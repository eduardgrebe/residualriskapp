#!/usr/bin/env bash
# audit_deps.sh — audit Python + Go project dependencies for known
# vulnerabilities and supply-chain risk (compromised / malicious packages).
#
# Defence layers:
#   1. Lockfile integrity   — `uv lock --check` + `go mod verify`
#                             (detects lockfile drift + tampered module cache)
#   2. Python vulnerabilities — `osv-scanner` on uv.lock (OSV / NVD / GHSA
#                             + MAL-* malicious-release advisories)
#                             + `pip-audit` against the PyPI advisory DB
#                             (independent source; cross-check via `uvx`,
#                              no persistent install)
#   3. Go vulnerabilities   — `osv-scanner` on go/go.mod
#                             + `govulncheck` (reachable-only, call-graph)
#                             + `snyk test` on go/go.mod (independent DB)
#
# We deliberately use *two* scanners per language, backed by *different*
# vulnerability databases, so the two cross-check each other. Every scanner
# runs independently: missing tools SKIP, findings FAIL, but the script keeps
# going so one run produces a complete report. Exit code is non-zero if any
# scanner failed (not just skipped).
#
# Usage:
#   bash scripts/audit_deps.sh          # run all layers
#   bash scripts/audit_deps.sh python   # Python only (lockfile + vulns)
#   bash scripts/audit_deps.sh go       # Go only     (lockfile + vulns)
#   bash scripts/audit_deps.sh lockfile # lockfile integrity only

set -uo pipefail  # not -e: we want every scanner to run even if others fail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# ── Colour helpers (match scripts/run_tests.sh) ───────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

pass() { echo -e "${GREEN}${BOLD}✔  $*${NC}"; }
fail() { echo -e "${RED}${BOLD}✘  $*${NC}"; }
info() { echo -e "${CYAN}▶  $*${NC}"; }
warn() { echo -e "${YELLOW}⚠  $*${NC}"; }
rule() { echo -e "${BOLD}────────────────────────────────────────────────────${NC}"; }

MODE="${1:-all}"
case "$MODE" in
    all|python|go|lockfile) ;;
    *) echo "Usage: $0 [all|python|go|lockfile]"; exit 1 ;;
esac

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S %Z')

# Scratch dir for generated artifacts (e.g. requirements.txt for pip-audit,
# pip-audit's own on-disk cache). Cleared on exit — never leaves the repo dirty.
TMPDIR_AUDIT=$(mktemp -d)
trap 'rm -rf "$TMPDIR_AUDIT"' EXIT

echo
rule
echo -e "${BOLD} Dependency audit — $TIMESTAMP${NC}"
echo -e "${BOLD} Repo:  $REPO_ROOT${NC}"
echo -e "${BOLD} Mode:  $MODE${NC}"
rule

# Per-section result codes: 0=PASS, 1=FAIL, 2=SKIPPED (never assign 2 to a
# variable that isn't guarded in the summary — see block at the bottom).
LOCK_UV=2
LOCK_GO=2
PY_OSV=2
PY_AUDIT=2
GO_OSV=2
GO_VULN=2
GO_SNYK=2

have() { command -v "$1" &>/dev/null; }

# ── Layer 1: lockfile integrity ───────────────────────────────────────────────
run_lockfile_checks() {
    rule
    info "Lockfile integrity"
    rule

    # uv lock --check verifies uv.lock is consistent with pyproject.toml. It
    # does NOT re-verify already-downloaded artifact hashes — for that, use
    # `uv sync --locked --reinstall` (heavy; save for CI / release gating).
    echo
    info "→ uv lock --check   (Python: uv.lock ⇔ pyproject.toml)"
    if ! have uv; then
        warn "uv not on PATH — skipping"
    elif uv lock --check; then
        pass "uv.lock is in sync with pyproject.toml"
        LOCK_UV=0
    else
        fail "uv.lock is out of sync with pyproject.toml (run 'uv lock' and review the diff)"
        LOCK_UV=1
    fi

    # go mod verify checks that downloaded modules in the local module cache
    # match the hashes recorded in go.sum. Detects tampered / corrupt cache.
    echo
    info "→ go mod verify     (Go: module cache ⇔ go.sum)"
    if ! have go; then
        warn "go not on PATH — skipping"
    else
        pushd "$REPO_ROOT/go" >/dev/null
        if go mod verify; then
            pass "Go module cache matches go.sum"
            LOCK_GO=0
        else
            fail "Go module cache does NOT match go.sum (tampered / corrupt cache)"
            LOCK_GO=1
        fi
        popd >/dev/null
    fi
}

# ── Layer 2: Python vulnerability + supply-chain scans ────────────────────────
run_python_scans() {
    rule
    info "Python — vulnerability + supply-chain scans"
    rule

    # osv-scanner reads uv.lock natively. OSV covers PyPI advisories (PYSEC,
    # GHSA) AND MAL-* advisories for known-malicious / compromised releases.
    echo
    info "→ osv-scanner       (OSV / NVD / GHSA + MAL-* malicious releases)"
    if ! have osv-scanner; then
        warn "osv-scanner not on PATH — skipping   (install: brew install osv-scanner)"
    elif osv-scanner scan source --lockfile "$REPO_ROOT/uv.lock"; then
        pass "osv-scanner: no known vulnerabilities in Python dependencies"
        PY_OSV=0
    else
        fail "osv-scanner reported vulnerabilities in Python dependencies (see table above)"
        PY_OSV=1
    fi

    # pip-audit (PyPA) using PyPI's advisory DB as an independent cross-check
    # against osv-scanner (which uses OSV). Run via `uvx` so no persistent
    # install is needed. `--no-deps` because uv.lock is already fully flattened
    # (no need for pip to re-resolve). `--cache-dir` inside our tmpdir keeps
    # us working in restricted environments where ~/Library/Caches is blocked.
    echo
    info "→ pip-audit         (PyPA — PyPI advisory DB, independent of OSV)"
    if ! have uv; then
        warn "uv not on PATH — can't run 'uvx pip-audit' — skipping"
    else
        REQ_TXT="$TMPDIR_AUDIT/requirements.txt"
        if ! uv export --format requirements-txt --no-hashes --no-emit-project \
                --quiet -o "$REQ_TXT" 2>&1; then
            fail "'uv export' failed — cannot run pip-audit on Python deps"
            PY_AUDIT=1
        elif uvx --quiet pip-audit \
                -r "$REQ_TXT" \
                --no-deps \
                --disable-pip \
                --cache-dir="$TMPDIR_AUDIT/pip-audit-cache" \
                --vulnerability-service=pypi; then
            pass "pip-audit: no vulnerabilities in Python dependencies (PyPI DB)"
            PY_AUDIT=0
        else
            fail "pip-audit reported vulnerabilities in Python dependencies (see output above)"
            PY_AUDIT=1
        fi
    fi
}

# ── Layer 3: Go vulnerability + supply-chain scans ────────────────────────────
run_go_scans() {
    rule
    info "Go — vulnerability + supply-chain scans"
    rule

    # osv-scanner reads go.mod / go.sum. Catches all known CVEs regardless of
    # reachability from our code.
    echo
    info "→ osv-scanner       (all CVEs affecting Go deps, reachable or not)"
    if ! have osv-scanner; then
        warn "osv-scanner not on PATH — skipping   (install: brew install osv-scanner)"
    elif osv-scanner scan source --lockfile "$REPO_ROOT/go/go.mod"; then
        pass "osv-scanner: no known vulnerabilities in Go dependencies"
        GO_OSV=0
    else
        fail "osv-scanner reported vulnerabilities in Go dependencies (see table above)"
        GO_OSV=1
    fi

    # govulncheck uses static call-graph analysis to filter to CVEs the binary
    # actually reaches — dramatically less noise than osv-scanner's raw list.
    # Uses `go run` if not installed → transient download, no GOBIN write.
    echo
    info "→ govulncheck       (Go-source call-graph analysis; reachable-only)"
    if have govulncheck; then
        GOVULNCHECK_CMD=(govulncheck ./...)
    elif have go; then
        info "  govulncheck not installed — using 'go run golang.org/x/vuln/cmd/govulncheck@latest'"
        GOVULNCHECK_CMD=(go run golang.org/x/vuln/cmd/govulncheck@latest ./...)
    else
        GOVULNCHECK_CMD=()
    fi
    if [[ ${#GOVULNCHECK_CMD[@]} -eq 0 ]]; then
        warn "neither govulncheck nor go on PATH — skipping"
    else
        pushd "$REPO_ROOT/go" >/dev/null
        if "${GOVULNCHECK_CMD[@]}"; then
            pass "govulncheck: no reachable vulnerabilities in Go source"
            GO_VULN=0
        else
            fail "govulncheck reported reachable vulnerabilities in Go source (see output above)"
            GO_VULN=1
        fi
        popd >/dev/null
    fi

    # Snyk also scans Go modules; independent DB → cross-check.
    echo
    info "→ snyk test         (Snyk vulnerability + supply-chain DB)"
    if ! have snyk; then
        warn "snyk not on PATH — skipping   (install: brew install snyk-cli && snyk auth)"
    else
        pushd "$REPO_ROOT/go" >/dev/null
        if snyk test --file=go.mod; then
            pass "snyk: no vulnerabilities in Go dependencies"
            GO_SNYK=0
        else
            rc=$?
            if [[ $rc -eq 1 ]]; then
                fail "snyk reported vulnerabilities in Go dependencies (see output above)"
            else
                fail "snyk failed (exit $rc) — check auth ('snyk auth') and network"
            fi
            GO_SNYK=1
        fi
        popd >/dev/null
    fi
}

# ── Dispatch ──────────────────────────────────────────────────────────────────
case "$MODE" in
    lockfile)
        run_lockfile_checks
        ;;
    python)
        run_lockfile_checks
        run_python_scans
        ;;
    go)
        run_lockfile_checks
        run_go_scans
        ;;
    all)
        run_lockfile_checks
        run_python_scans
        run_go_scans
        ;;
esac

# ── Summary ───────────────────────────────────────────────────────────────────
echo
rule
echo -e "${BOLD}Summary${NC}"
rule

report() {
    local label="$1" code="$2"
    case "$code" in
        0) pass "$label: PASSED" ;;
        1) fail "$label: FAILED" ;;
        2) warn "$label: SKIPPED" ;;
    esac
}

if [[ "$MODE" == "all" || "$MODE" == "python" || "$MODE" == "go" || "$MODE" == "lockfile" ]]; then
    report "uv lock --check         " $LOCK_UV
    report "go mod verify           " $LOCK_GO
fi
if [[ "$MODE" == "all" || "$MODE" == "python" ]]; then
    report "Python osv-scanner      " $PY_OSV
    report "Python pip-audit        " $PY_AUDIT
fi
if [[ "$MODE" == "all" || "$MODE" == "go" ]]; then
    report "Go     osv-scanner      " $GO_OSV
    report "Go     govulncheck      " $GO_VULN
    report "Go     snyk test        " $GO_SNYK
fi

rule
echo -e "${BOLD}Reminders${NC}"
echo "  • SKIPPED sections mean the scanner isn't installed — install and re-run"
echo "    to close the coverage gap (see the 'install:' hints above)."
echo "  • This audit checks lockfile integrity + published vulnerability DBs."
echo "    It does NOT re-verify artifact hashes for already-downloaded packages."
echo "    For release gating, additionally run:  uv sync --locked --reinstall"
rule

# Exit non-zero if any executed scanner actually FAILED (code 1).
# SKIPPED (code 2) does not fail the run — reported clearly in the summary.
FAILED=0
for code in $LOCK_UV $LOCK_GO $PY_OSV $PY_AUDIT $GO_OSV $GO_VULN $GO_SNYK; do
    [[ $code -eq 1 ]] && FAILED=1
done
exit $FAILED
