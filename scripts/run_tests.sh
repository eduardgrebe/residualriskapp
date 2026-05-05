#!/usr/bin/env bash
# run_tests.sh — Run the full test suite (Go + Python) with verbose output.
#
# Usage:
#   bash scripts/run_tests.sh           # all tests (default)
#   bash scripts/run_tests.sh go        # Go tests only
#   bash scripts/run_tests.sh python    # Python tests only
#   bash scripts/run_tests.sh fast      # Go + Python tests that don't need
#                                         multiprocessing (safe in all envs)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# ── Colour helpers ────────────────────────────────────────────────────────────
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

GO_OK=0
PY_OK=0
GO_SKIPPED=0
PY_SKIPPED=0

# ── Go tests ──────────────────────────────────────────────────────────────────
run_go_tests() {
    rule
    info "Go tests  (go test -v ./...)"
    rule
    if ! command -v go &>/dev/null; then
        warn "go not found on PATH — skipping Go tests"
        GO_SKIPPED=1
        return
    fi
    pushd "$REPO_ROOT/go" >/dev/null
    if go test -v ./...; then
        GO_OK=1
    fi
    popd >/dev/null
}

# ── Python tests ──────────────────────────────────────────────────────────────
# Tests that use ProcessPoolExecutor fail in sandboxed environments with
# PermissionError on semaphore creation. They are identified by class/method
# names containing: Python, Agreement, Bootstrap, or agree_with_python.
SANDBOX_FILTER='not (Python or Agreement or Bootstrap or agree_with_python)'

run_python_tests() {
    local extra_args=("$@")
    rule
    info "Python tests  (uv run pytest -v${extra_args:+ ${extra_args[*]}})"
    rule
    if ! command -v uv &>/dev/null; then
        warn "uv not found — skipping Python tests"
        PY_SKIPPED=1
        return
    fi
    if uv run pytest tests/ -v "${extra_args[@]+"${extra_args[@]}"}"; then
        PY_OK=1
    fi
}

# ── Dispatch ──────────────────────────────────────────────────────────────────
case "$MODE" in
    go)
        run_go_tests
        ;;
    python)
        run_python_tests
        ;;
    fast)
        # Go tests are always safe; Python filtered to exclude ProcessPoolExecutor tests
        run_go_tests
        rule
        info "Skipping tests that require ProcessPoolExecutor (sandbox-incompatible)"
        info "Filter: ${SANDBOX_FILTER}"
        run_python_tests -k "$SANDBOX_FILTER"
        ;;
    all)
        run_go_tests
        run_python_tests
        ;;
    *)
        echo "Usage: $0 [all|go|python|fast]"
        exit 1
        ;;
esac

# ── Summary ───────────────────────────────────────────────────────────────────
rule
echo -e "${BOLD}Summary${NC}"
rule

if [[ "$MODE" == "all" || "$MODE" == "go" || "$MODE" == "fast" ]]; then
    if   [[ $GO_SKIPPED -eq 1 ]]; then warn "Go tests:     SKIPPED (go not found)"
    elif [[ $GO_OK      -eq 1 ]]; then pass "Go tests:     PASSED"
    else                               fail "Go tests:     FAILED"; fi
fi

if [[ "$MODE" == "all" || "$MODE" == "python" || "$MODE" == "fast" ]]; then
    if   [[ $PY_SKIPPED -eq 1 ]]; then warn "Python tests: SKIPPED (uv not found)"
    elif [[ $PY_OK      -eq 1 ]]; then pass "Python tests: PASSED"
    else                               fail "Python tests: FAILED (see output above)"; fi
fi

if [[ "$MODE" == "all" || "$MODE" == "python" ]]; then
    warn "Note: tests using ProcessPoolExecutor fail in sandboxed environments (PermissionError) — this is expected"
fi

rule

# Exit non-zero if any suite actually failed (not just skipped)
FAILED=0
if [[ ("$MODE" == "all" || "$MODE" == "go" || "$MODE" == "fast") \
      && $GO_SKIPPED -eq 0 && $GO_OK -eq 0 ]]; then
    FAILED=1
fi
if [[ ("$MODE" == "all" || "$MODE" == "python" || "$MODE" == "fast") \
      && $PY_SKIPPED -eq 0 && $PY_OK -eq 0 ]]; then
    FAILED=1
fi
exit $FAILED
