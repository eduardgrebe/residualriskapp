#!/usr/bin/env bash
# Build the residualrisk Go binary at go/bin/riskdays_go. Idempotent.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
make -C "$REPO_ROOT/go" build
echo "Built: $REPO_ROOT/go/bin/riskdays_go"
