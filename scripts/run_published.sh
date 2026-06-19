#!/usr/bin/env bash
# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version. See <https://www.gnu.org/licenses/>.
#
# Pull a published Residual HIV-TT Risk Estimator image from ghcr.io and run it.
# (For running a locally *built* image instead, see docker/run.sh.)

set -euo pipefail

IMAGE="${IMAGE:-ghcr.io/eduardgrebe/residualrisk}"
PORT="${PORT:-8501}"

usage() {
    cat <<EOF
Pull and run a published Residual HIV-TT Risk Estimator container.

Usage: $(basename "$0") [VERSION]

  VERSION   Image tag to run (default: latest). A leading "v" is accepted and
            stripped, so both "v1.1.0a4" and "1.1.0a4" work — git tags are
            vX.Y.Z, the published image tags drop the "v".

Examples:
  $(basename "$0")               # :latest  (current production release)
  $(basename "$0") v1.1.0a4      # :1.1.0a4 (a specific pre-release)
  PORT=8600 $(basename "$0")     # publish on a different host port

Environment overrides:
  IMAGE   container repository (default: ghcr.io/eduardgrebe/residualrisk)
  PORT    host port to publish on (default: 8501)

The app is served at http://localhost:\$PORT (bound to localhost only).
Press Ctrl-C to stop; the container is removed automatically on exit.
EOF
}

case "${1:-}" in
    -h|--help) usage; exit 0 ;;
esac

VERSION="${1:-latest}"
# git tags are vX.Y.Z; published image tags drop the leading "v".
if [[ "$VERSION" =~ ^v[0-9] ]]; then
    TAG="${VERSION#v}"
else
    TAG="$VERSION"
fi
REF="${IMAGE}:${TAG}"

# Pre-flight: Docker CLI present and daemon reachable.
if ! command -v docker >/dev/null 2>&1; then
    echo "error: 'docker' not found on PATH. Install Docker, Colima, or OrbStack." >&2
    exit 1
fi
if ! docker info >/dev/null 2>&1; then
    echo "error: cannot reach the Docker daemon — is it running?" >&2
    echo "       (e.g. 'colima start', or launch OrbStack / Docker Desktop)" >&2
    exit 1
fi

echo "Pulling ${REF} ..."
if ! docker pull "$REF"; then
    echo "error: failed to pull ${REF}." >&2
    echo "       Check the tag exists and the package is public" >&2
    echo "       (private packages need: docker login ghcr.io)." >&2
    exit 1
fi

echo
echo "Starting ${REF}"
echo "  -> open http://localhost:${PORT} in your browser"
echo "  -> press Ctrl-C to stop"
echo
exec docker run --rm -p "127.0.0.1:${PORT}:8501" "$REF"
