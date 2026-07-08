#!/usr/bin/env bash
# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# deploy.sh — pull the pinned image, redeploy, and reclaim disk.
#
# Keep this next to the server's docker-compose.yml. The image tag is pinned in
# that file (`latest` does not track the 1.1.0aX pre-releases), so bump the tag
# there first, then run:
#
#     ./deploy.sh            # or: sudo ./deploy.sh
#
# Order is deliberate and safe: pull -> up -d -> prune. The new container is
# running (and therefore references the new image) before the prune, so
# `docker image prune -af` only removes images no container uses. `set -euo
# pipefail` aborts before deploying on a failed pull, and before pruning on a
# failed up.

set -euo pipefail

GREEN='\033[0;32m'
NC='\033[0m'
step() { echo -e "${GREEN}==> $*${NC}"; }

# Operate from this script's own directory (where docker-compose.yml lives),
# so it works no matter where you invoke it from.
cd "$(dirname "$(readlink -f "$0")")"

if [[ ! -f docker-compose.yml ]]; then
  echo "error: docker-compose.yml not found in $(pwd)" >&2
  echo "Place this script alongside the server's docker-compose.yml." >&2
  exit 1
fi

step "Pulling image(s) pinned in docker-compose.yml ..."
sudo docker compose pull

step "Redeploying (detached) ..."
sudo docker compose up -d

step "Reclaiming disk ..."
# -a also removes unused *tagged* images — the older 1.1.0aX builds that no
# running container references — which is what frees space as alphas pile up.
# The live container keeps its image; to roll back, re-pin the old tag and pull.
# Volumes and data are never touched. For a deeper clean (stopped containers,
# unused networks, build cache) swap this for: sudo docker system prune -af
sudo docker image prune -af

step "Done. Status and disk usage:"
sudo docker compose ps
sudo docker system df
