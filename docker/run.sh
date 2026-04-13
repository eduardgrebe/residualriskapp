#!/bin/bash
# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# Quick run script for Docker container

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Configuration
IMAGE_NAME="${IMAGE_NAME:-residualrisk:latest}"
CONTAINER_NAME="${CONTAINER_NAME:-residualrisk_app}"
PORT="${PORT:-8501}"

# Check if container is already running
if docker ps --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo -e "${YELLOW}Container ${CONTAINER_NAME} is already running${NC}"
    echo "Stop it with: docker stop ${CONTAINER_NAME}"
    exit 1
fi

# Remove existing stopped container with same name
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo -e "${YELLOW}Removing stopped container: ${CONTAINER_NAME}${NC}"
    docker rm "${CONTAINER_NAME}"
fi

echo -e "${GREEN}Starting Residual Risk Estimator...${NC}"
echo "Image: ${IMAGE_NAME}"
echo "Port: ${PORT}"
echo ""

# Run container
docker run -d \
    --name "${CONTAINER_NAME}" \
    -p "127.0.0.1:${PORT}:8501" \
    --restart unless-stopped \
    --health-cmd="curl -f http://localhost:8501/_stcore/health || exit 1" \
    --health-interval=30s \
    --health-timeout=10s \
    --health-retries=3 \
    --health-start-period=40s \
    "${IMAGE_NAME}"

echo ""
echo -e "${GREEN}Container started successfully!${NC}"
echo ""
echo "Access the application at: http://localhost:${PORT}"
echo ""
echo "Useful commands:"
echo "  View logs:    docker logs -f ${CONTAINER_NAME}"
echo "  Stop:         docker stop ${CONTAINER_NAME}"
echo "  Restart:      docker restart ${CONTAINER_NAME}"
echo "  Remove:       docker rm -f ${CONTAINER_NAME}"
echo ""

# Show logs
echo "Showing startup logs (Ctrl+C to exit log view)..."
sleep 2
docker logs -f "${CONTAINER_NAME}"
