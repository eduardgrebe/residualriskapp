#!/bin/bash
# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# Docker build script with multi-architecture support

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
IMAGE_NAME="${IMAGE_NAME:-residualrisk}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
REGISTRY="${REGISTRY:-}"  # e.g., docker.io/username or ghcr.io/username

# Parse command line arguments
PLATFORMS="linux/amd64,linux/arm64"
PUSH=false
LOAD=false
AUTO_PLATFORM=false
BUILD_TYPE="local"

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -p, --platform PLATFORMS   Comma-separated platforms (default: linux/amd64,linux/arm64)"
    echo "  -t, --tag TAG              Image tag (default: latest)"
    echo "  -r, --registry REGISTRY    Registry to push to (e.g., docker.io/username)"
    echo "  --push                     Push to registry after build"
    echo "  --load                     Load image to local Docker (auto-detects platform)"
    echo "  -h, --help                 Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Build for local use (simplest - auto-detects your architecture)"
    echo "  $0 --load"
    echo ""
    echo "  # Build multi-arch and push to Docker Hub"
    echo "  $0 --registry docker.io/username --push"
    echo ""
    echo "  # Build for specific architecture"
    echo "  $0 --platform linux/amd64 --load"
    echo "  $0 --platform linux/arm64 --load"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--platform)
            PLATFORMS="$2"
            shift 2
            ;;
        -t|--tag)
            IMAGE_TAG="$2"
            shift 2
            ;;
        -r|--registry)
            REGISTRY="$2"
            shift 2
            ;;
        --push)
            PUSH=true
            shift
            ;;
        --load)
            LOAD=true
            AUTO_PLATFORM=true
            shift
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo -e "${RED}Error: Unknown option $1${NC}"
            print_usage
            exit 1
            ;;
    esac
done

# Validate conflicting options
if [ "$PUSH" = true ] && [ "$LOAD" = true ]; then
    echo -e "${RED}Error: Cannot use --push and --load together${NC}"
    exit 1
fi

# Auto-detect platform for --load if not specified
if [ "$LOAD" = true ] && [ "$AUTO_PLATFORM" = true ]; then
    # Detect current architecture
    ARCH=$(uname -m)
    case "$ARCH" in
        x86_64)
            PLATFORMS="linux/amd64"
            echo -e "${GREEN}Auto-detected platform: linux/amd64${NC}"
            ;;
        aarch64|arm64)
            PLATFORMS="linux/arm64"
            echo -e "${GREEN}Auto-detected platform: linux/arm64${NC}"
            ;;
        *)
            echo -e "${RED}Error: Unsupported architecture: $ARCH${NC}"
            echo "Please specify platform manually with --platform"
            exit 1
            ;;
    esac
fi

if [ "$LOAD" = true ] && [[ "$PLATFORMS" == *","* ]]; then
    echo -e "${RED}Error: --load requires a single platform${NC}"
    echo "You specified: $PLATFORMS"
    echo "Either:"
    echo "  1. Use --load without --platform (auto-detects your architecture)"
    echo "  2. Specify a single platform: --platform linux/amd64 --load"
    exit 1
fi

# Construct full image name
if [ -n "$REGISTRY" ]; then
    FULL_IMAGE_NAME="${REGISTRY}/${IMAGE_NAME}:${IMAGE_TAG}"
else
    FULL_IMAGE_NAME="${IMAGE_NAME}:${IMAGE_TAG}"
fi

echo -e "${GREEN}=== Docker Build Configuration ===${NC}"
echo "Image: ${FULL_IMAGE_NAME}"
echo "Platforms: ${PLATFORMS}"
echo "Push: ${PUSH}"
echo "Load: ${LOAD}"
echo ""

# Check if buildx is available
if ! docker buildx version > /dev/null 2>&1; then
    echo -e "${RED}Error: docker buildx is not available${NC}"
    echo "Install it with: docker buildx install"
    exit 1
fi

# Create/use buildx builder
BUILDER_NAME="residualrisk-builder"
if ! docker buildx inspect "$BUILDER_NAME" > /dev/null 2>&1; then
    echo -e "${YELLOW}Creating new buildx builder: ${BUILDER_NAME}${NC}"
    docker buildx create --name "$BUILDER_NAME" --driver docker-container --use
else
    echo -e "${GREEN}Using existing buildx builder: ${BUILDER_NAME}${NC}"
    docker buildx use "$BUILDER_NAME"
fi

# Bootstrap builder if needed
docker buildx inspect --bootstrap

# Build command
BUILD_CMD="docker buildx build"
BUILD_CMD="$BUILD_CMD --platform ${PLATFORMS}"
BUILD_CMD="$BUILD_CMD --tag ${FULL_IMAGE_NAME}"

if [ "$PUSH" = true ]; then
    BUILD_CMD="$BUILD_CMD --push"
    echo -e "${YELLOW}Image will be pushed to registry${NC}"
fi

if [ "$LOAD" = true ]; then
    BUILD_CMD="$BUILD_CMD --load"
    echo -e "${YELLOW}Image will be loaded to local Docker${NC}"
fi

BUILD_CMD="$BUILD_CMD ."

echo -e "${GREEN}=== Starting Docker Build ===${NC}"
echo "Command: $BUILD_CMD"
echo ""

# Execute build
eval $BUILD_CMD

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}=== Build Successful ===${NC}"
    echo "Image: ${FULL_IMAGE_NAME}"

    if [ "$LOAD" = true ]; then
        echo ""
        echo "Run with: docker run -p 127.0.0.1:8501:8501 ${FULL_IMAGE_NAME}"
        echo "Or use docker-compose: docker-compose up"
    fi

    if [ "$PUSH" = true ]; then
        echo ""
        echo "Image pushed to: ${FULL_IMAGE_NAME}"
        echo "Pull with: docker pull ${FULL_IMAGE_NAME}"
    fi
else
    echo -e "${RED}=== Build Failed ===${NC}"
    exit 1
fi
