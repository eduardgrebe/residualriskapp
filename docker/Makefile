# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# Makefile for Docker operations
# Run from repository root: make -f docker/Makefile [target]

.PHONY: help build build-multiarch build-dev run stop logs clean push test

# Configuration
IMAGE_NAME ?= residualrisk
IMAGE_TAG ?= latest
REGISTRY ?=
PLATFORMS ?= linux/amd64,linux/arm64

help: ## Show this help message
	@echo "Residual Risk Estimator - Docker Commands"
	@echo ""
	@echo "Usage: make -f docker/Makefile [target]"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  %-20s %s\n", $$1, $$2}'

build: ## Build Docker image for local architecture
	./docker/build.sh --load

build-multiarch: ## Build for multiple architectures (amd64, arm64)
	./docker/build.sh --platform $(PLATFORMS)

build-dev: ## Build for local development with current architecture
	docker buildx build --platform $$(uname -m | sed 's/x86_64/linux\/amd64/;s/aarch64/linux\/arm64/') --tag $(IMAGE_NAME):dev --load .

run: ## Run container using docker-compose
	docker-compose up -d

stop: ## Stop running container
	docker-compose down

restart: ## Restart container
	docker-compose restart

logs: ## Show container logs
	docker-compose logs -f

logs-tail: ## Show last 100 lines of logs
	docker-compose logs --tail=100

shell: ## Open shell in running container
	docker-compose exec residualrisk /bin/bash

stats: ## Show container resource usage
	docker stats $$(docker-compose ps -q) --no-stream

ps: ## Show running containers
	docker-compose ps

health: ## Check container health status
	@docker inspect --format='{{.State.Health.Status}}' residualrisk_app 2>/dev/null || echo "Container not running"

clean: ## Remove containers and images
	docker-compose down --rmi local -v

clean-all: ## Remove everything including volumes
	docker-compose down --rmi all -v
	docker buildx prune -f

push: ## Build multi-arch and push to registry
ifndef REGISTRY
	$(error REGISTRY is not set. Use: make push REGISTRY=docker.io/username)
endif
	./docker/build.sh --platform $(PLATFORMS) --registry $(REGISTRY) --push

push-latest: ## Push latest tag to registry
ifndef REGISTRY
	$(error REGISTRY is not set. Use: make push-latest REGISTRY=docker.io/username)
endif
	./docker/build.sh --platform $(PLATFORMS) --registry $(REGISTRY) --tag latest --push

push-ghcr: ## Build multi-arch and push to ghcr.io/vitalant-research-institute
	./docker/build.sh --platform $(PLATFORMS) --registry ghcr.io/vitalant-research-institute --push

push-ghcr-tag: ## Build multi-arch and push a specific version tag to ghcr.io (usage: make push-ghcr-tag IMAGE_TAG=1.2.3)
ifndef IMAGE_TAG
	$(error IMAGE_TAG is not set. Use: make push-ghcr-tag IMAGE_TAG=1.2.3)
endif
	./docker/build.sh --platform $(PLATFORMS) --registry ghcr.io/vitalant-research-institute --tag $(IMAGE_TAG) --push

test: ## Test the container
	@echo "Building and testing container..."
	docker-compose up -d
	@echo "Waiting for container to be healthy..."
	@timeout 60 sh -c 'until docker inspect --format="{{.State.Health.Status}}" residualrisk_app | grep -q healthy; do sleep 2; done' || (echo "Container failed to become healthy" && docker-compose logs && exit 1)
	@echo "Testing health endpoint..."
	curl -f http://localhost:8501/_stcore/health || (echo "Health check failed" && docker-compose logs && exit 1)
	@echo "✓ Container is healthy and responding"
	docker-compose down

rebuild: ## Rebuild without cache
	docker-compose build --no-cache
	docker-compose up -d

update: ## Pull latest changes and rebuild
	git pull
	$(MAKE) rebuild

version: ## Show version information
	@echo "Image: $(IMAGE_NAME):$(IMAGE_TAG)"
	@echo "Registry: $(REGISTRY)"
	@echo "Platforms: $(PLATFORMS)"
	@docker --version
	@docker-compose --version
	@docker buildx version

# Development helpers
dev-up: build-dev ## Build and run in development mode
	docker run -d \
		--name residualrisk_dev \
		-p 127.0.0.1:8501:8501 \
		-v $$(pwd):/app \
		$(IMAGE_NAME):dev

dev-down: ## Stop development container
	docker rm -f residualrisk_dev

# Quick shortcuts
up: run ## Alias for run
down: stop ## Alias for stop
