# Docker Quick Reference

Quick commands for common Docker operations with the Residual Risk Estimator.

> **New to Docker?** See [BUILD-QUICK-START.md](BUILD-QUICK-START.md) for a beginner-friendly guide.

## 🚀 Quick Start

```bash
# Fastest way to get started
docker-compose up -d
```

Access at: http://localhost:8501

## 📦 Building

```bash
# Build for your machine (simplest - auto-detects your architecture)
./docker/build.sh --load

# Build for multiple architectures (amd64 + arm64) - cannot load locally
./docker/build.sh

# Build and push to Docker Hub (multi-arch)
./docker/build.sh --registry docker.io/username --push

# Build for specific platform (if you need to cross-compile)
./docker/build.sh --platform linux/amd64 --load   # Intel/AMD
./docker/build.sh --platform linux/arm64 --load   # ARM/Apple Silicon
```

## ▶️ Running

```bash
# Using docker-compose (recommended)
docker-compose up -d

# Using helper script
./docker/run.sh

# Using docker directly
docker run -d -p 127.0.0.1:8501:8501 --name residualrisk_app residualrisk:latest

# Custom port
docker run -d -p 127.0.0.1:8080:8501 residualrisk:latest
```

## 🔍 Monitoring

```bash
# View logs
docker logs -f residualrisk_app
docker-compose logs -f

# Check health
docker inspect --format='{{.State.Health.Status}}' residualrisk_app

# Resource usage
docker stats residualrisk_app
```

## 🛑 Stopping

```bash
# Stop container
docker stop residualrisk_app
docker-compose down

# Stop and remove
docker rm -f residualrisk_app
docker-compose down -v

# Restart
docker restart residualrisk_app
docker-compose restart
```

## 🔧 Troubleshooting

```bash
# View full logs
docker logs residualrisk_app

# Check if Go binary is present
docker exec residualrisk_app ls -lh go/bin/riskdays_go

# Test health endpoint
curl http://localhost:8501/_stcore/health

# Shell into container
docker exec -it residualrisk_app /bin/bash

# Rebuild from scratch
docker-compose build --no-cache
```

## 🌐 Multi-Architecture

```bash
# Check supported architectures
docker buildx imagetools inspect residualrisk:latest

# Build for ARM64 on AMD64 machine (or vice versa)
./docker/build.sh --platform linux/arm64 --registry docker.io/user --push

# Pull and run multi-arch image (automatic platform detection)
docker pull docker.io/user/residualrisk:latest
docker run -d -p 127.0.0.1:8501:8501 docker.io/user/residualrisk:latest
```

## 📝 Configuration

```bash
# With environment variables
docker run -d -p 127.0.0.1:8501:8501 \
  -e STREAMLIT_SERVER_MAX_UPLOAD_SIZE=500 \
  residualrisk:latest

# With custom config file
docker run -d -p 127.0.0.1:8501:8501 \
  -v $(pwd)/config.toml:/app/.streamlit/config.toml:ro \
  residualrisk:latest

# With resource limits (NOT recommended - reduces performance)
docker run -d -p 127.0.0.1:8501:8501 \
  --cpus="8" --memory="8g" \
  residualrisk:latest

# No limits (recommended - uses all CPUs)
docker run -d -p 127.0.0.1:8501:8501 residualrisk:latest
```

## 🏷️ Tagging & Publishing

```bash
# Tag image
docker tag residualrisk:latest username/residualrisk:v1.0.0

# Push to Docker Hub
docker push username/residualrisk:v1.0.0

# Push to GitHub Container Registry
docker tag residualrisk:latest ghcr.io/username/residualrisk:latest
docker push ghcr.io/username/residualrisk:latest
```

## 🧹 Cleanup

```bash
# Remove stopped container
docker rm residualrisk_app

# Remove image
docker rmi residualrisk:latest

# Remove all (containers + images + volumes)
docker-compose down --rmi all -v

# Clean build cache
docker buildx prune -f
docker system prune -a -f
```

## 📊 Using Makefile

If you prefer using Make:

```bash
# Show all available commands
make -f docker/Makefile help

# Common operations
make -f docker/Makefile build       # Build for local
make -f docker/Makefile run         # Start with docker-compose
make -f docker/Makefile logs        # View logs
make -f docker/Makefile stop        # Stop container
make -f docker/Makefile test        # Build and test
make -f docker/Makefile clean       # Remove everything

# Push to registry
make -f docker/Makefile push REGISTRY=docker.io/username
```

## 🔗 Useful Links

- Full documentation: [docker/README.md](README.md)
- nginx reverse proxy config: [docker/nginx/conf.d/app.conf](nginx/conf.d/app.conf)
- Caddy reverse proxy config: [docker/Caddyfile](Caddyfile)
- Application README: [../README.md](../README.md)

## 💡 Tips

1. **First time setup**: Use `docker-compose up -d` - it's the simplest
2. **Production**: Build multi-arch and push to a registry
3. **Development**: Mount local code with volumes for hot reload
4. **Performance**:
   - Container has access to ALL host CPUs by default (see CPU-PERFORMANCE.md)
   - Ensure the Go implementation is being used (check UI)
   - Use 250K+ simulations to see full multiprocessing benefits
5. **Security**: Container runs as non-root user (appuser, UID 1000)

## ⚡ One-Liners

```bash
# Build, run, and open browser
./docker/build.sh --load && ./docker/run.sh && open http://localhost:8501

# Quick restart
docker restart residualrisk_app

# Quick rebuild
docker-compose up -d --build

# Check if running
docker ps | grep residualrisk

# Stream logs until error
docker logs -f residualrisk_app 2>&1 | grep -i error

# Get container IP
docker inspect -f '{{range.NetworkSettings.Networks}}{{.IPAddress}}{{end}}' residualrisk_app
```