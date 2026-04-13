# Docker Build Quick Start

**TL;DR for building locally:**

```bash
./docker/build.sh --load
```

That's it! The script auto-detects your architecture and builds the right image.

---

## What This Does

When you run `./docker/build.sh --load`:

1. **Detects your architecture** (amd64 or arm64)
2. **Builds Docker image** for that platform only
3. **Compiles Go binary** for your architecture
4. **Loads image** into your local Docker

On your system (arm64/Apple Silicon):
- Builds for `linux/arm64`
- Go binary compiled for ARM64
- ~3-5 minutes to build

---

## Common Commands

### Local Development
```bash
# Build for your machine
./docker/build.sh --load

# Run the container
./docker/run.sh

# Or use docker-compose (builds automatically)
docker-compose up -d
```

### Multi-Architecture (for registries)
```bash
# Build for both amd64 and arm64, push to Docker Hub
./docker/build.sh --registry docker.io/username --push

# Build for both amd64 and arm64, push to GitHub Container Registry
./docker/build.sh --registry ghcr.io/username --push
```

### Specific Platform (cross-compilation)
```bash
# Build for Intel/AMD (even if you're on Apple Silicon)
./docker/build.sh --platform linux/amd64 --load

# Build for ARM (even if you're on Intel/AMD)
./docker/build.sh --platform linux/arm64 --load
```

---

## Understanding the Error Message

If you see:
```
Error: --load requires a single platform
```

**Cause:** You tried to build for multiple platforms with `--load`
- `--load` saves image to local Docker
- Docker can only load ONE architecture at a time
- Multi-arch images need to be pushed to a registry

**Solution:** Either:
1. Use `--load` without specifying `--platform` (auto-detects)
2. Use `--load` with a single platform: `--platform linux/arm64`
3. Push multi-arch to registry: `--push` (don't use `--load`)

---

## Why Auto-Detection?

Docker buildx has a limitation: you can only load single-platform images locally.

**Before (confusing):**
```bash
./docker/build.sh --load
# Error: --load requires a single platform
# Specify a single platform with --platform (e.g., --platform linux/amd64)
```

**Now (simple):**
```bash
./docker/build.sh --load
# Auto-detected platform: linux/arm64
# ✓ Building...
```

---

## Full Help

```bash
./docker/build.sh --help
```

Shows all options and examples.
