# Docker Deployment Guide

This guide covers building and running the Residual Risk Estimator webapp in Docker containers with multi-architecture support.

> **Just want to build and run?** See [BUILD-QUICK-START.md](BUILD-QUICK-START.md) for the simplest instructions.

## Quick Start

### Using Docker Compose (Recommended for Local Development)

```bash
# Build and start the container
docker-compose up -d

# View logs
docker-compose logs -f

# Stop the container
docker-compose down
```

Access the application at: http://localhost:8501

**Note:** The container has access to **all host CPU cores by default** for maximum multiprocessing performance. See [CPU-PERFORMANCE.md](CPU-PERFORMANCE.md) for details.

### Using Helper Scripts

```bash
# Build for local architecture (auto-detects your system)
./docker/build.sh --load

# Run the container
./docker/run.sh
```

**Note:** `--load` automatically detects your architecture (amd64 or arm64) and builds for that platform only.

## Multi-Architecture Builds

The Docker image supports both AMD64 and ARM64 architectures, making it deployable on:
- x86_64 servers (Intel/AMD)
- ARM64 servers (AWS Graviton, Apple Silicon, etc.)

### Prerequisites

1. **Docker with buildx** (included in Docker Desktop and recent Docker Engine versions)
   ```bash
   docker buildx version
   ```

2. **For multi-arch builds**: QEMU emulation (usually pre-installed)
   ```bash
   docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
   ```

### Build Options

#### 1. Local Development (Single Architecture)

**Simplest method - auto-detects your architecture:**

```bash
./docker/build.sh --load
```

This will:
- Detect if you're on amd64 (Intel/AMD) or arm64 (Apple Silicon/ARM)
- Build for that platform only
- Load the image to your local Docker

**Or manually specify platform:**
```bash
./docker/build.sh --platform linux/amd64 --load   # For Intel/AMD
./docker/build.sh --platform linux/arm64 --load   # For ARM/Apple Silicon
```

**Or use docker buildx directly:**
```bash
docker buildx build --platform linux/amd64 --tag residualrisk:latest --load .
```

#### 2. Multi-Architecture Build (No Push)

Build for both AMD64 and ARM64:

```bash
./docker/build.sh
```

**Important:** Multi-arch builds create a manifest but **cannot be loaded to local Docker** (Docker buildx limitation). Use this for testing the build process before pushing to a registry.

To test locally, use `--load` (builds for your architecture only).

#### 3. Build and Push to Registry

Push to Docker Hub:
```bash
./docker/build.sh --registry docker.io/yourusername --push
```

Push to GitHub Container Registry:
```bash
./docker/build.sh --registry ghcr.io/yourusername --push
```

#### 4. Build Specific Architecture

```bash
# AMD64 only
./docker/build.sh --platform linux/amd64 --load

# ARM64 only
./docker/build.sh --platform linux/arm64 --load
```

### Advanced Build Options

```bash
# Custom tag
./docker/build.sh --tag v1.0.0 --load

# Multiple platforms with custom registry
./docker/build.sh \
  --platform linux/amd64,linux/arm64 \
  --registry ghcr.io/vitalant \
  --tag latest \
  --push
```

## Running the Container

### Option 1: Docker Compose

Edit `docker-compose.yml` to customize settings, then:

```bash
docker-compose up -d
```

### Option 2: Helper Script

```bash
./docker/run.sh
```

### Option 3: Docker Run Command

```bash
docker run -d \
  --name residualrisk_app \
  -p 127.0.0.1:8501:8501 \
  --restart unless-stopped \
  residualrisk:latest
```

### With Custom Port

```bash
PORT=8080 ./docker/run.sh
```

Or:
```bash
docker run -d -p 127.0.0.1:8080:8501 residualrisk:latest
```

## Configuration

### Environment Variables

Configure the container using environment variables:

```bash
docker run -d \
  -p 127.0.0.1:8501:8501 \
  -e STREAMLIT_SERVER_MAX_UPLOAD_SIZE=500 \
  -e STREAMLIT_SERVER_ENABLE_CORS=true \
  residualrisk:latest
```

Available variables:
- `STREAMLIT_SERVER_PORT` (default: 8501)
- `STREAMLIT_SERVER_ADDRESS` (default: 0.0.0.0)
- `STREAMLIT_SERVER_MAX_UPLOAD_SIZE` (default: 200 MB)
- `STREAMLIT_SERVER_ENABLE_CORS` (default: false)
- `STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION` (default: true)

### Custom Streamlit Configuration

Mount a custom config file:

```bash
docker run -d \
  -p 127.0.0.1:8501:8501 \
  -v $(pwd)/custom-config.toml:/app/.streamlit/config.toml:ro \
  residualrisk:latest
```

### Resource Limits

**Important:** By default, containers have access to ALL host CPUs for maximum performance.

Limit CPU and memory if needed:

```bash
docker run -d \
  -p 127.0.0.1:8501:8501 \
  --cpus="8" \
  --memory="8g" \
  residualrisk:latest

# Or no limit (default)
docker run -d -p 127.0.0.1:8501:8501 residualrisk:latest
```

**Warning:** Limiting CPUs reduces multiprocessing performance. See [CPU-PERFORMANCE.md](CPU-PERFORMANCE.md) for details.

## Container Management

### View Logs

```bash
# Follow logs
docker logs -f residualrisk_app

# Last 100 lines
docker logs --tail 100 residualrisk_app

# With docker-compose
docker-compose logs -f
```

### Health Checks

The container includes a health check that monitors the Streamlit service:

```bash
# Check health status
docker inspect --format='{{.State.Health.Status}}' residualrisk_app

# View health check logs
docker inspect --format='{{json .State.Health}}' residualrisk_app | jq
```

### Stop and Remove

```bash
# Stop
docker stop residualrisk_app

# Remove
docker rm residualrisk_app

# Stop and remove
docker rm -f residualrisk_app

# With docker-compose
docker-compose down
```

### Restart

```bash
docker restart residualrisk_app

# Or with docker-compose
docker-compose restart
```

## Production Deployment

### SSL/HTTPS via Reverse Proxy (Required)

The Docker container does **not** handle SSL/TLS. For production, run a reverse proxy in front of the container that terminates SSL. The container binds only to `127.0.0.1:8501` so it is not directly reachable from the network.

A reference nginx configuration is provided at `docker/nginx/conf.d/app.conf`. Key requirements for the reverse proxy:

- **WebSocket proxying** — Streamlit requires `Upgrade` / `Connection: upgrade` headers
- **Long timeouts** — simulations can run for several minutes; set `proxy_read_timeout` to at least 300 seconds
- **Disable buffering** — set `proxy_buffering off` to avoid stalling the streaming UI

#### nginx example

```nginx
upstream streamlit {
    server 127.0.0.1:8501;
}

server {
    listen 80;
    server_name example.com;
    return 301 https://$host$request_uri;
}

server {
    listen 443 ssl http2;
    server_name example.com;

    ssl_certificate     /path/to/fullchain.pem;
    ssl_certificate_key /path/to/privkey.pem;

    location / {
        proxy_pass http://streamlit;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 300;
        proxy_send_timeout 300;
        proxy_buffering off;
        proxy_cache off;
    }
}
```

#### Caddy example

Caddy is the simplest option: it automatically obtains and renews SSL certificates via ACME (Let's Encrypt / ZeroSSL), handles HTTP→HTTPS redirection, and proxies WebSocket connections out of the box — no additional configuration required.

A reference `Caddyfile` is provided at `docker/Caddyfile`.

```
residualrisk.yourdomain.com {
    reverse_proxy 127.0.0.1:8501 {
        transport http {
            # Simulations can take several minutes — keep connections open
            read_timeout  5m
            write_timeout 5m
        }
    }
}
```

**If running Caddy as a systemd service** (recommended for Linux production deployments), copy the Caddyfile to your Caddy configuration directory and reload the service:

```bash
sudo cp docker/Caddyfile /etc/caddy/Caddyfile
sudo systemctl reload caddy
```

If Caddy is not yet running:

```bash
sudo systemctl enable --now caddy
```

**If running Caddy manually**, save as `Caddyfile` and run:

```bash
caddy run --config Caddyfile
# or to run in the background:
caddy start --config Caddyfile
```

Caddy handles SSL certificate acquisition and renewal automatically. No additional SSL configuration is needed.

#### Traefik example

```yaml
labels:
  - "traefik.enable=true"
  - "traefik.http.routers.residualrisk.rule=Host(`example.com`)"
  - "traefik.http.routers.residualrisk.entrypoints=websecure"
  - "traefik.http.routers.residualrisk.tls.certresolver=letsencrypt"
  - "traefik.http.services.residualrisk.loadbalancer.server.port=8501"
```

### Security Considerations

1. **Localhost-only binding**: The container exposes port 8501 only on `127.0.0.1` — not reachable from the network except through the reverse proxy
2. **Non-root user**: The container runs as `appuser` (UID 1000)
3. **Read-only filesystem**: Enabled in docker-compose.yml with tmpfs for temporary writes
4. **Resource limits**: Set appropriate CPU and memory limits for your environment
5. **Regular updates**: Rebuild with latest base images to pick up security patches

### Persistent Storage (Optional)

```bash
docker run -d \
  -p 127.0.0.1:8501:8501 \
  -v residualrisk_data:/app/data \
  residualrisk:latest
```

## Troubleshooting

### Container Won't Start

Check logs:
```bash
docker logs residualrisk_app
```

Common issues:
- Port 8501 already in use: Change with `-p 127.0.0.1:8080:8501`
- Permission issues: Check file ownership
- Resource constraints: Increase Docker resource limits

### Go Binary Not Found

The Go binary is compiled during the Docker build. If you see errors about missing `riskdays_go`:

1. Rebuild the image: `./docker/build.sh --load`
2. Check build logs for Go compilation errors
3. Verify platform matches: `docker inspect residualrisk:latest | grep Architecture`

### Performance Issues

1. **Increase CPU allocation**: `--cpus="8"`
2. **Increase memory**: `--memory="8g"`
3. **Verify Go binary is being used**: Check webapp UI shows "Go" as implementation
4. **Check system load**: `docker stats residualrisk_app`

### Health Check Failing

```bash
# Check if Streamlit is responding
docker exec residualrisk_app curl -f http://localhost:8501/_stcore/health

# Check if port is accessible
curl http://localhost:8501/_stcore/health
```

## Building from Scratch

If you need to customize the Dockerfile:

```bash
# Clean build (no cache)
docker buildx build --no-cache --platform linux/amd64 --tag residualrisk:latest --load .

# View build logs
docker buildx build --progress=plain --platform linux/amd64 --tag residualrisk:latest --load .
```

## Image Information

### Size

- Final image: ~500-700 MB (varies by architecture)
- Go binary: ~3-5 MB
- Python + dependencies: ~400-500 MB

### Layers

The multi-stage build minimizes layers:
1. Go builder stage (discarded)
2. Python base image
3. System dependencies
4. Python packages (uv)
5. Application code

### Verify Multi-Arch Support

```bash
# Inspect manifest
docker buildx imagetools inspect residualrisk:latest

# Should show both linux/amd64 and linux/arm64
```

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Build and Push Docker Image

on:
  push:
    branches: [main]
    tags: ['v*']

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v2

      - name: Login to GitHub Container Registry
        uses: docker/login-action@v2
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Build and push
        uses: docker/build-push-action@v4
        with:
          context: .
          platforms: linux/amd64,linux/arm64
          push: true
          tags: ghcr.io/${{ github.repository }}:latest
```

### GitLab CI Example

```yaml
build:
  image: docker:latest
  services:
    - docker:dind
  before_script:
    - docker login -u $CI_REGISTRY_USER -p $CI_REGISTRY_PASSWORD $CI_REGISTRY
  script:
    - docker buildx create --use
    - docker buildx build --platform linux/amd64,linux/arm64 --push -t $CI_REGISTRY_IMAGE:latest .
```

## License

This Docker configuration is part of the Residual Risk Estimator and is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).
