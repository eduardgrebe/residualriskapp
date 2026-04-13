# Residual HIV Transfusion Transmission Risk Estimation Tool
# Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
# Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
#
# Multi-architecture Docker image for the Streamlit web application
# Supports: linux/amd64, linux/arm64

# =============================================================================
# Stage 1: Build Go binary for high-performance simulations
# =============================================================================
FROM --platform=${BUILDPLATFORM} golang:1.26-alpine AS go-builder

# Build arguments automatically set by Docker buildx
ARG TARGETOS
ARG TARGETARCH

# Install build dependencies
RUN apk add --no-cache git make

# Set working directory
WORKDIR /build

# Copy Go module files
COPY go/go.mod go/go.sum ./

# Download dependencies
RUN go mod download

# Copy Go source code
COPY go/ ./

# Build the binary for target architecture
# CGO_ENABLED=0 creates a static binary (no C dependencies)
RUN CGO_ENABLED=0 GOOS=${TARGETOS} GOARCH=${TARGETARCH} \
    go build -ldflags="-w -s" -o riskdays_go main.go

# Verify the binary was built and is executable
RUN ls -lh riskdays_go && \
    test -x riskdays_go && \
    echo "✓ Go binary built successfully: $(ls -lh riskdays_go | awk '{print $5}')"

# =============================================================================
# Stage 2: Runtime image with Python and Streamlit
# =============================================================================
FROM python:3.14-slim

# Metadata
LABEL org.opencontainers.image.title="Residual Risk Estimator"
LABEL org.opencontainers.image.description="HIV transfusion transmission risk estimation tool"
LABEL org.opencontainers.image.vendor="Vitalant Research Institute"
LABEL org.opencontainers.image.licenses="AGPL-3.0"

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    # Streamlit configuration
    STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_ENABLE_CORS=false \
    STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION=true \
    # Application
    APP_HOME=/app

# Install system dependencies (curl needed for health checks)
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install uv for fast Python package management using pip
RUN pip install --no-cache-dir uv && \
    uv --version

# Create non-root user for security
RUN useradd -m -u 1000 -s /bin/bash appuser && \
    mkdir -p ${APP_HOME} && \
    chown -R appuser:appuser ${APP_HOME}

# Set working directory
WORKDIR ${APP_HOME}

# Copy Python dependency files
COPY --chown=appuser:appuser pyproject.toml uv.lock ./

# Install Python dependencies using uv
# This creates a virtual environment at .venv
RUN uv sync --frozen

# Copy compiled Go binary from builder stage
COPY --from=go-builder --chown=appuser:appuser /build/riskdays_go ./go/bin/riskdays_go
RUN chmod +x ./go/bin/riskdays_go

# Copy application code
COPY --chown=appuser:appuser *.py ./
COPY --chown=appuser:appuser pages/ ./pages/
COPY --chown=appuser:appuser docs/ ./docs/
COPY --chown=appuser:appuser static/ ./static/
COPY --chown=appuser:appuser LICENSE ./

# Copy Streamlit configuration (optional - will use defaults if not present)
COPY --chown=appuser:appuser .streamlit/ ./.streamlit/

# Switch to non-root user
USER appuser

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Run Streamlit application using uv's virtual environment
CMD [".venv/bin/streamlit", "run", "app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--server.fileWatcherType=none", \
     "--browser.gatherUsageStats=false"]
