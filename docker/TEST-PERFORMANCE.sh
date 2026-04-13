#!/bin/bash
# Quick performance test for Docker container

set -e

echo "=== Residual Risk Estimator - Performance Test ==="
echo ""

# Check if container is running
if ! docker ps | grep -q residualrisk_app; then
    echo "Error: Container 'residualrisk_app' is not running"
    echo "Start it with: docker-compose up -d"
    exit 1
fi

echo "1. Checking CPU access..."
CONTAINER_CPUS=$(docker exec residualrisk_app python -c "import multiprocessing; print(multiprocessing.cpu_count())")
HOST_CPUS=$(python3 -c "import multiprocessing; print(multiprocessing.cpu_count())" 2>/dev/null || echo "unknown")

echo "   Container CPUs: $CONTAINER_CPUS"
echo "   Host CPUs:      $HOST_CPUS"

if [ "$CONTAINER_CPUS" == "$HOST_CPUS" ]; then
    echo "   ✓ Container has full CPU access"
elif [ "$CONTAINER_CPUS" -lt "$HOST_CPUS" ]; then
    echo "   ⚠ Container has limited CPU access (performance will be reduced)"
    echo "   Check docker-compose.yml for CPU limits"
else
    echo "   ✓ CPU detection working"
fi

echo ""
echo "2. Checking Go binary..."
if docker exec residualrisk_app test -f go/bin/riskdays_go; then
    GO_SIZE=$(docker exec residualrisk_app ls -lh go/bin/riskdays_go | awk '{print $5}')
    echo "   ✓ Go binary exists (${GO_SIZE})"
else
    echo "   ✗ Go binary missing - performance will be reduced"
    echo "   Rebuild with: ./docker/build.sh --load"
fi

echo ""
echo "3. Current resource usage..."
docker stats residualrisk_app --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}"

echo ""
echo "=== Performance Test Complete ==="
echo ""
echo "Expected behavior during simulation:"
echo "  - CPU % should spike to ${CONTAINER_CPUS}00%+ (when using all cores)"
echo "  - Monitor with: docker stats residualrisk_app"
echo ""
echo "For detailed performance information, see:"
echo "  docker/CPU-PERFORMANCE.md"
