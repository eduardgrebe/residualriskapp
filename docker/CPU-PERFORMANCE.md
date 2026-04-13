# CPU Configuration and Performance

The Residual Risk Estimator uses multiprocessing for high-performance simulations. Proper CPU configuration is critical for performance.

## Default Configuration

**The Docker container has access to ALL host CPU cores by default.**

This ensures maximum performance for the multiprocessing architecture used by both the Python and Go implementations.

## CPU Detection

The application automatically detects available CPUs:

**Python implementation:**
```python
import multiprocessing
cores = multiprocessing.cpu_count()  # Detects all available cores
threads = cores - 1                   # Default: uses n-1 cores
```

**Go implementation:**
```go
import "runtime"
cores := runtime.NumCPU()            // Detects all available cores
```

**In the Streamlit UI:**
- Sidebar shows detected CPU count
- User can select number of threads to use (default: all cores - 1)

## Testing CPU Access

Verify the container can see all CPUs:

```bash
# Start container
docker-compose up -d

# Check CPU count inside container
docker exec residualrisk_app python -c "import multiprocessing; print(f'CPUs: {multiprocessing.cpu_count()}')"

# Compare with host
python3 -c "import multiprocessing; print(f'Host CPUs: {multiprocessing.cpu_count()}')"
```

## Limiting CPU Resources (Optional)

If you need to limit CPU usage (e.g., on shared servers), uncomment the `deploy` section in `docker-compose.yml`:

```yaml
deploy:
  resources:
    limits:
      cpus: '16'        # Maximum CPUs (use '0' for unlimited)
      memory: 16G       # Maximum memory
    reservations:
      cpus: '4'         # Minimum guaranteed
      memory: 4G
```

**Warning:** Limiting CPUs will reduce performance proportionally.

### Using docker run

```bash
# Limit to 8 CPUs
docker run -d --cpus="8" -p 127.0.0.1:8501:8501 residualrisk:latest

# Limit to 50% of host CPUs
docker run -d --cpus="0.5" -p 127.0.0.1:8501:8501 residualrisk:latest

# No limit (default)
docker run -d -p 127.0.0.1:8501:8501 residualrisk:latest
```

## Performance Expectations

Performance scales approximately linearly with CPU cores:

| CPUs | Simulations/sec | 250K sims | 1M sims |
|------|-----------------|-----------|---------|
| 1    | ~50-100         | ~42 min   | ~3 hrs  |
| 4    | ~200-400        | ~11 min   | ~42 min |
| 8    | ~400-800        | ~5 min    | ~21 min |
| 16   | ~800-1600       | ~3 min    | ~10 min |
| 32   | ~1600-3200      | ~1.5 min  | ~5 min  |

*Note: Actual performance varies by CPU architecture. Times shown are for Go implementation.*

### Python vs Go Performance

- **Python**: ~10-50 simulations/sec (single-threaded portions limit scaling)
- **Go**: ~100-500 simulations/sec (scales well with cores)
- **Speedup**: 10-50× faster with Go

## Monitoring CPU Usage

### Real-time monitoring:

```bash
# CPU usage of running container
docker stats residualrisk_app

# Detailed resource usage
docker exec residualrisk_app top -b -n 1
```

### Check if multiprocessing is working:

```bash
# Watch CPU usage while running simulation
docker stats residualrisk_app --no-stream

# Should show CPU% > 100% if using multiple cores
# e.g., 800% = using 8 cores at 100%
```

## Troubleshooting

### Container shows fewer CPUs than host

**Possible causes:**

1. **CPU limits set in docker-compose.yml**
   - Check: `grep -A 5 "resources:" docker-compose.yml`
   - Fix: Comment out or increase limits

2. **Docker Desktop resource limits** (Mac/Windows)
   - Check: Docker Desktop → Settings → Resources → CPUs
   - Fix: Increase allocated CPUs

3. **Container orchestration limits** (Kubernetes, Swarm)
   - Check: Deployment/pod resource limits
   - Fix: Update resource requests/limits

### Performance not scaling with CPUs

1. **Using Python implementation instead of Go**
   - Check: Streamlit UI shows "Python" or "Go"
   - Fix: Ensure Go binary exists: `docker exec residualrisk_app ls -lh go/bin/riskdays_go`
   - Rebuild if needed: `./docker/build.sh --load`

2. **Low number of simulations**
   - Multiprocessing overhead dominates for n_bs < 5,000
   - Use ≥ 25,000 simulations to see scaling benefits

3. **Thread count set too low in UI**
   - Check: Streamlit sidebar "Select number of CPU cores to use"
   - Should default to: CPU count - 1

## Production Recommendations

1. **No CPU limits** - Let container use all available cores
2. **Use Go implementation** - 10-50× faster than Python
3. **Monitor resources** - Use `docker stats` to verify multicore usage
4. **Memory scaling** - Allow ~100-200 MB per CPU core for large simulations

## Example Configurations

### Development (local, resource-constrained):
```yaml
deploy:
  resources:
    limits:
      cpus: '4'
      memory: 4G
```

### Production (dedicated server):
```yaml
# No limits - use all resources
# (deploy section commented out or removed)
```

### Shared environment (fair usage):
```yaml
deploy:
  resources:
    limits:
      cpus: '8'        # Fair share of 32-core server
      memory: 8G
    reservations:
      cpus: '2'        # Guaranteed minimum
      memory: 2G
```

## Verification Checklist

After starting the container, verify:

- [ ] Container can see all host CPUs
- [ ] Go binary is present and working
- [ ] Streamlit UI shows correct CPU count
- [ ] Simulations use multiple cores (check `docker stats`)
- [ ] Performance scales with CPU count

```bash
# Quick verification script
docker exec residualrisk_app python -c "
import multiprocessing
import os
print(f'CPUs detected: {multiprocessing.cpu_count()}')
print(f'Go binary exists: {os.path.exists(\"go/bin/riskdays_go\")}')
"

docker stats residualrisk_app --no-stream
```

Expected output:
```
CPUs detected: 16
Go binary exists: True

CONTAINER ID   CPU %    MEM USAGE / LIMIT
abc123...      0.5%     250MiB / 16GiB
```

When running simulations, CPU % should spike to 800%+ (for 8+ cores).
