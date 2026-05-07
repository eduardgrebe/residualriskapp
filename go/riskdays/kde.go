// Residual HIV Transfusion Transmission Risk Estimation Tool
// Copyright (C) 2025  Vitalant and Eduard Grebe Consulting
// Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

package riskdays

import (
	"math"
	"math/rand"
	"runtime"
	"sort"
	"sync"
)

// KDEModeLog estimates the mode of a positive, right-skewed distribution
// via KDE on the log scale.  It matches the Python _kde_mode_log() exactly.
//
//   - data:    positive-valued samples
//   - nGrid:   number of log-spaced grid points; 0 → auto (max 100k, ≤200k)
//   - cap:     maximum number of samples to use; 0 means no cap
//   - threads: number of parallel goroutines; ≤ 0 uses runtime.NumCPU()
//
// Algorithm:
//   1.  log-transform data
//   2.  Silverman bandwidth on log scale
//   3.  Gaussian KDE evaluated on a log-spaced grid
//   4.  Change-of-variables: f(k) = f_logk(log k) / k
//   5.  Return the grid value at maximum density
func KDEModeLog(data []float64, nGrid int, cap int, threads int) float64 {
	if len(data) == 0 {
		return 0
	}

	// ---- cap large inputs ----
	work := data
	if cap > 0 && len(data) > cap {
		work = make([]float64, cap)
		// Reservoir-style random subset
		rng := rand.New(rand.NewSource(42))
		for i := 0; i < cap; i++ {
			j := rng.Intn(len(data))
			work[i] = data[j]
		}
	}

	// ---- auto grid size when nGrid <= 0 ----
	if nGrid <= 0 {
		nGrid = len(work)
		if nGrid < 100_000 {
			nGrid = 100_000
		} else if nGrid > 200_000 {
			nGrid = 200_000
		}
	} else if nGrid < 100 {
		nGrid = 100
	}

	// ---- log transform & descriptive stats ----
	n := len(work)
	logData := make([]float64, n)
	minVal := math.MaxFloat64
	maxVal := -math.MaxFloat64
	sum := 0.0
	for i, v := range work {
		lv := math.Log(v)
		logData[i] = lv
		sum += lv
		if lv < minVal {
			minVal = lv
		}
		if lv > maxVal {
			maxVal = lv
		}
	}
	mean := sum / float64(n)

	// Standard deviation of log data
	var ssq float64
	for _, lv := range logData {
		d := lv - mean
		ssq += d * d
	}
	std := math.Sqrt(ssq / float64(n-1))

	// ---- Silverman bandwidth ----
	if std <= 0 || n <= 1 {
		// Single value or all identical: mode is the (only) value
		return work[0]
	}
	bw := 1.06 * std * math.Pow(float64(n), -0.2)

	// ---- log-spaced grid: work directly in log space ----
	// grid[i] = Exp(logGridVal[i]) so log(grid[i]) == logGridVal[i] — no need
	// to store the exp'd grid or call Log() inside the hot loop.
	step := (maxVal - minVal) / float64(nGrid-1)

	// ---- pre-sort logData so the 5σ cutoff can use a two-pointer scan ----
	sorted := make([]float64, n)
	copy(sorted, logData)
	sort.Float64s(sorted)

	// ---- evaluate KDE on grid in parallel ----
	// Pre-compute reciprocals to replace divisions with multiplications in the
	// inner loop (division is ~3–5× slower than multiplication on modern CPUs).
	invBw := 1.0 / bw
	invBwSqNegHalf := -0.5 * invBw * invBw
	normFactor := float64(n) * bw * math.Sqrt(2*math.Pi)
	invNormFactor := 1.0 / normFactor
	cutoff := 5.0 * bw // beyond ±5σ the Gaussian kernel contributes < 3.7e-6

	density := make([]float64, nGrid)
	var wg sync.WaitGroup
	workers := threads
	if workers <= 0 {
		workers = runtime.NumCPU()
	}
	chunkSize := (nGrid + workers - 1) / workers
	for c := 0; c < workers; c++ {
		start := c * chunkSize
		end := start + chunkSize
		if end > nGrid {
			end = nGrid
		}
		if start >= end {
			continue
		}
		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			// Two-pointer bounds into sorted logData for the 5σ window.
			// As logGridVal increases monotonically with i, lo/hi only move right.
			lo, hi := 0, 0
			for i := s; i < e; i++ {
				logGridVal := minVal + step*float64(i)
				loBound := logGridVal - cutoff
				hiBound := logGridVal + cutoff
				// Advance lower bound
				for lo < n && sorted[lo] < loBound {
					lo++
				}
				// Advance upper bound
				for hi < n && sorted[hi] <= hiBound {
					hi++
				}
				var sumK float64
				for j := lo; j < hi; j++ {
					diff := logGridVal - sorted[j]
					sumK += math.Exp(diff * diff * invBwSqNegHalf)
				}
				// Change-of-variables: f(k) = f_logk(log k) / k
				density[i] = (sumK * invNormFactor) / math.Exp(logGridVal)
			}
		}(start, end)
	}
	wg.Wait()

	// ---- argmax → return original-scale mode ----
	bestIdx := 0
	bestVal := density[0]
	for i := 1; i < nGrid; i++ {
		if density[i] > bestVal {
			bestVal = density[i]
			bestIdx = i
		}
	}
	return math.Exp(minVal + step*float64(bestIdx))
}
