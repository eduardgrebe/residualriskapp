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
	"sync"
)

// KDEModeLog estimates the mode of a positive, right-skewed distribution
// via KDE on the log scale.  It matches the Python _kde_mode_log() exactly.
//
//   - data:  positive-valued samples
//   - nGrid: number of log-spaced grid points (default 100_000)
//   - cap:   maximum number of samples to use; 0 means no cap
//
// Algorithm:
//   1.  log-transform data
//   2.  Silverman bandwidth on log scale
//   3.  Gaussian KDE evaluated on a log-spaced grid
//   4.  Change-of-variables: f(k) = f_logk(log k) / k
//   5.  Return the grid value at maximum density
func KDEModeLog(data []float64, nGrid int, cap int) float64 {
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

	// ---- log-spaced grid (natural log) ----
	grid := make([]float64, nGrid)
	step := (maxVal - minVal) / float64(nGrid-1)
	for i := 0; i < nGrid; i++ {
		grid[i] = math.Exp(minVal + step*float64(i))
	}

	// ---- evaluate KDE on grid in parallel ----
	normFactor := float64(n) * bw * math.Sqrt(2*math.Pi)

	density := make([]float64, nGrid)
	var wg sync.WaitGroup
	// Split grid into chunks for goroutines
	workers := runtime.NumCPU()
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
			for i := s; i < e; i++ {
				logGrid := math.Log(grid[i])
				var sumK float64
				for _, lv := range logData {
					z := (logGrid - lv) / bw
					sumK += math.Exp(-0.5 * z * z)
				}
				// Density on log scale, then change-of-variables to original
				density[i] = (sumK / normFactor) / grid[i]
			}
		}(start, end)
	}
	wg.Wait()

	// ---- argmax ----
	bestIdx := 0
	bestVal := density[0]
	for i := 1; i < nGrid; i++ {
		if density[i] > bestVal {
			bestVal = density[i]
			bestIdx = i
		}
	}
	return grid[bestIdx]
}
