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

	"gonum.org/v1/gonum/dsp/fourier"
)

// KDEModeLog estimates the mode of a positive, right-skewed distribution
// via KDE on the log scale.
//
//   - data:    positive-valued samples
//   - nGrid:   number of log-spaced grid points; 0 → auto (max 100k, ≤200k)
//   - cap:     maximum number of samples to use; 0 means no cap
//   - threads: unused (kept for API compatibility); FFT convolution is single-pass
//
// Algorithm (FFT-accelerated, O(N + M log M)):
//  1. Log-transform data
//  2. Silverman bandwidth on log scale
//  3. Linear binning of log-data onto a regular grid (O(N))
//  4. Gaussian kernel with wraparound for circular convolution (O(M))
//  5. FFT convolution: IFFT(FFT(counts) ⊙ FFT(kernel)) (O(M log M))
//  6. Change-of-variables: f(k) = f_logk(log k) / k
//  7. Return the grid value at maximum density
func KDEModeLog(data []float64, nGrid int, cap int, threads int) float64 {
	_ = threads

	if len(data) == 0 {
		return 0
	}

	// ---- cap large inputs ----
	work := data
	if cap > 0 && len(data) > cap {
		work = make([]float64, cap)
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
	std := math.Sqrt(ssq / float64(n - 1))

	// ---- Silverman bandwidth ----
	if std <= 0 || n <= 1 {
		return work[0]
	}
	bw := 1.06 * std * math.Pow(float64(n), -0.2)

	// ---- regular grid in log space ----
	step := (maxVal - minVal) / float64(nGrid-1)

	// ---- Step 1: Linear binning of logData onto the grid (O(N)) ----
	// Each data point distributes its unit mass to the two nearest grid points
	// proportionally to the distance.
	counts := make([]float64, nGrid)
	invStep := 1.0 / step
	for _, lv := range logData {
		pos := (lv - minVal) * invStep
		j := int(pos)
		if j < 0 {
			counts[0] += 1.0
		} else if j >= nGrid-1 {
			counts[nGrid-1] += 1.0
		} else {
			frac := pos - float64(j)
			counts[j] += 1.0 - frac
			counts[j+1] += frac
		}
	}

	// ---- Step 2: Gaussian kernel (symmetric, wrapped for circular convolution) ----
	// kernel[d] = exp(-0.5 * (d * step / bw)^2) for the positive half,
	// mirrored into the upper half of the array for the negative offsets
	// (circular convolution wraparound).
	kernel := make([]float64, nGrid)
	stepOverBw := step / bw
	halfM := nGrid / 2
	for d := 0; d <= halfM; d++ {
		z := float64(d) * stepOverBw
		if z > 5.0 {
			break // Gaussian < 3.7e-6 beyond 5σ
		}
		kval := math.Exp(-0.5 * z * z)
		kernel[d] = kval
		if d > 0 && nGrid-d > halfM {
			kernel[nGrid-d] = kval // wrap negative offsets
		}
	}

	// ---- Step 3: FFT convolution (O(M log M)) ----
	fftObj := fourier.NewFFT(nGrid)
	countsFFT := fftObj.Coefficients(nil, counts)
	kernelFFT := fftObj.Coefficients(nil, kernel)
	for i := range countsFFT {
		countsFFT[i] *= kernelFFT[i]
	}
	conv := fftObj.Sequence(nil, countsFFT)

	// ---- Step 4: Normalize and change-of-variables, then argmax ----
	// f_log(g_i) = conv[i] / (N * bw * sqrt(2π))
	// f(k_i) = f_log(g_i) / exp(g_i)
	invNormFactor := 1.0 / (float64(n) * bw * math.Sqrt(2*math.Pi))
	bestIdx := 0
	bestVal := -math.MaxFloat64
	for i := 0; i < nGrid; i++ {
		logGridVal := minVal + step*float64(i)
		d := conv[i] * invNormFactor / math.Exp(logGridVal)
		if d > bestVal {
			bestVal = d
			bestIdx = i
		}
	}

	return math.Exp(minVal + step*float64(bestIdx))
}
