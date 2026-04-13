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
	"sort"

	"gonum.org/v1/gonum/stat"
)

// Concentration calculates viral concentration using exponential growth model
// C(t) = C0 * 2^(t / doubling_time)
func Concentration(C0, doublingTime, t float64) float64 {
	return C0 * math.Pow(2, t/doublingTime)
}

// ProbInfectiousCopies calculates probability of infection using dose-response model
// Standard single-hit exponential dose-response: prob = 1 - exp(-k * n_copies)
// Result is clamped to [0, 1] to guard against floating-point edge cases.
func ProbInfectiousCopies(nCopies, k float64) float64 {
	prob := 1.0 - math.Exp(-k*nCopies)
	return math.Max(0.0, math.Min(1.0, prob))
}

// Quantile calculates the quantile (percentile) of a sorted slice
func Quantile(data []float64, q float64) float64 {
	sorted := make([]float64, len(data))
	copy(sorted, data)
	sort.Float64s(sorted)
	return stat.Quantile(q, stat.Empirical, sorted, nil)
}

// Min returns the minimum value in a slice
func Min(data []float64) float64 {
	if len(data) == 0 {
		return 0
	}
	min := data[0]
	for _, v := range data[1:] {
		if v < min {
			min = v
		}
	}
	return min
}

// Max returns the maximum value in a slice
func Max(data []float64) float64 {
	if len(data) == 0 {
		return 0
	}
	max := data[0]
	for _, v := range data[1:] {
		if v > max {
			max = v
		}
	}
	return max
}

// Mean calculates the mean of a slice
func Mean(data []float64) float64 {
	if len(data) == 0 {
		return 0
	}
	sum := 0.0
	for _, v := range data {
		sum += v
	}
	return sum / float64(len(data))
}

// Median calculates the median of a slice
func Median(data []float64) float64 {
	return Quantile(data, 0.5)
}

// ModeRounded calculates the mode of data rounded to specified precision
// Matches scipy.stats.mode(np.array(list).round(precision)).mode
func ModeRounded(data []float64, precision int) float64 {
	if len(data) == 0 {
		return 0
	}

	// Round all values to specified precision
	rounded := make([]float64, len(data))
	multiplier := math.Pow(10, float64(precision))
	for i, v := range data {
		rounded[i] = math.Round(v*multiplier) / multiplier
	}

	// Count frequencies
	freq := make(map[float64]int)
	for _, v := range rounded {
		freq[v]++
	}

	// Find mode (most frequent value)
	maxCount := 0
	var mode float64
	for v, count := range freq {
		if count > maxCount {
			maxCount = count
			mode = v
		}
	}

	return mode
}

// Round rounds a value to specified decimal places
func Round(value float64, precision int) float64 {
	multiplier := math.Pow(10, float64(precision))
	return math.Round(value*multiplier) / multiplier
}
