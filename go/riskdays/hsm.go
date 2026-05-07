// Residual HIV Transfusion Transmission Risk Estimation Tool
// Copyright (C) 2025-2026  Vitalant and Eduard Grebe Consulting
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
	"sort"
)

// HalfSampleMode estimates the mode of a continuous distribution using the
// Half-Sample Mode (HSM) algorithm (Bickel & Frühwirth, 2006; Robertson &
// Cryer, 1974).
//
// The algorithm iteratively finds the shortest interval ("shortest half")
// containing ⌈n/2⌉ of the sorted data points, then recurses on that interval
// until 1–3 points remain.
//
// Properties:
//   - Bandwidth-free: no kernel width or bin size to choose
//   - Robust to outliers: heavy tails are pruned away iteratively
//   - O(N log N) due to sorting; the recursion itself is O(N) per level with
//     O(log N) levels
//   - Not suitable for U-shaped or J-shaped distributions
//
// Parameters:
//   - data: positive-valued samples (will not be modified)
//
// Returns the estimated mode, or 0 if data is empty.
func HalfSampleMode(data []float64) float64 {
	n := len(data)
	if n == 0 {
		return 0
	}
	if n == 1 {
		return data[0]
	}

	// Sort a copy to avoid mutating the caller's slice.
	sorted := make([]float64, n)
	copy(sorted, data)
	sort.Float64s(sorted)

	return hsmRecurse(sorted)
}

// hsmRecurse operates on an already-sorted slice.
func hsmRecurse(sorted []float64) float64 {
	n := len(sorted)
	if n <= 2 {
		// Base case: return midpoint of the 1 or 2 remaining values.
		if n == 1 {
			return sorted[0]
		}
		return (sorted[0] + sorted[1]) / 2.0
	}
	if n == 3 {
		// Three points: pick the pair with the smallest gap, return midpoint.
		d01 := sorted[1] - sorted[0]
		d12 := sorted[2] - sorted[1]
		if d01 <= d12 {
			return (sorted[0] + sorted[1]) / 2.0
		}
		return (sorted[1] + sorted[2]) / 2.0
	}

	// Half-width: number of points in the shortest-half window.
	halfN := (n + 1) / 2 // ⌈n/2⌉

	// Find the window of halfN consecutive sorted values with the smallest range.
	bestStart := 0
	bestWidth := sorted[halfN-1] - sorted[0]
	for i := 1; i+halfN-1 < n; i++ {
		w := sorted[i+halfN-1] - sorted[i]
		if w < bestWidth {
			bestWidth = w
			bestStart = i
		}
	}

	// Recurse on the shortest-half window.
	return hsmRecurse(sorted[bestStart : bestStart+halfN])
}
