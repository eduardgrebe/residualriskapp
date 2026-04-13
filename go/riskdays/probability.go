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
	"fmt"
	"math"
)

// ProbInfectious calculates probability of infection at time t
// Corresponds to Python _prob_infectious()
func ProbInfectious(t, C0, doublingTime, volumeTransfused, k float64, copiesPerVirion int) float64 {
	C := Concentration(C0, doublingTime, t)
	nCopies := C * float64(copiesPerVirion) * volumeTransfused
	return ProbInfectiousCopies(nCopies, k)
}

// ProbPosInit calculates probability of initial positive test
// Corresponds to Python _prob_pos_init()
func ProbPosInit(C, doublingTime float64, poolSize int, lod50, lod95LOD50Ratio, z float64) (float64, error) {
	if poolSize < 1 {
		return 0, fmt.Errorf("pool_size must be at least 1")
	}

	// X = z * log10(C / (pool_size * lod50)) / log10(lod95_lod50_ratio)
	ratio := C / (float64(poolSize) * lod50)
	if ratio <= 0 {
		return 0, nil
	}

	X := z * math.Log10(ratio) / math.Log10(lod95LOD50Ratio)

	// prob = norm.cdf(X)
	prob := NormalCDF(X)
	return prob, nil
}

// ProbNegRetest calculates probability all retests are negative
// Corresponds to Python _prob_neg_retest()
func ProbNegRetest(C, doublingTime float64, poolSize int, lod50, lod95LOD50Ratio float64, retests int, z float64) (float64, error) {
	if poolSize < 1 {
		return 0, fmt.Errorf("pool_size must be at least 1")
	}
	if retests < 0 {
		return 0, fmt.Errorf("retests must be non-negative")
	}
	if retests == 0 {
		return 0, nil
	}

	// X = z * log10(C / lod50) / log10(lod95_lod50_ratio)
	ratio := C / lod50
	if ratio <= 0 {
		return 1, nil // If concentration is 0 or negative, all retests are negative
	}

	X := z * math.Log10(ratio) / math.Log10(lod95LOD50Ratio)

	// prob = (1 - norm.cdf(X))^retests
	probNeg := 1 - NormalCDF(X)
	prob := math.Pow(probNeg, float64(retests))
	return prob, nil
}

// ProbNondetection calculates probability of non-detection by testing
// Corresponds to Python _prob_nondetection()
func ProbNondetection(t, C0, doublingTime float64, copiesPerVirion int, poolSize int,
	lod50, lod95LOD50Ratio float64, retests int, z float64) (float64, error) {

	Cv := Concentration(C0, doublingTime, t)
	Cc := float64(copiesPerVirion) * Cv

	pPosInit, err := ProbPosInit(Cc, doublingTime, poolSize, lod50, lod95LOD50Ratio, z)
	if err != nil {
		return 0, err
	}

	pNegRetest, err := ProbNegRetest(Cc, doublingTime, poolSize, lod50, lod95LOD50Ratio, retests, z)
	if err != nil {
		return 0, err
	}

	prob := 1 - pPosInit*(1-pNegRetest)
	return prob, nil
}

// ProbInfectiousNondetection calculates the product of infection and non-detection probabilities
// This is the integrand used in risk days calculation
// Corresponds to Python _prob_infectious_nondetection()
func ProbInfectiousNondetection(t float64, params RiskDaysInnerParams) (float64, error) {
	pInfectious := ProbInfectious(t, params.C0, params.DoublingTime, params.VolumeTransfused,
		params.K, params.CopiesPerVirion)

	pNondetection, err := ProbNondetection(t, params.C0, params.DoublingTime, params.CopiesPerVirion,
		params.PoolSize, params.LOD50, params.LOD95LOD50Ratio, params.Retests, params.Z)
	if err != nil {
		return 0, err
	}

	return pInfectious * pNondetection, nil
}
