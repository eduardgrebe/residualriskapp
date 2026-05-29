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
)

// SinVaried computes the sinusoidal set-point oscillation.
// Returns offset + a*sin(b*t).
// Corresponds to Python _sin_varied().
func SinVaried(t, a, b, offset float64) float64 {
	return offset + a*math.Sin(b*t)
}

// FindTcrit computes the time at which exponential viral growth first reaches
// the set-point. This is the analytic solution to:
//
//	C0 * 2^((t - eclipse) / doubling_time) = set_point
//
// Solving: tcrit = eclipse + doubling_time * log2(set_point / C0)
//
// The Python implementation uses a grid search (np.arange(0, 265, 0.1)) which
// gives ~0.01-day resolution. The analytic formula is exact and faster.
func FindTcrit(eclipse, C0, doublingTime, setPoint float64) float64 {
	return eclipse + doublingTime*math.Log2(setPoint/C0)
}

// VLPostBT computes the PrEP breakthrough viral load at time t.
// Three phases:
//  1. t < eclipse: VL = 0 (eclipse phase)
//  2. eclipse <= t <= tcrit: VL = C0 * 2^((t-eclipse)/doubling_time) (exponential growth)
//  3. t > tcrit: VL = set_point * SinVaried(t-tcrit, a, b, offset) (oscillating plateau)
//
// tcrit must be precomputed via FindTcrit().
// Corresponds to Python _vl_postbt().
func VLPostBT(t, eclipse, C0, doublingTime, setPoint, a, b, offset, tcrit float64) float64 {
	if t < eclipse {
		return 0.0
	}
	if t <= tcrit {
		return C0 * math.Pow(2, (t-eclipse)/doublingTime)
	}
	// Modelled viral load can dip below zero when the sinusoidal set-point
	// oscillation amplitude exceeds its offset (a > offset); clamp to a
	// physical floor of zero.
	return math.Max(0.0, setPoint*SinVaried(t-tcrit, a, b, offset))
}

// DrugEffectFactor returns the antiretroviral transmissibility-reduction factor
// at time t — the multiplier applied to the per-time infection probability in
// ProbInfectiousPrep. drugEffect is a scalar in (0, 1] (1.0 = no reduction).
//
// It is currently CONSTANT in t, so it factors straight out of the RDE integral
// (multiplying here is numerically identical to scaling the final RDE), and the
// default 1.0 leaves the integrand bit-for-bit unchanged.
//
// t is taken deliberately as a placeholder for a future time-varying drug
// effect: breakthrough infections on long-acting injectable PrEP typically
// occur as the drug washes out, so the factor should relax toward 1.0 across
// the window as drug concentration decays. Returning a function of t here (e.g.
// an exponential wash-out from the last-injection time) is then the only
// correct placement, since it would no longer factor out of the integral.
// Corresponds to Python _drug_effect().
func DrugEffectFactor(t, drugEffect float64) float64 {
	// Placeholder: constant in t. Replace with a t-dependent expression
	// (e.g. decaying long-acting drug concentration) to model PrEP wash-out;
	// expected to matter most for injectable PrEP.
	return drugEffect
}

// ProbInfectiousPrep calculates the probability that a transfusion is infectious
// given PrEP breakthrough viral dynamics.
// Corresponds to Python _prob_infectious_prep().
func ProbInfectiousPrep(t float64, params PrepInnerParams) float64 {
	tcrit := FindTcrit(params.Eclipse, params.C0, params.DoublingTime, params.SetPoint)
	C := VLPostBT(t, params.Eclipse, params.C0, params.DoublingTime,
		params.SetPoint, params.A, params.B, params.Offset, tcrit)
	nCopies := C * float64(params.CopiesPerVirion) * params.VolumeTransfused
	// Drug effect is a linear scalar on the realized infection probability.
	return DrugEffectFactor(t, params.DrugEffect) * ProbInfectiousCopies(nCopies, params.K)
}

// ProbNondetectionSerology computes the probability that serology (antibody)
// testing fails to detect infection at time t.
//
// Returns:
//   - 1.0 if t < min (before seroconversion window)
//   - 0.0 if t > max (after max detection time)
//   - exp(-((t - min) / alpha)^beta) otherwise (Weibull-like decay)
//
// Corresponds to Python _prob_nondetection_serology_prep().
func ProbNondetectionSerology(t, min, max, alpha, beta float64) float64 {
	if t < min {
		return 1.0
	}
	if t > max {
		return 0.0
	}
	return math.Exp(-math.Pow((t-min)/alpha, beta))
}

// ProbNondetectionPrep calculates the probability that NAT testing fails to
// detect infection given PrEP breakthrough viral dynamics.
// Reuses ProbPosInit and ProbNegRetest from probability.go.
// Corresponds to Python _prob_nondetection_prep().
func ProbNondetectionPrep(t float64, params PrepInnerParams) float64 {
	tcrit := FindTcrit(params.Eclipse, params.C0, params.DoublingTime, params.SetPoint)
	Cv := VLPostBT(t, params.Eclipse, params.C0, params.DoublingTime,
		params.SetPoint, params.A, params.B, params.Offset, tcrit)
	Cc := float64(params.CopiesPerVirion) * Cv

	if Cc == 0 {
		return 1.0
	}

	pPosInit, _ := ProbPosInit(Cc, params.DoublingTime, params.PoolSize,
		params.LOD50, params.LOD95LOD50Ratio, params.Z)
	pNegRetest, _ := ProbNegRetest(Cc, params.DoublingTime, params.PoolSize,
		params.LOD50, params.LOD95LOD50Ratio, params.Retests, params.Z)

	return 1 - pPosInit*(1-pNegRetest)
}

// ProbInfectiousNondetectionPrep computes the integrand for the PrEP risk days
// calculation: the product of infectivity, NAT non-detection, and serology
// non-detection probabilities.
//
// This is a 3-product integrand (vs 2-product for baseline):
//
//	P_infectious × P_NAT_nondetection × P_serology_nondetection
//
// Corresponds to Python _prob_infectious_nondetection_prep().
func ProbInfectiousNondetectionPrep(t float64, params PrepInnerParams) float64 {
	// Precompute tcrit once for all three components
	tcrit := FindTcrit(params.Eclipse, params.C0, params.DoublingTime, params.SetPoint)

	// Viral load at time t
	C := VLPostBT(t, params.Eclipse, params.C0, params.DoublingTime,
		params.SetPoint, params.A, params.B, params.Offset, tcrit)

	// Infectivity (drug effect is a linear scalar on the infection probability;
	// see DrugEffectFactor — constant in t today, placeholder for wash-out)
	nCopies := C * float64(params.CopiesPerVirion) * params.VolumeTransfused
	pInfectious := DrugEffectFactor(t, params.DrugEffect) * ProbInfectiousCopies(nCopies, params.K)

	// NAT non-detection
	Cc := float64(params.CopiesPerVirion) * C
	var pNATNondet float64
	if Cc == 0 {
		pNATNondet = 1.0
	} else {
		pPosInit, _ := ProbPosInit(Cc, params.DoublingTime, params.PoolSize,
			params.LOD50, params.LOD95LOD50Ratio, params.Z)
		pNegRetest, _ := ProbNegRetest(Cc, params.DoublingTime, params.PoolSize,
			params.LOD50, params.LOD95LOD50Ratio, params.Retests, params.Z)
		pNATNondet = 1 - pPosInit*(1-pNegRetest)
	}

	// Serology non-detection
	pSeroNondet := ProbNondetectionSerology(t, params.SerMin, params.SerMax,
		params.SerAlpha, params.SerBeta)

	return pInfectious * pNATNondet * pSeroNondet
}
