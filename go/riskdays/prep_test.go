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
	"testing"
)

// Default PrEP test parameters matching Python defaults
var defaultPrepParams = PrepInnerParams{
	CopiesPerVirion:  2,
	C0:               0.00025,
	DoublingTime:     0.8542,
	SetPoint:         336,
	Eclipse:          7.0,
	A:                0.7,
	B:                0.6,
	Offset:           1.0,
	DrugEffect:       1.0,
	VolumeTransfused: 200.0,
	K:                0.000673,
	PoolSize:         16,
	LOD50:            2.73,
	LOD95LOD50Ratio:  3.5,
	Retests:          1,
	SerMin:           10,
	SerMax:           500,
	SerAlpha:         9.1,
	SerBeta:          5.2,
	Z:                1.6449,
	LimitMin:         -100,
	LimitMax:         500,
}

// productionPrepParams mirrors defaultPrepParams but with the production
// serology defaults used by risk_days_prep_bs / the _go.py bridge / the app.
// These give a much wider, slower-decaying active integration window (~[10, 169]
// days) than the narrow defaultPrepParams window (~[8.7, 22.5] days).
var productionPrepParams = func() PrepInnerParams {
	p := defaultPrepParams
	p.SerMin = 28.7
	p.SerMax = 250
	p.SerAlpha = 50.49434
	p.SerBeta = 1.15062
	return p
}()

// --- SinVaried ---

func TestSinVaried_Zero(t *testing.T) {
	// At t=0, sin(0)=0, so result = offset
	got := SinVaried(0, 0.7, 0.6, 1.0)
	if got != 1.0 {
		t.Errorf("SinVaried(0, 0.7, 0.6, 1.0) = %v, want 1.0", got)
	}
}

func TestSinVaried_Peak(t *testing.T) {
	// At t=π/(2*0.6), sin(0.6*t)=sin(π/2)=1, result = 1.0 + 0.7*1 = 1.7
	tVal := math.Pi / (2 * 0.6)
	got := SinVaried(tVal, 0.7, 0.6, 1.0)
	if math.Abs(got-1.7) > 1e-10 {
		t.Errorf("SinVaried(π/1.2, 0.7, 0.6, 1.0) = %v, want 1.7", got)
	}
}

// --- FindTcrit ---

func TestFindTcrit(t *testing.T) {
	// Analytic: eclipse + dt * log2((sp/copiesPerVirion) / C0)
	// = 7 + 0.8542 * log2((336/2) / 0.00025) = 7 + 0.8542 * 19.3355... ≈ 23.54
	tcrit := FindTcrit(7.0, 0.00025, 0.8542, 336, 2.0)
	expected := 23.5357 // From Python verification (336 copies/mL -> 168 virions/mL)
	if math.Abs(tcrit-expected) > 0.001 {
		t.Errorf("FindTcrit(7, 0.00025, 0.8542, 336) = %v, want ~%v", tcrit, expected)
	}
}

// --- VLPostBT ---

func TestVLPostBT_Eclipse(t *testing.T) {
	// During eclipse (t < 7), VL should be 0
	tcrit := FindTcrit(7.0, 0.00025, 0.8542, 336, 2.0)
	for _, tVal := range []float64{0, 3, 5, 6.999} {
		vl := VLPostBT(tVal, 7.0, 0.00025, 0.8542, 336, 0.7, 0.6, 1.0, tcrit, 2.0)
		if vl != 0 {
			t.Errorf("VLPostBT(%v) during eclipse = %v, want 0", tVal, vl)
		}
	}
}

func TestVLPostBT_ExponentialPhase(t *testing.T) {
	// Between eclipse and tcrit, VL = C0 * 2^((t-eclipse)/dt)
	tcrit := FindTcrit(7.0, 0.00025, 0.8542, 336, 2.0)
	tVal := 10.0
	got := VLPostBT(tVal, 7.0, 0.00025, 0.8542, 336, 0.7, 0.6, 1.0, tcrit, 2.0)
	expected := 0.00025 * math.Pow(2, (10-7)/0.8542)
	if math.Abs(got-expected)/expected > 1e-10 {
		t.Errorf("VLPostBT(10) exponential = %v, want %v", got, expected)
	}
	// Cross-validate with Python: VL(10) = 0.0028521662
	if math.Abs(got-0.0028521662) > 1e-6 {
		t.Errorf("VLPostBT(10) = %v, want ~0.0028521662", got)
	}
}

func TestVLPostBT_SetPointPhase(t *testing.T) {
	// After tcrit, VL = (set_point/copiesPerVirion) * SinVaried(t-tcrit, a, b, offset)
	tcrit := FindTcrit(7.0, 0.00025, 0.8542, 336, 2.0)
	tVal := 50.0
	got := VLPostBT(tVal, 7.0, 0.00025, 0.8542, 336, 0.7, 0.6, 1.0, tcrit, 2.0)
	expected := (336.0 / 2.0) * SinVaried(50-tcrit, 0.7, 0.6, 1.0)
	if math.Abs(got-expected) > 1e-6 {
		t.Errorf("VLPostBT(50) set-point = %v, want %v", got, expected)
	}
}

// --- ProbNondetectionSerology ---

func TestProbNondetectionSerology_BeforeMin(t *testing.T) {
	// Before ser_min, should return 1.0
	got := ProbNondetectionSerology(5, 10, 500, 9.1, 5.2)
	if got != 1.0 {
		t.Errorf("ProbNondetectionSerology(5) = %v, want 1.0", got)
	}
}

func TestProbNondetectionSerology_AfterMax(t *testing.T) {
	// After ser_max, should return 0.0
	got := ProbNondetectionSerology(501, 10, 500, 9.1, 5.2)
	if got != 0.0 {
		t.Errorf("ProbNondetectionSerology(501) = %v, want 0.0", got)
	}
}

func TestProbNondetectionSerology_InWindow(t *testing.T) {
	// At t=15, Python gives 0.9565472068
	got := ProbNondetectionSerology(15, 10, 500, 9.1, 5.2)
	if math.Abs(got-0.9565472068) > 1e-8 {
		t.Errorf("ProbNondetectionSerology(15) = %v, want ~0.9565472068", got)
	}
}

func TestProbNondetectionSerology_Monotonic(t *testing.T) {
	// Should decrease monotonically between min and max
	prev := 1.0
	for tVal := 11.0; tVal < 500; tVal += 1.0 {
		got := ProbNondetectionSerology(tVal, 10, 500, 9.1, 5.2)
		if got > prev {
			t.Errorf("ProbNondetectionSerology(%v) = %v > prev %v (not monotonic)", tVal, got, prev)
		}
		prev = got
	}
}

// --- ProbInfectiousPrep ---

func TestProbInfectiousPrep_DuringEclipse(t *testing.T) {
	// During eclipse, VL=0, so infectivity should be 0
	p := ProbInfectiousPrep(3.0, defaultPrepParams)
	if p != 0 {
		t.Errorf("ProbInfectiousPrep(3) during eclipse = %v, want 0", p)
	}
}

func TestProbInfectiousPrep_DuringGrowth(t *testing.T) {
	// During exponential growth, infectivity should be positive and increasing
	p10 := ProbInfectiousPrep(10.0, defaultPrepParams)
	p15 := ProbInfectiousPrep(15.0, defaultPrepParams)
	if p10 <= 0 {
		t.Errorf("ProbInfectiousPrep(10) = %v, want > 0", p10)
	}
	if p15 <= p10 {
		t.Errorf("ProbInfectiousPrep(15) = %v should be > ProbInfectiousPrep(10) = %v", p15, p10)
	}
}

// --- ProbNondetectionPrep ---

func TestProbNondetectionPrep_DuringEclipse(t *testing.T) {
	// During eclipse, Cc=0, so non-detection should be 1.0
	p := ProbNondetectionPrep(3.0, defaultPrepParams)
	if p != 1.0 {
		t.Errorf("ProbNondetectionPrep(3) during eclipse = %v, want 1.0", p)
	}
}

func TestProbNondetectionPrep_Decreases(t *testing.T) {
	// Non-detection should decrease as VL increases
	p10 := ProbNondetectionPrep(10.0, defaultPrepParams)
	p15 := ProbNondetectionPrep(15.0, defaultPrepParams)
	if p10 <= p15 {
		t.Errorf("ProbNondetectionPrep(10) = %v should be > ProbNondetectionPrep(15) = %v", p10, p15)
	}
}

// --- ProbInfectiousNondetectionPrep ---

func TestProbInfectiousNondetectionPrep_DuringEclipse(t *testing.T) {
	// During eclipse, integrand should be 0 (infectivity is 0)
	p := ProbInfectiousNondetectionPrep(3.0, defaultPrepParams)
	if p != 0 {
		t.Errorf("ProbInfectiousNondetectionPrep(3) = %v, want 0", p)
	}
}

func TestProbInfectiousNondetectionPrep_CrossValidate(t *testing.T) {
	// At t=15, Python gives 4.153462802811983e-02
	p := ProbInfectiousNondetectionPrep(15.0, defaultPrepParams)
	expected := 4.153462802811983e-02
	relErr := math.Abs(p-expected) / expected
	// Allow 1% tolerance due to analytic vs grid tcrit
	if relErr > 0.01 {
		t.Errorf("ProbInfectiousNondetectionPrep(15) = %e, want ~%e (relErr=%v)", p, expected, relErr)
	}
}

// --- RiskDaysPrep ---

func TestRiskDaysPrep_GoldenValue(t *testing.T) {
	// defaultPrepParams uses a NARROW serology window (ser_alpha=9.1, ser_beta=5.2):
	// the integrand has compact support over ~[8.7, 22.5] days. Cross-validated
	// against Python Simpson integration (100k points): risk_days_prep ≈ 1.0086.
	//
	// This is precisely the regime where scipy's adaptive quad silently misses
	// the peak and returns ~0 (≈5e-18); Go's fixed Gauss-Legendre — and Python's
	// "gauss-legendre" default — integrate it correctly. The golden is therefore
	// validated against the Simpson/GL truth, NOT against Python quad.
	// See TestRiskDaysPrep_GoldenValue_Production for the wide-window case.
	rd := RiskDaysPrep(defaultPrepParams)
	expected := 1.0086
	relErr := math.Abs(rd-expected) / expected
	// Allow 1% tolerance (fixed quadrature vs Simpson)
	if relErr > 0.01 {
		t.Errorf("RiskDaysPrep = %v, want ~%v (relErr=%v)", rd, expected, relErr)
	}
}

func TestRiskDaysPrep_GoldenValue_Production(t *testing.T) {
	// Production serology defaults (wide active window), set_point 336 copies/mL.
	// After the set-point units fix (a clinical copies/mL value is converted to the
	// model's virions/mL, ÷copiesPerVirion), Go GL and Python GL agree at 4.289877
	// — was 3.091868 with the pre-fix 2×-high plateau. Python↔Go parity verified to
	// ~1e-6 (see the fix-prep-setpoint-units branch validation).
	rd := RiskDaysPrep(productionPrepParams)
	expected := 4.289877
	relErr := math.Abs(rd-expected) / expected
	if relErr > 0.001 {
		t.Errorf("RiskDaysPrep(production) = %v, want ~%v (relErr=%v)", rd, expected, relErr)
	}
}

func TestRiskDaysPrep_Positive(t *testing.T) {
	rd := RiskDaysPrep(defaultPrepParams)
	if rd <= 0 {
		t.Errorf("RiskDaysPrep = %v, want > 0", rd)
	}
}

func TestRiskDaysPrep_HigherK(t *testing.T) {
	// Higher k → higher risk days (more infectious per virion)
	params1 := defaultPrepParams
	params1.K = 0.0001
	params2 := defaultPrepParams
	params2.K = 0.01

	rd1 := RiskDaysPrep(params1)
	rd2 := RiskDaysPrep(params2)
	if rd2 <= rd1 {
		t.Errorf("RiskDaysPrep(k=0.01)=%v should be > RiskDaysPrep(k=0.0001)=%v", rd2, rd1)
	}
}

func TestRiskDaysPrep_NoEclipse(t *testing.T) {
	// With eclipse=0, VL starts immediately → higher risk
	params1 := defaultPrepParams
	params1.Eclipse = 0
	params2 := defaultPrepParams
	params2.Eclipse = 14

	rd1 := RiskDaysPrep(params1)
	rd2 := RiskDaysPrep(params2)
	if rd1 <= rd2 {
		t.Errorf("RiskDaysPrep(eclipse=0)=%v should be > RiskDaysPrep(eclipse=14)=%v", rd1, rd2)
	}
}

// --- Bootstrap: RiskDaysBSPrep ---

func prepInput(nbs int, seed int64) RiskDaysInput {
	alpha := 2.0
	beta := 0.002019
	return RiskDaysInput{
		K:                   0.000673,
		DoublingTime:        0.8542,
		DoublingTimeNormSD:  0.2813,
		LOD50:               2.73,
		LOD50SD:             0.53,
		LOD95LOD50Ratio:     3.5,
		VolumeTransfused:    200,
		VolumeTransfusedMin: 100,
		VolumeTransfusedMax: 340,
		PoolSize:            16,
		Retests:             1,
		NBS:                 nbs,
		Seed:                seed,
		Threads:             2,
		PointEstimate:       "median",
		PrepMode:            true,
		SetPoint:            336,
		Eclipse:             7,
		A:                   0.7,
		B:                   0.6,
		Offset:              1.0,
		SerMin:              10,
		SerMax:              500,
		SerAlpha:            9.1,
		SerBeta:             5.2,
		KInvGammaAlpha:      &alpha,
		KInvGammaBeta:       &beta,
		ReturnParams:        true,
	}
}

func TestRiskDaysBSPrep_Sanity(t *testing.T) {
	input := prepInput(100, 42)
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatalf("RiskDaysBS (PrEP) error: %v", err)
	}
	if len(out.Simulations) != 100 {
		t.Errorf("got %d simulations, want 100", len(out.Simulations))
	}
	if out.PointEstimate <= 0 {
		t.Errorf("PE = %v, want > 0", out.PointEstimate)
	}
	if out.CredibleInterval[0] > out.PointEstimate || out.PointEstimate > out.CredibleInterval[1] {
		t.Errorf("PE %v not within CrI [%v, %v]", out.PointEstimate, out.CredibleInterval[0], out.CredibleInterval[1])
	}
}

func TestRiskDaysBSPrep_Reproducible(t *testing.T) {
	input := prepInput(50, 999)
	out1, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatal(err)
	}
	out2, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatal(err)
	}
	if out1.PointEstimate != out2.PointEstimate {
		t.Errorf("same seed → different PE: %v vs %v", out1.PointEstimate, out2.PointEstimate)
	}
	for i := 0; i < len(out1.Simulations); i++ {
		if out1.Simulations[i] != out2.Simulations[i] {
			t.Errorf("sim[%d] differs: %v vs %v", i, out1.Simulations[i], out2.Simulations[i])
			break
		}
	}
}

func TestRiskDaysBSPrep_DifferentSeeds(t *testing.T) {
	out1, _ := RiskDaysBS(prepInput(50, 1), nil)
	out2, _ := RiskDaysBS(prepInput(50, 2), nil)
	if out1.PointEstimate == out2.PointEstimate {
		t.Error("different seeds produced identical PE")
	}
}

func TestRiskDaysBSPrep_ReturnParams(t *testing.T) {
	input := prepInput(50, 42)
	input.ReturnParams = true
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatal(err)
	}
	if len(out.Ks) != 50 {
		t.Errorf("Ks len = %d, want 50", len(out.Ks))
	}
	if len(out.SetPoints) != 50 {
		t.Errorf("SetPoints len = %d, want 50", len(out.SetPoints))
	}
	if len(out.Eclipses) != 50 {
		t.Errorf("Eclipses len = %d, want 50", len(out.Eclipses))
	}
}

func TestRiskDaysBSPrep_PosteriorSample(t *testing.T) {
	input := prepInput(50, 42)
	input.KInvGammaAlpha = nil
	input.KInvGammaBeta = nil
	// Create a small posterior sample
	rng := NewRandomGenerator(99)
	input.KPosteriorSample = rng.GenerateGamma(2.0, 0.001, 200)
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatal(err)
	}
	if out.PointEstimate <= 0 {
		t.Errorf("PE = %v, want > 0", out.PointEstimate)
	}
}

func TestRiskDaysBSPrep_LnMix(t *testing.T) {
	input := prepInput(50, 42)
	input.KInvGammaAlpha = nil
	input.KInvGammaBeta = nil
	w := 0.90
	mu1 := -7.2403
	sigma1 := 0.3241
	mu2 := -3.7423
	sigma2 := 0.5258
	input.KLnMixW = &w
	input.KLnMixMu1 = &mu1
	input.KLnMixSigma1 = &sigma1
	input.KLnMixMu2 = &mu2
	input.KLnMixSigma2 = &sigma2
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatal(err)
	}
	if out.PointEstimate <= 0 {
		t.Errorf("PE = %v, want > 0", out.PointEstimate)
	}
}

func TestRiskDaysBSPrep_SetPointDist(t *testing.T) {
	// With uniform set-point distribution, sampled values should vary
	input := prepInput(100, 42)
	input.SetPointDistUniform = [2]float64{200, 500}
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatal(err)
	}
	// Check that set-points were sampled (not all the same)
	allSame := true
	for _, sp := range out.SetPoints[1:] {
		if sp != out.SetPoints[0] {
			allSame = false
			break
		}
	}
	if allSame {
		t.Error("SetPoints are all identical — uniform sampling not working")
	}
}

func TestRiskDaysBSPrep_Progress(t *testing.T) {
	input := prepInput(50, 42)
	var lastCompleted int
	progress := func(completed, total int) {
		if completed < lastCompleted {
			t.Errorf("progress went backwards: %d → %d", lastCompleted, completed)
		}
		lastCompleted = completed
	}
	_, err := RiskDaysBS(input, progress)
	if err != nil {
		t.Fatal(err)
	}
	if lastCompleted != 50 {
		t.Errorf("final progress = %d, want 50", lastCompleted)
	}
}

func TestRiskDaysBSPrep_PrimaryParamsPE(t *testing.T) {
	input := prepInput(50, 42)
	input.PointEstimate = "primary parameters"
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatal(err)
	}
	// The PE from primary params should match a single-call RiskDaysPrep
	pe := RiskDaysPrep(defaultPrepParams)
	// Allow some tolerance since bootstrap uses different LOD50 etc. for PE
	// but the primary params PE should match the single-call exactly
	if math.Abs(out.PointEstimate-pe)/pe > 0.01 {
		t.Errorf("PE (primary params) = %v, single-call = %v", out.PointEstimate, pe)
	}
}
