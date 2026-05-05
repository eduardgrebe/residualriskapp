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

// approxEqual returns true if a and b agree within the given relative tolerance,
// or within absTol when both values are near zero.
func approxEqual(a, b, relTol, absTol float64) bool {
	if math.Abs(a-b) <= absTol {
		return true
	}
	ref := math.Max(math.Abs(a), math.Abs(b))
	return math.Abs(a-b)/ref <= relTol
}

// Default parameter values matching app.py UI defaults.
const (
	defaultC0              = 0.00025
	defaultDoublingTime    = 20.5 / 24.0 // hours → days
	defaultVolumeTransfused = 20.0
	defaultK               = 0.013
	defaultLOD50           = 2.73
	defaultLOD95           = 12.33
	defaultZ               = 1.6449
	defaultCopiesPerVirion = 2
	defaultPoolSize        = 16
	defaultRetests         = 1
)

var defaultLOD95LOD50Ratio = defaultLOD95 / defaultLOD50

func defaultInnerParams() RiskDaysInnerParams {
	return RiskDaysInnerParams{
		CopiesPerVirion:  defaultCopiesPerVirion,
		C0:               defaultC0,
		DoublingTime:     defaultDoublingTime,
		VolumeTransfused: defaultVolumeTransfused,
		K:                defaultK,
		PoolSize:         defaultPoolSize,
		LOD50:            defaultLOD50,
		LOD95LOD50Ratio:  defaultLOD95LOD50Ratio,
		Retests:          defaultRetests,
		Z:                defaultZ,
		LimitMin:         -100,
		LimitMax:         500,
	}
}

// ---------------------------------------------------------------------------
// Concentration
// ---------------------------------------------------------------------------

func TestConcentration_AtT0_EqualsC0(t *testing.T) {
	got := Concentration(defaultC0, defaultDoublingTime, 0)
	if !approxEqual(got, defaultC0, 1e-9, 0) {
		t.Errorf("Concentration at t=0: got %v, want %v", got, defaultC0)
	}
}

func TestConcentration_AtOneDoublingTime_Doubles(t *testing.T) {
	got := Concentration(defaultC0, defaultDoublingTime, defaultDoublingTime)
	want := 2 * defaultC0
	if !approxEqual(got, want, 1e-9, 0) {
		t.Errorf("Concentration at t=doubling_time: got %v, want %v", got, want)
	}
}

func TestConcentration_AtTwoDoublingTimes_Quadruples(t *testing.T) {
	got := Concentration(defaultC0, defaultDoublingTime, 2*defaultDoublingTime)
	want := 4 * defaultC0
	if !approxEqual(got, want, 1e-9, 0) {
		t.Errorf("Concentration at t=2*doubling_time: got %v, want %v", got, want)
	}
}

func TestConcentration_IncreasesWithT(t *testing.T) {
	times := []float64{0, 1, 5, 10, 20}
	prev := Concentration(defaultC0, defaultDoublingTime, times[0])
	for _, ts := range times[1:] {
		curr := Concentration(defaultC0, defaultDoublingTime, ts)
		if curr <= prev {
			t.Errorf("Concentration not increasing: C(%v)=%v <= C(prev)=%v", ts, curr, prev)
		}
		prev = curr
	}
}

// ---------------------------------------------------------------------------
// ProbInfectiousCopies
// ---------------------------------------------------------------------------

func TestProbInfectiousCopies_ZeroCopies_GivesZero(t *testing.T) {
	got := ProbInfectiousCopies(0, 0.013)
	if !approxEqual(got, 0.0, 0, 1e-9) {
		t.Errorf("ProbInfectiousCopies(0, k): got %v, want 0", got)
	}
}

func TestProbInfectiousCopies_StandardSingleHitFormula(t *testing.T) {
	// 1 - exp(-0.01 * 100) = 1 - exp(-1) ≈ 0.6321
	want := 1.0 - math.Exp(-0.01*100)
	got := ProbInfectiousCopies(100, 0.01)
	if !approxEqual(got, want, 1e-9, 0) {
		t.Errorf("ProbInfectiousCopies(100, 0.01): got %v, want %v", got, want)
	}
}

func TestProbInfectiousCopies_LargeCopies_ApproachesOne(t *testing.T) {
	got := ProbInfectiousCopies(1_000_000, 0.013)
	if !approxEqual(got, 1.0, 0, 1e-6) {
		t.Errorf("ProbInfectiousCopies(1e6, k): got %v, want ~1.0", got)
	}
}

func TestProbInfectiousCopies_ResultBounded(t *testing.T) {
	for _, n := range []float64{0, 1, 10, 100, 1000} {
		p := ProbInfectiousCopies(n, 0.013)
		if p < 0 || p > 1 {
			t.Errorf("ProbInfectiousCopies(%v, k) = %v, outside [0, 1]", n, p)
		}
	}
}

func TestProbInfectiousCopies_IncreasesWithCopies(t *testing.T) {
	copies := []float64{0, 1, 10, 100, 1000}
	prev := ProbInfectiousCopies(copies[0], 0.013)
	for _, n := range copies[1:] {
		curr := ProbInfectiousCopies(n, 0.013)
		if curr <= prev {
			t.Errorf("ProbInfectiousCopies not increasing: P(%v)=%v <= P(prev)=%v", n, curr, prev)
		}
		prev = curr
	}
}

// ---------------------------------------------------------------------------
// ProbPosInit
// ---------------------------------------------------------------------------

func TestProbPosInit_AtPoolLOD50_GivesHalf(t *testing.T) {
	// When C = poolSize * lod50 the pooled concentration equals lod50,
	// so the initial detection probability should be 0.5.
	C := float64(defaultPoolSize) * defaultLOD50
	got, err := ProbPosInit(C, defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !approxEqual(got, 0.5, 1e-6, 0) {
		t.Errorf("ProbPosInit at pool*lod50: got %v, want 0.5", got)
	}
}

func TestProbPosInit_IncreasesWithConcentration(t *testing.T) {
	concentrations := []float64{defaultLOD50, 10 * defaultLOD50, 100 * defaultLOD50, 1000 * defaultLOD50}
	prev, err := ProbPosInit(concentrations[0], defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	for _, C := range concentrations[1:] {
		curr, err := ProbPosInit(C, defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, defaultZ)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if curr <= prev {
			t.Errorf("ProbPosInit not increasing with concentration: P(%v)=%v <= P(prev)=%v", C, curr, prev)
		}
		prev = curr
	}
}

func TestProbPosInit_InvalidPoolSize_ReturnsError(t *testing.T) {
	_, err := ProbPosInit(10.0, defaultDoublingTime, 0, defaultLOD50, defaultLOD95LOD50Ratio, defaultZ)
	if err == nil {
		t.Error("expected error for pool_size=0, got nil")
	}
}

// ---------------------------------------------------------------------------
// ProbNegRetest
// ---------------------------------------------------------------------------

func TestProbNegRetest_ZeroRetests_ReturnsZero(t *testing.T) {
	got, err := ProbNegRetest(100.0, defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, 0, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != 0 {
		t.Errorf("ProbNegRetest with retests=0: got %v, want 0", got)
	}
}

func TestProbNegRetest_HighConcentration_NearZero(t *testing.T) {
	got, err := ProbNegRetest(1e8, defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, 1, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !approxEqual(got, 0.0, 0, 1e-6) {
		t.Errorf("ProbNegRetest at very high concentration: got %v, want ~0", got)
	}
}

func TestProbNegRetest_DecreasesWithConcentration(t *testing.T) {
	concentrations := []float64{defaultLOD50, 10 * defaultLOD50, 1000 * defaultLOD50}
	prev, err := ProbNegRetest(concentrations[0], defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, 1, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	for _, C := range concentrations[1:] {
		curr, err := ProbNegRetest(C, defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, 1, defaultZ)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if curr >= prev {
			t.Errorf("ProbNegRetest not decreasing with concentration: P(%v)=%v >= P(prev)=%v", C, curr, prev)
		}
		prev = curr
	}
}

func TestProbNegRetest_NegativeRetests_ReturnsError(t *testing.T) {
	_, err := ProbNegRetest(10.0, defaultDoublingTime, defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, -1, defaultZ)
	if err == nil {
		t.Error("expected error for retests=-1, got nil")
	}
}

// ---------------------------------------------------------------------------
// ProbNondetection
// ---------------------------------------------------------------------------

func TestProbNondetection_VeryEarlyTime_IsOne(t *testing.T) {
	got, err := ProbNondetection(-50, defaultC0, defaultDoublingTime, defaultCopiesPerVirion,
		defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, defaultRetests, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !approxEqual(got, 1.0, 0, 1e-6) {
		t.Errorf("ProbNondetection at t=-50: got %v, want ~1.0", got)
	}
}

func TestProbNondetection_VeryLateTime_IsZero(t *testing.T) {
	got, err := ProbNondetection(100, defaultC0, defaultDoublingTime, defaultCopiesPerVirion,
		defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, defaultRetests, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !approxEqual(got, 0.0, 0, 1e-6) {
		t.Errorf("ProbNondetection at t=100: got %v, want ~0", got)
	}
}

func TestProbNondetection_DecreasesOverTime(t *testing.T) {
	times := []float64{-20, -10, 0, 10, 20}
	prev, err := ProbNondetection(times[0], defaultC0, defaultDoublingTime, defaultCopiesPerVirion,
		defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, defaultRetests, defaultZ)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	for _, ts := range times[1:] {
		curr, err := ProbNondetection(ts, defaultC0, defaultDoublingTime, defaultCopiesPerVirion,
			defaultPoolSize, defaultLOD50, defaultLOD95LOD50Ratio, defaultRetests, defaultZ)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if curr > prev {
			t.Errorf("ProbNondetection not decreasing: P(%v)=%v > P(prev)=%v", ts, curr, prev)
		}
		prev = curr
	}
}

// ---------------------------------------------------------------------------
// RiskDays (deterministic integral)
// ---------------------------------------------------------------------------

func TestRiskDays_GoldenValue_DefaultParams(t *testing.T) {
	got, err := RiskDays(defaultInnerParams())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	// Reference value matches the Python _risk_days() output to within 0.5%.
	want := 3.7207
	if !approxEqual(got, want, 0.005, 0) {
		t.Errorf("RiskDays default params: got %v, want ~%v", got, want)
	}
}

func TestRiskDays_IDNAT_GoldenValue(t *testing.T) {
	params := defaultInnerParams()
	params.PoolSize = 1
	params.Retests = 0
	got, err := RiskDays(params)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	want := 0.9001
	if !approxEqual(got, want, 0.005, 0) {
		t.Errorf("RiskDays ID-NAT: got %v, want ~%v", got, want)
	}
}

func TestRiskDays_HigherK_MoreRiskDays(t *testing.T) {
	pLow := defaultInnerParams()
	pLow.K = 0.005
	pHigh := defaultInnerParams()
	pHigh.K = 0.05

	rdLow, err := RiskDays(pLow)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	rdHigh, err := RiskDays(pHigh)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if rdHigh <= rdLow {
		t.Errorf("higher k should give more risk days: k=0.05 → %v, k=0.005 → %v", rdHigh, rdLow)
	}
}

func TestRiskDays_LargerVolume_MoreRiskDays(t *testing.T) {
	pSmall := defaultInnerParams()
	pSmall.VolumeTransfused = 5
	pLarge := defaultInnerParams()
	pLarge.VolumeTransfused = 100

	rdSmall, err := RiskDays(pSmall)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	rdLarge, err := RiskDays(pLarge)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if rdLarge <= rdSmall {
		t.Errorf("larger volume should give more risk days: v=100 → %v, v=5 → %v", rdLarge, rdSmall)
	}
}

func TestRiskDays_IDNATLessThanMinipool(t *testing.T) {
	pMinipool := defaultInnerParams()
	pIDNAT := defaultInnerParams()
	pIDNAT.PoolSize = 1
	pIDNAT.Retests = 0

	rdMinipool, err := RiskDays(pMinipool)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	rdIDNAT, err := RiskDays(pIDNAT)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if rdIDNAT >= rdMinipool {
		t.Errorf("ID-NAT should give fewer risk days than minipool: ID-NAT=%v, minipool=%v", rdIDNAT, rdMinipool)
	}
}

func TestRiskDays_ResultIsPositive(t *testing.T) {
	got, err := RiskDays(defaultInnerParams())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got <= 0 {
		t.Errorf("RiskDays should be positive, got %v", got)
	}
}

// ---------------------------------------------------------------------------
// RiskDaysInput: SetDefaults and Validate
// ---------------------------------------------------------------------------

func TestSetDefaults_FillsMissingValues(t *testing.T) {
	input := RiskDaysInput{}
	input.SetDefaults()

	if input.C0 != 0.00025 {
		t.Errorf("C0 default: got %v, want 0.00025", input.C0)
	}
	if input.CopiesPerVirion != 2 {
		t.Errorf("CopiesPerVirion default: got %v, want 2", input.CopiesPerVirion)
	}
	if input.Alpha != 0.05 {
		t.Errorf("Alpha default: got %v, want 0.05", input.Alpha)
	}
	if input.Z != 1.6449 {
		t.Errorf("Z default: got %v, want 1.6449", input.Z)
	}
	if input.NBS != 10000 {
		t.Errorf("NBS default: got %v, want 10000", input.NBS)
	}
	if input.PointEstimate != "primary parameters" {
		t.Errorf("PointEstimate default: got %q, want %q", input.PointEstimate, "primary parameters")
	}
}

func TestValidate_InvalidNBS_ReturnsError(t *testing.T) {
	kSample := []float64{0.013}
	input := RiskDaysInput{NBS: 0, PoolSize: 16, Retests: 1, KPosteriorSample: kSample}
	if err := input.Validate(); err == nil {
		t.Error("expected error for n_bs=0, got nil")
	}
}

func TestValidate_InvalidPoolSize_ReturnsError(t *testing.T) {
	kSample := []float64{0.013}
	input := RiskDaysInput{NBS: 100, PoolSize: 0, Retests: 1, KPosteriorSample: kSample}
	if err := input.Validate(); err == nil {
		t.Error("expected error for pool_size=0, got nil")
	}
}

func TestValidate_MissingKDistribution_ReturnsError(t *testing.T) {
	input := RiskDaysInput{NBS: 100, PoolSize: 16, Retests: 1}
	if err := input.Validate(); err == nil {
		t.Error("expected error when no k distribution provided, got nil")
	}
}

// ---------------------------------------------------------------------------
// RiskDaysBS (bootstrap simulation)
// ---------------------------------------------------------------------------

func defaultBSInput() RiskDaysInput {
	kSample := make([]float64, 1000)
	for i := range kSample {
		kSample[i] = 0.006 + float64(i%3)*0.007 // simple synthetic distribution
	}
	input := RiskDaysInput{
		K:                   defaultK,
		DoublingTime:        defaultDoublingTime,
		DoublingTimeNormSD:  1.33 / 24.0,
		LOD50:               defaultLOD50,
		LOD50SD:             0.193,
		LOD95LOD50Ratio:     defaultLOD95LOD50Ratio,
		VolumeTransfused:    defaultVolumeTransfused,
		VolumeTransfusedMin: 15,
		VolumeTransfusedMax: 30,
		PoolSize:            defaultPoolSize,
		Retests:             defaultRetests,
		KPosteriorSample:    kSample,
		NBS:                 200,
		Seed:                42,
		Threads:             2,
		PointEstimate:       "primary parameters",
	}
	input.SetDefaults()
	return input
}

func TestRiskDaysBS_ReturnsResults(t *testing.T) {
	out, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	if out == nil {
		t.Fatal("RiskDaysBS returned nil output")
	}
}

func TestRiskDaysBS_SimulationCount(t *testing.T) {
	input := defaultBSInput()
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	if len(out.Simulations) != input.NBS {
		t.Errorf("simulation count: got %v, want %v", len(out.Simulations), input.NBS)
	}
}

func TestRiskDaysBS_AllSimulationsPositive(t *testing.T) {
	out, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	for i, v := range out.Simulations {
		if v <= 0 {
			t.Errorf("simulation[%d] = %v, want > 0", i, v)
		}
	}
}

func TestRiskDaysBS_PointEstimatePositive(t *testing.T) {
	out, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	if out.PointEstimate <= 0 {
		t.Errorf("point estimate: got %v, want > 0", out.PointEstimate)
	}
}

func TestRiskDaysBS_CrIIsOrdered(t *testing.T) {
	out, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	if out.CredibleInterval[0] >= out.CredibleInterval[1] {
		t.Errorf("CrI not ordered: [%v, %v]", out.CredibleInterval[0], out.CredibleInterval[1])
	}
}

func TestRiskDaysBS_RangeContainsCrI(t *testing.T) {
	out, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	if out.Range[0] > out.CredibleInterval[0] {
		t.Errorf("range min %v > CrI lower %v", out.Range[0], out.CredibleInterval[0])
	}
	if out.CredibleInterval[1] > out.Range[1] {
		t.Errorf("CrI upper %v > range max %v", out.CredibleInterval[1], out.Range[1])
	}
}

func TestRiskDaysBS_PointEstimate_MatchesRiskDays(t *testing.T) {
	// With point_estimate="primary parameters" the pe must equal RiskDays at primary params.
	out, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	expected, err := RiskDays(defaultInnerParams())
	if err != nil {
		t.Fatalf("RiskDays failed: %v", err)
	}
	if !approxEqual(out.PointEstimate, expected, 1e-6, 0) {
		t.Errorf("point estimate: got %v, want %v", out.PointEstimate, expected)
	}
}

func TestRiskDaysBS_Reproducible(t *testing.T) {
	r1, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("first run failed: %v", err)
	}
	r2, err := RiskDaysBS(defaultBSInput(), nil)
	if err != nil {
		t.Fatalf("second run failed: %v", err)
	}
	if r1.PointEstimate != r2.PointEstimate {
		t.Errorf("point estimates differ across runs: %v vs %v", r1.PointEstimate, r2.PointEstimate)
	}
	for i := range r1.Simulations {
		if r1.Simulations[i] != r2.Simulations[i] {
			t.Errorf("simulation[%d] differs: %v vs %v", i, r1.Simulations[i], r2.Simulations[i])
			break
		}
	}
}

func TestRiskDaysBS_DifferentSeeds_DifferentResults(t *testing.T) {
	i1 := defaultBSInput()
	i1.Seed = 42
	i2 := defaultBSInput()
	i2.Seed = 99999

	r1, err := RiskDaysBS(i1, nil)
	if err != nil {
		t.Fatalf("first run failed: %v", err)
	}
	r2, err := RiskDaysBS(i2, nil)
	if err != nil {
		t.Fatalf("second run failed: %v", err)
	}
	// It would be astronomically unlikely for all values to match with different seeds.
	allSame := true
	for i := range r1.Simulations {
		if r1.Simulations[i] != r2.Simulations[i] {
			allSame = false
			break
		}
	}
	if allSame {
		t.Error("different seeds produced identical simulation results")
	}
}

func TestRiskDaysBS_InvalidNBS_ReturnsError(t *testing.T) {
	input := defaultBSInput()
	input.NBS = -1 // 0 is replaced by SetDefaults; use -1 to reach Validate
	_, err := RiskDaysBS(input, nil)
	if err == nil {
		t.Error("expected error for n_bs=-1, got nil")
	}
}

// ---------------------------------------------------------------------------
// GenerateInvGamma
// ---------------------------------------------------------------------------

func TestGenerateInvGamma_AllPositive(t *testing.T) {
	rng := NewRandomGenerator(42)
	samples := rng.GenerateInvGamma(2.0, 0.002019, 1000)
	for i, v := range samples {
		if v <= 0 {
			t.Errorf("sample[%d] = %v, want > 0", i, v)
		}
	}
}

func TestGenerateInvGamma_MeanApproximate(t *testing.T) {
	// InvGamma(alpha, beta): mean = beta / (alpha - 1) for alpha > 1.
	// With alpha=2.0, beta=0.002019: mean = 0.002019 / 1 = 0.002019.
	rng := NewRandomGenerator(42)
	n := 100_000
	alpha, beta := 2.0, 0.002019
	samples := rng.GenerateInvGamma(alpha, beta, n)
	sum := 0.0
	for _, v := range samples {
		sum += v
	}
	mean := sum / float64(n)
	expectedMean := beta / (alpha - 1)
	// Allow 3% relative tolerance for a large sample.
	if !approxEqual(mean, expectedMean, 0.03, 0) {
		t.Errorf("GenerateInvGamma mean: got %v, want ~%v", mean, expectedMean)
	}
}

func TestGenerateInvGamma_ModeApproximate(t *testing.T) {
	// InvGamma(alpha=2, beta=0.002019): mode = beta / (alpha + 1) = 0.002019 / 3 ≈ 0.000673.
	// Verify that the sample mode (via histogram peak) is in the right ballpark.
	rng := NewRandomGenerator(42)
	samples := rng.GenerateInvGamma(2.0, 0.002019, 100_000)
	// Rough check: most samples should be below the mean (right-skewed distribution).
	mean := 0.002019
	belowMean := 0
	for _, v := range samples {
		if v < mean {
			belowMean++
		}
	}
	fraction := float64(belowMean) / float64(len(samples))
	// For InvGamma(2, ...) most mass is below the mean; expect > 60%.
	if fraction < 0.60 {
		t.Errorf("expected >60%% of samples below mean, got %.1f%%", fraction*100)
	}
}

func TestGenerateInvGamma_Reproducible(t *testing.T) {
	r1 := NewRandomGenerator(42).GenerateInvGamma(2.0, 0.002019, 500)
	r2 := NewRandomGenerator(42).GenerateInvGamma(2.0, 0.002019, 500)
	for i := range r1 {
		if r1[i] != r2[i] {
			t.Errorf("sample[%d] differs across identical seeds: %v vs %v", i, r1[i], r2[i])
			break
		}
	}
}

// ---------------------------------------------------------------------------
// RiskDaysBS with InvGamma k distribution
// ---------------------------------------------------------------------------

func invGammaBSInput() RiskDaysInput {
	alpha := 2.0
	beta := 0.002019
	input := RiskDaysInput{
		K:                   0.000673, // InvGamma mode (beta / (alpha + 1))
		DoublingTime:        defaultDoublingTime,
		DoublingTimeNormSD:  1.33 / 24.0,
		LOD50:               defaultLOD50,
		LOD50SD:             0.193,
		LOD95LOD50Ratio:     defaultLOD95LOD50Ratio,
		VolumeTransfused:    defaultVolumeTransfused,
		VolumeTransfusedMin: 15,
		VolumeTransfusedMax: 30,
		PoolSize:            defaultPoolSize,
		Retests:             defaultRetests,
		KInvGammaAlpha:      &alpha,
		KInvGammaBeta:       &beta,
		NBS:                 200,
		Seed:                42,
		Threads:             2,
		PointEstimate:       "primary parameters",
	}
	input.SetDefaults()
	return input
}

func TestRiskDaysBS_InvGamma_Sanity(t *testing.T) {
	out, err := RiskDaysBS(invGammaBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS with InvGamma failed: %v", err)
	}
	if out == nil {
		t.Fatal("RiskDaysBS returned nil output")
	}
	if len(out.Simulations) != 200 {
		t.Errorf("simulation count: got %v, want 200", len(out.Simulations))
	}
	for i, v := range out.Simulations {
		if v <= 0 {
			t.Errorf("simulation[%d] = %v, want > 0", i, v)
		}
	}
	if out.PointEstimate <= 0 {
		t.Errorf("point estimate: got %v, want > 0", out.PointEstimate)
	}
	if out.CredibleInterval[0] >= out.CredibleInterval[1] {
		t.Errorf("CrI not ordered: [%v, %v]", out.CredibleInterval[0], out.CredibleInterval[1])
	}
}

func TestRiskDaysBS_InvGamma_Reproducible(t *testing.T) {
	r1, err := RiskDaysBS(invGammaBSInput(), nil)
	if err != nil {
		t.Fatalf("first run failed: %v", err)
	}
	r2, err := RiskDaysBS(invGammaBSInput(), nil)
	if err != nil {
		t.Fatalf("second run failed: %v", err)
	}
	if r1.PointEstimate != r2.PointEstimate {
		t.Errorf("point estimates differ: %v vs %v", r1.PointEstimate, r2.PointEstimate)
	}
	for i := range r1.Simulations {
		if r1.Simulations[i] != r2.Simulations[i] {
			t.Errorf("simulation[%d] differs: %v vs %v", i, r1.Simulations[i], r2.Simulations[i])
			break
		}
	}
}

func TestRiskDaysBS_InvGamma_PointEstimateMatchesRiskDays(t *testing.T) {
	// With point_estimate="primary parameters", pe must equal RiskDays(K=mode).
	input := invGammaBSInput()
	out, err := RiskDaysBS(input, nil)
	if err != nil {
		t.Fatalf("RiskDaysBS failed: %v", err)
	}
	params := defaultInnerParams()
	params.K = input.K // K = mode of InvGamma
	expected, err := RiskDays(params)
	if err != nil {
		t.Fatalf("RiskDays failed: %v", err)
	}
	if !approxEqual(out.PointEstimate, expected, 1e-6, 0) {
		t.Errorf("point estimate: got %v, want %v", out.PointEstimate, expected)
	}
}

// ---------------------------------------------------------------------------
// GenerateLogNormalMixture
// ---------------------------------------------------------------------------

// Default mixture parameters: 90% human + 10% animal (Recommendation B)
const (
	defaultLnMixW      = 0.90
	defaultLnMixMu1    = -7.2403
	defaultLnMixSigma1 = 0.3241
	defaultLnMixMu2    = -3.7423
	defaultLnMixSigma2 = 0.5258
)

func TestGenerateLogNormalMixture_AllPositive(t *testing.T) {
	rng := NewRandomGenerator(42)
	samples := rng.GenerateLogNormalMixture(
		defaultLnMixW, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, 1000)
	for i, v := range samples {
		if v <= 0 {
			t.Errorf("sample[%d] = %v, want > 0", i, v)
		}
	}
}

func TestGenerateLogNormalMixture_Reproducible(t *testing.T) {
	r1 := NewRandomGenerator(42).GenerateLogNormalMixture(
		defaultLnMixW, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, 500)
	r2 := NewRandomGenerator(42).GenerateLogNormalMixture(
		defaultLnMixW, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, 500)
	for i := range r1 {
		if r1[i] != r2[i] {
			t.Errorf("sample[%d] differs across identical seeds: %v vs %v", i, r1[i], r2[i])
			break
		}
	}
}

func TestGenerateLogNormalMixture_DifferentSeedsAreIndependent(t *testing.T) {
	r1 := NewRandomGenerator(42).GenerateLogNormalMixture(
		defaultLnMixW, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, 100)
	r2 := NewRandomGenerator(99).GenerateLogNormalMixture(
		defaultLnMixW, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, 100)
	allSame := true
	for i := range r1 {
		if r1[i] != r2[i] {
			allSame = false
			break
		}
	}
	if allSame {
		t.Error("samples from different seeds are identical (should differ)")
	}
}

func TestGenerateLogNormalMixture_ComponentIsolation_W1(t *testing.T) {
	// w=1 → all samples from component 1 (human posterior)
	// Component 1 median = exp(mu1) ≈ exp(-7.2403) ≈ 0.000715
	rng := NewRandomGenerator(42)
	samples := rng.GenerateLogNormalMixture(1.0, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, 10000)
	sorted := make([]float64, len(samples))
	copy(sorted, samples)
	// Compute median
	n := len(sorted)
	sum := 0.0
	for _, v := range sorted {
		sum += math.Log(v)
	}
	logMean := sum / float64(n)
	expectedLogMean := defaultLnMixMu1
	// Allow 5% relative tolerance
	if !approxEqual(logMean, expectedLogMean, 0.05, 0.001) {
		t.Errorf("w=1: log-mean of samples = %v, want ~%v (mu1)", logMean, expectedLogMean)
	}
}

func TestGenerateLogNormalMixture_ComponentIsolation_W0(t *testing.T) {
	// w=0 → all samples from component 2 (animal posterior)
	// Component 2 median = exp(mu2) ≈ exp(-3.7423) ≈ 0.0237
	rng := NewRandomGenerator(42)
	samples := rng.GenerateLogNormalMixture(0.0, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, 10000)
	sum := 0.0
	for _, v := range samples {
		sum += math.Log(v)
	}
	logMean := sum / float64(len(samples))
	expectedLogMean := defaultLnMixMu2
	if !approxEqual(logMean, expectedLogMean, 0.05, 0.001) {
		t.Errorf("w=0: log-mean of samples = %v, want ~%v (mu2)", logMean, expectedLogMean)
	}
}

func TestGenerateLogNormalMixture_MedianPlausible(t *testing.T) {
	// Expected mixture median ≈ 0.000750 (from companion analysis)
	rng := NewRandomGenerator(42)
	n := 100_000
	samples := rng.GenerateLogNormalMixture(
		defaultLnMixW, defaultLnMixMu1, defaultLnMixSigma1,
		defaultLnMixMu2, defaultLnMixSigma2, n)
	// Approximate median via counting
	count := 0
	threshold := 0.000750
	for _, v := range samples {
		if v < threshold {
			count++
		}
	}
	fraction := float64(count) / float64(n)
	// Should be close to 0.5 (within ±10%)
	if fraction < 0.40 || fraction > 0.60 {
		t.Errorf("fraction below expected median threshold %v: got %.3f, want ~0.50", threshold, fraction)
	}
}

// ---------------------------------------------------------------------------
// RiskDaysBS with lognormal mixture k distribution
// ---------------------------------------------------------------------------

func lnMixBSInput() RiskDaysInput {
	w := defaultLnMixW
	mu1 := defaultLnMixMu1
	s1 := defaultLnMixSigma1
	mu2 := defaultLnMixMu2
	s2 := defaultLnMixSigma2
	input := RiskDaysInput{
		K:                   0.000649, // approximate mixture mode
		DoublingTime:        defaultDoublingTime,
		DoublingTimeNormSD:  1.33 / 24.0,
		LOD50:               defaultLOD50,
		LOD50SD:             0.193,
		LOD95LOD50Ratio:     defaultLOD95LOD50Ratio,
		VolumeTransfused:    defaultVolumeTransfused,
		VolumeTransfusedMin: 15,
		VolumeTransfusedMax: 30,
		PoolSize:            defaultPoolSize,
		Retests:             defaultRetests,
		KLnMixW:             &w,
		KLnMixMu1:          &mu1,
		KLnMixSigma1:       &s1,
		KLnMixMu2:          &mu2,
		KLnMixSigma2:       &s2,
		NBS:                 200,
		Seed:                42,
		Threads:             2,
		PointEstimate:       "primary parameters",
	}
	input.SetDefaults()
	return input
}

func TestRiskDaysBS_LnMix_Sanity(t *testing.T) {
	out, err := RiskDaysBS(lnMixBSInput(), nil)
	if err != nil {
		t.Fatalf("RiskDaysBS with LnMix failed: %v", err)
	}
	if out == nil {
		t.Fatal("RiskDaysBS returned nil output")
	}
	if len(out.Simulations) != 200 {
		t.Errorf("simulation count: got %v, want 200", len(out.Simulations))
	}
	for i, v := range out.Simulations {
		if v <= 0 {
			t.Errorf("simulation[%d] = %v, want > 0", i, v)
		}
	}
	if out.PointEstimate <= 0 {
		t.Errorf("point estimate: got %v, want > 0", out.PointEstimate)
	}
	if out.CredibleInterval[0] >= out.CredibleInterval[1] {
		t.Errorf("CrI not ordered: [%v, %v]", out.CredibleInterval[0], out.CredibleInterval[1])
	}
}

func TestRiskDaysBS_LnMix_Reproducible(t *testing.T) {
	r1, err := RiskDaysBS(lnMixBSInput(), nil)
	if err != nil {
		t.Fatalf("first run failed: %v", err)
	}
	r2, err := RiskDaysBS(lnMixBSInput(), nil)
	if err != nil {
		t.Fatalf("second run failed: %v", err)
	}
	if r1.PointEstimate != r2.PointEstimate {
		t.Errorf("point estimates differ: %v vs %v", r1.PointEstimate, r2.PointEstimate)
	}
	for i := range r1.Simulations {
		if r1.Simulations[i] != r2.Simulations[i] {
			t.Errorf("simulation[%d] differs: %v vs %v", i, r1.Simulations[i], r2.Simulations[i])
			break
		}
	}
}
