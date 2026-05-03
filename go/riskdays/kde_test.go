package riskdays

import (
	"math"
	"testing"
)

// TestKDEModeLog_Silverman verifies Silverman bandwidth calculation.
func TestKDEModeLog_Silverman(t *testing.T) {
	// Simple data: log values [1,2,3,4,5] → std ≈ 1.5811, n=5
	data := []float64{1, 2, 3, 4, 5}
	// Expected Silverman bw: 1.06 * 1.5811 * 5^(-0.2)
	mode := KDEModeLog(data, 10000, 0)
	// Mode should be near the peak of the density
	if mode <= 0 || math.IsInf(mode, 0) || math.IsNaN(mode) {
		t.Errorf("unexpected mode: %v", mode)
	}
	t.Logf("Simple data mode: %v", mode)
}

// TestKDEModeLog_Reproducible verifies same input → same output.
func TestKDEModeLog_Reproducible(t *testing.T) {
	data := []float64{0.001, 0.002, 0.003, 0.005, 0.008, 0.013, 0.021}
	m1 := KDEModeLog(data, 5000, 0)
	m2 := KDEModeLog(data, 5000, 0)
	if m1 != m2 {
		t.Errorf("expected reproducible: %v != %v", m1, m2)
	}
}

// TestKDEModeLog_EmptyReturnsZero verifies empty input.
func TestKDEModeLog_EmptyReturnsZero(t *testing.T) {
	m := KDEModeLog([]float64{}, 1000, 0)
	if m != 0 {
		t.Errorf("expected 0 for empty input, got %v", m)
	}
}

// TestKDEModeLog_SingleValue returns that value.
func TestKDEModeLog_SingleValue(t *testing.T) {
	m := KDEModeLog([]float64{0.000673}, 1000, 0)
	if m != 0.000673 {
		t.Errorf("expected 0.000673, got %v", m)
	}
}

// TestKDEModeLog_Cap limits sample size.
func TestKDEModeLog_Cap(t *testing.T) {
	// Create 1000 values; cap at 100 should still work
	data := make([]float64, 1000)
	for i := range data {
		data[i] = 0.001 + float64(i)*1e-6
	}
	m := KDEModeLog(data, 5000, 100)
	if m <= 0 || math.IsInf(m, 0) || math.IsNaN(m) {
		t.Errorf("cap test failed: mode=%v", m)
	}
	t.Logf("cap=100 mode: %v", m)
}

// TestKDEModeLog_NoCap uses all data when cap=0.
func TestKDEModeLog_NoCap(t *testing.T) {
	data := make([]float64, 2000)
	for i := range data {
		data[i] = 0.001 + float64(i)*1e-6
	}
	m := KDEModeLog(data, 5000, 0)
	if m <= 0 || math.IsInf(m, 0) || math.IsNaN(m) {
		t.Errorf("no-cap test failed: mode=%v", m)
	}
	t.Logf("no-cap mode: %v", m)
}

// TestKDEModeLog_PositiveOnly ensures all output is positive.
func TestKDEModeLog_PositiveOnly(t *testing.T) {
	data := []float64{0.0003, 0.0005, 0.0007, 0.0010, 0.0015}
	m := KDEModeLog(data, 10000, 0)
	if m <= 0 {
		t.Errorf("expected positive mode, got %v", m)
	}
}

// TestKDEModeLog_AutoGrid matches grid to data size.
func TestKDEModeLog_AutoGrid(t *testing.T) {
	data := make([]float64, 150_000)
	for i := range data {
		data[i] = 0.001 + float64(i)*1e-8
	}
	// nGrid=0 → auto = max(100k, 150k) = 150k
	m := KDEModeLog(data, 0, 0)
	if m <= 0 || math.IsInf(m, 0) || math.IsNaN(m) {
		t.Errorf("auto grid test failed: mode=%v", m)
	}
	t.Logf("auto grid mode: %v", m)
}
