package riskdays

import "testing"

// ModeRounded documents itself as matching scipy.stats.mode, which returns the
// SMALLEST of the tied values. Go randomises map iteration order, so the previous
// implementation (a bare `count > maxCount`) returned an arbitrary tied value —
// different from run to run. Repeat the call many times: a nondeterministic
// implementation fails this with overwhelming probability.
func TestModeRounded_TieBreaksToSmallestDeterministically(t *testing.T) {
	// 1, 2 and 3 each appear twice — a three-way tie. scipy returns 1.
	data := []float64{3, 3, 1, 1, 2, 2}
	for i := 0; i < 200; i++ {
		if got := ModeRounded(data, 6); got != 1 {
			t.Fatalf("iteration %d: ModeRounded = %v, want 1 (smallest tied value)", i, got)
		}
	}
}

// A clear (untied) mode must still win, regardless of map ordering.
func TestModeRounded_ClearModeWins(t *testing.T) {
	data := []float64{1, 2, 2, 2, 3, 3}
	if got := ModeRounded(data, 6); got != 2 {
		t.Fatalf("ModeRounded = %v, want 2", got)
	}
}

// Guards the tie-break against the float64 zero value of `mode`: with all-negative
// data the initial mode==0 must never win a comparison.
func TestModeRounded_NegativeTieBreaksToSmallest(t *testing.T) {
	data := []float64{-2, -2, -5, -5}
	for i := 0; i < 100; i++ {
		if got := ModeRounded(data, 6); got != -5 {
			t.Fatalf("iteration %d: ModeRounded = %v, want -5 (smallest tied value)", i, got)
		}
	}
}

func TestModeRounded_EmptyReturnsZero(t *testing.T) {
	if got := ModeRounded(nil, 6); got != 0 {
		t.Fatalf("ModeRounded(nil) = %v, want 0", got)
	}
}
