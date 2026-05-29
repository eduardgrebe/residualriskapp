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
	"fmt"
	"sync"
)

// ProgressCallback is a function type for progress updates
type ProgressCallback func(completed, total int)

// RiskDaysBS performs bootstrap simulation to calculate risk days with uncertainty
// This is the main entry point and corresponds to Python risk_days_bs()
func RiskDaysBS(input RiskDaysInput, progressCallback ProgressCallback) (*RiskDaysOutput, error) {
	// Set defaults and validate
	input.SetDefaults()
	if err := input.Validate(); err != nil {
		return nil, err
	}

	// Initialize random generator
	rng := NewRandomGenerator(input.Seed)

	// Generate random samples for all parameters
	var ks []float64
	if input.KPosteriorSample != nil {
		ks = rng.BootstrapChoice(input.KPosteriorSample, input.NBS)
	} else if input.KGammaShape != nil && input.KGammaScale != nil {
		ks = rng.GenerateGamma(*input.KGammaShape, *input.KGammaScale, input.NBS)
	} else if input.KInvGammaAlpha != nil && input.KInvGammaBeta != nil {
		ks = rng.GenerateInvGamma(*input.KInvGammaAlpha, *input.KInvGammaBeta, input.NBS)
	} else if input.KLnMixW != nil && input.KLnMixMu1 != nil && input.KLnMixSigma1 != nil &&
		input.KLnMixMu2 != nil && input.KLnMixSigma2 != nil {
		ks = rng.GenerateLogNormalMixture(
			*input.KLnMixW, *input.KLnMixMu1, *input.KLnMixSigma1,
			*input.KLnMixMu2, *input.KLnMixSigma2, input.NBS,
		)
	} else {
		return nil, fmt.Errorf("no valid k distribution specified")
	}

	doublingTimes := rng.GenerateTruncatedNormal(input.DoublingTime, input.DoublingTimeNormSD, input.NBS)
	lod50s := rng.GenerateTruncatedNormal(input.LOD50, input.LOD50SD, input.NBS)
	volumesTransfused := rng.GenerateUniform(input.VolumeTransfusedMin, input.VolumeTransfusedMax, input.NBS)

	if input.PrepMode {
		return riskDaysBSPrep(input, rng, ks, doublingTimes, lod50s, volumesTransfused, progressCallback)
	}
	return riskDaysBSBaseline(input, ks, doublingTimes, lod50s, volumesTransfused, progressCallback)
}

// riskDaysBSBaseline runs the standard (non-PrEP) bootstrap simulation.
func riskDaysBSBaseline(input RiskDaysInput, ks, doublingTimes, lod50s, volumesTransfused []float64, progressCallback ProgressCallback) (*RiskDaysOutput, error) {

	// Prepare args list for parallel execution
	argsList := make([]RiskDaysInnerParams, input.NBS)
	for i := 0; i < input.NBS; i++ {
		argsList[i] = RiskDaysInnerParams{
			CopiesPerVirion:  input.CopiesPerVirion,
			C0:               input.C0,
			DoublingTime:     doublingTimes[i],
			VolumeTransfused: volumesTransfused[i],
			K:                ks[i],
			PoolSize:         input.PoolSize,
			LOD50:            lod50s[i],
			LOD95LOD50Ratio:  input.LOD95LOD50Ratio,
			Retests:          input.Retests,
			Z:                input.Z,
			LimitMin:         -100,
			LimitMax:         500,
		}
	}

	// Parallel execution using worker pool pattern
	rdests := make([]float64, input.NBS)

	// Create job and result channels
	jobs := make(chan int, input.NBS)
	results := make(chan struct {
		index int
		value float64
		err   error
	}, input.NBS)

	// Start worker goroutines
	var wg sync.WaitGroup
	for w := 0; w < input.Threads; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				rd, err := RiskDays(argsList[i])
				results <- struct {
					index int
					value float64
					err   error
				}{i, rd, err}
			}
		}()
	}

	// Send jobs
	go func() {
		for i := 0; i < input.NBS; i++ {
			jobs <- i
		}
		close(jobs)
	}()

	// Close results channel when all workers are done
	go func() {
		wg.Wait()
		close(results)
	}()

	// Collect results with progress tracking
	completed := 0
	for result := range results {
		if result.err != nil {
			return nil, fmt.Errorf("simulation %d failed: %w", result.index, result.err)
		}
		rdests[result.index] = result.value
		completed++

		// Call progress callback if provided (callback decides when to send updates)
		if progressCallback != nil {
			progressCallback(completed, input.NBS)
		}
	}

	// Final progress update
	if progressCallback != nil {
		progressCallback(input.NBS, input.NBS)
	}

	// Calculate statistics
	rdRange := [2]float64{Min(rdests), Max(rdests)}
	rdCrI := [2]float64{
		Quantile(rdests, input.Alpha/2),
		Quantile(rdests, 1-input.Alpha/2),
	}

	// Calculate point estimate based on method
	var rdPE float64
	switch input.PointEstimate {
	case "primary parameters":
		// Use primary parameters (not sampled values)
		primaryParams := RiskDaysInnerParams{
			CopiesPerVirion:  input.CopiesPerVirion,
			C0:               input.C0,
			DoublingTime:     input.DoublingTime,
			VolumeTransfused: input.VolumeTransfused,
			K:                input.K,
			PoolSize:         input.PoolSize,
			LOD50:            input.LOD50,
			LOD95LOD50Ratio:  input.LOD95LOD50Ratio,
			Retests:          input.Retests,
			Z:                input.Z,
			LimitMin:         -100,
			LimitMax:         500,
		}
		var err error
		rdPE, err = RiskDays(primaryParams)
		if err != nil {
			return nil, fmt.Errorf("failed to calculate primary parameters estimate: %w", err)
		}
	case "median":
		rdPE = Median(rdests)
	case "mean":
		rdPE = Mean(rdests)
	case "mode":
		rdPE = KDEModeLog(rdests, 1_000_000, 0, input.Threads)
	default:
		return nil, fmt.Errorf("unknown point estimate method: %s", input.PointEstimate)
	}

	// Return results
	out := &RiskDaysOutput{
		Version:          Version,
		PointEstimate:    rdPE,
		CredibleInterval: rdCrI,
		Range:            rdRange,
		Simulations:      rdests,
	}
	if input.ReturnParams {
		out.Ks = ks
		out.DoublingTimes = doublingTimes
		out.LOD50s = lod50s
		out.VolumesTransfused = volumesTransfused
	}
	return out, nil
}

// riskDaysBSPrep runs the PrEP breakthrough infection bootstrap simulation.
func riskDaysBSPrep(input RiskDaysInput, rng *RandomGenerator, ks, doublingTimes, lod50s, volumesTransfused []float64, progressCallback ProgressCallback) (*RiskDaysOutput, error) {
	// Generate PrEP-specific samples
	var setPoints []float64
	if input.SetPointDistUniform[0] != 0 || input.SetPointDistUniform[1] != 0 {
		setPoints = rng.GenerateUniform(input.SetPointDistUniform[0], input.SetPointDistUniform[1], input.NBS)
	} else {
		setPoints = make([]float64, input.NBS)
		for i := range setPoints {
			setPoints[i] = input.SetPoint
		}
	}

	var eclipses []float64
	if input.EclipseDistUniform[0] != 0 || input.EclipseDistUniform[1] != 0 {
		eclipses = rng.GenerateUniform(input.EclipseDistUniform[0], input.EclipseDistUniform[1], input.NBS)
	} else {
		eclipses = make([]float64, input.NBS)
		for i := range eclipses {
			eclipses[i] = input.Eclipse
		}
	}

	// Sinusoidal oscillation params: fixed at the scalar A/B unless a uniform
	// range is given ([0,0] = fixed). offset is never varied.
	var aVals []float64
	if input.ADistUniform[0] != 0 || input.ADistUniform[1] != 0 {
		aVals = rng.GenerateUniform(input.ADistUniform[0], input.ADistUniform[1], input.NBS)
	} else {
		aVals = make([]float64, input.NBS)
		for i := range aVals {
			aVals[i] = input.A
		}
	}

	var bVals []float64
	if input.BDistUniform[0] != 0 || input.BDistUniform[1] != 0 {
		bVals = rng.GenerateUniform(input.BDistUniform[0], input.BDistUniform[1], input.NBS)
	} else {
		bVals = make([]float64, input.NBS)
		for i := range bVals {
			bVals[i] = input.B
		}
	}

	// Build PrEP args
	argsList := make([]PrepInnerParams, input.NBS)
	for i := 0; i < input.NBS; i++ {
		argsList[i] = PrepInnerParams{
			CopiesPerVirion:  input.CopiesPerVirion,
			C0:               input.C0,
			DoublingTime:     doublingTimes[i],
			SetPoint:         setPoints[i],
			Eclipse:          eclipses[i],
			A:                aVals[i],
			B:                bVals[i],
			Offset:           input.Offset,
			VolumeTransfused: volumesTransfused[i],
			K:                ks[i],
			PoolSize:         input.PoolSize,
			LOD50:            lod50s[i],
			LOD95LOD50Ratio:  input.LOD95LOD50Ratio,
			Retests:          input.Retests,
			SerMin:           input.SerMin,
			SerMax:           input.SerMax,
			SerAlpha:         input.SerAlpha,
			SerBeta:          input.SerBeta,
			Z:                input.Z,
			LimitMin:         -100,
			LimitMax:         500,
		}
	}

	// Parallel execution
	rdests := make([]float64, input.NBS)

	jobs := make(chan int, input.NBS)
	results := make(chan struct {
		index int
		value float64
	}, input.NBS)

	var wg sync.WaitGroup
	for w := 0; w < input.Threads; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				rd := RiskDaysPrep(argsList[i])
				results <- struct {
					index int
					value float64
				}{i, rd}
			}
		}()
	}

	go func() {
		for i := 0; i < input.NBS; i++ {
			jobs <- i
		}
		close(jobs)
	}()

	go func() {
		wg.Wait()
		close(results)
	}()

	completed := 0
	for result := range results {
		rdests[result.index] = result.value
		completed++
		if progressCallback != nil {
			progressCallback(completed, input.NBS)
		}
	}

	if progressCallback != nil {
		progressCallback(input.NBS, input.NBS)
	}

	// Statistics
	rdRange := [2]float64{Min(rdests), Max(rdests)}
	rdCrI := [2]float64{
		Quantile(rdests, input.Alpha/2),
		Quantile(rdests, 1-input.Alpha/2),
	}

	// Point estimate
	var rdPE float64
	switch input.PointEstimate {
	case "primary parameters":
		primaryParams := PrepInnerParams{
			CopiesPerVirion:  input.CopiesPerVirion,
			C0:               input.C0,
			DoublingTime:     input.DoublingTime,
			SetPoint:         input.SetPoint,
			Eclipse:          input.Eclipse,
			A:                input.A,
			B:                input.B,
			Offset:           input.Offset,
			VolumeTransfused: input.VolumeTransfused,
			K:                input.K,
			PoolSize:         input.PoolSize,
			LOD50:            input.LOD50,
			LOD95LOD50Ratio:  input.LOD95LOD50Ratio,
			Retests:          input.Retests,
			SerMin:           input.SerMin,
			SerMax:           input.SerMax,
			SerAlpha:         input.SerAlpha,
			SerBeta:          input.SerBeta,
			Z:                input.Z,
			LimitMin:         -100,
			LimitMax:         500,
		}
		rdPE = RiskDaysPrep(primaryParams)
	case "median":
		rdPE = Median(rdests)
	case "mean":
		rdPE = Mean(rdests)
	case "mode":
		rdPE = KDEModeLog(rdests, 1_000_000, 0, input.Threads)
	default:
		return nil, fmt.Errorf("unknown point estimate method: %s", input.PointEstimate)
	}

	out := &RiskDaysOutput{
		Version:          Version,
		PointEstimate:    rdPE,
		CredibleInterval: rdCrI,
		Range:            rdRange,
		Simulations:      rdests,
	}
	if input.ReturnParams {
		out.Ks = ks
		out.DoublingTimes = doublingTimes
		out.LOD50s = lod50s
		out.VolumesTransfused = volumesTransfused
		out.SetPoints = setPoints
		out.Eclipses = eclipses
		out.As = aVals
		out.Bs = bVals
	}
	return out, nil
}
