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
	} else {
		return nil, fmt.Errorf("no valid k distribution specified")
	}

	doublingTimes := rng.GenerateTruncatedNormal(input.DoublingTime, input.DoublingTimeNormSD, input.NBS)
	lod50s := rng.GenerateTruncatedNormal(input.LOD50, input.LOD50SD, input.NBS)
	volumesTransfused := rng.GenerateUniform(input.VolumeTransfusedMin, input.VolumeTransfusedMax, input.NBS)

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
		rdPE = KDEModeLog(rdests, 0, 0, input.Threads)
	default:
		return nil, fmt.Errorf("unknown point estimate method: %s", input.PointEstimate)
	}

	// Return results
	return &RiskDaysOutput{
		Version:          Version,
		PointEstimate:    rdPE,
		CredibleInterval: rdCrI,
		Range:            rdRange,
		Simulations:      rdests,
	}, nil
}
