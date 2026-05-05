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

import "fmt"

// RiskDaysInput represents all input parameters for the risk days bootstrap calculation
type RiskDaysInput struct {
	// Primary parameters
	K                   float64   `json:"k"`
	DoublingTime        float64   `json:"doubling_time"`
	DoublingTimeNormSD  float64   `json:"doubling_time_norm_sd"`
	LOD50               float64   `json:"lod50"`
	LOD50SD             float64   `json:"lod50_sd"`
	LOD95LOD50Ratio     float64   `json:"lod95_lod50_ratio"`
	VolumeTransfused    float64   `json:"volume_transfused"`
	VolumeTransfusedMin float64   `json:"volume_transfused_min"`
	VolumeTransfusedMax float64   `json:"volume_transfused_max"`
	PoolSize            int       `json:"pool_size"`
	Retests             int       `json:"retests"`

	// Optional parameters with defaults
	C0               float64   `json:"c0"`                // default: 0.00025
	CopiesPerVirion  int       `json:"copies_per_virion"` // default: 2
	Alpha            float64   `json:"alpha"`             // default: 0.05
	Z                float64   `json:"z"`                 // default: 1.6449

	// K distribution parameters (one of these groups must be provided)
	KPosteriorSample []float64 `json:"k_posterior_sample,omitempty"`
	KGammaShape      *float64  `json:"k_gamma_shape,omitempty"`  // Deprecated: use KInvGammaAlpha/Beta
	KGammaScale      *float64  `json:"k_gamma_scale,omitempty"`  // Deprecated: use KInvGammaAlpha/Beta
	KInvGammaAlpha   *float64  `json:"k_invgamma_alpha,omitempty"` // shape (α)
	KInvGammaBeta    *float64  `json:"k_invgamma_beta,omitempty"`  // scale (β, same as scipy's scale)

	// Lognormal mixture k distribution
	// Two-component mixture: w * LN(mu1, sigma1) + (1-w) * LN(mu2, sigma2)
	// mu/sigma are log-scale parameters (scipy: lognorm(s=sigma, scale=exp(mu)))
	KLnMixW      *float64  `json:"k_lnmix_w,omitempty"`       // weight of component 1 (0–1)
	KLnMixMu1    *float64  `json:"k_lnmix_mu1,omitempty"`     // log-mean of component 1
	KLnMixSigma1 *float64  `json:"k_lnmix_sigma1,omitempty"`  // log-sd of component 1
	KLnMixMu2    *float64  `json:"k_lnmix_mu2,omitempty"`     // log-mean of component 2
	KLnMixSigma2 *float64  `json:"k_lnmix_sigma2,omitempty"`  // log-sd of component 2

	// Simulation parameters
	NBS           int    `json:"n_bs"`            // default: 10000
	Seed          int64  `json:"seed"`            // default: 126887
	Threads       int    `json:"threads"`         // default: num_cpu - 1
	PointEstimate string `json:"point_estimate"`  // "primary parameters", "median", "mean", "mode"
	ModePrecision int    `json:"mode_precision"`  // default: 2
}

// RiskDaysOutput represents the output of the risk days calculation
type RiskDaysOutput struct {
	Version           string     `json:"version"`
	PointEstimate     float64    `json:"point_estimate"`
	CredibleInterval  [2]float64 `json:"credible_interval"` // [lower, upper]
	Range             [2]float64 `json:"range"`             // [min, max]
	Simulations       []float64  `json:"simulations"`
}

// ProgressMessage represents a progress update during calculation
type ProgressMessage struct {
	Type      string  `json:"type"`      // "progress"
	Completed int     `json:"completed"`
	Total     int     `json:"total"`
	Percent   float64 `json:"percent"`
}

// RiskDaysInnerParams contains all parameters needed for a single risk days calculation
type RiskDaysInnerParams struct {
	CopiesPerVirion  int
	C0               float64
	DoublingTime     float64
	VolumeTransfused float64
	K                float64
	PoolSize         int
	LOD50            float64
	LOD95LOD50Ratio  float64
	Retests          int
	Z                float64
	LimitMin         float64
	LimitMax         float64
}

// SetDefaults sets default values for optional parameters
func (input *RiskDaysInput) SetDefaults() {
	if input.C0 == 0 {
		input.C0 = 0.00025
	}
	if input.CopiesPerVirion == 0 {
		input.CopiesPerVirion = 2
	}
	if input.Alpha == 0 {
		input.Alpha = 0.05
	}
	if input.Z == 0 {
		input.Z = 1.6449
	}
	if input.NBS == 0 {
		input.NBS = 10000
	}
	if input.Seed == 0 {
		input.Seed = 126887
	}
	if input.PointEstimate == "" {
		input.PointEstimate = "primary parameters"
	}
	if input.ModePrecision == 0 {
		input.ModePrecision = 2
	}
}

// Validate checks that input parameters are valid
func (input *RiskDaysInput) Validate() error {
	if input.NBS <= 0 {
		return fmt.Errorf("n_bs must be greater than zero")
	}
	if input.PoolSize < 1 {
		return fmt.Errorf("pool_size must be at least 1")
	}
	if input.Retests < 0 {
		return fmt.Errorf("retests must be non-negative")
	}
	if input.KPosteriorSample == nil &&
		(input.KGammaShape == nil || input.KGammaScale == nil) &&
		(input.KInvGammaAlpha == nil || input.KInvGammaBeta == nil) &&
		(input.KLnMixW == nil || input.KLnMixMu1 == nil || input.KLnMixSigma1 == nil ||
			input.KLnMixMu2 == nil || input.KLnMixSigma2 == nil) {
		return fmt.Errorf("a k distribution must be provided: k_posterior_sample, both k_gamma parameters, both k_invgamma parameters, or all k_lnmix parameters")
	}
	return nil
}
