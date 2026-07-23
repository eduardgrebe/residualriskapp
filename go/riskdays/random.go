// Residual HIV Transfusion Transmission Risk Estimator
// Copyright (C) 2025-2026 Vitalant and Eduard Grebe Consulting
// Author: Eduard Grebe <egrebe@vitalant.org> <eduard@grebe.consulting>

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

package riskdays

import (
	"math"
	"math/rand"

	exprand "golang.org/x/exp/rand"
	"gonum.org/v1/gonum/stat/distuv"
)

// RandomGenerator holds the random number generator state
type RandomGenerator struct {
	rng    exprand.Source
	stdRng *rand.Rand // For methods that need standard library rand
}

// NewRandomGenerator creates a new random generator with specified seed
func NewRandomGenerator(seed int64) *RandomGenerator {
	return &RandomGenerator{
		rng:    exprand.NewSource(uint64(seed)),
		stdRng: rand.New(rand.NewSource(seed)),
	}
}

// GenerateTruncatedNormal generates samples from a normal distribution truncated at 0
// (positive values only). Equivalent to
// scipy.stats.truncnorm.rvs(a=-mean/sd, b=inf, loc=mean, scale=sd, size=n) — scipy's
// a/b bounds are standardised (SDs from loc), so a lower bound of 0 maps to
// a=(0-mean)/sd, NOT a=0 (which would truncate at the mean).
// This implementation uses rejection sampling
func (rg *RandomGenerator) GenerateTruncatedNormal(mean, sd float64, n int) []float64 {
	samples := make([]float64, n)

	// Degenerate distribution (sd<=0): there is no spread, so distuv.Normal.Rand()
	// always returns `mean` and the rejection loop below would spin forever when
	// mean<=0. Mirror the Python _sample_positive_normal helper and fill with the
	// mean. (Non-positive lod50 / doubling_time should be rejected upstream by
	// Validate(); see the core.py:737 follow-up.)
	if sd <= 0 {
		for i := range samples {
			samples[i] = mean
		}
		return samples
	}

	normal := distuv.Normal{
		Mu:    mean,
		Sigma: sd,
		Src:   rg.rng,
	}

	// Backstop against a non-terminating rejection loop when mean is far below 0
	// (positive draws virtually impossible). For valid inputs (mean>0) the first
	// draw is positive, so this is draw-for-draw identical to the original loop.
	const maxAttempts = 10000
	for i := 0; i < n; i++ {
		sample := normal.Rand()
		for attempt := 1; sample <= 0; attempt++ { // truncated at 0
			if attempt >= maxAttempts {
				sample = math.Abs(mean)
				break
			}
			sample = normal.Rand()
		}
		samples[i] = sample
	}

	return samples
}

// GenerateUniform generates samples from a uniform distribution
// Equivalent to np.random.uniform(low, high, n)
func (rg *RandomGenerator) GenerateUniform(low, high float64, n int) []float64 {
	samples := make([]float64, n)
	uniform := distuv.Uniform{
		Min: low,
		Max: high,
		Src: rg.rng,
	}

	for i := 0; i < n; i++ {
		samples[i] = uniform.Rand()
	}

	return samples
}

// GenerateInvGamma generates samples from an Inverse Gamma distribution.
// alpha is the shape parameter; beta is the scale parameter (equivalent to
// scipy's invgamma(a=alpha, scale=beta)).
//
// Uses the relationship: if X ~ Gamma(alpha, rate=beta) then 1/X ~ InvGamma(alpha, scale=beta).
// In Gonum's Gamma, Beta is the rate (1/scale), so setting Beta=beta gives the correct Gamma.
func (rg *RandomGenerator) GenerateInvGamma(alpha, beta float64, n int) []float64 {
	samples := make([]float64, n)
	// Gonum Gamma: Beta = rate = 1/scale. Setting rate=beta yields scale=1/beta.
	// Then 1/X ~ InvGamma(alpha, scale=beta).
	gamma := distuv.Gamma{
		Alpha: alpha,
		Beta:  beta, // rate parameter in Gonum = 1/scale
		Src:   rg.rng,
	}
	for i := 0; i < n; i++ {
		samples[i] = 1.0 / gamma.Rand()
	}
	return samples
}

// GenerateLogNormalMixture generates samples from a two-component lognormal mixture.
//
// Each sample is drawn from component 1 (weight w, LN(mu1, sigma1)) with probability w,
// or from component 2 (weight 1-w, LN(mu2, sigma2)) otherwise.
//
// mu and sigma are the log-scale mean and standard deviation (parameterisation matches
// scipy.stats.lognorm(s=sigma, scale=exp(mu)) and numpy.random.lognormal(mean=mu, sigma)).
// Gonum's LogNormal uses Mu=mu, Sigma=sigma — same convention.
func (rg *RandomGenerator) GenerateLogNormalMixture(w, mu1, sigma1, mu2, sigma2 float64, n int) []float64 {
	samples := make([]float64, n)
	comp1 := distuv.LogNormal{Mu: mu1, Sigma: sigma1, Src: rg.rng}
	comp2 := distuv.LogNormal{Mu: mu2, Sigma: sigma2, Src: rg.rng}
	uniform := distuv.Uniform{Min: 0, Max: 1, Src: rg.rng}
	for i := 0; i < n; i++ {
		if uniform.Rand() < w {
			samples[i] = comp1.Rand()
		} else {
			samples[i] = comp2.Rand()
		}
	}
	return samples
}

// GenerateGamma generates samples from a gamma distribution
// Equivalent to np.random.gamma(shape, scale, n)
func (rg *RandomGenerator) GenerateGamma(shape, scale float64, n int) []float64 {
	samples := make([]float64, n)

	// Gonum's Gamma uses (alpha, beta) where beta = 1/scale
	gamma := distuv.Gamma{
		Alpha: shape,
		Beta:  1.0 / scale,
		Src:   rg.rng,
	}

	for i := 0; i < n; i++ {
		samples[i] = gamma.Rand()
	}

	return samples
}

// BootstrapChoice performs bootstrap sampling (sampling with replacement)
// Equivalent to np.random.choice(array, size=n, replace=True)
func (rg *RandomGenerator) BootstrapChoice(array []float64, n int) []float64 {
	samples := make([]float64, n)
	arrLen := len(array)

	for i := 0; i < n; i++ {
		idx := rg.stdRng.Intn(arrLen)
		samples[i] = array[idx]
	}

	return samples
}

// NormalCDF computes the cumulative distribution function of the standard normal distribution
// Equivalent to scipy.stats.norm.cdf(x)
func NormalCDF(x float64) float64 {
	return 0.5 * (1 + math.Erf(x/math.Sqrt2))
}

// NormalCDFWithParams computes the CDF of a normal distribution with given mean and std
func NormalCDFWithParams(x, mean, std float64) float64 {
	if std == 0 {
		if x < mean {
			return 0
		}
		return 1
	}
	z := (x - mean) / std
	return NormalCDF(z)
}
