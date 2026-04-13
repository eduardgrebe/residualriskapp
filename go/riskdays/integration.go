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
	"gonum.org/v1/gonum/integrate/quad"
)

// RiskDays performs numerical integration to calculate risk days for a single set of parameters
// Corresponds to Python _risk_days()
func RiskDays(params RiskDaysInnerParams) (float64, error) {
	// Create integrand function
	integrand := func(t float64) float64 {
		prob, err := ProbInfectiousNondetection(t, params)
		if err != nil {
			// In case of error, return 0
			// This shouldn't happen in normal operation
			return 0
		}
		return prob
	}

	// Perform numerical integration using adaptive quadrature
	// Equivalent to scipy.integrate.quad(integrand, limits[0], limits[1])
	result := quad.Fixed(integrand, params.LimitMin, params.LimitMax, 1000, nil, 0)

	return result, nil
}
