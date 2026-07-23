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

// PrepInnerParams contains all parameters needed for a single PrEP risk days calculation.
// This parallels RiskDaysInnerParams but includes PrEP-specific fields for
// the 3-phase viral dynamics and serology non-detection.
type PrepInnerParams struct {
	CopiesPerVirion  int
	C0               float64
	DoublingTime     float64
	SetPoint         float64
	Eclipse          float64
	A                float64
	B                float64
	Offset           float64
	DrugEffect       float64
	VolumeTransfused float64
	K                float64
	PoolSize         int
	LOD50            float64
	LOD95LOD50Ratio  float64
	Retests          int
	SerMin           float64
	SerMax           float64
	SerAlpha         float64
	SerBeta          float64
	Z                float64
	LimitMin         float64
	LimitMax         float64
}
