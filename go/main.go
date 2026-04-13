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

package main

import (
	"encoding/json"
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/vitalant-research-institute/residualrisk/riskdays"
)

const helpText = `riskdays_go - Residual HIV Transfusion Transmission Risk Estimation Tool

Usage:
  riskdays_go [input.json]   Read parameters from a JSON file
  riskdays_go                Read parameters from stdin (JSON)

Options:
  -h, --help    Show this help message and exit

The tool accepts a JSON object with simulation parameters and writes results
to stdout as JSON. Progress updates are written to stderr.

Example (stdin):
  echo '{"doubling_time": 0.85, "lod50": 2.73, "pool_size": 16, "n_bs": 10000}' | riskdays_go

Example (file):
  riskdays_go input.json

See README.md for the full JSON parameter schema.
`

func main() {
	// Read input from stdin or file
	var inputData []byte
	var err error

	if len(os.Args) > 1 {
		arg := os.Args[1]
		if arg == "--help" || arg == "-h" {
			fmt.Print(helpText)
			os.Exit(0)
		}
		// Read from file if provided as argument
		inputData, err = os.ReadFile(arg)
		if err != nil {
			fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to read input file: %v"}`+"\n", err)
			os.Exit(1)
		}
	} else {
		// Read from stdin
		inputData, err = io.ReadAll(os.Stdin)
		if err != nil {
			fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to read input from stdin: %v"}`+"\n", err)
			os.Exit(1)
		}
	}

	// Parse input JSON
	var input riskdays.RiskDaysInput
	if err := json.Unmarshal(inputData, &input); err != nil {
		fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to parse input JSON: %v"}`+"\n", err)
		os.Exit(1)
	}

	// Set default threads if not specified
	if input.Threads == 0 {
		input.Threads = runtime.NumCPU() - 1
		if input.Threads < 1 {
			input.Threads = 1
		}
	}

	// Create progress callback that writes to stderr
	// Only send updates when percentage changes to reduce output
	var lastPercent int
	progressCallback := func(completed, total int) {
		newPercent := (completed * 100) / total
		// Update only when percentage changes (reduces warnings in Python)
		if completed == 1 || newPercent > lastPercent {
			lastPercent = newPercent
			progress := riskdays.ProgressMessage{
				Type:      "progress",
				Completed: completed,
				Total:     total,
				Percent:   float64(completed) / float64(total),
			}
			progressJSON, _ := json.Marshal(progress)
			fmt.Fprintf(os.Stderr, "%s\n", progressJSON)
		}
	}

	// Run risk days calculation
	output, err := riskdays.RiskDaysBS(input, progressCallback)
	if err != nil {
		fmt.Fprintf(os.Stderr, `{"type": "error", "message": "calculation failed: %v"}`+"\n", err)
		os.Exit(1)
	}

	// Write output to stdout as JSON
	outputJSON, err := json.Marshal(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to marshal output: %v"}`+"\n", err)
		os.Exit(1)
	}

	fmt.Println(string(outputJSON))
}
