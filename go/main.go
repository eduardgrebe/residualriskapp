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
	"encoding/binary"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/vitalant-research-institute/residualrisk/riskdays"
)

const helpText = `riskdays_go - Residual HIV Transfusion Transmission Risk Estimation Tool

Usage:
  riskdays_go [input.json]        Read parameters from a JSON file
  riskdays_go                     Read parameters from stdin (JSON)
  riskdays_go --kde-mode [f.json] KDE mode estimation (see below)
  riskdays_go --hsm-mode [f.json] Half-Sample Mode estimation (see below)

Options:
  -h, --help       Show this help message and exit
  -v, --version    Print version and exit
  --kde-mode       Run KDE mode estimation instead of bootstrap simulation
  --hsm-mode       Run Half-Sample Mode estimation instead of bootstrap simulation

KDE mode:
  Reads JSON {"data":[...], "n_grid":5000, "cap":5000, "threads":0} from
  stdin (or a file argument) and writes {"mode": X} to stdout.
  n_grid=0 → auto-size; cap=0 → no cap; threads=0 → use all CPU cores.

HSM mode:
  Reads JSON {"data":[...]} from stdin (or a file argument) and writes
  {"mode": X} to stdout. Uses the Half-Sample Mode algorithm (bandwidth-free,
  outlier-robust). No tuning parameters needed.

Example (stdin):
  echo '{"doubling_time": 0.85, "lod50": 2.73, "pool_size": 16, "n_bs": 10000}' | riskdays_go

Example (file):
  riskdays_go input.json

See README.md for the full JSON parameter schema.
`

// kdeModeInput is the JSON input schema for the --kde-mode subcommand.
type kdeModeInput struct {
	Data    []float64 `json:"data"`
	NGrid   int       `json:"n_grid"`
	Cap     int       `json:"cap"`
	Threads int       `json:"threads"`
}

// kdeModeOutput is the JSON output schema for the --kde-mode subcommand.
type kdeModeOutput struct {
	Mode float64 `json:"mode"`
}

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
		if arg == "--version" || arg == "-v" {
			fmt.Println(riskdays.Version)
			os.Exit(0)
		}
		if arg == "--kde-mode" {
			// KDE mode subcommand: read JSON from file arg or stdin
			if len(os.Args) > 2 {
				inputData, err = os.ReadFile(os.Args[2])
				if err != nil {
					fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to read kde-mode input file: %v"}`+"\n", err)
					os.Exit(1)
				}
			} else {
				inputData, err = io.ReadAll(os.Stdin)
				if err != nil {
					fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to read kde-mode input from stdin: %v"}`+"\n", err)
					os.Exit(1)
				}
			}
			var kdeInput kdeModeInput
			if err := json.Unmarshal(inputData, &kdeInput); err != nil {
				fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to parse kde-mode JSON: %v"}`+"\n", err)
				os.Exit(1)
			}
			mode := riskdays.KDEModeLog(kdeInput.Data, kdeInput.NGrid, kdeInput.Cap, kdeInput.Threads)
			out, _ := json.Marshal(kdeModeOutput{Mode: mode})
			fmt.Println(string(out))
			os.Exit(0)
		}
		if arg == "--hsm-mode" {
			// Half-Sample Mode subcommand: read JSON from file arg or stdin
			if len(os.Args) > 2 {
				inputData, err = os.ReadFile(os.Args[2])
				if err != nil {
					fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to read hsm-mode input file: %v"}`+"\n", err)
					os.Exit(1)
				}
			} else {
				inputData, err = io.ReadAll(os.Stdin)
				if err != nil {
					fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to read hsm-mode input from stdin: %v"}`+"\n", err)
					os.Exit(1)
				}
			}
			var hsmInput struct {
				Data []float64 `json:"data"`
			}
			if err := json.Unmarshal(inputData, &hsmInput); err != nil {
				fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to parse hsm-mode JSON: %v"}`+"\n", err)
				os.Exit(1)
			}
			mode := riskdays.HalfSampleMode(hsmInput.Data)
			out, _ := json.Marshal(kdeModeOutput{Mode: mode})
			fmt.Println(string(out))
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

	if input.ReturnParams {
		// Binary wire format:
		//   [8 bytes]  uint64 LE — byte length of JSON header
		//   [N bytes]  JSON header with summary stats + column metadata
		//   [rest]     column-major float64 LE arrays:
		//              iwp, k, doubling_time, lod50, volume_transfused
		columns := []string{"iwp", "k", "doubling_time", "lod50", "volume_transfused"}
		type binaryHeader struct {
			Version          string     `json:"version"`
			PointEstimate    float64    `json:"point_estimate"`
			CredibleInterval [2]float64 `json:"credible_interval"`
			Range            [2]float64 `json:"range"`
			NBS              int        `json:"n_bs"`
			Columns          []string   `json:"columns"`
		}
		header := binaryHeader{
			Version:          output.Version,
			PointEstimate:    output.PointEstimate,
			CredibleInterval: output.CredibleInterval,
			Range:            output.Range,
			NBS:              len(output.Simulations),
			Columns:          columns,
		}
		headerBytes, err := json.Marshal(header)
		if err != nil {
			fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to marshal binary header: %v"}`+"\n", err)
			os.Exit(1)
		}
		// Write header length as uint64 LE
		if err := binary.Write(os.Stdout, binary.LittleEndian, uint64(len(headerBytes))); err != nil {
			fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to write header length: %v"}`+"\n", err)
			os.Exit(1)
		}
		// Write JSON header
		if _, err := os.Stdout.Write(headerBytes); err != nil {
			fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to write header: %v"}`+"\n", err)
			os.Exit(1)
		}
		// Write each column contiguously (column-major)
		for _, col := range [][]float64{
			output.Simulations,
			output.Ks,
			output.DoublingTimes,
			output.LOD50s,
			output.VolumesTransfused,
		} {
			if err := binary.Write(os.Stdout, binary.LittleEndian, col); err != nil {
				fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to write binary column: %v"}`+"\n", err)
				os.Exit(1)
			}
		}
	} else {
		// Standard JSON output (backward compatible)
		outputJSON, err := json.Marshal(output)
		if err != nil {
			fmt.Fprintf(os.Stderr, `{"type": "error", "message": "failed to marshal output: %v"}`+"\n", err)
			os.Exit(1)
		}
		fmt.Println(string(outputJSON))
	}
}
