#!/bin/bash

template_name=$1
num_cores=$2
seq_start=${3:-40}
seq_step=${4:-2}
seq_end=${5:-60}

skip_existing=false

# Scan all passed arguments for the "--skip" flag
for arg in "$@"; do
  if [ "$arg" == "--skip" ]; then
    skip_existing=true
  fi
done

if [ -z "$template_name" ] || [ -z "$num_cores" ]; then
  echo "Error: Missing arguments."
  echo "Usage: $0 <template_file.in> <number_of_cores> [start] [step] [end] [--skip]"
  exit 1
fi

echo "Starting convergence test: ecutwfc from $seq_start to $seq_end (step: $seq_step)"

if [ "$skip_existing" = "true" ]; then
  echo "Mode: SKIP existing files."
fi
echo "--------------------------------------------------"

for ecut in $(seq "$seq_start" "$seq_step" "$seq_end"); do
  base_name="${template_name%.*}"
  output_file="${base_name}_ecut_${ecut}.out"
  input_file="${base_name}_ecut_${ecut}.in"

  # The skip check
  if [ "$skip_existing" = "true" ] && [ -f "$output_file" ]; then
    echo "Skipping ecutwfc = $ecut (output already exists)."
    continue
  fi

  # PAW Multiplier
  # ecutrho=$((ecut * 8))

  echo "Running ecutwfc = $ecut Ry | ecutrho = $ecutrho Ry"

  sed -e "s/ecutwfc[[:space:]]*=[[:space:]]*[0-9]*/ecutwfc = $ecut/g" \
    "$template_name" >"$input_file" # -e "s/ecutrho[[:space:]]*=[[:space:]]*[0-9]*/ecutrho = $ecutrho/g" \

  mpirun --use-hwthread-cpus -np "$num_cores" pw.x <"$input_file" >"$output_file"
done

echo "--------------------------------------------------"
echo "All jobs done!"
