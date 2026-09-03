#!/bin/bash
# Script to run convergence for K-points

template_name=$1
num_cores=$2
seq_start=${3:-6}
seq_step=${4:-2}
seq_end=${5:-24}
skip_existing=false

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

echo "Starting K-point convergence test: from $seq_start to $seq_end (step: $seq_step)"

if [ "$skip_existing" = "true" ]; then
  echo "Mode: SKIP existing files."
fi
echo "--------------------------------------------------"

for kpoint in $(seq "$seq_start" "$seq_step" "$seq_end"); do
  base_name="${template_name%.*}"
  output_file="${base_name}_kpoint_${kpoint}.out"
  input_file="${base_name}_kpoint_${kpoint}.in"

  # The skip check
  if [ "$skip_existing" = "true" ] && [ -f "$output_file" ]; then
    echo "Skipping kpoint grid = $kpoint x $kpoint x 1 (output already exists)."
    continue
  fi

  kpoint_line="${kpoint} ${kpoint} 1 0 0 0"
  echo "Running kpoint grid = $kpoint x $kpoint x 1"

  sed "s/^[[:space:]]*[0-9]\+[[:space:]]\+[0-9]\+[[:space:]]\+1 0 0 0.*/$kpoint_line/" "$template_name" >"$input_file"

  mpirun --use-hwthread-cpus -np "$num_cores" pw.x <"$input_file" >"$output_file"
done

echo "--------------------------------------------------"
echo "All jobs done!"
