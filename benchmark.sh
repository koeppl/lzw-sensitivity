#!/usr/bin/env bash

function die {
	echo "$1" >&2
	exit 1
}

if [[ -z "$1" ]]; then
	echo "Usage: $0 <output_path>"
	exit 1
fi
kOutputPath="$1"
[[ -d "$kOutputPath" ]] || mkdir -p "$kOutputPath" || die "Failed to create output directory: $kOutputPath"

kCompressors=(
	"lzmw"
	"lzw"
	"lzd"
)

for order in $(seq 3 15); do
	# Run the Python script to generate the output
	for compressor in "${kCompressors[@]}"; do
		echo "Running compressor: $compressor"
		python3 ${compressor}_generate.py -o "$kOutputPath/${compressor}_${order}" -k "$order" || die "Failed to run the Python script for compressor: $compressor"
		python3 ${compressor}_compress.py -i "$kOutputPath/${compressor}_${order}" || die "Failed to run the Python script for compressor: $compressor"
	done
done
