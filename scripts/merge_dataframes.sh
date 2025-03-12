#!/bin/bash

# Function to display usage instructions
usage() {
    >&2 echo "Usage: $0 [delimiter] [output_file] <file1> <file2> [...fileN]"
    >&2 echo "If no delimiter is specified, tab is used by default."
    exit 1
}

# Check if at least one input file is provided
if [[ $# -lt 1 ]]; then
    usage
fi

# Check if the first argument is a valid delimiter
delimiter='\\t'  # Default delimiter (tab)
if [[ "$1" != /* && ! -f "$1" ]]; then
    delimiter="$1"
    shift
fi

# Check if the next argument is an output file or input file
if [[ -f "$1" ]]; then
    output="/dev/stdout"  # No output file specified, use stdout
else
    output="$1"
    shift
fi

# Ensure there is at least one input file remaining
if [[ $# -lt 1 ]]; then
    usage
fi

# Temporary file to store merged output
temp_file=$(mktemp)

header=$(head -1 "$1")
header="sample_id""$delimiter""$header"
echo -e "$header" > "$temp_file"

# Loop over each input file and prepend the file path as the first column
for file in "$@"; do
    sid=$(basename "$file")
    sid=${sid%_*}
    awk -v FS="$delimiter" -v OFS="$delimiter" -v sid="$sid" '
        NR == 1 { next }  # Skip header for all files
        { print sid, $0 } # Add sample id as the first column
    ' "$file" >> "$temp_file"
done 

# Move the result to the output file or print to stdout
mv "$temp_file" "$output"
