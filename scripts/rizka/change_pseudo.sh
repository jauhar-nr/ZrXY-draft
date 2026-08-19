#!/bin/bash

# Check if new directory argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <new_pseudo_dir>"
  exit 1
fi

NEW_DIR="$1"

echo "Updating pseudo_dir in all .in files recursively to: $NEW_DIR"

# Find all .in files and update the pseudo_dir variable
find . -type f -name "*.in" | while read -r file; do
  sed -i -E "s|(pseudo_dir[[:space:]]*=[[:space:]]*['\"]).*(['\"])|\1$NEW_DIR\2|g" "$file"
  echo "  Updated: $file"
done

echo "=== Done! ==="
