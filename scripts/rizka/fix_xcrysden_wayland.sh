#!/bin/bash

# ==============================================================================
# XCrysDen Wayland/Xwayland Fix
# ==============================================================================
# Problem:
# When running XCrysDen on modern Linux systems (especially Wayland/Xwayland)
# or via containers like Distrobox, you may encounter the following error:
#
#   "Couldn't configure togl widget"
#
# Cause:
# By default, the Togl (OpenGL Tcl/Tk) widget in XCrysDen requests an
# "Accumulation Buffer" when setting up the visual context. Modern Wayland
# compositors and Xwayland often do not provide this buffer, causing the widget
# initialization to fail instantly and the application to crash.
#
# Fix:
# This script disables the accumulation buffer requirement by creating or
# modifying the XCrysDen user configuration file (~/.xcrysden/custom-definitions)
# and setting `toglOpt(accum)` to `false`.
# ==============================================================================

echo "Fixing XCrysDen Togl widget error for Wayland/Xwayland..."

# Create the user configuration directory if it doesn't exist
mkdir -p ~/.xcrysden

CONFIG_FILE="$REPO_ROOT/.xcrysden/custom-definitions"

# If the file exists, we ensure the line is uncommented or present
if [ -f "$CONFIG_FILE" ]; then
  # Check if the line is commented out and uncomment it
  if grep -q "#set toglOpt(accum)  false" "$CONFIG_FILE"; then
    sed -i 's/#set toglOpt(accum)  false/set toglOpt(accum)  false/' "$CONFIG_FILE"
    echo "Updated existing configuration in $CONFIG_FILE."
  elif ! grep -q "set toglOpt(accum)[[:space:]]*false" "$CONFIG_FILE"; then
    # If it's not in the file at all, append it
    echo "" >>"$CONFIG_FILE"
    echo "# Fix for Wayland Togl widget crash" >>"$CONFIG_FILE"
    echo "set toglOpt(accum) false" >>"$CONFIG_FILE"
    echo "Appended configuration to $CONFIG_FILE."
  else
    echo "Fix is already applied in $CONFIG_FILE."
  fi
else
  # If the file does not exist, check if we can copy it from a system install
  if [ -f "/usr/share/xcrysden/Tcl/custom-definitions" ]; then
    cp "/usr/share/xcrysden/Tcl/custom-definitions" "$CONFIG_FILE"
    sed -i 's/#set toglOpt(accum)  false/set toglOpt(accum)  false/' "$CONFIG_FILE"
    echo "Copied system configuration and applied fix."
  # Or from a local directory install
  elif [ -f "./Tcl/custom-definitions" ]; then
    cp "./Tcl/custom-definitions" "$CONFIG_FILE"
    sed -i 's/#set toglOpt(accum)  false/set toglOpt(accum)  false/' "$CONFIG_FILE"
    echo "Copied local configuration and applied fix."
  else
    # Otherwise, just create a minimal config file
    cat <<'EOF' >"$CONFIG_FILE"
# Custom settings to fix Wayland/Xwayland Togl widget crash
set toglOpt(accum) false
EOF
    echo "Created minimal configuration file with fix at $CONFIG_FILE."
  fi
fi

echo "Done!"
