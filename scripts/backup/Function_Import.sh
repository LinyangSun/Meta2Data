#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source all .sh files in this folder
for f in "$SCRIPT_DIR"/*.sh; do
    if file "$f" | grep -q "CRLF"; then
         dos2unix "$f"
    fi
    chmod +x "$f"
    [ "$f" != "$BASH_SOURCE" ] && source "$f"
done

# Make py_16s.py executable and add to PATH
dos2unix "$SCRIPT_DIR/py_16s.py"
chmod +x "$SCRIPT_DIR/py_16s.py"
export PATH="$SCRIPT_DIR:$PATH"