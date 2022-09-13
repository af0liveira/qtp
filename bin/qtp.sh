#!/bin/bash

# Set up the correct Python interpreter, if necessary
py3="python3"

# Set up the path to the QTP package directory
srcpath="$HOME/qtp/qtp"

# Run the QTP program.
# !!! You should not change it. !!!
exec $py3 "$srcpath"/qtp.py "$@"
