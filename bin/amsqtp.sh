#!/bin/bash

# If AMS is correctly configured, this doesn't need to be changed
py3="$AMSBIN/amspython"

# Set up the path to the QTP package directory
srcpath="$HOME/qtp/qtp"

# Run the QTP program.
# !!! You should not change it. !!!
exec $py3 "$srcpath"/qtp.py "$@"
