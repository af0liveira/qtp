#!/bin/bash

# Set up the correct Python interpreter.
# py3="$HOME/opt/anaconda3/bin/python3"
py3="python3"

# Set up the path to the QTP packager directory
srcpath="$HOME/00-PARA/02-Areas/Code/qtp/qtp"

# Run the QTP program.
# !!! You should not change it. !!!
exec $py3 "$srcpath"/qtp.py "$@"
