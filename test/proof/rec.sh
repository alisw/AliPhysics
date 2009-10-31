#!/bin/bash

[ -z $1 ] && { echo "Usage: rec.sh <run_number>"; exit 1; }

aliroot -q run.C\($1\)
