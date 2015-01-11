#!/bin/bash
cd "$( dirname "$0" )"
find "$(cd ../..;pwd)" -name '*.h' | xargs -l dirname | sort -u | xargs -l -- echo -I > includes.txt
cat includes.txt
