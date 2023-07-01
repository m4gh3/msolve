#!/bin/bash

file=F4SAT-31

$(pwd)/msolve -f input_files/$file.ms -o test/diff/$file.res -S -g 2
if [ $? -gt 0 ]; then
    exit 1
fi

diff test/diff/$file.res output_files/$file.res
if [ $? -gt 0 ]; then
    exit 2
fi

rm test/diff/$file.res
