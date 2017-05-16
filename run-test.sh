#!/bin/bash

#
# run-test.sh
#
# Tester for randomness_test Library
# by snovvcrash
# 04.2017
#

# shopt -s nullglob

GOOD=0
BAD=0

if [[ "$#" -ne 2 ]]; then
	echo "help: ./run-test.sh <./executable> <path_to_folder>"
	exit
fi

for f in "$2"/*.txt; do
	echo -n "Processing $f... "
	RET=$($1 $f)
	echo "$RET"
	if [[ "$RET" == "good" ]]; then
		((GOOD++))
	else
		((BAD++))
	fi
done

echo "GOOD: $GOOD"
echo "BAD:  $BAD"
