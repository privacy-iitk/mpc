#!/bin/bash

N=$1
if [ -z "$N" ]; then
  echo "Usage: $0 <parameter>"
  exit 1
fi

./remise p0 "$N" &
PID0=$!

./remise p1 "$N" &
PID1=$!

wait $PID0 $PID1

