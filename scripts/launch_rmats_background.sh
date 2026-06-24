#!/bin/bash
# launch_rmats_background.sh
nohup bash scripts/10_rmats.sh > logs/rmats_final.log 2>&1 &
echo "rMATS launched with PID $!"
