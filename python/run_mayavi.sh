#! /usr/bin/env bash
for i in `seq 0 22`
do
    python mayavi_plot.py $i
done
