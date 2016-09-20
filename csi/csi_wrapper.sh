#!/bin/bash
set -e

#tag replicates
python3 /scripts/csi/expression_filter.py --Input $1 --TagReps

#run CSI
python3 /scripts/csi/main.py --csv=csi_out.csv --hdf5=csi_out.h5 FilteredOutput.csv "${@:2}"

#kick tagged replicates file
rm FilteredOutput.csv

#create webapp folder inside the results folder
cp -r /scripts/html/ html/

#make webapp port of results
python3 /scripts/csi/webapp-json.py csi_out.h5 > html/results.json

#parse networks (MAP and marginal)
python3 /scripts/csi/getnetwork.py --map csi_out.h5 > csi_MAP.csv
python3 /scripts/csi/getnetwork.py --marginal csi_out.h5 > csi_marginal.csv