#!/bin/bash
set -e

#run CSI
python3 /scripts/csi/main.py --csv=csi_out.csv --hdf5=csi_out.h5 "${@:1}"

#create webapp folder inside the results folder
cp -r /scripts/html/ html/

#make webapp port of results
python3 /scripts/csi/webapp-json.py csi_out.h5 > html/results.json

#parse networks (MAP and marginal)
python3 /scripts/csi/getnetwork.py --map csi_out.h5 > csi_MAP.csv
python3 /scripts/csi/getnetwork.py --marginal csi_out.h5 > csi_marginal.csv