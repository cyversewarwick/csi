#!/bin/bash
set -e

#mark start
sleep 5
touch tempfile
sleep 5

#run thing
bash /scripts/csi/csi_wrapper.sh "${@:1}"

#wrap up output and kick out tempfile
find . -mindepth 1 -newer tempfile -exec tar -rf FullOutput.tar {} \;
rm tempfile