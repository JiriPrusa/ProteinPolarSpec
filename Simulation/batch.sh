#!/bin/bash

for i in {41..50}
do
  mkdir $i
  cd $i
  cp ../chainSelector.py .
  cp ../dipoleReporter.py .
  cp ../3_production_meta.py .
  cp ../jesslys_let.pdb .
  cp ../state_number_${i}.xml .
  run_openmm_dip.sh $i 168
  cd ..
done
