#!/bin/bash

./phbuilder.py

while read initialLambda; do

    sed -i "/initial_lambda/c\initial_lambda        = ${initialLambda}" constant_ph_input.dat
    ./run.sh &> /dev/null
    ./calibrate.py ${initialLambda}

done < initials.txt
