#!/bin/bash
rm main
rm -r data/*
cd build
cmake ..
time make
cd ..
time aprun ./main config.lua
#time mpirun -np 2 ./main
#time mpirun -np 2 ./main config.lua

