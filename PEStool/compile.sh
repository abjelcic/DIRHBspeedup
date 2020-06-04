#!/bin/bash


#generate adequate dirhb.par file
cd ./src/pesgenerator
g++ -Wall --pedantic-errors -Wpedantic pargenerator.cc -o run
./run
rm -f run
cd ..
cd ..


#compile dirhbt code
cd ./src/dirhbt
make clean
make run
rm -f *.par
cd ..
cd ..


#generate (beta,gamma) folders
cd ./src/pesgenerator
g++ -Wall --pedantic-errors -Wpedantic taskgenerator.cc -lboost_system -lboost_filesystem -o run
./run
rm -f run
cd ..
cd ..


#cleaning phase
cd ./src/dirhbt
make clean
cd ..
cd ..
