#!/bin/bash

cd ./src/pesgenerator
g++ -std=c++17 -Wall --pedantic-errors -Wpedantic getpes.cc -lboost_system -lboost_filesystem -o run
./run
rm -f run
cd ..
cd ..


