#!/bin/bash
g++ -g -o GE11.o GE11.cpp $(root-config --cflags --glibs --ldflags) -lMinuit
./GE11.o
