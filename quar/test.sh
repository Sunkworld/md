#!/bin/bash
gfortran random.f90 mol.phys.13.f90 -o test.x; ./test.x
gfortran correlation.f90 -o l.x; ./l.x
notifyme.py
