#!/bin/bash
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.1d0/h=0.2d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.2d0/h=0.3d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.3d0/h=0.4d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.4d0/h=0.5d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.5d0/h=0.6d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.6d0/h=0.7d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.7d0/h=0.8d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.8d0/h=0.9d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x
sed -i '' 's/h=0.9d0/h=1.0d0/g' mol.phys.13.f90
gfortran random.f90 mol.phys.13.f90 -o test.x;./test.x

