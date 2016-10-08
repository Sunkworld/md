gfortran random.f90 mol.phys.13.f90 -o test.x; ./test.x
gfortran correlation.f90 -o l.x; ./l.x
mv cor_energy.dat cor_energy_2.dat
rm *.txt
touch cortime
awk '{sum+=$3} END{print sum}' <cor_energy_0.5.dat>>cortime
awk '{sum+=$3} END{print sum}' <cor_energy_1.dat>>cortime
awk '{sum+=$3} END{print sum}' <cor_energy_2.dat>>cortime
awk '{sum+=$3} END{print sum}' <cor_energy_3.dat>>cortime
awk '{sum+=$3} END{print sum}' <cor_energy_5.dat>>cortime
awk '{sum+=$3} END{print sum}' <cor_energy_10.dat>>cortime
awk '{sum+=$3} END{print sum}' <cor_energy_15.dat>>cortime
paste cor_energy_0.5.dat cor_energy_1.dat cor_energy_2.dat cor_energy_3.dat cor_energy_5.dat cor_energy_10.dat cor_energy_15.dat > tmp
awk '{print $2,$3,$7,$11,$15,$19,$23,$27}' <tmp >all
