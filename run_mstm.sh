

rm *.mod *.o
#rm amn_temp.dat
#gfortran-mp-6 -o mstm.x mpidefs-serial.f90 intrinsics.f90 modules.f90 main.f90
#
./mstm-scan.x test3.inp
