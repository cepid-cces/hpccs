MobcalHPC Configuration file description

All parameters must be splited by tabs in one single line

1 - Number of Conformations -> I need to specify how many conformations the input file contains.
2 - Number of complete cycles for average mobility calculation in MOBIL2. Default value is 10 if 
    this field use a number < 1 or is not set.
3 - Number of points in velocity integration in MOBIL2. Default value is 20 if 
    this field use a number < 1 or is not set.
4 - Number of points in Monte Carlo integrations of impact parameter
    and orientation in MOBIL2. Default value is 500 if 
    this field use a number < 1 or is not set.
5 - Number of rotations in average potential calculation. If ipr=0
    an unrotated configuration is used. Otherwise ipr random rotations
    are employed. Default value is 1000 if it is not set.
6-  Temperature -> default value is 300K if it is not set.
7 - Which gas 1 for Helium, 2 - For Nitrogen, 3 - for CO2

=============================Compiling==================================
On source forder just enter make to compile.

Works with Intel and GCC.

========================================================================


=================Running===================

To run just enter ./hpccs filnename.pqr

===========================================

