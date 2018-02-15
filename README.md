High Performance Collision Cross Section Calculation – HPCCS

This program is licensed granted by STATE UNIVERSITY OF CAMPINAS – UNICAMP (the “University”) for use of HIGH PERFORMANCE COLLISION CROSS SECTION - HPCCS (“the Software”) through this website https://github.com/cepid-cces/hpccs (the ”Website”). By downloading the Software through the Website, you (the “Licensee”) are confirming that you agree that your use of the Software is subject to the academic license terms. For more information about HPCCS please contact: skaf@iqm.unicamp.br (Munir Skaf) or leandro.zanotto@gmail.com (Leandro Zanotto).

Citation: Zanotto, L., Heerdt, G., Souza, P. C. T., Araujo, G. & Skaf, M. S. (in press). High Performance Collision Cross Section Calculation - HPCCS. Journal of Computational Chemistry, 2018, DOI: 10.1002/jcc.25199.

•	How To Compile HPCCS:
Modify the HPCCS source directory in the Makefile (line 3):
INC_DIR = $HOME/HPCCS
On source folder just enter make to compile.
Works with Intel and GCC compilers. Its highly recommended to use GCC 7.2 for better performance due to memory alignments.

•	Configuration:
In HPCCS configuration file description, all parameters must be splitted by tabs in one single line (config/config.in).

1. Number of Conformations -> I need to specify how many conformations the input file contains (generally 1). 
2. Number of complete cycles for average mobility calculation in MOBIL2. Default value is 10 if this field use a number < 1 or is not set.
3. Number of points in velocity integration in MOBIL2. Default value is 20 if     this field use a number < 1 or is not set.
4. Number of points in Monte Carlo integrations of impact parameter     and orientation in MOBIL2. Default value is 500 if this field use a number < 1 or is not set.
5. Number of rotations in average potential calculation. If ipr=0 an unrotated configuration is used. Otherwise ipr random rotations are employed. Default value is 1000 if it is not set.
6. Temperature, default value is 300 K if it is not set.
7. Which gas: 1 for Helium, 2 - For Nitrogen.

The file config.in will be like this:
1       10      20      500     1000    298.0   2



•	Running:
To run just enter ./hpccs filnename.pqr

It is not possible to move the hpccs file to another directory because it needs the configuration and other files in the source directory.

•	Making the PQR File
PDB2PQR is a Python software package that provides a utility for converting protein files in PDB format to PQR format. For more information about PDB2PQR tool you can visit: http://pdb2pqr.sourceforge.net/.

H++ is a server that given a (PDB) structure file on input, H++ outputs the completed structure in several common formats (PDB, PQR, AMBER inpcrd/prmtop) and provides a set of tools for analysis of electrostatic-related molecular properties. For more information about H++ server see: http://biophysics.cs.vt.edu/.



