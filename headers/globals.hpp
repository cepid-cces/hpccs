/*
 * This program is licensed granted by STATE UNIVERSITY OF CAMPINAS – UNICAMP (the “University”)
 * for use of HIGH PERFORMANCE COLLISION CROSS SECTION - HPCCS (“the Software”) through this website
 * https://github.com/cepid-cces/hpccs (the ”Website”).
 *
 * By downloading the Software through the Website, you (the “Licensee”) are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about HPCCS please contact: skaf@iqm.unicamp.br (Munir Skaf)
 * or leandro.zanotto@gmail.com (Leandro Zanotto).
 *
 */

#include <headers/molecule.hpp>

#ifndef SRC_HEADERS_GLOBALS_HPP_
#define SRC_HEADERS_GLOBALS_HPP_

// Mass constants
extern double m2;

// Mobility constants modified once a long the software
extern double mconst;
extern  double mu;

//Atoms
extern Molecule *molecule;

extern  double romax;

//Angles
extern double theta, phi, agamma;

extern int ifailc;

extern int samples;

extern double mob;

extern double cs;

extern  double hvar, hcvar;

extern  unsigned int numberOfAtoms;

// Temperature
extern double temperature;

// Number of complete cycles for average mobility calculation in
// MOBIL2. Default value is 10.
extern int itn;

// Number of points in velocity integration in MOBIL2. Default
// value is 20.
extern int inp;

// Number of points in Monte Carlo integrations of impact parameter
// and orientation in MOBIL2. Default value is 500.
extern  int imp;

//  Define parameters for POTENT
//  Number of rotations in average potential calculation. If ipr=0
//  an unrotated configuration is used. Otherwise ipr random rotations
//  are employed.
extern int ipr;

//Which gas 1 for Helium, 2 - For Nitrogen, 3 - for CO2
extern  int gas;

// Constant for ion-induced dipole potential
// Mass of helium atom = 4.0026d0
// Polarizability of helium = 0.204956d-30 m3
// xeo is permitivity of vacuum, 8.854187817d-12 F.m-1
extern  double dipol;

// Mass constants
extern double m1;

extern  int trajlost;

extern double sdevpc;

extern double fxo, fyo, fzo;

#endif /* SRC_HEADERS_GLOBALS_HPP_ */
