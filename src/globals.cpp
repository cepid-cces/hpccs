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

#include <headers/globals.hpp>

// Mass constants
double m2;

// Mobility constants modified once a long the software
double mconst;
double mu;

//Atoms
Molecule *molecule;

double romax;

//Angles
double theta, phi, agamma;

double pot, dpotx, dpoty, dpotz, dmax;

int ifailc;

int samples;

double mob;

double cs;

double hvar, hcvar;

unsigned int numberOfAtoms;

double temperature;

int itn;

int inp;

int imp;

int ipr;

int gas;

double dipol;

// Mass constants
double m1;

int trajlost;

double sdevpc;

double fxo, fyo, fzo;
