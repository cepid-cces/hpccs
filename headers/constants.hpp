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

#ifndef SRC_HEADERS_CONTANTS_HPP_
#define SRC_HEADERS_CONTANTS_HPP_

//  Constants from Handbook of Chemistry and Physics, 70th Edition

const double pi = 3.14159265358979323846;
const double cang = (180.0 / pi);

const double xe = 1.60217733e-19;
const double xk = 1.380658e-23;
const double xn = 6.0221367e23;
const double xeo = 8.854187817e-12;
const double xmv = 0.02241410;

//  Minimum value of (1-cosX). This quantity determines the maximum
//  impact parameter at each velocity. Default value is 0.0005.
const double cmin = 0.0005;
const double dtsf1 = 0.5;
const double dtsf2 = 0.1;
const double sw1 = 0.00005;
const double sw2 = 0.005;

// Lennard-Jones scaling parameters
const double eo = 1.34e-03 * xe;
const double ro = 3.043 * 1.0e-10;
const double ro2 = ro * ro;

// Mobility constant
const double dens = xn / xmv;

/*
 *  inwr is the number of integration steps before the program tests
 *  to see if the trajectory is done or lost. ifail is the number of
 *  failed trajectories that are permitted (a failed trajectory is one
 *  that does not conserve energy to within 1%. Default values are:
 *  inwr = 1        ifail = 100
 *
 */

const int ifail = 100;

const int inwr = 1;

const int EMPTY = -1;
#endif /* SRC_HEADERS_CONTANTS_HPP_ */
