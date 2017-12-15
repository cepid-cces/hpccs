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


#ifndef HEADERS_POTENTIALN2_HPP_
#define HEADERS_POTENTIALN2_HPP_

double dljpotN2(const double x, const double y, const double z, double *fx, double *fy, double *fz);

double dljpotN2(const double x, const double y, const double z);

void dljpotN2(const double x, const double y, const double z, double *pot, double *dmax, double *dpotx, double *dpoty, double *dpotz);

void dljpotN2(const double x, const double y, const double z, double *pot,double *dmax, double *dpotx, double *dpoty, double *dpotz,  double *fx, double *fy, double *fz);


#endif /* HEADERS_POTENTIALN2_HPP_ */
