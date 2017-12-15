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


#ifndef SRC_HEADERS_DIFFEQ_DERIV_HPP_
#define SRC_HEADERS_DIFFEQ_DERIV_HPP_

double diffeq(int l, double tim, double dt, double w[6], double dw[10],  double array[6][6], double q[6], double *fx, double *fy, double *fz, double *phvar, double *phcvar);

void deriv(double *w,double *dw, double *fx, double *fy, double *fz);

double diffeq(int l, double tim, double dt, double w[6], double dw[10], double array[6][6], double q[6]);

void deriv(double *w, double *dw);

#endif /* SRC_HEADERS_DIFFEQ_DERIV_HPP_ */
