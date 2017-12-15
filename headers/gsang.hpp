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


#ifndef SRC_HEADERS_GSANG_HPP_
#define SRC_HEADERS_GSANG_HPP_

double gsangHe(double v,double x,double rnt, double rnp, double rng, unsigned int totalAtoms);

double gsangN2(double v,double x,double rnt, double rnp, double rng, unsigned int totalAtoms);

double gsangHe(const double v,const double x);

double gsangN2(const double v,const double x);
#endif /* SRC_HEADERS_GSANG_HPP_ */
