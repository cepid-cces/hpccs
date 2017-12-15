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

#ifndef SRC_HEADERS_ATOM_HPP_
#define SRC_HEADERS_ATOM_HPP_

class Molecule {


	public:
		double *fx __attribute__((aligned(16)));
		double *fy __attribute__((aligned(16)));
		double *fz __attribute__((aligned(16)));
		double *ox __attribute__((aligned(16)));
		double *oy __attribute__((aligned(16)));
		double *oz __attribute__((aligned(16)));
		double *charge __attribute__((aligned(16)));
		double *xmass __attribute__((aligned(16)));
		double *eox4 __attribute__((aligned(16)));
		double *ro6lj __attribute__((aligned(16)));
		double *ro12lj __attribute__((aligned(16)));
		double *dro6 __attribute__((aligned(16)));
		double *dro12 __attribute__((aligned(16)));
		//Constructor

		Molecule(unsigned int size);
		~Molecule();

};


#endif /* SRC_HEADERS_ATOM_HPP_ */
