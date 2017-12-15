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
#if defined(__INTEL_COMPILER)
#include <aligned_new>
#else
#include <new>
#endif
Molecule::Molecule(unsigned int size){
			fx = new double[size];
			fy = new double[size];
			fz = new double[size];
			ox = new double[size];
			oy = new double[size];
			oz = new double[size];
			charge = new double[size];
			xmass = new double[size];
			eox4 = new double[size];
			ro6lj = new double[size];
			ro12lj = new double[size];
			dro6 = new double[size];
			dro12 = new double[size];
		};

		Molecule::~Molecule(){
			delete[] fx;
			delete[] fy;
			delete[] fz;
			delete[] ox;
			delete[] oy;
			delete[] oz;
			delete[] charge;
			delete[] xmass;
			delete[] eox4;
			delete[] ro6lj;
			delete[] ro12lj;
			delete[] dro6;
			delete[] dro12;
		};
