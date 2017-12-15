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

#include <string>
#include <headers/atomMLJ.hpp>
using namespace std;

AtomMLJ::AtomMLJ(){ }

AtomMLJ::~AtomMLJ(){ }

// Atom member functions
void AtomMLJ::setAtomName(string s_atomName)
{
	atomName = s_atomName;
}

void AtomMLJ::setMass(double s_mass)
{
	mass = s_mass;
}

void AtomMLJ::setSigma(double s_sigma){

	sigma = s_sigma;

}

void AtomMLJ::setEpsilon(double s_epsilon){

	epsilon = s_epsilon;

}

void AtomMLJ::setEox4(double s_eox4)
{
	eox4 = s_eox4;
}

void AtomMLJ::setRo6lj(double s_ro6lj)
{
	ro6lj = s_ro6lj;
}

void AtomMLJ::setRo12lj(double s_ro12lj)
{
	ro12lj = s_ro12lj;
}

void AtomMLJ::setDro6(double s_dro6)
{
	dro6 = s_dro6;
}

void AtomMLJ::setDro12(double s_dro12)
{
	dro12 = s_dro12;
}
