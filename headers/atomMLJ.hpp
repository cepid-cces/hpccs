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


#ifndef SRC_HEADERS_ATOMMLJ_HPP_
#define SRC_HEADERS_ATOMMLJ_HPP_

#include <string>
using namespace std;
class AtomMLJ {

	private:
		string atomName;
		double mass, eox4, sigma, epsilon;
		double ro6lj,ro12lj,dro6,dro12;


	public:
		//Constructor
		AtomMLJ();

		~AtomMLJ();

		//Getters and Setters
		void setAtomName(string atomName);
		string  getAtomName(){ return atomName;}

		void setMass(double s_mass);
		double getMass(){ return mass;}

		void setSigma(double s_sigma);
		double getSigma(){ return sigma;}

		void setEpsilon(double s_epsilon);
		double getEpsilon(){ return epsilon;}

		void setEox4(double s_eox4);
		double  getEox4(){ return eox4;}

		void setRo6lj(double s_ro6lj);
		double  getRo6lj(){ return ro6lj;}

		void setRo12lj(double s_ro12lj);
		double  getRo12lj(){ return ro12lj;}

		void setDro6(double s_dro6);
		double  getDro6(){ return dro6;}

		void setDro12(double s_dro12);
		double  getDro12(){ return dro12;}

};

#endif /* SRC_HEADERS_ATOMMLJ_HPP_ */
