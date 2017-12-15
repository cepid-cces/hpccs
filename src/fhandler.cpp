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

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <headers/constants.hpp>
#include <headers/fhandler.hpp>
#include <headers/atomMLJ.hpp>
#include <headers/globals.hpp>
#include <headers/molecule.hpp>
#include <cmath>
#if defined(__INTEL_COMPILER)
#include <aligned_new>
#else
#include <new>
#endif
using namespace std;

int fcoordinates(const char *filename){

	/*
	 * This function will read the PQR file and fill the Atom
	 * objects creating a list with ATOMs.
	 * Status = 0 - File not Found
	 * Status = 1 - Cant Open File
	 * Status = 2 - Success the vectors are filled
	 */

	const int atomsMLJCount = 110;
	std::string string, recordName, atomName,residueName;
	ifstream infile;
	AtomMLJ *atomMLJArray __attribute__((aligned(16)));
	int status, serial, chainId, j = 0;
	int numberOfConformations;
	double x,y, z, charges, radius;
	double mass, epsilon, sigma;
	unsigned int i = 0;
	atomMLJArray = new AtomMLJ[atomsMLJCount];
	m2 = 0.0;

	/*
	 * Read the AtomsMLJ.csv file and create the
	 * list of atom with mass and LJ parameters.
	 *
	 */

	infile.open("config/config.in");

	if (infile.is_open()){
		while(getline(infile,string)){
			istringstream sstream(string);
			sstream >> numberOfConformations >> itn >> inp >> imp >> ipr >> temperature >> gas;
		}
		infile.close();
	}

	if(gas == 2){
		infile.open("config/AtomsMLJN2.csv");
	}else {
		infile.open("config/AtomsMLJHe.csv");
	}


	if (infile.is_open()){
		while(getline(infile,string)){
			istringstream sstream(string);
			sstream >> atomName >> mass >> epsilon >> sigma;
			atomMLJArray[i].setAtomName(atomName);
			atomMLJArray[i].setMass(mass);
			atomMLJArray[i].setSigma(sigma);
			atomMLJArray[i].setEpsilon(epsilon);
			i++;
		}
		infile.clear();
		infile.close();
	}

	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
	for (int atoms = 0; atoms < i; atoms++){
		atomMLJArray[atoms].setEox4(4.0 * atomMLJArray[atoms].getEpsilon());
		atomMLJArray[atoms].setRo6lj(pow(atomMLJArray[atoms].getSigma(),6));
		atomMLJArray[atoms].setRo12lj(pow(atomMLJArray[atoms].getSigma(),12));
		atomMLJArray[atoms].setDro6(6.0 * atomMLJArray[atoms].getRo6lj());
		atomMLJArray[atoms].setDro12(12.0 * atomMLJArray[atoms].getRo12lj());
	}


	/* Check the number of atoms */
	i = 0;
	infile.open (filename);
			if (infile.is_open()){
				while(getline(infile,string)) // To get you all the lines.
				{
					istringstream sstream(string);
					sstream >> recordName;

					if (recordName == "ATOM"){
						++i;
					}
				}
				//infile.clear();
				infile.close();
			}

	numberOfAtoms = i;
	molecule = new Molecule(numberOfAtoms);

	i = 0;
	romax = 0.0;
	infile.open(filename);
		if (infile.is_open()){
			while(getline(infile,string)) // To get you all the lines.
			{
				istringstream sstream(string);
				sstream >> recordName >> serial >> atomName >> residueName >> chainId >> x >> y >> z >> charges >> radius;

				if (recordName == "ATOM"){
							atomName.erase(std::remove_if(atomName.begin(), atomName.end(), [](char x){return std::isdigit(x);}), atomName.end());
							molecule->fx[i] = x;
							molecule->fy[i] = y;
							molecule->fz[i] = z;
							molecule->charge[i] = charges;
							for (j = 0; j < atomsMLJCount; j++){
								if (atomName.length() > 1 && isupper(atomName.at(1))){
										atomName = atomName.at(0);
								}
								if (atomMLJArray[j].getAtomName() == atomName){
									molecule->xmass[i] = atomMLJArray[j].getMass();
									molecule->eox4[i] = atomMLJArray[j].getEox4();
									molecule->ro6lj[i] = atomMLJArray[j].getRo6lj();
									molecule->ro12lj[i] = atomMLJArray[j].getRo12lj();
									molecule->dro6[i] = atomMLJArray[j].getDro6();
									molecule->dro12[i] = atomMLJArray[j].getDro12();
									if(atomMLJArray[j].getSigma() > romax) romax = atomMLJArray[j].getSigma();
									break;
								}
							}
							++i;
				}
			}
			infile.clear();
			infile.close();
			status = 2;
		}else{
			status = 1;
		}

		if (status != 1){
		#if defined(__INTEL_COMPILER)
    	#pragma vector aligned
		#endif
		for (i = 0; i < numberOfAtoms; i++){
			fxo += molecule->fx[i] * molecule->xmass[i];
			fyo += molecule->fy[i] * molecule->xmass[i];
			fzo += molecule->fz[i] * molecule->xmass[i];
			m2 += molecule->xmass[i];
		}

		fxo /= m2;
		fyo /= m2;
		fzo /= m2;

		#if defined(__INTEL_COMPILER)
    	#pragma vector aligned
		#endif
		for (i = 0; i < numberOfAtoms; i++){
			molecule->fx[i] = (molecule->fx[i] - fxo) * 1.0e-10;
			molecule->ox[i] = molecule->fx[i];
			molecule->fy[i] = (molecule->fy[i] - fyo) * 1.0e-10;
			molecule->oy[i] = molecule->fy[i];
			molecule->fz[i] = (molecule->fz[i] - fzo) * 1.0e-10;
			molecule->oz[i] = molecule->fz[i];
		}
		}

	delete[] atomMLJArray;
	return status;
}


