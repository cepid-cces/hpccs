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
#include <headers/globals.hpp>
#include <headers/constants.hpp>
#include <headers/fhandler.hpp>
#include <headers/mobil2.hpp>
#include <headers/molecule.hpp>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cmath>
using namespace std;
//Paramenter 1 -> Trajectorie Number
//Parameter 2 -> File name
int main (int argc, char **argv){

	double start = omp_get_wtime();
	int status, k;
	double percent;

	samples = 1;

	// Constant for ion-induced dipole potential
	// Mass of helium atom = 4.0026d0
	// Polarizability of helium = 0.204956d-30 m3
	// xeo is permitivity of vacuum, 8.854187817d-12 F.m-1

	if (gas == 2){
		dipol = (1.641e-30 / (2.0 * 4.0 * pi * xeo)) * (xe * xe);
		// Mass constants
		m1 = 28.0134;

	}else if (gas == 3){

	}else{
		dipol = (0.204956e-30 / (2.0 * 4.0 * pi * xeo)) * (xe * xe);
		// Mass constants
		m1 = 4.0026;
	}

	/*
	 * Call the fhandler to Fill the atoms list
	 */
	status = fcoordinates(argv[1]);
		/*
		 * Read the user input PQR file and populate
		 * the Atom list with the right values
		 * for Mass] Sigma] Epsilon] Charge] x] y] z.
		 *
		 */

	mconst = sqrt(18.0 * pi) / 16.0;
	mconst *= sqrt(xn * 1.0e3) * sqrt((1.0 / m1) + (1.0 / m2));
	mconst *= xe / sqrt(xk);
	mconst /= dens;
	mu = ((m1 * m2) / (m1 + m2)) / (xn * 1.0e3);

    if (status == 2){
  	  for (int i = 0; i < samples; i++){

  		  mobil2();
  		  //print out summary



  		  percent = (trajlost * 100.0) / (itn * inp * imp);

  	//	  cout << "Input file name -> " << "\t" << 	filename;
  	//	  cout << "Using a uniform charge distribution" << "\n";
  	//	  cout << "Using a calculated (non-uniform) charge distribution" << "\n";
  	//	  cout << "Using no charge - only LJ interactions" << "\n";
  	//	  cout << "Temperature" << t << "\n";
  	//	  cout << "Using Mersene Twister with seed integer" << "\n";
  	//	  cout << "Mobility calculation by MOBIL2 (trajectory method)" << "\n";
  	//	  cout << sw1 << "\t" << sw2 << "\t" << dtsf1 << "\t" << dtsf2 << "\n";
  	//	  cout << "Number of complete cycles (itn) = " << "\t" << itn << "\n";
  	//	  cout << "Number of velocity points (inp) = " << "\t" << inp << "\n";
  	//	  cout << "Number of random points (imp) = " << "\t" << imp << "\n";
  	//	  cout << "Total number of points = " << "\t" << itn * inp * imp << "\n";
  		  cout << "\n" << "=================================================" << "\n";
  		  cout << "===   Zanotto, L. et al. (2017) submitted to  ===" << "\n";
  		  cout << "===   Journal of Computational Chemistry.     ===" << "\n";
  		  cout << "=================================================" << "\n";
  		  cout << "coordinate set = " << i << "\n";
  		  cout << "Average (second order) TM mobility = " << "\t" << scientific << mob << "\n";
  		  cout << "Inverse average (second order) # TM mobility = " << "\t" << scientific << 1.0/mob << "\n";
  		  cout << "Average TM cross section = " << scientific << cs * 1.e20 << "\n";
  		  cout << "Standard deviation (percent) = " << fixed << sdevpc << "\n";
  		  cout << "Total number of Trajectories " << itn * inp * imp << "\n";
  		  cout << fixed << "Number of lost trajectories  = " << trajlost << "\t" << setprecision(2) << " = " << percent << "%" << "\n";
  		  cout << "Summary" << "Program -> " << argv[1] <<  "\n";
  		  double end = omp_get_wtime();
  		  cout << "Total Time: " << (end - start) << "s";
  		 ofstream myfile;
  		 myfile.open ("output.out");

  		 myfile << "=================================================" << "\n";
  		 myfile << "===   Zanotto, L. et al. (2017) submitted to  ===" << "\n";
  		 myfile << "===   Journal of Computational Chemistry.     ===" << "\n";
  		 myfile << "=================================================" << "\n";
  		 myfile << "coordinate set = " << i << "\n";
  		 myfile << "Average (second order) TM mobility = " << "\t" << scientific << mob << "\n";
  		 myfile << "Inverse average (second order) # TM mobility = " << "\t" << scientific << 1.0/mob << "\n";
  		 myfile << "Average TM cross section = " << scientific << cs * 1.e20 << "\n";
  		 myfile << "Standard deviation (percent) = " << fixed << sdevpc << "\n";
  		 myfile << "Total number of Trajectories " << itn * inp * imp << "\n";
  		 myfile << fixed << "Number of lost trajectories  = " << trajlost << "\t" << setprecision(2) << " = " << percent << "%" << "\n";
  		 myfile << "Total Time: " <<  setprecision(2) << end - start << "s";
  		 myfile.close();
  	  }
    }

    delete molecule;

    return 0;
  }



