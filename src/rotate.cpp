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

#include <headers/constants.hpp>
#include <headers/globals.hpp>
#include <cmath>

void rotate(double rnt, double rnp, double rng, double *fx, double *fy, double *fz){

//  Rotates the cluster/molecule.
    double ogamma = 0.0, ophi = 0.0, otheta = 0.0, ntheta = 0.0, ngamma = 0.0, nphi = 0.0, rxy, rtheta, rphi, rgamma;
    unsigned int iatom;

	rtheta = rnt * 2.0 * pi;
	rphi = asin((rnp * 2.0)-1.0) + (pi/2.0);
	rgamma = rng * 2.0 * pi;

	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (iatom = 0; iatom < numberOfAtoms; iatom++){

    	rxy = sqrt((molecule->ox[iatom] * molecule->ox[iatom]) + (molecule->oy[iatom] * molecule->oy[iatom]));
        if(rxy == 0.0) {
        	fx[iatom] = cos(ntheta) * rxy;
        	fy[iatom] = sin(ntheta) * rxy;
        } else {
            otheta = acos(molecule->ox[iatom] / rxy);
            if(molecule->oy[iatom] < 0.0){
                otheta = (2.0 * pi) - otheta;
            }
            ntheta = otheta + rtheta;
            fx[iatom] = cos(ntheta) * rxy;
            fy[iatom] = sin(ntheta) * rxy;
        }

        rxy = sqrt((molecule->oz[iatom] * molecule->oz[iatom]) + (fy[iatom] * fy[iatom]));
        if(rxy == 0.0){
        	 fz[iatom] = cos(nphi) * rxy;
        	 fy[iatom] = sin(nphi) * rxy;
        }else {
           ophi = acos(molecule->oz[iatom] / rxy);
           if(fy[iatom] < 0.0){
              ophi=(2.0 * pi) - ophi;
           }
            nphi = ophi + rphi;
            fz[iatom] = cos(nphi) * rxy;
            fy[iatom] = sin(nphi) * rxy;
        }

        rxy = sqrt((fx[iatom] *  fx[iatom]) + (fy[iatom] *  fy[iatom]));
        if(rxy == 0.0){
          fx[iatom] = cos(ngamma) * rxy;
          fy[iatom] = sin(ngamma) * rxy;
        } else {
           ogamma = acos(fx[iatom] / rxy);
           if(fy[iatom] < 0.0){
               ogamma = (2.0 * pi) - ogamma;
    	   }
           ngamma = ogamma + rgamma;
           fx[iatom] = cos(ngamma) * rxy;
           fy[iatom] = sin(ngamma) * rxy;
        }
    }
}


void rotate(){

//  Rotates the cluster/molecule.
    double ogamma, ophi, otheta, ntheta, ngamma, nphi, rxy;
    double rnt, rnp, rng;
    unsigned int iatom;

    ophi = 0.0;
    nphi = 0.0;
    ngamma = 0.0;
    ntheta = 0.0;

	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (iatom = 0; iatom < numberOfAtoms; iatom++){

    	rxy = sqrt((molecule->ox[iatom] * molecule->ox[iatom]) + (molecule->oy[iatom] * molecule->oy[iatom]));

        if(rxy == 0.0) {
        	molecule->fx[iatom] = cos(ntheta) * rxy;
        	molecule->fy[iatom] = sin(ntheta) * rxy;
        } else {
            otheta = acos(molecule->ox[iatom] / rxy);
            if(molecule->oy[iatom] < 0.0){
                otheta = (2.0 * pi) - otheta;
            }
            ntheta = otheta + theta;
            molecule->fx[iatom] = cos(ntheta) * rxy;
            molecule->fy[iatom] = sin(ntheta) * rxy;
        }

        rxy = sqrt((molecule->oz[iatom] * molecule->oz[iatom]) + (molecule->fy[iatom] * molecule->fy[iatom]));

        if(rxy == 0.0){
        	molecule->fz[iatom] = cos(nphi) * rxy;
        	molecule->fy[iatom] = sin(nphi) * rxy;
        }else {
           ophi = acos(molecule->oz[iatom] / rxy);
           if(molecule->fy[iatom] < 0.0){
              ophi=(2.0 * pi) - ophi;
           }
            nphi = ophi + phi;
            molecule->fz[iatom] = (cos(nphi) * rxy);
            molecule->fy[iatom] = (sin(nphi) * rxy);
        }

        rxy = sqrt((molecule->fx[iatom] * molecule->fx[iatom]) + (molecule->fy[iatom] * molecule->fy[iatom]));

        if(rxy == 0.0){
          molecule->fx[iatom] = (cos(ngamma) * rxy);
          molecule->fy[iatom] = (sin(ngamma) * rxy);
        } else {
           ogamma = acos(molecule->fx[iatom] / rxy);
           if(molecule->fy[iatom] < 0.0){
               ogamma = (2.0 * pi) - ogamma;
    	   }
           ngamma = ogamma + agamma;
           molecule->fx[iatom] = (cos(ngamma) * rxy);
           molecule->fy[iatom] = (sin(ngamma) * rxy);
        }
    }
}
