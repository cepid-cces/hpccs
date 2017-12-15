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
#if defined(__INTEL_COMPILER)
#include <aligned_new>
#else
#include <new>
#endif
double dljpotHe(const double x,const double y,const double z){

//  Subroutine to calculate L-J + ion-dipole potential.
    double rxyz[7] __attribute__((aligned(16)));
    double rx, ry, rz, e00, pot;
    double xx, xx2, yy, yy2, zz, zz2, rxyz3i;

    rx = 0.0;
    ry = 0.0;
    rz = 0.0;
    e00 = 0.0;
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int iatom = 0; iatom < numberOfAtoms; iatom++){
        xx = x - molecule->fx[iatom];
        xx2 = xx * xx;
        yy = y - molecule->fy[iatom];
        yy2 = yy * yy;
        zz = z - molecule->fz[iatom];
        zz2 = zz * zz;
        rxyz[0] = xx2 + yy2 + zz2;
        rxyz[1] = sqrt(rxyz[0]);
        rxyz[2] = rxyz[0] * rxyz[1];
        rxyz[3] = rxyz[2] * rxyz[0];
        rxyz[4] = rxyz[3] * rxyz[1];
        rxyz[5] = rxyz[3] * rxyz[2];
        rxyz[6] = rxyz[4] * rxyz[4];

        // LJ potential
        e00 += (molecule->eox4[iatom] * ((molecule->ro12lj[iatom] / rxyz[6]) - (molecule->ro6lj[iatom] / rxyz[4])));

		// ion-induced dipole potential
		rxyz3i = molecule->charge[iatom] / rxyz[2];
		rx += xx * rxyz3i;
		ry += yy * rxyz3i;
		rz += zz * rxyz3i;
    }
    pot =  e00 - (dipol * ((rx * rx) + (ry * ry) + (rz * rz)));

    return pot;
}


double dljpotHe(const double x,const double y,const double z, double *fx, double *fy, double *fz){

//  Subroutine to calculate L-J + ion-dipole potential.
	double rxyz[7] __attribute__((aligned(16)));
	double rx,ry,rz;
	double xx,xx2,yy, yy2, zz,zz2,rxyz3i;
	double pot, e00 = 0.0;;
    rx = 0.0;
    ry = 0.0;
    rz = 0.0;

	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int iatom = 0; iatom < numberOfAtoms; ++iatom){
        xx = x - fx[iatom];
    	xx2 = xx * xx;
    	yy = y - fy[iatom];
    	yy2 = yy * yy;
    	zz = z - fz[iatom];
    	zz2 = zz * zz;
        rxyz[0] = xx2 + yy2 + zz2;
        rxyz[1] = sqrt(rxyz[0]);
        rxyz[2] = rxyz[0] * rxyz[1];
        rxyz[3] = rxyz[2] * rxyz[0];
        rxyz[4] = rxyz[3] * rxyz[1];
        rxyz[5] = rxyz[3] * rxyz[2];
        rxyz[6] = rxyz[4] * rxyz[4];


        // LJ potential
        e00 += (molecule->eox4[iatom] * ((molecule->ro12lj[iatom] / rxyz[6]) - (molecule->ro6lj[iatom] / rxyz[4])));

		// ion-induced dipole potential
		rxyz3i = molecule->charge[iatom] / rxyz[2];
		rx += xx * rxyz3i;
		ry += yy * rxyz3i;
		rz += zz * rxyz3i;

    }
    pot =  e00 - (dipol * ((rx * rx) + (ry * ry) + (rz * rz)));
    return pot;
}


void dljpotHe(const double x,const double y,const double z, double *pot, double *dmax, double *dpotx, double *dpoty, double *dpotz, double *fx, double *fy, double *fz){

//  Subroutine to calculate L-J + ion-dipole potential.

	double rxyz[8] __attribute__((aligned(16)));
	double  de00, rxyz3i, rxyz5i;
   	double xx, xx2, yy, yy2, zz, zz2;
    double e00, de00x, de00y, de00z, dmaxx, eox4;
    double sum1, sum2, sum3, sum4, sum5, sum6;
    double *rxyz_vec __attribute__((aligned(16)));

    double rx = 0.0, ry = 0.0, rz = 0.0;
    e00 = 0.0;
    de00x = 0.0;
    de00y = 0.0;
    de00z = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    sum6 = 0.0;
    dmaxx = 2.0 * romax;
    rxyz_vec = new double[numberOfAtoms];

	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int iatom = 0; iatom < numberOfAtoms; ++iatom){

		xx = x - fx[iatom];
		xx2 = xx * xx;
		yy = y - fy[iatom];
		yy2 = yy * yy;
		zz = z - fz[iatom];
		zz2 = zz * zz;
        rxyz[0] = xx2 + yy2 + zz2;
        rxyz[1] = sqrt(rxyz[0]);
        rxyz_vec[iatom] = rxyz[1];
        rxyz[2] = rxyz[0] * rxyz[1];
        rxyz[3] = rxyz[2] * rxyz[0];
        rxyz[4] = rxyz[3] * rxyz[1];
        rxyz[5] = rxyz[3] * rxyz[2];
        rxyz[6] = rxyz[4] * rxyz[4];
		rxyz[7] = rxyz[6] * rxyz[0];
		eox4 = molecule->eox4[iatom];
		// LJ potential
		e00 += (eox4 * ((molecule->ro12lj[iatom] / rxyz[6]) - (molecule->ro6lj[iatom] / rxyz[4])));

        // LJ derivative
        de00 = eox4 * ((molecule->dro6[iatom] / rxyz[5]) - (molecule->dro12[iatom] / rxyz[7]));
        de00x += de00 * xx;
        de00y += de00 * yy;
        de00z += de00 * zz;

        // ion-induced dipole potential
		rxyz3i = molecule->charge[iatom] / rxyz[2];
		rxyz5i = -3.0 * molecule->charge[iatom] / rxyz[3];
		rx += xx * rxyz3i;
		ry += yy * rxyz3i;
		rz += zz * rxyz3i;
		//  ion-induced dipole derivative
		sum1 += rxyz3i + (xx2 * rxyz5i);
		sum2 += xx * yy * rxyz5i;
		sum3 += xx * zz * rxyz5i;
		sum4 += rxyz3i + (yy2 * rxyz5i);
		sum5 += yy * zz * rxyz5i;
		sum6 += rxyz3i + (zz2 * rxyz5i);
}
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int iatom = 0; iatom < numberOfAtoms; ++iatom){
    	if(rxyz_vec[iatom] < dmaxx)  dmaxx = rxyz_vec[iatom];
    }
    *dmax = dmaxx;
    delete[] rxyz_vec;
    *pot = e00 - (dipol * ((rx * rx) + (ry * ry) + (rz * rz)));
    *dpotx = de00x -(dipol * ((2.0 * rx * sum1) + (2.0 * ry * sum2) + (2.0 * rz * sum3)));
    *dpoty = de00y -(dipol * ((2.0 * rx * sum2) + (2.0 * ry * sum4) + (2.0 * rz * sum5)));
    *dpotz = de00z -(dipol * ((2.0 * rx * sum3) + (2.0 * ry * sum5) + (2.0 * rz * sum6)));
}

void dljpotHe(const double x,const double y,const double z, double *pot, double *dmax, double *dpotx, double *dpoty, double *dpotz){

//  Subroutine to calculate L-J + ion-dipole potential.

	double rxyz[8] __attribute__((aligned(16)));
    double rx, ry, rz, e00, de00, de00x, de00y, de00z;
    double sum1, sum2, sum3, sum4, sum5, sum6, dmaxx, eox4;
    double xx, xx2, yy, yy2, zz, zz2, rxyz3i, rxyz5i;
    double *rxyz_vec __attribute__((aligned(16)));

    rx = 0.0;
    ry = 0.0;
    rz = 0.0;
    e00 = 0.0;
    de00x = 0.0;
    de00y = 0.0;
    de00z = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    sum6 = 0.0;
    dmaxx = 2.0 * romax;
    rxyz_vec = new double[numberOfAtoms];
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int iatom = 0; iatom < numberOfAtoms; iatom++){
    	xx = x - molecule->fx[iatom];
		xx2 = xx * xx;
		yy = y - molecule->fy[iatom];
		yy2 = yy * yy;
		zz = z - molecule->fz[iatom];
		zz2 = zz * zz;
        rxyz[0] = xx2 + yy2 + zz2;
        rxyz[1] = sqrt(rxyz[0]);
        rxyz_vec[iatom] = rxyz[1];
        rxyz[2] = rxyz[0] * rxyz[1];
        rxyz[3] = rxyz[2] * rxyz[0];
        rxyz[4] = rxyz[3] * rxyz[1];
        rxyz[5] = rxyz[3] * rxyz[2];
        rxyz[6] = rxyz[4] * rxyz[4];
		rxyz[7] = rxyz[6] * rxyz[0];

		// LJ potential
		eox4 = molecule->eox4[iatom];
		e00 += (eox4 * ((molecule->ro12lj[iatom] / rxyz[6]) - (molecule->ro6lj[iatom] / rxyz[4])));

        // LJ derivative
        de00 = eox4 * ((molecule->dro6[iatom] / rxyz[5]) - (molecule->dro12[iatom] / rxyz[7]));
        de00x += de00 * xx;
        de00y += de00 * yy;
        de00z += de00 * zz;

        // ion-induced dipole potential
		rxyz3i = molecule->charge[iatom] / rxyz[2];
		rxyz5i = -3.0 * molecule->charge[iatom] / rxyz[3];
		rx += xx * rxyz3i;
		ry += yy * rxyz3i;
		rz += zz * rxyz3i;
		//  ion-induced dipole derivative
		sum1 += rxyz3i + (xx2 * rxyz5i);
		sum2 += xx * yy * rxyz5i;
		sum3 += xx * zz * rxyz5i;
		sum4 += rxyz3i + (yy2 * rxyz5i);
		sum5 += yy * zz * rxyz5i;
		sum6 += rxyz3i + (zz2 * rxyz5i);


}
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int iatom = 0; iatom < numberOfAtoms; ++iatom){
        	if(rxyz_vec[iatom] < dmaxx)  dmaxx = rxyz_vec[iatom];
    }
    *dmax = dmaxx;
    delete[] rxyz_vec;
    *pot = e00 - (dipol * ((rx * rx) + (ry * ry) + (rz * rz)));
    *dpotx = de00x -(dipol * ((2.0 * rx * sum1) + (2.0 * ry * sum2) + (2.0 * rz * sum3)));
    *dpoty = de00y -(dipol * ((2.0 * rx * sum2) + (2.0 * ry * sum4) + (2.0 * rz * sum5)));
    *dpotz = de00z -(dipol * ((2.0 * rx * sum3) + (2.0 * ry * sum5) + (2.0 * rz * sum6)));

}
