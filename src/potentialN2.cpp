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

double dljpotN2(const double x, const double y, const double z, double *fx, double *fy, double *fz){

	//     Subroutine to calculate L-J + ion-dipole potential.

		double rx, ry, rz, e00, qpol, pot;
		double xc, yc, zc, dpolx, dpoly, dpolz;
		double bond, Ptfn, xkT, pc, pc_center;
		double xx_center, xx, xx_center2, xx2;
		double yy_center, yy, yy_center2, yy2;
		double zz_center, zz, zz_center2, zz2;
		double dipolxx, dipolzz, pot_min;
		double rxyz_center2, rxyz2, rxyz_center, rxyz, rxyz3;
		double rxyz5, rxyz6, rxyz8, rxyz12, rxyz14, rxyz_center3;
		double rxyz_center5, const_k, rxyz3i, rxyz5i, temp_pot;
		double tpot, weight;
	    double pottry[2][3], pot_mol[3];
	    double i_pot, i_dpotx, i_dpoty, i_dpotz;
	    int isamp;
	//     nitrogen : five charge model (Allen Tildesley, 14 page)
	//     Every data from B3LYP//aug-cc-pVDZ
	      bond = 1.0976e-10;
	      Ptfn = 0.0;
	      xkT = 500.0 * xk;
	      pc = -0.4825;
	      pc_center = -(pc);
	      dipolzz = 1.710e-30 / (2.0 * 4.0 * pi * xeo);
	      dipolzz *= xe * xe;
	      //dipolxx = 1.710e-30 / (2.0 * 4.0 * pi * xeo);
	      dipolxx = dipolzz;
	      pot_min = 1.0e+8;

	      for(isamp = 0; isamp < 3; isamp++){
	    	  for(int ibatom = 0; ibatom < 2; ibatom++){
				  rx = 0.0;
				  ry = 0.0;
				  rz = 0.0;
				  e00 = 0.0;
				  qpol = 0.0;
				  #if defined(__INTEL_COMPILER)
				  #pragma vector aligned
				  #endif
				  for(int iatom = 0; iatom < numberOfAtoms; iatom++){

					  if(isamp == 0){
						  xc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
						  yc = 0.0;
						  zc = 0.0;
						  dpolx = dipolzz;
						  dpoly = dipolxx;
						  dpolz = dipolxx;
					  }
					  if(isamp == 1){
						  xc = 0.0;
						  yc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
						  zc = 0.0;
						  dpolx = dipolxx;
						  dpoly = dipolzz;
						  dpolz = dipolxx;
					  }
					  if(isamp == 2){
						  xc = 0.0;
						  yc = 0.0;
						  zc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
						  dpolx = dipolxx;
						  dpoly = dipolxx;
						  dpolz = dipolzz;
					  }

					  xx_center = x - fx[iatom];
					  xx = xx_center + xc;
					  xx_center2 = xx_center * xx_center;
					  xx2 = xx * xx;

					  yy_center = y - fy[iatom];
					  yy = yy_center + yc;
					  yy_center2 = yy_center * yy_center;
					  yy2 = yy * yy;

					  zz_center = z - fz[iatom];
					  zz = zz_center + zc;
					  zz_center2 = zz_center * zz_center;
					  zz2 = zz * zz;

					  rxyz_center2 = xx_center2 + yy_center2 + zz_center2;
					  rxyz2 = xx2 + yy2 + zz2;

					  rxyz_center = sqrt(rxyz_center2);
					  rxyz = sqrt(rxyz2);
					  rxyz3 = rxyz2 * rxyz;
					  rxyz5 = rxyz3 * rxyz2;
					  rxyz6 = rxyz5 * rxyz;
					  rxyz12 = rxyz6 * rxyz6;

					  rxyz_center3 = rxyz_center2 * rxyz_center;

	//			      LJ potential
					  e00 += (molecule->eox4[iatom] * ((molecule->ro12lj[iatom]/rxyz12) - (molecule->ro6lj[iatom]/rxyz6)));

	//			      ion-induced dipole potential

					  rxyz3i = molecule->charge[iatom] / rxyz_center3;
					  rx += xx_center * rxyz3i;
					  ry += yy_center * rxyz3i;
					  rz += zz_center * rxyz3i;

	//			     ion-partial charge coulomb potential(quadrupole)
					  const_k = molecule->charge[iatom] * (xe * xe) / (4.0 * pi * xeo);
					  qpol += pc_center * const_k/rxyz_center;
					  qpol += pc * const_k / rxyz;
				  }
			 pottry[ibatom][isamp] = e00 - 0.5 * (((dpolx * rx * rx) + (dpoly * ry * ry) + (dpolz * rz * rz))) + qpol;

				  }
			  pot_mol[isamp] = pottry[0][isamp] + pottry[1][isamp];
			  if(pot_min >= pot_mol[isamp]) pot_min = pot_mol[isamp];
	      }

	      for(isamp = 0; isamp < 3; isamp++){
	    	  temp_pot = pot_mol[isamp] - pot_min;
	    	  Ptfn += exp(-temp_pot / xkT);
		}

	      pot = 0.0;
	      for(isamp = 0; isamp < 3; isamp++){
	    	  temp_pot = pot_mol[isamp] - pot_min;
	    	  weight = exp(-temp_pot/xkT)/Ptfn;
			  pot += weight * pot_mol[isamp];
	      }

	      return pot;
}

double dljpotN2(const double x, const double y, const double z){

//     Subroutine to calculate L-J + ion-dipole potential.

	double rx, ry, rz, e00, qpol, pot;
	double xc, yc, zc, dpolx, dpoly, dpolz;
	double bond, Ptfn, xkT, pc, pc_center;
	double xx_center, xx, xx_center2, xx2;
	double yy_center, yy, yy_center2, yy2;
	double zz_center, zz, zz_center2, zz2;
	double dipolxx, dipolzz, pot_min;
	double rxyz_center2, rxyz2, rxyz_center, rxyz, rxyz3;
	double rxyz5, rxyz6, rxyz8, rxyz12, rxyz14, rxyz_center3;
	double rxyz_center5, const_k, rxyz3i, rxyz5i, temp_pot;
	double tpot, weight;
    double pottry[2][3], pot_mol[3];
    double i_pot, i_dpotx, i_dpoty, i_dpotz;
    int isamp;

//     nitrogen : five charge model (Allen Tildesley, 14 page)
//     Every data from B3LYP//aug-cc-pVDZ
      bond = 1.0976e-10;
      Ptfn = 0.0;
      xkT = 500.0 * xk;
      pc = -0.4825;
      pc_center = -(pc);
      dipolzz = 1.710e-30 / (2.0 * 4.0 * pi * xeo);
      dipolzz *= xe * xe;
      //dipolxx = 1.710e-30 / (2.0 * 4.0 * pi * xeo);
      dipolxx = dipolzz;
      pot_min = 1.0e+8;

      for(isamp = 0; isamp < 3; isamp++){
    	  for(int ibatom = 0; ibatom < 2; ibatom++){
			  rx = 0.0;
			  ry = 0.0;
			  rz = 0.0;
			  e00 = 0.0;
			  qpol = 0.0;

			  #if defined(__INTEL_COMPILER)
              #pragma vector aligned
			  #endif
			  for(int iatom = 0; iatom < numberOfAtoms; iatom++){

				  if(isamp == 0){
					  xc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  yc = 0.0;
					  zc = 0.0;
					  dpolx = dipolzz;
					  dpoly = dipolxx;
					  dpolz = dipolxx;
				  }
				  if(isamp == 1){
					  xc = 0.0;
					  yc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  zc = 0.0;
					  dpolx = dipolxx;
					  dpoly = dipolzz;
					  dpolz = dipolxx;
				  }
				  if(isamp == 2){
					  xc = 0.0;
					  yc = 0.0;
					  zc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  dpolx = dipolxx;
					  dpoly = dipolxx;
					  dpolz = dipolzz;
				  }

				  xx_center = x - molecule->fx[iatom];
				  xx = xx_center + xc;
				  xx_center2 = xx_center * xx_center;
				  xx2 = xx * xx;

				  yy_center = y - molecule->fy[iatom];
				  yy = yy_center + yc;
				  yy_center2 = yy_center * yy_center;
				  yy2 = yy * yy;

				  zz_center = z - molecule->fz[iatom];
				  zz = zz_center + zc;
				  zz_center2 = zz_center * zz_center;
				  zz2 = zz * zz;

				  rxyz_center2 = xx_center2 + yy_center2 + zz_center2;
				  rxyz2 = xx2 + yy2 + zz2;

				  rxyz_center = sqrt(rxyz_center2);
				  rxyz = sqrt(rxyz2);
				  rxyz3 = rxyz2 * rxyz;
				  rxyz5 = rxyz3 * rxyz2;
				  rxyz6 = rxyz5 * rxyz;
				  rxyz12 = rxyz6 * rxyz6;

				  rxyz_center3 = rxyz_center2 * rxyz_center;

//			      LJ potential
				  e00 += (molecule->eox4[iatom] * ((molecule->ro12lj[iatom]/rxyz12) - (molecule->ro6lj[iatom]/rxyz6)));

//			      ion-induced dipole potential

				  rxyz3i = molecule->charge[iatom] / rxyz_center3;
				  rx += xx_center * rxyz3i;
				  ry += yy_center * rxyz3i;
				  rz += zz_center * rxyz3i;

//			     ion-partial charge coulomb potential(quadrupole)
				  const_k = molecule->charge[iatom] * (xe * xe) / (4.0 * pi * xeo);
				  qpol += pc_center * const_k/rxyz_center;
				  qpol += pc * const_k / rxyz;
			  }
		 pottry[ibatom][isamp] = e00 - 0.5 * (((dpolx * rx * rx) + (dpoly * ry * ry) + (dpolz * rz * rz))) + qpol;

			  }
		  pot_mol[isamp] = pottry[0][isamp] + pottry[1][isamp];
		  if(pot_min >= pot_mol[isamp]) pot_min = pot_mol[isamp];
      }

      for(isamp = 0; isamp < 3; isamp++){
    	  temp_pot = pot_mol[isamp] - pot_min;
    	  Ptfn += exp(-temp_pot / xkT);
	}

      pot = 0.0;
      for(isamp = 0; isamp < 3; isamp++){
    	  temp_pot = pot_mol[isamp] - pot_min;
    	  weight = exp(-temp_pot/xkT)/Ptfn;
		  pot += weight * pot_mol[isamp];
      }

      return pot;
}

void dljpotN2(const double x, const double y, const double z, double *pot, double *dmax, double *dpotx, double *dpoty, double *dpotz){

//     Subroutine to calculate L-J + ion-dipole potential.

	double rx, ry, rz, e00, de00, de00x, de00y, de00z;
	double sum1, sum2, sum3, sum4, sum5, sum6;
	double qpol, dqpolx, dqpoly, dqpolz;
	double xc, yc, zc, dpolx, dpoly, dpolz;
	double bond, Ptfn, xkT, pc, pc_center;
	double xx_center, xx, xx_center2, xx2;
	double yy_center, yy, yy_center2, yy2;
	double zz_center, zz, zz_center2, zz2;
	double dipolxx, dipolzz, pot_min;
	double rxyz_center2, rxyz2, rxyz_center, rxyz, rxyz3;
	double rxyz5, rxyz6, rxyz8, rxyz12, rxyz14, rxyz_center3;
	double rxyz_center5, const_k, rxyz3i, rxyz5i, temp_pot;
	double tpot, weight;
    double pottry[2][3], dpotxtry[2][3], dpotytry[2][3], dpotztry[2][3] ;
    double pot_mol[3], dpotx_mol[3], dpoty_mol[3], dpotz_mol[3];
    double i_pot, i_dpotx, i_dpoty, i_dpotz;
    int isamp;
//     nitrogen : five charge model (Allen Tildesley, 14 page)
//     Every data from B3LYP//aug-cc-pVDZ
      *dmax = 2.0 * romax;
      bond = 1.0976e-10;
      Ptfn = 0.0;
      xkT = 500.0 * xk;
      pc = -0.4825;
      pc_center = -(pc);
      dipolzz = 1.710e-30 / (2.0 * 4.0 * pi * xeo);
      dipolzz *= xe * xe;
      dipolxx = dipolzz;
      pot_min = 1.0e+8;

      for(isamp = 0; isamp < 3; isamp++){
    	  for(int ibatom = 0; ibatom < 2; ibatom++){
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
			  qpol = 0.0;
			  dqpolx = 0.0;
			  dqpoly = 0.0;
			  dqpolz = 0.0;

			  #if defined(__INTEL_COMPILER)
    		  #pragma vector aligned
			  #endif
			  for(int iatom = 0; iatom < numberOfAtoms; iatom++){

				  if(isamp == 0){
					  xc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  yc = 0.0;
					  zc = 0.0;
					  dpolx = dipolzz;
					  dpoly = dipolxx;
					  dpolz = dipolxx;
				  }
				  if(isamp == 1){
					  xc = 0.0;
					  yc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  zc = 0.0;
					  dpolx = dipolxx;
					  dpoly = dipolzz;
					  dpolz = dipolxx;
				  }
				  if(isamp == 2){
					  xc = 0.0;
					  yc = 0.0;
					  zc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  dpolx = dipolxx;
					  dpoly = dipolxx;
					  dpolz = dipolzz;
				  }

				  xx_center = x - molecule->fx[iatom];
				  xx = xx_center + xc;
				  xx_center2 = xx_center * xx_center;
				  xx2 = xx * xx;

				  yy_center = y - molecule->fy[iatom];
				  yy = yy_center + yc;
				  yy_center2 = yy_center * yy_center;
				  yy2 = yy * yy;

				  zz_center = z - molecule->fz[iatom];
				  zz = zz_center + zc;
				  zz_center2 = zz_center * zz_center;
				  zz2 = zz * zz;

				  rxyz_center2 = xx_center2 + yy_center2 + zz_center2;
				  rxyz2 = xx2 + yy2 + zz2;

				  rxyz_center = sqrt(rxyz_center2);
				  rxyz = sqrt(rxyz2);
				  if(rxyz < *dmax) *dmax = rxyz;
				  rxyz3 = rxyz2 * rxyz;
				  rxyz5 = rxyz3 * rxyz2;
				  rxyz6 = rxyz5 * rxyz;
				  rxyz8 = rxyz5 * rxyz3;
				  rxyz12 = rxyz6 * rxyz6;
				  rxyz14 = rxyz12 * rxyz2;
				  rxyz_center3 = rxyz_center2 * rxyz_center;
				  rxyz_center5 = rxyz_center3 * rxyz_center2;
//			      LJ potential
				  e00 += (molecule->eox4[iatom] * ((molecule->ro12lj[iatom]/rxyz12) - (molecule->ro6lj[iatom]/rxyz6)));
//			      LJ derivative
				  de00 = molecule->eox4[iatom] * ((molecule->dro6[iatom] / rxyz8) - (molecule->dro12[iatom]/rxyz14));
				  de00x += de00 * xx;
				  de00y += de00 * yy;
				  de00z += de00 * zz;
//			      ion-induced dipole potential

				  rxyz3i = molecule->charge[iatom] / rxyz_center3;
				  rxyz5i = -3.0 * molecule->charge[iatom] / rxyz_center5;
				  rx += xx_center * rxyz3i;
				  ry += yy_center * rxyz3i;
				  rz += zz_center * rxyz3i;
//			     ion-induced dipole derivative
				  sum1 += rxyz3i + (xx_center2 * rxyz5i);
				  sum2 += xx_center * yy_center * rxyz5i;
				  sum3 += xx_center * zz_center * rxyz5i;
				  sum4 += rxyz3i + (yy_center2 * rxyz5i);
				  sum5 += yy_center * zz_center * rxyz5i;
				  sum6 += rxyz3i + (zz_center2 * rxyz5i);
//			     ion-partial charge coulomb potential(quadrupole)
				  const_k = molecule->charge[iatom] * (xe * xe) / (4.0 * pi * xeo);
				  qpol += (pc_center * const_k/rxyz_center);
				  qpol += (pc * const_k / rxyz);
//			     ion-partial charge coulomb derivative(quadrupole)
				  dqpolx -= (pc_center * const_k / rxyz_center3) * (xx_center);
				  dqpoly -= (pc_center * const_k / rxyz_center3) * (yy_center);
				  dqpolz -= (pc_center * const_k / rxyz_center3) * (zz_center);
				  dqpolx -= (pc * const_k / rxyz3) * xx;
				  dqpoly -= (pc * const_k / rxyz3) * yy;
				  dqpolz -= (pc * const_k / rxyz3) * zz;

			  }
		 pottry[ibatom][isamp] = e00 - 0.5 * (((dpolx * rx * rx) + (dpoly * ry * ry) + (dpolz * rz * rz))) + qpol;
		 dpotxtry[ibatom][isamp] = de00x - 0.5 * ((dpolx  *2.0 * rx * sum1)
				 + (dpoly * 2.0 * ry *sum2)+(dpolz*2.0*rz*sum3))+dqpolx;
		 dpotytry[ibatom][isamp] = de00y - 0.5 * ((dpolx * 2.0  *rx * sum2)
				 +(dpoly * 2.0 * ry * sum4) + (dpolz * 2.0 * rz * sum5)) + dqpoly;
		 dpotztry[ibatom][isamp] = de00z - 0.5 * ((dpolx * 2.0 * rx * sum3)
				 +(dpoly * 2.0 * ry * sum5) + (dpolz * 2.0 * rz * sum6)) + dqpolz;
			  }
		  pot_mol[isamp] = pottry[0][isamp] + pottry[1][isamp];
		  //tpot = pot_mol[isamp];
		  if(pot_min >= pot_mol[isamp]) pot_min = pot_mol[isamp];
		  dpotx_mol[isamp] = dpotxtry[0][isamp] + dpotxtry[1][isamp];
		  dpoty_mol[isamp] = dpotytry[0][isamp] + dpotytry[1][isamp];
		  dpotz_mol[isamp] = dpotztry[0][isamp] + dpotztry[1][isamp];
      }

      for(isamp = 0; isamp < 3; isamp++){
    	  temp_pot = pot_mol[isamp] - pot_min;
    	  Ptfn += exp(-temp_pot / xkT);
	}

      i_pot = 0.0;
      i_dpotx = 0.0;
      i_dpoty = 0.0;
      i_dpotz = 0.0;
      for(isamp = 0; isamp < 3; isamp++){
    	  temp_pot = pot_mol[isamp] - pot_min;
    	  weight = exp(-temp_pot/xkT)/Ptfn;
		  i_pot += weight * pot_mol[isamp];
    	  i_dpotx += weight * dpotx_mol[isamp];
    	  i_dpoty += weight * dpoty_mol[isamp];
    	  i_dpotz += weight * dpotz_mol[isamp];
      }

      *pot = i_pot;
      *dpotx = i_dpotx;
      *dpoty = i_dpoty;
      *dpotz = i_dpotz;
}

void dljpotN2(const double x, const double y, const double z, double *pot,double *dmax, double *dpotx, double *dpoty, double *dpotz, double *fx, double *fy, double *fz){

//     Subroutine to calculate L-J + ion-dipole potential.

	double rx, ry, rz, e00, de00, de00x, de00y, de00z;
	double sum1, sum2, sum3, sum4, sum5, sum6;
	double qpol, dqpolx, dqpoly, dqpolz;
	double xc, yc, zc, dpolx, dpoly, dpolz;
	double bond, Ptfn, xkT, pc, pc_center;
	double xx_center, xx, xx_center2, xx2;
	double yy_center, yy, yy_center2, yy2;
	double zz_center, zz, zz_center2, zz2;
	double dipolxx, dipolzz, pot_min;
	double rxyz_center2, rxyz2, rxyz_center, rxyz, rxyz3;
	double rxyz5, rxyz6, rxyz8, rxyz12, rxyz14, rxyz_center3;
	double rxyz_center5, const_k, rxyz3i, rxyz5i, temp_pot;
	double tpot, weight;
    double pottry[2][3], dpotxtry[2][3], dpotytry[2][3], dpotztry[2][3] ;
    double pot_mol[3], dpotx_mol[3], dpoty_mol[3], dpotz_mol[3];
    double i_pot, i_dpotx, i_dpoty, i_dpotz;
    int isamp;
//     nitrogen : five charge model (Allen Tildesley, 14 page)
//     Every data from B3LYP//aug-cc-pVDZ
      *dmax = 2.0 * romax;
      bond = 1.0976e-10;
      Ptfn = 0.0;
      xkT = 500.0 * xk;
      pc = -0.4825;
      pc_center = -(pc);
      dipolzz = 1.710e-30 / (2.0 * 4.0 * pi * xeo);
      dipolzz *= xe * xe;
      //dipolxx = 1.710e-30 / (2.0 * 4.0 * pi * xeo);
      dipolxx = dipolzz;
      pot_min = 1.0e+8;

      for(isamp = 0; isamp < 3; isamp++){
    	  for(int ibatom = 0; ibatom < 2; ibatom++){
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
			  qpol = 0.0;
			  dqpolx = 0.0;
			  dqpoly = 0.0;
			  dqpolz = 0.0;

			  #if defined(__INTEL_COMPILER)
    		  #pragma vector aligned
			  #endif
			  for(int iatom = 0; iatom < numberOfAtoms; iatom++){

				  if(isamp == 0){
					  xc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  yc = 0.0;
					  zc = 0.0;
					  dpolx = dipolzz;
					  dpoly = dipolxx;
					  dpolz = dipolxx;
				  }
				  if(isamp == 1){
					  xc = 0.0;
					  yc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  zc = 0.0;
					  dpolx = dipolxx;
					  dpoly = dipolzz;
					  dpolz = dipolxx;
				  }
				  if(isamp == 2){
					  xc = 0.0;
					  yc = 0.0;
					  zc = (bond / 2.0) * (2.0 * (ibatom + 1) - 3.0);
					  dpolx = dipolxx;
					  dpoly = dipolxx;
					  dpolz = dipolzz;
				  }

				  xx_center = x - fx[iatom];
				  xx = xx_center + xc;
				  xx_center2 = xx_center * xx_center;
				  xx2 = xx * xx;

				  yy_center = y - fy[iatom];
				  yy = yy_center + yc;
				  yy_center2 = yy_center * yy_center;
				  yy2 = yy * yy;

				  zz_center = z - fz[iatom];
				  zz = zz_center + zc;
				  zz_center2 = zz_center * zz_center;
				  zz2 = zz * zz;

				  rxyz_center2 = xx_center2 + yy_center2 + zz_center2;
				  rxyz2 = xx2 + yy2 + zz2;

				  rxyz_center = sqrt(rxyz_center2);
				  rxyz = sqrt(rxyz2);
				  if(rxyz < *dmax) *dmax = rxyz;
				  rxyz3 = rxyz2 * rxyz;
				  rxyz5 = rxyz3 * rxyz2;
				  rxyz6 = rxyz5 * rxyz;
				  rxyz8 = rxyz5 * rxyz3;
				  rxyz12 = rxyz6 * rxyz6;
				  rxyz14 = rxyz12 * rxyz2;
				  rxyz_center3 = rxyz_center2 * rxyz_center;
				  rxyz_center5 = rxyz_center3 * rxyz_center2;
//			      LJ potential
				  e00 += (molecule->eox4[iatom] * ((molecule->ro12lj[iatom]/rxyz12) - (molecule->ro6lj[iatom]/rxyz6)));
//			      LJ derivative
				  de00 = molecule->eox4[iatom] * ((molecule->dro6[iatom] / rxyz8) - (molecule->dro12[iatom]/rxyz14));
				  de00x += de00 * xx;
				  de00y += de00 * yy;
				  de00z += de00 * zz;
//			      ion-induced dipole potential

				  rxyz3i = molecule->charge[iatom] / rxyz_center3;
				  rxyz5i = -3.0 * molecule->charge[iatom] / rxyz_center5;
				  rx += xx_center * rxyz3i;
				  ry += yy_center * rxyz3i;
				  rz += zz_center * rxyz3i;
//			     ion-induced dipole derivative
				  sum1 += rxyz3i + (xx_center2 * rxyz5i);
				  sum2 += xx_center * yy_center * rxyz5i;
				  sum3 += xx_center * zz_center * rxyz5i;
				  sum4 += rxyz3i + (yy_center2 * rxyz5i);
				  sum5 += yy_center * zz_center * rxyz5i;
				  sum6 += rxyz3i + (zz_center2 * rxyz5i);
//			     ion-partial charge coulomb potential(quadrupole)
				  const_k = molecule->charge[iatom] * (xe * xe) / (4.0 * pi * xeo);
				  qpol += (pc_center * const_k/rxyz_center);
				  qpol += (pc * const_k / rxyz);
//			     ion-partial charge coulomb derivative(quadrupole)
				  dqpolx -= (pc_center * const_k / rxyz_center3) * (xx_center);
				  dqpoly -= (pc_center * const_k / rxyz_center3) * (yy_center);
				  dqpolz -= (pc_center * const_k / rxyz_center3) * (zz_center);
				  dqpolx -= (pc * const_k / rxyz3) * xx;
				  dqpoly -= (pc * const_k / rxyz3) * yy;
				  dqpolz -= (pc * const_k / rxyz3) * zz;

			  }
		 pottry[ibatom][isamp] = e00 - 0.5 * (((dpolx * rx * rx) + (dpoly * ry * ry) + (dpolz * rz * rz))) + qpol;
		 dpotxtry[ibatom][isamp] = de00x - 0.5 * ((dpolx  *2.0 * rx * sum1)
				 + (dpoly * 2.0 * ry *sum2)+(dpolz*2.0*rz*sum3))+dqpolx;
		 dpotytry[ibatom][isamp] = de00y - 0.5 * ((dpolx * 2.0  *rx * sum2)
				 +(dpoly * 2.0 * ry * sum4) + (dpolz * 2.0 * rz * sum5)) + dqpoly;
		 dpotztry[ibatom][isamp] = de00z - 0.5 * ((dpolx * 2.0 * rx * sum3)
				 +(dpoly * 2.0 * ry * sum5) + (dpolz * 2.0 * rz * sum6)) + dqpolz;
			  }
		  pot_mol[isamp] = pottry[0][isamp] + pottry[1][isamp];
		  //tpot = pot_mol[isamp];
		  if(pot_min >= pot_mol[isamp]) pot_min = pot_mol[isamp];
		  dpotx_mol[isamp] = dpotxtry[0][isamp] + dpotxtry[1][isamp];
		  dpoty_mol[isamp] = dpotytry[0][isamp] + dpotytry[1][isamp];
		  dpotz_mol[isamp] = dpotztry[0][isamp] + dpotztry[1][isamp];
      }

      for(isamp = 0; isamp < 3; isamp++){
    	  temp_pot = pot_mol[isamp] - pot_min;
    	  Ptfn += exp(-temp_pot / xkT);
	}

      i_pot = 0.0;
      i_dpotx = 0.0;
      i_dpoty = 0.0;
      i_dpotz = 0.0;
      for(isamp = 0; isamp < 3; isamp++){
    	  temp_pot = pot_mol[isamp] - pot_min;
    	  weight = exp(-temp_pot/xkT)/Ptfn;
		  i_pot += weight * pot_mol[isamp];
    	  i_dpotx += weight * dpotx_mol[isamp];
    	  i_dpoty += weight * dpoty_mol[isamp];
    	  i_dpotz += weight * dpotz_mol[isamp];
      }

      *pot = i_pot;
      *dpotx = i_dpotx;
      *dpoty = i_dpoty;
      *dpotz = i_dpotz;
}

