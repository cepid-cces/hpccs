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

#include <headers/globals.hpp>
#include <headers/constants.hpp>
#include <headers/rotate.hpp>
#include <headers/gsang.hpp>
#include <headers/mt.h>
#include <headers/potentialHe.hpp>
#include <headers/potentialN2.hpp>
#include <iostream>
#include <omp.h>
#include <cmath>
#if defined(__INTEL_COMPILER)
#include <aligned_new>
#else
#include <new>
#endif
using namespace std;

int mobil2(){


    double *q1st, *q2st, *pgst, *wgst, *b2max;
    double *om11st, *om12st, *om13st, *om22st;
    double pot, ayst,term, w, r, rnb,ddd;
    double mom11st, mom12st, mom13st, mom22st;
	double rxy, rzy, rmax, rmaxx,sdom11st, best, cest;
	double dbst2, b, bst2, delta, gst, ang;
    double tst, tst3, u2, v, sum, sum1, sum2,temp, temp1;
    double f, dgst, dbst22, gst2, gstt, hold, hold1, hold2;
    int i,ihold, irn, ig, im, ic, ibst;
    unsigned int iatom, total;
    double cosx[2500] __attribute__((aligned(16)));
    double top,dt1,dt2,e0, qerqw;
    double *v_vec __attribute__((aligned(16)));
    double *temp_vec __attribute__((aligned(16)));
    double *temp1_vec __attribute__((aligned(16)));
    double *temp_vec2 __attribute__((aligned(16)));
    double *temp1_vec2 __attribute__((aligned(16)));
    double *b_vec __attribute__((aligned(16)));
    double *b2max_vec __attribute__((aligned(16)));
    int itn_inp = itn * inp, id;

    trajlost = 0;
    total =  itn * inp * imp;
	// Initialize a Mersenne Twister
    q1st = new double [inp]();
    q2st = new double [inp]();
    pgst = new double [inp]();
    wgst = new double [inp]();
    b2max = new double [inp]();

    om11st = new double [itn]();
    om12st = new double [itn]();
    om13st = new double [itn]();
    om22st = new double [itn]();


    ihold = 0;
    rmax = 0.0;
    rmaxx = 0.0;


//     determine maximum extent and orientate along x axis
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
	for (iatom = 0; iatom < numberOfAtoms; iatom++){
        r = sqrt((molecule->ox[iatom] * molecule->ox[iatom])
        		+ (molecule->oy[iatom] * molecule->oy[iatom])
        		+ (molecule->oz[iatom] * molecule->oz[iatom]));

        if(r > rmax){
            rmax = r;
            ihold = iatom;
        }
	}

    rzy = sqrt((molecule->oy[ihold] * molecule->oy[ihold])
    		+ (molecule->oz[ihold] * molecule->oz[ihold]));
    phi = acos(molecule->oz[ihold] / rzy);
    phi += (pi / 2.0);
    if(molecule->oy[ihold] < 0.0) phi = (2.0 * pi) - phi;
    phi = (2.0 * pi) - phi;
    theta = 0.0;
    agamma = 0.0;
    rotate();


    rxy = sqrt((molecule->fx[ihold] * molecule->fx[ihold])
		+ (molecule->fy[ihold] * molecule->fy[ihold]));
    agamma = acos(molecule->fx[ihold] / rxy);
    if(molecule->fy[ihold] < 0.0) agamma = (2.0 * pi) - agamma;
    agamma = (2.0 * pi) - agamma;
    rotate();

//*************************************************************************************************************************************
    hold = molecule->fx[ihold] / rmax;
    if(hold < 0.9999999999 || hold > 1.0000000001
    || molecule->fy[ihold] > 1.0e-20 || molecule->fz[ihold] > 1.0e-20
    || molecule->fy[ihold] < -1.0e-20 || molecule->fz[ihold] < -1.0e-20){
    	cout << "Problem orientating along x axis " << "\n";

    #pragma vector aligned
    for (iatom = 0; iatom < numberOfAtoms; iatom++){
        hold = sqrt((molecule->fx[iatom] * molecule->fx[iatom])
        		+ (molecule->fy[iatom] * molecule->fy[iatom])
        		+ (molecule->fz[iatom] * molecule->fz[iatom]));
        cout << iatom << "\t" << molecule->fx[iatom] << "\t" << molecule->fy[iatom] << "\t" <<  molecule->fz[iatom] << "\t" << hold << "\n";
    }
    	return 0;
	}

//   Determine rmax, emax, and r00 along x, y, and z directions
    irn = 1000;
    ddd = (rmax + romax) / float(irn);

//   Set-up integration over gst

    tst = xk * temperature / eo;
    tst3 = pow(tst,3);
    dgst = 5.0e-7 * 6.0 * sqrt(tst);
    gst = dgst;
    sum = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;

    for (i = 1; i <= inp; ++i){
        sum1 += sqrt(float(i));
    }

     for (i = 0; i < inp; ++i){
        hold1 = sqrt(float(i + 1));
        hold2 = sqrt(float(i));
        sum2 += hold2;
        wgst[i] = hold1 / sum1;
        gstt = tst3 * (sum2 + (hold1 / 2.0)) /sum1;
        do{
        	sum += exp(-gst * gst/tst) * pow(gst,5) * dgst;
        	gst += dgst;
        	if(sum > gstt) pgst[i] = gst - (dgst / 2.0);
        }while(sum < gstt);
     }

//   determine b2max
    dbst2 = 1.0;
    dbst22 = dbst2 / 10.0;


	if (gas == 2){
	for (ig = inp -1; ig >= 0; --ig){
		gst2 = pgst[ig] * pgst[ig];
		v = sqrt((gst2 * eo) /(0.5 * mu));
		ibst = int(rmaxx / ro) - 6;
		if(ig < inp - 1) ibst = int(b2max[ig + 1] / dbst2) - 6;
		if(ibst < 0) ibst = 0;

		while(true){
			while(true){
				bst2 = dbst2 * float(ibst);
				b = ro * sqrt(bst2);
				ang = gsangN2(v,b);
				cosx[ibst] = 1.0 - cos(ang);
				if(ibst < 5){
				   ibst++;
				}else{break;}
			}

			if(cosx[ibst] < cmin && cosx[ibst-1] < cmin &&
				cosx[ibst-2] < cmin && cosx[ibst-3] < cmin &&
				cosx[ibst-4] < cmin){
					b2max[ig] = float(ibst-5) * dbst2;
					do {
						b2max[ig] += dbst22;
						b = ro * sqrt(b2max[ig]);
						ang = gsangN2(v,b);
					}while (1.0 - cos(ang) > cmin);

					break;
				}else{
					ibst++;
					if(ibst > 2500) {
						cout << "ibst greater than 500" << "\n";
						return 1;
					}
				}
			}
		}
	}else{
	   for (ig = inp -1; ig >= 0; --ig){
				gst2 = pgst[ig] * pgst[ig];
				v = sqrt((gst2 * eo) /(0.5 * mu));
				ibst = int(rmaxx / ro) - 6;
				if(ig < inp - 1) ibst = int(b2max[ig + 1] / dbst2) - 6;
				if(ibst < 0) ibst = 0;

				while(true){
					while(true){
						bst2 = dbst2 * float(ibst);
						b = ro * sqrt(bst2);
						ang = gsangHe(v,b);
						cosx[ibst] = 1.0 - cos(ang);
						if(ibst < 5){
						   ibst++;
						}else{break;}
					}

					if(cosx[ibst] < cmin && cosx[ibst-1] < cmin &&
						cosx[ibst-2] < cmin && cosx[ibst-3] < cmin &&
						cosx[ibst-4] < cmin){
							b2max[ig] = float(ibst-5) * dbst2;
							do {
								b2max[ig] += dbst22;
								b = ro * sqrt(b2max[ig]);
								ang = gsangHe(v,b);
							}while (1.0 - cos(ang) > cmin);

							break;
						}else{
							ibst++;
							if(ibst > 2500) {
								cout << "ibst greater than 500" << "\n";
								return 1;
							}
						}
					}
				}
			}


/*
 *  Calculate Omega(1,1)*, Omega(1,2)*, Omega(1,3)*, and Omega(2,2)*
 *  by integrating Q(1)* or Q(2)* over all orientations, and initial
 *  relative velocities.
 *
 */


	MersenneTwister mt, mt_rotatex, mt_rotatey, mt_rotatez;

	// Initialize a Mersenne Twister

double *rng1,*rng2,*rng3,*rng4;
rng1 = new double[total];
rng2 = new double[total];
rng3 = new double[total];
rng4 = new double[total];

// Initialize a Mersenne Twister
//
//
    mt.init_genrand(5013486);
    mt_rotatex.init_genrand(100000);
    mt_rotatey.init_genrand(200000);
    mt_rotatez.init_genrand(300000);

#if defined(__INTEL_COMPILER)
#pragma vector aligned
#endif
for (int ia = 0; ia < total; ++ia){
	rng1[ia] = mt.genrand_res53();
	rng2[ia] = mt_rotatex.genrand_res53();
	rng3[ia] = mt_rotatey.genrand_res53();
	rng4[ia] = mt_rotatez.genrand_res53();
}

/*******  Monte Carlo Integration ******/
    temp_vec = new double[total]();
    temp1_vec = new double[total]();
    v_vec = new double[total]();
    b_vec = new double[total]();
    b2max_vec = new double[total]();
    temp_vec2 = new double[itn_inp]();
    temp1_vec2 = new double[itn_inp]();

id=0;
	for (ic = 0; ic < itn; ++ic){
		for (ig = 0; ig < inp; ++ig){
			v = sqrt((pgst[ig] * pgst[ig] * eo) / (0.5 * mu));
			#if defined(__INTEL_COMPILER)
    		#pragma vector aligned
			#endif
			for (im = 0; im < imp; ++im){
				v_vec[id] = v;
				b_vec[id] = ro * sqrt(rng1[id] * (b2max[ig]));
				b2max_vec[id] = b2max[ig];
				id++;
			}
		}
	}

	if (gas == 1){
		#pragma omp parallel for private(id,ang) schedule(dynamic)
		for (id = 0; id < total; ++id){
			ang = gsangHe(v_vec[id],b_vec[id],rng2[id],rng3[id],rng4[id],numberOfAtoms);
			temp_vec[id] = ((1.0 - cos(ang)) * b2max_vec[id] / float(imp));
			temp1_vec[id] = (1.5 * (sin(ang) * sin(ang)) * b2max_vec[id] / float(imp));
		}
	}else if (gas == 2){
		#pragma omp parallel for private(id,ang) schedule(dynamic)
		for (id = 0; id < total; ++id){
			ang = gsangN2(v_vec[id],b_vec[id],rng2[id],rng3[id],rng4[id],numberOfAtoms);
			temp_vec[id] = ((1.0 - cos(ang)) * b2max_vec[id] / float(imp));
			temp1_vec[id] = (1.5 * (sin(ang) * sin(ang)) * b2max_vec[id] / float(imp));
		}
	}

	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
	for(int i = 0; i < total / imp; i++){
	  double count = 0, count1 = 0;
	  #if defined(__INTEL_COMPILER)
      #pragma vector aligned
	  #endif
	  for (int j = 0; j < imp; j++){
	     count += temp_vec[i * imp + j];
	     count1 += temp1_vec[i * imp + j];
	  }
	  temp_vec2[i] = count;
	  temp1_vec2[i] = count1;
	}

id= 0;

	for (ic = 0; ic < itn; ++ic){
		#if defined(__INTEL_COMPILER)
    	#pragma vector aligned
		#endif
		for (ig = 0; ig < inp; ++ig){
			om11st[ic] += temp_vec2[id] * wgst[ig];
			om12st[ic] += temp_vec2[id] * pgst[ig] * pgst[ig] * wgst[ig] * (1.0 / (3.0 * tst));
			om13st[ic] += temp_vec2[id] * (pow(pgst[ig],4)) * wgst[ig] * (1.0 / (12.0 * tst * tst));
			om22st[ic] += temp1_vec2[id] * pgst[ig] * pgst[ig] * wgst[ig] * (1.0 / (3.0 * tst));
			q1st[ig] += temp_vec2[id];
			q2st[ig] += temp1_vec2[id];
			id++;
		}
	}
/* End of Monte Carlo */

    if(ifailc < ifail){
//   Calculate running averages

    hold1 = 0.0;
    hold2 = 0.0;
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (ic = 0; ic < itn; ++ic){
        temp = 1.0 / (mconst / (sqrt(temperature) * om11st[ic] * pi * ro * ro));
        hold1 += om11st[ic];
        hold2 += temp;
    }

    mom11st = 0.0;
    mom12st = 0.0;
    mom13st = 0.0;
    mom22st = 0.0;
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (ic = 0; ic < itn; ++ic){
        mom11st += om11st[ic];
        mom12st += om12st[ic];
        mom13st += om13st[ic];
        mom22st += om22st[ic];
    }

    mom11st /= float(itn);
    mom12st /= float(itn);
    mom13st /= float(itn);
    mom22st /= float(itn);
    sdom11st = 0.0;
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (ic = 0; ic < itn; ++ic){
        hold = mom11st - om11st[ic];
        sdom11st += (hold * hold);
    }

    sdom11st = sqrt(sdom11st / float(itn));
    //sterr = sdom11st / sqrt(float(itn));

    cs = mom11st * pi * ro * ro;
    sdevpc = (100.0 * sdom11st)/mom11st;

//   Use omegas to obtain higher order correction factor to mobility
     ayst = mom22st / mom11st;
     best = ((5.0 * mom12st) - (4.0 * mom13st)) / mom11st;
     cest = mom12st/mom11st;
     term = ((4.0 * ayst) / (15.0)) + (0.5 * (pow((m2-m1),2)) / (m1 * m2));
     u2 = term - (0.08333 * (2.4 * best + 1.0) * (m1 / m2));
     w = (m1 / m2);
     delta = ((pow(((6.0 * cest)-5.0),2)) * w) / (60.0 * (1.0 + u2));
     f = 1.0 / (1.0 - delta);
     mob = (mconst * f)/(sqrt(temperature) * cs);
    }

      //Free allocated variables
      delete[] q1st;
	  delete[] q2st;
	  delete[] pgst;
	  delete[] wgst;
	  delete[] b2max;
	  delete[] om11st;
	  delete[] om12st;
	  delete[] om13st;
	  delete[] om22st;

	  delete[] rng1;
	  delete[] rng2;
	  delete[] rng3;
	  delete[] rng4;

	  delete[] temp_vec;
	  delete[] temp1_vec;
	  delete[] v_vec;
	  delete[] b_vec;
	  delete[] b2max_vec;
	  delete[] temp_vec2;
	  delete[] temp1_vec2;

    return 0;

}

