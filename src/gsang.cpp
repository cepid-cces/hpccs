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
#include <headers/diffeq_deriv.hpp>
#include <headers/potentialHe.hpp>
#include <headers/potentialN2.hpp>
#include <headers/rotate.hpp>
#include <iostream>
#include <cmath>
#if defined(__INTEL_COMPILER)
#include <aligned_new>
#else
#include <new>
#endif

using namespace std;

double gsangHe(double v,double x, double rnt, double rnp, double rng, unsigned int totalAtoms){

	/*
	 *    Calculates trajectory. Adapted from code written by George Schatz
	 *    A fixed step integrator is used which combines a runge-kutta-gill
	 *    initiator with an adams-moulton predictor-corrector propagator.
	 *
	*/
		int ns,nw,l,id2,iymax;
	    double w[6] __attribute__((aligned(16)));
	    double dw[10] __attribute__((aligned(16)));
	    double array[6][6]__attribute__((aligned(16)));
	    double q[6] __attribute__((aligned(16)));
	    double y,z, pot, dmax;
	    double vy, vx, vz, dt, e;
	    double den, iymin, top, num, ymax, etot, ymin;
	    double tim, ang, erat, dt1, dt2, e0;
	    double *fx __attribute__((aligned(16)));
	    double *fy __attribute__((aligned(16)));
	    double *fz __attribute__((aligned(16)));
	    double phvar, phcvar;
	    fx = new double[totalAtoms];
	    fy = new double[totalAtoms];
	    fz = new double[totalAtoms];

	    rotate(rnt,rnp,rng,fx,fy,fz);

	    top = (v / 95.2381) - 0.5;
		if(v >= 1000.0 && v < 2000.0){
			top = 10.0;
		} else if(v >= 2000.0 && v < 3000.0){
			top = 10.0 - ((v - 2000.0) * 7.5e-3);
		} else if(v >= 3000.0){
			top = 2.5;
		}
		dt1 = top * dtsf1 * 1.0e-11 / v;
		dt2 = dt1 * dtsf2;
		e0 = 0.5 * mu * v * v;

	    phvar = 0.0;
	    phcvar = 0.0;
	    vy = -v;
	    vx = 0.0;
	    vz = 0.0;
	    ang = 0.0;

	//  determine time step
	    dt = dt1;
	   //  determine trajectory start position
	    z = 0.0;
	    ymin = 0.0;
	    ymax = 0.0;

	    ymax /= 1.0e-10;
	    ymin /= 1.0e-10;
	    iymin = int(ymin) -1;
	    iymax = int(ymax) +1;

	    id2 = iymax;
	    y = float(id2) * 1.0e-10;
	    pot = dljpotHe(x,y,z,fx,fy,fz);
	    if(abs(pot/e0) > sw1){
	      do{
	      id2 += 10;
	      y = float(id2) * 1.0e-10;
	      pot = dljpotHe(x,y,z,fx,fy,fz);
	      } while(abs(pot / e0) > sw1);
	      do{
	         id2--;
	         y = float(id2) * 1.0e-10;
	         pot = dljpotHe(x,y,z,fx,fy,fz);
	      }while(abs(pot/e0) < sw1);
	    }else{
	    	do{
	    		id2--;
	    	    y = float(id2) * 1.0e-10;
	    	     pot = dljpotHe(x,y,z,fx,fy,fz);
	    	    if(id2 < iymin){
	    	       cout << "trajectory not started - potential too small" << "\n";
	    	       ang = 0.0;
	    	       erat = 1.000;
	    	       delete[] fx;
	    	       delete[] fy;
	    	       delete[] fz;
	    	       return ang;
	    	    }

				//pot = dljpotHe(x,y,z,fx,fy,fz);

	    	}while(fabs(pot / e0) < sw1);
	    }
	    etot = e0 + pot;

	//  initial coordinates and momenta
	    w[0] = x;
	    w[1] = vx * mu;
	    w[2] = y;
	    w[3] = vy * mu;
	    w[4] = z;
	    w[5] = vz * mu;
	    tim = 0.0;

	//  initialize the time derivatives of the coordinates and momenta
	    deriv(w,dw,fx,fy,fz);
	    ns = 0;
	    nw = 0;
	    l = -1;

	    do{
	    	tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
	    	l = dw[8];
	    	dt = dw[9];
	    	nw++;
	    }while (nw != inwr);



	    pot = dw[6];
	    dmax = dw[7];
	    ns += nw;
	    nw = 0;

	//  check if trajectory has become "lost" (too many steps)
	    if(ns > 29999){
	        ang = pi / 2.0;
	        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	        erat = (e + pot) / etot;
	        delete[] fx;
	        delete[] fy;
	        delete[] fz;
	        return ang;
	    }

	//  check if the trajectory is finished
	    while (dmax < romax){
			do{
				tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
				l = dw[8];
				dt = dw[9];
				nw++;
			}while (nw != inwr);
				pot = dw[6];
				dmax = dw[7];
				ns += nw;
				nw = 0;
			//  check if trajectory has become "lost" (too many steps)
				if(ns > 29999){
					ang = pi / 2.0;
					e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
					erat = (e + pot) / etot;
					delete[] fx;
					delete[] fy;
					delete[] fz;
					return ang;
				}
	    }

	   if(fabs(pot / e0) > sw2 && dt == dt1){
	        dt = dt2;
	        l = -1;
	    }

	    if(fabs(pot / e0) < sw2 && dt == dt2){
	        dt = dt1;
	        l = -1;
	    }

	    while(fabs(pot / e0) > sw1 || ns < 50){
	    	  do{
	    		  tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
	    	    	l = dw[8];
	    	    	dt = dw[9];
	    	    	nw++;
	    	    }while (nw != inwr);

	    	    pot = dw[6];
	    	    dmax = dw[7];
	    	    ns += nw;
	    	    nw = 0;

	    	//  check if trajectory has become "lost" (too many steps)
	    	    if(ns > 29999){
	    	        ang = pi / 2.0;
	    	        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	    	        erat = (e + pot) / etot;
	    	        delete[] fx;
	    	        delete[] fy;
	    	        delete[] fz;
	    	        return ang;
	    	    }

	    	//  check if the trajectory is finished
	    	    while (dmax < romax){
	    	    		do{
	    	    			tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
	    	    		    l = dw[8];
	    	    		    dt = dw[9];
	    	    		   	nw++;
	    	    		}while (nw != inwr);
	    	    			pot = dw[6];
	    	    		    dmax = dw[7];
	    	    		    ns += nw;
	    	    		    nw = 0;
	    	    		//  check if trajectory has become "lost" (too many steps)
	    	    		    if(ns > 29999){
	    	    		        ang = pi / 2.0;
	    	    		        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	    	    		        erat = (e + pot) / etot;
	    	    		        delete[] fx;
	    	    		        delete[] fy;
	    	    		        delete[] fz;
	    	    		        return ang;
	    	    		    }
	    	    }
	    	   if(fabs(pot / e0) > sw2 && dt == dt1){
	    	        dt = dt2;
	    	        l = -1;
	    	    }

	    	    if(fabs(pot / e0) < sw2 && dt == dt2){
	    	    	dt = dt1;
	    	        l = -1;
	    	    }
	    }

	//  determine scattering angle

	    if(dw[0] > 0.0){
	        num = dw[2] * (-v);
	        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	        ang  = acos(num / den);
	    }else if(dw[0] < 0.0){
	        num = dw[2] * (-v);
	        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	        ang = (-acos(num / den));
	    }

	//  check for energy conservation
	    e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	    erat = (e + pot) / etot;
	    //cout << scientific << erat << "\n";
	    if(erat < 1.01 && erat > 0.99) {
	        delete[] fx;
	        delete[] fy;
	        delete[] fz;
	    	return ang;
	    }
		#pragma omp atomic
	    trajlost++;
	    delete[] fx;
	    delete[] fy;
	    delete[] fz;
	    return ang;
}

double gsangN2(double v,double x,double rnt, double rnp, double rng, unsigned int totalAtoms){

	/*
	 *    Calculates trajectory. Adapted from code written by George Schatz
	 *    A fixed step integrator is used which combines a runge-kutta-gill
	 *    initiator with an adams-moulton predictor-corrector propagator.
	 *
	*/
		int ns,nw,l,id2,iymax;
	    double w[6] = {0}, dw[10] = {0.0}, array[6][6]={0}, q[6]={0};
	    double y,z, pot, dmax;
	    double vy, vx, vz, dt, e;
	    double den, iymin, top, num, ymax, etot, ymin;
	    double tim, ang, erat,dt1,dt2,e0;
	    double *fx __attribute__((aligned(16)));
	    double *fy __attribute__((aligned(16)));
	    double *fz __attribute__((aligned(16)));
	    double phvar, phcvar;

	    fx = new double[totalAtoms];
	    fy = new double[totalAtoms];
	    fz = new double[totalAtoms];

	    rotate(rnt,rnp,rng,fx,fy,fz);

	    phvar = 0.0;
	    phcvar = 0.0;
	    vy = -v;
	    vx = 0.0;
	    vz = 0.0;
	    ang = 0.0;
	    top = (v / 95.2381) - 0.5;
		if(v >= 1000.0 && v < 2000.0){
			top = 10.0;
		} else if(v >= 2000.0 && v < 3000.0){
			top = 10.0 - ((v - 2000.0) * 7.5e-3);
		} else if(v >= 3000.0){
			top = 2.5;
		}
		dt1 = top * dtsf1 * 1.0e-11 / v;
		dt2 = dt1 * dtsf2;
		e0 = 0.5 * mu * v * v;
	//  determine time step
	    dt = dt1;
	   //  determine trajectory start position
	    z = 0.0;
	    ymin = 0.0;
	    ymax = 0.0;

	    ymax /= 1.0e-10;
	    ymin /= 1.0e-10;
	    iymin = int(ymin) -1;
	    iymax = int(ymax) +1;

	    id2 = iymax;
	    y = float(id2) * 1.0e-10;

	    pot = dljpotN2(x,y,z,fx,fy,fz);


	    if(fabs(pot/e0) > sw1){
	      do{
	      id2 += 10;
	      y = float(id2) * 1.0e-10;

	      pot = dljpotN2(x,y,z,fx,fy,fz);

	      } while(fabs(pot / e0) > sw1);

	      do{
	         id2--;
	         y = float(id2) * 1.0e-10;

	         pot = dljpotN2(x,y,z,fx,fy,fz);


	      }while(fabs(pot/e0) < sw1);
	    }else{
	    	do{
	    		id2--;
	    	    y = float(id2) * 1.0e-10;

	    	     pot = dljpotN2(x,y,z,fx,fy,fz);

	    	    if(id2 < iymin){
	    	       cout << "trajectory not started - potential too small" << "\n";
	    	       ang = 0.0;
	    	       erat = 1.000;
	    	       delete[] fx;
	    	       delete[] fy;
	    	       delete[] fz;
	    	       return ang;
	    	    }

				pot = dljpotN2(x,y,z,fx,fy,fz);

	    	}while(fabs(pot / e0) < sw1);
	    }
	    etot = e0 + pot;

	//  initial coordinates and momenta
	    w[0] = x;
	    w[1] = vx * mu;
	    w[2] = y;
	    w[3] = vy * mu;
	    w[4] = z;
	    w[5] = vz * mu;
	    tim = 0.0;

	//  initialize the time derivatives of the coordinates and momenta
	    deriv(w,dw,fx,fy,fz);
	    ns = 0;
	    nw = 0;
	    l = -1;

	    do{
	    	tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
	    	l = dw[8];
	    	dt = dw[9];
	    	nw++;
	    }while (nw != inwr);



	    pot = dw[6];
	    dmax = dw[7];
	    ns += nw;
	    nw = 0;

	//  check if trajectory has become "lost" (too many steps)
	    if(ns > 29999){
	        ang = pi / 2.0;
	        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	        erat = (e + pot) / etot;
	        delete[] fx;
	        delete[] fy;
	        delete[] fz;
	        return ang;
	    }

	//  check if the trajectory is finished
	    while (dmax < romax){
			do{
				tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
				l = dw[8];
				dt = dw[9];
				nw++;
			}while (nw != inwr);
				pot = dw[6];
				dmax = dw[7];
				ns += nw;
				nw = 0;
			//  check if trajectory has become "lost" (too many steps)
				if(ns > 29999){
					ang = pi / 2.0;
					e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
					erat = (e + pot) / etot;
					delete[] fx;
					delete[] fy;
					delete[] fz;
					return ang;
				}
	    }

	   if(fabs(pot / e0) > sw2 && dt == dt1){
	        dt = dt2;
	        l = -1;
	    }

	    if(fabs(pot / e0) < sw2 && dt == dt2){
	        dt = dt1;
	        l = -1;
	    }

	    while(fabs(pot / e0) > sw1 || ns < 50){
	    	  do{
	    		  tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
	    	    	l = dw[8];
	    	    	dt = dw[9];
	    	    	nw++;
	    	    }while (nw != inwr);

	    	    pot = dw[6];
	    	    dmax = dw[7];
	    	    ns += nw;
	    	    nw = 0;

	    	//  check if trajectory has become "lost" (too many steps)
	    	    if(ns > 29999){
	    	        ang = pi / 2.0;
	    	        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	    	        erat = (e + pot) / etot;
	    	        delete[] fx;
	    	        delete[] fy;
	    	        delete[] fz;
	    	        return ang;
	    	    }

	    	//  check if the trajectory is finished
	    	    while (dmax < romax){
	    	    		do{
	    	    			tim = diffeq(l,tim,dt,w,dw,array, q, fx, fy, fz, &phvar, &phcvar);
	    	    		    l = dw[8];
	    	    		    dt = dw[9];
	    	    		   	nw++;
	    	    		}while (nw != inwr);
	    	    			pot = dw[6];
	    	    		    dmax = dw[7];
	    	    		    ns += nw;
	    	    		    nw = 0;
	    	    		//  check if trajectory has become "lost" (too many steps)
	    	    		    if(ns > 29999){
	    	    		        ang = pi / 2.0;
	    	    		        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	    	    		        erat = (e + pot) / etot;
	    	    		        delete[] fx;
	    	    		        delete[] fy;
	    	    		        delete[] fz;
	    	    		        return ang;
	    	    		    }
	    	    }
	    	   if(fabs(pot / e0) > sw2 && dt == dt1){
	    	        dt = dt2;
	    	        l = -1;
	    	    }

	    	    if(fabs(pot / e0) < sw2 && dt == dt2){
	    	    	dt = dt1;
	    	        l = -1;
	    	    }
	    }

	//  determine scattering angle

	    if(dw[0] > 0.0){
	        num = dw[2] * (-v);
	        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	        ang  = acos(num / den);
	    }else if(dw[0] < 0.0){
	        num = dw[2] * (-v);
	        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	        ang = (-acos(num / den));
	    }

	//  check for energy conservation
	    e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
	    erat = (e + pot) / etot;
	    //cout << scientific << erat << "\n";
	    if(erat < 1.01 && erat > 0.99) {
	        delete[] fx;
	        delete[] fy;
	        delete[] fz;
	    	return ang;
	    }
		#pragma omp atomic
	    trajlost++;
	    delete[] fx;
	    delete[] fy;
	    delete[] fz;
	    return ang;
}




double gsangHe(const double v,const double x){

/*
 *    Calculates trajectory. Adapted from code written by George Schatz
 *    A fixed step integrator is used which combines a runge-kutta-gill
 *    initiator with an adams-moulton predictor-corrector propagator.
 *
*/
	//int it;
	int ns,nw,l,id2,iymax;
	double w[6] __attribute__((aligned(16)));
	double dw[10] __attribute__((aligned(16)));
	double array[6][6]__attribute__((aligned(16)));
	double q[6] __attribute__((aligned(16)));
    double y,z, pot, dmax;
    double vy, vx, vz, dt, dt1, dt2, e, e0;
    double den, iymin, top, num, ymax, etot, ymin;
    double tim, ang, erat;

    hvar = 0.0;
    hcvar = 0.0;
    vy = -v;
    vx = 0.0;
    vz = 0.0;
    ang = 0.0;

//  determine time step
    top = (v / 95.2381) - 0.5;
    if(v >= 1000.0 && v < 2000.0){
        top = 10.0;
    } else if(v >= 2000.0 && v < 3000.0){
        top = 10.0 - ((v - 2000.0) * 7.5e-3);
    } else if(v >= 3000.0){
        top = 2.5;
    }

    dt1 = top * dtsf1 * 1.0e-11 / v;
    dt2 = dt1 * dtsf2;
    dt = dt1;

//  determine trajectory start position

    e0 = 0.5 * mu * v * v;
    z = 0.0;
    ymin = 0.0;
    ymax = 0.0;
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int i=0; i < numberOfAtoms; i++){
        if(molecule->fy[i] > ymax) ymax = molecule->fy[i];
        if(molecule->fy[i] < ymin) ymin = molecule->fy[i];
    }

    ymax /= 1.0e-10;
    ymin /= 1.0e-10;
    iymin = int(ymin) - 1;
    iymax = int(ymax) + 1;
    id2 = iymax;
    y = float(id2) * 1.0e-10;
    pot = dljpotHe(x,y,z);

    if(fabs(pot/e0) > sw1){
      do{
      id2 += 10;
      y = float(id2) * 1.0e-10;

      	pot = dljpotHe(x,y,z);

      } while(fabs(pot / e0) > sw1);

      do{
         id2--;
         y = float(id2) * 1.0e-10;
         	pot = dljpotHe(x,y,z);

      }while(fabs(pot/e0) < sw1);
    }else{
    	do{
    		id2--;
    	    y = float(id2) * 1.0e-10;
    	    	pot = dljpotHe(x,y,z);

    	    if(id2 < iymin){
    	       cout << "trajectory not started - potential too small" << "\n";
    	       ang = 0.0;
    	       erat = 1.000;
    	    }
				pot = dljpotHe(x,y,z);
    	}while(fabs(pot / e0) < sw1);
    }

    etot = e0 + pot;

//  initial coordinates and momenta

    w[0] = x;
    w[1] = vx * mu;
    w[2] = y;
    w[3] = vy * mu;
    w[4] = z;
    w[5] = vz * mu;
    tim = 0.0;

//  initialize the time derivatives of the coordinates and momenta
    deriv(w,dw);
    ns = 0;
    nw = 0;
    l = -1;

    do{
    	tim = diffeq(l,tim,dt,w,dw,array, q);
    	l = dw[8];
    	dt = dw[9];
    	nw++;
    }while (nw != inwr);

    pot = dw[6];
    dmax = dw[7];
    ns += nw;
    nw = 0;
//  check if trajectory has become "lost" (too many steps)
    if(ns > 29999){
        ang = pi / 2.0;
        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
        erat = (e + pot) / etot;
        return ang;
    }

//  check if the trajectory is finished
    while (dmax < romax){
    		do{
    		    tim = diffeq(l,tim,dt,w,dw,array,q);
    		    l = dw[8];
    		    dt = dw[9];
    		   	nw++;
    		}while (nw != inwr);
    			pot = dw[6];
    		    dmax = dw[7];
    		    ns += nw;
    		    nw = 0;
    		//  check if trajectory has become "lost" (too many steps)
    		    if(ns > 29999){
    		        ang = pi / 2.0;
    		        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    		        erat = (e + pot) / etot;
    		        return ang;
    		    }
    }
   if(fabs(pot / e0) > sw2 && dt == dt1){
        dt = dt2;
        l = -1;
    }

    if(fabs(pot / e0) < sw2 && dt == dt2){
        dt = dt1;
        l = -1;
    }

    while(fabs(pot / e0) > sw1 || ns < 50){
    	  do{
    	    	tim = diffeq(l,tim,dt,w,dw,array,q);
    	    	l = dw[8];
    	    	dt = dw[9];
    	    	nw++;
    	    }while (nw != inwr);

    	    pot = dw[6];
    	    dmax = dw[7];
    	    ns += nw;
    	    nw = 0;

    	//  check if trajectory has become "lost" (too many steps)
    	    if(ns > 29999){
    	        ang = pi / 2.0;
    	        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    	        erat = (e + pot) / etot;
    	        return ang;
    	    }

    	//  check if the trajectory is finished
    	    while (dmax < romax){
    	    		do{
    	    		    tim = diffeq(l,tim,dt,w,dw,array,q);
    	    		    l = dw[8];
    	    		    dt = dw[9];
    	    		   	nw++;
    	    		}while (nw != inwr);
    	    			pot = dw[6];
    	    		    dmax = dw[7];
    	    		    ns += nw;
    	    		    nw = 0;
    	    		//  check if trajectory has become "lost" (too many steps)
    	    		    if(ns > 29999){
    	    		        ang = pi / 2.0;
    	    		        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    	    		        erat = (e + pot) / etot;
    	    		        return ang;
    	    		    }
    	    }

    	   if(fabs(pot / e0) > sw2 && dt == dt1){
    	        dt = dt2;
    	        l = -1;
    	    }

    	    if(fabs(pot / e0) < sw2 && dt == dt2){
    	    	dt = dt1;
    	        l = -1;
    	    }

    }

//  determine scattering angle

    if(dw[0] > 0.0){
        num = dw[2] * (-v);
        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
        ang  = acos(num / den);
    }else if(dw[0] < 0.0){
        num = dw[2] * (-v);
        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
        ang = (-acos(num / den));
    }

//  check for energy conservation
    e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    erat = (e + pot) / etot;
    if(erat < 1.01 && erat > 0.99) {
    	return ang;
    }

    trajlost++;

    return ang;

}




double gsangN2(const double v,const double x){

/*
 *    Calculates trajectory. Adapted from code written by George Schatz
 *    A fixed step integrator is used which combines a runge-kutta-gill
 *    initiator with an adams-moulton predictor-corrector propagator.
 *
*/
	//int it;
	int ns,nw,l,id2,iymax;
    double w[6] __attribute__((aligned(16)));
    double dw[10] __attribute__((aligned(16)));
    double array[6][6]__attribute__((aligned(16)));
    double q[6] __attribute__((aligned(16)));
    double y,z, pot, dmax;
    double vy, vx, vz, dt, dt1, dt2, e, e0;
    double den, iymin, top, num, ymax, etot, ymin;
    double tim, ang, erat;

    hvar = 0.0;
    hcvar = 0.0;
    vy = -v;
    vx = 0.0;
    vz = 0.0;
    ang = 0.0;

//  determine time step
    top = (v / 95.2381) - 0.5;
    if(v >= 1000.0 && v < 2000.0){
        top = 10.0;
    } else if(v >= 2000.0 && v < 3000.0){
        top = 10.0 - ((v - 2000.0) * 7.5e-3);
    } else if(v >= 3000.0){
        top = 2.5;
    }

    dt1 = top * dtsf1 * 1.0e-11 / v;
    dt2 = dt1 * dtsf2;
    dt = dt1;

//  determine trajectory start position

    e0 = 0.5 * mu * v * v;
    z = 0.0;
    ymin = 0.0;
    ymax = 0.0;
	#if defined(__INTEL_COMPILER)
    #pragma vector aligned
	#endif
    for (unsigned int i=0; i < numberOfAtoms; i++){
        if(molecule->fy[i] > ymax) ymax = molecule->fy[i];
        if(molecule->fy[i] < ymin) ymin = molecule->fy[i];
    }

    ymax /= 1.0e-10;
    ymin /= 1.0e-10;
    iymin = int(ymin) - 1;
    iymax = int(ymax) + 1;
    id2 = iymax;
    y = float(id2) * 1.0e-10;
    pot = dljpotN2(x,y,z);

    if(fabs(pot/e0) > sw1){
      do{
      id2 += 10;
      y = float(id2) * 1.0e-10;
      pot = dljpotN2(x,y,z);
      } while(fabs(pot / e0) > sw1);

      do{
         id2--;
         y = float(id2) * 1.0e-10;
         pot = dljpotN2(x,y,z);
      }while(fabs(pot/e0) < sw1);
    }else{
    	do{
    		id2--;
    	    y = float(id2) * 1.0e-10;
    	   	pot = dljpotN2(x,y,z);

    	    if(id2 < iymin){
    	       cout << "trajectory not started - potential too small" << "\n";
    	       ang = 0.0;
    	       erat = 1.000;
    	    }
			pot = dljpotN2(x,y,z);
    	}while(fabs(pot / e0) < sw1);
    }

    etot = e0 + pot;

//  initial coordinates and momenta

    w[0] = x;
    w[1] = vx * mu;
    w[2] = y;
    w[3] = vy * mu;
    w[4] = z;
    w[5] = vz * mu;
    tim = 0.0;

//  initialize the time derivatives of the coordinates and momenta
    deriv(w,dw);
    ns = 0;
    nw = 0;
    l = -1;

    do{
    	tim = diffeq(l,tim,dt,w,dw,array, q);
    	l = dw[8];
    	dt = dw[9];
    	nw++;
    }while (nw != inwr);

    pot = dw[6];
    dmax = dw[7];
    ns += nw;
    nw = 0;
//  check if trajectory has become "lost" (too many steps)
    if(ns > 29999){
        ang = pi / 2.0;
        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
        erat = (e + pot) / etot;
        return ang;
    }

//  check if the trajectory is finished
    while (dmax < romax){
    		do{
    		    tim = diffeq(l,tim,dt,w,dw,array,q);
    		    l = dw[8];
    		    dt = dw[9];
    		   	nw++;
    		}while (nw != inwr);
    			pot = dw[6];
    		    dmax = dw[7];
    		    ns += nw;
    		    nw = 0;
    		//  check if trajectory has become "lost" (too many steps)
    		    if(ns > 29999){
    		        ang = pi / 2.0;
    		        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    		        erat = (e + pot) / etot;
    		        return ang;
    		    }
    }
   if(fabs(pot / e0) > sw2 && dt == dt1){
        dt = dt2;
        l = -1;
    }

    if(fabs(pot / e0) < sw2 && dt == dt2){
        dt = dt1;
        l = -1;
    }

    while(fabs(pot / e0) > sw1 || ns < 50){
    	  do{
    	    	tim = diffeq(l,tim,dt,w,dw,array,q);
    	    	l = dw[8];
    	    	dt = dw[9];
    	    	nw++;
    	    }while (nw != inwr);

    	    pot = dw[6];
    	    dmax = dw[7];
    	    ns += nw;
    	    nw = 0;

    	//  check if trajectory has become "lost" (too many steps)
    	    if(ns > 29999){
    	        ang = pi / 2.0;
    	        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    	        erat = (e + pot) / etot;
    	        return ang;
    	    }

    	//  check if the trajectory is finished
    	    while (dmax < romax){
    	    		do{
    	    		    tim = diffeq(l,tim,dt,w,dw,array,q);
    	    		    l = dw[8];
    	    		    dt = dw[9];
    	    		   	nw++;
    	    		}while (nw != inwr);
    	    			pot = dw[6];
    	    		    dmax = dw[7];
    	    		    ns += nw;
    	    		    nw = 0;
    	    		//  check if trajectory has become "lost" (too many steps)
    	    		    if(ns > 29999){
    	    		        ang = pi / 2.0;
    	    		        e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    	    		        erat = (e + pot) / etot;
    	    		        return ang;
    	    		    }
    	    }

    	   if(fabs(pot / e0) > sw2 && dt == dt1){
    	        dt = dt2;
    	        l = -1;
    	    }

    	    if(fabs(pot / e0) < sw2 && dt == dt2){
    	    	dt = dt1;
    	        l = -1;
    	    }

    }

//  determine scattering angle

    if(dw[0] > 0.0){
        num = dw[2] * (-v);
        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
        ang  = acos(num / den);
    }else if(dw[0] < 0.0){
        num = dw[2] * (-v);
        den = v * sqrt(dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
        ang = (-acos(num / den));
    }

//  check for energy conservation
    e = 0.5 * mu * (dw[0] * dw[0] + dw[2] * dw[2] + dw[4] * dw[4]);
    erat = (e + pot) / etot;
    if(erat < 1.01 && erat > 0.99) {
    	return ang;
    }

    trajlost++;

    return ang;

}



