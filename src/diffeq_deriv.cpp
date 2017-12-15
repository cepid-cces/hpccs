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
#include <cmath>

double diffeq(int l, double tim, double dt, double w[6], double dw[10],  double array[6][6], double q[6]){

/*
 * Integration subroutine - uses 5th order runge-kutta-gill to
 * initiate and 5th order adams-moulton predictor-corrector to
 * propagate. Parameter l is initially set to zero and then
 * incremented to tell the subroutine when to switch between
 * integration methods. DIFFEQ calls subroutine DERIV to define
 * the equations of motion to be integrated.
 */

    double savw[6], savdw[6], r;
    static double a[4] = {0.50, 0.292893218814, 1.70710678118, 0.1666666666667};
    static double b[4] = {2.0, 1.0, 1.0, 2.0};
    static double c[4] = {-0.5,-0.292893218814,-1.70710678118,-0.5};
    static double ampc[5] = {-0.111059153612,0.672667757774,-1.70633621697,2.33387888707,-1.8524668225};
    static double amcc[4] = {0.0189208128941,-0.121233356692,0.337771548703,-0.55921513665};
    static double acst = 0.332866152768;
    static double var = 2.97013888888, cvar = 0.990972222222;

    int i, j, k;

    if (l < -1){

/*  This is the adams-moulton predictor-corrector part. */

    	for (j = 0; j < 6; j++){
            savw[j] = w[j];
            savdw[j] = dw[j];
            array[5][j] = savdw[j];
			#if defined(__INTEL_COMPILER)
    		#pragma vector aligned
			#endif
			for (int i = 0; i < 5; i++){
                array[5][j] += ampc[i] * array[i][j];
			}
            w[j] += (array[5][j] * hvar);
        }
    	tim += dt;
        deriv(w,dw);
        for (j = 0; j < 6; j++){
            array[5][j] = acst * dw[j];
            for (i = 0; i < 4; i++){
                array[i][j] = array[i+1][j];
                array[5][j] += array[i][j] * amcc[i];
            }
            array[4][j] = savdw[j];
            w[j] = savw[j] + hcvar * (array[4][j] + array[5][j]);
        }
        deriv(w,dw);
    }else if (l > -1){
           l++;


    /*    This is the runge-kutta-gill part...the steps are broken up into
     *     half steps to improve accuracy.
     */

       for (k = 0; k < 2; k++){
            for (j = 0; j < 4; j++){
                if (pow(-1,j+1) > 0) tim += 0.5 * dt;
                deriv(w,dw);
				#if defined(__INTEL_COMPILER)
				#pragma vector aligned
				#endif
                for (i = 0; i < 6; i++){
                    dw[i] *= dt;
                }
				#if defined(__INTEL_COMPILER)
    			#pragma vector aligned
				#endif
                for (i = 0; i < 6; i++){
                    r = a[j] * (dw[i] -b[j] * q[i]);
                    w[i] += r;
                    q[i] += 3.0 * r + c[j] * dw[i];
                }
            }
            deriv(w,dw);
       }


        if (l-6 < -1){
            for (j = 0; j < 6; j++){
              array[l][j] = dw[j];
            }
        }else{
                l = -2;
                dt *= 2.0;
        }
        dw[8] = l;
        dw[9] = dt;
    }else{//l = 0;
        	for (j = 0; j < 6; j++){
                q[j] = 0.0;
        	}
            hvar = dt * var;
            hcvar = dt * cvar;
            dt = 0.5 * dt;
            l++;

			/*
			 *   This is the runge-kutta-gill part...the steps are broken up into
			 *   half steps to improve accuracy.
             */

			for (k = 0; k < 2; k++){
				for (j = 0; j < 4; j++){
						if (pow(-1,j+1) > 0) tim += 0.5 * dt;
						deriv(w,dw);
						#if defined(__INTEL_COMPILER)
    					#pragma vector aligned
						#endif
						for (i = 0; i < 6; i++){
							dw[i] *= dt;
						}
						#if defined(__INTEL_COMPILER)
    					#pragma vector aligned
						#endif
						for (i = 0; i < 6; i++){
							r = a[j] * (dw[i] - b[j] * q[i]);
							w[i] += r;
							q[i] += 3.0 * r + c[j] * dw[i];
						}
				}
					deriv(w,dw);
			}

				if (l - 6 < -1){
					 for (j = 0; j < 6; j++){
						 array[l][j] = dw[j];
					 }
				} else{
						l = -2;
						dt = 2.0 * dt;
				}

				dw[8] = l;
				dw[9] = dt;
        	}

		return tim;
	}

void deriv(double *w,double *dw){

	double pot, dmax, dpotx, dpoty, dpotz;


/*
 *   Defines Hamilton's equations of motion as the time derivatives
 *   of the coordinates and momenta.
 *   x = w[0]; y = w[2]; z = w[4];
 *   From Hamilton's equations, the time derivatives of the coordinates
 *   are the conjugates divided by the mass.
 *
 */

     dw[0] = w[1] / mu;
     dw[2] = w[3] / mu;
     dw[4] = w[5] / mu;

/*
 *    Hamilton's equations for the time derivatives of the momenta
 *    evaluated by using the coordinate derivatives together with the
 *    chain rule.
 *
 */

/*   These are analytical derivatives. */
    if (gas == 2) dljpotN2(w[0],w[2],w[4], &pot, &dmax, &dpotx, &dpoty, &dpotz);
    if (gas == 1) dljpotHe(w[0],w[2],w[4], &pot, &dmax, &dpotx, &dpoty, &dpotz);

     dw[1] = -dpotx;
     dw[3] = -dpoty;
     dw[5] = -dpotz;
     dw[6] = pot;
     dw[7] = dmax;


}

double diffeq(int l, double tim, double dt, double w[6], double dw[10],  double array[6][6], double q[6], double *fx, double *fy, double *fz, double *phvar, double *phcvar){

/*
 * Integration subroutine - uses 5th order runge-kutta-gill to
 * initiate and 5th order adams-moulton predictor-corrector to
 * propagate. Parameter l is initially set to zero and then
 * incremented to tell the subroutine when to switch between
 * integration methods. DIFFEQ calls subroutine DERIV to define
 * the equations of motion to be integrated.
 */

    double savw[6], savdw[6], r;
    static double a[4] = {0.50, 0.292893218814, 1.70710678118, 0.1666666666667};
    static double b[4] = {2.0, 1.0, 1.0, 2.0};
    static double c[4] = {-0.5,-0.292893218814,-1.70710678118,-0.5};
    static double ampc[5] = {-0.111059153612,0.672667757774,-1.70633621697,2.33387888707,-1.8524668225};
    static double amcc[4] = {0.0189208128941,-0.121233356692,0.337771548703,-0.55921513665};
    static double acst = 0.332866152768;
    static double var = 2.97013888888, cvar = 0.990972222222;

    int i, j, k;

    if (l < -1){
//*  This is the adams-moulton predictor-corrector part. */
    	for (j = 0; j < 6; j++){
            savw[j] = w[j];
            savdw[j] = dw[j];
            array[5][j] = savdw[j];
			#if defined(__INTEL_COMPILER)
    		#pragma vector aligned
			#endif
			for (int i = 0; i < 5; i++){
                array[5][j] += ampc[i] * array[i][j];
			}
            w[j] += (array[5][j] * *phvar);
        }
    	tim += dt;
        deriv(w,dw,fx,fy,fz);
        for (j = 0; j < 6; j++){
            array[5][j] = acst * dw[j];
            for (i = 0; i < 4; i++){
                array[i][j] = array[i+1][j];
                array[5][j] = array[i][j] * amcc[i] + array[5][j];
            }
            array[4][j] = savdw[j];
            w[j] = savw[j] + *phcvar * (array[4][j] + array[5][j]);
        }
        deriv(w,dw,fx,fy,fz);
    }else if (l > -1){
           l++;
//
//
//    /*    This is the runge-kutta-gill part...the steps are broken up into
//     *     half steps to improve accuracy.
//     */
//
       for (k = 0; k < 2; k++){
            for (j = 0; j < 4; j++){
                if (pow(-1,j+1) > 0) tim += 0.5 * dt;
                deriv(w,dw,fx,fy,fz);
				#if defined(__INTEL_COMPILER)
    			#pragma vector aligned
				#endif
                for (i = 0; i < 6; i++){
                    dw[i] *= dt;

                }
				#if defined(__INTEL_COMPILER)
    			#pragma vector aligned
				#endif
                for (i = 0; i < 6; i++){
                	r = a[j] * (dw[i] -b[j] * q[i]);
                	w[i] += r;
                	q[i] += 3.0 * r + c[j] * dw[i];
                }
            }
            deriv(w,dw,fx,fy,fz);
       }


        if (l-6 < -1){
            for (j = 0; j < 6; j++){
              array[l][j] = dw[j];
            }

        }else{
                l = -2;
                dt *= 2.0;
        }
        dw[8] = l;
        dw[9] = dt;
    }else{//l = 0;
        	for (j = 0; j < 6; j++){
                q[j] = 0.0;
        	}
            *phvar = dt * var;
            *phcvar = dt * cvar;
            dt = 0.5 * dt;
            l++;

			/*
			 *   This is the runge-kutta-gill part...the steps are broken up into
			 *   half steps to improve accuracy.
             */

			for (k = 0; k < 2; k++){
				for (j = 0; j < 4; j++){
						if (pow(-1,j+1) > 0) tim += 0.5 * dt;
						deriv(w,dw,fx,fy,fz);
						#if defined(__INTEL_COMPILER)
    					#pragma vector aligned
						#endif
						for (i = 0; i < 6; i++){
							dw[i] *= dt;
						}
						#if defined(__INTEL_COMPILER)
    					#pragma vector aligned
						#endif
						for (i = 0; i < 6; i++){
							r = a[j] * (dw[i] - b[j] * q[i]);
							w[i] += r;
							q[i] += 3.0 * r + c[j] * dw[i];
						}
				}
				 deriv(w,dw,fx,fy,fz);
			}

				if (l - 6 < -1){
					 for (j = 0; j < 6; j++){
						 array[l][j] = dw[j];
					 }
				} else{
						l = -2;
						dt = 2.0 * dt;
				}

				dw[8] = l;
				dw[9] = dt;
        	}

		return tim;
	}

void deriv(double *w,double *dw, double *fx, double *fy, double *fz){

	double pot, dmax, dpotx, dpoty, dpotz;


/*
 *   Defines Hamilton's equations of motion as the time derivatives
 *   of the coordinates and momenta.
 *   x = w[0]; y = w[2]; z = w[4];
 *   From Hamilton's equations, the time derivatives of the coordinates
 *   are the conjugates divided by the mass.
 *
 */

     dw[0] = w[1] / mu;
     dw[2] = w[3] / mu;
     dw[4] = w[5] / mu;

/*
 *    Hamilton's equations for the time derivatives of the momenta
 *    evaluated by using the coordinate derivatives together with the
 *    chain rule.
 *
 */

/*   These are analytical derivatives. */
     if (gas == 2) dljpotN2(w[0],w[2],w[4], &pot, &dmax, &dpotx, &dpoty, &dpotz, fx, fy, fz);
     if (gas == 1) dljpotHe(w[0],w[2],w[4], &pot, &dmax, &dpotx, &dpoty, &dpotz, fx, fy, fz);


     dw[1] = -dpotx;
     dw[3] = -dpoty;
     dw[5] = -dpotz;
     dw[6] = pot;
     dw[7] = dmax;


}
