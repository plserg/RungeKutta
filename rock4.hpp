  //    Numerical solution of a (mildly) stiff system
  //    of first order differential equations. ROCK4 is
  //    based on a family of 4th order explicit Runge-
  //    Kutta methods with nearly optimal stability domain
  //    on the negative real axis. The numerical method is
  //    based on a three-term recursion relation, and a
  //    composition of two sub-methods.
  //    The size (along the negative axis) of the stability
  //    domains increases quadratically with the stage number.
  //
  //    Intended for problems of large dimensions with
  //    eigenvalues of the Jacobian close to the negative
  //    real axis. Typically for problems originating from
  //    parabolic PDEs. 
  //       C++ code is based on work of A. Abdulle.
  //     The analysis of the ROCK4 method is described in:
  //
  //    [ ] A. Abdulle
  //        Fourth order Chebyshev methods with
  //        recurrence relation

  /*
     send bug reports and comments to sergey.plyasunov@gmail.com
   */

#include <iostream>
#include <cstring>

#include "rk_utils.hpp"
#include "rock4_data.h"

namespace rk {

  template < typename TRHSFun > void rock4_step(const int neqn, 
						const double t,
						const double *y, 
						const double h, 
						const TRHSFun  fun,   //void (*fun)(const int,const  double,const double*, double*),
                                                int mdeg,
                                                double *yn,
                                                double &err,
                                                const double atol,
                                                const double rtol, double *work
                                                /*temporal array */ ) {

    double ci1, ci2, ci3, temp1, temp2, temp3, temp4, temp5;
    int mr, mz, i, j;

  /*-- initializations --*/
     err = 0.0e0;
  /*-- select entries in the data arrays--*/
     mz = 0;
     mr = 1;
    for (i = 0; i < 50; ++i) {

      if (ms[i] >= mdeg) {
        mdeg = ms[i];
        mz = (i + 1);
        break;
      } else {
        mr += ms[i] * 2 - 1;
      }
    }
  /*--allocate temporary arrays*/
    double *yjm1 = work;        //new double[neqn];
    double *yjm2 = work + neqn; //new double[neqn];
    double *yjm3 = work + 2 * neqn;     //new double[neqn];
    double *yjm4 = work + 3 * neqn;     //new double[neqn];
    double *fn = work + 4 * neqn;       //new double[neqn];
    double *fnt = work + 5 * neqn;      //new double[neqn];

    /* -------- First stage.-------- */

    std::memcpy(yn, y, sizeof(double) * neqn);

    temp1 = h * recf[mr - 1];
    ci1 = t + temp1;
    ci2 = t + temp1;
    ci3 = t;
    fun(neqn, t, y, fn);

    for (i = 0; i < neqn; ++i) {

      yjm2[i] = y[i];
      yjm1[i] = y[i] + temp1 * fn[i];

      if (mdeg < 2) {
        yn[i] = yjm1[i];
      }

    }
    //--------- Stage for j=2..mdeg.--------
    for (i = 2; i <= mdeg; ++i) {

      temp1 = h * recf[mr + 2 * (i - 2) + 1 - 1];
      temp3 = -recf[mr + 2 * (i - 2) + 2 - 1];
      temp2 = 1.0e0 - temp3;
      fun(neqn, ci1, yjm1, fn);
      ci1 = temp1 + temp2 * ci2 + temp3 * ci3;

      for (j = 0; j < neqn; ++j) {
        yn[j] = temp1 * fn[j] + temp2 * yjm1[j] + temp3 * yjm2[j];
        // -------- Shift the value "y" for the next stage.--------
        if (i < mdeg) {
          yjm2[j] = yjm1[j];
          yjm1[j] = yn[j];
        }
      }
      ci3 = ci2;
      ci2 = ci1;
    }
    // -------- The finishing procedure (4-stage method).
    //-------- Stage 1.--------  
    temp1 = h * fpa[mz - 1][1 - 1];
    fun(neqn, ci1, yn, yjm1);

    for (j = 0; j < neqn; ++j) {
      yjm3[j] = yn[j] + temp1 * yjm1[j];
    }
    //-------- Stage 2.--------    
    ci2 = ci1 + temp1;
    temp1 = h * fpa[mz - 1][2 - 1];
    temp2 = h * fpa[mz - 1][3 - 1];
    fun(neqn, ci2, yjm3, yjm2);

    for (j = 0; j < neqn; ++j) {
      yjm4[j] = yn[j] + temp1 * yjm1[j] + temp2 * yjm2[j];
    }
    //-------- Stage 3.--------
    ci2 = ci1 + temp1 + temp2;
    temp1 = h * fpa[mz - 1][4 - 1];
    temp2 = h * fpa[mz - 1][5 - 1];
    temp3 = h * fpa[mz - 1][6 - 1];
    fun(neqn, ci2, yjm4, yjm3);

    for (j = 0; j < neqn; ++j) {
      fnt[j] = yn[j] + temp1 * yjm1[j] + temp2 * yjm2[j] + temp3 * yjm3[j];
    }
    // -------- Stage 4.--------
    ci2 = ci1 + temp1 + temp2 + temp3;
    temp1 = h * fpb[mz - 1][1 - 1];
    temp2 = h * fpb[mz - 1][2 - 1];
    temp3 = h * fpb[mz - 1][3 - 1];
    temp4 = h * fpb[mz - 1][4 - 1];
    fun(neqn, ci2, fnt, yjm4);

    for (j = 0; j < neqn; ++j) {
      yn[j] =
        yn[j] + temp1 * yjm1[j] + temp2 * yjm2[j] + temp3 * yjm3[j] +
        temp4 * yjm4[j];
    }
    // -------- Error evaluation (embedded method of order 3).--------
    temp1 = h * fpbe[mz - 1][1 - 1] - temp1;
    temp2 = h * fpbe[mz - 1][2 - 1] - temp2;
    temp3 = h * fpbe[mz - 1][3 - 1] - temp3;
    temp4 = h * fpbe[mz - 1][4 - 1] - temp4;
    temp5 = h * fpbe[mz - 1][5 - 1];
    fun(neqn, t + h, yn, fnt);
    err = 0.00e0;
    // -------- Atol and rtol are scalar.--------
    for (j = 0; j < neqn; ++j) {
      err =
        err +
        std::pow((temp1 * yjm1[j] + temp2 * yjm2[j] + temp3 * yjm3[j] +
                  temp4 * yjm4[j] + temp5 * fnt[j]) / (atol +
                                                       std::abs(y[j]) *
                                                       rtol), 2.0);
    }

    err = std::sqrt(err / double (neqn));

         /*--deallocate--*/
    // delete [] yjm1;
    // delete [] yjm2;
    // delete [] yjm3;
    // delete [] yjm4;
    // delete [] fn;
    // delete [] fnt;
    return;
  }
// *----------------------------------------------
//    Core integrator for ROCK4.
// ----------------------------------------------*
  template < typename TRHSFun > void rock4_solve(const int neqn,
						 const double t_begin, 
						 const double *y0, 
						 const double t_end,
						 const double h0,
						 const TRHSFun  fun,      //void (*fun)(const int,const double,const double*, double*),
                                                 double *y_out,
                                                 const double atol = 1e-4,
                                                 const double rtol = 1e-3,
                                                 const bool verbose = false) {

    /*local vars */
    const double hmin =
      1.0e2 * EPSILON * std::max(std::abs(t_begin), std::abs(t_end));
    const double safety = 0.8e0;
    double rho = 0.0e0;
    double t_step = h0, t_step_old = h0;
    double t_now = t_begin;
    int mdeg = 5;
    double fac = 1.0e0;
    const double facmax = 1.5e0;
    double scale = 1.0e0;
    double err_loc = 0.0e0, err_loc_old = 0.0e0;
    double temp1 = 1.0e0, temp2 = 1.0e0;
    int nrej = 0;
    int naccp = 0;
    int nrho = 0;
    bool spectrOK = false;

    //
    double *y1 = new double[neqn];
    double *y2 = new double[neqn];
    double *work = new double[6 * neqn];

        /*--initial sptep size --*/

    t_step = std::min(t_step, t_end - t_now);

    std::memcpy(y1, y0, sizeof(double) * neqn);
    std::memcpy(y2, y0, sizeof(double) * neqn);

    spectrOK = spectralRadius(neqn, t_now, y1, fun, rho);
    t_step = std::max(0.25 / rho, t_step);
    t_step = std::min(t_step, t_end - t_now);
    t_step_old = t_step;

    while (t_now < t_end) {

      if (nrho % 5 == 0) {
        spectrOK = spectralRadius(neqn, t_now, y1, fun, rho);
      }
      nrho++;
      //calculate the eigmax: call spectral radius
      mdeg = 1 + int (std::sqrt((3.e0 + t_step * rho) / 0.353e0));

      if (mdeg > 152) {
        t_step = safety * (std::pow(152.e0, 2.0e0) * 0.353e0 - 3.e0) / rho;
        mdeg = 152;
      }
      mdeg = std::max(mdeg, 5) - 4;
      t_step = std::min(t_step, t_end - t_now);
      if (verbose) {
        std::cerr << "_rock4: mdeg-4=" << mdeg << "\t" << "t_step=" <<
          t_step << std::endl;
      }
      //rk::rock4_step(neqn,t_now,y1,t_step,fun, mdeg, y2, err_loc, atol, rtol);//perform the step
      rk::rock4_step(neqn, t_now, y1, t_step, fun, mdeg, y2, err_loc, atol, rtol, work);        //perform the step

      if (err_loc >= 1.0e0) {
        nrej++;
        t_step = safety * t_step / std::pow(err_loc, 1.0 / 5.0);
        t_step = std::max(t_step, hmin);
      } else {
        naccp++;
                  /*--determine new step size--*/
        t_now += t_step;
        for (int i = 0; i < neqn; ++i) {
          y1[i] = y2[i];
        }

        fac = facmax;

        if (naccp == 1) {
          temp2 = std::pow(err_loc, 1.0 / 5.0);
          if (safety < fac * temp2) {
            fac = safety / temp2;
          }
        } else {
          temp1 = safety * t_step * std::pow(err_loc_old, 1.0 / 5.0);
          temp2 = t_step_old * std::pow(err_loc, 2.0 / 5.0);
          if (temp1 < fac * temp2) {
            fac = temp1 / temp2;
          }
        }
        t_step_old = t_step;
        err_loc_old = err_loc;
        scale = std::max(0.1e0, fac);
        t_step *= scale;
        t_step = std::max(hmin, t_step);
      }

    }                           //while(t<t_end)

    std::memcpy(y_out, y2, sizeof(double) * neqn);

    if (verbose) {
      std::cerr << "_rock4: nrej=" << nrej << "\t" << "naccpt=" << naccp
        << std::endl;
    }
    delete [] y1;
    delete [] y2;
    delete [] work;

    return;
  }


};                              //namespace
