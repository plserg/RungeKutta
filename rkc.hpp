/*
 Stabilized Runge-Kuuta (2nd order)
 based on 
 "explicit solver for parabolic PDEs"
  by Sommeijer, Shampine, Verwer, J. Comp. Appl. Math., 1997  

  send commnets to sergey.plyasunov@gmail.com
*/

#pragma once

#include <iostream>
#include <cstring>
#include "rk_utils.hpp"
 
namespace rk{

//second order RKC method  
template<typename TRHSFun>
void
rkc_step(const int     neqn,
         const double  t,
         const double  *yn,
         const double  t_step,
         const TRHSFun fun,
         int nstages,   
         double *y_out,  
         double& err, 
         const double atol,
         const double rtol,
	 const double eps, /*=2.0/13.0*/
	 double *work){
 
  double *y= work;//new double[neqn];
  double *yjm1= work+neqn;//new double[neqn];
  double *yjm2= work+2*neqn;//new double[neqn];
  double *rhs=  work+3*neqn;//new double[neqn];
  double *dye_0= work+4*neqn;//new double[neqn];
  double *dye_jm1= work+5*neqn;//new double[neqn];
  
  double nu,mu,mus,mus1;
  double w0,w1; 
  double ajm1,bj,bjm1,bjm2,dzj,dzjm1,dzjm2,d2zj,d2zjm1,d2zjm2,thj,thjm1,thjm2;
  double zj,zjm1,zjm2,gamma;
  double arg,temp1,temp2;
  long   NfunE=0;
  double abs_t_step = std::abs(t_step);

  /*-compute basic arguments--*/
  //w0    = 1.0e0 + 2.0e0/(13.0e0*nstages*nstages);
  w0    = 1.0 + eps/(nstages*nstages);
  temp1 = w0 * w0 - 1.0;
  temp2 = std::sqrt(temp1);
  arg = nstages * std::log(w0 + temp2);
  w1 = std::sinh(arg) * temp1/( std::cosh(arg) * nstages * temp2 - w0 * std::sinh(arg) );
  bjm1 = 1.0/w0;
  bjm2 = 1.0/(4.0*w0*w0);
  /*---call y0----*/
  std::memcpy(y, yn, sizeof(double)*neqn );
  //
  fun(neqn,t,yn,dye_0);
  NfunE++;
  //
  std::memcpy(yjm2, yn, sizeof(double)*neqn );
  mus = w1 * bjm1;
  mus1 = mus;
  //
  for(int i=0;i<neqn;++i){
    rhs[i] = y[i] + abs_t_step * mus1 * dye_0[i];
  }
  //
  std::memcpy(yjm1, rhs, sizeof(double)*neqn );
  std::memcpy(y, rhs, sizeof(double)*neqn );
        
  thjm2 = 0.0;
  thjm1 = mus;
  zjm1  = w0;
  zjm2  = 1.0;
  dzjm1 = 1.0;
  dzjm2 = 0.0;
  d2zjm1= 0.0;
  d2zjm2= 0.0;
 
  for(int s=2;s<=nstages;++s){
 
        zj  =  2.0 * w0 * zjm1 - zjm2;
        dzj =  2.0 * w0 * dzjm1 - dzjm2 + 2 * zjm1;
        d2zj = 2.0 * w0 * d2zjm1 - d2zjm2  + 4 * dzjm1;
        bj = d2zj/(dzj*dzj);
        ajm1 = 1.0 - zjm1*bjm1;
        mu = 2.0 * w0 *bj/bjm1;
        nu = -bj/bjm2;
        mus = mu * w1/w0;
        gamma = -ajm1*mus;
        
        fun(neqn,t + t_step * thjm1, yjm1, dye_jm1);   
        NfunE++;

        for(int i=0;i<neqn;++i){
           rhs[i] = (1.0 - mu -nu) * yn[i] + mu * yjm1[i] + nu * yjm2[i]  + abs_t_step * mus * dye_jm1[i] + abs_t_step * gamma * dye_0[i]; 
        }

        thj = mu *thjm1 + nu * thjm2 + mus*(1.0 - ajm1);
        /* explcit step */
        std::memcpy(y, rhs, sizeof(double)*neqn );
        
        if(s<nstages){
 
                std::memcpy(yjm2, yjm1, sizeof(double)*neqn );
                std::memcpy(yjm1, y, sizeof(double)*neqn );
                thjm2 = thjm1;
                thjm1 = thj;
                bjm2 = bjm1;
                bjm1 = bj;
                zjm2 = zjm1;
                zjm1 = zj;
                dzjm2 = dzjm1;
                dzjm1 = dzj;
                d2zjm2 = d2zjm1;
                d2zjm1 = d2zj;  
        }       
  }
  //assign the y_out
  std::memcpy(y_out, y, sizeof(double)*neqn );
  //calcaulte local erorr:
  
  fun(neqn,t+t_step,y_out,dye_jm1);
  err =0.0e0;
 
  for(int i=0;i<neqn;++i){
    err = err + std::pow( 0.5e0 * abs_t_step *(dye_jm1[i]-dye_0[i])/(atol + rtol*std::abs(yn[i]) ), 2.0e0 );
  }
  err = std::sqrt(err/neqn);
  //clean up
 /* delete [] y;
  delete [] yjm1;
  delete [] yjm2;
  delete [] rhs;
  delete [] dye_0;
  delete [] dye_jm1;
 */ 
}
 
template<typename TRHSFun>
void 
rkc_solve(const int neqn,
         const double t_begin,
         const double *y0,
         const double t_end,
         const double h0,
         const TRHSFun& fun,// void (*fun)(const int,const double, const double*, double*),
         double *y_out,
         const  double atol=1e-4,
         const  double rtol=1e-3,
         const  bool verbose=false,
	 const  double eps=2.0/13.0){
  
        const int dir = (t_end>=t_begin)? +1:-1;
        const int  max_num_stages=30;
        double rho=0;
        double t_step,t_step_old;
        double t_now = t_begin;     
        int num_stages = 2;
        double fac = 10.0e0;
        //double scale =1.0e0;
        double err_loc =0.0e0,err_loc_old=0.0e0;
        double temp1,temp2;
        const double hmin = 1.0e2 *EPSILON * std::max( std::abs(t_begin),std::abs(t_end) );        
        const double safety =0.8e0;
        long nrej =0;
        long naccpt=0;
        long nrho =0;
 
        double *y1 = new double[neqn];
        double *y2 = new double[neqn];
	double *work = new double[6*neqn];
 
       /*--estimate the initial  spectral  radius--*/
        bool spektrOk = spectralRadius( neqn, t_begin, y0, fun, rho);
        rho = std::max(rho,EPSILON);
       /*--initial sptep size --*/
        t_step =  std::max( 0.25/rho, h0);
        t_step = std::min( t_step , std::abs(t_end - t_now));
        /*-one stage step -*/
        std::memcpy(y1, y0, sizeof(double)*neqn );
        std::memcpy(y2, y0, sizeof(double)*neqn );
 
	while((t_end-t_now)*dir>1.0e-7){
       
          nrho++;      
         if(nrho % 5 ==0){
              spektrOk = spectralRadius(neqn, t_now, y1, fun, rho);
         }
 
	 const double inv_beta = 1.0/( (2.0/3.0) * (1.0 - 2.0/15.0 *eps) );//~1.54e0
         //num_stages = 1  + (int)std::sqrt(1.0e0 + 1.54e0 * t_step *rho );
         num_stages = 1  + (int)std::sqrt(1.0e0 + inv_beta * t_step *rho );  
         
	if(num_stages > max_num_stages){  
                        num_stages = max_num_stages;
                        t_step = (num_stages*num_stages - 1.0e0 )/(inv_beta*rho);
          }
 
	 t_step = std::min( t_step , std::abs(t_end - t_now));
          if(verbose){
                        std::cerr<< "_rkc: spectral radius="<< rho << ", num_stages="<< num_stages <<std::endl;
          } 
 
           /*rk::rkc_step(neqn,t_now,y1,dir*t_step,fun,num_stages,y2,err_loc,atol,rtol,eps);*/
	   rk::rkc_step(neqn,t_now,y1,dir*t_step,fun,num_stages,y2,err_loc,atol,rtol,eps,work);
 
           if(err_loc >= 1.0e0){
 
                        nrej++;
                        t_step = safety * t_step/std::pow( err_loc, 0.33e0 );
                        t_step = std::max(t_step , hmin );              
                                
           }else{
                       naccpt++;
                        /*--new step size--*/
                        t_now += t_step*dir;
                        std::memcpy(y1, y2, sizeof(double)*neqn );
                        fac = 10.0e0;
         
                        if(naccpt==1){
                           temp2 = std::pow( err_loc , 0.33e0 );
                           if(safety < fac * temp2 ) {
                                        fac = safety/temp2;     
                        }            
                        }else{
                          temp1 = safety  * t_step * std::pow(err_loc_old, 0.33e0 );
                          temp2 = t_step_old * std::pow(err_loc, 0.66e0);
                          if(temp1 < fac *temp2){
                                        fac = temp1/temp2;
                          }               
                        }
 
                        t_step_old = t_step; 
                        err_loc_old = err_loc;
                        t_step *= std::max( 0.1e0, fac );
                        t_step = std::max( hmin, t_step );
                
                     
         }
        
        }//while 
 
       std::memcpy(y_out, y2, sizeof(double)*neqn );
       
       delete [] work; 
       delete [] y1;
       delete [] y2;
}
 
template<typename TRHSFun,typename TBoundFun>
void
rkc_lcp_step(const int  neqn,
         const double  t,
         const double *yn,
        const double *lmN,
         const double  t_step,
         const TRHSFun& fun,//void (*fun)(const int,const  double,const double*, double*),
         const TBoundFun& bound_fun,//void (*bound_fun)(const int,const double, double*),/*provides a bound vector*/
         int nstages,   
         double *y_out, 
         double *lm_OUT,
         double& err, 
         const double atol=1e-4,
         const double rtol=1e-3,
	 const double eps=2.0/13.0){
 
  double *y= new double[neqn];
  double *yjm1= new double[neqn];
  double *yjm2= new double[neqn];
  double *rhs=new double[neqn];
  double *dye_0=new double[neqn];
  double *dye_jm1= new double[neqn];
  double *g = new double[neqn];
  
  double nu,mu,nuS,mus,mus1;
  double w0,w1; 
  double ajm1,bj,bjm1,bjm2,dzj,dzjm1,dzjm2,d2zj,d2zjm1,d2zjm2,thj,thjm1,thjm2;
  double zj,zjm1,zjm2,gamma;
  double arg,temp1,temp2;
  long   NfunE=0;
  double abs_t_step = std::abs(t_step);
  /*-compute basic arguments--*/
  w0    = 1.0e0 + eps/(nstages*nstages);
  temp1 = w0 * w0 - 1.0e0;
  temp2 = std::sqrt(temp1);
  arg = nstages * std::log(w0 + temp2);
  w1 = std::sinh(arg) * temp1/( std::cosh(arg) * nstages * temp2 - w0 * std::sinh(arg) );
  bjm1 = 1.0e0/w0;
  bjm2 = 1.0e0/(4.0e0*w0*w0);
  /*---call y0----*/
  std::memcpy(y, yn, sizeof(double)*neqn );
  //
  fun(neqn,t,yn,dye_0);
  NfunE++;
  //
  std::memcpy(yjm2, yn, sizeof(double)*neqn );
  mus = w1 * bjm1;
  mus1 = mus;
 
  for(int i=0;i<neqn;++i){
    rhs[i] = y[i] + abs_t_step * mus1 * (dye_0[i] + lmN[i]);
  }
  std::memcpy(yjm1, rhs, sizeof(double)*neqn );
  std::memcpy(y, rhs, sizeof(double)*neqn );
        
  thjm2 = 0.0e0;
  thjm1 = mus;
  zjm1  = w0;
  zjm2  = 1.0e0;
  dzjm1 = 1.0e0;
  dzjm2 = 0.0e0;
  d2zjm1= 0.0e0;
  d2zjm2= 0.0e0;
 
  for(int s=2;s<=nstages;++s){
 
        zj  =  2.0e0 * w0 * zjm1 - zjm2;
        dzj =  2.0e0 * w0 * dzjm1 - dzjm2 + 2 * zjm1;
        d2zj = 2.0e0 * w0 * d2zjm1 - d2zjm2  + 4 * dzjm1;
        bj = d2zj/(dzj*dzj);
        ajm1 = 1.0e0 - zjm1*bjm1;
        mu = 2.0e0 * w0 *bj/bjm1;
        nu = -bj/bjm2;
        mus = mu * w1/w0;
        gamma = -ajm1*mus;
        
        fun(neqn,t + t_step * thjm1, yjm1, dye_jm1);   
        NfunE++;
  
       for(int i=0;i<neqn;++i){
            rhs[i] = (1.0e0 - mu -nu) * yn[i] + mu * yjm1[i] + nu * yjm2[i]  + abs_t_step * mus * (dye_jm1[i] + lmN[i]) + abs_t_step * gamma * (dye_0[i] + lmN[i] ); 
       }

        thj = mu *thjm1 + nu * thjm2 + mus*(1.0e0 - ajm1);
        /* totally explcit*/
        std::memcpy(y, rhs, sizeof(double)*neqn );
        
        if(s<nstages){
                std::memcpy(yjm2, yjm1, sizeof(double)*neqn );
                std::memcpy(yjm1, y, sizeof(double)*neqn );
                thjm2 = thjm1;
                thjm1 = thj;
                bjm2 = bjm1;
                bjm1 = bj;
                zjm2 = zjm1;
                zjm1 = zj;
                dzjm2 = dzjm1;
                dzjm1 = dzj;
                d2zjm2 = d2zjm1;
                d2zjm1 = d2zj;  
        }       
  }
 
 
  bound_fun(neqn, t, g);//what is the constraint at the current step?
 
  for(int i=0;i<neqn;++i){
    lm_OUT[i] = std::max((g[i]-y[i])/abs_t_step + lmN[i],0.0e0);
    y[i] = y[i] + ( lm_OUT[i] - lmN[i] )*abs_t_step;
  }
  //assign the y_out
  std::memcpy(y_out, y, sizeof(double)*neqn );
  //calcaulte local erorr:
  
  fun(neqn,t+t_step,y_out,dye_jm1);
  err = 0.0e0;
 
  for(int i=0;i<neqn;++i){
    err = err + std::pow( 0.5e0*abs_t_step*(dye_jm1[i]-dye_0[i])/(atol + rtol*std::abs(yn[i]) ), 2.0e0 );
  }
  err = std::sqrt(err/neqn);
  //clean up
  delete [] y;
  delete [] yjm1;
  delete [] yjm2;
  delete [] rhs;
  delete [] dye_0;
  delete [] dye_jm1;
  delete [] g;
 
}
 
template<typename TRHSFun,typename TBoundFun>
void 
rkc_lcp_solve(const int neqn,
         const double t_begin,
         const double *y0,
         const double t_end,
         const double h0,
         const TRHSFun& fun,//void (*fun)(const int,const double, const double*, double*),
         const TBoundFun& bound_fun,//void (*bound_fun)(const int,const double, double*),
         double *y_out,
         const  double atol=1e-4,
         const  double rtol=1e-3,
         const  bool verbose=false,
	 const  double eps=2.0/13.0){
 
        const int dir = (t_end>=t_begin)? +1:-1;
        const int  max_num_stages=30;
        double rho=0.0e0;
        double t_step,t_step_old;
        double t_now = t_begin;     
        int num_stages = 2;
        double fac = 10.0e0;
        double scale =1.0e0;
        double err_loc =0.0e0,err_loc_old=0.0e0;
        double temp1,temp2;
        const double hmin = 1.0e2 *EPSILON * std::max( std::abs(t_begin),std::abs(t_end) );        
        const double safety =0.8e0;
        long nrej =0;
        long naccpt=0;
        long nrho =0;
 
        double *y1 = new double[neqn];
        double *y2 = new double[neqn];
        double *lm1= new double[neqn];
        double *lm2= new double[neqn];
 
        /*--estimate the initial  spectral  radius--*/
        bool spektrOk= spectralRadius( neqn, t_begin, y0, fun, rho);
        rho = std::max(rho,EPSILON);
        /*--initial sptep size --*/
        t_step =  std::max( 0.25/rho, h0);
        t_step = std::min( t_step , std::abs(t_end - t_now));
        /*-one stage step -*/
 
        std::memcpy(y1, y0, sizeof(double)*neqn );
        std::memcpy(y2, y0, sizeof(double)*neqn );
        std::memset(lm1, 0, sizeof(double)*neqn );
        std::memset(lm2, 0, sizeof(double)*neqn );
 
	while((t_end-t_now)*dir>1.0e-7){
       
         nrho++;      
         if(nrho % 20 ==0){
              spektrOk = spectralRadius(neqn, t_now, y1, fun, rho);
         }
 	 
 	 const double inv_beta = 1.0/( (2.0/3.0) *(1.0 - 2.0/15.0*eps));//1.54 
         //num_stages = 1 + (int)std::sqrt(1.0e0 + 1.54e0 * t_step *rho );
         num_stages = 1  + (int)std::sqrt(1.0e0 + inv_beta * t_step *rho );
         
	if(num_stages > max_num_stages){  
                        num_stages = max_num_stages;
                        t_step = (num_stages*num_stages - 1.0e0 )/(inv_beta*rho);
          }
 
         t_step = std::min( t_step , std::abs(t_end - t_now));
 
          if(verbose){
                        std::cerr<< "_rkc_lcp: spectral radius="<< rho << ", num_stages="<< num_stages <<std::endl;
          } 
 
          rkc_lcp_step(neqn,t_now,y1,lm1,dir*t_step,fun,bound_fun,num_stages,y2,lm2,err_loc,atol,rtol,eps);
 
           if(err_loc >= 1.0e0){
 
                        nrej++;
                        t_step = safety * t_step/std::pow( err_loc, 0.33e0 );
                        t_step = std::max(t_step , hmin );              
                                
           }else{
                       naccpt++;
                        /*--new step size--*/
                        t_now += t_step*dir;
 
                        std::memcpy(y1, y2, sizeof(double)*neqn );
                        std::memcpy(lm1, lm2, sizeof(double)*neqn );
 
                        fac = 10.0e0;
         
                        if(naccpt==1){
                         temp2 = std::pow( err_loc , 0.33e0 );
                           if(safety < fac * temp2 ) {
                                        fac = safety/temp2;     
                        }            
                        }else{
                          temp1 = safety  * t_step * std::pow(err_loc_old, 0.33e0 );
                          temp2 = t_step_old * std::pow(err_loc, 0.66e0);
                          if(temp1 < fac *temp2){
                                        fac = temp1/temp2;
                          }               
                        }
 
                        t_step_old = t_step; 
                        err_loc_old = err_loc;
                        t_step *= std::max( 0.1e0, fac );
                        t_step = std::max( hmin, t_step );
                
                     
         }
        
        }//while 
 
        std::memcpy(y_out, y2, sizeof(double)*neqn );
 
       delete [] y1;
       delete [] y2;
       delete [] lm1;
       delete [] lm2;
}
 
};//namespace 
