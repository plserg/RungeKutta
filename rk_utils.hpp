    

#ifndef _RK_UTILS_
#define _RK_UTILS_
/*
  basic utilities for Runge-Kutta methods
  send comments to sergey.plyasunov@gmail.com  
 */
 
#include <cmath>
#include <algorithm>
#include <cstring>
 
namespace rk{
 
const static double EPSILON = 1.0e-16;
//return the L2-norm
inline
double norm(const int neqn,const double* v){
        double sum=0.0e0;
        for(int i=0;i<neqn;++i)
                sum += v[i]*v[i];
        
        return std::sqrt(sum);
}
   
template<typename TFun>
bool
spectralRadius( const int neqn,
                const double t,
                const double* yn,
                const TFun& fun,//void (*fun)(const int,const double,const double*,double*),
                double& rho){
 
        const double epsilon = EPSILON;
        bool converged = false;
        double ynrm,vnrm,dynrm,dfunnrm;
        double lambda=0.0e0;
        double lambdaold=0.0e0;
        int i=0;
        int itmax =30;
        long nfune =0;
 
        /* --temporaly arrays-- */
        double *ev = new double[neqn];
        double *v  = new double[neqn];
        double *funv = new double[neqn];
        double *funn = new double[neqn];
        
        /*------*/      
        fun(neqn,t,yn,funn);
        nfune++;
 
        std::memcpy(v, funn, sizeof(double)*neqn );
 
        vnrm = norm(neqn,v);
        ynrm = norm(neqn,yn);
 
       if(std::abs(ynrm)>0 && std::abs(vnrm) >0 ){
                dynrm  = ynrm * std::sqrt(epsilon);
 
                for(i=0;i<neqn;++i){
                        v[i] = yn[i] + (dynrm/vnrm) * v[i];
                }      
        }else{
                if(std::abs(ynrm)>0){
                        dynrm = ynrm *std::sqrt(epsilon);
                        for(i=0;i<neqn;++i)
                          v[i]  =(1.0 + std::sqrt(epsilon))*yn[i];
                }
                if(std::abs(vnrm)>0){
                        dynrm = epsilon;
                        for(i=0;i<neqn;++i)
                          v[i] = (dynrm/vnrm) *v[i];
                }
                if(std::abs(ynrm)==0 && std::abs(vnrm) ==0 ){
                         dynrm = epsilon;
                         v[0] = dynrm;
                }
        }       
        
 
        std::memset(ev, 0, sizeof(double)*neqn );
       lambda =0.0e0;
        //iterative loop
        for(int iter=0;iter<itmax;++iter){
            
                fun(neqn, t,v,funv);
                nfune++;
                dfunnrm=0.0e0;
                for(i=0;i<neqn;++i) {
                  dfunnrm += (funv[i]-funn[i])*(funv[i]-funn[i]);
                }
                dfunnrm = std::sqrt( dfunnrm );
                lambdaold = lambda;
                lambda = dfunnrm/dynrm;
                rho = 1.1e0*lambda; 
        
                if( iter>=1 && std::abs(lambda - lambdaold) <= 0.01*std::max(lambda , lambdaold) ){
 
                        for(i=0;i<neqn;++i) {
                            ev[i] = v[i] - yn[i]; 
                        }         
                        converged = true;
                        break;
                }       
                if( std::abs(dfunnrm) >0  ){
                        for(i=0;i<neqn;++i) {         
                             v[i] = yn[i] +( funv[i]-funn[i] ) *(dynrm/dfunnrm);
                        }
                }else{
                        int index = 0 + (iter % neqn);
                        v[index] = yn[index] - (v[index]-yn[index]);
                }
        }
        //clean up 
        
        delete [] funn;
        delete [] funv;
        delete [] v;
        delete [] ev;
 
        return converged;
} 
};//namespace 
#endif
