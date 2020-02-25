//
// example of 2D RKC solver 
//
#include <cstdio>
#include <cmath>
#include <memory>
#include <vector>
#include <initializer_list>
#include "rkc.hpp"


class GridFunction
{
    public:
        GridFunction(std::initializer_list<int> szs): 
        strides_(szs),
        x_grid_(std::make_unique<double[]>(strides_[0])),
        y_grid_(std::make_unique<double[]>(strides_[1]))
        {
            size_ = 1;
            for(int i=0;i<strides_.size(); ++i)  
                size_ *=strides_[i];
            const double x_min = 0.000;
            const double x_max = 1.000;
            const double y_min = 0.000;
            const double y_max = 0.999;
            double dx = (x_max - x_min)/(strides_[0] -1);
            double dy = (y_max - y_min)/(strides_[1] -1);

            for(int i=0; i<strides_[0];++i)
            {
                x_grid_[i] = x_min + dx * i;
            }
            for(int i=0;i<strides_[1]; ++i)
            {
                y_grid_[i] = y_min +dy *i;
            }

        }
        int size() const noexcept { return size_; }

        void operator()(const int neqn, const double t, const  double* u, double* dudt) const noexcept 
         {
        
            for(int i=0;i < neqn; ++i)
            {
                dudt[i] = - 1.0 * u[i]; 
            }
        }

        double a_x(const double t, const double x, const double y) const noexcept
        {
            return 0.00;
        }
        
        double a_y(const double t, const double x, const double y) const noexcept
        {
            return 0.00;
        }

        double b_xx(const double t, const double x, const double y) const noexcept
        {
            return 0.00;
        }
        
        double b_xy(const double t, const double x, const double y) const noexcept
        {
            return 0.00;
        }

        double b_yy(const double t, const double x, const double y) const noexcept
        {
            return 0.00;
        }

        ~GridFunction()
        {
                std::printf("GridFunction: destructor called\n");
        }

    private:
        std::vector<int> strides_;
        std::unique_ptr<double[]> x_grid_;
        std::unique_ptr<double[]> y_grid_;
        int size_;
};

int main(int argc,char **argv)
{
   const int nx = 100; 
   const int ny = 10;
   GridFunction gridf( {nx,ny} );
   std::printf("size of the grid: %d\n", gridf.size());

   return 0;
}
