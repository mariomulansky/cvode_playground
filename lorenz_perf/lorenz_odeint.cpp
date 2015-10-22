#include <iostream>
#include <array>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

const int M = 100;   // repetitions

const double t_start = 0.0;
const double t_end = 100.0;
const double dt = 1.0;


typedef boost::array< double , 3 > state_type;

struct lorenz
{
    double sigma, R, b;
    int calls;

    lorenz(double sigma, double R, double b)
        : sigma(sigma), R(R), b(b), calls(0)
    {}

    void operator()( const state_type &x , state_type &dxdt , double t )
    {
        calls++;
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }

};

int main(int argc, char **argv)
{
    double x0_final = 0.0;
    lorenz rhs(sigma, R, b);
    for(int m=0; m<M; ++m)
    {
        state_type x = {{ 10.0 , 1.0 , 1.0 }}; // initial conditions
        integrate_const(make_dense_output(1E-6, 1E-6, runge_kutta_dopri5<state_type>()),
                        std::ref(rhs), x, t_start, t_end, dt);
        x0_final += x[0];
    }
    std::cout << x0_final << std::endl;
    std::cout << "rhs_calls: " << rhs.calls << std::endl;
}
