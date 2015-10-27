#include <stdio.h>

#include <cvode/cvode.h>             /* main integrator header file */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_types.h> /* definition of realtype */

#define SIGMA RCONST(10.0)
#define R RCONST(28.0)
#define B RCONST(8.0/3.0)

#define N 3       // dimensionality
#define M 100     // repetitions

int rhs_calls = 0;

int lorenz(realtype t, N_Vector y, N_Vector y_dot, void *params)
{
    rhs_calls++;

    // lorenz equations
    NV_Ith_S(y_dot, 0) = SIGMA * ( NV_Ith_S(y, 1) - NV_Ith_S(y, 0) );
    NV_Ith_S(y_dot, 1) = R * NV_Ith_S(y, 0) - NV_Ith_S(y, 1) - NV_Ith_S(y, 0) * NV_Ith_S(y, 2);
    NV_Ith_S(y_dot, 2) = -B * NV_Ith_S(y, 2) + NV_Ith_S(y, 0) * NV_Ith_S(y, 1);

    return 0;
}

int main()
{
    const realtype reltol = 1E-6;
    const realtype abstol = 1E-6;
    const realtype t_start = 0.0;
    const realtype t_end = 100.0;
    const realtype dt = 1.0;
    const int steps = 100;
    int m, t_iter;
    realtype t;
    realtype y0_tot = 0;
    N_Vector y;

    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);

    y = N_VNew_Serial(N);

    NV_Ith_S(y, 0) = 10.0;
    NV_Ith_S(y, 1) = 1.0;
    NV_Ith_S(y, 2) = 1.0;

    for(m=0; m<M; m++)
    {
        if( m == 0 )
        {
            CVodeInit(cvode_mem, lorenz, t_start, y);
            CVodeSStolerances(cvode_mem, reltol, abstol);
        } else
        {
            CVodeReInit(cvode_mem, t_start, y);
        }

        realtype t_out = t_start+dt;
        // main step iteration
        for(t_iter=0; t_iter < steps; t_iter++, t_out += dt)
        {
            CVode(cvode_mem, t_out, y, &t, CV_NORMAL);
        }

        y0_tot += NV_Ith_S(y, 0);
    }
    CVodeFree(&cvode_mem);
    N_VDestroy_Serial(y);

    printf("%.5f\n", y0_tot);
    printf("rhs_calls: %d\n", rhs_calls);

    return 0;
}
