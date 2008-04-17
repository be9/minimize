#include "nm.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric::ublas;

namespace minimize { namespace nelder_mead {

// Nelder-Mead parameters
const double RHO   = 1;
const double CHI   = 2;
const double PSI   = 0.5;
const double SIGMA = 0.5;

// sort so simplex column 0 has the lowest function value
static void sort_simplex(values_vec &FV, simplex_mat &V);
static int change_simplex(values_vec &FV, simplex_mat &V, cost_function func, report_details &det);

struct BreakException {};

void run(const params_vector &x0, cost_function func, output &out, report_func reporter,
        int maxiter, int maxfunevals, double tolx, double tolfun)
{
    int n = x0.size();

    if (maxiter < 0)
        maxiter *= -n;

    if (maxfunevals < 0)
        maxfunevals *= -n;
  
    // Set up a simplex near the initial guess.
    simplex_mat V(n, n+1);
    values_vec FV(n+1);

    // Place input guess in the simplex! (credit L.Pfeffer at Stanford)
    column(V, 0) = x0;
    FV[0] = func(x0);
  
    int func_evals = 1;
    int itercount = 0;

    // Initial simplex setup continues later
    try {
        if (reporter && 
                (reporter(INIT, NOTHING, itercount, func_evals, FV, V) || 
                 reporter(ITER, NOTHING, itercount, func_evals, FV, V)))
            throw BreakException();

        // Continue setting up the initial simplex.
        // Following improvement suggested by L.Pfeffer at Stanford

        const double USUAL_DELTA = 0.05;          // 5 percent deltas for non-zero terms
        const double ZERO_TERM_DELTA = 0.00025;   // Even smaller delta for zero elements of x

        for (int j = 0; j < n; ++j) {
            params_vector y = x0;

            if (y[j] != 0.0)
                y[j] *= 1.0 + USUAL_DELTA;
            else
                y[j] = ZERO_TERM_DELTA;

            column(V, j+1) = y;
            FV[j+1] = func(y);
        }

        sort_simplex(FV, V);

        ++itercount;
        func_evals = n+1;

        if (reporter && reporter(ITER, INITIAL_SIMPLEX, itercount, func_evals, FV, V))
            throw BreakException();
  
        // Main algorithm
        // Iterate until the diameter of the simplex is less than tolx
        // AND the function values differ from the min by less than tolf,
        // or the max function evaluations are exceeded. (Cannot use OR instead of
        //  AND.)
        while (func_evals < maxfunevals && itercount < maxiter) {
            double FVdiff = norm_inf(subrange(FV, 1, n) - scalar_vector<double>(n-1, FV[0]));
            
            simplex_mat V2 = subrange(V, 0,n, 1,n+1);

            for (int j = 0; j < n; ++j)
                column(V2, j) -= column(V, 0);

            double Vdiff = norm_inf(V2);

            if (FVdiff <= tolfun && Vdiff <= tolx)
                break;

            report_details det;
            func_evals += change_simplex(FV, V, func, det);

            sort_simplex(FV, V);

            ++itercount;

            if (reporter && reporter(ITER, det, itercount, func_evals, FV, V))
                throw BreakException();
    
        }

        if (reporter && reporter(DONE, NOTHING, itercount, func_evals, FV, V)) 
            ;


        if (func_evals >= maxfunevals)
            out.code = MAX_FUNC_EVALS;
        else if (itercount >= maxiter)
            out.code = MAX_ITER;
        else
            out.code = SUCCESS;
    }

    catch (BreakException) {
        out.code = ABORT;
    }
    
    out.last_x           = column(V, 0);
    out.last_f           = FV[0];
    out.iterations       = itercount;
    out.func_evaluations = func_evals;
}

static void sort_simplex(values_vec &FV, simplex_mat &V)
{
    //indexes = range(shape(simplex)[1])
    //indexes.sort(key=lambda i: values[i])
    //return simplex[:, indexes], values[indexes]
}

static int change_simplex(values_vec &FV, simplex_mat &V, cost_function func, report_details &det)
{
    return 0;
}

} }
