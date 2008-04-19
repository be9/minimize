#ifndef __nm_h
#define __nm_h

#include "minimize.h"
#include <boost/numeric/ublas/matrix.hpp>

namespace minimize { namespace nelder_mead {

// function return codes
enum return_code {
    SUCCESS,
    ABORT,
    MAX_FUNC_EVALS,
    MAX_ITER
};

struct output {
    return_code   code;
    params_vector last_x;
    double        last_f;
    int           iterations;
    int           func_evaluations;
};

// reporter function
enum report_stage {
    INIT, ITER, DONE
};

enum report_details {
    NOTHING,
    INITIAL_SIMPLEX,
    REFLECT,
    EXPAND,
    CONTRACT_INSIDE,
    CONTRACT_OUTSIDE,
    SHRINK
};

typedef boost::numeric::ublas::vector<double> values_vec;
typedef boost::numeric::ublas::matrix<double> simplex_mat;

// if report function returns true, process is aborted
typedef bool (*report_func)(report_stage stage, report_details details, 
        int iter, int func_evals, 
        const values_vec &values, const simplex_mat &simplex);

// main function
void run(
        const params_vector &x0,                // initial guess
        cost_function       func,               // cost function
        void                *farg,				// cost function arg
        output              &out,               // output results
        report_func         reporter = 0,       // reporter function
        int                 maxiter = -200,     // max. iterations (mult. by params count if negative)
        int                 maxfunevals = -200, // max. func. evals (mult. by params count if negative)
        double              tolx = 1e-4,        // simplex diameter tolerance
        double              tolfun = 1e-4);     // function values tolerance
} }

#endif
