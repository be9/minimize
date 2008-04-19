#ifndef __hj_h
#define __hj_h

#include "minimize.h"

namespace minimize { namespace hooke_jeeves {

// function return codes
enum return_code {
    SUCCESS,
    ABORT,
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

// if report function returns true, process is aborted
typedef bool (*report_func)(report_stage stage, int iter, int func_evals, 
		const params_vector &last_x, double last_f);

// main function
void run(
        const params_vector &x0,                // initial guess
        cost_function       func,               // cost function
        output              &out,               // output results
        report_func         reporter = 0,       // reporter function
        int                 maxiter = 5000,   	// max. iterations
        double              rho = 0.5,       	// simplex diameter tolerance
        double              epsilon = 1e-6);	// function values tolerance
} }

#endif
