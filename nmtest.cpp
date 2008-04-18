#include <assert.h>
#include <iostream>
#include <iomanip>
#include "nm.h"
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace minimize;
using namespace minimize::nelder_mead;

#define WOODS

#ifndef WOODS
double func(const params_vector &params)
{
    assert(params.size() == 2);

    double x = params[0], y = params[1];

    return (1.0-x)*(1.0-x) + 100 * (y-x*x)* (y-x*x);
}

params_vector initial()
{
    params_vector x0(2);

    x0[0] = 100.0;
    x0[1] = 100.0;

    return x0;
}
#else
/* Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566) */

double func(const params_vector &x)
{
	double s1, s2, s3, t1, t2, t3, t4, t5;
	s1 = x[1] - x[0]*x[0];
	s2 = 1 - x[0];
	s3 = x[1] - 1;
	t1 = x[3] - x[2]*x[2];
	t2 = 1 - x[2];
	t3 = x[3] - 1;
	t4 = s3 + t3;
	t5 = s3 - t3;
	return 100*(s1*s1) + s2*s2 + 90*(t1*t1) + t2*t2 + 10*(t4*t4) + t5*t5/10.;
}

params_vector initial()
{
    params_vector x0(4);
	
    x0[0] = -3;
	x0[1] = -1;
	x0[2] = -3;
	x0[3] = -1;

    return x0;
}

#endif

bool reporter(report_stage stage, report_details details, int iter, int func_evals, 
        const values_vec &values, const simplex_mat &simplex)
{
    static const char *text_stages[] = {
        "init", "iter", "done"
    };
    
    static const char *text_details[] = {
        "nothing", "initial simplex", 
        "reflect", "expand", 
        "contract inside", "contract outside", "shrink"
    };

    cout << "===================================================" << endl;
    cout << "Iteration " << iter << " (" << func_evals << " func. evals): ";
    cout << text_stages[stage] << " / " << text_details[details] << endl;
   
    cout << "Simplex: " << simplex << endl;
    cout << "Values:  " << values << endl;
    cout << "---------------------------------------------------" << endl;
                                           
    return false; 
}

int main()
{
    output out;

    run(initial(), func, out, reporter);

    static const char *return_codes[] = {
        "SUCCESS", "ABORT", "MAX_FUNC_EVALS", "MAX_ITER"
    };

    cout << ">>>>>>>>>>>> Nelder-Mead exit <<<<<<<<<<<<" << endl;
    cout << "Return code: " << return_codes[out.code] << endl;
    cout << "Iterations:  " << out.iterations << endl;
    cout << "Func. evals: " << out.func_evaluations << endl;
    cout << "Reached min: " << out.last_f << endl;
    cout << "Params vec:  " << out.last_x << endl;

    return 0;
}
