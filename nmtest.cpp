#include <assert.h>
#include <iostream>
#include <iomanip>
#include "nm.h"
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace minimize;
using namespace minimize::nelder_mead;

double rosenbrocks(const params_vector &params)
{
    assert(params.size() == 2);

    double x = params[0], y = params[1];

    return (1.0-x)*(1.0-x) + 100 * (y-x*x)* (y-x*x);
}

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
    params_vector x0(2);

    x0[0] = 100.0;
    x0[1] = 100.0;
    
    output out;

    run(x0, rosenbrocks, out, reporter);

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
