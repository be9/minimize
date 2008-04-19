#include "nm.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <vector>

using namespace boost::numeric::ublas;

namespace minimize { namespace nelder_mead {

// Nelder-Mead parameters
const double RHO   = 1;
const double CHI   = 2;
const double PSI   = 0.5;
const double SIGMA = 0.5;

// sort so simplex column 0 has the lowest function value
static void sort_simplex(values_vec &FV, simplex_mat &V);
static int change_simplex(values_vec &FV, simplex_mat &V, 
		cost_function func, void *farg, report_details &det);

struct BreakException {};

void run(const params_vector &x0, cost_function func, void *farg, 
		output &out, report_func reporter,
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
    FV[0] = func(x0, farg);
  
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
            FV[j+1] = func(y, farg);
        }

        sort_simplex(FV, V);

        ++itercount;
        func_evals = n+1;

        if (reporter && reporter(ITER, INITIAL_SIMPLEX, itercount, func_evals, FV, V))
            throw BreakException();
  
        // Main algorithm
        // Iterate until the diameter of the simplex is less than tolx AND the
        // function values differ from the min by less than tolf, or the max
        // function evaluations are exceeded. (Cannot use OR instead of AND.)
        while (func_evals < maxfunevals && itercount < maxiter) {
            double FVdiff = norm_inf(subrange(FV, 1, n) - scalar_vector<double>(n-1, FV[0]));
            
            simplex_mat V2 = subrange(V, 0,n, 1,n+1);

            for (int j = 0; j < n; ++j)
                column(V2, j) -= column(V, 0);

            double Vdiff = norm_inf(V2);

            if (FVdiff <= tolfun && Vdiff <= tolx)
                break;

            report_details det;
            func_evals += change_simplex(FV, V, func, farg, det);

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

class SortByCostFunctionValue {
public:
	SortByCostFunctionValue(const values_vec &_FV) : FV(_FV) {}

	bool operator()(int a, int b) {
		return FV[a] < FV[b];
	}
private:
	const values_vec &FV;
};

// sort so simplex column 0 has the lowest function value
static void sort_simplex(values_vec &FV, simplex_mat &V)
{
	values_vec FVcopy = FV;
	simplex_mat Vcopy = V;
	std::vector<int> indexes(FV.size());
	int j;

	for (j = 0; j < indexes.size(); ++j)
		indexes[j] = j;

	sort(indexes.begin(), indexes.end(), SortByCostFunctionValue(FV));

	for (j = 0; j < indexes.size(); ++j) {
		FV[j] = FVcopy[indexes[j]];
		column(V, j) = column(Vcopy, indexes[j]);
	}
}

#define REPLACE_AND_RETURN(newx, newf, detail) \
	do { \
		column(V, n) = newx; \
		FV[n] = newf; \
		det = detail; \
		return func_evals; \
	} while(0)

#define SHRINK_AND_RETURN(blah) \
	do { \
		func_evals += perform_shrink(FV, V, func, farg); \
		det = SHRINK; \
		return func_evals; \
	} while(0)

static int perform_shrink(values_vec &FV, simplex_mat &V, cost_function func, void *farg);

static int change_simplex(values_vec &FV, simplex_mat &V, 
		cost_function func, void *farg, report_details &det)
{
    int func_evals = 0;
	int j, n = FV.size()-1;
	
  	// xbar = average of the n (NOT n+1) best points
	params_vector xbar = zero_vector<double>(n), worst = column(V, n);

	for (j = 0; j < n; ++j)
		xbar += column(V, j);

	xbar /= double(n);
 
	// Compute the reflection point
 	params_vector xr = (1 + RHO)*xbar - RHO*worst;
	double fxr = func(xr, farg); ++func_evals;
	
  	/***** Check the reflection point against our current best *****/
  	if (fxr < FV[0]) {
    	// Calculate the expansion point
    	params_vector xe = (1.0 + RHO*CHI)*xbar - RHO*CHI*worst;
	    double fxe = func(xe, farg); ++func_evals;

	    if (fxe < fxr)
			REPLACE_AND_RETURN(xe, fxe, EXPAND);
		else
			REPLACE_AND_RETURN(xr, fxr, REFLECT);
	}
  
	/***** Continuing, the reflection point is worse than the current best *****/
	if (fxr < FV[n-1])
	    /* ... but at least it's better than the second worst point from the end */
		REPLACE_AND_RETURN(xr, fxr, REFLECT);
  
	/***** Worse than the second worse point from the end *****/

  	// Perform contraction
	if (fxr < FV[n]) {
		/***** better than the worst point we had *****/
    	
		// Perform an outside contraction
		params_vector xc = (1 + PSI*RHO)*xbar - PSI*RHO*worst;
		double fxc = func(xc, farg); ++func_evals;
          
		if (fxc <= fxr)
			REPLACE_AND_RETURN(xc, fxc, CONTRACT_OUTSIDE);
		else
			SHRINK_AND_RETURN();
	}

	/***** Reflection bad, really bad, even worse than the worst one *****/

	// Perform an inside contraction
	params_vector xcc = (1-PSI)*xbar + PSI*worst;
	double fxcc = func(xcc, farg); ++func_evals;

	if (fxcc < FV[n])
		REPLACE_AND_RETURN(xcc, fxcc, CONTRACT_INSIDE);

	SHRINK_AND_RETURN();
}

static int perform_shrink(values_vec &FV, simplex_mat &V, cost_function func, void *farg)
{
	params_vector best = column(V, 0);

	for (int j = 1; j < FV.size(); ++j) {
		column(V, j) = best + SIGMA * (column(V, j) - best);
		FV[j] = func(column(V, j), farg);
	}
	
	// number of function evaluations
	return FV.size() - 1;
}

} }
