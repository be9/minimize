/* 
Copyright notice from the original hooke.c source:

The author of this software is M.G. Johnson. Permission to use, copy, modify,
and distribute this software for any purpose without fee is hereby granted,
provided that this entire notice is included in all copies of any software
which is or includes a copy or modification of this software and in all copies
of the supporting documentation for such software. THIS SOFTWARE IS BEING
PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED WARRANTY. IN PARTICULAR,
NEITHER THE AUTHOR NOR AT&T MAKE ANY REPRESENTATION OR WARRANTY OF ANY KIND
CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS FITNESS FOR ANY
PARTICULAR PURPOSE.
*/
#include <cmath>
#include <algorithm>
#include "hj.h"

namespace minimize { namespace hooke_jeeves {

static int best_nearby(params_vector &delta, params_vector &point, double &pval, double prevbest, cost_function func);

static params_vector abs(const params_vector &x)
{
    params_vector result(x.size());

    for (int i = 0; i < x.size(); ++i)
        result[i] = std::abs(x[i]);

    return result;
}

inline double max(const params_vector &vec)
{
    return *std::max_element(vec.begin(), vec.end());
}

void run(const params_vector &x0, cost_function func, output &out, report_func reporter,
        int maxiter, double rho, double epsilon)
{
    int i, nvars = x0.size();
    params_vector delta(nvars);

    for (i = 0; i < nvars; ++i) {
        double d = std::abs(x0[i] * rho);
        delta[i] = (d == 0.0) ? rho : d;
    }

    int iters = 0;
    double steplength = rho;

    params_vector xbefore = x0;
    double fbefore = func(xbefore);

    params_vector newx = x0;
    double newf = fbefore;

    int func_evals = 1;
    
    reporter(INIT, iters, func_evals, newx, newf);

    bool abort = false;

    while (iters < maxiter && steplength > epsilon) {
        ++iters;

        if (reporter(ITER, iters, func_evals, newx, newf)) {
            abort = true;
            break;
        }

        newx = xbefore;
        func_evals += best_nearby(delta, newx, newf, fbefore, func);
		
        // if we made some improvements, pursue that direction
        while (newf < fbefore) {
            // firstly, arrange the sign of delta[]
            for (i = 0; i < nvars; ++i) 
                delta[i] = (newx[i] <= xbefore[i]) ? -std::abs(delta[i]) : std::abs(delta[i]);

            // now, move further in this direction
            params_vector tmp = newx;
            newx = 2.0*newx - xbefore;
            xbefore = tmp;
            fbefore = newf;

            func_evals += best_nearby(delta, newx, newf, fbefore, func);

            // if the further (optimistic) move was bad....
            if (newf >= fbefore)
                break;

            // make sure that the differences between the new and the old
            // points are due to actual displacements; beware of roundoff
            // errors that might cause newf < fbefore
            if (max(abs(newx-xbefore) - 0.5*abs(delta)) <= 0.0)
               break;
        }

        if (steplength >= epsilon && newf >= fbefore) {
            steplength *= rho;
            delta *= rho;
        }
    }

    if (abort)
        out.code = ABORT;
    else if (iters >= maxiter)
        out.code = MAX_ITER;
    else
        out.code = SUCCESS;

    out.last_x = xbefore;
    out.last_f = fbefore;
    out.iterations = iters;
    out.func_evaluations = func_evals;
}

static int best_nearby(params_vector &delta, params_vector &point, double &pval, double prevbest, cost_function func)
{
    double minf = prevbest;
    int nvars = point.size();
    int fevals = 0;

    params_vector z = point;

    for (int i = 0; i < nvars; ++i) {
        // try positive step
        z[i] = point[i] + delta[i];
        
        double ftmp = func(z); ++fevals;

        if (ftmp < minf) {
            minf = ftmp;
            continue;
        }

        // try negative step
        delta[i] = -delta[i];
        z[i] = point[i] + delta[i];
        
        ftmp = func(z); ++fevals;

        if (ftmp < minf) {
            minf = ftmp;
            continue;
        }

        // leave as is
        z[i] = point[i];
    }

    point = z;
    pval = minf;

    return fevals;
}

} }
