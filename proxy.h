#ifndef __minimize_proxy_h
#define __minimize_proxy_h

#include <vector>
#include "minimize.h"

// Minimize proxy
namespace minimize { namespace proxy {

enum constraint_state { 
    LOOSE, CONSTRAINED, FROZEN
};

struct constraint {
    constraint_state state;

    double min, max;        // valid when state == CONSTRAINED
};

typedef std::vector<constraint> constraints;

class proxy {
public:
    proxy(const constraints &cons, const params_vector &initial, 
        cost_function func, void *farg);

    // converts "real" vector into a "raw" vector used by minimization algos
    params_vector to_raw(const params_vector &real) const;
    
    // converts "raw" back to "real"
    params_vector from_raw(const params_vector &raw) const;

    // expects proxy * as arg!
    static double costfunc(const params_vector &vec, void *arg);

private:
    const constraints cons;
    cost_function func;
    void *farg;
    
    int live_params;
    params_vector frozen;
};

} }

#endif
