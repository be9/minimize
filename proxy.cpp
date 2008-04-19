#include <cmath>
#include "proxy.h"

namespace minimize { namespace proxy {

proxy::proxy(const constraints &_cons, const params_vector &initial, 
        cost_function _func, void *_farg)
    : cons(_cons), func(_func), farg(_farg)
{
    live_params = 0;

    for (int i = 0; i < cons.size(); ++i) {
        assert(cons[i].state != CONSTRAINED || cons[i].min < cons[i].max);
    
        if (cons[i].state != FROZEN)
            ++live_params;
    }

    frozen = initial;
}

double proxy::costfunc(const params_vector &vec, void *arg)
{
    const proxy *prox = reinterpret_cast<const proxy *>(arg);

    return prox->func(prox->from_raw(vec), prox->farg);
}

const double PI    = 3.14159265358979323846;    // pi
const double _1_PI = 0.31830988618379067154;    // 1/pi

// converts [-inf, +inf) -> [a, b)
inline double infinite2range(double value, double a, double b)
{
    return 0.5*(b+a) + (std::atan(value) * _1_PI) * (b-a);
}

// converts [a, b) -> [-inf, +inf)
inline double range2infinite(double value, double a, double b)
{
    assert(a <= value && value < b);

    return std::tan( (value - 0.5*(b+a)) * PI / (b-a) );
}

params_vector proxy::from_raw(const params_vector &raw) const
{
    assert(raw.size() == live_params);
    
    int i, j;
    params_vector real(cons.size());

    for (i = j = 0; i < real.size(); ++i) {
        switch (cons[i].state) {
            case FROZEN:
                real[i] = frozen[i];
                break;

            case LOOSE:
                real[i] = raw[j++];
                break;

            case CONSTRAINED:
                real[i] = infinite2range(raw[j++], cons[i].min, cons[i].max);
                break;
        }
    }

    return real;
}

params_vector proxy::to_raw(const params_vector &real) const
{
    assert(real.size() == cons.size());
    
    int i, j;
    params_vector raw(live_params);

    for (i = j = 0; i < real.size(); ++i) {
        switch (cons[i].state) {
            case FROZEN:
                continue;

            case LOOSE:
                raw[j++] = real[i];
                break;

            case CONSTRAINED:
                raw[j++] = range2infinite(real[i], cons[i].min, cons[i].max);
                break;
        }
    }

    return raw;
}

} }
