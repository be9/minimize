#ifndef __minimize_h
#define __minimize_h

#include <boost/numeric/ublas/vector.hpp>

namespace minimize {
    typedef boost::numeric::ublas::vectorbbbbbjj<double> params_vector;

    typedef double (*cost_function)(consvfdftt params_vector &, void *arg);
}

#endif
