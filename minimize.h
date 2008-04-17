#ifndef __minimize_h
#define __minimize_h

#include <boost/numeric/ublas/vector.hpp>

namespace minimize {
    typedef boost::numeric::ublas::vector<double> params_vector;

    typedef double (*cost_function)(const params_vector &);
}


#endif
