#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "../proxy.h"
#include <iostream>
#include <iomanip>

BOOST_AUTO_TEST_SUITE(proxy_suite);

using namespace std;
using namespace minimize::proxy;

#define PV minimize::params_vector
#define CS constraint_state
#define K BOOST_CHECK

#define BOOST_CHECK_CLOSE_MESSAGE( L, R, T, M ) \
    BOOST_CHECK_WITH_ARGS_IMPL( ::boost::test_tools::check_is_close, M, CHECK, CHECK_CLOSE, \
        (L)(R)(::boost::test_tools::percent_tolerance(T)) )

inline bool operator==(const PV &v1, const PV &v2)
{
    if (v1.size() != v2.size())
        return false;

    for (int i = 0; i < v1.size(); ++i)
        if (v1[i] != v2[i])
            return false;

    return true;
}

PV vec(double a) { PV v(1); v[0] = a; return v; }
PV vec(double a, double b) { PV v(2); v[0] = a; v[1] = b; return v; }
PV vec(double a, double b, double c) { PV v(3); v[0] = a; v[1] = b; v[2] = c; return v; }

constraints cons(CS a, CS b, CS c) 
{ 
    constraints cs(3); 
    
    cs[0].state = a;
    cs[1].state = b;
    cs[2].state = c;

    return cs;
}

static void check_vec(const constraints &cs, const PV &v1, const PV &v2)
{
    proxy p(cs, v1, 0, 0);

    PV raw = p.to_raw(v1);

    K(raw == v2);
    K(p.from_raw(raw) == v1);
}

BOOST_AUTO_TEST_CASE(proxy_passthru)
{
    check_vec(cons(LOOSE, LOOSE, LOOSE), vec(1,2,3), vec(1,2,3));
}

BOOST_AUTO_TEST_CASE(proxy_frozen)
{
    check_vec(cons(FROZEN, LOOSE,  LOOSE),  vec(1,2,3), vec(2,3));
    check_vec(cons(LOOSE,  FROZEN, LOOSE),  vec(1,2,3), vec(1,3));
    check_vec(cons(LOOSE,  LOOSE,  FROZEN), vec(1,2,3), vec(1,2));
    check_vec(cons(FROZEN, FROZEN, LOOSE),  vec(1,2,3), vec(3));
    check_vec(cons(FROZEN, LOOSE,  FROZEN), vec(1,2,3), vec(2));
    check_vec(cons(LOOSE,  FROZEN, FROZEN), vec(1,2,3), vec(1));
}

static PV was_vec;
void *was_arg;
static int called;

static double myfunc(const PV &vec, void *arg)
{
    was_vec = vec;
    was_arg = arg;

    ++called;

    return 1.2345;
}

BOOST_AUTO_TEST_CASE(proxy_fun_call)
{
    PV vec123 = vec(1,2,3);
    proxy p(cons(LOOSE, FROZEN, LOOSE), vec123, myfunc, &vec123);

    called = 0;
    K(proxy::costfunc(vec(1,3), &p) == 1.2345);

    K(called == 1);
    K(was_vec == vec123);
    K(was_arg == &vec123);
}

static void check_conv(double min, double max, double v1, double v2, double tol = 1e-4)
{
    constraints cs(1);
    cs[0].state = CONSTRAINED;
    cs[0].min = min;
    cs[0].max = max;

    PV vec1 = vec(v1), vec2 = vec(v2);
    proxy p(cs, vec1, 0, 0);

    PV raw = p.to_raw(vec1);

    BOOST_CHECK_CLOSE(raw[0], vec2[0], tol); 

    PV real = p.from_raw(vec2);
  
    BOOST_CHECK_CLOSE(real[0], vec1[0], tol); 
}

BOOST_AUTO_TEST_CASE(proxy_cons)
{
    check_conv(1, 2, 1.001, -318, 1);
    check_conv(1, 2, 1.01, -31.8205);
    check_conv(1, 2, 1.5,  0);
    check_conv(1, 2, 1.99, 31.8205);
    check_conv(1, 2, 1.999, 318.205, 1);
}

BOOST_AUTO_TEST_SUITE_END();
