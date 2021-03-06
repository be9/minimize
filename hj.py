# Nonlinear Optimization using the algorithm of Hooke and Jeeves
#	12 February 1994	author: Mark G. Johnson

# Find a point X where the nonlinear function f(X) has a local minimum. X is an
# n-vector and f(X) is a scalar. In mathematical notation f: R^n -> R^1.  The
# objective function f() is not required to be continuous. Nor does f() need to
# be differentiable. The program does not use or require derivatives of f().

# The software user supplies three things: a subroutine that computes f(X), an
# initial "starting guess" of the minimum point X, and values for the algorithm
# convergence parameters. Then the program searches for a local minimum,
# beginning from the starting guess, using the Direct Search algorithm of Hooke
# and Jeeves.

# This C program is adapted from the Algol pseudocode found in "Algorithm 178:
# Direct Search" by Arthur F. Kaupe Jr., Communications of the ACM, Vol 6.
# p.313 (June 1963).  It includes the improvements suggested by Bell and Pike
# (CACM v.9, p. 684, Sept 1966) and those of Tomlin and Smith, "Remark on
# Algorithm 178" (CACM v.12). The original paper, which I don't recommend as
# highly as the one by A. Kaupe, is: R. Hooke and T. A. Jeeves, "Direct Search
# Solution of Numerical and Statistical Problems", Journal of the ACM, Vol. 8,
# April 1961, pp. 212-229.

# Calling sequence:
#  int hooke(nvars, startpt, endpt, rho, epsilon, itermax)
#
#     nvars	   {an integer}  This is the number of dimensions
#		   in the domain of f().  It is the number of
#		   coordinates of the starting point (and the
#		   minimum point.)
#     startpt	   {an array of doubles}  This is the user-
#		   supplied guess at the minimum.
#     endpt	   {an array of doubles}  This is the location of
#		   the local minimum, calculated by the program
#     rho	   {a double}  This is a user-supplied convergence
#		   parameter (more detail below), which should be
#		   set to a value between 0.0 and 1.0.	Larger
#		   values of rho give greater probability of
#		   convergence on highly nonlinear functions, at a
#		   cost of more function evaluations.  Smaller
#		   values of rho reduces the number of evaluations
#		   (and the program running time), but increases
#		   the risk of nonconvergence.	See below.
#     epsilon	   {a double}  This is the criterion for halting
#		   the search for a minimum.  When the algorithm
#		   begins to make less and less progress on each
#		   iteration, it checks the halting criterion: if
#		   the stepsize is below epsilon, terminate the
#		   iteration and return the current best estimate
#		   of the minimum.  Larger values of epsilon (such
#		   as 1.0e-4) give quicker running time, but a
#		   less accurate estimate of the minimum.  Smaller
#		   values of epsilon (such as 1.0e-7) give longer
#		   running time, but a more accurate estimate of
#		   the minimum.
#     itermax	   {an integer}  A second, rarely used, halting
#		   criterion.  If the algorithm uses >= itermax
#		   iterations, halt.


# The user-supplied objective function f(x,n) should return a C
# "double".  Its  arguments are  x -- an array of doubles, and
# n -- an integer.  x is the point at which f(x) should be
# evaluated, and n is the number of coordinates of x.	That is,
# n is the number of coefficients being fitted.

# rho, the algorithm convergence control
#	The algorithm works by taking "steps" from one estimate of
#    a minimum, to another (hopefully better) estimate.  Taking
#    big steps gets to the minimum more quickly, at the risk of
#    "stepping right over" an excellent point.  The stepsize is
#    controlled by a user supplied parameter called rho.  At each
#    iteration, the stepsize is multiplied by rho  (0 < rho < 1),
#    so the stepsize is successively reduced.
#	Small values of rho correspond to big stepsize changes,
#    which make the algorithm run more quickly.  However, there
#    is a chance (especially with highly nonlinear functions)
#    that these big changes will accidentally overlook a
#    promising search vector, leading to nonconvergence.
#	Large values of rho correspond to small stepsize changes,
#    which force the algorithm to carefully examine nearby points
#    instead of optimistically forging ahead.	This improves the
#    probability of convergence.
#	The stepsize is reduced until it is equal to (or smaller
#    than) epsilon.  So the number of iterations performed by
#    Hooke-Jeeves is determined by rho and epsilon:
#	    rho**(number_of_iterations) = epsilon
#	In general it is a good idea to set rho to an aggressively
#    small value like 0.5 (hoping for fast convergence).  Then,
#    if the user suspects that the reported minimum is incorrect
#    (or perhaps not accurate enough), the program can be run
#    again with a larger value of rho such as 0.85, using the
#    result of the first minimization as the starting guess to
#    begin the second minimization.

# Normal use: (1) Code your function f() in the C language
#	       (2) Install your starting guess {or read it in}
#	       (3) Run the program
#	       (4) {for the skeptical}: Use the computed minimum
#		      as the starting point for another run

# Data Fitting:
#	Code your function f() to be the sum of the squares of the
#	errors (differences) between the computed values and the
#	measured values.  Then minimize f() using Hooke-Jeeves.
#	EXAMPLE: you have 20 datapoints (ti, yi) and you want to
#	find A,B,C such that  (A*t*t) + (B*exp(t)) + (C*tan(t))
#	fits the data as closely as possible.  Then f() is just
#	f(x) = SUM (measured_y[i] - ((A*t[i]*t[i]) + (B*exp(t[i]))
#				  + (C*tan(t[i]))))^2
#	where x[] is a 3-vector consisting of {A, B, C}.

#
#  The author of this software is M.G. Johnson.
#  Permission to use, copy, modify, and distribute this software
#  for any purpose without fee is hereby granted, provided that
#  this entire notice is included in all copies of any software
#  which is or includes a copy or modification of this software
#  and in all copies of the supporting documentation for such
#  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT
#  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE
#  AUTHOR NOR AT&T MAKE ANY REPRESENTATION OR WARRANTY OF ANY
#  KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS
#  FITNESS FOR ANY PARTICULAR PURPOSE.
#

from numpy import array, zeros, shape, ones, arange

SUCCESS = 0
ITERMAX = 1
BREAK = 2

def dumbreporter(place, itercount, func_evals, newx, newf):
  print "Iteration %d (%d func. evals): %s\n" % (itercount, func_evals, place)
  print "X = ", newx, "; F = ", newf
  return False

# given a point, look for a better one nearby, one coord at a time
def _best_nearby(delta, point, prevbest, func):
  minf = prevbest
  nvars = shape(point)[0]
  fevals = 0

  z = point.copy()

  for i in range(nvars):
    # try positive step
    z[i] = point[i] + delta[i]
    
    ftmp, fevals = func(z), fevals+1

    if ftmp < minf:
      minf = ftmp
      continue

    # try negative step
    delta[i] = -delta[i]
    z[i] = point[i] + delta[i]

    ftmp, fevals = func(z), fevals+1

    if ftmp < minf:
      minf = ftmp
      continue

    # leave as is
    z[i] = point[i]

  point[:] = z

  return minf, fevals

def hooke_jeeves(x0, func, outputfunc = None, rho=0.5, epsilon=1e-6, itermax=5000):
  nvars = shape(x0)[0]
  
  delta = abs(x0 * rho)

  for i in range(nvars):
    if delta[i] == 0.0:
      delta[i] = rho

  iters = 0
  steplength = rho

  xbefore = x0
  fbefore = func(xbefore)
  
  newx = x0
  newf = fbefore

  func_evals = 1

  def call_output(place):
    return outputfunc and outputfunc(place, iters, func_evals, newx, newf)

  call_output('init')

  while iters < itermax and steplength > epsilon:
    iters += 1

    if call_output('iter'): 
      return BREAK, newx, newf

    # find best new point, one coord at a time
    newx = xbefore.copy()
    newf, fevals = _best_nearby(delta, newx, fbefore, func)
    func_evals += fevals

		# if we made some improvements, pursue that direction
    while newf < fbefore:
      # firstly, arrange the sign of delta[]
      for i in range(nvars):
        if newx[i] <= xbefore[i]:
          delta[i] = -abs(delta[i])
        else:
          delta[i] = abs(delta[i])
        
      # now, move further in this direction
      newx, xbefore = 2*newx-xbefore, newx

      fbefore = newf
      newf, fevals = _best_nearby(delta, newx, fbefore, func)
      func_evals += fevals
      
      # if the further (optimistic) move was bad....
      if newf >= fbefore:
        break

      # make sure that the differences between the new and the old points are
      # due to actual displacements; beware of roundoff errors that might cause
      # newf < fbefore
      if max(abs(newx-xbefore) - 0.5*abs(delta)) <= 0.0:
        break
    
    if steplength >= epsilon and newf >= fbefore:
      steplength *= rho
      delta *= rho

  reason = SUCCESS
  if iters >= itermax:
    reason = ITERMAX

  return reason, xbefore, fbefore
