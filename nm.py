from numpy import array, zeros, shape, ones, arange

# return codess for neldermead
SUCCESS = 0
BREAK   = 1
MAX_FUNC_EVALS = 2
MAX_ITER = 3

# Nelder-Mead algorithm constants
RHO   = 1
CHI   = 2
PSI   = 0.5
SIGMA = 0.5

def dumbreporter(what, itercount, func_evals, details, xval, fval, v, fv):
  pass
  
# _change_simplex finds the direction where to go and transforms the simplex
#
# v    - an NxN+1-sized matrix with the simplex 
# fv   - an N+1-sized vector. fv[i] = func(v[:,i])
# func - cost function.

# returns a tuple ('transformation', number of function evaluations made)
def _change_simplex(v, fv, func):
  # handy func to replace the worst point with a new one
  def replace_worst(newx, newf):
    v[:,-1] = newx
    fv[-1] = newf

  evals = 0

  # evaluate function and calculate the number of evaluations
  def evaluate(x):
    evals += 1
    return func(x)

  n = numpy.shape(v)[0]
  
  # perform a shrink
  def shrink():
    for j in range(1, n+1):
      v[:,j] = v[:,0] + SIGMA * (v[:,j] - v[:,0])
      fv[j] = evaluate(v[:,j])

  # xbar = average of the n (NOT n+1) best points
  xbar = v[:, 0:n].mean(1)
  worst = v[:,-1]

  # Compute the reflection point
  xr = (1 + RHO)*xbar - RHO*worst
  fxr = evaluate(xr)

  # >>> Check the reflection point against our current best
  if fxr < fv[0]:
    # Calculate the expansion point
    xe = (1 + RHO*CHI)*xbar - RHO*CHI*worst
    fxe = evaluate(xe)

    if fxe < fxr:
      replace_worst(xe, fxe)
      how = 'expand'
    else:
      replace_worst(xr, fxr)
      how = 'reflect'

    return how, evals

  # >>> Continuing, the reflection point is worse than the current best
  if fxr < fv[-2]:
    # >>> but at least it's better than the second worst point from the end
    replace_worst(xr, fxr)

    return 'reflect', evals

  # >>> Worse than the second worse point from the end
  
  # Perform contraction
  if fxr < fv[-1]:
    # >>> better than the worst point we had

    # Perform an outside contraction
    xc = (1 + PSI*RHO)*xbar - PSI*RHO*worst
    fxc = evaluate(xc)
          
    if fxc <= fxr:
        replace_worst(xc, fxc)
        return 'contract outside', evals
    else:
        # perform a shrink
        shrink()
        return 'shrink', evals
    
  # Reflection bad, really bad, even worse than the worst one

  # Perform an inside contraction
  xcc = (1-PSI)*xbar + PSI*worst
  fxcc = evaluate(xcc)
          
  if fxcc < fv[-1]:
    replace_worst(xcc, fxcc)
    return 'contract inside', evals

  shrink()
  return 'shrink', evals

def neldermead(x0, func, outputfunc = None, maxiter = -200, maxfunevals = -200, \
               tolx = 1e-4, tolfun = 1e-4, funvalcheck = False):
  n = shape(x0)[0]
  x = x0.copy()

  if maxiter < 0:     maxiter *= n
  if maxfunevals < 0: maxfunevals *= n

  # Set up a simplex near the initial guess.
  xin = x.copy()
  v = zeros((n, n+1))
  fv = zeros(n+1)
 
  # Place input guess in the simplex! (credit L.Pfeffer at Stanford)
  v[:,0] = xin
  fv[0] = func(xin)
  
  func_evals = 1
  itercount = 0

  # Initial simplex setup continues later

  # calls output and returns true if we need to break
  def call_output(place, details):
    return outputfunc and outputfunc(place, itercount, func_evals, details, v, fv) 

  if call_output('init', ''): return BREAK

  # 0th iteration, initial f(x)
  if call_output('iter', ''): return BREAK

  # Continue setting up the initial simplex.
  # Following improvement suggested by L.Pfeffer at Stanford
  
  usual_delta = 0.05              # 5 percent deltas for non-zero terms
  zero_term_delta = 0.00025       # Even smaller delta for zero elements of x

  for j in range(n):
    y = xin.copy()

    if y[j] != 0.0:
      y[j] *= 1 + usual_delta
    else:
      y[j] = zero_term_delta

    v[:,j+1] = y
    fv[j+1] = func(y)

  # sort so v[0,:] has the lowest function value
  def sort_points():
    indexes = range(shape(v)[1])
    indexes.sort(key=lambda i: fv[i])
    fv = fv[indexes]
    v = v[:, indexes]

  sort_points()

  itercount += 1
  func_evals = n+1
  
  if call_output('iter', 'initial simplex'): return BREAK

  # Main algorithm
  # Iterate until the diameter of the simplex is less than tolx
  #   AND the function values differ from the min by less than tolf,
  #   or the max function evaluations are exceeded. (Cannot use OR instead of
  #   AND.)
  while func_evals < maxfun and itercount < maxiter:
    if (fv[0] - fv[1:n+1]).abs().max() <= tolf and (v[:,1:n+1] - v[:,0]).abs().max() <= tolx:
      break

    what, evals = _change_simplex(v, fv)
    func_evals += evals

    sort_points()
    itercount += 1

    if call_output('iter', what): 
      return BREAK

  call_output('done', '')

  if func_evals >= maxfun:
    return MAX_FUNC_EVALS

  if itercount >= maxiter:
    return MAX_ITER

  return SUCCESS
