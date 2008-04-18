#!/usr/bin/env python
import nm
from numpy import array

def rosenbrocks(vec):
  x,y = vec[0], vec[1]

  return (1-x)**2 + 100*((y-x*x)**2)

def rosenbrock3(vec):
  x,y,z = vec[0], vec[1], vec[2]

  return (1-x)**2 + 100*((y-x*x)**2) + (1-y)**2 + 100*((z-y*y)**2)

x0 = array([5.0, 5.0])
ret = nm.neldermead(x0, rosenbrocks, nm.dumbreporter)
#x0 = array([5.0, 5.0, 5.0])
#ret = nm.neldermead(x0, rosenbrock3, nm.dumbreporter)
print ret
