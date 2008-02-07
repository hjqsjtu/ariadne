#!/usr/bin/python
#
# integration example
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
import math


# y' = z & z' = -y
dyn=AffineVectorField(Matrix([[0,0,0],[0,0,1],[0,-1,0]]),Vector([1,0,0]))

grid=Grid(Vector([0.02,0.01,0.01]))
block=Box([[-0.1,7.1],[-1.2,1.2],[-1.2,1.2]])
fgrid=FiniteGrid(grid,block)
print fgrid

time=7

print "Creating reference set"
reference_set=GridMaskSet(fgrid)
for i in range(71):
	x = 0.1*i
	y = math.cos(x)
	print x, y
        if y <= 1 and y >= -1 :
	  reference_set.adjoin_outer_approximation(Box([[x,x],[y,y],[0.001,0.001]]))
	  	
print "reference_set.size(),capacity()=",reference_set.size(),reference_set.capacity()

print "Creating inital set"
initial_set=RectangularSet([[0.001,0.002],[0.999,1],[0.001,0.002]])

# Evolution parameters
par=EvolutionParameters()
par.set_maximum_step_size(0.5);
par.set_lock_to_grid_time(time);
par.set_grid_length(0.01);


integrator=AffineIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing upper reach set with Affine Integrator..."
reach_set=evolver.lower_reach(dyn,initial_set,Rational(time))

print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("cos-affine.eps",block,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done.\n"

integrator=LohnerIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing chainreach sets with Lohner Integrator..."
chainreach_set=evolver.upper_reach(dyn,initial_set,Rational(time))

print "Exporting to postscript output...",
eps.open("cos-lohner.eps",block,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=KuhnIntegrator(3,3)
evolver=VectorFieldEvolver(par,integrator)

print "Computing chainreach sets with Kuhn Integrator..."
chainreach_set=evolver.upper_reach(dyn,initial_set,Rational(time))

eps.open("cos-kuhn.eps",block,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

  	

