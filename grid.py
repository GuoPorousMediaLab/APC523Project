#!/usr/bin/python
import functionsIO as wf
from time import clock

print "Read time"
print clock()
# write BoSlice.txt files to be used in VESA

gridName = 'grid'

init_condition = []

# density [kg/m^3]
init_condition.append(1)
# velocity in the x direction [m/s]
init_condition.append(1)
# velocity in the y direction [m/s]
init_condition.append(0)
# pressure [N/m^2]
init_condition.append(100000)


dx1 = 1
dx2 = 0.5
dy1 = 1
dy2 = 0.5
Nx1 = 1
Nx2 = 2
Ny1 = 2
Ny2 = 4

wf.writeGrid(gridName, init_condition, Nx1, Ny1, Nx2, Ny2, dx1, dy1, dx2, dy2)

print "Write Time"
print clock()

