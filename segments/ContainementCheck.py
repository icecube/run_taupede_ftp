# IceCube imports
from icecube import dataio, icetray, dataclasses
from icecube import millipede, MuonGun
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants
from icecube.dataclasses import I3RecoPulse, I3RecoPulseSeriesMap, I3RecoPulseSeriesMapMask
from icecube.icetray import I3Units, I3Frame, I3ConditionalModule, traysegment
from I3Tray import I3Tray

# Python system imports
import sys, datetime, os
from glob import glob
from optparse import OptionParser
import numpy as n





"""
containment selection
"""
def ispointinpolygon(x, y,outeredge_x,outeredge_y):


    x = n.asarray(x)
    y = n.asarray(y)
    cx = n.asarray(outeredge_x); cx = cx.reshape((cx.size, 1))
    cy = n.asarray(outeredge_y); cy = cy.reshape((cy.size, 1))

    # coordindates of the "next" vertex. Since the polygon is closed, the last vertex is next to the first.
    nx = n.roll(cx, 1, axis=0)
    ny = n.roll(cy, 1, axis=0)

    # draw a horizontal line at y. Which edges does it intersect?
    crossings = ((cy <= y)&(y < ny))|((ny <= y)&(y < cy))

    # now, cast a ray to the right. Which edges does it intersect?
    crossings &= (x < (nx-cx)*(y-cy)/(ny-cy) + cx)

    # count the number of crossings.
    inside = (crossings.sum(axis=0) % 2) != 0
    if inside.size == 1:
        return inside[0]
    else:
        return inside

"""
closest distance to outer detector boundary
"""
def getclosestdistance(x, y,outeredge_x,outeredge_y):

    

    x = n.asarray(x)
    y = n.asarray(y)
    cx = n.asarray(outeredge_x)
    cy = n.asarray(outeredge_y)

    # find two closest points
    mindex1, mindex2 = n.argsort(n.sqrt((cx-x)**2+(cy-y)**2))[0], n.argsort(n.sqrt((cx-x)**2+(cy-y)**2))[1]
    string1x, string1y = cx[mindex1], cy[mindex1]
    string2x, string2y = cx[mindex2], cy[mindex2]

    # project point onto connecting vector
    vec_a = (x-string1x, y-string1y)
    vec_b = (string2x-string1x, string2y-string1y)
    projection_c = (vec_a[0]*vec_b[0]+vec_a[1]*vec_b[1])/(vec_b[0]**2+vec_b[1]**2)
    vec_c = (projection_c*vec_b[0], projection_c*vec_b[1])
    ref_pos = (string1x+vec_c[0], string1y+vec_c[1])

    # calculate orthogonal distance to edge
    closestdistance = n.sqrt((x-ref_pos[0])**2+(y-ref_pos[1])**2)

    return closestdistance
"""
define vertex containment with respect to a boundary
"""
def iscontained(x, y, z, boundary,outeredge_x,outeredge_y):
    
    contained = True
    if (n.abs(z) > boundary) or (not ispointinpolygon(x, y,outeredge_x, outeredge_y) and (getclosestdistance(x, y,outeredge_x, outeredge_y) > boundary - 500)):
        contained = False
    return contained

