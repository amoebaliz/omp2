# function [dist,phaseangle] = sw_dist(lat,lon,units)

#  SW_DIST    Distance between two lat,lon coordinates
# ===================================================================
#  SW_DIST  $Id: sw_dist.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
#           Copyright (C) CSIRO, Phil Morgan & Steve Rintoul 1992.
#
# USAGE:  [dist,phaseangle] = sw_dist(lat,lon {,units} )
#
# DESCRIPTION:
#   Calculate distance between two positions on glode using the "Plane
#   Sailing" method.  Also uses simple geometry to calculate the bearing of
#   the path between position pairs.
#
# INPUT:
#    lat      = decimal degrees (+ve N, -ve S) [- 90.. +90]
#    lon      = decimal degrees (+ve E, -ve W) [-180..+180]
#    units    = optional string specifing units of distance
#               'nm'  = nautical miles (default)
#               'km'  = kilometres
#
# OUTPUT:
#    dist        = distance between positions in units
#    phaseangle  = angle of line between stations with x axis (East).
#                  Range of values are -180..+180. (E=0, N=90, S=-90)
#
# AUTHOR:   Phil Morgan and Steve Rintoul 92-02-10
#
# DISCLAIMER:
#   This software is provided "as is" without warranty of any kind.
#   See the file sw_copy.m for conditions of use and licence.
#
# REFERENCE:
#    The PLANE SAILING method as descriibed in "CELESTIAL NAVIGATION" 1989 by
#    Dr. P. Gormley. The Australian Antartic Division.
# ==================================================================

# Modifications
# 99-06-25. Lindsay Pender, Function name change from distance to sw_dist.
# 99-06-25. Lindsay Pender, Fixed transpose of row vectors.

# CALLER:   general purpose
# CALLEE:   none

# ----------------------
#  CHECK INPUT ARGUMENTS
# ----------------------
def sw_dist(lat,lon,units):
    import numpy as np
    import math
    #if nargin > 3:
    #   raise ValueError('sw_dist.m: No more than 3 arguments allowed')
    #elif nargin==3:
    #   if not isinstance(units, str):
    #      raise ValueError('sw_dist.m: units argument must be string')
    #elif nargin==2:
    #   units = 'nm'  # default units
    #else:
    #   raise ValueError('sw_dist.m: wrong number of arguments')

    #mlat,nlat = lat.size
    #if ((mlat!=1) & (nlat!=1)):
    #   raise ValueError('sw_dist.m: lat, lon must be vectors.  No matrices allowed')

    if len(lat)!=len(lon):
       raise ValueError('sw_dist.m: lat and lon must have same number of elements')

    # -----------------
    # DEFINE CONSTANTS
    # -----------------

    DEG2RAD = (2*math.pi/360)
    RAD2DEG = 1/DEG2RAD
    DEG2MIN = 60
    DEG2NM  = 60
    NM2KM   = 1.8520    # Defined in Pond & Pickard p303.

    # BEGIN
    npositions = len(lat)
    # ind=1:npositions-1;     % index to first of position pairs
    # ind = range(npositions-1) # index to first of position pairs
    dlon = np.diff(lon.squeeze())
    if any(abs(dlon)>180):
       flag = np.where(abs(dlon)>180)[0]
       # for ii in range(len(flag)):
             # dlon(flag(ii))= -sign(dlon(flag(ii))) * (360 - abs(dlon(flag(ii))) );
       for ii in flag:
           dlon[ii]= -np.sign(dlon[ii]) * (360 - abs(dlon[ii]))
       #end for
    #end if
    latrad = abs(lat*DEG2RAD)
    temp_vals = (latrad[1:]+latrad[:-1])/2.
    dep    = np.cos(temp_vals.squeeze()) * dlon
    dlat   = np.diff(lat.squeeze())
   
    dist   = DEG2NM*np.sqrt(dlat**2 + dep**2)  # in n.miles

    if units == 'km':   # defaults to n.miles
       dist = dist * NM2KM
    # end %if

    # CALCUALTE ANGLE TO X AXIS
    phaseangle  = np.angle(dep+dlat*np.sqrt(-1))*RAD2DEG
    return dist, phaseangle
