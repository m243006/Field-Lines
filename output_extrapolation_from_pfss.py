"""
GONG PFSS extrapolation
=======================

Calculating PFSS solution for a GONG synoptic magnetic field map.

To save to a JSON file, possibly use Binary JSON to reduce the size of the file

bjson.org
msgpack.org

"""
import json
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord

from sunpy.coordinates import frames, get_earth, transform_with_sun_center

import os

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import pfsspy.utils
import pfsspy

from pfsspy import coords, tracing
from pfsspy.sample_data import get_gong_map

###############################################################################
# Load a GONG magnetic field map
gong_fname = get_gong_map()
print(gong_fname)
gong_map = sunpy.map.Map(gong_fname)

###############################################################################
# Load HMI magnetic feild 
#hmi attempt ------------------------------------------
time = a.Time('2022/12/09', '2023/12/09')
series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')
crot = a.jsoc.PrimeKey('CAR_ROT', 2210)
result = Fido.search(time, series, crot,
                     a.jsoc.Notify('m243006@usna.edu'))

#downloads the files 
print(result)
files = Fido.fetch(result)

gong_map = sunpy.map.Map(files[0])
pfsspy.utils.fix_hmi_meta(gong_map)
#hmi_map.peek()

#-----------------------------------------------------------------

gong_fname = get_gong_map()
print(gong_fname)
#gong_map = sunpy.map.Map(gong_fname)


###############################################################################
# The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
# rho = ln(r), and r is the standard spherical radial coordinate. We need to
# define the number of rho grid points, and the source surface radius.
nrho = 35
rss = 2.5

###############################################################################
# From the boundary condition, number of radial grid points, and source
# surface, we now construct an Input object that stores this information
pfss_in = pfsspy.Input(gong_map, nrho, rss)


def set_axes_lims(ax):
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)


###############################################################################
# Using the Input object, plot the input field
#m = pfss_in.map
#fig = plt.figure()
#ax = plt.subplot(projection=m)
#m.plot()
#m.plot_settings["norm"]
#plt.colorbar()
#ax.set_title('Input field')
#set_axes_lims(ax)

###############################################################################
# Now calculate the PFSS solution
print('works')
pfss_out = pfsspy.pfss(pfss_in)
print('works')
###############################################################################
# Using the Output object we can plot the source surface field, and the
# polarity inversion line.

# Plot the map. Since are not interested in the exact map coordinates, we can
# simply use :meth:`~matplotlib.Axes.imshow`.
norm = gong_map.plot_settings['norm']
norm.vmin, norm.vmax = np.percentile(gong_map.data, [1, 99.9])
ax.imshow(gong_map.data,
          norm=norm,
          cmap=gong_map.plot_settings['cmap'],
          origin="lower")
fig.show()

ss_br = pfss_out.source_surface_br
# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=ss_br)

# Plot the source surface map
ss_br.plot()
# Plot the polarity inversion line
ax.plot_coord(pfss_out.source_surface_pils[0])
# Plot formatting
plt.colorbar()
ax.set_title('Source surface magnetic field')
set_axes_lims(ax)

###############################################################################
# It is also easy to plot the magnetic field at an arbitrary height within
# the PFSS solution.

# Get the radial magnetic field at a given height
ridx = 15
br = pfss_out.bc[0][:, :, ridx]
# Create a sunpy Map object using output WCS
br = sunpy.map.Map(br.T, pfss_out.source_surface_br.wcs)
# Get the radial coordinate
r = np.exp(pfss_out.grid.rc[ridx])

# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=br)

# Plot the source surface map
br.plot(cmap='RdBu')
# Plot formatting
plt.colorbar()
ax.set_title('$B_{r}$ ' + f'at r={r:.2f}' + '$r_{\\odot}$')
set_axes_lims(ax)


###############################################################################
# Finally, using the 3D magnetic field solution we can trace some field lines.
# In this case 64 points equally gridded in theta and phi are chosen and
# traced from the source surface outwards.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

tracer = tracing.FortranTracer()
r = 1.2 * const.R_sun
lat = np.linspace(-np.pi / 2, np.pi / 2, 8, endpoint=False)
lon = np.linspace(0, 2 * np.pi, 8, endpoint=False)
lat, lon = np.meshgrid(lat, lon, indexing='ij')
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad

seeds = SkyCoord(lon, lat, r, frame=pfss_out.coordinate_frame)

field_lines = tracer.trace(seeds, pfss_out)

# HGS observer time = 2018-08-11 00:00:00
# Representation is x, y, z
# Distance units in solar radii

# Set the observer an measurement units as expected by helios
field_line_observer = get_earth("2018-08-11 00:00:00")
field_line_observer.representation_type='cartesian'
d_unit = u.solRad
b_unit = u.G

# Get the observer information and put it in a dictionary
observer = {"obstime": {"value": field_line_observer.obstime.value,
                        "scale": field_line_observer.obstime.scale,
                        "format": field_line_observer.obstime.format},
            "x": {"value": field_line_observer.x.to(d_unit).value,
                    "unit": str(d_unit)},
            "y": {"value": field_line_observer.y.to(d_unit).value,
                    "unit": str(d_unit)},
            "z": {"value": field_line_observer.z.to(d_unit).value,
                       "unit": str(d_unit)},
            "frame": field_line_observer.name}

# Create the fieldlines dictionary
fieldlines = {"frame": {"x_unit": str(d_unit),
                                        "y_unit": str(d_unit),
                                        "z_unit": str(d_unit),
                                        "coordinate_system": "Heliographic Stonyhurst",
                            "source_map_obstime": {"value": gong_map.date.value,
                                            "scale": gong_map.date.scale,
                                            "format": gong_map.date.format}},
              "field_description": {"b_unit": str(b_unit),
                                    "bx_value": "x component of field vector",
                                    "by_value": "y component of field vector",
                                    "bz_value": "z component of field vector",
                                    "b_mag": "magnitude of field"},
              "lines": []}
#------------------------------------------------------------------------------------------------------------
# Go through the field lines and extract information 
this_field_line = -1
for field_line in field_lines:
    flc = field_line.coords
    if len(flc) > 0:
        this_field_line = this_field_line + 1
        with transform_with_sun_center():
            hgs = field_line.coords.transform_to(field_line_observer).transform_to(frames.HeliographicStonyhurst)
            hgs.representation_type='cartesian'
            b_along_fline = field_line.b_along_fline.to(b_unit).value
            fieldlines["lines"].append({"x":hgs.x.to(d_unit).value.tolist(),
                                            "y": hgs.y.to(d_unit).value.tolist(),
                                            "z": hgs.z.to(d_unit).value.tolist(),
                                            "polarity":field_line.polarity,
                                            "bx_value": b_along_fline[:,0].tolist(),
                                            "by_value": b_along_fline[:,1].tolist(),
                                            "bz_value": b_along_fline[:,2].tolist(),
                                            "b_mag": np.sqrt(b_along_fline[:,0]**2 + b_along_fline[:,1]**2 + b_along_fline[:,2]**2).tolist()})

            color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
            ax.plot(hgs.x.to(d_unit),
                        hgs.y.to(d_unit),
                        hgs.z.to(d_unit), color=color, linewidth=1)

#d unit is the solar radius 
ax.set_title('PFSS solution')
plt.show()


# Write out the PFSS field line information
output = {'observer': observer, 'fieldlines': fieldlines}
filename = "test_pfsspy2.json"
with open(filename, "w") as outfile:
    json.dump(output, outfile)

