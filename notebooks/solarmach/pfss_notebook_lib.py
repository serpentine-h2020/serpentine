'''
A library intended for importing into pfsspy_notebook.ibynp. It contains the necessary functions 
for seeking the footpoints of IMF field lines connecting back to the photosphere and plotting them.

@Author: Christian Palmroos
         <chospa@utu.fi>

Last updated: 2022-03-02
'''

#imports:
import astropy.constants as const
import astropy.units as units
import astropy.coordinates
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.colors as mcolors
import cmasher as cmr
import numpy as np
import pandas as pd
import pfsspy
#import scipy
import pickle
import warnings
import sunpy.map
#import sunpy.io.fits

#from astropy.time import Time
from astropy.coordinates import SkyCoord
#from pfsspy import coords
from pfsspy import tracing
#from pfsspy.sample_data import get_gong_map
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
from sunpy.net import Fido, attrs as a
#from sunpy.coordinates import frames

import os

#some matplotlib settings:
plt.rcParams['axes.linewidth'] = 1.5
# plt.rcParams['font.size'] = 20
plt.rcParams['agg.path.chunksize'] = 20000

#living on the edge:
warnings.simplefilter('ignore')

# ----------------------------------------------------------------------------------------

# functions:

def get_color(key: str = None) -> str:
    '''
    Returns the right color for an object according to SOLAR-MACH tool.
    '''
    if key is not None:
        key = key.lower()
    else:
        key = 'default'

    color_dict = {'mercury': 'darkturquoise',
                'venus': 'darkorchid',
                'earth': 'green',
                'mars': 'maroon',
                'jupiter': 'navy',
                'stereo a': 'red',
                'stereo-a': 'red',
                'stereo b': 'b',
                'stereo-b': 'b',
                'soho': 'darkgreen',
                'solo': 'dodgerblue',
                'solar orbiter': 'dodgerblue',
                'psp': 'purple',
                'parker solar probe': 'purple',
                'bepicolombo': 'orange',
                'maven': 'brown',
                'mars express': 'darkorange',
                'messenger': 'olivedrab',
                'juno': 'orangered',
                'cassini': 'mediumvioletred',
                'rosetta': 'blueviolet',
                'pioneer 10': 'teal',
                'pioneer 11': 'darkblue',
                'ulysses': 'dimgray',
                'voyager 1': 'darkred',
                'voyager 2': 'midnightblue',
                'default': 'k'}


    return color_dict[key]

# ----------------------------------------------------------------------------------------

def get_sc_data(csvfile: str):
    '''
    Reads the contents of solar-mach produced csv file, and returns lists
    of necessary data to run pfss field line tracing analysis.
    
    csvfile: str, the name of the csv file one wants to read
    '''
    import pandas as pd

    if type(csvfile) is not str:
        raise TypeError("File name is not a string.")

    csvdata = pd.read_csv(csvfile)
    
    names = list(csvdata['Spacecraft/Body'])
    lons = list(csvdata['Carrington Longitude (°)'])
    lats = list(csvdata['Latitude (°)'])
    dist = au_to_km(list(csvdata['Heliocentric Distance (AU)']))
    sw = list(csvdata['Vsw'])
    
    return names, sw, dist, lons, lats

# ----------------------------------------------------------------------------------------

def field_line_accuracy(flines):
    '''
    Calculates at least now the central point, average distance from the central point and 
    the standard deviation of the photospheric footpoints for a set of field lines
    '''

    if isinstance(flines[0], list):
        flines = flatten(flines)

    footpoints = []
    for fline in flines:

        r, lon, lat = get_coord_values(fline)

        #index 0 of coordinates corresponds to photospheric coordinates, index -1 to pfss
        footpoint = [lon[0], lat[0]]
        footpoints.append(footpoint)

    #declare longitudes and latitudes of footpoints
    footp_lons = [pair[0] for pair in footpoints]
    footp_lats = [pair[1] for pair in footpoints]

    #if separation in longitudes is over half a circle, then points are probably
    #scattered near the 0-deg line -> transfer them far enough from that line
    #while calculations are ran
    if max(footp_lons) > min(footp_lons)+180.0:

        #transfer of longitudes
        footp_lons = shift_longitudes(footp_lons)

        #calculate the central point
        c_point = [np.mean(footp_lons), np.mean(footp_lats)]

        #standard deviation of longitudes and latitudes
        lon_std = np.std(footp_lons)
        lat_std = np.std(footp_lats)

        #calculate mean distance from the central point
        dist_sum = 0
        for i in range(len(footp_lons)):

            lon1, lat1 = footp_lons[i], footp_lats[i]
            angular_separation = orthodrome(lon1, lat1, c_point[0], c_point[1])

            #distance is in solar radii
            distance_rs = angular_separation
            dist_sum = dist_sum + distance_rs

        #transfer lons and the central point back the same amount that the longitudes were moved
        footp_lons = shift_longitudes(footp_lons, shift=-180.0)
        c_point[0] = shift_longitudes([c_point[0]], shift=-180.0)[0]

    else:
        
        #calculate the central point
        c_point = [np.mean(footp_lons), np.mean(footp_lats)]

        #standard deviation of longitudes and latitudes
        lon_std = np.std(footp_lons)
        lat_std = np.std(footp_lats)

        #calculate mean distance from the central point
        dist_sum = 0
        for i in range(len(footp_lons)):

            lon1, lat1 = footp_lons[i], footp_lats[i]
            angular_separation = orthodrome(lon1, lat1, c_point[0], c_point[1])

            #distance is in solar radii
            distance_rs = angular_separation
            dist_sum = dist_sum + distance_rs

    avg_dist = dist_sum/len(footp_lons)


    return footpoints, c_point, avg_dist, [lon_std, lat_std]

# ----------------------------------------------------------------------------------------

def shift_longitudes(lons, shift=180.0):
    '''
    Shifts the longitudes by <shift> amount
    '''
    if shift>0:
        lons = [lon+shift for lon in lons]
        lons = [lon-360.0 if lon>360 else lon for lon in lons]
    else:
        #if shift is negative, points are likely being moved back to their 
        #original place
        lons = [lon+shift for lon in lons]
        lons = [lon+360.0 if lon<0 else lon for lon in lons]

    return lons

# ----------------------------------------------------------------------------------------

def map_on_surface(fps, c_point, avg_d, shift=None, zoom=None, show_avg_d=False):
    '''
    Plots a really simple 2d representation of fieldline objects' footpoints.
    '''
    
    
    import matplotlib.patches as mpatch

    centre = np.array(c_point)
    fpslons = [item[0] for item in fps]
    fpslats = [item[1] for item in fps]

    if shift is not None:
        fpslons = shift_longitudes(fpslons, shift=shift)
        centre[0] = c_point[0]+shift
        
    fig_tuple = (16,7)
    if zoom is not None:
        fig_tuple = (12,12)

    fig = plt.figure(figsize=fig_tuple)
    ax = plt.subplot()

    ax.scatter(fpslons[0], fpslats[0], c='navy', s=60, label="original footpoint")
    ax.scatter(fpslons[1:], fpslats[1:], c='C0', label="dummy footpoints")
    ax.scatter(centre[0], centre[1], c='r', label="avg(lons,lats)")
    
    if show_avg_d:
        avg_d_deg = np.rad2deg(avg_d)
        ax.add_patch(mpatch.Circle((centre[0],centre[1]), avg_d_deg, color='r', lw=0.8, ls='--', fill=False))

    plt.ylim(-90,90)
    plt.xlim(0,360)
    if zoom is not None:
        plt.ylim(centre[1]-zoom/2, centre[1]+zoom/2)
        plt.xlim(centre[0]-zoom/2, centre[0]+zoom/2)

    plt.legend()
    plt.grid("True")
    plt.show()

# ----------------------------------------------------------------------------------------

def orthodrome(lon1,lat1, lon2,lat2, rad=False) -> float:
    '''
    calculates the othodrome (total angular separtion) between two coordinates defined by their lon/lat positions
    '''
    import numpy as np

    if(rad == False):
        lon1 = np.deg2rad(lon1)
        lon2 = np.deg2rad(lon2)
        lat1 = np.deg2rad(lat1)
        lat2 = np.deg2rad(lat2)

    ortho = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))

    return ortho

#----------------------------------------------------------------------------------------

def ortho_to_points(lon1, lat1, orthodrome, rad=False):
    '''
    Calculates a lon/lat pair from a central point and an orthodrome(=total angular separation between points)
    '''

    if rad == False:
        lon1 = np.deg2rad(lon1)
        lat1 = np.deg2rad(lat1)

    lon2, lat2 = np.cos(lon1+orthodrome), np.sin(lat1+orthodrome)

    return lon2,lat2

# ----------------------------------------------------------------------------------------

def arcsec_to_carrington(arc_x, arc_y, time):

    from astropy.coordinates import SkyCoord#sky_coordinate
    import astropy.units as u
    from sunpy.coordinates import frames
    from sunpy.coordinates import get_horizons_coord

    """
    transforms arcsec on the solar disk (as seen from Earth) to Carrington longitude & latitude    Parameters
    ----------
    arc_x : float
        Helioprojective x coordinate in arcsec
    arc_y : float
        Helioprojective y coordinate in arcsec
    time : string
        date and time, for example: '2021-04-17 16:00'

    Returns
    -------
    lon, lat : float 
        longitude and latitude in Carrington coordinates
    """

    c = SkyCoord(arc_x*u.arcsec, arc_y*u.arcsec, frame=frames.Helioprojective, obstime=time, observer="earth")
    Carr = c.transform_to(frames.HeliographicCarrington(observer="Sun"))
    lon = Carr.lon.value
    lat = Carr.lat.value

    return lon, lat

# ----------------------------------------------------------------------------------------

def get_pfss_hmimap(filepath, email, carrington_rot, date, rss=2.5, nrho=35):
    '''
    downloading hmi map or calculating the PFSS solution
    '''

    time = a.Time(date, date)
    pfname =  filepath+"PFSS_output_"+ str(time.start.datetime.date())+'_CR'+str(carrington_rot)+'_SS'+str(rss)+'_nrho'+str(nrho)+'.p'

    #check if PFSS file already exists locally:
    try:
        with open(pfname, 'rb') as f:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            output = u.load()
        print('Found pickled PFSS file')

    #if not, then download MHI mag, calc. PFSS, and save as picle file for next time
    except FileNotFoundError:
        series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')
        crot = a.jsoc.PrimeKey('CAR_ROT', carrington_rot)
        result = Fido.search(time, series, crot, a.jsoc.Notify(email))
        files = Fido.fetch(result)
        hmi_map = sunpy.map.Map(files[0])
        pfsspy.utils.fix_hmi_meta(hmi_map)
        print('Data shape: ', hmi_map.data.shape)

        hmi_map = hmi_map.resample([360, 180]*units.pix)
        print('New shape: ', hmi_map.data.shape)

        pfss_input = pfsspy.Input(hmi_map, nrho, rss)
        output = pfsspy.pfss(pfss_input)
        with open(pfname, 'wb') as f:
            pickle.dump(output, f)

    return output

# ----------------------------------------------------------------------------------------

def circle_around(x,y,n,r=0.1):
    '''
    Produces new points around a (x,y) point in a circle.
    
    x,y: coordinates of the original point
    n: the amount of new points around the origin
    r: the radius of the circle at which new points are placed (in radians)
    
    returns:
    pointlist: list of new points (tuples) around the original point in a circle, placed
                at equal intervals
                
    At the moment does not work perfectly in the immediate vicinity of either pole.
    '''
 
    origin = (x,y)

    x_coords = np.array([])
    y_coords = np.array([])
    for i in range(0,n):

        theta = (2*i*np.pi)/n
        newx = origin[0] + r*np.cos(theta)
        newy = origin[1] + r*np.sin(theta)

        if newx >= 2*np.pi:
            newx = newx - 2*np.pi

        if newy > np.pi/2:
            overflow = newy - np.pi/2
            newy = newy - 2*overflow

        if newy < -np.pi/2:
            overflow = newy + np.pi/2
            newy = newy + 2*overflow

        x_coords = np.append(x_coords,newx)
        y_coords = np.append(y_coords,newy)

    pointlist = np.array([x_coords,y_coords])

    return pointlist

# ----------------------------------------------------------------------------------------

def vary_flines(lon, lat, hmimap, n_varies):
    '''
    Finds a set of sub-pfss fieldlines connected to or very near a single footpoint on the pfss.
    
    lon: longitude of the footpoint [rad]
    lat: latitude of the footpoint [rad]
    
    n_varies:   tuple that holds the amount of circles and the number of dummy flines per circle
                if type(n_varies)=int, consider that as the amount of circles, and set the 
                amount of dummy flines per circle to 16

    n_circles:  the amount of circles of fieldlines traced
    n_flines:   the number of dummy fieldlines per one circle of fieldlines
    '''
    #field lines per n_circles (circle)
    if isinstance(n_varies,list):
        n_circles = n_varies[0]
        n_flines = n_varies[1]
    else:
        n_circles = n_varies
        n_flines = 16

    #first produce new points around the given lonlat_pair
    lons,lats= np.array([lon]), np.array([lat])
    increments = np.array([0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29])
    for circle in range(n_circles):

        newlons,newlats = circle_around(lon,lat,n_flines,r=increments[circle])
        lons, lats = np.append(lons,newlons), np.append(lats,newlats)


    pointlist = np.array([lons,lats])

    #trace fieldlines from all of these points
    varycoords, varyflines = get_field_line_coords(pointlist[0],pointlist[1],hmimap)

    #because the original fieldlines and the varied ones are all in the same arrays,
    #extract the varied ones to their own arrays
    coordlist, flinelist = [], []

    #total amount of flines = 1 + (circles) * (fieldlines_per_circle)
    total_per_fp = n_flines*n_circles+1
    erased_indices = []
    for i in range(len(varycoords)):
        #n_flines*n_circles = the amount of extra field lines between each "original" field line
        if i%(total_per_fp)==0:
            erased_indices.append(i)
            #pop(i) removes the ith element from the list and returns it
            #-> we append it to the list of original footpoint fieldlines
            coordlist.append(varycoords[i]) #.pop(i)
            flinelist.append(varyflines[i])

    #really ugly quick fix to erase values from varycoords and varyflines
    for increment, index in enumerate(erased_indices):
        varycoords.pop(index-increment)
        varyflines.pop(index-increment)

    return coordlist, flinelist, varycoords, varyflines

# ----------------------------------------------------------------------------------------

def get_coord_values(field_line):
    '''
    Gets the coordinate values from FieldLine object and makes sure that they are in the right order.
    '''
    
    #first check that the field_line object is oriented correctly (start on photosphere and end at pfss)
    fl_coordinates = field_line.coords
    fl_coordinates = check_field_line_alignment(fl_coordinates)
    
    fl_r = fl_coordinates.radius.value / const.R_sun.value
    fl_lon = fl_coordinates.lon.value
    fl_lat = fl_coordinates.lat.value
    
        
    return fl_r, fl_lon, fl_lat

# ----------------------------------------------------------------------------------------

def get_field_line_coords(longitude, latitude, hmimap):
    '''
    Returns triplets of open magnetic field line coordinates, and the field line object itself
    
    longitude and latitude are given in radians
    '''

    #the amount of coordinate triplets we are going to trace
    try: 
        coord_triplets = len(latitude)
    except TypeError:
        coord_triplets = 1
        latitude = [latitude]
        longitude = [longitude]

    #the loop in which we trace the field lines and collect them to the coordlist
    coordlist = []
    flinelist = []
    for i in range(coord_triplets):

        increment = 1
        sign_switch = 1
        init_lat0 = latitude[i]
        init_lon0 = longitude[i]

        #keep tracing the field line until a valid one is found
        while(True):

            #trace a field line downward from the point lon,lat on the pfss
            fline = trace_field_line(longitude[i], latitude[i], hmimap)

            radius0 = fline.coords.radius[0].value
            radius9 = fline.coords.radius[-1].value
            bool_key = (radius0==radius9)

            #if fline is not a valid field line, then alter lat a little and try again
            #also check if this is a null line (all coordinates identical)
            if( (len(fline.coords) < 10) or bool_key ):

                latitude[i] = init_lat0 + increment*sign_switch*0.0001
                sign_switch = sign_switch*(-1)
                if(sign_switch>0):
                    increment += 1

            #skip closed lines
            elif fline.polarity == 0:

                longitude[i] = init_lon0 + increment*sign_switch*0.0001
                sign_switch = sign_switch*(-1)
                if(sign_switch>0):
                    increment += 1

            #check that we are not too far from the original coordinate
            elif(increment > 500):
                raise Exception('no field line footpoint found on the given coordinate')

            #if there was nothing wrong, break the loop and proceed with the traced field line
            else:
                #print("increment:",increment)
                break

        #get the field line coordinate values in the correct order
        #start on the photopshere, end at the pfss
        fl_r, fl_lon, fl_lat   = get_coord_values(fline)

        #fill in the lists
        triplet = [fl_r, fl_lon, fl_lat]
        coordlist.append(triplet)
        flinelist.append(fline)


    return coordlist, flinelist

# ----------------------------------------------------------------------------------------

def au_to_km(distlist):
    
    for i in range(len(distlist)):
        distlist[i] = distlist[i]*(150e6)

    return distlist

# ----------------------------------------------------------------------------------------

def multicolorline(x, y, cvals, ax, vmin=-90, vmax=90):
    '''
    original example from:
    https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
    '''

    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(vmin, vmax)

    cmrmap = cmr.redshift

    #sample the colormaps that you want to use. Use 90 from each so there is one
    #color for each degree
    colors_pos = cmrmap(np.linspace(0.0, 0.30, 45))
    colors_neg = cmrmap(np.linspace(0.70, 1.0, 45))

    #combine them and build a new colormap
    colors = np.vstack((colors_pos, colors_neg))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    #establish the linecollection object
    lc = LineCollection(segments, cmap=mymap, norm=norm)

    #set the values used for colormapping
    lc.set_array(cvals)

    #set the width of line 
    lc.set_linewidth(3)

    #this we want to return 
    line = ax.add_collection(lc)

    return line

# ----------------------------------------------------------------------------------------

def plot3d(field_lines):
    
    if not isinstance(field_lines, list):
        field_lines = [field_lines]
        
    if isinstance(field_lines[0],list):
        field_lines = flatten(field_lines)

    fig, axarr = plt.subplots(subplot_kw={"projection": "3d"})

    axarr.set_box_aspect((1, 1, 1))

    r_sun, r_ss = 1.0, 2.5
    
    #Draw the Sun
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    axarr.plot_wireframe(x*r_sun, y*r_sun, z*r_sun, color="darkorange")
    axarr.set_xlim(-2,2)
    axarr.set_ylim(-2,2)
    axarr.set_zlim(-2,2)
    
    for field_line in field_lines:
            coords = field_line.coords
            coords.representation = 'cartesian'
            color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
            axarr.plot(coords.x / const.R_sun,
            coords.y / const.R_sun,
            coords.z / const.R_sun,
            color=color, linewidth=1)
    
    try:
        axarr.set_aspect('equal', adjustable='box')
    except NotImplementedError:
        axarr.set_aspect('auto', adjustable='box')

# ----------------------------------------------------------------------------------------

def draw_fieldlines(field_lines, frame='yz', save=False):

    import matplotlib.patches as mpatch
    
    #check if there's a list inside a list, if there is -> flatten
    if isinstance(field_lines[0],list):
        field_lines = flatten(field_lines)

    fig, ax = plt.subplots(figsize=[10,10])
    r_ss = 2.5

    ax.set_aspect('equal')

    if(isinstance(field_lines,list)):

        for field_line in field_lines:
            coords = field_line.coords
            coords.representation = 'cartesian'
            color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
            if(frame=='yz'):
                ax.plot(coords.cartesian.y / const.R_sun, coords.cartesian.z / const.R_sun, color=color)
                projection = 'POV: Carrington longitude 0'
            elif(frame=='xy'):
                ax.plot(coords.cartesian.x / const.R_sun, coords.cartesian.y / const.R_sun, color=color)
                projection = 'POV: North'
            elif(frame=='xz'):
                ax.plot(coords.cartesian.x / const.R_sun, coords.cartesian.z / const.R_sun, color=color)
                projection = 'POV: Carrington longitude 270'
            else:
                raise Exception("Invalid frame")

    else:

        field_line = field_lines
        coords = field_line.coords
        coords.representation = 'cartesian'
        color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
        if(frame=='yz'):
            ax.plot(coords.cartesian.y / const.R_sun, coords.cartesian.z / const.R_sun, color=color)
            projection = 'POV: Carrington longitude 0'
        elif(frame=='xy'):
            ax.plot(coords.cartesian.x / const.R_sun, coords.cartesian.y / const.R_sun, color=color)
            projection = 'POV: North'
        elif(frame=='xz'):
            ax.plot(coords.cartesian.x / const.R_sun, coords.cartesian.z / const.R_sun, color=color)
            projection = 'POV: Carrington longitude 270'
        else:
            raise Exception("Invalid frame")

    ax.add_patch(mpatch.Circle((0,0), 1, color='darkorange', lw=2.0, fill=False))
    ax.add_patch(mpatch.Circle((0,0), r_ss, color='k', linestyle='--', fill=False))

    plt.title(projection)

    if save:
        plt.savefig('overhead.png')

    plt.show()

# ----------------------------------------------------------------------------------------

def flatten(l):
    '''
    Flattens a list of lists to a single list.
    '''
    flat = []
    for item in l:
        try:
            for subitem in item:
                flat.append(subitem)
        except TypeError:
            return l

    return flat

# ----------------------------------------------------------------------------------------

def check_field_line_alignment(coordinates):
    '''
    Checks that a field line object is oriented such that it starts from
    the photpshere and ends at the pfss. If that is not the case, then
    flips the field line coordinates over and returns the flipped object.
    '''
    
    fl_r = coordinates.radius.value
    
    if fl_r[0] > fl_r[-1]:
        coordinates = np.flip(coordinates)

    return coordinates

# ----------------------------------------------------------------------------------------

def trace_field_line(lon0, lat0, hmimap, rad=True):
    '''
    Traces a single open magnetic field line at coordinates (lon0,lat0) on the pfss down
    to the photosphere
    '''

    #if lat0 and lon0 are given in deg for some reason, transform them to rad
    if not rad:
        lat0 = np.deg2rad(lat0)
        lon0 = np.deg2rad(lon0)

    #start tracing from the pfss height
    height = 2.5*const.R_sun
    tracer = tracing.PythonTracer()

    #add unit to longitude and latitude, so that SkyCoord understands them
    lon, lat = lon0*units.rad, lat0*units.rad

    #seed the starting coordinate at the desired coordinates
    seed = SkyCoord(lon, lat, height, frame=hmimap.coordinate_frame)

    #trace the field line from the seed point given the hmi map
    field_lines = tracer.trace(seed, hmimap)

    #field_lines is a list of len=1, because there's only one seed point given to the tracer
    field_line = field_lines[0]

    return field_line

# ----------------------------------------------------------------------------------------

def parker_spiral(sw, distance, longitude, resolution, endpoint=2.5, backtrack=True):
    '''
    construct one magnetic parker spiral arm

    INPUT
    sw: solar wind speed in km/s
    distance: distance to the object in km
    longitude: angular coordinate of the object (stellar longitude?) in deg
    resolution: resolution of the curve (amount of points making the curve)
    endpoint: the point at which one wants to stop tracing back the spiral arm in solar radii

    RETURNS
    phi: array of angular coordinates in rad
    r: array of radial coordinates in solar radii
    '''
    #parker spiral solution e.g. here:
    #http://www.physics.usyd.edu.au/~cairns/teaching/2010/lecture8_2010.pdf

    omega = 2.694e-6 #rad/s

    r = np.linspace(endpoint, distance, resolution) #linspace to get even spacing

    #backtracking means going from sc to source surface
    if backtrack:
        phi = longitude + (omega)*(distance-r)/sw
    else:
        phi = longitude + (omega*r)/sw

    return phi, r

# ----------------------------------------------------------------------------------------

def symlog_pspiral(sw, distance, longitude, latitude, hmimap, names=None, title='', r_s=2.5, \
                    vary=False, n_varies=1, save=False):
    '''
    Produces a figure of the heliosphere in polar coordinates with logarithmic r-axis outside the pfss.
    Also tracks an open field line down to photosphere given a point on the pfss.

    sw = solar wind in km/s, type: int, float or list
    distance = distance to the object in km, type: int, float or list
    longitude = angular coordinate of the object in deg (stellar longitude?), type: int, float or list
    latitude = see longitude
    r_s = source surface height of the potential field, type: float
    resolution = resolution of the spiral curve, type: int
    save = boolean value to save the figure
    
    --RETURNS---
    either,
        fline_objects: a list holding all the field line objects that were calculated for the plot
    or if vary,
        [fline objects, varyfline_objects]: a list containing a list of field line objects, and also
                                            lists of all the varied field line objects
    
    '''

    if names == None:
        names = [None]*len(sw)

    sun_radius = const.R_sun.value  #/ 10**6 #m

    #treat lons and lats as radians in the function
    longitude = np.deg2rad(longitude)
    latitude = np.deg2rad(latitude)

    #normalize variables to solar radii
    if isinstance(sw,list):
        sw_norm = [(u*1000)/sun_radius for u in sw]
        distance_norm = [(d*1000)/sun_radius for d in distance]

    else:
        sw_norm = (sw*1000)/sun_radius
        distance_norm =(distance*1000)/sun_radius

    #projection of the objects on the plane of ecliptic
    projection = np.cos(latitude)

    #calculate parker spiral for given objects
    if isinstance(sw,list):
        phis, rs = [], []
        for i in range(len(sw)):
            phi, r = parker_spiral(sw_norm[i], distance_norm[i], longitude[i], resolution=1000)
            phis.append(phi)
            rs.append(r)

        sc_footphis = [phi[0] for phi in phis]
        sc_footrs = [r[0] for r in rs]
        sc_footpoint = [sc_footphis,sc_footrs]

    else:
        phi, r = parker_spiral(sw_norm, distance_norm, longitude, resolution=1000)
        sc_footpoint = [phi[0],r[0]]

    #----------------------------------------------------------
    #tracing the closest field line to sc_footpoint down to photosphere:

    #acquire an array of (r,lon,lat) coordinates of the open field lines under the pfss
    #based on the footpoint(s) of the sc
    if vary:
        #if there is more than one objects being both traced and varied, we have to vary them one at a time
        if len(sc_footpoint[0]) > 1:
            fline_triplets, fline_objects, varyfline_triplets, varyfline_objects = [], [], [], []
            for i, footpoint in enumerate(sc_footpoint[0]):
                #Append doesn't work here, but a simple + does. I wonder why.
                tmp_triplets, tmp_objects, varytmp_triplets, varytmp_objects = vary_flines(footpoint, latitude[i], hmimap, n_varies)
                fline_triplets = fline_triplets + tmp_triplets
                fline_objects = fline_objects + tmp_objects
                varyfline_triplets = varyfline_triplets + varytmp_triplets
                varyfline_objects = varyfline_objects + varytmp_objects

        #if only a single object, then just run vary_flines for it
        else:
            fline_triplets, fline_objects, varyfline_triplets, varyfline_objects = vary_flines(sc_footpoint[0], latitude, hmimap, n_varies)


    #if no varying, then just get one field line from get_field_line_coords()
    else:
        fline_triplets, fline_objects = get_field_line_coords(sc_footpoint[0], latitude, hmimap)

    #we need fl_r, fl_lon, fl_lat
    #they are located in fline_triplets[i][j]

    #source surface:
    theta = 2*np.pi*np.linspace(0,1,200)
    ss = np.ones(200)*r_s

    #Plotting commands------------------------------------------------------------->
    fig, ax = plt.subplots(figsize = [19,17], subplot_kw={'projection': 'polar'})

    #plot the source_surface and solar surface
    ax.plot(theta, ss, c='k', ls='--', zorder=1)
    ax.plot(theta, np.ones(200), c='darkorange', lw=2.5, zorder=1)
    
    #plot the 30 and 60 deg lines on the Sun
    ax.plot(theta, np.ones(len(theta))*0.866, c='darkgray', lw=0.5, zorder=1) #cos(30deg) = 0.866(O)
    ax.plot(theta, np.ones(len(theta))*0.500, c='darkgray', lw=0.5, zorder=1) #cos(60deg) = 0.5

    #plot the spiral
    if isinstance(sw,list):
        for i in range(len(sw)):
            ax.plot(phis[i], projection[i]*rs[i], c=get_color(names[i]), label='sw={} km/s'.format(int(sw[i])), zorder=1)
            ax.scatter(sc_footpoint[0][i], projection[i]*sc_footpoint[1][i], s=115, c=get_color(names[i]), marker='D', zorder=3)
            ax.scatter(phis[i][-1], projection[i]*rs[i][-1], c=get_color(names[i]), alpha=1.0, s=175, zorder=3)
    else:
        ax.plot(phi, projection*r, c='k', label='sw={} km/s'.format(int(sw)), zorder=1)
        #mark the footpoint of the spiral arm on the source surface
        ax.scatter(sc_footpoint[0], projection*sc_footpoint[1], s=140, c='C0', marker='x', zorder=3)
        #mark the object
        ax.scatter(phi[-1], projection*r[-1], c='C0', s=175, zorder=3)


    #plot the field line(s) connecting ss_footpoint and the solar surface and collect relevant 
    #points (footpoints) to an array
    display_fl_footpoints = []
    display_fl_sourcepoints = []
    for fline in fline_triplets:
        
        fl_r   = fline[0]
        fl_lon = fline[1]
        fl_lat = fline[2]
        
        #remember the foot and source points
        display_fl_sourcepoints.append([np.round(fl_r[0],1), np.round(fl_lon[0],1), np.round(fl_lat[0],1)])
        display_fl_footpoints.append([np.round(fl_r[-1],1), np.round(fl_lon[-1],1), np.round(fl_lat[-1],1)])
        
        #plot the color coded field line
        fieldline = multicolorline(np.deg2rad(fl_lon), np.cos(np.deg2rad(fl_lat))*fl_r, ax=ax, cvals=fl_lat, vmin=-90, vmax=90)

    if vary:
        for fline in varyfline_triplets:
        
            fl_r   = fline[0]
            fl_lon = fline[1]
            fl_lat = fline[2]

            #plot the color coded varied field line
            fieldline = multicolorline(np.deg2rad(fl_lon), np.cos(np.deg2rad(fl_lat))*fl_r, ax=ax, cvals=fl_lat, vmin=-90, vmax=90)

    #Settings--------------------------------->

    r_max = 500e9 / sun_radius
    ax.set_rmax(r_max)
    ax.set_rscale('symlog', linthresh=r_s)

    ax.grid(True)
    #ax.set_thetamin(0)
    #ax.set_thetamax(135)
    ax.set_rticks([1.0, 2.5, 10.0, 100.0])
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    rlabels = ['1', '2.5', r'$10^1$', r'$10^2$']
    ax.set_yticklabels(rlabels)

    #plt.legend(loc=10, bbox_to_anchor=(-0.1, 0.1))

    #create the colorbar displaying values of the last fieldline plotted
    cb = fig.colorbar(fieldline)
    
    #@TODO: shrink the colorbar andmove it to the top right corner
    
    #colorbar is the last object created -> it is the final index in the list of axes
    cb_ax = fig.axes[-1]
    cb_ax.set_ylabel('latitude [deg]')
    
    #before adding txtboxes, make sure that sc_footpoint is of proper shape
    if(isinstance(sc_footpoint[0],float)):
        display_sc_footpoint = [[sc_footpoint[0]],[sc_footpoint[1]]]
        latitude = [latitude]

    #also add magnetic polarity info to display:
    display_polarities = []
    for fline in fline_objects:
        display_polarities.append(int(fline.polarity))
    
    
    #txtbox stuff----------->
    txtsize = 12
    
    #Make legend for abbreviations:
    legendlabel = "ABBREVIATIONS \nPS: Photosphere \nSS: Source Surface \nFP: Footpoint \nP: Polarity"
    legendbox = AnchoredText(f"{legendlabel}", prop=dict(size=txtsize+1), frameon=True, loc=(10), bbox_to_anchor=(92,880)) #float('1.{}'.format(i))
    legendbox.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
    legendbox.patch.set_linewidth(3.0)
    legendbox.patch.set_edgecolor(get_color(None))
    ax.add_artist(legendbox)
    
    #add textbox for footpoint coordinates
    if(names is not None):
        for i in range(len(display_fl_footpoints)):
            plabel = AnchoredText(f"{names[i]}, sw={int(sw[i])} km/s\nPS FP = ({display_fl_sourcepoints[i][1]},{display_fl_sourcepoints[i][2]}) \nSS FP = ({display_fl_footpoints[i][1]},{display_fl_footpoints[i][2]}) \nP: {display_polarities[i]}", 
                                  prop=dict(size=txtsize+1), frameon=True, loc=(2), bbox_to_anchor=(15, (840-i*75) )) #float('1.{}'.format(i))
            plabel.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            plabel.patch.set_linewidth(3.0)
            plabel.patch.set_edgecolor(get_color(names[i]))
            ax.add_artist(plabel)
    else:
        for i in range(len(display_fl_footpoints)):
            plabel = AnchoredText(f"PS FP = ({display_fl_sourcepoints[i][1]},{display_fl_sourcepoints[i][2]}) \n SS FP = ({display_fl_footpoints[i][1]},{display_fl_footpoints[i][2]}) \n SPACECRAFT FOOTPOINT \n ({np.round(np.rad2deg(sc_footpoint[0][i]),1)},{np.round(latitude[i],1)})", 
                                  prop=dict(size=txtsize+1), frameon=True, loc=(10), bbox_to_anchor=(100, (800-i*75) )) #float('1.{}'.format(i))
            plabel.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            plabel.patch.set_linewidth(1.5)
            plabel.patch.set_edgecolor(get_color(names[i]))
            ax.add_artist(plabel)

    #set the title of the figure
    plt.title(title)

    if(save):
        plt.savefig('testfig.png', transparent=False, facecolor='white', bbox_inches='tight')

    plt.show()

    #the function returns all the calculated field line objects, which include many attributes
    #of the field lines such as coordinates, polarity, and wether they are open or closed
    if vary:
        return [fline_objects, varyfline_objects]
    else:
        return fline_objects

# ---------------------------------------------------------------------------

def field_line_info_to_df(flines, names):
    '''
    Takes a list of field line objects and names, and assembles
    a pd dataframe that includes the footpoint and head of each magnetic
    field line and polarity corresponding to each sc.
    '''

    if not isinstance(names, np.ndarray):
        names = np.array(names)

    if len(names) != len(flines):

        for name in names:
            if(name is not None):
                raise Exception("The list of names does not match the list of field line objects.")

        names = np.array([None]*len(flines))


    #init arrays to store into dataframe
    polarities = np.array([])
    photosphere_lons = np.array([])
    photosphere_lats = np.array([])
    pfss_lons = np.array([])
    pfss_lats = np.array([])

    for line in flines:

        polarities = np.append(polarities,line.polarity)
        
        coordinates = check_field_line_alignment(line.coords)

        photospheric_footpoint = (coordinates.lon.value[0], coordinates.lat.value[0])
        photosphere_lons = np.append(photosphere_lons,photospheric_footpoint[0])
        photosphere_lats = np.append(photosphere_lats,photospheric_footpoint[1])

        pfss_footpoint = (coordinates.lon.value[-1], coordinates.lat.value[-1])
        pfss_lons = np.append(pfss_lons, pfss_footpoint[0])
        pfss_lats = np.append(pfss_lats, pfss_footpoint[1])
        

    data_dict = {'names': names,
                'footpoint lon': photosphere_lons,
                'footpoint lat': photosphere_lats,
                'pfss lon': pfss_lons,
                'pfss lat': pfss_lats,
                'polarity': polarities}
    
    df = pd.DataFrame(data=data_dict)
    
    return df
    
# ---------------------------------------------------------------------------

def df_to_file(df, filename: str):
    '''
    Writes a dataframe in to a csv file
    '''
    
    if not isinstance(filename,str):
        raise Exception("The file name must be string.")
    
    current_directory = os.getcwd()
    
    filestr = f"{current_directory}/{filename}.csv"
    
    df.to_csv(filestr)
    print(f"Created file {filename}.csv to {current_directory}/")

# ---------------------------------------------------------------------------

def write_info_to_csv(flines, names=[None], filename='magnetic_info'):
    
    #first assemble a dataframe from the field line object
    coords_df = field_line_info_to_df(flines, names)
    
    #write the df into csv
    df_to_file(coords_df, filename)

# ---------------------------------------------------------------------------

def get_sc_data(csvfile: str):
    '''
    Reads the contents of solar-mach produced csv file, and returns lists
    of necessary data to run pfss field line tracing analysis.
    
    csvfile: str, the name of the csv file one wants to read
    '''
    import pandas as pd

    if type(csvfile) is not str:
        raise TypeError("File name is not a string.")

    csvdata = pd.read_csv(csvfile)
    
    names = list(csvdata['Spacecraft/Body'])
    lons = list(csvdata['Carrington Longitude (°)'])
    lats = list(csvdata['Latitude (°)'])
    dist = au_to_km(list(csvdata['Heliocentric Distance (AU)']))
    sw = list(csvdata['Vsw'])
    
    return names, sw, dist, lons, lats

#============================================================================================================================

#============================================================================================================================
#this will be ran when importing
#set Jupyter notebook cells to 100% screen size:
# from IPython.core.display import display, HTML
# display(HTML(data="""<style> div#notebook-container { width: 99%; } div#menubar-container { width: 85%; } div#maintoolbar-container { width: 99%; } </style>"""))
