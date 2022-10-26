'''
A library intended for importing into pfsspy_notebook.ibynp. It contains the necessary functions
for seeking the footpoints of IMF field lines connecting back to the photosphere and plotting them.

@Author: Christian Palmroos
         <chospa@utu.fi>

Last updated: 2022-10-26
'''

# imports:
import os
import numpy as np
import pandas as pd

import astropy.constants as const
import astropy.units as units
# import astropy.coordinates
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatch
import cmasher as cmr
import pfsspy
import pickle
import warnings
import sunpy.map

from matplotlib.collections import LineCollection

# from astropy.time import Time
from astropy.coordinates import SkyCoord
from pfsspy import tracing
# from pfsspy.sample_data import get_gong_map
# from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
from sunpy.net import Fido, attrs as a
from sunpy.coordinates import frames

# some matplotlib settings:
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.size'] = 20
plt.rcParams['agg.path.chunksize'] = 20000

# living on the edge:
warnings.simplefilter('ignore')

# ----------------------------------------------------------------------------------------


def save_pickle(obj, file_name):
    """
    Saves an object as a python pickle file.
    
    Parameters:
    -----------
    obj : object
    file_name : str
    """

    with open(file_name, 'wb') as file:
        pickle.dump(obj, file)


def load_pickle(file_name):
    """
    Loads a python pickle and returns it.

    Parameters:
    -----------
    file_name : str

    Returns:
    ----------
    obj : object
    """

    with open(file_name, 'rb') as file:
        obj = pickle.load(file)

    if isinstance(obj,list):
        obj = np.asarray(obj, dtype=object)

    return obj


def get_color(key: str = None) -> str:
    """
    Returns the right color for an object according to SOLAR-MACH tool.

    params
    -------
    key: str (default = None)
            The name of the spacecraft or planet.

    returns
    -------
    color: str (default = 'k')
            The color identifier that matplotlib understands
    """

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

    try:
        return color_dict[key]
    except KeyError:
        return color_dict['default']


def get_sc_data(csvfile: str):
    """
    Reads the contents of solar-mach produced csv file, and returns lists
    of necessary data to run pfss field line tracing analysis.

    params
    ------
    csvfile: str
            The name of the csv file one wants to read
    
    returns
    -------
    names: list[str]
            List of names
    sw: list[int/float]
            List of solar winds in units of km/s
    dist: list[int/float]
            List of heliocentric distances in units of AU
    lons: list[float]
            List of Carrington longitudes of the corresponding objects
    lats: list[float]
            List of Carrington latitudes of the corresponding objects
    """

    if type(csvfile) is not str:
        raise TypeError("File name is not a string.")

    csvdata = pd.read_csv(csvfile)

    names = list(csvdata['Spacecraft/Body'])
    lons = list(csvdata['Carrington Longitude (°)'])

    # for some reason the output of Streamlit Solar-MACH is different from notebook-produced csv
    try:
        lats = list(csvdata['Latitude (°)'])
    except:
        lats = list(csvdata['Carrington Latitude (°)'])

    dist = au_to_km(list(csvdata['Heliocentric Distance (AU)']))
    sw = list(csvdata['Vsw'])

    return names, sw, dist, lons, lats


def field_line_accuracy(flines):
    """
    Calculates at least now the central point, average distance from the central point and
    the standard deviation of the photospheric footpoints for a set of field lines
    """

    if isinstance(flines[0], list):
        flines = flatten(flines)

    footpoints = []
    for fline in flines:

        r, lon, lat = get_coord_values(fline)

        # Index 0 of coordinates corresponds to photospheric coordinates, index -1 to pfss
        footpoint = [lon[0], lat[0]]
        footpoints.append(footpoint)

    # Declare longitudes and latitudes of footpoints
    footp_lons = [pair[0] for pair in footpoints]
    footp_lats = [pair[1] for pair in footpoints]

    # If separation in longitudes is over half a circle, then points are probably
    # scattered near the 0-deg line -> transfer them far enough from that line
    # while calculations are ran
    if max(footp_lons) > min(footp_lons)+180.0:

        # Transfer of longitudes
        footp_lons = shift_longitudes(footp_lons)

        # Calculate the central point
        c_point = [np.mean(footp_lons), np.mean(footp_lats)]

        # Standard deviation of longitudes and latitudes
        lon_std = np.std(footp_lons)
        lat_std = np.std(footp_lats)

        # Calculate mean distance from the central point
        dist_sum = 0
        for i in range(len(footp_lons)):

            lon1, lat1 = footp_lons[i], footp_lats[i]
            angular_separation = orthodrome(lon1, lat1, c_point[0], c_point[1])

            # Distance is in solar radii
            distance_rs = angular_separation
            dist_sum = dist_sum + distance_rs

        # Transfer lons and the central point back the same amount that the longitudes were moved
        footp_lons = shift_longitudes(footp_lons, shift=-180.0)
        c_point[0] = shift_longitudes([c_point[0]], shift=-180.0)[0]

    else:

        # Calculate the central point
        c_point = [np.mean(footp_lons), np.mean(footp_lats)]

        # Standard deviation of longitudes and latitudes
        lon_std = np.std(footp_lons)
        lat_std = np.std(footp_lats)

        # Calculate mean distance from the central point
        dist_sum = 0
        for i in range(len(footp_lons)):

            lon1, lat1 = footp_lons[i], footp_lats[i]
            angular_separation = orthodrome(lon1, lat1, c_point[0], c_point[1])

            # Distance is in solar radii
            distance_rs = angular_separation
            dist_sum = dist_sum + distance_rs

    avg_dist = dist_sum/len(footp_lons)

    return footpoints, c_point, avg_dist, [lon_std, lat_std]


def shift_longitudes(lons, shift=180.0):
    """
    Shifts the longitudes by <shift> amount
    """

    if shift>0:
        lons = [lon+shift for lon in lons]
        lons = [lon-360.0 if lon>360 else lon for lon in lons]
    else:
        # If shift is negative, points are likely being moved back to their
        # original place
        lons = [lon+shift for lon in lons]
        lons = [lon+360.0 if lon<0 else lon for lon in lons]

    return lons


def map_on_surface(fps, c_point, avg_d, shift=None, zoom=None, show_avg_d=False):
    """
    Plots a really simple 2d representation of fieldline objects' footpoints.
    """

    import matplotlib.patches as mpatch

    centre = np.array(c_point)
    fpslons = [item[0] for item in fps]
    fpslats = [item[1] for item in fps]

    if shift is not None:
        fpslons = shift_longitudes(fpslons, shift=shift)
        centre[0] = c_point[0]+shift

    fig_tuple = (16, 7)
    if zoom is not None:
        fig_tuple = (12, 12)

    fig = plt.figure(figsize=fig_tuple)
    ax = plt.subplot()

    ax.scatter(fpslons[0], fpslats[0], c='navy', s=60, label="original footpoint")
    ax.scatter(fpslons[1:], fpslats[1:], c='C0', label="dummy footpoints")
    ax.scatter(centre[0], centre[1], c='r', label="avg(lons,lats)")

    if show_avg_d:
        avg_d_deg = np.rad2deg(avg_d)
        ax.add_patch(mpatch.Circle((centre[0], centre[1]), avg_d_deg, color='r', lw=0.8, ls='--', fill=False))

    plt.ylim(-90, 90)
    plt.xlim(0, 360)
    if zoom is not None:
        plt.ylim(centre[1]-zoom/2, centre[1]+zoom/2)
        plt.xlim(centre[0]-zoom/2, centre[0]+zoom/2)

    plt.legend()
    plt.grid("True")
    plt.show()


def orthodrome(lon1, lat1, lon2, lat2, rad=False) -> float:
    """
    Calculates the othodrome (total angular separtion) between two coordinates defined by their lon/lat positions

    params
    -------
    lon1: int/float
            Longitude of the first coordinate point
    lat1: int/float
            Latitude of the first coordinate point
    lon2: int/float
            See lon1
    lat2: int/float
            See lat1
    rad: bool (default = False)
            Is the input given as radians? If not, treat them as degrees
    
    returns
    -------
    ortho: float
            The total angular separation between (lon1,lat1) and (lon2,lat2)
    """

    if(rad == False):
        lon1 = np.deg2rad(lon1)
        lon2 = np.deg2rad(lon2)
        lat1 = np.deg2rad(lat1)
        lat2 = np.deg2rad(lat2)

    ortho = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))

    return ortho


def ortho_to_points(lon1, lat1, orthodrome, rad=False):
    """
    Calculates a lon/lat pair from a reference point given orthodrome away.
    """

    if rad == False:
        lon1 = np.deg2rad(lon1)
        lat1 = np.deg2rad(lat1)

    lon2, lat2 = np.cos(lon1+orthodrome), np.sin(lat1+orthodrome)

    return lon2, lat2


def _isstreamlit():
    """
    Function to check whether python code is run within streamlit

    Returns
    -------
    use_streamlit : boolean
        True if code is run within streamlit, else False
    """
    # https://discuss.streamlit.io/t/how-to-check-if-code-is-run-inside-streamlit-and-not-e-g-ipython/23439
    try:
        from streamlit.scriptrunner import get_script_run_ctx
        if not get_script_run_ctx():
            use_streamlit = False
        else:
            use_streamlit = True
    except ModuleNotFoundError:
        use_streamlit = False
    return use_streamlit


def null_decorator(a):
    # decorator that does nothing
    return a


# if run inside streamlit, use streamlit's caching decorator on get_pfss_hmimap()
if _isstreamlit():
    import streamlit as st
    st_cache_decorator = st.cache(persist=True, allow_output_mutation=True)
else:
    st_cache_decorator = null_decorator


@st_cache_decorator
def get_pfss_hmimap(filepath, email, carrington_rot, date, rss=2.5, nrho=35):
    """
    Downloading hmi map or calculating the PFSS solution

    params
    -------
    filepath: str
            Path to the hmimap, if exists.
    email: str
            The email address of a registered user
    carrington_rot: int
            The Carrington rotation corresponding to the hmi map
    date: str
            The date of the map. Format = 'YYYY/MM/DD'
    rss: float (default = 2.5)
            The height of the potential field source surface for the solution.
    nrho: int (default = 35)
            The resolution of the PFSS-solved field line objects
    
    returns
    -------
    output: hmi_synoptic_map object
            The PFSS-solution
    """

    time = a.Time(date, date)
    pfname =  f"PFSS_output_{str(time.start.datetime.date())}_CR{str(carrington_rot)}_SS{str(rss)}_nrho{str(nrho)}.p"

    # Check if PFSS file already exists locally:
    print(f"Searching for PFSS file from {filepath}")
    try:
        with open(f"{filepath}/{pfname}", 'rb') as f:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            output = u.load()
        print("Found pickled PFSS file!")

    # If not, then download MHI mag, calc. PFSS, and save as picle file for next time
    except FileNotFoundError:
        print("PFSS file not found.\nDownloading...")
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


def circle_around(x, y, n, r=0.1):
    """
    Produces new points around a (x,y) point in a circle.
    At the moment does not work perfectly in the immediate vicinity of either pole.

    params
    -------
    x,y: int/float
            Coordinates of the original point
    n: int
            The amount of new points around the origin
    r: int/float (default = 0.1)
            The radius of the circle at which new points are placed

    returns
    -------
    pointlist: list[float]
            List of new points (tuples) around the original point in a circle, placed at equal intervals
    """

    origin = (x, y)

    x_coords = np.array([])
    y_coords = np.array([])
    for i in range(0, n):

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

        x_coords = np.append(x_coords, newx)
        y_coords = np.append(y_coords, newy)

    pointlist = np.array([x_coords, y_coords])

    return pointlist


def vary_flines(lon, lat, hmimap, n_varies, rss):
    """
    Finds a set of sub-pfss fieldlines connected to or very near a single footpoint on the pfss.
    
    lon: longitude of the footpoint [rad]
    lat: latitude of the footpoint [rad]
    
    n_varies:   tuple that holds the amount of circles and the number of dummy flines per circle
                if type(n_varies)=int, consider that as the amount of circles, and set the 
                amount of dummy flines per circle to 16

    params
    -------
    lon: int/float
            The longitude of the footpoint in radians
    lat: int/float
            The latitude of the footpoint in radians
    hmimap: hmi_synoptic_map object
            The pfss-solution used to calculate the field lines
    n_varies: list[int,int] or int 
            A list that holds the amount of circles and the number of dummy flines per circle
            if type(n_varies)=int, consider that as the amount of circles, and set the
            amount of dummy flines per circle to 16
    rss: float
            Heliocentric height of the source surface

    returns
    -------
    coordlist: list[float,float,float]
            List of coordinate triplets of the original field lines (lon,lat,height)
    flinelist: list[FieldLine-object]
            List of Fieldline objects of the original field lines
    varycoords: list[float,float,float]
            List of coordinate triplets of the varied field lines
    varyflines: list[FieldLine-object]
            List of Fieldline objects of the varied field lines
    """

    # Field lines per n_circles (circle)
    if isinstance(n_varies,list):
        n_circles = n_varies[0]
        n_flines = n_varies[1]
    else:
        n_circles = n_varies
        n_flines = 16

    # First produce new points around the given lonlat_pair
    lons, lats= np.array([lon]), np.array([lat])
    increments = np.array([0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29])
    for circle in range(n_circles):

        newlons, newlats = circle_around(lon,lat,n_flines,r=increments[circle])
        lons, lats = np.append(lons,newlons), np.append(lats,newlats)

    pointlist = np.array([lons, lats])

    # Trace fieldlines from all of these points
    varycoords, varyflines = get_field_line_coords(pointlist[0],pointlist[1],hmimap, rss)

    # Because the original fieldlines and the varied ones are all in the same arrays,
    # Extract the varied ones to their own arrays
    coordlist, flinelist = [], []

    # Total amount of flines = 1 + (circles) * (fieldlines_per_circle)
    total_per_fp = n_flines*n_circles+1
    erased_indices = []
    for i in range(len(varycoords)):
        # n_flines*n_circles = the amount of extra field lines between each "original" field line
        if i%(total_per_fp)==0:
            erased_indices.append(i)
            # pop(i) removes the ith element from the list and returns it
            # -> we append it to the list of original footpoint fieldlines
            coordlist.append(varycoords[i]) #.pop(i)
            flinelist.append(varyflines[i])

    # Really ugly quick fix to erase values from varycoords and varyflines
    for increment, index in enumerate(erased_indices):
        varycoords.pop(index-increment)
        varyflines.pop(index-increment)

    return coordlist, flinelist, varycoords, varyflines


def get_coord_values(field_line):
    """
    Gets the coordinate values from FieldLine object and makes sure that they are in the right order.

    params
    -------
    field_line: FieldLine object

    returns
    -------
    fl_r: list[float]
            The list of heliocentric distances of each segment of the field line
    fl_lon: list[float]
            The list of Carrington longitudes of each field line segment
    fl_lat: list[float]
            The list of Carrington latitudes of each field line segment
    """

    # first check that the field_line object is oriented correctly (start on photosphere and end at pfss)
    fl_coordinates = field_line.coords
    fl_coordinates = check_field_line_alignment(fl_coordinates)

    fl_r = fl_coordinates.radius.value / const.R_sun.value
    fl_lon = fl_coordinates.lon.value
    fl_lat = fl_coordinates.lat.value

    return fl_r, fl_lon, fl_lat


def get_field_line_coords(longitude, latitude, hmimap, seedheight):
    """
    Returns triplets of open magnetic field line coordinates, and the field line object itself

    params
    -------
    longitude: int/float
            Carrington longitude of the seeding point for the FieldLine tracing
    latitude: int/float
            Carrington latitude of the seeding point for the FieldLine tracing
    hmimap: hmi_synoptic_map object
    seedheight: float
            Heliocentric height of the seeding point

    returns
    -------
    coordlist: list[list[float,float,float]]
            The list of lists  of all coordinate triplets that correspond to the FieldLine objects traced
    flinelist: list[FieldLine]
            List of all FieldLine objects traced
    """

    # The amount of coordinate triplets we are going to trace
    try:
        coord_triplets = len(latitude)
    except TypeError:
        coord_triplets = 1
        latitude = [latitude]
        longitude = [longitude]

    # The loop in which we trace the field lines and collect them to the coordlist
    coordlist = []
    flinelist = []
    for i in range(coord_triplets):

        # Inits for finding another seed point if we hit null or closed line
        turn = 'lon'
        sign_switch = 1
        
        # Steps to the next corner, steps taken
        corner_tracker = [1,0]
        
        init_lon, init_lat = longitude[i], latitude[i]

        # Keep tracing the field line until a valid one is found
        while(True):

            # Trace a field line downward from the point lon,lat on the pfss
            fline = trace_field_line(longitude[i], latitude[i], hmimap, seedheight=seedheight)

            radius0 = fline.coords.radius[0].value
            radius9 = fline.coords.radius[-1].value
            bool_key = (radius0==radius9)

            # If fline is not a valid field line, then spiral out from init_point and try again
            # Also check if this is a null line (all coordinates identical)
            # Also check if polarity is 0, meaning that the field line is NOT open
            if( (len(fline.coords) < 10) or bool_key or fline.polarity==0): #fline.polarity==0 

                longitude[i], latitude[i], sign_switch, corner_tracker, turn = spiral_out(longitude[i], latitude[i], sign_switch, corner_tracker, turn)

            # If there was nothing wrong, break the loop and proceed with the traced field line
            else:
                break

            # Check that we are not too far from the original coordinate
            if(corner_tracker[0] >= 10):
                print(f"no open field line found in the vicinity of {init_lon},{init_lat}")
                break

        # Get the field line coordinate values in the correct order
        # Start on the photopshere, end at the pfss
        fl_r, fl_lon, fl_lat = get_coord_values(fline)

        # Fill in the lists
        triplet = [fl_r, fl_lon, fl_lat]
        coordlist.append(triplet)
        flinelist.append(fline)

    return coordlist, flinelist


def spiral_out(lon, lat, sign_switch, corner_tracker, turn):
    """
    Moves the seeding point in an outward spiral.
    
    Parameters
    ---------
    lon, lat: float
            the carrington coordinates on a surface of a sphere (sun or pfss)
    sign_switch: int
            -1 or 1, dictates the direction in which lon or lat is 
            incremented
    corner_tracker: tuple
            first entry is steps_unti_corner, int that tells how many steps to the next corner of a spiral
            the second entry is steps taken on a given spiral turn
    turn: str
            indicates which is to be incremented, lon or lat
            
    returns
    -----------
    lon, lat: float
            new coordinate pair
    """
    
    # In radians, 1 rad \approx 57.3 deg
    step = 0.01
    
    # Keeps track of how many steps until it's time to turn
    steps_until_corner, steps_moved = corner_tracker[0], corner_tracker[1]

    if turn=='lon':
        
        lon = lon + step*sign_switch
        lat = lat
        
        steps_moved += 1
        
        # We have arrived in a corner, time to move in lat direction
        if steps_until_corner == steps_moved:
            steps_moved = 0
            turn = 'lat'

        return lon, lat, sign_switch, [steps_until_corner, steps_moved], turn


    if turn=='lat':
        
        lon = lon
        lat = lat + step*sign_switch
        
        steps_moved += 1

        # Hit a corner; start moving in the lon direction
        if steps_until_corner == steps_moved:
            steps_moved = 0
            steps_until_corner += 1
            turn = 'lon'
            sign_switch = sign_switch*(-1)
        
        return lon, lat, sign_switch, [steps_until_corner, steps_moved], turn


def au_to_km(list_of_distances):
    """
    Takes a list of distances in AU and converts them to kilometers
    """
   
    # Convert to numpy array for performance
    list_of_distances = np.asarray(list_of_distances)

    # Add units of astronomical units
    list_of_distances = list_of_distances * u.au
    
    # Convert to units of kilometers
    list_of_distances = list_of_distances.to(u.km)

    return list_of_distances.value


def multicolorline(x, y, cvals, ax, vmin=-90, vmax=90):
    """
    Constructs a line object, with each segemnt of the line color coded
    Original example from: https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html

    params
    -------
    x, y: float
    cvals: str
    ax: Figure.Axes object
    vmin, vmax: int (default = -90, 90)

    returns
    -------
    line: LineCollection object
    """

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(vmin, vmax)

    cmrmap = cmr.redshift

    # sample the colormaps that you want to use. Use 90 from each so there is one
    # color for each degree
    colors_pos = cmrmap(np.linspace(0.0, 0.30, 45))
    colors_neg = cmrmap(np.linspace(0.70, 1.0, 45))

    # combine them and build a new colormap
    colors = np.vstack((colors_pos, colors_neg))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    # establish the linecollection object
    lc = LineCollection(segments, cmap=mymap, norm=norm)

    # set the values used for colormapping
    lc.set_array(cvals)

    # set the width of line
    lc.set_linewidth(3)

    # this we want to return
    line = ax.add_collection(lc)

    return line


def plot3d(field_lines, names, color_code='polarity'):
    """
    Creates an interactive 3D plot that the user is free to rotate and zoom in a Jupyter notebook.

    params
    -------
    field_lines : FieldLine object, or a list of them

    names : str
            names of the objects corresponding to the tracked field lines
    
    color_code : str
            either 'polarity' or 'object', dictates the color coding of the field lines
    """

    if not isinstance(field_lines, list):
        field_lines = [field_lines]

    if isinstance(field_lines[0], list):
        modulator = len(field_lines[1])//len(names) #modulates the order of field lines and the colors picked from c_list
        field_lines = flatten_flines(field_lines, modulator)

    if color_code=='object':
        num_objects = len(names)
        if len(field_lines) % num_objects != 0:
            raise Exception("Names do not match field lines")
        c_list = [get_color(x) for x in names]
    elif color_code=='polarity':
        colors = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}
    else:
        raise Exception("Choose either 'polarity' or 'object' as the color coding.")

    fig, axarr = plt.subplots(subplot_kw={"projection": "3d"})

    axarr.set_box_aspect((1, 1, 1))

    axarr.set_xlabel(r"x / $R_{\odot}$", labelpad=20)
    axarr.set_ylabel(r"y / $R_{\odot}$", labelpad=20)
    axarr.set_zlabel(r"z / $R_{\odot}$", labelpad=20)

    # Draw the Sun
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    axarr.plot_wireframe(x*1.0, y*1.0, z*1.0, color="darkorange")
    axarr.set_xlim(-2, 2)
    axarr.set_ylim(-2, 2)
    axarr.set_zlim(-2, 2)

    for i, field_line in enumerate(field_lines):
        coords = field_line.coords
        coords.representation = 'cartesian'

        if color_code=='polarity':
            color = colors.get(field_line.polarity)
        if color_code=='object':
            color = c_list[i//(modulator+1)]
        axarr.plot(coords.cartesian.x / const.R_sun,
                   coords.cartesian.y / const.R_sun,
                   coords.cartesian.z / const.R_sun,
                   color=color, linewidth=1)

    try:
        axarr.set_aspect('equal', adjustable='box')
    except NotImplementedError:
        axarr.set_aspect('auto', adjustable='box')


def draw_fieldlines(field_lines, rss=2.5, frame='yz', color_code='polarity', names=[], save=False):

    # check if there's a list inside a list, if there is -> flatten
    # NOTICE: flatten() messes up the order of the field lines, putting the original flines
    # first and then the varied field lines
    if isinstance(field_lines[0], list):
        modulator = len(field_lines[1])//len(names)  # modulates the order of field lines and the colors picked from c_list
        field_lines = flatten_flines(field_lines, modulator) # flatten_flines() conserves the correct ordering of field line objects

    fig, ax = plt.subplots(figsize=[10, 10])
    ax.set_aspect('equal')

    ax.set_xlabel(frame[0]+r" / $R_{\odot}$")
    ax.set_ylabel(frame[1]+r" / $R_{\odot}$")

    # set color coding to field lines
    if color_code=='object':
        num_objects = len(names)
        if len(field_lines) % num_objects != 0:
            raise Exception("Names do not match field lines")
        c_list = [get_color(x) for x in names]
    elif color_code=='polarity':
        colors = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}
    else:
        raise Exception("Choose either 'polarity' or 'object' as the color coding.")

    if(isinstance(field_lines, list)):

        for i, field_line in enumerate(field_lines):
            coords = field_line.coords
            coords.representation = 'cartesian'

            if color_code == 'polarity':
                color = colors.get(field_line.polarity)
            if color_code == 'object':
                color = c_list[i//(modulator+1)]

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

    ax.add_patch(mpatch.Circle((0, 0), 1, color='darkorange', lw=2.0, fill=False))
    ax.add_patch(mpatch.Circle((0, 0), rss, color='k', linestyle='--', fill=False))

    plt.title(projection)

    if save:
        plt.savefig('overhead.png')

    plt.show()


def flatten_flines(field_lines, modulator):
    '''
    Flattens a list of field_lines in such a way that the order of
    corresponding field lines and their varied counterparts is preserved.
    '''
    flines0, flines1 = field_lines[0], field_lines[1]
    flines_all = []
    # append the original field line, and then all the varied field lines corresponding to
    # that original field line in order
    for i in range(len(flines0)):
        flines_all.append(flines0[i])
        for j in range(modulator):
            index = i*modulator+j #index takes into account that there are a modulator amount of dummy flines
            flines_all.append(flines1[index])

    return flines_all


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


def trace_field_line(lon0, lat0, hmimap, seedheight, rad=True):
    '''
    Traces a single open magnetic field line at coordinates (lon0,lat0) on the pfss down
    to the photosphere

    Parameters:
    -----------
    lon0, lat0: float
            Longitude and latitude of the seedpoint
    hmimap: hmimap-object

    seedheight: float
            The height at which field line tracing is started (in solar radii)
    rad: bool, (default True)
            Wether or not input coordinates are in radians. If False, consider them degrees
    
    Returns:
    --------
    field_lines: FieldLine or list[FieldLine]
            A FieldLine object, or a list of them, if input coordinates were a list
    
    '''

    # if lat0 and lon0 are given in deg for some reason, transform them to rad
    if not rad:
        lat0 = np.deg2rad(lat0)
        lon0 = np.deg2rad(lon0)

    # Start tracing from a given height
    height = seedheight*const.R_sun
    tracer = tracing.PythonTracer()

    # Add unit to longitude and latitude, so that SkyCoord understands them
    lon, lat = lon0*units.rad, lat0*units.rad

    # Seed the starting coordinate at the desired coordinates
    seed = SkyCoord(lon, lat, height, frame=hmimap.coordinate_frame)

    # Trace the field line from the seed point given the hmi map
    field_lines = tracer.trace(seed, hmimap)

    # Field_lines could be list of len=1, because there's only one seed point given to the tracer
    if len(field_lines) == 1:
        return field_lines[0]
    else:
        return field_lines


def parker_spiral(sw, distance, longitude, resolution=1000, endpoint=2.5, backtrack=True):
    '''
    construct one magnetic parker spiral arm

    params
    ---------------------
    sw: solar wind speed in km/s
    distance: distance to the object in km
    longitude: angular coordinate of the object (stellar longitude?) in deg
    resolution: resolution of the curve (amount of points making the curve)
    endpoint: the point at which one wants to stop tracing back the spiral arm in solar radii

    returns
    ---------------------
    phi: array of angular coordinates in rad, type = [float]
    r: array of radial coordinates in solar radii, type = [float]
    '''

    # parker spiral solution e.g. here:
    # http://www.physics.usyd.edu.au/~cairns/teaching/2010/lecture8_2010.pdf

    omega = 2.694e-6 #rad/s

    r = np.linspace(endpoint, distance, resolution) #linspace to get even spacing

    # backtracking means going from sc to source surface
    if backtrack:
        phi = longitude + (omega)*(distance-r)/sw
    # if not backtracking, solve parker spiral from pfss outward all the way to max distance
    else:
        phi = longitude - (omega*r)/sw

    return phi, r


def symlog_pspiral(sw, distance, longitude, latitude, hmimap, names=None, title='', rss=2.5,
                   vary=False, n_varies=1, reference_longitude=None, save=False):
    '''
    Produces a figure of the heliosphere in polar coordinates with logarithmic r-axis outside the pfss.
    Also tracks an open field line down to photosphere given a point on the pfss.

    sw = solar wind in km/s, type: int, float or list
    distance = distance to the object in km, type: int, float or list
    longitude = angular coordinate of the object in deg (stellar longitude?), type: int, float or list
    latitude = see longitude
    rss = source surface height of the potential field, type: float
    reference_longitude : draws an arrow pointing radially outward, degree type: int or floar
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

    # Treat lons and lats as radians in the function
    longitude = np.deg2rad(longitude)
    latitude = np.deg2rad(latitude)

    # Normalize variables to solar radii
    if isinstance(sw, list):
        sw_norm = [(u*1000)/sun_radius for u in sw]
        distance_norm = [(d*1000)/sun_radius for d in distance]

    else:
        sw_norm = (sw*1000)/sun_radius
        distance_norm =(distance*1000)/sun_radius

    # Projection of the objects on the plane of ecliptic
    projection = np.cos(latitude)

    # Calculate parker spiral for given objects
    if isinstance(sw, list):
        phis, rs = [], []
        for i in range(len(sw)):
            phi, r = parker_spiral(sw_norm[i], distance_norm[i], longitude[i], resolution=1000, endpoint=rss)
            phis.append(phi)
            rs.append(r)

        sc_footphis = [phi[0] for phi in phis]
        sc_footrs = [r[0] for r in rs]
        sc_footpoint = [sc_footphis, sc_footrs]

    else:
        phi, r = parker_spiral(sw_norm, distance_norm, longitude, resolution=1000, endpoint=rss)
        sc_footpoint = [phi[0], r[0]]

    # calculate parker spiral for ref longitude, if one is inputted
    if reference_longitude is not None:
        default_sw = 400e3 / sun_radius
        max_distance = 480e9 / sun_radius
        ref_phi, ref_r = parker_spiral(default_sw, max_distance, np.deg2rad(reference_longitude), backtrack=False)

    # ----------------------------------------------------------
    # tracing the closest field line to sc_footpoint down to photosphere:

    # acquire an array of (r,lon,lat) coordinates of the open field lines under the pfss
    # based on the footpoint(s) of the sc
    if vary:

        # If there is more than one objects being both traced and varied, we have to vary them one at a time
        if len(sc_footpoint[0]) > 1:

            fline_triplets, fline_objects, varyfline_triplets, varyfline_objects = [], [], [], []

            for i, footpoint in enumerate(sc_footpoint[0]):
                
                #Append doesn't work here, but a simple + does. I wonder why.
                tmp_triplets, tmp_objects, varytmp_triplets, varytmp_objects = vary_flines(footpoint, latitude[i], hmimap, n_varies, rss)
                fline_triplets = fline_triplets + tmp_triplets
                fline_objects = fline_objects + tmp_objects
                varyfline_triplets = varyfline_triplets + varytmp_triplets
                varyfline_objects = varyfline_objects + varytmp_objects

        # If only a single object, then just run vary_flines for it
        else:
            fline_triplets, fline_objects, varyfline_triplets, varyfline_objects = vary_flines(sc_footpoint[0], latitude, hmimap, n_varies, rss)


    # If no varying, then just get one field line from get_field_line_coords()
    else:
        fline_triplets, fline_objects = get_field_line_coords(sc_footpoint[0], latitude, hmimap, rss)

    # We need fl_r, fl_lon, fl_lat,
    # they are located in fline_triplets[i][j]

    # Source surface:
    theta = 2*np.pi*np.linspace(0, 1, 200)
    ss = np.ones(200)*rss

    # Plotting commands ->
    fig, ax = plt.subplots(figsize = [23,17], subplot_kw={'projection': 'polar'})

    # Plot the source_surface and solar surface
    ax.plot(theta, ss, c='k', ls='--', zorder=1)
    ax.plot(theta, np.ones(200), c='darkorange', lw=2.5, zorder=1)
    
    # Plot the 30 and 60 deg lines on the Sun
    ax.plot(theta, np.ones(len(theta))*0.866, c='darkgray', lw=0.5, zorder=1) #cos(30deg) = 0.866(O)
    ax.plot(theta, np.ones(len(theta))*0.500, c='darkgray', lw=0.5, zorder=1) #cos(60deg) = 0.5(0)

    # Plot the spiral
    if isinstance(sw, list):
        for i in range(len(sw)):
            ax.plot(phis[i], projection[i]*rs[i], c=get_color(names[i]), label='sw={} km/s'.format(int(sw[i])), zorder=1)
            ax.scatter(sc_footpoint[0][i], projection[i]*sc_footpoint[1][i], s=115, c=get_color(names[i]), marker='D', zorder=3)
            ax.scatter(phis[i][-1], projection[i]*rs[i][-1], c=get_color(names[i]), alpha=1.0, s=175, zorder=3)

    else:
        ax.plot(phi, projection*r, c='k', label='sw={} km/s'.format(int(sw)), zorder=1)
        # Mark the footpoint of the spiral arm on the source surface
        ax.scatter(sc_footpoint[0], projection*sc_footpoint[1], s=140, c='C0', marker='x', zorder=3)
        # Mark the object
        ax.scatter(phi[-1], projection*r[-1], c='C0', s=175, zorder=3)

    # The reference longitude arrow and parker spiral
    if reference_longitude is not None:

        ref_lon_rad = np.deg2rad(reference_longitude)
        arrowprops = {'width':0.8, 'headwidth':9, 'color':'black'}

        ax.annotate(text=None, xy=(ref_lon_rad, 2.5), 
                    xytext=(ref_lon_rad, 1.0),
                    arrowprops=arrowprops)

        ax.plot(ref_phi, ref_r, c='black', ls='--', label='reference_longitude')

    # it seems that all other plotting must be done before this line
    # the exact reason for this still eludes me (2022-04-29, Vappu Eve) maybe something to do with multicolorline()?

    # plot the field line(s) connecting ss_footpoint and the solar surface and collect relevant
    # points (footpoints) to an array
    display_fl_footpoints = []
    display_fl_sourcepoints = []
    for fline in fline_triplets:

        fl_r   = fline[0]
        fl_lon = fline[1]
        fl_lat = fline[2]

        # Remember the foot and source points
        display_fl_sourcepoints.append([np.round(fl_r[0],1), np.round(fl_lon[0],1), np.round(fl_lat[0],1)])
        display_fl_footpoints.append([np.round(fl_r[-1],1), np.round(fl_lon[-1],1), np.round(fl_lat[-1],1)])

        # Plot the color coded field line
        fieldline = multicolorline(np.deg2rad(fl_lon), np.cos(np.deg2rad(fl_lat))*fl_r, ax=ax, cvals=fl_lat, vmin=-90, vmax=90)

    if vary:

        for fline in varyfline_triplets:

            fl_r   = fline[0]
            fl_lon = fline[1]
            fl_lat = fline[2]

            # Plot the color coded varied field line
            fieldline = multicolorline(np.deg2rad(fl_lon), np.cos(np.deg2rad(fl_lat))*fl_r, ax=ax, cvals=fl_lat, vmin=-90, vmax=90)


    r_max = 500e9 / sun_radius
    ax.set_rmax(r_max)
    ax.set_rscale('symlog', linthresh=rss)

    ax.grid(True)

    # These commands allow for plotting just a specific sector of the whole circle:
    # ax.set_thetamin(0)
    # ax.set_thetamax(135)
    ax.set_rticks([1.0, rss, 10.0, 100.0])
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    rlabels = ['1', str(rss), r'$10^1$', r'$10^2$']
    ax.set_yticklabels(rlabels)

    # Create the colorbar displaying values of the last fieldline plotted
    cb = fig.colorbar(fieldline)
    
    # @TODO: shrink the colorbar and move it to the top right corner
    
    # Colorbar is the last object created -> it is the final index in the list of axes
    cb_ax = fig.axes[-1]
    cb_ax.set_ylabel('Heliographic latitude [deg]')
    
    # Before adding txtboxes, make sure that sc_footpoint is of proper shape
    if(isinstance(sc_footpoint[0],float)):
        display_sc_footpoint = [[sc_footpoint[0]],[sc_footpoint[1]]]
        latitude = [latitude]

    # Also add magnetic polarity info to display:
    display_polarities = []
    for fline in fline_objects:
        display_polarities.append(int(fline.polarity))

    # Txtbox stuff----------->
    txtsize = 12
    abbrv_height = 880
    
    # Make legend for abbreviations:
    legendlabel = "ABBREVIATIONS \nPS: Photosphere \nSS: Source Surface \nFP: Footpoint \nP: Polarity"
    legendbox = AnchoredText(f"{legendlabel}", prop=dict(size=txtsize+1), frameon=True, loc=(10), bbox_to_anchor=(92,abbrv_height)) #bbox_to_anchor=(92,880)
    legendbox.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
    legendbox.patch.set_linewidth(3.0)
    legendbox.patch.set_edgecolor(get_color(None))
    ax.add_artist(legendbox)


    # Add textbox for footpoint coordinates
    if(names is not None):
        for i in range(len(display_fl_footpoints)):
            plabel = AnchoredText(f"{names[i]}, sw={int(sw[i])} km/s\nPS FP = ({display_fl_sourcepoints[i][1]},{display_fl_sourcepoints[i][2]}) \nSS FP = ({display_fl_footpoints[i][1]},{display_fl_footpoints[i][2]}) \nP: {display_polarities[i]}", 
                                  prop=dict(size=txtsize+1), frameon=True, loc=(2), bbox_to_anchor=(15, (abbrv_height-40-i*75) )) #bbox_to_anchor=(15, (840-i*75) )
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
    
    # Set the title of the figure
    plt.title(title)

    if(save):
        plt.savefig('pfss_constellation_figure.png', transparent=False, facecolor='white', bbox_inches='tight')

    # if using streamlit, send plot to streamlit output, else call plt.show()
    if _isstreamlit():
        import streamlit as st
        st.pyplot(fig)  # , dpi=200)
    else:
        plt.show()

    # The function returns all the calculated field line objects, which include many attributes
    # of the field lines such as coordinates, polarity, and wether they are open or closed
    if vary:
        return [fline_objects, varyfline_objects]
    else:
        return fline_objects


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

    # Init arrays to store into dataframe
    polarities = np.array([])
    photosphere_lons = np.array([])
    photosphere_lats = np.array([])
    pfss_lons = np.array([])
    pfss_lats = np.array([])

    for line in flines:

        polarities = np.append(polarities, line.polarity)

        coordinates = check_field_line_alignment(line.coords)

        photospheric_footpoint = (coordinates.lon.value[0], coordinates.lat.value[0])
        photosphere_lons = np.append(photosphere_lons, photospheric_footpoint[0])
        photosphere_lats = np.append(photosphere_lats, photospheric_footpoint[1])

        pfss_footpoint = (coordinates.lon.value[-1], coordinates.lat.value[-1])
        pfss_lons = np.append(pfss_lons, pfss_footpoint[0])
        pfss_lats = np.append(pfss_lats, pfss_footpoint[1])

    data_dict = {'Observing object': names,
                 'Solar surface \nfootpoint lon': photosphere_lons,
                 'Solar surface \nfootpoint lat': photosphere_lats,
                 'PFSS \nfootpoint lon': pfss_lons,
                 'PFSS \nfootpoint lat': pfss_lats,
                 'Magnetic Polarity': polarities}

    df = pd.DataFrame(data=data_dict)

    return df


def df_to_file(df, filename: str):
    '''
    Writes a dataframe in to a csv file
    '''

    if not isinstance(filename, str):
        raise Exception("The file name must be string.")

    current_directory = os.getcwd()

    filestr = f"{current_directory}/{filename}.csv"

    df.to_csv(filestr)
    print(f"Created file {filename}.csv to {current_directory}/")


def write_info_to_csv(flines, names=[None], filename='magnetic_info'):

    # first assemble a dataframe from the field line object
    coords_df = field_line_info_to_df(flines, names)

    # write the df into csv
    df_to_file(coords_df, filename)


# ============================================================================================================================

# Coordinate transformations:

def heeq2hgc(xyz_list, obstimes, observer='earth', unit=None, returns='objects'):
    """
    Takes a list of cartesian xyz coordinates in HEEQ (Heliospheric Earth Equatorial) frame and transforms them into HGC (HelioGraphic Carrington) coordinates.
    """

    from astropy.coordinates import SkyCoord, CartesianRepresentation
    # from sunpy.coordinates import HeliographicCarrington
    import astropy.units as units
    from sunpy.coordinates import frames

    if unit is None:
        unit = units.AU

    if returns != 'objects' and returns != 'coordinates':
        raise Exception("Choose either 'objects' or 'coordinates' for the function return.")

    coordlist = []
    for i, xyz in enumerate(xyz_list):

        # First seek a spherical stonyhurst representation for the HEEQ cartesian representation
        c_stonyhurst = SkyCoord(CartesianRepresentation(xyz[0]*unit, xyz[1]*unit, xyz[2]*unit),
                                obstime=obstimes[i],
                                frame="heliographic_stonyhurst")

        # Convert stonyhurst coordinates (effectively just longitude) to carrington coordinates
        c_carrington = c_stonyhurst.transform_to(frames.HeliographicCarrington(observer='earth'))

        coordlist.append(c_carrington)

    if returns=='objects':
        return coordlist
    if returns=='coordinates':
        return [(x.lon.value, x.lat.value, x.radius.value) for x in coordlist]


def arcsec_to_carrington(arc_x, arc_y, time):

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


def _isnotebook():
    # https://stackoverflow.com/a/39662359/2336056
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False


# ============================================================================================================================
# this will be ran when importing
# set Jupyter notebook cells to 100% screen size:
# if _isnotebook():
#     from IPython.core.display import display, HTML
#     display(HTML(data="""<style> div#notebook-container { width: 99%; } div#menubar-container { width: 85%; } div#maintoolbar-container { width: 99%; } </style>"""))
