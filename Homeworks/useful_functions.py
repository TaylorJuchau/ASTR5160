import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
from astropy.table import Table, vstack, Column, Row
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree
import subprocess
import requests


def decode_sweep_name(sweepname, nside=None, inclusive=True, fact=4):
    """Retrieve RA/Dec edges from a full directory path to a sweep file

    Parameters
    ----------
    sweepname : :class:`str`
        Full path to a sweep file, e.g., /a/b/c/sweep-350m005-360p005.fits
    nside : :class:`int`, optional, defaults to None
        (NESTED) HEALPixel nside
    inclusive : :class:`book`, optional, defaults to ``True``
        see documentation for `healpy.query_polygon()`
    fact : :class:`int`, optional defaults to 4
        see documentation for `healpy.query_polygon()`

    Returns
    -------
    :class:`list` (if nside is None)
        A 4-entry list of the edges of the region covered by the sweeps file
        in the form [RAmin, RAmax, DECmin, DECmax]
        For the above example this would be [350., 360., -5., 5.]
    :class:`list` (if nside is not None)
        A list of HEALPixels that touch the  files at the passed `nside`
        For the above example this would be [16, 17, 18, 19]
    """
    # ADM extract just the file part of the name.
    sweepname = os.path.basename(sweepname)

    # ADM the RA/Dec edges.
    ramin, ramax = float(sweepname[6:9]), float(sweepname[14:17])
    decmin, decmax = float(sweepname[10:13]), float(sweepname[18:21])

    # ADM flip the signs on the DECs, if needed.
    if sweepname[9] == 'm':
        decmin *= -1
    if sweepname[17] == 'm':
        decmax *= -1

    if nside is None:
        return [ramin, ramax, decmin, decmax]

    pixnum = hp_in_box(nside, [ramin, ramax, decmin, decmax],
                       inclusive=inclusive, fact=fact)

    return pixnum


def is_it_in(objs, box):
    """Determine which of an array of objects are inside an RA, Dec box. Trimmed version of ADAM's function with inclusive boundaries

    Parameters
    -------------
    objs :  type = astropy.table or list - array of objects, must include columns titled "RA" and "DEC", or be SkyCoord objects
    box :  type = list - list of floats corresponding to the boundaries of the box
                *must be in format [ra_min, ra_max, dec_min, dec_max] with boundaries in degrees*
    
    Returns
    -------------
    original objs replaced with boolean values stating whether they are in the box (True) or not (False)
        
    """

    ramin, ramax, decmin, decmax = box
    if type(objs)==SkyCoord:
        ii = ((objs.ra.deg >= ramin) & (objs.ra.deg <= ramax) & (objs.dec.deg >= decmin) & (objs.dec.deg <= decmax))
    else:
        ii = ((objs["RA"] >= ramin) & (objs["RA"] <= ramax) & (objs["DEC"] >= decmin) & (objs["DEC"] <= decmax))

    return ii #TJ return true or false array


def find_sweep_file(directory, obj):
    '''searches a directory with multiple .fits files and returns which file contains this one object

    imports all needed packages automatically
    -------------


    Parameters
    -------------
    directory :  type = str - string representing a directory to search, this directory should contain at least 1 .fits file
    obj : type = skycoord object
                        
    
    Returns
    -------------
    File name in directory '''

    fits_files = glob.glob(os.path.join(directory, "*.fits")) #TJ find only fits files in directory
    boxes = [decode_sweep_name(file) for file in fits_files]
    all_files = np.vstack([is_it_in(obj, box) for box in boxes])
    needed_index = np.where(all_files.any(axis=1))[0][0]
    needed_file = fits_files[needed_index]
    
    return needed_file



def find_sweep_files_with(directory, objects_list):
    '''searches a directory with multiple .fits files and returns which file(s) contain at least one of your supplied objects

    imports all needed packages automatically
    -------------


    Parameters
    -------------
    directory :  type = str - string representing a directory to search, this directory should contain at least 1 .fits file
    ra :  type = float - right ascension of desired object in degrees
    dec : type = float - declination of desired object in degrees
                        
    
    Returns
    -------------
    list of file names in directory, each of these files have objects that have BOTH larger and smaller ra's AND dec's'''

    fits_files = sorted(glob.glob(os.path.join(directory, "*.fits"))) #TJ find only fits files in directory
    boxes = [decode_sweep_name(file) for file in fits_files] #TJ extract boxes that each fits file contains
    all_files = np.atleast_2d([is_it_in(objects_list, box) for box in boxes])
    #TJ for whether or not that object is in that file's box. An array of all Falses would mean no objects are in that file's box
    needed_indices = np.where(all_files.any(axis=1))[0] #TJ this grabs the indices corresponding to all the sublists where at least one True exists
    needed_files = [fits_files[i] for i in needed_indices] #TJ this uses those indices to grab just the files where at least one object is inside
    
    return needed_files


def read_partial_fits(filepath, kept_column_names):
    '''read in table and only save necessary columns this saves memory and time with further manipulations compared to saving EVERY column
    -------------

    Parameters
    -------------
    filepath :  type = str - string representing path to .fits files you want to read in
    kept_column_names :  type = list - list of strings corresponding to which columns need to be read in. All other columns are ignored
    
    Returns
    -------------
    table with only the columns that we wanted to keep
    '''
    
    with fits.open(filepath) as hdul:
        data = hdul[1].data  #TJ read in all data
        filtered_data = {col: data[col] for col in kept_column_names} #TJ filter data to only keep necessary columns
        return Table(filtered_data)


def add_mag_column(original_table, fluxes):
    '''Add one new column in table for every column with name in fluxes (example: FLUX_x). New column will be named "x_mag".
    *Assumes fluxes are given in nanomaggy
    -------------

    Parameters
    -------------
    table :  type = astropy.table - table with at least one column name in "fluxes"
    fluxes : type = list - list of strings corresponding to which columns of fluxes you want to add a corresponding magnitude column

    Returns
    -------------
    Original table, with new columns for magnitudes
    '''

    for flux_col in fluxes:
        band = flux_col.split('_')[1] #TJ extract column name, "FLUX_R" will have the "R" pulled out and used to name new column "R_mag"
        new_col_name = f"{band}_mag"  #TJ define new column e.g., 'G_mag' for original "FLUX_G"

        original_table[new_col_name] = 22.5 - 2.5*np.log10(original_table[flux_col]) #TJ convert flux to magnitude
    return original_table




def get_dots_in_box(box, num_dots=1000000, radians=False):
    '''populates a sphere evenly with dots and selects only those within the box provided
    
    imports neeeded general packages automatically
    requires following functions from parent module "ASTRO5160/Homeworks/HW_2.py":
                        "find_latlong_area()"
                        "evenly_populate_sphere()"
                        "is_it_in()"
    
    Parameters
    -------------
    box : type = list - given in the format [ra_min, ra_max, dec_min, dec_max] with location in degrees
                                                           
    num_dots (optional, defaults to 1,000,000) : type = int - number of dots on entire sphere.     
                                               
    radians (optional, defaults to False): type = boolean - change to True if locations are in radians
    
    Returns
    -------------
    list of [ra,dec] for every dot that is inside the box
    '''
    
    if radians == False: #TJ if box isnt already in radians, convert to radians.
        box = np.radians(box)

    #TJ generate random points on a sphere
    ra, dec = evenly_populate_sphere(num_dots)
    all_dots = np.column_stack((ra, dec))  #TJ store as a NumPy array for easier indexing
    
    #TJ generate list of dots that are in box
    in_box = all_dots[is_it_in(all_dots, box)]
    
    
    return in_box



def aitoff_plot_region(path, show_plot=True, save_plot=None, radians = False, name = 'bounded_region'):
    '''Plot regions bounded by path with vertices given in ra, dec in degrees (or radians if specified)
        example: box_path = [[0,0], [0,15], [20,15], [10,0]] will plot an assymmetrical trapezoidal box
    
    imports needed packages automatically
    
    Parameters
    -------------
    path : type = list - contains entries for each corner of a closed path
                        all coordinates given in degrees (unless radians = True is specified)
                    
    show_plot (optional, defaults to True) : type = boolean - False will suppress plot output
                string multiple iterations of this function with show_plot = False, then call plt.show() to show 
                all regions in single plot
                    
    save_plot (optional, defaults to None) : type = str - path to directory where plots will be saved to
                             if directory is set to None, output plot will not be saved (can still print to screen)
    
    name (optional, defaults to 'bounded_region'): type = str - name of .png file to save as (if save_plot == True)
                                                                This is also the name of the legend label for the region

    radians (optional, defaults to False): type = boolean - specify whether or not corners are defined in radians or degrees
    
                For example:
                    path1 = [[0,0], [0,15], [20,15], [10,0]]
                    path2 = [[30,30], [30,45], [45,30]]
                    aitoff_plot_region(path1, show_plot=False, save_plot=None, radians = False, name = 'region1')
                    aitoff_plot_region(path2, show_plot=False, save_plot=None, radians = False, name = 'region2')
                    aitoff_plot_region([], show_plot=True, save_plot='C:/Users/tj360/ASTR5160', radians = False, name = 'testing_function')
                    

    
    Returns
    -------------
    Nothing, just used to plot regions
    
    '''

    #TJ the way this function is written makes it difficult to save a plot of multiple regions with an informative title
    #TJ to save a plot like this, call the function with an empty path, and save_plot = True, with name = figure_title
    if path == []: #TJ check that path is not empty, if it is, this is likely just to save a set of plots
        if not save_plot: #TJ this should never happem, print error message
            print('function called with a blank path but save_plot set to False')
            return None
        else: 
            os.makedirs(save_plot, exist_ok=True)  #TJ create the directory if it doesn't exist
            save_path = os.path.join(save_plot, f'{name}.png') #TJ define the full path for saving plot
            plt.savefig(save_path, dpi=300, bbox_inches='tight') #TJ save plot to specified directory
            print(f'output plot saved to {save_plot} as filename = {name}.png')
    else: #TJ if path is not empty, continue with plotting function
        if not radians: #TJ convert to radians for aitoff projection if needed
            path = np.radians(path)
            ra_path = [p[0] for p in path] #TJ define right ascension of corners
            dec_path = [p[1] for p in path] #TJ define declination of corners
            ra_path.append(path[0][0]) #TJ append with first value so the path is closed
            dec_path.append(path[0][1])
        elif radians:
            ra_path = [p[0] for p in path]
            dec_path = [p[1] for p in path]
            ra_path.append(path[0][0])
            dec_path.append(path[0][1])

        if not plt.get_fignums(): #If number of stored figures is zero, create new figure with appropriate axes
            plt.figure(figsize=(8, 5))
            ax = plt.subplot(111, projection="aitoff")
            ax.set_xlabel("Right Ascension", fontsize=12)
            ax.set_ylabel("Declination", fontsize=12)
            ax.set_title("Regions on Sphere", fontsize=12)
            ax.grid(color='black', linestyle='--', linewidth=0.5)
        else:
            ax = plt.gca() #TJ if a plot already exists, then just use that current plot and only add things to it

        #TJ plot the path through all ra,dec corners
        ax.plot(ra_path, dec_path, linestyle='-', linewidth=2, label = name) 
        ax.legend()


        if save_plot != None:
            os.makedirs(save_plot, exist_ok=True)  #TJ create the directory if it doesn't exist
            save_path = os.path.join(save_plot, f'{name}.png') #TJ define the full path for saving plot
            plt.savefig(save_path, dpi=300, bbox_inches='tight') #TJ save plot to specified directory
            print(f'output plot saved to {save_plot} with filename {name}.png') #TJ print if file was saved
        if show_plot:
            plt.show()
    return None


def aitoff_plot_rectangles(boxes, directory, name = 'lat-long box', radians=False):
    '''Plot lat-long rectangles bounded by ra_min, ra_max, dec_min, dec_max, where:
        boxes = [[ra_min1, ra_max1, dec_min1, dec_max1], [ra_min2, ra_max2, dec_min2, dec_max2]]
    
    imports needed packages automatically
    
    Parameters
    -------------
    boxes : type = list - contains entries for ra and dec boundaries of the lat-long rectangle
                        all coordinates given in degrees (unless radians = True is specified)
                        
    directory : type = str - path to directory where plots will be saved to
                             if directory is set to None, output plot will not be saved (will still print)
    
    name (optional, defaults to 'lat-long box'): type = str - name of .png file

    radians (optional, defaults to False): type = boolean - Whether or not boundaries are in radians not degrees
        
    Returns
    -------------
    Nothing, prints plot to screen and saves a .png file into working directory if that is specified
    '''

    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="aitoff") #TJ set projection format
    if isinstance(boxes[0], int) or isinstance(boxes[0], float):
        #TJ if there is only one box, perform function on the single box
        if radians == False:
            ra_min, ra_max, dec_min, dec_max = np.radians(boxes) #TJ this makes sure ra/dec are in radians
        elif radians == True:
            ra_min, ra_max, dec_min, dec_max = boxes
        #TJ make a path that circles around the whole box perimeter
        ra_path = [ra_min, ra_max, ra_max, ra_min, ra_min]
        dec_path = [dec_min, dec_min, dec_max, dec_max, dec_min]
        
        ax.plot(ra_path, dec_path, color='red',linestyle = '-', linewidth=2, 
                label = 'lat-long box')
    else: #TJ if the first entry in boxes is a list or array, iterate through each entry as individual boxes
        colors = cm.rainbow(np.linspace(1, 0, len(boxes))) #TJ this is prep for each box to be a different color
        areas = [] #TJ initialize area list
        for i, b in enumerate(boxes): #TJ iterate through each box and do what is described above to each box
            if radians == False:
                ra_min, ra_max = np.radians([b[0], b[1]])
                dec_min, dec_max = np.radians([b[2], b[3]])
            else:
                ra_min, ra_max = [b[0], b[1]]
                dec_min, dec_max = [b[2], b[3]]
            ra_path = [ra_min, ra_max, ra_max, ra_min, ra_min]
            dec_path = [dec_min, dec_min, dec_max, dec_max, dec_min]
            area = ((ra_max - ra_min)*(np.sin(dec_max) - np.sin(dec_min)))*((180/np.pi)**2)
            areas.append(area) #TJ this isnt optimal, but the "boxes" array will never be large, so it should be fine
            ax.plot(ra_path, dec_path, color=colors[i],linestyle='-', linewidth=2, 
                    label = f"box {i+1}")
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'] #TJ set right ascension tick labels
    ax.set_xticklabels(xlab, weight=300)
    ax.grid(color='black', linestyle='--', linewidth=0.5) #TJ make grid lines pretty thin so box borders show up better
    ax.set_xlabel("Right Ascension", fontsize=12)
    ax.set_ylabel("Declination", fontsize=12)
    ax.set_title(f"{name}", fontsize=12)
    
    ax.legend() #TJ I couldnt get the legend to appear at a consistent place outside of the plot

    

    if directory != None:
        os.makedirs(directory, exist_ok=True)  #TJ create the directory if it doesn't exist
        save_path = os.path.join(directory, name) #TJ define the full path for saving plot
        plt.savefig(save_path, dpi=300, bbox_inches='tight') #TJ save plot to specified directory
        print(f'output plot saved to {directory}')
    plt.show()
    return None


def cross_match_objects(match_all_these, find_matches_from_here, max_sep=1/3600):
    '''find objects in the find_matches_from_here table that are the closest RA,DEC from the objects in the match_all_these table

    -------------

    Parameters
    -------------
    match_all_these :  type = astropy.table - table with columns for "RA" and "DEC"
    find_matches_from_here :  type = astropy.table - table with columns for "RA" and "DEC"
    max_sep (optional, defaults to 1/3600) : type = float - fail to match objects if no object is closer than this many degrees away from it                  
    
    Returns
    -------------
    Table with only rows selected from the match_all_these table that found matches in the find_matches_from_here table
    Table from the find_matches_from_here table that are the valid matches to entries in the match_all_these table
    Table with only rows from the match_all_these table that were NOT found in the find_matches_here table
    '''
    #TJ create skycoord object for all objects in survey table, convert to cartesian for accurate ckd searching
    #TJ ckdtree uses euclidiean distance formula, so leaving it in ra and dec would distort near poles
    if type(find_matches_from_here) == SkyCoord:
        survey_coords = find_matches_from_here
    else:
        survey_coords = SkyCoord(find_matches_from_here['RA'], find_matches_from_here['DEC'], unit=u.deg)
    survey_xyz = np.vstack(survey_coords.cartesian.xyz).T 
    
    #TJ create skycoord object for all objects in observation table, convert to cartesian to remove dimensional distortion for cdk tree structuring
    if type(find_matches_from_here) == SkyCoord:
        obs_coords = match_all_these
    else:
        obs_coords = SkyCoord(match_all_these['RA'], match_all_these['DEC'], unit=u.deg)
    obs_xyz = np.vstack(obs_coords.cartesian.xyz).T
    tree = cKDTree(survey_xyz) #TJ build ckd tree for efficient searching

    #TJ query the tree for the nearest neighbor of each object in match_all_these table
    _, matched_indices = tree.query(obs_xyz) #TJ first output is the matched entry's euclidean distance from source, not easy to convert to deg
    
    separations = obs_coords.separation(survey_coords[matched_indices]).deg  #TJ calculate separation from each "closest match" object

    #TJ apply filters requested in HW3 prompts
    valid_matches = separations < max_sep

    #TJ extract the objects that satisfy all the above constraints, return BOTH tables' filtered rows
    filtered_match_all_these = match_all_these[valid_matches]
    filtered_matches_table = find_matches_from_here[matched_indices[valid_matches]]
    unmatched_match_all_these = match_all_these[~valid_matches]
    
    return filtered_match_all_these, filtered_matches_table, unmatched_match_all_these


def find_matching_object(object_to_match, find_match_from_here, max_sep=1/3600):
    '''find closest object in the find_matches_from_here table that is the closest RA,DEC to the object_to_match

    -------------

    Parameters
    -------------
    object_to_match :  type = SkyCoord object
    find_matches_from_here :  type = astropy.table - table with columns for "RA" and "DEC"
    max_sep (optional, defaults to 1/3600) : type = float - fail to match objects if no object is closer than this many degrees away from it                  
    
    Returns
    -------------
    Table with a single row that is the object that was closest. Includes all columns from find_matches_here table
    '''

    #TJ create SkyCoord object for the survey table and convert to cartesian
    survey_coords = SkyCoord(find_match_from_here['RA'], find_match_from_here['DEC'], unit=u.deg)
    survey_xyz = np.vstack(survey_coords.cartesian.xyz).T
    
    #TJ convert original object to cartesian
    obj_xyz = np.vstack(object_to_match.cartesian.xyz).T
    tree = cKDTree(survey_xyz)

    #TJ query for closest match
    _, matched_index = tree.query(obj_xyz[0])
    matched_obj = find_match_from_here[matched_index]

    #TJ compute separation and apply quality filters
    separation = object_to_match.separation(survey_coords[matched_index]).deg

    if separation < max_sep:
        return matched_obj  #TJ return matched object
    else:
        print(f'no matches were found within max_sep of {max_sep}')
        return None


def get_sdss_mags(ra, dec, radius = 1.2/3600):
    '''query SDSS navigator for objects within given radius of ra, dec location

    Note that the website will decline the 61st query made in under 1 minute
    This automatic throttling means we dont need to sleep for 1 second between queries if the total numebr is under 60 ;)

    -------------

    Parameters
    -------------
    ra :  type = float - right ascension in degrees of query object
    dec :  type = float - declination in degrees of query object
    radius (optional, defaults to 1.2 arcseconds) : 
    
    Returns
    -------------
    list of magnitudes in order u, g, r, i, z
    '''
    #TJ define SDSS query url

    url = f"http://skyserver.sdss.org/dr16/SkyServerWS/SearchTools/SqlSearch?cmd=SELECT+TOP+1+u,g,r,i,z,ra,dec+FROM+PhotoPrimary+WHERE+ra BETWEEN {ra-radius} AND {ra+radius}+AND+dec BETWEEN {dec-radius} AND {dec+radius}&format=json"
    response = requests.get(url).json() #TJ perform request and parse result
    if response[0]['Rows'] != []: #TJ if an object is returned, function returns the u and i magnitudes
        return response[0]['Rows'][0]['u'], response[0]['Rows'][0]['g'], response[0]['Rows'][0]['r'], response[0]['Rows'][0]['i'], response[0]['Rows'][0]['z']
    else:
        return np.nan, np.nan, np.nan, np.nan, np.nan

def what_is_it(obj, m = 0.8862876254180598, b= -0.652173913043478):
    '''Use a dividing line with slope = m and y-intercept = b so determine if the object is a star or a QSO
    defaults here were from a dataset in Week10.Class18_notebook. Currently only detects either QSO or star

    -------------

    (star["g_mag"]-star["z_mag"], star["r_mag"]-star["W1_mag"])

    Parameters
    -------------
    object :  type = astropy.table.row.Row - row with magnitude values for "g_mag", "z_mag", "r_mag", "W1_mag"
    m (optional, defaults to 0.886) :  type = float - slope of the dividing line
    b (optional, defaults to -0.6521) : type = float - y-intercept of dividing line                   
    
    Returns
    -------------
    string representing object classification
    'QSO' if object is above the dividing line
    'star' if object is below the dividing line
    '''
    gmz = obj["g_mag"] - obj["z_mag"]
    rmW1 = obj["r_mag"]-obj["W1_mag"]
    if rmW1 > (m * gmz + b):
        return 'QSO'
    elif rmW1 < (m * gmz + b):
        return 'star'
    else:
        print('something went wrong classifying this object')
        return None


if __name__ == "__main__":
    print('running test environment...')





