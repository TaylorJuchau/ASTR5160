#!/usr/bin/env python
# coding: utf-8

# In[23]:


def find_latlongrect_area(box, radians = False, vertices = False):
    '''find area enclosed in a lat-long box defined as box = [ra_min, ra_max, dec_min, dec_max] in square degrees
    
    *imports needed packages automatically
    
    Parameters
    -------------
    box : type = list - list of ra_min, ra_max, dec_min, dec_max
    
    radians (optional, defaults to False) : type = boolean - true if boundaries are defined in radians
    
    vertices (optional, defaults to False) : type = boolean - specifies whether the box is given as a set of boundaries
                                            or as a list of vertices 

    Returns
    -------------
    area contained within the lat-long box in square degrees
    
    '''

    import numpy as np
    
    if not radians:
        box = np.radians(box)
    if vertices: #TJ if box is given as a set of vertices, convert to format [ra_min, ra_max, dec_min, dec_max]
        ra_min = min([b[0] for b in box])
        ra_max = max([b[0] for b in box])
        dec_min = min([b[1] for b in box])
        dec_max = max([b[1] for b in box])
    if not vertices: #TJ if box is given as boundaries, assign ra/dec _ min/max accordingly
        ra_min, ra_max, dec_min, dec_max = box
    #TJ compute area and convert to square degrees
    return ((ra_max - ra_min)*(np.sin(dec_max) - np.sin(dec_min)))*((180/np.pi)**2) 


import numpy as np


def boundaries_to_vertices(box, radians = False):
    '''Convert a box defined as bounded by [ra_min, ra_max, dec_min, dec_max] to a box defined by its vertices
    new box will be given as set of vertices sorted to be counter-clockwise:
                        [[ra_min, dec_min], [ra_max, dec_min], [ra_max, dec_max], [ra_min, dec_max]]
    
    *imports needed packages automatically
    
    Parameters
    -------------
    dot : type = list - boundaries of lat-long box in the format [ra_min, ra_max, dec_min, dec_max]
    
    radians (optional, defaults to False) : type = boolean - True if coordinates of vertices are in radians

    Returns
    -------------
    vertices : type = list - vertices of lat-long box in format: 
    [[ra_min, dec_min], [ra_max, dec_min], [ra_max, dec_max], [ra_min, dec_max]] in degrees
    '''
    
    if radians:
        import numpy as np
        box = np.degrees(box)
    return [[box[0],box[2]], [box[1], box[2]], [box[1],box[3]], [box[0], box[3]]]
    
    

    
def cartesian(dot, radians = True):
    '''Convert (ra, dec) to Cartesian coordinates x, y, z.
    
    *imports needed packages automatically
    
    Parameters
    -------------
    dot : type = list - coordinates for a dot in right ascension and declination in format [ra,dec]
    
    radians (optional, defaults to True) : type = boolean - True if coordinates of vertices are in radians

    Returns
    -------------
    3-vector [x, y, z] representing the same point on a unit sphere as the given dot
    
    '''
    
    if not radians: #TJ convert to radians if needed
        dot = np.radians(dot)
    ra, dec = dot
    return [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
    


def find_polygon_area(vertices, radians = True):
    '''find area enclosed by a polygon in square degrees, made from vertices connected by sections of great circles.
    
    **Important note** lat-lon boxes are NOT bounded by sections of great circles!!!
    **Use find_latlongrect_area() instead for these regions.
    
    *imports needed packages automatically
    *needs
    Parameters
    -------------
    vertices : type = list - list of vertices [[ra1,dec1], [ra2,dec2], [ra3,dec3], etc]
    
    radians (optional, defaults to True) : type = boolean - True if coordinates of vertices are in radians
    
    Returns
    -------------
    area contained within polygon in square degrees
    
    '''
    
    import numpy as np
    
    if not radians: #TJ convert to radians if needed
        vertices = np.radians(vertices)

    num_corners = len(vertices) #TJ calculate the number of corners
    total_excess = 0 #TJ initialize total excess angle

    for i in range(num_corners):
        A = np.array(cartesian(vertices[i - 1]))
        B = np.array(cartesian(vertices[i]))
        C = np.array(cartesian(vertices[(i + 1) % num_corners]))

        #TJ compute angle at vertex B using dot product of great-circle vectors
        #TJ see Girard’s Theorem (Spherical Excess) via Cartesian Vectors for derivation
        #TJ np.linalg.norm calculates the magnitude of a vector
        numerator = np.dot(np.cross(A, B), np.cross(C, B)) 
        denominator = np.linalg.norm(np.cross(A, B)) * np.linalg.norm(np.cross(C, B))
        angle = np.arccos(numerator / denominator)
        
        total_excess += angle

    #TJ spherical excess formula
    E = total_excess - (num_corners - 2) * np.pi

    #TJ area on a sphere, converted to square degrees
    return E * (180/np.pi)**2


def evenly_populate_sphere(n_dots):
    '''generates n_dots number of ra and dec values (in radians) such that they will be evenly distributed over a sphere.

    *imports needed packages automatically
    
    Parameters
    -------------
    n_dots : type = int - integer number of dots you want to appear on the sphere
    

    Returns
    -------------
    ra : list of right ascension values in radians
    dec : list of declination values in radians
    '''
    
    import numpy as np
    from numpy.random import random
    ra = 2*np.pi*(random(n_dots)) #TJ generates random RA's measured in radians between 0 and 2pi 
    dec = np.arcsin(1.-random(n_dots)*2.) #TJ creates random declinations in radians between +pi/2 and -pi/2
    return ra, dec



def is_it_in(dots, lat_long_box, radians = False):
    '''determines whether a coordinate (or list of coordinates) is within a given lat-long box
    
    Parameters
    -------------
    dots : type = list - list of ra, dec coordinates in radians, corresponding to dots on a sphere
                
    lat_long_box : type = list - given in the format [ra_min, ra_max, dec_min, dec_max]
    
    radians (optional, defaults to False) : type = boolean - units for BOTH boundaries of box AND coordinates of dots
    

    Returns
    -------------
    boolean value True if coordinate is within box (inclusive of end values on both ends), False if not.
    '''
    if not radians:
        import numpy as np
        lat_long_box = np.radians(lat_long_box)
    ii = [((dot[0] >= lat_long_box[0]) & (dot[0] <= lat_long_box[1]) & 
          (dot[1] >= lat_long_box[2]) & (dot[1] <= lat_long_box[3])) for dot in dots]

    return ii #TJ return true or false array corresponding to True if a dot is inside the box

def populate_box(box, num_dots = 1000000, radians = True):
    '''populates a lat-long box with random dots by populating entire sphere, then filtering the ones outside of the box
    
    imports all standard packages automatically
    *requires the following functions from parent module "ASTRO5160/Homeworks/HW_2.py":
                        "evenly_populate_sphere()"
                        "is_it_in()"
    
    Parameters
    -------------
    box : type = list - describes box in format [ra_min, ra_max, dec_min, dec_max]
                
    num_dots (optional, defaults to 1 million) : type = int - number of dots on full sphere (not number of dots inside box)
    
    radians (optional, defaults to true) : type = boolean - False if output array of locations should be in degrees
    
    Returns
    -------------
    ra,dec for all the dots inside the lat-long box
    '''
    import numpy as np

    ra, dec = evenly_populate_sphere(num_dots) #TJ populate entire sphere with num_dots
    all_dots = np.column_stack((ra, dec))  #TJ store as a NumPy array for easier indexing

    #TJ generate list of dots that are in box of format specified in homework task 2
    in_box = all_dots[is_it_in(all_dots, box)]
    ra = [x[0] for x in in_box]
    dec = [x[1] for x in in_box]
    if not radians:
        return np.degrees(ra), np.degrees(dec)
    else:
        return ra,dec
    
    

def aitoff_plot_region(path, show_plot=True, save_plot=None, radians = False, name = 'bounded_region', 
                       label_with_area = False, title = "regions on sphere"):
    '''Plot regions bounded by path with vertices given in ra, dec in degrees (or radians if specified)
        example: box_path = [[0,0], [0,15], [20,15], [10,0]] will plot an assymmetrical trapezoidal box
    
    imports needed packages automatically
    *will need to manually import os before calling if save_plot = f'{os.getcwd()}'
    
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
    
    label_with_area (optional, defaults to False) : type = boolean - choose to label each box with the area in square degrees
                                                    **NOTE** area only works if this is a rectangular lat-long box.
    
    title (optional, defaults to "regions on sphere") : type = str - plot title
                For example, try:
                    import os
                    path1 = [[0,0], [0,15], [20,15], [10,0]]
                    path2 = [[30,30], [30,45], [45,30]]
                    aitoff_plot_region(path1, show_plot=False, save_plot=None, radians = False, 
                    name = 'region1', title = 'test plot')
                    aitoff_plot_region(path2, show_plot=False, save_plot=None, radians = False, 
                    name = 'region2')
                    aitoff_plot_region([], show_plot=True, save_plot=f'{os.getcwd()}', radians = False,
                    name = 'testing_function')
                    

    
    Returns
    -------------
    Nothing, just used to plot regions
    
    '''
        
    import matplotlib.pyplot as plt
    import numpy as np
    import os
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
            global num_plots
            num_plots = 0
            plt.figure(figsize=(8, 5))
            ax = plt.subplot(111, projection="aitoff")
            ax.set_xlabel("Right Ascension", fontsize=12)
            ax.set_ylabel("Declination", fontsize=12)
            ax.set_title(f"{title}", fontsize=12)
            ax.grid(color='black', linestyle='--', linewidth=0.5)
        num_plots+=1
        ax = plt.gca() #TJ if a plot already exists, then just use that current plot and only add things to it

        #TJ plot the path through all ra,dec corners
        if label_with_area:
            name = f'box {num_plots} area = {round(find_latlongrect_area(box),3)} sq deg'
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


if __name__ == "__main__":
    import argparse
    import os
    import numpy as np
    #TJ add argparse object and description
    parser = argparse.ArgumentParser(description="Populate sphere evenly, then return dots that are in given box.") 
    #TJ assign argument type and help note
    parser.add_argument("directory", type=str, help="directory to save output plots to", default = './') 

    args = parser.parse_args() #TJ assign arugments
    
    directory = args.directory #TJ assign variable to the provided argument
    
    #TJ create boxes that are progressively higher in declination and one thats a full hemisphere
    boxes = [[(5*15), (8*15), (0), (15)], [(5*15), (8*15), (20), (35)], 
         [(5*15), (8*15), (40), (55)], [(5*15), (8*15), (60), (75)],
         [0, 360, 0, 90]]
    for box in boxes: #TJ plot all boxes in same graph
        aitoff_plot_region(boundaries_to_vertices(box), show_plot=False, label_with_area=True, title = 'task1')
    #TJ name the final .png file and output to screen and save to directory
    aitoff_plot_region([], show_plot=True, save_plot = directory, name = 'progressively_higher_declination_boxes')
    
    box = [(5*15), (8*15), (20), (35)] #TJ define box to populate with dots
    num_dots = 1000000 #TJ set number of dots to populate entire sphere with
    ra, dec = populate_box(box, num_dots = num_dots) #TJ generate ra,dec values all over the sphere
    surface_area = 4*np.pi*(180/np.pi)**2 #TJ define entire surface area of sphere in square degrees
    print(f'expected number of dots within box based on area within : {round(num_dots*find_latlongrect_area(box)/(surface_area))}')
    print(f'actual number of dots observed in box : {len(ra)}')
    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




