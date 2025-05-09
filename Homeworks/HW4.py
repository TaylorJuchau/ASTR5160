#!/usr/bin/env python
# coding: utf-8

# In[4]:


from astropy.table import Table
import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..'))) #TJ change directory to include the entire ASTRO5160 directory
from Homeworks.useful_functions import *

def splendid_function(table):
    """Determines which objects in table are likely to be qsos. Expected ~95% of qsos will be correctly identified as qsos,
    expected ~1% of non-qsos to be erroneously selected as a qso. Note for LARGE tables, it is likely that MOST of these objects
    are not qsos, but MOST of the qsos in the table will be correctly identified as qsos. In other words, false positives are 
    significantly more likely than false negatives.

    Parameters
    -------------
    table :  type = astropy.table - table of objects
    
    Returns
    -------------
    True/False array with length equal to rows in table, True means the object was identified as a qso
        
    """

    cut_values= [1.67838290e+12, 3.11068380e+12, 1.42528763e+01, 1.98868179e+01,
           1.29941151e-01] #TJ these cuts were determined by a machine learning algorithm I wrote
    cut_columns = ['RA_IVAR', 'DEC_IVAR', 'FLUX_W1', 'FLUX_W2', 
           'GAIA_ASTROMETRIC_SIGMA5D_MAX'] 
    #TJ these columns were selected by that same machine learning algorithm that tried many columns and combinations of columns
    objects = trim_high_R_mag(table) #TJ remove objects that are dim in R_band
    kept_objects = trim_columns(objects,cut_columns) #TJ trim the table of all the extra columns so it is faster to work with
    tally_array = [
        (np.array([obj[col] for col in cut_columns]) > cut_values).astype(int).tolist()
        for obj in kept_objects
    ] #TJ this is an array of [1, 1, 1, 1, 1] if it is flagged as a qso for all columns
    return [sum(row)>4 for row in tally_array] #TJ return true only if all of the columns passed the benchmark

def trim_columns(table, kept_columns):
    """Trim a table with a large number of columns to only include columns we need, this is for faster processing.

    Parameters
    -------------
    table :  type = astropy.table - table of objects, must include columns with titles in the kept_columns argument
    kept_columns :  type = list - list of strings corresponding to the column names
    
    Returns
    -------------
    original table with only columns in the kept_columns list
        
    """
    return table[kept_columns]

def trim_high_R_mag(table):
    """Trim a table of objects to only include objects with R_mag of less than 19

    Parameters
    -------------
    table :  type = astropy.table - table of objects, must include "FLUX_R" column
    
    Returns
    -------------
    original table with any objects with R_magnitude higher than 19 removed 
    """
    new_table = add_mag_column(table, ["FLUX_R"]) #TJ add the magnitude column for R
    return new_table[new_table['R_mag'] <19] #TJ delete any rows that have R_magnitudes larger than 19


if __name__ == "__main__":
    import argparse
    from astropy.table import Table
    #TJ add argparse object and description
    parser = argparse.ArgumentParser(description="""Uses cuts in the following columns (RA_IVAR, DEC_IVAR, FLUX_W1, FLUX_W2, 
    and GAIA_ASTROMETRIC_SIGMA5D_MAX) to determine which of the objects in a given fits file are determined to be qsos. These 
    cuts and  columns were determined by a machine learning algorithm that can be found by asking Taylor (it's stored locally 
    on his computer and was not pushed to github. Note for LARGE tables, it is likely that MOST of these objects
    are not qsos, but MOST of the qsos in the table will be correctly identified as qsos. In other words, false positives are 
    significantly more likely than false negatives.""") 
    #TJ assign argument type and help note
    parser.add_argument("path", type=str, help="path to fits file") 

    args = parser.parse_args() #TJ assign arugments
    
    path = args.path #TJ assign variable to the provided argument
    print(f'reading in data from {path}.')
    table = Table.read(path) #TJ read in fits file
    print("Approximate number of qsos in fits file:", sum(splendid_function(table))) #TJ print the total number of objects flagged as qso

