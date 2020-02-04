#-------------------------------------------------------------------------------
# Name:        readsbet.py
# Version 	   1.1
# Purpose:     This script parses the sbet and uses the event time from the exiflog to figure out the exact position in the event
# Input:       sbet.out and exiflog.csv
#
# Requirements:  Anaconda Python 2.7 (standard library), pyproj
#
# Inputs:	
# Author:       J. Heath Harwood
#        
# References:   
#
# Created:     09/21/2018
# Copyright:   (c) USACE 2018
# Licence:     Public
#
# Change Log:
#       H. Harwood; V 1.1 Script is functional
#
#-------------------------------------------------------------------------------

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.IMPORT STATEMENTs.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

import numpy as np
import os, glob, sys, easygui, math, csv
import pandas as pd

# change to the cwd and get cwf
cwd = os.getcwd()
os.chdir(cwd)
cwf = os.path.basename(cwd)


def find_nearest(array,value):
    value0 = value - 0.0025
    value1 = value + 0.0025
    #print value0, value1
    idx = (np.abs(array['time'] - value)).argmin()
    idx0 = (np.abs(array['time'] - value0)).argmin()
    idx1 = (np.abs(array['time'] - value1)).argmin()
    #print idx, idx0, idx1
    if idx > idx0 and idx == idx1:
        recordAfter = array[idx]
        recordBefore = array[idx0]
        #print array[idx]
        #print array[idx0]
        return recordBefore,recordAfter
    else:
        recordBefore = array[idx]
        recordAfter = array[idx1]
        #print array[idx]
        #print array[idx1]
        return recordBefore,recordAfter

###########################################################################
# Get a list of the SBET record types
# This is the definition of a SBET record
###########################################################################
def sbet_record_types():
    """Function sbet_record_types
       Get a list of the sbet record types
       Arguments:
       Returns: list of the data types in a sbet record
    """
    return [ ("time", np.float64),
             ("lat", np.float64),
             ("lon", np.float64),
             ("alt", np.float64),
             ("ewspeed", np.float64),
             ("nsspeed", np.float64),
             ("vertspeed", np.float64),
             ("roll", np.float64),
             ("pitch", np.float64),
             ("heading", np.float64),
             ("wander", np.float64),
             ("ewacc", np.float64),
             ("nsacc", np.float64),
             ("vertacc", np.float64),
             ("xacc", np.float64),
             ("yacc", np.float64),
             ("zacc", np.float64) ]

###########################################################################
# Read a sbet file into a numpy array.
# This is the function that reads the SBET file and returns the data as a
# numpy array
###########################################################################
def readSbet(filename):
    """Function readSbet
       Read an sbet file into a numpy array.
       Arguments:
                filename: string of filename to read into a numpy array
       Returns: 2-d numpy array of sbet data
    """
    if not isinstance(filename, str):
        raise TypeError("argument 1 to readSbet must be a string")
    sbetData = np.fromfile(filename, dtype=np.dtype(sbet_record_types()))
    return sbetData

