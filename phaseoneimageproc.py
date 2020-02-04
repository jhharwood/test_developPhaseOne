#-------------------------------------------------------------------------------
# Name:         phaseoneimageproc.py
# Version 	    2.0
# Purpose:      This script contains functions to manipulate ExifLog information to produce
#		        kml/kmz, geoevents, copy template files, processes imagery, and processing file in general
#
#
# Requirements:  Anaconda Python 3.7.4 (standard library), pyproj
#
# Inputs:	
# Author:       J. Heath Harwood; Parts have been modified from Jon Sellars at NOAA DSS_quickThumbs64bit.py
#               See below changelog.
# References:   
#
# Created:     09/21/2018
# Copyright:   (c) USACE 2018
# Licence:     Public
#
# Change Log:
#       H. Harwood; V 1.0 Script is functional
# 	    H. Harwood; v 2.0
#		- added geoevent html with geojson tags for viewing in browser	
#-------------------------------------------------------------------------------

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.IMPORT STATEMENTs.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###


import os, re, glob, easygui, time, queue, threading, fnmatch
import numpy as np
import traceback
import shutil
from math import *
from pyproj import Proj
import zipfile


### Usage Constants
###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.CONST.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###
WEEK_OFFSET = 7.0 * 86400.0

imgtype = 'C'
imgelev = -23.5
imgyear = 2019

FL = 0.053  # FL= Focal Length in cm
CCD_X = 0.0540  # CCD_X = along track meters in m 0.006/9000
CCD_Y = 0.0404  # CCD_Y = across track meters Left Wing + in m
PP_X = 0.0078  # PP_X = Principle Point x(super)i from Terrestrial Cal Report/and or Boresite in m
PP_Y = -0.0042  # PP_Y = Principle Point y(super)i from Terrestrial Cal Report/and or Boresite in m
PP_Pix_Y = 3366  # PP_Pix_X = Principle Point x(super)P off from upper left in pixel
PP_Pix_X = 4500  # PP_Pix_Y = Principle Point y(super)P off from uppper left in pixel
CCD_XY = 0.0000060  # CCD_XY pixel size on CCD array meters in m

PP_Pix_X = PP_Pix_X
PP_Pix_Y = PP_Pix_Y

Geoid = 0  # Geoid ->Estimate for project area
Radius = 6378137  # Radius ->SemiMajor Radius of the Datum

# PhaseOne coarse navigation parameters
BS_P = 5
BS_R = 8
BS_H = -0.5

# Scale factor 1 for full scale
scale = 1
scale = scale + 0.00000
# SF = 1/scale

# Level of detail to turn layer on
minlod = 1
# set this so that management types can't zoom in beyond 1:1
maxlod = -1
clip = 0.0  # A percent used to clip the region to limit the number of images that pop up


###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.CONST.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.FUNCTIONS.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

### Define function for scanning flight data folder
### Function gets flight date folder as input and scans
### for the state/lake/area and grabs the block names in the area folder
def scanBlock(blk):
    """
    Function: Scans the block folder to get flight date information
    :param blk:
    :return:
            rawBlkPath: fields defined in GUI D:\2018_NCMP_GL\processing\_camera\LO\raw\20180908_2\blk_631
            blkName: block named assigned blk_631 or fnp_751
            areaName: area name ex. LO or MI
            flightDateNum: flight data folder ex. 20180908_2
            exifLog = exif log name ex ExifLog.csv
            camSysModel: ex. YC030284/MM000174
    """

    rawBlkPth = blk.replace('/', '\\')
    numEls = len(rawBlkPth.split('\\'))

    if numEls == 8:
        blkName = os.path.basename(rawBlkPth)
        areaName = rawBlkPth.split('\\')[4:5]
        flightDateNum = rawBlkPth.split('\\')[6:7]
        exifLogs = [os.path.join(rawBlkPth, exif)
                    for dirpath, dirnames, files in os.walk(rawBlkPth)
                    for exif in fnmatch.filter(files, 'ExifLog.csv')]
        sbetFiles = [os.path.join(rawBlkPth, sbet)
                     for dirpath, dirnames, files in os.walk(rawBlkPth)
                     for sbet in fnmatch.filter(files, 'sbet*.out')]

        exifLog = exif
        camSysModel = dirpath.split('\\')[8:9]

        if not sbetFiles:
            print ("No sbet file found.")
            return rawBlkPth, blkName, areaName, flightDateNum, exifLog, camSysModel
        else:
            sbetFile = sbet
            print ("Found sbet file " + sbetFile)
            return rawBlkPth, blkName, areaName, flightDateNum, exifLog, camSysModel, sbetFile

    else:
        # prints the associate error
        easygui.msgbox("Check raw path for correct setup, ex. D:\\2018_NCMP_GL\\processing\\_camera\LO\\raw\\20180908_2\\blk_631\n\n"
                       "'drive:\\effort\\processing\\_camera\\area\\raw\\YYYYMMDD_#\\blk_###'\n\n"
                       "Area needs to be two letters\n\n"
                       "Blocks can be either blk_* or fnp_* (always lowercase)",
                       "File Path Read Error")

### Define function to scan exif file to get first event time
### also checks for test events because we don't want the DC directory name
### to be based off test event times
def getEvt1Time(pth2exif):
    """Function readExif
           Read an exif file into a numpy array.
           Arguments:
                    filename: string of filename to read into a numpy array
           Returns: 2-d numpy array of exif data
           # Filename,GPS Event,GPS Time,Weeks:Seconds,Longitude,Latitude,Altitude,Pitch,Roll,Yaw,\
             Yaw Type,ISO,Aperture,Shutter Speed 1/#
        """
    try:
        print (pth2exif)
        if not isinstance(pth2exif, str):
            raise TypeError("argument 1 to readexif must be a string")

        data = np.genfromtxt(pth2exif, skip_header=1, dtype=None,
                             usecols=(2),
                             delimiter=',',
                             names=['gps_time'],
                             encoding='ascii',
                             invalid_raise=False)
        for event in data:
            #print event
            firstEvtTime = event[0].split(':')[0:2]
            break
        return firstEvtTime


    except Exception as e:
        print (e)
        # prints the associate error
        traceback.print_exc()
        # prints the associate error
        easygui.msgbox("An ExifLog.csv Error has occurred! Check log for any broken rows", "ExifLog Read Error")

def buildDSDir(dsDir):
    """
    Function to read in the dataset directory path and create a DC folder

    :param dsDir: Name of the dataset directory ex. D:\2018_NCMP_GL\processing\_camera\LO\fnp_751
    :return: Returns nothing
    """
    if os.path.isdir(dsDir):
        print ("Dataset directory already exists, moving on to next step\n")
        pass
    else:
        os.makedirs(dsDir)

def buildAreaDir(stDir):
    """
    Function to read in the dataset directory path and create a Area folder

    :param stDir: Name of the dataset directory ex. DC_DS_P_180908_1640
    :return: Returns nothing
    """
    if os.path.isdir(stDir):
        print ("Area directory already exists, moving on to next step\n")
        pass
    else:
        os.makedirs(stDir)

def buildBlkDir(blkDir):
    """
    Function to read in the dataset directory path and create a block folder

    :param stDir: Name of the dataset directory ex. D:\2018_NCMP_GL\processing\_camera\LO
    :return: Returns nothing
    """
    if os.path.isdir(blkDir):
        print ("Block directory already exists, moving on to next step\n")
        pass
    else:
        os.makedirs(blkDir)

### Define Function for retrieving EO Data
def readExif(filename):
    """Function readExif
       Read an exif file into a numpy array.
       Arguments:
                filename: string of filename to read into a numpy array
       Returns: 2-d numpy array of exif data
       # Filename,GPS Event,GPS Time,Weeks:Seconds,Longitude,Latitude,Altitude,Pitch,Roll,Yaw,\
         Yaw Type,ISO,Aperture,Shutter Speed 1/#
    """
    try:
        if not isinstance(filename, str):
            raise TypeError("argument 1 to readexif must be a string")

        data = np.genfromtxt(filename, skip_header=1, dtype=None,
                             usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                             delimiter=',',
                             names=['image', 'event_num', 'gps_time', 'gpsWeekSec', 'long', 'lat', 'alt', 'pitch',
                                    'roll', 'yaw'],
                             encoding='ascii',
                             invalid_raise=False)
        return data

    except Exception as e:
        print (e)
        # prints the associate error
        traceback.print_exc()
        easygui.msgbox("An ExifLog.csv Error has occurred! Check log for any broken rows", "ExifLog Read Error")

### Function to get picture time
def getPictureTime(yr, mth, day, gps_seconds, prevTime):
    """
    :Pass the gps seconds of the image to get the unix time (picture time):
    :return: pic time in Unix time
    """

    try:

        # Define flight date using the folder name and strip the time to create a time struct
        # flightDate = time.strptime((("%2s/%2s/%4s")%(mth,day,yr)), "%m/%d/%Y")
        flightDate = time.strptime((("%2s/%2s/%4s") % (mth, day, yr)), "%m/%d/%Y")

        # Get seconds from the epoch (01-01-1970) for the date in the filename. This will also give us the day of the
        # week for the GPS seconds of week calculation.
        tv_sec = time.mktime(flightDate)
        # print "GPS seconds from the epoch (01-01-1970) for the date in the filename :", tv_sec

        # Get flight date weekday
        flightDate_weekday = flightDate.tm_wday + 1
        if flightDate_weekday >= 7:
            flightDate_weekday = flightDate_weekday - 7
        # print "The flight day weekday is: ", flightDate_weekday
        # Subtract the number of days since Saturday midnight (Sunday morning) in seconds.
        tv_sec = tv_sec - (flightDate_weekday * 86400)
        start_week = tv_sec
        # print "Subtract the number of days since Saturday midnight (Sunday morning) in seconds. :", start_week
        # print "The number of GPS Seconds is: ", gps_seconds

        picture_time = (start_week + gps_seconds) * 1000000.0
        # print "The new Picture time is: ", str(int(picture_time))
        # print "The previous Picture time is: ", str(int(prevTime))

        # Test for GPS Week Rollover
        if int(picture_time) < int(prevTime):
            picture_time += WEEK_OFFSET * 1000000.0
            # print "The picture time after week rollover is: %s" % str(int(picture_time))
            return str(int(picture_time))
        else:
            # print "The picture time before week rollover is: %s" % str(int(picture_time))
            return str(int(picture_time))


    except IOError as e:  # prints the associate error
        print ("I/O Error %d: %s" % (e.errno, e.strerror))
        raise e

# Function for creating IPAS Files
def createIPAS(eoData, year, month, day, evtFile, photoID, photoIDPix4D, frameSer, coarseDat):
    """
    Function to create the eo data file for Ipas CO+

    :param eoData: accepts the 2-D numpy array
    :param year: 4 Digit year ex 2018
    :param month: 2 Digit month 09
    :param day: 2 Digit Day 08
    :param evtFile: name and path of the event file for Ipas CO+
    :param photoID: name and path of the photo id file for Ipas CO+
    :param frameSer: name and path of the frameser file for Ipas CO+
    :param coarseDat: name and path of the coarse dat file for Ipas CO+ Boresight Computation
    :return: Returns nothing
    """
    # run line by line to get image name and gps time
    try:

        for row in eoData:
            fileName = row[0].strip('.IIQ')
            iName = os.path.join(fileName + '.tif')
            iNameJpg = os.path.join(fileName + '.jpg')
            gpsWeek = row[3].split(':')[0]
            gpsTime = row[3].split(':')[1]
            nameTime = iName.split('_')[2:3]
            lon = float(row[4])
            lat = float(row[5])
            alt = row[6].split(' ')[0]
            roll = row[8]
            pitch = row[7]
            yaw = row[9]
            # print 'Name Time is ', nameTime
            timeImage = ''.join(nameTime)
            timePattern = timeImage.replace('-', '').replace('.', '')
            # print 'Image time is ', timePattern

            # Create data for evt, photo id, and frame serial number files
            evtLines = ('%s 62065%s%s%s%s \n') % (gpsTime, year, month, day, timePattern)
            # print evtLines
            photoIDLines = ('%s %s \n') % (gpsTime, iName)
            photoIDPix4DLines = ('%s %s \n') % (gpsTime, iNameJpg)
            # print photoIDLines
            frameSerLines = ('%s 62065 62065%s%s%s%s %s \n') % (gpsTime, year, month, day, timePattern, iName)
            # print frameSerLines
            coarseDatLines = ('62065%s%s%s%s %s %s %s %s %s %s %s \n') % (year, month, day, timePattern, iName, lon, lat, alt, roll, pitch, yaw)
            #print coarseDatLines

            # Write data to files
            evtFile.writelines(evtLines)
            photoID.writelines(photoIDLines)
            photoIDPix4D.writelines(photoIDPix4DLines)
            frameSer.writelines(frameSerLines)
            coarseDat.writelines(coarseDatLines)

        evtFile.close()
        photoID.close()
        photoIDPix4D.close()
        frameSer.close()
        coarseDat.close()

    except IOError as e:  # prints the associate error
        print ("I/O Error %d: %s" % (e.errno, e.strerror))
        raise e
 # Create frame and camera tables
def createFCTables(eoData,camTableCsv):
    """
        Function to create the eo data file for Ipas CO+

        :param eoData: accepts the 2-D numpy array
        :param tifPath: filepath to the DC directory/tifs
        :param camTableCsv: name and path of the camera (IO file) table for ArcMap/ArcGIS Pro Frame Camera Table script or OrthoMapping
        :return: Returns nothing

        camera table format
        CameraID,CameraMaker,CameraModel,FocalLength,PrincipalX,PrincipalY,NRows,NCols,PixelSize,FilmFiducials,\
        DistortionType,Radial,Tangential,RadialDistances,RadialDistortions
        """
    # run line by line to get image name and gps time
    try:

        # Create the header list for the camera csv file
        camerafieldnames = ['OBJECTID','CameraID','CameraMaker','CameraModel','FocalLength','PrincipalX','PrincipalY',
                            'NRows','NCols','PixelSize','FilmFiducials','DistortionType','Radial','Tangential',
                            'RadialDistances','RadialDistortions']

        with open(camTableCsv, 'wb') as fc:
            fcWriter = csv.DictWriter(fc, fieldnames=camerafieldnames)
            fcWriter.writeheader()

            i = 0
            for row in eoData:
                objectID = i
                cameraID = row[1]
                cameraMaker = 'Phase One'
                cameraModel = 'iXM-RS150F'
                focalLength = '51371.2'
                PrincipalX = '58.4'
                PrincipalY = '-23.9'
                NRows = '10652'
                NCols = '14204'
                PixelSize = '3.76'
                FilmFiucials = ''
                DistortionType = 'DistortionModel'
                Radial = '0.0;1.605600e-5;-6.189660e-9;1.313290e-12'
                Tangential = '-2.0419e-6;-4.51202e-6'
                RadialDistance = ''
                RadialDistortion = ''

                camDataTable = [objectID,cameraID,cameraMaker,cameraModel,focalLength,PrincipalX,PrincipalY,NRows,
                                NCols,PixelSize,FilmFiucials,DistortionType,Radial,Tangential,RadialDistance,
                                RadialDistortion]
                fcWriter.writerow(dict(zip(camerafieldnames,camDataTable)))

                i+=1

        fc.close()
        print ("Done writing the formated camera table for ArcGIS")


    except IOError as e:  # prints the associate error
        print ("I/O Error %d: %s" % (e.errno, e.strerror))
        raise e

# Function for create the KML and KMZ files
def processKML(eoData, name2, kml, year, month, day, prev_time, zoneNum, zoneHem, cameraSyncR):
    """
    Function to create the KML/KMZ from the EO data

    :param eoData: accepts the 2-D numpy array
    :param name2: Dataset name ex DC_DS_P_180908_1359
    :param kml: name and path of Dataset KML file
    :param year: 4 Digit year ex 2018
    :param month: 2 Digit month 09
    :param day: 2 Digit Day 08
    :param prev_time: unix time stamp of previous image
    :param zoneNum: 2 digit zone number string
    :param zoneHem: 1 letter zone hemisphere string
    :param cameraSyncR: name and path of Dataset camera sync file
    :return: Returns nothing
    """
    # Write out the KML header
    khead = """<?xml version="1.0" encoding="UTF-8"?>
	<kml xmlns="http://earth.google.com/kml/2.2">
	<Document>
	  <Style id="FEATURES_LABELS">
		<IconStyle>
		  <color>FFFFFFFF</color>
		  <Icon>
			<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>
		  </Icon>
		</IconStyle>
		<LabelStyle>
		  <color>FFFFFFFF</color>
		</LabelStyle>
	  </Style>
	  <Folder>
		<name>""" + name2 + """ Thumbnails</name>\n"""

    kml.writelines(khead)

    # Initialize a line counter
    counter = 0
    # Begin looping through

    for row in eoData:

        # parse the exiflog to make the event, photoid, frameserial, kml, and camera sync files
        pid = row[1]
        iiqName = row[0]
        stripName = row[0].strip('.IIQ')
        tifName = os.path.join(stripName + '.tif')
        newFileName = re.sub(r'[-|.]', r'_', stripName)
        iName = os.path.join(newFileName + '_scaled.jpeg')
        lat = float(row[5])
        latF = format(lat, '.11f')
        lon = float(row[4])
        lonF = format(lon, '.11f')
        alt = float(row[6].split(' ')[0])
        altF = format(alt, '.4f')
        pitch = float(row[7])
        pitchF = format(pitch, '.13f')
        roll = float(row[8])
        rollF = format(roll, '.13f')
        yaw = float(row[9])
        yawF = format(yaw, '.13f')
        gpsWeek = row[3].split(':')[0]
        gpsSeconds = float(row[3].split(':')[1])
        gpsSecondsF = format(gpsSeconds, '.5f')

        nameTime = iName.split('_')[2:3]
        # print 'Name Time is ', nameTime
        timeImage = ''.join(nameTime)
        timePattern = timeImage.replace('-', '').replace('.', '')
        # print 'Image time is ', timePattern

        # Convert Lat and Lon to radians
        Latrad = radians(lat)
        Lonrad = radians(lon)

        # Convert Lat and Lon to meters
        X = (Radius * (cos(Latrad))) * Lonrad
        Y = Radius * Latrad

        # Calculate pixel size based on flying height
        FH = alt
        PixSz = (FH * CCD_XY) / FL

        # Calculate pixel size in degrees for X and Y diminsions
        PixSzY = (degrees(PixSz / Radius))
        PixSzX = degrees(PixSz / (Radius * (cos(Latrad))))

        # Convert degress to radians (using * !DTOR <pi/180 ~= 0.01745>) ****Note H2 is set above******
        Prad = radians(pitch + BS_P)
        Rrad = radians(roll + BS_R)
        Hrad = radians(yaw + BS_H)

        # Calculate s_inv
        s_inv = 1 / (FL / FH)

        # Create terms for the M matrix
        M11 = cos(Prad) * sin(Hrad)
        M12 = -cos(Hrad) * cos(Rrad) - sin(Hrad) * sin(Prad) * sin(Rrad)
        M13 = cos(Hrad) * sin(Rrad) - sin(Hrad) * sin(Prad) * cos(Rrad)
        M21 = cos(Prad) * cos(Hrad)
        M22 = sin(Hrad) * cos(Rrad) - (cos(Hrad) * sin(Prad) * sin(Rrad))
        M23 = (-sin(Hrad) * sin(Rrad)) - (cos(Hrad) * sin(Prad) * cos(Rrad))
        M31 = sin(Prad)
        M32 = cos(Prad) * sin(Rrad)
        M33 = cos(Prad) * cos(Rrad)
        # print M11, M12, M13, M21, M22, M23, M31, M32, M33

        # Define p matrix using PP offsets (image center) (rotate to +X direction of flight along track, +Y left wing across track, -FL)
        Xi = PP_Y
        Yi = -1.00 * PP_X
        FLneg = -1 * FL

        # s_inv * M * p + T(GPSxyz)
        CP_X = (s_inv * (M11 * Xi + M12 * Yi + M13 * FLneg)) + X
        CP_Y = (s_inv * (M21 * Xi + M22 * Yi + M23 * FLneg)) + Y
        CP_Z = (s_inv * (M31 * Xi + M32 * Yi + M33 * FLneg)) + FH

        # Calculate Upper left corner (from center) in mapping space (DIR FLT = +Y, Right Wing = +X), rotate, apply to center coords in mapping space
        ULX = (PixSz * PP_Pix_X) + CP_X
        ULY = (PixSz * PP_Pix_Y) + CP_Y
        LRX = CP_X - (PixSz * PP_Pix_X)
        LRY = CP_Y - (PixSz * PP_Pix_Y)

        # Convert CP_X and CP_Y to Degrees Lat and - Long
        # East and West appear switched but it is the way X and Y
        # are handeled in the image reference frame above
        East = (degrees(ULX / (Radius * cos(Latrad))))  # - 360
        North = (degrees(ULY / Radius))
        West = (degrees(LRX / (Radius * cos(Latrad))))  # - 360
        South = (degrees(LRY / Radius))

        # Calculated center point of lat and Long convereed to degrees
        CenLat = (degrees(CP_Y / Radius))
        CenLon = (degrees(CP_X / (Radius * cos(Latrad))))

        # Calculate Rotation
        if imgtype == "C":
            Rot = yaw
            if Rot >= 180: Rotation = 360.00 - Rot
            if Rot < 180: Rotation = -1.00 * (Rot * 1.00)
        if imgtype == "R":
            Rot = yaw + 180.0
            if Rot > 360: Rot = Rot - 360
            if Rot < 180: Rotation = 360 - Rot
            if Rot >= 180: Rotation = -1 * (Rot * 1)

        # Calculate the clip
        nsClip = (North - South) * clip
        ewClip = (West - East) * clip

        # Write to the master KML
        kbody = """  	<Document>
		<name>""" + iName + """</name>
		<Region>
			<LatLonAltBox>
								<north>""" + str(North - nsClip) + """</north>
								<south>""" + str(South + nsClip) + """</south>
								<east>""" + str(East + ewClip) + """</east>
								<west>""" + str(West - ewClip) + """</west>
				<minAltitude>0</minAltitude>
				<maxAltitude>0</maxAltitude>
			</LatLonAltBox>
			<Lod>
				<minLodPixels>""" + str(minlod) + """</minLodPixels>
				<maxLodPixels>""" + str(maxlod) + """</maxLodPixels>
				<minFadeExtent>0</minFadeExtent>
				<maxFadeExtent>0</maxFadeExtent>
			</Lod>
		</Region>
		<Style>
			<ListStyle id="hideChildren">
				<listItemType>checkHideChildren</listItemType>
				<bgColor>00ffffff</bgColor>
			</ListStyle>
		</Style>
		<GroundOverlay>
			<name>r08569772</name>
			<Icon>
				<href>""" + iName + """</href>
			</Icon>
			<LatLonBox>
						<north>""" + str(North) + """</north>
						<south>""" + str(South) + """</south>
						<east>""" + str(East) + """</east>
						<west>""" + str(West) + """</west>
						<rotation>""" + str(Rotation) + """</rotation>
			</LatLonBox>
		</GroundOverlay>
	</Document>\n"""
        # print kbody
        kml.writelines(kbody)

        # Increment counter
        counter += 1

    kmiddle = """</Folder>
		<Folder>
			<visibility>0</visibility>
		<name>""" + name2 + """ Labels</name>\n"""
    kml.writelines(kmiddle)

    # # Write state to log
    # kmlLoopLabelInfoMsg = datetime.now().strftime('%m\%d\%Y %H:%M:%d')+"    Starting KML Label Loop...\n"
    # logFile.writelines(kmlLoopLabelInfoMsg)
    # print kmlLoopLabelInfoMsg

    counter = 0
    for row in eoData:

        # parse the exiflog to make the event, photoid, frameserial, kml, and camera sync files
        pid = row[1]
        iiqName = row[0]
        stripName = row[0].strip('.IIQ')
        tifName = os.path.join(stripName + '.tif')
        newFileName = re.sub(r'[-|.]', r'_', stripName)
        iName = os.path.join(newFileName + '_scaled.jpeg')
        lat = float(row[5])
        latF = format(lat, '.11f')
        lon = float(row[4])
        lonF = format(lon, '.11f')
        alt = float(row[6].split(' ')[0])
        altF = format(alt, '.4f')
        pitch = float(row[7])
        pitchF = format(pitch, '.13f')
        roll = float(row[8])
        rollF = format(roll, '.13f')
        yaw = float(row[9])
        yawF = format(yaw, '.13f')
        gpsWeek = row[3].split(':')[0]
        gpsSeconds = float(row[3].split(':')[1])
        gpsSecondsF = format(gpsSeconds, '.5f')
        picTime = getPictureTime(year, month, day, gpsSeconds, prev_time)
        prev_time = picTime

        nameTime = iName.split('_')[2:3]
        # print 'Name Time is ', nameTime
        timeImage = ''.join(nameTime)
        timePattern = timeImage.replace('-', '').replace('.', '')
        # print 'Image time is ', timePattern

        # Convert Lat and Lon to radians
        Latrad = radians(lat)
        Lonrad = radians(lon)

        # Convert Lat and Lon to meters
        X = (Radius * (cos(Latrad))) * Lonrad
        Y = Radius * Latrad

        # Calculate pixel size based on flying height
        FH = alt
        PixSz = (FH * CCD_XY) / FL

        # Calculate pixel size in degrees for X and Y diminsions
        PixSzY = (degrees(PixSz / Radius))  # * 1.10 WHAT WAS THIS FOR?????
        PixSzYneg = 1.00 * (degrees(PixSz / Radius))
        PixSzX = degrees(PixSz / (Radius * (cos(Latrad))))

        # Convert degress to radians (using * !DTOR <pi/180 ~= 0.01745>) ****Note H2 is set above******
        Prad = radians(pitch + BS_P)
        Rrad = radians(roll + BS_R)
        Hrad = radians(yaw + BS_H)

        # Calculate s_inv
        s_inv = 1 / (FL / FH)

        # Create terms for the M matrix
        M11 = cos(Prad) * sin(Hrad)
        M12 = -cos(Hrad) * cos(Rrad) - sin(Hrad) * sin(Prad) * sin(Rrad)
        M13 = cos(Hrad) * sin(Rrad) - sin(Hrad) * sin(Prad) * cos(Rrad)
        M21 = cos(Prad) * cos(Hrad)
        M22 = sin(Hrad) * cos(Rrad) - (cos(Hrad) * sin(Prad) * sin(Rrad))
        M23 = (-sin(Hrad) * sin(Rrad)) - (cos(Hrad) * sin(Prad) * cos(Rrad))
        M31 = sin(Prad)
        M32 = cos(Prad) * sin(Rrad)
        M33 = cos(Prad) * cos(Rrad)

        # Define p matrix using PP offsets (image center) (rotate to +X direction of flight along track, +Y left wing across track, -FL)
        Xi = PP_Y
        Yi = -1.00 * PP_X
        FLneg = -1 * FL

        CP_X = (s_inv * (M11 * Xi + M12 * Yi + M13 * FLneg)) + X
        CP_Y = (s_inv * (M21 * Xi + M22 * Yi + M23 * FLneg)) + Y
        CP_Z = (s_inv * (M31 * Xi + M32 * Yi + M33 * FLneg)) + FH

        # Calculate Upper left corner (from center) in mapping space (DIR FLT = +Y, Right Wing = +X), rotate, apply to center coords in mapping space
        ULX = CP_X - (PixSz * PP_Pix_X)
        ULY = (PixSz * PP_Pix_Y) + CP_Y
        LRX = CP_X - (PixSz * PP_Pix_X)
        LRY = CP_Y - (PixSz * PP_Pix_Y)

        # Convert CP_X and CP_Y to Degrees Lat and - Long
        # East and West appear switched but it is the way X and Y
        # are handeled in the image reference frame above
        ULXLon = (degrees(ULX / (Radius * cos(Latrad))))  # - 360
        ULYLat = (degrees(ULY / Radius))
        LRXLon = (degrees(LRX / (Radius * cos(Latrad))))  # - 360
        LRYLat = (degrees(LRY / Radius))

        CenLat = (degrees(CP_Y / Radius))
        CenLon = (degrees(CP_X / (Radius * cos(Latrad))))

        # Convert center Lat, Long to UTM
        p = Proj(proj='utm', zone=zoneNum, ellps='WGS84')
        cenEast, cenNorth = p(lon, lat)
        cenEastF = format(cenEast, '.11f')
        cenNorthF = format(cenNorth, '.11f')

        # Convert center Lat, Long to UTM
        p = Proj(proj='utm', zone=zoneNum, ellps='WGS84')
        ulE, ulN = p(ULXLon, ULYLat)

        # Use rotation to get params for 0 to 180 Degrees of yaw
        Rot = 0 - yaw
        ulEp = cenEast + ((ulE - cenEast) * cos(radians(Rot)) - (((ulN - cenNorth)) * sin(radians(Rot))))
        ulNp = cenNorth + ((((ulE - cenEast)) * sin(radians(Rot))) + (((ulN - cenNorth)) * cos(radians(Rot))))

        Easting = ulEp
        Northing = ulNp

        # Calculate the center of the upper left pixel
        cenULEst = cenEast - (PP_Pix_Y * PixSz)
        cenULNrt = cenNorth + ((PP_Pix_X * PixSz) * 2)

        # Calculate Rotation
        Rotation = 0
        if imgtype == "C":
            Rot = yaw
            # print "The Rotation for the nav file is: " + str(Rot)
            if Rot > 0 and Rot < 180:
                Rotation = Rot
                # print "Rotation is from 0 to 180: " + str(Rotation)
                # Calculate Rotation for pixels of upper left for tfw
                PixSzNeg = -1 * PixSz
                lineA = PixSz * cos((pi / 180) * Rotation)
                lineD = PixSzNeg * sin((pi / 180) * Rotation)
                lineB = PixSzNeg * sin((pi / 180) * Rotation)
                lineE = PixSzNeg * cos((pi / 180) * Rotation)

            else:
                Rot = 360 + yaw
                # print "The Correct Rotation for the nav file is: " + str(Rot)
                Rotation = Rot
                # Calculate Rotation for pixels of upper left for tfw
                PixSzNeg = -1 * PixSz
                lineA = PixSz * cos((pi / 180) * Rotation)
                lineD = PixSzNeg * sin((pi / 180) * Rotation)
                lineB = PixSzNeg * sin((pi / 180) * Rotation)
                lineE = PixSzNeg * cos((pi / 180) * Rotation)

        if imgtype == "R":
            Rot = yaw + 180.0
            if Rot > 360: Rot = Rot - 360
            if Rot < 180: Rotation = 360 - Rot
            if Rot >= 180: Rotation = -1 * (Rot * 1)

        # Calculate the clip
        nsClip = (North - South) * clip
        ewClip = (West - East) * clip

        # Write to the master KML
        kbody2 = """  	<Placemark>
				<visibility>0</visibility>
				<name>""" + iName + """</name>
				<styleUrl>#FEATURES_LABELS</styleUrl>
				<Point>
					<extrude>0</extrude>
					<altitudeMode>absolute</altitudeMode>
					<coordinates>""" + str(CenLon) + """,""" + str(CenLat) + """,0</coordinates>
				</Point>
				</Placemark>\n"""
        kml.writelines(kbody2)

        # Contructs the string list for use in the Original Camera Sync file associated with HF and
        # the RCD30 Camera Sync file with that same formatting; 0 file can be used in
        # camLines0 = (' %-20s  %60s  %-13s  %-13s  %-10s  -1  -1  -1  %s  %s  %-14s  %-15s  %-13s  %16s  %16s  %-14s  \n')%\
        # (pid, iName, cenEastF, cenNorthF, altF, zoneNum, zoneHem, latF, lonF, gpsSecondsF, rollF, pitchF, yawF)
        camLinesR = (
                        ' %-20s  %60s  %-13s  %-13s  %-10s  -1  -1  -1  %s  %s  %-14s  %-15s  %-13s  %16s  %16s  %-14s  %-16s  \n') % \
                    (pid, iName, cenEastF, cenNorthF, altF, zoneNum, zoneHem, latF, lonF, gpsSecondsF, rollF, pitchF,
                     yawF, picTime)
        #print camLinesR

        # Write the new strings to the new file
        cameraSyncR.writelines(camLinesR)

        # Increment counter
        counter += 1

    # Finish up kml creation
    kfoot = """</Folder>
	</Document>
	</kml>"""
    kml.writelines(kfoot)

def createGeoJson(eoData, tifPath, dcDirName):
    """
    Function to create the geojson events file in html and returns the geolocation of the block

    :param eoData: 2D numpy array of eodata
    :param tifPath: name and path of DC directory ex. D:\2018_NCMP_GL\processing\_camera\LO\blk_631\DC_DS_P_180908_1359
    :param dcDirName: Dataset directory name ex. DC_DS_P_180908_1359
    :return: lat and long of the block of data
    """
    outFileHandle = open(tifPath + '\\' + dcDirName + '.geojson', 'w')

    # the head of the geojson file
    header = \
        ''' 
{ 
   "type": "FeatureCollection",
   "features": ['''

    outFileHandle.write(header)

    # the template. where data from the csv will be formatted to geojson
    template = '''
  { 
    "type": "Feature",
    "geometry": {
       "type": "Point",
       "coordinates": [ %s, %s ]
    },
    "properties": {
    "Filename":"%s",
    "Latitude":"%s",
    "Longitude":"%s",
    "GPS Event":%s,
    "GPS Time":"%s",
    "Weeks:Seconds":"%s",
    "Pitch":%s,
    "Roll":%s,
    "Yaw":%s
    }
  }%s'''
    numRecs = len(eoData)
    comma = ','
    space = ' '
    iter = 0
    for row in eoData:
        filename = str(row[0])
        gps_event = str(row[1])
        gps_time = row[2]
        week_sec = row[3]
        long = float(row[4])
        lat = float(row[5])
        pitch = row[7]
        roll = row[8]
        yaw = row[9]
        if iter < numRecs - 1:
            output = template % (
                str(long), str(lat), str(filename), str(lat), str(long), str(gps_event), str(gps_time), str(week_sec), str(pitch), str(roll),
                str(yaw), comma)
            # print output
            outFileHandle.write(output)
        else:
            output = template % (
                str(long), str(lat), str(filename), str(lat), str(long), str(gps_event), str(gps_time), str(week_sec), str(pitch), str(roll),
                str(yaw), space)
            # print output
            outFileHandle.write(output)
        iter += 1

    # the tail of the geojson file
    footer = ''' 
]
}'''

    # opens an geoJSON file to write the output
    outFileHandle.write(footer)
    outFileHandle.close()

    #return outFileHandle
    return lat, long
    # print outFileHandle

def getNumImgs(imgDir):
    """
    Function to get the number of images in the raw image/camera directory
    :param imgDir: Name and path to the raw image directory ex.
    D:\2018_NCMP_GL\processing\_camera\LO\raw\20180908_2\blk_631\YC030284
    :return: Returns the total number of files in the images directory
    """
    # Get images from image directory
    # imgs = glob.glob('*.IIQ')
    imgs = [os.path.join(imgDir, f)
            for dirpath, dirnames, files in os.walk(imgDir)
            for f in fnmatch.filter(files, '*.IIQ')]
    totFiles = len(imgs)
    print (totFiles)
    return totFiles

def getImgs(imgDir):
    """
    Function to get the raw images (IIQ) files in the raw image directory for processing
    :param imgDir: name and path to the raw image direcotry ex.
    D:\2018_NCMP_GL\processing\_camera\LO\raw\20180908_2\blk_631\YC030284
    :return: Returns a list of raw images to be processed for the block
    """
    # Get images from image directory
    # imgs = glob.glob('*.IIQ')
    imgs = [os.path.join(imgDir, f)
            for dirpath, dirnames, files in os.walk(imgDir)
            for f in fnmatch.filter(files, '*.IIQ')]
    return imgs

def processImgs(imgs, tifPath, copeSwCmd, irfSwCmd, nCopeProc, nIrfanProc, kml, name2, cameraSyncR, evtFile, photoID,
                frameSer,choice,cope,irfanview):
    """
    Function to process the images in COPE and Irfanview
    :param imgs: Raw image files to be processed
    :param tifPath: name and path to the DC directory ex.
    D:\2018_NCMP_GL\processing\_camera\LO\blk_631\DC_DS_P_180908_1359
    :param copeSwCmd: list to store cope commands in
    :param irfSwCmd: list to store irfanview commands in
    :param nCopeProc: interger for number of CPUs to use for processing
    :param nIrfanProc: interger for number of CPUs to use for processing
    :param kml: name of KML file to store points in
    :param name2: DC directory name ex. DC_DS_P_180908_1359
    :param cameraSyncR: name and file path of camera sync file
    :param evtFile: name and file path of event file
    :param photoID: name and file path of photoid file
    :param frameSer: name and file path of frameser file
    :param choice: Choice from GUI that could be one of 4 choices: 'Tifs and Scaled Jpegs',
    'Tifs Only,'Scaled Jpegs Only','Rerun the EO/KMZ'
    :return: Returns nothing
    """

    if choice == "Tifs and Scaled Jpegs":
        # Runs the jpeg scaling based on choices above
        # Yes it will run the jpegs, no it passes and goes ahead and creates the kmz (though there may be nothing
        # in the kmz, best to only run if the jpegs are already developed.

        # Begin processing loop
        for img in imgs:
            imgName = os.path.basename(img)
            imgOut = imgName.strip('.IIQ') + '.tif'
            newJpeg = re.sub(r'[-|.]', r'_', imgName)
            jpgOut = newJpeg.strip('_IIQ') + '_scaled.jpeg'
            #print jpgOut

            procCopeCmd = (r'%s %s %s\%s -outputformat=tif -resolution=300 -resolutionunit=inch -bits=8 -Rotation=0.0 -brightness=0.0 -contrast=0.0 -saturation=0.0 -levelHighlight=1.0 -levelShadow=0.0 -targetHighlight=1.0 -targetShadow=0.0 -highlightRecovery=0.0 -shadowRecovery=0.0 -clarity=0.0 -enableOpenCL=1') \
                          % (cope, img, tifPath, imgOut)
            irfSw = (r'%s %s\%s /resize=(800,600) /resample /aspectratio /jpeg=100 /convert=%s\%s') \
                    % (irfanview, tifPath, imgOut, tifPath, jpgOut)
            #print procCopeCmd
            #print irfSw
            copeSwCmd.append(procCopeCmd)
            irfSwCmd.append(irfSw)

        queue = Queue.Queue()

        class threadJobs(threading.Thread):
            def __init__(self, queue):
                threading.Thread.__init__(self)
                self.queue = queue

            def run(self):
                while True:
                    # get job from queue
                    myJob = self.queue.get()
                    inCmd = myJob
                    os.system(inCmd)

                    # signal queue job is done
                    self.queue.task_done()

        def cope():
            # spawn a pool of processes, and let the queue manage them
            for i in range(nCopeProc):
                t = threadJobs(queue)
                t.setDaemon(True)
                t.start()
            # populate queue with jobs
            for cmd in copeSwCmd:
                queue.put(cmd)
            # wait on the queue until everything has been processed
            queue.join()
            # wait for processing to finish
            while not queue.empty():
                time.sleep(1)
            return 'Done'

        def irfan():
            # spawn a pool of processes, and let the queue manage them
            for i in range(nIrfanProc):
                t = threadJobs(queue)
                t.setDaemon(True)
                t.start()
            # populate queue with jobs
            for cmd in irfSwCmd:
                queue.put(cmd)
            # wait on the queue until everything has been processed
            queue.join()
            # wait for processing to finish
            while not queue.empty():
                time.sleep(1)
            return 'Done'

        cope()
        irfan()

    elif choice == "Jpgs and Scaled Jpegs":
        # Runs the jpeg scaling based on choices above
        # Yes it will run the jpegs, no it passes and goes ahead and creates the kmz (though there may be nothing
        # in the kmz, best to only run if the jpegs are already developed.

        # Begin processing loop
        for img in imgs:
            imgName = os.path.basename(img)
            imgOut = imgName.strip('.IIQ') + '.jpg'
            newJpeg = re.sub(r'[-|.]', r'_', imgName)
            jpgOut = newJpeg.strip('_IIQ') + '_scaled.jpeg'
            #print jpgOut

            procCopeCmd = (r'%s %s %s\%s -outputformat=jpg -jpgquality=75 -resolution=300 -resolutionunit=inch -bits=8 -Rotation=0.0 -brightness=0.0 -contrast=0.0 -saturation=0.0 -levelHighlight=1.0 -levelShadow=0.0 -targetHighlight=1.0 -targetShadow=0.0 -highlightRecovery=0.0 -shadowRecovery=0.0 -clarity=0.0 -enableOpenCL=1') \
                          % (cope, img, tifPath, imgOut)
            irfSw = (r'%s %s\%s /resize=(800,600) /resample /aspectratio /jpeg=100 /convert=%s\%s') \
                    % (irfanview, tifPath, imgOut, tifPath, jpgOut)
            #print procCopeCmd
            #print irfSw
            copeSwCmd.append(procCopeCmd)
            irfSwCmd.append(irfSw)

        queue = Queue.Queue()

        class threadJobs(threading.Thread):
            def __init__(self, queue):
                threading.Thread.__init__(self)
                self.queue = queue

            def run(self):
                while True:
                    # get job from queue
                    myJob = self.queue.get()
                    inCmd = myJob
                    os.system(inCmd)

                    # signal queue job is done
                    self.queue.task_done()

        def cope():
            # spawn a pool of processes, and let the queue manage them
            for i in range(nCopeProc):
                t = threadJobs(queue)
                t.setDaemon(True)
                t.start()
            # populate queue with jobs
            for cmd in copeSwCmd:
                queue.put(cmd)
            # wait on the queue until everything has been processed
            queue.join()
            # wait for processing to finish
            while not queue.empty():
                time.sleep(1)
            return 'Done'

        def irfan():
            # spawn a pool of processes, and let the queue manage them
            for i in range(nIrfanProc):
                t = threadJobs(queue)
                t.setDaemon(True)
                t.start()
            # populate queue with jobs
            for cmd in irfSwCmd:
                queue.put(cmd)
            # wait on the queue until everything has been processed
            queue.join()
            # wait for processing to finish
            while not queue.empty():
                time.sleep(1)
            return 'Done'

        cope()
        irfan()

    elif choice == "Tifs Only":
        # Runs the jpeg scaling based on choices above
        # Yes it will run the jpegs, no it passes and goes ahead and creates the kmz (though there may be nothing
        # in the kmz, best to only run if the jpegs are already developed.

        # Write state to log
        #tifFilesInfoMsg = datetime.now().strftime('%m\%d\%Y %H:%M:%d') + "    Started process for tifs.\n"
        #logFile.writelines(tifFilesInfoMsg)
        #print tifFilesInfoMsg

        # Begin processing loop
        for img in imgs:
            imgName = os.path.basename(img)
            imgOut = imgName.strip('.IIQ') + '.tif'

            procCopeCmd = (r'%s %s %s\%s -outputformat=tif -resolution=300 -resolutionunit=inch -bits=8 -Rotation=0.0 -brightness=0.0 -contrast=0.0 -saturation=0.0 -levelHighlight=1.0 -levelShadow=0.0 -targetHighlight=1.0 -targetShadow=0.0 -highlightRecovery=0.0 -shadowRecovery=0.0 -clarity=0.0 -enableOpenCL=1') \
                          % (cope, img, tifPath, imgOut)
            #print procCopeCmd
            copeSwCmd.append(procCopeCmd)

        queue = Queue.Queue()

        class threadJobs(threading.Thread):
            def __init__(self, queue):
                threading.Thread.__init__(self)
                self.queue = queue

            def run(self):
                while True:
                    # get job from queue
                    myJob = self.queue.get()
                    inCmd = myJob
                    os.system(inCmd)

                    # signal queue job is done
                    self.queue.task_done()

        def cope():
            # spawn a pool of processes, and let the queue manage them
            for i in range(nCopeProc):
                t = threadJobs(queue)
                t.setDaemon(True)
                t.start()
            # populate queue with jobs
            for cmd in copeSwCmd:
                queue.put(cmd)
            # wait on the queue until everything has been processed
            queue.join()
            # wait for processing to finish
            while not queue.empty():
                time.sleep(1)
            return 'Done'

        cope()
        #tifFinishedInfoMsg = datetime.now().strftime('%m\%d\%Y %H:%M:%d') + "    Finshed Processing Tif files.\n"
        #logFile.writelines(tifFinishedInfoMsg)

    elif choice == "Scaled Jpegs Only":
        # Runs the jpeg scaling based on choices above
        # Yes it will run the jpegs, no it passes and goes ahead and creates the kmz (though there may be nothing
        # in the kmz, best to only run if the jpegs are already developed.

        # Write state to log
        #jpegFilesInfoMsg = datetime.now().strftime('%m\%d\%Y %H:%M:%d') + "    Started process for jpegs.\n"
        #logFile.writelines(jpegFilesInfoMsg)
        #print jpegFilesInfoMsg

        # Begin processing loop
        for img in imgs:
            imgName = os.path.basename(img)
            imgOut = imgName.strip('.IIQ') + '.tif'
            newJpeg = re.sub(r'[-|.]', r'_', imgName)
            jpgOut = newJpeg.strip('_IIQ') + '_scaled.jpeg'

            irfSw = (r'%s %s\%s /resize=(800,600) /resample /aspectratio /jpeg=100 /convert=%s\%s') \
                    % (irfanview, tifPath, imgOut, tifPath, jpgOut)
            #print irfSw
            irfSwCmd.append(irfSw)

        queue = Queue.Queue()

        class threadJobs(threading.Thread):
            def __init__(self, queue):
                threading.Thread.__init__(self)
                self.queue = queue

            def run(self):
                while True:
                    # get job from queue
                    myJob = self.queue.get()
                    inCmd = myJob
                    os.system(inCmd)

                    # signal queue job is done
                    self.queue.task_done()

        def irfan():
            # spawn a pool of processes, and let the queue manage them
            for i in range(nIrfanProc):
                t = threadJobs(queue)
                t.setDaemon(True)
                t.start()
            # populate queue with jobs
            for cmd in irfSwCmd:
                queue.put(cmd)
            # wait on the queue until everything has been processed
            queue.join()
            # wait for processing to finish
            while not queue.empty():
                time.sleep(1)
            return 'Done'

        irfan()
        #irFanFinishedInfoMsg = datetime.now().strftime(
        #    '%m\%d\%Y %H:%M:%d') + "    Finshed Processing scaled jpeg files.\n"
        #logFile.writelines(irFanFinishedInfoMsg)

    elif choice == "Rerun the EO/KMZ":
        pass
        # msg = "Do you want to continue?  Have the tifs been created for this dataset?"
        # title = "Please Confirm"
        # choices = ["Yes", "No"]
        # if easygui.buttonbox(msg, title, choices=choices):
        #     pass
        # else:
        #     sys.exit(0)

    kml.flush()
    kml.close()

    kmzName = os.path.join(tifPath + "\\" + '_' + name2 + "_thumbs.kmz")
    #print kmzName
    kmz = zipfile.ZipFile(kmzName, "w")

    if os.name == 'nt':
        # change directories to scoop up the scaled jpegs and kml
        os.chdir(tifPath)
        for name in glob.glob("*_scaled.jpeg"):
            kmz.write(name, os.path.basename(name), zipfile.ZIP_DEFLATED)
        for name in glob.glob("*_thumbs.kml"):
            kmz.write(name, os.path.basename(name), zipfile.ZIP_DEFLATED)

    kmz.close()
    cameraSyncR.flush()
    cameraSyncR.close()
    evtFile.close()
    photoID.close()
    frameSer.close()

def copySysFiles(camInstall,tifPath,supportPath):

    if os.path.exists(supportPath) and camInstall == 'MM000174_20190325':
       shutil.copy2(supportPath + r'\phase_one\parameters\CZ04_Install_MM000174_20190325_20190514\_ipas_files\Installation_settings_RCD30_62065.ini',\
                    '%s\\Installation_settings_RCD30_62065.ini' % tifPath)
       shutil.copy2(supportPath + r'\phase_one\parameters\CZ04_Install_MM000174_20190325_20190514\\_ipas_params\_CZ04_Install_MM000174_20190325_20190514_utm16n_wgs84_template.ptu',\
                    '%s\\_CZ04_Install_MM000174_20190325_20190514_utm16n_wgs84_template.ptu' % tifPath)
       shutil.copy2(supportPath + r'\phase_one\parameters\CZ04_Install_MM000174_20190325_20190514\\_ipas_params\_CZ04_Install_MM000174_20190325_20190514_utm17n_wgs84_template.ptu', \
           '%s\\_CZ04_Install_MM000174_20190325_20190514_utm17n_wgs84_template.ptu' % tifPath)
       shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.bat','%s\\_eo2c3d_v1.2.bat' % tifPath)
       shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.py', '%s\\_eo2c3d_v1.2.py' % tifPath)
       shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.bat', '%s\\_buildC3D_v1.1.bat' % tifPath)
       shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.py', '%s\\_buildC3D_v1.1.py' % tifPath)

    elif os.path.exists(supportPath) and camInstall == 'MM000174_20190607':
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_MM00147_20190607_20190904\_ipas_files\Installation_settings_RCD30_62065.ini', \
            '%s\\Installation_settings_RCD30_62065.ini' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_MM00147_20190607_20190904\_ipas_params\_CZ04_Install_MM00147_20190607_20190904_utm16n_wgs84_template.ptu', \
            '%s\\_CZ04_Install_MM00147_20190607_20190904_utm16n_wgs84_template.ptu' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_MM00147_20190607_20190904\_ipas_params\_CZ04_Install_MM00147_20190607_20190904_utm17n_wgs84_template.ptu', \
            '%s\\_CZ04_Install_MM00147_20190607_20190904_utm17n_wgs84_template.ptu' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.bat', '%s\\_eo2c3d_v1.2.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.py', '%s\\_eo2c3d_v1.2.py' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.bat', '%s\\_buildC3D_v1.1.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.py', '%s\\_buildC3D_v1.1.py' % tifPath)
        
    elif os.path.exists(supportPath) and camInstall == 'MM000134_20190905':
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_MM00134_20190905_20191111\_ipas_files\Installation_settings_RCD30_62065.ini', \
            '%s\\Installation_settings_RCD30_62065.ini' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_MM00134_20190905_20191111\_ipas_params\_CZ04_Install_MM00134_20190905_20191111_utm17n_wgs84_template.ptu', \
            '%s\\_CZ04_Install_MM00134_20190905_20191111_utm17n_wgs84_template.ptu' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_MM00134_20190905_20191111\_ipas_params\_CZ04_Install_MM00134_20190905_20191111_utm18n_wgs84_template.ptu', \
            '%s\\_CZ04_Install_MM00134_20190905_20191111_utm18n_wgs84_template.ptu' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.bat', '%s\\_eo2c3d_v1.2.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.py', '%s\\_eo2c3d_v1.2.py' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.bat', '%s\\_buildC3D_v1.1.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.py', '%s\\_buildC3D_v1.1.py' % tifPath)

    elif os.path.exists(supportPath) and camInstall == 'YC030333_20181115':
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ02_Install_YC30333_20181115_20190608\_ipas_files\Installation_settings_RCD30_62065.ini', \
            '%s\\Installation_settings_RCD30_62065.ini' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ02_Install_YC30333_20181115_20190608\\_ipas_params\_CZ02_Install_YC30333_20181115_20190608_utm16n_wgs84_template.ptu', \
            '%s\\_CZ02_Install_YC30333_20181115_20190608_utm16n_wgs84_template.ptu' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ02_Install_YC30333_20181115_20190608\\_ipas_params\_CZ02_Install_YC30333_20181115_20190608_utm17n_wgs84_template.ptu', \
            '%s\\_CZ02_Install_YC30333_20181115_20190608_utm17n_wgs84_template.ptu' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.bat', '%s\\_eo2c3d_v1.2.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.py', '%s\\_eo2c3d_v1.2.py' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.bat', '%s\\_buildC3D_v1.1.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.py', '%s\\_buildC3D_v1.1.py' % tifPath)

    elif os.path.exists(supportPath) and camInstall == 'YC030333_20200113':
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ02_Install_YC30333_20200113_Pres\_ipas_files\Installation_settings_RCD30_62065.ini', \
            '%s\\Installation_settings_RCD30_62065.ini' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ02_Install_YC30333_20200113_Pres\\_ipas_params\_CZ02_Install_YC30333_20200113_Pres_utm17n_wgs84_template.ptu', \
            '%s\\_CZ02_Install_YC30333_20200113_Pres_utm17n_wgs84_template.ptu' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ02_Install_YC30333_20200113_Pres\\_ipas_params\_CZ02_Install_YC30333_20200113_Pres_utm18n_wgs84_template.ptu', \
            '%s\\_CZ02_Install_YC30333_20200113_Pres_utm18n_wgs84_template.ptu' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.bat', '%s\\_eo2c3d_v1.2.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.py', '%s\\_eo2c3d_v1.2.py' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.bat', '%s\\_buildC3D_v1.1.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.py', '%s\\_buildC3D_v1.1.py' % tifPath)

    elif os.path.exists(supportPath) and camInstall == 'YC030284_20190514':
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_YC030284_20190514_Pres\_ipas_files\Installation_settings_RCD30_62065.ini', \
            '%s\\Installation_settings_RCD30_62065.ini' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_YC030284_20190514_Pres\_ipas_params\_CZ04_Install_YC030284_20190514_Pres_utm16n_wgs84_template.ptu', \
            '%s\\_CZ04_Install_YC030284_20190514_Pres_utm16n_wgs84_template.ptu' % tifPath)
        shutil.copy2(
            supportPath + r'\phase_one\parameters\CZ04_Install_YC030284_20190514_Pres\_ipas_params\_CZ04_Install_YC030284_20190514_Pres_utm17n_wgs84_template.ptu', \
            '%s\\_CZ04_Install_YC030284_20190514_Pres_utm17n_wgs84_template.ptu' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.bat', '%s\\_eo2c3d_v1.2.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\eo2c3d\_eo2c3d_v1.2.py', '%s\\_eo2c3d_v1.2.py' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.bat', '%s\\_buildC3D_v1.1.bat' % tifPath)
        shutil.copy2(supportPath + r'\scripts\buildC3D\_buildC3D_v1.1.py', '%s\\_buildC3D_v1.1.py' % tifPath)

    else:
        # prints the associate error
        easygui.msgbox("Please check connection to the Support folder, This script relies on paths: \n"
                       "D:\\Support\\phase_one and D:\\Support\\scripts", "File Path Read Error")
