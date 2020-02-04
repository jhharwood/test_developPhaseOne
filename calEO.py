#-------------------------------------------------------------------------------
# Name:        calEO.py
# Version:     1.5
# Purpose:     This script called from the phaseoneimageryproc.py.  
# Input:        ExifLog and Sbet.out
#
#
# Requirements:  Anaconda Python 2.7 (standard library), read_sbet.py, pyproj
#
#
# Author:       J. Heath Harwood; Modified from Jon Sellars DSS_quickThumbs64bit.py
#               See above changelog.
# References:   - New Calibration and Computing Method for Direct Georeferencing of Image and Scanner Data Using
#               the Position and Angular Data of an Hybrid Inertial Navigation System; Baumker adn Heimes, 2001
#               https://www.hochschule-bochum.de/fileadmin/media/fb_v/labore/photogrammetrie/Artikel/Veroeffentlichungen/OEEPE_Symposium_Hannover_2001.PDF
#               - ESTIMATION OF ANGLE ELEMENTS OF EXTERIOR ORIENTATION FOR UAV
#                 IMAGES BASED ON INS DATA AND AERIAL TRIANGULATION PROCESSING; Wierzbicki 1990
#               http://www.tf.llu.lv/conference/proceedings2018/Papers/N054.pdf
#               - The ellipsoid and the Transverse Mercator projection; Geodetic Information Paper No 1 2/1998 (version 2.2)
#               http://fgg-web.fgg.uni-lj.si/~/mkuhar/Zalozba/TM_projection.pdf
#               - http://mathworld.wolfram.com/EulerAngles.html
#               - http://www.cmlab.csie.ntu.edu.tw/~jsyeh/vision2ta/chapter14/
#               - New Calibration and Computing Method for Direct Georeferencing of Image and Scanner Data Using (1).pdf
#
# Created:     09/21/2018
# Copyright:   (c) USACE 2018
# Licence:     Public
#
# Change Log:
#       H. Harwood; Script is functional
#-------------------------------------------------------------------------------

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.IMPORT STATEMENTs.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###
import re
import os
from math import *
from datetime import *
import time
from pyproj import Proj
import numpy as np
import read_sbet

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.CONSTANTS.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

#Geoid ->Estimate for project area
Geoid = 0
# Scale factor on Central meridian
F0 = 0.99963209
#Semi-Major Axis-A Equitorial Radius of the WGS84 Datum in meters
a = 6378137.0
aF0 = a * F0
#Semi-Minor Axis-B Polar Radius of the WGS84 Datum in meters
b = 6356752.3142
bF0 = b * F0
# radius of earth
radius = 6371000 # meters
# eccentricity
e2 = (aF0**2 - bF0**2) / aF0**2
#e2 = (a**2 - b**2) / a**2

#E = (1 - (b**2/a**2))**0.5
#e = sqrt((a**2 - b**2) / a**2)
#e2 = e**2
e2sqr = e2**2
#print "Eccentricity " + str(e2sqr)
# Ellipsoiodal constant ratio
n = (a - b) / (a + b)
#print str(n)
# Define the utm zone
#zone = 19



###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.FUNCTIONS.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

# Function to calculate the event lat/long based on time, speed and distance between the two sbet events
def calc_coords(a,lat1,lat2,lon1,lon2,evtTime,t_sbet1,t_sbet2):

    # convert decimal degrees to radians
    #lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # calculate the bearing between to the two points from sbet
    yPos = sin(lon2 - lon1) * cos(lat2)
    xPos = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon2 - lon1)
    bearing = atan2(yPos, xPos)
    # Calculate Bearing
    # bearingRad = math.atan2(math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(lon2 - lon1),
    #                         math.sin(lon2 - lon1) * math.cos(lat2))
    # bearing = math.degrees(bearingRad)
    # print "This is the bearing " + str(bearing)

    # calculate the distance between the two points from sbet using haversine
    lat = lat2 - lat1
    lon = lon2 - lon1
    hav = sin(lat * 0.5)**2 + cos(lat1) * cos(lat2) * sin(lon * 0.5)**2
    pt_sbet_dist = 2 * a * asin(sqrt(hav))
    # print "This is distance between sbet events " + str(pt_sbet_dist) + " meters"
    # # New haversine formula
    # dlat = lat2 - lat1
    # dlon = lon2 - lon1
    #
    # hav = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    # dist = 2 * asin(sqrt(hav))
    # pt_sbet_dist = a * dist
    #print "This is distance between sbet events " + str(pt_sbet_dist) + " meters"

    # calculate the speed it takes to get from sbet point 1 to sbet point 2

    curSpd = pt_sbet_dist / (t_sbet2 - t_sbet1)
    #print "This is the current speed " + str(curSpd) + " m/s"

    # calculate the distance to the event point based on current speed and event time
    dist2pt = curSpd * (evtTime - t_sbet1)
    #print "This is the distance to the event point " + str(dist2pt) + " meters"

    # get the lat and long of the event based on the bearing, lat1, lon2, and distance to the event
    newLat = asin(sin(lat1)*cos(dist2pt/a) + cos(lat1)*sin(dist2pt/a)*cos(bearing))
    newLon =  lon1 + atan2(sin(bearing)*sin(dist2pt/a)*cos(lat1), cos(dist2pt/a) - sin(lat1)*sin(lat2))
    #print newLat, newLon

    return newLat,newLon

### Define function for reading EO data and sbet data
### and created an ouput eo file with the correcte omage, phi, kappa values

def calEO(exifData, sbetFile, tifPath, dcDirName, zone):
    """

    :param exifData:
    :param sbetFile:
    :param tifPath:
    :param dcDirName:
    :return:
    """
    # Create EO out file
    eoTxtFileName = tifPath + '\\' + dcDirName + '_eo_for_c3d_GRS80.txt'
    eoTxtFile = open(eoTxtFileName, 'w')

    # run line by line to get gps time
    try:
        start = time.time()
        count_run_times = 0
        # Get the sbet and read the data into a numpy array - want to do this here so we can have access
        # to the array and not in the loop so we are only ready the sbet array in to memory once
        recordData = read_sbet_v1_1.readSbet(sbetFile)
        for row in exifData:
            #gpsWeek = row[3].split(':')[0]
            imageName = row[0].strip('.IIQ') + '.tif'
            imageNum = row[1]
            gpsTime = row[3].split(':')[1]
            gpsTimeFlt = np.float64(gpsTime)
            #print gpsTimeFlt
            record1,record2 = read_sbet_v1_1.find_nearest(recordData, gpsTimeFlt)
            #print "This is record before GPS Time: " + str(record1)
            #print "This is record after GPS Time: " + str(record2)
            #print "Got records"

            # Definitions for Records 1 and 2
            t_sbet1 = record1[0]
            #print t_sbet1
            lat1 = degrees(record1[1])
            lon1 = degrees(record1[2])
            #print lat1,lon1
            alt1 = record1[3]
            roll1 = record1[7]
            pitch1 = record1[8]
            heading1 = record1[9]

            t_sbet2 = record2[0]
            #print t_sbet2
            lat2 = degrees(record2[1])
            lon2 = degrees(record2[2])
            #print lat2,lon2
            alt2 = record2[3]
            roll2 = record2[7]
            pitch2 = record2[8]
            heading2 = record2[9]

            # Interpolate the positions from sbet times
            dif_t_ratio = (t_sbet2 - t_sbet1 / t_sbet2 - t_sbet1)
            #print "Time diff ratio " + str(dif_t_ratio)
            deltTime = (t_sbet2 - t_sbet1) / 2
            #print deltTime
            newTime = t_sbet1 + deltTime
            #print newTime
            deltaLat = lat2 - lat1
            deltaLon = lon2 - lon1
            deltaRoll = roll2 - roll1
            deltaPitch = pitch2 - pitch1
            deltaHeading = heading2 - heading1
            deltaAlt = alt2 - alt1

            oldLat = lat1 + ((deltaLat/deltTime) * ((t_sbet2 - t_sbet1)/(dif_t_ratio/0.65)))
            oldLon = lon1 + ((deltaLon/deltTime) * ((t_sbet2 - t_sbet1)/(dif_t_ratio/0.65)))
            #print oldLat, oldLon

            newCoords = calc_coords(aF0, radians(lat1), radians(lat2), radians(lon1), radians(lon2), gpsTimeFlt, t_sbet1, t_sbet2)
            newLat = degrees(newCoords[0])
            newLon = degrees(newCoords[1])
            #print newLat, newLon

            newRoll = roll1 + deltaRoll #* dif_t_ratio
            #newRoll = ((roll2 - roll1)/(t_sbet2 - t_sbet1)) * ((t_sbet2 - t_sbet1)/dif_t_ratio) + roll1
            #print newRoll
            newPitch = pitch1 + deltaPitch #* dif_t_ratio
            #newPitch = ((pitch2 - pitch1) / (t_sbet2 - t_sbet1)) * ((t_sbet2 - t_sbet1)/dif_t_ratio) + pitch1
            newHeading = heading1 + deltaHeading #* dif_t_ratio
            #newHeading = ((heading2 - heading1) / (t_sbet2 - t_sbet1)) * ((t_sbet2- t_sbet1)/dif_t_ratio) + heading1
            newAlt = alt1 + deltaAlt #* dif_t_ratio
            #newAlt = (alt1 + alt2) * abs(dif_t_ratio) / 2
            #newAlt = ((alt2 - alt1) / (t_sbet2 - t_sbet1)) * ((t_sbet2- t_sbet1)/dif_t_ratio) + alt1
            #print newRoll, newPitch, newHeading, newAlt

            # Lever Arms in Degrees
            laDegX = -0.743
            laDegY = 0.07
            laDegZ = 0.071
            # Boresight Angles
            misBx = radians(0.0341121)
            misBy = radians(0.12524213)
            misBz = radians(-0.01981901)

            # Lever Arms in Radians
            laRadX = laDegX
            laRadY = laDegY
            laRadZ = laDegZ

            # Test settings in Radians
            latRad = radians(newLat)
            lonRad = radians(newLon)
            rollRad = newRoll + misBx
            pitchRad = newPitch + misBy
            headingRad = newHeading + misBz
            #print latRad, lonRad, rollRad, pitchRad, headingRad

            # Latitude angles
            sinLatRad = sin(latRad)
            cosLatRad = cos(latRad)
            tanLatRad = tan(latRad)

            # Longitude angles
            sinLonRad = sin(lonRad)
            cosLonRad = cos(lonRad)
            tanLonRad = tan(lonRad)

            # TM_Projections by Ordance Survey Geodetic information paper No. 1 2/1998 (version 2.2) document equations
            # v = radius of curvature at latitude  perpendicular to a meridian
            nu = aF0 / sqrt(1 - (e2sqr * sinLatRad ** 2))
            #nu = aF0 / (1 - e2sqr * sinLatRad**2)**0.5
            # p (rho) = radius of curvature of a meridian at latitude phi
            rho = aF0 * (1 - e2sqr) / (1 - (e2sqr * sinLatRad ** 2)) ** 1.5
            # n2 (eta2) = v/p - 1 or nu over rho - 1
            eta2 = (nu / rho) - 1
            # P = longitude at point measured east (+) or west (-) of Greenwich minus
            #     longitude of true origin of central meridian
            cm = (abs(zone - 30) * 6) + 3
            # cm = (zone - 31) * 6 + 3
            #print "Central Meridian " + str(cm)
            P = lonRad - radians(cm)
            #print "P equals " + str(P)

            # Convergence calculation C from Phi (latitdue and Lamda (longitude)
            # TM_Projections by Ordance Survey Geodetic information paper No. 1 2/1998 (version 2.2) document equations
            # XIII = sinLatRad
            XIV = ((-sinLatRad * cosLatRad ** 2) / 3) * (1 + 3 * eta2 + 2 * eta2 ** 4)
            # print degrees(XIV)
            XV = ((-sinLatRad * cosLatRad ** 4) / 15) * (2 - tanLatRad ** 2)
            # print degrees(XV )

            # Convergence from phi and lambda
            # TM_Projections by Ordance Survey Geodetic information paper No. 1 2/1998 (version 2.2) document equations
            C = radians((P * (sinLatRad)) + (P ** 3 * XIV) + (P ** 5 * XV))
            #print "Convergence: " + str(C)
            headingCor = headingRad + C
            #print "Heading Correction " + str(headingCor)

            pitchCos = cos(pitchRad)
            pitchSin = sin(pitchRad)
            rollCos = cos(rollRad)
            rollSin = sin(rollRad)
            headCorCos = cos(headingCor)
            headCorSin = sin(headingCor)

            # Matrix m terms rph to wpk
            m11 = headCorCos * rollCos
            # m12 = headCorSin * pitchCos
            # m13 = -pitchSin
            m21 = pitchSin * rollSin * headCorCos + pitchCos * headCorSin
            # m22 = pitchSin * rollSin * headCorSin + pitchCos * headCorCos
            # m23 = rollCos * pitchSin
            m31 = pitchCos * rollSin * headCorCos - pitchSin * headCorSin
            m32 = -pitchCos * rollSin * headCorSin - pitchSin * headCorCos
            m33 = pitchCos * rollCos

            # Euler Angles
            phi = asin(m31)
            omega = atan(-m32 / m33)
            kappa = atan2(-m21, m11)

            # Print test values
            #print ("Omega, Phi, Kappa from IPAS CO+ EO: 0.125100000000,1.90880000000,30.49918000000")
            #print ("Omega, Phi, Kappa from Cam Sync EO: 0.251720000000,1.9283760000,30.27496600000")
            #print ("Omega, Phi, Kappa from Cam DatM EO: 0.225442000000,2.05513100000,30.46123300000")
            #print ("Omega, Phi, Kappa from Spt File EO: ") + str(degrees(omega)), str(degrees(phi)), str(degrees(kappa))

            # # Add lever arm transforamtion to roll, pitch and heading
            # # From ALB Blue Book II Manual eq. 4.1.4
            rollRad_laX = laRadX * (cos(pitchRad) * cos(headingRad)) + laRadY * (
            sin(rollRad) * sin(pitchRad) * cos(headingRad) - cos(rollRad) * sin(headingRad)) + laRadZ * (
            cos(rollRad) * sin(pitchRad) * cos(headingRad) + sin(rollRad) * sin(headingRad))
            pitchRad_laY = laRadX * (cos(pitchRad) * sin(headingRad)) + laRadY * (
            sin(headingRad) * sin(pitchRad) * sin(headingRad) + cos(rollRad) * cos(headingRad)) + laRadZ * (
            cos(rollRad) * sin(pitchRad) * sin(headingRad) - sin(rollRad) * cos(headingRad))
            headingRad_laZ = -laRadX * (sin(pitchRad)) + laRadY * (sin(rollRad) * cos(pitchRad)) + laRadZ * (
            cos(rollRad) * cos(pitchRad))
            # # print "Roll, Pitch, Heading with LA in radians: " + str(rollRad_laX),str(pitchRad_laY),str(headingRad_laZ)
            # # print "Roll, Pitch, Heading with LA in degrees: " + str(degrees(rollRad_laX)),str(degrees(pitchRad_laY)),str(degrees(headingRad_laZ))

            # # Get Roll Pitch and Heading from GPS/IMU to ECEF Ref Frame
            # # From ALB Blue Book II Manual eq. 4.1.11
            rollRad_laX_gpsimu_ecef = -rollRad_laX * (sinLatRad * cosLonRad) - pitchRad_laY * (
            sinLonRad) - headingRad_laZ * (cosLatRad * cosLonRad)
            pitchRad_laY_gpsimu_ecef = -rollRad_laX * (sinLatRad * sinLonRad) + pitchRad_laY * (
            cosLonRad) - headingRad_laZ * (cosLatRad * sinLonRad)
            headingRad_laZ_gpsimu_ecef = rollRad_laX * (cosLatRad) - headingRad_laZ * (sinLatRad)
            # print "Roll, Pitch, Heading with LA from IMU/GPS in radians: " + str(rollRad_laX_gpsimu_ecef),str(pitchRad_laY_gpsimu_ecef),str(headingRad_laZ_gpsimu_ecef)
            # print "Roll, Pitch, Heading with LA from IMU/GPS in degrees: " + str(degrees(rollRad_laX_gpsimu_ecef)),str(degrees(pitchRad_laY_gpsimu_ecef)),str(degrees(headingRad_laZ_gpsimu_ecef))

            # Get geodetic position of the ECEF coordinate of the ellipsoid height (altitude as input)
            # From ALB Blue Book II Manual eq. 4.1.12
            # Meridian radius of curvature M
            # M = a*(1 - E**2) / ((1 - E**2 * sinLatRad**2)**1.5)
            # print "M = " + str(M)
            # Tranverse Radius N
            N = a/(1 - e2**2) / ((1 - e2**2 * sinLatRad**2)**0.5)
            rN = aF0 / sqrt(1 - (e2sqr * sinLatRad**2))
            #print "rN = " + str(rN)
            # e2sqr = e2**2

            # print e2sqr
            # N = 0.9
            rollRad_laX_gpsimu_ecef_hgt = (rN + newAlt) * (cosLatRad * cosLonRad) + rollRad_laX_gpsimu_ecef
            pitchRad_laY_gpsimu_ecef_hgt = (rN + newAlt) * (cosLatRad * sinLonRad) + pitchRad_laY_gpsimu_ecef
            headingRad_laZ_gpsimu_ecef_hgt = (rN * (1 - e2sqr) + newAlt) * sinLatRad + headingRad_laZ_gpsimu_ecef
            # print ("XYZ ") + str(rollRad_laX_gpsimu_ecef_hgt),str(pitchRad_laY_gpsimu_ecef_hgt),str(headingRad_laZ_gpsimu_ecef_hgt)
            # print "Roll, Pitch, Heading with LA from IMU/GPS above ellipsoid in radians: " + str(rollRad_laX_gpsimu_ecef_hgt),str(pitchRad_laY_gpsimu_ecef_hgt),str(headingRad_laZ_gpsimu_ecef_hgt)

            # Compute the coordinates of the camera centroid latitude
            # Reference http://www.oc.nps.edu/oc2902w/coord/coordcvt.pdf, Section B ECEF xyz to Lat, Lon, Height
            p = sqrt(rollRad_laX_gpsimu_ecef_hgt ** 2 + pitchRad_laY_gpsimu_ecef_hgt ** 2)
            z = headingRad_laZ_gpsimu_ecef_hgt
            r = sqrt(rollRad_laX_gpsimu_ecef_hgt ** 2 + pitchRad_laY_gpsimu_ecef_hgt ** 2 + headingRad_laZ_gpsimu_ecef_hgt** 2)
            # Compute the coordinates of the camera centroid longitude
            #outLon = atan(pitchRad_laY_gpsimu_ecef_hgt / rollRad_laX_gpsimu_ecef_hgt)
            outLon = atan2(pitchRad_laY_gpsimu_ecef_hgt,rollRad_laX_gpsimu_ecef_hgt)
            # print "Out Longitude " + str(degrees(outLon))
            inLat = atan(z / p)  #* (1/sqrt(1 - (e2sqr * sinLatRad**2)))
            #inLat = atan2(p,z)  # * (1/sqrt(1 - (e2sqr * sinLatRad**2)))
            # print "In Latitude " + str(degrees(inLat))
            nowLat = inLat
            # print "Now Latitude " + str(degrees(nowLat))
            #RN = aF0 / sqrt(1 - (e2sqr * sin(nowLat) ** 2))
            #RN = a / sqrt(1 - (e2sqr * sin(nowLat) ** 2))
            loop = 0
            while loop <= 3:
                h = (p / cos(nowLat)) - (rN * nowLat)
                # print "Height is equal to " + str(h)
                #nextLat = atan( (z/p) * (1 - e2sqr * ( RN / (RN + h ) ) ) )
                nextLat = atan((z / p) * (1 + (aF0 * sin(nowLat) * e2sqr) / z * sqrt(1 - e2sqr * sin(nowLat)**2)))
                #nextLat = atan((headingRad_laZ_gpsimu_ecef_hgt/sqrt(rollRad_laX_gpsimu_ecef_hgt**2 + pitchRad_laY_gpsimu_ecef_hgt**2)) * (1/sqrt(1 - (e2sqr * sinLatRad**2))))
                # print "Next Latitude is " + str(degrees(nextLat))
                loop += 1
                #print nextLat
                nowLat = nextLat
            outLat = nextLat
            # Compute the corrected ellipsoid height based on position in ECEF
            #alt_adj = p/cos(outLat) -  (aF0 / sqrt(1 - e2sqr * sin(outLat)**2) + 0.1)
            alt_adj = p / cos(outLat) - rN
            #alt_adj = ((1 - e2)*RN + h) * sin(outLat)
            #print alt_adj

            # x1 = (r + alt_adj)*cos(outLat)*cos(outLon)
            # y1 = (r + alt_adj)*cos(outLat)*sin(outLon)
            # z1 = (r + alt_adj)*sin(outLat)
            #
            # x2 = (RN + alt_adj) * cos(outLat) * cos(outLon)
            # y2 = (RN + alt_adj) * cos(outLat) * sin(outLon)
            # z2 = ((1-e2sqr)*RN + alt_adj) * sin(outLat)

            #print ("CameraSync Lat/Long/Atl 41.7385661690,-70.1425501650 353.285805000")
            #print ("IPAS C0+   Lat/Long/Atl 41.7385661500,-70.1425501500 353.286000000")
            #print ("Calculated Lat/Long/Atl ") + str(degrees(outLat)), str(degrees(outLon)), str(alt_adj)
            # print x1,y1,z1
            # print x2,y2,z2

            #Had this in here to calculate the difference between the Ipas CO+ image locations but not necessary here
            # ipLat = radians(41.73856615)
            # ipLon = radians(-70.14255015)
            # dLat = (outLat - ipLat)
            # dLon = (outLon - ipLon)
            # a2 = sin(dLat / 2) * sin(dLat / 2) + cos(ipLat) * cos(outLat) * sin(dLon / 2) * sin(dLon / 2)
            # c = 2 * atan2(sqrt(a2), sqrt(1 - a2))
            # d = radius * c
            # print "Distance is equal to " + str(d)

            # Convert center Lat, Long to UTM
            p = Proj(proj='utm', zone=zone, ellps='GRS80')
            cenEast, cenNorth = p(degrees(outLon),degrees(outLat))

            cenEastF = format(cenEast, '.11f')
            cenNorthF = format(cenNorth, '.11f')

            count_run_times += 1
            eoLines = ('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n') % (imageName,imageNum,cenEastF,cenNorthF,alt_adj,degrees(omega),degrees(phi),degrees(kappa),degrees(outLat),degrees(outLon))
            #print eoLines
            eoTxtFile.writelines(eoLines)

        eoTxtFile.close()
        print ("Number of times run: " + str(count_run_times))
        print ('\n\nProcessing completed in: %s' % ((time.time() - start) / 60), 'minutes\n')

    except IOError as e:  # prints the associate error
        print ("I/O Error %d: %s" % (e.errno, e.strerror))
        raise e
