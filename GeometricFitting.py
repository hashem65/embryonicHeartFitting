#!/usr/bin/env python

#> \file 
#> \author Hashem Yousefi 
#> \This code creates the stage meshes of the developing heart from the linear
#> heart tube to the deformed immature s-looped heart. The aim of this study 
#> is to provide a micro-mechanical model of the heart for the stages of the
#> development. Two different atlases have been combined to provide this data
#> in 5 different stages. We used the Kaufman histology atlas and Edinburgh 
#> mouse anatomy atlas in the early stages of heart formation. The aim is to 
#> track the recognised regions of the heart partitions in each developing 
#> stage and provide a consequence for them. It will help to understand the  
#> kinematic of the heart tissue motion, and it will be useful for further 
#> studies. The final aim is to study the growth process of forming heart.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>
#<

import sys, os
import exfile
import numpy
from numpy import linalg
import math
import random

# Intialise OpenCMISS/iron 
from opencmiss.iron import iron

# defining the output file to be written in the ExDataFile
def writeExdataFile(filename,dataPointLocations,dataErrorVector,dataErrorDistance,offset):
    "Writes data points to an exdata file"

    numberOfDimensions = dataPointLocations[1].shape[0]
    try:
        f = open(filename,"w")
        if numberOfDimensions == 1:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=2, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=3, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 2:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=3, #Derivatives=0, #Versions=1
  y.  Value index=4, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=5, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 3:
             header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
  x.  Value index=3, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=4, #Derivatives=0, #Versions=1
  y.  Value index=5, #Derivatives=0, #Versions=1
  z.  Value index=6, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=7, #Derivatives=0, #Versions=1
'''
        f.write(header)

        numberOfDataPoints = len(dataPointLocations)
        for i in range(numberOfDataPoints):
            line = " Node: " + str(offset+i+1) + '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataPointLocations[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataErrorVector[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            line = ' ' + str(dataErrorDistance[i])
            f.write(line)
            line = '\n'
            f.write(line)
        f.close()
            
    except IOError:
        print ('Could not open file: ' + filename)



#=================================================================#
#       Which stage to show from Kaufman 10, 11, 12, 13, 14       #     
#                                                                 #
#       which surface to fit  inner surface / outer surface       #
#=================================================================#
# one of the five below lines should be True 
Kaufman10 = False  
Kaufman11 = False  
Kaufman12 = False 
Kaufman13 = False
Kaufman14 = True 

KaufmanNumber = 10 
if (Kaufman10):
    KaufmanNumber = 10 
if (Kaufman11):
    KaufmanNumber = 11 
if (Kaufman12):
    KaufmanNumber = 12 
if (Kaufman13):
    KaufmanNumber = 13 
if (Kaufman14):
    KaufmanNumber = 14 

Epi = True 
#Epi = False

if (Kaufman10):
    if (Epi):
        numberOfDataPoints = 1024
    else: 
        numberOfDataPoints = 473

if (Kaufman11):
    if (Epi):
        numberOfDataPoints = 351
    else: 
        numberOfDataPoints = 92

if (Kaufman12):
    if (Epi):
        numberOfDataPoints = 551
    else: 
        numberOfDataPoints = 256

if (Kaufman13):
    if (Epi):
        numberOfDataPoints = 854
    else: 
        numberOfDataPoints = 532

if (Kaufman14):
    if (Epi):
        numberOfDataPoints = 2691
    else: 
        numberOfDataPoints = 1169


#=================================================================#
#                                                                 #
#                       Control Panel                             #
#                                                                 #
#=================================================================#

# set the number of elements and the number of nodes for the cylinder 
numberOfDimensions = 3
numberOfGaussXi = 3 
numberOfCircumfrentialElementsPerQuarter = 2
numberOfCircumfrentialElements = 4*numberOfCircumfrentialElementsPerQuarter
numberOfCircumfrentialNodes = numberOfCircumfrentialElements
numberOfLengthElements = 8
numberOfLengthNodes = numberOfLengthElements+1
numberOfWallElements = 1
numberOfWallNodes = numberOfWallElements+1
origin = [0.0,0.0,0.0]
meshOrigin = [0.0,0.0,0.0]
print "mesh resolution and parameters fixed"

# The number of data points which are digitised from the heart segments 
# fix interior nodes so that fitting only applies to surface
# If start iteration > 1, read in geometry from a previous fit iteration
# Set Sobolev smoothing parameters

tau = 0.00001
kappa = 0.00001
zeroTolerance = 0.00001
fixInterior = True
hermite = True
startpoint = 0
startIteration = 0
if startIteration > 1:
    exfileMesh = True
    exnode = exfile.Exnode("DeformedGeometry" + str(startIteration-1) + ".part0.exnode")
    exelem = exfile.Exelem("UndeformedGeometry.part0.exelem")
else:
    exfileMesh = False
print "other parameters setted up "


#=================================================================#
#                                                                 #
#        CS      ,        Region     ,      Basis                 #
#                                                                 #
#=================================================================#
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    independentFieldUserNumber,
    dataPointFieldUserNumber,
    materialFieldUserNumber,
    analyticFieldUserNumber,
    dependentDataFieldUserNumber,
    dataPointsUserNumber,
    dataProjectionUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,19)

# Get the computational nodes information
#print dir(iron),'\n\n'
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC CS
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "FittingRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# define a basis 
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = numberOfDimensions
if hermite:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
else:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()
print "CS, Region and basis setted up"



#=================================================================#
#                                                                 #
#                         M e s h                                 #
#                                                                 #
#=================================================================#

# creating the number of elements and the mesh origins ... and/or
# Start the creation of a manually generated mesh in the region
numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements

print "numberOfElements = ", numberOfElements
print "numberOfNodes = ", numberOfNodes


if (exfileMesh):
    # Read previous mesh
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, numberOfDimensions)
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(exelem.num_elements)
    # Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region, exnode.num_nodes)
    nodes.CreateFinish()
    # Define elements for the mesh
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    for elem in exelem.elements:
        elements.NodesSet(elem.number, elem.nodes)
    elements.CreateFinish()
    mesh.CreateFinish()
else:
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber,region,3)
    mesh.origin = meshOrigin
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(numberOfElements)
# Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region,numberOfNodes)
    nodes.CreateFinish()
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    elementNumber = 0
    for wallElementIdx in range(1,numberOfWallElements+1):
       for lengthElementIdx in range(1,numberOfLengthElements+1):
            for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                elementNumber = elementNumber + 1
                localNode1 = circumfrentialElementIdx + (lengthElementIdx - 1)*numberOfCircumfrentialElements + \
                    (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                if circumfrentialElementIdx == numberOfCircumfrentialElements:
                    localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                        (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                else: 
                    localNode2 = localNode1 + 1
                localNode3 = localNode1 + numberOfCircumfrentialNodes
                localNode4 = localNode2 + numberOfCircumfrentialNodes
                localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
#		print "Element Number = ",elementNumber
#         	print "Node numbers of the element", localNode1, localNode2, localNode3, localNode4, localNode5, localNode6, localNode7, localNode8 
                elements.NodesSet(elementNumber,localNodes)  
    elements.CreateFinish()
    mesh.CreateFinish() 


# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

print "mesh decomposition finished"

#=================================================================#
#                                                                 #
#            p o s i t i o n   of the   n o d e s                 #      
#                                                                 #
#=================================================================#
# the location of nodes for the mesh of the each of the stages are different  
# for stage 8dpc equivalent to Kaufman10 inner surface  
if (Kaufman10): 
    manualNodePoints = numpy.zeros((numberOfLengthNodes,8,3,2))

    manualNodePoints[0,0,:,0] = [528,318,-5]
    manualNodePoints[0,1,:,0] = [521,321,-5]
    manualNodePoints[0,2,:,0] = [510,324,-5]
    manualNodePoints[0,3,:,0] = [502,322,-5]
    manualNodePoints[0,4,:,0] = [497,317,-5]
    manualNodePoints[0,5,:,0] = [500,311,-5]
    manualNodePoints[0,6,:,0] = [508,305,-5]
    manualNodePoints[0,7,:,0] = [522,308,-5]

    manualNodePoints[1,0,:,0] = [529,314,25]
    manualNodePoints[1,1,:,0] = [519,319,25]
    manualNodePoints[1,2,:,0] = [509,322,25]
    manualNodePoints[1,3,:,0] = [497,319,25]
    manualNodePoints[1,4,:,0] = [494,311,25]
    manualNodePoints[1,5,:,0] = [501,306,25]
    manualNodePoints[1,6,:,0] = [511,306,25]
    manualNodePoints[1,7,:,0] = [523,307,25]

    manualNodePoints[2,0,:,0] = [548,300,130]
    manualNodePoints[2,1,:,0] = [541,312,127]
    manualNodePoints[2,2,:,0] = [526,330,131]
    manualNodePoints[2,3,:,0] = [512,311,128]
    manualNodePoints[2,4,:,0] = [505,296,130]
    manualNodePoints[2,5,:,0] = [520,290,130]
    manualNodePoints[2,6,:,0] = [531,287,130]
    manualNodePoints[2,7,:,0] = [544,285,130]

    manualNodePoints[3,0,:,0] = [596,270,255]
    manualNodePoints[3,1,:,0] = [573,302,245]
    manualNodePoints[3,2,:,0] = [546,332,217]
    manualNodePoints[3,3,:,0] = [498,296,242]
    manualNodePoints[3,4,:,0] = [467,262,255]
    manualNodePoints[3,5,:,0] = [491,246,255]
    manualNodePoints[3,6,:,0] = [534,253,255]
    manualNodePoints[3,7,:,0] = [564,252,255]

    manualNodePoints[4,0,:,0] = [625,230,380]
    manualNodePoints[4,1,:,0] = [591,250,380]
    manualNodePoints[4,2,:,0] = [530,223,380]
    manualNodePoints[4,3,:,0] = [488,267,380]
    manualNodePoints[4,4,:,0] = [413,265,380]
    manualNodePoints[4,5,:,0] = [455,230,380]
    manualNodePoints[4,6,:,0] = [547,194,380]
    manualNodePoints[4,7,:,0] = [587,201,380]

    manualNodePoints[5,0,:,0] = [672,250,490]
    manualNodePoints[5,1,:,0] = [610,262,490]
    manualNodePoints[5,2,:,0] = [564,255,490]
    manualNodePoints[5,3,:,0] = [524,257,490]
    manualNodePoints[5,4,:,0] = [471,296,490]
    manualNodePoints[5,5,:,0] = [495,245,490]
    manualNodePoints[5,6,:,0] = [523,213,490]
    manualNodePoints[5,7,:,0] = [590,232,490]

    manualNodePoints[6,0,:,0] = [637,249,580]
    manualNodePoints[6,1,:,0] = [627,270,580]
    manualNodePoints[6,2,:,0] = [600,294,581]
    manualNodePoints[6,3,:,0] = [584,298,580]
    manualNodePoints[6,4,:,0] = [562,302,580]
    manualNodePoints[6,5,:,0] = [557,271,580]
    manualNodePoints[6,6,:,0] = [582,260,580]
    manualNodePoints[6,7,:,0] = [599,231,580]

    manualNodePoints[7,0,:,0] = [609,307,610]
    manualNodePoints[7,1,:,0] = [595,336,610]
    manualNodePoints[7,2,:,0] = [576,318,610]
    manualNodePoints[7,3,:,0] = [570,300,610]
    manualNodePoints[7,4,:,0] = [563,282,610]
    manualNodePoints[7,5,:,0] = [574,273,610]
    manualNodePoints[7,6,:,0] = [585,272,610]
    manualNodePoints[7,7,:,0] = [604,285,610]
    
    manualNodePoints[8,0,:,0] = [609,307,635]
    manualNodePoints[8,1,:,0] = [595,336,635]
    manualNodePoints[8,2,:,0] = [576,318,635]
    manualNodePoints[8,3,:,0] = [570,300,635]
    manualNodePoints[8,4,:,0] = [563,282,635]
    manualNodePoints[8,5,:,0] = [574,273,635]
    manualNodePoints[8,6,:,0] = [585,272,635]
    manualNodePoints[8,7,:,0] = [602,284,635]
    
    # node positions of the outer surface ... 

    manualNodePoints[0,0,:,1] = [600,311,-5]
    manualNodePoints[0,1,:,1] = [578,331,-5]
    manualNodePoints[0,2,:,1] = [518,345,-5]
    manualNodePoints[0,3,:,1] = [461,340,-5]
    manualNodePoints[0,4,:,1] = [437,311,-5]
    manualNodePoints[0,5,:,1] = [455,288,-5]
    manualNodePoints[0,6,:,1] = [520,278,-5]
    manualNodePoints[0,7,:,1] = [574,290,-5]

    manualNodePoints[1,0,:,1] = [611,285,25]
    manualNodePoints[1,1,:,1] = [590,325,25]
    manualNodePoints[1,2,:,1] = [520,350,25]
    manualNodePoints[1,3,:,1] = [435,355,25]
    manualNodePoints[1,4,:,1] = [395,325,25]
    manualNodePoints[1,5,:,1] = [430,285,25]
    manualNodePoints[1,6,:,1] = [513,252,25]
    manualNodePoints[1,7,:,1] = [575,260,25]

    manualNodePoints[2,0,:,1] = [661,280,135]
    manualNodePoints[2,1,:,1] = [656,367,158]
    manualNodePoints[2,2,:,1] = [534,397,158]
    manualNodePoints[2,3,:,1] = [391,373,159]
    manualNodePoints[2,4,:,1] = [331,300,125]
    manualNodePoints[2,5,:,1] = [390,245,135]
    manualNodePoints[2,6,:,1] = [495,220,135]
    manualNodePoints[2,7,:,1] = [640,235,135]
    
    manualNodePoints[3,0,:,1] = [706,216,260]
    manualNodePoints[3,1,:,1] = [673,293,270]
    manualNodePoints[3,2,:,1] = [530,366,252]
    manualNodePoints[3,3,:,1] = [392,347,260]
    manualNodePoints[3,4,:,1] = [315,290,260]
    manualNodePoints[3,5,:,1] = [416,203,260]
    manualNodePoints[3,6,:,1] = [526,176,260]
    manualNodePoints[3,7,:,1] = [675,185,260]
    
    manualNodePoints[4,0,:,1] = [775,210,380]
    manualNodePoints[4,1,:,1] = [750,280,380]
    manualNodePoints[4,2,:,1] = [595,310,380]
    manualNodePoints[4,3,:,1] = [403,325,380]
    manualNodePoints[4,4,:,1] = [313,270,380]
    manualNodePoints[4,5,:,1] = [350,190,380]
    manualNodePoints[4,6,:,1] = [500,140,380]
    manualNodePoints[4,7,:,1] = [700,160,380]
    
    manualNodePoints[5,0,:,1] = [805,223,490]
    manualNodePoints[5,1,:,1] = [755,323,499]
    manualNodePoints[5,2,:,1] = [590,320,490]
    manualNodePoints[5,3,:,1] = [441,365,492]
    manualNodePoints[5,4,:,1] = [305,305,490]
    manualNodePoints[5,5,:,1] = [324,205,490]
    manualNodePoints[5,6,:,1] = [465,180,490]
    manualNodePoints[5,7,:,1] = [675,150,490]
    
    manualNodePoints[6,0,:,1] = [845,240,590]
    manualNodePoints[6,1,:,1] = [780,312,590]
    manualNodePoints[6,2,:,1] = [603,363,590]
    manualNodePoints[6,3,:,1] = [425,340,590]
    manualNodePoints[6,4,:,1] = [330,290,590]
    manualNodePoints[6,5,:,1] = [385,205,590]
    manualNodePoints[6,6,:,1] = [524,155,590]
    manualNodePoints[6,7,:,1] = [717,164,590]
    
    manualNodePoints[7,0,:,1] = [840,267,615]
    manualNodePoints[7,1,:,1] = [785,337,615]
    manualNodePoints[7,2,:,1] = [603,363,615]
    manualNodePoints[7,3,:,1] = [427,353,615]
    manualNodePoints[7,4,:,1] = [342,312,615]
    manualNodePoints[7,5,:,1] = [385,205,615]
    manualNodePoints[7,6,:,1] = [524,155,615]
    manualNodePoints[7,7,:,1] = [717,164,615]
    
    manualNodePoints[8,0,:,1] = [827,248,640]
    manualNodePoints[8,1,:,1] = [785,348,640]
    manualNodePoints[8,2,:,1] = [603,363,640]
    manualNodePoints[8,3,:,1] = [431,365,640]
    manualNodePoints[8,4,:,1] = [307,333,640]
    manualNodePoints[8,5,:,1] = [385,205,640]
    manualNodePoints[8,6,:,1] = [524,155,640]
    manualNodePoints[8,7,:,1] = [717,164,640]

# for stage Kaufman11 inner surface 
if (Kaufman11): 
    manualNodePoints = numpy.zeros((numberOfLengthNodes,8,3,2))

    manualNodePoints[0,0,:,0] = [748,235,-38]
    manualNodePoints[0,1,:,0] = [700,260,-40]
    manualNodePoints[0,2,:,0] = [650,270,-31]
    manualNodePoints[0,3,:,0] = [617,236,-35]
    manualNodePoints[0,4,:,0] = [605,205,-25]
    manualNodePoints[0,5,:,0] = [635,175,-18]
    manualNodePoints[0,6,:,0] = [695,165,-37]
    manualNodePoints[0,7,:,0] = [735,185,-32]
    
    manualNodePoints[1,0,:,0] = [718,235,60]
    manualNodePoints[1,1,:,0] = [670,260,60]
    manualNodePoints[1,2,:,0] = [620,270,60]
    manualNodePoints[1,3,:,0] = [587,236,60]
    manualNodePoints[1,4,:,0] = [575,205,60]
    manualNodePoints[1,5,:,0] = [605,175,60]
    manualNodePoints[1,6,:,0] = [665,165,60]
    manualNodePoints[1,7,:,0] = [705,185,60]
    
    manualNodePoints[2,0,:,0] = [714,220,203]
    manualNodePoints[2,1,:,0] = [701,279,197]
    manualNodePoints[2,2,:,0] = [603,295,201]
    manualNodePoints[2,3,:,0] = [531,307,193]
    manualNodePoints[2,4,:,0] = [486,266,207]
    manualNodePoints[2,5,:,0] = [520,232,205]
    manualNodePoints[2,6,:,0] = [579,201,206]
    manualNodePoints[2,7,:,0] = [654,191,210]
    
    manualNodePoints[3,0,:,0] = [679,201,380]
    manualNodePoints[3,1,:,0] = [638,227,380]
    manualNodePoints[3,2,:,0] = [553,256,380]
    manualNodePoints[3,3,:,0] = [473,252,378]
    manualNodePoints[3,4,:,0] = [400,207,380]
    manualNodePoints[3,5,:,0] = [456,162,382]
    manualNodePoints[3,6,:,0] = [565,160,380]
    manualNodePoints[3,7,:,0] = [637,160,384]
    
    manualNodePoints[4,0,:,0] = [738,303,540]
    manualNodePoints[4,7,:,0] = [695,207,540]
    manualNodePoints[4,6,:,0] = [573,126,540]
    manualNodePoints[4,5,:,0] = [478,118,540]
    manualNodePoints[4,4,:,0] = [452,180,540]
    manualNodePoints[4,3,:,0] = [537,247,540]
    manualNodePoints[4,2,:,0] = [592,266,540]
    manualNodePoints[4,1,:,0] = [670,304,540]
    
    manualNodePoints[5,0,:,0] = [600,303,720]
    manualNodePoints[5,1,:,0] = [531,354,720]
    manualNodePoints[5,2,:,0] = [453,319,720]
    manualNodePoints[5,3,:,0] = [382,235,720]
    manualNodePoints[5,4,:,0] = [422,174,720]
    manualNodePoints[5,5,:,0] = [496,159,720]
    manualNodePoints[5,6,:,0] = [562,216,720]
    manualNodePoints[5,7,:,0] = [615,233,715]
    
    manualNodePoints[6,0,:,0] = [631,338,897]
    manualNodePoints[6,1,:,0] = [611,412,890]
    manualNodePoints[6,2,:,0] = [498,398,890]
    manualNodePoints[6,3,:,0] = [381,370,888]
    manualNodePoints[6,4,:,0] = [318,305,890]
    manualNodePoints[6,5,:,0] = [422,251,888]
    manualNodePoints[6,6,:,0] = [472,262,890]
    manualNodePoints[6,7,:,0] = [591,312,890]
    
    manualNodePoints[7,0,:,0] = [656,330,977]
    manualNodePoints[7,1,:,0] = [628,416,982]
    manualNodePoints[7,2,:,0] = [503,398,960]
    manualNodePoints[7,3,:,0] = [393,367,968]
    manualNodePoints[7,4,:,0] = [338,305,970]
    manualNodePoints[7,5,:,0] = [380,257,965]
    manualNodePoints[7,6,:,0] = [455,286,963]
    manualNodePoints[7,7,:,0] = [572,326,970]
    
    manualNodePoints[8,0,:,0] = [652,329,1050]
    manualNodePoints[8,1,:,0] = [605,416,1037]
    manualNodePoints[8,2,:,0] = [510,405,1045]
    manualNodePoints[8,3,:,0] = [413,361,1044]
    manualNodePoints[8,4,:,0] = [358,305,1050]
    manualNodePoints[8,5,:,0] = [392,284,1059]
    manualNodePoints[8,6,:,0] = [480,321,1053]
    manualNodePoints[8,7,:,0] = [575,335,1047]

    # node positions of the outer surface ... 
    
    manualNodePoints[0,0,:,1] = [855,220,-70]
    manualNodePoints[0,1,:,1] = [881,346,-50]
    manualNodePoints[0,2,:,1] = [705,393,-54]
    manualNodePoints[0,3,:,1] = [528,387,-40]
    manualNodePoints[0,4,:,1] = [442,315,-30]
    manualNodePoints[0,5,:,1] = [452,123,-30]
    manualNodePoints[0,6,:,1] = [614,58,-40]
    manualNodePoints[0,7,:,1] = [797,76,-50]
    
    manualNodePoints[1,0,:,1] = [925,248,60]
    manualNodePoints[1,1,:,1] = [899,402,63]
    manualNodePoints[1,2,:,1] = [733,433,60]
    manualNodePoints[1,3,:,1] = [545,412,75]
    manualNodePoints[1,4,:,1] = [371,278,63]
    manualNodePoints[1,5,:,1] = [452,170,79]
    manualNodePoints[1,6,:,1] = [604,123,80]
    manualNodePoints[1,7,:,1] = [850,107,65]

    manualNodePoints[2,0,:,1] = [965,313,216]
    manualNodePoints[2,1,:,1] = [904,413,211]
    manualNodePoints[2,2,:,1] = [696,434,220]
    manualNodePoints[2,3,:,1] = [436,407,179]
    manualNodePoints[2,4,:,1] = [275,186,180]
    manualNodePoints[2,5,:,1] = [422,79,182]
    manualNodePoints[2,6,:,1] = [615,60,213]
    manualNodePoints[2,7,:,1] = [884,150,216]
    
    manualNodePoints[3,0,:,1] = [966,352,397]
    manualNodePoints[3,1,:,1] = [876,410,375]
    manualNodePoints[3,2,:,1] = [705,389,369]
    manualNodePoints[3,3,:,1] = [407,347,356]
    manualNodePoints[3,4,:,1] = [264,173,360]
    manualNodePoints[3,5,:,1] = [430,29,380]
    manualNodePoints[3,6,:,1] = [656,55,390]
    manualNodePoints[3,7,:,1] = [915,158,382]

    manualNodePoints[4,0,:,1] = [1047,336,541]
    manualNodePoints[4,1,:,1] = [911,408,525]
    manualNodePoints[4,2,:,1] = [679,397,528]
    manualNodePoints[4,3,:,1] = [395,377,526]
    manualNodePoints[4,4,:,1] = [275,236,535]
    manualNodePoints[4,5,:,1] = [460,79,550]
    manualNodePoints[4,6,:,1] = [665,63,537]
    manualNodePoints[4,7,:,1] = [919,150,549]

    manualNodePoints[5,0,:,1] = [1013,407,717]
    manualNodePoints[5,1,:,1] = [897,487,693]
    manualNodePoints[5,2,:,1] = [669,429,715]
    manualNodePoints[5,3,:,1] = [364,392,718]
    manualNodePoints[5,4,:,1] = [178,194,741]
    manualNodePoints[5,5,:,1] = [406,35,751]
    manualNodePoints[5,6,:,1] = [624,41,738]
    manualNodePoints[5,7,:,1] = [845,158,718]

    manualNodePoints[6,0,:,1] = [970,443,861]
    manualNodePoints[6,1,:,1] = [891,567,837]
    manualNodePoints[6,2,:,1] = [630,465,892]
    manualNodePoints[6,3,:,1] = [313,436,865]
    manualNodePoints[6,4,:,1] = [165,253,894]
    manualNodePoints[6,5,:,1] = [295,110,876]
    manualNodePoints[6,6,:,1] = [570,89,851]
    manualNodePoints[6,7,:,1] = [776,145,871]

    manualNodePoints[7,0,:,1] = [921,501,952]
    manualNodePoints[7,1,:,1] = [811,551,938]
    manualNodePoints[7,2,:,1] = [535,460,978]
    manualNodePoints[7,3,:,1] = [246,507,961]
    manualNodePoints[7,4,:,1] = [72,437,981]
    manualNodePoints[7,5,:,1] = [261,194,976]
    manualNodePoints[7,6,:,1] = [546,166,979]
    manualNodePoints[7,7,:,1] = [772,281,965]

    manualNodePoints[8,0,:,1] = [897,789,1048]
    manualNodePoints[8,1,:,1] = [706,588,1050]
    manualNodePoints[8,2,:,1] = [461,464,1048]
    manualNodePoints[8,3,:,1] = [226,561,1040]
    manualNodePoints[8,4,:,1] = [75,626,1048]
    manualNodePoints[8,5,:,1] = [176,230,1048]
    manualNodePoints[8,6,:,1] = [505,277,1048]
    manualNodePoints[8,7,:,1] = [828,339,1044]

# for stage Kaufman12 inner surface 
if (Kaufman12): 
    manualNodePoints = numpy.zeros((numberOfLengthNodes,8,3,2))
    
    manualNodePoints[0,0,:,0] = [646,322,-38]
    manualNodePoints[0,1,:,0] = [630,373,-40]
    manualNodePoints[0,2,:,0] = [581,391,-31]
    manualNodePoints[0,3,:,0] = [536,381,-35]
    manualNodePoints[0,4,:,0] = [491,347,-25]
    manualNodePoints[0,5,:,0] = [505,307,-18]
    manualNodePoints[0,6,:,0] = [561,282,-37]
    manualNodePoints[0,7,:,0] = [641,275,-32]

    manualNodePoints[1,0,:,0] = [841,308,57]
    manualNodePoints[1,1,:,0] = [812,353,60]
    manualNodePoints[1,2,:,0] = [783,400,72]
    manualNodePoints[1,3,:,0] = [747,366,70]
    manualNodePoints[1,4,:,0] = [742,332,82]
    manualNodePoints[1,5,:,0] = [696,255,65]
    manualNodePoints[1,6,:,0] = [790,250,70]
    manualNodePoints[1,7,:,0] = [834,260,58]
    
    manualNodePoints[2,0,:,0] = [1106,259,253]
    manualNodePoints[2,1,:,0] = [1095,327,250]
    manualNodePoints[2,2,:,0] = [1019,362,260]
    manualNodePoints[2,3,:,0] = [949,340,255]
    manualNodePoints[2,4,:,0] = [890,300,250]
    manualNodePoints[2,5,:,0] = [920,272,243]
    manualNodePoints[2,6,:,0] = [974,262,238]
    manualNodePoints[2,7,:,0] = [1038,256,239]
    
    manualNodePoints[3,0,:,0] = [1050,300,455]
    manualNodePoints[3,1,:,0] = [1018,357,460]
    manualNodePoints[3,2,:,0] = [933,386,449]
    manualNodePoints[3,3,:,0] = [851,370,440]
    manualNodePoints[3,4,:,0] = [780,337,431]
    manualNodePoints[3,5,:,0] = [834,315,431]
    manualNodePoints[3,6,:,0] = [945,280,441]
    manualNodePoints[3,7,:,0] = [1035,275,450]
    
    manualNodePoints[4,0,:,0] = [1138,303,728]
    manualNodePoints[4,7,:,0] = [1045,207,685]
    manualNodePoints[4,6,:,0] = [873,126,621]
    manualNodePoints[4,5,:,0] = [728,118,601]
    manualNodePoints[4,4,:,0] = [656,225,561]
    manualNodePoints[4,3,:,0] = [707,247,600]
    manualNodePoints[4,2,:,0] = [842,266,636]
    manualNodePoints[4,1,:,0] = [970,304,660]
    
    manualNodePoints[5,0,:,0] = [681,372,812]
    manualNodePoints[5,1,:,0] = [570,370,800]
    manualNodePoints[5,2,:,0] = [474,328,774]
    manualNodePoints[5,3,:,0] = [382,235,736]
    manualNodePoints[5,4,:,0] = [422,174,727]
    manualNodePoints[5,5,:,0] = [558,167,740]
    manualNodePoints[5,6,:,0] = [640,236,757]
    manualNodePoints[5,7,:,0] = [732,265,783]
    
    manualNodePoints[6,0,:,0] = [722,329,927]
    manualNodePoints[6,1,:,0] = [674,419,917]
    manualNodePoints[6,2,:,0] = [563,412,925]
    manualNodePoints[6,3,:,0] = [371,397,945]
    manualNodePoints[6,4,:,0] = [318,305,911]
    manualNodePoints[6,5,:,0] = [379,233,915]
    manualNodePoints[6,6,:,0] = [472,262,930]
    manualNodePoints[6,7,:,0] = [591,312,930]
    
    manualNodePoints[7,0,:,0] = [800,342,1090]
    manualNodePoints[7,1,:,0] = [680,370,1080]
    manualNodePoints[7,2,:,0] = [490,360,1080]
    manualNodePoints[7,3,:,0] = [315,355,1080]
    manualNodePoints[7,4,:,0] = [206,503,1101]
    manualNodePoints[7,5,:,0] = [175,375,1080]
    manualNodePoints[7,6,:,0] = [461,276,1082]
    manualNodePoints[7,7,:,0] = [645,318,1084]
    
    manualNodePoints[8,0,:,0] = [703,368,1211]
    manualNodePoints[8,1,:,0] = [482,403,1211]
    manualNodePoints[8,2,:,0] = [315,383,1213]
    manualNodePoints[8,3,:,0] = [228,544,1211]
    manualNodePoints[8,4,:,0] = [130,580,1163]
    manualNodePoints[8,5,:,0] = [115,465,1199]
    manualNodePoints[8,6,:,0] = [248,316,1220]
    manualNodePoints[8,7,:,0] = [487,322,1208]

    # node positions of the outer surface ... 
    
    manualNodePoints[0,0,:,1] = [844,231,-32]
    manualNodePoints[0,7,:,1] = [792,88,-51]
    manualNodePoints[0,6,:,1] = [630,39,-48]
    manualNodePoints[0,5,:,1] = [450,138,-41]
    manualNodePoints[0,4,:,1] = [325,365,-30]
    manualNodePoints[0,3,:,1] = [501,481,-20]
    manualNodePoints[0,2,:,1] = [700,490,-40]
    manualNodePoints[0,1,:,1] = [833,350,-50]

    manualNodePoints[1,0,:,1] = [944,285,61]
    manualNodePoints[1,1,:,1] = [870,440,65]
    manualNodePoints[1,2,:,1] = [757,542,125]
    manualNodePoints[1,3,:,1] = [500,517,153]
    manualNodePoints[1,4,:,1] = [384,386,120]
    manualNodePoints[1,5,:,1] = [452,173,85]
    manualNodePoints[1,6,:,1] = [706,93,73]
    manualNodePoints[1,7,:,1] = [896,114,45]

    manualNodePoints[2,0,:,1] = [1323,317,206]
    manualNodePoints[2,1,:,1] = [1140,430,230]
    manualNodePoints[2,2,:,1] = [832,489,267]
    manualNodePoints[2,3,:,1] = [493,400,263]
    manualNodePoints[2,4,:,1] = [265,183,260]
    manualNodePoints[2,5,:,1] = [438,51,246]
    manualNodePoints[2,6,:,1] = [762,41,237]
    manualNodePoints[2,7,:,1] = [1083,122,228]

    manualNodePoints[3,0,:,1] = [1290,360,450]
    manualNodePoints[3,7,:,1] = [1153,181,459]
    manualNodePoints[3,6,:,1] = [804,59,465]
    manualNodePoints[3,5,:,1] = [435,42,452]
    manualNodePoints[3,4,:,1] = [265,230,450]
    manualNodePoints[3,3,:,1] = [466,419,455]
    manualNodePoints[3,2,:,1] = [810,445,468]
    manualNodePoints[3,1,:,1] = [1120,450,450]

    manualNodePoints[4,0,:,1] = [1456,315,702]
    manualNodePoints[4,1,:,1] = [1078,427,697]
    manualNodePoints[4,2,:,1] = [740,430,710]
    manualNodePoints[4,3,:,1] = [400,400,720]
    manualNodePoints[4,4,:,1] = [204,334,692]
    manualNodePoints[4,5,:,1] = [430,110,701]
    manualNodePoints[4,6,:,1] = [814,17,709]
    manualNodePoints[4,7,:,1] = [1155,151,676]

    manualNodePoints[5,0,:,1] = [1185,203,877]
    manualNodePoints[5,1,:,1] = [1047,402,833]
    manualNodePoints[5,2,:,1] = [691,432,834]
    manualNodePoints[5,3,:,1] = [382,382,824]
    manualNodePoints[5,4,:,1] = [227,302,831]
    manualNodePoints[5,5,:,1] = [397,90,836]
    manualNodePoints[5,6,:,1] = [797,-11,841]
    manualNodePoints[5,7,:,1] = [1050,67,831]

    manualNodePoints[6,0,:,1] = [1039,249,969]
    manualNodePoints[6,1,:,1] = [955,480,950]
    manualNodePoints[6,2,:,1] = [615,442,952]
    manualNodePoints[6,3,:,1] = [300,400,950]
    manualNodePoints[6,4,:,1] = [220,303,949]
    manualNodePoints[6,5,:,1] = [342,145,950]
    manualNodePoints[6,6,:,1] = [733,29,951]
    manualNodePoints[6,7,:,1] = [930,74,970]

    manualNodePoints[7,7,:,1] = [890,110,1100]
    manualNodePoints[7,6,:,1] = [663,99,1110]
    manualNodePoints[7,5,:,1] = [254,219,1080]
    manualNodePoints[7,4,:,1] = [95,450,1070]
    manualNodePoints[7,3,:,1] = [235,575,1070]
    manualNodePoints[7,2,:,1] = [501,399,1071]
    manualNodePoints[7,1,:,1] = [860,480,1070]
    manualNodePoints[7,0,:,1] = [1029,277,1081]

    manualNodePoints[8,0,:,1] = [1123,421,1224]
    manualNodePoints[8,1,:,1] = [795,565,1225]
    manualNodePoints[8,2,:,1] = [427,465,1225]
    manualNodePoints[8,3,:,1] = [200,638,1225]
    manualNodePoints[8,4,:,1] = [-83,703,1206]
    manualNodePoints[8,5,:,1] = [181,295,1242]
    manualNodePoints[8,6,:,1] = [626,134,1235]
    manualNodePoints[8,7,:,1] = [860,199,1237]

# for stage Kaufman13 inner surface 
if (Kaufman13): 
    manualNodePoints = numpy.zeros((numberOfLengthNodes,8,3,2))
    manualNodePoints[0,0,:,0] = [1020,510,150]
    manualNodePoints[0,1,:,0] = [920,510,150]
    manualNodePoints[0,2,:,0] = [810,505,150]
    manualNodePoints[0,3,:,0] = [690,515,150]
    manualNodePoints[0,4,:,0] = [610,523,150]
    manualNodePoints[0,5,:,0] = [700,450,150]
    manualNodePoints[0,6,:,0] = [800,410,150]
    manualNodePoints[0,7,:,0] = [930,435,150]

    manualNodePoints[1,0,:,0] = [1280,250,440]
    manualNodePoints[1,1,:,0] = [1250,315,440]
    manualNodePoints[1,2,:,0] = [1170,315,440]
    manualNodePoints[1,3,:,0] = [1080,300,440]
    manualNodePoints[1,4,:,0] = [1020,285,440]
    manualNodePoints[1,5,:,0] = [1090,240,440]
    manualNodePoints[1,6,:,0] = [1180,210,440]
    manualNodePoints[1,7,:,0] = [1275,210,440]
    
    manualNodePoints[2,0,:,0] = [1360,265,770]
    manualNodePoints[2,1,:,0] = [1320,320,755]
    manualNodePoints[2,2,:,0] = [1210,365,740]
    manualNodePoints[2,3,:,0] = [1112,380,725]
    manualNodePoints[2,4,:,0] = [1072,325,710]
    manualNodePoints[2,5,:,0] = [1098,254,725]
    manualNodePoints[2,6,:,0] = [1200,270,740]
    manualNodePoints[2,7,:,0] = [1280,235,755]

    manualNodePoints[3,0,:,0] = [1260,335,1040]
    manualNodePoints[3,1,:,0] = [1220,376,995]
    manualNodePoints[3,2,:,0] = [1177,397,950]
    manualNodePoints[3,3,:,0] = [1142,410,897]
    manualNodePoints[3,4,:,0] = [1076,319,852]
    manualNodePoints[3,5,:,0] = [1090,208,897]
    manualNodePoints[3,6,:,0] = [1155,237,943]
    manualNodePoints[3,7,:,0] = [1213,273,996]
    
    manualNodePoints[4,0,:,0] = [1123,270,1275]
    manualNodePoints[4,1,:,0] = [1040,378,1053]
    manualNodePoints[4,2,:,0] = [1000,350,960]
    manualNodePoints[4,3,:,0] = [964,324,875]
    manualNodePoints[4,4,:,0] = [977,260,741]
    manualNodePoints[4,5,:,0] = [954,175,875]
    manualNodePoints[4,6,:,0] = [995,200,960]
    manualNodePoints[4,7,:,0] = [1035,238,1053]
    
    manualNodePoints[5,0,:,0] = [860,320,1053]
    manualNodePoints[5,1,:,0] = [850,394,950]
    manualNodePoints[5,2,:,0] = [838,370,849]
    manualNodePoints[5,3,:,0] = [766,390,716]
    manualNodePoints[5,4,:,0] = [680,240,585]
    manualNodePoints[5,5,:,0] = [680,240,710]
    manualNodePoints[5,6,:,0] = [770,307,830]
    manualNodePoints[5,7,:,0] = [770,307,950]
    
    manualNodePoints[6,0,:,0] = [700,400,937]
    manualNodePoints[6,1,:,0] = [667,480,934]
    manualNodePoints[6,2,:,0] = [820,610,900]
    manualNodePoints[6,3,:,0] = [650,650,860]
    manualNodePoints[6,4,:,0] = [570,620,860]
    manualNodePoints[6,5,:,0] = [530,500,877]
    manualNodePoints[6,6,:,0] = [560,400,877]
    manualNodePoints[6,7,:,0] = [620,340,917]
    
    manualNodePoints[7,0,:,0] = [800,600,1140]
    manualNodePoints[7,1,:,0] = [970,625,1170]
    manualNodePoints[7,2,:,0] = [830,680,1140]
    manualNodePoints[7,3,:,0] = [650,680,1110]
    manualNodePoints[7,4,:,0] = [450,770,1080]
    manualNodePoints[7,5,:,0] = [330,715,1050]
    manualNodePoints[7,6,:,0] = [500,600,1080]
    manualNodePoints[7,7,:,0] = [650,550,1110]
    
    manualNodePoints[8,0,:,0] = [1080,480,1390]
    manualNodePoints[8,1,:,0] = [1270,600,1390]
    manualNodePoints[8,2,:,0] = [1070,620,1390]
    manualNodePoints[8,3,:,0] = [820,590,1390]
    manualNodePoints[8,4,:,0] = [600,625,1390]
    manualNodePoints[8,5,:,0] = [410,650,1390]
    manualNodePoints[8,6,:,0] = [550,490,1390]
    manualNodePoints[8,7,:,0] = [810,460,1390]

    # node positions of the outer surface ... 
    
    manualNodePoints[0,0,:,1] = [1187,507,150]
    manualNodePoints[0,1,:,1] = [1050,663,150]
    manualNodePoints[0,2,:,1] = [821,584,150]
    manualNodePoints[0,3,:,1] = [580,634,150]
    manualNodePoints[0,4,:,1] = [487,520,150]
    manualNodePoints[0,5,:,1] = [626,331,150]
    manualNodePoints[0,6,:,1] = [880,225,150]
    manualNodePoints[0,7,:,1] = [1206,238,150]
    
    manualNodePoints[1,0,:,1] = [1571,260,453]
    manualNodePoints[1,1,:,1] = [1333,394,440]
    manualNodePoints[1,2,:,1] = [1104,483,440]
    manualNodePoints[1,3,:,1] = [924,420,440]
    manualNodePoints[1,4,:,1] = [697,334,440]
    manualNodePoints[1,5,:,1] = [903,131,440]
    manualNodePoints[1,6,:,1] = [1159,78,440]
    manualNodePoints[1,7,:,1] = [1385,94,440]

    manualNodePoints[2,0,:,1] = [1600,251,771]
    manualNodePoints[2,1,:,1] = [1433,408,744]
    manualNodePoints[2,2,:,1] = [1210,465,740]
    manualNodePoints[2,3,:,1] = [1045,489,625]
    manualNodePoints[2,4,:,1] = [900,350,572]
    manualNodePoints[2,5,:,1] = [1021,206,606]
    manualNodePoints[2,6,:,1] = [1173,114,738]
    manualNodePoints[2,7,:,1] = [1323,135,766]

    manualNodePoints[3,0,:,1] = [1419,326,1160]
    manualNodePoints[3,1,:,1] = [1285,440,1010]
    manualNodePoints[3,2,:,1] = [1137,522,876]
    manualNodePoints[3,3,:,1] = [964,499,621]
    manualNodePoints[3,4,:,1] = [853,309,593]
    manualNodePoints[3,5,:,1] = [976,100,579]
    manualNodePoints[3,6,:,1] = [1140,80,921]
    manualNodePoints[3,7,:,1] = [1314,157,1046]

    manualNodePoints[4,0,:,1] = [963,218,1541]
    manualNodePoints[4,1,:,1] = [1074,470,1081]
    manualNodePoints[4,2,:,1] = [1011,531,886]
    manualNodePoints[4,3,:,1] = [865,503,555]
    manualNodePoints[4,4,:,1] = [629,199,471]
    manualNodePoints[4,5,:,1] = [781,91,711]
    manualNodePoints[4,6,:,1] = [992,96,974]
    manualNodePoints[4,7,:,1] = [1002,111,1271]

    manualNodePoints[5,0,:,1] = [815,367,1315]
    manualNodePoints[5,1,:,1] = [927,556,1040]
    manualNodePoints[5,2,:,1] = [1012,577,954]
    manualNodePoints[5,3,:,1] = [721,569,576]
    manualNodePoints[5,4,:,1] = [359,380,396]
    manualNodePoints[5,5,:,1] = [426,198,717]
    manualNodePoints[5,6,:,1] = [659,140,907]
    manualNodePoints[5,7,:,1] = [722,183,1138]

    manualNodePoints[6,0,:,1] = [734,452,1085]
    manualNodePoints[6,1,:,1] = [730,527,1056]
    manualNodePoints[6,2,:,1] = [1021,638,954]
    manualNodePoints[6,3,:,1] = [763,650,867]
    manualNodePoints[6,4,:,1] = [461,657,758]
    manualNodePoints[6,5,:,1] = [259,492,739]
    manualNodePoints[6,6,:,1] = [400,339,908]
    manualNodePoints[6,7,:,1] = [598,310,1094]

    manualNodePoints[7,0,:,1] = [800,500,1145]
    manualNodePoints[7,1,:,1] = [1137,493,1166]
    manualNodePoints[7,2,:,1] = [1244,665,1117]
    manualNodePoints[7,3,:,1] = [857,705,1075]
    manualNodePoints[7,4,:,1] = [423,832,1015]
    manualNodePoints[7,5,:,1] = [230,715,1080]
    manualNodePoints[7,6,:,1] = [383,448,1109]
    manualNodePoints[7,7,:,1] = [650,450,1150]

    manualNodePoints[8,0,:,1] = [1114,446,1390]
    manualNodePoints[8,1,:,1] = [1517,545,1390]
    manualNodePoints[8,2,:,1] = [1212,693,1397]
    manualNodePoints[8,3,:,1] = [823,652,1390]
    manualNodePoints[8,4,:,1] = [477,739,1390]
    manualNodePoints[8,5,:,1] = [186,765,1390]
    manualNodePoints[8,6,:,1] = [448,493,1390]
    manualNodePoints[8,7,:,1] = [807,440,1390]

# for stage 8.5 dpc equivalent to Kaufman14 inner surface 
if (Kaufman14): 
    manualNodePoints = numpy.zeros((numberOfLengthNodes,8,3,2))
    manualNodePoints[0,0,:,0] = [1190,685,-35]
    manualNodePoints[0,1,:,0] = [1090,730,-35]
    manualNodePoints[0,2,:,0] = [920,745,-35]
    manualNodePoints[0,3,:,0] = [795,715,-35]
    manualNodePoints[0,4,:,0] = [745,625,-35]
    manualNodePoints[0,5,:,0] = [700,550,-35]
    manualNodePoints[0,6,:,0] = [840,535,-35]
    manualNodePoints[0,7,:,0] = [1010,585,-35]
    
    manualNodePoints[1,0,:,0] = [1615,725,375]
    manualNodePoints[1,1,:,0] = [1560,775,375]
    manualNodePoints[1,2,:,0] = [1430,820,375]
    manualNodePoints[1,3,:,0] = [1310,788,375]
    manualNodePoints[1,4,:,0] = [1270,735,375]
    manualNodePoints[1,5,:,0] = [1375,675,375]
    manualNodePoints[1,6,:,0] = [1475,650,375]
    manualNodePoints[1,7,:,0] = [1590,680,375]
    
    manualNodePoints[2,7,:,0] = [1466,558,710]
    manualNodePoints[2,6,:,0] = [1401,539,681]
    manualNodePoints[2,5,:,0] = [1342,538,656]
    manualNodePoints[2,4,:,0] = [1281,563,630]
    manualNodePoints[2,3,:,0] = [1309,625,661]
    manualNodePoints[2,2,:,0] = [1373,673,684]
    manualNodePoints[2,1,:,0] = [1438,708,702]
    manualNodePoints[2,0,:,0] = [1524,671,736]
    
    manualNodePoints[3,7,:,0] = [1419,461,1039]
    manualNodePoints[3,6,:,0] = [1369,375,902]
    manualNodePoints[3,5,:,0] = [1322,409,804]
    manualNodePoints[3,4,:,0] = [1261,470,722]
    manualNodePoints[3,3,:,0] = [1207,563,744]
    manualNodePoints[3,2,:,0] = [1197,579,853]
    manualNodePoints[3,1,:,0] = [1245,558,931]
    manualNodePoints[3,0,:,0] = [1312,557,1040]
    
    manualNodePoints[4,7,:,0] = [1113,306,1051]
    manualNodePoints[4,6,:,0] = [1217,354,925]
    manualNodePoints[4,5,:,0] = [1181,380,825]
    manualNodePoints[4,4,:,0] = [1143,415,760]
    manualNodePoints[4,3,:,0] = [1074,557,749]
    manualNodePoints[4,2,:,0] = [1024,573,829]
    manualNodePoints[4,1,:,0] = [1043,514,950]
    manualNodePoints[4,0,:,0] = [1025,421,1013]
    
    manualNodePoints[5,6,:,0] = [780,255,672]
    manualNodePoints[5,5,:,0] = [902,245,500]
    manualNodePoints[5,4,:,0] = [917,290,400]
    manualNodePoints[5,3,:,0] = [952,438,560]
    manualNodePoints[5,2,:,0] = [865,500,700]
    manualNodePoints[5,1,:,0] = [850,570,784]
    manualNodePoints[5,0,:,0] = [710,435,900]
    manualNodePoints[5,7,:,0] = [766,223,840]
    
    manualNodePoints[6,0,:,0] = [496,526,555]
    manualNodePoints[6,1,:,0] = [595,605,600]
    manualNodePoints[6,2,:,0] = [700,632,575]
    manualNodePoints[6,3,:,0] = [588,592,512]
    manualNodePoints[6,4,:,0] = [510,544,470]
    manualNodePoints[6,5,:,0] = [357,492,435]
    manualNodePoints[6,6,:,0] = [256,412,431]
    manualNodePoints[6,7,:,0] = [385,463,517]
    
    manualNodePoints[7,1,:,0] = [570,850,855]
    manualNodePoints[7,2,:,0] = [515,900,855]
    manualNodePoints[7,3,:,0] = [350,885,855]
    manualNodePoints[7,4,:,0] = [160,800,855]
    manualNodePoints[7,5,:,0] = [48.0,700,855]
    manualNodePoints[7,6,:,0] = [110,670,855]
    manualNodePoints[7,7,:,0] = [270,700,855]
    manualNodePoints[7,0,:,0] = [450,775,855]
    
    manualNodePoints[8,1,:,0] = [1030,560,1210]
    manualNodePoints[8,2,:,0] = [770,610,1210]
    manualNodePoints[8,3,:,0] = [600,670,1210]
    manualNodePoints[8,4,:,0] = [230,750,1210]
    manualNodePoints[8,5,:,0] = [170,600,1210]
    manualNodePoints[8,6,:,0] = [290,450,1210]
    manualNodePoints[8,7,:,0] = [658,438,1210]
    manualNodePoints[8,0,:,0] = [900,475,1210]

    # node positions of the outer surface ... 

    manualNodePoints[0,7,:,1] = [1180,310,-35]
    manualNodePoints[0,6,:,1] = [865,365,-35]
    manualNodePoints[0,5,:,1] = [660,535,-35]
    manualNodePoints[0,4,:,1] = [650,650,-35]
    manualNodePoints[0,3,:,1] = [740,760,-35]
    manualNodePoints[0,2,:,1] = [860,810,-35]
    manualNodePoints[0,1,:,1] = [1075,810,-35]
    manualNodePoints[0,0,:,1] = [1300,700,-35]
    
    manualNodePoints[1,7,:,1] = [1754,513,381]
    manualNodePoints[1,6,:,1] = [1403,476,388]
    manualNodePoints[1,5,:,1] = [1100,560,375]
    manualNodePoints[1,4,:,1] = [813,701,368]
    manualNodePoints[1,3,:,1] = [1000,900,375]
    manualNodePoints[1,2,:,1] = [1365,880,375]
    manualNodePoints[1,1,:,1] = [1670,863,375]
    manualNodePoints[1,0,:,1] = [1880,710,375]
    
    manualNodePoints[2,7,:,1] = [1720,485,710]
    manualNodePoints[2,6,:,1] = [1572,407,640]
    manualNodePoints[2,5,:,1] = [1357,540,571]
    manualNodePoints[2,4,:,1] = [1160,568,561]
    manualNodePoints[2,3,:,1] = [1079,697,604]
    manualNodePoints[2,2,:,1] = [1332,787,718]
    manualNodePoints[2,1,:,1] = [1541,846,710]
    manualNodePoints[2,0,:,1] = [1825,616,862]
    
    manualNodePoints[3,7,:,1] = [1582,282,1071]
    manualNodePoints[3,6,:,1] = [1413,211,772]
    manualNodePoints[3,5,:,1] = [1294,181,475]
    manualNodePoints[3,4,:,1] = [1091,522,546]
    manualNodePoints[3,3,:,1] = [979,696,587]
    manualNodePoints[3,2,:,1] = [1133,695,848]
    manualNodePoints[3,1,:,1] = [1298,723,1070]
    manualNodePoints[3,0,:,1] = [1728,642,1239]

    manualNodePoints[4,7,:,1] = [1135,107,1058]
    manualNodePoints[4,6,:,1] = [1020,105,704]
    manualNodePoints[4,5,:,1] = [1055,166,464]
    manualNodePoints[4,4,:,1] = [990,414,354]
    manualNodePoints[4,3,:,1] = [892,678,533]
    manualNodePoints[4,2,:,1] = [970,759,731]
    manualNodePoints[4,1,:,1] = [1083,571,1048]
    manualNodePoints[4,0,:,1] = [1261,270,1383]

    manualNodePoints[5,7,:,1] = [723,73,1113]
    manualNodePoints[5,6,:,1] = [554,53,394]
    manualNodePoints[5,5,:,1] = [569,131,131]
    manualNodePoints[5,4,:,1] = [764,275,113]
    manualNodePoints[5,3,:,1] = [780,655,400]
    manualNodePoints[5,2,:,1] = [853,785,606]
    manualNodePoints[5,1,:,1] = [820,718,872]
    manualNodePoints[5,0,:,1] = [574,445,993]

    manualNodePoints[6,7,:,1] = [181,263,450]
    manualNodePoints[6,6,:,1] = [-39,384,209]
    manualNodePoints[6,5,:,1] = [263,576,175]
    manualNodePoints[6,4,:,1] = [515,772,381]
    manualNodePoints[6,3,:,1] = [653,804,507]
    manualNodePoints[6,2,:,1] = [707,774,661]
    manualNodePoints[6,1,:,1] = [608,586,714]
    manualNodePoints[6,0,:,1] = [369,426,614]

    manualNodePoints[7,0,:,1] = [545,667,855]
    manualNodePoints[7,1,:,1] = [655,872,855]
    manualNodePoints[7,2,:,1] = [525,940,855]
    manualNodePoints[7,3,:,1] = [320,950,855]
    manualNodePoints[7,4,:,1] = [75,840,855]
    manualNodePoints[7,5,:,1] = [-70,710,855]
    manualNodePoints[7,6,:,1] = [90,550,855]
    manualNodePoints[7,7,:,1] = [242,590,855]

    manualNodePoints[8,7,:,1] = [675,400,1210]
    manualNodePoints[8,6,:,1] = [305,370,1210]
    manualNodePoints[8,5,:,1] = [70,580,1210]
    manualNodePoints[8,4,:,1] = [210,795,1210]
    manualNodePoints[8,3,:,1] = [635,700,1210]
    manualNodePoints[8,2,:,1] = [834,660,1210]
    manualNodePoints[8,1,:,1] = [1150,620,1210]
    manualNodePoints[8,0,:,1] = [940,430,1210]

#=========================================================#
#                                                         #
#   Derivatives of the whole models use the same fomula   #
#                                                         #
#=========================================================#

#calculating the derivatives 
difference = numpy.zeros((numberOfLengthNodes,8,3,2))
differenceAverage = numpy.zeros((numberOfLengthNodes,8,3,2))
circumDeriv = numpy.zeros((numberOfLengthNodes,8,3,2))
directDeriv = numpy.zeros((numberOfLengthNodes,8,3,2))
lengthDeriv = numpy.zeros((numberOfLengthNodes,8,3,2))
#circumferential derivative to be calculated 
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            if (i<7):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i+1,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,0,m,k]-manualNodePoints[j,7,m,k]
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            if (i<7):
                for m in range (3):
                    differenceAverage[j,i+1,m,k]=(difference[j,i+1,m,k]+difference[j,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,0,m,k]=(difference[j,0,m,k]+difference[j,7,m,k])/2
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            for m in range (3):
                circumDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# derivative of the length direction
for k in range (2):
    for i in range (8):
        for j in range (numberOfLengthNodes):
            if (j<numberOfLengthNodes-1):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j+1,i,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i,m,k]-manualNodePoints[j-1,i,m,k]
for k in range (2):
    for i in range (8):
        for j in range (numberOfLengthNodes):
            if (j == 0):
                for m in range (3): 
                    differenceAverage[j,i,m,k]=difference[j,i,m,k]
            if (j<numberOfLengthNodes-1):
                for m in range (3):
                    differenceAverage[j+1,i,m,k]=(difference[j,i,m,k]+difference[j+1,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,i,m,k]=difference[j-1,i,m,k]
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            for m in range (3):
                lengthDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# the derivatives of the wall direction is defined in the below lines ... 
for i in range (8):
    for j in range (numberOfLengthNodes):
        for m in range (3):
            for k in range (2):
                difference[j,i,m,k] = manualNodePoints[j,i,m,1] - manualNodePoints[j,i,m,0]
for i in range (8):
    for j in range (numberOfLengthNodes):

        for k in range (2):
            for m in range (3):
                directDeriv[j,i,m,k] = difference[j,i,m,k]/math.sqrt(math.pow(difference[j,i,0,k],2) + math.pow(difference[j,i,1,k],2) + math.pow(difference[j,i,2,k],2))

#=================================================================#
#                                                                 #
#             G e o m e t r i c      F i e l d                    #      
#                                                                 #
#=================================================================#

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
for dimension in range(3):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,dimension+1,1)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()

# Get nodes
nodes = iron.Nodes()
region.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes

# Get or calculate geometric parameters
if (exfileMesh):
    # Read the geometric field from the exnode file
    for node_num in range(1, exnode.num_nodes + 1):
        for derivative in range(1,9):
            version = 1
            for component in range(1, numberOfDimensions + 1):
                component_name = ["x", "y", "z"][component - 1]
                value = exnode.node_value("Coordinate", component_name, node_num, derivative)
                geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                      version, derivative, node_num, component, value)
    geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
else:
    # Create the geometric field
    for wallNodeIdx in range(1,numberOfWallNodes+1):
        for lengthNodeIdx in range(1,numberOfLengthNodes+1):
            for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
                x = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                y = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                z = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                xtangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                ytangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                ztangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                xnormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                ynormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                znormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                zxnormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                zynormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                zznormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,1,x)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,2,y)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,3,z)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,ztangent)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,zxnormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,zynormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,zznormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,znormal)
    # Update the geometric field
    geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    # Export undeformed mesh geometry
    print("Writing undeformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("Kaufman" + str(KaufmanNumber) + "UndeformedGeometry","FORTRAN")
    fields.ElementsExport("Kaufman" + str(KaufmanNumber) + "UndeformedGeometry","FORTRAN")
    fields.Finalise()


#=================================================================#
#                                                                 #
#                D a t a       P o i n t s                        #
#                                                                 #
#=================================================================#

# Create the data points
dataPoints = iron.DataPoints()
dataPoints.CreateStart(dataPointsUserNumber,region,numberOfDataPoints)
dataPointLocations = numpy.zeros((numberOfDataPoints,3))
print("Number of data points: " + str(numberOfDataPoints))
# reading from a text file containing the point clouds   
# reading from stage Kaufman 10 
if (Kaufman10): 
    if (Epi):
        with open("Kaufman10EpiDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)
    else: 
        with open("Kaufman10EndoDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)

# reading from stage Kaufman 11
if (Kaufman11): 
    if (Epi):
        with open("Kaufman11EpiDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)
    else: 
        with open("Kaufman11EndoDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)

# reading from stage Kaufman 12
if (Kaufman12): 
    if (Epi):
        with open("Kaufman12EpiDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)
    else: 
        with open("Kaufman12EndoDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)

# reading from stage Kaufman 13
if (Kaufman13): 
    if (Epi):
        with open("Kaufman13EpiDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)
    else: 
        with open("Kaufman13EndoDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)

# reading from stage Kaufman 14
if (Kaufman14): 
    if (Epi):
        with open("Kaufman14EpiDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)
    else: 
        with open("Kaufman14EndoDataPoints.txt", "r") as ins:
	        arrayOfInputData = []
	        for line in ins:
		        arrayOfInputData.append(line)

x = 0.0
y = 0.0
z = 0.0
#for i in range (startpoint, numberOfDataPoints + startpoint):
for i in range (numberOfDataPoints):
	for j in range (5):
		sample = arrayOfInputData[i*5 + j]
#		sample = arrayOfInputData[i*25 + j]
		if (math.fmod(j,5) == 1):
			x = float (sample[12:25])				
		elif (math.fmod(j,5) == 2):
			y = float (sample[12:25])
		elif (math.fmod(j,5) == 3):
			z = float (sample[12:17])
#		dataPointLocations[i - startpoint,:] = [x,y,z]
		dataPointLocations[i,:] = [x,y,z]


# Set up data points with geometric values
for dataPoint in range(numberOfDataPoints):
    dataPointId = dataPoint + 1
    dataList = dataPointLocations[dataPoint,:]
    dataPoints.PositionSet(dataPointId,dataList)
dataPoints.CreateFinish()

#=================================================================#
#                                                                 #
# D a t a   P r o j e c t i o n  on  G e o m e t r i c  F i e l d #
#                                                                 #
#=================================================================#
print("Projecting data points onto geometric field")
candidateElements = range(1,numberOfElements+1)
if (Epi): 
    candidateFaceNormals = iron.ElementNormalXiDirections.PLUS_XI3*numpy.ones(numberOfElements,dtype=numpy.int32)
else: 
    candidateFaceNormals = iron.ElementNormalXiDirections.MINUS_XI3*numpy.ones(numberOfElements,dtype=numpy.int32)
# Set up data projection
dataProjection = iron.DataProjection()
dataProjection.CreateStart(dataProjectionUserNumber,dataPoints,geometricField,iron.FieldVariableTypes.U)
#dataProjection.projectionType = iron.DataProjectionProjectionTypes.ALL_ELEMENTS
dataProjection.projectionType = iron.DataProjectionProjectionTypes.BOUNDARY_FACES
dataProjection.ProjectionCandidateFacesSet(candidateElements,candidateFaceNormals)
#dataProjection.ProjectionDataCandidateFacesSet([1,2,3],[1,2],[iron.ElementNormalXiDirections.PLUS_XI3,iron.ElementNormalXiDirections.PLUS_XI3])
dataProjection.CreateFinish()

# Evaluate data projection based on geometric field
dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
# Create mesh topology for data projection
mesh.TopologyDataPointsCalculateProjection(dataProjection)
# Create decomposition topology for data projection
decomposition.TopologyDataProjectionCalculate()

# Cancel some projections
dataProjection.ProjectionCancelByDistance(iron.DataProjectionDistanceRelations.GREATER_EQUAL,20.0)

# Output data projection results
dataProjection.ResultAnalysisOutput("ProjectionAnalysis")

rmsError=dataProjection.ResultRMSErrorGet()
print("RMS error = "+ str(rmsError))

# Output the .exdata file.                                           
dataErrorVector = numpy.zeros((numberOfDataPoints,3))
dataErrorDistance = numpy.zeros(numberOfDataPoints)
for elementIdx in range(1,numberOfElements+1):
    numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementIdx)
    for dataPointIdx in range(1,numberOfProjectedDataPoints+1):
        dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementIdx,dataPointIdx)
        errorVector = dataProjection.ResultProjectionVectorGet(dataPointNumber,3)
        dataErrorVector[dataPointNumber-1,0]=errorVector[0]
        dataErrorVector[dataPointNumber-1,1]=errorVector[1]
        dataErrorVector[dataPointNumber-1,2]=errorVector[2]
        errorDistance = dataProjection.ResultDistanceGet(dataPointNumber)
        dataErrorDistance[dataPointNumber-1]=errorDistance
 
# write data points to exdata file for CMGUI
offset = 0
writeExdataFile("DataPoints"+".Kaufman"+str(KaufmanNumber)+".part"+str(computationalNodeNumber)+".exdata",dataPointLocations,dataErrorVector,dataErrorDistance,offset)
print("Projection complete")
#exit(0)

#=================================================================#
#                                                                 #
#             E q u a t i o n  s         S e t                    #
#                                                                 #
#=================================================================#

# Create vector fitting equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                             iron.EquationsSetTypes.DATA_FITTING_EQUATION,
                             iron.EquationsSetSubtypes.DATA_POINT_FITTING, 
 			     iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

#=================================================================#
#                                                                 #
#           D e p e n d e n  t       F i e l d                    #
#                                                                 #   
#=================================================================#

# Create dependent field (will be deformed fitted values based on data point locations)
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,numberOfDimensions)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions)
equationsSet.DependentCreateFinish()
# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

# Initialise dependent field to undeformed geometric field
for component in range (1,numberOfDimensions+1):
    geometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            component, dependentField, iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, component)


#=================================================================#
#                                                                 #
#           I n d e p e n d e n t        F i e l d                #
#                                                                 #
#=================================================================#

# Create data point field (independent field, with vector values stored at the data points)
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
independentField.VariableLabelSet(iron.FieldVariableTypes.U,"data point vector")
independentField.VariableLabelSet(iron.FieldVariableTypes.V,"data point weight")
independentField.DataProjectionSet(dataProjection)
equationsSet.IndependentCreateFinish()
# Initialise data point vector field to 0
#independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
# Initialise data point weight field to 1
#independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,1.0)
# loop over each element's data points and set independent field values to data point locations on surface of the sphere
for element in range(numberOfElements):
    elementId = element + 1
    elementDomain = decomposition.ElementDomainGet(elementId)
    if (elementDomain == computationalNodeNumber):
        numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementId)
        for dataPoint in range(numberOfProjectedDataPoints):
            dataPointId = dataPoint + 1
            dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementId,dataPointId)
            dataList = dataPoints.PositionGet(dataPointNumber,3)
            # set data point field values
            for component in range(numberOfDimensions):
                componentId = component + 1
                dataPointNumberIndex = dataPointNumber - 1
                value = dataList[component]
                independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementId,dataPointId,componentId,value)

#=================================================================
# Material Field
#=================================================================
# Create material field (Sobolev parameters)
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Smoothing Parameters")
equationsSet.MaterialsCreateFinish()
# Set kappa and tau - Sobolev smoothing parameters
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

#=================================================================
# Equations
#=================================================================
# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.FULL
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#=================================================================
# Problem setup
#=================================================================
# Create fitting problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FITTING,
                        iron.ProblemTypes.DATA_FITTING,
                        iron.ProblemSubtypes.STATIC_FITTING]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.NONE # NONE / MATRIX
#solver.outputType = iron.SolverOutputTypes.MATRIX # NONE / MATRIX
solver.linearType = iron.LinearSolverTypes.ITERATIVE
#solver.linearType = iron.LinearSolverTypes.DIRECT
#solver.LibraryTypeSet(iron.SolverLibraries.UMFPACK) # UMFPACK/SUPERLU
#solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
solver.linearIterativeAbsoluteTolerance = 1.0E-10
solver.linearIterativeRelativeTolerance = 1.0E-05
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
#solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.FULL
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#=================================================================
# Boundary Conditions
#=================================================================


# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
'''
for nodeIdx in range(73,73+32):
    for componentIdx in range(1,4):
        for derivativeIdx in range(1,9):
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
for nodeIdx in range(144-32,145):
    for componentIdx in range(1,4):
        for derivativeIdx in range(1,9):
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
'''
# for nodeIdx in range(numberOfLengthNodes*numberOfCircumfrentialNodes+1,numberOfLengthNodes*numberOfCircumfrentialNodes*2+1):
if (Epi):
    for nodeIdx in range(1,numberOfLengthNodes*numberOfCircumfrentialNodes+1):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
else: 
    for nodeIdx in range(numberOfLengthNodes*numberOfCircumfrentialNodes+1,numberOfLengthNodes*numberOfCircumfrentialNodes*numberOfWallNodes+1):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
#for nodeIdx in range(5,9):
#    for componentIdx in range(1,4):
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeIdx,2,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeIdx,1,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeIdx,2,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeIdx,3,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
       
solverEquations.BoundaryConditionsCreateFinish()


#=================================================================
# S o l v e    a n d    E x p o r t    D a t a
#=================================================================
derivativeVector=[0.0,0.0,0.0,0.0]
numberOfIterations = 5
for iteration in range (startIteration,startIteration+numberOfIterations+1):
    # Solve the problem
    print("Solving fitting problem, iteration: " + str(iteration))
    problem.Solve()
    # Normalise derivatives
    for nodeIdx in range(1,numberOfNodes+1):
      for derivativeIdx in [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]:
          length=0.0
          for componentIdx in range(1,4):
              derivativeVector[componentIdx]=dependentField.ParameterSetGetNode(iron.FieldVariableTypes.U,
                                                                                iron.FieldParameterSetTypes.VALUES,
                                                                                1,derivativeIdx,nodeIdx,componentIdx)
              length=length + derivativeVector[componentIdx]*derivativeVector[componentIdx]
          length=math.sqrt(length)
          for componentIdx in range(1,4):
              value=derivativeVector[componentIdx]/length
              dependentField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,derivativeIdx,nodeIdx,componentIdx,value)
    # Copy dependent field to geometric 
    for componentIdx in range(1,numberOfDimensions+1):
        dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES,
                                                                componentIdx,geometricField,
                                                                iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                componentIdx)
    # Reproject
    dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
    rmsError=dataProjection.ResultRMSErrorGet()
    print("RMS error = "+ str(rmsError))
    # Export fields
    print("Writing deformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("Kaufman" + str(KaufmanNumber) + "DeformedGeometry" + str(iteration),"FORTRAN")
    fields.ElementsExport("Kaufman" + str(KaufmanNumber) + "DeformedGeometry" + str(iteration),"FORTRAN")
    fields.Finalise()
#-----------------------------------------------------------------


iron.Finalise()
