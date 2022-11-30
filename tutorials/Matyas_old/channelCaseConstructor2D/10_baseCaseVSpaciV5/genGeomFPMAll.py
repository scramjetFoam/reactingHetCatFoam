#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  genGeomFPM.py
#  
#  Copyright 2020 Tomáš Hlavatý <hlavatyo@vscht.cz>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software  
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# This script creates isOnWallCoating and isInWallCoating for diffReactSolver
# note different settings for :
# "basic geometry" - without overlay of the coating along the channel
# "non-basic geometry$ - with overlay of the coating along the channel

import os
import numpy as np

coatOnWallFrac = 0.99
coatInWallFrac = 0.0

basicGeometry = 1														# 1 - non-overlay, 0 - overlay
if basicGeometry == 1:
	print ('Creating isOnWall and isInWall without overlay\n')
else:
	print ('Creating isOnWall and isInWall with overlay\n')
	
# -- spliting mesh into regions for in wall and on wall fields 
os.system("splitMeshRegions -detectOnly -cellZones > log.splitMesh")
with open('log.splitMesh','r') as fl:
	lines = fl.readlines()
wallReg = "-1\n"
freeReg = "-1\n"
stInd = -1
endInd = -1
for lnInd in range(len(lines)):
	if "Region\tCells\n" == lines[lnInd]:
		stInd = lnInd+2
	if "Region\tZone\tName\n" == lines[lnInd]:
		endInd = lnInd-2
nCells = np.empty((0,2),float)
for i in range(stInd,endInd):
	nCells = np.append(nCells,np.array([lines[i].replace('\n','').split('\t')]).astype(float).reshape(1,2),axis=0)

sort = (np.argsort(-nCells[:,1]))
# print(sort)
# print(nCells[sort][0,0])
freeReg = int(nCells[sort][0,0])
wallReg = int(nCells[sort][1,0])
print(freeReg,wallReg)

fSpMesh  = open('0/cellToRegion','r')
onWall = fSpMesh.readlines()
fSpMesh.close()
inWall = onWall[:]
wallFracs = inWall[:]

def isInt(val):
    try:
        int(val)
        return True
    except:
        return False
#creating fields isOnWallCoat and is InWallCoat
#NOTE: this has to be checked if new zones are added
if basicGeometry == 1:
	for ind in range(len(wallFracs)):
		if wallFracs[ind] == '%d\n'%wallReg:        
			wallFracs[ind] = '%g\n'%coatInWallFrac
		elif wallFracs[ind] == '%d\n'%freeReg:
			wallFracs[ind] = '0\n'
		elif isInt(wallFracs[ind].replace('\n','')):
			if int(wallFracs[ind].replace('\n','')) < 10000: 
				wallFracs[ind] = '%g\n'%coatOnWallFrac
	for ind in range(len(inWall)):
		if inWall[ind] == '%d\n'%wallReg:       
			inWall[ind] = '1\n'
		elif isInt(inWall[ind].replace('\n','')):
			if int(inWall[ind].replace('\n','')) < 10000: 
				inWall[ind] = '0\n'
	for ind in range(len(onWall)):
		if onWall[ind] == '%d\n'%wallReg or onWall[ind] == '%d\n'%freeReg:
			onWall[ind] = '0\n'
		elif isInt(onWall[ind].replace('\n','')):                    
			if int(onWall[ind].replace('\n','')) < 10000: 
				onWall[ind] = '1\n'
# else:
# 	for ind in range(len(inWall)):
# 		if inWall[ind] == '2\n' or inWall[ind] == '3\n'  or inWall[ind] == '4\n' or inWall[ind] == '5\n' or inWall[ind] == '6\n':       
# 			inWall[ind] = '0\n'
# 	for ind in range(len(onWall)):
# 		if onWall[ind] == '1\n' or onWall[ind] == '4\n':                                                              
# 			onWall[ind] = '0\n'
# 		if onWall[ind] == '3\n' or onWall[ind] == '2\n' or onWall[ind] == '6\n'  or onWall[ind] == '5\n':                    
# 			onWall[ind] = '1\n'
with open('0/coatWallFrac','w+') as isInWall:
    isInWall.writelines(wallFracs)
isInWall.close()
with open('0/isOnWallCoat','w+') as isOnWall:
    isOnWall.writelines(onWall)
isOnWall.close()
with open('0/isInWallCoat','w+') as isOnWall:
    isOnWall.writelines(inWall)
isOnWall.close()
os.system ('cp 0/coatWallFrac 0.org/')
os.system ('cp 0/isOnWallCoat 0.org/')
os.system ('cp 0/isInWallCoat 0.org/')
os.system ('cp 0/cellToRegion 0.org/')
