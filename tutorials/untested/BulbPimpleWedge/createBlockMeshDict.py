
import math
from blockMeshDictClassV8 import *


Vbulb = 77.99*1e-6

r1 = (Vbulb / 4 * 3 / math.pi) ** (1./3)

d1 = 2 * r1

# parameters
l1 = 0.025
l2 = 85.9e-3
l3 = 0.150

r0 = 0.00208

# d0 = 0.01
# d1 = 0.01044842266

# rAdd = d0*0.2

dY = 0.5e-3
dX = 1e-3
dZ = 1.0e-3

print(d1+l2/2)
print("d1", d1)

nCZ = 1

x0 = y0 = z0 = 0.0
grX = grY = grZ = "1.0"

yNas=5
yNas0 = 1.5
xNas = 2

gradX1 = "0.1"
gradX1I = "10"

wAng = 1.0 # wedge angle
yMax = dZ/math.atan(wAng/180*math.pi*0.5)

# create mesh
fvMesh = mesh()

### BLOCKS ###
# -- first block
xC, yC = x0, y0
xE, yE = xC+d1, yC+r0

# vertices
vertices = [
        [xC, yC, z0-yC/yMax*dZ],
        [xE, yC, z0-yC/yMax*dZ],
        [xE, yE, z0-yE/yMax*dZ],
        [xC, yE, z0-yE/yMax*dZ],
        [xC, yC, z0+yC/yMax*dZ],
        [xE, yC, z0+yC/yMax*dZ],
        [xE, yE, z0+yE/yMax*dZ],
        [xC, yE, z0+yE/yMax*dZ],
    ]

# neighbouring blocks
neighbours = []

# number of cells
nCX = int(round(abs(xE-xC)/dX))
nCY = int(round(abs(yE-yC)/dY))
nCells = [nCX, nCY*yNas0, nCZ]

# grading
grading = [gradX1, grY, grZ]

# create the block
firstMid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

# -- first top block
xC, yC = xC, yE
xE, yE = xE, yC+1e-3

# vertices
vertices = [
        [xC, yC, z0-yC/yMax*dZ],
        [xE, yC, z0-yC/yMax*dZ],
        [xE, yE, z0-yE/yMax*dZ],
        [xC, yE, z0-yE/yMax*dZ],
        [xC, yC, z0+yC/yMax*dZ],
        [xE, yC, z0+yC/yMax*dZ],
        [xE, yE, z0+yE/yMax*dZ],
        [xC, yE, z0+yE/yMax*dZ],
    ]

# neighbouring blocks
neighbours = [firstMid]

# number of cells
nCX = int(round(abs(xE-xC)/dX))
nCY = int(round(abs(yE-yC)/dY))
nCells = [nCX, nCY*yNas, nCZ]

# grading
grading = [gradX1, grY, grZ]

# create the block
firstTop = fvMesh.addBlock(vertices, neighbours, nCells, grading)

fvMesh.addEdge("arc", firstTop.retEYEZ0(), [((xE+xC)*0.5, yE+r1, z0-(yE+r1)/yMax*dZ)])
fvMesh.addEdge("arc", firstTop.retEYEZE(), [((xE+xC)*0.5, yE+r1, z0+(yE+r1)/yMax*dZ)])

# -- first top block
xC, yC = xE, y0
xE, yE = xC+l2, yC+r0

# vertices
vertices = [
        [xC, yC, z0-yC/yMax*dZ],
        [xE, yC, z0-yC/yMax*dZ],
        [xE, yE, z0-yE/yMax*dZ],
        [xC, yE, z0-yE/yMax*dZ],
        [xC, yC, z0+yC/yMax*dZ],
        [xE, yC, z0+yC/yMax*dZ],
        [xE, yE, z0+yE/yMax*dZ],
        [xC, yE, z0+yE/yMax*dZ],
    ]

# neighbouring blocks
neighbours = [firstMid]

# number of cells
nCX = int(round(abs(xE-xC)/dX))
nCY = int(round(abs(yE-yC)/dY))
nCells = [nCX*xNas, nCY*yNas0, nCZ]

# grading
grading = [grX, grY, grZ]

# create the block
secondMid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

# -- first top block
xC, yC = xE, y0
xE, yE = xC+d1, yC+r0

# vertices
vertices = [
        [xC, yC, z0-yC/yMax*dZ],
        [xE, yC, z0-yC/yMax*dZ],
        [xE, yE, z0-yE/yMax*dZ],
        [xC, yE, z0-yE/yMax*dZ],
        [xC, yC, z0+yC/yMax*dZ],
        [xE, yC, z0+yC/yMax*dZ],
        [xE, yE, z0+yE/yMax*dZ],
        [xC, yE, z0+yE/yMax*dZ],
    ]

# neighbouring blocks
neighbours = [secondMid]

# number of cells
nCX = int(round(abs(xE-xC)/dX))
nCY = int(round(abs(yE-yC)/dY))
nCells = [nCX, nCY*yNas0, nCZ]

# grading
grading = [gradX1I, grY, grZ]

# create the block
thirdMid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

# -- first top block
xC, yC = xC, yE
xE, yE = xE, yC+1e-3

# vertices
vertices = [
        [xC, yC, z0-yC/yMax*dZ],
        [xE, yC, z0-yC/yMax*dZ],
        [xE, yE, z0-yE/yMax*dZ],
        [xC, yE, z0-yE/yMax*dZ],
        [xC, yC, z0+yC/yMax*dZ],
        [xE, yC, z0+yC/yMax*dZ],
        [xE, yE, z0+yE/yMax*dZ],
        [xC, yE, z0+yE/yMax*dZ],
    ]

# neighbouring blocks
neighbours = [thirdMid]

# number of cells
nCX = int(round(abs(xE-xC)/dX))
nCY = int(round(abs(yE-yC)/dY))
nCells = [nCX, nCY*yNas, nCZ]

# grading
grading = [gradX1I, grY, grZ]

# create the block
thirdTop = fvMesh.addBlock(vertices, neighbours, nCells, grading)

fvMesh.addEdge("arc", thirdTop.retEYEZ0(), [((xE+xC)*0.5, yE+r1, z0-(yE+r1)/yMax*dZ)])
fvMesh.addEdge("arc", thirdTop.retEYEZE(), [((xE+xC)*0.5, yE+r1, z0+(yE+r1)/yMax*dZ)])

# # -- second block
# xC, yC = xE+r1, y0
# xE, yE = xC+l2, yC+r0

# # vertices
# vertices = [
#         [xC, yC, z0-yC/yMax*dZ],
#         [xE, yC, z0-yC/yMax*dZ],
#         [xE, yE, z0-yE/yMax*dZ],
#         [xC, yE, z0-yE/yMax*dZ],
#         [xC, yC, z0+yC/yMax*dZ],
#         [xE, yC, z0+yC/yMax*dZ],
#         [xE, yE, z0+yE/yMax*dZ],
#         [xC, yE, z0+yE/yMax*dZ],
#     ]

# # neighbouring blocks
# neighbours = [firstMid]

# # number of cells
# nCX = int(round(abs(xE-xC)/dX))
# nCY = int(round(abs(yE-yC)/dY))
# nCells = [nCX, nCY, nCZ]

# # grading
# grading = [grX, grY, grZ]

# # create the block
# secondMid = fvMesh.addBlock(vertices, neighbours, nCells, grading)



# # -- third block
# xC, yC = xE, y0
# xE, yE = x0+l3, yC+r0

# # vertices
# vertices = [
#         [xC, yC, z0-yC/yMax*dZ],
#         [xE, yC, z0-yC/yMax*dZ],
#         [xE, yE, z0-yE/yMax*dZ],
#         [xC, yE, z0-yE/yMax*dZ],
#         [xC, yC, z0+yC/yMax*dZ],
#         [xE, yC, z0+yC/yMax*dZ],
#         [xE, yE, z0+yE/yMax*dZ],
#         [xC, yE, z0+yE/yMax*dZ],
#     ]

# # neighbouring blocks
# neighbours = [secondMid]

# # number of cells
# nCX = int(round(abs(xE-xC)/dX))
# nCY = int(round(abs(yE-yC)/dY))
# nCells = [nCX, nCY, nCZ]

# # grading
# grading = [grX, grY, grZ]

# # create the block
# thirdMid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

### PATCHES ###
# -- wedge patches
wedgeZ0 = list()
for block in fvMesh.blocks:
    wedgeZ0.append(block.retFXY0())

fvMesh.addPatch("wedgeZ0", "wedge", wedgeZ0)

wedgeZE = list()
for block in fvMesh.blocks:
    wedgeZE.append(block.retFXYE())

fvMesh.addPatch("wedgeZE", "wedge", wedgeZE)

# # -- inlet
# inlet = list()
# inlet.append(firstMid.retFYZ0())

# fvMesh.addPatch("inlet", "patch", inlet)

# -- outlet
# outlet = list()
# outlet.append(thirdMid.retFYZE())

# fvMesh.addPatch("outlet", "patch", outlet)

# -- walls
walls = list()
walls.append(firstMid.retFYZ0())
walls.append(firstTop.retFYZ0())
walls.append(firstTop.retFYZE())
# walls.append(firstMid.retFXZE())

walls.append(firstTop.retFXZE())

walls.append(secondMid.retFXZE())

walls.append(thirdMid.retFYZE())
walls.append(thirdTop.retFYZ0())
walls.append(thirdTop.retFYZE())
walls.append(thirdTop.retFXZE())

# walls.append(secondTop.retFYZ0())
# walls.append(secondTop.retFXZE())
# walls.append(secondTop.retFYZE())

# walls.append(thirdMid.retFXZ0())
# walls.append(thirdMid.retFXZE())

fvMesh.addPatch("walls", "wall", walls)

### WRITE ###
fvMesh.writeBMD("./system/")
