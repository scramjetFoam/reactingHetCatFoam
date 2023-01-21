from os import system, path
from sys import exit

solver = False
np = False
# np = '4'
# solver = 'reactingHetCatSimpleFoam'

with open('system/controlDict', mode='r') as f1:
    lines = f1.readlines()
    for line in lines:
        if 'application' in line and line[0] != '/':
            solver = line.split(' ')[-1][:-2]
            # print(solver)

with open('system/decomposeParDict', mode='r') as f2:
    lines = f2.readlines()
    for line in lines:
        if 'numberOfSubdomains' in line and line[0] != '/':
            np = line.split('\t')[-1][:-2]
            # print(np)

if (not np) or (not solver):
    print('solver or np not found')
    exit()

system('rm -rf 0')
system('cp -rf 0.org 0')

print('1/3 generating mesh')
system('blockMesh > log.blockMesh')
system('decomposePar > log.decomposePar1')
system(f'mpirun -np {np} snappyHexMesh -overwrite -parallel > log.snappyHexMesh')
system('reconstructParMesh -constant -time 0 > log.reconstructParMesh1')


print('2/3 running solver')
system('decomposePar -force > log.decomposePar2')
system(f'mpirun -np {np} renumberMesh -overwrite -parallel > log.renumberMesh')
system(f'mpirun -np {np} {solver} -parallel > log.{solver}1')
if not path.isdir('processor0/100'):
    system(f'mpirun -np {np} {solver} -parallel > log.{solver}2')
system('reconstructPar -latestTime > log.reconstructPar')

print('3/3 postprocessing')
system('postProcess -func integrace -latestTime > log.integrace')
system('intSrcSphere > log.intSrcSphere')
system('rm -rf processor*')