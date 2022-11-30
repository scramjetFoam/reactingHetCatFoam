#!/bin/bash
# Name of the job
#SBATCH -J SpaciN 
# Partition to use
#SBATCH -p Mlong
# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed.
# time: 24hours
#SBATCH --time=14-00:00:0
# Number of processes
#SBATCH -n12
#SBATCH -N1
# Do NOT use some nodes
#SBATCH --exclude=kraken-m[1]
#### SBATCH --nodelist=kraken-m8

# run setup
# load the modules
module load openfoam-org/6-10.3.0
source /opt/modules/spack_installs/linux-centos7-broadwell/gcc-10.3.0/openfoam-org-6-vqdvitbtcrbdg66aa2x2zpf4jtpcvqks/etc/bashrc

module load python
module load py-numpy

rm -rf 0

blockMesh >> log.blockMesh

echo 'stitching the mesh'
chmod 777 stitchMeshSc.sh 
bash stitchMeshSc.sh 

# Note: stitchMesh now returns warning, but it all seems to end well

createPatch -overwrite >> log.createPatch1

#~ runApplication refineWallLayer -overwrite "(walls)" 0.5
#~ mv log.refineWallLayer log.refineWallLayer_0
#~ runApplication refineWallLayer -overwrite "(walls)" 0.5
#~ mv log.refineWallLayer log.refineWallLayer_1
#~ runApplication refineWallLayer -overwrite "(walls)" 0.5
#~ mv log.refineWallLayer log.refineWallLayer_2
decomposePar >> log.decomposePar1
srun -n12 snappyHexMesh -overwrite -parallel >> log.snappyHexMesh
reconstructParMesh -constant -mergeTol 1 >> log.reconstructParMesh
topoSet >> log.topoSet1
createPatch -dict system/createPatchDict.2 -overwrite >> log.createPatch2

#~ runApplication refineMesh -overwrite

paraFoam -touch

mkdir 0
cp -rf 0.org/* 0

rm -rf 0/meshPhi

# python genGeomFPMAll.py

topoSet -dict system/topoSetDict.3 >> log.topoSet2

rm -rf processor*
decomposePar >> log.decomposePar2

## runApplication renumberMesh -overwrite -frontWidth


srun -n12 nonIsoFlowDCRV2 -parallel >> log.nonIsoFlowDCRV2

# application=`getApplication`
# runParallel $application

reconstructPar -latestTime >> log.reconstructPar
postProcess -func 'patchAverage(name=inletCyl,CO)' -latestTime > log.postProcess


