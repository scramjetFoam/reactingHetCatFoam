#!/bin/bash
# Name of the job
#SBATCH -J rashV7
# Partition to use
#SBATCH -p Mlong
#SBATCH -o log.slurm
# Time limit. Often not needed as there are defaults, 
# but you might have to specify it to get the maximum allowed.
# time: 24hours
#SBATCH --time=14-00:00:0
# Number of processes
#SBATCH -n16
#SBATCH -N1
# Do NOT use some nodes
#SBATCH --exclude=kraken-m[1]
#! ## SBATCH --nodelist=kraken-m9

# run setup
# load the modules
module load openfoam-org/6-10.3.0
source /opt/modules/spack_installs/linux-centos7-broadwell/gcc-10.3.0/openfoam-org-6-vqdvitbtcrbdg66aa2x2zpf4jtpcvqks/etc/bashrc


#prepare Geom
# blockMesh >> log.blockMesh

# decomposePar >> log.decomposePar

# srun -n16 snappyHexMesh -parallel -overwrite >> log.snappyHexMesh
# reconstructParMesh -mergeTol 1 -constant

# rm -rf 0

# paraFoam -touch

# mkdir 0
# cp -rf 0.org/* 0

# rm -rf 0/meshPhi

# python genGeomFPMAll.py

# rm -rf processor*


# decomposePar >> log.decomposePar

# run SIm
srun -n16 nonIsoFlowDCR -parallel >> log.nonIsoFlowDCR



# blockMesh >> log.blockMesh

# echo 'stitching the mesh'
# chmod 777 stitchMeshSc.sh 
# bash stitchMeshSc.sh 

# # Note: stitchMesh now returns warning, but it all seems to end well

# createPatch -overwrite >> log.createPatch1

# #~ runApplication refineWallLayer -overwrite "(walls)" 0.5
# #~ mv log.refineWallLayer log.refineWallLayer_0
# #~ runApplication refineWallLayer -overwrite "(walls)" 0.5
# #~ mv log.refineWallLayer log.refineWallLayer_1
# #~ runApplication refineWallLayer -overwrite "(walls)" 0.5
# #~ mv log.refineWallLayer log.refineWallLayer_2
# decomposePar >> log.decomposePar1
# srun -n16 snappyHexMesh -overwrite -parallel >> log.snappyHexMesh
# reconstructParMesh -constant -mergeTol 1 >> log.reconstructParMesh
# topoSet >> log.topoSet1
# createPatch -dict system/createPatchDict.2 -overwrite >> log.createPatch2

# #~ runApplication refineMesh -overwrite

# paraFoam -touch

# mkdir 0
# cp -rf 0.org/* 0

# rm -rf 0/meshPhi

# python genGeomFPMAll.py

# topoSet -dict system/topoSetDict.3 >> log.topoSet2

# rm -rf processor*
# decomposePar >> log.decomposePar2

# ## runApplication renumberMesh -overwrite -frontWidth


#srun -n16 nonIsoFlowDCR -parallel >> log.nonIsoFlowDCR

# application=`getApplication`
# runParallel $application



