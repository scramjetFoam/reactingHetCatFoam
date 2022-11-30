#!/bin/sh

stitchMesh -perfect -overwrite masterPerp slavePerp >> log.stitchMesh_0
stitchMesh -perfect -overwrite masterHor slaveHor >> log.stitchMesh_1
stitchMesh -perfect -overwrite masterVer slaveVer >> log.stitchMesh_2
