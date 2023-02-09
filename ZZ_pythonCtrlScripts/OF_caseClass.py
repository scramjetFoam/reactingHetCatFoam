# -- Python class to store, change and write OpenFOAM case data

# --    usage:
# --        1) load OpenFOAMCase from baseCase or from the different OpenFOAMCase
# --            -- loadOFCaseFromBaseCase(baseDir)
# --                    baseDir -- base case directory
# --            -- loadOFCaseFromOFCase(parent)
# --                    parent -- OpenFOAMCase you want to load from
# --                NOTE: does not work so far
# --        2) (optional) change OpenFOAMCase.dir (where the changes will be made), default -- baseDir
# --            -- changeOFCaseDir(dir)
# --        3) (optional) copy OpenFOAMCase.basedir to OpenFOAMCase.dir
# --            -- copyBaseCase()
# --        4) specify what you want to modify in openFoam source files:
# --            options:            
# --                a) replace(in, what, by)
# --                        inFl -- openFoam source file you want to change
# --                        what -- what you want to replace (string)
# --                        by -- new value (string)
# --        5) run commands in OpenFOAMCase.dir
# --            -- runCommands(commands)

# --    additional functions:
# --        -- updateTimes() 
# --                update OpenFOAMCase.times variable (includes float folders in OpenFOAMCase.dir)
# --                and OpenFOAMCase.latestTime variable (max of OpenFOAMCase.times)

# --    TODO:
# --        1) connect with Lucka blockMeshDict class

# -- imports 
import numpy as np
import os
import shutil as sh
from myAddFcs import *

class OpenFOAMCase:
    def __init__(self):
        """initialization of the new OpenFOAMCase"""
        # -- save starting directory
        self.whereIStart = os.getcwd()
        print("New OpenFOAMCase has been initialized, starting directory is:\n\t%s"%self.whereIStart)
    
    def loadOFCaseFromBaseCase(self, baseDir):
        """load OpenFOAMCase from base case directory"""
        print("OpenFOAMCase.baseCase has been set to :\n\t%s" % baseDir)
        self.baseDir = baseDir
        # -- OpenFOAMCase directory is initialized as its baseDir, can be changed by changeOFCaseDir(dir)
        self.dir = baseDir
    
    def loadOFCaseFromOFCase(self, parent):
        """load OpenFOAMCase from the different OpenFOAMCase"""
        print("OpenFOAMCase has been loaded from parent OpenFOAMCase.")
    
    def changeOFCaseDir(self,dir):
        """change OpenFOAMCase directory"""
        self.dir = dir
        print("OpenFOAMCase directory has been changed from:\n\t%s to %s"%(self.baseDir,dir))

    def copyBaseCase(self):
        """copy OpenFOAMCase.baseDir to OpenFOAMCase.dir"""
        # -- if OpenFOAMCase.dir exists, remove it
        if os.path.isdir(self.dir): 
            sh.rmtree(self.dir)

        # -- copy baseDir to dir
        sh.copytree(self.baseDir,self.dir)
        print("OpenFOAMCase.baseDir has been copied to OpenFOAMCase.dir:\n\t %s --> %s"%(self.baseDir,self.dir))

    def replace(self, replaces):
        """replace option -- replace = [in inFl (file), what (string), by (string)]"""
        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        # -- make the replaces
        for replace in replaces:
            inFl, what, by = replace
            with open(inFl, 'r') as fl:
                linesInFl = fl.readlines()
            for lnI in range(len(linesInFl)):
                for whatI in range(len(what)):
                    whatTu = what[whatI]
                    byTu = by[whatI]
                    if whatTu in linesInFl[lnI]:
                        linesInFl[lnI] = linesInFl[lnI].replace(whatTu,byTu)
                        print("In %s, I have replaced %s by %s on line %d." % (inFl, whatTu, byTu, lnI))
            with open(inFl, 'w') as fl:
                for lnI in range(len(linesInFl)):
                    fl.writelines(linesInFl[lnI])
        
        # -- move back where I start
        os.chdir(self.whereIStart)
    
    def setParameter(self, inParVal):
        """setParameter option -- inParVal = [in inFl (file), par (string), val (string), inSubDict (string)]"""
        """inSubDict option enables to find parameter in subDictionary, if "" then does nothing """
        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        # -- make the replaces
        inFl, par, val, inSubDict = inParVal
        with open(inFl, 'r') as fl:
            linesInFl = fl.readlines()
        
        # -- if inSubDict speciefied find where it is
        if not inSubDict == "":
            subDictSt, subDictEnd = -1, -1
            for lnI in range(len(linesInFl)):
                if inSubDict in linesInFl[lnI]:
                    subDictSt = lnI
                if subDictSt != -1 and subDictEnd == -1 and "}" in linesInFl[lnI]:
                    subDictEnd = lnI
            if subDictSt == -1:
                print("I could not find subDict %s in file %s."%(inSubDict,inFl))
        for lnI in range(len(linesInFl)):
            if par in linesInFl[lnI] and (inSubDict == "" or (lnI >= subDictSt and lnI <= subDictEnd)):
                try:
                    oldVal = linesInFl[lnI].split(par)[1].replace(' ','').replace('\n','').replace(';','').replace('\t','')
                    if not oldVal == '':
                        linesInFl[lnI] = linesInFl[lnI].replace(oldVal,val)
                        print("In %s, I have set parameter %s in subDictionary '%s' from value %s to value %s on line %d." % (inFl, par, inSubDict,oldVal, val, lnI))
                    else: print("I could not replaced parameter %s on line %d:\n\t'%s'"%(par,lnI,linesInFl[lnI].replace('\n','')))
                except: print("I could not replaced parameter %s on line %d:\n\t'%s'"%(par,lnI,linesInFl[lnI].replace('\n','')))
        with open(inFl, 'w') as fl:
            for lnI in range(len(linesInFl)):
                fl.writelines(linesInFl[lnI])

        # -- move back where I start
        os.chdir(self.whereIStart)
    
    def runCommands(self,commands):
        """run commands in OpenFOAMCase dir"""
        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        # -- runCommand
        print("I will execute commands in %s:"%self.dir)
        for commandI in range(len(commands)):
            command = commands[commandI]
            print("\t%d. %s"%(commandI+1,command))
            os.system(command)
        
        # -- move back where I start
        os.chdir(self.whereIStart)
    
    def updateTimes(self):
        """update the OpenFOAMCase.times variable"""

        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        lstDir = os.listdir()
        flLstDir = [float(flDir) for flDir in lstDir if isFloat(flDir)]
        self.times = flLstDir
        self.latestTime = max(flLstDir)
        print("Found OpenFOAMCase times are:\n\t%s,\nlatestTime:\n\t%g"%(str(self.times),self.latestTime))

        # -- move back where I start
        os.chdir(self.whereIStart)  

    
