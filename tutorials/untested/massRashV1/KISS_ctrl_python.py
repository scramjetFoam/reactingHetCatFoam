# control for reactor sim

import os

def isInt(value):
    try:
        int(value)
        return True
    except:
        return False

def latestTime():
    lstDir = os.listdir('./')
    lstInt = [lstdir for lstdir in lstDir if isInt(lstdir)]
    return int(max((lstInt)))

def changeInControlDict(endT):
    with open('system/controlDict','r') as fl:
        lines = fl.readlines()
    for i in range(len(lines)):
        if 'endlineLine' in lines[i]:
            lines[i] = 'endTime	       %d;  //endlineLine\n'%endT
        if 'writeInterval' in lines[i]:
            lines[i] = 'writeInterval	       %d;\n'%endT
    with open('system/controlDict','w') as fl:
        for i in range(len(lines)):
            fl.writelines(lines[i])        



N = 3
# simDir = 'rashigsV1'
# simDir = 'testMappedBC'
simDir = './'
scalarfields = ['T','ethylene']
vectorfields = ['U']
fields = "'(T ethylene U)'"
nIt = 2
# N = 33

os.chdir(simDir)

for sim in range(N):
    print('Running sim. %d %g -- %g m'%(sim,(sim)*0.1,(sim+1)*0.1))
    os.system('foamJob -parallel -screen reactingHetCatSimpleFoam > log.rash%d'%sim)
    os.system('reconstructPar -latestTime -fields %s> log.reconstruct%d'%(fields,sim))
    os.system("postProcess -func 'sample' > log.sample%d"%sim)
    lT = latestTime()
    print('Result in %d folder'%lT)
    os.system('mkdir constant/boundaryData')
    os.system('mkdir constant/boundaryData/inlet')
    os.system('mkdir constant/boundaryData/inlet/%d'%lT)
    os.system('cp postProcessing/sample/%d/outlet_field/points constant/boundaryData/inlet/%d/'%(lT,lT))
    for field in scalarfields:
        os.system('cp postProcessing/sample/%d/outlet_field/scalarField/%s constant/boundaryData/inlet/%d/'%(lT,field,lT))
    for field in vectorfields:
        os.system('cp postProcessing/sample/%d/outlet_field/vectorField/%s constant/boundaryData/inlet/%d/'%(lT,field,lT))

    # os.system('decomposePar -latestTime -fields > log.decompose%d'%(sim))
    # changeInControlDict(nIt+lT)
    if sim == 0:
        for field in scalarfields:
            with open('%d/%s'%(lT,field), 'r') as fl:
                lineBC = fl.readlines()
            indInlet = -1
            indZavorka = -1
            uz = False
            for i in range(len(lineBC)):
                if 'inlet' in lineBC[i]:
                    indInlet = i
                    uz = True
                if uz and '}' in lineBC[i]:
                    indZavorka = i
                    uz = False
            newBC = []
            print('inlet at', indInlet,indZavorka)
            for i in range(indInlet+2):
                newBC.append(lineBC[i])
            newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset 0;\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
            for i in range(indZavorka,len(lineBC)):
                newBC.append(lineBC[i])
            with open('%d/%s'%(lT,field), 'w') as fl:
                for i in range(len(newBC)):
                    fl.writelines(newBC[i])
        
        for field in vectorfields:
            with open('%d/%s'%(lT,field), 'r') as fl:
                lineBC = fl.readlines()
            indInlet = -1
            indZavorka = -1
            uz = False
            for i in range(len(lineBC)):
                if ('inlet' in lineBC[i]) and not uz:
                    indInlet = i
                    uz = True
                if uz and '}' in lineBC[i]:
                    indZavorka = i
                    break
            print('inlet at', indInlet,indZavorka)
            newBC = []
            for i in range(indInlet+2):
                newBC.append(lineBC[i])
            newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset (0 0 0);\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
            for i in range(indZavorka,len(lineBC)):
                newBC.append(lineBC[i])
            with open('%d/%s'%(lT,field), 'w') as fl:
                for i in range(len(newBC)):
                    fl.writelines(newBC[i])
    # os.system('')
    os.system('decomposePar -latestTime -fields > log.decompose%d'%(sim))
    changeInControlDict(nIt+lT)
    # change controldict