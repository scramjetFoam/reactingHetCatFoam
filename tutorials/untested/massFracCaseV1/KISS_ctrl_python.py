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
    lstIntT = [int(lstInt1) for lstInt1 in lstInt] 
    return (max((lstIntT)))

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

def changeK0(value):
    with open('constant/reactiveProperties','r') as fl:
        lines = fl.readlines()
    for i in range(len(lines)):
        if 'thisline' in lines[i]:
            lines[i] = '\t\tk0  k0  [0 0 -1 0 0 0 0]   %g; // thisline\n'%value
    with open('constant/reactiveProperties','w') as fl:
        for i in range(len(lines)):
            fl.writelines(lines[i])     

def p(kde,value):
    with open('%s/p'%kde,'r') as fl:
        lines = fl.readlines()
    for i in range(len(lines)):
        if 'p0' in lines[i]:
            lines[i] = '\tp0\t\tuniform %.12g;\n'%value
        if 'value' in lines[i]:
            lines[i] = '\tvalue\t\tuniform %.12g;\n'%value
    with open('%s/p'%kde,'w') as fl:
        for i in range(len(lines)):
            fl.writelines(lines[i])   



N = 5
Nzacatek = 0
simDir = 'rashigsV10'
simDir = 'rashigsV12'
simDir = 'rashigsV13'
simDir = 'rashigsV14'
simDir = 'rashigsV15'
simDir = 'rashigsV16'
simDir = 'rashigsV17'
simDir = 'rashigsV18'
simDir = 'rashigsV19'
simDir = 'rashigsV20'
simDir = 'rashigsV21'
simDir = 'rashigsV22'
# simDir = 'reformaxV1'
# simDir = 'reformaxV2'
# simDir = 'vagWheelV1'
# simDir = 'vagWheelV2'
# simDir = 'testMappedBC'
simDir = './'
# simDir = '../testMappedBC'
scalarfields = ['T','ethylene','prod','N2']
vectorfields = ['U']
fields = "'(T ethylene U p prod N2)'"
# nIt = 850
nIt = 850

whenChange = [0,7,22, 26]
values = [55000*1.1,75000*1.1,91700*1.1]
values = [55000*1.1,75000*1.1,91700*1.1,100000*1.1]

# whenChange = [0,7,13,20,26]
# values = [55000,69623,93694,94400,221000]
# values = [55000*1.5,69623*1.5,75000*1.5,94400*1.5,221000*1.5]


os.chdir(simDir)
pOut = 501325

if Nzacatek > 0:
    lT = latestTime()
    changeInControlDict(nIt+lT)
else:
    changeInControlDict(nIt)

for sim in range(Nzacatek,N):
    print('Running sim. %d %g -- %g m'%(sim,(sim)*0.1,(sim+1)*0.1))
    for j in range(len(whenChange)):
        if whenChange[j] == sim:
            changeK0(values[j])
    os.system('foamJob -parallel -screen reactingHetCatSimpleFoamFC > log.rash%d'%sim)
    os.system('reconstructPar -latestTime -fields %s > log.reconstruct%d'%(fields,sim))
    os.system("postProcess -func 'sample' > log.sample%d"%sim)
    lT = latestTime()
    print('Result in %d folder'%lT)
    if sim == 0:
        os.system('mkdir constant/boundaryData')
        os.system('mkdir constant/boundaryData/inlet')
    
    os.system('mkdir constant/boundaryData/inlet/%d'%lT)
        # os.system('cp postProcessing/sample/%d/outlet_field/points constant/boundaryData/inlet/%d/'%(lT,lT))
    os.system('cp postProcessing/sample/%d/outlet_field/points constant/boundaryData/inlet/'%(lT))

    
    os.system('rm -f constant/boundaryData/inlet/%d/*'%(lT))
    for field in scalarfields:
        os.system('cp postProcessing/sample/%d/outlet_field/scalarField/%s constant/boundaryData/inlet/%d/'%(lT,field,lT))

    os.system("postProcess -func 'patchAverage(p,name=inlet,patch=inlet)' > log.postProcessP%d"%sim)
    with open ('postProcessing/patchAverage(p,name=inlet,patch=inlet)/0/surfaceFieldValue.dat', 'r') as fl:
        lines = fl.readlines()
    print(lines[-1].split("\t")[-1].replace('\n',''))
    dp = float(lines[-1].split("\t")[-1].replace('\n',''))-pOut
    # print('nastavuji nove pOut %g'%(pOut-dp)) 
    # p(lT,pOut-dp)
    # pOut = pOut-dp

    for field in vectorfields:
        os.system('cp postProcessing/sample/%d/outlet_field/vectorField/%s constant/boundaryData/inlet/%d/'%(lT,field,lT))

    # os.system('decomposePar -latestTime -fields > log.decompose%d'%(sim))
    # changeInControlDict(nIt+lT)
    # if sim == 0 and not pozdejsiSim:
    for field in scalarfields:
        with open('%d/%s'%(lT,field), 'r') as fl:
            lineBC = fl.readlines()
        indInlet = -1
        indZavorka = -1
        uz = 0
        for i in range(len(lineBC)):
            if 'inlet' in lineBC[i]:
                indInlet = i
                uz = 1
            if uz == 1 and '}' in lineBC[i]:
                if sim == 0:
                    indZavorka = i
                    uz = 0
                else:
                    uz = 2
            elif uz == 2 and '}' in lineBC[i]:
                if sim != 0:
                    indZavorka = i
                uz = 0
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
        # with open('%d/%s_n'%(lT,field), 'w') as fl:
        #     for i in range(len(newBC)):
        #         fl.writelines(newBC[i])
            

    for field in vectorfields:
        with open('%d/%s'%(lT,field), 'r') as fl:
            lineBC = fl.readlines()
        indInlet = -1
        indZavorka = -1
        uz = 0
        # for i in range(len(lineBC)):
        #     if ('inlet' in lineBC[i]) and not uz:
        #         indInlet = i
        #         uz = True
        #     if uz and '}' in lineBC[i]:
        #         indZavorka = i
        #         break
        for i in range(len(lineBC)):
            if 'inlet' in lineBC[i] and uz == 0:
                indInlet = i
                uz = 1
            if uz == 1 and '}' in lineBC[i]:
                if sim == 0:
                    indZavorka = i
                    uz = 3
                else:
                    uz = 2
            elif uz == 2 and '}' in lineBC[i]:
                if sim != 0:
                    indZavorka = i
                uz = 3
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

    for field in scalarfields:
        # repaire here
        os.system("postProcess -func 'patchAverage(%s,name=outlet,patch=outlet)' > log.postProcess%s%d"%(field,field,sim))
        with open ('postProcessing/patchAverage(%s,name=outlet,patch=outlet)/0/surfaceFieldValue.dat'%(field), 'r') as fl:
            lines = fl.readlines()
            fCorr = float(lines[-1].split("\t")[-1].replace('\n',''))
            # print(fCorr)

        os.system("postProcess -func 'patchAverage(%s,name=inlet,patch=inlet)' > log.postProcess_o%s%d"%(field,field,sim))
        with open ('postProcessing/patchAverage(%s,name=inlet,patch=inlet)/0/surfaceFieldValue.dat'%(field), 'r') as fl:
            lines = fl.readlines()
            # print(lines[-1].split("\t")[-1].replace('\n',''))
            fsouc = float(lines[-1].split("\t")[-1].replace('\n',''))
        print(fCorr,fsouc,fCorr/fsouc) 

        with open ('constant/boundaryData/inlet/%d/%s'%(lT,field) , 'r') as fl:
            lines = fl.readlines()
        for i in range(3,len(lines)-1):
            # f = float(lines[i].replace('\n',''))
            f = float(lines[i].replace('\n',''))*float(fCorr)/float(fsouc)
            lines[i] = '%.15g\n'%(f)
        with open ('constant/boundaryData/inlet/%d/%s'%(lT,field) , 'w') as fl:
            for i in range(len(lines)):
                fl.writelines(lines[i])
        os.system("postProcess -func 'patchAverage(%s,name=inlet,patch=inlet)' > log.postProcess_n%s%d"%(field,field,sim))

    # repaire
    os.system("postProcess -func 'patchAverage(U,name=outlet,patch=outlet)' > log.postProcessU%d"%sim)
    with open ('postProcessing/patchAverage(U,name=outlet,patch=outlet)/0/surfaceFieldValue.dat', 'r') as fl:
        lines = fl.readlines()
    Ucorr = (lines[-1].split("\t")[-1].replace('\n','').replace('(','').replace(')','').split(' '))
    # print(Ucorr)

    os.system("postProcess -func 'patchAverage(U,name=inlet,patch=inlet)' > log.postProcessU%d"%sim)
    with open ('postProcessing/patchAverage(U,name=inlet,patch=inlet)/0/surfaceFieldValue.dat', 'r') as fl:
        lines = fl.readlines()
    # print(lines[-1].split("\t")[-1].replace('\n',''))
    Usouc = (lines[-1].split("\t")[-1].replace('\n','').replace('(','').replace(')','').split(' '))
    print(Ucorr,Usouc,float(Ucorr[0])/float(Usouc[0]),float(Ucorr[1])/float(Usouc[1]),float(Ucorr[2])/float(Usouc[2])) 
    
    with open ('constant/boundaryData/inlet/%d/U'%(lT) , 'r') as fl:
        lines = fl.readlines()
    for i in range(3,len(lines)-1):
        Ux = float(lines[i].split(' ')[0].replace('(',''))*float(Ucorr[0])/float(Usouc[0])#*pOut/(pOut+dp)
        Uy = float(lines[i].split(' ')[1])*float(Ucorr[1])/float(Usouc[1])*pOut/(pOut+dp)*0
        Uz = float(lines[i].split(' ')[2].replace(')\n',''))*float(Ucorr[2])/float(Usouc[2])*pOut/(pOut+dp)*0
        # Ux = float(lines[i].split(' ')[0].replace('(',''))#*float(Ucorr[0])/float(Usouc[0])#*pOut/(pOut+dp)
        # Uy = float(lines[i].split(' ')[1])
        # Uz = float(lines[i].split(' ')[2].replace(')\n',''))
        lines[i] = '(%.15g %.15g %.15g)\n'%(Ux,Uy, Uz)
    with open ('constant/boundaryData/inlet/%d/U'%(lT) , 'w') as fl:
        for i in range(len(lines)):
            fl.writelines(lines[i])
    # os.system('')
    os.system('decomposePar -latestTime -fields > log.decompose%d'%(sim))
    changeInControlDict(nIt+lT)
    # change controldict
