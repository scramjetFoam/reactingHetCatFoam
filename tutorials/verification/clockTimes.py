from os import walk
from os.path import isfile


verCases = [dir for dir in next(walk('.'))[1] if 'baseCase' not in dir and 'inner' in dir]
# verCases = [dir for dir in next(walk('.'))[1] if 'baseCase' not in dir and 'flow' in dir]
# verCases = [dir for dir in next(walk('.'))[1] if 'baseCase' not in dir]
# print(verCases)
for verCase in verCases:
    ZZ_dir = verCase+'/ZZ_cases'
    solver = 'reactingHetCatSimpleFoam'
    if 'Mass' in verCase: solver += 'M'

    caseDirs = next(walk(ZZ_dir))[1]
    # print(verCase + ' -> '+str(caseDirs[0])+', ...')
    times = []
    timessum = 0
    timesnum = 0
    for caseDir in caseDirs:
        if caseDir != 'protoMesh':
            caseDir = ZZ_dir+'/'+caseDir+'/'
            log1 = caseDir+'log1.'+solver
            log2 = caseDir+'log2.'+solver
            time = 'not found'

            if isfile(log2): log = log2
            else: log = log1

            with open(log, 'r') as logfile:
                lines = logfile.readlines()
            time = int(lines[-5].split(' ')[-2])
            times.append(time)
            timessum += time
            timesnum += 1

    print('clocktimes for '+verCase+'\n\t'+str(times))
    print('average clocktime for '+verCase+' (in s):\n\t'+str(timessum/timesnum))