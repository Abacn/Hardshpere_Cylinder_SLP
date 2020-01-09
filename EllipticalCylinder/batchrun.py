from __future__ import print_function
from __future__ import division

import sys, os

import iomodule
import lp

# Data folder structure:
# p/lambda/n/resultConfig_0.dat
# Ex: 1/1/16/resultConfig_0.dat


if __name__ == '__main__':
    # default file name of initial and final configuration
    default_initfile = "%s/initialConfig_%d.dat"
    default_finalfile = "%s/resultConfig_%d.dat"
    # default parameter file
    default_parafile = "config.json"
    
    # set parameters
    # phoas = [2, 2.5, 3, 4]
    phoas = [4]
    aratios = [1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.30]
    nparts = range(10, 31)
    # phoas = [1]; aratios = [1]; nparts = [10, 12, 15]
    repeatruns = 10
    runstart = 0
    # read in parameters
    parameter = iomodule.readParameter()
    for p in phoas:
        for ratio in aratios:
            for npart in nparts:
                parameter["NParticle"] = npart
                parameter["aspectRatio"] = ratio
                parameter["z"] = parameter["Ra"]*npart/p
                foldername = '%g/%g/%d' % (p, ratio, npart)
                try:
                    os.makedirs(foldername)
                except OSError:
                    pass
                for i in range(runstart,runstart+repeatruns):
                    lpobj = lp.LPCore(parameter)
                    iomodule.writeCoords(lpobj.coords, lpobj.boxNow, lpobj.diam, default_initfile%(foldername, i))
                    # print parameters
                    for k in parameter:
                        if not k.startswith('_'):
                            print(k, ': ', parameter[k])
                    print("Optimization begin")
                    status = lpobj.optimize()
                    print("Optimization ended code ", status)
                    iomodule.writeCoords(lpobj.coords, lpobj.boxNow, lpobj.diam, default_finalfile%(foldername, i))