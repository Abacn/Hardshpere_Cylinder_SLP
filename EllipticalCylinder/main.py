from __future__ import print_function

import sys

import iomodule
import lp

if __name__ == '__main__':
    # default file name of initial and final configuration
    default_initfile = "initialConfig_%d.dat"
    default_finalfile = "resultConfig_%d.dat"
    # default parameter file
    default_parafile = "config.json"
    
    # parse in parameters
    if 1 == len(sys.argv):
        # do not have parse in parameters, run one sample
        repeatrun = 1
    else:
        try:
            repeatrun = int(sys.argv[1])
        except ValueError:
            print('Usage: main.py [repeatrun]')
            exit()
    # read in parameters
    parameter = iomodule.readParameter()
    for i in range(repeatrun):
        lpobj = lp.LPCore(parameter)
        iomodule.writeCoords(lpobj.coords, lpobj.boxNow, lpobj.diam, default_initfile%(i))
        # print parameters
        for k in parameter:
            if not k.startswith('_'):
                print(k, ': ', parameter[k])
        print("Optimization begin")
        status = lpobj.optimize()
        print("Optimization ended code ", status)
        iomodule.writeCoords(lpobj.coords, lpobj.boxNow, lpobj.diam, default_finalfile%(i))
