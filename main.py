import iomodule, randomconfig, lp

if __name__ == '__main__':
    parameter = iomodule.readParameter()
    lpobj = lp.LPCore(parameter)
    iomodule.writeCoords(lpobj.coords, lpobj.boxNow, "initialConfig.dat")
    status = lpobj.optimize()
    print("Optimization ended code ", status)
    iomodule.writeCoords(lpobj.coords, lpobj.boxNow, "resultConfig.dat")
