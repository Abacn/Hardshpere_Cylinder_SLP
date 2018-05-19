import iomodule, randomconfig, lp

if __name__ == '__main__':
    parameter = iomodule.readParameter()
    lpobj = lp.LPCore(parameter)
    iomodule.writeCoords(lpobj.coords, lpobj.boxNow, "initialConfig.dat")
    lpobj.optimize()
    iomodule.writeCoords(lpobj.coords, lpobj.boxNow, "resultConfig.dat")
