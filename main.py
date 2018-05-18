import iomodule, randomconfig, lp

if __name__ == '__main__':
    parameter = iomodule.readParameter()
    lpobj = lp.LPCore(parameter)
    lpobj.optimize()
    rststr = iomodule.writeCoords(lpobj.coords, lpobj.boxNow)
    print(rststr)