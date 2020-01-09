import os, glob
import numpy as np
import matplotlib.pyplot as plt

from packingfraction import getPackingFraction

if __name__ == '__main__':
    rhos = ['1', '2', '2.5', '3', '4']
    saveResults = True  # option for saving results or not 
    results = {}
    for rho in rhos:
        nrho = float(rho)
        data = []
        fnames = glob.glob(rho+'/*/resultConfid_max.dat')
        for fname in fnames:
            aratio = os.path.basename(os.path.dirname(fname))
            nratio = float(aratio)
            eta = getPackingFraction(fname)
            data.append([nratio, eta])
        arr = np.array(data)
        arr.sort(axis=0)
        results[rho] = arr
    plt.figure()
    for rho in rhos:
        plt.plot(results[rho][:,0], results[rho][:,1], label=rho)
    plt.xlabel('$\lambda$')
    plt.ylabel('$\eta$')
    plt.legend()
    plt.show()
    if saveResults:
        with open('packingfractions.txt', 'w') as fout:
            fout.write("rho*a\taspectratio\tpackingFraction\n")
            for rho in rhos:
                for rp in range(results[rho].shape[0]):
                    fout.write("%s\t%f\t%f\n" % (rho, results[rho][rp][0], results[rho][rp][1]))
        