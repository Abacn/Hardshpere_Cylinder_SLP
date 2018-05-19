import numpy as np
from tools import pbc_dz

class NeiList:
    """Neighbor List"""
    def genNeiList(self, coords):
        """Generate neighbor list"""
        self.N = coords.shape[0]
        self.neighNum = np.zeros(self.N)
        self.neighIdx = [list([]) for i in range(self.N)]
        for i in range(self.N):
            for j in range(i+1, self.N):
                dz = coords[i, 2] - coords[j, 2]
                dz = pbc_dz(dz)
                if abs(np.abs(dz) < self.sqrNNLextraDist):
                    self.neighNum[i] += 1
                    self.neighIdx[i].append(j)

    def __init__(self, parameter, coords=None):
        self.parameter = parameter
        # parameters for neighborlist is optional
        if "NNLextraDist" not in parameter:
            self.parameter["NNLextraDist"] = 0.1
        self.sqrNNLextraDist = (1+self.parameter["NNLextraDist"])*(1+self.parameter["NNLextraDist"])
        if "NNLupdFreq" not in parameter:
            self.parameter["NNLupdFreq"] = 20
        self.updInterval = self.parameter["NNLupdFreq"]
        if coords is not None:
            self.genNeiList(coords)

