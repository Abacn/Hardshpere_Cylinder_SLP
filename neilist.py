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
        self.resetIter()

    def resetIter(self):
        """Reset iterator"""
        self.iter_n = 0
        self.iter_idx = 0

    def __init__(self, parameter, coords=None):
        self.parameter = parameter
        # parameters for neighborlist is optional
        if "NNLextraDist" not in parameter:
            self.parameter["NNLextraDist"] = 0.1
        self.sqrNNLextraDist = self.parameter["NNLextraDist"]*self.parameter["NNLextraDist"]
        if "NNLupdFreq" not in parameter:
            self.parameter["NNLupdFreq"] = 20
        self.updInterval = self.parameter["NNLupdFreq"]
        if coords is not None:
            self.genNeiList(coords)

    def __iter__(self):
        return self

    def __next__(self):
        """Iteration method"""
        while self.iter_idx >= self.neighNum[self.iter_n]:
            # meet the end of neilist and skip the blank neilist of such particle
            self.iter_idx = 0
            self.iter_n += 1
            if self.iter_n >= self.N:
                raise StopIteration
        pairs =(self.iter_n, self.neighIdx[self.iter_n][self.iter_idx])
        self.iter_idx += 1
        return pairs