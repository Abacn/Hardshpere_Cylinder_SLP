import glob, os, shutil
from packingfraction import getPackingFraction

def densestconfig(folder, nsubfolder=0):
    """Get the filename of the densest configuration"""
    namestr = folder + '/*'*nsubfolder + '/resultConfig_*.dat'
    filenames = glob.glob(namestr)
    maxeta = 0.  # maximal packing fraction
    maxfname = None
    for fn in filenames:
        eta = getPackingFraction(fn)
        if maxeta < eta:
            maxeta = eta
            maxfname = fn
    return maxfname
    
if __name__ == '__main__':
    folders = filter(os.path.isdir, glob.glob('*/*'))
    for fdname in folders:
        maxfname = densestconfig(fdname, nsubfolder=1)
        cpname = os.path.normpath(fdname + '/resultConfid_max.dat')
        if maxfname is not None and cpname != os.path.normpath(maxfname): # robust fix
            print(maxfname)
            shutil.copyfile(maxfname, cpname)
            
        