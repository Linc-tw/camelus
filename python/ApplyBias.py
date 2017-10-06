import numpy as np
import matplotlib.pyplot as plt
import sys
import linecache
import os

# set path to Camelus outputs

def ComputeDensity(galcat):
    """Compute density delta(x,y) from a Camelus-generated galaxy catalog.

        Parameters
        ----------
        galcat : np.ndarray
            Camelus-generated galaxy catalog as numpy array.

        Returns
        -------
        delta : tuple
            Density in numpy.histogram2d convention, that is, delta[0][i,j] is the 
            density for galaxies at position x=delta[1][i], y=delta[2][j].
    """
    x_bins = range(int(np.min(galcat[:,0])),int(np.max(galcat[:,0]))+1)
    y_bins = range(int(np.min(galcat[:,1])),int(np.max(galcat[:,1]))+1)
    delta = np.histogram2d(galcat[:,0], galcat[:,1], bins=[x_bins,y_bins])
    return delta

def ApplyBias(galcat, delta, n, b= 0.00856):
    """Compute and apply multiplicative bias to shear.

        Parameters
        ----------
        galcat : np.ndarray
            Camelus-generated galaxy catalog as numpy array.
            
        delta : tuple
            Density in numpy.histogram2d convention.
            
        a : float
            Slope for the linear relationship between bias and density.
            Default: value from Hoekstra et al., 2016, Figure 5/Section 5.2.
            
        b : float
            Intercept for the linear relationship between bias and density.
            Default: value from Hoekstra et al., 2016, Figure 5/Section 5.2.
        
        Returns
        -------
        galcat_b : np.ndarray
            Galaxy catalog with multiplicative bias applied to shear measurements.
    """
    a=-b/n
    galcat_b = np.copy(galcat)
    xs, ys = galcat[:,0].astype(int), galcat[:,1].astype(int)
    # add maximum position to last bin
    xs[xs == np.max(delta[1])] = np.max(delta[1])-1
    ys[ys == np.max(delta[2])] = np.max(delta[2])-1
    delta_pos = np.array([delta[0][np.where(delta[1]==x)[0],np.where(delta[2]==y)[0]] 
                            for x,y in zip(xs,ys)])
    m = a + b*delta_pos
    galcat_b[:,-2:] *= (1+m)
    return galcat_b


def read_n_mean(galcat_path) :
    line = linecache.getline(galcat_path,3)
    return np.float(line.split()[3])
    
    

    
def main():
    galcat_dir, galcat_name = sys.argv[1], sys.argv[2]
    galcat_files = [catfile for catfile in os.listdir(galcat_dir) 
                    if galcat_name in catfile]
    galcats = np.array([np.loadtxt(galcat_dir+catfile) for catfile in galcat_files])
    nmeans = np.array([read_n_mean(galcat_dir+catfile) for catfile in galcat_files])
    # compute density from galaxy catalogs
    deltas = np.array([ComputeDensity(galcat) for galcat in galcats])
    # apply bias
    galcats_b = np.array([ApplyBias(galcat, delta, nmean) for galcat, delta, nmean in zip(galcats, deltas, nmeans)])
    # save biased galaxy catalogs
    for galcat_b, filename in zip(galcats_b, galcat_files):
        idnb = '_'
        for digit in [char for char in filename[-3:] if char.isdigit()]:
            idnb += digit
        np.savetxt(galcat_dir+sys.argv[3]+idnb, galcat_b)
    
if __name__ == "__main__":
    main()
