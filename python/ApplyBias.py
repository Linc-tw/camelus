import numpy as np
import sys
import linecache
import os

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

def ApplyBias(galcat, delta, n, b=-0.00856):
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
    """Read mean density from Camelus-generated galaxy catalog.

        Parameters
        ----------
        galcat_path : str
            Path to galaxy catalog.

        Returns
        -------
        n_mean : float
            Average density in catalog.
    """
    line = linecache.getline(galcat_path,3)
    return np.float(line.split()[3])
    
def ConvertCats(galcat_dir, filename, savestub):
    """Read, compute local densities and apply bias to given galaxy catalog.

        Parameters
        ----------
        galcat_dir : str
            Path to folder containing galaxy catalogs. Note the biased catalogs
            will also be saved there.
            
        filename : str
            Name of the file containing the galaxy catalog.
            
        savestub : str
            What biased galaxy catalogs should be called. Will be followed by an
            underscore (_) and the up-to-3 digit(s) id of the read catalog where
            applicable.
                
    """
    print 'Applying bias to galaxy catalog {}.'.format(filename)
    galcat = np.loadtxt(galcat_dir+filename)
    # read average density
    nmean = read_n_mean(galcat_dir+filename)
    # compute local densities from galaxy catalogs
    delta = ComputeDensity(galcat)
    # apply bias
    galcat_b = ApplyBias(galcat, delta, nmean)
    # save biased galaxy catalog
    idnb = '_'
    for digit in [char for char in filename[-3:] if char.isdigit()]:
        idnb += digit
    np.savetxt(galcat_dir+savestub+idnb, galcat_b)
    
def main():
    """Apply bias to Camelus-generated galaxy catalogs. Syntax:
    
    > python ApplyBias path/to/catalogfolders/ galcat bgalcat
    
    Where galcat_names and biased_galcat are the stubs of the filenames for
    catalogs to be read and (biased) ones to be written, e.g. with the run
    above, if files galcat_000, galcat_001 and galcat_002 are present in
    catalogfolders, ApplyBias will generate and save bgalcat_000, bgalcat_001
    and bgalcat_002 to catalogfolders.
    
    """
    galcat_dir, galcat_name = sys.argv[1], sys.argv[2]
    stub_length = len(galcat_name)
    galcat_files = [catfile for catfile in os.listdir(galcat_dir) 
                    if galcat_name == catfile[:stub_length]]
    _ = [ConvertCats(galcat_dir, filename, sys.argv[3]) for filename in galcat_files]
    
if __name__ == "__main__":
    main()
