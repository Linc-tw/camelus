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

def CamelusNz(z, alpha=2., beta=1., z_0=.5):
    x = z/z_0
    return x**alpha * np.exp(-x**beta)

def CutOff(fullzs, bin_edges, nobj, dz=None):
    """apply cut off
    """
    if bin_edges is None:
        bin_edges = np.arange(0,np.max(fullzs)+dz,dz)
    select_idx = []
    zmid = [rf-(rf-lf)/2 for lf, rf in zip(bin_edges, bin_edges[1:])]
    nz = np.array([CamelusNz(z) for z in zmid])
    nz = (nz * nobj/np.sum(nz)).astype(int)
    for nb_rand, (lf, rf) in zip(nz,zip(bin_edges, bin_edges[1:])):
        print '   > Working on z bin [{},{}]'.format(lf,rf)
        zidx = np.where((fullzs>lf) & (fullzs<=rf))[0]
        nb_hod = len(zidx)
        if nb_hod < nb_rand:
            print '/!\/!\ Not enough sources for this redshift bin!\t{} v {} for HOD and random cats. /!\/!\ '.format(nb_hod,nb_rand)
            nb_rand = nb_hod
        this_selection = np.random.choice(zidx, nb_rand, False)
        select_idx += list(this_selection)
    print '   > Total number of sources : {}'.format(len(select_idx))
    return select_idx
    
def ConvertCats(galcat_dir, filename, bin_edges, savestub, savestub_b, nobj, dz=None):
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
    print ' > Applying bias to galaxy catalog {}.'.format(filename)
    galcat = np.loadtxt(galcat_dir+filename)
    # compute local densities from galaxy catalogs
    delta = ComputeDensity(galcat)
    # compute average density
    nmean = np.mean(delta[0])
    # apply bias
    galcat_b = ApplyBias(galcat, delta, nmean)
    print ' > Bias applied.'
    # apply cutoff
    print ' > Applying cutoff to galaxy catalogs {}.'.format(filename)
    select_idx = CutOff(galcat[:,2], bin_edges, nobj, dz)
    cutoff = galcat[select_idx,:]
    cutoff_b = galcat_b[select_idx,:]
    print ' > Cutoff performed.'.format(filename)
    # save cut off and biased galaxy catalog
    idnb = '_'
    for digit in [char for char in filename[-3:] if char.isdigit()]:
        idnb += digit
    np.savetxt(galcat_dir+savestub+idnb, cutoff)
    np.savetxt(galcat_dir+savestub_b+idnb, cutoff_b)
    
def main():
    """Apply bias to Camelus-generated galaxy catalogs. Syntax:
    
    > python cutoff.py path/to/catalogfolders/ galcat nobj dz output_cutoff output_cutoff_b [zmin] [zmax]
    
    """
    galcat_dir, galcat_name, nobj = sys.argv[1], sys.argv[2], float(sys.argv[3])
    dz = float(sys.argv[4])
    savestub, savestub_b = sys.argv[5], sys.argv[6]
    if len(sys.argv) > 7:
        zmin = float(sys.argv[7])
        zmax = float(sys.argv[8])
        bin_edges = np.arange(zmin,zmax+dz,dz)
    else:
        bin_edges = None
    stub_length = len(galcat_name)
    galcat_files = [catfile for catfile in os.listdir(galcat_dir) 
                    if galcat_name == catfile[:stub_length]]
    _ = [ConvertCats(galcat_dir, filename, bin_edges, savestub, savestub_b, nobj, dz) for 
         filename in galcat_files]
    
if __name__ == "__main__":
    main()
