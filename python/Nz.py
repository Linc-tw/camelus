import numpy as np
import sys
import os
import matplotlib.pyplot as plt



def ComputeNzFromHalo(filename, dz):
    halcat = np.loadtxt(filename)
    zs = halcat[:,3]
    Ngals = halcat[:,7]
    bin_edges = np.arange(0,np.max(zs),dz)
    counts = []
    for lf, rf in zip(bin_edges, bin_edges[1:]):
        idx = np.where((zs > lf) & (zs <=rf))
        counts += [np.sum(Ngals[idx])]
    counts = np.array(counts)
    return np.array([counts, bin_edges])
    
    

def ComputeNz(filename, dz):
    galcat = np.loadtxt(filename)
    zs = galcat[:,2]
    bin_edges = np.arange(0,np.max(zs),dz)
    Nz = np.histogram(zs, bin_edges)
    return Nz
    
#%%
def main():
    """Read Camelus-generated halo catalog and populate them following a HOD. Syntax:
    
    > python Nz.py path/to/catalogfolders/ halcat galcat dz
    
    Where halcat and galcat are the stubs of the filenames for halo and galaxy
    catalogs to be read, respectively, and dz is the redshift bin width.
    
    """
    cats_dir, halcat_name, galcat_name = sys.argv[1], sys.argv[2], sys.argv[3]
    dz = float(sys.argv[4])
    stub_length = len(halcat_name)
    halcat_files = [catfile for catfile in os.listdir(cats_dir) 
                    if halcat_name == catfile[:stub_length]]
    stub_length = len(galcat_name)
    galcat_files = [catfile for catfile in os.listdir(cats_dir) 
                    if galcat_name == catfile[:stub_length]]
    Nzs_hal = np.array([ComputeNzFromHalo(cats_dir+filename, dz) 
                    for filename in halcat_files])
    Nzs_gal = np.array([ComputeNz(cats_dir+filename, dz) 
                    for filename in galcat_files])
    plt.errorbar(Nzs_hal[0,1][1:]-dz/2, np.mean(Nzs_hal[:,0]), np.std(Nzs_hal[:,0]),
                 label='HOD population')
    plt.errorbar(Nzs_gal[0,1][1:]-dz/2, np.mean(Nzs_gal[:,0]), np.std(Nzs_gal[:,0]),
                 label='Camelus random population')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$N(z)$')
    plt.legend()
    #plt.show()
    plt.savefig(cats_dir+'Nz_plot.png')
    
    
if __name__ == "__main__":
    main()
