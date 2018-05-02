import numpy as np
import sys
import os
import matplotlib.pyplot as plt



def ComputeNzFromHalo(filename, dz):
    halcat = np.loadtxt(filename)
    zs = halcat[:,3]
    Ngals = halcat[:,5] * (1+halcat[:,6])
    bin_edges = np.arange(0,np.max(zs)+dz,dz)
    counts = []
    for lf, rf in zip(bin_edges, bin_edges[1:]):
        idx = np.where((zs > lf) & (zs <=rf))
        counts += [np.sum(Ngals[idx])]
    counts = np.array(counts)
    return np.array([counts, bin_edges])

def CamelusNz(z, alpha=2., beta=1., z_0=.5):
    x = z/z_0
    #return z**2 / (2*z_0**3) * np.exp(-x**beta)
    return x**alpha * np.exp(-x**beta)

def ComputeNz(filename, dz):
    galcat = np.loadtxt(filename)
    zs = galcat[:,2]
    bin_edges = np.arange(0,np.max(zs)+dz,dz)
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
    print 'Average number of generated galaxies (HOD):\t{}'.format(np.mean([np.sum(count) for count in Nzs_hal[:,0]]))
    Nzs_gal = np.array([ComputeNz(cats_dir+filename, dz) 
                    for filename in galcat_files])
    print 'Average number of generated galaxies (galaxy catalog):\t{}'.format(np.mean([np.sum(count) for count in Nzs_gal[:,0]]))
    np.save(cats_dir+'Nzs_hal.npy', Nzs_hal)
    np.save(cats_dir+'Nzs_gal.npy', Nzs_gal)
    zs = Nzs_hal[0,1][1:]-dz/2
    HOD_n = np.mean(Nzs_hal[:,0])
    HOD_n /= np.sum(HOD_n)
    plt.errorbar(zs, HOD_n, np.std(HOD_n),
                 label='HOD population')
    CamZs = [CamelusNz(zee) for zee in zs]
    CamZs /= np.sum(CamZs)
    plt.plot(zs, CamZs, label='Camelus n(z)')
    gal_n = np.mean(Nzs_gal[:,0]).astype('float')
    gal_n /= np.sum(gal_n)
    plt.errorbar(zs, gal_n, np.std(gal_n),
                 label='Camelus random population')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$n(z)$')
    plt.legend(loc=5, bbox_to_anchor=(1.6,.5))
    #plt.show()
    plt.savefig(cats_dir+'Nz_plot.png', bbox_inches='tight')
    
    
if __name__ == "__main__":
    main()
