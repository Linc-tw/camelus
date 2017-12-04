import numpy as np
import sys
import os
import matplotlib.pyplot as plt


def PopulateGalaxiesFromHalo(filename):
    halocat = np.genfromtxt(filename, names=True, skip_header = 10)
    First = True
    for halo in halocat :
        if First :
            gal = PopulateOneHalo(halo)
            First = False
        else :
            gal = np.vstack((gal,PopulateOneHalo(halo)))
            

    return gal

def PopulateOneHalo(halo) :

    First = True
    if np.random.random() < halo['Ngal_c'] :
        gal = np.array([halo['theta_x'], halo['theta_y'], halo['z']])
        First = False
        
    for i in xrange(np.int(halo['Ngal_s']+0.5)) :
        R_false = True
        while R_false :
            xtest = np.random.random()
            if NFW(5*xtest)>np.random.random() :
                r = halo['Rv'] * xtest
                R_false = False
        theta = np.random.random() * 2 * np.pi
        phi = np.arccos(np.random.random() * 2 - 1)
        x = np.cos(theta) * np.sin(phi) * r
        y = np.sin(theta) * np.sin(phi) * r
        if First :
            gal = np.array([halo['theta_x'] + x, halo['theta_y'] + y, halo['z']])
        else :
            gal = np.vstack((gal,[halo['theta_x'] + x, halo['theta_y'] + y, halo['z']]))

    return gal


def NFW(x) :
    return 1/(x*(1+np.power(x,2)))
    



def main():

    cats_dir, halcat_name, galcat_name = sys.argv[1], sys.argv[2], sys.argv[3]
    stub_length = len(halcat_name)
    halcat_files = [catfile for catfile in os.listdir(cats_dir) 
                    if halcat_name == catfile[:stub_length]]
    for halcat in halcat_files :
        galcat = PopulateGalaxiesFromHalo(cats_dir+halcat)
        idnb = '_'
        for digit in [char for char in halcat[-3:] if char.isdigit()]:
            idnb += digit
            np.savetxt(cats_dir + galcat_name + idnb, galcat)
    
    
if __name__ == "__main__":
   main()
