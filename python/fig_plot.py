### git clone 
 
from __future__ import print_function
import scipy.integrate as integrate
from astropy.io import ascii
import sys
import math
import numpy as np
from scipy import interpolate
import pylab as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

import scipy.integrate as integrate
import scipy.special as special
from scipy.interpolate import interp1d
import matplotlib.image as mpimg
import matplotlib.colors as colors

from matplotlib.mlab import griddata

plt.rc('text',usetex=True)
plt.rc('font', family='serif', size=12, serif='cz00')

def nz_multi(catHalo,nhalo):
	plt.close('all')
	plt.figure(1)
	nnz2=10
	zmin=0.0001 #min(z)
	zmax=2.
	zz2=np.linspace(zmin,zmax,nnz2)
	nz2=np.linspace(zmin,zmax,nnz2)*0
	dz=zz2[1]-zz2[0]

	for iii in range(nhalo):
		print('test {0}'.format(iii))
		catHalo2="{0}{1:02d}".format(catHalo,iii+1)
		dat = ascii.read(catHalo2)
		z =dat['col4']
		print("z : {0} {1} ".format(min(z),max(z)))
		Ngal_c =dat['col6']
		Ngal_s =dat['col7']
		Ntot=Ngal_c*(1.0+Ngal_s) #=ngc*(1+ngs)
		print("tot gal : {0} ".format(sum(Ntot)))
		for j in range(len(z)):
			ii=int(((z[j]-zmin)/dz)-0.5)
			nz2[ii]=nz2[ii]+Ntot[j]		
		nz2=nz2/(sum(nz2))
		plt.plot(zz2,nz2,'+')
		plt.plot(zz2,nz2,color='b')

	#plt.hist(z,bins=10,normed=True,color='b',alpha=0.5)
	zz=np.linspace(0.01,2+0.1,100)
	nn2=nnz(zz)
	nn2=nn2/(sum(nn2))
	print("tot gal : {0} ".format(sum(nn2)))
	plt.plot(zz,nn2,color='r')
	plt.show()
	return

def nz(catHalo):
	plt.close('all')
	dat = ascii.read(catHalo)
	z =dat['col4']
	print("z : {0} {1} ".format(min(z),max(z)))
	Ngal_c =dat['col6']
	Ngal_s =dat['col7']
	Ntot=Ngal_c*(1.0+Ngal_s) #=ngc*(1+ngs)
	print("tot gal : {0} ".format(sum(Ntot)))
	nnz2=30
	zmin=min(z)
	zz2=np.linspace(zmin,max(z),nnz2)
	nz2=np.linspace(zmin,max(z),nnz2)*0
	dz=zz2[1]-zz2[0]
	for j in range(len(z)):
		ii=int(((z[j]-zmin)/dz)-0.5)
		nz2[ii]=nz2[ii]+Ntot[j]		

	plt.figure(1)
	#nz2=nz2 #/(sum(nz2))
	#plt.hist(z,bins=10,normed=True,color='b',alpha=0.5)
	#zz=np.linspace(0.01,max(z)+0.1,len(zz2))
	#nn2=nnz(zz)
	#nn2=nn2*12.*180.*180. #/(sum(nn2))
	dat = ascii.read('../build/g1')
	z1 =dat['col3']
	nz1=np.linspace(zmin,max(z),nnz2)*0
	for j in range(len(z1)):
		ii=int(((z1[j]-zmin)/dz)-0.5)
		nz1[ii]=nz1[ii]+1	
	plt.semilogy(zz2,nz1,color='r')
	plt.semilogy(zz2,nz2,color='b')
	plt.show()
	return	
	
def nzg(catHalo):
	plt.close('all')
	dat = ascii.read(catHalo)
	z =dat['col3']
	plt.hist(z,bins=10,normed=True,color='b',alpha=0.5)
	zz=np.linspace(0.0001,max(z),10)
	nn2=nnz(zz)
	#nn2=nn2/(sum(nn2))
	plt.plot(zz,nn2,color='r')*10.*6.
	#plt.plot(zz2,nz2,'+')
	#plt.plot(zz2,nz2,color='b')
	plt.show()
	return	

def nnz(z):
	z0=0.5
	alpha=2.
	beta=1.
	x=z/z0
	return  x**alpha * np.exp(-x**beta)

##############
def mass_function_halo(z):
		fich='../build/massFct_z{0:.3f}'.format(z)
		dat = ascii.read(fich)
		m =dat['col1']
		dn=dat['col2']
		plt.close('all')
		plt.figure(1)
		plt.xlabel('log(M) $M_{odot}/h$')
		plt.ylabel(' dn/dlogM   $Mpc/h ^{-3} $')	
		plt.semilogy(np.log10(m),dn)
		plt.title(' Halo mass function ')	
		return

def z_function_halo(fich,zbin):
		dat = ascii.read(fich)
		x= dat['col1']
		y= dat['col2']
		z= dat['col3']
		zz= dat['col4']
		mass= dat['col5']

		plt.close('all')

		zmin=min(zz)
		zmax=max(zz)
		dz=(zmax-zmin)/float(N)
		for i in range(len(zz)):
			iz = int((z-zmin)/dz)
			

		plt.figure(1)
		plt.xlabel('log(M) $M_{odot}/h$')
		plt.ylabel(' dn/dlogM   $Mpc/h ^{-3} $')	
		plt.semilogy(np.log10(m),dn)
		plt.title(' Halo mass function ')	
		return

def histogram_bias(fich,fichbias,N):
	plt.close('all')

	
	if (N==1):
		dat = ascii.read(fich)
		xmin =dat['col1']
		xmax=dat['col2']
		nn=dat['col3']
		nbins = xmax

		plt.figure(1)
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.step(xmin,nn,where='post',color='b',alpha=1,label='no bias')
		plt.title('Peak abundance histogram ')

		dat = ascii.read(fichbias)
		xmin =dat['col1']
		xmax=dat['col2']
		nn=dat['col3']
		nbins = xmax
		plt.step(xmin,nn,where='post',color='r',alpha=1,label='bias')
		plt.legend()
		plt.show()
		return

	if(N!=1):
		fich2=fich+str(1)

		dat = np.loadtxt(fich2)
		xmin =dat[:,0]
		xmax=dat[:,1]
		nn=dat[:,2]

		nbins = xmax
		mean_snr = np.copy(nn)
		mean_snr_error = np.copy(nn)*0
		
		for ii in range(2,N):
			fich2=fich+str(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			nbins = xmax
			mean_snr[:]=mean_snr[:]+nn[:]
		
		mean_snr[:]=mean_snr[:]/float(N)

		for ii in range(1,N):
			fich2=fich+str(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			mean_snr_error = mean_snr_error +(nn-mean_snr)**2.

		mean_snr_error =(1./(N-1.)*mean_snr_error)**(0.5)
		nbin_snr=9
		data=np.linspace(1,5,nbin_snr)

		plt.errorbar(data[:(nbin_snr-1)]+0.25, mean_snr[:(nbin_snr-1)], yerr=mean_snr_error[:(nbin_snr-1)],fmt='+',color='crimson',alpha=0.7)
		plt.step(data,mean_snr,where='post',color='crimson',alpha=1,label="NoBias")


		fich2=fichbias+str(1)
		dat = np.loadtxt(fich2)
		xmin =dat[:,0]
		xmax=dat[:,1]
		nn=dat[:,2]
		nbins = xmax
		mean_snr = np.copy(nn)
		mean_snr_error = np.copy(nn)*0
		
		for ii in range(2,N):
			fich2=fichbias+str(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			nbins = xmax
			mean_snr=mean_snr+nn
		mean_snr=mean_snr/float(N)

		for ii in range(1,N):
			fich2=fichbias+str(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			mean_snr_error = mean_snr_error +(nn-mean_snr)**2.

		mean_snr_error =(1./(N-1.)*mean_snr_error)**(0.5)
		nbin_snr=9
		data=np.linspace(1,5,nbin_snr)

		plt.errorbar(data[:(nbin_snr-1)]+0.25, mean_snr[:(nbin_snr-1)], yerr=mean_snr_error[:(nbin_snr-1)],fmt='+',color='deepskyblue',alpha=0.7)
		plt.step(data,mean_snr,where='post',color='deepskyblue',alpha=1,label="Bias")
		plt.title('Peak abundance histogram (averaged over {0} realizations)'.format(N))
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.legend()
		plt.ion()
		return


##################################

###################################
def histogram_snr(N):
	plt.close('all')

	if (N==1):
		fich='../build/peakHist'
		dat = ascii.read(fich)
		xmin =dat['col1']
		xmax=dat['col2']
		nn=dat['col3']
		nbins = xmax

		plt.figure(1)
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.step(xmin,nn,where='post',color='b',alpha=1)
		plt.title('Peak abundance histogram (only one realization)')
		plt.show()
	if(N!=1):
		fich='../build/dataMat_nbFilters2_N'+str(N)
		dat = np.loadtxt(fich)
		nbin_snr=9
		mean_dat=np.ones(nbin_snr)
		err_dat=np.ones(nbin_snr)*0.
		snr=np.linspace(1,5,nbin_snr)
		for i in range(nbin_snr):
			mean_dat[i]=np.mean(dat[:,i])
			for j in range(N):
				err_dat[i]= err_dat[i]+(dat[j,i]-mean_dat[i])**2.
				if(i==1): plt.step(snr,dat[j,0:nbin_snr],where='post',color='crimson',alpha=0.05)
			err_dat[i]=(1./((np.size(dat[:,i])-1.))*err_dat[i])**(1./2.)

		plt.figure(1)
		plt.errorbar(snr[:(nbin_snr-1)]+0.25, mean_dat[:(nbin_snr-1)],  yerr=err_dat[:(nbin_snr-1)],fmt='+',color='crimson',alpha=0.7)
		plt.step(snr,mean_dat,where='post',color='crimson',alpha=1)
		plt.title('Peak abundance histogram (averaged over {0} realizations)'.format(N))
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.show()
	return


##################################

def map_lensing(fich,opt):
	plt.close('all')
	dat = np.loadtxt(fich)
	if(opt == 1):
		x = dat[:,0]
		y = dat[:,1]
		z = dat[:,2]
		N = int(len(z)**.5)
		x1=np.linspace(0.5,179.5,N)
		y1=np.linspace(0.5,179.5,N)
		X1,Y1 = np.meshgrid(x1,y1)  
  		z2=np.reshape(z,(N,N))

		plt.figure(1)
		plt.pcolor(X1,Y1,(z2),cmap=cm.terrain,vmin=z2.min(),vmax=z2.max())    
		cbar =plt.colorbar()
		plt.title( "Convergence map" )
		plt.xlabel(r' $\theta_x$ (arcmin)')
		plt.ylabel(r' $\theta_y$ (arcmin)')

	if(opt == 2):
		x = dat[:,0]
		y = dat[:,1]
		g1 = dat[:,2]
		g2 = dat[:,3]
		N = int(len(x)**.5)
		x1=np.linspace(0.5,179.5,N)
		y1=np.linspace(0.5,179.5,N)
		X1,Y1 = np.meshgrid(x1,y1)  
  		z1=np.reshape(g1,(N,N))
  		z2=np.reshape(g2,(N,N))

		plt.figure(1)
		plt.subplot(1,2,1)
		plt.pcolor(X1,Y1,(z1),cmap=cm.OrRd,vmin=z1.min(),vmax=z1.max())    
		cbar =plt.colorbar()
		plt.title( '$\\rm \gamma_1$ Map')
		plt.xlabel(r' $\theta_x$ (arcmin)')
		plt.ylabel(r' $\theta_y$ (arcmin)')

		plt.subplot(1,2,2)
		plt.pcolor(X1,Y1,(z2),cmap=cm.OrRd,vmin=z2.min(),vmax=z2.max())    
		cbar =plt.colorbar()
		plt.title( '$\\rm \gamma_2$ Map')
		plt.xlabel(r' $\theta_x$ (arcmin)')
		plt.ylabel(r' $\theta_y$ (arcmin)')

	plt.show()
	return;
#####################################

def map_peak(fich,fich2,lim):
	
	dat = np.loadtxt(fich2)
	snr = dat[:,0]
	x = dat[:,1]
	y = dat[:,2]
	xx=x[snr>lim]
	yy=y[snr>lim]

	dat = np.loadtxt(fich)
	x2 = dat[:,0]
	y2 = dat[:,1]
	z = dat[:,2]

	plt.close('all')

	N = int(len(z)**.5)
	z2=np.reshape(z,(N,N))
	x1=np.linspace(min(x2),max(x2),N)
	y1=np.linspace(min(y2),max(y2),N)
	X1,Y1 = np.meshgrid(x1,y1)     
    
	plt.figure(9)
	plt.pcolor(X1,Y1,(z2),cmap=cm.terrain,vmin=z2.min(),vmax=z2.max())    
	cbar =plt.colorbar()
	cbar.set_label(r'$\kappa$', rotation=-270)
 	plt.scatter(xx,yy, marker='o',s=200, edgecolor='red', alpha=1,zorder=10,facecolor='None')
	plt.title( 'Convergence Map and peak detection for SNR $>$ {0:.1f}'.format(lim))
	plt.xlabel(r' $\theta_x$ (arcmin)')
	plt.ylabel(r' $\theta_y$ (arcmin)')
	plt.axis([X1.min(), X1.max(), Y1.min(),Y1.max()]) 
	plt.show()

	return;

##########################################

def read_catalogue_galaxies(fich):
	dat = np.loadtxt(fich)
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	red = dat[:,2]
	k  = dat[:,3]
	g1 = dat[:,4]
	g2 = dat[:,5]
	N=np.size(k)

	plt.close('all')
	plt.figure(1)

	
	d1 =np.reshape(red, (1,np.product(red.shape)))[0]
	plt.hist(d1,bins=100,color = 'limegreen', edgecolor = 'green',alpha=0.60, normed=1)
	plt.yscale('log')
	plt.title( 'Redshift histogramm of galaxy sources')
	plt.ylabel(r'\LARGE{$\rm N_{sources}$}')
	plt.xlabel(r'Redshift')

	return ;
########################################

def read_catalog_halo(fich,opt):
	dat = np.loadtxt(fich)
	plt.close('all')
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	ww = dat[:,2]
	zz = dat[:,3]
	mm = dat[:,4]
	plt.figure(1)

	if(opt==1):
		MIN=min(mm)
		MAX=max(mm)
		nbins =50
		nbins = np.linspace(np.log10(MIN), np.log10(MAX), nbins)
		plt.title( 'Halo mass histogramm')
		res =plt.hist(np.log10(mm), bins = nbins,color = 'limegreen', edgecolor = 'green',alpha=0.60, stacked=True)
		plt.yscale('log')
		ylab=r'{ $ \rm N_{halo} $}'
		xlab=r'{ $ \rm log( M_{halo} ( \rm M_{\odot}/h ) ) $}'
		axes = plt.gca()
		plt.xlabel(xlab)
		plt.ylabel(ylab)
	if(opt==2):
		res =plt.hist(zz, bins = 30,color = 'limegreen', edgecolor = 'green',alpha=0.60, stacked=True)
		plt.yscale('log')
		plt.title( 'Halo redshift histogramm')	
		ylab=r' $ \rm N_{halo} $'
		xlab=r'redshift'
		plt.xlabel(xlab)
		plt.ylabel(ylab)
	return ;
######################

def rewrite_Galaxy_catalogue(fich2,fich3):
	dat = np.loadtxt(fich2)
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	red = dat[:,2]
	k  = dat[:,3]
	g1 = dat[:,4]
	g2 = dat[:,5]
	nn=np.size(k)
	m=0 #0.1 ##bias used routine   m= -0.00856*d  ## champ densite   k mass/surface 
	g1n =np.copy(g1)
	g2n =np.copy(g2)
	file3 = open(fich3,"w") 
	file2 = open(fich2,"r") 
	aa=0
	for i in range(nn):
		x = theta_x[i]
		for j in range(nn):
			y = theta_y[j]


	for line in file2: 
		file3.write(line)
		aa=aa+1
		if(aa==12): break

	file2.close()

	for i in range(nn):
		file3.write('\t {0:2.3f} \t {1:3.3f} \t {2:4.5f} \t {3:4.5f} \t {4:4.5f} \t {5:4.5f} \n'.format(dat[i,0],dat[i,1],dat[i,2],dat[i,3],g1[i],g2[i]))
  	file3.close()
	return ;

##########################
######################

def redefine_Galaxy_catalogue(fich2,fich3,zz,dz):
	datg = np.loadtxt(fich2)
	xg = datg[:,0]
	yg = datg[:,1]
	zg = datg[:,2]
	k  = datg[:,3]
	g1 = datg[:,4]
	g2 = datg[:,5]
	ng=np.size(k)

	dath = np.loadtxt(fich3)
	xh = dath[:,0]
	yh = dath[:,1]
	wh = dath[:,2]
	zh = dath[:,3]
	mh = dath[:,4]
	nh=np.size(zh)

	xh2=xh[((zh>zz-dz)&(zh<zz+dz))]
	yh2=yh[((zh>zz-dz)&(zh<zz+dz))]

	xg2=xg[((zg>zz-dz)&(zg<zz+dz))]
	yg2=yg[((zg>zz-dz)&(zg<zz+dz))]

	print("\nnb galaxies : {0}, density nb/arcmin2 n={1} \n".format(ng,ng/(180.**2)))
	print("nb halo : {0}, density nb/arcmin2 n={1} \n".format(nh,nh/(180.**2)))

	plt.close('all')

	plt.figure(1)
	plt.xlim([-10,190])
	plt.ylim([-10,190])
	plt.plot(xh2,yh2,'.',color='crimson',alpha=0.2)
	plt.show()

	

	plt.figure(2)
	plt.xlim([-10,190])
	plt.ylim([-10,190])
	plt.plot(xg2,yg2,'.',color='deepskyblue',alpha=0.2)
	plt.show()
	return

