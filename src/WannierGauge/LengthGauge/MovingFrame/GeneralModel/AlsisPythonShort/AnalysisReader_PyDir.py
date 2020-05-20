#!/usr/bin/env python3
"""
Reader and visual. of harmonic emission from solids, Chacon model """
import matplotlib.pyplot as plt 
from matplotlib.colors import BoundaryNorm 
from matplotlib.ticker import MaxNLocator 
import numpy as np
import os 
import sys
import argparse
import shutil 
import cmath
import errno
from scipy.fftpack import fft, ifft
import io, libconf



##########################################
#Functions 
def masking_dipole(t,f,ta,tb,asigma,bsigma):
  Nt=len(t) 
  amask=f 
  for j in range(0,Nt):
    if (t[j] <= ta):
      amask[j] = f[j]*np.exp(- ((t[j]-ta)/asigma)**2 ) 
      
    if (t[j] >= tb):
      amask[j] = f[j]*np.exp(- ((t[j]-tb)/bsigma)**2 ) 
      
  return amask 


parser = argparse.ArgumentParser(description='Analysis HHG or HSG outputs.')
parser.add_argument('datapath', nargs='?', default='', type=str,help='path to folder contains data and config file')
parser.add_argument('nparam', nargs='?', default='', type=str,help='additional name to dat file')
parser.add_argument('--noshow', dest='isShow', action='store_const', const=False, default=True, help = 'Do not show the plot (default: show)')
parser.add_argument('--f', dest='forceOverwrite', action='store_const', const=True, default=False, help = 'Force overwrite (default: ask)')

args = parser.parse_args()

mparam = 'xy'
lparam = args.datapath
nparam = args.nparam

#Creating path and file name for reading files 
BasicPath	        = os.getcwd()

FolderName              = "/" + lparam
ProjPath  	        = BasicPath + FolderName

FileNameIntraInter      = '/interband_dipole_full_evol.dat'
FileNameIntraInter1     = '/intraband_current_full_evol.dat'
FileNameLaser 	        = '/outlaserdata.dat'
set_DataName 	        = '/SetData0'

FileArgSave = 'ListAnalysis.txt'

#####################################################
FigureDir 	        = '/' + mparam + 'Figure'
IntraInterPath 	        = ProjPath  + FileNameIntraInter
IntraInterPath1         = ProjPath  + FileNameIntraInter1
LaserPath 	        = ProjPath  + FileNameLaser



try:
  os.mkdir( ProjPath + FigureDir )
except OSError as exc:
  if (exc.errno != errno.EEXIST):
    raise exc  
  pass  

try:
    os.mkdir( BasicPath + set_DataName )
except OSError as exc:
    if (exc.errno != errno.EEXIST):
        raise exc
    pass


#########################
# LOADING DATA # 
InterC 	        = np.loadtxt( IntraInterPath )
IntraC          = np.loadtxt( IntraInterPath1 )
Laser 	        = np.loadtxt( LaserPath )


with io.open(ProjPath + '/inputParam.cfg') as f:
  config = libconf.load(f)

Npulses = len(config['laser']['pulses'])
print("\n\n+=+++++++++++++++++++++++++=+\nLaser features:") 
for i in range(Npulses):
  print(("Pulse ", i))
  print("E0        = ", config['laser']['pulses'][i]['E0'], " a.u.") 
  print("w0        = ", config['laser']['pulses'][i]['w0'], " a.u.") 
  print("N.O.C.    = ", config['laser']['pulses'][i]['ncycles'], " ") 
  print("Ellip     = ", config['laser']['pulses'][i]['ellip'], " ") 



#####################################################


#####
# paramters changing for each plots
# width and height is applied for 'every' plots including current
width = 11
hight = width/1.62

axesSym = ['x', 'y', 'z']
axesColor = ['red', 'green', 'blue']
# these are for plotting spectrums
hzero           = .98e-25

refPulse = 0
w0		= config['laser']['pulses'][i]['w0']
T0		= 2.*np.pi/w0

spectrum_xmin = 0
spectrum_xmax = 37
spectrum_ymin = np.log10(hzero)
spectrum_ymax = 1
xticks0  = np.arange(1,50,4)
yticks0  = np.arange(-20,6,5)

print("Reference pulse = ", refPulse)
print("Period, T0  = ", T0, " a.u.")
print("Mean-freq   = ", w0, " a.u.")

#####################################################
#Time axis and electric field 
t		    = Laser[:,0]
Efield      = []
for i in range(3):
  Efield.append(Laser[:,i+1])
Enorm = np.sqrt(Efield[0]**2 + Efield[1]**2 + Efield[2]**2)

Nt          = len(t)
dt 	        = t[1]-t[0]


print("\n\n\n+++++++++++++++++++++++++++++")
print("Time-step dt             = ", dt, " a.u. ") 
print("Total No. of time-steps  = ", Nt, " ") 
print('+=+++++++++++++++++++++++=+\n\n')

#####################################################
#Frequency axis 
wmax 	    	= np.pi/dt
dw		        = 2.*wmax/float(Nt)
wmin  	    	= -wmax
wmax 	    	= +wmax
w           	= np.linspace( wmin, wmax-dw, num = Nt )


#####################################################
if (Nt != len(w)):
    print("\nsNtime have to be equal to length (w)\n")
    exit()
#####################################################


print("\nNomega-Shape   = ", w.shape)
print("dw             = ", dw)
print('wmax           = {0:.2f}'.format(wmax)) 
print('+=++++++++++++++++++++++++=+\n')


#####################################################
#width = 11
#hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
for i in range(3):
  plt.plot( t, Efield[i], axesColor[i], lw=2 ,label='E'+axesSym[i])
plt.plot(t, Enorm, 'black', lw=2, label='|E|')

plt.xlabel('time (a.u.) ', fontsize=18)
plt.ylabel('Efield (a.u.) ', fontsize=18)
plt.legend()
plt.tick_params(labelsize=18)

xaxmin      = t.min()     #
xaxmax      = t.max()     # 

plt.xlim(xaxmin, xaxmax)
plt.tight_layout()

filename0           = FigureDir + '/LaserField.png'
fileNamePicture     = ProjPath  + filename0
plt.savefig( fileNamePicture ) 



#####################################################
print("\n\nLen of interC = ", InterC.shape)


#####################################################
#Momentum Integration for the "remains" of MPI code
xinterC 	= InterC[:,0] - InterC[0,0]
yinterC     = InterC[:,1] - InterC[0,1]

xJintra  	= IntraC[:,0] - IntraC[0,0]
yJintra     = IntraC[:,1] - IntraC[0,1]


#####################################################
#Computing inter-current contribution from dipole polarization InterC1
xJinter		   = np.zeros( xinterC.shape, np.float )
yJinter        = np.zeros( yinterC.shape, np.float )



xJinter 	   = xinterC    #np.diff( xinterC )/dt;
yJinter        = yinterC    #np.diff( yinterC )/dt;



#####################################################
#Ploting the current oscillations
#width = 11
#hight = width/1.62


fig = plt.figure(figsize=(width,hight) )


p1, = plt.plot( t, Enorm/max(Enorm)*max(xJinter), 'red', lw=2 )
p2, = plt.plot( t, xJinter, axesColor[0], lw = 1.5 )
p3, = plt.plot( t, yJinter, axesColor[1],  lw = 1.5 )



plt.legend([p1, p2, p3], ['$E_{L}$', '$J_{x,er}$', '$J_{y,er}$']) 



plt.xlabel('time (a.u.) ', fontsize=18)
plt.ylabel('Jx/Jy (a.u.) ', fontsize=18)
plt.tick_params(labelsize=18)



xaxmin      = t.min()     
xaxmax      = t.max()     #
plt.xlim( xaxmin, xaxmax )


plt.tight_layout()




filename0           = FigureDir + '/CurrentOscillations.pdf' 
fileNamePicture     = ProjPath + filename0  
plt.savefig( fileNamePicture )


#############################################################
## Computing harmonic spectra, inter and intra contribution
i 	= cmath.sqrt(-1)#complex(0,1);
print("complex number, imag. base = ", i)


#############################################################
#filtering dipole and current oscillations by means of 
#applying a smoth time mask over the beginning and end of pulses, this will avoid 
#high esporeous frequencies...

xJinterMasked = np.blackman(Nt)*xJinter
yJinterMasked = np.blackman(Nt)*yJinter

xJintraMasked = np.blackman(Nt)*xJintra
yJintraMasked = np.blackman(Nt)*yJintra


#####################################################
fig = plt.figure(figsize=(width,hight) )


p1, = plt.plot( t, Enorm/max(Enorm)*max(xJintra), 'red', lw=2 )
p2, = plt.plot( t, xJintraMasked, 'green', lw = 1.5 )
p3, = plt.plot( t, yJintraMasked, 'blue',  lw = 1.5 )
#plt.title(title_name,fontsize=18);



plt.legend([p1, p2, p3], ['$E_{L}$', '$J_{x,ra}$', '$J_{y,ra}$'])



plt.xlabel('time (a.u.) ', fontsize=18)
plt.ylabel('Jx/Jy (a.u.) ', fontsize=18)
plt.tick_params(labelsize=18)



xaxmin      = t.min()     #
xaxmax      = t.max()     #
plt.xlim( xaxmin, xaxmax )


plt.tight_layout()




filename0           = FigureDir + '/IntraCurrentOscillations.pdf' 
fileNamePicture     = ProjPath + filename0 
plt.savefig( fileNamePicture ) 



#####################################################
#Ploting the current oscillations
#width   = 11
#hight   = width/1.62
fig     = plt.figure( figsize=(width,hight) )

p2,= plt.plot( t, xJinterMasked, 'green', lw = 1 ) 
p3,= plt.plot( t, yJinterMasked, 'blue',  lw = 1 ) 
plt.legend([p2, p3], ['$J_{x,er}$', '$J_{y,er}$'], fontsize = 18 ) 

plt.xlabel( 'time (a.u.) ', fontsize = 18 ) 
plt.ylabel( 'Filter Currents-Mask (a.u.) ', fontsize = 18 ) 
plt.tick_params( labelsize = 18 ) 

xaxmin      = t.min() 
xaxmax      = t.max() 

plt.xlim(xaxmin, xaxmax) 
plt.tight_layout() 


#########################################
###   Saving picture   ###
filename0           = FigureDir + '/CurrentOscillationsMF.pdf' 
fileNamePicture     = ProjPath + filename0  
plt.savefig( fileNamePicture ) 



#####################################################
##Calculating FFT 
xFFT_Jinter 	= fft( xJinterMasked )*dt 
yFFT_Jinter 	= fft( yJinterMasked )*dt 

xFFT_Jintra     = fft( xJintraMasked )*dt 
yFFT_Jintra     = fft( yJintraMasked )*dt 

xFFT_Jinter      = -np.fft.fftshift( xFFT_Jinter )*(1j)*w 
yFFT_Jinter      = -np.fft.fftshift( yFFT_Jinter )*(1j)*w 

xFFT_Jintra  	 = -np.fft.fftshift( xFFT_Jintra )*(1j)*w 
yFFT_Jintra  	 = -np.fft.fftshift( yFFT_Jintra )*(1j)*w 

xFullRadiation 	 = xFFT_Jinter   +  xFFT_Jintra 
yFullRadiation   = yFFT_Jinter   +  yFFT_Jintra 

# RCP, LCP
dJ_p            = xFullRadiation - (1j)*yFullRadiation
dJ_m            = xFullRadiation + (1j)*yFullRadiation


############################
############################

xSinter         = np.log10( abs( xFFT_Jinter )**2 + hzero ) 
ySinter         = np.log10( abs( yFFT_Jinter )**2 + hzero ) 
Sinter          = np.log10( abs( xFFT_Jinter )**2 + abs( yFFT_Jinter )**2 + hzero ) 


xSintra         = np.log10( abs( xFFT_Jintra )**2 + hzero ) 
ySintra         = np.log10( abs( yFFT_Jintra )**2 + hzero ) 
Sintra 	        = np.log10( abs( xFFT_Jintra )**2 + abs( yFFT_Jintra )**2 + hzero) 


xSpectrum 	    = np.log10( abs( xFullRadiation )**2 + hzero ) 
ySpectrum       = np.log10( abs( yFullRadiation )**2 + hzero ) 
Spectrum        = np.log10( abs( xFullRadiation )**2 + abs( yFullRadiation )**2 + hzero ) 


a_dJ_p          = np.log10( abs( dJ_p )**2 + hzero ) 
a_dJ_m          = np.log10( abs( dJ_m )**2 + hzero ) 
#####################################################


print("\n\n+=+++++++++++++++++++++=+\nSpectra shape       = ",  xSinter.shape)
print("Frequency axis shape = ", w.shape, "\n")

## Note that this assumes that HSG uses only two pulses...
#if (spectrumType == 'HSG'):
#  w -= LParams[(refPulse+1)%2, 2]


w.shape             = (Nt,1)
xSpectrum.shape     = (Nt,1)
ySpectrum.shape     = (Nt,1)
Spectrum.shape      = (Nt,1)

xFullRadiation.shape = (Nt,1) 
yFullRadiation.shape = (Nt,1) 



Nthalf          = int(np.floor(Nt/2)) 
#temp0           = np.array([Nthalf]) 
#temp0.shape     = (1,1)
#w_output        = np.concatenate( ( temp0, w[Nthalf:Nt-1]/w0 ),axis=0 ) 
w_output        = w[Nthalf:Nt-1]/w0


#temp0           = np.array([Phi0]) 
#temp0.shape     = (1,1)
#real_jxoutput   = np.concatenate( ( temp0, np.real(xFullRadiation[Nthalf:Nt-1]) ),axis=0 ) 
real_jxoutput   = np.real(xFullRadiation[Nthalf:Nt-1]) 


#temp0           = np.array([M0t2]) 
#temp0.shape     = (1,1)
#imag_jxoutput    = np.concatenate( ( temp0, np.imag(xFullRadiation[Nthalf:Nt-1]) ),axis=0 ) 
imag_jxoutput    = np.imag(xFullRadiation[Nthalf:Nt-1]) 



#temp0           = np.array([ChernNo]) 
#temp0.shape     = (1,1)
#real_jyoutput   = np.concatenate( ( temp0, np.real( yFullRadiation[Nthalf:Nt-1]) ),axis=0 ) 
real_jyoutput   = np.real( yFullRadiation[Nthalf:Nt-1]) 


#temp0           = np.array([Eg]) 
#temp0.shape     = (1,1)
#imag_jyoutput   = np.concatenate( ( temp0, np.imag( yFullRadiation[Nthalf:Nt-1] ) ),axis=0 ) 
imag_jyoutput   = np.imag( yFullRadiation[Nthalf:Nt-1] ) 


#temp0           = np.array([dt]) 
#temp0.shape     = (1,1)
#spectrum_output = np.concatenate( ( temp0, pow( 10., Spectrum[Nthalf:Nt-1] ) ),axis=0 )
spectrum_output = pow( 10., Spectrum[Nthalf:Nt-1] )


OutputData      = w_output 
OutputData      = np.concatenate( (OutputData, real_jxoutput   ), axis=1 ) 
OutputData      = np.concatenate( (OutputData, imag_jxoutput   ), axis=1 ) 
OutputData      = np.concatenate( (OutputData, real_jyoutput   ), axis=1 ) 
OutputData      = np.concatenate( (OutputData, imag_jyoutput   ), axis=1 ) 
OutputData      = np.concatenate( (OutputData, spectrum_output ), axis=1 ) 

name_laser_offset = ("w0__" + str("%.3f"%config['laser']['pulses'][refPulse]['w0']) + "__E0__" + str("%.3f"%config['laser']['pulses'][refPulse]['E0']) 
  + "__e0__" + str("%.3f"%config['laser']['pulses'][refPulse]['ellip']) + "__" + nparam  + "__")
ofname          = "/HHG__" + name_laser_offset + '.dat'

out                = ProjPath + ofname
np.savetxt(out, OutputData , fmt='%1.16e')


NewDir = '.' + set_DataName + ofname
if (not os.path.exists(NewDir) or args.forceOverwrite):
  shutil.copyfile( out, NewDir )
else:
  print("Warning file already exists - check whether you make typo")
  inputkey = input("If you want to replace the file and proceed, press y: ")
  if (inputkey == 'y'):
    shutil.copyfile( out, NewDir )
  else:
    sys.exit()

#####################################################
############################
#Ploting the harmonic current radiations-oscillations


fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

p1, = plt.plot( w/w0, xSpectrum, 'r-', lw = 2., label='$J_{x}$' )
p2, = plt.plot( w/w0, ySpectrum, 'b',  lw = 1.5, label='$J_{y}$' )
p3, = plt.plot( w/w0, Spectrum, 'g',  lw = 1.2, label='$J_{t}$' )

plt.legend([p1, p2], ['$J_{x}$', '$J_{y}$'],fontsize=18)

plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 )
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 )
plt.tick_params(labelsize = 28 )


plt.xticks( xticks0 )
plt.yticks( yticks0 )

plt.grid(True)

plt.xlim(spectrum_xmin, spectrum_xmax)
plt.ylim(spectrum_ymin, spectrum_ymax)

fname='/'+ mparam+'LinearHarmonicSpectrum__'+ name_laser_offset + '.pdf'


filename1           = FigureDir + fname
fileNamePicture     = ProjPath  + filename1  


plt.savefig(fileNamePicture, dpi = 300)
shutil.copyfile( fileNamePicture, BasicPath + set_DataName + fname )

#####################################################
############################
#Ploting the harmonic current radiations-oscillations
#width = 11
#hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])


p11, = plt.plot( w[Nthalf:Nt-1]/w0, Sinter[Nthalf:Nt-1], 'b',  lw = 1.2, label='$J_{er}$' ) 
p22, = plt.plot( w[Nthalf:Nt-1]/w0, Sintra[Nthalf:Nt-1], 'r',  lw = 1.2, label='$J_{ra}$' ) 

plt.legend([p11, p22], [ '$J_{er}$', '$J_{ra}$'],fontsize=18)

plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 ) 
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 ) 
plt.tick_params(labelsize = 28 ) 

plt.xticks( xticks0 )
plt.yticks( yticks0 )

plt.grid(True)

plt.xlim(spectrum_xmin, spectrum_xmax)
plt.ylim(spectrum_ymin, spectrum_ymax)

plt.grid(True)

fname='/'+str(mparam)+'InterIntraHarmonicSpectrum__'+ name_laser_offset + '.pdf'


filename1           = FigureDir + fname 
fileNamePicture     = ProjPath  + filename1  


plt.savefig(fileNamePicture, dpi = 300) 
shutil.copyfile( fileNamePicture, BasicPath + set_DataName + fname ) 





#####################################################
############################
#Ploting the harmonic current radiations-oscillations
#width = 11
#hight = width/1.62


fig  = plt.figure(figsize=(width,hight) )
ax1  = fig.add_axes([0.2, 0.15, 0.75, 0.75])
p1,  = plt.plot( w[Nthalf:Nt-1]/w0, a_dJ_p[Nthalf:Nt-1], 'r-', lw = 2., label='$J_{rcp}$' ) 
p2,  = plt.plot( w[Nthalf:Nt-1]/w0, a_dJ_m[Nthalf:Nt-1], 'b',  lw = 1.5, label='$J_{lcp}$' ) 
p3,  = plt.plot( w[Nthalf:Nt-1]/w0, Spectrum[Nthalf:Nt-1], 'g',  lw = 1.2, label='$J_{t}$' ) 

plt.legend([p1, p2,p3], [r'$J_{rcp}$', r'$J_{lcp}$', r'$J_{tot}'],fontsize=18) 

plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 ) 
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 ) 
plt.tick_params(labelsize = 28 )

plt.xticks( xticks0 )
plt.yticks( yticks0 )

plt.grid(True)

plt.xlim(spectrum_xmin, spectrum_xmax)
plt.ylim(spectrum_ymin, spectrum_ymax)


fname='/'+str(mparam)+'CircularHarmonicSpectrum__' + name_laser_offset + '.pdf'


filename1           = FigureDir + fname 
fileNamePicture     = ProjPath  + filename1  


plt.savefig(fileNamePicture, dpi = 300) 
shutil.copyfile( fileNamePicture, BasicPath + set_DataName + fname ) 

if (args.isShow):
  plt.show()

argsave = open(FileArgSave, 'a+')
argsave.write(f'{lparam}  {nparam} \n')
argsave.close()
#####################################################
