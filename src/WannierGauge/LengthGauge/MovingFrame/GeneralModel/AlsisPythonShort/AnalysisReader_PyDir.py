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

#####################################################
FigureDir 	        = '/' + 'Figure'
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
InterC 	        = np.transpose( np.loadtxt( IntraInterPath ) )
IntraC          = np.transpose( np.loadtxt( IntraInterPath1 ) )
Laser 	        = np.transpose( np.loadtxt( LaserPath ) )

# memory order align
InterC = np.ascontiguousarray(InterC, dtype=np.float64)
IntraC = np.ascontiguousarray(IntraC, dtype=np.float64)

with io.open(ProjPath + '/inputParam.cfg') as f:
  config = libconf.load(f)

Npulses = len(config['laser']['pulses'])
print("\n\n==============================")
print(f"Total {Npulses} pulses") 
for i in range(Npulses):
  print('------------------------------')
  print("Pulse ", i+1)
  print(fr"E0 = {config['laser']['pulses'][i]['E0']} a.u. , w0 = {config['laser']['pulses'][i]['w0']} a.u. ,"
    fr"ellip = {config['laser']['pulses'][i]['ellip']}, ncycles = {config['laser']['pulses'][i]['ncycles']}, "
    fr"cep = {config['laser']['pulses'][i]['cep']} rad")
  print(fr"t0 = {config['laser']['pulses'][i]['t0']} a.u. , "
    fr"phix = {config['laser']['pulses'][i]['phix']} rad, thetaz = {config['laser']['pulses'][i]['thetaz']} rad, "
    fr"phiz = {config['laser']['pulses'][i]['phiz']} rad, envelope = {config['laser']['pulses'][i]['env_name']}")
print("==============================")

#####################################################

#####
# paramters changing for each plots
# width and height is applied for 'every' plots including current
width = 11
hight = width/1.62

axesSym = ['x', 'y', 'z']
axesColor = ['b', 'g', 'r']
#axesColor = {'x':'red', 'y':'blue', 'z':'green', 'total':'black', 'inter':'blue', 'intra':'green', 'rcp':'red', 'lcp':'blue'}
# these are for plotting spectrums
hzero           = .98e-25

refPulse = 0
w0		= config['laser']['pulses'][refPulse]['w0']
T0		= 2.*np.pi/w0

spectrum_xmin = 0
spectrum_xmax = 37
spectrum_ymin = np.log10(hzero)
spectrum_ymax = 1
xticks0  = np.arange(1,50,4)
yticks0  = np.arange(-20,6,5)

name_laser_offset = ("w0__" + str("%.3f"%config['laser']['pulses'][refPulse]['w0']) + "__E0__" + str("%.3f"%config['laser']['pulses'][refPulse]['E0']) 
  + "__e0__" + str("%.3f"%config['laser']['pulses'][refPulse]['ellip']))

name_material_offset = config['target']

name_offset = name_material_offset + '__' + name_laser_offset + "__" + nparam

print("Reference pulse = ", refPulse)
print("Period, T0  = ", T0, " a.u.")
print("Mean-freq   = ", w0, " a.u.")

Dim = 2

#####################################################
#Time axis and electric field 
t		    = Laser[0, :]
Efield = Laser[refPulse*4 + 1:refPulse*4 + 3, :]

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
plt.plot( t, Efield[0, :], 'b', lw=2 ,label="$E_{x'}$")
plt.plot( t, Efield[1, :], 'g', lw=2 ,label="$E_{y'}$")

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
#Ploting the current oscillations
#width = 11
#hight = width/1.62


fig = plt.figure(figsize=(width,hight) )


plt.plot( t, Efield[0, :]/max(Efield[0, :])*np.amax(InterC), 'k', lw=2 , label="$E_{x'}$")
for i in range(Dim):
  plt.plot( t, InterC[i, :], axesColor[i], lw=1.5, label="$J_{"+axesSym[i]+",er}$")


#plt.legend([p1, p2, p3], ['$E_{L}$', '$J_{x,er}$', '$J_{y,er}$']) 
plt.legend()

plt.xlabel('time (a.u.) ', fontsize=18)
plt.ylabel('$J_{x,y,z}$ (a.u.) ', fontsize=18)
plt.tick_params(labelsize=18)

xaxmin      = t.min()     
xaxmax      = t.max()    
plt.xlim( xaxmin, xaxmax )

plt.tight_layout()

filename0           = FigureDir + '/CurrentOscillations.pdf' 
fileNamePicture     = ProjPath + filename0  
plt.savefig( fileNamePicture )

#############################################################
#filtering dipole and current oscillations by means of 
#applying a smoth time mask over the beginning and end of pulses, this will avoid 
#high esporeous frequencies...

InterCMasked = np.zeros(InterC.shape, dtype=np.float64)
IntraCMasked = np.zeros(IntraC.shape, dtype=np.float64)
for i in range(Dim):
  InterCMasked[i, :] = np.blackman(Nt) * InterC[i, :]
  IntraCMasked[i, :] = np.blackman(Nt) * IntraC[i, :]

#####################################################
fig = plt.figure(figsize=(width,hight) )


plt.plot( t, Efield[0, :]/max(Efield[0, :])*np.amax(IntraCMasked), 'red', lw=2 , label=r'$E_{x}$')
for i in range(Dim):
  plt.plot( t, IntraCMasked[i, :], axesColor[i], lw=1.5, label='$J_{'+axesSym[i]+',ra}$')
#plt.title(title_name,fontsize=18);
#plt.legend([p1, p2, p3], ['$E_{L}$', '$J_{x,ra}$', '$J_{y,ra}$'])
plt.legend()

plt.xlabel('time (a.u.) ', fontsize=18)
plt.ylabel('$J_{x,y,z}$ (a.u.) ', fontsize=18)
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

for i in range(Dim):
  plt.plot( t, InterCMasked[i, :], axesColor[i], lw=1, label='$J_{'+axesSym[i]+',er}$')

#plt.legend([p2, p3], ['$J_{x,er}$', '$J_{y,er}$'], fontsize = 18 ) 

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

FFT_InterC = np.zeros(InterC.shape, dtype=np.complex128)
FFT_IntraC = np.zeros(InterC.shape, dtype=np.complex128)

FullRadiation = np.zeros(InterC.shape, dtype=np.complex128)
FullRadiationRotated = np.zeros(InterC.shape, dtype=np.complex128)

for i in range(Dim):
  FFT_InterC[i, :] = fft( InterCMasked[i, :] ) * dt
  FFT_IntraC[i, :] = fft( IntraCMasked[i, :] ) * dt

  FFT_InterC[i, :] = -1j*w*np.fft.fftshift( FFT_InterC[i, :] )
  FFT_IntraC[i, :] = -1j*w*np.fft.fftshift( FFT_IntraC[i, :] )

  FullRadiation[i, :] = FFT_InterC[i, :] + FFT_IntraC[i, :]

Sinter = FFT_InterC[0, :] + FFT_InterC[1, :] #+ FFT_InterC[2, :]
Sintra = FFT_IntraC[0, :] + FFT_IntraC[1, :] #+ FFT_IntraC[2, :]
# rotate Rz(phiz)Ry(thetaz)
thetaz = config['laser']['pulses'][refPulse]['thetaz']
phiz = config['laser']['pulses'][refPulse]['phiz']
FullRadiationRotated[0, :] = (np.cos(thetaz)*np.cos(phiz)    * FullRadiation[0, :]
                              + np.cos(thetaz)*np.sin(phiz)  * FullRadiation[1, :] )
                              #- np.sin(thetaz)               * FullRadiation[2, :] )
FullRadiationRotated[1, :] = (-np.sin(phiz)    * FullRadiation[0, :]
                              + np.cos(phiz)  * FullRadiation[1, :] )
# FullRadiationRotated[2, :] = (np.sin(thetaz)*np.cos(phiz)    * FullRadiation[0, :]
#                               + np.sin(thetaz)*np.sin(phiz)  * FullRadiation[1, :]
#                               + np.cos(thetaz)               * FullRadiation[2, :] )

# RCP, LCP
dJ_p            = FullRadiationRotated[0, :] + (1j)*FullRadiationRotated[1, :]
dJ_m            = FullRadiationRotated[0, :] - (1j)*FullRadiationRotated[1, :]


############################
############################

Spectrum = np.log10( abs(FullRadiation[0, :])**2 +  abs(FullRadiation[1, :])**2 + hzero )#+ abs(FullRadiation[2, :])**2 + hzero)
axisSpectrum = np.zeros(FullRadiation.shape, dtype=np.float64)
for i in range(Dim):
  axisSpectrum[i, :] = np.log10( abs(FullRadiation[i, :])**2 + hzero )

a_dJ_p          = np.log10( abs( dJ_p )**2 + hzero ) 
a_dJ_m          = np.log10( abs( dJ_m )**2 + hzero ) 
#####################################################


print("\n\n+=+++++++++++++++++++++=+\nSpectra shape       = ",  Spectrum.shape)
print("Frequency axis shape = ", w.shape, "\n")

Nthalf          = int(np.floor(Nt/2)) 
w_output        = w[Nthalf:Nt-1]/w0

spectrum_output = pow( 10., Spectrum[Nthalf:Nt-1] )


OutputData      = w_output[:, None] 
for i in range(Dim):
  OutputData = np.concatenate( (OutputData, np.real(FullRadiation[i, Nthalf:Nt-1])[:, None] ), axis=1 )
  OutputData = np.concatenate( (OutputData, np.imag(FullRadiation[i, Nthalf:Nt-1])[:, None] ), axis=1 )
OutputData      = np.concatenate( (OutputData, spectrum_output[:, None] ), axis=1 ) 

ofname          = "/HHG__" + name_offset + '.dat'

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

for i in range(Dim):
  plt.plot(w/w0, axisSpectrum[i], axesColor[i], lw=1.5, label='$J_{'+axesSym[i]+'}$' )
plt.plot( w/w0, Spectrum, 'k',  lw = 1.2, label='$J_{tot}$' )

plt.legend(fontsize=18)
#plt.legend([p1, p2], ['$J_{x}$', '$J_{y}$'],fontsize=18)

plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 )
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 )
plt.tick_params(labelsize = 28 )

plt.xticks( xticks0 )
plt.yticks( yticks0 )

plt.grid(True)

plt.xlim(spectrum_xmin, spectrum_xmax)
plt.ylim(spectrum_ymin, spectrum_ymax)

fname='/LinearHarmonicSpectrum__'+ name_offset + '.pdf'


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

fname='/InterIntraHarmonicSpectrum__'+ name_offset + '.pdf'


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


fname='/CircularHarmonicSpectrum__' + name_offset + '.pdf'


filename1           = FigureDir + fname 
fileNamePicture     = ProjPath  + filename1  


plt.savefig(fileNamePicture, dpi = 300) 
shutil.copyfile( fileNamePicture, BasicPath + set_DataName + fname ) 

if (args.isShow):
  plt.show()

#####################################################
