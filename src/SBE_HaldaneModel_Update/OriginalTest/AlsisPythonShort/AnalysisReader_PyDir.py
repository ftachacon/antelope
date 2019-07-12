#!/usr/bin/python
"""
Reader and visual. of harmonic emission from solids, Chacon model """
import matplotlib.pyplot as plt 
from matplotlib.colors import BoundaryNorm 
from matplotlib.ticker import MaxNLocator 
import numpy as np
import os 
import sys 
import shutil 
import cmath
import errno
from scipy.fftpack import fft, ifft
import parameters as cs



##########################################
#Functions 
def masking_dipole(t,f,ta,tb,asigma,bsigma):
  Nt=len(t);
  amask=f;
  for j in range(0,Nt):
    if (t[j] <= ta):
      amask[j] = f[j]*np.exp(- ((t[j]-ta)/asigma)**2 );
      
    if (t[j] >= tb):
      amask[j] = f[j]*np.exp(- ((t[j]-tb)/bsigma)**2 );
      
  return amask;



#Clear
#plt.cla()
#plt.clf()

#Getting parameters
set_of_params   = sys.argv;

iparam 	    	=  int( set_of_params[1] );  #intensity param
kparam          =  int( set_of_params[2] );  # Chern-No.
lparam      	=  set_of_params[3] ;       #dir-name
mparam      	=  set_of_params[4] ;       #xdirection





#Creating path and file name for reading files 
BasicPath	        = os.getcwd();

FolderName              = "/" + lparam;
ProjPath  	        = BasicPath + FolderName;


FileNameIntraInter      = '/interband_dipole_full_evol.dat';
FileNameIntraInter1     = '/intraband_current_full_evol.dat';
FileNameLaser 	        = '/outlaserdata.dat';
FileParameters 		= '/setOfparameters.dat';
set_DataName 	        = '/SetData0';



#####################################################
FigureDir 	        = '/' + mparam + 'Figure';
IntraInterPath 	        = ProjPath  + FileNameIntraInter;
IntraInterPath1         = ProjPath  + FileNameIntraInter1;
LaserPath 	        = ProjPath  + FileNameLaser;
ParamPath 	        = ProjPath  + FileParameters;



try:
  os.mkdir( ProjPath + FigureDir );
except OSError as exc:
  if (exc.errno != errno.EEXIST):
    raise exc; 
  pass  



# LOADING DATA # 
InterC 	        = np.loadtxt( IntraInterPath );
IntraC          = np.loadtxt( IntraInterPath1 );
Laser 	        = np.loadtxt( LaserPath );
Params 		= np.loadtxt( ParamPath );

print "Params = \n" , Params


#####################################################
#shutil.copyfile( out, NewDir );
Phi0        = Params[2,0] ;
M0          = Params[2,1] ;
t1          = Params[1,5] ;
t2          = Params[1,6] ;
M0t2        = M0/t2;
ChernNo     = Params[2,4] ;
Eg          = Params[2,2]; 



#####################################################
print( "\nEg= ", Eg, " Chern No.= ", ChernNo, "; Phi= ", Phi0 );
print('\n............');
#####################################################


#print "size dim of ByHand", ByHand.shape
#Arranging/organizing data and parameters 
#Nt			= cs.Ntime;



#frequency and laser period
w0		= Params[0,2];
T0		= 2.*np.pi/w0;



print( "period: ", T0, " a.u.")
print( "mean-freq: ", w0, " a.u.")


#####################################################
#Laser intensity
E0          = np.sqrt( Params[0,1]/3.5e16);


#####################################################
#Time axis and electric field 
t		= Laser[:,0];
Efield		= Laser[:,1];

Nt          	= len(t)
dt 	        = t[1]-t[0]; 


print( "dt: ", dt, " a.u." );
print( "Ntime: ", Nt );


#####################################################
#Frequency axis 
wmax 	    	= np.pi/dt;
dw		= 2.*wmax/float(Nt);
wmin  	    	= -wmax;
wmax 	    	= +wmax;
w		= np.arange( wmin, wmax-dw, dw );
w           	= np.linspace( wmin, wmax-dw, num = Nt )




#####################################################
if (Nt != len(w)):
    print ("\n Ntime have to be equal to length (w)\n")
    exit()
#####################################################




print( "Lengths of omega axis, Nomega: ", w.shape)
print( "dw= ", dw )
print( '\nwmax={0:.2f}'.format(wmax) );




#####################################################
width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) );
plt.plot( t, Efield, 'red', lw=2 );

plt.xlabel('time (a.u.) ', fontsize=18);
plt.ylabel('Efield-MIR (a.u.) ', fontsize=18);
plt.tick_params(labelsize=18);

xaxmin      = t.min();    #
xaxmax      = t.max();    # 

plt.xlim(xaxmin, xaxmax);
plt.tight_layout();

filename0           = FigureDir + '/LaserField.png';
fileNamePicture     = ProjPath  + filename0;
plt.savefig( fileNamePicture ) ;



#####################################################
print( "\nlen of interC: ", InterC.shape)
print('\n---')




#####################################################
#Momentum Integration for the "remains" of MPI code
xinterC 	= InterC[:,0] - InterC[0,0];
yinterC     	= InterC[:,1] - InterC[0,1];

xJintra  	= IntraC[:,0] - IntraC[0,0];
yJintra     	= IntraC[:,1] - IntraC[0,1];




#####################################################
#Computing inter-current contribution from dipole polarization InterC1
xJinter		   = np.zeros( xinterC.shape, np.float );
yJinter        = np.zeros( yinterC.shape, np.float );



xJinter 	   = xinterC    #np.diff( xinterC )/dt;
yJinter        = yinterC    #np.diff( yinterC )/dt;



#####################################################
#Ploting the current oscillations
width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) )


p1, = plt.plot( t, Efield/max(Efield)*max(xJinter), 'red', lw=2 );
p2, = plt.plot( t, xJinter, 'green', lw = 1.5 );
p3, = plt.plot( t, yJinter, 'blue',  lw = 1.5 );



plt.legend([p1, p2, p3], ['$E_{L}$', '$J_{x,er}$', '$J_{y,er}$']);



plt.xlabel('time (a.u.) ', fontsize=18);
plt.ylabel('Jx/Jy (a.u.) ', fontsize=18);
plt.tick_params(labelsize=18);



xaxmin      = t.min();    #
xaxmax      = t.max();    # 
plt.xlim( xaxmin, xaxmax );


plt.tight_layout();




filename0           = FigureDir + '/CurrentOscillations.pdf';
fileNamePicture     = ProjPath + filename0; 
plt.savefig( fileNamePicture );






#####################################################
#Populations
width = 11
hight = width/1.62


#############################################################
## Computing harmonic spectra, inter and intra contribution
i 	= cmath.sqrt(-1)#complex(0,1);
print( "\ncomplex number: ", i )


#############################################################
#filtering dipole and current oscillations by means of 
#applying a smoth time mask over the beginning and end of pulses, this will avoid 
#high esporeous frequencies...
ta 	    = t[0]  + T0*5.0;#T0*3.0;
tb 	    = t[-1] - T0*5.0;#T0*3.0;

asigma  = T0*1.5
bsigma  = T0*1.5


xJinterMasked = masking_dipole(t, xJinter, ta, tb, asigma, bsigma);
yJinterMasked = masking_dipole(t, yJinter, ta, tb, asigma, bsigma);



ta     = t[0]  + T0*5.0;
tb     = t[-1] - T0*5.0;


asigma  = T0*1.50;
bsigma  = T0*1.50;

xJintraMasked = 1.*masking_dipole(t, xJintra, ta, tb, asigma, bsigma);
yJintraMasked = 1.*masking_dipole(t, yJintra, ta, tb, asigma, bsigma);



#####################################################
fig = plt.figure(figsize=(width,hight) )


p1, = plt.plot( t, Efield/max(Efield)*max(xJintra), 'red', lw=2 );
p2, = plt.plot( t, xJintraMasked, 'green', lw = 1.5 );
p3, = plt.plot( t, yJintraMasked, 'blue',  lw = 1.5 );
#plt.title(title_name,fontsize=18);



plt.legend([p1, p2, p3], ['$E_{L}$', '$J_{x,ra}$', '$J_{y,ra}$']);



plt.xlabel('time (a.u.) ', fontsize=18);
plt.ylabel('Jx/Jy (a.u.) ', fontsize=18);
plt.tick_params(labelsize=18);



xaxmin      = t.min();    #
xaxmax      = t.max();    #
plt.xlim( xaxmin, xaxmax );


plt.tight_layout();




filename0           = FigureDir + '/IntraCurrentOscillations.pdf';
fileNamePicture     = ProjPath + filename0;
plt.savefig( fileNamePicture );



#####################################################
#Ploting the current oscillations
width   = 11
hight   = width/1.62
fig     = plt.figure( figsize=(width,hight) )

p2,= plt.plot( t, xJinterMasked, 'green', lw = 1 );
p3,= plt.plot( t, yJinterMasked, 'blue',  lw = 1 );
plt.legend([p2, p3], ['$J_{x,er}$', '$J_{y,er}$'], fontsize = 18 );

plt.xlabel( 'time (a.u.) ', fontsize = 18 );
plt.ylabel( 'Filter Currents-Mask (a.u.) ', fontsize = 18 );
plt.tick_params( labelsize = 18 );

xaxmin      = t.min();
xaxmax      = t.max();

plt.xlim(xaxmin, xaxmax);
plt.tight_layout();


#########################################
###   Saving picture   ###
filename0           = FigureDir + '/CurrentOscillationsMF.pdf';
fileNamePicture     = ProjPath + filename0; 
plt.savefig( fileNamePicture );



#####################################################
##Calculating FFT 
xFFT_Jinter 	= fft( xJinterMasked )*dt;
yFFT_Jinter 	= fft( yJinterMasked )*dt;

xFFT_Jintra     = fft( xJintraMasked )*dt;
yFFT_Jintra     = fft( yJintraMasked )*dt;

xFFT_Jinter      = np.fft.fftshift( xFFT_Jinter );
yFFT_Jinter      = np.fft.fftshift( yFFT_Jinter );

xFFT_Jintra  	 = np.fft.fftshift( xFFT_Jintra );
yFFT_Jintra  	 = np.fft.fftshift( yFFT_Jintra );

xFullRadiation 	 = 1.*xFFT_Jinter   +  1.*xFFT_Jintra;
yFullRadiation   = 1.*yFFT_Jinter   +  1.*yFFT_Jintra;


############################
############################
hzero           = 9e-16
xSinter 	    = np.log10( abs( xFFT_Jinter)**2 + hzero );
ySinter 	    = np.log10( abs( yFFT_Jinter)**2 + hzero );

xSintra         = np.log10( abs( xFFT_Jintra)**2 + hzero );
ySintra         = np.log10( abs( yFFT_Jintra)**2 + hzero );

xSpectrum 	    = np.log10( abs( xFullRadiation )**2 + hzero );
ySpectrum       = np.log10( abs( yFullRadiation )**2 + hzero );
#####################################################




print( "Spectra shape= ",  xSinter.shape )
print( ";     frequency axis shape= ", w.shape )


w.shape             = (Nt,1)
xSpectrum.shape     = (Nt,1)
ySpectrum.shape     = (Nt,1)

OutputData      = w/w0  #[]#np.array([Nt])

OutputData      = np.concatenate( (OutputData,w/w0),    axis=1 );

OutputData      = np.concatenate( (OutputData,xSpectrum), axis=1 );

OutputData      = np.concatenate( (OutputData,ySpectrum), axis=1 );

ofname          = "/HHG__CNo" + str("%.1f"%ChernNo) +"__E0index__" +str("%.4d"%iparam) + "__.dat"

out                = ProjPath + ofname

np.savetxt(out, OutputData , fmt='%1.10e')

#outfile            = open( out, 'w' );
#for i in range(0,len(OutputData)):
#    outfile.write( str("%.12e"%OutputData[i,:])+"\n" );
#outfile.close()



NewDir = "./SetData0" + ofname
shutil.copyfile( out, NewDir );




#####################################################
############################
#Ploting the harmonic current radiations-oscillations
width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])


p1, = plt.plot( w/w0, xSpectrum, 'r-', lw = 2., label='$J_{x}$' );
p2, = plt.plot( w/w0, ySpectrum, 'b',  lw = 1.5, label='$J_{y}$' );
plt.legend([p1, p2], ['$J_{x}$', '$J_{y}$'],fontsize=18);


xp = 1.3
yp = -1.1


plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 );
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 );
plt.tick_params(labelsize = 28 );



xaxmin      = np.log10(hzero);      # controling harmonic-order
xaxmax      = +5;                   # controling harmonic-order

plt.ylim(xaxmin, xaxmax);
xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );

yticks0  = [-20., -15, -10.0, -5, 0, 5 ]#10., 15, 5.0];
plt.yticks( yticks0 );


if (mparam=='x'):
    plt.ylim([-12, .2])
if (mparam=='y'):
    plt.ylim([-10, .2])

plt.grid(True)


xaxmin      = 0;    # controling harmonic-order axis limits, down
xaxmax      = 37;   # controling harmonic-order axis limits, u


plt.xlim(xaxmin, xaxmax);
plt.ylim( -15,5 )

fname='/'+str(mparam)+'PaperHarmonicSpectrumCNo'+ str("%.1f"%ChernNo) + 'M_' + str("%.2f"%M0) + 'Phi_' + str("%.2f"%Phi0) + '.pdf'


filename1           = FigureDir + fname;
fileNamePicture     = ProjPath  + filename1; 


plt.savefig(fileNamePicture, dpi = 300);
shutil.copyfile( fileNamePicture, BasicPath + set_DataName + fname );



print( 'Eg = ',Eg, 'M0 = ',M0, ' phi0 = ', Phi0, 'C = ',ChernNo, 'mparam =', mparam, '-direction')

plt.show();
#####################################################
