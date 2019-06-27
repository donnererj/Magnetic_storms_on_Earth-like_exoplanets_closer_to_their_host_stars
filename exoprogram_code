import sys
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sunpy.time
import time
import pickle
import seaborn as sns
import os
import copy
import pdb
from matplotlib import rc
import scipy.io
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import math
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.dates import DayLocator, DateFormatter, HourLocator
import re
import decimal
import bisect 
from sys import exit
from matplotlib.ticker import AutoMinorLocator
from itertools import repeat
import pickle
from scipy.interpolate import interp1d
from scipy import interpolate
import copy
from scipy import stats
from scipy.stats import t


###get data###
def getcat(filename):
  #print('reading CAT')
  cat=scipy.io.readsav(filename, verbose='true')  
  #print('done CAT')
  return cat  
   
def decode_array(bytearrin):
 #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
 #make list of python lists with arbitrary length
 bytearrout= ['' for x in range(len(bytearrin))]
 for i in range(0,len(bytearrin)-1):
  bytearrout[i]=bytearrin[i].decode()
 #has to be np array so to be used with numpy "where"
 bytearrout=np.array(bytearrout)
 return bytearrout  
  
def time_to_num_cat(time_in):  

  j=0
  #time_str=np.empty(np.size(time_in),dtype='S19')
  time_str= ['' for x in range(len(time_in))]
  #=np.chararray(np.size(time_in),itemsize=19)
  time_num=np.zeros(np.size(time_in))
  
  for i in time_in:

   #convert from bytes (output of scipy.readsav) to string
   time_str[j]=time_in[j][0:16].decode()+':00'
   year=int(time_str[j][0:4])
   time_str[j]
   #convert time to sunpy friendly time and to matplotlibdatetime
   #only for valid times so 9999 in year is not converted
   #pdb.set_trace()
   if year < 2100:
    	  time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
   j=j+1  
   #the date format in matplotlib is e.g. 735202.67569444
   #this is time in days since 0001-01-01 UTC, plus 1.
   
   #return time_num which is already an array and convert the list of strings to an array
  return time_num, np.array(time_str)

###### duration values from other program (cloud1.py)
rh=np.array([1.00,4.80,0.87,1.00,1.00,1.06,1.08,0.62,1.00,1.58,1.59,0.94,1.00,1.98,0.87,2.49,0.51,1.00,0.73,0.96,0.47,1.08,0.73,1.08,0.72,0.96])
durationxs=np.array([28.05,24,7,15.06,21,23.49,25,9,19,14.03,17.06,22,32,24,28,50.57,9.02,10.3,24.4,30.47,20.28,30,11,14.42,30.3,30.55])

filename_icmecat='ICMECAT/HELCATS_ICMECAT_v10_SCEQ.sav'     #take your own file path 
i=getcat(filename_icmecat)


#access each element of the array see http://docs.scipy.org/doc/numpy/user/basics.rec.html
#i.icmecat['id']
#print(a.arrcat.dtype)
#print(i.icmecat.dtype)
#get spacecraft and planet positions
#pos=getcat('DATACAT/positions_2007_2018_HEEQ_6hours.sav')
#pos_time_num=time_to_num_cat(pos.time)[0]


####get all parameters from ICMECAT
iid=i.icmecat['id']

#need to decode all strings
iid=decode_array(iid)

###ID for events
isc=i.icmecat['sc_insitu'] #string
isc=decode_array(isc)
icme_start_time=i.icmecat['ICME_START_TIME']
[icme_start_time_num,icme_start_time_str]=time_to_num_cat(icme_start_time)
mo_start_time=i.icmecat['MO_START_TIME']
[mo_start_time_num,mo_start_time_str]=time_to_num_cat(mo_start_time)
mo_end_time=i.icmecat['MO_END_TIME']
[mo_end_time_num,mo_end_time_str]=time_to_num_cat(mo_end_time)
icme_end_time=i.icmecat['ICME_END_TIME']
[icme_end_time_num,icme_end_time_str]=time_to_num_cat(icme_end_time)
sc_heliodistance=i.icmecat['SC_HELIODISTANCE']
sc_long_heeq=i.icmecat['SC_LONG_HEEQ']
sc_lat_heeq=i.icmecat['SC_LAT_HEEQ']
mo_bmax=i.icmecat['MO_BMAX']
mo_bmean=i.icmecat['MO_BMEAN']
mo_bstd=i.icmecat['MO_BSTD']
mo_bzmean=i.icmecat['MO_BZMEAN']
mo_bzmin=i.icmecat['MO_BZMIN']
mo_duration=i.icmecat['MO_DURATION']
mo_mva_axis_long=i.icmecat['MO_MVA_AXIS_LONG']
mo_mva_axis_lat=i.icmecat['MO_MVA_AXIS_LAT']
mo_mva_ratio=i.icmecat['MO_MVA_RATIO']
sheath_speed=i.icmecat['SHEATH_SPEED']
sheath_speed_std=i.icmecat['SHEATH_SPEED_STD']
mo_speed=i.icmecat['MO_SPEED']
mo_speed_st=i.icmecat['MO_SPEED_STD']
sheath_density=i.icmecat['SHEATH_DENSITY']
sheath_density_std=i.icmecat['SHEATH_DENSITY_STD']
mo_density=i.icmecat['MO_DENSITY']
mo_density_std=i.icmecat['MO_DENSITY_STD']
sheath_temperature=i.icmecat['SHEATH_TEMPERATURE']
sheath_temperature_std=i.icmecat['SHEATH_TEMPERATURE_STD']
mo_temperature=i.icmecat['MO_TEMPERATURE']
mo_temperature_std=i.icmecat['MO_TEMPERATURE_STD']

ivexind=np.where(isc == 'VEX')
istaind=np.where(isc == 'STEREO-A')
istbind=np.where(isc == 'STEREO-B')
iwinind=np.where(isc == 'Wind')
imesind=np.where(isc == 'MESSENGER')
iulyind=np.where(isc == 'ULYSSES')

############################################
#get current directory
#os.system('pwd')
os.getcwd()
#define global variables from OMNI2 dataset
#see http://omniweb.gsfc.nasa.gov/html/ow_data.html
dataset=473376;
#global Variables
spot=np.zeros(dataset) 
btot=np.zeros(dataset) #floating points
bx=np.zeros(dataset) #floating points
by=np.zeros(dataset) #floating points
bz=np.zeros(dataset) #floating points
bzgsm=np.zeros(dataset) #floating points
bygsm=np.zeros(dataset) #floating points
speed=np.zeros(dataset) #floating points
speedx=np.zeros(dataset) #floating points
speed_phi=np.zeros(dataset) #floating points
speed_theta=np.zeros(dataset) #floating points
dst=np.zeros(dataset) #float
kp=np.zeros(dataset) #float
temp=np.zeros(dataset)
den=np.zeros(dataset) #float
pdyn=np.zeros(dataset) #float
year=np.zeros(dataset)
day=np.zeros(dataset)
hour=np.zeros(dataset)
t=np.zeros(dataset) #index time
times1=np.zeros(dataset) #datetime time
mo_start=np.zeros(dataset)

#distanceentry=float(input("choose value between 0.3 and 1.5 AU:"))
#print('distanceentry',distanceentry)  

#read in data from file 
def getdata():
 
 #FORMAT(2I4,I3,I5,2I3,2I4,14F6.1,F9.0,F6.1,F6.0,2F6.1,F6.3,F6.2, F9.0,F6.1,F6.0,2F6.1,F6.3,2F7.2,F6.1,I3,I4,I6,I5,F10.2,5F9.2,I3,I4,2F6.1,2I6,F5.1)
 #1963   1  0 1771 99 99 999 999 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 9999999. 999.9 9999. 999.9 999.9 9.999 99.99 9999999. 999.9 9999. 999.9 999.9 9.999 999.99 999.99 999.9  7  23    -6  119 999999.99 99999.99 99999.99 99999.99 99999.99 99999.99  0   3 999.9 999.9 99999 99999 99.9

 
 j=0
 #print('start reading variables from file')
 with open (os.path.dirname(os.path.abspath("omni2_all_years.dat"))+"/omni2_all_years.dat") as f:   #file path 
  for line in f:
   line = line.split() # to deal with blank 
   #print line #41 is Dst index, in nT
   dst[j]=line[40]
   kp[j]=line[38]
   
   if dst[j] == 99999: dst[j]=np.NaN
   #40 is sunspot number
   spot[j]=line[39]
   #if spot[j] == 999: spot[j]=NaN

   #25 is bulkspeed F6.0, in km/s
   speed[j]=line[24]
   if speed[j] == 9999: speed[j]=np.NaN
 
   #get speed angles F6.1
   speed_phi[j]=line[25]
   if speed_phi[j] == 999.9: speed_phi[j]=np.NaN

   speed_theta[j]=line[26]
   if speed_theta[j] == 999.9: speed_theta[j]=np.NaN
   #convert speed to GSE x see OMNI website footnote
   speedx[j] = - speed[j] * np.cos(np.radians(speed_theta[j])) * np.cos(np.radians(speed_phi[j]))

   #9 is total B  F6.1 also fill ist 999.9, in nT
   btot[j]=line[9]
   if btot[j] == 999.9: btot[j]=np.NaN

   #GSE components from 13 to 15, so 12 to 14 index, in nT
   bx[j]=line[12]
   if bx[j] == 999.9: bx[j]=np.NaN
   by[j]=line[13]
   if by[j] == 999.9: by[j]=np.NaN
   bz[j]=line[14]
   if bz[j] == 999.9: bz[j]=np.NaN
 
   #GSM
   bygsm[j]=line[15]
   if bygsm[j] == 999.9: bygsm[j]=np.NaN
 
   bzgsm[j]=line[16]
   if bzgsm[j] == 999.9: bzgsm[j]=np.NaN 	
 
 
   #24 in file, index 23 proton density /ccm
   den[j]=line[23]
   if den[j] == 999.9: den[j]=np.NaN
 
   #29 in file, index 28 Pdyn, F6.2, fill values are 99.99, in nPa
   pdyn[j]=line[28]
   if pdyn[j] == 99.99: pdyn[j]=np.NaN 		
 
   year[j]=line[0]
   day[j]=line[1]
   hour[j]=line[2]
   
#    temp[j]=line[24]
#    if temp[j]==9999999.: temp[j]=np.NaN
   
   j=j+1     

############################################################# 


def converttime():
 #http://docs.sunpy.org/en/latest/guide/time.html
 #http://matplotlib.org/examples/pylab_examples/date_demo2.html

 #print('convert time start')
 for index in range(0,dataset):
      #first to datetimeobject 
      timedum=datetime.datetime(int(year[index]), 1, 1) + datetime.timedelta(day[index] - 1) +datetime.timedelta(hours=hour[index])
      #then to matlibplot dateformat:
      times1[index] = matplotlib.dates.date2num(timedum)
      #mdates.date2num(sunpy.time.parse_time('2011-06-30 11:00'))
      #print time
      #print year[index], day[index], hour[index]
 #print('convert time done')   #for time conversion
  
  
 
def make_3dcore_dst(bz_in,v_in,pdyn_in,time_in):

#make_3dcore_dst(btoti,bxi,  bygsmi, bzgsmi, speedi, deni, timesi)

 #this makes from synthetic or observed solar wind the Dst index	
 #btot_in IMF total field, in nT
 #by_in - the IMF By field in nT
 #bz_in - the IMF Bz field in nT

 #v_in - the speed in km/s
 #vx_in - the solar wind speed x component (GSE or GSM?) in km/s
 #time_in - the time in matplotlib date format


 #add clock angle and total field for temerin and li

 #define variables
 Ey=np.zeros(len(bz_in))
 #dynamic pressure is constant
 pdyn1=pdyn_in
 dststar2=np.zeros(len(bz_in))
 dstcalc2=np.zeros(len(bz_in))
 
 #set all fields above 0 to 0 
 bz_in_negind=np.where(bz_in > 0)  
 
 #important: make a deepcopy because you manipulate the input variable
 bzneg=copy.deepcopy(bz_in)
 bzneg[bz_in_negind]=0

 #define interplanetary electric field 
 Ey=v_in*abs(bzneg)*1e-3; #now Ey is in mV/m
 #in EY steckt bz drin  #make 3D core
 
 ######## OBrien and McPherron 2000 ##########
 #constants
 Ec=0.49
 b=7.26  
 c=11  #nT
 for i in range(len(bz_in)-1):
  if Ey[i] > Ec:            #Ey in mV m
   Q=-4.4*(Ey[i]-Ec) 
  else: Q=0
  tau=2.4*np.exp(9.74/(4.69+Ey[i])) #tau in hours
  #OBrien abstract: Dst=Dst*+7.26P^1/2 -11
  #this is the ring current Dst
  #dststar2[i+1]=(Q-dststar2[i]/tau)+dststar2[i] #t is pro stunde, time intervall ist auch 1h
  deltat_hours=(time_in[i+1]-time_in[i])*24 #time_in is in days - convert to hours
  #deltat=dststar
  dststar2[i+1]=((Q-dststar2[i]/tau))*deltat_hours+dststar2[i] #t is pro stunde, time intervall ist auch 1h
  #this is the Dst of ring current and magnetopause currents 
  dstcalc2[i+1]=dststar2[i+1]+b*np.sqrt(pdyn1[i+1])-c; 
   #OBrien & McPherron
  #bei zeitpunkt 0 mit 0 anfangen, steigert sich, wie ändert solar wind input den jetzigen wert vom dst, bei d0=0 anfangen, input = solar wind, Ey=convection electric field=speed*Bz
  #
 return dstcalc2
   
#read in data from omni file -> 1 , from save_file -> 0
data_from_omni_file = 1 #hier 1 setzen beim ersten maL

if data_from_omni_file == 1:
 getdata()
 converttime()
 pickle.dump([spot,btot,bx,by,bz,bygsm,bzgsm,speed,speedx, dst,kp,den,pdyn,year,day,hour,times1], open(os.path.dirname(os.path.abspath("omni3save_july2016.p"))+"/omni3save_july2016.p", 'wb') ) #enter your own file path
else: [spot,btot,bx,by,bz,bygsm, bzgsm,speed,speedx, dst,kp,den,pdyn,year,day,hour,times1]= pickle.load(open(os.path.dirname(os.path.abspath("omni3save_july2016.p"))+"/omni3save_july2016.p", 'rb' ) ) #enter your own file path


#############################################
#startvalue
start='2007-Jun-1'     
ndays=2928

btoti=np.zeros(24*ndays)
bxi=np.zeros(24*ndays)
bygsmi=np.zeros(24*ndays)
bzgsmi=np.zeros(24*ndays)
speedi=np.zeros(24*ndays)
speedix=np.zeros(24*ndays)
timesi=np.zeros(24*ndays)
dsti=np.zeros(24*ndays)
kpi=np.zeros(24*ndays)
# tempi=np.zeros(24*ndays)
deni=np.zeros(24*ndays)
pdyni=np.zeros(24*ndays)
s=mdates.date2num(sunpy.time.parse_time(start))
#"i" stands for interval

ind=int(np.where(s==times1)[0])
btoti=btot[ind:ind+24*ndays]
bxi=bx[ind:ind+24*ndays]
bygsei=by[ind:ind+24*ndays]
bzgsei=bz[ind:ind+24*ndays]
bygsmi=bygsm[ind:ind+24*ndays]
bzgsmi=bzgsm[ind:ind+24*ndays]
speedi=speed[ind:ind+24*ndays]
speedix=speedx[ind:ind+24*ndays]
deni=den[ind:ind+24*ndays]
timesi=times1[ind:ind+24*ndays]
pdyni=pdyn[ind:ind+24*ndays]
dsti=dst[ind:ind+24*ndays]
kpi=kp[ind:ind+24*ndays]

#IMF clock angle
thetai=np.arctan2(bygsmi,bzgsmi)
thetai_deg=thetai*180/np.pi

#deni[np.isnan(deni)]=5*1**-2	
speedi[np.isnan(speedi)]=400	#400 slow solar wind 

bzgsmi[np.isnan(bzgsmi)]=5    #5 nT for background solar wind

# #p[nanopascal]=2.0*1e-6*n
# newsp=speednj**2
# pdyn_n=2.0*1e-6*deni*speedi	






#def make_3dcore_dst(btot_in,bx_in, by_in,bz_in,v_in,vx_in,density_in,time_in):
dst_obrien=make_3dcore_dst(bzgsmi,speedi, pdyni,timesi)
#voriger 3dcore (btoti,bxi,bygsmi,bzgsmi,speedi,speedix,deni,timesi)  
 



###define start & end time for MO's within time range###
pfaffe=timesi.astype(int)					#take values as integers in order to calculate with it 
maffe=mo_start_time_num[0:164].astype(int)  #intervall 0:164 to take only wind values 
klaffe=mo_end_time_num[0:164].astype(int)

khaki=np.where(np.in1d(maffe,pfaffe))   	
sharon=np.where(np.in1d(klaffe,pfaffe)) 

jau=mo_start_time_num[khaki]
joi=mo_end_time_num[sharon]

#define start value for MO 
start1=[]
for i in range(len(jau)):
	aleu=bisect.bisect(timesi,jau[i])
	#print('index aleu', aleu)
	start1.append(aleu)
st=np.asarray(start1)   

#define end value from MO 
end1=[]
for i in range(len(joi)):
	gau=bisect.bisect(timesi,joi[i])
	end1.append(gau)
en=np.asarray(end1)

if len(st)>len(en):
	print('\x1b[1;31m'+'Must have equal number of start & end dates, choose more ndays'+'\x1b[0m')
	
### B scale factor with bmean ### 
Bq=[]	
#buri=bur*(1.0**bmeanfit[0])		#take fit at 1 AU to compare with other data and get factor for each distance 
buffalo=8.76*1**(-1.72)
facto=[]
for r in np.arange(0.3,1.6,0.1):		 
	Bq=8.76*r**(-1.72)				 
	factor1=Bq/buffalo
	facto.append(factor1)
bfacmean=np.asarray(facto) 


### B scale factor with bmax ###
Br=[]							
buffalo2=12.13*1**(-1.76)
bfact=[]
for r in np.arange(0.3,1.6,0.1):    #define distances between 0.3AU-1.5AU
	Br=12.13*r**(-1.76)
	fac=Br/buffalo2
	bfact.append(fac)
bfac=np.asarray(bfact)   			#B scale factor  


### duration scale factor ###					
fadr=[]
durf=22.7*1+5.9						#duration fit at 1AU for comparison in order to get scale factor, here we take a linear fit
for r in np.arange(0.3,1.6,0.1):   
	dr=22.7*r+5.9
	facdr=dr/durf
	fadr.append(facdr)	
facdur=np.asarray(fadr)  			#scale factor for duration  


### define intervalls for speed scalation ###
sp=[]      							#create factor for multiplication with speed 
speedy1=np.linspace(2,1,8)    		#from 0.3 AU until 1 AU starting with factor being=2 and ending at 1 AU with factor being=1
speedy2=np.ones(5)            		#after 1 AU (1.1AU until 1,5AU) factor stays =1
sp.extend(speedy1)            		#using extend instead append to get correct array
sp.extend(speedy2)
spd=np.transpose(sp)		  		#need to transpose to get it in horizontal array 


### density scale factor ###
dens1=6.5*1**(-2.4)    				#leitner et. al. takes 6.5, 6.63;  bothmer & schwenn takes 6.47
dens3=[]
for r in np.arange(0.3,1.6,0.1):
	dens2=6.5*r**(-2.4) 
	facdens1=dens2/dens1
	dens3.append(facdens1)    
facdens=np.asarray(dens3)   		#scale factor for density 

#define lists for value calculation in function
bxbx=[]
byby=[]
bzbz=[]
btotb=[]
tnewj=[]
tai=[]
tufi=[]
dr=[]
timex=[]
bey=[]
bez=[]	
betot=[]
facden=[]
dens4=[]
spni=[]
spnj=[]
bxu=[]
bxnew=[]

idefix=range(len(st))
timeaxis=[]
for a in idefix:       
	time2=timesi[st[a]:en[a]]-timesi[st[a]]			
	dt2=time2[2]-time2[1]						#definine time step 

	#t3=np.zeros(len(time2))
	t2=[]
	
	for i in range(len(time2)):
		t3=i*(dt2*24.*facdur)   					 
		t2.append(t3)				
		t4=np.asarray(t2) 	
		t5=t4.transpose()
	timeaxis.append(t5)	

####create function for MO scaling 
def scalation(v1,v2):
	bko=[]
	for q in range(len(st)):  	
		ix=np.arange(len(v1[st[q]:en[q]]))		#v1=bxi
		bex=[]
		for i in range(len(ix)): 				#scale between starting point of bxi from MO and end point of MO for bxi 
			bxj=v1[i+st[q]]*v2					#v2=bfac			 
			bex.append(bxj)
			bexi=np.asarray(bex) 
			baxi=bexi.transpose()		
		bko.append(baxi)
	return bko	
	

###create lists for scalation function####		
bx_mean_scaled=[]
bx_max_scaled=[]
by_mean_scaled=[]
by_max_scaled=[]
bz_mean_scaled=[]
bz_max_scaled=[]
btot_mean_scaled=[]
btot_max_scaled=[]
density_scaled=[]
speed_scaled=[]

### copy data ###
bxr=copy.copy(bxi)
byr=copy.copy(bygsmi)
bzr=copy.copy(bzgsmi)
btotr=copy.copy(btoti)
densr=copy.copy(deni)
speedr=copy.copy(speedi)
maxbxr=copy.copy(bxi)
maxbyr=copy.copy(bygsmi)
maxbzr=copy.copy(bzgsmi)
maxbtotr=copy.copy(btoti)

### apply scalation function ###
bxmax1=scalation(bxr,bfac)
bxmean1=scalation(bxr,bfacmean)
bymax1=scalation(byr,bfac)
bymean1=scalation(byr,bfacmean)
bzmax1=scalation(bzr,bfac)
bzmean1=scalation(bzr,bfacmean)
btotmax1=scalation(btotr,bfac)
btotmean1=scalation(btotr,bfacmean)
dens1=scalation(densr,facdens)
speed1=scalation(speedr,spd)

### attach parameters to lists ###
bx_mean_scaled.append(bxmean1)
bx_max_scaled.append(bxmax1)
by_mean_scaled.append(bymean1)
by_max_scaled.append(bymax1)
bz_mean_scaled.append(bzmean1)
bz_max_scaled.append(bzmax1)
btot_mean_scaled.append(btotmean1)
btot_max_scaled.append(btotmax1)	
density_scaled.append(dens1)
speed_scaled.append(speed1)
###output = lists with the scaled values per MO 

### duplicate original arrays ###
bxro=copy.copy(bxi)
byro=copy.copy(bygsmi)
bzro=copy.copy(bzgsmi)
btotro=copy.copy(btoti)
densro=copy.copy(deni)
speedro=copy.copy(speedi)
maxbxro=copy.copy(bxi)
maxbyro=copy.copy(bygsmi)
maxbzro=copy.copy(bzgsmi)
maxbtotro=copy.copy(btoti)

bachi=bxro
bychi=byro
bzchi=bzro
btotchi=btotro
denschi=densro
speedchi=speedro
maxbx=maxbxro
maxby=maxbyro
maxbz=maxbzro
maxbtot=maxbtotro

### exchange scaled MO values within original MO values in array ###
fau=[x/24 for x in timeaxis]
tas=copy.copy(timesi)
ais=np.arange(0,13,1)

### create function for MO scaling to different distances ###
def scalationAU(AUx): 
### scaled for 0.3 AU until 1.5 AU ###
	ri=np.arange(0.3,1.6,0.1)
	for s in range(len(st)):
		bachi[st[s]:en[s]]=bx_mean_scaled[0][s][AUx]
		bychi[st[s]:en[s]]=by_mean_scaled[0][s][AUx]
		bzchi[st[s]:en[s]]=bz_mean_scaled[0][s][AUx]
		btotchi[st[s]:en[s]]=btot_mean_scaled[0][s][AUx]
		denschi[st[s]:en[s]]=density_scaled[0][s][AUx]
		speedchi[st[s]:en[s]]=speed_scaled[0][s][AUx]
		maxbx[st[s]:en[s]]=bx_max_scaled[0][s][AUx]
		maxby[st[s]:en[s]]=by_max_scaled[0][s][AUx]
		maxbz[st[s]:en[s]]=bz_max_scaled[0][s][AUx]
		maxbtot[st[s]:en[s]]=btot_max_scaled[0][s][AUx]
	
	### define new time line with MO scalation ###


	for p in range(len(st)):
		oluj=[tas[st[p]]+x for x in fau[p][AUx]]
		tas[st[p]:en[p]]=oluj

	tn=copy.copy(tas)
	asd=[]
	nas=[]
	xul=[]
	lako=[]
	lakn=[]
	nasi=[]
	csd=[]	

	### start to fill with nans and therefore start and end time change ### 
	for p in range(len(st)):
		
	
		### first st & en are still the right ones ###
		if p==0:
			if (tn[en[0]])>(tn[en[0]-1]+0.04):    					#make sure the attached values aren't bigger than the next values 
				dag=np.arange(tn[en[0]-1],tn[en[0]],0.04)			#create array with time entries between last value from MO and first next value after MO
				ful=dag[1:]											#dont take 1 value because it would be double entry 
				nas.append(ful)
				if [x<tn[en[0]] for x in ful]:
					#tn[en[0]:en[0]]=nas[0]    
					tnx=np.insert(tn,en[0],[x for x in nas[0]])		#insert values between last value from MO and first next value after MO
					lau=np.empty(len(nas[0]))
					lau[:]=np.nan	
					asd.append(lau)									#length of this new array gives length for nans to insert in b components 

					### insert nans at specified index -> nan array has length from new array so it is always dependent on time ### 	
					newbx=(np.random.random(len(asd[0]))*5-2.5)*2
					newby=(np.random.random(len(asd[0]))*5-2.5)*2
					newbz=(np.random.random(len(asd[0]))*5-2.5)*2
					newbtot=np.sqrt(newbx**2+newby**2+newbz**2)
					newdens=np.full(len(asd[0]),(5*ri[AUx]**-2),dtype=float) #5protonen/ccm at 1 AU 
					newspeed= np.full(len(asd[0]),400,dtype=int)
					###creates random numbers for the lengths of the time axis section  
					###the new background solar wind  
					
					buchi=np.insert(bachi,en[0],[x for x in newbx])
					buychi=np.insert(bychi,en[0],[x for x in newby])
					buzchi=np.insert(bzchi,en[0],[x for x in newbz])
					butotchi=np.insert(btotchi,en[0],[x for x in newbtot])
					buspeedchi=np.insert(speedchi,en[0],[x for x in newspeed])
					budensitychi=np.insert(denschi,en[0],[x for x in newdens])
					maxbxi=np.insert(maxbx,en[0],[x for x in newbx])
					maxbyi=np.insert(maxby,en[0],[x for x in newby])
					maxbzi=np.insert(maxbz,en[0],[x for x in newbz])
					maxbtoti=np.insert(maxbtot,en[0],[x for x in newbtot])
				
		### st and en change =>creation of new st(st_x) and en(en_x) ###
		if p!=0:
			if p==1:
				rst=np.where(tnx==timesi[st[p]])  	     
			if p!=0 and p!=1:
				rst=np.where(tnxi==timesi[st[p]]) 
			#print('p',p)
			#print('rst',rst)
			lako.extend(rst)
			st_x=np.asarray(lako)
		
			ren=st_x[p-1]+len(fau[p][0])   			
			lakn.extend(ren)
			en_x=np.asarray(lakn)
		
			### extra case in order to exchange values for tnxi ###
			if p==1: 
				if (tnx[en_x[0]])>(tnx[en_x[0]-1]+0.04):
					dagi=np.arange(tnx[en_x[0]-1],tnx[en_x[0]],0.04)				
					fuli=dagi[1:]
					nasi.append(fuli)
			
					if [x<tnx[en_x[0]] for x in fuli]:   
						tnxi=np.insert(tnx,en_x[0],fuli)
						#tn[en_x[p-1]:en_x[p-1]]=fuli
					
						laui=np.empty(len(fuli))
						laui[:]=np.nan	
						csd.append(laui)
						
						newbx=(np.random.random(len(csd[0]))*5-2.5)*2
						newby=(np.random.random(len(csd[0]))*5-2.5)*2
						newbz=(np.random.random(len(csd[0]))*5-2.5)*2
						newbtot=np.sqrt(newbx**2+newby**2+newbz**2)
						newdens=np.full(len(csd[0]),(5*ri[AUx]**-2),dtype=float) #5protonen/ccm at 1 AU 
						newspeed= np.full(len(csd[0]),400,dtype=int)
						
						bxni=np.insert(buchi,en_x[0],[x for x in newbx])
						byni=np.insert(buychi,en_x[0],[x for x in newby])
						bzni=np.insert(buzchi,en_x[0],[x for x in newbz])
						btotni=np.insert(butotchi,en_x[0],[x for x in newbtot])
						speedni=np.insert(buspeedchi,en_x[0],[x for x in newspeed])
						densityni=np.insert(budensitychi,en_x[0],[x for x in newdens])   #csd[p-1]
						maxbxn=np.insert(maxbxi,en_x[0],[x for x in newbx])
						maxbyn=np.insert(maxbyi,en_x[0],[x for x in newby])
						maxbzn=np.insert(maxbzi,en_x[0],[x for x in newbz])
						maxbtotn=np.insert(maxbtoti,en_x[0],[x for x in newbtot])

			### output=finalized scaled arrays ###
			if p!=0 and p!=1:
				if (tnxi[en_x[p-1]])>(tnxi[en_x[p-1]-1]+0.04):
					dagi=np.arange(tnxi[en_x[p-1]-1],tnxi[en_x[p-1]],0.04)				
					fuli=dagi[1:]
				#print('ful',ful)
					#print('fuli',fuli)
					nasi.append(fuli)
			
					if [x<tnxi[en_x[p-1]] for x in fuli]:   

						tnxi=np.insert(tnxi,en_x[p-1],fuli)

						laui=np.empty(len(fuli))
						laui[:]=np.nan	
						csd.append(laui)
						
						
						newbx=(np.random.random(len(csd[p-1]))*5-2.5)*2
						newby=(np.random.random(len(csd[p-1]))*5-2.5)*2
						newbz=(np.random.random(len(csd[p-1]))*5-2.5)*2
						newbtot=np.sqrt(newbx**2+newby**2+newbz**2)
						
						newdens=np.full(len(csd[p-1]),(5*ri[AUx]**-2),dtype=float) #5protonen/ccm at 1 AU 
						newspeed= np.full(len(csd[p-1]),400,dtype=int)
						
						
						bxni=np.insert(bxni,en_x[p-1],[x for x in newbx])
						byni=np.insert(byni,en_x[p-1],[x for x in newby])
						bzni=np.insert(bzni,en_x[p-1],[x for x in newbz])
						btotni=np.insert(btotni,en_x[p-1],[x for x in newbtot])
						speedni=np.insert(speedni,en_x[p-1],[x for x in newspeed])
						densityni=np.insert(densityni,en_x[p-1],[x for x in newdens])
						maxbxn=np.insert(maxbxn,en_x[p-1],[x for x in newbx])
						maxbyn=np.insert(maxbyn,en_x[p-1],[x for x in newby])
						maxbzn=np.insert(maxbzn,en_x[p-1],[x for x in newbz])
						maxbtotn=np.insert(maxbtotn,en_x[p-1],[x for x in newbtot])
			
				
		if len(st)==1:
			bxni=copy.copy(buchi)	
			byni=copy.copy(buychi)	
			bzni=copy.copy(buzchi)	
			btotni=copy.copy(butotchi)	
			speedni=copy.copy(buspeedchi)	
			densityni=copy.copy(budensitychi)	
			maxbxn=copy.copy(maxbxi)	
			maxbyn=copy.copy(maxbyi)	
			maxbzn=copy.copy(maxbzi)	
			maxbtotn=copy.copy(maxbtoti)	
			tnxi=copy.copy(tnx)			
	
	
### create distance AU array ###
	AU=np.arange(0.3,1.6,0.1)
### create start and end values ###
	st_n=[]
	en_n=[]
	en_nan=[]

	### new start values for scaled arrays ###
	st_n.append(st[0])
	if len(st)>1:
		st_n.extend(st_x)

	### new end values for scaled arrays  ==startvalues for nans, endvaluebefore nans==en_n[x]-1 ###
	en_n.append(en[0])
	if len(en)>1:
		en_n.extend(en_x)

	#end values after nans for scaled arrays
	en_nan.append(en[0]+len(nas[0]))
	if len(en)>1:
		for x in range(len(en_x)):
			en_nan.append(en_x[x]+len(nasi[x]))
	
	
	ri=np.arange(0.3,1.6,0.1)
	
	bxj=copy.copy(bxni)
	byj=copy.copy(byni)
	bzj=copy.copy(bzni)
	btotj=copy.copy(btotni)
	densj=copy.copy(densityni)
	speedj=copy.copy(speedni)
	maxbxj=copy.copy(maxbxn)
	maxbyj=copy.copy(maxbyn)
	maxbzj=copy.copy(maxbzn)
	maxbtotj=copy.copy(maxbtotn)
	tnxj=copy.copy(tnxi)
	
	
	for a in range(len(st_n)):	
	
		if a==0 and len(st)>1:
			bxj[0:st_n[0]]=[x*ri[AUx]**-2 for x in bxj[0:st_n[0]]]
			byj[0:st_n[0]]=[x*ri[AUx]**-2 for x in byj[0:st_n[0]]]
			bzj[0:st_n[0]]=[x*ri[AUx]**-2 for x in bzj[0:st_n[0]]]
			btotj[0:st_n[0]]=[x*ri[AUx]**-2 for x in btotj[0:st_n[0]]]
			# speedj[0:st_n[0]]=[x*ri[AUx]**-2 for x in speedj[0:st_n[0]]]
# 			densj[0:st_n[0]]=[x*ri[AUx]**-2 for x in densj[0:st_n[0]]]
			maxbxj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbxj[0:st_n[0]]]
			maxbyj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbyj[0:st_n[0]]]
			maxbzj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbzj[0:st_n[0]]]
			maxbtotj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbtotj[0:st_n[0]]]
 
			bxj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in bxj[int(en_n[a]):int(st_n[a+1])]]
			byj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in byj[int(en_n[a]):int(st_n[a+1])]]
			bzj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in bzj[int(en_n[a]):int(st_n[a+1])]]
			btotj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in btotj[int(en_n[a]):int(st_n[a+1])]]
			# speedj[en_n[a]:st_n[a+1]]=[x*ri[AUx]**-2 for x in speedj[en_n[a]:st_n[a+1]]]
# 			densj[en_n[a]:st_n[a+1]]=[x*ri[AUx]**-2 for x in densj[en_n[a]:st_n[a+1]]]
			maxbxj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbxj[int(en_n[a]):int(st_n[a+1])]]
			maxbyj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbyj[int(en_n[a]):int(st_n[a+1])]]
			maxbzj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbzj[int(en_n[a]):int(st_n[a+1])]]
			maxbtotj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbtotj[int(en_n[a]):int(st_n[a+1])]]	
			
			
		if a==0 and len(st)==1:
			bxj[0:st_n[0]]=[x*ri[AUx]**-2 for x in bxj[0:st_n[0]]]
			byj[0:st_n[0]]=[x*ri[AUx]**-2 for x in byj[0:st_n[0]]]
			bzj[0:st_n[0]]=[x*ri[AUx]**-2 for x in bzj[0:st_n[0]]]
			btotj[0:st_n[0]]=[x*ri[AUx]**-2 for x in btotj[0:st_n[0]]]
			# speedj[0:st_n[0]]=[x*ri[AUx]**-2 for x in speedj[0:st_n[0]]]
# 			densj[0:st_n[0]]=[x*ri[AUx]**-2 for x in densj[0:st_n[0]]]
			maxbxj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbxj[0:st_n[0]]]
			maxbyj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbyj[0:st_n[0]]]
			maxbzj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbzj[0:st_n[0]]]
			maxbtotj[0:st_n[0]]=[x*ri[AUx]**-2 for x in maxbtotj[0:st_n[0]]]
		
			########## 
		if a!=0 and a!=range(len(st_n))[-1]:  
			bxj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in bxj[int(en_n[a]):int(st_n[a+1])]]
			byj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in byj[int(en_n[a]):int(st_n[a+1])]]
			bzj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in bzj[int(en_n[a]):int(st_n[a+1])]]
			btotj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in btotj[int(en_n[a]):int(st_n[a+1])]]
			# speedj[en_n[a]:st_n[a+1]]=[x*ri[AUx]**-2 for x in speedj[en_n[a]:st_n[a+1]]]
# 			densj[en_n[a]:st_n[a+1]]=[x*ri[AUx]**-2 for x in densj[en_n[a]:st_n[a+1]]]
			maxbxj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbxj[int(en_n[a]):int(st_n[a+1])]]
			maxbyj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbyj[int(en_n[a]):int(st_n[a+1])]]
			maxbzj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbzj[int(en_n[a]):int(st_n[a+1])]]
			maxbtotj[int(en_n[a]):int(st_n[a+1])]=[x*ri[AUx]**-2 for x in maxbtotj[int(en_n[a]):int(st_n[a+1])]]
			
	
	if len(st)==1:  #in case there is only 1 MO 
		bxj[en_n[0]:]=[x*ri[AUx]**-2 for x in bxj[en_n[0]:]]
		byj[en_n[0]:]=[x*ri[AUx]**-2 for x in byj[en_n[0]:]]
		bzj[en_n[0]:]=[x*ri[AUx]**-2 for x in bzj[en_n[0]:]]
		btotj[en_n[0]:]=[x*ri[AUx]**-2 for x in btotj[en_n[0]:]]
		# speedj[en_n[-1]:]=[x*ri[AUx]**-2 for x in speedj[en_n[-1]:]]
	# 	densj[en_n[-1]:]=[x*ri[AUx]**-2 for x in densj[en_n[-1]:]]
		maxbxj[en_n[0]:]=[x*ri[AUx]**-2 for x in maxbxj[en_n[0]:]]
		maxbyj[en_n[0]:]=[x*ri[AUx]**-2 for x in maxbyj[en_n[0]:]]
		maxbzj[en_n[0]:]=[x*ri[AUx]**-2 for x in maxbzj[en_n[0]:]]
		maxbtotj[en_n[0]:]=[x*ri[AUx]**-2 for x in maxbtotj[en_n[0]:]]
	
	if len(st)>1: #in case there is more than 1 MO 
		bxj[en_n[-1]:]=[x*ri[AUx]**-2 for x in bxj[en_n[-1]:]]
		byj[en_n[-1]:]=[x*ri[AUx]**-2 for x in byj[en_n[-1]:]]
		bzj[en_n[-1]:]=[x*ri[AUx]**-2 for x in bzj[en_n[-1]:]]
		btotj[en_n[-1]:]=[x*ri[AUx]**-2 for x in btotj[en_n[-1]:]]
		# speedj[en_n[-1]:]=[x*ri[AUx]**-2 for x in speedj[en_n[-1]:]]
	# 	densj[en_n[-1]:]=[x*ri[AUx]**-2 for x in densj[en_n[-1]:]]
		maxbxj[en_n[-1]:]=[x*ri[AUx]**-2 for x in maxbxj[en_n[-1]:]]
		maxbyj[en_n[-1]:]=[x*ri[AUx]**-2 for x in maxbyj[en_n[-1]:]]
		maxbzj[en_n[-1]:]=[x*ri[AUx]**-2 for x in maxbzj[en_n[-1]:]]
		maxbtotj[en_n[-1]:]=[x*ri[AUx]**-2 for x in maxbtotj[en_n[-1]:]]
	
	bxnj=copy.copy(bxj)	
	bynj=copy.copy(byj)
	bznj=copy.copy(bzj)
	btotnj=copy.copy(btotj)
	speednj=copy.copy(speedj)
	densnj=copy.copy(densj)
	maxbxnj=copy.copy(maxbxj)
	maxbynj=copy.copy(maxbyj)
	maxbznj=copy.copy(maxbzj)
	maxbtotnj=copy.copy(maxbtotj)
		
		
	### NaN's for density are being replaced with 5protons/cm^-3 -> newdens=np.full(len(csd[p-1]),(5*ri[AUx]**-2),dtype=float) ### 	
	densnj[np.isnan(densnj)]=5*ri[AUx]**-2	   #johnstone et. al 2015
	speednj[np.isnan(speednj)]=400
	#p[nanopascal]=2.0*1e-6*n
	newsp=speednj**2
	pdyn_n=2.0*1e-6*densnj*newsp	
	####### create dst values #######
	dstmean=make_3dcore_dst(bznj,speednj,pdyn_n,tnxj)  
	dstmax=make_3dcore_dst(maxbznj,speednj,pdyn_n,tnxj)  
	
	startneu=np.asarray(st_n)
	stneu=startneu.astype(int)
	
	endeneu=np.asarray(en_n)
	enneu=endeneu.astype(int)
	
	endenanneu=np.asarray(en_nan)
	ennanneu=endenanneu.astype(int)
	
	return (bxnj,bynj,bznj,btotnj,speednj,densnj,pdyn_n,maxbxnj,maxbynj,maxbznj,maxbtotnj,tnxj,stneu,enneu,ennanneu,dstmean,dstmax)
	
### create arrays for different distances ###
bxni03,byni03,bzni03,btotni03,speedni03,densityni03,pdyn_n03,maxbxn03,maxbyn03,maxbzn03,maxbtotn03,tnxi03,st_n03,en_n03,en_nan03,dstmean03,dstmax03=scalationAU(ais[0])	
bxni04,byni04,bzni04,btotni04,speedni04,densityni04,pdyn_n04,maxbxn04,maxbyn04,maxbzn04,maxbtotn04,tnxi04,st_n04,en_n04,en_nan04,dstmean04,dstmax04=scalationAU(ais[1])		
bxni05,byni05,bzni05,btotni05,speedni05,densityni05,pdyn_n05,maxbxn05,maxbyn05,maxbzn05,maxbtotn05,tnxi05,st_n05,en_n05,en_nan05,dstmean05,dstmax05=scalationAU(ais[2])	
bxni06,byni06,bzni06,btotni06,speedni06,densityni06,pdyn_n06,maxbxn06,maxbyn06,maxbzn06,maxbtotn06,tnxi06,st_n06,en_n06,en_nan06,dstmean06,dstmax06=scalationAU(ais[3])		
bxni07,byni07,bzni07,btotni07,speedni07,densityni07,pdyn_n07,maxbxn07,maxbyn07,maxbzn07,maxbtotn07,tnxi07,st_n07,en_n07,en_nan07,dstmean07,dstmax07=scalationAU(ais[4])		
bxni08,byni08,bzni08,btotni08,speedni08,densityni08,pdyn_n08,maxbxn08,maxbyn08,maxbzn08,maxbtotn08,tnxi08,st_n08,en_n08,en_nan08,dstmean08,dstmax08=scalationAU(ais[5])		
bxni09,byni09,bzni09,btotni09,speedni09,densityni09,pdyn_n09,maxbxn09,maxbyn09,maxbzn09,maxbtotn09,tnxi09,st_n09,en_n09,en_nan09,dstmean09,dstmax09=scalationAU(ais[6])

r_i=np.arange(0.3,1.6,0.1)

### calculates maximum Dst-factor ###
def neuden(dsm,dn):
	#dsn=[]	
	#fdsne=copy.copy(dn)
	factod=[]
	#dsn=copy.copy(dsm)				 
	factor12=np.nanmean(dsm)/np.nanmean(dn)
	factod.append(factor12)
	dstfac=np.asarray(factod) 
	return(dstfac)

### Dst-index at 1 AU ###		
dst1n=copy.copy(dsti)

### Dst-index compared with Obrien & McPherron model ###
dst03mneu=neuden(dstmean03,dst_obrien)
dst04mneu=neuden(dstmean04,dst_obrien)
dst05mneu=neuden(dstmean05,dst_obrien)
dst06mneu=neuden(dstmean06,dst_obrien)
dst07mneu=neuden(dstmean07,dst_obrien)
dst08mneu=neuden(dstmean08,dst_obrien)
dst09mneu=neuden(dstmean09,dst_obrien)

### Dst-index compared with original Dst observed at 1 AU ###
dst3h=neuden(dstmean03,dst1n)
dst4h=neuden(dstmean04,dst1n)
dst5h=neuden(dstmean05,dst1n)
dst6h=neuden(dstmean06,dst1n)
dst7h=neuden(dstmean07,dst1n)
dst8h=neuden(dstmean08,dst1n)
dst9h=neuden(dstmean09,dst1n)

### Table with Dst-index values ###
dstftex=pd.DataFrame([["%.1f"%dst3h[0],"%.1f"%dst4h[0],"%.1f"%dst5h[0],"%.1f"%dst6h[0],"%.1f"%dst7h[0],"%.1f"%dst8h[0],"%.1f"%dst9h[0]]],
	index=['Dst-index'],columns=['0.3AU','0.4AU','0.5AU','0.6AU','0.7AU','0.8AU','0.9AU'])

with open('dst_factors.tex','w') as tf:
	tf.write(dstftex.to_latex())

### create dictionary / pickle file ###
### exoprogram file ### 
dict = {'time_normal': timesi,'bx_normal': bxi,'by_normal': bygsmi,'bz_normal': bzgsmi,'btot_normal': btoti, 
'speed_normal': speedi,'density_normal': deni,'pdyn_normal': pdyni,'dst_normal':dsti,'dst_obrien_normal':dst_obrien,
'time_i_scaled_03AU': tnxi03,'time_i_scaled_04AU': tnxi04,'time_i_scaled_05AU': tnxi05,'time_i_scaled_06AU': tnxi06,'time_i_scaled_07AU': tnxi07,'time_i_scaled_08AU': tnxi08,'time_i_scaled_09AU': tnxi09,   
'bx_i_mean_scaled_03AU': bxni03,'bx_i_mean_scaled_04AU': bxni04,'bx_i_mean_scaled_05AU': bxni05,'bx_i_mean_scaled_06AU': bxni06,'bx_i_mean_scaled_07AU': bxni07,'bx_i_mean_scaled_08AU': bxni08,'bx_i_mean_scaled_09AU': bxni09, 
'by_i_mean_scaled_03AU': byni03,'by_i_mean_scaled_04AU': byni04,'by_i_mean_scaled_05AU': byni05,'by_i_mean_scaled_06AU': byni06,'by_i_mean_scaled_07AU': byni07,'by_i_mean_scaled_08AU': byni08,'by_i_mean_scaled_09AU': byni09,
'bz_i_mean_scaled_03AU': bzni03,'bz_i_mean_scaled_04AU': bzni04,'bz_i_mean_scaled_05AU': bzni05,'bz_i_mean_scaled_06AU': bzni06,'bz_i_mean_scaled_07AU': bzni07,'bz_i_mean_scaled_08AU': bzni08,'bz_i_mean_scaled_09AU': bzni09,
'btot_i_mean_scaled_03AU': btotni03,'btot_i_mean_scaled_04AU': btotni04,'btot_i_mean_scaled_05AU': btotni05,'btot_i_mean_scaled_06AU': btotni06,'btot_i_mean_scaled_07AU': btotni07,'btot_i_mean_scaled_08AU': btotni08,'btot_i_mean_scaled_09AU': btotni09,
'bx_i_max_scaled_03AU': maxbxn03,'bx_i_max_scaled_04AU': maxbxn04,'bx_i_max_scaled_05AU': maxbxn05,'bx_i_max_scaled_06AU': maxbxn06,'bx_i_max_scaled_07AU': maxbxn07,'bx_i_max_scaled_08AU': maxbxn08,'bx_i_max_scaled_09AU': maxbxn09,
'by_i_max_scaled_03AU': maxbyn03,'by_i_max_scaled_04AU': maxbyn04,'by_i_max_scaled_05AU': maxbyn05,'by_i_max_scaled_06AU': maxbyn06,'by_i_max_scaled_07AU': maxbyn07,'by_i_max_scaled_08AU': maxbyn08,'by_i_max_scaled_09AU': maxbyn09,
'bz_i_max_scaled_03AU': maxbzn03,'bz_i_max_scaled_04AU': maxbzn04,'bz_i_max_scaled_05AU': maxbzn05,'bz_i_max_scaled_06AU': maxbzn06,'bz_i_max_scaled_07AU': maxbzn07,'bz_i_max_scaled_08AU': maxbzn08,'bz_i_max_scaled_09AU': maxbzn09,
'btot_i_max_scaled_03AU': maxbtotn03,'btot_i_max_scaled_04AU': maxbtotn04,'btot_i_max_scaled_05AU': maxbtotn05,'btot_i_max_scaled_06AU': maxbtotn06,'btot_i_max_scaled_07AU': maxbtotn07,'btot_i_max_scaled_08AU': maxbtotn08,'btot_i_max_scaled_09AU': maxbtotn09,
'speed_i_scaled_03AU': speedni03,'speed_i_scaled_04AU': speedni04,'speed_i_scaled_05AU': speedni05,'speed_i_scaled_06AU': speedni06,'speed_i_scaled_07AU': speedni07,'speed_i_scaled_08AU': speedni08,'speed_i_scaled_09AU': speedni09,
'density_i_scaled_03AU': densityni03,'density_i_scaled_04AU': densityni04,'density_i_scaled_05AU': densityni05,'density_i_scaled_06AU': densityni06,'density_i_scaled_07AU': densityni07,'density_i_scaled_08AU': densityni08,'density_i_scaled_09AU': densityni09,
'pdyn_i_scaled_03AU': pdyn_n03,'pdyn_i_scaled_04AU': pdyn_n04,'pdyn_i_scaled_05AU': pdyn_n05,'pdyn_i_scaled_06AU': pdyn_n06,'pdyn_i_scaled_07AU': pdyn_n07,'pdyn_i_scaled_08AU': pdyn_n08,'pdyn_i_scaled_09AU': pdyn_n09,
'dst_i_mean_scaled_03AU': dstmean03,'dst_i_mean_scaled_04AU': dstmean04,'dst_i_mean_scaled_05AU': dstmean05,'dst_i_mean_scaled_06AU': dstmean06,'dst_i_mean_scaled_07AU': dstmean07,'dst_i_mean_scaled_08AU': dstmean08,'dst_i_mean_scaled_09AU': dstmean09,
'dst_i_max_scaled_03AU': dstmax03,'dst_i_max_scaled_04AU': dstmax04,'dst_i_max_scaled_05AU': dstmax05,'dst_i_max_scaled_06AU': dstmax06,'dst_i_max_scaled_07AU': dstmax07,'dst_i_max_scaled_08AU': dstmax08,'dst_i_max_scaled_09AU': dstmax09,
'st_ni03': st_n03,'st_ni04': st_n04,'st_ni05': st_n05,'st_ni06': st_n06,'st_ni07': st_n07,'st_ni08': st_n08,'st_ni09': st_n09,
'en_ni03': en_n03,'en_ni04': en_n04,'en_ni05': en_n05,'en_ni06': en_n06,'en_ni07': en_n07,'en_ni08': en_n08,'en_ni09': en_n09,
'en_nan03': en_nan03,'en_nan04': en_nan04,'en_nan05': en_nan05,'en_nan06': en_nan06,'en_nan07': en_nan07,'en_nan08': en_nan08,'en_nan09': en_nan09,
'start':start,'distance':r_i,'st_normal':st,'en_normal':en}
exofile = open(os.path.dirname(os.path.abspath("exoprogram.txt"))+"/exoprogram.txt", 'wb')     #enter your own file path
pickle.dump(dict, exofile)
exofile.close()


### create dictionary with scale_factors ###
dict3={'Bmean_3':bfacmean[0],'Bmean_4':bfacmean[1],'Bmean_5':bfacmean[2],'Bmean_6':bfacmean[3],'Bmean_':bfacmean[4],'Bmean_8':bfacmean[5],'Bmean_9':bfacmean[6],
	'Bmax3':bfac[0],'Bmax4':bfac[1],'Bmax5':bfac[2],'Bmax6':bfac[3],'Bmax7':bfac[4],'Bmax8':bfac[5],'Bmax9':bfac[6],
	'Duration3':facdur[0],'Duration4':facdur[1],'Duration5':facdur[2],'Duration6':facdur[3],'Duration7':facdur[4],'Duration8':facdur[5],'Duration9':facdur[6],
	'Density3':facdens[0],'Density4':facdens[1],'Density5':facdens[2],'Density6':facdens[3],'Density7':facdens[4],'Density8':facdens[5],'Density9':facdens[6],
	'Speed3':spd[0],'Speed4':spd[1],'Speed5':spd[2],'Speed6':spd[3],'Speed7':spd[4],'Speed8':spd[5],'Speed9':spd[6]}
exofile3 = open(os.path.dirname(os.path.abspath("scale_factors.txt"))+"/scale_factors.txt", 'wb')     #enter your own file path
pickle.dump(dict3, exofile3)
exofile3.close()

### create tabe with scale factors ### 
oweli=pd.DataFrame([["%.1f"%bfacmean[0],"%.1f"%bfacmean[1],"%.1f"%bfacmean[2],"%.1f"%bfacmean[3],"%.1f"%bfacmean[4],"%.1f"%bfacmean[5],"%.1f"%bfacmean[6]],
	["%.1f"%bfac[0],"%.1f"%bfac[1],"%.1f"%bfac[2],"%.1f"%bfac[3],"%.1f"%bfac[4],"%.1f"%bfac[5],"%.1f"%bfac[6]],
	["%.1f"%facdur[0],"%.1f"%facdur[1],"%.1f"%facdur[2],"%.1f"%facdur[3],"%.1f"%facdur[4],"%.1f"%facdur[5],"%.1f"%facdur[6]],
	["%.1f"%facdens[0],"%.1f"%facdens[1],"%.1f"%facdens[2],"%.1f"%facdens[3],"%.1f"%facdens[4],"%.1f"%facdens[5],"%.1f"%facdens[6]],
	["%.1f"%spd[0],"%.1f"%spd[1],"%.1f"%spd[2],"%.1f"%spd[3],"%.1f"%spd[4],"%.1f"%spd[5],"%.1f"%spd[6]]],
	index=['Bmean','Bmax','Duration','Density','Speed'],columns=['0.3AU','0.4AU','0.5AU','0.6AU','0.7AU','0.8AU','0.9AU'])

with open('factor_table.tex','w') as tf:
	tf.write(oweli.to_latex())
	
sys.exit()

