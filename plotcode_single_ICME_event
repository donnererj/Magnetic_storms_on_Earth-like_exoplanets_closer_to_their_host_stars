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
from matplotlib import rc
import scipy.io
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
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
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas as pd
from datetime import timedelta
from time import mktime 
from numpy import median
from scipy import stats

#open file 
file = open(os.path.dirname(os.path.abspath("exoprogram.txt"))+"/exoprogram.txt", 'rb')
dict = pickle.load(file)


#get current directory 
current_directory = os.getcwd()
output1= os.path.join(current_directory, r"plots_single_ICME_event")
if not os.path.exists(output1):
   os.makedirs(output1)

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/solar_cycle") == False: os.mkdir(output1+"/solar_cycle")
output1= os.path.abspath(output1+"/solar_cycle")

#create folders for plots 
if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/Bfield_comparison_solar") == False: os.mkdir(output1+"/Bfield_comparison_solar")
Bfield_comparison_solar= os.path.abspath(output1+"/Bfield_comparison_solar")

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/Dst_comparison_solar") == False: os.mkdir(output1+"/Dst_comparison_solar")
Dst_comparison_solar= os.path.abspath(output1+"/Dst_comparison_solar")

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/Observed_flux_solar") == False: os.mkdir(output1+"/Observed_flux_solar")
Observed_flux_solar= os.path.abspath(output1+"/Observed_flux_solar")

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/Scaled_fluxrope_solar") == False: os.mkdir(output1+"/Scaled_fluxrope_solar")
Scaled_fluxrope_solar= os.path.abspath(output1+"/Scaled_fluxrope_solar")



def startfunc(inputdate):
#define start and end dates for ICME event
	api=mdates.num2date(inputdate)
	obi0=[x.day==5 and x.month==7 and x.year==2013 for x in api]    #type in first day of minimum
	wobi0=np.where(obi0)
	start2008=wobi0[0][0]

	#last day/entry from 2008 
	obi=[x.day==7 and x.month==7 and x.year==2013 for x in api]		#type in last day of minimum
	wobi=np.where(obi)
	end2008=wobi[0][-1]

	#first day/entry from June 2013
	obi1=[x.day==16 and x.month==6 and x.year==2012 for x in api]	#type in first day of maximum
	wobi1=np.where(obi1)
	start2013=wobi1[0][0]

	#first day/entry from June 2014 
	obi2=[x.day==18 and x.month==6 and x.year==2012 for x in api]	#type in last day of maximum
	wobi2=np.where(obi2)
	end2014=wobi2[0][0]
	return (start2008, end2008, start2013, end2014)

start2008_03, end2008_03, start2013_03, end2014_03=startfunc(dict['time_i_scaled_03AU'])
start2008_04, end2008_04, start2013_04, end2014_04=startfunc(dict['time_i_scaled_04AU'])
start2008_05, end2008_05, start2013_05, end2014_05=startfunc(dict['time_i_scaled_05AU'])
start2008_06, end2008_06, start2013_06, end2014_06=startfunc(dict['time_i_scaled_06AU'])
start2008_07, end2008_07, start2013_07, end2014_07=startfunc(dict['time_i_scaled_07AU'])
start2008_08, end2008_08, start2013_08, end2014_08=startfunc(dict['time_i_scaled_08AU'])
start2008_09, end2008_09, start2013_09, end2014_09=startfunc(dict['time_i_scaled_09AU'])
start2008_1, end2008_1, start2013_1, end2014_1=startfunc(dict['time_normal'])


time_normal_date=mdates.num2date(dict['time_normal'])
start_n=time_normal_date[0].strftime('%Y-%m-%d')

#range for 0.3 AU until 0.9 AU 
ron=np.arange(3,10,1)



#create plot function 
def plotcomparison(start1,end1,start2,end2,start1n,end1n,start2n,end2n,timen,bxvalue,byvalue,bzvalue,btotvalue,dstvalue,dstobrienvalue,times,bxsvalue,bysvalue,bzsvalue,btotsvalue,dstscaledvalue,savevalue):
	
	stmax=mdates.num2date(timen[start2n]).strftime('%Y-%m-%d')
	stmin=mdates.num2date(timen[start1n]).strftime('%Y-%m-%d') 
	
	stmaxs=mdates.num2date(times[start2]).strftime('%Y-%m-%d')
	stmins=mdates.num2date(times[start1]).strftime('%Y-%m-%d') 
	
	###1### fluxrope comparison for 1 AU and chosen heliospheric distance - solar minimum 
	f1=plt.figure(41,figsize=(16,10))
	sns.set_style("whitegrid")  
	a1=f1.add_subplot(211)
	a1.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a1.transAxes,fontsize=17)
	#use timeaxis & B-data for 1 AU 
	plt.plot_date(timen[start1n:end1n],bxvalue[start1n:end1n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 		#data for Bx-component 
	plt.plot_date(timen[start1n:end1n],byvalue[start1n:end1n],'k',color='green',label=r'$\mathrm{By\,component}$')		#data for By-component 
	plt.plot_date(timen[start1n:end1n],bzvalue[start1n:end1n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')		#data for Bz-component 
	plt.plot_date(timen[start1n:end1n],btotvalue[start1n:end1n],'k',color='black',label=r'$\mathrm{B\,total}$')			#data for Btotal-component 
	#get date & plot it on the ticks 
	a1.xaxis_date()
	a1.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a1.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a1.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a1.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	#write text into plot
	plt.text(0.1,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a1.transAxes)
	plt.text(0.1,0.9, '$\mathrm{%s}$'%(stmin),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a1.transAxes)
	#make grey shade for duration of ICME event 
	for p in range(len(dict['st_normal'])):
			a1.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5, color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#select axis-range dependent on the chosen heliospheric distance 
	if savevalue==ron[0]:
		plt.ylim(-110,150)
	if savevalue==ron[1]:
		plt.ylim(-100,100)	
	if savevalue==ron[2]:
		plt.ylim(-50,60)
	if savevalue==ron[3]:
		plt.ylim(-50,50)
	if savevalue==ron[4]:
		plt.ylim(-30,40)	
	if savevalue==ron[5]:
		plt.ylim(-20,30)
	if savevalue==ron[6]:
		plt.ylim(-20,25)	
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 
	
	a2=f1.add_subplot(212)  
	a2.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a2.transAxes,fontsize=17)
	#use timeaxis & B-data for chosen heliospheric distance 
	line1=plt.plot_date(times[start1:end1],bxsvalue[start1:end1],'k',color='red',label=r'$\mathrm{Bx\,component}$')     
	line2=plt.plot_date(times[start1:end1],bysvalue[start1:end1],'k',color='green',label=r'$\mathrm{By\,component}$')
	line3=plt.plot_date(times[start1:end1],bzsvalue[start1:end1],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	line4=plt.plot_date(times[start1:end1],btotsvalue[start1:end1],'k',color='black',label=r'$\mathrm{B\,total}$')
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
			a2.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a2.xaxis_date()
	a2.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a2.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a2.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a2.xaxis.set_tick_params(which='minor', labelsize=16)
	plt.minorticks_on()
	plt.text(0.1,0.8,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=16,color='black',horizontalalignment='center',verticalalignment='center',transform = a2.transAxes)
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-110,150)
	if savevalue==ron[1]:
		plt.ylim(-100,100)	
	if savevalue==ron[2]:
		plt.ylim(-50,60)
	if savevalue==ron[3]:
		plt.ylim(-50,50)
	if savevalue==ron[4]:
		plt.ylim(-30,40)	
	if savevalue==ron[5]:
		plt.ylim(-20,30)
	if savevalue==ron[6]:
		plt.ylim(-20,25)
	plt.xlim(times[start1],times[end1])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 
	plt.savefig(Bfield_comparison_solar+"/fluxropecomparison_0.%dAU_solarmin.pdf"%(savevalue),dpi=300)
	plt.savefig(Bfield_comparison_solar+"/fluxropecomparison_0.%dAU_solarmin.png"%(savevalue),dpi=300)
	plt.close()
	
	
	###2### fluxrope comparison for 1 AU and chosen heliospheric distance - solar maximum 
	f2=plt.figure(42,figsize=(16,10))
	a3=f2.add_subplot(211)
	a3.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a3.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],bxvalue[start2n:end2n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start2n:end2n],byvalue[start2n:end2n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start2n:end2n],bzvalue[start2n:end2n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start2n:end2n],btotvalue[start2n:end2n],'k',color='black',label=r'$\mathrm{B\,total}$')
	a3.xaxis_date()
	a3.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a3.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a3.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a3.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a3.transAxes)
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmax),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a3.transAxes)
	for p in range(len(dict['st_normal'])):
			a3.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5, color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=18,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-190,350)
	if savevalue==ron[1]:
		plt.ylim(-100,200)	
	if savevalue==ron[2]:
		plt.ylim(-100,150)
	if savevalue==ron[3]:
		plt.ylim(-50,100)
	if savevalue==ron[4]:
		plt.ylim(-40,90)	
	if savevalue==ron[5]:
		plt.ylim(-40,60)
	if savevalue==ron[6]:
		plt.ylim(-20,50)
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a4=f2.add_subplot(212)  
	a4.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a4.transAxes,fontsize=17)
	line1=plt.plot_date(times[start2:end2],bxsvalue[start2:end2],'k',color='red',label=r'$\mathrm{Bx\,component}$')     
	line2=plt.plot_date(times[start2:end2],bysvalue[start2:end2],'k',color='green',label=r'$\mathrm{By\,component}$')
	line3=plt.plot_date(times[start2:end2],bzsvalue[start2:end2],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	line4=plt.plot_date(times[start2:end2],btotsvalue[start2:end2],'k',color='black',label=r'$\mathrm{B\,total}$')
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
			a4.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a4.xaxis_date()
	a4.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a4.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a4.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a4.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a4.transAxes)
	plt.ylabel('$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-190,350)
	if savevalue==ron[1]:
		plt.ylim(-100,200)	
	if savevalue==ron[2]:
		plt.ylim(-100,150)
	if savevalue==ron[3]:
		plt.ylim(-50,100)
	if savevalue==ron[4]:
		plt.ylim(-40,90)	
	if savevalue==ron[5]:
		plt.ylim(-40,60)
	if savevalue==ron[6]:
		plt.ylim(-20,50)
	plt.xlim(times[start2],times[end2])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 	

	plt.savefig(Bfield_comparison_solar+"/fluxropecomparison_0.%dAU_solarmax.pdf"%(savevalue),dpi=300)
	plt.savefig(Bfield_comparison_solar+"/fluxropecomparison_0.%dAU_solarmax.png"%(savevalue),dpi=300)
	plt.close()


    ###3### Dst-index comparison for 1 AU and chosen heliospheric distance - solar minimum 
	f3=plt.figure(43,figsize=(16,10))
	a5 = f3.add_subplot(211)
	a5.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a5.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],dstvalue[start1n:end1n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start1n:end1n],dstobrienvalue[start1n:end1n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a5.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a5.xaxis_date()
	a5.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a5.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a5.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a5.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a5.transAxes)
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmin),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a5.transAxes)
	plt.ylabel('$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')    ####obere zeile soll start auch als text wiedergegeben werden 
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=18,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-300,50)
	if savevalue==ron[1]:
		plt.ylim(-200,50)	
	if savevalue==ron[2]:
		plt.ylim(-150,50)
	if savevalue==ron[3]:
		plt.ylim(-100,40)
	if savevalue==ron[4]:
		plt.ylim(-80,40)	
	if savevalue==ron[5]:
		plt.ylim(-60,40)
	if savevalue==ron[6]:
		plt.ylim(-50,40)
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)

	a6=f3.add_subplot(212)
	#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
	a6.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a6.transAxes,fontsize=17)
	plt.plot_date(times[start1:end1],dstscaledvalue[start1:end1],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a6.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a6.xaxis_date()
	a6.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a6.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a6.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a6.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.text(0.9,0.9,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a6.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-300,50)
	if savevalue==ron[1]:
		plt.ylim(-200,50)	
	if savevalue==ron[2]:
		plt.ylim(-150,50)
	if savevalue==ron[3]:
		plt.ylim(-100,40)
	if savevalue==ron[4]:
		plt.ylim(-80,40)	
	if savevalue==ron[5]:
		plt.ylim(-60,40)
	if savevalue==ron[6]:
		plt.ylim(-50,40)
	plt.xlim(times[start1],times[end1])
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Dst_comparison_solar+"/Dst_comparison_0.%dAU_solarmin.pdf"%(savevalue),dpi=300)
	plt.savefig(Dst_comparison_solar+"/Dst_comparison_0.%dAU_solarmin.png"%(savevalue),dpi=300)
	plt.close()
	
	
	###4### Dst-index comparison for 1 AU and chosen heliospheric distance - solar maximum 
	f4=plt.figure(44,figsize=(16,10))
	a7 = f4.add_subplot(211)
	a7.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a7.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],dstvalue[start2n:end2n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start2n:end2n],dstobrienvalue[start2n:end2n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a7.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a7.xaxis_date()
	a7.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a7.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a7.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a7.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a7.transAxes)
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmax),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a7.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')    ####obere zeile soll start auch als text wiedergegeben werden 
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-900,200)
	if savevalue==ron[1]:
		plt.ylim(-600,180)	
	if savevalue==ron[2]:
		plt.ylim(-400,100)
	if savevalue==ron[3]:
		plt.ylim(-300,100)
	if savevalue==ron[4]:
		plt.ylim(-210,100)	
	if savevalue==ron[5]:
		plt.ylim(-200,100)
	if savevalue==ron[6]:
		plt.ylim(-150,100)
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	
	a8=f4.add_subplot(212)
	#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
	a8.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a8.transAxes,fontsize=17)
	plt.plot_date(times[start2:end2],dstscaledvalue[start2:end2],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a8.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a8.xaxis_date()
	a8.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a8.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a8.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a8.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.text(0.9,0.9,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a8.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-900,200)
	if savevalue==ron[1]:
		plt.ylim(-600,180)	
	if savevalue==ron[2]:
		plt.ylim(-400,100)
	if savevalue==ron[3]:
		plt.ylim(-300,100)
	if savevalue==ron[4]:
		plt.ylim(-210,100)	
	if savevalue==ron[5]:
		plt.ylim(-200,100)
	if savevalue==ron[6]:
		plt.ylim(-150,100)
	plt.xlim(times[start2],times[end2])
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Dst_comparison_solar+"/Dst_comparison_0.%dAU_solarmax.pdf"%(savevalue),dpi=300)
	plt.savefig(Dst_comparison_solar+"/Dst_comparison_0.%dAU_solarmax.png"%(savevalue),dpi=300)
	plt.close()
	
	
	###5### observed fluxrope adjusted to specific heliospheric distance for solar maximum 
	f5=plt.figure(45,figsize=(16,10))
	sns.set_style("whitegrid")  
	a9=f5.add_subplot(211)
	a9.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a9.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],bxvalue[start2n:end2n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start2n:end2n],byvalue[start2n:end2n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start2n:end2n],bzvalue[start2n:end2n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start2n:end2n],btotvalue[start2n:end2n],'k',color='black',label=r'$\mathrm{B\,total}$')	
	a9.xaxis_date()
	a9.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a9.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a9.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a9.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a9.transAxes)
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmax),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a9.transAxes)
	#plt.text(0.9,0.9,start_n,fontsize=14,color='black',horizontalalignment='center',verticalalignment='center',transform = a9.transAxes)
	for p in range(len(dict['st_normal'])):
		a9.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5, color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-190,350)
	if savevalue==ron[1]:
		plt.ylim(-100,200)	
	if savevalue==ron[2]:
		plt.ylim(-100,150)
	if savevalue==ron[3]:
		plt.ylim(-50,100)
	if savevalue==ron[4]:
		plt.ylim(-40,90)	
	if savevalue==ron[5]:
		plt.ylim(-40,60)
	if savevalue==ron[6]:
		plt.ylim(-20,50)
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a10 = f5.add_subplot(212)
	a10.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a10.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],dstvalue[start2n:end2n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start2n:end2n],dstobrienvalue[start2n:end2n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a10.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a10.xaxis_date()
	a10.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a10.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a10.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a10.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(1.18,0.16,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a10.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')     
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-900,200)
	if savevalue==ron[1]:
		plt.ylim(-600,180)	
	if savevalue==ron[2]:
		plt.ylim(-400,100)
	if savevalue==ron[3]:
		plt.ylim(-300,100)
	if savevalue==ron[4]:
		plt.ylim(-210,100)	
	if savevalue==ron[5]:
		plt.ylim(-200,100)
	if savevalue==ron[6]:
		plt.ylim(-150,100)
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Observed_flux_solar+"/obs_adjusted_1AU_%s_0.%dAU_solarmax.pdf"%(start_n,savevalue),dpi=300)	
	plt.savefig(Observed_flux_solar+"/obs_adjusted_1AU_%s_0.%dAU_solarmax.png"%(start_n,savevalue),dpi=300)	
	plt.close()
	
	
	###6### observed fluxrope adjusted to specific heliospheric distance for solar minimum 
	f6=plt.figure(46,figsize=(16,10))
	sns.set_style("whitegrid")  
	a11=f6.add_subplot(211)
	a11.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a11.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],bxvalue[start1n:end1n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start1n:end1n],byvalue[start1n:end1n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start1n:end1n],bzvalue[start1n:end1n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start1n:end1n],btotvalue[start1n:end1n],'k',color='black',label=r'$\mathrm{B\,total}$')	
	a11.xaxis_date()
	a11.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a11.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a11.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a11.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a11.transAxes)
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmin),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a11.transAxes)
	for p in range(len(dict['st_normal'])):
		a11.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5, color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=18,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-110,150)
	if savevalue==ron[1]:
		plt.ylim(-100,100)	
	if savevalue==ron[2]:
		plt.ylim(-50,60)
	if savevalue==ron[3]:
		plt.ylim(-50,50)
	if savevalue==ron[4]:
		plt.ylim(-30,40)	
	if savevalue==ron[5]:
		plt.ylim(-20,30)
	if savevalue==ron[6]:
		plt.ylim(-20,25)
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a12 = f6.add_subplot(212)
	a12.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a12.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],dstvalue[start1n:end1n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start1n:end1n],dstobrienvalue[start1n:end1n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a12.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a12.xaxis_date()
	a12.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a12.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a12.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a12.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(1.18,0.16,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a12.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')     
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-300,50)
	if savevalue==ron[1]:
		plt.ylim(-200,50)	
	if savevalue==ron[2]:
		plt.ylim(-150,50)
	if savevalue==ron[3]:
		plt.ylim(-100,40)
	if savevalue==ron[4]:
		plt.ylim(-80,40)	
	if savevalue==ron[5]:
		plt.ylim(-60,40)
	if savevalue==ron[6]:
		plt.ylim(-50,40)
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.yticks(fontsize=16) 
	plt.xticks(fontsize=16)
	plt.tight_layout()
	plt.savefig(Observed_flux_solar+"/obs_adjusted_1AU_%s_0.%dAU_solarmin.pdf"%(start_n,savevalue),dpi=300)	
	plt.savefig(Observed_flux_solar+"/obs_adjusted_1AU_%s_0.%dAU_solarmin.png"%(start_n,savevalue),dpi=300)	
	plt.close()
	
	
	###7### observed fluxrope at 1 AU for solar maximum 
	f51=plt.figure(49,figsize=(16,10))
	sns.set_style("whitegrid")  
	a91=f51.add_subplot(211)
	a91.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a91.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],bxvalue[start2n:end2n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start2n:end2n],byvalue[start2n:end2n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start2n:end2n],bzvalue[start2n:end2n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start2n:end2n],btotvalue[start2n:end2n],'k',color='black',label=r'$\mathrm{B\,total}$')	
	a91.xaxis_date()
	a91.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a91.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a91.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a91.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a91.transAxes)
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmax),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a91.transAxes)
	#plt.text(0.9,0.9,start_n,fontsize=14,color='black',horizontalalignment='center',verticalalignment='center',transform = a9.transAxes)
	for p in range(len(dict['st_normal'])):
		a91.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5, color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=18,fontweight='bold')
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a101 = f51.add_subplot(212)
	a101.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a101.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],dstvalue[start2n:end2n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start2n:end2n],dstobrienvalue[start2n:end2n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a101.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a101.xaxis_date()
	a101.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a101.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a101.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a101.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(1.18,0.16,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a101.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')     
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Observed_flux_solar+"/observedfluxrope_1AU_%s_1AU_solarmax.pdf"%(start_n),dpi=300)	
	plt.savefig(Observed_flux_solar+"/observedfluxrope_1AU_%s_1AU_solarmax.png"%(start_n),dpi=300)	
	plt.close()
	
	
	###8### observed fluxrope at 1 AU for solar minimum 
	f16=plt.figure(50,figsize=(16,10))
	sns.set_style("whitegrid")  
	a111=f16.add_subplot(211)
	a111.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a111.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],bxvalue[start1n:end1n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start1n:end1n],byvalue[start1n:end1n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start1n:end1n],bzvalue[start1n:end1n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start1n:end1n],btotvalue[start1n:end1n],'k',color='black',label=r'$\mathrm{B\,total}$')	
	a111.xaxis_date()
	a111.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a111.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a111.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a111.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a111.transAxes)
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmin),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a111.transAxes)
	for p in range(len(dict['st_normal'])):
		a111.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5, color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a121 = f16.add_subplot(212)
	a121.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a121.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],dstvalue[start1n:end1n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start1n:end1n],dstobrienvalue[start1n:end1n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a121.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a121.xaxis_date()
	a121.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a121.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a121.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a121.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(1.18,0.16,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a121.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')     
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Observed_flux_solar+"/observedfluxrope_1AU_%s_1AU_solarmin.pdf"%(start_n),dpi=300)	
	plt.savefig(Observed_flux_solar+"/observedfluxrope_1AU_%s_1AU_solarmin.png"%(start_n),dpi=300)	
	plt.close()
	
	
	###9### fluxrope & Dst-index scaled to chosen heliospheric distance - solar maximum 
	f7=plt.figure(47,figsize=(16,10)) 
	a13=f7.add_subplot(211)  
	a13.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a13.transAxes,fontsize=17)
	line1=plt.plot_date(times[start2:end2],bxsvalue[start2:end2],'k',color='red',label=r'$\mathrm{Bx\,component}$')     
	line2=plt.plot_date(times[start2:end2],bysvalue[start2:end2],'k',color='green',label=r'$\mathrm{By\,component}$')
	line3=plt.plot_date(times[start2:end2],bzsvalue[start2:end2],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	line4=plt.plot_date(times[start2:end2],btotsvalue[start2:end2],'k',color='black',label=r'$\mathrm{B\,total}$')
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
			a13.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a13.xaxis_date()
	a13.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a13.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a13.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a13.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmaxs),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a13.transAxes)
	plt.text(0.9,0.8,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a13.transAxes)
	plt.ylabel('$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-190,350)
	if savevalue==ron[1]:
		plt.ylim(-100,200)	
	if savevalue==ron[2]:
		plt.ylim(-100,150)
	if savevalue==ron[3]:
		plt.ylim(-50,100)
	if savevalue==ron[4]:
		plt.ylim(-40,90)	
	if savevalue==ron[5]:
		plt.ylim(-40,60)
	if savevalue==ron[6]:
		plt.ylim(-20,50)
	plt.xlim(times[start2],times[end2])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 	
	
	
	a14=f7.add_subplot(212)
	a14.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a14.transAxes,fontsize=17)
	plt.plot_date(times[start2:end2],dstscaledvalue[start2:end2],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a14.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a14.xaxis_date()
	a14.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a14.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a14.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a14.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-900,200)
	if savevalue==ron[1]:
		plt.ylim(-600,180)	
	if savevalue==ron[2]:
		plt.ylim(-400,100)
	if savevalue==ron[3]:
		plt.ylim(-300,100)
	if savevalue==ron[4]:
		plt.ylim(-210,100)	
	if savevalue==ron[5]:
		plt.ylim(-200,100)
	if savevalue==ron[6]:
		plt.ylim(-150,100)
	plt.xlim(times[start2],times[end2])
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Scaled_fluxrope_solar+"/scaled_fluxrope_0.%dAU_solarmax.pdf"%(savevalue),dpi=300)
	plt.savefig(Scaled_fluxrope_solar+"/scaled_fluxrope_0.%dAU_solarmax.png"%(savevalue),dpi=300)
	plt.close()
	
	
	###10### fluxrope & Dst-index scaled to chosen heliospheric distance - solar minimum 
	f8=plt.figure(48,figsize=(16,10)) 
	a15=f8.add_subplot(211)  
	a15.text(-0.07, 1.0,'$\mathrm{a)}$',transform=a15.transAxes,fontsize=17)
	line1=plt.plot_date(times[start1:end1],bxsvalue[start1:end1],'k',color='red',label=r'$\mathrm{Bx\,component}$')     
	line2=plt.plot_date(times[start1:end1],bysvalue[start1:end1],'k',color='green',label=r'$\mathrm{By\,component}$')
	line3=plt.plot_date(times[start1:end1],bzsvalue[start1:end1],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	line4=plt.plot_date(times[start1:end1],btotsvalue[start1:end1],'k',color='black',label=r'$\mathrm{B\,total}$')
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
			a15.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a15.xaxis_date()
	a15.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a15.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a15.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a15.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmins),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a15.transAxes)
	plt.text(0.9,0.8,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a15.transAxes)
	plt.ylabel('$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-110,150)
	if savevalue==ron[1]:
		plt.ylim(-100,100)	
	if savevalue==ron[2]:
		plt.ylim(-50,60)
	if savevalue==ron[3]:
		plt.ylim(-50,50)
	if savevalue==ron[4]:
		plt.ylim(-30,40)	
	if savevalue==ron[5]:
		plt.ylim(-20,30)
	if savevalue==ron[6]:
		plt.ylim(-20,25)
	plt.xlim(times[start1],times[end1])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 	
	
	
	a16=f8.add_subplot(212)
	a16.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a16.transAxes,fontsize=17)
	plt.plot_date(times[start1:end1],dstscaledvalue[start1:end1],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a16.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a16.xaxis_date()
	a16.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))   	# normal '%H:%M'
	a16.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a16.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a16.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=15)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-300,50)
	if savevalue==ron[1]:
		plt.ylim(-200,50)	
	if savevalue==ron[2]:
		plt.ylim(-150,50)
	if savevalue==ron[3]:
		plt.ylim(-100,40)
	if savevalue==ron[4]:
		plt.ylim(-80,40)	
	if savevalue==ron[5]:
		plt.ylim(-60,40)
	if savevalue==ron[6]:
		plt.ylim(-50,40)
	plt.xlim(times[start1],times[end1])
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Scaled_fluxrope_solar+"/scaled_fluxrope_0.%dAU_solarmin.pdf"%(savevalue),dpi=300)
	plt.savefig(Scaled_fluxrope_solar+"/scaled_fluxrope_0.%dAU_solarmin.png"%(savevalue),dpi=300)
	plt.close()
	
	
	return (f1,f2,f3,f4,f5,f6,f7,f8,f16,f51)

#make scaled plots from 0.3 AU until 0.9 AU	
#start1,end1,start2,end2,start1n,end1n,start2n,end2n,timen,bxvalue,byvalue,bzvalue,btotvalue,dstvalue,dstobrienvalue,times,bxsvalue,bysvalue,bzsvalue,btotsvalue,dstscaledvalue,savevalue		
plotcomparison(start2008_03,end2008_03,start2013_03,end2014_03,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_03AU'],dict['bx_i_mean_scaled_03AU'],dict['by_i_mean_scaled_03AU'],dict['bz_i_mean_scaled_03AU'],dict['btot_i_mean_scaled_03AU'],dict['dst_i_mean_scaled_03AU'],ron[0])		
plotcomparison(start2008_04,end2008_04,start2013_04,end2014_04,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_04AU'],dict['bx_i_mean_scaled_04AU'],dict['by_i_mean_scaled_04AU'],dict['bz_i_mean_scaled_04AU'],dict['btot_i_mean_scaled_04AU'],dict['dst_i_mean_scaled_04AU'],ron[1])		
plotcomparison(start2008_05,end2008_05,start2013_05,end2014_05,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_05AU'],dict['bx_i_mean_scaled_05AU'],dict['by_i_mean_scaled_05AU'],dict['bz_i_mean_scaled_05AU'],dict['btot_i_mean_scaled_05AU'],dict['dst_i_mean_scaled_05AU'],ron[2])		
plotcomparison(start2008_06,end2008_06,start2013_06,end2014_06,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_06AU'],dict['bx_i_mean_scaled_06AU'],dict['by_i_mean_scaled_06AU'],dict['bz_i_mean_scaled_06AU'],dict['btot_i_mean_scaled_06AU'],dict['dst_i_mean_scaled_06AU'],ron[3])		
plotcomparison(start2008_07,end2008_07,start2013_07,end2014_07,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_07AU'],dict['bx_i_mean_scaled_07AU'],dict['by_i_mean_scaled_07AU'],dict['bz_i_mean_scaled_07AU'],dict['btot_i_mean_scaled_07AU'],dict['dst_i_mean_scaled_07AU'],ron[4])		
plotcomparison(start2008_08,end2008_08,start2013_08,end2014_08,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_08AU'],dict['bx_i_mean_scaled_08AU'],dict['by_i_mean_scaled_08AU'],dict['bz_i_mean_scaled_08AU'],dict['btot_i_mean_scaled_08AU'],dict['dst_i_mean_scaled_08AU'],ron[5])		
plotcomparison(start2008_09,end2008_09,start2013_09,end2014_09,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_09AU'],dict['bx_i_mean_scaled_09AU'],dict['by_i_mean_scaled_09AU'],dict['bz_i_mean_scaled_09AU'],dict['btot_i_mean_scaled_09AU'],dict['dst_i_mean_scaled_09AU'],ron[6])		


##get maximum & minimum values during ICME event for chosen heliospheric distances:
##1AU
#define start and end date of ICME event
al=np.where([x for x in dict['st_normal']]>start2013_1)	
bl=np.where([x for x in dict['st_normal']]<end2014_1)	
aku=np.where(np.in1d(bl,al))							
allo=int(dict['st_normal'][aku])						
ab=np.where([x for x in dict['en_normal']]>start2013_1)
bb=np.where([x for x in dict['en_normal']]<end2014_1)
abd=np.where(np.in1d(bb,ab))
ollo=int(dict['en_normal'][abd])

#get maximum/minimum data during ICME event
toi1=max(dict['btot_normal'][allo:ollo])		 		#Maximum Bfield
toi2=min(dict['btot_normal'][allo:ollo])				#Minimum Bfield		
toi3=min(dict['dst_normal'][allo:ollo])					#Minimum Dst
toi5=max(dict['density_normal'][allo:ollo])				#Maximum Density  
toi5a=min(dict['density_normal'][allo:ollo])			#Minimum Density 
toi6=max(dict['speed_normal'][allo:ollo])				#Maximum Speed 

#get duration of ICME event 
mori1=mdates.num2date(dict['time_normal'][allo])		#start date 
mori21=mdates.num2date(dict['time_normal'][ollo])		#end date 
kata1=mori21-mori1										#time between start & end 
#kata1/3600												#datetime object 
secs1 = kata1.total_seconds()							#get duration in seconds 
toi4 = int(secs1 / 3600)								#duration of ICME event at 1 AU in hours 


###same procedure for chosen heliospheric distances:
##0.3AU
n3=np.where([x for x in dict['st_ni03']]>start2013_03)
m3=np.where([x for x in dict['st_ni03']]<end2014_03)
mn3=np.where(np.in1d(m3,n3))
mollo=int(dict['st_ni03'][mn3])
bn3=np.where([x for x in dict['en_ni03']]>start2013_03)
bm3=np.where([x for x in dict['en_ni03']]<end2014_03)
apf=np.where(np.in1d(bm3,bn3))
was=int(dict['en_ni03'][apf])

tois1=max(dict['btot_i_mean_scaled_03AU'][mollo:was])		 
tois2=min(dict['btot_i_mean_scaled_03AU'][mollo:was])		
tois3=min(dict['dst_i_mean_scaled_03AU'][mollo:was])					
tois5=max(dict['density_i_scaled_03AU'][mollo:was])					  
tois5a=min(dict['density_i_scaled_03AU'][mollo:was])	
tois6=max(dict['speed_i_scaled_03AU'][mollo:was])

mori=mdates.num2date(dict['time_i_scaled_03AU'][mollo])
mori2=mdates.num2date(dict['time_i_scaled_03AU'][was])
kata=mori2-mori
#kata/3600
secs = kata.total_seconds()
tois4 = int(secs / 3600)


##0.4 AU
n34=np.where([x for x in dict['st_ni04']]>start2013_04)
m34=np.where([x for x in dict['st_ni04']]<end2014_04)
mn34=np.where(np.in1d(m34,n34))
mollo4=int(dict['st_ni04'][mn34])
bn34=np.where([x for x in dict['en_ni04']]>start2013_04)
bm34=np.where([x for x in dict['en_ni04']]<end2014_04)
apf4=np.where(np.in1d(bm34,bn34))
was4=int(dict['en_ni04'][apf4])
 
toisj1=max(dict['btot_i_mean_scaled_04AU'][mollo4:was4])		 			
toisj2=min(dict['btot_i_mean_scaled_04AU'][mollo4:was4])						
toisj3=min(dict['dst_i_mean_scaled_04AU'][mollo4:was4])					
toisj5=max(dict['density_i_scaled_04AU'][mollo4:was4])					 
toisj5a=min(dict['density_i_scaled_04AU'][mollo4:was4])
toisj6=max(dict['speed_i_scaled_04AU'][mollo4:was4])

mori4=mdates.num2date(dict['time_i_scaled_04AU'][mollo4])
mori24=mdates.num2date(dict['time_i_scaled_04AU'][was4])
kata4=mori24-mori4
#kata4/3600
secs4 = kata4.total_seconds()
toisj4 = int(secs4 / 3600)


##0.5AU
n35=np.where([x for x in dict['st_ni05']]>start2013_05)
m35=np.where([x for x in dict['st_ni05']]<end2014_05)
mn35=np.where(np.in1d(m35,n35))
mollo5=int(dict['st_ni05'][mn35])
bn35=np.where([x for x in dict['en_ni05']]>start2013_05)
bm35=np.where([x for x in dict['en_ni05']]<end2014_05)
apf5=np.where(np.in1d(bm35,bn35))
was5=int(dict['en_ni05'][apf5])
 
toisw1=max(dict['btot_i_mean_scaled_05AU'][mollo5:was5])		 			
toisw2=min(dict['btot_i_mean_scaled_05AU'][mollo5:was5])							
toisw3=min(dict['dst_i_mean_scaled_05AU'][mollo5:was5])					 
toisw5=max(dict['density_i_scaled_05AU'][mollo5:was5])					 
toisw5a=min(dict['density_i_scaled_05AU'][mollo5:was5])	
toisw6=max(dict['speed_i_scaled_05AU'][mollo5:was5])

mori5=mdates.num2date(dict['time_i_scaled_05AU'][mollo5])
mori25=mdates.num2date(dict['time_i_scaled_05AU'][was5])
kata5=mori25-mori5
#kata5/3600
secs5 = kata5.total_seconds()
toisw4 = int(secs5 / 3600)


##0.6AU
n36=np.where([x for x in dict['st_ni06']]>start2013_06)
m36=np.where([x for x in dict['st_ni06']]<end2014_06)
mn36=np.where(np.in1d(m36,n36))
mollo6=int(dict['st_ni06'][mn36])
bn36=np.where([x for x in dict['en_ni06']]>start2013_06)
bm36=np.where([x for x in dict['en_ni06']]<end2014_06)
apf6=np.where(np.in1d(bm36,bn36))
was6=int(dict['en_ni06'][apf6])

toisc1=max(dict['btot_i_mean_scaled_06AU'][mollo6:was6])		 			
toisc2=min(dict['btot_i_mean_scaled_06AU'][mollo6:was6])							
toisc3=min(dict['dst_i_mean_scaled_06AU'][mollo6:was6])					
toisc5=max(dict['density_i_scaled_06AU'][mollo6:was6])					 
toisc5a=min(dict['density_i_scaled_06AU'][mollo6:was6])	
toisc6=max(dict['speed_i_scaled_06AU'][mollo6:was6])

mori6=mdates.num2date(dict['time_i_scaled_06AU'][mollo6])
mori26=mdates.num2date(dict['time_i_scaled_06AU'][was6])
kata6=mori26-mori6
#kata6/3600
secs6 = kata6.total_seconds()
toisc4 = int(secs6 / 3600)


##0.7AU 
n37=np.where([x for x in dict['st_ni07']]>start2013_07)
m37=np.where([x for x in dict['st_ni07']]<end2014_07)
mn37=np.where(np.in1d(m37,n37))
mollo7=int(dict['st_ni07'][mn37])
bn37=np.where([x for x in dict['en_ni07']]>start2013_07)
bm37=np.where([x for x in dict['en_ni07']]<end2014_07)
apf7=np.where(np.in1d(bm37,bn37))
was7=int(dict['en_ni07'][apf7])

toise1=max(dict['btot_i_mean_scaled_07AU'][mollo7:was7])		 		
toise2=min(dict['btot_i_mean_scaled_07AU'][mollo7:was7])							
toise3=min(dict['dst_i_mean_scaled_07AU'][mollo7:was7])					
toise5=max(dict['density_i_scaled_07AU'][mollo7:was7])					  
toise5a=min(dict['density_i_scaled_07AU'][mollo7:was7])
toise6=max(dict['speed_i_scaled_07AU'][mollo7:was7])

mori7=mdates.num2date(dict['time_i_scaled_07AU'][mollo7])
mori27=mdates.num2date(dict['time_i_scaled_07AU'][was7])
kata7=mori27-mori7
#kata7/3600
secs7 = kata7.total_seconds()
toise4 = int(secs7 / 3600)


##0.8AU
n38=np.where([x for x in dict['st_ni08']]>start2013_08)
m38=np.where([x for x in dict['st_ni08']]<end2014_08)
mn38=np.where(np.in1d(m38,n38))
mollo8=int(dict['st_ni08'][mn38])
bn38=np.where([x for x in dict['en_ni08']]>start2013_08)
bm38=np.where([x for x in dict['en_ni08']]<end2014_08)
apf8=np.where(np.in1d(bm38,bn38))
was8=int(dict['en_ni08'][apf8])

toisf1=max(dict['btot_i_mean_scaled_08AU'][mollo8:was8])		 		
toisf2=min(dict['btot_i_mean_scaled_08AU'][mollo8:was8])					
toisf3=min(dict['dst_i_mean_scaled_08AU'][mollo8:was8])						
toisf5=max(dict['density_i_scaled_08AU'][mollo8:was8])					
toisf5a=min(dict['density_i_scaled_08AU'][mollo8:was8])
toisf6=max(dict['speed_i_scaled_08AU'][mollo8:was8])

mori8=mdates.num2date(dict['time_i_scaled_08AU'][mollo8])
mori28=mdates.num2date(dict['time_i_scaled_08AU'][was8])
kata8=mori28-mori8
#kata8/3600
secs8 = kata8.total_seconds()
toisf4 = int(secs8 / 3600)


##0.9 AU
n39=np.where([x for x in dict['st_ni09']]>start2013_09)
m39=np.where([x for x in dict['st_ni09']]<end2014_09)
mn39=np.where(np.in1d(m39,n39))
mollo9=int(dict['st_ni09'][mn39])
bn39=np.where([x for x in dict['en_ni09']]>start2013_09)
bm39=np.where([x for x in dict['en_ni09']]<end2014_09)
apf9=np.where(np.in1d(bm39,bn39))
was9=int(dict['en_ni09'][apf9])
# tui9=mdates.num2date(dict['time_i_scaled_09AU'][mollo9])
# zuik9=mdates.num2date(dict['time_i_scaled_09AU'][was9])
# adf9=time.mktime(zuik9.timetuple()) - time.mktime(tui9.timetuple())  #get seconds between start & end 

toisr1=max(dict['btot_i_mean_scaled_09AU'][mollo9:was9])		 			
toisr2=min(dict['btot_i_mean_scaled_09AU'][mollo9:was9])						
toisr3=min(dict['dst_i_mean_scaled_09AU'][mollo9:was9])					
toisr5=max(dict['density_i_scaled_09AU'][mollo9:was9])					
toisr5a=min(dict['density_i_scaled_09AU'][mollo9:was9])	
toisr6=max(dict['speed_i_scaled_09AU'][mollo9:was9])

mori9=mdates.num2date(dict['time_i_scaled_09AU'][mollo9])
mori29=mdates.num2date(dict['time_i_scaled_09AU'][was9])
kata9=mori29-mori9
#kata9/3600
secs9 = kata9.total_seconds()
toisr4 = int(secs9 / 3600)



table2=pd.DataFrame([["%.1f"%tois1,"%.1f"%toisj1,"%.1f"%toisw1,"%.1f"%toisc1,"%.1f"%toise1,"%.1f"%toisf1,"%.1f"%toisr1,"%.1f"%toi1],
	["%.1f"%tois3,"%.1f"%toisj3,"%.1f"%toisw3,"%.1f"%toisc3,"%.1f"%toise3,"%.1f"%toisf3,"%.1f"%toisr3,"%.1f"%toi3],
	["%.1f"%tois4,"%.1f"%toisj4,"%.1f"%toisw4,"%.1f"%toisc4,"%.1f"%toise4,"%.1f"%toisf4,"%.1f"%toisr4,"%.1f"%toi4],
	["%.1f"%tois5,"%.1f"%toisj5,"%.1f"%toisw5,"%.1f"%toisc5,"%.1f"%toise5,"%.1f"%toisf5,"%.1f"%toisr5,"%.1f"%toi5],
	["%.1f"%tois6,"%.1f"%toisj6,"%.1f"%toisw6,"%.1f"%toisc6,"%.1f"%toise6,"%.1f"%toisf6,"%.1f"%toisr6,"%.1f"%toi6]],
	index=['Max. B-field [nT]','Min. Dst-index [nT]','Duration [hours]','Max. Density [N cm**-3]','Max. Speed [km s**-1]'],columns=['0.3AU','0.4AU','0.5AU','0.6AU','0.7AU','0.8AU','0.9AU','1AU'])

with open('examples_for_minimum_and_maximum_scaled_values.tex','w') as tf:
	tf.write(table2.to_latex())	


#minimum Dst-index at chosen heliospheric distances compared to 1 AU 	
ap3=tois3/toi3
ap4=toisj3/toi3
ap5=toisw3/toi3
ap6=toisc3/toi3
ap7=toise3/toi3
ap8=toisf3/toi3
ap9=toisr3/toi3

dstftex2=pd.DataFrame([["%.1f"%ap3,"%.1f"%ap4,"%.1f"%ap5,"%.1f"%ap6,"%.1f"%ap7,"%.1f"%ap8,"%.1f"%ap9]],
	index=['Dst-index'],columns=['0.3AU','0.4AU','0.5AU','0.6AU','0.7AU','0.8AU','0.9AU'])

with open('dst_neu_fac.tex','w') as tf:
	tf.write(dstftex2.to_latex())






sys.exit()