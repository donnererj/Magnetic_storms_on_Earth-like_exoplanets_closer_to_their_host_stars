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
from astropy.table import Table, Column
import pdb
import pandas as pd

#open file
file = open(os.path.dirname(os.path.abspath("exoprogram.txt"))+"/exoprogram.txt", 'rb')
dict = pickle.load(file)

#get current directory
current_directory = os.getcwd()
output1= os.path.join(current_directory, r"plots_timeseries_histograms_for_ICME_events")
if not os.path.exists(output1):
   os.makedirs(output1)
   
#create directories
if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/solar_cycle") == False: os.mkdir(output1+"/solar_cycle")
output1= os.path.abspath(output1+"/solar_cycle")

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/histogram_solar") == False: os.mkdir(output1+"/histogram_solar")
histo_solar= os.path.abspath(output1+"/histogram_solar")

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/histofit_solar") == False: os.mkdir(output1+"/histofit_solar")
histof_solar= os.path.abspath(output1+"/histofit_solar")

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/fits_solar") == False: os.mkdir(output1+"/fits_solar")
fits_solar= os.path.abspath(output1+"/fits_solar")

if os.path.isdir(output1) == False: os.mkdir(output1)
if os.path.isdir(output1+"/fits_histo") == False: os.mkdir(output1+"/fits_histo")
fits_histo= os.path.abspath(output1+"/fits_histo")


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

plt.close()



def startfunc(inputdate):
#define start and end dates for icme events
	api=mdates.num2date(inputdate)
	obi0=[x.day==1 and x.month==1 and x.year==2008 for x in api]    #first day of minimum
	wobi0=np.where(obi0)
	start2008=wobi0[0][0]

	#last day/entry from 2008 
	obi=[x.day==31 and x.month==12 and x.year==2008 for x in api]	#last day of minimum
	wobi=np.where(obi)
	end2008=wobi[0][-1]

	#first day/entry from June 2013
	obi1=[x.day==1 and x.month==6 and x.year==2013 for x in api]	#first day of maximum
	wobi1=np.where(obi1)
	start2013=wobi1[0][0]

	#first day/entry from June 2014 
	obi2=[x.day==1 and x.month==6 and x.year==2014 for x in api]	#last day of maximum
	wobi2=np.where(obi2)
	end2014=wobi2[0][0]
	return (start2008, end2008, start2013, end2014)

#define start and end dates for time range according to chosen heliospheric distance
start2008_03, end2008_03, start2013_03, end2014_03=startfunc(dict['time_i_scaled_03AU'])
start2008_04, end2008_04, start2013_04, end2014_04=startfunc(dict['time_i_scaled_04AU'])
start2008_05, end2008_05, start2013_05, end2014_05=startfunc(dict['time_i_scaled_05AU'])
start2008_06, end2008_06, start2013_06, end2014_06=startfunc(dict['time_i_scaled_06AU'])
start2008_07, end2008_07, start2013_07, end2014_07=startfunc(dict['time_i_scaled_07AU'])
start2008_08, end2008_08, start2013_08, end2014_08=startfunc(dict['time_i_scaled_08AU'])
start2008_09, end2008_09, start2013_09, end2014_09=startfunc(dict['time_i_scaled_09AU'])
start2008_1, end2008_1, start2013_1, end2014_1=startfunc(dict['time_normal'])
#hereby solar minimum = year 2008 
#		solar maximum = 1. June 2013 until 1. June 2014

ron=np.arange(3,10,1)


time_normal_date=mdates.num2date(dict['time_normal'])
start_n=time_normal_date[0].strftime('%Y-%m-%d')



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
	plt.plot_date(timen[start1n:end1n],bxvalue[start1n:end1n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start1n:end1n],byvalue[start1n:end1n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start1n:end1n],bzvalue[start1n:end1n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start1n:end1n],btotvalue[start1n:end1n],'k',color='black',label=r'$\mathrm{B\,total}$')
	a1.xaxis_date()
	a1.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a1.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a1.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a1.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.1,0.8,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a1.transAxes)
	plt.text(0.1,0.9, '$\mathrm{%s}$'%(stmin),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a1.transAxes)
	for p in range(len(dict['st_normal'])):
			a1.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5, color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-250,250)
	if savevalue==ron[1]:
		plt.ylim(-150,150)	
	if savevalue==ron[2]:
		plt.ylim(-100,100)
	if savevalue==ron[3]:
		plt.ylim(-70,70)
	if savevalue==ron[4]:
		plt.ylim(-60,60)	
	if savevalue==ron[5]:
		plt.ylim(-50,50)
	if savevalue==ron[6]:
		plt.ylim(-30,30)	
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 
	
	a2=f1.add_subplot(212)  
	a2.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a2.transAxes,fontsize=17)
	line1=plt.plot_date(times[start1:end1],bxsvalue[start1:end1],'k',color='red',label=r'$\mathrm{Bx\,component}$')     
	line2=plt.plot_date(times[start1:end1],bysvalue[start1:end1],'k',color='green',label=r'$\mathrm{By\,component}$')
	line3=plt.plot_date(times[start1:end1],bzsvalue[start1:end1],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	line4=plt.plot_date(times[start1:end1],btotsvalue[start1:end1],'k',color='black',label=r'$\mathrm{B\,total}$')
	a2.xaxis_date()
	a2.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a2.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a2.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a2.xaxis.set_tick_params(which='minor', labelsize=16)
	plt.minorticks_on()
	plt.text(0.1,0.9,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a2.transAxes)
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
			a2.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	plt.ylabel(r'$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-250,250)
	if savevalue==ron[1]:
		plt.ylim(-150,150)	
	if savevalue==ron[2]:
		plt.ylim(-100,100)
	if savevalue==ron[3]:
		plt.ylim(-70,70)
	if savevalue==ron[4]:
		plt.ylim(-60,60)	
	if savevalue==ron[5]:
		plt.ylim(-50,50)
	if savevalue==ron[6]:
		plt.ylim(-30,30)
	plt.xlim(times[start1],times[end1])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
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
	a3.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
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
		plt.ylim(-290,290)
	if savevalue==ron[1]:
		plt.ylim(-150,150)	
	if savevalue==ron[2]:
		plt.ylim(-100,100)
	if savevalue==ron[3]:
		plt.ylim(-60,60)
	if savevalue==ron[4]:
		plt.ylim(-45,45)	
	if savevalue==ron[5]:
		plt.ylim(-40,40)
	if savevalue==ron[6]:
		plt.ylim(-30,30)
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
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
	a4.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a4.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a4.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a4.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.9,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a4.transAxes)
	plt.ylabel('$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-290,290)
	if savevalue==ron[1]:
		plt.ylim(-150,150)	
	if savevalue==ron[2]:
		plt.ylim(-100,100)
	if savevalue==ron[3]:
		plt.ylim(-60,60)
	if savevalue==ron[4]:
		plt.ylim(-45,45)	
	if savevalue==ron[5]:
		plt.ylim(-40,40)
	if savevalue==ron[6]:
		plt.ylim(-30,30)
	plt.xlim(times[start2],times[end2])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
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
	a5.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a5.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a5.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a5.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.3,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a5.transAxes)
	plt.text(0.9,0.4, '$\mathrm{%s}$'%(stmin),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a5.transAxes)
	plt.ylabel('$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')    
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=18,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-600,70)
	if savevalue==ron[1]:
		plt.ylim(-400,70)	
	if savevalue==ron[2]:
		plt.ylim(-300,50)
	if savevalue==ron[3]:
		plt.ylim(-200,40)
	if savevalue==ron[4]:
		plt.ylim(-150,40)	
	if savevalue==ron[5]:
		plt.ylim(-120,40)
	if savevalue==ron[6]:
		plt.ylim(-100,40)
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)

	a6=f3.add_subplot(212)
	#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
	a6.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a6.transAxes,fontsize=17)
	plt.plot_date(times[start1:end1],dstscaledvalue[start1:end1],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a6.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a6.xaxis_date()
	a6.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a6.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a6.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a6.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.text(0.9,0.3,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a6.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-600,90)
	if savevalue==ron[1]:
		plt.ylim(-400,90)	
	if savevalue==ron[2]:
		plt.ylim(-300,50)
	if savevalue==ron[3]:
		plt.ylim(-200,40)
	if savevalue==ron[4]:
		plt.ylim(-150,40)	
	if savevalue==ron[5]:
		plt.ylim(-120,40)
	if savevalue==ron[6]:
		plt.ylim(-100,40)
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
	a7.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a7.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a7.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a7.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.3,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a7.transAxes)
	plt.text(0.9,0.4, '$\mathrm{%s}$'%(stmax),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a7.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')    
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
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	
	a8=f4.add_subplot(212)
	#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
	a8.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a8.transAxes,fontsize=17)
	plt.plot_date(times[start2:end2],dstscaledvalue[start2:end2],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a8.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a8.xaxis_date()
	a8.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a8.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a8.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a8.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.text(0.9,0.3,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a8.transAxes)
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
	a9.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
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
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a10 = f5.add_subplot(212)
	a10.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a10.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],dstvalue[start2n:end2n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start2n:end2n],dstobrienvalue[start2n:end2n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a10.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a10.xaxis_date()
	a10.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
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
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
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
	a11.text(-0.06, 1.0,'$\mathrm{a)}$',transform=a11.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],bxvalue[start1n:end1n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start1n:end1n],byvalue[start1n:end1n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start1n:end1n],bzvalue[start1n:end1n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start1n:end1n],btotvalue[start1n:end1n],'k',color='black',label=r'$\mathrm{B\,total}$')	
	a11.xaxis_date()
	a11.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
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
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a12 = f6.add_subplot(212)
	a12.text(-0.06, 1.0,'$\mathrm{b)}$',transform=a12.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],dstvalue[start1n:end1n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start1n:end1n],dstobrienvalue[start1n:end1n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a12.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a12.xaxis_date()
	a12.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
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
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
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
	a91.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
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
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a101 = f51.add_subplot(212)
	a101.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a101.transAxes,fontsize=17)
	plt.plot_date(timen[start2n:end2n],dstvalue[start2n:end2n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start2n:end2n],dstobrienvalue[start2n:end2n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a101.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a101.xaxis_date()
	a101.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a101.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a101.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a101.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(1.18,0.16,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a101.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')    
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	plt.xlim(timen[start2n],timen[end2n])
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
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
	a111.text(-0.06, 1.0,'$\mathrm{a)}$',transform=a111.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],bxvalue[start1n:end1n],'k',color='red',label=r'$\mathrm{Bx\,component}$') 
	plt.plot_date(timen[start1n:end1n],byvalue[start1n:end1n],'k',color='green',label=r'$\mathrm{By\,component}$')
	plt.plot_date(timen[start1n:end1n],bzvalue[start1n:end1n],'k',color='blue',label=r'$\mathrm{Bz\,component}$')
	plt.plot_date(timen[start1n:end1n],btotvalue[start1n:end1n],'k',color='black',label=r'$\mathrm{B\,total}$')	
	a111.xaxis_date()
	a111.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
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
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 

	a121 = f16.add_subplot(212)
	a121.text(-0.06, 1.0,'$\mathrm{b)}$',transform=a121.transAxes,fontsize=17)
	plt.plot_date(timen[start1n:end1n],dstvalue[start1n:end1n],'k',color='lightcoral',label='$\mathrm{Observed\,hourly\,Dst}$',linestyle='--')
	plt.plot_date(timen[start1n:end1n],dstobrienvalue[start1n:end1n],'k',color='darkmagenta', label='$\mathrm{OBrien\,McPherron\,2000}$')
	for p in range(len(dict['st_normal'])):
		a121.axvspan(timen[dict['st_normal'][p]],timen[dict['en_normal'][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a121.xaxis_date()
	a121.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a121.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a121.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a121.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(1.18,0.16,'$\mathrm{@1AU}$',fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a121.transAxes)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')    
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	plt.xlim(timen[start1n],timen[end1n])
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
	plt.legend(loc=4,ncol=3, fancybox=True,fontsize=16)
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
	a13.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a13.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a13.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a13.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmaxs),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a13.transAxes)
	plt.text(0.9,0.8,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a13.transAxes)
	plt.ylabel('$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-290,290)
	if savevalue==ron[1]:
		plt.ylim(-150,150)	
	if savevalue==ron[2]:
		plt.ylim(-100,100)
	if savevalue==ron[3]:
		plt.ylim(-60,60)
	if savevalue==ron[4]:
		plt.ylim(-45,45)	
	if savevalue==ron[5]:
		plt.ylim(-40,40)
	if savevalue==ron[6]:
		plt.ylim(-30,30)
	plt.xlim(times[start2],times[end2])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 	
	
	
	a14=f7.add_subplot(212)
	a14.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a14.transAxes,fontsize=17)
	plt.plot_date(times[start2:end2],dstscaledvalue[start2:end2],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a14.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a14.xaxis_date()
	a14.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a14.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a14.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a14.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
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
	a15.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a15.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a15.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a15.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.text(0.9,0.9, '$\mathrm{%s}$'%(stmins),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a15.transAxes)
	plt.text(0.9,0.8,'$\mathrm{@0.%dAU}$'%(savevalue),fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform = a15.transAxes)
	plt.ylabel('$\mathrm{B\,\,[nT]}$',fontsize=20,fontweight='bold')
	#plt.xlabel('$\mathrm{time\,[days]}$',fontsize=14,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-250,250)
	if savevalue==ron[1]:
		plt.ylim(-150,150)	
	if savevalue==ron[2]:
		plt.ylim(-100,100)
	if savevalue==ron[3]:
		plt.ylim(-70,70)
	if savevalue==ron[4]:
		plt.ylim(-60,60)	
	if savevalue==ron[5]:
		plt.ylim(-50,50)
	if savevalue==ron[6]:
		plt.ylim(-30,30)
	plt.xlim(times[start1],times[end1])
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=4, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17) 	
	
	
	a16=f8.add_subplot(212)
	a16.text(-0.07, 1.0,'$\mathrm{b)}$',transform=a16.transAxes,fontsize=17)
	plt.plot_date(times[start1:end1],dstscaledvalue[start1:end1],'k',color='darkmagenta', label='$\mathrm{scaled\,Dst}$') 
	for p in range(len(dict['st_ni0%d'%(savevalue)])):
		a16.axvspan(times[dict['st_ni0%d'%(savevalue)][p]],times[dict['en_ni0%d'%(savevalue)][p]],alpha=0.5,color='grey',label='Magnetic Obstacle'if p == 0 else '')
	a16.xaxis_date()
	a16.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))   	# normal '%H:%M'
	a16.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))   #normal '%Y.%m.%d'
	a16.xaxis.set_tick_params(tickdir='out',length=5,color='black',pad=15, labelsize=17)
	a16.xaxis.set_tick_params(which='minor', labelsize=16)	
	plt.minorticks_on()
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.legend(loc=3, ncol=3, fancybox=True, shadow=True,fontsize=16)
	plt.ylabel(r'$\mathrm{Dst\,\,[nT]}$',fontsize=20,fontweight='bold')
	plt.xlabel('$\mathrm{time\,[days]}$',fontsize=20,fontweight='bold')
	if savevalue==ron[0]:
		plt.ylim(-600,90)
	if savevalue==ron[1]:
		plt.ylim(-400,90)	
	if savevalue==ron[2]:
		plt.ylim(-300,50)
	if savevalue==ron[3]:
		plt.ylim(-200,40)
	if savevalue==ron[4]:
		plt.ylim(-150,40)	
	if savevalue==ron[5]:
		plt.ylim(-120,40)
	if savevalue==ron[6]:
		plt.ylim(-100,40)
	plt.xlim(times[start1],times[end1])
	plt.yticks(fontsize=17) 
	plt.xticks(fontsize=17)
	plt.tight_layout()
	plt.savefig(Scaled_fluxrope_solar+"/scaled_fluxrope_0.%dAU_solarmin.pdf"%(savevalue),dpi=300)
	plt.savefig(Scaled_fluxrope_solar+"/scaled_fluxrope_0.%dAU_solarmin.png"%(savevalue),dpi=300)
	plt.close()
	
	
	return (f1,f2,f3,f4,f5,f6,f7,f8,f16,f51)

#make scaled plots from 0.3 AU until 0.9 AU	
#start1,end1,start2,end2,timen,bxvalue,byvalue,bzvalue,btotvalue,dstvalue,dstobrienvalue,times,bxsvalue,bysvalue,bzsvalue,btotsvalue,dstscaledvalue,savevalue		
plotcomparison(start2008_03,end2008_03,start2013_03,end2014_03,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_03AU'],dict['bx_i_mean_scaled_03AU'],dict['by_i_mean_scaled_03AU'],dict['bz_i_mean_scaled_03AU'],dict['btot_i_mean_scaled_03AU'],dict['dst_i_mean_scaled_03AU'],ron[0])		
plotcomparison(start2008_04,end2008_04,start2013_04,end2014_04,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_04AU'],dict['bx_i_mean_scaled_04AU'],dict['by_i_mean_scaled_04AU'],dict['bz_i_mean_scaled_04AU'],dict['btot_i_mean_scaled_04AU'],dict['dst_i_mean_scaled_04AU'],ron[1])		
plotcomparison(start2008_05,end2008_05,start2013_05,end2014_05,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_05AU'],dict['bx_i_mean_scaled_05AU'],dict['by_i_mean_scaled_05AU'],dict['bz_i_mean_scaled_05AU'],dict['btot_i_mean_scaled_05AU'],dict['dst_i_mean_scaled_05AU'],ron[2])		
plotcomparison(start2008_06,end2008_06,start2013_06,end2014_06,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_06AU'],dict['bx_i_mean_scaled_06AU'],dict['by_i_mean_scaled_06AU'],dict['bz_i_mean_scaled_06AU'],dict['btot_i_mean_scaled_06AU'],dict['dst_i_mean_scaled_06AU'],ron[3])		
plotcomparison(start2008_07,end2008_07,start2013_07,end2014_07,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_07AU'],dict['bx_i_mean_scaled_07AU'],dict['by_i_mean_scaled_07AU'],dict['bz_i_mean_scaled_07AU'],dict['btot_i_mean_scaled_07AU'],dict['dst_i_mean_scaled_07AU'],ron[4])		
plotcomparison(start2008_08,end2008_08,start2013_08,end2014_08,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_08AU'],dict['bx_i_mean_scaled_08AU'],dict['by_i_mean_scaled_08AU'],dict['bz_i_mean_scaled_08AU'],dict['btot_i_mean_scaled_08AU'],dict['dst_i_mean_scaled_08AU'],ron[5])		
plotcomparison(start2008_09,end2008_09,start2013_09,end2014_09,start2008_1,end2008_1,start2013_1,end2014_1,dict['time_normal'],dict['bx_normal'],dict['by_normal'],dict['bz_normal'],dict['btot_normal'],dict['dst_normal'],dict['dst_obrien_normal'],dict['time_i_scaled_09AU'],dict['bx_i_mean_scaled_09AU'],dict['by_i_mean_scaled_09AU'],dict['bz_i_mean_scaled_09AU'],dict['btot_i_mean_scaled_09AU'],dict['dst_i_mean_scaled_09AU'],ron[6])		

		
		
#make histogram fits
#formula from Wanliss et al.(2005), doi: 10.1029/2004JA010996
def pdi(x,a1,b1,c1,a2,b2,c2):
	return (a1*np.exp(-((np.log(x)-b1)**2)/(c1**2))+a2*np.exp(-((np.log(x)-b2)**2)/(c2**2)))

#make histograms 
#### 1 AU solar minimum 
yhist1_2008,xhist1_2008=np.histogram(dict['dst_normal'][start2008_1:end2008_1],bins=np.arange(-100,50,5))
yhist1_2008=yhist1_2008/len(dict['dst_normal'][start2008_1:end2008_1])
xhist1new_2008=-xhist1_2008+100
init_vals_2008 = [1e-2, 1, 1,1e-2, 1, 1]    # for [a1, b1, c1]  [1e-3,1,1,1e-3,1,1]
best_vals1_2008, covar1_2008 = curve_fit(pdi, xhist1new_2008[1:],yhist1_2008, p0=init_vals_2008) 
xwerte1_2008=np.arange(-100,50,0.1)
xwertenew1_2008=-xwerte1_2008+100
fit1_2008=pdi(xwertenew1_2008,best_vals1_2008[0],best_vals1_2008[1],best_vals1_2008[2],best_vals1_2008[3],best_vals1_2008[4],best_vals1_2008[5])
f1=plt.figure(31,figsize=(16,10))
plt.bar(xhist1_2008[1:],yhist1_2008,color='orangered',edgecolor='orangered',linewidth=3,fill=True,label='histogram for 1AU - solar min')
plt.plot(xwerte1_2008, fit1_2008, color='red',linewidth=3, label='fit for 1AU - solar min')
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
#plt.text(,'$\mathrm{period: }$',fontsize=14,color='black',horizontalalignment='center',verticalalignment='center',transform = .transAxes)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histo_solar+"/1AU_solarmin.pdf",dpi=300)
plt.savefig(histo_solar+"/1AU_solarmin.png",dpi=300)	

#### 1 AU solar maximum 
yhist1_2013,xhist1_2013=np.histogram(dict['dst_normal'][start2013_1:end2014_1],bins=np.arange(-100,50,5))
yhist1_2013=yhist1_2013/len(dict['dst_normal'][start2013_1:end2014_1])
xhist1new_2013=-xhist1_2013+100
init_vals_2013 = [1e-2, 1, 1,1e-2, 1, 1]    # for [a1, b1, c1]  [1e-3,1,1,1e-3,1,1]
best_vals1_2013, covar1_2013 = curve_fit(pdi, xhist1new_2013[1:],yhist1_2013, p0=init_vals_2013) 
xwerte1_2013=np.arange(-100,50,0.1)
xwertenew1_2013=-xwerte1_2013+100
fit1_2013=pdi(xwertenew1_2013,best_vals1_2013[0],best_vals1_2013[1],best_vals1_2013[2],best_vals1_2013[3],best_vals1_2013[4],best_vals1_2013[5])
fio=plt.figure(32,figsize=(16,10))
plt.bar(xhist1_2013[1:],yhist1_2013,color='dodgerblue',edgecolor='dodgerblue',linewidth=3,fill=True,label='histogram for 1AU - solar max ')
plt.plot(xwerte1_2013, fit1_2013, color='blue',label='fit for 1AU - solar max')
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)
plt.savefig(histo_solar+"/1AU_solarmax.pdf",dpi=300)
plt.savefig(histo_solar+"/1AU_solarmax.png",dpi=300)	

############ 0.3 AU solar minimum 
yhist3_2008,xhist3_2008=np.histogram(dict['dst_i_mean_scaled_03AU'][start2008_03:end2008_03],bins=np.arange(-500,0,5))
yhist3_2008=yhist3_2008/len(dict['dst_i_mean_scaled_03AU'][start2008_03:end2008_03])
xhist3n_2008=-xhist3_2008
init_vals3_2008 = [1e-2, 1, 1,1e-2,1, 1]
best_vals3_2008, covar3_2008 = curve_fit(pdi, xhist3n_2008[1:],yhist3_2008, p0=init_vals3_2008) 
xwerte3_2008=np.arange(-500,0,0.1)
xwertenew3_2008=-xwerte3_2008
fit3_2008=pdi(xwertenew3_2008,best_vals3_2008[0],best_vals3_2008[1],best_vals3_2008[2],best_vals3_2008[3],best_vals3_2008[4],best_vals3_2008[5])
figur=plt.figure(22,figsize=(16,10))
plt.bar(xhist3_2008[1:],yhist3_2008,color='green',edgecolor='green',linewidth=3,fill=True,label='histogram for 0.3AU - solar min')
plt.plot(xwerte3_2008, fit3_2008, color='saddlebrown',label='fit for 0.3AU - solar min')
plt.xlim(-500,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histo_solar+"/03AU_solarmin.pdf",dpi=300)
plt.savefig(histo_solar+"/03AU_solarmin.png",dpi=300)

####### flipped histogram ############
yhist3r_2008,xhist3r_2008=np.histogram(dict['dst_i_mean_scaled_03AU'][start2008_03:end2008_03],bins=np.arange(-500,0,5))
yhist3r_2008=yhist3r_2008/len(dict['dst_i_mean_scaled_03AU'][start2008_03:end2008_03])
xhist3nr_2008=-xhist3r_2008
init_vals3r_2008 = [1e-2, 1, 1,1e-2,1, 1]
best_vals3r_2008, covar3r_2008 = curve_fit(pdi, xhist3nr_2008[1:],yhist3r_2008, p0=init_vals3r_2008) 
xwerte3r_2008=np.arange(-500,0,0.1)
xwertenew3r_2008=-xwerte3r_2008
fit3r_2008=pdi(xwertenew3r_2008,best_vals3r_2008[0],best_vals3r_2008[1],best_vals3r_2008[2],best_vals3r_2008[3],best_vals3r_2008[4],best_vals3r_2008[5])
figurk=plt.figure(3567,figsize=(16,10))
plt.bar(-xhist3r_2008[1:],yhist3r_2008,color='orange',edgecolor='orange',linewidth=3,fill=True,label='@0.3AU - histogram with positive values')
plt.plot(xwertenew3r_2008, fit3r_2008, color='red',linewidth=3,label='@0.3AU - fit for positive histogram')
plt.bar(xhist3r_2008[1:],yhist3r_2008,color='green',edgecolor='green',linewidth=3,fill=True,label='@0.3AU - histogram reversed to negative axis')
plt.plot(xwerte3r_2008, fit3r_2008, color='saddlebrown',linewidth=3,label='@0.3AU - fit reversed to negative axis')
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.legend(loc=2,fontsize=15)
plt.legend(loc=2,fontsize=15)
plt.legend(loc=2,fontsize=15)
plt.legend(loc=2,fontsize=15)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histo_solar+"/flippedhisto_solarmin.pdf",dpi=300)
plt.savefig(histo_solar+"/flippedhisto_solarmin.png",dpi=300)


#### 0.9 AU solar minimum 
yhist9_2008,xhist9_2008=np.histogram(dict['dst_i_mean_scaled_09AU'][start2008_09:end2008_09],bins=np.arange(-100,50,5))
yhist9_2008=yhist9_2008/len(dict['dst_i_mean_scaled_09AU'][start2008_09:end2008_09])
xhist9n_2008=-xhist9_2008+100
init_vals9_2008 = [1e-2, 1, 1,1e-2,1, 1]
best_vals9_2008, covar9_2008 = curve_fit(pdi, xhist9n_2008[1:],yhist9_2008, p0=init_vals9_2008) 
xwerte9_2008=np.arange(-100,50,0.1)
xwertenew9_2008=-xwerte9_2008+100
fit9_2008=pdi(xwertenew9_2008,best_vals9_2008[0],best_vals9_2008[1],best_vals9_2008[2],best_vals9_2008[3],best_vals9_2008[4],best_vals9_2008[5])
figur9=plt.figure(292,figsize=(16,10))
plt.bar(xhist9_2008[1:],yhist9_2008,color='green',edgecolor='green',linewidth=3,fill=True,label='histogram for 0.9AU - solar min')
plt.plot(xwerte9_2008, fit9_2008, color='saddlebrown',linewidth=3,label='fit for 0.9AU - solar min')
plt.xlim(-100,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histo_solar+"/09AU_solarmin.pdf",dpi=300)
plt.savefig(histo_solar+"/09AU_solarmin.png",dpi=300)

#### 0.9 AU solar maximum 
yhist9_2013,xhist9_2013=np.histogram(dict['dst_i_mean_scaled_09AU'][start2013_09:end2014_09],bins=np.arange(-100,50,5))
yhist9_2013=yhist9_2013/len(dict['dst_i_mean_scaled_09AU'][start2013_09:end2014_09])
xhist9n_2013=-xhist9_2013+100
init_vals9_2013 = [1e-2, 1, 1,1e-2,1, 1]
best_vals9_2013, covar9_2013 = curve_fit(pdi, xhist9n_2013[1:],yhist9_2013, p0=init_vals9_2013) 
xwerte9_2013=np.arange(-100,50,0.1)
xwertenew9_2013=-xwerte9_2013+100
fit9_2013=pdi(xwertenew9_2013,best_vals9_2013[0],best_vals9_2013[1],best_vals9_2013[2],best_vals9_2013[3],best_vals9_2013[4],best_vals9_2013[5])
fegu9=plt.figure(239,figsize=(16,10))
plt.bar(xhist9_2013[1:],yhist9_2013,color='darkmagenta',edgecolor='darkmagenta',linewidth=3,fill=True,label='histogram for 0.9AU - solar max')
plt.plot(xwerte9_2013, fit9_2013, color='crimson',linewidth=3,label='fit for 0.9AU - solar max')
plt.xlim(-100,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histo_solar+"/09AU_solarmax.pdf",dpi=300)
plt.savefig(histo_solar+"/09AU_solarmax.png",dpi=300)


#create histogram plot function 
def histo(dstAU,st1,en1,st2,en2,savevalue1):
	fa=list()
	#############
	yhist3_2013,xhist3_2013=np.histogram(dstAU[st2:en2],bins=np.arange(-500,0,5))
	yhist3_2013=yhist3_2013/len(dstAU[st2:en2])
	xhist3n_2013=-xhist3_2013
	init_vals3_2013 = [1e-2, 1, 1,1e-2,1, 1]
	best_vals3_2013, covar3_2013 = curve_fit(pdi, xhist3n_2013[1:],yhist3_2013, p0=init_vals3_2013) 
	xwerte3_2013=np.arange(-500,0,0.1)
	xwertenew3_2013=-xwerte3_2013
	fit3_2013=pdi(xwertenew3_2013,best_vals3_2013[0],best_vals3_2013[1],best_vals3_2013[2],best_vals3_2013[3],best_vals3_2013[4],best_vals3_2013[5])
	fo=plt.figure(35,figsize=(16,10))
	plt.bar(xhist3_2013[1:],yhist3_2013,color='darkmagenta',edgecolor='darkmagenta',linewidth=3,fill=True,label='histogram for 0.%dAU - solar maximum'%(savevalue1))
	plt.plot(xwerte3_2013, fit3_2013, color='crimson',linewidth=3,label='fit for 0.%dAU - solar max'%(savevalue1))
	plt.xlim(-500,50)
	plt.ylabel('Probability density',fontsize=28,labelpad=15)    
	plt.xlabel('Dst values [nT]',fontsize=28)
	plt.legend(loc=2,fontsize=21)
	plt.legend(loc=2,fontsize=21)
	plt.yticks(fontsize=22) 
	plt.xticks(fontsize=22)
	plt.savefig(fits_histo+"/0%dAU_auto_solarmax.pdf"%(savevalue1),dpi=300)
	plt.savefig(fits_histo+"/0%dAU_auto_solarmax.png"%(savevalue1),dpi=300)
	plt.close()
	
	############
	yhist3_2008,xhist3_2008=np.histogram(dstAU[st1:en1],bins=np.arange(-500,0,5))
	yhist3_2008=yhist3_2008/len(dstAU[st1:en1])
	xhist3n_2008=-xhist3_2008
	init_vals3_2008 = [1e-2, 1, 1,1e-2,1, 1]
	best_vals3_2008, covar3_2008 = curve_fit(pdi, xhist3n_2008[1:],yhist3_2008, p0=init_vals3_2008) 
	xwerte3_2008=np.arange(-500,0,0.1)
	xwertenew3_2008=-xwerte3_2008
	fit3_2008=pdi(xwertenew3_2008,best_vals3_2008[0],best_vals3_2008[1],best_vals3_2008[2],best_vals3_2008[3],best_vals3_2008[4],best_vals3_2008[5])
	fop=plt.figure(36,figsize=(16,10))
	plt.bar(xhist3_2008[1:],yhist3_2008,color='green',edgecolor='green',linewidth=3,fill=True,label='histogram for 0.%dAU - solar min'%(savevalue1))
	plt.plot(xwerte3_2008, fit3_2008, color='saddlebrown',linewidth=3,label='fit for 0.%dAU - solar min'%(savevalue1))
	plt.xlim(-500,50)
	plt.ylabel('Probability density',fontsize=28,labelpad=15)    
	plt.xlabel('Dst values [nT]',fontsize=28)
	plt.legend(loc=2,fontsize=21)
	plt.legend(loc=2,fontsize=21)
	plt.yticks(fontsize=22) 
	plt.xticks(fontsize=22)
	plt.savefig(fits_histo+"/0%dAU_auto_solarmin.pdf"%(savevalue1),dpi=300)
	plt.savefig(fits_histo+"/0%dAU_auto_solarmin.png"%(savevalue1),dpi=300)
	plt.close()
	fa.append(best_vals3_2008[0])
	fa.append(best_vals3_2008[1])
	fa.append(best_vals3_2008[2])
	fa.append(best_vals3_2008[3])
	fa.append(best_vals3_2008[4])
	fa.append(best_vals3_2008[5])
	fa.append(best_vals3_2013[0])
	fa.append(best_vals3_2013[1])
	fa.append(best_vals3_2013[2])
	fa.append(best_vals3_2013[3])
	fa.append(best_vals3_2013[4])
	fa.append(best_vals3_2013[5])
	return (fa,fo,fop)
	
d3=histo(dict['dst_i_mean_scaled_03AU'],start2008_03,end2008_03,start2013_03,end2014_03,ron[0])
d4=histo(dict['dst_i_mean_scaled_04AU'],start2008_04,end2008_04,start2013_04,end2014_04,ron[1])
d5=histo(dict['dst_i_mean_scaled_05AU'],start2008_05,end2008_05,start2013_05,end2014_05,ron[2])
d6=histo(dict['dst_i_mean_scaled_06AU'],start2008_06,end2008_06,start2013_06,end2014_06,ron[3])
d7=histo(dict['dst_i_mean_scaled_07AU'],start2008_07,end2008_07,start2013_07,end2014_07,ron[4])
d8=histo(dict['dst_i_mean_scaled_08AU'],start2008_08,end2008_08,start2013_08,end2014_08,ron[5])



#### 0.3 AU solar maximum 
yhist3_2013,xhist3_2013=np.histogram(dict['dst_i_mean_scaled_03AU'][start2013_03:end2014_03],bins=np.arange(-500,0,5))
yhist3_2013=yhist3_2013/len(dict['dst_i_mean_scaled_03AU'][start2013_03:end2014_03])
xhist3n_2013=-xhist3_2013
init_vals3_2013 = [1e-2, 1, 1,1e-2,1, 1]
best_vals3_2013, covar3_2013 = curve_fit(pdi, xhist3n_2013[1:],yhist3_2013, p0=init_vals3_2013) 
xwerte3_2013=np.arange(-500,0,0.1)
xwertenew3_2013=-xwerte3_2013
fit3_2013=pdi(xwertenew3_2013,best_vals3_2013[0],best_vals3_2013[1],best_vals3_2013[2],best_vals3_2013[3],best_vals3_2013[4],best_vals3_2013[5])
fegu=plt.figure(23,figsize=(16,10))
plt.bar(xhist3_2013[1:],yhist3_2013,color='darkmagenta',edgecolor='darkmagenta',linewidth=3,fill=True,label='histogram for 0.3AU - solar max')
plt.plot(xwerte3_2013, fit3_2013, color='crimson',label='fit for 0.3AU - solar max')
plt.xlim(-500,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histo_solar+"/03AU_solarmax.pdf",dpi=300)
plt.savefig(histo_solar+"/03AU_solarmax.png",dpi=300)

#### 0.3 AU solar minimum 
fiii=plt.figure(24,figsize=(16,10))
plt.bar(xhist3_2008[1:],yhist3_2008,color='green',edgecolor='green',linewidth=3,fill=True,label='histogram for 0.3AU - solar min')
plt.plot(xwerte3_2008, fit3_2008, color='saddlebrown',linewidth=3, label='fit for 0.3AU - solar min')
plt.bar(xhist3_2013[1:],yhist3_2013,color='darkmagenta',edgecolor='darkmagenta',linewidth=3,fill=True,label='histogram for 0.3AU - solar max')
plt.plot(xwerte3_2013, fit3_2013, color='crimson',linewidth=3, label='fit for 0.3AU - solar max')
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.xlim(-500,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histof_solar+"/03histograms_solarcycle.pdf",dpi=300)
plt.savefig(histof_solar+"/03histograms_solarcycle.png",dpi=300)

#### 0.3 AU comparison of solar minimum & maximum 
fixi=plt.figure(26,figsize=(16,10))
plt.bar(xhist3_2008[1:],yhist3_2008,color='green',edgecolor='green',linewidth=3,fill=True,label='histogram for 0.3AU - solar min')
plt.plot(xwerte3_2008, fit3_2008, color='saddlebrown',linewidth=3,label='fit for 0.3AU - solar min')
plt.bar(xhist3_2013[1:],yhist3_2013,color='darkmagenta',edgecolor='darkmagenta',linewidth=3,fill=True,label='histogram for 0.3AU - solar max')
plt.plot(xwerte3_2013, fit3_2013, color='crimson',linewidth=3,label='fit for 0.3AU - solar max')
plt.bar(xhist1_2008[1:],yhist1_2008,color='orangered',edgecolor='orangered',linewidth=3,fill=True,label='histogram for 1AU - solar min')
plt.plot(xwerte1_2008, fit1_2008, color='red',linewidth=3,label='fit for 1AU - solar min')
plt.bar(xhist1_2013[1:],yhist1_2013,color='dodgerblue',edgecolor='dodgerblue',linewidth=3,fill=True,label='histogram for 1AU - solar max')
plt.plot(xwerte1_2013, fit1_2013, color='blue',linewidth=3,label='fit for 1AU - solar max')
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.xlim(-500,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histof_solar+"/allhistograms_solarcycle1AU.pdf",dpi=300)
plt.savefig(histof_solar+"/allhistograms_solarcycle1AU.png",dpi=300)

#### 1 AU comparison of solar minimum & maximum 
fixin=plt.figure(29,figsize=(16,10))
plt.bar(xhist1_2008[1:],yhist1_2008,color='orangered',edgecolor='orangered',linewidth=3,fill=True,label='histogram for 1AU - solar min')
plt.plot(xwerte1_2008, fit1_2008, color='red',linewidth=3, label='fit for 1AU - solar min')
plt.bar(xhist1_2013[1:],yhist1_2013,color='dodgerblue',edgecolor='dodgerblue',linewidth=3,fill=True,label='histogram for 1AU - solar max')
plt.plot(xwerte1_2013, fit1_2013, color='blue',linewidth=3, label='fit for 1AU - solar max')
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.xlim(-110,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(histof_solar+"/1AUhistograms_solarcycle1AU.pdf",dpi=300)
plt.savefig(histof_solar+"/1AUhistograms_solarcycle1AU.png",dpi=300)

#### 0.3 AU fit comparison for solar minimum & maximum 
fee=plt.figure(25,figsize=(16,10))
plt.plot(xwerte3_2008, fit3_2008, color='saddlebrown',linewidth=3, label='fit for 0.3AU - solar min')
plt.plot(xwerte3_2013, fit3_2013, color='crimson',linewidth=3, label='fit for 0.3AU - solar max')
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.xlim(-500,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(fits_solar+"/03fits_solarcycle.pdf",dpi=300)
plt.savefig(fits_solar+"/03fits_solarcycle.png",dpi=300)

#### 1 AU fit comparison for solar minimum & maximum 
fexei=plt.figure(28,figsize=(16,10))
plt.plot(xwerte1_2008, fit1_2008, color='red',linewidth=3,label='fit for 1AU - solar min')
plt.plot(xwerte1_2013, fit1_2013, color='blue',linewidth=3,label='fit for 1AU - solar max')
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.xlim(-110,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(fits_solar+"/1AUfits_solarcycle1AU.pdf",dpi=300)
plt.savefig(fits_solar+"/1AUfits_solarcycle1AU.png",dpi=300)

#### fit comparison for solar minimum & maximum for 1 AU & 0.3 AU
fexe=plt.figure(27,figsize=(16,10))
plt.plot(xwerte3_2008, fit3_2008, color='saddlebrown',linewidth=3,label='fit for 0.3AU - solar min')
plt.plot(xwerte3_2013, fit3_2013, color='green',linewidth=3,label='fit for 0.3AU - solar max')
plt.plot(xwerte1_2008, fit1_2008, color='red',linewidth=3,label='fit for 1AU - solar min')
plt.plot(xwerte1_2013, fit1_2013, color='blue',linewidth=3,label='fit for 1AU - solar max')
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.legend(loc=2,fontsize=21)
plt.xlim(-500,50)
plt.ylabel('Probability density',fontsize=28,labelpad=15)    
plt.xlabel('Dst values [nT]',fontsize=28)
plt.yticks(fontsize=22) 
plt.xticks(fontsize=22)
plt.savefig(fits_solar+"/allfits_solarcycle1AU.pdf",dpi=300)
plt.savefig(fits_solar+"/allfits_solarcycle1AU.png",dpi=300)

#create table of fit values 
table1=pd.DataFrame([["%.4f"%best_vals3_2008[0],"%.4f"%best_vals3_2008[1],"%.4f"%best_vals3_2008[2],"%.4f"%best_vals3_2008[3],"%.4f"%best_vals3_2008[4],"%.4f"%best_vals3_2008[5]],
	["%.4f"%best_vals3_2013[0],"%.4f"%best_vals3_2013[1],"%.4f"%best_vals3_2013[2],"%.4f"%best_vals3_2013[3],"%.4f"%best_vals3_2013[4],"%.4f"%best_vals3_2013[5]],
	["%.4f"%d4[0][0],"%.4f"%d4[0][1],"%.4f"%d4[0][2],"%.4f"%d4[0][3],"%.4f"%d4[0][4],"%.4f"%d4[0][5]],["%.4f"%d4[0][6],"%.4f"%d4[0][7],"%.4f"%d4[0][8],"%.4f"%d4[0][9],"%.4f"%d4[0][10],"%.4f"%d4[0][11]],
	["%.4f"%d5[0][0],"%.4f"%d5[0][1],"%.4f"%d5[0][2],"%.4f"%d5[0][3],"%.4f"%d5[0][4],"%.4f"%d5[0][5]],["%.4f"%d5[0][6],"%.4f"%d5[0][7],"%.4f"%d5[0][8],"%.4f"%d5[0][9],"%.4f"%d5[0][10],"%.4f"%d5[0][11]],
	["%.4f"%d6[0][0],"%.4f"%d6[0][1],"%.4f"%d6[0][2],"%.4f"%d6[0][3],"%.4f"%d6[0][4],"%.4f"%d6[0][5]],["%.4f"%d6[0][6],"%.4f"%d6[0][7],"%.4f"%d6[0][8],"%.4f"%d6[0][9],"%.4f"%d6[0][10],"%.4f"%d6[0][11]],
	["%.4f"%d7[0][0],"%.4f"%d7[0][1],"%.4f"%d7[0][2],"%.4f"%d7[0][3],"%.4f"%d7[0][4],"%.4f"%d7[0][5]],["%.4f"%d7[0][6],"%.4f"%d7[0][7],"%.4f"%d7[0][8],"%.4f"%d7[0][9],"%.4f"%d7[0][10],"%.4f"%d7[0][11]],
	["%.4f"%d8[0][0],"%.4f"%d8[0][1],"%.4f"%d8[0][2],"%.4f"%d8[0][3],"%.4f"%d8[0][4],"%.4f"%d8[0][5]],["%.4f"%d8[0][6],"%.4f"%d8[0][7],"%.4f"%d8[0][8],"%.4f"%d8[0][9],"%.4f"%d8[0][10],"%.4f"%d8[0][11]],
	["%.4f"%best_vals9_2008[0],"%.4f"%best_vals9_2008[1],"%.4f"%best_vals9_2008[2],"%.4f"%best_vals9_2008[3],"%.4f"%best_vals9_2008[4],"%.4f"%best_vals9_2008[5]],
	["%.4f"%best_vals9_2013[0],"%.4f"%best_vals9_2013[1],"%.4f"%best_vals9_2013[2],"%.4f"%best_vals9_2013[3],"%.4f"%best_vals9_2013[4],"%.4f"%best_vals9_2013[5]],
	["%.4f"%best_vals1_2008[0],"%.4f"%best_vals1_2008[1],"%.4f"%best_vals1_2008[2],"%.4f"%best_vals1_2008[3],"%.4f"%best_vals1_2008[4],"%.4f"%best_vals1_2008[5]],
	["%.4f"%best_vals1_2013[0],"%.4f"%best_vals1_2013[1],"%.4f"%best_vals1_2013[2],"%.4f"%best_vals1_2013[3],"%.4f"%best_vals1_2013[4],"%.4f"%best_vals1_2013[5]]],
	index=['0.3 AU solar minimum','0.3 AU solar maximum','0.4 AU solar minimum','0.4 AU solar maximum','0.5 AU solar minimum','0.5 AU solar maximum',
	'0.6 AU solar minimum','0.6 AU solar maximum','0.7 AU solar minimum','0.7 AU solar maximum','0.8 AU solar minimum','0.8 AU solar maximum',
	'0.9 AU solar minimum','0.9 AU solar maximum','1AU solar minimum','1AU solar maximum'],
	columns=['A1','b1','c1','A2','b2','c2'])

with open('table_solar_new.tex','w') as tf:
	tf.write(table1.to_latex())



sys.exit()