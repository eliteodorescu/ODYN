#!/usr/bin/python
import matplotlib as mpl
import matplotlib.dates as mpd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogFormatterMathtext
import matplotlib.ticker as plticker


#######  Axes FORMATTING 
locbt = mpd.HourLocator(interval=1) #interval=1)  # every 20 minutes
fmtbt = mpd.DateFormatter('%H:%M:%S')
locb = mpd.MinuteLocator(interval=45) #interval=1)  # every 20 minutes
fmtb = mpd.DateFormatter('%H:%M')
strf = FormatStrFormatter('%.2f')
dpr = dict(c='#ee4000', marker='o',ls=':', markeredgecolor='#ee4000', markerfacecolor='none')
#legprops = dict(loc=2, numpoints=1, frameon=False, borderaxespad=0.)
majgridprops=dict(which='major', color=[0.3, 0.3, 0.3],ls=':')
mingridprops=dict(which='minor', color=[0.5, 0.5, 0.5],ls=':')
pdfpr=dict(marker='d', ls=':', lw=1.5, markeredgecolor='none')
mpr=dict(ms=10, mfc='None', mew=2.)
mpr2=dict(ms=8, mec='None')
lblpr={'size':'12'}
mpl.rc('font',family='serif')#, weight='bold')
plt.rc('legend',**{'fontsize':10})#, weight='bold')
legprops = dict(bbox_to_anchor=(0.999, 0.5, 0.01, 0.05), loc=6, frameon=False, borderaxespad=0.)
tickprops=dict(pad=8, labelsize=12)
txtpr={'size':'11'}
mym=['o',  'd', 'D','v', '^', 'x', 'd', 'p', 's', 'o',  'd', 'D','v', '^', 'x', 'd', 'p', 's','o',  'd', 'D','v', '^', 'x', 'd', 'p', 's']
mym=1000*['o','d','D','v','^','x','d','p','s','+','*']
myc=1000*['b','g','r','m','c','k','y', '#8A2BE2', '#7FFF00', '#8FBC8F', '#FF1493', '#FFD700', '#FF69B4', '#7CFC00', '#ADD8E6'] 
