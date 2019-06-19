import sys, os, urllib
from spacepy import pycdf
import matplotlib as mpl
import numpy as np
import pylab
from scipy import signal, fftpack, optimize
from scipy.stats import mode
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.dates as mpd
from datetime import datetime, date, time, timedelta
import bisect
import random
from sklearn.metrics import mutual_info_score
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogFormatterMathtext
import math, copy
import re


from dateutil import parser
sys.path.append('../AnalysisMethods/')
import draw_attr as D
os.getcwd()



################################################################
######      READ DATA ALGORITHMS        ########################
################################################################

def CLUSTER_DATA_CDF(default_data_path,input_file_name,satellite,probe):
    input_file = [default_data_path+satellite+'/'+input_file_name]
    #input_file = [input_file_name]    
    if not input_file:
        print ('DATA file path or name incorrect')
    else:
        ifile=pycdf.CDF(str(input_file)[2:-2])
    
    ######## DATA retrieving: vars and attrs
    m_Time = ifile['time_tags__'+probe+'_CP_FGM_FULL'][...]
    mag_data = ifile['B_vec_xyz_gse__'+probe+'_CP_FGM_FULL'][...]
    mag_btot  = ifile['B_mag__'+probe+'_CP_FGM_FULL'][...]

    ######## DATA in a fixed format to be used for further analysis
    m_Data=np.c_[mag_data, mag_btot]
    m_Resolution=float(ifile.attrs['TIME_RESOLUTION'][0]) ## data resolution in seconds
    
    return m_Data,m_Time,m_Resolution
    
def CLUSTER_DATA_TXT(default_data_path,input_file_name,satellite,probe,year):
    input_file = [default_data_path+satellite+'/'+input_file_name]
    
    with open(input_file[0]) as fb: 
        tb,bx,by,bz,bxmfa,bymfa,bzmfa,bperpmfa=[],[],[],[],[],[],[],[]
        for lb in fb.readlines():
            if lb.startswith(year):
                #print ('lb: ', lb)
                lb.strip()    
                lb=re.split(r'[;,\r\t\s]\s*',lb)
                tb.append(str(lb[0])+'-'+str(lb[1])+'-'+str(lb[2])+'T'+str(lb[3])+':'+str(lb[4])+':'+str(lb[5]))
                #print ('tb: ', tb)
                bx.append(float(lb[6]))
                by.append(float(lb[7]))
                bz.append(float(lb[8]))
                #mag_tot.append(float(lb[9]))
                bxmfa.append(float(lb[10]))
                bymfa.append(float(lb[11]))
                bzmfa.append(float(lb[12]))
                bperpmfa.append(np.sqrt(float(lb[10])**2+float(lb[11])**2)) 
        fb.close()

        m_Time=[]
        for ti in range(len(tb)):
            m_Time.append(parser.parse(tb[ti]))     

        m_Data=np.c_[bx, by, bz, bxmfa, bymfa, bzmfa, bperpmfa]
            
    return m_Data,m_Time            
    
def VEX_DATA_TXT(default_data_path,input_file_name,satellite,year):
    input_file = [default_data_path+satellite+'/'+input_file_name]
    
    with open(input_file[0]) as fb:
        ll=fb.readlines()
        lf=(ll[-1]).strip()
        lf=re.split(r'[;,\r\t\s]\s*',lf)
        time_list,var_list=[],[[] for il in range(len(lf)-1)]
        for lb in ll:
                if lb.startswith(year):   
                    lb.strip()    
                    lb=re.split(r'[;,\r\t\s]\s*',lb)
                    time_list.append(str(lb[0]))
                    for ili in range(len(lf)-1):
                        var_list[ili].append(float(lb[ili+1]))
    fb.close()                    
    m_Time=[]
    for ti in range(len(time_list)):
        m_Time.append(parser.parse(time_list[ti]))     
    m_Data=np.c_[tuple(var_list)]
    return m_Data,m_Time

def ULYSSES_DATA_ASC(default_data_path,input_file_name,satellite):
    input_file = [default_data_path+satellite+'/'+input_file_name]
    
    with open(input_file[0]) as fb:
        year, day, hour, minute, sec=[],[],[],[],[]
        ll=fb.readlines()
        lf=(ll[-1]).strip()
        lf=re.split(r'[;,\r\t\s]\s*',lf)
        var_list=[[] for il in range(len(lf)-5)]
        for lb in ll:
            if not lb.startswith('#'):
                lb.strip() 
                lb=re.split(r'[;,\r\t\s]\s*',lb)
                if float(lb[5])>=60.00:
                    lb[4]=int(lb[4])+1
                    lb[5]=00.00
                year.append(str(lb[1]))
                day.append(str(lb[2]))
                hour.append(str(lb[3]))
                minute.append(str(lb[4]))
                sec.append(str(lb[5]))
                for ili in range(len(lf)-5):
                    var_list[ili].append(float(lb[ili+6]))
    fb.close()
    mtime=[]
    for it in range(len(year)):
        if int(year[it])>=90:
            mtime.append('19'+str(year[it])+'-'+str(day[it])+'T'+str(hour[it])+':'+str(minute[it])+':'+str(sec[it]))
        elif int(year[it])<10:
            mtime.append('20'+str(year[it])+'-'+str(day[it])+'T'+str(hour[it])+':'+str(minute[it])+':'+str(sec[it]))


    m_Time=[]
    for ti in range(len(mtime)):
        m_Time.append(datetime.strptime(mtime[ti],"%Y-%jT%H:%M:%S.%f"))     

    m_Data=np.c_[tuple(var_list)]
    return m_Data,m_Time


def USER_DATA_TXT(default_data_path,input_file_name,satellite):
    input_file = [default_data_path+satellite+'/'+input_file_name]
    
    with open(input_file[0]) as fb:
        ll=fb.readlines()
        lf=(ll[-1]).strip()
        lf=re.split(r'[;,\r\t\s]\s*',lf)
        time_list,var_list=[],[[] for il in range(len(lf)-1)]
        print (len(lf))
        for lb in ll:
                if not lb.startswith('#'):
                    lb.strip()    
                    lb=re.split(r'[;,\r\t\s]\s*',lb)
                    time_list.append(float(lb[0]))
                    for ili in range(len(lf)-1):
                        var_list[ili].append(float(lb[ili+1]))
    fb.close()                    
    m_Time=[]
    for ti in range(len(time_list)):
        m_Time.append(datetime.now()+timedelta(seconds=time_list[ti]) )    
    m_Data=np.c_[tuple(var_list)]
    return m_Data,m_Time
###############################################################


################################################################
########    PLOT DATA ALGORITHMS      ##########################
################################################################

def PLOT_DATA(DATA, time, spacecraft, labels=None):
    if labels==None:
        labels=[]
        if len(DATA.shape)>1:
            for il in range(DATA.shape[1]):
                labels.append('var'+str(il+1))
        elif len(DATA.shape)==1:        
            labels.append('var0')
                
                
    fig = plt.figure('DATA PLOT '+str(random.randint(1,100))+str(random.randint(1,100)),figsize = (10,3))#dpi=200,
    m1=fig.add_axes([0.,0.,1.,1.])
    m1.axes.get_xaxis().set_visible(False)
    m1.axes.get_yaxis().set_visible(False)
    mp=fig.add_axes([0.12, 0.12, 0.7, 0.78])
    if len(DATA.shape)>1:
        for i in range(DATA.shape[1]):
             mp.plot(time, DATA[:,i], '+-', label=labels[i])
    elif len(DATA.shape)==1:
        mp.plot(time, DATA, '+-', label=labels)
    else:
        print ('shape of DATA array is not appropriate for this function')
    
    l1=mp.legend(**D.legprops)
    f1=l1.get_frame()
    f1.set_color('None')
    mp.tick_params(**D.tickprops)
    mp.xaxis.set_major_formatter(D.fmtb)
    mp.set_title(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))

    plt.show()
    
    return fig
    
def PLOT_PLASMAVAR(time, vel_data, den, temp):    
    figVDT = plt.figure('PLASMAVAR'+str(random.randint(1,100))+str(random.randint(1,100)),figsize = (10,3))
    m1 = figVDT.add_subplot(311)
    #mp.plot(vel_time, vel_data, '+', label='v')
    m1.plot(time, vel_data, '+', label=r'$\mathrm{v_{tot}}$')
    m1.set_ylabel('v [km/s]',**D.lblpr)
    l11=m1.legend(**D.legprops)
    f11=l11.get_frame()
    f11.set_color('None')
    m1.tick_params(**D.tickprops)

    m2 = figVDT.add_subplot(312)
    m2.plot(time, den, '+', label='n')
    m2.set_ylabel(r'$\mathrm{n [cm^{-3}]}$',**D.lblpr)
    l2=m2.legend(**D.legprops)
    f2=l2.get_frame()
    f2.set_color('None')
    m2.tick_params(**D.tickprops)

    m3 = figVDT.add_subplot(313)
    m3.plot(time, temp, '+', label='T')
    m3.set_ylabel('T [MK]',**D.lblpr)
    l3=m3.legend(**D.legprops)
    f3=l3.get_frame()
    f3.set_color('None')
    m3.tick_params(**D.tickprops)

    figVDT.suptitle(str(probe)+(str(vel_time[0]))[:22]+'--'+str(vel_time[-1])[11:22]+'\n'+r'$\mathrm{\overline{v}}$, $\mathrm{\overline{n}}$,$\mathrm{\overline{T}}$ = %.2f %s, %.2f %s, %.2f %s'%(average_vel, '[m/s]', average_den, '$[\mathrm{cm^{-3}}]$', average_temp, '[MK]'), fontsize=12,  fontweight='bold')

    plt.show()
    #fig_name=str(output_path)+'VDTTS_'+str(probe)+'_'+str(vel_time[0])[0:10]+'_'+str(vel_time[0])[11:13]+str(vel_time[0])[14:16]+str(vel_time[0])[17:19]+'-'+str(vel_time[-1])[11:13]+str(vel_time[-1])[14:16]+str(vel_time[-1])[17:19]+'.png'
    #figVDT.savefig(fig_name)
    
def DRAW_ALL_PSD(freqs,PSDs,time,data):
    figPSD = plt.figure('ALL_PSD'+str(random.randint(1,100))+str(random.randint(1,100)),figsize = (10,7))
    p0=figPSD.add_axes([0.,0.,1.,1.])
    p0.axes.get_xaxis().set_visible(False)
    p0.axes.get_yaxis().set_visible(False)
    p = figPSD.add_axes([0.2,0.75,0.5,0.15])
    p.plot(time,data)
    p.xaxis.set_major_formatter(D.fmtb)
    p1 = figPSD.add_axes([0.2,0.2,0.5,0.5])
    for i in range(len(PSDs.keys())):
            p1.loglog(freqs,PSDs[i], ms=5, c=D.myc[i],label='var'+str(i+1))
            l=p1.legend()
    p1.set_xlabel('S/C frame frequency [Hz]'+'\n'+str(time[0])[0:10]+'\n'+str(time[0])[11:19]+'-'+str(time[-1])[11:19], **D.lblpr)  
    p1.set_ylabel(r'$\mathrm{PSD}$', **D.lblpr)
    p1.xaxis.grid(True, **D.mingridprops)
    p1.yaxis.grid(True, **D.majgridprops)
    p1.xaxis.grid(True, **D.majgridprops)
    #p.yaxis.grid(True, **D.mingridprops)
    plt.show()

    return figPSD


def DRAW_PDF(PDF_var, time, scl, frequency, bins_array, pdfs, flatness,spacecraft):    
    #Gauss Parameters
    #gp1=[1/(np.sqrt(2*np.pi)*STD_TS[0]), MEAN_TS[0], STD_TS[0]]
    gp1=[1/(np.sqrt(2*np.pi)), 0., 1.]
    #PDF x limits
    xlimit=max(abs(np.nanmin(PDF_var)),np.nanmax(PDF_var))
    #Gauss Bins
    gbi=np.arange(-xlimit, xlimit, 0.005)        

    drsc=[0,int(0.4*len(scl))+1,int(0.6*len(scl))+1,len(scl)-1]

    c=['b','m','g','y']
    col={}
    for ic in range(len(drsc)):
      col[drsc[ic]]=c[ic]
    fig=plt.figure('PDF'+str(random.randint(1,100))+str(random.randint(1,100)),figsize = (10,12)) 
    mf=fig.add_axes([0.,0.,1.,1.])
    mf.axes.get_xaxis().set_visible(False)
    mf.axes.get_yaxis().set_visible(False)
    m0=fig.add_axes([0.15, 0.8, 0.75,0.1])
    m0.plot(time, np.where(PDF_var>9999., -3.*np.mean(PDF_var[PDF_var<999.]), PDF_var))
    m0.set_ylabel('var',**D.lblpr)
    m0.xaxis.set_major_formatter(D.fmtb)

    m1=fig.add_axes([0.15, 0.35, 0.75,0.4])    
    for ih in drsc: 
      m1.semilogy(bins_array[ih], pdfs[ih], '+-', ms=12, mew=2., c=col[ih], label=r'$\tau$=%.1f [s]'%(scl[ih]/frequency))
    m1.semilogy(gbi, gauss(gbi, *gp1), 'k', lw=1.5)
    m1.set_xlabel(r'$\mathrm{\Delta}$', **D.lblpr)
    m1.set_ylabel('PDF', **D.lblpr)
    m1.set_xlim(-1.5*xlimit,1.5*xlimit)
    m1.set_ylim(ymin=1.e-2*np.nanmin(np.where(pdfs[0]==0., np.nan, pdfs[0])))
    l1=m1.legend(**D.legprops)
    l1.set_bbox_to_anchor((0.35,0.17))
    m1.xaxis.grid(True, **D.mingridprops)
    m1.yaxis.grid(True, **D.majgridprops)
    m1.xaxis.grid(True, **D.majgridprops)
    m1.yaxis.grid(True, **D.mingridprops)
    f1=l1.get_frame()
    f1.set_color('None')
    m1.tick_params(**D.tickprops)

    m2=fig.add_axes([0.15, 0.1, 0.75,0.17])
    m2.semilogx(np.divide(scl,frequency), flatness, 'pr-',**D.mpr2)
    m2.semilogx([scl[0]/frequency,scl[-1]/frequency],[3.,3.], 'k--',lw=2.)
    m2.set_xlabel(r'$\mathrm{\tau}$ [s]',**D.lblpr)
    m2.set_ylabel(r'Flatness',**D.lblpr)
    m2.set_xlim(scl[0]/frequency,scl[-1]/frequency)
    m2.set_ylim(2.,1.2*max(flatness))
    m2.tick_params(**D.tickprops)
    m2.yaxis.grid(True, **D.majgridprops)
    m2.xaxis.grid(True, **D.majgridprops)
    m2.xaxis.grid(True, **D.mingridprops) 
    fig.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))    
    plt.show()

    return fig
        
def DRAW_SF(PDF_var, time, scl, frequency, rank, structure_functions, flatness,spacecraft):    

    fig=plt.figure('SF'+str(random.randint(1,100))+str(random.randint(1,100)),figsize = (10,13)) 
    mf=fig.add_axes([0.,0.,1.,1.])
    mf.axes.get_xaxis().set_visible(False)
    mf.axes.get_yaxis().set_visible(False)        
    m0=fig.add_axes([0.15, 0.8, 0.65,0.1])
    m0.plot(time, np.where(PDF_var>9999., -3.*np.mean(PDF_var[PDF_var<999.]), PDF_var))
    m0.set_ylabel('var',**D.lblpr)
    m0.xaxis.set_major_formatter(D.fmtb)

    m1=fig.add_axes([0.15, 0.35, 0.65,0.4])    
    for ih in range(len(rank)): 
      m1.loglog(np.divide(scl,frequency), structure_functions[ih], 's-', ms=6, mec='None', label=r'$\mathrm{q}$=%i'%(rank[ih]))
    m1.set_xlabel(r'$\mathrm{\tau}$ [s]', **D.lblpr)
    m1.set_ylabel(r'$\mathrm{SF_q}$', **D.lblpr)
    m1.set_xlim(scl[0]/frequency,scl[-1]/frequency)
    l1=m1.legend(**D.legprops)
    #l1.set_bbox_to_anchor((0.35,0.17))
    m1.xaxis.grid(True, **D.mingridprops)
    m1.yaxis.grid(True, **D.majgridprops)
    m1.xaxis.grid(True, **D.majgridprops)
    m1.yaxis.grid(True, **D.mingridprops)
    f1=l1.get_frame()
    f1.set_color('None')
    m1.tick_params(**D.tickprops)

    m2=fig.add_axes([0.15, 0.1, 0.65,0.17])
    m2.semilogx(np.divide(scl,frequency), flatness, 'pr-',**D.mpr2)
    m2.semilogx([scl[0]/frequency,scl[-1]/frequency],[3.,3.], 'k--',lw=2.)
    m2.set_xlabel(r'$\mathrm{\tau}$ [s]',**D.lblpr)
    m2.set_ylabel(r'Flatness',**D.lblpr)
    m2.set_xlim(scl[0]/frequency,scl[-1]/frequency)
    m2.set_ylim(2.,1.2*max(flatness))
    m2.tick_params(**D.tickprops)
    m2.yaxis.grid(True, **D.majgridprops)
    m2.xaxis.grid(True, **D.majgridprops)
    m2.xaxis.grid(True, **D.mingridprops)     
    fig.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))        
    plt.show()

    return fig

def DRAW_PF(PF_var,time,scl,chosen_scl,frequency,rank,partition_functions,amplitude,TAUq,Dq,alpha,falpha,alpha_model,falpha_model,spacecraft):
    
    fig=plt.figure('PF'+str(random.randint(1,100))+str(random.randint(1,100)),figsize = (10,10)) 
    mf=fig.add_axes([0.,0.,1.,1.])
    mf.axes.get_xaxis().set_visible(False)
    mf.axes.get_yaxis().set_visible(False)
    m0=fig.add_axes([0.15, 0.83, 0.65,0.07])
    m0.plot(time, PF_var)
    m0.xaxis.set_major_formatter(D.fmtb)    
    m1=fig.add_axes([0.15, 0.45, 0.65,0.3])    
    for ih in range(0,len(rank),1): 
        m1.plot(np.log2(np.divide(scl,frequency)), np.log2(partition_functions[ih]), 's-', ms=6, mec='None', label=r'$\mathrm{q}$=%i'%(rank[ih]))
    for iu in range(0,len(rank),1):
        m1.plot(np.log2(np.divide(chosen_scl,frequency)), TAUq[iu]*np.log2(np.divide(chosen_scl,frequency))+amplitude[iu], 'k', linewidth=2.)                        
    m1.set_xlabel(r'$\mathrm{log_2 (l)}$', **D.lblpr)
    m1.set_ylabel(r'$\mathrm{log_2 \chi (q,l)}$', **D.lblpr)
    l1=m1.legend(**D.legprops)
    m2=fig.add_axes([0.15, 0.15, 0.25,0.2])    
    for ihi in range(len(rank)): 
        m2.plot(rank[ihi],TAUq[ihi], 'o-', ms=3, mec='r',mfc='None', mew=2)
    m2.set_xlabel('q', **D.lblpr)
    m2.set_ylabel(r'$\mathrm{\tau _q}$', **D.lblpr)
    m22=fig.add_axes([0.55, 0.15, 0.25,0.2])    
    for ihi in range(len(rank)): 
        m22.plot(rank[ihi],Dq[ihi], 's-', ms=3, mec='k',mfc='None', mew=2)            
    m22.set_xlabel('q', **D.lblpr)
    m22.set_ylabel(r'$\mathrm{D_q}$', **D.lblpr)
    fig.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))

    fig2=plt.figure('PFSpectrum',figsize = (7,7)) 
    mf=fig2.add_axes([0.,0.,1.,1.])
    mf.axes.get_xaxis().set_visible(False)
    mf.axes.get_yaxis().set_visible(False)
    m0=fig2.add_axes([0.15, 0.83, 0.65,0.07])
    m0.plot(time, PF_var)
    m0.xaxis.set_major_formatter(D.fmtb) 
    m1=fig2.add_axes([0.15, 0.15, 0.65,0.6])        
    m1.plot(alpha,falpha, '+', ms=3, mec='k',mfc='None', mew=2)
    m1.plot(alpha_model,falpha_model, 'r-')
    m1.set_xlabel(r'$\mathrm{\alpha (q)}$', **D.lblpr)
    m1.set_ylabel(r'$\mathrm{f(\alpha )}$', **D.lblpr)
    fig2.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))
    plt.show()
    
    return fig, fig2
    
    
def YInterval(deltaY, lowlim, highlim):
  DeltaY=np.logical_and(np.greater(deltaY,lowlim), np.less(deltaY,highlim))
  return DeltaY

def DRAW_ROMA(ROMA_TS,time,PDFs,PDFbn,R_PDF,R_PDFbn,DY,RL_SF_DY,slopes_DY,SF_ROMASol,ROMASols,chscl,data_freq,q_all,s_all,dy_ind,spacecraft):
    s_ind=list(s_all).index(ROMASols[dy_ind][0])
    
    DF_DY=np.empty((len(chscl)),dtype=object)
    binss=np.empty((len(chscl)),dtype=object)
    for iscl in range(len(chscl)):
        DF=YInterval(PDFbn[iscl], DY[dy_ind][0]*((chscl[iscl]/chscl[0])**ROMASols[dy_ind][0]), DY[dy_ind][-1]*((chscl[iscl]/chscl[0])**ROMASols[dy_ind][0]))
        binss[iscl]=np.extract(DF,PDFbn[iscl])
        DF_DY[iscl]=np.extract(DF,PDFs[iscl])             
        
    ###Draw fit for the first and last order q
    fitq1s1=np.polyfit(np.log10(chscl),np.log10(RL_SF_DY[dy_ind][0][s_ind]),1)
    polq1s1=np.poly1d(fitq1s1)
    logfitq1s1=lambda x: np.power(10, polq1s1(np.log10(x)))

    fitq2s1=np.polyfit(np.log10(chscl),np.log10(RL_SF_DY[dy_ind][-1][s_ind]),1)
    polq2s1=np.poly1d(fitq2s1)
    logfitq2s1=lambda x: np.power(10, polq2s1(np.log10(x)))
    
    fig=plt.figure('ROMA_PDF'+str(random.randint(1,100))+str(random.randint(1,100)),figsize=(10,10))
    mf=fig.add_axes([0.,0.,1.,1.])
    mf.axes.get_xaxis().set_visible(False)
    mf.axes.get_yaxis().set_visible(False)
    m0=fig.add_axes([0.1, 0.55, 0.33, 0.33])
    for i in range(len(chscl)):
        m0.semilogy(PDFbn[i], PDFs[i], 's-', ms=2.1,lw=2., color=D.myc[i], label=r'$\tau$=%.2f s'%(chscl[i]/data_freq))
    m0.set_xlabel(r'$\mathrm{\delta B}$ [nT]', **D.lblpr)
    m0.set_ylabel(r'$\mathrm{P(\delta B,\tau)}$', **D.lblpr)
    m0.yaxis.grid(True, **D.mingridprops)
    m0.xaxis.grid(True, **D.majgridprops)
    legpr = dict(loc='best', frameon=True, borderaxespad=0.)    
    lp1=m0.legend(**legpr)
    lp1.get_frame().set_color('None')
    m0.set_title('(a)')        

        
    m1=fig.add_axes([0.6, 0.55, 0.33, 0.33])
    for iplot in range(len(q_all)):
        m1.loglog(np.divide(chscl,data_freq), (RL_SF_DY[dy_ind][iplot][s_ind]), 'D--', markersize=8, label='q=%.1f'%q_all[iplot])
    m1.loglog(np.divide(chscl,data_freq), logfitq1s1(chscl))
    m1.loglog(np.divide(chscl,data_freq), logfitq2s1(chscl))
    m1.set_xlabel(r'$\tau$ [s]', **D.lblpr)
    m1.set_ylabel(r'SF($\tau$,q)', **D.lblpr)
    m1.tick_params(axis='both', which='major', labelsize=12)
    l1=m1.legend(loc='best',frameon=False,borderaxespad=0.)
    m1.set_title(r'$\mathrm{(b)\ \ \Delta Y_{%i}=[%.2f,%.2f]}$'%(dy_ind+1,DY[dy_ind][0],DY[dy_ind][-1]))        

    m2=fig.add_axes([0.1, 0.12, 0.33, 0.33])
    for qid in range(len(q_all)):
        m2.plot(s_all, slopes_DY[dy_ind][qid], 'o', ms=4, mec=D.myc[qid] , c=D.myc[qid], label='q=%.1f'%q_all[qid])
        m2.plot(s_all, q_all[qid]*s_all, c=D.myc[qid])
    m2.set_xlabel('s', **D.lblpr)
    m2.set_ylabel(r'$\zeta$(q)', **D.lblpr)
    l2=m2.legend(**D.legprops)
    l2.get_frame().set_color('None') 
    m2.set_title(r'$\mathrm{(c)\ \ \Delta Y_{%i}=[%.2f,%.2f]}$'%(dy_ind+1,DY[dy_ind][0],DY[dy_ind][-1]))    

    m3=fig.add_axes([0.6, 0.12, 0.33, 0.33])      
    for jplot in range(len(DY)):
        if ROMASols[jplot][0]>0.0:
            m3.plot(q_all, SF_ROMASol[jplot],  marker=D.mym[jplot], c=D.myc[jplot], lw=0., label=r'$\mathrm{\Delta}$ Y=[%.2f,%.2f], s=%.2f'%(DY[jplot][0],DY[jplot][-1],ROMASols[jplot][0]))
            m3.plot(q_all, q_all*ROMASols[jplot][0], ls='--', c=D.myc[jplot])
    m3.set_xlabel('q',**D.lblpr)
    m3.set_ylabel(r'$\zeta_{q}$',**D.lblpr)
    m3.tick_params(axis='both', which='major', labelsize=12)        
    l3=m3.legend(loc=2, frameon=False, borderaxespad=0.)
    m3.set_title('(d)')    
    fig.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))        
    
    figROMA=plt.figure('RescaledPDFs'+str(random.randint(1,100))+str(random.randint(1,100)),figsize=(10,9))    
    mr=figROMA.add_axes([0.,0.,1.,1.])
    mr.axes.get_xaxis().set_visible(False)
    mr.axes.get_yaxis().set_visible(False)    
    p4=figROMA.add_axes([0.1,0.8,0.8,0.1])
    p4.plot(time,ROMA_TS)
    p4.xaxis.set_major_formatter(D.fmtb)
    p4.set_title('(a)')
    
    p2=figROMA.add_axes([0.1, 0.48, 0.8, 0.25])    
    for jplot in range(len(DY)):
        if ROMASols[jplot]!=None:
            for isol in range(len(ROMASols[jplot])):
                p2.plot(DY[jplot],[ROMASols[jplot][isol], ROMASols[jplot][isol]], c='k', lw=2.5)
                p2.plot(DY[jplot],[ROMASols[jplot][isol], ROMASols[jplot][isol]], c=D.myc[isol], lw=1.9)
                p2.plot(DY[jplot][0], ROMASols[jplot][isol], marker='|', markeredgecolor='k', markerfacecolor='k')
                p2.plot(DY[jplot][-1], ROMASols[jplot][isol], marker='|', markeredgecolor='k', markerfacecolor='k')
    p2.set_xlim(DY[0][0], DY[-1][-1])
    p2.set_ylim(0.,1.)
    p2.set_xlabel('Y [nT]', fontsize=12)
    p2.set_ylabel('s(Y)', fontsize=12)
    p2.tick_params(axis='both', which='major', labelsize=12) 
    p2.set_title('(b)')
    
    m11=figROMA.add_axes([0.1,0.15,0.35,0.25])
    for i in range(0,len(chscl),1):    
        m11.semilogy(PDFbn[i], PDFs[i],'s-',lw=0.3, ms=3.,mfc='None',mew=1.2, mec=D.myc[i],color=D.myc[i], label=r'$\tau$=%.2f s'%(chscl[i]/data_freq))
        m11.semilogy(binss[i], DF_DY[i],lw=4., c=D.myc[dy_ind])
    m11.xaxis.grid(True, **D.majgridprops)   
    m11.yaxis.grid(True, **D.majgridprops)
    m11.set_xlim(xmin=0.)    
    m11.set_xlabel(r'$\mathrm{\delta B}$ [nT]', **D.lblpr)
    m11.set_ylabel(r'$\mathrm{P(\delta B,\tau)}$', **D.lblpr)
    l11=m11.legend(loc='best', frameon=True, borderaxespad=0.)
    l11.get_frame().set_color('None')    
    m11.set_title('(c)')
    
    p3=figROMA.add_axes([0.55,0.15,0.35,0.25])
    for ix in range(len(DY)):  
        for iy in range(len(chscl)):
            p3.semilogy(R_PDFbn[ix][iy], R_PDF[ix][iy], marker=D.mym[iy], mew=2.,ms=3, mfc='None',lw=0.,mec=D.myc[ix])
        p3.semilogy([DY[ix][-1],DY[ix][-1]], [1.e-3,10.], '--', lw=1.2, color=D.myc[ix])
        p3.semilogy([-DY[ix][-1],-DY[ix][-1]], [1.e-3,10.], '--', lw=1.2, color=D.myc[ix])
    p3.set_xlabel('Y [nT]', fontsize=12)
    p3.set_ylabel(r'$\mathrm{P_s(Y) [nT^{-1}]}$', fontsize=12)
    llg=[]
    for i in range(len(DY)):
        llg.append('$\mathrm{s_{\Delta Y%i}}$=%3.2f'%(i+1,ROMASols[i][0])) 
    titl1 =figROMA.add_axes([0.56, 0.16, 0.35, 0.25])
    for i in range(len(DY)):
        titl1.text(0.7, 0.85-0.1*i, llg[i], color=D.myc[i], **D.txtpr)
    titl1.axis('off')
    p3.set_title('(d)')    
    figROMA.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))        
    plt.show()   
    
    return fig, figROMA

def DRAW_MI(MI,MI_Baseline,scales,resolution,compBaseline,time,spacecraft):
    if compBaseline:
        v1 = np.logspace(np.log10(MI.min()), np.log10(MI.max()), 5, endpoint=True)
        vb = np.logspace(np.log10(np.min(np.nanmin(MI_Baseline))), np.log10(np.max(np.nanmax(MI_Baseline))), 5, endpoint=True)

        fig = plt.figure('MI and Baseline'+str(random.randint(1,100))+str(random.randint(1,100)), figsize = (10,4))
        mf=fig.add_axes([0.,0.,1.,1.])
        mf.axes.get_xaxis().set_visible(False)
        mf.axes.get_yaxis().set_visible(False)
        m0=fig.add_axes([0.15, 0.25, 0.33, 0.6])
        im0=m0.imshow(MI, interpolation='bicubic', origin='lower', norm=mpl.colors.LogNorm(vmin=MI.min(), vmax=MI.max()), cmap=plt.cm.jet)
        cbar0=fig.colorbar(im0,ticks=v1,format='$%.2f$')

        m0.set_ylabel(r'$\tau$ [s]', **D.lblpr)
        m0.set_xlabel(r'$\tau$ [s]', **D.lblpr)
        cbar0.set_label('MI',**D.lblpr)
        ticks=np.arange(0,len(scales))
        plt.xticks(ticks)
        plt.yticks(ticks)
        labels = [item.get_text() for item in m0.get_xticklabels()]
        for il in range(len(labels)):
            labels[il]=scales[il]*resolution
        m0.set_xticklabels(['%.2f'% a for a in labels], rotation=90, fontsize=12)
        m0.set_yticklabels(['%.2f'% a for a in labels], fontsize=12)


        m1=fig.add_axes([0.55, 0.25, 0.33, 0.6])
        im1=m1.imshow(list(MI_Baseline), interpolation='bicubic', origin='lower', norm=mpl.colors.LogNorm(vmin=np.min(np.nanmin(MI_Baseline)), vmax=np.max(np.nanmax(MI_Baseline))), cmap=plt.cm.jet)
        cbar1=fig.colorbar(im1,ticks=vb,format='$%.2f$')
        cbar1.set_label(r'$\mathrm{\frac{MI-\overline{MI}}{\sigma_{MI}}}$',**D.lblpr)
        ticks=np.arange(0,len(scales))
        plt.xticks(ticks)
        plt.yticks(ticks)
        labels = [item.get_text() for item in m1.get_xticklabels()]
        for il in range(len(labels)):
            labels[il]=scales[il]*resolution
        m1.set_xticklabels(['%.2f'% a for a in labels], rotation=90, fontsize=12)
        m1.set_yticklabels(['%.2f'% a for a in labels], fontsize=12) 
        fig.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))
    else:
        v1 = np.logspace(np.log10(MI.min()), np.log10(MI.max()), 5, endpoint=True)

        fig = plt.figure('MI', figsize = (6,6))
        mf=fig.add_axes([0.,0.,1.,1.])
        mf.axes.get_xaxis().set_visible(False)
        mf.axes.get_yaxis().set_visible(False)
        m0=fig.add_axes([0.25, 0.25, 0.65, 0.65])
        im0=m0.imshow(MI, interpolation='bicubic', origin='lower', norm=mpl.colors.LogNorm(vmin=MI.min(), vmax=MI.max()), cmap=plt.cm.jet)
        cbar0=fig.colorbar(im0,ticks=v1,format='$%.2f$')

        m0.set_ylabel(r'$\tau$ [s]', **D.lblpr)
        m0.set_xlabel(r'$\tau$ [s]', **D.lblpr)
        cbar0.set_label('MI',**D.lblpr)
        ticks=np.arange(0,len(scales))
        plt.xticks(ticks)
        plt.yticks(ticks)
        labels = [item.get_text() for item in m0.get_xticklabels()]
        for il in range(len(labels)):
            labels[il]=scales[il]*resolution
        m0.set_xticklabels(['%.2f'% a for a in labels], rotation=90, fontsize=12)
        m0.set_yticklabels(['%.2f'% a for a in labels], fontsize=12)
        fig.suptitle(spacecraft+'  '+str(time[0])+' - '+str(time[-1]))
    plt.show()

    return fig
################################################################


################################################################
####        TRIM DATA - SELECT DATA SUBSET      ################
################################################################
def CHOP_TS(ti, tf, DATA_TIME, DATA):
    d=DATA_TIME[0].date()
    #if str(ti)[0]!=0 and str(tf)[0]!=0:
    t_start=time(int(str(ti)[:2]),int(str(ti)[2:4]), int(str(ti)[4:6]))
    t_stop=time(int(str(tf)[:2]),int(str(tf)[2:4]), int(str(tf)[4:6]))

    start=datetime.combine(DATA_TIME[0].date(),t_start)
    stop=datetime.combine(DATA_TIME[0].date(),t_stop)
    start_ind = bisect.bisect_left(DATA_TIME, start)
    stop_ind = bisect.bisect_left(DATA_TIME, stop) 
    
    if stop_ind-start_ind > 0:    
        return DATA_TIME[start_ind:stop_ind], DATA[start_ind:stop_ind,:]
    else:
        print ('WARNING!!! time series could not be trimmed, limits out of bounds')
        return [],[]
###############################################################


################################################################
########    PRE_PROCESS DATA ALGORITHMS      ###################
################################################################    
    
def PRINT_FLAGS(flag, time, DATA):
    ### Just printing the erroneous data (flagged, nans) and whether there are any data gaps in time
    flag_indices=np.transpose(np.where((DATA>=flag) | (np.isnan(DATA))))
    print ('flagged_indices', flag_indices)
    ######## DATA GAPS
    mag_freq_array=np.subtract(time[1:],time[:-1])
    ts=[]
    for it in range(len(mag_freq_array)):
        ts.append(mag_freq_array[it].total_seconds())
    print ('mean time resolution (arithmetic mean): ', np.nanmean(ts), '\nmedian time resolution (middle value): ', np.median(ts), '\nmost frequent time resolution: ', mode(ts))
    idx_arr=np.where(mag_freq_array>timedelta(seconds=0.045))[0]
    print ('time index with resolution larger than nominal: ', idx_arr, len(idx_arr))

    
def MASK_ERRDATA(flag, time, DATA):    
    MSK_MAG_DATA=np.ma.masked_greater(DATA,flag)
    MSK_MAG_DATA=np.ma.masked_invalid(MSK_MAG_DATA)
    #FILLED_MAG_DATA=FILLED_MAG_DATA.filled(np.nan)

    mask=np.ma.getmaskarray(MSK_MAG_DATA[:,0])

    DATA_MASKED=[]#np.empty((DATA.shape[1]),dtype=object)
    for i in range(DATA.shape[1]):
      DATA_MASKED.append(MSK_MAG_DATA[:,i][~mask])

    DATA_MASKED=np.array(DATA_MASKED)
    DATA_MASKED=np.swapaxes(DATA_MASKED,0,1)
    t_MASKED=np.array(time)[~mask]
    
    ts=[]
    mtdiff = np.array(t_MASKED) - datetime(1970,1,1)
    for isec in range(len(t_MASKED)):
      ts.append(mtdiff[isec].total_seconds())

    return t_MASKED, ts, DATA_MASKED 


def SELF_FLAGG(resolution, time_in_sec, DATA_MASKED):

    ##PSEUDO - INTERPOLATION - fill any data-gaps with flagg_vales
    #PSEUDOINTERP_MAG_DATA=MAG_DATA_MASKED
    t_pseudointerp=[]
    var_pseudointerp=[ [] for i in range(DATA_MASKED.shape[1]) ]
    
    jj=0
    for iint in range(len(time_in_sec)-1):
      if abs(time_in_sec[iint+1]-time_in_sec[iint])<=resolution+1.e-3:
        t_pseudointerp.append(time_in_sec[jj])
        for ii in range(DATA_MASKED.shape[1]):
            var_pseudointerp[ii].append(DATA_MASKED[jj,ii])
      else:
        t_pseudointerp.append(time_in_sec[jj])  
        for ii in range(DATA_MASKED.shape[1]):
            var_pseudointerp[ii].append(DATA_MASKED[jj,ii])
        ss=1
        for nb in range(int(math.ceil(abs(time_in_sec[iint+1]-time_in_sec[iint])/resolution))-1):
          t_pseudointerp.append(time_in_sec[iint]+ss*resolution)
          for ii in range(DATA_MASKED.shape[1]):
              var_pseudointerp[ii].append(1.e5)
          ss=ss+1
      jj=jj+1

    tparsed=[]
    for tpars in range(len(t_pseudointerp)):
      tparsed.append(datetime.utcfromtimestamp(t_pseudointerp[tpars]))
    var_pseudointerp=np.array(var_pseudointerp)
    var_pseudointerp=np.swapaxes(var_pseudointerp,0,1)
    
    return t_pseudointerp, tparsed, var_pseudointerp   


def STATIONARITY_TEST(flagged_field,iter_param):
    field={}
    for jf in flagged_field.keys():
        field[jf] = np.where(flagged_field[jf]>=9999., np.nan, flagged_field[jf])
    partial_mean_components, total_mean_components, partial_std_components, total_std_components={},{},{},{}
    for ik in field.keys():
        partial_mean,partial_std=[],[]
        total_mean_components[ik]=np.nanmean(field[ik])
        total_std_components[ik]=np.nanstd(field[ik])
        #print ('component mean and std: ', ik,mean_component_field,std_component_field)
        for ipar in iter_param:
            partial_mean.append(np.nanmean(field[ik][0:int(ipar*len(field[ik]))])/total_mean_components[ik])
            partial_std.append(np.nanstd(field[ik][0:int(ipar*len(field[ik]))])/total_std_components[ik])
        partial_mean_components[ik]=partial_mean
        partial_std_components[ik]=partial_std
    return  partial_mean_components, partial_std_components, total_mean_components, total_std_components   
###############################################################


################################################################
########    ANALYSIS METHODS ALGORITHMS      ###################
################################################################       
def PSD(DATA_PSEUDOINTERP, freq, window, segment_magnitude, overlap):
    #segm_length_minutes=segment_magnitude/(freq*60.)
    PSDs={}
    
    for i in range(DATA_PSEUDOINTERP.shape[1]):
        frequencies,PSD = signal.welch(DATA_PSEUDOINTERP[:,i], freq, window, segment_magnitude, overlap)
        PSDs[i]=PSD
    return frequencies, PSDs
 
    
def PSEUDOINTERP_PSD(frequency, masked_time_in_sec, time_of_flagged_data, DATA_MASKED, window, segment_magnitude, overlap_percent,draw_PSD=False, lbls=None):
    if lbls==None:
        lbls=[]
        if len(DATA_MASKED.shape)>1:
            for il in range(DATA_MASKED.shape[1]):
                lbls.append('var'+str(il+1))
        elif len(DATA_MASKED.shape)==1:        
            lbls.append('var0')    
    
    soverlap=int(overlap_percent*segment_magnitude)
    segm_length_min=segment_magnitude/(frequency*60.)
    print ('Welch segment length = %.2f minutes'%segm_length_min)

    #INTERPOLATE ONLY DATA GAPS LARGER THAN 2 TIMES THE RESOLUTION IN THE TIME SERIES
    DATA_PSEUDOINTERP = [[] for i in range(DATA_MASKED.shape[1])]        
    for i in range(DATA_MASKED.shape[1]):
        interp_fct = interp1d(masked_time_in_sec, DATA_MASKED[:,i])
        interp_data = interp_fct(time_of_flagged_data)
        DATA_PSEUDOINTERP[i]=interp_data
    DATA_PSEUDOINTERP=np.array(DATA_PSEUDOINTERP)
    DATA_PSEUDOINTERP=np.swapaxes(DATA_PSEUDOINTERP,0,1)
        
    #PSD CALCULATION
    freqs, power_spectral_density = PSD(DATA_PSEUDOINTERP, frequency, window, segment_magnitude, soverlap)  

    #DRAW
    if draw_PSD:
        figs=[]
        for ifig in range(len(lbls)):
                figPSD = plt.figure(lbls[ifig],figsize = (10,7))
                p0=figPSD.add_axes([0.,0.,1.,1.])
                p0.axes.get_xaxis().set_visible(False)
                p0.axes.get_yaxis().set_visible(False)
                p = figPSD.add_axes([0.2, 0.2, 0.7, 0.7])
                p.loglog(freqs,power_spectral_density[ifig],'b', label='PSD '+ lbls[ifig])
                l=p.legend()
                p.set_xlabel('S/C frame frequency [Hz]', **D.lblpr)  
                p.set_ylabel(r'$\mathrm{PSD}$', **D.lblpr)
                #p.set_ylabel(r'$\mathrm{PSD \ [nT^2 Hz^{-1}]}$', **D.lblpr) 
                p.xaxis.grid(True, **D.mingridprops)
                p.yaxis.grid(True, **D.majgridprops)
                p.xaxis.grid(True, **D.majgridprops)
                #p.yaxis.grid(True, **D.mingridprops)
                figs.append(figPSD)
        #plt.show()        
        return DATA_PSEUDOINTERP, freqs, power_spectral_density, figs
    else:
        return DATA_PSEUDOINTERP, freqs, power_spectral_density        

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def maxscale(ts_length):
    for i in range(0,100):
      if not 2**i<ts_length:
        break
    return i-1 

def Scales(s1, s2, fineRes):
  s = np.arange(s1, s2+1, 1, dtype=int)
  pow2_scl = 2**s

  if fineRes:
    intermediate_scl=np.add(2**s[:-2], 2**s[1:-1])
    final_scales=np.concatenate((pow2_scl, intermediate_scl), axis=0)
    final_scales.sort()
    fin_scales=[]
    for i in range(len(final_scales)):
        fin_scales.append(int(final_scales[i]))
    return fin_scales
  else:
    fin_scales=[]
    for i in range(len(pow2_scl)):
        fin_scales.append(pow2_scl[i])
    #return pow2_scl
    return fin_scales

def PDF(PDF_var, rank, scl, bin_no, compute_SF):
    DifferenceTS, AbsDifferenceTS, histos, binc_array={},{}, {}, {}
    Flatness=[]
    SFTOT=np.empty((len(rank),len(scl)),dtype=float)
    for s in range(len(scl)):
        ### COMPUTE TIME SERIES OF DIFFERENCES
        DifferenceTS[s]= np.subtract(PDF_var[scl[s]:len(PDF_var)],PDF_var[0:len(PDF_var)-scl[s]])
        DifferenceTS[s]=DifferenceTS[s][(np.abs(DifferenceTS[s])<9999.)& (np.abs(DifferenceTS[s])!=0.)]
        ###NON_STANDARDIZED
        AbsDifferenceTS[s]=np.abs(DifferenceTS[s])
        ###STANDARDIZED
        #AbsDifferenceTS[s]=np.divide(np.abs(DifferenceTS[s]), np.nanstd(DifferenceTS[s]))
        DifferenceTS[s]=np.divide(np.subtract(DifferenceTS[s],np.nanmean(DifferenceTS[s])), np.nanstd(DifferenceTS[s]))
        #DifferenceTS[s][abs(DifferenceTS[s])>9999.]=np.nan
        Flatness.append(np.nanmean(np.power(DifferenceTS[s],4))/(np.nanmean(np.power(DifferenceTS[s],2)))**2)
        
        ### COMPUTE HISTOGRAMS OF THE DIFFERENCES
        minbin=(np.floor(min(DifferenceTS[s])*10))/10
        maxbin=(np.ceil(max(DifferenceTS[s])*10))/10
        hist, bin_edges = np.histogram(DifferenceTS[s], bins=bin_no, range=(minbin,maxbin), density=False)
        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
        hist=hist/(sum(hist)*(bin_centers[2]-bin_centers[1]))
        histos[s]=hist
        binc_array[s]=bin_centers
        
        ### COMPUTE STRUCTURE FUNCTIONS
        if compute_SF:
            for q in rank:
                SFTOT[list(rank).index(q),s]=np.nanmean(np.power(AbsDifferenceTS[s],q))    

    if compute_SF:
        return histos, binc_array, Flatness, SFTOT
    elif not compute_SF:
        return histos, binc_array, Flatness, None      

#def Wavelet()  
####PARTITION FUNCTIONS####
def PartitionFunctions(PDF_var, frequency, rank, scl, choose_scl): 
    #COMPUTE THE PARTITION FUNCTIONS
    print ('COMPUTING PARTITION FUNCTIONS: IT MAY TAKE A WHILE')
    DifferenceTS, AbsDifferenceTS,partial_probability,sum_D={},{},{},{}
    partition_functions=np.zeros((len(rank),len(scl)),dtype=float)
    for s in range(len(scl)):
        ### COMPUTE TIME SERIES OF DIFFERENCES
        DifferenceTS[s]= np.subtract(PDF_var[scl[s]:len(PDF_var)],PDF_var[0:len(PDF_var)-scl[s]])
        DifferenceTS[s]=DifferenceTS[s][(np.abs(DifferenceTS[s])<9999.)& (np.abs(DifferenceTS[s])!=0.)]
        ###NON_STANDARDIZED
        AbsDifferenceTS[s]=np.abs(DifferenceTS[s])
        sum_D[s]=np.nansum(AbsDifferenceTS[s])
        lens=len(AbsDifferenceTS[s])
        for qj in range(len(rank)):
            prob_q,pfq=0,[]        
            for j in range(0,math.floor(lens/scl[s])):
                part_prob=np.nansum(np.divide(AbsDifferenceTS[s][j*scl[s]:(j+1)*scl[s]],sum_D[s]))               
                prob_q=part_prob**rank[qj]
                pfq.append(prob_q)
            partition_functions[qj,s]=np.nansum(pfq) 
            
    #PARTITION FUNCTIONS FIT        
    chosen_fit_scl=scl[choose_scl[0]:choose_scl[-1]]
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y: (y - fitfunc(p, x))
    pinit = [1.0, -1.0]
    powerlaw = lambda x, amp, index: amp * (x**index) 

    #FIT THE PARTITION FUNCTIONS AND COMPUTE TAUq
    out, param_final, tau_q, amp={}, {}, [], []
    for iq in range(0,len(rank)):
        out[iq] = optimize.leastsq(errfunc, pinit, args=(np.log2(np.divide(chosen_fit_scl,frequency)), np.log2(partition_functions[iq][choose_scl[0]:choose_scl[-1]])), full_output=1)
        param_final[iq] = out[iq][0]
        amp.append(param_final[iq][0])
        tau_q.append(param_final[iq][1])

    #COMPUTE THE GENERALIZED DIMENSIONS
    Dq={}
    for tq in range(0,len(rank)):
            Dq[tq]=(tau_q[tq]/(rank[tq]-1))
    for key, value in Dq.items():
        if value>9999.:
            Dq[key] = np.nan

    #COMPUTE THE MULTIFRACTAL SPECTRUM f(alpha)
    alpha=np.diff(tau_q)/np.diff(rank)
    f_alpha=np.multiply(rank[:-1],alpha)-tau_q[:-1]            

    #SUPERPOSE p-model ON THE f(alpha) PLOT
    p=0.4
    rank2=np.arange(-10, 10, 0.1)
    alpha_model, Dq_model= [],[]
    for iqm in range(len(rank2)):
        alpha_model.append(((p**rank2[iqm]) * np.log2(p) + ((1-p)**rank2[iqm])*np.log2(1-p)) /  (p**rank2[iqm] + (1-p)**rank2[iqm])*np.log2(0.5))
        #if rank[iqm]!=1:
        Dq_model.append((np.log2((p**rank2[iqm])+ (1-p)**rank2[iqm] ))/((rank2[iqm]-1)*np.log2(0.5)))
    falpha_model=np.subtract(np.multiply(rank2,alpha_model),np.multiply((np.subtract(rank2,1)),Dq_model))

    
    
    return partition_functions,tau_q,amp,Dq,alpha,f_alpha,alpha_model,falpha_model
    

#####ROMA####
def Fluctuations(bvec, scl,rank):
    DY_all=np.empty((len(scl)), dtype=list)
    bin_no=200
    mx=[]
    for i in range(0,len(scl)):
        Dy=[]
        for j in range(0,len(bvec)-scl[i]):
            if abs(bvec[j+scl[i]])>9999. or abs(bvec[j])>9999. or (abs(bvec[j+scl[i]])>9999. and abs(bvec[j])>9999.):
                continue
            else:
                Dy.append(abs(bvec[j+scl[i]]-bvec[j]))    
                DY_all[i]=Dy
        mx.append(np.nanmax(DY_all[i]))
    return DY_all, mx

def DY_limits(DYmax):
    limits=[0.05, 0.15*DYmax]
    while 1.2*limits[-1]<DYmax:
        limits.append(limits[-1]+1.2*(limits[-1]-limits[-2]))
    if DYmax-limits[-1]<0.3*(limits[-2]-limits[-3]):
        limits[-1]=DYmax
    elif DYmax-limits[-1]>=0.3*(limits[-2]-limits[-3]) and DYmax-limits[-1]<0.6*(limits[-2]-limits[-3]):
        limits.append(2*DYmax-limits[-1]) 
    else:
        limits.append(DYmax)

    return limits



def RangeLimited_SF(DY_all, DeltaY, q_rank, scl, s_roma):
  print ('this is the function which calculates the Range-Limited Structure Functions, SF(tau,q)')
  SFq=[]#np.empty((len(q_rank)))
  RankFluct=np.empty((len(s_roma),len(scl)), dtype=object)
  for q in q_rank:
    SF=np.empty((len(s_roma),len(scl))) ##if array is 'empty' all entries must be filled!!, careful with this!
    for iscl in range(0,len(scl)):
        n=0
        for sr in s_roma:
          FluctDY1=[]
          SF_scl_q_s=0
          DeltaY1=YInterval(np.divide(DY_all[iscl], (scl[iscl]/scl[0])**sr), DeltaY[0], DeltaY[-1])
          FluctDY1=np.extract(DeltaY1,DY_all[iscl])
          SF_scl_q_s=(np.nanmean(np.power(FluctDY1,q)))
          SF[n,iscl]=SF_scl_q_s
          n=n+1
    SFq.append(SF) ##array with q*s*tau dimensions

  return SFq

def FindSolutions(slope, q, s_roma):
  differences=np.subtract(slope, q*s_roma)
  j=0
  solROMA=[]
  for di in range(len(differences)-1):
    solutions=0
    if (differences[di]*differences[di+1]) > 0:
      continue
    elif differences[di]*differences[di+1] <= 0:
      if abs(differences[di])<abs(differences[di+1]):
        solutions=s_roma[di]
      elif abs(differences[di])>=abs(differences[di+1]):
        solutions=s_roma[di+1]
    solROMA.append(solutions)
  #print 'solutions are: ', solutions, solROMA[0]
  return solROMA

def PDF_DY(bvec, scl, sol_roma):
  #minbin=-17.#(floor(min(Dy)*10))/10
  #maxbin= 17.#(ceil(max(Dy)*10))/10
  bin_no=200
  DY_all=np.empty((len(scl)), dtype=list)
  PDFDY_all=np.empty((len(scl)), dtype=object)
  histos=np.empty((len(scl)), dtype=object)
  binc_array=np.empty((len(scl)), dtype=object)
  for i in range(0,len(scl)):
      #print 'SCALE IS:', i , scl[i]
      Dy=[]
      mx=[]
      for j in range(0,len(bvec)-scl[i]):
        if abs(bvec[j+scl[i]])>9999. or abs(bvec[j])>9999. or (abs(bvec[j+scl[i]])>9999. and abs(bvec[j])>9999.):
          continue
        else:
          Dy.append(bvec[j+scl[i]]-bvec[j])
          #print 'len(Dy)', Dy      
          DY_all[i]=Dy

      minbin=(np.floor(min(Dy)*10))/10
      maxbin=(np.ceil(max(Dy)*10))/10
     
      
      PDFDY_all[i]=np.divide(DY_all[i], (scl[i]/scl[0])**sol_roma)   
      hist, bin_edges = np.histogram(PDFDY_all[i], bins=bin_no, range=(minbin,maxbin), density=False)
      bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
      #hist=np.where(hist>=30.,hist,0.) # valorile mai mici de 30 sunt inlocuite cu 0
      hist=hist/(sum(hist)*(bin_centers[2]-bin_centers[1]))
      histos[i]=hist
      binc_array[i]=bin_centers
      
  return DY_all, histos, binc_array 

def ROMA(ROMA_timeseries,chscl,binno,s_all,q_all,given_Y_bins):
    ScaleFluctuations, MaxScaleFluctuations= Fluctuations(ROMA_timeseries,chscl,q_all)  

    '''Define Y bins, either automatically or as they have been chosen '''    
    if given_Y_bins==None:
        max_Yedge = MaxScaleFluctuations[0]
        DYlim=DY_limits(max_Yedge)
        Y_bins=[]
        for di in range(len(DYlim)-1):
            Y_bins.append([DYlim[di], DYlim[di+1]])
    else:
        Y_bins=given_Y_bins
            
    '''RANGE-LIMITED STRUCTURE FUNCTION CALCULATION'''
    RL_SF_allY=np.empty((len(Y_bins)),dtype=list)
    #RankFluct=np.empty((len(Y_bins)),dtype=object)
    for idy in range(len(Y_bins)):
        #RangeLimited_SF=[]
        RL_SF=RangeLimited_SF(ScaleFluctuations, Y_bins[idy], q_all, chscl, s_all)
        RL_SF_allY[idy]=RL_SF
    print ('done computing range limited structure functions')   
    
    '''Identify all possible solutions'''
    slROMA=np.empty((len(q_all),len(Y_bins)), dtype=list)
    slopes_allq=np.empty((len(q_all),len(s_all)), dtype=list)
    RL_SF_slopes_allDY=np.empty((len(Y_bins)),dtype=list)

    for dly in range(len(Y_bins)):
        slopes_ALLQ=[]
        #resid_ALLQ=[]
        for q_index in range(len(q_all)):
            q_toDraw=q_all[q_index]
            for fi in range(len(s_all)):
                fitq1,resid=None,None
                slopes_allq[q_index,fi]=None
                fitq1,resid,_,_,_=np.polyfit(np.log10(chscl),np.log10(RL_SF_allY[dly][q_index][fi]),1, full=True)
                slopes_allq[q_index,fi]=fitq1[0]
            slopes_ALLQ.append(list(slopes_allq[q_index]))
            slROMA[q_index,dly]=FindSolutions(slopes_allq[q_index], q_toDraw, s_all)
        RL_SF_slopes_allDY[dly]=slopes_ALLQ    

    all_solutionsDY=np.empty((len(Y_bins)), dtype=list) ##all solutions for all q for each DY
    sROMA_setsY=np.empty((len(Y_bins)), dtype=object) ## sets of independent solutions defined above
    for isy in range(len(Y_bins)):
        all_s_DY=[]
        for isq in range(len(q_all)):
            for iss in range(len(slROMA[isq][isy])):
                if slROMA[isq][isy][iss]>0.:
                    all_s_DY.append(slROMA[isq][isy][iss])
            all_solutionsDY[isy]=all_s_DY ## all solutions for one DY
            sROMA_setsY[isy]=set(all_solutionsDY[isy]) ##set of independent solutions for one DY

    all_sROMA_ind=np.empty((len(Y_bins)), dtype=list)
    for idY in range(len(Y_bins)):
        all_s_index=[]
        for isi in range(len(sROMA_setsY[idY])):
            all_s_index.append(list(s_all).index(list(sROMA_setsY[idY])[isi]))
        all_sROMA_ind[idY]=all_s_index

    print ('solutions sets for each DY', sROMA_setsY)    
  
    '''Find the best solution'''
    print ('Find the Best ROMA Solution for each DY')
    sfq_allsolDY=np.empty((len(Y_bins)), dtype=object)
    chi_2_DY = np.empty((len(Y_bins)), dtype=object) 
    for iDY in range(len(Y_bins)):
        sfq_allsol=np.empty((len(all_sROMA_ind[iDY])), dtype=list)
        chi_2 = np.empty((len(all_sROMA_ind[iDY])), dtype=list)
        fit_sol=np.empty((len(all_sROMA_ind[iDY])), dtype=list)
        for jas in range(len(all_sROMA_ind[iDY])): #there are as many solutions as I have q
            sfq=[]
            for qrk in range(len(q_all)):
                sfq.append(RL_SF_slopes_allDY[iDY][qrk][all_sROMA_ind[iDY][jas]])
            sfq_allsol[jas]=sfq
            fitsfq=np.polyfit(q_all,sfq_allsol[jas],1)
            residuals=np.subtract(sfq_allsol[jas], np.multiply(q_all,list(sROMA_setsY[iDY])[jas]))
            sol_err=np.std(sfq_allsol[jas])
            chi_2[jas]=np.nansum(np.power(residuals,2))
        sfq_allsolDY[iDY]=sfq_allsol
        chi_2_DY[iDY]=chi_2

    ROMASol=np.empty((len(Y_bins)), dtype=object)
    for idly in range(len(Y_bins)):
        if chi_2_DY[idly].size:
            GoodSol=np.array(chi_2_DY[idly]==np.nanmin(np.array(chi_2_DY[idly])))
            chosen_sol=np.array(list(sROMA_setsY[idly]))[GoodSol]
            ROMASol[idly]=chosen_sol
        else:
            print ("CHI2 IS EMPTY: ", chi_2_DY[idly])  
            ROMASol[idly]=np.array([0.])#chosen_sol#np.array([0.])
    ok_sol=(ROMASol==ROMASol)
    ROMASol=ROMASol[ok_sol]
    print ('ROMASol: ', ROMASol)    

    '''Choose only those Structure Functions that correspond to the best ROMA solutions'''
    RL_SF_ROMASol={}
    for dyi in range(len(Y_bins)):
        if ROMASol[dyi][0]>0.0:
            RL_SF_ROMASol[dyi]=sfq_allsolDY[dyi][list(sROMA_setsY[dyi]).index(ROMASol[dyi][0])]    
    
    '''SCALED PDFs'''
    _,PDF,PDF_bins=PDF_DY(ROMA_timeseries,chscl, 0.)    
    RangeLimited_PDF=np.empty((len(Y_bins)),dtype=object)
    bincen=np.empty((len(Y_bins)),dtype=object)
    RL_PDF=np.empty((len(Y_bins), len(chscl)), dtype=object)
    RL_PDF_bins=np.empty((len(Y_bins), len(chscl)), dtype=object)
    for idypdf in range(len(Y_bins)):
        _,RangeLimited_PDF[idypdf], bincen[idypdf]=PDF_DY(ROMA_timeseries,chscl, ROMASol[idypdf][0])
        for ipdf in range(len(chscl)):
            OK_PDF=YInterval(abs(bincen[idypdf][ipdf]), Y_bins[idypdf][0], Y_bins[idypdf][-1])
            RL_PDF[idypdf, ipdf]=np.extract(OK_PDF,RangeLimited_PDF[idypdf][ipdf])
            RL_PDF_bins[idypdf, ipdf]=np.extract(OK_PDF,bincen[idypdf][ipdf])

    return PDF,PDF_bins,RL_PDF,RL_PDF_bins,Y_bins,RL_SF_allY,RL_SF_slopes_allDY,RL_SF_ROMASol,ROMASol   
##########END ROMA   
    
########## MI ###############    
def EqualLengthFluctuations(fluctuations_TS, scl):
    DifferenceTS={}
    for s in range(len(scl)):
        ### COMPUTE TIME SERIES OF DIFFERENCES
        DifferenceTS[s]=np.subtract(fluctuations_TS[int(scl[-1]/2)+int(scl[s]/2):-scl[-1]+int(scl[s]/2)], fluctuations_TS[int(scl[-1]/2)-int(scl[s]/2):-scl[-1]-int(scl[s]/2)])
        DifferenceTS[s]=np.abs(DifferenceTS[s])
        DifferenceTS[s][abs(DifferenceTS[s])>9999.]=np.nan
    return DifferenceTS
    
def MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    c_xy = np.rot90(c_xy)
    c_xy = np.flipud(c_xy)
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi

def MI_Baseline(diffVector, scales, ref_scales, binnb, noScr=3):
    scrVect=copy.deepcopy(diffVector)
    Mean_MI=np.empty((len(ref_scales)), dtype=object)
    Std_MI=np.empty((len(ref_scales)), dtype=object)
    for refs in range(len(scales)):
        MI_tautau=np.empty((noScr, len(scales)))
        Lambda_tautau=np.empty((noScr, len(scales)))
        for ifms in range(0,noScr):
            for imi in range(len(scales)):
                np.random.shuffle(scrVect[imi])
                MI_tautau[ifms,imi]=MI(diffVector[refs], scrVect[imi], binnb)
        Mean_MI[refs]=list(np.mean(MI_tautau, axis=0))
        Std_MI[refs]=list(np.std(MI_tautau, axis=0))
    return Mean_MI, Std_MI
    
###############################################################  
