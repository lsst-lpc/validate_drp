# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.

from __future__ import print_function, division

import matplotlib.pylab as plt
import numpy as np
import scipy.stats
from .check import fitExp, fitAstromErrModel, fitPhotErrModel, expModel, astromErrModel, photErrModel


# Plotting defaults
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 20
# plt.rcParams['figure.titlesize'] = 30

color = {'all': 'grey', 'bright': 'blue',
         'iqr': 'green', 'rms': 'red'}


def plotOutlinedLinesHorizontal(ax, *args, **kwargs):
    """Plot horizontal lines outlined in white.

    The motivation is to let horizontal lines stand out clearly
    even against a cluttered background.
    """
    plotOutlinedLines(ax.axhline, *args, **kwargs)


def plotOutlinedLinesVertical(ax, *args, **kwargs):
    """Plot vertical lines outlined in white.

    The motivation is to let horizontal lines stand out clearly
    even against a cluttered background.
    """
    plotOutlinedLines(ax.axvline, *args, **kwargs)


def plotOutlinedLines(ax_plot, x1, x2, x1_color=color['all'], x2_color=color['bright']):
    """Plot horizontal lines outlined in white.

    The motivation is to let horizontal lines stand out clearly
    even against a cluttered background.
    """
    ax_plot(x1, color='white', linewidth=4)
    ax_plot(x2, color='white', linewidth=4)
    ax_plot(x1, color=x1_color, linewidth=3)
    ax_plot(x2, color=x2_color, linewidth=3)



def plotVisitVsTime(goodMatches,
                    outputPrefix=""):
    mjd=[]
    visit=[]
    for group in goodMatches.groups:
        mjd += list(group.get('MJD-OBS'))
        visit += list(group.get('visit'))
    plt.figure(figsize=(10,8))
    time= mjd-min(mjd)
    plt.scatter(time, visit )
    plt.xlabel('t-tmin (mjd)')
    plt.ylabel('visit')
    plt.title('Visit in function of Time (mjd)')
    plotPath = outputPrefix + 'VisitVsTime.png'
    plt.savefig(plotPath, format="png")

    plt.figure(figsize=(12,10))
    plt.scatter(mjd, visit )
    plt.xlabel('t (mjd)')
    plt.ylabel('visit')
    plt.title('Visit in function of Time (mjd)')
    plotPath = outputPrefix + 'VisitVsTime_mjd.png'
    plt.savefig(plotPath, format="png")
     #plt.show()



def plotAstromPhotRMSvsTimeCcd(dist, mag, snr, goodMatches, mmagrms,
                               fit_params=None, brightSnr=100, srcFluxField='base_PsfFlux',
                               outputPrefix="",
                               zoom=False):

    sourceFluxField=srcFluxField

    if zoom:
        outputPrefix=outputPrefix+"_zoom_"

    sizelegend=12 # taille des legendes
    digits=1000. # precision des valeurs dans les legendes des histos
  
  #  print("LONGUEUR MAGRMS", len(mmagrms))
  #  print("LONGUEUR dist", len(dist))
  #  print(" MAGRMS", mmagrms)
    compt = 0
    goodmjd = []
    radToDeg = 180./np.pi
    degToArcs = 3600.
    radToArcs = radToDeg* degToArcs

    deltaRAcosdecs = [] #goodMatches.aggregate(np.ones,'coord_ra')
    deltaDecs = [] # goodMatchesaggregate(np.ones,'coord_dec')

    ccds = []
    sourcemag = []
    sourcedmag = []
    sourcesnr = []
    #dist=[]
    dist_calc = []
    FluxMag0s = []
    FluxMag0Errs = []

    grpMeanRAcosdec = []
    grpMeanDec = []
    groupRMSracosdec = []
    groupRMSdec = []
    posRMS = []

    medsnr = []
    medsnrlong = []
    Psf_fwhm = []
    grpMeanPsf_fwhm = []
   # grpMeanShapex  = []### test shape

  #  BRIGHTS=['bright', 'brightmean']

  #  for suffix in BRIGHTS:
  #      exec('psf_fwhm'+'_'+suffix +"=[]")


    for group in goodMatches.groups:
      #  group_schema=group.getSchema()
      #  print('group_schema=group.getSchema()',group_schema.getOrderedNames())
      #  print('')
        
        
        snr = group.get((sourceFluxField+"_snr"))
        sourcesnr += list(snr)
        medsnr.append(np.median(snr))
        medsnrlong += list(np.ones(len(snr))*np.mean(snr))
        # print('medsnrlong',medsnrlong)
        RA = group.get('coord_ra')
        Dec = group.get('coord_dec')
        
    #    brightmean,= np.where(np.asarray(medsnr) > brightSnr)
     #   print(' brightmean,', brightmean)

   #     bright, = np.where(np.asarray(snr) > brightSnr)
   #     print(' bright,', bright)

     #   print('test base_SdssShape_x', group.get('base_SdssShape_x'))
     #   Shapex = group.get('base_SdssShape_x')
      #  MeanShapex = np.mean( Shapex)
      #  grpMeanShapex.append(MeanShapex)
 #  print('test base_SdssShape_y', group.get('base_SdssShape_y'))


        psf_fwhm = group.get('PSF-FWHM')
        Psf_fwhm+=list(psf_fwhm)
        MeanPsf_fwhm = np.mean(psf_fwhm)
        grpMeanPsf_fwhm.append(MeanPsf_fwhm)

      #  print(' psf_fwhm', psf_fwhm)
        MeanRA = np.mean(RA)
        MeanDec = np.mean(Dec)

        dra = ((RA -  MeanRA) * np.cos(MeanDec)) *radToArcs*1000
        ddec = (Dec - MeanDec) *radToArcs*1000
        
        deltaRAcosdecs += list(dra)
        deltaDecs += list(ddec)
        
        grpMeanRAcosdec.append(MeanRA  * np.cos(MeanDec))
        grpMeanDec.append(MeanDec)
        
        groupRMSracosdec.append( np.std(dra))
        groupRMSdec.append( np.std(ddec))
        posRMS.append(np.sqrt((np.std(dra))**2 + (np.std(ddec))**2))
        
        goodmjd += list(group.get('MJD-OBS'))
        sourcemag += list((group.get(sourceFluxField+"_mag"))*1000)
        sourcedmag += list((group.get(sourceFluxField+"_mag") - np.mean(group.get(sourceFluxField+"_mag")))*1000)

        FluxMag0s += list(group.get('FLUXMAG0'))
        FluxMag0Errs += list(group.get('FLUXMAG0ERR'))

        ccds += list(group.get('ccd'))


        # for suffix in BRIGHTS:
        #     var_dra='dra'+'_'+suffix
        #     var_ddec='ddec'+'_'+suffix
        #     var_grpMeanRAcosdec=' grpMeanRAcosdec'+'_'+suffix
        #     var_grpMeanDec='grpMeanDec'+'_'+suffix
      


        # test grandes snr
       # if medsnr>100:
       #     nbn=5
       #     plt.figure()
       #     plt.hist( dra, bins=nbn)
       #     plt.title('(ra-ramean)cosdec')
       #     plt.figure()
       #     plt.hist(ddec, bins=nbn)
       #     plt.title('(dec-decmean)')
       #     plt.figure()
       #     plt.hist(snr, bins=nbn)
       #     plt.title('snr')
       #     plt.figure()
       #     plt.hist(psf_fwhm, bins=nbn)
       #     plt.title(' psf_fwhm')
       #     plt.show()
            
    

    
    bright, = np.where(np.asarray(medsnr) > brightSnr)
   # bright2, = np.where(np.asarray(snr) > brightSnr) # ???  # vient de goodPsfSnr = goodMatches.aggregate(np.median, field=psfSnrKey) dans validate.py# 

    groupRMSracosdec_bright = np.array(groupRMSracosdec)[bright]
    groupRMSdec_bright = np.array(groupRMSdec)[bright]
    posRMS_bright = np.array(posRMS)[bright]

    nb_sigma=5
    bright_outliers, = np.where((np.asarray(posRMS_bright-np.median(posRMS_bright)) < nb_sigma*np.std(posRMS_bright))) and  np.where((np.asarray(groupRMSracosdec_bright-np.median(groupRMSracosdec_bright)) < nb_sigma*np.std(groupRMSracosdec_bright))) and np.where((np.asarray(groupRMSdec_bright-np.median(groupRMSdec_bright)) < nb_sigma*np.std(groupRMSdec_bright)))

    posRMS_bright_outliers = posRMS_bright[bright_outliers]
    groupRMSracosdec_bright_outliers = groupRMSracosdec_bright[bright_outliers]
    groupRMSdec_bright_outliers = groupRMSdec_bright[bright_outliers]

    grpMeanRAcosdec=np.array(grpMeanRAcosdec)

    deltaRAcosdecs = np.array(deltaRAcosdecs)
    deltaDecs = np.array(deltaDecs)
    Psf_fwhm = np.array(Psf_fwhm)

    #brightallsnr, = np.where(np.asarray(sourcesnr) > brightSnr)
    brightallsnr, = np.where(np.asarray(medsnrlong) > brightSnr) #pour avoir le meme cut sur les snr medianes


    plt.close('all')

    plt.figure(figsize=(12,12))
    plt.title('PSF FWHM')
    plt.hist(Psf_fwhm, histtype ='stepfilled', alpha=0.8, color='b')
    histB = plt.hist(Psf_fwhm[brightallsnr] ,histtype ='stepfilled',alpha=0.8,color='r')
    plt.axvline(np.median(Psf_fwhm), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(Psf_fwhm)*digits)/digits)+'as')
    plt.axvline(np.median(Psf_fwhm[brightallsnr]), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(Psf_fwhm[brightallsnr] )*digits)/digits)+'as')
 
    plt.xlabel('psf fwhm (as)')
    plt.ylabel('# / bin')
    plt.legend(prop={'size':sizelegend})
    if zoom:
        plt.xlim(0.,55.)
        plt.ylim(0., max( histB[0]))
    plotPath = outputPrefix+'PsfFwhmvshist.png'
    plt.savefig(plotPath, format="png")





    plt.figure(figsize=(12,12))
    plt.title('PSF FWHM vs deltaRAcosdecs')
    plt.scatter(deltaRAcosdecs, Psf_fwhm, color=color['all'], label='all')
    plt.scatter(deltaRAcosdecs[brightallsnr], Psf_fwhm[brightallsnr], color=color['bright'], label='bright')
    plt.legend()
    plt.xlabel('DeltaRaCosDec (mas)')
    plt.ylabel('psf fwhm (as)')
    if zoom:
        plt.xlim(-75.,75.)
        plt.ylim(min(Psf_fwhm),max(Psf_fwhm))
    plotPath = outputPrefix+'PsfFwhmvsdeltaRacosDec.png'
    plt.savefig(plotPath, format="png")

    plt.figure(figsize=(12,12))
    plt.title('PSF FWHM vs deltaDecs')
    plt.scatter(deltaDecs, Psf_fwhm, color=color['all'], label='all')
    plt.scatter(deltaDecs[brightallsnr], Psf_fwhm[brightallsnr], color=color['bright'], label='bright')
    plt.legend()
    plt.xlabel('DeltaDec (mas)')
    plt.ylabel('psf fwhm (as)')
    if zoom:
        plt.xlim(-75.,75.)
        plt.ylim(min(Psf_fwhm),max(Psf_fwhm))
    plotPath = outputPrefix+'PsfFwhmvsdeltaDec.png'
    plt.savefig(plotPath, format="png")

    #plt.show()

    plt.close('all')
    plt.figure()
    plt.title('racosdec')
    plt.scatter(grpMeanRAcosdec, groupRMSracosdec, color=color['all'], label='all')
    plt.scatter(grpMeanRAcosdec[bright], groupRMSdec_bright, color=color['bright'], label='bright')
    plt.scatter((grpMeanRAcosdec[bright])[bright_outliers], (groupRMSdec_bright)[bright_outliers], color='r', label='bright outliers')
    plt.legend()
    plt.xlabel('mean ra cosdec')
    plt.ylabel('rms ra cos dec')
    plotPath = outputPrefix+'RAcosDecRMSvsRacosDecMean.png'
    plt.savefig(plotPath, format="png")

    plt.figure()
    plt.title('racosdec')
    plt.scatter(grpMeanRAcosdec[bright], groupRMSdec_bright, color=color['bright'], label='bright')
    plt.scatter((grpMeanRAcosdec[bright])[bright_outliers], (groupRMSdec_bright)[bright_outliers], color='r', label='bright outliers')

    plt.legend()
    plt.xlabel('mean ra cosdec')
    plt.ylabel('rms ra cos dec')

    plotPath = outputPrefix+'RAcosDecRMSvsRacosDecMean_bright.png'
    plt.savefig(plotPath, format="png")


# grpMeanShapex

   # plt.figure()
   # plt.title('test corr shape_x /RA')
   # plt.scatter(groupRMSracosdec,grpMeanShapex, color=color['all'], label='all')
    # #marche pas plt.scatter(groupRMSracosdec[bright], grpMeanShapex[bright], color=color['bright'], label='bright')#
    #plt.legend()
   # plt.xlabel('mean ra cosdec')
   # plt.ylabel('Shape_x')
    """
    plt.figure()
    plt.title('racosdec')
    plt.scatter(grpMeanDec, groupRMSracosdec, color=color['all'], label='all')
    plt.scatter(grpMeanDec[bright], groupRMSdec_bright, color=color['bright'], label='bright')
    plt.scatter((grpMeanRAcosdec[bright])[bright_outliers], (groupRMSdec_bright)[bright_outliers], color='r', label='bright outliers')
    plt.legend()

    plt.legend()
    plt.xlabel('dec')
    plt.ylabel('rms racos')
    """
  #  plt.show()
    
    #plot(x,y,"k.")
  #  y_av = movingaverage(groupRMSracosdec, 100)
   # plt.plot(grpMeanRAcosdec, y_av,"r")
  #  plt.xlim(0,1000)
   # plt.grid(True)
   # plt.show()


    plt.figure(figsize=(12,12))
    plt.title('RMS racosdec vs RMS dec')
    plt.scatter(groupRMSracosdec, groupRMSdec, alpha=0.8, color='b', label='all')
    plt.scatter(groupRMSracosdec_bright, groupRMSdec_bright, alpha=0.8, color='r', label='bright')
    plt.scatter(groupRMSracosdec_bright_outliers, groupRMSdec_bright_outliers, alpha=0.8, color='g', label='bright + 5 sigma outliers')
    plt.xlabel('RMS RAcosDec (mas)')
    plt.ylabel('RMS Dec (mas)')
    plt.xlim(0,max(groupRMSracosdec))
    plt.ylim(0,max(groupRMSdec))
    if zoom:
        plt.xlim(0.,55.)
        plt.ylim(0.,55.)
    plt.legend(prop={'size':sizelegend})
    plotPath = outputPrefix+'RAcosDecRMSvsDecRMS.png'
    plt.savefig(plotPath, format="png")
    #plt.show()


    plt.figure()
    plt.title('racosdec')
    plt.hist(groupRMSracosdec, bins=50, histtype ='stepfilled', alpha=0.8, color='b')# label='RMS='+str(int(np.std(groupRMSracosdec)*digits)/digits)+'mAcs\nMean='+str(int(np.mean(groupRMSracosdec)*digits)/digits)+'mAcs',alpha=0.5)#,histtype ='stepfilled',alpha=0.8,color='r')
    histB = plt.hist(groupRMSracosdec_bright,histtype ='stepfilled',alpha=0.8,color='r')
    plt.axvline(np.median(groupRMSracosdec), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(groupRMSracosdec)*digits)/digits)+'mas')
    plt.axvline(np.median(groupRMSracosdec_bright), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(groupRMSracosdec_bright)*digits)/digits)+'mas')
    plt.hist(groupRMSracosdec_bright_outliers ,histtype ='stepfilled',alpha=0.5,color='g')
    plt.axvline(np.median( groupRMSracosdec_bright_outliers), 0, 1, linewidth=2, color='green', label='Median 5sigm clipp='+str(int(np.median( groupRMSracosdec_bright_outliers )*digits)/digits)+'mas')
    plt.xlabel('RMS RAcosDec')
    plt.ylabel('#/bin')
    plt.legend(prop={'size':sizelegend})
    if zoom:
        plt.xlim(0.,55.)
        plt.ylim(0., max( histB[0]))
    plotPath = outputPrefix+'RAcosDecRMS.png'
    plt.savefig(plotPath, format="png")

    plt.figure()
    plt.title('dec')
    plt.hist(groupRMSdec, bins=50, histtype ='stepfilled', alpha=0.8, color='b')
    histB = plt.hist(groupRMSdec_bright,histtype ='stepfilled', alpha=0.8,color='r')
    plt.hist(groupRMSdec_bright_outliers ,histtype ='stepfilled',alpha=0.5,color='g')
    plt.axvline(np.median(groupRMSdec), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(groupRMSdec)*digits)/digits)+'mas')
    plt.axvline(np.median(groupRMSdec_bright), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(groupRMSdec_bright)*digits)/digits)+'mas')
    plt.axvline(np.median( groupRMSdec_bright_outliers), 0, 1, linewidth=2, color='green', label='Median 5sigm clipp='+str(int(np.median( groupRMSdec_bright_outliers )*digits)/digits)+'mas')
    plt.xlabel('RMS Dec')
    plt.ylabel('#/bin')
    plt.legend(prop={'size':sizelegend})
    if zoom:
        plt.xlim(0.,55.)
        plt.ylim(0., max( histB[0]))
    plotPath = outputPrefix+'DecRMS.png'
    plt.savefig(plotPath, format="png")

    plt.figure()
    plt.title('RMS distances (equiv plotAstrometry)')
    plt.hist(posRMS, bins=50, histtype ='stepfilled', alpha=0.8, color='b')
    histB = plt.hist(posRMS_bright,histtype ='stepfilled',alpha=0.8,color='r')
    plt.hist(posRMS_bright_outliers ,histtype ='stepfilled',alpha=0.5,color='g')
    plt.axvline(np.median(posRMS), 0, 1, linewidth=2, color='blue', label='Median='+str(int(np.median(posRMS)*digits)/digits)+'mas')
    plt.axvline(np.median(posRMS_bright), 0, 1, linewidth=2, color='red', label='Median bright='+str(int(np.median(posRMS_bright)*digits)/digits)+'mas')
    plt.axvline(np.median( posRMS_bright_outliers), 0, 1, linewidth=2, color='green', label='Median 5sigm clipp='+str(int(np.median( posRMS_bright_outliers )*digits)/digits)+'mas')

    plt.xlabel('RMS distance')
    plt.ylabel('#/bin')
    plt.legend(prop={'size':sizelegend})
    if zoom:
        plt.xlim(0.,55.)
        plt.ylim(0., max( histB[0]))
    plotPath = outputPrefix+'DistanceRMS.png'
    plt.savefig(plotPath, format="png")
   # plt.show()
    print('lenccd', len(ccds))
    print('lenddec',len(deltaDecs))
    print('len fluxmag0',len(FluxMag0s))
    #     FluxMag0Errs
   # plt.figure() #test
   # plt.plot( FluxMag0s)
  #  plt.title('FluxMag0')
  #  plt.show()
    if zoom:
        return

  #  print('ccd',ccds)
    print("LONGUEUR DELTARAS", len(deltaRAcosdecs))
    print("LONGUEUR SOURCEMAG", len(sourcemag))
    print("LONGUEUR sourceSNR", len(sourcesnr))
    #faire un script pour sortir la med et rms en fonction du temps

    times=[]
    dic_time={}
    for i in range(len(deltaRAcosdecs)):
         dic_time[goodmjd[i]]={}
         dic_time[goodmjd[i]]['dra']=[]
         dic_time[goodmjd[i]]['ddec']=[]
         dic_time[goodmjd[i]]['sourcedmag']=[]
         dic_time[goodmjd[i]]['sourcesnr']=[]

         dic_time[goodmjd[i]]['FluxMag0s']=[]
         dic_time[goodmjd[i]]['FluxMag0Errs']=[]


    for i in range(len(deltaRAcosdecs)):
         dic_time[goodmjd[i]]['dra'].append(deltaRAcosdecs[i])
         dic_time[goodmjd[i]]['ddec'].append(deltaDecs[i])
         dic_time[goodmjd[i]]['sourcedmag'].append(sourcedmag[i])
         dic_time[goodmjd[i]]['sourcesnr'].append(sourcesnr[i])

         dic_time[goodmjd[i]]['FluxMag0s'].append(FluxMag0s[i])
         dic_time[goodmjd[i]]['FluxMag0Errs'].append(FluxMag0Errs[i])


    for time in dic_time.keys():
         times.append(time)
         dic_time[time]['dra_med'] = np.median(dic_time[time]['dra'])
         dic_time[time]['dra_rms'] = np.std(dic_time[time]['dra'])

         dic_time[time]['ddec_med'] = np.median(dic_time[time]['ddec'])
         dic_time[time]['ddec_rms'] = np.std(dic_time[time]['ddec'])

         dic_time[time]['sourcedmag_med'] = np.median(dic_time[time]['sourcedmag'])
         dic_time[time]['sourcedmag_rms'] = np.std(dic_time[time]['sourcedmag'])

         dic_time[time]['sourcesnr_med'] = np.median(dic_time[time]['sourcesnr'])
         dic_time[time]['sourcesnr_rms'] = np.std(dic_time[time]['sourcesnr']) # barre err inutile
         dic_time[time]['FluxMag0s_mean'] = np.mean(dic_time[time]['FluxMag0s'])
         dic_time[time]['FluxMag0s_med'] = np.median(dic_time[time]['FluxMag0s'])
         dic_time[time]['FluxMag0s_rms'] = np.std(dic_time[time]['FluxMag0s'])
         dic_time[time]['FluxMag0Errs_mean'] = np.mean(dic_time[time]['FluxMag0Errs'])
         dic_time[time]['FluxMag0Errs_med'] = np.median(dic_time[time]['FluxMag0Errs'])
         dic_time[time]['FluxMag0Errs_rms'] = np.std(dic_time[time]['FluxMag0Errs'])
 
       #  plt.errorbar(time,  dic_time[time]['dra_med'], xerr=None, yerr=dic_time[time]['dra_rms'],linestyle='',alpha=0.75, marker='o', color='b')
     #    plt.errorbar(time,  dic_time[time]['ddec_med'], xerr=None, yerr=dic_time[time]['ddec_rms'],linestyle='',alpha=0.75, marker='o', color='r') 
    #plt.figure(figsize=(10,8))
 
    sizelegend=12
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False, figsize=(16,9))
    # f.subplots_adjust(hspace=0) # pour enlever l'espace entre les graphiques (vertical)
    # plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    for time in dic_time.keys():
         ax1.scatter(time, dic_time[time]['dra_rms'],alpha=0.85, marker='o', color='b')
         ax1.scatter(time, dic_time[time]['ddec_rms'],alpha=0.85, marker='o', color='r')
    ax1.scatter([], [] , color='b',marker='o',label='$\delta(RA) \cos(dec)$')
    ax1.scatter([], [] , color='r',marker='o',label='$\delta(Dec)$')
    ax1.legend(prop={'size':sizelegend})
    ax1.set_xlim(min(goodmjd)-(max(goodmjd)-min(goodmjd))/30., max(goodmjd)+(max(goodmjd)-min(goodmjd))/30.)
    ax1.set_ylabel('RMS (mas)')
   # ax1.set_xlabel('time (mjd)')
    ax1.set_title('RMS Dra/Ddec in function of time')
    for time in dic_time.keys():
         ax2.scatter(time, dic_time[time]['sourcedmag_rms'],alpha=0.85, marker='o', color='g') #
    ax2.scatter([], [], marker='o', color='g',label='source $\Delta$Magnitude RMS') #
    ax2.set_xlim(min(goodmjd)-(max(goodmjd)-min(goodmjd))/30., max(goodmjd)+(max(goodmjd)-min(goodmjd))/30.)
    ax2.set_ylabel('RMS $\Delta$Mag (mmag)')
    #plt.xlabel('time (mjd)')
    ax2.set_title('RMS mag in function of time')
    ax2.legend(prop={'size':sizelegend})

    for time in dic_time.keys():
         ax3.scatter(time, dic_time[time]['sourcesnr_med'],alpha=0.85, marker='o', color='k') #
        # ax3.errorbar(time,  dic_time[time]['sourcesnr_med'], xerr=None, yerr=dic_time[time]['sourcesnr_rms'],linestyle='',alpha=0.75, marker='o', color='b')
    ax3.scatter([], [], marker='o', color='b',label='source SNR Median') #
    ax3.set_xlim(min(goodmjd)-(max(goodmjd)-min(goodmjd))/30., max(goodmjd)+(max(goodmjd)-min(goodmjd))/30.)
    ax3.set_ylabel('SNR')
    ax3.set_xlabel('time (mjd)')
    ax3.set_title('Median SNR in function of time')
    ax3.legend(prop={'size':sizelegend})
    plotPath = outputPrefix+'AstrometryVsTime.png'
    plt.savefig(plotPath, format="png")

    #### astrometry and photometry vs CCD ####
    dic_ccd={}
    for i in range(len(deltaRAcosdecs)):
         dic_ccd[ccds[i]]={}
         dic_ccd[ccds[i]]['dra']=[]
         dic_ccd[ccds[i]]['ddec']=[]
         dic_ccd[ccds[i]]['sourcedmag']=[]
         dic_ccd[ccds[i]]['sourcesnr']=[]
         dic_ccd[ccds[i]]['FluxMag0s']=[]
         dic_ccd[ccds[i]]['FluxMag0Errs']=[]

    for i in range(len(deltaRAcosdecs)):
         dic_ccd[ccds[i]]['dra'].append(deltaRAcosdecs[i])
         dic_ccd[ccds[i]]['ddec'].append(deltaDecs[i])
         dic_ccd[ccds[i]]['sourcedmag'].append(sourcedmag[i])
         dic_ccd[ccds[i]]['sourcesnr'].append(sourcesnr[i])
         dic_ccd[ccds[i]]['FluxMag0s'].append(FluxMag0s[i])
         dic_ccd[ccds[i]]['FluxMag0Errs'].append(FluxMag0Errs[i])

    for ccd in dic_ccd.keys():
         ccds.append(ccd)
         dic_ccd[ccd]['dra_med'] = np.median(dic_ccd[ccd]['dra'])
         dic_ccd[ccd]['dra_rms'] = np.std(dic_ccd[ccd]['dra'])

         dic_ccd[ccd]['ddec_med'] = np.median(dic_ccd[ccd]['ddec'])
         dic_ccd[ccd]['ddec_rms'] = np.std(dic_ccd[ccd]['ddec'])

         dic_ccd[ccd]['sourcedmag_med'] = np.median(dic_ccd[ccd]['sourcedmag'])
         dic_ccd[ccd]['sourcedmag_rms'] = np.std(dic_ccd[ccd]['sourcedmag'])

         dic_ccd[ccd]['sourcesnr_med'] = np.median(dic_ccd[ccd]['sourcesnr'])
         dic_ccd[ccd]['sourcesnr_rms'] = np.std(dic_ccd[ccd]['sourcesnr'])# barre err totalement inutile
 
         dic_ccd[ccd]['FluxMag0s_mean'] = np.mean(dic_ccd[ccd]['FluxMag0s'])
         dic_ccd[ccd]['FluxMag0s_med'] = np.median(dic_ccd[ccd]['FluxMag0s'])
         dic_ccd[ccd]['FluxMag0s_rms'] = np.std(dic_ccd[ccd]['FluxMag0s'])
         dic_ccd[ccd]['FluxMag0Errs_mean'] = np.mean(dic_ccd[ccd]['FluxMag0Errs'])
         dic_ccd[ccd]['FluxMag0Errs_med'] = np.median(dic_ccd[ccd]['FluxMag0Errs'])
         dic_ccd[ccd]['FluxMag0Errs_rms'] = np.std(dic_ccd[ccd]['FluxMag0Errs'])
 
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False, figsize=(16,9))

    for ccd in dic_ccd.keys():
         ax1.scatter(ccd, dic_ccd[ccd]['dra_rms'],alpha=0.85, marker='o', color='b')
         ax1.scatter(ccd, dic_ccd[ccd]['ddec_rms'],alpha=0.85, marker='o', color='r')
    ax1.scatter([], [] , color='b',marker='o',label='$\delta(RA) \cos(dec)$')
    ax1.scatter([], [] , color='r',marker='o',label='$\delta(Dec)$')
    ax1.legend(prop={'size':sizelegend})
    ax1.set_xlim(min(ccds)-(max(ccds)-min(ccds))/30., max(ccds)+(max(ccds)-min(ccds))/30.)
    ax1.set_ylabel('RMS (mas)')
    ax1.set_title('RMS Dra/Ddec in function of CCD')

    for ccd in dic_ccd.keys():
         ax2.scatter(ccd, dic_ccd[ccd]['sourcedmag_rms'],alpha=0.85, marker='o', color='g') #
    ax2.scatter([], [], marker='o', color='g',label='source $\Delta$Magnitude RMS') #
    ax2.set_ylabel('RMS $\Delta$Mag (mmag)')
    ax2.set_title('RMS mag in function of CCD')
    ax2.legend(prop={'size':sizelegend})

    for ccd in dic_ccd.keys():
         ax3.scatter(ccd, dic_ccd[ccd]['sourcesnr_med'],alpha=0.85, marker='o', color='k') #
         # ax3.errorbar(ccd,  dic_ccd[ccd]['sourcesnr_med'], xerr=None, yerr=dic_ccd[ccd]['sourcesnr_rms'],linestyle='',alpha=0.75, marker='o', color='b')
    ax3.scatter([], [], marker='o', color='b',label='SNR Median') #
    ax3.set_ylabel('SNR')
    ax3.set_xlabel('CCD #')
    ax3.set_title('Median SNR in function of CCD')
    ax3.legend(prop={'size':sizelegend})
    plotPath = outputPrefix+'AstrometryVsCcd.png'
    plt.savefig(plotPath, format="png")
    s=100
    plt.figure( figsize=(14,5))
    plt.title('FluxMag0s vs ccd')
    plt.xlabel('ccd')
    plt.ylabel('FluxMag0s')
    for ccd in dic_ccd.keys():
        plt.scatter(ccd, dic_ccd[ccd]['FluxMag0s_med'],s=s,alpha=0.95, marker='+', color='k') 
        plt.scatter(ccd, dic_ccd[ccd]['FluxMag0s_mean'],s=s,alpha=0.95, marker='+', color='r') 
        #plt.scatter(ccd, dic_ccd[ccd]['FluxMag0Errs_med'],alpha=0.85, marker='o', color='k') 

    plt.scatter([], [], marker='+',s=s, color='k',label='Median')
    plt.scatter([], [], marker='+',s=s, color='r',label='Mean')
    plt.legend()
    plotPath = outputPrefix+'ZpVsCCD.png'
    plt.savefig(plotPath, format="png")

    plt.figure( figsize=(14,5))
    plt.title('FluxMag0s vs time')
    plt.xlabel('time')
    plt.ylabel('FluxMag0s')
    for time in dic_time.keys():
         plt.scatter(time, dic_time[time]['FluxMag0s_med'],s=s,alpha=0.95, marker='+', color='k') 
         plt.scatter(time, dic_time[time]['FluxMag0s_mean'],s=s,alpha=0.95, marker='+', color='r')
    plt.scatter([], [], marker='+', color='k',s=s,label='Median')
    plt.scatter([], [], marker='+', color='r',s=s,label='Mean')
    plt.legend()
    plotPath = outputPrefix+'ZpVsTime.png'
    plt.savefig(plotPath, format="png")
    # plt.show()


def plotAstrometry(dist, mag, snr, fit_params=None, brightSnr=100,
                   outputPrefix=""):
    """Plot angular distance between matched sources from different exposures.

    Creates a file containing the plot with a filename beginning with `outputPrefix`.

    Parameters
    ----------
    dist : list or numpy.array
        Separation from reference [mas]
    mag : list or numpy.array
        Mean magnitude of PSF flux
    snr : list or numpy.array
        Median SNR of PSF flux
    fit_params : list or numpy.array
        Fit parameters to display
    brightSnr : float, optional
        Minimum SNR for a star to be considered "bright".
    outputPrefix : str, optional
        Prefix to use for filename of plot file.  Will also be used in plot titles.
        E.g., outputPrefix='Cfht_output_r_' will result in a file named
           'Cfht_output_r_check_astrometry.png'
    """
    bright, = np.where(np.asarray(snr) > brightSnr)

    numMatched = len(dist)
    dist_median = np.median(dist)
    bright_dist_median = np.median(np.asarray(dist)[bright])

    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(18, 12))

    ax[0].hist(dist, bins=100, color=color['all'],
               histtype='stepfilled', orientation='horizontal')
    ax[0].hist(np.asarray(dist)[bright], bins=100, color=color['bright'],
               histtype='stepfilled', orientation='horizontal',
               label='SNR > %.0f' % brightSnr)

    ax[0].set_ylim([0., 500.])
    ax[0].set_ylabel("Distance [mas]")
    ax[0].set_title("Median : %.1f, %.1f mas" %
                    (bright_dist_median, dist_median),
                    x=0.55, y=0.88)
    plotOutlinedLinesHorizontal(ax[0], dist_median, bright_dist_median)

    ax[1].scatter(snr, dist, s=10, color=color['all'], label='All')
    ax[1].scatter(np.asarray(snr)[bright], np.asarray(dist)[bright], s=10,
                  color=color['bright'],
                  label='SNR > %.0f' % brightSnr)
    ax[1].set_xlabel("SNR")
    ax[1].set_xscale("log")
    ax[1].set_ylim([0., 500.])
    ax[1].set_title("# of matches : %d, %d" % (len(bright), numMatched))

    w, = np.where(dist < 200)
    plotAstromErrModelFit(snr[w], dist[w], fit_params=fit_params, ax=ax[1])

    ax[1].legend(loc='upper right')
    ax[1].axvline(brightSnr, color='red', linewidth=4, linestyle='dashed')
    plotOutlinedLinesHorizontal(ax[1], dist_median, bright_dist_median)

    plt.suptitle("Astrometry Check : %s" % outputPrefix.rstrip('_'), fontsize=30)
    plotPath = outputPrefix+"check_astrometry.png"
    plt.savefig(plotPath, format="png")
    plt.close(fig)


def plotExpFit(x, y, y_err, fit_params=None, deg=2, ax=None, verbose=False):
    """Plot an exponential quadratic fit to x, y, y_err.

    Parameters
    ----------
    fit_params : list or numpy.array
        Fit parameters to display
        If None, then will be fit.
    """

    if ax is None:
        ax = plt.figure()
        xlim = [1, 1e4]
    else:
        xlim = ax.get_xlim()

    if fit_params is None:
        fit_params = fitExp(x, y, y_err, deg=2)

    x_model = np.linspace(*xlim, num=100)
    fit_model = expModel(x_model, *fit_params)
    label = '%.4g exp(mag/%.4g) + %.4g' % \
            (fit_params[0], fit_params[2], fit_params[1])
    if verbose:
        print(fit_params)
        print(label)

    ax.plot(x_model, fit_model, color='red',
            label=label)


def plotAstromErrModelFit(snr, dist, fit_params=None,
                          color='red', ax=None, verbose=True):
    """Plot model of photometric error from LSST Overview paper
    http://arxiv.org/abs/0805.2366v4

    Astrometric Errors
    error = C * theta / SNR

    Parameters
    ----------
    snr : list or numpy.array
        S/N of photometric measurements
    dist : list or numpy.array
        Separation from reference [mas]
    """
    if ax is None:
        ax = plt.figure()
        xlim = [10, 30]
    else:
        xlim = ax.get_xlim()

    if fit_params is None:
        fit_params = fitAstromErrModel(snr, dist)

    x_model = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), num=100)
    fit_model_mas_err = astromErrModel(x_model, **fit_params)
    label = r'$C, \theta, \sigma_{\rm sys}$ =' + '\n' + \
            '{C:.2g}, {theta:.4g}, {sigmaSys:.4g} [mas]'.format(**fit_params)

    if verbose:
        print(fit_params)
        print(label)

    ax.plot(x_model, fit_model_mas_err,
            color=color, linewidth=2,
            label=label)
    # Set the x limits back to their original values.
    ax.set_xlim(xlim)


def plotPhotErrModelFit(mag, mmag_err, fit_params=None, color='red', ax=None, verbose=True):
    """Plot model of photometric error from LSST Overview paper (Eq. 4 & 5)

    Parameters
    ----------
    mag : list or numpy.array
        Magnitude
    mmag_err : list or numpy.array
        Magnitude uncertainty or variation in *mmag*.
    fit_params : list or numpy.array
        Fit parameters to display
    ax : matplotlib.Axis, optional
        The Axis object to plot to.
    verbose : bool, optional
        Produce extra output to STDOUT
    """

    if ax is None:
        ax = plt.figure()
        xlim = [10, 30]
    else:
        xlim = ax.get_xlim()

    if fit_params is None:
        fit_params = fitPhotErrModel(mag, mmag_err)

    x_model = np.linspace(*xlim, num=100)
    fit_model_mag_err = photErrModel(x_model, **fit_params)
    fit_model_mmag_err = 1000*fit_model_mag_err
    labelFormatStr = r'$\sigma_{{\rm sys}} {{\rm [mmag]}}$, $\gamma$, $m_5 {{\rm [mag]}}$=' + '\n' + \
                     r'{sigmaSysMmag:6.4f}, {gamma:6.4f}, {m5:6.3f}'
    label = labelFormatStr.format(sigmaSysMmag=1000*fit_params['sigmaSys'],
                                  **fit_params)

    if verbose:
        print(fit_params)
        print(label)

    ax.plot(x_model, fit_model_mmag_err,
            color=color, linewidth=2,
            label=label)

    return fit_params


def plotMagerrFit(*args, **kwargs):
    plotExpFit(*args, **kwargs)


def plotPhotometry(mag, snr, mmagerr, mmagrms, brightSnr=100,
                   fit_params=None,
                   filterName='Magnitude',
                   outputPrefix=""):
    """Plot photometric RMS for matched sources.

    Parameters
    ----------
    snr : list or numpy.array
        Median SNR of PSF flux
    mag : list or numpy.array
        Average Magnitude
    mmagerr : list or numpy.array
        Average Magnitude uncertainty [millimag]
    mmagrms ; list or numpy.array
        Magnitude RMS across visits [millimag]
    fit_params : list or numpy.array
        Fit parameters for photometry error model
    brightSnr : float, optional
        Minimum SNR for a star to be considered "bright".
    filterName : str, optional
        Name of the observed filter to use on axis labels.
    outputPrefix : str, optional
        Prefix to use for filename of plot file.  Will also be used in plot titles.
        E.g., outputPrefix='Cfht_output_r_' will result in a file named
           'Cfht_output_r_check_photometry.png'
    """

    bright, = np.where(np.asarray(snr) > brightSnr)

    numMatched = len(mag)
    mmagrms_median = np.median(mmagrms)
    bright_mmagrms_median = np.median(np.asarray(mmagrms)[bright])

    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(18, 16))
    ax[0][0].hist(mmagrms, bins=100, range=(0, 500), color=color['all'], label='All',
                  histtype='stepfilled', orientation='horizontal')
    ax[0][0].hist(np.asarray(mmagrms)[bright], bins=100, range=(0, 500), color=color['bright'],
                  label='SNR > %.0f' % brightSnr,
                  histtype='stepfilled', orientation='horizontal')
    ax[0][0].set_ylim([0, 500])
    ax[0][0].set_ylabel("RMS [mmag]")
    ax[0][0].set_title("Median : %.1f, %.1f mmag" %
                       (bright_mmagrms_median, mmagrms_median),
                       x=0.55, y=0.88)
    plotOutlinedLinesHorizontal(ax[0][0], mmagrms_median, bright_mmagrms_median)

    ax[0][1].scatter(mag, mmagrms, s=10, color=color['all'], label='All')
    ax[0][1].scatter(np.asarray(mag)[bright], np.asarray(mmagrms)[bright],
                     s=10, color=color['bright'],
                     label='SNR > %.0f' % brightSnr)

    ax[0][1].set_xlabel("%s [mag]" % filterName)
    ax[0][1].set_ylabel("RMS [mmag]")
    ax[0][1].set_xlim([17, 24])
    ax[0][1].set_ylim([0, 500])
    ax[0][1].set_title("# of matches : %d, %d" % (len(bright), numMatched))
    ax[0][1].legend(loc='upper left')
    plotOutlinedLinesHorizontal(ax[0][1], mmagrms_median, bright_mmagrms_median)

    ax[1][0].scatter(mmagrms, mmagerr, s=10, color=color['all'], label=None)
    ax[1][0].scatter(np.asarray(mmagrms)[bright], np.asarray(mmagerr)[bright],
                     s=10, color=color['bright'],
                     label=None)
    ax[1][0].set_xscale('log')
    ax[1][0].set_yscale('log')
    ax[1][0].plot([0, 1000], [0, 1000],
                  linestyle='--', color='black', linewidth=2)
    ax[1][0].set_xlabel("RMS [mmag]")
    ax[1][0].set_ylabel("Median Reported Magnitude Err [mmag]")
    brightSnrInMmag = 2.5*np.log10(1 + (1/brightSnr)) * 1000

    ax[1][0].axhline(brightSnrInMmag, color='red', linewidth=4, linestyle='dashed',
                     label=r'SNR > %.0f = $\sigma_{\rm mmag} < $ %0.1f' % (brightSnr, brightSnrInMmag))
    ax[1][0].set_xlim([1, 500])
    ax[1][0].set_ylim([1, 500])
    ax[1][0].legend(loc='upper center')

    ax[1][1].scatter(mag, mmagerr, color=color['all'], label=None)
    ax[1][1].set_yscale('log')
    ax[1][1].scatter(np.asarray(mag)[bright], np.asarray(mmagerr)[bright],
                     s=10, color=color['bright'],
                     label=None,
                     )
    ax[1][1].set_xlabel("%s [mag]" % filterName)
    ax[1][1].set_ylabel("Median Reported Magnitude Err [mmag]")
    ax[1][1].set_xlim([17, 24])
    ax[1][1].set_ylim([1, 500])
    ax[1][1].axhline(brightSnrInMmag, color='red', linewidth=4, linestyle='dashed',
                     label=r'$\sigma_{\rm mmag} < $ %0.1f' % (brightSnrInMmag))

    ax2 = ax[1][1].twinx()
    ax2.scatter(mag, snr,
                color=color['all'], facecolor='none',
                marker='.', label=None)
    ax2.scatter(np.asarray(mag)[bright], np.asarray(snr)[bright],
                color=color['bright'], facecolor='none',
                marker='.', label=None)
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel("SNR")
    ax2.axhline(brightSnr, color='red', linewidth=2, linestyle='dashed',
                label=r'SNR > %.0f' % (brightSnr))

    w, = np.where(mmagerr < 200)
    plotPhotErrModelFit(mag[w], mmagerr[w], fit_params=fit_params, ax=ax[1][1])
    ax[1][1].legend(loc='upper left')

    plt.suptitle("Photometry Check : %s" % outputPrefix.rstrip('_'), fontsize=30)
    plotPath = outputPrefix+"check_photometry.png"
    plt.savefig(plotPath, format="png")
    plt.close(fig)


def plotPA1(pa1, outputPrefix=""):
    """Plot the results of calculating the LSST SRC requirement PA1.

    Creates a file containing the plot with a filename beginning with `outputPrefix`.

    Parameters
    ----------
    pa1 : pipeBase.Struct
        Must contain:
        rms, iqr, magMean, magDiffs
        rmsUnits, iqrUnits, magDiffsUnits
    outputPrefix : str, optional
        Prefix to use for filename of plot file.  Will also be used in plot titles.
        E.g., outputPrefix='Cfht_output_r_' will result in a file named
           'Cfht_output_r_AM1_D_5_arcmin_17.0-21.5.png'
        for an AMx.name=='AM1' and AMx.magRange==[17, 21.5]
    """
    diffRange = (-100, +100)

    fig = plt.figure(figsize=(18, 12))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.scatter(pa1.magMean, pa1.magDiffs, s=10, color=color['bright'], linewidth=0)
    ax1.axhline(+pa1.rms, color=color['rms'], linewidth=3)
    ax1.axhline(-pa1.rms, color=color['rms'], linewidth=3)
    ax1.axhline(+pa1.iqr, color=color['iqr'], linewidth=3)
    ax1.axhline(-pa1.iqr, color=color['iqr'], linewidth=3)

    ax2 = fig.add_subplot(1, 2, 2, sharey=ax1)
    ax2.hist(pa1.magDiffs, bins=25, range=diffRange,
             orientation='horizontal', histtype='stepfilled',
             normed=True, color=color['bright'])
    ax2.set_xlabel("relative # / bin")

    yv = np.linspace(diffRange[0], diffRange[1], 100)
    ax2.plot(scipy.stats.norm.pdf(yv, scale=pa1.rms), yv,
             marker='', linestyle='-', linewidth=3, color=color['rms'],
             label="PA1(RMS) = %4.2f %s" % (pa1.rms, pa1.rmsUnits))
    ax2.plot(scipy.stats.norm.pdf(yv, scale=pa1.iqr), yv,
             marker='', linestyle='-', linewidth=3, color=color['iqr'],
             label="PA1(IQR) = %4.2f %s" % (pa1.iqr, pa1.iqrUnits))
    ax2.set_ylim(*diffRange)
    ax2.legend()
#    ax1.set_ylabel(u"12-pixel aperture magnitude diff (mmag)")
#    ax1.set_xlabel(u"12-pixel aperture magnitude")
    ax1.set_xlabel("psf magnitude")
    ax1.set_ylabel("psf magnitude diff (%s)" % pa1.magDiffsUnits)
    for label in ax2.get_yticklabels():
        label.set_visible(False)

    plt.suptitle("PA1: %s" % outputPrefix.rstrip('_'))
    plotPath = "%s%s" % (outputPrefix, "PA1.png")
    plt.savefig(plotPath, format="png")
    plt.close(fig)


def plotAM1(*args, **kwargs):
    return plotAMx(*args, x=1, **kwargs)


def plotAM2(*args, **kwargs):
    return plotAMx(*args, x=2, **kwargs)


def plotAM3(*args, **kwargs):
    return plotAMx(*args, x=3, **kwargs)


def plotAMx(AMx, outputPrefix=""):
    """Plot a histogram of the RMS in relative distance between pairs of stars.

    Creates a file containing the plot with a filename beginning with `outputPrefix`.

    Parameters
    ----------
    AMx : pipeBase.Struct
        Must contain:
        AMx, rmsDistMas, fractionOver, annulus, magRange, x, level,
        AMx_spec, AFx_spec, ADx_spec
    outputPrefix : str, optional
        Prefix to use for filename of plot file.  Will also be used in plot titles.
        E.g., outputPrefix='Cfht_output_r_' will result in a file named
           'Cfht_output_r_AM1_D_5_arcmin_17.0-21.5.png'
        for an AMx.name=='AM1' and AMx.magRange==[17, 21.5]
    """

    percentOver = 100*AMx.fractionOver

    AMxAsDict = AMx.getDict()
    AMxAsDict['AMxADx'] = AMxAsDict['AMx_spec']+AMxAsDict['ADx_spec']
    AMxAsDict['percentOver'] = percentOver

    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.hist(AMx.rmsDistMas, bins=25, range=(0.0, 100.0),
             histtype='stepfilled',
             label='D: %.1f-%.1f %s\nMag Bin: %.1f-%.1f' %
                   (AMx.annulus[0], AMx.annulus[1], AMx.annulusUnits, AMx.magRange[0], AMx.magRange[1]))
    ax1.axvline(AMx.AMx, 0, 1, linewidth=2, color='black',
                label='median RMS of relative\nseparation: %.2f %s' % (AMx.AMx, AMx.amxUnits))
    ax1.axvline(AMx.AMx_spec, 0, 1, linewidth=2, color='red',
                label='%s: %.0f %s' % (AMx.name, AMx.AMx_spec, AMx.amxUnits))
    formatStr = 'AM{x:d}+AD{x:d}: {AMxADx:.0f} {amxUnits:s}\n' + \
                'AF{x:d}: {AFx_spec:.2f}{afxUnits:s} > AM{x:d}+AD{x:d} = ' + \
                '{percentOver:.2f}%'
    ax1.axvline(AMx.AMx_spec+AMx.ADx_spec, 0, 1, linewidth=2, color='green',
                label=formatStr.format(**AMxAsDict))

    ax1.set_title('The %d stars separated by D = [%.2f, %.2f] %s' %
                  (len(AMx.rmsDistMas), AMx.annulus[0], AMx.annulus[1], AMx.annulusUnits))
    ax1.set_xlim(0.0, 100.0)
    ax1.set_xlabel('rms Relative Separation (%s)' % AMx.rmsUnits)
    ax1.set_ylabel('# pairs / bin')

    ax1.legend(loc='upper right', fontsize=16)

    figName = outputPrefix+'%s_D_%d_%s_%.1f-%.1f.png' % \
        (AMx.name, int(sum(AMx.annulus)/2), AMx.DUnits.upper(), AMx.magRange[0], AMx.magRange[1])

    plt.savefig(figName, dpi=300)
    plt.close(fig)
