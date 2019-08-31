import numpy as np
import matplotlib
matplotlib.use("Qt4agg")
import matplotlib.pyplot as plt
import pdb
import kplr
import fnmatch
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.stats import LombScargle
from scipy.signal import savgol_filter as savgol
from astropy import units as u
import lightkurve as lkk
import matplotlib.lines as mlines
from echelle import plot_echelle


# Want to debug? Use pdb.set_trace() 

#Action Items:
        ## identify delta v and make echelle DONE
        ## make echelle better - 
        ## plot 1D echelle diagram
        ## plot on Jen's graphs 
        ## tell if it is l=0 etc.DOING 
        ## make the echelle diagram DONE
        ## identify peaks (might need to decrease smoothening)
        ## do the asteroseismic age for the star using the small seperations (not vmax)
        ## look through all again and identify good ones
        
# subroutine to perform rough sigma clipping
def sigclip(x,y,subs,sig):
    keep = np.zeros_like(x)
    start=0
    end=subs
    nsubs=int((len(x)/subs)+1)
    for i in range(0,nsubs):        
        me=np.mean(y[start:end])
        sd=np.std(y[start:end])
        good=np.where((y[start:end] > me-sig*sd) & (y[start:end] < me+sig*sd))[0]
        keep[start:end][good]=1
        start=start+subs
        end=end+subs
    return keep


# main program starts here
if __name__ == '__main__':
    
#load gaia data
    gaiad = np.loadtxt('./observed_sc_targets_Q0-Q17_prob_gaia.txt',skiprows=1)
    
#load rotation period data : Table 1 of http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJS/211/24
    keepkics = np.loadtxt('./keepgyr.txt',skiprows=1) #kicids with rot periods # only one col so no need to sep into cols
    
#load prev measured asteroseismic detections. : Table 1 of http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJS/210/1
    ignorekics = np.loadtxt('./ignoreprev.txt',skiprows=1)#kicids with prev periods astero detections
    
    kicidlistuf = gaiad[:,0]
    problistuf = gaiad[:,5]
    tlengthlistuf = gaiad[:,4] #unit is days
    tefflistuf = gaiad[:,1]
    radlistuf = gaiad[:,2]
    kepmaglistuf = gaiad[:,3]

    kicidlist = np.array([]) #creating empty array to append it to
    problist = np.array([])
    tlengthlist = np.array([])
    tefflist = np.array([])
    radlist = np.array([])
    kepmaglist = np.array([])

    for i in range(0,len(kicidlistuf)):
        if (problistuf[i] <= 0.89): # selecting *s with prob of asteroseismic detection > 0.89
            continue
        if (tlengthlistuf[i] != 30):# selecting *s with 30 day long exposures
            continue
        if (radlistuf[i] > 10): # selecting *s smaller than 10*solar radius
            continue
        
        igkic = np.where(ignorekics==kicidlistuf[i])[0]
        kekic = np.where(keepkics==kicidlistuf[i])[0]
        
        if (len(igkic)==1): # selecting *s that do not have previously detected asteroseismic measurements
            continue
        if (len(kekic)==0): # selecting *s that have a rotation period 
            continue
        
        kicid = kicidlistuf[i]
        
        ## putting in st coz got cut off
        if (int(kicid)!=11029516): #the vvg one!
            continue
        ##
        
        prob = problistuf[i]
        tlength = tlengthlistuf[i]
        kicrow = list(kicidlistuf).index(kicid) # get row value where particular kicid valid
        
        teff = gaiad[kicrow,1]
        rad = gaiad[kicrow,2]
        kepmag = gaiad[kicrow,3]

        kicidlist = np.append(kicidlist,kicid) # append to create a list
        tlengthlist = np.append(tlengthlist,tlength)
        problist = np.append(problist,prob)
        tefflist = np.append(tefflist,teff)
        radlist = np.append(radlist,rad)
        kepmaglist = np.append(kepmaglist,kepmag)

    for kicid in kicidlist:
        kicrow = list(kicidlistuf).index(kicid) # get row value where particular kicid valid
        teff = gaiad[kicrow,1]
        rad = gaiad[kicrow,2]
        kepmag = gaiad[kicrow,3]

        client = kplr.API()
        star = client.star(kicid)
        lcs = star.get_light_curves(short_cadence=True)
        if (len(lcs)==0):
            continue
        time, flux, qual = [], [], []
        for lc in lcs:
            # if data is long-cadence, skip this entry
            if (fnmatch.fnmatch(str(lc),'*LC*')):
                continue
            with lc.open() as f:
                # The lightcurve data are in the first FITS HDU.
                print(lc)
                hdu_data = f[1].data
                time.append(hdu_data["time"])
                flux.append(hdu_data["pdcsap_flux"])
                qual.append(hdu_data["sap_quality"])
        # convert to simple arrays
        if (time==[]):
            continue
        time=time[0]
        flux=flux[0]
        qual=qual[0]
    
        # only keep data with good quality flags
        good=np.where(qual == 0)[0]
        time=time[good]
        flux=flux[good]

        # plot the light curve
        plt.ion()
        # plt.clf()
        # plt.subplot(3,1,1)
        # plt.plot(time,flux)
        # plt.xlabel('Time (Days)')
        # plt.ylabel('Flux (counts)')

        # sigma-clip outliers from the light curve and overplot it
        res=sigclip(time,flux,50,3)
        good = np.where(res == 1)[0]
        time=time[good]
        flux=flux[good]
        #plt.plot(time,flux)
        
####################################################
     # plotting
         # shift for a certain star flux 
        m_ms=1.2 #m/msun = 1.2
        vmaxs = 3100 #microHz
        tsun = 5778 #K
        vmax = (m_ms*((teff/tsun)**3.5)*vmaxs)/((rad**2)*((teff/tsun)**4))
        width= (1/(0.1*vmax*(10**(-6))))/86400
        
        print("width=",width,"days \n"
              "teff=",teff,"K \n"
              "radius=",rad,"R_sun", "\n"
              "nu_max=",vmax,"microHz")

        #plt.title('some text') 
        #plt.title(str(numax))

            # next, run a filter through the data to remove long-periodic (low frequency) variations
            # let's pick a 5 day width for the filter

            #5  ### tie to vmax cox T = 1/f but width is in days
        ###### WIDTH has to be an odd number ####
        #pdb.set_trace()
            #86400 s in a day
            ## if estimated vmax is very low coould be correct
            ## title with estimated vmax
            
            
        boxsize= width/(1./60./24.)
        box_kernel = Box1DKernel(boxsize)
        intboxsize = int(boxsize)
        if intboxsize%2 == 0:
            intboxsize = intboxsize+1
        smoothed_flux = savgol(flux,intboxsize,1,mode='mirror')
            # overplot this smoothed version, and then divide the light curve through it
        #plt.plot(time,smoothed_flux)
        #plt.title("kicid is "+ str(int(kicid)) +" & numax is "+ str(int(vmax))+" $\mu$Hz")

        flux=flux/(smoothed_flux)

            # plot the filtered light curve
        # plt.subplot(3,1,2)
        # plt.plot(time,flux)
        # plt.xlabel('Time (Days)')
        # plt.ylabel('Relative flux')

            # now let's calculate the fourier transform. the nyquist frequency is:
        nyq=1./(1./60./24.)
            ## expected to be 40 microHz set an appropriate filter
            ## red giant stars long cadence not short cadenc 

            # FT magic
        freq, amp = LombScargle(time,flux).autopower(method='fast',samples_per_peak=50,maximum_frequency=nyq)#10

            # unit conversions
        freq = 1000.*freq/86.4
        bin = freq[1]-freq[0]
        amp = 2.*amp*np.var(flux*1e6)/(np.sum(amp)*bin)
        gauss_kernel = Gaussian1DKernel(26)
        pssm = convolve(amp, gauss_kernel)

        #     # plot the power spectrum log scale
        # plt.subplot(3,2,5)
        # plt.loglog(freq,amp)
        # plt.loglog(freq,pssm)
        # plt.axvline(x=vmax,linewidth=2, color='r')
        # plt.xlabel('Frequency ($\mu$Hz)')
        # plt.ylabel('Power Density')
        # plt.xlim([100,8000]) ## plot all 10-8000
        # plt.tight_layout()

        # #data artifact 300-400 microhz
        #     # plot the power spectrum regular ### low freq no point in plotting. 
        # plt.subplot(3,2,6)
        # plt.plot(freq,amp) 
        # plt.plot(freq,pssm)
        # plt.axvline(x=vmax,linewidth=2, color='r')
        # plt.xlabel('Frequency ($\mu$Hz)')
        # plt.ylabel('Power Density')
        # plt.xlim([vmax-0.5*vmax,vmax+0.5*vmax])
        
        # plt.tight_layout()
        
        # save the output as png
        #plt.savefig(str(input("Is it good(g) or bad(b)?"))+'_'+str(int(kicid))+'.png',dpi=200)

#_______________________________________________________________________________
#                             Echelle
#_______________________________________________________________________________
   #Generate lightkurve 
        lk = lkk.LightCurve(time=time,flux=flux)
       
   #Do Fourier transform
        pg = lk.to_periodogram(method='lombscargle',normalization='psd',minimum_frequency=1000,maximum_frequency=3000)
       
   #Get numax, seis  
        snr = pg.flatten()
        seis = snr.to_seismology()
        numax = seis.estimate_numax()
       
   #Change deltanu
        adelenu= np.float64(0.9)#*u.uHz #
        #ax = plt.figure()#.add_suplot(111)
        #for adelenu in np.arange(0.8,1.1,0.02):
        #adelenu=np.float64(input("add or subtract how much from dnu?"))*u.uHz
        seis.deltanu=seis.estimate_deltanu()+(adelenu*u.uHz)
            #seis.plot_echelle(deltanu=seis.deltanu,numax=numax,smooth_filter_width=3.,scale='log',cmap='viridis')
        #plt.savefig('final_Echelle'+str(seis.deltanu)+'_'+'.png',dpi=200)
            #plt.clear()
            #plt.pause(0.01)
            
        ##l=1
        l1freq =np.array([2316.36,2220.47,2122.40,2024.33,1925.68,1829.04])
        #l1freqog = np.array([2418.48,2316.36,2220.47,2122.40,2024.33,1925.68,1829.04,1733.90,1550.97])
        l1freqmod = l1freq%(seis.deltanu/u.uHz)
        #print(l1freqmod)
        
        ##l=2
        l2freq = np.array([2264.08,2166.27,1971.82,1873.42])
        #l2freqog = np.array([2264.08,2166.27,2076.13,1971.82,1873.42])
        l2freqmod = l2freq%(seis.deltanu/u.uHz)
        #print(l2freqmod)

        ##l=0
        l0freq = np.array([2269.92,2173.08,2076.13,1977.90,1881.65])
        #l0freqog = np.array([2269.92,2173.08,2084.47,1977.90,1881.65])
        l0freqmod = l0freq%(seis.deltanu/u.uHz)
        #print(l0freqmod)

### Plot modes on power spectrum
        # plt.figure()
        # plt.plot(time,flux)
        # plt.plot(freq,pssm,'black')
        # #plt.xlim([1750,2400])
        # plt.ylim([-0.5,21])
        # nu_max = plt.axvline(x=vmax,linewidth=5, color='yellow',alpha=0.5,label='nu_max') #vmaxline -- note vmax more accurate than numax

        # for l1freaks in l1freq:
        #     l1 =plt.axvline(x=l1freaks,linewidth=1.1,color='fuchsia',ls='--',label='l1')
        # for l2freaks in l2freq:
        #     l2 = plt.axvline(x=l2freaks,linewidth=1.1,color='springgreen',ls='--',label='l2')
        # for l0freaks in l0freq:
        #     l0 = plt.axvline(x=l0freaks,linewidth=1.1,color='dodgerblue',ls='--',label='l0')
        # plt.legend(handles=(nu_max,l1,l2,l0))

        # plt.xlabel('Frequency ($\mu$Hz)') 
        # plt.ylabel('Power Density')
        # plt.xlim([vmax-0.5*vmax,vmax+0.5*vmax])

#### Plot Echelle
        # seis.plot_echelle(deltanu=seis.deltanu,numax=numax,smooth_filter_width=3.,scale='log',cmap='gray')
        # plt.plot(l1freqmod,l1freq,'-ok',color='fuchsia')
        # plt.plot(l2freqmod,l2freq,'-ok',color='springgreen')
        # plt.plot(l0freqmod,l0freq,'-ok',color='dodgerblue')
        # l1ech = mlines.Line2D([], [], color='fuchsia', marker='.',
        #                   markersize=10, label='l = 1')
        # l2ech = mlines.Line2D([], [], color='springgreen', marker='.',
        #                   markersize=10, label='l = 2')
        # l0ech = mlines.Line2D([], [], color='dodgerblue', marker='.',
        #                   markersize=10, label='l = 0')
        # plt.legend(handles=(l1ech,l2ech,l0ech))
        # plt.show()

        ###
        #print(freq) 
        # freqmodech = freq%(seis.deltanu/u.uHz)
        # plt.figure()
        # plt.scatter(freqmodech,freq)
        # plt.yscale('log')#,cmap='viridis')
        #plt.imshow
        

### danhey ech
        #plt.figure()
        ## plot_echelle(freq,pssm,97.73,fmin=1700,fmax=2500)
        from echelle import echelle
        X,Y,Z = echelle(freq,pssm,97.73) ## X = freqmod, Y= pssm, Z=freq        
        
        #plt.ylim([1600,2400])
        
        # plt.plot(l1freqmod,l1freq,'-ok',color='fuchsia')
        # plt.plot(l2freqmod,l2freq,'-ok',color='springgreen')
        # plt.plot(l0freqmod,l0freq,'-ok',color='dodgerblue') 
        # l1ech = mlines.Line2D([], [], color='fuchsia', marker='.',
        #                   markersize=2, label='l = 1')
        # l2ech = mlines.Line2D([], [], color='springgreen', marker='.',
        #                   markersize=2, label='l = 2')
        # l0ech = mlines.Line2D([], [], color='dodgerblue', marker='.',
        #                   markersize=2, label='l = 0')
        # plt.legend(handles=(l1ech,l2ech,l0ech))
        # plt.show()
        #pdb.set_trace()
        plt.plot(X, np.sum(Z, axis=0), 'k', linewidth=0.7)
        plt.axvline(x=np.mean(l1freqmod),linewidth=8, color='fuchsia',alpha=0.4,label='l1')
        plt.axvline(x=np.mean(l2freqmod),linewidth=8, color='springgreen',alpha=0.4,label='l2')
        plt.axvline(x=np.mean(l0freqmod),linewidth=8, color='dodgerblue',alpha=0.4,label='l0')
        l1ech = mlines.Line2D([], [], color='fuchsia', marker='s',markersize=8,alpha=0.4, label='l = 1')
        l2ech = mlines.Line2D([], [], color='springgreen', marker='s',markersize=8,alpha=0.4, label='l = 2')
        l0ech = mlines.Line2D([], [], color='dodgerblue', marker='s',markersize=8,alpha=0.4, label='l = 0')
        plt.legend(handles=(l1ech,l2ech,l0ech))
        plt.xlabel('Frequency mod 97.73')
        plt.ylabel('Relative Power') #Ampl in danhey code

# # fig, axes = plt.subplots(2,1, figsize=[10,10])

# # dnu = 7

# # ax = axes[0]
# # echelle.plot_echelle(pg.frequency.value, np.sqrt(pg.power.value), dnu, fmin=10, fmax=70, ax=ax, cmap='Blues')
# # ax.set_xlabel('')
# # ax.set_xticks([])

# # ax = axes[1]
# # X, Y, Z = echelle.echelle(pg.frequency.value, pg.power.value, dnu)
# # ax.plot(X, np.sum(Z, axis=0), 'k', linewidth=0.7)
# # ax.set_xlabel('Frequency mod 5')
# # ax.set_ylabel('Amplitude')
# # ax.set_xlim(X[0], X[-1])
# # ax.set_ylim(0, None)

# # plt.subplots_adjust(hspace=0.)

        # redfreq = freq.flatten()
        #redfreqmod = redfreq%97.73
        #freqmod=freq%97.73
        # l1y1 =l1freqmod/10.0
        # l2y2 =l2freqmod/10.0
        # l0y0 =l0freqmod/10.0
        # #plt.plot(freqmod,amp)
        # plt.plot(freqmod,pssm-amp)
        # plt.plot(l1freqmod,l1y1,'-ok',color='fuchsia')
        # plt.plot(l2freqmod,l2y2,'-ok',color='springgreen')
        # plt.plot(l0freqmod,l0y0,'-ok',color='dodgerblue')
        # # plt.legend(handles=(l1ech,l2ech,l0ech))
        # plt.show()

        

        #imagine.flatten()
        #print(freqmod)
        #print((np.sum(freq,axis=0)).shape())
        #print((np.sum(freq),shape()))
        #plt.plot(freqmod,np.sum(freq,axis=0),'k')
        #from echelle import echelle
       # X, Y, Z = echelle(freq, pssm, 97.73)
        #Yf = Y.flatten()
        #imagine = plt.contourf(X, Y, Z)
        #imagine.flatten()
             #    pssm)

## calc mass
        msun = 1.9891 * (10**30) #kg
        deltanusun = 135 * (10**-6) #Hz
        numaxsun = vmaxs * (10**-6) #Hz
        teffsun = tsun #K
        rsun = 6.957 * (10**8) #m

        teffstar = teff #K
        radstar = rad * rsun #m
        deltanustar = (seis.deltanu/u.uHz) * (10**-6) #Hz
        #print(deltanustar)
        numaxstar = vmax * (10**-6)
        #print(numaxstar)
        
        mstarnum = msun * (numaxstar**3) * (deltanusun**4) * (teffstar**1.5)
        mstardenom = (numaxsun**3) * (deltanustar**4) * (teffsun**1.5)
        mstar = mstarnum/mstardenom ##kg

        mstarsol = mstar/msun #in Msun

        print("mstar = ",mstar, " kg")
        print("mstar = ",mstarsol, " m_sun")
        #)#*tsun#*(vmaxs*(10**-6))



        ###### binrubbish
   #seis.plot_echelle(deltanu=seis.deltanu,numax=numax,smooth_filter_width=3.,scale='log',cmap='viridis')
        #from pprint import pprint 
        #pprint(seis)

        #integer rounding #nrows, cols
        
        input(':')

###### meeting notes ##################
        ###
        # hack the large sep to get a better diag
        #over sample the power spectrum
        # make freq grid much smaller

        # find a way to use your power specturm but lk's echelle diag
        # look through others

        # possibly collaps eon one axis
        ## do that if not obvious on echelle diagream

        # rough estimate on large sep and put
        # next week - remind about metallicity bc all plots are assuming some metallicity.

        ## good goal : dv02 spacing and other plot white pper dv02 and and put jen's
        ## error bar micqillan paper . 
        

        ### tried to show intitiative
        ### getting upto speed with codings
        ## read up about the physics behind it 
        
        ###good taking notes
        ###good showing motivations
        ###good ask qs about the physics behind it


        ## look at power spectrum which peaks are significant DONE
        ## identify peaks DONE
        ## decrease smoothening DONE
        ## echelle flatten get dv02 WAITING ON LIGHTKURVE PPL 

        # coding questions
        ## good 

        ## measure radius
        ## dan hey freq range <-
        ## plotting issue shifted echeel

        # lorentizaian profile
        ## gaussian has a width of sigman, and a mean.
        ## lorientiazian also have a freq and a width and height
        ## lorentizian -- osc exciticed by convection, dame out - lorentxian diecribes it
        ## gamma prop 1/mode  lifetime
        ## short mode lifetimes, narraw peak
        
