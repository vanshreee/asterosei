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
#import echelle #don't need this unless using danhey's 
from astropy import units as u
import lightkurve as lkk

# Want to debug? Use pdb.set_trace() 

#Action Items:
        ## identify delta v and make echelle DONE
        ## plot 1D echelle diagram
        ## plot on Jen's graphs 
        ## tell if it is l=0 etc.DOING 
        ## make the echelle diagram DONE
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
        plt.clf()
        plt.subplot(3,1,1)
        plt.plot(time,flux)
        plt.xlabel('Time (Days)')
        plt.ylabel('Flux (counts)')

        # sigma-clip outliers from the light curve and overplot it
        res=sigclip(time,flux,50,3)
        good = np.where(res == 1)[0]
        time=time[good]
        flux=flux[good]
        plt.plot(time,flux)
        
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
        plt.plot(time,smoothed_flux)
        plt.title("kicid is "+ str(int(kicid)) +" & numax is "+ str(int(vmax))+" $\mu$Hz")

        flux=flux/(smoothed_flux)

            # plot the filtered light curve
        plt.subplot(3,1,2)
        plt.plot(time,flux)
        plt.xlabel('Time (Days)')
        plt.ylabel('Relative flux')

            # now let's calculate the fourier transform. the nyquist frequency is:
        nyq=1./(1./60./24.)
            ## expected to be 40 microHz set an appropriate filter
            ## red giant stars long cadence not short cadenc 

            # FT magic
        freq, amp = LombScargle(time,flux).autopower(method='fast',samples_per_peak=10,maximum_frequency=nyq)

            # unit conversions
        freq = 1000.*freq/86.4
        bin = freq[1]-freq[0]
        amp = 2.*amp*np.var(flux*1e6)/(np.sum(amp)*bin)
        gauss_kernel = Gaussian1DKernel(26)
        pssm = convolve(amp, gauss_kernel)

            # plot the power spectrum log scale
        plt.subplot(3,2,5)
        plt.loglog(freq,amp)
        plt.loglog(freq,pssm)
        plt.axvline(x=vmax,linewidth=2, color='r')
        plt.xlabel('Frequency ($\mu$Hz)')
        plt.ylabel('Power Density')
        plt.xlim([100,8000]) ## plot all 10-8000
        plt.tight_layout()

        #data artifact 300-400 microhz
            # plot the power spectrum regular ### low freq no point in plotting. 
        plt.subplot(3,2,6)
        plt.plot(freq,amp) 
        plt.plot(freq,pssm)
        plt.axvline(x=vmax,linewidth=2, color='r')
        plt.xlabel('Frequency ($\mu$Hz)')
        plt.ylabel('Power Density')
        plt.xlim([vmax-0.5*vmax,vmax+0.5*vmax]) # x range goes to vmax +- 0.5 vmax
        plt.tight_layout()

## also was just plotting power specturen seperaterly here: 
        #plt.figure()
        #plt.plot(time,flux)
        # plt.plot(freq,pssm)     
        # plt.axvline(x=vmax,linewidth=2, color='r')
        # plt.xlabel('Frequency ($\mu$Hz)') 
        # plt.ylabel('Power Density')
        # plt.xlim([vmax-0.5*vmax,vmax+0.5*vmax])
        # x range goes to vmax +- 0.5 vmaxm

       import warnings
        warnings.filterwarnings('ignore')
        import lightkurve as lkk

    #Getting Data from MAST
        datalist = lkk.search_lightcurvefile('KIC11029516',  cadence='short')
        data = datalist[:].download_all()
        lk = data.stitch

    #Plot Normalized flux over Time
        lk=lk().normalize().remove_outliers().remove_nans()
        #lk.plot();

    #Plot Fourier transform with smoothened one &numax
        pg = lk.to_periodogram(method='lombscargle',normalization='psd',minimum_frequency=1000,maximum_frequency=3000)
        #ax = pg.plot()
        #ax.axvline(pg.frequency_at_max_power.value,lw=5,ls='-.')

        #pg.smooth(method= 'boxkernel', filter_width=1.).plot(ax=ax,label='Smoothed',c='red',lw=2)
        #pg.plot(scale = 'log')
        #pg.show_properties()

    #Plot Signal-to-Noise and use seismology module
        snr = pg.flatten()
        #snr.plot()
        seis = snr.to_seismology()
        #seis
        #seis.periodogram.plot()
        
        numax = seis.estimate_numax()
        #deltanu = seis.estimate_deltanu()
        #print(np.dtype(deltanu))
        
    #Plot that measured deltanu
        #ax = seis.diagnose_deltanu()

        
        # save the output as png
        #plt.savefig(str(input("Is it good(g) or bad(b)?"))+'_'+str(int(kicid))+'.png',dpi=200)

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
