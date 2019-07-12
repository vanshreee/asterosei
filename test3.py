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

# subroutine to perform rough sigma clipping
def sigclip(x,y,subs,sig):
    keep = np.zeros_like(x)
    #pdb.set_trace()
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
#
    gaiad = np.loadtxt('./observed_sc_targets_Q0-Q17_prob_gaia.txt',skiprows=1)
    
#load rotation period data : Table 1 of http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJS/211/24
#
    gyrkics = np.loadtxt('./keepthese.txt',skiprows=1)#kicids with rot periods
    
#load prev measured asteroseismic detections. : Spec table of http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJS/210/1
#
    prevkics = np.loadtxt('./ignorethese.txt',skiprows=1)#kicids with rot periods
    
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
        if (problistuf[i] > 0.89):
            if (tlengthlistuf[i] == 30):
                kicid = kicidlistuf[i]
                prob = problistuf[i]
                tlength = tlengthlistuf[i]
                kicrow = list(kicidlistuf).index(kicid) # get row value where particular kicid valid
                #print(kicrow)
                teff = gaiad[kicrow,1]
                rad = gaiad[kicrow,2]
                kepmag = gaiad[kicrow,3]
                

                kicidlist = np.append(kicidlist,kicid)
                tlengthlist = np.append(tlengthlist,tlength)
                problist = np.append(problist,prob)
                tefflist = np.append(tefflist,teff)
                radlist = np.append(radlist,rad)
                kepmaglist = np.append(kepmaglist,kepmag)
                
    for kicid in kicidlist:
        ## highlight promising ones
        ## radius filter 10 solar mass
        ## plt.savefig('test.png') 
        #if kicid != 2013883:
        #    continue
        kicrow = list(kicidlistuf).index(kicid) # get row value where particular kicid valid
        #print(kicrow)
        teff = gaiad[kicrow,1]
        rad = gaiad[kicrow,2]
        kepmag = gaiad[kicrow,3]

        client = kplr.API()
        star = client.star(kicid)
        lcs = star.get_light_curves(short_cadence=True)

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
## print lcs and check if it has anything. if it doesnt, skip that star.
        #pdb.set_trace()
        ### error st ## mayb not finding sdata for star.
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
##### limit on star radius.. < 10 solar rad
##### shift for a certain star flux 
        m_ms=1.2 #m/msun = 1.2
        vmaxs = 3100 #microHz
        tsun = 5778 #K
        vmax = (m_ms*((teff/tsun)**3.5)*vmaxs)/((rad**2)*((teff/tsun)**4))
        width= (1/(0.1*vmax*(10**(-6))))/86400
        print("width=",width,"days \n"
              "teff=",teff,"K \n"
              "radius=",rad,"R_sun", "\n"
              "nu_max=",vmax,"microHz")
# plt.title 
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
            
            #
        boxsize= width/(1./60./24.)
        box_kernel = Box1DKernel(boxsize)
        intboxsize = int(boxsize)
        if intboxsize%2 == 0:
            intboxsize = intboxsize+1
        smoothed_flux = savgol(flux,intboxsize,1,mode='mirror')
            # overplot this smoothed version, and then divide the light curve through it
        plt.plot(time,smoothed_flux)

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
        plt.xlim([vmax-0.5*vmax,vmax+0.5*vmax])
        plt.tight_layout()


        # plot the power spectrum regular
        plt.subplot(3,2,6)
        plt.plot(freq,amp) ####### change x range to vmax + 0.5 vmax
        plt.plot(freq,pssm)
        plt.axvline(x=vmax,linewidth=2, color='r')
        plt.xlabel('Frequency ($\mu$Hz)')
        plt.ylabel('Power Density')
        plt.xlim([vmax-0.5*vmax,vmax+0.5*vmax])
        plt.tight_layout()

        # save the output as png
        plt.savefig('fig.png',dpi=200)

        input(':')

        ### sap flux : aperture added flux
        ###pdc flux has corrections by kepler team


        ## rotating star having asteroseismic detection really valuable. Sending 2 catalogs with rotation periods## skip the ones with asteroseismic detection and prioriitize those with a rotation period ## make sample manageable. From those, look at high probablity stars.
        # keep the code as you have, download the catalogs.. if my star in this catalog, skip it etc.. bunch of if statement s to implement that. 
        
