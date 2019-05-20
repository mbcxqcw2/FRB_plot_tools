import numpy as np
import matplotlib.gridspec as gs
from matplotlib import pyplot as plt
import os
from sigpyproc.Readers import FilReader as fr
from astropy.convolution import convolve,Box1DKernel
import sys

#input parameters

a=sys.argv

n_inputs=9

if len(a)!=n_inputs+1:
    print "Error: {0} arguments required. {1} provided.".format(n_inputs,len(a))

#search algorithm used to create candidate file
#options: presto, destroy, dedisperse_all, seek, astroaccelerate, heimdall

#currently having issues with astroAccelerate as it doesn't recognise liam's filterbank file as an 8-bit file

searchtype = str(a[1])
if searchtype not in ['presto','seek','destroy','dedisperse_all','astroaccelerate','heimdall']:
    print 'Warning: invalid search algorithm. Must be presto, seek, destroy, dedisperse_all, astroaccelerate or heimdall.'

print 'Candidates files generated using {0} will be plotted.'.format(searchtype)

#filterbank file used in searching

filfile=str(a[2])

print 'The filterbank file used to generate the candidates was: {0}'.format(filfile)

#dm range of candidates to plot (pc/cc)

loDM = np.float(a[3])
hiDM = np.float(a[4])

print 'Candidates between {0} and {1} pc/cc will be plotted.'.format(loDM, hiDM)

#time range of candidates to plot (seconds)

loTime = np.float(a[5])
hiTime = np.float(a[6])

print 'Candidates between {0} and {1} seconds will be plotted.'.format(loTime,hiTime)

#snr threshold to plot

snrThresh = np.float(a[7])

print 'Candidates above the S/N threshold: {0} will be plotted'.format(snrThresh)

#candidate file location

#if search type was presto, input the directory containing the single_pulse_search.py output
#'.singlepulse' candidate files, as it splits them up by DM.

#if search type was destroy, input the directory containing the output '.pls'
#candidate files. Dstroy must be run multiple times on multiple DM .tim files.

#if search type was seek, rename output '.pls' candidate file extension '.s' to
#avoid confusion with destroy. Seek must be run multiple times on multiple DM .tim
#files.

#if search type was dedisperse_all, ensure the output cand file was given
#a '.dd' file extension. Input its directory location.

candfile_loc = a[8]

print 'Candidate files within directory: {0} will be plotted.'.format(candfile_loc)

#boxcar width to convolve data with when plotting timeseries/dynamic spectra

boxcar = int(a[9])

print 'When plotting dynamic spectra, data will be convolved with a boxcar of width {0} bins'.format(boxcar)

#########################################################################################

class Click():
    def __init__(self, ax, func, button=1):
        self.ax=ax
        self.func=func
        self.button=button
        self.press=False
        self.move = False
        self.c1=self.ax.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2=self.ax.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.c3=self.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmove)

    def onclick(self,event):
        if event.inaxes == self.ax:
            if event.button == self.button:
                self.func(event, self.ax,boxcar)
    def onpress(self,event):
        self.press=True
    def onmove(self,event):
        if self.press:
            self.move=True
    def onrelease(self,event):
        if self.press and not self.move:
            self.onclick(event)
        self.press=False; self.move=False


def func(event, ax,boxcar):
    """
    Function to plot dynamic spectra when clicked

    """
    #extract time and DM info from plot click
    time = event.xdata
    DM = event.ydata
    print 'Clicked DM, time were: {0} pc/cc, {1} s'.format(DM,time)
    print 'Filterbank file details...'
    print 'Number of samples: {0}, Sampling time: {1}, top frequency channel: {2}, bottom frequency channel {3}'.format(nsamples,
                                                                                                                   samptime,
                                                                                                                   topchan,
                                                                                                                   botchan)


    #find extent of DM sweep between band top and bottom for pulse
    #according to DM
    dt = (1./2.41e-16)*(((botchan*1e6)**(-2)-((topchan*1e6)**(-2))))*DM
    print '\nCalculated DM sweep across the band: {0}s'.format(dt)

    #find nearest sample, sample time in fb file to the FRB
    nearest_sample = (np.abs(timesamples-time)).argmin()

    #find dm sweep in timesamples
    dt_samps = int(np.ceil(np.abs(dt)/samptime))

    #get the sample range to read
    minsamp = nearest_sample-(dt_samps)
    nsamps = 2*int(round(dt_samps/np.float(boxcar))*np.int(boxcar)) #make sure it is a factor of boxcar so data can be downsampled

    #read filterbank
    block=fil.readBlock(minsamp,nsamps)
    block_dedisp = block.dedisperse(DM)

    #create dispersed downsampled image
    for i in range(block.shape[0]):
        block[i,:] = convolve(block[i,:],Box1DKernel(int(boxcar)))


    #downsample by convolving with boxcar
    dedisp_timeseries = block_dedisp.sum(axis=0)
    dedisp_timeseries_downsamp = convolve(dedisp_timeseries,Box1DKernel(int(boxcar)))
    dedisp_timeseries_downsamp = dedisp_timeseries_downsamp[boxcar-1:-(boxcar-1+1)]#remove the first downsamp-1 and last downsamp-1 bins as convolution pads with zeros.
                                                                                   #the additional +1 added to the end sample is as python is zero-indexed. -1 removes no bins,... etc
    for i in range(block_dedisp.shape[0]):
        block_dedisp[i,:] = convolve(block_dedisp[i,:],Box1DKernel(int(boxcar)))
    
    #plot downsampled dedispersed pulse
    grid=gs.GridSpec(5,8)
    fig=plt.figure(figsize=(5,8))

    ax1=fig.add_subplot(grid[1:3,0:6])
    #ax1.imshow(block_dedisp,origin='lower',aspect='auto',vmax=np.mean(block_dedisp[:,0:((block_dedisp.shape[1]/2)-int(boxcar))]),vmin=4*np.std(block_dedisp[:,0:((block_dedisp.shape[1]/2)-int(boxcar))]),cmap='hot',extent=[timesamples[minsamp],timesamples[minsamp+nsamps],topchan,botchan])
    ax1.imshow(block_dedisp,origin='lower',aspect='auto',cmap='hot_r',extent=[timesamples[minsamp],timesamples[minsamp+nsamps],topchan,botchan])
    ax1.set_xlabel('Seconds since {0} {1}'.format(fil.header.obs_date,fil.header.obs_time),fontsize=8)
    ax1.set_ylabel('Frequency (MHz)',fontsize=8)
    ax1.set_title('Dedispersed Dynamic Spectrum')

    #plot downsampled timeseries
    ax2=fig.add_subplot(grid[0:1,0:6])
    ax2.set_title('Dedispersed Timeseries')
    ax2.set_xlim([timesamples[minsamp+(boxcar-1)],timesamples[minsamp+nsamps-(boxcar-1)]])
    ax2.plot(timesamples[minsamp+(boxcar-1):minsamp+nsamps-(boxcar-1+1)],dedisp_timeseries_downsamp)

    #plot dispersed pulse
    ax3=fig.add_subplot(grid[3:5,0:6])
    ax3.imshow(block,origin='lower',aspect='auto',cmap='gray_r',extent=[timesamples[minsamp],timesamples[minsamp+nsamps],topchan,botchan])
    ax3.set_xlabel('Seconds since {0} {1}'.format(fil.header.obs_date,fil.header.obs_time),fontsize=8)
    ax3.set_ylabel('Frequency (MHz)',fontsize=8)
    ax3.set_title('Dispersed Dynamic Spectrum')

    #plot dispersed bandpass
    ax4=fig.add_subplot(grid[3:5,6:8])
    ax4.set_title('Bandpass')
    ax4.set_xticklabels([])
    ax4.set_yticklabels([])
    ax4.plot(block.sum(axis=1),np.arange(block.shape[0]))

    #plot dedispersed bandpass
    ax5=fig.add_subplot(grid[1:3,6:8])
    ax5.set_title('Bandpass')
    ax5.set_xticklabels([])
    ax5.set_yticklabels([])
    ax5.plot(block_dedisp.sum(axis=1),np.arange(block.shape[0]))
    plt.tight_layout()
    plt.show()




#define file extension based on search algorithm
if searchtype=='presto':
    candfile_ext = '.singlepulse'
elif searchtype=='destroy':
    candfile_ext = '.pls'
elif searchtype=='dedisperse_all':
    candfile_ext = '.dd'
elif searchtype=='seek':
    candfile_ext = '.s'
elif searchtype=='astroaccelerate':
    candfile_ext = '.aa'
elif searchtype=='heimdall':
    candfile_ext = '.h'


#extract filterbank header
fil=fr(filfile)
filhead=fil.header
nsamples=filhead.nsamples
samptime=filhead.tsamp
topchan=filhead.ftop
botchan=filhead.fbottom
#times od samples in filterbank
timesamples=np.arange(nsamples)*samptime

#extract all candidate files
files = os.listdir(candfile_loc)
candfiles = [i for i in files if i[-len(candfile_ext):]==candfile_ext]

print '... {0} candidate files with file extension: {1} found.'.format(len(candfiles),candfile_ext)

#case: presto
if searchtype=='presto':
    
    #initialise presto candidate array (contains five columns)
    cands=np.empty((1,5))
    for i in range(len(candfiles)):
        #load first file
        check=np.loadtxt(candfile_loc+candfiles[i])

        if check.shape!=(0,): #if cand file is not empty
            if check.shape==(5,): #if only one candidate, reshape so appending can happen
                check=check.reshape(1,5)
            #append candidates
            cands=np.concatenate((cands,check),axis=0)
            
    #reassign candidates to arrays
    times = cands[:,2]
    dms = cands[:,0]
    snrs = cands[:,1]
    sample = cands[:,3]
    downsamp = cands[:,4]

#case: destroy
if searchtype=='destroy':
    
    #initialise destroy candidate array (contains four columns)
    cands=np.empty((1,4))
    for i in range(len(candfiles)):
        #load first file
        check=np.loadtxt(candfile_loc+candfiles[i])
        if check.shape!=(0,): #if cand file is not empty
            if check.shape==(4,): #if only one candidate, reshape to allow appending
                check=check.reshape(1,4)
            cands=np.concatenate((cands,check),axis=0) #append candidates
    #reassign candidates to arrays
    dms = cands[:,0]
    downsamp = cands[:,1] #known as nsmoothings in destroy files. check this is boxcar width in bins.
    sample = cands[:,2]
    snrs = cands[:,3]
    times = sample*filhead.tsamp

#case: dedisperse_all
if searchtype=='dedisperse_all':
    
    #initialise dedisperse_all cand. array (contains 11 columns)
    cands=np.empty((1,11))
    for i in range(len(candfiles)):
        check=np.loadtxt(candfile_loc+candfiles[i],dtype=str)
        if check.shape!=(0,): #if cand file is not empty
            if check.shape==(11,):#if only one candidate, reshape to allow appending
                check=check.reshape(1,11)
            cands=np.concatenate((cands,check),axis=0) #append candidates
    #reassign candidates to arrays #need to check with evan that downsamp and sample are the columns I think...
    snrs = np.array(cands[:,1],dtype='float') #1 or 2???
    dms = np.array(cands[:,7],dtype='float')
    downsamp = np.array(cands[:,3],dtype='float')
    sample = np.array(cands[:,10],dtype='float') #4 or 5 or 10???
    times = sample*filhead.tsamp
    
#case: seek
if searchtype=='seek':
    
    #initialise seek cand. array (contains 5 columns)
    cands=np.empty((1,5))
    for i in range(len(candfiles)):
        check=np.loadtxt(candfile_loc+candfiles[i],skiprows=1)
        if check.shape!=(0,):#if candfile is not empty
            if check.shape==(5,):#if only one candidate, reshape to allow appending
                check=check.reshape(1,5)
            cands=np.concatenate((cands,check),axis=0)
    #reassign cands to arrays
    dms = cands[:,0]
    downsamp = cands[:,1]
    sample=cands[:,2]
    snrs=cands[:,3]
    times=sample*filhead.tsamp

#case: AstroAccelerate
if searchtype=='astroaccelerate':
    
    #initialise aa cand. array (contains 4 columns)
    cands=np.empty((1,4))
    for i in range(len(candfiles)):
        #open binary .dat file
        binfile=open(candfile_loc+candfiles[i])
        #load candidates
        check=np.fromfile(binfile,dtype='f4')
        if check.shape!=(0,):#if candfile is not empty
            if check.shape==(4,):#if only one candidate, reshape to allow appending
                check=check.reshape(1,4)
            else:#else reshape to correct shape
                check=check.reshape((check.shape[0]/4,4))
            cands=np.concatenate((cands,check),axis=0)
        #close binary file        
        binfile.close()
    #reassign cands to arrays
    dms      = cands[:,0]
    downsamp = cands[:,3]
    snrs     = cands[:,2]
    times    = cands[:,1]
    sample   = np.zeros_like(times) #astro_accelerate doesn't provide sample number

#case: heimdall
if searchtype=='heimdall':

    #initialise heimdall cand. array (contains 9 columns)
    cands=np.empty((1,9))
    for i in range(len(candfiles)):
        check=np.loadtxt(candfile_loc+candfiles[i],skiprows=0)
        if check.shape!=(0,):#if candidate file is not empty
            if check.shape==(9,0):#if only one candidate, reshape to allow appending
                check=check.reshape(1,9)
            cands=np.concatenate((cands,check),axis=0)
    #reassign cands to arrays
    dms      = cands[:,5]
    downsamp = cands[:,3]
    sample   = cands[:,1]#heimdall returns highest frequency channel sample number
    snrs     = cands[:,0]
    times    = cands[:,2]
    
print 'Total number of candidates found by {0}: {1}'.format(searchtype,cands.shape[0])
    
#mask all cands outside of DM range
DMmask=np.where((dms>=loDM) & (dms <=hiDM))

print 'Candidates masked by DM range: {0}'.format(cands.shape[0]-DMmask[0].shape[0])

times=times[DMmask]
dms=dms[DMmask]
snrs=snrs[DMmask]
sample=sample[DMmask]
downsamp=downsamp[DMmask]

#mask all cands outside of time range
timemask=np.where((times>=loTime) & (times<=hiTime))

print 'Candidates masked by time range: {0}'.format(cands.shape[0]-timemask[0].shape[0])

times=times[timemask]
dms=dms[timemask]
snrs=snrs[timemask]
sample=sample[timemask]
downsamp=downsamp[timemask]

#mask all cands below snr threshold
snrmask=np.where((snrs>snrThresh))

print 'Candidates masked by S/N threshold: {0}'.format(cands.shape[0]-snrmask[0].shape[0])

times=times[snrmask]
dms=dms[snrmask]
snrs=snrs[snrmask]
sample=sample[snrmask]
downsamp=downsamp[snrmask]

print 'Plotting...'


#plots

grid=gs.GridSpec(9,9,wspace=12,hspace=12)
fig=plt.figure(figsize=(9,9))

#Time-DM plot
ax1=fig.add_subplot(grid[3:9,0:9])
ax1.scatter(times,dms,s=np.array(snrs,dtype=int)**2,facecolors='none',edgecolors='k')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('DM (pc/cc)')

#S/N-DM plot
ax2=fig.add_subplot(grid[0:3,6:9])
plt.scatter(dms,snrs,facecolors='none',edgecolors='k')
plt.xlabel('DM (pc/cc)')
plt.ylabel('S/N')

#S/N histogram
ax3=fig.add_subplot(grid[0:3,3:6])
ax3.hist(snrs,facecolor='none',edgecolor='k')
ax3.set_yscale('log')
ax3.set_xlabel('S/N')
ax3.set_ylabel('Number of pulses')

#DM histogram
ax4=fig.add_subplot(grid[0:3,0:3])
ax4.hist(dms,facecolor='none',edgecolor='k')
ax4.set_xlabel('DM (pc/cc)')
ax4.set_ylabel('Number of pulses')


ax3.set_title('Single pulse results for "{0}"\
\nAlgorithm used: {1}\
\n\n Source: {2}   RA (J2000): {3}   N samples: {4}\
\n Telescope ID: {5}   Dec (J2000): {6}   Sampling time: {7} ms\
\n Machine ID: {8}    tstart: {9}   Freq_{{ctr}}: {10}\n'.format(filhead.filename.split('/')[-1],
                        searchtype,
                        filhead.source_name,
                        filhead.ra,
                        filhead.nsamples,
                        filhead.telescope_id,
                        filhead.dec,
                        filhead.tsamp,
                        filhead.machine_id,
                        filhead.tstart,
                        filhead.fcenter),fontsize=8)
fig.tight_layout()
#plt.savefig('single_pulse_cands_{0}.png'.format(searchtype))

click = Click(ax1, func, button=1)
plt.show()
