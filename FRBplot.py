from matplotlib import pyplot as plt
from sigpyproc.Readers import FilReader as fr
import numpy as np
from astropy.convolution import convolve,Box1DKernel
import matplotlib.gridspec as gs
import sys

#input parameters
a=sys.argv
n_inputs=4
if len(a)!=n_inputs+1:
    print "Insufficient arguments. {0} needed".format(n_inputs)

#filterbank file
filfile=str(a[1])

#expected FRB time (s)
time=float(a[2])

#expected frb DM (pc/cc)
DM=float(a[3])

#boxcar downsampling factor (ideally FRB width) in bins
downsamp=int(a[4])



print 'FRB details...'
print 'Expected FRB time: {0} s, Expected FRB DM: {1}, Expected FRB width/downsampling factor {2}'.format(time,DM,downsamp)

#read filterbank file header
fil=fr(filfile)
nsamples=fil.header.nsamples
samptime=fil.header.tsamp
topchan=fil.header.ftop
botchan=fil.header.fbottom
#times od samples in filterbank
timesamples=np.arange(nsamples)*samptime

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
nsamps = 2*int(round(dt_samps/np.float(downsamp))*np.int(downsamp)) #make sure it is a factor of downsamp so data can be downsampled

#read filterbank
block=fil.readBlock(minsamp,nsamps)
block_dedisp = block.dedisperse(DM)

#create dispersed downsampled image
for i in range(block.shape[0]):
    block[i,:] = convolve(block[i,:],Box1DKernel(int(downsamp)))


#downsample by convolving with boxcar
dedisp_timeseries = block_dedisp.sum(axis=0)
dedisp_timeseries_downsamp = convolve(dedisp_timeseries,Box1DKernel(int(downsamp)))
dedisp_timeseries_downsamp = dedisp_timeseries_downsamp[downsamp-1:-(downsamp-1)]#remove the first downsamp-1 and last downsamp-1 bins as convolution pads with zeros
for i in range(block_dedisp.shape[0]):
    block_dedisp[i,:] = convolve(block_dedisp[i,:],Box1DKernel(int(downsamp)))
    
#plot downsampled dedispersed pulse
grid=gs.GridSpec(5,8)
fig=plt.figure(figsize=(5,8))

ax1=fig.add_subplot(grid[1:3,0:6])
ax1.imshow(block_dedisp,origin='lower',aspect='auto',cmap='gray_r',extent=[timesamples[minsamp],timesamples[minsamp+nsamps],topchan,botchan])
ax1.set_xlabel('Seconds since {0} {1}'.format(fil.header.obs_date,fil.header.obs_time),fontsize=8)
ax1.set_ylabel('Frequency (MHz)',fontsize=8)
ax1.set_title('Dedispersed Dynamic Spectrum')

#plot downsampled timeseries
ax2=fig.add_subplot(grid[0:1,0:6])
ax2.set_title('Dedispersed Timeseries')
print minsamp,minsamp+nsamps
ax2.set_xlim([timesamples[minsamp+(downsamp-1)],timesamples[minsamp+nsamps-(downsamp-1)]])
ax2.plot(timesamples[minsamp+(downsamp-1):minsamp+nsamps-(downsamp-1)],dedisp_timeseries_downsamp)

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
