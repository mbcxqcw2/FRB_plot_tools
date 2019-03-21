import numpy as np
import matplotlib.gridspec as gs
from matplotlib import pyplot as plt
import os
from sigpyproc.Readers import FilReader as fr
import sys

#input parameters

a=sys.argv

n_inputs=8

if len(a)!=n_inputs+1:
    print "Error: {0} arguments required. {1} provided.".format(n_inputs,len(a))

#search algorithm used to create candidate file
#options: presto, destroy, dedisperse_all, seek

#currently having issues with astroAccelerate as it doesn't recognise liam's filterbank file as an 8-bit file

searchtype = str(a[1])
if searchtype not in ['presto','seek','destroy','dedisperse_all']:
    print 'Warning: invalid search algorithm. Must be presto, seek, destroy or dedisperse_all.'

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

#########################################################################################


#define file extension based on search algorithm
if searchtype=='presto':
    candfile_ext = '.singlepulse'
elif searchtype=='destroy':
    candfile_ext = '.pls'
elif searchtype=='dedisperse_all':
    candfile_ext = '.dd'
elif searchtype=='seek':
    candfile_ext = '.s'


#extract filterbank header
fil=fr(filfile)
filhead=fil.header

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
                check=check.reshape(1,4)
            cands=np.concatenate((cands,check),axis=0)
    #reassign cands to arrays
    dms = cands[:,0]
    downsamp = cands[:,1]
    sample=cands[:,2]
    snrs=cands[:,3]
    times=sample*filhead.tsamp
    
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
plt.show()
