import numpy
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable

data = []
for i in range(1,34):
    data.append(loadmat('/home/andrewstier/Downloads/sfi_gains_subak/volatility_shock_run_kl_model_seed_%d.mat'%i))
params = numpy.vstack(data[0]['params'][0]).astype(float)
params[:,1] = params[:,1]/10
### plot the power law with distance parameter over the grid
### insets blowing up one square with high dependence and one with lower dependence include variability there
### panel 2 is the plots of how many farms fail overall

### code for plotting the boundary/area and failure relationship
distance_pdfs = []
boundary_pdfs = []
distance_dependence = []
distance = numpy.array(range(21))
biggest_spike_failure = []
biggest_spike_failure_gtz = []
for tf in range(1,11):
    pdfs = []
    dp = []
    bsf = []
    bpdfs = []
    for sg in range(1,11):
        idx = numpy.argwhere((params[:, 0] == tf) & (params[:, 1] == sg / 10)).flatten()[0]
        distance_failure_pdfs = data[0]['failure_pdfs_distance_all'][0][idx]
        distance_pdf = numpy.nansum(distance_failure_pdfs, 0) / max(numpy.nansum(distance_failure_pdfs),1)
        bsf.append(numpy.diff(data[0]['dts_all'][0][idx]).max())
        bpdf = data[0]['failure_pdfs_psize_all'][0][idx]
        bpdfs.append(numpy.nansum(bpdf,0)/numpy.nansum(bpdf))
        power = numpy.nan
        try:
            power = curve_fit(lambda x, a, b: b * (x + 1) ** a, distance, distance_pdf, [1, 1])[0][0]
        except:
            pass
        dp.append(power)
        pdfs.append(distance_pdf)
    distance_pdfs.append(pdfs)
    boundary_pdfs.append(bpdfs)
    distance_dependence.append(dp)
    biggest_spike_failure.append(bsf)
    biggest_spike_failure_gtz.append((numpy.array(bsf)>0)*1.)
distance_dependence = numpy.vstack(distance_dependence)
for i in range(1,len(data)):
    dd = []
    for tf in range(1, 11):
        pdfs = []
        dp = []
        bsf = []
        bpdfs = []
        for sg in range(1, 11):
            idx = numpy.argwhere((params[:, 0] == tf) & (params[:, 1] == sg / 10)).flatten()[0]
            distance_failure_pdfs = data[i]['failure_pdfs_distance_all'][0][idx]
            distance_pdf = numpy.nansum(distance_failure_pdfs, 0) / max(
                numpy.nansum(distance_failure_pdfs), 1)
            pdfs.append(distance_pdf)
            bsf.append(numpy.diff(data[i]['dts_all'][0][idx]).max())
            bpdf = data[i]['failure_pdfs_psize_all'][0][idx]
            bpdfs.append(numpy.nansum(bpdf, 0) / numpy.nansum(bpdf))
            power = numpy.nan
            try:
                power = curve_fit(lambda x, a, b: b * (x + 1) ** a, distance, distance_pdf, [1, 1])[0][
                    0]
            except:
                pass
            dp.append(power)
        distance_pdfs[tf-1] = numpy.vstack(distance_pdfs[tf-1]) + numpy.vstack(pdfs)
        biggest_spike_failure[tf - 1] = numpy.hstack(biggest_spike_failure[tf - 1]) + numpy.hstack(bsf)
        biggest_spike_failure_gtz[tf - 1] = numpy.hstack(biggest_spike_failure_gtz[tf - 1]) + (numpy.hstack(bsf)>0)*1.
        boundary_pdfs[tf-1] = numpy.vstack(boundary_pdfs[tf - 1]) + numpy.vstack(bpdfs)
        dd.append(dp)
    dd = numpy.vstack(dd)
    distance_dependence = numpy.dstack((distance_dependence,dd))

plt.clf()
fig = plt.figure(figsize=(9,8))
ax = plt.subplot(111)
axin1 = ax.inset_axes([0, .2, 0.45, .75])
dd = numpy.nanmean(distance_dependence,2)
dd[dd==1] = numpy.nan
im = axin1.imshow(dd,cmap='RdBu',vmin=numpy.nanmin(dd),vmax=-numpy.nanmin(dd))
# add color bar below chart
divider = make_axes_locatable(axin1)
cax = divider.new_vertical(size='5%', pad=0.6, pack_start = True)
fig.add_axes(cax)
c=fig.colorbar(im, cax = cax, orientation = 'horizontal')
c.set_label(r'$\beta$')
axin1.set_ylabel('time to failure')
axin1.set_yticks(range(10),range(1,11))
axin1.set_xlabel('shock (fraction of water stress)')
axin1.set_xticks(range(10),numpy.array(range(1,11))/10)
axin1.text(.1,9,r'Pr(failure,$r_{periphery}$)'+'\n         =\n'+r'    $b*r_{periphery}^\beta$')
ax.set_axis_off()
#############################
axin2 = ax.inset_axes([0, .89, 0.15, .125])
axin3 = ax.inset_axes([.22, .89, 0.15, .125])
axin4 = ax.inset_axes([.43, .89, 0.15, .125])
##############################
idx = numpy.argwhere((params[:, 0] == 1) & (params[:, 1] == 10 / 10)).flatten()[0]
dy = [numpy.nansum(data[i]['failure_pdfs_distance_all'][0][idx],0) for i in range(len(data))]
dy = numpy.vstack([x/max(numpy.nansum(x),1) for x in dy]).mean(0)
axin2.plot(distance,dy)
axin2.set_ylim(0,1)
axin2.set_xlim(0,10)
axin2.set_ylabel('Pr(failure)')
axin2.set_xlabel('distance to periphery')
################################
idx3 = numpy.argwhere((params[:, 0] == 3) & (params[:, 1] == 8 / 10)).flatten()[0]
dy = [numpy.nansum(data[i]['failure_pdfs_distance_all'][0][idx3],0) for i in range(len(data))]
dy = numpy.vstack([x/numpy.nansum(x) for x in dy]).mean(0)
axin3.plot(distance,dy)
axin3.set_ylim(0,1)
axin3.set_xlim(0,10)
axin3.set_xlabel('distance to periphery')
################################
idx4 = numpy.argwhere((params[:, 0] == 5) & (params[:, 1] == 6 / 10)).flatten()[0]
dy = [numpy.nansum(data[i]['failure_pdfs_distance_all'][0][idx4],0) for i in range(len(data))]
dy = numpy.vstack([x/max(numpy.nansum(x),1) for x in dy]).mean(0)
axin4.plot(distance,dy)
axin4.set_ylim(0,1)
axin4.set_xlim(0,10)
axin4.set_xlabel('distance to periphery')
axin4.set_yticks([])
axin3.set_yticks([])
##################################
axin5 = ax.inset_axes([.5, .2, 0.45, .75])
bsf = numpy.vstack(biggest_spike_failure)/numpy.vstack(biggest_spike_failure_gtz)
bsf[bsf==0]= numpy.nan
im = axin5.imshow(bsf)
# add color bar below chart
divider = make_axes_locatable(axin5)
cax = divider.new_vertical(size='5%', pad=0.6, pack_start = True)
fig.add_axes(cax)
c=fig.colorbar(im, cax = cax, orientation = 'horizontal')
c.set_label('largest spike in failures')
axin5.set_yticks(range(10),range(1,11))
axin5.set_xlabel('shock (fraction of water stress)')
axin5.set_xticks(range(10),numpy.array(range(1,11))/10)
axin5.set_yticks([])
#########################################
axin6 = ax.inset_axes([.38, .1, 0.15, .125])
axin7 = ax.inset_axes([.61, .1, 0.15, .125])
axin8 = ax.inset_axes([.83, .1, 0.15, .125])
##############################
idx = numpy.argwhere((params[:, 0] == 1) & (params[:, 1] == 10 / 10)).flatten()[0]
dy = data[0]['dts_all'][0][idx].flatten()
axin6.plot(range(1,len(dy)+1),dy)
#####################################
idx = numpy.argwhere((params[:, 0] == 3) & (params[:, 1] == 8 / 10)).flatten()[0]
dy = data[0]['dts_all'][0][idx].flatten()
axin7.plot(range(1,len(dy)+1),dy)
######################################
idx = numpy.argwhere((params[:, 0] == 5) & (params[:, 1] == 6 / 10)).flatten()[0]
dy = data[0]['dts_all'][0][idx].flatten()
axin8.plot(range(1,len(dy)+1),dy)
axin8.set_yticks([0,1e-3],['0',r'$10^{-3}$'])
axin6.set_ylabel('failed farms')
axin6.set_xlabel('model time')
axin7.set_xlabel('model time')
axin8.set_xlabel('model time')
plt.tight_layout()
plt.savefig('/home/andrewstier/Downloads/sfi_gains_subak/figures/h1_all.png',dpi=350)


plt.clf()
fig = plt.figure(figsize=(10,10))
c=1
for i in range(10):
    for j in range(10):
        ax = plt.subplot(10,10,c)
        p = boundary_pdfs[i][j,:][1:]
        if numpy.nansum(p)>0:
            plt.plot(numpy.array(range(len(p)))/len(p),p,alpha=.75)
            plt.yticks([])
        else:
            ax.set_axis_off()
        c+=1
plt.tight_layout()
plt.savefig('/home/andrewstier/Downloads/sfi_gains_subak/figures/shape_boundary.png',dpi=350)
