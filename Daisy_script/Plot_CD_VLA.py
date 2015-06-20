'''
Look for .fits in ../PriorImage/

Purpose:
--------
plots CD, and VLA,

Notes:
------

'''

from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from FgContIm import *

path = '../PriorImage/'
Plotpath = '../OverlayFig/'

label = dict.fromkeys(['VLA', 'mom0', 'cont'])
for k in label.iterkeys():
    files = glob.glob(path+'*'+k+'*.fits')
    label[k] = files

print label

#import sys; sys.exit()

########################################
# user define area
########################################
# 13:35:42.245, 30:04:06.854
VLA_ra = 203.927
VLA_dec = 30.06857
sigma_VLA = 7.1885e-6
vla_min = -3e-5
vla_max = 1.08e58-4

ra_center = VLA_ra
dec_center = VLA_dec

sigma_line = 0.5161
sigma_cont = 2.444e-4

sizep = 0.0035
sizePLine = 0.0040
sizeCont = 0.007

# # positions for crosses, centering on radio core
ra_cross, dec_cross = ra_center, dec_center

figCD = plt.figure(1, figsize=(12, 5))
fig = plt.figure(4, figsize=(12, 5))

row_a = 0.10
width = 0.35
x_gap = 0.05
x0 = 0.10
dy = 0.90


########################################
# intialize base figure
########################################
fvlaCD = aplpy.FITSFigure(label['VLA'][0], \
        figure=figCD, subplot=[x0,row_a,width,dy])
fvlaCD.show_grayscale() #stretch='log', vmin=vla_min, vmax=vla_max, vmid=vla_min-0.1)

flinCD = aplpy.FITSFigure(label['mom0'][0], \
        figure=figCD, subplot=[x0+width+2*x_gap, row_a, width, dy])
flinCD.show_colorscale(cmap=mpl.cm.jet)

figCont = aplpy.FITSFigure(label['cont'][0], \
        figure=fig, subplot=[x0,row_a,width,dy])
figCont.show_colorscale(cmap=mpl.cm.jet)
figContzoom = aplpy.FITSFigure(label['cont'][0], \
        figure=fig, subplot=[x0+width+2*x_gap, row_a, width, dy])
figContzoom.show_colorscale(cmap=mpl.cm.jet)


########################################
# Contours
########################################

fvlaCD.show_contour(label['VLA'][0], colors="blue", levels=sigma_contour_array(sigma_VLA), linewidths=2, layer='fg')
fvlaCD.show_contour(label['mom0'][0], colors='red', levels=sigma_contour_CARMA(sigma_line), linewidths=2, layer='fg_cont')
flinCD.show_contour(label['mom0'][0], colors="white", levels=sigma_contour_CARMA(sigma_line), linewidths=2, layer='mol')
figCont.show_contour(label['cont'][0], colors='red', levels=sigma_contour_CARMA(sigma_cont), linewidths=1.5, layer='fg_cont')
figContzoom.show_contour(label['cont'][0], colors='red', levels=sigma_contour_CARMA(sigma_cont), linewidths=1.5, layer='fg_cont')
########################################
# beam
########################################

setup_beam(fvlaCD)
setup_beam(flinCD)
setup_beam(figCont)
setup_beam(figContzoom)


########################################
# scale bar
########################################
lg_1arcsec = 1./3600
# lg_20kpc_fg = lg_1arcsec * 20./scale_radio
# lg_20kpc_bg = lg_1arcsec * 20./scale_SMG
# setup_scalebar(fvla, lg_20kpc_fg, str('20kpc'))
# setup_scalebar(fcont, lg_20kpc_fg, str('20kpc'))
# setup_scalebar(flin, lg_20kpc_bg, str('20kpc'))
# setup_scalebar(fSMA, lg_20kpc_bg, str('20kpc'))

########################################
# axes
########################################
standard_plot_setup(fvlaCD, ra_center, dec_center, sizep)
standard_plot_setup(flinCD, ra_center, dec_center, sizePLine)
standard_plot_setup(figCont, ra_center, dec_center, sizeCont)
standard_plot_setup(figContzoom, ra_center, dec_center, sizePLine)
# fcont.tick_labels.hide()
# fcont.axis_labels.hide()
flinCD.axis_labels.hide()
figContzoom.axis_labels.hide()


########################################
# markers
########################################
markers_cross(fvlaCD, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(flinCD, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(figCont, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(figContzoom, ra_cross, dec_cross, layer='marker_set_1')


########################################
# Labels
########################################
put_label(fvlaCD, 0.25, 0.95, 'VLA 6GHz', 'titleBand')
put_label(fvlaCD, 0.2625, 0.9, 'NA.v1.489', 'titleObj')
put_label(flinCD, 0.31, 0.95, 'CARMA CO(3-2)', 'titleBand', c='black')
put_label(flinCD, 0.31, 0.9, 'NA.v1.489', 'titleObj', c='black')
put_label(figCont, 0.45, 0.95, 'CARMA Continuum', 'titleBand', c='black')
put_label(figContzoom, 0.5, 0.95, 'CARMA Continuum Zoomed in', 'titleBand', c='black')


labsize = 'xx-large'
labc = 'white'
#put_label(fvla, 0.80, 0.925, '(a)', 'ref', c=labc, s=labsize)
#put_label(fcont, 0.80, 0.925, '(b)', 'ref', c=labc, s=labsize)
#put_label(flin, 0.80, 0.925, '(c)', 'ref', c=labc, s=labsize)
#put_label(fSMA, 0.80, 0.925, '(d)', 'ref', c=labc, s=labsize)
########################################
# Colorbar
########################################
# axisflin = fig.add_axes([0.92,0.19,0.02,0.68])
# normflin = mpl.colors.Normalize(vmin=min_line, vmax=max_line)
# cbflin = mpl.colorbar.ColorbarBase(axisflin, cmap=mpl.cm.jet, norm=normflin, orientation='vertical')
# cbflin.set_label('mJy')
# fig.canvas.draw()
# fig_line.canvas.draw()
# plt.show()

if __name__ == '__main__':
    """
    run script.py True
    sys.argv[1] determines whether to save all the plots or not
    """
    import sys, os
    if len(sys.argv) < 2:
        errmsg = "Invalid number of arguments: {0:d}\n  run script.py Save_True"
        raise IndexError(errmsg.format(len(sys.argv)))
    saveFig = True if sys.argv[1].lower() == 'true' else False
    if saveFig == True:
        os.system('rm -rf ' + 'CO32mom0.png' + ' ' + 'CO32mom0.eps')
#        figCD.savefig(Plotpath + CD + '.eps', dpi=600)
        print ("saving to {:s}").format(Plotpath + 'CO32mom0.png')
        figCD.savefig(Plotpath + 'CO32mom0.png', dpi=300)
        fig.savefig(Plotpath + 'Cont.png', dpi=300)
    else:
        plt.show()
