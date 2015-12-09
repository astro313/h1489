"""
Author: Daisy Leung

Last edited: Dec 09 2015

History:
Dec 09 2015
- changed path and Plotpath as I moved this directory to inside CARMA/

Jun 19 2015
- Plot HST, Keck NIR images with CO(3-2)
- fix alignment of HST

"""


from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from FgContIm import *

path = '/Users/admin/Research/hatlas_Cooray/PriorImage/'
Plotpath = '/Users/admin/Research/hatlas_Cooray/OverlayFig/'
# z = '2.685.'
z = '2.685_arc.bin2'
z2 = '2.688.'
zN = '2.688.NReg.'

sym = z        # change the sigma
#CD, C, D = 'CD_USB_HST'+sym, 'C_USB_HST'+sym, 'D_USB_HST'+sym
CD = 'CD_USB_'+sym

label = dict.fromkeys(['J133543', '-H', '-Ks', CD])
for k in label.iterkeys():
    files = glob.glob(path+'*'+k+'*.fits')
    label[k] = files

print label

figC = plt.figure(3, figsize=(12, 5))
figD = plt.figure(2, figsize=(12, 5))
figCD = plt.figure(1, figsize=(12, 5))

########################################
# user define area
########################################
# 13:35:42.245, 30:04:06.854
VLA_ra = 203.927
VLA_dec = 30.06857

ra_center = VLA_ra
dec_center = VLA_dec

# # # NReg:
# sigma_C = 0.655
# sigma_D = 0.76
# sigma_CD = 0.57

# 2.685
sigma_C = 0.65
sigma_D = 0.973
sigma_CD = 0.675

# 2.688
# sigma_C = 0.5
# sigma_D = 0.998
# sigma_CD = 0.615

sizep = 0.0035
sizePLine = 0.0040

# # positions for crosses, centering on radio core
ra_cross, dec_cross = ra_center, dec_center
row_a = 0.10
width = 0.35
x_gap = 0.05
x0 = 0.10
dy = 0.90


########################################
# intialize base figure
########################################
fvlaC = aplpy.FITSFigure(label['-Ks'][0], \
        figure=figC, subplot=[x0,row_a,width,dy])
fvlaC.show_grayscale() #stretch='log', vmin=vla_min, vmax=vla_max, vmid=vla_min-0.1)

flinC = aplpy.FITSFigure(label[CD][0], \
        figure=figC, subplot=[x0+width+2*x_gap, row_a, width, dy])
flinC.show_colorscale(cmap=mpl.cm.jet) #, vmin=min_line, vmax=max_line, vmid=min_line-0.5, stretch='log')

fvlaD = aplpy.FITSFigure(label['J133543'][0], \
        figure=figD, subplot=[x0,row_a,width,dy], north=True)
fvlaD.show_grayscale() #stretch='log', vmin=vla_min, vmax=vla_max, vmid=vla_min-0.1)

flinD = aplpy.FITSFigure(label[CD][0], \
        figure=figD, subplot=[x0+width+2*x_gap, row_a, width, dy])
flinD.show_colorscale(cmap=mpl.cm.jet)


fvlaCD = aplpy.FITSFigure(label['-H'][0], \
        figure=figCD, subplot=[x0,row_a,width,dy])
fvlaCD.show_grayscale() #stretch='log', vmin=vla_min, vmax=vla_max, vmid=vla_min-0.1)

flinCD = aplpy.FITSFigure(label[CD][0], \
        figure=figCD, subplot=[x0+width+2*x_gap, row_a, width, dy])
flinCD.show_colorscale(cmap=mpl.cm.jet)


########################################
# Contours
########################################
#fvlaC.show_contour(label['-Ks.'][0], colors="lime", levels=sigma_contour_array(sigma_VLA), linewidths=2, layer='fg')
fvlaC.show_contour(label[CD][0], colors='red', levels=sigma_contour_CARMA(sigma_CD), linewidths=2, layer='fg_cont')
flinC.show_contour(label[CD][0], colors="white", levels=sigma_contour_CARMA(sigma_CD), linewidths=2, layer='mol')

#fvlaD.show_contour(label['J133543'][0], colors="lime", levels=sigma_contour_array(sigma_VLA), linewidths=2, layer='fg')
fvlaD.show_contour(label[CD][0], colors='red', levels=sigma_contour_CARMA(sigma_CD), linewidths=2, layer='fg_cont')
flinD.show_contour(label[CD][0], colors="white", levels=sigma_contour_CARMA(sigma_CD), linewidths=2, layer='mol')

#fvlaCD.show_contour(label['-H.'][0], colors="lime", levels=sigma_contour_array(sigma_VLA), linewidths=2, layer='fg')
fvlaCD.show_contour(label[CD][0], colors='red', levels=sigma_contour_CARMA(sigma_CD), linewidths=2, layer='fg_cont')
flinCD.show_contour(label[CD][0], colors="white", levels=sigma_contour_CARMA(sigma_CD), linewidths=2, layer='mol')


########################################
# beam
########################################
#setup_beam(fvlaC)
#setup_beam(fvlaD)
#setup_beam(fvlaCD)
setup_beam(flinC)
setup_beam(flinD)
setup_beam(flinCD)

########################################
# Compass, bug in APLpy
########################################
# flinD.compass = aplpy.overlays.Compass()
# import pdb;pdb.set_trace()
# flinD.compass.show_compass()

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
standard_plot_setup(fvlaC, ra_center, dec_center, sizep)
standard_plot_setup(flinC, ra_center, dec_center, sizePLine, latex=False)
standard_plot_setup(fvlaD, ra_center, dec_center, sizep)
standard_plot_setup(flinD, ra_center, dec_center, sizePLine, latex=False)
standard_plot_setup(fvlaCD, ra_center, dec_center, sizep)
standard_plot_setup(flinCD, ra_center, dec_center, sizePLine)
# fcont.tick_labels.hide()
# fcont.axis_labels.hide()
# flin.tick_labels.hide()
# flin.axis_labels.hide()


########################################
# markers
########################################
markers_cross(fvlaC, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(flinC, ra_cross, dec_cross, layer='marker_set_1')

markers_cross(fvlaD, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(flinD, ra_cross, dec_cross, layer='marker_set_1')

markers_cross(fvlaCD, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(flinCD, ra_cross, dec_cross, layer='marker_set_1')

########################################
# Labels
########################################
if '_' in sym[:-1]: symf = sym.replace('_', ' ')
put_label(fvlaC, 0.25, 0.95, 'Keck Ks', 'titleBand')
put_label(fvlaC, 0.2625, 0.9, 'NA.v1.489', 'titleObj')
put_label(flinC, 0.5, 0.95, 'CARMA CO(3-2) CD data'+symf[:-1], 'titleBand', c='black')
put_label(flinC, 0.31, 0.9, 'NA.v1.489', 'titleObj', c='black')

put_label(fvlaD, 0.25, 0.95, 'J133543', 'titleBand')
put_label(fvlaD, 0.2625, 0.9, 'NA.v1.489', 'titleObj')
put_label(flinD, 0.5, 0.95, 'CARMA CO(3-2) CD data'+symf[:-1], 'titleBand', c='black')
put_label(flinD, 0.31, 0.9, 'NA.v1.489', 'titleObj', c='black')

put_label(fvlaCD, 0.25, 0.95, 'Keck H', 'titleBand')
put_label(fvlaCD, 0.2625, 0.9, 'NA.v1.489', 'titleObj')
put_label(flinCD, 0.5, 0.95, 'CARMA CO(3-2) CD data'+symf[:-1], 'titleBand', c='black')
put_label(flinCD, 0.31, 0.9, 'NA.v1.489', 'titleObj', c='black')

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
#        os.system('rm -rf ' + C[:-1] + '.png' + ' ' + C[:-1] + '.eps')
#        os.system('rm -rf ' + D[:-1] + '.png' + ' ' + D[:-1] + '.eps')
#        os.system('rm -rf ' + CD[:-1] + '.png' + ' ' + CD[:-1] + '.eps')
#        figC.savefig(Plotpath + C[:-1] + '.eps', dpi=600)
        figC.savefig(Plotpath + CD[:-1] + 'KeckKs.png', dpi=300)
#        figD.savefig(Plotpath + D[:-1] + '.eps', dpi=600)
        figD.savefig(Plotpath + CD[:-1] + 'HST.png', dpi=300)
#        figCD.savefig(Plotpath + CD[:-1] + '.eps', dpi=600)
        figCD.savefig(Plotpath + CD[:-1] + 'KeckH.png', dpi=300)
    else:
        plt.show()

