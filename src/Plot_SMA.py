'''
Last Modified: Dec 09 2015

Author: Daisy Leung

History:
Dec 09 2015
- created this script to overplot SMA with CO
- need cube /Users/admin/Research/hatlas_Cooray/OverlayFig/062015/sendthem/CD_USB_2.685_S.bin5.co32.mom0.fits for CO and /Users/admin/Research/hatlas_Cooray/SMA/NA.V1.489.SMA.SUB.DSB.FITS
- need to use case-insenstive glob because DR's file are uppercase

'''
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from FgContIm import *


def insensitive_glob(pattern):
    '''
        Case insensitive file pattern match
    '''
    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(either, pattern)))


COpath = '/Users/admin/Research/hatlas_Cooray/OverlayFig/062015/sendthem/'
SMApath = '/Users/admin/Research/hatlas_Cooray/SMA/'
Plotpath = '../OverlayFig/Dec0915/'


# 'S.bin5.co32.mom0'    # 0.680
# 'bin5.co32.mom0', sigma = 0.753
# 'bin5.32_43co32.mom0', sigma = 0.928 (Personally prefer this over bin5.co32.mom0, i.e. the entry above)
COitermomName = 'S.bin5.co32.mom0'


label = dict.fromkeys([COitermomName, 'SMA'])
path = [COpath, SMApath]

for i, k in enumerate(label.iterkeys()):
    print path[i]
    files = insensitive_glob(path[i] + '*' + k + '*.fits')
    label[k] = files

print label

# import sys
#sys.exit()

########################################
# user define area
########################################
# 13:35:42.245, 30:04:06.854
ra_center = 203.927
dec_center = 30.06857
# positions for crosses
ra_cross, dec_cross = ra_center, dec_center

sigma_CO = 0.680  # 0.928      #0.753      #0.680   per beam
sigma_SMA = 2.9E-04

sizep = 0.007

fig = plt.figure(1, figsize=(6, 6))

row_a = 0.10
width = 0.35
x_gap = 0.05
x0 = 0.10
dy = 0.90


########################################
# intialize base figure
########################################
figCont = aplpy.FITSFigure(label['SMA'][0],
                           figure=fig, subplot=[x0, row_a, 0.8, dy])
figCont.set_system_latex(True)
# figCont.show_grayscale()


########################################
# Contours
########################################

figCont.show_contour(label['SMA'][0], colors='red', levels=sigma_contour_array(
    sigma_SMA), linewidths=1.5, layer='Cont')
figCont.show_contour(label[COitermomName][0], colors="blue", levels=sigma_contour_CARMA(
    sigma_CO), linewidths=1.5, layer='mol')

########################################
# beam
########################################

# setup_beam(figCont)
# for some reason auto setup doesnt for SMA data

# the following hack didn't work either
# rad2deg = 57.2958
# BMAJ = 1.2278E-03 * rad2deg
# BMIN = 8.1474E-04 * rad2deg
# BPA = 43.82

# figCont.add_beam(0, 0, 0)
# figCont.beam.show()
# figCont.beam.set_major(BMAJ)  # degrees
# figCont.beam.set_minor(BMIN)  # degrees
# figCont.beam.set_angle(BPA)  # degrees

# fig.beam.set_corner('bottom left')
# fig.beam.set_frame(False)
# fig.beam.set_alpha(0.5)
# fig.beam.set_color('black')
# fig.beam.set_edgecolor('white')
# fig.beam.set_facecolor('green')
# fig.beam.set_linestyle('dashed')
# fig.beam.set_linewidth(2)  # points

########################################
# scale bar
########################################
lg_1arcsec = 1. / 3600
# lg_20kpc_fg = lg_1arcsec * 20./scale_radio
# lg_20kpc_bg = lg_1arcsec * 20./scale_SMG
# setup_scalebar(fvla, lg_20kpc_fg, str('20kpc'))
# setup_scalebar(fcont, lg_20kpc_fg, str('20kpc'))
# setup_scalebar(flin, lg_20kpc_bg, str('20kpc'))
# setup_scalebar(fSMA, lg_20kpc_bg, str('20kpc'))

########################################
# axes
########################################
standard_plot_setup(figCont, ra_center, dec_center, sizep)

# fcont.tick_labels.hide()
# fcont.axis_labels.hide()

########################################
# markers
########################################
# markers_cross(figCont, ra_cross, dec_cross, layer='marker_set_1')


########################################
# Labels
########################################
labsize = 'xx-large'
put_label(figCont, 0.25, 0.9, 'NA.v1.489', 'titleObj', c='black', s=labsize)
put_label(figCont, 0.25, 0.95, 'SMA 228 GHz', 'titleBand', c='red', s=labsize)
put_label(figCont, 0.25, 0.85, 'CARMA CO(3-2)', 'titleBand2', c='blue', s=labsize)

put_label(figCont, 0.75, 0.95, 'SMA: +/-2^n' + r'$\times\sigma$', 'titleerrSMA', c='black', w=1000)
put_label(figCont, 0.75, 0.85, 'CO: +/-2,3,...' + r'$\times\sigma$', 'titleerrCO', c='black', w=1000)

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
    import sys
    import os
    if len(sys.argv) < 2:
        errmsg = "Invalid number of arguments: {0:d}\n  run script.py Save_True"
        raise IndexError(errmsg.format(len(sys.argv)))
    saveFig = True if sys.argv[1].lower() == 'true' else False
    if saveFig == True:
        os.system('rm -rf ' + itermomName + '.png' + ' ' + itermomName + '.eps')
#        figCD.savefig(Plotpath + CD + '.eps', dpi=600)
        print ("saving to {:s}").format(Plotpath + itermomName + '.png')
        figCD.savefig(Plotpath + itermomName + '.png', dpi=300)
        fig.savefig(Plotpath + 'Cont.png', dpi=300)
    else:
        plt.show()
