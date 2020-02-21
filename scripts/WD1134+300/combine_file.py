import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'/Users/Ryan/research/pipeline/')
import utilities as uti
import process as pro
import pre_process as pre
import numpy as np
from astropy.io import fits

master_path_1 = '/Users/Ryan/research/data/WD1134+300/01-05-18/'
path_final_1 = master_path_1+'final/'
master_path_2 = '/Users/Ryan/research/data/WD1134+300/02-08-17/'
path_final_2 = master_path_2+'final/'
return_path = '/Users/Ryan/research/data/WD1134+300/'

##median psf subtracted
#file_one_msm = fits.open('{}/median_psf_subtraction/median_sub_master.fits'.format(path_final_1))[0].data.astype(float)
#file_two_msm = fits.open('{}/median_psf_subtraction/median_sub_master.fits'.format(path_final_2))[0].data.astype(float)
#
#combine_msm = np.median([file_one_msm,file_two_msm],axis=0)
#
#out_fits = fits.HDUList(fits.PrimaryHDU(combine_msm))
#out_fits.writeto('{}median_sub_comb.fits'.format(return_path), overwrite = True)
#
##deroto
#file_one_de = fits.open('{}/median_psf/master_deroto.fits'.format(path_final_1))[0].data.astype(float)
#file_two_de = fits.open('{}/median_psf/master_deroto.fits'.format(path_final_2))[0].data.astype(float)
#
#combine_de = np.median([file_one_de,file_two_de],axis=0)
#
#out_fits = fits.HDUList(fits.PrimaryHDU(combine_de))
#out_fits.writeto('{}deroto_comb.fits'.format(return_path), overwrite = True)
#
##no deroto
#file_one_nd = fits.open('{}/median_psf/master_no_deroto.fits'.format(path_final_1))[0].data.astype(float)
#file_two_nd = fits.open('{}/median_psf/master_no_deroto.fits'.format(path_final_2))[0].data.astype(float)
#
#combine_nd = np.median([file_one_nd,file_two_nd],axis=0)
#
#out_fits = fits.HDUList(fits.PrimaryHDU(combine_nd))
#out_fits.writeto('{}no_deroto_comb.fits'.format(return_path), overwrite = True)

star_mag = 13.105
plot_radius = 330
image_center = 331

uti.contrast_curves(return_path,'deroto_comb',return_path,image_center,star_mag,'g',plot_radius,find_star_counts=True)
uti.contrast_curves(return_path,'no_deroto_comb',return_path,image_center,star_mag,'b',plot_radius)
uti.contrast_curves(return_path,'median_sub_comb',return_path,image_center,star_mag,'b',plot_radius)

