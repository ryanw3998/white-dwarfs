#!/usr/bin/python
import os
import math
import numpy as np
import vip_hci as vip
from photutils import CircularAperture as ca
from photutils import aperture_photometry as ap
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.filters import gaussian_filter
from vip_hci.preproc import cube_recenter_2dfit
import pickle
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'/Users/Ryan/research/pipeline/')
from utilities import data_opener

def nod_subtract(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,test_type='',luci=False):
    #We take data on different parts of the detector, and then subtract one set of images from the 
    # median of the other set and vise versa. This eliminates the majority of bad pixels and background noise

    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    nods = []
    fnames,fname_prefixes = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print,luci)
    for fname_prefix in fname_prefixes:
        if fits.open('{}raw_data/{}.fits'.format(master_path,fname_prefix))[0].header['FLAG'] == 'NOD_A':
            nods.append('NOD_A')
        else: 
            nods.append('NOD_B')
    fnames_np = np.array(fnames, dtype=np.float)
    fname_prefixes_np = np.array(fname_prefixes)
    line = 0
    end = False
    while end == False:
        chnod = False
        chnod_2 = False
        fnames_nod_A = []
        fname_prefixes_nod_A = []
        while chnod == False:
            fnames_nod_A.append(fnames_np[line])
            fname_prefixes_nod_A.append(fname_prefixes_np[line])
            try:
                if nods[line+1] != nods[line]:
                    chnod = True
            except:
                IndexError
                chnod = True
            line += 1
        fnames_nod_B = []
        fname_prefixes_nod_B = []
        while chnod_2 == False:
            fnames_nod_B.append(fnames_np[line])
            fname_prefixes_nod_B.append(fname_prefixes_np[line])
            try:
                if nods[line+1] != nods[line]:
                    chnod_2 = True
            except:
                IndexError
                chnod_2 = True
            line += 1
        fnames_nod_A_med = np.median(fnames_nod_A,axis=0)
        fnames_nod_B_med = np.median(fnames_nod_B,axis=0)
        for i in range(len(fnames_nod_A)):
            fname_sub_A = fnames_nod_A[i] - fnames_nod_B_med
            out_fits = fits.HDUList(fits.PrimaryHDU(fname_sub_A))
            out_fits.writeto('{}{}.fits'.format(return_path,fname_prefixes_nod_A[i]), overwrite = True)
            if progress_print == True:
                print ('{} sci subtracted'.format(fname_prefixes_nod_A[i]))
        for i in range(len(fnames_nod_B)):
            fname_sub_B = fnames_nod_B[i] - fnames_nod_A_med
            out_fits = fits.HDUList(fits.PrimaryHDU(fname_sub_B))
            out_fits.writeto('{}{}.fits'.format(return_path,fname_prefixes_nod_B[i]), overwrite = True)
            if progress_print == True:
                print ('{} sci subtracted'.format(fname_prefixes_nod_B[i]))
        if line == len(nods):
            end = True
    print('All frames nod subtracted')
    del fnames_nod_A
    del fname_prefixes_nod_A
    del fnames_nod_B
    del fname_prefixes_nod_B
    del fnames
    del fname_prefixes
    del nods
    del fnames_np
    del fname_prefixes_np

def star_center(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,star_location,crop_parameter,test_type=''):
    #Roughly centers star based on pixel coordinates of brightest pixel

    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    nods = []
    fnames,fname_prefixes = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print)
    for fname_prefix in fname_prefixes:
        if fits.open('{}raw_data/{}.fits'.format(master_path,fname_prefix))[0].header['FLAG'] == 'NOD_A': #changed to fits.gz
            nods.append('NOD_A')
        else: 
            nods.append('NOD_B')
    fnames_np = np.array(fnames, dtype=np.float)
    fname_prefixes_np = np.array(fname_prefixes)
    line = 0
    end = False
    while end == False:
        chnod = False
        fnames_nod = []
        fname_prefixes_nod = []
        while chnod == False:
            fnames_nod.append(fnames_np[line])
            fname_prefixes_nod.append(fname_prefixes_np[line])
            try:
                if nods[line+1] != nods[line]:
                    chnod = True
            except:
                IndexError
                chnod = True
            line += 1
        fnames_med = np.median(fnames_nod,axis=0)
        #out_fits = fits.HDUList(fits.PrimaryHDU(fnames_med))
        #out_fits.writeto('{}nod_med_comp_{}.fits'.format(return_path,line), overwrite = True)
        fnames_blur = gaussian_filter(fnames_med, 7)
        star = np.unravel_index(np.argmax(fnames_blur), fnames_blur.shape)
        if star_location == True:
            print('star',star)
        for i in range(len(fnames_nod)):
            fname_crop = fnames_nod[i][star[0]-crop_parameter:star[0]+crop_parameter+1,star[1]-crop_parameter:star[1]+crop_parameter+1] 
            out_fits = fits.HDUList(fits.PrimaryHDU(fname_crop))
            out_fits.writeto('{}{}.fits'.format(return_path,fname_prefixes_nod[i]), overwrite = True)
            if progress_print == True:
                print ('{} star centered'.format(fname_prefixes_nod[i]))
        if line == len(nods):
            end = True
    print('All frames star centered')
    del nods
    del fnames_np
    del fname_prefixes_np
    del fnames_nod
    del fname_prefixes_nod

def stripe_cleaner(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,test_type='',mask_size=6000,luci=False):
    #Cleans detector stripes in images by subtracting the meadian of the entire column from each pixel
    # in the column

    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        try:
            if aux_path == 'raw_data':
                fname = fits.open('{}{}/{}.fits'.format(master_path,aux_path,fname_prefix))
                fname_float = fname[0].data.astype(float)
            else:
                fname = fits.open('{}reduced_data{}/{}/{}.fits'.format(master_path,test_type,aux_path,fname_prefix))
                fname_float = fname[0].data.astype(float)
            if progress_print == True:
                print('{} loaded'.format(fname_prefix))
        except:
            FileNotFoundError
            print('File Not Found Error')
            continue
        #fname_uncrcted = fname.copy()
        fname_star = fname_float.copy()
        y_dim,x_dim = fname_float.shape
        fname_blur = gaussian_filter(fname_float, 7)
        peak = np.unravel_index(np.argmax(fname_blur), fname_blur.shape) #finding star in blurred frame
        x = [i for i in range(x_dim)] #setting up mask
        y = [i for i in range(y_dim)]
        X,Y = np.meshgrid(x,y)
        x_0 = X[peak[0],peak[1]]
        y_0 = Y[peak[0],peak[1]]
        Mask = ((X-x_0)**2+(Y-y_0)**2)<mask_size #making circular mask around star (star is masked out)
        Mask_2 = ((X-x_0)**2+(Y-y_0)**2)>mask_size #area around star is masked out
        fname_float[Mask] = np.nan #star is replaced by nans
        fname_star[Mask_2] = np.nan #area surrounding stars replaced by nans
        one_d_array = np.nanmedian(fname_float, axis=0) #medians columns, ignoring nans
        two_d_array = np.tile(one_d_array,(y_dim,1))
        fname_stripe_crct = fname_float - two_d_array
        fname_star_stripe_crct = fname_star - two_d_array 
        fname_stripe_crct[Mask] = fname_star_stripe_crct[Mask]
        out_fits = fits.HDUList(fits.PrimaryHDU(fname_stripe_crct))
        out_fits.writeto('{}{}.fits'.format(return_path,fname_prefix), overwrite = True)
        if progress_print == True:
            print('{} stripe corrected'.format(fname_prefix))
        del Mask
        del x
        del y
        fname.close()
    print('All frames stripe corrected')

def frame_filter(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,fwhm_tolerance,peak_tolerance,luci=False,test_type=''):
    #Filters frame based on FWHM size and peak brightness 
    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    try:
        os.remove('{}reduced_data/fwhm_median.p'.format(master_path))
    except:
        FileNotFoundError
    cube,fname_prefixes = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print,luci)
    for fname_prefix in fname_prefixes: #Removing exsisting frames from directory might need this for functions called after frame filter
        try:
            os.remove('{}{}.fits'.format(return_path,fname_prefix))
        except:
            FileNotFoundError
    cube = np.array(cube,dtype=np.float)
    fname_prefixes = np.array(fname_prefixes)
    fwhms = []
    fname_prefixes_fwhm_filtered = []
    #Filtering frames based on fwhm size
    for fname,fname_prefix in zip(cube,fname_prefixes):
        fwhm = vip.var.fit_2dgaussian(fname, crop=True, cropsize=9, full_output=True, debug=False)
        fwhm_mean = np.mean([fwhm.loc[0,'fwhm_x'],fwhm.loc[0,'fwhm_y']])
        if progress_print == True:
            print('{} fwhm_mean = {}'.format(fname_prefix,fwhm_mean))
        fwhms.append(fwhm_mean)
    fwhms = np.array(fwhms)
    #Finding smallest fwhms 
    fwhm_tolerance = fwhm_tolerance*0.01
    print(len(cube))
    g = int(np.ceil(len(cube)*fwhm_tolerance)) #make comment describing the effects of ceiling 
    #g = int(np.ceil((end_frame_num-start_frame_num)*fwhm_tolerance)) #Kept for refrence
    if progress_print == True:
        print('{} frames within fwhm tolerance'.format(g))
    fwhm_filter = np.argpartition(fwhms, g)[:g] #returns indicies of smallest fwhms
    cube_fwhm_filtered = cube[fwhm_filter]
    fname_prefixes_fwhm_filtered = fname_prefixes[fwhm_filter] #keeping track of prefixes
    #Finding peak values in frames that passed fwhm filtering
    peak_list = []
    for fname,fname_prefix in zip(cube_fwhm_filtered,fname_prefixes_fwhm_filtered): 
        peak = np.amax(fname)
        peak_list.append(peak)
        if progress_print == True:
            print('{} peak = {}'.format(fname_prefix,peak))
    peak_list = np.array(peak_list)
    #Finding frames with highest peaks
    peak_tolerance = peak_tolerance*0.01
    k = int(np.ceil((len(peak_list))*peak_tolerance))
    if progress_print == True:
        print('{} frames within peak tolerance'.format(k))
    peak_filter = np.argpartition(peak_list, -k)[-k:]
    cube_peak_filtered = cube_fwhm_filtered[peak_filter]
    fname_prefixes_peak_filtered = fname_prefixes_fwhm_filtered[peak_filter]
    #Frames are now peak filtered
    for fname,fname_prefix in zip(cube_peak_filtered,fname_prefixes_peak_filtered):
        out_fits = fits.HDUList(fits.PrimaryHDU(fname))
        out_fits.writeto('{}{}.fits'.format(return_path,fname_prefix), overwrite = True)
    #Finding the median peak value across all frames
    cube_fwhm_median = np.median(fwhms)
    cube_fwhm_median_array = np.array([cube_fwhm_median])
    np.savetxt('{}reduced_data{}/fwhm_median.txt'.format(master_path,test_type),cube_fwhm_median_array)
    #pickle.dump(cube_fwhm_median, open('{}reduced_data{}/fwhm_median.p'.format(master_path,test_type),'wb'))
    print('{} frames passed filtering'.format(len(cube_peak_filtered)))
    
    del cube
    del peak_list
    del fname_prefixes
    del cube_peak_filtered
    del fname_prefixes_peak_filtered
    del cube_fwhm_filtered
    del fwhms
    del fname_prefixes_fwhm_filtered

def star_ultra_center(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,crop_parameter,luci=False,test_type=''):
    #Centers each image via interpolation
    if os.path.isdir(return_path):
        print('Deleting directory')
    else:
        os.mkdir(return_path)
        print('Making directory')
    filelist = [f for f in os.listdir(return_path) if f.endswith('.fits')]
    for f in filelist:
        os.remove(os.path.join(return_path,f))
    fwhm_load = np.loadtxt('{}reduced_data{}/fwhm_median.txt'.format(master_path,test_type))
    fwhm_load = float(fwhm_load)
    print(fwhm_load)
    cube,fname_prefixes = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print,luci)
    cube_np_array = np.array(cube, dtype=np.float)
    cube_ultra_center = cube_recenter_2dfit(cube_np_array, (crop_parameter+1,crop_parameter+1), fwhm = fwhm_load,  nproc=1, subi_size=20, model='gauss', negative=False, full_output=False, debug=False)
    for fname_ultra_centered,fname_prefix_ultra_centered in zip(cube_ultra_center,fname_prefixes):
        out_fits = fits.HDUList(fits.PrimaryHDU(fname_ultra_centered))
        out_fits.writeto('{}{}.fits'.format(return_path,fname_prefix_ultra_centered), overwrite = True)
        if progress_print == True:
            print('{} ultra centered'.format(fname_prefix_ultra_centered))
    print('All frames ultra centered')
    del cube
    del fname_prefixes
    del cube_ultra_center

def llsg_psf_subtraction(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,frame_set,rank_input,thresh_input,max_iter_input,luci=False,test_type=''):
    #PSF subtraction based off LLSG algorithm published in Gomez et al 2018 (https://arxiv.org/pdf/1602.08381.pdf)
    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    cube,fname_prefixes,parangs = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print,return_parang=True,luci=luci)
    fwhm_load = np.loadtxt('{}reduced_data{}/fwhm_median.txt'.format(master_path,test_type))
    #fwhm_load = pickle.load(open('{}reduced_data{}/fwhm_median.p'.format(master_path,test_type),'rb'))
    cube_np_array = np.array(cube, dtype=np.float)
    parangs_np = np.array(parangs, dtype=np.float)
    psf_subtraction = vip.llsg.llsg(cube_np_array, parangs_np, fwhm_load, rank=rank_input, thresh=thresh_input,max_iter=max_iter_input,random_seed=10,full_output=False)
    out_fits = fits.HDUList(fits.PrimaryHDU(psf_subtraction))
    out_fits.writeto('{}llsg_sub_{}.fits'.format(return_path,frame_set), overwrite = True)
    print('Star evaluated using llsg method')
    del cube
    del parangs
    del fname_prefixes

def median_psf_subtraction(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,frame_set,luci=False,test_type=''):
    #Conventional psf subtraction algorithm
    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    cube,fname_prefixes,parangs = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print,return_parang=True,luci=luci)

    a_frames = []
    a_frames_pa = []
    b_frames = []
    b_frames_pa = []
    for i,fname_prefix,parang in zip(cube,fname_prefixes,parangs):
        if fits.open('{}raw_data/{}.fits'.format(master_path,fname_prefix))[0].header['FLAG'] == 'NOD_A':
            a_frames.append(i)
            a_frames_pa.append(parang)
            print('a')
        else:
            b_frames.append(i)
            b_frames_pa.append(parang)
            print('b')
    print(len(a_frames),len(b_frames))
    a_frames_np = np.array(a_frames, dtype=np.float)
    a_frames_pa_np = np.array(a_frames_pa, dtype=np.float)
    b_frames_np = np.array(b_frames, dtype=np.float)
    b_frames_pa_np =np.array(b_frames_pa, dtype=np.float)
    fwhm_load = np.loadtxt('{}reduced_data{}/fwhm_median.txt'.format(master_path,test_type))


    psf_subtraction_a = vip.medsub.median_sub(a_frames_np, a_frames_pa_np, fwhm=fwhm_load, mode='fullfr',full_output=False)
    out_fits = fits.HDUList(fits.PrimaryHDU(psf_subtraction_a))
    out_fits.writeto('{}median_sub_{}_a.fits'.format(return_path,frame_set), overwrite = True)

    psf_subtraction_b = vip.medsub.median_sub(b_frames_np, b_frames_pa_np, fwhm=fwhm_load, mode='fullfr',full_output=False)
    out_fits = fits.HDUList(fits.PrimaryHDU(psf_subtraction_b))
    out_fits.writeto('{}median_sub_{}_b.fits'.format(return_path,frame_set), overwrite = True)

    
    ##print(parangs)
    #fwhm_load = np.loadtxt('{}reduced_data{}/fwhm_median.txt'.format(master_path,test_type))
    ##fwhm_load = pickle.load(open('{}reduced_data{}/fwhm_median.p'.format(master_path,test_type),'rb'))
    #cube_np_array = np.array(cube, dtype=np.float)
    #parangs_np = np.array(parangs, dtype=np.float)
    #psf_subtraction = vip.medsub.median_sub(cube_np_array, parangs_np, fwhm=fwhm_load, mode='fullfr',full_output=False)
    #out_fits = fits.HDUList(fits.PrimaryHDU(psf_subtraction))
    #out_fits.writeto('{}median_sub_{}.fits'.format(return_path,frame_set), overwrite = True)
    #print('Star evaluated using median psf subtraction method')
    #del cube
    #del parangs
    #del fname_prefixes

def psf_stacker(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,frame_set,luci=False,test_type=''):
    #Median combines images together

    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    fnames,fname_prefixes,parangs = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print,return_parang=True,luci=luci)
    fnames_median = np.median(fnames,axis=0)
    fnames_deroto = []
    for fname,parang,fname_prefix in zip(fnames,parangs,fname_prefixes):
        fname_deroto = rotate(fname, angle = -parang, cval = np.nan, reshape=False)
        if progress_print == True:
            print('{} derotated'.format(fname_prefix), parang)
        fnames_deroto.append(fname_deroto)
    fnames_median_deroto = np.median(fnames_deroto,axis=0)
    out_fits = fits.HDUList(fits.PrimaryHDU(fnames_median))
    out_fits.writeto('{}{}_no_deroto.fits'.format(return_path,frame_set), overwrite = True)
    out_fits = fits.HDUList(fits.PrimaryHDU(fnames_median_deroto))
    out_fits.writeto('{}{}_deroto.fits'.format(return_path,frame_set), overwrite = True)
    print('Star median combined')
    print('Star derotated and median combined')
    del fnames
    del fname_prefixes
    del parangs

def contrast_curve_test(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,frame_set,luci=False,test_type='',rank=1):
    #Improved Contrast curve using algorithm from VIP, under construction

    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    fnames,fname_prefixes,parangs = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print,return_parang=True,luci=luci)
    del fname_prefixes

    parangs = np.array(parangs,dtype=np.float)
    fnames_array = np.array(fnames,dtype=np.float)
    fnames_med = np.median(fnames_array,axis=0)
    fwhm = vip.var.fit_2dgaussian(fnames_med, crop=True, full_output=True)
    fwhm_mean = np.mean([fwhm.loc[0,'fwhm_x'],fwhm.loc[0,'fwhm_y']])
    #print(fwhm)
    print('FWHM = {}'.format(fwhm_mean))

    image_center = int((fnames_med.shape)[0]/2 + .5) #Finding image center
    print('Image Center {}'.format(image_center))
    flux = vip.metrics.aperture_flux(fnames_med,[image_center],[image_center],fwhm=fwhm_mean)
    print('Star Flux {}'.format(flux[0]))

    if luci == True:
        plate_scale = 0.015
    else:
        plate_scale = 10.7/1000.
    print('Plate Scale {}'.format(plate_scale))

    path = '{}final/median_psf_subtraction/contrast_curve_test.png'.format(master_path)
    #path = '{}final/llsg_psf_subtraction/contrast_curve_rank_{}.png'.format(master_path,rank)
    
    data = vip.metrics.contrast_curve(fnames_array,parangs,fnames_med,fwhm_mean,plate_scale,flux[0],vip.medsub.median_sub,plot=True,save_plot=path,debug=True,full_output=True)
    print(data[0])
    print(data[1].shape)
    print(data[2].shape)
    print(data[3].shape)
    #print(type(data[4]))
    #vip.metrics.contrast_curve(fnames,parangs,psf,fwhm,plate_scale,fwhm_flux[0],vip.llsg.llsg,plot=True,save_plot=path,debug=True)
    
    plt.show()