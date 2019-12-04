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

def edge_crop(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,file_type,test_type='',crop_bottom=80,crop_top=1200,crop_left=650,crop_right=2042):
    #Function for cropping images
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
        fname_edge_crop = fname_float[crop_bottom:crop_top,crop_left:crop_right]
        out_fits = fits.HDUList(fits.PrimaryHDU(fname_edge_crop))
        out_fits.writeto('{}{}.fits'.format(return_path,fname_prefix), overwrite = True)
        if progress_print == True:
            print('{} edge cropped'.format(fname_prefix))
        fname.close()
    print('{} frames edge cropping complete'.format(file_type))

def stacker(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,file_type,test_type=''):
    #Function for median stacking images
    frames,fname_prefixes = data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type,progress_print)
    frames_array = np.array(frames, dtype=np.float)
    stacked_frames = np.median(frames_array, axis=0)
    out_fits = fits.HDUList(fits.PrimaryHDU(stacked_frames))
    out_fits.writeto('{}stacked_{}.fits'.format(return_path,file_type), overwrite = True)
    del frames
    del fname_prefixes
    print('Stacked {} frame created'.format(file_type))
    return(stacked_frames)

def bp_mask_mkr(master_path,return_path,sigma_coe,test_type=''):
    #This function creates a mask of bad pixels from a median combined dark frame

    flt_drk = fits.open('{}reduced_data{}/flt_drk.fits'.format(master_path,test_type))
    #flt_drk = fits.open('{}reduced_data{}/stacked_drk.fits'.format(master_path,test_type))
    flt_drk_float = flt_drk[0].data.astype(float)
    flt_drk_blur = gaussian_filter(flt_drk_float, 7)
    flt_drk_blur_sub = flt_drk_float - flt_drk_blur
    bp_list_neg = np.where(flt_drk_blur_sub <= -sigma_coe*np.std(flt_drk_blur_sub))
    bp_list_pos = np.where(flt_drk_blur_sub >= sigma_coe*np.std(flt_drk_blur_sub))
    bp_list_negx = list(bp_list_neg[0])
    bp_list_posx = list(bp_list_pos[0])
    bp_list_negy = list(bp_list_neg[1])
    bp_list_posy = list(bp_list_pos[1])
    bp_list_totx = bp_list_posx + bp_list_negx
    bp_list_toty = bp_list_posy + bp_list_negy
    bp_list = [bp_list_totx,bp_list_toty]
    zero_mask = np.zeros(np.shape(flt_drk_float))
    zero_mask[bp_list] = 1.0
    bp_mask = zero_mask
    out_fits = fits.HDUList(fits.PrimaryHDU(bp_mask))
    out_fits.writeto('{}bp_mask.fits'.format(return_path), overwrite = True)
    out_fits = fits.HDUList(fits.PrimaryHDU(bp_list))
    out_fits.writeto('{}bp_list.fits'.format(return_path), overwrite = True)
    flt_drk.close()
    print('Bad pixel mask created')

def flat_field_mkr(master_path,prefix,aux_path_drk,aux_path_flt,start_frame_num_drk,end_frame_num_drk,start_frame_num_flt,end_frame_num_flt,return_path,progress_print,test_type=''):
    #Function takes indivudual flat frames and median combines them into a master flat frame

    master_drk_float = stacker(master_path,prefix,aux_path_drk,start_frame_num_drk,end_frame_num_drk,return_path,progress_print,file_type='drk',test_type=test_type)
    master_flt_float = stacker(master_path,prefix,aux_path_flt,start_frame_num_flt,end_frame_num_flt,return_path,progress_print,file_type='flt',test_type=test_type)
    #master_drk_float = master_drk[0].data.astype(float)
    #master_flt_float = master_flt[0].data.astype(float)
    flt_drk = master_flt_float - master_drk_float
    flt_drk_median = flt_drk/np.median(flt_drk)
    out_fits = fits.HDUList(fits.PrimaryHDU(flt_drk))
    out_fits.writeto('{}flt_drk.fits'.format(return_path), overwrite = True)
    out_fits = fits.HDUList(fits.PrimaryHDU(flt_drk_median))
    out_fits.writeto('{}flt_drk_median.fits'.format(return_path), overwrite = True)
    #master_drk.close()
    #master_flt.close()
    print('Flat corrected dark frame created')

def bp_crct(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,file_type,test_type='',iterations=3,blurring=7): 
    #Custom function for eliminating bad pixels. It uses a mask to find bad pixels, gaussian blurs them,
    # and then replaces the original bad pixels with blured ones. Its more effective at eliminating
    # bad pixels than replacing bad pixels with the average of its surrounding pixels.

    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    bp_mask = fits.open('{}reduced_data{}/bp_mask.fits'.format(master_path,test_type))[0].data.astype(int)
    bp_list = fits.open('{}reduced_data{}/bp_list.fits'.format(master_path,test_type))[0].data.astype(int)
    bp_list = tuple(bp_list)
    if file_type == 'sci': #bad pixel correcting sci frames
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
            for j in range(0,iterations):
                fname_blur = gaussian_filter(fname_float, blurring)
                fname_float[bp_list] = fname_blur[bp_list]
            fname_float = vip.preproc.badpixremoval.frame_fix_badpix_isolated(fname_float,bp_mask,verbose=False)
            out_fits = fits.HDUList(fits.PrimaryHDU(fname_float))
            out_fits.writeto('{}{}.fits'.format(return_path,fname_prefix), overwrite = True)
            if progress_print == True:
                print('{} bad pixel corrected'.format(fname_prefix))
            fname.close()
        print('{} frames bad pixel correcting complete'.format(file_type))
    else: #this is for DRK - FLT frame
        fname = fits.open('{}reduced_data{}/flt_drk_median.fits'.format(master_path,test_type))
        fname_float = fname[0].data.astype(float)
        for j in range(0,iterations):
            fname_blur = gaussian_filter(fname_float, blurring)
            fname_float[bp_list] = fname_blur[bp_list]
        #fname_float = vip.preproc.badpixremoval.frame_fix_badpix_isolated(fname_float,bp_mask,verbose=False)
        out_fits = fits.HDUList(fits.PrimaryHDU(fname_float))
        out_fits.writeto('{}{}_bp_crct.fits'.format(return_path,file_type), overwrite = True)
        fname.close()
        print('Master {} frame bad pixel correcting complete'.format(file_type))
        fname.close

def flat_fielder(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print,test_type=''): 
    #function for flat fielding images
    
    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')
    flt_drk_median = fits.open('{}reduced_data{}/flt_drk_median_bp_crct.fits'.format(master_path,test_type))
    flt_drk_median_float = flt_drk_median[0].data.astype(float)
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
        fname_med = fname_float/flt_drk_median_float
        out_fits = fits.HDUList(fits.PrimaryHDU(fname_med))
        out_fits.writeto('{}{}.fits'.format(return_path,fname_prefix), overwrite = True)
        if progress_print == True:
            print('{} flat fielded'.format(fname_prefix))
        fname.close()
    flt_drk_median.close()
    print('All sci frames flat fielded')

