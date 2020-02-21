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
from dewarp_package.astrom_lmircam_soln import *
from dewarp_package.astrom_lmircam_soln import dewarp


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

def dewarper(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,crop_bottom,crop_top,crop_left,crop_right,date,side,test_type='',progress_print=False,luci=False):
    if os.path.isdir(return_path):
        print('Directory exists')
    else:
        os.mkdir(return_path)
        print('Making directory')

    if date == '2016B' or '2017A':
        #dewarp_file = '/Users/Ryan/research/pipeline/dewarp_package/coeff_lmir_double_161122.dat'
        Kx = [[ -4.74621436e+00,   9.99369200e-03,  -4.69741638e-06,   4.11937105e-11],
              [  1.01486148e+00,  -2.84066638e-05,   2.10787962e-08,  -3.90558311e-12],
              [ -1.61139243e-05,   2.24876212e-08,  -2.29864156e-11,   6.59792237e-15],
              [  8.88888428e-09,  -1.03720381e-11,   1.05406782e-14,  -3.06854175e-18]]
        Ky = [[  9.19683947e+00,   9.84613002e-01,  -1.28813904e-06,   6.26844974e-09],
              [ -7.28218373e-03,  -1.06359740e-05,   2.43203662e-09,  -1.17977589e-12],
              [  9.48872623e-06,   1.03410741e-08,  -2.38036199e-12,   1.17914143e-15],
              [  3.56510910e-10,  -5.62885797e-13,  -5.67614656e-17,  -4.21794191e-20]]

    elif date == '2017B' or '2018A':
        if side == 'dx':
            #dewarp_file = '/Users/Ryan/research/pipeline/dewarp_package/coeff_lmir_dx_171108.dat'
            Kx = [[ -1.34669677e+01,   2.25398365e-02,  -7.39846082e-06,  -8.00559920e-11],
                  [  1.03267422e+00,  -1.10283816e-05,   5.30280579e-09,  -1.18715846e-12],
                  [ -2.60199694e-05,  -3.04570646e-09,   1.12558669e-12,   1.40993647e-15],
                  [  8.14712290e-09,   9.36542070e-13,  -4.20847687e-16,  -3.46570596e-19]]
            Ky = [[  1.43440109e+01,   9.90752231e-01,  -3.52171557e-06,   7.17391873e-09],
                  [ -2.43926351e-02,  -1.76691374e-05,   5.69247088e-09,  -2.86064608e-12],
                  [  1.06635297e-05,   8.63408955e-09,  -2.66504801e-12,   1.47775242e-15],
                  [ -1.10183664e-10,  -1.67574602e-13,   2.66154718e-16,  -1.13635710e-19]]
        elif side == 'sx':
            #dewarp_file = '/Users/Ryan/research/pipeline/dewarp_package/coeff_lmir_sx_171108.dat'
            Kx = [[ -1.97674665e+01,   2.26890756e-02,  -1.06483884e-05,   1.33174251e-09],
                  [  1.04269459e+00,  -2.68457747e-05,   2.08519317e-08,  -4.74786541e-12],
                  [ -3.30919802e-05,   9.48855944e-09,  -1.00804780e-11,   3.45894384e-15],
                  [  1.00196745e-08,  -2.58289058e-12,   2.58960909e-15,  -8.74827083e-19]]
            Ky = [[  1.05428609e+01,   9.91877631e-01,  -1.30947328e-06,   5.98620664e-09],
                  [ -2.65330464e-02,  -6.14857421e-06,  -1.56796197e-08,   6.61213303e-12],
                  [  1.50777505e-05,  -8.14931285e-09,   2.28968428e-11,  -9.37645995e-15],
                  [ -1.54162134e-09,   5.25556977e-12,  -7.46189515e-15,   3.04540450e-18]]
        else:
            #dewarp_file = '/Users/Ryan/research/pipeline/dewarp_package/coeff_lmir_double_171108.dat'
            Kx = [[ -1.28092986e+01,   2.37843229e-02,  -8.09328998e-06,  -1.28301567e-10],
                  [  1.02650326e+00,  -2.15600674e-05,   1.42563655e-08,  -3.38050198e-12],
                  [ -1.94200903e-05,   3.51059572e-09,  -4.80126009e-12,   2.96495549e-15],
                  [  6.26453161e-09,  -4.65825179e-13,   9.00810316e-16,  -7.28975474e-19]]
            Ky = [[  1.47852158e+01,   9.92480569e-01,  -6.05953133e-06,   8.04550607e-09],
                  [ -3.22153904e-02,  -1.77238376e-05,   1.09565607e-08,  -5.03884071e-12],
                  [  1.90672398e-05,   3.96777274e-09,  -2.84884006e-12,   1.98359721e-15],
                  [ -2.57086429e-09,   1.87084793e-12,  -3.96206006e-16,  -7.81470678e-20]]

    elif date == '2018B' or '2019A' or '2019B' or '2020A':
        if side == 'dx':
            #dewarp_file = '/Users/Ryan/research/pipeline/dewarp_package/coeff_lmir_dx_20200102.rtf'
            Kx = [[-1.13034544e+01,   1.45852226e-02,  -1.13372175e-05,   1.32602063e-09],
 	              [ 1.03220943e+00,  -1.93058980e-05,   1.55798844e-08,  -3.86115281e-12],
 	              [-2.57352199e-05,   2.70371257e-09,  -6.62650011e-12,   3.25017078e-15],
 	              [ 8.02004325e-09,  -5.59037685e-13,   1.60256679e-15,  -8.18749145e-19]]
            Ky = [[-9.37023860e-01,   9.89973161e-01,  -4.32284634e-06,   7.08000564e-09],
             	  [-1.59438718e-02,  -1.95099497e-05,   1.55801035e-09,   2.13743170e-13],
             	  [ 9.34251811e-06,   1.26473736e-08,  -3.71968684e-12,   5.88384784e-16],
             	  [ 3.41071678e-10,  -1.82060569e-12,   1.59690189e-15,  -3.14810895e-19]]
        elif side == 'sx':
            #dewarp_file = '/Users/Ryan/research/pipeline/dewarp_package/coeff_lmir_sx_191119.dat'
            Kx = [[  1.53484821e+00,   1.02242631e-02,  -7.62529074e-06,   3.01397583e-10],
                  [  1.01609380e+00,  -1.67032792e-05,   8.83706127e-09,  -9.33817389e-14],
                  [ -2.30147408e-05,  -1.56357707e-09,   2.07411145e-12,  -8.65727769e-16],
                  [  7.13358606e-09,   1.49417672e-12,  -1.61619128e-15,   5.24106088e-19]]
            Ky = [[  1.51620772e+01,   9.77986282e-01,  -3.28698979e-06,   6.96761931e-09],
                  [ -1.56775963e-02,  -1.41085603e-05,  -4.89279783e-09,   1.75661648e-12],
                  [  1.15874971e-05,   5.20977887e-10,   8.83909100e-12,  -2.78133098e-15],
                  [ -5.66801592e-10,   2.66518618e-12,  -3.08650185e-15,   9.23903266e-19]]
        else:
            #dewarp_file = '/Users/Ryan/research/pipeline/dewarp_package/coeff_lmir_sx_20200102.rtf'
            Kx = [[-7.96016109e+00,   9.26584096e-03,  -7.93676069e-06,   5.13414639e-10],	  
                  [ 1.02925896e+00,  -1.59974177e-05,   9.57769272e-09,  -1.14409822e-12],
 	              [-2.30169348e-05,  -3.45351550e-09,   1.89621010e-12,   6.72971750e-17],
                  [ 7.07041647e-09,   2.20511200e-12,  -1.66433082e-15,   2.55463711e-19]]
            Ky = [[-2.26409123e+00,   9.93175401e-01,  -6.67169688e-06,   8.08275391e-09],
             	  [-1.38521740e-02,  -2.27910031e-05,   4.72367613e-09,  -1.64738716e-12],
             	  [ 8.17060299e-06,   1.35240460e-08,  -5.52374318e-12,   2.14966954e-15],
             	  [ 9.25982725e-10,  -2.50607186e-12,   2.58675626e-15,  -9.82463036e-19]]

    # this is just to get the shape; image content doesn't matter
    hdul = fits.open('/Users/Ryan/research/pipeline/dewarp_package/pinholes_lmir_double_190919.fits')
    imagePinholes = hdul[0].data.copy()

    # map the coordinates that define the entire image plane
    # (transposition due to a coefficient definition change btwn Python and IDL)
    #print(np.array(Kx).T)
    dewarp_coords = dewarp.make_dewarp_coordinates(imagePinholes.shape,np.array(Kx).T,np.array(Ky).T)
    '''
    dewarp_coords = dewarp.make_dewarp_coordinates(imagePinholes.shape,
                                                   np.array(Kx).T,
                                                   np.array(Ky).T)
    '''
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        try:
            if aux_path == 'raw_data':
                imageAsterism, header = fits.getdata('{}{}/{}.fits'.format(master_path,aux_path,fname_prefix),0,header=True)
                imageAsterism_junk, header = fits.getdata('{}{}/raw_data/{}.fits'.format(master_path,test_type,fname_prefix),0,header=True)
            else:
                imageAsterism, header_junk = fits.getdata('{}reduced_data{}/{}/{}.fits'.format(master_path,test_type,aux_path,fname_prefix),0,header=True)
                imageAsterism_junk, header = fits.getdata('{}{}/raw_data/{}.fits'.format(master_path,test_type,fname_prefix),0,header=True)

            if progress_print == True:
                print('{} loaded'.format(fname_prefix))
        except:
            FileNotFoundError
            print('File Not Found Error')
            continue
        image_y = imageAsterism.shape[0]
        image_x = imageAsterism.shape[1]
        #print('XSTART = ', header['XSTART'])
        #print('XSTOP = ', header['XSTOP'])
        #print('YSTART = ', header['YSTART'])
        #print('YSTOP = ', header['YSTOP'])

        try:
            xstart = header['XSTART']
            xstop = header['XSTOP']
            ystart = header['YSTART']
            ystop = header['YSTOP']
        except:
            KeyError
            xstart = 0
            xstop = 2047
            ystart = 0
            ystop = 2047
        print('Unpadded = (y,x) {}'.format(imageAsterism.shape))

        before_1 = ((2048-ystop) + crop_bottom)-1
        after_1 = (2048 - (before_1+image_y))
        before_2 = xstart + crop_left
        after_2 = 2048-(xstart+crop_right)
        #print(before_1,after_1,before_2,after_2)
        imageAsterism_pad = np.pad(imageAsterism,((before_1,after_1),(before_2,after_2)),'constant',constant_values=0)
        #fits.writeto('{}{}_pad_test.fits'.format(return_path,fname_prefix),np.squeeze(imageAsterism_pad),header,overwrite=True)
        if imageAsterism_pad.shape[0] != 2048 or imageAsterism_pad.shape[1] != 2048:
            print('Incorrect pad dimensions!')
            break
        print('Padded = (y,x) {}'.format(imageAsterism_pad.shape))
        dewarpedAsterism = dewarp.dewarp_with_precomputed_coords(imageAsterism_pad,dewarp_coords,order=3)
        fits.writeto('{}{}.fits'.format(return_path,fname_prefix),np.squeeze(dewarpedAsterism),header,overwrite=True)
        if progress_print == True:
            print('{} dewarped'.format(fname_prefix))
    print('All frames dewarped')