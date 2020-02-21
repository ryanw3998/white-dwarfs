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

def data_opener(master_path,prefix,aux_path,start_frame_num,end_frame_num,test_type='',progress_print=False,return_parang=False,luci=False):
    #Function for opening multiple fits files at once

    fnames = []
    fname_prefixes = []
    parangs = []
    if luci == True:
        parang_header = 'POSANGLE'
    else:
        parang_header = 'LBT_PARA'
    if return_parang == True:
        try:
            parangs_load = np.loadtxt('{}reduced_data{}/parangs.txt'.format(master_path,test_type))
        except:
            FileNotFoundError
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        #print(fname_prefix)
        try:
            if aux_path == 'raw_data':
                fname = fits.open('{}{}/{}.fits'.format(master_path,aux_path,fname_prefix))
                fname_float_data = fname[0].data.astype(float)
            else:
                fname = fits.open('{}reduced_data{}/{}/{}.fits'.format(master_path,test_type,aux_path,fname_prefix))
                fname_float_data = fname[0].data.astype(float)
                #print(fname_float_data.shape)
            if progress_print == True:
                print('{} loaded'.format(fname_prefix))
            if return_parang == True:
                try:
                    parang = fits.open('{}raw_data/{}.fits'.format(master_path,fname_prefix))[0].header['{}'.format(parang_header)]
                    parangs.append(parang)
                except:
                    parangs.append(parangs_load[i])
                    print('No parangs in header')
            fnames.append(fname_float_data)
            fname_prefixes.append(fname_prefix)
            fname.close()
        except:
            FileNotFoundError
            print('File Not Found Error')
    if return_parang == True:
        return fnames,fname_prefixes,parangs
    else:
        return fnames,fname_prefixes

def nod_appender_A(master_path,prefix,aux_path,start_frame_num,end_frame_num,progress_print,luci=False):
    #Function for writing information to fits headers that is essential to the pipeline
    # This information is not included with files from the LUCI imager on LBT

    if luci == True:
        fname_end = '.fits.gz'
    else:
        fname_end = '.fits'
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        try:
            hdulist = fits.open('{}raw_data/{}{}'.format(master_path,fname_prefix,fname_end))
        except:
            FileNotFoundError
            print('File Not Found Error')
            continue
        header = hdulist[0].header
        header.set('FLAG','NOD_A')
        hdulist.writeto('{}raw_data/{}{}'.format(master_path,fname_prefix,fname_end),overwrite=True)
        if progress_print == True:
            print('{} NOD_A appended'.format(fname_prefix))
        hdulist.close()
        print(fits.open('{}raw_data/{}{}'.format(master_path,fname_prefix,fname_end))[0].header['FLAG'])

def nod_appender_B(master_path,prefix,aux_path,start_frame_num,end_frame_num,progress_print,luci=False):
    #Function for writing information to fits headers that is essential to the pipeline
    # This information is not included with files from the LUCI imager on LBT

    if luci == True:
        fname_end = '.fits.gz'
    else:
        fname_end = '.fits'
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        try:
            hdulist = fits.open('{}raw_data/{}{}'.format(master_path,fname_prefix,fname_end))
        except:
            FileNotFoundError
            print('File Not Found Error')
            continue
        header = hdulist[0].header
        header.set('FLAG','NOD_B')
        hdulist.writeto('{}raw_data/{}{}'.format(master_path,fname_prefix,fname_end),overwrite=True)
        if progress_print == True:
            print('{} NOD_B appended'.format(fname_prefix))
        hdulist.close()
        print(fits.open('{}raw_data/{}{}'.format(master_path,fname_prefix,fname_end))[0].header['FLAG'])

def parang_calculator(master_path,prefix,aux_path,start_frame_num,end_frame_num,progress_print,luci=False):
    #calculates parallactic angles for image derotation

    if luci == True:
        dec_kywrd = 'TELDEC'
        ra_kywrd = 'TELRA'
        date_obs_kywrd = 'MJD-OBS'
        parang_header = 'POSANGLE'
    else:
        dec_kywrd = 'LBT_DEC'
        ra_kywrd = 'LBT_RA'
        date_obs_kywrd = 'LBT_LST'
        acqtime_kywrd = 'EXPTIME'
        parang_header = 'LBT_PARA'
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        try:
            hdulist = fits.open('{}raw_data/{}.fits'.format(master_path,fname_prefix))
        except:
            FileNotFoundError
            print('File Not Found Error')
            continue
        header = hdulist[0].header
        #parang_true = fits.open('{}raw_data/{}.fits'.format(master_path,fname_prefix))[0].header['{}'.format(parang_header)]
        parang = np.round(vip.preproc.compute_paral_angles(header,32.70131,ra_kywrd,dec_kywrd,date_obs_kywrd,acqtime_kywrd),5)
        header.set('LBT_PARA',parang)
        hdulist.writeto('{}raw_data/{}.fits'.format(master_path,fname_prefix),overwrite=True)
        hdulist.close()
    print('All parangs calculated')

def contrast_curves(master_path,file,return_path,image_center,mag_star,graph_color,x_lim,test_type='',find_star_counts=False,r_o=2.7,sigma=5,luci=False,save_table=False): #table_type,
    #Contrast curve image quality metric

    print('Creating Contrast Curves')
    d_o = 2*r_o
    frame = fits.open('{}{}.fits'.format(return_path,file))[0].data.astype(float)
    fwhm_mean = np.loadtxt('{}reduced_data{}/fwhm_median.txt'.format(master_path,test_type))
    fwhm_mean = float(fwhm_mean)
    #print(fwhm_mean)
    quality = fwhm_mean/d_o #comparing measured fwhm to theoretical fwhm
    print('Meas FWHM/Theo FWHM', quality)
    if find_star_counts == True: #finding value of star to divide each radii by
        try:
            os.remove('{}reduced_data{}/star_counts.p'.format(master_path,test_type))
        except:
            FileNotFoundError
        star_ap = ca((image_center,image_center),r_o)
        star_qtable = ap(frame,star_ap)
        star_counts = star_qtable['aperture_sum'].quantity
        np.savetxt('{}reduced_data{}/star_counts.txt'.format(master_path,test_type),star_counts)
        #pickle.dump(star_counts, open('{}reduced_data{}/star_counts.p'.format(master_path,test_type),'wb'))
    else:
        star_counts = np.loadtxt('{}reduced_data{}/star_counts.txt'.format(master_path,test_type))
        #star_counts = pickle.load(open('{}reduced_data{}/star_counts.p'.format(master_path,test_type),'rb'))
    #print(fwhm_mean)
    x,y = np.shape(frame)
    #print(int(x),int(y))
    #x = 2*x_lim+1
    #numrad = np.round(x/(2*fwhm_mean)) #number of concentric circles from central star
    #print(numrad)
    numrad = np.round(x/(2*d_o))-1
    #print(numrad)
    xpos_array = []
    counts_array = []
    flux_array = []
    mag_array = []
    for i in range(1,int(numrad)): #for each circle, creating ring of circle aperatures
        #ap_rad = np.round(i * fwhm_mean,2)
        ap_rad = np.round(i * d_o,2)
        #xpos_array.append(ap_rad) #appending circle's radius for graph
        num_ap = np.round(2*np.pi*ap_rad/d_o)
        deg_ap = np.pi/num_ap#no 2* for rotational overlap. Also in radians
        x_array = []
        y_array = []
        for j in range(int(2*num_ap)):#generating coordinates for each annuli
            x = np.round(ap_rad*np.cos(j*deg_ap),5)+image_center 
            y = np.round(ap_rad*np.sin(j*deg_ap),5)+image_center
            x_array.append(x)
            y_array.append(y)
        positions = list(zip(x_array,y_array))
        apertures = ca(positions,r_o)
        qtable = ap(frame,apertures)#generating photometry
        std = np.std(qtable['aperture_sum'].quantity)*sigma 
        counts_array.append(std)
        std_star_crcted = std/star_counts
        flux_array.append(std_star_crcted)
        mag = mag_star - 2.5*math.log10(std_star_crcted)
        xpos_array.append(ap_rad) #appending circle's radius for graph
        mag_array.append(mag)
    print('Done')
    del xpos_array[i-1]
    del counts_array[i-1]
    del flux_array[i-1]
    del mag_array[i-1]


    #if save_table == True:
    #    data = np.array([np.around(xpos_array,2),mag_array])
    #    data = data.T
    #    np.savetxt('{}magvsradius_cleaned.txt'.format(return_path),data,header='Mag vs radii for WD0208+396 (cleaned) \nLeft row is the radii and right row is the mag at that radii')

    #if save_table == True:
    #    if table_type == 'mag':
    #        array = mag_array
    #    elif table_type == 'flux':
    #        array = flux_array
    #    elif table_type == 'counts':
    #        array = counts_array
    #    data = np.array([np.around(xpos_array,2),array])
    #    data = data.T
    #    np.savetxt('{}{}vsradius_{}{}.txt'.format(return_path,table_type,file,test_type),data)
    plt.figure(1)
    plt.title('Counts vs radius for {}'.format(file))
    plt.xlabel('Radii (pixels)')
    plt.ylabel('Counts')
    plt.xlim(0,x_lim)
    counts_plot = plt.plot(xpos_array,counts_array)
    plt.savefig('{}countsvsradius_{}.png'.format(return_path,file))

    plt.figure(2)
    plt.title('Flux vs radius for {}'.format(file))
    plt.xlabel('Radii (pixels)')
    plt.ylabel('Flux')
    plt.xlim(0,x_lim)
    flux_plot = plt.plot(xpos_array,flux_array)
    plt.savefig('{}fluxvsradius_{}.png'.format(return_path,file))

    plt.figure(3)
    plt.title('Magnitude vs radius for {}'.format(file))
    plt.xlabel('Radii (pixels)')
    plt.ylabel('Magnitude')
    plt.xlim(0,x_lim)
    plt.gca().invert_yaxis()
    mag_plot = plt.plot(xpos_array,mag_array)
    plt.savefig('{}magvsradius_{}.png'.format(return_path,file))
    
    plt.show()

    return(xpos_array,counts_array,flux_array,mag_array)
    #return(counts_plot,flux_plot,mag_plot) #plotting 

def multiplot(xpos_array,array1,array2,array3,array4,array5,array6,array7,array8,array9,path_final):
    #For displaying multiple plots at once, not used anymore

    plt.figure(1)
    plt.plot(xpos_array,array1)
    plt.plot(xpos_array,array2)
    plt.plot(xpos_array,array3)

    plt.legend(('master_deroto','median_sub_master','llsg_sub_master'),loc='upper right')
    plt.title('Counts vs radius')
    plt.xlabel('Radii (pixels)')
    plt.ylabel('Counts')
    plt.savefig('{}countsvsradius.png'.format(path_final))

    plt.figure(2)
    plt.plot(xpos_array,array4)
    plt.plot(xpos_array,array5)
    plt.plot(xpos_array,array6)

    plt.legend(('master_deroto','median_sub_master','llsg_sub_master'),loc='upper right')
    plt.title('Flux vs radius')
    plt.xlabel('Radii (pixels)')
    plt.ylabel('Flux')
    plt.savefig('{}fluxvsradius.png'.format(path_final))

    plt.figure(3)
    plt.plot(xpos_array,array7)
    plt.plot(xpos_array,array8)
    plt.plot(xpos_array,array9)

    plt.legend(('master_deroto','median_sub_master','llsg_sub_master'),loc='upper right')
    plt.title('Magnitude vs radius')
    plt.xlabel('Radii (pixels)')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    plt.savefig('{}magvsradius.png'.format(path_final))

    plt.show()

def decompressor(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,progress_print=False,luci=False):
    #Decompresses fit.gz images

    if luci == True:
        fname_end = '.fits.gz'
    else:
        fname_end = '.fits'
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        try:
            if aux_path == 'raw_data':
                fname = fits.open('{}{}/{}{}'.format(master_path,aux_path,fname_prefix,fname_end))
            else:
                fname = fits.open('{}reduced_data/{}/{}{}'.format(master_path,aux_path,fname_prefix,fname_end))
            if progress_print == True:
                print('{} loaded'.format(fname_prefix))
        except:
            FileNotFoundError
            print('File Not Found Error')
        data = fname[0].data
        header = fname[0].header
        fits.writeto('{}{}.fits'.format(return_path,fname_prefix),data=data, header=header, overwrite = True)
        if progress_print == True:
            print('{} uncompressed'.format(fname_prefix))
        fname.close()
    print('All frames uncompressed')


def biassub(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,test_type='',progress_print=False,luci=False):
    for i in range(start_frame_num, end_frame_num, 1):
        fname_prefix = '%s%03d'%(prefix,i)
        print(fname_prefix)
        try:
            if aux_path == 'raw_data':
                fname = fits.open('{}{}/{}.fits'.format(master_path,aux_path,fname_prefix))
                fname_float = fname[0].data.astype(float)
            elif aux_path == 'raw_data_bias':
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
        #print(type(fname))
        #print(fname_float.shape)
        fname_bias_sub = fname_float[1] - fname_float[0]
        #fname_bias_sub = fname_float[1]
        header = fname[0].header
        fits.writeto('{}{}.fits'.format(return_path,fname_prefix),data=fname_bias_sub, header=header, overwrite = True)
        print(fname_bias_sub.shape)

def rms(master_path,prefix,aux_path,start_frame_num,end_frame_num,return_path,points,width,wd,tag='',test_type='',progress_print=False,luci=False):
    all_rms = []
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
            all_rms.append([0,0,0])
            continue
        rms_list = []
        print('Frame =', i)
        for j in range(3):
            cb = points[j][1] - width
            ct = points[j][1] + width
            cl = points[j][0] - width
            cr = points[j][0] + width
            #print(cb,ct,cl,cr)
            section = fname_float[cb:ct,cl:cr]
            #print('Frame = ', i)
            #print('Median = ', np.median(section))
            #print('Mean = ', np.mean(section) )
            #print('Standard deviation =', np.std(section))
            bp_list_pos = np.where(section >= np.median(section)+5*np.std(section))
            bp_list_neg = np.where(section <= np.median(section)-5*np.std(section))
            bp_list_posx = list(bp_list_pos[0])
            bp_list_posy = list(bp_list_pos[1])
            bp_list_negx = list(bp_list_neg[0])
            bp_list_negy = list(bp_list_neg[1])
            bp_list_totx = bp_list_posx + bp_list_negx
            bp_list_toty = bp_list_posy + bp_list_negy
            bp_list = [bp_list_totx,bp_list_toty]
            #print('Number of bad pixels = ', len(bp_list[0]))
            section[tuple(bp_list)] = 0
            #fits.writeto('section_test{}.fits'.format(j),np.squeeze(section),overwrite=True)
            #fits.writeto('{}section_upperleft_post5_7124.fits'.format(return_path),data=section, overwrite = True)
            array_size = (section.shape[0])*(section.shape[1])-len(bp_list[0])
            rms = np.sqrt(np.sum(np.square(section))/(array_size))
            rms_list.append(rms)
            print('Post RMS = ',rms)
            #print()
        all_rms.append(rms_list)
    all_rms = np.array (all_rms)
    #print(all_rms[:,0])
    plt.plot(np.arange(start_frame_num,end_frame_num),all_rms[:,0],label=points[0])
    plt.plot(np.arange(start_frame_num,end_frame_num),all_rms[:,1],label=points[1])
    plt.plot(np.arange(start_frame_num,end_frame_num),all_rms[:,2],label=points[2])
    plt.legend(fontsize=10)
    plt.title('RMS per frame {}'.format(wd))
    plt.xlabel('Frame Number ({}+)'.format(prefix))
    plt.ylabel('RMS')
    plt.savefig('{}rms_{}{}.png'.format(return_path,wd,tag))
    plt.show()