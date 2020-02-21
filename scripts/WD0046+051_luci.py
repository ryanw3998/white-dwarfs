import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'/Users/Ryan/research/pipeline/')
import utilities as uti
import process as pro
import pre_process as pre
data_set = 'luci2.20181126.0'
data_set_flt_drk = 'luci2.20181123.0'
wd = 'WD0046+051'
master_path = '/Users/Ryan/research/data/WD0046+051_luci/'

###### PATHS FOR ROUTINE ######
path_main_directory = master_path+'reduced_data/'
path_drk_edge_crop = path_main_directory+'drk_edge_cropped/'
path_flt_edge_crop = path_main_directory+'flt_edge_cropped/'
path_sci_edge_crop = path_main_directory+'sci_edge_cropped/'
path_bp_crct = path_main_directory+'bp_crcted/'
path_flat_fielded = path_main_directory+'flat_fielded/'
path_subtracted = path_main_directory+'subtracted/'
path_centered = path_main_directory+'centered/'
path_stripe_crct = path_main_directory+'stripe_crcted/'
path_filtered = path_main_directory+'filtered/'
path_ultra_center = path_main_directory+'ultra_centered/'
path_final = master_path+'final/'
path_llsg_psf_subtraction = path_final+'llsg_psf_subtraction/'
path_median_psf_subtraction = path_final+'median_psf_subtraction/'
path_median_psf = path_final+'median_psf/'
path_raw_data = master_path+'raw_data/'


start_frame = 180
start_frame_drk = 85
start_frame_flt = 90

end_frame = 307
end_frame_drk = 90
end_frame_flt = 95

crop_parameter = 350
star_mag = 11.572
#points = [(,),(,),(,)]
width = 100

cb = 200
ct = 1850
cl = 160
cr = 1815

#De-compressing images
#uti.decompressor(master_path,data_set,'raw_data',start_frame,end_frame,path_raw_data,progress_print=True,luci=True)

##Nod appending images
#uti.nod_appender_A(master_path,data_set,'raw_data',start_frame,199,progress_print=True)
#uti.nod_appender_B(master_path,data_set,'raw_data',199,219,progress_print=True)
#uti.nod_appender_A(master_path,data_set,'raw_data',219,239,progress_print=True)
#uti.nod_appender_B(master_path,data_set,'raw_data',239,259,progress_print=True)
#uti.nod_appender_A(master_path,data_set,'raw_data',259,279,progress_print=True)
#uti.nod_appender_B(master_path,data_set,'raw_data',279,299,progress_print=True)
#uti.nod_appender_A(master_path,data_set,'raw_data',299,end_frame,progress_print=True)

##cropping
#pre.edge_crop(master_path,data_set_flt_drk,'raw_data',start_frame_drk,end_frame_drk,path_drk_edge_crop,crop_bottom=cb,crop_top=ct,crop_left=cl,crop_right=cr,progress_print=True,file_type='drk')
#pre.edge_crop(master_path,data_set_flt_drk,'raw_data',start_frame_flt,end_frame_flt,path_flt_edge_crop,crop_bottom=cb,crop_top=ct,crop_left=cl,crop_right=cr,progress_print=True,file_type='flt')
#pre.edge_crop(master_path,data_set,'raw_data',start_frame,end_frame,path_sci_edge_crop,crop_bottom=cb,crop_top=ct,crop_left=cl,crop_right=cr,progress_print=True,file_type='sci')
#
##creating flat fielding frame
#pre.flat_field_mkr(master_path,data_set_flt_drk,'drk_edge_cropped','flt_edge_cropped',start_frame_drk,end_frame_drk,start_frame_flt,end_frame_flt,path_main_directory,progress_print=True)
#
##creating bp mask
#pre.bp_mask_mkr(master_path,path_main_directory,sigma_coe=3)
#
##nod sub science frames
#pro.nod_subtract(master_path,data_set,'sci_edge_cropped',start_frame,299,path_subtracted,progress_print=True)
#pro.nod_subtract(master_path,data_set,'sci_edge_cropped',291,end_frame,path_subtracted,progress_print=True)
#
##bad pixel correcting flt_drk_median frame
#pre.bp_crct(master_path,'N/A','N/A','N/A','N/A',path_main_directory,progress_print=True,file_type='flt_drk_median')
#
##bad pixel correcting science frames
#pre.bp_crct(master_path,data_set,'subtracted',start_frame,end_frame,path_bp_crct,progress_print=True,file_type='sci')
##
###Testing RMS of data set
##uti.rms(master_path,data_set,'bp_crcted',start_frame,end_frame,path_main_directory,points,width,wd,luci=True)
##
##Flat fielding induces more bad pixels. There is also no obvious vingetting so the are ultimatley unessicary
##flat fielding science frames 
#pre.flat_fielder(master_path,data_set,'bp_crcted',start_frame,end_frame,path_flat_fielded,progress_print=True)
#
##Cleaning sensor readout
#pro.stripe_cleaner(master_path,data_set,'flat_fielded',start_frame,end_frame,path_stripe_crct,progress_print=True,mask_size=6000) #changed mask size
#
##Centering science frames
#pro.star_center(master_path,data_set,'stripe_crcted',start_frame,end_frame,path_centered,progress_print=True,star_location=True,crop_parameter=crop_parameter)
#
##Bad frame detection
#pro.frame_filter(master_path,data_set,'centered',start_frame,end_frame,path_filtered,progress_print=True,fwhm_tolerance=95,peak_tolerance=95)

#Ultra centering science frames
pro.star_ultra_center(master_path,data_set,'centered',start_frame,end_frame,path_ultra_center,progress_print=True,crop_parameter=crop_parameter,luci=True)

#LLSG PSF subtraction method
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,frame_set='master_all_frames',rank_input=2,thresh_input=1,max_iter_input=10,luci=True)

#Median PSF subtraction method
pro.median_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf_subtraction,progress_print=True,frame_set='master_all_frames',luci=True)

#Median combining PSF with and without rotation
pro.psf_stacker(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,wd=wd,progress_print=True,frame_set='master_all_frames',luci=True)

#Contrast curve for median psf subtraction method
pro.contrast_curve_vip(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,progress_print=True,method='median_psf_subtraction',luci=True,tag='_all_frames')

#Contrast curve for llsg psf subtraction method
pro.contrast_curve_vip(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,method='llsg_psf_subtraction',luci=True,tag='_all_frames')
