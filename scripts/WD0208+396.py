import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'/Users/Ryan/research/pipeline/')
import utilities as uti
import process as pro
import pre_process as pre

data_set = 'lm_161012_'
data_set_flt_drk = 'lm_161013_00'
wd = 'WD0208+396'
master_path = '/Users/Ryan/research/data/WD0208+396/'

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
path_dewarped = path_main_directory+'dewarped/'
path_ultra_center = path_main_directory+'ultra_centered/'
path_final = master_path+'final/'
path_llsg_psf_subtraction = path_final+'llsg_psf_subtraction/'
path_median_psf_subtraction = path_final+'median_psf_subtraction/'
path_median_psf = path_final+'median_psf/'

start_frame = 32126
start_frame_drk = 21
start_frame_flt = 1

end_frame = 32190
end_frame_drk = 31
end_frame_flt = 11

crop_parameter = 150
star_mag = 13.67
points = [(1300,600),(380,320),(400,950)]
width = 100

cb = 459
ct = 1816
cl = 40
cr = 1764

##cropping
#pre.edge_crop(master_path,data_set_flt_drk,'raw_data',start_frame_drk,end_frame_drk,path_drk_edge_crop,crop_bottom=cb,crop_top=ct,crop_left=cl,crop_right=cr,progress_print=True,file_type='drk')
#pre.edge_crop(master_path,data_set_flt_drk,'raw_data',start_frame_flt,end_frame_flt,path_flt_edge_crop,crop_bottom=cb,crop_top=ct,crop_left=cl,crop_right=cr,progress_print=True,file_type='flt')
#pre.edge_crop(master_path,data_set,'raw_data',start_frame,end_frame,path_sci_edge_crop,crop_bottom=cb,crop_top=ct,crop_left=cl,crop_right=cr,progress_print=True,file_type='sci')
#
##creating flat fielding frame
#pre.flat_field_mkr(master_path,data_set_flt_drk,'drk_edge_cropped','flt_edge_cropped',start_frame_drk,end_frame_drk,start_frame_flt,end_frame_flt,path_main_directory,progress_print=True)
#
##creating bp mask from flt_drk or flt_drk_median frame
#pre.bp_mask_mkr(master_path,path_main_directory,sigma_coe=3)
#
##nod sub science frames
#pro.nod_subtract(master_path,data_set,'sci_edge_cropped',start_frame,end_frame,path_subtracted,progress_print=True)
#
##bad pixel correcting flt_drk_median frame
#pre.bp_crct(master_path,'N/A','N/A','N/A','N/A',path_main_directory,progress_print=True,file_type='flt_drk_median')
#
##bad pixel correcting science frames
#pre.bp_crct(master_path,data_set,'subtracted',start_frame,end_frame,path_bp_crct,progress_print=True,file_type='sci')
#
##Testing noise change across data set
#uti.rms(master_path,data_set,'bp_crcted',start_frame,end_frame,path_main_directory,points,width,wd)
#
##flat fielding science frames 
#pre.flat_fielder(master_path,data_set,'bp_crcted',start_frame,end_frame,path_flat_fielded,progress_print=True)
#
##Cleaning sensor readout
#pro.stripe_cleaner(master_path,data_set,'flat_fielded',start_frame,end_frame,path_stripe_crct,progress_print=True,mask_size=6000) #changed mask size
#
#Dewarping science images
#pre.dewarper(master_path,data_set,'stripe_crcted',start_frame,end_frame,path_dewarped,crop_bottom=cb,crop_top=ct,crop_left=cl,crop_right=cr,date='2016B',side='',progress_print=True)
#
##Centering science frames
#pro.star_center(master_path,data_set,'dewarped',start_frame,32145,path_centered,progress_print=True,star_location=True,crop_parameter=crop_parameter)
#pro.star_center(master_path,data_set,'dewarped',32145,end_frame,path_centered,progress_print=True,star_location=True,crop_parameter=crop_parameter)
#
##Bad frame detection
#pro.frame_filter(master_path,data_set,'centered',start_frame,end_frame,path_filtered,progress_print=True,fwhm_tolerance=100,peak_tolerance=100)
#
##Ultra centering science frames
#pro.star_ultra_center(master_path,data_set,'filtered',start_frame,end_frame,path_ultra_center,progress_print=True,crop_parameter=crop_parameter)
#
##LLSG PSF subtraction method
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,frame_set='master_all_frames',rank_input=2,thresh_input=1,max_iter_input=10)
#
##Median PSF subtraction method
#pro.median_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf_subtraction,wd=wd,progress_print=True,frame_set='master_all_frames')
#
##Median combining PSF with and without rotation
#pro.psf_stacker(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,wd=wd,progress_print=True,frame_set='master_all_frames')
#
##Contrast curve for median psf subtraction method
#pro.contrast_curve_vip(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,wd=wd,star_mag=star_mag,progress_print=True,method='median_psf_subtraction',luci=False,tag='_all_frames')
#
#Contrast curve for llsg psf subtraction method
#pro.contrast_curve_vip(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,wd=wd,star_mag=star_mag,rank=20,progress_print=True,method='llsg_psf_subtraction',luci=False,tag='_20_all_frames')
#