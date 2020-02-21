import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'/Users/Ryan/research/pipeline/')
import utilities as uti
import process as pro
import pre_process as pre

data_set = 'luci2.20181126.0'
data_set_flt_drk = 'luci2.20181123.0'
master_path = '/Users/Ryan/research/data/WD2336-079_luci/'

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

start_frame = 12
start_frame_drk = 85
start_frame_flt = 90

end_frame = 172
end_frame_drk = 90
end_frame_flt = 95

crop_parameter = 330
star_mag = 13.341
#image_center = ***
#plot_radius = ***

cb = 200
ct = 1850
cl = 160
cr = 1815

#De-compressing images
#uti.decompressor(master_path,data_set,'raw_data',start_frame,end_frame,path_raw_data,progress_print=True,luci=True)

#Nod appending images
#uti.nod_appender_A(master_path,data_set,'raw_data',32,33,progress_print=True)
#uti.nod_appender_B(master_path,data_set,'raw_data',33,52,progress_print=True)
#uti.nod_appender_A(master_path,data_set,'raw_data',52,72,progress_print=True)
#uti.nod_appender_B(master_path,data_set,'raw_data',72,92,progress_print=True)
#uti.nod_appender_A(master_path,data_set,'raw_data',92,112,progress_print=True)
#uti.nod_appender_B(master_path,data_set,'raw_data',112,132,progress_print=True)
#uti.nod_appender_A(master_path,data_set,'raw_data',132,152,progress_print=True)
#uti.nod_appender_B(master_path,data_set,'raw_data',152,172,progress_print=True)

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
#nod sub science frames
#pro.nod_subtract(master_path,data_set,'sci_edge_cropped',start_frame,end_frame,path_subtracted,progress_print=True)
#
##bad pixel correcting flt_drk_median frame
#pre.bp_crct(master_path,'N/A','N/A','N/A','N/A',path_main_directory,progress_print=True,file_type='flt_drk_median')
#
##bad pixel correcting science frames
#pre.bp_crct(master_path,data_set,'subtracted',start_frame,end_frame,path_bp_crct,progress_print=True,file_type='sci')
#
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
#pro.frame_filter(master_path,data_set,'centered',start_frame,end_frame,path_filtered,progress_print=True,fwhm_tolerance=85,peak_tolerance=85)
#
##Ultra centering science frames
#pro.star_ultra_center(master_path,data_set,'filtered',start_frame,end_frame,path_ultra_center,progress_print=True,crop_parameter=crop_parameter)
#
#LLSG PSF subtraction method
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank3',rank_input=3,thresh_input=1,max_iter_input=10)
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank4',rank_input=4,thresh_input=1,max_iter_input=10)
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank5',rank_input=5,thresh_input=1,max_iter_input=10)
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank6',rank_input=6,thresh_input=1,max_iter_input=10)
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank7',rank_input=7,thresh_input=1,max_iter_input=10)
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank8',rank_input=8,thresh_input=1,max_iter_input=10)
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank9',rank_input=9,thresh_input=1,max_iter_input=10)
#pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank10',rank_input=10,thresh_input=1,max_iter_input=10)

#Apply this to kelt 4 data to see how different terms change over time. Also find out how to extract direct tables/arrays from contrast curve data
#is "Student" discussed in Mawet paper?
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank11',rank_input=11,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank12',rank_input=12,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank13',rank_input=13,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank14',rank_input=14,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank15',rank_input=15,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank16',rank_input=16,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank17',rank_input=17,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank18',rank_input=18,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank19',rank_input=19,thresh_input=1,max_iter_input=10)
pro.llsg_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_llsg_psf_subtraction,progress_print=True,luci=True,frame_set='master_rank20',rank_input=20,thresh_input=1,max_iter_input=10)


#pro.contrast_curve_test(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,progress_print=True,frame_set='master',luci=True,rank=2)
#
##Median PSF subtraction method
#pro.median_psf_subtraction(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf_subtraction,progress_print=True,frame_set='master',luci=True)
#
##Median combining PSF with and without rotation
#pro.psf_stacker(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,progress_print=True,frame_set='master',luci=True)
#
##plotting
#xpos_array,counts_array_md,flux_array_md,mag_array_md=uti.contrast_curves(master_path,'master_deroto',path_median_psf,image_center,star_mag,'g',plot_radius,find_star_counts=True)
#plt.legend(('master_deroto','median_sub_master','llsg_sub_master'),loc='upper right')
#
#xpos_array,counts_array_msm,flux_array_msm,mag_array_msm=uti.contrast_curves(master_path,'median_sub_master',path_median_psf_subtraction,image_center,star_mag,'b',plot_radius)
#plt.legend(('master_deroto','median_sub_master','llsg_sub_master'),loc='upper right')
#
##xpos_array,counts_array_lsm,flux_array_lsm,mag_array_lsm=uti.contrast_curves(master_path,'llsg_sub_master',path_llsg_psf_subtraction,image_center,star_mag,'r',plot_radius)
##plt.legend(('master_deroto','median_sub_master','llsg_sub_master'),loc='upper right')
#
#uti.multiplot(xpos_array,counts_array_md,counts_array_msm,counts_array_lsm,flux_array_md,flux_array_msm,flux_array_lsm,mag_array_md,mag_array_msm,mag_array_lsm,path_final)

#pro.contrast_curve_test(master_path,data_set,'ultra_centered',start_frame,end_frame,path_median_psf,progress_print=True,frame_set='master',luci=True)
