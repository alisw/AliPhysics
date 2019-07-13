
Commands for V0M in fill 7483

//-------------------------------------------------------
// separations
//-------------------------------------------------------
.L Create_nominal_separation_file_Fill_7483.C+g
Create_nominal_separation_file_Fill_7483()

//-------------------------------------------------------
// intensity
//-------------------------------------------------------

.L Create_beam_intensity_file.C+g
Create_beam_intensity_file(7483,0)
Create_beam_intensity_file(7483,1)

.L QA_beam_intensity.C+g
QA_beam_intensity(7483,0,0) // fbct, scan0
QA_beam_intensity(7483,0,1) // fbct, scan1
QA_beam_intensity(7483,1,0) // bptx, scan0
QA_beam_intensity(7483,1,1) // bptx, scan1

.L QA_normalisation_histograms.C+g
QA_normalisation_histograms(7483,0) // fbct nomalisation
QA_normalisation_histograms(7483,1) // bptx nomalisation

.L QA_get_intensities.C+g
QA_get_intensities(7483,0,10) // fbct, bc=10
QA_get_intensities(7483,1,10) // bptx, bc=10

.L Create_intensity_correction_file.C+g
Create_intensity_correction_file(7483,1) // bptx
Create_intensity_correction_file(7483,2) // fbct

//-------------------------------------------------------
// rates
//-------------------------------------------------------

.L Create_raw_rate_file.C+g
Create_raw_rate_file(7483,"V0M")

.L Create_bkgd_correction_file_V0T0.C+g
Create_bkgd_correction_file_V0T0(7483,1)

.L Create_bkgd_corrected_rate_file.C+g
Create_bkgd_corrected_rate_file(7483,"V0M","V0")

.L Create_pileup_corrected_rate_file.C+g
Create_pileup_corrected_rate_file(7483,"V0M",0,0)

.L Create_intensity_corrected_rate_file.C+g
Create_intensity_corrected_rate_file(7483,"V0M","BPTX")
Create_intensity_corrected_rate_file(7483,"V0M","FBCT")

.L QA_rate_vs_sep.C+g
QA_rate_vs_sep(7483,"V0M","Raw","Nom",1,0,0) //1=xscan,0=scan,0=bc
QA_rate_vs_sep(7483,"V0M","Raw","Nom",1,0,100)
QA_rate_vs_sep(7483,"V0M","IntensityCorrBPTX","Nom",1,0,100)

.L QA_corr_sep.C+g
QA_corr_vs_sep(7483,"V0M",0,2,0) // 0=scan,2=y scan, 0=bc

//-------------------------------------------------------
// Cross section
//-------------------------------------------------------


.L Create_hxhy_file.C+g
Create_hxhy_file(7483,"V0M","IntensityCorrBPTX","Nom",0) // 0=gauss+p2 model
Create_hxhy_file(7483,"V0M","IntensityCorrFBCT","Nom",0)

.L Create_xs_file.C+g
Create_xs_file(7483,"V0M","IntensityCorrBPTX","Nom","BPTX",0,1,1,1,1)
Create_xs_file(7483,"V0M","IntensityCorrFBCT","Nom","FBCT",0,1,1,1,1)

.L QA_xs.C+g
QA_xs(7483,"V0M","IntensityCorrFBCT","Nom","FBCT",0,0)
QA_xs(7483,"V0M","IntensityCorrFBCT","Nom","FBCT",0,1)

//-------------------------------------------------------
// Cross section ODC
//-------------------------------------------------------

.L Create_ODC_separation_file.C+g
Create_ODC_separation_file(7483)

.L Create_hxhy_file.C+g
Create_hxhy_file(7483,"V0M","IntensityCorrBPTX","ODC",0) // 0=gauss+p2 model
Create_hxhy_file(7483,"V0M","IntensityCorrFBCT","ODC",0)

.L Create_xs_file.C+g
Create_xs_file(7483,"V0M","IntensityCorrBPTX","ODC","BPTX",0,1,1,1,1)
Create_xs_file(7483,"V0M","IntensityCorrFBCT","ODC","FBCT",0,1,1,1,1)

.L QA_xs.C+g
QA_xs(7483,"V0M","IntensityCorrFBCT","ODC","FBCT",0,0)
QA_xs(7483,"V0M","IntensityCorrFBCT","ODC","FBCT",0,1)

//-------------------------------------------------------
// Numerical integration
//-------------------------------------------------------
.L Create_hxhy_file.C+g
Create_hxhy_file(7483,"V0M","IntensityCorrBPTX","NOM",3)
Create_hxhy_file(7483,"V0M","IntensityCorrBPTX","ODC",3)
Create_hxhy_file(7483,"V0M","IntensityCorrFBCT","NOM",3)
Create_hxhy_file(7483,"V0M","IntensityCorrFBCT","ODC",3)

.L Create_xs_file.C+g
Create_xs_file(7483,"V0M","IntensityCorrBPTX","NOM","BPTX",3,1,1,1,1)
Create_xs_file(7483,"V0M","IntensityCorrBPTX","ODC","BPTX",3,1,1,1,1)
Create_xs_file(7483,"V0M","IntensityCorrFBCT","NOM","FBCT",3,1,1,1,1)
Create_xs_file(7483,"V0M","IntensityCorrFBCT","ODC","FBCT",3,1,1,1,1)

.L QA_xs.C+g
QA_xs(7483,"V0M","IntensityCorrFBCT","ODC","FBCT",3,0)
QA_xs(7483,"V0M","IntensityCorrFBCT","ODC","FBCT",3,1)


//-------------------------------------------------------
// Creation of bbd files for Ivan
//-------------------------------------------------------
.L Create_bbd_file.C+g
Create_bbd_file(7483,0,"V0M","IntensityCorrFBCT","NOM","FBCT")
Create_bbd_file(7483,1,"V0M","IntensityCorrFBCT","NOM","FBCT")
Create_bbd_file(7483,0,"V0M","IntensityCorrFBCT","ODC","FBCT")
Create_bbd_file(7483,1,"V0M","IntensityCorrFBCT","ODC","FBCT")

//-------------------------------------------------------
// Use bbd files from Ivan
//-------------------------------------------------------

.L Create_BBD_separation_file.C+g
Create_BBD_separation_file(7483)

.L Create_hxhy_file.C+g
Create_hxhy_file(7483,"V0M","IntensityCorrFBCT","NomBBD",3)
Create_hxhy_file(7483,"V0M","IntensityCorrFBCT","ODCBBD",3)

.L Create_xs_file.C+g
Create_xs_file(7483,"V0M","IntensityCorrFBCT","NomBBD","FBCT",3,1,1,1,1)
Create_xs_file(7483,"V0M","IntensityCorrFBCT","ODCBBD","FBCT",3,1,1,1,1)

.L QA_xs.C+g
QA_xs(7483,"V0M","IntensityCorrFBCT","NomBBD","FBCT",3,0)
QA_xs(7483,"V0M","IntensityCorrFBCT","NomBBD","FBCT",3,1)
QA_xs(7483,"V0M","IntensityCorrFBCT","ODCBBD","FBCT",3,0)
QA_xs(7483,"V0M","IntensityCorrFBCT","ODCBBD","FBCT",3,1)
