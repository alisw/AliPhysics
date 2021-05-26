#/bin/bash -f

export Fill=${1:-6864}
export TypeFlux=${2:-0} # 0 for FBCT, 1 for BPTX (* flux = beam intensity)
export TypeRate=${3:-0} # 0 for V0 (VBAandVBC), 1 for T0 (TVX)
export TypeFit=${4:-1}  # GP2 (index 0), GP6 {1), G(2), NI (3), and DG (4)

# Syst. check or QA
export SystBG=0    # 0 for false, 1 for true - enable systematic err option in Create_bkgd_correction_file_V0T0.C
export SystODC=0   # 0 for Nom, 1 for ODC
export SystBBD=0   # 0 for false, 1 for true
export SystOptic=0 # 0 for false, 1 for true
export BXing=10    # Bunch Crossing # (in LHC scheme) for QA

#--------------------------------------------------------------------

[ $TypeFlux -eq 0 ] && nFlux="FBCT" || nFlux="BPTX"     # prefix "n" stands for name
[ $TypeFlux -eq 0 ] && TypeFluxCorr=2 || TypeFluxCorr=1 # Exclusively for Create_intensity_correction_file.C
[ $TypeRate -eq 0 ] && nRate="VBAandVBC" || nRate="TVX"
[ $TypeRate -eq 0 ] && nRateCorr="V0" || nRateCorr="T0" #Create_bkgd_corrected_rate_file.C
[ $TypeRate -eq 0 ] && TypeRateCorr=1 || TypeRateCorr=2 # Exclusively for Create_bkgd_correction_file_V0T0.C
[ $SystODC -eq 0 ] && nSep="Nom" || nSep="ODC"          # Start application from hxhy calculation
[ $SystBBD -eq 0 ] && nBBD="" || nBBD="BBD"
[ $SystOptic -eq 0 ] && nOptic="" || nOptic="Optical"

# Pileup correction factors (V0/T0)
if [ $Fill -eq 4937 ]; then
	[ $TypeRate -eq 0 ] && RatioA=0.0755 || RatioA=0.4459 #err: 0.0002, 0.0008
	[ $TypeRate -eq 0 ] && RatioC=0.0611 || RatioC=0.3911 #err: 0.0002, 0.0007
elif [ $Fill -eq 6012 ]; then
	[ $TypeRate -eq 0 ] && RatioA=0.07703 || RatioA=0.499  #err: 0.00004, 0.000?
	[ $TypeRate -eq 0 ] && RatioC=0.06216 || RatioC=0.3933 #err: 0.00004, 0.0002
elif [ $Fill -eq 6864 ]; then
	[ $TypeRate -eq 0 ] && RatioA=0.07684 || RatioA=0.49   #err: 0.00004, 0.00?
	[ $TypeRate -eq 0 ] && RatioC=0.06193 || RatioC=0.4005 #err: 0.00004, 0.0002
else
	echo "Unknown Fill number: stop." 
	exit
fi

#--------------------------------------------------------------------

# Check requirements
if	[ ! -f "Create_nominal_separation_file.C" ] ||
	[ ! -f "Create_ODC_separation_file.C" ] ||
	[ ! -f "Create_BBD_separation_file.C" ] ||
	[ ! -f "Create_beam_normalisation_tree.C" ] ||
	[ ! -f "Create_beam_intensity_file.C" ] ||
	[ ! -f "Create_intensity_correction_file.C" ] ||
	[ ! -f "Create_bkgd_correction_file_V0T0.C" ] ||
	[ ! -f "Create_raw_rate_file.C" ] ||
	[ ! -f "Create_bkgd_corrected_rate_file.C" ] ||
	[ ! -f "Create_pileup_corrected_rate_file.C" ] ||
	[ ! -f "Create_intensity_corrected_rate_file.C" ] ||
	[ ! -f "Create_optical_corrected_rate_file.C" ] ||
	[ ! -f "Create_hxhy_file.C" ] ||
	[ ! -f "Create_xs_file.C" ] ||
	[ ! -f "QA_xs.C" ] ||
	[ ! -f "QA_xs_N1N2.C" ]	; then
	echo "Cannot find one of the requirements: stop." 
	exit
fi

export OutPath="../Fill-$Fill"
mkdir -p $OutPath
mkdir -p $OutPath/QA_intensity
mkdir -p $OutPath/QA_rate
mkdir -p $OutPath/QA_fits
mkdir -p $OutPath/QA_xs

echo ""
echo " - Output directory:" $OutPath
echo " - Separation:" $nSep
echo " - Intensity:" $nFlux
echo " - Rate:" $nRate
echo " - Fit model:" $TypeFit
echo " - Pileup factors:" $RatioA "(A)," $RatioC "(C)"
[ $SystODC -eq 1 ] && echo " - ODC on"
[ $SystBBD -eq 1 ] && echo " - BBD on"
[ $SystOptic -eq 1 ] && echo " - Optical rate correction on"

#-------------------------------------------------------------------------

export rlbq="root -l -b -q"

#: '
#Separation, Beam intensity (flux)
$rlbq "Create_nominal_separation_file.C+($Fill)"
$rlbq "Create_beam_normalisation_tree.C+($Fill, $TypeFlux)"
$rlbq "Create_beam_intensity_file.C+($Fill, $TypeFlux)"
$rlbq "Create_intensity_correction_file.C+($Fill, $TypeFluxCorr)"

#$rlbq "QA_beam_intensity.C+($Fill, $TypeFlux, 0)"
#$rlbq "QA_beam_intensity.C+($Fill, $TypeFlux, 1)"
#$rlbq "QA_normalisation_histograms.C+($Fill, $TypeFlux)"
#$rlbq "QA_get_intensities.C+($Fill, $TypeFlux, $BXing)" 
#mv c1*.png ${OutPath}/QA_intensity
#'

#: '
#Rate
$rlbq "Create_raw_rate_file.C+($Fill, \"$nRate\")"
$rlbq "Create_bkgd_correction_file_V0T0.C+($Fill, $TypeRateCorr, $SystBG)"
$rlbq "Create_bkgd_corrected_rate_file.C+($Fill, \"$nRate\", \"$nRateCorr\", $SystBG)"
$rlbq "Create_pileup_corrected_rate_file.C+($Fill, \"$nRate\", $RatioA, $RatioC, $SystBG)"
$rlbq "Create_intensity_corrected_rate_file.C+($Fill, \"$nRate\", \"$nFlux\", $SystBG)"

#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"Raw\", \"Nom\", 1, $TypeFlux, $BXing)"
#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"Raw\", \"Nom\", 2, $TypeFlux, $BXing)"
#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"IntensityCorr$nFlux\", \"Nom\", 1, $TypeFlux, $BXing)"
#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"IntensityCorr$nFlux\", \"Nom\", 2, $TypeFlux, $BXing)"
#$rlbq "QA_corr_vs_sep.C+($Fill, \"$nRate\", $TypeFlux, 1, $BXing)"
#$rlbq "QA_corr_vs_sep.C+($Fill, \"$nRate\", $TypeFlux, 2, $BXing)"
#mv c2*.png $OutPath/QA_rate
#'

#: '
#hxhy and cross section
[ $SystODC -eq 1 ] && $rlbq "Create_ODC_separation_file.C+($Fill)" #ODC file generation
[ $SystBBD -eq 1 ] && $rlbq "Create_BBD_separation_file.C+($Fill)" #BBD file generation
[ $SystOptic -eq 1 ] && $rlbq "Create_optical_corrected_rate_file.C+($Fill, \"$nRate\")" #Optical rate correction

$rlbq "Create_hxhy_file.C+($Fill, \"$nRate\", \"${nOptic}IntensityCorr$nFlux\", \"$nSep$nBBD\", $TypeFit, $SystBG)"
$rlbq "Create_xs_file.C+($Fill, \"$nRate\", \"${nOptic}IntensityCorr$nFlux\", \"$nSep$nBBD\", \"$nFlux\", $TypeFit, 1,1,1,1, $SystBG)"

$rlbq "QA_xs.C+($Fill, \"$nRate\",\"${nOptic}IntensityCorr$nFlux\", \"$nSep$nBBD\",\"$nFlux\",$TypeFit, 0, $SystBG)"
$rlbq "QA_xs.C+($Fill, \"$nRate\",\"${nOptic}IntensityCorr$nFlux\", \"$nSep$nBBD\",\"$nFlux\",$TypeFit, 1, $SystBG)"
$rlbq "QA_xs_N1N2.C+($Fill, \"$nRate\", \"${nOptic}IntensityCorr$nFlux\", \"$nSep$nBBD\", \"$nFlux\", $TypeFit, 0)"
$rlbq "QA_xs_N1N2.C+($Fill, \"$nRate\", \"${nOptic}IntensityCorr$nFlux\", \"$nSep$nBBD\", \"$nFlux\", $TypeFit, 1)"
mv c3*.png $OutPath/QA_xs
#'

