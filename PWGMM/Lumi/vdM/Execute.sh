#/bin/bash -f

export Fill=${1:-6864}
export BXing=${2:-10} # Bunch Crossing for QA

# Options: intensity, rate, fit model
export TypeFlux=0 # 0 for FBCT, 1 for BPTX (* flux = beam intensity)
export TypeRate=0 # 0 for V0 (VBAandVBC), 1 for T0 (TVX)
export TypeFit=1 # GP2 (index 0), GP6 {1), G(2), NI (3), and DG (4)
export SystChk=0 # 0 for false, 1 for true - enable systematic option from Create_bkgd_correction_file_V0T0.C
export SystODC=0 # 0 for Nom, 1 for ODC

[ $TypeFlux -eq 0 ] && NameFlux="FBCT" || NameFlux="BPTX"
[ $TypeFlux -eq 0 ] && TypeFluxCorr=2 || TypeFluxCorr=1 # Exclusively for Create_intensity_correction_file.C
[ $TypeRate -eq 0 ] && NameRate="VBAandVBC" || NameRate="TVX"
[ $TypeRate -eq 0 ] && NameRateCorr="V0" || NameRateCorr="T0" #Create_bkgd_corrected_rate_file.C
[ $TypeRate -eq 0 ] && TypeRateCorr=1 || TypeRateCorr=2 # Exclusively for Create_bkgd_correction_file_V0T0.C
[ $SystODC -eq 0 ] && NameSep="Nom" || NameSep="ODC" # Start application from hxhy calculation

# Pileup correction factors (V0/T0)
if [ $Fill -eq 4937 ]; then
	[ $TypeRate -eq 0 ] && RatioA=0.0755 || RatioA=0.4459 #err: 0.0002, 0.0008
	[ $TypeRate -eq 0 ] && RatioC=0.0611 || RatioC=0.3911 #err: 0.0002, 0.0007
elif [ $Fill -eq 6012 ]; then
	[ $TypeRate -eq 0 ] && RatioA=0.0755 || RatioA=0.4459 #err: 0.0002, 0.0008
	[ $TypeRate -eq 0 ] && RatioC=0.0611 || RatioC=0.3911 #err: 0.0002, 0.0007
elif [ $Fill -eq 6864 ]; then
	[ $TypeRate -eq 0 ] && RatioA=0.07703 || RatioA=0.4990 #err: 0.00004, 0.0002
	[ $TypeRate -eq 0 ] && RatioC=0.06216 || RatioC=0.3933 #err: 0.00004, 0.0002
else
	echo "Unknown Fill number: stop." 
	exit
fi

# Check requirements
if	[ ! -f "Checks.C" ] ||
	[ ! -f "Create_nominal_separation_file.C" ] || #!
	[ ! -f "Create_ODC_separation_file.C" ] ||
	[ ! -f "Create_beam_normalisation_tree.C" ] ||
	[ ! -f "Create_beam_intensity_file.C" ] ||
	[ ! -f "Create_intensity_correction_file.C" ] ||
	[ ! -f "QA_beam_intensity.C" ] ||
	[ ! -f "QA_normalisation_histograms.C" ] ||
	[ ! -f "QA_get_intensities.C" ] ||	#!
	[ ! -f "Create_bkgd_correction_file_V0T0.C" ] ||
	[ ! -f "Create_raw_rate_file.C" ] ||
	[ ! -f "Create_bkgd_corrected_rate_file.C" ] ||
	[ ! -f "Create_pileup_corrected_rate_file.C" ] ||
	[ ! -f "Create_intensity_corrected_rate_file.C" ] ||
	[ ! -f "QA_rate_vs_sep.C" ] ||
	[ ! -f "QA_corr_vs_sep.C" ]	|| #!
	[ ! -f "Create_hxhy_file.C" ] ||
	[ ! -f "Create_xs_file.C" ] ||
	[ ! -f "QA_xs.C" ] ||
	[ ! -f "QA_xs_N1N2.C" ]	; then
	echo "Cannot find one of the requirements: stop." 
	exit
fi

export OutPath="../Fill-$Fill"
mkdir -p $OutPath
mkdir -p $OutPath/Fits
mkdir -p $OutPath/QA_intensity
mkdir -p $OutPath/QA_rate
mkdir -p $OutPath/QA_xs

echo ""
#echo "Procced with following setup: "
echo " - Output directory:" $OutPath
echo " - Separation (from hxhy):" $NameSep
echo " - Intensity type:" $NameFlux
echo " - Rate type:" $NameRate
echo " - Pileup factors:" $RatioA"," $RatioC

#-------------------------------------------------------------------------

export rlbq="root -l -b -q"

: '
#Separation, Beam intensity (flux)
#$rlbq "Checks.C+($Fill, 0)"
#$rlbq "Checks.C+($Fill, 1)"
$rlbq "Create_nominal_separation_file.C+($Fill)"
#$rlbq "Create_ODC_separation_file.C+($Fill)"
$rlbq "Create_beam_normalisation_tree.C+($Fill, $TypeFlux)"
$rlbq "Create_beam_intensity_file.C+($Fill, $TypeFlux)"
$rlbq "Create_intensity_correction_file.C+($Fill, $TypeFluxCorr)"

$rlbq "QA_beam_intensity.C+($Fill, $TypeFlux, 0)"
$rlbq "QA_beam_intensity.C+($Fill, $TypeFlux, 1)"
$rlbq "QA_normalisation_histograms.C+($Fill, $TypeFlux)"
$rlbq "QA_get_intensities.C+($Fill, $TypeFlux, $BXing)" 
mv c1*.png ${OutPath}/QA_intensity
'

: '
#Rate
$rlbq "Create_raw_rate_file.C+($Fill, \"$NameRate\")"
$rlbq "Create_bkgd_correction_file_V0T0.C+($Fill, $TypeRateCorr, $SystChk)"
$rlbq "Create_bkgd_corrected_rate_file.C+($Fill, \"$NameRate\", \"$NameRateCorr\", $SystChk)"
$rlbq "Create_pileup_corrected_rate_file.C+($Fill, \"$NameRate\", $RatioA, $RatioC, $SystChk)"
$rlbq "Create_intensity_corrected_rate_file.C+($Fill, \"$NameRate\", \"$NameFlux\", $SystChk)"

$rlbq "QA_rate_vs_sep.C+($Fill, \"$NameRate\", \"Raw\", \"Nom\", 1, $TypeFlux, $BXing)"
$rlbq "QA_rate_vs_sep.C+($Fill, \"$NameRate\", \"Raw\", \"Nom\", 2, $TypeFlux, $BXing)"
$rlbq "QA_rate_vs_sep.C+($Fill, \"$NameRate\", \"IntensityCorr$NameFlux\", \"Nom\", 1, $TypeFlux, $BXing)"
$rlbq "QA_rate_vs_sep.C+($Fill, \"$NameRate\", \"IntensityCorr$NameFlux\", \"Nom\", 2, $TypeFlux, $BXing)"
$rlbq "QA_corr_vs_sep.C+($Fill, \"$NameRate\", $TypeFlux, 1, $BXing)"
$rlbq "QA_corr_vs_sep.C+($Fill, \"$NameRate\", $TypeFlux, 2, $BXing)"
mv c2*.png $OutPath/QA_rate
'

: '
#hxhy and cross section
$rlbq "Create_hxhy_file.C+($Fill, \"$NameRate\", \"IntensityCorr$NameFlux\", \"$NameSep\", $TypeFit, $SystChk)"
$rlbq "Create_xs_file.C+($Fill, \"$NameRate\",\"IntensityCorr$NameFlux\", \"$NameSep\", \"$NameFlux\", $TypeFit, 1,1,1,1, $SystChk)"

$rlbq "QA_xs.C+($Fill, \"$NameRate\", \"IntensityCorr$NameFlux\", \"$NameSep\", \"$NameFlux\", $TypeFit, 0, $SystChk)"
$rlbq "QA_xs.C+($Fill, \"$NameRate\", \"IntensityCorr$NameFlux\", \"$NameSep\", \"$NameFlux\", $TypeFit, 1, $SystChk)"
$rlbq "QA_xs_N1N2.C+($Fill, \"$NameRate\", \"IntensityCorr$NameFlux\", \"$NameSep\", \"$NameFlux\", $TypeFit, 0)"
$rlbq "QA_xs_N1N2.C+($Fill, \"$NameRate\", \"IntensityCorr$NameFlux\", \"$NameSep\", \"$NameFlux\", $TypeFit, 1)"
mv c3*.png $OutPath/QA_xs
'
