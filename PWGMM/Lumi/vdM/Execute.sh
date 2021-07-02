#/bin/bash -f

export Fill=${1:-4937}
export TypeFlux=${2:-0} # 0 for FBCT (default), 1 for BPTX (* flux = beam intensity)
export TypeFit=${3:-1}  # GP2 (index 0), GP6 {1, default), G(2), NI (3), and DG (4)

# Systematic err options
export SystBG=${4:-0}     # 0 for false (default), 1 for true: enable syst. err in Create_bkgd_correction_file_V0T0.C
export SystODC=${5:-1}    # 0 for Nom, 1 for ODC (default)
export SystBBD=${6:-1}    # 0 for false, 1 for true (default), +Q (2), -Q (3), +xi (4), and -xi (5)
export SystOptic=${7:-1}  # 0 for false, 1 for true (default)
export SystPileup=${8:-0} # 0 for false (default), ++ (1, +err on RatioA/C), +- (2), -+ (3), and -- (4)

# QA related
export BXing=10  # Bunch Crossing # (in LHC scheme) for QA
export FForm=eps # File format of plots being printed out: ex. eps

[ $TypeFlux   -eq 0 ] && TypeFluxCorr=2 || TypeFluxCorr=1   # Exclusively for Create_intensity_correction_file.C
[ $TypeFlux   -eq 0 ] && nFlux="FBCT"   || nFlux="BPTX"     # Prefix "n" stands for name
[ $SystBG     -eq 0 ] && nBG=""         || nBG="SystBG"     # Effective from bkgd_correction_file_V0T0
[ $SystBBD    -eq 0 ] && nBBD=""        || nBBD="BBD"       # Effective from hxhy calculation part
[ $SystODC    -eq 0 ] && nSep="Nom"     || nSep="ODC"       # Effective from hxhy calculation part
[ $SystOptic  -eq 0 ] && nOptic=""      || nOptic="Optical" # Effective from hxhy calculation part
[ $SystPileup -eq 0 ] && nPileup=""     || nPileup="PU${SystPileup}"

#For BBD effect's syst. err
if [ $SystBBD -le 1 ]; then
	nBBDopt=""
elif [ $SystBBD -eq 2 ]; then
	nBBDopt="+Q"
elif [ $SystBBD -eq 3 ]; then
	nBBDopt="-Q"
elif [ $SystBBD -eq 4 ]; then
	nBBDopt="+xi"
elif [ $SystBBD -eq 5 ]; then
	nBBDopt="-xi"
else
	echo "Invalid BBD option! Stop."
	exit
fi

#-------------------------------------------------------------------------

# Check essential requirements
if	[ ! -f "Create_nominal_separation_file.C" ] ||
	[ ! -f "Create_beam_normalisation_tree.C" ] ||
	[ ! -f "Create_beam_intensity_file.C" ] ||
	[ ! -f "Create_intensity_correction_file.C" ] ||
	[ ! -f "Create_raw_rate_file.C" ] ||
	[ ! -f "Create_bkgd_correction_file_V0T0.C" ] ||
	[ ! -f "Create_bkgd_corrected_rate_file.C" ] ||
	[ ! -f "Create_pileup_corrected_rate_file.C" ] ||
	[ ! -f "Create_intensity_corrected_rate_file.C" ] ||
	[ ! -f "Create_ODC_separation_file.C" ] ||
	[ ! -f "Create_BBD_separation_file.C" ] ||
	[ ! -f "Create_optical_corrected_rate_file.C" ] ||
	[ ! -f "Create_hxhy_file.C" ] ||
	[ ! -f "Create_xs_file.C" ] ||
	[ ! -f "QA_xs.C" ] ||
	[ ! -f "QA_xs_N1N2.C" ]	; then
	echo "Cannot find one of essential requirements: stop." 
	exit
fi

export OutPath="../Fill-$Fill"
mkdir -p $OutPath
mkdir -p $OutPath/QA_intensity
mkdir -p $OutPath/QA_rate
mkdir -p $OutPath/QA_fits
mkdir -p $OutPath/QA_xs

# Link vdM input
export InPath="../Input"
echo ""
if compgen -G "$OutPath/vdm*.root" > /dev/null; then
    echo " vdM data files found: "
	ls -1 $OutPath/vdm*.root
else
	ln -s $InPath/vdm_*$Fill*.root $OutPath
	echo " Link files $InPath/vdm*$Fill*.root to $OutPath "
fi

# Link BBD/Optcial correction input - temporary for now (June 2021)
if [ $SystBBD -eq 1 ] || [ $SystOptic -eq 1 ]; then
	ln -s $InPath/Corr-$Fill-sys $OutPath
fi

echo ""
echo " - Output directory:" $OutPath
echo " - Separation:" $nSep
echo " - Intensity:" $nFlux
echo " - Fit model:" $TypeFit
echo ""
[ $SystBG -eq 1 ] && echo " - Syst. err mode on for BG rate correction"
[ $SystODC -eq 1 ] && echo " - ODC (orbit drift correction) on"
[ $SystBBD -ne 0 ] && echo " - BBD (beam-beam dynamics correction) on $nBBDopt"
[ $SystOptic -eq 1 ] && echo " - Optical distortion correction on"
[ $SystPileup -ne 0 ] && echo " - Pileup correction is being tested in mode ${SystPileup}"

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
#mv c1*.$FForm ${OutPath}/QA_intensity
#'

#Loop over rate types: V0 and T0
for TypeRate in {0..1}
do

	[ $TypeRate -eq 0 ] && nRate="VBAandVBC" || nRate="TVX"
	[ $TypeRate -eq 0 ] && nRateCorr="V0"    || nRateCorr="T0" # Create_bkgd_corrected_rate_file.C
	[ $TypeRate -eq 0 ] && TypeRateCorr=1    || TypeRateCorr=2 # Exclusively for Create_bkgd_correction_file_V0T0.C

	# 2016-2018 pileup correction factors
	if [ $Fill -eq 4937 ]; then
		[ $TypeRate -eq 0 ] && RatioAVal=0.0755 || RatioAVal=0.4459 #V0 (left), T0 (right)
		[ $TypeRate -eq 0 ] && RatioAErr=0.0002 || RatioAErr=0.0008
		[ $TypeRate -eq 0 ] && RatioCVal=0.0611 || RatioCVal=0.3911
		[ $TypeRate -eq 0 ] && RatioCErr=0.0002 || RatioCErr=0.0007
	elif [ $Fill -eq 6012 ]; then
		[ $TypeRate -eq 0 ] && RatioAVal=0.07703 || RatioAVal=0.499  
		[ $TypeRate -eq 0 ] && RatioAErr=0.00004 || RatioAErr=0.0002
		[ $TypeRate -eq 0 ] && RatioCVal=0.06216 || RatioCVal=0.3933 
		[ $TypeRate -eq 0 ] && RatioCErr=0.00004 || RatioCErr=0.0002
	elif [ $Fill -eq 6864 ]; then
		[ $TypeRate -eq 0 ] && RatioAVal=0.07684 || RatioAVal=0.49
		[ $TypeRate -eq 0 ] && RatioAErr=0.00004 || RatioAErr=0.0002
		[ $TypeRate -eq 0 ] && RatioCVal=0.06193 || RatioCVal=0.4005
		[ $TypeRate -eq 0 ] && RatioCErr=0.00004 || RatioCErr=0.0002
	else
		echo "Unknown Fill number: stop." 
		exit
	fi

	if [ $SystPileup -eq 0 ]; then
		RatioA=$RatioAVal 
		RatioC=$RatioCVal
	elif [ $SystPileup -eq 1 ]; then
		RatioA=$RatioAVal+$RatioAErr
		RatioC=$RatioCVal+$RatioCErr
	elif [ $SystPileup -eq 2 ]; then
		RatioA=$RatioAVal+$RatioAErr
		RatioC=$RatioCVal-$RatioCErr
	elif [ $SystPileup -eq 3 ]; then
		RatioA=$RatioAVal-$RatioAErr
		RatioC=$RatioCVal+$RatioCErr
	elif [ $SystPileup -eq 4 ]; then
		RatioA=$RatioAVal-$RatioAErr
		RatioC=$RatioCVal-$RatioCErr
	fi

	: '
	# Vary pileup factors
	PUVSide=1     # 0 for A, 1 for C
	PUVFactor=1.1 # Pileup factor vary by 10%
	[ $PUVSide -eq 0 ] && RatioA=$RatioAVal*$PUVFactor || RatioC=$RatioCVal*$PUVFactor
	'

	echo ""
	echo " =========================================================="
	echo " - Rate: $nRateCorr"
	echo " - Pileup factors:" $RatioA "(A)," $RatioC "(C)"
	echo " =========================================================="

	#: '
	#Rate
	$rlbq "Create_raw_rate_file.C+($Fill, \"$nRate\")"
	$rlbq "Create_bkgd_correction_file_V0T0.C+($Fill, $TypeRateCorr, $SystBG)"
	$rlbq "Create_bkgd_corrected_rate_file.C+($Fill, \"$nRate\", \"$nRateCorr$nBG\")"
	$rlbq "Create_pileup_corrected_rate_file.C+($Fill, \"$nRate\", $RatioA, $RatioC)"
	$rlbq "Create_intensity_corrected_rate_file.C+($Fill, \"$nRate\", \"$nFlux\")"

	#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"Raw\", \"Nom\", 1, $TypeFlux, $BXing)"
	#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"Raw\", \"Nom\", 2, $TypeFlux, $BXing)"
	#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"IntensityCorr$nFlux\", \"Nom\", 1, $TypeFlux, $BXing)"
	#$rlbq "QA_rate_vs_sep.C+($Fill, \"$nRate\", \"IntensityCorr$nFlux\", \"Nom\", 2, $TypeFlux, $BXing)"
	#$rlbq "QA_corr_vs_sep.C+($Fill, \"$nRate\", $TypeFlux, 1, $BXing)"
	#$rlbq "QA_corr_vs_sep.C+($Fill, \"$nRate\", $TypeFlux, 2, $BXing)"
	#mv c2*.$FForm $OutPath/QA_rate
	#'

	#: '
	# ODC, BBD, and Optical correciton
	[ $SystODC -eq 1 ] && $rlbq "Create_ODC_separation_file.C+($Fill)"
	[ $SystBBD -ne 0 ] && $rlbq "Create_BBD_separation_file.C+($Fill, \"$nBBDopt\")"
	[ $SystOptic -eq 1 ] && $rlbq "Create_optical_corrected_rate_file.C+($Fill,\"$nFlux\",\"$nRate\",\"$nBBDopt\")"

	# hxhy and cross section
	RATETYPE=${nOptic}Intensity${nBBDopt}Corr${nFlux} # This is getting too much...
	SEPTYPE=${nSep}${nBBD}${nBBDopt}

	$rlbq "Create_hxhy_file.C+($Fill, \"$nRate\", \"$RATETYPE\", \"$SEPTYPE\", $TypeFit)"
	$rlbq "Create_xs_file.C+($Fill, \"$nRate\", \"$RATETYPE\", \"$SEPTYPE\", \"$nFlux\", $TypeFit, 1,1,1,1)"

	#$rlbq "QA_xs.C+($Fill, \"$nRate\", \"$RATETYPE\", \"$SEPTYPE\", \"$nFlux\", $TypeFit, 0)"
	#$rlbq "QA_xs.C+($Fill, \"$nRate\", \"$RATETYPE\", \"$SEPTYPE\", \"$nFlux\", $TypeFit, 1)"
	$rlbq "QA_xs_N1N2.C+($Fill, \"$nRate\", \"$RATETYPE\", \"$SEPTYPE\", \"$nFlux\", $TypeFit, 0)"
	$rlbq "QA_xs_N1N2.C+($Fill, \"$nRate\", \"$RATETYPE\", \"$SEPTYPE\", \"$nFlux\", $TypeFit, 1)"
	#$rlbq "QA_headonrate.C+($Fill, \"$nRate\", \"$RATETYPE\", \"$SEPTYPE\", \"$nFlux\",$TypeFit, true)"
	#[ $TypeRate -eq 1 ] && $rlbq "QA_hxhy_V0T0.C+($Fill, \"$RATETYPE\", \"$SEPTYPE\", \"$nFlux\", $TypeFit, true)"
	mv c3*.$FForm $OutPath/QA_xs
	#'

done

#-------------------------------------------------------------------------

echo ""
echo " Done: cleaning house..."
rm *.d *.so *.pcm $OutPath/Corr-$Fill-sys $OutPath/vdm_*$Fill*.root
mv ../Fill-$Fill ../Fill${Fill}_${nSep}${nBBD}${nBBDopt}_${nFlux}_${nBG}${nOptic}${nPileup}_Fit${TypeFit}
