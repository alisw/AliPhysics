#/bin/bash -f

export Fill=${1:-6864}

: '
# DO NOT ENABLE THIS PART: just check which argument is what
# ----------------------------------------------------------
# export TypeFlux   = ${2:-0} # 0 for FBCT, 1 for BPTX (* flux = beam intensity)
# export TypeFit    = ${3:-1} # GP2 (index 0), GP6 {1), G(2), NI (3), and DG (4)
# export SystBG     = ${4:-0} # 0 for false (default), 1 for true: enable syst. err in Create_bkgd...file_V0T0.C
# export SystODC    = ${5:-1} # 0 for Nom, 1 for ODC (default)
# export SystBBD    = ${6:-1} # 0 for false, 1 for true (default), +Q (2), -Q (3), +xi (4), and -xi (5)
# export SystOptic  = ${7:-1} # 0 for false, 1 for true (default)
# export SystPileup = ${8:-0} # 0 for false, ++ (1, +err on RatioA/C), +- (2), -+ (3), and -- (4)
# ----------------------------------------------------------
'

#: '
# Default, Intensity syst, and Fit model syst
./Execute.sh $Fill 0 1 > log_${Fill}_00_default.txt
./Execute.sh $Fill 1 1 > log_${Fill}_01_systBPTX.txt
./Execute.sh $Fill 0 3 > log_${Fill}_02_systNUM.txt

# BG corresction (on rate) syst, ODC syst, and BBD syst
./Execute.sh $Fill 0 1 1     > log_${Fill}_03_systBG.txt    # Syst. err of BG correction on rate
./Execute.sh $Fill 0 1 0 0   > log_${Fill}_04_systNoODC.txt # Disable ODC
./Execute.sh $Fill 0 1 0 1 0 > log_${Fill}_05_systNoBBD.txt # Disable BBD

# BBD syst
./Execute.sh $Fill 0 1 0 1 2 > log_${Fill}_10_systBBD+Q.txt  # BBD+Q
./Execute.sh $Fill 0 1 0 1 3 > log_${Fill}_11_systBBD-Q.txt  # BBD-Q
./Execute.sh $Fill 0 1 0 1 4 > log_${Fill}_12_systBBD+Xi.txt # BBD+Q
./Execute.sh $Fill 0 1 0 1 5 > log_${Fill}_13_systBBD-Xi.txt # BBD-Q

# Pileup syst
./Execute.sh $Fill 0 1 0 1 1 1 1 > log_${Fill}_20_systPU1.txt # ++
./Execute.sh $Fill 0 1 0 1 1 1 2 > log_${Fill}_21_systPU2.txt # +-
./Execute.sh $Fill 0 1 0 1 1 1 3 > log_${Fill}_22_systPU3.txt # -+
./Execute.sh $Fill 0 1 0 1 1 1 4 > log_${Fill}_23_systPU4.txt # --

echo ""
echo "Finished!"
mv log_$Fill*.txt ../
#'
