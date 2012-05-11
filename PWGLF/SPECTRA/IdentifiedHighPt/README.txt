The idea is to document the different steps needed to make spectra.

STEP 1: GRID CODE
=================

In the directory grid you will find the code used to produce the trees. For
info on how to use have a look at runAAF.C

STEP 2: EXTRACT TREES
=====================

The file produced on the grid contains the trees in the list. This means that
one cannot direcrtly chain tghem. Therefore we use the code in extract_code to
make new files with the trees only.

Example:
./merge.sh aortizve/Trees_LHC10b_Pass3/files/
and 
./merge.sh aortizve/Trees_LHC10b_Pass3/files/ HighPtDeDxV0

STEP 3: COMPILE LIBRARY
=======================

cd lib
make clean
make


STEP 4: DETERMINE RATIOS
========================

This is the biggest step.

mkdir ratios_7tevb
cd ratios_7tevb

First we need to create the text files we want to analyze.
Example:
find /home/pchristi/work/analysis/7tev/ | grep HighPtDeDx_Tree | grep new | grep 117059 > 7tev_b_test.dat
find /home/pchristi/work/analysis/7tev/ | grep HighPtDeDx_Tree | grep new > 7tev_b.dat

find /home/pchristi/work/analysis/7tev/ | grep HighPtDeDxV0_Tree | grep new | grep 117059 > 7tev_b_test_v0.dat
find /home/pchristi/work/analysis/7tev/ | grep HighPtDeDxV0_Tree | grep new > 7tev_b_v0.dat

ln -s ../macros/run_code.C .
ln -s ../macros/calibrate_de_dx.C .

cp ../macros/drawText.C .
Edit the text here. This macro is used to tag the pictures.

Follow the example in the macro run_code.C

Step 1-5 is about determining the dE/dx calibrations and the input data to the
fits in pT.

Now it is time for extracting the uncorrected ratios.

ln -s ../macros/fit_yields_final.C .

This is documented in fit_yields_final.C


There is still things missing:
- Option to generate tree when generating the data.
- The code to estimate systematic errors
- Efficiency code
- Corrected fractions
- Spectra and RAA code

But that will come soon - hopefully next week.




To zip:

ls -1 README.txt lib/Makefile lib/*.cxx lib/*.h macros/*.C > tozip.txt
tar -hcvzf analysis.tgz -T tozip.txt

