****************************************************************
* Instructions for running multiplicity distribution analysis  *
*                    ------------                              *
* Casper Nygaard, nygaard.casper@gmail.com, July 21st, 2012    *
****************************************************************

In these instructions the abbrevation $FMD corresponds to the
$ALICE_ROOT/PWGLF/FORWARD/analysis2/ directory.

For running the multiplicity distribution analysis the following files
are used (in addition to some general AliRoot FORWARD code):

$FMD/MakeMultAOD.C
$FMD/MakeMultiplicityDistribution.C
$FMD/CreateResponseMatrices.C
$FMD/AliForwardMultiplicityDistribution.cxx/.h
$FMD/AliForwardCreateResponseMatrices.cxx/.h
$FMD/AddTaskMultiplicity.C
$FMD/AddTaskCreateResponseMatrices.C

$FMD/scripts/doUnfolding.C
$FMD/scripts/unfoldBase.C
$FMD/scripts/unfoldChi2Method.C

Additionally the unfolding ROOT framework RooUnfold should be
installed (http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html).


1) CREATE AODs
------------------------------------------------------------
Create AODs on the GRID, similarly to the standard
FORWARD/analysis2 method. 

One can use $FMD/MakeMultAOD.C macro to create the AODs, specifying
which runs/ESD pass/etc. to use.


2) CREATE BASIC MULTIPLICITY FILES
------------------------------------------------------------
Run over the (local) AODs to create basic multiplicity distribution
analysis files, and response matrices.

The basic idea is to run over data AODs to create raw multiplicity
distrbutions in a number of eta bins, and run over MC AODs to create
response matrices over the same eta bins.

The output file is a .root file with a folder for each eta-bin,
containing 'raw' multiplicity distributions or response matrices.


Multiplicity distributions:
aliroot -l 
.L MakeMultiplicityDistribution.C
MakeMultiplicityDistributions(const char* aoddir, Int_t nEvents, 
			      const char* trig, 
 			      Double_t vzMin , Double_t vzMax, 
			      Int_t lowCent, Int_t highCent, 
			      char* output, Int_t nBins)

aoddir   = directory of the local AOD files
nEvents  = number of events to run over, -1 for all
trig     = trigger word to analyse, "NSD" for Non-Single-Diffractive
           or "INEL" for inelastic
vzMin    = lower bound for z-vertex
vzMax    = upper bound for z-vertex. To have full continuous
           eta-coverage between the SPD and FMD: -4<vz<4
lowCent  = lower bound for centrality
highCent = upper bound for centrality. Only relevant for PbPb. In pp
           set both to zero.
output   = output filename
nBins    = max multiplicity. The multiplicity axis runs from 
	   [-0.5 ; nBins-0.5] in nBins bins. For pp 500 suffices, for PbPb the most
           central collisions need nBins=30000. 

response matrices:

aliroot -l
.L CreateResponseMatrices.C
CreateResponseMatrices(const char* aoddir, Int_t nEvents, 
		       const char* trig, 
		       Double_t vzMin, Double_t vzMax, 
		       char* output")

Parameters are the same as for
MakeMultiplicityDistribution(). lowCent, highCent, nBins are not
included since at the moment unfolding is not attempted for PbPb.


3) UNFOLDING
-------------------------------------------------------------
The 'raw' distributions must be unfolded to correct for detector
effects. For pp the unfolding can be done in a number of ways. I
haveused two methods; a Single Value Decomposition (SVD) method and a
Bayesen Iterative method. 

The SVD method is implemented in AliRoot by Jan Fiete, and the
Bayesian method uses the ROOT framework RooUnfold.

doUnfolding.C is a small steering routine, which calls unfoldBase.C,
which handles the unfolding itself. 

aliroot -l
.L unfoldBase.C
unfoldBase(const Char_t* outputFile, Method method, 
  	   const Char_t* responseFileName, const Char_t* dataFileName)

outputFile        = name of unfolded distributions file
method            = unfolding method, possible values are kBayes or kSvd
responseFilenName = response matrix filename
dataFileName      = uncorrected multiplicity filename


