/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/***********************************************

  Lambda Analysis Module 
  ----------------------

This Class enables the analysis of output files 
created with AliAnalysisTaskExtractV0 and
AliAnalysisTaskExtractPerformanceV0 grid tasks.

It constructs corrected Lambda, AntiLambda and
K0Short spectra. 

This version: 27th April 2012 

To be compiled as a library 
(.L AliV0Module.cxx++ or something to that effect)

--- David Dobrigkeit Chinellato
    daviddc@ifi.unicamp.br

***********************************************/

//--- For C++ ----
#include <TApplication.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <THashList.h>
#include <TDirectoryFile.h>
#include <cstdlib>
using namespace std;

//--- For ROOT --- 
#include "AliV0Module.h"

AliV0Module::AliV0Module()
{ 
  //Dummy Constructor that sets some default values
  //so that nothing too horrible will happen (hopefully)
  fWhichParticle           = "Lambda"; //Default
  fCINT1BoverINELratio     = 1;  
  fRapidityBoundary        = 0.5;
  fRealDataFile            = "";
  fMCDataFile              = "";
  fFeedDownDataFile        = "";
  fOutputDataFile          = "";

  //Selections
  fCutV0Radius                     = -1; 
  fCutDCANegToPV                   = -1;
  fCutDCAPosToPV                   = -1; 
  fCutDCAV0Daughters               = 1000;
  fCutV0CosPA                      = -2;
  fCutProperLifetime               = 1e+6;
  fCutTPCPIDNSigmas                = 1e+6; 
  fCutNSigmasForSignalExtraction   = 5;
  fCutLeastNumberOfCrossedRows     = 70;
  fCutDaughterEta                  = 0.8;
  fCutCompetingV0Rejection         = -1;

  //Set Feeddown Treatment
  fFDSwitch = "UseMCRatio";

  //Default: Use bin counting
  fFitBackgroundSwitch = kFALSE;

  //Pt Bins: undefined
  fptbinnumb = -1;
}

AliV0Module::AliV0Module(TString fParticleType)
{ 
  // Allows definition of Particle Type in analysis. 
  // Possible Options are "Lambda", "AntiLambda" and "K0Short". 
  // If some other string is given, this constructor will 
  // default to "Lambda". 
  fWhichParticle = fParticleType;
  if(fWhichParticle!="Lambda"&&fWhichParticle!="AntiLambda"&&fWhichParticle!="K0Short"){
    cout<<"Particle Type "<<fParticleType<<" unknown. Set to lambda."<<endl;
    fWhichParticle = "Lambda"; 
  }
  fCINT1BoverINELratio     = 1;  
  fRapidityBoundary        = 0.5;
  fRealDataFile            = "";
  fMCDataFile              = "";
  fFeedDownDataFile        = "";
  fOutputDataFile          = "";

  //Selections
  fCutV0Radius                     = -1; 
  fCutDCANegToPV                   = -1;
  fCutDCAPosToPV                   = -1; 
  fCutDCAV0Daughters               = 1000;
  fCutV0CosPA                      = -2;
  fCutProperLifetime               = 1e+6;
  fCutTPCPIDNSigmas                = 1e+6; 
  fCutNSigmasForSignalExtraction   = 5;
  fCutLeastNumberOfCrossedRows     = 70;
  fCutDaughterEta                  = 0.8;
  fCutCompetingV0Rejection         = -1;

  //Set Feeddown Treatment
  fFDSwitch = "UseMCRatio";

  //Default: Use bin counting
  fFitBackgroundSwitch = kFALSE;

  //Pt Bins: undefined
  fptbinnumb = -1;
}

/***********************************************
 --- Setters For Configuration ---
***********************************************/

// Filename Setters
void AliV0Module::SetRealDataFile    ( TString RealDataFilename     ){
  //Set root file containing real data candidates. 
  fRealDataFile = RealDataFilename; }
void AliV0Module::SetMCDataFile      ( TString MCDataFilename       ){
  //Set root file containing Monte Carlo data (for efficiency computation).
  fMCDataFile   = MCDataFilename;   }
void AliV0Module::SetFeedDownDataFile( TString FeedDownDataFilename ){
  //Set root file containing Monte Carlo data (for feeddown computation).
  fFeedDownDataFile = FeedDownDataFilename; }
void AliV0Module::SetOutputFile      ( TString OutputFilename       ){
  //Set root filename for the analysis output.
  fOutputDataFile = OutputFilename;
  cout<<"[AliV0Module] Set output file to \'"<<OutputFilename<<"\'."<<endl;
}

// Bin Limit Setter
void AliV0Module::SetPtBinLimits(Long_t got_fptbinnumb,const Double_t *got_fptbinlimits){
  //Function to set pt binning. First argument is the number of pt bins, second is
  //an array with bin limits. 
  fptbinnumb = got_fptbinnumb;
	for(int ix = 0;ix<fptbinnumb+1;ix++){
		fptbinlimits[ix] = got_fptbinlimits[ix];
	}
	for(int ix = 0;ix<fptbinnumb;ix++){
		fptX[ix] = (fptbinlimits[ix+1] + fptbinlimits[ix])/2.;
	}
  cout<<"[AliV0Module] Received "<<fptbinnumb<<" pt bins, set accordingly."<<endl;
}

// Rapidity Window Setter
void AliV0Module::SetRapidityWindow(Double_t got_fRapidityBoundary){
  //Set Rapidity Boundary used in analysis. 
  //Value provided will be used as upper limit for the modulus of y (|y|< given value)
  fRapidityBoundary = got_fRapidityBoundary;
  cout<<"[AliV0Module] Received "<<got_fRapidityBoundary<<" as rapidity limits, set accordingly."<<endl;
}

// CINT1B/INEL Setter for normalization to yields
void AliV0Module::SetCINT1BoverINEL(Double_t got_fCINT1BoverINEL){
  //Set CINT1B/INEL ratio (determined by collaboration).
  fCINT1BoverINELratio = got_fCINT1BoverINEL;
  cout<<"[AliV0Module] Received CINT1B/INEL = "<<got_fCINT1BoverINEL<<" to normalize with, set accordingly."<<endl;
}

// Topological Selection Setters
void AliV0Module::SetCutV0Radius(Double_t cut){
  //Set minimum decay radius for the V0 in centimeters.
  //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
  fCutV0Radius = cut;
  cout<<"[AliV0Module] Received V0 Radius (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutDCANegToPV(Double_t cut){
  //Set minimum distance of closest approach between V0 negative daughter
  //track and primary vertex (in centimeters).
  fCutDCANegToPV = cut;
  cout<<"[AliV0Module] Received DCA Negative track to PV (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutDCAPosToPV(Double_t cut){
  //Set minimum distance of closest approach between V0 positive daughter
  //track and primary vertex (in centimeters).
  fCutDCAPosToPV = cut;
  cout<<"[AliV0Module] Received DCA Positive track to PV (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutDCAV0Daughters(Double_t cut){
  //Set minimum distance of closest approach between V0 daughter
  //tracks, in sigmas. This is not in centimeters because if the 
  //tracks have been determined with ITS refit the resolution is
  //greatly improved; thus, the cut may be tighter in that case.
  //Using sigmas will therefore take the tracking method in 
  //consideration. 
  fCutDCAV0Daughters = cut;
  cout<<"[AliV0Module] Received DCA V0 Daughters (max value) = "<<cut<<endl;
}
void AliV0Module::SetCutV0CosPA(Double_t cut){
  //Set minimum value for the cosine of pointing angle of the V0.
  fCutV0CosPA = cut;
  cout<<"[AliV0Module] Received V0 Cosine of Pointing Angle (min value) = "<<cut<<endl;
}

// Other Selection Setters
void AliV0Module::SetCutProperLifetime(Double_t cut){
  //Set maximum value for m*L/p variable for the V0.
  //This is the "proper lifetime" selection and is usually called a 
  //"c*tau cut". Should be set to a value larger than the c*tau for
  //the V0 considered. 
  fCutProperLifetime = cut;
  cout<<"[AliV0Module] Received proper lifetime cut (max value) = "<<cut<<endl;
}
void AliV0Module::SetCutTPCPIDNSigmas(Double_t cut){
  //Set maximum deviation from the expected energy loss in 
  //the TPC, in multiples of sigmas as computed from the AliPIDResponse
  //object. Selection is only used in real data and should thus 
  //be very loose to ensure negligible signal loss.  
  fCutTPCPIDNSigmas = cut;
  cout<<"[AliV0Module] Received TPC N-sigmas selection (max dist from BB curve) = "<<cut<<endl;
}
void AliV0Module::SetCutSigmaForSignalExtraction(Double_t cut){
  //Set number of sigmas for the signal extraction method. The value 
  //provided is the number of sigmas from the peak position used for 
  //the peak sampling, i.e. the peak will be from [<m>-cut*sigma,<m>+cut*sigma]
  //while the background sampling regions will be from 
  //[<m>-2*cut*sigma,<m>+2*cut*sigma.
  fCutNSigmasForSignalExtraction = cut;
  cout<<"[AliV0Module] Received N-sigmas for sig. ext.: peak is (-"<<cut<<",+"<<cut<<") in sigmas"<<endl;
}
void AliV0Module::SetCutLeastNumberOfCrossedRows(Double_t cut){
  //Set smallest allowed number of TPC clusters for the V0 daughter tracks. 
  fCutLeastNumberOfCrossedRows = cut;
  cout<<"[AliV0Module] Received Least Nbr of crossed rows (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutLeastNumberOfCrossedRowsOverFindable(Double_t cut){
  //Set smallest allowed number of TPC clusters for the V0 daughter tracks. 
  fCutLeastNumberOfCrossedRowsOverFindable = cut;
  cout<<"[AliV0Module] Received Least Nbr of crossed rows over findable clusters (min value) = "<<cut<<endl;
}

void AliV0Module::SetCutDaughterEta(Double_t cut){
  //Set maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
  fCutDaughterEta = cut;
  cout<<"[AliV0Module] Received Daughter |eta| cut (max value) = "<<cut<<endl;
}

void AliV0Module::SetCutCompetingV0Rejection(Double_t cut){
  //Set rejection window around invariant mass of competing V0 species. 
  //If negative, no rejection will occur. If analysis is for Lambdas, K0s 
  //will be rejected and vice-versa. Typical values revolve around 
  //0.003 - 0.010 GeV/c^2. 
  fCutCompetingV0Rejection = cut;
  cout<<"[AliV0Module] Received Competing V0 Rejection window of +/- (in GeV/c^2) = "<<cut<<endl;
}

void AliV0Module::SetFeeddownTreatment ( TString fFDMethod ){
  //Set method used to compute charged and neutral Xi baryon feeddown to Lambda.
  //Methods allowed are:
  //"NoFD" - No Feeddown subtraction performed. 
  //"DoubleChargedXi" - Feeddown is computed for the charged Xi, and is then 
  //subtracted twice from real data. Assumes charged and neutral Xi enter 
  //the analysis in the same way and are produced at a ratio 1:1. 
  //"UseMCRatio" - The Feeddown matrix F_{ij} is filled not only with Lambdas
  //coming from charged Xis but also from neutral ones. The scaling is performed
  //according to the measured charged Xis, and since the matrix is more populated
  //a larger subtraction will occurr. This method assumes that MC correctly
  //reproduces the ratio between charged and neutral Xi baryons. Warning: this
  //should not be used with charged Xi triggered Monte Carlo data. 
  fFDSwitch = fFDMethod;
  cout<<"[AliV0Module] Received Feeddown treatment method: "<<fFDMethod<<endl;
  if( fFDMethod == "NoFD" ) 
    cout<<"[AliV0Module] ---> No Feeddown correction will be performed."<<endl;
  if( fFDMethod == "DoubleChargedXi" ) 
    cout<<"[AliV0Module] ---> Feeddown performed by doubling charged Xi correction."<<endl;
  if( fFDMethod == "UseMCRatio" ) 
    cout<<"[AliV0Module] ---> Feeddown performed by using MC neutral/charged Xi."<<endl;
}

void AliV0Module::SetFitBackground ( Bool_t fitBgSwitch ){
  //Turns on background fitting for signal extraction instead of pure 
  //bin counting. Useful, among other things, for systematics. 
  fFitBackgroundSwitch = fitBgSwitch;
}

void AliV0Module::SetDefaultCuts(){ 
   //Sets Default cuts for analysis. (adjusted for adequate pp analysis)
  cout<<"[AliV0Module] Setting default cuts for particle species: "<<fWhichParticle<<endl;
  //Set Cuts - topological
  SetCutV0Radius       (0.500);
  SetCutDCANegToPV     (0.060);
  SetCutDCAPosToPV     (0.060);
  SetCutDCAV0Daughters (1.000);
  if ( fWhichParticle == "K0Short")
    SetCutV0CosPA      (0.970);
  if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
    SetCutV0CosPA      (0.995);


  //Set Cuts - other
  SetCutProperLifetime                          ( 30);
  SetCutTPCPIDNSigmas                           (  5);
  SetCutSigmaForSignalExtraction                (  5);
  SetCutLeastNumberOfCrossedRows                ( 70);
  SetCutLeastNumberOfCrossedRowsOverFindable    (0.8);
  SetCutDaughterEta                             (0.8);
  if ( fWhichParticle == "K0Short")
    SetCutCompetingV0Rejection      (0.005);
  if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
    SetCutCompetingV0Rejection      (0.010);
}

TString AliV0Module::IntToString(int input){
  //Integer to TString Converter
	char dummyChar[50];
	sprintf(dummyChar, "%d", input);
	TString outputstring = dummyChar;
	return outputstring;
}

TString AliV0Module::DoubleToString(double input){
  //Double to TString Converter
	char dummyChar[50];
	sprintf(dummyChar, "%.3f", input);
	TString outputstring = dummyChar;
	return outputstring;
}

Double_t AliV0Module::ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ){
  //Error in a Ratio
        if(B!=0){
                Double_t errorfromtop = Aerr*Aerr / (B*B) ;
                Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
                return TMath::Sqrt( errorfromtop + errorfrombottom );
        }
        return 1;
}

Double_t AliV0Module::MyGeant3FlukaCorrectionForProtons(const Double_t *x, const Double_t *par){
        //Parametrization Used for Geant3/Fluka Correction for protons
        //Credit: Antonin Maire
        return 1 - par[0]*TMath::Exp(par[1]*x[0]) + par[2];        
}


Double_t AliV0Module::MyGeant3FlukaCorrectionForAntiProtons(const Double_t *x, const Double_t *par){
        // Parametrization Used for Geant3/Fluka Correction for antiprotons
        // Credit: Antonin Maire
        //
        // La fonction A*TMath::Erf( B.x ) ne marche pas à bas pt.
        // Différentes fonctions ont été testées.
        // La fonction suivante semble donner de bons résultats
        // On peut jouer sur la puissance n du terme "1/x^n" pour repousser le pt auquel la fonction prend la valeur 1.0
        // Ici, pour n = 0.2, la fonction prend la valeur 0.9990 en pt = 10 GeV/c
        
        return 1 - par[0]*TMath::Exp(par[1]*x[0]) + par[2] + par[3]*1/TMath::Power(x[0], 0.2)*TMath::Log(x[0]);       
}

Double_t AliV0Module::MyLevyPtXi(const Double_t *pt, const Double_t *par)
{
  //Levy Fit Function
  Double_t lMass  = 1.32171; //Xi Mass
  Double_t ldNdy  = par[0];
  Double_t lTemp = par[1];
  Double_t lPower = par[2];

  Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
  Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);

  return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
  //return ldNdy * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}

Double_t AliV0Module::MyBgPol1(const Double_t *x, const Double_t *par)
{
  //Function for background fitting, rejects peak region
        if ( x[0] > par[2] && x[0] < par[3]) {
                TF1::RejectPoint();
                return 0;
        }
        return par[0] + par[1]*x[0];
}

Double_t AliV0Module::MyBgPolToEval1(const Double_t *x, const Double_t *par)
{
  //Just a plain linear function.
        return par[0] + par[1]*x[0];
}

Double_t AliV0Module::RoundToThousandth( const Double_t lToRound ){
  //Round any number to a hundredth... 
  //Well, within machine precision...
  return TMath::Nint( 1000 * lToRound ) * 0.001;
}

void AliV0Module::DoAnalysis(){
  //----------------------
  // V0 Analysis Function
  //----------------------
  //
  //Consists of the following steps: 
  // 
  // (1) Loop over Real data candidates, acquire peak position and widths.
  // (2) Loop over Real data candidates, extract signal with variable extraction 
  //     areas. Two loops: totally avoids binning data, allowing for sub-MeV/c^2 
  //     granularity in signal extraction areas. 
  // (3) Loop over MC data reconstructed candidates, find associated-to-MC primary 
  //     candidates for efficiency numerator.
  // (4) (if Lambda or AntiLambda) Perform Geant3/Fluka correction.
  // (5) Get generated primary V0 histograms from MC file for efficiency denominator. 
  // (6) (if Lambda or AntiLambda + FD correction enabled) Open MC file for feeddown 
  //     subtraction. Loop over candidates, find Lambda associated with primary Xi 
  //     to fill feeddown matrix. Scale the feeddown contribution according to real-life
  //     measured Xi- production under specified feeddown subtraction scheme (see 
  //     function SetFeeddownTreatment). Perform subtraction of raw Lambda or AntiLambda
  //     count estimated to be coming from charged or neutral Xi. 
  // (7) Perform detection efficiency correction and compute final corrected spectra.
  //
  //
  // Normalization: 
  //  --- Number of Inelastic Events estimated to be:
  //      <Number of CINT1B triggers> / <CINT1B/INEL ratio>
  //  
  //  --- Efficiency denominator: filled at pre-physics selection stage. All signal 
  //      losses are therefore estimated from Monte Carlo. 
  //
  //  --- Pileup: not included in MC. Fraction of events removed by IsPileupFromSPD 
  //      is artificially removed from Number of Inelastic Events as well.    
  //
  // Output: written to specified file (with SetOutputFile(...)). 
  //   --- includes spectra and efficiencies in base directory. 
  //   --- includes a number of different subdirectories storing plots such as
  //       such as invariant mass histos, signal extraction, pt resolution, G3/F 
  //       related distributions, feeddown related distributions, and so on.

  //Set Batch Mode: Ignore all canvases in this section 
  gROOT->SetBatch (kTRUE);

  cout<<"======================================"<<endl;
  cout<<" -------- Spectra Extraction -------- "<<endl;
  cout<<"======================================"<<endl;
  cout<<endl;
  if(fptbinnumb == -1){
    cout<<"[AliV0Module] It's unclear what kind of pt binning you want."<<endl;
    cout<<"[AliV0Module] Most likely you forgot to set it using SetPtBinLimits..."<<endl;
    cout<<"[AliV0Module] Analysis will NOT be done. Returning."<<endl;
    return;
  }
  cout<<"--------------- Configuration --------------------------"<<endl;
  cout<<" Analysed Particle.............: "<<fWhichParticle<<endl;
  cout<<" Rapidity Window...............: "<<fRapidityBoundary<<endl;
  cout<<" CINT1B/INEL Ratio used........: "<<fCINT1BoverINELratio<<endl;
  cout<<" V0 Decay Radius...............: "<<fCutV0Radius<<endl;
  cout<<" DCA Negative track to PV......: "<<fCutDCANegToPV<<endl;
  cout<<" DCA Positive track to PV......: "<<fCutDCAPosToPV<<endl;
  cout<<" DCA V0 Daughters..............: "<<fCutDCAV0Daughters<<endl;
  cout<<" Cosine of Pointing Angle V0...: "<<fCutV0CosPA<<endl;
  cout<<" Proper Lifetime cut (cm)......: "<<fCutProperLifetime<<endl;
  cout<<" TPC dE/dx sigmas cut (Real)...: "<<fCutTPCPIDNSigmas<<endl;
  cout<<" CutNSigmasForSignalExtraction.: "<<fCutNSigmasForSignalExtraction<<endl;
  cout<<" Least # of Crossed Rows.......: "<<fCutLeastNumberOfCrossedRows<<endl;
  cout<<" Least # of C. R. over findable: "<<fCutLeastNumberOfCrossedRowsOverFindable<<endl;
  cout<<" Daughter Track |eta| < .......: "<<fCutDaughterEta<<endl;
  cout<<" Competing V0 Reject. (GeV/c^2): "<<fCutCompetingV0Rejection<<endl;
  cout<<"--------------- File Names -----------------------------"<<endl;
  cout<<" Real Data File................: "<<fRealDataFile<<endl;
  cout<<" MC File.......................: "<<fMCDataFile<<endl;
  cout<<" MC File for feeddown..........: "<<fFeedDownDataFile<<endl;
  cout<<" Analysis output file..........: "<<fOutputDataFile<<endl;
  cout<<"--------------------------------------------------------"<<endl;
  cout<<endl;

  //=== Real data loop 1: Acquire peak positions, widths ===
  //--- Preparing... ---
  //defining helping histogram - only used to FindBin Index!========
	TH1F* fHistPt 		= new TH1F("fHistPt","Dummy;p_{T} (GeV/c);Counts",fptbinnumb,fptbinlimits);
  //Peak Position Histograms 
	TH1F* fHistPeakPosition	= new TH1F("fHistPeakPosition","Peak Position (Real);p_{T} (GeV/c);Peak Position",fptbinnumb,fptbinlimits);
	TH1F* fHistPeakPositionMC	= new TH1F("fHistPeakPositionMC","Peak Position (MC);p_{T} (GeV/c);Peak Position",fptbinnumb,fptbinlimits);
  //Signal to Noise Histograms
	TH1F* fHistSigToNoise	  = new TH1F("fHistSigToNoise","Signal to Noise Ratio (real);p_{T} (GeV/c);Sig / Bg",fptbinnumb,fptbinlimits);
	TH1F* fHistSigToNoiseMC	= new TH1F("fHistSigToNoiseMC","Signal to Noise Ratio (MC);p_{T} (GeV/c);Sig / Bg",fptbinnumb,fptbinlimits);

  //Signal Extraction Range Histogram
  TH1F* fHistSignalExtractionRange	= new TH1F("fHistSignalExtractionRange","Sig. Ext. Range;p_{T} (GeV/c);Range",fptbinnumb,fptbinlimits);

   //Resolution Histogram (filled with MC)
	TH2F* f2dHistPtResolution 		= new TH2F("f2dHistPtResolution","p_{t} Resolution;p_{t} (reco);p_{t} (mc)",fptbinnumb,fptbinlimits, fptbinnumb, fptbinlimits);

  /********************************************************

    ---> Let's Remember the limits of the data we're analyzing! 
    ---> Important so that we don't try signal extraction outside 
    ---> the bundaries of available data in the Tree object.

    From AliAnalysisTaskExtractV0.cxx 

    //Second Selection: rough 20-sigma band, parametric. 
    //K0Short: Enough to parametrize peak broadening with linear function.    
    Double_t UpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*tree_lPt; 
    Double_t LowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*tree_lPt;

    //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
    //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
    Double_t UpperLimitLambda = (1.13688e+00) + (5.27838e-03)*tree_lPt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*tree_lPt); 
    Double_t LowerLimitLambda = (1.09501e+00) - (5.23272e-03)*tree_lPt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*tree_lPt);

  ********************************************************/

  TF1 *fKDataUpper = new TF1("fKDataUpper","[0]+[1]*x",0,20);
  TF1 *fKDataLower = new TF1("fKDataLower","[0]+[1]*x",0,20);
  TF1 *fLDataUpper = new TF1("fLDataUpper","[0]+[1]*x+[2]*TMath::Exp([3]*x)",0,20);
  TF1 *fLDataLower = new TF1("fLDataLower","[0]+[1]*x+[2]*TMath::Exp([3]*x)",0,20);

  fKDataUpper->SetParameter(0, 5.63707e-01);
  fKDataUpper->SetParameter(1, 1.14979e-02);
  fKDataLower->SetParameter(0, 4.30006e-01);
  fKDataLower->SetParameter(1,-1.10029e-02);
  
  fLDataUpper->SetParameter(0, 1.13688e+00);
  fLDataUpper->SetParameter(1, 5.27838e-03);
  fLDataUpper->SetParameter(2, 8.42220e-02);
  fLDataUpper->SetParameter(3,-3.80595e+00);

  fLDataLower->SetParameter(0, 1.09501e+00);
  fLDataLower->SetParameter(1,-5.23272e-03);
  fLDataLower->SetParameter(2,-7.52690e-02);
  fLDataLower->SetParameter(3,-3.46339e+00);

  Int_t lWeAreAtBin = 0;
  Double_t lParticleMass = -1;          //needed for proper lifetime selection
  Double_t lCompetingParticleMass = -1; //needed for Competing V0 Rejection
  Double_t lHistoLowerBoundary = -1;
  Double_t lHistoUpperBoundary = -1;
  Long_t lHistoNBins = -1;
  if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){ 
    lParticleMass          = 1.115683;
    lCompetingParticleMass = 0.4976;
    lHistoLowerBoundary = 1.116 - 0.1; 
    lHistoUpperBoundary = 1.116 + 0.1; 
    lHistoNBins = 200; //1MeV/c^2 binning
  }
  if ( fWhichParticle == "K0Short" ){ 
    lParticleMass          = 0.4976;
    lCompetingParticleMass = 1.115683;
    lHistoLowerBoundary = 0.498 - 0.15; 
    lHistoUpperBoundary = 0.498 + 0.15;
    lHistoNBins = 300; //1MeV/c^2 binning 
  }
	//Setting Up Base Histograms to use===============================
  TH1F* lHistoFullV0[100];
  TH1F* lHistoFullV0MC[100];
  TCanvas* lCanvasHistoFullV0[100];
  TCanvas* lCanvasHistoFullV0MC[100];
  char histname[80]; 	TString bindescription = "";
  for(Int_t ihist=0;ihist<100;ihist++){
    //Histo For Real Data
    sprintf(histname,"lHistoFullV0%i",ihist);
    if(fWhichParticle == "Lambda")     bindescription="#Lambda, bin #";
    if(fWhichParticle == "AntiLambda") bindescription="#bar{#Lambda}, bin #";
    if(fWhichParticle == "K0Short")    bindescription="K^{0}_{S}, bin #";
		bindescription.Append(IntToString( ihist ));
    if ( ihist < fptbinnumb ){
		  bindescription.Append(", ");
		  bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
		  bindescription.Append("-");
		  bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
		  bindescription.Append("GeV/c");
    }
    lHistoFullV0[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
    lHistoFullV0[ihist]->SetTitle(bindescription);
    sprintf(histname,"lCanvasHistoFullV0%i",ihist);
    lCanvasHistoFullV0[ihist] = new TCanvas(histname, bindescription, 800, 600);
    lCanvasHistoFullV0[ihist] -> SetFillColor(kWhite);
    lCanvasHistoFullV0[ihist] -> SetLeftMargin( 0.175 );
    lCanvasHistoFullV0[ihist] -> SetBottomMargin( 0.175 );
    

    //Histo for MC 
    sprintf(histname,"lHistoFullV0MC%i",ihist);
    if(fWhichParticle == "Lambda")     bindescription="#Lambda, MC, bin #";
    if(fWhichParticle == "AntiLambda") bindescription="#bar{#Lambda}, MC, bin #";
    if(fWhichParticle == "K0Short")    bindescription="K^{0}_{S}, MC, bin #";
		bindescription.Append(IntToString( ihist ));
    if ( ihist < fptbinnumb ){
		  bindescription.Append(", ");
		  bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
		  bindescription.Append("-");
		  bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
		  bindescription.Append("GeV/c");
    }
    lHistoFullV0MC[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
    lHistoFullV0MC[ihist]->SetTitle(bindescription);
    sprintf(histname,"lCanvasHistoFullV0MC%i",ihist);
    lCanvasHistoFullV0MC[ihist] = new TCanvas(histname, bindescription, 800, 600);
    lCanvasHistoFullV0MC[ihist] -> SetFillColor(kWhite);
    lCanvasHistoFullV0MC[ihist] -> SetLeftMargin( 0.175 );
    lCanvasHistoFullV0MC[ihist] -> SetBottomMargin( 0.175 );
  }
	//================================================================	
  cout<<endl;


  cout<<"--------------- Open Real Data File --------------------"<<endl;
	TFile* file = TFile::Open(fRealDataFile, "READ");
	file->cd("PWG2CheckLambda_PP");
	TList* v0list  = (TList*)file->FindObjectAny("clistV0");
  TTree* lTree; 
  TH1F* fHistV0MultiplicityForTrigEvt;
  TH1F* fHistV0MultiplicityForSelEvtNoTPCOnly;
  TH1F* fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup;
  lTree = (TTree*)file->FindObjectAny("fTree");
  fHistV0MultiplicityForTrigEvt = (TH1F*)v0list->FindObject("fHistV0MultiplicityForTrigEvt");
  fHistV0MultiplicityForSelEvtNoTPCOnly = (TH1F*)v0list->FindObject("fHistV0MultiplicityForSelEvtNoTPCOnly");
  fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup = (TH1F*)v0list->FindObject("fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup");

	Long_t lNCandidates = lTree->GetEntries();
	Long_t lNEvents = fHistV0MultiplicityForTrigEvt->GetEntries();
  Double_t lEvtFracAfterPileup = 
    ((double)(fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup->GetEntries()))/
    ((double)(fHistV0MultiplicityForSelEvtNoTPCOnly->GetEntries()));

  Double_t lNInelasticEvents = lEvtFracAfterPileup * (lNEvents / fCINT1BoverINELratio);

  cout<<" Number of CINT1B triggers.....: "<<lNEvents<<endl;
	cout<<" Pileup Fraction...............: "<<1 - lEvtFracAfterPileup<<endl;
  cout<<" Estimated N_{INEL}............: "<<lNInelasticEvents<<endl;
	cout<<" V0 Candidates.................: "<<lNCandidates<<endl;
  cout<<"--------------------------------------------------------"<<endl;

	//Variable Definition=============================================
	//Kinematic
	Float_t lPt, lRap, lPtMC, lNegEta, lPosEta;
  //Invariant Masses
  Float_t lInvariantMass, lInvariantMassCompetingOne, lInvariantMassCompetingTwo; //wildcard
	//DCA Variables
	Float_t lDcaV0Daughters;
	Float_t lDcaPosToPrimVertex,  lDcaNegToPrimVertex; 
	//Cosine of Pointing Angle variable
	Float_t lV0CosinePointingAngle; 
  //Decay Radius and distance over total momentum
  Float_t lV0Radius, lDistOverTotMom;
  //Least Number of TPC Clusters
  Int_t lLeastNbrCrossedRows;
  Float_t lLeastNbrCrossedRowsOverFindable;
  //TPC dE/dx acquired with AliPIDResponse class
  Float_t lNSigmasPosProton,lNSigmasNegProton,lNSigmasPosPion,lNSigmasNegPion;
	//================================================================

	//Linking to Tree=================================================

  //--- Base Variables ----------------------------------------------
	lTree->SetBranchAddress("fTreeVariablePosEta",&lPosEta);
	lTree->SetBranchAddress("fTreeVariableNegEta",&lNegEta);
	lTree->SetBranchAddress("fTreeVariablePt",&lPt);
  if ( fWhichParticle == "Lambda"      )  lTree->SetBranchAddress("fTreeVariableInvMassLambda",&lInvariantMass);
  if ( fWhichParticle == "AntiLambda"  ) 	lTree->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMass);
  if ( fWhichParticle == "K0Short"     ) 	lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMass);
  if ( fWhichParticle == "Lambda"      ){ //For symmetry of computation...
    lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
    lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingTwo);
  }
  if ( fWhichParticle == "AntiLambda"  ){ //For symmetry of computation...
    lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
    lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingTwo);
  }
  if ( fWhichParticle == "K0Short"     ){
    lTree->SetBranchAddress("fTreeVariableInvMassLambda"    ,&lInvariantMassCompetingOne);
    lTree->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMassCompetingTwo);
  }
  if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
  	lTree->SetBranchAddress("fTreeVariableRapLambda",&lRap);
  if ( fWhichParticle == "K0Short" )
  	lTree->SetBranchAddress("fTreeVariableRapK0Short",&lRap);
	lTree->SetBranchAddress("fTreeVariableDistOverTotMom",&lDistOverTotMom);
	lTree->SetBranchAddress("fTreeVariableLeastNbrCrossedRows",&lLeastNbrCrossedRows);
	lTree->SetBranchAddress("fTreeVariableLeastRatioCrossedRowsOverFindable",&lLeastNbrCrossedRowsOverFindable);
  //--- TPC dEdx Variables ------------------------------------------
	lTree->SetBranchAddress("fTreeVariableNSigmasPosProton",&lNSigmasPosProton);
	lTree->SetBranchAddress("fTreeVariableNSigmasNegProton",&lNSigmasNegProton);
	lTree->SetBranchAddress("fTreeVariableNSigmasPosPion",&lNSigmasPosPion);
	lTree->SetBranchAddress("fTreeVariableNSigmasNegPion",&lNSigmasNegPion);
  //--- Topological selection variables -----------------------------
	lTree->SetBranchAddress("fTreeVariableV0Radius",&lV0Radius);
	lTree->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",&lDcaNegToPrimVertex);
	lTree->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",&lDcaPosToPrimVertex);
	lTree->SetBranchAddress("fTreeVariableDcaV0Daughters",&lDcaV0Daughters);
	lTree->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",&lV0CosinePointingAngle);
	//================================================================


  cout<<"--------------- Real Data File Loop 1 ------------------"<<endl;
  Long_t lOneTenthOfNCandidates = ((double)(lNCandidates) / 10. );
  for(Long_t icand = 0;icand<lNCandidates;icand++){
		lTree->GetEntry(icand);
    if( icand % lOneTenthOfNCandidates == 0 ) 
    	cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidates<<" ( "<<(long)(((double)(icand)/(double)(lNCandidates))*(100.+1e-3))<<"% )"<<endl;
		if(TMath::Abs(lRap)<fRapidityBoundary &&
            TMath::Abs(lNegEta)       <= fCutDaughterEta               &&                   
            TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
            lV0Radius                 >= fCutV0Radius                  &&
            lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
            lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
            lDcaV0Daughters           <= fCutDCAV0Daughters            &&
            lV0CosinePointingAngle    >= fCutV0CosPA                   && 
            lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
            lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
            lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
            TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
            TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
          ( //official response code
          ( fWhichParticle == "Lambda"
          && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas
          && TMath::Abs(lNSigmasPosProton) <= fCutTPCPIDNSigmas) || 
          ( fWhichParticle == "AntiLambda"
          && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
          && TMath::Abs(lNSigmasNegProton) <= fCutTPCPIDNSigmas) || 
          ( fWhichParticle == "K0Short"
          && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
          && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas)
          )
      ){
		  lWeAreAtBin = fHistPt->FindBin(lPt)-1;
		  if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment
		  lHistoFullV0[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass
		}
	}
  cout<<"--------------- Loop Completed -------------------------"<<endl;
  cout<<endl;

  cout<<"--------------- Peak Finding (gauss+linear) ------------"<<endl;
 //== Variables for holding peak position, width ==============
  Double_t lPeakPosition[100];
  Double_t lPeakWidth[100]; 
  Double_t lLeftBgLeftLimit[100];
  Double_t lLeftBgRightLimit[100];
  Double_t lPeakLeftLimit[100];
  Double_t lPeakRightLimit[100];
  Double_t lRightBgLeftLimit[100];
  Double_t lRightBgRightLimit[100];
  //May be needed if bg areas are different in the future
  //Double_t lScaleFactor[100];

  TLine *lLineLeftMost[100];
  TLine *lLineLeft[100];
  TLine *lLineRight[100];
  TLine *lLineRightMost[100];

  TLine *lLineLeftMostMC[100];
  TLine *lLineLeftMC[100];
  TLine *lLineRightMC[100];
  TLine *lLineRightMostMC[100];

  char fgausname[100]; 
  TF1 *fgausPt[100];
  for(Int_t ibin = 0; ibin<fptbinnumb; ibin++){ 
    cout<<"---> Peak Finding, bin #"<<ibin<<"..."<<endl;
    sprintf(fgausname,"fGausPt%i",ibin);
   	if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){
      fgausPt[ibin]= new TF1(fgausname,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]",1.116-0.040,1.116+0.040);
  	  fgausPt[ibin]->SetParameter(1,1.116);
  	  fgausPt[ibin]->SetParameter(2,0.0025); 
      fgausPt[ibin]->SetParLimits(1,1,1.2);
      fgausPt[ibin]->SetParLimits(2,0.001,0.02);
    }
   	if ( fWhichParticle == "K0Short"){ 
      fgausPt[ibin]= new TF1(fgausname,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]",0.498-0.06,0.498+0.060);
  	  fgausPt[ibin]->SetParameter(1,0.498);
  	  fgausPt[ibin]->SetParameter(2,0.004); 
    }
    fgausPt[ibin]->SetParameter(0,lHistoFullV0[ibin]->GetMaximum() * 0.9);


	  fgausPt[ibin]->SetParameter(3,0); 
	  fgausPt[ibin]->SetParameter(4,lHistoFullV0[ibin]->GetMaximum() * 0.1); 
    lHistoFullV0[ibin]->Fit(fgausname,"QREM0");
    lPeakPosition[ibin] = fgausPt[ibin]->GetParameter(1);
    lPeakWidth[ibin] = fgausPt[ibin]->GetParameter(2);
    cout<<"---> ["<<fptbinlimits[ibin]<<" - "<<fptbinlimits[ibin+1]<<" GeV/c]\tPeak at: "<<lPeakPosition[ibin]<<", sigma = "<<lPeakWidth[ibin]<<endl;
    //Find Corresponding Limits In this bin
    lLeftBgLeftLimit[ibin]  = lPeakPosition[ibin] - 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin]; 
    lLeftBgRightLimit[ibin] = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin]; 
    lPeakLeftLimit[ibin]    = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin]; 
    lPeakRightLimit[ibin]   = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin]; 
    lRightBgLeftLimit[ibin] = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin]; 
    lRightBgRightLimit[ibin]= lPeakPosition[ibin] + 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin]; 
    //lScaleFactor[ibin] = (lPeakRightLimit[ibin] - lPeakLeftLimit[ibin]) / ( (lLeftBgRightLimit[ibin]-lLeftBgLeftLimit[ibin]) + (lRightBgRightLimit[ibin] - lRightBgLeftLimit[ibin]) );
    fHistPeakPosition->SetBinContent(ibin+1, lPeakPosition[ibin]);
    fHistPeakPosition->SetBinError(ibin+1, lPeakWidth[ibin]);
    //Create Signal Extraction Range Histogram
    fHistSignalExtractionRange->SetBinContent(ibin+1, lPeakPosition[ibin]);
    fHistSignalExtractionRange->SetBinError(ibin+1, 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin] );
    //Create appropriate TLine Objects for Canvases 
    lLineLeftMost[ibin]  = new TLine( lLeftBgLeftLimit[ibin],   0, lLeftBgLeftLimit[ibin],   lHistoFullV0[ibin]->GetMaximum() * 0.95 );
    lLineLeft[ibin]      = new TLine( lLeftBgRightLimit[ibin],  0, lLeftBgRightLimit[ibin],  lHistoFullV0[ibin]->GetMaximum() * 0.95 );
    lLineRight[ibin]     = new TLine( lRightBgLeftLimit[ibin],  0, lRightBgLeftLimit[ibin],  lHistoFullV0[ibin]->GetMaximum() * 0.95 );
    lLineRightMost[ibin] = new TLine( lRightBgRightLimit[ibin], 0, lRightBgRightLimit[ibin], lHistoFullV0[ibin]->GetMaximum() * 0.95 );

    //Preparing Canvas for storing...
    lCanvasHistoFullV0[ibin]->cd();
    lHistoFullV0[ibin]->Draw();
    lLineLeftMost[ibin]->Draw();
    lLineLeft[ibin]->Draw();
    lLineRight[ibin]->Draw();
    lLineRightMost[ibin]->Draw();
    fgausPt[ibin]->Draw("same");
  
  }
  cout<<"--------------- Peak Finding Finished ------------------"<<endl;

  //Defining Signal holding variables===============================
  Double_t lSigRealV0[100];
  Double_t lSigErrRealV0[100];
  Double_t lSigMCV0[100];
  Double_t lSigErrMCV0[100];
  //================================================================
	Double_t lLeftPlusRightBgV0[100];
	Double_t lSigPlusCenterBgV0[100];
	Double_t lLeftPlusRightBgV0MC[100];
	Double_t lSigPlusCenterBgV0MC[100];
  for(Long_t n1 = 0; n1<100; n1++){
     lLeftPlusRightBgV0[n1]=0;
     lSigPlusCenterBgV0[n1]=0;
     lLeftPlusRightBgV0MC[n1]=0;
     lSigPlusCenterBgV0MC[n1]=0;
  }
  //================================================================
  cout<<endl;
  cout<<"--------------- Real Data File Loop 2 ------------------"<<endl;
  for(Long_t icand = 0;icand<lNCandidates;icand++){
		lTree->GetEntry(icand);
    if( icand % lOneTenthOfNCandidates == 0 ) 
    	cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidates<<" ( "<<(long)(((double)(icand)/(double)(lNCandidates))*(100.+1e-3))<<"% )"<<endl;
		if(TMath::Abs(lRap)<fRapidityBoundary &&
            TMath::Abs(lNegEta)       <= fCutDaughterEta               &&                   
            TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
            lV0Radius                 >= fCutV0Radius                  &&
            lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
            lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
            lDcaV0Daughters           <= fCutDCAV0Daughters            &&
            lV0CosinePointingAngle    >= fCutV0CosPA                   && 
            lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
            lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
            lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
            TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
            TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
          ( //official response code
          ( fWhichParticle == "Lambda"
          && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas
          && TMath::Abs(lNSigmasPosProton) <= fCutTPCPIDNSigmas) || 
          ( fWhichParticle == "AntiLambda"
          && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
          && TMath::Abs(lNSigmasNegProton) <= fCutTPCPIDNSigmas) || 
          ( fWhichParticle == "K0Short"
          && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
          && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas)
          )
      ){ // Start Entry Loop
		  lWeAreAtBin = fHistPt->FindBin(lPt)-1;
		  if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment

      //Extract left and right background areas and peak 
      //--- Left Background Sample
			if(lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]    && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]   ){ lLeftPlusRightBgV0[lWeAreAtBin]++; }
			if(lInvariantMass>lRightBgLeftLimit[lWeAreAtBin]   && lInvariantMass<lRightBgRightLimit[lWeAreAtBin]  ){ lLeftPlusRightBgV0[lWeAreAtBin]++; }
      //--- Peak Region
			if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin]      && lInvariantMass<lPeakRightLimit[lWeAreAtBin]     ){ lSigPlusCenterBgV0[lWeAreAtBin]++; }
		} // End Entry Loop
	}
  cout<<"--------------- Loop Completed -------------------------"<<endl;
  cout<<endl;
  cout<<"--------------- Memory Cleanup -------------------------"<<endl;
  if ( v0list ) { 
    v0list->Delete();
    delete v0list;
  } 

  file->Close("R");
  file->Delete();
  delete file; 
  cout<<endl;

  TF1 *lfitNoise[50];
  TF1 *lSampleNoise[50];
  char lFitNameOne[50];

  if( fFitBackgroundSwitch ){ 
    cout<<"--------------- Backgrounds: fitting with linear -------"<<endl;
    for(long i=0; i<fptbinnumb; i++){
      //Define Function to Fit Background 
      sprintf(lFitNameOne,"lfitNoise%i",((int)(i)));
      //cout<<"creating fitnoise, named "<<lFitNameOne<<endl;
      lfitNoise[i] = new TF1(lFitNameOne, this, &AliV0Module::MyBgPol1, 
        RoundToThousandth ( lLeftBgLeftLimit[i]   ),
        RoundToThousandth ( lRightBgRightLimit[i] ), 4 , "AliV0Module", "MyBgPol1");
      lfitNoise[i] -> FixParameter(2,RoundToThousandth ( lLeftBgRightLimit[i] ) );
      lfitNoise[i] -> FixParameter(3,RoundToThousandth ( lRightBgLeftLimit[i] ) );
      lfitNoise[i] -> SetParameter(0,lLeftPlusRightBgV0[i] * lHistoFullV0[i]->GetBinWidth(5) / (lRightBgLeftLimit[i]-lLeftBgRightLimit[i] + 1e-6 ) );
      lfitNoise[i] -> SetParameter(1,0);
      cout<<"Guessed Parameter 0 fot "<<lFitNameOne<<" to be "<<lLeftPlusRightBgV0[i] * lHistoFullV0[i]->GetBinWidth(5) / (lRightBgLeftLimit[i]-lLeftBgRightLimit[i] + 1e-6 )<<endl;
      sprintf(lFitNameOne,"lSampleNoise%i",((int)(i)));
      
      //Define Function to Sample Background
      //cout<<"creating sample "<<i<<endl;
      lSampleNoise[i] = new TF1(lFitNameOne, this, &AliV0Module::MyBgPolToEval1, 
        RoundToThousandth ( lLeftBgLeftLimit[i]   ),
        RoundToThousandth ( lRightBgRightLimit[i] ), 2 , "AliV0Module", "MyBgPolToEval1");
    }    
    for(long i=0; i<fptbinnumb; i++){
      //cout<<"Fitting function for bin "<<i<<", get name = "<<lfitNoise[i]->GetName()<<endl;
      sprintf(lFitNameOne,"lfitNoise%i",((int)(i)));
      lHistoFullV0[i] -> Fit(lFitNameOne,"LLrie+0");
      lSampleNoise[i]->SetParameter(0, lfitNoise[i]->GetParameter(0) );
      lSampleNoise[i]->SetParameter(1, lfitNoise[i]->GetParameter(1) );
    }
    for(long i=0; i<fptbinnumb; i++){
      cout<<"Overriding Background info: Was "<<lLeftPlusRightBgV0[i]<<", is now "<<lSampleNoise[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0[i]->GetBinWidth(5)<<endl;
      lLeftPlusRightBgV0[i] = lSampleNoise[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0[i]->GetBinWidth(5);
    }
    cout<<"--------------- Fitting Finished! ----------------------"<<endl;
  }

  //=============================================================
  // Compute Signal + Sig to noise
  for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
      //Signal computation
      lSigRealV0[ipoint] = lSigPlusCenterBgV0[ipoint] - lLeftPlusRightBgV0[ipoint];
      lSigErrRealV0[ipoint] = TMath::Sqrt(lSigPlusCenterBgV0[ipoint]+lLeftPlusRightBgV0[ipoint]);
      //Error Formula: Equivalent to Sqrt(S+B+B) = Sqrt(S+2B)

      //Signal-to-noise computation
      if( lLeftPlusRightBgV0[ipoint] != 0 ){
        fHistSigToNoise->SetBinContent(ipoint+1, (lSigPlusCenterBgV0[ipoint] - lLeftPlusRightBgV0[ipoint])/lLeftPlusRightBgV0[ipoint] );
      }else{
        fHistSigToNoise->SetBinContent(ipoint+1, -1) ; //-1 means: no background
      }
  }
  //=============================================================
 
  cout<<"--------------- Extracted Signal -----------------------"<<endl;
  for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
      cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<lSigRealV0[ipoint]<<" +/- "<<lSigErrRealV0[ipoint]<<endl;
  }
  cout<<"--------------------------------------------------------"<<endl;

  //=============================================================
  // Preparations for MC loop 
  //Extra info available only in MC
  Int_t lPID = 0; 
  Float_t lNegTransvMomentumMC, lPosTransvMomentumMC;
  Int_t lPIDPositive = 0;
  Int_t lPIDNegative = 0;
  Int_t lPIDMother = -1; Float_t lPtMother = -1;
  Int_t lPrimaryStatus = 0; 
  Int_t lPrimaryStatusMother = 0; 
  //=============================================================
  cout<<endl;
  //Prepare containers for Geant3/fluka
  //Will only be used and saved if we're dealing with Lambdas
  TH1F *lProtonMomentum[100]; 
  char lNameOne[100];
  for(Int_t ibin=0;ibin<100;ibin++){ 
    sprintf(lNameOne,"lProtonMomentumBin%i",ibin);
    lProtonMomentum[ibin] = new TH1F(lNameOne,"",800,0,20);
    if(fWhichParticle == "Lambda")     bindescription="#Lambda, bin #";
    if(fWhichParticle == "AntiLambda") bindescription="#bar{#Lambda}, bin #";
    if(fWhichParticle == "K0Short")    bindescription="K^{0}_{S}, bin #";
		bindescription.Append(IntToString( ibin ));
    if ( ibin < fptbinnumb ){
		  bindescription.Append(", ");
		  bindescription.Append(DoubleToString(fptbinlimits[ ibin ]));
		  bindescription.Append("-");
		  bindescription.Append(DoubleToString(fptbinlimits[ ibin+1 ]));
		  bindescription.Append("GeV/c");
    }
    lProtonMomentum[ibin]->SetTitle(bindescription);
  }

	//MC File Acquisition=============================================
  cout<<"--------------- Opening MC file ------------------------"<<endl;
	TFile* fileMC = TFile::Open(fMCDataFile, "READ");
	fileMC->cd("PWG2CheckPerformanceLambda_PP_MC");
	TList* v0listMC  = (TList*)fileMC->FindObjectAny("clistV0MC");
  TTree *lTreeMC;
  TH1F* fHistMultiplicityBeforeTrigSelMC;
  lTreeMC = (TTree*)fileMC->FindObjectAny("fTree");
  fHistMultiplicityBeforeTrigSelMC = (TH1F*)v0listMC->FindObject("fHistMultiplicityBeforeTrigSel");
 
	Long_t lNCandidatesMC = lTreeMC->GetEntries();
	Long_t lNEventsMC = fHistMultiplicityBeforeTrigSelMC->GetEntries();
	cout<<" MC Events (before trig sel)...: "<<lNEventsMC<<endl;
	cout<<" Cascade MC Candidates.........: "<<lNCandidatesMC<<endl;
  cout<<"--------------------------------------------------------"<<endl;
	//================================================================

  //Linking to Tree=================================================
  //--- Base Variables ----------------------------------------------
	lTreeMC->SetBranchAddress("fTreeVariablePosEta",&lPosEta);
	lTreeMC->SetBranchAddress("fTreeVariableNegEta",&lNegEta);
	lTreeMC->SetBranchAddress("fTreeVariablePrimaryStatus",&lPrimaryStatus);
	lTreeMC->SetBranchAddress("fTreeVariablePrimaryStatusMother",&lPrimaryStatusMother);
	lTreeMC->SetBranchAddress("fTreeVariablePt",&lPt);
	lTreeMC->SetBranchAddress("fTreeVariablePtXiMother",&lPtMother);
	lTreeMC->SetBranchAddress("fTreeVariablePtMC",&lPtMC);
  if ( fWhichParticle == "Lambda"      )  lTreeMC->SetBranchAddress("fTreeVariableInvMassLambda",&lInvariantMass);
  if ( fWhichParticle == "AntiLambda"  ) 	lTreeMC->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMass);
  if ( fWhichParticle == "K0Short"     ) 	lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMass);
  if ( fWhichParticle == "Lambda"      ){ //For symmetry of computation...
    lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
    lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingTwo);
  }
  if ( fWhichParticle == "AntiLambda"  ){ //For symmetry of computation...
    lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
    lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingTwo);
  }
  if ( fWhichParticle == "K0Short"     ){
    lTreeMC->SetBranchAddress("fTreeVariableInvMassLambda"    ,&lInvariantMassCompetingOne);
    lTreeMC->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMassCompetingTwo);
  }
  //if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
  //	lTreeMC->SetBranchAddress("fTreeVariableRapLambda",&lRap);
  //if ( fWhichParticle == "K0Short" )
  //	lTreeMC->SetBranchAddress("fTreeVariableRapK0Short",&lRap);

  lTreeMC->SetBranchAddress("fTreeVariableRapMC",&lRap);
	lTreeMC->SetBranchAddress("fTreeVariableLeastNbrCrossedRows",&lLeastNbrCrossedRows);
	lTreeMC->SetBranchAddress("fTreeVariableLeastRatioCrossedRowsOverFindable",&lLeastNbrCrossedRowsOverFindable);
  //--- Topological selection variables -----------------------------
	lTreeMC->SetBranchAddress("fTreeVariableV0Radius",&lV0Radius);
	lTreeMC->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",&lDcaNegToPrimVertex);
	lTreeMC->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",&lDcaPosToPrimVertex);
	lTreeMC->SetBranchAddress("fTreeVariableDcaV0Daughters",&lDcaV0Daughters);
	lTreeMC->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",&lV0CosinePointingAngle);
	lTreeMC->SetBranchAddress("fTreeVariableDistOverTotMom",&lDistOverTotMom);
  //---- Extra  Information Available in MC only --------------------
	lTreeMC->SetBranchAddress("fTreeVariableNegTransvMomentumMC",&lNegTransvMomentumMC);
	lTreeMC->SetBranchAddress("fTreeVariablePosTransvMomentumMC",&lPosTransvMomentumMC);
	lTreeMC->SetBranchAddress("fTreeVariablePID",&lPID);
	lTreeMC->SetBranchAddress("fTreeVariablePIDMother",&lPIDMother);
	lTreeMC->SetBranchAddress("fTreeVariablePIDPositive",&lPIDPositive);
	lTreeMC->SetBranchAddress("fTreeVariablePIDNegative",&lPIDNegative);
	//================================================================

  //================================================================
  cout<<endl;
  cout<<"--------------- MC Data File Loop  ---------------------"<<endl;
  Long_t lOneTenthOfNCandidatesMC = ((double)(lNCandidatesMC) / 10. );
  for(Long_t icand = 0;icand<lNCandidatesMC;icand++){
		lTreeMC->GetEntry(icand);
    if( icand % lOneTenthOfNCandidatesMC == 0 ) 
    	cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidatesMC<<" ( "<<(long)(((double)(icand)/(double)(lNCandidatesMC))*(100.+1e-3))<<"% )"<<endl;
		if(TMath::Abs(lRap)<fRapidityBoundary &&
            TMath::Abs(lNegEta)       <= fCutDaughterEta               &&                   
            TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
            lV0Radius                 >= fCutV0Radius                  &&
            lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
            lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
            lDcaV0Daughters           <= fCutDCAV0Daughters            &&
            lV0CosinePointingAngle    >= fCutV0CosPA                   && 
            lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
            lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
            lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
            TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
            TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
          ( //perfect PID association, IsPhysicalPrimary association
          ( fWhichParticle == "Lambda"
          && lPID           == 3122 //V0 is a Lambda
          && lPIDPositive   == 2212 //Pos Daughter is p
          && lPIDNegative   == -211 //Neg Daughter is pi-
          && lPrimaryStatus == 1
          ) || 
          ( fWhichParticle == "AntiLambda"
          && lPID           == -3122  //V0 is an AntiLambda
          && lPIDPositive   == 211    //Pos Daughter is pi+
          && lPIDNegative   == -2212  //Neg Daughter is antiproton
          && lPrimaryStatus == 1 
          ) || 
          ( fWhichParticle == "K0Short"
          && lPID           == 310   //V0 is an AntiLambda
          && lPIDPositive   == 211   //Pos Daughter is pi+
          && lPIDNegative   == -211  //Neg Daughter is pi-
          && lPrimaryStatus == 1     //K0Short is a primary          
          )
          )
      ){ // Start Entry Loop
		  lWeAreAtBin = fHistPt->FindBin(lPtMC)-1;
		  if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment
		  lHistoFullV0MC[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass
        f2dHistPtResolution->Fill(lPt,lPtMC);
      //Extract left and right background areas and peak 
      //--- Left Background Sample
			if(lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]    && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]   ){ lLeftPlusRightBgV0MC[lWeAreAtBin]++; }
			if(lInvariantMass>lRightBgLeftLimit[lWeAreAtBin]   && lInvariantMass<lRightBgRightLimit[lWeAreAtBin]  ){ lLeftPlusRightBgV0MC[lWeAreAtBin]++; }
      //--- Peak Region
			if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin]      && lInvariantMass<lPeakRightLimit[lWeAreAtBin]     ){ lSigPlusCenterBgV0MC[lWeAreAtBin]++; }
      //--- Info may be needed for geant3/fluka 
      if ( fWhichParticle == "Lambda"     ) lProtonMomentum[lWeAreAtBin]->Fill( lPosTransvMomentumMC );
      if ( fWhichParticle == "AntiLambda" ) lProtonMomentum[lWeAreAtBin]->Fill( lNegTransvMomentumMC );
		} // End Entry Loop
	}
  cout<<"--------------- Loop Completed -------------------------"<<endl;
  cout<<endl;

  cout<<"--------------- X-Check Peak Finding (gauss+linear) ----"<<endl;
 //== Variables for holding peak position, width ==============
  Double_t lPeakPositionMC[100];
  Double_t lPeakWidthMC[100]; 

  char fgausnameMC[100]; 
  TF1 *fgausPtMC[100];
  for(Int_t ibin = 0; ibin<fptbinnumb; ibin++){ 
    cout<<"---> Peak Finding, bin #"<<ibin<<" (perfect count = "<<lHistoFullV0MC[ibin]->GetEntries()<<")..."<<endl;
    sprintf(fgausnameMC,"fGausPtMC%i",ibin);
   	if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){
      fgausPtMC[ibin]= new TF1(fgausnameMC,"[0]*TMath::Gaus(x,[1],[2])",1.116-0.040,1.116+0.040);
  	  fgausPtMC[ibin]->SetParameter(1,1.116);
  	  fgausPtMC[ibin]->SetParameter(2,0.0025); 
  	  fgausPtMC[ibin]->SetParLimits(1,1.1,1.2);
  	  fgausPtMC[ibin]->SetParLimits(2,0.001,0.02);
    }
   	if ( fWhichParticle == "K0Short"){ 
      fgausPtMC[ibin]= new TF1(fgausnameMC,"[0]*TMath::Gaus(x,[1],[2])",0.498-0.06,0.498+0.060);
  	  fgausPtMC[ibin]->SetParameter(1,0.498);
  	  fgausPtMC[ibin]->SetParameter(2,0.004); 
    }
    fgausPtMC[ibin]->SetParameter(0,lHistoFullV0MC[ibin]->GetMaximum() * 0.9);
	  //fgausPtMC[ibin]->SetParameter(3,0); 
	  //fgausPtMC[ibin]->SetParameter(4,HistoFullV0MC[ibin]->GetMaximum() * 0.1); 
    lHistoFullV0MC[ibin]->Fit(fgausnameMC,"QREM0");
    lPeakPositionMC[ibin] = fgausPtMC[ibin]->GetParameter(1);
    lPeakWidthMC[ibin] = fgausPtMC[ibin]->GetParameter(2);
    cout<<"---> ["<<fptbinlimits[ibin]<<" - "<<fptbinlimits[ibin+1]<<" GeV/c]\tPeak at: "<<lPeakPositionMC[ibin]<<", sigma = "<<lPeakWidthMC[ibin]<<endl;
    fHistPeakPositionMC->SetBinContent(ibin+1, lPeakPositionMC[ibin]);
    fHistPeakPositionMC->SetBinError(ibin+1, lPeakWidthMC[ibin]);

    //Create appropriate TLine Objects for Canvases 
    lLineLeftMostMC[ibin]  = new TLine( lLeftBgLeftLimit[ibin],   0, lLeftBgLeftLimit[ibin],   lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );
    lLineLeftMC[ibin]      = new TLine( lLeftBgRightLimit[ibin],  0, lLeftBgRightLimit[ibin],  lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );
    lLineRightMC[ibin]     = new TLine( lRightBgLeftLimit[ibin],  0, lRightBgLeftLimit[ibin],  lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );
    lLineRightMostMC[ibin] = new TLine( lRightBgRightLimit[ibin], 0, lRightBgRightLimit[ibin], lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );

    //Preparing Canvas for storing...
    lCanvasHistoFullV0MC[ibin]->cd();
    lHistoFullV0MC[ibin]->Draw();
    lLineLeftMostMC[ibin]->Draw();
    lLineLeftMC[ibin]->Draw();
    lLineRightMC[ibin]->Draw();
    lLineRightMostMC[ibin]->Draw();
    fgausPtMC[ibin]->Draw("same");

  }
  cout<<"--------------- Peak Finding Finished (in MC) ----------"<<endl;
  cout<<endl;

  TF1 *lfitNoiseMC[50];
  TF1 *lSampleNoiseMC[50];
  char lFitNameOneMC[50];

  if( fFitBackgroundSwitch ){ 
    cout<<"--------------- Backgrounds: fitting with linear -------"<<endl;
    for(long i=0; i<fptbinnumb; i++){
      //Define Function to Fit Background 
      sprintf(lFitNameOneMC,"lfitNoiseMC%i",(int)i);
      cout<<"creating fitnoise, named "<<lFitNameOneMC<<endl;
      lfitNoiseMC[i] = new TF1(lFitNameOneMC, this, &AliV0Module::MyBgPol1, 
        RoundToThousandth ( lLeftBgLeftLimit[i]   ),
        RoundToThousandth ( lRightBgRightLimit[i] ), 4 , "AliV0Module", "MyBgPol1");
      lfitNoiseMC[i] -> FixParameter(2,RoundToThousandth ( lLeftBgRightLimit[i] ) );
      lfitNoiseMC[i] -> FixParameter(3,RoundToThousandth ( lRightBgLeftLimit[i] ) );
      lfitNoiseMC[i] -> SetParameter(0,lLeftPlusRightBgV0MC[i]*lHistoFullV0MC[i]->GetMaximum() / (lSigPlusCenterBgV0MC[i]+1e-6) );
      lfitNoiseMC[i] -> SetParameter(1,0);
      sprintf(lFitNameOneMC,"lSampleNoiseMC%i",(int)i);
      
      //Define Function to Sample Background
      cout<<"creating sample "<<i<<endl;
      lSampleNoiseMC[i] = new TF1(lFitNameOneMC, this, &AliV0Module::MyBgPolToEval1, 
        RoundToThousandth ( lLeftBgLeftLimit[i]   ),
        RoundToThousandth ( lRightBgRightLimit[i] ), 2, "AliV0Module", "MyBgPolToEval1" );
    }    
    for(long i=0; i<fptbinnumb; i++){
      cout<<"Fitting function for bin "<<i<<", get name = "<<lfitNoiseMC[i]->GetName()<<endl;
      sprintf(lFitNameOneMC,"lfitNoiseMC%i",(int)i);
      lHistoFullV0MC[i] -> Fit(lFitNameOneMC,"LLrie+0");
      lSampleNoiseMC[i]->SetParameter(0, lfitNoiseMC[i]->GetParameter(0) );
      lSampleNoiseMC[i]->SetParameter(1, lfitNoiseMC[i]->GetParameter(1) );
    }
    for(long i=0; i<fptbinnumb; i++){
      cout<<"Overriding Background info: Was "<<lLeftPlusRightBgV0MC[i]<<", is now "<<lSampleNoiseMC[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0MC[i]->GetBinWidth(5)<<endl;
      lLeftPlusRightBgV0MC[i] = lSampleNoiseMC[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0MC[i]->GetBinWidth(5);
    }
    cout<<"--------------- Fitting Finished! ----------------------"<<endl;
  }


  //=============================================================  
  // Compute Signal 
  for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
      lSigMCV0[ipoint] = lSigPlusCenterBgV0MC[ipoint] - lLeftPlusRightBgV0MC[ipoint];
      lSigErrMCV0[ipoint] = TMath::Sqrt(lSigPlusCenterBgV0MC[ipoint]+lLeftPlusRightBgV0MC[ipoint]);
      //Error Formula: Equivalent to Sqrt(S+B+B) = Sqrt(S+2B)

      //Signal-to-noise computation
      if( lLeftPlusRightBgV0MC[ipoint] != 0 ){
        fHistSigToNoiseMC->SetBinContent(ipoint+1, (lSigPlusCenterBgV0MC[ipoint] - lLeftPlusRightBgV0MC[ipoint])/lLeftPlusRightBgV0MC[ipoint] );
      }else{
        fHistSigToNoiseMC->SetBinContent(ipoint+1, -1) ; //-1 means: no background
      }
  }
  //=============================================================
  cout<<"--------------- Extracted Signal (MC) ------------------"<<endl;
  for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
      cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<lSigMCV0[ipoint]<<" +/- "<<lSigErrMCV0[ipoint]<<endl;
  }
  cout<<"--------------------------------------------------------"<<endl;
  cout<<endl;


  cout<<"--------------- Process Generated V0s ------------------"<<endl;
	//Filling histogram with original MC particles====================
  TH3F* f3dHistGenPtVsYVsMultV0 = 0x0;
 	  if(fWhichParticle == "Lambda") 
      f3dHistGenPtVsYVsMultV0 			= (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultLambda");
	  if(fWhichParticle == "AntiLambda") 
      f3dHistGenPtVsYVsMultV0 			= (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultAntiLambda");
	  if(fWhichParticle == "K0Short") 
      f3dHistGenPtVsYVsMultV0 			= (TH3F*)v0listMC->FindObject("f3dHistPrimRawPtVsYVsMultK0Short");
	
	TH1D *fHistDummyV0 = f3dHistGenPtVsYVsMultV0->ProjectionX("fHistDummyV0",
		f3dHistGenPtVsYVsMultV0->GetYaxis()->FindBin(-fRapidityBoundary+1e-2),  //avoid taking previous bin
		f3dHistGenPtVsYVsMultV0->GetYaxis()->FindBin(+fRapidityBoundary-1e-2), //avoid taking next bin
		f3dHistGenPtVsYVsMultV0->GetZaxis()->FindBin(-1), //Not interested in Multiplicity so far, integrate
		f3dHistGenPtVsYVsMultV0->GetZaxis()->FindBin(101) //Not interested in Multiplicity so far, integrate
	);

	TH1F *fHistMCCountbyptV0	= new TH1F("fHistMCCountbyptV0","V0 MC count;p_{T} (GeV/c);Counts",fptbinnumb,fptbinlimits);

	Double_t temppt;
	for(long i = 1; i<fHistDummyV0->GetNbinsX()+1;i++){
		temppt = fHistDummyV0->GetXaxis()->GetBinCenter(i);
		for(long filling = 0; filling<fHistDummyV0->GetBinContent(i);filling++){
			fHistMCCountbyptV0->Fill(temppt);
		}	
	}
  cout<<"--------------- Generated V0 Dump ----------------------"<<endl;
  for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
      cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<fHistMCCountbyptV0->GetBinContent(ipoint+1)<<" +/- "<<TMath::Sqrt(fHistMCCountbyptV0->GetBinContent(ipoint+1))<<endl;
  }
  cout<<"--------------------------------------------------------"<<endl;

    //=============================================================  
  // Compute Efficiency 
  Double_t lEfficiency[100];
  Double_t lEfficiencyError[100]; 
  for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
    lEfficiency[ipoint] = lSigMCV0[ipoint] / fHistMCCountbyptV0->GetBinContent(ipoint+1);
    lEfficiencyError[ipoint] = ErrorInRatio(
      lSigMCV0[ipoint],
      lSigErrMCV0[ipoint],
      fHistMCCountbyptV0->GetBinContent(ipoint+1),
      TMath::Sqrt(fHistMCCountbyptV0->GetBinContent(ipoint+1))
    );
  }
  //=============================================================
  cout<<endl;
  cout<<"--------------- Memory Cleanup -------------------------"<<endl;
  fHistDummyV0->Delete();
  fHistMCCountbyptV0->Delete();
  v0listMC->Delete();
  delete v0listMC;
  fileMC->Close("R");
  fileMC->Delete();
  delete fileMC; 
  cout<<endl;
  cout<<"--------------- Pure Efficiency Numbers ----------------"<<endl;
	TH1F* fHistPureEfficiency 		= new TH1F("fHistPureEfficiency","Pure Efficiency;p_{T} (GeV/c);Pure Efficiency",fptbinnumb,fptbinlimits);
	TH1F* fHistEfficiency 		    = new TH1F("fHistEfficiency","Efficiency;p_{T} (GeV/c);Efficiency",fptbinnumb,fptbinlimits);
  if( fWhichParticle == "K0Short"      ) fHistPureEfficiency->SetTitle("K^{0}_{S} Efficiency");
  if( fWhichParticle == "Lambda"       ) fHistPureEfficiency->SetTitle("#Lambda Efficiency (no G3/F corr.)");
  if( fWhichParticle == "AntiLambda"   ) fHistPureEfficiency->SetTitle("#bar{#Lambda} Efficiency (no G3/F corr.)");

  if( fWhichParticle == "Lambda"       ) fHistEfficiency->SetTitle("#Lambda Efficiency (with G3/F corr.)");
  if( fWhichParticle == "AntiLambda"   ) fHistEfficiency->SetTitle("#bar{#Lambda} Efficiency (with G3/F corr.)");

  for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
    cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiency[ipoint]<<" +/- "<<lEfficiencyError[ipoint]<<endl;
    fHistPureEfficiency->SetBinContent(ipoint+1, lEfficiency[ipoint]);
    fHistPureEfficiency->SetBinError(ipoint+1, lEfficiencyError[ipoint]);
  }
  cout<<"--------------------------------------------------------"<<endl;


  cout<<endl;
  //If it's Lambda or AntiLambda, do Geant3/fluka correction
  //Keep fit function outside 'if' scope
  TF1 *fitGeant3FlukaCorr = 0x0;
  TH1D* h1dPanosCorrections =0x0;
  if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){ 
    cout<<"--------------- Geant3/Fluka Correction ----------------"<<endl;

    //Credit for the original version of this code: Antonin Maire (thanks!)    

			TFile *fileCorrGeanT3FlukaPanos = new TFile( "aCorrectionForPPbarCrossSections_FromPanos.root", "READ");
      TH2D* h2dCorrForCrossSectionProtons      = (TH2D*)  ( fileCorrGeanT3FlukaPanos->Get("gHistCorrectionForCrossSectionProtons")      );
      TH2D* h2dCorrCrossSectionAntiProtons     = (TH2D*)  ( fileCorrGeanT3FlukaPanos->Get("gHistCorrectionForCrossSectionAntiProtons")  );
			//Find bins...

      //TH2D* h2dAsMCPtBaryonVsPtCasc  = 0x0;
      TH2D* h2dCorrCrossSection      = 0x0;

			if( fWhichParticle == "Lambda" )     h2dCorrCrossSection   = h2dCorrForCrossSectionProtons;  
      if( fWhichParticle == "AntiLambda" ) h2dCorrCrossSection   = h2dCorrCrossSectionAntiProtons; 

      Int_t myFirstYBinPanos = h2dCorrCrossSection ->GetXaxis()->FindBin(-fRapidityBoundary);
      Int_t myLastYBinPanos  = h2dCorrCrossSection ->GetXaxis()->FindBin(+fRapidityBoundary-1e-3);

			cout<<"Check limiting bins: "<<endl;		
        Printf(" - Bin %d = %f to %f", myFirstYBinPanos, h2dCorrCrossSection ->GetXaxis()->GetBinLowEdge(myFirstYBinPanos), h2dCorrCrossSection ->GetXaxis()->GetBinUpEdge(myFirstYBinPanos) );
        Printf(" - Bin %d = %f to %f", myLastYBinPanos,  h2dCorrCrossSection ->GetXaxis()->GetBinLowEdge(myLastYBinPanos ), h2dCorrCrossSection ->GetXaxis()->GetBinUpEdge(myLastYBinPanos) );

        // B.1.a.
        h1dPanosCorrections = h2dCorrCrossSection->ProjectionY("h1dPanosCorrections", myFirstYBinPanos, myLastYBinPanos, "e" ); // ~mid-rapidity
        h1dPanosCorrections->SetDirectory(0); //means manual cleanup later
        h1dPanosCorrections->Scale( (Double_t) (1.0/(myLastYBinPanos - myFirstYBinPanos +1.0)) );
        h1dPanosCorrections->SetYTitle("#epsilon_{GEANT3} / #epsilon_{FLUKA}");

        //TF1 *fitGeant3FlukaCorr = 0x0; declared outside this scope
        if(fWhichParticle == "Lambda") fitGeant3FlukaCorr = new TF1("fitGeant3FlukaCorr", this, &AliV0Module::MyGeant3FlukaCorrectionForProtons,     0.25, 1.1, 3, "AliV0Module","MyGeant3FlukaCorrectionForProtons");
        if(fWhichParticle == "AntiLambda") fitGeant3FlukaCorr = new TF1("fitGeant3FlukaCorr", this, &AliV0Module::MyGeant3FlukaCorrectionForAntiProtons, 0.25, 1.1, 4, "AliV0Module", "MyGeant3FlukaCorrectionForAntiProtons");

        h1dPanosCorrections->Fit("fitGeant3FlukaCorr","rime+0"); // Chi2

        Printf("Test fit function...");
        Printf(" - fitGeant3FlukaCorr for pt =  .25 GeV/c : %f", fitGeant3FlukaCorr->Eval( .25) );
        Printf(" - fitGeant3FlukaCorr for pt =  .5  GeV/c : %f", fitGeant3FlukaCorr->Eval( .5) );
        Printf(" - fitGeant3FlukaCorr for pt =  1   GeV/c : %f", fitGeant3FlukaCorr->Eval( 1.0) );
        Printf(" - fitGeant3FlukaCorr for pt =  2   GeV/c : %f", fitGeant3FlukaCorr->Eval( 2.0) );
        Printf(" - fitGeant3FlukaCorr for pt =  5   GeV/c : %f", fitGeant3FlukaCorr->Eval( 5.0) );
        Printf(" - fitGeant3FlukaCorr for pt =  7   GeV/c : %f", fitGeant3FlukaCorr->Eval( 7.0) );
        Printf(" - fitGeant3FlukaCorr for pt =  8   GeV/c : %f", fitGeant3FlukaCorr->Eval( 8.0) );
        Printf(" - fitGeant3FlukaCorr for pt = 10   GeV/c : %f", fitGeant3FlukaCorr->Eval(10.0) );
    cout<<"--------------- Embedding G3/F in Efficiencies ---------"<<endl;
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
		  lEfficiency[ipoint]         /= fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint]->GetMean() ) ;
      lEfficiencyError[ipoint]    /= fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint]->GetMean() ) ;
    }
    cout<<"--------------- G3/F Corrected Efficiencies ------------"<<endl;
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
      cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiency[ipoint]<<" +/- "<<lEfficiencyError[ipoint]<<endl;
      fHistEfficiency->SetBinContent(ipoint+1, lEfficiency[ipoint]);
      fHistEfficiency->SetBinError(ipoint+1, lEfficiencyError[ipoint]);
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;

    fileCorrGeanT3FlukaPanos->Close();
    fileCorrGeanT3FlukaPanos->Delete();
    delete fileCorrGeanT3FlukaPanos;

  } //end Geant3/fluka correction if

  //Declaring outside Feeddown scope to keep in memory
  TH2F *f2dFeedDownMatrix = 0x0;
  TH2F *f2dFeedDownEfficiency = 0x0;
  TH2F *f2dFeedDownEfficiencyGFCorrected = 0x0;
  TH1F *fHistFeeddownSubtraction = 0x0;

  //=========================================================================
  // --- Feeddown correction section 
  if( fFDSwitch != "NoFD" && fWhichParticle != "K0Short"){ 
    cout<<"--------------- Feeddown Correction --------------------"<<endl;
	  //MC File Acquisition=============================================
    cout<<"--------------- Opening MC file for Feeddown -----------"<<endl;
	  TFile* fileMCFD = TFile::Open(fFeedDownDataFile, "READ");
	  fileMCFD->cd("PWG2CheckPerformanceLambda_PP_MC");
	  TList* v0listMCFD  = (TList*)fileMCFD->FindObjectAny("clistV0MC");
    TTree* lTreeMCFD;
    TH1F* fHistMultiplicityBeforeTrigSelMCFD;
    lTreeMCFD = (TTree*)fileMCFD->FindObjectAny("fTree");
    fHistMultiplicityBeforeTrigSelMCFD = (TH1F*)v0listMCFD->FindObject("fHistMultiplicityBeforeTrigSel");
 
	  Long_t lNCandidatesMCFD = lTreeMCFD->GetEntries();
	  Long_t lNEventsMCFD = fHistMultiplicityBeforeTrigSelMCFD->GetEntries();
	  cout<<" MC Events (before trig sel)...: "<<lNEventsMCFD<<endl;
	  cout<<" Cascade MC Candidates.........: "<<lNCandidatesMCFD<<endl;
    cout<<"--------------------------------------------------------"<<endl;
	  //================================================================

    //Linking to Tree=================================================
    //--- Base Variables ----------------------------------------------
	  lTreeMCFD->SetBranchAddress("fTreeVariablePosEta",&lPosEta);
	  lTreeMCFD->SetBranchAddress("fTreeVariableNegEta",&lNegEta);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePrimaryStatus",&lPrimaryStatus);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePrimaryStatusMother",&lPrimaryStatusMother);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePt",&lPt);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePtMC",&lPtMC);
  	lTreeMCFD->SetBranchAddress("fTreeVariablePtXiMother",&lPtMother);
    if ( fWhichParticle == "Lambda"      )  lTreeMCFD->SetBranchAddress("fTreeVariableInvMassLambda",&lInvariantMass);
    if ( fWhichParticle == "AntiLambda"  ) 	lTreeMCFD->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMass);
    if ( fWhichParticle == "K0Short"     ) 	lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMass);
    if ( fWhichParticle == "Lambda"      ){ //For symmetry of computation...
      lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
      lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingTwo);
    }
    if ( fWhichParticle == "AntiLambda"  ){ //For symmetry of computation...
      lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
      lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingTwo);
    }
    if ( fWhichParticle == "K0Short"     ){
      lTreeMCFD->SetBranchAddress("fTreeVariableInvMassLambda"    ,&lInvariantMassCompetingOne);
      lTreeMCFD->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMassCompetingTwo);
    }
    //if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
    //	lTreeMCFD->SetBranchAddress("fTreeVariableRapLambda",&lRap);
    //if ( fWhichParticle == "K0Short" )
    //	lTreeMCFD->SetBranchAddress("fTreeVariableRapK0Short",&lRap);

    lTreeMCFD->SetBranchAddress("fTreeVariableRapMC",&lRap);
	  lTreeMCFD->SetBranchAddress("fTreeVariableLeastNbrCrossedRows",&lLeastNbrCrossedRows);
	  lTreeMCFD->SetBranchAddress("fTreeVariableLeastRatioCrossedRowsOverFindable",&lLeastNbrCrossedRowsOverFindable);
    //--- Topological selection variables -----------------------------
	  lTreeMCFD->SetBranchAddress("fTreeVariableV0Radius",&lV0Radius);
	  lTreeMCFD->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",&lDcaNegToPrimVertex);
	  lTreeMCFD->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",&lDcaPosToPrimVertex);
	  lTreeMCFD->SetBranchAddress("fTreeVariableDcaV0Daughters",&lDcaV0Daughters);
	  lTreeMCFD->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",&lV0CosinePointingAngle);
	  lTreeMCFD->SetBranchAddress("fTreeVariableDistOverTotMom",&lDistOverTotMom);
    //---- Extra  Information Available in MC only --------------------
	  lTreeMCFD->SetBranchAddress("fTreeVariableNegTransvMomentumMC",&lNegTransvMomentumMC);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePosTransvMomentumMC",&lPosTransvMomentumMC);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePID",&lPID);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePIDMother",&lPIDMother);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePIDPositive",&lPIDPositive);
	  lTreeMCFD->SetBranchAddress("fTreeVariablePIDNegative",&lPIDNegative);
	  //================================================================

    //================================================================
    // Define Feeddown matrix 
    Double_t lFeedDownMatrix [100][100]; 
    //       FeedDownMatrix [Lambda Bin][Xi Bin];
    for(Int_t ilb = 0; ilb<100; ilb++){
      for(Int_t ixb = 0; ixb<100; ixb++){ 
        lFeedDownMatrix[ilb][ixb]=0;
      }
    }
    //Define Xi Binning 
	  Double_t xibinlimits[] = {
       0.00,  0.20,  0.40,  0.60,  0.80,  0.90,
       1.00,  1.10,  1.20,  1.30,  1.40,  1.50,
       1.70,  1.90,  2.20,  2.60,  3.10,  3.90,
       4.90,  6.00,  7.20,  8.50 ,10.00, 12.00 }; 
    Long_t xibinnumb = sizeof(xibinlimits)/sizeof(Double_t) - 1;

  	TH1F* fHistPtXiReference =
      new TH1F("fHistPtXiReference","#Xi candidates uncorrected p_{T};p_{T} (GeV/c);Counts",xibinnumb,xibinlimits);
  
    //Feeddown: will use a double-coordinate system: 
    // ( lambda bin, xi bin ) all the time! 
    lWeAreAtBin = 0;         // Lambda bin 
    Int_t lWeAreAtXiBin = 0; // Xi Bin 

    cout<<"--------------- MC Data File Loop, Feeddown  -----------"<<endl;
    Long_t lOneTenthOfNCandidatesMCFD = ((double)(lNCandidatesMCFD) / 10. );
    for(Long_t icand = 0;icand<lNCandidatesMCFD;icand++){
		  lTreeMCFD->GetEntry(icand);
      if( icand % lOneTenthOfNCandidatesMCFD == 0 ) 
      	cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidatesMCFD<<" ( "<<(long)(((double)(icand)/(double)(lNCandidatesMCFD))*(100.+1e-3))<<"% )"<<endl;
		  if(TMath::Abs(lRap)<fRapidityBoundary &&
              TMath::Abs(lNegEta)       <= fCutDaughterEta               &&                   
              TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
              lV0Radius                 >= fCutV0Radius                  &&
              lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
              lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
              lDcaV0Daughters           <= fCutDCAV0Daughters            &&
              lV0CosinePointingAngle    >= fCutV0CosPA                   && 
              lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
              lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
              lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
              TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
              TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
            ( //perfect PID association, IsPhysicalPrimary association
            ( fWhichParticle == "Lambda"
            && lPID           == 3122 //V0 is a Lambda
            && lPIDPositive   == 2212 //Pos Daughter is p
            && lPIDNegative   == -211 //Neg Daughter is pi-
            && (lPIDMother     == 3312 || (lPIDMother == 3322 && fFDSwitch=="UseMCRatio") ) 
            && lPrimaryStatusMother == 1 //Xi is actually a primary (should matter little)
            ) || 
            ( fWhichParticle == "AntiLambda"
            && lPID           == -3122  //V0 is an AntiLambda
            && lPIDPositive   == 211    //Pos Daughter is pi+
            && lPIDNegative   == -2212  //Neg Daughter is antiproton
            && (lPIDMother     ==-3312 || (lPIDMother ==-3322 && fFDSwitch=="UseMCRatio") ) 
            && lPrimaryStatusMother == 1 //Xi is actually a primary (should matter little) 
            )
            )
        ){ // Start Entry Loop
		    lWeAreAtBin   = fHistPt->FindBin(lPtMC)-1;
		    if(lWeAreAtBin == -1) lWeAreAtBin   = 99; //UnderFlow, special treatment
        lWeAreAtXiBin = fHistPtXiReference->FindBin(lPtMother)-1; 
		    if(lWeAreAtXiBin == -1) lWeAreAtXiBin = 99; //UnderFlow, special treatment
        //Populate Feeddown Matrix
        //cout<<" Populate at coordinates "<<WeAreAtBin<<" vs "<<WeAreAtXiBin<<" ....."<<endl;
        lFeedDownMatrix[lWeAreAtBin][lWeAreAtXiBin]++;
		  } // End Entry Loop
	  }
    cout<<"--------------- Loop Completed -------------------------"<<endl;
    cout<<endl;


	  //Filling histogram with original MC particles====================
    TH3F* f3dHistGenPtVsYVsMultV0FeedDown = 0x0;

	    if (fWhichParticle == "Lambda" )
        f3dHistGenPtVsYVsMultV0FeedDown = (TH3F*)v0listMCFD->FindObject("f3dHistGenPtVsYVsMultXiMinus");
	    if (fWhichParticle == "AntiLambda" )
        f3dHistGenPtVsYVsMultV0FeedDown = (TH3F*)v0listMCFD->FindObject("f3dHistGenPtVsYVsMultXiPlus");
	
	  TH1D *fHistDummyV0FeedDown = f3dHistGenPtVsYVsMultV0FeedDown->ProjectionX("fHistDummyV0FeedDown",
		  f3dHistGenPtVsYVsMultV0FeedDown->GetYaxis()->FindBin(-fRapidityBoundary+1e-2),  //avoid taking previous bin
		  f3dHistGenPtVsYVsMultV0FeedDown->GetYaxis()->FindBin(+fRapidityBoundary-1e-2), //avoid taking next bin
		  f3dHistGenPtVsYVsMultV0FeedDown->GetZaxis()->FindBin(-1),
		  f3dHistGenPtVsYVsMultV0FeedDown->GetZaxis()->FindBin(101)
	  );

	  TH1F *fHistMCCountbyptXiFeedDown	= new TH1F("fHistMCCountbyptXiFeedDown","#Xi MC count;p_{T} (GeV/c);Counts",xibinnumb,xibinlimits);
	  //Double_t temppt, y; --- declared before
    for(long i = 1; i<fHistDummyV0FeedDown->GetNbinsX()+1;i++){
      temppt = fHistDummyV0FeedDown->GetXaxis()->GetBinCenter(i);
      for(long filling = 0; filling<fHistDummyV0FeedDown->GetBinContent(i);filling++){
        fHistMCCountbyptXiFeedDown->Fill(temppt);
      }	
    }
    if(fWhichParticle == "Lambda") 
      cout<<"--------------- Generated Xi- Dump ---------------------"<<endl;
    if(fWhichParticle == "AntiLambda") 
      cout<<"--------------- Generated Xi+ Dump ---------------------"<<endl;
    for(Int_t ipoint = 0; ipoint<xibinnumb; ipoint++){
        cout<<"---> ["<<xibinlimits[ipoint]<<" - "<<xibinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<fHistMCCountbyptXiFeedDown->GetBinContent(ipoint+1)<<" +/- "<<TMath::Sqrt(fHistMCCountbyptXiFeedDown->GetBinContent(ipoint+1))<<endl;
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;
    cout<<"Computing actual feeddown matrix"<<endl;
    Double_t lFeedDownEfficiency[100][100];
    Double_t lFeedDownEfficiencyGFCorrected[100][100];
    Double_t lFeedDownEfficiencyError[100][100];
    Double_t lFeedDownEfficiencyGFCorrectedError[100][100];

    for(Int_t ilb = 0; ilb<fptbinnumb; ilb++){
      for(Int_t ixb = 0; ixb<xibinnumb; ixb++){ 
        if( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1)!= 0 ){ //avoid div by zero
          lFeedDownEfficiency[ilb][ixb]=((double)(lFeedDownMatrix[ilb][ixb])) / ((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) )) ;
          lFeedDownEfficiencyError[ilb][ixb]=ErrorInRatio( 
            ((double)(lFeedDownMatrix[ilb][ixb])),
            TMath::Sqrt((double)(lFeedDownMatrix[ilb][ixb])),
            ((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) )),
            TMath::Sqrt((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) ))
          );
        }else{ 
          lFeedDownEfficiency[ilb][ixb] = 0;
          lFeedDownEfficiencyError[ilb][ixb] = 0;
        }
      }
    }

    //Beware: the feeddown correction efficiency is computed with MC,
    //Which has the geant3/fluka problem. Thus, we actually geant3/fluka
    //correct the Feeddown efficiencies, too...
  	for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
      for(Int_t ixb = 0;ixb<xibinnumb; ixb++){
        lFeedDownEfficiencyGFCorrected[ipoint][ixb] = lFeedDownEfficiency[ipoint][ixb]/fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint]->GetMean() ) ;
        lFeedDownEfficiencyGFCorrectedError[ipoint][ixb] = lFeedDownEfficiencyError[ipoint][ixb]/fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint]->GetMean() ) ;
      } 
    }

    //Create FD TH2Fs for storing
    f2dFeedDownMatrix = new TH2F("f2dFeedDownMatrix","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
    f2dFeedDownEfficiency = new TH2F("f2dFeedDownEfficiency","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
    f2dFeedDownEfficiencyGFCorrected = new TH2F("f2dFeedDownEfficiencyGFCorrected","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
    f2dFeedDownMatrix->SetDirectory(0);
    f2dFeedDownEfficiency->SetDirectory(0);
    f2dFeedDownEfficiencyGFCorrected->SetDirectory(0);
    for(Int_t ilb = 0; ilb<fptbinnumb; ilb++){
      for(Int_t ixb = 0; ixb<xibinnumb; ixb++){ 
        f2dFeedDownMatrix->SetBinContent(ilb+1,ixb+1,lFeedDownMatrix[ilb][ixb]);
        f2dFeedDownEfficiency->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiency[ilb][ixb]);
        f2dFeedDownEfficiency->SetBinError(ilb+1,ixb+1,lFeedDownEfficiencyError[ilb][ixb]);
        f2dFeedDownEfficiencyGFCorrected->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiencyGFCorrected[ilb][ixb]);
        f2dFeedDownEfficiencyGFCorrected->SetBinError(ilb+1,ixb+1,lFeedDownEfficiencyGFCorrectedError[ilb][ixb]);
      }
    }

    Double_t lProducedXi[100]; 
    
	  TF1 *lLevyFitXiFeedDown = new TF1("LevyFitXiFeedDown",this ,&AliV0Module::MyLevyPtXi, 0.0,15,3, "AliV0Module","MyLevyPtXi");
    
    //Set Parameters Adequately, as in paper
    //FIXME: These are the 7TeV Xi- parameters!!! 
    //FIXME: Use points vs fit function (should matter little if Xi- fit is good)

    if(fWhichParticle == "Lambda"){ 
      lLevyFitXiFeedDown->SetParameter(0, 7.98e-03);
      lLevyFitXiFeedDown->SetParameter(1, 3.4429e-1);
      lLevyFitXiFeedDown->SetParameter(2, 1.0787e+1);
    }
    if(fWhichParticle == "AntiLambda"){ 
      lLevyFitXiFeedDown->SetParameter(0, 7.79e-03);
      lLevyFitXiFeedDown->SetParameter(1, 3.3934e-1);
      lLevyFitXiFeedDown->SetParameter(2, 1.0428e+1);
    }
    //If you want me to double charged Xi feeddown, I'll do it here
    if( fFDSwitch == "DoubleChargedXi" ){ 
    lLevyFitXiFeedDown->SetParameter(0, lLevyFitXiFeedDown->GetParameter(0)*2 );
    }
    for(Int_t ixb = 0; ixb<xibinnumb; ixb++){ 
      lProducedXi[ixb] = lNInelasticEvents * lLevyFitXiFeedDown->Integral(xibinlimits[ixb],xibinlimits[ixb+1]); 
    }
    if(fWhichParticle == "Lambda") 
      cout<<"--------------- Generated Xi- Dump (real-corrected) ----"<<endl;
    if(fWhichParticle == "AntiLambda") 
      cout<<"--------------- Generated Xi+ Dump (real-corrected) ----"<<endl;
    for(Int_t ixb = 0; ixb<xibinnumb; ixb++){ 
        cout<<"Xi bin "<<ixb<<"\t"<<lProducedXi[ixb]<<endl;
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;

    //=== Computing actual feeddown numbers... ===
    Double_t lFeedDownToSubtract[100];
    Double_t lFeedDownToSubtractError[100];
    for(Int_t ilb = 0; ilb < fptbinnumb; ilb++){ 
      lFeedDownToSubtract[ilb] = 0; 
      lFeedDownToSubtractError[ilb] = 0; 
      for(Int_t ixb = 0; ixb < xibinnumb; ixb++){ 
        lFeedDownToSubtract[ilb] += lProducedXi[ixb]*lFeedDownEfficiencyGFCorrected[ilb][ixb]; 
        lFeedDownToSubtractError[ilb] += TMath::Power( (lProducedXi[ixb]*lFeedDownEfficiencyGFCorrectedError[ilb][ixb]), 2 ); 
      }
      lFeedDownToSubtractError[ilb] = TMath::Sqrt(lFeedDownToSubtractError[ilb]);
    }
  	fHistFeeddownSubtraction 		= new TH1F("fHistFeeddownSubtraction","FD Subtraction;p_{T} (GeV/c);Ratio removed",fptbinnumb,fptbinlimits);
    fHistFeeddownSubtraction->SetDirectory(0);
    for(Int_t ilb = 0; ilb<fptbinnumb; ilb++){ 
      fHistFeeddownSubtraction->SetBinContent( ilb+1 , ((double)(lFeedDownToSubtract[ilb])) / ((double)(lSigRealV0[ilb]) ) ); 
      fHistFeeddownSubtraction->SetBinError  ( ilb+1 , 
        ErrorInRatio( 
          ((double)(lFeedDownToSubtract[ilb])), 
          ((double)(lFeedDownToSubtractError[ilb])), 
          ((double)(lSigRealV0[ilb]) ), 
          ((double)(lSigErrRealV0[ilb]) ) ) 
        ); 
    }



    cout<<"--------------- FD Subtraction Fraction ----------------"<<endl;
    for(Int_t ilb = 0; ilb<fptbinnumb; ilb++){ 
      cout<<"---> ["<<fptbinlimits[ilb]<<" - "<<fptbinlimits[ilb+1]<<" GeV/c]\t"<<fHistFeeddownSubtraction->GetBinContent(ilb+1)<<"\t+/-\t"<<fHistFeeddownSubtraction->GetBinError(ilb+1)<<endl;
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<" Performing Actual Correction... "<<endl;
    for(Int_t ilb = 0; ilb<fptbinnumb; ilb++){
      lSigRealV0   [ilb] = lSigRealV0[ilb] - lFeedDownToSubtract[ilb];
      lSigErrRealV0[ilb] = TMath::Sqrt( TMath::Power(lSigErrRealV0[ilb],2) + TMath::Power(lFeedDownToSubtractError[ilb],2) );
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;
    cout<<"--------------- Memory Cleanup -------------------------"<<endl;
	  lLevyFitXiFeedDown->Delete();
    v0listMCFD->Delete();
    delete v0listMCFD;
    fileMCFD->Close("R");
    fileMCFD->Delete();
    delete fileMCFD;
    cout<<endl;
  }//--- End Feeddown Correction section
  //=========================================================================

  //=========================================================================
  //---> At this stage, everthing's just ready for the actual spectrum 
  //---> computation to occur! 
  Double_t lSpectrum[100];
  Double_t lSpectrumError[100];

  for(Int_t ibin=0;ibin<fptbinnumb;ibin++){ 
    lSpectrum[ibin] = lSigRealV0[ibin] / lEfficiency [ibin]; 
    lSpectrumError[ibin] = ErrorInRatio( 
      lSigRealV0    [ibin],
      lSigErrRealV0 [ibin],
      lEfficiency      [ibin],
      lEfficiencyError [ibin]);
  }

  //Divide by: Bin Width, Rapidity Window, N_{Inel} 
  for(Int_t ibin=0;ibin<fptbinnumb;ibin++){ 
    lSpectrum[ibin] /= (fptbinlimits[ibin+1]-fptbinlimits[ibin]);
    lSpectrum[ibin] /= lNInelasticEvents; 
    lSpectrum[ibin] /= 2*fRapidityBoundary; 
    lSpectrumError[ibin] /= (fptbinlimits[ibin+1]-fptbinlimits[ibin]);
    lSpectrumError[ibin] /= lNInelasticEvents; 
    lSpectrumError[ibin] /= 2*fRapidityBoundary; 
  }

  TH1F* fHistPtLambda = new TH1F("fHistPtLambda","#Lambda Corrected Spectrum;p_{T} (GeV/c);1/N_{INEL} #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);
  TH1F* fHistPtAntiLambda = new TH1F("fHistPtAntiLambda","#bar{#Lambda} Corrected Spectrum;p_{T} (GeV/c);1/N_{INEL} #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);
  TH1F* fHistPtK0Short = new TH1F("fHistPtK0Short","K^{0}_{S} Corrected Spectrum;p_{T} (GeV/c);1/N_{INEL} #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);

  //Copy to Histogram
  for(Int_t ibin=0;ibin<fptbinnumb;ibin++){ 
    if(fWhichParticle == "Lambda" ){
      fHistPtLambda->SetBinContent( ibin+1, lSpectrum[ibin]      );
      fHistPtLambda->SetBinError  ( ibin+1, lSpectrumError[ibin] );
    }
    if(fWhichParticle == "AntiLambda" ){
      fHistPtAntiLambda->SetBinContent( ibin+1, lSpectrum[ibin]      );
      fHistPtAntiLambda->SetBinError  ( ibin+1, lSpectrumError[ibin] );
    }
    if(fWhichParticle == "K0Short" ){
      fHistPtK0Short->SetBinContent( ibin+1, lSpectrum[ibin]      );
      fHistPtK0Short->SetBinError  ( ibin+1, lSpectrumError[ibin] );
    }    
  }
  

  //=========================================================================

  cout<<"--------------- Result Output --------------------------"<<endl;
  cout<<" ---> Writing information to "<<fOutputDataFile<<endl;
  // Open an output file
  TFile* lResultsFile = TFile::Open(fOutputDataFile, "RECREATE");
  if (!lResultsFile || !lResultsFile->IsOpen()){ 
    cout<<"Error! Couldn't open file!"<<endl;
    return;
  }
  
  //Preparing Signal Extraction Range Canvas 
  TCanvas *cSigExtRange = new TCanvas("cSigExtRange","Extraction Range",900,600);
  cSigExtRange->SetFillColor(kWhite);
  cSigExtRange->SetLeftMargin(0.17);
  cSigExtRange->SetRightMargin(0.17);
  cSigExtRange->cd();
  fHistSignalExtractionRange->SetFillColor(18);
  fHistSignalExtractionRange->SetMarkerStyle(20);
  if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){ 
    fHistSignalExtractionRange->GetYaxis()->SetRangeUser(1.115-0.08,1.115+0.08);
  }
  if( fWhichParticle == "K0Short" ){ 
    fHistSignalExtractionRange->GetYaxis()->SetRangeUser(0.498-0.15,0.498+0.15);
  }
  fHistSignalExtractionRange->Draw("E2");
  if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){ 
    fLDataUpper->Draw("same");
    fLDataLower->Draw("same");
  }
  if( fWhichParticle == "K0Short" ){ 
    fKDataUpper->Draw("same");
    fKDataLower->Draw("same");
  }
  //Saving Invariant Mass Plots (real data)
  lResultsFile->cd();
  TDirectoryFile *lInvMassReal = new TDirectoryFile("lInvMassReal","Invariant Mass Plots (Real Data)");
  lInvMassReal->cd();
  cSigExtRange->Write();
  for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
    lCanvasHistoFullV0[ibin] -> Write();
  }

  TDirectoryFile *lInvMassRealRawData = new TDirectoryFile("lInvMassRealRawData","Objects for Inv Mass Plots (Real Data)");
  lInvMassRealRawData->cd();
  fHistSignalExtractionRange->Write();
  fHistPeakPosition->Write();
  fHistSigToNoise->Write();
  for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
    lHistoFullV0[ibin] -> Write();
    fgausPt[ibin]         -> Write();
    if( fFitBackgroundSwitch ) lfitNoise[ibin] -> Write();
  }
  if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){ 
    fLDataUpper->Write();
    fLDataLower->Write();
  }
  if( fWhichParticle == "K0Short" ){ 
    fKDataUpper->Write();
    fKDataLower->Write();
  }

  //Saving Invariant Mass Plots (MC)
  lResultsFile->cd();
  TDirectoryFile *lInvMassMC = new TDirectoryFile("lInvMassMC","Invariant Mass Plots (Monte Carlo)");
  lInvMassMC->cd();

  for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
    lCanvasHistoFullV0MC[ibin] -> Write();
  }

  TDirectoryFile *lInvMassMCRawData = new TDirectoryFile("lInvMassMCRawData","Objects for Inv Mass Plots (MC)");
  lInvMassMCRawData->cd();
  fHistPeakPositionMC->Write();
  fHistSigToNoiseMC->Write();
   f2dHistPtResolution->Write();
  for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
    lHistoFullV0MC[ibin] ->Write();
    fgausPtMC[ibin]         ->Write();
    if( fFitBackgroundSwitch ) lfitNoiseMC[ibin] -> Write();
  }

  //Saving Geant3/Fluka Correction Data (MC)
  if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){
    lResultsFile->cd();
    TDirectoryFile *lGFCorrection = new TDirectoryFile("lGFCorrection","Geant3/Fluka Correction Histograms");
    lGFCorrection->cd();
    h1dPanosCorrections->Write();
    fitGeant3FlukaCorr->Write();
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
      lProtonMomentum[ibin]  ->Write();
    }
  }

  //Saving Feeddown Correction information, if needed
  if( (fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda")&&fFDSwitch!="NoFD" ){
    lResultsFile->cd();
    TDirectoryFile *lFeeddown = new TDirectoryFile("lFeeddown","Feeddown Subtraction Information");
    lFeeddown->cd();
    f2dFeedDownMatrix->Write();
    f2dFeedDownEfficiency->Write();
    f2dFeedDownEfficiencyGFCorrected->Write();
    fHistFeeddownSubtraction->Write();
  }

  lResultsFile->cd();

  if( fWhichParticle == "K0Short") fHistPureEfficiency->Write();
  if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){
    fHistPureEfficiency->Write();
    fHistEfficiency->Write();
  }
  if(fWhichParticle == "Lambda" ){
    fHistPtLambda     ->  Write();
  }
  if(fWhichParticle == "AntiLambda" ){
    fHistPtAntiLambda -> Write();
  }
  if(fWhichParticle == "K0Short" ){
    fHistPtK0Short    -> Write();
  }    

  lResultsFile->Write();
  lResultsFile->Close();
  delete lResultsFile;

  //================================================
  //Manual Cleanup of all created Histograms 
  //================================================
  // needed if you want to re-run the whole thing without 
  // memory leaks (systematics, etc)
  
  //switch on if you want large amounts of printout
  Bool_t lDebugCleaningProcess = kFALSE; 

  if (lDebugCleaningProcess) cout<<"fHistPt->Delete()"<<endl;
  fHistPt->Delete();
  if (lDebugCleaningProcess) cout<<"fHistPeakPosition->Delete()"<<endl;
  fHistPeakPosition->Delete();
  if (lDebugCleaningProcess) cout<<"fHistPeakPositionMC->Delete()"<<endl;
  fHistPeakPositionMC->Delete();
  if (lDebugCleaningProcess) cout<<"fHistSigToNoise->Delete()"<<endl;
  fHistSigToNoise->Delete();
  if (lDebugCleaningProcess) cout<<"fHistSigToNoiseMC->Delete()"<<endl;
  fHistSigToNoiseMC->Delete();
  if (lDebugCleaningProcess) cout<<"fHistSignalExtractionRange->Delete()"<<endl;
  fHistSignalExtractionRange->Delete();
  if (lDebugCleaningProcess) cout<<"fKData*->Delete()"<<endl;
  fKDataUpper->Delete();
  fKDataLower->Delete();
  fLDataUpper->Delete();
  fLDataLower->Delete(); 
   if(lDebugCleaningProcess) cout<<"f2dHistPtResolution->Delete()"<<endl;   
   f2dHistPtResolution->Delete();

  if(lDebugCleaningProcess) cout<<"lfitNoise*[*]->Delete(); lSampleNoise*[*]->Delete()"<<endl;   
  if( fFitBackgroundSwitch ){ 
    for(long i=0; i<fptbinnumb; i++){
      lfitNoise[i]    -> Delete();
      lSampleNoise[i] -> Delete();
      lfitNoiseMC[i]    -> Delete();
      lSampleNoiseMC[i] -> Delete();
    }
  }


  //pt-by-pt histos
  if (lDebugCleaningProcess) cout<<"lHistoFullV0*[*]->Delete()"<<endl;
  for(Int_t ihist=0;ihist<100;ihist++){
    lHistoFullV0[ihist]->Delete();
    lHistoFullV0MC[ihist]->Delete();
  }

  if (lDebugCleaningProcess) cout<<"lLine*[*]->Delete()"<<endl;
  //gaussian fit functions, drawing lines
  for(Int_t ibin = 0; ibin<fptbinnumb; ibin++){ 
    lLineLeftMost[ibin] ->Delete(); 
    lLineLeft[ibin]     ->Delete();
    lLineRight[ibin]    ->Delete();
    lLineRightMost[ibin]->Delete();
    lLineLeftMostMC[ibin] ->Delete(); 
    lLineLeftMC[ibin]     ->Delete();
    lLineRightMC[ibin]    ->Delete();
    lLineRightMostMC[ibin]->Delete();

   	if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){
      fgausPt[ibin]->Delete();
      fgausPtMC[ibin]->Delete();
    }
   	if ( fWhichParticle == "K0Short"){ 
      fgausPt[ibin]->Delete();
      fgausPtMC[ibin]->Delete();
    }
  }
  if (lDebugCleaningProcess) cout<<"lCanvasHistoFullV0*[*]->Delete()"<<endl;
  for(Int_t ihist=0;ihist<100;ihist++){
    lCanvasHistoFullV0[ihist]->Close();
    lCanvasHistoFullV0MC[ihist]->Close();
    delete lCanvasHistoFullV0[ihist];
    delete lCanvasHistoFullV0MC[ihist];
  }

  if (lDebugCleaningProcess) cout<<"cSigExtRange->Delete()"<<endl;
  cSigExtRange->Close();  
  delete cSigExtRange;

  if (lDebugCleaningProcess) cout<<"lProtonMomentum[*]->Delete()"<<endl;
  //histograms for G3/F correction
  for(Int_t ibin=0;ibin<100;ibin++){ 
    lProtonMomentum[ibin]->Delete();
  }

  if (lDebugCleaningProcess) cout<<"fHist*Efficiency->Delete()"<<endl;
  //MC info, efficiencies
	fHistPureEfficiency->Delete();
	fHistEfficiency->Delete();
  //data for G3/F, continued
  if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ){ 
      h1dPanosCorrections->Delete();
      fitGeant3FlukaCorr->Delete();
  }

  //histograms, feeddown 
  if (lDebugCleaningProcess) cout<<"f2dFeedDown*->Delete()"<<endl;
  if( fFDSwitch != "NoFD" && fWhichParticle != "K0Short"){ 
    //fHistPtXiReference->Delete();
    //fHistMCCountbyptXiFeedDown->Delete();
    f2dFeedDownMatrix->Delete();
    f2dFeedDownEfficiency->Delete();
    f2dFeedDownEfficiencyGFCorrected->Delete();
    fHistFeeddownSubtraction->Delete();
  }
  if (lDebugCleaningProcess) cout<<"fHistPt*->Delete()"<<endl;
  //Corrected Spectra Histograms
  fHistPtLambda->Delete();
  fHistPtAntiLambda->Delete();
  fHistPtK0Short->Delete();

  //Exit Batch Mode
  gROOT->SetBatch (kFALSE);

  cout<<"--------------------------------------------------------"<<endl;
  cout<<" There, done! "<<endl;
  cout<<endl;




}
