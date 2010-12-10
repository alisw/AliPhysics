// author: Ionut Cristian Arsene
// Date:  01/09/2010

#include <iostream>
#include <fstream>
using namespace std;

#include <TObjArray.h>
#include <TNamed.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>

#include "AliCFContainer.h"
#include "AliDielectronCFdraw.h"

/*
This macro makes projections and saves histograms from a list of CF containers
generated with the dielectron package. These histograms can later be used to
calculate efficiencies.
To use it, the following modifications are needed:
1) Modify the global variables listed below according to your needs.
2) Make projections by applying as many cut sets as needed on the CF containers.
   Call the FillHistograms() after each cut set.
3) To extract efficiencies use the ExtractEfficiencies() function where one needs to 
   specify the indexes for the nominator and denominator histograms. The
   indexes are based on the CF step number, cut set number and histogram number
   as defined in the global variables below.
*/

// The CF variable indexes ------------------------------------------------
enum Variables {    // create an enumeration item for every variable from your CF container
  kNothing = -1,    // kNothing should be always here
  kPt = 0,
  kY,
  kPairType,
  kThetaCS,
  kThetaHE,
  kM,
  kLeg1_Eta,
  //  kLeg1_TPC_nSigma_Electrons,
  kLeg1_Pt,
  //  kLeg1_P,
  kLeg1_NclsTPC,
  kLeg2_Eta,
  //  kLeg2_TPC_nSigma_Electrons,
  kLeg2_Pt,
  //  kLeg2_P,
  kLeg2_NclsTPC,
  kNVariables       // kNVariables should be always here!
};
const Char_t* gkVarNames[kNVariables] = {     // variable names to be put on histograms axes and titles
  "p_{T} [GeV/c]",
  "y",
  "Pair type (0=++; 1=+-; 2=--)",
  "cos #theta^{*}_{CS}",
  "cos #theta^{*}_{HE}",
  "M [GeV/c^{2}]",
  "#eta^{leg1}",
  //  "TPC n #sigma electrons (leg1)",
  "p_{T}^{leg1} [GeV/c]",
  //  "P^{leg1} [GeV/c]",
  "# TPC clusters (leg1)",
  "#eta^{leg2}",
  //  "TPC n #sigma electrons (leg2)",
  "p_{T}^{leg2} [GeV/c]",
  //  "P^{leg2} [GeV/c]",
  "# TPC clusters (leg2)"
};
Int_t gNbins[kNVariables];           // number of bins for every variable --> filled automatically
Double_t* gBinLimits[kNVariables];   // bin limits for every variable --> filled automatically
// ------------------------------------------------------------------------

// Put here all the CF steps of interest ----------------------------------
enum Steps {        // step indexes in the CF containers to be analyzed
  kPureMC = 0,
  kESDSPDany = 2,
  kESDSPDfirst = 4,
  kFullSPDany = 6,
  kFullSPDfirst = 8,
  kNSteps = 5        // total number of steps (the number of steps above)
};
const Int_t gkStepNumbers[kNSteps] = {   // array with step indexes (all from the enumeration above)
  kPureMC, 
  kESDSPDany, kESDSPDfirst,
  kFullSPDany, kFullSPDfirst
};
const Char_t* gkStepNames[kNSteps][2] = {// names for each CF step
  {"PureMC",         "Pure MC"},         // NOTE: short names go to histo names, long names go to titles
  {"ESDSPDany",      "ESD track cuts, SPD any, TPCnclus>90"}, 
  {"ESDSPDfirst",    "ESD track cuts, SPD first, TPCnclus>90"},
  {"FullSPDany",     "All track cuts (with SPD any) and TPC-PID"},
  {"FullSPDfirst",   "All track cuts (with SPD first) and TPC-PID"}
};
//------------------------------------------------------------------------

// Put here info about the cut sets for which projections will be made ---
const Int_t gkNCutSets = 18;     // number of cut sets for which histos will be filled
const Char_t* gkCutSetNames[gkNCutSets][2] = {   // short and long names for all the cut sets
  // baseline
  {"Ycut",                 "|y_{J/#Psi}|<0.88"},
  {"YcutPt1",              "|y_{J/#Psi}|<0.88 & 0.0<p_{T J/#Psi}<0.8"},
  {"YcutPt2",              "|y_{J/#Psi}|<0.88 & 0.8<p_{T J/#Psi}<1.4"},
  {"YcutPt3",              "|y_{J/#Psi}|<0.88 & 1.4<p_{T J/#Psi}<2.8"},
  {"YcutPt4",              "|y_{J/#Psi}|<0.88 & 2.8<p_{T J/#Psi}<5.0"},
  {"YcutPt5",              "|y_{J/#Psi}|<0.88 & 5.0<p_{T J/#Psi}<10.0"},
  // pure kinematic acceptance
  {"YcutLegsEtaPt1",       "|y_{J/#Psi}|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c"},
  {"YcutPt1LegsEtaPt1",    "|y_{J/#Psi}|<0.88 & 0.0<p_{T J/#Psi}<0.8 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c"},
  {"YcutPt2LegsEtaPt1",    "|y_{J/#Psi}|<0.88 & 0.8<p_{T J/#Psi}<1.4 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c"},
  {"YcutPt3LegsEtaPt1",    "|y_{J/#Psi}|<0.88 & 1.4<p_{T J/#Psi}<2.8 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c"},
  {"YcutPt4LegsEtaPt1",    "|y_{J/#Psi}|<0.88 & 2.8<p_{T J/#Psi}<5.0 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c"},
  {"YcutPt5LegsEtaPt1",    "|y_{J/#Psi}|<0.88 & 5.0<p_{T J/#Psi}<10.0 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c"},
  // track quality
  {"YcutLegsEtaPt1TPC90",     "|y_{J/#Psi}|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90"},
  {"YcutPt1LegsEtaPt1TPC90",  "|y_{J/#Psi}|<0.88 & 0.0<p_{T J/#Psi}<0.8 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90"},
  {"YcutPt2LegsEtaPt1TPC90",  "|y_{J/#Psi}|<0.88 & 0.8<p_{T J/#Psi}<1.4 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90"},
  {"YcutPt3LegsEtaPt1TPC90",  "|y_{J/#Psi}|<0.88 & 1.4<p_{T J/#Psi}<2.8 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90"},
  {"YcutPt4LegsEtaPt1TPC90",  "|y_{J/#Psi}|<0.88 & 2.8<p_{T J/#Psi}<5.0 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90"},
  {"YcutPt5LegsEtaPt1TPC90",  "|y_{J/#Psi}|<0.88 & 5.0<p_{T J/#Psi}<10.0 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90"}
  // track quality + PID
  //  {"YcutLegsEtaPt1TPC90PID",     "|y_{J/#Psi}|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90  & TPC pid"},
  //  {"YcutPt1LegsEtaPt1TPC90PID",  "|y_{J/#Psi}|<0.88 & 0.0<p_{T J/#Psi}<0.8 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90  & TPC pid"},
  //  {"YcutPt2LegsEtaPt1TPC90PID",  "|y_{J/#Psi}|<0.88 & 0.8<p_{T J/#Psi}<1.4 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90  & TPC pid"},
  //  {"YcutPt3LegsEtaPt1TPC90PID",  "|y_{J/#Psi}|<0.88 & 1.4<p_{T J/#Psi}<2.8 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90  & TPC pid"},
  //  {"YcutPt4LegsEtaPt1TPC90PID",  "|y_{J/#Psi}|<0.88 & 2.8<p_{T J/#Psi}<5.0 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90  & TPC pid"},
  //  {"YcutPt5LegsEtaPt1TPC90PID",  "|y_{J/#Psi}|<0.88 & 5.0<p_{T J/#Psi}<10.0 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c TPCnclus>90  & TPC pid"}
};
// -----------------------------------------------------------------------

// Put here info about the histograms to be filled -----------------------
const Int_t gkNhistos = 6;                        // how many histograms for every (step,cut set) combination
const Char_t* gkHistoNames[gkNhistos][2] = {      // short and long names of the histograms
  {"pt","p_{T}(J/#Psi)"},                         // NOTE: short names go to the histo name, long name goes to title
  {"y","y(J/#Psi)"}, 
  {"pty","p_{T} vs y(J/#Psi)"}, 
  //  {"m","e^{+}e^{-} invariant mass"}, 
  {"ThetaCS","cos #theta^{*}_{CS}"}, 
  {"ThetaHE","cos #theta^{*}_{HE}"},
  {"Minv", "Invariant mass"}
};
const Int_t gkDims[gkNhistos][4] = {      // dimensions and variables for histograms
// ndim  xVar      yVar      zVar
  {1,    kPt,      kNothing, kNothing},   // pt dependence
  {1,    kY,       kNothing, kNothing},   // y dependence
  {2,    kY,       kPt,      kNothing},   // pt,y dependence
  //  {1,    kM,       kNothing, kNothing},   // inv. mass dependence
  {1,    kThetaCS, kNothing, kNothing},   // cos theta* CS dependence
  {1,    kThetaHE, kNothing, kNothing},    // cos theta* HE dependence
  {1,    kM,       kNothing, kNothing}    // invariant mass
};
// -----------------------------------------------------------------------

// ******************************************************************************** 
// Define here all the efficiencies you want (if any)
// Efficiency = (nominator histogram)/(denominator histogram)
// A histogram is defined by its step,cut set, and number defined according the
// global variables above
// ********************************************************************************
const Int_t gkNeffs = 210;
const Int_t gkEffs[gkNeffs][6] = {
  //nominator: Step  Cut  Histo | denominator: Step  Cut  Histo     comment
  // full corrections, pt dependence
             { 3,    12,   0,                   0,    0,   0     },     // full correction, SPD any
             { 4,    12,   0,                   0,    0,   0     },     // full correction, SPD first
             { 3,    13,   0,                   0,    1,   0     },     // full correction, SPD any,  0.0<pt<0.8
             { 4,    13,   0,                   0,    1,   0     },     // full correction, SPD first,0.0<pt<0.8 
             { 3,    14,   0,                   0,    2,   0     },     // full correction, SPD any,  0.8<pt<1.4
             { 4,    14,   0,                   0,    2,   0     },     // full correction, SPD first,0.8<pt<1.4
             { 3,    15,   0,                   0,    3,   0     },     // full correction, SPD any,  1.4<pt<2.8
             { 4,    15,   0,                   0,    3,   0     },     // full correction, SPD first,1.4<pt<2.8
             { 3,    16,   0,                   0,    4,   0     },     // full correction, SPD any,  2.8<pt<5.0
             { 4,    16,   0,                   0,    4,   0     },     // full correction, SPD first,2.8<pt<5.0
             { 3,    17,   0,                   0,    5,   0     },     // full correction, SPD any,  5.0<pt<10.0
             { 4,    17,   0,                   0,    5,   0     },     // full correction, SPD first,5.0<pt<10.0
	     // full corrections, y dependence
	     { 3,    12,   1,                   0,    0,   1     },     // full correction, SPD any
             { 4,    12,   1,                   0,    0,   1     },     // full correction, SPD first
             { 3,    13,   1,                   0,    1,   1     },     // full correction, SPD any,  0.0<pt<0.8
             { 4,    13,   1,                   0,    1,   1     },     // full correction, SPD first,0.0<pt<0.8 
             { 3,    14,   1,                   0,    2,   1     },     // full correction, SPD any,  0.8<pt<1.4
             { 4,    14,   1,                   0,    2,   1     },     // full correction, SPD first,0.8<pt<1.4
             { 3,    15,   1,                   0,    3,   1     },     // full correction, SPD any,  1.4<pt<2.8
             { 4,    15,   1,                   0,    3,   1     },     // full correction, SPD first,1.4<pt<2.8
             { 3,    16,   1,                   0,    4,   1     },     // full correction, SPD any,  2.8<pt<5.0
             { 4,    16,   1,                   0,    4,   1     },     // full correction, SPD first,2.8<pt<5.0
             { 3,    17,   1,                   0,    5,   1     },     // full correction, SPD any,  5.0<pt<10.0
             { 4,    17,   1,                   0,    5,   1     },     // full correction, SPD first,5.0<pt<10.0
	     // full corrections, pt,y dependence
	     { 3,    12,   2,                   0,    0,   2     },     // full correction, SPD any
             { 4,    12,   2,                   0,    0,   2     },     // full correction, SPD first
             { 3,    13,   2,                   0,    1,   2     },     // full correction, SPD any,  0.0<pt<0.8
             { 4,    13,   2,                   0,    1,   2     },     // full correction, SPD first,0.0<pt<0.8 
             { 3,    14,   2,                   0,    2,   2     },     // full correction, SPD any,  0.8<pt<1.4
             { 4,    14,   2,                   0,    2,   2     },     // full correction, SPD first,0.8<pt<1.4
             { 3,    15,   2,                   0,    3,   2     },     // full correction, SPD any,  1.4<pt<2.8
             { 4,    15,   2,                   0,    3,   2     },     // full correction, SPD first,1.4<pt<2.8
             { 3,    16,   2,                   0,    4,   2     },     // full correction, SPD any,  2.8<pt<5.0
             { 4,    16,   2,                   0,    4,   2     },     // full correction, SPD first,2.8<pt<5.0
             { 3,    17,   2,                   0,    5,   2     },     // full correction, SPD any,  5.0<pt<10.0
             { 4,    17,   2,                   0,    5,   2     },     // full correction, SPD first,5.0<pt<10.0
	     // full corrections, cos* Theta CS dependence
	     { 3,    12,   3,                   0,    0,   3     },     // full correction, SPD any
             { 4,    12,   3,                   0,    0,   3     },     // full correction, SPD first
             { 3,    13,   3,                   0,    1,   3     },     // full correction, SPD any,  0.0<pt<0.8
             { 4,    13,   3,                   0,    1,   3     },     // full correction, SPD first,0.0<pt<0.8 
             { 3,    14,   3,                   0,    2,   3     },     // full correction, SPD any,  0.8<pt<1.4
             { 4,    14,   3,                   0,    2,   3     },     // full correction, SPD first,0.8<pt<1.4
             { 3,    15,   3,                   0,    3,   3     },     // full correction, SPD any,  1.4<pt<2.8
             { 4,    15,   3,                   0,    3,   3     },     // full correction, SPD first,1.4<pt<2.8
             { 3,    16,   3,                   0,    4,   3     },     // full correction, SPD any,  2.8<pt<5.0
             { 4,    16,   3,                   0,    4,   3     },     // full correction, SPD first,2.8<pt<5.0
             { 3,    17,   3,                   0,    5,   3     },     // full correction, SPD any,  5.0<pt<10.0
             { 4,    17,   3,                   0,    5,   3     },     // full correction, SPD first,5.0<pt<10.0
	     // full corrections, cos* Theta HE dependence
	     { 3,    12,   4,                   0,    0,   4     },     // full correction, SPD any
             { 4,    12,   4,                   0,    0,   4     },     // full correction, SPD first
             { 3,    13,   4,                   0,    1,   4     },     // full correction, SPD any,  0.0<pt<0.8
             { 4,    13,   4,                   0,    1,   4     },     // full correction, SPD first,0.0<pt<0.8 
             { 3,    14,   4,                   0,    2,   4     },     // full correction, SPD any,  0.8<pt<1.4
             { 4,    14,   4,                   0,    2,   4     },     // full correction, SPD first,0.8<pt<1.4
             { 3,    15,   4,                   0,    3,   4     },     // full correction, SPD any,  1.4<pt<2.8
             { 4,    15,   4,                   0,    3,   4     },     // full correction, SPD first,1.4<pt<2.8
             { 3,    16,   4,                   0,    4,   4     },     // full correction, SPD any,  2.8<pt<5.0
             { 4,    16,   4,                   0,    4,   4     },     // full correction, SPD first,2.8<pt<5.0
             { 3,    17,   4,                   0,    5,   4     },     // full correction, SPD any,  5.0<pt<10.0
             { 4,    17,   4,                   0,    5,   4     },     // full correction, SPD first,5.0<pt<10.0

	     // acceptance corrections, pt dependence
	     { 0,     6,   0,                   0,    0,   0     },     // acc. correction, full pt range
	     { 0,     7,   0,                   0,    1,   0     },     // acc. correction, 0.0<pt<0.8
	     { 0,     8,   0,                   0,    2,   0     },     // acc. correction, 0.8<pt<1.4
	     { 0,     9,   0,                   0,    3,   0     },     // acc. correction, 1.4<pt<2.8
	     { 0,    10,   0,                   0,    4,   0     },     // acc. correction, 2.8<pt<5.0
	     { 0,    11,   0,                   0,    5,   0     },     // acc. correction, 5.0<pt<10.0
	     // acceptance corrections, y dependence
	     { 0,     6,   1,                   0,    0,   1     },     // acc. correction, full pt range
	     { 0,     7,   1,                   0,    1,   1     },     // acc. correction, 0.0<pt<0.8
	     { 0,     8,   1,                   0,    2,   1     },     // acc. correction, 0.8<pt<1.4
	     { 0,     9,   1,                   0,    3,   1     },     // acc. correction, 1.4<pt<2.8
	     { 0,    10,   1,                   0,    4,   1     },     // acc. correction, 2.8<pt<5.0
	     { 0,    11,   1,                   0,    5,   1     },     // acc. correction, 5.0<pt<10.0
	     // acceptance corrections, pt,y dependence
	     { 0,     6,   2,                   0,    0,   2     },     // acc. correction, full pt range
	     { 0,     7,   2,                   0,    1,   2     },     // acc. correction, 0.0<pt<0.8
	     { 0,     8,   2,                   0,    2,   2     },     // acc. correction, 0.8<pt<1.4
	     { 0,     9,   2,                   0,    3,   2     },     // acc. correction, 1.4<pt<2.8
	     { 0,    10,   2,                   0,    4,   2     },     // acc. correction, 2.8<pt<5.0
	     { 0,    11,   2,                   0,    5,   2     },     // acc. correction, 5.0<pt<10.0
	     // acceptance corrections, cos Theta* CS dependence
	     { 0,     6,   3,                   0,    0,   3     },     // acc. correction, full pt range
	     { 0,     7,   3,                   0,    1,   3     },     // acc. correction, 0.0<pt<0.8
	     { 0,     8,   3,                   0,    2,   3     },     // acc. correction, 0.8<pt<1.4
	     { 0,     9,   3,                   0,    3,   3     },     // acc. correction, 1.4<pt<2.8
	     { 0,    10,   3,                   0,    4,   3     },     // acc. correction, 2.8<pt<5.0
	     { 0,    11,   3,                   0,    5,   3     },     // acc. correction, 5.0<pt<10.0
	     // acceptance corrections, cos Theta* HE dependence
	     { 0,     6,   4,                   0,    0,   4     },     // acc. correction, full pt range
	     { 0,     7,   4,                   0,    1,   4     },     // acc. correction, 0.0<pt<0.8
	     { 0,     8,   4,                   0,    2,   4     },     // acc. correction, 0.8<pt<1.4
	     { 0,     9,   4,                   0,    3,   4     },     // acc. correction, 1.4<pt<2.8
	     { 0,    10,   4,                   0,    4,   4     },     // acc. correction, 2.8<pt<5.0
	     { 0,    11,   4,                   0,    5,   4     },     // acc. correction, 5.0<pt<10.0


	     // tracking corrections with SPD any and TPCncls>90, pt dependence
	     { 1,    12,   0,                   0,    6,   0     },     // tracking correction, full pt range
	     { 1,    13,   0,                   0,    7,   0     },     // tracking correction, 0.0<pt<0.8
	     { 1,    14,   0,                   0,    8,   0     },     // tracking correction, 0.8<pt<1.4
	     { 1,    15,   0,                   0,    9,   0     },     // tracking correction, 1.4<pt<2.8
	     { 1,    16,   0,                   0,   10,   0     },     // tracking correction, 2.8<pt<5.0
	     { 1,    17,   0,                   0,   11,   0     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD any and TPCncls>90, y dependence
	     { 1,    12,   1,                   0,    6,   1     },     // tracking correction, full pt range
	     { 1,    13,   1,                   0,    7,   1     },     // tracking correction, 0.0<pt<0.8
	     { 1,    14,   1,                   0,    8,   1     },     // tracking correction, 0.8<pt<1.4
	     { 1,    15,   1,                   0,    9,   1     },     // tracking correction, 1.4<pt<2.8
	     { 1,    16,   1,                   0,   10,   1     },     // tracking correction, 2.8<pt<5.0
	     { 1,    17,   1,                   0,   11,   1     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD any and TPCncls>90, pt,y dependence
	     { 1,    12,   2,                   0,    6,   2     },     // tracking correction, full pt range
	     { 1,    13,   2,                   0,    7,   2     },     // tracking correction, 0.0<pt<0.8
	     { 1,    14,   2,                   0,    8,   2     },     // tracking correction, 0.8<pt<1.4
	     { 1,    15,   2,                   0,    9,   2     },     // tracking correction, 1.4<pt<2.8
	     { 1,    16,   2,                   0,   10,   2     },     // tracking correction, 2.8<pt<5.0
	     { 1,    17,   2,                   0,   11,   2     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD any and TPCncls>90, cos Theta* CS dependence
	     { 1,    12,   3,                   0,    6,   3     },     // tracking correction, full pt range
	     { 1,    13,   3,                   0,    7,   3     },     // tracking correction, 0.0<pt<0.8
	     { 1,    14,   3,                   0,    8,   3     },     // tracking correction, 0.8<pt<1.4
	     { 1,    15,   3,                   0,    9,   3     },     // tracking correction, 1.4<pt<2.8
	     { 1,    16,   3,                   0,   10,   3     },     // tracking correction, 2.8<pt<5.0
	     { 1,    17,   3,                   0,   11,   3     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD any and TPCncls>90, cos Theta* HE dependence
	     { 1,    12,   4,                   0,    6,   4     },     // tracking correction, full pt range
	     { 1,    13,   4,                   0,    7,   4     },     // tracking correction, 0.0<pt<0.8
	     { 1,    14,   4,                   0,    8,   4     },     // tracking correction, 0.8<pt<1.4
	     { 1,    15,   4,                   0,    9,   4     },     // tracking correction, 1.4<pt<2.8
	     { 1,    16,   4,                   0,   10,   4     },     // tracking correction, 2.8<pt<5.0
	     { 1,    17,   4,                   0,   11,   4     },     // tracking correction, 5.0<pt<10.0

	     // tracking corrections with SPD first and TPCncls>90, pt dependence
	     { 2,    12,   0,                   0,    6,   0     },     // tracking correction, full pt range
	     { 2,    13,   0,                   0,    7,   0     },     // tracking correction, 0.0<pt<0.8
	     { 2,    14,   0,                   0,    8,   0     },     // tracking correction, 0.8<pt<1.4
	     { 2,    15,   0,                   0,    9,   0     },     // tracking correction, 1.4<pt<2.8
	     { 2,    16,   0,                   0,   10,   0     },     // tracking correction, 2.8<pt<5.0
	     { 2,    17,   0,                   0,   11,   0     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD first and TPCncls>90, y dependence
	     { 2,    12,   1,                   0,    6,   1     },     // tracking correction, full pt range
	     { 2,    13,   1,                   0,    7,   1     },     // tracking correction, 0.0<pt<0.8
	     { 2,    14,   1,                   0,    8,   1     },     // tracking correction, 0.8<pt<1.4
	     { 2,    15,   1,                   0,    9,   1     },     // tracking correction, 1.4<pt<2.8
	     { 2,    16,   1,                   0,   10,   1     },     // tracking correction, 2.8<pt<5.0
	     { 2,    17,   1,                   0,   11,   1     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD first and TPCncls>90, pt,y dependence
	     { 2,    12,   2,                   0,    6,   2     },     // tracking correction, full pt range
	     { 2,    13,   2,                   0,    7,   2     },     // tracking correction, 0.0<pt<0.8
	     { 2,    14,   2,                   0,    8,   2     },     // tracking correction, 0.8<pt<1.4
	     { 2,    15,   2,                   0,    9,   2     },     // tracking correction, 1.4<pt<2.8
	     { 2,    16,   2,                   0,   10,   2     },     // tracking correction, 2.8<pt<5.0
	     { 2,    17,   2,                   0,   11,   2     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD first and TPCncls>90, cos Theta*CS dependence
	     { 2,    12,   3,                   0,    6,   3     },     // tracking correction, full pt range
	     { 2,    13,   3,                   0,    7,   3     },     // tracking correction, 0.0<pt<0.8
	     { 2,    14,   3,                   0,    8,   3     },     // tracking correction, 0.8<pt<1.4
	     { 2,    15,   3,                   0,    9,   3     },     // tracking correction, 1.4<pt<2.8
	     { 2,    16,   3,                   0,   10,   3     },     // tracking correction, 2.8<pt<5.0
	     { 2,    17,   3,                   0,   11,   3     },     // tracking correction, 5.0<pt<10.0
	     // tracking corrections with SPD first and TPCncls>90, cos Theta*HE dependence
	     { 2,    12,   4,                   0,    6,   4     },     // tracking correction, full pt range
	     { 2,    13,   4,                   0,    7,   4     },     // tracking correction, 0.0<pt<0.8
	     { 2,    14,   4,                   0,    8,   4     },     // tracking correction, 0.8<pt<1.4
	     { 2,    15,   4,                   0,    9,   4     },     // tracking correction, 1.4<pt<2.8
	     { 2,    16,   4,                   0,   10,   4     },     // tracking correction, 2.8<pt<5.0
	     { 2,    17,   4,                   0,   11,   4     },     // tracking correction, 5.0<pt<10.0

	     // PID corrections with SPD any, pt dependence
	     { 3,    12,   0,                   1,   12,   0     },     // PID correction, full pt range
	     { 3,    13,   0,                   1,   13,   0     },     // PID correction, 0.0<pt<0.8
	     { 3,    14,   0,                   1,   14,   0     },     // PID correction, 0.8<pt<1.4
	     { 3,    15,   0,                   1,   15,   0     },     // PID correction, 1.4<pt<2.8
	     { 3,    16,   0,                   1,   16,   0     },     // PID correction, 2.8<pt<5.0
	     { 3,    17,   0,                   1,   17,   0     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD any, y dependence
	     { 3,    12,   1,                   1,   12,   1     },     // PID correction, full pt range
	     { 3,    13,   1,                   1,   13,   1     },     // PID correction, 0.0<pt<0.8
	     { 3,    14,   1,                   1,   14,   1     },     // PID correction, 0.8<pt<1.4
	     { 3,    15,   1,                   1,   15,   1     },     // PID correction, 1.4<pt<2.8
	     { 3,    16,   1,                   1,   16,   1     },     // PID correction, 2.8<pt<5.0
	     { 3,    17,   1,                   1,   17,   1     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD any, pt,y dependence
	     { 3,    12,   2,                   1,   12,   2     },     // PID correction, full pt range
	     { 3,    13,   2,                   1,   13,   2     },     // PID correction, 0.0<pt<0.8
	     { 3,    14,   2,                   1,   14,   2     },     // PID correction, 0.8<pt<1.4
	     { 3,    15,   2,                   1,   15,   2     },     // PID correction, 1.4<pt<2.8
	     { 3,    16,   2,                   1,   16,   2     },     // PID correction, 2.8<pt<5.0
	     { 3,    17,   2,                   1,   17,   2     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD any, cos Theta*CS dependence
	     { 3,    12,   3,                   1,   12,   3     },     // PID correction, full pt range
	     { 3,    13,   3,                   1,   13,   3     },     // PID correction, 0.0<pt<0.8
	     { 3,    14,   3,                   1,   14,   3     },     // PID correction, 0.8<pt<1.4
	     { 3,    15,   3,                   1,   15,   3     },     // PID correction, 1.4<pt<2.8
	     { 3,    16,   3,                   1,   16,   3     },     // PID correction, 2.8<pt<5.0
	     { 3,    17,   3,                   1,   17,   3     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD any, cos Theta*HE dependence
	     { 3,    12,   4,                   1,   12,   4     },     // PID correction, full pt range
	     { 3,    13,   4,                   1,   13,   4     },     // PID correction, 0.0<pt<0.8
	     { 3,    14,   4,                   1,   14,   4     },     // PID correction, 0.8<pt<1.4
	     { 3,    15,   4,                   1,   15,   4     },     // PID correction, 1.4<pt<2.8
	     { 3,    16,   4,                   1,   16,   4     },     // PID correction, 2.8<pt<5.0
	     { 3,    17,   4,                   1,   17,   4     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD first, pt dependence
	     { 4,    12,   0,                   2,   12,   0     },     // PID correction, full pt range
	     { 4,    13,   0,                   2,   13,   0     },     // PID correction, 0.0<pt<0.8
	     { 4,    14,   0,                   2,   14,   0     },     // PID correction, 0.8<pt<1.4
	     { 4,    15,   0,                   2,   15,   0     },     // PID correction, 1.4<pt<2.8
	     { 4,    16,   0,                   2,   16,   0     },     // PID correction, 2.8<pt<5.0
	     { 4,    17,   0,                   2,   17,   0     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD first, y dependence
	     { 4,    12,   1,                   2,   12,   1     },     // PID correction, full pt range
	     { 4,    13,   1,                   2,   13,   1     },     // PID correction, 0.0<pt<0.8
	     { 4,    14,   1,                   2,   14,   1     },     // PID correction, 0.8<pt<1.4
	     { 4,    15,   1,                   2,   15,   1     },     // PID correction, 1.4<pt<2.8
	     { 4,    16,   1,                   2,   16,   1     },     // PID correction, 2.8<pt<5.0
	     { 4,    17,   1,                   2,   17,   1     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD first, pt,y dependence
	     { 4,    12,   2,                   2,   12,   2     },     // PID correction, full pt range
	     { 4,    13,   2,                   2,   13,   2     },     // PID correction, 0.0<pt<0.8
	     { 4,    14,   2,                   2,   14,   2     },     // PID correction, 0.8<pt<1.4
	     { 4,    15,   2,                   2,   15,   2     },     // PID correction, 1.4<pt<2.8
	     { 4,    16,   2,                   2,   16,   2     },     // PID correction, 2.8<pt<5.0
	     { 4,    17,   2,                   2,   17,   2     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD first, cos Theta*CS dependence
	     { 4,    12,   3,                   2,   12,   3     },     // PID correction, full pt range
	     { 4,    13,   3,                   2,   13,   3     },     // PID correction, 0.0<pt<0.8
	     { 4,    14,   3,                   2,   14,   3     },     // PID correction, 0.8<pt<1.4
	     { 4,    15,   3,                   2,   15,   3     },     // PID correction, 1.4<pt<2.8
	     { 4,    16,   3,                   2,   16,   3     },     // PID correction, 2.8<pt<5.0
	     { 4,    17,   3,                   2,   17,   3     },     // PID correction, 5.0<pt<10.0
	     // PID corrections with SPD first, cos Theta*HE dependence
	     { 4,    12,   4,                   2,   12,   4     },     // PID correction, full pt range
	     { 4,    13,   4,                   2,   13,   4     },     // PID correction, 0.0<pt<0.8
	     { 4,    14,   4,                   2,   14,   4     },     // PID correction, 0.8<pt<1.4
	     { 4,    15,   4,                   2,   15,   4     },     // PID correction, 1.4<pt<2.8
	     { 4,    16,   4,                   2,   16,   4     },     // PID correction, 2.8<pt<5.0
	     { 4,    17,   4,                   2,   17,   4     }      // PID correction, 5.0<pt<10.0

};
// custom names and titles for efficiency histograms
const Char_t* gkEffNames[gkNeffs][2] = {
  // full corrections, pt dependence
  {"fullCorrectionSPDany_pt",      "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionSPDfirst_pt",    "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionPt1SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt1SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt2SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt2SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt3SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt3SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt4SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt4SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt5SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  {"fullCorrectionPt5SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  // full corrections, y dependence
  {"fullCorrectionSPDany_y",      "Eff. vs. y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionSPDfirst_y",    "Eff. vs. y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionPt1SPDany_y",   "Eff. vs. y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt1SPDfirst_y", "Eff. vs. y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt2SPDany_y",   "Eff. vs. y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt2SPDfirst_y", "Eff. vs. y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt3SPDany_y",   "Eff. vs. y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt3SPDfirst_y", "Eff. vs. y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt4SPDany_y",   "Eff. vs. y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt4SPDfirst_y", "Eff. vs. y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt5SPDany_y",   "Eff. vs. y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  {"fullCorrectionPt5SPDfirst_y", "Eff. vs. y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  // full corrections, pt-y dependence
  {"fullCorrectionSPDany_pty",      "Eff. vs. pt-y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionSPDfirst_pty",    "Eff. vs. pt-y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionPt1SPDany_pty",   "Eff. vs. pt-y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt1SPDfirst_pty", "Eff. vs. pt-y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt2SPDany_pty",   "Eff. vs. pt-y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt2SPDfirst_pty", "Eff. vs. pt-y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt3SPDany_pty",   "Eff. vs. pt-y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt3SPDfirst_pty", "Eff. vs. pt-y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt4SPDany_pty",   "Eff. vs. pt-y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt4SPDfirst_pty", "Eff. vs. pt-y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt5SPDany_pty",   "Eff. vs. pt-y, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  {"fullCorrectionPt5SPDfirst_pty", "Eff. vs. pt-y, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  // full corrections, cos Theta*CS dependence
  {"fullCorrectionSPDany_ThetaCS",      "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionSPDfirst_ThetaCS",    "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionPt1SPDany_ThetaCS",   "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt1SPDfirst_ThetaCS", "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt2SPDany_ThetaCS",   "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt2SPDfirst_ThetaCS", "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt3SPDany_ThetaCS",   "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt3SPDfirst_ThetaCS", "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt4SPDany_ThetaCS",   "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt4SPDfirst_ThetaCS", "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt5SPDany_ThetaCS",   "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  {"fullCorrectionPt5SPDfirst_ThetaCS", "Eff. vs. cos #theta^{*}_{CS}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  // full corrections, cos Theta*HE dependence
  {"fullCorrectionSPDany_ThetaHE",      "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionSPDfirst_ThetaHE",    "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88)"},
  {"fullCorrectionPt1SPDany_ThetaHE",   "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt1SPDfirst_ThetaHE", "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0<pt<0.8)"},
  {"fullCorrectionPt2SPDany_ThetaHE",   "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt2SPDfirst_ThetaHE", "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 0.8<pt<1.4)"},
  {"fullCorrectionPt3SPDany_ThetaHE",   "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt3SPDfirst_ThetaHE", "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 1.4<pt<2.8)"},
  {"fullCorrectionPt4SPDany_ThetaHE",   "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt4SPDfirst_ThetaHE", "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 2.8<pt<5.0)"},
  {"fullCorrectionPt5SPDany_ThetaHE",   "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  {"fullCorrectionPt5SPDfirst_ThetaHE", "Eff. vs. cos #theta^{*}_{HE}, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.88 & 5.0<pt<10.0)"},
  // acceptance corrections, pt dependence
  {"accCorrection_pt", "Kinematic correction vs pt, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt1_pt", "Kinematic correction vs pt, 0.0<p_{T}(J/#Psi)<0.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt2_pt", "Kinematic correction vs pt, 0.8<p_{T}(J/#Psi)<1.4, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt3_pt", "Kinematic correction vs pt, 1.4<p_{T}(J/#Psi)<2.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt4_pt", "Kinematic correction vs pt, 2.8<p_{T}(J/#Psi)<5.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt5_pt", "Kinematic correction vs pt, 5.0<p_{T}(J/#Psi)<10.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  // acceptance corrections, y dependence
  {"accCorrection_y", "Kinematic correction vs y, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt1_y", "Kinematic correction vs y, 0.0<p_{T}(J/#Psi)<0.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt2_y", "Kinematic correction vs y, 0.8<p_{T}(J/#Psi)<1.4, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt3_y", "Kinematic correction vs y, 1.4<p_{T}(J/#Psi)<2.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt4_y", "Kinematic correction vs y, 2.8<p_{T}(J/#Psi)<5.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt5_y", "Kinematic correction vs y, 5.0<p_{T}(J/#Psi)<10.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  // acceptance corrections, pt-y dependence
  {"accCorrection_pty", "Kinematic correction vs pt-y, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt1_pty", "Kinematic correction vs pt-y, 0.0<p_{T}(J/#Psi)<0.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt2_pty", "Kinematic correction vs pt-y, 0.8<p_{T}(J/#Psi)<1.4, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt3_pty", "Kinematic correction vs pt-y, 1.4<p_{T}(J/#Psi)<2.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt4_pty", "Kinematic correction vs pt-y, 2.8<p_{T}(J/#Psi)<5.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt5_pty", "Kinematic correction vs pt-y, 5.0<p_{T}(J/#Psi)<10.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  // acceptance corrections, cos Theta*CS dependence
  {"accCorrection_ThetaCS", "Kinematic correction vs cos #theta^{*}_{CS}, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt1_ThetaCS", "Kinematic correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<0.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt2_ThetaCS", "Kinematic correction vs cos #theta^{*}_{CS}, 0.8<p_{T}(J/#Psi)<1.4, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt3_ThetaCS", "Kinematic correction vs cos #theta^{*}_{CS}, 1.4<p_{T}(J/#Psi)<2.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt4_ThetaCS", "Kinematic correction vs cos #theta^{*}_{CS}, 2.8<p_{T}(J/#Psi)<5.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt5_ThetaCS", "Kinematic correction vs cos #theta^{*}_{CS}, 5.0<p_{T}(J/#Psi)<10.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  // acceptance corrections, cos Theta*HE dependence
  {"accCorrection_ThetaHE", "Kinematic correction vs cos #theta^{*}_{HE}, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt1_ThetaHE", "Kinematic correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<0.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt2_ThetaHE", "Kinematic correction vs cos #theta^{*}_{HE}, 0.8<p_{T}(J/#Psi)<1.4, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt3_ThetaHE", "Kinematic correction vs cos #theta^{*}_{HE}, 1.4<p_{T}(J/#Psi)<2.8, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt4_ThetaHE", "Kinematic correction vs cos #theta^{*}_{HE}, 2.8<p_{T}(J/#Psi)<5.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  {"accCorrectionPt5_ThetaHE", "Kinematic correction vs cos #theta^{*}_{HE}, 5.0<p_{T}(J/#Psi)<10.0, (J/#Psi in |y|<0.88 & |#eta_{legs}|<0.88 & p_{Tlegs}>1.0 GeV/c)/(J/#Psi in |y|<0.88)"},
  // tracking corrections with SPD any and TPCncls>90, pt dependence
  {"trackingCorrectionSPDany_pt", "Tracking correction vs pt, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt1_pt", "Tracking correction vs pt, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt2_pt", "Tracking correction vs pt, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt3_pt", "Tracking correction vs pt, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt4_pt", "Tracking correction vs pt, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt5_pt", "Tracking correction vs pt, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // tracking corrections with SPD any and TPCncls>90, y dependence
  {"trackingCorrectionSPDany_y", "Tracking correction vs y, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt1_y", "Tracking correction vs y, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt2_y", "Tracking correction vs y, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt3_y", "Tracking correction vs y, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt4_y", "Tracking correction vs y, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt5_y", "Tracking correction vs y, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // tracking corrections with SPD any and TPCncls>90, pt-y dependence
  {"trackingCorrectionSPDany_pty", "Tracking correction vs pt-y, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt1_pty", "Tracking correction vs pt-y, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt2_pty", "Tracking correction vs pt-y, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt3_pty", "Tracking correction vs pt-y, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt4_pty", "Tracking correction vs pt-y, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt5_pty", "Tracking correction vs pt-y, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // tracking corrections with SPD any and TPCncls>90, cos Theta*CS dependence
  {"trackingCorrectionSPDany_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt1_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt2_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt3_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt4_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt5_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // tracking corrections with SPD any and TPCncls>90, cos Theta*HE dependence
  {"trackingCorrectionSPDany_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt1_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt2_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt3_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt4_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"trackingCorrectionSPDanyPt5_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},

  // tracking corrections with SPD first and TPCncls>90, pt dependence
  {"trackingCorrectionSPDfirst_pt", "Tracking correction vs pt, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt1_pt", "Tracking correction vs pt, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt2_pt", "Tracking correction vs pt, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt3_pt", "Tracking correction vs pt, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt4_pt", "Tracking correction vs pt, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt5_pt", "Tracking correction vs pt, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // tracking corrections with SPD first and TPCncls>90, y dependence
  {"trackingCorrectionSPDfirst_y", "Tracking correction vs y, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt1_y", "Tracking correction vs y, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt2_y", "Tracking correction vs y, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt3_y", "Tracking correction vs y, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt4_y", "Tracking correction vs y, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt5_y", "Tracking correction vs y, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // tracking corrections with SPD first and TPCncls>90, pt-y dependence
  {"trackingCorrectionSPDfirst_pty", "Tracking correction vs pt-y, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt1_pty", "Tracking correction vs pt-y, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt2_pty", "Tracking correction vs pt-y, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt3_pty", "Tracking correction vs pt-y, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt4_pty", "Tracking correction vs pt-y, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt5_pty", "Tracking correction vs pt-y, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // tracking corrections with SPD first and TPCncls>90, cos Theta*CS dependence
  {"trackingCorrectionSPDfirst_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt1_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt2_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt3_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt4_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt5_ThetaCS", "Tracking correction vs cos #theta^{*}_{CS}, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // tracking corrections with SPD first and TPCncls>90, cos Theta*HE dependence
  {"trackingCorrectionSPDfirst_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt1_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt2_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt3_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt4_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"trackingCorrectionSPDfirstPt5_ThetaHE", "Tracking correction vs cos #theta^{*}_{HE}, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // PID corrections with SPD any and TPCncls>90, pt dependence
  {"pidCorrectionSPDany_pt", "PID correction vs pt, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt1_pt", "PID correction vs pt, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt2_pt", "PID correction vs pt, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt3_pt", "PID correction vs pt, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt4_pt", "PID correction vs pt, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt5_pt", "PID correction vs pt, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // PID corrections with SPD any and TPCncls>90, y dependence
  {"pidCorrectionSPDany_y", "PID correction vs y, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt1_y", "PID correction vs y, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt2_y", "PID correction vs y, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt3_y", "PID correction vs y, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt4_y", "PID correction vs y, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt5_y", "PID correction vs y, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // PID corrections with SPD any and TPCncls>90, pt-y dependence
  {"pidCorrectionSPDany_pty", "PID correction vs pt-y, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt1_pty", "PID correction vs pt-y, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt2_pty", "PID correction vs pt-y, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt3_pty", "PID correction vs pt-y, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt4_pty", "PID correction vs pt-y, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt5_pty", "PID correction vs pt-y, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // PID corrections with SPD any and TPCncls>90, cos Theta*CS dependence
  {"pidCorrectionSPDany_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt1_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt2_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt3_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt4_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt5_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // PID corrections with SPD any and TPCncls>90, cos Theta*HE dependence
  {"pidCorrectionSPDany_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt1_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<0.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt2_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 0.8<p_{T}(J/#Psi)<1.4, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt3_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 1.4<p_{T}(J/#Psi)<2.8, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt4_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 2.8<p_{T}(J/#Psi)<5.0, with SPD any and TPCncls>90"},
  {"pidCorrectionSPDanyPt5_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 5.0<p_{T}(J/#Psi)<10.0, with SPD any and TPCncls>90"},
  // PID corrections with SPD first and TPCncls>90, pt dependence
  {"pidCorrectionSPDfirst_pt", "PID correction vs pt, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt1_pt", "PID correction vs pt, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt2_pt", "PID correction vs pt, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt3_pt", "PID correction vs pt, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt4_pt", "PID correction vs pt, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt5_pt", "PID correction vs pt, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // PID corrections with SPD first and TPCncls>90, y dependence
  {"pidCorrectionSPDfirst_y", "PID correction vs y, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt1_y", "PID correction vs y, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt2_y", "PID correction vs y, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt3_y", "PID correction vs y, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt4_y", "PID correction vs y, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt5_y", "PID correction vs y, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // PID corrections with SPD first and TPCncls>90, pt-y dependence
  {"pidCorrectionSPDfirst_pty", "PID correction vs pt-y, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt1_pty", "PID correction vs pt-y, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt2_pty", "PID correction vs pt-y, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt3_pty", "PID correction vs pt-y, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt4_pty", "PID correction vs pt-y, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt5_pty", "PID correction vs pt-y, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // PID corrections with SPD first and TPCncls>90, cos Theta*CS dependence
  {"pidCorrectionSPDfirst_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt1_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt2_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt3_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt4_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt5_ThetaCS", "PID correction vs cos #theta^{*}_{CS}, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  // PID corrections with SPD first and TPCncls>90, cos Theta*HE dependence
  {"pidCorrectionSPDfirst_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt1_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 0.0<p_{T}(J/#Psi)<0.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt2_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 0.8<p_{T}(J/#Psi)<1.4, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt3_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 1.4<p_{T}(J/#Psi)<2.8, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt4_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 2.8<p_{T}(J/#Psi)<5.0, with SPD first and TPCncls>90"},
  {"pidCorrectionSPDfirstPt5_ThetaHE", "PID correction vs cos #theta^{*}_{HE}, 5.0<p_{T}(J/#Psi)<10.0, with SPD first and TPCncls>90"}
};


// Function prototypes ---------------------------------------------------
// The user must modify the DefineHistograms() and FillHistograms() functions
// according to need
Double_t* GetBinning(AliCFContainer* cont, Int_t variable, Int_t& nBins);
void FillHistograms(TObjArray* histosArray, AliCFContainer* cont, Int_t currentRangeStep, Bool_t firstTime);
void GetBinLimits(AliCFContainer* cont);
void DefineHistograms(TObjArray* objArray, Int_t iCutSet);
void AddHistogram(TObjArray* objArray, Int_t ndim, 
		  const Char_t* name, const Char_t* title, 
		  Int_t nbinsx, Double_t* binsx, const Char_t* xLabel = "",
		  Int_t nbinsy=0, Double_t* binsy=0, const Char_t* yLabel = "",
		  Int_t nbinsz=0, Double_t* binsz=0, const Char_t* zLabel = "");
void ProjectManyRuns(const Char_t* runList, Int_t howMany=1, Int_t offset = 0);
void ProjectAll(const Char_t* inputList, const Char_t* outfilename="HistosFromCFs.root", 
		Int_t howMany=1, Int_t offset=0);
void ExtractEfficienciesMany(const Char_t* runList, Int_t howMany=1, Int_t offset=0);
void ExtractEfficiencies(const Char_t* inputFile, const Char_t* outfilename="Efficiencies.root", const Char_t* numbersFile="");
TH1* DivideHists(TH1* nominator, TH1* denominator);
//-------------------------------------------------------------------------


//_______________________________________________________________________________________
void ProjectManyRuns(const Char_t* runList, Int_t howMany, Int_t offset) {
  //
  //
  //
  
  // loop over all runs -----------------------
  ifstream input; input.open(runList);
  Int_t runCounter = 0;
  while(input.good()) {
    Char_t readString[256];
    input.getline(readString, 256, '\n');  // get a chunk
    TString runStr = readString;
    Int_t run = runStr.Atoi();
    if(run<=0) continue;

    if(runCounter<offset) {
      runCounter++;
      continue;
    }
    if(runCounter>=offset+howMany) 
      break;

    cout << "=================== run " << run << " ============================" << endl;

    ProjectAll(Form("LHC10f7a/invMass_BB1/%s/listCF.txt",readString), 
	       Form("LHC10f7a/invMass_BB1/%s/Projections.root",readString), 
	       100, 0);
    runCounter++;
  }
}


//_______________________________________________________________________________________
void ProjectAll(const Char_t* inputList, 
		const Char_t* outfilename, 
		Int_t howMany, Int_t offset) {
  //
  //  Main function for making projections from a list of CF containers (inputList).
  //  The resulting histograms are placed in the ROOT file specified by outfilename
  //
  //  Modify the global variables above to match your requirements
  //

  // open the output file
  TFile *outFile = new TFile(outfilename,"RECREATE");
  // -----------------------------------------------------------------------------

  // copy the current ExtractEfficiency macro in the same dir as the output file
  TString outStr = "";
  outStr += outfilename;
  outStr.ReplaceAll(".root", "_ExtractEfficienciesMacro.C");
  gSystem->Exec(Form("cp ExtractEfficiencies.C %s", outStr.Data()));
  // ---------------------------------------------------------------------------
  
  // create the container for all the histograms ---------------------
  TObjArray *histoArray=new TObjArray();
  histoArray->SetOwner();
  //------------------------------------------------------------------

  
  // loop over all CF files, project and merge -----------------------
  ifstream input; input.open(inputList);
  Int_t currentFile=0;
  Bool_t firstTime = kTRUE;
  while(input.good()) {
    Char_t readString[256];
    input.getline(readString, 256, '\n');  // get a chunk
    TString readStringString = readString;
    if(readStringString[0]!='/') continue;
    if(!readStringString.Contains(".root")) continue;

    if(currentFile<offset) {
      currentFile++;
      continue;
    }
    if(currentFile>=offset+howMany) 
      break;

    cout << "file: " << readString << endl;

    AliDielectronCFdraw *cf=new AliDielectronCFdraw(readString);
    AliCFContainer* cont=cf->GetCFContainer();

    // ****************************************************************************
    // Below apply all your cut sets then call the FillHistograms() function
    // Don't forget to increment the "currentCutSet" variable after every cut set
    // ****************************************************************************

    // pair type (0 ++, 1 +-, 2 --) ----------------------------------
    cf->SetRangeUser("PairType", 1, 1);
    // Pair rapidity cut
    cf->SetRangeUser("Y", -0.899, 0.899);
    Int_t currentCutSet = 0;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 0<pt<0.8 ----------------------------------------------------
    cf->SetRangeUser("Pt", 0.001, 0.799);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 0.8<pt<1.4 ----------------------------------------------------
    cf->SetRangeUser("Pt", 0.801, 1.399);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 1.4<pt<2.8 ----------------------------------------------------
    cf->SetRangeUser("Pt", 1.401, 2.799);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 2.8<pt<5.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 2.801, 4.999);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 5.0<pt<10.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 5.001, 9.999);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // Leg pseudo-rapidity cut ---------------------------------------
    cf->SetRangeUser("Pt", 0.001, 9.999);
    cf->SetRangeUser("Leg1_Eta", -0.899, 0.899);
    cf->SetRangeUser("Leg2_Eta", -0.899, 0.899);
    cf->SetRangeUser("Leg1_Pt", 0.801, 10.0);
    cf->SetRangeUser("Leg2_Pt", 0.801, 10.0);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 0<pt<0.8 ----------------------------------------------------
    cf->SetRangeUser("Pt", 0.001, 0.799);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 0.8<pt<1.4 ----------------------------------------------------
    cf->SetRangeUser("Pt", 0.801, 1.399);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 1.4<pt<2.8 ----------------------------------------------------
    cf->SetRangeUser("Pt", 1.401, 2.799);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 2.8<pt<5.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 2.801, 4.999);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 5.0<pt<10.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 5.001, 9.999);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // |LegEta|<0.88 & Leg_Pt>1.0 & TPCncls>90 ------------------------------------------
    cf->SetRangeUser("Pt", 0.001, 9.999);
    cf->SetRangeUser("Leg1_NclsTPC", 90.1, 160.0);
    cf->SetRangeUser("Leg2_NclsTPC", 90.1, 160.0);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 0<pt<0.8 ----------------------------------------------------
    cf->SetRangeUser("Pt", 0.001, 0.799);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 0.8<pt<1.4 ----------------------------------------------------
    cf->SetRangeUser("Pt", 0.801, 1.399);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 1.4<pt<2.8 ----------------------------------------------------
    cf->SetRangeUser("Pt", 1.401, 2.799);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 2.8<pt<5.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 2.801, 4.999);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);

    // j/psi 5.0<pt<10.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 5.001, 9.999);
    currentCutSet++;
    FillHistograms(histoArray,cont,currentCutSet, firstTime);
    

    currentFile++;
    firstTime = kFALSE;
    delete cont;
    delete cf;
  }  // end loop over CF files

  outFile->cd();  
  histoArray->Write();
  outFile->Close();
  delete histoArray;
  return;
}

//_______________________________________________________________________________________
void ExtractEfficienciesMany(const Char_t* runList, Int_t howMany, Int_t offset) {
  //
  //
  //

  // loop over all runs -----------------------
  ifstream input; input.open(runList);
  Int_t runCounter = 0;

  TGraphErrors* trends[gkNeffs];
  Double_t weightedEffs[gkNeffs];
  Double_t weightedErrs[gkNeffs];
  Double_t nTotalEvents = 0;
  for(Int_t iTrend=0; iTrend<gkNeffs; iTrend++) {
    trends[iTrend] = new TGraphErrors();
    trends[iTrend]->SetName(gkEffNames[iTrend][0]);
    trends[iTrend]->SetTitle(gkEffNames[iTrend][1]);
    weightedEffs[iTrend] = 0.0; weightedErrs[iTrend] = 0.0;
  }
  TFile* file=0x0;
  TFile* normalizationFile=0x0;
  TNamed* object;
  while(input.good()) {
    Char_t readString[256];
    input.getline(readString, 256, '\n');  // get a chunk
    TString runStr = readString;
    Int_t run = runStr.Atoi();
    if(run<=0) continue;

    if(runCounter<offset) {
      runCounter++;
      continue;
    }
    if(runCounter>=offset+howMany) 
      break;

    cout << "=================== run " << run << " ============================" << endl;
    TString periodStr;
    if(run<=117222) periodStr = "LHC10b.pass2";
    if(run>117222 && run<=120829) periodStr = "LHC10c.pass2";
    if(run>=122374 && run<=126437) periodStr = "LHC10d.pass1";
    Double_t nPhysicsEvents = 0;
    normalizationFile = TFile::Open(Form("/u/iarsene/work/ALICE/normalization/2010-10-03_0628.3293/%s/%d/iarsene_normalization.root", periodStr.Data(), run));
    if(normalizationFile) {
      TObjArray *histos=(TObjArray*)normalizationFile->Get("iarsene_normalization");
      TH1I* triggers=(TH1I*)histos->FindObject("TriggersHistogram");
      nPhysicsEvents = triggers->GetBinContent(2);   // PHYSICS events
      normalizationFile->Close();
      nTotalEvents += nPhysicsEvents;
    }

    ExtractEfficiencies(Form("LHC10f7a/iter10_BB2/%s/Projections_iter10.root",readString), 
			Form("LHC10f7a/iter10_BB2/%s/Efficiencies_iter10.root",readString));
    file = TFile::Open(Form("LHC10f7a/iter10_BB2/%s/Efficiencies_iter10.root",readString));
    if(file && !file->IsZombie()) {   
      for(Int_t iTrend=0; iTrend<gkNeffs; iTrend++) {
	object = (TNamed*)file->Get(Form("%s_value",gkEffNames[iTrend][0]));
	if(!object) continue;
	Float_t eff = (TString(object->GetTitle())).Atof();
	trends[iTrend]->SetPoint(runCounter, run, eff);
	object = (TNamed*)file->Get(Form("%s_error",gkEffNames[iTrend][0]));
	if(!object) continue;
	Float_t err = (TString(object->GetTitle())).Atof();
	trends[iTrend]->SetPointError(runCounter, 0.0, err);
	weightedEffs[iTrend] += nPhysicsEvents*eff;
	weightedErrs[iTrend] += nPhysicsEvents*nPhysicsEvents*err*err;
      }
      file->Close();
    }
    if(normalizationFile)
      normalizationFile->Close();

    runCounter++;
  }

  TFile *saveTrend = new TFile(Form("%s.trend_iter10.root", runList), "RECREATE");
  TNamed *weightedFactors;
  TNamed *weightedErrors;
  for(Int_t iTrend=0; iTrend<gkNeffs; iTrend++) {
    trends[iTrend]->Write();
    weightedEffs[iTrend] /= nTotalEvents;
    weightedErrs[iTrend] = TMath::Sqrt(weightedErrs[gkNeffs]/nTotalEvents/nTotalEvents);
    weightedFactors = new TNamed(Form("%s_weighted", gkEffNames[iTrend][0]),
				 Form("%f", weightedEffs[iTrend]));
    weightedErrors = new TNamed(Form("%s_weightedErr", gkEffNames[iTrend][0]),
				 Form("%f", weightedErrs[iTrend]));
    weightedFactors->Write();
    weightedErrors->Write();
  }
  weightedFactors = new TNamed("TotalEvents", Form("%f",nTotalEvents));
  weightedFactors->Write();

  saveTrend->Close();
}

//_______________________________________________________________________________________
void ExtractEfficiencies(const Char_t* inputFilename,
			 const Char_t* outfilename,
			 const Char_t* numbersFile) {
  //
  // Main function to extract efficiencies
  //

  // Open the output file
  TFile *output = new TFile(outfilename, "RECREATE");
  // ---------------------------------------------------------------------------

  // copy the current ExtractEfficiency macro in the same dir as the output file
  TString outStr = "";
  outStr += outfilename;
  outStr.ReplaceAll(".root", "_ExtractEfficienciesMacro.C");
  gSystem->Exec(Form("cp ExtractEfficiencies.C %s", outStr.Data()));
  // ---------------------------------------------------------------------------

  // open the input file and read out all the histograms
  TFile *file = TFile::Open(inputFilename);
  if(!file || file->IsZombie()) return;
  
  TObjArray *effArray = new TObjArray();
  effArray->SetOwner();
  TH1* nominator;
  TH1* denominator;
  TNamed* effValue;
  TNamed* effError;
  ofstream asciiOut;
  if(numbersFile[0]!='\0') {
    asciiOut.open(numbersFile);
    asciiOut << "#Format:  Name  |  Value   |   Abs. Error" << endl;
  }
  for(Int_t iEff=0; iEff<gkNeffs; iEff++) {
    cout << gkEffNames[iEff][0] << " (" << gkEffNames[iEff][1] << " )" << endl;
    if(gkDims[gkEffs[iEff][2]][0]==1) {      // 1-dim histos
      nominator = (TH1D*)(file->Get(Form("%s_%s_%s",gkStepNames[gkEffs[iEff][0]][0],
					  gkCutSetNames[gkEffs[iEff][1]][0],
					  gkHistoNames[gkEffs[iEff][2]][0])));
      denominator = (TH1D*)(file->Get(Form("%s_%s_%s",gkStepNames[gkEffs[iEff][3]][0],
					    gkCutSetNames[gkEffs[iEff][4]][0],
					    gkHistoNames[gkEffs[iEff][5]][0])));
      if(!nominator) continue;
      if(!denominator) continue;
      nominator->GetYaxis()->SetTitle("efficiency");
    }
    if(gkDims[gkEffs[iEff][2]][0]==2) {      // 2-dim histos
      nominator = (TH2D*)(file->Get(Form("%s_%s_%s",gkStepNames[gkEffs[iEff][0]][0],
					  gkCutSetNames[gkEffs[iEff][1]][0],
					  gkHistoNames[gkEffs[iEff][2]][0])));
      denominator = (TH2D*)(file->Get(Form("%s_%s_%s",gkStepNames[gkEffs[iEff][3]][0],
					    gkCutSetNames[gkEffs[iEff][4]][0],
					    gkHistoNames[gkEffs[iEff][5]][0])));
      if(!nominator) continue;
      if(!denominator) continue;
      nominator->GetZaxis()->SetTitle("efficiency");
    }
    if(gkDims[gkEffs[iEff][2]][0]==3) {      // 3-dim histos
      nominator = (TH3D*)(file->Get(Form("%s_%s_%s",gkStepNames[gkEffs[iEff][0]][0],
					  gkCutSetNames[gkEffs[iEff][1]][0],
					  gkHistoNames[gkEffs[iEff][2]][0])));
      denominator = (TH3D*)(file->Get(Form("%s_%s_%s",gkStepNames[gkEffs[iEff][3]][0],
					    gkCutSetNames[gkEffs[iEff][4]][0],
					    gkHistoNames[gkEffs[iEff][5]][0])));
      if(!nominator) continue;
      if(!denominator) continue;
    }
    Double_t nomIntegral = nominator->Integral();
    Double_t denomIntegral = denominator->Integral();
    Double_t eff = (denomIntegral>0 ? nomIntegral/denomIntegral : 0);
    Double_t error = (nomIntegral>0 && denomIntegral>0 ? eff*TMath::Sqrt(1.0/nomIntegral + 1.0/denomIntegral) : 0);
    //nominator->Divide(denominator);
    TH1* ratio = DivideHists(nominator, denominator);
    cout << "efficiency = " << nomIntegral << " / " << denomIntegral << " = "
	 << eff << " +/- " << error << endl;
    TString title = gkEffNames[iEff][1];
    title += Form(", integrated eff. = %f #pm %f", eff, error);
    ratio->SetTitle(title.Data());
    ratio->SetName(gkEffNames[iEff][0]);
    effArray->Add(ratio);

    if(numbersFile[0]!='\0') {
      asciiOut << gkEffNames[iEff][0] << "\t" << eff << "\t" << error << endl;
    }
    effValue = new TNamed(Form("%s_value", gkEffNames[iEff][0]),
			  Form("%f", eff));
    effArray->Add(effValue);
    effError = new TNamed(Form("%s_error", gkEffNames[iEff][0]),
			  Form("%f", error));
    effArray->Add(effError);
  }

  output->cd();
  effArray->Write();
  output->Close();
  file->Close();
  asciiOut.close();
}


//________________________________________________________________________________________
void DefineHistograms(TObjArray* histoArray, Int_t iCutSet) {
  //
  // Define the histograms to be filled for every step and a given cut set
  // This function is called by the FillHistograms() if the firstTime flag is set

  for(Int_t iStep=0; iStep<kNSteps; iStep++) {
    for(Int_t iHisto = 0; iHisto<gkNhistos; iHisto++) {
      AddHistogram(histoArray, gkDims[iHisto][0],
		   Form("%s_%s_%s", gkStepNames[iStep][0], gkCutSetNames[iCutSet][0], gkHistoNames[iHisto][0]),
		   Form("%s, %s, %s", gkHistoNames[iHisto][1], gkStepNames[iStep][1], gkCutSetNames[iCutSet][1]),
		   gNbins[gkDims[iHisto][1]], gBinLimits[gkDims[iHisto][1]], gkVarNames[gkDims[iHisto][1]],
		   (gkDims[iHisto][2]!=kNothing ? gNbins[gkDims[iHisto][2]] : 0), 
		   (gkDims[iHisto][2]!=kNothing ? gBinLimits[gkDims[iHisto][2]] : 0),
		   (gkDims[iHisto][2]!=kNothing ? gkVarNames[gkDims[iHisto][2]] : ""),
		   (gkDims[iHisto][3]!=kNothing ? gNbins[gkDims[iHisto][3]] : 0), 
		   (gkDims[iHisto][3]!=kNothing ? gBinLimits[gkDims[iHisto][3]] : 0),
		   (gkDims[iHisto][3]!=kNothing ? gkVarNames[gkDims[iHisto][3]] : ""));
      
    }  // end loop over histos
  }   // end loop over steps
}


//_________________________________________________________________________________________
void AddHistogram(TObjArray* histoArray, Int_t ndim, 
		  const Char_t* name, const Char_t* title, 
		  Int_t nbinsx, Double_t* binsx, const Char_t* xLabel,
		  Int_t nbinsy, Double_t* binsy, const Char_t* yLabel,
		  Int_t nbinsz, Double_t* binsz, const Char_t* zLabel) {
  //
  // Create a 1,2 or 3 - dimensional histogram and add it to the object array
  //
  if(ndim<1 || ndim>3) return;
  TH1* histo;
  if(ndim==1) {
    histo = new TH1D(name, title, nbinsx, binsx);
    histo->Sumw2();
    histo->GetXaxis()->SetTitle(xLabel);
  }
  if(ndim==2) {
    histo = new TH2D(name, title, nbinsx, binsx, nbinsy, binsy);
    histo->Sumw2();
    histo->GetXaxis()->SetTitle(xLabel);
    histo->GetYaxis()->SetTitle(yLabel);
  }
  if(ndim==3) {
    histo = new TH3D(name, title, nbinsx, binsx, nbinsy, binsy, nbinsz, binsz);
    histo->Sumw2();
    histo->GetXaxis()->SetTitle(xLabel);
    histo->GetYaxis()->SetTitle(yLabel);
    histo->GetZaxis()->SetTitle(zLabel);
  }
  histoArray->Add(histo);
}

//__________________________________________________________________________________________
void FillHistograms(TObjArray* histosArray, AliCFContainer* cont, Int_t currentCutSet, Bool_t firstTime) {
  //
  // Fill the user defined histograms for a given cut set
  // 
  // If the firstTime flag is on then update the bin limits and call DefineHistograms()
  if(firstTime) {
    GetBinLimits(cont);
    DefineHistograms(histosArray, currentCutSet);
  }

  TH1* histo;
  for(Int_t iStep=0; iStep<kNSteps; ++iStep) {  // loop over CF container steps
    for(Int_t iHisto=0; iHisto<gkNhistos; iHisto++) {
      // fill 1-dim histos
      if(gkDims[iHisto][0]==1) {
	histo = (TH1D*)histosArray->FindObject(Form("%s_%s_%s",gkStepNames[iStep][0],
						    gkCutSetNames[currentCutSet][0],
						    gkHistoNames[iHisto][0]));
	histo->Add(cont->Project(gkDims[iHisto][1],gkStepNumbers[iStep]));
      }
      // fill 2-dim histos
      if(gkDims[iHisto][0]==2) {
	histo = (TH2D*)histosArray->FindObject(Form("%s_%s_%s",gkStepNames[iStep][0],
						    gkCutSetNames[currentCutSet][0],
						    gkHistoNames[iHisto][0]));
	histo->Add(cont->Project(gkDims[iHisto][1], gkDims[iHisto][2], gkStepNumbers[iStep]));
      }
      // fill 3-dim histos
      if(gkDims[iHisto][0]==3) {
	histo = (TH3D*)histosArray->FindObject(Form("%s_%s_%s",gkStepNames[iStep][0],
						    gkCutSetNames[currentCutSet][0],
						    gkHistoNames[iHisto][0]));
	histo->Add(cont->Project(gkDims[iHisto][1], gkDims[iHisto][2], gkDims[iHisto][3], gkStepNumbers[iStep]));
      }
    }   // end loop over histos
  }  // end loop over steps
}

//____________________________________________________________________________________________
void GetBinLimits(AliCFContainer* cont) {
  //
  // Extract the bin limits from the CF container
  //
  cout << "********* New cut set ****************" << endl;
  for(Int_t iVar=0; iVar<kNVariables; iVar++) {
    gNbins[iVar] = 0;
    gBinLimits[iVar] = GetBinning(cont, iVar, gNbins[iVar]);
    cout << "n bins on " << cont->GetVarTitle(iVar) << " : " << gNbins[iVar];
    cout << "; range = " << gBinLimits[iVar][0] << " --> " << gBinLimits[iVar][gNbins[iVar]] << endl;
  }
}

//________________________________________________________________________________________
Double_t* GetBinning(AliCFContainer* cont, Int_t variable, 
		     Int_t& nBins) {
  //
  // Get the number of bins and the bin limits for the projection of a given variable
  //
  TH1D* tempHist = cont->Project(variable, kPureMC);
  nBins = tempHist->GetXaxis()->GetNbins();
  Double_t* binLimits = new Double_t[nBins+1];
  for(Int_t i=1; i<=nBins; i++)
    binLimits[i-1]=tempHist->GetXaxis()->GetBinLowEdge(i);
  binLimits[nBins] = tempHist->GetXaxis()->GetBinLowEdge(nBins) + 
    tempHist->GetXaxis()->GetBinWidth(nBins);
  return binLimits;
}

//________________________________________________________________________________________
TH1* DivideHists(TH1* nominator, TH1* denominator) {
  //
  // divide 2 histograms with error propagation
  //
  TH1* ratio;
  if(nominator->InheritsFrom("TH3")) {
    Int_t nBinsXNom = nominator->GetXaxis()->GetNbins();
    Int_t nBinsXDenom = denominator->GetXaxis()->GetNbins();
    Int_t nBinsYNom = nominator->GetYaxis()->GetNbins();
    Int_t nBinsYDenom = denominator->GetYaxis()->GetNbins();
    Int_t nBinsZNom = nominator->GetZaxis()->GetNbins();
    Int_t nBinsZDenom = denominator->GetZaxis()->GetNbins();
    if(nBinsXNom!=nBinsXDenom || nBinsYNom!=nBinsYDenom || nBinsZNom!=nBinsZDenom) {
      cout << "Trying to divide histograms with different number of bins" << endl;
      return 0x0;
    }
    ratio = (TH3D*)nominator->Clone("ratio");
    ratio->Reset();
    for(Int_t iXbin=1; iXbin<=nBinsXNom; ++iXbin) {
      for(Int_t iYbin=1; iYbin<=nBinsYNom; ++iYbin) {
	for(Int_t iZbin=1; iZbin<=nBinsZNom; ++iZbin) {
	  Double_t countsN = nominator->GetBinContent(iXbin, iYbin, iZbin);
	  Double_t countsD = denominator->GetBinContent(iXbin, iYbin, iZbin);
	  if(countsN<1 || countsD<1) continue;    // zero entry bins
	  ratio->SetBinContent(iXbin, iYbin, iZbin, countsN/countsD);
	  ratio->SetBinError(iXbin, iYbin, iZbin, (countsN/countsD)*TMath::Sqrt(1.0/countsN)+(1.0/countsD));
	}
      }
    }
    return ratio;
  }

  if(nominator->InheritsFrom("TH2")) {
    Int_t nBinsXNom = nominator->GetXaxis()->GetNbins();
    Int_t nBinsXDenom = denominator->GetXaxis()->GetNbins();
    Int_t nBinsYNom = nominator->GetYaxis()->GetNbins();
    Int_t nBinsYDenom = denominator->GetYaxis()->GetNbins();
    if(nBinsXNom!=nBinsXDenom || nBinsYNom!=nBinsYDenom) {
      cout << "Trying to divide histograms with different number of bins" << endl;
      return 0x0;
    }
    ratio = (TH2D*)nominator->Clone("ratio");
    ratio->Reset();
    for(Int_t iXbin=1; iXbin<=nBinsXNom; ++iXbin) {
      for(Int_t iYbin=1; iYbin<=nBinsYNom; ++iYbin) {
	Double_t countsN = nominator->GetBinContent(iXbin, iYbin);
	Double_t countsD = denominator->GetBinContent(iXbin, iYbin);
	if(countsN<1 || countsD<1) continue;    // zero entry bins
	ratio->SetBinContent(iXbin, iYbin, countsN/countsD);
	ratio->SetBinError(iXbin, iYbin, (countsN/countsD)*TMath::Sqrt(1.0/countsN)+(1.0/countsD));
      }
    }
    return ratio;
  }

  if(nominator->InheritsFrom("TH1")) {
    Int_t nBinsXNom = nominator->GetXaxis()->GetNbins();
    Int_t nBinsXDenom = denominator->GetXaxis()->GetNbins();
    if(nBinsXNom!=nBinsXDenom) {
      cout << "Trying to divide histograms with different number of bins" << endl;
      return 0x0;
    }
    ratio = (TH1D*)nominator->Clone("ratio");
    ratio->Reset();
    for(Int_t iXbin=1; iXbin<=nBinsXNom; ++iXbin) {
      Double_t countsN = nominator->GetBinContent(iXbin);
      Double_t countsD = denominator->GetBinContent(iXbin);
      if(countsN<1 || countsD<1) continue;    // zero entry bins
      ratio->SetBinContent(iXbin, countsN/countsD);
      ratio->SetBinError(iXbin, (countsN/countsD)*TMath::Sqrt(1.0/countsN)+(1.0/countsD));
    }
    return ratio;
  }
    
  return 0x0;
}
