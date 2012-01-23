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
  kThetaCS,
  kThetaHE,
  kM,
  kPairType,
  //  kPhi,
  kLeg1_Pt,
  kLeg1_NclsTPC,
  kLeg1_Eta,
  //kLeg1_Phi,
  //  kLeg1_TPC_nSigma_Electrons,
  //  kLeg1_P,
  //kLeg2_Phi,
  //  kLeg2_TPC_nSigma_Electrons,
  kLeg2_Pt,
  //  kLeg2_P,
  kLeg2_NclsTPC,
  kLeg2_Eta,
  kNVariables       // kNVariables should be always here!
};
const Char_t* gkVarNames[kNVariables] = {     // variable names to be put on histograms axes and titles
  "p_{T} [GeV/c]",
  "y",
  "cos #theta^{*}_{CS}",
  "cos #theta^{*}_{HE}",
  "M [GeV/c^{2}]",
  "Pair type (0=++; 1=+-; 2=--)",
  // "#phi(J/#Psi) [rad.]",
  //"#phi^{leg1}",
  //  "TPC n #sigma electrons (leg1)",
  "p_{T}^{leg1} [GeV/c]",
  //  "P^{leg1} [GeV/c]",
  "# TPC clusters (leg1)",
  "#eta^{leg1}",
  //"#phi^{leg2}",
  //  "TPC n #sigma electrons (leg2)",
  "p_{T}^{leg2} [GeV/c]",
  //  "P^{leg2} [GeV/c]",
  "# TPC clusters (leg2)"
  "#eta^{leg2}",
};
Int_t gNbins[kNVariables];           // number of bins for every variable --> filled automatically
Double_t* gBinLimits[kNVariables];   // bin limits for every variable --> filled automatically
// ------------------------------------------------------------------------

// Put here all the CF steps of interest ----------------------------------
enum Steps {        // step indexes in the CF containers to be analyzed
  kPureMC = 0,
  kESDSPDany = 2,
  kESDSPDfirst = 4,
  kESDv0SPDany = 6,
  kESDv0SPDfirst = 8,
  kFullSPDany = 10,
  kFullSPDfirst = 12,
  kNSteps = 7        // total number of steps (the number of steps above)
};
const Int_t gkStepNumbers[kNSteps] = {   // array with step indexes (all from the enumeration above)
  kPureMC, 
  kESDSPDany, kESDSPDfirst,
  kESDv0SPDany, kESDv0SPDfirst,
  kFullSPDany, kFullSPDfirst
};
const Char_t* gkStepNames[kNSteps][2] = {// names for each CF step
  {"PureMC",         "Pure MC"},         // NOTE: short names go to histo names, long names go to titles
  {"ESDSPDany",      "ESD track cuts, SPD any"}, 
  {"ESDSPDfirst",    "ESD track cuts, SPD first"},
  {"ESDv0SPDany",    "ESD track cuts, conv. cuts, SPD any"}, 
  {"ESDv0SPDfirst",  "ESD track cuts, conv. cuts, SPD first"},
  {"FullSPDany",     "All track cuts (with SPD any) and TPC-PID"},
  {"FullSPDfirst",   "All track cuts (with SPD first) and TPC-PID"}
};
//------------------------------------------------------------------------

// Put here info about the cut sets for which projections will be made ---
const Int_t gkNCutSets = 10*3;     // number of cut sets for which histos will be filled
const Char_t* gkCutSetNames[gkNCutSets][2] = {   // short and long names for all the cut sets
  // baseline
  {"Ycut",                 "|y_{J/#Psi}|<0.9 & 0<pt_{J/#Psi}<10.0"},
  {"Ycut1",                "|y_{J/#Psi}|<0.3 & 0<pt_{J/#Psi}<10.0"},
  {"Ycut2",                "0.3<y_{J/#Psi}<0.9 & 0<pt_{J/#Psi}<10.0"},
  {"Ycut3",                "-0.9<y_{J/#Psi}<-0.3 & 0<pt_{J/#Psi}<10.0"},
  {"YcutPt1",              "|y_{J/#Psi}|<0.9 & 0.0<p_{T J/#Psi}<1.0"},
  {"YcutPt2",              "|y_{J/#Psi}|<0.9 & 1.0<p_{T J/#Psi}<2.0"},
  {"YcutPt3",              "|y_{J/#Psi}|<0.9 & 2.0<p_{T J/#Psi}<3.0"},
  {"YcutPt4",              "|y_{J/#Psi}|<0.9 & 3.0<p_{T J/#Psi}<5.0"},
  {"YcutPt5",              "|y_{J/#Psi}|<0.9 & 5.0<p_{T J/#Psi}<7.0"},
  {"YcutPt6",              "|y_{J/#Psi}|<0.9 & 7.0<p_{T J/#Psi}<10.0"},
  // track cuts
  {"YcutLegsFull",       "|y_{J/#Psi}|<0.9 & 0<pt_{J/#Psi}<10.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"Ycut1LegsFull",      "|y_{J/#Psi}|<0.3 & 0<pt_{J/#Psi}<10.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"Ycut2LegsFull",      "0.3<y_{J/#Psi}<0.9 & 0<pt_{J/#Psi}<10.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"Ycut3LegsFull",      "-0.9<y_{J/#Psi}<-0.3 & 0<pt_{J/#Psi}<10.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt1LegsFull",    "|y_{J/#Psi}|<0.9 & 0<pt_{J/#Psi}<1.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt2LegsFull",    "|y_{J/#Psi}|<0.9 & 1.0<pt_{J/#Psi}<2.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt3LegsFull",    "|y_{J/#Psi}|<0.9 & 2.0<pt_{J/#Psi}<3.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt4LegsFull",    "|y_{J/#Psi}|<0.9 & 3.0<pt_{J/#Psi}<5.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt5LegsFull",    "|y_{J/#Psi}|<0.9 & 5.0<pt_{J/#Psi}<7.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt6LegsFull",    "|y_{J/#Psi}|<0.9 & 7.0<pt_{J/#Psi}<10.0 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  // track cuts + cut on signal integration range
  {"YcutMcutLegsFull",       "|y_{J/#Psi}|<0.9 & 0<pt_{J/#Psi}<10.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"Ycut1McutLegsFull",      "|y_{J/#Psi}|<0.3 & 0<pt_{J/#Psi}<10.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"Ycut2McutLegsFull",      "0.3<y_{J/#Psi}<0.9 & 0<pt_{J/#Psi}<10.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"Ycut3McutLegsFull",      "-0.9<y_{J/#Psi}<-0.3 & 0<pt_{J/#Psi}<10.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt1McutLegsFull",    "|y_{J/#Psi}|<0.9 & 0<pt_{J/#Psi}<1.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt2McutLegsFull",    "|y_{J/#Psi}|<0.9 & 1.0<pt_{J/#Psi}<2.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt3McutLegsFull",    "|y_{J/#Psi}|<0.9 & 2.0<pt_{J/#Psi}<3.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt4McutLegsFull",    "|y_{J/#Psi}|<0.9 & 3.0<pt_{J/#Psi}<5.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt5McutLegsFull",    "|y_{J/#Psi}|<0.9 & 5.0<pt_{J/#Psi}<7.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"},
  {"YcutPt6McutLegsFull",    "|y_{J/#Psi}|<0.9 & 7.0<pt_{J/#Psi}<10.0 & 2.92<M_{inv}<3.16 & |#eta_{legs}|<0.9 & p_{Tlegs}>1.0 GeV/c & Ncls_TPC>70"}
  
};
// -----------------------------------------------------------------------

// Put here info about the histograms to be filled -----------------------
const Int_t gkNhistos = 3;                        // how many histograms for every (step,cut set) combination
const Char_t* gkHistoNames[gkNhistos][2] = {      // short and long names of the histograms
  {"pt","p_{T}(J/#Psi)"},                         // NOTE: short names go to the histo name, long name goes to title
  //  {"y","y(J/#Psi)"}, 
  // {"pty","p_{T} vs y(J/#Psi)"},
  //  {"phi","#phi(J/#Psi) [rad.]"}, 
  //{"m","e^{+}e^{-} invariant mass"}, 
  {"ThetaCS","cos #theta^{*}_{CS}"}, 
  {"ThetaHE","cos #theta^{*}_{HE}"}
  // {"Minv", "Invariant mass"}
};
const Int_t gkDims[gkNhistos][4] = {      // dimensions and variables for histograms
// ndim  xVar      yVar      zVar
  {1,    kPt,      kNothing, kNothing},   // pt dependence
  //  {1,    kY,       kNothing, kNothing},   // y dependence
  //  {2,    kY,       kPt,      kNothing},   // pt,y dependence
  //  {1,    kPhi,     kNothing, kNothing},   // phi dependence
  //  {1,    kM,       kNothing, kNothing},   // inv. mass dependence
  {1,    kThetaCS, kNothing, kNothing},   // cos theta* CS dependence
  {1,    kThetaHE, kNothing, kNothing}    // cos theta* HE dependence
  //{1,    kM,       kNothing, kNothing}    // invariant mass
};
// -----------------------------------------------------------------------

// ******************************************************************************** 
// Define here all the efficiencies you want (if any)
// Efficiency = (nominator histogram)/(denominator histogram)
// A histogram is defined by its step,cut set, and number defined according the
// global variables above
// ********************************************************************************
const Int_t gkNeffs = 20*3*5;
const Int_t gkEffs[gkNeffs][6] = {
  //nominator: Step  Cut  Histo | denominator: Step  Cut  Histo     comment
  // full corrections, pt dependence
             { 5,    10,   0,                   0,    0,   0     },     // full correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 6,    10,   0,                   0,    0,   0     },     // full correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 5,    11,   0,                   0,    1,   0     },     // full correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    11,   0,                   0,    1,   0     },     // full correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    12,   0,                   0,    2,   0     },     // full correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    12,   0,                   0,    2,   0     },     // full correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    13,   0,                   0,    3,   0     },     // full correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    13,   0,                   0,    3,   0     },     // full correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    14,   0,                   0,    4,   0     },     // full correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 6,    14,   0,                   0,    4,   0     },     // full correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 5,    15,   0,                   0,    5,   0     },     // full correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 6,    15,   0,                   0,    5,   0     },     // full correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 5,    16,   0,                   0,    6,   0     },     // full correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 6,    16,   0,                   0,    6,   0     },     // full correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 5,    17,   0,                   0,    7,   0     },     // full correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 6,    17,   0,                   0,    7,   0     },     // full correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 5,    18,   0,                   0,    8,   0     },     // full correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 6,    18,   0,                   0,    8,   0     },     // full correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 5,    19,   0,                   0,    9,   0     },     // full correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 6,    19,   0,                   0,    9,   0     },     // full correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0

	     // full corrections, cosThetaCS dependence
             { 5,    10,   1,                   0,    0,   1     },     // full correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 6,    10,   1,                   0,    0,   1     },     // full correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 5,    11,   1,                   0,    1,   1     },     // full correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    11,   1,                   0,    1,   1     },     // full correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    12,   1,                   0,    2,   1     },     // full correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    12,   1,                   0,    2,   1     },     // full correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    13,   1,                   0,    3,   1     },     // full correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    13,   1,                   0,    3,   1     },     // full correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    14,   1,                   0,    1,   1     },     // full correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 6,    14,   1,                   0,    1,   1     },     // full correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 5,    15,   1,                   0,    2,   1     },     // full correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 6,    15,   1,                   0,    2,   1     },     // full correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 5,    16,   1,                   0,    3,   1     },     // full correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 6,    16,   1,                   0,    3,   1     },     // full correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 5,    17,   1,                   0,    4,   1     },     // full correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 6,    17,   1,                   0,    4,   1     },     // full correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 5,    18,   1,                   0,    5,   1     },     // full correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 6,    18,   1,                   0,    5,   1     },     // full correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 5,    19,   1,                   0,    6,   1     },     // full correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 6,    19,   1,                   0,    6,   1     },     // full correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0

	     // full corrections, cosThetaHE dependence
	     { 5,    10,   2,                   0,    0,   2     },     // full correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 6,    10,   2,                   0,    0,   2     },     // full correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 5,    11,   2,                   0,    1,   2     },     // full correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    11,   2,                   0,    1,   2     },     // full correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    12,   2,                   0,    2,   2     },     // full correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    12,   2,                   0,    2,   2     },     // full correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    13,   2,                   0,    3,   2     },     // full correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    13,   2,                   0,    3,   2     },     // full correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    14,   2,                   0,    4,   2     },     // full correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 6,    14,   2,                   0,    4,   2     },     // full correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 5,    15,   2,                   0,    5,   2     },     // full correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 6,    15,   2,                   0,    5,   2     },     // full correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 5,    16,   2,                   0,    6,   2     },     // full correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 6,    16,   2,                   0,    6,   2     },     // full correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 5,    17,   2,                   0,    7,   2     },     // full correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 6,    17,   2,                   0,    7,   2     },     // full correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 5,    18,   2,                   0,    8,   2     },     // full correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 6,    18,   2,                   0,    8,   2     },     // full correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 5,    19,   2,                   0,    9,   2     },     // full correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 6,    19,   2,                   0,    9,   2     },     // full correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0             


	     // tracking corrections, pt dependence
             { 1,    10,   0,                   0,    0,   0     },     // tracking correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 2,    10,   0,                   0,    0,   0     },     // tracking correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 1,    11,   0,                   0,    1,   0     },     // tracking correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 2,    11,   0,                   0,    1,   0     },     // tracking correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 1,    12,   0,                   0,    2,   0     },     // tracking correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 2,    12,   0,                   0,    2,   0     },     // tracking correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 1,    13,   0,                   0,    3,   0     },     // tracking correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 2,    13,   0,                   0,    3,   0     },     // tracking correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 1,    14,   0,                   0,    4,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 2,    14,   0,                   0,    4,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 1,    15,   0,                   0,    5,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 2,    15,   0,                   0,    5,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 1,    16,   0,                   0,    6,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 2,    16,   0,                   0,    6,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 1,    17,   0,                   0,    7,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 2,    17,   0,                   0,    7,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 1,    18,   0,                   0,    8,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 2,    18,   0,                   0,    8,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 1,    19,   0,                   0,    9,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 2,    19,   0,                   0,    9,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0

	     // tracking corrections, cosThetaCS dependence
             { 1,    10,   1,                   0,    0,   1     },     // tracking correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 2,    10,   1,                   0,    0,   1     },     // tracking correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 1,    11,   1,                   0,    1,   1     },     // tracking correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 2,    11,   1,                   0,    1,   1     },     // tracking correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 1,    12,   1,                   0,    2,   1     },     // tracking correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 2,    12,   1,                   0,    2,   1     },     // tracking correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 1,    13,   1,                   0,    3,   1     },     // tracking correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 2,    13,   1,                   0,    3,   1     },     // tracking correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 1,    14,   1,                   0,    4,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 2,    14,   1,                   0,    4,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 1,    15,   1,                   0,    5,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 2,    15,   1,                   0,    5,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 1,    16,   1,                   0,    6,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 2,    16,   1,                   0,    6,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 1,    17,   1,                   0,    7,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 2,    17,   1,                   0,    7,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 1,    18,   1,                   0,    8,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 2,    18,   1,                   0,    8,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 1,    19,   1,                   0,    9,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 2,    19,   1,                   0,    9,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0

	     // tracking corrections, cosThetaHE dependence
             { 1,    10,   2,                   0,    0,   2     },     // tracking correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 2,    10,   2,                   0,    0,   2     },     // tracking correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 1,    11,   2,                   0,    1,   2     },     // tracking correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 2,    11,   2,                   0,    1,   2     },     // tracking correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 1,    12,   2,                   0,    2,   2     },     // tracking correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 2,    12,   2,                   0,    2,   2     },     // tracking correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 1,    13,   2,                   0,    3,   2     },     // tracking correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 2,    13,   2,                   0,    3,   2     },     // tracking correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 1,    14,   2,                   0,    4,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 2,    14,   2,                   0,    4,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 1,    15,   2,                   0,    5,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 2,    15,   2,                   0,    5,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 1,    16,   2,                   0,    6,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 2,    16,   2,                   0,    6,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 1,    17,   2,                   0,    7,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 2,    17,   2,                   0,    7,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 1,    18,   2,                   0,    8,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 2,    18,   2,                   0,    8,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 1,    19,   2,                   0,    9,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 2,    19,   2,                   0,    9,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0


	     // tracking corrections with V0 rejection cut, pt dependence
             { 3,    10,   0,                   0,    0,   0     },     // tracking correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 4,    10,   0,                   0,    0,   0     },     // tracking correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 3,    11,   0,                   0,    1,   0     },     // tracking correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 4,    11,   0,                   0,    1,   0     },     // tracking correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 3,    12,   0,                   0,    2,   0     },     // tracking correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 4,    12,   0,                   0,    2,   0     },     // tracking correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 3,    13,   0,                   0,    3,   0     },     // tracking correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 4,    13,   0,                   0,    3,   0     },     // tracking correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 3,    14,   0,                   0,    4,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 4,    14,   0,                   0,    4,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 3,    15,   0,                   0,    5,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 4,    15,   0,                   0,    5,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 3,    16,   0,                   0,    6,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 4,    16,   0,                   0,    6,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 3,    17,   0,                   0,    7,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 4,    17,   0,                   0,    7,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 3,    18,   0,                   0,    8,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 4,    18,   0,                   0,    8,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 3,    19,   0,                   0,    9,   0     },     // tracking correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 4,    19,   0,                   0,    9,   0     },     // tracking correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0

	     // tracking corrections with V0 rejection cut, cosThetaCS dependence
             { 3,    10,   1,                   0,    0,   1     },     // tracking correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 4,    10,   1,                   0,    0,   1     },     // tracking correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 3,    11,   1,                   0,    1,   1     },     // tracking correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 4,    11,   1,                   0,    1,   1     },     // tracking correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 3,    12,   1,                   0,    2,   1     },     // tracking correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 4,    12,   1,                   0,    2,   1     },     // tracking correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 3,    13,   1,                   0,    3,   1     },     // tracking correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 4,    13,   1,                   0,    3,   1     },     // tracking correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 3,    14,   1,                   0,    4,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 4,    14,   1,                   0,    4,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 3,    15,   1,                   0,    5,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 4,    15,   1,                   0,    5,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 3,    16,   1,                   0,    6,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 4,    16,   1,                   0,    6,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 3,    17,   1,                   0,    7,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 4,    17,   1,                   0,    7,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 3,    18,   1,                   0,    8,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 4,    18,   1,                   0,    8,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 3,    19,   1,                   0,    9,   1     },     // tracking correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 4,    19,   1,                   0,    9,   1     },     // tracking correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0

	     // tracking corrections with V0 rejection cut, cosThetaHE dependence
             { 3,    10,   2,                   0,    0,   2     },     // tracking correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 4,    10,   2,                   0,    0,   2     },     // tracking correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 3,    11,   2,                   0,    1,   2     },     // tracking correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 4,    11,   2,                   0,    1,   2     },     // tracking correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 3,    12,   2,                   0,    2,   2     },     // tracking correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 4,    12,   2,                   0,    2,   2     },     // tracking correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 3,    13,   2,                   0,    3,   2     },     // tracking correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 4,    13,   2,                   0,    3,   2     },     // tracking correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 3,    14,   2,                   0,    4,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 4,    14,   2,                   0,    4,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 3,    15,   2,                   0,    5,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 4,    15,   2,                   0,    5,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 3,    16,   2,                   0,    6,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 4,    16,   2,                   0,    6,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 3,    17,   2,                   0,    7,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 4,    17,   2,                   0,    7,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 3,    18,   2,                   0,    8,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 4,    18,   2,                   0,    8,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 3,    19,   2,                   0,    9,   2     },     // tracking correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 4,    19,   2,                   0,    9,   2     },     // tracking correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0


	     // PID corrections, pt dependence
             { 5,    10,    0,                   1,    10,    0     },     // PID correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 6,    10,    0,                   2,    10,    0     },     // PID correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 5,    11,    0,                   1,    11,    0     },     // PID correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    11,    0,                   2,    11,    0     },     // PID correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    12,    0,                   1,    12,    0     },     // PID correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    12,    0,                   2,    12,    0     },     // PID correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    13,    0,                   1,    13,    0     },     // PID correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    13,    0,                   2,    13,    0     },     // PID correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    14,    0,                   1,    14,    0     },     // PID correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 6,    14,    0,                   2,    14,    0     },     // PID correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 5,    15,    0,                   1,    15,    0     },     // PID correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 6,    15,    0,                   2,    15,    0     },     // PID correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 5,    16,    0,                   1,    16,    0     },     // PID correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 6,    16,    0,                   2,    16,    0     },     // PID correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 5,    17,    0,                   1,    17,    0     },     // PID correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 6,    17,    0,                   2,    17,    0     },     // PID correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 5,    18,    0,                   1,    18,    0     },     // PID correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 6,    18,    0,                   2,    18,    0     },     // PID correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 5,    19,    0,                   1,    19,    0     },     // PID correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 6,    19,    0,                   2,    19,    0     },     // PID correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0

	     // PID corrections, cosThetaCS dependence
	     { 5,    10,    1,                   1,    10,    1     },     // PID correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 6,    10,    1,                   2,    10,    1     },     // PID correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 5,    11,    1,                   1,    11,    1     },     // PID correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    11,    1,                   2,    11,    1     },     // PID correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    12,    1,                   1,    12,    1     },     // PID correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    12,    1,                   2,    12,    1     },     // PID correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    13,    1,                   1,    13,    1     },     // PID correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    13,    1,                   2,    13,    1     },     // PID correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    14,    1,                   1,    14,    1     },     // PID correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 6,    14,    1,                   2,    14,    1     },     // PID correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 5,    15,    1,                   1,    15,    1     },     // PID correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 6,    15,    1,                   2,    15,    1     },     // PID correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 5,    16,    1,                   1,    16,    1     },     // PID correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 6,    16,    1,                   2,    16,    1     },     // PID correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 5,    17,    1,                   1,    17,    1     },     // PID correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 6,    17,    1,                   2,    17,    1     },     // PID correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 5,    18,    1,                   1,    18,    1     },     // PID correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 6,    18,    1,                   2,    18,    1     },     // PID correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 5,    19,    1,                   1,    19,    1     },     // PID correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 6,    19,    1,                   2,    19,    1     },     // PID correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0


	     // PID corrections, cosThetaHE dependence
             { 5,    10,    2,                   1,    10,    2     },     // PID correction, SPD any   -0.9<y<+0.9  0.0<pt<10
             { 6,    10,    2,                   2,    10,    2     },     // PID correction, SPD first -0.9<y<+0.9  0.0<pt<10
	     { 5,    11,    2,                   1,    11,    2     },     // PID correction, SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    11,    2,                   2,    11,    2     },     // PID correction, SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    12,    2,                   1,    12,    2     },     // PID correction, SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    12,    2,                   2,    12,    2     },     // PID correction, SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    13,    2,                   1,    13,    2     },     // PID correction, SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    13,    2,                   2,    13,    2     },     // PID correction, SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    14,    2,                   1,    14,    2     },     // PID correction, SPD any,  -0.9<y<+0.9  0.0<pt<1.0
             { 6,    14,    2,                   2,    14,    2     },     // PID correction, SPD first,-0.9<y<+0.9  0.0<pt<1.0
             { 5,    15,    2,                   1,    15,    2     },     // PID correction, SPD any,  -0.9<y<+0.9  1.0<pt<2.0
             { 6,    15,    2,                   2,    15,    2     },     // PID correction, SPD first,-0.9<y<+0.9  1.0<pt<2.0
	     { 5,    16,    2,                   1,    16,    2     },     // PID correction, SPD any,  -0.9<y<+0.9  2.0<pt<3.0
             { 6,    16,    2,                   2,    16,    2     },     // PID correction, SPD first,-0.9<y<+0.9  2.0<pt<3.0
	     { 5,    17,    2,                   1,    17,    2     },     // PID correction, SPD any,  -0.9<y<+0.9  3.0<pt<5.0
             { 6,    17,    2,                   2,    17,    2     },     // PID correction, SPD first,-0.9<y<+0.9  3.0<pt<5.0
	     { 5,    18,    2,                   1,    18,    2     },     // PID correction, SPD any,  -0.9<y<+0.9  5.0<pt<7.0
             { 6,    18,    2,                   2,    18,    2     },     // PID correction, SPD first,-0.9<y<+0.9  5.0<pt<7.0
	     { 5,    19,    2,                   1,    19,    2     },     // PID correction, SPD any,  -0.9<y<+0.9  7.0<pt<10.0
             { 6,    19,    2,                   2,    19,    2     },     // PID correction, SPD first,-0.9<y<+0.9  7.0<pt<10.0
             

	     // correction for the invariant mass integration range, pt dependence
	     { 5,    20,   0,                   5,   10,   0     },     // SPD any, -0.9<y<+0.9  0.0<pt<10
	     { 6,    20,   0,                   6,   10,   0     },     // SPD first, -0.9<y<+0.9  0.0<pt<10
	     { 5,    21,   0,                   5,   11,   0     },     // SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    21,   0,                   6,   11,   0     },     // SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    22,   0,                   5,   12,   0     },     // SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    22,   0,                   6,   12,   0     },     // SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    23,   0,                   5,   13,   0     },     // SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    23,   0,                   6,   13,   0     },     // SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    24,   0,                   5,   14,   0     },     // SPD any, -0.9<y<+0.9  0.0<pt<1.0
	     { 6,    24,   0,                   6,   14,   0     },     // SPD first, -0.9<y<+0.9  0.0<pt<1.0
	     { 5,    25,   0,                   5,   15,   0     },     // SPD any, -0.9<y<+0.9  1.0<pt<2.0
	     { 6,    25,   0,                   6,   15,   0     },     // SPD first, -0.9<y<+0.9  1.0<pt<2.0
	     { 5,    26,   0,                   5,   16,   0     },     // SPD any, -0.9<y<+0.9  2.0<pt<3.0
	     { 6,    26,   0,                   6,   16,   0     },     // SPD first, -0.9<y<+0.9  2.0<pt<3.0
	     { 5,    27,   0,                   5,   17,   0     },     // SPD any, -0.9<y<+0.9  3.0<pt<5.0
	     { 6,    27,   0,                   6,   17,   0     },     // SPD first, -0.9<y<+0.9  3.0<pt<5.0
	     { 5,    28,   0,                   5,   18,   0     },     // SPD any, -0.9<y<+0.9  5.0<pt<7
	     { 6,    28,   0,                   6,   18,   0     },     // SPD first, -0.9<y<+0.9  5.0<pt<7
	     { 5,    29,   0,                   5,   19,   0     },     // SPD any, -0.9<y<+0.9  7<pt<10
	     { 6,    29,   0,                   6,   19,   0     },      // SPD first, -0.9<y<+0.9  7.0<pt<10


	     // correction for the invariant mass integration range, cosThetaCS dependence
	     { 5,    20,   1,                   5,   10,   1     },     // SPD any, -0.9<y<+0.9  0.0<pt<10
	     { 6,    20,   1,                   6,   10,   1     },     // SPD first, -0.9<y<+0.9  0.0<pt<10
	     { 5,    21,   1,                   5,   11,   1     },     // SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    21,   1,                   6,   11,   1     },     // SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    22,   1,                   5,   12,   1     },     // SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    22,   1,                   6,   12,   1     },     // SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    23,   1,                   5,   13,   1     },     // SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    23,   1,                   6,   13,   1     },     // SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    24,   1,                   5,   14,   1     },     // SPD any, -0.9<y<+0.9  0.0<pt<1.0
	     { 6,    24,   1,                   6,   14,   1     },     // SPD first, -0.9<y<+0.9  0.0<pt<1.0
	     { 5,    25,   1,                   5,   15,   1     },     // SPD any, -0.9<y<+0.9  1.0<pt<2.0
	     { 6,    25,   1,                   6,   15,   1     },     // SPD first, -0.9<y<+0.9  1.0<pt<2.0
	     { 5,    26,   1,                   5,   16,   1     },     // SPD any, -0.9<y<+0.9  2.0<pt<3.0
	     { 6,    26,   1,                   6,   16,   1     },     // SPD first, -0.9<y<+0.9  2.0<pt<3.0
	     { 5,    27,   1,                   5,   17,   1     },     // SPD any, -0.9<y<+0.9  3.0<pt<5.0
	     { 6,    27,   1,                   6,   17,   1     },     // SPD first, -0.9<y<+0.9  3.0<pt<5.0
	     { 5,    28,   1,                   5,   18,   1     },     // SPD any, -0.9<y<+0.9  5.0<pt<7
	     { 6,    28,   1,                   6,   18,   1     },     // SPD first, -0.9<y<+0.9  5.0<pt<7
	     { 5,    29,   1,                   5,   19,   1     },     // SPD any, -0.9<y<+0.9  7<pt<10
	     { 6,    29,   1,                   6,   19,   1     },      // SPD first, -0.9<y<+0.9  7.0<pt<10
	     

	     // correction for the invariant mass integration range, cosThetaHE dependence
	     { 5,    20,   2,                   5,   10,   2     },     // SPD any, -0.9<y<+0.9  0.0<pt<10
	     { 6,    20,   2,                   6,   10,   2     },     // SPD first, -0.9<y<+0.9  0.0<pt<10
	     { 5,    21,   2,                   5,   11,   2     },     // SPD any   -0.3<y<+0.3  0.0<pt<10
             { 6,    21,   2,                   6,   11,   2     },     // SPD first -0.3<y<+0.3  0.0<pt<10
	     { 5,    22,   2,                   5,   12,   2     },     // SPD any    0.3<y<+0.9  0.0<pt<10
             { 6,    22,   2,                   6,   12,   2     },     // SPD first  0.3<y<+0.9  0.0<pt<10
	     { 5,    23,   2,                   5,   13,   2     },     // SPD any   -0.9<y<-0.3  0.0<pt<10
             { 6,    23,   2,                   6,   13,   2     },     // SPD first -0.9<y<-0.3  0.0<pt<10
	     { 5,    24,   2,                   5,   14,   2     },     // SPD any, -0.9<y<+0.9  0.0<pt<1.0
	     { 6,    24,   2,                   6,   14,   2     },     // SPD first, -0.9<y<+0.9  0.0<pt<1.0
	     { 5,    25,   2,                   5,   15,   2     },     // SPD any, -0.9<y<+0.9  1.0<pt<2.0
	     { 6,    25,   2,                   6,   15,   2     },     // SPD first, -0.9<y<+0.9  1.0<pt<2.0
	     { 5,    26,   2,                   5,   16,   2     },     // SPD any, -0.9<y<+0.9  2.0<pt<3.0
	     { 6,    26,   2,                   6,   16,   2     },     // SPD first, -0.9<y<+0.9  2.0<pt<3.0
	     { 5,    27,   2,                   5,   17,   2     },     // SPD any, -0.9<y<+0.9  3.0<pt<5.0
	     { 6,    27,   2,                   6,   17,   2     },     // SPD first, -0.9<y<+0.9  3.0<pt<5.0
	     { 5,    28,   2,                   5,   18,   2     },     // SPD any, -0.9<y<+0.9  5.0<pt<7
	     { 6,    28,   2,                   6,   18,   2     },     // SPD first, -0.9<y<+0.9  5.0<pt<7
	     { 5,    29,   2,                   5,   19,   2     },     // SPD any, -0.9<y<+0.9  7<pt<10
	     { 6,    29,   2,                   6,   19,   2     }      // SPD first, -0.9<y<+0.9  7.0<pt<10
	     

};
// custom names and titles for efficiency histograms
const Char_t* gkEffNames[gkNeffs][3] = {
  // full corrections, pt dependence
  {"fullCorrectionSPDany_pt",      "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       "any,-0.9,0.9,0.0,10.0"},
  {"fullCorrectionSPDfirst_pt",    "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     "first,-0.9,0.9,0.0,10.0"},
  {"fullCorrectionY1SPDany_pt",    "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       "any,-0.3,0.3,0.0,10.0"},
  {"fullCorrectionY1SPDfirst_pt",  "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     "first,-0.3,0.3,0.0,10.0"},
  {"fullCorrectionY2SPDany_pt",    "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     "any,0.3,0.9,0.0,10.0"},
  {"fullCorrectionY2SPDfirst_pt",  "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   "first,0.3,0.9,0.0,10.0"},
  {"fullCorrectionY3SPDany_pt",    "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   "any,-0.9,-0.3,0.0,10.0"},
  {"fullCorrectionY3SPDfirst_pt",  "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", "first,-0.9,-0.3,0.0,10.0"},
  {"fullCorrectionPt1SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              "any,-0.9,0.9,0.0,1.0"},
  {"fullCorrectionPt1SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            "first,-0.9,0.9,0.0,1.0"},
  {"fullCorrectionPt2SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            "any,-0.9,0.9,1.0,2.0"},
  {"fullCorrectionPt2SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          "first,-0.9,0.9,1.0,2.0"},
  {"fullCorrectionPt3SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            "any,-0.9,0.9,2.0,3.0"},
  {"fullCorrectionPt3SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          "first,-0.9,0.9,2.0,3.0"},
  {"fullCorrectionPt4SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            "any,-0.9,0.9,3.0,5.0"},
  {"fullCorrectionPt4SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          "first,-0.9,0.9,3.0,5.0"},
  {"fullCorrectionPt5SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            "any,-0.9,0.9,5.0,7.0"},
  {"fullCorrectionPt5SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          "first,-0.9,0.9,5.0,7.0"},
  {"fullCorrectionPt6SPDany_pt",   "Eff. vs. pt, (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           "any,-0.9,0.9,7.0,10.0"},
  {"fullCorrectionPt6SPDfirst_pt", "Eff. vs. pt, (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         "first,-0.9,0.9,7.0,10.0"},
  // full corrections, cosThetaCS dependence
  {"fullCorrectionSPDany_ThetaCS",      "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"fullCorrectionSPDfirst_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"fullCorrectionY1SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"fullCorrectionY1SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"fullCorrectionY2SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"fullCorrectionY2SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"fullCorrectionY3SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"fullCorrectionY3SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"fullCorrectionPt1SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"fullCorrectionPt1SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"fullCorrectionPt2SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"fullCorrectionPt2SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"fullCorrectionPt3SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"fullCorrectionPt3SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"fullCorrectionPt4SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"fullCorrectionPt4SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"fullCorrectionPt5SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"fullCorrectionPt5SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"fullCorrectionPt6SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"fullCorrectionPt6SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // full corrections, cosThetaHE dependence
  {"fullCorrectionSPDany_ThetaHE",      "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"fullCorrectionSPDfirst_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"fullCorrectionY1SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"fullCorrectionY1SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"fullCorrectionY2SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"fullCorrectionY2SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"fullCorrectionY3SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"fullCorrectionY3SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"fullCorrectionPt1SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"fullCorrectionPt1SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"fullCorrectionPt2SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"fullCorrectionPt2SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"fullCorrectionPt3SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"fullCorrectionPt3SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"fullCorrectionPt4SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"fullCorrectionPt4SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"fullCorrectionPt5SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"fullCorrectionPt5SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"fullCorrectionPt6SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (full cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"fullCorrectionPt6SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (full cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // tracking corrections, pt dependence
  {"trackingCorrectionSPDany_pt",      "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionSPDfirst_pt",    "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY1SPDany_pt",    "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionY1SPDfirst_pt",  "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2SPDany_pt",    "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2SPDfirst_pt",  "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3SPDany_pt",    "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3SPDfirst_pt",  "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"trackingCorrectionPt1SPDany_pt",   "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"trackingCorrectionPt1SPDfirst_pt", "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"trackingCorrectionPt2SPDany_pt",   "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"trackingCorrectionPt2SPDfirst_pt", "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"trackingCorrectionPt3SPDany_pt",   "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"trackingCorrectionPt3SPDfirst_pt", "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"trackingCorrectionPt4SPDany_pt",   "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"trackingCorrectionPt4SPDfirst_pt", "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"trackingCorrectionPt5SPDany_pt",   "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"trackingCorrectionPt5SPDfirst_pt", "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"trackingCorrectionPt6SPDany_pt",   "Eff. vs. pt, (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"trackingCorrectionPt6SPDfirst_pt", "Eff. vs. pt, (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // tracking corrections, cosThetaCS dependence
  {"trackingCorrectionSPDany_ThetaCS",      "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionSPDfirst_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY1SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionY1SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"trackingCorrectionPt1SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"trackingCorrectionPt1SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"trackingCorrectionPt2SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"trackingCorrectionPt2SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"trackingCorrectionPt3SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"trackingCorrectionPt3SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"trackingCorrectionPt4SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"trackingCorrectionPt4SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"trackingCorrectionPt5SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"trackingCorrectionPt5SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"trackingCorrectionPt6SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"trackingCorrectionPt6SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // tracking corrections, cosThetaHE dependence
  {"trackingCorrectionSPDany_ThetaHE",      "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionSPDfirst_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY1SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionY1SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"trackingCorrectionPt1SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"trackingCorrectionPt1SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"trackingCorrectionPt2SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"trackingCorrectionPt2SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"trackingCorrectionPt3SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"trackingCorrectionPt3SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"trackingCorrectionPt4SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"trackingCorrectionPt4SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"trackingCorrectionPt5SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"trackingCorrectionPt5SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"trackingCorrectionPt6SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"trackingCorrectionPt6SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // tracking corrections with V0 cuts, pt dependence
  {"trackingCorrectionV0SPDany_pt",      "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionV0SPDfirst_pt",    "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY1V0SPDany_pt",    "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionY1V0SPDfirst_pt",  "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2V0SPDany_pt",    "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2V0SPDfirst_pt",  "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3V0SPDany_pt",    "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3V0SPDfirst_pt",  "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"trackingCorrectionV0Pt1SPDany_pt",   "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"trackingCorrectionV0Pt1SPDfirst_pt", "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"trackingCorrectionV0Pt2SPDany_pt",   "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"trackingCorrectionV0Pt2SPDfirst_pt", "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"trackingCorrectionV0Pt3SPDany_pt",   "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"trackingCorrectionV0Pt3SPDfirst_pt", "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"trackingCorrectionV0Pt4SPDany_pt",   "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"trackingCorrectionV0Pt4SPDfirst_pt", "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"trackingCorrectionV0Pt5SPDany_pt",   "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"trackingCorrectionV0Pt5SPDfirst_pt", "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"trackingCorrectionV0Pt6SPDany_pt",   "Eff. vs. pt, (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"trackingCorrectionV0Pt6SPDfirst_pt", "Eff. vs. pt, (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // tracking corrections with V0 cuts, cosThetaCS dependence
  {"trackingCorrectionV0SPDany_ThetaCS",      "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionV0SPDfirst_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY1V0SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionY1V0SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2V0SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2V0SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3V0SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3V0SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"trackingCorrectionV0Pt1SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"trackingCorrectionV0Pt1SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"trackingCorrectionV0Pt2SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"trackingCorrectionV0Pt2SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"trackingCorrectionV0Pt3SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"trackingCorrectionV0Pt3SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"trackingCorrectionV0Pt4SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"trackingCorrectionV0Pt4SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"trackingCorrectionV0Pt5SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"trackingCorrectionV0Pt5SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"trackingCorrectionV0Pt6SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"trackingCorrectionV0Pt6SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // tracking corrections with V0 cuts, cosThetaHE dependence
  {"trackingCorrectionV0SPDany_ThetaHE",      "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionV0SPDfirst_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY1V0SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"trackingCorrectionY1V0SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2V0SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"trackingCorrectionY2V0SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3V0SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"trackingCorrectionY3V0SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"trackingCorrectionV0Pt1SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"trackingCorrectionV0Pt1SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"trackingCorrectionV0Pt2SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"trackingCorrectionV0Pt2SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"trackingCorrectionV0Pt3SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"trackingCorrectionV0Pt3SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"trackingCorrectionV0Pt4SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"trackingCorrectionV0Pt4SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"trackingCorrectionV0Pt5SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"trackingCorrectionV0Pt5SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"trackingCorrectionV0Pt6SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"trackingCorrectionV0Pt6SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (tracking +V0 cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},

  // pid corrections, pt dependence
  {"pidCorrectionSPDany_pt",      "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"pidCorrectionSPDfirst_pt",    "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY1SPDany_pt",    "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"pidCorrectionY1SPDfirst_pt",  "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY2SPDany_pt",    "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY2SPDfirst_pt",  "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"pidCorrectionY3SPDany_pt",    "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"pidCorrectionY3SPDfirst_pt",  "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"pidCorrectionPt1SPDany_pt",   "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"pidCorrectionPt1SPDfirst_pt", "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"pidCorrectionPt2SPDany_pt",   "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"pidCorrectionPt2SPDfirst_pt", "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"pidCorrectionPt3SPDany_pt",   "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"pidCorrectionPt3SPDfirst_pt", "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"pidCorrectionPt4SPDany_pt",   "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"pidCorrectionPt4SPDfirst_pt", "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"pidCorrectionPt5SPDany_pt",   "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"pidCorrectionPt5SPDfirst_pt", "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"pidCorrectionPt6SPDany_pt",   "Eff. vs. pt, (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"pidCorrectionPt6SPDfirst_pt", "Eff. vs. pt, (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // pid corrections, cosThetaCS dependence
  {"pidCorrectionSPDany_ThetaCS",      "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"pidCorrectionSPDfirst_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY1SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"pidCorrectionY1SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY2SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY2SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"pidCorrectionY3SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"pidCorrectionY3SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"pidCorrectionPt1SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"pidCorrectionPt1SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"pidCorrectionPt2SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"pidCorrectionPt2SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"pidCorrectionPt3SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"pidCorrectionPt3SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"pidCorrectionPt4SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"pidCorrectionPt4SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"pidCorrectionPt5SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"pidCorrectionPt5SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"pidCorrectionPt6SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"pidCorrectionPt6SPDfirst_ThataCS", "Eff. vs. cos(#theta_{CS}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // pid corrections, cosThetaHE dependence
  {"pidCorrectionSPDany_ThetaHE",      "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",       ""},
  {"pidCorrectionSPDfirst_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY1SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",       ""},
  {"pidCorrectionY1SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.3 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY2SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",     ""},
  {"pidCorrectionY2SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in 0.3<y<0.9 & 0<pt<10.0 GeV/c)",   ""},
  {"pidCorrectionY3SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)",   ""},
  {"pidCorrectionY3SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c)", ""},
  {"pidCorrectionPt1SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",              ""},
  {"pidCorrectionPt1SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 0<pt<1.0)",            ""},
  {"pidCorrectionPt2SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",            ""},
  {"pidCorrectionPt2SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 1.0<pt<2.0)",          ""},
  {"pidCorrectionPt3SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",            ""},
  {"pidCorrectionPt3SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 2.0<pt<3.0)",          ""},
  {"pidCorrectionPt4SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",            ""},
  {"pidCorrectionPt4SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 3.0<pt<5.0)",          ""},
  {"pidCorrectionPt5SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",            ""},
  {"pidCorrectionPt5SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 5.0<pt<7.0)",          ""},
  {"pidCorrectionPt6SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD any, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",           ""},
  {"pidCorrectionPt6SPDfirst_ThataHE", "Eff. vs. cos(#theta_{HE}), (pid cuts, SPD first, TPC PID)/(J/#Psi in |y|<0.9 & 7.0<pt<10.0)",         ""},
  // correction for the invariant mass integration range, pt dependence
  {"signalFractionSPDany_pt",      "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 0<pt<10.0 GeV/c",      "any,-0.9,0.9,0.0,10.0"},
  {"signalFractionSPDfirst_pt",    "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 0<pt<10.0 GeV/c",    "first,-0.9,0.9,0.0,10.0"},
  {"signalFractionY1SPDany_pt",    "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.3<y<+0.3 & 0<pt<10.0 GeV/c",      "any,-0.3,0.3,0.0,10.0"},
  {"signalFractionY1SPDfirst_pt",  "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.3<y<+0.3 & 0<pt<10.0 GeV/c",    "first,-0.3,0.3,0.0,10.0"},
  {"signalFractionY2SPDany_pt",    "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in 0.3<y<+0.9 & 0<pt<10.0 GeV/c",       "any,0.3,0.9,0.0,10.0"},
  {"signalFractionY2SPDfirst_pt",  "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in 0.3<y<+0.9 & 0<pt<10.0 GeV/c",     "first,0.3,0.9,0.0,10.0"},
  {"signalFractionY3SPDany_pt",    "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c",      "any,-0.9,-0.3,0.0,10.0"},
  {"signalFractionY3SPDfirst_pt",  "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c",    "first,-0.9,-0.3,0.0,10.0"},
  {"signalFractionPt1SPDany_pt",   "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 0<pt<1.0 GeV/c",       "any,-0.9,0.9,0.0,1.0"},
  {"signalFractionPt1SPDfirst_pt", "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 0<pt<1.0 GeV/c",     "first,-0.9,0.9,0.0,1.0"},
  {"signalFractionPt2SPDany_pt",   "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 1.0<pt<2.0 GeV/c",     "any,-0.9,0.9,1.0,2.0"},
  {"signalFractionPt2SPDfirst_pt", "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 1.0<pt<2.0 GeV/c",   "first,-0.9,0.9,1.0,2.0"},
  {"signalFractionPt3SPDany_pt",   "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 2.0<pt<3.0 GeV/c",     "any,-0.9,0.9,2.0,3.0"},
  {"signalFractionPt3SPDfirst_pt", "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 2.0<pt<3.0 GeV/c",   "first,-0.9,0.9,2.0,3.0"},
  {"signalFractionPt4SPDany_pt",   "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 3.0<pt<5.0 GeV/c",     "any,-0.9,0.9,3.0,5.0"},
  {"signalFractionPt4SPDfirst_pt", "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 3.0<pt<5.0 GeV/c",   "first,-0.9,0.9,3.0,5.0"},
  {"signalFractionPt5SPDany_pt",   "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 5.0<pt<7.0 GeV/c",    "any,-0.9,0.9,5.0,7.0"},
  {"signalFractionPt5SPDfirst_pt", "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 5.0<pt<7.0 GeV/c",  "first,-0.9,0.9,5.0,7.0"},
  {"signalFractionPt6SPDany_pt",   "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 7.0<pt<10.0 GeV/c",    "any,-0.9,0.9,7.0,10.0"},
  {"signalFractionPt6SPDfirst_pt", "Eff. vs. pt, (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 7.0<pt<10.0 GeV/c",  "first,-0.9,0.9,7.0,10.0"},
  // correction for the invariant mass integration range, cosThetaCS dependence
  {"signalFractionSPDany_ThetaCS",      "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 0<pt<10.0 GeV/c",      ""},
  {"signalFractionSPDfirst_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 0<pt<10.0 GeV/c",    ""},
  {"signalFractionY1SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.3<y<+0.3 & 0<pt<10.0 GeV/c",      ""},
  {"signalFractionY1SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.3<y<+0.3 & 0<pt<10.0 GeV/c",    ""},
  {"signalFractionY2SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in 0.3<y<+0.9 & 0<pt<10.0 GeV/c",       ""},
  {"signalFractionY2SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in 0.3<y<+0.9 & 0<pt<10.0 GeV/c",     ""},
  {"signalFractionY3SPDany_ThetaCS",    "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c",      ""},
  {"signalFractionY3SPDfirst_ThetaCS",  "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c",    ""},
  {"signalFractionPt1SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 0<pt<1.0 GeV/c",       ""},
  {"signalFractionPt1SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 0<pt<1.0 GeV/c",     ""},
  {"signalFractionPt2SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 1.0<pt<2.0 GeV/c",     ""},
  {"signalFractionPt2SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 1.0<pt<2.0 GeV/c",   ""},
  {"signalFractionPt3SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 2.0<pt<3.0 GeV/c",     ""},
  {"signalFractionPt3SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 2.0<pt<3.0 GeV/c",   ""},
  {"signalFractionPt4SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 3.0<pt<5.0 GeV/c",     ""},
  {"signalFractionPt4SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 3.0<pt<5.0 GeV/c",   ""},
  {"signalFractionPt5SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 5.0<pt<7.0 GeV/c",    ""},
  {"signalFractionPt5SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 5.0<pt<7.0 GeV/c",  ""},
  {"signalFractionPt6SPDany_ThetaCS",   "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 7.0<pt<10.0 GeV/c",    ""},
  {"signalFractionPt6SPDfirst_ThetaCS", "Eff. vs. cos(#theta_{CS}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 7.0<pt<10.0 GeV/c",  ""},  

  // correction for the invariant mass integration range, cosThetaHE dependence
  {"signalFractionSPDany_ThetaHE",      "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 0<pt<10.0 GeV/c",      ""},
  {"signalFractionSPDfirst_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 0<pt<10.0 GeV/c",    ""},
  {"signalFractionY1SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.3<y<+0.3 & 0<pt<10.0 GeV/c",      ""},
  {"signalFractionY1SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.3<y<+0.3 & 0<pt<10.0 GeV/c",    ""},
  {"signalFractionY2SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in 0.3<y<+0.9 & 0<pt<10.0 GeV/c",       ""},
  {"signalFractionY2SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in 0.3<y<+0.9 & 0<pt<10.0 GeV/c",     ""},
  {"signalFractionY3SPDany_ThetaHE",    "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c",      ""},
  {"signalFractionY3SPDfirst_ThetaHE",  "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<-0.3 & 0<pt<10.0 GeV/c",    ""},
  {"signalFractionPt1SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 0<pt<1.0 GeV/c",       ""},
  {"signalFractionPt1SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 0<pt<1.0 GeV/c",     ""},
  {"signalFractionPt2SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 1.0<pt<2.0 GeV/c",     ""},
  {"signalFractionPt2SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 1.0<pt<2.0 GeV/c",   ""},
  {"signalFractionPt3SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 2.0<pt<3.0 GeV/c",     ""},
  {"signalFractionPt3SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 2.0<pt<3.0 GeV/c",   ""},
  {"signalFractionPt4SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 3.0<pt<5.0 GeV/c",     ""},
  {"signalFractionPt4SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 3.0<pt<5.0 GeV/c",   ""},
  {"signalFractionPt5SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 5.0<pt<7.0 GeV/c",    ""},
  {"signalFractionPt5SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 5.0<pt<7.0 GeV/c",  ""},
  {"signalFractionPt6SPDany_ThetaHE",   "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD any, J/#Psi in -0.9<y<+0.9 & 7.0<pt<10.0 GeV/c",    ""},
  {"signalFractionPt6SPDfirst_ThetaHE", "Eff. vs. cos(#theta_{HE}), (2.92<M_{inv}<3.16) / (0.0<M_{inv}<5.0), SPD first, J/#Psi in -0.9<y<+0.9 & 7.0<pt<10.0 GeV/c",  ""}

};


// Function prototypes ---------------------------------------------------
Double_t* GetBinning(AliCFContainer* cont, Int_t variable, Int_t& nBins);
void FillHistograms(TObjArray* histosArray, AliCFContainer* cont, Int_t currentRangeStep, Bool_t firstTime);
void GetBinLimits(AliCFContainer* cont);
void DefineHistograms(TObjArray* objArray, Int_t iCutSet);
void AddHistogram(TObjArray* objArray, Int_t ndim, 
		  const Char_t* name, const Char_t* title, 
		  Int_t nbinsx, Double_t* binsx, const Char_t* xLabel = "",
		  Int_t nbinsy=0, Double_t* binsy=0, const Char_t* yLabel = "",
		  Int_t nbinsz=0, Double_t* binsz=0, const Char_t* zLabel = "");
void ProjectManyRuns(const Char_t* runList, const Char_t* pattern, Int_t howMany=1, Int_t offset = 0);
void ProjectAll(const Char_t* inputList, const Char_t* outfilename="HistosFromCFs.root", 
		Int_t howMany=1, Int_t offset=0);
void ExtractEfficienciesMany(const Char_t* runList, const Char_t* pattern, const Char_t* outAscii, Int_t howMany=1, Int_t offset=0);
void ExtractEfficiencies(const Char_t* inputFile, const Char_t* outfilename="Efficiencies.root", const Char_t* numbersFile="");
TH1* DivideHists(TH1* nominator, TH1* denominator, Int_t dimension);
//-------------------------------------------------------------------------


//_______________________________________________________________________________________
void ProjectManyRuns(const Char_t* runList,
		     const Char_t* pattern, 
		     Int_t howMany, Int_t offset) {
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

    ProjectAll(Form("%s/%s/listCF.txt", pattern, readString), 
	       Form("%s/%s/Projections.root",pattern, readString), 
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
    cf->SetRangeUser("PairType", 1.1, 1.9);
    // Pair rapidity cut
    cf->SetRangeUser("Y", -0.899, 0.899);
    cf->SetRangeUser("Pt", 0.001, 9.999);
    Int_t currentCutSet = 0;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi -0.3<y<0.3 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.299, 0.299);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 0.3<y<0.9 ----------------------------------------------------
    cf->SetRangeUser("Y", 0.301, 0.899);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi -0.9<y<-0.3 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.899, -0.301);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 0<pt<1.0 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.899, 0.899);
    cf->SetRangeUser("Pt", 0.001, 0.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 1.0<pt<2.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 1.001, 1.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 2.0<pt<3.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 2.001, 2.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 3.0<pt<5.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 3.001, 4.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 5.0<pt<7.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 5.001, 6.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 7.0<pt<10.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 7.001, 9.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // Leg pseudo-rapidity cut ---------------------------------------
    cf->SetRangeUser("Pt", 0.001, 9.999);
    cf->SetRangeUser("Y", -0.899, 0.899);
    cf->SetRangeUser("Leg1_Eta", -0.899, 0.899);
    cf->SetRangeUser("Leg2_Eta", -0.899, 0.899);
    cf->SetRangeUser("Leg1_Pt", 1.001, 10.0);
    cf->SetRangeUser("Leg2_Pt", 1.001, 10.0);
    cf->SetRangeUser("Leg1_NclsTPC", 70.1, 160.0);
    cf->SetRangeUser("Leg2_NclsTPC", 70.1, 160.0);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi -0.3<y<0.3 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.299, 0.299);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 0.3<y<0.9 ----------------------------------------------------
    cf->SetRangeUser("Y", 0.301, 0.899);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi -0.9<y<-0.3 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.899, -0.301);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 0<pt<1.0 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.899, 0.899);
    cf->SetRangeUser("Pt", 0.001, 0.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 1.0<pt<2.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 1.001, 1.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 2.0<pt<3.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 2.001, 2.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 3.0<pt<5.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 3.001, 4.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 5.0<pt<7.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 5.001, 6.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 7.0<pt<10.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 7.001, 9.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    
    // j/psi 2.92<M<3.16
    cf->SetRangeUser("Pt", 0.001, 9.999);
    cf->SetRangeUser("Y", -0.899, 0.899);
    cf->SetRangeUser("M", 2.9201, 3.1599);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi -0.3<y<0.3 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.299, 0.299);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 0.3<y<0.9 ----------------------------------------------------
    cf->SetRangeUser("Y", 0.301, 0.899);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi -0.9<y<-0.3 ----------------------------------------------------
    cf->SetRangeUser("Y", -0.899, -0.301);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 0<pt<1.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 0.001, 0.999);
    cf->SetRangeUser("Y", -0.899, 0.899);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 1.0<pt<2.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 1.001, 1.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 2.0<pt<3.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 2.001, 2.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 3.0<pt<5.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 3.001, 4.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 5.0<pt<7.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 5.001, 6.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);

    // j/psi 7.0<pt<10.0 ----------------------------------------------------
    cf->SetRangeUser("Pt", 7.001, 9.999);
    currentCutSet++;
    FillHistograms(histoArray, cont, currentCutSet, firstTime);
    

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
void ExtractEfficienciesMany(const Char_t* runList, const Char_t* pattern, const Char_t* outAscii, Int_t howMany, Int_t offset) {
  //
  //
  //

  // loop over all runs -----------------------
  ifstream input; input.open(runList);
  Int_t runCounter = 0;

  TGraphErrors* trends[gkNeffs];
  Double_t weightedEffs[gkNeffs];
  Double_t weightedErrs[gkNeffs];
  Int_t nPoints[gkNeffs];
  Double_t nTotalEvents[gkNeffs];
  for(Int_t iTrend=0; iTrend<gkNeffs; iTrend++) {
    trends[iTrend] = new TGraphErrors();
    trends[iTrend]->SetName(gkEffNames[iTrend][0]);
    trends[iTrend]->SetTitle(gkEffNames[iTrend][1]);
    weightedEffs[iTrend] = 0.0; weightedErrs[iTrend] = 0.0;
    nPoints[iTrend] = 0;
    nTotalEvents[iTrend] = 0;
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
    if(run>121000 && run<=126437) periodStr = "LHC10d.pass2";
    Double_t nPhysicsEvents = 0;
    normalizationFile = TFile::Open(Form("/lustre/alice/train/V006.pp/2011-03-18_2242.6024/mergedRuns/pp/7TeV/%s/%d.ana/iarsene_normalization.root", periodStr.Data(), run));
    cout << "# physics events = ";
    if(normalizationFile) {
      TObjArray *histos=(TObjArray*)normalizationFile->Get("iarsene_normalization");
      TH1I* triggers=(TH1I*)histos->FindObject("TriggersHistogram");
      nPhysicsEvents = triggers->GetBinContent(2);   // PHYSICS events
      cout << nPhysicsEvents;
      normalizationFile->Close();
    }
    else
      cout << " NOT FOUND";
    cout << endl;

    ExtractEfficiencies(Form("%s/%s/Projections.root",pattern,readString), 
			Form("%s/%s/Efficiencies.root",pattern,readString));
    file = TFile::Open(Form("%s/%s/Efficiencies.root",pattern,readString));
    if(file && !file->IsZombie()) {   
      for(Int_t iTrend=0; iTrend<gkNeffs; iTrend++) {
	object = (TNamed*)file->Get(Form("%s_value",gkEffNames[iTrend][0]));
	if(!object) continue;
	Float_t eff = (TString(object->GetTitle())).Atof();
	trends[iTrend]->SetPoint(nPoints[iTrend], run, eff);
	object = (TNamed*)file->Get(Form("%s_error",gkEffNames[iTrend][0]));
	if(!object) continue;
	Float_t err = (TString(object->GetTitle())).Atof();
	trends[iTrend]->SetPointError(nPoints[iTrend], 0.0, err);
	weightedEffs[iTrend] += nPhysicsEvents*eff;
	weightedErrs[iTrend] += nPhysicsEvents*nPhysicsEvents*err*err;
	nTotalEvents[iTrend] += nPhysicsEvents;
	//	cout << "trend " << iTrend << "; eff = " << eff << endl;
	nPoints[iTrend]+=1;
      }
      file->Close();
    }
    if(normalizationFile)
      normalizationFile->Close();
    runCounter++;
  }

  // write the efficiencies also in an ascii file
  ofstream asciiOut;
  asciiOut.open(outAscii);

  TFile *saveTrend = new TFile(Form("%s.trend.root", runList), "RECREATE");
  TNamed *weightedFactors;
  TNamed *weightedErrors;
  TNamed *nEventsObject;
  for(Int_t iTrend=0; iTrend<gkNeffs; iTrend++) {
    trends[iTrend]->Write();
    weightedEffs[iTrend] /= nTotalEvents[iTrend];
    weightedErrs[iTrend] = TMath::Sqrt(weightedErrs[gkNeffs])/nTotalEvents[iTrend];
    weightedFactors = new TNamed(Form("%s_weighted", gkEffNames[iTrend][0]),
				 Form("%f", weightedEffs[iTrend]));
    weightedErrors = new TNamed(Form("%s_weightedErr", gkEffNames[iTrend][0]),
				 Form("%f", weightedErrs[iTrend]));
    weightedFactors->Write();
    weightedErrors->Write();
    nEventsObject = new TNamed(Form("TotalEvents_%s",gkEffNames[iTrend][0]), Form("%f",nTotalEvents[iTrend]));
    nEventsObject->Write();
    // write a table into an ascii file ---------------------------------
    TString effNameStr(gkEffNames[iTrend][2]);
    // The array contains y,pt rapidity intervals for this efficiency, and/or other variables
    // Always, the first element of the array should be "any" or "first"
    TObjArray* array = effNameStr.Tokenize(",");
    if(array->GetEntries()<=1) continue;
    asciiOut << iTrend << "\t";
    if(((TObjString*)array->At(0))->GetString()=="any")
      asciiOut << "1";
    else if(((TObjString*)array->At(0))->GetString()=="first")
      asciiOut << "2";
    for(Int_t iStr=1; iStr<array->GetEntries(); iStr++)
      asciiOut << "\t" << ((TObjString*)array->At(iStr))->GetString().Data();
    asciiOut << "\t " << weightedEffs[iTrend] << endl;
  }
  // At the end of the ascii file write the format and brief explanations
  asciiOut << endl;
  asciiOut << "# Format:  efficiencyId    SPD(1-any/2-first)     yLow       yHigh      ptLow      ptHigh      efficiency" << endl;
  asciiOut << "# Efficiency descriptions (based on efficiencyId) : " << endl;
  for(Int_t iTrend=0; iTrend<gkNeffs; iTrend++) {
    TString effNameStr(gkEffNames[iTrend][2]);
    TObjArray* array = effNameStr.Tokenize(",");
    if(array->GetEntries()<=1) continue;
    asciiOut << "# " << iTrend << " - " << gkEffNames[iTrend][1] << endl;
  }

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
    //    cout << gkEffNames[iEff][0] << " (" << gkEffNames[iEff][1] << " )" << endl;
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
    // Error calculation: take into account that nominator and denominator are correlated.
    // The nominator is a subset of denominator
    Double_t error = eff;
    if(nomIntegral>0 && denomIntegral>0)
      error = eff*TMath::Sqrt(TMath::Abs(denomIntegral-nomIntegral)/nomIntegral/denomIntegral);
    //nominator->Divide(denominator);
    TH1* ratio = DivideHists(nominator, denominator, gkDims[gkEffs[iEff][2]][0]);
    //    cout << "efficiency = " << nomIntegral << " / " << denomIntegral << " = "
    //	 << eff << " +/- " << error << endl;
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
	//histo->Add(cont->Project(gkDims[iHisto][1],gkStepNumbers[iStep]));

	//cout << "Histo: " << Form("%s_%s_%s",gkStepNames[iStep][0],gkCutSetNames[currentCutSet][0],gkHistoNames[iHisto][0]) << endl;
	//cout << "bin lims x: ";
	//for(Int_t iBinx=1; iBinx<=histo->GetXaxis()->GetNbins(); iBinx++)
	//  cout << histo->GetXaxis()->GetBinLowEdge(iBinx) << "  ";
	//cout << histo->GetXaxis()->GetBinUpEdge(histo->GetXaxis()->GetNbins()) << endl;
	histo->Add(cont->Project(gkStepNumbers[iStep],gkDims[iHisto][1]));
      }
      // fill 2-dim histos
      if(gkDims[iHisto][0]==2) {
	histo = (TH2D*)histosArray->FindObject(Form("%s_%s_%s",gkStepNames[iStep][0],
						    gkCutSetNames[currentCutSet][0],
						    gkHistoNames[iHisto][0]));
	//histo->Add(cont->Project(gkDims[iHisto][1], gkDims[iHisto][2], gkStepNumbers[iStep]));
	//cout << "Histo: " << Form("%s_%s_%s",gkStepNames[iStep][0],gkCutSetNames[currentCutSet][0],gkHistoNames[iHisto][0]) << endl;
	//cout << "bin lims x: ";
	//for(Int_t iBinx=1; iBinx<=histo->GetXaxis()->GetNbins(); iBinx++)
	//  cout << histo->GetXaxis()->GetBinLowEdge(iBinx) << "  ";
	//cout << histo->GetXaxis()->GetBinUpEdge(histo->GetXaxis()->GetNbins()) << endl;
	//cout << "bin lims y: ";
	//for(Int_t iBiny=1; iBiny<=histo->GetYaxis()->GetNbins(); iBiny++)
	//  cout << histo->GetYaxis()->GetBinLowEdge(iBiny) << "  ";
	//cout << histo->GetYaxis()->GetBinUpEdge(histo->GetYaxis()->GetNbins()) << endl;
	histo->Add(cont->Project(gkStepNumbers[iStep], gkDims[iHisto][1], gkDims[iHisto][2]));
      }
      // fill 3-dim histos
      if(gkDims[iHisto][0]==3) {
	histo = (TH3D*)histosArray->FindObject(Form("%s_%s_%s",gkStepNames[iStep][0],
						    gkCutSetNames[currentCutSet][0],
						    gkHistoNames[iHisto][0]));
	//histo->Add(cont->Project(gkDims[iHisto][1], gkDims[iHisto][2], gkDims[iHisto][3], gkStepNumbers[iStep]));
	histo->Add(cont->Project(gkStepNumbers[iStep], gkDims[iHisto][1], gkDims[iHisto][2], gkDims[iHisto][3]));
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
    //cout << "bin limits = ";
    //for(Int_t iBin=0; iBin<=gNbins[iVar]; iBin++) 
    //  cout << gBinLimits[iVar][iBin] << "  ";
    //cout << endl;
  }
}

//________________________________________________________________________________________
Double_t* GetBinning(AliCFContainer* cont, Int_t variable, 
		     Int_t& nBins) {
  //
  // Get the number of bins and the bin limits for the projection of a given variable
  //
  //TH1D* tempHist = cont->Project(variable, kPureMC);
  TH1* tempHist = cont->Project(kPureMC, variable);
  nBins = tempHist->GetXaxis()->GetNbins();
  Double_t* binLimits = new Double_t[nBins+1];
  for(Int_t i=1; i<=nBins; i++)
    binLimits[i-1]=tempHist->GetXaxis()->GetBinLowEdge(i);
  binLimits[nBins] = tempHist->GetXaxis()->GetBinLowEdge(nBins) + 
    tempHist->GetXaxis()->GetBinWidth(nBins);
  return binLimits;
}

//________________________________________________________________________________________
TH1* DivideHists(TH1* nominator, TH1* denominator, Int_t dimension) {
  //
  // divide 2 histograms with error propagation
  //
  TH1* ratio;
  if(dimension==3) {
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
	  Double_t eff = countsN/countsD;
	  Double_t error = eff*TMath::Sqrt(TMath::Abs(countsD-countsN)/countsN/countsD);
	  ratio->SetBinContent(iXbin, iYbin, iZbin, eff);
	  ratio->SetBinError(iXbin, iYbin, iZbin, error);
	}
      }
    }
    return ratio;
  }

  if(dimension==2) {
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
	Double_t eff = countsN/countsD;
	Double_t error = eff*TMath::Sqrt(TMath::Abs(countsD-countsN)/countsN/countsD);
	ratio->SetBinContent(iXbin, iYbin, eff);
	ratio->SetBinError(iXbin, iYbin, error);
      }
    }
    return ratio;
  }

  if(dimension==1) {
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
      Double_t eff = countsN/countsD;
      Double_t error = eff*TMath::Sqrt(TMath::Abs(countsD-countsN)/countsN/countsD);
      ratio->SetBinContent(iXbin, eff);
      ratio->SetBinError(iXbin, error);
    }
    return ratio;
  }
    
  return 0x0;
}
