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

// c++ headers
#include <iostream>
#include <fstream>
// #include <vector>
// #include <algorithm>


// root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1.h"
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"


// aliroot headers
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliAODVertex.h"


// my headers
#include "AliAnalysisTaskMatchTriggerForward.h"


// headers for MC
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"



class AliAnalysisTaskMatchTriggerForward;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMatchTriggerForward) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskMatchTriggerForward::AliAnalysisTaskMatchTriggerForward()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fLumiPerRun(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      // fVectorGoodRunNumbers(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0),
      fEfficiencyPerRunH(0),
      fMCEfficiencyPerRunH(0),
      fEfficiencyPerRunWithTriggeringH(0),
      fSingleMuonPtDistributionH(0),
      fZNCEnergyAgainstEntriesH(0),
      fZNAEnergyAgainstEntriesH(0),
      fZNCEnergyBeforeTimingSelectionH(0),
      fZNAEnergyBeforeTimingSelectionH(0),
      fZNCEnergyCalibratedH(0),
      fZNAEnergyCalibratedH(0),
      fZNCEnergyUncalibratedH(0),
      fZNAEnergyUncalibratedH(0),
      fZNCEnergyCalibratedHigherGainH(0),
      fZNAEnergyCalibratedHigherGainH(0),
      fZNCTimeAgainstEntriesH(0),
      fZNATimeAgainstEntriesH(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskMatchTriggerForward::AliAnalysisTaskMatchTriggerForward( const char* name )
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fLumiPerRun(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      // fVectorGoodRunNumbers(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0),
      fEfficiencyPerRunH(0),
      fMCEfficiencyPerRunH(0),
      fEfficiencyPerRunWithTriggeringH(0),
      fSingleMuonPtDistributionH(0),
      fZNCEnergyAgainstEntriesH(0),
      fZNAEnergyAgainstEntriesH(0),
      fZNCEnergyBeforeTimingSelectionH(0),
      fZNAEnergyBeforeTimingSelectionH(0),
      fZNCEnergyCalibratedH(0),
      fZNAEnergyCalibratedH(0),
      fZNCEnergyUncalibratedH(0),
      fZNAEnergyUncalibratedH(0),
      fZNCEnergyCalibratedHigherGainH(0),
      fZNAEnergyCalibratedHigherGainH(0),
      fZNCTimeAgainstEntriesH(0),
      fZNATimeAgainstEntriesH(0)
{
    // FillGoodRunVector(fVectorGoodRunNumbers);
    for( Int_t iRun = 0; iRun < 60000; iRun++) {
      fCounterGeneratedLevel[iRun] = 0;
    }

    for( Int_t iRun = 0; iRun < 364; iRun++) {
      fDeadZoneEtaVsPhiPerRunH[iRun]               = 0x0;
      fDeadZoneEtaVsPhiPerRunWithTriggeringH[iRun] = 0x0;
      fZNCEnergyPerRunH[iRun]                      = 0x0;
      fZNAEnergyPerRunH[iRun]                      = 0x0;
    }


    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskMatchTriggerForward::~AliAnalysisTaskMatchTriggerForward()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
void AliAnalysisTaskMatchTriggerForward::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //

  //muon track cuts
  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonTrackCuts->SetFilterMask(    AliMuonTrackCuts::kMuEta     |
                                    AliMuonTrackCuts::kMuThetaAbs|
                                    AliMuonTrackCuts::kMuPdca    //|
                                    // AliMuonTrackCuts::kMuMatchLpt
                                    );
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->Print("mask");


  fOutputList = new TList();          // this is a list which will contain all
                                      // of your histograms at the end of the
                                      // analysis, the contents of this list
                                      // are written to the output file

  fOutputList->SetOwner(kTRUE);       // memory management: the list is owner
                                      // of all objects it contains and will
                                      // delete them if requested

  //_______________________________
  // - Adding histograms
  fNumberMuonsH = new TH1F("fNumberMuonsH", "fNumberMuonsH", 12, -0.5, 11.5);
  fOutputList->Add(fNumberMuonsH);    // don't forget to add it to the list!

  fCounterH = new TH1F("fCounterH", "fCounterH", 24, -0.5, 23.5);
  fOutputList->Add(fCounterH);

  fEtaMuonH = new TH1F("fEtaMuonH", "fEtaMuonH", 90, -2, -5);
  fOutputList->Add(fEtaMuonH);

  fRAbsMuonH = new TH1F("fRAbsMuonH", "fRAbsMuonH", 100, 0, 100);
  fOutputList->Add(fRAbsMuonH);

  fEntriesAgainstRunNumberH = new TH1F("fEntriesAgainstRunNumberH", "fEntriesAgainstRunNumberH", 10000, 290000, 300000);
  fOutputList->Add(fEntriesAgainstRunNumberH);

  /* - Trying to reproduce the histogram for the RunNumbers as they always
     - show it, properly labelled. Inspiration has come from the website:
     - https://root.cern.ch/doc/master/hlabels1_8C.html
     - Let us see if it works properly.
     -
   */
  fEntriesAgainstRunNumberProperlyH = new TH1F("fEntriesAgainstRunNumberProperlyH", "fEntriesAgainstRunNumberProperlyH", 10000, 290000, 300000);
  fEntriesAgainstRunNumberProperlyH->SetStats(0);
  fEntriesAgainstRunNumberProperlyH->SetFillColor(38);
  fEntriesAgainstRunNumberProperlyH->LabelsDeflate();
  fOutputList->Add(fEntriesAgainstRunNumberProperlyH);

  fEfficiencyPerRunH = new TH1F("fEfficiencyPerRunH", "fEfficiencyPerRunH", 3, 0, 3);
  fEfficiencyPerRunH->SetStats(0);
  fEfficiencyPerRunH->SetFillColor(38);
  fEfficiencyPerRunH->LabelsDeflate();
  fOutputList->Add(fEfficiencyPerRunH);

  fMCEfficiencyPerRunH = new TH1F("fMCEfficiencyPerRunH", "fMCEfficiencyPerRunH", 3, 0, 3);
  fMCEfficiencyPerRunH->SetStats(0);
  fMCEfficiencyPerRunH->SetFillColor(38);
  fMCEfficiencyPerRunH->LabelsDeflate();
  fOutputList->Add(fMCEfficiencyPerRunH);

  fEfficiencyPerRunWithTriggeringH = new TH1F("fEfficiencyPerRunWithTriggeringH", "fEfficiencyPerRunWithTriggeringH", 3, 0, 3);
  fEfficiencyPerRunWithTriggeringH->SetStats(0);
  fEfficiencyPerRunWithTriggeringH->SetFillColor(38);
  fEfficiencyPerRunWithTriggeringH->LabelsDeflate();
  fOutputList->Add(fEfficiencyPerRunWithTriggeringH);

  /* - Eta vs Phi dead zones per Run.
   * -
   */
  /* - [0] refers to Eta.
   * - I am plotting from -5.0 to -2.0,
   * - hence 150 bins are reasonable...  (REBIN 10x)
   * - [1] refers to Phi.
   * - To avoid problems related to TMath::Pi(),
   * - I am plotting from 0 to 8.
   * - Hence I had thought of 200 bins... (REBIN 10x)
   */
  Int_t listOfGoodRunNumbers[]       = { 295585, 295586, 295587, 295588, 295589, 295612,
                                         295615, 295665, 295666, 295667, 295668, 295671,
                                         295673, 295675, 295676, 295677, 295714, 295716,
                                         295717, 295718, 295719, 295723, 295725, 295753,
                                         295754, 295755, 295758, 295759, 295762, 295763,
                                         295786, 295788, 295791, 295816, 295818, 295819,
                                         295822, 295825, 295826, 295829, 295831, 295854,
                                         295855, 295856, 295859, 295860, 295861, 295863,
                                         295881, 295908, 295909, 295910, 295913, 295936,
                                         295937, 295941, 295942, 295943, 295945, 295947,
                                         296061, 296062, 296063, 296065, 296066, 296068,
                                         296123, 296128, 296132, 296133, 296134, 296135,
                                         296142, 296143, 296191, 296192, 296194, 296195,
                                         296196, 296197, 296198, 296241, 296242, 296243,
                                         296244, 296246, 296247, 296269, 296270, 296273,
                                         296279, 296280, 296303, 296304, 296307, 296309,
                                         296312, 296376, 296377, 296378, 296379, 296380,
                                         296381, 296383, 296414, 296419, 296420, 296423,
                                         296424, 296433, 296472, 296509, 296510, 296511,
                                         296514, 296516, 296547, 296548, 296549, 296550,
                                         296551, 296552, 296553, 296615, 296616, 296618,
                                         296619, 296622, 296623,
                                         296690, 296691, 296694, 296749, 296750, 296781,
                                         296784, 296785, 296786, 296787, 296791, 296793,
                                         296794, 296799, 296836, 296838, 296839, 296848,
                                         296849, 296850, 296851, 296852, 296890, 296894,
                                         296899, 296900, 296903, 296930, 296931, 296932,
                                         296934, 296935, 296938, 296941, 296966, 296967,
                                         296968, 296969, 296971, 296975, 296976, 296977,
                                         296979, 297029, 297031, 297035, 297085, 297117,
                                         297118, 297119, 297123, 297124, 297128, 297129,
                                         297132, 297133, 297193, 297194, 297196, 297218,
                                         297219, 297221, 297222, 297278, 297310, 297312,
                                         297315, 297317, 297363, 297366, 297367, 297372,
                                         297379, 297380, 297405, 297408, 297413, 297414,
                                         297415, 297441, 297442, 297446, 297450, 297451,
                                         297452, 297479, 297481, 297483, 297512, 297537,
                                         297540, 297541, 297542, 297544, 297558, 297588,
                                         297590, 297595,/*, 297623, 297624*/
                                         244918, 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
                                         245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
                                         245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
                                         245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
                                         245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
                                         245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
                                         246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
                                         246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
                                         246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
                                         246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
                                         246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
                                         246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
                                         246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
                                         246949, 246980, 246982, 246984, 246989, 246991, 246994
                                       };
  for( Int_t iRuns = 0; iRuns < 364; iRuns++ ) {
    fDeadZoneEtaVsPhiPerRunH[iRuns] = new TH2F( Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                                Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                                150, -5.0, -2.0,
                                                // 200, -4.0,  4.0
                                                200,  0.0,  8.0
                                                );
    fOutputList->Add(fDeadZoneEtaVsPhiPerRunH[iRuns]);
  }

  for( Int_t iRuns = 0; iRuns < 364; iRuns++ ) {
    fDeadZoneEtaVsPhiPerRunWithTriggeringH[iRuns] = new TH2F( Form( "fDeadZoneEtaVsPhiPerRunWithTriggeringH_%d", listOfGoodRunNumbers[iRuns] ),
                                                Form( "fDeadZoneEtaVsPhiPerRunWithTriggeringH_%d", listOfGoodRunNumbers[iRuns] ),
                                                150, -5.0, -2.0,
                                                // 200, -4.0,  4.0
                                                200,  0.0,  8.0
                                                );
    fOutputList->Add(fDeadZoneEtaVsPhiPerRunWithTriggeringH[iRuns]);
  }

  /* - ZDC energy spectra for calibration.
   * - Needed for XNXN analysis.
   * -
   */
  fZNCEnergyAgainstEntriesH = new TH1F("fZNCEnergyAgainstEntriesH", "fZNCEnergyAgainstEntriesH", 20000, -10000, 400000);
  fOutputList->Add(fZNCEnergyAgainstEntriesH);

  fZNAEnergyAgainstEntriesH = new TH1F("fZNAEnergyAgainstEntriesH", "fZNAEnergyAgainstEntriesH", 20000, -10000, 400000);
  fOutputList->Add(fZNAEnergyAgainstEntriesH);

  fZNCEnergyBeforeTimingSelectionH = new TH1F("fZNCEnergyBeforeTimingSelectionH", "fZNCEnergyBeforeTimingSelectionH", 20000, -10000, 400000);
  fOutputList->Add(fZNCEnergyBeforeTimingSelectionH);

  fZNAEnergyBeforeTimingSelectionH = new TH1F("fZNAEnergyBeforeTimingSelectionH", "fZNAEnergyBeforeTimingSelectionH", 20000, -10000, 400000);
  fOutputList->Add(fZNAEnergyBeforeTimingSelectionH);

  fZNCEnergyCalibratedH = new TH1F("fZNCEnergyCalibratedH", "fZNCEnergyCalibratedH", 20000, -10000, 400000);
  fOutputList->Add(fZNCEnergyCalibratedH);

  fZNAEnergyCalibratedH = new TH1F("fZNAEnergyCalibratedH", "fZNAEnergyCalibratedH", 20000, -10000, 400000);
  fOutputList->Add(fZNAEnergyCalibratedH);

  fZNCEnergyUncalibratedH = new TH1F("fZNCEnergyUncalibratedH", "fZNCEnergyUncalibratedH", 20000, -10000, 400000);
  fOutputList->Add(fZNCEnergyUncalibratedH);

  fZNAEnergyUncalibratedH = new TH1F("fZNAEnergyUncalibratedH", "fZNAEnergyUncalibratedH", 20000, -10000, 400000);
  fOutputList->Add(fZNAEnergyUncalibratedH);

  fZNCEnergyCalibratedHigherGainH = new TH1F("fZNCEnergyCalibratedHigherGainH", "fZNCEnergyCalibratedHigherGainH", 20000, -80000, 3200000);
  fOutputList->Add(fZNCEnergyCalibratedHigherGainH);

  fZNAEnergyCalibratedHigherGainH = new TH1F("fZNAEnergyCalibratedHigherGainH", "fZNAEnergyCalibratedHigherGainH", 20000, -80000, 3200000);
  fOutputList->Add(fZNAEnergyCalibratedHigherGainH);

  for( Int_t iRuns = 0; iRuns < 364; iRuns++ ) {
    fZNCEnergyPerRunH[iRuns] = new TH1F( Form( "fZNCEnergyPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                         Form( "fZNCEnergyPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                         20000, -10000, 400000
                                         );
    fOutputList->Add(fZNCEnergyPerRunH[iRuns]);
  }

  for( Int_t iRuns = 0; iRuns < 364; iRuns++ ) {
    fZNAEnergyPerRunH[iRuns] = new TH1F( Form( "fZNAEnergyPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                         Form( "fZNAEnergyPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                         20000, -10000, 400000
                                         );
    fOutputList->Add(fZNAEnergyPerRunH[iRuns]);
  }

  fZNCTimeAgainstEntriesH = new TH1F("fZNCTimeAgainstEntriesH", "fZNCTimeAgainstEntriesH", 6000, -150, 150);
  fOutputList->Add(fZNCTimeAgainstEntriesH);

  fZNATimeAgainstEntriesH = new TH1F("fZNATimeAgainstEntriesH", "fZNATimeAgainstEntriesH", 6000, -150, 150);
  fOutputList->Add(fZNATimeAgainstEntriesH);

  fSingleMuonPtDistributionH = new TH1F("fSingleMuonPtDistributionH", "fSingleMuonPtDistributionH", 4000, 0, 20);
  fOutputList->Add(fSingleMuonPtDistributionH);

  //_______________________________
  // - End of the function
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskMatchTriggerForward::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskMatchTriggerForward::UserExec(Option_t *)
{
  /* - This iSelectionCounter is used as a token. So at every passing step it is
     - increased by one. All the events are supposed to pass the first step
     - obviously, but only a few get to the end. This is effect is clearly
     - noticeable in fCounterH event with the small trial local sample.
     - Almost 160k possible events at the 0-th step, while only 2k at the 4th step.
   */
  Int_t iSelectionCounter = 0; // no selection applied yet
  fCounterH->Fill(iSelectionCounter); // entering UserExec
  iSelectionCounter++;

  // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
      PostData(1, fOutputList);
      return;
  }
  // fMCEvent = MCEvent();
  // if(!fMCEvent) {
  //     PostData(1, fOutputList);
  //     return;
  // }
  // if(fMCEvent) {
  //   fRunNum    = fAOD->GetRunNumber();
  //   SetLuminosityCap();
  //   fCounterGeneratedLevel[ fRunNum - 240000 ] += 1;
  //   // cout << "fCounterGeneratedLevel[ " << (fRunNum - 240000) << " ] = " << fCounterGeneratedLevel[ fRunNum - 240000 ] << endl;
  //   // if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( (Int_t)fLumiPerRun * (Int_t)40000 ) ) {
  //   if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( fLumiPerRun * 40000 ) ) {
  //         PostData(1, fOutputList);
  //         return;
  //   }
  //   ProcessMCParticles(fMCEvent);
  //   fMCEfficiencyPerRunH->Fill( Form("%d", fRunNum) , 1 );
  // }
  /* - Trigger selection:
   - here we verify which trigger was the event selected upon.
   - The useful triggers are only those needed for the CTRUE
   - analysis. Hence, if we cannot find them we return...
   -
 */

  Int_t fCtrue = -1;
  TString trigger = fAOD->GetFiredTriggerClasses();
  if (trigger.Contains("CINT7-B-NOPF-MUFAST")) fCtrue = 1;
  if (trigger.Contains("CINT7ZAC-B-NOPF-CENTNOTRD")) fCtrue = 1;
  // if (trigger.Contains("CTRUE-B")) fCtrue = 1;
  // if (trigger.Contains("CTRUE-A")) fCtrue = 2;
  // if (trigger.Contains("CTRUE-C")) fCtrue = 3;
  // if (trigger.Contains("CTRUE-E")) fCtrue = 4;
  if ( fCtrue == -1 ) {
    PostData(1, fOutputList);
    return;
  }
  // if (    !(trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
  //           trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
  //           trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
  //           trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
  //           trigger.Contains("CMUP13-B-NOPF-MUFAST")  )
  //         )  {
  //                 PostData(1, fOutputList);
  //                 return;
  //             }
  fCounterH->Fill(3);

  /* - We are now checking if there were any tracks. If there were at least one,
     - then the histogram gets filled again. If not we are returning. There
     - would be no point in going further.
   */
  Int_t nTracks(fAOD->GetNumberOfTracks());
  if(nTracks<1) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter); // At least one track
  iSelectionCounter++;


  //_______________________________
  // EVENT DATA EXTRACTION
  /* - Eugeny Krishen's event data extraction. I am trying to implement it.
     - The only thing I am a bit worried about is whether it should go before or
     - after the "nTracks<1" check... I will try and switch it if it sounds
     - better. These data are used for the event selection and maybe later on
     - for track selection, but I did not get to that part yet. If after all of
     - this I remember to do so, I will come back to this point and correct this
     - statement. If you find this part, please, keep in mind to check the
     - following.
   */

  /* - Event information:
     - Run Number, maybe to select the GOOD Runs and discard the others;
     - Number of Tracklets, these are in this case the SPD tracklets, so the
     - almost unit vector roughly 2 cm between two pixels of the SPD in different
     - layers.
   */
  fRunNum    = fAOD->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();

  /* - Trigger Inputs:
     - L0: ..... ;
     - L1: ..... .
   */
  // fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  // fL1inputs = fAOD->GetHeader()->GetL1TriggerInputs();

  /* - Past-future protection maps:
     - IR1: .... ;
     - IR2: .... .
   */
  // fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  // fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();

  /* - ZDC: we try to find the ZDC object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without ZDC's information. Or at least, this
     - is my impression of it (filling fCounterH). ZDC information:
     - fZem1Energy:
     - fZem2Energy:
     - fZNAEnergy:
     - fZNCEnergy:
     - fZPAEnergy:
     - fZPCEnergy:
     - fZNATime:
     - fZNCTime:
     - fZNATDC[i]:
     - fZNCTDC[i]:
     - fZPATDC[i]:
     - fZPCTDC[i]:
   */
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fZem1Energy = dataZDC->GetZEM1Energy();
  fZem2Energy = dataZDC->GetZEM2Energy();

  /* - Reset Event information.
   * -
   */
  fZNAEnergy  = -8999;
  fZNCEnergy  = -8999;
  fZPAEnergy  = -8999;
  fZPCEnergy  = -8999;

  fZNAEnergy  = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy  = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy  = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy  = dataZDC->GetZPCTowerEnergy()[0];

  fZNATime    = dataZDC->GetZNATime();
  fZNCTime    = dataZDC->GetZNCTime();

  /* - Reset Event information.
   * -
   */
  for (Int_t i=0;i<4;i++) fZNATDC[i] = -999;
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = -999;
  for (Int_t i=0;i<4;i++) fZPATDC[i] = -999;
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = -999;

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);

  /* - These lines are the calibration for the ZDC as provided by Evgeny Kryshen.
     -
   */
  // Bool_t calibrated = 0;
  // if ( fRunNum <  295726 ) calibrated = 1;
  // if ( fRunNum == 296509 ) calibrated = 1;
  // if ( fRunNum >  296689 ) calibrated = 1;
  // if ( fRunNum >  296695 ) calibrated = 0;
  // if ( fRunNum == 297219 ) calibrated = 1;
  // if ( fRunNum == 297221 ) calibrated = 1;
  // if ( fRunNum == 297415 ) calibrated = 1;
  //
  // if ( !calibrated ) {
  //   fZNAEnergy *= (2500./190.);
  //   fZNCEnergy *= (2500./190.);
  // }

  // V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  // Reset event info
  fBBCFlags = 0;
  fBGCFlags = 0;
  fBBAFlags = 0;
  fBGAFlags = 0;
  for (Int_t i=0; i<64; i++) {
    // get array of fired cells
    fBBFlag[i] = dataVZERO->GetBBFlag(i);
    fBGFlag[i] = dataVZERO->GetBGFlag(i);
  }

  for (Int_t i=0; i<32; i++){ // loop over cells
    fBBCFlags+=fBBFlag[i];
    fBGCFlags+=fBGFlag[i];
    fBBAFlags+=fBBFlag[i+32];
    fBGAFlags+=fBGFlag[i+32];
  }


  //_____________________________________
  // RUN SELECTION
  /* - This part is the run selection. We call the std::find() method of the
     - <algorithm> library STL. This returns a kTRUE if you find the value you
     - are looking for inside the vector. So, what happens is that I look for
     - the Run Numbers inside the vector containing them. If I cannot find
     - them I move on to the next event.
     -
   */
  // auto findRunNumber = std::find(  std::begin(fVectorGoodRunNumbers),
  //                                  std::end(fVectorGoodRunNumbers),
  //                                  fRunNum
  //                                  );
  // if (findRunNumber != std::end(fVectorGoodRunNumbers)) {
  //     // std::cout << "fVectorGoodRunNumbers DOES     contain: " << fRunNum << std::endl;
  // } else {
  //     PostData(1, fOutputList);
  //     return;
  // }

  fCounterH->Fill(15);
  Int_t listOfGoodRunNumbersLHC18q[] = { 295585, 295586, 295587, 295588, 295589, 295612,
                                         295615, 295665, 295666, 295667, 295668, 295671,
                                         295673, 295675, 295676, 295677, 295714, 295716,
                                         295717, 295718, 295719, 295723, 295725, 295753,
                                         295754, 295755, 295758, 295759, 295762, 295763,
                                         295786, 295788, 295791, 295816, 295818, 295819,
                                         295822, 295825, 295826, 295829, 295831, 295854,
                                         295855, 295856, 295859, 295860, 295861, 295863,
                                         295881, 295908, 295909, 295910, 295913, 295936,
                                         295937, 295941, 295942, 295943, 295945, 295947,
                                         296061, 296062, 296063, 296065, 296066, 296068,
                                         296123, 296128, 296132, 296133, 296134, 296135,
                                         296142, 296143, 296191, 296192, 296194, 296195,
                                         296196, 296197, 296198, 296241, 296242, 296243,
                                         296244, 296246, 296247, 296269, 296270, 296273,
                                         296279, 296280, 296303, 296304, 296307, 296309,
                                         296312, 296376, 296377, 296378, 296379, 296380,
                                         296381, 296383, 296414, 296419, 296420, 296423,
                                         296424, 296433, 296472, 296509, 296510, 296511,
                                         296514, 296516, 296547, 296548, 296549, 296550,
                                         296551, 296552, 296553, 296615, 296616, 296618,
                                         296619, 296622, 296623 };
  Int_t listOfGoodRunNumbersLHC18r[] = { 296690, 296691, 296694, 296749, 296750, 296781,
                                         296784, 296785, 296786, 296787, 296791, 296793,
                                         296794, 296799, 296836, 296838, 296839, 296848,
                                         296849, 296850, 296851, 296852, 296890, 296894,
                                         296899, 296900, 296903, 296930, 296931, 296932,
                                         296934, 296935, 296938, 296941, 296966, 296967,
                                         296968, 296969, 296971, 296975, 296976, 296977,
                                         296979, 297029, 297031, 297035, 297085, 297117,
                                         297118, 297119, 297123, 297124, 297128, 297129,
                                         297132, 297133, 297193, 297194, 297196, 297218,
                                         297219, 297221, 297222, 297278, 297310, 297312,
                                         297315, 297317, 297363, 297366, 297367, 297372,
                                         297379, 297380, 297405, 297408, 297413, 297414,
                                         297415, 297441, 297442, 297446, 297450, 297451,
                                         297452, 297479, 297481, 297483, 297512, 297537,
                                         297540, 297541, 297542, 297544, 297558, 297588,
                                         297590, 297595/*, 297623, 297624*/ };
  /* - This good run number list has been taken from the analysis
     - note of Kay's talk for DIS 2017, see:
     - https://alice-notes.web.cern.ch/system/files/notes/analysis/596/2017-Feb-08-analysis_note-2017-Feb-08-analysis-note.pdf
     -
   */
  Int_t listOfGoodRunNumbersLHC15o[] = { 244918, 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
                                         245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
                                         245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
                                         245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
                                         245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
                                         245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
                                         246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
                                         246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
                                         246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
                                         246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
                                         246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
                                         246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
                                         246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
                                         246949, 246980, 246982, 246984, 246989, 246991, 246994
                                       };
  Int_t listOfRunNumbersZDC[] = { 296244, 296750, 296849, 297219, 297481 };
  Bool_t checkIfGoodRun = kFALSE;
  for( Int_t iRunLHC18q = 0; iRunLHC18q < 129; iRunLHC18q++){
    if( fRunNum == listOfGoodRunNumbersLHC18q[iRunLHC18q] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC18r = 0; iRunLHC18r <  98; iRunLHC18r++){
    if( fRunNum == listOfGoodRunNumbersLHC18r[iRunLHC18r] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC15o = 0; iRunLHC15o < 137; iRunLHC15o++){
    if( fRunNum == listOfGoodRunNumbersLHC15o[iRunLHC15o] ) checkIfGoodRun = kTRUE;
  }
  // for( Int_t iRunZDC = 0; iRunZDC < 5; iRunZDC++){
  //   if( fRunNum == listOfRunNumbersZDC[iRunZDC] )           checkIfGoodRun = kTRUE;
  // }
  // cout << "fRunNum = " << fRunNum << "   and   checkIfGoodRun = " << checkIfGoodRun << endl;
  if(checkIfGoodRun != 1) {
       PostData(1, fOutputList);
       // cout << "OPS!" << endl;
       return;
  }
  fCounterH->Fill(17);

  // END RUN SELECTION
  //_____________________________________



  /* - We have to get the number of fired V0C cells. So firstly, we get the
     - boolean information about the hit cells for all V0. This is done through
     - the GetBBFlag(i) method, where 0<i<32 stands for the V0C cells and
     - 32<i<64 for the V0A cells. Then I thought the easiest way to check
     - whether the number of fired V0C cells is above 2 is just to add up the
     - boolean numbers for 0<i<32. Let's see.
   */
  fV0TotalNCells = 0;
  for(Int_t iV0Hits = 0; iV0Hits < 64; iV0Hits++) {
        fV0Hits[iV0Hits] = dataVZERO->GetBBFlag(iV0Hits);
        if(fV0Hits[iV0Hits] == kTRUE) {
              if(iV0Hits < 32) fV0TotalNCells += 1;
        }
  }

  // AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(dataAD) {
        fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();

        // Reset event info
        fBBCFlagsAD = 0;
        fBGCFlagsAD = 0;
        fBBAFlagsAD = 0;
        fBGAFlagsAD = 0;
        for(Int_t i=0; i<16; i++) {
          // get array of fired pads
          fBBFlagAD[i] = dataAD->GetBBFlag(i);
          fBGFlagAD[i] = dataAD->GetBGFlag(i);
        }

        for(Int_t i=0; i<4; i++) { // loop over pairs of pads
          if ( fBBFlagAD[i]   && fBBFlagAD[i+4]  ) fBBCFlagsAD++;
          if ( fBGFlagAD[i]   && fBGFlagAD[i+4]  ) fBGCFlagsAD++;
          if ( fBBFlagAD[i+8] && fBBFlagAD[i+12] ) fBBAFlagsAD++;
          if ( fBGFlagAD[i+8] && fBGFlagAD[i+12] ) fBGAFlagsAD++;
        }
  }
  // END EVENT DATA EXTRACTION
  //_______________________________
  // APPLY TRIGGER MC!
  // if(!IsTriggered()) {
  //   cout << "Ehm" ;
  //   PostData(1, fOutputList);
  //   return;
  // }
  fCounterH->Fill(iSelectionCounter); // right trigger found
  iSelectionCounter++;
  //_______________________________
  // EVENT SELECTION
  /* - This is Eugeny Krishen's event selection from the talk in 14/1/2019 for
     - the PWG-UD (UPC oriented) meeting. The event selection requires:
     - CMUP11-B triggers;
     - Maximum 2 V0C cells fired;
     - Empty V0A decision;
     - Empty ADA decision;
     - Empty ADC decision;
     - 0 tracklets in SPD;
     - Exactly 2 unlike-sign muons;
   */
  /* - Empty V0A decision
     - Empty ADA decision
     - Empty ADC decision
   */
  // if(fV0ADecision != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  // if(fADADecision != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  // if(fADCDecision != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  // /* - 0 tracklets in SPD
  //  */
  // if(fTracklets != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  // /* - Maximum 2 V0C cells fired.
  //  */
  // if( fV0TotalNCells > 2 ) {
  //      PostData(1, fOutputList);
  //      return;
  // }




  // loop over tracks and select good muons
  Int_t nGoodMuons = 0;
  for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
    // get track
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) return;

    // is it a good muon track?
    if(!track->IsMuonTrack()) continue;
    if(!fMuonTrackCuts->IsSelected(track)) continue;

    // increase counter
    nGoodMuons++;

    // fill muon info
    fEtaMuonH ->Fill(track->Eta());
    fRAbsMuonH->Fill(track->GetRAtAbsorberEnd());
    fEntriesAgainstRunNumberH->Fill(fRunNum);
    /* - This is the last part of my try to obtain a proper RunNumbers histogram...
       -
     */
    fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
    fEfficiencyPerRunH               ->Fill( Form("%d", fRunNum) , 1 );
    if( track->GetMatchTrigger() != 0 ) {
      fEfficiencyPerRunWithTriggeringH->Fill( Form("%d", fRunNum) , 1 );
      ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunWithTriggeringH_%d", fRunNum )) )->Fill( track->Eta(), track->Phi() );
    }
    ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d",                 fRunNum )) )->Fill( track->Eta(), track->Phi() );

    fSingleMuonPtDistributionH->Fill( track->Pt() );

  }




  /* - ZDC plots for calibration of the energy spectra.
   * -
   * -
   */
  Bool_t isZNAfired = kFALSE;
  Bool_t isZNCfired = kFALSE;
  Bool_t isZNAfiredStrict = kFALSE;
  Bool_t isZNCfiredStrict = kFALSE;
  Int_t  counterZNA = 0;
  Int_t  counterZNC = 0;
  /* - Note that in C++ the && and || operators "short-circuit". That means that
     - they only evaluate a parameter if required. If the first parameter to &&
     - is false, or the first to || is true, the rest will not be evaluated.
     - That means that writing:
     - if ( (isZNAfired == 0) && (...) )
     - should mean effectively
     - if ( isZNAfired != 0 ) continue;
     - hence it should be *at least* one hit!!!
     -
   */
  for(Int_t iZDC = 0; iZDC < 4 ; iZDC++) {
    if ( (isZNAfired == 0) && (fZNATDC[iZDC] > -2.) && (fZNATDC[iZDC] < 2.) ) {
      isZNAfired = kTRUE;
      /* - After mail with Chiara Oppedisano, it seems like the best way
         - to proceed is to firstly call the IsZNAfired() and then filling...
         -
         - If this doesn't appear in later pulls it is because this
         - doesn't seem to suit my case...
         -
       */
      if( dataZDC->IsZNAfired() ) fZNATimeAgainstEntriesH->Fill(fZNATDC[iZDC]);
      // fCounterZNAH->Fill(counterZNA);
    }
    if ( (isZNCfired == 0) && (fZNCTDC[iZDC] > -2.) && (fZNCTDC[iZDC] < 2.) ) {
      isZNCfired = kTRUE;
      if( dataZDC->IsZNCfired() ) fZNCTimeAgainstEntriesH->Fill(fZNCTDC[iZDC]);
      // fCounterZNCH->Fill(counterZNC);
    }
    counterZNA++;
    counterZNC++;
  }

  if ( isZNCfired != 0 ) {
    fZNCEnergyAgainstEntriesH->Fill(fZNCEnergy);
    // if ( calibrated == 0 ) fZNCEnergyUncalibratedH->Fill(fZNCEnergy);
    // if ( calibrated == 1 ) {
    //   fZNCEnergyCalibratedH          ->Fill( fZNCEnergy );
    //   fZNCEnergyCalibratedHigherGainH->Fill( dataZDC->GetZNCTowerEnergyLR()[0] );
    // }
    ((TH1F*) fOutputList->FindObject(Form( "fZNCEnergyPerRunH_%d", fRunNum )) )->Fill( fZNCEnergy );
  }
  fZNCEnergyBeforeTimingSelectionH->Fill(fZNCEnergy);
  if ( isZNAfired != 0 ) {
    fZNAEnergyAgainstEntriesH->Fill(fZNAEnergy);
    // if ( calibrated == 0 ) fZNAEnergyUncalibratedH->Fill(fZNAEnergy);
    // if ( calibrated == 1 ) {
    //   fZNAEnergyCalibratedH          ->Fill( fZNAEnergy );
    //   fZNAEnergyCalibratedHigherGainH->Fill( dataZDC->GetZNATowerEnergyLR()[0] );
    // }
    ((TH1F*) fOutputList->FindObject(Form( "fZNAEnergyPerRunH_%d", fRunNum )) )->Fill( fZNAEnergy );
  }
  fZNAEnergyBeforeTimingSelectionH->Fill(fZNAEnergy);







  fNumberMuonsH->Fill(nGoodMuons);



  // post the data
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskMatchTriggerForward::IsTriggered()
{
  /* - This function IS roughly speaking the trigger for the MC.
     - This code has been greatly inspired by David Horak's work.
     - It returns kTRUE if the trigger has been fired.
     -
   */
  Bool_t is0VBAfired = fBBAFlags > 0;
  Bool_t is0VBCfired = fBBCFlags > 0;
  Bool_t is0UBAfired = fBBAFlagsAD > 0;
  Bool_t is0UBCfired = fBBCFlagsAD > 0;
  // cout << "is0VBAfired = ( fBBAFlags   = " << fBBAFlags   << " ) > 0 => " << is0VBAfired << endl;
  // cout << "is0VBCfired = ( fBBCFlags   = " << fBBCFlags   << " ) > 0 => " << is0VBCfired << endl;
  // cout << "is0UBAfired = ( fBBAFlagsAD = " << fBBAFlagsAD << " ) > 0 => " << is0UBAfired << endl;
  // cout << "is0UBCfired = ( fBBCFlagsAD = " << fBBCFlagsAD << " ) > 0 => " << is0UBCfired << endl;
  if (!is0VBAfired && !is0UBAfired && !is0UBCfired ) return kTRUE;
  else return kFALSE;

}
//_____________________________________________________________________________
void AliAnalysisTaskMatchTriggerForward::ProcessMCParticles(AliMCEvent* fMCEventArg)
{


}
//_____________________________________________________________________________
void AliAnalysisTaskMatchTriggerForward::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
/* - There were plenty of ways to do this...
 * - However, recently the STL libraries were
 * - creating confusion on the LEGO framework
 * - (they didn't fire at all).
 * - This problem was not found on local, where
 * - things were working properly...
 * - So I am using the most barbaric C-style
 * - arrays/for...
 */
void AliAnalysisTaskMatchTriggerForward::SetLuminosityCap()
{
  fLumiPerRun = 0;
  /* - Here I am rounding up the number for 10k,
   * - so that I don't have to do tedious conversions...
   * - I am adding 1 entry to the number obtained by 40k,
   * - so that I am not missing anything...
   * -
   */
  if      ( fRunNum == 244980 ) { fLumiPerRun = 0.0505; }
  else if ( fRunNum == 244982 ) { fLumiPerRun = 0.0761; }
  else if ( fRunNum == 244983 ) { fLumiPerRun = 0.0291; }
  else if ( fRunNum == 245064 ) { fLumiPerRun = 0.1643; }
  else if ( fRunNum == 245066 ) { fLumiPerRun = 0.0236; }
  else if ( fRunNum == 245068 ) { fLumiPerRun = 0.0202; }
  else if ( fRunNum == 245145 ) { fLumiPerRun = 1.2115; }
  else if ( fRunNum == 245146 ) { fLumiPerRun = 1.3773; }
  else if ( fRunNum == 245151 ) { fLumiPerRun = 0.1469; }
  else if ( fRunNum == 245152 ) { fLumiPerRun = 0.1655; }
  else if ( fRunNum == 245231 ) { fLumiPerRun = 0.3084; }
  else if ( fRunNum == 245232 ) { fLumiPerRun = 1.0146; }
  else if ( fRunNum == 245233 ) { fLumiPerRun = 0.2373; }
  else if ( fRunNum == 245253 ) { fLumiPerRun = 0.3067; }
  else if ( fRunNum == 245259 ) { fLumiPerRun = 0.4893; }
  else if ( fRunNum == 245343 ) { fLumiPerRun = 0.7006; }
  else if ( fRunNum == 245345 ) { fLumiPerRun = 2.2153; }
  else if ( fRunNum == 245346 ) { fLumiPerRun = 0.2785; }
  else if ( fRunNum == 245347 ) { fLumiPerRun = 1.1752; }
  else if ( fRunNum == 245353 ) { fLumiPerRun = 1.6505; }
  else if ( fRunNum == 245401 ) { fLumiPerRun = 0.7485; }
  else if ( fRunNum == 245407 ) { fLumiPerRun = 2.0625; }
  else if ( fRunNum == 245409 ) { fLumiPerRun = 0.8705; }
  else if ( fRunNum == 245410 ) { fLumiPerRun = 0.1819; }
  else if ( fRunNum == 245446 ) { fLumiPerRun = 0.1261; }
  else if ( fRunNum == 245450 ) { fLumiPerRun = 0.2621; }
  else if ( fRunNum == 245496 ) { fLumiPerRun = 1.06;   }
  else if ( fRunNum == 245501 ) { fLumiPerRun = 1.334;  }
  else if ( fRunNum == 245504 ) { fLumiPerRun = 0.6492; }
  else if ( fRunNum == 245505 ) { fLumiPerRun = 0.3624; }
  else if ( fRunNum == 245507 ) { fLumiPerRun = 1.6192; }
  else if ( fRunNum == 245535 ) { fLumiPerRun = 1.3612; }
  else if ( fRunNum == 245540 ) { fLumiPerRun = 0.7121; }
  else if ( fRunNum == 245542 ) { fLumiPerRun = 1.1181; }
  else if ( fRunNum == 245543 ) { fLumiPerRun = 2.0169; }
  else if ( fRunNum == 245554 ) { fLumiPerRun = 1.7248; }
  else if ( fRunNum == 245683 ) { fLumiPerRun = 4.0406; }
  else if ( fRunNum == 245692 ) { fLumiPerRun = 1.9090; }
  else if ( fRunNum == 245700 ) { fLumiPerRun = 1.1167; }
  else if ( fRunNum == 245705 ) { fLumiPerRun = 0.3239; }
  else if ( fRunNum == 245729 ) { fLumiPerRun = 1.1548; }
  else if ( fRunNum == 245731 ) { fLumiPerRun = 3.3932; }
  else if ( fRunNum == 245738 ) { fLumiPerRun = 1.9485; }
  else if ( fRunNum == 245752 ) { fLumiPerRun = 1.2497; }
  else if ( fRunNum == 245759 ) { fLumiPerRun = 1.3785; }
  else if ( fRunNum == 245766 ) { fLumiPerRun = 1.1429; }
  else if ( fRunNum == 245775 ) { fLumiPerRun = 1.7326; }
  else if ( fRunNum == 245785 ) { fLumiPerRun = 0.5102; }
  else if ( fRunNum == 245793 ) { fLumiPerRun = 0.7093; }
  else if ( fRunNum == 245829 ) { fLumiPerRun = 1.958;  }
  else if ( fRunNum == 245831 ) { fLumiPerRun = 1.9939; }
  else if ( fRunNum == 245833 ) { fLumiPerRun = 0.3559; }
  else if ( fRunNum == 245949 ) { fLumiPerRun = 0.5652; }
  else if ( fRunNum == 245952 ) { fLumiPerRun = 3.0759; }
  else if ( fRunNum == 245954 ) { fLumiPerRun = 1.9965; }
  else if ( fRunNum == 245963 ) { fLumiPerRun = 2.2815; }
  else if ( fRunNum == 245996 ) { fLumiPerRun = 0.4644; }
  else if ( fRunNum == 246001 ) { fLumiPerRun = 3.5684; }
  else if ( fRunNum == 246003 ) { fLumiPerRun = 0.5803; }
  else if ( fRunNum == 246012 ) { fLumiPerRun = 0.7302; }
  else if ( fRunNum == 246036 ) { fLumiPerRun = 0.2143; }
  else if ( fRunNum == 246037 ) { fLumiPerRun = 1.7465; }
  else if ( fRunNum == 246042 ) { fLumiPerRun = 4.8713; }
  else if ( fRunNum == 246048 ) { fLumiPerRun = 0.3835; }
  else if ( fRunNum == 246049 ) { fLumiPerRun = 3.2666; }
  else if ( fRunNum == 246053 ) { fLumiPerRun = 1.7691; }
  else if ( fRunNum == 246087 ) { fLumiPerRun = 14.184; }
  else if ( fRunNum == 246089 ) { fLumiPerRun = 0.3296; }
  else if ( fRunNum == 246113 ) { fLumiPerRun = 1.4761; }
  else if ( fRunNum == 246115 ) { fLumiPerRun = 0.4514; }
  else if ( fRunNum == 246148 ) { fLumiPerRun = 5.3175; }
  else if ( fRunNum == 246151 ) { fLumiPerRun = 3.0605; }
  else if ( fRunNum == 246152 ) { fLumiPerRun = 0.4734; }
  else if ( fRunNum == 246153 ) { fLumiPerRun = 4.6676; }
  else if ( fRunNum == 246178 ) { fLumiPerRun = 0.8156; }
  else if ( fRunNum == 246181 ) { fLumiPerRun = 2.7526; }
  else if ( fRunNum == 246182 ) { fLumiPerRun = 2.2047; }
  else if ( fRunNum == 246217 ) { fLumiPerRun = 3.4663; }
  else if ( fRunNum == 246220 ) { fLumiPerRun = 0.682;  }
  else if ( fRunNum == 246222 ) { fLumiPerRun = 3.6826; }
  else if ( fRunNum == 246225 ) { fLumiPerRun = 1.2534; }
  else if ( fRunNum == 246272 ) { fLumiPerRun = 5.5294; }
  else if ( fRunNum == 246275 ) { fLumiPerRun = 1.242;  }
  else if ( fRunNum == 246276 ) { fLumiPerRun = 0.5871; }
  else if ( fRunNum == 246390 ) { fLumiPerRun = 0.0448; }
  else if ( fRunNum == 246391 ) { fLumiPerRun = 0.1446; }
  else if ( fRunNum == 246392 ) { fLumiPerRun = 0.1765; }
  else if ( fRunNum == 246424 ) { fLumiPerRun = 2.866;  }
  else if ( fRunNum == 246428 ) { fLumiPerRun = 0.4417; }
  else if ( fRunNum == 246431 ) { fLumiPerRun = 1.7836; }
  else if ( fRunNum == 246433 ) { fLumiPerRun = 0.4164; }
  else if ( fRunNum == 246434 ) { fLumiPerRun = 4.103;  }
  else if ( fRunNum == 246487 ) { fLumiPerRun = 0.7286; }
  else if ( fRunNum == 246488 ) { fLumiPerRun = 7.5895; }
  else if ( fRunNum == 246493 ) { fLumiPerRun = 1.3534; }
  else if ( fRunNum == 246495 ) { fLumiPerRun = 0.4100; }
  else if ( fRunNum == 246675 ) { fLumiPerRun = 2.3469; }
  else if ( fRunNum == 246676 ) { fLumiPerRun = 0.4794; }
  else if ( fRunNum == 246750 ) { fLumiPerRun = 2.0756; }
  else if ( fRunNum == 246751 ) { fLumiPerRun = 2.0419; }
  else if ( fRunNum == 246755 ) { fLumiPerRun = 1.4197; }
  else if ( fRunNum == 246757 ) { fLumiPerRun = 0.59;   }
  else if ( fRunNum == 246758 ) { fLumiPerRun = 1.626;  }
  else if ( fRunNum == 246759 ) { fLumiPerRun = 0.3335; }
  else if ( fRunNum == 246760 ) { fLumiPerRun = 1.1753; }
  else if ( fRunNum == 246763 ) { fLumiPerRun = 0.549;  }
  else if ( fRunNum == 246765 ) { fLumiPerRun = 0.3274; }
  else if ( fRunNum == 246804 ) { fLumiPerRun = 1.0208; }
  else if ( fRunNum == 246805 ) { fLumiPerRun = 3.1925; }
  else if ( fRunNum == 246806 ) { fLumiPerRun = 2.5555; }
  else if ( fRunNum == 246807 ) { fLumiPerRun = 2.5962; }
  else if ( fRunNum == 246808 ) { fLumiPerRun = 0.3101; }
  else if ( fRunNum == 246809 ) { fLumiPerRun = 2.4707; }
  else if ( fRunNum == 246844 ) { fLumiPerRun = 0.7657; }
  else if ( fRunNum == 246845 ) { fLumiPerRun = 1.4355; }
  else if ( fRunNum == 246846 ) { fLumiPerRun = 0.8986; }
  else if ( fRunNum == 246847 ) { fLumiPerRun = 1.7064; }
  else if ( fRunNum == 246851 ) { fLumiPerRun = 1.2170; }
  else if ( fRunNum == 246855 ) { fLumiPerRun = 1.3014; }
  else if ( fRunNum == 246859 ) { fLumiPerRun = 1.2397; }
  else if ( fRunNum == 246864 ) { fLumiPerRun = 2.4832; }
  else if ( fRunNum == 246865 ) { fLumiPerRun = 0.8111; }
  else if ( fRunNum == 246867 ) { fLumiPerRun = 1.5019; }
  else if ( fRunNum == 246871 ) { fLumiPerRun = 0.8713; }
  else if ( fRunNum == 246930 ) { fLumiPerRun = 0.5641; }
  else if ( fRunNum == 246937 ) { fLumiPerRun = 0.699;  }
  else if ( fRunNum == 246942 ) { fLumiPerRun = 1.0555; }
  else if ( fRunNum == 246945 ) { fLumiPerRun = 2.1676; }
  else if ( fRunNum == 246948 ) { fLumiPerRun = 0.8855; }
  else if ( fRunNum == 246949 ) { fLumiPerRun = 2.8978; }
  else if ( fRunNum == 246980 ) { fLumiPerRun = 7.1999; }
  else if ( fRunNum == 246982 ) { fLumiPerRun = 0.5146; }
  else if ( fRunNum == 246984 ) { fLumiPerRun = 4.143;  }
  else if ( fRunNum == 246989 ) { fLumiPerRun = 3.8342; }
  else if ( fRunNum == 246991 ) { fLumiPerRun = 0.4368; }
  else if ( fRunNum == 246994 ) { fLumiPerRun = 1.2329; }
  else if ( fRunNum == 295585 ) { fLumiPerRun = 0.0793; }
  else if ( fRunNum == 295586 ) { fLumiPerRun = 0.2386; }
  else if ( fRunNum == 295587 ) { fLumiPerRun = 0.1095; }
  else if ( fRunNum == 295588 ) { fLumiPerRun = 0.1358; }
  else if ( fRunNum == 295589 ) { fLumiPerRun = 0.2819; }
  else if ( fRunNum == 295612 ) { fLumiPerRun = 0.449;  }
  else if ( fRunNum == 295615 ) { fLumiPerRun = 0.0566; }
  else if ( fRunNum == 295665 ) { fLumiPerRun = 0.3349; }
  else if ( fRunNum == 295666 ) { fLumiPerRun = 0.3239; }
  else if ( fRunNum == 295667 ) { fLumiPerRun = 0.0970; }
  else if ( fRunNum == 295668 ) { fLumiPerRun = 0.1303; }
  else if ( fRunNum == 295671 ) { fLumiPerRun = 0.3259; }
  else if ( fRunNum == 295673 ) { fLumiPerRun = 0.3128; }
  else if ( fRunNum == 295675 ) { fLumiPerRun = 0.132;  }
  else if ( fRunNum == 295676 ) { fLumiPerRun = 0.3213; }
  else if ( fRunNum == 295677 ) { fLumiPerRun = 0.2652; }
  else if ( fRunNum == 295714 ) { fLumiPerRun = 0.3456; }
  else if ( fRunNum == 295716 ) { fLumiPerRun = 0.3389; }
  else if ( fRunNum == 295717 ) { fLumiPerRun = 0.2880; }
  else if ( fRunNum == 295718 ) { fLumiPerRun = 0.2567; }
  else if ( fRunNum == 295719 ) { fLumiPerRun = 0.2947; }
  else if ( fRunNum == 295723 ) { fLumiPerRun = 0.5064; }
  else if ( fRunNum == 295725 ) { fLumiPerRun = 0.8890; }
  else if ( fRunNum == 295753 ) { fLumiPerRun = 0.3846; }
  else if ( fRunNum == 295754 ) { fLumiPerRun = 0.7055; }
  else if ( fRunNum == 295755 ) { fLumiPerRun = 0.7585; }
  else if ( fRunNum == 295758 ) { fLumiPerRun = 1.8934; }
  else if ( fRunNum == 295759 ) { fLumiPerRun = 0.5331; }
  else if ( fRunNum == 295762 ) { fLumiPerRun = 0.2749; }
  else if ( fRunNum == 295763 ) { fLumiPerRun = 1.0282; }
  else if ( fRunNum == 295786 ) { fLumiPerRun = 0.7490; }
  else if ( fRunNum == 295788 ) { fLumiPerRun = 3.0237; }
  else if ( fRunNum == 295791 ) { fLumiPerRun = 0.8580; }
  else if ( fRunNum == 295816 ) { fLumiPerRun = 1.2056; }
  else if ( fRunNum == 295818 ) { fLumiPerRun = 0.1453; }
  else if ( fRunNum == 295819 ) { fLumiPerRun = 2.7474; }
  else if ( fRunNum == 295822 ) { fLumiPerRun = 2.2529; }
  else if ( fRunNum == 295825 ) { fLumiPerRun = 0.2558; }
  else if ( fRunNum == 295826 ) { fLumiPerRun = 1.5814; }
  else if ( fRunNum == 295829 ) { fLumiPerRun = 0.9351; }
  else if ( fRunNum == 295831 ) { fLumiPerRun = 0.7762; }
  else if ( fRunNum == 295854 ) { fLumiPerRun = 1.3119; }
  else if ( fRunNum == 295855 ) { fLumiPerRun = 1.7466; }
  else if ( fRunNum == 295856 ) { fLumiPerRun = 1.4700; }
  else if ( fRunNum == 295859 ) { fLumiPerRun = 1.0510; }
  else if ( fRunNum == 295860 ) { fLumiPerRun = 0.8341; }
  else if ( fRunNum == 295861 ) { fLumiPerRun = 1.0670; }
  else if ( fRunNum == 295863 ) { fLumiPerRun = 0.7279; }
  else if ( fRunNum == 295881 ) { fLumiPerRun = 0.7115; }
  else if ( fRunNum == 295908 ) { fLumiPerRun = 2.9261; }
  else if ( fRunNum == 295909 ) { fLumiPerRun = 0.7875; }
  else if ( fRunNum == 295910 ) { fLumiPerRun = 3.1843; }
  else if ( fRunNum == 295913 ) { fLumiPerRun = 3.1294; }
  else if ( fRunNum == 295936 ) { fLumiPerRun = 1.4736; }
  else if ( fRunNum == 295937 ) { fLumiPerRun = 0.4057; }
  else if ( fRunNum == 295941 ) { fLumiPerRun = 1.6767; }
  else if ( fRunNum == 295942 ) { fLumiPerRun = 1.9237; }
  else if ( fRunNum == 295943 ) { fLumiPerRun = 1.6747; }
  else if ( fRunNum == 295945 ) { fLumiPerRun = 2.0370; }
  else if ( fRunNum == 295947 ) { fLumiPerRun = 2.6337; }
  else if ( fRunNum == 296061 ) { fLumiPerRun = 1.2968; }
  else if ( fRunNum == 296062 ) { fLumiPerRun = 1.8083; }
  else if ( fRunNum == 296063 ) { fLumiPerRun = 2.6876; }
  else if ( fRunNum == 296065 ) { fLumiPerRun = 2.4473; }
  else if ( fRunNum == 296066 ) { fLumiPerRun = 0.7337; }
  else if ( fRunNum == 296068 ) { fLumiPerRun = 1.9812; }
  else if ( fRunNum == 296123 ) { fLumiPerRun = 0.5065; }
  else if ( fRunNum == 296128 ) { fLumiPerRun = 0.4455; }
  else if ( fRunNum == 296132 ) { fLumiPerRun = 1.312;  }
  else if ( fRunNum == 296133 ) { fLumiPerRun = 1.7321; }
  else if ( fRunNum == 296134 ) { fLumiPerRun = 3.9104; }
  else if ( fRunNum == 296135 ) { fLumiPerRun = 2.3412; }
  else if ( fRunNum == 296142 ) { fLumiPerRun = 1.7893; }
  else if ( fRunNum == 296143 ) { fLumiPerRun = 0.5340; }
  else if ( fRunNum == 296191 ) { fLumiPerRun = 5.0507; }
  else if ( fRunNum == 296192 ) { fLumiPerRun = 0.4974; }
  else if ( fRunNum == 296194 ) { fLumiPerRun = 2.8725; }
  else if ( fRunNum == 296195 ) { fLumiPerRun = 0.7376; }
  else if ( fRunNum == 296196 ) { fLumiPerRun = 2.352;  }
  else if ( fRunNum == 296197 ) { fLumiPerRun = 2.0691; }
  else if ( fRunNum == 296198 ) { fLumiPerRun = 0.8140; }
  else if ( fRunNum == 296241 ) { fLumiPerRun = 0.8459; }
  else if ( fRunNum == 296242 ) { fLumiPerRun = 0.9517; }
  else if ( fRunNum == 296243 ) { fLumiPerRun = 1.5674; }
  else if ( fRunNum == 296244 ) { fLumiPerRun = 8.3722; }
  else if ( fRunNum == 296246 ) { fLumiPerRun = 1.8351; }
  else if ( fRunNum == 296247 ) { fLumiPerRun = 1.1765; }
  else if ( fRunNum == 296269 ) { fLumiPerRun = 3.8392; }
  else if ( fRunNum == 296270 ) { fLumiPerRun = 1.5116; }
  else if ( fRunNum == 296273 ) { fLumiPerRun = 7.2237; }
  else if ( fRunNum == 296279 ) { fLumiPerRun = 0.4057; }
  else if ( fRunNum == 296280 ) { fLumiPerRun = 1.5066; }
  else if ( fRunNum == 296303 ) { fLumiPerRun = 2.006;  }
  else if ( fRunNum == 296304 ) { fLumiPerRun = 6.0965; }
  else if ( fRunNum == 296307 ) { fLumiPerRun = 2.9023; }
  else if ( fRunNum == 296309 ) { fLumiPerRun = 2.1026; }
  else if ( fRunNum == 296312 ) { fLumiPerRun = 2.1228; }
  else if ( fRunNum == 296377 ) { fLumiPerRun = 6.0666; }
  else if ( fRunNum == 296378 ) { fLumiPerRun = 5.3897; }
  else if ( fRunNum == 296379 ) { fLumiPerRun = 2.0969; }
  else if ( fRunNum == 296380 ) { fLumiPerRun = 2.8820; }
  else if ( fRunNum == 296381 ) { fLumiPerRun = 1.4418; }
  else if ( fRunNum == 296383 ) { fLumiPerRun = 1.5136; }
  else if ( fRunNum == 296414 ) { fLumiPerRun = 4.8766; }
  else if ( fRunNum == 296419 ) { fLumiPerRun = 2.7523; }
  else if ( fRunNum == 296420 ) { fLumiPerRun = 1.4132; }
  else if ( fRunNum == 296423 ) { fLumiPerRun = 1.5981; }
  else if ( fRunNum == 296424 ) { fLumiPerRun = 0.3864; }
  else if ( fRunNum == 296433 ) { fLumiPerRun = 4.0456; }
  else if ( fRunNum == 296472 ) { fLumiPerRun = 0.8632; }
  else if ( fRunNum == 296509 ) { fLumiPerRun = 2.9592; }
  else if ( fRunNum == 296510 ) { fLumiPerRun = 9.0673; }
  else if ( fRunNum == 296511 ) { fLumiPerRun = 2.5666; }
  else if ( fRunNum == 296514 ) { fLumiPerRun = 0.4898; }
  else if ( fRunNum == 296516 ) { fLumiPerRun = 0.6134; }
  else if ( fRunNum == 296547 ) { fLumiPerRun = 1.0834; }
  else if ( fRunNum == 296548 ) { fLumiPerRun = 1.3771; }
  else if ( fRunNum == 296549 ) { fLumiPerRun = 4.8645; }
  else if ( fRunNum == 296550 ) { fLumiPerRun = 3.9901; }
  else if ( fRunNum == 296551 ) { fLumiPerRun = 2.0214; }
  else if ( fRunNum == 296552 ) { fLumiPerRun = 0.4842; }
  else if ( fRunNum == 296553 ) { fLumiPerRun = 0.7091; }
  else if ( fRunNum == 296615 ) { fLumiPerRun = 1.5676; }
  else if ( fRunNum == 296616 ) { fLumiPerRun = 0.5399; }
  else if ( fRunNum == 296618 ) { fLumiPerRun = 1.7014; }
  else if ( fRunNum == 296619 ) { fLumiPerRun = 1.5613; }
  else if ( fRunNum == 296622 ) { fLumiPerRun = 0.7064; }
  else if ( fRunNum == 296623 ) { fLumiPerRun = 2.1442; }
  else if ( fRunNum == 296690 ) { fLumiPerRun = 6.8615; }
  else if ( fRunNum == 296691 ) { fLumiPerRun = 0.6511; }
  else if ( fRunNum == 296694 ) { fLumiPerRun = 5.1826; }
  else if ( fRunNum == 296749 ) { fLumiPerRun = 9.2413; }
  else if ( fRunNum == 296750 ) { fLumiPerRun = 8.2161; }
  else if ( fRunNum == 296781 ) { fLumiPerRun = 0.8179; }
  else if ( fRunNum == 296784 ) { fLumiPerRun = 2.98;   }
  else if ( fRunNum == 296785 ) { fLumiPerRun = 1.9085; }
  else if ( fRunNum == 296786 ) { fLumiPerRun = 0.7537; }
  else if ( fRunNum == 296787 ) { fLumiPerRun = 3.2190; }
  else if ( fRunNum == 296791 ) { fLumiPerRun = 0.7573; }
  else if ( fRunNum == 296793 ) { fLumiPerRun = 1.3317; }
  else if ( fRunNum == 296794 ) { fLumiPerRun = 3.1335; }
  else if ( fRunNum == 296799 ) { fLumiPerRun = 2.7149; }
  else if ( fRunNum == 296836 ) { fLumiPerRun = 1.5116; }
  else if ( fRunNum == 296838 ) { fLumiPerRun = 0.5432; }
  else if ( fRunNum == 296839 ) { fLumiPerRun = 2.9424; }
  else if ( fRunNum == 296848 ) { fLumiPerRun = 2.1628; }
  else if ( fRunNum == 296849 ) { fLumiPerRun = 11.469; }
  else if ( fRunNum == 296850 ) { fLumiPerRun = 2.7979; }
  else if ( fRunNum == 296851 ) { fLumiPerRun = 0.1392; }
  else if ( fRunNum == 296852 ) { fLumiPerRun = 0.9565; }
  else if ( fRunNum == 296890 ) { fLumiPerRun = 8.0545; }
  else if ( fRunNum == 296894 ) { fLumiPerRun = 4.6472; }
  else if ( fRunNum == 296899 ) { fLumiPerRun = 2.1355; }
  else if ( fRunNum == 296900 ) { fLumiPerRun = 2.7833; }
  else if ( fRunNum == 296903 ) { fLumiPerRun = 1.0391; }
  else if ( fRunNum == 296930 ) { fLumiPerRun = 1.4575; }
  else if ( fRunNum == 296931 ) { fLumiPerRun = 0.5292; }
  else if ( fRunNum == 296932 ) { fLumiPerRun = 1.1863; }
  else if ( fRunNum == 296934 ) { fLumiPerRun = 2.5917; }
  else if ( fRunNum == 296935 ) { fLumiPerRun = 4.4039; }
  else if ( fRunNum == 296938 ) { fLumiPerRun = 1.6678; }
  else if ( fRunNum == 296941 ) { fLumiPerRun = 2.9181; }
  else if ( fRunNum == 296966 ) { fLumiPerRun = 3.3611; }
  else if ( fRunNum == 296967 ) { fLumiPerRun = 0.8051; }
  else if ( fRunNum == 296968 ) { fLumiPerRun = 3.1905; }
  else if ( fRunNum == 296969 ) { fLumiPerRun = 1.8878; }
  else if ( fRunNum == 296971 ) { fLumiPerRun = 0.6907; }
  else if ( fRunNum == 296975 ) { fLumiPerRun = 7.3683; }
  else if ( fRunNum == 296976 ) { fLumiPerRun = 1.1175; }
  else if ( fRunNum == 296979 ) { fLumiPerRun = 1.0995; }
  else if ( fRunNum == 297029 ) { fLumiPerRun = 7.2370; }
  else if ( fRunNum == 297031 ) { fLumiPerRun = 6.0499; }
  else if ( fRunNum == 297035 ) { fLumiPerRun = 0.5705; }
  else if ( fRunNum == 297085 ) { fLumiPerRun = 0.9774; }
  else if ( fRunNum == 297117 ) { fLumiPerRun = 2.3096; }
  else if ( fRunNum == 297118 ) { fLumiPerRun = 2.43;   }
  else if ( fRunNum == 297119 ) { fLumiPerRun = 2.6870; }
  else if ( fRunNum == 297123 ) { fLumiPerRun = 3.2804; }
  else if ( fRunNum == 297124 ) { fLumiPerRun = 0.6395; }
  else if ( fRunNum == 297128 ) { fLumiPerRun = 2.411;  }
  else if ( fRunNum == 297129 ) { fLumiPerRun = 2.8300; }
  else if ( fRunNum == 297132 ) { fLumiPerRun = 2.8179; }
  else if ( fRunNum == 297133 ) { fLumiPerRun = 1.1454; }
  else if ( fRunNum == 297193 ) { fLumiPerRun = 7.5602; }
  else if ( fRunNum == 297194 ) { fLumiPerRun = 8.8428; }
  else if ( fRunNum == 297196 ) { fLumiPerRun = 2.1255; }
  else if ( fRunNum == 297218 ) { fLumiPerRun = 6.42;   }
  else if ( fRunNum == 297219 ) { fLumiPerRun = 10.531; }
  else if ( fRunNum == 297221 ) { fLumiPerRun = 2.8309; }
  else if ( fRunNum == 297222 ) { fLumiPerRun = 1.7175; }
  else if ( fRunNum == 297278 ) { fLumiPerRun = 0.6019; }
  else if ( fRunNum == 297310 ) { fLumiPerRun = 0.6701; }
  else if ( fRunNum == 297312 ) { fLumiPerRun = 2.4002; }
  else if ( fRunNum == 297315 ) { fLumiPerRun = 7.8271; }
  else if ( fRunNum == 297317 ) { fLumiPerRun = 4.3148; }
  else if ( fRunNum == 297363 ) { fLumiPerRun = 1.9122; }
  else if ( fRunNum == 297366 ) { fLumiPerRun = 2.1293; }
  else if ( fRunNum == 297367 ) { fLumiPerRun = 3.1548; }
  else if ( fRunNum == 297372 ) { fLumiPerRun = 3.2003; }
  else if ( fRunNum == 297379 ) { fLumiPerRun = 6.8050; }
  else if ( fRunNum == 297380 ) { fLumiPerRun = 1.5488; }
  else if ( fRunNum == 297405 ) { fLumiPerRun = 0.6007; }
  else if ( fRunNum == 297408 ) { fLumiPerRun = 4.1021; }
  else if ( fRunNum == 297413 ) { fLumiPerRun = 2.9907; }
  else if ( fRunNum == 297414 ) { fLumiPerRun = 2.2140; }
  else if ( fRunNum == 297415 ) { fLumiPerRun = 6.8227; }
  else if ( fRunNum == 297441 ) { fLumiPerRun = 5.0556; }
  else if ( fRunNum == 297442 ) { fLumiPerRun = 1.9878; }
  else if ( fRunNum == 297446 ) { fLumiPerRun = 8.1326; }
  else if ( fRunNum == 297450 ) { fLumiPerRun = 1.9518; }
  else if ( fRunNum == 297451 ) { fLumiPerRun = 1.3327; }
  else if ( fRunNum == 297452 ) { fLumiPerRun = 1.1512; }
  else if ( fRunNum == 297479 ) { fLumiPerRun = 7.7463; }
  else if ( fRunNum == 297481 ) { fLumiPerRun = 10.645; }
  else if ( fRunNum == 297483 ) { fLumiPerRun = 1.9505; }
  else if ( fRunNum == 297512 ) { fLumiPerRun = 1.5848; }
  else if ( fRunNum == 297537 ) { fLumiPerRun = 1.8096; }
  else if ( fRunNum == 297540 ) { fLumiPerRun = 0.6286; }
  else if ( fRunNum == 297541 ) { fLumiPerRun = 4.0120; }
  else if ( fRunNum == 297542 ) { fLumiPerRun = 1.5362; }
  else if ( fRunNum == 297544 ) { fLumiPerRun = 7.2900; }
  else if ( fRunNum == 297558 ) { fLumiPerRun = 0.4783; }
  else if ( fRunNum == 297588 ) { fLumiPerRun = 5.2912; }
  else if ( fRunNum == 297590 ) { fLumiPerRun = 3.06;   }

}
