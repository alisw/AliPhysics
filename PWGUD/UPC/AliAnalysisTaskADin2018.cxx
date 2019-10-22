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
#include "AliAODVertex.h"         // My addition, to use Eugeny Krishen's format


// my headers
#include "AliAnalysisTaskADin2018.h"



class AliAnalysisTaskADin2018;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskADin2018) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskADin2018::AliAnalysisTaskADin2018()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fADcheck(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fRunNumberTriggerCMUP11ClassH(0),
      fRunNumberTriggerCMUP11ClassProperlyH(0),
      fRunNumberTriggerCMUP26ClassH(0),
      fRunNumberTriggerCMUP26ClassProperlyH(0),
      fRunNumberTriggerCMUP6ClassH(0),
      fRunNumberTriggerCMUP6ClassProperlyH(0),
      fRunNumberTriggerCMUP10ClassH(0),
      fRunNumberTriggerCMUP10ClassProperlyH(0),
      fRunNumberTriggerCMUP13ClassH(0),
      fRunNumberTriggerCMUP13ClassProperlyH(0),
      fTriggersVsRunH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionCoherentShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentShiftMinusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusTwoH(0),
      fDimuonPtDistributionH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
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
      fV0TotalNCells(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskADin2018::AliAnalysisTaskADin2018(const char* name)
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fADcheck(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fRunNumberTriggerCMUP11ClassH(0),
      fRunNumberTriggerCMUP11ClassProperlyH(0),
      fRunNumberTriggerCMUP26ClassH(0),
      fRunNumberTriggerCMUP26ClassProperlyH(0),
      fRunNumberTriggerCMUP6ClassH(0),
      fRunNumberTriggerCMUP6ClassProperlyH(0),
      fRunNumberTriggerCMUP10ClassH(0),
      fRunNumberTriggerCMUP10ClassProperlyH(0),
      fRunNumberTriggerCMUP13ClassH(0),
      fRunNumberTriggerCMUP13ClassProperlyH(0),
      fTriggersVsRunH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0 },
      fInvariantMassDistributionCoherentShiftMinusTwoH(0),
      fInvariantMassDistributionCoherentShiftMinusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusOneH(0),
      fInvariantMassDistributionCoherentShiftPlusTwoH(0),
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentShiftMinusTwoH(0),
      fInvariantMassDistributionIncoherentShiftMinusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusOneH(0),
      fInvariantMassDistributionIncoherentShiftPlusTwoH(0),
      fDimuonPtDistributionH(0),
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
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
      fV0TotalNCells(0)
{

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
AliAnalysisTaskADin2018::~AliAnalysisTaskADin2018()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::UserCreateOutputObjects()
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
                                    AliMuonTrackCuts::kMuPdca    |
                                    AliMuonTrackCuts::kMuMatchLpt   );
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

  fInvariantMassDistributionH = new TH1F("fInvariantMassDistributionH", "fInvariantMassDistributionH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionH);

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
  // fEntriesAgainstRunNumberProperlyH->SetCanExtend(TH1::kAllAxes);
  fEntriesAgainstRunNumberProperlyH->LabelsDeflate();
  fOutputList->Add(fEntriesAgainstRunNumberProperlyH);

  fRunNumberTriggerCMUP11ClassH = new TH1F("fRunNumberTriggerCMUP11ClassH", "fRunNumberTriggerCMUP11ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP11ClassH);

  fRunNumberTriggerCMUP11ClassProperlyH = new TH1F("fRunNumberTriggerCMUP11ClassProperlyH", "fRunNumberTriggerCMUP11ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP11ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP11ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP11ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP11ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP11ClassProperlyH);

  fRunNumberTriggerCMUP26ClassH = new TH1F("fRunNumberTriggerCMUP26ClassH", "fRunNumberTriggerCMUP26ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP26ClassH);

  fRunNumberTriggerCMUP26ClassProperlyH = new TH1F("fRunNumberTriggerCMUP26ClassProperlyH", "fRunNumberTriggerCMUP26ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP26ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP26ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP26ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP26ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP26ClassProperlyH);

  fRunNumberTriggerCMUP6ClassH = new TH1F("fRunNumberTriggerCMUP6ClassH", "fRunNumberTriggerCMUP6ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP6ClassH);

  fRunNumberTriggerCMUP6ClassProperlyH = new TH1F("fRunNumberTriggerCMUP6ClassProperlyH", "fRunNumberTriggerCMUP6ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP6ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP6ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP6ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP6ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP6ClassProperlyH);

  fRunNumberTriggerCMUP10ClassH = new TH1F("fRunNumberTriggerCMUP10ClassH", "fRunNumberTriggerCMUP10ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP10ClassH);

  fRunNumberTriggerCMUP10ClassProperlyH = new TH1F("fRunNumberTriggerCMUP10ClassProperlyH", "fRunNumberTriggerCMUP10ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP10ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP10ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP10ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP10ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP10ClassProperlyH);

  fRunNumberTriggerCMUP13ClassH = new TH1F("fRunNumberTriggerCMUP13ClassH", "fRunNumberTriggerCMUP13ClassH", 10000, 290000, 300000);
  fOutputList->Add(fRunNumberTriggerCMUP13ClassH);

  fRunNumberTriggerCMUP13ClassProperlyH = new TH1F("fRunNumberTriggerCMUP13ClassProperlyH", "fRunNumberTriggerCMUP13ClassProperlyH", 10000, 290000, 300000);
  fRunNumberTriggerCMUP13ClassProperlyH->SetStats(0);
  fRunNumberTriggerCMUP13ClassProperlyH->SetFillColor(38);
  // fRunNumberTriggerCMUP13ClassProperlyH->SetCanExtend(TH1::kAllAxes);
  fRunNumberTriggerCMUP13ClassProperlyH->LabelsDeflate();
  fOutputList->Add(fRunNumberTriggerCMUP13ClassProperlyH);

  fTriggersVsRunH = new TH2F("fTriggersVsRunH","",5,0,5,60000,240000,300000);
  fOutputList->Add(fTriggersVsRunH);

  fInvariantMassDistributionCoherentH = new TH1F("fInvariantMassDistributionCoherentH", "fInvariantMassDistributionCoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentH);

  for( Int_t iRapidity = 0; iRapidity < 6; iRapidity++ ) {
    fInvariantMassDistributionCoherentRapidityBinsH[iRapidity]
            = new TH1F( Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidity),
                        Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidity),
                        2000, 0, 20);
    fOutputList->Add(fInvariantMassDistributionCoherentRapidityBinsH[iRapidity]);
  }


  fInvariantMassDistributionCoherentShiftMinusTwoH = new TH1F("fInvariantMassDistributionCoherentShiftMinusTwoH", "fInvariantMassDistributionCoherentShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftMinusTwoH);

  fInvariantMassDistributionCoherentShiftMinusOneH = new TH1F("fInvariantMassDistributionCoherentShiftMinusOneH", "fInvariantMassDistributionCoherentShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftMinusOneH);

  fInvariantMassDistributionCoherentShiftPlusOneH = new TH1F("fInvariantMassDistributionCoherentShiftPlusOneH", "fInvariantMassDistributionCoherentShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftPlusOneH);

  fInvariantMassDistributionCoherentShiftPlusTwoH = new TH1F("fInvariantMassDistributionCoherentShiftPlusTwoH", "fInvariantMassDistributionCoherentShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentShiftPlusTwoH);

  fInvariantMassDistributionIncoherentH = new TH1F("fInvariantMassDistributionIncoherentH", "fInvariantMassDistributionIncoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentH);

  fInvariantMassDistributionIncoherentShiftMinusTwoH = new TH1F("fInvariantMassDistributionIncoherentShiftMinusTwoH", "fInvariantMassDistributionIncoherentShiftMinusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftMinusTwoH);

  fInvariantMassDistributionIncoherentShiftMinusOneH = new TH1F("fInvariantMassDistributionIncoherentShiftMinusOneH", "fInvariantMassDistributionIncoherentShiftMinusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftMinusOneH);

  fInvariantMassDistributionIncoherentShiftPlusOneH = new TH1F("fInvariantMassDistributionIncoherentShiftPlusOneH", "fInvariantMassDistributionIncoherentShiftPlusOneH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftPlusOneH);

  fInvariantMassDistributionIncoherentShiftPlusTwoH = new TH1F("fInvariantMassDistributionIncoherentShiftPlusTwoH", "fInvariantMassDistributionIncoherentShiftPlusTwoH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentShiftPlusTwoH);

  fDimuonPtDistributionH = new TH1F("fDimuonPtDistributionH", "fDimuonPtDistributionH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionH);

  fDimuonPtDistributionShiftPlusOneH = new TH1F("fDimuonPtDistributionShiftPlusOneH", "fDimuonPtDistributionShiftPlusOneH", 4000, 0.02, 20.02);
  fOutputList->Add(fDimuonPtDistributionShiftPlusOneH);


  /* - These histograms have an EXTENDED range (0,20)->(0,40)
     -
   */
  fInvariantMassDistributionExtendedH = new TH1F("fInvariantMassDistributionExtendedH", "fInvariantMassDistributionExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionExtendedH);

  fInvariantMassDistributionCoherentExtendedH = new TH1F("fInvariantMassDistributionCoherentExtendedH", "fInvariantMassDistributionCoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionCoherentExtendedH);

  fInvariantMassDistributionIncoherentExtendedH = new TH1F("fInvariantMassDistributionIncoherentExtendedH", "fInvariantMassDistributionIncoherentExtendedH", 4000, 0, 40);
  fOutputList->Add(fInvariantMassDistributionIncoherentExtendedH);





  //_______________________________
  // - End of the function
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing your output to file
  // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::UserExec(Option_t *)
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
  fCounterH->Fill(iSelectionCounter); // AOD event found
  iSelectionCounter++;

  /* - Is it the right trigger?
   * -
   * -
   */
  TString trigger = fAOD->GetFiredTriggerClasses();
  if ( !( trigger.Contains("CMUP6-B-NOPF-MUFAST") )  )  {
          PostData(1, fOutputList);
          return;
  }
  fCounterH->Fill(iSelectionCounter); // right trigger found
  iSelectionCounter++;


  /* - The following lines concern the LUMI computation.
   * - What is being done is that we fill the histograms
   * - with the number of events which pass the relative
   * - trigger conditions.
   * -
   */
  // Int_t  counterForTrigger = 0;
  fRunNum    = fAOD->GetRunNumber();
  if ( trigger.Contains("CMUP6-B-NOPF-MUFAST") )  {
    fTriggersVsRunH->Fill( 2.5, fRunNum );
    fRunNumberTriggerCMUP6ClassH        ->Fill(fRunNum);
    fRunNumberTriggerCMUP6ClassProperlyH->Fill( Form("%d", fRunNum) , 1 );
    // counterForTrigger++;
  }


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
  // fRunNum    = fAOD->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();

  /* - Trigger Inputs:
     - L0: ..... ;
     - L1: ..... .
   */
  fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  fL1inputs = fAOD->GetHeader()->GetL1TriggerInputs();

  /* - Past-future protection maps:
     - IR1: .... ;
     - IR2: .... .
   */
  fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();

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
  Bool_t calibrated = 0;
  // if ( fRunNum <= 245068 ) {
  //   calibrated = 1;
  // } else if ( ( fRunNum > 245068 ) && ( fRunNum <  246995 ) ){
  //   calibrated = 0;
  // } else {
  //   calibrated = 1;
  // }

  // if ( calibrated == 0 ) {
  //   if( fRunNum <= 246994 ) {
  //     fZNAEnergy *= (2500./250.);
  //     fZNCEnergy *= (2500./250.);
  //   }
  // }

  /* - V0: we try to find the V0 object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without V0's information. Or at least, this
     - is my impression of it (filling fCounterH). V0 information:
     - fV0ADecision: ..... ;
     - fV0CDecision: ..... .
     -
     -
     -
     - Plot the V0 variables to try to understand whether it is cells we are
     - talking about or boolean responses or something else altogether.
  */
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
        PostData(1, fOutputList);
        return;
  }
  fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;
  fCounterH->Fill(12);


  fV0ADecision = dataVZERO->GetV0ADecision();
  fCounterH->Fill(13);
  fV0CDecision = dataVZERO->GetV0CDecision();
  fCounterH->Fill(14);


  //_____________________________________
  // RUN SELECTION
  /* - NOTE: total run selection.
   * -
   */
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
                                         296312, /*296376,*/ 296377, 296378, 296379, 296380,
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
                                         296968, 296969, 296971, 296975, 296976, /*296977,*/
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
  Int_t listOfGoodRunNumbersLHC15o[] = { /*244918,*/ 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
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
  Bool_t checkIfGoodRun = kFALSE;
  for( Int_t iRunLHC18q = 0; iRunLHC18q < 128; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 129; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 125; iRunLHC18q++){
    if( fRunNum == listOfGoodRunNumbersLHC18q[iRunLHC18q] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC18r = 0; iRunLHC18r <  97; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  98; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  82; iRunLHC18r++){
    if( fRunNum == listOfGoodRunNumbersLHC18r[iRunLHC18r] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC15o = 0; iRunLHC15o < 136/*137*/; iRunLHC15o++){
  // for( Int_t iRunLHC15o = 0; iRunLHC15o < 134; iRunLHC15o++){
    if( fRunNum == listOfGoodRunNumbersLHC15o[iRunLHC15o] ) checkIfGoodRun = kTRUE;
  }
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
     -
     - Weird fact: this doesn't seem to work... I have changed it so that if
     - the single cell has recorded a signal (kTRUE) then it adds up to the
     - total number of cells. Hope for the best.
     -
     - I am an idiot!!!!!! I have to reset the variable everytime!!!!
   */
  fV0TotalNCells = 0;
  for(Int_t iV0Hits = 0; iV0Hits < 64; iV0Hits++) {
        fV0Hits[iV0Hits] = dataVZERO->GetBBFlag(iV0Hits);
        if(fV0Hits[iV0Hits] == kTRUE) {
              // if(iV0Hits < 32) fV0TotalNCells += fV0Hits[iV0Hits];
              if(iV0Hits < 32) fV0TotalNCells += 1;
        }
        // std::cout << "fV0Hits[iV0Hits = " << iV0Hits << ", fRunNum=" << fRunNum << "] = " << fV0Hits[iV0Hits] << endl;
        // std::cout << "fV0TotalNCells (fRunNum = " << fRunNum << ") = " << fV0TotalNCells << endl;
  }
  fCounterH->Fill(18);

  /* - AD: we try to find the AD object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without AD's information. Or at least, this
     - is my impression of it (filling fCounterH). AD information:
     - fADADecision: small detector in ALICE, ADA and ADC at large distances;
     - fADCDecision: again, maybe check whether it is cells or boolean, same as V0.
  */
  // AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  fCounterH->Fill(19);
  if(dataAD) {
        fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;
        fCounterH->Fill(20);

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();
        fCounterH->Fill(21);
  }
  fCounterH->Fill(22);

  // END EVENT DATA EXTRACTION
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
  /* - CMUP11-B triggers: I have to check with my supervisor, but this requirement
     - may have already been satisfied with the requirements for the trigger info
   */
  /* - Maximum 2 V0C cells fired:
     -
   */
  /* - Empty V0A decision
     - Empty ADA decision
     - Empty ADC decision
   */
  if(fV0ADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADADecision != 0) {
       PostData(1, fOutputList);
       return;
  }

  /**
   * - This is the AD check.
   * - This is selected at the
   * - level of train set up...
   * - Needed for both the neutron
   * - emission analysis and the
   * - p-Pb analysis!!
   * -
   */
  if( fADcheck != 0){
      if(fADCDecision != 0) {
           PostData(1, fOutputList);
           return;
      }
  }
  /* - Empty V0C decision
   * - or at least in beam timing.
   */
  if( !(fV0CDecision == 0 || fV0CDecision == 1) ) {
       PostData(1, fOutputList);
       return;
  }
  /* - 0 tracklets in SPD
     - Is it like this?? Not too sure what fTracklets was!
   */
  // if(fTracklets != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  /* - Maximum 2 V0C cells fired.
     -
     - Trying a more readable and immediate approach.
   */
  // if( !(fV0TotalNCells < 2) ) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  if( fV0TotalNCells > 2 ) {
       PostData(1, fOutputList);
       return;
  }



  /* - We are finally at the starting point. We loop over the tracks and select
     - the good muons. Later on everything should happen in this loop. Let us
     - see what the future has in hold.
     -
     - Saturday: I moved the creation of the AliAODTrack* track outside of the
     - loop as it would have been otherwise created for each single iteration.
     - This could have caused massive memory issues especially to grid. I have
     - added a second AliAODTrack* track[2] to hold the second supposed muon.
     - Now this is ready to send the information to two TLorentzVectors to
     - obtain the invariant mass of the J/Psi through the Mag() method of the
     - class. Hope for the best.
   */
  Int_t nGoodMuons = 0;
  AliAODTrack* track[2];
  track[0]         = 0x0;
  track[1]         = 0x0;
  for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
    /* - This should be another form of event selection.
       - I am basically requesting the presence of TWO good muons only.
       - Later I will be checking whether of they are likesign or unlikesign.
     */
    if(nGoodMuons > 2) {
         PostData(1, fOutputList);
         return;
    }
    // get track
    // AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    // if(!track) return;
    //
    // // is it a good muon track?
    // if(!track->IsMuonTrack()) continue;
    // if(!fMuonTrackCuts->IsSelected(track)) continue;

    track[nGoodMuons] = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track[nGoodMuons]) return;

    // is it a good muon track?
    if(!track[nGoodMuons]->IsMuonTrack()) {
        // track[nGoodMuons] = 0x0;
        continue;
    }
    if(!fMuonTrackCuts->IsSelected(track[nGoodMuons])) {
        // track[nGoodMuons] = 0x0;
        continue;
    }

    // MUON SELECTION
    /* - This is Eugeny Krishen's MUON selection from the talk in 14/1/2019 for
       - the PWG-UD (UPC oriented) meeting. The event selection requires:
       - Muon trigger matching >=2 (1 GeV/c threshold);
       - (-4) < eta < (-2.5);
       - (17.5 cm) < R_{abs} < (89.5 cm);
       - p \times DCA cut;
    */

    // increase counter
    nGoodMuons++;

    // fill muon info
    // fEtaMuonH->Fill(track->Eta());
    // fRAbsMuonH->Fill(track->GetRAtAbsorberEnd());

  }
  /* - We need EXACTLY 2 good muons !!!!!
     -
   */
  if( nGoodMuons != 2 ) {
        PostData(1, fOutputList);
        return;
  }
  /* - Implementing the track cut on the unlike muons
   * -
   */
  if( (track[0]->Charge()) == (track[1]->Charge()) ) {
        PostData(1, fOutputList);
        return;
  }
  for(Int_t iFilling = 0; iFilling < nGoodMuons; iFilling++) {
        fEtaMuonH ->Fill(track[iFilling]->Eta());
        fRAbsMuonH->Fill(track[iFilling]->GetRAtAbsorberEnd());
  }
  // store muons
  fNumberMuonsH->Fill(nGoodMuons);
  fEntriesAgainstRunNumberH->Fill(fRunNum);
  /* - This is the last part of my try to obtain a proper RunNumbers histogram...
     -
   */
  fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
  if (nGoodMuons>0) fCounterH->Fill(iSelectionCounter); // At least one good muon
  iSelectionCounter++;



  /* - Finally the core!!!
   * - What will be happening is that we will instantiate TLorentzVectors to
   * - obtain the invariant mass of the dimuon system. If everything goes fine
   * - after this we should be able to obtain the peak of the J/Psi. But
   * - things never go as expected, so who knows!
   */
  TLorentzVector muons[2];
  TLorentzVector possibleJPsi;
  Double_t       chargeOfMuons[2];
  for(int indexMuon = 0; indexMuon < 2; indexMuon++) {
        muons[indexMuon].SetPtEtaPhiM(   track[indexMuon]->Pt(),
                                         track[indexMuon]->Eta(),
                                         track[indexMuon]->Phi(),
                                         TDatabasePDG::Instance()->GetParticle(13)->Mass()
                                       );
        possibleJPsi += muons[indexMuon];
        chargeOfMuons[indexMuon] = track[indexMuon]->Charge();
  }
  fInvariantMassDistributionH->Fill(possibleJPsi.Mag());
  fInvariantMassDistributionExtendedH->Fill(possibleJPsi.Mag());


  /* - Now we are evaluating the pt of the dimuon pair. Generally speaking,
     - if such a pt is less than 0.25 GeV/c then it fills the coherent
     - component, otherwise the incoherent component. At this point we may fill
     - even the dimuon pt distribution histogram and see if it looks like Kay's.
     -
   */
  Double_t ptOfTheDimuonPair = possibleJPsi.Pt();
  if( ptOfTheDimuonPair < 0.25) {
        fInvariantMassDistributionCoherentH->Fill(possibleJPsi.Mag());
        fInvariantMassDistributionCoherentExtendedH->Fill(possibleJPsi.Mag());
        if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.75 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[0]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.75 && possibleJPsi.Rapidity() <= -3.50 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[1]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.25 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[2]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.25 && possibleJPsi.Rapidity() <= -3.00 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[3]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.75 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[4]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -2.75 && possibleJPsi.Rapidity() <= -2.50 ) {
          fInvariantMassDistributionCoherentRapidityBinsH[5]->Fill(possibleJPsi.Mag());
        }
  } else {
        fInvariantMassDistributionIncoherentH->Fill(possibleJPsi.Mag());
        fInvariantMassDistributionIncoherentExtendedH->Fill(possibleJPsi.Mag());
  }
  // fDimuonPtDistributionH->Fill(ptOfTheDimuonPair);
  if ( (possibleJPsi.Mag() > 2.8) && (possibleJPsi.Mag() < 3.3) ) {
    fDimuonPtDistributionH            ->Fill(ptOfTheDimuonPair);
    fDimuonPtDistributionShiftPlusOneH->Fill(ptOfTheDimuonPair);
  }

  if( ptOfTheDimuonPair < 0.200 ) {
    fInvariantMassDistributionCoherentShiftMinusTwoH->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftMinusTwoH->Fill(possibleJPsi.Mag());
  }
  if( ptOfTheDimuonPair < 0.225 ) {
    fInvariantMassDistributionCoherentShiftMinusOneH->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftMinusOneH->Fill(possibleJPsi.Mag());
  }
  if( ptOfTheDimuonPair < 0.275 ) {
    fInvariantMassDistributionCoherentShiftPlusOneH ->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftPlusOneH ->Fill(possibleJPsi.Mag());
  }
  if( ptOfTheDimuonPair < 0.300 ) {
    fInvariantMassDistributionCoherentShiftPlusTwoH ->Fill(possibleJPsi.Mag());
  } else {
    fInvariantMassDistributionIncoherentShiftPlusTwoH ->Fill(possibleJPsi.Mag());
  }




  /* -
     - Now we are ordering the muons. The first muon will always be positive.
     - This is useful for the histograms...
   */
  TLorentzVector muonsCopy[2];
  TLorentzVector muonsCopy2[2];
  TLorentzVector possibleJPsiCopy;
  if( chargeOfMuons[0] > 0 ){
    muonsCopy[0]     = muons[0];
    muonsCopy[1]     = muons[1];
  } else if( chargeOfMuons[0] < 0 ){
    muonsCopy[0]     = muons[1];
    muonsCopy[1]     = muons[0];
  }
  muonsCopy2[0] = muonsCopy[0];
  muonsCopy2[1] = muonsCopy[1];
  possibleJPsiCopy = possibleJPsi;
  /* - If we are in the J/Psi peak, hence 2.8 < M < 3.3 GeV/c, AND if we are
     - in the coherent regime, so if the Pt < 0.25 GeV/c, we fill the plots.
     -
     - In the following note that the rapidity is well computed, so we are
     - dealing with negative values... -4.0 < Y < -2.5 !!!
     -
   */
  Double_t possibleJPsiCopyMag =possibleJPsiCopy.Mag();





  // post the data
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskADin2018::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
