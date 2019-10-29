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
#include "THnSparse.h"
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
#include "AliAnalysisTaskForMCpPb.h"


// headers for MC
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"



class AliAnalysisTaskForMCpPb;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskForMCpPb) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskForMCpPb::AliAnalysisTaskForMCpPb()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fInvariantMassDistributionRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fDimuonPtDistributionH(0),
      fTemplatePtDistributionH(0),
      fTemplatePtDistributionRapidityH{ 0, 0, 0 },
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
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
      fMCpdgCodesH(0),
      fMCpdgCodesOnlyPrimaryH(0),
      fMCphiGeneratedTruthH(0),
      fMCetaGeneratedTruthH(0),
      fMCpseudorapidityGeneratedTruthH(0),
      fMCptGeneratedTruthH(0),
      fMCphiDimuonGeneratedTruthH(0),
      fMCetaDimuonGeneratedTruthH(0),
      fMCpseudorapidityDimuonGeneratedTruthH(0),
      fMCptDimuonGeneratedTruthSingleMuonsH(0),
      fMCptDimuonGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH(0),
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
      fCounterUPCevent(0),
      fEfficiencyPerRunH(0),
      fMCEfficiencyPerRunH(0),
      fEtaAndPhi(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskForMCpPb::AliAnalysisTaskForMCpPb( const char* name )
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fMCEvent(0),
      fNumberMuonsH(0),
      fCounterH(0),
      fEtaMuonH(0),
      fRAbsMuonH(0),
      fInvariantMassDistributionH(0),
      fInvariantMassDistributionRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fEntriesAgainstRunNumberH(0),
      fEntriesAgainstRunNumberProperlyH(0),
      fInvariantMassDistributionCoherentH(0),
      fInvariantMassDistributionCoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fInvariantMassDistributionIncoherentH(0),
      fInvariantMassDistributionIncoherentRapidityBinsH{ 0, 0, 0, 0, 0, 0},
      fDimuonPtDistributionH(0),
      fTemplatePtDistributionH(0),
      fTemplatePtDistributionRapidityH{ 0, 0, 0 },
      fInvariantMassDistributionExtendedH(0),
      fInvariantMassDistributionCoherentExtendedH(0),
      fInvariantMassDistributionIncoherentExtendedH(0),
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
      fMCpdgCodesH(0),
      fMCpdgCodesOnlyPrimaryH(0),
      fMCphiGeneratedTruthH(0),
      fMCetaGeneratedTruthH(0),
      fMCpseudorapidityGeneratedTruthH(0),
      fMCptGeneratedTruthH(0),
      fMCphiDimuonGeneratedTruthH(0),
      fMCetaDimuonGeneratedTruthH(0),
      fMCpseudorapidityDimuonGeneratedTruthH(0),
      fMCptDimuonGeneratedTruthSingleMuonsH(0),
      fMCptDimuonGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiGeneratedTruthH(0),
      fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH(0),
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
      fCounterUPCevent(0),
      fEfficiencyPerRunH(0),
      fMCEfficiencyPerRunH(0),
      fEtaAndPhi(0)
{
    // FillGoodRunVector(fVectorGoodRunNumbers);
    for( Int_t iRun = 0; iRun < 60000; iRun++) {
      fCounterGeneratedLevel[iRun]   = 0;
      // fDeadZoneEtaVsPhiPerRunH[iRun] = 0x0;
    }
    // fDeadZoneEtaVsPhiPerRunH[60000] = 0x0;
    for( Int_t iRun = 0; iRun < 134; iRun++) {
      fDeadZoneEtaVsPhiPerRunH[iRun] = 0x0;
    }
    fEtaAndPhi = new Double_t[2];

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
AliAnalysisTaskForMCpPb::~AliAnalysisTaskForMCpPb()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
void AliAnalysisTaskForMCpPb::UserCreateOutputObjects()
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

  for( Int_t iRapidityBin = 0; iRapidityBin < 6; iRapidityBin++ ){
    fInvariantMassDistributionRapidityBinsH[iRapidityBin]
          = new TH1F( Form("fInvariantMassDistributionRapidityBinsH_%d", iRapidityBin),
                      Form("fInvariantMassDistributionRapidityBinsH_%d", iRapidityBin),
                      2000, 0, 20
                      );
    fOutputList->Add(fInvariantMassDistributionRapidityBinsH[iRapidityBin]);
  }

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

  fInvariantMassDistributionCoherentH = new TH1F("fInvariantMassDistributionCoherentH", "fInvariantMassDistributionCoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionCoherentH);

  for( Int_t iRapidityBin = 0; iRapidityBin < 6; iRapidityBin++ ){
    fInvariantMassDistributionCoherentRapidityBinsH[iRapidityBin]
          = new TH1F( Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidityBin),
                      Form("fInvariantMassDistributionCoherentRapidityBinsH_%d", iRapidityBin),
                      2000, 0, 20
                      );
    fOutputList->Add(fInvariantMassDistributionCoherentRapidityBinsH[iRapidityBin]);
  }

  fInvariantMassDistributionIncoherentH = new TH1F("fInvariantMassDistributionIncoherentH", "fInvariantMassDistributionIncoherentH", 2000, 0, 20);
  fOutputList->Add(fInvariantMassDistributionIncoherentH);

  for( Int_t iRapidityBin = 0; iRapidityBin < 6; iRapidityBin++ ){
    fInvariantMassDistributionIncoherentRapidityBinsH[iRapidityBin]
          = new TH1F( Form("fInvariantMassDistributionIncoherentRapidityBinsH_%d", iRapidityBin),
                      Form("fInvariantMassDistributionIncoherentRapidityBinsH_%d", iRapidityBin),
                      2000, 0, 20
                      );
    fOutputList->Add(fInvariantMassDistributionIncoherentRapidityBinsH[iRapidityBin]);
  }

  fDimuonPtDistributionH = new TH1F("fDimuonPtDistributionH", "fDimuonPtDistributionH", 4000, 0, 20);
  fOutputList->Add(fDimuonPtDistributionH);

  fTemplatePtDistributionH = new TH1F("fTemplatePtDistributionH", "fTemplatePtDistributionH", 4000, 0, 20);
  fOutputList->Add(fTemplatePtDistributionH);

  for( Int_t iRapidityBin = 0; iRapidityBin < 3; iRapidityBin++ ){
    fTemplatePtDistributionRapidityH[iRapidityBin] =
          new TH1F( Form( "fTemplatePtDistributionRapidityH_%d", iRapidityBin),
                    Form( "fTemplatePtDistributionRapidityH_%d", iRapidityBin),
                    4000, 0, 20
                    );
    fOutputList->Add(fTemplatePtDistributionRapidityH[iRapidityBin]);
  }

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
  // - MC-only plots
  fMCpdgCodesH = new TH1F("fMCpdgCodesH", "fMCpdgCodesH", 3, 0, 3);
  fMCpdgCodesH ->SetStats(0);
  fMCpdgCodesH ->SetFillColor(38);
  fMCpdgCodesH->LabelsDeflate();
  fOutputList->Add(fMCpdgCodesH);

  fMCpdgCodesOnlyPrimaryH = new TH1F("fMCpdgCodesOnlyPrimaryH", "fMCpdgCodesOnlyPrimaryH", 3, 0, 3);
  fMCpdgCodesOnlyPrimaryH ->SetStats(0);
  fMCpdgCodesOnlyPrimaryH ->SetFillColor(38);
  fMCpdgCodesOnlyPrimaryH->LabelsDeflate();
  fOutputList->Add(fMCpdgCodesOnlyPrimaryH);

  fMCinvariantMassDistrJPsiGeneratedTruthH = new TH1F("fMCinvariantMassDistrJPsiGeneratedTruthH", "fMCinvariantMassDistrJPsiGeneratedTruthH", 2000000, 0, 20);
  fOutputList->Add(fMCinvariantMassDistrJPsiGeneratedTruthH);

  fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH = new TH1F("fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH", "fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH", 20000, 0, 20);
  fOutputList->Add(fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH);

  fMCphiGeneratedTruthH = new TH1F("fMCphiGeneratedTruthH", "fMCphiGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCphiGeneratedTruthH);

  fMCetaGeneratedTruthH = new TH1F("fMCetaGeneratedTruthH", "fMCetaGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCetaGeneratedTruthH);

  fMCpseudorapidityGeneratedTruthH = new TH1F("fMCpseudorapidityGeneratedTruthH", "fMCpseudorapidityGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCpseudorapidityGeneratedTruthH);

  fMCptGeneratedTruthH = new TH1F("fMCptGeneratedTruthH", "fMCptGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCptGeneratedTruthH);

  fMCphiDimuonGeneratedTruthH = new TH1F("fMCphiDimuonGeneratedTruthH", "fMCphiDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCphiDimuonGeneratedTruthH);

  fMCetaDimuonGeneratedTruthH = new TH1F("fMCetaDimuonGeneratedTruthH", "fMCetaDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCetaDimuonGeneratedTruthH);

  fMCpseudorapidityDimuonGeneratedTruthH = new TH1F("fMCpseudorapidityDimuonGeneratedTruthH", "fMCpseudorapidityDimuonGeneratedTruthH", 2000, 0, 20);
  fOutputList->Add(fMCpseudorapidityDimuonGeneratedTruthH);

  fMCptDimuonGeneratedTruthSingleMuonsH = new TH1F("fMCptDimuonGeneratedTruthSingleMuonsH", "fMCptDimuonGeneratedTruthSingleMuonsH", 4000, 0, 20);
  fOutputList->Add(fMCptDimuonGeneratedTruthSingleMuonsH);

  fMCptDimuonGeneratedTruthH = new TH1F("fMCptDimuonGeneratedTruthH", "fMCptDimuonGeneratedTruthH", 4000, 0, 20);
  fOutputList->Add(fMCptDimuonGeneratedTruthH);

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

  /* - Eta vs Phi dead zones per Run.
   * -
   */
  Int_t*    bins = new Int_t[2];
  Double_t* xmin = new Double_t[2];
  Double_t* xmax = new Double_t[2];
  /* - [0] refers to Eta.
   * - I am plotting from -5.0 to -2.0,
   * - hence 150 bins are reasonable...  (REBIN 10x)
   * - [1] refers to Phi.
   * - To avoid problems related to TMath::Pi(),
   * - I am plotting from -4 to 4.
   * - Hence I had thought of 200 bins... (REBIN 10x)
   */
  bins[0] = 150;
  bins[1] = 200;
  xmin[0] = -5.0;
  xmin[1] = -4.0;
  xmax[0] = -2.0;
  xmax[1] = +4.0;
  Int_t listOfGoodRunNumbers[]       = {  266318, 266316, 266312, 266305, 266304, 266300, 266299, 266296, 266235, 266234,
                                          266208, 266197, 266196, 266193, 266190, 266189, 266187, 266117, 266086, 266085,
                                          266084, 266081, 266076, 266074, 266034, 266025, 266023, 266022, 265841, 265840,
                                          265797, 265795, 265792, 265789, 265788, 265787, 265785, 265756, 265754, 265746,
                                          265744, 265742, 265741, 265740, 265714, 265713, 265709, 265701, 265700, 265698,
                                          265697, 265696, 265694, 265691, 265607, 265596, 265594,
                                          267131, 267130, 267110, 267109, 267077, 267072, 267070, 267067, 267063, 267062,
                                          267022, 267020, 266998, 266997, 266994, 266993, 266988, 266944, 266943, 266942,
                                          266940, 266915, 266912, 266886, 266885, 266883, 266882, 266880, 266878, 266857,
                                          266807, 266805, 266800, 266776, 266775, 266708, 266706, 266703, 266702, 266676,
                                          266674, 266669, 266668, 266665, 266659, 266658, 266657, 266630, 266621, 266618,
                                          266615, 266614, 266613, 266595, 266593, 266591, 266588, 266587, 266584, 266549,
                                          266543, 266539, 266534, 266533, 266525, 266523, 266522, 266520, 266518, 266516,
                                          266514, 266487, 266480, 266479, 266472, 266441, 266439
                                       };
  for( Int_t iRuns = 0; iRuns < 134; iRuns++ ) {
    fDeadZoneEtaVsPhiPerRunH[iRuns] = new TH2F( Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                                Form( "fDeadZoneEtaVsPhiPerRunH_%d", listOfGoodRunNumbers[iRuns] ),
                                                150, -5.0, -2.0,
                                                200,  0.0,  8.0
                                                );
    fOutputList->Add(fDeadZoneEtaVsPhiPerRunH[iRuns]);
  }



  //_______________________________
  // - End of the function
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskForMCpPb::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskForMCpPb::UserExec(Option_t *)
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
  fMCEvent = MCEvent();
  if(!fMCEvent) {
      PostData(1, fOutputList);
      return;
  }
  if(fMCEvent) {
    fRunNum    = fAOD->GetRunNumber();
    SetLuminosityCap();
    fCounterGeneratedLevel[ fRunNum - 240000 ] += 1;
    // cout << "fCounterGeneratedLevel[ " << (fRunNum - 240000) << " ] = " << fCounterGeneratedLevel[ fRunNum - 240000 ] << endl;
    // if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( (Int_t)fLumiPerRun * (Int_t)40000 ) ) {
    if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( fLumiPerRun * 40000 ) ) {
          PostData(1, fOutputList);
          return;
    }
    ProcessMCParticles(fMCEvent);
    fCounterUPCevent += 1;
    fMCEfficiencyPerRunH->Fill( Form("%d", fRunNum) , 1 );
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
  fZNAEnergy  = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy  = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy  = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy  = dataZDC->GetZPCTowerEnergy()[0];

  fZNATime    = dataZDC->GetZNATime();
  fZNCTime    = dataZDC->GetZNCTime();

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);


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
  /* - NOTE: total run selection.
   * -
   */
  fCounterH->Fill(15);
  Int_t listOfGoodRunNumbersLHC16r[]  = { 266318, 266316, 266312, 266305, 266304, 266300, 266299, 266296, 266235, 266234,
                                          266208, 266197, 266196, 266193, 266190, 266189, 266187, 266117, 266086, 266085,
                                          266084, 266081, 266076, 266074, 266034, 266025, 266023, 266022, 265841, 265840,
                                          265797, 265795, 265792, 265789, 265788, 265787, 265785, 265756, 265754, 265746,
                                          265744, 265742, 265741, 265740, 265714, 265713, 265709, 265701, 265700, 265698,
                                          265697, 265696, 265694, 265691, 265607, 265596, 265594 };
  Int_t listOfGoodRunNumbersLHC16s[]  = { 267131, 267130, 267110, 267109, 267077, 267072, 267070, 267067, 267063, 267062,
                                          267022, 267020, 266998, 266997, 266994, 266993, 266988, 266944, 266943, 266942,
                                          266940, 266915, 266912, 266886, 266885, 266883, 266882, 266880, 266878, 266857,
                                          266807, 266805, 266800, 266776, 266775, 266708, 266706, 266703, 266702, 266676,
                                          266674, 266669, 266668, 266665, 266659, 266658, 266657, 266630, 266621, 266618,
                                          266615, 266614, 266613, 266595, 266593, 266591, 266588, 266587, 266584, 266549,
                                          266543, 266539, 266534, 266533, 266525, 266523, 266522, 266520, 266518, 266516,
                                          266514, 266487, 266480, 266479, 266472, 266441, 266439 };
  Bool_t checkIfGoodRun = kFALSE;
  for( Int_t iRunLHC16r = 0; iRunLHC16r <  57; iRunLHC16r++){
    if( fRunNum == listOfGoodRunNumbersLHC16r[iRunLHC16r] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC16s = 0; iRunLHC16s <  77; iRunLHC16s++){
    if( fRunNum == listOfGoodRunNumbersLHC16s[iRunLHC16s] ) checkIfGoodRun = kTRUE;
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
  if(!IsTriggered()) {
    cout << "Ehm" ;
    PostData(1, fOutputList);
    return;
  }
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
  if(fV0ADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADADecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  if(fADCDecision != 0) {
       PostData(1, fOutputList);
       return;
  }
  /* - 0 tracklets in SPD
   */
  if(fTracklets != 0) {
       PostData(1, fOutputList);
       return;
  }
  /* - Maximum 2 V0C cells fired.
   */
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
  fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );
  fEfficiencyPerRunH               ->Fill( Form("%d", fRunNum) , 1 );
  if (nGoodMuons>0) fCounterH->Fill(iSelectionCounter); // At least one good muon
  iSelectionCounter++;

  /* - Filling the fDeadZoneEtaVsPhiPerRunH.
   * -
   */
  ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( track[0]->Eta(), track[0]->Phi() );
  ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( track[1]->Eta(), track[1]->Phi() );

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
  if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.75 ) {
    fInvariantMassDistributionRapidityBinsH[0]->Fill(possibleJPsi.Mag());
  } else if ( possibleJPsi.Rapidity() > -3.75 && possibleJPsi.Rapidity() <= -3.50 ) {
    fInvariantMassDistributionRapidityBinsH[1]->Fill(possibleJPsi.Mag());
  } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.25 ) {
    fInvariantMassDistributionRapidityBinsH[2]->Fill(possibleJPsi.Mag());
  } else if ( possibleJPsi.Rapidity() > -3.25 && possibleJPsi.Rapidity() <= -3.00 ) {
    fInvariantMassDistributionRapidityBinsH[3]->Fill(possibleJPsi.Mag());
  } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.75 ) {
    fInvariantMassDistributionRapidityBinsH[4]->Fill(possibleJPsi.Mag());
  } else if ( possibleJPsi.Rapidity() > -2.75 && possibleJPsi.Rapidity() <= -2.50 ) {
    fInvariantMassDistributionRapidityBinsH[5]->Fill(possibleJPsi.Mag());
  }
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
        if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.75 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[0]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.75 && possibleJPsi.Rapidity() <= -3.50 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[1]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.25 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[2]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.25 && possibleJPsi.Rapidity() <= -3.00 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[3]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.75 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[4]->Fill(possibleJPsi.Mag());
        } else if ( possibleJPsi.Rapidity() > -2.75 && possibleJPsi.Rapidity() <= -2.50 ) {
          fInvariantMassDistributionIncoherentRapidityBinsH[5]->Fill(possibleJPsi.Mag());
        }
  }
  fDimuonPtDistributionH->Fill(ptOfTheDimuonPair);
  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
    fTemplatePtDistributionH->Fill(ptOfTheDimuonPair);
    if (        possibleJPsi.Rapidity() > -4.0  && possibleJPsi.Rapidity() <= -3.50 ) {
      fTemplatePtDistributionRapidityH[0]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.50 && possibleJPsi.Rapidity() <= -3.00 ) {
      fTemplatePtDistributionRapidityH[1]->Fill(possibleJPsi.Mag());
    } else if ( possibleJPsi.Rapidity() > -3.00 && possibleJPsi.Rapidity() <= -2.50 ) {
      fTemplatePtDistributionRapidityH[2]->Fill(possibleJPsi.Mag());
    }
  }


  /* - Filling the J/Psi's polarization plots.
     -
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
  muonsCopy2[0]      = muonsCopy[0];
  muonsCopy2[1]      = muonsCopy[1];
  possibleJPsiCopy = possibleJPsi;


  /* - If we are in the J/Psi peak, hence 2.8 < M < 3.3 GeV/c, AND if we are
     - in the coherent regime, so if the Pt < 0.25 GeV/c, we fill the plots.
     -
     - In the following note that the rapidity is well computed, so we are
     - dealing with negative values... -4.0 < Y < -2.5 !!!
     -
   */
  if ( (possibleJPsiCopy.Mag() > 2.85) && (possibleJPsiCopy.Mag() < 3.35) && (possibleJPsiCopy.Pt() < 0.25) ) {










  }







  // post the data
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskForMCpPb::IsTriggered()
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
void AliAnalysisTaskForMCpPb::ProcessMCParticles(AliMCEvent* fMCEventArg)
{
  // Loop over the MC tracks
  for(int iPart = 0; iPart < (fMCEventArg->GetNumberOfTracks()); iPart++) {
    AliAODMCParticle *mcParticle  = (AliAODMCParticle*)fMCEventArg->GetTrack(iPart);
    if (!mcParticle) {
      AliError(Form("Could not receive track %d", iPart));
      continue;
    }
    if (  mcParticle->Charge()    == 0) continue;
    Double_t pT     = mcParticle->Pt();
    Double_t eta    = mcParticle->Eta();
    Double_t pseudo = mcParticle->Y();
    Double_t phi    = mcParticle->Phi();
    Short_t  charge = mcParticle->Charge();
    Int_t    pdg    = mcParticle->PdgCode();
    fMCpdgCodesH->Fill( Form("%d", pdg) , 1 );
    if ( mcParticle->IsPrimary() ) fMCpdgCodesOnlyPrimaryH->Fill( Form("%d", pdg) , 1 );
    fMCphiGeneratedTruthH           ->Fill(phi);
    fMCetaGeneratedTruthH           ->Fill(eta);
    fMCpseudorapidityGeneratedTruthH->Fill(pseudo);
    fMCptGeneratedTruthH            ->Fill(pT);
  }

  /* - This is loop where we record the MC TRUTH stuff BEFORE both
     - the event selection and the track selection.
     - It is nonetheless strange that we can't directly use the
     - J/Psi and that we have to o through such a roundabout
     - road like reconstructing the invariant mass...
     -
   */
  Int_t nGoodMuonsMC = 0;
  TLorentzVector muonsMC[2];
  TLorentzVector possibleJPsiMC;
  Double_t pT[2];
  Double_t eta[2];
  Double_t pseudo[2];
  Double_t phi[2];
  Short_t  charge[2];
  Int_t    pdg[2];
  // Loop over the MC tracks
  for(int iPart = 0; iPart < (fMCEventArg->GetNumberOfTracks()); iPart++) {
    if(nGoodMuonsMC > 2) continue;
    AliAODMCParticle *mcParticle  = (AliAODMCParticle*)fMCEventArg->GetTrack(iPart);
    if (!mcParticle) {
      AliError(Form("Could not receive track %d", iPart));
      continue;
    }
    if (  !mcParticle->IsPrimary()                 ) continue;
    if (   mcParticle->Charge()              == 0  ) continue;
    if (   TMath::Abs(mcParticle->PdgCode()) == 13 ) {
      if ( nGoodMuonsMC < 2 ) {
        // cout << "Ok" << nGoodMuonsMC << endl;
        pT[nGoodMuonsMC]     = mcParticle->Pt();
        eta[nGoodMuonsMC]    = mcParticle->Eta();
        pseudo[nGoodMuonsMC] = mcParticle->Y();
        phi[nGoodMuonsMC]    = mcParticle->Phi();
        charge[nGoodMuonsMC] = mcParticle->Charge();
        pdg[nGoodMuonsMC]    = mcParticle->PdgCode();
        muonsMC[nGoodMuonsMC].SetPtEtaPhiM( pT[nGoodMuonsMC],
                                            eta[nGoodMuonsMC],
                                            phi[nGoodMuonsMC],
                                            TDatabasePDG::Instance()->GetParticle(13)->Mass()
                                           );
        nGoodMuonsMC++;
      }
    }
  }
  // cout << "Tracks in this event:  " << fMCEventArg->GetNumberOfTracks() << endl;
  /* - We need EXACTLY 2 good muons !!!!!
     -
   */
  if( nGoodMuonsMC == 2 ) {
    /* - Unlike-sign muons
       -
     */
    if( charge[0] != charge[1] ) {
      TLorentzVector muonsMCcopy[2];
      muonsMCcopy[0] = muonsMC[0];
      muonsMCcopy[1] = muonsMC[1];
      for( Int_t iMuonsMC = 0; iMuonsMC < 2; iMuonsMC++ ) {
        possibleJPsiMC+=muonsMC[iMuonsMC];
        fMCphiDimuonGeneratedTruthH           ->Fill(phi[iMuonsMC]);
        fMCetaDimuonGeneratedTruthH           ->Fill(eta[iMuonsMC]);
        fMCpseudorapidityDimuonGeneratedTruthH->Fill(pseudo[iMuonsMC]);
        fMCptDimuonGeneratedTruthSingleMuonsH ->Fill(pT[iMuonsMC]);
      }
      fMCinvariantMassDistrJPsiGeneratedTruthH->Fill(possibleJPsiMC.Mag());
      fMCptDimuonGeneratedTruthH->Fill(possibleJPsiMC.Pt());
      if ( (possibleJPsiMC.Mag() > 2.8) && (possibleJPsiMC.Mag() < 3.3) && (possibleJPsiMC.Pt() < 0.25) ) {
          if( charge[0] > 0 ) {
                  /* - This means that [0] is the positive muon while [1]
                   * - is the negative muon!
                   * -
                   */



          } else  {
                  /* - This means that [0] is the negative muon while [1]
                   * - is the positive muon!
                   * -
                   */


          }
      }
    }
  }

}
//_____________________________________________________________________________
void AliAnalysisTaskForMCpPb::Terminate(Option_t *)
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
void AliAnalysisTaskForMCpPb::SetLuminosityCap()
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
  else                          { fLumiPerRun = 1.00;   }

}
