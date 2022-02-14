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

#include <bitset>

// my headers
#include "AliAnalysisTaskForMCpPb.h"


// headers for MC
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"



class AliAnalysisTaskForMCpPb;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'
typedef std::bitset<32> IntBits;

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
      fEtaDimuonH(0),
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
      fEfficiencyPerRunRestrictedRapidityH(0),
      fMCEfficiencyPerRunRestrictedRapidityH(0),
      fEfficiencyPerRunRestrictedRapidity36to31H(0),
      fMCEfficiencyPerRunRestrictedRapidity36to31H(0),
      fEfficiencyPerRunRestrictedRapidity31to26H(0),
      fMCEfficiencyPerRunRestrictedRapidity31to26H(0),
      fEfficiencyPerRunWithRunTwoSettings(0),
      fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH{0,0},
      fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH{0,0,0},
      fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH{0,0,0,0},
      fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH{0,0,0,0,0},
      fMCEfficiencyPerRunWithRunTwoSettings(0),
      fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH{0,0},
      fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH{0,0,0},
      fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH{0,0,0,0},
      fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH{0,0,0,0,0},
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
      fEtaDimuonH(0),
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
      fEfficiencyPerRunRestrictedRapidityH(0),
      fMCEfficiencyPerRunRestrictedRapidityH(0),
      fEfficiencyPerRunRestrictedRapidity36to31H(0),
      fMCEfficiencyPerRunRestrictedRapidity36to31H(0),
      fEfficiencyPerRunRestrictedRapidity31to26H(0),
      fMCEfficiencyPerRunRestrictedRapidity31to26H(0),
      fEfficiencyPerRunWithRunTwoSettings(0),
      fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH{0,0},
      fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH{0,0,0},
      fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH{0,0,0,0},
      fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH{0,0,0,0,0},
      fMCEfficiencyPerRunWithRunTwoSettings(0),
      fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH{0,0},
      fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH{0,0,0},
      fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH{0,0,0,0},
      fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH{0,0,0,0,0},
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

  fEtaMuonH = new TH1F("fEtaMuonH", "fEtaMuonH", 160, -1, -5);
  fOutputList->Add(fEtaMuonH);

  fEtaDimuonH = new TH1F("fEtaDimuonH", "fEtaDimuonH", 160, -1, -5);
  fOutputList->Add(fEtaDimuonH);

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

  fEfficiencyPerRunRestrictedRapidityH = new TH1F("fEfficiencyPerRunRestrictedRapidityH", "fEfficiencyPerRunRestrictedRapidityH", 3, 0, 3);
  fEfficiencyPerRunRestrictedRapidityH->SetStats(0);
  fEfficiencyPerRunRestrictedRapidityH->SetFillColor(38);
  fEfficiencyPerRunRestrictedRapidityH->LabelsDeflate();
  fOutputList->Add(fEfficiencyPerRunRestrictedRapidityH);

  fMCEfficiencyPerRunRestrictedRapidityH = new TH1F("fMCEfficiencyPerRunRestrictedRapidityH", "fMCEfficiencyPerRunRestrictedRapidityH", 3, 0, 3);
  fMCEfficiencyPerRunRestrictedRapidityH->SetStats(0);
  fMCEfficiencyPerRunRestrictedRapidityH->SetFillColor(38);
  fMCEfficiencyPerRunRestrictedRapidityH->LabelsDeflate();
  fOutputList->Add(fMCEfficiencyPerRunRestrictedRapidityH);

  fEfficiencyPerRunRestrictedRapidity36to31H = new TH1F("fEfficiencyPerRunRestrictedRapidity36to31H", "fEfficiencyPerRunRestrictedRapidity36to31H", 3, 0, 3);
  fEfficiencyPerRunRestrictedRapidity36to31H->SetStats(0);
  fEfficiencyPerRunRestrictedRapidity36to31H->SetFillColor(38);
  fEfficiencyPerRunRestrictedRapidity36to31H->LabelsDeflate();
  fOutputList->Add(fEfficiencyPerRunRestrictedRapidity36to31H);

  fMCEfficiencyPerRunRestrictedRapidity36to31H = new TH1F("fMCEfficiencyPerRunRestrictedRapidity36to31H", "fMCEfficiencyPerRunRestrictedRapidity36to31H", 3, 0, 3);
  fMCEfficiencyPerRunRestrictedRapidity36to31H->SetStats(0);
  fMCEfficiencyPerRunRestrictedRapidity36to31H->SetFillColor(38);
  fMCEfficiencyPerRunRestrictedRapidity36to31H->LabelsDeflate();
  fOutputList->Add(fMCEfficiencyPerRunRestrictedRapidity36to31H);

  fEfficiencyPerRunRestrictedRapidity31to26H = new TH1F("fEfficiencyPerRunRestrictedRapidity31to26H", "fEfficiencyPerRunRestrictedRapidity31to26H", 3, 0, 3);
  fEfficiencyPerRunRestrictedRapidity31to26H->SetStats(0);
  fEfficiencyPerRunRestrictedRapidity31to26H->SetFillColor(38);
  fEfficiencyPerRunRestrictedRapidity31to26H->LabelsDeflate();
  fOutputList->Add(fEfficiencyPerRunRestrictedRapidity31to26H);

  fMCEfficiencyPerRunRestrictedRapidity31to26H = new TH1F("fMCEfficiencyPerRunRestrictedRapidity31to26H", "fMCEfficiencyPerRunRestrictedRapidity31to26H", 3, 0, 3);
  fMCEfficiencyPerRunRestrictedRapidity31to26H->SetStats(0);
  fMCEfficiencyPerRunRestrictedRapidity31to26H->SetFillColor(38);
  fMCEfficiencyPerRunRestrictedRapidity31to26H->LabelsDeflate();
  fOutputList->Add(fMCEfficiencyPerRunRestrictedRapidity31to26H);

  fEfficiencyPerRunWithRunTwoSettings = new TH1F("fEfficiencyPerRunWithRunTwoSettings", "fEfficiencyPerRunWithRunTwoSettings", 3, 0, 3);
  fEfficiencyPerRunWithRunTwoSettings->SetStats(0);
  fEfficiencyPerRunWithRunTwoSettings->SetFillColor(38);
  fEfficiencyPerRunWithRunTwoSettings->LabelsDeflate();
  fOutputList->Add(fEfficiencyPerRunWithRunTwoSettings);

  for( Int_t iLoop = 0; iLoop < 2; iLoop++ ) {
    fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop] = new TH1F(
      Form( "fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH_%d", iLoop),
      Form( "fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]->SetStats(0);
    fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]->SetFillColor(38);
    fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]);
  }

  for( Int_t iLoop = 0; iLoop < 3; iLoop++ ) {
    fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop] = new TH1F(
      Form( "fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH_%d", iLoop),
      Form( "fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]->SetStats(0);
    fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]->SetFillColor(38);
    fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]);
  }

  for( Int_t iLoop = 0; iLoop < 4; iLoop++ ) {
    fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop] = new TH1F(
      Form( "fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH_%d", iLoop),
      Form( "fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]->SetStats(0);
    fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]->SetFillColor(38);
    fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]);
  }

  for( Int_t iLoop = 0; iLoop < 5; iLoop++ ) {
    fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop] = new TH1F(
      Form( "fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH_%d", iLoop),
      Form( "fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]->SetStats(0);
    fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]->SetFillColor(38);
    fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]);
  }





  fMCEfficiencyPerRunWithRunTwoSettings = new TH1F("fMCEfficiencyPerRunWithRunTwoSettings", "fMCEfficiencyPerRunWithRunTwoSettings", 3, 0, 3);
  fMCEfficiencyPerRunWithRunTwoSettings->SetStats(0);
  fMCEfficiencyPerRunWithRunTwoSettings->SetFillColor(38);
  fMCEfficiencyPerRunWithRunTwoSettings->LabelsDeflate();
  fOutputList->Add(fMCEfficiencyPerRunWithRunTwoSettings);

  for( Int_t iLoop = 0; iLoop < 2; iLoop++ ) {
    fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop] = new TH1F(
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH_%d", iLoop),
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]->SetStats(0);
    fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]->SetFillColor(38);
    fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[iLoop]);
  }

  for( Int_t iLoop = 0; iLoop < 3; iLoop++ ) {
    fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop] = new TH1F(
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH_%d", iLoop),
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]->SetStats(0);
    fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]->SetFillColor(38);
    fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[iLoop]);
  }

  for( Int_t iLoop = 0; iLoop < 4; iLoop++ ) {
    fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop] = new TH1F(
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH_%d", iLoop),
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]->SetStats(0);
    fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]->SetFillColor(38);
    fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[iLoop]);
  }

  for( Int_t iLoop = 0; iLoop < 5; iLoop++ ) {
    fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop] = new TH1F(
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH_%d", iLoop),
      Form( "fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH_%d", iLoop),
      3, 0, 3
      );
    fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]->SetStats(0);
    fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]->SetFillColor(38);
    fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]->LabelsDeflate();
    fOutputList->Add(fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[iLoop]);
  }




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
  fRunNum    = fAOD->GetRunNumber();
  if(fMCEvent) {
    fRunNum    = fAOD->GetRunNumber();
    SetLuminosityCap();
    fCounterGeneratedLevel[ fRunNum - 240000 ] += 1;
    // cout << "fCounterGeneratedLevel[ " << (fRunNum - 240000) << " ] = " << fCounterGeneratedLevel[ fRunNum - 240000 ] << endl;
    // if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( (Int_t)fLumiPerRun * (Int_t)40000 ) ) {
    if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( fLumiPerRun * 100.0 ) ) {
          PostData(1, fOutputList);
          return;
    }
    ProcessMCParticles(fMCEvent);
    fCounterUPCevent += 1;
    fMCEfficiencyPerRunH->Fill( Form("%d", fRunNum) , 1 );
  }
  if( fCounterGeneratedLevel[ fRunNum - 240000 ] > ( fLumiPerRun * 100.0 ) ) {
        PostData(1, fOutputList);
        return;
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
                                          265697, /*265696,*/ 265694, 265691, 265607, 265596, 265594 };
  Int_t listOfGoodRunNumbersLHC16s[]  = { 267131, 267130, 267110, 267109, 267077, 267072, 267070, 267067, 267063, 267062,
                                          267022, 267020, 266998, 266997, 266994, 266993, 266988, 266944, 266943, 266942,
                                          266940, 266915, 266912, 266886, 266885, 266883, 266882, 266880, 266878, 266857,
                                          266807, 266805, 266800, 266776, 266775, 266708, 266706, 266703, 266702, 266676,
                                          266674, 266669, 266668, 266665, 266659, 266658, 266657, 266630, 266621, 266618,
                                          /*266615,*/ 266614, 266613, 266595, 266593, 266591, 266588, 266587, 266584, 266549,
                                          266543, 266539, 266534, 266533, 266525, 266523, 266522, 266520, 266518, 266516,
                                          266514, 266487, 266480, 266479, 266472, 266441, 266439, 295585 };
  Bool_t checkIfGoodRun = kFALSE;
  // cout << "OK4" << endl;
  // for( Int_t iRunLHC16r = 0; iRunLHC16r <  56; iRunLHC16r++){
  //   if( fRunNum == listOfGoodRunNumbersLHC16r[iRunLHC16r] ) checkIfGoodRun = kTRUE;
  // }
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
  // AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  // if(dataAD) {
  //       fCounterH->Fill(iSelectionCounter);
  //       iSelectionCounter++;
  //
  //       fADADecision = dataAD->GetADADecision();
  //       fADCDecision = dataAD->GetADCDecision();
  //
  //       // Reset event info
  //       fBBCFlagsAD = 0;
  //       fBGCFlagsAD = 0;
  //       fBBAFlagsAD = 0;
  //       fBGAFlagsAD = 0;
  //       for(Int_t i=0; i<16; i++) {
  //         // get array of fired pads
  //         fBBFlagAD[i] = dataAD->GetBBFlag(i);
  //         fBGFlagAD[i] = dataAD->GetBGFlag(i);
  //       }
  //
  //       for(Int_t i=0; i<4; i++) { // loop over pairs of pads
  //         if ( fBBFlagAD[i]   && fBBFlagAD[i+4]  ) fBBCFlagsAD++;
  //         if ( fBGFlagAD[i]   && fBGFlagAD[i+4]  ) fBGCFlagsAD++;
  //         if ( fBBFlagAD[i+8] && fBBFlagAD[i+12] ) fBBAFlagsAD++;
  //         if ( fBGFlagAD[i+8] && fBGFlagAD[i+12] ) fBGAFlagsAD++;
  //       }
  // }




  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  fCounterH->Fill(19);
  Int_t is_ADA_set = -9;
  Int_t is_ADC_set = -9;
  Double_t ADmultiplicities[16]   = { -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1 };
  Double_t ADmultiplicitiesTotal  = 0;
  Double_t ADAmultiplicitiesTotal = 0;
  Double_t ADCmultiplicitiesTotal = 0;

  Int_t ADAPastFutureBeamBeamFlags[8][21];
  Int_t ADCPastFutureBeamBeamFlags[8][21];


  Int_t ADAPastFutureBoolean = 0;
  Int_t ADCPastFutureBoolean = 0;

  if(dataAD) {
        fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;
        fCounterH->Fill(20);

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();
        fCounterH->Fill(21);


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



        is_ADA_set = IntBits( dataAD->GetTriggerBits() ).test(12);
        is_ADC_set = IntBits( dataAD->GetTriggerBits() ).test(13);
        // cout << "is_ADA_set = " << is_ADA_set << endl;
        // cout << "is_ADC_set = " << is_ADC_set << endl;
        // cout << "is_ADA_set = " << IntBits( dataAD->GetTriggerBits() ) << endl;
        // cout << "is_ADC_set = " << dataAD->GetTriggerBits() << endl;
        for( Int_t iChannel = 0; iChannel < 16; iChannel++ ){
          ADmultiplicities[iChannel] = dataAD->GetMultiplicity(iChannel);
          ADmultiplicitiesTotal     += dataAD->GetMultiplicity(iChannel);
          if ( iChannel < 8 ) {
            ADCmultiplicitiesTotal  += dataAD->GetMultiplicity(iChannel);
          } else {
            ADAmultiplicitiesTotal  += dataAD->GetMultiplicity(iChannel);
          }
        }


        for(   Int_t iChannel = 0; iChannel < 8; iChannel++ ){
          for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
            ADAPastFutureBeamBeamFlags[iChannel][iClock] = 0;
            ADCPastFutureBeamBeamFlags[iChannel][iClock] = 0;
          }
        }

        for(   Int_t iChannel = 0; iChannel < 8; iChannel++ ){
          for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
            ADAPastFutureBeamBeamFlags[iChannel][iClock] = dataAD->GetPFBBFlag(iChannel + 8, iClock);
            ADCPastFutureBeamBeamFlags[iChannel][iClock] = dataAD->GetPFBBFlag(iChannel, iClock);
          }
        }

        for(   Int_t iChannel = 0; iChannel < 8;  iChannel++ ){
          for( Int_t iClock   = 0; iClock   < 21; iClock++   ){
            if( dataAD->GetPFBBFlag(iChannel + 8, iClock) != 0 ) ADAPastFutureBoolean = 1;
            if( dataAD->GetPFBBFlag(iChannel, iClock)     != 0 ) ADCPastFutureBoolean = 1;
          }
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
  // if(fTracklets != 0) {
  //      PostData(1, fOutputList);
  //      return;
  // }
  /* - Maximum 2 V0C cells fired.
   */
  if( fV0TotalNCells > 2 ) {
       PostData(1, fOutputList);
       return;
  }





  //_______________________________
  /* -
   * - ADC multiplicity cut
   * -
   */
  // if( ADCmultiplicitiesTotal != 0 ) {
  //      PostData(1, fOutputList);
  //      return;
  // }

  //_______________________________
  /* -
   * - ADA multiplicity cut
   * -
   */
  // if( ADAmultiplicitiesTotal != 0 ) {
  //      PostData(1, fOutputList);
  //      return;
  // }



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


    /* -
     * - Compatibility with Run 1 analysis.
     * -
     */
    // if ( !( (track[nGoodMuons]->Eta() < -2.5) && (track[nGoodMuons]->Eta() > -3.7) ) ) {
    //   continue;
    // }


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
  if ( fRunNum < 290000 ) {
    ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( track[0]->Eta(), track[0]->Phi() );
    ((TH2F*) fOutputList->FindObject(Form( "fDeadZoneEtaVsPhiPerRunH_%d", fRunNum )) )->Fill( track[1]->Eta(), track[1]->Phi() );
  }

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
  fEtaDimuonH ->Fill( possibleJPsi.Rapidity() );
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

  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
      if ( (possibleJPsi.Rapidity() > -3.6) && (possibleJPsi.Rapidity() < -2.6) ) {
          fEfficiencyPerRunRestrictedRapidityH->Fill( Form("%d", fRunNum) , 1 );
          if ( possibleJPsi.Rapidity() < -3.1 ) {
            fEfficiencyPerRunRestrictedRapidity36to31H->Fill( Form("%d", fRunNum) , 1 );
          } else {
            fEfficiencyPerRunRestrictedRapidity31to26H->Fill( Form("%d", fRunNum) , 1 );
          }
          // cout << "OK1" << endl;
      }
  }



  /* -
   * - Efficiency with Run 2 settings.
   * -
   */
  if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ) {
      if ( (possibleJPsi.Rapidity() > -4.0) && (possibleJPsi.Rapidity() < -2.5) ) {
          fEfficiencyPerRunWithRunTwoSettings->Fill( Form("%d", fRunNum) , 1 );
          //
          //
          // Two Bins
          if ( possibleJPsi.Rapidity() <= -3.25 ) {
            fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
          } else {
            fEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
          }
          //
          //
          // Three Bins
          if (         possibleJPsi.Rapidity() <= -3.5 ) {
            fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -3.5) && (possibleJPsi.Rapidity() <= -3.0) ) {
            fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -3.0) && (possibleJPsi.Rapidity() <= -2.5) ) {
            fEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[2]->Fill( Form("%d", fRunNum) , 1 );
          }
          //
          //
          // Four Bins
          if (         possibleJPsi.Rapidity() <= -3.625 ) {
            fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -3.625) && (possibleJPsi.Rapidity() <= -3.250) ) {
            fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -3.250) && (possibleJPsi.Rapidity() <= -2.875) ) {
            fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[2]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -2.875) && (possibleJPsi.Rapidity() <= -2.500) ) {
            fEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[3]->Fill( Form("%d", fRunNum) , 1 );
          }
          //
          //
          // Five Bins
          if (         possibleJPsi.Rapidity() <= -3.7 ) {
            fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -3.7) && (possibleJPsi.Rapidity() <= -3.4) ) {
            fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -3.4) && (possibleJPsi.Rapidity() <= -3.1) ) {
            fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[2]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -3.1) && (possibleJPsi.Rapidity() <= -2.8) ) {
            fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[3]->Fill( Form("%d", fRunNum) , 1 );
          } else if ( (possibleJPsi.Rapidity() >  -2.8) && (possibleJPsi.Rapidity() <= -2.5) ) {
            fEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[4]->Fill( Form("%d", fRunNum) , 1 );
          }



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
  Bool_t is0VC5fired = fBBCFlags > 4;
  Bool_t is0VGAfired = fBGAFlags > 0;
  Bool_t is0UBAfired = fBBAFlagsAD > 0;
  Bool_t is0UBCfired = fBBCFlagsAD > 0;
  Bool_t is0UGCfired = fBGCFlagsAD > 0;
  // cout << "is0VBAfired = ( fBBAFlags   = " << fBBAFlags   << " ) > 0 => " << is0VBAfired << endl;
  // cout << "is0VBCfired = ( fBBCFlags   = " << fBBCFlags   << " ) > 0 => " << is0VBCfired << endl;
  // cout << "is0UBAfired = ( fBBAFlagsAD = " << fBBAFlagsAD << " ) > 0 => " << is0UBAfired << endl;
  // cout << "is0UBCfired = ( fBBCFlagsAD = " << fBBCFlagsAD << " ) > 0 => " << is0UBCfired << endl;
  // if (!is0VBAfired && !is0UBAfired && !is0UBCfired ) return kTRUE;
  if (!is0VBAfired && !is0VGAfired && !is0VC5fired && !is0UBCfired && !is0UGCfired ) return kTRUE;
  // if (!is0VBAfired && !is0UBAfired ) return kTRUE;
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
    if (  !mcParticle->IsPrimary()                                ) continue;
    if (   mcParticle->Charge()              == 0                 ) continue;
    // if ( !(mcParticle->Eta() < -2.5 && mcParticle->Eta() > -3.7 ) ) continue;
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
      if ( (possibleJPsiMC.Mag() > 2.85) && (possibleJPsiMC.Mag() < 3.35) ) {
          // cout << "OK2" << endl;
          if ( (possibleJPsiMC.Rapidity() > -3.6) && (possibleJPsiMC.Rapidity() < -2.6) ) {
              fMCEfficiencyPerRunRestrictedRapidityH->Fill( Form("%d", fRunNum) , 1 );
              if ( possibleJPsiMC.Rapidity() < -3.1 ) {
                fMCEfficiencyPerRunRestrictedRapidity36to31H->Fill( Form("%d", fRunNum) , 1 );
              } else {
                fMCEfficiencyPerRunRestrictedRapidity31to26H->Fill( Form("%d", fRunNum) , 1 );
              }
          }
          // cout << "OK3" << endl;
      }



      /* -
       * - Efficiency with Run 2 settings.
       * -
       */
      if ( (possibleJPsiMC.Mag() > 2.85) && (possibleJPsiMC.Mag() < 3.35) ) {
          if ( (possibleJPsiMC.Rapidity() > -4.0) && (possibleJPsiMC.Rapidity() < -2.5) ) {
              fMCEfficiencyPerRunWithRunTwoSettings->Fill( Form("%d", fRunNum) , 1 );
              //
              //
              // Two Bins
              if ( possibleJPsiMC.Rapidity() <= -3.25 ) {
                fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
              } else {
                fMCEfficiencyPerRunWithRunTwoSettingsTwoBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
              }
              //
              //
              // Three Bins
              if (         possibleJPsiMC.Rapidity() <= -3.5 ) {
                fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -3.5) && (possibleJPsiMC.Rapidity() <= -3.0) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -3.0) && (possibleJPsiMC.Rapidity() <= -2.5) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsThreeBinsRapidityH[2]->Fill( Form("%d", fRunNum) , 1 );
              }
              //
              //
              // Four Bins
              if (         possibleJPsiMC.Rapidity() <= -3.625 ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -3.625) && (possibleJPsiMC.Rapidity() <= -3.250) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -3.250) && (possibleJPsiMC.Rapidity() <= -2.875) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[2]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -2.875) && (possibleJPsiMC.Rapidity() <= -2.500) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFourBinsRapidityH[3]->Fill( Form("%d", fRunNum) , 1 );
              }
              //
              //
              // Five Bins
              if (         possibleJPsiMC.Rapidity() <= -3.7 ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[0]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -3.7) && (possibleJPsiMC.Rapidity() <= -3.4) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[1]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -3.4) && (possibleJPsiMC.Rapidity() <= -3.1) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[2]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -3.1) && (possibleJPsiMC.Rapidity() <= -2.8) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[3]->Fill( Form("%d", fRunNum) , 1 );
              } else if ( (possibleJPsiMC.Rapidity() >  -2.8) && (possibleJPsiMC.Rapidity() <= -2.5) ) {
                fMCEfficiencyPerRunWithRunTwoSettingsFiveBinsRapidityH[4]->Fill( Form("%d", fRunNum) , 1 );
              }



          }
      }


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
  // if      ( fRunNum == 267131 ) { fLumiPerRun = 359.819; }
  // else if ( fRunNum == 267130 ) { fLumiPerRun = 163.229; }
  // else if ( fRunNum == 267110 ) { fLumiPerRun = 348.615; }
  // else if ( fRunNum == 267109 ) { fLumiPerRun = 618.29;  }
  // else if ( fRunNum == 267077 ) { fLumiPerRun = 64.7591; }
  // else if ( fRunNum == 267072 ) { fLumiPerRun = 132.145; }
  // else if ( fRunNum == 267070 ) { fLumiPerRun = 49.2531; }
  // else if ( fRunNum == 267067 ) { fLumiPerRun = 142.775; }
  // else if ( fRunNum == 267063 ) { fLumiPerRun = 104.584; }
  // else if ( fRunNum == 267062 ) { fLumiPerRun = 82.6746; }
  // else if ( fRunNum == 267022 ) { fLumiPerRun = 113.37;  }
  // else if ( fRunNum == 267020 ) { fLumiPerRun = 543.514; }
  // else if ( fRunNum == 266998 ) { fLumiPerRun = 26.5804; }
  // else if ( fRunNum == 266997 ) { fLumiPerRun = 26.7698; }
  // else if ( fRunNum == 266994 ) { fLumiPerRun = 68.9501; }
  // else if ( fRunNum == 266993 ) { fLumiPerRun = 30.8271; }
  // else if ( fRunNum == 266988 ) { fLumiPerRun = 374.406; }
  // else if ( fRunNum == 266944 ) { fLumiPerRun = 208.054; }
  // else if ( fRunNum == 266943 ) { fLumiPerRun = 200.326; }
  // else if ( fRunNum == 266942 ) { fLumiPerRun = 83.6478; }
  // else if ( fRunNum == 266940 ) { fLumiPerRun = 104.177; }
  // else if ( fRunNum == 266915 ) { fLumiPerRun = 70.9138; }
  // else if ( fRunNum == 266912 ) { fLumiPerRun = 198.575; }
  // else if ( fRunNum == 266886 ) { fLumiPerRun = 204.031; }
  // else if ( fRunNum == 266885 ) { fLumiPerRun = 136.839; }
  // else if ( fRunNum == 266883 ) { fLumiPerRun = 160.311; }
  // else if ( fRunNum == 266882 ) { fLumiPerRun = 91.8889; }
  // else if ( fRunNum == 266880 ) { fLumiPerRun = 41.6682; }
  // else if ( fRunNum == 266878 ) { fLumiPerRun = 530.376; }
  // else if ( fRunNum == 266857 ) { fLumiPerRun = 58.3583; }
  // else if ( fRunNum == 266807 ) { fLumiPerRun = 70.9539; }
  // else if ( fRunNum == 266805 ) { fLumiPerRun = 167.801; }
  // else if ( fRunNum == 266800 ) { fLumiPerRun = 327.175; }
  // else if ( fRunNum == 266776 ) { fLumiPerRun = 420.595; }
  // else if ( fRunNum == 266775 ) { fLumiPerRun = 709.504; }
  // else if ( fRunNum == 266708 ) { fLumiPerRun = 91.9383; }
  // else if ( fRunNum == 266706 ) { fLumiPerRun = 235.296; }
  // else if ( fRunNum == 266703 ) { fLumiPerRun = 61.2925; }
  // else if ( fRunNum == 266702 ) { fLumiPerRun = 106.221; }
  // else if ( fRunNum == 266676 ) { fLumiPerRun = 32.0462; }
  // else if ( fRunNum == 266674 ) { fLumiPerRun = 52.7185; }
  // else if ( fRunNum == 266669 ) { fLumiPerRun = 254.209; }
  // else if ( fRunNum == 266668 ) { fLumiPerRun = 68.8008; }
  // else if ( fRunNum == 266665 ) { fLumiPerRun = 131.907; }
  // else if ( fRunNum == 266659 ) { fLumiPerRun = 190.677; }
  // else if ( fRunNum == 266658 ) { fLumiPerRun = 58.6213; }
  // else if ( fRunNum == 266657 ) { fLumiPerRun = 593.284; }
  // else if ( fRunNum == 266630 ) { fLumiPerRun = 49.7633; }
  // else if ( fRunNum == 266621 ) { fLumiPerRun = 127.857; }
  // else if ( fRunNum == 266618 ) { fLumiPerRun = 277;     }
  // else if ( fRunNum == 266614 ) { fLumiPerRun = 457.478; }
  // else if ( fRunNum == 266613 ) { fLumiPerRun = 386.435; }
  // else if ( fRunNum == 266595 ) { fLumiPerRun = 329.874; }
  // else if ( fRunNum == 266593 ) { fLumiPerRun = 201.278; }
  // else if ( fRunNum == 266591 ) { fLumiPerRun = 69.7556; }
  // else if ( fRunNum == 266588 ) { fLumiPerRun = 70.7593; }
  // else if ( fRunNum == 266587 ) { fLumiPerRun = 59.1278; }
  // else if ( fRunNum == 266584 ) { fLumiPerRun = 121.033; }
  // else if ( fRunNum == 266549 ) { fLumiPerRun = 56.9813; }
  // else if ( fRunNum == 266543 ) { fLumiPerRun = 190.57;  }
  // else if ( fRunNum == 266539 ) { fLumiPerRun = 170.217; }
  // else if ( fRunNum == 266534 ) { fLumiPerRun = 80.2625; }
  // else if ( fRunNum == 266533 ) { fLumiPerRun = 98.5803; }
  // else if ( fRunNum == 266525 ) { fLumiPerRun = 84.4293; }
  // else if ( fRunNum == 266523 ) { fLumiPerRun = 55.9541; }
  // else if ( fRunNum == 266522 ) { fLumiPerRun = 135.808; }
  // else if ( fRunNum == 266520 ) { fLumiPerRun = 81.0349; }
  // else if ( fRunNum == 266518 ) { fLumiPerRun = 371.843; }
  // else if ( fRunNum == 266516 ) { fLumiPerRun = 24.7453; }
  // else if ( fRunNum == 266514 ) { fLumiPerRun = 25.1691; }
  // else if ( fRunNum == 266487 ) { fLumiPerRun = 36.0792; }
  // else if ( fRunNum == 266480 ) { fLumiPerRun = 940.67;  }
  // else if ( fRunNum == 266479 ) { fLumiPerRun = 214.998; }
  // else if ( fRunNum == 266472 ) { fLumiPerRun = 105.209; }
  // else if ( fRunNum == 266441 ) { fLumiPerRun = 177.071; }
  // else if ( fRunNum == 266439 ) { fLumiPerRun = 41.7182; }
  // else if ( fRunNum == 266318 ) { fLumiPerRun = 66.2315; }
  // else if ( fRunNum == 266316 ) { fLumiPerRun = 9.58106; }
  // else if ( fRunNum == 266312 ) { fLumiPerRun = 157.511; }
  // else if ( fRunNum == 266305 ) { fLumiPerRun = 198.184; }
  // else if ( fRunNum == 266304 ) { fLumiPerRun = 85.8848; }
  // else if ( fRunNum == 266300 ) { fLumiPerRun = 129.924; }
  // else if ( fRunNum == 266299 ) { fLumiPerRun = 133.485; }
  // else if ( fRunNum == 266296 ) { fLumiPerRun = 116.474; }
  // else if ( fRunNum == 266235 ) { fLumiPerRun = 414.245; }
  // else if ( fRunNum == 266234 ) { fLumiPerRun = 199.209; }
  // else if ( fRunNum == 266208 ) { fLumiPerRun = 149.731; }
  // else if ( fRunNum == 266197 ) { fLumiPerRun = 118.298; }
  // else if ( fRunNum == 266196 ) { fLumiPerRun = 68.0719; }
  // else if ( fRunNum == 266193 ) { fLumiPerRun = 88.4665; }
  // else if ( fRunNum == 266190 ) { fLumiPerRun = 79.6631; }
  // else if ( fRunNum == 266189 ) { fLumiPerRun = 32.3559; }
  // else if ( fRunNum == 266187 ) { fLumiPerRun = 127.149; }
  // else if ( fRunNum == 266117 ) { fLumiPerRun = 153.453; }
  // else if ( fRunNum == 266086 ) { fLumiPerRun = 116.072; }
  // else if ( fRunNum == 266085 ) { fLumiPerRun = 63.4458; }
  // else if ( fRunNum == 266084 ) { fLumiPerRun = 22.4614; }
  // else if ( fRunNum == 266081 ) { fLumiPerRun = 47.3174; }
  // else if ( fRunNum == 266076 ) { fLumiPerRun = 195.23;  }
  // else if ( fRunNum == 266074 ) { fLumiPerRun = 230.263; }
  // else if ( fRunNum == 266034 ) { fLumiPerRun = 86.9949; }
  // else if ( fRunNum == 266025 ) { fLumiPerRun = 658.516; }
  // else if ( fRunNum == 266023 ) { fLumiPerRun = 133.836; }
  // else if ( fRunNum == 266022 ) { fLumiPerRun = 340.669; }
  // else if ( fRunNum == 265841 ) { fLumiPerRun = 210.894; }
  // else if ( fRunNum == 265840 ) { fLumiPerRun = 3.65278; }
  // else if ( fRunNum == 265797 ) { fLumiPerRun = 41.2853; }
  // else if ( fRunNum == 265795 ) { fLumiPerRun = 65.9944; }
  // else if ( fRunNum == 265792 ) { fLumiPerRun = 34.2686; }
  // else if ( fRunNum == 265789 ) { fLumiPerRun = 168.22;  }
  // else if ( fRunNum == 265788 ) { fLumiPerRun = 145.194; }
  // else if ( fRunNum == 265787 ) { fLumiPerRun = 251.743; }
  // else if ( fRunNum == 265785 ) { fLumiPerRun = 316.091; }
  // else if ( fRunNum == 265756 ) { fLumiPerRun = 51.7427; }
  // else if ( fRunNum == 265754 ) { fLumiPerRun = 96.0737; }
  // else if ( fRunNum == 265746 ) { fLumiPerRun = 366.01;  }
  // else if ( fRunNum == 265744 ) { fLumiPerRun = 168.605; }
  // else if ( fRunNum == 265742 ) { fLumiPerRun = 184.746; }
  // else if ( fRunNum == 265741 ) { fLumiPerRun = 80.0317; }
  // else if ( fRunNum == 265740 ) { fLumiPerRun = 72.2736; }
  // else if ( fRunNum == 265714 ) { fLumiPerRun = 46.912;  }
  // else if ( fRunNum == 265713 ) { fLumiPerRun = 41.2605; }
  // else if ( fRunNum == 265709 ) { fLumiPerRun = 48.43;   }
  // else if ( fRunNum == 265701 ) { fLumiPerRun = 124.259; }
  // else if ( fRunNum == 265700 ) { fLumiPerRun = 33.3219; }
  // else if ( fRunNum == 265698 ) { fLumiPerRun = 140.203; }
  // else if ( fRunNum == 265697 ) { fLumiPerRun = 19.0271; }
  // else if ( fRunNum == 265694 ) { fLumiPerRun = 577.183; }
  // else if ( fRunNum == 265691 ) { fLumiPerRun = 351.54;  }
  // else if ( fRunNum == 265607 ) { fLumiPerRun = 0.647854;}
  // else if ( fRunNum == 265596 ) { fLumiPerRun = 4.23135; }

  //_______________________________
  /* -
   * - Asking for AD decisions.
   * -
   */
  if      ( fRunNum == 267131 ) { fLumiPerRun = 299.488; }
  else if ( fRunNum == 267130 ) { fLumiPerRun = 133.813; }
  else if ( fRunNum == 267110 ) { fLumiPerRun = 293.588; }
  else if ( fRunNum == 267109 ) { fLumiPerRun = 506.868; }
  else if ( fRunNum == 267077 ) { fLumiPerRun = 52.7582; }
  else if ( fRunNum == 267072 ) { fLumiPerRun = 108.475; }
  else if ( fRunNum == 267070 ) { fLumiPerRun = 40.3772; }
  else if ( fRunNum == 267067 ) { fLumiPerRun = 117.626; }
  else if ( fRunNum == 267063 ) { fLumiPerRun = 122.597; }
  else if ( fRunNum == 267062 ) { fLumiPerRun = 67.7758; }
  else if ( fRunNum == 267022 ) { fLumiPerRun = 95.3051; }
  else if ( fRunNum == 267020 ) { fLumiPerRun = 446.61; }
  else if ( fRunNum == 266998 ) { fLumiPerRun = 22.3174; }
  else if ( fRunNum == 266997 ) { fLumiPerRun = 22.4019; }
  else if ( fRunNum == 266994 ) { fLumiPerRun = 58.4044; }
  else if ( fRunNum == 266993 ) { fLumiPerRun = 26.0976; }
  else if ( fRunNum == 266988 ) { fLumiPerRun = 307.71; }
  else if ( fRunNum == 266944 ) { fLumiPerRun = 174.388; }
  else if ( fRunNum == 266943 ) { fLumiPerRun = 164.8; }
  else if ( fRunNum == 266942 ) { fLumiPerRun = 68.5736; }
  else if ( fRunNum == 266940 ) { fLumiPerRun = 110.686; }
  else if ( fRunNum == 266915 ) { fLumiPerRun = 58.1344; }
  else if ( fRunNum == 266912 ) { fLumiPerRun = 163.755; }
  else if ( fRunNum == 266886 ) { fLumiPerRun = 167.216; }
  else if ( fRunNum == 266885 ) { fLumiPerRun = 111.504; }
  else if ( fRunNum == 266883 ) { fLumiPerRun = 132.063; }
  else if ( fRunNum == 266882 ) { fLumiPerRun = 76.5574; }
  else if ( fRunNum == 266880 ) { fLumiPerRun = 34.1592; }
  else if ( fRunNum == 266878 ) { fLumiPerRun = 444.15; }
  else if ( fRunNum == 266857 ) { fLumiPerRun = 47.3076; }
  else if ( fRunNum == 266807 ) { fLumiPerRun = 58.4043; }
  else if ( fRunNum == 266805 ) { fLumiPerRun = 137.953; }
  else if ( fRunNum == 266800 ) { fLumiPerRun = 277.81; }
  else if ( fRunNum == 266776 ) { fLumiPerRun = 347.019; }
  else if ( fRunNum == 266775 ) { fLumiPerRun = 587.711; }
  else if ( fRunNum == 266708 ) { fLumiPerRun = 77.0559; }
  else if ( fRunNum == 266706 ) { fLumiPerRun = 255.623; }
  else if ( fRunNum == 266703 ) { fLumiPerRun = 47.0861; }
  else if ( fRunNum == 266702 ) { fLumiPerRun = 121.264; }
  else if ( fRunNum == 266676 ) { fLumiPerRun = 25.0438; }
  else if ( fRunNum == 266674 ) { fLumiPerRun = 43.2181; }
  else if ( fRunNum == 266669 ) { fLumiPerRun = 202.677; }
  else if ( fRunNum == 266668 ) { fLumiPerRun = 151.058; }
  else if ( fRunNum == 266665 ) { fLumiPerRun = 102.909; }
  else if ( fRunNum == 266659 ) { fLumiPerRun = 147.023; }
  else if ( fRunNum == 266658 ) { fLumiPerRun = 44.1322; }
  else if ( fRunNum == 266657 ) { fLumiPerRun = 451.944; }
  else if ( fRunNum == 266630 ) { fLumiPerRun = 40.7954; }
  else if ( fRunNum == 266621 ) { fLumiPerRun = 138.065; }
  else if ( fRunNum == 266618 ) { fLumiPerRun = 228.025; }
  else if ( fRunNum == 266614 ) { fLumiPerRun = 372.09; }
  else if ( fRunNum == 266613 ) { fLumiPerRun = 310.359; }
  else if ( fRunNum == 266595 ) { fLumiPerRun = 253.44; }
  else if ( fRunNum == 266593 ) { fLumiPerRun = 158.782; }
  else if ( fRunNum == 266591 ) { fLumiPerRun = 57.1858; }
  else if ( fRunNum == 266588 ) { fLumiPerRun = 60.8278; }
  else if ( fRunNum == 266587 ) { fLumiPerRun = 123.611; }
  else if ( fRunNum == 266584 ) { fLumiPerRun = 169.477; }
  else if ( fRunNum == 266549 ) { fLumiPerRun = 45.4535; }
  else if ( fRunNum == 266543 ) { fLumiPerRun = 151.764; }
  else if ( fRunNum == 266539 ) { fLumiPerRun = 139.022; }
  else if ( fRunNum == 266534 ) { fLumiPerRun = 65.8065; }
  else if ( fRunNum == 266533 ) { fLumiPerRun = 80.7648; }
  else if ( fRunNum == 266525 ) { fLumiPerRun = 67.4575; }
  else if ( fRunNum == 266523 ) { fLumiPerRun = 45.6672; }
  else if ( fRunNum == 266522 ) { fLumiPerRun = 111.334; }
  else if ( fRunNum == 266520 ) { fLumiPerRun = 59.4779; }
  else if ( fRunNum == 266518 ) { fLumiPerRun = 304.703; }
  else if ( fRunNum == 266516 ) { fLumiPerRun = 20.191; }
  else if ( fRunNum == 266514 ) { fLumiPerRun = 164.321; }
  else if ( fRunNum == 266487 ) { fLumiPerRun = 28.6638; }
  else if ( fRunNum == 266480 ) { fLumiPerRun = 738.723; }
  else if ( fRunNum == 266479 ) { fLumiPerRun = 164.735; }
  else if ( fRunNum == 266472 ) { fLumiPerRun = 86.2489; }
  else if ( fRunNum == 266441 ) { fLumiPerRun = 143.532; }
  else if ( fRunNum == 266439 ) { fLumiPerRun = 33.9646; }
  //_______________________________



  //_______________________________
  /* -
   * - pPb analysis
   * -
   */
  else if ( fRunNum == 266318 ) { fLumiPerRun = 66.2315; }
  else if ( fRunNum == 266316 ) { fLumiPerRun = 9.58106; }
  else if ( fRunNum == 266312 ) { fLumiPerRun = 157.511; }
  else if ( fRunNum == 266305 ) { fLumiPerRun = 198.184; }
  else if ( fRunNum == 266304 ) { fLumiPerRun = 85.8848; }
  else if ( fRunNum == 266300 ) { fLumiPerRun = 129.924; }
  else if ( fRunNum == 266299 ) { fLumiPerRun = 133.485; }
  else if ( fRunNum == 266296 ) { fLumiPerRun = 116.474; }
  else if ( fRunNum == 266235 ) { fLumiPerRun = 414.245; }
  else if ( fRunNum == 266234 ) { fLumiPerRun = 199.209; }
  else if ( fRunNum == 266208 ) { fLumiPerRun = 149.731; }
  else if ( fRunNum == 266197 ) { fLumiPerRun = 118.298; }
  else if ( fRunNum == 266196 ) { fLumiPerRun = 68.0719; }
  else if ( fRunNum == 266193 ) { fLumiPerRun = 88.4665; }
  else if ( fRunNum == 266190 ) { fLumiPerRun = 79.6631; }
  else if ( fRunNum == 266189 ) { fLumiPerRun = 32.3559; }
  else if ( fRunNum == 266187 ) { fLumiPerRun = 127.149; }
  else if ( fRunNum == 266117 ) { fLumiPerRun = 153.453; }
  else if ( fRunNum == 266086 ) { fLumiPerRun = 116.072; }
  else if ( fRunNum == 266085 ) { fLumiPerRun = 63.4458; }
  else if ( fRunNum == 266084 ) { fLumiPerRun = 22.4614; }
  else if ( fRunNum == 266081 ) { fLumiPerRun = 47.3174; }
  else if ( fRunNum == 266076 ) { fLumiPerRun = 195.23; }
  else if ( fRunNum == 266074 ) { fLumiPerRun = 230.263; }
  else if ( fRunNum == 266034 ) { fLumiPerRun = 86.9949; }
  else if ( fRunNum == 266025 ) { fLumiPerRun = 658.516; }
  else if ( fRunNum == 266023 ) { fLumiPerRun = 133.836; }
  else if ( fRunNum == 266022 ) { fLumiPerRun = 340.669; }
  else if ( fRunNum == 265841 ) { fLumiPerRun = 210.894; }
  else if ( fRunNum == 265840 ) { fLumiPerRun = 3.65278; }
  else if ( fRunNum == 265797 ) { fLumiPerRun = 41.2853; }
  else if ( fRunNum == 265795 ) { fLumiPerRun = 65.9944; }
  else if ( fRunNum == 265792 ) { fLumiPerRun = 34.2686; }
  else if ( fRunNum == 265789 ) { fLumiPerRun = 168.22; }
  else if ( fRunNum == 265788 ) { fLumiPerRun = 145.194; }
  else if ( fRunNum == 265787 ) { fLumiPerRun = 251.743; }
  else if ( fRunNum == 265785 ) { fLumiPerRun = 316.091; }
  else if ( fRunNum == 265756 ) { fLumiPerRun = 51.7427; }
  else if ( fRunNum == 265754 ) { fLumiPerRun = 96.0737; }
  else if ( fRunNum == 265746 ) { fLumiPerRun = 366.01; }
  else if ( fRunNum == 265744 ) { fLumiPerRun = 168.605; }
  else if ( fRunNum == 265742 ) { fLumiPerRun = 184.746; }
  else if ( fRunNum == 265741 ) { fLumiPerRun = 80.0317; }
  else if ( fRunNum == 265740 ) { fLumiPerRun = 72.2736; }
  else if ( fRunNum == 265714 ) { fLumiPerRun = 46.912; }
  else if ( fRunNum == 265713 ) { fLumiPerRun = 41.2605; }
  else if ( fRunNum == 265709 ) { fLumiPerRun = 48.43; }
  else if ( fRunNum == 265701 ) { fLumiPerRun = 124.259; }
  else if ( fRunNum == 265700 ) { fLumiPerRun = 33.3219; }
  else if ( fRunNum == 265698 ) { fLumiPerRun = 140.203; }
  else if ( fRunNum == 265697 ) { fLumiPerRun = 19.0271; }
  else if ( fRunNum == 265694 ) { fLumiPerRun = 577.183; }
  else if ( fRunNum == 265691 ) { fLumiPerRun = 351.54; }
  else if ( fRunNum == 265607 ) { fLumiPerRun = 0.647854; }
  else if ( fRunNum == 265596 ) { fLumiPerRun = 4.23135; }
  //_______________________________



  else                          { fLumiPerRun = 1.00;    }

}
