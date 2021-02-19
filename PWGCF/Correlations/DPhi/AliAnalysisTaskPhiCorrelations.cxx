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

/* $Id:$ */

#include <TROOT.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRandom.h>
#include <TParameter.h>

#include "AliAnalysisTaskPhiCorrelations.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliUEHistograms.h"
#include "AliUEHist.h"

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliVParticle.h"
#include "AliCFContainer.h"
#include "AliEventplane.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliStack.h"
#include "AliAODMCHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenHepMCEventHeader.h"

#include "AliEventPoolManager.h"
#include "AliBasicParticle.h"
#include "AliVHeader.h"

#include "AliESDZDC.h"
#include "AliESDtrackCuts.h"

#include "AliEmcalJet.h"

#include "AliMultSelection.h"

#include "AliHelperPID.h"
#include "AliAnalysisUtils.h"
#include "TMap.h"

////////////////////////////////////////////////////////////////////////
//
// Analysis class for azimuthal correlation studies
// Based on the UE task from Sara Vallero and Jan Fiete
//
// This class needs input AODs.
// The output is a list of analysis-specific containers.
//
// The AOD can be either connected to the InputEventHandler
// for a chain of AOD files
// or
// to the OutputEventHandler
// for a chain of ESD files,
// in this case the class should be in the train after the jet-finder
//
//    Authors:
//    Jan Fiete Grosse-Oetringhaus
//
////////////////////////////////////////////////////////////////////////


ClassImp( AliAnalysisTaskPhiCorrelations )

//____________________________________________________________________
AliAnalysisTaskPhiCorrelations::AliAnalysisTaskPhiCorrelations(const char* name):
AliAnalysisTask(name,""),
// general configuration
fDebug(0),
fMode(0),
fReduceMemoryFootprint(kFALSE),
fFillMixed(kTRUE),
fMixingTracks(50000),
fTwoTrackEfficiencyStudy(kFALSE),
fTwoTrackEfficiencyCut(0),
fTwoTrackCutMinRadius(0.8),
fUseVtxAxis(kFALSE),
fCourseCentralityBinning(kFALSE),
fSkipTrigger(kFALSE),
fInjectedSignals(kFALSE),
fRandomizeReactionPlane(kFALSE),
fV0CL1PileUp(0),
fESDTPCTrackPileUp(0),
fTPCITSTOFPileUp(0),
fHelperPID(0x0),
fAnalysisUtils(0x0),
fMap(0x0),
// pointers to UE classes
fAnalyseUE(0x0),
fHistos(0x0),
fHistosMixed(0),
fEfficiencyCorrectionTriggers(0),
fEfficiencyCorrectionAssociated(0),
fCentralityWeights(0),
fCentralityMCGen_V0M(0),
fCentralityMCGen_CL1(0),
// handlers and events
fAOD(0x0),
fESD(0x0),
fArrayMC(0x0),
fInputHandler(0x0),
fMcEvent(0x0),
fMcHandler(0x0),
fPoolMgr(0x0),
// histogram settings
fListOfHistos(0x0),
// event QA
fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default)
fZVertex(7.),
fAcceptOnlyMuEvents(kFALSE),
fCentralityMethod("V0M"),
// track cuts
fTrackEtaCut(0.8),
fTrackEtaCutMin(-1.),
fTrackPhiCutEvPlMin(0.),
fTrackPhiCutEvPlMax(0.),
fOnlyOneEtaSide(0),
fOnlyOneAssocEtaSide(0),
fPtMin(0.5),
fDCAXYCut(0),
fSharedClusterCut(-1),
fCrossedRowsCut(-1),
fFoundFractionCut(-1),
fFilterBit(0xFF),
fTrackStatus(0),
fSelectBit(AliVEvent::kMB|AliVEvent::kUserDefined),
fUseChargeHadrons(kFALSE),
fParticleSpeciesTrigger(-1),
fParticleSpeciesAssociated(-1),
fCheckMotherPDG(kTRUE),
fTrackletDphiCut(9999999.),
fSelectCharge(0),
fTriggerSelectCharge(0),
fAssociatedSelectCharge(0),
fTriggerRestrictEta(-1),
fEtaOrdering(kFALSE),
fCutConversionsV(-1),
fCutResonancesV(-1),
fCutOnCustomMass(-1),
fCutOnCustomFirst(-1),
fCutOnCustomSecond(-1),
fCutOnCustomV(-1),
fCutOnPhi(kFALSE),
fCutOnPhiV(-1),
fCutOnRho(kFALSE),
fCutOnRhoV(-1),
fCutOnLambdaV(-1),
fCutOnK0sV(-1),
fRejectResonanceDaughters(-1),
fFillOnlyStep0(kFALSE),
fSkipStep6(kFALSE),
fSkipStep9(kTRUE),
fRejectCentralityOutliers(kFALSE),
fRejectZeroTrackEvents(kFALSE),
fRemoveWeakDecays(kFALSE),
fRemoveDuplicates(kFALSE),
fSkipFastCluster(kFALSE),
fWeightPerEvent(kFALSE),
fCustomBinning(),
fPtOrder(kTRUE),
fTriggersFromDetector(0),
fAssociatedFromDetector(0),
fUseUncheckedCentrality(kFALSE),
fCheckCertainSpecies(-1),
fRemoveWeakDecaysInMC(kFALSE),
fFillYieldRapidity(kFALSE),
fFillCorrelationsRapidity(kFALSE),
fUseDoublePrecision(kFALSE),
fUseNewCentralityFramework(kFALSE),
fFillpT(kFALSE),
fJetBranchName("clustersAOD_ANTIKT04_B1_Filter00768_Cut00150_Skip00"),
fTrackEtaMax(.9),
fJetEtaMax(.9),
fJetPtMin(5.),
fJetConstMin(0),
fExclusionRadius(-1.),
fCustomParticlesA(""),
fCustomParticlesB(""),
fEventPoolOutputList(),
fUsePtBinnedEventPool(0),
fCheckEventNumberInMixedEvent(kFALSE)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
}

AliAnalysisTaskPhiCorrelations::~AliAnalysisTaskPhiCorrelations()
{
  // destructor

  if (fListOfHistos  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
    delete fListOfHistos;
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::ConnectInputData(Option_t* /*option*/)
{

  // Connect the input data
  if (fDebug > 1) AliInfo("ConnectInputData() ");

  // Since AODs can either be connected to the InputEventHandler
  // or to the OutputEventHandler ( the AOD is created by a previus task in the train )
  // we need to get the pointer to the AODEvent correctly.

  // Delta AODs are also accepted.

  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

  if ( handler && handler->InheritsFrom("AliAODInputHandler") )
  { // input AOD
    fAOD = ((AliAODInputHandler*)handler)->GetEvent();
    if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODInputHandler");
  }
  else
  {  //output AOD
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if (handler && handler->InheritsFrom("AliAODHandler") )
    {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODHandler");
    }
    else
    {  // no AOD
      AliWarning("I can't get any AOD Event Handler");
    }
  }

  if (handler && handler->InheritsFrom("AliESDInputHandler") )
  { // input ESD
    // pointer received per event in ::Exec
    if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliESDInputHandler");
  }

  // Initialize common pointers
  Initialize();
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::CreateOutputObjects()
{
  // Create the output container

  if (fDebug > 1) AliInfo("CreateOutputObjects()");

  // Initialize class with main algorithms, event and track selection.
  fAnalyseUE = new AliAnalyseLeadingTrackUE();
  fAnalyseUE->SetParticleSelectionCriteria(fFilterBit, fUseChargeHadrons, fTrackEtaCut, fTrackEtaCutMin, fPtMin);
  fAnalyseUE->SetDCAXYCut(fDCAXYCut);
  fAnalyseUE->SetSharedClusterCut(fSharedClusterCut);
  fAnalyseUE->SetCrossedRowsCut(fCrossedRowsCut);
  fAnalyseUE->SetFoundFractionCut(fFoundFractionCut);
  fAnalyseUE->SetTrackStatus(fTrackStatus);
  fAnalyseUE->SetCheckMotherPDG(fCheckMotherPDG);
  fAnalyseUE->SetDebug(fDebug);
  fAnalyseUE->DefineESDCuts(fFilterBit);
  fAnalyseUE->SetEventSelection(fSelectBit);
  fAnalyseUE->SetHelperPID(fHelperPID);
  if (fTrackPhiCutEvPlMax > 0.0001)
    fAnalyseUE->SetParticlePhiCutEventPlane(fTrackPhiCutEvPlMin,fTrackPhiCutEvPlMax);
  if ((fParticleSpeciesTrigger != -1 || fParticleSpeciesAssociated != -1) && !fHelperPID)
    AliFatal("HelperPID object should be set in the steering macro");

  // Initialize output list of containers
  if (fListOfHistos != NULL) {
    delete fListOfHistos;
    fListOfHistos = NULL;
  }
  if (!fListOfHistos) {
    fListOfHistos = new TList();
    fListOfHistos->SetOwner(kTRUE);
  }

  // Initialize class to handle histograms
  TString histType = "4R";
  if (fUseVtxAxis == 1)
    histType = "5R";
  else if (fUseVtxAxis == 2)
    histType = "6R";
  if (fCourseCentralityBinning)
    histType += "C";
  if (fUseDoublePrecision)
    histType += "D";
  fHistos = new AliUEHistograms("AliUEHistogramsSame", histType, fCustomBinning);
  fHistosMixed = new AliUEHistograms("AliUEHistogramsMixed", histType, fCustomBinning);

  // On demand, check event number before correlating tracks in mixed events
  // To avoid same event contributions in mixed event when importing event pool
  fHistosMixed->SetCheckEventNumberInCorrelation(fCheckEventNumberInMixedEvent);

  fHistos->SetSelectCharge(fSelectCharge);
  fHistosMixed->SetSelectCharge(fSelectCharge);

  fHistos->SetSelectTriggerCharge(fTriggerSelectCharge);
  fHistosMixed->SetSelectTriggerCharge(fTriggerSelectCharge);

  fHistos->SetSelectAssociatedCharge(fAssociatedSelectCharge);
  fHistosMixed->SetSelectAssociatedCharge(fAssociatedSelectCharge);

  fHistos->SetTriggerRestrictEta(fTriggerRestrictEta);
  fHistosMixed->SetTriggerRestrictEta(fTriggerRestrictEta);

  fHistos->SetOnlyOneEtaSide(fOnlyOneEtaSide);
  fHistosMixed->SetOnlyOneEtaSide(fOnlyOneEtaSide);

  fHistos->SetOnlyOneAssocEtaSide(fOnlyOneAssocEtaSide);
  fHistosMixed->SetOnlyOneAssocEtaSide(fOnlyOneAssocEtaSide);

  fHistos->SetEtaOrdering(fEtaOrdering);
  fHistosMixed->SetEtaOrdering(fEtaOrdering);

  fHistos->SetPairCuts(fCutConversionsV, fCutResonancesV);
  fHistosMixed->SetPairCuts(fCutConversionsV, fCutResonancesV);

  fHistos->SetCustomCut(fCutOnCustomMass, fCutOnCustomFirst, fCutOnCustomSecond, fCutOnCustomV);
  fHistosMixed->SetCustomCut(fCutOnCustomMass, fCutOnCustomFirst, fCutOnCustomSecond, fCutOnCustomV);

  fHistos->SetCutOnPhi(fCutOnPhi);
  fHistosMixed->SetCutOnPhi(fCutOnPhi);

  fHistos->SetCutOnPhi(fCutOnPhiV);
  fHistosMixed->SetCutOnPhi(fCutOnPhiV);

  fHistos->SetCutOnRho(fCutOnRho);
  fHistosMixed->SetCutOnRho(fCutOnRho);

  fHistos->SetCutOnRho(fCutOnRhoV);
  fHistosMixed->SetCutOnRho(fCutOnRhoV);

  fHistos->SetCutOnK0s(fCutOnK0sV);
  fHistosMixed->SetCutOnK0s(fCutOnK0sV);

  fHistos->SetCutOnLambda(fCutOnLambdaV);
  fHistosMixed->SetCutOnLambda(fCutOnLambdaV);

  fHistos->SetRejectResonanceDaughters(fRejectResonanceDaughters);
  fHistosMixed->SetRejectResonanceDaughters(fRejectResonanceDaughters);

  fHistos->SetTrackEtaCut(fTrackEtaCut);
  fHistosMixed->SetTrackEtaCut(fTrackEtaCut);

  fHistos->SetWeightPerEvent(fWeightPerEvent);
  fHistosMixed->SetWeightPerEvent(fWeightPerEvent);

  fHistos->SetPtOrder(fPtOrder);
  fHistosMixed->SetPtOrder(fPtOrder);

  fHistos->SetTwoTrackCutMinRadius(fTwoTrackCutMinRadius);
  fHistosMixed->SetTwoTrackCutMinRadius(fTwoTrackCutMinRadius);

  if (fEfficiencyCorrectionTriggers) {
    fHistos->SetEfficiencyCorrectionTriggers(fEfficiencyCorrectionTriggers);
    fHistosMixed->SetEfficiencyCorrectionTriggers((THnF*) fEfficiencyCorrectionTriggers->Clone());
  }
  if (fEfficiencyCorrectionAssociated) {
    fHistos->SetEfficiencyCorrectionAssociated(fEfficiencyCorrectionAssociated);
    fHistosMixed->SetEfficiencyCorrectionAssociated((THnF*) fEfficiencyCorrectionAssociated->Clone());
  }

  // add histograms to list
  fListOfHistos->Add(fHistos);
  fListOfHistos->Add(fHistosMixed);
  // add HelperPID to list
  if (fHelperPID)
    fListOfHistos->Add(fHelperPID);
  // add TMap to list
  if (fMap)
    fListOfHistos->Add(fMap);

  fListOfHistos->Add(new TH2F("processIDs", ";#Delta#phi;process id", 100, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), kPNoProcess + 1, -0.5, kPNoProcess + 0.5));
  fListOfHistos->Add(new TH1F("eventStat", ";;events", 4, -0.5, 3.5));
  fListOfHistos->Add(new TH2F("mixedDist", ";centrality;tracks;events", 101, 0, 101, 200, 0, fMixingTracks * 1.5));
  fListOfHistos->Add(new TH2F("mixedDist2", ";centrality;events;events", 101, 0, 101, 100, -0.5, 99.5));
  fListOfHistos->Add(new TH2F("referenceMultiplicity", ";centrality;tracks;events", 101, 0, 101, 200, -0.5, 199.5));
  if (fCentralityMethod == "V0A_MANUAL") {
    fListOfHistos->Add(new TH2F("V0AMult", "V0A multiplicity;V0A multiplicity;V0A multiplicity (scaled)", 1000, -.5, 999.5, 1000, -.5, 999.5));
    fListOfHistos->Add(new TH2F("V0AMultCorrelation", "V0A multiplicity;V0A multiplicity;SPD tracklets", 1000, -.5, 999.5, 1000, -.5, 999.5));
  }
  if (fTriggersFromDetector == 1 || fTriggersFromDetector == 2 || fAssociatedFromDetector == 1 || fAssociatedFromDetector == 2)
    fListOfHistos->Add(new TH1F("V0SingleCells", "V0 single cell multiplicity;multiplicity;events", 100, -0.5, 99.5));
  if (fTriggersFromDetector == 3 || fAssociatedFromDetector == 3)
    fListOfHistos->Add(new TH1F("DphiTrklets", "tracklets Dphi;#Delta#phi,trklets (mrad);entries", 100, -100, 100));
  if (fCheckCertainSpecies > 0)
    fListOfHistos->Add(new TH2F("checkSpecies", ";eta;pt;particles", 20, -1, 1, 40, 0, 10));

  if (fCentralityMethod == "ZNAC")
    fListOfHistos->Add(new TH1D("ZNA+C_energy", "ZNA+C_energy", 4100, -100, 4000));
  if (fCentralityMethod == "MCGen_V0M")
    fListOfHistos->Add(new TH2D("Mult_MCGen_V0M", "Mult_MCGen_V0M", 1010, -9.5, 1000.5, 1010, -9.5, 1000.5));
  if (fCentralityMethod == "MCGen_CL1")
    fListOfHistos->Add(new TH2D("Mult_MCGen_CL1", "Mult_MCGen_CL1", 1010, -9.5, 1000.5, 1010, -9.5, 1000.5));
  if (fV0CL1PileUp) {
    fListOfHistos->Add(new TH2I("fHistGlobalvsV0BeforePileUpCuts", "fHistGlobalvsV0BeforePileUpCuts;V0;CL1", 100, 0, 100, 100, 0, 100));
    fListOfHistos->Add(new TH2I("fHistGlobalvsV0AfterPileUpCuts", "fHistGlobalvsV0AfterPileUpCuts;V0;CL1", 100, 0, 100, 100, 0, 100));
  }
  if (fESDTPCTrackPileUp) {
    fListOfHistos->Add(new TH2I("fHistGlobalvsESDBeforePileUpCuts", "fHistGlobalvsESDBeforePileUpCuts;nTracks;multESD", 100, 0, 30000, 100, 0, 30000));
    fListOfHistos->Add(new TH2I("fHistGlobalvsESDAfterPileUpCuts", "fHistGlobalvsESDAfterPileUpCuts;nTracks;multESD", 100, 0, 30000, 100, 0, 30000));
  }
  if (fTPCITSTOFPileUp) {
    fListOfHistos->Add(new TH2I("fHistV0MvsTPCoutBeforePileUpCuts", "fHistV0MvsTPCoutBeforePileUpCuts;ntrkTPCout;multVZERO", 100, 0, 40000, 100, 0, 40000));
    fListOfHistos->Add(new TH2I("fHistV0MvsTPCoutAfterPileUpCuts", "fHistV0MvsTPCoutAfterPileUpCuts;ntrkTPCout;multVZERO", 100, 0, 40000, 100, 0, 40000));
  }
  Int_t nCentralityBins = fHistos->GetUEHist(2)->GetEventHist()->GetNBins(1);
  Double_t* centralityBins = (Double_t*) fHistos->GetUEHist(2)->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray();

  if (fFillYieldRapidity) {
    const Int_t nPtBins = 200;
    Double_t ptBins[nPtBins+1];
    for (int i=0; i<=nPtBins; i++)
      ptBins[i] = 20.0 / nPtBins * i;

    const Int_t nyBins = 200;
    Double_t yBins[nyBins+1];
    for (int i=0; i<=nyBins; i++)
      yBins[i] = -10.0 + 20.0 / nyBins * i;

    fListOfHistos->Add(new TH3F("yieldsRapidity", ";centrality;pT;y", nCentralityBins, centralityBins, nPtBins, ptBins, nyBins, yBins));
  }

  PostData(0,fListOfHistos);

  // Add task configuration to output list
  AddSettingsTree();

  // event mixing
  Int_t poolsize = -1; // Maximum number of events, -1 means no limit

  const Int_t kNZvtxBins = 10+(1+10)*4;
  // bins for further buffers are shifted by 100 cm
  Double_t vertexBins[kNZvtxBins+1] = {-10, -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,   10, 
                                        90,  92,  94,  96,  98,  100, 102, 104, 106, 108, 110, 
                                        190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 
                                        290, 292, 294, 296, 298, 300, 302, 304, 306, 308, 310, 
                                        390, 392, 394, 396, 398, 400, 402, 404, 406, 408, 410 };

  Int_t nZvtxBins = kNZvtxBins;
  Double_t* zvtxbin = vertexBins;

  if (fMode == 0 && fHistos->GetUEHist(2)->GetEventHist()->GetNVar() > 2) {
    nZvtxBins = fHistos->GetUEHist(2)->GetEventHist()->GetNBins(2);
    zvtxbin = (Double_t*) fHistos->GetUEHist(2)->GetEventHist()->GetAxis(2, 0)->GetXbins()->GetArray();
  }

  Int_t nPsiBins = 1;
  Double_t psibins[2] = {-999.,999.};

  Int_t nPtBins = 1;
  Double_t defaultPtBins[2] = {-9999., 9999.};
  Double_t* ptbins = defaultPtBins;

  // Retrieve event pool pt binning only if pool should be devided in bins of pt
  if (fUsePtBinnedEventPool)
    if (fHistos->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetNBins(1)) {
      nPtBins  = fHistos->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetNBins(1);
      ptbins = (Double_t*) fHistos->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetAxis(1, 0)->GetXbins()->GetArray();
    }

  // Create default event pool in case no external pool is given
  if (!fPoolMgr) {
    fPoolMgr = new AliEventPoolManager(poolsize, fMixingTracks, nCentralityBins, centralityBins, nZvtxBins, zvtxbin, nPsiBins, psibins, nPtBins, ptbins);
    fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
  }

  // Check binning of pool manager (basic dimensional check for the time being)
  if ( (fPoolMgr->GetNumberOfMultBins() != nCentralityBins) || (fPoolMgr->GetNumberOfZVtxBins() != nZvtxBins) || (fPoolMgr->GetNumberOfPtBins() != nPtBins) )
    AliFatal("Binning of given pool manager not compatible with binning of correlation task!");

  // If some bins of the pool should be saved, fEventPoolOutputList must be given
  // using AddEventPoolToOutput()
  // Note that this is in principle also possible, if an external poolmanager was given
  for (UInt_t i = 0; i < fEventPoolOutputList.size(); i++) {
    Double_t minCent = fEventPoolOutputList[i][0];
    Double_t maxCent = fEventPoolOutputList[i][1];
    Double_t minZvtx = fEventPoolOutputList[i][2];
    Double_t maxZvtx = fEventPoolOutputList[i][3];
    Double_t minPt   = fEventPoolOutputList[i][4];
    Double_t maxPt   = fEventPoolOutputList[i][5];

    fPoolMgr->SetSaveFlag(minCent, maxCent, minZvtx, maxZvtx, 0, 0, minPt, maxPt);
  }

  // Basic checks and printing of pool properties
  fPoolMgr->Validate();

  // save to output if requested
  if (fEventPoolOutputList.size())
    fListOfHistos->Add(fPoolMgr);
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::Exec(Option_t */*option*/)
{
  // exec (per event)
  fAnalyseUE->NextEvent();

  // receive ESD pointer if we are not running AOD analysis
  if (!fAOD) {
    AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    if (handler && handler->InheritsFrom("AliESDInputHandler"))
      fESD = (AliESDEvent*)((AliESDInputHandler*)handler)->GetEvent();
  }

  if (fMode) {
    // correction mode

    if (fMcHandler)
      fMcEvent = fMcHandler->MCEvent();

    if (fAOD) {
      // array of MC particles
      fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
        AliFatal("No array of MC particles found !!!");
    }

    AnalyseCorrectionMode();
  }
  else AnalyseDataMode();
}

/******************** ANALYSIS METHODS *****************************/

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  settingsTree->Branch("fnTracksVertex", &fnTracksVertex,"nTracksVertex/I");
  settingsTree->Branch("fZVertex", &fZVertex,"ZVertex/D");
  settingsTree->Branch("fAcceptOnlyMuEvents", &fAcceptOnlyMuEvents,"AcceptOnlyMuEvents/O");
  //settingsTree->Branch("fCentralityMethod", fCentralityMethod.Data(),"CentralityMethod/C");
  settingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  settingsTree->Branch("fTrackEtaCutMin", &fTrackEtaCutMin, "TrackEtaCutMin/D");
  settingsTree->Branch("fTrackPhiCutEvPlMin", &fTrackPhiCutEvPlMin, "TrackPhiCutEvPlMin/D");
  settingsTree->Branch("fTrackPhiCutEvPlMax", &fTrackPhiCutEvPlMax, "TrackPhiCutEvPlMax/D");
  settingsTree->Branch("fOnlyOneEtaSide", &fOnlyOneEtaSide,"OnlyOneEtaSide/I");
  settingsTree->Branch("fOnlyOneAssocEtaSide", &fOnlyOneAssocEtaSide,"OnlyOneAssocEtaSide/I");
  settingsTree->Branch("fPtMin", &fPtMin, "PtMin/D");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fSharedClusterCut", &fSharedClusterCut,"SharedClusterCut/D");
  settingsTree->Branch("fCrossedRowsCut", &fCrossedRowsCut,"CrossedRowsCut/I");
  settingsTree->Branch("fFoundFractionCut", &fFoundFractionCut,"FoundFractionCut/D");
  settingsTree->Branch("fTrackStatus", &fTrackStatus,"TrackStatus/I");
  settingsTree->Branch("fSelectBit", &fSelectBit,"EventSelectionBit/I");
  settingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  settingsTree->Branch("fParticleSpeciesTrigger", &fParticleSpeciesTrigger,"ParticleSpeciesTrigger/I");
  settingsTree->Branch("fParticleSpeciesAssociated", &fParticleSpeciesAssociated,"ParticleSpeciesAssociated/I");
  settingsTree->Branch("fCheckMotherPDG", &fCheckMotherPDG,"CheckMotherPDG/I");
  settingsTree->Branch("fSelectCharge", &fSelectCharge,"SelectCharge/I");
  settingsTree->Branch("fTriggerSelectCharge", &fTriggerSelectCharge,"TriggerSelectCharge/I");
  settingsTree->Branch("fAssociatedSelectCharge", &fAssociatedSelectCharge,"fAssociatedSelectCharge/I");
  settingsTree->Branch("fTriggerRestrictEta", &fTriggerRestrictEta,"TriggerRestrictEta/D");
  settingsTree->Branch("fEtaOrdering", &fEtaOrdering,"EtaOrdering/O");
  settingsTree->Branch("fCutConversionsV", &fCutConversionsV,"CutConversionsV/D");
  settingsTree->Branch("fCutResonancesV", &fCutResonancesV,"CutResonancesV/D");
  settingsTree->Branch("fCutOnCustomMass", &fCutOnCustomMass,"CutOnCustomMass/D");
  settingsTree->Branch("fCutOnCustomFirst", &fCutOnCustomFirst,"CutOnCustomFirst/D");
  settingsTree->Branch("fCutOnCustomSecond", &fCutOnCustomSecond,"CutOnCustomSecond/D");
  settingsTree->Branch("fCutOnCustomV", &fCutOnCustomV,"CutOnCustomV/D");
  settingsTree->Branch("fCutOnPhi", &fCutOnPhi,"CutOnPhi/O");
  settingsTree->Branch("fCutOnPhiV", &fCutOnPhiV,"CutOnPhiV/D");
  settingsTree->Branch("fCutOnRho", &fCutOnRho,"CutOnRho/O");
  settingsTree->Branch("fCutOnRhoV", &fCutOnRhoV,"CutOnRhoV/D");
  settingsTree->Branch("fCutOnLambdaV", &fCutOnLambdaV,"CutOnLambdaV/D");
  settingsTree->Branch("fCutOnK0sV", &fCutOnK0sV,"CutOnK0sV/D");
  settingsTree->Branch("fRejectResonanceDaughters", &fRejectResonanceDaughters,"RejectResonanceDaughters/I");
  settingsTree->Branch("fFillpT", &fFillpT,"FillpT/O");
  settingsTree->Branch("fMixingTracks", &fMixingTracks,"MixingTracks/I");
  settingsTree->Branch("fSkipTrigger", &fSkipTrigger,"SkipTrigger/O");
  settingsTree->Branch("fInjectedSignals", &fInjectedSignals,"InjectedSignals/O");
  settingsTree->Branch("fRandomizeReactionPlane", &fRandomizeReactionPlane,"RandomizeReactionPlane/O");
  settingsTree->Branch("fRejectCentralityOutliers", &fRejectCentralityOutliers,"RejectCentralityOutliers/O");
  settingsTree->Branch("fRejectZeroTrackEvents", &fRejectZeroTrackEvents,"RejectZeroTrackEvents/O");
  settingsTree->Branch("fRemoveWeakDecays", &fRemoveWeakDecays,"RemoveWeakDecays/O");
  settingsTree->Branch("fRemoveDuplicates", &fRemoveDuplicates,"RemoveDuplicates/O");
  settingsTree->Branch("fSkipFastCluster", &fSkipFastCluster,"SkipFastCluster/O");
  settingsTree->Branch("fWeightPerEvent", &fWeightPerEvent,"WeightPerEvent/O");
  settingsTree->Branch("fPtOrder", &fPtOrder,"PtOrder/O");
  settingsTree->Branch("fTriggersFromDetector", &fTriggersFromDetector,"TriggersFromDetector/I");
  settingsTree->Branch("fAssociatedFromDetector", &fAssociatedFromDetector,"AssociatedFromDetector/I");
  settingsTree->Branch("fUseUncheckedCentrality", &fUseUncheckedCentrality,"UseUncheckedCentrality/O");
  settingsTree->Branch("fCheckCertainSpecies", &fCheckCertainSpecies,"fCheckCertainSpecies/I");
  settingsTree->Branch("fRemoveWeakDecaysInMC", &fRemoveWeakDecaysInMC,"RemoveWeakDecaysInMC/O");
  settingsTree->Branch("fFillYieldRapidity", &fFillYieldRapidity,"fFillYieldRapidity/O");
  settingsTree->Branch("fFillCorrelationsRapidity", &fFillCorrelationsRapidity,"fFillCorrelationsRapidity/O");
  settingsTree->Branch("fUseNewCentralityFramework", &fUseNewCentralityFramework,"fUseNewCentralityFramework/O");
  settingsTree->Branch("fTwoTrackEfficiencyCut", &fTwoTrackEfficiencyCut,"TwoTrackEfficiencyCut/D");
  settingsTree->Branch("fTwoTrackCutMinRadius", &fTwoTrackCutMinRadius,"TwoTrackCutMinRadius/D");

  //fCustomBinning

  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AnalyseCorrectionMode()
{
  // Run the analysis on MC to get the correction maps
  //

  if ( fDebug > 3 )
    AliInfo( " Processing event in Corrections mode ..." );

  TObject* mc = fArrayMC;
  if (!mc)
    mc = fMcEvent;

  // Support for ESD and AOD based analysis
  AliVEvent* inputEvent = fAOD;
  if (!inputEvent)
    inputEvent = fESD;

  Double_t centrality = GetCentrality(inputEvent, mc);

  Float_t bSign = 0;

  if (inputEvent) {
    fHistos->SetRunNumber(inputEvent->GetRunNumber());
    bSign = (inputEvent->GetMagneticField() > 0) ? 1 : -1;
  }

  // count all events
  fHistos->FillEvent(centrality, -1);

  if (centrality < 0)
    return;

  // Only consider MC events within the vtx-z region used also as cut on the reconstructed vertex
  TObject* vertexSupplier = fMcEvent;
  if (fAOD) // AOD
    vertexSupplier = fAOD->FindListObject(AliAODMCHeader::StdBranchName());

  if (!fAnalyseUE->VertexSelection(vertexSupplier, 0, fZVertex))
    return;

  Float_t zVtx = 0;
  if (fAOD)
    zVtx = ((AliAODMCHeader*) vertexSupplier)->GetVtxZ();
  else
    zVtx = fMcEvent->GetPrimaryVertex()->GetZ();
  Float_t weight = 1;
  if (fFillpT)
    weight = -1;

  Double_t evtPlanePhi = -999.; //A value outside [-pi/2,pi/2] will be ignored
  if (fTrackPhiCutEvPlMax > 0.0001)
    if (!InitiateEventPlane(evtPlanePhi, inputEvent)) return; //Reject event if plane is not available

  // For productions with injected signals, figure out above which label to skip particles/tracks
  Int_t skipParticlesAbove = 0;
  if (fInjectedSignals) {
    AliGenEventHeader* eventHeader = 0;
    Int_t headers = 0;

    if (fMcEvent) {
      // ESD
      AliHeader* header = (AliHeader*) fMcEvent->Header();
      if (!header)
        AliFatal("fInjectedSignals set but no MC header found");

      AliGenCocktailEventHeader* cocktailHeader = dynamic_cast<AliGenCocktailEventHeader*> (header->GenEventHeader());
      if (!cocktailHeader) {
        header->Dump();
        AliFatal("fInjectedSignals set but no MC cocktail header found");
      }

      headers = cocktailHeader->GetHeaders()->GetEntries();
      eventHeader = dynamic_cast<AliGenEventHeader*> (cocktailHeader->GetHeaders()->First());

      if (fDebug > 4) {
        for (Int_t i=0; i<cocktailHeader->GetHeaders()->GetEntries(); i++) {
          AliGenEventHeader* headerTmp = dynamic_cast<AliGenEventHeader*> (cocktailHeader->GetHeaders()->At(i));
          if (headerTmp)
            Printf("%d particles in header:", headerTmp->NProduced());
          cocktailHeader->GetHeaders()->At(i)->Dump();
        }
      }
    }
    else {
      // AOD
      AliAODMCHeader* header = (AliAODMCHeader*) fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      if (!header)
        AliFatal("fInjectedSignals set but no MC header found");

      headers = header->GetNCocktailHeaders();
      eventHeader = header->GetCocktailHeader(0);
    }

    if (!eventHeader) {
      // We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
      // (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
      AliError("First event header not found. Skipping this event.");
      fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
      return;
    }

    skipParticlesAbove = eventHeader->NProduced();
    AliInfo(Form("Injected signals in this event (%d headers). Keeping particles/tracks of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove));
  }

  if (fCentralityWeights && !AcceptEventCentralityWeight(centrality)) {
    AliInfo(Form("Rejecting event because of centrality weighting: %f", centrality));
    fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
    return;
  }

  // debug for certain species
  if (fCheckCertainSpecies > 0) {
    // need to get also neutral particles here
    TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, fParticleSpeciesTrigger, kTRUE, kTRUE, evtPlanePhi, kFALSE);
    for (Int_t i=0; i<tmpList->GetEntriesFast(); i++) {
      AliMCParticle* particle = dynamic_cast<AliMCParticle*> (tmpList->UncheckedAt(i));
      if (!particle)
        continue;
      if (TMath::Abs(particle->PdgCode()) != fCheckCertainSpecies)
        continue;
      ((TH2F*) fListOfHistos->FindObject("checkSpecies"))->Fill(particle->Eta(), particle->Pt());
    }
    delete tmpList;
  }

  // Get MC primaries
  // triggers
  TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, fParticleSpeciesTrigger, kTRUE, kTRUE, evtPlanePhi);
  CleanUp(tmpList, mc, skipParticlesAbove);
  if (fFillYieldRapidity) {
    for (Int_t i=0; i<tmpList->GetEntriesFast(); i++) {
      AliVParticle* particle = dynamic_cast<AliVParticle*> (tmpList->UncheckedAt(i));
      if (!particle)
        continue;
      ((TH3F*) fListOfHistos->FindObject("yieldsRapidity"))->Fill(centrality, particle->Pt(), particle->Y());
    }
  }
  TObjArray* tracksMC = CloneAndReduceTrackList(tmpList);
  delete tmpList;

  // associated
  TObjArray* tracksCorrelateMC = tracksMC;
  if (fParticleSpeciesAssociated != fParticleSpeciesTrigger || fTrackPhiCutEvPlMax > 0.0001) {
    tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, fParticleSpeciesAssociated, kTRUE);
    CleanUp(tmpList, mc, skipParticlesAbove);
    tracksCorrelateMC = CloneAndReduceTrackList(tmpList);
    delete tmpList;
  }

  if (fRandomizeReactionPlane) {
    Double_t centralityDigits = centrality*1000. - (Int_t)(centrality*1000.);
    Double_t angle = TMath::TwoPi() * centralityDigits;
    AliInfo(Form("Shifting phi of all tracks by %f (digits %f)", angle, centralityDigits));
    ShiftTracks(tracksMC, angle);
    if (tracksCorrelateMC != tracksMC)
      ShiftTracks(tracksCorrelateMC, angle);
  }

  if (fFillOnlyStep0)
    zVtx = 0;

  // Event selection based on number of number of MC particles
  if (fRejectZeroTrackEvents && tracksMC->GetEntriesFast() == 0) {
    AliInfo(Form("Rejecting event due to kinematic selection: %f %d", centrality, tracksMC->GetEntriesFast()));
    fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
    if (tracksMC != tracksCorrelateMC)
      delete tracksCorrelateMC;
    delete tracksMC;
    return;
  }

  // (MC-true all particles)
  // STEP 0
  fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepAll, tracksMC, tracksCorrelateMC, weight);

  // mixed event
  if (fFillMixed) {
    for (Int_t iPool=0; iPool<fPoolMgr->GetNumberOfPtBins(); iPool++) {
      AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx, 0., iPool);
      if (fFillOnlyStep0) {
        ((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
        ((TH2F*) fListOfHistos->FindObject("mixedDist2"))->Fill(centrality, pool->GetCurrentNEvents());
      }
      if (pool->IsReady())
        for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++)
          fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepAll, tracksMC, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix == 0));
      pool->UpdatePool(CloneAndReduceTrackList(tracksCorrelateMC, pool->GetPtMin(), pool->GetPtMax()));
    }
  }

//   Printf("trigger: %d", ((AliInputEventHandler*)fInputHandler)->IsEventSelected());

  // Trigger selection ************************************************
  if (!fFillOnlyStep0 && (fSkipTrigger || fAnalyseUE->TriggerSelection(fInputHandler))) {
    // (MC-true all particles)
    // STEP 1
    if (!fReduceMemoryFootprint)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTriggered, tracksMC, tracksCorrelateMC, weight);
    else
      fHistos->FillEvent(centrality, AliUEHist::kCFStepTriggered);

    if (!inputEvent) {
      AliFatal("UNEXPECTED: inputEvent is 0. Trigger selection should have failed");
      return;
    }

    // Vertex selection *************************************************
    if (fAnalyseUE->VertexSelection(inputEvent, fnTracksVertex, fZVertex)) {
      // fill here for tracking efficiency
      // loop over particle species

      for (Int_t particleSpecies = 0; particleSpecies < 4; particleSpecies++) {
        TObjArray* primMCParticles = fAnalyseUE->GetAcceptedParticles(mc, 0x0, kTRUE, particleSpecies, kTRUE, kTRUE, evtPlanePhi);
        TObjArray* primRecoTracksMatched = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, particleSpecies, kTRUE, kFALSE, evtPlanePhi);
        TObjArray* allRecoTracksMatched  = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, particleSpecies, kTRUE, kFALSE, evtPlanePhi);
        TObjArray* primRecoTracksMatchedPID = 0;
        TObjArray* allRecoTracksMatchedPID  = 0;

        if (fHelperPID) {
          primRecoTracksMatchedPID = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, particleSpecies, kTRUE, kTRUE, evtPlanePhi);
          allRecoTracksMatchedPID  = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, particleSpecies, kTRUE, kTRUE, evtPlanePhi);
        }

        CleanUp(primMCParticles, mc, skipParticlesAbove);
        CleanUp(primRecoTracksMatched, mc, skipParticlesAbove);
        CleanUp(allRecoTracksMatched, mc, skipParticlesAbove);
        CleanUp(primRecoTracksMatchedPID, mc, skipParticlesAbove);
        CleanUp(allRecoTracksMatchedPID, mc, skipParticlesAbove);

        // select charges
        if (fTriggerSelectCharge != 0) {
          SelectCharge(primMCParticles);
          SelectCharge(primRecoTracksMatched);
          SelectCharge(allRecoTracksMatched);
          SelectCharge(primRecoTracksMatchedPID);
          SelectCharge(allRecoTracksMatchedPID);
        }

        fHistos->FillTrackingEfficiency(primMCParticles, primRecoTracksMatched, allRecoTracksMatched, primRecoTracksMatchedPID, allRecoTracksMatchedPID, 0, particleSpecies, centrality, zVtx);

//      Printf("%d --> %d %d %d", particleSpecies, primMCParticles->GetEntries(), primRecoTracksMatched->GetEntries(), allRecoTracksMatched->GetEntries());

        delete primMCParticles;
        delete primRecoTracksMatched;
        delete allRecoTracksMatched;
      }

      TObjArray* fakeParticles = fAnalyseUE->GetFakeParticles(inputEvent, mc, kFALSE, -1, kTRUE);
      CleanUp((TObjArray*) fakeParticles->At(0), mc, skipParticlesAbove);
      CleanUp((TObjArray*) fakeParticles->At(1), mc, skipParticlesAbove);

      fHistos->FillTrackingEfficiency(0, 0, 0, 0, 0, (TObjArray*) fakeParticles->At(2), 0, centrality, zVtx);
      fHistos->FillFakePt(fakeParticles, centrality);
//       Printf(">>>>> %d %d %d fakes", ((TObjArray*) fakeParticles->At(0))->GetEntriesFast(), ((TObjArray*) fakeParticles->At(1))->GetEntriesFast(), ((TObjArray*) fakeParticles->At(2))->GetEntriesFast());
      delete fakeParticles;

      // (MC-true all particles)
      // STEP 2
      if (!fReduceMemoryFootprint)
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepVertex, tracksMC, tracksCorrelateMC, weight);
      else
        fHistos->FillEvent(centrality, AliUEHist::kCFStepVertex);

      // Get MC primaries that match reconstructed track
      // triggers
      tmpList = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, fParticleSpeciesTrigger, kTRUE, kTRUE, evtPlanePhi);
      CleanUp(tmpList, mc, skipParticlesAbove);
      TObjArray* tracksRecoMatchedPrim = CloneAndReduceTrackList(tmpList);
      delete tmpList;

      // associated
      TObjArray* tracksCorrelateRecoMatchedPrim = tracksRecoMatchedPrim;
      if (fParticleSpeciesAssociated != fParticleSpeciesTrigger || fTrackPhiCutEvPlMax > 0.0001) {
        tmpList = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, fParticleSpeciesAssociated, kTRUE);
        CleanUp(tmpList, mc, skipParticlesAbove);
        tracksCorrelateRecoMatchedPrim = CloneAndReduceTrackList(tmpList);
        delete tmpList;
      }

      // (RECO-matched (quantities from MC particle) primary particles)
      // STEP 4
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTrackedOnlyPrim, tracksRecoMatchedPrim, tracksCorrelateRecoMatchedPrim, weight);

      // mixed event
      if (fFillMixed) {
        for (Int_t iPool=0; iPool<fPoolMgr->GetNumberOfPtBins(); iPool++) {
          AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx + 200, 0., iPool);
          if (pool->IsReady())
            for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++)
              fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTrackedOnlyPrim, tracksRecoMatchedPrim, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix == 0));
          pool->UpdatePool(CloneAndReduceTrackList(tracksCorrelateRecoMatchedPrim, pool->GetPtMin(), pool->GetPtMax()));
        }
      }

      // Get MC primaries + secondaries that match reconstructed track
      // triggers
      tmpList = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, fParticleSpeciesTrigger, kTRUE, kTRUE, evtPlanePhi);
      CleanUp(tmpList, mc, skipParticlesAbove);
      TObjArray* tracksRecoMatchedAll = CloneAndReduceTrackList(tmpList);
      delete tmpList;

      // associated
      TObjArray* tracksCorrelateRecoMatchedAll = tracksRecoMatchedAll;
      if (fParticleSpeciesAssociated != fParticleSpeciesTrigger || fTrackPhiCutEvPlMax > 0.0001) {
        tmpList = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, fParticleSpeciesAssociated, kTRUE);
        CleanUp(tmpList, mc, skipParticlesAbove);
        tracksCorrelateRecoMatchedAll = CloneAndReduceTrackList(tmpList);
        delete tmpList;
      }

      // (RECO-matched (quantities from MC particle) all particles)
      // STEP 5
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTracked, tracksRecoMatchedAll, tracksCorrelateRecoMatchedAll, weight);

      // mixed event
      if (fFillMixed) {
        for (Int_t iPool=0; iPool<fPoolMgr->GetNumberOfPtBins(); iPool++) {
          AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx + 300, 0., iPool);
          if (pool->IsReady())
            for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++)
              fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTracked, tracksRecoMatchedAll, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix == 0));
          pool->UpdatePool(CloneAndReduceTrackList(tracksCorrelateRecoMatchedAll, pool->GetPtMin(), pool->GetPtMax()));
        }
      }

      // Get RECO tracks
      // triggers
      tmpList = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, fParticleSpeciesTrigger, kTRUE, kTRUE, evtPlanePhi);
      CleanUp(tmpList, mc, skipParticlesAbove);
      TObjArray* tracks = CloneAndReduceTrackList(tmpList);
      delete tmpList;

      // associated
      TObjArray* tracksCorrelate = tracks;
      if (fParticleSpeciesAssociated != fParticleSpeciesTrigger || fTrackPhiCutEvPlMax > 0.0001) {
        tmpList = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, fParticleSpeciesAssociated, kTRUE);
        CleanUp(tmpList, mc, skipParticlesAbove);
        tracksCorrelate = CloneAndReduceTrackList(tmpList);
        delete tmpList;
      }

      // (RECO all tracks)
      // STEP 6
      if (!fSkipStep6)
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, tracksCorrelate, weight);

      // two track cut, STEP 8
      if (fTwoTrackEfficiencyCut > 0)
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracks, tracksCorrelate, weight, kTRUE, kTRUE, bSign, fTwoTrackEfficiencyCut);

      // apply correction efficiency, STEP 9 and 10
      if (fEfficiencyCorrectionTriggers || fEfficiencyCorrectionAssociated) {

        // all two track cuts disabled, STEP 9
        if (!fSkipStep9)
          fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy2, tracks, tracksCorrelate, weight, kTRUE, kFALSE, bSign, .0, kTRUE);

        // STEP 10
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepCorrected, tracks, tracksCorrelate, weight, kTRUE, kTRUE, bSign, fTwoTrackEfficiencyCut, kTRUE);
      }

      // mixed event
      if (fFillMixed) {
        for (Int_t iPool=0; iPool<fPoolMgr->GetNumberOfPtBins(); iPool++) {
          AliEventPool* pool2 = fPoolMgr->GetEventPool(centrality, zVtx + 100, 0., iPool);
          ((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool2->NTracksInPool());
          ((TH2F*) fListOfHistos->FindObject("mixedDist2"))->Fill(centrality, pool2->GetCurrentNEvents());
          if (pool2->IsReady()) {
            for (Int_t jMix=0; jMix<pool2->GetCurrentNEvents(); jMix++) {
              // STEP 6
              if (!fSkipStep6)
                fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, pool2->GetEvent(jMix), 1.0 / pool2->GetCurrentNEvents(), (jMix == 0));

              // two track cut, STEP 8
              if (fTwoTrackEfficiencyCut > 0)
                fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracks, pool2->GetEvent(jMix), 1.0 / pool2->GetCurrentNEvents(), (jMix == 0), kTRUE, bSign, fTwoTrackEfficiencyCut);

              // apply correction efficiency, STEP 9 and 10
              if (fEfficiencyCorrectionTriggers || fEfficiencyCorrectionAssociated) {

                // all two track cuts disabled, STEP 9
                if (!fSkipStep9)
                  fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy2, tracks, pool2->GetEvent(jMix), 1.0 / pool2->GetCurrentNEvents(), (jMix == 0), kFALSE, bSign, .0, kTRUE);

                // STEP 10
                fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepCorrected, tracks, pool2->GetEvent(jMix), 1.0 / pool2->GetCurrentNEvents(), (jMix == 0), kTRUE, bSign, fTwoTrackEfficiencyCut, kTRUE);
              }
            }
          }
          pool2->UpdatePool(CloneAndReduceTrackList(tracksCorrelate, pool2->GetPtMin(), pool2->GetPtMax()));
        }
      }

      if (0 && !fReduceMemoryFootprint) {
        // make list of secondaries (matched with MC)
        TObjArray* tracksRecoMatchedSecondaries = new TObjArray;
        for (Int_t i=0; i<tracksRecoMatchedAll->GetEntriesFast(); i++)
          if (((AliAODMCParticle*)tracksRecoMatchedAll->At(i))->IsPhysicalPrimary() == kFALSE)
            tracksRecoMatchedSecondaries->Add(tracksRecoMatchedAll->At(i));

        // Study: Use only secondaries as trigger particles and plot the correlation vs. all particles; store in step 9
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy2, tracksRecoMatchedSecondaries, tracksRecoMatchedAll, weight);

        // Study: Use only primaries as trigger particles and plot the correlation vs. secondaries; store in step 8
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracksRecoMatchedPrim, tracksRecoMatchedSecondaries, weight);

        // plot delta phi vs process id of secondaries
        // trigger particles: primaries in 4 < pT < 10
        // associated particles: secondaries in 1 < pT < 10

        for (Int_t i=0; i<tracksRecoMatchedPrim->GetEntriesFast(); i++) {
          AliVParticle* triggerParticle = (AliVParticle*) tracksRecoMatchedPrim->At(i);

          if (triggerParticle->Pt() < 4 || triggerParticle->Pt() > 10)
            continue;

          for (Int_t j=0; j<tracksRecoMatchedSecondaries->GetEntriesFast(); j++) {
            AliAODMCParticle* particle = (AliAODMCParticle*) tracksRecoMatchedSecondaries->At(j);

            if (particle->Pt() < 1 || particle->Pt() > 10)
              continue;

            if (particle->Pt() > triggerParticle->Pt())
              continue;

            Double_t deltaPhi = triggerParticle->Phi() - particle->Phi();
            if (deltaPhi > 1.5 * TMath::Pi())
              deltaPhi -= TMath::TwoPi();
            if (deltaPhi < -0.5 * TMath::Pi())
              deltaPhi += TMath::TwoPi();

            Int_t processID = fMcEvent->Stack()->Particle(particle->GetLabel())->GetUniqueID();

            ((TH2F*) fListOfHistos->FindObject("processIDs"))->Fill(deltaPhi, processID);
          }
        }

        delete tracksRecoMatchedSecondaries;
      }

      if (tracksCorrelateRecoMatchedPrim != tracksRecoMatchedPrim)
        delete tracksCorrelateRecoMatchedPrim;
      delete tracksRecoMatchedPrim;

      if (tracksCorrelateRecoMatchedAll != tracksRecoMatchedAll)
        delete tracksCorrelateRecoMatchedAll;
      delete tracksRecoMatchedAll;

      if (tracksCorrelate != tracks)
        delete tracksCorrelate;
      delete tracks;
    }
  }

  if (tracksMC != tracksCorrelateMC)
    delete tracksCorrelateMC;
  delete tracksMC;
}

//____________________________________________________________________
AliGenEventHeader* AliAnalysisTaskPhiCorrelations::GetFirstHeader()
{
  // get first MC header from either ESD/AOD (including cocktail header if available)

  if (fMcEvent) {
    // ESD
    AliHeader* header = (AliHeader*) fMcEvent->Header();
    if (!header)
      return 0;

    AliGenCocktailEventHeader* cocktailHeader = dynamic_cast<AliGenCocktailEventHeader*> (header->GenEventHeader());
    if (cocktailHeader)
      return dynamic_cast<AliGenEventHeader*> (cocktailHeader->GetHeaders()->First());

    return dynamic_cast<AliGenEventHeader*> (header->GenEventHeader());
  }
  else {
    // AOD
    AliAODMCHeader* header = (AliAODMCHeader*) fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!header)
      return 0;

    return header->GetCocktailHeader(0);
  }
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::AnalyseDataMode()
{

  // Run the analysis on DATA or MC to get raw distributions

  if ( fDebug > 3 )
    AliInfo( " Processing event in Data mode ..." );

  ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(0);

  if (!fInputHandler)
    return;

  // skip not selected events here (the AOD is not updated for those)
  if (!fSkipTrigger && !(fInputHandler->IsEventSelected() & fSelectBit))
    return;

  // skip fast cluster events here if requested
  if (fSkipFastCluster && (fInputHandler->IsEventSelected() & AliVEvent::kFastOnly))
    return;

  // Support for ESD and AOD based analysis
  AliVEvent* inputEvent = fAOD;
  if (!inputEvent)
    inputEvent = fESD;

  Double_t centrality = GetCentrality(inputEvent, 0);

  Float_t bSign = (inputEvent->GetMagneticField() > 0) ? 1 : -1;

  fHistos->SetRunNumber(inputEvent->GetRunNumber());

  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(centrality, AliUEHist::kCFStepAll);

  // Trigger selection ************************************************
  if (!fSkipTrigger && !fAnalyseUE->TriggerSelection(fInputHandler)) return;

  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(centrality, AliUEHist::kCFStepTriggered);

  // Pileup selection ************************************************
  if (fAnalysisUtils && fAnalysisUtils->IsPileUpEvent(inputEvent)) {
    // count the removed events
    fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);

    return;
  }

  if (fV0CL1PileUp) {
    AliMultSelection *multSelection = (AliMultSelection*) inputEvent->FindListObject("MultSelection");
    ((TH2I*)fListOfHistos->FindObject("fHistGlobalvsV0BeforePileUpCuts"))->Fill(multSelection->GetMultiplicityPercentile("V0M"),multSelection->GetMultiplicityPercentile("CL1"));
    if (TMath::Abs(multSelection->GetMultiplicityPercentile("V0M") - multSelection->GetMultiplicityPercentile("CL1")) > 7.5) {
      fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
      return;
    }
    ((TH2I*)fListOfHistos->FindObject("fHistGlobalvsV0AfterPileUpCuts"))->Fill(multSelection->GetMultiplicityPercentile("V0M"),multSelection->GetMultiplicityPercentile("CL1"));
  }

  if (fESDTPCTrackPileUp) {
    const Int_t nTracks = inputEvent->GetNumberOfTracks();
    Int_t multEsd = ((AliAODHeader*)inputEvent->GetHeader())->GetNumberOfESDTracks();
    ((TH2D*)fListOfHistos->FindObject("fHistGlobalvsESDBeforePileUpCuts"))->Fill(nTracks,multEsd);
    Int_t multTPC = 0;
    for (Int_t it = 0; it < nTracks; it++) {
      AliAODTrack* AODTrk = (AliAODTrack*)inputEvent->GetTrack(it);
     if (!AODTrk){ delete AODTrk; continue; }
     if (AODTrk->TestFilterBit(128)) {multTPC++;}
    } // end of for (Int_t it = 0; it < nTracks; it++)
    double fPileupLHC15oSlope = 3.38;
    double fPileupLHC15oOffset = 15000;
    if ((multEsd - fPileupLHC15oSlope*multTPC) > fPileupLHC15oOffset) {
      fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
      return;
    }
    ((TH2I*)fListOfHistos->FindObject("fHistGlobalvsESDAfterPileUpCuts"))->Fill(nTracks,multEsd);
  }

  if (fTPCITSTOFPileUp) {
    Int_t ntrkTPCout = 0;
    for (int it = 0; it < inputEvent->GetNumberOfTracks(); it++) {
      AliAODTrack* AODTrk = (AliAODTrack*)inputEvent->GetTrack(it);
      if ((AODTrk->GetStatus() & AliAODTrack::kTPCout) && AODTrk->GetID() > 0)
        ntrkTPCout++;
    }

    Double_t multVZERO =0;
    AliVVZERO *vzero = (AliVVZERO*)inputEvent->GetVZEROData();
    if (vzero) {
      for (int ich=0; ich < 64; ich++)
        multVZERO += vzero->GetMultiplicity(ich);
    }


    ((TH2I*)fListOfHistos->FindObject("fHistV0MvsTPCoutBeforePileUpCuts"))->Fill(ntrkTPCout, multVZERO);

    if (multVZERO < (-2200 + 2.5*ntrkTPCout + 1.2e-5*ntrkTPCout*ntrkTPCout)) {
      fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
      return;
    }
    ((TH2I*)fListOfHistos->FindObject("fHistV0MvsTPCoutAfterPileUpCuts"))->Fill(ntrkTPCout, multVZERO);
  }

  // Reject events without a muon in the muon arm ************************************************
  if (fAcceptOnlyMuEvents && !IsMuEvent())
    return;

  // Vertex selection *************************************************
  if (!fAnalyseUE->VertexSelection(inputEvent, fnTracksVertex, fZVertex))
    return;

  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(centrality, AliUEHist::kCFStepVertex);

  // fill V0 control histograms
  if (fCentralityMethod == "V0A_MANUAL") {
    ((TH2F*) fListOfHistos->FindObject("V0AMult"))->Fill(inputEvent->GetVZEROData()->GetMTotV0A(), centrality);
    if (fAOD)
      ((TH2F*) fListOfHistos->FindObject("V0AMultCorrelation"))->Fill(inputEvent->GetVZEROData()->GetMTotV0A(), fAOD->GetTracklets()->GetNumberOfTracklets());
  }

  // optimization
  if (centrality < 0)
    return;

  TObjArray* tracks = 0;

  Double_t evtPlanePhi = -999.; //A value outside [-pi/2,pi/2] will be ignored
  if (fTrackPhiCutEvPlMax > 0.0001)
    if (!InitiateEventPlane(evtPlanePhi, inputEvent)) return; //Reject event if plane is not available

  if (fTriggersFromDetector == 0)
    tracks = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, fParticleSpeciesTrigger, kTRUE, kTRUE, evtPlanePhi);
  else if (fTriggersFromDetector <= 7)
    tracks=GetParticlesFromDetector(inputEvent,fTriggersFromDetector);
  else
    AliFatal(Form("Invalid setting for fTriggersFromDetector: %d", fTriggersFromDetector));

  //Printf("Accepted %d tracks", tracks->GetEntries());

  // check for outlier in centrality vs number of tracks (rough constants extracted from correlation histgram)
  Bool_t reject = kFALSE;
  if (fRejectCentralityOutliers) {
    if (centrality > 40 && centrality <= 50 && tracks->GetEntriesFast() > 1160)
      reject = kTRUE;
    if (centrality > 50 && centrality <= 60 && tracks->GetEntriesFast() > 650)
      reject = kTRUE;
    if (centrality > 60 && centrality <= 70 && tracks->GetEntriesFast() > 370)
      reject = kTRUE;
    if (centrality > 70 && centrality <= 80 && tracks->GetEntriesFast() > 220)
      reject = kTRUE;
    if (centrality > 80 && centrality <= 90 && tracks->GetEntriesFast() > 130)
      reject = kTRUE;
    if (centrality > 90 && tracks->GetEntriesFast() > 75)
      reject = kTRUE;
  }

  if (reject) {
    AliInfo(Form("Rejecting event due to centrality vs tracks correlation: %f %d", centrality, tracks->GetEntriesFast()));
    fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
    delete tracks;
    return;
  }

  if (fRejectZeroTrackEvents && tracks->GetEntriesFast() == 0) {
    AliInfo(Form("Rejecting event because it has no tracks: %f %d", centrality, tracks->GetEntriesFast()));
    fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
    delete tracks;
    return;
  }

  if (fCentralityWeights && !AcceptEventCentralityWeight(centrality)) {
    AliInfo(Form("Rejecting event because of centrality weighting: %f", centrality));
    fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
    delete tracks;
    return;
  }

  // correlate particles with...
  TObjArray* tracksCorrelate = 0;
  if (fAssociatedFromDetector==0) {
    if (fParticleSpeciesAssociated != fParticleSpeciesTrigger || fTriggersFromDetector > 0 || fTrackPhiCutEvPlMax > 0.0001)
      tracksCorrelate = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, fParticleSpeciesAssociated, kTRUE);
  }
  else if (fAssociatedFromDetector <= 7) {
    if (fAssociatedFromDetector != fTriggersFromDetector)
      tracksCorrelate = GetParticlesFromDetector(inputEvent,fAssociatedFromDetector);
  }
  else
    AliFatal(Form("Invalid setting for fAssociatedFromDetector: %d", fAssociatedFromDetector));

  // reference multiplicity
  Int_t referenceMultiplicity = -1;
  if (fESD)
    referenceMultiplicity = AliESDtrackCuts::GetReferenceMultiplicity(fESD);
  else if (fAOD)
    referenceMultiplicity = tracks->GetEntriesFast(); // TODO to be replaced by the estimator once available in the AOD
//    referenceMultiplicity = ((AliVAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb05();

  ((TH2F*) fListOfHistos->FindObject("referenceMultiplicity"))->Fill(centrality, referenceMultiplicity);

  const AliVVertex* vertex = inputEvent->GetPrimaryVertex();
  Double_t zVtx = vertex->GetZ();

  Float_t weight = 1;
  if (fFillpT)
    weight = -1;

  // Fill containers at STEP 6 (reconstructed)
  if (centrality >= 0) {
    if (!fSkipStep6)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, tracksCorrelate, weight, kTRUE, kTRUE, 0, -1, kTRUE);

    ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(1);

    if (fTwoTrackEfficiencyCut > 0)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracks, tracksCorrelate, weight, kTRUE, kTRUE, bSign, fTwoTrackEfficiencyCut, kTRUE);

    if (!fSkipStep9)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy2, tracks, tracksCorrelate, weight, kTRUE, kFALSE, bSign, fTwoTrackEfficiencyCut, kTRUE);

  }

  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = CloneAndReduceTrackList(tracks);
  delete tracks;

  if (fFillMixed) {
    // event mixing

    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    for (Int_t iPool=0; iPool<fPoolMgr->GetNumberOfPtBins(); iPool++) {
      AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx, 0., iPool);

      if (!pool)
        AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, zVtx));

  //     pool->SetDebug(1);

      if (pool->IsReady()) {
        Int_t nMix = pool->GetCurrentNEvents();
  //       cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;

        ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
        ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(3, nMix);
        ((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
        ((TH2F*) fListOfHistos->FindObject("mixedDist2"))->Fill(centrality, nMix);

        // Fill mixed-event histos here
        for (Int_t jMix=0; jMix<nMix; jMix++) {
          TObjArray* bgTracks = pool->GetEvent(jMix);

          if (!fSkipStep6)
            fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracksClone, bgTracks, 1.0 / nMix, (jMix == 0), kTRUE, 0, -1, kTRUE);

          if (fTwoTrackEfficiencyCut > 0)
            fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracksClone, bgTracks, 1.0 / nMix, (jMix == 0), kTRUE, bSign, fTwoTrackEfficiencyCut, kTRUE);

          if (!fSkipStep9)
            fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy2, tracksClone, bgTracks, 1.0 / nMix, (jMix == 0), kFALSE, bSign, fTwoTrackEfficiencyCut, kTRUE);

        }
      }

      if (!pool->GetLockFlag()) {
        TObjArray* storeInMixed = 0;
        if (tracksCorrelate)
          storeInMixed = CloneAndReduceTrackList(tracksCorrelate, pool->GetPtMin(), pool->GetPtMax());
        else
          storeInMixed = CloneAndReduceTrackList(tracksClone, pool->GetPtMin(), pool->GetPtMax());

        // ownership is with the pool now
        pool->UpdatePool(storeInMixed);
      }
      //pool->PrintInfo();
    }
  }

  delete tracksClone;
  if (tracksCorrelate)
    delete tracksCorrelate;
}

Double_t AliAnalysisTaskPhiCorrelations::GetCentrality(AliVEvent* inputEvent, TObject* mc)
{
  // return centrality

  if (fCentralityMethod.Length() == 0)
    return 0;

  Double_t centrality = 0;

  if (fUseNewCentralityFramework)  {
    AliMultSelection *multSelection = (AliMultSelection*) inputEvent->FindListObject("MultSelection");
    if (!multSelection)
      AliFatal("MultSelection not found in input event");

    if (fUseUncheckedCentrality)
      centrality = multSelection->GetMultiplicityPercentile(fCentralityMethod, kFALSE);
    else
      centrality = multSelection->GetMultiplicityPercentile(fCentralityMethod, kTRUE);

    // error handling
    if (centrality > 100)
      centrality = -1;
  }
  else {
    AliCentrality *centralityObj = 0;

    if (fCentralityMethod == "ZNA_MANUAL") {
      Bool_t zna = kFALSE;
      for (Int_t j = 0; j < 4; ++j) {
        if (fESD->GetZDCData()->GetZDCTDCData(12,j) != 0) {
          zna = kTRUE;
        }
      }

//       Printf("%d %f", zna, fZNAtower[0]);
      if (zna) {
        // code from Chiara O (23.10.12)
        const Double_t *fZNAtower = fESD->GetZDCData()->GetZN2TowerEnergy();
        Float_t znacut[4] = {681., 563., 413., 191.};

        if (fZNAtower[0]>znacut[0])
          centrality = 1;
        else if (fZNAtower[0]>znacut[1])
          centrality = 21;
        else if (fZNAtower[0]>znacut[2])
          centrality = 41;
        else if (fZNAtower[0]>znacut[3])
          centrality = 61;
        else
          centrality = 81;
      }
      else
        centrality = -1;
    }
    else if (fCentralityMethod == "ZNAC") { // pp
      // values from Cvetan
      const Double_t *towZNA = fAOD->GetZDCData()->GetZNATowerEnergy();
      const Double_t *towZNC = fAOD->GetZDCData()->GetZNCTowerEnergy();

      Double_t enZNA = 1.13 * towZNA[0];
      Double_t enZNC = towZNC[0];
      Double_t enZN = enZNA + enZNC;

      Float_t znaccut[8] = {1e9, 797.743, 526.879, 238.595, 107.325, 39.214, 7.633, -100};

      if (enZN > znaccut[1] && enZN < znaccut[0])
        centrality = 96;
      else if (enZN > znaccut[2])
        centrality = 91;
      else if (enZN > znaccut[3])
        centrality = 81;
      else if (enZN > znaccut[4])
        centrality = 71;
      else if (enZN > znaccut[5])
        centrality = 61;
      else if (enZN > znaccut[6])
        centrality = 51;
      else if (enZN > znaccut[7])
        centrality = 1;
      else
        centrality = -1;

      ((TH1D*) fListOfHistos->FindObject("ZNA+C_energy"))->Fill(enZN);
    }
    else if (fCentralityMethod == "TRACKS_MANUAL") {
      // for pp
      TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, -1, kTRUE);
      centrality = tracks->GetEntriesFast();
      if (centrality > 40)
        centrality = 41;
//       Printf("%d %f", tracks->GetEntriesFast(), centrality);
      delete tracks;
    }
    else if (fCentralityMethod == "V0A_MANUAL") {
      // for pp

      //Total multiplicity in the VZERO A detector
      Float_t MV0A=inputEvent->GetVZEROData()->GetMTotV0A();
      Float_t MV0AScaled=0.;
      if (fMap) {
        TParameter<float>* sf=(TParameter<float>*)fMap->GetValue(Form("%d",inputEvent->GetRunNumber()));
        if (sf)
          MV0AScaled=MV0A*sf->GetVal();
      }

      if (MV0AScaled > 0)
        centrality = MV0AScaled;
      else
        centrality = -1;
    }
    else if (fCentralityMethod == "nano") {
      centrality = ((AliNanoAODHeader*) fAOD->GetHeader())->GetCentrality();
    }
    else if (fCentralityMethod.BeginsWith("Nano.")) {
      AliNanoAODHeader* nanoHeader = dynamic_cast<AliNanoAODHeader*>(fAOD->GetHeader());
      if (!nanoHeader)
        AliFatal("Nano Header not found");

      static TString nanoField = fCentralityMethod(5, fCentralityMethod.Length());
      static const Int_t kField = nanoHeader->GetVarIndex(nanoField);

      centrality = nanoHeader->GetVar(kField);
    }
    else if (fCentralityMethod == "PPVsMultUtils") {
      if (fAnalysisUtils) centrality = fAnalysisUtils->GetMultiplicityPercentile(inputEvent);
      else centrality = -1;
    }
    else if (fCentralityMethod == "MC_b") {
      AliGenEventHeader* eventHeader = GetFirstHeader();
      if (!eventHeader) {
        // We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
        // (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
        AliError("Event header not found. Skipping this event.");
        fHistos->FillEvent(0, AliUEHist::kCFStepAnaTopology);
        return -1;
      }

      AliCollisionGeometry* collGeometry = dynamic_cast<AliCollisionGeometry*> (eventHeader);
      AliGenHepMCEventHeader* hepMCHeader = dynamic_cast<AliGenHepMCEventHeader*> (eventHeader);
      if (!collGeometry && !hepMCHeader) {
        eventHeader->Dump();
        AliFatal("Asking for MC_b centrality, but event header has no collision geometry information");
      }

      if (collGeometry)
        centrality = collGeometry->ImpactParameter();
      else if (hepMCHeader)
        centrality = hepMCHeader->impact_parameter();
    }
    else if (fCentralityMethod == "MCGen_V0M") {
//      TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, -1, kFALSE, kFALSE, -999.,kTRUE);
      TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kFALSE, -1, kFALSE, kFALSE, -999.,kTRUE);
      Float_t MultV0M=0.;
      Float_t dNchdeta = 0.;
      Float_t INEL0 = 0.;
      for (Int_t i=0; i<tmpList->GetEntriesFast(); i++) {
        AliMCParticle* particle = dynamic_cast<AliMCParticle*> (tmpList->UncheckedAt(i));
        if (!particle)
          continue;
        Int_t pdgabs=TMath::Abs(particle->PdgCode());
        if (pdgabs==9902210)
          return -1; //no diffractive protons
        if (pdgabs!=211 && pdgabs!=321 && pdgabs!=2212)
          continue; //only charged pi+K+p
        if (!(((AliMCEvent*)mc)->IsPhysicalPrimary(particle->GetLabel())))
          continue; // only primaries (probably not necessary)
        if ( particle->Pt() < 0.001 || particle->Pt() > 50. )
          continue;
        Float_t eta=particle->Eta();
        if ((eta > 2.8 && eta < 5.1) || (eta > -3.7 && eta < -1.7))
          MultV0M += 1.;
        if (eta < 0.5 && eta > -0.5)
          dNchdeta += 1.;
        if (eta < 1. && eta > -1.)
          INEL0 += 1.;
      }
      if (INEL0 < 0.5)
        centrality = -1.; // INEL>0 cut
      else {
        ((TH2D*)fListOfHistos->FindObject("Mult_MCGen_V0M"))->Fill(MultV0M,dNchdeta);
        if (fCentralityMCGen_V0M)
          centrality = fCentralityMCGen_V0M->GetBinContent(fCentralityMCGen_V0M->GetXaxis()->FindBin(MultV0M));
        else
          centrality=-1.;
      }
    }
    else if (fCentralityMethod == "MCGen_CL1") {
//      TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, -1, kFALSE, kFALSE, -999.,kTRUE);
      TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kFALSE, -1, kFALSE, kFALSE, -999.,kTRUE);
      Float_t MultCL1=0.;
      Float_t dNchdeta = 0.;
      Float_t INEL0 = 0.;
      for (Int_t i=0; i<tmpList->GetEntriesFast(); i++) {
        AliMCParticle* particle = dynamic_cast<AliMCParticle*> (tmpList->UncheckedAt(i));
        if (!particle)
          continue;
        Int_t pdgabs = TMath::Abs(particle->PdgCode());
        if (pdgabs==9902210)
          return -1; //no diffractive protons
        if (pdgabs!=211 && pdgabs!=321 && pdgabs!=2212)
          continue; //only charged pi+K+p
        if (!(((AliMCEvent*)mc)->IsPhysicalPrimary(particle->GetLabel())))
          continue; // only primaries (probably not necessary)
        if ( particle->Pt() < 0.001 || particle->Pt() > 50. )
          continue;
        Float_t eta = particle->Eta();
        if (eta < 1.4 && eta > -1.4)
          MultCL1 += 1.;
        if (eta < 0.5 && eta > -0.5)
          dNchdeta += 1.;
        if (eta < 1. && eta > -1.)
          INEL0 += 1.;
      }
      if (INEL0 < 0.5)
        centrality = -1.; // INEL>0 cut
      else {
        ((TH2D*)fListOfHistos->FindObject("Mult_MCGen_CL1"))->Fill(MultCL1,dNchdeta);
        if (fCentralityMCGen_CL1)
          centrality = fCentralityMCGen_CL1->GetBinContent(fCentralityMCGen_CL1->GetXaxis()->FindBin(MultCL1));
        else
          centrality=-1.;
      }
    }
    else {
      if (fAOD)
        centralityObj = fAOD->GetCentrality();
      else if (fESD)
        centralityObj = fESD->GetCentrality();

      if (centralityObj) {
        if (fUseUncheckedCentrality)
          centrality = centralityObj->GetCentralityPercentileUnchecked(fCentralityMethod);
        else
          centrality = centralityObj->GetCentralityPercentile(fCentralityMethod);
      }
      else
        centrality = -1;

      if (fAOD) {
        // remove outliers
        if (centrality == 0) {
          if (fAOD->GetVZEROData()) {
            Float_t multV0 = 0;
            for (Int_t i=0; i<64; i++)
              multV0 += fAOD->GetVZEROData()->GetMultiplicity(i);
            if (multV0 < 19500) {
              centrality = -1;
              AliInfo("Rejecting event due to too small V0 multiplicity");
            }
          }
        }
      }
    }
  }
  AliInfo(Form("Centrality is %f", centrality));

  return centrality;
}

//____________________________________________________________________
TObjArray* AliAnalysisTaskPhiCorrelations::CloneAndReduceTrackList(TObjArray* tracks, Double_t minPt, Double_t maxPt)
{
  // clones a track list by using AliBasicParticle which uses much less memory (used for event mixing)
  // Clone only a certain pt bin on demand

  // Check if we already have a reduced track list. In that case a simple Clone is enough
  // Only possible for inclusive pt
  if (maxPt-minPt < 0)
    if (tracks->GetEntriesFast() == 0 || tracks->UncheckedAt(0)->InheritsFrom("AliBasicParticle"))
      return (TObjArray*) tracks->Clone();

  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
    AliVParticle* particle = (AliVParticle*) tracks->UncheckedAt(i);
    AliBasicParticle* copy = 0;

    if ( (maxPt-minPt > 0) && ((particle->Pt()<minPt) || (particle->Pt()>=maxPt)) )
      continue;

    if (fFillCorrelationsRapidity)
      copy = new AliBasicParticle(particle->Y(), particle->Phi(), particle->Pt(), particle->Charge());
    else
      copy = new AliBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
    copy->SetUniqueID(particle->GetUniqueID());

    // Set unique event index if input tracks are AliBasicParticles
    AliBasicParticle* particleBasic = dynamic_cast<AliBasicParticle*>(particle);
    if (particleBasic)
      copy->SetEventIndex(particleBasic->GetEventIndex());

    tracksClone->Add(copy);
  }

  return tracksClone;
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::Initialize()
{
  // input handler
  fInputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  // MC handler
  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::RemoveDuplicates(TObjArray* tracks)
{
  // remove particles with the same label

  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i) {
    AliVParticle* part = (AliVParticle*) tracks->At(i);

    for (Int_t j=i+1; j<before; ++j) {
      AliVParticle* part2 = (AliVParticle*) tracks->At(j);

      if (part->GetLabel() == part2->GetLabel()) {
        Printf("Removing %d with label %d (duplicated in %d)", i, part->GetLabel(), j); part->Dump(); part2->Dump();
        TObject* object = tracks->RemoveAt(i);
        if (tracks->IsOwner())
          delete object;
        break;
      }
    }
  }

  tracks->Compress();

  if (before > tracks->GetEntriesFast())
    AliInfo(Form("Reduced from %d to %d", before, tracks->GetEntriesFast()));
}

void AliAnalysisTaskPhiCorrelations::CleanUp(TObjArray* tracks, TObject* mcObj, Int_t maxLabel)
{
  // calls RemoveInjectedSignals, RemoveWeakDecays and RemoveDuplicates

  if (!tracks)
    return;

  if (fInjectedSignals)
    fAnalyseUE->RemoveInjectedSignals(tracks, mcObj, maxLabel);
  if (fRemoveWeakDecays)
    fAnalyseUE->RemoveWeakDecays(tracks, mcObj);
  if (fRemoveDuplicates)
    RemoveDuplicates(tracks);
  if (fRemoveWeakDecaysInMC)
    RemoveWeakDecaysInMC(tracks, mcObj);
}

void AliAnalysisTaskPhiCorrelations::RemoveWeakDecaysInMC(TObjArray* tracks, TObject* mcObj)
{
  // Exclude weak decay products (if not done by IsPhysicalPrimary)
  // In order to prevent analyzing daughters from weak decays
  // - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it

  const Int_t kNWeakParticles = 7;
  const Int_t kWeakParticles[kNWeakParticles] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
                                                  3122, 3112, // Lambda0 Sigma+-
                                                  130, 310 }; // K_L0 K_S0

  AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(mcObj);
  if (!mcEvent)
    return;

  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i) {
    AliMCParticle* mcParticle = dynamic_cast<AliMCParticle*> (tracks->UncheckedAt(i));
    if (!mcParticle)
      continue;

    TParticle *particle = mcParticle->Particle();
    if (!particle)
      continue;

    Int_t motherIndex = particle->GetMother(0);
    if (motherIndex == -1)
      continue;

    AliMCParticle* motherTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(motherIndex));
    if (!motherTrack)
      continue;

    TParticle *motherParticle = motherTrack->Particle();
    if (!motherParticle)
      continue;

    Int_t pdgcode = TMath::Abs(motherParticle->GetPdgCode());

    for (Int_t j=0; j != kNWeakParticles; ++j) {
      if (kWeakParticles[j] == pdgcode) {
        AliDebug(1, Form("Removing particle %d (pdg code %d; mother %d)", i, particle->GetPdgCode(), motherParticle->GetPdgCode()));
        TObject* object = tracks->RemoveAt(i);
        if (tracks->IsOwner())
          delete object;
        break;
      }
    }
  }

  tracks->Compress();
  if (before > tracks->GetEntriesFast())
    AliInfo(Form("Reduced from %d to %d", before, tracks->GetEntriesFast()));
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::SelectCharge(TObjArray* tracks)
{
  // remove particles with charge not selected (depending on fTriggerSelectCharge)

  if (!tracks)
    return;

  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i) {
    AliVParticle* part = (AliVParticle*) tracks->At(i);

    if (part->Charge() * fTriggerSelectCharge < -1) {
//       Printf("Removing %d with charge %d", i, part->Charge());
      TObject* object = tracks->RemoveAt(i);
      if (tracks->IsOwner())
        delete object;
    }
  }

  tracks->Compress();

  if (before > tracks->GetEntriesFast())
    AliInfo(Form("Reduced from %d to %d", before, tracks->GetEntriesFast()));
}

//____________________________________________________________________
Bool_t AliAnalysisTaskPhiCorrelations::AcceptEventCentralityWeight(Double_t centrality)
{
  // rejects "randomly" events such that the centrality gets flat
  // uses fCentralityWeights histogram

  // TODO code taken and adapted from AliRDHFCuts; waiting for general class AliCentralityFlattening

  Double_t weight = fCentralityWeights->GetBinContent(fCentralityWeights->FindBin(centrality));
  Double_t centralityDigits = centrality*100. - (Int_t)(centrality*100.);

  Bool_t result = kFALSE;
  if (centralityDigits < weight)
    result = kTRUE;

  AliInfo(Form("Centrality: %f; Digits: %f; Weight: %f; Result: %d", centrality, centralityDigits, weight, result));

  return result;
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::ShiftTracks(TObjArray* tracks, Double_t angle)
{
  // shifts the phi angle of all tracks by angle
  // 0 <= angle <= 2pi

  for (Int_t i=0; i<tracks->GetEntriesFast(); ++i) {
    AliBasicParticle* part = (AliBasicParticle*) tracks->At(i);
    Double_t newAngle = part->Phi() + angle;
    if (newAngle >= TMath::TwoPi())
      newAngle -= TMath::TwoPi();

    part->SetPhi(newAngle);
  }
}

//____________________________________________________________________
TObjArray* AliAnalysisTaskPhiCorrelations::GetParticlesFromDetector(AliVEvent* inputEvent, Int_t idet)
{
  //1 = VZERO_A; 2 = VZERO_C; 3 = SPD tracklets;
  //4 = muon tracks; 5 = global tracks w/o jets
  //6 = custom particles A; 7 = custom particles B
  TObjArray* obj = new TObjArray;
  obj->SetOwner(kTRUE);

  if (idet == 0 || idet > 9)
    AliFatal("Value of idet causes problem with uniqueID to avoid double counting");

  if (idet == 1 || idet == 2) {
    AliVVZERO* vZero = inputEvent->GetVZEROData();

    const Int_t vZeroStart = (idet == 1) ? 32 : 0;

    TH1F* singleCells = (TH1F*) fListOfHistos->FindObject("V0SingleCells");
    for (Int_t i=vZeroStart; i<vZeroStart+32; i++) {
      Float_t weight = vZero->GetMultiplicity(i);
      singleCells->Fill(weight);

      // rough estimate of multiplicity
      for (Int_t j=0; j<TMath::Nint(weight); j++) {
        AliBasicParticle* particle = new AliBasicParticle((AliVVZERO::GetVZEROEtaMax(i) + AliVVZERO::GetVZEROEtaMin(i)) / 2, AliVVZERO::GetVZEROAvgPhi(i), 1.1, 0); // fix pT = 1.1 and charge = 0
        particle->SetUniqueID((fAnalyseUE->GetEventCounter() * 50000 + j + i * 1000) * 10 + idet);
        particle->SetEventIndex(GetUniqueEventID(inputEvent));

        obj->Add(particle);
      }
    }
  }
  else if (idet == 3) {
    if (!fAOD)
      AliFatal("Tracklets only available on AOD");
    AliAODTracklets* trklets=(AliAODTracklets*)fAOD->GetTracklets();
    if (!trklets)
      AliFatal("AliAODTracklets not found");
    for (Int_t itrklets=0;itrklets<trklets->GetNumberOfTracklets();itrklets++) {
      Float_t eta=-TMath::Log(TMath::Tan(trklets->GetTheta(itrklets)/2));
      if (TMath::Abs(eta)>fTrackEtaCut)
        continue;
      Float_t pT=1000*TMath::Abs(trklets->GetDeltaPhi(itrklets));//in mrad
      if (pT>fTrackletDphiCut)
        continue;
      TH1F* DphiTrklets = (TH1F*)fListOfHistos->FindObject("DphiTrklets");
      DphiTrklets->Fill(1000*trklets->GetDeltaPhi(itrklets)); //in mrad
      Float_t phi=trklets->GetPhi(itrklets);
      phi+=trklets->GetDeltaPhi(itrklets)*39./34.; //correction dphi*39./34. (Dphi in rad)
      if (phi<0)
        phi+=TMath::TwoPi();
      if (phi>TMath::TwoPi())
        phi-=TMath::TwoPi();

      AliBasicParticle* particle = new AliBasicParticle(eta,phi, pT, 0); // pT = TMath::Abs(trklets->GetDeltaPhi(itrklets)) in mrad and charge = 0
      particle->SetUniqueID((fAnalyseUE->GetEventCounter() * 50000 + itrklets) * 10 + idet);
      particle->SetEventIndex(GetUniqueEventID(inputEvent));

      obj->Add(particle);
    }
  }
  else if (idet == 4) {
    if (!fAOD)
      AliFatal("Muon selection only implemented on AOD");//FIXME to be implemented also for ESDs as in AliAnalyseLeadingTrackUE::GetAcceptedPArticles
    for (Int_t iTrack = 0; iTrack < fAOD->GetNumberOfTracks(); iTrack++) {
      AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
      if (!track)
        AliFatal("Not a standard AOD");
      if (!track->IsMuonTrack())
        continue;
      //Float_t dca    = track->DCA();
      //Float_t chi2   = track->Chi2perNDF();
      Float_t rabs   = track->GetRAtAbsorberEnd();
      Float_t eta    = track->Eta();
      Int_t   matching   = track->GetMatchTrigger();
      if (rabs < 17.6 || rabs > 89.5)
        continue;
      if (eta < -4 || eta > -2.5)
        continue;
      if (matching < 2)
        continue;

      AliBasicParticle* particle = new AliBasicParticle(eta,track->Phi(), track->Pt(), track->Charge()); 
      particle->SetUniqueID((fAnalyseUE->GetEventCounter() * 50000 + iTrack) * 10 + idet);
      particle->SetEventIndex(GetUniqueEventID(inputEvent));

      obj->Add(particle);
    }
  }
  else if (idet == 5) {
    if (!fAOD)
      AliFatal("Cannot access jets without AOD");

    // retrieve jet array
    TClonesArray *jetArray = 0x0;
    if (fJetBranchName.Length() > 0) {
      jetArray = dynamic_cast<TClonesArray*> (fAOD->FindListObject(fJetBranchName));
      if (!jetArray) {
        fAOD->Print();
        AliFatal(Form("Cannot access jet branch <%s>", fJetBranchName.Data()));
      }
    }

    const Int_t nJets = jetArray ? jetArray->GetEntries() : 0;
    const Int_t nTracks = fAOD->GetNumberOfTracks();

    for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
      AliAODTrack *track = dynamic_cast<AliAODTrack*> (fAOD->GetTrack(iTrack));

      if (!track)
        AliFatal("Not a standard AOD");

      // track cuts
      if (!track->IsHybridGlobalConstrainedGlobal() ||
          (TMath::Abs(track->Eta()) > fTrackEtaMax))
        continue;

      // check if track is part of a jet
      Bool_t rejectTrack = kFALSE;
      for (Int_t iJet = 0; iJet < nJets; ++iJet) {
        AliEmcalJet *jet = dynamic_cast<AliEmcalJet*> (jetArray->At(iJet));

        // skip jets outside of acceptance and too low pt
        if ((TMath::Abs(jet->Eta()) > fJetEtaMax) ||
            (jet->Pt() < fJetPtMin) ||
            (jet->GetNumberOfConstituents() < fJetConstMin))
          continue;

        Float_t dEta = track->Eta() - jet->Eta();
        Float_t dPhi = TVector2::Phi_mpi_pi(track->Phi() - jet->Phi());
        Bool_t trackInCone = (TMath::Power(dPhi, 2.) + TMath::Power(dEta, 2.)) < TMath::Power(fExclusionRadius, 2.);
        Bool_t trackInJet  = jet->ContainsTrack(track, fAOD->GetTracks()) >= 0;

        if (trackInJet != trackInCone)
          AliDebug(2, Form("track %15s - %15s",
                           trackInJet ? "in jet" : "not in jet",
                           trackInCone ? "in cone" : "not in cone"));

        // exclude track if in cone around jet or assigned to jet
        rejectTrack = fExclusionRadius > 0. ? trackInCone : trackInJet;

        if (rejectTrack)
          break;
      }

      // skip track if it is associated to a selected jet
      // (either by jet finder or by cone radius)
      if (rejectTrack)
        continue;

      // add particle to array
      AliBasicParticle* particle =
        new AliBasicParticle(track->Eta(), track->Phi(), track->Pt(), track->Charge());
      // NOTE "+ idet" is missing here on purpose as the tracks are the same as in the case of idet == 0
      particle->SetUniqueID((fAnalyseUE->GetEventCounter() * 50000 + iTrack) * 10);
      particle->SetEventIndex(GetUniqueEventID(inputEvent));
      obj->Add(particle);
    }
  }
  else if ( (idet == 6) || (idet == 7) ) { // custom particle arrays
    TClonesArray* clonesArray = 0;
    if (idet == 6) {
      if (fCustomParticlesA=="")
        AliFatal("Custom particle array A demanded. Use SetCustomParticlesA() to give particle array's name");

      // Get custom particle array
      clonesArray = dynamic_cast<TClonesArray*> (inputEvent->FindListObject(fCustomParticlesA.Data()));
      if (!clonesArray) {
        inputEvent->Print();
        AliFatal(Form("Cannot find custom particle array A <%s>", fCustomParticlesA.Data()));
      }
    }
    else if (idet == 7) {
      if (fCustomParticlesB=="")
        AliFatal("Custom particle array B demanded. Use SetCustomParticlesB() to give particle array's name");

      // Get custom particle array
      clonesArray = dynamic_cast<TClonesArray*> (inputEvent->FindListObject(fCustomParticlesB.Data()));
      if (!clonesArray) {
        inputEvent->Print();
        AliFatal(Form("Cannot find custom particle array B <%s>", fCustomParticlesB.Data()));
      }
    }

    for (Int_t iParticle = 0; iParticle < clonesArray->GetEntries(); iParticle++) {
      AliVParticle* customParticle = dynamic_cast<AliVParticle*> (clonesArray->At(iParticle));

      if (!customParticle) {
        AliError(Form("Custom particle in array %s cannot be casted to AliVParticle", ( (idet==6) ? "A" : "B" ) ));
        continue;
      }

      // add particle to object array
      AliBasicParticle* particle = new AliBasicParticle(customParticle->Eta(), customParticle->Phi(), customParticle->Pt(), customParticle->Charge());

      // Set custom ID. Should be the same if custom particle A and B are the same
      if (fCustomParticlesA == fCustomParticlesB)
        particle->SetUniqueID((fAnalyseUE->GetEventCounter() * 50000 + iParticle) * 10 + 6);
      else
        particle->SetUniqueID((fAnalyseUE->GetEventCounter() * 50000 + iParticle) * 10 + idet);

      particle->SetEventIndex(GetUniqueEventID(inputEvent));

      obj->Add(particle);
    }
  }
  else
    AliFatal(Form("GetParticlesFromDetector: Invalid idet value: %d", idet));

  return obj;
}

//____________________________________________________________________
Bool_t AliAnalysisTaskPhiCorrelations::IsMuEvent()
{
  if (!fAOD)
    AliFatal("Muon selection only implemented on AOD");//FIXME to be implemented also for ESDs as in AliAnalyseLeadingTrackUE::GetAcceptedPArticles
  for (Int_t iTrack = 0; iTrack < fAOD->GetNumberOfTracks(); iTrack++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if (!track)
      AliFatal("Not a standard AOD");
    if (!track->IsMuonTrack())
      continue;
    //Float_t dca    = track->DCA();
    //Float_t chi2   = track->Chi2perNDF();
    Float_t rabs   = track->GetRAtAbsorberEnd();
    Float_t eta    = track->Eta();
    Int_t   matching   = track->GetMatchTrigger();
    if (rabs < 17.6 || rabs > 89.5)
      continue;
    if (eta < -4 || eta > -2.5)
      continue;
    if (matching < 2)
      continue;
    return kTRUE;
  }
  return kFALSE;

}

//____________________________________________________________________
Bool_t AliAnalysisTaskPhiCorrelations::InitiateEventPlane(Double_t& evtPlanePhi, AliVEvent* inputEvent)
{
  AliEventplane* evtPlane = inputEvent->GetEventplane();
  Double_t qx = 0; Double_t qy = 0;
  if (evtPlane) {
    evtPlanePhi = evtPlane->CalculateVZEROEventPlane(inputEvent, 10, 2, qx, qy);
    return 1;
  }
  else
    return 0;
}

//____________________________________________________________________
Long64_t AliAnalysisTaskPhiCorrelations::GetUniqueEventID(AliVEvent* inputEvent)
{
  // Get event ID from header
  AliVHeader* eventIDHeader = inputEvent->GetHeader();
  if (eventIDHeader)
    return eventIDHeader->GetEventIdAsLong();
  else
    return 0;
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt)
{
  std::vector<Double_t> binVec;
  binVec.push_back(minCent);
  binVec.push_back(maxCent);
  binVec.push_back(minZvtx);
  binVec.push_back(maxZvtx);
  binVec.push_back(minPt);
  binVec.push_back(maxPt);
  fEventPoolOutputList.push_back(binVec);
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::FinishTaskOutput()
{
  // Clear unnecessary pools before saving
  fPoolMgr->ClearPools();
}
