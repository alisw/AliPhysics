/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables
// 6 Steps introduced: MC, MC Acc, Reco, Reco Acc, Reco Acc + ITS Cl,
// Reco Acc + ITS Cl + PPR cuts
// 12 variables used: pt, y, cosThetaStar, ptPi, ptK, ct,
// dca, d0Pi, d0K, d0Pixd0K, cosPointingAngle, phi
//
//-----------------------------------------------------------------------
// Author : C. Zampolli, CERN
//          D. Caffarri, Univ & INFN Padova  caffarri@pd.infn.it
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Base class for HF Unfolding (pt and eta)
// correlation matrix filled at Acceptance and PPR level
// Author: A.Grelli ,  Utrecht - agrelli@uu.nl
//-----------------------------------------------------------------------
#include <TCanvas.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TProfile.h>
#include <TH1I.h>
#include <TStyle.h>
#include <TFile.h>
#include <TF1.h>

#include "AliCFTaskVertexingHF.h"
#include "AliMCEvent.h"
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliLog.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TH2D.h"
#include "AliESDtrackCuts.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliCFVertexingHF2Prong.h"
#include "AliCFVertexingHF3Prong.h"
#include "AliCFVertexingHFCascade.h"
#include "AliCFVertexingHFLctoV0bachelor.h"
#include "AliCFVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisVertexingHF.h"
#include "AliPIDResponse.h"

//__________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF() :
  AliAnalysisTaskSE(),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
  fListProfiles(0),
  fCountMC(0),
  fCountGenLimAcc(0),
  fCountGenLimAccNoAcc(0),
  fCountAcc(0),
  fCountVertex(0),
  fCountRefit(0),
  fCountReco(0),
  fCountRecoAcc(0),
  fCountRecoITSClusters(0),
  fCountRecoPPR(0),
  fCountRecoPID(0),
  fEvents(0),
  fDecayChannel(0),
  fFillFromGenerated(kFALSE),
  fOriginDselection(0),
  fAcceptanceUnf(kTRUE),
  fCuts(0),
  fUseWeight(kFALSE),
  fWeight(1.),
  fUseFlatPtWeight(kFALSE),
  fUseZWeight(kFALSE),
  fUseNchWeight(kFALSE),
  fUseTrackletsWeight(kFALSE),
  fUseMultRatioAsWeight(kFALSE),
  fNvar(0),
  fPartName(""),
  fDauNames(""),
  fSign(2),
  fCentralitySelection(kTRUE),
  fFakeSelection(0),
  fRejectIfNoQuark(kTRUE),
  fUseMCVertex(kFALSE),
  fDsOption(1),
  fGenDsOption(3),
  fConfiguration(kCheetah), // by default, setting the fast configuration
  fFuncWeight(0x0),
  fHistoPtWeight(0x0),
  fHistoMeasNch(0x0),
  fHistoMCNch(0x0),
  fResonantDecay(0),
  fLctoV0bachelorOption(1),
  fGenLctoV0bachelorOption(0),
  fUseSelectionBit(kTRUE),
  fPDGcode(0),
  fMultiplicityEstimator(kNtrk10),
  fRefMult(9.26),
  fZvtxCorrectedNtrkEstimator(kFALSE),
  fIsPPData(kFALSE),
  fIsPPbData(kFALSE),
  fUseAdditionalCuts(kFALSE),
  fUseCutsForTMVA(kFALSE),
  fUseCascadeTaskForLctoV0bachelor(kFALSE),
  fCutOnMomConservation(0.00001)
{
  //
  //Default ctor
  //
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}
//___________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF(const Char_t* name, AliRDHFCuts* cuts, TF1* func) :
  AliAnalysisTaskSE(name),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
  fListProfiles(0),
  fCountMC(0),
  fCountGenLimAcc(0),
  fCountGenLimAccNoAcc(0),
  fCountAcc(0),
  fCountVertex(0),
  fCountRefit(0),
  fCountReco(0),
  fCountRecoAcc(0),
  fCountRecoITSClusters(0),
  fCountRecoPPR(0),
  fCountRecoPID(0),
  fEvents(0),
  fDecayChannel(0),
  fFillFromGenerated(kFALSE),
  fOriginDselection(0),
  fAcceptanceUnf(kTRUE),
  fCuts(cuts),
  fUseWeight(kFALSE),
  fWeight(1.),
  fUseFlatPtWeight(kFALSE),
  fUseZWeight(kFALSE),
  fUseNchWeight(kFALSE),
  fUseTrackletsWeight(kFALSE),
  fUseMultRatioAsWeight(kFALSE),
  fNvar(0),
  fPartName(""),
  fDauNames(""),
  fSign(2),
  fCentralitySelection(kTRUE),
  fFakeSelection(0),
  fRejectIfNoQuark(kTRUE),
  fUseMCVertex(kFALSE),
  fDsOption(1),
  fGenDsOption(3),
  fConfiguration(kCheetah),  // by default, setting the fast configuration
  fFuncWeight(func),
  fHistoPtWeight(0x0),
  fHistoMeasNch(0x0),
  fHistoMCNch(0x0),
  fResonantDecay(0),
  fLctoV0bachelorOption(1),
  fGenLctoV0bachelorOption(0),
  fUseSelectionBit(kTRUE),
  fPDGcode(0),
  fMultiplicityEstimator(kNtrk10),
  fRefMult(9.26),
  fZvtxCorrectedNtrkEstimator(kFALSE),
  fIsPPData(kFALSE),
  fIsPPbData(kFALSE),
  fUseAdditionalCuts(kFALSE),
  fUseCutsForTMVA(kFALSE),
  fUseCascadeTaskForLctoV0bachelor(kFALSE),
  fCutOnMomConservation(0.00001)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());
  DefineOutput(3,THnSparseD::Class());
  DefineOutput(4,AliRDHFCuts::Class());
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
  DefineOutput(5,TList::Class()); // slot #5 keeps the zvtx Ntrakclets correction profiles

  fCuts->PrintAll();
}

//___________________________________________________________________________
AliCFTaskVertexingHF& AliCFTaskVertexingHF::operator=(const AliCFTaskVertexingHF& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fCFManager  = c.fCFManager;
    fHistEventsProcessed = c.fHistEventsProcessed;
    fCuts = c.fCuts;
    fFuncWeight = c.fFuncWeight;
    fHistoPtWeight = c.fHistoPtWeight;
    fHistoMeasNch = c.fHistoMeasNch;
    fHistoMCNch = c.fHistoMCNch;
    for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=c.fMultEstimatorAvg[i];
  }
  return *this;
}

//___________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF(const AliCFTaskVertexingHF& c) :
  AliAnalysisTaskSE(c),
  fCFManager(c.fCFManager),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fCorrelation(c.fCorrelation),
  fListProfiles(c.fListProfiles),
  fCountMC(c.fCountMC),
  fCountGenLimAcc(c.fCountGenLimAcc),
  fCountGenLimAccNoAcc(c.fCountGenLimAccNoAcc),
  fCountAcc(c.fCountAcc),
  fCountVertex(c.fCountVertex),
  fCountRefit(c.fCountRefit),
  fCountReco(c.fCountReco),
  fCountRecoAcc(c.fCountRecoAcc),
  fCountRecoITSClusters(c.fCountRecoITSClusters),
  fCountRecoPPR(c.fCountRecoPPR),
  fCountRecoPID(c.fCountRecoPID),
  fEvents(c.fEvents),
  fDecayChannel(c.fDecayChannel),
  fFillFromGenerated(c.fFillFromGenerated),
  fOriginDselection(c.fOriginDselection),
  fAcceptanceUnf(c.fAcceptanceUnf),
  fCuts(c.fCuts),
  fUseWeight(c.fUseWeight),
  fWeight(c.fWeight),
  fUseFlatPtWeight(c.fUseFlatPtWeight),
  fUseZWeight(c.fUseZWeight),
  fUseNchWeight(c.fUseNchWeight),
  fUseTrackletsWeight(c.fUseTrackletsWeight),
  fUseMultRatioAsWeight(c.fUseMultRatioAsWeight),
  fNvar(c.fNvar),
  fPartName(c.fPartName),
  fDauNames(c.fDauNames),
  fSign(c.fSign),
  fCentralitySelection(c.fCentralitySelection),
  fFakeSelection(c.fFakeSelection),
  fRejectIfNoQuark(c.fRejectIfNoQuark),
  fUseMCVertex(c.fUseMCVertex),
  fDsOption(c.fDsOption),
  fGenDsOption(c.fGenDsOption),
  fConfiguration(c.fConfiguration),
  fFuncWeight(c.fFuncWeight),
  fHistoPtWeight(c.fHistoPtWeight),
  fHistoMeasNch(c.fHistoMeasNch),
  fHistoMCNch(c.fHistoMCNch),
  fResonantDecay(c.fResonantDecay),
  fLctoV0bachelorOption(c.fLctoV0bachelorOption),
  fGenLctoV0bachelorOption(c.fGenLctoV0bachelorOption),
  fUseSelectionBit(c.fUseSelectionBit),
  fPDGcode(c.fPDGcode),
  fMultiplicityEstimator(c.fMultiplicityEstimator),
  fRefMult(c.fRefMult),
  fZvtxCorrectedNtrkEstimator(c.fZvtxCorrectedNtrkEstimator),
  fIsPPData(c.fIsPPData),
  fIsPPbData(c.fIsPPbData),
  fUseAdditionalCuts(c.fUseAdditionalCuts),
  fUseCutsForTMVA(c.fUseCutsForTMVA),
  fUseCascadeTaskForLctoV0bachelor(c.fUseCascadeTaskForLctoV0bachelor),
  fCutOnMomConservation(c.fCutOnMomConservation)
{
  //
  // Copy Constructor
  //
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=c.fMultEstimatorAvg[i];
}

//___________________________________________________________________________
AliCFTaskVertexingHF::~AliCFTaskVertexingHF()
{
  //
  //destructor
  //
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fCorrelation)         delete fCorrelation ;
  if (fListProfiles)        delete fListProfiles;
  if (fCuts)                delete fCuts;
  if (fFuncWeight)          delete fFuncWeight;
  if (fHistoPtWeight)       delete fHistoPtWeight;
  if (fHistoMeasNch)        delete fHistoMeasNch;
  if (fHistoMCNch)          delete fHistoMCNch;
  for(Int_t i=0; i<4; i++) { if(fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i]; }
}

//_________________________________________________________________________-
void AliCFTaskVertexingHF::Init()
{
  //
  // Initialization
  //

  if (fDebug>1) printf("AliCFTaskVertexingHF::Init()");
  if(fUseWeight && fUseZWeight) { AliFatal("Can not use at the same time pt and z-vtx weights, please choose"); return; }
  if(fUseWeight && fUseNchWeight) { AliInfo("Beware, using at the same time pt and Nch weights, please check"); }
  if(fUseNchWeight && !fHistoMCNch) { AliFatal("Need to pass the MC Nch distribution to use Nch weights"); return; }
  if(fUseNchWeight && !fHistoMeasNch) CreateMeasuredNchHisto();

  AliRDHFCuts *copyfCuts = 0x0;
  if (!fCuts){
    AliFatal("No cuts defined - Exiting...");
    return;
  }

  switch (fDecayChannel){
  case 2:{
    fPDGcode = 421;
    copyfCuts = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
    switch (fConfiguration) {
    case kSnail:  // slow configuration: all variables in
      fNvar = 16;
      break;
    case kCheetah:// fast configuration: only pt_candidate, y, phi, ct, fake, z_vtx, centrality, multiplicity will be filled
      fNvar = 8;
      break;
    }
    fPartName="D0";
    fDauNames="K+pi";
    break;
  }
  case 21:{
    fPDGcode = 413;
    copyfCuts = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
    switch (fConfiguration) {
    case kSnail:  // slow configuration: all variables in
      fNvar = 16;
      break;
    case kCheetah:// fast configuration: only pt_candidate, y, phi, ct, fake, z_vtx, centrality, multiplicity will be filled
      fNvar = 8;
      break;
    }
    fPartName="Dstar";
    fDauNames="K+pi+pi";
    break;
  }
  case 22:{
    fPDGcode = 4122;
    copyfCuts = new AliRDHFCutsLctoV0(*(static_cast<AliRDHFCutsLctoV0*>(fCuts)));
    switch (fConfiguration) {
    case kSnail:  // slow configuration: all variables in
      fNvar = 16;
      break;
    case kCheetah:// fast configuration: only pt_candidate, y, phi, ct, fake, z_vtx, centrality, multiplicity will be filled
      fNvar = 8;
      break;
    }
    fPartName="Lambdac";
    fDauNames="V0+bachelor";
    break;
  }
  case 31:{
    fPDGcode = 411;
    copyfCuts = new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fCuts)));
    switch (fConfiguration) {
    case kSnail:  // slow configuration: all variables in
      fNvar = 14;
      break;
    case kCheetah:// fast configuration: only pt_candidate, y, phi, ct, fake, z_vtx, centrality, multiplicity will be filled
      fNvar = 8;
      break;
    }
    fPartName="Dplus";
    fDauNames="K+pi+pi";
    break;
  }
  case 32:{
    fPDGcode = 4122;
    copyfCuts = new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fCuts)));
    switch (fConfiguration) {
    case kSnail:  // slow configuration: all variables in
      fNvar = 18;
      break;
    case kCheetah:// fast configuration: only pt_candidate, y, phi, ct, fake, z_vtx, centrality, multiplicity will be filled
      fNvar = 8;
      break;
    }
    fPartName="Lambdac";
    fDauNames="p+K+pi";
    break;
  }
  case 33:{
    fPDGcode = 431;
    copyfCuts = new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fCuts)));
    switch (fConfiguration) {
    case kSnail:  // slow configuration: all variables in
      fNvar = 14;
      break;
    case kCheetah:// fast configuration: only pt_candidate, y, phi, ct, fake, z_vtx, centrality, multiplicity will be filled
      fNvar = 8;
      break;
    }
    fPartName="Ds";
    fDauNames="K+K+pi";
    break;
  }
  case 4:{
    fPDGcode = 421;
    copyfCuts = new AliRDHFCutsD0toKpipipi(*(static_cast<AliRDHFCutsD0toKpipipi*>(fCuts)));
    switch (fConfiguration) {
    case kSnail:  // slow configuration: all variables in
      fNvar = 16;
      break;
    case kCheetah:// fast configuration: only pt_candidate, y, phi, ct, fake, z_vtx, centrality, multiplicity will be filled
      fNvar = 8;
      break;
    }
    fPartName="D0";
    fDauNames="K+pi+pi+pi";
    break;
  }
  default:
    AliFatal("The decay channel MUST be defined according to AliCFVertexing::DecayChannel - Exiting...");
    break;
  }

  const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
  if (copyfCuts){
    copyfCuts->SetName(nameoutput);

    //Post the data
    PostData(4, copyfCuts);
  }
  else{
    AliFatal("Failing initializing AliRDHFCuts object - Exiting...");
  }

  fListProfiles = new TList();
  fListProfiles->SetOwner();
  TString period[4];
  Int_t nProfiles=4;

  if (fIsPPbData) { //if pPb, use only two estimator histos
    period[0] = "LHC13b"; period[1] = "LHC13c";
    nProfiles = 2;
  } else {        // else assume pp (four histos for LHC10)
    period[0] = "LHC10b"; period[1] = "LHC10c"; period[2] = "LHC10d"; period[3] = "LHC10e";
    nProfiles = 4;
  }

  for(Int_t i=0; i<nProfiles; i++){
    if(fMultEstimatorAvg[i]){
      TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
      hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
      fListProfiles->Add(hprof);
    }
  }

  // Save also the weight functions or histograms
  if(fFuncWeight) fListProfiles->Add(fFuncWeight);
  if(fHistoPtWeight) fListProfiles->Add(fHistoPtWeight);
  if(fHistoMeasNch) fListProfiles->Add(fHistoMeasNch);
  if(fHistoMCNch) fListProfiles->Add(fHistoMCNch);

  PostData(5,fListProfiles);

  return;
}

//_________________________________________________
void AliCFTaskVertexingHF::UserExec(Option_t *)
{
  //
  // Main loop function
  //

  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fCorrelation) ;

  AliDebug(3,Form("*** Processing event %d\n", fEvents));

  if (fFillFromGenerated){
    AliWarning("Flag to fill container with generated value ON ---> dca, d0pi, d0K, d0xd0, cosPointingAngle will be set as dummy!");
  }

  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  TClonesArray *arrayBranch=0;

  if(!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();

      switch (fDecayChannel){
      case 2:{
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
        break;
      }
      case 21:{
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
        break;
      }
      case 22:{
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
        break;
      }
      case 31:
      case 32:
      case 33:{
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
        break;
      }
      case 4:{
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm4Prong");
        break;
      }
      default:
        break;
      }
    }
  }
  else {
    switch (fDecayChannel){
    case 2:{
      arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
      break;
    }
    case 21:{
      arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
      break;
    }
    case 22:{
      arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("CascadesHF");
      break;
    }
    case 31:
    case 32:
    case 33:{
      arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("Charm3Prong");
      break;
    }
    case 4:{
      arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("Charm4Prong");
      break;
    }
    default:
      break;
    }
  }

  AliAODVertex *aodVtx = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!aodVtx) {
    AliDebug(3, "The event was skipped due to missing vertex");
    return;
  }

  if (!arrayBranch) {
    AliError("Could not find array of HF vertices");
    return;
  }

  fEvents++;

  fCFManager->SetRecEventInfo(aodEvent);
  fCFManager->SetMCEventInfo(aodEvent);

  //******** DEFINE number of variables of the container***** for now set at 13, in the future in the config macro.

  //loop on the MC event

  TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) {
    AliError("Could not find Monte-Carlo in AOD");
    return;
  }
  Int_t icountMC = 0;
  Int_t icountGenLimAcc = 0;
  Int_t icountGenLimAccNoAcc = 0;
  Int_t icountAcc = 0;
  Int_t icountReco = 0;
  Int_t icountVertex = 0;
  Int_t icountRefit = 0;
  Int_t icountRecoAcc = 0;
  Int_t icountRecoITSClusters = 0;
  Int_t icountRecoPPR = 0;
  Int_t icountRecoPID = 0;
  Int_t cquarks = 0;

  AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  if (!mcHeader) {
    AliError("Could not find MC Header in AOD");
    return;
  }

  fHistEventsProcessed->Fill(0.5);

  Double_t* containerInput = new Double_t[fNvar];
  Double_t* containerInputMC = new Double_t[fNvar];


  AliCFVertexingHF* cfVtxHF=0x0;
  switch (fDecayChannel){
  case 2:{
    cfVtxHF = new AliCFVertexingHF2Prong(mcArray, fOriginDselection);
    break;
  }
  case 21:{
    cfVtxHF = new AliCFVertexingHFCascade(mcArray, fOriginDselection);
    ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGcascade(413);
    ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGbachelor(211);
    ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaugh(421);
    ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaughForMC(421);
    ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaughPositive(211);
    ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaughNegative(321);
    ((AliCFVertexingHFCascade*)cfVtxHF)->SetPrimaryVertex(aodVtx);
    break;
  }
  case 22:{
    // Lc ->  K0S+proton
    if (fUseCascadeTaskForLctoV0bachelor){
      cfVtxHF = new AliCFVertexingHFCascade(mcArray, fOriginDselection);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGcascade(4122);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGbachelor(2212);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaugh(310);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaughForMC(311);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaughPositive(211);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetPDGneutrDaughNegative(211);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetPrimaryVertex(aodVtx);
      ((AliCFVertexingHFCascade*)cfVtxHF)->SetCutOnMomConservation(fCutOnMomConservation);
      if (fUseAdditionalCuts) ((AliCFVertexingHFCascade*)cfVtxHF)->SetUseCutsForTMVA(fUseCutsForTMVA);
    }
    else {
      cfVtxHF = new AliCFVertexingHFLctoV0bachelor(mcArray, fOriginDselection,fGenLctoV0bachelorOption);
    }
    break;
  }
  case 31:
    //	case 32:
  case 33:{
    cfVtxHF = new AliCFVertexingHF3Prong(mcArray, fOriginDselection, fDecayChannel);
    if(fDecayChannel==33){
      ((AliCFVertexingHF3Prong*)cfVtxHF)->SetGeneratedDsOption(fGenDsOption);
    }
    break;
  }
  case 32:{
    cfVtxHF = new AliCFVertexingHF3Prong(mcArray, fOriginDselection, fDecayChannel,fResonantDecay);
  }
  case 4:{
    //cfVtxHF = new AliCFVertexingHF4Prong(mcArray, originDselection);  // not there yet
    break;
  }
  default:
    break;
  }
  if (!cfVtxHF){
    AliError("No AliCFVertexingHF initialized");
    delete[] containerInput;
    delete[] containerInputMC;
    return;
  }

  Double_t zPrimVertex = aodVtx ->GetZ();
  Double_t zMCVertex = mcHeader->GetVtxZ();
  Int_t runnumber = aodEvent->GetRunNumber();

  // Multiplicity definition with tracklets
  Double_t nTracklets = 0;
  Int_t nTrackletsEta10 = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));
  Int_t nTrackletsEta16 = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.6,1.6));
  nTracklets = (Double_t)nTrackletsEta10;
  if(fMultiplicityEstimator==kNtrk10to16) { nTracklets = (Double_t)(nTrackletsEta16 - nTrackletsEta10); }

  // Apply the Ntracklets z-vtx data driven correction
  if(fZvtxCorrectedNtrkEstimator) {
    TProfile* estimatorAvg = GetEstimatorHistogram(aodEvent);
    if(estimatorAvg) {
      Int_t nTrackletsEta10Corr = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nTrackletsEta10,zPrimVertex,fRefMult));
      Int_t nTrackletsEta16Corr = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nTrackletsEta16,zPrimVertex,fRefMult));
      nTracklets = (Double_t)nTrackletsEta10Corr;
      if(fMultiplicityEstimator==kNtrk10to16) { nTracklets = (Double_t)(nTrackletsEta16Corr - nTrackletsEta10Corr); }
    }
  }


  fWeight=1.;
  if(fUseZWeight) fWeight *= GetZWeight(zMCVertex,runnumber);
  if(fUseNchWeight){
    Int_t nChargedMCPhysicalPrimary=AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(mcArray,-1.0,1.0);
    if(!fUseTrackletsWeight) fWeight *= GetNchWeight(nChargedMCPhysicalPrimary);
    else fWeight *= GetNchWeight(static_cast<Int_t>(nTracklets));
    AliDebug(2,Form("Using Nch weights, Mult=%d Weight=%f\n",nChargedMCPhysicalPrimary,fWeight));
  }
  Double_t eventWeight=fWeight;

  if (TMath::Abs(zMCVertex) > fCuts->GetMaxVtxZ()){
    AliDebug(3,Form("z coordinate of MC vertex = %f, it was required to be within [-%f, +%f], skipping event", zMCVertex, fCuts->GetMaxVtxZ(), fCuts->GetMaxVtxZ()));
    delete[] containerInput;
    delete[] containerInputMC;
    delete cfVtxHF;
    return;
  }

  if(aodEvent->GetTriggerMask()==0 &&
     (runnumber>=195344 && runnumber<=195677)){
    AliDebug(3,"Event rejected because of null trigger mask");
    delete[] containerInput;
    delete[] containerInputMC;
    delete cfVtxHF;
    return;
  }

  AliESDtrackCuts** trackCuts = new AliESDtrackCuts*[cfVtxHF->GetNProngs()];
  if (fDecayChannel == 21){
    // for the D*, setting the third element of the array of the track cuts to those for the soft pion
    for (Int_t iProng = 0; iProng<cfVtxHF->GetNProngs()-1; iProng++){
      trackCuts[iProng]=fCuts->GetTrackCuts();
    }
    trackCuts[2] = fCuts->GetTrackCutsSoftPi();
  }
  else if (fDecayChannel == 22) {
    // for the Lc->V0+bachelor, setting the second and third elements of the array of the track cuts to those for the V0 daughters
    trackCuts[0]=fCuts->GetTrackCuts();
    trackCuts[1]=fCuts->GetTrackCutsV0daughters();
    trackCuts[2]=fCuts->GetTrackCutsV0daughters();
  }
  else {
    for (Int_t iProng = 0; iProng<cfVtxHF->GetNProngs(); iProng++){
      trackCuts[iProng]=fCuts->GetTrackCuts();
    }
  }

  //General settings: vertex, feed down and fill reco container with generated values.
  cfVtxHF->SetRecoPrimVertex(zPrimVertex);
  cfVtxHF->SetMCPrimaryVertex(zMCVertex);
  cfVtxHF->SetFillFromGenerated(fFillFromGenerated);
  cfVtxHF->SetNVar(fNvar);
  cfVtxHF->SetFakeSelection(fFakeSelection);
  cfVtxHF->SetRejectCandidateIfNotFromQuark(fRejectIfNoQuark);
  cfVtxHF->SetConfiguration(fConfiguration);

  // switch-off the trigger class selection (doesn't work for MC)
  fCuts->SetTriggerClass("");

  // MC vertex, to be used, in case, for pp
  if (fUseMCVertex) fCuts->SetUseMCVertex();

  if (fCentralitySelection){ // keep only the requested centrality

    if(fCuts->IsEventSelectedInCentrality(aodEvent)!=0) {
      delete[] containerInput;
      delete[] containerInputMC;
      delete [] trackCuts;
      delete cfVtxHF;
      return;
    }
  }  else { // keep all centralities
    fCuts->SetMinCentrality(0.);
    fCuts->SetMaxCentrality(100.);
  }

  Float_t centValue = 0.;
  if(!fIsPPData) centValue = fCuts->GetCentrality(aodEvent);
  cfVtxHF->SetCentralityValue(centValue);

  // multiplicity estimator with VZERO
  Double_t vzeroMult=0;
  AliAODVZERO *vzeroAOD = (AliAODVZERO*)aodEvent->GetVZEROData();
  if(vzeroAOD) vzeroMult = vzeroAOD->GetMTotV0A() +  vzeroAOD->GetMTotV0C();

  Double_t multiplicity = nTracklets; // set to the Ntracklet estimator
  if(fMultiplicityEstimator==kVZERO) { multiplicity = vzeroMult; }

  cfVtxHF->SetMultiplicity(multiplicity);

  //  printf("Multiplicity estimator %d, value %2.2f\n",fMultiplicityEstimator,multiplicity);

  for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) {
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!mcPart){
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }

    //counting c quarks
    cquarks += cfVtxHF->MCcquarkCounting(mcPart);

    // check the MC-level cuts, must be the desidered particle
    if (!fCFManager->CheckParticleCuts(0, mcPart)) {
      AliDebug(2,"Check the MC-level cuts - not desidered particle");
      continue;  // 0 stands for MC level
    }
    else {
      AliDebug(3, Form("\n\n---> COOL! we found a particle (particle %d)!!! with PDG code = %d \n\n", iPart, mcPart->GetPdgCode()));
    }
    cfVtxHF->SetMCCandidateParam(iPart);


    if (!(cfVtxHF->SetLabelArray())){
      AliDebug(2,Form("Impossible to set the label array for particle %d (decaychannel = %d)", iPart, fDecayChannel));
      continue;
    }

    //check the candiate family at MC level
    if (!(cfVtxHF->CheckMCPartFamily(mcPart, mcArray))) {
      AliDebug(2,Form("Check on the family wrong for particle %d!!! (decaychannel = %d)", iPart, fDecayChannel));
      continue;
    }
    else{
      AliDebug(2,Form("Check on the family OK for particle %d!!! (decaychannel = %d)", iPart, fDecayChannel));
    }

    //Fill the MC container
    Bool_t mcContainerFilled = cfVtxHF -> FillMCContainer(containerInputMC);
    AliDebug(2, Form("particle = %d mcContainerFilled = %d", iPart, mcContainerFilled));
    if (mcContainerFilled) {
      if (fUseWeight){
        if (fHistoPtWeight) { // using an histogram as weight function
          AliDebug(2,"Using Histogram as Pt weight function");
          fWeight = eventWeight*GetPtWeightFromHistogram(containerInputMC[0]);
        }
        else if (fFuncWeight){ // using user-defined function
          AliDebug(2,"Using function");
          fWeight = eventWeight*fFuncWeight->Eval(containerInputMC[0]);
        }
        else{ // using FONLL
          AliDebug(2,"Using FONLL");
          fWeight = eventWeight*GetWeight(containerInputMC[0]);
        }
        AliDebug(2,Form("pt = %f, weight = %f",containerInputMC[0], fWeight));
      }
      if (!fCuts->IsInFiducialAcceptance(containerInputMC[0],containerInputMC[1])) {
        AliDebug(3, Form("Not in limited acceptance, containerInputMC[0] = %f, containerInputMC[1] = %f", containerInputMC[0], containerInputMC[1]));
        continue;
      }
      else{
        AliDebug(3, Form("YES!! in limited acceptance, containerInputMC[0] = %f, containerInputMC[1] = %f", containerInputMC[0],containerInputMC[1]));
      }

      //MC Limited Acceptance
      if (TMath::Abs(containerInputMC[1]) < 0.5) {
        fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGeneratedLimAcc, fWeight);
        icountGenLimAcc++;
	AliDebug(3,"MC Lim Acc container filled\n");
	if (fCFManager->GetParticleContainer()->GetNStep() == kStepGenLimAccNoAcc+1) {
	  if (!(cfVtxHF-> MCAcceptanceStep())) {
	    fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGenLimAccNoAcc, fWeight);
	    icountGenLimAccNoAcc++;
	    AliDebug(3,"MC Lim Acc No Acc container filled\n");
	  }
	}
      }

      //MC
      fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepGenerated, fWeight);
      icountMC++;
      AliDebug(3,"MC container filled \n");

      // MC in acceptance
      // check the MC-Acceptance level cuts
      // since standard CF functions are not applicable, using Kine Cuts on daughters
      Bool_t mcAccepStep = cfVtxHF-> MCAcceptanceStep();
      if (mcAccepStep){
        fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepAcceptance, fWeight);
        AliDebug(3,"MC acceptance cut passed\n");
        icountAcc++;

        //MC Vertex step
        if (fCuts->IsEventSelected(aodEvent)){
          // filling the container if the vertex is ok
          fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepVertex, fWeight) ;
          AliDebug(3,"Vertex cut passed and container filled\n");
          icountVertex++;

          //mc Refit requirement
          Bool_t mcRefitStep = cfVtxHF->MCRefitStep(aodEvent, &trackCuts[0]);
          if (mcRefitStep){
            fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepRefit, fWeight);
            AliDebug(3,"MC Refit cut passed and container filled\n");
            icountRefit++;
          }
          else{
            AliDebug(3,"MC Refit cut not passed\n");
            continue;
          }
        }
        else{
          AliDebug (3, "MC vertex step not passed\n");
          continue;
        }
      }
      else{
        AliDebug (3, "MC in acceptance step not passed\n");
        continue;
      }
    }
    else {
      AliDebug (3, "MC container not filled\n");
    }
  }

  if (cquarks<2) AliDebug(2,Form("Event with %d c-quarks", cquarks));
  AliDebug(2,Form("Found %i MC particles that are %s!!",icountMC,fPartName.Data()));
  AliDebug(2,Form("Found %i MC particles that are %s and satisfy Acc cuts!!",icountAcc,fPartName.Data()));
  AliDebug(2,Form("Found %i MC particles that are %s and satisfy Vertex cuts!!",icountVertex,fPartName.Data()));
  AliDebug(2,Form("Found %i MC particles that are %s and satisfy Refit cuts!!",icountRefit,fPartName.Data()));

  // Now go to rec level
  fCountMC += icountMC;
  fCountGenLimAcc += icountGenLimAcc;
  fCountGenLimAccNoAcc += icountGenLimAccNoAcc;
  fCountAcc += icountAcc;
  fCountVertex+= icountVertex;
  fCountRefit+= icountRefit;

  AliDebug(2,Form("Found %d vertices for decay channel %d",arrayBranch->GetEntriesFast(),fDecayChannel));
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();

  for(Int_t iCandid = 0; iCandid<arrayBranch->GetEntriesFast();iCandid++){
    AliAODRecoDecayHF* charmCandidate=0x0;
    switch (fDecayChannel){
    case 2:{
      charmCandidate = (AliAODRecoDecayHF2Prong*)arrayBranch->At(iCandid);
      if(!vHF->FillRecoCand(aodEvent,(AliAODRecoDecayHF2Prong*)charmCandidate)) continue;
			break;
    }
    case 21:{
      charmCandidate = (AliAODRecoCascadeHF*)arrayBranch->At(iCandid);
      if(!vHF->FillRecoCasc(aodEvent,((AliAODRecoCascadeHF*)charmCandidate),kTRUE)) continue; //DStar
			break;
    }
    case 22:{
      charmCandidate = (AliAODRecoCascadeHF*)arrayBranch->At(iCandid);
      if(!vHF->FillRecoCasc(aodEvent,((AliAODRecoCascadeHF*)charmCandidate),kFALSE)) continue; //Cascade
			break;
    }
    case 31:
    case 32:
    case 33:{
      charmCandidate = (AliAODRecoDecayHF3Prong*)arrayBranch->At(iCandid);
      if(!vHF->FillRecoCand(aodEvent,(AliAODRecoDecayHF3Prong*)charmCandidate)) continue;
			break;
    }
    case 4:{
      charmCandidate = (AliAODRecoDecayHF4Prong*)arrayBranch->At(iCandid);
      break;
    }
    default:
      break;
    }

    Bool_t unsetvtx=kFALSE;
    if(!charmCandidate->GetOwnPrimaryVtx()) {
      charmCandidate->SetOwnPrimaryVtx(aodVtx); // needed to compute all variables
      unsetvtx=kTRUE;
    }

    Bool_t signAssociation = cfVtxHF->SetRecoCandidateParam((AliAODRecoDecayHF*)charmCandidate);
    if (!signAssociation){
      if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
      continue;
    }

    Int_t isPartOrAntipart = cfVtxHF->CheckReflexion(fSign);
    if (isPartOrAntipart == 0){
      AliDebug(2, Form("The candidate pdg code doesn't match the requirement set in the task (fSign = %d)",fSign));
      if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
      continue;
    }

    AliDebug(3,Form("iCandid=%d - signAssociation=%d, isPartOrAntipart=%d",iCandid, signAssociation, isPartOrAntipart));

    Bool_t recoContFilled = cfVtxHF->FillRecoContainer(containerInput);
    AliDebug(3, Form("CF task: RecoContFilled for candidate %d is %d", iCandid, (Int_t)recoContFilled));
    if (recoContFilled){

      // weight according to pt
      if (fUseWeight){
        if (fHistoPtWeight) {
          AliDebug(2,"Using Histogram as Pt weight function");
          fWeight = eventWeight*GetPtWeightFromHistogram(containerInput[0]);
        }
        else if (fFuncWeight){ // using user-defined function
          AliDebug(2, "Using function");
          fWeight = eventWeight*fFuncWeight->Eval(containerInput[0]);
        }
        else{ // using FONLL
          AliDebug(2, "Using FONLL");
          fWeight = eventWeight*GetWeight(containerInput[0]);
        }
        AliDebug(2, Form("pt = %f, weight = %f",containerInput[0], fWeight));
      }

      if (!fCuts->IsInFiducialAcceptance(containerInput[0],containerInput[1])){
        if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
        continue;
      }
      //Reco Step
      Bool_t recoStep = cfVtxHF->RecoStep();
      if (recoStep) AliDebug(2, Form("particle = %d --> CF task: Reco step for candidate %d is %d", iCandid, iCandid, (Int_t)recoStep));
      Bool_t vtxCheck = fCuts->IsEventSelected(aodEvent);


      // Selection on the filtering bit
      Bool_t isBitSelected = kTRUE;
      if(fDecayChannel==2) {
        if(fUseSelectionBit && !charmCandidate->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)) isBitSelected = kFALSE;
      }else if(fDecayChannel==31){
        if(fUseSelectionBit && !charmCandidate->HasSelectionBit(AliRDHFCuts::kDplusCuts)) isBitSelected = kFALSE;
      }else if(fDecayChannel==33){
        if(fUseSelectionBit && !charmCandidate->HasSelectionBit(AliRDHFCuts::kDsCuts)) isBitSelected = kFALSE;
      }else if(fDecayChannel==32){
        if(fUseSelectionBit && !charmCandidate->HasSelectionBit(AliRDHFCuts::kLcCuts)) isBitSelected = kFALSE;
      }else if(fDecayChannel==22){
        if(!((dynamic_cast<AliAODRecoCascadeHF*>(charmCandidate))->CheckCascadeFlags())) isBitSelected = kFALSE; // select only Lc among cascade candidates
      }
      if(!isBitSelected){
        if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
        continue;
      }


      if (recoStep && recoContFilled && vtxCheck){
        fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed, fWeight) ;
        icountReco++;
        AliDebug(3,"Reco step  passed and container filled\n");

        //Reco in the acceptance -- take care of UNFOLDING!!!!
        Bool_t recoAcceptanceStep = cfVtxHF->RecoAcceptStep(&trackCuts[0]);
        if (recoAcceptanceStep) {
          fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoAcceptance, fWeight) ;
          icountRecoAcc++;
          AliDebug(3,"Reco acceptance cut passed and container filled\n");

          if(fAcceptanceUnf){
            Double_t fill[4]; //fill response matrix
            Bool_t bUnfolding = cfVtxHF -> FillUnfoldingMatrix(fPDGcode,fill);
            if (bUnfolding) fCorrelation->Fill(fill);
          }

          //Number of ITS cluster requirements
          Int_t recoITSnCluster = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kTracks, aodEvent);
          if (recoITSnCluster){
            fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoITSClusters, fWeight) ;
            icountRecoITSClusters++;
            AliDebug(3,"Reco n ITS cluster cut passed and container filled\n");

            Bool_t iscutsusingpid = fCuts->GetIsUsePID();
            Int_t recoAnalysisCuts = -1, recoPidSelection = -1;
            fCuts->SetUsePID(kFALSE);
            recoAnalysisCuts = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kCandidate, aodEvent);

            if (fDecayChannel==33){ // Ds case, where more possibilities are considered
              Bool_t keepDs=ProcessDs(recoAnalysisCuts);
              if(keepDs) recoAnalysisCuts=3;
            }
            else if (fDecayChannel==22){ // Lc->V0+bachelor case, where more possibilities are considered
              Bool_t keepLctoV0bachelor = ProcessLctoV0Bachelor(recoAnalysisCuts);
              if (keepLctoV0bachelor) recoAnalysisCuts=3;
              AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
              AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
              AliPIDResponse* pidResponse = inputHandler->GetPIDResponse();
              if (fUseAdditionalCuts){
                if (!((AliCFVertexingHFCascade*)cfVtxHF)->CheckAdditionalCuts(pidResponse)) recoAnalysisCuts = -1;
              }
            }

            fCuts->SetUsePID(iscutsusingpid); //restoring usage of the PID from the cuts object
            Bool_t tempAn=(recoAnalysisCuts == 3 || recoAnalysisCuts == isPartOrAntipart);
            if (fDecayChannel == 22) tempAn = (recoAnalysisCuts == 3);
            if (fDecayChannel == 32) tempAn=(recoAnalysisCuts >0 || recoAnalysisCuts == isPartOrAntipart);

            if (tempAn){
              fCFManager->GetParticleContainer()->Fill(containerInput, kStepRecoPPR, fWeight);
              icountRecoPPR++;
              AliDebug(3,"Reco Analysis cuts passed and container filled \n");
              //pid selection
              //recoPidSelection = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kPID);
              //if((fCuts->CombineSelectionLevels(3,recoAnalysisCuts,recoPidSelection)==isPartOrAntipart)||(fCuts->CombineSelectionLevels(3,recoAnalysisCuts,recoPidSelection)==3)){
              recoPidSelection = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kCandidate, aodEvent);

              if (fDecayChannel==33){ // Ds case, where more possibilities are considered
                Bool_t keepDs=ProcessDs(recoPidSelection);
                if(keepDs) recoPidSelection=3;
              } else if (fDecayChannel==22){ // Lc->V0+bachelor case, where more possibilities are considered
                Bool_t keepLctoV0bachelor=ProcessLctoV0Bachelor(recoPidSelection);
                if (keepLctoV0bachelor) recoPidSelection=3;
              }

              Bool_t tempPid=(recoPidSelection == 3 || recoPidSelection == isPartOrAntipart);
              if (fDecayChannel == 22) tempPid = (recoPidSelection == 3);
              if (fDecayChannel == 32) tempPid=(recoPidSelection >0 || recoPidSelection == isPartOrAntipart);

              if (tempPid){
                Double_t weigPID = 1.;
                if (fDecayChannel == 2){ // D0 with Bayesian PID using weights
                  if(((AliRDHFCutsD0toKpi*)fCuts)->GetCombPID() && (((AliRDHFCutsD0toKpi*)fCuts)->GetBayesianStrategy() == AliRDHFCutsD0toKpi::kBayesWeight || ((AliRDHFCutsD0toKpi*)fCuts)->GetBayesianStrategy() == AliRDHFCutsD0toKpi::kBayesWeightNoFilter)){
                    if (isPartOrAntipart == 1){
                      weigPID = ((AliRDHFCutsD0toKpi*)fCuts)->GetWeightsNegative()[AliPID::kKaon] * ((AliRDHFCutsD0toKpi*)fCuts)->GetWeightsPositive()[AliPID::kPion];
                    }else if (isPartOrAntipart == 2){
                      weigPID = ((AliRDHFCutsD0toKpi*)fCuts)->GetWeightsPositive()[AliPID::kKaon] * ((AliRDHFCutsD0toKpi*)fCuts)->GetWeightsNegative()[AliPID::kPion];
                    }
                    if ((weigPID  < 0) || (weigPID > 1)) weigPID = 0.;
                  }
                }else if (fDecayChannel == 33){ // Ds with Bayesian PID using weights
                  if(((AliRDHFCutsDstoKKpi*)fCuts)->GetPidOption()==AliRDHFCutsDstoKKpi::kBayesianWeights){
                    Int_t labDau0=((AliAODTrack*)charmCandidate->GetDaughter(0))->GetLabel();
                    AliAODMCParticle* firstDau=(AliAODMCParticle*)mcArray->UncheckedAt(TMath::Abs(labDau0));
                    if(firstDau){
                      Int_t pdgCode0=TMath::Abs(firstDau->GetPdgCode());
                      if(pdgCode0==321){
                        weigPID=((AliRDHFCutsDstoKKpi*)fCuts)->GetWeightForKKpi();
                      }else if(pdgCode0==211){
                        weigPID=((AliRDHFCutsDstoKKpi*)fCuts)->GetWeightForpiKK();
                      }
                      if ((weigPID  < 0) || (weigPID > 1)) weigPID = 0.;
                    }else{
                      weigPID=0.;
                    }
                  }
                }
                fCFManager->GetParticleContainer()->Fill(containerInput, kStepRecoPID, fWeight*weigPID);
                icountRecoPID++;
                AliDebug(3,"Reco PID cuts passed and container filled \n");
                if(!fAcceptanceUnf){
                  Double_t fill[4]; //fill response matrix
                  Bool_t bUnfolding = cfVtxHF -> FillUnfoldingMatrix(fPDGcode,fill);
                  if (bUnfolding) fCorrelation->Fill(fill);
                }
              }
              else {
                AliDebug(3, "Analysis Cuts step not passed \n");
                if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
                continue;
              }
            }
            else {
              AliDebug(3, "PID selection not passed \n");
              if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
              continue;
            }
          }
          else{
            AliDebug(3, "Number of ITS cluster step not passed\n");
            if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
            continue;
          }
        }
        else{
          AliDebug(3, "Reco acceptance step not passed\n");
          if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
          continue;
        }
      }
      else {
        AliDebug(3, "Reco step not passed\n");
        if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
        continue;
      }
    }

    if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
  } // end loop on candidate
  delete vHF;
  fCountReco+= icountReco;
  fCountRecoAcc+= icountRecoAcc;
  fCountRecoITSClusters+= icountRecoITSClusters;
  fCountRecoPPR+= icountRecoPPR;
  fCountRecoPID+= icountRecoPID;

  fHistEventsProcessed->Fill(1.5);

  delete[] containerInput;
  delete[] containerInputMC;
  delete cfVtxHF;
  if (trackCuts){
    //	for (Int_t i=0; i<cfVtxHF->GetNProngs(); i++){
    //		delete [] trackCuts[i];
    //	}
    delete [] trackCuts;
  }


}

//___________________________________________________________________________
void AliCFTaskVertexingHF::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliAnalysisTaskSE::Terminate();

  AliInfo(Form("Found %i MC particles that are %s in MC, in %d events",fCountMC,fPartName.Data(),fEvents));
  AliInfo(Form("Found %i MC particles that are %s in MC in limited acceptance, in %d events",fCountGenLimAcc,fPartName.Data(),fEvents));
  AliInfo(Form("Found %i MC particles that are %s in MC in limited acceptance and don't satisfy Acc cuts, in %d events",fCountGenLimAccNoAcc,fPartName.Data(),fEvents));
  AliInfo(Form("Found %i MC particles that are %s in MC and satisfy Acc cuts, in %d events",fCountAcc,fPartName.Data(),fEvents));
  AliInfo(Form("Found %i MC particles that are %s in MC and satisfy Acc cuts, and satisfy Vertex requirement in %d events",fCountVertex,fPartName.Data(),fEvents));
  AliInfo(Form("Found %i MC particles that are %s in MC and satisfy Acc cuts, and satisfy ITS+TPC refit requirementin %d events",fCountRefit,fPartName.Data(),fEvents));
  AliInfo(Form("Found %i reco %s that are decaying in %s, in %d events",fCountReco,fPartName.Data(),fDauNames.Data(),fEvents));
  AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and are in the requested acceptance, in %d events",fCountRecoAcc,fPartName.Data(),fDauNames.Data(),fEvents));
  AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and have at least 5 clusters in ITS, in %d events",fCountRecoITSClusters,fPartName.Data(),fDauNames.Data(),fEvents));
  AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and satisfy PPR cuts, in %d events",fCountRecoPPR,fPartName.Data(),fDauNames.Data(),fEvents));
  AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and satisfy PPR+PID cuts, in %d events",fCountRecoPID,fPartName.Data(),fDauNames.Data(),fEvents));

  // draw some example plots....
  AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
  if(!cont) {
    printf("CONTAINER NOT FOUND\n");
    return;
  }
  // projecting the containers to obtain histograms
  // first argument = variable, second argument = step

  TH1D** h = new TH1D*[3];
  Int_t nvarToPlot = 0;
  if (fConfiguration == kSnail){
    //h = new TH1D[3][12];
    for (Int_t ih = 0; ih<3; ih++){
      if(fDecayChannel==22){
        nvarToPlot = 16;
      } else {
        nvarToPlot = 12;
      }
      h[ih] = new TH1D[nvarToPlot];
    }
    for(Int_t iC=1;iC<nvarToPlot; iC++){
      // MC-level
      h[0][iC] =   *(cont->ShowProjection(iC,0));
      // MC-Acceptance level
      h[1][iC] =   *(cont->ShowProjection(iC,1));
      // Reco-level
      h[2][iC] =   *(cont->ShowProjection(iC,4));
    }
  }
  else {
    //h = new TH1D[3][12];
    nvarToPlot = 8;
    for (Int_t ih = 0; ih<3; ih++){
      h[ih] = new TH1D[nvarToPlot];
    }
    for(Int_t iC=0;iC<nvarToPlot; iC++){
      // MC-level
      h[0][iC] =   *(cont->ShowProjection(iC,0));
      // MC-Acceptance level
      h[1][iC] =   *(cont->ShowProjection(iC,1));
      // Reco-level
      h[2][iC] =   *(cont->ShowProjection(iC,4));
    }
  }
  TString* titles;
  //Int_t nvarToPlot = 0;
  if (fConfiguration == kSnail){
    if(fDecayChannel==31){
      //nvarToPlot = 12;
      titles = new TString[nvarToPlot];
      titles[0]="pT_Dplus (GeV/c)";
      titles[1]="rapidity";
      titles[2]="phi (rad)";
      titles[3]="cT (#mum)";
      titles[4]="cosPointingAngle";
      titles[5]="pT_1 (GeV/c)";
      titles[6]="pT_2 (GeV/c)";
      titles[7]="pT_3 (GeV/c)";
      titles[8]="d0_1 (#mum)";
      titles[9]="d0_2 (#mum)";
      titles[10]="d0_3 (#mum)";
      titles[11]="zVertex (cm)";
    } else if (fDecayChannel==22) {
      //nvarToPlot = 16;
      titles = new TString[nvarToPlot];
      titles[0]="p_{T}(#Lambda_{c}) [GeV/c]";
      titles[1]="y(#Lambda_{c})";
      titles[2]="#varphi(#Lambda_{c}) [rad]";
      titles[3]="onTheFlyStatusV0";
      titles[4]="z_{vtx} [cm]";
      titles[5]="centrality";
      titles[6]="fake";
      titles[7]="multiplicity";
      //titles[8]="pT(bachelor) [GeV/c]";
      titles[8]="p(bachelor) [GeV/c]";
      titles[9]="p_{T}(V0) [GeV/c]";
      titles[10]="y(V0)";
      titles[11]="#varphi(V0) [rad]";
      titles[12]="m_{inv}(#pi^{+}#pi^{+}) [GeV/c^{2}]";
      titles[13]="dcaV0 (nSigma)";
      titles[14]="cosine pointing angle (V0)";
      titles[15]="cosine pointing angle (#Lambda_{c})";
      //titles[16]="c#tauV0 (#mum)";
      //titles[17]="c#tau (#mum)";
    } else {
      //nvarToPlot = 12;
      titles = new TString[nvarToPlot];
      titles[0]="pT_D0 (GeV/c)";
      titles[1]="rapidity";
      titles[2]="cosThetaStar";
      titles[3]="pT_pi (GeV/c)";
      titles[4]="pT_K (Gev/c)";
      titles[5]="cT (#mum)";
      titles[6]="dca (#mum)";
      titles[7]="d0_pi (#mum)";
      titles[8]="d0_K (#mum)";
      titles[9]="d0xd0 (#mum^2)";
      titles[10]="cosPointingAngle";
      titles[11]="phi (rad)";
    }
  }
  else {
    //nvarToPlot = 8;
    titles = new TString[nvarToPlot];
    if (fDecayChannel==22) {
      titles[0]="p_{T}(#Lambda_{c}) [GeV/c]";
      titles[1]="y(#Lambda_{c})";
      titles[2]="#varphi(#Lambda_{c}) [rad]";
      titles[3]="onTheFlyStatusV0";
      titles[4]="z_{vtx} [cm]";
      titles[5]="centrality";
      titles[6]="fake";
      titles[7]="multiplicity";
    } else {
      titles[0]="pT_candidate (GeV/c)";
      titles[1]="rapidity";
      titles[2]="cT (#mum)";
      titles[3]="phi";
      titles[4]="z_{vtx}";
      titles[5]="centrality";
      titles[6]="fake";
      titles[7]="multiplicity";
    }
  }

  Int_t markers[16]={20,24,21,25,27,28,
                     20,24,21,25,27,28,
                     20,24,21,25};
  Int_t colors[3]={2,8,4};
  for(Int_t iC=0;iC<nvarToPlot; iC++){
    for(Int_t iStep=0;iStep<3;iStep++){
      h[iStep][iC].SetTitle(titles[iC].Data());
      h[iStep][iC].GetXaxis()->SetTitle(titles[iC].Data());
      Double_t maxh=h[iStep][iC].GetMaximum();
      h[iStep][iC].GetYaxis()->SetRangeUser(0,maxh*1.2);
      h[iStep][iC].SetMarkerStyle(markers[iC]);
      h[iStep][iC].SetMarkerColor(colors[iStep]);
    }
  }

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  // drawing in 2 separate canvas for a matter of clearity
  TCanvas * c1 =new TCanvas(Form("c1New_%d",fDecayChannel),"Vars 0, 1, 2, 3",1100,1200);
  c1->Divide(3,4);
  Int_t iPad=1;
  for(Int_t iVar=0; iVar<4; iVar++){
    c1->cd(iPad++);
    h[0][iVar].DrawCopy("p");
    c1->cd(iPad++);
    h[1][iVar].DrawCopy("p");
    c1->cd(iPad++);
    h[2][iVar].DrawCopy("p");
  }

  TCanvas * c2 =new TCanvas(Form("c2New_%d",fDecayChannel),"Vars 4, 5, 6, 7",1100,1200);
  c2->Divide(3,4);
  iPad=1;
  for(Int_t iVar=4; iVar<8; iVar++){
    c2->cd(iPad++);
    h[0][iVar].DrawCopy("p");
    c2->cd(iPad++);
    h[1][iVar].DrawCopy("p");
    c2->cd(iPad++);
    h[2][iVar].DrawCopy("p");
  }

  if (fConfiguration == kSnail){
    TCanvas * c3 =new TCanvas(Form("c3New_%d",fDecayChannel),"Vars 8, 9, 10, 11",1100,1200);
    c3->Divide(3,4);
    iPad=1;
    for(Int_t iVar=8; iVar<12; iVar++){
      c3->cd(iPad++);
      h[0][iVar].DrawCopy("p");
      c3->cd(iPad++);
      h[1][iVar].DrawCopy("p");
      c3->cd(iPad++);
      h[2][iVar].DrawCopy("p");
    }
    if (fDecayChannel==22) {
      TCanvas * c4 =new TCanvas(Form("c4New_%d",fDecayChannel),"Vars 12, 13, 14, 15",1100,1200);
      c4->Divide(3,4);
      iPad=1;
      for(Int_t iVar=12; iVar<16; iVar++){
        c4->cd(iPad++);
        h[0][iVar].DrawCopy("p");
        c4->cd(iPad++);
        h[1][iVar].DrawCopy("p");
        c4->cd(iPad++);
        h[2][iVar].DrawCopy("p");
      }
    }
  }

  /*
    THnSparseD* hcorr = dynamic_cast<THnSparseD*> (GetOutputData(3));

    TH2D* corr1 = hcorr->Projection(0,2);
    TH2D* corr2 = hcorr->Projection(1,3);

    TCanvas * c7 =new TCanvas(Form("c7New_%d",fDecayChannel),"",800,400);
    c7->Divide(2,1);
    c7->cd(1);
    corr1->DrawCopy("text");
    c7->cd(2);
    corr2->DrawCopy("text");
  */
  TFile* file_projection = new TFile("CFtaskHFprojectionNew.root","RECREATE");

  //	corr1->Write();
  //	corr2->Write();

  for(Int_t iC=0;iC<nvarToPlot; iC++){
    for(Int_t iStep=0;iStep<3;iStep++){
      h[iStep][iC].Write(Form("Step%d_%s",iStep,titles[iC].Data()));
    }
  }
  file_projection->Close();
  for (Int_t ih = 0; ih<3; ih++) delete [] h[ih];
  delete [] h;
  delete [] titles;

}

//___________________________________________________________________________
void AliCFTaskVertexingHF::UserCreateOutputObjects()
{
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());

  //slot #1
  OpenFile(1);
  const char* nameoutput=GetOutputSlot(1)->GetContainer()->GetName();
  fHistEventsProcessed = new TH1I(nameoutput,"",2,0,2) ;
  fHistEventsProcessed->GetXaxis()->SetBinLabel(1,"Events processed (all)");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(2,"Events analyzed (after selection)");

  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fCorrelation) ;

}


//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL276overLHC12a17a(){
  // ad-hoc weight function from ratio of
  // pt spectra from FONLL 2.76 TeV and
  // pt spectra from MC production LHC12a17a (PYTHIA Perugia0 with pthard bins)
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","[0]+[1]*TMath::Exp(-[2]*x)",0.,50.);
  fFuncWeight->SetParameter(0,4.63891e-02);
  fFuncWeight->SetParameter(1,1.51674e+01);
  fFuncWeight->SetParameter(2,4.09941e-01);
  fUseWeight=kTRUE;
}
//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromDataPbPb276overLHC12a17a(){
  // ad-hoc weight function from ratio of
  // D0 pt spectra in PbPb 2011 0-10% centrality and
  // pt spectra from MC production LHC12a17a (PYTHIA Perugia0 with pthard bins)
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","[0]+[1]/TMath::Power(x,[2])",0.05,50.);
  fFuncWeight->SetParameter(0,1.43116e-02);
  fFuncWeight->SetParameter(1,4.37758e+02);
  fFuncWeight->SetParameter(2,3.08583);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL276overLHC12a17b(){
  // weight function from the ratio of the LHC12a17b MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,50.);
  fFuncWeight->SetParameters(1.92381e+01, 5.05055e+00, 1.05314e+01, 2.5, 1.88214e-03, 3.44871e+00, -9.74325e-02, 1.97671e+00, -3.21278e-01);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL276andBAMPSoverLHC12a17b(){
  // weight function from the ratio of the LHC12a17b MC
  // and FONLL calculations for pp data
  // corrected by the BAMPS Raa calculation for 30-50% CC
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,50.);
  fFuncWeight->SetParameters(6.10443e+00, 1.53487e+00, 1.99474e+00, 2.5, 5.51172e-03, 5.86590e+00, -5.46963e-01, 9.41201e-02, -1.64323e-01);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5anddataoverLHC16i2a(){
    // weight function from the ratio of the LHC16i2a MC
    // 1.5-14 GeV/c using data and 1-1.5, 14-50 GeV/c using FONLL calculations
    
    if(fHistoPtWeight) delete fHistoPtWeight;
    fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
    Float_t binc[500]={ 1.695705, 1.743693, 1.790289, 1.835410, 1.878981, 1.920938, 1.961223, 1.999787, 2.036589, 2.071597, 2.104784, 2.136132, 2.165629, 2.193270, 2.219057, 2.174545, 2.064698, 1.959489, 1.858770, 1.762396, 1.670224, 1.582115, 1.497931, 1.417541, 1.340814, 1.267622, 1.197842, 1.131352, 1.068033, 1.007770, 0.950450, 0.895963, 0.844202, 0.795062, 0.748441, 0.704241, 0.662363, 0.622715, 0.585204, 0.549742, 0.516242, 0.484620, 0.454795, 0.426686, 0.400217, 0.375314, 0.351903, 0.329915, 0.309281, 0.289936, 0.271816, 0.254860, 0.239007, 0.224201, 0.210386, 0.197508, 0.185516, 0.174360, 0.163992, 0.154366, 0.145438, 0.137166, 0.129508, 0.122426, 0.115882, 0.109840, 0.104266, 0.099128, 0.094395, 0.090036, 0.086023, 0.082331, 0.078933, 0.075805, 0.072925, 0.070271, 0.067823, 0.065562, 0.063471, 0.061532, 0.059730, 0.058051, 0.056481, 0.055007, 0.053619, 0.052306, 0.051059, 0.049867, 0.048725, 0.047624, 0.046558, 0.045522, 0.044511, 0.043521, 0.042548, 0.041590, 0.040643, 0.039706, 0.038778, 0.037857, 0.036944, 0.036039, 0.035141, 0.034251, 0.033370, 0.032500, 0.031641, 0.030796, 0.029966, 0.029153, 0.028359, 0.027587, 0.026837, 0.026113, 0.025416, 0.024748, 0.024111, 0.023507, 0.022937, 0.022402, 0.021904, 0.021443, 0.021020, 0.020634, 0.020286, 0.019974, 0.019698, 0.019455, 0.019244, 0.019062, 0.018905, 0.018770, 0.018652, 0.018545, 0.018444, 0.018342, 0.018231, 0.018102, 0.017947, 0.017755, 0.017536, 0.017327, 0.017120, 0.016915, 0.016713, 0.016514, 0.016317, 0.016122, 0.015929, 0.015739, 0.015551, 0.015366, 0.015182, 0.015001, 0.014822, 0.014645, 0.014470, 0.014297, 0.014127, 0.013958, 0.013791, 0.013627, 0.013464, 0.013303, 0.013145, 0.012988, 0.012833, 0.012679, 0.012528, 0.012378, 0.012231, 0.012085, 0.011940, 0.011798, 0.011657, 0.011518, 0.011380, 0.011244, 0.011110, 0.010978, 0.010846, 0.010717, 0.010589, 0.010463, 0.010338, 0.010214, 0.010092, 0.009972, 0.009853, 0.009735, 0.009619, 0.009504, 0.009391, 0.009279, 0.009168, 0.009058, 0.008950, 0.008843, 0.008738, 0.008633, 0.008530, 0.008429, 0.008328, 0.008229, 0.008130, 0.008033, 0.007937, 0.007843, 0.007749, 0.007656, 0.007565, 0.007475, 0.007385, 0.007297, 0.007210, 0.007124, 0.007039, 0.006955, 0.006872, 0.006790, 0.006709, 0.006629, 0.006550, 0.006471, 0.006394, 0.006318, 0.006242, 0.006168, 0.006094, 0.006022, 0.005950, 0.005879, 0.005808, 0.005739, 0.005671, 0.005603, 0.005536, 0.005470, 0.005405, 0.005340, 0.005276, 0.005213, 0.005151, 0.005090, 0.005029, 0.004969, 0.004909, 0.004851, 0.004793, 0.004736, 0.004679, 0.004623, 0.004568, 0.004514, 0.004460, 0.004406, 0.004354, 0.004302, 0.004251, 0.004200, 0.004150, 0.004100, 0.004051, 0.004003, 0.003955, 0.003908, 0.003861, 0.003815, 0.003770, 0.003725, 0.003680, 0.003636, 0.003593, 0.003550, 0.003507, 0.003466, 0.003424, 0.003383, 0.003343, 0.003303, 0.003264, 0.003225, 0.003186, 0.003148, 0.003110, 0.003073, 0.003037, 0.003000, 0.002965, 0.002929, 0.002894, 0.002860, 0.002826, 0.002792, 0.002758, 0.002726, 0.002693, 0.002661, 0.002629, 0.002598, 0.002567, 0.002536, 0.002506, 0.002476, 0.002446, 0.002417, 0.002388, 0.002360, 0.002332, 0.002304, 0.002276, 0.002249, 0.002222, 0.002196, 0.002169, 0.002144, 0.002118, 0.002093, 0.002068, 0.002043, 0.002019, 0.001995, 0.001971, 0.001947, 0.001924, 0.001901, 0.001878, 0.001856, 0.001834, 0.001812, 0.001790, 0.001769, 0.001748, 0.001727, 0.001706, 0.001686, 0.001666, 0.001646, 0.001626, 0.001607, 0.001588, 0.001569, 0.001550, 0.001531, 0.001513, 0.001495, 0.001477, 0.001460, 0.001442, 0.001425, 0.001408, 0.001391, 0.001374, 0.001358, 0.001342, 0.001326, 0.001310, 0.001294, 0.001279, 0.001264, 0.001249, 0.001234, 0.001219, 0.001204, 0.001190, 0.001176, 0.001162, 0.001148, 0.001134, 0.001121, 0.001107, 0.001094, 0.001081, 0.001068, 0.001055, 0.001043, 0.001030, 0.001018, 0.001006, 0.000994, 0.000982, 0.000970, 0.000959, 0.000947, 0.000936, 0.000925, 0.000914, 0.000903, 0.000892, 0.000881, 0.000871, 0.000860, 0.000850, 0.000840, 0.000830, 0.000820, 0.000810, 0.000801, 0.000791, 0.000782, 0.000772, 0.000763, 0.000754, 0.000745, 0.000736, 0.000727, 0.000719, 0.000710, 0.000702, 0.000693, 0.000685, 0.000677, 0.000669, 0.000661, 0.000653, 0.000645, 0.000637, 0.000630, 0.000622, 0.000615, 0.000607, 0.000600, 0.000593, 0.000586, 0.000579, 0.000572, 0.000565, 0.000558, 0.000552, 0.000545, 0.000539, 0.000532, 0.000526, 0.000520, 0.000513, 0.000507, 0.000501, 0.000495, 0.000489, 0.000483, 0.000478, 0.000472, 0.000466, 0.000461, 0.000455, 0.000450, 0.000444, 0.000439, 0.000434, 0.000429, 0.000424, 0.000419, 0.000414, 0.000409, 0.000404, 0.000399, 0.000394, 0.000389, 0.000385, 0.000380, 0.000376, 0.000371, 0.000367, 0.000362, 0.000358, 0.000354, 0.000350, 0.000345, 0.000341, 0.000337, 0.000333, 0.000329, 0.000325, 0.000321, 0.000318, 0.000314, 0.000310, 0.000306, 0.000303, 0.000299, 0.000295, 0.000292, 0.000288, 0.000285, 0.000282, 0.000278, 0.000275, 0.000272, 0.000268, 0.000265, 0.000262, 0.000259, 0.000256, 0.000253, 0.000250, 0.000247, 0.000244, 0.000241, 0.000238, 0.000235};
    for(Int_t i=0; i<500; i++){
        fHistoPtWeight->SetBinContent(i+1,binc[i]);
    }
    //SetWeightHistogram();
    fUseWeight=kTRUE;
}
//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5overLHC16i2abc(){
  // weight function from the ratio of the LHC16i2a+b+c MC
  // and FONLL calculations for pp data

  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",400,0.,40.);
  Float_t binc[400]={1.118416, 1.003458, 0.935514, 0.907222, 0.904359, 0.913668, 0.933906, 0.963898, 0.996388, 1.031708, 1.066404, 1.099683, 1.125805, 1.145181, 1.165910, 1.181905, 1.193425, 1.203891, 1.204726, 1.209411, 1.209943, 1.204763, 1.205291, 1.198912, 1.197390, 1.182005, 1.184194, 1.175994, 1.167881, 1.158348, 1.147190, 1.139833, 1.126940, 1.123322, 1.108389, 1.102199, 1.089464, 1.075874, 1.061964, 1.051429, 1.038113, 1.026668, 1.011441, 0.998567, 0.987658, 0.972434, 0.950068, 0.940758, 0.916880, 0.911931, 0.894512, 0.878691, 0.860589, 0.848025, 0.830774, 0.819399, 0.801134, 0.775276, 0.766382, 0.750495, 0.736935, 0.717529, 0.702637, 0.689152, 0.671334, 0.652030, 0.635696, 0.621365, 0.608362, 0.599019, 0.576024, 0.562136, 0.550938, 0.533587, 0.516410, 0.509744, 0.501655, 0.487402, 0.476469, 0.463762, 0.445979, 0.438088, 0.422214, 0.417467, 0.404357, 0.391450, 0.379996, 0.371201, 0.361497, 0.352912, 0.343189, 0.329183, 0.327662, 0.310783, 0.304525, 0.301007, 0.293306, 0.278332, 0.274419, 0.267361, 0.261459, 0.255514, 0.249293, 0.241129, 0.237600, 0.231343, 0.221982, 0.216872, 0.211094, 0.206954, 0.202333, 0.196572, 0.193274, 0.188240, 0.181817, 0.178364, 0.173614, 0.167135, 0.166055, 0.163423, 0.156557, 0.155821, 0.151985, 0.144909, 0.145062, 0.139720, 0.138873, 0.131892, 0.129969, 0.126509, 0.126978, 0.120451, 0.117661, 0.116300, 0.115604, 0.112215, 0.109237, 0.107720, 0.106419, 0.102050, 0.102777, 0.097406, 0.098447, 0.095964, 0.093868, 0.092430, 0.089329, 0.088249, 0.085881, 0.084417, 0.085498, 0.082444, 0.079151, 0.079565, 0.077811, 0.077293, 0.075218, 0.072445, 0.073054, 0.071545, 0.070279, 0.068046, 0.067854, 0.068092, 0.065378, 0.064405, 0.062060, 0.063391, 0.061718, 0.059616, 0.058913, 0.058895, 0.058311, 0.056320, 0.056527, 0.055349, 0.053701, 0.054735, 0.052264, 0.051277, 0.051554, 0.050545, 0.048995, 0.049507, 0.048466, 0.048156, 0.046809, 0.047600, 0.046078, 0.044801, 0.044113, 0.043700, 0.043530, 0.043396, 0.042556, 0.041048, 0.041657, 0.040394, 0.041314, 0.040720, 0.039656, 0.038478, 0.039276, 0.038777, 0.037730, 0.036918, 0.036466, 0.035827, 0.035285, 0.035963, 0.034371, 0.034757, 0.033205, 0.033666, 0.033266, 0.032583, 0.033570, 0.032102, 0.032107, 0.031464, 0.032160, 0.030091, 0.030564, 0.029464, 0.029613, 0.029626, 0.029512, 0.029324, 0.028607, 0.027628, 0.027251, 0.027072, 0.027077, 0.026724, 0.026961, 0.026303, 0.026237, 0.025454, 0.025133, 0.025365, 0.026014, 0.024807, 0.023901, 0.023459, 0.023405, 0.023654, 0.023981, 0.023675, 0.022493, 0.022781, 0.021801, 0.021704, 0.022372, 0.021189, 0.020681, 0.020779, 0.021324, 0.020558, 0.020901, 0.020586, 0.020808, 0.019276, 0.019516, 0.019706, 0.018935, 0.018632, 0.018516, 0.019187, 0.018916, 0.018039, 0.018208, 0.018045, 0.017628, 0.017916, 0.017711, 0.017838, 0.017222, 0.016565, 0.015733, 0.016264, 0.015826, 0.016090, 0.016622, 0.015802, 0.016621, 0.015441, 0.015309, 0.014860, 0.014935, 0.014968, 0.014443, 0.014485, 0.015136, 0.014078, 0.014414, 0.013908, 0.014071, 0.014078, 0.013766, 0.013436, 0.013507, 0.013480, 0.013224, 0.013041, 0.013935, 0.012885, 0.012453, 0.012528, 0.012492, 0.012225, 0.012542, 0.012706, 0.012136, 0.011902, 0.011560, 0.011448, 0.011861, 0.011271, 0.011831, 0.011159, 0.011171, 0.010966, 0.011311, 0.011002, 0.011130, 0.010995, 0.010450, 0.010663, 0.010678, 0.010492, 0.009861, 0.010507, 0.009916, 0.010121, 0.010029, 0.010046, 0.009370, 0.009647, 0.010104, 0.009282, 0.009830, 0.009403, 0.009148, 0.009172, 0.008893, 0.009158, 0.009019, 0.008780, 0.008579, 0.009063, 0.008634, 0.008988, 0.008265, 0.008581, 0.008575, 0.008690, 0.008181, 0.008352, 0.008150, 0.008430, 0.008256, 0.008119, 0.008453, 0.008447, 0.008021, 0.007938, 0.008025, 0.007718, 0.008127, 0.007651, 0.007590, 0.007316, 0.007839, 0.007504, 0.007341, 0.007527, 0.007263, 0.007668, 0.007306, 0.007271, 0.006910, 0.007257, 0.007260, 0.006810, 0.006967, 0.006887, 0.006867, 0.007202, 0.006829, 0.006370, 0.006710, 0.006417, 0.006361, 0.006800, 0.006410, 0.006323, 0.006790, 0.006322, 0.006673, 0.006547};
  for(Int_t i=0; i<400; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  //SetWeightHistogram();
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5andBAMPSoverLHC16i2abc(){
  // weight function from the ratio of the LHC16i2a+b+c MC
  // and FONLL calculations for pp data
  // corrected by the BAMPS Raa calculation for 30-50% CC
  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",400,0.,40.);
  Float_t binc[400]={2.166180, 1.866117, 1.667595, 1.547176, 1.486661, 1.457891, 1.426949, 1.399055, 1.383278, 1.349383, 1.317009, 1.282321, 1.234257, 1.181136, 1.136655, 1.087523, 1.037912, 0.993256, 0.944746, 0.900948, 0.865869, 0.827193, 0.794424, 0.757723, 0.733020, 0.700164, 0.682189, 0.659872, 0.637918, 0.615749, 0.593020, 0.574402, 0.556158, 0.542663, 0.525494, 0.516038, 0.503629, 0.490980, 0.479143, 0.469005, 0.457749, 0.447668, 0.436803, 0.427073, 0.418282, 0.407867, 0.395093, 0.387861, 0.374742, 0.369462, 0.360146, 0.351991, 0.342990, 0.336259, 0.327730, 0.322382, 0.314602, 0.303874, 0.299820, 0.293049, 0.287539, 0.280329, 0.274866, 0.269939, 0.263299, 0.256057, 0.249215, 0.242170, 0.235704, 0.230709, 0.220529, 0.213921, 0.208394, 0.202424, 0.196700, 0.194943, 0.192620, 0.187894, 0.184411, 0.180204, 0.172915, 0.169077, 0.162201, 0.159636, 0.153904, 0.148296, 0.143282, 0.139306, 0.135561, 0.132342, 0.128696, 0.123444, 0.122873, 0.116544, 0.114197, 0.112878, 0.110018, 0.104547, 0.103222, 0.100707, 0.098622, 0.096513, 0.094295, 0.091334, 0.090122, 0.087870, 0.084894, 0.083729, 0.082265, 0.081404, 0.080323, 0.078750, 0.078132, 0.076781, 0.074823, 0.074050, 0.072614, 0.070093, 0.069828, 0.068907, 0.066189, 0.066054, 0.064600, 0.061757, 0.061986, 0.059862, 0.059656, 0.056807, 0.055956, 0.054386, 0.054507, 0.051629, 0.050358, 0.049702, 0.049331, 0.047814, 0.046476, 0.045762, 0.045142, 0.043224, 0.043484, 0.041282, 0.041794, 0.040809, 0.039985, 0.039439, 0.038181, 0.037782, 0.036831, 0.036264, 0.036790, 0.035535, 0.034173, 0.034409, 0.033659, 0.033308, 0.032290, 0.030981, 0.031121, 0.030361, 0.029708, 0.028653, 0.028461, 0.028449, 0.027208, 0.026697, 0.025623, 0.026069, 0.025279, 0.024332, 0.024341, 0.024629, 0.024677, 0.024117, 0.024490, 0.024257, 0.023804, 0.024537, 0.023692, 0.023502, 0.023888, 0.023673, 0.023193, 0.023684, 0.023429, 0.023521, 0.023014, 0.023346, 0.022544, 0.021866, 0.021477, 0.021224, 0.021089, 0.020972, 0.020515, 0.019739, 0.019982, 0.019328, 0.019719, 0.019387, 0.018833, 0.018227, 0.018558, 0.018276, 0.017738, 0.017460, 0.017365, 0.017178, 0.017033, 0.017478, 0.016817, 0.017119, 0.016463, 0.016802, 0.016711, 0.016475, 0.017083, 0.016441, 0.016548, 0.016320, 0.016786, 0.015804, 0.016153, 0.015668, 0.015843, 0.015810, 0.015651, 0.015454, 0.014981, 0.014376, 0.014089, 0.013906, 0.013818, 0.013549, 0.013580, 0.013160, 0.013040, 0.012566, 0.012324, 0.012353, 0.012582, 0.011915, 0.011401, 0.011112, 0.011008, 0.011046, 0.011119, 0.010954, 0.010439, 0.010604, 0.010179, 0.010163, 0.010507, 0.009981, 0.009771, 0.009846, 0.010134, 0.009798, 0.009991, 0.009869, 0.010005, 0.009295, 0.009438, 0.009557, 0.009210, 0.009088, 0.009057, 0.009412, 0.009306, 0.008899, 0.009009, 0.008952, 0.008764, 0.008926, 0.008842, 0.008924, 0.008634, 0.008322, 0.007920, 0.008205, 0.008000, 0.008151, 0.008438, 0.008037, 0.008472, 0.007886, 0.007835, 0.007621, 0.007675, 0.007707, 0.007452, 0.007489, 0.007841, 0.007308, 0.007497, 0.007248, 0.007348, 0.007367, 0.007227, 0.007097, 0.007179, 0.007209, 0.007115, 0.007059, 0.007588, 0.007058, 0.006862, 0.006945, 0.006965, 0.006856, 0.007075, 0.007209, 0.006925, 0.006830, 0.006672, 0.006645, 0.006923, 0.006615, 0.006982, 0.006622, 0.006666, 0.006579, 0.006823, 0.006673, 0.006786, 0.006740, 0.006440, 0.006606, 0.006650, 0.006568, 0.006206, 0.006646, 0.006305, 0.006468, 0.006442, 0.006486, 0.006080, 0.006291, 0.006622, 0.006113, 0.006506, 0.006254, 0.006114, 0.006161, 0.006002, 0.006211, 0.006146, 0.006012, 0.005902, 0.006264, 0.005996, 0.006271, 0.005793, 0.006043, 0.006067, 0.006177, 0.005842, 0.005991, 0.005872, 0.006102, 0.006003, 0.005930, 0.006201, 0.006224, 0.005937, 0.005901, 0.005992, 0.005788, 0.006121, 0.005787, 0.005766, 0.005582, 0.006006, 0.005774, 0.005672, 0.005841, 0.005660, 0.006000, 0.005741, 0.005737, 0.005475, 0.005773, 0.005799, 0.005462, 0.005610, 0.005569, 0.005574, 0.005871, 0.005589, 0.005234, 0.005535, 0.005314, 0.005288, 0.005676, 0.005371, 0.005319, 0.005734, 0.005360, 0.005679, 0.005593};
  for(Int_t i=0; i<400; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5andTAMUoverLHC16i2abc(){
  // weight function from the ratio of the LHC16i2a+b+c MC
  // and FONLL calculations for pp data
  // corrected by the TAMU Raa calculation for 0-10% CC (not available in 30-50% CC)
  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",400,0.,40.);
  Float_t binc[400]={1.179906, 1.091249, 1.047774, 1.045579, 1.071679, 1.112413, 1.167414, 1.236240, 1.310301, 1.390289, 1.471711, 1.553389, 1.626886, 1.692115, 1.760647, 1.813658, 1.850817, 1.886699, 1.907671, 1.934832, 1.955433, 1.966727, 1.987262, 1.996316, 2.013326, 1.973926, 1.931144, 1.871654, 1.812942, 1.752718, 1.690846, 1.635303, 1.572611, 1.523510, 1.459790, 1.402510, 1.331908, 1.261575, 1.192241, 1.127915, 1.061798, 0.998830, 0.933514, 0.871774, 0.812936, 0.762844, 0.719340, 0.686587, 0.644108, 0.615714, 0.579512, 0.545254, 0.510508, 0.479884, 0.447423, 0.426154, 0.408934, 0.388264, 0.376424, 0.361389, 0.347757, 0.331685, 0.318029, 0.305285, 0.290922, 0.278523, 0.269807, 0.262025, 0.254878, 0.249325, 0.238179, 0.230899, 0.224792, 0.216253, 0.207879, 0.204465, 0.201153, 0.195373, 0.190926, 0.185773, 0.178589, 0.175371, 0.168959, 0.167004, 0.161705, 0.156809, 0.152788, 0.149806, 0.146429, 0.143478, 0.140037, 0.134813, 0.134679, 0.128205, 0.126078, 0.125038, 0.122214, 0.116329, 0.115044, 0.112427, 0.110279, 0.108098, 0.105784, 0.102628, 0.101429, 0.099101, 0.095464, 0.093631, 0.091491, 0.090045, 0.088374, 0.086188, 0.085067, 0.083168, 0.080636, 0.079414, 0.077610, 0.075013, 0.074825, 0.073932, 0.071106, 0.071050, 0.069574, 0.066593, 0.066924, 0.064876, 0.065064, 0.062345, 0.061980, 0.060859, 0.061616, 0.058952, 0.058079, 0.057894, 0.058031, 0.056604, 0.055180, 0.054490, 0.053909, 0.051768, 0.052210, 0.049552, 0.050152, 0.048955, 0.047953, 0.047224, 0.045588, 0.044985, 0.043728, 0.042934, 0.043434, 0.041834, 0.040118, 0.040281, 0.039348, 0.038987, 0.037793, 0.036258, 0.036420, 0.035528, 0.034761, 0.033524, 0.033296, 0.033280, 0.031825, 0.031351, 0.030329, 0.031103, 0.030401, 0.029481, 0.029247, 0.029352, 0.029174, 0.028286, 0.028500, 0.028017, 0.027293, 0.027932, 0.026779, 0.026379, 0.026628, 0.026211, 0.025508, 0.025877, 0.025433, 0.025328, 0.024636, 0.025069, 0.024282, 0.023625, 0.023278, 0.023074, 0.023000, 0.022943, 0.022514, 0.021767, 0.022180, 0.021594, 0.022175, 0.021944, 0.021456, 0.020901, 0.021419, 0.021230, 0.020738, 0.020322, 0.020055, 0.019686, 0.019371, 0.019725, 0.018835, 0.019029, 0.018163, 0.018398, 0.018163, 0.017719, 0.018126, 0.017208, 0.017086, 0.016622, 0.016865, 0.015663, 0.015791, 0.015108, 0.015069, 0.015033, 0.015006, 0.014940, 0.014604, 0.014133, 0.013968, 0.013904, 0.013934, 0.013780, 0.013930, 0.013727, 0.013940, 0.013763, 0.013826, 0.014192, 0.014801, 0.014347, 0.014048, 0.014009, 0.014197, 0.014571, 0.014999, 0.015030, 0.014491, 0.014891, 0.014456, 0.014596, 0.015256, 0.014648, 0.014492, 0.014756, 0.015344, 0.014986, 0.015433, 0.015394, 0.015756, 0.014778, 0.015145, 0.015478, 0.015051, 0.014986, 0.015067, 0.015793, 0.015748, 0.015188, 0.015502, 0.015533, 0.015340, 0.015759, 0.015745, 0.016026, 0.015635, 0.015194, 0.014579, 0.015225, 0.014963, 0.015365, 0.016030, 0.015387, 0.016341, 0.015327, 0.015340, 0.015030, 0.015246, 0.015420, 0.015015, 0.015195, 0.016021, 0.015034, 0.015528, 0.015114, 0.015423, 0.015564, 0.015348, 0.015107, 0.015314, 0.015411, 0.015243, 0.015154, 0.016324, 0.015215, 0.014823, 0.015030, 0.015104, 0.014896, 0.015400, 0.015721, 0.015131, 0.014951, 0.014630, 0.014597, 0.015235, 0.014583, 0.015418, 0.014648, 0.014769, 0.014601, 0.015167, 0.014857, 0.015134, 0.015053, 0.014405, 0.014800, 0.014921, 0.014760, 0.013966, 0.014979, 0.014230, 0.014620, 0.014581, 0.014701, 0.013799, 0.014299, 0.015071, 0.013931, 0.014846, 0.014290, 0.013988, 0.014113, 0.013767, 0.014263, 0.014131, 0.013840, 0.013604, 0.014456, 0.013853, 0.014505, 0.013416, 0.014010, 0.014081, 0.014352, 0.013589, 0.013952, 0.013690, 0.014241, 0.014024, 0.013868, 0.014517, 0.014587, 0.013927, 0.013857, 0.014084, 0.013619, 0.014417, 0.013644, 0.013607, 0.013185, 0.014200, 0.013665, 0.013437, 0.013849, 0.013431, 0.014252, 0.013648, 0.013652, 0.013039, 0.013761, 0.013836, 0.013043, 0.013408, 0.013319, 0.013344, 0.014065, 0.013400, 0.012560, 0.013294, 0.012773, 0.012721, 0.013663, 0.012939, 0.012823, 0.013835, 0.012942, 0.013723, 0.013525};
  for(Int_t i=0; i<400; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5overLHC13d3(){
  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
  fFuncWeight->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL7overLHC10f6a(){
  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
  fFuncWeight->SetParameters(2.41522e+01,4.92146e+00,6.72495e+00,2.5,6.15361e-03,4.78995e+00,-4.29135e-01,3.99421e-01,-1.57220e-02);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL7overLHC12a12(){
  // weight function from the ratio of the LHC12a12 MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,50.);
  fFuncWeight->SetParameters(1.31497e+01,3.30503e+00,3.45594e+00,2.5,2.28642e-02,1.42372e+00,2.32892e-04,5.21986e-02,-2.14236e-01,3.86200e+00);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL7overLHC12a12bis(){
  // weight function from the ratio of the LHC12a12bis MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,50.);
  fFuncWeight->SetParameters(6.54519e+00,2.74007e+00,2.48325e+00,2.5,1.61113e-01,-5.32546e-01,-3.75916e-04,2.38189e-01,-2.17561e-01,2.35975e+00);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL7overLHC13e2fix(){
  // weight function from the ratio of the LHC13e2fix MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,50.);
  fFuncWeight->SetParameters(1.85862e+01,2.48171e+00,3.39356e+00,2.5,1.70426e-02,2.28453e+00,-4.57749e-02,5.84585e-02,-3.19719e-01,4.16789e+00);
  fUseWeight=kTRUE;
}

//________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5overLHC10f6a(){
  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
  fFuncWeight->SetParameters(2.77730e+01,4.78942e+00,7.45378e+00,2.5,9.86255e-02,2.30120e+00,-4.16435e-01,3.43770e-01,-2.29380e-02);
  fUseWeight=kTRUE;
}

//________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL276overLHC10f6a(){
  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,40.);
  fFuncWeight->SetParameters(1.34412e+01,3.20068e+00,5.14481e+00,2.5,7.59405e-04,7.51821e+00,-3.93811e-01,2.16849e-02,-3.37768e-02,2.40308e+00);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5overLHC13d3Lc(){
	// weight function from the ratio of the LHC13d3 MC
	// and FONLL calculations for pPb data
	// using Lc simulated pt distribution
	if(fFuncWeight) delete fFuncWeight;
	fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,20.);
	fFuncWeight->SetParameters(5.94428e+01,1.63585e+01,9.65555e+00,6.71944e+00,8.88338e-02,2.40477e+00,-4.88649e-02,-6.78599e-01,-2.10951e-01);
	fUseWeight=kTRUE;
}
//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5overLHC16i6a(){
	// weight function from the ratio of the LHC16i6a MC
	// and FONLL calculations for pp data
	// using D0 simulated pt distribution
	if(fFuncWeight) delete fFuncWeight;
	fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.2,20.);
	fFuncWeight->SetParameters(2.70821e+03,2.98122e+00,8.67776e+01,1.16611e+00,1.17276e-02,5.41670e+00,6.01099e-02,-2.04524e+00,6.69208e-02);
	fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL7overLHC11b2Lc(){
	// weight function from the ratio of the LHC11b2 MC
	// and FONLL calculations for pp data
	// using Lc simulated pt distribution
	if(fFuncWeight) delete fFuncWeight;
	fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",1.,20.);
	fFuncWeight->SetParameters(2.11879e+02,3.73290e+00,2.01235e+01,1.41508e+00,1.06268e-01,1.86285e+00,-4.52956e-02,-9.90631e-01,-1.31615e+00);
	fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL7overLHC10f7aLc(){
	// weight function from the ratio of the LHC10f7a MC
	// and FONLL calculations for pp data
	// using Lc simulated pt distribution
	if(fFuncWeight) delete fFuncWeight;
	fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",1.,20.);
	fFuncWeight->SetParameters(2.84268e+02,2.18850e+01,2.36298e+01,7.46144e+00,1.69747e-01,1.66993e+00,-5.54726e-02,-1.53869e+00,-1.18404e+00);
	fUseWeight=kTRUE;
}


//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL8overLHC15l2a2(){
	// weight function from the ratio of the LHC15l2a2 MC
	// and FONLL calculations for pp data
	// using D0 simulated pt distribution
	if(fFuncWeight) delete fFuncWeight;
	fFuncWeight=new TF1("funcWeight","(x<28.5)*(([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x))+(x>=28.5)*0.140023",0.15,40);
	fFuncWeight->SetParameters(4.66092e+01, 4.27321e+01, 4.46858e+01, 1.67788e+00, -2.88457e-02, 4.40656e+00, -7.31064e-01, 2.96431e+00, -2.79976e-01);
	fUseWeight=kTRUE;
}


//_________________________________________________________________________
Double_t AliCFTaskVertexingHF::GetWeight(Float_t pt)
{
  //
  // calculating the weight to fill the container
  //

  // FNOLL central:
  // p0 = 1.63297e-01 --> 0.322643
  // p1 = 2.96275e+00
  // p2 = 2.30301e+00
  // p3 = 2.50000e+00

  // PYTHIA
  // p0 = 1.85906e-01 --> 0.36609
  // p1 = 1.94635e+00
  // p2 = 1.40463e+00
  // p3 = 2.50000e+00

  Double_t func1[4] = {0.322643,2.96275,2.30301,2.5};
  Double_t func2[4] = {0.36609,1.94635,1.40463,2.5};

  Double_t dndpt_func1 = dNdptFit(pt,func1);
  if(fUseFlatPtWeight) dndpt_func1 = 1./30.;
  Double_t dndpt_func2 = dNdptFit(pt,func2);
  AliDebug(2,Form("pt = %f, FONLL = %f, Pythia = %f, ratio = %f",pt,dndpt_func1,dndpt_func2,dndpt_func1/dndpt_func2));
  return dndpt_func1/dndpt_func2;
}

//__________________________________________________________________________________________________
Double_t AliCFTaskVertexingHF::dNdptFit(Float_t pt, Double_t* par)
{
  //
  // calculating dNdpt
  //

  Double_t denom =  TMath::Power((pt/par[1]), par[3] );
  Double_t dNdpt = par[0]*pt/TMath::Power(1.+denom, par[2]);

  return dNdpt;
}

//_________________________________________________________________________
Double_t AliCFTaskVertexingHF::GetPtWeightFromHistogram(Float_t pt)
{
  //
  // Using an histogram as weight function
  //  weight = 0 in the range outside the function
  //
  Double_t weight = 0.0;
  Int_t histoNbins = fHistoPtWeight->GetNbinsX();
  Int_t bin2 = fHistoPtWeight->FindBin(pt);
  if( (bin2>0) && (bin2<=histoNbins) ) {
    Int_t bin1=bin2-1;
    Int_t bin3=bin2+1;
    if(bin2==1) bin1=bin2+2;
    if(bin2==histoNbins) bin3=bin2-2;
    Float_t x_1=fHistoPtWeight->GetXaxis()->GetBinCenter(bin1);
    Float_t x_2=fHistoPtWeight->GetXaxis()->GetBinCenter(bin2);
    Float_t x_3=fHistoPtWeight->GetXaxis()->GetBinCenter(bin3);
    Float_t y_1=fHistoPtWeight->GetBinContent(bin1);
    Float_t y_2=fHistoPtWeight->GetBinContent(bin2);
    Float_t y_3=fHistoPtWeight->GetBinContent(bin3);
    Double_t a=( (y_3-y_2)*(x_1-x_2) - (y_1-y_2)*(x_3-x_2) )/( (x_3*x_3-x_2*x_2)*(x_1-x_2) - (x_1*x_1-x_2*x_2)*(x_3-x_2) );
    Double_t b=((y_1-y_2)-a*(x_1*x_1-x_2*x_2))/(x_1-x_2);
    Double_t c=y_3-a*(x_3*x_3)-b*x_3;
    weight = a*pt*pt+b*pt+c;
  }
  return weight;
}

//__________________________________________________________________________________________________
Double_t AliCFTaskVertexingHF::GetZWeight(Float_t z, Int_t runnumber){
  //
  //  calculating the z-vtx weight for the given run range
  //

  if(runnumber>146824 || runnumber<146803) return 1.0;

  Double_t func1[3] = {1.0, -0.5, 6.5 };
  Double_t func2[3] = {1.0, -0.5, 5.5 };

  Double_t dzFunc1 = DodzFit(z,func1);
  Double_t dzFunc2 = DodzFit(z,func2);

  return dzFunc1/dzFunc2;
}

//__________________________________________________________________________________________________
Double_t AliCFTaskVertexingHF::DodzFit(Float_t z, Double_t* par) {

  //
  // Gaussian z-vtx shape
  //
  //gaussian = [0]/TMath::Sqrt(2.*TMath::Pi())/[2]*exp[-(x-[1])*(x-[1])/(2*[2]*[2])]

  Double_t value =  par[0]/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(z-par[1])*(z-par[1])/2./par[2]/par[2]);

  return value;
}
//__________________________________________________________________________________________________
Double_t AliCFTaskVertexingHF::GetNchWeight(Int_t nch){
  //
  //  calculates the Nch weight using the measured and generateed Nch distributions
  //
  if(nch<=0) return 0.;
  if(!fHistoMeasNch || !fHistoMCNch) { AliError("Input histos to evaluate Nch weights missing"); return 0.; }
  Double_t pMeas=fHistoMeasNch->GetBinContent(fHistoMeasNch->FindBin(nch));
  Double_t pMC=fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(nch));
  Double_t weight = pMC>0 ? pMeas/pMC : 0.;
  if(fUseMultRatioAsWeight)  weight = pMC;
  return weight;
}
//__________________________________________________________________________________________________
void AliCFTaskVertexingHF::CreateMeasuredNchHisto(){
  // creates historgam with measured multiplcity distribution in pp 7 TeV collisions (from Eur. Phys. J. C (2010) 68: 345â€“354)
  //
  // for Nch  > 70 the points were obtained with a double NBD distribution fit
  // TF1 *fit1 = new TF1("fit1","[0]*(TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])))*(TMath::Power(([2]/[1]),x))*(TMath::Power((1+([2]/[1])),-x-[1]))"); fit1->SetParameter(0,1.);// normalization constant
  // fit1->SetParameter(1,1.63); // k parameter
  // fit1->SetParameter(2,12.8); // mean multiplicity
  //
  Double_t nchbins[82]={0.50,1.50,2.50,3.50,4.50,5.50,6.50,7.50,8.50,9.50,
                        10.50,11.50,12.50,13.50,14.50,15.50,16.50,17.50,18.50,19.50,
                        20.50,21.50,22.50,23.50,24.50,25.50,26.50,27.50,28.50,29.50,
                        30.50,31.50,32.50,33.50,34.50,35.50,36.50,37.50,38.50,39.50,
                        40.50,41.50,42.50,43.50,44.50,45.50,46.50,47.50,48.50,49.50,
                        50.50,51.50,52.50,53.50,54.50,55.50,56.50,57.50,58.50,59.50,
                        60.50,62.50,64.50,66.50,68.50,70.50,72.50,74.50,76.50,78.50,
                        80.50,82.50,84.50,86.50,88.50,90.50,92.50,94.50,96.50,98.50,
                        100.50,102.50};
  Double_t pch[81]={0.062011,0.072943,0.070771,0.067245,0.062834,0.057383,0.051499,0.04591,0.041109,0.036954,
                    0.03359,0.030729,0.028539,0.026575,0.024653,0.0229,0.021325,0.019768,0.018561,0.017187,
                    0.01604,0.014836,0.013726,0.012576,0.011481,0.010393,0.009502,0.008776,0.008024,0.007452,
                    0.006851,0.006428,0.00594,0.005515,0.005102,0.00469,0.004162,0.003811,0.003389,0.003071,
                    0.002708,0.002422,0.002184,0.001968,0.00186,0.00165,0.001577,0.001387,0.001254,0.001118,
                    0.001037,0.000942,0.000823,0.000736,0.000654,0.000579,0.000512,0.00049,0.00045,0.000355,
                    0.000296,0.000265,0.000193,0.00016,0.000126,0.0000851, 0.0000676,0.0000537,0.0000426, 0.0000338,
                    0.0000268,0.0000213,0.0000166,0.0000133,0.0000106,0.00000837,0.00000662, 0.00000524,0.00000414, 0.00000327,
                    0.00000258};

  if(fHistoMeasNch) delete fHistoMeasNch;
  fHistoMeasNch=new TH1F("hMeaseNch","",81,nchbins);
  for(Int_t i=0; i<81; i++){
    fHistoMeasNch->SetBinContent(i+1,pch[i]);
    fHistoMeasNch->SetBinError(i+1,0.);
  }
}

//__________________________________________________________________________________________________
Bool_t AliCFTaskVertexingHF::ProcessDs(Int_t recoAnalysisCuts) const{
  // processes output of Ds is selected
  Bool_t keep=kFALSE;
  if(recoAnalysisCuts > 0){
    Int_t isKKpi=recoAnalysisCuts&1;
    Int_t ispiKK=recoAnalysisCuts&2;
    Int_t isPhiKKpi=recoAnalysisCuts&4;
    Int_t isPhipiKK=recoAnalysisCuts&8;
    Int_t isK0starKKpi=recoAnalysisCuts&16;
    Int_t isK0starpiKK=recoAnalysisCuts&32;
    if(fDsOption==1){
      if(isKKpi && isPhiKKpi) keep=kTRUE;
      if(ispiKK && isPhipiKK) keep=kTRUE;
    }
    else if(fDsOption==2){
      if(isKKpi && isK0starKKpi) keep=kTRUE;
      if(ispiKK && isK0starpiKK) keep=kTRUE;
    }
    else if(fDsOption==3)keep=kTRUE;
  }
  return keep;
}
//__________________________________________________________________________________________________
Bool_t AliCFTaskVertexingHF::ProcessLctoV0Bachelor(Int_t recoAnalysisCuts) const{

  // processes output of Lc->V0+bnachelor is selected

  Bool_t keep=kFALSE;

  if (recoAnalysisCuts > 0){

    Int_t isK0Sp = recoAnalysisCuts&1;
    Int_t isLambdaBarpi = recoAnalysisCuts&2;
    Int_t isLambdapi = recoAnalysisCuts&4;

    if(fLctoV0bachelorOption == 1){
      if(isK0Sp) keep=kTRUE;
    }
    else if(fLctoV0bachelorOption == 2){
      if(isLambdaBarpi) keep=kTRUE;
    }
    else if(fLctoV0bachelorOption == 4){
      if(isLambdapi) keep=kTRUE;
    }
    else if(fLctoV0bachelorOption == 7) {
      if (isK0Sp || isLambdaBarpi || isLambdapi) keep=kTRUE;
    }
  }
  return keep;
}

//____________________________________________________________________________
TProfile* AliCFTaskVertexingHF::GetEstimatorHistogram(const AliVEvent* event){
  // Get Estimator Histogram from period event->GetRunNumber();
  //
  // If you select SPD tracklets in |eta|<1 you should use type == 1
  //

  Int_t runNo  = event->GetRunNumber();
  Int_t period = -1;   // pp:  0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
                       // pPb: 0-LHC13b, 1-LHC13c

  if (fIsPPbData) {    // setting run numbers for LHC13 if pPb
    if (runNo>195343 && runNo<195484) period = 0;
    if (runNo>195528 && runNo<195678) period = 1;
    if (period<0 || period>1) return 0;
  } else {             //else assume pp 2010
    if(runNo>114930 && runNo<117223) period = 0;
    if(runNo>119158 && runNo<120830) period = 1;
    if(runNo>122373 && runNo<126438) period = 2;
    if(runNo>127711 && runNo<130841) period = 3;
    if(period<0 || period>3) return 0;
  }

  return fMultEstimatorAvg[period];
}
