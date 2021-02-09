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
#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"

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
  fIsPP13TeVData(kFALSE),
  fUseAdditionalCuts(kFALSE),
  fUseCutsForTMVA(kFALSE),
  fUseCascadeTaskForLctoV0bachelor(kFALSE),
  fFillMinimumSteps(kFALSE),
  fCutOnMomConservation(0.00001),
  fMinLeadPtRT(6.0),
  fAODProtection(0),
  fRejectOOBPileUpEvents(kFALSE),
  fKeepOnlyOOBPileupEvents(kFALSE)
{
  //
  //Default ctor
  //
  for(Int_t i=0; i<33; i++) fMultEstimatorAvg[i]=0;
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
  fIsPP13TeVData(kFALSE),
  fUseAdditionalCuts(kFALSE),
  fUseCutsForTMVA(kFALSE),
  fUseCascadeTaskForLctoV0bachelor(kFALSE),
  fFillMinimumSteps(kFALSE),
  fCutOnMomConservation(0.00001),
  fMinLeadPtRT(6.0),
  fAODProtection(0),
  fRejectOOBPileUpEvents(kFALSE),
  fKeepOnlyOOBPileupEvents(kFALSE)
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
  for(Int_t i=0; i<33; i++) fMultEstimatorAvg[i]=0;
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
    for(Int_t i=0; i<33; i++) fMultEstimatorAvg[i]=c.fMultEstimatorAvg[i];
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
  fIsPP13TeVData(c.fIsPP13TeVData),
  fUseAdditionalCuts(c.fUseAdditionalCuts),
  fUseCutsForTMVA(c.fUseCutsForTMVA),
  fUseCascadeTaskForLctoV0bachelor(c.fUseCascadeTaskForLctoV0bachelor),
  fFillMinimumSteps(c.fFillMinimumSteps),
  fCutOnMomConservation(c.fCutOnMomConservation),
  fMinLeadPtRT(c.fMinLeadPtRT),
  fAODProtection(c.fAODProtection),
  fRejectOOBPileUpEvents(c.fRejectOOBPileUpEvents),
  fKeepOnlyOOBPileupEvents(c.fKeepOnlyOOBPileupEvents)
{
  //
  // Copy Constructor
  //
  for(Int_t i=0; i<33; i++) fMultEstimatorAvg[i]=c.fMultEstimatorAvg[i];
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
  for(Int_t i=0; i<33; i++) { if(fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i]; }
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
    case kFalcon:// super fast configuration: only pt_candidate, y, centrality
      fNvar = 4;
      break;
    case kESE:// configuration with variables for ESE analysis (pt,y,centrality,mult,q2)
      fNvar = 6;
      break;
    case kRT: // config with variables for RT analysis (pt,y,mult,RT,deltaphi)
      fNvar = 5;
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
    case kFalcon:// super fast configuration: only pt_candidate, y, centrality
      fNvar = 4;
      break;
    case kESE:// configuration with variables for ESE analysis (pt,y,centrality,mult,q2)
      fNvar = 6;
      break;
    case kRT: // config with variables for RT analysis (pt,y,mult,RT,deltaphi)
      fNvar = 5;
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
    case kFalcon:// super fast configuration: only pt_candidate, y, centrality
      fNvar = 4;
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
    case kFalcon:// super fast configuration: only pt_candidate, y, centrality
      fNvar = 4;
      break;
    case kESE:// configuration with variables for ESE analysis (pt,y,centrality,mult,q2)
      fNvar = 6;
      break;
    case kRT: // config with variables for RT analysis (pt,y,mult,RT,deltaphi)
      fNvar = 5;
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
    case kFalcon:// super fast configuration: only pt_candidate, y, centrality
      fNvar = 4;
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
    case kFalcon:// super fast configuration: only pt_candidate, y, centrality
      fNvar = 4;
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
    case kFalcon:// super fast configuration: only pt_candidate, y, centrality
      fNvar = 4;
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
  TString period[33];
  Int_t nProfiles=33;

  if (fIsPP13TeVData) { //if pp at 13 TeV, use 33 estimator histos
    period[0] = "LHC16d"; period[1] = "LHC16e"; period[2] = "LHC16g"; period[3] = "LHC16h"; period[4] = "LHC16j"; period[5] = "LHC16k"; period[6] = "LHC16l"; period[7] = "LHC16o"; period[8] = "LHC16p";
    period[9] = "LHC17c"; period[10] = "LHC17e"; period[11] = "LHC17f"; period[12] = "LHC17h"; period[13] = "LHC17i"; period[14] = "LHC17j"; period[15] = "LHC17k"; period[16] = "LHC17l"; period[17] = "LHC17m"; period[18] = "LHC17o"; period[19] = "LHC17r";
    period[20] = "LHC18b"; period[21] = "LHC18d"; period[22] = "LHC18e"; period[23] = "LHC18f"; period[24] = "LHC18g"; period[25] = "LHC18h"; period[26] = "LHC18i"; period[27] = "LHC18k"; period[28] = "LHC18l"; period[29] = "LHC18m"; period[30] = "LHC18n"; period[31] = "LHC18o"; period[32] = "LHC18p";
    nProfiles = 33;
  } else if (fIsPPbData) { //if pPb, use only two estimator histos
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

  if(fAODProtection>=0){
        //   Protection against different number of events in the AOD and deltaAOD
        //   In case of discrepancy the event is rejected.
        Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
            // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
            fHistEventsProcessed->Fill(6.5);
            return;
        }
        fHistEventsProcessed->Fill(7.5);
  }

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

  // reject / keep events with simulated OOB pileup
  Bool_t isHijing = kFALSE;
  // check if Hijing is among the generators used in the MC simulation
  TList *lgen = mcHeader->GetCocktailHeaders();
  if(lgen)
  {
    for(Int_t i=0;i<lgen->GetEntries();i++){
      AliGenEventHeader* gh=(AliGenEventHeader*)lgen->At(i);
      TString genname=gh->GetName();
      if(genname.Contains("Hijing"))
      {
        isHijing = kTRUE;
        break;
      }
    }
  }

  if(isHijing) {
    Bool_t isPileUp = AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader, "Hijing");
    if(isPileUp && fRejectOOBPileUpEvents)
      return;
    if(!isPileUp && fKeepOnlyOOBPileupEvents)
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
  if(fUseNchWeight) {
    Int_t nChargedMCPhysicalPrimary=AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(mcArray,-1.0,1.0);
    if(!fUseTrackletsWeight) {
    	fWeight *= GetNchWeight(nChargedMCPhysicalPrimary);
    	AliDebug(2,Form("Using Nch weights, Mult=%d Weight=%f\n",nChargedMCPhysicalPrimary,fWeight));
    } else {
    	fWeight *= GetNchWeight(static_cast<Int_t>(nTracklets));
    	AliDebug(2,Form("Using Nch weights (with tracklets), TrklMult=%d Weight=%f\n",nTracklets,fWeight));
    }
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

  Double_t q2=0;
  if(fConfiguration==kESE) {
    //set q2 in case of kESE configuration
    q2=ComputeTPCq2(aodEvent,mcHeader,-0.8,0.8,0.2,5.);
    cfVtxHF->Setq2Value(q2);

    //set track array in case of kESE configuration
    cfVtxHF->SetTrackArray(aodEvent->GetTracks());
  }

  Double_t rtval=-1.;
  if (fConfiguration==kRT) {
     //do RT determination if RT analysis
     rtval = CalculateRTValue(aodEvent,mcHeader,cfVtxHF);
     cfVtxHF->SetRTValue(rtval);
  }

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

    // PILEUP protection for PbPb2018: remove particles from pileup events in efficiency computation
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iPart, mcHeader, mcArray)) {
      AliDebug(2, Form("Check on the out-of-bunch pile-up wrong for particle %d!!!", iPart));
      fHistEventsProcessed->Fill(8.5);
      continue;
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
      if(!fFillMinimumSteps)
        fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGenLimAccNoAcc, fWeight);
	    icountGenLimAccNoAcc++;
	    AliDebug(3,"MC Lim Acc No Acc container filled\n");
	  }
	}
      }

      //MC
      if(!fFillMinimumSteps)
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
          if(!fFillMinimumSteps)
            fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepVertex, fWeight) ;
          AliDebug(3,"Vertex cut passed and container filled\n");
          icountVertex++;

          //mc Refit requirement
          Bool_t mcRefitStep = cfVtxHF->MCRefitStep(aodEvent, &trackCuts[0]);
          if (mcRefitStep){
            if(!fFillMinimumSteps)
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
    fHistEventsProcessed->Fill(2.5);
    AliAODRecoDecayHF* charmCandidate=0x0;
    switch (fDecayChannel){
    case 2:{
      charmCandidate = (AliAODRecoDecayHF2Prong*)arrayBranch->At(iCandid);
      if(charmCandidate->GetIsFilled()!=0) fHistEventsProcessed->Fill(3.5);
      if(!vHF->FillRecoCand(aodEvent,(AliAODRecoDecayHF2Prong*)charmCandidate)){
	fHistEventsProcessed->Fill(5.5);
	continue;
      }else{
	fHistEventsProcessed->Fill(4.5);
      }
      break;
    }
    case 21:{
      charmCandidate = (AliAODRecoCascadeHF*)arrayBranch->At(iCandid);
      if(charmCandidate->GetIsFilled()!=0) fHistEventsProcessed->Fill(3.5);
      if(!vHF->FillRecoCasc(aodEvent,((AliAODRecoCascadeHF*)charmCandidate),kTRUE)){
	fHistEventsProcessed->Fill(5.5);
	continue; //DStar
      }else{
	fHistEventsProcessed->Fill(4.5);
      }
      break;
    }
    case 22:{
      charmCandidate = (AliAODRecoCascadeHF*)arrayBranch->At(iCandid);
      if(charmCandidate->GetIsFilled()!=0) fHistEventsProcessed->Fill(3.5);
      if(!vHF->FillRecoCasc(aodEvent,((AliAODRecoCascadeHF*)charmCandidate),kFALSE)){
	fHistEventsProcessed->Fill(5.5);
	continue; //Cascade
      }else{
	fHistEventsProcessed->Fill(4.5);
      }
      break;
    }
    case 31:
    case 32:
    case 33:{
      charmCandidate = (AliAODRecoDecayHF3Prong*)arrayBranch->At(iCandid);
      if(charmCandidate->GetIsFilled()!=0) fHistEventsProcessed->Fill(3.5);
      if(!vHF->FillRecoCand(aodEvent,(AliAODRecoDecayHF3Prong*)charmCandidate)){
	fHistEventsProcessed->Fill(5.5);
	continue;
      }else{
	fHistEventsProcessed->Fill(4.5);
      }
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
        if(!fFillMinimumSteps)
          fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed, fWeight) ;
        icountReco++;
        AliDebug(3,"Reco step  passed and container filled\n");

        //Reco in the acceptance -- take care of UNFOLDING!!!!
        Bool_t recoAcceptanceStep = cfVtxHF->RecoAcceptStep(&trackCuts[0]);
        if (recoAcceptanceStep) {
          if(!fFillMinimumSteps)
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
            if(!fFillMinimumSteps)
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
  else if(fConfiguration == kCheetah){
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
  else if(fConfiguration == kFalcon){
    //h = new TH1D[3][12];
    nvarToPlot = 4;
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
  else if(fConfiguration == kESE){
    //h = new TH1D[3][12];
    nvarToPlot = 6;
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
  else if(fConfiguration == kRT) {
     nvarToPlot = 5;
     for (Int_t ih = 0; ih < 3; ih++) {
        h[ih] = new TH1D[nvarToPlot];
     }
     for (Int_t iC=0;iC<nvarToPlot;iC++){
        //MC-level
        h[0][iC] = *(cont->ShowProjection(iC,0));
        //MC-Acceptance level
        h[1][iC] = *(cont->ShowProjection(iC,1));
        //Reco-level
        h[2][iC] = *(cont->ShowProjection(iC,4));
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
  else if(fConfiguration == kCheetah){
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
  else if(fConfiguration == kFalcon){
    //nvarToPlot = 4;
    titles = new TString[nvarToPlot];
    if (fDecayChannel==22) {
      titles[0]="p_{T}(#Lambda_{c}) [GeV/c]";
      titles[1]="y(#Lambda_{c})";
      titles[2]="centrality";
      titles[3]="multiplicity";
    } else {
      titles[0]="pT_candidate (GeV/c)";
      titles[1]="rapidity";
      titles[2]="centrality";
      titles[3]="multiplicity";
    }
  }
  else if(fConfiguration == kESE){
    //nvarToPlot = 4;
    titles = new TString[nvarToPlot];
    titles[0]="pT_candidate (GeV/c)";
    titles[1]="rapidity";
    titles[2]="centrality";
    titles[3]="multiplicity";
    titles[4]="N_{tracks} (R<0.4)";
    titles[5]="q_{2}";
  }
  else if(fConfiguration == kRT) {
     //nvarToPlot =  5;
     titles = new TString[nvarToPlot];
     titles[0]="pT_candidate (GeV/c)";
     titles[1]="rapidity";
     titles[2]="multiplicity";
     titles[3]="RT";
     titles[4]="deltaPhi";
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
    if(iVar>=fNvar) continue;
    c1->cd(iPad++);
    h[0][iVar].DrawCopy("p");
    c1->cd(iPad++);
    h[1][iVar].DrawCopy("p");
    c1->cd(iPad++);
    h[2][iVar].DrawCopy("p");
  }

  if (fConfiguration == kSnail || fConfiguration == kCheetah){
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
  fHistEventsProcessed = new TH1I(nameoutput,"",9,0,9) ;
  fHistEventsProcessed->GetXaxis()->SetBinLabel(1,"Events processed (all)");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(2,"Events analyzed (after selection)");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(3,"Candidates processed (all)");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(4,"Candidates already filled");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(5,"Candidates OK in FillRecoCand");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(6,"Candidates failing in FillRecoCand");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(7,"AOD/dAOD mismatch");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(8,"AOD/dAOD #events ok");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(9,"Candidates from OOB pile-up");

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
void AliCFTaskVertexingHF::SetPtWeightsFDFromFONLL5anddataoverLHC19c3a(){
    // non-prompt D0 weight function from the ratio of the LHC19c3a MC 0-10%
    // 0-1.5 GeV/c using FONLL*Raa and 1.5-50 GeV/c using data

    if(fHistoPtWeight) delete fHistoPtWeight;
    fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
    Float_t binc[500]={ 0.527663 ,0.616222 ,0.694634 ,0.761106 ,0.815208 ,0.857483 ,0.889068 ,0.911394 ,0.925968 ,0.934232 ,0.937491 ,0.936874 ,0.882565 ,0.81223 ,0.754265 ,0.705777 ,0.664635 ,0.629247 ,0.598402 ,0.571173 ,0.546839 ,0.524839 ,0.504729 ,0.486156 ,0.468842 ,0.452563 ,0.437141 ,0.422434 ,0.408328 ,0.394732 ,0.381576 ,0.368804 ,0.356371 ,0.344245 ,0.332401 ,0.320819 ,0.309487 ,0.298396 ,0.287539 ,0.276914 ,0.266519 ,0.256354 ,0.246421 ,0.236721 ,0.227257 ,0.21803 ,0.209043 ,0.200297 ,0.191795 ,0.183538 ,0.175526 ,0.16776 ,0.16024 ,0.152964 ,0.145933 ,0.139143 ,0.132593 ,0.12628 ,0.1202 ,0.114352 ,0.108729 ,0.105356 ,0.104384 ,0.103368 ,0.102312 ,0.103163 ,0.100649 ,0.0981966 ,0.0958041 ,0.0934698 ,0.0911924 ,0.0889705 ,0.0868027 ,0.0846878 ,0.0826243 ,0.0806112 ,0.0786471 ,0.0767308 ,0.0748613 ,0.0730373 ,0.0712577 ,0.0695215 ,0.0678276 ,0.066175 ,0.0645627 ,0.0629896 ,0.0614548 ,0.0599575 ,0.0584966 ,0.0570713 ,0.0556808 ,0.0543241 ,0.0530005 ,0.0517092 ,0.0504493 ,0.0492201 ,0.0480208 ,0.0468508 ,0.0457093 ,0.0445956 ,0.043509 ,0.0424489 ,0.0414146 ,0.0404055 ,0.0394211 ,0.0384606 ,0.0375235 ,0.0366092 ,0.0357172 ,0.034847 ,0.0339979 ,0.0331696 ,0.0323614 ,0.0315729 ,0.0308036 ,0.0300531 ,0.0293208 ,0.0286064 ,0.0279094 ,0.0271751 ,0.0271134 ,0.0270518 ,0.0269903 ,0.0269289 ,0.0268677 ,0.026661 ,0.0263976 ,0.0261369 ,0.0258786 ,0.025623 ,0.0253698 ,0.0251192 ,0.024871 ,0.0246253 ,0.024382 ,0.0241411 ,0.0239026 ,0.0236665 ,0.0234327 ,0.0232012 ,0.022972 ,0.022745 ,0.0225203 ,0.0222978 ,0.0220775 ,0.0218594 ,0.0216435 ,0.0214296 ,0.0212179 ,0.0210083 ,0.0208007 ,0.0205952 ,0.0203918 ,0.0201903 ,0.0199908 ,0.0197933 ,0.0195978 ,0.0194042 ,0.0192125 ,0.0190227 ,0.0188347 ,0.0186487 ,0.0184644 ,0.018282 ,0.0181014 ,0.0179225 ,0.0177455 ,0.0175702 ,0.0173966 ,0.0172247 ,0.0170545 ,0.0168861 ,0.0167192 ,0.016554 ,0.0163905 ,0.0162286 ,0.0160682 ,0.0159095 ,0.0157523 ,0.0155967 ,0.0154426 ,0.01529 ,0.015139 ,0.0149894 ,0.0148413 ,0.0146947 ,0.0145495 ,0.0144058 ,0.0142635 ,0.0141226 ,0.013983 ,0.0138449 ,0.0137081 ,0.0135727 ,0.0134386 ,0.0133058 ,0.0131744 ,0.0130442 ,0.0129153 ,0.0127877 ,0.0126614 ,0.0125363 ,0.0124125 ,0.0122898 ,0.0121684 ,0.0120482 ,0.0119292 ,0.0118113 ,0.0116946 ,0.0115791 ,0.0114647 ,0.0113514 ,0.0112393 ,0.0111282 ,0.0110183 ,0.0109095 ,0.0108017 ,0.010695 ,0.0105893 ,0.0104847 ,0.0103811 ,0.0102785 ,0.010177 ,0.0100764 ,0.0099769 ,0.00987833 ,0.00978074 ,0.00968411 ,0.00958844 ,0.00949371 ,0.00939991 ,0.00930705 ,0.0092151 ,0.00912406 ,0.00903392 ,0.00894467 ,0.0088563 ,0.0087688 ,0.00868217 ,0.0085964 ,0.00851147 ,0.00842738 ,0.00834412 ,0.00826169 ,0.00818007 ,0.00809925 ,0.00801923 ,0.00794001 ,0.00786157 ,0.0077839 ,0.007707 ,0.00763086 ,0.00755547 ,0.00748082 ,0.00740692 ,0.00733374 ,0.00726129 ,0.00718955 ,0.00711852 ,0.00704819 ,0.00697856 ,0.00690962 ,0.00684135 ,0.00677376 ,0.00670684 ,0.00664058 ,0.00657498 ,0.00651002 ,0.0064457 ,0.00638202 ,0.00631897 ,0.00625654 ,0.00619473 ,0.00613353 ,0.00607294 ,0.00601294 ,0.00595353 ,0.00589472 ,0.00583648 ,0.00577882 ,0.00572173 ,0.0056652 ,0.00560923 ,0.00555381 ,0.00549894 ,0.00544462 ,0.00539083 ,0.00533757 ,0.00528484 ,0.00523263 ,0.00518093 ,0.00512974 ,0.00507907 ,0.00502889 ,0.0049792 ,0.00493001 ,0.00488131 ,0.00483308 ,0.00478533 ,0.00473806 ,0.00469125 ,0.0046449 ,0.00459901 ,0.00455358 ,0.00450859 ,0.00446405 ,0.00441994 ,0.00437628 ,0.00433304 ,0.00429023 ,0.00424785 ,0.00420588 ,0.00416433 ,0.00412319 ,0.00408245 ,0.00404212 ,0.00400219 ,0.00396265 ,0.0039235 ,0.00388474 ,0.00384636 ,0.00380836 ,0.00377073 ,0.00373348 ,0.0036966 ,0.00366007 ,0.00362391 ,0.00358811 ,0.00355266 ,0.00351757 ,0.00348281 ,0.00344841 ,0.00341434 ,0.00338061 ,0.00334721 ,0.00331414 ,0.0032814 ,0.00324898 ,0.00321688 ,0.0031851 ,0.00315363 ,0.00312247 ,0.00309163 ,0.00306108 ,0.00303084 ,0.0030009 ,0.00297125 ,0.0029419 ,0.00291283 ,0.00288405 ,0.00285556 ,0.00282735 ,0.00279942 ,0.00277176 ,0.00274438 ,0.00271726 ,0.00269042 ,0.00266384 ,0.00263752 ,0.00261146 ,0.00258566 ,0.00256012 ,0.00253483 ,0.00250978 ,0.00248499 ,0.00246044 ,0.00243613 ,0.00241206 ,0.00238823 ,0.00236464 ,0.00234128 ,0.00231815 ,0.00229524 ,0.00227257 ,0.00225012 ,0.00222789 ,0.00220588 ,0.00218408 ,0.00216251 ,0.00214114 ,0.00211999 ,0.00209904 ,0.00207831 ,0.00205777 ,0.00203744 ,0.00201732 ,0.00199739 ,0.00197765 ,0.00195811 ,0.00193877 ,0.00191961 ,0.00190065 ,0.00188187 ,0.00186328 ,0.00184487 ,0.00182665 ,0.0018086 ,0.00179073 ,0.00177304 ,0.00175552 ,0.00173818 ,0.00172101 ,0.001704 ,0.00168717 ,0.0016705 ,0.001654 ,0.00163766 ,0.00162148 ,0.00160546 ,0.0015896 ,0.00157389 ,0.00155834 ,0.00154295 ,0.00152771 ,0.00151261 ,0.00149767 ,0.00148287 ,0.00146822 ,0.00145372 ,0.00143936 ,0.00142513 ,0.00141106 ,0.00139711 ,0.00138331 ,0.00136965 ,0.00135611 ,0.00134272 ,0.00132945 ,0.00131632 ,0.00130331 ,0.00129044 ,0.00127769 ,0.00126506 ,0.00125257 ,0.00124019 ,0.00122794 ,0.00121581 ,0.0012038 ,0.0011919 ,0.00118013 ,0.00116847 ,0.00115693 ,0.0011455 ,0.00113418 ,0.00112297 ,0.00111188 ,0.00110089 ,0.00109002 ,0.00107925 ,0.00106859 ,0.00105803 ,0.00104758 ,0.00103723 ,0.00102698 ,0.00101683 ,0.00100679 ,0.000996842 ,0.000986994 ,0.000977243 ,0.000967588 ,0.000958029 ,0.000948564 ,0.000939193 ,0.000929914 ,0.000920727 ,0.000911631 ,0.000902624 ,0.000893707 ,0.000884877 ,0.000876135 ,0.000867479 ,0.000858909 ,0.000850424 ,0.000842022 ,0.000833703 ,0.000825467 ,0.000817311 ,0.000809237 ,0.000801242 ,0.000793326 ,0.000785488 ,0.000777728 ,0.000770045 ,0.000762437 ,0.000754905 ,0.000747447 ,0.000740062 ,0.000732751 ,0.000725512 ,0.000718344 ,0.000711247 ,0.00070422 ,0.000697263 ,0.000690374 ,0.000683554 ,0.000676801 ,0.000670114 ,0.000663494 ,0.000656939 ,0.000650449 };
    for(Int_t i=0; i<500; i++){
        fHistoPtWeight->SetBinContent(i+1,binc[i]);
    }
    //SetWeightHistogram();
    fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFDFromFONLL5anddataoverLHC19c3b(){
    // non-prompt D0 weight function from the ratio of the LHC19c3a MC 30-50%
    // 0-1.5 GeV/c using FONLL*Raa and 1.5-50 GeV/c using data

    if(fHistoPtWeight) delete fHistoPtWeight;
    fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
    Float_t binc[500]={0.60507345, 0.72757414, 0.83929434, 0.93633835, 1.0169410, 1.0810175, 1.1296302, 1.1645060, 1.1876600, 1.2011356, 1.2068452, 1.2064864, 1.2015105, 1.1931229, 1.1823027, 1.1698298, 1.1084997, 1.0566777, 1.0117197, 0.97224002, 0.93715379, 0.90560414, 0.87690916, 0.85052290, 0.82600624, 0.80300486, 0.78123236, 0.76045723, 0.74049255, 0.72118792, 0.70242292, 0.68410187, 0.66614952, 0.64850757, 0.63113170, 0.61398918, 0.59705677, 0.58031908, 0.56376703, 0.54739672, 0.53120830, 0.51520515, 0.49939307, 0.48377973, 0.46837404, 0.45318577, 0.43822516, 0.42350256, 0.40902826, 0.39481218, 0.47709550, 0.46266798, 0.44867675, 0.43510861, 0.42195079, 0.40919086, 0.39681679, 0.38481692, 0.37317993, 0.36189484, 0.35095103, 0.34033815, 0.33004621, 0.32006551, 0.31038662, 0.30100043, 0.29189808, 0.28307098, 0.27451082, 0.26620952, 0.25815926, 0.25035244, 0.24278170, 0.23543990, 0.22832012, 0.22141565, 0.21471996, 0.20822676, 0.20192992, 0.19582349, 0.18990172, 0.18415903, 0.17859000, 0.17318938, 0.16795208, 0.16287315, 0.15794781, 0.15317142, 0.14853947, 0.14404758, 0.13969154, 0.13546722, 0.13137065, 0.12739795, 0.12354540, 0.11980934, 0.11618627, 0.11267276, 0.10926550, 0.10596127, 0.10275697, 0.099649565, 0.096636130, 0.093713821, 0.090879885, 0.088131647, 0.085466517, 0.082881981, 0.080375603, 0.077945018, 0.076174685, 0.076159299, 0.076143916, 0.076128537, 0.076113160, 0.076097787, 0.076082417, 0.076067050, 0.076051686, 0.076036325, 0.076020967, 0.076005612, 0.075990261, 0.075974912, 0.075959567, 0.075944225, 0.075928886, 0.075913549, 0.075898216, 0.075882887, 0.071260953, 0.070592453, 0.069930223, 0.069274206, 0.068624344, 0.067980577, 0.067342850, 0.066711105, 0.066085287, 0.065465340, 0.064851208, 0.064242837, 0.063640174, 0.063043164, 0.062451755, 0.061865894, 0.061285528, 0.060710608, 0.060141080, 0.059576895, 0.059018003, 0.058464354, 0.057915899, 0.057372588, 0.056834375, 0.056301211, 0.055773048, 0.055249840, 0.054731540, 0.054218102, 0.053709481, 0.053205631, 0.052706508, 0.052212067, 0.051722264, 0.051237057, 0.050756401, 0.050280254, 0.049808574, 0.049341319, 0.048878447, 0.048419917, 0.047965689, 0.047515722, 0.047069976, 0.046628411, 0.046190989, 0.045757670, 0.045328417, 0.044903190, 0.044481952, 0.044064666, 0.043651294, 0.043241801, 0.042836148, 0.042434302, 0.042036225, 0.041641882, 0.041251238, 0.040864260, 0.040480911, 0.040101159, 0.039724969, 0.039352308, 0.038983144, 0.038617442, 0.038255171, 0.037896298, 0.037540792, 0.037188621, 0.036839754, 0.036494160, 0.036151807, 0.035812666, 0.035476707, 0.035143899, 0.034814213, 0.034487621, 0.034164091, 0.033843597, 0.033526110, 0.033211601, 0.032900042, 0.032591406, 0.032285665, 0.031982793, 0.031682762, 0.031385545, 0.031091117, 0.030799450, 0.030510520, 0.030224300, 0.029940765, 0.029659890, 0.029381650, 0.029106020, 0.028832976, 0.028562493, 0.028294548, 0.028029116, 0.027766175, 0.027505699, 0.027247668, 0.026992057, 0.026738844, 0.026488006, 0.026239522, 0.025993368, 0.025749524, 0.025507967, 0.025268676, 0.025031630, 0.024796808, 0.024564189, 0.024333751, 0.024105476, 0.023879342, 0.023655329, 0.023433418, 0.023213589, 0.022995822, 0.022780097, 0.022566397, 0.022354701, 0.022144991, 0.021937248, 0.021731455, 0.021527591, 0.021325641, 0.021125584, 0.020927405, 0.020731084, 0.020536606, 0.020343951, 0.020153104, 0.019964048, 0.019776764, 0.019591238, 0.019407452, 0.019225391, 0.019045037, 0.018866375, 0.018689389, 0.018514064, 0.018340383, 0.018168331, 0.017997894, 0.017829055, 0.017661800, 0.017496115, 0.017331983, 0.017169392, 0.017008325, 0.016848770, 0.016690711, 0.016534135, 0.016379028, 0.016225376, 0.016073165, 0.015922383, 0.015773014, 0.015625048, 0.015478469, 0.015333265, 0.015189423, 0.015046931, 0.014905775, 0.014765944, 0.014627424, 0.014490204, 0.014354271, 0.014219613, 0.014086219, 0.013954076, 0.013823172, 0.013693497, 0.013565038, 0.013437784, 0.013311724, 0.013186846, 0.013063140, 0.012940595, 0.012819199, 0.012698941, 0.012579812, 0.012461801, 0.012344896, 0.012229089, 0.012114367, 0.012000722, 0.011888143, 0.011776620, 0.011666144, 0.011556703, 0.011448290, 0.011340893, 0.011234504, 0.011129113, 0.011024710, 0.010921287, 0.010818834, 0.010717343, 0.010616803, 0.010517207, 0.010418544, 0.010320808, 0.010223988, 0.010128077, 0.010033065, 0.0099389446, 0.0098457071, 0.0097533443, 0.0096618480, 0.0095712100, 0.0094814223, 0.0093924769, 0.0093043658, 0.0092170814, 0.0091306158, 0.0090449613, 0.0089601103, 0.0088760553, 0.0087927888, 0.0087103035, 0.0086285920, 0.0085476470, 0.0084674613, 0.0083880279, 0.0083093396, 0.0082313895, 0.0081541707, 0.0080776762, 0.0080018994, 0.0079268334, 0.0078524716, 0.0077788074, 0.0077058342, 0.0076335457, 0.0075619352, 0.0074909965, 0.0074207233, 0.0073511094, 0.0072821485, 0.0072138345, 0.0071461613, 0.0070791230, 0.0070127136, 0.0069469272, 0.0068817580, 0.0068172000, 0.0067532477, 0.0066898954, 0.0066271373, 0.0065649680, 0.0065033819, 0.0064423736, 0.0063819375, 0.0063220684, 0.0062627610, 0.0062040099, 0.0061458099, 0.0060881559, 0.0060310428, 0.0059744655, 0.0059184189, 0.0058628981, 0.0058078981, 0.0057534141, 0.0056994412, 0.0056459746, 0.0055930096, 0.0055405415, 0.0054885655, 0.0054370772, 0.0053860719, 0.0053355450, 0.0052854921, 0.0052359088, 0.0051867907, 0.0051381333, 0.0050899323, 0.0050421836, 0.0049948827, 0.0049480256, 0.0049016081, 0.0048556260, 0.0048100753, 0.0047649519, 0.0047202518, 0.0046759710, 0.0046321056, 0.0045886517, 0.0045456055, 0.0045029631, 0.0044607207, 0.0044188745, 0.0043774210, 0.0043363563, 0.0042956769, 0.0042553790, 0.0042154592, 0.0041759139, 0.0041367396, 0.0040979327, 0.0040594899, 0.0040214078, 0.0039836829, 0.0039463119, 0.0039092914, 0.0038726183, 0.0038362892, 0.0038003009, 0.0037646502, 0.0037293339, 0.0036943489, 0.0036596922, 0.0036253605, 0.0035913509, 0.0035576604, 0.0035242859, 0.0034912245, 0.0034584733, 0.0034260292, 0.0033938896, 0.0033620514, 0.0033305120, 0.0032992684, 0.0032683179, 0.0032376577, 0.0032072852, 0.0031771976, 0.0031473922, 0.0031178664, 0.0030886177, 0.0030596433, 0.0030309407, 0.0030025074, 0.0029743408, 0.0029464385, 0.0029187978, 0.0028914165, 0.0028642921, 0.0028374221, 0.0028108042, 0.0027844360, 0.0027583151, 0.0027324393, 0.0027068063, 0.0026814136, 0.0026562592, 0.0026313408, 0.0026066562, 0.0025822031, 0.0025579794, 0.0025339829, 0.0025102116, 0.0024866632, 0.0024633358, 0.0024402272, 0.0024173353, 0.0023946583, 0.0023721939, 0.0023499403, 0.0023278955, 0.0023060575, 0.0022844243, 0.0022629940, 0.0022417649, 0.0022207348, 0.0021999020};
    for(Int_t i=0; i<500; i++){
        fHistoPtWeight->SetBinContent(i+1,binc[i]);
    }
    //SetWeightHistogram();
    fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5andLBToverLHC16i2a(){
    // weight function from the ratio of the LHC16i2a MC
    // using FONLL*Raa from Xu.Cao.Bass(LBT model)

    if(fHistoPtWeight) delete fHistoPtWeight;
    fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
    Float_t binc[500]={ -3.491009, -2.151267, -1.169817, -0.462113, 0.038852, 0.385901, 0.620406, 0.774498, 0.872869, 0.934235, 0.972536, 0.997905, 1.017465, 1.035981, 1.056389, 1.080233, 1.108019, 1.139503, 1.173928, 1.210210, 1.247090, 1.283254, 1.317421, 1.348413, 1.375204, 1.396948, 1.413000, 1.422920, 1.426469, 1.423596, 1.414424, 1.399224, 1.378394, 1.352429, 1.321897, 1.287417, 1.249627, 1.209171, 1.166675, 1.122735, 1.077901, 1.032675, 0.987496, 0.942748, 0.898751, 0.855766, 0.813999, 0.773607, 0.734700, 0.697350, 0.661596, 0.627448, 0.594895, 0.563911, 0.534455, 0.506477, 0.479923, 0.454732, 0.430845, 0.408201, 0.386738, 0.366399, 0.347127, 0.328866, 0.311564, 0.295172, 0.279641, 0.264928, 0.250988, 0.237781, 0.248506, 0.240198, 0.232287, 0.224750, 0.217563, 0.210707, 0.204160, 0.197906, 0.191927, 0.186209, 0.180735, 0.175494, 0.170472, 0.165657, 0.161039, 0.156606, 0.152350, 0.148261, 0.144331, 0.140552, 0.136916, 0.133416, 0.130046, 0.126800, 0.123671, 0.120654, 0.117744, 0.114937, 0.112226, 0.109608, 0.107080, 0.104636, 0.102273, 0.099988, 0.097778, 0.095638, 0.093567, 0.091561, 0.089618, 0.087735, 0.085910, 0.084141, 0.082424, 0.080759, 0.079143, 0.077574, 0.076051, 0.074571, 0.073134, 0.071737, 0.070379, 0.069059, 0.067775, 0.066526, 0.065311, 0.064128, 0.062977, 0.061856, 0.060764, 0.059700, 0.058664, 0.057654, 0.056670, 0.055710, 0.054775, 0.053862, 0.052972, 0.052103, 0.051255, 0.050428, 0.049620, 0.048831, 0.048061, 0.047308, 0.046573, 0.045855, 0.045152, 0.044466, 0.043795, 0.043139, 0.042498, 0.041870, 0.041256, 0.040656, 0.040068, 0.039493, 0.038930, 0.038379, 0.037839, 0.037310, 0.036793, 0.036286, 0.035789, 0.035302, 0.034825, 0.034358, 0.033900, 0.033451, 0.033010, 0.032578, 0.032155, 0.031740, 0.031332, 0.030932, 0.030540, 0.030155, 0.029778, 0.029407, 0.029043, 0.028686, 0.028335, 0.027991, 0.027653, 0.027320, 0.026994, 0.026674, 0.026359, 0.026049, 0.025745, 0.025447, 0.025153, 0.024864, 0.024581, 0.024302, 0.024027, 0.023758, 0.023492, 0.023232, 0.022975, 0.022723, 0.022474, 0.022230, 0.021990, 0.021754, 0.021521, 0.021292, 0.021066, 0.020845, 0.020626, 0.020411, 0.020199, 0.019991, 0.019786, 0.019584, 0.019385, 0.019189, 0.018995, 0.018805, 0.018618, 0.018433, 0.018251, 0.018072, 0.017895, 0.017721, 0.017549, 0.017380, 0.017214, 0.017049, 0.016887, 0.016727, 0.016570, 0.016415, 0.016262, 0.016111, 0.015962, 0.015815, 0.015670, 0.015527, 0.015386, 0.015247, 0.015110, 0.014974, 0.014841, 0.014709, 0.014579, 0.014450, 0.014324, 0.014199, 0.014075, 0.013953, 0.013833, 0.013714, 0.013597, 0.013481, 0.013367, 0.013254, 0.013143, 0.013033, 0.012924, 0.012817, 0.012711, 0.012606, 0.012503, 0.012401, 0.012300, 0.012200, 0.012102, 0.012004, 0.011908, 0.011813, 0.011719, 0.011626, 0.011535, 0.011444, 0.011355, 0.011266, 0.011179, 0.011092, 0.011007, 0.010922, 0.010839, 0.010756, 0.010675, 0.010594, 0.010514, 0.010435, 0.010357, 0.010280, 0.010204, 0.010128, 0.010053, 0.009980, 0.009907, 0.009834, 0.009763, 0.009692, 0.009622, 0.009553, 0.009485, 0.009417, 0.009350, 0.009284, 0.009218, 0.009153, 0.009089, 0.009026, 0.008963, 0.008901, 0.008839, 0.008778, 0.008718, 0.008658, 0.008599, 0.008540, 0.008483, 0.008425, 0.008368, 0.008312, 0.008257, 0.008202, 0.008147, 0.008093, 0.008040, 0.007987, 0.007934, 0.007883, 0.007831, 0.007780, 0.007730, 0.007680, 0.007631, 0.007582, 0.007533, 0.007485, 0.007438, 0.007390, 0.007344, 0.007298, 0.007252, 0.007206, 0.007161, 0.007117, 0.007073, 0.007029, 0.006986, 0.006943, 0.006900, 0.006858, 0.006817, 0.006775, 0.006734, 0.006694, 0.006653, 0.006613, 0.006574, 0.006535, 0.006496, 0.006457, 0.006419, 0.006381, 0.006344, 0.006307, 0.006270, 0.006233, 0.006197, 0.006161, 0.006126, 0.006090, 0.006055, 0.006021, 0.005986, 0.005952, 0.005918, 0.005885, 0.005852, 0.005819, 0.005786, 0.005754, 0.005722, 0.005690, 0.005658, 0.005627, 0.005596, 0.005565, 0.005534, 0.005504, 0.005474, 0.005444, 0.005415, 0.005385, 0.005356, 0.005327, 0.005299, 0.005270, 0.005242, 0.005214, 0.005186, 0.005159, 0.005132, 0.005104, 0.005078, 0.005051, 0.005025, 0.004998, 0.004972, 0.004946, 0.004921, 0.004895, 0.004870, 0.004845, 0.004820, 0.004796, 0.004771, 0.004747, 0.004723, 0.004699, 0.004675, 0.004651, 0.004628, 0.004605, 0.004582, 0.004559, 0.004536, 0.004514, 0.004491, 0.004469, 0.004447, 0.004425, 0.004404, 0.004382, 0.004361, 0.004339, 0.004318, 0.004297, 0.004277, 0.004256, 0.004235, 0.004215, 0.004195, 0.004175, 0.004155, 0.004135, 0.004116, 0.004096, 0.004077, 0.004058, 0.004039, 0.004020, 0.004001, 0.003982, 0.003964, 0.003945, 0.003927, 0.003909, 0.003891, 0.003873, 0.003855, 0.003837, 0.003820, 0.003802, 0.003785, 0.003768, 0.003751, 0.003734, 0.003717, 0.003700, 0.003684, 0.003667, 0.003651, 0.003634, 0.003618, 0.003602, 0.003586, 0.003570, 0.003555, 0.003539, 0.003523, 0.003508, 0.003493, 0.003477, 0.003462, 0.003447, 0.003432, 0.003417, 0.003403, 0.003388, 0.003373, 0.003359, 0.003345, 0.003330, 0.003316, 0.003302, 0.003288, 0.003274, 0.003260, 0.003246, 0.003233, 0.003219, 0.003206, 0.003192};
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
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5andDplusdataoverLHC16i2a(){
  // weight function from the ratio of the LHC16i2a MC
  // 2.5-14 GeV/c using data and 0-2.5, 14-50 GeV/c using FONLL calculations

  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
  Float_t binc[500]={ 1.977386, 1.092380, 0.982069, 0.990772, 1.031869, 1.075578, 1.121814, 1.160738, 1.194926, 1.223790, 1.248790, 1.265487,
    1.280118, 1.292664, 1.302578, 1.308750, 1.316914, 1.321336, 1.325862, 1.335865, 1.330973, 1.338854, 1.347197, 1.355202,
    1.345868, 1.306835, 1.217713, 1.150661, 1.082328, 1.016336, 0.961361, 0.909377, 0.859229, 0.808615, 0.770923, 0.726400,
    0.690444, 0.656639, 0.626977, 0.601282, 0.565640, 0.541549, 0.517613, 0.492610, 0.467891, 0.450866, 0.427611, 0.405115,
    0.393432, 0.367097, 0.355256, 0.344176, 0.327230, 0.308159, 0.299667, 0.288000, 0.278949, 0.263341, 0.250445, 0.243064,
    0.232859, 0.221800, 0.214465, 0.207126, 0.197129, 0.188575, 0.186439, 0.179715, 0.170882, 0.164423, 0.158682, 0.154019,
    0.150277, 0.144272, 0.137206, 0.135669, 0.128931, 0.125974, 0.121932, 0.119566, 0.113806, 0.112343, 0.107325, 0.106804,
    0.104220, 0.098339, 0.095757, 0.095992, 0.092098, 0.089896, 0.087719, 0.086209, 0.085378, 0.082050, 0.078781, 0.078904,
    0.077688, 0.077198, 0.075539, 0.071257, 0.071342, 0.070415, 0.069445, 0.067272, 0.064635, 0.065670, 0.063817, 0.062011,
    0.060226, 0.058354, 0.058203, 0.057580, 0.055477, 0.055900, 0.054993, 0.053452, 0.052578, 0.051701, 0.050199, 0.049844,
    0.047829, 0.047402, 0.046873, 0.044984, 0.045633, 0.044068, 0.043317, 0.042371, 0.042091, 0.042391, 0.041007, 0.039471,
    0.038575, 0.038065, 0.039097, 0.037241, 0.035288, 0.035440, 0.035264, 0.033790, 0.033902, 0.033853, 0.032486, 0.031943,
    0.031140, 0.030939, 0.030288, 0.029462, 0.029969, 0.028097, 0.028703, 0.028627, 0.027646, 0.027564, 0.026505, 0.026356,
    0.025599, 0.026124, 0.025300, 0.024600, 0.024092, 0.024186, 0.023315, 0.023520, 0.023087, 0.022226, 0.022345, 0.021455,
    0.020492, 0.021519, 0.020416, 0.020084, 0.020114, 0.019717, 0.019508, 0.018334, 0.018877, 0.018658, 0.018210, 0.017702,
    0.018155, 0.017342, 0.016846, 0.017155, 0.016306, 0.016247, 0.016075, 0.015419, 0.015835, 0.016264, 0.015473, 0.015212,
    0.014994, 0.014338, 0.015062, 0.014039, 0.014740, 0.014412, 0.014111, 0.013556, 0.013372, 0.012951, 0.012977, 0.013086,
    0.012727, 0.012118, 0.012176, 0.012056, 0.012350, 0.011744, 0.011840, 0.011740, 0.012028, 0.011028, 0.010922, 0.011132,
    0.011492, 0.010939, 0.011151, 0.010446, 0.010507, 0.009898, 0.010788, 0.010553, 0.010305, 0.010137, 0.009395, 0.009406,
    0.009338, 0.009635, 0.009477, 0.009497, 0.008915, 0.009116, 0.009080, 0.009294, 0.008765, 0.008222, 0.008483, 0.008666,
    0.008659, 0.008378, 0.008100, 0.008121, 0.008041, 0.007868, 0.008032, 0.008100, 0.008079, 0.007581, 0.007421, 0.007118,
    0.007416, 0.007828, 0.007152, 0.007088, 0.007122, 0.007133, 0.006934, 0.007103, 0.006555, 0.007113, 0.006829, 0.006883,
    0.006288, 0.006291, 0.006312, 0.006090, 0.006457, 0.006410, 0.006449, 0.005728, 0.006082, 0.006032, 0.005872, 0.005989,
    0.005934, 0.005962, 0.005698, 0.005851, 0.005615, 0.005442, 0.005838, 0.005589, 0.005422, 0.005286, 0.005219, 0.005153,
    0.005222, 0.005520, 0.004971, 0.005049, 0.005134, 0.004912, 0.004870, 0.004990, 0.005039, 0.004782, 0.004498, 0.004962,
    0.004755, 0.004995, 0.004437, 0.004686, 0.004753, 0.004703, 0.004486, 0.004614, 0.004152, 0.004107, 0.004131, 0.004349,
    0.004254, 0.004112, 0.004103, 0.004045, 0.004095, 0.004148, 0.004417, 0.004104, 0.003765, 0.003868, 0.004080, 0.003937,
    0.004091, 0.003892, 0.003849, 0.003681, 0.003834, 0.003530, 0.003934, 0.003552, 0.003714, 0.003282, 0.003898, 0.003909,
    0.003685, 0.003800, 0.003469, 0.003311, 0.003483, 0.003242, 0.003083, 0.003145, 0.003158, 0.003073, 0.003331, 0.003388,
    0.003156, 0.002935, 0.003256, 0.003261, 0.002972, 0.002966, 0.003312, 0.003253, 0.003122, 0.003075, 0.002978, 0.003029,
    0.002883, 0.002793, 0.003039, 0.002905, 0.002511, 0.002779, 0.002964, 0.002755, 0.002747, 0.002863, 0.002777, 0.002514,
    0.002775, 0.002845, 0.002427, 0.002540, 0.002406, 0.002490, 0.002592, 0.002457, 0.002560, 0.002223, 0.002582, 0.002725,
    0.002456, 0.002258, 0.002419, 0.002630, 0.002263, 0.002158, 0.002518, 0.002364, 0.002034, 0.002077, 0.002132, 0.002602,
    0.002309, 0.002357, 0.002333, 0.002194, 0.002213, 0.002270, 0.002202, 0.002308, 0.002042, 0.002417, 0.002220, 0.002359,
    0.002305, 0.002082, 0.002036, 0.002075, 0.002176, 0.002291, 0.002174, 0.002057, 0.001983, 0.001972, 0.001997, 0.001875,
    0.001958, 0.001838, 0.001895, 0.001868, 0.001996, 0.001807, 0.001765, 0.002055, 0.001933, 0.001797, 0.001908, 0.001816,
    0.001790, 0.001774, 0.001731, 0.002141, 0.001511, 0.001603, 0.002065, 0.001799, 0.001736, 0.001954, 0.001836, 0.001729,
    0.001817, 0.001671, 0.001766, 0.001652, 0.001706, 0.001519, 0.001649, 0.001518, 0.001505, 0.001552, 0.001531, 0.001730,
    0.001717, 0.001442, 0.001626, 0.001668, 0.001591, 0.001509, 0.001313, 0.001555, 0.001393, 0.001587, 0.001363, 0.001360,
    0.001785, 0.001410, 0.001363, 0.001278, 0.001427, 0.001407, 0.001102, 0.001358, 0.001408, 0.001509, 0.001399, 0.001411,
    0.001391, 0.001284, 0.001641, 0.001205, 0.001172, 0.001277, 0.001163, 0.001321, 0.001243, 0.001284, 0.001353, 0.001219,
    0.001338, 0.001204, 0.001206, 0.001250, 0.001127, 0.001366, 0.001240, 0.001124 };

  for(Int_t i=0; i<500; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromD0Cent080dataoverLHC16i2abc(){
  //
  // weight function from the ratio of the D0 0-80% over lhc16i2abc
  // This is for Lc/D0 analysis, so the functional form is derived using only 3-16 GeV/c
  //
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]*x+[6])",3.,16.);
  fFuncWeight->SetParameters(-93.179,2.15711,6.45459,1.33679,0.373436,-0.070921,2.52372);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromD0Cent080dataModOhoverLHC16i2abc(){
  //
  // weight function from the ratio of the D0 0-80% over lhc16i2abc
  // Modified taking into account Lc/D0 ratio given in PRC 79, 044905 (2009)
  // This is for Lc/D0 analysis, so the functional form is derived using only 3-16 GeV/c
  //
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","(([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]*x+[6]))*[7]*TMath::Gaus(x,[8],[9])",3.,16.);
  fFuncWeight->SetParameters(-93.179,2.15711,6.45459,1.33679,0.373436,-0.070921,2.52372,1.02369,1.41492,3.09519);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromD0Cent080dataModMartinezoverLHC16i2abc(){
  //
  // weight function from the ratio of the D0 0-80% over lhc16i2abc
  // Modified taking into account Lc/D0 ratio given in PLB 663, 55-60 (2008)
  // This is for Lc/D0 analysis, so the functional form is derived using only 3-16 GeV/c
  //
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","(([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]*x+[6]))*[7]*TMath::Gaus(x,[8],[9])",3.,16.);
  fFuncWeight->SetParameters(-93.179,2.15711,6.45459,1.33679,0.373436,-0.070921,2.52372,0.9,5.0,2.9);
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
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL13overLHC17c3a12(){
	// weight function from the ratio of the LHC17c3a1-LHC17c3a2 MC
	// and FONLL calculations for pp data
	// using D0 simulated pt distribution
	if(fFuncWeight) delete fFuncWeight;
	fFuncWeight=new TF1("funcWeight","(x<10.71)*(([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x))+(x>=10.71)*2.192",0.15,40);
	fFuncWeight->SetParameters(-8.28294e+04,1.04367e+00,1.43041e+02,1.40041e+00,-7.23198e-03,4.71629e+00,-3.15545e+00,2.36495e+00,-7.08928e-03);
	fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5overLHC18a4a2() {
  // Weight function from the ratio of the D0 spectrum in LHC18a4a2 (fast+cent) MC
  // and FONLL calculations for pp data at 5 TeV

  if(fHistoPtWeight) delete fHistoPtWeight;
  fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
  fHistoPtWeight->Sumw2();
  Float_t binc[500]={0.7836, 0.7203, 0.6840, 0.6749, 0.6840, 0.7038, 0.7337, 0.7691, 0.8098, 0.8520, 0.8917, 0.9331, 0.9673, 0.9985, 1.0297, 1.0569, 1.0781, 1.0993, 1.1166, 1.1285, 1.1425, 1.1520, 1.1598, 1.1685, 1.1762, 1.1791, 1.1858, 1.1899, 1.1944, 1.1986, 1.1968, 1.2041, 1.2072, 1.2031, 1.2068, 1.2090, 1.2079, 1.2089, 1.2081, 1.2047, 1.2088, 1.2093, 1.2085, 1.2105, 1.2099, 1.2125, 1.2108, 1.2090, 1.2071, 1.2066, 1.2122, 1.2126, 1.2131, 1.2136, 1.2141, 1.2145, 1.2150, 1.2155, 1.2160, 1.2165, 1.2169, 1.2174, 1.2179, 1.2184, 1.2188, 1.2193, 1.2198, 1.2203, 1.2207, 1.2212, 1.2217, 1.2222, 1.2227, 1.2231, 1.2236, 1.2241, 1.2246, 1.2250, 1.2255, 1.2260, 1.2265, 1.2269, 1.2274, 1.2279, 1.2284, 1.2289, 1.2293, 1.2298, 1.2303, 1.2308, 1.2312, 1.2317, 1.2322, 1.2327, 1.2331, 1.2336, 1.2341, 1.2346, 1.2351, 1.2355, 1.2360, 1.2365, 1.2370, 1.2374, 1.2379, 1.2384, 1.2389, 1.2393, 1.2398, 1.2403, 1.2408, 1.2413, 1.2417, 1.2422, 1.2427, 1.2432, 1.2436, 1.2441, 1.2446, 1.2451, 1.2455, 1.2460, 1.2465, 1.2470, 1.2475, 1.2479, 1.2484, 1.2489, 1.2494, 1.2498, 1.2503, 1.2508, 1.2513, 1.2517, 1.2522, 1.2527, 1.2532, 1.2537, 1.2541, 1.2546, 1.2551, 1.2556, 1.2560, 1.2565, 1.2570, 1.2575, 1.2579, 1.2584, 1.2589, 1.2594, 1.2599, 1.2603, 1.2608, 1.2613, 1.2618, 1.2622, 1.2627, 1.2632, 1.2637, 1.2641, 1.2646, 1.2651, 1.2656, 1.2661, 1.2665, 1.2670, 1.2675, 1.2680, 1.2684, 1.2689, 1.2694, 1.2699, 1.2703, 1.2708, 1.2713, 1.2718, 1.2723, 1.2727, 1.2732, 1.2737, 1.2742, 1.2746, 1.2751, 1.2756, 1.2761, 1.2765, 1.2770, 1.2775, 1.2780, 1.2785, 1.2789, 1.2794, 1.2799, 1.2804, 1.2808, 1.2813, 1.2818, 1.2823, 1.2827, 1.2832, 1.2837, 1.2842, 1.2847, 1.2851, 1.2856, 1.2861, 1.2866, 1.2870, 1.2875, 1.2880, 1.2885, 1.2889, 1.2894, 1.2899, 1.2904, 1.2909, 1.2913, 1.2918, 1.2923, 1.2928, 1.2932, 1.2937, 1.2942, 1.2947, 1.2951, 1.2956, 1.2961, 1.2966, 1.2971, 1.2975, 1.2980, 1.2985, 1.2990, 1.2994, 1.2999, 1.3004, 1.3009, 1.3013, 1.3018, 1.3023, 1.3028, 1.3033, 1.3037, 1.3042, 1.3047, 1.3052, 1.3056, 1.3061, 1.3066, 1.3071, 1.3075, 1.3080, 1.3085, 1.3090, 1.3095, 1.3099, 1.3104, 1.3109, 1.3114, 1.3118, 1.3123, 1.3128, 1.3133, 1.3137, 1.3142, 1.3147, 1.3152, 1.3157, 1.3161, 1.3166, 1.3171, 1.3176, 1.3180, 1.3185, 1.3190, 1.3195, 1.3199, 1.3204, 1.3209, 1.3214, 1.3219, 1.3223, 1.3228, 1.3233, 1.3238, 1.3242, 1.3247, 1.3252, 1.3257, 1.3262, 1.3266, 1.3271, 1.3276, 1.3281, 1.3285, 1.3290, 1.3295, 1.3300, 1.3304, 1.3309, 1.3314, 1.3319, 1.3324, 1.3328, 1.3333, 1.3338, 1.3343, 1.3347, 1.3352, 1.3357, 1.3362, 1.3366, 1.3371, 1.3376, 1.3381, 1.3386, 1.3390, 1.3395, 1.3400, 1.3405, 1.3409, 1.3414, 1.3419, 1.3424, 1.3428, 1.3433, 1.3438, 1.3443, 1.3448, 1.3452, 1.3457, 1.3462, 1.3467, 1.3471, 1.3476, 1.3481, 1.3486, 1.3490, 1.3495, 1.3500, 1.3505, 1.3510, 1.3514, 1.3519, 1.3524, 1.3529, 1.3533, 1.3538, 1.3543, 1.3548, 1.3552, 1.3557, 1.3562, 1.3567, 1.3572, 1.3576, 1.3581, 1.3586, 1.3591, 1.3595, 1.3600, 1.3605, 1.3610, 1.3614, 1.3619, 1.3624, 1.3629, 1.3634, 1.3638, 1.3643, 1.3648, 1.3653, 1.3657, 1.3662, 1.3667, 1.3672, 1.3676, 1.3681, 1.3686, 1.3691, 1.3696, 1.3700, 1.3705, 1.3710, 1.3715, 1.3719, 1.3724, 1.3729, 1.3734, 1.3738, 1.3743, 1.3748, 1.3753, 1.3758, 1.3762, 1.3767, 1.3772, 1.3777, 1.3781, 1.3786, 1.3791, 1.3796, 1.3800, 1.3805, 1.3810, 1.3815, 1.3820, 1.3824, 1.3829, 1.3834, 1.3839, 1.3843, 1.3848, 1.3853, 1.3858, 1.3862, 1.3867, 1.3872, 1.3877, 1.3882, 1.3886, 1.3891, 1.3896, 1.3901, 1.3905, 1.3910, 1.3915, 1.3920, 1.3924, 1.3929, 1.3934, 1.3939, 1.3944, 1.3948, 1.3953, 1.3958, 1.3963, 1.3967, 1.3972, 1.3977, 1.3982, 1.3986, 1.3991, 1.3996, 1.4001, 1.4006, 1.4010, 1.4015, 1.4020, 1.4025, 1.4029, 1.4034, 1.4039, 1.4044, 1.4048, 1.4053, 1.4058, 1.4063, 1.4068, 1.4072, 1.4077, 1.4082, 1.4087, 1.4091, 1.4096, 1.4101, 1.4106, 1.4110, 1.4115, 1.4120, 1.4125, 1.4130, 1.4134, 1.4139, 1.4144, 1.4149, 1.4153, 1.4158, 1.4163, 1.4168, 1.4172, 1.4177, 1.4182, 1.4187, 1.4192, 1.4196, 1.4201, 1.4206, 1.4211, 1.4215, 1.4220, 1.4225, 1.4230, 1.4234, 1.4239, 1.4244, 1.4249, 1.4254, 1.4258, 1.4263};
  for(Int_t i=0; i<500; i++){
    fHistoPtWeight->SetBinContent(i+1,binc[i]);
  }
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5anddataoverLHC20g2a(){
    // weight function from the ratio of the LHC20g2a MC
    // 1-50 GeV/c using data (0-10% centrality) and 0-1 GeV/c using FONLL calculations

    if(fHistoPtWeight) delete fHistoPtWeight;
    fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
    Float_t binc[500]={ 1.275279, 1.178278, 1.097561, 1.037290, 0.995188, 0.969488, 0.958207, 0.958101, 0.967583, 0.989210, 1.026397, 1.073075, 1.117044, 1.154229, 1.185649, 1.211302, 1.231190, 1.245313, 1.253670, 1.256261, 1.253087, 1.244147, 1.229442, 1.208970, 1.182734, 1.150731, 1.112963, 1.069430, 1.032005, 1.003455, 0.973970, 0.942495, 0.899813, 0.848803, 0.801111, 0.756539, 0.714863, 0.675874, 0.639391, 0.605245, 0.573278, 0.543344, 0.515306, 0.489036, 0.464417, 0.441337, 0.419694, 0.399391, 0.380339, 0.362455, 0.345661, 0.329884, 0.315057, 0.301118, 0.288007, 0.275670, 0.264056, 0.253117, 0.242810, 0.233092, 0.223926, 0.215276, 0.207108, 0.199391, 0.192096, 0.185196, 0.178666, 0.172482, 0.166621, 0.161064, 0.155792, 0.150786, 0.146030, 0.141509, 0.137207, 0.133112, 0.129210, 0.125490, 0.121941, 0.118553, 0.115316, 0.112220, 0.109259, 0.106423, 0.103706, 0.101100, 0.098600, 0.096199, 0.093892, 0.091674, 0.089539, 0.087483, 0.085502, 0.083592, 0.081749, 0.079969, 0.078250, 0.076588, 0.074981, 0.073424, 0.071917, 0.070457, 0.069041, 0.067668, 0.066335, 0.065040, 0.063782, 0.062560, 0.061371, 0.060214, 0.059088, 0.057992, 0.056923, 0.055882, 0.054868, 0.053878, 0.052912, 0.051969, 0.051049, 0.050151, 0.049273, 0.048415, 0.047577, 0.046757, 0.045955, 0.045171, 0.044404, 0.043653, 0.042918, 0.042198, 0.041494, 0.040804, 0.040128, 0.039465, 0.038816, 0.038180, 0.037556, 0.036945, 0.036346, 0.035758, 0.035182, 0.034616, 0.034062, 0.033518, 0.032984, 0.032460, 0.031946, 0.031442, 0.030947, 0.030461, 0.029983, 0.029515, 0.029055, 0.028604, 0.028160, 0.027725, 0.027298, 0.026878, 0.026465, 0.026060, 0.025663, 0.025272, 0.024888, 0.024511, 0.024141, 0.023777, 0.023419, 0.023068, 0.022723, 0.022384, 0.022051, 0.021724, 0.021402, 0.021086, 0.020775, 0.020470, 0.020170, 0.019876, 0.019586, 0.019301, 0.019022, 0.018747, 0.018476, 0.018211, 0.017950, 0.017693, 0.017441, 0.017193, 0.016950, 0.016710, 0.016475, 0.016244, 0.016016, 0.015793, 0.015573, 0.015357, 0.015145, 0.014936, 0.014731, 0.014529, 0.014331, 0.014136, 0.013945, 0.013757, 0.013571, 0.013389, 0.013210, 0.013035, 0.012862, 0.012692, 0.012525, 0.012360, 0.012199, 0.012040, 0.011884, 0.011730, 0.011579, 0.011431, 0.011285, 0.011142, 0.011001, 0.010862, 0.010726, 0.010592, 0.010460, 0.010330, 0.010203, 0.010078, 0.009955, 0.009834, 0.009715, 0.009597, 0.009482, 0.009369, 0.009258, 0.009148, 0.009041, 0.008935, 0.008831, 0.008729, 0.008628, 0.008529, 0.008432, 0.008336, 0.008242, 0.008149, 0.008058, 0.007969, 0.007881, 0.007794, 0.007709, 0.007626, 0.007543, 0.007462, 0.007383, 0.007305, 0.007228, 0.007152, 0.007077, 0.007004, 0.006932, 0.006861, 0.006792, 0.006723, 0.006656, 0.006589, 0.006524, 0.006460, 0.006397, 0.006335, 0.006274, 0.006214, 0.006154, 0.006096, 0.006039, 0.005983, 0.005928, 0.005873, 0.005820, 0.005767, 0.005715, 0.005664, 0.005614, 0.005565, 0.005516, 0.005468, 0.005421, 0.005375, 0.005330, 0.005285, 0.005241, 0.005198, 0.005155, 0.005113, 0.005072, 0.005031, 0.004991, 0.004952, 0.004913, 0.004875, 0.004838, 0.004801, 0.004764, 0.004729, 0.004693, 0.004659, 0.004625, 0.004591, 0.004558, 0.004526, 0.004494, 0.004462, 0.004431, 0.004401, 0.004371, 0.004341, 0.004312, 0.004284, 0.004255, 0.004228, 0.004200, 0.004173, 0.004147, 0.004121, 0.004095, 0.004070, 0.004045, 0.004020, 0.003996, 0.003972, 0.003949, 0.003926, 0.003903, 0.003881, 0.003859, 0.003837, 0.003816, 0.003795, 0.003774, 0.003754, 0.003733, 0.003714, 0.003694, 0.003675, 0.003656, 0.003637, 0.003619, 0.003601, 0.003583, 0.003565, 0.003548, 0.003531, 0.003514, 0.003497, 0.003481, 0.003465, 0.003449, 0.003433, 0.003418, 0.003403, 0.003388, 0.003373, 0.003358, 0.003344, 0.003330, 0.003316, 0.003302, 0.003288, 0.003275, 0.003262, 0.003249, 0.003236, 0.003223, 0.003211, 0.003198, 0.003186, 0.003174, 0.003162, 0.003151, 0.003139, 0.003128, 0.003116, 0.003105, 0.003094, 0.003084, 0.003073, 0.003063, 0.003052, 0.003042, 0.003032, 0.003022, 0.003012, 0.003002, 0.002993, 0.002983, 0.002974, 0.002965, 0.002955, 0.002946, 0.002938, 0.002929, 0.002920, 0.002912, 0.002903, 0.002895, 0.002886, 0.002878, 0.002870, 0.002862, 0.002854, 0.002847, 0.002839, 0.002831, 0.002824, 0.002816, 0.002809, 0.002802, 0.002795, 0.002787, 0.002780, 0.002773, 0.002767, 0.002760, 0.002753, 0.002746, 0.002740, 0.002733, 0.002727, 0.002721, 0.002714, 0.002708, 0.002702, 0.002696, 0.002690, 0.002684, 0.002678, 0.002672, 0.002666, 0.002661, 0.002655, 0.002649, 0.002644, 0.002638, 0.002633, 0.002627, 0.002622, 0.002617, 0.002611, 0.002606, 0.002601, 0.002596, 0.002591, 0.002586, 0.002581, 0.002576, 0.002571, 0.002566, 0.002562, 0.002557, 0.002552, 0.002547, 0.002543, 0.002538, 0.002534, 0.002529, 0.002525, 0.002520, 0.002516, 0.002511, 0.002507, 0.002503, 0.002499, 0.002494, 0.002490, 0.002486, 0.002482, 0.002478, 0.002474, 0.002470, 0.002466, 0.002462, 0.002458, 0.002454, 0.002450, 0.002446, 0.002442, 0.002438, 0.002434, 0.002431, 0.002427, 0.002423, 0.002420, 0.002416, 0.002412, 0.002409, 0.002405, 0.002401, 0.002398, 0.002394, 0.002391};
    for(Int_t i=0; i<500; i++){
        fHistoPtWeight->SetBinContent(i+1,binc[i]);
    }
    //SetWeightHistogram();
    fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHF::SetPtWeightsFromFONLL5anddataoverLHC20g2b(){
    // weight function from the ratio of the LHC20g2b MC
    // 1-50 GeV/c using data (30-50% centrality) and 0-1 GeV/c using FONLL calculations

    if(fHistoPtWeight) delete fHistoPtWeight;
    fHistoPtWeight = new TH1F("histoWeight","histoWeight",500,0.,50.);
    Float_t binc[500]={ 1.244760, 1.150068, 1.071265, 1.012415, 0.971299, 0.946193, 0.935163, 0.935049, 0.944369, 0.964778, 0.997505, 1.036911, 1.073968, 1.105850, 1.133293, 1.156296, 1.174859, 1.188983, 1.198668, 1.203912, 1.204718, 1.201084, 1.193010, 1.180497, 1.163544, 1.142152, 1.116320, 1.086049, 1.051338, 1.017632, 0.990162, 0.964726, 0.938994, 0.914297, 0.890133, 0.866480, 0.843481, 0.821118, 0.799372, 0.778227, 0.757666, 0.737672, 0.718229, 0.699323, 0.680939, 0.663061, 0.645675, 0.628769, 0.612328, 0.596340, 0.580792, 0.565672, 0.550967, 0.536667, 0.522759, 0.509234, 0.496080, 0.483287, 0.470845, 0.458744, 0.446975, 0.435528, 0.424395, 0.413567, 0.403035, 0.392791, 0.382827, 0.373135, 0.363708, 0.354539, 0.345619, 0.336943, 0.328503, 0.320293, 0.312306, 0.304537, 0.296978, 0.289626, 0.282473, 0.275513, 0.268743, 0.262156, 0.255748, 0.249513, 0.243446, 0.237544, 0.231801, 0.226213, 0.220775, 0.215484, 0.210335, 0.205325, 0.200450, 0.195705, 0.191088, 0.186594, 0.182221, 0.177965, 0.173823, 0.169791, 0.165866, 0.162047, 0.158329, 0.154709, 0.151186, 0.147757, 0.144418, 0.141168, 0.138003, 0.134923, 0.131923, 0.129003, 0.126159, 0.123390, 0.120694, 0.118068, 0.115512, 0.113022, 0.110597, 0.108235, 0.105934, 0.103694, 0.101511, 0.099385, 0.097314, 0.095297, 0.093331, 0.091416, 0.089551, 0.087733, 0.085961, 0.084235, 0.082553, 0.080914, 0.079316, 0.077759, 0.076241, 0.074762, 0.073320, 0.071914, 0.070543, 0.069207, 0.067904, 0.066634, 0.065395, 0.064187, 0.063008, 0.061859, 0.060738, 0.059645, 0.058578, 0.057538, 0.056523, 0.055532, 0.054565, 0.053622, 0.052702, 0.051803, 0.050927, 0.050071, 0.049235, 0.048419, 0.047623, 0.046845, 0.046086, 0.045344, 0.044620, 0.043912, 0.043221, 0.042546, 0.041887, 0.041242, 0.040612, 0.039997, 0.039395, 0.038807, 0.038233, 0.037671, 0.037122, 0.036584, 0.036059, 0.035545, 0.035043, 0.034552, 0.034071, 0.033600, 0.033140, 0.032690, 0.032249, 0.031818, 0.031395, 0.030982, 0.030577, 0.030181, 0.029792, 0.029412, 0.029040, 0.028675, 0.028318, 0.027967, 0.027624, 0.027288, 0.026958, 0.026635, 0.026318, 0.026007, 0.025703, 0.025404, 0.025111, 0.024823, 0.024541, 0.024265, 0.023993, 0.023726, 0.023465, 0.023208, 0.022956, 0.022708, 0.022465, 0.022226, 0.021992, 0.021761, 0.021535, 0.021312, 0.021094, 0.020879, 0.020667, 0.020460, 0.020255, 0.020055, 0.019857, 0.019663, 0.019472, 0.019284, 0.019099, 0.018917, 0.018738, 0.018561, 0.018388, 0.018217, 0.018049, 0.017883, 0.017720, 0.017559, 0.017401, 0.017245, 0.017091, 0.016939, 0.016790, 0.016643, 0.016498, 0.016355, 0.016214, 0.016075, 0.015937, 0.015802, 0.015669, 0.015537, 0.015407, 0.015279, 0.015153, 0.015028, 0.014904, 0.014783, 0.014663, 0.014544, 0.014427, 0.014311, 0.014197, 0.014084, 0.013973, 0.013863, 0.013754, 0.013647, 0.013540, 0.013435, 0.013332, 0.013229, 0.013128, 0.013027, 0.012928, 0.012830, 0.012733, 0.012637, 0.012543, 0.012449, 0.012356, 0.012264, 0.012173, 0.012083, 0.011995, 0.011907, 0.011819, 0.011733, 0.011648, 0.011563, 0.011480, 0.011397, 0.011315, 0.011234, 0.011154, 0.011074, 0.010995, 0.010917, 0.010840, 0.010763, 0.010687, 0.010612, 0.010538, 0.010464, 0.010391, 0.010318, 0.010247, 0.010175, 0.010105, 0.010035, 0.009966, 0.009897, 0.009829, 0.009762, 0.009695, 0.009629, 0.009563, 0.009498, 0.009433, 0.009369, 0.009305, 0.009242, 0.009180, 0.009118, 0.009057, 0.008996, 0.008935, 0.008875, 0.008816, 0.008757, 0.008698, 0.008640, 0.008583, 0.008525, 0.008469, 0.008413, 0.008357, 0.008301, 0.008246, 0.008192, 0.008138, 0.008084, 0.008031, 0.007978, 0.007925, 0.007873, 0.007822, 0.007770, 0.007719, 0.007669, 0.007619, 0.007569, 0.007520, 0.007471, 0.007422, 0.007373, 0.007325, 0.007278, 0.007231, 0.007184, 0.007137, 0.007091, 0.007045, 0.006999, 0.006954, 0.006909, 0.006864, 0.006820, 0.006776, 0.006732, 0.006689, 0.006646, 0.006603, 0.006560, 0.006518, 0.006476, 0.006435, 0.006393, 0.006352, 0.006311, 0.006271, 0.006231, 0.006191, 0.006151, 0.006112, 0.006073, 0.006034, 0.005995, 0.005957, 0.005919, 0.005881, 0.005843, 0.005806, 0.005769, 0.005732, 0.005696, 0.005659, 0.005623, 0.005587, 0.005552, 0.005517, 0.005481, 0.005447, 0.005412, 0.005378, 0.005343, 0.005309, 0.005276, 0.005242, 0.005209, 0.005176, 0.005143, 0.005110, 0.005078, 0.005046, 0.005014, 0.004982, 0.004950, 0.004919, 0.004888, 0.004857, 0.004826, 0.004795, 0.004765, 0.004735, 0.004705, 0.004675, 0.004646, 0.004616, 0.004587, 0.004558, 0.004529, 0.004501, 0.004472, 0.004444, 0.004416, 0.004388, 0.004360, 0.004333, 0.004305, 0.004278, 0.004251, 0.004224, 0.004197, 0.004171, 0.004145, 0.004118, 0.004093, 0.004067, 0.004041, 0.004016, 0.003990, 0.003965, 0.003940, 0.003915, 0.003891, 0.003866, 0.003842, 0.003817, 0.003793, 0.003769, 0.003746, 0.003722, 0.003699, 0.003675, 0.003652, 0.003629, 0.003606, 0.003584, 0.003561, 0.003539, 0.003516, 0.003494, 0.003472, 0.003450, 0.003429, 0.003407, 0.003386, 0.003364, 0.003343, 0.003322, 0.003301, 0.003280, 0.003260, 0.003239, 0.003219, 0.003199, 0.003178, 0.003158, 0.003139, 0.003119, 0.003099, 0.003080, 0.003060};
    for(Int_t i=0; i<500; i++){
        fHistoPtWeight->SetBinContent(i+1,binc[i]);
    }
    //SetWeightHistogram();
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
  if(fUseMultRatioAsWeight)  weight = pMC; //in this case, fHistoMCNch is already the ratio of data/MC!
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
  Int_t period = -1;   // pp 7 TeV:  0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
                       // pPb: 0-LHC13b, 1-LHC13c
  					   // pp 13 TeV: 0-32, see below

  if (fIsPP13TeVData) {    // setting run numbers for LHC16-17-18 at 13 TeV if pp 13 TeV
  //2016
    if(runNo>=252235 && runNo<=252375)period = 0;//16d
    if(runNo>=252603 && runNo<=253591)period = 1;//16e
    if(runNo>=254124 && runNo<=254332)period = 2;//16g
    if(runNo>=254378 && runNo<=255469)period = 3;//16h
    if(runNo>=256146 && runNo<=256420)period = 4;//16j
    if(runNo>=256504 && runNo<=258537)period = 5;//16k
    if(runNo>=258883 && runNo<=260187)period = 6;//16l
    if(runNo>=262395 && runNo<=264035)period = 7;//16o
    if(runNo>=264076 && runNo<=264347)period = 8;//16p
  //2017
    if(runNo>=270531 && runNo<=270667)period = 9;//17c
    if(runNo>=270822 && runNo<=270830)period = 10;//17e
    if(runNo>=270854 && runNo<=270865)period = 11;//17f
    if(runNo>=271868 && runNo<=273103)period = 12;//17h
    if(runNo>=273591 && runNo<=274442)period = 13;//17i
    if(runNo>=274593 && runNo<=274671)period = 14;//17j
    if(runNo>=274690 && runNo<=276508)period = 15;//17k
    if(runNo>=276551 && runNo<=278216)period = 16;//17l
    if(runNo>=278914 && runNo<=280140)period = 17;//17m
    if(runNo>=280282 && runNo<=281961)period = 18;//17o
    if(runNo>=282504 && runNo<=282704)period = 19;//17r
  //2018
    if(runNo>=284706 && runNo<=285447)period = 20;//18b
    if(runNo>=285978 && runNo<=286350)period = 21;//18d
    if(runNo>=286380 && runNo<=286937)period = 22;//18e
    if(runNo>=287000 && runNo<=287977)period = 23;//18f
    if(runNo>=288619 && runNo<=288750)period = 24;//18g
    if(runNo>=288804 && runNo<=288806)period = 25;//18h
    if(runNo>=288861 && runNo<=288909)period = 26;//18i
    if(runNo>=289165 && runNo<=289201)period = 27;//18k
    if(runNo>=289240 && runNo<=289971)period = 28;//18l
    if(runNo>=290222 && runNo<=292839)period = 29;//18m
    if(runNo>=293357 && runNo<=293359)period = 30;//18n
    if(runNo>=293368 && runNo<=293898)period = 31;//18o
    if(runNo>=294009 && runNo<=294925)period = 32;//18p
    if (period<0 || period>32) return 0;
  } else if (fIsPPbData) {    // setting run numbers for LHC13 if pPb
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

//________________________________________________________________________
Double_t AliCFTaskVertexingHF::ComputeTPCq2(AliAODEvent* aod, AliAODMCHeader* mcHeader, Double_t etamin, Double_t etamax, Double_t ptmin, Double_t ptmax) const {
  /// Compute the q2 for ESE starting from TPC tracks

  Int_t nTracks=aod->GetNumberOfTracks();
  Double_t nHarmonic=2.;
  Double_t q2Vec[2] = {0.,0.};
  Int_t multQvec=0;

  for(Int_t it=0; it<nTracks; it++){
    AliAODTrack* track=(AliAODTrack*)aod->GetTrack(it);
    if(!track) continue;
    TString genname = AliVertexingHFUtils::GetGenerator(track->GetLabel(),mcHeader);
    if(genname.Contains("Hijing")) {
      if(track->TestFilterBit(BIT(8))||track->TestFilterBit(BIT(9))) {
        Double_t pt=track->Pt();
        Double_t eta=track->Eta();
        Double_t phi=track->Phi();
        Double_t qx=TMath::Cos(nHarmonic*phi);
        Double_t qy=TMath::Sin(nHarmonic*phi);
        if(eta<etamax && eta>etamin && pt>ptmin && pt<ptmax) {
          q2Vec[0]+=qx;
          q2Vec[1]+=qy;
          multQvec++;
        }
      }
    }
  }

  Double_t q2 = 0.;
  if(multQvec>0) q2 = TMath::Sqrt(q2Vec[0]*q2Vec[0]+q2Vec[1]*q2Vec[1])/TMath::Sqrt(multQvec);

  return q2;
}
Double_t AliCFTaskVertexingHF::CalculateRTValue(AliAODEvent* esdEvent, AliAODMCHeader *mcHeader, AliCFVertexingHF* cf)
{
   ///! TODO: check MM RT task for MC header version of calc

   /// Calculate RT value for input event (method ported from AliAnalysisTaskUeSpectraRT)
   /// also sets phi of leading particle (fPhiLeading)
   Int_t runNumber = esdEvent->GetRunNumber();
   Int_t eventId = 0;
   Double_t trackRTval = -1;
   if (esdEvent->GetHeader()) eventId = GetEventIdAsLong(esdEvent->GetHeader());
   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   AliESDtrackCuts* esdCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
   trackFilter->AddCuts(esdCutsTPC);


   const Int_t nESDTracks = esdEvent->GetNumberOfTracks();
   TObjArray *fCTSTracks = new TObjArray();
   Double_t nRecTracks = 0;
   AliVParticle* part = 0x0;
   Double_t eta, pt, LeadingPt = -1;
   for (Int_t iT = 0; iT < nESDTracks; iT++)
   {
      part = esdEvent->GetTrack(iT);
      eta = part->Eta();
      pt  = part->Pt();
      if (TMath::Abs(eta) > 1.5) continue; //temporary, hardcoded eta cut to default value
      if (!(TMath::Abs(pt) > 0.15)) continue;

      // Default track filter (to be checked)
      for ( Int_t i = 0; i < 1; i++)
      {
         UInt_t selectDebug = 0;
         if (trackFilter)
         {
            selectDebug = trackFilter->IsSelected(part);
            if (!selectDebug)
            {
               continue;
            }
            /// fill tracks array
            fCTSTracks->Add(part);
            if (!part) continue;
         }
      }
   }

   //find leading object
   TObjArray *LeadingTrackReco = FindLeading(fCTSTracks);
   AliVParticle* LeadingReco = 0;
   TObjArray *regionSortedParticlesReco = 0;
   TObjArray *regionsMinMaxReco = 0;
   TList *listMax = 0;
   TList *listMin = 0;
   if (LeadingTrackReco) {
      LeadingReco = (AliVParticle*)LeadingTrackReco->At(0);
      LeadingPt = LeadingReco->Pt();
      cf->SetPhiLeading(LeadingReco->Phi());
      if (LeadingPt > fMinLeadPtRT && LeadingPt < 300. ) {// calculate only if leading pt is in acceptable range
         //Sorting
         regionSortedParticlesReco = SortRegionsRT((AliVParticle*)LeadingTrackReco->At(0), fCTSTracks);
         // Transverse regions
         regionsMinMaxReco = GetMinMaxRegionRT((TList*)regionSortedParticlesReco->At(2),(TList*)regionSortedParticlesReco->At(3));
         listMax = (TList*)regionsMinMaxReco->At(0);
         listMin = (TList*)regionsMinMaxReco->At(1);

         trackRTval = (listMax->GetEntries() + listMin->GetEntries()) / cf->GetAveMultiInTrans(); //sum of transverse regions / average
      }

   }
  // clean up trackFilter object (else leak)
  if (regionSortedParticlesReco) delete regionSortedParticlesReco;
  if (regionsMinMaxReco) delete regionsMinMaxReco;
  if (LeadingTrackReco) delete LeadingTrackReco;
  if(trackFilter) delete trackFilter;
  return trackRTval;
}


ULong64_t AliCFTaskVertexingHF::GetEventIdAsLong(AliVHeader* header)
{
//	unique ID for each event
   return ((ULong64_t)header->GetBunchCrossNumber() +
           (ULong64_t)header->GetOrbitNumber()*3564 +
           (ULong64_t)header->GetPeriodNumber()*16777215*3564);

}
TObjArray *AliCFTaskVertexingHF::FindLeading(TObjArray *array)
{
   if (!array) return 0;
   Int_t nTracks = array->GetEntriesFast();
   if (!nTracks) return 0;
   AliVParticle *part = 0x0;
   TObjArray *tracks = new TObjArray(nTracks);
   for (Int_t ipart = 0; ipart < nTracks; ipart++) {
      part = (AliVParticle*)(array->At(ipart));
      if(!part) continue;
      tracks->AddLast(part);
   }
   QSortTracks(*tracks, 0, tracks->GetEntriesFast());
   nTracks = tracks->GetEntriesFast();
   if (!nTracks) return 0;
   return tracks;
}

void AliCFTaskVertexingHF::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
   //Sort array by pT
  static TObject *tmp;
  static int i;           // "static" to save stack space
  int j;

  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;

      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;
    } else {
      QSortTracks(a, j + 1, last);
      last = j;
    }
  }
}



TObjArray *AliCFTaskVertexingHF::SortRegionsRT(const AliVParticle* leading, TObjArray *array)
{
   if (!array) return 0;
   static const Double_t k60rad = 60.*TMath::Pi()/180.;
   static const Double_t k120rad = 120.*TMath::Pi()/180.;

   // define output lists of particles
   TList *toward = new TList();
   TList *away = new TList();
   TList *transverse1 = new TList();
   TList *transverse2 = new TList();
   TObjArray *regionParticles = new TObjArray;
   regionParticles->SetOwner();

   regionParticles->AddLast(toward);
   regionParticles->AddLast(away);
   regionParticles->AddLast(transverse1);
   regionParticles->AddLast(transverse2);
   if (!leading) return regionParticles;

   TVector3 leadVect(leading->Px(),leading->Py(),leading->Pz());
   Int_t nTracks = array->GetEntriesFast();
   if (!nTracks) return 0;

   //loop over tracks
   AliVParticle* part = 0x0;
   for (Int_t ipart = 0; ipart < nTracks; ipart++) {
      part = (AliVParticle*)(array->At(ipart));
      if(!part) continue;
      //vector notation for particles
      TVector3 partVect(part->Px(), part->Py(), part->Pz());
      Int_t region = 0;
      Float_t deltaPhi = leadVect.DeltaPhi(partVect);
      if (deltaPhi <= -TMath::PiOver2()) deltaPhi+= TMath::TwoPi();
      if (deltaPhi > 3*TMath::PiOver2()) deltaPhi-= TMath::TwoPi();
      Double_t fUeDeltaPhiMinCut = TMath::DegToRad()*60.;
      Double_t fUeDeltaPhiMaxCut = TMath::DegToRad()*120.;

      //transverse regions

      if((deltaPhi<-fUeDeltaPhiMinCut) || (deltaPhi >2*fUeDeltaPhiMaxCut))region = -1; //left
      if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi <fUeDeltaPhiMaxCut)) region = 1;   //right

      if(deltaPhi > -fUeDeltaPhiMinCut && deltaPhi < fUeDeltaPhiMinCut) region = 2;    //forward
      if(deltaPhi > fUeDeltaPhiMaxCut && deltaPhi < 2*fUeDeltaPhiMaxCut) region = -2;  //backward

      // skip leading particle
      if(leading == part) continue;
      if(part->Pt() >= leading->Pt()) continue;
      if(!region)continue;

      if(region == 1) transverse1->Add(part);
      if(region == -1) transverse2->Add(part);
      if(region == 2) toward->Add(part);
      if(region == -2) away->Add(part);
   }//end loop on tracks

   return regionParticles;
}


TObjArray* AliCFTaskVertexingHF::GetMinMaxRegionRT(TList *transv1, TList *transv2)
{
// Returns two lists of particles, one for MIN and one for MAX region
  Double_t sumpT1 = 0.;
  Double_t sumpT2 = 0.;

  Int_t particles1 = transv1->GetEntries();
  Int_t particles2 = transv2->GetEntries();
  AliVParticle *part = 0x0;
// Loop on transverse region 1
  for(Int_t i=0; i<particles1; i++){
   part = (AliVParticle*)transv1->At(i);
   sumpT1 +=  part->Pt();
   }

// Loop on transverse region 2
  for(Int_t i=0; i<particles2; i++){
   part = (AliVParticle*)transv2->At(i);
   sumpT2 +=  part->Pt();
   }

  TObjArray *regionParticles = new TObjArray;
  if(sumpT2 >= sumpT1){
   regionParticles->AddLast(transv1); // MIN
   regionParticles->AddLast(transv2); // MAX
   }
  else{
   regionParticles->AddLast(transv2); // MIN
   regionParticles->AddLast(transv1); // MAX
   }

  return regionParticles;
}
