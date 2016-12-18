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
//-----------------------------------------------------------------------
// Modified by for analysis of D0 candidates using the cut-variation
// feed-down subtraction method
//  by
// Andrea Rossi <Andrea.Rossi@cern.ch>
// Felix Reidt  <Felix.Reidt@cern.ch>
//-----------------------------------------------------------------------
#include "AliCFTaskVertexingHFCutVarFDSub.h"

#include "TCanvas.h"
#include "TProfile.h"
#include "TH1I.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"

#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliLog.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"

#include "THnSparse.h"
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
#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliHFsubtractBFDcuts.h"

/// \cond CLASSIMP
ClassImp(AliCFTaskVertexingHFCutVarFDSub);
/// \endcond

//__________________________________________________________________________
AliCFTaskVertexingHFCutVarFDSub::AliCFTaskVertexingHFCutVarFDSub() :
  AliAnalysisTaskSE(),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
  fListProfiles(0),
  fCountMC(0),
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
  fHistoMotherPtWeight(0x0),
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
  fCutOnMomConservation(0.00001),
  fObjSpr(0x0),
  fTHnAnalysis(0x0),
  fTHnGenerator(0x0),
  fhBptCutVar(0x0),
  fhBquarkPt(0x0),
  fhCquarkPt(0x0),
  fListBdecays(0x0),
  fQAHists(0x0)
{
  //
  //Default ctor
  //
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}
//___________________________________________________________________________
AliCFTaskVertexingHFCutVarFDSub::AliCFTaskVertexingHFCutVarFDSub(const Char_t* name, AliRDHFCuts* cuts, TF1* func) :
  AliAnalysisTaskSE(name),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
  fListProfiles(0),
  fCountMC(0),
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
  fHistoMotherPtWeight(0x0),
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
  fCutOnMomConservation(0.00001),
  fObjSpr(0x0),
  fTHnAnalysis(0x0),
  fTHnGenerator(0x0),
  fhBptCutVar(0x0),
  fhBquarkPt(0x0),
  fhCquarkPt(0x0),
  fListBdecays(0x0),
  fQAHists(0x0)
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
  DefineOutput(6,THnSparseF::Class());
  DefineOutput(7,TH3F::Class());
  DefineOutput(8,TH1F::Class());
  DefineOutput(9,TH1F::Class());
  DefineOutput(10,TH1F::Class());
  DefineOutput(11,TList::Class());
  DefineOutput(12,TList::Class());

  fCuts->PrintAll();
}

//___________________________________________________________________________
AliCFTaskVertexingHFCutVarFDSub& AliCFTaskVertexingHFCutVarFDSub::operator=(const AliCFTaskVertexingHFCutVarFDSub& c)
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
    fHistoMotherPtWeight = c.fHistoMotherPtWeight;
    fHistoMeasNch = c.fHistoMeasNch;
    fHistoMCNch = c.fHistoMCNch;
    for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=c.fMultEstimatorAvg[i];
    fObjSpr=c.fObjSpr;
    fTHnAnalysis=c.fTHnAnalysis;
    fTHnGenerator=c.fTHnGenerator;
    fhBptCutVar=c.fhBptCutVar;
    fhBquarkPt=c.fhBquarkPt;
    fhCquarkPt=c.fhCquarkPt;
    fListBdecays=c.fListBdecays;  // FIXME: TList copy contructor not implemented
    fQAHists=c.fQAHists;  // FIXME: TList copy contructor not implemented
  }
  return *this;
}

//___________________________________________________________________________
AliCFTaskVertexingHFCutVarFDSub::AliCFTaskVertexingHFCutVarFDSub(const AliCFTaskVertexingHFCutVarFDSub& c) :
  AliAnalysisTaskSE(c),
  fCFManager(c.fCFManager),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fCorrelation(c.fCorrelation),
  fListProfiles(c.fListProfiles),
  fCountMC(c.fCountMC),
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
  fHistoMotherPtWeight(c.fHistoMotherPtWeight),
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
  fCutOnMomConservation(c.fCutOnMomConservation),
  fObjSpr(c.fObjSpr),
  fTHnAnalysis(c.fTHnAnalysis),
  fTHnGenerator(c.fTHnGenerator),
  fhBptCutVar(c.fhBptCutVar),
  fhBquarkPt(c.fhBquarkPt),
  fhCquarkPt(c.fhCquarkPt),
  fListBdecays(c.fListBdecays), // FIXME: TList copy contructor not implemented
  fQAHists(c.fQAHists) // FIXME: TList copy contructor not implemented
{
  //
  // Copy Constructor
  //
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=c.fMultEstimatorAvg[i];
}

//___________________________________________________________________________
AliCFTaskVertexingHFCutVarFDSub::~AliCFTaskVertexingHFCutVarFDSub()
{
  //
  //destructor
  //
  if (fCFManager)           { delete fCFManager;           fCFManager=0x0;           }
  if (fHistEventsProcessed) { delete fHistEventsProcessed; fHistEventsProcessed=0x0; }
  if (fCorrelation)         { delete fCorrelation;         fCorrelation=0x0;         }
  if (fListProfiles)        { delete fListProfiles;        fListProfiles=0x0;        }
  if (fCuts)                { delete fCuts;                fCuts=0x0;                }
  if (fFuncWeight)          { delete fFuncWeight;          fFuncWeight=0x0;          }
  if (fHistoPtWeight)       { delete fHistoPtWeight;       fHistoPtWeight=0x0;       }
  if (fHistoMotherPtWeight) { delete fHistoMotherPtWeight; fHistoMotherPtWeight=0x0; }
  if (fHistoMeasNch)        { delete fHistoMeasNch;        fHistoMeasNch=0x0;        }
  if (fHistoMCNch)          { delete fHistoMCNch;          fHistoMCNch=0x0;          }
  if (fObjSpr)              { delete fObjSpr;              fObjSpr=0x0;              }
  if (fTHnAnalysis)         { delete fTHnAnalysis;         fTHnAnalysis=0x0;         }
  if (fTHnGenerator)        { delete fTHnGenerator;        fTHnGenerator=0x0;        }
  if (fhBptCutVar)          { delete fhBptCutVar;          fhBptCutVar=0x0;          }
  if (fhBquarkPt)           { delete fhBquarkPt;           fhBquarkPt=0x0;           }
  if (fhCquarkPt)           { delete fhCquarkPt;           fhCquarkPt=0x0;           }
  if (fListBdecays)         { delete fListBdecays;         fListBdecays=0x0;         }
  if (fQAHists)             { delete fQAHists;             fQAHists=0x0;             }
  for(Int_t i=0; i<4; i++) {
    if(fMultEstimatorAvg[i]) { delete fMultEstimatorAvg[i]; fMultEstimatorAvg[i]=0x0; }
  }
}

//__________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::Init()
{
  //
  // Initialization
  //

  if (fDebug>1) printf("AliCFTaskVertexingHFCutVarFDSub::Init()");
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

  /// Save also the weight functions or histograms
  if(fFuncWeight)          fListProfiles->Add(fFuncWeight);
  if(fHistoPtWeight)       fListProfiles->Add(fHistoPtWeight);
  if(fHistoMotherPtWeight) fListProfiles->Add(fHistoMotherPtWeight);
  if(fHistoMeasNch)        fListProfiles->Add(fHistoMeasNch);
  if(fHistoMCNch)          fListProfiles->Add(fHistoMCNch);

  PostData(5,fListProfiles);

  return;
}

//_________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::UserExec(Option_t *)
{
  //
  // Main loop function
  //

  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fCorrelation);

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

    // Obtain B0 pt spectrum
    if ((((mcPart->GetPdgCode()%1000)/100==5||(mcPart->GetPdgCode()%1000)/100==-5)) &&  //||(mcPart->GetPdgCode()%10000)/1000==-5)||(mcPart->GetPdgCode()%10000)/1000==-5)) &&
        (mcPart->Y()<1.)&&(mcPart->Y()>-1.)) {
      fhBptCutVar->Fill(mcPart->Pt());
    }

    // Obtain the pt spectrum of the injected quarks
    if (TMath::Abs(mcPart->GetPdgCode())==5&&TMath::Abs(mcPart->Y())<1.) fhBquarkPt->Fill(mcPart->Pt());
    if (TMath::Abs(mcPart->GetPdgCode())==4&&TMath::Abs(mcPart->Y())<1.) fhCquarkPt->Fill(mcPart->Pt());

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
        AliDebug(3,"MC Lim Acc container filled\n");
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
        if(fTHnGenerator){fObjSpr->FillGenStep(mcPart, mcPart->Pt(), fWeight, mcArray, mcHeader);}

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
  fCountAcc += icountAcc;
  fCountVertex+= icountVertex;
  fCountRefit+= icountRefit;

  AliDebug(2,Form("Found %d vertices for decay channel %d",arrayBranch->GetEntriesFast(),fDecayChannel));

  for(Int_t iCandid = 0; iCandid<arrayBranch->GetEntriesFast();iCandid++){
    AliAODRecoDecayHF* charmCandidate=0x0;
    switch (fDecayChannel){
    case 2:{
      charmCandidate = (AliAODRecoDecayHF2Prong*)arrayBranch->At(iCandid);
      break;
    }
    case 21:
    case 22:{
      charmCandidate = (AliAODRecoCascadeHF*)arrayBranch->At(iCandid);
      break;
    }
    case 31:
    case 32:
    case 33:{
      charmCandidate = (AliAODRecoDecayHF3Prong*)arrayBranch->At(iCandid);
      break;
    }
    case 4:{
      charmCandidate = (AliAODRecoDecayHF4Prong*)arrayBranch->At(iCandid);
      break;
    }
    default:
      break;
    }

    Double_t motherPtWeight=1.;
    if (fOriginDselection>0 && fUseMotherPtWeight) {
      if (fHistoMotherPtWeight) { /// using an histogram as weight function
        AliDebug(2,"Using Histogram as mother pt weight function");
        Double_t motherPt = GetMotherPtFromRecoDecay(charmCandidate,mcArray);
        motherPtWeight = GetMotherPtWeightFromHistogram(motherPt);
        AliDebug(2,Form("mother pt = %f, weight = %f", motherPt, fWeight));
      }
      else AliDebug(2, "Couldn't find mother pt dependend weights");
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
          Int_t recoITSnCluster = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kTracks);
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

                AliAODRecoDecayHF2Prong *d0toKpi = (AliAODRecoDecayHF2Prong*)charmCandidate;
                fObjSpr->SetFillMC(kTRUE);
                fObjSpr->FillSparses(d0toKpi, recoAnalysisCuts, d0toKpi->Pt(),d0toKpi->InvMassD0(), d0toKpi->InvMassD0bar(), fWeight*weigPID, mcArray, aodEvent, mcHeader);

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
void AliCFTaskVertexingHFCutVarFDSub::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliAnalysisTaskSE::Terminate();

  AliInfo(Form("Found %i MC particles that are %s in MC, in %d events",fCountMC,fPartName.Data(),fEvents));
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
void AliCFTaskVertexingHFCutVarFDSub::UserCreateOutputObjects()
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

  if(!fObjSpr){
    fObjSpr=new AliHFsubtractBFDcuts("feedDownCuts","feedDownCuts");
    fObjSpr->InitHistos();
    fTHnAnalysis=fObjSpr->GetSparseMC();
    fTHnGenerator=fObjSpr->GetSparseMCgen();
    fListBdecays=fObjSpr->GetDecayStrings();
    fQAHists=fObjSpr->GetQAhists();
  }

  fhBptCutVar = new TH1F("hBptCutVar", "B meson #it{p}_{T} spectrum;#it{p}_{T};Counts (a.u.)",201,0.,50.25);
  fhBquarkPt = new TH1F("hBquarkPt", "B quark #it{p}_{T} spectrum;#it{p}_{T};Counts (a.u.)",201,0.,50.25);
  fhCquarkPt = new TH1F("hCquarkPt", "C quark #it{p}_{T} spectrum;#it{p}_{T};Counts (a.u.)",201,0.,50.25);

  PostData( 1,fHistEventsProcessed) ;
  PostData( 2,fCFManager->GetParticleContainer()) ;
  PostData( 3,fCorrelation) ;
  PostData( 6,fTHnAnalysis);
  PostData( 7,fTHnGenerator);
  PostData( 8,fhBptCutVar);
  PostData( 9,fhBquarkPt);
  PostData(10,fhCquarkPt);
  PostData(11,fListBdecays);
  PostData(12,fQAHists);
}


//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL276overLHC12a17a(){
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
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromDataPbPb276overLHC12a17a(){
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
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL276overLHC12a17b(){
  // weight function from the ratio of the LHC12a17b MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,50.);
  fFuncWeight->SetParameters(1.92381e+01, 5.05055e+00, 1.05314e+01, 2.5, 1.88214e-03, 3.44871e+00, -9.74325e-02, 1.97671e+00, -3.21278e-01);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL276andBAMPSoverLHC12a17b(){
  // weight function from the ratio of the LHC12a17b MC
  // and FONLL calculations for pp data
  // corrected by the BAMPS Raa calculation for 30-50% CC
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,50.);
  fFuncWeight->SetParameters(6.10443e+00, 1.53487e+00, 1.99474e+00, 2.5, 5.51172e-03, 5.86590e+00, -5.46963e-01, 9.41201e-02, -1.64323e-01);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL5overLHC13d3(){
  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
  fFuncWeight->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetMotherPtWeightsFromFONLL5overLHC13d3(){
  /// weight function from the ratio of the LHC13d3 MC
  /// and FONLL calculations for pp data
  if(fHistoMotherPtWeight) delete fHistoMotherPtWeight;

  const Int_t nBins=48;

  Double_t weights[nBins] = { 0.958796, 0.94262, 1.00746, 1.00264, 1.0103, 1.00033, 0.970773, 0.925759, 0.873651, 0.827537, 0.793838, 0.770711, 0.753224, 0.740698, 0.729354, 0.720129, 0.71487, 0.711469, 0.706371, 0.704551, 0.705003, 0.702572, 0.704875, 0.706797, 0.70574, 0.710273, 0.714451, 0.716691, 0.720462, 0.719736, 0.722672, 0.724752, 0.731169, 0.734225, 0.741758, 0.749343, 0.741116, 0.747457, 0.746682, 0.754226, 0.763797, 0.765175, 0.762687, 0.772638, 0.767542, 0.779808, 0.782834, 0.781626 };

  Double_t xMin=0.;
  Double_t xMax=24.;
  fHistoMotherPtWeight = new TH1F("hMotherPtWeight",";Mother #it{p}_T;Weight",nBins,xMin,xMax);
  for (Int_t ix=1;ix<nBins+1;++ix) fHistoMotherPtWeight->SetBinContent(ix, weights[ix-1]);
  fUseMotherPtWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL7overLHC10f6a(){
  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
  fFuncWeight->SetParameters(2.41522e+01,4.92146e+00,6.72495e+00,2.5,6.15361e-03,4.78995e+00,-4.29135e-01,3.99421e-01,-1.57220e-02);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL7overLHC12a12(){
  // weight function from the ratio of the LHC12a12 MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,50.);
  fFuncWeight->SetParameters(1.31497e+01,3.30503e+00,3.45594e+00,2.5,2.28642e-02,1.42372e+00,2.32892e-04,5.21986e-02,-2.14236e-01,3.86200e+00);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL7overLHC12a12bis(){
  // weight function from the ratio of the LHC12a12bis MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,50.);
  fFuncWeight->SetParameters(6.54519e+00,2.74007e+00,2.48325e+00,2.5,1.61113e-01,-5.32546e-01,-3.75916e-04,2.38189e-01,-2.17561e-01,2.35975e+00);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL7overLHC13e2fix(){
  // weight function from the ratio of the LHC13e2fix MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,50.);
  fFuncWeight->SetParameters(1.85862e+01,2.48171e+00,3.39356e+00,2.5,1.70426e-02,2.28453e+00,-4.57749e-02,5.84585e-02,-3.19719e-01,4.16789e+00);
  fUseWeight=kTRUE;
}

//________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL5overLHC10f6a(){
  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
  fFuncWeight->SetParameters(2.77730e+01,4.78942e+00,7.45378e+00,2.5,9.86255e-02,2.30120e+00,-4.16435e-01,3.43770e-01,-2.29380e-02);
  fUseWeight=kTRUE;
}

//________________________________________________________________
void AliCFTaskVertexingHFCutVarFDSub::SetPtWeightsFromFONLL276overLHC10f6a(){
  // weight function from the ratio of the LHC10f6a MC
  // and FONLL calculations for pp data
  if(fFuncWeight) delete fFuncWeight;
  fFuncWeight=new TF1("funcWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x+[9])",0.15,40.);
  fFuncWeight->SetParameters(1.34412e+01,3.20068e+00,5.14481e+00,2.5,7.59405e-04,7.51821e+00,-3.93811e-01,2.16849e-02,-3.37768e-02,2.40308e+00);
  fUseWeight=kTRUE;
}

//_________________________________________________________________________
Double_t AliCFTaskVertexingHFCutVarFDSub::GetWeight(Float_t pt)
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
Double_t AliCFTaskVertexingHFCutVarFDSub::dNdptFit(Float_t pt, Double_t* par)
{
  //
  // calculating dNdpt
  //

  Double_t denom =  TMath::Power((pt/par[1]), par[3] );
  Double_t dNdpt = par[0]*pt/TMath::Power(1.+denom, par[2]);

  return dNdpt;
}

//_________________________________________________________________________
Double_t AliCFTaskVertexingHFCutVarFDSub::GetPtWeightFromHistogram(Float_t pt)
{
  //
  // Using an histogram as weight function
  //  weight = 0 in the range outside the function
  //
  Double_t weight = 0.0;
  Int_t histoNbins = fHistoPtWeight->GetNbinsX();
  Int_t histobin = fHistoPtWeight->FindBin(pt);
  if( (histobin>0) && (histobin<=histoNbins) ) {
    weight = fHistoPtWeight->GetBinContent(histobin);
  }

  return weight;
}

//_________________________________________________________________________
Double_t AliCFTaskVertexingHFCutVarFDSub::GetMotherPtWeightFromHistogram(Float_t motherPt)
{
  //
  /// Using an histogram as weight function
  ///  weight = 0 in the range outside the function
  //
  Double_t weight = 0.0;
  Int_t histoNbins = fHistoMotherPtWeight->GetNbinsX();
  Int_t histobin = fHistoMotherPtWeight->FindBin(motherPt);
  if( (histobin>0) && (histobin<=histoNbins) ) {
    weight = fHistoMotherPtWeight->GetBinContent(histobin);
  }

  return weight;
}

//_________________________________________________________________________
Double_t AliCFTaskVertexingHFCutVarFDSub::GetMotherPtFromRecoDecay(AliAODRecoDecayHF* charmCandidate, TClonesArray* mcArray)
{
  //
  /// Get mother (B hadron) pt from charm candidate
  //
  Double_t motherPt = -1.;
  Int_t labDau0=((AliAODTrack*)charmCandidate->GetDaughter(0))->GetLabel();
  if (labDau0<0) labDau0*=-1;
  AliAODMCParticle* firstDau=(AliAODMCParticle*)mcArray->UncheckedAt(labDau0);
  Int_t labMother=firstDau->GetMother();
  if (labMother<0) labMother*=-1;
  AliAODMCParticle* mother=(AliAODMCParticle*)mcArray->UncheckedAt(labMother);
  Int_t pdgMother=mother->GetPdgCode();
  while ((pdgMother%1000)/100==4 || (pdgMother%10000)/1000==4) {
    /// follow chained decay up to a B hadron
    labMother=mother->GetMother();
    if (labMother<0) labMother*=-1;
    mother=(AliAODMCParticle*)mcArray->UncheckedAt(labMother);
    pdgMother=mother->GetPdgCode();
  }
  if ((pdgMother%1000)/100!=5 && (pdgMother%10000)/1000!=5) {
    AliDebug(3, "Found strange decay, expected the mother to be a beauty hadron!");
  }
  else {
    motherPt=mother->Pt();
  }
  return motherPt;
}


//__________________________________________________________________________________________________
Double_t AliCFTaskVertexingHFCutVarFDSub::GetZWeight(Float_t z, Int_t runnumber){
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
Double_t AliCFTaskVertexingHFCutVarFDSub::DodzFit(Float_t z, Double_t* par) {

  //
  // Gaussian z-vtx shape
  //
  //gaussian = [0]/TMath::Sqrt(2.*TMath::Pi())/[2]*exp[-(x-[1])*(x-[1])/(2*[2]*[2])]

  Double_t value =  par[0]/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(z-par[1])*(z-par[1])/2./par[2]/par[2]);

  return value;
}
//__________________________________________________________________________________________________
Double_t AliCFTaskVertexingHFCutVarFDSub::GetNchWeight(Int_t nch){
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
void AliCFTaskVertexingHFCutVarFDSub::CreateMeasuredNchHisto(){
  // creates historgam with measured multiplcity distribution in pp 7 TeV collisions (from Eur. Phys. J. C (2010) 68: 345–354)
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
Bool_t AliCFTaskVertexingHFCutVarFDSub::ProcessDs(Int_t recoAnalysisCuts) const{
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
Bool_t AliCFTaskVertexingHFCutVarFDSub::ProcessLctoV0Bachelor(Int_t recoAnalysisCuts) const{
  // processes output of Lc->V0+bnachelor is selected
  Bool_t keep=kFALSE;
  if(recoAnalysisCuts > 0){

    Int_t isK0Sp = recoAnalysisCuts&1;
    Int_t isLambdaBarpi = recoAnalysisCuts&2;
    Int_t isLambdapi = recoAnalysisCuts&4;

    if(fLctoV0bachelorOption==1){
      if(isK0Sp) keep=kTRUE;
    }
    else if(fLctoV0bachelorOption==2){
      if(isLambdaBarpi) keep=kTRUE;
    }
    else if(fLctoV0bachelorOption==4){
      if(isLambdapi) keep=kTRUE;
    }
    else if(fLctoV0bachelorOption==7) {
      if (isK0Sp || isLambdaBarpi || isLambdapi) keep=kTRUE;
    }
  }
  return keep;
}

//____________________________________________________________________________
TProfile* AliCFTaskVertexingHFCutVarFDSub::GetEstimatorHistogram(const AliVEvent* event){
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
