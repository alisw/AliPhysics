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
#include <TH1I.h>
#include <TStyle.h>
#include <TFile.h>
#include <TF1.h>

#include "AliCFTaskVertexingHF.h"
#include "AliStack.h"
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
#include "AliPIDResponse.h"

//__________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF() :
  AliAnalysisTaskSE(),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
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
  fHistoMeasNch(0x0),
  fHistoMCNch(0x0),
  fResonantDecay(0),
  fLctoV0bachelorOption(1),
  fGenLctoV0bachelorOption(0),
  fUseSelectionBit(kTRUE),
  fPDGcode(0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF(const Char_t* name, AliRDHFCuts* cuts, TF1* func) :
  AliAnalysisTaskSE(name),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0),
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
  fHistoMeasNch(0x0),
  fHistoMCNch(0x0),
  fResonantDecay(0),
  fLctoV0bachelorOption(1),
  fGenLctoV0bachelorOption(0),
  fUseSelectionBit(kTRUE),
  fPDGcode(0)
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
    fHistoMeasNch = c.fHistoMeasNch;
    fHistoMCNch = c.fHistoMCNch;
  }
  return *this;
}

//___________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF(const AliCFTaskVertexingHF& c) :
  AliAnalysisTaskSE(c),
  fCFManager(c.fCFManager),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fCorrelation(c.fCorrelation),
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
  fHistoMeasNch(c.fHistoMeasNch),
  fHistoMCNch(c.fHistoMCNch),
  fResonantDecay(c.fResonantDecay),
  fLctoV0bachelorOption(c.fLctoV0bachelorOption),
  fGenLctoV0bachelorOption(c.fGenLctoV0bachelorOption),
  fUseSelectionBit(c.fUseSelectionBit),
  fPDGcode(c.fPDGcode)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFTaskVertexingHF::~AliCFTaskVertexingHF() 
{
  //
  //destructor
  //
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fCorrelation)	  delete fCorrelation ;
  if (fCuts)                delete fCuts;
  if (fFuncWeight)          delete fFuncWeight;
  if (fHistoMeasNch)        delete fHistoMeasNch;
  if (fHistoMCNch)          delete fHistoMCNch;
}

//_________________________________________________________________________-
void AliCFTaskVertexingHF::Init()
{
  //
  // Initialization
  //
	
  if (fDebug>1) printf("AliCFTaskVertexingHF::Init()");
  if(fUseWeight && fUseZWeight) { AliFatal("Can not use at the same time pt and z-vtx weights, please choose"); return; }
  if(fUseWeight && fUseNchWeight) { AliFatal("Can not use at the same time pt and Nch weights, please choose"); return; }
  if(fUseNchWeight && !fHistoMCNch) { AliFatal("Need to pass the MC Nch distribution to use Nch weights"); return; }
  if(fUseNchWeight) CreateMeasuredNchHisto();

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
  if (!aodVtx) return;
	
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
    break;
  }
  case 22:{
    cfVtxHF = new AliCFVertexingHFLctoV0bachelor(mcArray, fOriginDselection,fGenLctoV0bachelorOption); // Lc -> K0S+proton
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
  fWeight=1.;
  if(fUseZWeight) fWeight *= GetZWeight(zMCVertex,runnumber);
  if(fUseNchWeight){
    Int_t nChargedMCPhysicalPrimary=AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(mcArray,-1.0,1.0);
    fWeight *= GetNchWeight(nChargedMCPhysicalPrimary);
    AliDebug(2,Form("Using Nch weights, Mult=%d Weight=%f\n",nChargedMCPhysicalPrimary,fWeight));	
  }

  if (TMath::Abs(zMCVertex) > fCuts->GetMaxVtxZ()){
    AliDebug(3,Form("z coordinate of MC vertex = %f, it was required to be within [-%f, +%f], skipping event", zMCVertex, fCuts->GetMaxVtxZ(), fCuts->GetMaxVtxZ()));
    delete[] containerInput;
    delete[] containerInputMC;
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
      return;
    }    
  }  else { // keep all centralities
    fCuts->SetMinCentrality(0.);
    fCuts->SetMaxCentrality(100.);
  }
	
  Float_t centValue = fCuts->GetCentrality(aodEvent);
  cfVtxHF->SetCentralityValue(centValue);  
	
  // number of tracklets - multiplicity estimator
  Double_t multiplicity = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.)); // casted to double because the CF is filled with doubles
  cfVtxHF->SetMultiplicity(multiplicity);
	
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
    cfVtxHF->SetMCCandidateParam(iPart);
	 
	  
    if (!(cfVtxHF->SetLabelArray())){
      AliDebug(2,Form("Impossible to set the label array (decaychannel = %d)",fDecayChannel));
      continue;
    }		   

    //check the candiate family at MC level
    if (!(cfVtxHF->CheckMCPartFamily(mcPart, mcArray))) {
      AliDebug(2,Form("Check on the family wrong!!! (decaychannel = %d)",fDecayChannel));
      continue;
    }
    else{
      AliInfo(Form("Check on the family OK!!! (decaychannel = %d)",fDecayChannel));
    }
		
    //Fill the MC container
    Bool_t mcContainerFilled = cfVtxHF -> FillMCContainer(containerInputMC);
    AliDebug(2,Form("mcContainerFilled = %d)",mcContainerFilled));
    if (mcContainerFilled) {
      if (fUseWeight){
	if (fFuncWeight){ // using user-defined function
	  AliDebug(2,"Using function");
	  fWeight = fFuncWeight->Eval(containerInputMC[0]);				     
	}
	else{ // using FONLL
	  AliDebug(2,"Using FONLL");
	  fWeight = GetWeight(containerInputMC[0]);
	}
	AliDebug(2,Form("pt = %f, weight = %f",containerInputMC[0], fWeight));
      }
      if (!fCuts->IsInFiducialAcceptance(containerInputMC[0],containerInputMC[1])) continue;
      //MC Limited Acceptance
      if (TMath::Abs(containerInputMC[1]) < 0.5) {
	fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGeneratedLimAcc, fWeight);
	AliDebug(3,"MC Lim Acc container filled\n");
      }	    
			
      //MC 
      fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepGenerated, fWeight);
      icountMC++;
      AliDebug(3,"MC cointainer filled \n");
			
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
		
    Bool_t unsetvtx=kFALSE;
    if(!charmCandidate->GetOwnPrimaryVtx()) {
      charmCandidate->SetOwnPrimaryVtx(aodVtx); // needed to compute all variables
      unsetvtx=kTRUE;
    }
		
    Bool_t signAssociation = cfVtxHF->SetRecoCandidateParam((AliAODRecoDecayHF*)charmCandidate);
    if (!signAssociation){
      charmCandidate = 0x0;
      continue;
    }

    Int_t isPartOrAntipart = cfVtxHF->CheckReflexion(fSign);
    if (isPartOrAntipart == 0){
      AliDebug(2, Form("The candidate pdg code doesn't match the requirement set in the task (fSign = %d)",fSign));
      continue;
    }

    AliDebug(3,Form("iCandid=%d - signAssociation=%d, isPartOrAntipart=%d",iCandid, signAssociation, isPartOrAntipart));

    Bool_t recoContFilled = cfVtxHF->FillRecoContainer(containerInput);
    if (recoContFilled){

      // weight according to pt
      if (fUseWeight){
	if (fFuncWeight){ // using user-defined function
	  AliDebug(2, "Using function");
	  fWeight = fFuncWeight->Eval(containerInput[0]);
	}
	else{ // using FONLL
	  AliDebug(2, "Using FONLL");
	  fWeight = GetWeight(containerInput[0]);
	}
	AliDebug(2, Form("pt = %f, weight = %f",containerInput[0], fWeight));
      }

      if (!fCuts->IsInFiducialAcceptance(containerInput[0],containerInput[1])) continue;
			
      //Reco Step
      Bool_t recoStep = cfVtxHF->RecoStep();
      Bool_t vtxCheck = fCuts->IsEventSelected(aodEvent);
			

      // Selection on the filtering bit
      Bool_t isBitSelected = kTRUE;
      if(fDecayChannel==2) {
	if(fUseSelectionBit && !charmCandidate->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)) isBitSelected = kFALSE;
      }else if(fDecayChannel==31){
	if(fUseSelectionBit && !charmCandidate->HasSelectionBit(AliRDHFCuts::kDplusCuts)) isBitSelected = kFALSE;
      }else if(fDecayChannel==33){
	if(fUseSelectionBit && !charmCandidate->HasSelectionBit(AliRDHFCuts::kDsCuts)) isBitSelected = kFALSE;
      }
      if(!isBitSelected) continue;



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
	      Bool_t keepLctoV0bachelor=ProcessLctoV0Bachelor(recoAnalysisCuts);
	      if (keepLctoV0bachelor) recoAnalysisCuts=3;
	    }

						    

	    fCuts->SetUsePID(iscutsusingpid); //restoring usage of the PID from the cuts object	
	    Bool_t tempAn=(recoAnalysisCuts == 3 || recoAnalysisCuts == isPartOrAntipart);
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
	      if (fDecayChannel == 32) tempPid=(recoPidSelection >0 || recoPidSelection == isPartOrAntipart);

	      if (tempPid){
		fCFManager->GetParticleContainer()->Fill(containerInput, kStepRecoPID, fWeight);
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
		continue;
	      }
	    }
	    else {
	      AliDebug(3, "PID selection not passed \n");
	      continue;
	    }
	  }
	  else{
	    AliDebug(3, "Number of ITS cluster step not passed\n");
	    continue;
	  }
	}
	else{
	  AliDebug(3, "Reco acceptance step not passed\n");
	  continue;
	}
      }
      else {
	AliDebug(3, "Reco step not passed\n");
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
void AliCFTaskVertexingHF::Terminate(Option_t*)
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
      titles[0]="pT_Lc (GeV/c)";
      titles[1]="rapidity";
      titles[2]="phi (rad)";
      titles[3]="cosPAV0";
      titles[4]="onTheFlyStatusV0";
      titles[5]="centrality";
      titles[6]="fake";
      titles[7]="multiplicity";
      titles[8]="pT_bachelor (GeV/c)";
      titles[9]="pT_V0pos (GeV/c)";
      titles[10]="pT_V0neg (GeV/c)";
      titles[11]="invMassV0 (GeV/c2)";
      titles[12]="dcaV0 (nSigma)";
      titles[13]="c#tauV0 (#mum)";
      titles[14]="c#tau (#mum)";
      titles[15]="cosPA";
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
      titles[0]="pT_candidate (GeV/c)";
      titles[1]="rapidity";
      titles[2]="phi (rad)";
      titles[3]="cosPAV0";
      titles[4]="onTheFlyStatusV0";
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

	
}

//___________________________________________________________________________
void AliCFTaskVertexingHF::UserCreateOutputObjects() 
{
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  AliPIDResponse *localPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();

  if (fCuts->GetIsUsePID() && fDecayChannel==22) {
	  
    fCuts->GetPidHF()->SetPidResponse(localPIDResponse);
    fCuts->GetPidHF()->SetOldPid(kFALSE);
    AliRDHFCutsLctoV0* lcv0Cuts=dynamic_cast<AliRDHFCutsLctoV0*>(fCuts);
    if(lcv0Cuts){
      lcv0Cuts->GetPidV0pos()->SetPidResponse(localPIDResponse);
      lcv0Cuts->GetPidV0neg()->SetPidResponse(localPIDResponse);
      lcv0Cuts->GetPidV0pos()->SetOldPid(kFALSE);
      lcv0Cuts->GetPidV0neg()->SetOldPid(kFALSE);
    }
  }

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
  Double_t pMeas=fHistoMeasNch->GetBinContent(fHistoMeasNch->FindBin(nch));
  Double_t pMC=fHistoMCNch->GetBinContent(fHistoMCNch->FindBin(nch));
  return pMeas/pMC;
}
//__________________________________________________________________________________________________
void AliCFTaskVertexingHF::CreateMeasuredNchHisto(){
  // creates historgam with measured multiplcity distribution in pp 7 TeV collisions (from Eur. Phys. J. C (2010) 68: 345–354)
  Double_t nchbins[66]={0.50,1.50,2.50,3.50,4.50,5.50,6.50,7.50,8.50,9.50,
			10.50,11.50,12.50,13.50,14.50,15.50,16.50,17.50,18.50,19.50,
			20.50,21.50,22.50,23.50,24.50,25.50,26.50,27.50,28.50,29.50,
			30.50,31.50,32.50,33.50,34.50,35.50,36.50,37.50,38.50,39.50,
			40.50,41.50,42.50,43.50,44.50,45.50,46.50,47.50,48.50,49.50,
			50.50,51.50,52.50,53.50,54.50,55.50,56.50,57.50,58.50,59.50,
			60.50,62.50,64.50,66.50,68.50,70.50};
  Double_t pch[65]={0.062011,0.072943,0.070771,0.067245,0.062834,0.057383,0.051499,0.04591,0.041109,0.036954,
		    0.03359,0.030729,0.028539,0.026575,0.024653,0.0229,0.021325,0.019768,0.018561,0.017187,
		    0.01604,0.014836,0.013726,0.012576,0.011481,0.010393,0.009502,0.008776,0.008024,0.007452,
		    0.006851,0.006428,0.00594,0.005515,0.005102,0.00469,0.004162,0.003811,0.003389,0.003071,
		    0.002708,0.002422,0.002184,0.001968,0.00186,0.00165,0.001577,0.001387,0.001254,0.001118,
		    0.001037,0.000942,0.000823,0.000736,0.000654,0.000579,0.000512,0.00049,0.00045,0.000355,
		    0.000296,0.000265,0.000193,0.00016,0.000126};

  if(fHistoMeasNch) delete fHistoMeasNch;
  fHistoMeasNch=new TH1F("hMeaseNch","",65,nchbins);
  for(Int_t i=0; i<65; i++){
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
