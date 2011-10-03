/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// Base class for cuts on AOD reconstructed heavy-flavour decay
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////
#include <Riostream.h>

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliAODRecoDecayHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"

ClassImp(AliRDHFCuts)

//--------------------------------------------------------------------------
AliRDHFCuts::AliRDHFCuts(const Char_t* name, const Char_t* title) : 
AliAnalysisCuts(name,title),
fMinVtxType(3),
fMinVtxContr(1),
fMaxVtxRedChi2(1e6),
fMaxVtxZ(10.),
fMinSPDMultiplicity(0),
fTriggerMask(0),
fTriggerClass("CINT1"),
fTrackCuts(0),
fnPtBins(1),
fnPtBinLimits(1),
fPtBinLimits(0),
fnVars(1),
fVarNames(0),
fnVarsForOpt(0),
fVarsForOpt(0),
fGlobalIndex(1),
fCutsRD(0),
fIsUpperCut(0),
fUsePID(kFALSE),
fUseAOD049(kFALSE),
fPidHF(0),
fWhyRejection(0),
fEvRejectionBits(0),
fRemoveDaughtersFromPrimary(kFALSE),
fUseMCVertex(kFALSE),
fOptPileup(0),
fMinContrPileup(3),
fMinDzPileup(0.6),
fUseCentrality(0),
fMinCentrality(0.),
fMaxCentrality(100.),
fFixRefs(kFALSE),
fIsSelectedCuts(0),
fIsSelectedPID(0),
fMinPtCand(-1.),
fMaxPtCand(100000.),
fKeepSignalMC(kFALSE)
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliRDHFCuts::AliRDHFCuts(const AliRDHFCuts &source) :
  AliAnalysisCuts(source),
  fMinVtxType(source.fMinVtxType),
  fMinVtxContr(source.fMinVtxContr),
  fMaxVtxRedChi2(source.fMaxVtxRedChi2),
  fMaxVtxZ(source.fMaxVtxZ),
  fMinSPDMultiplicity(source.fMinSPDMultiplicity),
  fTriggerMask(source.fTriggerMask),
  fTriggerClass(source.fTriggerClass),
  fTrackCuts(0),
  fnPtBins(source.fnPtBins),
  fnPtBinLimits(source.fnPtBinLimits),
  fPtBinLimits(0),
  fnVars(source.fnVars),
  fVarNames(0),
  fnVarsForOpt(source.fnVarsForOpt),
  fVarsForOpt(0),
  fGlobalIndex(source.fGlobalIndex),
  fCutsRD(0),
  fIsUpperCut(0),
  fUsePID(source.fUsePID),
  fUseAOD049(source.fUseAOD049),
  fPidHF(0),
  fWhyRejection(source.fWhyRejection),
  fEvRejectionBits(source.fEvRejectionBits),
  fRemoveDaughtersFromPrimary(source.fRemoveDaughtersFromPrimary),
  fUseMCVertex(source.fUseMCVertex),
  fOptPileup(source.fOptPileup),
  fMinContrPileup(source.fMinContrPileup),
  fMinDzPileup(source.fMinDzPileup),
  fUseCentrality(source.fUseCentrality),
  fMinCentrality(source.fMinCentrality),
  fMaxCentrality(source.fMaxCentrality),
  fFixRefs(source.fFixRefs),
  fIsSelectedCuts(source.fIsSelectedCuts),
  fIsSelectedPID(source.fIsSelectedPID),
  fMinPtCand(source.fMinPtCand),
  fMaxPtCand(source.fMaxPtCand),
  fKeepSignalMC(source.fKeepSignalMC)
{
  //
  // Copy constructor
  //
  cout<<"Copy constructor"<<endl;
  if(source.GetTrackCuts()) AddTrackCuts(source.GetTrackCuts());
  if(source.fPtBinLimits) SetPtBins(source.fnPtBinLimits,source.fPtBinLimits);
  if(source.fVarNames) SetVarNames(source.fnVars,source.fVarNames,source.fIsUpperCut);
  if(source.fCutsRD) SetCuts(source.fGlobalIndex,source.fCutsRD);
  if(source.fVarsForOpt) SetVarsForOpt(source.fnVarsForOpt,source.fVarsForOpt);
  if(source.fPidHF) SetPidHF(source.fPidHF);
  PrintAll();

}
//--------------------------------------------------------------------------
AliRDHFCuts &AliRDHFCuts::operator=(const AliRDHFCuts &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAnalysisCuts::operator=(source);

  fMinVtxType=source.fMinVtxType;
  fMinVtxContr=source.fMinVtxContr;
  fMaxVtxRedChi2=source.fMaxVtxRedChi2;
  fMaxVtxZ=source.fMaxVtxZ;
  fMinSPDMultiplicity=source.fMinSPDMultiplicity;
  fTriggerMask=source.fTriggerMask;
  fTriggerClass=source.fTriggerClass;
  fnPtBins=source.fnPtBins;
  fnVars=source.fnVars;
  fGlobalIndex=source.fGlobalIndex;
  fnVarsForOpt=source.fnVarsForOpt;
  fUsePID=source.fUsePID;
  fUseAOD049=source.fUseAOD049;
  SetPidHF(source.GetPidHF());
  fWhyRejection=source.fWhyRejection;
  fEvRejectionBits=source.fEvRejectionBits;
  fRemoveDaughtersFromPrimary=source.fRemoveDaughtersFromPrimary;
  fUseMCVertex=source.fUseMCVertex;
  fOptPileup=source.fOptPileup;
  fMinContrPileup=source.fMinContrPileup;
  fMinDzPileup=source.fMinDzPileup;
  fUseCentrality=source.fUseCentrality;
  fMinCentrality=source.fMinCentrality;
  fMaxCentrality=source.fMaxCentrality;
  fFixRefs=source.fFixRefs;
  fIsSelectedCuts=source.fIsSelectedCuts;
  fIsSelectedPID=source.fIsSelectedPID;
  fMinPtCand=source.fMinPtCand;
  fMaxPtCand=source.fMaxPtCand;
  fKeepSignalMC=source.fKeepSignalMC;

  if(source.GetTrackCuts()) AddTrackCuts(source.GetTrackCuts());
  if(source.fPtBinLimits) SetPtBins(source.fnPtBinLimits,source.fPtBinLimits);
  if(source.fVarNames) SetVarNames(source.fnVars,source.fVarNames,source.fIsUpperCut);
  if(source.fCutsRD) SetCuts(source.fGlobalIndex,source.fCutsRD);
  if(source.fVarsForOpt) SetVarsForOpt(source.fnVarsForOpt,source.fVarsForOpt);
  PrintAll();

  return *this;
}
//--------------------------------------------------------------------------
AliRDHFCuts::~AliRDHFCuts() {
  //  
  // Default Destructor
  //
  if(fPtBinLimits) {delete [] fPtBinLimits; fPtBinLimits=0;}
  if(fVarNames) {delete [] fVarNames; fVarNames=0;}
  if(fVarsForOpt) {delete [] fVarsForOpt; fVarsForOpt=0;}
  if(fCutsRD) {
    delete [] fCutsRD;
    fCutsRD=0;
  }
  if(fIsUpperCut) {delete [] fIsUpperCut; fIsUpperCut=0;}
  if(fPidHF){ 
    delete fPidHF;
    fPidHF=0;
  }
}
//---------------------------------------------------------------------------
Int_t AliRDHFCuts::IsEventSelectedInCentrality(AliVEvent *event) {
  //
  // Centrality selection
  //
  if(fUseCentrality<kCentOff||fUseCentrality>=kCentInvalid){    
    AliWarning("Centrality estimator not valid");    
    return 3;  
  }else{    
    Float_t centvalue=GetCentrality((AliAODEvent*)event);          
    if (centvalue<-998.){//-999 if no centralityP
      return 0;    
    }else{      
      if (centvalue<fMinCentrality || centvalue>fMaxCentrality){
	return 2;      
      }    
    } 
  }  
  return 0;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsEventSelected(AliVEvent *event) {
  //
  // Event selection
  // 
  //if(fTriggerMask && event->GetTriggerMask()!=fTriggerMask) return kFALSE;

  fWhyRejection=0;
  fEvRejectionBits=0;
  Bool_t accept=kTRUE;

  // check if it's MC
  Bool_t isMC=kFALSE;
  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(mcArray) {isMC=kTRUE;fUseAOD049=kFALSE;}

  // settings for the TPC dE/dx BB parameterization
  if(fPidHF && fPidHF->GetOldPid()) {
    // pp, from LHC10d onwards
    if((event->GetRunNumber()>121693 && event->GetRunNumber()<136851) ||
       event->GetRunNumber()>139517) fPidHF->SetOnePad(kTRUE);
    // pp, 2011 low energy run
    if((event->GetRunNumber()>=146686 && event->GetRunNumber()<=146860)){
      fPidHF->SetppLowEn2011(kTRUE);
      fPidHF->SetOnePad(kFALSE);
    }
    // PbPb LHC10h
    if(event->GetRunNumber()>=136851 && event->GetRunNumber()<=139517) fPidHF->SetPbPb(kTRUE);
    // MC
    if(isMC) fPidHF->SetMC(kTRUE);
    if(isMC && (event->GetRunNumber()>=146686 && event->GetRunNumber()<=146860))
      fPidHF->SetMClowenpp2011(kTRUE);
  } 
  else if(fPidHF && !fPidHF->GetOldPid()) {
   if(fPidHF->GetPidResponse()==0x0){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
    fPidHF->SetPidResponse(pidResp);
   }
  }


  // trigger class
  TString firedTriggerClasses=((AliAODEvent*)event)->GetFiredTriggerClasses();
  // don't do for MC and for PbPb 2010 data
  if(!isMC && (event->GetRunNumber()<136851 || event->GetRunNumber()>139517)) {
    if(!firedTriggerClasses.Contains(fTriggerClass.Data())) {
      fWhyRejection=5;
      fEvRejectionBits+=1<<kNotSelTrigger;
      accept=kFALSE;
    }
  }

  // TEMPORARY FIX FOR REFERENCES
  // Fix references to daughter tracks
  //  if(fFixRefs) {
  //    AliAnalysisVertexingHF *fixer = new AliAnalysisVertexingHF();
  //    fixer->FixReferences((AliAODEvent*)event);
  //    delete fixer;
  //  }
  //



  // vertex requirements
   
  const AliVVertex *vertex = event->GetPrimaryVertex();

  if(!vertex){
    accept=kFALSE;
    fEvRejectionBits+=1<<kNoVertex;
  }else{
    TString title=vertex->GetTitle();
    if(title.Contains("Z") && fMinVtxType>1){
      accept=kFALSE;
      fEvRejectionBits+=1<<kNoVertex;
    }
    else if(title.Contains("3D") && fMinVtxType>2){
      accept=kFALSE;
      fEvRejectionBits+=1<<kNoVertex;
    }
    if(vertex->GetNContributors()<fMinVtxContr){
      accept=kFALSE;
      fEvRejectionBits+=1<<kTooFewVtxContrib;
    }
    if(TMath::Abs(vertex->GetZ())>fMaxVtxZ) {
      fEvRejectionBits+=1<<kZVtxOutFid;
      if(accept) fWhyRejection=6;
      accept=kFALSE;
    } 
  }


  // pile-up rejection
  if(fOptPileup==kRejectPileupEvent){
    Int_t cutc=(Int_t)fMinContrPileup;
    Double_t cutz=(Double_t)fMinDzPileup;
    if(event->IsPileupFromSPD(cutc,cutz,3.,2.,10.)) {
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kPileupSPD;
      accept=kFALSE;
    }
  }

  // centrality selection
  if (fUseCentrality!=kCentOff) {  
    Int_t rejection=IsEventSelectedInCentrality(event);    
    if(rejection>1){      
      if(accept) fWhyRejection=rejection;      
      fEvRejectionBits+=1<<kOutsideCentrality;
      accept=kFALSE;
    }
  }

  return accept;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::AreDaughtersSelected(AliAODRecoDecayHF *d) const {
  //
  // Daughter tracks selection
  // 
  if(!fTrackCuts) return kTRUE;
 
  Int_t ndaughters = d->GetNDaughters();
  AliAODVertex *vAOD = d->GetPrimaryVtx();
  Double_t pos[3],cov[6];
  vAOD->GetXYZ(pos);
  vAOD->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  
  Bool_t retval=kTRUE;
  
  for(Int_t idg=0; idg<ndaughters; idg++) {
    AliAODTrack *dgTrack = (AliAODTrack*)d->GetDaughter(idg);
    if(!dgTrack) {retval = kFALSE; continue;}
    //printf("charge %d\n",dgTrack->Charge());
    if(dgTrack->Charge()==0) continue; // it's not a track, but a V0

    if(!IsDaughterSelected(dgTrack,&vESD,fTrackCuts)) retval = kFALSE;
  }
  
  return retval;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsDaughterSelected(AliAODTrack *track,const AliESDVertex *primary,AliESDtrackCuts *cuts) const {
  //
  // Convert to ESDtrack, relate to vertex and check cuts
  //
  if(!cuts) return kTRUE;

  Bool_t retval=kTRUE;

  // convert to ESD track here
  AliESDtrack esdTrack(track);
  // set the TPC cluster info
  esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(track->GetTPCNclsF());
  // needed to calculate the impact parameters
  esdTrack.RelateToVertex(primary,0.,3.); 

  if(!cuts->IsSelected(&esdTrack)) retval = kFALSE;
 
  if(fOptPileup==kRejectTracksFromPileupVertex){
    // to be implemented
    // we need either to have here the AOD Event, 
    // or to have the pileup vertex object
  }
  return retval; 
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetPtBins(Int_t nPtBinLimits,Float_t *ptBinLimits) {
  // Set the pt bins

  if(fPtBinLimits) {
    delete [] fPtBinLimits;
    fPtBinLimits = NULL;
    printf("Changing the pt bins\n");
  }

  if(nPtBinLimits != fnPtBins+1){
    cout<<"Warning: ptBinLimits dimention "<<nPtBinLimits<<" != nPtBins+1 ("<<fnPtBins+1<<")\nSetting nPtBins to "<<nPtBinLimits-1<<endl;
    SetNPtBins(nPtBinLimits-1);
  }

  fnPtBinLimits = nPtBinLimits;
  SetGlobalIndex();
  //cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
  fPtBinLimits = new Float_t[fnPtBinLimits];
  for(Int_t ib=0; ib<nPtBinLimits; ib++) fPtBinLimits[ib]=ptBinLimits[ib];

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetVarNames(Int_t nVars,TString *varNames,Bool_t *isUpperCut){
  // Set the variable names

  if(fVarNames) {
    delete [] fVarNames;
    fVarNames = NULL;
    //printf("Changing the variable names\n");
  }
  if(nVars!=fnVars){
    printf("Wrong number of variables: it has to be %d\n",fnVars);
    return;
  }
  //fnVars=nVars;
  fVarNames = new TString[nVars];
  fIsUpperCut = new Bool_t[nVars];
  for(Int_t iv=0; iv<nVars; iv++) {
    fVarNames[iv] = varNames[iv];
    fIsUpperCut[iv] = isUpperCut[iv];
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetVarsForOpt(Int_t nVars,Bool_t *forOpt) {
  // Set the variables to be used for cuts optimization

  if(fVarsForOpt) {
    delete [] fVarsForOpt;
    fVarsForOpt = NULL;
    //printf("Changing the variables for cut optimization\n");
  }
  
  if(nVars==0){//!=fnVars) {
    printf("%d not accepted as number of variables: it has to be %d\n",nVars,fnVars);
    return;
  } 
  
  fnVarsForOpt = 0;
  fVarsForOpt = new Bool_t[fnVars];
  for(Int_t iv=0; iv<fnVars; iv++) {
    fVarsForOpt[iv]=forOpt[iv];
    if(fVarsForOpt[iv]) fnVarsForOpt++;
  }

  return;
}

//---------------------------------------------------------------------------
void AliRDHFCuts::SetUseCentrality(Int_t flag) {
  //
  // set centrality estimator  
  //
  fUseCentrality=flag;
  if(fUseCentrality<kCentOff||fUseCentrality>=kCentInvalid) AliWarning("Centrality estimator not valid");
 
  return;
}


//---------------------------------------------------------------------------
void AliRDHFCuts::SetCuts(Int_t nVars,Int_t nPtBins,Float_t **cutsRD) {
  //
  // store the cuts
  //
  if(nVars!=fnVars) {
    printf("Wrong number of variables: it has to be %d\n",fnVars);
    AliFatal("exiting");
  } 
  if(nPtBins!=fnPtBins) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBins);
    AliFatal("exiting");
  } 

  if(!fCutsRD)  fCutsRD = new Float_t[fGlobalIndex];
  

  for(Int_t iv=0; iv<fnVars; iv++) {

    for(Int_t ib=0; ib<fnPtBins; ib++) {

      //check
      if(GetGlobalIndex(iv,ib)>=fGlobalIndex) {
	cout<<"Overflow, exit..."<<endl;
	return;
      }

      fCutsRD[GetGlobalIndex(iv,ib)] = cutsRD[iv][ib];

    }
  }
  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetCuts(Int_t glIndex,Float_t* cutsRDGlob){
  //
  // store the cuts
  //
  if(glIndex != fGlobalIndex){
    cout<<"Wrong array size: it has to be "<<fGlobalIndex<<endl;
    AliFatal("exiting");
  }
  if(!fCutsRD)  fCutsRD = new Float_t[fGlobalIndex];

  for(Int_t iGl=0;iGl<fGlobalIndex;iGl++){
    fCutsRD[iGl] = cutsRDGlob[iGl];
  }
  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::PrintAll() const {
  //
  // print all cuts values
  // 

  printf("Minimum vtx type %d\n",fMinVtxType);
  printf("Minimum vtx contr %d\n",fMinVtxContr);
  printf("Max vtx red chi2 %f\n",fMaxVtxRedChi2);
  printf("Min SPD mult %d\n",fMinSPDMultiplicity);
  printf("Use PID %d\n",(Int_t)fUsePID);
  printf("Remove daughters from vtx %d\n",(Int_t)fRemoveDaughtersFromPrimary);
  printf("Pileup rejection: %s\n",(fOptPileup > 0) ? "Yes" : "No");
  if(fOptPileup==1) printf(" -- Reject pileup event");
  if(fOptPileup==2) printf(" -- Reject tracks from pileup vtx");
  if(fUseCentrality>0) {
    TString estimator="";
    if(fUseCentrality==1) estimator = "V0";
    if(fUseCentrality==2) estimator = "Tracks";
    if(fUseCentrality==3) estimator = "Tracklets";
    if(fUseCentrality==4) estimator = "SPD clusters outer"; 
    printf("Centrality class considered: %.1f-%.1f, estimated with %s",fMinCentrality,fMaxCentrality,estimator.Data());
  }

  if(fVarNames){
    cout<<"Array of variables"<<endl;
    for(Int_t iv=0;iv<fnVars;iv++){
      cout<<fVarNames[iv]<<"\t";
    }
    cout<<endl;
  }
  if(fVarsForOpt){
    cout<<"Array of optimization"<<endl;
    for(Int_t iv=0;iv<fnVars;iv++){
      cout<<fVarsForOpt[iv]<<"\t";
    }
    cout<<endl;
  }
  if(fIsUpperCut){
    cout<<"Array of upper/lower cut"<<endl;
   for(Int_t iv=0;iv<fnVars;iv++){
     cout<<fIsUpperCut[iv]<<"\t";
   }
   cout<<endl;
  }
  if(fPtBinLimits){
    cout<<"Array of ptbin limits"<<endl;
    for(Int_t ib=0;ib<fnPtBinLimits;ib++){
      cout<<fPtBinLimits[ib]<<"\t";
    }
    cout<<endl;
  }
  if(fCutsRD){
    cout<<"Matrix of cuts"<<endl;
   for(Int_t iv=0;iv<fnVars;iv++){
     for(Int_t ib=0;ib<fnPtBins;ib++){
       cout<<"fCutsRD["<<iv<<"]["<<ib<<"] = "<<fCutsRD[GetGlobalIndex(iv,ib)]<<"\t";
     } 
     cout<<endl;
   }
   cout<<endl;
  }
  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::GetCuts(Float_t**& cutsRD) const{
  //
  // get the cuts
  //

  //cout<<"Give back a "<<fnVars<<"x"<<fnPtBins<<" matrix."<<endl;


  Int_t iv,ib;
  if(!cutsRD) {
    //cout<<"Initialization..."<<endl;
    cutsRD=new Float_t*[fnVars];
    for(iv=0; iv<fnVars; iv++) {
      cutsRD[iv] = new Float_t[fnPtBins];
    }
  }
  
  for(Int_t iGlobal=0; iGlobal<fGlobalIndex; iGlobal++) {
    GetVarPtIndex(iGlobal,iv,ib);
    cutsRD[iv][ib] = fCutsRD[iGlobal];
  }

  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCuts::GetGlobalIndex(Int_t iVar,Int_t iPtBin) const{
  //
  // give the global index from variable and pt bin
  //
  return iPtBin*fnVars+iVar;
}

//---------------------------------------------------------------------------
void AliRDHFCuts::GetVarPtIndex(Int_t iGlob, Int_t& iVar, Int_t& iPtBin) const {
  //
  //give the index of the variable and of the pt bin from the global index
  //
  iPtBin=(Int_t)iGlob/fnVars;
  iVar=iGlob%fnVars;

  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCuts::PtBin(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fPtBinLimits[0])return ptbin;
  for (Int_t i=0;i<fnPtBins;i++){
    if(pt<fPtBinLimits[i+1]) {
      ptbin=i;
      break;
    }
  }
  return ptbin;
}
//-------------------------------------------------------------------
Float_t AliRDHFCuts::GetCutValue(Int_t iVar,Int_t iPtBin) const {
  // 
  // Give the value of cut set for the variable iVar and the pt bin iPtBin
  //
  if(!fCutsRD){
    cout<<"Cuts not iniziaisez yet"<<endl;
    return 0;
  }
  return fCutsRD[GetGlobalIndex(iVar,iPtBin)];
}
//-------------------------------------------------------------------
Float_t AliRDHFCuts::GetCentrality(AliAODEvent* aodEvent,AliRDHFCuts::ECentrality estimator) {
  //
  // Get centrality percentile
  //

  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)aodEvent)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(mcArray) {fUseAOD049=kFALSE;}

  AliAODHeader *header=aodEvent->GetHeader();
  AliCentrality *centrality=header->GetCentralityP();
  Float_t cent=-999.;
  Bool_t isSelRun=kFALSE;
  Int_t selRun[5]={138364, 138826, 138828, 138836, 138871};
  if(!centrality) return cent;
  else{
    if (estimator==kCentV0M){
      cent=(Float_t)(centrality->GetCentralityPercentile("V0M"));
      if(cent<0){
	Int_t quality = centrality->GetQuality();
	if(quality<=1){
	  cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
	}else{
	  Int_t runnum=aodEvent->GetRunNumber();
	  for(Int_t ir=0;ir<5;ir++){
	    if(runnum==selRun[ir]){
	      isSelRun=kTRUE;
	      break;
	    }
	  }
	  if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
	}
      }

      //temporary fix for AOD049 outliers
      if(fUseAOD049&&cent>=0){
	Float_t v0=0;
	AliAODVZERO* aodV0 = aodEvent->GetVZEROData();
	v0+=aodV0->GetMTotV0A();
	v0+=aodV0->GetMTotV0C();
	if(cent==0&&v0<19500)return -1;//filtering issue
	Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
	Float_t val= 1.30552 +  0.147931 * v0;
	Float_t tklSigma[101]={176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86, 120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654, 92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334, 68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224, 51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255, 37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398, 26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235, 19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504, 12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544};
	if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )return -1;//outlier
      }
    }
    else {
      if (estimator==kCentTRK) {
	cent=(Float_t)(centrality->GetCentralityPercentile("TRK"));
	if(cent<0){
	  Int_t quality = centrality->GetQuality();
	  if(quality<=1){
	    cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TRK");
	  }else{
	    Int_t runnum=aodEvent->GetRunNumber();
	    for(Int_t ir=0;ir<5;ir++){
	      if(runnum==selRun[ir]){
		isSelRun=kTRUE;
		break;
	      }
	    }
	    if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TRK");
	  }
	}
      }
      else{
	if (estimator==kCentTKL){
	  cent=(Float_t)(centrality->GetCentralityPercentile("TKL"));
	  if(cent<0){
	    Int_t quality = centrality->GetQuality();
	    if(quality<=1){
	      cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TKL");
	    }else{
	      Int_t runnum=aodEvent->GetRunNumber();
	      for(Int_t ir=0;ir<5;ir++){
		if(runnum==selRun[ir]){
		  isSelRun=kTRUE;
		  break;
	    }
	      }
	      if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TKL");
	    }   
	  }
	}
	else{
	  if (estimator==kCentCL1){
	    cent=(Float_t)(centrality->GetCentralityPercentile("CL1"));
	    if(cent<0){
	      Int_t quality = centrality->GetQuality();
	      if(quality<=1){
		cent=(Float_t)centrality->GetCentralityPercentileUnchecked("CL1");
	      }else{
		Int_t runnum=aodEvent->GetRunNumber();
		for(Int_t ir=0;ir<5;ir++){
		  if(runnum==selRun[ir]){
		    isSelRun=kTRUE;
		    break;
		  }
		}
		if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("CL1");
	      }
	    }
	  }
	  else {
	    AliWarning("Centrality estimator not valid");
	    
	  }
	}
      }
    } 
  }
  return cent;
}
//-------------------------------------------------------------------
Bool_t AliRDHFCuts::CompareCuts(const AliRDHFCuts *obj) const {
  //
  // Compare two cuts objects
  //

  Bool_t areEqual=kTRUE;

  if(fMinVtxType!=obj->fMinVtxType) { printf("Minimum vtx type %d  %d\n",fMinVtxType,obj->fMinVtxType); areEqual=kFALSE;}

  if(fMinVtxContr!=obj->fMinVtxContr) { printf("Minimum vtx contr %d  %d\n",fMinVtxContr,obj->fMinVtxContr); areEqual=kFALSE;}

  if(TMath::Abs(fMaxVtxRedChi2-obj->fMaxVtxRedChi2)>1.e-10) {   printf("Max vtx red chi2 %f  %f\n",fMaxVtxRedChi2,obj->fMaxVtxRedChi2);areEqual=kFALSE;}

  if(fMinSPDMultiplicity!=obj->fMinSPDMultiplicity) {  printf("Min SPD mult %d\n  %d",fMinSPDMultiplicity,obj->fMinSPDMultiplicity);areEqual=kFALSE;}

  if(fUsePID!=obj->fUsePID) { printf("Use PID %d  %d\n",(Int_t)fUsePID,(Int_t)obj->fUsePID); areEqual=kFALSE;}

  if(fRemoveDaughtersFromPrimary!=obj->fRemoveDaughtersFromPrimary) {printf("Remove daughters from vtx %d  %d\n",(Int_t)fRemoveDaughtersFromPrimary,(Int_t)obj->fRemoveDaughtersFromPrimary); areEqual=kFALSE;}
  if(fTrackCuts){
    if(fTrackCuts->GetMinNClusterTPC()!=obj->fTrackCuts->GetMinNClusterTPC()) {printf("MinNClsTPC %d  %d\n",fTrackCuts->GetMinNClusterTPC(),obj->fTrackCuts->GetMinNClusterTPC()); areEqual=kFALSE;}

    if(fTrackCuts->GetMinNClustersITS()!=obj->fTrackCuts->GetMinNClustersITS()) {printf("MinNClsITS %d  %d\n",fTrackCuts->GetMinNClustersITS(),obj->fTrackCuts->GetMinNClustersITS()); areEqual=kFALSE;}

    if(TMath::Abs(fTrackCuts->GetMaxChi2PerClusterTPC()-obj->fTrackCuts->GetMaxChi2PerClusterTPC())>1.e-10) {printf("MaxChi2ClsTPC %f  %f\n",fTrackCuts->GetMaxChi2PerClusterTPC(),obj->fTrackCuts->GetMaxChi2PerClusterTPC()); areEqual=kFALSE;}

    if(fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)!=obj->fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)) {printf("ClusterReq SPD %d  %d\n",fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD),obj->fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)); areEqual=kFALSE;}
  }

  if(fCutsRD) {
   for(Int_t iv=0;iv<fnVars;iv++) {
     for(Int_t ib=0;ib<fnPtBins;ib++) {
       if(TMath::Abs(fCutsRD[GetGlobalIndex(iv,ib)]-obj->fCutsRD[GetGlobalIndex(iv,ib)])>1.e-10) {
	 cout<<"fCutsRD["<<iv<<"]["<<ib<<"] = "<<fCutsRD[GetGlobalIndex(iv,ib)]<<"   "<<obj->fCutsRD[GetGlobalIndex(iv,ib)]<<"\n";
	 areEqual=kFALSE;
       }
     }
   }
  }

  return areEqual;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::MakeTable() const {
  //
  // print cuts values in table format
  // 

	TString ptString = "pT range";
	if(fVarNames && fPtBinLimits && fCutsRD){
		TString firstLine(Form("*       %-15s",ptString.Data()));
		for (Int_t ivar=0; ivar<fnVars; ivar++){
			firstLine+=Form("*    %-15s  ",fVarNames[ivar].Data());
			if (ivar == fnVars){
				firstLine+="*\n";
			}
		}
		Printf("%s",firstLine.Data());
		
		for (Int_t ipt=0; ipt<fnPtBins; ipt++){
			TString line;
			if (ipt==fnPtBins-1){
				line=Form("*  %5.1f < pt < inf    ",fPtBinLimits[ipt]);
			}
			else{
				line=Form("*  %5.1f < pt < %4.1f   ",fPtBinLimits[ipt],fPtBinLimits[ipt+1]);
			}
			for (Int_t ivar=0; ivar<fnVars; ivar++){
				line+=Form("*     %-15f ",fCutsRD[GetGlobalIndex(ivar,ipt)]);
			}
			Printf("%s",line.Data());
		}

	}


  return;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFCuts::RecalcOwnPrimaryVtx(AliAODRecoDecayHF *d,
					AliAODEvent *aod) const
{
  //
  // Recalculate primary vertex without daughters
  //

  if(!aod) {
    AliError("Can not remove daughters from vertex without AOD event");
    return 0;
  }   

  AliAODVertex *recvtx=d->RemoveDaughtersFromPrimaryVtx(aod);
  if(!recvtx){
    AliDebug(2,"Removal of daughter tracks failed");
    return kFALSE;
  }


  //set recalculed primary vertex
  d->SetOwnPrimaryVtx(recvtx);
  delete recvtx;

  return kTRUE;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFCuts::SetMCPrimaryVtx(AliAODRecoDecayHF *d,AliAODEvent *aod) const
{
  //
  // Recalculate primary vertex without daughters
  //

  if(!aod) {
    AliError("Can not get MC vertex without AOD event");
    return kFALSE;
  }   

  // load MC header
  AliAODMCHeader *mcHeader = 
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    AliError("Can not get MC vertex without AODMCHeader event");
    return kFALSE;
  }
  Double_t pos[3];
  Double_t covmatrix[6]={0.,0.,0.,0.,0.,0.};
  mcHeader->GetVertex(pos);
  AliAODVertex *recvtx=new AliAODVertex(pos,covmatrix);

  if(!recvtx){
    AliDebug(2,"Removal of daughter tracks failed");
    return kFALSE;
  }

  //set recalculed primary vertex
  d->SetOwnPrimaryVtx(recvtx);

  d->RecalculateImpPars(recvtx,aod);

  delete recvtx;

  return kTRUE;
}
//--------------------------------------------------------------------------
void AliRDHFCuts::CleanOwnPrimaryVtx(AliAODRecoDecayHF *d,
				     AliAODEvent *aod,
				     AliAODVertex *origownvtx) const
{
  //
  // Clean-up own primary vertex if needed
  //

  if(fRemoveDaughtersFromPrimary || fUseMCVertex) {
    d->UnsetOwnPrimaryVtx();
    if(origownvtx) {
      d->SetOwnPrimaryVtx(origownvtx);
      delete origownvtx; origownvtx=NULL;
    }
    d->RecalculateImpPars(d->GetPrimaryVtx(),aod);
  } else {
    if(origownvtx) {
      delete origownvtx; origownvtx=NULL;
    }
  }
  return;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsSignalMC(AliAODRecoDecay *d,AliAODEvent *aod,Int_t pdg) const 
{
  //
  // Checks if this candidate is matched to MC signal
  // 

  if(!aod) return kFALSE;

  // get the MC array
  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)aod)->GetList()->FindObject(AliAODMCParticle::StdBranchName());

  if(!mcArray) return kFALSE;

  // try to match  
  Int_t label = d->MatchToMC(pdg,mcArray);
  
  if(label>=0) {
    //printf("MATCH!\n");
    return kTRUE;
  }

  return kFALSE;
}


