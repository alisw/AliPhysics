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
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"

ClassImp(AliRDHFCuts)

//--------------------------------------------------------------------------
AliRDHFCuts::AliRDHFCuts(const Char_t* name, const Char_t* title) : 
AliAnalysisCuts(name,title),
fMinVtxType(3),
fMinVtxContr(1),
fMaxVtxRedChi2(1e6),
fMinSPDMultiplicity(0),
fTriggerMask(0),
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
fUsePID(kFALSE)
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
  fMinSPDMultiplicity(source.fMinSPDMultiplicity),
  fTriggerMask(source.fTriggerMask),
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
  fUsePID(source.fUsePID)
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
  fMinSPDMultiplicity=source.fMinSPDMultiplicity;
  fTriggerMask=source.fTriggerMask;
  fnPtBins=source.fnPtBins;
  fnVars=source.fnVars;
  fGlobalIndex=source.fGlobalIndex;
  fnVarsForOpt=source.fnVarsForOpt;
  fUsePID=source.fUsePID;

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

}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsEventSelected(AliVEvent *event) const {
  //
  // Event selection
  // 
  //if(fTriggerMask && event->GetTriggerMask()!=fTriggerMask) return kFALSE;

  // multiplicity cuts no implemented yet

  const AliVVertex *vertex = event->GetPrimaryVertex();

  if(!vertex) return kFALSE;

  TString title=vertex->GetTitle();
  if(title.Contains("Z") && fMinVtxType>1) return kFALSE; 
  if(title.Contains("3D") && fMinVtxType>2) return kFALSE; 

  if(vertex->GetNContributors()<fMinVtxContr) return kFALSE; 

  return kTRUE;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::AreDaughtersSelected(AliAODRecoDecayHF *d) const {
  //
  // Daughter tracks selection
  // 
  if(!fTrackCuts) return kTRUE;

  Int_t ndaughters = d->GetNDaughters();
  AliESDtrack* esdTrack=0;
  AliAODVertex *vAOD = d->GetPrimaryVtx();
  Double_t pos[3],cov[6];
  vAOD->GetXYZ(pos);
  vAOD->GetCovarianceMatrix(cov);
  const AliESDVertex *vESD = new AliESDVertex(pos,cov,100.,100);

  Bool_t retval=kTRUE;

  for(Int_t idg=0; idg<ndaughters; idg++) {
    AliAODTrack *dgTrack = (AliAODTrack*)d->GetDaughter(idg);
    if(!dgTrack) retval = kFALSE;
    //printf("charge %d\n",dgTrack->Charge());
    if(dgTrack->Charge()==0) continue; // it's not a track, but a V0
    // convert to ESD track here
    esdTrack=new AliESDtrack(dgTrack);
    // needed to calculate the impact parameters
    esdTrack->RelateToVertex(vESD,0.,3.); 
    if(!fTrackCuts->IsSelected(esdTrack)) retval = kFALSE;
    delete esdTrack; esdTrack=0;
  }

  delete vESD; vESD=0;

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
  cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
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
    printf("Changing the variable names\n");
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
    printf("Changing the variables for cut optimization\n");
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
void AliRDHFCuts::SetCuts(Int_t nVars,Int_t nPtBins,Float_t **cutsRD) {
  //
  // store the cuts
  //
  if(nVars!=fnVars) {
    printf("Wrong number of variables: it has to be %d\n",fnVars);
    return;
  } 
  if(nPtBins!=fnPtBins) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBins);
    return;
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
    return;
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

