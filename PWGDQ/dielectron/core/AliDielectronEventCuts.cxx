/*************************************************************************
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron EventCuts                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include <AliTriggerAnalysis.h>
#include <AliESDVertex.h>
#include <AliAODVertex.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMultiplicity.h>
#include <AliCentrality.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>

#include "AliDielectronVarManager.h"
#include "AliDielectronEventCuts.h"

ClassImp(AliDielectronEventCuts)

const char* AliDielectronEventCuts::fgkVtxNames[AliDielectronEventCuts::kVtxTracksOrSPD+1] = {"Tracks", "SPD", "TPC", "Any", "TracksOrSPD"};

AliDielectronEventCuts::AliDielectronEventCuts() :
  AliAnalysisCuts(),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fRun(),
  fVtxZmin(0.),
  fVtxZmax(0.),
  fRequireVtx(kFALSE),
  fMinVtxContributors(0),
  fMultITSTPC(kFALSE),
  fCentMin(1.),
  fCentMax(0.),
  fVtxType(kVtxTracks),
  fRequire13sel(kFALSE),
  fUtils(),
  fRequireV0and(0),
  fTriggerAnalysis(0x0),
  fkVertex(0x0),
  fkVertexAOD(0x0),
  fparMean(0x0),
  fparSigma(0x0),
  fcutSigma(3.),
  fparMinVtxContributors(0x0),
  fparMaxVtxContributors(0x0)
{
  //
  // Default Constructor
  //
  for (Int_t icut=0; icut<5; ++icut){
    fCorrCutMin[icut]=0x0;
    fCorrCutMax[icut]=0x0;
  }
}

//______________________________________________
AliDielectronEventCuts::AliDielectronEventCuts(const char* name, const char* title) :
  AliAnalysisCuts(name, title),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fRun(),
  fVtxZmin(0.),
  fVtxZmax(0.),
  fRequireVtx(kFALSE),
  fMinVtxContributors(0),
  fMultITSTPC(kFALSE),
  fCentMin(1.),
  fCentMax(0.),
  fVtxType(kVtxTracks),
  fRequire13sel(kFALSE),
  fUtils(),
  fRequireV0and(0),
  fTriggerAnalysis(0x0),
  fkVertex(0x0),
  fkVertexAOD(0x0),
  fparMean(0x0),
  fparSigma(0x0),
  fcutSigma(3.),
  fparMinVtxContributors(0x0),
  fparMaxVtxContributors(0x0)
{
  //
  // Named Constructor
  //
  for (Int_t icut=0; icut<5; ++icut){
    fCorrCutMin[icut]=0x0;
    fCorrCutMax[icut]=0x0;
  }
}

//______________________________________________
AliDielectronEventCuts::~AliDielectronEventCuts()
{
  //
  // Default Destructor
  //
  if (fUsedVars) delete fUsedVars;
  if (fTriggerAnalysis) delete fTriggerAnalysis;
}

//______________________________________________
Bool_t AliDielectronEventCuts::IsSelected(TObject* event)
{
  //
  // check the cuts
  //
  
  if(event->IsA() == AliESDEvent::Class()) return IsSelectedESD(event);
  else if(event->IsA() == AliAODEvent::Class()) return IsSelectedAOD(event);
  else return kFALSE;
}
//____________________________________________________________________
Bool_t AliDielectronEventCuts::IsSelectedESD(TObject* event)
{
  //
  // check the cuts
  //
  
  AliESDEvent *ev=dynamic_cast<AliESDEvent*>(event);
  if (!ev) return kFALSE;

  if (fCentMin<fCentMax){
    AliCentrality *centrality=ev->GetCentrality();
    Double_t centralityF=-1;
    if (centrality) centralityF = centrality->GetCentralityPercentile("V0M");
    if (centralityF<fCentMin || centralityF>=fCentMax) return kFALSE;
  }

  fkVertex=0x0;

  switch(fVtxType){
  case kVtxTracks:
  case kVtxTracksOrSPD:
    fkVertex=ev->GetPrimaryVertexTracks();
    break;
  case kVtxSPD:    fkVertex=ev->GetPrimaryVertexSPD(); break;
  case kVtxTPC:    fkVertex=ev->GetPrimaryVertexTPC(); break;
  case kVtxAny:    fkVertex=ev->GetPrimaryVertex(); break;
  }

  if ((fRequireVtx||fVtxZmin<fVtxZmax||fMinVtxContributors>0)&&!fkVertex) return kFALSE;
  

  if (fMinVtxContributors>0){
    Int_t nCtrb = fkVertex->GetNContributors();
    if (nCtrb<fMinVtxContributors){
      if (fVtxType==kVtxTracksOrSPD){
        fkVertex=ev->GetPrimaryVertexSPD();
        nCtrb = fkVertex->GetNContributors();
        if (nCtrb<fMinVtxContributors) return kFALSE;
      } else {
        return kFALSE;
      }
    }
  }

  if (fVtxZmin<fVtxZmax){
    Double_t zvtx=fkVertex->GetZ();
    if (zvtx<fVtxZmin||zvtx>fVtxZmax) return kFALSE;
  }

  if(fRequire13sel){
    if(!fUtils.IsVertexSelected2013pA(ev)) return kFALSE;
    if(fUtils.IsFirstEventInChunk(ev)) return kFALSE;
  }

  if (fRequireV0and){
    if (!fTriggerAnalysis) fTriggerAnalysis=new AliTriggerAnalysis;
    Bool_t v0AND = kFALSE;
    if (fRequireV0and==1){
      Bool_t v0A       = fTriggerAnalysis->IsOfflineTriggerFired(ev, AliTriggerAnalysis::kV0A);
      Bool_t v0C       = fTriggerAnalysis->IsOfflineTriggerFired(ev, AliTriggerAnalysis::kV0C);
      v0AND = v0A && v0C;
    }

    if (fRequireV0and==2){
      Bool_t v0AHW     = (fTriggerAnalysis->V0Trigger(ev, AliTriggerAnalysis::kASide, kTRUE) == AliTriggerAnalysis::kV0BB);
      Bool_t v0CHW     = (fTriggerAnalysis->V0Trigger(ev, AliTriggerAnalysis::kCSide, kTRUE) == AliTriggerAnalysis::kV0BB);
      v0AND = v0AHW && v0CHW;
    }

    if (!v0AND) return kFALSE;
  }

  if (fMultITSTPC){
    const AliESDVertex *vtxESDTPC=ev->GetPrimaryVertexTPC();
    const AliMultiplicity *multESD = ev->GetMultiplicity();
    if ( vtxESDTPC && multESD && vtxESDTPC->GetNContributors() < (-10.+0.25*multESD->GetNumberOfITSClusters(0)) )
      return kFALSE;
  }

  if(fparMean && fparSigma) {
    Double_t nTrks  = ev->GetNumberOfTracks();
    Double_t multV0 = 0.0;
    for(Int_t j=0; j<64; j++) multV0 += ev->GetVZEROData()->GetMultiplicity(j);
    Double_t mV0 = fparMean->Eval(nTrks);
    Double_t sV0 = fparSigma->Eval(nTrks);
    if(multV0 > mV0+fcutSigma*sV0 || multV0 < mV0-fcutSigma*sV0) return kFALSE;
  }

  // cut on the number of vertex contributors using TPC versus global vertex
  if(fparMinVtxContributors && fparMaxVtxContributors) {
    const AliESDVertex *vtxTPC = ev->GetPrimaryVertexTPC();
    const AliESDVertex *vtxGbl = ev->GetPrimaryVertex();
    Double_t nContribTPC = (vtxTPC ? vtxTPC->GetNContributors() : 0);
    Double_t nContribGbl = (vtxGbl ? vtxGbl->GetNContributors() : 0);
    Double_t minCut = fparMinVtxContributors->Eval(nContribGbl);
    Double_t maxCut = fparMaxVtxContributors->Eval(nContribGbl);
    if(nContribTPC > maxCut || nContribTPC < minCut) return kFALSE;
  }

  // run rejection
  Int_t run = ev->GetRunNumber();
  if(fRun.GetNrows()) {
    for(Int_t irun=0; irun<fRun.GetNrows(); irun++) {
      if(fRun(irun)==run) return kFALSE;
    }
  }

  return kTRUE;
}
//______________________________________________
Bool_t AliDielectronEventCuts::IsSelectedAOD(TObject* event)
{
  //
  // check the cuts
  //
  
  AliAODEvent *ev=dynamic_cast<AliAODEvent*>(event);
  if (!ev) return kFALSE;

  //Fill values
  Double_t values[AliDielectronVarManager::kNMaxValues];
  if(fUsedVars->CountBits()) {
    AliDielectronVarManager::SetFillMap(fUsedVars);
    AliDielectronVarManager::Fill(ev,values);

    // correlation cuts
    for(Int_t i=0; i<5; i++) {
      if(fCorrCutMin[i]) {
	Double_t varx = values[fCorrCutMin[i]->GetXaxis()->GetUniqueID()];
	Double_t vary = values[fCorrCutMin[i]->GetYaxis()->GetUniqueID()];
	Double_t min  = ((TF1*)fCorrCutMin[i]->GetListOfFunctions()->At(0))->Eval(varx);
	//      printf("coor cut %d: varx %f -> eval %f > %f \n",i,varx,min,vary);
	if(vary<min) return kFALSE;
      }
      if(fCorrCutMax[i]) {
	Double_t varx = values[fCorrCutMax[i]->GetXaxis()->GetUniqueID()];
	Double_t vary = values[fCorrCutMax[i]->GetYaxis()->GetUniqueID()];
	Double_t max  = ((TF1*)fCorrCutMax[i]->GetListOfFunctions()->At(0))->Eval(varx);
	if(vary>max) return kFALSE;
      }
    }
  }

  // run rejection
  Int_t run = ev->GetRunNumber();
  if(fRun.GetNrows()) {
    for(Int_t irun=0; irun<fRun.GetNrows(); irun++) {
      if(fRun(irun)==run) return kFALSE;
    }
  }

  if (fCentMin<fCentMax){
    AliCentrality *centrality=ev->GetCentrality();
    Double_t centralityF=-1;
    if (centrality) centralityF = centrality->GetCentralityPercentile("V0M");
    if (centralityF<fCentMin || centralityF>=fCentMax) return kFALSE;
  }

  fkVertexAOD=0x0;

  switch(fVtxType){
  case kVtxTracks:         fkVertexAOD=0x0;                       break;
  case kVtxTPC:            fkVertexAOD=AliDielectronVarManager::GetVertex(ev, AliAODVertex::kMainTPC);   break;
  case kVtxSPD:
  case kVtxTracksOrSPD:    fkVertexAOD=ev->GetPrimaryVertexSPD(); break;
  case kVtxAny:            fkVertexAOD=ev->GetPrimaryVertex();    break;
  }

  if ((fRequireVtx||fVtxZmin<fVtxZmax||fMinVtxContributors>0)&&!fkVertexAOD) return kFALSE;
  
  if (fMinVtxContributors>0){
    Int_t nCtrb = fkVertexAOD->GetNContributors();
    if (nCtrb<fMinVtxContributors){
      // if (fVtxType==kVtxTracksOrSPD){
      //   fkVertexAOD=ev->GetVertex(AliAODVertex::kPrimary);
      //   nCtrb = fkVertexAOD->GetNContributors();
      //   if (nCtrb<fMinVtxContributors) return kFALSE;
      //      } else {
      return kFALSE;
      //}
    }
  }
  

  if (fVtxZmin<fVtxZmax){
    Double_t zvtx=fkVertexAOD->GetZ();
    if (zvtx<fVtxZmin||zvtx>fVtxZmax) return kFALSE;
  }

  if(fRequire13sel){
    if(!fUtils.IsVertexSelected2013pA(ev)) return kFALSE;
//     if(fUtils.IsFirstEventInChunk(ev)) return kFALSE;
  }

  /*
  if (fRequireV0and){
    //    if (!fTriggerAnalysis) fTriggerAnalysis=new AliTriggerAnalysis;
    Bool_t v0AND = kFALSE;
    if (fRequireV0and==1){
      Bool_t v0A       = fTriggerAnalysis->IsOfflineTriggerFired(ev, AliTriggerAnalysis::kV0A);
      Bool_t v0A       = header->GetOfflineTrigger(); //TODO
      Bool_t v0C       = fTriggerAnalysis->IsOfflineTriggerFired(ev, AliTriggerAnalysis::kV0C);
      v0AND = v0A && v0C;
    }

    if (fRequireV0and==2){
      Bool_t v0AHW     = (fTriggerAnalysis->V0Trigger(ev, AliTriggerAnalysis::kASide, kTRUE) == AliTriggerAnalysis::kV0BB);
      Bool_t v0CHW     = (fTriggerAnalysis->V0Trigger(ev, AliTriggerAnalysis::kCSide, kTRUE) == AliTriggerAnalysis::kV0BB);
      v0AND = v0AHW && v0CHW;
    }

    if (!v0AND) return kFALSE;
  }
  */
  /*  if (fMultITSTPC){
    const AliESDVertex *vtxESDTPC=ev->GetPrimaryVertexTPC();
    const AliMultiplicity *multESD = ev->GetMultiplicity();
    if ( vtxESDTPC && multESD && vtxESDTPC->GetNContributors() < (-10.+0.25*multESD->GetNumberOfITSClusters(0)) )
      return kFALSE;
  }
  */
  
  // correlation cut Ntrks vs. multV0
  if(fparMean && fparSigma) {
    Double_t nTrks  = ev->GetNumberOfTracks();
    Double_t multV0 = 0.0;
    for(Int_t j=0; j<64; j++) multV0 += ev->GetVZEROData()->GetMultiplicity(j);
    Double_t mV0 = fparMean->Eval(nTrks);
    Double_t sV0 = fparSigma->Eval(nTrks);
    if(multV0 > mV0+fcutSigma*sV0 || multV0 < mV0-fcutSigma*sV0) return kFALSE;
  }

  // cut on the number of vertex contributors using TPC versus global vertex
  if(fparMinVtxContributors && fparMaxVtxContributors) {
    const AliAODVertex *vtxTPC = ev->GetVertex(AliAODVertex::kMainTPC);
    const AliAODVertex *vtxGbl = ev->GetPrimaryVertex();
    Double_t nContribTPC = (vtxTPC ? vtxTPC->GetNContributors() : 0);
    Double_t nContribGbl = (vtxGbl ? vtxGbl->GetNContributors() : 0);
    Double_t minCut = fparMinVtxContributors->Eval(nContribGbl);
    Double_t maxCut = fparMaxVtxContributors->Eval(nContribGbl);
    if(nContribTPC > maxCut || nContribTPC < minCut) return kFALSE;
  }

  return kTRUE;
}

//______________________________________________
void AliDielectronEventCuts::SetMinCorrCutFunction(TF1 *fun, UInt_t varx, UInt_t vary)
{
  //
  // add correlation cut using a TF1
  //
  fUsedVars->SetBitNumber(varx,kTRUE);
  fUsedVars->SetBitNumber(vary,kTRUE);
  // store variables
  fun->GetXaxis()->SetUniqueID(varx);
  fun->GetYaxis()->SetUniqueID(vary);

  Int_t i=0;
  for(i=0; i<5; i++) {
    if(fCorrCutMin[i]) continue;
    else {
      TString key=GetName(); key+=Form("Min%d",i);
      // clone temporare histogram since otherwise it will not be streamed to file!
      fCorrCutMin[i] = (TH1*)fun->GetHistogram()->Clone(key.Data());
      fCorrCutMin[i]->GetListOfFunctions()->AddAt(fun,0);
      break;
    }
  }
  //printf("-----> corr cut added to %d %p \n",i,fCorrCutMin[i]);
  //  fCorrCutMin[i]->Print();
  //fCorrCutMin[i]->GetListOfFunctions()->ls();
}

//______________________________________________
void AliDielectronEventCuts::SetMaxCorrCutFunction(TF1 *fun, UInt_t varx, UInt_t vary)
{
  //
  // add correlation cut using a TF1
  //
  fUsedVars->SetBitNumber(varx,kTRUE);
  fUsedVars->SetBitNumber(vary,kTRUE);
  // store variables
  fun->GetXaxis()->SetUniqueID(varx);
  fun->GetYaxis()->SetUniqueID(vary);

  Int_t i=0;
  for(i=0; i<5; i++) {
    if(fCorrCutMax[i]) continue;
    else {
      TString key=GetName(); key+=Form("Max%d",i);
      // clone temporare histogram since otherwise it will not be streamed to file!
      fCorrCutMax[i] = (TH1*)fun->GetHistogram()->Clone(key.Data());
      fCorrCutMax[i]->GetListOfFunctions()->AddAt(fun,0);
      break;
    }
  }
  //printf("-----> corr cut added to %d %p \n",i,fCorrCutMax[i]);
  //  fCorrCutMax[i]->Print();
  //  fCorrCutMax[i]->GetListOfFunctions()->ls();
}

//________________________________________________________________________
void AliDielectronEventCuts::Print(const Option_t* /*option*/) const
{
  //
  // Print cuts and the range
  //
  printf("cut ranges for '%s'\n",GetTitle());
  printf("All Cuts have to be fulfilled\n");

  Int_t iCut=0;
  if(fRequireVtx) {
    printf("Cut %02d: vertex required \n",iCut);                                   iCut++; }
  printf("Cut %02d: vertex type: %s \n", iCut, fgkVtxNames[fVtxType]);             iCut++;
  if(fMinVtxContributors) {
    printf("Cut %02d: vertex contributors >= %d \n", iCut, fMinVtxContributors);   iCut++; }
  if(fVtxZmin<fVtxZmax) {
    printf("Cut %02d: %f < %s < %f\n",   iCut, fVtxZmin, "Zvtx", fVtxZmax);        iCut++; }
  if(fCentMin<fCentMax) {
    printf("Cut %02d: %f < %s < %f\n",   iCut, fCentMin, "V0centrality", fCentMax);iCut++; }
  if(fMultITSTPC) {
    printf("Cut %02d: cut on multiplcity ITS vs. TPC \n", iCut);                   iCut++; }
  if(fparMean&&fparSigma) {
    printf("Cut %02d: multplicity vs. #tracks correlation +-%.1f sigma inclusion \n", iCut, fcutSigma); iCut++; }
  if(fRequire13sel){
    printf("Cut %02d: vertex and event selection for 2013 pPb data taking required \n",iCut);   iCut++; } 
  if(fRequireV0and) {
    printf("Cut %02d: require V0and type: %c \n", iCut, fRequireV0and);            iCut++; }

}

