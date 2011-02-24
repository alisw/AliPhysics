/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//*************************************************************************
// Class AliAnalysisTaskMEVertexingHF
// AliAnalysisTaskME for event mixing, building the background for 
// heavy-flavour decay candidates
// Author: R.Romita, r.romita@gsi.de
//*************************************************************************



#include "TH1F.h"
#include "TObjArray.h"
#include "TList.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"

#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisVertexingHF.h"
#include "AliMixedEvent.h"
#include "AliAnalysisTaskMEVertexingHF.h"
#include "AliAnalysisManager.h"
#include "AliMultiEventInputHandler.h"

ClassImp(AliAnalysisTaskMEVertexingHF)

//________________________________________________________________________
AliAnalysisTaskMEVertexingHF::AliAnalysisTaskMEVertexingHF(const char *name) : 
AliAnalysisTaskME(name), 
fvHF(0), 
fMixedEvent(),
fVerticesHFTClArr(0),
fD0toKpiTClArr(0), 
fJPSItoEleTClArr(0),
fCharm3ProngTClArr(0),
fCharm4ProngTClArr(0),
fDstarTClArr(0),
fCascadesTClArr(0),
fLikeSign2ProngTClArr(0),
fLikeSign3ProngTClArr(0)
{
  // Constructor
}
//________________________________________________________________________
void AliAnalysisTaskMEVertexingHF::Init()
{
 // Initialization
 // Instanciates vHF and loads its parameters
 // Some parameters are changed
 
  if(gROOT->LoadMacro("ConfigVertexingHF.C")) {
    printf("AnalysisTaskMEVertexingHF::Init() \n Using $ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C\n");
    gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");
  }
  fvHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  fvHF->SetMixEventOn();
  fvHF->SetInputAOD();
  fvHF->PrintStatus();
  if(fvHF->GetLikeSign()) {
    printf("WARNING: fLikeSign will be switched off!");
    fvHF->SetLikeSignOff();
  }
  if(fvHF->GetRecoPrimVtxSkippingTrks() || fvHF->GetRmTrksFromPrimVtx()){
    fvHF->UnsetRecoPrimVtxSkippingTrks();
    printf("WARNING: if on, fRecoPrimVtxSkippingTrks and fRmTrksFromPrimVtx  will be switched off!\n");
  }
  
  AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile("AliAOD.VertexingHF.root");

  return;
}
//________________________________________________________________________
void AliAnalysisTaskMEVertexingHF::UserCreateOutputObjects()
{  
// Create the output container


  if (!AODEvent()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }
  
  if(!fvHF) {
    printf("AnalysisTaskMEVertexingHF::UserCreateOutPutData() \n ERROR! no fvHF!\n");
    return;
  }
  fVerticesHFTClArr = new TClonesArray("AliAODVertex", 0);
  fVerticesHFTClArr->SetName("VerticesHF");
  AddAODBranch("TClonesArray", &fVerticesHFTClArr);
  if(fvHF->GetD0toKpi()) {
    fD0toKpiTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fD0toKpiTClArr->SetName("D0toKpi");
    AddAODBranch("TClonesArray", &fD0toKpiTClArr);
  }
  if(fvHF->GetJPSItoEle()) {
    fJPSItoEleTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fJPSItoEleTClArr->SetName("JPSItoEle");
    AddAODBranch("TClonesArray", &fJPSItoEleTClArr);
  }
  if(fvHF->Get3Prong()) {
    fCharm3ProngTClArr = new TClonesArray("AliAODRecoDecayHF3Prong", 0);
    fCharm3ProngTClArr->SetName("Charm3Prong");
    AddAODBranch("TClonesArray", &fCharm3ProngTClArr);
  }
  if(fvHF->Get4Prong()) {
    fCharm4ProngTClArr = new TClonesArray("AliAODRecoDecayHF4Prong", 0);
    fCharm4ProngTClArr->SetName("Charm4Prong");
    AddAODBranch("TClonesArray", &fCharm4ProngTClArr);
  }
  
  if(fvHF->GetDstar()) {
    fDstarTClArr = new TClonesArray("AliAODRecoCascadeHF", 0);
    fDstarTClArr->SetName("Dstar");
    AddAODBranch("TClonesArray", &fDstarTClArr);
  }
  
  if(fvHF->GetCascades()){
    fCascadesTClArr = new TClonesArray("AliAODRecoCascadeHF", 0);
    fCascadesTClArr->SetName("CascadesHF");
    AddAODBranch("TClonesArray", &fCascadesTClArr);
  }

  if(fvHF->GetLikeSign()) {
    fLikeSign2ProngTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fLikeSign2ProngTClArr->SetName("LikeSign2Prong");
    AddAODBranch("TClonesArray", &fLikeSign2ProngTClArr);
  }
 
  if(fvHF->GetLikeSign() && fvHF->Get3Prong()) {
    fLikeSign3ProngTClArr = new TClonesArray("AliAODRecoDecayHF3Prong", 0);
    fLikeSign3ProngTClArr->SetName("LikeSign3Prong");
    AddAODBranch("TClonesArray", &fLikeSign3ProngTClArr);
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskMEVertexingHF::UserExec(Option_t *) 
{
  // Execute analysis for current event:
  // first build the mixed event, compute the new primary vtx 
  // then heavy flavor vertexing 

  Int_t nev = fInputHandler->GetBufferSize();
  fMixedEvent = new AliMixedEvent();
  fMixedEvent->Reset();
  TString primTitle;
  TString primTitleFirst;
  AliAODVertex *vtxCopy=0;

  TObjArray *vertices=new TObjArray(nev);
  for (Int_t iev = 0; iev < nev; iev++) {
    AliAODEvent *evt = (AliAODEvent*)GetEvent(iev);
    if(!evt) {delete vertices;return;}
    AliAODVertex *evtVtx=(AliAODVertex*)evt->GetPrimaryVertex();
    if(!evtVtx) {delete vertices;return;}
    primTitle = evtVtx->GetTitle();
    Int_t nContrib=evtVtx->GetNContributors();
    if(!primTitle.Contains("VertexerTracks") || nContrib<=0) {
      delete vertices;
      return;
    }

    vtxCopy=new AliAODVertex(*evtVtx);
    primTitleFirst=evtVtx->GetTitle();
        

    fMixedEvent->AddEvent(evt);

    vertices->AddLast(vtxCopy);
  }


  fMixedEvent->Init();
  Double_t vtxPos[3]={0.,0.,0.},vtxSigma[3]={0.,0.,0.};
  Int_t nContributors[1]={0};
  Double_t chi2=0;
  Bool_t primaryOk=fMixedEvent->ComputeVtx(vertices,vtxPos,vtxSigma,nContributors);
  if(!primaryOk) {
    delete vertices;
    delete vtxCopy;
    vtxCopy=NULL;
    return;
  }
  Int_t contribCopy=nContributors[0];
  Double_t vtxCov[6]={vtxSigma[0]*vtxSigma[0],0,vtxSigma[1]*vtxSigma[1],0,0,vtxSigma[2]*vtxSigma[2]};
  AliVVertex* newVertex=new AliESDVertex(vtxPos,vtxCov,chi2,contribCopy);
  newVertex->SetTitle(primTitleFirst.Data());
  fMixedEvent->SetPrimaryVertex(newVertex);

  delete vertices;
  delete vtxCopy;
  vtxCopy=NULL;

  fvHF->FindCandidates(fMixedEvent,
		       fVerticesHFTClArr,
		       fD0toKpiTClArr,
		       fJPSItoEleTClArr,
		       fCharm3ProngTClArr,
		       fCharm4ProngTClArr,
		       fDstarTClArr,
		       fCascadesTClArr, 
		       fLikeSign2ProngTClArr,
		       fLikeSign3ProngTClArr);

  delete newVertex;
  return;
}      

//________________________________________________________________________
void AliAnalysisTaskMEVertexingHF::Terminate(Option_t *) 
{
  // Terminate analysis
}
