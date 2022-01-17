/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
// AliAnalysisTaskSE for the reconstruction of heavy flavor
// decays, using the class AliAnalysisVertexingHF (for conversions to AO2D).
//
// Author: F.Prino, prino@to.infn.it
/////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>
#include <TTree.h>

#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliESDEvent.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSEVertexingHFRun3Conversion.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDUtils.h"

#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEVertexingHFRun3Conversion);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEVertexingHFRun3Conversion::AliAnalysisTaskSEVertexingHFRun3Conversion():
AliAnalysisTaskSE(),
fVHF(0),
fMakeReducedCandidates(kTRUE),
fResetTreeAtEachEv(kTRUE),
fD0CandidateTree(0x0),
f3ProngCandidateTree(0x0),
fDstarCandidateTree(0x0),
fCascadeCandidateTree(0x0),
fEventIndex(-1),
fD0track0(-1),
fD0track1(-1),
f3ptrack0(-1),
f3ptrack1(-1),
f3ptrack2(-1),
fDstD0tr0(-1),
fDstD0tr1(-1),
fDstSofPi(-1),
fCasV0ind(-1),
fCasV0tr0(-1),
fCasV0tr1(-1),
fCasBachl(-1),
fListOfCuts(0),
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
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEVertexingHFRun3Conversion::AliAnalysisTaskSEVertexingHFRun3Conversion(const char *name):
AliAnalysisTaskSE(name),
fVHF(0),
fMakeReducedCandidates(kTRUE),
fResetTreeAtEachEv(kTRUE),
fD0CandidateTree(0x0),
f3ProngCandidateTree(0x0),
fDstarCandidateTree(0x0),
fCascadeCandidateTree(0x0),
fEventIndex(-1),
fD0track0(-1),
fD0track1(-1),
f3ptrack0(-1),
f3ptrack1(-1),
f3ptrack2(-1),
fDstD0tr0(-1),
fDstD0tr1(-1),
fDstSofPi(-1),
fCasV0ind(-1),
fCasV0tr0(-1),
fCasV0tr1(-1),
fCasBachl(-1),
fListOfCuts(0),
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
  // Standard constructor

  DefineOutput(1,TList::Class()); // analysis cuts
  DefineOutput(2,TTree::Class());  // TTree D0
  DefineOutput(3,TTree::Class());  // TTree 3prong
  DefineOutput(4,TTree::Class());  // TTree D*
  DefineOutput(5,TTree::Class());  // TTree Cascade
}

//________________________________________________________________________
AliAnalysisTaskSEVertexingHFRun3Conversion::~AliAnalysisTaskSEVertexingHFRun3Conversion()
{
  // Destructor

  if(fListOfCuts) {
    delete fListOfCuts;
    fListOfCuts=NULL;
  }
  delete fD0CandidateTree;
  delete f3ProngCandidateTree;
  delete fDstarCandidateTree;
  delete fCascadeCandidateTree;
}  

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHFRun3Conversion::Init()
{
  // Initialization
  // Instanciates vHF and loads its parameters

  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF::Init() \n");

  if(gROOT->LoadMacro("ConfigVertexingHF.C")) {
    printf("AnalysisTaskSEVertexingHF::Init() \n Using $ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF.C\n");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF.C");
  }

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  if(fMakeReducedCandidates) fVHF->SetMakeReducedRHF(kTRUE);
  fVHF->PrintStatus();


  // write the objects AliRDHFCuts to a list to store in the output

  fListOfCuts = fVHF->FillListOfCuts();

  PostData(1,fListOfCuts);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHFRun3Conversion::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEVertexingHF::UserCreateOutPutData() \n");
  if(!fVHF) {
    printf("AnalysisTaskSEVertexingHF::UserCreateOutPutData() \n ERROR! no fvHF!\n");
    return;
  }
  
  fD0CandidateTree = new TTree("candD0Tree", "Tree of D0->Kpi candidates");
  if(!fResetTreeAtEachEv) fD0CandidateTree->Branch("eventId",&fEventIndex,"eventId/I");
  fD0CandidateTree->Branch("trackId0",&fD0track0,"trackId0/I");
  fD0CandidateTree->Branch("trackId1",&fD0track1,"trackId1/I");

  f3ProngCandidateTree = new TTree("cand3ProngTree", "Tree of c->3prong candidates");
  if(!fResetTreeAtEachEv) f3ProngCandidateTree->Branch("eventId",&fEventIndex,"eventId/I");
  f3ProngCandidateTree->Branch("trackId0",&f3ptrack0,"trackId0/I");
  f3ProngCandidateTree->Branch("trackId1",&f3ptrack1,"trackId1/I");
  f3ProngCandidateTree->Branch("trackId2",&f3ptrack2,"trackId2/I");
  
  fDstarCandidateTree = new TTree("candDstarTree", "Tree of D*->D0pi candidates");
  if(!fResetTreeAtEachEv) fDstarCandidateTree->Branch("eventId",&fEventIndex,"eventId/I");
  fDstarCandidateTree->Branch("trackD0Dau0",&fDstD0tr0,"trackD0Dau0/I");
  fDstarCandidateTree->Branch("trackD0Dau1",&fDstD0tr1,"trackD0Dau1/I");
  fDstarCandidateTree->Branch("trackSoftPi",&fDstSofPi,"trackSoftPi/I");

  fCascadeCandidateTree = new TTree("candCascade", "Tree of Lc->V0+bach candidates");
  if(!fResetTreeAtEachEv) fCascadeCandidateTree->Branch("eventId",&fEventIndex,"eventId/I");
  fCascadeCandidateTree->Branch("v0index",&fCasV0ind,"v0index/I");
  fCascadeCandidateTree->Branch("trackV0Dau0",&fCasV0tr0,"trackV0Dau0/I");
  fCascadeCandidateTree->Branch("trackV0Dau1",&fCasV0tr1,"trackV0Dau1/I");
  fCascadeCandidateTree->Branch("trackBachel",&fCasBachl,"trackBachel/I");

  fVerticesHFTClArr = new TClonesArray("AliAODVertex", 0);
  fVerticesHFTClArr->SetName("VerticesHF");
  
  if(fVHF->GetD0toKpi()) {
    fD0toKpiTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fD0toKpiTClArr->SetName("D0toKpi");
  }

  if(fVHF->GetJPSItoEle()) {
    fJPSItoEleTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fJPSItoEleTClArr->SetName("JPSItoEle");
  }

  if(fVHF->Get3Prong()) {
    fCharm3ProngTClArr = new TClonesArray("AliAODRecoDecayHF3Prong", 0);
    fCharm3ProngTClArr->SetName("Charm3Prong");
  }

  if(fVHF->Get4Prong()) {
    fCharm4ProngTClArr = new TClonesArray("AliAODRecoDecayHF4Prong", 0);
    fCharm4ProngTClArr->SetName("Charm4Prong");
  }

  if(fVHF->GetDstar()) {
    fDstarTClArr = new TClonesArray("AliAODRecoCascadeHF", 0);
    fDstarTClArr->SetName("Dstar");
  }

  if(fVHF->GetCascades()){
    fCascadesTClArr = new TClonesArray("AliAODRecoCascadeHF", 0);
    fCascadesTClArr->SetName("CascadesHF");
  }

  if(fVHF->GetLikeSign()) {                      
    fLikeSign2ProngTClArr = new TClonesArray("AliAODRecoDecayHF2Prong", 0);
    fLikeSign2ProngTClArr->SetName("LikeSign2Prong");
  }

  if(fVHF->GetLikeSign() && fVHF->Get3Prong()) {                      
    fLikeSign3ProngTClArr = new TClonesArray("AliAODRecoDecayHF3Prong", 0);
    fLikeSign3ProngTClArr->SetName("LikeSign3Prong");
  }


  
  //--------------------------------------------------------------  

  PostData(2,fD0CandidateTree);
  PostData(3,f3ProngCandidateTree);
  PostData(4,fDstarCandidateTree);
  PostData(5,fCascadeCandidateTree);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHFRun3Conversion::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor vertexing
  
  AliVEvent *event = dynamic_cast<AliVEvent*> (InputEvent());
  // In case there is an AOD handler writing a standard AOD, use the AOD 
  // event in memory rather than the input (ESD) event. (A.G. 27/04/09)
  if (AODEvent() && IsStandardAOD()) event = dynamic_cast<AliVEvent*> (AODEvent());


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  fVHF->SetPidResponse(pidResp);

  // heavy flavor vertexing
  fVHF->FindCandidates(event,
                       fVerticesHFTClArr,
                       fD0toKpiTClArr,
                       fJPSItoEleTClArr,
                       fCharm3ProngTClArr,
                       fCharm4ProngTClArr,
                       fDstarTClArr,
                       fCascadesTClArr,
                       fLikeSign2ProngTClArr,
                       fLikeSign3ProngTClArr);
  
  if(fResetTreeAtEachEv){
    fD0CandidateTree->Reset();
    f3ProngCandidateTree->Reset();
    fDstarCandidateTree->Reset();
    fCascadeCandidateTree->Reset();
  }
  
  fEventIndex=event->GetEventNumberInFile();
  Int_t nD0=fD0toKpiTClArr->GetEntriesFast();
  for(Int_t iD0=0; iD0<nD0; iD0++){
    AliAODRecoDecayHF2Prong* dcand=(AliAODRecoDecayHF2Prong*)fD0toKpiTClArr->At(iD0);
    fD0track0=dcand->GetProngID(0);
    fD0track1=dcand->GetProngID(1);
    fD0CandidateTree->Fill();
    //    printf("Event %d cand %d  tracks %d %d\n",fEventIndex,iD0,fD0track0,fD0track1);
  }
  
  Int_t n3p=fCharm3ProngTClArr->GetEntriesFast();
  for(Int_t i3p=0; i3p<n3p; i3p++){
    AliAODRecoDecayHF3Prong* dcand=(AliAODRecoDecayHF3Prong*)fCharm3ProngTClArr->At(i3p);
    f3ptrack0=dcand->GetProngID(0);
    f3ptrack1=dcand->GetProngID(1);
    f3ptrack2=dcand->GetProngID(2);
    f3ProngCandidateTree->Fill();
  }

  Int_t nDst=fDstarTClArr->GetEntriesFast();
  for(Int_t iDst=0; iDst<nDst; iDst++){
    AliAODRecoCascadeHF* dcand=(AliAODRecoCascadeHF*)fDstarTClArr->At(iDst);
    fDstSofPi=dcand->GetProngID(0);
    Int_t iD0=dcand->GetProngID(1);
    AliAODRecoDecayHF2Prong* d0cand=(AliAODRecoDecayHF2Prong*)fD0toKpiTClArr->At(iD0);
    fDstD0tr0=d0cand->GetProngID(0);
    fDstD0tr1=d0cand->GetProngID(1);
    //    printf("Event %d D0 from D* id %d  dau %d %d\n",fEventIndex,iD0,fDstD0tr0,fDstD0tr1);
    fDstarCandidateTree->Fill();
  }

  Int_t nCas=fCascadesTClArr->GetEntriesFast();
  for(Int_t iCas=0; iCas<nCas; iCas++){
    AliAODRecoCascadeHF* dcand=(AliAODRecoCascadeHF*)fCascadesTClArr->At(iCas);
    fCasBachl=dcand->GetProngID(0);
    fCasV0ind=dcand->GetProngID(1);
    if(event->IsA()->InheritsFrom("AliAODEvent")){
      AliAODv0 *v0=((AliAODEvent*)event)->GetV0(fCasV0ind);
      fCasV0tr0=v0->GetPosID();
      fCasV0tr1=v0->GetNegID();
    }else{
      AliESDv0* v0=((AliESDEvent*)event)->GetV0(fCasV0ind);
      fCasV0tr0=v0->GetPindex();
      fCasV0tr1=v0->GetNindex();
    }
    fCascadeCandidateTree->Fill();
  }
  

  PostData(2,fD0CandidateTree);
  PostData(3,f3ProngCandidateTree);
  PostData(4,fDstarCandidateTree);
  PostData(5,fCascadeCandidateTree);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEVertexingHFRun3Conversion::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEVertexingHFRun3Conversion: Terminate() \n");
}
