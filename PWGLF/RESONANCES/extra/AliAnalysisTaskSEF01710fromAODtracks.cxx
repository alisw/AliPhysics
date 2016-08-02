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
 * appeuear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//               F01710 analysis code
//
//  Input: AOD
//  Output: TTree or THnSparse
//
//-------------------------------------------------------------------------
//
//                 Authors: Y.S Watanabe(a)
//  (a) CNS, the University of Tokyo
//  Contatcs: wyosuke@cns.s.u-tokyo.ac.jp
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <THnSparse.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliAODRecoDecay.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEF01710fromAODtracks.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNeutralTrackParam.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"
#include "AliVertexerTracks.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEF01710fromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEF01710fromAODtracks::AliAnalysisTaskSEF01710fromAODtracks() : 
  AliAnalysisTaskSE(),
  fAnalysisType(0),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fCEvents(0),
  fCEventsNorm(0),
  fHTrigger(0),
  fHCentrality(0),
  fHVtxZ(0),
  fHMultiplicityV0A(0),
  fHMultiplicityV0C(0),
  fIsEventSelected(kFALSE),
  fWhyRejection(0),
  fWriteVariableTree(kFALSE),
  fVariablesTree(0),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0.),
  fBzkG(0),
  fCentrality(0),
  fTriggerCheck(0),
  fMultiplicityVZEROA(0),
  fMultiplicityVZEROC(0),
  fPIDResponse(0),
  fHistoF01710Mass(0),
  fHistoF01710MassMix(0),
  fHistoF01710MassMCGen(0),
  fHistoF01710MassMCS(0),
  fHistoF21525MassMCGen(0),
  fHistoF21525MassMCS(0),
  fHistoA21320MassMCGen(0),
  fHistoA21320MassMCS(0),
  fHistoF21270MassMCGen(0),
  fHistoF21270MassMCS(0),
  fHistoF01710ChargedMass(0),
  fHistoF01710ChargedLikeMass(0),
  fHistoF01710ChargedMassMix(0),
  fHistoF01710ChargedLikeMassMix(0),
  fHistoF01710ChargedMassMCGen(0),
  fHistoF01710ChargedMassMCS(0),
  fHistoF21525ChargedMassMCGen(0),
  fHistoF21525ChargedMassMCS(0),
  fHistoA21320ChargedMassMCGen(0),
  fHistoA21320ChargedMassMCS(0),
  fHistoF21270ChargedMassMCGen(0),
  fHistoF21270ChargedMassMCS(0),
  fHistonEvtvsRunNumber(0),
  fHistonKaonvsRunNumber(0),
  fHistonK0vsRunNumber(0),
  fDoEventMixing(0),
  fNumberOfEventsForMixing		(5),
  fNzVtxBins					(0), 
  fNCentBins					(0),
  fNOfPools(1),
  fEventBuffer(0x0),
  fEventInfo(new TObjString("")),
  fK0Tracks(0x0),
  fK0CutVarsArray(0x0),
  fKaonTracks(0x0),
  fKaonCutVarsArray(0x0)
{
  //
  // Default Constructor. 
  //
}

//___________________________________________________________________________
AliAnalysisTaskSEF01710fromAODtracks::AliAnalysisTaskSEF01710fromAODtracks(const Char_t* name,
    Bool_t writeVariableTree) :
  AliAnalysisTaskSE(name),
  fAnalysisType(0),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fCEvents(0),
  fCEventsNorm(0),
  fHTrigger(0),
  fHCentrality(0),
  fHVtxZ(0),
  fHMultiplicityV0A(0),
  fHMultiplicityV0C(0),
  fIsEventSelected(kFALSE),
  fWhyRejection(0),
  fWriteVariableTree(writeVariableTree),
  fVariablesTree(0),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0.),
  fBzkG(0),
  fCentrality(0),
  fTriggerCheck(0),
  fMultiplicityVZEROA(0),
  fMultiplicityVZEROC(0),
  fPIDResponse(0),
  fHistoF01710Mass(0),
  fHistoF01710MassMix(0),
  fHistoF01710MassMCGen(0),
  fHistoF01710MassMCS(0),
  fHistoF21525MassMCGen(0),
  fHistoF21525MassMCS(0),
  fHistoA21320MassMCGen(0),
  fHistoA21320MassMCS(0),
  fHistoF21270MassMCGen(0),
  fHistoF21270MassMCS(0),
  fHistoF01710ChargedMass(0),
  fHistoF01710ChargedLikeMass(0),
  fHistoF01710ChargedMassMix(0),
  fHistoF01710ChargedLikeMassMix(0),
  fHistoF01710ChargedMassMCGen(0),
  fHistoF01710ChargedMassMCS(0),
  fHistoF21525ChargedMassMCGen(0),
  fHistoF21525ChargedMassMCS(0),
  fHistoA21320ChargedMassMCGen(0),
  fHistoA21320ChargedMassMCS(0),
  fHistoF21270ChargedMassMCGen(0),
  fHistoF21270ChargedMassMCS(0),
  fHistonEvtvsRunNumber(0),
  fHistonKaonvsRunNumber(0),
  fHistonK0vsRunNumber(0),
  fDoEventMixing(0),
  fNumberOfEventsForMixing		(5),
  fNzVtxBins					(0), 
  fNCentBins					(0),
  fNOfPools(1),
  fEventBuffer(0x0),
  fEventInfo(new TObjString("")),
  fK0Tracks(0x0),
  fK0CutVarsArray(0x0),
  fKaonTracks(0x0),
  fKaonCutVarsArray(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEF01710fromAODtracks","Calling Constructor");

  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
  DefineOutput(3,TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskSEF01710fromAODtracks::~AliAnalysisTaskSEF01710fromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEF01710fromAODtracks","Calling Destructor");


  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }

  if (fVariablesTree) {
    delete fVariablesTree;
    fVariablesTree = 0;
  }

  if(fCandidateVariables){
    delete fCandidateVariables;
    fCandidateVariables = 0;
  }

  if(fK0Tracks) fK0Tracks->Delete();
  delete fK0Tracks;
  if(fKaonTracks) fKaonTracks->Delete();
  delete fKaonTracks;

  if(fK0CutVarsArray) fK0CutVarsArray->Delete();
  delete fK0CutVarsArray;
  if(fKaonCutVarsArray) fKaonCutVarsArray->Delete();
  delete fKaonCutVarsArray;

  if(fEventBuffer){
    for(Int_t i=0; i<fNOfPools; i++) delete fEventBuffer[i];
    delete fEventBuffer;
  }
  delete fEventInfo;

}

//_________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::Init() {
  //
  // Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  return;
}

//_________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::UserExec(Option_t *)
{
  //
  // UserExec code
  //

  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);


  fCEvents->Fill(1);
  //------------------------------------------------
  // First check if the event has proper vertex and B
  //------------------------------------------------
  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  AliKFParticle::SetField(fBzkG);
  if (TMath::Abs(fBzkG)<0.001) {
    delete fV1;
    return;
  }
  fCEvents->Fill(2);

  fIsEventSelected = IsEventSelected(aodEvent); 
  if(fWhyRejection==0) fCEventsNorm->Fill(0);
  if(fIsEventSelected){
    fCEventsNorm->Fill(3);
  }else{
    if(fWhyRejection==0) fCEventsNorm->Fill(1);
    if(fWhyRejection==6) {
      fCEventsNorm->Fill(2);
      fCEventsNorm->Fill(3);
    }
  }

  //------------------------------------------------
  // MC analysis setting
  //------------------------------------------------
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;
  if (fUseMCInfo) {
    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fCEvents->Fill(6); // in case of MC events

    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: MC header branch not found!\n");
      return;
    }
    fCEvents->Fill(7); // in case of MC events

    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > 10.) {
      AliDebug(2,Form("Event rejected: abs(zVtxMC)=%f > 10.=%f",zMCVertex,10.));
      return;
    } else {
      fCEvents->Fill(17); // in case of MC events
    }
    if ((TMath::Abs(zMCVertex) < 10.) && (IsEventAcceptedPhysicsSelection(aodEvent)) && (IsEventAcceptedTrigger(aodEvent))) {
      MakeMCAnalysis(mcArray);
    }
  }

  //------------------------------------------------
  // Event selection 
  //------------------------------------------------
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) return;

  Double_t pos[3],cov[6];
  fVtx1->GetXYZ(pos);
  fVtx1->GetCovarianceMatrix(cov);
  fV1 = new AliESDVertex(pos,cov,100.,100,fVtx1->GetName());
  fVtxZ = pos[2];

  Bool_t fIsTriggerOK = IsEventAcceptedTrigger(aodEvent);
  if(fIsTriggerOK) fCEvents->Fill(3);
  if(!fIsEventSelected) {
    delete fV1;
    return;
  }
  fCEvents->Fill(4);

  fIsMB=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  fIsSemi=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  fIsCent=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral); 
  fIsINT7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);  
  fIsEMC7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);   
  fTriggerCheck = fIsMB+2*fIsSemi+4*fIsCent+8*fIsINT7+16*fIsEMC7;
  if(fIsMB) fHTrigger->Fill(1);
  if(fIsSemi) fHTrigger->Fill(2);
  if(fIsCent) fHTrigger->Fill(3);
  if(fIsINT7) fHTrigger->Fill(4);
  if(fIsEMC7) fHTrigger->Fill(5);
  if(fIsMB|fIsSemi|fIsCent) fHTrigger->Fill(7);
  if(fIsINT7|fIsEMC7) fHTrigger->Fill(8);
  if(fIsMB&fIsSemi) fHTrigger->Fill(10);
  if(fIsMB&fIsCent) fHTrigger->Fill(11);
  if(fIsINT7&fIsEMC7) fHTrigger->Fill(12);

  AliCentrality *cent = aodEvent->GetCentrality();
  //fCentrality = cent->GetCentralityPercentile("V0M");
  fCentrality = 1.;
  fHCentrality->Fill(fCentrality);
  fHVtxZ->Fill(fVtxZ);
  fCEventsNorm->Fill(4);

  AliAODVZERO* aodV0 = aodEvent->GetVZEROData();
  fMultiplicityVZEROA = aodV0->GetMTotV0A();
  fMultiplicityVZEROC = aodV0->GetMTotV0C();
  fHMultiplicityV0A->Fill(fMultiplicityVZEROA);
  fHMultiplicityV0C->Fill(fMultiplicityVZEROC);

  fCEvents->Fill(5);

  Int_t runnumber_offset = 0;
  Int_t runnumber = aodEvent->GetRunNumber();
  if(runnumber<=131000&&runnumber>=114000){
    runnumber_offset = 114000;//lhc10bcde
  }else if(runnumber<=196000&&runnumber>=195000){
    runnumber_offset = 195000;//lhc13bc
  }else if(runnumber<=170593&&runnumber>=167902){
    runnumber_offset = 167902;//lhc11h
  }
  fHistonEvtvsRunNumber->Fill(runnumber-runnumber_offset,1.);

  //------------------------------------------------
  // Main analysis done in this function
  //------------------------------------------------
  MakeAnalysis(aodEvent, mcArray);

  PostData(1,fOutput);
  PostData(2,fVariablesTree);
  PostData(3,fOutputAll);

  fIsEventSelected=kFALSE;

  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSEF01710fromAODtracks::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    AliError("fOutput not available");
    return;
  }

  fOutputAll = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputAll) {     
    AliError("fOutputAll not available");
    return;
  }

  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::UserCreateOutputObjects() 
{ 
  //
  // UserCreateOutputObjects
  //

  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

  //------------------------------------------------
  // output object setting
  //------------------------------------------------
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");
  DefineGeneralHistograms(); // define general histograms
  PostData(1,fOutput);

  DefineTreeVariables();
  PostData(2,fVariablesTree);

  fOutputAll = new TList();
  fOutputAll->SetOwner();
  fOutputAll->SetName("anahisto");
  DefineAnalysisHistograms(); // define general histograms
  PostData(3,fOutputAll);

  if(fDoEventMixing){
    fK0Tracks = new TObjArray();
    fK0Tracks->SetOwner();
    fK0CutVarsArray = new TObjArray();
    fK0CutVarsArray->SetOwner();
    fKaonTracks = new TObjArray();
    fKaonTracks->SetOwner();
    fKaonCutVarsArray = new TObjArray();
    fKaonCutVarsArray->SetOwner();

    fNOfPools=fNCentBins*fNzVtxBins;
    fEventBuffer = new TTree*[fNOfPools];
    for(Int_t i=0; i<fNOfPools; i++){
      fEventBuffer[i]=new TTree(Form("EventBuffer_%d",i), "Temporary buffer for event mixing");
      fEventBuffer[i]->Branch("zVertex", &fVtxZ);
      fEventBuffer[i]->Branch("centrality", &fCentrality);
      fEventBuffer[i]->Branch("eventInfo", "TObjString",&fEventInfo);
      fEventBuffer[i]->Branch("k0array", "TObjArray", &fK0Tracks);
      fEventBuffer[i]->Branch("k0varsarray", "TObjArray", &fK0CutVarsArray);
      fEventBuffer[i]->Branch("kaarray", "TObjArray", &fKaonTracks);
      fEventBuffer[i]->Branch("kavarsarray", "TObjArray", &fKaonCutVarsArray);
    }
  }


  return;
}

//________________________________________________________________________
  void AliAnalysisTaskSEF01710fromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  // Main analysis part called from UserExec
  //
  if(fDoEventMixing){
    if(fK0Tracks) fK0Tracks->Delete();
    if(fK0CutVarsArray) fK0CutVarsArray->Delete();
    if(fKaonTracks) fKaonTracks->Delete();
    if(fKaonCutVarsArray) fKaonCutVarsArray->Delete();
  }

  //------------------------------------------------
  // Select good track before hand to save time
  //------------------------------------------------
  Int_t nTracks= aodEvent->GetNumberOfTracks();
  Int_t  seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags,mcArray);

  Int_t nV0s= aodEvent->GetNumberOfV0s();
  Bool_t  seleV0Flags[nV0s];
  Int_t     nSeleV0=0;
  SelectV0(aodEvent,nV0s,nSeleV0,seleV0Flags,mcArray);

  //------------------------------------------------
  // v0 loop 
  //------------------------------------------------
  for (Int_t iv01 = 0; iv01<nV0s-1; iv01++) {
    if(!seleV0Flags[iv01]) continue;
    AliAODv0 *v01 = aodEvent->GetV0(iv01);
    if(!v01) continue;

    AliAODTrack *cptrack1 =  (AliAODTrack*)(v01->GetDaughter(0));
    AliAODTrack *cntrack1 =  (AliAODTrack*)(v01->GetDaughter(1));
    Int_t cpid1 = cptrack1->GetID();
    Int_t cnid1 = cntrack1->GetID();
    //------------------------------------------------
    // v0 loop 
    //------------------------------------------------
    for (Int_t iv02 = iv01+1; iv02<nV0s; iv02++) {
      if(!seleV0Flags[iv02]) continue;
      AliAODv0 *v02 = aodEvent->GetV0(iv02);
      if(!v02) continue;

      AliAODTrack *cptrack2 =  (AliAODTrack*)(v02->GetDaughter(0));
      AliAODTrack *cntrack2 =  (AliAODTrack*)(v02->GetDaughter(1));
      Int_t cpid2 = cptrack2->GetID();
      Int_t cnid2 = cntrack2->GetID();

      if(cpid1==cpid2) continue;
      if(cnid1==cnid2) continue;

      FillROOTObjects(aodEvent,v01,v02,mcArray);
    }
  }

  for (Int_t itrk1 = 0; itrk1<nTracks-1; itrk1++) {
    if(seleTrkFlags[itrk1]!=1) continue;
    AliAODTrack *trk1 = (AliAODTrack*)aodEvent->GetTrack(itrk1);

    for (Int_t itrk2 = itrk1+1; itrk2<nTracks; itrk2++) {
      if(seleTrkFlags[itrk2]!=1) continue;
      AliAODTrack *trk2 = (AliAODTrack*)aodEvent->GetTrack(itrk2);

      Int_t id1= trk1->GetID();
      Int_t id2= trk2->GetID();
      if(id1==id2) continue;

      FillROOTObjects_ChargedKaon(aodEvent,trk1,trk2,mcArray);

    }
  }

  if(fDoEventMixing){
    fEventInfo->SetString(Form("Ev%d_esd%d_V%d_K%d",AliAnalysisManager::GetAnalysisManager()->GetNcalls(),((AliAODHeader*)aodEvent->GetHeader())->GetEventNumberESDFile(),fK0Tracks->GetEntries(),fKaonTracks->GetEntries()));
    Int_t ind=GetPoolIndex(fVtxZ,fCentrality);
    if(ind>=0 && ind<fNOfPools){
      if(fEventBuffer[ind]->GetEntries() >= fNumberOfEventsForMixing){
        DoEventMixingWithPools(ind);
        if(fEventBuffer[ind]->GetEntries() >= 1000*fNumberOfEventsForMixing){
          ResetPool(ind);
        }
      }
      fEventBuffer[ind]->Fill();
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::FillROOTObjects(AliAODEvent *aod,AliAODv0 *v01, AliAODv0 *v02, TClonesArray *mcArray) 
{
  //
  // Fill histogram or Tree depending on fWriteVariableTree flag
  //
  Double_t mk0s_pdg = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  TLorentzVector vf0,vk0s1,vk0s2;
  vk0s1.SetXYZM(v01->Px(),v01->Py(),v01->Pz(),mk0s_pdg);
  vk0s2.SetXYZM(v02->Px(),v02->Py(),v02->Pz(),mk0s_pdg);
  vf0 = vk0s1+vk0s2;

  Double_t posVtx[3];
  fVtx1->GetXYZ(posVtx);

  for(Int_t i=0;i<26;i++)
    fCandidateVariables[i] = -9999.;

  AliAODMCParticle *mcf01710 = 0;
  AliAODMCParticle *mck01 = 0;
  AliAODMCParticle *mck02 = 0;
  Int_t mclabf01710 = 0;
  Int_t mcpdgk01_array[100];
  Int_t mcpdgk02_array[100];
  Int_t mclabelk01_array[100];
  Int_t mclabelk02_array[100];
  Int_t mcngen_k01=-9999;
  Int_t mcngen_k02=-9999;
  if(fUseMCInfo){
    mclabf01710 =  MatchToMC(v01,v02,mcArray,mcpdgk01_array, mcpdgk02_array,mclabelk01_array,mclabelk02_array,mcngen_k01,mcngen_k02);
    if(mclabf01710>-1){
      mcf01710 = (AliAODMCParticle*) mcArray->At(mclabf01710);
    }
  }

  fCandidateVariables[ 0] = vf0.M();
  fCandidateVariables[ 1] = vf0.Px();
  fCandidateVariables[ 2] = vf0.Py();
  fCandidateVariables[ 3] = vf0.Pz();
  fCandidateVariables[ 4] = v01->MassK0Short();
  fCandidateVariables[ 5] = v01->Px();
  fCandidateVariables[ 6] = v01->Py();
  fCandidateVariables[ 7] = v01->Pz();
  fCandidateVariables[ 8] = v02->MassK0Short();
  fCandidateVariables[ 9] = v02->Px();
  fCandidateVariables[10] = v02->Py();
  fCandidateVariables[11] = v02->Pz();

  fCandidateVariables[12] = v01->CosPointingAngle(posVtx);
  fCandidateVariables[13] = v02->CosPointingAngle(posVtx);

  fCandidateVariables[14] = (Float_t) aod->GetNumberOfTracks();
  fCandidateVariables[15] = 0;
  fCandidateVariables[16] = fMultiplicityVZEROA;
  fCandidateVariables[17] = fMultiplicityVZEROC;
  if(fUseMCInfo){
    if(mcf01710){
      fCandidateVariables[18] = mcf01710->GetPdgCode();
      fCandidateVariables[19] = mcf01710->Px();
      fCandidateVariables[20] = mcf01710->Py();
      fCandidateVariables[21] = mcf01710->Pz();

      fCandidateVariables[22] = mclabelk01_array[0];
      fCandidateVariables[23] = mclabelk02_array[0];
      fCandidateVariables[24] = mclabf01710;
    }
  }
  fCandidateVariables[25] = fCentrality;
  fCandidateVariables[26] = 0;

  if(fWriteVariableTree)
    fVariablesTree->Fill();

  Double_t cont[4];
  cont[0] = vf0.M();
  cont[1] = vf0.Pt();
  if(v01->Pt()>v02->Pt()){
    cont[2] = v02->Pt();
  }else{
    cont[2] = v01->Pt();
  }
  cont[3] = fCentrality;

  fHistoF01710Mass->Fill(cont);

  if(fUseMCInfo){
    if(mcf01710){
      Int_t pdgcode = mcf01710->GetPdgCode();
      if(abs(pdgcode)==10331){
        fHistoF01710MassMCS->Fill(cont);
      }
      if(abs(pdgcode)==335){
        fHistoF21525MassMCS->Fill(cont);
      }
      if(abs(pdgcode)==115){
        fHistoA21320MassMCS->Fill(cont);
      }
      if(abs(pdgcode)==225){
        fHistoF21270MassMCS->Fill(cont);
      }
    }
  }


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::FillROOTObjects_ChargedKaon(AliAODEvent *aod,AliAODTrack *trk1, AliAODTrack *trk2, TClonesArray *mcArray) 
{
  //
  // Fill histogram or Tree depending on fWriteVariableTree flag
  //
  for(Int_t i=0;i<27;i++)
    fCandidateVariables[i] = -9999.;

  TLorentzVector vka1,vka2;
  vka1.SetXYZM(trk1->Px(),trk1->Py(),trk1->Pz(),0.493677);
  vka2.SetXYZM(trk2->Px(),trk2->Py(),trk2->Pz(),0.493677);
  TLorentzVector vf0 = vka1+vka2;

  AliAODMCParticle *mcf01710 = 0;
  AliAODMCParticle *mcka1 = 0;
  AliAODMCParticle *mcka2 = 0;
  Int_t mclabf01710 = 0;
  Int_t mcpdgka1_array[100];
  Int_t mcpdgka2_array[100];
  Int_t mclabelka1_array[100];
  Int_t mclabelka2_array[100];
  Int_t mcngen_ka1=-9999;
  Int_t mcngen_ka2=-9999;
  if(fUseMCInfo){
    mclabf01710 =  MatchToMC_ChargedKaon(trk1,trk2,mcArray,mcpdgka1_array, mcpdgka2_array,mclabelka1_array,mclabelka2_array,mcngen_ka1,mcngen_ka2);
    if(mclabf01710>-1){
      mcf01710 = (AliAODMCParticle*) mcArray->At(mclabf01710);
    }
  }

  fCandidateVariables[ 0] = vf0.M();
  fCandidateVariables[ 1] = vf0.Px();
  fCandidateVariables[ 2] = vf0.Py();
  fCandidateVariables[ 3] = vf0.Pz();
  fCandidateVariables[ 4] = 0.493677;
  fCandidateVariables[ 5] = trk1->Px();
  fCandidateVariables[ 6] = trk1->Py();
  fCandidateVariables[ 7] = trk1->Pz();
  fCandidateVariables[ 8] = 0.493677;
  fCandidateVariables[ 9] = trk2->Px();
  fCandidateVariables[10] = trk2->Py();
  fCandidateVariables[11] = trk2->Pz();

  fCandidateVariables[14] = (Float_t) aod->GetNumberOfTracks();
  fCandidateVariables[15] = 0;
  fCandidateVariables[16] = fMultiplicityVZEROA;
  fCandidateVariables[17] = fMultiplicityVZEROC;
  if(fUseMCInfo){
    if(mcf01710){
      fCandidateVariables[18] = mcf01710->GetPdgCode();
      fCandidateVariables[19] = mcf01710->Px();
      fCandidateVariables[20] = mcf01710->Py();
      fCandidateVariables[21] = mcf01710->Pz();

      fCandidateVariables[22] = mclabelka1_array[0];
      fCandidateVariables[23] = mclabelka2_array[0];
      fCandidateVariables[24] = mclabf01710;
    }
  }
  fCandidateVariables[25] = fCentrality;
  fCandidateVariables[26] = 1;

  if(vf0.Pt()>4.){
    fVariablesTree->Fill();
  }


  Double_t cont[4];
  cont[0] = vf0.M();
  cont[1] = vf0.Pt();
  if(trk1->Pt()>trk2->Pt()){
    cont[2] = trk2->Pt();
  }else{
    cont[2] = trk1->Pt();
  }
  cont[3] = fCentrality;

  if(trk1->Charge()*trk2->Charge()<0){
    fHistoF01710ChargedMass->Fill(cont);
  }else{
    fHistoF01710ChargedLikeMass->Fill(cont);
  }

  if(fUseMCInfo){
    if(mcf01710){
      Int_t pdgcode = mcf01710->GetPdgCode();
      if(abs(pdgcode)==10331){
        fHistoF01710ChargedMassMCS->Fill(cont);
      }
      if(abs(pdgcode)==335){
        fHistoF21525ChargedMassMCS->Fill(cont);
      }
      if(abs(pdgcode)==115){
        fHistoA21320ChargedMassMCS->Fill(cont);
      }
      if(abs(pdgcode)==225){
        fHistoF21270ChargedMassMCS->Fill(cont);
      }

    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::FillMixROOTObjects(TLorentzVector *v01, TLorentzVector *v02, TVector *v01vars, TVector *v02vars) 
{
  //
  // Fill Mix root object
  //

  for(Int_t i=0;i<26;i++)
    fCandidateVariables[i] = -9999.;

  Double_t mk0s_pdg = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  TLorentzVector vf0,vk0s1,vk0s2;
  vk0s1.SetXYZM(v01->Px(),v01->Py(),v01->Pz(),mk0s_pdg);
  vk0s2.SetXYZM(v02->Px(),v02->Py(),v02->Pz(),mk0s_pdg);
  vf0 = vk0s1+vk0s2;

  fCandidateVariables[ 0] = vf0.M();
  fCandidateVariables[ 1] = vf0.Px();
  fCandidateVariables[ 2] = vf0.Py();
  fCandidateVariables[ 3] = vf0.Pz();
  fCandidateVariables[ 4] = v01->M();
  fCandidateVariables[ 5] = v01->Px();
  fCandidateVariables[ 6] = v01->Py();
  fCandidateVariables[ 7] = v01->Pz();
  fCandidateVariables[ 8] = v02->M();
  fCandidateVariables[ 9] = v02->Px();
  fCandidateVariables[10] = v02->Py();
  fCandidateVariables[11] = v02->Pz();
  fCandidateVariables[12] = (*v01vars)[0];
  fCandidateVariables[13] = (*v02vars)[0];
  fCandidateVariables[14] = -9999;
  fCandidateVariables[15] = 1;
  fCandidateVariables[16] = fMultiplicityVZEROA;
  fCandidateVariables[17] = fMultiplicityVZEROC;

  fCandidateVariables[25] = fCentrality;
  fCandidateVariables[26] = 0;

  if(fWriteVariableTree)
    fVariablesTree->Fill();

  Double_t cont[4];
  cont[0] = vf0.M();
  cont[1] = vf0.Pt();
  if(v01->Pt()>v02->Pt()){
    cont[2] = v02->Pt();
  }else{
    cont[2] = v01->Pt();
  }
  cont[3] = fCentrality;
  fHistoF01710MassMix->Fill(cont);
}
//________________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::DefineTreeVariables() 
{
  //
  // This is to define tree variables
  //
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 27;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="InvMassF0";
  fCandidateVariableNames[ 1]="F0Px";
  fCandidateVariableNames[ 2]="F0Py";
  fCandidateVariableNames[ 3]="F0Pz";
  fCandidateVariableNames[ 4]="KaMass1";
  fCandidateVariableNames[ 5]="KaPx1";
  fCandidateVariableNames[ 6]="KaPy1";
  fCandidateVariableNames[ 7]="KaPz1";
  fCandidateVariableNames[ 8]="KaMass2";
  fCandidateVariableNames[ 9]="KaPx2";
  fCandidateVariableNames[10]="KaPy2";
  fCandidateVariableNames[11]="KaPz2";
  fCandidateVariableNames[12]="V0CosPA1";
  fCandidateVariableNames[13]="V0CosPA2";
  fCandidateVariableNames[14]="nTracks";
  fCandidateVariableNames[15]="Mixing";
  fCandidateVariableNames[16]="MultVZEROA";
  fCandidateVariableNames[17]="MultVZEROC";
  fCandidateVariableNames[18]="PdgCode";
  fCandidateVariableNames[19]="F0PxMC";
  fCandidateVariableNames[20]="F0PyMC";
  fCandidateVariableNames[21]="F0PzMC";
  fCandidateVariableNames[22]="MatchedLabelKa1";
  fCandidateVariableNames[23]="MatchedLabelKa2";
  fCandidateVariableNames[24]="MatchedLabel";
  fCandidateVariableNames[25]="Centrality";
  fCandidateVariableNames[26]="ChargedMode";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

//__________________________________________________________________________
void  AliAnalysisTaskSEF01710fromAODtracks::DefineGeneralHistograms() {
  //
  // This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",18,-0.5,17.5);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(1,"X1");
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"TriggerOK");
  fCEvents->GetXaxis()->SetBinLabel(5,"IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(6,"CascadesHF exists");
  fCEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
  fCEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
  fCEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fCEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fCEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fCEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fCEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fCEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=10cm"));
  fCEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=10cm"));
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fCEventsNorm = new TH1F("fCEventsNorm","conter",5,-0.5,4.5);
  fCEventsNorm->SetStats(kTRUE);
  fCEventsNorm->GetXaxis()->SetBinLabel(1,"CountForNorm");
  fCEventsNorm->GetXaxis()->SetBinLabel(2,"NoPrimaryV");
  fCEventsNorm->GetXaxis()->SetBinLabel(3,"ZvtxGT10");
  fCEventsNorm->GetXaxis()->SetBinLabel(4,"PrimaryV");
  fCEventsNorm->GetXaxis()->SetBinLabel(5,"Analyzed events");
  fCEventsNorm->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger","counter",18,-0.5,17.5);
  fHTrigger->SetStats(kTRUE);
  fHTrigger->GetXaxis()->SetBinLabel(1,"X1");
  fHTrigger->GetXaxis()->SetBinLabel(2,"kMB");
  fHTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(4,"kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(5,"kINT7");
  fHTrigger->GetXaxis()->SetBinLabel(6,"kEMC7");
  //fHTrigger->GetXaxis()->SetBinLabel(7,"Space");
  fHTrigger->GetXaxis()->SetBinLabel(8,"kMB|kSemiCentral|kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(9,"kINT7|kEMC7");
  fHTrigger->GetXaxis()->SetBinLabel(11,"kMB&kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(12,"kMB&kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(13,"kINT7&kEMC7");

  fHCentrality = new TH1F("fHCentrality","conter",100,0.,100.);
  fHVtxZ = new TH1F("fHVtxZ","conter",100,-15.,15.);
  fHMultiplicityV0A = new TH1F("fHMultiplicityV0A","",1000,-0.5,999.5);
  fHMultiplicityV0C = new TH1F("fHMultiplicityV0C","",1000,-0.5,999.5);

  fOutput->Add(fCEvents);
  fOutput->Add(fCEventsNorm);
  fOutput->Add(fHTrigger);
  fOutput->Add(fHCentrality);
  fOutput->Add(fHVtxZ);
  fOutput->Add(fHMultiplicityV0A);
  fOutput->Add(fHMultiplicityV0C);


  return;
}
//__________________________________________________________________________
void  AliAnalysisTaskSEF01710fromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define histograms
  //

  //------------------------------------------------
  // Basic histograms
  //------------------------------------------------
  Int_t bins_base[4]=		{400,40,20,1};
  Double_t xmin_base[4]={0.9,0,0.,0.00};
  Double_t xmax_base[4]={4.9,20.,10.,100};
  fHistoF01710Mass = new THnSparseF("fHistoF01710Mass","",4,bins_base,xmin_base,xmax_base);
  fHistoF01710MassMix = new THnSparseF("fHistoF01710MassMix","",4,bins_base,xmin_base,xmax_base);
  fHistoF01710ChargedMass = new THnSparseF("fHistoF01710ChargedMass","",4,bins_base,xmin_base,xmax_base);
  fHistoF01710ChargedLikeMass = new THnSparseF("fHistoF01710ChargedLikeMass","",4,bins_base,xmin_base,xmax_base);
  fHistoF01710ChargedMassMix = new THnSparseF("fHistoF01710ChargedMassMix","",4,bins_base,xmin_base,xmax_base);
  fHistoF01710ChargedLikeMassMix = new THnSparseF("fHistoF01710ChargedLikeMassMix","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoF01710Mass);
  fOutputAll->Add(fHistoF01710MassMix);
  fOutputAll->Add(fHistoF01710ChargedMass);
  fOutputAll->Add(fHistoF01710ChargedLikeMass);
  fOutputAll->Add(fHistoF01710ChargedMassMix);
  fOutputAll->Add(fHistoF01710ChargedLikeMassMix);

  Int_t bins_mcs[4]=		{100,200,20,1};
  Double_t xmin_mcs[4]={1.1,0,0.,0.00};
  Double_t xmax_mcs[4]={1.9,20.,10.,100};
  fHistoF01710MassMCS = new THnSparseF("fHistoF01710MassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fHistoF21525MassMCS = new THnSparseF("fHistoF21525MassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fHistoA21320MassMCS = new THnSparseF("fHistoA21320MassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fHistoF21270MassMCS = new THnSparseF("fHistoF21270MassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fHistoF01710ChargedMassMCS = new THnSparseF("fHistoF01710ChargedMassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fHistoF21525ChargedMassMCS = new THnSparseF("fHistoF21525ChargedMassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fHistoA21320ChargedMassMCS = new THnSparseF("fHistoA21320ChargedMassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fHistoF21270ChargedMassMCS = new THnSparseF("fHistoF21270ChargedMassMCS","",4,bins_mcs,xmin_mcs,xmax_mcs);
  fOutputAll->Add(fHistoF01710MassMCS);
  fOutputAll->Add(fHistoF21525MassMCS);
  fOutputAll->Add(fHistoA21320MassMCS);
  fOutputAll->Add(fHistoF21270MassMCS);
  fOutputAll->Add(fHistoF01710ChargedMassMCS);
  fOutputAll->Add(fHistoF21525ChargedMassMCS);
  fOutputAll->Add(fHistoA21320ChargedMassMCS);
  fOutputAll->Add(fHistoF21270ChargedMassMCS);

  Int_t bins_mcgen[3]={200, 20, 1};
  Double_t xmin_mcgen[3]={0., -1., 0.00};
  Double_t xmax_mcgen[3]={20., 1., 100};
  fHistoF01710MassMCGen = new THnSparseF("fHistoF01710MassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fHistoF21525MassMCGen = new THnSparseF("fHistoF21525MassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fHistoA21320MassMCGen = new THnSparseF("fHistoA21320MassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fHistoF21270MassMCGen = new THnSparseF("fHistoF21270MassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fHistoF01710ChargedMassMCGen = new THnSparseF("fHistoF01710ChargedMassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fHistoF21525ChargedMassMCGen = new THnSparseF("fHistoF21525ChargedMassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fHistoA21320ChargedMassMCGen = new THnSparseF("fHistoA21320ChargedMassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fHistoF21270ChargedMassMCGen = new THnSparseF("fHistoF21270ChargedMassMCGen","",3,bins_mcgen,xmin_mcgen,xmax_mcgen);
  fOutputAll->Add(fHistoF01710MassMCGen);
  fOutputAll->Add(fHistoF21525MassMCGen);
  fOutputAll->Add(fHistoA21320MassMCGen);
  fOutputAll->Add(fHistoF21270MassMCGen);
  fOutputAll->Add(fHistoF01710ChargedMassMCGen);
  fOutputAll->Add(fHistoF21525ChargedMassMCGen);
  fOutputAll->Add(fHistoA21320ChargedMassMCGen);
  fOutputAll->Add(fHistoF21270ChargedMassMCGen);

  Int_t bins_rundep[3]={20000, 20, 10};
  Double_t xmin_rundep[3]={-0.5, 0.48, 0.00};
  Double_t xmax_rundep[3]={19999.5, 0.52, 10.};
  fHistonEvtvsRunNumber=new TH1D("fHistonEvtvsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonEvtvsRunNumber);
  fHistonKaonvsRunNumber=new THnSparseF("fHistonKaonvsRunNumber","",3,bins_rundep,xmin_rundep,xmax_rundep);
  fOutputAll->Add(fHistonKaonvsRunNumber);
  fHistonK0vsRunNumber=new THnSparseF("fHistonK0vsRunNumber","",3,bins_rundep,xmin_rundep,xmax_rundep);
  fOutputAll->Add(fHistonK0vsRunNumber);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::SelectV0( const AliVEvent *event,Int_t nV0s,Int_t &nSeleV0, Bool_t *seleV0Flags, TClonesArray *mcArray)
{
  //
  // Select good V0  and return the array of their ids
  //

  nSeleV0 = 0;
  Double_t posVtx[3];
  fVtx1->GetXYZ(posVtx);

  for(Int_t iv0=0;iv0<nV0s;iv0++)
  {
    seleV0Flags[iv0] = kFALSE;
    AliAODv0 *v0 = ((AliAODEvent*)event)->GetV0(iv0);

    if(SingleV0Cuts(v0,fVtx1)){
      seleV0Flags[iv0] = kTRUE;
      nSeleV0++;

      Int_t runnumber_offset = 0;
      Int_t runnumber = event->GetRunNumber();
      if(runnumber<=131000&&runnumber>=114000){
        runnumber_offset = 114000;//lhc10bcde
      }else if(runnumber<=196000&&runnumber>=195000){
        runnumber_offset = 195000;//lhc13bc
      }else if(runnumber<=170593&&runnumber>=167902){
        runnumber_offset = 167902;//lhc11h
      }
      Double_t cont_rundep[3];
      cont_rundep[0] = runnumber-runnumber_offset;
      cont_rundep[1] = v0->MassK0Short();
      cont_rundep[2] = v0->Pt();
      fHistonK0vsRunNumber->Fill(cont_rundep);

      if(fDoEventMixing){
        TVector *varvec = new TVector(1);
        (*varvec)[0] = v0->CosPointingAngle(posVtx);

        Double_t Ek0s = sqrt(pow(v0->P(),2)+pow(v0->MassK0Short(),2));
        fK0Tracks->AddLast(new TLorentzVector(v0->Px(),v0->Py(),v0->Pz(),Ek0s));
        fK0CutVarsArray->AddLast(varvec);
      }
    }
  }
}

//_________________________________________________________________
Int_t AliAnalysisTaskSEF01710fromAODtracks::GetPoolIndex(Double_t zvert, Double_t mult){
  //
  // check in which of the pools the current event falls
  //

  Int_t theBinZ=TMath::BinarySearch(fNzVtxBins,fZvtxBins,zvert);
  if(theBinZ<0 || theBinZ>=fNzVtxBins) return -1;
  Int_t theBinM=TMath::BinarySearch(fNCentBins,fCentBins,mult);
  if(theBinM<0 || theBinM>=fNCentBins) return -1;
  return fNCentBins*theBinZ+theBinM;
}
//_________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::ResetPool(Int_t poolIndex){
  //
  // delete the contets of the pool
  //
  if(poolIndex<0 || poolIndex>=fNOfPools) return;
  delete fEventBuffer[poolIndex];
  fEventBuffer[poolIndex]=new TTree(Form("EventBuffer_%d",poolIndex), "Temporary buffer for event mixing");

  fEventBuffer[poolIndex]->Branch("zVertex", &fVtxZ);
  fEventBuffer[poolIndex]->Branch("centrality", &fCentrality);
  fEventBuffer[poolIndex]->Branch("eventInfo", "TObjString",&fEventInfo);
  fEventBuffer[poolIndex]->Branch("k0array", "TObjArray", &fK0Tracks);
  fEventBuffer[poolIndex]->Branch("k0varsarray", "TObjArray", &fK0CutVarsArray);
  fEventBuffer[poolIndex]->Branch("kaarray", "TObjArray", &fKaonTracks);
  fEventBuffer[poolIndex]->Branch("kavarsarray", "TObjArray", &fKaonCutVarsArray);

  return;
}
//_________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::DoEventMixingWithPools(Int_t poolIndex)
{
  //
  // perform mixed event analysis
  //

  if(poolIndex<0 || poolIndex>fNzVtxBins*fNCentBins) return;
  if(fEventBuffer[poolIndex]->GetEntries()<fNumberOfEventsForMixing) return;

  Int_t nK0 = fK0Tracks->GetEntries();
  Int_t nKaon = fKaonTracks->GetEntries();
  Int_t nEvents=fEventBuffer[poolIndex]->GetEntries();

  TObjArray* k0array=0x0;
  TObjArray* k0varsarray=0x0;
  TObjArray* kaarray=0x0;
  TObjArray* kavarsarray=0x0;
  Float_t zVertex,cent;
  TObjString* eventInfo=0x0;
  fEventBuffer[poolIndex]->SetBranchAddress("eventInfo",&eventInfo);
  fEventBuffer[poolIndex]->SetBranchAddress("zVertex", &zVertex);
  fEventBuffer[poolIndex]->SetBranchAddress("centrality", &cent);
  fEventBuffer[poolIndex]->SetBranchAddress("k0array", &k0array);
  fEventBuffer[poolIndex]->SetBranchAddress("k0varsarray", &k0varsarray);
  fEventBuffer[poolIndex]->SetBranchAddress("kaarray", &kaarray);
  fEventBuffer[poolIndex]->SetBranchAddress("kavarsarray", &kavarsarray);

  for (Int_t i=0; i<nK0; i++)
  {
    TLorentzVector* trk=(TLorentzVector*) fK0Tracks->At(i);
    if(!trk)continue;
    TVector *varsarray1 = (TVector*)fK0CutVarsArray->At(i);

    for(Int_t iEv=0; iEv<fNumberOfEventsForMixing; iEv++){
      fEventBuffer[poolIndex]->GetEvent(iEv + nEvents - fNumberOfEventsForMixing);

      Int_t nV01=k0array->GetEntries();
      for(Int_t iTr1=0; iTr1<nV01; iTr1++){
        TLorentzVector* v0=(TLorentzVector*)k0array->At(iTr1);
        if(!v0) continue;
        TVector *varsarray2 = (TVector*) k0varsarray->At(iTr1);
        FillMixROOTObjects(trk,v0,varsarray1,varsarray2);
      }//v0 loop
    }//event loop
  }//track loop

  for (Int_t i=0; i<nK0; i++)
  {
    TLorentzVector* trk=(TLorentzVector*) fKaonTracks->At(i);
    if(!trk)continue;
    TVector *varsarray1 = (TVector*)fKaonCutVarsArray->At(i);

    for(Int_t iEv=0; iEv<fNumberOfEventsForMixing; iEv++){
      fEventBuffer[poolIndex]->GetEvent(iEv + nEvents - fNumberOfEventsForMixing);

      Int_t nV01=kaarray->GetEntries();
      for(Int_t iTr1=0; iTr1<nV01; iTr1++){
        TLorentzVector* v0=(TLorentzVector*)kaarray->At(iTr1);
        if(!v0) continue;
        TVector *varsarray2 = (TVector*) kavarsarray->At(iTr1);

        Double_t mk0s_pdg = TDatabasePDG::Instance()->GetParticle(310)->Mass();
        TLorentzVector vf0,vka1,vka2;
        vka1.SetXYZM(trk->Px(),trk->Py(),trk->Pz(),0.493677);
        vka2.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),0.493677);
        vf0 = vka1+vka2;

        Double_t cont[4];
        cont[0] = vf0.M();
        cont[1] = vf0.Pt();
        if(vka1.Pt()>vka2.Pt()){
          cont[2] = vka2.Pt();
        }else{
          cont[2] = vka1.Pt();
        }
        cont[3] = fCentrality;

        if(trk->T()*v0->T()<0){
          fHistoF01710ChargedMassMix->Fill(cont);
        }else{
          fHistoF01710ChargedLikeMassMix->Fill(cont);
        }


      }//v0 loop
    }//event loop
  }//track loop

  delete eventInfo;
  if(k0array) k0array->Delete();
  delete k0array;
  if(k0varsarray) k0varsarray->Delete();
  delete k0varsarray;
  if(kaarray) kaarray->Delete();
  delete kaarray;
  if(kavarsarray) kavarsarray->Delete();
  delete kavarsarray;

  fEventBuffer[poolIndex]->SetBranchAddress("zVertex", &fVtxZ);
  fEventBuffer[poolIndex]->SetBranchAddress("centrality", &fCentrality);
  fEventBuffer[poolIndex]->SetBranchAddress("eventInfo",&fEventInfo);
  fEventBuffer[poolIndex]->SetBranchAddress("k0array", &fK0Tracks);
  fEventBuffer[poolIndex]->SetBranchAddress("k0varsarray", &fK0CutVarsArray);
  fEventBuffer[poolIndex]->SetBranchAddress("kaarray", &fKaonTracks);
  fEventBuffer[poolIndex]->SetBranchAddress("kavarsarray", &fKaonCutVarsArray);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEF01710fromAODtracks::SingleV0Cuts(AliAODv0 *v0, AliAODVertex *primVert)
{
  //
  // Single V0 Cut to be applied before object creation
  //
  Bool_t fUseOnTheFlyV0 = kFALSE;
  Int_t fProdV0DaughterTPCClusterMin = 70;
  Double_t fProdV0MassTolK0s = 0.01;
  Double_t fProdV0MassRejLambda = 0.005;
  Double_t fProdV0MassRejPhoton = 0.3;
  Double_t fProdV0DcaDaughtersMax = 1.5;
  Double_t fProdV0CosPointingAngleToPrimVtxMin = 0.97;
  Double_t fProdV0PtMin = 0.0;
  Double_t fProdV0DaughterEtaRange = 0.8;
  Double_t fProdV0DaughterPtMin = 0.10;

  Bool_t onFlyV0 = v0->GetOnFlyStatus(); // on-the-flight V0s
  if ( onFlyV0 && !fUseOnTheFlyV0 ) return kFALSE;

  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
  if(!cptrack || !cntrack) return kFALSE;
  if ( cptrack->Charge() == cntrack->Charge() ) return kFALSE;
  if(!(cptrack->GetStatus() & AliESDtrack::kTPCrefit) ||
      !(cntrack->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
  AliAODVertex *maybeKinkPos = (AliAODVertex*)cptrack->GetProdVertex();
  AliAODVertex *maybeKinkNeg = (AliAODVertex*)cntrack->GetProdVertex();
  if (maybeKinkPos->GetType()==AliAODVertex::kKink || maybeKinkNeg->GetType()==AliAODVertex::kKink) 
    return kFALSE;

  if ( ( ( cptrack->GetTPCClusterInfo(2,1) ) < (Float_t)fProdV0DaughterTPCClusterMin ) || 
      ( ( cntrack->GetTPCClusterInfo(2,1) ) < (Float_t)fProdV0DaughterTPCClusterMin) ) return kFALSE;

  Double_t massK0S = v0->MassK0Short();
  Double_t mk0sPDG   = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  if(fabs(massK0S-mk0sPDG)>fProdV0MassTolK0s) return kFALSE;

  Double_t massLambda = v0->MassLambda();
  Double_t massAntiLambda = v0->MassAntiLambda();
  Double_t mlamPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  if((fabs(massAntiLambda-mlamPDG)<fProdV0MassRejLambda) || (fabs(massLambda-mlamPDG)<fProdV0MassRejLambda)) return kFALSE;

  Double_t pxe1 = v0->MomPosX();
  Double_t pye1 = v0->MomPosY();
  Double_t pze1 = v0->MomPosZ();
  Double_t Ee1 = sqrt(pxe1*pxe1+pye1*pye1+pze1*pze1+0.000511*0.000511);
  Double_t pxe2 = v0->MomNegX();
  Double_t pye2 = v0->MomNegY();
  Double_t pze2 = v0->MomNegZ();
  Double_t Ee2 = sqrt(pxe2*pxe2+pye2*pye2+pze2*pze2+0.000511*0.000511);
  Double_t mphoton = sqrt(pow(Ee1+Ee2,2)-pow(pxe1+pxe2,2)-pow(pye1+pye2,2)-pow(pze1+pze2,2));
  if(mphoton<fProdV0MassRejPhoton) return kFALSE;


  if(TMath::Abs(v0->DcaV0Daughters())>fProdV0DcaDaughtersMax) return kFALSE;
  Double_t posVtx[3] = {0.,0.,0.};
  primVert->GetXYZ(posVtx);
  Double_t cospav0 = v0->CosPointingAngle(posVtx); 
  if(cospav0<fProdV0CosPointingAngleToPrimVtxMin) return kFALSE;
  if(v0->Pt()<fProdV0PtMin) return kFALSE;
  if(fabs(cptrack->Eta())>fProdV0DaughterEtaRange) return kFALSE;
  if(fabs(cntrack->Eta())>fProdV0DaughterEtaRange) return kFALSE;
  if(cptrack->Pt()<fProdV0DaughterPtMin) return kFALSE;
  if(cntrack->Pt()<fProdV0DaughterPtMin) return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEF01710fromAODtracks::IsEventSelected(AliVEvent *event) {
  //
  // Event selection 
  // 
  Bool_t accept=kTRUE;
  fWhyRejection = 0;

  //----------------------- Setup -----------------------------
  TString  TriggerClass[2]; 
  TriggerClass[0] = "";
  TriggerClass[1] = "";
  ULong64_t TriggerMask = AliVEvent::kAnyINT;
  Int_t MinVtxType = 3;
  Int_t MinVtxContr = 1;   
  Float_t MaxVtxZ = 10.;
  Int_t cutc=3;
  Double_t cutz=0.6;
  Bool_t DoPileUpRejection = kTRUE;

  if(fAnalysisType==0){
    TriggerClass[0] = "CINT1";
    TriggerMask = AliVEvent::kAnyINT;
    DoPileUpRejection = kTRUE;
  }

  // check if it's MC
  Bool_t isMC=kFALSE;
  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(mcArray) {isMC=kTRUE;}

  // trigger class
  TString firedTriggerClasses=((AliAODEvent*)event)->GetFiredTriggerClasses();
  if(!isMC){
    if(!firedTriggerClasses.Contains(TriggerClass[0].Data()) && 
        (TriggerClass[1].CompareTo("")==0 || !firedTriggerClasses.Contains(TriggerClass[1].Data())) ) {
      fWhyRejection = 5;
      accept=kFALSE;
    }
  }

  // physics selection requirements
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & TriggerMask);
  if(!isSelected) {
    if(accept) fWhyRejection = 7;
    accept=kFALSE;
  }

  // vertex requirements
  const AliVVertex *vertex = event->GetPrimaryVertex();
  if(!vertex){
    accept=kFALSE;
  }else{
    TString title=vertex->GetTitle();
    if(title.Contains("Z") && MinVtxType>1){
      accept=kFALSE;
    }
    else if(title.Contains("3D") && MinVtxType>2){
      accept=kFALSE;
    }
    if(vertex->GetNContributors()<MinVtxContr){
      accept=kFALSE;
    }
    if(TMath::Abs(vertex->GetZ())>MaxVtxZ) {
      if(accept) fWhyRejection = 6;
      accept=kFALSE;
    } 
  }

  // pile-up rejection
  if(!isMC && DoPileUpRejection){
    if(event->IsPileupFromSPD(cutc,cutz,3.,2.,10.)) {
      if(accept)fWhyRejection = 1;
      accept=kFALSE;
    }
  }

  return accept;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEF01710fromAODtracks::IsEventAcceptedPhysicsSelection(AliVEvent *event) {
  //
  // Event selection 
  // 
  Bool_t accept=kTRUE;

  //----------------------- Setup -----------------------------
  ULong64_t TriggerMask = AliVEvent::kAnyINT;
  if(fAnalysisType==0){
    TriggerMask = AliVEvent::kAnyINT;
  }

  // physics selection requirements
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & TriggerMask);
  if(!isSelected) {
    accept=kFALSE;
  }

  return accept;
}

Bool_t AliAnalysisTaskSEF01710fromAODtracks::IsEventAcceptedTrigger(AliVEvent *event) {
  //
  // Event selection 
  // 
  Bool_t accept=kTRUE;

  //----------------------- Setup -----------------------------
  TString  TriggerClass[2]; 
  TriggerClass[0] = "";
  TriggerClass[1] = "";
  if(fAnalysisType==0){
    TriggerClass[0] = "CINT1";
  }

  // check if it's MC
  Bool_t isMC=kFALSE;
  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(mcArray) {isMC=kTRUE;}

  // trigger class
  TString firedTriggerClasses=((AliAODEvent*)event)->GetFiredTriggerClasses();
  if(!isMC){
    if(!firedTriggerClasses.Contains(TriggerClass[0].Data()) && 
        (TriggerClass[1].CompareTo("")==0 || !firedTriggerClasses.Contains(TriggerClass[1].Data())) ) {
      accept=kFALSE;
    }
  }

  return accept;
}
//_________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::MakeMCAnalysis(TClonesArray *mcArray)
{
  //
  // Analyze AliAODmcparticle
  //
  Int_t nmcpart = mcArray->GetEntriesFast();

  for(Int_t i=0;i<nmcpart;i++)
  {
    AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(i);
    if((TMath::Abs(mcpart->GetPdgCode())==10331) || (TMath::Abs(mcpart->GetPdgCode())==335)|| (TMath::Abs(mcpart->GetPdgCode())==115)|| (TMath::Abs(mcpart->GetPdgCode())==225)){
      Bool_t k01_flag = kFALSE;
      Bool_t k02_flag = kFALSE;
      Bool_t ka1_flag = kFALSE;
      Bool_t ka2_flag = kFALSE;
      AliAODMCParticle *mck01 = 0;
      AliAODMCParticle *mck02 = 0;
      AliAODMCParticle *mcka1 = 0;
      AliAODMCParticle *mcka2 = 0;
      Int_t ndau = mcpart->GetNDaughters();
      Int_t nk0s = 0;
      Int_t nkas = 0;

      if(ndau!=2) continue;


      for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
      {
        if(idau<0) break;
        AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
        if(!mcdau) continue;
        if(TMath::Abs(mcdau->GetPdgCode())==310 && nk0s == 0){
          k01_flag = kTRUE;
          mck01 = mcdau;
          nk0s ++;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==310 && nk0s == 1){
          k02_flag = kTRUE;
          mck02 = mcdau;
          nk0s ++;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==321 && nkas == 0){
          ka1_flag = kTRUE;
          mcka1 = mcdau;
          nkas ++;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==321 && nkas == 1){
          ka2_flag = kTRUE;
          mcka2 = mcdau;
          nkas ++;
        }
      }

      Int_t pdgcode = mcpart->GetPdgCode();
      if(k01_flag && k02_flag){
        Double_t contmc[3];
        contmc[0] = mcpart->Pt();
        contmc[1] = mcpart->Y();
        contmc[2] = fCentrality;
        if(abs(pdgcode)==10331){
          fHistoF01710MassMCGen->Fill(contmc);
        }
        if(abs(pdgcode)==335){
          fHistoF21525MassMCGen->Fill(contmc);
        }
        if(abs(pdgcode)==115){
          fHistoA21320MassMCGen->Fill(contmc);
        }
        if(abs(pdgcode)==225){
          fHistoF21270MassMCGen->Fill(contmc);
        }
      }
      if(ka1_flag && ka2_flag){
        Double_t contmc[3];
        contmc[0] = mcpart->Pt();
        contmc[1] = mcpart->Y();
        contmc[2] = fCentrality;
        if(abs(pdgcode)==10331){
          fHistoF01710ChargedMassMCGen->Fill(contmc);
        }
        if(abs(pdgcode)==335){
          fHistoF21525ChargedMassMCGen->Fill(contmc);
        }
        if(abs(pdgcode)==115){
          fHistoA21320ChargedMassMCGen->Fill(contmc);
        }
        if(abs(pdgcode)==225){
          fHistoF21270ChargedMassMCGen->Fill(contmc);
        }
      }
    }
  }

  return;
}
//________________________________________________________________________
Int_t AliAnalysisTaskSEF01710fromAODtracks::MatchToMC(AliAODv0 *v01, AliAODv0 *v02, TClonesArray *mcArray, Int_t *pdgarray_v01, Int_t *pdgarray_v02, Int_t *labelarray_v01, Int_t *labelarray_v02, Int_t &ngen_v01, Int_t &ngen_v02)
{
  //
  // Match to MC
  //
  for(Int_t i=0;i<100;i++){
    pdgarray_v01[i] = -9999;
    labelarray_v01[i] = -9999;
    pdgarray_v02[i] = -9999;
    labelarray_v02[i] = -9999;
  }
  ngen_v01 = 0;
  ngen_v02 = 0;

  // K0s matching
  Int_t pdgdgv0[2]={211,211};
  Int_t labV01 = v01->MatchToMC(310,mcArray,2,pdgdgv0); // the V0
  if(labV01<0) return -1;
  AliAODMCParticle *mcv01 = (AliAODMCParticle*)mcArray->At(labV01);
  if(!mcv01) return -1;
  labelarray_v01[0] = labV01;
  pdgarray_v01[0] = mcv01->GetPdgCode();
  ngen_v01 ++;

  AliAODMCParticle *mcprimv01=0;
  mcprimv01 = mcv01;
  while(mcprimv01->GetMother()>=0) {
    Int_t labprim_v01=mcprimv01->GetMother();
    AliAODMCParticle *tmcprimv01 = (AliAODMCParticle*)mcArray->At(labprim_v01);
    if(!tmcprimv01) {
      break;
    }

    mcprimv01 = tmcprimv01;
    pdgarray_v01[ngen_v01] = mcprimv01->GetPdgCode();
    labelarray_v01[ngen_v01] = labprim_v01;
    ngen_v01 ++;
    if(ngen_v01==100) break;
  }

  Int_t labV02 = v02->MatchToMC(310,mcArray,2,pdgdgv0); // the V0
  if(labV02<0) return -1;
  AliAODMCParticle *mcv02 = (AliAODMCParticle*)mcArray->At(labV02);
  if(!mcv02) return -1;
  labelarray_v02[0] = labV02;
  pdgarray_v02[0] = mcv02->GetPdgCode();
  ngen_v02 ++;

  AliAODMCParticle *mcprimv02=0;
  mcprimv02 = mcv02;
  while(mcprimv02->GetMother()>=0) {
    Int_t labprim_v02=mcprimv02->GetMother();
    AliAODMCParticle *tmcprimv02 = (AliAODMCParticle*)mcArray->At(labprim_v02);
    if(!tmcprimv02) {
      break;
    }

    mcprimv02 = tmcprimv02;
    pdgarray_v02[ngen_v02] = mcprimv02->GetPdgCode();
    labelarray_v02[ngen_v02] = labprim_v02;
    ngen_v02 ++;
    if(ngen_v02==100) break;
  }

  //Check whether they are coming from the same particle
  Bool_t same_flag = kFALSE;
  Int_t matchedlabel=-9999;
  Int_t matchedpdg=-9999;

  for(Int_t iv01=0;iv01<ngen_v01;iv01++){
    for(Int_t iv02=0;iv02<ngen_v02;iv02++){
      if(labelarray_v01[iv01]==labelarray_v02[iv02]){
        matchedpdg = pdgarray_v01[iv01];
        matchedlabel = labelarray_v01[iv01];
        same_flag = kTRUE;
        break;
      }//if matched  (k0-pi)
    }//pi1 loop
    if(same_flag) break;
  }//k0 loop

  return matchedlabel;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSEF01710fromAODtracks::MatchToMC_ChargedKaon(AliAODTrack *trk1, AliAODTrack *trk2, TClonesArray *mcArray, Int_t *pdgarray_t1, Int_t *pdgarray_t2, Int_t *labelarray_t1, Int_t *labelarray_t2, Int_t &ngen_t1, Int_t &ngen_t2)
{
  //
  // Match to MC
  //
  for(Int_t i=0;i<100;i++){
    pdgarray_t1[i] = -9999;
    labelarray_t1[i] = -9999;
    pdgarray_t2[i] = -9999;
    labelarray_t2[i] = -9999;
  }
  ngen_t1 = 0;
  ngen_t2 = 0;

  // Ka 1
  Int_t labT1 = trk1->GetLabel();
  if(labT1<0) return -1;
  AliAODMCParticle *mctrk1 = (AliAODMCParticle*)mcArray->At(labT1);
  if(!mctrk1) return -1;
  if(abs(mctrk1->GetPdgCode())!=321) return -1;
  labelarray_t1[0] = labT1;
  pdgarray_t1[0] = mctrk1->GetPdgCode();
  ngen_t1 ++;

  AliAODMCParticle *mcprim1=0;
  mcprim1 = mctrk1;
  while(mcprim1->GetMother()>=0) {
    Int_t labprim_t1=mcprim1->GetMother();
    AliAODMCParticle *tmcprim1 = (AliAODMCParticle*)mcArray->At(labprim_t1);
    if(!tmcprim1) {
      break;
    }

    mcprim1 = tmcprim1;
    pdgarray_t1[ngen_t1] = mcprim1->GetPdgCode();
    labelarray_t1[ngen_t1] = labprim_t1;
    ngen_t1 ++;
    if(ngen_t1==100) break;
  }

  // Ka 2
  Int_t labT2 = trk2->GetLabel();
  if(labT2<0) return -1;
  AliAODMCParticle *mctrk2 = (AliAODMCParticle*)mcArray->At(labT2);
  if(!mctrk2) return -1;
  if(abs(mctrk2->GetPdgCode())!=321) return -1;
  labelarray_t2[0] = labT2;
  pdgarray_t2[0] = mctrk2->GetPdgCode();
  ngen_t2 ++;

  AliAODMCParticle *mcprim2=0;
  mcprim2 = mctrk2;
  while(mcprim2->GetMother()>=0) {
    Int_t labprim_t2=mcprim2->GetMother();
    AliAODMCParticle *tmcprim2 = (AliAODMCParticle*)mcArray->At(labprim_t2);
    if(!tmcprim2) {
      break;
    }

    mcprim2 = tmcprim2;
    pdgarray_t2[ngen_t2] = mcprim2->GetPdgCode();
    labelarray_t2[ngen_t2] = labprim_t2;
    ngen_t2 ++;
    if(ngen_t2==100) break;
  }

  //Check whether they are coming from the same particle
  Bool_t same_flag = kFALSE;
  Int_t matchedlabel=-9999;
  Int_t matchedpdg=-9999;

  for(Int_t it1=0;it1<ngen_t1;it1++){
    for(Int_t it2=0;it2<ngen_t2;it2++){
      if(labelarray_t1[it1]==labelarray_t2[it2]){
        matchedpdg = pdgarray_t1[it1];
        matchedlabel = labelarray_t1[it1];
        same_flag = kTRUE;
        break;
      }//if matched  (k0-pi)
    }//pi1 loop
    if(same_flag) break;
  }//k0 loop

  return matchedlabel;
}
//________________________________________________________________________
void AliAnalysisTaskSEF01710fromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Int_t *seleFlags, TClonesArray *mcArray)
{
  //
  // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //

  if(trkEntries==0) return;

  if(fPIDResponse==0x0){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    fPIDResponse =inputHandler->GetPIDResponse();
  }

  nSeleTrks=0;
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = 0;

    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);

    //if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;


    AliAODTrack *aodt = (AliAODTrack*)track;
    if(!aodt->TestFilterMask(BIT(4))) continue;

    if(aodt->Pt()<0.2) continue;

    Double_t nsigma_tpcka = -9999;
    Double_t nsigma_tofka = -9999;
    nsigma_tpcka = fPIDResponse->NumberOfSigmasTPC(aodt,AliPID::kKaon);
    nsigma_tofka = fPIDResponse->NumberOfSigmasTOF(aodt,AliPID::kKaon);


    if(fabs(nsigma_tpcka)<2 && fabs(nsigma_tofka)<2.){
      seleFlags[i]=1;
      nSeleTrks++;

      Int_t runnumber_offset = 0;
      Int_t runnumber = event->GetRunNumber();
      if(runnumber<=131000&&runnumber>=114000){
        runnumber_offset = 114000;//lhc10bcde
      }else if(runnumber<=196000&&runnumber>=195000){
        runnumber_offset = 195000;//lhc13bc
      }else if(runnumber<=170593&&runnumber>=167902){
        runnumber_offset = 167902;//lhc11h
      }
      Double_t cont_rundep[3];
      cont_rundep[0] = runnumber-runnumber_offset;
      cont_rundep[1] = 0.493677;
      cont_rundep[2] = aodt->Pt();
      fHistonKaonvsRunNumber->Fill(cont_rundep);


      if(fDoEventMixing){
        TVector *varvec = new TVector(1);
        (*varvec)[0] = -9999;

        fKaonTracks->AddLast(new TLorentzVector(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->Charge()));
        fKaonCutVarsArray->AddLast(varvec);
      }
    }
  } // end loop on tracks
}
