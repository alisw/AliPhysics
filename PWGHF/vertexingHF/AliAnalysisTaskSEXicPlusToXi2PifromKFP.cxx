/**************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
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
// Author: Jianhui Zhu (1,2)
// (1) Central China Normal University
// (2) GSI Helmholtz Centre for Heavy Ion Research
// E-mail: zjh@ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <TDatabasePDG.h>
#include <vector>
#include <TVector3.h>
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLine.h"
#include "TList.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisTaskSEXicPlusToXi2PifromKFP.h"
#include "AliPIDResponse.h"

#include "AliAODMCParticle.h"

// includes added to play with KFParticle
#ifndef HomogeneousField
#define HomogeneousField 
#endif

using std::cout;
using std::endl;

class AliAnalysisTaskSEXicPlusToXi2PifromKFP;    // your analysis class

ClassImp(AliAnalysisTaskSEXicPlusToXi2PifromKFP) // classimp: necessary for root

AliAnalysisTaskSEXicPlusToXi2PifromKFP::AliAnalysisTaskSEXicPlusToXi2PifromKFP() :
  AliAnalysisTaskSE(),
  fIsMC(kFALSE),
  fPID(0),
  fAnaCuts(0),
  fpVtx(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
  fAodTrackInd(0),
  fOutputList(0),
  fListCuts(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_XicPlus(0),
  fVar_XicPlus(0),
  fTree_XicPlusMCGen(0),
  fVar_XicPlusMCGen(0),
  fCounter(0),
  fHistEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fWriteXicPlusTree(kFALSE),
  fWriteXicPlusMCGenTree(kFALSE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSEXicPlusToXi2PifromKFP::AliAnalysisTaskSEXicPlusToXi2PifromKFP(const char* name, AliRDHFCutsKFP* cuts) :
  AliAnalysisTaskSE(name),
  fIsMC(kFALSE),
  fPID(0),
  fAnaCuts(cuts),
  fpVtx(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
  fAodTrackInd(0),
  fOutputList(0),
  fListCuts(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_XicPlus(0),
  fVar_XicPlus(0),
  fTree_XicPlusMCGen(0),
  fVar_XicPlusMCGen(0),
  fCounter(0),
  fHistEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fWriteXicPlusTree(kFALSE),
  fWriteXicPlusMCGenTree(kFALSE)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
  DefineOutput(2, AliNormalizationCounter::Class());
  DefineOutput(3, TTree::Class()); // event
  DefineOutput(4, TTree::Class()); // XicPlus
  DefineOutput(5, TTree::Class()); // XicPlus MCGen
  DefineOutput(6, TList::Class()); // XicPlus event & trigger

}
//_____________________________________________________________________________
AliAnalysisTaskSEXicPlusToXi2PifromKFP::~AliAnalysisTaskSEXicPlusToXi2PifromKFP()
{
    // destructor
    if (fOutputList) {
      delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
      fOutputList = 0;
    }

    if (fListCuts) {
      delete fListCuts;
      fListCuts = 0;
    }

    if (fAnaCuts) {
      delete fAnaCuts;
      fAnaCuts = 0;
    }

    if (fTree_Event) {
      delete fTree_Event;
      fTree_Event = 0;
    }

    if (fVar_Event) {
      delete fVar_Event;
      fVar_Event = 0;
    }

    if (fTree_XicPlus) {
      delete fTree_XicPlus;
      fTree_XicPlus = 0;
    }

    if (fVar_XicPlus) {
      delete fVar_XicPlus;
      fVar_XicPlus = 0;
    }

    if (fTree_XicPlusMCGen) {
      delete fTree_XicPlusMCGen;
      fTree_XicPlusMCGen = 0;
    }

    if (fVar_XicPlusMCGen) {
      delete fVar_XicPlusMCGen;
      fVar_XicPlusMCGen = 0;
    }

    if (fCounter) {
      delete fCounter;
      fCounter = 0;
    }


}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::Init()
{
  // Initialization

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsKFP(*fAnaCuts));
  PostData(1, fListCuts);

  return;

}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
  fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

  DefineAnaHist(); // define analysis histograms

    // example of a histogram
  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 18, 0.5, 18.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"AliAODVertex exists");
  fHistEvents->GetXaxis()->SetBinLabel(3,"TriggerOK");
  fHistEvents->GetXaxis()->SetBinLabel(4,"IsEventSelected");
  fHistEvents->GetXaxis()->SetBinLabel(5,"V0 exists");
  fHistEvents->GetXaxis()->SetBinLabel(6,"Cascade exists");

  fHistEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
  fHistEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
  fHistEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fHistEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fHistEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fHistEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fHistEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fHistEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fHistEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm", fAnaCuts->GetMaxVtxZ()));
  fHistEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fHistEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fHistEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnaCuts->GetMaxVtxZ()));
  fHistEvents->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger", "counter", 18, -0.5, 17.5);                                      
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

  fHCentrality = new TH1F("fHCentrality", "counter", 100, 0., 100.);

  fOutputList->Add(fHistEvents); // don't forget to add it to the list! the list will be written to file, so if you want
  fOutputList->Add(fHTrigger);

  // Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
  if(cont) normName = (TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(2, fCounter);
  DefineEvent();
  PostData(3, fTree_Event);  // postdata will notify the analysis manager of changes / updates to the 

  DefineTreeRecXicPlus();
  PostData(4, fTree_XicPlus);

  DefineTreeGenXicPlus();
  PostData(5, fTree_XicPlusMCGen);

  PostData(6, fOutputList);

  return;
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::UserExec(Option_t *)
{
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you 
  // have access to the current event. 
  // once you return from the UserExec function, the manager will retrieve the next event from the chain

  if (!fInputEvent) { // if the event is empty (getting it failed) skip this event
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* AODEvent = dynamic_cast<AliAODEvent*>(fInputEvent);    // get an event (called AODEvent) from the input file
                                                        // there's another event format (ESD) which works in a similar way
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's

  fHistEvents->Fill(1);

  //--------------------------------------------------------------
  // First check if the event has magnetic field and proper vertex
  //--------------------------------------------------------------

  fBzkG = (Double_t)AODEvent->GetMagneticField();
  if (TMath::Abs(fBzkG)<0.001) return;
  KFParticle::SetField(fBzkG);

  fpVtx = (AliAODVertex*)AODEvent->GetPrimaryVertex();
  if (!fpVtx) return;
  fHistEvents->Fill(2);

  fCounter->StoreEvent(AODEvent,fAnaCuts,fIsMC);

  //------------------------------------------------
  // MC analysis setting                                                                    
  //------------------------------------------------

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if (fIsMC) {
    fMCEvent = MCEvent(); // get the corresponding MC event fMCEvent
    if (!fMCEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(AODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !mcArray ) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fHistEvents->Fill(7); // number of MC array exist

    // load MC header
    mcHeader = (AliAODMCHeader*)AODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if ( !mcHeader ) {
      AliError("AliAnalysisTaskSEXicPlusToXi2PifromKFP::UserExec: MC header branch not found!\n");
      return;
    }
    fHistEvents->Fill(8); // number of MC header exist

    Double_t zMCvtx = mcHeader->GetVtxZ();
    if ( TMath::Abs(zMCvtx) > fAnaCuts->GetMaxVtxZ() ) {
      AliDebug(2,Form("Event rejected: fabs(zVtxMC)=%f > fAnaCuts->GetMaxVtxZ()=%f", zMCvtx, fAnaCuts->GetMaxVtxZ()));
      return;
    } else {
      fHistEvents->Fill(18);
    }
    if ((TMath::Abs(zMCvtx) < fAnaCuts->GetMaxVtxZ()) && (!fAnaCuts->IsEventRejectedDuePhysicsSelection()) && (!fAnaCuts->IsEventRejectedDueToTrigger())) {
      Bool_t selevt = MakeMCAnalysis(mcArray);
      if(!selevt) return;
    }
  }

  //------------------------------------------------
  // Event selection
  //------------------------------------------------
  Bool_t IsTriggerNotOK = fAnaCuts->IsEventRejectedDueToTrigger();
  Bool_t IsPhysSelNotOK = fAnaCuts->IsEventRejectedDuePhysicsSelection();
  Bool_t IsNoVertex = fAnaCuts->IsEventRejectedDueToNotRecoVertex();
  if( !IsTriggerNotOK && !IsPhysSelNotOK && !IsNoVertex && fabs(fpVtx->GetZ())<fAnaCuts->GetMaxVtxZ() ) fHistEvents->Fill(3);

  Bool_t IsEventSelected = fAnaCuts->IsEventSelected(AODEvent);
  if(!IsEventSelected) {
//    cout<<"Why: "<<fAnaCuts->GetWhyRejection()<<endl;
    return;
  }
  fHistEvents->Fill(4);


  Bool_t IsMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  Bool_t IsSemi = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  Bool_t IsCent = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral);
  Bool_t IsINT7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);
  Bool_t IsEMC7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);
  if(IsMB) fHTrigger->Fill(1);
  if(IsSemi) fHTrigger->Fill(2);
  if(IsCent) fHTrigger->Fill(3);
  if(IsINT7) fHTrigger->Fill(4);
  if(IsEMC7) fHTrigger->Fill(5);
  if(IsMB||IsSemi||IsCent) fHTrigger->Fill(7);
  if(IsINT7||IsEMC7) fHTrigger->Fill(8);
  if(IsMB&&IsSemi) fHTrigger->Fill(10);
  if(IsMB&&IsCent) fHTrigger->Fill(11);
  if(IsINT7&&IsEMC7) fHTrigger->Fill(12);

//  AliCentrality *cent = AODEvent->GetCentrality();
//  Float_t Centrality = cent->GetCentralityPercentile("V0M");
//  fHCentrality->Fill(Centrality);

  //------------------------------------------------
  // Check if the event has v0 candidate
  //------------------------------------------------
  Int_t num_v0 = AODEvent->GetNumberOfV0s();
  if (num_v0>0) fHistEvents->Fill(5);

  //------------------------------------------------
  // Check if the event has cascade candidate
  //------------------------------------------------
  Int_t num_casc = AODEvent->GetNumberOfCascades();
  if (num_casc<=0) return;
  fHistEvents->Fill(6);

  // set primary vertex
  KFPVertex pVertex;
  Double_t pos[3],cov[6];
  fpVtx->GetXYZ(pos);
  if ( fabs(pos[2])>10. ) return; // vertex cut on z-axis direction
  fpVtx->GetCovarianceMatrix(cov);
  pVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
  Float_t covF[6];
  for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }
  pVertex.SetCovarianceMatrix(covF);
  pVertex.SetChi2(fpVtx->GetChi2());
  pVertex.SetNDF(fpVtx->GetNDF());
  pVertex.SetNContributors(fpVtx->GetNContributors());

  KFParticle PV(pVertex);

  if(!fAnaCuts) return;

  FillEventROOTObjects();

//------------------------------------------------
// Main analysis done in this function
//------------------------------------------------
  
  fPID = fInputHandler->GetPIDResponse();
  MakeAnaXicPlusFromCasc(AODEvent, mcArray, PV);

  PostData(2, fCounter);
  PostData(3, fTree_Event);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
  PostData(4, fTree_XicPlus);
  PostData(5, fTree_XicPlusMCGen);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
/*
    TCanvas *c1 = new TCanvas();
    TLine *lLPDG = new TLine(7.89, 0, 7.89, 1e10);
    lLPDG->SetLineColor(2);
    lLPDG->Draw();

    TCanvas *c2 = new TCanvas();
    TLine *lXiPDG = new TLine(4.91, 0, 4.91, 1e10);
    lXiPDG->SetLineColor(2);
    lXiPDG->Draw();
*/
    return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MakeMCAnalysis(TClonesArray *mcArray)
{
  // Analyse AliAODMCParticle
  
  for(Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) {
    AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!mcPart) {
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }
    Int_t pdg = mcPart->GetPdgCode();
    if ( TMath::Abs(pdg)!=4232 ) {
      AliDebug(2, Form("MC particle %d is not a Xic+: its pdg code is %d", iPart, pdg));
      continue;
    }
    AliDebug(2, Form("Step 0 ok: MC particle %d is a Xic+: its pdg code is %d", iPart, pdg));
    if ( mcPart->GetNDaughters()!=3 ) {
      AliDebug(2, "Xic+ does not have 3 daughters");
      continue;
    }
    Int_t index_FirstDau = mcPart->GetDaughterFirst();
    if (index_FirstDau<0) continue;
    AliAODMCParticle* mcXicPlusDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(index_FirstDau));
    AliAODMCParticle* mcXicPlusDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(index_FirstDau+1));
    AliAODMCParticle* mcXicPlusDau_2 = dynamic_cast<AliAODMCParticle*>(mcArray->At(index_FirstDau+2));
    if ( TMath::Abs(mcXicPlusDau_0->GetPdgCode())!=211 && TMath::Abs(mcXicPlusDau_0->GetPdgCode())!=3312 ) continue;
    if ( TMath::Abs(mcXicPlusDau_1->GetPdgCode())!=211 && TMath::Abs(mcXicPlusDau_1->GetPdgCode())!=3312 ) continue;
    if ( TMath::Abs(mcXicPlusDau_2->GetPdgCode())!=211 && TMath::Abs(mcXicPlusDau_2->GetPdgCode())!=3312 ) continue;

    // -------------------------------- First daughter is Xi --------------------------------
    if ( TMath::Abs(mcXicPlusDau_0->GetPdgCode())==3312 && TMath::Abs(mcXicPlusDau_1->GetPdgCode())==211 && TMath::Abs(mcXicPlusDau_2->GetPdgCode())==211 ) { // First is Xi
      if (mcXicPlusDau_1->GetPdgCode()!=mcXicPlusDau_2->GetPdgCode()) continue; // Check pion charge is same
      if ( (mcXicPlusDau_0->GetPdgCode()>0&&mcXicPlusDau_1->GetPdgCode()<0) || (mcXicPlusDau_0->GetPdgCode()<0&&mcXicPlusDau_1->GetPdgCode()>0) ) continue; // Check charge of pion and Xi is different
      // Start to check Xi decay
      if ( mcXicPlusDau_0->GetNDaughters()!=2 ) continue;
      Bool_t pifromXi_flag = kFALSE;
      Bool_t v0_flag = kFALSE;
      Bool_t pifromLam_flag = kFALSE;
      Bool_t prfromLam_flag = kFALSE;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcXicPlusDau_0->GetDaughterFirst();idau<=mcXicPlusDau_0->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = dynamic_cast<AliAODMCParticle*>(mcArray->At(idau));
        if(TMath::Abs(mcdau->GetPdgCode())==211){ // 211: pion
          pifromXi_flag = kTRUE;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==2) { // 3122: Lambda
          v0_flag = kTRUE;
          mcv0part = mcdau;
          for(Int_t jdau=mcv0part->GetDaughterFirst();jdau<=mcv0part->GetDaughterLast();jdau++) {
            if (jdau<0) break;
            AliAODMCParticle *mcDau_Lam = (AliAODMCParticle*) mcArray->At(jdau);
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==211) pifromLam_flag = kTRUE;
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==2212) prfromLam_flag = kTRUE;
          }
        }
      }
      if ( pifromXi_flag && v0_flag && pifromLam_flag && prfromLam_flag ) {
        if ( TMath::Abs(mcPart->Y()) < 0.8 ) { // Acceptance
          Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
          FillTreeGenXicPlus(mcPart, CheckOrigin);
        }
      }

      /*
      // Start to check Xi decay (method 2)
      if ( mcXicPlusDau_0->GetNDaughters()!=2 ) continue;
      Int_t label_XiDau_0 = mcXicPlusDau_0->GetDaughterLabel(0);
      Int_t label_XiDau_1 = mcXicPlusDau_0->GetDaughterLabel(1);
      AliAODMCParticle* mcXiDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_XiDau_0));
      AliAODMCParticle* mcXiDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_XiDau_1));
      if ( !mcXiDau_0 || !mcXiDau_1 ) {
        AliDebug(2, "Could not access Xi daughters, continuing...");
        continue;
      }
      Int_t pdgXiDau_0 = mcXiDau_0->GetPdgCode();
      Int_t pdgXiDau_1 = mcXiDau_1->GetPdgCode();
      if ( TMath::Abs(pdgXiDau_0)!=3122 || TMath::Abs(pdgXiDau_0)!=211 ) continue;
      if ( TMath::Abs(pdgXiDau_1)!=3122 || TMath::Abs(pdgXiDau_1)!=211 ) continue;
      if ( TMath::Abs(pdgXiDau_0)==3122 && TMath::Abs(pdgXiDau_1)==211) { // First is Lam
        if ( mcXiDau_0->GetNDaughters()!=2 ) continue;
        Int_t label_LamDau_0 = mcXiDau_0->GetDaughterLabel(0);
        Int_t label_LamDau_1 = mcXiDau_0->GetDaughterLabel(1);
        AliAODMCParticle* mcLamDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_0));
        AliAODMCParticle* mcLamDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_1));
        if ( !mcLamDau_0 || !mcLamDau_1 ) {
          AliDebug(2, "Could not access Lambda daughters, continuing...");
          continue;
        }
        Int_t pdgLamDau_0 = mcLamDau_0->GetPdgCode();
        Int_t pdgLamDau_1 = mcLamDau_1->GetPdgCode();
        if ( (TMath::Abs(pdgLamDau_0)==2212 && TMath::Abs(pdgLamDau_1)==211) || (TMath::Abs(pdgLamDau_1)==2212 && TMath::Abs(pdgLamDau_0)==211) ) {
          if ( TMath::Abs(mcPart->Y()) < 0.8 ) {
            Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
            FillTreeGenXicPlus(mcPart, CheckOrigin);
          }
        }
      } // First is Lam
      if ( TMath::Abs(pdgXiDau_1)==3122 && TMath::Abs(pdgXiDau_0)==211) { // Second is Lam
        if ( mcXiDau_1->GetNDaughters()!=2 ) continue;
        Int_t label_LamDau_0 = mcXiDau_1->GetDaughterLabel(0);
        Int_t label_LamDau_1 = mcXiDau_1->GetDaughterLabel(1);
        AliAODMCParticle* mcLamDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_0));
        AliAODMCParticle* mcLamDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_1));
        if ( !mcLamDau_0 || !mcLamDau_1 ) {
          AliDebug(2, "Could not access Lambda daughters, continuing...");
          continue;
        }
        Int_t pdgLamDau_0 = mcLamDau_0->GetPdgCode();
        Int_t pdgLamDau_1 = mcLamDau_1->GetPdgCode();
        if ( (TMath::Abs(pdgLamDau_0)==2212 && TMath::Abs(pdgLamDau_1)==211) || (TMath::Abs(pdgLamDau_1)==2212 && TMath::Abs(pdgLamDau_0)==211) ) {
          if ( TMath::Abs(mcPart->Y()) < 0.8 ) {
            Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
            FillTreeGenXicPlus(mcPart, CheckOrigin);
          }
        }
      } // Second is Lam
      */
    } // First is Xi

    // -------------------------------- Second daughter is Xi --------------------------------
    if ( TMath::Abs(mcXicPlusDau_1->GetPdgCode())==3312 && TMath::Abs(mcXicPlusDau_0->GetPdgCode())==211 && TMath::Abs(mcXicPlusDau_2->GetPdgCode())==211 ) { // Second is Xi
      if (mcXicPlusDau_0->GetPdgCode()!=mcXicPlusDau_2->GetPdgCode()) continue; // Check pion charge is same
      if ( (mcXicPlusDau_1->GetPdgCode()>0&&mcXicPlusDau_0->GetPdgCode()<0) || (mcXicPlusDau_1->GetPdgCode()<0&&mcXicPlusDau_0->GetPdgCode()>0) ) continue; // Check charge of pion and Xi is different
      // Start to check Xi decay
      if ( mcXicPlusDau_1->GetNDaughters()!=2 ) continue;
      Bool_t pifromXi_flag = kFALSE;
      Bool_t v0_flag = kFALSE;
      Bool_t pifromLam_flag = kFALSE;
      Bool_t prfromLam_flag = kFALSE;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcXicPlusDau_1->GetDaughterFirst();idau<=mcXicPlusDau_1->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = dynamic_cast<AliAODMCParticle*>(mcArray->At(idau));
        if(TMath::Abs(mcdau->GetPdgCode())==211){ // 211: pion
          pifromXi_flag = kTRUE;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==2) { // 3122: Lambda
          v0_flag = kTRUE;
          mcv0part = mcdau;
          for(Int_t jdau=mcv0part->GetDaughterFirst();jdau<=mcv0part->GetDaughterLast();jdau++) {
            if (jdau<0) break;
            AliAODMCParticle *mcDau_Lam = (AliAODMCParticle*) mcArray->At(jdau);
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==211) pifromLam_flag = kTRUE;
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==2212) prfromLam_flag = kTRUE;
          }
        }
      }
      if ( pifromXi_flag && v0_flag && pifromLam_flag && prfromLam_flag ) {
        if ( TMath::Abs(mcPart->Y()) < 0.8 ) { // Acceptance
          Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
          FillTreeGenXicPlus(mcPart, CheckOrigin);
        }
      }
    } // Second is Xi

    // -------------------------------- Third daughter is Xi --------------------------------
    if ( TMath::Abs(mcXicPlusDau_2->GetPdgCode())==3312 && TMath::Abs(mcXicPlusDau_0->GetPdgCode())==211 && TMath::Abs(mcXicPlusDau_1->GetPdgCode())==211 ) { // Third is Xi
      if (mcXicPlusDau_0->GetPdgCode()!=mcXicPlusDau_1->GetPdgCode()) continue; // Check pion charge is same
      if ( (mcXicPlusDau_2->GetPdgCode()>0&&mcXicPlusDau_0->GetPdgCode()<0) || (mcXicPlusDau_2->GetPdgCode()<0&&mcXicPlusDau_0->GetPdgCode()>0) ) continue; // Check charge of pion and Xi is different
      // Start to check Xi decay
      if ( mcXicPlusDau_2->GetNDaughters()!=2 ) continue;
      Bool_t pifromXi_flag = kFALSE;
      Bool_t v0_flag = kFALSE;
      Bool_t pifromLam_flag = kFALSE;
      Bool_t prfromLam_flag = kFALSE;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcXicPlusDau_2->GetDaughterFirst();idau<=mcXicPlusDau_2->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = dynamic_cast<AliAODMCParticle*>(mcArray->At(idau));
        if(TMath::Abs(mcdau->GetPdgCode())==211){ // 211: pion
          pifromXi_flag = kTRUE;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==2) { // 3122: Lambda
          v0_flag = kTRUE;
          mcv0part = mcdau;
          for(Int_t jdau=mcv0part->GetDaughterFirst();jdau<=mcv0part->GetDaughterLast();jdau++) {
            if (jdau<0) break;
            AliAODMCParticle *mcDau_Lam = (AliAODMCParticle*) mcArray->At(jdau);
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==211) pifromLam_flag = kTRUE;
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==2212) prfromLam_flag = kTRUE;
          }
        }
      }
      if ( pifromXi_flag && v0_flag && pifromLam_flag && prfromLam_flag ) {
        if ( TMath::Abs(mcPart->Y()) < 0.8 ) { // Acceptance
          Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
          FillTreeGenXicPlus(mcPart, CheckOrigin);
        }
      }
    } // Third is Xi
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillTreeGenXicPlus(AliAODMCParticle *mcpart, Int_t CheckOrigin)
{
  // Fill histograms or tree depending

  for(Int_t i=0;i<4;i++){
    fVar_XicPlusMCGen[i] = -9999.;
  }

  fVar_XicPlusMCGen[0] = mcpart->Y();
  fVar_XicPlusMCGen[1] = mcpart->Pt();
  fVar_XicPlusMCGen[2] = CheckOrigin;
  fVar_XicPlusMCGen[3] = mcpart->GetPdgCode();

  if (fWriteXicPlusMCGenTree && fVar_XicPlusMCGen[1]>0.9999) fTree_XicPlusMCGen->Fill();

//  fVar_XicPlusMCGen[ 0] = fCentrality;
//  fVar_XicPlusMCGen[ 1] = decaytype;
//  if (mcpart->IsPrimary() && (!mcpart->IsPhysicalPrimary())) fVar_XicPlusMCGen[2] = 1;
//  if (mcpart->IsPhysicalPrimary()) fVar_XicPlusMCGen[2] = 2;
//  if (mcpart->IsSecondaryFromWeakDecay()) fVar_XicPlusMCGen[2] = 3;
//  if (mcpart->IsSecondaryFromMaterial()) fVar_XicPlusMCGen[2] = 4;
//  if (mcpart->IsFromSubsidiaryEvent()) fVar_XicPlusMCGen[2] = 5;
//  fVar_XicPlusMCGen[ 3] = mcpart->Eta();
//  fVar_XicPlusMCGen[ 4] = mcpart->Y();
//  fVar_XicPlusMCGen[ 5] = mcpart->Px();
//  fVar_XicPlusMCGen[ 6] = mcpart->Py();
//  fVar_XicPlusMCGen[ 7] = mcpart->Pz();
//  fVar_XicPlusMCGen[ 8] = mcpipart->Px();
//  fVar_XicPlusMCGen[ 9] = mcpipart->Py();
//  fVar_XicPlusMCGen[10] = mcpipart->Pz();
//  fVar_XicPlusMCGen[11] = mccascpart->Px();
//  fVar_XicPlusMCGen[12] = mccascpart->Py();
//  fVar_XicPlusMCGen[13] = mccascpart->Pz();
//  fVar_XicPlusMCGen[14] = mcpart->GetPdgCode();
//  fVar_XicPlusMCGen[15] = mcpipart->GetPdgCode();
//  fVar_XicPlusMCGen[16] = mccascpart->GetPdgCode();
//  fVar_XicPlusMCGen[17] = fRunNumber;
//  fVar_XicPlusMCGen[18] = fEvNumberCounter;
//  fVar_XicPlusMCGen[19] = CheckOrigin;

  /*
  const Double_t massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t massXi   = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t pipx = mcpipart->Px();
  Double_t pipy = mcpipart->Py();
  Double_t pipz = mcpipart->Pz();
//  Double_t piE  = sqrt(pipx*pipx+pipy*pipy+pipz*pipz+0.000511*0.000511);
  Double_t piE  = sqrt(pipx*pipx+pipy*pipy+pipz*pipz+massPion*massPion);
  Double_t cascpx = mccascpart->Px();
  Double_t cascpy = mccascpart->Py();
  Double_t cascpz = mccascpart->Pz();
  Double_t cascE  = sqrt(cascpx*cascpx+cascpy*cascpy+cascpz*cascpz+massXi*massXi);

  Double_t InvMassPiXi = sqrt(pow(piE+cascE,2)-pow(pipx+cascpx,2)-pow(pipy+cascpy,2)-pow(pipz+cascpz,2));

  Double_t contXicPlusMC[3];
  contXicPlusMC[0] = mcpart->Pt();
  contXicPlusMC[1] = mcpart->Y();
  contXicPlusMC[2] = fCentrality;

  Double_t contPionMC[3];
  contPionMC[0] = mcpipart->Pt();
  contPionMC[1] = mcpipart->Eta();
  contPionMC[2] = fCentrality;

  Double_t contXiMC[3];
  contXiMC[0] = mccascpart->Pt();
  contXiMC[1] = mccascpart->Y();
  contXiMC[2] = fCentrality;

  Double_t contPiXiMassMCGen[3];
  contPiXiMassMCGen[0] = InvMassPiXi;
  contPiXiMassMCGen[1] = mcpart->Pt();
  contPiXiMassMCGen[2] = fCentrality;

  Double_t contPiXiMassvsPiPtMCGen[3];
  contPiXiMassvsPiPtMCGen[0] = InvMassPiXi;
  contPiXiMassvsPiPtMCGen[1] = mcpipart->Pt();
  contPiXiMassvsPiPtMCGen[2] = fCentrality;

  if (decaytype==0) {
    if (fabs(mcpipart->Eta())<fAnaCuts->GetProdTrackEtaRange()) {
    }
  }
  */
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::MakeAnaXicPlusFromCasc(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV)
{
  // Main analysis called from "UserExec"

//  std::cout.setf(std::ios::fixed);
//  std::cout.setf(std::ios::showpoint);
//  std::cout.precision(3);

  // set the magnetic field
  KFParticle::SetField(fBzkG);

  const UInt_t nCasc = AODEvent->GetNumberOfCascades();

  Double_t xyzP[3], xyzN[3];
  Double_t xvyvzvP[3], xvyvzvN[3];
  Double_t covP[21], covN[21], covB[21];
  const Float_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Float_t massXi     = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Float_t massK0S    = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  // select good candidates for pion
  const UInt_t nTracks = AODEvent->GetNumberOfTracks();
  AliAODTrack *trackP[nTracks], *trackN[nTracks];
  Int_t flag_trkP = 0, flag_trkN = 0;

  for (UInt_t itrk=0; itrk<nTracks; itrk++) {
    AliAODTrack *trk = static_cast<AliAODTrack*>(AODEvent->GetTrack(itrk));
    Double_t covtest[21];
    if ( !trk || trk->GetID()<0 || !trk->GetCovarianceXYZPxPyPz(covtest) || !AliVertexingHFUtils::CheckAODtrackCov(trk) ) continue;

    if  (trk->Charge() > 0 ) {
      trackP[flag_trkP] = trk;
      flag_trkP++;
    }
    if  (trk->Charge() < 0 ) {
      trackN[flag_trkN] = trk;
      flag_trkN++;
    }
  }

  for (UInt_t iCasc=0; iCasc<nCasc; iCasc++) {
    AliAODcascade *casc = AODEvent->GetCascade(iCasc);
    // cascade cut
    if ( !fAnaCuts->SingleCascCuts(casc, kFALSE) ) continue;

    AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0));
    AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1));
    AliAODTrack *btrack = (AliAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0));

    if ( !ptrack||!ntrack||!btrack ) continue;

    // check charge of the first daughter, if negative, define it as the second one
    if ( ptrack->Charge()<0 ) {
      ptrack = (AliAODTrack*) (casc->GetDaughter(1));
      ntrack = (AliAODTrack*) (casc->GetDaughter(0));
    }

    if ( !ptrack->GetCovarianceXYZPxPyPz(covP) || !ntrack->GetCovarianceXYZPxPyPz(covN) || !btrack->GetCovarianceXYZPxPyPz(covB) ) continue;

    if ( !AliVertexingHFUtils::CheckAODtrackCov(ptrack) || !AliVertexingHFUtils::CheckAODtrackCov(ntrack) || !AliVertexingHFUtils::CheckAODtrackCov(btrack) ) continue;

    KFParticle kfpProton     = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack, 2212);
    KFParticle kfpPionMinus  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack, -211);
    KFParticle kfpAntiProton = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack, -2212);
    KFParticle kfpPionPlus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack, 211);

    KFParticle kfpElePlus    = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack, -11);
    KFParticle kfpEleMinus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack, 11);

    // === K0S ===
    KFParticle kfpK0Short;
    const KFParticle *vk0sDaughters[2]  = {&kfpPionPlus, &kfpPionMinus};
    kfpK0Short.Construct(vk0sDaughters, 2);
    // ============
    // === Gamma ===
    KFParticle kfpGamma;
    const KFParticle *vGammaDaughters[2]  = {&kfpElePlus, &kfpEleMinus};
    kfpGamma.Construct(vGammaDaughters, 2);
    // =============

    if ( btrack->Charge()<0 ) { // Xi^-

      const KFParticle *vDaughters[2] = {&kfpProton, &kfpPionMinus};

      KFParticle kfpLambda;
      kfpLambda.Construct(vDaughters, 2);
      Float_t massLambda_Rec, err_massLambda;
      kfpLambda.GetMass(massLambda_Rec, err_massLambda);

      // check rapidity of lambda
      if ( TMath::Abs(kfpLambda.GetE())<=TMath::Abs(kfpLambda.GetPz()) ) continue;

      // chi2>0 && NDF>0 for selecting Lambda
      if ( (kfpLambda.GetNDF()<=0 || kfpLambda.GetChi2()<=0) ) continue;

      // check cov. of Lambda
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLambda) ) continue;

      // err_mass>0 of Lambda
      if ( err_massLambda<=0 ) continue;

      // Chi2geo cut of Lambda
      if ( (kfpLambda.GetChi2()/kfpLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue; 

      //************************** calculate l/Δl for Lambda *************************************
      Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
      Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
      Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
      Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
      Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
      if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
      dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
      Double_t nErr_l_Lambda = l_Lambda/dl_Lambda;
      //***************************************************************************************

      // l/Deltal cut of Lambda
      if ( nErr_l_Lambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Lambda
      if ( TMath::Abs(massLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      KFParticle kfpLambda_m = kfpLambda;
      kfpLambda_m.SetNonlinearMassConstraint(massLambda);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLambda_m) || TMath::Abs(kfpLambda_m.GetE()) <= TMath::Abs(kfpLambda_m.GetPz()) ) continue;

      KFParticle kfpPionOrKaon;
      kfpPionOrKaon = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, -211); // pion-
      KFParticle kfpXiMinus;
      const KFParticle *vXiDs[2] = {&kfpPionOrKaon, &kfpLambda_m};
      kfpXiMinus.Construct(vXiDs, 2);

      // check rapidity of Xi-
      if ( TMath::Abs(kfpXiMinus.GetE())<=TMath::Abs(kfpXiMinus.GetPz()) ) continue;

      // err_massXi > 0
      Float_t massXiMinus_Rec, err_massXiMinus;
      kfpXiMinus.GetMass(massXiMinus_Rec, err_massXiMinus);
      if ( err_massXiMinus<=0 ) continue;

      // chi2>0 && NDF>0
      if ( kfpXiMinus.GetNDF()<=0 || kfpXiMinus.GetChi2()<=0 ) continue;

      // Prefilter
      if ( kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

      // check covariance matrix
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus) ) continue;

      // mass window cut of Xi-
      if ( TMath::Abs(massXiMinus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi()) ) continue;

      KFParticle kfpXiMinus_m = kfpXiMinus;
      kfpXiMinus_m.SetNonlinearMassConstraint(massXi);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus_m) || TMath::Abs(kfpXiMinus_m.GetE()) <= TMath::Abs(kfpXiMinus_m.GetPz()) ) continue;

      for (Int_t itrkBP_trk1=0; itrkBP_trk1<flag_trkP-1; itrkBP_trk1++) { // Loop for first bachelor pion+
        for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkP; itrkBP_trk2++) { // Loop for second bachelor pion+

          if ( trackP[itrkBP_trk1]->GetID() == ptrack->GetID() || trackP[itrkBP_trk2]->GetID() == ptrack->GetID() || trackP[itrkBP_trk1]->GetID() == trackP[itrkBP_trk2]->GetID() ) continue;

          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk1]) ) continue;
          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk2]) ) continue;

          KFParticle kfpBP_trk1  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackP[itrkBP_trk1], 211);
          KFParticle kfpBP_trk2   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackP[itrkBP_trk2], 211);

          // reconstruct XicPlus
          KFParticle kfpXicPlus;
          const KFParticle *vXicPlusDs[3] = {&kfpXiMinus_m, &kfpBP_trk1, &kfpBP_trk2};
          kfpXicPlus.Construct(vXicPlusDs, 3);

          // chi2>0 && NDF>0
          if ( kfpXicPlus.GetNDF()<=0 || kfpXicPlus.GetChi2()<=0 ) continue;

          // Prefilter
          if ( kfpXicPlus.GetChi2()/kfpXicPlus.GetNDF() >= fAnaCuts->GetKFPXicPlus_Chi2geoMax() ) continue;
          if ( kfpXicPlus.GetPt() < fAnaCuts->GetPtMinXicPlus() ) continue;

          // check rapidity of XicPlus
          if ( TMath::Abs(kfpXicPlus.GetE())<=TMath::Abs(kfpXicPlus.GetPz()) ) continue;

          // check covariance matrix
          if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus) ) continue;

          // err_massXicPlus > 0
          Float_t massXicPlus_Rec, err_massXicPlus;
          kfpXicPlus.GetMass(massXicPlus_Rec, err_massXicPlus);
          if ( err_massXicPlus<=0 ) continue;

          if (fWriteXicPlusTree) {
            Int_t lab_XicPlus = -9999;
            if (fIsMC) {
              lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
              if (lab_XicPlus>=0) FillTreeRecXicPlusFromCasc(kfpXicPlus, trackP[itrkBP_trk1], kfpBP_trk1, kfpXiMinus, kfpXiMinus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpLambda, kfpLambda_m, ptrack, ntrack, trackP[itrkBP_trk2], kfpBP_trk2, kfpProton, kfpPionMinus, PV, mcArray, lab_XicPlus);
            }
            if (!fIsMC) {
             FillTreeRecXicPlusFromCasc(kfpXicPlus, trackP[itrkBP_trk1], kfpBP_trk1, kfpXiMinus, kfpXiMinus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpLambda, kfpLambda_m, ptrack, ntrack, trackP[itrkBP_trk2], kfpBP_trk2, kfpProton, kfpPionMinus, PV, mcArray, lab_XicPlus);
            }
          }
          kfpXicPlus.Clear();
          kfpBP_trk2.Clear();
          kfpBP_trk1.Clear();
        } // Loop for second bachelor pion+
      } // Loop for first bachelor pion+
      kfpXiMinus_m.Clear();
      kfpXiMinus.Clear();
      kfpPionOrKaon.Clear();
      kfpLambda_m.Clear();
      kfpLambda.Clear();
    }

    if ( btrack->Charge()>0 ) { // Xi^+

      const KFParticle *vAntiDaughters[2] = {&kfpPionPlus, &kfpAntiProton};

      KFParticle kfpAntiLambda;
      kfpAntiLambda.Construct(vAntiDaughters, 2);
      Float_t massAntiLambda_Rec, err_massAntiLambda;
      kfpAntiLambda.GetMass(massAntiLambda_Rec, err_massAntiLambda);

      // check rapidity of Anti-Lambda
      if ( TMath::Abs(kfpAntiLambda.GetE())<=TMath::Abs(kfpAntiLambda.GetPz()) ) continue;

      // chi2>0 && NDF>0 for selecting Anti-Lambda
      if ( kfpAntiLambda.GetNDF()<=0 || kfpAntiLambda.GetChi2()<=0 ) continue;

      // check cov. of Anti-Lambda
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda) ) continue;

      // err_mass>0 of Anti-Lambda
      if ( err_massAntiLambda<=0 ) continue;

      // Chi2geo cut of Anti-Lambda
      if ( (kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue;

      //************************** calculate l/Δl for Anti-Lambda *************************************
      Double_t dx_AntiLambda = PV.GetX()-kfpAntiLambda.GetX();
      Double_t dy_AntiLambda = PV.GetY()-kfpAntiLambda.GetY();
      Double_t dz_AntiLambda = PV.GetZ()-kfpAntiLambda.GetZ();
      Double_t l_AntiLambda = TMath::Sqrt(dx_AntiLambda*dx_AntiLambda + dy_AntiLambda*dy_AntiLambda + dz_AntiLambda*dz_AntiLambda);
      Double_t dl_AntiLambda = (PV.GetCovariance(0)+kfpAntiLambda.GetCovariance(0))*dx_AntiLambda*dx_AntiLambda + (PV.GetCovariance(2)+kfpAntiLambda.GetCovariance(2))*dy_AntiLambda*dy_AntiLambda + (PV.GetCovariance(5)+kfpAntiLambda.GetCovariance(5))*dz_AntiLambda*dz_AntiLambda + 2*( (PV.GetCovariance(1)+kfpAntiLambda.GetCovariance(1))*dx_AntiLambda*dy_AntiLambda + (PV.GetCovariance(3)+kfpAntiLambda.GetCovariance(3))*dx_AntiLambda*dz_AntiLambda + (PV.GetCovariance(4)+kfpAntiLambda.GetCovariance(4))*dy_AntiLambda*dz_AntiLambda );
      if ( fabs(l_AntiLambda)<1.e-8f ) l_AntiLambda = 1.e-8f;
      dl_AntiLambda = dl_AntiLambda<0. ? 1.e8f : sqrt(dl_AntiLambda)/l_AntiLambda;
      Double_t nErr_l_AntiLambda = l_AntiLambda/dl_AntiLambda;
      //***************************************************************************************

      // l/Deltal cut of Anti-Lambda
      if ( nErr_l_AntiLambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Anti-Lambda
      if ( TMath::Abs(massAntiLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      KFParticle kfpAntiLambda_m = kfpAntiLambda;
      kfpAntiLambda_m.SetNonlinearMassConstraint(massLambda);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda_m) || TMath::Abs(kfpAntiLambda_m.GetE()) <= TMath::Abs(kfpAntiLambda_m.GetPz()) ) continue;

      KFParticle kfpPionOrKaon;
      kfpPionOrKaon = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, 211); // pion+
      KFParticle kfpXiPlus;
      const KFParticle *vXiDs[2] = {&kfpPionOrKaon, &kfpAntiLambda_m};
      kfpXiPlus.Construct(vXiDs, 2);

      // check rapidity of Xi+
      if ( TMath::Abs(kfpXiPlus.GetE())<=TMath::Abs(kfpXiPlus.GetPz()) ) continue;

      // err_massXi > 0
      Float_t massXiPlus_Rec, err_massXiPlus;
      kfpXiPlus.GetMass(massXiPlus_Rec, err_massXiPlus);
      if ( err_massXiPlus<=0 ) continue;

      // chi2>0 && NDF>0
      if ( kfpXiPlus.GetNDF()<=0 || kfpXiPlus.GetChi2()<=0 ) continue;

      // Prefilter
      if ( kfpXiPlus.GetChi2()/kfpXiPlus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

      // check covariance matrix
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus) ) continue;

      // mass window cut of Xi+
      if ( TMath::Abs(massXiPlus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi()) ) continue;

      KFParticle kfpXiPlus_m = kfpXiPlus;
      kfpXiPlus_m.SetNonlinearMassConstraint(massXi);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus_m) || TMath::Abs(kfpXiPlus_m.GetE()) <= TMath::Abs(kfpXiPlus_m.GetPz()) ) continue;

      for (Int_t itrkBP_trk1=0; itrkBP_trk1<flag_trkN-1; itrkBP_trk1++) { // Loop for first bachelor pion-
        for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkN; itrkBP_trk2++) { // Loop for second bachelor pion-

          if ( trackN[itrkBP_trk1]->GetID() == ntrack->GetID() || trackN[itrkBP_trk2]->GetID() == ntrack->GetID() || trackN[itrkBP_trk1]->GetID() == trackN[itrkBP_trk2]->GetID() ) continue;

          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk1]) ) continue;
          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk2]) ) continue;

          KFParticle kfpBP_trk1  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackN[itrkBP_trk1], -211);
          KFParticle kfpBP_trk2   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackN[itrkBP_trk2], -211);

          // reconstruct Anti-XicPlus
          KFParticle kfpAntiXicPlus;
          const KFParticle *vXicPlusDs[3] = {&kfpXiPlus_m, &kfpBP_trk1, &kfpBP_trk2};
          kfpAntiXicPlus.Construct(vXicPlusDs, 3);

          // chi2>0 && NDF>0
          if ( kfpAntiXicPlus.GetNDF()<=0 || kfpAntiXicPlus.GetChi2()<=0 ) continue;

          // Prefilter
          if ( kfpAntiXicPlus.GetChi2()/kfpAntiXicPlus.GetNDF() >= fAnaCuts->GetKFPXicPlus_Chi2geoMax() ) continue;
          if ( kfpAntiXicPlus.GetPt() < fAnaCuts->GetPtMinXicPlus() ) continue;

          // check rapidity of Anti-XicPlus
          if ( TMath::Abs(kfpAntiXicPlus.GetE())<=TMath::Abs(kfpAntiXicPlus.GetPz()) ) continue;

          // check covariance matrix
          if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiXicPlus) ) continue;

          // err_massAntiXicPlus > 0
          Float_t massAntiXicPlus_Rec, err_massAntiXicPlus;
          kfpAntiXicPlus.GetMass(massAntiXicPlus_Rec, err_massAntiXicPlus);
          if ( err_massAntiXicPlus<=0 ) continue;

          if (fWriteXicPlusTree) {
            Int_t lab_AntiXicPlus = -9999.;
            if (fIsMC) {
              lab_AntiXicPlus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
              if (lab_AntiXicPlus>=0) FillTreeRecXicPlusFromCasc(kfpAntiXicPlus, trackN[itrkBP_trk1], kfpBP_trk1, kfpXiPlus, kfpXiPlus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, trackN[itrkBP_trk2], kfpBP_trk2, kfpAntiProton, kfpPionPlus, PV, mcArray, lab_AntiXicPlus);
            }
            if (!fIsMC) {
              FillTreeRecXicPlusFromCasc(kfpAntiXicPlus, trackN[itrkBP_trk1], kfpBP_trk1, kfpXiPlus, kfpXiPlus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, trackN[itrkBP_trk2], kfpBP_trk2, kfpAntiProton, kfpPionPlus, PV, mcArray, lab_AntiXicPlus);
            }
          }
          kfpAntiXicPlus.Clear();
          kfpBP_trk2.Clear();
          kfpBP_trk1.Clear();
        } // Loop for second bachelor pion-
      } // Loop for first bachelor pion-
      kfpXiPlus_m.Clear();
      kfpXiPlus.Clear();
      kfpPionOrKaon.Clear();
      kfpAntiLambda_m.Clear();
      kfpAntiLambda.Clear();
    }

    kfpK0Short.Clear();
    kfpPionPlus.Clear();
    kfpAntiProton.Clear();
    kfpPionMinus.Clear();
    kfpProton.Clear();

  }

  return;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCXicPlus(AliAODTrack *trackProton, AliAODTrack *trackPionMinus_2, AliAODTrack *trackPionMinus_1, AliAODTrack *trackPionPlus_trk1, AliAODTrack *trackPionPlus_trk2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle
  
  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));

  Int_t labelPionMinus_2  = fabs(trackPionMinus_2->GetLabel());
  if (labelPionMinus_2<0) return -1;
  AliAODMCParticle* mcPionMinus_2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_2));

  Int_t labelPionMinus_1  = fabs(trackPionMinus_1->GetLabel());
  if (labelPionMinus_1<0) return -1;
  AliAODMCParticle* mcPionMinus_1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_1));

  Int_t labelPionPlus_trk1  = fabs(trackPionPlus_trk1->GetLabel());
  if (labelPionPlus_trk1<0) return -1;
  AliAODMCParticle* mcPionPlus_trk1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_trk1));

  Int_t labelPionPlus_trk2  = fabs(trackPionPlus_trk2->GetLabel());
  if (labelPionPlus_trk2<0) return -1;
  AliAODMCParticle* mcPionPlus_trk2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_trk2));

  if ( mcProton->GetPdgCode() != 2212 || mcPionMinus_2->GetPdgCode() != -211 || mcPionMinus_1->GetPdgCode() != -211 || mcPionPlus_trk1->GetPdgCode() != 211 || mcPionPlus_trk2->GetPdgCode() != 211 ) return -1; // check pdg

  // === check Lambda ===
  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPionMinus_2->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda
  // ====================

  // === check Xi- ===
  IndexMother[0] = mcMother->GetMother(); // mother of lambda
  IndexMother[1] = mcPionMinus_1->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi-
  // =================

  // === check XicPlus ===
  Int_t IndexMother_XicPlusDau[3];
  IndexMother_XicPlusDau[0] = mcMother->GetMother(); // mother of Xi-
  IndexMother_XicPlusDau[1] = mcPionPlus_trk1->GetMother();
  IndexMother_XicPlusDau[2] = mcPionPlus_trk2->GetMother();
  if ( IndexMother_XicPlusDau[0]<0 || IndexMother_XicPlusDau[1]<0 || IndexMother_XicPlusDau[2]<0 ) return -1; // check mother exist
  if ( IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[1] || IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[2] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[0]));
  if ( mcMother->GetPdgCode() != 4232 || mcMother->GetNDaughters()!=3 ) return -1; // check mother is XicPlus
  // ==================

//  if ( mcMother->IsPrimary() ) return 1;
//  if ( mcMother->IsPhysicalPrimary() ) return 2;
//  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
//  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
//  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCAntiXicPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackPionPlus_2, AliAODTrack *trackPionPlus_1, AliAODTrack *trackPionMinus_trk1, AliAODTrack *trackPionMinus_trk2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle
  
  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));

  Int_t labelPionPlus_2  = fabs(trackPionPlus_2->GetLabel());
  if (labelPionPlus_2<0) return -1;
  AliAODMCParticle* mcPionPlus_2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_2));

  Int_t labelPionPlus_1  = fabs(trackPionPlus_1->GetLabel());
  if (labelPionPlus_1<0) return -1;
  AliAODMCParticle* mcPionPlus_1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_1));

  Int_t labelPionMinus_trk1  = fabs(trackPionMinus_trk1->GetLabel());
  if (labelPionMinus_trk1<0) return -1;
  AliAODMCParticle* mcPionMinus_trk1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_trk1));

  Int_t labelPionMinus_trk2  = fabs(trackPionMinus_trk2->GetLabel());
  if (labelPionMinus_trk2<0) return -1;
  AliAODMCParticle* mcPionMinus_trk2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_trk2));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcPionPlus_2->GetPdgCode() != 211 || mcPionPlus_1->GetPdgCode() != 211 || mcPionMinus_trk1->GetPdgCode() != -211 || mcPionMinus_trk2->GetPdgCode() != -211 ) return -1; // check pdg

  // === check Anti-Lambda ===
  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcPionPlus_2->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Anti-Lambda

  // === check Xi+ ===
  IndexMother[0] = mcMother->GetMother(); // mother of Anti-Lambda
  IndexMother[1] = mcPionPlus_1->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi+
  // =================

  // === check Anti-XicPlus ===
  Int_t IndexMother_XicPlusDau[3];
  IndexMother_XicPlusDau[0] = mcMother->GetMother(); // mother of Xi+
  IndexMother_XicPlusDau[1] = mcPionMinus_trk1->GetMother();
  IndexMother_XicPlusDau[2] = mcPionMinus_trk2->GetMother();
  // =======================
  if ( IndexMother_XicPlusDau[0]<0 || IndexMother_XicPlusDau[1]<0 || IndexMother_XicPlusDau[2]<0 ) return -1; // check mother exist
  if ( IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[1] || IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[2] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[0]));
  if ( mcMother->GetPdgCode() != -4232 || mcMother->GetNDaughters()!=3 ) return -1; // check mother is Anti-XicPlus

//  if ( mcMother->IsPrimary() ) return 1;
//  if ( mcMother->IsPhysicalPrimary() ) return 2;
//  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
//  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
//  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags)
{
  // Select good tracks using fAnaCuts (AliRDHFCuts object)
  if(trkEntries==0) return;

  nSeleTrks=0;                                                                                                 
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = kFALSE;
    
    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);
    
//    if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;
    
//    AliAODTrack *aodt = (AliAODTrack*)track;

/*
    if(!fAnaCuts) continue;
    if(fAnaCuts->SingleTrkCuts(aodt)){
      seleFlags[i]=kTRUE;
      nSeleTrks++;
//      fHistoPiPtRef->Fill(aodt->Pt());
    }
*/
  } // end loop on tracks
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCXiMinus(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));
  Int_t labelPion2  = fabs(trackPion2->GetLabel());
  if (labelPion2<0) return -1;
  AliAODMCParticle* mcPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion2));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 || mcPion2->GetPdgCode() != -211 ) return -1; // check pdg

  Int_t IndexMother[4];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda

  IndexMother[2] = mcMother->GetMother(); // mother of lambda
  IndexMother[3] = mcPion2->GetMother();
  if ( IndexMother[2]<0 || IndexMother[3]<0 ) return -1; // check mother exist
  if ( IndexMother[2] != IndexMother[3] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[2]));
  if ( mcMother->GetPdgCode() != 3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi-

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[2];
  return 1;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCXiPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));
  Int_t labelAntiPion2  = fabs(trackAntiPion2->GetLabel());
  if (labelAntiPion2<0) return -1;
  AliAODMCParticle* mcAntiPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion2));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 || mcAntiPion2->GetPdgCode() != 211 ) return -1; // check pdg

  Int_t IndexMother[4];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Anti-lambda

  IndexMother[2] = mcMother->GetMother(); // mother of lambda
  IndexMother[3] = mcAntiPion2->GetMother();
  if ( IndexMother[2]<0 || IndexMother[3]<0 ) return -1; // check mother exist
  if ( IndexMother[2] != IndexMother[3] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[2]));
  if ( mcMother->GetPdgCode() != -3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi+

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[2];
  return 1;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCLambda(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<=0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<=0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda and only have two daughters

//  AliAODMCParticle* mcMother2 = static_cast<AliAODMCParticle*>(mcArray->At(fabs(mcMother->GetMother())));

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[0];
  return 1;

}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCLambdaFromXi(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<=0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<=0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 ) return -1; // check mother is lambda
  
  Int_t Index_Mother_Xi = mcMother->GetMother();
  if ( Index_Mother_Xi<=0 ) return -1;
  AliAODMCParticle* mcMother_Xi = static_cast<AliAODMCParticle*>(mcArray->At(Index_Mother_Xi));
  if ( mcMother_Xi->GetPdgCode() != 3312 ) return -1;

//  if ( !mcMother_Xi->IsPhysicalPrimary() ) return -1; // check IsPhysicalPrimary()

  return IndexMother[0];

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCAntiLambda(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<=0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<=0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Anti-lambda and only have two daughters

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[0];
  return 1;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCAntiLambdaFromXi(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<=0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<=0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 ) return -1; // check mother is Anti-lambda

  Int_t Index_Mother_Xi = mcMother->GetMother();
  if ( Index_Mother_Xi<=0 ) return -1;
  AliAODMCParticle* mcMother_Xi = static_cast<AliAODMCParticle*>(mcArray->At(Index_Mother_Xi));
  if ( mcMother_Xi->GetPdgCode() != -3312 ) return -1;

//  if ( !mcMother_Xi->IsPhysicalPrimary() ) return -1; // check IsPhysicalPrimary()

  return IndexMother[0];

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCPion(AliAODTrack *track, TClonesArray *mcArray)
{
  Int_t labelPion = fabs(track->GetLabel());
  if (labelPion<=0) return -1;
  AliAODMCParticle* mcPion = static_cast<AliAODMCParticle*>(mcArray->At(labelPion));
  if ( TMath::Abs(mcPion->GetPdgCode()) != 211 ) return -1;

  return labelPion;
}

/*
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToXicPlusMC(TClonesArray *mcArray, Int_t PDGXicPlus, const Int_t nDaughters, const Int_t *daughterIndex, const Int_t *daughterPDG)
{
  ///
  /// Check if this candidate is matched to a MC signal XicPlus
  /// If yes, return Index (>=0) of the AliAODMCParticle
  /// If no, return -1
  ///

  Int_t IndexMom[10] = {0};
  Int_t IndexMother=-1;
  Bool_t pdgUsed[10] = {0};

  // loop on daughter Index
  for(Int_t i=0; i<nDaughters; i++) {
    IndexMom[i]=-1;
    Int_t Index = daughterIndex[i];
    if(Index<0) {
      printf("daughter with negative index %d\n", Index);
      return -1;
    }
    AliAODMCParticle *part = (AliAODMCParticle*)mcArray->At(Index);
    if(!part) { 
      printf("no MC particle\n");
      return -1;
    }

    // check the PDG of the daughter
    Int_t pdgPart = part->GetPdgCode();
    for(Int_t j=0; j<nDaughters; j++) {
      if(!pdgUsed[j] && pdgPart==daughterPDG[j]) {
        pdgUsed[j]=kTRUE;
        break;
      }
    }

    AliAODMCParticle *mother = part;
    while ( mother->GetMother()>=0 ) {
      IndexMother = mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(IndexMother);
      if (!mother) {
        printf("no MC mother particle\n");
        break;
      }
      Int_t pdgMother = mother->GetPdgCode();
      if ( pdgMother==PDGXicPlus ) { // check mother is XicPlus
        IndexMom[i]=IndexMother;
        break;
      } else if( pdgMother>PDGXicPlus || pdgMother<10 ) {
        break;
      }
    }

    if( IndexMom[i]==-1 ) return -1; // mother PDG not ok for this daughter

  } // end loop on daughters

  IndexMother=IndexMom[0];
  for(Int_t i=0; i<nDaughters; i++) {
    // all Index have to be the same and !=-1
    if(IndexMom[i]==-1)          return -1;
    if(IndexMom[i]!=IndexMother) return -1;
    // check that all daughter PDGs are matched
    if(pdgUsed[i]==kFALSE)       return -1;
  }

  return IndexMother;
}
*/

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineEvent()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(3)->GetContainer()->GetName();
  fTree_Event = new TTree(nameoutput, "Event");
  Int_t nVar = 7;
  fVar_Event = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "centrality";
  fVarNames[1]  = "z_vtx_reco";
  fVarNames[2]  = "n_vtx_contributors";
  fVarNames[3]  = "n_tracks";
  fVarNames[4]  = "is_ev_rej";
  fVarNames[5]  = "run_number";
  fVarNames[6]  = "ev_id";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Event->Branch(fVarNames[ivar].Data(), &fVar_Event[ivar], Form("%s/F", fVarNames[ivar].Data()));
  }

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeRecXicPlus()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fTree_XicPlus = new TTree(nameoutput, "XicPlus variables tree");
  Int_t nVar = 47;
  fVar_XicPlus = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "nSigmaTPC_Pi0FromXicPlus"; // TPC nsigma for pion_0 coming from XicPlus
  fVarNames[1]  = "nSigmaTOF_Pi0FromXicPlus"; // TOF nsigma for pion_0 coming from XicPlus
  fVarNames[2]  = "nSigmaTPC_Pi1FromXicPlus"; // TPC nsigma for pion_1 coming from XicPlus
  fVarNames[3]  = "nSigmaTOF_Pi1FromXicPlus"; // TOF nsigma for pion_1 coming from XicPlus
  fVarNames[4]  = "nSigmaTPC_PiFromXi"; // TPC nsigma for pion coming from Xi
  fVarNames[5]  = "nSigmaTOF_PiFromXi"; // TOF nsigma for pion coming from Xi
  fVarNames[6]  = "nSigmaTPC_PiFromLam"; // TPC nsigma for pion coming from Lambda
  fVarNames[7]  = "nSigmaTPC_PrFromLam"; // TPC nsigma for proton coming from Lambda

  fVarNames[8] = "chi2geo_Lam"; // chi2_geometry of Lambda (without mass constraint)
  fVarNames[9]  = "ldl_Lam"; // l/dl of Lambda
  fVarNames[10]  = "chi2topo_LamToPV"; // chi2_topo of Lambda (with mass constraint) to PV

  fVarNames[11] = "chi2geo_Xi"; // chi2_geometry of Xi (with Lambda mass const.)
  fVarNames[12] = "ldl_Xi"; // l/dl of Xi (with Lambda mass const.)
  fVarNames[13] = "chi2topo_XiToPV"; // chi2_topo of Xi (with mass constraint) to PV
  fVarNames[14] = "chi2MassConst_Xi"; // chi2_MassConst of Xi

  fVarNames[15] = "DecayLxy_Lam"; // decay length of Lambda in x-y plane
  fVarNames[16] = "ct_Lam"; // life time of Lambda
  fVarNames[17] = "DecayLxy_Xi"; // decay length of Xi in x-y plane
  fVarNames[18] = "ct_Xi"; // life time of Xi

  fVarNames[19] = "PA_LamToXi"; // pointing angle of Lmabda (pointing back to Xi)
  fVarNames[20] = "PA_LamToPV"; // pointing angle of Lambda (pointing back to PV)
  fVarNames[21] = "PA_XiToPV"; // pointing angle of Xi (pointing back to PV)

  fVarNames[22] = "mass_Lam"; // mass of Lambda (without mass const.)
  fVarNames[23] = "mass_Xi"; // mass of Xi (without mass const.)

  fVarNames[24] = "pt_Pi0FromXicPlus"; // pt of pion_0 coming from XicPlus
  fVarNames[25] = "pt_Pi1FromXicPlus"; // pt of pion_1 coming from XicPlus

  fVarNames[26] = "pt_XicPlus"; // pt of XicPlus
  fVarNames[27] = "rap_XicPlus"; // rapidity of XicPlus
  fVarNames[28] = "mass_XicPlus"; // mass of XicPlus
  fVarNames[29] = "chi2geo_XicPlus"; // chi2_geometry of XicPlus

  fVarNames[30] = "chi2prim_Pi0FromXicPlus"; // chi2_topo of pion_0 to PV
  fVarNames[31] = "chi2prim_Pi1FromXicPlus"; // chi2_topo of pion_1 to PV
  fVarNames[32] = "DCAxy_Pi0FromXicPlusToPV_KF"; // DCA of pion_0 coming from XicPlus in x-y plane
  fVarNames[33] = "DCAxy_Pi1FromXicPlusToPV_KF"; // DCA of pion_1 coming from XicPlus in x-y plane
  fVarNames[34] = "mass_K0S"; // mass of Ks0
  fVarNames[35] = "mass_Gamma"; // mass of e+e-

  fVarNames[36] = "DCAxy_LamDau"; // DCA of Lam's daughters
  fVarNames[37] = "DCAxy_XiDau"; // DCA of Xi's daughters (calculated from KF after Lambda mass constraint)
  fVarNames[38] = "DCAxy_XiToPV"; // DCA of Xi to PV in x-y plane (calculated from KF after Xi mass constraint)
  fVarNames[39] = "DCAxy_PiToPi"; // DCA of pi to pi in x-y plane
  fVarNames[40] = "DCA_PiToPi"; // DCA of pi to pi
  fVarNames[41] = "DCAxy_Pi0ToXi"; // DCA of pion_0 to Xi in x-y plane
  fVarNames[42] = "DCAxy_Pi1ToXi"; // DCA of pion_1 to Xi in x-y plane
  fVarNames[43] = "PA_XicPlusToPV"; // pointing angle of Xic (pointing back to PV)
  fVarNames[44] = "DecayLxy_XicPlus"; // decay length of XicPlus in x-y plane
  fVarNames[45] = "chi2topo_XicPlus"; // chi2_topo of XicPlus to PV

  fVarNames[46] = "Source_XicPlus"; // flag for XicPlus MC truth (“4” prompt, "5" feed-down, “<0” background)

//  fVarNames[26] = "CosThetaStar_PiFromXicPlus"; // CosThetaStar of pion coming from XicPlus
//  fVarNames[27] = "CosThetaStar_Xi"; // CosThetaStar of Xi coming from XicPlus
//  fVarNames[41] = "DCA_XiDau_Cascade"; // DCA of Xi's daughters (calculated from AOD cascade)
//  fVarNames[33] = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
//  fVarNames[34] = "chi2topo_XiToXicPlus"; // chi2_topo of Xi to XicPlus
//  fVarNames[35] = "DecayLxy_XicPlus"; // decay length of XicPlus in x-y plane
//  fVarNames[36] = "ct_XicPlus"; // life time of XicPlus
//  fVarNames[38] = "DCA_LamDau"; // DCA of Lambda's daughters (calculated from AOD cascade)
//  fVarNames[40] = "DCA_XicPlusDau_KF"; // DCA of XicPlus's daughters (calculated from KF after Xi mass constraint)

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_XicPlus->Branch(fVarNames[ivar].Data(), &fVar_XicPlus[ivar], Form("%s/F", fVarNames[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeGenXicPlus()
{
  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fTree_XicPlusMCGen = new TTree(nameoutput,"XicPlus MC variables tree");
  Int_t nVar = 4;
  fVar_XicPlusMCGen = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];
  
  fVarNames[0] = "rap_XicPlus";
  fVarNames[1] = "pt_XicPlus";
  fVarNames[2] = "Source_XicPlus";
  fVarNames[3] = "PDG_XicPlus";

  /*
  fVarNames[ 0]="Centrality";
  fVarNames[ 1]="DecayType";
  fVarNames[ 2]="XicSource";
  fVarNames[ 3]="XicEta";
  fVarNames[ 4]="XicY";
  fVarNames[ 5]="XicPx";
  fVarNames[ 6]="XicPy";
  fVarNames[ 7]="XicPz";
  fVarNames[ 8]="PiPx";
  fVarNames[ 9]="PiPy";
  fVarNames[10]="PiPz";
  fVarNames[11]="CascPx";
  fVarNames[12]="CascPy";
  fVarNames[13]="CascPz";
  fVarNames[14]="XicPdgCode";
  fVarNames[15]="PiPdgCode";
  fVarNames[16]="CascPdgCode";
  fVarNames[17]="RunNumber";
  fVarNames[18]="EvNumber";
  fVarNames[19]="IsPrompt";
  */

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_XicPlusMCGen->Branch(fVarNames[ivar].Data(),&fVar_XicPlusMCGen[ivar],Form("%s/F",fVarNames[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2)
{
  
  Double_t mass1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
  Double_t mass2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
  Double_t E1 = TMath::Sqrt(mass1*mass1 + trk1->P()*trk1->P());
  Double_t E2 = TMath::Sqrt(mass2*mass2 + trk2->P()*trk2->P());
  Double_t mass = TMath::Sqrt( (E1+E2)*(E1+E2) - (trk1->Px()+trk2->Px())*(trk1->Px()+trk2->Px()) - (trk1->Py()+trk2->Py())*(trk1->Py()+trk2->Py()) - (trk1->Pz()+trk2->Pz())*(trk1->Pz()+trk2->Pz()) );

  return mass;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineAnaHist()
{
  // Define analysis histograms
  
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillEventROOTObjects()
{

  for (Int_t i=0; i<7; i++) {
    fVar_Event[i] = 0.;
  }

  Double_t pos[3];
  fpVtx->GetXYZ(pos);

  fVar_Event[1] = pos[2];

  fTree_Event->Fill();

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillTreeRecXicPlusFromCasc(KFParticle kfpXicPlus, AliAODTrack *trackPiFromXicPlus_trk1, KFParticle kfpBP_trk1, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, KFParticle kfpPionOrKaon, AliAODTrack *trackPiFromXiOrKaonFromOmega, KFParticle kfpK0Short, KFParticle kfpGamma, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, AliAODTrack *trackPiFromXicPlus_trk2, KFParticle kfpBP_trk2, KFParticle kfpProtonFromLam, KFParticle kfpPionFromLam, KFParticle PV, TClonesArray *mcArray, Int_t lab_XicPlus)
{

  for (Int_t i=0; i<47; i++) {
    fVar_XicPlus[i] = -9999.;
  }

  Float_t nSigmaTPC_PiFromXicPlus_trk1 = fPID->NumberOfSigmasTPC(trackPiFromXicPlus_trk1,AliPID::kPion);
  Float_t nSigmaTOF_PiFromXicPlus_trk1 = fPID->NumberOfSigmasTOF(trackPiFromXicPlus_trk1,AliPID::kPion);
  Float_t nSigmaTPC_PiFromXicPlus_trk2  = fPID->NumberOfSigmasTPC(trackPiFromXicPlus_trk2,AliPID::kPion);
  Float_t nSigmaTOF_PiFromXicPlus_trk2  = fPID->NumberOfSigmasTOF(trackPiFromXicPlus_trk2,AliPID::kPion);

  Float_t nSigmaTPC_PrFromLam  = fPID->NumberOfSigmasTPC(trkProton,AliPID::kProton);
  Float_t nSigmaTPC_PiFromLam  = fPID->NumberOfSigmasTPC(trkPion,AliPID::kPion);

  Float_t nSigmaTPC_PiFromXi = -9999., nSigmaTOF_PiFromXi = -9999.;
  nSigmaTPC_PiFromXi = fPID->NumberOfSigmasTPC(trackPiFromXiOrKaonFromOmega,AliPID::kPion);
  nSigmaTOF_PiFromXi = fPID->NumberOfSigmasTOF(trackPiFromXiOrKaonFromOmega,AliPID::kPion);

  if ( fabs(nSigmaTPC_PiFromXicPlus_trk1)>=4. || fabs(nSigmaTPC_PiFromXicPlus_trk2)>=4. || fabs(nSigmaTPC_PiFromXi)>=4. || fabs(nSigmaTPC_PrFromLam)>=4. || fabs(nSigmaTPC_PiFromLam)>=4. ) return;

//  AliAODTrack *trk0 = (AliAODTrack*) (v0->GetDaughter(0));

//  Double_t alpha_FirstDaugPos = (v0->MomPosAlongV0() - v0->MomNegAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());

//  if (trk0->Charge()<0) {
//    alpha_FirstDaugPos = (v0->MomNegAlongV0() - v0->MomPosAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());
//  }

  KFParticle kfpXicPlus_PV = kfpXicPlus;
  kfpXicPlus_PV.SetProductionVertex(PV);

  KFParticle kfpXiMinus_XicPlus = kfpXiMinus_m;
  kfpXiMinus_XicPlus.SetProductionVertex(kfpXicPlus);
  KFParticle kfpXiMinus_PV = kfpXiMinus_m;
  kfpXiMinus_PV.SetProductionVertex(PV);

  KFParticle kfpLambda_Xi = kfpLambda_m;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus);
  KFParticle kfpLambda_PV = kfpLambda_m;
  kfpLambda_PV.SetProductionVertex(PV);
//  KFParticle kfpXiMinus_pv = kfpXiMinus;
//  kfpXiMinus_pv.SetProductionVertex(PV);

  KFParticle kfpBP_XicPlus = kfpBP_trk1;
  kfpBP_XicPlus.SetProductionVertex(kfpXicPlus);

  kfpBP_trk1.GetDistanceFromVertexXY(PV);
  kfpBP_trk2.GetDistanceFromVertexXY(PV);

  //************************** calculate l/Δl for Lambda *************************************
  Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
  Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
  Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
  Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
  Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
  if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
  dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
  if ( dl_Lambda<=0 ) return;
  //***************************************************************************************
  //************************** calculate l/Δl for Xi- *************************************
  Double_t dx_Xi = PV.GetX()-kfpXiMinus.GetX();
  Double_t dy_Xi = PV.GetY()-kfpXiMinus.GetY();
  Double_t dz_Xi = PV.GetZ()-kfpXiMinus.GetZ();
  Double_t l_Xi = TMath::Sqrt(dx_Xi*dx_Xi + dy_Xi*dy_Xi + dz_Xi*dz_Xi);
  Double_t dl_Xi = (PV.GetCovariance(0)+kfpXiMinus.GetCovariance(0))*dx_Xi*dx_Xi + (PV.GetCovariance(2)+kfpXiMinus.GetCovariance(2))*dy_Xi*dy_Xi + (PV.GetCovariance(5)+kfpXiMinus.GetCovariance(5))*dz_Xi*dz_Xi + 2*( (PV.GetCovariance(1)+kfpXiMinus.GetCovariance(1))*dx_Xi*dy_Xi + (PV.GetCovariance(3)+kfpXiMinus.GetCovariance(3))*dx_Xi*dz_Xi + (PV.GetCovariance(4)+kfpXiMinus.GetCovariance(4))*dy_Xi*dz_Xi );
  if ( fabs(l_Xi)<1.e-8f ) l_Xi = 1.e-8f;
  dl_Xi = dl_Xi<0. ? 1.e8f : sqrt(dl_Xi)/l_Xi;
  if ( dl_Xi<=0 ) return;
  //***************************************************************************************

  if ( kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF() <= fAnaCuts->GetKFPLam_Chi2topoMin() ) return;
  if ( kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF() >= fAnaCuts->GetKFPXi_Chi2topoMax() ) return;
  if ( l_Xi/dl_Xi <= fAnaCuts->GetKFPXi_lDeltalMin() ) return;

  const Float_t PDGmassXicPlus = TDatabasePDG::Instance()->GetParticle(4232)->Mass();
  Float_t mass_XicPlus_PV, err_mass_XicPlus_PV;
  kfpXicPlus_PV.GetMass(mass_XicPlus_PV, err_mass_XicPlus_PV);
  fVar_XicPlus[28] = mass_XicPlus_PV; // mass of XicPlus

  if ( fabs(mass_XicPlus_PV-PDGmassXicPlus) > fAnaCuts->GetProdMassTolXicPlus() ) return;


  fVar_XicPlus[0]  = nSigmaTPC_PiFromXicPlus_trk1; // TPC nsigma for pion_0 coming from XicPlus
  fVar_XicPlus[1]  = nSigmaTOF_PiFromXicPlus_trk1; // TOF nsigma for pion_0 coming from XicPlus
  fVar_XicPlus[2]  = nSigmaTPC_PiFromXicPlus_trk2; // TPC nsigma for pion_1 coming from XicPlus
  fVar_XicPlus[3]  = nSigmaTOF_PiFromXicPlus_trk2; // TOF nsigma for pion_1 coming from XicPlus
  fVar_XicPlus[4]  = nSigmaTPC_PiFromXi; // TPC nsigma for pion coming from Xi
  fVar_XicPlus[5]  = nSigmaTOF_PiFromXi; // TOF nsigma for pion coming from Xi
  fVar_XicPlus[6]  = nSigmaTPC_PiFromLam; // TPC nsigma for pion coming from Lambda
  fVar_XicPlus[7]  = nSigmaTPC_PrFromLam; // TPC nsigma for proton coming from Lambda

  fVar_XicPlus[8] = kfpLambda.GetChi2()/kfpLambda.GetNDF(); // chi2geo_Lam
  fVar_XicPlus[9] = l_Lambda/dl_Lambda; // l/dl of Lambda
  fVar_XicPlus[10] = kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF(); // chi2_topo of Lambda (with mass constraint) to PV

  fVar_XicPlus[11] = kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF(); // chi2_geometry of Xi (with Lambda mass const.)
  fVar_XicPlus[12] = l_Xi/dl_Xi; // l/dl of Xi (with Lambda mass const.)
  fVar_XicPlus[13] = kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF(); // chi2_topo of Xi (with mass constraint) to PV
  fVar_XicPlus[14] = kfpXiMinus_m.GetChi2()/kfpXiMinus_m.GetNDF(); // chi2_MassConst of Xi

  Float_t DecayLxy_Lam, err_DecayLxy_Lam;
  kfpLambda_Xi.GetDecayLengthXY(DecayLxy_Lam, err_DecayLxy_Lam);
  fVar_XicPlus[15] = DecayLxy_Lam; // decay length of Lambda in x-y plane
  Float_t ct_Lam=0., err_ct_Lam=0.;
  kfpLambda_Xi.GetLifeTime(ct_Lam, err_ct_Lam);
  fVar_XicPlus[16] = ct_Lam; // life time of Lambda

  Float_t DecayLxy_Xi, err_DecayLxy_Xi;
  kfpXiMinus_PV.GetDecayLengthXY(DecayLxy_Xi, err_DecayLxy_Xi);
  fVar_XicPlus[17] = DecayLxy_Xi; // decay length of Xi in x-y plane
  Float_t ct_Xi=0., err_ct_Xi=0.;
  kfpXiMinus_PV.GetLifeTime(ct_Xi, err_ct_Xi);
  fVar_XicPlus[18] = ct_Xi; // life time of Xi

  // calculate CosPointingAngle
  Double_t cosPA_v0toXi = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLambda_m, kfpXiMinus);
  Double_t cosPA_v0toPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLambda_m, PV);
  Double_t cosPA_XiToPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXiMinus_m, PV);
  fVar_XicPlus[19] = TMath::ACos(cosPA_v0toXi); // pointing angle of Lmabda (pointing back to Xi)
  fVar_XicPlus[20] = TMath::ACos(cosPA_v0toPV); // pointing angle of Lambda (pointing back to PV)
  fVar_XicPlus[21] = TMath::ACos(cosPA_XiToPV); // pointing angle of Xi (pointing back to PV)

  Float_t mass_Lam, err_mass_Lam;
  kfpLambda.GetMass(mass_Lam, err_mass_Lam);
  fVar_XicPlus[22] = mass_Lam; // mass of Lambda (without mass const.)

  Float_t mass_Xi, err_mass_Xi;
  kfpXiMinus.GetMass(mass_Xi, err_mass_Xi);
  fVar_XicPlus[23] = mass_Xi; // mass of Xi (without mass const.)

  fVar_XicPlus[24] = trackPiFromXicPlus_trk1->Pt(); // pt of pion_0 coming from XicPlus
  fVar_XicPlus[25] = trackPiFromXicPlus_trk2->Pt(); // pt of pion_1 coming from XicPlus

  fVar_XicPlus[26] = kfpXicPlus_PV.GetPt(); // pt of XicPlus
  if ( TMath::Abs(kfpXicPlus_PV.GetE())>TMath::Abs(kfpXicPlus_PV.GetPz()) ) {
    fVar_XicPlus[27] = kfpXicPlus_PV.GetRapidity(); // rapidity of XicPlus
  }

  fVar_XicPlus[29] = kfpXicPlus.GetChi2()/kfpXicPlus.GetNDF(); // chi2_geometry of XicPlus

  // --- chi2_prim of Pion to PV ---
  KFParticle kfpBP_trk1_PV = kfpBP_trk1;
  kfpBP_trk1_PV.SetProductionVertex(PV);
  fVar_XicPlus[30] = kfpBP_trk1_PV.GetChi2()/kfpBP_trk1_PV.GetNDF(); // chi2_topo of pion_0 to PV
  KFParticle kfpBP_trk2_PV = kfpBP_trk2;
  kfpBP_trk2_PV.SetProductionVertex(PV);
  fVar_XicPlus[31] = kfpBP_trk2_PV.GetChi2()/kfpBP_trk2_PV.GetNDF(); // chi2_topo of pion_1 to PV
  // -------------------------------
  
  // --- DCA of Pion to PV ---
  fVar_XicPlus[32] = kfpBP_trk1.GetDistanceFromVertexXY(PV); // DCA of pion_0 coming from XicPlus in x-y plane
  fVar_XicPlus[33] = kfpBP_trk2.GetDistanceFromVertexXY(PV); // DCA of pion_1 coming from XicPlus in x-y plane
  // -------------------------

  Float_t massK0S_Rec, err_massK0S;
  kfpK0Short.GetMass(massK0S_Rec, err_massK0S);
  fVar_XicPlus[34] = massK0S_Rec; // mass of Ks0
  Float_t massGamma_Rec, err_massGamma;
  kfpGamma.GetMass(massGamma_Rec, err_massGamma);
  fVar_XicPlus[35] = massGamma_Rec; // mass of e+e-

  fVar_XicPlus[36] = kfpPionFromLam.GetDistanceFromParticleXY(kfpProtonFromLam); // DCA of Lam's daughters
  fVar_XicPlus[37] = kfpPionOrKaon.GetDistanceFromParticleXY(kfpLambda_m); // DCA of Xi's daughters (calculated from KF after Lambda mass constraint)
  fVar_XicPlus[38] = kfpXiMinus_m.GetDistanceFromVertexXY(PV); // DCA of Xi to PV in x-y plane (calculated from KF after Xi mass constraint)
  fVar_XicPlus[39] = kfpBP_trk1.GetDistanceFromParticleXY(kfpBP_trk2); // DCA of pi to pi in x-y plane
  fVar_XicPlus[40] = kfpBP_trk1.GetDistanceFromParticle(kfpBP_trk2); // DCA of pi to pi
  fVar_XicPlus[41] = kfpBP_trk1.GetDistanceFromParticleXY(kfpXiMinus_m); // DCA of pi_trk1 to Xi in x-y plane
  fVar_XicPlus[42] = kfpBP_trk2.GetDistanceFromParticleXY(kfpXiMinus_m); // DCA of pi_trk2 to Xi in x-y plane

  // --- CosPointingAngle ---
  Double_t cosPA_XicPlusToPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXicPlus, PV);
  fVar_XicPlus[43] = TMath::ACos(cosPA_XicPlusToPV); // pointing angle of XicPlus (pointing back to PV)
  // -------------------------
  Float_t DecayLxy_XicPlus, err_DecayLxy_XicPlus;
  kfpXicPlus_PV.GetDecayLengthXY(DecayLxy_XicPlus, err_DecayLxy_XicPlus);
  fVar_XicPlus[44] = DecayLxy_XicPlus; // decay length of XicPlus in x-y plane
  fVar_XicPlus[45] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF(); // chi2_topo of XicPlus to PV

//  fVar_XicPlus[26] = AliVertexingHFUtils::CosThetaStarFromKF(0, 4132, 211, 3312, kfpXicPlus, kfpBP_XicPlus, kfpXiMinus_XicPlus);
//  fVar_XicPlus[27] = AliVertexingHFUtils::CosThetaStarFromKF(1, 4132, 211, 3312, kfpXicPlus, kfpBP_XicPlus, kfpXiMinus_XicPlus);


//  fVar_XicPlus[33] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();
//  fVar_XicPlus[34] = kfpXiMinus_XicPlus.GetChi2()/kfpXiMinus_XicPlus.GetNDF();

//  Float_t ct_XicPlus=0., err_ct_XicPlus=0.;
//  kfpXicPlus_PV.GetLifeTime(ct_XicPlus, err_ct_XicPlus);
//  fVar_XicPlus[36] = ct_XicPlus;
//  fVar_XicPlus[38] = casc->DcaV0Daughters(); // DCA_LamDau
//  fVar_XicPlus[40] = kfpBP_trk1.GetDistanceFromParticle(kfpXiMinus_m); // DCA_XicPlusDau_KF

  if (fIsMC) {
    fVar_XicPlus[46] = lab_XicPlus;
    // === weight ===
    /*
    if (lab_XicPlus>0) {
      Int_t labelPion_trk1 = fabs(trackPiFromXicPlus_trk1->GetLabel());
      AliAODMCParticle* mcPion_trk1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion_trk1));
      AliAODMCParticle* mcXicPlus  = static_cast<AliAODMCParticle*>(mcArray->At(mcPion_trk1->GetMother()));
    }
    */
  }

  // pt(XicPlus)>=1
  if (fVar_XicPlus[26]>0.9999) fTree_XicPlus->Fill();



//  fVar_XicPlus[10] = casc->DcaXiDaughters(); // DCA_XiDau
//  fVar_XicPlus[32] = AliVertexingHFUtils::CosThetaStarFromKF(0, 4132, 211, 3312, kfpXicPlus, kfpBP_trk1, kfpXiMinus_m);
//  fVar_XicPlus[33] = AliVertexingHFUtils::CosThetaStarFromKF(1, 4132, 211, 3312, kfpXicPlus, kfpBP_trk1, kfpXiMinus_m);

//  fVar_XicPlus[3]  = TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF());
//  fVar_XicPlus[7]  = kfpLambda.GetPz();
//  fVar_XicPlus[8]  = kfpLambda.GetX();
//  fVar_XicPlus[9]  = kfpLambda.GetY();
//  fVar_XicPlus[10] = kfpLambda.GetZ();
//  fVar_XicPlus[8] = kfpLambda_pv.GetChi2()/kfpLambda_pv.GetNDF();
//  fVar_XicPlus[9] = TMath::Prob(kfpLambda_pv.GetChi2(), kfpLambda_pv.GetNDF());
//  fVar_XicPlus[15] = alpha_FirstDaugPos;
//  fVar_XicPlus[16] = v0->PtArmV0();

//  fVar_XicPlus[17] = trackPiFromXi->Charge();
//  fVar_XicPlus[20] = TMath::Prob(kfpXiMinus.GetChi2(), kfpXiMinus.GetNDF());
//  fVar_XicPlus[21] = kfpXiMinus_pv.GetChi2()/kfpXiMinus_pv.GetNDF();
//  fVar_XicPlus[22] = TMath::Prob(kfpXiMinus_pv.GetChi2(), kfpXiMinus_pv.GetNDF());
//  Float_t mass_Xi_PV, err_mass_Xi_PV;
//  kfpXiMinus_pv.GetMass(mass_Xi_PV, err_mass_Xi_PV);
//  fVar_XicPlus[23] = mass_Xi_PV;
//  if ( TMath::Abs(kfpXiMinus_pv.GetE())>TMath::Abs(kfpXiMinus_pv.GetPz()) ) {
//    fVar_XicPlus[24] = kfpXiMinus_pv.GetRapidity();
//  }
//  fVar_XicPlus[25] = kfpXiMinus_pv.GetPt();
//  fVar_XicPlus[26] = kfpXiMinus_pv.GetPz();
//  fVar_XicPlus[27] = kfpXiMinus_pv.GetX();
//  fVar_XicPlus[28] = kfpXiMinus_pv.GetY();
//  fVar_XicPlus[29] = kfpXiMinus_pv.GetZ();
//  fVar_XicPlus[33] = kfpXiMinus.GetPz();
//  fVar_XicPlus[34] = kfpXiMinus.GetX();
//  fVar_XicPlus[35] = kfpXiMinus.GetY();
//  fVar_XicPlus[36] = kfpXiMinus.GetZ();

//  fVar_XicPlus[39] = trackPiFromXicPlus_trk1->Charge();
//  fVar_XicPlus[41] = trackPiFromXicPlus_trk1->Pz();

//  fVar_XicPlus[20] = kfpXicPlus.GetChi2()/kfpXicPlus.GetNDF();
//  fVar_XicPlus[45] = TMath::Prob(kfpXicPlus.GetChi2(), kfpXicPlus.GetNDF());
//  fVar_XicPlus[21] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF();
//  fVar_XicPlus[47] = TMath::Prob(kfpXicPlus_PV.GetChi2(), kfpXicPlus_PV.GetNDF());
//  fVar_XicPlus[51] = kfpXicPlus_PV.GetPz();
//  fVar_XicPlus[52] = kfpXicPlus_PV.GetX();
//  fVar_XicPlus[53] = kfpXicPlus_PV.GetY();
//  fVar_XicPlus[54] = kfpXicPlus_PV.GetZ();
//  Float_t mass_Xic, err_mass_Xic;
//  kfpXicPlus.GetMass(mass_Xic, err_mass_Xic);
//  fVar_XicPlus[55] = mass_Xic;
//  fVar_XicPlus[56] = kfpXicPlus.GetRapidity();
//  fVar_XicPlus[57] = kfpXicPlus.GetPt();
//  fVar_XicPlus[58] = kfpXicPlus.GetPz();
//  fVar_XicPlus[59] = kfpXicPlus.GetX();
//  fVar_XicPlus[60] = kfpXicPlus.GetY();
//  fVar_XicPlus[61] = kfpXicPlus.GetZ();

//  fVar_XicPlus[27] = cosPA_XiToXic;


//  Float_t CT_Lam, err_CT_Lam;
//  Float_t CT_Xi, err_CT_Xi;
//  Float_t CT_XicPlus, err_CT_XicPlus;
//  kfpLambda_Xi.GetLifeTime(CT_Lam, err_CT_Lam);
//  kfpXiMinus_pv.GetLifeTime(CT_Xi, err_CT_Xi);
//  kfpXicPlus_PV.GetLifeTime(CT_XicPlus, err_CT_XicPlus);
//  fVar_XicPlus[38] = CT_Lam;
//  fVar_XicPlus[39] = CT_Xi;
//  fVar_XicPlus[40] = CT_XicPlus;

//  if (fIsMC) {
//    fVar_XicPlus[78] = lab_XicPlus;
//  }

  return;
}
