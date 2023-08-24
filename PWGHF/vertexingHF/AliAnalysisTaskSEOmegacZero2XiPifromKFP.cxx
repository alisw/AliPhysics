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
// Author: Federica Zanone (1)
// (1) Ruprecht-Karls-Universität Heidelberg
// E-mail: federica.zanone@cern.ch
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
#include "AliAnalysisTaskSEOmegacZero2XiPifromKFP.h"
#include "AliPIDResponse.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCutsKFP.h"
#include "AliVertexingHFUtils.h"

//includes added to play with KFParticle
#ifndef HomogeneousField
#define HomogeneousField
#endif

using std::cout;
using std::endl;

class AliAnalysisTaskSEOmegacZero2XiPifromKFP;    //your analysis class

ClassImp(AliAnalysisTaskSEOmegacZero2XiPifromKFP) //classimp: necessary for root

AliAnalysisTaskSEOmegacZero2XiPifromKFP::AliAnalysisTaskSEOmegacZero2XiPifromKFP() :
  AliAnalysisTaskSE(),
  fIsMC(kFALSE),
  fPID(0),
  fAnaCuts(0),
  fpVtx(0),
  fBzkG(0),
  fOutputList(0),
  fListCuts(0),
  fTree_Omegac0(0),
  fVar_Omegac0(0),
  fCounter(0),
  fHistEvents(0),
  fHistMult(0),
  fHistCheckKF(0),
  fEvCount(0),
  fHistMCPdgCode(0),
  fHistMCOmegacDauNumber(0),
  fHistMCOmegacPt(0),
  fHistMCOmegacEta(0),
  fHistMCOmegacSign(0),
  fHistMCOmegacCTau(0),
  fHistMCOmegacM(0),
  fHistMCDecayChain(0),
  fHistMCOrigin(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSEOmegacZero2XiPifromKFP::AliAnalysisTaskSEOmegacZero2XiPifromKFP(const char* name, AliRDHFCutsKFP* cuts, Bool_t isMC) :
  AliAnalysisTaskSE(name),
  fIsMC(isMC),
  fPID(0),
  fAnaCuts(cuts),
  fpVtx(0),
  fBzkG(0),
  fOutputList(0),
  fListCuts(0),
  fTree_Omegac0(0),
  fVar_Omegac0(0),
  fCounter(0),
  fHistEvents(0),
  fHistMult(0),
  fHistCheckKF(0),
  fEvCount(0),
  fHistMCPdgCode(0),
  fHistMCOmegacDauNumber(0),
  fHistMCOmegacPt(0),
  fHistMCOmegacEta(0),
  fHistMCOmegacSign(0),
  fHistMCOmegacCTau(0),
  fHistMCOmegacM(0),
  fHistMCDecayChain(0),
  fHistMCOrigin(0)
{
    // constructor
  DefineInput(0, TChain::Class());  // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
  DefineOutput(1, TList::Class());  // define the ouptut of the analysis: in this case it's a list of histograms
  //list of cuts by CutObject           // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
                                        //meaning postdata: PostData will notify the client tasks of the data container that that the data pointer has changed compared to the previous post
  DefineOutput(2, AliNormalizationCounter::Class());
  DefineOutput(3, TTree::Class()); // Omegac0
  DefineOutput(4, TList::Class()); // Event hist
  // the number defines a slot (see AddTask)
}
//_____________________________________________________________________________
AliAnalysisTaskSEOmegacZero2XiPifromKFP::~AliAnalysisTaskSEOmegacZero2XiPifromKFP()
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

    if (fTree_Omegac0) {
      delete fTree_Omegac0;
      fTree_Omegac0 = 0;
    }

    if (fVar_Omegac0) {
      delete fVar_Omegac0;
      fVar_Omegac0 = 0;
    }

    if (fCounter) {
      delete fCounter;
      fCounter = 0;
    }

}

//______________________________________________________________________________
//______________________________________________________________________________

void AliAnalysisTaskSEOmegacZero2XiPifromKFP::Init()
{
  // Initialization (of the analysis task - to use the cut objct, too)

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner(); //By calling TList::SetOwner(true), we transfer ownership of all memory allocated by the list items to the TList itself, this means that in our destructor, we can simply call delete list to delete all our list items, rather than calling delete for all items individually
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsKFP(*fAnaCuts)); //AliRDHFCutsKFP is the CutObject class
  PostData(1, fListCuts); //PostData will notify the client tasks of the data container that that the data pointer has changed compared to the previous post

  return;

}

//______________________________________________________________________________
//______________________________________________________________________________

void AliAnalysisTaskSEOmegacZero2XiPifromKFP::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //

  fEvCount=0;

  fOutputList = new TList();          // this is a list which will contain all of your histograms, at the end of the analysis, the contents of this list are written to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested

  //histogram for general checks on number of events
  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0.5, 11.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events"); //all the events
  fHistEvents->GetXaxis()->SetBinLabel(2,"AliAODVertex exists"); //events with magnetic field B and primary vertex PV
  fHistEvents->GetXaxis()->SetBinLabel(3,"TriggerOK"); //events with B, PV, trigger selected, physicsselection, PV reconstructed, |z|<10
  fHistEvents->GetXaxis()->SetBinLabel(4,"IsEventSelected"); //events con B, PV, physicsselection ok, pileup ok, vertex reconstructed, trigger class selected
  fHistEvents->GetXaxis()->SetBinLabel(5,"V0 exists"); //events with B, PV, at least one V0, physicsselection ok, pileup ok, vertex reconstructed, trigger class selected
  fHistEvents->GetXaxis()->SetBinLabel(6,"Cascade exists"); //events with B, PV, at least one cascade, physicsselection ok, pileup ok, vertex reconstructed, trigger class selected
  fHistEvents->GetXaxis()->SetBinLabel(7,"TriggerMB"); //events with B, PV, trigger MB, PhysicsSelection
  fHistEvents->GetXaxis()->SetBinLabel(8,"TriggerkINT7zCut"); //events with B, PV, trigger kINT7, physicsselection, |z|<10
  fHistEvents->GetXaxis()->SetBinLabel(9,"TriggerkINT7"); //events with B, PV, trigger kINT7, physicsselection
  fHistEvents->GetXaxis()->SetBinLabel(10,"kINT7&B"); //events with B, trigger kINT7, (for rejection factor studies)
  fHistEvents->GetXaxis()->SetBinLabel(11,"OnlykINT7"); //events with trigger kINT7, (for rejection factor studies)


  fOutputList->Add(fHistEvents); // don't forget to add it to the list! the list will be written to file

  //histogram multiplicity distribution
  fHistMult = new TH1F("fHistMult", "fHistMult", 2500, 0., 5000);
  fHistMult->GetXaxis()->SetTitle("N_{tracks}");

  fOutputList->Add(fHistMult);

  //histogram check on KF failure (E<p_z)
  fHistCheckKF = new TH1F("fHistCheckKF", "fHistCheckKF", 4, 0.5, 4.5);
  fHistCheckKF->GetXaxis()->SetBinLabel(1,"KFOkayBeforeMC"); //KF does not fail before setting mass constraint
  fHistCheckKF->GetXaxis()->SetBinLabel(2,"KFFailsBeforeMC"); //KF fails before setting mass constraint
  fHistCheckKF->GetXaxis()->SetBinLabel(3,"KFOkayAfterMC"); //KF does not fail after setting mass constraint
  fHistCheckKF->GetXaxis()->SetBinLabel(4,"KFFailsAfterMC"); //KF fails after setting mass constraint

  fOutputList->Add(fHistCheckKF);

  fHistMCPdgCode= new TH1F("fHistMCPdgCode","fHistMCPdgCode",10000,-5000,5000);
  fOutputList->Add(fHistMCPdgCode);
  fHistMCOmegacDauNumber= new TH1F("fHistMCOmegacDauNumber","fHistMCOmegacDauNumber",10,0,10);
  fOutputList->Add(fHistMCOmegacDauNumber);
  fHistMCOmegacPt= new TH1F("fHistMCOmegacPt","fHistMCOmegacPt",20,0,20);
  fOutputList->Add(fHistMCOmegacPt);
  fHistMCOmegacEta= new TH1F("fHistMCOmegacEta","fHistMCOmegacEta",20,-1,1);
  fOutputList->Add(fHistMCOmegacEta);
  fHistMCOmegacSign= new TH1F("fHistMCOmegacSign","fHistMCOmegacSign",10,-5,5);
  fOutputList->Add(fHistMCOmegacSign);
  fHistMCOmegacCTau= new TH1F("fHistMCOmegacCTau","fHistMCOmegacCTau",3000,0.,1500.);
  fOutputList->Add(fHistMCOmegacCTau);
  fHistMCOmegacM= new TH1F("fHistMCOmegacM","fHistMCOmegacM",500,2.5,2.8);
  fOutputList->Add(fHistMCOmegacM);
  fHistMCDecayChain= new TH1F("fHistMCDecayChain","fHistMCDecayChain",10,-5,5);
  fOutputList->Add(fHistMCDecayChain);
  fHistMCOrigin= new TH1F("fHistMCOrigin","fHistMCOrigin",20,-10,10);
  fOutputList->Add(fHistMCOrigin);

  //counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
  if(cont) normName = (TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(2, fCounter); // postdata will notify the analysis manager of changes/updates to the data pointer (it has changed compared to the previous post)

  //Omegac0 tree
  DefineTreeOmegac0();
  PostData(3, fTree_Omegac0);

  //histogram list for general checks
  PostData(4, fOutputList); //fOutputList object - the manager will in the end take care of writing your output to file, so it needs to know what's in the output

  return;

}

//______________________________________________________________________________
//______________________________________________________________________________

void AliAnalysisTaskSEOmegacZero2XiPifromKFP::UserExec(Option_t *)
{
  // user exec: function that is`called once for each event. The manager will take care of reading the events from file, and with the static function InputEvent() you have access to the current event.
  // Once you return from the UserExec function, the manager will retrieve the next event from the chain.

  if (!fInputEvent) { // if the event is empty (getting it failed) skip this event
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* AODEvent = dynamic_cast<AliAODEvent*>(fInputEvent);    // get an event (called AODEvent) from the input file (there's another event format (ESD) which works in a similar way, but is more cpu/memory unfriendly -> for now, we'll stick with aod)

  fHistMult->Fill(AODEvent->GetNumberOfTracks());
  fHistEvents->Fill(1);

  //----------------------------------------------------------------------------
  //verify CutObject exists
  //----------------------------------------------------------------------------
  if(!fAnaCuts) return;

  //----------------------------------------------------------------------------
  // First check if the event has magnetic field and proper vertex - selecting events with B e PV (and check on kINT7 trigger for rejection factor studies)
  //----------------------------------------------------------------------------
  Bool_t IsINT7RejFac = (  ( (AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7 )==(AliVEvent::kINT7);

  if(IsINT7RejFac){
    fHistEvents->Fill(11);
  }

  fBzkG = (Double_t)AODEvent->GetMagneticField();
  if (TMath::Abs(fBzkG)<0.001) return;
  KFParticle::SetField(fBzkG);

  if(IsINT7RejFac){
    fHistEvents->Fill(10);
  }

  fpVtx = (AliAODVertex*)AODEvent->GetPrimaryVertex();
  if (!fpVtx) return;
  fHistEvents->Fill(2);

  fCounter->StoreEvent(AODEvent,fAnaCuts,fIsMC);

  //----------------------------------------------------------------------------
  // MC analysis setting
  //----------------------------------------------------------------------------
  if(fIsMC){

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(AODEvent->FindListObject(AliAODMCParticle::StdBranchName())); //findlistobject(nameobject)=return the pointer to the object with the given name.
    if ( !mcArray ) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }

    // load MC header
    mcHeader = (AliAODMCHeader*)AODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()); //getlist=returns a list of AOD objects
    if ( !mcHeader ) {
      AliError("AliAnalysisTaskSEOmegacZero2XiPifromKFP::UserExec: MC header branch not found!\n");
      return;
    }

    Double_t zMCvtx = mcHeader->GetVtxZ();
    if ( TMath::Abs(zMCvtx) > fAnaCuts->GetMaxVtxZ() ) {
      AliDebug(2,Form("Event rejected: fabs(zVtxMC)=%f > fAnaCuts->GetMaxVtxZ()=%f", zMCvtx, fAnaCuts->GetMaxVtxZ()));
      return;
    }
    if ((TMath::Abs(zMCvtx) < fAnaCuts->GetMaxVtxZ()) && (!fAnaCuts->IsEventRejectedDuePhysicsSelection()) && (!fAnaCuts->IsEventRejectedDueToTrigger())) {
      Bool_t selevt = MakeMCCheck(mcArray);
      if(!selevt) return;
    }
  }

  //----------------------------------------------------------------------------
  // Event selection - selecting events kINT7, PhysicsSelection ok, PilueUp ok, con PV reconstructed and trigger selection ok
  //----------------------------------------------------------------------------

  Bool_t IsTriggerNotOK = fAnaCuts->IsEventRejectedDueToTrigger();
  Bool_t IsPhysSelNotOK = fAnaCuts->IsEventRejectedDuePhysicsSelection();
  Bool_t IsNoVertex = fAnaCuts->IsEventRejectedDueToNotRecoVertex();

  if( !IsTriggerNotOK && !IsPhysSelNotOK && !IsNoVertex && fabs(fpVtx->GetZ())<fAnaCuts->GetMaxVtxZ() ) fHistEvents->Fill(3);

  Bool_t IsMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  if(IsMB && !IsPhysSelNotOK){
    fHistEvents->Fill(7);
  }

  Bool_t IsINT7 = (  ( (AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7 )==(AliVEvent::kINT7);
  if(IsINT7 && !IsPhysSelNotOK){
    fHistEvents->Fill(9);
  }

  if(!IsPhysSelNotOK && fabs(fpVtx->GetZ())<fAnaCuts->GetMaxVtxZ() && IsINT7){
    fHistEvents->Fill(8);
  }


  //see AliRDHFCuts::IsEventSelected(AliVEvent *event)
  Bool_t IsEventSelected = fAnaCuts->IsEventSelected(AODEvent);
  if(!IsEventSelected) {
    return;
  }
  fHistEvents->Fill(4);

  //----------------------------------------------------------------------------
  // Check if the event has v0 candidate
  //----------------------------------------------------------------------------
  Int_t num_v0 = AODEvent->GetNumberOfV0s();
  if (num_v0>0) fHistEvents->Fill(5);

  //----------------------------------------------------------------------------
  // Check if the event has cascade candidate - selecting events with at least one cascade
  //----------------------------------------------------------------------------
  Int_t num_casc = AODEvent->GetNumberOfCascades();
  if (num_casc<=0) return;
  fHistEvents->Fill(6);

  //----------------------------------------------------------------------------
  // set primary vertex - selecting events with |z|<10cm
  //----------------------------------------------------------------------------
  KFPVertex pVertex;
  Double_t pos[3],cov[6];
  fpVtx->GetXYZ(pos);
  if ( fabs(pos[2])>fAnaCuts->GetMaxVtxZ() ) return; // vertex cut on z-axis direction (cut implemented in CutObject |z|<10cm)
  fpVtx->GetCovarianceMatrix(cov);
  pVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
  Float_t covF[6];
  for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }
  pVertex.SetCovarianceMatrix(covF);
  pVertex.SetChi2(fpVtx->GetChi2());
  pVertex.SetNDF(fpVtx->GetNDF());
  pVertex.SetNContributors(fpVtx->GetNContributors());

  KFParticle PV(pVertex);

  fPID = fInputHandler->GetPIDResponse();

  //----------------------------------------------------------------------------
  // Main analysis done in this function
  //----------------------------------------------------------------------------

  MakeAnaOmegacZero(AODEvent, PV);

  fEvCount=fEvCount+1;

  /*if(fEvCount==100){
    fEvCount=0;
  }*/

  // stream the results the analysis of this event to the output manager which will take care of writing it to a file
  PostData(2, fCounter);
  PostData(3, fTree_Omegac0);
  PostData(4, fOutputList);

  return;
}

//______________________________________________________________________________
//______________________________________________________________________________

void AliAnalysisTaskSEOmegacZero2XiPifromKFP::Terminate(Option_t *)
{
    // terminate: function called at the END of the analysis (when all events are processed)
    return;
}








//------------------------------------------------------------------------------
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// NEW FUNCTIONS ANALYSISTASK
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//------------------------------------------------------------------------------

// Main analysis called from "UserExec", applied for all the analyzed events
void AliAnalysisTaskSEOmegacZero2XiPifromKFP::MakeAnaOmegacZero(AliAODEvent *AODEvent, KFParticle PV)
{

  // set the magnetic field
  KFParticle::SetField(fBzkG);

  //get number of cascades in the event
  const UInt_t nCasc = AODEvent->GetNumberOfCascades();

  Double_t covP[21], covN[21], covB[21];
  const Int_t NDaughters = 2;

  //get mass values from PDG
  const Float_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); //in GeV
  const Float_t massXi  = TDatabasePDG::Instance()->GetParticle(3312)->Mass();

  //separating positive and negative tracks in the events (plus checks)
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

//==================== loop on cascades in the event =========================

  for (UInt_t iCasc=0; iCasc<nCasc; iCasc++) {

    AliAODcascade *casc = AODEvent->GetCascade(iCasc);

   //check whether cascade is reconstructed offline or on the fly
   // cout << "Cascade is reconstructed on the fly: " << casc->GetOnFlyStatus() << endl; //GetOnFlyStatus returns kTRUE (=1) if the particle is recontructed on fly during the tracking (kFALSE (=0)->reconstructed offline). Here 0 --> cascades reconstructed offline
   //Nota: V0 sreconstructed on the fly (see function AliCascadeVertexer::V0sTracks2CascadeVertices)

    // cascade cut
    if ( !fAnaCuts->SingleCascCuts(casc, kFALSE) ) continue; //cuts implemented in the corresponding class, kFALSE --> (see AliRDHFCutsKFP.cxx - the only thing the boolean affects is the particle species assigned to the charged cascade's daughter, if the boolean is false it is assumed to be a pion, if it's true a kaon)

    AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0)); //ptrack stands for positive track, i.e. track of the daughter with positive charge (see how the trak bends in the magnetic field) - first daughter of lambda<-cascade
    AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1)); //ntrack stands for negative track, i.e. track of the daughter with negative charge (see how the trak bends in the magnetic field) - second daughter of lambda<-cascade
    AliAODTrack *btrack = (AliAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0)); //bachelor track, pi- <--- cascade-

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

    if ( btrack->Charge()<0 ) {

//^^^^^^^^^^^^^^^^^^^^^^^^  LAMBDA   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      //reconstruct lambda with KFP
      const KFParticle *vDaughters[2] = {&kfpProton, &kfpPionMinus}; //ricostruisco la lambda a partire da protone e pione negativo
      KFParticle kfpLambda;
      kfpLambda.Construct(vDaughters, NDaughters);
      Float_t massLambda_Rec, err_massLambda;
      kfpLambda.GetMass(massLambda_Rec, err_massLambda);

      //checks
      if ( TMath::Abs(kfpLambda.GetE())<=TMath::Abs(kfpLambda.GetPz()) ){
        fHistCheckKF->Fill(2);
        continue;
      }
      fHistCheckKF->Fill(1);

      if ( (kfpLambda.GetNDF()<=0 || kfpLambda.GetChi2()<=0) ) continue;
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLambda) ) continue;
      if ( err_massLambda<=0 ) continue;
      if ( (kfpLambda.GetChi2()/kfpLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue;

      //calculate l/Δl for Lambda (distance between PV of the event and the decay point of the lamda, normalized to its error)
      Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
      Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
      Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
      Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
      Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
      if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
      dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
      Double_t nErr_l_Lambda = l_Lambda/dl_Lambda;

      // l/Deltal cut of Lambda
      if ( nErr_l_Lambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Lambda
      if ( TMath::Abs(massLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      KFParticle kfpLambda_m = kfpLambda;
      kfpLambda_m.SetNonlinearMassConstraint(massLambda);

      //checks
      if ( TMath::Abs(kfpLambda_m.GetE()) <= TMath::Abs(kfpLambda_m.GetPz()) ){
        fHistCheckKF->Fill(4);
        continue;
      }
      fHistCheckKF->Fill(3);

      if(!AliVertexingHFUtils::CheckKFParticleCov(kfpLambda_m)) continue;

//^^^^^^^^^^^^^^^^^^^  PI- <-- CASCADE-  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      //construct il pion (pi<-xi) KFP
      KFParticle kfpPionFromXi;
      kfpPionFromXi = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, -211); // pion-


//^^^^^^^^^^^^^^^^^^^  CASCADE  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      //reconstruct cascade with KFP
      KFParticle kfpXiMinus;
      const KFParticle *vXiDs[2] = {&kfpPionFromXi, &kfpLambda_m};
      kfpXiMinus.Construct(vXiDs, NDaughters);

      // checks on cascade
      if ( TMath::Abs(kfpXiMinus.GetE())<=TMath::Abs(kfpXiMinus.GetPz()) ){
        fHistCheckKF->Fill(2);
        continue;
      }
      fHistCheckKF->Fill(1);

      // err_mass_cascade > 0
      Float_t massXiMinus_Rec, err_massXiMinus;
      kfpXiMinus.GetMass(massXiMinus_Rec, err_massXiMinus);
      if ( err_massXiMinus<=0 ) continue;

      // chi2>0 && NDF>0
      if ( kfpXiMinus.GetNDF()<=0 || kfpXiMinus.GetChi2()<=0 ) continue;

      // Prefilter
      if ( kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

      // check covariance matrix
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus) ) continue;

      // mass window cut of Omega-
      if ( (TMath::Abs(massXiMinus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi())) ) continue;

      // mass constraint
      KFParticle kfpXiMinus_m = kfpXiMinus;
      kfpXiMinus_m.SetNonlinearMassConstraint(massXi);

      //checks
      if ( TMath::Abs(kfpXiMinus_m.GetE()) <= TMath::Abs(kfpXiMinus_m.GetPz()) ){
        fHistCheckKF->Fill(4);
        continue;
      }
      fHistCheckKF->Fill(3);

      if(!AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus_m)) continue;


//=============== loop on pi+<-Omegac0 ==================================

      for (Int_t itrkBP=0; itrkBP<flag_trkP; itrkBP++) { // Loop for first bachelor pion+ (loop on all positive tracks in the event)

        if ( trackP[itrkBP]->GetID() == ptrack->GetID() ) continue; //Id of the considered track different from ID of positive daughter of lambda

        //standard quality checks on track
        if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP]) ) continue;

        KFParticle kfpBP = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackP[itrkBP], 211);

        //reconstruct Omegac0 with KFP
        KFParticle kfpOmegac0;
        const KFParticle *vOmegacZeroDs[2] = {&kfpBP, &kfpXiMinus_m};
        kfpOmegac0.Construct(vOmegacZeroDs, NDaughters);

        // chi2>0 && NDF>0
        if ( kfpOmegac0.GetNDF()<=0 || kfpOmegac0.GetChi2()<=0 ) continue;

        // Prefilter
        if ( kfpOmegac0.GetChi2()/kfpOmegac0.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
        if ( kfpOmegac0.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;

        //check
        if ( TMath::Abs(kfpOmegac0.GetE())<=TMath::Abs(kfpOmegac0.GetPz()) ){
          fHistCheckKF->Fill(2);
          continue;
        }
        fHistCheckKF->Fill(1);

        //check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpOmegac0) ) continue;

        //err_massOmegac0 > 0
        Float_t massOmegac0_Rec, err_massOmegac0;
        kfpOmegac0.GetMass(massOmegac0_Rec, err_massOmegac0);
        if ( err_massOmegac0<=0 ) continue;


       //fill the tree
       FillTreeOmegac0(1, kfpOmegac0, trackP[itrkBP], kfpBP, kfpXiMinus, kfpXiMinus_m, kfpPionFromXi, btrack, casc, kfpLambda, kfpLambda_m, ptrack, ntrack, PV);

       kfpOmegac0.Clear();
       kfpBP.Clear();

      } //end loop for first bachelor pion+


      //LIKE-SIGN PAIRS (ximinus&pi-  ---> for combinatorial background studies)

      //loop over all the negative tracks of the event
      for (Int_t itrkBP_LS=0; itrkBP_LS<flag_trkN; itrkBP_LS++) {

        //check the negative track I'm considering is different from both the negative track of the final state (pi- <-- lambda) and the bachelor track (pi- <-- xi-)
        if ( trackN[itrkBP_LS]->GetID()==ntrack->GetID() || trackN[itrkBP_LS]->GetID()==btrack->GetID() ) continue;

        //quality checks on track
        if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_LS]) ) continue;

        //create pi- from the track
        KFParticle kfpBP_LS = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackN[itrkBP_LS], -211);

        //reconstruct Omegac0 bkg (Like-sign pairs)
        KFParticle kfpOmegac0_LS;
        const KFParticle *vOmegac0Ds_LS[2] = {&kfpBP_LS, &kfpXiMinus_m};
        kfpOmegac0_LS.Construct(vOmegac0Ds_LS, NDaughters);

        //chi2>0 && NDF>0
        if ( kfpOmegac0_LS.GetNDF()<=0 || kfpOmegac0_LS.GetChi2()<=0 ) continue;

        //Prefilter
        if ( kfpOmegac0_LS.GetChi2()/kfpOmegac0_LS.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
        if ( kfpOmegac0_LS.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;
        if ( TMath::Abs(kfpOmegac0_LS.GetE())<=TMath::Abs(kfpOmegac0_LS.GetPz()) ) continue;

        // check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpOmegac0_LS) ) continue;

        // err_massXic0 > 0
        Float_t massOmegac0_LS_Rec, err_massOmegac0_LS;
        kfpOmegac0_LS.GetMass(massOmegac0_LS_Rec, err_massOmegac0_LS);
        if ( err_massOmegac0_LS<=0 ) continue;

        //fill the tree
        FillTreeOmegac0(0, kfpOmegac0_LS, trackN[itrkBP_LS], kfpBP_LS, kfpXiMinus, kfpXiMinus_m, kfpPionFromXi, btrack, casc, kfpLambda, kfpLambda_m, ptrack, ntrack, PV);

        //clear
        kfpOmegac0_LS.Clear();
        kfpBP_LS.Clear();


      }

      kfpXiMinus_m.Clear();
      kfpXiMinus.Clear();
      kfpPionFromXi.Clear();
      kfpLambda_m.Clear();
      kfpLambda.Clear();
    }


    //ANTIPARTICLES
    //repeat same process (but ptrack <---> ntrack), then depending on pi<--cascade charge you see if you are dealing either with omega_c^0 or anti-omega_c^0
    if ( btrack->Charge()>0 ) {

//^^^^^^^^^^^^^^^^^^^   ANTI-LAMBDA  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      const KFParticle *vAntiDaughters[2] = {&kfpPionPlus, &kfpAntiProton};
      KFParticle kfpAntiLambda;
      kfpAntiLambda.Construct(vAntiDaughters, NDaughters);
      Float_t massAntiLambda_Rec, err_massAntiLambda;
      kfpAntiLambda.GetMass(massAntiLambda_Rec, err_massAntiLambda);

      // check
      if ( TMath::Abs(kfpAntiLambda.GetE())<=TMath::Abs(kfpAntiLambda.GetPz()) ){
        fHistCheckKF->Fill(2);
        continue;
      }
      fHistCheckKF->Fill(1);

      // chi2>0 && NDF>0 for selecting Anti-Lambda
      if ( kfpAntiLambda.GetNDF()<=0 || kfpAntiLambda.GetChi2()<=0 ) continue;

      // check cov. of Anti-Lambda
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda) ) continue;

      // err_mass>0 of Anti-Lambda
      if ( err_massAntiLambda<=0 ) continue;

      // Chi2geo cut of Anti-Lambda
      if ( (kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue;

      //calculate l/Δl for Anti-Lambda
      Double_t dx_AntiLambda = PV.GetX()-kfpAntiLambda.GetX();
      Double_t dy_AntiLambda = PV.GetY()-kfpAntiLambda.GetY();
      Double_t dz_AntiLambda = PV.GetZ()-kfpAntiLambda.GetZ();
      Double_t l_AntiLambda = TMath::Sqrt(dx_AntiLambda*dx_AntiLambda + dy_AntiLambda*dy_AntiLambda + dz_AntiLambda*dz_AntiLambda);
      Double_t dl_AntiLambda = (PV.GetCovariance(0)+kfpAntiLambda.GetCovariance(0))*dx_AntiLambda*dx_AntiLambda + (PV.GetCovariance(2)+kfpAntiLambda.GetCovariance(2))*dy_AntiLambda*dy_AntiLambda + (PV.GetCovariance(5)+kfpAntiLambda.GetCovariance(5))*dz_AntiLambda*dz_AntiLambda + 2*( (PV.GetCovariance(1)+kfpAntiLambda.GetCovariance(1))*dx_AntiLambda*dy_AntiLambda + (PV.GetCovariance(3)+kfpAntiLambda.GetCovariance(3))*dx_AntiLambda*dz_AntiLambda + (PV.GetCovariance(4)+kfpAntiLambda.GetCovariance(4))*dy_AntiLambda*dz_AntiLambda );
      if ( fabs(l_AntiLambda)<1.e-8f ) l_AntiLambda = 1.e-8f;
      dl_AntiLambda = dl_AntiLambda<0. ? 1.e8f : sqrt(dl_AntiLambda)/l_AntiLambda;
      Double_t nErr_l_AntiLambda = l_AntiLambda/dl_AntiLambda;

      // l/Deltal cut of Anti-Lambda
      if ( nErr_l_AntiLambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Anti-Lambda
      if ( TMath::Abs(massAntiLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      KFParticle kfpAntiLambda_m = kfpAntiLambda;
      kfpAntiLambda_m.SetNonlinearMassConstraint(massLambda);

      if (TMath::Abs(kfpAntiLambda_m.GetE()) <= TMath::Abs(kfpAntiLambda_m.GetPz()) ){
        fHistCheckKF->Fill(4);
        continue;
      }
      fHistCheckKF->Fill(3);

      if(!AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda_m)) continue;

//^^^^^^^^^^^^^^^^^^^^^^^ PI+ <-- CASCADE+  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      KFParticle kfpPionFromXi;
      kfpPionFromXi = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, 211); // pion+

//^^^^^^^^^^^^^^^^^^^^^^ ANTI-CASCADE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      KFParticle kfpXiPlus;
      const KFParticle *vXiDs[2] = {&kfpPionFromXi, &kfpAntiLambda_m};
      kfpXiPlus.Construct(vXiDs, NDaughters);

      // check
      if ( TMath::Abs(kfpXiPlus.GetE())<=TMath::Abs(kfpXiPlus.GetPz()) ){
        fHistCheckKF->Fill(2);
        continue;
      }
      fHistCheckKF->Fill(1);

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
      if ( (TMath::Abs(massXiPlus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi())) ) continue;

      KFParticle kfpXiPlus_m = kfpXiPlus;
      kfpXiPlus_m.SetNonlinearMassConstraint(massXi);

     if ( TMath::Abs(kfpXiPlus_m.GetE()) <= TMath::Abs(kfpXiPlus_m.GetPz()) ){
        fHistCheckKF->Fill(4);
        continue;
      }
      fHistCheckKF->Fill(3);

      if( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus_m) ) continue;


      //=============== loop on pi-<-AntiOmegac0 ==================================

      for (Int_t itrkBP=0; itrkBP<flag_trkN; itrkBP++) { // Loop for first bachelor pion-

        if ( trackN[itrkBP]->GetID() == ntrack->GetID() ) continue;

        if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP]) ) continue; //see the function in file AliRDHFCutsKFP.cxx  (here I set FilterBit 4)

        KFParticle kfpBP = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackN[itrkBP], -211);

        // reconstruct Anti-Omegac0
        KFParticle kfpAntiOmegac0;
        const KFParticle *vOmegac0Ds[2] = {&kfpBP, &kfpXiPlus_m};
        kfpAntiOmegac0.Construct(vOmegac0Ds, NDaughters);

        // chi2>0 && NDF>0
        if ( kfpAntiOmegac0.GetNDF()<=0 || kfpAntiOmegac0.GetChi2()<=0 ) continue;

        // Prefilter
        if ( kfpAntiOmegac0.GetChi2()/kfpAntiOmegac0.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
        if ( kfpAntiOmegac0.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;

        // check rapidity of Anti-Omegac0
        if ( TMath::Abs(kfpAntiOmegac0.GetE())<=TMath::Abs(kfpAntiOmegac0.GetPz()) ){
          fHistCheckKF->Fill(2);
          continue;
        }
        fHistCheckKF->Fill(1);

        // check covariance matrix (quality check, checks whether the uncertainties contained in the covariance matrix are reasonable)
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiOmegac0) ) continue;

        // err_massAntiOmegac0 > 0
        Float_t massAntiOmegac0_Rec, err_massAntiOmegac0;
        kfpAntiOmegac0.GetMass(massAntiOmegac0_Rec, err_massAntiOmegac0);
        if ( err_massAntiOmegac0<=0 ) continue;

        //fill the tree
        FillTreeOmegac0(1, kfpAntiOmegac0, trackN[itrkBP], kfpBP, kfpXiPlus, kfpXiPlus_m, kfpPionFromXi, btrack, casc, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV);

        kfpAntiOmegac0.Clear();
        kfpBP.Clear();
      } //end loop for first bachelor pion-

      //LIKE-SIGN PAIRS (xiplus&pi+  --->  for combinatorial background studies)

      //loop over all the positive tracks of the event
      for (Int_t itrkBP_LS=0; itrkBP_LS<flag_trkP; itrkBP_LS++) {

        //check the positive track I'm considering is different from both the positive track of the final state (pi+ <-- lambda) and the bachelor track (pi+ <-- xi+)
        if ( trackP[itrkBP_LS]->GetID()==ptrack->GetID() || trackP[itrkBP_LS]->GetID()==btrack->GetID() ) continue;

        //quality checks on track
        if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_LS]) ) continue;

        //create pi+ from the track
        KFParticle kfpBP_LS = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackP[itrkBP_LS], 211);

        //reconstruct anti-Omegac0 bkg (Like-sign pairs)
        KFParticle kfpAntiOmegac0_LS;
        const KFParticle *vAntiOmegac0Ds_LS[2] = {&kfpBP_LS, &kfpXiPlus_m};
        kfpAntiOmegac0_LS.Construct(vAntiOmegac0Ds_LS, NDaughters);

        //chi2>0 && NDF>0
        if ( kfpAntiOmegac0_LS.GetNDF()<=0 || kfpAntiOmegac0_LS.GetChi2()<=0 ) continue;

        //Prefilter
        if ( kfpAntiOmegac0_LS.GetChi2()/kfpAntiOmegac0_LS.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
        if ( kfpAntiOmegac0_LS.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;
        if ( TMath::Abs(kfpAntiOmegac0_LS.GetE())<=TMath::Abs(kfpAntiOmegac0_LS.GetPz()) ) continue;

        // check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiOmegac0_LS) ) continue;

        // err_massOmegac0 > 0
        Float_t massAntiOmegac0_LS_Rec, err_massAntiOmegac0_LS;
        kfpAntiOmegac0_LS.GetMass(massAntiOmegac0_LS_Rec, err_massAntiOmegac0_LS);
        if ( err_massAntiOmegac0_LS<=0 ) continue;

        //fill the tree
        FillTreeOmegac0(0, kfpAntiOmegac0_LS, trackP[itrkBP_LS], kfpBP_LS, kfpXiPlus, kfpXiPlus_m, kfpPionFromXi, btrack, casc, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV);

        //clear
        kfpAntiOmegac0_LS.Clear();
        kfpBP_LS.Clear();
      }

      kfpXiPlus_m.Clear();
      kfpXiPlus.Clear();
      kfpPionFromXi.Clear();
      kfpAntiLambda_m.Clear();
      kfpAntiLambda.Clear();
    }

    kfpPionPlus.Clear();
    kfpAntiProton.Clear();
    kfpPionMinus.Clear();
    kfpProton.Clear();

  }//end loop cascades

  return;
}




//______________________________________________________________________________
void AliAnalysisTaskSEOmegacZero2XiPifromKFP::DefineTreeOmegac0()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(3)->GetContainer()->GetName();
  fTree_Omegac0 = new TTree(nameoutput, "Omegac0 variables tree");
  Int_t nVar = 46;
  fVar_Omegac0 = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "nSigmaTPC_PiFromOmegac0"; // TPC nsigma for pion coming from Omegac0
  fVarNames[1]  = "nSigmaTOF_PiFromOmegac0"; // TOF nsigma for pion coming from Omegac0
  fVarNames[2]  = "nSigmaTPC_PiFromXi"; // TPC nsigma for pion coming from Xi (cascade)
  fVarNames[3]  = "nSigmaTOF_PiFromXi"; // TOF nsigma for pion coming from Xi (cascade)
  fVarNames[4]  = "nSigmaTPC_PiFromLam"; // TPC nsigma for pion coming from Lambda
  fVarNames[5]  = "nSigmaTPC_PrFromLam"; // TPC nsigma for proton coming from Lambda
  fVarNames[6]  = "DCA_XiDau"; // DCA of Cascade's daughters (calculated from AOD cascade)
  fVarNames[7]  = "chi2geo_Lam"; // chi2_geometry of Lambda
  fVarNames[8]  = "ldl_Lam"; // l/dl of Lambda
  fVarNames[9]  = "chi2topo_LamToPV"; // chi2_topo of Lambda to PV (after using SetNonlinearMassConstraint on Lambda)
  fVarNames[10] = "chi2geo_Xi"; // chi2_geometry of Xi (with Lambda mass const.)
  fVarNames[11] = "ldl_Xi"; // l/dl of Xi (with Lambda mass const.)
  fVarNames[12] = "chi2topo_XiToPV"; // chi2_topo of Xi to PV (after using SetNonlinearMassConstraint on Omega-)
  fVarNames[13] = "DecayLxy_Lam"; // decay length of Lambda in x-y plane
  fVarNames[14] = "DecayLxy_Xi"; // decay length of Xi in x-y plane
  fVarNames[15] = "PA_LamToXi"; // pointing angle of Lambda (pointing back to Xi)
  fVarNames[16] = "PA_LamToPV"; // pointing angle of Lambda (pointing back to PV)
  fVarNames[17] = "PA_XiToOmegac0"; // pointing angle of Xi (pointing back to Omegac0)
  fVarNames[18] = "PA_XiToPV"; // pointing angle of Xi (pointing back to PV)
  fVarNames[19] = "mass_Lam"; // mass of Lambda (without mass const.)
  fVarNames[20] = "mass_Xi"; // mass of Xi (without mass const.)
  fVarNames[21] = "pt_PiFromOmegac0"; // pt of pion coming from Omegac0
  fVarNames[22] = "pt_Omegac0"; // pt of Omegac0
  fVarNames[23] = "rap_Omegac0"; // rapidity of Omegac0
  fVarNames[24] = "mass_Omegac0"; // mass of Omegac0
  fVarNames[25] = "CosThetaStar_PiFromOmegac0"; // CosThetaStar of pion coming from Omegac0
  fVarNames[26] = "chi2prim_PiFromOmegac0"; // chi2_topo of pion to PV
  fVarNames[27] = "DCAxy_PiFromOmegac0_KF"; // DCA of pion coming from Omegac0 in x-y plane
  fVarNames[28] = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi (after using SetNonlinearMassConstraint on Lambda)
  fVarNames[29] = "chi2topo_XiToOmegac0"; // chi2_topo of Xi to Omegac0 (after using SetNonlinearMassConstraint on Xi-)
  fVarNames[30] = "DecayLxy_Omegac0"; // decay length of Omegac0 in x-y plane
  fVarNames[31] = "chi2geo_Omegac0"; // chi2_geometry of Omegac0
  fVarNames[32] = "DCA_LamDau"; // DCA of Lambda's daughters (calculated from AOD cascade)
  fVarNames[33] = "DCA_XiDau_KF"; // DCA of Xi's daughters (calculated from KF after Lambda mass constraint)
  fVarNames[34] = "DCA_Omegac0Dau_KF"; // DCA of Omegac0's daughters (calculated from KF after Omega mass constraint)
  fVarNames[35] = "DCAxy_XiToPV_KF"; // DCA of Xi to PV in x-y plane (calculated from KF after Omega mass constraint)
  fVarNames[36] = "mass_Omega"; // mass of Omega
  fVarNames[37] = "chi2topo_Omegac0ToPV"; // chi2_topo of Omegac0 to PV
  fVarNames[38] = "PA_Omegac0ToPV"; // pointing angle of Omegac0 (pointing back to PV)
  fVarNames[39] = "ldl_Omegac0"; // ldl of Omegac0
  fVarNames[40] = "ct_Omegac0"; // lifetime of Omegac0
  fVarNames[41] = "chi2mass_Lam"; // chi2 Lambda after using SetNonlinearMassConstraint (before fitting in PV/mother decay point)
  fVarNames[42] = "chi2mass_Xi"; // chi2 Xi after using SetNonlinearMassConstraint (before fitting in PV/mother decay point)
  fVarNames[43] = "flag_UnlikeOrLike_Sign"; // flag of unlike sign or like sign pair (0=likesign, 1=unlikesign) , to study the combinatorial bck (likesign)
  fVarNames[44] = "EventCounter"; //int that counts the number of event
  fVarNames[45] = "pseudorap_Omegac0"; //pseudorapidity of Omegac0

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Omegac0->Branch(fVarNames[ivar].Data(), &fVar_Omegac0[ivar], Form("%s/F", fVarNames[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEOmegacZero2XiPifromKFP::FillTreeOmegac0(Int_t flagUSorLS, KFParticle kfpOmegac0, AliAODTrack *trackPiFromOmegac0, KFParticle kfpBP, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, KFParticle kfpPionFromXi, AliAODTrack *trackPionFromXi, AliAODcascade *casc, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, KFParticle PV)
{

//fill Omegac0 TTree

  for (Int_t i=0; i<46; i++) {
    fVar_Omegac0[i] = -9999.;
  }

  //PID with nSigma criterion
  Float_t nSigmaTOF_PiFromOmegac0 = fPID->NumberOfSigmasTOF(trackPiFromOmegac0,AliPID::kPion); //calculate the number of sigmas (between the measurement expected for the indicated particle and the performed measurement) in the TOF
  Float_t nSigmaTPC_PiFromOmegac0 = fPID->NumberOfSigmasTPC(trackPiFromOmegac0,AliPID::kPion); //calculate the number of sigmas (between the measurement expected for the indicated particle and the performed measurement) in the TPC
  Float_t nSigmaTPC_PrFromLam  = fPID->NumberOfSigmasTPC(trkProton,AliPID::kProton);
  Float_t nSigmaTPC_PiFromLam  = fPID->NumberOfSigmasTPC(trkPion,AliPID::kPion);
  Float_t nSigmaTPC_PiFromXi = -9999., nSigmaTOF_PiFromXi = -9999.;
  nSigmaTPC_PiFromXi = fPID->NumberOfSigmasTPC(trackPionFromXi,AliPID::kPion);
  nSigmaTOF_PiFromXi = fPID->NumberOfSigmasTOF(trackPionFromXi,AliPID::kPion);

  if ( fabs(nSigmaTPC_PiFromOmegac0)>=fAnaCuts->GetPidPiFromXic0()->GetSigma(0) || fabs(nSigmaTPC_PiFromXi)>=fAnaCuts->GetPidPiFromXi()->GetSigma(0) || fabs(nSigmaTPC_PrFromLam)>=fAnaCuts->GetPidPrFromV0()->GetSigma(0) || fabs(nSigmaTPC_PiFromLam)>=fAnaCuts->GetPidPiFromV0()->GetSigma(0) ) return; //here I use nSigma criterion to do PID (see CutObject for the values)

  //fit Omegac0 track in PV
  KFParticle kfpOmegac0_PV = kfpOmegac0;
  kfpOmegac0_PV.SetProductionVertex(PV);

  //fit lXi- track in Omegac0 decay vertex
  KFParticle kfpXiMinus_Omegac0 = kfpXiMinus_m;
  kfpXiMinus_Omegac0.SetProductionVertex(kfpOmegac0);

  //fit Xi- track in PV
  KFParticle kfpXiMinus_PV = kfpXiMinus_m;
  kfpXiMinus_PV.SetProductionVertex(PV);

  //fit lambda track in Xi- decay vertex
  KFParticle kfpLambda_Xi = kfpLambda_m;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus);

  //fit lambda trck in PV
  KFParticle kfpLambda_PV = kfpLambda_m;
  kfpLambda_PV.SetProductionVertex(PV);

  //fit pi<-Omegac0 track in Omegac0 decay vertex
  KFParticle kfpBP_Omegac0 = kfpBP;
  kfpBP_Omegac0.SetProductionVertex(kfpOmegac0);

  //calculate l/Δl for Lambda
  Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
  Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
  Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
  Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
  Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
  if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
  dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
  if ( dl_Lambda<=0 ) return;

  //calculate l/Δl for Xi-
  Double_t dx_Xi = PV.GetX()-kfpXiMinus.GetX();
  Double_t dy_Xi = PV.GetY()-kfpXiMinus.GetY();
  Double_t dz_Xi = PV.GetZ()-kfpXiMinus.GetZ();
  Double_t l_Xi = TMath::Sqrt(dx_Xi*dx_Xi + dy_Xi*dy_Xi + dz_Xi*dz_Xi);
  Double_t dl_Xi = (PV.GetCovariance(0)+kfpXiMinus.GetCovariance(0))*dx_Xi*dx_Xi + (PV.GetCovariance(2)+kfpXiMinus.GetCovariance(2))*dy_Xi*dy_Xi + (PV.GetCovariance(5)+kfpXiMinus.GetCovariance(5))*dz_Xi*dz_Xi + 2*( (PV.GetCovariance(1)+kfpXiMinus.GetCovariance(1))*dx_Xi*dy_Xi + (PV.GetCovariance(3)+kfpXiMinus.GetCovariance(3))*dx_Xi*dz_Xi + (PV.GetCovariance(4)+kfpXiMinus.GetCovariance(4))*dy_Xi*dz_Xi );
  if ( fabs(l_Xi)<1.e-8f ) l_Xi = 1.e-8f;
  dl_Xi = dl_Xi<0. ? 1.e8f : sqrt(dl_Xi)/l_Xi;
  if ( dl_Xi<=0 ) return;

  //calculate l/Δl for Omegac0
  Double_t dx_Omegac0 = PV.GetX()-kfpOmegac0.GetX();
  Double_t dy_Omegac0 = PV.GetY()-kfpOmegac0.GetY();
  Double_t dz_Omegac0 = PV.GetZ()-kfpOmegac0.GetZ();
  Double_t l_Omegac0 = TMath::Sqrt(dx_Omegac0*dx_Omegac0 + dy_Omegac0*dy_Omegac0 + dz_Omegac0*dz_Omegac0);
  Double_t dl_Omegac0 = (PV.GetCovariance(0)+kfpOmegac0.GetCovariance(0))*dx_Omegac0*dx_Omegac0 + (PV.GetCovariance(2)+kfpOmegac0.GetCovariance(2))*dy_Omegac0*dy_Omegac0 + (PV.GetCovariance(5)+kfpOmegac0.GetCovariance(5))*dz_Omegac0*dz_Omegac0 + 2*( (PV.GetCovariance(1)+kfpOmegac0.GetCovariance(1))*dx_Omegac0*dy_Omegac0 + (PV.GetCovariance(3)+kfpOmegac0.GetCovariance(3))*dx_Omegac0*dz_Omegac0 + (PV.GetCovariance(4)+kfpOmegac0.GetCovariance(4))*dy_Omegac0*dz_Omegac0 );
  if ( fabs(l_Omegac0)<1.e-8f ) l_Omegac0 = 1.e-8f;
  dl_Omegac0 = dl_Omegac0<0. ? 1.e8f : sqrt(dl_Omegac0)/l_Omegac0;
  if ( dl_Omegac0<=0 ) return;

  //apply cuts contained in CutObject
  if ( kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF() <= fAnaCuts->GetKFPLam_Chi2topoMin() ) return;
  if ( kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF() >= fAnaCuts->GetKFPXi_Chi2topoMax() ) return;
  if ( l_Xi/dl_Xi <= fAnaCuts->GetKFPXi_lDeltalMin() ) return;

  //mass of Omegac0 -->  checks and save
  const Float_t PDGmassOmegac0 = 2.6952;	//value from PDG
  Float_t mass_Omegac0_PV, err_mass_Omegac0_PV;
  kfpOmegac0_PV.GetMass(mass_Omegac0_PV, err_mass_Omegac0_PV);
  fVar_Omegac0[24] = mass_Omegac0_PV; // mass of Omegac0
  if ( (fabs(mass_Omegac0_PV-PDGmassOmegac0) > fAnaCuts->GetProdMassTolXic0()) ) return;


  fVar_Omegac0[0]  = nSigmaTPC_PiFromOmegac0;
  fVar_Omegac0[1]  = nSigmaTOF_PiFromOmegac0;
  fVar_Omegac0[2]  = nSigmaTPC_PiFromXi;
  fVar_Omegac0[3]  = nSigmaTOF_PiFromXi;
  fVar_Omegac0[4]  = nSigmaTPC_PiFromLam;
  fVar_Omegac0[5]  = nSigmaTPC_PrFromLam;
  fVar_Omegac0[6]  = casc->DcaXiDaughters(); // DCA_XiDau
  fVar_Omegac0[7]  = kfpLambda.GetChi2()/kfpLambda.GetNDF(); // chi2geo_Lam
  fVar_Omegac0[8]  = l_Lambda/dl_Lambda;
  fVar_Omegac0[9]  = kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF();
  fVar_Omegac0[10] = kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF();
  fVar_Omegac0[11] = l_Xi/dl_Xi;
  fVar_Omegac0[12] = kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF();

  //decay lenghts
  Float_t DecayLxy_Lam, err_DecayLxy_Lam;
  kfpLambda_Xi.GetDecayLengthXY(DecayLxy_Lam, err_DecayLxy_Lam);
  fVar_Omegac0[13] = DecayLxy_Lam;

  Float_t DecayLxy_Xi, err_DecayLxy_Xi;
  kfpXiMinus_PV.GetDecayLengthXY(DecayLxy_Xi, err_DecayLxy_Xi);
  fVar_Omegac0[14] = DecayLxy_Xi;

  // calculate CosPointingAngle
  Double_t cosPA_v0toXi = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLambda_m, kfpXiMinus);
  Double_t cosPA_v0toPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLambda_m, PV);
  Double_t cosPA_XiToOmegac0 = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXiMinus_m, kfpOmegac0);
  Double_t cosPA_XiToPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXiMinus_m, PV);
  fVar_Omegac0[15] = TMath::ACos(cosPA_v0toXi); // PA_LamToXi
  fVar_Omegac0[16] = TMath::ACos(cosPA_v0toPV); // PA_LamToPV
  fVar_Omegac0[17] = TMath::ACos(cosPA_XiToOmegac0); // PA_XiToOmegac0
  fVar_Omegac0[18] = TMath::ACos(cosPA_XiToPV); // PA_XiToPV

  //masses
  Float_t mass_Lam, err_mass_Lam;
  kfpLambda.GetMass(mass_Lam, err_mass_Lam);
  fVar_Omegac0[19] = mass_Lam;

  Float_t mass_Xi, err_mass_Xi;
  kfpXiMinus.GetMass(mass_Xi, err_mass_Xi);
  fVar_Omegac0[20] = mass_Xi;

  //pt
  fVar_Omegac0[21] = trackPiFromOmegac0->Pt();

  //rapidity
  fVar_Omegac0[22] = kfpOmegac0_PV.GetPt();
  if ( TMath::Abs(kfpOmegac0_PV.GetE())>TMath::Abs(kfpOmegac0_PV.GetPz()) ) {
    fVar_Omegac0[23] = kfpOmegac0_PV.GetRapidity();
    fVar_Omegac0[45] = kfpOmegac0_PV.GetEta();
  }

  //angle
  fVar_Omegac0[25] = AliVertexingHFUtils::CosThetaStarFromKF(0, 4332, 211, 3312, kfpOmegac0, kfpBP_Omegac0, kfpXiMinus_Omegac0);


  //fit pi<-Omegac0 track in PV and save chi2 (topo)
  KFParticle kfpBP_PV = kfpBP;
  kfpBP_PV.SetProductionVertex(PV);
  fVar_Omegac0[26] = kfpBP_PV.GetChi2()/kfpBP_PV.GetNDF();

  //DCA of Pion<-Omegac0 to PV
  fVar_Omegac0[27] = kfpBP.GetDistanceFromVertexXY(PV); // DCA of pion in x-y

  //chi2
  fVar_Omegac0[28] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();
  fVar_Omegac0[29] = kfpXiMinus_Omegac0.GetChi2()/kfpXiMinus_Omegac0.GetNDF();

  //decay lengths
  Float_t DecayLxy_Omegac0, err_DecayLxy_Omegac0;
  kfpOmegac0_PV.GetDecayLengthXY(DecayLxy_Omegac0, err_DecayLxy_Omegac0);
  fVar_Omegac0[30] = DecayLxy_Omegac0;

  //chi2
  fVar_Omegac0[31] = kfpOmegac0.GetChi2()/kfpOmegac0.GetNDF();

  //DCAs
  fVar_Omegac0[32] = casc->DcaV0Daughters(); // DCA_LamDau
  fVar_Omegac0[33] = kfpPionFromXi.GetDistanceFromParticle(kfpLambda_m); // DCA_XiDau_KF - DCA of cascade-'s daughters (nel mio caso Xi-'s daughters)
  fVar_Omegac0[34] = kfpBP.GetDistanceFromParticle(kfpXiMinus_m); // DCA_Omegac0Dau_KF
  fVar_Omegac0[35] = kfpXiMinus_m.GetDistanceFromVertexXY(PV); // DCA of Xi in x-y

  //check rej
  KFParticle kfpPionOrKaon_Rej;
  if(trackPionFromXi->Charge()<0){
    kfpPionOrKaon_Rej = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPionFromXi, -321); // -321 K-
  }
  if(trackPionFromXi->Charge()>0){
    kfpPionOrKaon_Rej = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPionFromXi, 321); // 321 K+
  }
  KFParticle kfpCasc_Rej;
  const KFParticle *vCasc_Rej_Ds[2] = {&kfpPionOrKaon_Rej, &kfpLambda_m};
  kfpCasc_Rej.Construct(vCasc_Rej_Ds, 2);
  Float_t massCasc_Rej, err_massCasc_Rej;
  kfpCasc_Rej.GetMass(massCasc_Rej, err_massCasc_Rej);
  fVar_Omegac0[36] = massCasc_Rej;

  //chi2topo
  fVar_Omegac0[37] = kfpOmegac0_PV.GetChi2()/kfpOmegac0_PV.GetNDF(); // chi2topo_Omegac0ToPV

  //pointing angle
  Double_t cosPA_Omegac0ToPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpOmegac0, PV);
  fVar_Omegac0[38] = TMath::ACos(cosPA_Omegac0ToPV); // PA_Omegac0ToPV

  //l/dl
  fVar_Omegac0[39] = l_Omegac0/dl_Omegac0; // ldl_Omegac0

  //lifetime
  Float_t ct_Omegac0=0., err_ct_Omegac0=0.;
  kfpOmegac0_PV.GetLifeTime(ct_Omegac0, err_ct_Omegac0);
  fVar_Omegac0[40] = ct_Omegac0; // lifetime of Omegac0

  //chi2_mass
  fVar_Omegac0[41]  = kfpLambda_m.GetChi2()/kfpLambda_m.GetNDF(); // chi2mass_lambda before fitting in PV/mother decay vertex
  fVar_Omegac0[42]  = kfpXiMinus_m.GetChi2()/kfpXiMinus_m.GetNDF(); // chi2mass_omega before fitting in PV/mother decay vertex

  //flag like/unlike sign
  fVar_Omegac0[43]=flagUSorLS;

  //flag like/unlike sign
  fVar_Omegac0[44]=fEvCount;

  //fill tree
  fTree_Omegac0->Fill();

  return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEOmegacZero2XiPifromKFP::MakeMCCheck(TClonesArray *mcArray)
{
    // Analysis AliAODMCParticle
    
    Double_t Num_Omegac0 = 0;
    Int_t    nmcpart = mcArray -> GetEntriesFast();
    Int_t    NDaughtersOmegac = 2;
    Int_t    NDaughtersXi = 2;
    Int_t    NDaughtersLambda = 2;
    Int_t    sign = 99;
    
    for (Int_t i=0; i<nmcpart; i++){
        AliAODMCParticle *mcpart = NULL;
        mcpart = (AliAODMCParticle*)mcArray->At(i);
                
        if(TMath::Abs(mcpart->GetPdgCode())==4332 && mcpart->GetNDaughters()==NDaughtersOmegac){  // 4332: Omegac0
          
            Int_t  CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcpart,kTRUE);
            Bool_t pi_flag = kFALSE;
            Bool_t xi_flag = kFALSE;
            Bool_t lambda_flag = kFALSE;
            Bool_t pifromxi_flag = kFALSE;
            Bool_t proton_flag = kFALSE;
            Bool_t pifromlambda_flag = kFALSE;
            Int_t    correctDecayChannel = 0;
            AliAODMCParticle *mcpipart = NULL;
            AliAODMCParticle *mcxipart = NULL;
            AliAODMCParticle *mclambdapart = NULL;
            AliAODMCParticle *mcpifromxipart = NULL;
            AliAODMCParticle *mcpifromlambdapart = NULL;
            AliAODMCParticle *mcprotonpart = NULL;

            for(Int_t idau = mcpart->GetDaughterFirst(); idau<mcpart->GetDaughterLast()+1; idau++){
                if(idau <0) break;
                AliAODMCParticle *mcdau = (AliAODMCParticle*)mcArray->At(idau);
                if(TMath::Abs(mcdau->GetPdgCode())==3312 && mcdau->GetNDaughters()==NDaughtersXi){  // 3312: xi
                    xi_flag = kTRUE;
                    mcxipart = mcdau;
                }
                if(TMath::Abs(mcdau->GetPdgCode())==211){   // 211: pion
                    pi_flag = kTRUE;
                    mcpipart = mcdau;
                }
            }

            if(pi_flag && xi_flag){
              for(Int_t idau = mcxipart->GetDaughterFirst(); idau<mcxipart->GetDaughterLast()+1; idau++){
                if(idau <0) break;
                AliAODMCParticle *mcdau = (AliAODMCParticle*)mcArray->At(idau);
                if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==NDaughtersLambda){  // 3122: lambda
                    lambda_flag = kTRUE;
                    mclambdapart = mcdau;
                }
                if(TMath::Abs(mcdau->GetPdgCode())== 211){   // 211: pion
                    pifromxi_flag = kTRUE;
                    mcpifromxipart = mcdau;
                }
              }
            }

            if(pi_flag && xi_flag && pifromxi_flag && lambda_flag){
              for(Int_t idau = mclambdapart->GetDaughterFirst(); idau<mclambdapart->GetDaughterLast()+1; idau++){
                if(idau <0) break;
                AliAODMCParticle *mcdau = (AliAODMCParticle*)mcArray->At(idau);
                if(TMath::Abs(mcdau->GetPdgCode())==2212){  // 2212: proton
                    proton_flag = kTRUE;
                    mcprotonpart = mcdau;
                }
                if(TMath::Abs(mcdau->GetPdgCode())== 211){   // 211: pion
                    pifromlambda_flag = kTRUE;
                    mcpifromlambdapart = mcdau;
                }
              }
            }

            if(pi_flag && xi_flag && pifromxi_flag && lambda_flag && proton_flag && pifromlambda_flag){
              correctDecayChannel = 2;
            } else {
              correctDecayChannel = -2;
            }

            if(pi_flag && xi_flag && pifromxi_flag && lambda_flag && proton_flag && pifromlambda_flag){
            if(mcxipart->Charge() < 0 && mcxipart->Charge() != -99){
              sign = 1; //omegac0
            } else if (mcxipart->Charge() > 0){
              sign = -1; //anti-omegac0
            } else {
              sign = 4;
            }
            }

            if(pi_flag && xi_flag && pifromxi_flag && lambda_flag && proton_flag && pifromlambda_flag){
              fHistMCDecayChain->Fill(3);
            } else {
              fHistMCDecayChain->Fill(-3);
            }

            AliAODMCParticle *mcdau_0 = (AliAODMCParticle*)mcArray->At(mcpart->GetDaughterFirst());
            Double_t MLOverP = sqrt( pow(mcpart->Xv() - mcdau_0->Xv(),2.) +  pow(mcpart->Yv() - mcdau_0->Yv(),2.) +  pow(mcpart->Zv() - mcdau_0->Zv(),2.)) * mcpart-> M() / mcpart->P()*1.e4; //cosi' espresso in micrometri

              
            fHistMCPdgCode->Fill(mcpart->GetPdgCode()); //pdgcode of omegac0 (netries=number of generated omegac0)
            fHistMCOmegacDauNumber->Fill(mcpart->GetNDaughters()); //number of daughters of omegac0 (netries=number of generated omegac0)
            fHistMCOmegacPt->Fill(mcpart->Pt()); //Pt distribution of generated omegac0 (netries=number of generated omegac0)
            fHistMCOmegacEta->Fill(mcpart->Eta()); //pseudorapidity distribution of simulated omegac0 (netries=number of generated omegac0)
            fHistMCOmegacSign->Fill(sign); //omegac0 or antiomegac0 (netries=number of generated omegac0)
            fHistMCOmegacCTau->Fill(MLOverP); //ctau distribution of simulated omegac0 (netries=number of generated omegac0)
            fHistMCOmegacM->Fill(mcpart->M()); //mass distribution of simulated omegac0 (netries=number of generated omegac0)
            fHistMCOrigin->Fill(CheckOrigin); //0 no quark found, 4 prompt, 5 non prompt (from beauty)

            } //close loop over omegac0      
    }  // close loop over nmcpart
    
    return kTRUE;
    
}//close MakeMCCheck
