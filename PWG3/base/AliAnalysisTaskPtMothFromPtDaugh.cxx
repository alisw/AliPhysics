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

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//             Class AnalysisTaskAliPtMothFromPtDaugh                   //
//   AnalysisTaskSE used for the reconstruction of mothers particles    //
//   spectra (pT and pTMin) starting from the pT-spectra of             //
//   daughters particles.                                               //
//                                                                      //
//             Authors: Giuseppe Bruno & Fiorella Fionda                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TChain.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliAnalysisTaskPtMothFromPtDaugh.h"
#include "AliPtMothFromPtDaugh.h"
#include "AliAnalysisTaskSE.h"

ClassImp(AliAnalysisTaskPtMothFromPtDaugh)

//_________________________________________________________________________________
AliAnalysisTaskPtMothFromPtDaugh::AliAnalysisTaskPtMothFromPtDaugh() :
  AliAnalysisTaskSE(),
  fPtMothDaugh(0),
  fDecayKine(0),
  fReadKineFromNtupla(0),
  fFileNtuplaName(0),
  fList(0)
{
  //
  // Default Constructor
  //
  AliInfo("Default Constructor!\n");  
}

//_________________________________________________________________________________
AliAnalysisTaskPtMothFromPtDaugh::AliAnalysisTaskPtMothFromPtDaugh(Bool_t IsNtuplaCreated) :
  AliAnalysisTaskSE("TaskAliPtMothFromDaugh"),
  fPtMothDaugh(0),
  fDecayKine(0),
  fReadKineFromNtupla(IsNtuplaCreated),
  fFileNtuplaName(0),
  fList(0)
{
    //
    // Basic AnalysisTaskSE Constructor  
    // Basic input of the analysis: TChain of galice.root files (optional)
    // Basic ouput: TList of mothers histograms (pT and pTMin) (standard)
    //              TNtupla with kinematics informations of mothers 
    //              and daughters (optional)
    // If Ntupla already exists loop on kinematics is not needed
    // therefore input is not defined
    //
    if(!IsNtuplaCreated){
           DefineInput(0,TChain::Class()); DefineOutput(2,TNtuple::Class());}
    else {AliInfo("Ntupla is already created! Loop on events is skipped!\n");}
    DefineOutput(1,TList::Class());
}

//___________________________________________________________________________________
AliAnalysisTaskPtMothFromPtDaugh::~AliAnalysisTaskPtMothFromPtDaugh()
   {
    //
    // Destructor
    //
    if(fPtMothDaugh)
     { delete fPtMothDaugh; fPtMothDaugh=0;}
    if(fDecayKine)
     { delete fDecayKine; fDecayKine=0;}
    if(fFileNtuplaName)
     { delete fFileNtuplaName; fFileNtuplaName=0;}
    if (fList)
     { delete fList; fList = 0; }
   }



//_________________________________________________________________________________
void AliAnalysisTaskPtMothFromPtDaugh::UserCreateOutputObjects()
{
  //
  // Initialise the framework objects
  //
   //
    // Initialise the framework objects. OutPut objects created are:
    // 1) TList with mothers histograms
    // 2) TNtuple with kinematics informations of mothers and daughters (optional)  
    //
    fList = new TList();
    fList->SetOwner();
    if(fReadKineFromNtupla) return;

    if(!AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) {
     Fatal("UserCreateOutputObjects", "This task needs a MC handler");
     return;
     }
    fDecayKine=new TNtuple("DecayKine","Decay kinematics","pdgM:pxM:pyM:pzM:yM:etaM:pdgD:pxD:pyD:pzD:yD:etaD");
    return;
}

//_________________________________________________________________________________
void AliAnalysisTaskPtMothFromPtDaugh::UserExec(Option_t */*option*/)
{
    //
    // Main loop. Called for every event
    // This loop fill a TNtuple with kinematics informations for 
    // mother and daughter particles. TNtupla contains:  
    //                   pdg of Mothers
    //                    px      "
    //                    py      "
    //                    pz      "
    //                  rapidity  "
    //                   eta      "
    //                   pdg of Daughter
    //                    px      "
    //                    py      "
    //                    pz      "
    //                  rapidity  "
    //                   eta      "         
    //
    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcHandler) {AliError("Could not get MC event handler!"); return; }
    mcHandler->SetReadTR(kFALSE);
    AliMCEvent* mcEvent = mcHandler->MCEvent();
    if(!mcEvent){ AliError("Could not get MC event!"); return; }
    AliStack* stack = mcEvent->Stack();
    if(!stack){ AliError("Could not get stack!"); return; }
    Int_t nPrims = stack->GetNprimary();
    float *inf = new float[12];
    for (Int_t iTrack = 0; iTrack < nPrims; ++iTrack) {
    TParticle *part = stack->Particle(iTrack);
    Int_t pdg=TMath::Abs(part->GetPdgCode());
    Double_t y,y2;
    //check if particle is in mothers list
    Int_t labelDaugh;
    if(fPtMothDaugh->IsMothers(pdg))
           {
           //check if mother particle has a selected daugh
           if(!fPtMothDaugh->IsSelectedDaugh(part,labelDaugh,stack)) continue;
           TParticle *pDaugh=stack->Particle(labelDaugh);
           fPtMothDaugh->Rapidity(part,y);
           fPtMothDaugh->Rapidity(pDaugh,y2);

           inf[0]=part->GetPdgCode();
           inf[1]=part->Px();
           inf[2]=part->Py();
           inf[3]=part->Pz();
           inf[4]=y;
           inf[5]=part->Eta();
           inf[6]=pDaugh->GetPdgCode();
           inf[7]=pDaugh->Px();
           inf[8]=pDaugh->Py();
           inf[9]=pDaugh->Pz();
           inf[10]=y2;
           inf[11]=pDaugh->Eta();
           fDecayKine->Fill(inf);
         } //close if statement for mothers particles
        } //end of tracks loop
    PostData(2,fDecayKine);
    PostData(1,fList);
    delete [] inf;
    return;
}

//___________________________________________________________________________________
void AliAnalysisTaskPtMothFromPtDaugh::Terminate(Option_t */*option*/)
   {
    //
    // Terminate method called at the end of the events loop.
    // Get the Ntupla with kineamtics informations after the
    // events loop or read it from the file fFileNtuplaName
    // if the Ntupla is already created. Then use it to
    // evaluate pT and pTMin spectra for mothers by means
    // of the method implemented in AliPtMothFromPtDaugh
    //
    if(fReadKineFromNtupla) fDecayKine = ReadNtuplaFromFile(fFileNtuplaName);
    else{ 
         fDecayKine = dynamic_cast<TNtuple*>(GetOutputData(2));
         fList = dynamic_cast<TList*>(GetOutputData(1));
        }
    if(!fDecayKine) { AliInfo("TNtupla not available!\n"); return; }
    if(!fList) { AliInfo("TList not availble!\n"); return; }
    fPtMothDaugh->SetDecayNtupla(fDecayKine);
    fPtMothDaugh->EvaluatePtMoth();
    TH1F *fHistoPt = (TH1F*)fPtMothDaugh->GetHistoPtMother()->Clone();
    TH1F *fHistoPtMin = (TH1F*)fPtMothDaugh->GetHistoPtMinMother()->Clone();
    fList->Add(fHistoPt);
    fList->Add(fHistoPtMin);
    PostData(1,fList);
    return;
   }

//___________________________________________________________________________________
TNtuple *AliAnalysisTaskPtMothFromPtDaugh::ReadNtuplaFromFile(char *inFileName)
   {
    // 
    // Get Ntupla from the file inFileName
    // after the it is created. 
    // Input: name of the file - Output: TNtupla
    //
    TFile *f = new TFile(inFileName,"READ");
    if(!f) {AliError(Form("File %s with TNtupla doesn't exist!",inFileName)); return 0x0;}
    TNtuple *DecayKine=(TNtuple*)f->Get("DecayKine");
    if(!DecayKine) { AliError("The TNtupla doesn't exist!\n"); return 0x0;}
    return DecayKine;
   }
