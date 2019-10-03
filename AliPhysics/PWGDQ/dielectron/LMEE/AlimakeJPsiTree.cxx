/*************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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
//                                                                       //
//    Class to produce charm cocktail                  (analysis task)   //
//                                        (description in .h file)       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// c++ includes
#include <iostream>
// ROOT includes
#include "TParticle.h"
#include "TTree.h"
#include "TChain.h"
#include "TList.h"
// AliRoot includes
#include "AlimakeJPsiTree.h" // this task
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"


ClassImp(AlimakeJPsiTree)

//______________________________| Default Constructor
AlimakeJPsiTree::AlimakeJPsiTree():
AliAnalysisTaskSE(),
  fMcEvent(0x0),
  fMcHandler(0x0),
  ftree(0x0),
  fmother(0x0),
  fdaughter1(0x0),
  fdaughter2(0x0),
  fOutputList(0x0)
{

}

//______________________________| Specific Constructor
AlimakeJPsiTree::AlimakeJPsiTree(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fMcEvent(0x0),
  fMcHandler(0x0),
  ftree(0x0),
  fmother(0x0),
  fdaughter1(0x0),
  fdaughter2(0x0),
  fOutputList(0x0)
{
  Info("AlimakeJPsiTree","Calling Constructor");
  // Output slot #1 writes into a TList container (nevents histogram)
  DefineInput(0, TChain::Class());
  DefineOutput(1,TList::Class());
}

//______________________________| Destructor
AlimakeJPsiTree::~AlimakeJPsiTree()
{
  // Destructor
  Info("~AlimakeJPsiTree","Calling Destructor");
  if (fOutputList) delete fOutputList;
 
}

//______________________________| User Output
void AlimakeJPsiTree::UserCreateOutputObjects()
{
  Info("AlimakeJPsiTree","CreateOutputObjects of task %s", GetName());
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  ftree = new TTree("JPsi_ee","JPsi_and_daughters");
  ftree->Branch("mother",&fmother);
  ftree->Branch("daughter1",&fdaughter1);
  ftree->Branch("daughter2",&fdaughter2);
  fOutputList->Add(ftree);

  PostData(1, fOutputList);
}

//______________________________| Init
void AlimakeJPsiTree::Init()
{
  if(fDebug > 1) printf("AlimakeJPsiTree::Init() \n");
  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}


//______________________________| User Exec
void AlimakeJPsiTree::UserExec(Option_t *)
{


  //
  //Nb event
  //
   
  if(fDebug > 1) printf("AlimakeJPsiTree::UserExec \n");
  Init();

  if(fMcHandler){
    fMcEvent = fMcHandler->MCEvent();
  }else{
    if(fDebug > 1) printf("AlimakeJPsiTree::Handler() fMcHandler=NULL\n");
    return;
  }
  if(!fMcEvent){
    if(fDebug > 1) printf("AlimakeJPsiTree::UserExec()   fMcEvent=NULL \n");
    return;
  }


  int nparticles = fMcEvent->GetNumberOfPrimaries();
  if(fDebug > 1) printf("Number of particles %d and primary %d\n",fMcEvent->GetNumberOfTracks(),fMcEvent->GetNumberOfPrimaries());


  //// loop over particles
  //============================ single loop =================================
  
  
  for(int iparticle=0; iparticle<nparticles;iparticle++){
    
    AliAODMCParticle *p = (AliAODMCParticle *) fMcEvent->GetTrack(iparticle);
    if (!p) continue;
    int pdg = TMath::Abs( p->PdgCode() );
    //printf("pdg of particles %d\n",pdg);

    // find JPsi
    if(pdg==443) {
      fmother = p;
      // loop over daughters
      int k1 = p->GetDaughterFirst();
      int k2 = p->GetDaughterLast();
      Bool_t found_positron = kFALSE;
      Bool_t found_electron = kFALSE;
      for(int d=k1; d <= k2; d++) {
	if(d>0){
	  AliAODMCParticle *decay = (AliAODMCParticle *) fMcEvent->GetTrack(d);
	  int pdg_daughter = decay->PdgCode();
	  if(fDebug > 1) printf("pdg of daughter %d\n",pdg_daughter);
	  if(pdg_daughter==11){
	    // found positron
	    fdaughter1 = decay;
	    found_positron = kTRUE;
	  }
	  if(pdg_daughter==-11){
	    // found electron
	    fdaughter2 = decay;
	    found_electron = kTRUE;
	  }
	}
      }// loop over daughters
      if(found_positron && found_electron){
	if(fDebug > 1) printf("Fill tree\n");
	ftree->Fill();
      } 
    } // find JPsi    
  }

  
  PostData(1, fOutputList);

  return;
}


//______________________________| Terminate
void AlimakeJPsiTree::Terminate(Option_t*)
{
  Info("Terminate","Start and end of Method");
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    
    ftree = (TTree *) fOutputList->FindObject("JPsi_ee");
    if (!ftree) {
      Printf("ERROR: tree not available");
      return;
    }
  }
  return;
}
