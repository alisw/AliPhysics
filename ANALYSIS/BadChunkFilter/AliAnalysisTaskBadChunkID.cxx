/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// MC bad chunk identifier 
// --- david.dobrigkeit.chinellato@cern.ch
//
// Loops over all chunks and fills a TTree object with a TString locating 
// chunk name for each event and a "number of global tracks" variable. 
//
// TTree is filled event-by-event but has only very few data members, 
// so memory consumption should still be reasonable. 
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
//#include "AliMCEventHandler.h"
//#include "AliMCEvent.h"
//#include "AliStack.h"
#include "AliAnalysisTaskBadChunkID.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskBadChunkID)

AliAnalysisTaskBadChunkID::AliAnalysisTaskBadChunkID() 
: AliAnalysisTaskSE(), fList(0), fTree(0),
  fHistNEvents(0),
  fRunNumber(0),
  fFileName(0),
  fNGlobalTracks(0),
  fNTracks(0)
{
  // Dummy Constructor

}

AliAnalysisTaskBadChunkID::AliAnalysisTaskBadChunkID(const char *name) 
  : AliAnalysisTaskSE(name), fList(0), fTree(0),
    fHistNEvents(0),
    fRunNumber(0),
    fFileName(0),
    fNGlobalTracks(0),
    fNTracks(0)
{
  // Constructor
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


AliAnalysisTaskBadChunkID::~AliAnalysisTaskBadChunkID()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fList){
      delete fList;
      fList = 0x0;
   }
   /*
   if (fTree){
      delete fTree;
      fTree = 0x0;
   }
   */
}

//________________________________________________________________________
void AliAnalysisTaskBadChunkID::UserCreateOutputObjects()
{

   
   // Create histograms
   OpenFile(1);
   fList = new TList();
   fList->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

   if(! fHistNEvents) {
      fHistNEvents = new TH1F("fHistNEvents", 
         "NumberOfEvents", 
         1, 0, 1); 		
      fList->Add(fHistNEvents);
   }

   OpenFile(2);	
   // Called once

//------------------------------------------------

   fTree = new TTree("fTree","V0Candidates");

//------------------------------------------------
// fTree Branch definitions - V0 Tree
//------------------------------------------------

//-----------BASIC-INFO---------------------------
   /*1*/ fTree->Branch("fNGlobalTracks",&fNGlobalTracks,"fNGlobalTracks/I");	
   /*2*/ fTree->Branch("fNTracks",&fNTracks,"fNTracks/I");	
   /*3*/ fTree->Branch("fRunNumber",&fRunNumber,"fRunNumber/I");	
   /*4*/ fTree->Branch("fFileName",&fFileName,16000,0);	

   //List of Histograms: Normal
   PostData(1, fList);

   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTree);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskBadChunkID::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliESDEvent *lESDevent = 0x0;
   //AliMCEvent  *lMCevent  = 0x0;
   //AliStack    *lMCstack  = 0x0;

  // Connect to the InputEvent   
  // After these lines, we should have an ESD/AOD event

   // Appropriate for ESD analysis! 
      
   lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
   if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
   }

/* Anyhow not needed for cross-check
 
   lMCevent = MCEvent();
   if (!lMCevent) {
      Printf("ERROR: Could not retrieve MC event \n");
      cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;   
      return;
   }

   lMCstack = lMCevent->Stack();
   if (!lMCstack) {
      Printf("ERROR: Could not retrieve MC stack \n");
      cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
      return;
   }
 
*/
//------------------------------------------------
// Track Multiplicity Information Acquistion
//------------------------------------------------

   Int_t lNTracks = -1; 
   Int_t lNGlobalTracks = 0;

   lNTracks = lESDevent->GetNumberOfTracks();

//----- Loop on Tracks --------------------------------------------------------------
   for (Int_t iCurrentTrack = 0; iCurrentTrack < lNTracks; iCurrentTrack++) 
   {// This is the begining of the loop on tracks
      AliESDtrack *lThisTrack=((AliESDEvent*)lESDevent)->GetTrack(iCurrentTrack);

      // kITSrefit refit condition: Global Track
      if( !(lThisTrack->GetStatus() & AliESDtrack::kITSrefit)) continue; 
      lNGlobalTracks++;
   }
//----- End Loop on Tracks ----------------------------------------------------------

  fHistNEvents->Fill(0.5); // valid event
  fNGlobalTracks  = lNGlobalTracks;
  fNTracks        = lNTracks;
  fRunNumber      = lESDevent->GetRunNumber();

  //Gymnastics to get the chunk number...
  TString lFileName = CurrentFileName();
  //  TObjArray *lLocationArray = lFileName.Tokenize("/");
  fFileName = fInputHandler->GetTree()->GetCurrentFile()->GetName();
  //  delete lLocationArray;
  fTree->Fill();

  //Printf("%i in run %i, Nglob = %i, Ntrack = %i \n",fChunkNumber,fRunNumber,fNGlobalTracks,fNTracks);

   // Post output data.
  PostData(1, fList);
  PostData(2, fTree);
}

//________________________________________________________________________
void AliAnalysisTaskBadChunkID::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   // Not interesting at this point.
}
