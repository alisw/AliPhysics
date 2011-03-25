// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This includes a TList of AliLRCProcess objects that are processing LRC analysis
// for a given Eta window 
// This task is worcking with ESD data only

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.5
// Version 3.5.5

#include "TChain.h"
#include "TTree.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskLRC.h"
#include <AliPID.h>

ClassImp(AliAnalysisTaskLRC)

//________________________________________________________________________
AliAnalysisTaskLRC::AliAnalysisTaskLRC(const char *name,Bool_t runKine) 
  : AliAnalysisTaskSE(name),fMaxPtLimit(5.0),fMinPtLimit(0.0),fDropKineE(kFALSE), fLRCproc(0),fOutList(0),fRunKine(0)
{
  //Init
  
  fRunKine=runKine;
  
  // Constructor

  // Output slot #1 writes into a TList container for all histogramms
  DefineOutput(1, TList::Class());
  
}


// ---------------------------------------  Setters ------------------

 
//________________________________________________________________________
void AliAnalysisTaskLRC::UserCreateOutputObjects()
{
  // Create histograms
  
  //LRC processors init
  Int_t lLrcNum=fLRCproc.GetEntries();
  
  for(Int_t i=0; i < lLrcNum; i++)
  {
    AliLRCProcess *p = (dynamic_cast<AliLRCProcess*> (fLRCproc.At(i)));
    if(p) p->InitDataMembers();
    else continue;
  }
  

    // --------- Output list
  
  fOutList = new TList();
 
   
  
  // ---------- Adding data members to output list
  
    
   for(Int_t i=0; i < lLrcNum; i++)
  {
  fOutList->Add((dynamic_cast<AliLRCProcess*> (fLRCproc.At(i)))->CreateOutput());
  }
  
  
    
}

//________________________________________________________________________
void AliAnalysisTaskLRC::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliVEvent *event = InputEvent();
   
   if( fRunKine ) 
   {
   	AliMCEvent *eventMC =    MCEvent(); 
  	event = eventMC;
   }

   
   
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }

   Int_t lLrcNum=fLRCproc.GetEntries();

    
  //Track variables
  double lPt;   // Temp Pt
  double lEta;	  // Temp ETA	
  
  
  for(Int_t i=0; i < lLrcNum; i++)
  {
    AliLRCProcess *prc = dynamic_cast<AliLRCProcess*> (fLRCproc.At(i));
    if(!prc) continue;
    else prc->StartEvent();
    //(dynamic_cast<AliLRCProcess*> (fLRCproc.At(i)))->StartEvent();
  }
 
    
    
    // Track loop 
  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
    AliVParticle* track = event->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    
    lPt=track->Pt();
    lEta=track->Eta();
    
    if(lPt<fMinPtLimit)continue;     // Dropping traks with lo Pt
    if(lPt>fMaxPtLimit)continue;     // Dropping traks with hi Pt
    
    if(fRunKine&&(track->Charge()*track->Charge()<0.5))continue; // Dropping neutral in Kine
    //if(fDropKineE&&fRunKine&&(track->PID()==AliPID::kElectron))continue; // Drop all e+- in Kine;
      
      for(Int_t i=0; i < lLrcNum; i++)
  {
    AliLRCProcess *p = (dynamic_cast<AliLRCProcess*> (fLRCproc.At(i)));
    if(p) p->AddTrackPtEta(lPt,lEta);
    else continue;
  }

    
        
  } //end of track loop 
  
  
  
   for(Int_t i=0; i < lLrcNum; i++)
  {
  (dynamic_cast<AliLRCProcess*> (fLRCproc.At(i)))->FinishEvent();
  }

  
   // Post output data.
  
  PostData(1, fOutList);
}      

//________________________________________________________________________
void AliAnalysisTaskLRC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
fOutList = dynamic_cast<TList*> (GetOutputData(0));
 
}


void AliAnalysisTaskLRC::AddLRCProcess(AliLRCProcess *newProc)
{
if(!newProc)
{Printf("ERROR:No AliLRCProcess object -  NULL pointer!");
return;
}
fLRCproc.Add(newProc);
return ;
}

