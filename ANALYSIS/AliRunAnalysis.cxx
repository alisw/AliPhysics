#include "AliRunAnalysis.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliRunAnalysis
//
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TString.h>
#include <TObjString.h>
#include <TClass.h>
#include <TFile.h>
#include <TKey.h>
#include <TObjArray.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliESDtrack.h>
#include <AliESD.h>


#include "AliEventCut.h"
#include "AliReader.h"
#include "AliVAODParticle.h"


ClassImp(AliRunAnalysis)
AliRunAnalysis::AliRunAnalysis():
 TTask("RunAnalysis","Alice Analysis Manager"),
 fAnalysies(10),
 fReader(0x0),
 fEventCut(0x0),
 fCutOnSim(kFALSE),
 fCutOnRec(kTRUE)
{
  //ctor
}
/*********************************************************/

AliRunAnalysis::~AliRunAnalysis()
{
  //dtor
  delete fReader;
  delete fEventCut;
}
/*********************************************************/

Int_t AliRunAnalysis::Run()
{
 //makes analysis

 if (fReader == 0x0)
  {
    Error("Run","Reader is not set");
    return 1;
  }
 TDirectory* cwd = gDirectory; 
 Int_t nanal = fAnalysies.GetEntries();
 if (AliVAODParticle::GetDebug()) Info("Run","There is %d analysies",nanal);
 /******************************/ 
 /*  Init Event                */ 
 /******************************/ 
 if (AliVAODParticle::GetDebug()) Info("Run","Intializing analyses...");
 for (Int_t an = 0; an < nanal; an++)
  {   
      if (AliVAODParticle::GetDebug()) Info("Run","Intializing analysis %d", an);
      AliAnalysis* analysis = (AliAnalysis*)fAnalysies.At(an);
      if (AliVAODParticle::GetDebug()) 
       { 
         Info("Run","Intializing analysis %d address %#x", an, analysis);
         Info("Run","Intializing analysis %d name %d", an, analysis->GetName());
         Info("Run","Intializing analysis %d: Calling Init...", an);
       } 
      analysis->Init();
      if (AliVAODParticle::GetDebug()) Info("Run","Intializing analysis %d: Calling Init... Done");
  }
 if (AliVAODParticle::GetDebug()) Info("Run","Intializing analyses... Done.");
  
 while (fReader->Next() == kFALSE)
  {
     AliAOD* eventrec = fReader->GetEventRec();
     AliAOD* eventsim = fReader->GetEventSim();

     /******************************/ 
     /*  Event Cut                 */ 
     /******************************/ 
     if ( Rejected(eventrec,eventsim) )
      {
        if (AliVAODParticle::GetDebug()) Info("Run","Event rejected by Event Cut");
        continue; //Did not pass the 
      }
      
     /******************************/ 
     /*  Process Event             */ 
     /******************************/ 
     if (AliVAODParticle::GetDebug())  Info("Run","There is %d analyses",fAnalysies.GetEntries());
     
     for (Int_t an = 0; an < fAnalysies.GetEntries(); an++)
      {
          AliAnalysis* analysis = (AliAnalysis*)fAnalysies.At(an);
          analysis->ProcessEvent(eventrec,eventsim);
      }
    
  }//end of loop over events

 /******************************/ 
 /*  Finish Event              */ 
 /******************************/ 
 if (AliVAODParticle::GetDebug()) Info("Run","Finishing analyses... ");
 if (AliVAODParticle::GetDebug()) Info("Run","There is %d anlyses",fAnalysies.GetEntries());
 if (cwd) cwd->cd();
 for (Int_t an = 0; an < fAnalysies.GetEntries(); an++)
  {
      if (AliVAODParticle::GetDebug()) Info("Run","Finishing analysis %d", an);
      AliAnalysis* analysis = (AliAnalysis*)fAnalysies.At(an);
      if (AliVAODParticle::GetDebug()) 
       { 
         Info("Run","Finishing analysis %d address %#x", an, analysis);
         Info("Run","Finishing analysis %d name %d", an, analysis->GetName());
         Info("Run","Finishing analysis %d: Calling Finish...",an);
       } 
      analysis->Finish();
      if (AliVAODParticle::GetDebug()) Info("Run","Finishing analysis %d: Calling Finish... Done");
  }
 if (AliVAODParticle::GetDebug()) Info("Run","Finishing done");

 return 0;   
}
/*********************************************************/

void  AliRunAnalysis::Add(AliAnalysis* a)
{
  //adds a to the list of analysis
  fAnalysies.Add(a);
}
/*********************************************************/

void AliRunAnalysis::SetEventCut(AliEventCut* evcut)
{
//Sets event -  makes a private copy
  delete fEventCut;
  if (evcut) fEventCut = (AliEventCut*)evcut->Clone();
  else fEventCut = 0x0;
}

/*********************************************************/

Bool_t AliRunAnalysis::Rejected(AliAOD* recevent, AliAOD* simevent)
{
  //checks the event cut
  if (fEventCut == 0x0) return kFALSE;
  
  if (fCutOnRec)
    if (fEventCut->Rejected(recevent)) return kTRUE;
    
  if (fCutOnSim)
    if (fEventCut->Rejected(simevent)) return kTRUE;
  
  return kFALSE;
}
