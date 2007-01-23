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

/* $Id$ */

//********************************************************
// class AliRunAnalysis                                  *
// Analysis manager                                      *
// Author: Piotr.Skowronski@cern.ch                      *
//********************************************************

#include <TDirectory.h>

#include "AliRunAnalysis.h"
#include "AliLog.h"
#include "AliAnalysis.h"
#include "AliEventCut.h"
#include "AliReader.h"


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
    AliError("Reader is not set");
    return 1;
  }
 TDirectory* cwd = gDirectory; 
 Int_t nanal = fAnalysies.GetEntries();
 AliDebug(1,Form("There are %d analyses",nanal));
 /******************************/ 
 /*  Init Event                */ 
 /******************************/ 
 AliDebug(1,"Intializing analyses...");
 for (Int_t an = 0; an < nanal; an++)
  {   
      AliAnalysis* analysis = (AliAnalysis*)fAnalysies.At(an);
      AliDebug(1,Form("Intializing analysis %d,address=%#x, name=%s", 
		      an, analysis, analysis->GetName()));
      analysis->Init();
      AliDebug(1,Form("Init done for analysis %d",an));
  }
 AliDebug(1,"Intializing analyses... Done.");
  
 while (fReader->Next() == kFALSE)
  {
     AliAOD* eventrec = fReader->GetEventRec();
     AliAOD* eventsim = fReader->GetEventSim();

     /******************************/ 
     /*  Event Cut                 */ 
     /******************************/ 
     if ( Rejected(eventrec,eventsim) )
      {
        AliDebug(1,"Event rejected by Event Cut");
        continue; //Did not pass the 
      }
      
     /******************************/ 
     /*  Process Event             */ 
     /******************************/ 
     AliDebug(1,Form("There is %d analyses",fAnalysies.GetEntries()));
     
     for (Int_t an = 0; an < fAnalysies.GetEntries(); an++)
      {
          AliAnalysis* analysis = (AliAnalysis*)fAnalysies.At(an);
          analysis->ProcessEvent(eventrec,eventsim);
      }
    
  }//end of loop over events

 /******************************/ 
 /*  Finish Event              */ 
 /******************************/ 
 AliDebug(1,Form("Finishing analyses...\n There are %d anlyses",fAnalysies.GetEntries()));
 if (cwd) cwd->cd();
 for (Int_t an = 0; an < fAnalysies.GetEntries(); an++)
  {
      AliAnalysis* analysis = (AliAnalysis*)fAnalysies.At(an);
      AliDebug(1,Form("Calling Finish for analysis %d address %#x name=%s", 
		      an, analysis,analysis->GetName()));
      analysis->Finish();
      AliDebug(1,Form("Called Finish for analysis %d",an));
  }
 AliDebug(1,"Finishing done");

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
