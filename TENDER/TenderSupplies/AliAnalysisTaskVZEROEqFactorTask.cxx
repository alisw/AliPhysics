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
//  AliAnalysisTaskVZEROEqFactorTask.cxx, February 12th 2014
//   --- David Dobrigkeit Chinellato
//
// This task is meant to set correct VZERO equalization factors in AliESDRun 
// so that AliCentrality makes use of the correct values. NB This task has to
// be executed prior to AliCentrality for this to work properly! It is meant
// to be used as a Tender. 
//
//    Comments, Suggestions, Bug reports: Please send them to: 
//          --- daviddc@ifi.unicamp.br
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#include "AliLog.h"

#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDRun.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"

#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskVZEROEqFactorTask.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"

#include "AliVZEROCalibData.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskVZEROEqFactorTask)

AliAnalysisTaskVZEROEqFactorTask::AliAnalysisTaskVZEROEqFactorTask()
: AliAnalysisTaskSE(), fListHist(0), fEqFactors(0), fCalibData(0), fRunNumber(0), fHistEventCounter(0), fisAOD(kFALSE)
//------------------------------------------------
// Tree Variables 
{
  // Dummy Constructor  
}

AliAnalysisTaskVZEROEqFactorTask::AliAnalysisTaskVZEROEqFactorTask(const char *name) 
  : AliAnalysisTaskSE(name), fListHist(0), fEqFactors(0), fCalibData(0), fRunNumber(0), fHistEventCounter(0), fisAOD(kFALSE)
{
  // Constructor
  DefineOutput(1, TList::Class());
}


AliAnalysisTaskVZEROEqFactorTask::~AliAnalysisTaskVZEROEqFactorTask()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fListHist){
      delete fListHist;
      fListHist = 0x0;
   }
}

//________________________________________________________________________
void AliAnalysisTaskVZEROEqFactorTask::UserCreateOutputObjects()
{
//------------------------------------------------
// Output: Empty at the moment 
//------------------------------------------------

   fListHist = new TList();
   fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

   if(! fHistEventCounter ) {
    //Histogram Output: Event-by-Event
    // --- Single "Events Processed" Counter at this stage 
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",4,0,4); 
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "Has ESD");
    fHistEventCounter->GetXaxis()->SetBinLabel(3, "Has ESDRun");
    fHistEventCounter->GetXaxis()->SetBinLabel(4, "Rewrote EqFactors");
    fListHist->Add(fHistEventCounter); 
   }

   //List of Histograms
   PostData(1, fListHist);
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskVZEROEqFactorTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliESDEvent *lESDevent = 0x0;
   AliAODEvent *lAODevent = 0x0;
  // Connect to the InputEvent   
  // After these lines, we should have an ESD/AOD event + the number of V0s in it.

   // Appropriate for ESD analysis! 
   //Count Processed Events 
   fHistEventCounter->Fill(0.5);

   if(fisAOD) {
     lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
     if (!lAODevent) {
       AliError("AOD event not available \n");
       return;
     }
   } else {
     lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
     if (!lESDevent) {
       AliError("ESD event not available \n");
       return;
     }
   }
   fHistEventCounter->Fill(1.5);
       
   
   //Acquire ESDRun object - Will be needed to invoke AliESDEvent::SetVZEROEqFactors
   Int_t runNumber=-1;
   // const AliESDRun *lESDRun;
   if(fisAOD) {
     runNumber = lAODevent->GetRunNumber();
       } else {
   //   lESDRun = lESDevent->GetESDRun();
   //   if (!lESDRun) {
   //     AliError("ERROR: lESDRun not available, won't be able to write Equalization Factors! Exiting. \n");
   //     return;
   //   }
     runNumber = lESDevent->GetRunNumber();
   }
   fHistEventCounter->Fill(2.5);

   //CDB Processing only necessary if Run Number changed! Check for change (no need to redo this every event) 
   if( runNumber != fRunNumber ){
      AliWarning("Run Changed! Reloading CDB values!");
      //Load CDB Entries - Mirroring AliVZEROReconstructor
      AliCDBManager *cdbmanager = AliCDBManager::Instance();
      cdbmanager->SetDefaultStorage("raw://");
      cdbmanager->SetRun(runNumber);
      if(!cdbmanager) AliFatal("No CDB Manager !");
      AliCDBEntry *entry7 = cdbmanager->Get("VZERO/Calib/EqualizationFactors");
      if (!entry7) AliFatal("VZERO equalization factors are not found in OCDB !");
      fEqFactors = (TH1F*)entry7->GetObject();

      //Load Calibration object fCalibData
      fCalibData = GetCalibData(); // Mirror AliVZEROReconstructor Functionality 
      fRunNumber = runNumber; //New Run
   }
   if(!fCalibData) AliFatal("No VZERO CalibData Object found!");


   Float_t factors[64];
   Float_t factorSum = 0;
   for(Int_t i = 0; i < 64; ++i) {
      factors[i] = fEqFactors->GetBinContent(i+1)*fCalibData->GetMIPperADC(i);
      factorSum += factors[i];
   }
   for(Int_t i = 0; i < 64; ++i) { 
      factors[i] *= (64./factorSum);
   }

   // Set the equalized factors
   if(fisAOD) {
     lAODevent->SetVZEROEqFactors(factors);
   } else {
     lESDevent->SetVZEROEqFactors(factors);
   }
   fHistEventCounter->Fill(3.5);

   // Post output data.
   PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskVZEROEqFactorTask::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskVZEROEqFactorTask : ouput data container list not available\n");
      return;
   }	
	
   fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
   if (!fHistEventCounter) {
      Printf("ERROR - AliAnalysisTaskVZEROEqFactorTask : fHistEventCounter not available");
      return;
   }
  
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskVZEROEqFactorTask","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();

   fHistEventCounter->SetMarkerStyle(22);
   fHistEventCounter->DrawCopy("E");
}

//_____________________________________________________________________________
AliVZEROCalibData* AliAnalysisTaskVZEROEqFactorTask::GetCalibData() const
{
  // Gets calibration object for VZERO set
 
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry=0;
  entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibdata = 0;
  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");
  return calibdata;
}
