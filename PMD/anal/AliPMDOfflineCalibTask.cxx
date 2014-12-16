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

/************************************
 * AliPMD offline calibration task
 *
 * Satyajit Jena, IIT Bombay
 * sjena@cern.ch
 * Fri Feb 12 13:30:19 IST 2010
 *
 ************************************/

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "TObjArray.h"

#include "AliPMDOfflineCalibTask.h"

ClassImp(AliPMDOfflineCalibTask)

//________________________________________________________________________
AliPMDOfflineCalibTask::AliPMDOfflineCalibTask(const char *name) 
: AliAnalysisTaskSE(name),
  fListOfHistos(0),
  fPmdCalibAdcP(0),
  fPmdCalibAdcC(0),
  fPmdCalibEntP(0),
  fPmdCalibEntC(0),
  fNEvents(0),
  fSelectedTrigger(new TObjArray()),
  fRejected(kTRUE)
{
  // Constructor
  
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliPMDOfflineCalibTask::UserCreateOutputObjects()
{
  fListOfHistos = new TList();

  fNEvents = new TH1I("hNEvents","Events Statistic", 3, 0, 5);

  fPmdCalibAdcP = new TH1F("fPmdCalibAdcP", "PMD Calibration ADC fill for PRE", 234795, 0, 234795);
  fPmdCalibAdcP->GetYaxis()->SetTitle("ADC");
  fPmdCalibAdcP->GetXaxis()->SetTitle("Cells in SMN-ROW-COL");

  fPmdCalibAdcC = new TH1F("fPmdCalibAdcC", "PMD Calibration ADC fill for CPV", 234795, 0, 234795);
  fPmdCalibAdcC->GetYaxis()->SetTitle("ADC");
  fPmdCalibAdcC->GetXaxis()->SetTitle("Cells in SMN-ROW-COL");

  fPmdCalibEntP = new TH1F("fPmdCalibEntP", "PMD Calibration Entries fill for PRE", 234795, 0, 234795);
  fPmdCalibEntP->GetYaxis()->SetTitle("Frequescy");
  fPmdCalibEntP->GetXaxis()->SetTitle("Cells in SMN-ROW-COL");

  fPmdCalibEntC = new TH1F("fPmdCalibEntC", "PMD Calibration Entries fill for CPV", 234795, 0, 234795);
  fPmdCalibEntC->GetYaxis()->SetTitle("Frequescy ");
  fPmdCalibEntC->GetXaxis()->SetTitle("Cells in SMN-ROW-COL");
 
  fListOfHistos->Add(fPmdCalibAdcP);
  fListOfHistos->Add(fPmdCalibAdcC);
  fListOfHistos->Add(fPmdCalibEntP);
  fListOfHistos->Add(fPmdCalibEntC);
  fListOfHistos->Add(fNEvents);
}

//________________________________________________________________________
void AliPMDOfflineCalibTask::UserExec(Option_t *) 
{

  AliVEvent *event = InputEvent();
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }

  fNEvents->Fill(1);

  AliESDEvent* pmdesd = dynamic_cast<AliESDEvent*>(event);
  if(!pmdesd) return; 
  fNEvents->Fill(2);

  Bool_t pass = kTRUE;
  Int_t numberOfTriggerSelected = fSelectedTrigger->GetEntriesFast();
  
  if(fRejected) 
    {
      pass = kTRUE;
      for(Int_t k = 0; k < numberOfTriggerSelected; k++)
	{
	  const TObjString *const obString = (TObjString*)fSelectedTrigger->At(k);
	
	  const TString tString = obString->GetString();
	  if(pmdesd->IsTriggerClassFired((const char*)tString)) 
	    {
	      pass = kFALSE;
	    }
	}
    }
  else 
    {
      pass = kFALSE;
      for(Int_t k = 0; k < numberOfTriggerSelected; k++)
	{
	  const TObjString *const obString=(TObjString*)fSelectedTrigger->At(k);
	  
	  const TString tString = obString->GetString();
	  if(pmdesd->IsTriggerClassFired((const char*)tString)) 
	    {
	      pass = kTRUE;
	    }
	}
    }
  
  if(!pass) 
    {
      PostData(1, fListOfHistos);
      return;
    }

  fNEvents->Fill(3);
  
  Int_t npmdcl = pmdesd->GetNumberOfPmdTracks();
  if(npmdcl < 1) fNEvents->Fill(4);

  while (npmdcl--)
    {
      AliESDPmdTrack *pmdtr = pmdesd->GetPmdTrack(npmdcl);
      
      Int_t   det     = pmdtr->GetDetector();
      Int_t   smn     = pmdtr->GetSmn();
      Float_t adc     = pmdtr->GetClusterADC();
      Float_t sTag    = pmdtr->GetClusterSigmaX();
      Float_t sRowCol = pmdtr->GetClusterSigmaY();
      
      
      Float_t rc      = smn*10000 + sRowCol;
      if(adc > 1200.) continue;
      
      if(sTag > 999. && sTag < 1000.) 
	{
	  if(det == 0) 
	    {
	      fPmdCalibAdcP->Fill(rc,adc);
	      fPmdCalibEntP->Fill(rc);
	    }
	  else if(det == 1) 
	    {
	      fPmdCalibAdcC->Fill(rc,adc);
	      fPmdCalibEntC->Fill(rc);
	    }
	}
    }
  

  PostData(1, fListOfHistos);
}      

//________________________________________________________________________
void AliPMDOfflineCalibTask::Terminate(Option_t *) 
{
  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
    
    fPmdCalibAdcP = dynamic_cast<TH1F*>(fListOfHistos->At(0));
    if(!fPmdCalibAdcP) {
      printf("ERROR: No ADC File Generated for PMD-PRE ");
      return;
    }
    
    fPmdCalibAdcC = dynamic_cast<TH1F*>(fListOfHistos->At(1));  
    if(!fPmdCalibAdcC) {
      printf("ERROR: No ADC File Generated for PMD-CPV ");
      return;
    }
    
    fPmdCalibEntP = dynamic_cast<TH1F*>(fListOfHistos->At(2));
    if(!fPmdCalibEntP) {
      printf("ERROR: No Entry File Generated for PMD-PRE ");
      printf(" No fhXyPRE ");
      return;
    }

    fPmdCalibEntC = dynamic_cast<TH1F*>(fListOfHistos->At(3));
    if(!fPmdCalibEntC) {
      printf("ERROR: No Entry File Generated for PMD-CPV ");
      return;
    }
   
    fNEvents = dynamic_cast<TH1I*>(fListOfHistos->At(4));
    if(!fNEvents) {
      printf("ERROR: No Entry File Generated for Event Counter ");
      return;
    }

 
    Info("AliPMDOfflineCalibTask","PMD offline Calibration Successfully finished");


}
