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

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliESDPmdTrack.h"

#include "AliPmdCustomCalibTask.h"

/*------------------------------------------------------------------------
  .  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  . 
  .              PMD Offiline Calibration Task to run Over
  .		 ESDs - can be used in user task 
  .		      - runs for single module in both plan
  .		      - Centrality is not needed for PbPb
  .                   - (Let me know if you need any updates)
  .
  .		     Satyajit Jena, IIT Bombay
  .		     sjena@cern.ch
  .		     3/8/2011
  .   
  .		     
  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  ------------------------------------------------------------------------*/

ClassImp(AliPmdCustomCalibTask)

//________________________________________________________________________
AliPmdCustomCalibTask::AliPmdCustomCalibTask() : AliAnalysisTaskSE(),
  fOutput(0),
  fSmn(10), 
  fHistAdcPre(0), 
  fHistAdcCpv(0),
  fHistClusterXYPre(0),           
  fHistClusterXYCpv(0) {
  for(Int_t i =0;i<48;i++) {
    for(Int_t j=0;j<96;j++)  {
      fHistAdcPreRC[i][j] = NULL;
      fHistAdcCpvRC[i][j] = NULL;
    }
  }
}

//________________________________________________________________________
AliPmdCustomCalibTask::AliPmdCustomCalibTask(const char *name) 
  :AliAnalysisTaskSE(name),
   fOutput(0), 
   fSmn(10),
   fHistAdcPre(0), 
   fHistAdcCpv(0),
   fHistClusterXYPre(0),           
   fHistClusterXYCpv(0) {
  for(Int_t i =0;i<48;i++) {
    for(Int_t j=0;j<96;j++)  {
      fHistAdcPreRC[i][j] = NULL;
      fHistAdcCpvRC[i][j] = NULL;
    }
  }
  DefineOutput(1, TList::Class());                                       }

//________________________________________________________________________
AliPmdCustomCalibTask::~AliPmdCustomCalibTask() {
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutput;
}

//________________________________________________________________________
void AliPmdCustomCalibTask::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetOwner(kTRUE);  
  
  // Module wise ADC
  fHistAdcPre = new TH1F("fHistAdcPre", "ADC Distribution for PRE module", 120, 0, 1200);
  fHistAdcCpv = new TH1F("fHistAdcCpv", "ADC Distribution for CPV module", 120, 0, 1200);
  
  // XY QA
  fHistClusterXYPre = new TH2F("fHistClusterXYPre", "ClusterX-Y Plot for PRE module", 400,-200,200,400,-200,200);
  fHistClusterXYCpv = new TH2F("fHistClusterXYCpv", "ClusterX-Y Plot for CPV module", 400,-200,200,400,-200,200);

  fOutput->Add(fHistAdcPre);
  fOutput->Add(fHistAdcCpv);
  fOutput->Add(fHistClusterXYPre);
  fOutput->Add(fHistClusterXYCpv);
  
  Char_t *name  = new Char_t[100];
  Char_t *title = new Char_t[200];

  // Cell-wise ADC
  for(Int_t i = 0; i < 48; i++)  {
    for(Int_t j = 0; j < 96; j++) {
      sprintf(name,"fHistPreAdcRaw%02dCol%02d",i,j);
      sprintf(title,"ADC Distribution for PRE | SM %d - Raw%02dCol%02d",fSmn,i,j);
      fHistAdcPreRC[i][j] = new TH1F(name,title,120,0,1200);
      fOutput->Add(fHistAdcPreRC[i][j]);
      
      sprintf(name,"fHistCpvAdcRaw%02dCol%02d",i,j);
      sprintf(title,"ADC Distribution for CPV | SM %d - Raw%02dCol%02d",fSmn,i,j);
      fHistAdcCpvRC[i][j] = new TH1F(name,title,120,0,1200);
      fOutput->Add(fHistAdcCpvRC[i][j]);
    }
  }	
 
  PostData(1, fOutput);   
}

//________________________________________________________________________
void AliPmdCustomCalibTask::UserExec(Option_t *) {
  AliVEvent *event = InputEvent();
  
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  if (!esd) {
    AliError("Cannot get the ESD event");
    return;
  }  
  
   

  Int_t npmdcl = esd->GetNumberOfPmdTracks();
  
  while (npmdcl--) {
    AliESDPmdTrack *pmdtr = esd->GetPmdTrack(npmdcl);

    Int_t   smn   = pmdtr->GetSmn();

    if (smn != fSmn) continue; // if not the SM go back

    Float_t adc   = pmdtr->GetClusterADC();
    if(adc>1200)  continue;

    Float_t sTag  = pmdtr->GetClusterSigmaX();
    if(sTag < 999. || sTag > 1000.) continue; // Isolated Cell Section
    
    Int_t   det   = pmdtr->GetDetector();
    
    Float_t clsX  = pmdtr->GetClusterX();
    Float_t clsY  = pmdtr->GetClusterY();

    Float_t sRowCol = pmdtr->GetClusterSigmaY();
    Int_t row = (Int_t)sRowCol/100;
    Int_t col = (Int_t)sRowCol - row*100;
    
    if(det == 0) {
      fHistAdcPre->Fill(adc);
      fHistAdcPreRC[row][col]->Fill(adc);  
      fHistClusterXYPre->Fill(clsX,clsY);    
    } else if(det == 1) {
      fHistAdcCpv->Fill(adc);	  
      fHistAdcCpvRC[row][col]->Fill(adc);
      fHistClusterXYCpv->Fill(clsX,clsY);    
    }
  }
  
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliPmdCustomCalibTask::Terminate(Option_t *) {
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
}
 
		  
