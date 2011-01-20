
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena      .                                           *
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

/*-------------------------------------------------------------------------
 *
 *             AliEbyEFluctuationAnalysis base class
 *           This class deals with the Multiplicity   fluctuation 
 *                Origin: Satyajit Jena <sjena@cern.ch>
 *
 ------------------------------------------------------------------------*/

#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include "TH1F.h"
#include "TH2F.h"
#include <TH3I.h>

#include <TParticle.h>
#include <TList.h>

#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliMCParticle.h>
#include <AliPID.h>
#include <AliStack.h>
#include <AliCFContainer.h>
#include <AliCFEffGrid.h>
#include <AliCFDataGrid.h>
#include <AliTPCPIDResponse.h>
#include <AliESDpid.h>

#include "AliESDtrackCuts.h"
#include "AliEbyEEventBase.h"
#include "AliEbyEFluctuationAnalysis.h"

ClassImp(AliEbyEFluctuationAnalysis)

//____________________________________________________________________//
AliEbyEFluctuationAnalysis::AliEbyEFluctuationAnalysis() : 
  TObject(), 
  fListMeasure(0),
  fEbyEBase(0),
  fESDtrackCuts(0)
{
  //Default constructor
  InitHistos();
}

//____________________________________________________________________//
AliEbyEFluctuationAnalysis::~AliEbyEFluctuationAnalysis() {
  //Default destructor
  if(fListMeasure) delete fListMeasure;
  if(fEbyEBase) delete fEbyEBase;
  if(fESDtrackCuts) delete fESDtrackCuts;
}

//____________________________________________________________________//
void AliEbyEFluctuationAnalysis::InitHistos() {
  fListMeasure = new TList();
  fListMeasure->SetName("MeasuresList");
  Int_t cent = 50;

  // cent = fEbyEBase->GetCentralityBin();  

  TH1F *fEvtStat = new TH1F("hEventStatistuc","Analysis Event Counter", 10,0,10);
  fListMeasure->Add(fEvtStat); // 
  
  TH2F *fCvsM = new TH2F("hCentralityVsMultiplisity", "Correlation Plots; Centrality; Multiplicity",
					     6000,0,6000,cent,0,(Double_t)cent);
  fListMeasure->Add(fCvsM);  // Correlation Histo centrality and Multiplicity
  
  TH2F *fCharge[cent];
  
  char hname[40]; char htitle[100]; 
  Int_t binm = 100/cent;
  
  for(Int_t i = 0; i < cent; i++ ) {  // Preparing Histograms for Each Centrality Bin
    
    sprintf(hname,"hNPinCent%02d",i);
    sprintf(htitle,"N_{+} vs N_{-} in Centrality %02d%% - %02d%%; N_{+}; N_{-}\n",i*binm, (i+1)*binm);
    
    fCharge[i] = new TH2F(hname,htitle,3000,0.5,3000.5,3000,0.5,3000.5);
    fListMeasure->Add(fCharge[i]); 
  }
  
  
}


//____________________________________________________________________//
void AliEbyEFluctuationAnalysis::Analyze(AliESDEvent* esd) {

  // Analysis from ESD
  ((TH1F*)(fListMeasure->At(0)))->Fill(1);  
  
  Int_t kCent = 1;
   kCent =  fEbyEBase->GetCentrality(esd);
   
  TObjArray *array = new TObjArray();
  Int_t itr = esd->GetNumberOfTracks();
  
  for (Int_t iTracks = 0; iTracks < itr; iTracks++) 
    {
      AliESDtrack* track = esd->GetTrack(iTracks);
      if (!track) continue;
      if(!fESDtrackCuts->AcceptTrack(track)) continue;
      array->Add(track);
    }

  Calculate(array,kCent);

  delete array;
}


//____________________________________________________________________//
void AliEbyEFluctuationAnalysis::Analyze(AliAODEvent* aod) {
  // Analysis from AOD
 
  printf("%d\n",aod->GetNTracks());
  

}


//____________________________________________________________________//
void AliEbyEFluctuationAnalysis::Analyze(AliStack* stack) {
  // Analysis from MC stack


  printf("%d \n",stack->GetNtrack());

}

//____________________________________________________________________________//
void AliEbyEFluctuationAnalysis::Calculate(TObjArray *gTrackArray, Int_t cent){

  ((TH1F*)(fListMeasure->At(0)))->Fill(1,1); 
  if(cent < 0) return;
  ((TH1F*)(fListMeasure->At(0)))->Fill(2,1); 
  Int_t gNtrack = gTrackArray->GetEntries();
  if(gNtrack < 1) return;
  ((TH1F*)(fListMeasure->At(0)))->Fill(3,1); 

  AliVParticle* track = 0;
  Int_t pCharge = 0;
  Int_t nCharge = 0;
  ((TH1F*)(fListMeasure->At(1)))->Fill(3,1); 
  
  TString gAnalysisLevel = fEbyEBase->GetAnalysisLevel();  
  for(Int_t i = 0; i < gNtrack; i++) {
    if(gAnalysisLevel == "ESD")
      track = dynamic_cast<AliESDtrack*>(gTrackArray->At(i));
    Short_t charge = track->Charge();
    if(charge > 0) pCharge += 1; 
    if(charge < 0) nCharge += 1; 
  }
  
  ((TH2F*)(fListMeasure->At(2)))->Fill(cent,pCharge, nCharge); //+ vs -
 
  
}
