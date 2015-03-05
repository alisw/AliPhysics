/**************************************************************************
 * Author: Michal Broz                                               *
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskADQA class
//            This task is for QAing the AD data from ESD/AOD
//              Origin:michal.broz@cern.ch
//-----------------------------------------------------------------
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"

#include "AliAnalysisTaskADQA.h"

ClassImp(AliAnalysisTaskADQA)

//________________________________________________________________________
AliAnalysisTaskADQA::AliAnalysisTaskADQA() 
  : AliAnalysisTaskSE(),fListHist(0),fHistTotalChargePerEventADA(0),fHistTotalChargePerEventADC(0),fHistChargePerPM(0),fHistTimePerPM(0),fHistWidthPerPM(0),fHistTimeVsCharge(0),fHistWidthVsCharge(0),fHistNBBflagsADA(0),fHistNBBflagsADC(0),fHistNBBflagsADAVsADC(0),
    fHistChargeNoFlag(0),fHistTimeNoFlag(0), fHistChargeNoTime(0),
    fHistNCoincidencesADA(0),fHistNCoincidencesADC(0),fHistNCoincidencesADAVsADC(0)
{
  // Dummy constructor
}
//________________________________________________________________________
AliAnalysisTaskADQA::AliAnalysisTaskADQA(const char *name) 
  : AliAnalysisTaskSE(name),fListHist(0),fHistTotalChargePerEventADA(0),fHistTotalChargePerEventADC(0),fHistChargePerPM(0),fHistTimePerPM(0),fHistWidthPerPM(0),fHistTimeVsCharge(0),fHistWidthVsCharge(0),fHistNBBflagsADA(0),fHistNBBflagsADC(0),fHistNBBflagsADAVsADC(0),
    fHistChargeNoFlag(0),fHistTimeNoFlag(0),fHistChargeNoTime(0),
    fHistNCoincidencesADA(0),fHistNCoincidencesADC(0),fHistNCoincidencesADAVsADC(0)
{
  // Constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskADQA::~AliAnalysisTaskADQA(){
  // Destructor
  if (fListHist) { delete fListHist; fListHist = 0x0; }
}
//________________________________________________________________________
void AliAnalysisTaskADQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fListHist = new TList();
  fListHist->SetOwner(); 

if (!fHistTotalChargePerEventADA) {
    fHistTotalChargePerEventADA = CreateHist1D("fHistTotalChargePerEventADA","Total Charge in ADA per event",5000,0,5000,"ADC counts","Entries");
    fListHist->Add(fHistTotalChargePerEventADA);
  }
if (!fHistTotalChargePerEventADC) {
    fHistTotalChargePerEventADC = CreateHist1D("fHistTotalChargePerEventADC","Total Charge in ADC per event",5000,0,5000,"ADC counts","Entries");
    fListHist->Add(fHistTotalChargePerEventADC);
  }        
if (!fHistChargePerPM) {
    fHistChargePerPM = CreateHist2D("fHistChargePerPM","Charge per PM",16,-0.5,15.5,3000,0,3000,"PM number","ADC counts");
    fListHist->Add(fHistChargePerPM);
  } 
if (!fHistTimePerPM) {
    fHistTimePerPM = CreateHist2D("fHistTimePerPM","Time per PM",16,-0.5,15.5,1024, 240.039062, 340.039062,"PM number","Leading time [ns]");
    fListHist->Add(fHistTimePerPM);
  }  
if (!fHistWidthPerPM) {
    fHistWidthPerPM = CreateHist2D("fHistWidthPerPM","Width per PM",16,-0.5,15.5,128,0,100,"PM number","Time width [ns]");
    fListHist->Add(fHistWidthPerPM);
  }  
if (!fHistTimeVsCharge) {
    fHistTimeVsCharge = CreateHist2D("fHistTimeVsCharge","Time vs Charge",1024, 240.039062, 340.039062,3000,0,3000,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsCharge);
  }
if (!fHistWidthVsCharge) {
    fHistWidthVsCharge = CreateHist2D("fHistWidthVsCharge","Width vs Charge",128,0,100,3000,0,3000,"Time width [ns]","ADC counts");
    fListHist->Add(fHistWidthVsCharge);
  } 
if (!fHistNBBflagsADA) {
    fHistNBBflagsADA = CreateHist1D("fHistNBBflagsADA","Number of BB flags in ADA",9,-0.5,8.5,"Number of BB flags in ADA per event","Entries");
    fListHist->Add(fHistNBBflagsADA);
  } 
if (!fHistNBBflagsADC) {
    fHistNBBflagsADC = CreateHist1D("fHistNBBflagsADC","Number of BB flags in ADC",9,-0.5,8.5,"Number of BB flags in ADC per event","Entries");
    fListHist->Add(fHistNBBflagsADC);
  }
if (!fHistNBBflagsADAVsADC) {
    fHistNBBflagsADAVsADC = CreateHist2D("fHistNBBflagsADAVsADC","Number of BB flags",9,-0.5,8.5,9,-0.5,8.5,"Number of BB flags in ADA","Number of BB flags in ADC");
    fListHist->Add(fHistNBBflagsADAVsADC);
  }
  
if (!fHistNCoincidencesADA) {
    fHistNCoincidencesADA = CreateHist1D("fHistNCoincidencesADA","Number of BB coincidences in ADA",5,-0.5,4.5,"Number of BB coincidences in ADA per event","Entries");
    fListHist->Add(fHistNCoincidencesADA);
  } 
if (!fHistNCoincidencesADC) {
    fHistNCoincidencesADC = CreateHist1D("fHistNCoincidencesADC","Number of BB coincidences in ADC",5,-0.5,4.5,"Number of BB coincidences in ADC per event","Entries");
    fListHist->Add(fHistNCoincidencesADC);
  } 
if (!fHistNCoincidencesADAVsADC) {
    fHistNCoincidencesADAVsADC = CreateHist2D("fHistNCoincidencesADAVsADC","Number of BB coincidences",5,-0.5,4.5, 5,-0.5,4.5,"Number of BB coincidences in ADA","Number of BB coincidences in ADC");
    fListHist->Add(fHistNCoincidencesADAVsADC);
  }   
if (!fHistChargeNoFlag) {
    fHistChargeNoFlag = CreateHist1D("fHistChargeNoFlag","Charge in PM without BB flag",300,0,300,"Charge","Entries");
    fListHist->Add(fHistChargeNoFlag);
  } 
if (!fHistTimeNoFlag) {
    fHistTimeNoFlag = CreateHist1D("fHistTimeNoFlag","Time in PM without BB flag",400,0,400,"Time","Entries");
    fListHist->Add(fHistTimeNoFlag);
  }   
if (!fHistChargeNoTime) {
    fHistChargeNoTime = CreateHist1D("fHistChargeNoTime","Charge in PM without time measurement",300,0,300,"Charge","Entries");
    fListHist->Add(fHistChargeNoTime);
  } 

  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskADQA::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  AliVEvent* fEvent = InputEvent();
  if (!fEvent) {
    Printf("ERROR: Event not available");
    return;
  }
  AliESDEvent* fESD = (AliESDEvent*)fEvent;
  if (!fESD) {
    Printf("ERROR: ESD not available");
    return;
  }
  AliESDAD* esdAD = fESD->GetADData();
  if (!esdAD) {
    Printf("ERROR: No ESD AD");
    return;
  }
  /*/
  AliESDfriend *fESDfriend = fESD->FindFriend();
  if (!fESDfriend) {
    Printf("ERROR: No ESD friend");
    return;
  }
  
  AliESDADfriend* esdADfriend = fESDfriend->GetADfriend();
  if (!esdADfriend) {
    Printf("ERROR: No ESD AD friend");
    return;
  }
  /*/

  Float_t totChargeADA = 0;
  Float_t totChargeADC = 0;
  Int_t nBBflagsADA = 0;
  Int_t nBBflagsADC = 0;
  
  for(Int_t i=0; i<16; i++){
  	if(i<8)totChargeADC += esdAD->GetAdc(i);
	if(i>7)totChargeADA += esdAD->GetAdc(i);
	if(i<8 && esdAD->GetBBFlag(i)) nBBflagsADC++;
	if(i>7 && esdAD->GetBBFlag(i)) nBBflagsADA++;
	
	fHistChargePerPM->Fill(i,esdAD->GetAdc(i));
	fHistTimePerPM->Fill(i,esdAD->GetTime(i));
	fHistTimeVsCharge->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
	fHistWidthPerPM->Fill(i,esdAD->GetWidth(i));
	fHistWidthVsCharge->Fill(esdAD->GetWidth(i),esdAD->GetAdc(i));
	if(!esdAD->GetBBFlag(i)) {
		fHistChargeNoFlag->Fill(esdAD->GetAdc(i));
		fHistTimeNoFlag->Fill(esdAD->GetTime(i));
		}
	if(esdAD->GetTime(i) < 1e-6) fHistChargeNoTime->Fill(esdAD->GetAdc(i));	
	}
	
  Int_t nCoincidencesADA = 0;
  Int_t nCoincidencesADC = 0;
  for(Int_t i=0; i<4; i++){
  	if(esdAD->GetBBFlag(i) && esdAD->GetBBFlag(i+4)) nCoincidencesADC++;
	if(esdAD->GetBBFlag(i+8) && esdAD->GetBBFlag(i+12)) nCoincidencesADA++;
  	}
	
  fHistTotalChargePerEventADA->Fill(totChargeADA);
  fHistTotalChargePerEventADC->Fill(totChargeADC);
  fHistNBBflagsADA->Fill(nBBflagsADA);
  fHistNBBflagsADC->Fill(nBBflagsADC);
  fHistNBBflagsADAVsADC->Fill(nBBflagsADA,nBBflagsADC);
  fHistNCoincidencesADA->Fill(nCoincidencesADA);
  fHistNCoincidencesADC->Fill(nCoincidencesADC);
  fHistNCoincidencesADAVsADC->Fill(nCoincidencesADA,nCoincidencesADC);
  

  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskADQA::Terminate(Option_t *) 
{
 // Draw result to the screen
  // Called once at the end of the query

}
//________________________________________________________________________
TH1F * AliAnalysisTaskADQA::CreateHist1D(const char* name, const char* title,Int_t nBins, 
				    Double_t xMin, Double_t xMax,
				    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH1F* result = new TH1F(name, title, nBins, xMin, xMax);
  //result->SetOption("E");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  //result->SetMarkerStyle(kFullCircle);
  return result;
}
//________________________________________________________________________
TH2F * AliAnalysisTaskADQA::CreateHist2D(const char* name, const char* title,Int_t nBinsX, 
				    Double_t xMin, Double_t xMax,
				    Int_t nBinsY,
				    Double_t yMin, Double_t yMax,
				    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH2F* result = new TH2F(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  return result;
}

