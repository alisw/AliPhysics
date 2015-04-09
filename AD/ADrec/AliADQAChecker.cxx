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


/*
  Checks the quality assurance. Under construction. 
  By comparing with reference data

*/

// --- ROOT system ---
#include <TClass.h>
#include <TH2F.h> 
#include <TH1D.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliADQAChecker.h"
#include "AliADQADataMakerRec.h"


ClassImp(AliADQAChecker)

//__________________________________________________________________
AliADQAChecker::AliADQAChecker() : AliQACheckerBase("AD","AD Quality Assurance Data Checker"),
  fLowEventCut(1000),
  fORvsANDCut(0.2),
  fBGvsBBCut(0.2)
{
  // Default constructor
  // Nothing else here
}

//__________________________________________________________________
void AliADQAChecker::Check(Double_t * check, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/) 
{
  // Main check function: Depending on the TASK, different checks will be applied
  // Check for missing channels and check on the trigger type for raw data
  // Check for missing disk or rings for esd (to be redone)

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    check[specie] = 1.0;
    // no check on calibration events
    if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib)continue;
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue;
    if (index == AliQAv1::kRAW) {
      check[specie] =  CheckRaws(list[specie]);
    } else if (index == AliQAv1::kESD) {
      check[specie] =  CheckEsds(list[specie]);
    } 
  }
}

//_________________________________________________________________
Double_t AliADQAChecker::CheckRaws(TObjArray * list) const
{

  //  Check on the QA histograms on the raw-data input list:
  //  Two things are checked: the presence of data in all channels and
  //  the ratio between different trigger types

  Double_t test = 1.0;
  if (list->GetEntries() == 0){  
    AliWarning("There are no histograms to be checked");
  }
  else {  
    TH2F *hChargeEoIInt0 = (TH2F*)list->At(AliADQADataMakerRec::kChargeEoIInt0);
    if (!hChargeEoIInt0) {
      AliWarning("ChargeEoIInt0 histogram is not found");
    }
    else {
    	if(hChargeEoIInt0->GetListOfFunctions()->GetEntries()<1) hChargeEoIInt0->GetListOfFunctions()->Add(new TPaveText(0.5,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hChargeEoIInt0->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hChargeEoIInt0->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hChargeEoIInt0->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hChargeSlice;
			Int_t NbadChannels = 0;
			TString badChannels = " ";
			
    			for(Int_t i=0; i<16; i++){
      				hChargeSlice = hChargeEoIInt0->ProjectionY("hChargeSlice",i+1,i+1);
				Double_t saturation = hChargeSlice->Integral(1000,1025);
				if(saturation > 1000){
					test = 0.1;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kRed);
						QAbox->AddText("Saturation on channel ");
						}
					badChannels += i;
					badChannels += ", ";
					NbadChannels++;
					}
				if(NbadChannels != 0 && i==15) QAbox->AddText(badChannels.Data());
				}
			if(NbadChannels == 0){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("OK");
				}		
    			}
		}
    	}
    TH2F *hChargeEoIInt1 = (TH2F*)list->At(AliADQADataMakerRec::kChargeEoIInt1);
    if (!hChargeEoIInt1) {
      AliWarning("ChargeEoIInt1 histogram is not found");
    }
    else {
    	if(hChargeEoIInt1->GetListOfFunctions()->GetEntries()<1) hChargeEoIInt1->GetListOfFunctions()->Add(new TPaveText(0.5,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hChargeEoIInt1->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hChargeEoIInt1->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hChargeEoIInt1->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hChargeSlice;
			Int_t NbadChannels = 0;
			TString badChannels = " ";
			
    			for(Int_t i=0; i<16; i++){
      				hChargeSlice = hChargeEoIInt1->ProjectionY("hChargeSlice",i+1,i+1);
				Double_t saturation = hChargeSlice->Integral(1000,1025);
				if(saturation > 1000){
					test = 0.1;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kRed);
						QAbox->AddText("Saturation on channel ");
						}
					badChannels += i;
					badChannels += ", ";
					NbadChannels++;
					}
				if(NbadChannels != 0 && i==15) QAbox->AddText(badChannels.Data());
				}
			if(NbadChannels == 0){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("OK");
				}		
    			}
		}
    	}  
    TH2F *hBBFlagVsClock = (TH2F*)list->At(AliADQADataMakerRec::kBBFlagVsClock);
    if (!hBBFlagVsClock) {
      AliWarning("BBFlagVsClock histogram is not found");
    }
    else {
    	if(hBBFlagVsClock->GetListOfFunctions()->GetEntries()<1) hBBFlagVsClock->GetListOfFunctions()->Add(new TPaveText(0.5,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hBBFlagVsClock->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hBBFlagVsClock->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hBBFlagVsClock->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hClockSlice;
			Bool_t notSynch = kFALSE;

    			for(Int_t i=0; i<16; i++){
      				hClockSlice = hBBFlagVsClock->ProjectionY("hClockSlice",i+1,i+1);
				Double_t center = hClockSlice->GetBinContent(11);
				Double_t around = hClockSlice->Integral(0,10) + hClockSlice->Integral(12,21);
				if(center == 0){
					test = 0.1;
					QAbox->SetFillColor(kRed);
					QAbox->AddText("AD not synchronized");
					notSynch = kTRUE;
					break;
					} 
				if(around/center > 0.5) {
					test = 0.1;
					test = 0.1;
					QAbox->SetFillColor(kRed);
					QAbox->AddText("AD not synchronized");
					notSynch = kTRUE;
					break;
					}
				}
			if(!notSynch){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("OK");
				}		
    			}
		}
    	} 
    TH2F *hBGFlagVsClock = (TH2F*)list->At(AliADQADataMakerRec::kBGFlagVsClock);
    if (!hBGFlagVsClock) {
      AliWarning("BGFlagVsClock histogram is not found");
    }
    else {
    	if(hBGFlagVsClock->GetListOfFunctions()->GetEntries()<1) hBGFlagVsClock->GetListOfFunctions()->Add(new TPaveText(0.5,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hBGFlagVsClock->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hBGFlagVsClock->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hBGFlagVsClock->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hClockSlice;
			Bool_t notSynch = kFALSE;

    			for(Int_t i=0; i<16; i++){
      				hClockSlice = hBGFlagVsClock->ProjectionY("hClockSlice",i+1,i+1);
				Double_t center = hClockSlice->GetBinContent(11);
				Double_t around = hClockSlice->Integral(0,10) + hClockSlice->Integral(12,21);
				if(center == 0){
					test = 0.1;
					QAbox->SetFillColor(kRed);
					QAbox->AddText("AD not synchronized");
					notSynch = kTRUE;
					break;
					} 
				if(around/center > 0.5) {
					test = 0.1;
					test = 0.1;
					QAbox->SetFillColor(kRed);
					QAbox->AddText("AD not synchronized");
					notSynch = kTRUE;
					break;
					}
				}
			if(!notSynch){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("OK");
				}		
    			}
		}
    	}  
     
    TH2F *hPedestalDiffInt0  = (TH2F*)list->At(AliADQADataMakerRec::kPedestalDiffInt0);
    if (!hPedestalDiffInt0) {
      AliWarning("PedestalInt0 histogram is not found");
    }
    else {
    	if(hPedestalDiffInt0->GetListOfFunctions()->GetEntries()<1) hPedestalDiffInt0->GetListOfFunctions()->Add(new TPaveText(0.5,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hPedestalDiffInt0->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hPedestalDiffInt0->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hPedestalDiffInt0->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hPedestalSlice;
			Int_t NbadChannels = 0;
			TString badChannels = " ";
    			for(Int_t i=0; i<16; i++){
      				hPedestalSlice = hPedestalDiffInt0->ProjectionY("hPedestalSlice",i+1,i+1);
				Double_t mean = hPedestalSlice->GetMean();
				if(TMath::Abs(mean)>1) {
					test = 0.1;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kRed);
						QAbox->AddText("Bad pedestal for channel ");
						}
					badChannels += i;
					badChannels += ", ";
					NbadChannels++;
					}
				if(NbadChannels != 0 && i==15) QAbox->AddText(badChannels.Data());
				}
			if(NbadChannels == 0){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("OK");
				}		
    			}
		}
    	}
    TH2F *hPedestalDiffInt1  = (TH2F*)list->At(AliADQADataMakerRec::kPedestalDiffInt1);
    if (!hPedestalDiffInt1) {
      AliWarning("PedestalInt1 histogram is not found");
    }
    else {
    	if(hPedestalDiffInt1->GetListOfFunctions()->GetEntries()<1) hPedestalDiffInt1->GetListOfFunctions()->Add(new TPaveText(0.5,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hPedestalDiffInt1->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hPedestalDiffInt1->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hPedestalDiffInt1->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hPedestalSlice;
			Int_t NbadChannels = 0;
			TString badChannels = " ";
    			for(Int_t i=0; i<16; i++){
      				hPedestalSlice = hPedestalDiffInt1->ProjectionY("hPedestalSlice",i+1,i+1);
				Double_t mean = hPedestalSlice->GetMean();
				if(TMath::Abs(mean)>1) {
					test = 0.1;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kRed);
						QAbox->AddText("Bad pedestal for channel ");
						}
					badChannels += i;
					badChannels += ", ";
					NbadChannels++;
					}
				if(NbadChannels != 0 && i==15) QAbox->AddText(badChannels.Data());
				}
			if(NbadChannels == 0){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("OK");
				}		
    			}
		}
    	}
  } 
  return test ;  
}  

//_________________________________________________________________
Double_t AliADQAChecker::CheckEsds(TObjArray * list) const
{
  
//  check the ESDs for missing disk or ring
//  printf(" Number of entries in ESD list = %d\n", list->GetEntries()); 
//  list->Print();

  Double_t test     = 1.0;     // initialisation to OK

  return test ; 
} 

//______________________________________________________________________________
void AliADQAChecker::Init(const AliQAv1::DETECTORINDEX_t det) 
{
  // intialises QA and QA checker settings
  AliQAv1::Instance(det) ; 
  Float_t * hiValue = new Float_t[AliQAv1::kNBIT] ; 
  Float_t * lowValue = new Float_t[AliQAv1::kNBIT] ;
  lowValue[AliQAv1::kINFO]      = 0.5   ; 
  hiValue[AliQAv1::kINFO]       = 1.0 ; 
  lowValue[AliQAv1::kWARNING]   = 0.2 ; 
  hiValue[AliQAv1::kWARNING]    = 0.5 ; 
  lowValue[AliQAv1::kERROR]     = 0.0   ; 
  hiValue[AliQAv1::kERROR]      = 0.2 ; 
  lowValue[AliQAv1::kFATAL]     = -1.0   ; 
  hiValue[AliQAv1::kFATAL]      = 0.0 ; 
  SetHiLo(hiValue, lowValue) ; 
  delete [] hiValue;
  delete [] lowValue;
}

//______________________________________________________________________________
void AliADQAChecker::SetQA(AliQAv1::ALITASK_t index, Double_t * value) const
{
// sets the QA word according to return value of the Check
  AliQAv1 * qa = AliQAv1::Instance(index);
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    qa->UnSet(AliQAv1::kFATAL, specie);
    qa->UnSet(AliQAv1::kWARNING, specie);
    qa->UnSet(AliQAv1::kERROR, specie);
    qa->UnSet(AliQAv1::kINFO, specie);
    if ( ! value ) { // No checker is implemented, set all QA to Fatal
      qa->Set(AliQAv1::kFATAL, specie) ; 
    } else {
      if ( value[specie] >= fLowTestValue[AliQAv1::kFATAL] && value[specie] < fUpTestValue[AliQAv1::kFATAL] ) 
        qa->Set(AliQAv1::kFATAL, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kERROR] && value[specie] <= fUpTestValue[AliQAv1::kERROR]  )
        qa->Set(AliQAv1::kERROR, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kWARNING] && value[specie] <= fUpTestValue[AliQAv1::kWARNING]  )
        qa->Set(AliQAv1::kWARNING, specie) ;
      else if ( value[specie] > fLowTestValue[AliQAv1::kINFO] && value[specie] <= fUpTestValue[AliQAv1::kINFO] ) 
        qa->Set(AliQAv1::kINFO, specie) ; 	
    }
  }
}

//____________________________________________________________________________ 
void AliADQAChecker::MakeImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) 
{
  // makes the QA image for sim and rec
  TObjArray tmpArr;  // array to store flat version of original array (which may contain clones)
  //
  for (Int_t esIndex = 0; esIndex < AliRecoParam::kNSpecies; esIndex++) {
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(esIndex)) || list[esIndex]->GetEntries() == 0) continue;
    Int_t nImages = 0;
    TIter next(list[esIndex]);
    TObject* hdata = NULL;
    tmpArr.Clear();
    while ( (hdata=(next())) ) { // count histos and transfere to flat array
      if (hdata->InheritsFrom(TH1::Class()) && hdata->TestBit(AliQAv1::GetImageBit()) ) {  // histo, not cloned
	nImages++; 
	tmpArr.AddLast(hdata); 
	continue;
      }
      if (!hdata->TestBit(AliQAv1::GetClonedBit())) continue;  // not an array of clones, unknown object
      TIter nextCl((TObjArray*)hdata);   // array of histo clones
      TObject* hcl = 0;
      while ((hcl=nextCl())) if (hcl->InheritsFrom(TH1::Class()) && hcl->TestBit(AliQAv1::GetImageBit())) {tmpArr.AddLast(hcl); nImages++;}
    }
    //
    if ( nImages == 0 ) {
      AliDebug(AliQAv1::GetQADebugLevel(), Form("No histogram will be plotted for %s %s %s\n", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(esIndex)));  
      continue;
    }
    AliDebug(AliQAv1::GetQADebugLevel(), Form("%d histograms will be plotted for %s %s %s\n", nImages, GetName(), AliQAv1::GetTaskName(task).Data(),AliRecoParam::GetEventSpecieName(esIndex)));  
    //        
    const Char_t * title = Form("QA_%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(esIndex)); 
    //
    if ( !fImage[esIndex] ) fImage[esIndex] = new TCanvas(title, title,1000,1500);
    //
    fImage[esIndex]->Clear(); 
    fImage[esIndex]->SetTitle(title); 
    fImage[esIndex]->cd(); 
    TPaveText someText(0.015, 0.015, 0.98, 0.98) ;
    someText.AddText(title) ;
    someText.SetFillColor(0);
    someText.SetFillStyle(0);
    someText.SetBorderSize(0);
    someText.SetTextColor(kRed+1);
    someText.Draw() ; 
    TString outName(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()));
    fImage[esIndex]->Print(outName, "ps") ; 

    // Now set some parameters on the canvas 
    fImage[esIndex]->Clear(); 
    
    const char* topT = Form("%s, %s, Run: %d", AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(esIndex),AliQAChecker::Instance()->GetRunNumber());
    TLatex* topText = new TLatex(0.5, 0.99, topT);
    topText->SetTextAlign(23);
    topText->SetTextSize(.038);
    topText->SetTextFont(42);
    topText->SetTextColor(kBlue+3);
    topText->SetNDC();
    topText->Draw();
    
    if(esIndex != AliRecoParam::kCalib){

    	TVirtualPad* pCharge = 0;
    	TVirtualPad* pTime  = 0;
    	TVirtualPad* pClock  = 0;
    	TVirtualPad* pCoinc  = 0;
    	TVirtualPad* pPed  = 0;
    	TVirtualPad* pMaxCh  = 0;
	TVirtualPad* pFlagCorr  = 0;
    	TVirtualPad* pTrigger  = 0;

    	TPad* pCh = new TPad("Charge", "Charge Pad", 0, 0.83, 1.0, 0.95);
    	fImage[esIndex]->cd();
    	pCh->Draw();
    	pCharge = pCh;

    	TPad* pT = new TPad("Time", "Time Pad", 0, 0.59, 1.0, 0.83);
    	fImage[esIndex]->cd();
    	pT->Draw();
    	pTime = pT;
    
    	TPad* pCl = new TPad("Clock", "Clock Pad", 0, 0.47, 1.0, 0.59);
    	fImage[esIndex]->cd();
    	pCl->Draw();
    	pClock = pCl;
    
    	TPad* pCo = new TPad("Coincidences", "Coincidences Pad", 0, 0.23, 1.0, 0.35);
    	fImage[esIndex]->cd();
    	pCo->Draw();
    	pCoinc = pCo;
    
    	TPad* pP = new TPad("Pedestal", "Pedestal Pad", 0, 0.35, 0.5, 0.47);
    	fImage[esIndex]->cd();
    	pP->Draw();
    	pPed = pP;
    
    	TPad* pM = new TPad("Max Charge", "Max Charge Pad", 0.5, 0.35, 1.0, 0.47);
    	fImage[esIndex]->cd();
    	pM->Draw();
    	pMaxCh = pM;
	
	TPad* pFC = new TPad("Flag correlations", "Flag correlations Pad", 0.25, 0.11, 0.75, 0.23);
    	fImage[esIndex]->cd();
    	pFC->Draw();
    	pFlagCorr = pFC;
	
	TPad* pTr = new TPad("Triggers", "Triggers Pad", 0, 0.0, 1.0, 0.11);
    	fImage[esIndex]->cd();
    	pTr->Draw();
    	pTrigger = pTr;

    	pCharge->Divide(3, 1);
    	pTime->Divide(4, 1);
    	pClock->Divide(4, 1);
    	pCoinc->Divide(4, 1);
    	pPed->Divide(2, 1);
    	pMaxCh->Divide(2, 1);
	pFlagCorr->Divide(2, 1);
	
  
    	TIter nexthist(&tmpArr);
    	Int_t npad = 1; 
    	TH1* histo = 0;
    	while ( npad < 23) { // tmpArr is guaranteed to contain only plottable histos, no checks needed
      		histo=(TH1*)nexthist();
      		histo->SetStats(kFALSE);
      		TVirtualPad* pad = 0;
       
      		if(npad<4) pad = pCharge->cd(npad);
      		if(npad>3 && npad<8) pad = pTime->cd(npad-3);
      		if(npad>7 && npad<12) pad = pClock->cd(npad-7);
      		if(npad>11 && npad<16) pad = pCoinc->cd(npad-11);
      		if(npad>15 && npad<18) pad = pPed->cd(npad-15);
      		if(npad>17 && npad<20) pad = pMaxCh->cd(npad-17);
		if(npad>19 && npad<22) pad = pFlagCorr->cd(npad-19);
		if(npad == 22) pad = pTrigger->cd();
      
      		pad->SetRightMargin(0.10);
      		pad->SetLeftMargin(0.10);
      		pad->SetBottomMargin(0.10);
      
      		if(npad ==1 || npad==2 || npad==12 || npad==13 || npad==14 || npad==15) gPad->SetLogy();
      		if(npad ==3 || npad ==18 || npad ==19 || npad ==20 || npad ==21) gPad->SetLogz();
      		histo->DrawCopy("colz");  
     
      		npad++; 
    	}
    }
    else{
	TVirtualPad* pPed  = 0;
	
    	TPad* pP = new TPad("Pedestal", "Pedestal Pad", 0.0, 0.1, 1.0, 0.95);
    	fImage[esIndex]->cd();
    	pP->Draw();
    	pPed = pP;
	pPed->Divide(1, 2);
	
	
	TIter nexthist(&tmpArr);
    	Int_t npad = 1; 
    	TH1* histo = 0;
    	while ( npad < 25) { // tmpArr is guaranteed to contain only plottable histos, no checks needed
      		histo=(TH1*)nexthist();
      		histo->SetStats(kFALSE);
      		TVirtualPad* pad = 0;
        	if((npad>22 && npad<25)) {
			pad = pPed->cd(npad-22);
      			pad->SetRightMargin(0.10);
      			pad->SetLeftMargin(0.10);
      			pad->SetBottomMargin(0.10);
      			histo->DrawCopy("colz"); 
			} 
      		npad++; 
    	}
    }
    

    fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps"); 
  }
}
 
