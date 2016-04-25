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
#include <TLegend.h>
#include <TStyle.h>
#include <TSpline.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliADQAChecker.h"
#include "AliADQADataMakerRec.h"
#include "AliADQAParam.h"


ClassImp(AliADQAChecker)

//__________________________________________________________________
AliADQAChecker::AliADQAChecker() : AliQACheckerBase("AD","AD Quality Assurance Data Checker"),
  fLowEventCut(1000),
  fORvsANDCut(0.2),
  fBGvsBBCut(0.2),
  fSatMed(0.1),
  fSatHigh(0.3),
  fSatHuge(0.5),
  fMaxPedDiff(1),
  fMaxPedWidth(1.5),
  fChargeChannelZoomMin(0),
  fChargeChannelZoomMax(50),
  fTimeRatioBBZoomMin(170),
  fTimeRatioBBZoomMax(210),
  fTimeRatioBGZoomMin(50),
  fTimeRatioBGZoomMax(90),
  fChargeTrendMin(0),
  fChargeTrendMax(1000),
  fMaxNoTimeRate(10e-3),
  fMaxNoFlagRate(10e-2),
  fMaxBBVariation(5),
  fMaxBGVariation(5),
  fAsynchronBB(0.5),
  fAsynchronBG(0.5)
{
  fQAParam = (AliADQAParam*)GetQAParam();
  fSatMed = fQAParam->GetSatMed();
  fSatHigh = fQAParam->GetSatHigh();
  fSatHuge = fQAParam->GetSatHuge();
  fMaxPedDiff = fQAParam->GetMaxPedDiff();
  fMaxPedWidth = fQAParam->GetMaxPedWidth();
  fChargeChannelZoomMin = fQAParam->GetChargeChannelZoomMin();
  fChargeChannelZoomMax = fQAParam->GetChargeChannelZoomMax();
  fTimeRatioBBZoomMin =  fQAParam->GetTdcTimeMinBBFlag();
  fTimeRatioBBZoomMax =  fQAParam->GetTdcTimeMaxBBFlag();
  fTimeRatioBGZoomMin =  fQAParam->GetTdcTimeMinBGFlag();
  fTimeRatioBGZoomMax =  fQAParam->GetTdcTimeMaxBGFlag();
  fChargeTrendMin = fQAParam->GetChargeTrendMin();
  fChargeTrendMax = fQAParam->GetChargeTrendMax();
  fMaxNoTimeRate = fQAParam->GetMaxNoTimeRate();
  fMaxNoFlagRate = fQAParam->GetMaxNoFlagRate();
  fMaxBBVariation = fQAParam->GetMaxBBVariation();
  fMaxBGVariation = fQAParam->GetMaxBGVariation();
  fAsynchronBB = fQAParam->GetAsynchronBB();
  fAsynchronBG = fQAParam->GetAsynchronBG();
  
}

//____________________________________________________________________________
AliADQAParam* AliADQAChecker::GetQAParam() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("AD/Calib/QAParam");
  if(!entry)AliWarning("Load of QA param from default storage failed!");
 
  // Retrieval of data in directory AD/Calib/QA:

  AliADQAParam *QAParam = 0;

  if (entry) QAParam = (AliADQAParam*) entry->GetObject();
  if (!QAParam)  AliFatal("No QA param from calibration database !");

  return QAParam;
}
//__________________________________________________________________
void AliADQAChecker::Check(Double_t * check, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/) 
{
  // Main check function: Depending on the TASK, different checks will be applied
  // Check for missing channels and check on the trigger type for raw data
  // Check for missing disk or rings for esd (to be redone)

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    check[specie] = 1.0;    
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue;
    if (index == AliQAv1::kRAW) {
      if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib) check[specie] =  CheckPedestals(list[specie]);
      else check[specie] =  CheckRaws(list[specie]);
    } 
    else if (index == AliQAv1::kESD) {
      check[specie] =  CheckEsds(list[specie]);
    } 
  }
}
//_________________________________________________________________
Double_t AliADQAChecker::CheckPedestals(TObjArray * list) const
{
  //  Check on the QA histograms on the raw-data input list:

Double_t test = 1.0;
  if (list->GetEntries() == 0){  
    AliWarning("There are no histograms to be checked");
  }
  else {
    TH2F *hPedestalInt0  = (TH2F*)list->At(AliADQADataMakerRec::kPedestalInt0);
    if (!hPedestalInt0) {
      AliWarning("PedestalInt0 histogram is not found");
    }
    else {
    	if(hPedestalInt0->GetListOfFunctions()->GetEntries()<1) hPedestalInt0->GetListOfFunctions()->Add(new TPaveText(0.15,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hPedestalInt0->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hPedestalInt0->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hPedestalInt0->GetListOfFunctions()->At(i);
      			QAbox->Clear();
			if(hPedestalInt0->GetEntries() == 0){
				QAbox->Clear();
        			QAbox->SetFillColor(kYellow);
        			QAbox->AddText("Is this AD pedestal run?");
				}
			else{
    				TH1D *hPedestalSlice;
				Int_t NbadChannels = 0;
				TString badChannels = " ";
    				for(Int_t i=0; i<16; i++){
      					hPedestalSlice = hPedestalInt0->ProjectionY("hPedestalSlice",i+1,i+1);
					Double_t RMS = hPedestalSlice->GetRMS();
					if(RMS > fMaxPedWidth) {
						test = 0.1;
						if(NbadChannels == 0){
							QAbox->SetFillColor(kRed);
							QAbox->AddText("Noisy pedestal for channel ");
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
    TH2F *hPedestalInt1  = (TH2F*)list->At(AliADQADataMakerRec::kPedestalInt1);
    if (!hPedestalInt1) {
      AliWarning("PedestalInt1 histogram is not found");
    }
    else {
    	if(hPedestalInt1->GetListOfFunctions()->GetEntries()<1) hPedestalInt1->GetListOfFunctions()->Add(new TPaveText(0.15,0.63,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hPedestalInt1->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hPedestalInt1->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hPedestalInt1->GetListOfFunctions()->At(i);
      			QAbox->Clear();
			if(hPedestalInt1->GetEntries() == 0){
				QAbox->Clear();
        			QAbox->SetFillColor(kYellow);
        			QAbox->AddText("Is this AD pedestal run?");
				}
			else{
    				TH1D *hPedestalSlice;
				Int_t NbadChannels = 0;
				TString badChannels = " ";
    				for(Int_t i=0; i<16; i++){
      					hPedestalSlice = hPedestalInt1->ProjectionY("hPedestalSlice",i+1,i+1);
					Double_t RMS = hPedestalSlice->GetRMS();
					if(RMS > fMaxPedWidth) {
						test = 0.1;
						if(NbadChannels == 0){
							QAbox->SetFillColor(kRed);
							QAbox->AddText("Noisy pedestal for channel ");
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
  }
  
  return test ; 
}

//_________________________________________________________________
Double_t AliADQAChecker::CheckRaws(TObjArray * list) const
{

  //  Check on the QA histograms on the raw-data input list:

  Double_t test = 1.0;
  if (list->GetEntries() == 0){  
    AliWarning("There are no histograms to be checked");
  }
  else {
    Float_t nEvents;
    TH1F *hChargeADA = (TH1F*)list->At(AliADQADataMakerRec::kChargeADA);
    if (!hChargeADA) {
      AliWarning("ChargeADA histogram is not found");
    }
    else {
    	nEvents = hChargeADA->Integral(-1,-1);
    	if(hChargeADA->GetListOfFunctions()->GetEntries()<1) hChargeADA->GetListOfFunctions()->Add(new TPaveText(0.43,0.70,0.66,0.83,"NDC"));
    	for(Int_t i=0; i<hChargeADA->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hChargeADA->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hChargeADA->GetListOfFunctions()->At(i);
      			
			TString meanChargeADA = " Mean charge ADA = ";
			meanChargeADA += TString::Format("%4.1f ", hChargeADA->GetMean());
			
			QAbox->Clear();
			QAbox->SetFillColor(kWhite);
        		QAbox->SetTextColor(kBlue);
			QAbox->AddText(meanChargeADA.Data());				
    			}
		}
    	}
	
    TH1F *hChargeADC = (TH1F*)list->At(AliADQADataMakerRec::kChargeADC);
    if (!hChargeADC) {
      AliWarning("ChargeADC histogram is not found");
    }
    else {
    	if(hChargeADC->GetListOfFunctions()->GetEntries()<1) hChargeADC->GetListOfFunctions()->Add(new TPaveText(0.43,0.55,0.66,0.68,"NDC"));
    	for(Int_t i=0; i<hChargeADC->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hChargeADC->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hChargeADC->GetListOfFunctions()->At(i);
      			
			TString meanChargeADC = " Mean charge ADC = ";
			meanChargeADC += TString::Format("%4.1f ", hChargeADC->GetMean());
			
			QAbox->Clear();
			QAbox->SetFillColor(kWhite);
        		QAbox->SetTextColor(kRed);
			QAbox->AddText(meanChargeADC.Data());				
    			}
		}
    	}
    TH1F *hTriggers = (TH1F*)list->At(AliADQADataMakerRec::kTriggers);
    if (!hTriggers) {
      AliWarning("Triggers histogram is not found");
    }
    else {
    	if(hTriggers->GetListOfFunctions()->GetEntries()<1) hTriggers->GetListOfFunctions()->Add(new TPaveText(0.56,0.73,0.72,0.85,"NDC"));
    	for(Int_t i=0; i<hTriggers->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hTriggers->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hTriggers->GetListOfFunctions()->At(i);
      			
			TString nEventsText = " Number of events = ";
			nEventsText += TString::Format("%.0f ", nEvents);
			
			QAbox->Clear();
			QAbox->SetFillColor(kWhite);
        		QAbox->SetTextColor(kAzure-8);
			QAbox->AddText(nEventsText.Data());				
    			}
		}
    	}
    TH1F *hFlagNoTime = (TH1F*)list->At(AliADQADataMakerRec::kFlagNoTime);
    if (!hFlagNoTime) {
      AliWarning("FlagNoTime histogram is not found");
    }
    else {
    	if(hFlagNoTime->GetListOfFunctions()->GetEntries()<1) hFlagNoTime->GetListOfFunctions()->Add(new TPaveText(0.15,0.67,0.85,0.76,"NDC"));
    	for(Int_t i=0; i<hFlagNoTime->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hFlagNoTime->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hFlagNoTime->GetListOfFunctions()->At(i);
      			
			TH1F *histoRate = (TH1F*)hFlagNoTime->Clone("histoRate");
			histoRate->Sumw2();
			if(nEvents != 0) histoRate->Scale(1/nEvents);
			Int_t NbadChannels = 0;
			TString badChannels = "No time rate too high, ch:";
			
			for(Int_t i=0; i<16; i++){
				if(histoRate->GetBinContent(i+1)>fMaxNoTimeRate) {
					test = 0.3;
					badChannels += i;
					badChannels += ", ";
					NbadChannels++;
					}
				}
			if(NbadChannels != 0){
				QAbox->Clear();
				QAbox->SetFillColor(kRed); 
				QAbox->AddText(badChannels.Data());
				}
			else{
				QAbox->Clear();
				QAbox->SetFillColor(kGreen);
				QAbox->AddText("No time rate OK");				
    				}
			}
    		}
    	}
    TH1F *hTimeNoFlag = (TH1F*)list->At(AliADQADataMakerRec::kTimeNoFlag);
    if (!hTimeNoFlag) {
      AliWarning("TimeNoFlag histogram is not found");
    }
    else {
    	if(hTimeNoFlag->GetListOfFunctions()->GetEntries()<1) hTimeNoFlag->GetListOfFunctions()->Add(new TPaveText(0.15,0.56,0.85,0.65,"NDC"));
    	for(Int_t i=0; i<hTimeNoFlag->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hTimeNoFlag->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hTimeNoFlag->GetListOfFunctions()->At(i);
      			
			TH1F *histoRate = (TH1F*)hTimeNoFlag->Clone("histoRate");
			histoRate->Sumw2();
			if(nEvents != 0) histoRate->Scale(1/nEvents);
			Int_t NbadChannels = 0;
			TString badChannels = "No flag rate too high, ch: ";
			
			for(Int_t i=0; i<16; i++){
				if(histoRate->GetBinContent(i+1)>fMaxNoFlagRate) {
					test = 0.3;
					badChannels += i;
					badChannels += ", ";
					NbadChannels++;
					}
				}
			if(NbadChannels != 0){
				QAbox->Clear();
				QAbox->SetFillColor(kRed); 
				QAbox->AddText(badChannels.Data());
				}
			else{
				QAbox->Clear();
				QAbox->SetFillColor(kGreen);
				QAbox->AddText("No flag rate OK");				
    				}
			}
    		}
    	}
    TH1F *hNEventsBBFlag = (TH1F*)list->At(AliADQADataMakerRec::kNEventsBBFlag);
    if (!hNEventsBBFlag) {
      AliWarning("NEventsBBFlag histogram is not found");
    }
    else {
    	if(hNEventsBBFlag->GetListOfFunctions()->GetEntries()<1) hNEventsBBFlag->GetListOfFunctions()->Add(new TPaveText(0.15,0.67,0.85,0.76,"NDC"));
    	for(Int_t i=0; i<hNEventsBBFlag->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hNEventsBBFlag->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hNEventsBBFlag->GetListOfFunctions()->At(i);
      			
			TH1F *histoRate = (TH1F*)hNEventsBBFlag->Clone("histoRate");
			histoRate->Sumw2();
			if(nEvents != 0) histoRate->Scale(1/nEvents);
			
			Float_t meanRateADA  = 0;
			Float_t meanRateADC  = 0;
			for(Int_t i=1; i<=8; i++){
				meanRateADC += histoRate->GetBinContent(i);
				meanRateADA += histoRate->GetBinContent(i+8);
				}
			meanRateADA = meanRateADA/8;
			meanRateADC = meanRateADC/8;
			Bool_t highVar = kFALSE;
			
			for(Int_t i=1; i<=8; i++){
				if(((TMath::Abs(histoRate->GetBinContent(i)-meanRateADC))>fMaxBBVariation) || 
				   ((TMath::Abs(histoRate->GetBinContent(i+8)-meanRateADA))>fMaxBBVariation)){
					test = 0.7;
					highVar = kTRUE;
					}
				}
			if(highVar){
				QAbox->Clear();
				QAbox->SetFillColor(kYellow); 
				QAbox->AddText("BB rate variation too high");
				}
			else{
				QAbox->Clear();
				QAbox->SetFillColor(kGreen);
				QAbox->AddText("BB rate variation OK");				
    				}
			}
    		}
    	}                                 
    TH1F *hNEventsBGFlag = (TH1F*)list->At(AliADQADataMakerRec::kNEventsBGFlag);
    if (!hNEventsBGFlag) {
      AliWarning("NEventsBGFlag histogram is not found");
    }
    else {
    	if(hNEventsBGFlag->GetListOfFunctions()->GetEntries()<1) hNEventsBGFlag->GetListOfFunctions()->Add(new TPaveText(0.15,0.56,0.85,0.65,"NDC"));
    	for(Int_t i=0; i<hNEventsBGFlag->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hNEventsBGFlag->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hNEventsBGFlag->GetListOfFunctions()->At(i);
      			
			TH1F *histoRate = (TH1F*)hNEventsBGFlag->Clone("histoRate");
			histoRate->Sumw2();
			if(nEvents != 0) histoRate->Scale(1/nEvents);
			
			Float_t meanRateADA  = 0;
			Float_t meanRateADC  = 0;
			for(Int_t i=1; i<=8; i++){
				meanRateADC += histoRate->GetBinContent(i);
				meanRateADA += histoRate->GetBinContent(i+8);
				}
			meanRateADA = meanRateADA/8;
			meanRateADC = meanRateADC/8;
			Bool_t highVar = kFALSE;
			
			for(Int_t i=1; i<=8; i++){
				if(((TMath::Abs(histoRate->GetBinContent(i)-meanRateADC))>fMaxBGVariation) || 
				   ((TMath::Abs(histoRate->GetBinContent(i+8)-meanRateADA))>fMaxBGVariation)){
					test = 0.7;
					highVar = kTRUE;
					}
				}
			if(highVar){
				QAbox->Clear();
				QAbox->SetFillColor(kYellow); 
				QAbox->AddText("BG rate variation too high");
				}
			else{
				QAbox->Clear();
				QAbox->SetFillColor(kGreen);
				QAbox->AddText("BG rate variation OK");				
    				}
			}
    		}
    	}                                                      
    TH2F *hChargeSaturation = (TH2F*)list->At(AliADQADataMakerRec::kChargeSaturation);
    if (!hChargeSaturation) {
      AliWarning("ChargeSaturation histogram is not found");
    }
    else {
    	if(hChargeSaturation->GetListOfFunctions()->GetEntries()<1) hChargeSaturation->GetListOfFunctions()->Add(new TPaveText(0.11,0.40,0.89,0.65,"NDC"));
    	for(Int_t i=0; i<hChargeSaturation->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hChargeSaturation->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hChargeSaturation->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hChargeSlice;
			TString satText = "    ";
			Char_t satValue[5];
			Bool_t	medSat = kFALSE; 
			Bool_t	highSat = kFALSE;
			Bool_t	hugeSat = kFALSE;
			
    			for(Int_t i=0; i<16; i++){
      				hChargeSlice = hChargeSaturation->ProjectionY("hChargeSlice",i+1,i+1);
				Double_t saturation;
				if(hChargeSlice->Integral() != 0) saturation = hChargeSlice->Integral(1000,1025)/hChargeSlice->Integral();
				satText += TString::Format("%1.3f ", saturation);
				satText +=" "; 
				if(saturation > fSatMed && saturation < fSatHigh){
					test = 0.7;
					medSat = kTRUE;
					}
				if(saturation > fSatHigh && saturation < fSatHuge){
					test = 0.3;
					highSat = kTRUE;
					}
				if(saturation > fSatHuge){
					test = 0.1;
					hugeSat = kTRUE;
					}
				}
			if(!medSat && !highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("Saturation");
				QAbox->AddText(satText.Data());
				}
			if(medSat && !highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kYellow);
        			QAbox->AddText("Saturation");
				QAbox->AddText(satText.Data());
				}
			if(highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kOrange);
        			QAbox->AddText("Saturation");
				QAbox->AddText(satText.Data());
				}
			if(hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kRed);
        			QAbox->AddText("Saturation");
				QAbox->AddText(satText.Data());
				}		
    			}
		}
    	}    
    TH2F *hBBFlagVsClock = (TH2F*)list->At(AliADQADataMakerRec::kBBFlagVsClock);
    if (!hBBFlagVsClock) {
      AliWarning("BBFlagVsClock histogram is not found");
    }
    else {
    	if(hBBFlagVsClock->GetListOfFunctions()->GetEntries()<1) hBBFlagVsClock->GetListOfFunctions()->Add(new TPaveText(0.30,0.15,0.70,0.37,"NDC"));
    	for(Int_t i=0; i<hBBFlagVsClock->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hBBFlagVsClock->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hBBFlagVsClock->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hClockSlice;
			Bool_t notConfgADA = kFALSE;
			Bool_t notConfgADC = kFALSE;
			Bool_t notSynchADA = kFALSE;
			Bool_t notSynchADC = kFALSE;

    			for(Int_t i=0; i<16; i++){
      				hClockSlice = hBBFlagVsClock->ProjectionY("hClockSlice",i+1,i+1);
				Double_t center = hClockSlice->GetBinContent(11);
				Double_t flagArray[21];
				for(Int_t iClock = 0; iClock<21; iClock++)flagArray[iClock] = hClockSlice->GetBinContent(iClock+1);
				Int_t maxClock = TMath::LocMax(21,flagArray);
				if(center == 0){
					if(i>7)notConfgADA = kTRUE;
					if(i<8)notConfgADC = kTRUE;
					continue;
					} 
				if(maxClock != 10) {
					if(i>7)notSynchADA = kTRUE;
					if(i<8)notSynchADC = kTRUE;
					}
				}
			if(notConfgADA || notConfgADC){
				QAbox->Clear();
        			QAbox->SetFillColor(kViolet);
        			if(notConfgADA)QAbox->AddText("ADA: dead channels!");
				if(!notConfgADA)QAbox->AddText("ADA ok");
				if(notConfgADC)QAbox->AddText("ADC: dead channels!");
				if(!notConfgADC)QAbox->AddText("ADC ok");
				if(notConfgADA && notConfgADC)QAbox->SetFillColor(kViolet);
				}
				
			else if(notSynchADA || notSynchADC){
				QAbox->Clear();
        			QAbox->SetFillColor(kRed);
        			if(notSynchADA)QAbox->AddText("ADA misconfigured!");
				if(!notSynchADA)QAbox->AddText("ADA ok");
				if(notSynchADC)QAbox->AddText("ADC misconfigured!");
				if(!notSynchADC)QAbox->AddText("ADC ok");
				if(notSynchADA && notSynchADC)QAbox->SetFillColor(kRed);
				}
			else {
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("ADA ok");
				QAbox->AddText("ADC ok");
				}				
    			}
		}
    	} 

    TH2F *hBBFlagVsClock_ADOR = (TH2F*)list->At(AliADQADataMakerRec::kBBFlagVsClock_ADOR);
    if (!hBBFlagVsClock_ADOR) {
      AliWarning("BBFlagVsClock_ADOR histogram is not found");
    }
    else {
    	if(hBBFlagVsClock_ADOR->GetListOfFunctions()->GetEntries()<1) hBBFlagVsClock_ADOR->GetListOfFunctions()->Add(new TPaveText(0.30,0.15,0.70,0.37,"NDC"));
    	for(Int_t i=0; i<hBBFlagVsClock_ADOR->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hBBFlagVsClock_ADOR->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hBBFlagVsClock_ADOR->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hClockSlice;
			Bool_t notConfgADA = kFALSE;
			Bool_t notConfgADC = kFALSE;
			Bool_t notSynchADA = kFALSE;
			Bool_t notSynchADC = kFALSE;

    			for(Int_t i=0; i<16; i++){
      				hClockSlice = hBBFlagVsClock_ADOR->ProjectionY("hClockSlice",i+1,i+1);
				Double_t center = hClockSlice->GetBinContent(11);
				Double_t flagArray[19];
				for(Int_t iClock = 1; iClock<20; iClock++)flagArray[iClock-1] = hClockSlice->GetBinContent(iClock+1);
				Int_t maxClock = TMath::LocMax(19,flagArray);
				if(center == 0){
					if(i>7)notConfgADA = kTRUE;
					if(i<8)notConfgADC = kTRUE;
					continue;
					} 
				if(maxClock != 9) {
					if(i>7)notSynchADA = kTRUE;
					if(i<8)notSynchADC = kTRUE;
					}
				}
			if(notConfgADA || notConfgADC){
				QAbox->Clear();
        			QAbox->SetFillColor(kViolet);
        			if(notConfgADA)QAbox->AddText("ADA: dead channels!");
				if(!notConfgADA)QAbox->AddText("ADA ok");
				if(notConfgADC)QAbox->AddText("ADC: dead channels!");
				if(!notConfgADC)QAbox->AddText("ADC ok");
				if(notConfgADA && notConfgADC)QAbox->SetFillColor(kViolet);
				}
				
			else if(notSynchADA || notSynchADC){
				QAbox->Clear();
        			QAbox->SetFillColor(kRed);
        			if(notSynchADA)QAbox->AddText("ADA misconfigured!");
				if(!notSynchADA)QAbox->AddText("ADA ok");
				if(notSynchADC)QAbox->AddText("ADC misconfigured!");
				if(!notSynchADC)QAbox->AddText("ADC ok");
				if(notSynchADA && notSynchADC)QAbox->SetFillColor(kRed);
				}
			else {
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("ADA ok");
				QAbox->AddText("ADC ok");
				}				
    			}
		}
    	}     
    TH2F *hPedestalDiffInt0  = (TH2F*)list->At(AliADQADataMakerRec::kPedestalDiffInt0);
    if (!hPedestalDiffInt0) {
      AliWarning("PedestalInt0 histogram is not found");
    }
    else {
    	if(hPedestalDiffInt0->GetListOfFunctions()->GetEntries()<1) hPedestalDiffInt0->GetListOfFunctions()->Add(new TPaveText(0.15,0.63,0.85,0.85,"NDC"));
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
				if(TMath::Abs(mean)>fMaxPedDiff) {
					test = 0.3;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kYellow);
						QAbox->AddText("Unstable pedestal for channel ");
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
    	if(hPedestalDiffInt1->GetListOfFunctions()->GetEntries()<1) hPedestalDiffInt1->GetListOfFunctions()->Add(new TPaveText(0.15,0.63,0.85,0.85,"NDC"));
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
				if(TMath::Abs(mean)>fMaxPedDiff) {
					test = 0.3;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kYellow);
						QAbox->AddText("Unstable pedestal for channel ");
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
    if ( !fImage[esIndex] ) {
    	if(AliRecoParam::ConvertIndex(esIndex) != AliRecoParam::kCalib) fImage[esIndex] = new TCanvas(title, title,2500,5500);
	else fImage[esIndex] = new TCanvas(title, title,400,600);
	}
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
    
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(1.1,"xyz");
    gStyle->SetLabelFont(42,"xyz"); 
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetTitleFont(42,"xyz");  
    gStyle->SetTitleOffset(1.0,"xyz");  
    gStyle->SetTitleSize(1.1,"xyz");  
    
    if(AliRecoParam::ConvertIndex(esIndex) != AliRecoParam::kCalib){

    	TVirtualPad* pCharge = 0;
    	TVirtualPad* pTime  = 0;
    	TVirtualPad* pClockFg  = 0;
	TVirtualPad* pClockCh  = 0;
    	TVirtualPad* pCoinc  = 0;
    	TVirtualPad* pPed  = 0;
    	TVirtualPad* pMaxCh  = 0;
	TVirtualPad* pMeanTime  = 0;
    	TVirtualPad* pTrigger  = 0;
	TVirtualPad* pChargeZoom = 0;
	TVirtualPad* pTimeRatio = 0;
	TVirtualPad* pDecision = 0;
	TVirtualPad* pChargeTrend = 0;
	
	const UInt_t nRows = 11;
	const UInt_t nLargeRows = 1;
	Bool_t isLarge[] = {0,0,1,0,0,0,0,0,0,0,0};
	const Double_t timesLarger = 1.5;
	Double_t xRow[nRows];
	xRow[0] = 0.95;
	Double_t xStep = 0.95/((nRows-nLargeRows)+ timesLarger*nLargeRows);
	for(Int_t iRow = 1; iRow<nRows; iRow++){
		if(!isLarge[iRow-1]) xRow[iRow] = xRow[iRow-1]-xStep;
		else xRow[iRow] = xRow[iRow-1]-timesLarger*xStep;
		}
	
	UInt_t iRow = 0;	

    	TPad* pCh = new TPad("Charge", "Charge Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pCh->Draw();
    	pCharge = pCh;
	iRow++;
	
	TPad* pChZ = new TPad("ChargeZoom", "Charge Zoom Pad", 0.0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pChZ->Draw();
    	pChargeZoom = pChZ;
	iRow++;

    	TPad* pT = new TPad("Time", "Time Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pT->Draw();
    	pTime = pT;
	iRow++;
	
	TPad* pTR = new TPad("TimeRatio", "Time Ratio Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pTR->Draw();
    	pTimeRatio = pTR;
	iRow++;
	
	TPad* pMT = new TPad("Mean time", "Mean time Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pMT->Draw();
    	pMeanTime = pMT;
	iRow++;
    
    	TPad* pClFg = new TPad("ClockFg", "ClockFg Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pClFg->Draw();
    	pClockFg = pClFg;
	iRow++;
	
	TPad* pP = new TPad("Pedestal", "Pedestal Pad", 0, xRow[iRow+1], 0.5, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pP->Draw();
    	pPed = pP;
    
    	TPad* pM = new TPad("Max Charge", "Max Charge Pad", 0.5, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pM->Draw();
    	pMaxCh = pM;
	iRow++;
    
    	TPad* pCo = new TPad("Coincidences", "Coincidences Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pCo->Draw();
    	pCoinc = pCo;
	iRow++;
	
        TPad* pTr = new TPad("Triggers", "Triggers Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pTr->Draw();
    	pTrigger = pTr;
	iRow++;

	TPad* pD = new TPad("Decisions", "Decisions Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pD->Draw();
    	pDecision = pD;
	iRow++;
	
	TPad* pCT = new TPad("Charge Trend", "Charge Trend Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pCT->Draw();
    	pChargeTrend = pCT;

    	pCharge->Divide(2, 1);
	pChargeZoom->Divide(6, 1);
    	pTime->Divide(4, 1);
	pTimeRatio->Divide(4, 1);
    	pClockFg->Divide(5, 1);
    	pCoinc->Divide(4, 1);
    	pPed->Divide(2, 1);
	pMeanTime->Divide(4, 1);
	pDecision->Divide(3, 1);
	
	Float_t nEvents = 0;
    	TH1* histo = 0;
	TH1* histoBlue = 0;
	TH1* histoRed = 0;
	TH1* histoNumerator = 0;
	TH1* histoDenominator = 0;
	TH1* histoRatio = 0;
	TH2* histo2D = 0;
	TH2* histo2DCopy = 0;
      	TVirtualPad* pad = 0;
       		
	//Charge pad
      	pad = pCharge->cd(1);
	gPad->SetLogy();
	histoBlue=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeADA);
	nEvents = histoBlue->Integral(-1,-1);
	histoRed=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeADC);
	Float_t max[2];
	max[0] = histoBlue->GetBinContent(histoBlue->GetMaximumBin());
	max[1] = histoRed->GetBinContent(histoRed->GetMaximumBin());
	Float_t min[2];
	min[0] = histoBlue->GetBinContent(histoBlue->GetMinimumBin());
	min[1] = histoRed->GetBinContent(histoRed->GetMinimumBin());
	
	histoBlue->GetYaxis()->SetRangeUser(TMath::MinElement(2,min)+1 ,2*TMath::MaxElement(2,max));
	histoBlue->DrawCopy();
	histoRed->DrawCopy("same");
	
	TLegend *myLegend1 = new TLegend(0.70,0.67,0.97,0.82);
	myLegend1->SetTextFont(42);
  	myLegend1->SetBorderSize(0);
  	myLegend1->SetFillStyle(0);
  	myLegend1->SetFillColor(0);
  	myLegend1->SetMargin(0.25);
  	myLegend1->SetTextSize(0.05);
  	myLegend1->SetEntrySeparation(0.5);
	myLegend1->AddEntry(histoBlue,"ADA","l");
	myLegend1->AddEntry(histoRed,"ADC","l");
	myLegend1->Draw();
	pad = pCharge->cd(2);
	gPad->SetLogz();
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeEoIBB);
	histo->DrawCopy("COLZ");
	//Zoommed charge pad
	pad = pChargeZoom->cd(5);
	gPad->SetLogz();
	histoRatio = (TH1*)histo->Clone("histoRatio");
	histoRatio->GetYaxis()->SetRangeUser(fChargeChannelZoomMin,fChargeChannelZoomMax);
	histoRatio->DrawCopy("COLZ");
	pad = pChargeZoom->cd(4);
	gPad->SetLogz();
	histoDenominator=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeEoI);
	histoDenominator->GetYaxis()->SetRangeUser(fChargeChannelZoomMin,fChargeChannelZoomMax);
	histoDenominator->DrawCopy("COLZ");
	pad = pChargeZoom->cd(6);
	gPad->SetLogz();
	histoRatio->Divide(histoDenominator);
	histoRatio->SetTitle("Ratio: w_BB_Flag/All");
	histoRatio->DrawCopy("COLZ");
	
	for(Int_t iHist = 0; iHist<3; iHist++){
		pad = pChargeZoom->cd(iHist+1);
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeVsClockInt0+iHist);
		if(iHist==2)gPad->SetLogz();
		histo->DrawCopy("COLZ");
		}
	
	//Time pad
	for(Int_t iHist = 0; iHist<4; iHist++){
		pad = pTime->cd(iHist+1);
		if(iHist==3)gPad->SetLogz();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kHPTDCTime+iHist);
		histo->DrawCopy("COLZ");
		}
	//Time ratio pad
	pad = pTimeRatio->cd(1);
	//gPad->SetLogy();
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kFlagNoTime);
	histoBlue = (TH1*)histo->Clone("histoBlue");
	if(nEvents != 0) histoBlue->Scale(1/nEvents);
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kTimeNoFlag);
	histoRed = (TH1*)histo->Clone("histoRed");
	if(nEvents != 0) histoRed->Scale(1/nEvents);
	max[0] = histoBlue->GetBinContent(histoBlue->GetMaximumBin());
	max[1] = histoRed->GetBinContent(histoRed->GetMaximumBin());
	histoBlue->GetYaxis()->SetRangeUser(0,2*TMath::MaxElement(2,max));
	histoBlue->DrawCopy();
	histoRed->DrawCopy("same");
	TLegend *myLegend2 = new TLegend(0.15,0.78,0.85,0.88);
	myLegend2->SetTextFont(42);
  	myLegend2->SetBorderSize(0);
  	myLegend2->SetFillStyle(0);
  	myLegend2->SetFillColor(0);
  	myLegend2->SetMargin(0.25);
  	myLegend2->SetTextSize(0.04);
  	myLegend2->SetEntrySeparation(0.5);
	myLegend2->AddEntry(histoBlue,"Events with BB/BG flag but no time","l");
	myLegend2->AddEntry(histoRed,"Events with time but no BB/BG flag","l");
	myLegend2->Draw();

	histoDenominator=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kHPTDCTimeRebin);
		
	for(Int_t iHist = 1; iHist<3; iHist++){
		pad = pTimeRatio->cd(1+iHist);
		gPad->SetLogz();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kHPTDCTimeRebin + iHist);
		histoRatio = (TH1*)histo->Clone("histoRatio");
		histoRatio->Divide(histoDenominator);
		if(iHist == 1)histoRatio->GetYaxis()->SetRangeUser(fTimeRatioBBZoomMin,fTimeRatioBBZoomMax);
		if(iHist == 2)histoRatio->GetYaxis()->SetRangeUser(fTimeRatioBGZoomMin,fTimeRatioBGZoomMax);
		histoRatio->DrawCopy("COLZ");
		}
	pad = pTimeRatio->cd(4);
	//gPad->SetLogy();
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kNEventsBBFlag);
	histoBlue = (TH1*)histo->Clone("histoBlue");
	histoBlue->Sumw2();
	if(nEvents != 0) histoBlue->Scale(1/nEvents);
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kNEventsBGFlag);
	histoRed = (TH1*)histo->Clone("histoRed");
	histoRed->Sumw2();
	if(nEvents != 0) histoRed->Scale(1/nEvents);
	
	max[0] = histoBlue->GetBinContent(histoBlue->GetMaximumBin());
	max[1] = histoRed->GetBinContent(histoRed->GetMaximumBin());

	min[0] = histoBlue->GetBinContent(histoBlue->GetMinimumBin());
	min[1] = histoRed->GetBinContent(histoRed->GetMinimumBin());
	
	histoBlue->GetYaxis()->SetRangeUser(0.8*TMath::MinElement(2,min),1.8*TMath::MaxElement(2,max));
	histoBlue->DrawCopy("e");
	//histoBlue->DrawCopy("esame");
	//histoRed->DrawCopy("samehist");
	histoRed->DrawCopy("esame");
	TLegend *myLegend3 = new TLegend(0.15,0.78,0.85,0.88);
	myLegend3->SetTextFont(42);
  	myLegend3->SetBorderSize(0);
  	myLegend3->SetFillStyle(0);
  	myLegend3->SetFillColor(0);
  	myLegend3->SetMargin(0.25);
  	myLegend3->SetTextSize(0.04);
  	myLegend3->SetEntrySeparation(0.5);
	myLegend3->AddEntry(histoBlue,"Events with BB flag","l");
	myLegend3->AddEntry(histoRed,"Events with BG flag","l");
	myLegend3->Draw();
		
	//Clock Flags pad
	for(Int_t iHist = 0; iHist<5; iHist++){
		pad = pClockFg->cd(iHist+1);
		if(iHist>1)gPad->SetLogz();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kBBFlagVsClock+iHist);
		histo->DrawCopy("COLZ");
		}
	
	//Coincidences pad
	for(Int_t iHist = 0; iHist<2; iHist++){
		pad = pCoinc->cd(iHist+1);
		gPad->SetLogy();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kNBBCoincADC+2*iHist);
		histo->DrawCopy();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kNBBCoincADA+2*iHist);
		histo->DrawCopy("same");
		myLegend1->Draw();
		}
	for(Int_t iHist = 0; iHist<2; iHist++){
		pad = pCoinc->cd(iHist+3);
		gPad->SetLogz();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kNBBCoincCorr+iHist);
		histo->DrawCopy("COLZ");
		}
	//Pedestal monitoring pad
	for(Int_t iHist = 0; iHist<2; iHist++){
		pad = pPed->cd(iHist+1);
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kPedestalDiffInt0+iHist);
		histo->DrawCopy("COLZ");
		}
	//Saturation monitoring pad
	pad = pMaxCh->cd();
	gPad->SetLogz();
	histo2D=(TH2*)list[esIndex]->At(AliADQADataMakerRec::kChargeSaturation);
	histo2DCopy = (TH2*)histo2D->Clone("histoBlue");
	histo2DCopy->RebinY(64);
	histo2DCopy->DrawCopy("COLZ");

	//Trigger inputs pad
	pad = pTrigger->cd();
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kTriggers);
	gPad->SetLogy();
	histo->DrawCopy("HIST");
	histo->DrawCopy("TEXT0SAME");
	//Mean time pad
	pad = pMeanTime->cd(1);
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kMeanTimeADA);
	histoBlue = (TH1*)histo->Clone("histoBlue");
	
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kMeanTimeADC);
	histoRed = (TH1*)histo->Clone("histoRed");
	max[0] = histoBlue->GetBinContent(histoBlue->GetMaximumBin());
	max[1] = histoRed->GetBinContent(histoRed->GetMaximumBin());
	histoBlue->GetYaxis()->SetRangeUser(0,1.1*TMath::MaxElement(2,max));
	histoBlue->DrawCopy();
	histoRed->DrawCopy("same");
	myLegend1->Draw();
	for(Int_t iHist = 0; iHist<3; iHist++){
		pad = pMeanTime->cd(iHist+2);
		gPad->SetLogz();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kMeanTimeDiff+iHist);
		histo->DrawCopy("COLZ");
		}
	//Decisions pad
	TH1D *hTimeSlice = 0x0;
	TH1F *hMeanTimeVsChargeADA = new TH1F("hMeanTimeVsChargeADA","hMeanTimeVsChargeADA",200,-4,0);
	
	pad = pDecision->cd(1);
	gPad->SetLogz();
	histo2D=(TH2*)list[esIndex]->At(AliADQADataMakerRec::kTimeSlewingADA);
	histo2DCopy = (TH2*)histo2D->Clone("histo2DCopy");
	histo2DCopy->GetYaxis()->SetRangeUser(fTimeRatioBBZoomMin,fTimeRatioBBZoomMax);
	histo2DCopy->SetTitle("Time slewing ADA");
		
	for(Int_t j=0; j<200;j++){
  		hTimeSlice = histo2DCopy->ProjectionY("hTimeSlice",j+1,j+1);
		if(hTimeSlice->GetEntries()<100 && j>100){
			hMeanTimeVsChargeADA->SetBinContent(j+1,hMeanTimeVsChargeADA->GetBinContent(j));
			hMeanTimeVsChargeADA->SetBinError(j+1,hMeanTimeVsChargeADA->GetBinError(j));
			}
		else{
			hMeanTimeVsChargeADA->SetBinContent(j+1,hTimeSlice->GetMean());
        		if(hTimeSlice->GetEntries()>0)hMeanTimeVsChargeADA->SetBinError(j+1,1/TMath::Power(hTimeSlice->GetEntries(),2));
			}
  		}
  	TSpline3 *fTimeSlewingSplineADA = new TSpline3(hMeanTimeVsChargeADA);
	fTimeSlewingSplineADA->SetLineWidth(3);
	fTimeSlewingSplineADA->SetLineColor(kPink);
		
	histo2DCopy->DrawCopy("COLZ");
	fTimeSlewingSplineADA->Draw("Psame");
	delete hMeanTimeVsChargeADA;
	
	TH1F *hMeanTimeVsChargeADC = new TH1F("hMeanTimeVsChargeADC","hMeanTimeVsChargeADC",200,-4,0);
	pad = pDecision->cd(3);
	gPad->SetLogz();
	histo2D=(TH2*)list[esIndex]->At(AliADQADataMakerRec::kTimeSlewingADC);
	histo2DCopy = (TH2*)histo2D->Clone("histo2DCopy");
	histo2DCopy->GetYaxis()->SetRangeUser(fTimeRatioBBZoomMin,fTimeRatioBBZoomMax);
	histo2DCopy->SetTitle("Time slewing ADC");
		
	for(Int_t j=0; j<200;j++){
  		hTimeSlice = histo2DCopy->ProjectionY("hTimeSlice",j+1,j+1);
		if(hTimeSlice->GetEntries()<100 && j>100){
			hMeanTimeVsChargeADC->SetBinContent(j+1,hMeanTimeVsChargeADC->GetBinContent(j));
			hMeanTimeVsChargeADC->SetBinError(j+1,hMeanTimeVsChargeADC->GetBinError(j));
			}
		else{
			hMeanTimeVsChargeADC->SetBinContent(j+1,hTimeSlice->GetMean());
        		if(hTimeSlice->GetEntries()>0)hMeanTimeVsChargeADC->SetBinError(j+1,1/TMath::Power(hTimeSlice->GetEntries(),2));
			}
  		}
  	TSpline3 *fTimeSlewingSplineADC = new TSpline3(hMeanTimeVsChargeADC);
	fTimeSlewingSplineADC->SetLineWidth(3);
	fTimeSlewingSplineADC->SetLineColor(kPink);
		
	histo2DCopy->DrawCopy("COLZ");
	fTimeSlewingSplineADC->Draw("Psame");
	delete hMeanTimeVsChargeADC;
	

	pad = pDecision->cd(2);
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kDecisions);
	histo->Draw("COLZTEXT");
	
	pad = pChargeTrend->cd();
	histoBlue = (TH1*)list[esIndex]->At(AliADQADataMakerRec::kTrend_TriggerChargeQuantileADA);
	histoRed = (TH1*)list[esIndex]->At(AliADQADataMakerRec::kTrend_TriggerChargeQuantileADC);
	histoBlue->GetYaxis()->SetRangeUser(fChargeTrendMin,fChargeTrendMax);
	histoBlue->GetYaxis()->SetTitle("Quantile 0.9");
	histoBlue->GetXaxis()->SetRange(1,histoBlue->GetNbinsX()-1);
	histoBlue->DrawCopy();
	histoRed->DrawCopy("same");
	myLegend1->Draw();
	

    }
    else{
	TVirtualPad* pPed  = 0;
	TH1* histo = 0;
	
    	TPad* pP = new TPad("Pedestal", "Pedestal Pad", 0.0, 0.1, 1.0, 0.95);
    	fImage[esIndex]->cd();
    	pP->Draw();
    	pPed = pP;
	pPed->Divide(1, 2);
	
	TVirtualPad* pad = 0;
	
      	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kPedestalInt0);
      	histo->SetStats(kFALSE);
	pad = pPed->cd(1);
      	histo->DrawCopy("colz");
	
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kPedestalInt1);
      	histo->SetStats(kFALSE); 
	pad = pPed->cd(2);
	histo->DrawCopy("colz");
    }
    
    //fImage[esIndex]->SaveAs(Form("QAsummary_%d_%d.png",AliQAChecker::Instance()->GetRunNumber(),esIndex));
    //fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps"); 
  }
}
 
