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
    	if(hChargeEoIInt0->GetListOfFunctions()->GetEntries()<1) hChargeEoIInt0->GetListOfFunctions()->Add(new TPaveText(0.4,0.53,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hChargeEoIInt0->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hChargeEoIInt0->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hChargeEoIInt0->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hChargeSlice;
			
			TString messageText = " "; 
			TString satTextADA = " A-side = ";
			TString satTextADC = " C-side = ";
			Char_t satValue[5];
			Bool_t	medSat = kFALSE; 
			Bool_t	highSat = kFALSE;
			Bool_t	hugeSat = kFALSE;
			
    			for(Int_t i=0; i<2; i++){
      				hChargeSlice = hChargeEoIInt0->ProjectionY("hChargeSlice",8*i+1,8*i+9);
				Double_t saturation = hChargeSlice->Integral(1000,1025)/hChargeSlice->Integral();
				sprintf(satValue, "%1.3f", saturation);
				if(i == 0) satTextADC +=satValue; 
				if(i == 1) satTextADA +=satValue;
				if(saturation > 0.1 && saturation < 0.3){
					test = 0.7;
					medSat = kTRUE;
					}
				if(saturation > 0.3 && saturation < 0.5){
					test = 0.3;
					highSat = kTRUE;
					}
				if(saturation > 0.5){
					test = 0.1;
					hugeSat = kTRUE;
					}
				}
			if(!medSat && !highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("Saturation OK");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
				}
			if(medSat && !highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kYellow);
        			QAbox->AddText("Medium Saturation");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
				}
			if(highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kOrange);
        			QAbox->AddText("High Saturation");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
				}
			if(hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kRed);
        			QAbox->AddText("Very High Saturation");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
				}		
							
    			}
		}
    	}
    TH2F *hChargeEoIInt1 = (TH2F*)list->At(AliADQADataMakerRec::kChargeEoIInt1);
    if (!hChargeEoIInt1) {
      AliWarning("ChargeEoIInt1 histogram is not found");
    }
    else {
    	if(hChargeEoIInt1->GetListOfFunctions()->GetEntries()<1) hChargeEoIInt1->GetListOfFunctions()->Add(new TPaveText(0.4,0.53,0.85,0.85,"NDC"));
    	for(Int_t i=0; i<hChargeEoIInt1->GetListOfFunctions()->GetEntries(); i++){
     		TString funcName = hChargeEoIInt1->GetListOfFunctions()->At(i)->ClassName();
     		if(funcName.Contains("TPaveText")){
      			TPaveText *QAbox = (TPaveText*)hChargeEoIInt1->GetListOfFunctions()->At(i);
      			QAbox->Clear();

    			TH1D *hChargeSlice;
			
			TString messageText = " "; 
			TString satTextADA = " A-side = ";
			TString satTextADC = " C-side = ";
			Char_t satValue[5];
			Bool_t	medSat = kFALSE; 
			Bool_t	highSat = kFALSE;
			Bool_t	hugeSat = kFALSE;
			
    			for(Int_t i=0; i<2; i++){
      				hChargeSlice = hChargeEoIInt1->ProjectionY("hChargeSlice",8*i+1,8*i+9);
				Double_t saturation = hChargeSlice->Integral(1000,1025)/hChargeSlice->Integral();
				sprintf(satValue, "%1.3f", saturation);
				if(i == 0) satTextADC +=satValue; 
				if(i == 1) satTextADA +=satValue;
				if(saturation > 0.1 && saturation < 0.3){
					test = 0.7;
					medSat = kTRUE;
					}
				if(saturation > 0.3 && saturation < 0.5){
					test = 0.3;
					highSat = kTRUE;
					}
				if(saturation > 0.5){
					test = 0.1;
					hugeSat = kTRUE;
					}
				}
			if(!medSat && !highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kGreen);
        			QAbox->AddText("Saturation OK");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
				}
			if(medSat && !highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kYellow);
        			QAbox->AddText("Medium Saturation");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
				}
			if(highSat && !hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kOrange);
        			QAbox->AddText("High Saturation");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
				}
			if(hugeSat){
				QAbox->Clear();
        			QAbox->SetFillColor(kRed);
        			QAbox->AddText("Very High Saturation");
				QAbox->AddText(satTextADA.Data());
				QAbox->AddText(satTextADC.Data());
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
				if(TMath::Abs(mean)>1) {
					test = 0.3;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kOrange);
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
				if(TMath::Abs(mean)>1) {
					test = 0.3;
					if(NbadChannels == 0){
						QAbox->SetFillColor(kOrange);
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
    	if(esIndex != AliRecoParam::kCalib) fImage[esIndex] = new TCanvas(title, title,2000,3500);
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
    
    if(esIndex != AliRecoParam::kCalib){

    	TVirtualPad* pCharge = 0;
    	TVirtualPad* pTime  = 0;
    	TVirtualPad* pClock  = 0;
    	TVirtualPad* pCoinc  = 0;
    	TVirtualPad* pPed  = 0;
    	TVirtualPad* pMaxCh  = 0;
	TVirtualPad* pMeanTime  = 0;
    	TVirtualPad* pTrigger  = 0;
	
	const UInt_t nRows = 7;
	const UInt_t nLargeRows = 1;
	Bool_t isLarge[] = {0,1,0,0,0,0,0};
	const Double_t timesLarger = 2.0;
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

    	TPad* pT = new TPad("Time", "Time Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pT->Draw();
    	pTime = pT;
	iRow++;
	
	TPad* pMT = new TPad("Mean time", "Mean time Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pMT->Draw();
    	pMeanTime = pMT;
	iRow++;
    
    	TPad* pCl = new TPad("Clock", "Clock Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pCl->Draw();
    	pClock = pCl;
	iRow++;
    
    	TPad* pCo = new TPad("Coincidences", "Coincidences Pad", 0, xRow[iRow+1], 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pCo->Draw();
    	pCoinc = pCo;
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
	
	TPad* pTr = new TPad("Triggers", "Triggers Pad", 0, 0.0, 1.0, xRow[iRow]);
    	fImage[esIndex]->cd();
    	pTr->Draw();
    	pTrigger = pTr;

    	pCharge->Divide(2, 1);
    	pTime->Divide(4, 1);
    	pClock->Divide(4, 1);
    	pCoinc->Divide(4, 1);
    	pPed->Divide(2, 1);
    	pMaxCh->Divide(2, 1);
	pMeanTime->Divide(4, 1);
	
    	TH1* histo = 0;
	TH1* histoADA = 0;
	TH1* histoADC = 0;
      	TVirtualPad* pad = 0;
       		
	//Charge pad
      	pad = pCharge->cd(1);
	gPad->SetLogy();
	histoADA=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeADA);
	histoADA->DrawCopy();
	histoADC=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeADC);
	histoADC->DrawCopy("same");
	TLegend *myLegend1 = new TLegend(0.70,0.67,0.97,0.82);
	myLegend1->SetTextFont(42);
  	myLegend1->SetBorderSize(0);
  	myLegend1->SetFillStyle(0);
  	myLegend1->SetFillColor(0);
  	myLegend1->SetMargin(0.25);
  	myLegend1->SetTextSize(0.04);
  	myLegend1->SetEntrySeparation(0.5);
	myLegend1->AddEntry(histoADA,"ADA","l");
	myLegend1->AddEntry(histoADC,"ADC","l");
	myLegend1->Draw();
	pad = pCharge->cd(2);
	gPad->SetLogz();
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeEoI);
	histo->DrawCopy("COLZ");
	//Time pad
	for(Int_t iHist = 0; iHist<4; iHist++){
		pad = pTime->cd(iHist+1);
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kHPTDCTime+iHist);
		histo->DrawCopy("COLZ");
		}
	//Clock pad
	for(Int_t iHist = 0; iHist<4; iHist++){
		pad = pClock->cd(iHist+1);
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeVsClockInt0+iHist);
		histo->DrawCopy("COLZ");
		}
	//Coincidences pad
	for(Int_t iHist = 0; iHist<2; iHist++){
		pad = pCoinc->cd(iHist+1);
		gPad->SetLogy();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kNBBCoincADA+2*iHist);
		histo->DrawCopy();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kNBBCoincADC+2*iHist);
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
	for(Int_t iHist = 0; iHist<2; iHist++){
		pad = pMaxCh->cd(iHist+1);
		gPad->SetLogz();
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kChargeEoIInt0+iHist);
		histo->DrawCopy("COLZ");
		}
	//Trigger inputs pad
	pad = pTrigger->cd();
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kTriggers);
	gPad->SetLogy();
	histo->DrawCopy();
	//Mean time pad
	pad = pMeanTime->cd(1);
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kMeanTimeADA);
	histo->DrawCopy();
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kMeanTimeADC);
	histo->DrawCopy("same");
	myLegend1->Draw();
	for(Int_t iHist = 0; iHist<3; iHist++){
		pad = pMeanTime->cd(iHist+2);
		histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kMeanTimeDiff+iHist);
		histo->DrawCopy("COLZ");
		}	

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
	
	histo=(TH1*)list[esIndex]->At(AliADQADataMakerRec::kPedestalInt0);
      	histo->SetStats(kFALSE); 
	pad = pPed->cd(2);
	histo->DrawCopy("colz");
    }
    
    fImage[esIndex]->SaveAs(Form("QAcanvas%d.png",esIndex));
    fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps"); 
  }
}
 
