/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  P. Cerello Apr 2008
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include <TH1D.h>
#include <TH2.h>
// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliITSQASDDChecker.h"
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSgeomTGeo.h"


ClassImp(AliITSQASDDChecker)
//__________________________________________________________________
AliITSQASDDChecker& AliITSQASDDChecker::operator = (const AliITSQASDDChecker& qac ) 
{
  // Equal operator.
  this->~AliITSQASDDChecker();
  new(this) AliITSQASDDChecker(qac);
  return *this;
}

AliITSQASDDChecker::~AliITSQASDDChecker() 
{

  //destructor
  if(fStepBitSDD) 
    {
      delete[] fStepBitSDD ;
      fStepBitSDD = NULL;
    }
  if(fLowSDDValue)
    {
      delete[]fLowSDDValue;
      fLowSDDValue=NULL;
    }
  if(fHighSDDValue)
    { 
      delete[]fHighSDDValue;
      fHighSDDValue=NULL;
    }
  if(fCalibration)
    {
      delete fCalibration;
      fCalibration=NULL;
    }
} // dtor

//__________________________________________________________________
Double_t AliITSQASDDChecker::Check(AliQAv1::ALITASK_t index, const TObjArray * list, const AliDetectorRecoParam * /*recoparam*/) 
{
  //check histograms of the different lists  
  AliInfo(Form("AliITSQASDDChecker called with offset: %d\n", fSubDetOffset) );

  AliDebug(1,Form("AliITSQASDDChecker called with offset: %d\n", fSubDetOffset));

  Double_t SDDQACheckerValue = 0.;
  TH1 *hdata=NULL;
  Double_t entries=0.;
  Double_t entries2[2];
  for(Int_t i=0;i<2;i++)entries2[i]=0.;

  if(!fCalibration){
    AliCDBEntry *calibSDD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSDD");
    Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
    if(!calibSDD)
      {
	AliError("Calibration object retrieval failed! SDD will not be processed");
	fCalibration = NULL;
	SDDQACheckerValue= fHighSDDValue[AliQAv1::kWARNING];
      }
    fCalibration = (TObjArray *)calibSDD->GetObject();
    
    if(!cacheStatus)calibSDD->SetObject(NULL);
    calibSDD->SetOwner(kTRUE);
    if(!cacheStatus)
      {
	delete calibSDD;
      }
  }

  AliInfo("Calib SDD Created\n ");

  TIter next(list);

  switch(index) {
    case AliQAv1::kRAW:{
      AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index)));
//       if(fRawModulePattern) { delete fRawModulePattern; fRawModulePattern = 0; }
//       if(fRawL3Pattern) { delete fRawL3Pattern; fRawL3Pattern = 0; }
//       if(fRawL4Pattern) { delete fRawL4Pattern; fRawL4Pattern = 0; }
      if (list->GetEntries() == 0){SDDQACheckerValue += fHighSDDValue[AliQAv1::kFATAL];	break;}
      TH1 *hmodule=NULL;
      TH2 *hlayer[2]; 
      Int_t emptymodules[2], filledmodules[2],emptyladders[2],filledladders[2];
      for(Int_t i=0;i<2;i++){emptymodules[i]=0; filledmodules[i]=0; emptyladders[i]=0; filledladders[i]=0; }
      for(Int_t i=0;i<2;i++)hlayer[i]=NULL;   
      while( (hdata = dynamic_cast<TH1* >(next())) ){
	if (hdata){TString hname=hdata->GetName();
	  if(hname.Contains("SDDchargeMap"))continue;
	  if(hname.Contains("SDDModPattern")){
	    if(hname.Contains("NORM")) continue;
	    hmodule=(TH1*)hdata->Clone();
	    entries= hdata->GetEntries();
	    if(AliITSQADataMakerRec::AreEqual(entries,0.)){AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));SDDQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];}//endif entries
	    else{int modmax=hdata->GetNbinsX();
	      Int_t empty=0;
	      Int_t filled=0;
	      Double_t content=0;
	      for(Int_t i=1;i<=modmax;i++){content=hdata->GetBinContent(i);if(AliITSQADataMakerRec::AreEqual(content,0.)) empty++; else filled++; }//end for
	      AliInfo(Form(" %s : empty modules %i \t filled modules %i",hname.Data(), empty, filled));}//end else pattern entries !=0
	  } 		
	  if(hname.Contains("_RelativeOccupancy")) {
	    //fRawModulePattern = (TH1F *) hdata;
	    Float_t threshold = hdata->GetMean() + 4*hdata->GetRMS();
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 1: L3 mean: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 2: L4 mean: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    Int_t aboveThreshold = 0;
	    for(Int_t k=0; k<= hdata->GetNbinsX(); k++) {if(hdata->GetBinLowEdge(k) > threshold) aboveThreshold += (int)(hdata->GetBinContent(k));}
	    Float_t fractionAboveThreshold = ((Float_t) aboveThreshold)/hdata->GetEntries();
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 1, L3: Raw fractionAboveThreshold: %f",fractionAboveThreshold));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 2, L4: Raw fractionAboveThreshold: %f",fractionAboveThreshold));
	    if(fractionAboveThreshold > fThresholdForRelativeOccupancy) {SDDQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	      if(hname.Contains("L3")) AliInfo(Form("SDD check number 1: Set Warning (L3 Raw)"));
	      if(hname.Contains("L4")) AliInfo(Form("SDD check number 2: Set Warning (L4 Raw)")); } }
	  if(hname.Contains("SDDphizL3") || hname.Contains("SDDphizL4")){if(hname.Contains("NORM"))continue;
	    //if(hname.Contains("L3")) {fRawL3Pattern = (TH2F *) hdata;}
	    //if(hname.Contains("L4")) {fRawL4Pattern = (TH2F *) hdata;}
	    Int_t layer=0;
	    if(hname.Contains("3"))layer=0;
	    else  if(hname.Contains("4"))layer=1;
	    entries2[layer]=hdata->GetEntries();
	    if(entries2[layer]==0){AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));
	      SDDQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];}//end if getentries
	    else{
	      Int_t layer1=0;
	      if(hname.Contains("3"))layer1=0;
	      else  if(hname.Contains("4"))layer1=1;
	      TH2* htemp=dynamic_cast<TH2*>(hdata);
	      hlayer[layer1]=(TH2*)htemp->Clone();
	      char newname[50];
	      sprintf(newname,"%s_copy",hname.Data());
	      hlayer[layer1]->SetName(newname);
	      hlayer[layer1]->RebinX(2);
	      int modmay=hlayer[layer1]->GetNbinsY();
	      TH1D* hproj= hlayer[layer1]->ProjectionY();
	      Double_t ladcontent=0;
	      for(Int_t i=1;i<=modmay;i++) {//loop on the ladders
		ladcontent=hproj->GetBinContent(i);
		if(AliITSQADataMakerRec::AreEqual(ladcontent,0.)) emptyladders[layer1]++;
		else filledladders[layer1]++;}//end for
	      AliInfo(Form(" %s : empty ladders %i \t filled ladders %i\n",hname.Data(), emptyladders[layer], filledladders[layer]));//end else layer 3
	      delete hproj;
	      hproj=NULL;}//end else entries !=0
	  }//end check on phiz	      
	}//end if hdata	
      }//end while
      if(AliITSQADataMakerRec::AreEqual(entries,0.)&&AliITSQADataMakerRec::AreEqual(entries2[0],0.)&&AliITSQADataMakerRec::AreEqual(entries2[1],0.)) break;
      //else{
      if(hmodule || (hlayer[0] && hlayer[1])){
	Int_t excluded=0;
	Int_t active=0;
	Int_t exactive=0;//excluded but taking data
	for(Int_t imod=0;imod<fgknSDDmodules;imod++){
	  Int_t lay=0;
	  Int_t lad=0;
	  Int_t det=0;
	  Int_t module=0;
	  module=imod+fgkmodoffset;
	  AliITSCalibrationSDD * cal=(AliITSCalibrationSDD*)fCalibration->At(imod);
	  if(cal==0) { delete cal; continue;}
	  AliITSgeomTGeo::GetModuleId(module,lay,lad,det);
	  if (cal->IsBad()){
	    excluded++;
	    Double_t content=0.;
	    Double_t contentlayer[2];
	    for(Int_t i=0;i<2;i++)contentlayer[i]=0.;
	    if(hmodule)content=hmodule->GetBinContent(imod+1);//if expert bit is active the histogram has been created 
	    contentlayer[lay-3]=hlayer[lay-3]->GetBinContent(det,lad);
	    if(AliITSQADataMakerRec::AreEqual(content,0.)== kFALSE || AliITSQADataMakerRec::AreEqual(contentlayer[lay-3],0.)==kFALSE) {
	      filledmodules[lay-3]++;
	      AliWarning(Form("The module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data \n ",module,lay,lad,det));
	      exactive++;
	    } else if(AliITSQADataMakerRec::AreEqual(content,0.) && AliITSQADataMakerRec::AreEqual(contentlayer[lay-3],0.)) 
	      emptymodules[lay-3]++;
	  } else {
	    Double_t contentgood=0.;
	    active++;
	    contentgood=hlayer[lay-3]->GetBinContent(det,lad);
	    if(AliITSQADataMakerRec::AreEqual(contentgood,0.)) 
	      emptymodules[lay-3]++;
	    else 
	      filledmodules[lay-3]++;
	  }
	}//end for
	for(Int_t i=0;i<2;i++){AliInfo(Form("Layer %i \tempty modules %i \t filled modules %i\n", i+3,emptymodules[i], filledmodules[i]));}//end else layers
	if(exactive==0){
	  AliInfo(Form("All the active modules (%i) are in acquisition. The number of excluded modules are %i \n",active,excluded));
	  SDDQACheckerValue=fHighSDDValue[AliQAv1::kINFO];
	}
	if(exactive!=0){
	  AliWarning(Form("%i modules excluded from the acquisition took data. Active modules%i \n ",exactive,active));
	  SDDQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	}
	if(excluded==exactive){
	  AliWarning(Form("All the modules excluded from the acquisition (%d) took data!  Active modules %i\n",excluded,active));
	  SDDQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	}
	if(active==0){
	  AliWarning(Form("No modules took data: excluded %i \t exactive %i \n", excluded, exactive)); 
	  SDDQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	}
	
      }//end else 
      delete hmodule;
      hmodule=NULL;
      for(Int_t i=0;i<2;i++) {
	delete hlayer[i];
	hlayer[i]=NULL;
      }

    }
      
      break;
      
  case AliQAv1::kNULLTASK:{
    AliInfo(Form("No Check on %s\n",AliQAv1::GetAliTaskName(index))); 
    SDDQACheckerValue=1.;
  }
    break;
    
  case AliQAv1::kREC:
    {
      
      AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index))); 

      if (list->GetEntries() == 0){ //check if the list is empty
	//printf("SDDQACheckerValue = %f \t value %f\n",SDDQACheckerValue,fHighSDDValue[AliQAv1::kFATAL]);
	SDDQACheckerValue=fHighSDDValue[AliQAv1::kFATAL]; 
	break;			
      }//end if getentries
      
      while((hdata=dynamic_cast<TH1* >(next()))){
	if (hdata){
	  TString hname=hdata->GetName();
	  if(hname.Contains("_RelativeOccupancy")) {
	    Float_t threshold = hdata->GetMean() + 4*hdata->GetRMS();
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 3: L3 mean: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 4: L4 mean: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    Int_t aboveThreshold = 0;
	    for(Int_t k=0; k<= ((Int_t)hdata->GetNbinsX()); k++) {
	      if(hdata->GetBinLowEdge(k) > threshold) aboveThreshold += (Int_t)(hdata->GetBinContent(k));
	    }
	    Float_t fractionAboveThreshold = ((Float_t) aboveThreshold)/hdata->GetEntries();
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 3, L3: RecPoints fractionAboveThreshold: %f",fractionAboveThreshold));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 4, L4: RecPoints fractionAboveThreshold: %f",fractionAboveThreshold));
	    if(fractionAboveThreshold > fThresholdForRelativeOccupancy) { 
	      SDDQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	      if(hname.Contains("L3")) AliInfo(Form("SDD check number 3: Set Warning (L3 RecPoints)"));
	      if(hname.Contains("L4")) AliInfo(Form("SDD check number 4: Set Warning (L4 RecPoints)"));
	    }
	  }
	  if(hname.Contains("Rec2Raw") && !hname.Contains("2D")) {
	    //Float_t threshold = 0.;
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 5: L3 R2R mean: %f, rms: ,%f",((TH1F *) hdata)->GetMean(),((TH1F *) hdata)->GetRMS()));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 6: L4 R2R mean: %f, rms: ,%f",((TH1F *) hdata)->GetMean(),((TH1F *) hdata)->GetRMS()));
	    Int_t belowThreshold = 0;
	    for(Int_t k=0; k<=((TH1F *)hdata)->GetNbinsX(); k++) {
	      if(((TH1F *) hdata)->GetBinLowEdge(k) < fThresholdForRecToRawRatio) belowThreshold += ((Int_t)((TH1F *) hdata)->GetBinContent(k));
	    }
	    Double_t fractionBelowThreshold =0.;
	    Double_t entries3=((TH1F *)hdata)->GetEntries();
	    if(entries3>0.001)fractionBelowThreshold = ((Double_t)(belowThreshold))/entries3;
	    else{ AliWarning(Form("No entries on %s. The check will retuns zero.\n",hdata->GetName() )); }
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 5, L3: RecPoints2Raws fractionBelowThreshold: %f",fractionBelowThreshold));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 6, L4: RecPoints2Raws fractionBelowThreshold: %f",fractionBelowThreshold));
	    if(fractionBelowThreshold > fThresholdForRelativeOccupancy) { 
	      SDDQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	      if(hname.Contains("L3")) AliInfo(Form("SDD check number 5: Set Warning (L3 RecPoints2Raws)"));
	      if(hname.Contains("L4")) AliInfo(Form("SDD check number 6: Set Warning (L4 RecPoints2Raws)"));
	    }
	  }
	  if(hname.Contains("dedx")) {
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 7: L3 average charge: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 8: L4 average charge: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	  }
	}
      }				
      
      SDDQACheckerValue=1.;
    }
    break;
  case AliQAv1::kANA:
    {
      AliInfo(Form("===================> No Check on %s\n",AliQAv1::GetAliTaskName(index)));
      SDDQACheckerValue=1.; 
    }
    break;
  case AliQAv1::kESD:
    {
      AliInfo(Form("==================>  No Check on %s\n",AliQAv1::GetAliTaskName(index)));
      SDDQACheckerValue=1.;
    } 
    break;
  case AliQAv1::kNTASK:{
    AliInfo(Form("==================>  No Check on %s\n",AliQAv1::GetAliTaskName(index))); 
    SDDQACheckerValue=1.;
  }
    break;
  case AliQAv1::kSIM:{
    AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index))); 
    Int_t uid=list->GetUniqueID();
    if(uid==60) {
      //digits
      if (list->GetEntries() == 0){ 
	SDDQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	break;
      } else{
	
	while( (hdata = dynamic_cast<TH1* >(next())) ){
	  if (hdata){
	    if(hdata->GetEntries()==0)SDDQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];
	    else {
	      TString hname=hdata->GetName();
	      if(hname.Contains("SDDDIGITSModulePattern")) {
		//see raws
		
		SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDDAnodeDistribution")) {
		SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDDTbinDistribution")) {
		//to do as rp
		SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDDADCCountsDistribution")) {
		SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      }//end adc counts
	      
	    }//end entries !=0
	  }//end hdata
	}//end while
      }//end else
    } else if(uid==50) 
      {
	//hits
	if (list->GetEntries() == 0){ 
	  SDDQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	  break;
	} 
	else{
	  
	  while( (hdata = dynamic_cast<TH1* >(next())) ){
	    if (hdata){
	      if(hdata->GetEntries()==0)SDDQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];
	      else {
		TString hname=hdata->GetName();
		if(hname.Contains("SDDHITSModulePattern")) {
		  //to do as raws
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDHITlenghtalonglocalYCoord")) {
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDHITlenghtalonglocalYCoordZoom")) {
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDDepositedEnergyDistribution")) {
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		}//end deposited energy
		
	      }//end entries !=0
	    }//end hdata
	  }//end while
	}//end else
      } else if(uid==70) 
      {
	//sdigits
	if (list->GetEntries() == 0){ 
	  SDDQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	  break;
	} else{
	  
	  while( (hdata = dynamic_cast<TH1* >(next())) ){
	    if (hdata){
	      if(hdata->GetEntries()==0)SDDQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];
	      else {
		TString hname=hdata->GetName();
		if(hname.Contains("SDDSDIGITSModulePattern")) {
		  //to do as raws
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDAnodeDistribution")) {
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDTbinDistribution")) {
		  //to do as rp
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDADCCountsDistribution")) {
		  SDDQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		}//end adc counts bindistribution
	      }//end entries !=0
	    }//end hdata
	  }//end while
	}//end else
      }//end sdigits
    SDDQACheckerValue=1.;
  }
    break;
    
  }//end switch
  
  fCalibration=NULL;
  delete hdata;


  return SDDQACheckerValue;	
}

//__________________________________________________________________
void AliITSQASDDChecker::SetTaskOffset(Int_t taskoffset)
{
  //set the number of the histograms already present in the list before the SDD histograms
  fSubDetOffset = taskoffset;
}


//__________________________________________________________________
void AliITSQASDDChecker::SetStepBit(const Double_t *steprange)
{
  //set the values of the step bit for each QA bit range calculated in the AliITSQAChecker class
  fStepBitSDD = new Double_t[AliQAv1::kNBIT];
  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fStepBitSDD[bit]=steprange[bit];
    }
}

//__________________________________________________________________
void  AliITSQASDDChecker::SetSDDLimits(const Float_t *lowvalue, const Float_t * highvalue)
{
  //set the low and high values in for each QA bit range calculated in the AliITSQAChecker class
  fLowSDDValue = new Float_t[AliQAv1::kNBIT];
  fHighSDDValue= new Float_t[AliQAv1::kNBIT];

  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fLowSDDValue[bit]=lowvalue[bit];
      fHighSDDValue[bit]= highvalue[bit];
    }

}

