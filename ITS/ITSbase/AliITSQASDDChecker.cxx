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
//  first implementation: P. Cerello Apr 2008
//  last review: F. Prino Apr 2015
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include <TH1D.h>
#include <TH2.h>
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
//#include "TPaletteAxis.h"
// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliITSQASDDChecker.h"
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSgeomTGeo.h"
#include "AliQAManager.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliQACheckerBase.h"


ClassImp(AliITSQASDDChecker)



//_____________________________________________________________________

AliITSQASDDChecker::AliITSQASDDChecker():TObject(),
  fSubDetOffset(0),
  fStepBitSDD(NULL),
  fLowSDDValue(NULL),
  fHighSDDValue(NULL),
  fCalibration(NULL),
  fThresholdForRelativeOccupancy(0.01),
  fThresholdForRecToRawRatio(0.04),
  fImage(NULL),
  fESforCheck(0)
{
// Default constructor
  fStepBitSDD=new Double_t[AliQAv1::kNBIT];
  fLowSDDValue=new Float_t[AliQAv1::kNBIT];
  fHighSDDValue=new Float_t[AliQAv1::kNBIT];
  
  for(Int_t ibit=0;ibit<AliQAv1::kNBIT;ibit++){
    fStepBitSDD[ibit]=0.;
    fLowSDDValue[ibit]=0.;
    fHighSDDValue[ibit]=0.;
  }
  for(Int_t i=0;i<AliRecoParam::kNSpecies;i++) fPaveText[i] = NULL;
}          // ctor



AliITSQASDDChecker::~AliITSQASDDChecker() 
{
  //destructor
  delete[] fStepBitSDD ;
  delete[]fLowSDDValue;
  delete[]fHighSDDValue;
  delete fCalibration;
  delete []fImage; 

  for(Int_t i=0;i<AliRecoParam::kNSpecies;i++) {
    delete fPaveText[i]; 
  }
} // dtor

//__________________________________________________________________
Double_t AliITSQASDDChecker::Check(AliQAv1::ALITASK_t index, const TObjArray * list, const AliDetectorRecoParam * /*recoparam*/) 
{
  //check histograms of the different lists  
  AliInfo(Form("AliITSQASDDChecker called with offset: %d \t and specie %d", fSubDetOffset,fESforCheck));

  AliDebug(1,Form("AliITSQASDDChecker called with offset: %d", fSubDetOffset));

  Double_t sddQACheckerValue = 0.;
  TH1 *hdata=NULL;

  if(!fCalibration){
    AliCDBEntry *calibSDD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSDD");
    Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
    if(!calibSDD)
      {
	AliError("Calibration object retrieval failed! SDD will not be processed");
	fCalibration = NULL;
	sddQACheckerValue= fHighSDDValue[AliQAv1::kWARNING];
      }
    else{
      fCalibration = (TObjArray *)calibSDD->GetObject();
      
      if(!cacheStatus)calibSDD->SetObject(NULL);
      calibSDD->SetOwner(kTRUE);
      if(!cacheStatus)
	{
	  delete calibSDD;
	}
    }//end calibsdd 
  }//end fcalibration

  AliInfo("Calib SDD Created");

  TIter next(list);

  switch(index) {
  case AliQAv1::kRAW:{

    AliInfo(Form("Check on %s",AliQAv1::GetAliTaskName(index)));
    if (list->GetEntries() == 0){
      AliError("Raw List for SDD is empty");
      sddQACheckerValue += fHighSDDValue[AliQAv1::kFATAL];	
      break;
    }
    TH1 *hmodule=NULL;
    TH2 *hlayer[2]; 
    hdata=NULL;
    for(Int_t i=0;i<2;i++)hlayer[i]=NULL;

    //check counters
    Int_t counters[kNumOfSDDCheckerCounters][3];
    for(Int_t ic=0; ic<kNumOfSDDCheckerCounters; ic++){
      for(Int_t jl=0; jl<3; jl++){
	counters[ic][jl]=0;
      }
    }

    Int_t neventsraw=0;

    //take the number of events

    while( (hdata = dynamic_cast<TH1* >(next())) ){
      if (hdata){
	TString hname=hdata->GetName();	
	if(hname.Contains("SDDRawDataCheck")){ 
	  neventsraw=(Int_t)hdata->GetBinContent(1); 
	}
	if(hname.Contains("SDDModPattern")){
	  if(hname.Contains("NORM")) continue;
	  hmodule=(TH1*)hdata->Clone("hmoduleOcc_temp");
	  if(hdata->GetEntries()<0.1){
	    AliWarning(Form("===================>>>>>> No entries in  %s",hname.Data()));
	  }
	}
	if(hname.Contains("SDDphizL3") || hname.Contains("SDDphizL4")){
	  if(hname.Contains("NORM"))continue;
	  Int_t layer=0;
	  if(hname.Contains("3"))layer=0;
	  else  if(hname.Contains("4"))layer=1;
	  hlayer[layer]=(TH2*)hdata->Clone(Form("%s_copy",hname.Data()));
	  if(hdata->GetEntries()<0.1){
	    AliWarning(Form("===================>>>>>> No entries in  %s",hname.Data()));
	  }
	}//end check on phiz	      	
      }//end if hdata	
    }//end while

    if(hmodule || (hlayer[0] && hlayer[1])){
      FillCounters(counters,hmodule,hlayer[0],hlayer[1]);

      for(Int_t i=0;i<2;i++){
	AliInfo(Form("Layer %i   \tempty modules %i             \t filled modules %i", i+3,counters[kEmptyMod][i], counters[kFilledMod][i]));
	AliInfo(Form("Layer %i   \tempty single drift regions %i \t filled single drift regions %i",i+3,counters[kEmptyWing][i], counters[kFilledWing][i]));
      }
      
      for(Int_t ic=0; ic<kNumOfSDDCheckerCounters; ic++){
	counters[ic][2]=counters[ic][0]+counters[ic][1];
      }

      AliInfo(Form("In total %d modules and %d single drift regions took data. ",counters[kFilledMod][2], counters[kFilledWing][2]));
      AliInfo(Form("In total %d modules and %d single drift regions were empty",counters[kEmptyMod][2],counters[kEmptyWing][2]));
	    
      next.Begin();
      Int_t outputLay3=0,outputLay4=0,outputAll=0;
      while( (hdata=dynamic_cast<TH1* >(next())) ) {
	if (hdata){
	  TString hname=hdata->GetName();
	  if(hname.Contains("NORM"))continue;
	  if(hname.Contains("SDDModPattern")){
	    TPaveText *ptext = ((TPaveText *)hdata->GetListOfFunctions()->FindObject("TPave"));
	    outputAll=CheckCounters(counters,2,neventsraw,ptext);
	  }
	  if(hname.Contains("SDDphizL3")){
	    TPaveText *ptext = ((TPaveText *)hdata->GetListOfFunctions()->FindObject("TPave"));
	    outputLay3=CheckCounters(counters,0,neventsraw,ptext);
	  }else if(hname.Contains("SDDphizL4")){
	    TPaveText *ptext = ((TPaveText *)hdata->GetListOfFunctions()->FindObject("TPave"));
	    outputLay4=CheckCounters(counters,1,neventsraw,ptext);
	  } else if(hname.Contains("SDDRawDataCheck")) {
	    for(Int_t ic=0; ic<kNumOfSDDCheckerCounters; ic++){
	      ((TH1F*)hdata)->SetBinContent(5+ic,counters[ic][2]);
	      //layer 3
	      ((TH1F*)hdata)->SetBinContent(19+ic,counters[ic][0]);
	      //layer 4
	      ((TH1F*)hdata)->SetBinContent(33+ic,counters[ic][1]);
	    }
	  }
	}	
      }
      if(outputLay3==fHighSDDValue[AliQAv1::kERROR] || outputLay4==fHighSDDValue[AliQAv1::kERROR]) sddQACheckerValue=fHighSDDValue[AliQAv1::kERROR];
      else if(outputLay3==fHighSDDValue[AliQAv1::kWARNING] || outputLay4==fHighSDDValue[AliQAv1::kWARNING]) sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
      else sddQACheckerValue=fHighSDDValue[AliQAv1::kINFO];
    }

    delete hmodule;
    hmodule=NULL;
    for(Int_t i=0;i<2;i++) {
      delete hlayer[i];
      hlayer[i]=NULL;
    }
    
  }//end raw
    
    break;
    
  case AliQAv1::kNULLTASK:{
    AliInfo(Form("No Check on %s",AliQAv1::GetAliTaskName(index))); 
    sddQACheckerValue=1.;
  }
    break;
    
  case AliQAv1::kREC:{

    Int_t uidrec=list->GetUniqueID();
    AliInfo(Form("Check on %s",AliQAv1::GetAliTaskName(index))); 
    if(uidrec==20){
      //recpoints
      if (list->GetEntries() == 0){ //check if the list is empty
	//printf("sddQACheckerValue = %f \t value %f\n",sddQACheckerValue,fHighSDDValue[AliQAv1::kFATAL]);
	sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL]; 
	//break;			
      }//end if getentries
      
      
      TH1 *hmodule=NULL;
      TH2 *hlayer[2]; 
      hdata=NULL;
      for(Int_t i=0;i<2;i++)hlayer[i]=NULL;
      
      //check counters
      Int_t counters[kNumOfSDDCheckerCounters][3];
      for(Int_t ic=0; ic<kNumOfSDDCheckerCounters; ic++){
	for(Int_t jl=0; jl<3; jl++){
	  counters[ic][jl]=0;
	}
      }
      
      TString results1;
      TString results2;
      Int_t neventsrecpoints=0;

      while( (hdata = dynamic_cast<TH1* >(next())) ){
	if (hdata){
	  TString hname=hdata->GetName();	  
	  if(hname.Contains("SDDRecPointCheck")){ 
	    neventsrecpoints=(Int_t)hdata->GetBinContent(1); 
	  }
	  if(hname.Contains("SDDModPatternRP")){
	    if(hname.Contains("NORM")) continue;
	    hmodule=(TH1*)hdata->Clone();
	    if(hdata->GetEntries()<0.1){
	      AliWarning(Form("===================>>>>>> No entries in  %s ",hname.Data()));
	    }
	  }
	  if(hname.Contains("SDDModPatternL3RP") || hname.Contains("SDDModPatternL4RP")){
	    if(hname.Contains("NORM"))continue;
	    Int_t layer=0;
	    if(hname.Contains("3"))layer=0;
	    else  if(hname.Contains("4"))layer=1;
	    hlayer[layer]=(TH2*)hdata->Clone(Form("%s_copy",hname.Data()));
	    if(hdata->GetEntries()<0.1){
	      AliWarning(Form("===================>>>>>> No entries in  %s ",hname.Data()));
	    }
	  }//end check on phiz
	}//end if hdata
      }//end while
      
      if(hmodule || (hlayer[0] && hlayer[1])){
	FillCounters(counters,hmodule,hlayer[0],hlayer[1]);

	for(Int_t i=0;i<2;i++){
	  AliInfo(Form("Layer %i   \tempty modules %i             \t filled modules %i", i+3,counters[kEmptyMod][i], counters[kFilledMod][i]));
	  AliInfo(Form("Layer %i   \tempty single drift regions %i \t filled single drift regions %i",i+3,counters[kEmptyWing][i], counters[kFilledWing][i]));
	}
      
	for(Int_t ic=0; ic<kNumOfSDDCheckerCounters; ic++){
	  counters[ic][2]=counters[ic][0]+counters[ic][1];
	}

	AliInfo(Form("In total %d modules and %d single drift regions took data. ",counters[kFilledMod][2], counters[kFilledWing][2]));
	AliInfo(Form("In total %d modules and %d single drift regions were empty",counters[kEmptyMod][2],counters[kEmptyWing][2]));


	next.Begin();
	Int_t outputLay3=0,outputLay4=0;;
	while( (hdata=dynamic_cast<TH1* >(next())) ) {
	  if (hdata){
	    TString hname=hdata->GetName();
	    if(hname.Contains("NORM"))continue;
	    if(hname.Contains("SDDModPatternL3RP")){
	      TPaveText *ptext = ((TPaveText *)hdata->GetListOfFunctions()->FindObject("TPave"));
	      outputLay3=CheckCounters(counters,0,neventsrecpoints,ptext);
	    }else if(hname.Contains("SDDModPatternL4RP")){
	      TPaveText *ptext = ((TPaveText *)hdata->GetListOfFunctions()->FindObject("TPave"));
	      outputLay4=CheckCounters(counters,1,neventsrecpoints,ptext);
	    }else if(hname.Contains("SDDRecPointCheck")) {
	      for(Int_t ic=0; ic<kNumOfSDDCheckerCounters; ic++){
		((TH1F*)hdata)->SetBinContent(5+ic,counters[ic][2]);
		//layer 3
		((TH1F*)hdata)->SetBinContent(19+ic,counters[ic][0]);
		//layer 4
		((TH1F*)hdata)->SetBinContent(33+ic,counters[ic][1]);
	      }
	    }
	  }	
	}
      }

      delete hmodule;
      hmodule=NULL;
      for(Int_t i=0;i<2;i++) {
	delete hlayer[i];
	hlayer[i]=NULL;
      }

    }//end recpoint list uid = 20
    if(uidrec==40){
      //digitsr
      if (list->GetEntries() == 0){ 
	sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	break;
      } else{
	    
	while( (hdata = dynamic_cast<TH1* >(next())) ){
	  if (hdata){
	    if(hdata->GetEntries()==0)sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];
	    else {
	      TString hname=hdata->GetName();
	      if(hname.Contains("SDD DIGITS Module Pattern")) {
		//see raws
		
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDD Anode Distribution")) {
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDD Tbin Distribution")) {
		//to do as rp
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDD ADC Counts Distribution")) {
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      }//end adc counts
	      
	    }//end entries !=0
	  }//end hdata
	}//end while
      }//end else
      sddQACheckerValue=1.;
    }
    
  }
    break;
  case AliQAv1::kANA:
    {
      AliInfo(Form("===================> No Check on %s",AliQAv1::GetAliTaskName(index)));
      sddQACheckerValue=1.; 
    }
    break;
  case AliQAv1::kESD:
    {
      AliInfo(Form("==================>  No Check on %s",AliQAv1::GetAliTaskName(index)));
      sddQACheckerValue=1.;
    } 
    break;
  case AliQAv1::kNTASK:{
    AliInfo(Form("==================>  No Check on %s",AliQAv1::GetAliTaskName(index))); 
    sddQACheckerValue=1.;
  }
    break;
  case AliQAv1::kSIM:{
    AliInfo(Form("Check on %s",AliQAv1::GetAliTaskName(index))); 
    Int_t uid=list->GetUniqueID();
    if(uid==60) {
      //digits
      if (list->GetEntries() == 0){ 
	sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	break;
      } else{
	
	while( (hdata = dynamic_cast<TH1* >(next())) ){
	  if (hdata){
	    if(hdata->GetEntries()==0)sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];
	    else {
	      TString hname=hdata->GetName();
	      if(hname.Contains("SDDDIGITSModulePattern")) {
		//see raws
		
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDDAnodeDistribution")) {
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDDTbinDistribution")) {
		//to do as rp
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      } else if(hname.Contains("SDDADCCountsDistribution")) {
		sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
	      }//end adc counts
	      
	    }//end entries !=0
	  }//end hdata
	}//end while
      }//end else
    } else if(uid==50) 
      {
	//hits
	if (list->GetEntries() == 0){ 
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	  break;
	} 
	else{
	  
	  while( (hdata = dynamic_cast<TH1* >(next())) ){
	    if (hdata){
	      if(hdata->GetEntries()==0)sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];
	      else {
		TString hname=hdata->GetName();
		if(hname.Contains("SDDHITSModulePattern")) {
		  //to do as raws
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDHITlenghtalonglocalYCoord")) {
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDHITlenghtalonglocalYCoordZoom")) {
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDDepositedEnergyDistribution")) {
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		}//end deposited energy
		
	      }//end entries !=0
	    }//end hdata
	  }//end while
	}//end else
      } else if(uid==70) 
      {
	//sdigits
	if (list->GetEntries() == 0){ 
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	  break;
	} else{
	  
	  while( (hdata = dynamic_cast<TH1* >(next())) ){
	    if (hdata){
	      if(hdata->GetEntries()==0)sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];
	      else {
		TString hname=hdata->GetName();
		if(hname.Contains("SDDSDIGITSModulePattern")) {
		  //to do as raws
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDAnodeDistribution")) {
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDTbinDistribution")) {
		  //to do as rp
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		} else if(hname.Contains("SDDADCCountsDistribution")) {
		  sddQACheckerValue += fStepBitSDD[AliQAv1::kINFO];    
		}//end adc counts bindistribution
	      }//end entries !=0
			
	    }//end hdata
	  }//end while
	}//end else
      }//end sdigits
    sddQACheckerValue=1.;
  }
    break;
    
  }//end switch
  
  fCalibration=NULL;
  return sddQACheckerValue;	
}


//__________________________________________________________________
void AliITSQASDDChecker::FillCounters(Int_t counters[kNumOfSDDCheckerCounters][3], TH1* hmodule, TH2* hlay3, TH2* hlay4){
  // fill the counters of module filling
  for(Int_t imod=0;imod<fgknSDDmodules;imod++){
    Int_t lay=0;
    Int_t lad=0;
    Int_t det=0;
    Int_t module=imod+fgkmodoffset;
    AliITSCalibrationSDD * cal=(AliITSCalibrationSDD*)fCalibration->At(imod);
    if(cal==0) continue;
    AliITSgeomTGeo::GetModuleId(module,lay,lad,det);
    if (cal->IsBad()){
      counters[kExcludedMod][lay-3]++;
      Double_t content=0.;
      if(hmodule) content=hmodule->GetBinContent(imod+1);//if expert bit is active the histogram will be stored in the QA file otherwise the histogram will not be written on the logbook
      if(AliITSQADataMakerRec::AreEqual(content,0.)== kFALSE) {
	counters[kFilledMod][lay-3]++;
	counters[kExcludedButFilledMod][lay-3]++;
	AliError(Form("The module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data  ",module,lay,lad,det));
      }else{
	counters[kEmptyMod][lay-3]++; //it has to be empty
      }
    } else { // module is good
      Int_t sideFilled=0;
      Int_t sideActive=0;
      for(Int_t i=0;i<2;i++){
	if(cal->IsWingBad(i)==kFALSE) sideActive+=(1<<i);
	Double_t counts=0.;
	if(lay==3 && hlay3) counts=hlay3->GetBinContent(2*det+i-1,lad);
	else if(lay==4 && hlay4) counts=hlay4->GetBinContent(2*det+i-1,lad);
	if(counts>0.) sideFilled+=(1<<i);
      }
      if(sideActive==3) counters[kActiveMod][lay-3]++;
      else if(sideActive==0) counters[kExcludedMod][lay-3]++;
      else{
	counters[kActiveWing][lay-3]++;
	counters[kExcludedWing][lay-3]++;	    
      }
      if(sideFilled==3) counters[kFilledMod][lay-3]++;
      else if(sideFilled==0) counters[kEmptyMod][lay-3]++;
      else{
	counters[kEmptyWing][lay-3]++;
	counters[kFilledWing][lay-3]++;
      }      
      if(sideFilled!=sideActive){ // error cases
	if(sideActive==3){
	  if(sideFilled==0){
	    AliWarning(Form("The  module %d (layer %i, ladder %i det %i ) is in acquisition, but it didn't take data ",module, lay, lad, det)); 
	    counters[kActiveButEmptyMod][lay-3]++;
	  }
	  else if(sideFilled ==1 || sideFilled==2){
	    AliWarning(Form("The side %d of the module %d (layer %i, ladder %i det %i ) is in acquisition, but it didn't take data ",2-sideFilled,module, lay, lad, det));
	    counters[kActiveButEmptyWing][lay-3]++;
	  }
	}else if(sideActive==2 || sideActive==1){
	  if(sideFilled==0){
	    AliWarning(Form("The side %d of the module %d (layer %i, ladder %i det %i ) is in acquisition, but it didn't take data ",sideActive-1,module, lay, lad, det));
	    counters[kActiveButEmptyWing][lay-3]++;
	  }else if(sideFilled==1 || sideFilled==2){
	    AliWarning(Form("The side %d of the module %d (layer %i, ladder %i det %i ) is in acquisition, but it didn't take data ",sideActive-1,module, lay, lad, det));
	    AliError(Form("The side %d of the module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data",sideFilled-1,module,lay,lad,det));
	    counters[kExcludedButFilledWing][lay-3]++;
	    counters[kActiveButEmptyWing][lay-3]++;
	  }else if(sideFilled==3){
	    AliError(Form("The side %d of the module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data",2-sideActive,module,lay,lad,det));
	    counters[kExcludedButFilledWing][lay-3]++;		
	  }
	}else if(sideActive==0){
	  if(sideFilled==1 || sideFilled==2){
	    AliError(Form("The side %d of the module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data ",sideFilled-1,module,lay,lad,det));
	    counters[kExcludedButFilledWing][lay-3]++;
	  }else if(sideFilled==3){
	    AliError(Form("The module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data",module,lay,lad,det));
	    counters[kExcludedButFilledMod][lay-3]++;
	  }
	}
      }
    }
  }
}

//__________________________________________________________________
Int_t AliITSQASDDChecker::CheckCounters(Int_t counters[kNumOfSDDCheckerCounters][3], Int_t jl, Int_t nevents, TPaveText *ptext){
  // check the counters of module/hybrid filling

  Int_t numlimit=1000;
  if(AliRecoParam::ConvertIndex(GetEventSpecieForCheck())!=AliRecoParam::kCosmic) numlimit=10000;

  TString results1="";
  TString results2="";
  Int_t color=kBlack;
  Int_t textcolor=kBlack;
  Int_t retval=0;
  if(jl>=3) return retval;
  Bool_t warn=kFALSE;
  Bool_t err=kFALSE;
  if(jl<2) results1.Append(Form(" Layer %d:",jl+3));

  if(counters[kExcludedButFilledMod][jl]==0 && counters[kActiveButEmptyMod][jl]==0 && counters[kExcludedButFilledWing][jl]==0 && counters[kActiveButEmptyWing][jl]==0){
    AliInfo(Form("All the active modules (%i) and single drift regions (%i) are in acquisition. The number of excluded modules are %i and the excluded single drift regions are %i",counters[kActiveMod][jl],counters[kActiveWing][jl],counters[kExcludedMod][jl],counters[kExcludedWing][jl]));
    results1.Append(" OK.");
    results2.Append(" All active modules and drift regions in acquisition took data");
    color=kGreen;
    textcolor=kBlack;
    retval=fHighSDDValue[AliQAv1::kINFO];
  }else{
    if(counters[kActiveButEmptyMod][jl]>0){
      results1.Append(Form(" %i good module(s) didn't take data",counters[kActiveButEmptyMod][jl]));
      warn=kTRUE;
    }
    if(counters[kActiveButEmptyWing][jl]>0){
      results1.Append(Form(" %i good drift regions didn't take data",counters[kActiveButEmptyWing][jl]));
      warn=kTRUE;
    }
    if(warn){
      if(nevents<numlimit) {
	results2.Form(" Events %d .Too few events. DO NOT CALL the Expert ",nevents);
      }else{
	err=kTRUE;
 	results2.Form(" Events %d. If PHYSICS, follow the TWiki instruction and call the Expert ",nevents);
      }
    }
    if(counters[kExcludedButFilledMod][jl]!=0){
      results1.Append(Form(" %i modules excluded from the acquisition took data",counters[kExcludedButFilledMod][jl]));
      results2="Follow the TWiki instructions and Call the SDD expert ";
      err=kTRUE;
    }
    if(counters[kExcludedButFilledWing][jl]!=0){
      results1.Append(Form(" %i drift regions  excluded from the acquisition took data",counters[kExcludedButFilledWing][jl]));
      results2="Follow the TWiki instructions and Call the SDD expert ";
      err=kTRUE;
    }
    if(err){
      color=kRed;
      textcolor=kWhite;
      retval=fHighSDDValue[AliQAv1::kERROR];
      AliError(results1.Data());
    }  
    if(warn && !err){
      color=kYellow;
      textcolor=kBlack;
      retval=fHighSDDValue[AliQAv1::kWARNING];
      AliWarning(results1.Data());
    }
  } 
  if(ptext) {
    ptext->Clear();
    ptext->AddText(results1.Data());
    ptext->AddText(results2.Data());
    ptext->SetFillColor(color);
    ptext->SetFillStyle(1001);
    ptext->SetTextColor(textcolor);
  }

  return retval;
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
  //if(fStepBitSDD){/*delete fStepBitSDD;*/ fStepBitSDD=NULL;}
  //fStepBitSDD = new Double_t[AliQAv1::kNBIT];
  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fStepBitSDD[bit]=steprange[bit];
    }
}

//__________________________________________________________________
void  AliITSQASDDChecker::SetSDDLimits(const Float_t *lowvalue, const Float_t * highvalue)
{
  //set the low and high values in for each QA bit range calculated in the AliITSQAChecker class
  //  fLowSDDValue = new Float_t[AliQAv1::kNBIT];
  //  fHighSDDValue= new Float_t[AliQAv1::kNBIT];

  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fLowSDDValue[bit]=lowvalue[bit];
      fHighSDDValue[bit]= highvalue[bit];
    }

}
//__________________________________________________________________
Bool_t  AliITSQASDDChecker::MakeSDDImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode)
{
  //create the image for raws and recpoints. In the other case, the default methodof CheckerBase class will be used
  //
  Bool_t rval=kFALSE;
  fImage=(TCanvas**)AliQAChecker::Instance()->GetDetQAChecker(0)->GetImage();
  switch(task)
    {
    case AliQAv1::kRAWS:{
      rval=DrawHistos(list, task,mode); 
    }
      break;
    case AliQAv1::kRECPOINTS:{ 
      rval=DrawHistos(list, task,mode);
    }
      break;
    case AliQAv1::kHITS:; case AliQAv1::kESDS:; case AliQAv1::kDIGITS:;case AliQAv1::kDIGITSR:;case AliQAv1::kSDIGITS:;case AliQAv1::kTRACKSEGMENTS:;case AliQAv1::kRECPARTICLES:; default:
      {
	rval=kFALSE;
	//AliQAChecker::Instance()->GetDetQAChecker(0)->MakeImage(list,task,mode);
      }
      break;
    case AliQAv1::kNULLTASKINDEX:; case  AliQAv1::kNTASKINDEX: 
      {
	Int_t ts=(Int_t)task;
	AliWarning(Form("No histograms for this task number %d", ts)); 
	rval=kFALSE;
      }
      break;
    }
  return rval;  
}


//_______________________________________________________________________
Bool_t AliITSQASDDChecker::DrawHistos(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode)
{
  // MakeSDDRawsImage: raw data QA plots

  for (Int_t esIndex = 0 ; esIndex < AliRecoParam::kNSpecies ; esIndex++) {
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(esIndex)) || list[esIndex]->GetEntries() == 0)  {
      continue;
    } else {
      const Char_t * title = Form("QA_%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(esIndex)) ; 
      if ( !fImage[esIndex] ) {
	fImage[esIndex] = new TCanvas(title, title,1280,980) ;
      }
      fImage[esIndex]->Clear() ; 
      fImage[esIndex]->SetTitle(title) ; 
      fImage[esIndex]->cd();
      TPaveText someText(0.015, 0.015, 0.98, 0.98);
      someText.AddText(title);
      someText.Draw(); 
      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
      fImage[esIndex]->Clear() ; 
      if(task == AliQAv1::kRAWS) fImage[esIndex]->Divide(2,3);
      else fImage[esIndex]->Divide(2,6);
      TIter nexthist(list[esIndex]) ; 
      TH1* hist = NULL ;
      Int_t npad = 1 ; 
      fImage[esIndex]->cd(npad); 
      fImage[esIndex]->cd(npad)->SetBorderMode(0) ;
      while ( (hist=static_cast<TH1*>(nexthist())) ) {
	TString hname(hist->GetName());
	TString cln(hist->ClassName()) ; 
	if ( ! cln.Contains("TH") ) continue ;
	if(hist->TestBit(AliQAv1::GetImageBit())) {
	  hist->GetXaxis()->SetTitleSize(0.04);
	  hist->GetYaxis()->SetTitleSize(0.04);
	  hist->GetXaxis()->SetLabelSize(0.02);
	  hist->GetYaxis()->SetLabelSize(0.02);
	  if(cln.Contains("TH1") && task == AliQAv1::kRECPOINTS){
	    if(!hname.Contains("Check")) hist->SetFillColor(kOrange+7);
	  }
	  if(cln.Contains("TH2")) {
	    gPad->SetRightMargin(0.15);
	    gPad->SetLeftMargin(0.05);
	    hist->SetStats(0);
	    hist->SetOption("colz") ;
	  }
	  hist->DrawCopy(); 
	  fImage[esIndex]->cd(++npad) ; 
	  fImage[esIndex]->cd(npad)->SetBorderMode(0) ; 
	} 
      }
      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
    }
  }
  return kTRUE;
}
