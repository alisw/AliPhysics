
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
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>



// --- AliRoot header files ---
#include "AliITSQASDDChecker.h"
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"
#include "AliQACheckerBase.h"
#include "TSystem.h"
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
Double_t AliITSQASDDChecker::Check(AliQAv1::ALITASK_t index, TObjArray * list, const AliDetectorRecoParam * /*recoparam*/) 
{  
  AliInfo(Form("AliITSQASDDChecker called with offset: %d\n", fSubDetOffset) );

  AliDebug(1,Form("AliITSQASDDChecker called with offset: %d\n", fSubDetOffset));

  Double_t test = 0.;
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
	test= fHighSDDValue[AliQAv1::kWARNING];
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

  switch(index)
    {

    case AliQAv1::kRAW:
      AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index)));
      
      if (list->GetEntries() == 0.){ //check if the list is empty
	//printf("test = %f \t value %f\n",test,fHighSDDValue[AliQAv1::kFATAL]);
	test=test+fHighSDDValue[AliQAv1::kFATAL];
	break;
      }//end if getentries
      else{
	TIter next(list);
	Int_t offset = 0;
	for(offset =0;offset < fSubDetOffset; offset++){hdata = dynamic_cast<TH1*>(next());}//end for
	Int_t emptymodules[2];
	Int_t filledmodules[2];
	Int_t emptyladders[2];
	Int_t filledladders[2];
	for(Int_t i=0;i<2;i++){
	  emptymodules[i]=0;
	  filledmodules[i]=0;
	  emptyladders[i]=0;
	  filledladders[i]=0;
	}
	TH1 *hmodule=NULL;
	TH2 *hlayer[2];
	for(Int_t i=0;i<2;i++)hlayer[i]=NULL;    
	while( hdata = dynamic_cast<TH1* >(next()) ){
	  if (hdata){
	    TString hname=hdata->GetName();
	    if(hname.Contains("SDDchargeMap"))continue;
	    if(hname.Contains("SDDModPattern")){
	      hmodule=hdata;
	      entries= hdata->GetEntries();
	      if(entries==0){
		AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));
		//printf("test = %f \t value %f\n",test,fHighSDDValue[AliQAv1::kFATAL]);
		test=test+fStepBitSDD[AliQAv1::kFATAL];
	      }//endif entries
	      else{
		int modmax=hdata->GetNbinsX();
		Int_t empty=0;
		Int_t filled=0;
		Double_t content=0;
		for(Int_t i=1;i<=modmax;i++){
		  content=hdata->GetBinContent(i);
		  if(content==0.){empty++;}
		  else if(content!=0.){filled++;}
		}//end for
		AliInfo(Form(" %s : empty modules %i \t filled modules %i",hname.Data(), empty, filled));
	      }//end else pattern entries !=0	      
	    }//end if modpattern
	    else if(hname.Contains("SDDphizL3")||hname.Contains("SDDphizL4")){
	      Int_t layer=0;
	      if(hname.Contains("3"))layer=0;
	      else  if(hname.Contains("4"))layer=1;
	      entries2[layer]=hdata->GetEntries();
	      if(entries2[layer]==0){
		AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));
		//printf("test = %f \t value %f\n",test,fStepBitSDD[AliQAv1::kFATAL]);
		test=test+fStepBitSDD[AliQAv1::kFATAL];
		if(entries==0){ 
		  //return test; 
		  //break;
		}
	      }//end if getentries
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
		  if(ladcontent==0){emptyladders[layer1]++;}
		  else if(ladcontent!=0){filledladders[layer1]++;} 
		}//end for
		AliInfo(Form(" %s : empty ladders %i \t filled ladders %i\n",hname.Data(), emptyladders[layer], filledladders[layer]));//end else layer 3
		delete hproj;
		hproj=NULL;	
		//delete htemp;
		//htemp=NULL;
	      }//end else entries !=0	      
	    }//end if layer 3
	  }//end if hdata	
	}//end while
	if(entries==0.&&entries2[0]==0.&&entries2[1]==0.) break;
	else{
	  if(hmodule||(hlayer[0]&&hlayer[1])){
	    Int_t excluded=0;
	    Int_t active=0;
	    Int_t exactive=0;//excluded but taking data
	    //AliITSCalibrationSDD *cal;
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
		if(content!=0.||contentlayer[lay-3]!=0.)
		  {
		    filledmodules[lay-3]++;
		    AliWarning(Form("The module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data \n ",module,lay,lad,det));
		    exactive++;
		  }
		else if(content==0.&&contentlayer[lay-3]==0.)emptymodules[lay-3]++;
		//AliInfo(Form("The module %d (layer %i, ladder %i det %i ) is bad, content %f content layer %f  filled modules position %d ",module,lay,lad,det,contentlayer[lay-3],content,lay-3) );
	      }//end if bad
	      else
		{
		  Double_t contentgood=0.;
		  active++;
		  //printf("lay: %i\t det %i \t lad %i \n",lay,det,lad );
		  contentgood=hlayer[lay-3]->GetBinContent(det,lad);
		  if(contentgood==0.){emptymodules[lay-3]++;}
		  else if(contentgood!=0.){filledmodules[lay-3]++;}
		}
	      
	      //delete cal;
	      //cal=NULL;
	    }//end for
	    for(Int_t i=0;i<2;i++){AliInfo(Form("Layer %i \tempty modules %i \t filled modules %i\n", i+3,emptymodules[i], filledmodules[i]));}//end else layers
	    if(exactive==0){
	      AliInfo(Form("All the active modules (%i) are in acquisition. The number of excluded modules are %i \n",active,excluded));
	      test=fHighSDDValue[AliQAv1::kINFO];}
	    if(exactive!=0){
	      AliWarning(Form("%i modules excluded from the acquisition took data. Active modules%i \n ",exactive,active));
	      test=fHighSDDValue[AliQAv1::kWARNING];
	    }
	    if(excluded==exactive){
	      AliWarning(Form("All the modules exluded from the acquisition (%d) took data!  Active modules %i\n",excluded,active));
	      test=fHighSDDValue[AliQAv1::kWARNING];
	    }
	    if(active==0){
	      AliWarning(Form("No modules took data: excluded %i \t exactive %i \n", excluded, exactive)); 
	      test=fHighSDDValue[AliQAv1::kFATAL];
	    }
	    for(Int_t i=0;i<2;i++)
	      {
		delete hlayer[i];
		hlayer[i]=NULL;
	      }
	  }//end else 
	}
      }//end getentries !=0
      //delete calSDD;
      
      //delete calibSDD;
      //delete calSDD;
   
      break;
    case AliQAv1::kNULLTASK:
      AliInfo(Form("No Check on %s\n",AliQAv1::GetAliTaskName(index))); 
      test=1.;
      break;
    case AliQAv1::kREC:
      /*
      AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index))); 
      //TH1*hdata=NULL;
      if(list->GetUniqueID()==40){
	if (list->GetEntries() == 0.){ //check if the list is empty
	  //printf("test = %f \t value %f\n",test,fHighSDDValue[AliQAv1::kFATAL]);
	  test=fHighSDDValue[AliQAv1::kFATAL];
	  
	}//end if getentries
	else{
	  
	  TIter next(list);
	
	  while( hdata = dynamic_cast<TH1* >(next()) ){
	    if (hdata){
	      if(hdata->GetEntries()==0)test=test+fStepBitSDD[AliQAv1::kFATAL];
	      else
		{
		  TString hname=hdata->GetName();
		  if(hname.Contains("FSE"))continue;
		  else if(hname.Contains("SDDLay3TotCh")||hname.Contains("SDDLay4TotCh")){
		    Double_t meancharge=hdata->GetMean();
		    Double_t rmscharge=hdata->GetRMS();
		    AliInfo(Form("%s : Mean value:%f RMS value%f \n ",hname.Data(),meancharge,rmscharge));
		    test=test+fStepBitSDD[AliQAv1::kINFO];    
		  }//end if name charge
		  else if(hname.Contains("SDDGlobalCoordDistribYX" ))
		    {
		      test=test+fStepBitSDD[AliQAv1::kINFO];    
		    }//end xy
		  else if(hname.Contains("SDDGlobalCoordDistribRZ"))
		    {
		      
		      test=test+fStepBitSDD[AliQAv1::kINFO];    
		    } //end rz
		  else if(hname.Contains("SDDGlobalCoordDistribL3PHIZ" )||hname.Contains("SDDGlobalCoordDistribL3PHIZ"))
		    {    

		    }//end phiz
		  else if(hname.Contains("SDDModPatternRP"))
		    {
		      
		      //to do :se raws

		    }//modpattern
		  else if(hname.Contains("SDDModPatternL3RP")||hname.Contains("SDDModPatternL4RP") )
		    {
		      //to do: see raws
		    }//end ladpattern
		  else if(hname.Contains("SDDLocalCoordDistrib"))
		    {
		      test=test+fStepBitSDD[AliQAv1::kINFO];    
		    }//end local coord
		  else if(hname.Contains("SDDrdistrib_Layer3")||hname.Contains("SDDrdistrib_Layer4"))
		    {
		      
		    }//end r distribution
		  else if(hname.Contains("SDDphidistrib_Layer3")||hname.Contains("SDDphidistrib_Layer4"))
		    {
		      
		    }//end phi distribution
		  else if(hname.Contains("SDDdrifttime_Layer3")||hname.Contains("SDDdrifttime_Layer4"))
		    {
		      
		    }//end drift time
		}
	    }//end if hdata
	    
	  }//end while
	}//end else geentries
      }//end uniqueid
      */
      test=1.;
      break;
    case AliQAv1::kANA:
      AliInfo(Form("===================> No Check on %s\n",AliQAv1::GetAliTaskName(index)));
      test=1.; 
      break;
    case AliQAv1::kESD:
      AliInfo(Form("==================>  No Check on %s\n",AliQAv1::GetAliTaskName(index)));
      test=1.; 
      break;
    case AliQAv1::kNTASK:
      AliInfo(Form("==================>  No Check on %s\n",AliQAv1::GetAliTaskName(index))); 
      test=1.;
      break;
    case AliQAv1::kSIM:
      AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index))); 
      Int_t uid=list->GetUniqueID();
      if(uid==60)
	{
	  //digits
	  if (list->GetEntries() == 0.){ //check if the list is empty
	    //printf("test = %f \t value %f\n",test,fHighSDDValue[AliQAv1::kFATAL]);
	    test=fHighSDDValue[AliQAv1::kFATAL];
	    
	  }//end if getentries
	  else{
	    
	    TIter next(list);
	    
	    while( hdata = dynamic_cast<TH1* >(next()) ){
	      if (hdata){
		if(hdata->GetEntries()==0)test=test+fStepBitSDD[AliQAv1::kFATAL];
		else
		  {
		    TString hname=hdata->GetName();
		    if(hname.Contains("SDDDIGITSModulePattern"))
		      {
			//see raws

			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end modpattern
		    else if(hname.Contains("SDDAnodeDistribution"))
		      {
		       
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end anode distribution
		    else if(hname.Contains("SDDTbinDistribution"))
		      {

			//to do as rp
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end timebindistribution
		    else if(hname.Contains("SDDADCCountsDistribution"))
		      {
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end adc counts

		  }//end entries !=0
	      }//end hdata
	    }//end while
	  }//end else
	}//end digits
      else if(uid==50)
	{
	  //hits
	  if (list->GetEntries() == 0.){ //check if the list is empty
	    //printf("test = %f \t value %f\n",test,fHighSDDValue[AliQAv1::kFATAL]);
	    test=fHighSDDValue[AliQAv1::kFATAL];
	    
	  }//end if getentries
	  else{
	    
	    TIter next(list);
	    
	    while( hdata = dynamic_cast<TH1* >(next()) ){
	      if (hdata){
		if(hdata->GetEntries()==0)test=test+fStepBitSDD[AliQAv1::kFATAL];
		else
		  {
		    TString hname=hdata->GetName();
		    if(hname.Contains("SDDHITSModulePattern"))
		      {
			//to do as raws
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end modpattern
		    else if(hname.Contains("SDDHITlenghtalonglocalYCoord"))
		      {
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end hit lenght
		    else if(hname.Contains("SDDHITlenghtalonglocalYCoordZoom"))
		      {
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end hit lenght
		    else if(hname.Contains("SDDDepositedEnergyDistribution"))
		      {
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end deposited energy

		  }//end entries !=0
	      }//end hdata
	    }//end while
	  }//end else
	}//end hits
      else if(uid==70)
	{
	  //sdigits
	  if (list->GetEntries() == 0.){ //check if the list is empty
	    //printf("test = %f \t value %f\n",test,fHighSDDValue[AliQAv1::kFATAL]);
	    test=fHighSDDValue[AliQAv1::kFATAL];
	    
	  }//end if getentries
	  else{
	    
	    TIter next(list);
	    
	    while( hdata = dynamic_cast<TH1* >(next()) ){
	      if (hdata){
		if(hdata->GetEntries()==0)test=test+fStepBitSDD[AliQAv1::kFATAL];
		else
		  {
		    TString hname=hdata->GetName();
		    if(hname.Contains("SDDSDIGITSModulePattern"))
		      {
			//to do as raws
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end modpattern
		    else if(hname.Contains("SDDAnodeDistribution"))
		      {
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end anode bindistribution
		    else if(hname.Contains("SDDTbinDistribution"))
		      {
			//to do as rp
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end timebindistribution
		    else if(hname.Contains("SDDADCCountsDistribution"))
		      {
			test=test+fStepBitSDD[AliQAv1::kINFO];    
		      }//end adc counts bindistribution

		  }//end entries !=0
	      }//end hdata
	    }//end while
	  }//end else
	}//end sdigits
      test=1.;
      break;
      
    }//end switch

  fCalibration=NULL;
  delete hdata;
  return test;	
}
 
//__________________________________________________________________
void AliITSQASDDChecker::SetTaskOffset(Int_t taskoffset)
{
  fSubDetOffset = taskoffset;
}


//__________________________________________________________________
void AliITSQASDDChecker::SetStepBit(Double_t *steprange)
{

  fStepBitSDD = new Double_t[AliQAv1::kNBIT];
  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fStepBitSDD[bit]=steprange[bit];
    }
}

//__________________________________________________________________
void  AliITSQASDDChecker::SetSDDLimits(Float_t *lowvalue, Float_t * highvalue)
{

  fLowSDDValue = new Float_t[AliQAv1::kNBIT];
  fHighSDDValue= new Float_t[AliQAv1::kNBIT];

  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fLowSDDValue[bit]=lowvalue[bit];
      fHighSDDValue[bit]= highvalue[bit];
    }

}


