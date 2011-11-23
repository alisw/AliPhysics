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

AliITSQASDDChecker::AliITSQASDDChecker():
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
	for(Int_t ibit=0;ibit<AliQAv1::kNBIT;ibit++)
	  {
	    fStepBitSDD[ibit]=0.;
	    fLowSDDValue[ibit]=0.;
	    fHighSDDValue[ibit]=0.;
	  }

}          // ctor



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
  if(fImage)
    {
      delete []fImage; 
      fImage=NULL;
    }
} // dtor

//__________________________________________________________________
Double_t AliITSQASDDChecker::Check(AliQAv1::ALITASK_t index, const TObjArray * list, const AliDetectorRecoParam * /*recoparam*/) 
{
  //check histograms of the different lists  
  AliInfo(Form("AliITSQASDDChecker called with offset: %d\n", fSubDetOffset) );

  AliDebug(1,Form("AliITSQASDDChecker called with offset: %d\n", fSubDetOffset));

  Double_t sddQACheckerValue = 0.;
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
  }//end f calibration

  AliInfo("Calib SDD Created\n ");

  TIter next(list);
  TString results1;
  TString results2;
  Int_t color=1;

  switch(index) {
    case AliQAv1::kRAW:{
      AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index)));
//       if(fRawModulePattern) { delete fRawModulePattern; fRawModulePattern = 0; }
//       if(fRawL3Pattern) { delete fRawL3Pattern; fRawL3Pattern = 0; }
//       if(fRawL4Pattern) { delete fRawL4Pattern; fRawL4Pattern = 0; }
      if (list->GetEntries() == 0){AliError("Raw List for SDD is empty \n");sddQACheckerValue += fHighSDDValue[AliQAv1::kFATAL];	break;}
      TH1 *hmodule=NULL;
      TH2 *hlayer[2]; 
      hdata=NULL;
      for(Int_t i=0;i<2;i++)hlayer[i]=NULL;

      //check counters
      Int_t emptymodules[2], filledmodules[2],emptyladders[2],filledladders[2],emptydriftregion[2], filleddriftregion[2], excludedmoduleperlayer[2], excludeddrperlayer[2], activemoduleperlayer[2],activedrperlayer[2],exactivemoduleperlayer[2],exactivedrperlayer[2];
      Int_t excluded=0;            //excluded modules  
      Int_t excludeddriftregion=0; //excluded single drift region
      Int_t active=0;              //active modules  
      Int_t activedriftregion=0;   //active single drift region
      Int_t exactive=0;            //excluded modules but taking data
      Int_t exactivedriftregion=0; //excluded single drift region but taking data   
      Int_t empty=0;
      Int_t filled=0;
      Int_t emptydr=0;
      Int_t filleddr=0;
      //Int_t emptyactivemodule=0;
      Int_t emptyactivemoduleperlayer[2];
      //Int_t emptydractivemodule=0;
      Int_t emptyactivedrperlayer[2];
      Int_t emptysum=0;
      Int_t emptydiff=0;
      Int_t emptydrsum=0;
      Int_t emptydrdiff=0;

    for(Int_t i=0;i<2;i++)
      {
	emptymodules[i]=0; 
	filledmodules[i]=0; 
	emptyladders[i]=0; 
	filledladders[i]=0; 
	emptydriftregion[i]=0;
	filleddriftregion[i]=0; 
	excludedmoduleperlayer[i]=0;
	excludeddrperlayer[i]=0;
	activemoduleperlayer[i]=0;
	activedrperlayer[i]=0;
	exactivemoduleperlayer[i]=0;
	exactivedrperlayer[i]=0;
	emptyactivemoduleperlayer[i]=0;
	emptyactivedrperlayer[i]=0;
      }  

    Int_t neventsraw=0;

    //take the number of events

      while( (hdata = dynamic_cast<TH1* >(next())) ){
	if (hdata){
	  TString hname=hdata->GetName();

	  if(hname.Contains("SDDRawDataCheck"))
	    { 
	      neventsraw=(Int_t)hdata->GetBinContent(1); 
	      //break;
	    }
	  else continue;
	}//end if hdata
      }//end while

      next.Begin();
      while( (hdata = dynamic_cast<TH1* >(next())) ){
	if (hdata){TString hname=hdata->GetName();
	  if(hname.Contains("SDDchargeMap"))continue;
	  if(hname.Contains("SDDDDLPattern"))continue;
	  if(hname.Contains("SDDEventSize"))continue;
	  if(hname.Contains("_RelativeOccupancy"))continue;
	  if(hname.Contains("SDDModPattern")){
	    if(hname.Contains("NORM")) continue;
	    hmodule=(TH1*)hdata->Clone();
	    entries= hdata->GetEntries();
	    if(AliITSQADataMakerRec::AreEqual(entries,0.)){AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];}//endif entries
	    else{
	      int modmax=hdata->GetNbinsX();
	      Double_t content=0;
	      Int_t llay=0;
	      for(Int_t i=1;i<=modmax;i++)
		{
		  if(i<85)llay=0;
		  else llay=1;
		  content=hdata->GetBinContent(i);
		  if(AliITSQADataMakerRec::AreEqual(content,0.))
		    {
		    empty++;
		    emptymodules[llay]++;
		    } 
		  else 
		    {
		      filled++;
		      filledmodules[llay]++;
		    }
		}//end for
	      //output of the check at the level of the modules. Drift region in the following checks

	      AliInfo(Form(" %s : empty modules %i \t filled modules %i",hname.Data(), empty, filled));
	      AliInfo(Form(" %s : Layer 3 empty modules %i \t filled modules %i",hname.Data(), emptymodules[0], filledmodules[0]));
	      AliInfo(Form(" %s : Layer 4 empty modules %i \t filled modules %i",hname.Data(), emptymodules[1], filledmodules[1]));
	    }//end else pattern entries !=0
	  } //modpattern (1d histogram) 		
	  if(hname.Contains("_RelativeOccupancy")) {
	    //fRawModulePattern = (TH1F *) hdata;
	    Float_t threshold = hdata->GetMean() + 4*hdata->GetRMS();
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 1: L3 mean: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 2: L4 mean: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    Int_t aboveThreshold = 0;
	    for(Int_t k=0; k<= hdata->GetNbinsX(); k++) {if(hdata->GetBinLowEdge(k) > threshold) aboveThreshold += (int)(hdata->GetBinContent(k));}
	    Float_t fractionAboveThreshold=0.;
	    if(hdata->GetEntries()>0.)fractionAboveThreshold=((Float_t) aboveThreshold)/hdata->GetEntries();
	    if(hname.Contains("L3")) AliInfo(Form("SDD check number 1, L3: Raw fractionAboveThreshold: %f",fractionAboveThreshold));
	    if(hname.Contains("L4")) AliInfo(Form("SDD check number 2, L4: Raw fractionAboveThreshold: %f",fractionAboveThreshold));
	    if(fractionAboveThreshold > fThresholdForRelativeOccupancy) {sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	      if(hname.Contains("L3")) AliInfo(Form("SDD check number 1: Set Warning (L3 Raw)"));
	      if(hname.Contains("L4")) AliInfo(Form("SDD check number 2: Set Warning (L4 Raw)")); } }//relativeoccupancy

	  if(hname.Contains("SDDphizL3") || hname.Contains("SDDphizL4")){if(hname.Contains("NORM"))continue;
	    Int_t layer=0;
	    if(hname.Contains("3"))layer=0;
	    else  if(hname.Contains("4"))layer=1;
	    entries2[layer]=hdata->GetEntries();
	    if(entries2[layer]==0){AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));
	      sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];}//end if getentries
	    else{
	      Int_t layer1=0;
	      if(hname.Contains("3"))layer1=0;
	      else  if(hname.Contains("4"))layer1=1;
	      TH2* htemp=dynamic_cast<TH2*>(hdata);
	      if(htemp){
		hlayer[layer1]=(TH2*)htemp->Clone();
		hlayer[layer1]->SetName(Form("%s_copy",hname.Data()));
		//hlayer[layer1]->RebinX(2);
		int modmay=hlayer[layer1]->GetNbinsY();
		TH1D* hproj= hlayer[layer1]->ProjectionY();
		Double_t ladcontent=0;
		for(Int_t i=1;i<=modmay;i++) {//loop on the ladders
		  ladcontent=hproj->GetBinContent(i);
		  if(AliITSQADataMakerRec::AreEqual(ladcontent,0.)) emptyladders[layer1]++;
		  else filledladders[layer1]++;}//end for
		AliInfo(Form(" %s : empty ladders %i \t filled ladders %i\n",hname.Data(), emptyladders[layer1], filledladders[layer1]));//end else layer 3
		delete hproj;
		hproj=NULL;
	      }//end if htemp
	    }//end else entries !=0
	  }//end check on phiz	      
	}//end if hdata	
      }//end while
      for(Int_t ii=0;ii<2;ii++)
	{
	  filledmodules[ii]=0;
	  emptymodules[ii]=0;
	}
      filled=0;
      empty=0;
      if(AliITSQADataMakerRec::AreEqual(entries,0.)&& AliITSQADataMakerRec::AreEqual(entries2[0],0.)&& AliITSQADataMakerRec::AreEqual(entries2[1],0.)) break;
      //else{
      if(hmodule || (hlayer[0] && hlayer[1])){
	for(Int_t imod=0;imod<fgknSDDmodules;imod++){
	  Int_t lay=0;
	  Int_t lad=0;
	  Int_t det=0;
	  Int_t module=0;
	  module=imod+fgkmodoffset;
	  AliITSCalibrationSDD * cal=(AliITSCalibrationSDD*)fCalibration->At(imod);
	  if(cal==0) { continue;}
	  AliITSgeomTGeo::GetModuleId(module,lay,lad,det);
	  if (cal->IsBad()){
	    excluded++;
	    excludedmoduleperlayer[lay-3]++;
	    Double_t content=0.;
	    Double_t contentlayer[2];
	    for(Int_t i=0;i<2;i++)contentlayer[i]=0.;
	    if(hmodule)content=hmodule->GetBinContent(imod+1);//if expert bit is active the histogram will be stored in the QA file otherwise the histogram will not be written on the logbook
	    //	    if(hlayer[lay-3]) contentlayer[lay-3]=hlayer[lay-3]->GetBinContent(det,lad);
	    if(AliITSQADataMakerRec::AreEqual(content,0.)== kFALSE) {
	      filledmodules[lay-3]++;
	      filled++;
	      AliError(Form("The module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data \n ",module,lay,lad,det));
	      exactive++;
	      exactivemoduleperlayer[lay-3]++;
	    } else if(AliITSQADataMakerRec::AreEqual(content,0.)) 
	      {
		emptymodules[lay-3]++; //it has to be empty
		empty++;
	      }
	  } else {
	    Int_t totside=0;
	    Int_t totactiveside=0;
	    //Int_t totbadside=0;
	    Double_t contentgood=0.;

	    for(Int_t i=0;i<2;i++){
	      if(hlayer[lay-3]) contentgood=hlayer[lay-3]->GetBinContent(2*det+i-1,lad);
	      if(cal->IsWingBad(i))
		{ 
		  excludeddriftregion++; 
		  excludeddrperlayer[lay-3]++; 
		  if(AliITSQADataMakerRec::AreEqual(contentgood,0.)==kFALSE){
		    AliError(Form("The side %d of the module %d (layer %i, ladder %i det %i ) excluded from the acquisition, took data \n ",i,module,lay,lad,det));
		    exactivedriftregion++; 
		    exactivedrperlayer[lay-3]++;
		    filleddr++;
		    filleddriftregion[lay-3]++;
		  }
		}//end wingbad
	      else{
		if(AliITSQADataMakerRec::AreEqual(contentgood,0.)==kTRUE)
		  {
		    AliWarning(Form("The side %d of the module %d (layer %i, ladder %i det %i ) is in acquisition, but it didn't take data\n ",i,module, lay, lad, det));
		  }
		else 
		  {
		    totside++;
		  }
		totactiveside++;
	      }
	    }//end for
	    if(totside==0){
	    AliWarning(Form("The  module %d (layer %i, ladder %i det %i ) is in acquisition, but it didn't take data\n ",module, lay, lad, det));
	      emptymodules[lay-3]++;
	      empty++;

	    }
	      else 
		if(totside==2){
		filledmodules[lay-3]++;
		filled++;
		}
		else
		  if(totside==1)
		    {
		      //		      emptydr++;
		      //emptydriftregion[lay-3]++; //it has to be empty
		      emptydriftregion[lay-3]++; 
		      emptydr++;
		      filleddr++;
		      filleddriftregion[lay-3]++;
		    }
	    if(totactiveside==1)
	      {
		activedriftregion++;
		activedrperlayer[lay-3]++;
	      }else if(totactiveside==2)
	      {
		active++;
		activemoduleperlayer[lay-3]++;
	      }

	  }
	}//end for
	AliInfo(Form("In total %d modules and %d single drift regions took data.\n ",filled, filleddr));
	AliInfo(Form("In total %d modules and %d single drift regions were empty\n",empty, emptydr));
	for(Int_t i=0;i<2;i++)
	  {
	    AliInfo(Form("Layer %i   \tempty modules %i             \t filled modules %i\n", i+3,emptymodules[i], filledmodules[i]));
	    AliInfo(Form("Layer %i   \tempty single drift regions %i \t filled single drift regions %i\n",i+3,emptydriftregion[i], filleddriftregion[i]));
	  }//end else layers
	emptysum=emptymodules[0]+emptymodules[1];
	emptydiff=emptysum-excluded;
	emptyactivemoduleperlayer[0]=emptymodules[0]- excludedmoduleperlayer[0];
	emptyactivemoduleperlayer[1]=emptymodules[1]- excludedmoduleperlayer[1];

	emptydrsum=emptydriftregion[0]+emptydriftregion[1];
	emptydrdiff=emptydrsum-excludeddriftregion;
	emptyactivedrperlayer[0]=emptydriftregion[0]- excludeddrperlayer[0];
	emptyactivedrperlayer[1]=emptydriftregion[1]- excludeddrperlayer[1];

	Int_t numlimit=1000;

 	if(emptysum>excluded||emptydrsum>excludeddriftregion){
	  AliWarning(Form(" %i good module(s) and %i good single drift regions didn't take data! \n",emptydiff,emptydrdiff));
	  AliWarning(Form(" Layer 3: %i good module(s) and %i good single drift regions didn't take data! \n",emptyactivemoduleperlayer[0] ,emptyactivedrperlayer[0] ));
	  AliWarning(Form(" Layer 4: %i good module(s) and %i good single drift regions didn't take data! \n",emptyactivemoduleperlayer[1] ,emptyactivedrperlayer[1] ));
	  //printf("========================= %d",AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::kCosmic));
	  //	  if((AliQAv1::Instance(AliQAv1::GetDetIndex("ITS")))->IsEventSpecieSet(AliRecoParam::kCosmic)==kFALSE){     
	  if(AliRecoParam::Convert(GetEventSpecieForCheck())!=AliRecoParam::kCosmic){     
	    
	    results1.Form("%i good module(s) and %i good drift regions didn't take data!",emptydiff,emptydrdiff);
	    if(neventsraw<numlimit)
	      {
		results2.Form(" Events %d .Too few events.DO NOT CALL the Expert ",neventsraw);
		color=5;
		sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	      }
	    else
	      {
		results2.Form(" Events %d. If PHYISICS, follow the TWiki instruction and call the Expert ",neventsraw);
		color=2;
		sddQACheckerValue=fHighSDDValue[AliQAv1::kERROR];
	      }
	  }
	  else
	    //	    if((AliQAv1::Instance(AliQAv1::GetDetIndex("ITS")))->IsEventSpecieSet(AliRecoParam::kCosmic)==kTRUE)
	    if(AliRecoParam::Convert(GetEventSpecieForCheck())==AliRecoParam::kCosmic)
	      {     
		numlimit=10000;
		if(neventsraw<numlimit)
		  {
		    AliWarning(Form("This is a cosmic run. Some module and drift region are empty but all is OK. "));
		    results1.Form("OK. This is a cosmic run. you need a lot of events. Events %i",neventsraw);
		    results2.Form("%i good module(s) and %i good drift are empty! DO NOT CALL THE EXPERT",emptydiff,emptydrdiff);
		    color=5;
		    sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
		  }
		else
		  {
		    results1.Form("%i good module(s) and %i good drift region(s) have no data!",emptydiff,emptydrdiff);
		    results2.Form(" Cosmic Events %d .Follow the TWiki instruction and call the Expert ",neventsraw);
		    color=2;
		    sddQACheckerValue=fHighSDDValue[AliQAv1::kERROR];
		  }
	      }
	}
	
	if(exactive==0 && emptydiff==0 && exactivedriftregion==0 && emptydrdiff==0){
	  AliInfo(Form("All the active modules (%i) and single drift regions (%i) are in acquisition. The number of excluded modules are %i and the excluded single drift regions are %i\n",active,activedriftregion,excluded,excludeddriftregion));
	  results1.Form("OK.");
	  results2.Form(" All active modules and drift regions in acquisition");
	  color=3;
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kINFO];
	}
	if(exactive!=0||exactivedriftregion!=0){
	  AliError(Form("%i modules and %i single  drift regions excluded from the acquisition took data. Active modules%i single drift region %i \n ",exactive,exactivedriftregion,active,activedriftregion));
	  AliError(Form("Layer 3: %i modules and %i single  drift regions excluded from the acquisition took data. Active modules%i single drift region %i \n ",exactivemoduleperlayer[0],exactivedrperlayer[0],activemoduleperlayer[0],activedrperlayer[0]));
	  AliError(Form("Layer 3: %i modules and %i single  drift regions excluded from the acquisition took data. Active modules%i single drift region %i \n ",exactivemoduleperlayer[1],exactivedrperlayer[1],activemoduleperlayer[1],activedrperlayer[1]));
	  results1.Form("%i modules and %i drift region excluded from the acquisition took data",exactive,exactivedriftregion);
	  results2.Form("Follow the TWiki instructions and Call the SDD expert ");
	  color=2;	
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kERROR];
	}
	if(excluded==exactive||excludeddriftregion==exactivedriftregion){
	  AliError(Form("All the modules (%d) or single drift regions (%d) excluded from the acquisition  took data!\n  Active modules %i \t Active drfift regions %i\n",excluded,excludeddriftregion,active,activedriftregion));
	  results1.Form("All the modules (%d) or drift regions (%d) excluded from the acquisition took data!",excluded,excludeddriftregion );
	  results2.Form("Follow the TWiki instructions and Call the SDD expert ");
	  color=6;
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	}
	if(active==0||activedriftregion==0){
	  AliError(Form("No modules or single drift regions took data: excluded %i \t excluded active %i \n\t\t excluded single drift regions %i \t excluded active drift regions %i \n", excluded, exactive, excludeddriftregion, exactivedriftregion)); 
	  results1.Form("No modules or drift region took data: excluded modules %i  excluded drift regions %i ", excluded, excludeddriftregion );
	  results2.Form("Follow the TWiki instructions and Call the SDD expert ");
	  color=6;
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	}

	TPaveText *pave[2];
	next.Begin();

	while( (hdata=dynamic_cast<TH1* >(next())) )
	  {
	    if (hdata){
	      TString hname=hdata->GetName();
	      if(hname.Contains("SDDphizL3") || hname.Contains("SDDphizL4")){
		if(hname.Contains("NORM"))continue;
		//AliInfo("========================================Found histo 11\n");
		Int_t lay=0;
		if(hname.Contains("3"))lay=0;
		else if(hname.Contains("4"))lay=1;
		pave[lay]=new TPaveText(0.3,0.88,0.9,0.99,"NDC");
		pave[lay]->AddText(results1.Data());
		pave[lay]->AddText(results2.Data());
		pave[lay]->SetFillColor(color);
		pave[lay]->SetBorderSize(1);
		pave[lay]->SetLineWidth(1);
		hdata->GetListOfFunctions()->Add(pave[lay]);
	      }
	      else
		if(hname.Contains("SDDRawDataCheck"))
		  {
		    
		    //AliInfo("========================================Found histo\n");
		    ((TH1F*)hdata)->SetBinContent(5,active);
		    ((TH1F*)hdata)->SetBinContent(6,filled);
		    ((TH1F*)hdata)->SetBinContent(7,activedriftregion);
		    ((TH1F*)hdata)->SetBinContent(8,filleddr);
		    ((TH1F*)hdata)->SetBinContent(9,excluded);
		    ((TH1F*)hdata)->SetBinContent(10,empty);
		    ((TH1F*)hdata)->SetBinContent(11,excludeddriftregion);
		    ((TH1F*)hdata)->SetBinContent(12,emptydr);
		    ((TH1F*)hdata)->SetBinContent(13,exactive);
		    ((TH1F*)hdata)->SetBinContent(14,emptydiff);
		    ((TH1F*)hdata)->SetBinContent(15,exactivedriftregion);
		    ((TH1F*)hdata)->SetBinContent(16,emptydr);
		    
		    //layer 3
		    ((TH1F*)hdata)->SetBinContent(19,activemoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(20,filledmodules[0]);
		    ((TH1F*)hdata)->SetBinContent(21,activedrperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(22,filleddriftregion[0]);
		    ((TH1F*)hdata)->SetBinContent(23,excludedmoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(24,emptymodules[0]);
		    ((TH1F*)hdata)->SetBinContent(25,excludeddrperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(26,emptydriftregion[0]);
		    ((TH1F*)hdata)->SetBinContent(27,exactivemoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(28,emptyactivemoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(29,activedrperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(30,emptyactivedrperlayer[0]);
		    
		    //layer 4
		    ((TH1F*)hdata)->SetBinContent(33,activemoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(34,filledmodules[1]);
		    ((TH1F*)hdata)->SetBinContent(35,activedrperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(36,filleddriftregion[1]);
		    ((TH1F*)hdata)->SetBinContent(37,excludedmoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(38,emptymodules[1]);
		    ((TH1F*)hdata)->SetBinContent(39,excludeddrperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(40,emptydriftregion[1]);
		    ((TH1F*)hdata)->SetBinContent(41,exactivemoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(42,emptyactivemoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(43,activedrperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(44,emptyactivedrperlayer[1]);
		    //break; 
		  }
	    }//if hdata
	    
	  }//end while 
	
      }//end else 
      delete hmodule;
      hmodule=NULL;
      for(Int_t i=0;i<2;i++) {
	delete hlayer[i];
	hlayer[i]=NULL;
      }//end for

    }//end raw
      
      break;
      
  case AliQAv1::kNULLTASK:{
    AliInfo(Form("No Check on %s\n",AliQAv1::GetAliTaskName(index))); 
    sddQACheckerValue=1.;
  }
    break;
    
  case AliQAv1::kREC:
    {
      Int_t uidrec=list->GetUniqueID();
      AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index))); 
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
      Int_t emptymodules[2], filledmodules[2],emptyladders[2],filledladders[2],emptydriftregion[2], filleddriftregion[2], excludedmoduleperlayer[2], excludeddrperlayer[2], activemoduleperlayer[2],activedrperlayer[2],exactivemoduleperlayer[2],exactivedrperlayer[2];
      Int_t excluded=0;            //excluded modules  
      Int_t excludeddriftregion=0; //excluded single drift region
      Int_t active=0;              //active modules  
      Int_t activedriftregion=0;   //active single drift region
      Int_t exactive=0;            //excluded modules but taking data
      Int_t exactivedriftregion=0; //excluded single drift region but taking data   
      Int_t empty=0;
      Int_t filled=0;
      Int_t emptydr=0;
      Int_t filleddr=0;
      //Int_t emptyactivemodule=0;
      Int_t emptyactivemoduleperlayer[2];
      //Int_t emptydractivemodule=0;
      Int_t emptyactivedrperlayer[2];
      Int_t emptysum=0;
      Int_t emptydiff=0;
      Int_t emptydrsum=0;
      Int_t emptydrdiff=0;

    for(Int_t i=0;i<2;i++)
      {
	emptymodules[i]=0; 
	filledmodules[i]=0; 
	emptyladders[i]=0; 
	filledladders[i]=0; 
	emptydriftregion[i]=0;
	filleddriftregion[i]=0; 
	excludedmoduleperlayer[i]=0;
	excludeddrperlayer[i]=0;
	activemoduleperlayer[i]=0;
	activedrperlayer[i]=0;
	exactivemoduleperlayer[i]=0;
	exactivedrperlayer[i]=0;
	emptyactivemoduleperlayer[i]=0;
	emptyactivedrperlayer[i]=0;
      }  

    Int_t neventsrecpoints=0;

      while( (hdata = dynamic_cast<TH1* >(next())) ){
	if (hdata){
	  TString hname=hdata->GetName();

	  if(hname.Contains("SDDRecPointCheck"))
	    { 
	      neventsrecpoints=(Int_t)hdata->GetBinContent(1); 
	      //break;
	    }
	  else{continue;}
	}//end if hdata
      }//end while

      next.Begin();

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
	      Float_t fractionAboveThreshold=0.;
	      if(hdata->GetEntries()>0.) fractionAboveThreshold = ((Float_t) aboveThreshold)/hdata->GetEntries();
	      if(hname.Contains("L3")) AliInfo(Form("SDD check number 3, L3: RecPoints fractionAboveThreshold: %f",fractionAboveThreshold));
	      if(hname.Contains("L4")) AliInfo(Form("SDD check number 4, L4: RecPoints fractionAboveThreshold: %f",fractionAboveThreshold));
	      if(fractionAboveThreshold > fThresholdForRelativeOccupancy) { 
		sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
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
		sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
		if(hname.Contains("L3")) AliInfo(Form("SDD check number 5: Set Warning (L3 RecPoints2Raws)"));
		if(hname.Contains("L4")) AliInfo(Form("SDD check number 6: Set Warning (L4 RecPoints2Raws)"));
	      }
	    }
	    if(hname.Contains("dedx")) {
	      if(hname.Contains("L3")) AliInfo(Form("SDD check number 7: L3 average charge: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	      if(hname.Contains("L4")) AliInfo(Form("SDD check number 8: L4 average charge: %f, rms: ,%f",hdata->GetMean(),hdata->GetRMS()));
	    }//end if dedx
	  if(hname.Contains("SDDModPatternRP")){
	    if(hname.Contains("NORM")) continue;
	    hmodule=(TH1*)hdata->Clone();
	    entries= hdata->GetEntries();
	    if(AliITSQADataMakerRec::AreEqual(entries,0.)){AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];}//endif entries
	    else{
	      int modmax=hdata->GetNbinsX();
	      Double_t content=0;
	      Int_t llay=0;
	      for(Int_t i=1;i<=modmax;i++)
		{
		  if(i<85)llay=0;
		  else llay=1;
		  content=hdata->GetBinContent(i);
		  if(AliITSQADataMakerRec::AreEqual(content,0.))
		    {
		    empty++;
		    emptymodules[llay]++;
		    } 
		  else 
		    {
		      filled++;
		      filledmodules[llay]++;
		    }
		}//end for
	      //output of the check at the level of the modules. Drift region in the following checks

	      AliInfo(Form(" %s : empty modules %i \t filled modules %i",hname.Data(), empty, filled));
	      AliInfo(Form(" %s : Layer 3 empty modules %i \t filled modules %i",hname.Data(), emptymodules[0], filledmodules[0]));
	      AliInfo(Form(" %s : Layer 4 empty modules %i \t filled modules %i",hname.Data(), emptymodules[1], filledmodules[1]));
	    }//end else pattern entries !=0
	  } //modpattern (1d histogram) 		

	  if(hname.Contains("SDDModPatternL3RP") || hname.Contains("SDDModPatternL4RP")){if(hname.Contains("NORM"))continue;
	    Int_t layer=0;
	    if(hname.Contains("3"))layer=0;
	    else  if(hname.Contains("4"))layer=1;
	    entries2[layer]=hdata->GetEntries();
	    if(entries2[layer]==0){AliWarning(Form("===================>>>>>> No entries in  %s \n",hname.Data()));
	      sddQACheckerValue += fStepBitSDD[AliQAv1::kFATAL];}//end if getentries
	    else{
	      Int_t layer1=0;
	      if(hname.Contains("3"))layer1=0;
	      else  if(hname.Contains("4"))layer1=1;
	      TH2* htemp=dynamic_cast<TH2*>(hdata);
	      if(htemp){
		hlayer[layer1]=(TH2*)htemp->Clone();
		hlayer[layer1]->SetName(Form("%s_copy",hname.Data()));
		//hlayer[layer1]->RebinX(2);
		int modmay=hlayer[layer1]->GetNbinsY();
		TH1D* hproj= hlayer[layer1]->ProjectionY();
		Double_t ladcontent=0;
		for(Int_t i=1;i<=modmay;i++) {//loop on the ladders
		  ladcontent=hproj->GetBinContent(i);
		  if(AliITSQADataMakerRec::AreEqual(ladcontent,0.)) emptyladders[layer1]++;
		  else filledladders[layer1]++;}//end for
		AliInfo(Form(" %s : empty ladders %i \t filled ladders %i\n",hname.Data(), emptyladders[layer1], filledladders[layer1]));//end else layer 3
		delete hproj;
		hproj=NULL;
	      }//end if htemp
	    }//end else entries !=0
	  }//end check on phiz
	  }//end if hdata
	}//end while				
      for(Int_t ii=0;ii<2;ii++)
	{
	  filledmodules[ii]=0;
	  emptymodules[ii]=0;
	}
      filled=0;
      empty=0;
      if(AliITSQADataMakerRec::AreEqual(entries,0.)&& AliITSQADataMakerRec::AreEqual(entries2[0],0.)&& AliITSQADataMakerRec::AreEqual(entries2[1],0.)) break;
      //else{
      if(hmodule || (hlayer[0] && hlayer[1])){
	for(Int_t imod=0;imod<fgknSDDmodules;imod++){
	  Int_t lay=0;
	  Int_t lad=0;
	  Int_t det=0;
	  Int_t module=0;
	  module=imod+fgkmodoffset;
	  AliITSCalibrationSDD * cal=(AliITSCalibrationSDD*)fCalibration->At(imod);
	  if(cal==0) { continue;}
	  AliITSgeomTGeo::GetModuleId(module,lay,lad,det);
	  if (cal->IsBad()){
	    excluded++;
	    excludedmoduleperlayer[lay-3]++;
	    Double_t content=0.;
	    Double_t contentlayer[2];
	    for(Int_t i=0;i<2;i++)contentlayer[i]=0.;
	    if(hmodule)content=hmodule->GetBinContent(imod+1);//if expert bit is active the histogram will be stored in the QA file otherwise the histogram will not be written on the logbook
	    //	    if(hlayer[lay-3]) contentlayer[lay-3]=hlayer[lay-3]->GetBinContent(det,lad);
	    if(AliITSQADataMakerRec::AreEqual(content,0.)== kFALSE) {
	      filledmodules[lay-3]++;
	      filled++;
	      AliError(Form("The module %d (layer %i, ladder %i det %i ) excluded from the acquisition,has recpoints \n ",module,lay,lad,det));
	      exactive++;
	      exactivemoduleperlayer[lay-3]++;
	    } else if(AliITSQADataMakerRec::AreEqual(content,0.)) 
	      {
		emptymodules[lay-3]++; //it has to be empty
		empty++;
	      }
	  } else {
	    Int_t totside=0;
	    Int_t totactiveside=0;
	    //Int_t totbadside=0;
	    Double_t contentgood=0.;

	    for(Int_t i=0;i<2;i++){
	      if(hlayer[lay-3]) contentgood=hlayer[lay-3]->GetBinContent(2*det+i-1,lad);
	      if(cal->IsWingBad(i))
		{ 
		  excludeddriftregion++; 
		  excludeddrperlayer[lay-3]++; 
		  if(AliITSQADataMakerRec::AreEqual(contentgood,0.)==kFALSE){
		    AliError(Form("The side %d of the module %d (layer %i, ladder %i det %i ) excluded from the acquisition, has recpoints \n ",i,module,lay,lad,det));
		    exactivedriftregion++; 
		    exactivedrperlayer[lay-3]++;
		    filleddr++;
		    filleddriftregion[lay-3]++;
		  }
		}//end wingbad
	      else{
		if(AliITSQADataMakerRec::AreEqual(contentgood,0.)==kTRUE)
		  {
		    AliWarning(Form("The side %d of the module %d (layer %i, ladder %i det %i ) is in acquisition, but no recpoints are present\n ",i,module, lay, lad, det));
		  }
		else 
		  {
		    totside++;
		  }
		totactiveside++;
	      }
	    }//end for
	    if(totside==0){
	    AliWarning(Form("The  module %d (layer %i, ladder %i det %i ) is in acquisition, but no recpoints are present \n ",module, lay, lad, det));
	      emptymodules[lay-3]++;
	      empty++;

	    }
	      else 
		if(totside==2){
		filledmodules[lay-3]++;
		filled++;
		}
		else
		  if(totside==1)
		    {
		      //		      emptydr++;
		      //emptydriftregion[lay-3]++; //it has to be empty
		      emptydriftregion[lay-3]++; 
		      emptydr++;
		      filleddr++;
		      filleddriftregion[lay-3]++;
		    }
	    if(totactiveside==1)
	      {
		activedriftregion++;
		activedrperlayer[lay-3]++;
	      }else if(totactiveside==2)
	      {
		active++;
		activemoduleperlayer[lay-3]++;
	      }

	  }
	}//end for
	AliInfo(Form("In total %d modules and %d single drift regions have recpoints.\n ",filled, filleddr));
	AliInfo(Form("In total %d modules and %d single drift regions are empty\n",empty, emptydr));
	for(Int_t i=0;i<2;i++)
	  {
	    AliInfo(Form("Layer %i   \tempty modules %i             \t filled modules %i\n", i+3,emptymodules[i], filledmodules[i]));
	    AliInfo(Form("Layer %i   \tempty single drift regions %i \t filled single drift regions %i\n",i+3,emptydriftregion[i], filleddriftregion[i]));
	  }//end else layers
	emptysum=emptymodules[0]+emptymodules[1];
	emptydiff=emptysum-excluded;
	emptyactivemoduleperlayer[0]=emptymodules[0]- excludedmoduleperlayer[0];
	emptyactivemoduleperlayer[1]=emptymodules[1]- excludedmoduleperlayer[1];

	emptydrsum=emptydriftregion[0]+emptydriftregion[1];
	emptydrdiff=emptydrsum-excludeddriftregion;
	emptyactivedrperlayer[0]=emptydriftregion[0]- excludeddrperlayer[0];
	emptyactivedrperlayer[1]=emptydriftregion[1]- excludeddrperlayer[1];

	Int_t numlimit=1000;
 	if(emptysum>excluded||emptydrsum>excludeddriftregion){ 
	  AliWarning(Form(" %i good module(s) and %i good single drift regions have not recpoints! \n",emptydiff,emptydrdiff));
	  AliWarning(Form(" Layer 3: %i good module(s) and %i good single drift regions have not recpoints! \n",emptyactivemoduleperlayer[0] ,emptyactivedrperlayer[0] ));
	  AliWarning(Form(" Layer 4: %i good module(s) and %i good single drift regions have not recpoints! \n",emptyactivemoduleperlayer[1] ,emptyactivedrperlayer[1] ));
	  
	  //sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
	  Int_t numlimits=1000;	  

	  results1.Form("%i good module(s) and %i good drift region(s) have not recpoints!",emptydiff,emptydrdiff);
	  //	  if((AliQAv1::Instance(AliQAv1::GetDetIndex("ITS")))->IsEventSpecieSet(AliRecoParam::kCosmic)==kFALSE)
	  if(AliRecoParam::Convert(GetEventSpecieForCheck())!=AliRecoParam::kCosmic){     
	    
	      if(neventsrecpoints<numlimits)
		{
		  results2.Form(" Events %d .Too few events.DO NOT CALL the Expert ",neventsrecpoints);
		  color=5;
		  sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
		}
	      else
		{
		  results2.Form(" Events %d .If PHYISICS, follow the TWiki instruction and call the Expert ",neventsrecpoints);
		  color=2;
		  sddQACheckerValue=fHighSDDValue[AliQAv1::kERROR];
		}
	    }
	  else
	    //if((AliQAv1::Instance(AliQAv1::GetDetIndex("ITS")))->IsEventSpecieSet(AliRecoParam::kCosmic)==kTRUE)
	    if(AliRecoParam::Convert(GetEventSpecieForCheck())==AliRecoParam::kCosmic){     
		numlimit=10000;
		if( neventsrecpoints<numlimit)
		  {
		    AliWarning(Form("This is a cosmic run. Some module and drift region are empty but all is OK. "));
		    results1.Form("OK. Thi is a cosmic run. You need a lot of events. Events %i",neventsrecpoints);
		    results2.Form("%i good module(s) and %i good drift are empty! DO NOT CALL THE EXPERT",emptydiff,emptydrdiff);
		    color=5;
		    sddQACheckerValue=fHighSDDValue[AliQAv1::kWARNING];
		  }
		else
		  {
		    results1.Form("%i good module(s) and %i good drift region(s) have not recpoints!",emptydiff,emptydrdiff);
		    results2.Form("Cosmic Events %d .Follow the TWiki instruction and call the Expert ",neventsrecpoints);
		    color=2;
		    sddQACheckerValue=fHighSDDValue[AliQAv1::kERROR];
		  }
	      }
	}
	
      

	if(exactive==0 && emptydiff==0 && exactivedriftregion==0 && emptydrdiff==0){
	  AliInfo(Form("All the active modules (%i) and single drift regions (%i) are in acquisition. The number of excluded modules are %i and the excluded single drift regions are %i\n",active,activedriftregion,excluded,excludeddriftregion));
	  results1.Form("OK.");
	  results2.Form(" All active modules have recpoints");
	  color=3;
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kINFO];
	}
	if(exactive!=0||exactivedriftregion!=0){
	  AliError(Form("%i modules and %i single  drift regions excluded from the acquisition have recpoints. Active modules %i single drift region %i \n ",exactive,exactivedriftregion,active,activedriftregion));
	  AliError(Form("Layer 3: %i modules and %i single  drift regions excluded from the acquisition have recpoints. Active modules %i single drift region %i \n ",exactivemoduleperlayer[0],exactivedrperlayer[0],activemoduleperlayer[0],activedrperlayer[0]));
	  AliError(Form("Layer 3: %i modules and %i single  drift regions excluded from the acquisition have recpoints. Active modules %i single drift region %i \n ",exactivemoduleperlayer[1],exactivedrperlayer[1],activemoduleperlayer[1],activedrperlayer[1]));
	  results1.Form("%i modules and %i drift region excluded from the acquisition have recpoints",exactive,exactivedriftregion);
	  results2.Form("Follow the TWiki instructions and Call the SDD expert ");
	  color=2;	
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kERROR];
	}
	if(excluded==exactive||excludeddriftregion==exactivedriftregion){
	  AliError(Form("All the modules (%d) or single drift regions (%d) excluded from the acquisition have recpoints!\n  Active modules %i \t Active drfift regions %i\n",excluded,excludeddriftregion,active,activedriftregion));
	  results1.Form("All the modules (%d) or drift regions (%d) excluded from the acquisition have recpoints!",excluded,excludeddriftregion );
	  results2.Form("Follow the TWiki instructions and Call the SDD expert ");
	  color=6;
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	}
	if(active==0||activedriftregion==0){
	  AliError(Form("No modules or single drift regions have recpoints: excluded %i \t excluded active %i \n\t\t excluded single drift regions %i \t excluded active drift regions %i \n", excluded, exactive, excludeddriftregion, exactivedriftregion)); 
	  results1.Form("No modules or drift region have recpoints: excluded modules %i  excluded drift regions %i ", excluded, excludeddriftregion );
	  results2.Form("Follow the TWiki instructions and Call the SDD expert ");
	  color=6;
	  sddQACheckerValue=fHighSDDValue[AliQAv1::kFATAL];
	}

	TPaveText *pave[2];
	next.Begin();

	while( (hdata=dynamic_cast<TH1* >(next())) )
	  {
	    if (hdata){
	      TString hname=hdata->GetName();
	      if(hname.Contains("SDDModPatternL3RP") || hname.Contains("SDDModPatternL4RP")){
		if(hname.Contains("NORM"))continue;
		//AliInfo("========================================Found histo 11\n");
		Int_t lay=0;
		if(hname.Contains("3"))lay=0;
		else if(hname.Contains("4"))lay=1;
		pave[lay]=new TPaveText(0.3,0.88,0.9,0.99,"NDC");
		pave[lay]->AddText(results1.Data());
		pave[lay]->AddText(results2.Data());
		pave[lay]->SetFillColor(color);
		pave[lay]->SetBorderSize(1);
		pave[lay]->SetLineWidth(1);
		hdata->GetListOfFunctions()->Add(pave[lay]);
	      }
	      else
		if(hname.Contains("SDDRecPointCheck"))
		  {
		    
		    //AliInfo("========================================Found histo\n");
		    ((TH1F*)hdata)->SetBinContent(5,active);
		    ((TH1F*)hdata)->SetBinContent(6,filled);
		    ((TH1F*)hdata)->SetBinContent(7,activedriftregion);
		    ((TH1F*)hdata)->SetBinContent(8,filleddr);
		    ((TH1F*)hdata)->SetBinContent(9,excluded);
		    ((TH1F*)hdata)->SetBinContent(10,empty);
		    ((TH1F*)hdata)->SetBinContent(11,excludeddriftregion);
		    ((TH1F*)hdata)->SetBinContent(12,emptydr);
		    ((TH1F*)hdata)->SetBinContent(13,exactive);
		    ((TH1F*)hdata)->SetBinContent(14,emptydiff);
		    ((TH1F*)hdata)->SetBinContent(15,exactivedriftregion);
		    ((TH1F*)hdata)->SetBinContent(16,emptydr);
		    
		    //layer 3
		    ((TH1F*)hdata)->SetBinContent(19,activemoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(20,filledmodules[0]);
		    ((TH1F*)hdata)->SetBinContent(21,activedrperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(22,filleddriftregion[0]);
		    ((TH1F*)hdata)->SetBinContent(23,excludedmoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(24,emptymodules[0]);
		    ((TH1F*)hdata)->SetBinContent(25,excludeddrperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(26,emptydriftregion[0]);
		    ((TH1F*)hdata)->SetBinContent(27,exactivemoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(28,emptyactivemoduleperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(29,activedrperlayer[0]);
		    ((TH1F*)hdata)->SetBinContent(30,emptyactivedrperlayer[0]);
		    
		    //layer 4
		    ((TH1F*)hdata)->SetBinContent(35,activemoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(36,filledmodules[1]);
		    ((TH1F*)hdata)->SetBinContent(37,activedrperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(38,filleddriftregion[1]);
		    ((TH1F*)hdata)->SetBinContent(39,excludedmoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(40,emptymodules[1]);
		    ((TH1F*)hdata)->SetBinContent(41,excludeddrperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(42,emptydriftregion[1]);
		    ((TH1F*)hdata)->SetBinContent(43,exactivemoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(44,emptyactivemoduleperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(45,activedrperlayer[1]);
		    ((TH1F*)hdata)->SetBinContent(46,emptyactivedrperlayer[1]);
		    
		  }
	    }//if hadata
	    
	  }//end while 
	
      }//end else 
      delete hmodule;
      hmodule=NULL;
      for(Int_t i=0;i<2;i++) {
	delete hlayer[i];
	hlayer[i]=NULL;
      }//end for      


      //sddQACheckerValue=1.;
      }//end recpoint list uid = 20
      if(uidrec==40)
	{
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
      AliInfo(Form("===================> No Check on %s\n",AliQAv1::GetAliTaskName(index)));
      sddQACheckerValue=1.; 
    }
    break;
  case AliQAv1::kESD:
    {
      AliInfo(Form("==================>  No Check on %s\n",AliQAv1::GetAliTaskName(index)));
      sddQACheckerValue=1.;
    } 
    break;
  case AliQAv1::kNTASK:{
    AliInfo(Form("==================>  No Check on %s\n",AliQAv1::GetAliTaskName(index))); 
    sddQACheckerValue=1.;
  }
    break;
  case AliQAv1::kSIM:{
    AliInfo(Form("Check on %s\n",AliQAv1::GetAliTaskName(index))); 
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
  if(hdata) delete hdata;


  return sddQACheckerValue;	
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
      rval=MakeSDDRawsImage(list, task,mode);
    }
      break;
    case AliQAv1::kRECPOINTS:{ rval=MakeSDDRecPointsImage(list, task,mode); }
      break;
    case AliQAv1::kHITS:; case AliQAv1::kESDS:; case AliQAv1::kDIGITS:;case AliQAv1::kDIGITSR:;case AliQAv1::kSDIGITS:;case AliQAv1::kTRACKSEGMENTS:;case AliQAv1::kRECPARTICLES:; default:
    {
       rval=kFALSE;
       //AliQAChecker::Instance()->GetDetQAChecker(0)->MakeImage(list,task,mode);
    }
    break;
    case AliQAv1::kNULLTASKINDEX:; case  AliQAv1::kNTASKINDEX: 
      {AliWarning(Form("No histograms for this task ( %s ) \n", AliQAv1::GetTaskName(task).Data())); rval=kFALSE;}
      break;
    }
return rval;  
}


//_______________________________________________________________________
Bool_t AliITSQASDDChecker::MakeSDDRawsImage(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode )
{
  // MakeSDDRawsImage: raw data QA plots

    for (Int_t esIndex = 0 ; esIndex < AliRecoParam::kNSpecies ; esIndex++) {
      //printf("-------------------------> %i \n", esIndex);
      if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(esIndex)) || list[esIndex]->GetEntries() == 0) 
          {//printf ("Nothing for %s \n", AliRecoParam::GetEventSpecieName(esIndex));
	 continue;
	}
      else{
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
	Int_t nx =2; //TMath::Nint(TMath::Sqrt(nImages));
	Int_t ny =2; // nx  ; 
	//if (nx < TMath::Sqrt(nImages))
	//ny++ ;  
	fImage[esIndex]->Divide(nx, ny) ; 
	TIter nexthist(list[esIndex]) ; 
	TH1* hist = NULL ;
	Int_t npad = 1 ; 
	fImage[esIndex]->cd(npad); 
	fImage[esIndex]->cd(npad)->SetBorderMode(0) ;
	while ( (hist=static_cast<TH1*>(nexthist())) ) {
	  //gPad=fImage[esIndex]->cd(npad)->GetPad(npad);
	  TString cln(hist->ClassName()) ; 
	  if ( ! cln.Contains("TH") )
	    continue ;
	  
	  if(hist->TestBit(AliQAv1::GetImageBit())) {
	    hist->GetXaxis()->SetTitleSize(0.02);
	    hist->GetYaxis()->SetTitleSize(0.02);
	    hist->GetXaxis()->SetLabelSize(0.02);
	    hist->GetYaxis()->SetLabelSize(0.02);
	    if(cln.Contains("TH2"))
	      {
		gPad->SetRightMargin(0.15);
		gPad->SetLeftMargin(0.05);
		hist->SetStats(0);
		hist->SetOption("colz") ;
		//hist->GetListOfFunctions()->FindObject("palette")->SetLabelSize(0.025);
		//gPad->Update();
	      }
	    hist->DrawCopy() ; 
	    fImage[esIndex]->cd(++npad) ; 
	    fImage[esIndex]->cd(npad)->SetBorderMode(0) ; 
	  }
	}
	fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
      }
    }
   return kTRUE;
}




//_______________________________________________________________________
Bool_t AliITSQASDDChecker::MakeSDDRecPointsImage(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode )
{
  // MakeSDDRecPointsImage: rec point QA plots

    for (Int_t esIndex = 0 ; esIndex < AliRecoParam::kNSpecies ; esIndex++) {
      if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(esIndex)) || list[esIndex]->GetEntries() == 0) 
        {
	//printf ("Nothing for %s \n", AliQAv1::GetTaskName(task).Data()); 
	continue;
	}
      const Char_t * title = Form("QA_%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(esIndex)) ; 
      if ( !fImage[esIndex] ) {
        fImage[esIndex] = new TCanvas(title, title,1280,980) ;
      }
      fImage[esIndex]->Clear() ; 
      fImage[esIndex]->SetTitle(title) ; 
      fImage[esIndex]->cd();
      fImage[esIndex]->SetBorderMode(0) ;  
      TPaveText someText(0.015, 0.015, 0.98, 0.98);
      someText.AddText(title);
      someText.Draw(); 
      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
      fImage[esIndex]->Clear() ; 
      Int_t nx =2; //TMath::Nint(TMath::Sqrt(nImages));
      Int_t ny =4; // nx  ; 
      //if (nx < TMath::Sqrt(nImages))
      //ny++ ;  
      fImage[esIndex]->Divide(nx, ny) ; 
      TIter nexthist(list[esIndex]) ; 
      TH1* hist = NULL ;
      Int_t npad = 1 ; 
      fImage[esIndex]->cd(npad) ; 
      fImage[esIndex]->cd(npad)->SetBorderMode(0) ; 
      while ( (hist=static_cast<TH1*>(nexthist())) ) {
	//gPad=fImage[esIndex]->cd(npad)->GetPad(npad);
        TString cln(hist->ClassName()) ;
	//printf("=====================> Class name %s \n",cln.Data()); 
        if ( ! cln.Contains("TH") )
          continue ;
        if(hist->TestBit(AliQAv1::GetImageBit())) {
	    hist->GetXaxis()->SetTitleSize(0.02);
	    hist->GetYaxis()->SetTitleSize(0.02);
	    hist->GetXaxis()->SetLabelSize(0.02);
	    hist->GetYaxis()->SetLabelSize(0.02);
	  if(cln.Contains("TH1"))
	    {
	      hist->SetFillColor(kOrange+7);
	      //SetFrameFillColor(kAzure-9);
	      //hist->DrawCopy() ; 
	    }
	  if(cln.Contains("TH2"))
	    {
	      gPad->SetRightMargin(0.15);
	      gPad->SetLeftMargin(0.05);
	      hist->SetStats(0);
	      hist->SetOption("colz") ;
	      //	      TPaletteAxis *paletta =(TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
	      //paletta->SetLabelSize(0.025);
	      //gPad->Update(); 
	    }
	  hist->DrawCopy();
          fImage[esIndex]->cd(++npad) ; 
	  fImage[esIndex]->cd(npad)->SetBorderMode(0) ; 
        }
      }
      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
    }
    // }  
   return kTRUE;
}
