/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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


///////////////////////////////////////////////////////////////////////////////
// CPV Preprocessor class. It runs by Shuttle at the end of the run,
// calculates pedestals, calibration coefficients and dead/bad channels
// to be posted in OCDB.
//
// Author: Boris Polichtchouk, 25 January 2008
// Updated:Sergey Evdokimov, 28 Mar 2015
///////////////////////////////////////////////////////////////////////////////

#include "AliPHOSCpvPreprocessor.h"
#include "AliLog.h"
#include "AliCDBMetaData.h"
#include "AliPHOSCpvBadChannelsMap.h"
#include "AliPHOSCpvCalibData.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMap.h"
#include "TKey.h"
#include "TList.h"
#include "TObjString.h"
#include "AliCDBEntry.h"
ClassImp(AliPHOSCpvPreprocessor)

//_______________________________________________________________________________________
AliPHOSCpvPreprocessor::AliPHOSCpvPreprocessor() :
AliPreprocessor("CPV",0)
{
  //default constructor
}

//_______________________________________________________________________________________
AliPHOSCpvPreprocessor::AliPHOSCpvPreprocessor(AliShuttleInterface* shuttle):
AliPreprocessor("CPV",shuttle)
{
  // Constructor

  AddRunType("PHYSICS");
  AddRunType("PEDESTAL");

}

//_______________________________________________________________________________________
UInt_t AliPHOSCpvPreprocessor::Process(TMap* /*valueSet*/)
{
  // process data retrieved by the Shuttle

  // The fileName with the histograms which have been produced by
  // CPVPEDda.cxx, CPVBCMda.cxx or CPVGAINda.cxx.
  // It is a responsibility of the SHUTTLE framework to form the fileName.
  
  AliCDBMetaData md;
  md.SetResponsible("Sergey Evdokimov");
  
  TString runType = GetRunType();
  Log(Form("Run type: %s",runType.Data()));

  if(runType=="PHYSICS") {
    Bool_t BCM = false, GAIN = false;
    Bool_t storeOK_BCM=false, storeOK_GAIN=false;
    AliPHOSCpvCalibData* calib=0;
    AliPHOSCpvBadChannelsMap* badMap=0;
    //===================================================
    //=============== Bad Channels Map ==================
    //===================================================
    TList* list = GetFileSources(kDAQ, "CPVBADMAP");
    if(!list) {
      Log("CPVBADMAP sources list not found, moving to gain calibration.");
    }else BCM=true;
    
    if(BCM)
      {
	TIter iter(list);
	TObjString *source;
	badMap = new AliPHOSCpvBadChannelsMap();
	while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
	  AliInfo(Form("found source %s", source->String().Data()));
	  
	  TString fileName = GetFile(kDAQ, "CPVBADMAP", source->GetName());
	  AliInfo(Form("Got filename: %s",fileName.Data()));
	  
	  TFile f(fileName);
	  
	  if(!f.IsOpen()) {
	    Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
	    return 1;
	  }
	  for(Int_t iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL;iDDL+=2)
	    if(f.Get(Form("hBadChMap%d",iDDL))){
	      badMap->Reset(AliPHOSCpvParam::DDL2Mod(iDDL));
	      TH2* hBadMap = (TH2*)f.Get(Form("hBadChMap%d",iDDL));
	      for(Int_t iX = 0; iX<AliPHOSCpvParam::kPadPcX;iX++)
		for(Int_t iY = 0; iY<AliPHOSCpvParam::kPadPcY;iY++)
		  if(hBadMap->GetBinContent(iX+1,iY+1))
		    badMap->SetBadChannel(AliPHOSCpvParam::DDL2Mod(iDDL),iX+1,iY+1);
	    }
	  f.Close();
	}//while((source = dynamic_cast<TObjString *> (iter.Next())))
	Bool_t storeOK_BCM = Store("Calib", "CpvBadChannels", badMap, &md, 0, kTRUE);
      }//if(BCM)
    //===================================================
    //=============== GAIN calibration ==================
    //===================================================
    list = GetFileSources(kDAQ, "CPVAMPLITUDES");
    if(!list) {
      Log("Sources list not found, exit.");
    }else GAIN=true;

    if(GAIN){
      //Retrieve the last Cpv calibration & BadChMap objects
      AliCDBEntry* entryCalib = GetFromOCDB("Calib", "CpvGainPedestals");
      if(!entryCalib){
	Log(Form("Cannot find any AliCDBEntry for [Calib, CpvGainPedestals]!\nGoing to create new CpvGainPedestals object."));
	calib = new AliPHOSCpvCalibData();
      }
      else
	calib = (AliPHOSCpvCalibData*)entryCalib->GetObject();

      //I think that we don't need bad map because bad channels exclusion is done at CPVGAINda working time.
      // if(!storeOK_BCM){//check if we already prepared BCM in for current run
      // 	AliCDBEntry* entryBadMap = GetFromOCDB("Calib", "CpvBadChannels");
      // 	if(!entryBadMap){
      // 	  Log(Form("Cannot find any AliCDBEntry for [Calib, CpvBadChannels]!"));
      // 	  badMap = new AliPHOSCpvBadChannelsMap();//dummy bad channels map just to simplify following code
      // 	}else
      // 	  badMap = (AliPHOSCpvBadChannelsMap*)entryCalib->GetObject();
      // }
      
      TIter iter(list);
      TObjString *source;
  
      while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
	AliInfo(Form("found source %s", source->String().Data()));

	TString fileName = GetFile(kDAQ, "CPVAMPLITUDES", source->GetName());
	AliInfo(Form("Got filename: %s",fileName.Data()));

	TFile f(fileName);

	if(!f.IsOpen()) {
	  Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
	  return 1;
	}
	Int_t minimalStatistics = 150;// minimal statistics to calculate calibration coeff
	Float_t coeff;
	TF1* fitFunc = new TF1("fitFunc","landau",0.,4000.);
	fitFunc->SetParameters(1.,200.,60.);fitFunc->SetParLimits(1,10.,2000.);
	for(Int_t iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL;iDDL+=2){
	  TH2* entriesMap=(TH2*)f.Get(Form("hEntriesMap%d",iDDL));
	  if(!entriesMap)continue;
	  for(Int_t iX = 0; iX<AliPHOSCpvParam::kPadPcX;iX++)
	    for(Int_t iY = 0; iY<AliPHOSCpvParam::kPadPcY;iY++)
	      if(entriesMap->GetBinContent(iX+1,iY+1)){//we have some statistics for current channel
		fitFunc->SetParameters(1.,200.,60.);
		TH1* hAmpl = (TH1*)f.Get(Form("hAmplA0_DDL%d_iX%d_iY%d",iDDL,iX,iY));
		if(!hAmpl)continue;
		hAmpl->Rebin(4);
		if(hAmpl->Integral(20,2000)<minimalStatistics)continue;
		hAmpl->Fit(fitFunc,"QL0","",20.,2000.);
		coeff = 200./fitFunc->GetParameter(1);
		calib->SetADCchannelCpv(AliPHOSCpvParam::DDL2Mod(iDDL),iX+1,iY+1,coeff);
		hAmpl->Delete();
	      }
	}
	f.Close();
      }//while ((source =...
      
      //Store CPV calibration data
      Bool_t storeOK_GAIN = Store("Calib", "CpvGainPedestals", calib, &md, 0, kTRUE);
    }//if(GAIN)
    if(GAIN) Log("PHYSICS run: I successfully updated gain calibration coeffs!");
    if(BCM) Log("PHYSICS run: I successfully updated bad channel map!");
    if(!(GAIN||BCM)) Log("PHYSICS run: I did nothing this time");
    return 0;
  }//if(runType=="PHYSICS")

  if(runType=="PEDESTAL") {
    Bool_t PED = false;
    Bool_t storeOK_PED=false;
    AliPHOSCpvCalibData* calib=0;
    //===================================================
    //================== PEDESTALS ======================
    //===================================================
    TList* list = GetFileSources(kDAQ, "CPVPEDS");
    if(!list) {
      Log("Sources list not found, exit.");
    }else PED=true;

    if(PED){
      //Retrieve the last Cpv calibration object
      AliCDBEntry* entryCalib = GetFromOCDB("Calib", "CpvGainPedestals");
      if(!entryCalib){
	Log(Form("Cannot find any AliCDBEntry for [Calib, CpvGainPedestals]!\nGoing to create new CpvGainPedestals object."));
	calib = new AliPHOSCpvCalibData();
      }
      else
	calib = (AliPHOSCpvCalibData*)entryCalib->GetObject();
      TIter iter(list);
      TObjString *source;
  
      while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
	AliInfo(Form("found source %s", source->String().Data()));

	TString fileName = GetFile(kDAQ, "CPVPEDS", source->GetName());
	AliInfo(Form("Got filename: %s",fileName.Data()));

	TFile f(fileName);

	if(!f.IsOpen()) {
	  Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
	  return 1;
	}
      
	for(Int_t iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL;iDDL+=2){
	  TH2* pedMap=(TH2*)f.Get(Form("fPedMeanMap%d",iDDL));
	  if(!pedMap)continue;
	  for(Int_t iX = 0; iX<AliPHOSCpvParam::kPadPcX;iX++)
	    for(Int_t iY = 0; iY<AliPHOSCpvParam::kPadPcY;iY++)
	      calib->SetADCpedestalCpv(AliPHOSCpvParam::DDL2Mod(iDDL),iX+1,iY+1,pedMap->GetBinContent(iX+1,iY+1));
	}
      }//while ((source =...	
      //Store CPV calibration data
      Bool_t storeOK_PED = Store("Calib", "CpvGainPedestals", calib, &md, 0, kTRUE);
    }//if(PED)
    if(PED) Log("PEDESTAL run: I successfully updated pedestals!");
    if(!PED) Log("PEDESTAL run: I did nothing this time");
  }//if(runType=="PEDESTAL")
  Log(Form("Unknown or unused run type %s. Do nothing and return OK.",runType.Data()));
  return 0;
  
}

