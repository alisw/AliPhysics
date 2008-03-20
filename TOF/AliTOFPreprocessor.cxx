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

/* $Id$ */

//#include <Riostream.h>
//#include <stdio.h>
//#include <stdlib.h>

#include <TFile.h>
//#include <TH1.h>
//#include <TH1F.h>
//#include <TH1S.h>
#include <TH2S.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TTimeStamp.h>

#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliTOFChannelOnline.h"
#include "AliTOFChannelOnlineStatus.h"
#include "AliTOFDataDCS.h"
#include "AliTOFGeometry.h"
#include "AliTOFPreprocessor.h"
#include "AliTOFFEEReader.h"

//class TF1;
//class AliDCSValue;
//class AliTOFGeometry;

// TOF preprocessor class.
// It takes data from DCS and passes them to the class AliTOFDataDCS, which
// processes them. The result is then written to the CDB.
// analogously, it takes data form DAQ (both at Run level and inclusive - 
// of all the runs - level, processes them, and stores both Reference Data
// and Online Calibration files in the CDB. 
// Processing of Pulser/Noise Run data and of TOF FEE DCS map

// return codes:
// return=0 : all ok
// return=1 : no DCS input data Map
// return=2 : no DCS input data processing
// return=3 : no DCS processed data was stored in Ref Data
// return=4 : no DAQ input for Ref Data
// return=5 : failed to store DAQ Ref Data
// return=6 : failed to retrieve DAQ data for calibration 
// return=7 : problems in processing histos in the input DAQ file 
// return=8 : failed to store Online Delays
// return=9 : failed to store Reference Data for Pulser
// return=10: failed to retrieve Pulser data 
// return=11: failed to store Pulser map in OCDB
// return=12: failed to store Reference Data for Noise
// return=13: failed to retrieve Noise data 
// return=14: failed to store Noise map in OCDB
// return=15: failed to retrieve FEE data from FXS
// return=16: failed to retrieve FEE data from OCDB
// return=17: failed to store FEE data in OCDB
// return=18: failed to store FEE reference data in OCDB

ClassImp(AliTOFPreprocessor)

const Int_t    AliTOFPreprocessor::fgkBinRangeAve = 13;    // number of bins where to calculate the mean 
const Double_t AliTOFPreprocessor::fgkIntegralThr = 100;   // min number of entries to perform computation of delay per channel 
const Double_t AliTOFPreprocessor::fgkThrPar      = 0.013; // parameter used to trigger the calculation of the delay

//_____________________________________________________________________________

AliTOFPreprocessor::AliTOFPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TOF", shuttle),
  fData(0),
  fh2(0),
  fCal(0),
  fCalStatus(0),
  fFEEStatus(0),
  fNChannels(0),
  fStoreRefData(kTRUE),
  fFDRFlag(kTRUE)
{
  // constructor

}

//_____________________________________________________________________________

AliTOFPreprocessor::~AliTOFPreprocessor()
{
  // destructor
  if (fData){
    delete fData;
    fData = 0;
  }
  if (fh2){
    delete fh2;
    fh2 = 0;
  }
  if (fCal){
    //    fCal->Clear();
    delete fCal;
    fCal = 0;
  }
  if (fCalStatus){
    delete fCalStatus;
    fCalStatus = 0;
  }
  if (fFEEStatus){
    delete fFEEStatus;
    fFEEStatus = 0;
  }
}

//______________________________________________________________________________
void AliTOFPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliTOFDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fData = new AliTOFDataDCS(fRun, fStartTime, fEndTime);
	fData->SetFDRFlag(fFDRFlag);
	fh2 = 0x0;
	fNChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
	fCal = new TObjArray(fNChannels);
	fCal->SetOwner();
	for (Int_t ich = 0; ich<fNChannels; ich ++){
	  AliTOFChannelOnline * calChOnline = new AliTOFChannelOnline();
	  fCal->AddAt(calChOnline,ich);
	}
	fCalStatus = new TObjArray(fNChannels);
	fCalStatus->SetOwner();
	for (Int_t ich = 0; ich<fNChannels; ich ++){
	  AliTOFChannelOnlineStatus * calChOnlineStatus = new AliTOFChannelOnlineStatus();
	  fCalStatus->AddAt(calChOnlineStatus,ich);
	}
	fFEEStatus = new TObjArray(fNChannels);
	fFEEStatus->SetOwner();
	for (Int_t ich = 0; ich<fNChannels; ich ++){
	  AliTOFChannelOnlineStatus * calChOnlineStatus = new AliTOFChannelOnlineStatus();
	  fFEEStatus->AddAt(calChOnlineStatus,ich);
	}
}
//_____________________________________________________________________________
Bool_t AliTOFPreprocessor::ProcessDCS(){

  // check whether DCS should be processed or not...

  TString runType = GetRunType();
  Log(Form("RunType %s",runType.Data()));

  if (runType != "PHYSICS"){
    return kFALSE;
  }

  return kTRUE;
}
//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::ProcessDCSDataPoints(TMap* dcsAliasMap)
{
  // Fills data into a AliTOFDataDCS object


  Log("Processing DCS DP");
  TH1::AddDirectory(0);

  Bool_t resultDCSMap=kFALSE;
  Bool_t resultDCSStore=kFALSE;

  // processing DCS

  if (!dcsAliasMap){
    Log("No DCS map found: TOF exiting from Shuttle");
    return 1;// return error Code for DCS input data not found 
  }
  else {
  // The processing of the DCS input data is forwarded to AliTOFDataDCS
    resultDCSMap=fData->ProcessData(*dcsAliasMap);
    if(!resultDCSMap){
      Log("Some problems occurred while processing DCS data, TOF exiting from Shuttle");
      return 2;// return error Code for processed DCS data not stored 
    }
    else{
      AliCDBMetaData metaDataDCS;
      metaDataDCS.SetBeamPeriod(0);
      metaDataDCS.SetResponsible("Chiara Zampolli");
      metaDataDCS.SetComment("This preprocessor fills an AliTOFDataDCS object.");
      AliInfo("Storing DCS Data");
      resultDCSStore = StoreReferenceData("Calib","DCSData",fData, &metaDataDCS);
      if (!resultDCSStore){
	Log("Some problems occurred while storing DCS data results in Reference Data, TOF exiting from Shuttle");
	return 3;// return error Code for processed DCS data not stored 
	         // in reference data
	
      }
    }
  }
  return 0;
}
//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::ProcessOnlineDelays()
{
  // Processing data from DAQ for online calibration 

  Log("Processing DAQ delays");

  TH1::AddDirectory(0);

  Bool_t resultDAQRef=kFALSE;
  Bool_t resultTOFPP=kFALSE;
  // processing DAQ
  
  TFile * daqFile=0x0;
  
  if(fStoreRefData){
    //retrieving data at Run level
    TList* list = GetFileSources(kDAQ, "RUNLevel");
    if (list !=0x0 && list->GetEntries()!=0)
      {
	AliInfo("The following sources produced files with the id RUNLevel");
	list->Print();
	for (Int_t jj=0;jj<list->GetEntries();jj++){
	  TObjString * str = dynamic_cast<TObjString*> (list->At(jj));
	  AliInfo(Form("found source %s", str->String().Data()));
	  // file to be stored run per run
	  TString fileNameRun = GetFile(kDAQ, "RUNLevel", str->GetName());
	  if (fileNameRun.Length()>0){
	    AliInfo(Form("Got the file %s, now we can store the Reference Data for the current Run.", fileNameRun.Data()));
	    daqFile = new TFile(fileNameRun.Data(),"READ");
	    fh2 = (TH2S*) daqFile->Get("htof");
	    AliCDBMetaData metaDataHisto;
	    metaDataHisto.SetBeamPeriod(0);
	    metaDataHisto.SetResponsible("Chiara Zampolli");
	    metaDataHisto.SetComment("This preprocessor stores the array of histos object as Reference Data.");
	    AliInfo("Storing Reference Data");
	    resultDAQRef = StoreReferenceData("Calib","DAQData",fh2, &metaDataHisto);
	    if (!resultDAQRef){
	      Log("some problems occurred::No Reference Data stored, TOF exiting from Shuttle");
	      return 5;//return error code for failure in storing Ref Data 
	    }
	    daqFile->Close();
	    delete daqFile;
	  }
	  
	  else{
	    Log("The input data file from DAQ (run-level) was not found, TOF exiting from Shuttle "); 
	    return 4;//return error code for failure in retrieving Ref Data 
	  }
	}
	delete list;
      }
    else{
      Log("The input data file list from DAQ (run-level) was not found, TOF exiting from Shuttle "); 
      return 4;//return error code for failure in retrieving Ref Data 
    }	
  }


//Total files, with cumulative histos
  
  TList* listTot = GetFileSources(kDAQ, "DELAYS");
  if (listTot !=0x0 && listTot->GetEntries()!=0)
    {
      AliInfo("The following sources produced files with the id DELAYS");
      listTot->Print();
      for (Int_t jj=0;jj<listTot->GetEntries();jj++){
	TObjString * str = dynamic_cast<TObjString*> (listTot->At(jj));
	AliInfo(Form("found source %s", str->String().Data()));

	// file with summed histos, to extract calib params
	TString fileName = GetFile(kDAQ, "DELAYS", str->GetName());
	if (fileName.Length()>0){
	  AliInfo(Form("Got the file %s, now we can extract some values.", fileName.Data()));

	  daqFile = new TFile(fileName.Data(),"READ");
	  if (fh2) delete fh2;
	  fh2 = (TH2S*) daqFile->Get("htoftot");
	  if (!fh2){
	    Log("some problems occurred:: No histo retrieved, TOF exiting from Shuttle");
	    delete daqFile;
	    return 7; //return error code for histograms not existing/junky
	  }
	  else {
	    static const Int_t kSize=fh2->GetNbinsX();
	    static const Int_t kNBins=fh2->GetNbinsY();
	    static const Double_t kXBinmin=fh2->GetYaxis()->GetBinLowEdge(1);
	    if (kSize != fNChannels){
	      Log(" number of bins along x different from number of pads, found only a subset of the histograms, TOF exiting from Shuttle");
	      delete daqFile;
              return 7; //return error code for histograms not existing/junky
	    }
	    Int_t nNotStatistics = 0; // number of channel with not enough statistics
	    for (Int_t ich=0;ich<kSize;ich++){
	      TH1S *h1 = new TH1S("h1","h1",kNBins,kXBinmin-0.5,kNBins*1.+kXBinmin-0.5);
	      for (Int_t ibin=0;ibin<kNBins;ibin++){
		h1->SetBinContent(ibin+1,fh2->GetBinContent(ich+1,ibin+1));
	      }
	      if(h1->Integral()<fgkIntegralThr) {
		nNotStatistics++;
		if (!fFDRFlag)  Log(Form(" Not enough statistics for bin %i, skipping this channel",ich));  // printing message only if not in FDR runs
		delete h1;
		h1=0x0;
		continue;
	      }
	      if (!fFDRFlag) {  // not computing delays if in FDR runs
		Bool_t found=kFALSE; 
		Float_t minContent=h1->Integral()*fgkThrPar; 
		Int_t nbinsX = h1->GetNbinsX();
		Int_t startBin=1;
		for (Int_t j=1; j<=nbinsX; j++){
		  if ((
		       h1->GetBinContent(j) +     
		       h1->GetBinContent(j+1)+
		       h1->GetBinContent(j+2)+ 
		       h1->GetBinContent(j+3))>minContent){
		    found=kTRUE;
		    startBin=j;
		    break;
		  }
		}
		if(!found) AliInfo(Form("WARNING!!! no start of fit found for histo # %i",ich));
		// Now calculate the mean over the interval. 
		Double_t mean = 0;
		Double_t sumw2 = 0;
		Double_t nent = 0;
		for(Int_t k=0;k<fgkBinRangeAve;k++){
		  mean=mean+h1->GetBinCenter(startBin+k)*h1->GetBinContent(startBin+k);                 
		  nent=nent+h1->GetBinContent(startBin+k);                 
		  sumw2=sumw2+(h1->GetBinCenter(startBin+k))*(h1->GetBinCenter(startBin+k))*(h1->GetBinContent(startBin+k));
		}
		mean= mean/nent; //<x>
		sumw2=sumw2/nent; //<x^2>
		Double_t rmsmean= 0;
		rmsmean = TMath::Sqrt((sumw2-mean*mean)/nent);
		if (ich<fNChannels) {
		  AliTOFChannelOnline * ch = (AliTOFChannelOnline *)fCal->At(ich);
		  ch->SetDelay((Double_t)mean*AliTOFGeometry::TdcBinWidth()*1.E-3);  // delay in ns
		  // ch->SetStatus(1);  // calibrated channel, removed for the time being from AliTOFChannelOnline
		}
	      }
	    delete h1;
	    h1=0x0;
	    }
	    if (nNotStatistics!=0) Log(Form("Too little statistics for %d channels!",nNotStatistics)); 
	  }
	  daqFile->Close();
	  delete daqFile;
	}
	else{
	  Log("The Cumulative data file from DAQ does not exist, TOF exiting from Shuttle"); 
          return 6;//return error code for problems in retrieving DAQ data 
	}
      }
      delete listTot;
    }
  else{
    Log("Problem: no list for Cumulative data file from DAQ was found, TOF exiting from Shuttle");
    return 6; //return error code for problems in retrieving DAQ data 
  }

  daqFile=0;
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Chiara Zampolli");
  metaData.SetComment("This preprocessor fills a TObjArray object.");
  AliInfo("Storing Calibration Data");
  resultTOFPP = Store("Calib","ParOnline",fCal, &metaData,0,kTRUE);
  if(!resultTOFPP){
    Log("Some problems occurred while storing online object resulting from DAQ data processing");
    return 8;//return error code for problems in storing DAQ data 
  }

  return 0;
}
//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::ProcessPulserData()
{
  // Processing Pulser Run data for TOF channel status

  Log("Processing Pulser");

  TH1::AddDirectory(0);

  Bool_t resultPulserRef=kFALSE;
  Bool_t resultPulser=kFALSE;

  static const Int_t kSize = AliTOFGeometry::NPadXSector()*AliTOFGeometry::NSectors();
  TH1S * htofPulser = new TH1S("hTOFpulser","histo with signals on TOF during pulser", kSize,-0.5,kSize-0.5);
  for (Int_t ibin =1;ibin<=kSize;ibin++){
    htofPulser->SetBinContent(ibin,-1);
  }

  // retrieving last stored pulser object, and copying

  AliCDBEntry *cdbEntry = GetFromOCDB("Calib","Pulser");
  if (cdbEntry!=0x0){
    TObjArray *currentCalStatus = (TObjArray*)cdbEntry->GetObject();
    for (Int_t ich = 0; ich<fNChannels; ich ++){
      AliTOFChannelOnlineStatus * calChOnlineSt = (AliTOFChannelOnlineStatus*)fCalStatus->At(ich);
      calChOnlineSt->SetStatus(((AliTOFChannelOnlineStatus*)currentCalStatus->At(ich))->GetStatus());
    }
  }

  // processing pulser
  
  TFile * daqFile=0x0;
  TH1S *h1=0x0;
  
  //retrieving Pulser data 
  TList* listPulser = GetFileSources(kDAQ, "PULSER");
  if (listPulser !=0x0 && listPulser->GetEntries()!=0)
    {
      AliInfo("The following sources produced files with the id PULSER");
      listPulser->Print();
      for (Int_t jj=0;jj<listPulser->GetEntries();jj++){
	TObjString * str = dynamic_cast<TObjString*> (listPulser->At(jj));
	AliInfo(Form("found source %s", str->String().Data()));
	// file to be stored run per run
	TString fileNamePulser = GetFile(kDAQ, "PULSER", str->GetName());
	if (fileNamePulser.Length()>0){
	  // storing refernce data
	  AliInfo(Form("Got the file %s, now we can process pulser data.", fileNamePulser.Data()));
	  daqFile = new TFile(fileNamePulser.Data(),"READ");
	  h1 = (TH1S*) daqFile->Get("hTOFpulser");
	  for (Int_t ibin=0;ibin<kSize;ibin++){
	    if ((h1->GetBinContent(ibin+1))!=-1){
	      if ((htofPulser->GetBinContent(ibin+1))==-1){
		htofPulser->SetBinContent(ibin+1,h1->GetBinContent(ibin+1));
	      }
	      else {
		Log(Form("Something strange occurred during Pulser run, channel %i already read by another LDC, please check!",ibin));
	      }
	    }
	  }
	  
	  // elaborating infos
	  Double_t mean =0;
	  Int_t nread=0;
	  Int_t nreadNotEmpty=0;
	  for (Int_t ientry=1;ientry<=h1->GetNbinsX();ientry++){

	    /* check whether channel has been read out during current run.
	     * if the status is bad it means it has not been read out.
	     * in this case skip channel in order to not affect the mean */ 
	    if (((AliTOFChannelOnlineStatus *)fFEEStatus->At(ientry-1))->GetStatus() == AliTOFChannelOnlineStatus::kTOFHWBad)
		    continue;

	    if (h1->GetBinContent(ientry)==-1) continue;
	    else {
	      if (h1->GetBinContent(ientry)>0) {
		nreadNotEmpty++;
		AliDebug(1,Form(" channel %i is ok with entry = %f; so far %i channels added ",ientry-1,h1->GetBinContent(ientry),nreadNotEmpty));
	      }
	      mean+=h1->GetBinContent(ientry);
	      nread++;
	    }
	  }
	  if (nread!=0) {
		  mean/=nread;
		  AliDebug(1,Form(" nread =  %i , mean = %f",nread,mean));
		  for (Int_t ich =0;ich<fNChannels;ich++){
			  AliTOFChannelOnlineStatus * chSt = (AliTOFChannelOnlineStatus *)fCalStatus->At(ich);
			  if (h1->GetBinContent(ich+1)==-1) continue;
			  AliDebug(1,Form(" channel %i ",ich));
			  AliDebug(1,Form(" channel status before pulser = %i",(Int_t)chSt->GetStatus()));
			  
			  /* check whether channel has been read out during current run.
			   * if the status is bad it means it has not been read out.
			   * in this case skip channel in order to leave its status 
			   * unchanged */
       			  if (((AliTOFChannelOnlineStatus *)fFEEStatus->At(ich))->GetStatus() == AliTOFChannelOnlineStatus::kTOFHWBad)
			  continue;
			  
			  if (h1->GetBinContent(ich+1)<0.05*mean){
				  chSt->SetStatus(chSt->GetStatus()|AliTOFChannelOnlineStatus::kTOFPulserBad);  // bad status for pulser
				  AliDebug(1,Form(" channel status after pulser = %i",(Int_t)chSt->GetStatus()));
			  }
			  else {
				  chSt->SetStatus(chSt->GetStatus()|AliTOFChannelOnlineStatus::kTOFPulserOk);  // bad status for pulser
				  AliDebug(1,Form(" channel status after pulser = %i",(Int_t)chSt->GetStatus()));
			  }
		  }
	  }
	  else {
		  Log("No channels read!! No action taken, keeping old status");
	  }

	  daqFile->Close();
	  delete daqFile;
	  delete h1;
	}
	
	else{
	  Log("The input data file from DAQ (pulser) was not found, TOF exiting from Shuttle "); 
	  return 10;//return error code for failure in retrieving Ref Data 
	}
	
      }
      delete listPulser;
    }
  
  else{
    Log("The input data file list from DAQ (pulser) was not found, TOF exiting from Shuttle "); 
    return 10;//return error code for failure in retrieving Ref Data 
  }	

  //storing in OCDB  

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Chiara Zampolli");
  metaData.SetComment("This preprocessor fills a TObjArray object for Pulser data.");
  AliInfo("Storing Calibration Data from Pulser Run");
  resultPulser = Store("Calib","Pulser",fCalStatus, &metaData,0,kTRUE);
  if(!resultPulser){
    Log("Some problems occurred while storing online object resulting from Pulser data processing");
    return 11;//return error code for problems in storing Pulser data 
  }

  if(fStoreRefData){
    
    AliCDBMetaData metaDataHisto;
    metaDataHisto.SetBeamPeriod(0);
    metaDataHisto.SetResponsible("Chiara Zampolli");
    char comment[200];
    sprintf(comment,"This preprocessor stores the result of the pulser run");
    metaDataHisto.SetComment(comment);
    AliInfo("Storing Reference Data");
    resultPulserRef = StoreReferenceData("Calib","PulserData",htofPulser, &metaDataHisto);
    if (!resultPulserRef){
      Log("some problems occurred::No Reference Data for pulser stored, TOF exiting from Shuttle");
      return 9;//return error code for failure in storing Ref Data 
    }
  }
  
  daqFile=0;

  return 0;
}
//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::ProcessNoiseData()
{

  // Processing Noise Run data for TOF channel status

  Log("Processing Noise");

  TH1::AddDirectory(0);

  Bool_t resultNoiseRef=kFALSE;
  Bool_t resultNoise=kFALSE;

  static const Int_t kSize = AliTOFGeometry::NPadXSector()*AliTOFGeometry::NSectors();
  TH1F * htofNoise = new TH1F("hTOFnoise","histo with signals on TOF during pulser", kSize,-0.5,kSize-0.5);
  for (Int_t ibin =1;ibin<=kSize;ibin++){
    htofNoise->SetBinContent(ibin,-1);
  }

  // retrieving last stored noise object, and copying

  AliCDBEntry *cdbEntry = GetFromOCDB("Calib","Noise");
  if (cdbEntry!=0x0){
    TObjArray *currentCalStatus = (TObjArray*)cdbEntry->GetObject();
    for (Int_t ich = 0; ich<fNChannels; ich ++){
      AliTOFChannelOnlineStatus * calChOnlineSt = (AliTOFChannelOnlineStatus*)fCalStatus->At(ich);
      calChOnlineSt->SetStatus(((AliTOFChannelOnlineStatus*)currentCalStatus->At(ich))->GetStatus());
    }
  }

  // processing noise
  
  TFile * daqFile=0x0;
  TH1F * h1=0x0;
  
  //retrieving Noise data 
  TList* listNoise = GetFileSources(kDAQ, "NOISE");
  if (listNoise !=0x0 && listNoise->GetEntries()!=0)
    {
      AliInfo("The following sources produced files with the id NOISE");
      listNoise->Print();
      for (Int_t jj=0;jj<listNoise->GetEntries();jj++){
	TObjString * str = dynamic_cast<TObjString*> (listNoise->At(jj));
	AliInfo(Form("found source %s", str->String().Data()));
	// file to be stored run per run
	TString fileNameNoise = GetFile(kDAQ, "NOISE", str->GetName());
	if (fileNameNoise.Length()>0){
	  // storing refernce data
	  AliInfo(Form("Got the file %s, now we can process noise data.", fileNameNoise.Data()));
	  daqFile = new TFile(fileNameNoise.Data(),"READ");
	  h1 = (TH1F*) daqFile->Get("hTOFnoise");
	  for (Int_t ibin=0;ibin<kSize;ibin++){
	    if ((h1->GetBinContent(ibin+1))!=-1){
	      if ((htofNoise->GetBinContent(ibin+1))==-1){
		htofNoise->SetBinContent(ibin+1,h1->GetBinContent(ibin+1));
	      }
	      else {
		Log(Form("Something strange occurred during Noise run, channel %i already read by another LDC, please check!",ibin));
	      }
	    }
	  }
	  // elaborating infos
	  for (Int_t ich =0;ich<fNChannels;ich++){
	    AliTOFChannelOnlineStatus * chSt = (AliTOFChannelOnlineStatus *)fCalStatus->At(ich);
	    if (h1->GetBinContent(ich+1)==-1) continue;
	    AliDebug(1,Form( " channel %i",ich));
	    AliDebug(1,Form( " channel status before noise = %i",(Int_t)chSt->GetStatus()));

	    /* check whether channel has been read out during current run.
	     * if the status is bad it means it has not been read out.
	     * in this case skip channel in order to leave its status 
	     * unchanged */
	    if (((AliTOFChannelOnlineStatus *)fFEEStatus->At(ich))->GetStatus() == AliTOFChannelOnlineStatus::kTOFHWBad)
	      continue;

	    if (h1->GetBinContent(ich+1)>=1){  // setting limit for noise to 1 kHz
	      chSt->SetStatus(chSt->GetStatus()|AliTOFChannelOnlineStatus::kTOFNoiseBad);  // bad status for noise
	      AliDebug(1,Form( " channel status after noise = %i",(Int_t)chSt->GetStatus()));
	    }
	    else {
	      chSt->SetStatus(chSt->GetStatus()|AliTOFChannelOnlineStatus::kTOFNoiseOk);  // bad status for noise
	      AliDebug(1,Form(" channel status after noise = %i",(Int_t)chSt->GetStatus()));
	    }
	  }
 
	  daqFile->Close();
	  delete daqFile;
	  delete h1;

	}
	
	else{
	  Log("The input data file from DAQ (noise) was not found, TOF exiting from Shuttle "); 
	  return 13;//return error code for failure in retrieving Ref Data 
	}
	
      }
      delete listNoise;
    }
  else{
    Log("The input data file list from DAQ (noise) was not found, TOF exiting from Shuttle "); 
    return 13;//return error code for failure in retrieving Ref Data 
  }	
  
  daqFile=0;

  //storing in OCDB

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Chiara Zampolli");
  metaData.SetComment("This preprocessor fills a TObjArray object for Noise data.");
  AliInfo("Storing Calibration Data from Noise Run");
  resultNoise = Store("Calib","Noise",fCalStatus, &metaData,0,kTRUE);
  if(!resultNoise){
    Log("Some problems occurred while storing online object resulting from Noise data processing");
    return 14;//return error code for problems in storing Noise data 
  }

  if(fStoreRefData){
    
    AliCDBMetaData metaDataHisto;
    metaDataHisto.SetBeamPeriod(0);
    metaDataHisto.SetResponsible("Chiara Zampolli");
    char comment[200];
    sprintf(comment,"This preprocessor stores the result of the noise run, TOF exiting from Shuttle ");
    metaDataHisto.SetComment(comment);
    AliInfo("Storing Reference Data");
    resultNoiseRef = StoreReferenceData("Calib","NoiseData",htofNoise, &metaDataHisto);
    if (!resultNoiseRef){
      Log("some problems occurred::No Reference Data for noise stored");
      return 12;//return error code for failure in storing Ref Data 
    }
  }

  return 0;
}
//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::ProcessFEEData()
{
  // Processing Pulser Run data for TOF channel status
  // dummy for the time being

  Log("Processing FEE");

  Bool_t updateOCDB = kFALSE;
  AliTOFFEEReader feeReader;
  AliTOFChannelOnlineStatus *currentChannel = NULL, *storedChannel = NULL;
  TObjArray *currentFEE = fFEEStatus;
  
  TH1C hCurrentFEE("hCurrentFEE","histo with current FEE channel status", fNChannels, 0, fNChannels);
  
  /* load current TOF FEE config from DCS FXS, parse, 
   * fill current FEE histogram and set FEE status */
  
  const char * nameFile = GetFile(kDCS,"TofFeeMap",""); 
  AliInfo(Form("nameFile = %s",nameFile));
  if (nameFile == NULL) {
	  return 15;
  } 
  feeReader.LoadFEEConfig(nameFile);
  feeReader.ParseFEEConfig();
  /* loop over channels */
  for (Int_t iChannel = 0; iChannel < fNChannels; iChannel++) {
    currentChannel = (AliTOFChannelOnlineStatus *)currentFEE->At(iChannel);
    if (!currentChannel)
      continue;
    /* channel enabled. set FEE channel status ok */
    if (feeReader.IsChannelEnabled(iChannel)) {
      currentChannel->SetStatus(AliTOFChannelOnlineStatus::kTOFHWOk);
      hCurrentFEE.SetBinContent(iChannel + 1, 1);
    }
    /* channel disabled. set FEE channel status bad */
    else {
      currentChannel->SetStatus(AliTOFChannelOnlineStatus::kTOFHWBad);
    }
  }
  
  /* load stored TOF FEE from OCDB and compare it with current FEE.
   * if stored FEE is different from current FEE set update flag.
   * if there is no stored FEE in OCDB set update flag */
  
  AliCDBEntry *cdbEntry = GetFromOCDB("Calib","FEE");
  /* no CDB entry found. set update flag */
  if (cdbEntry == NULL) {
    updateOCDB = kTRUE;
  }
  /* CDB entry OK. loop over channels */
  else {
    TObjArray *storedFEE = (TObjArray *)cdbEntry->GetObject();
    for (Int_t iChannel = 0; iChannel < fNChannels; iChannel++){
      currentChannel = (AliTOFChannelOnlineStatus *)currentFEE->At(iChannel);
      storedChannel = (AliTOFChannelOnlineStatus *)storedFEE->At(iChannel);
      /* compare current FEE channel status with stored one 
       * if different set update flag and break loop */
      if (currentChannel->GetStatus() != storedChannel->GetStatus()) {
	updateOCDB = kTRUE;
	break;
      }
    }
  }

  /* check whether we don't have to store reference data.
   * in this case we return without errors. */
  if (fStoreRefData) {
	  /* store reference data */
	  AliCDBMetaData metaDataHisto;
	  metaDataHisto.SetBeamPeriod(0);
	  metaDataHisto.SetResponsible("Roberto Preghenella");
	  metaDataHisto.SetComment("This preprocessor stores the FEE referece data of the current run.");
	  AliInfo("Storing FEE reference data");
	  /* store FEE reference data */
	  if (!StoreReferenceData("Calib", "FEEData", &hCurrentFEE, &metaDataHisto)) {
		  /* failed */
		  Log("problems while storing FEE reference data");
		  return 18; /* error return code for problems while storing FEE reference data */
	  }
  }

  /* check whether we don't need to update OCDB.
   * in this case we can return without errors and
   * the current FEE is stored in the fFEEStatus TObjArray. */
  if (!updateOCDB) {
    AliInfo("TOF FEE config has not changed. Do not overwrite stored file.");
    return 0; /* return ok */
  }

  /* update the OCDB with the current FEE since even 
   * a little difference has been detected. */

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Roberto Preghenella");
  metaData.SetComment("This preprocessor fills a TObjArray object for FEE data.");
  AliInfo("Storing FEE data from current run");
  /* store FEE data */
  if (!Store("Calib", "FEE", fFEEStatus, &metaData, 0, kTRUE)) {
    /* failed */
    Log("problems while storing FEE data object");
    return 17; /* return error code for problems  while storing FEE data */
  }

  /* everything fine. return */

  return 0;

}

//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::Process(TMap* dcsAliasMap)
{
  //
  //
  //

  TString runType = GetRunType();
  Log(Form("RunType %s",runType.Data()));
  
  // processing 

  /* always process FEE data */
  Int_t iresultFEE = ProcessFEEData();
  if (iresultFEE != 0)
    return iresultFEE;

  if (runType == "PULSER") {
    Int_t iresultPulser = ProcessPulserData();
    return iresultPulser; 
  }
  
  if (runType == "NOISE") { // for the time being associating noise runs with pedestal runs; proper run type to be defined 
    Int_t iresultNoise = ProcessNoiseData();
    return iresultNoise; 
  }
  
  if (runType == "PHYSICS") {
    Int_t iresultDCS = ProcessDCSDataPoints(dcsAliasMap);
    if (iresultDCS != 0) {
      return iresultDCS;
    }
    else { 
      Int_t iresultDAQ = ProcessOnlineDelays();
      return iresultDAQ;
    }
  }

  // storing
  return 0;
}


