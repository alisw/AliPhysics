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

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1S.h>
#include <TH2S.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TTimeStamp.h>

#include "AliCDBMetaData.h"
#include "AliLog.h"
#include "AliTOFChannelOnline.h"
#include "AliTOFDataDCS.h"
#include "AliTOFGeometry.h"
#include "AliTOFPreprocessor.h"

class TF1;
class AliDCSValue;
class AliTOFGeometry;

// TOF preprocessor class.
// It takes data from DCS and passes them to the class AliTOFDataDCS, which
// processes them. The result is then written to the CDB.
// analogously, it takes data form DAQ (both at Run level and inclusive - 
// of all the runs - level, processes them, and stores both Reference Data
// and Online Calibration files in the CDB. 


ClassImp(AliTOFPreprocessor)

const Int_t AliTOFPreprocessor::fgkBinRangeAve    = 13; // number of bins where to calculate the mean 
const Double_t AliTOFPreprocessor::fgkThrPar    = 0.013; // parameter used to trigger the calculation of the delay

//_____________________________________________________________________________

AliTOFPreprocessor::AliTOFPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TOF", shuttle),
  fData(0),
  fh2(0),
  fCal(0),
  fNChannels(0),
  fStoreRefData(kTRUE)
{
  // constructor

}

//_____________________________________________________________________________

AliTOFPreprocessor::~AliTOFPreprocessor()
{
  // destructor
  delete fData;
  fData = 0;
  delete fh2;
  fh2 = 0;
  fCal->Clear();
  delete fCal;
  fCal = 0;
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
	fh2 = 0x0;
	fNChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
	fCal = new TObjArray(fNChannels);
	fCal->SetOwner();
	for (Int_t ich = 0; ich<fNChannels; ich ++){
	  AliTOFChannelOnline * calChOnline = new AliTOFChannelOnline();
	  fCal->AddAt(calChOnline,ich);
	}
}

//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliTOFDataDCS object
  // return codes:
  // return=0 : all ok
  // return=1 : no DCS input data Map
  // return=2 : no DCS input data processing
  // return=3 : no DCS processed data was stored
  // return=4 : no DAQ input for Ref Data
  // return=5 : failed to store Ref data
  // return=6 : failed to retrieve DAQ data for calibration 
  // return=7 : problems in histos in the input DAQ file 
  // return=8 : failed to store Online Delays

  TH1::AddDirectory(0);

  Bool_t resultDCSMap=kFALSE;
  Bool_t resultDCSStore=kFALSE;
  Bool_t resultDAQ=kFALSE;
  Bool_t resultDAQRef=kFALSE;

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
      resultDCSStore = Store("Calib","DCSData",fData, &metaDataDCS);
      if (!resultDCSStore){
	Log("Some problems occurred while storing DCS data results in OCDB, TOF exiting from Shuttle");
	return 3;// return error Code for processed DCS data not stored 
	
      }
    }
  }
  
  // processing DAQ
  
  TFile * daqFile=0x0;
  
  if(fStoreRefData){
    //retrieving data at Run level
    TList* list = GetFileSources(kDAQ, "RUNLevel");
    if (list)
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
      }
    else{
      Log("The input data file list from DAQ (run-level) was not found, TOF exiting from Shuttle "); 
      return 4;//return error code for failure in retrieving Ref Data 
    }	
  }


//Total files, with cumulative histos
  
  TList* listTot = GetFileSources(kDAQ, "DELAYS");
  if (listTot)
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
	    for (Int_t ich=0;ich<kSize;ich++){
	      TH1S *h1 = new TH1S("h1","h1",kNBins,kXBinmin-0.5,kNBins*1.+kXBinmin-0.5);
	      for (Int_t ibin=0;ibin<kNBins;ibin++){
		h1->SetBinContent(ibin+1,fh2->GetBinContent(ich+1,ibin+1));
	      }
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
	      }
	    delete h1;
	    h1=0x0;
	    }
	  }
	  daqFile->Close();
	  delete daqFile;
	  AliCDBMetaData metaData;
	  metaData.SetBeamPeriod(0);
	  metaData.SetResponsible("Chiara Zampolli");
	  metaData.SetComment("This preprocessor fills an AliTOFCal object.");
	  AliInfo("Storing Calibration Data");
	  resultDAQ = Store("Calib","OnlineDelay",fCal, &metaData);
          if(!resultDAQ){
	    Log("Some problems occurred while storing DAQ data processing results");
	    return 8;//return error code for problems in storing DAQ data 
	  }
	}
	else{
	  Log("The Cumulative data file from DAQ does not exist, TOF exiting from Shuttle"); 
          return 6;//return error code for problems in retrieving DAQ data 
	}
      }
    }
  else{
    Log("Problem: no list for Cumulative data file from DAQ was found, TOF exiting from Shuttle");
    return 6; //return error code for problems in retrieving DAQ data 
  }

  daqFile=0;

  return 0;
}


