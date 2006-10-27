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
$Log$
Revision 1.1  2006/10/26 09:09:29  arcelli
prototype for the TOF Shuttle preprocessor (C.Zampolli)

*/  

#include "AliTOFPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"
#include "AliTOFDataDCS.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliTOFCalOnline.h"
#include "AliTOFChannelOnline.h"
#include "AliTOFGeometryV5.h"
#include "TTimeStamp.h"
#include "TH1F.h"

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

//_____________________________________________________________________________

AliTOFPreprocessor::AliTOFPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  AliPreprocessor(detector, shuttle),
  fData(0),
  fArray(0),
  fCal(0),
  fTOFGeometry(0)
{
  // constructor
}

//_____________________________________________________________________________

AliTOFPreprocessor::~AliTOFPreprocessor()
{
  // destructor
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
	fArray = 0x0;
	fTOFGeometry = new AliTOFGeometryV5();
	fCal = new AliTOFCalOnline(fTOFGeometry);
	fCal->CreateArray();
}

//_____________________________________________________________________________

UInt_t AliTOFPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliTOFDataDCS object

  UInt_t resultDCS=0;
  UInt_t resultDAQ=0;
  UInt_t resultDAQRef=0;
  UInt_t result=0;

  // processing DCS

  if (!dcsAliasMap){
    AliInfo(Form("No DCS map found "));
  }
  else {
  // The processing of the DCS input data is forwarded to AliTOFDataDCS
    fData->ProcessData(*dcsAliasMap);
    AliCDBMetaData metaDataDCS;
    metaDataDCS.SetBeamPeriod(0);
    metaDataDCS.SetResponsible("Chiara Zampolli");
    metaDataDCS.SetComment("This preprocessor fills an AliTOFDataDCS object.");     resultDCS = Store("Calib","DCSData",fData, &metaDataDCS);
    result+=resultDCS;
    if (!resultDCS){
      AliInfo(Form("some problems occurred while storing DCS data processing results"));
    }
  }

  // processing DAQ

  TFile * daqFile=0x0;
  TList* list = GetFileSources(kDAQ, "DELAYS");

  if (list)
    {
      AliInfo("The following sources produced files with the id DELAYS");
      list->Print();
      for (Int_t jj=0;jj<list->GetEntries();jj++){
	TObjString * str = dynamic_cast<TObjString*> (list->At(jj));
	AliInfo(Form("found source %s", str->String().Data()));
	// file to be stored run per run
	const char* fileNameRun = GetFile(kDAQ, "RUNLevel", str->GetName());
	if (fileNameRun){
	  AliInfo(Form("Got the file %s, now we can store the Reference Data for the current Run.", fileNameRun));
	  daqFile = new TFile(fileNameRun,"READ");
	  fArray = (TObjArray*) daqFile->Get("ciccio");
	  AliCDBMetaData metaDataHisto;
	  metaDataHisto.SetBeamPeriod(0);
	  metaDataHisto.SetResponsible("Chiara Zampolli");
	  metaDataHisto.SetComment("This preprocessor stores the array of histos object as Reference Data.");
	  resultDAQRef = StoreReferenceData("Calib","DAQData",fArray, &metaDataHisto);
	  result+=resultDAQRef*2;
	  if (!resultDAQRef){
	    AliInfo(Form("some problems occurred::No Reference Data stored"));
	  }
	  daqFile->Close();
	  delete daqFile;
	}

	else{
	  AliError(Form("The file %s does not exist",fileNameRun)); 
	}
       
	// file with summed histos, to extract calib params
	const char *fileName = GetFile(kDAQ, "DELAYS", str->GetName());
	if (fileName){
	  AliInfo(Form("Got the file %s, now we can extract some values.", fileName));

	  daqFile = new TFile(fileName,"READ");
	  fArray = (TObjArray*) daqFile->Get("ciccio");
	  if (!fArray){
	    AliInfo(Form("some problems occurred:: No histo array retrieved"));
	  }
	  
	  else {
	    Int_t nentries=(Int_t)fArray->GetEntries();
	    AliInfo(Form(" il numero di entries e' = %i ", nentries)); 
	    Int_t npads = fCal->NPads();
	    for (Int_t i=0 ; i<nentries; i++){
	      Bool_t found=kFALSE; 
	      TH1F * h1 = (TH1F*)fArray->At(i);
	      Float_t minContent=h1->Integral()*0.01; //for the time being 
	      //we use integral() since we simulate landau distribution
	      //Float_t minContent=h1->GetEntries()*0.01;
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
	      if(!found) AliInfo(Form("WARNING!!! no start of fit found for histo # %i",i));
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
	      if (i<npads) {
		AliTOFChannelOnline * ch = fCal->GetChannel(i);
		ch->SetDelay(mean);
	      }
	    }
	  }
	  daqFile->Close();
	  delete daqFile;
	  AliCDBMetaData metaData;
	  metaData.SetBeamPeriod(0);
	  metaData.SetResponsible("Chiara Zampolli");
	  metaData.SetComment("This preprocessor fills an AliTOFCal object.");
	  resultDAQ = Store("Calib","OnlineDelay",fCal, &metaData);
	  result+=resultDAQ*2*2;
	}
	else{
	  AliError(Form("The file %s does not exist",fileName)); 
	}
      }
    }
  else{
    AliInfo(Form("Problem: no list found"));
  }

  delete list;
  list = 0;
  daqFile=0;
  delete fData;
  fData = 0;
  delete fCal;
  fCal = 0;
  delete fArray;
  fArray = 0;
  return result;
}


