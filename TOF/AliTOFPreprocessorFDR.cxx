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
Revision 1.4  2007/12/05 13:55:18  zampolli
Initialization bug fixed.

Revision 1.3  2007/11/27 13:12:30  zampolli
CDB object run range upper limit extended to AliCDBRunRange::Infinity()

Revision 1.2  2007/11/27 07:24:41  zampolli
Log used, fData member removed

Revision 1.1  2007/11/24 18:36:27  zampolli
TOF Preprocessor for FDR

*/

//#include <Riostream.h>
//#include <stdio.h>
//#include <stdlib.h>

//#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TObjArray.h>
//#include <TObjString.h>
#include <TTimeStamp.h>

#include "AliCDBMetaData.h"
#include "AliLog.h"

//#include "AliTOFDataDCS.h"
//#include "AliTOFGeometry.h"
#include "AliTOFPreprocessorFDR.h"
#include "AliTOFFormatDCS.h"
#include "AliDCSValue.h"

// TOF preprocessor class.
// It takes data from DCS and passes them to the class AliTOFDataDCS, which
// processes them. The result is then written to the CDB.
// analogously, it takes data form DAQ (both at Run level and inclusive - 
// of all the runs - level, processes them, and stores both Reference Data
// and Online Calibration files in the CDB. 


ClassImp(AliTOFPreprocessorFDR)

//_____________________________________________________________________________

AliTOFPreprocessorFDR::AliTOFPreprocessorFDR(AliShuttleInterface* shuttle) :
  AliPreprocessor("TOF", shuttle),
  fStoreRefData(kTRUE)
{
  // constructor
  AddRunType("PHYSICS");

}

//_____________________________________________________________________________

AliTOFPreprocessorFDR::~AliTOFPreprocessorFDR()
{
  // destructor
}

//______________________________________________________________________________
void AliTOFPreprocessorFDR::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliTOFDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

}

//_____________________________________________________________________________

UInt_t AliTOFPreprocessorFDR::ProcessDCSDataPoints(TMap* aliasMap)
{
  // Fills data into a AliTOFDataDCS object
  // return codes:
  // return=0 : all ok
  // return=1 : no DCS input data Map
  // return=2 : no DCS processed data was stored in Ref Data
  // return=3 : no DCS processed data was stored in OCDB

  if (!aliasMap){
    Log("No DCS map found: TOF exiting from Shuttle");
    return 1;// return error Code for DCS input data not found 
  }

  else {

    AliDCSValue* aValue;
    AliDCSValue* aValue1;
    Float_t timeMin = (Float_t)fStartTime;
    Float_t timeMax = (Float_t)fEndTime;
    Float_t val=0;
    Float_t val1=0;
    Float_t time=0; 
    Float_t delta[2];
    Float_t timedelta[2];

    TH1::AddDirectory(0);
    
    Bool_t resultDCSMap=kFALSE;
    Bool_t resultDCSStore=kFALSE;
    
    TString aliasDP[4]={"tof_lv_i48_02","tof_lv_v48_02","tof_lv_i33_02","tof_lv_v33_02"};
    
    TObjArray *array = new TObjArray(4);
    array->SetOwner();

    AliTOFFormatDCS *lvI4802 = new AliTOFFormatDCS();
    AliTOFFormatDCS *lvV4802 = new AliTOFFormatDCS();
    AliTOFFormatDCS *lvI3302 = new AliTOFFormatDCS();
    AliTOFFormatDCS *lvV3302 = new AliTOFFormatDCS();
    
    array->AddAt(lvI4802,0);
    array->AddAt(lvV4802,1);
    array->AddAt(lvI3302,2);
    array->AddAt(lvV3302,3);

    // processing DCS
    
    for (Int_t i=0;i<4;i++){
      for (Int_t idelta =0;idelta<2;idelta++){
	delta[idelta]=0;
	timedelta[idelta]=0;
      }
      TObjArray *aliasArr = (TObjArray*) aliasMap->GetValue(aliasDP[i].Data());
      
      if(!aliasArr){
	Log(Form("Alias %s not found!", aliasDP[i].Data()));
	return kFALSE;
      }
      
      if(aliasArr->GetEntries()<3){
	Log(Form("Alias %s has just %d entries!",
		      aliasDP[i].Data(),aliasArr->GetEntries()));
	continue;
      }
      
      TIter iterarray(aliasArr);
      
      Int_t nentries = aliasArr->GetEntries();
      Int_t deltaTimeStamp = (Int_t) nentries/3;
      Int_t deltaTimeStamp1 = (Int_t) nentries/2;
      
      // filling aliases with 10 floats+1 Usign
      Int_t index = 0;
      for (Int_t k=0;k<3;k++){
	index = deltaTimeStamp*k;
	if (k==0) {
	  index=0;
	}
	else if (k==1) {
	  index=deltaTimeStamp1;
	} 
	else if (k==2) {
	  index=nentries-1; 
	}
	aValue = (AliDCSValue*) aliasArr->At(index);
	val = aValue->GetFloat();
	time = (Float_t) (aValue->GetTimeStamp());
	if (i==0){
	  AliDebug(1,Form("tof_lv_i48_02: setting value %i to %f at %f",k,val,time));
	  lvI4802->SetFloat(k,val);
	  lvI4802->SetTimeStampFloat(k,time);
	}
	else if (i==1){
	  AliDebug(1,Form("tof_lv_v48_02: setting value %i to %f at %f",k,val,time));
	  lvV4802->SetFloat(k,val);
	  lvV4802->SetTimeStampFloat(k,time);
	}
	else if (i==2){
	  AliDebug(1,Form("tof_lv_i33_02: setting value %i to %f at %f",k,val,time));
	  lvI3302->SetFloat(k,val);
	  lvI3302->SetTimeStampFloat(k,time);
	}
	else if (i==3){
	  AliDebug(1,Form("tof_lv_v33_02: setting value %i to %f at %f",k,val,time));
	  lvV3302->SetFloat(k,val);
	  lvV3302->SetTimeStampFloat(k,time);
	}
      }
  
      // computing the most significant variations
      
      Int_t deltamin = (Int_t)(60/(timeMax-timeMin)*nentries);
      Int_t klast = nentries-deltamin;
      
      for (Int_t k=0;k<klast;k++){
	aValue = (AliDCSValue*) aliasArr->At(k);
	aValue1 = (AliDCSValue*) aliasArr->At(k+deltamin);
	val = aValue->GetFloat();
	val1 = aValue1->GetFloat();
	if (delta[0]<=TMath::Abs(val1-val)) {
	  delta[0]=TMath::Abs(val1-val);
	  timedelta[0] = (Float_t)k;
	}
	if (delta[1]<=delta[0]) {
	  Float_t temp = delta[1];
	  Float_t timetemp = timedelta[1];
	  delta[1]=delta[0];
	  delta[0]=temp;
	  timedelta[1]=timedelta[0];
	  timedelta[0]=timetemp;
	}
      }
      
      for (Int_t kk=0;kk<2;kk++){
	if (i==0){
	  AliDebug(1,Form("tof_lv_i48: setting variation %i to %f at %f",kk,delta[kk],timedelta[kk]));
	  lvI4802->SetDelta(kk,delta[kk]);
	  lvI4802->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (i==1){
	  AliDebug(1,Form("tof_lv_v48: setting variation %i to %f at %f",kk,delta[kk],timedelta[kk]));
	  lvV4802->SetDelta(kk,delta[kk]);
	  lvV4802->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (i==2){
	  AliDebug(1,Form("tof_lv_i33: setting variation %i to %f at %f",kk,delta[kk],timedelta[kk]));
	  lvI3302->SetDelta(kk,delta[kk]);
	  lvI3302->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
	else if (i==3){
	  AliDebug(1,Form("tof_lv_v33: setting variation %i to %f at %f",kk,delta[kk],timedelta[kk]));
	  lvV3302->SetDelta(kk,delta[kk]);
	  lvV3302->SetTimeStampDelta(kk,(Float_t)timedelta[kk]);
	}
      }
    }

    AliCDBMetaData metaDataDCS;
    metaDataDCS.SetBeamPeriod(0);
    metaDataDCS.SetResponsible("Chiara Zampolli");
    metaDataDCS.SetComment("This preprocessor fills an AliTOFDataDCS object.");
    Log("Storing DCS Data");
    resultDCSStore = StoreReferenceData("Calib","DCSData",array, &metaDataDCS);
    if (!resultDCSStore){
      Log("Some problems occurred while storing DCS data results in Reference Data, TOF exiting from Shuttle");
      return 2;// return error Code for processed DCS data not stored 
      // in reference data
    }
    
    Log("Storing DCS Data in OCDB");
    resultDCSMap = Store("Calib","DCSData",array, &metaDataDCS,0,kTRUE);
    if (!resultDCSMap){
      Log("Some problems occurred while storing DCS data results in OCDB, TOF exiting from Shuttle");
      return 3;// return error Code for processed DCS data not stored 
      // in reference data
    }
    
    delete array;
    return 0;

  }

}
//_____________________________________________________________________________

UInt_t AliTOFPreprocessorFDR::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliTOFDataDCS object
  // return codes:
  // return=0  : all ok
  // return=1  : no DCS input data Map
  // return=2  : no DCS processed data was stored in Ref Data
  // return=3  : no DCS processed data was stored in OCDB

  TH1::AddDirectory(0);

  // processing 

  Int_t iresultDCS = ProcessDCSDataPoints(dcsAliasMap);
  if ((iresultDCS == 1) || (iresultDCS == 2) || (iresultDCS == 3)) return iresultDCS; 
  
  return 0;
}


