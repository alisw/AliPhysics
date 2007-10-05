/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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


#include "AliTPCPreprocessor.h"
#include "AliShuttleInterface.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliTPCSensorTempArray.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalibPedestal.h"
#include "AliTPCCalibPulser.h"
#include "AliTPCCalibCE.h"
#include "TFile.h"
#include "TTree.h"
#include "TEnv.h"

#include <TTimeStamp.h>

const Int_t kValCutTemp = 100;               // discard temperatures > 100 degrees
const Int_t kDiffCutTemp = 5;	             // discard temperature differences > 5 degrees
const TString kPedestalRunType = "PEDESTAL_RUN";  // pedestal run identifier
const TString kPulserRunType = "PULSER_RUN";   // pulser run identifier

//
// This class is the SHUTTLE preprocessor for the TPC detector.
// It contains several components, this far the part containing
// temperatures is implemented
//

ClassImp(AliTPCPreprocessor)

//______________________________________________________________________________________________
AliTPCPreprocessor::AliTPCPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TPC",shuttle),
  fConfEnv(0), fTemp(0), fPressure(0), fConfigOK(kTRUE), fROC(0)
{
  // constructor
  fROC = AliTPCROC::Instance();
}
//______________________________________________________________________________________________
// AliTPCPreprocessor::AliTPCPreprocessor(const AliTPCPreprocessor& org) :
//   AliPreprocessor(org),
//   fConfEnv(0), fTemp(0), fPressure(0), fConfigOK(kTRUE)
// {
//   // copy constructor not implemented
//   //   -- missing underlying copy constructor in AliPreprocessor
//
//   Fatal("AliTPCPreprocessor", "copy constructor not implemented");
//
// //  fTemp = new AliTPCSensorTempArray(*(org.fTemp));
// }

//______________________________________________________________________________________________
AliTPCPreprocessor::~AliTPCPreprocessor()
{
  // destructor

  delete fTemp;
  delete fPressure;
}
//______________________________________________________________________________________________
AliTPCPreprocessor& AliTPCPreprocessor::operator = (const AliTPCPreprocessor& )
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}


//______________________________________________________________________________________________
void AliTPCPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliTestDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

  // Preprocessor configuration

	AliCDBEntry* entry = GetFromOCDB("Config", "Preprocessor");
        if (entry) fConfEnv = (TEnv*) entry->GetObject();
        if ( fConfEnv==0 ) {
           Log("AliTPCPreprocsessor: Preprocessor Config OCDB entry missing.\n");
	   fConfigOK = kFALSE;
           return;
        }

  // Temperature sensors

        TTree *confTree = 0;
	entry = GetFromOCDB("Config", "Temperature");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           Log("AliTPCPreprocsessor: Temperature Config OCDB entry missing.\n");
	   fConfigOK = kFALSE;
	   return;
        }
        fTemp = new AliTPCSensorTempArray(fStartTime, fEndTime, confTree);
	fTemp->SetValCut(kValCutTemp);
	fTemp->SetDiffCut(kDiffCutTemp);

  // Pressure sensors

        confTree=0;
	entry=0;
	entry = GetFromOCDB("Config", "Pressure");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           Log("AliTPCPreprocsessor: Pressure Config OCDB entry missing.\n");
	   fConfigOK = kFALSE;
	   return;
        }
	fPressure = new AliDCSSensorArray(fStartTime, fEndTime, confTree);

}

//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into TPC calibrations objects

  // Amanda servers provide information directly through dcsAliasMap

  if (!dcsAliasMap) return 9;
  if (dcsAliasMap->GetEntries() == 0 ) return 9;
  if (!fConfigOK) return 9;

  TString runType = GetRunType();

  // Temperature sensors are processed by AliTPCCalTemp


  UInt_t tempResult = MapTemperature(dcsAliasMap);
  UInt_t result=tempResult;

  // Pressure sensors

  UInt_t pressureResult = MapPressure(dcsAliasMap);
  result += pressureResult;

  // Other calibration information will be retrieved through FXS files
  //  examples:
  //    TList* fileSourcesDAQ = GetFile(AliShuttleInterface::kDAQ, "pedestals");
  //    const char* fileNamePed = GetFile(AliShuttleInterface::kDAQ, "pedestals", "LDC1");
  //
  //    TList* fileSourcesHLT = GetFile(AliShuttleInterface::kHLT, "calib");
  //    const char* fileNameHLT = GetFile(AliShuttleInterface::kHLT, "calib", "LDC1");

  // pedestal entries

  if(runType == kPedestalRunType) {
    Int_t pedestalSource = AliShuttleInterface::kDAQ;
    TString source = fConfEnv->GetValue("Pedestal","DAQ");
    source.ToUpper();
    if ( source == "HLT" ) pedestalSource = AliShuttleInterface::kHLT;
    if (!GetHLTStatus()) pedestalSource = AliShuttleInterface::kDAQ;
    UInt_t pedestalResult = ExtractPedestals(pedestalSource);
    result += pedestalResult;

  }

  // pulser trigger processing

  if( runType == kPulserRunType ) {
    Int_t pulserSource = AliShuttleInterface::kDAQ;
    TString source = fConfEnv->GetValue("Pulser","DAQ");
    source.ToUpper();
    if ( source == "HLT" ) pulserSource = AliShuttleInterface::kHLT;
    if (!GetHLTStatus()) pulserSource = AliShuttleInterface::kDAQ;
    UInt_t pulserResult = ExtractPulser(pulserSource);
    result += pulserResult;

  }

  // Central Electrode processing

  if( false ) {    // CE input file not generated yet
    Int_t ceSource = AliShuttleInterface::kDAQ;
    TString source = fConfEnv->GetValue("CE","DAQ");
    source.ToUpper();
    if ( source == "HLT" ) ceSource = AliShuttleInterface::kHLT;
    if (!GetHLTStatus()) ceSource = AliShuttleInterface::kDAQ;
    UInt_t ceResult = ExtractCE(ceSource);
    result += ceResult;

  }

  return result;
}
//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::MapTemperature(TMap* dcsAliasMap)
{

   // extract DCS temperature maps. Perform fits to save space

  UInt_t result=0;
  TMap *map = fTemp->ExtractDCS(dcsAliasMap);
  if (map) {
    fTemp->MakeSplineFit(map);
    AliInfo(Form("Temperature values extracted, fits performed.\n"));
  } else {
    Log("AliTPCPreprocsessor: no temperature map extracted. \n");
    result=9;
  }
  delete map;
  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "Temperature", fTemp, &metaData, 0, kFALSE);
        if ( !storeOK )  result=1;

   }

   return result;

}
//______________________________________________________________________________________________

UInt_t AliTPCPreprocessor::MapPressure(TMap* dcsAliasMap)
{

   // extract DCS temperature maps. Perform fits to save space

  UInt_t result=0;
  TMap *map = fPressure->ExtractDCS(dcsAliasMap);
  if (map) {
    fPressure->MakeSplineFit(map);
    AliInfo(Form("Pressure values extracted, fits performed.\n"));
  } else {
    Log("AliTPCPreprocsessor: no atmospheric pressure map extracted. \n");
    result=9;
  }
  delete map;
  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "Pressure", fPressure, &metaData, 0, 0);
        if ( !storeOK ) result=1;

   }

   return result;

}


//______________________________________________________________________________________________

UInt_t AliTPCPreprocessor::ExtractPedestals(Int_t sourceFXS)
{
 //
 //  Read pedestal file from file exchage server
 //  Keep original entry from OCDB in case no new pedestals are available
 //
 AliTPCCalPad *calPadPed=0;
 AliCDBEntry* entry = GetFromOCDB("Calib", "Pedestals");
 if (entry) calPadPed = (AliTPCCalPad*)entry->GetObject();
 if ( calPadPed==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC pedestal entry available.\n");
     calPadPed = new AliTPCCalPad("PedestalsMean","PedestalsMean");
 }

 AliTPCCalPad *calPadRMS=0;
 entry = GetFromOCDB("Calib", "Noise");
 if (entry) calPadRMS = (AliTPCCalPad*)entry->GetObject();
 if ( calPadRMS==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC noise entry available.\n");
     calPadRMS = new AliTPCCalPad("PedestalsRMS","PedestalsRMS");
 }


 UInt_t result=0;

 Int_t nSectors = fROC->GetNSectors();
 TList* list = GetFileSources(sourceFXS,"pedestals");
 
 if (list && list->GetEntries()>0) {

//  loop through all files from LDCs

    UInt_t index = 0;
    while (list->At(index)!=NULL) {
     TObjString* fileNameEntry = (TObjString*) list->At(index);
     if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "pedestals",
	                                 fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  Log ("Error opening pedestal file.");
	  result =2;
	  break;
	}
        AliTPCCalibPedestal *calPed;
	f->GetObject("AliTPCCalibPedestal",calPed);

        //  replace entries for the sectors available in the present file

        for (Int_t sector=0; sector<nSectors; sector++) {
           AliTPCCalROC *rocPed=calPed->GetCalRocPedestal(sector, kFALSE);
           if ( rocPed )  calPadPed->SetCalROC(rocPed,sector);
           AliTPCCalROC *rocRMS=calPed->GetCalRocRMS(sector, kFALSE);
           if ( rocRMS )  calPadRMS->SetCalROC(rocRMS,sector);
        }
      }
     ++index;
    }  // while(list)
//
//  Store updated pedestal entry to OCDB
//
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Haavard Helstrup");
    metaData.SetComment("Preprocessor AliTPC data base entries.");

    Bool_t storeOK = Store("Calib", "Pedestals", calPadPed, &metaData, 0, kTRUE);
    if ( !storeOK ) ++result;
    storeOK = Store("Calib", "PadNoise", calPadRMS, &metaData, 0, kTRUE);
    if ( !storeOK ) ++result;

  } else {
    Log ("Error: no entries!");
    result = 1;
  }

  return result;
}

//______________________________________________________________________________________________


UInt_t AliTPCPreprocessor::ExtractPulser(Int_t sourceFXS)
{
 //
 //  Read pulser calibration file from file exchage server
 //  Keep original entry from OCDB in case no new pulser calibration is available
 //
 TObjArray    *pulserObjects=0;
 AliTPCCalPad *pulserTmean=0;
 AliTPCCalPad *pulserTrms=0;
 AliTPCCalPad *pulserQmean=0;
 AliCDBEntry* entry = GetFromOCDB("Calib", "Pulser");
 if (entry) pulserObjects = (TObjArray*)entry->GetObject();
 if ( pulserObjects==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC pulser entry available.\n");
     pulserObjects = new TObjArray;    
 }

 pulserTmean = (AliTPCCalPad*)pulserObjects->FindObject("PulserTmean");
 if ( !pulserTmean ) {
    pulserTmean = new AliTPCCalPad("PulserTmean","PulserTmean");
    pulserObjects->Add(pulserTmean);
 }
 pulserTrms = (AliTPCCalPad*)pulserObjects->FindObject("PulserTrms");
 if ( !pulserTrms )  { 
    pulserTrms = new AliTPCCalPad("PulserTrms","PulserTrms");
    pulserObjects->Add(pulserTrms);
 }
 pulserQmean = (AliTPCCalPad*)pulserObjects->FindObject("PulserQmean");
 if ( !pulserQmean )  { 
    pulserQmean = new AliTPCCalPad("PulserQmean","PulserQmean");
    pulserObjects->Add(pulserQmean);
 }


 UInt_t result=0;

 Int_t nSectors = fROC->GetNSectors();
 TList* list = GetFileSources(sourceFXS,"pulser");
 
 if (list && list->GetEntries()>0) {

//  loop through all files from LDCs

    UInt_t index = 0;
    while (list->At(index)!=NULL) {
     TObjString* fileNameEntry = (TObjString*) list->At(index);
     if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "pulser",
	                                 fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  Log ("Error opening pulser file.");
	  result =2;
	  break;
	}
        AliTPCCalibPulser *calPulser;
	f->GetObject("AliTPCCalibPulser",calPulser);

        //  replace entries for the sectors available in the present file

        for (Int_t sector=0; sector<nSectors; sector++) {
           AliTPCCalROC *rocTmean=calPulser->GetCalRocT0(sector);
           if ( rocTmean )  pulserTmean->SetCalROC(rocTmean,sector);
           AliTPCCalROC *rocTrms=calPulser->GetCalRocRMS(sector);
           if ( rocTrms )  pulserTrms->SetCalROC(rocTrms,sector);
           AliTPCCalROC *rocQmean=calPulser->GetCalRocQ(sector);
           if ( rocQmean )  pulserQmean->SetCalROC(rocQmean,sector);
        }
      }
     ++index;
    }  // while(list)
//
//  Store updated pedestal entry to OCDB
//
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Haavard Helstrup");
    metaData.SetComment("Preprocessor AliTPC data base entries.");

    Bool_t storeOK = Store("Calib", "Pulser", pulserObjects, &metaData, 0, kTRUE);
    if ( !storeOK ) ++result;
    
  } else {
    Log ("Error: no entries!");
    result = 1;
  }

  return result;
}

UInt_t AliTPCPreprocessor::ExtractCE(Int_t sourceFXS)
{
 //
 //  Read Central Electrode file from file exchage server
 //  Keep original entry from OCDB in case no new CE calibration is available
 //
 TObjArray    *ceObjects=0;
 AliTPCCalPad *ceTmean=0;
 AliTPCCalPad *ceTrms=0;
 AliTPCCalPad *ceQmean=0;
 AliCDBEntry* entry = GetFromOCDB("Calib", "CE");
 if (entry) ceObjects = (TObjArray*)entry->GetObject();
 if ( ceObjects==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC central electrode entry available.\n");
     ceObjects = new TObjArray;    
 }

 ceTmean = (AliTPCCalPad*)ceObjects->FindObject("CETmean");
 if ( !ceTmean ) {
    ceTmean = new AliTPCCalPad("CETmean","CETmean");
    ceObjects->Add(ceTmean);
 }
 ceTrms = (AliTPCCalPad*)ceObjects->FindObject("CETrms");
 if ( !ceTrms )  { 
    ceTrms = new AliTPCCalPad("CETrms","CETrms");
    ceObjects->Add(ceTrms);
 }
 ceQmean = (AliTPCCalPad*)ceObjects->FindObject("CEQmean");
 if ( !ceQmean )  { 
    ceQmean = new AliTPCCalPad("CEQmean","CEQmean");
    ceObjects->Add(ceQmean);
 }


 UInt_t result=0;

 Int_t nSectors = fROC->GetNSectors();
 TList* list = GetFileSources(sourceFXS,"CE");
 
 if (list && list->GetEntries()>0) {

//  loop through all files from LDCs

    UInt_t index = 0;
    while (list->At(index)!=NULL) {
     TObjString* fileNameEntry = (TObjString*) list->At(index);
     if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "CE",
	                                 fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  Log ("Error opening central electrode file.");
	  result =2;
	  break;
	}
        AliTPCCalibCE *calCE;
	f->GetObject("AliTPCCalibCE",calCE);

        //  replace entries for the sectors available in the present file

        for (Int_t sector=0; sector<nSectors; sector++) {
           AliTPCCalROC *rocTmean=calCE->GetCalRocT0(sector);
           if ( rocTmean )  ceTmean->SetCalROC(rocTmean,sector);
           AliTPCCalROC *rocTrms=calCE->GetCalRocRMS(sector);
           if ( rocTrms )  ceTrms->SetCalROC(rocTrms,sector);
           AliTPCCalROC *rocQmean=calCE->GetCalRocQ(sector);
           if ( rocQmean )  ceQmean->SetCalROC(rocQmean,sector);
        }
      }
     ++index;
    }  // while(list)
//
//  Store updated pedestal entry to OCDB
//
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Haavard Helstrup");
    metaData.SetComment("Preprocessor AliTPC data base entries.");

    Bool_t storeOK = Store("Calib", "CE", ceObjects, &metaData, 0, kTRUE);
    if ( !storeOK ) ++result;
    
  } else {
    Log ("Error: no entries!");
    result = 1;
  }

  return result;
}
