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

///////////////////////////////////////////////////////////////////////////////
// EMCAL Preprocessor class. It runs by Shuttle at the end of the run,
// calculates stuff to be posted in OCDB
//
// Author: Boris Polichtchouk, 4 October 2006
// Adapted for EMCAL by Gustavo Conesa Balbastre, October 2006
// Updated by David Silvermyr May 2008, based on TPC code
///////////////////////////////////////////////////////////////////////////////

//Root
#include "TFile.h"
#include "TTree.h"
#include "TEnv.h"
#include "TParameter.h"

#include <TTimeStamp.h>

//AliRoot
#include "AliShuttleInterface.h"
#include "AliEMCALPreprocessor.h"
#include "AliLog.h"
#include "AliDCSValue.h"
#include "AliCDBMetaData.h"
#include "AliEMCALTriggerTRUDCSConfig.h"
#include "AliEMCALTriggerDCSConfig.h"
#include "AliCaloCalibPedestal.h"
#include "AliCaloCalibSignal.h"
#include "AliEMCALSensorTempArray.h"

const Int_t kValCutTemp = 100;               // discard temperatures > 100 degrees
const Int_t kDiffCutTemp = 5;	             // discard temperature differences > 5 degrees
const TString kPedestalRunType = "PEDESTAL";  // pedestal run identifier
const TString kPhysicsRunType = "PHYSICS";   // physics run identifier
const TString kStandAloneRunType = "STANDALONE_BC"; // standalone run identifier
const TString kAmandaTemp = "PT_%02d.Temperature"; // Amanda string for temperature entries
//const Double_t kFitFraction = 0.7;                 // Fraction of DCS sensor fits required 
const Double_t kFitFraction = -1.0;          // Don't require minimum number of fits during commissioning 

const TString kMetaResponsible = "David Silvermyr";
//legacy comments and return codes from TPC
const TString kMetaComment = "Preprocessor AliEMCAL data base entries.";
const int kReturnCodeNoInfo = 9;
const int kReturnCodeNoObject = 2;
const int kReturnCodeNoEntries = 1;

const int kNTRU = 30; // From 2011; 10 SuperModules (SM) * 3 TRU per SM

ClassImp(AliEMCALPreprocessor)
  
//_______________________________________________________________________________________
AliEMCALPreprocessor::AliEMCALPreprocessor() :
  AliPreprocessor("EMC",0),
  fConfEnv(0), 
  fTemp(0), 
  fConfigOK(kTRUE)
{
  //default constructor
}

//_______________________________________________________________________________________
AliEMCALPreprocessor::AliEMCALPreprocessor(AliShuttleInterface* shuttle):
  AliPreprocessor("EMC",shuttle),
  fConfEnv(0), 
  fTemp(0), 
  fConfigOK(kTRUE)
{
  // Constructor AddRunType(kPedestalRunType);
  
  // define run types to be processed
  AddRunType(kPedestalRunType);
  AddRunType(kPhysicsRunType);
}

//______________________________________________________________________________________________
AliEMCALPreprocessor::AliEMCALPreprocessor(const AliEMCALPreprocessor&  ) :
  AliPreprocessor("EMCAL",0),
  fConfEnv(0), fTemp(0), fConfigOK(kTRUE)
{
  Fatal("AliEMCALPreprocessor", "copy constructor not implemented");
}

// assignment operator; use copy ctor to make life easy.
//______________________________________________________________________________________________
AliEMCALPreprocessor& AliEMCALPreprocessor::operator = (const AliEMCALPreprocessor &source ) 
{
  // assignment operator; use copy ctor
  if (&source == this) return *this;
  
  new (this) AliEMCALPreprocessor(source);
  return *this;
}

//____________________________________________________________________________
AliEMCALPreprocessor::~AliEMCALPreprocessor()
{
  // destructor
  if (fTemp) delete fTemp;
}

//______________________________________________________________________________________________
void AliEMCALPreprocessor::Initialize(Int_t run, UInt_t startTime,
				      UInt_t endTime)
{
  // Creates AliTestDataDCS object -- start maps half an hour beforre actual run start
  UInt_t startTimeLocal = startTime-1800;
  AliPreprocessor::Initialize(run, startTimeLocal, endTime);
  
  AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
	       TTimeStamp((time_t)startTime,0).AsString(),
	       TTimeStamp((time_t)endTime,0).AsString()));
  
  // Preprocessor configuration
  AliCDBEntry* entry = GetFromOCDB("Config", "Preprocessor");
  if (entry) fConfEnv = (TEnv*) entry->GetObject();
  if ( fConfEnv==0 ) {
    Log("AliEMCALPreprocessor: Preprocessor Config OCDB entry missing.\n");
    fConfigOK = kFALSE;
    return;
  }
  
  // Temperature sensors
  TTree *confTree = 0;
  
  TString tempConf = fConfEnv->GetValue("Temperature","ON");
  tempConf.ToUpper();
  if (tempConf != "OFF" ) {
    entry = GetFromOCDB("Config", "Temperature");
    if (entry) confTree = (TTree*) entry->GetObject();
    if ( confTree==0 ) {
      Log("AliEMCALPreprocessor: Temperature Config OCDB entry missing.\n");
      fConfigOK = kFALSE;
      return;
    }
    fTemp = new AliEMCALSensorTempArray(startTimeLocal, fEndTime, confTree, kAmandaTemp);
    fTemp->SetValCut(kValCutTemp);
    fTemp->SetDiffCut(kDiffCutTemp);
  }
  
  return;
}

//______________________________________________________________________________________________
UInt_t AliEMCALPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into EMCAL calibrations objects
  // Amanda servers provide information directly through dcsAliasMap
  
  if (!fConfigOK) return kReturnCodeNoInfo;
  UInt_t result = 0;
  TObjArray *resultArray = new TObjArray();
  TString errorHandling = fConfEnv->GetValue("ErrorHandling","ON");
  errorHandling.ToUpper();
  TObject * status;
  
  UInt_t dcsResult=0;
  if (errorHandling == "OFF" ) {
    if (!dcsAliasMap) dcsResult = kReturnCodeNoEntries;
    else if (dcsAliasMap->GetEntries() == 0 ) dcsResult = kReturnCodeNoEntries;  
    status = new TParameter<int>("dcsResult",dcsResult);
    resultArray->Add(status);
  } 
  else {
    if (!dcsAliasMap) return kReturnCodeNoInfo;
    else if (dcsAliasMap->GetEntries() == 0 ) return kReturnCodeNoInfo;
  }
  
    
  TString runType = GetRunType();
  
  // Temperature sensors are processed by AliEMCALCalTemp
  TString tempConf = fConfEnv->GetValue("Temperature","ON");
  tempConf.ToUpper();
  if (tempConf != "OFF" && dcsAliasMap ) {
    UInt_t tempResult = MapTemperature(dcsAliasMap);
    result=tempResult;
    status = new TParameter<int>("tempResult",tempResult);
    resultArray->Add(status);
  }
  // Trigger configuration processing: only for Physics runs
  TString triggerConf = fConfEnv->GetValue("Trigger","ON");
  triggerConf.ToUpper();
  if( runType == kPhysicsRunType ) {
    //  if (triggerConf != "OFF" && dcsAliasMap ) {
    if ( dcsAliasMap ) {
      UInt_t triggerResult = MapTriggerConfig(dcsAliasMap);
      AliInfo(Form("triggerConf %s\n", triggerConf.Data()));
      result+=triggerResult;
      status = new TParameter<int>("triggerResult",triggerResult);
      resultArray->Add(status);
    }
  }
  
  // Other calibration information will be retrieved through FXS files
  //  examples:
  //    TList* fileSourcesDAQ = GetFile(AliShuttleInterface::kDAQ, "pedestals");
  //    const char* fileNamePed = GetFile(AliShuttleInterface::kDAQ, "pedestals", "LDC1");
  //
  //    TList* fileSourcesHLT = GetFile(AliShuttleInterface::kHLT, "calib");
  //    const char* fileNameHLT = GetFile(AliShuttleInterface::kHLT, "calib", "LDC1");
  
  // PEDESTAL ENTRIES:
  
  if ( runType == kPedestalRunType ) {
    Int_t numSources = 1;
    Int_t pedestalSource[2] = {AliShuttleInterface::kDAQ, AliShuttleInterface::kHLT} ;
    TString source = fConfEnv->GetValue("Pedestal","DAQ");
    source.ToUpper();
    if (source != "OFF" ) { 
      if ( source == "HLT") pedestalSource[0] = AliShuttleInterface::kHLT;
      if (!GetHLTStatus()) pedestalSource[0] = AliShuttleInterface::kDAQ;
      if (source == "HLTDAQ" ) {
        numSources=2;
        pedestalSource[0] = AliShuttleInterface::kHLT;
        pedestalSource[1] = AliShuttleInterface::kDAQ;
      }
      if (source == "DAQHLT" ) numSources=2;
      UInt_t pedestalResult=0;
      for (Int_t i=0; i<numSources; i++ ) {	
        pedestalResult = ExtractPedestals(pedestalSource[i]);
        if ( pedestalResult == 0 ) break;
      }
      result += pedestalResult;
      status = new TParameter<int>("pedestalResult",pedestalResult);
      resultArray->Add(status);
    }
  }
  
  // SIGNAL/LED ENTRIES:
  if( runType == kPhysicsRunType ) {
    Int_t numSources = 1;
    Int_t signalSource[2] = {AliShuttleInterface::kDAQ,AliShuttleInterface::kHLT} ;
    TString source = fConfEnv->GetValue("Signal","DAQ");
    source.ToUpper();
    if ( source != "OFF") { 
      if ( source == "HLT") signalSource[0] = AliShuttleInterface::kHLT;
      if (!GetHLTStatus()) signalSource[0] = AliShuttleInterface::kDAQ;
      if (source == "HLTDAQ" ) {
        numSources=2;
        signalSource[0] = AliShuttleInterface::kHLT;
        signalSource[1] = AliShuttleInterface::kDAQ;
      }
      if (source == "DAQHLT" ) numSources=2;
      UInt_t signalResult=0;
      for (Int_t i=0; i<numSources; i++ ) {	
        signalResult = ExtractSignal(signalSource[i]);
        if ( signalResult == 0 ) break;
      }
      result += signalResult;
      status = new TParameter<int>("signalResult",signalResult);
      resultArray->Add(status);
    }
  }
  
  
  // overall status at the end
  if (errorHandling == "OFF" ) {
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible(kMetaResponsible);
    metaData.SetComment("Preprocessor AliEMCAL status.");
    Bool_t storeOK = Store("Calib", "PreprocStatus", resultArray, &metaData, 0, kFALSE);
    resultArray->Delete();
    result = 0;
    if ( !storeOK )  result=1;
    return result;
  } 
  else { 
    return result;
  }
  
}
//______________________________________________________________________________________________
UInt_t AliEMCALPreprocessor::MapTemperature(TMap* dcsAliasMap)
{ // extract DCS temperature maps. Perform fits to save space
  UInt_t result=0;

  TMap *map = fTemp->ExtractDCS(dcsAliasMap);
  if (map) {
    fTemp->MakeSplineFit(map);
    Double_t fitFraction = 1.0*fTemp->NumFits()/fTemp->NumSensors(); 
    if (fitFraction > kFitFraction ) {
      AliInfo(Form("Temperature values extracted, fits performed.\n"));
    } 
    else { 
      Log ("Too few temperature maps fitted. \n");
      result = kReturnCodeNoInfo;
    }
  } 
  else {
    Log("No temperature map extracted. \n");
    result = kReturnCodeNoInfo;
  }
  delete map;
  // Now store the final CDB file
  
  if ( result == 0 ) { // some info was found
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible(kMetaResponsible);
    metaData.SetComment(kMetaComment);
    
    Bool_t storeOK = Store("Calib", "Temperature", fTemp, &metaData, 0, kFALSE);
    if ( !storeOK )  result=1;
  }
  
  return result;
}

//______________________________________________________________________________________________
UInt_t AliEMCALPreprocessor::MapTriggerConfig(TMap* dcsAliasMap)
{ // extract DCS trigger info
  AliInfo(Form("Get TRU info from DCS DPs.\n"));
  Int_t i, iTRU;
  char buf[100];

  AliDCSValue *dcsVal;
  TObjArray *arrL0ALGSEL, *arrPEAKFINDER, *arrGLOBALTHRESH, *arrCOSMTHRESH;
  TObjArray *arrMASK[6];

  // overall object to hold STU and DCS config info
  // DS comment: for now only holds TRU info, i.e. only partially filled
  // (STU info only in raw data header; unfortunately not also picked up via DCS DPs)
  AliEMCALTriggerDCSConfig *trigConfig = new AliEMCALTriggerDCSConfig();

  // loop through all TRUs
  bool debug = true; // debug flag for AliInfo printouts for each TRU
  for( iTRU = 0; iTRU < kNTRU; iTRU++){
    if (debug) AliInfo( Form("iTRU %d \n", iTRU) );
    // get the shuttled values
    sprintf( buf, "EMC_TRU%02d_L0ALGSEL", iTRU );
    arrL0ALGSEL = (TObjArray*) dcsAliasMap->GetValue( buf );
    sprintf( buf, "EMC_TRU%02d_PEAKFINDER", iTRU );
    arrPEAKFINDER = (TObjArray*) dcsAliasMap->GetValue( buf );
    sprintf( buf, "EMC_TRU%02d_GLOBALTHRESH", iTRU );
    arrGLOBALTHRESH = (TObjArray*) dcsAliasMap->GetValue( buf );
    sprintf( buf, "EMC_TRU%02d_COSMTHRESH", iTRU );
    arrCOSMTHRESH = (TObjArray*) dcsAliasMap->GetValue( buf );
    
    for( i = 0; i < 6; i++ ){
      sprintf( buf, "EMC_TRU%02d_MASK%d", iTRU, i );
      arrMASK[i] = (TObjArray*) dcsAliasMap->GetValue( buf );
    }
    
    // fill the objects
    AliEMCALTriggerTRUDCSConfig* truConfig = trigConfig->GetTRUDCSConfig(iTRU);
    if( ! truConfig ){
      AliWarning( Form("EMC DCS TRU%02d config not retrieved!\n", iTRU ));
    }

    // get last entries. fill the TRU object
    if( ! arrL0ALGSEL ){
      AliWarning( Form("EMC DCS TRU%02d L0ALGSEL alias not found!\n", iTRU ));
    }
    else{
      if (debug) AliInfo( Form("arrL0ALGSEL has %d entries \n", arrL0ALGSEL->GetEntries()) );
      if ( arrL0ALGSEL->GetEntries() > 0 ) {
	dcsVal = (AliDCSValue *) arrL0ALGSEL->At( arrL0ALGSEL->GetEntries() - 1 );
	if (dcsVal) truConfig->SetL0SEL( dcsVal->GetUInt() );
      }
    }

    if( ! arrPEAKFINDER ){
      AliWarning( Form("EMC DCS TRU%02d PEAKFINDER alias not found!\n", iTRU ));
    }
    else{
      if (debug) AliInfo( Form("arrPEAKFINDER has %d entries \n", arrPEAKFINDER->GetEntries()) );
      if ( arrPEAKFINDER->GetEntries() > 0 ) {
	dcsVal = (AliDCSValue *) arrPEAKFINDER->At( arrPEAKFINDER->GetEntries() - 1 );
	if (dcsVal) truConfig->SetSELPF( dcsVal->GetUInt() );
      }
    }

    if( ! arrGLOBALTHRESH ){
      AliWarning( Form("EMC DCS TRU%02d GLOBALTHRESH alias not found!\n", iTRU ));
    }
    else{
      if (debug) AliInfo( Form("arrGLOBALTHRESH has %d entries \n", arrGLOBALTHRESH->GetEntries()) );
      if ( arrGLOBALTHRESH->GetEntries() > 0 ) {
	dcsVal = (AliDCSValue *) arrGLOBALTHRESH->At( arrGLOBALTHRESH->GetEntries() - 1 );
	if (dcsVal) truConfig->SetGTHRL0( dcsVal->GetUInt() );
      }
    }

    if( ! arrCOSMTHRESH ){
      AliWarning( Form("EMC DCS TRU%02d COSMTHRESH alias not found!\n", iTRU ));
    }
    else{
      if (debug) AliInfo( Form("arrCOSMTHRESH has %d entries \n", arrCOSMTHRESH->GetEntries()) );
      if ( arrCOSMTHRESH->GetEntries() > 0 ) {
	dcsVal = (AliDCSValue *) arrCOSMTHRESH->At( arrCOSMTHRESH->GetEntries() - 1 );
	if (dcsVal) truConfig->SetL0COSM( dcsVal->GetUInt() );
      }
    }
    
    for( i = 0; i < 6; i++ ){
      if( ! arrMASK[i] ){
	AliWarning( Form("EMC DCS TRU%02d MASK%d alias not found!\n", iTRU, i ));
      }
      else{
	if (debug) AliInfo( Form("arrMASK[%d] has %d entries \n", i, arrMASK[i]->GetEntries()) );
	if ( arrMASK[i]->GetEntries() > 0 ) {
	  dcsVal = (AliDCSValue *) arrMASK[i]->At( arrMASK[i]->GetEntries() - 1 );
	  if (dcsVal) truConfig->SetMaskReg( dcsVal->GetUInt(), i );
	}
      }
    }
    
  } // TRUs
  AliInfo(Form("TRU info retrieved.\n"));
  // save the objects
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible(kMetaResponsible);
  metaData.SetComment(kMetaComment); 
      
  Bool_t retCode = Store("Calib", "Trigger", trigConfig, &metaData, 0, kFALSE);
  AliInfo(Form("TRU info stored.\n"));
  return retCode;
}

//______________________________________________________________________________________________
UInt_t AliEMCALPreprocessor::ExtractPedestals(Int_t sourceFXS)
{
  //  Read pedestal file from file exchange server
  //  Only store if new pedestal info is available
  //
  UInt_t result=0;

  AliCaloCalibPedestal *calibPed = new AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal);
  calibPed->Init();

  TList* list = GetFileSources(sourceFXS,"pedestals");
  if (list && list->GetEntries()>0) {
    
    //  loop through all files from LDCs

    int changes = 0;
    UInt_t index = 0;
    while (list->At(index)!=NULL) {
      TObjString* fileNameEntry = (TObjString*) list->At(index);
      if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "pedestals",
				   fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  Log ("Error opening pedestal file.");
	  result = kReturnCodeNoObject;
	  break;
	}
        AliCaloCalibPedestal *calPed;
	f->GetObject("emcCalibPedestal",calPed);
        if ( !calPed ) {
	  Log ("No pedestal calibration object in file.");
	  result = kReturnCodeNoObject;
	  break;
	}
	if ( calPed->GetNEvents()>0 && calPed->GetNChanFills()>0 ) {
	  // add info for the modules available in the present file
	  Bool_t status = calibPed->AddInfo(calPed);
	  if (status) { changes++; }
	}
	
        delete calPed; 
        f->Close();
      }
      index++;
    }  // while(list)
    
    //
    //  Store updated pedestal entry to OCDB
    //
    if (changes>0) {
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0);
      metaData.SetResponsible(kMetaResponsible);
      metaData.SetComment(kMetaComment); 
      
      Bool_t storeOK = StoreReferenceData("Calib", "Pedestals", calibPed, &metaData);
      if ( !storeOK ) result++;
    }
  } 
  else {
    Log ("Error: no entries in input file list!");
    result = kReturnCodeNoEntries;
  }
  
  return result;
}

//______________________________________________________________________________________________
UInt_t AliEMCALPreprocessor::ExtractSignal(Int_t sourceFXS)
{ //  Read signal file from file exchange server
  //  Only store if new signal info is available
  //
  UInt_t result=0;
  AliCaloCalibSignal *calibSig = new AliCaloCalibSignal(AliCaloCalibSignal::kEmCal); 
  
  TList* list = GetFileSources(sourceFXS,"signal");
  if (list && list->GetEntries()>0) {

    //  loop through all files from LDCs
    
    int changes = 0;
    UInt_t index = 0;
    while (list->At(index)!=NULL) {
      TObjString* fileNameEntry = (TObjString*) list->At(index);
      if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "signal",
				   fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  Log ("Error opening signal file.");
	  result = kReturnCodeNoObject;
	  break;
	}
	AliCaloCalibSignal *calSig;
	f->GetObject("emcCalibSignal",calSig);
        if ( !calSig ) {
	  Log ("No signal calibration object in file.");
	  result = kReturnCodeNoObject;
	  break;
	}
	if ( calSig->GetNEvents()>0 ) {
	  // add info for the modules available in the present file
	  Bool_t status = calibSig->AddInfo(calSig);
	  if (status) { changes++; }
	}
	
        delete calSig; 
        f->Close();
      }
      index++;
    }  // while(list)
    
    //
    //  Store updated signal entry to OCDB
    //
    if (changes>0) {
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0);
      metaData.SetResponsible(kMetaResponsible);
      metaData.SetComment(kMetaComment); 
      
      Bool_t storeOK = Store("Calib", "LED", calibSig, &metaData, 0, kFALSE);
      if ( !storeOK ) result++;
    }
  } 
  else {
    Log ("Error: no entries in input file list!");
    result = kReturnCodeNoEntries;
  }

  return result;
}


