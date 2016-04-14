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
#include "AliTPCCalibRaw.h"
#include "AliTPCdataQA.h"
#include "ARVersion.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h" 
#include "TEnv.h"
#include "TParameter.h"

#include <TTimeStamp.h>

const Int_t kValCutTemp = 100;               // discard temperatures > 100 degrees
const Int_t kDiffCutTemp = 5;	             // discard temperature differences > 5 degrees
const Double_t kHighVoltageDifference = 1e-4; // don't record High Voltage points 
                                             // differing by less than 1e-4 from
					     // previous point.
const Double_t kGasCompositionDifference = 1e-3; // don't record Gas Composition points 
                                             // differing by less than 1e-3 from
					     // previous point
const TString kPedestalRunType = "PEDESTAL";  // pedestal run identifier
const TString kPulserRunType = "PULSER";     // pulser run identifier
const TString kPhysicsRunType = "PHYSICS";   // physics run identifier
const TString kCosmicRunType = "COSMIC";     // cosmic run identifier
const TString kLaserRunType = "LASER";       // laser run identifier
const TString kDaqRunType = "DAQ"; 	     // DAQ run identifier
const TString kStandaloneRunType = "STANDALONE";  // STANDALONE run identifier
const TString kAmandaTemp = "TPC_PT_%d_TEMPERATURE"; // Amanda string for temperature entries
const TString kAmandaDDL = "DDL%d";   // Amanda string for list of active DDLs
const Int_t  kNumDDL = 216;           // number of TPC DDLs
const Int_t  kFirstDDL = 768;         // identifier of first DDL
//const Double_t kFitFraction = 0.7;                 // Fraction of DCS sensor fits required              
const Double_t kFitFraction = -1.0;          // Don't require minimum number of fits in commissioning run 
const Int_t   kNumPressureSensors = 3;    // number of pressure sensors
const char* kPressureSensorNames[kNumPressureSensors] = {
                   "CavernAtmosPressure",
                   "CavernAtmosPressure2",
                   "SurfaceAtmosPressure" };
const Int_t  kMinCESectors = 32;      // minimum number of sectors (per side)
                                      // to accept CE calibration
      

//
// This class is the SHUTTLE preprocessor for the TPC detector.
//

ClassImp(AliTPCPreprocessor)

//______________________________________________________________________________________________
AliTPCPreprocessor::AliTPCPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TPC",shuttle),
  fConfEnv(0), fTemp(0), fHighVoltage(0), fHighVoltageStat(0), fGoofie(0),
  fPressure(0), fGasComposition(0), fConfigOK(kTRUE), fROC(0),
  fForceSingleRun(kFALSE)
{
  // constructor
  fROC = AliTPCROC::Instance();

  // define run types to be processed
  
  AddRunType(kPedestalRunType);
  AddRunType(kPulserRunType);
  AddRunType(kPhysicsRunType);
  AddRunType(kCosmicRunType);
  AddRunType(kLaserRunType);
  AddRunType(kDaqRunType);
  AddRunType(kStandaloneRunType);
  
}
//______________________________________________________________________________________________
 AliTPCPreprocessor::AliTPCPreprocessor(const AliTPCPreprocessor&  ) :
   AliPreprocessor("TPC",0),
   fConfEnv(0), fTemp(0), fHighVoltage(0), fHighVoltageStat(0), fGoofie(0),
   fPressure(0), fGasComposition(0), fConfigOK(kTRUE), fROC(0),
   fForceSingleRun(kFALSE)
 {

   Fatal("AliTPCPreprocessor", "copy constructor not implemented");
//
// //  fTemp = new AliTPCSensorTempArray(*(org.fTemp));
 }

//______________________________________________________________________________________________
AliTPCPreprocessor::~AliTPCPreprocessor()
{
  // destructor

  delete fTemp;
  delete fHighVoltage;
  delete fHighVoltageStat;
  delete fGoofie;
  delete fGasComposition;
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

  AliPreprocessor::Initialize(run, startTime, endTime);

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp((time_t)startTime,0).AsString(),
		TTimeStamp((time_t)endTime,0).AsString()));

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

       TString tempConf = fConfEnv->GetValue("Temperature","ON");
       tempConf.ToUpper();
       if (tempConf != "OFF" ) {
	entry = GetFromOCDB("Config", "Temperature");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           Log("AliTPCPreprocsessor: Temperature Config OCDB entry missing.\n");
	   fConfigOK = kFALSE;
	   return;
        }
        fTemp = new AliTPCSensorTempArray(startTime, endTime, confTree, kAmandaTemp);
	fTemp->SetValCut(kValCutTemp);
	fTemp->SetDiffCut(kDiffCutTemp);
       }

  // High voltage measurements

      TString hvConf = fConfEnv->GetValue("HighVoltage","ON");
      hvConf.ToUpper();
      if (hvConf != "OFF" ) { 
        confTree=0;
        entry=0;
        entry = GetFromOCDB("Config", "HighVoltage");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           Log("AliTPCPreprocsessor: High Voltage Config OCDB entry missing.\n");
           fConfigOK = kFALSE;
           return;
        }
        time_t timeStart = (time_t)(((TString)GetRunParameter("DAQ_time_start")).Atoi());
	time_t timeEnd = (time_t)(((TString)GetRunParameter("DAQ_time_end")).Atoi());
        fHighVoltage = new AliDCSSensorArray (UInt_t(timeStart), 
	                                    UInt_t(timeEnd), confTree);
      }

   // High voltage status values
     
      TString hvStatConf = fConfEnv->GetValue("HighVoltageStat","ON");
      hvStatConf.ToUpper();
      if (hvStatConf != "OFF" ) { 
        confTree=0;
        entry=0;
        entry = GetFromOCDB("Config", "HighVoltageStat");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           Log("AliTPCPreprocsessor: High Voltage Status Config OCDB entry missing.\n");
           fConfigOK = kFALSE;
           return;
        }
        fHighVoltageStat = new AliDCSSensorArray(startTime, endTime, confTree);
      }

   // Goofie values
     
      TString goofieConf = fConfEnv->GetValue("Goofie","ON");
      goofieConf.ToUpper();
      if (goofieConf != "OFF" ) { 
        confTree=0;
        entry=0;
        entry = GetFromOCDB("Config", "Goofie");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           Log("AliTPCPreprocsessor: Goofie Config OCDB entry missing.\n");
           fConfigOK = kFALSE;
           return;
        }
        fGoofie = new AliDCSSensorArray(startTime, endTime, confTree);
      }


  // Gas composition measurements

      TString gasConf = fConfEnv->GetValue("GasComposition","OFF");
      gasConf.ToUpper();
      if (gasConf != "OFF" ) { 
        confTree=0;
        entry=0;
        entry = GetFromOCDB("Config", "GasComposition");
        if (entry) confTree = (TTree*) entry->GetObject();
        if ( confTree==0 ) {
           Log("AliTPCPreprocsessor: Gas Composition Config OCDB entry missing.\n");
           fConfigOK = kFALSE;
           return;
        }
        time_t timeStart = (time_t)(((TString)GetRunParameter("DAQ_time_start")).Atoi());
	time_t timeEnd = (time_t)(((TString)GetRunParameter("DAQ_time_end")).Atoi());
        fGasComposition = new AliDCSSensorArray (UInt_t(timeStart), 
	                                    UInt_t(timeEnd), confTree);
      }
   // Pressure values
     
       TString runType = GetRunType();

       if( runType == kPhysicsRunType || 
        runType == kLaserRunType || runType == kStandaloneRunType ) {    
       TString pressureConf = fConfEnv->GetValue("Pressure","ON");
       pressureConf.ToUpper();
       if (pressureConf != "OFF" ) { 
         TClonesArray * array = new TClonesArray("AliDCSSensor",kNumPressureSensors); 
         for(Int_t j = 0; j < kNumPressureSensors; j++) {
           AliDCSSensor * sens = new ((*array)[j])AliDCSSensor;
           sens->SetStringID(kPressureSensorNames[j]);
         }
         fPressure = new AliDCSSensorArray(startTime, endTime, array);
       }
     }
}

//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into TPC calibrations objects

  // Amanda servers provide information directly through dcsAliasMap

  
  if (!fConfigOK) return 9;
  UInt_t result = 0;
  TObjArray *resultArray = new TObjArray();
  TString errorHandling = fConfEnv->GetValue("ErrorHandling","ON");
  errorHandling.ToUpper();
  TObject * status;

  UInt_t dcsResult=0;
  if (!dcsAliasMap) { 
     dcsResult=1;
  } else {
     if (dcsAliasMap->GetEntries() == 0 ) dcsResult=1;
  }  
  status = new TParameter<int>("dcsResult",dcsResult);
  resultArray->Add(status);


  TString runType = GetRunType();

  if ( dcsResult == 0 ) {

  // Temperature sensors are processed by AliTPCCalTemp

    TString tempConf = fConfEnv->GetValue("Temperature","ON");
    tempConf.ToUpper();
    if (tempConf != "OFF" ) {
      UInt_t tempResult = MapTemperature(dcsAliasMap);
      if ( tempConf != "TRY") result+=tempResult;
      status = new TParameter<int>("tempResult",tempResult);
      resultArray->Add(status);
    }

  // High Voltage recordings


    TString hvConf = fConfEnv->GetValue("HighVoltage","ON");
    hvConf.ToUpper();
    if (hvConf != "OFF" ) { 
     UInt_t hvResult = MapHighVoltage(dcsAliasMap);
     if (hvConf != "TRY") result+=hvResult;
     status = new TParameter<int>("hvResult",hvResult);
     resultArray->Add(status);
    }

    // Goofie values


    TString goofieConf = fConfEnv->GetValue("Goofie","ON");
    goofieConf.ToUpper();
    if (goofieConf != "OFF" ) { 
     UInt_t goofieResult = MapGoofie(dcsAliasMap);
     if (goofieConf != "TRY") result+=goofieResult;
     status = new TParameter<int>("goofieResult",goofieResult);
     resultArray->Add(status);
    }


    // Gas composition recordings


    TString gasConf = fConfEnv->GetValue("GasComposition","OFF");
    gasConf.ToUpper();
    if (gasConf != "OFF" ) { 
     UInt_t gasResult = MapGasComposition(dcsAliasMap);
     if (gasConf != "TRY") result+=gasResult;
     status = new TParameter<int>("gasResult",gasResult);
     resultArray->Add(status);
    }

    // Pressure values

    if( runType == kPhysicsRunType || 
      runType == kLaserRunType || runType == kStandaloneRunType ) {    

      TString pressureConf = fConfEnv->GetValue("Pressure","ON");
      pressureConf.ToUpper();
      if (pressureConf != "OFF" ) { 
       UInt_t pressureResult = MapPressure(dcsAliasMap);
       status = new TParameter<int>("pressureResult",pressureResult);
       resultArray->Add(status);
      }
    }
  }
  // Other calibration information will be retrieved through FXS files
  //  examples:
  //    TList* fileSourcesDAQ = GetFile(AliShuttleInterface::kDAQ, "pedestals");
  //    const char* fileNamePed = GetFile(AliShuttleInterface::kDAQ, "pedestals", "LDC1");
  //
  //    TList* fileSourcesHLT = GetFile(AliShuttleInterface::kHLT, "calib");
  //    const char* fileNameHLT = GetFile(AliShuttleInterface::kHLT, "calib", "LDC1");

  // pedestal entries

  if(runType == kPedestalRunType) {
    Int_t numSources = 1;
    Int_t pedestalSource[2] = {AliShuttleInterface::kDAQ,AliShuttleInterface::kHLT} ;
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

  // pulser trigger processing

  if(runType == kPulserRunType) {
    Int_t numSources = 1;
    Int_t pulserSource[2] = {AliShuttleInterface::kDAQ,AliShuttleInterface::kHLT} ;
    TString source = fConfEnv->GetValue("Pulser","DAQ");
    source.ToUpper();
    if ( source != "OFF") { 
     if ( source == "HLT") pulserSource[0] = AliShuttleInterface::kHLT;
     if (!GetHLTStatus()) pulserSource[0] = AliShuttleInterface::kDAQ;
     if (source == "HLTDAQ" ) {
         numSources=2;
	 pulserSource[0] = AliShuttleInterface::kHLT;
 	 pulserSource[1] = AliShuttleInterface::kDAQ;
     }
     if (source == "DAQHLT" ) numSources=2;
     if (source == "TRY" ) numSources=1;
     UInt_t pulserResult=0;
     for (Int_t i=0; i<numSources; i++ ) {	
       pulserResult = ExtractPulser(pulserSource[i]);
       if ( pulserResult == 0 ) break;
     }
     if (source != "TRY") result += pulserResult;
     status = new TParameter<int>("pulserResult",pulserResult);
     resultArray->Add(status);
    }
  }


// raw calibration processing

  if(runType == kPhysicsRunType) {
    Int_t numSources = 1;
    Int_t rawSource[2] = {AliShuttleInterface::kDAQ,AliShuttleInterface::kHLT} ;
    TString source = fConfEnv->GetValue("Raw","DAQ");
    source.ToUpper();
    if ( source != "OFF") { 
     if ( source == "HLT") rawSource[0] = AliShuttleInterface::kHLT;
     if (!GetHLTStatus()) rawSource[0] = AliShuttleInterface::kDAQ;
     if (source == "HLTDAQ" ) {
         numSources=2;
	 rawSource[0] = AliShuttleInterface::kHLT;
 	 rawSource[1] = AliShuttleInterface::kDAQ;
     }
     if (source == "DAQHLT" ) numSources=2;
     if (source == "TRY" ) numSources=1;
     UInt_t rawResult=0;
     for (Int_t i=0; i<numSources; i++ ) {	
       rawResult = ExtractRaw(rawSource[i]);
       if ( rawResult == 0 ) break;
     }
     if (source != "TRY" )result += rawResult;
     status = new TParameter<int>("rawResult",rawResult);
     resultArray->Add(status);
    }
  }


  // Altro configuration


  TString altroConf = fConfEnv->GetValue("AltroConf","ON");
  altroConf.ToUpper();
  if (altroConf != "OFF" ) { 
   UInt_t altroResult = ExtractAltro(AliShuttleInterface::kDCS,dcsAliasMap);
   if (altroConf != "TRY" ) result+=altroResult;
   status = new TParameter<int>("altroResult",altroResult);
   resultArray->Add(status);
  }


  // Central Electrode processing

  if( runType == kPhysicsRunType || 
      runType == kLaserRunType ) {    

    Int_t numSources = 1;
    Int_t ceSource[2] = {AliShuttleInterface::kDAQ,AliShuttleInterface::kHLT} ;
    TString source = fConfEnv->GetValue("CE","DAQ");
    source.ToUpper();
    if ( source != "OFF" ) { 
     if ( source == "HLT") ceSource[0] = AliShuttleInterface::kHLT;
     if (!GetHLTStatus()) ceSource[0] = AliShuttleInterface::kDAQ;
     if (source == "HLTDAQ" ) {
        numSources=2;
	ceSource[0] = AliShuttleInterface::kHLT;
	ceSource[1] = AliShuttleInterface::kDAQ;
     }
     if (source == "DAQHLT" ) numSources=2;
     if (source == "TRY" ) numSources=1;
     UInt_t ceResult=0;
     for (Int_t i=0; i<numSources; i++ ) {	
       ceResult = ExtractCE(ceSource[i]);
       if ( ceResult == 0 ) break;
     }

   // only flag error if CE result is missing from LASER runs
   //    -- for PHYSICS run do CE processing if data available
   
     if ( runType == kLaserRunType && source != "TRY"  && ceResult<10 ) result += ceResult;
     status = new TParameter<int>("ceResult",ceResult);
     resultArray->Add(status);

    numSources = 1;
    Int_t qaSource[2] = {AliShuttleInterface::kDAQ,AliShuttleInterface::kHLT} ;
    source = fConfEnv->GetValue("QA","DAQ");
    source.ToUpper();
    if ( source != "OFF" ) { 
     if ( source == "HLT") qaSource[0] = AliShuttleInterface::kHLT;
     if (!GetHLTStatus()) qaSource[0] = AliShuttleInterface::kDAQ;
     if (source == "HLTDAQ" ) {
        numSources=2;
	qaSource[0] = AliShuttleInterface::kHLT;
	qaSource[1] = AliShuttleInterface::kDAQ;
     }
     if (source == "DAQHLT" ) numSources=2;
     if (source == "TRY" ) numSources=1;
     UInt_t qaResult=0;
     for (Int_t i=0; i<numSources; i++ ) {	
       qaResult = ExtractQA(qaSource[i]);
       if ( qaResult == 0 ) break;
     }
//     result += qaResult;
     if ( qaResult !=0 ) Log ("ExtractQA failed, no QA entry available.");
     status = new TParameter<int>("qaResult",qaResult);
     resultArray->Add(status);
    }
   }
  }
  
// Store component status to OCDB

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Haavard Helstrup");
  metaData.SetAliRootVersion(ALIROOT_VERSION);
  metaData.SetComment("Preprocessor AliTPC status.");
  Bool_t storeOK = Store("Calib", "PreprocStatus", resultArray, &metaData, 0,  kFALSE);
  if (!storeOK) Log ("Unable to store preprocessor status entry");
 
  resultArray->Delete();
  delete resultArray;

  if (errorHandling == "OFF" ) return 0;
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
    Double_t fitFraction = 1.0*fTemp->NumFits()/fTemp->NumSensors(); 
    if (fitFraction > kFitFraction ) {
      AliInfo(Form("Temperature values extracted, fits performed.\n"));
    } else { 
      Log ("Too few temperature maps fitted. \n");
      result = 9;
    }
  } else {
    Log("No temperature map extracted. \n");
    result=9;
  }
  delete map;
  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetAliRootVersion(ALIROOT_VERSION);
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "Temperature", fTemp, &metaData, 0, kFALSE);
        if ( !storeOK )  result=1;

   }

   return result;

}
//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::MapPressure(TMap* dcsAliasMap)
{

   // extract DCS pressure maps. Perform fits to save space

  UInt_t result=0;
  TMap *map = fPressure->ExtractDCS(dcsAliasMap);
  if (map) {
    fPressure->MakeSplineFit(map);
    Double_t fitFraction = 1.0*fPressure->NumFits()/fPressure->NumSensors(); 
    if (fitFraction > kFitFraction ) {
      AliInfo(Form("Pressure values extracted, fits performed.\n"));
    } else { 
      Log ("Too few pressure maps fitted. \n");
      result = 9;
    }
  } else {
    Log("No pressure map extracted. \n");
    result=9;
  }
  delete map;
  return result;
}

//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::MapHighVoltage(TMap* dcsAliasMap)
{

   // extract DCS HV maps. Perform fits to save space

  UInt_t result=0;
  TMap *map = fHighVoltage->ExtractDCS(dcsAliasMap);
  if (map) {
    fHighVoltage->ClearFit();
    fHighVoltage->SetGraph(map);
    fHighVoltage->RemoveGraphDuplicates(kHighVoltageDifference);
               // don't keep new point if too similar to previous one
  } else {
    Log("No high voltage recordings extracted. \n");
    result=9;
  }
  delete map;

  TString hvStatConf = fConfEnv->GetValue("HighVoltageStat","ON");
  hvStatConf.ToUpper();
  if (hvStatConf != "OFF" ) { 
    TMap *map2 = fHighVoltageStat->ExtractDCS(dcsAliasMap);
    if (map2) {
      fHighVoltageStat->ClearFit();
      fHighVoltageStat->SetGraph(map2);
    } else {
       Log("No high voltage status recordings extracted. \n");
      result=9;
    }
    delete map2;

    // add status maps to high voltage sensor array

    fHighVoltage->AddSensors(fHighVoltageStat);
   }
  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetAliRootVersion(ALIROOT_VERSION);
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "HighVoltage", fHighVoltage, &metaData, 0, kFALSE);
        if ( !storeOK )  result=1;

   }

   return result;

}

//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::MapGoofie(TMap* dcsAliasMap)
{

   // extract DCS Goofie maps. Do not perform fits (low update rate)

  UInt_t result=0;

  TMap *map = fGoofie->ExtractDCS(dcsAliasMap);
  if (map) {
    fGoofie->ClearFit();
    fGoofie->SetGraph(map);
  } else {
    Log("No Goofie recordings extracted. \n");
    result=9;
  }
  delete map;

  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
        metaData.SetAliRootVersion(ALIROOT_VERSION);
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "Goofie", fGoofie, &metaData, 0, kFALSE);
        if ( !storeOK )  result=1;

   }

   return result;

}

//______________________________________________________________________________________________
UInt_t AliTPCPreprocessor::MapGasComposition(TMap* dcsAliasMap)
{

   // extract DCS HV maps. Perform fits to save space

  UInt_t result=0;
  TMap *map = fGasComposition->ExtractDCS(dcsAliasMap);
  if (map) {
    fGasComposition->ClearFit();
    fGasComposition->SetGraph(map);
    fGasComposition->RemoveGraphDuplicates(kGasCompositionDifference);
               // don't keep new point if too similar to previous one
  } else {
    Log("No gas composition recordings extracted. \n");
    result=9;
  }
  delete map;

  // Now store the final CDB file

  if ( result == 0 ) {
        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Haavard Helstrup");
	metaData.SetAliRootVersion(ALIROOT_VERSION);
	metaData.SetComment("Preprocessor AliTPC data base entries.");

	Bool_t storeOK = Store("Calib", "GasComposition", fGasComposition, &metaData, 0, kFALSE);
        if ( !storeOK )  result=1;

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
 AliTPCCalPad *calPadPedOCDB=0;
 AliCDBEntry* entry = GetFromOCDB("Calib", "Pedestals");
 if (entry) calPadPedOCDB = (AliTPCCalPad*)entry->GetObject();
 if ( calPadPedOCDB==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC pedestal entry available.\n");
 }

 AliTPCCalPad *calPadRMSOCDB=0;
 entry = GetFromOCDB("Calib", "PadNoise");
 if (entry) calPadRMSOCDB = (AliTPCCalPad*)entry->GetObject();
 if ( calPadRMSOCDB==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC noise entry available.\n");
 }

  AliTPCCalPad* calPadPed = new AliTPCCalPad("PedestalsMean","PedestalsMean");
  AliTPCCalPad* calPadRMS = new AliTPCCalPad("PedestalsRMS","PedestalsRMS");


 UInt_t result=0;

 Int_t nSectors = fROC->GetNSectors();
 TVectorD foundSectorsPed(nSectors);
 foundSectorsPed=0;
 TVectorD foundSectorsRMS(nSectors);
 foundSectorsRMS=0;

 TList* list = GetFileSources(sourceFXS,"pedestals");
 
 if (list && list->GetEntries()>0) {

//  loop through all files from LDCs

    Bool_t changed=false;
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
	f->GetObject("tpcCalibPedestal",calPed);
        if ( !calPed ) {
	  Log ("No pedestal calibration object in file.");
	  result = 2;
	  break;
	}

        //  replace entries for the sectors available in the present file

        changed=true;
        for (Int_t sector=0; sector<nSectors; sector++) {
           AliTPCCalROC *rocPed=calPed->GetCalRocPedestal(sector, kFALSE);
           if ( rocPed ) {
	       AliTPCCalROC* roc=calPadPed->GetCalROC(sector);
	       roc->Add(rocPed,1);
	       foundSectorsPed[sector]++;
	   }
           AliTPCCalROC *rocRMS=calPed->GetCalRocRMS(sector, kFALSE);
           if ( rocRMS )  {
               AliTPCCalROC* roc=calPadRMS->GetCalROC(sector);
	       roc->Add(rocRMS,1);
	       foundSectorsRMS[sector]++;
	   } 
        }
        delete calPed; 
        f->Close();
      }
     ++index;
    }  // while(list)

//
//  Check if calibration is complete -- otherwise take from old OCDB entry

// inner sectors
    for (Int_t sector=0; sector<nSectors/2; sector++) {
       if (foundSectorsPed[sector] < 1 ) {
          if (calPadPedOCDB) {
            AliTPCCalROC *rocOCDB=calPadPedOCDB->GetCalROC(sector);
	    calPadPed->SetCalROC(rocOCDB,sector);
	  } else {
	    const int mess_length=100;
            char message[mess_length];
	    snprintf(message,mess_length,"Missing pedestals for sector %d - also not available from previous OCDB entry.\n",
	             sector);
	    Log (message);
	    result = 2;
	  }
       }
       if (foundSectorsRMS[sector] < 1 ) {
          if (calPadRMSOCDB) {
            AliTPCCalROC *rocOCDB=calPadRMSOCDB->GetCalROC(sector);
	    calPadRMS->SetCalROC(rocOCDB,sector);
	  } else {
	    const int mess_length=100;
            char message[mess_length];
	    snprintf(message,mess_length,"Missing pedestal RMS for sector %d - also not available from previous OCDB entry.\n",
	             sector);
	    Log (message);
	    result = 2;
	  }       
       }
    }

// outer sectors -- two updates needed

    for (Int_t sector=nSectors/2; sector<nSectors; sector++) {
       if (foundSectorsPed[sector] < 2 ) {
          if (calPadPedOCDB) {
            AliTPCCalROC *rocOCDB=calPadPedOCDB->GetCalROC(sector);
	    calPadPed->SetCalROC(rocOCDB,sector);
	  } else {
	    const int mess_length=100;
            char message[mess_length];
	    snprintf(message,mess_length,"Missing pedestals for sector %d - also not available from previous OCDB entry.\n",
	             sector);
	    Log (message);
	    result = 2;
	  }

       }
       if (foundSectorsRMS[sector] < 2 ) {
          if (calPadRMSOCDB) {
            AliTPCCalROC *rocOCDB=calPadRMSOCDB->GetCalROC(sector);
	    calPadRMS->SetCalROC(rocOCDB,sector);
	  } else {
	    const int mess_length=100;
            char message[mess_length];
	    snprintf(message,mess_length,"Missing pedestal RMS for sector %d - also not available from previous OCDB entry.\n",
	             sector);
	    Log (message);
	    result = 2;
	  }       

       }       
    }

    
//
//  Store updated pedestal entry to OCDB
//
    if (changed) {
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Haavard Helstrup");
     metaData.SetAliRootVersion(ALIROOT_VERSION);
     metaData.SetComment("Preprocessor AliTPC data base entries."); 
 
     Bool_t storeOK = Store("Calib", "Pedestals", calPadPed, &metaData, 0, kTRUE ^ fForceSingleRun);
     if ( !storeOK ) ++result;
     storeOK = Store("Calib", "PadNoise", calPadRMS, &metaData, 0, kTRUE ^ fForceSingleRun);
     if ( !storeOK ) ++result;
    }
  } else {
    Log ("Error: no entries in pedestal file list!");
    result = 1;
  }

  delete calPadPed;
  delete calPadRMS;
  delete calPadPedOCDB;
  delete calPadRMSOCDB;

  return result;
}

//______________________________________________________________________________________________


UInt_t AliTPCPreprocessor::ExtractPulser(Int_t sourceFXS)
{
 //
 //  Read pulser calibration file from file exchage server
 //  Keep original entry from OCDB in case no new pulser calibration is available
 //
 
 TObjArray *pulserObjects = new TObjArray;
 TObjArray *pulserObjectsOCDB=0;   
 pulserObjects->SetOwner(kTRUE);

  
 AliCDBEntry* entry = GetFromOCDB("Calib", "Pulser");
 if (entry) pulserObjectsOCDB = (TObjArray*)entry->GetObject();
 if ( pulserObjectsOCDB==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC pulser entry available.\n");
 } else {
     pulserObjectsOCDB->SetOwner(kTRUE);
 }

 AliTPCCalPad *pulserTmean = new AliTPCCalPad("PulserTmean","PulserTmean");
 pulserObjects->Add(pulserTmean);
 
 AliTPCCalPad *pulserTrms = new AliTPCCalPad("PulserTrms","PulserTrms");
 pulserObjects->Add(pulserTrms);
 
 AliTPCCalPad *pulserQmean = new AliTPCCalPad("PulserQmean","PulserQmean");
 pulserObjects->Add(pulserQmean);
 


 UInt_t result=0;

 Int_t nSectors = fROC->GetNSectors();
 TVectorD foundTmean(nSectors);
 foundTmean=0;
 TVectorD foundTrms(nSectors);
 foundTrms=0;
 TVectorD foundQmean(nSectors);
 foundQmean=0;


 TList* list = GetFileSources(sourceFXS,"pulser");
 
 if (list && list->GetEntries()>0) {

//  loop through all files from LDCs

    Bool_t changed=false;
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
	f->GetObject("tpcCalibPulser",calPulser);
        if ( !calPulser ) {
	  Log ("No pulser calibration object in file.");
	  result = 2;
	  break;
	}

        //  replace entries for the sectors available in the present file

        changed=true;
        for (Int_t sector=0; sector<nSectors; sector++) {
           AliTPCCalROC *rocTmean=calPulser->GetCalRocT0(sector);
           if ( rocTmean ) {
               AliTPCCalROC* roc=pulserTmean->GetCalROC(sector);
	       roc->Add(rocTmean,1);
	       foundTmean[sector]++;
           }
           AliTPCCalROC *rocTrms=calPulser->GetCalRocRMS(sector);
           if ( rocTrms ) {
               AliTPCCalROC* roc=pulserTrms->GetCalROC(sector);
	       roc->Add(rocTrms,1);
	       foundTrms[sector]++;
           }
           AliTPCCalROC *rocQmean=calPulser->GetCalRocQ(sector);
           if ( rocQmean ) {
               AliTPCCalROC* roc=pulserQmean->GetCalROC(sector);
	       roc->Add(rocQmean,1);
	       foundQmean[sector]++;
           }
        }
       delete calPulser;
       f->Close();
      }
     ++index;
    }  // while(list)




if (pulserObjectsOCDB) {
  AliTPCCalPad* pulserTmeanOCDB = (AliTPCCalPad*)pulserObjectsOCDB->FindObject("PulserTmean");
  AliTPCCalPad* pulserTrmsOCDB  = (AliTPCCalPad*)pulserObjectsOCDB->FindObject("PulserTrms");
  AliTPCCalPad* pulserQmeanOCDB = (AliTPCCalPad*)pulserObjectsOCDB->FindObject("PulserQmean");
 
//
//  Check if calibration is complete -- otherwise take from old OCDB entry

// inner sectors
    for (Int_t sector=0; sector<nSectors/2; sector++) {
       if (foundTmean[sector] < 1 ) {
         if (pulserTmeanOCDB) {
          AliTPCCalROC* rocOCDB = pulserTmeanOCDB->GetCalROC(sector);
          if ( rocOCDB ) pulserTmean->SetCalROC(rocOCDB,sector);
         }
       }
       if (foundTrms[sector] < 1 ) {
         if (pulserTrmsOCDB) {
          AliTPCCalROC* rocOCDB = pulserTrmsOCDB->GetCalROC(sector);
          if ( rocOCDB ) pulserTrms->SetCalROC(rocOCDB,sector);
         }
       }
       if (foundQmean[sector] < 1 ) {
         if (pulserQmeanOCDB) {
          AliTPCCalROC* rocOCDB = pulserQmeanOCDB->GetCalROC(sector);
          if ( rocOCDB ) pulserQmean->SetCalROC(rocOCDB,sector);
         }
       }
    }
  
// outer sectors -- two updates needed

    for (Int_t sector=0; sector<nSectors/2; sector++) {
       if (foundTmean[sector] < 2 ) {
         if (pulserTmeanOCDB) {
          AliTPCCalROC* rocOCDB = pulserTmeanOCDB->GetCalROC(sector);
          if ( rocOCDB ) pulserTmean->SetCalROC(rocOCDB,sector);
         }
       }
       if (foundTrms[sector] < 2 ) {
         if (pulserTrmsOCDB) {
          AliTPCCalROC* rocOCDB = pulserTrmsOCDB->GetCalROC(sector);
          if ( rocOCDB ) pulserTrms->SetCalROC(rocOCDB,sector);
         }
       }
       if (foundQmean[sector] < 2 ) {
         if (pulserQmeanOCDB) {
          AliTPCCalROC* rocOCDB = pulserQmeanOCDB->GetCalROC(sector);
          if ( rocOCDB ) pulserQmean->SetCalROC(rocOCDB,sector);
         }
       }

    }
}

//
//  Store updated pedestal entry to OCDB
//
    if (changed) {
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Haavard Helstrup");
     metaData.SetAliRootVersion(ALIROOT_VERSION);
     metaData.SetComment("Preprocessor AliTPC data base entries.");

     Bool_t storeOK = Store("Calib", "Pulser", pulserObjects, &metaData, 0, kTRUE ^ fForceSingleRun);
     if ( !storeOK ) ++result;
    }  
  } else {
    Log ("Error: no entries in pulser file list!");
    result = 1;
  }
  pulserObjects->Delete();
  delete pulserObjects;
  if (pulserObjectsOCDB) {
    pulserObjectsOCDB->Delete();
  }
  delete pulserObjectsOCDB;

  return result;
}

//______________________________________________________________________________________________


UInt_t AliTPCPreprocessor::ExtractRaw(Int_t sourceFXS)
{
 //
 //  Read Raw calibration file from file exchage server
 //
 
 UInt_t result=0;
 TObjArray* rawArray = new TObjArray;
 rawArray->SetOwner(kTRUE);

 TList* list = GetFileSources(sourceFXS,"tpcCalibRaw");
 
 if (list && list->GetEntries()>0) {

//  loop through all files

    UInt_t index = 0;
    while (list->At(index)!=NULL) {
     TObjString* fileNameEntry = (TObjString*) list->At(index);
     if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "tpcCalibRaw",
	                                 fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  Log ("Error opening raw file.");
	  result =2;
	  break;
	}
        AliTPCCalibRaw *calRaw;
	f->GetObject("tpcCalibRaw",calRaw);
        if ( !calRaw ) {
	  Log ("No raw calibration object in file.");
	  result = 2;
	  break;
	}
       rawArray->Add(calRaw);
       f->Close();
      }
     ++index;
    }  // while(list)
//
//  Store updated pedestal entry to OCDB
//
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Haavard Helstrup");
     metaData.SetAliRootVersion(ALIROOT_VERSION);
     metaData.SetComment("Preprocessor AliTPC data base entries.");

     Bool_t storeOK = Store("Calib", "Raw", rawArray, &metaData, 0, kTRUE ^ fForceSingleRun);
     if ( !storeOK ) ++result;
  } else {
    Log ("Error: no entries in raw file list!");
    result = 1;
  }
  
  rawArray->Delete();
  delete rawArray;

  return result;
}
//______________________________________________________________________________________________

UInt_t AliTPCPreprocessor::ExtractCE(Int_t sourceFXS)
{
 //
 //  Read Central Electrode file from file exchage server
 //
 //
  AliTPCCalPad *ceTmean=0;
  AliTPCCalPad *ceTrms=0;
  AliTPCCalPad *ceQmean=0;
  TObjArray    *rocTtime=0;
  TObjArray    *rocQtime=0;
  
  TObjArray    *ceObjects= new TObjArray;
  ceObjects->SetOwner(kTRUE);
  
  Int_t nSectors = fROC->GetNSectors();
  
  ceTmean = new AliTPCCalPad("CETmean","CETmean");
  ceObjects->Add(ceTmean);
  
  ceTrms = new AliTPCCalPad("CETrms","CETrms");
  ceObjects->Add(ceTrms);
  
  ceQmean = new AliTPCCalPad("CEQmean","CEQmean");
  ceObjects->Add(ceQmean);
  
  rocTtime = new TObjArray(nSectors+2);   // also make room for A and C side average
  rocTtime->SetName("rocTtime");
  ceObjects->Add(rocTtime);
  rocTtime->SetOwner(kTRUE);
  
  rocQtime = new TObjArray(nSectors);
  rocQtime->SetName("rocQtime");
  ceObjects->Add(rocQtime);
  rocQtime->SetOwner(kTRUE);  

  //=== new part
  TObjArray *arrFitGraphs=new TObjArray;
  arrFitGraphs->SetName("ceFitsDrift");
  ceObjects->Add(arrFitGraphs);
  arrFitGraphs->SetOwner(kTRUE);
  
// Temperature maps 
  
  if (fTemp) {
    AliTPCSensorTempArray *tempMap = new AliTPCSensorTempArray(*fTemp);
    tempMap->SetNameTitle("TempMap","TempMap");
    ceObjects->Add(tempMap);
  }
  
// Pressure maps
  
  if (fPressure) {
    AliDCSSensor *sensor=0, *sensorCopy=0;
    for (Int_t isensor=0; isensor<kNumPressureSensors; ++isensor ) {
      sensor = fPressure->GetSensor(kPressureSensorNames[isensor]);
      if (sensor) {
        sensorCopy = new AliDCSSensor(*sensor);
        sensorCopy->SetNameTitle(kPressureSensorNames[isensor],kPressureSensorNames[isensor]);
        ceObjects->Add(sensorCopy);
      }
    }
  }
  
  UInt_t result=0;
  
  TList* list = GetFileSources(sourceFXS,"CE");
  
  if (list && list->GetEntries()>0) {
    
//  loop through all files from LDCs
    
    UInt_t index = 0;
    while (list->At(index)!=NULL) {
      TObjString* fileNameEntry = (TObjString*) list->At(index);
      if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "CE",
                                   fileNameEntry->GetString().Data());
        AliTPCCalibCE *calCE=AliTPCCalibCE::ReadFromFile(fileName.Data());        
        
        if (!calCE) {
          Log ("No valid calibCE object.");
          result=2;
          break;
        }
        //  replace entries for the sectors available in the present file
        
        for (Int_t sector=0; sector<nSectors; sector++) {
          AliTPCCalROC *rocTmean=calCE->GetCalRocT0(sector);
          if ( rocTmean )  ceTmean->SetCalROC(rocTmean,sector);
          AliTPCCalROC *rocTrms=calCE->GetCalRocRMS(sector);
          if ( rocTrms )  ceTrms->SetCalROC(rocTrms,sector);
          AliTPCCalROC *rocQmean=calCE->GetCalRocQ(sector);
          if ( rocQmean )  ceQmean->SetCalROC(rocQmean,sector);
          TGraph *grT=calCE->MakeGraphTimeCE(sector,0,2); // T time graph
          if ( grT ) rocTtime->AddAt(grT,sector);
          TGraph *grQ=calCE->MakeGraphTimeCE(sector,0,3); // Q time graph
          if ( grQ ) rocQtime->AddAt(grQ,sector);
        }
        
        TGraph *grT=calCE->MakeGraphTimeCE(-1,0,2); // A side average
        if ( grT ) {
          rocTtime->AddAt(grT,nSectors);
        } else {
          result=10;
        }
        grT=calCE->MakeGraphTimeCE(-2,0,2); // C side average
        if ( grT ) {
          rocTtime->AddAt(grT,nSectors+1);
        } else {
          result=10;
        }

    //    calCE->Delete();
        delete calCE;
      }
      ++index;
    }  // while(list)
//
//   Check number of calibrated sectors per side
// 
    Int_t aside=0, cside=0;
    for (Int_t ind=0; ind<nSectors/4; ind++ ) {
        TGraph *grT=(TGraph*)rocTtime->At(ind);
	if (grT) aside++;
	grT=(TGraph*)rocTtime->At(ind+nSectors/2);
	if (grT) aside++;
	grT=(TGraph*)rocTtime->At(ind+nSectors/4);
	if (grT) cside++;
	grT=(TGraph*)rocTtime->At(ind+3*nSectors/4);
	if (grT) cside++;
     }
     if ( (aside<kMinCESectors) && (cside<kMinCESectors) ) {
        Log (Form("ExtractCE: Too few fitted sectors: Aside =%d, Cside=%d\n",
	     aside, cside)) ;
	result=10;
     }

    //
    //=== New CE part
    //    if it is validated this part needs to be modified again
    //    currently its only processed if there is a valid standard CE object
    //
    list = GetFileSources(sourceFXS,"CEnew");
    
    if (result==0 && list && list->GetEntries()>0) {
      
//  loop through all files from LDCs
      
      UInt_t index2 = 0;
      while (list->At(index2)!=NULL) {
        TObjString* fileNameEntry = (TObjString*) list->At(index2);
        if (fileNameEntry!=NULL) {
          TString fileName = GetFile(sourceFXS, "CEnew",
                                     fileNameEntry->GetString().Data());
          AliTPCCalibCE *calCE=AliTPCCalibCE::ReadFromFile(fileName.Data());        
          
          if (!calCE) {
            Log ("No valid new calibCE object.");
//             result=2;
            break;
          }

          TIter nextObj(calCE->GetArrFitGraphs());
          TObject *obj=0x0;
          while ( (obj=nextObj()) ){
            arrFitGraphs->Add(obj->Clone());
          }
//          calCE->Delete();
          delete calCE;
        }
        ++index2;
      }
    }
    
//
//  Store updated pedestal entry to OCDB
//
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Haavard Helstrup");
    metaData.SetAliRootVersion(ALIROOT_VERSION);
    metaData.SetComment("Preprocessor AliTPC data base entries.");
    
    if ( result == 0 ) {
      Bool_t storeOK = Store("Calib", "CE", ceObjects, &metaData, 0, kTRUE ^ fForceSingleRun);
      if ( !storeOK ) ++result;
    } else {
      Log ("Warning: Average time graphs not available - no OCDB entry written");
    }
  } else {
    Log ("Error: no CE entries available from FXS!");
    result = 1;
  }
  
  ceObjects->Delete();
  delete ceObjects;
  
  return result;
}
//______________________________________________________________________________________________

UInt_t AliTPCPreprocessor::ExtractQA(Int_t sourceFXS)
{
 //
 //  Read Quality Assurance file from file exchage server
 //
 
 UInt_t result=0;

 TList* list = GetFileSources(sourceFXS,"QA");
 
 if (list && list->GetEntries()>0) {

//  only one QA object should be available!

    AliTPCdataQA *calQA;

    UInt_t nentries = list->GetEntries();  
    UInt_t index=0;
    if ( nentries > 1) Log ( "More than one QA entry. First one processed");      
    TObjString* fileNameEntry = (TObjString*) list->At(index);
    if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, "QA",
	                                 fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  Log ("Error opening QA file.");
	  result =2;          
	} else {
   	  f->GetObject("tpcCalibQA",calQA);
          if ( calQA ) {      
//
//  Store updated pedestal entry to OCDB
//
           AliCDBMetaData metaData;
           metaData.SetBeamPeriod(0);
           metaData.SetResponsible("Haavard Helstrup");
	   metaData.SetAliRootVersion(ALIROOT_VERSION);
           metaData.SetComment("Preprocessor AliTPC data base entries.");

           Bool_t storeOK = Store("Calib", "QA", calQA, &metaData, 0, kFALSE);
           if ( !storeOK ) ++result;

           delete calQA;
	  }
        }
    } else {
    Log ("Error: no QA files on FXS!");
    result = 2;
    }
  } else {
    Log ("Error: no QA entries in FXS list!");
    result = 1;
  }
  return result;
}

//______________________________________________________________________________________________


UInt_t AliTPCPreprocessor::ExtractAltro(Int_t sourceFXS, TMap* dcsMap)
{
 //
 //  Read Altro configuration file from file exchage server
 //  Keep original entry from OCDB in case no new pulser calibration is available
 //
 TObjArray    *altroObjects=0;
 AliTPCCalPad *acqStart=0;
 AliTPCCalPad *zsThr=0;
 AliTPCCalPad *acqStop=0;
 AliTPCCalPad *FPED=0;
 AliTPCCalPad *masked=0;
 AliTPCCalPad *k1=0, *k2=0, *k3=0;
 AliTPCCalPad *l1=0, *l2=0, *l3=0;
 TMap *mapRCUconfig=0;

 AliCDBEntry* entry = GetFromOCDB("Calib", "AltroConfig");
 if (entry) altroObjects = (TObjArray*)entry->GetObject();
 if ( altroObjects==NULL ) {
     Log("AliTPCPreprocsessor: No previous TPC altro calibration entry available.\n");
     altroObjects = new TObjArray;    
 }
 altroObjects->SetOwner(kTRUE);

 acqStart = (AliTPCCalPad*)altroObjects->FindObject("AcqStart");
 if ( !acqStart ) {
    acqStart = new AliTPCCalPad("AcqStart","AcqStart");
    altroObjects->Add(acqStart);
 }
 zsThr = (AliTPCCalPad*)altroObjects->FindObject("ZsThr");
 if ( !zsThr )  { 
    zsThr = new AliTPCCalPad("ZsThr","ZsThr");
    altroObjects->Add(zsThr);
 }
 FPED = (AliTPCCalPad*)altroObjects->FindObject("FPED");
 if ( !FPED )  { 
    FPED = new AliTPCCalPad("FPED","FPED");
    altroObjects->Add(FPED);
 }
 acqStop = (AliTPCCalPad*)altroObjects->FindObject("AcqStop");
 if ( !acqStop ) {
    acqStop = new AliTPCCalPad("AcqStop","AcqStop");
    altroObjects->Add(acqStop);
 }
 masked = (AliTPCCalPad*)altroObjects->FindObject("Masked");
 if ( !masked )  { 
    masked = new AliTPCCalPad("Masked","Masked");
    altroObjects->Add(masked);
 }
 k1 = (AliTPCCalPad*)altroObjects->FindObject("K1");
 if ( !k1 )  { 
    k1 = new AliTPCCalPad("K1","K1");
    altroObjects->Add(k1);
 }
 k2 = (AliTPCCalPad*)altroObjects->FindObject("K2");
 if ( !k2 )  { 
    k2 = new AliTPCCalPad("K2","K2");
    altroObjects->Add(k2);
 }
 k3 = (AliTPCCalPad*)altroObjects->FindObject("K3");
 if ( !k3 )  { 
    k3 = new AliTPCCalPad("K3","K3");
    altroObjects->Add(k3);
 }
 l1 = (AliTPCCalPad*)altroObjects->FindObject("L1");
 if ( !l1 )  { 
    l1 = new AliTPCCalPad("L1","L1");
    altroObjects->Add(l1);
 }
 l2 = (AliTPCCalPad*)altroObjects->FindObject("L2");
 if ( !l2 )  { 
    l2 = new AliTPCCalPad("L2","L2");
    altroObjects->Add(l2);
 }
 l3 = (AliTPCCalPad*)altroObjects->FindObject("L3");
 if ( !l3 )  { 
    l3 = new AliTPCCalPad("L3","L3");
    altroObjects->Add(l3);
 }
 mapRCUconfig = (TMap*)altroObjects->FindObject("RCUconfig");
 if (!mapRCUconfig) {
    mapRCUconfig = new TMap();
    mapRCUconfig->SetName("RCUconfig");
    altroObjects->Add(mapRCUconfig);
 }


 UInt_t result=0;
 TString idFXS[2]={"AltroConfigA","AltroConfigC"};

 Int_t nSectors = fROC->GetNSectors();
 Bool_t changed=false;
 if (altroObjects == 0 ) altroObjects = new TObjArray;
 altroObjects->SetOwner(kTRUE);

// extract list of active DDLs

 if (dcsMap) {
  Bool_t found; 
  TString arrDDL(kNumDDL);
  arrDDL.Append('x',kNumDDL);
  for ( Int_t iDDL = 0; iDDL<kNumDDL; iDDL++ ) {
    TString stringID = Form (kAmandaDDL.Data(),iDDL+kFirstDDL);
    TPair *pair = (TPair*)dcsMap->FindObject(stringID.Data());
    found = false;
    if ( pair ) {
        TObjArray *valueSet=(TObjArray*)pair->Value();
        if ( valueSet) { 
	  AliDCSValue *val = (AliDCSValue*)valueSet->At(0);
	  if (val) { 
	      found = val->GetBool();
	      if (found){
               arrDDL[iDDL] = '1';
              } else { 
               arrDDL[iDDL] = '0';
              }
	  }    
	}
    } 
  }
  TObjString *ddlArray = new TObjString;
  ddlArray->SetString(arrDDL);
  TMap *activeDDL = new TMap;
  activeDDL->SetName("DDLArray");
  TObjString *key = new TObjString("DDLArray");
  activeDDL->Add(key,ddlArray);
  altroObjects->Add(activeDDL);
  changed=true;
 } else {
   Log ("ExtractAltro: No DCS map available. Active DDL list cannot be obtained.");
   result = 3;
 }        

// extract Altro configuration files

 for ( Int_t id=0; id<2; id++) {
   TList* list = GetFileSources(sourceFXS,idFXS[id].Data());
 
   if (list && list->GetEntries()>0) {

//  loop through all files from LDCs

    UInt_t index = 0;
    while (list->At(index)!=NULL) {
     TObjString* fileNameEntry = (TObjString*) list->At(index);
     if (fileNameEntry!=NULL) {
        TString fileName = GetFile(sourceFXS, idFXS[id].Data(),
	                                 fileNameEntry->GetString().Data());
        TFile *f = TFile::Open(fileName);
        if (!f) {
	  const int mess_length=40;
          char message[mess_length];
	  snprintf(message,mess_length,"Error opening Altro configuration file, id = %d",id);
	  Log (message);
	  result =2;
	  break;
	}
        TObjArray *altroFXS;
	f->GetObject("AltroConfig",altroFXS);
        if ( !altroFXS ) {
	  Log ("No Altro configuration object in file.");
	  result = 2;
	  break;
	}

        //  replace entries for the sectors available in the present file
        AliTPCCalPad *acqStartFXS=(AliTPCCalPad*)altroFXS->FindObject("AcqStart");
        AliTPCCalPad *zsThrFXS=(AliTPCCalPad*)altroFXS->FindObject("ZsThr");
        AliTPCCalPad *acqStopFXS=(AliTPCCalPad*)altroFXS->FindObject("AcqStop");
        AliTPCCalPad *FPEDFXS=(AliTPCCalPad*)altroFXS->FindObject("FPED");
        AliTPCCalPad *maskedFXS=(AliTPCCalPad*)altroFXS->FindObject("Masked");
        AliTPCCalPad *k1FXS=(AliTPCCalPad*)altroFXS->FindObject("K1");
        AliTPCCalPad *k2FXS=(AliTPCCalPad*)altroFXS->FindObject("K2");
        AliTPCCalPad *k3FXS=(AliTPCCalPad*)altroFXS->FindObject("K3");
        AliTPCCalPad *l1FXS=(AliTPCCalPad*)altroFXS->FindObject("L1");
        AliTPCCalPad *l2FXS=(AliTPCCalPad*)altroFXS->FindObject("L2");
        AliTPCCalPad *l3FXS=(AliTPCCalPad*)altroFXS->FindObject("L3");
        TMap *mapRCUconfigFXS = (TMap*)altroFXS->FindObject("RCUconfig");
        TIterator *mapFXSiter = mapRCUconfigFXS->MakeIterator();
	
        changed=true;
        for (Int_t sector=0; sector<nSectors; sector++) {
            
           if (acqStartFXS) {
	      AliTPCCalROC *rocAcqStart=acqStartFXS->GetCalROC(sector);
              if ( rocAcqStart )  acqStart->SetCalROC(rocAcqStart,sector);
	   }
	   if (zsThrFXS ) {
              AliTPCCalROC *rocZsThr=zsThrFXS->GetCalROC(sector);
              if ( rocZsThr )  zsThr->SetCalROC(rocZsThr,sector);
	   }
	   if (acqStopFXS) {
              AliTPCCalROC *rocAcqStop=acqStopFXS->GetCalROC(sector);
              if ( rocAcqStop )  acqStop->SetCalROC(rocAcqStop,sector);
	   }
	   if (FPEDFXS ) {
              AliTPCCalROC *rocFPED=FPEDFXS->GetCalROC(sector);
              if ( rocFPED )  FPED->SetCalROC(rocFPED,sector);
	   }
	   if (maskedFXS) {
              AliTPCCalROC *rocMasked=maskedFXS->GetCalROC(sector);
              if ( rocMasked )  masked->SetCalROC(rocMasked,sector);
	   }
	   if (k1FXS) {
              AliTPCCalROC *rocK1=k1FXS->GetCalROC(sector);
              if ( rocK1 )  k1->SetCalROC(rocK1,sector);
	   }
	   if (k2FXS) {
              AliTPCCalROC *rocK2=k2FXS->GetCalROC(sector);
              if ( rocK2 )  k2->SetCalROC(rocK2,sector);
	   }
	   if (k3FXS) {
              AliTPCCalROC *rocK3=k3FXS->GetCalROC(sector);
              if ( rocK3 )  k3->SetCalROC(rocK3,sector);
	   }
	   if (l1FXS) {
              AliTPCCalROC *rocL1=l1FXS->GetCalROC(sector);
              if ( rocL1 )  l1->SetCalROC(rocL1,sector);
	   }
	   if (l2FXS) {
              AliTPCCalROC *rocL2=l2FXS->GetCalROC(sector);
              if ( rocL2 )  l2->SetCalROC(rocL2,sector);
	   }
	   if (l3FXS) {
              AliTPCCalROC *rocL3=l3FXS->GetCalROC(sector);
              if ( rocL3 )  l3->SetCalROC(rocL3,sector);
	   }
	 }
	 if (mapRCUconfigFXS) {
          Int_t mapEntries = mapRCUconfigFXS->GetEntries();
          TObjString* keyFXS;
	  TVectorF* vecFXS;
	  TVectorF* vec;              // nSectors = 72  (total number of inner/outer sectors)
	  for (Int_t i=0; i<mapEntries; ++i) {
	    keyFXS=(TObjString*)mapFXSiter->Next();
            vecFXS=(TVectorF*)mapRCUconfigFXS->GetValue(keyFXS);
            vec=(TVectorF*)mapRCUconfig->GetValue(keyFXS);
	    if (!vec) {
	      vec = new TVectorF(3*nSectors);
	      *vec = -1;
	      mapRCUconfig->Add(keyFXS,vec);
	    }
	    if (vec->GetNoElements() != 3*nSectors ) {
	      vec->ResizeTo(3*nSectors);
            }
	    if (id==0) {                        // A side
	      vec->SetSub(0,vecFXS->GetSub(0,nSectors/2-1));
	      vec->SetSub(nSectors,vecFXS->GetSub(nSectors,2*nSectors-1));
	    } else {                             // C side
	      vec->SetSub(nSectors/2,vecFXS->GetSub(nSectors/2,nSectors-1));
	      vec->SetSub(2*nSectors,vecFXS->GetSub(2*nSectors,3*nSectors-1));
	    }
	  }
        }
       altroFXS->Delete();
       delete altroFXS;
       f->Close();
      }
     ++index;
     }  // while(list)
    } else {
      Log ("Error: no entries in AltroConfig file list!");
      result = 1;
    }

   }   // for - id
//
//  Store updated pedestal entry to OCDB
//
    if (changed) {
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Haavard Helstrup");
     metaData.SetAliRootVersion(ALIROOT_VERSION);
     metaData.SetComment("Preprocessor AliTPC data base entries.");

     Bool_t storeOK = Store("Calib", "AltroConfig", altroObjects, &metaData, 0, kFALSE);
     if ( !storeOK ) ++result;
    }  
    
  altroObjects->Delete();
  delete altroObjects;
  
  return result;
}
