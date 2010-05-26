////////////////////////////////////////////
// Author: Ionut Cristian Arsene          //          
// email:  iarsene@mail.cern.ch           //
////////////////////////////////////////////
// Use this macro to create ROOT trees with time dependent information from the TRD OCDB
//
// Usage:
//
//DumpOCDBtoTree(const Char_t* runListFilename, const Char_t* outFilename,
//    	         Int_t firstRun = -1, Int_t lastRun = -1,
//		 Bool_t getMonitoringInfo = kTRUE,
//               Bool_t getCalibrationInfo = kTRUE,
// 	         Bool_t getGoofieInfo = kTRUE,
//		 Bool_t getDCSInfo = kFALSE,
//		 const Char_t* storage = "local:///lustre/alice/alien/alice/data/2010/OCDB/")
//
//    * runListFilename   - name of an ascii file containing run numbers
//    * outFilename       - name of the root file where the TRD OCDB information tree to be stored
//    * firstRun, lastRun - lowest and highest run numbers (from the ascii file) to be dumped
//                          if these numbers are not specified (-1) all run numbers in the input ascii file will
//                          be used. If the run list file is not specified then all runs in this interval
//                          will be queried
//    * getMonitoringInfo - flag to switch on/off monitoring information (HV, temperatures, gas pressure,
//                          gas composition, ADC clock phase, chamber status)
//    * getCalibrationInfo- flag to switch on/off calibration information (gain, pedestal, T0, vdrift, pad status)
//    * getGoofieInfo     - flag to switch on/off goofie information (gain, HV, pressure, temperature, drift velocity,
//                          gas composition)
//    * getDCSInfo        - flag to switch on/off DCS information ()
//    * storage           - path of the OCDB database (if it is on alien, be sure to have a valid/active token)


#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TObjString.h"
#include "TTimeStamp.h"
#include "TH1.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TSystem.h"
#include "TGrid.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTRDcalibDB.h"
#include "AliGRPObject.h"
#include "AliDCSSensor.h"
#include "AliTRDSensorArray.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalROC.h"
#include "AliTRDCalPadStatus.h"
#include "AliTRDCalChamberStatus.h"
#include "AliTRDCalSingleChamberStatus.h"
#include "AliTRDCalDCS.h"
#include "AliTRDCalDCSFEE.h"
using namespace std;

// global variables
// histograms used for extracting the mean and RMS of calibration parameters
TH1F *gRunWiseHisto;
TH1F *gSuperModuleWiseHisto;
TH1F *gChamberWiseHisto;

// global constants
const Int_t gkSuperModuleStatus[18] = {1, 1, 0, 0, 0, 0, 0, 1, 1,   // (1-installed)
				       1, 1, 0, 0, 0, 0, 0, 0, 1};  

void MakeRunListFromOCDB(const Char_t* directory, const Char_t* outfile, Bool_t fromAlien=kFALSE);
void ProcessTRDSensorArray(AliTRDSensorArray*, TTimeStamp, TVectorD&);
void ProcessTRDCalibArray(AliTRDCalDet*, AliTRDCalPad*, TString, Double_t&, Double_t&,
			  TVectorD&, TVectorD&, TVectorD&, TVectorD&);
void ProcessTRDstatus(AliTRDCalChamberStatus*, AliTRDCalPadStatus*, Float_t&, TVectorD&, TVectorD&);
void ProcessTRDCalDCSFEE(AliTRDCalDCS*, AliTRDCalDCS*, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Bool_t&, 
			 TVectorD&, TVectorD&);

//__________________________________________________________________________________________
void DumpOCDBtoTree(const Char_t* runListFilename,
		    const Char_t* outFilename,
		    Int_t firstRun = -1, Int_t lastRun = -1,
		    Bool_t getMonitoringInfo = kTRUE,
                    Bool_t getCalibrationInfo = kTRUE,
		    Bool_t getGoofieInfo = kTRUE,
		    Bool_t getDCSInfo = kFALSE,
		    const Char_t* storage = "alien://folder=/alice/data/2010/OCDB/") {
  //
  // Main function to steer the extraction of TRD OCDB information
  //


  TTimeStamp jobStartTime;
  // if the storage is on alien than we need to do some extra stuff
  TString storageString(storage);
  if(storageString.Contains("alien://")) {
    TGrid::Connect("alien://");
  }
  // initialize the OCDB manager
  AliCDBManager *manager = AliCDBManager::Instance();
  manager->SetDefaultStorage(storage);
    
  // initialize the tree
  TTreeSRedirector *treeStreamer = new TTreeSRedirector(outFilename);
  
  // initialize the histograms used for extracting the mean and RMS
  gRunWiseHisto = new TH1F("runHisto", "runHisto", 200, -10.0, 10.0);
  gSuperModuleWiseHisto = new TH1F("smHisto", "smHisto", 200, -10.0, 10.0);
  gChamberWiseHisto = new TH1F("chamberHisto", "chamberHisto", 200, -10.0, 10.0);

  // open the ascii file with run numbers
  ifstream in;
  if(runListFilename[0]!='\0')
    in.open(runListFilename);

  // if a run list file was not specified then use the run range
  Int_t currRun;
  if(runListFilename[0]=='\0' && firstRun!=-1 && lastRun!=-1)
    currRun = firstRun-1;  

  TVectorD runs;
  TVectorD rejectedRuns;
  TTimeStamp loopStartTime;

  // loop over runs
  while(1) {
    // check if we still have run numbers in the file or provided range
    if(runListFilename[0]=='\0' && firstRun!=-1 && lastRun!=-1) {
      currRun++;
      if(currRun>lastRun) break;
    }
    if(runListFilename[0]!='\0') {
      if(in.eof()) break;
      if(!(in>>currRun)) continue;
      if(currRun < (firstRun==-1 ? 0 : firstRun) ||
	 currRun > (lastRun==-1 ? 999999999 : lastRun))
	continue;
    }

    cout << "run = " << currRun << endl;
    // check if the run was processed already
    Bool_t runProcessed = kFALSE;
    for(Int_t iRun=0; iRun<runs.GetNoElements(); iRun++) {
      if(runs[iRun]==currRun) runProcessed = kTRUE;
    }
    if(runProcessed) {
      cout << "Run processed already" << endl;
      continue;
    }

    manager->SetRun(currRun);

    // Get the GRP data. Only runs with a corresponding GRP entry in the OCDB 
    // will be processed.
    time_t startTime = 0;
    time_t endTime = 0;
    TObjString runType("UNKNOWN");
    AliDCSSensor *cavern_pressure;
    AliDCSSensor *surface_pressure;
    UInt_t detectorMask = 0;
    AliCDBEntry *entry = manager->Get("GRP/GRP/Data");
    if(entry) {
      entry->SetOwner(kFALSE);
    }
    AliGRPObject* grpObject = 0;
    if(entry) 
      grpObject = dynamic_cast<AliGRPObject*>(entry->GetObject());
    else {
      // add the run number to the list of rejected runs
      rejectedRuns.ResizeTo(rejectedRuns.GetNoElements()+1);
      rejectedRuns[rejectedRuns.GetNoElements()-1] = currRun;
      continue;
    }
    if(grpObject) {
      startTime = grpObject->GetTimeStart();
      endTime = grpObject->GetTimeEnd();
      runType = grpObject->GetRunType().Data();
      cavern_pressure = grpObject->GetCavernAtmosPressure();
      surface_pressure = grpObject->GetSurfaceAtmosPressure();
      detectorMask = grpObject->GetDetectorMask();
      TTimeStamp start(grpObject->GetTimeStart());
      TTimeStamp end(grpObject->GetTimeEnd());
      cout << "Start time: " << start.GetDate()/10000 << "/" 
	   << (start.GetDate()/100)-(start.GetDate()/10000)*100 << "/" 
	   << start.GetDate()%100 << "   "
	   << start.GetTime()/10000 << ":"
	   << (start.GetTime()/100)-(start.GetTime()/10000)*100 << ":" 
	   << start.GetTime()%100 << endl;
      cout << "End time: " << end.GetDate()/10000 << "/" 
	   << (end.GetDate()/100)-(end.GetDate()/10000)*100 << "/" 
	   << end.GetDate()%100 << "   "
	   << end.GetTime()/10000 << ":"
	   << (end.GetTime()/100)-(end.GetTime()/10000)*100 << ":"
	   << end.GetTime()%100 << endl;
      cout << "Run type = " << grpObject->GetRunType().Data() << endl;
    }
    else {
      // add the run number to the list of rejected runs
      rejectedRuns.ResizeTo(rejectedRuns.GetNoElements()+1);
      rejectedRuns[rejectedRuns.GetNoElements()-1] = currRun;
      continue;
    }

    // remove runs with zero time duration
    if(startTime==endTime) {
      if(grpObject) delete grpObject;
      // add the run number to the list of rejected runs
      rejectedRuns.ResizeTo(rejectedRuns.GetNoElements()+1);
      rejectedRuns[rejectedRuns.GetNoElements()-1] = currRun;
      continue;
    }

    // time step for time dependent information (change this if you need something else)
    UInt_t dTime = TMath::Max((endTime-startTime)/20, Long_t(5*60));

    // get monitoring information
    AliTRDSensorArray *anodeISensors = 0;
    AliTRDSensorArray *anodeUSensors = 0;
    AliTRDSensorArray *driftISensors = 0;
    AliTRDSensorArray *driftUSensors = 0;
    AliTRDSensorArray *temperatureSensors = 0;
    AliTRDSensorArray *chamberStatusSensors = 0;
    AliTRDSensorArray *overpressureSensors = 0;
    AliTRDSensorArray *gasCO2Sensors = 0;
    AliTRDSensorArray *gasH2OSensors = 0;
    AliTRDSensorArray *gasO2Sensors = 0;
  //  AliTRDSensorArray *adcClkPhaseSensors = 0;

    if(getMonitoringInfo) {
      // anode hv currents (per chamber)
      entry = manager->Get("TRD/Calib/trd_hvAnodeImon");
      if(entry) {
	entry->SetOwner(kTRUE);
	anodeISensors = (AliTRDSensorArray*)entry->GetObject();
      }

      // anode hv voltages (per chamber)
      entry = manager->Get("TRD/Calib/trd_hvAnodeUmon");
      if(entry) {
	entry->SetOwner(kTRUE);
	anodeUSensors = (AliTRDSensorArray*)entry->GetObject();
      }

      // drift hv currents (per chamber)
      entry = manager->Get("TRD/Calib/trd_hvDriftImon");
      if(entry) {
	entry->SetOwner(kTRUE);
	driftISensors = (AliTRDSensorArray*)entry->GetObject();
      }
    
      // drift hv voltages (per chamber)
      entry = manager->Get("TRD/Calib/trd_hvDriftUmon");
      if(entry) {
	entry->SetOwner(kTRUE);
	driftUSensors = (AliTRDSensorArray*)entry->GetObject();
      }
    
      // temperatures from chamber sensors (per chamber)
      entry = manager->Get("TRD/Calib/trd_envTemp");
      if(entry) {
	entry->SetOwner(kTRUE);
	temperatureSensors = (AliTRDSensorArray*)entry->GetObject();
      }
    
      // chamber status (from sensors)
      entry = manager->Get("TRD/Calib/trd_chamberStatus");
      if(entry) {
	entry->SetOwner(kTRUE);
	chamberStatusSensors = (AliTRDSensorArray*)entry->GetObject();
      }
    
      // gas overpressure (per whole TRD)
      entry = manager->Get("TRD/Calib/trd_gasOverpressure");
      if(entry) {
	entry->SetOwner(kTRUE);
	overpressureSensors = (AliTRDSensorArray*)entry->GetObject();
      }
    
      // gas CO2 fraction (whole TRD)
      entry = manager->Get("TRD/Calib/trd_gasCO2");
      if(entry) {
	entry->SetOwner(kTRUE);
	gasCO2Sensors = (AliTRDSensorArray*)entry->GetObject();
      }
    
      // gas H2O fraction (whole TRD)
      entry = manager->Get("TRD/Calib/trd_gasH2O");
      if(entry) {
	entry->SetOwner(kTRUE);
	gasH2OSensors = (AliTRDSensorArray*)entry->GetObject();
      }
      
      // gas O2 fraction (whole TRD)
      entry = manager->Get("TRD/Calib/trd_gasO2");
      if(entry) {
	entry->SetOwner(kTRUE);
	gasO2Sensors = (AliTRDSensorArray*)entry->GetObject();
      }
    
      // ADC Clk phase (whole TRD)
/*
      entry = manager->Get("TRD/Calib/trd_adcClkPhase");
      if(entry) {
	entry->SetOwner(kTRUE);
	adcClkPhaseSensors = (AliTRDSensorArray*)entry->GetObject();
      }
*/
    }  // end if getMonitoringInfo


    // get calibration information
    // process gains
    AliTRDCalDet *chamberGainFactor = 0;
    AliTRDCalPad *padGainFactor = 0;
    Double_t runMeanGain=0.0, runRMSGain=0.0;
    TVectorD chamberMeanGain(AliTRDcalibDB::kNdet);
    TVectorD chamberRMSGain(AliTRDcalibDB::kNdet);
    TVectorD smMeanGain(AliTRDcalibDB::kNsector);
    TVectorD smRMSGain(AliTRDcalibDB::kNsector);
    TString parName("Gain");
    if(getCalibrationInfo) {
      entry = manager->Get("TRD/Calib/ChamberGainFactor");
      if(entry) {
	entry->SetOwner(kFALSE);
	chamberGainFactor = (AliTRDCalDet*)entry->GetObject();
      }
      entry = manager->Get("TRD/Calib/LocalGainFactor");
      if(entry) {
	entry->SetOwner(kFALSE);
	padGainFactor = (AliTRDCalPad*)entry->GetObject();
      }
      ProcessTRDCalibArray(chamberGainFactor, padGainFactor, 
			   parName,
			   runMeanGain, runRMSGain,
			   chamberMeanGain, chamberRMSGain,
			   smMeanGain, smRMSGain);
    }

    // process pedestals
    AliTRDCalDet *chamberNoise = 0;
    AliTRDCalPad *padNoise = 0;
    Double_t runMeanNoise=0.0, runRMSNoise=0.0;
    TVectorD chamberMeanNoise(AliTRDcalibDB::kNdet);
    TVectorD chamberRMSNoise(AliTRDcalibDB::kNdet);
    TVectorD smMeanNoise(AliTRDcalibDB::kNsector);
    TVectorD smRMSNoise(AliTRDcalibDB::kNsector);
    parName = "Noise";
    if(getCalibrationInfo) {
      entry = manager->Get("TRD/Calib/DetNoise");
      if(entry) {
	entry->SetOwner(kFALSE);
	chamberNoise = (AliTRDCalDet*)entry->GetObject();
      }
      entry = manager->Get("TRD/Calib/PadNoise");
      if(entry) {
	entry->SetOwner(kFALSE);
	padNoise = (AliTRDCalPad*)entry->GetObject();
      }
      ProcessTRDCalibArray(chamberNoise, padNoise, 
			   parName,
			   runMeanNoise, runRMSNoise,
			   chamberMeanNoise, chamberRMSNoise,
			   smMeanNoise, smRMSNoise);
    }

    // process drift velocity
    AliTRDCalDet *chamberVdrift = 0;
    AliTRDCalPad *padVdrift = 0;
    Double_t runMeanVdrift=0.0, runRMSVdrift=0.0;
    TVectorD chamberMeanVdrift(AliTRDcalibDB::kNdet);
    TVectorD chamberRMSVdrift(AliTRDcalibDB::kNdet);
    TVectorD smMeanVdrift(AliTRDcalibDB::kNsector);
    TVectorD smRMSVdrift(AliTRDcalibDB::kNsector);
    parName = "Vdrift";
    if(getCalibrationInfo) {
      entry = manager->Get("TRD/Calib/ChamberVdrift");
      if(entry) {
	entry->SetOwner(kFALSE);
	chamberVdrift = (AliTRDCalDet*)entry->GetObject();
      }
      entry = manager->Get("TRD/Calib/LocalVdrift");
      if(entry) {
	entry->SetOwner(kFALSE);
	padVdrift = (AliTRDCalPad*)entry->GetObject();
      }
      ProcessTRDCalibArray(chamberVdrift, padVdrift, 
			   parName,
			   runMeanVdrift, runRMSVdrift,
			   chamberMeanVdrift, chamberRMSVdrift,
			   smMeanVdrift, smRMSVdrift);
    }
    
    // process T0
    AliTRDCalDet *chamberT0 = 0;
    AliTRDCalPad *padT0 = 0;
    Double_t runMeanT0=0.0, runRMST0=0.0;
    TVectorD chamberMeanT0(AliTRDcalibDB::kNdet);
    TVectorD chamberRMST0(AliTRDcalibDB::kNdet);
    TVectorD smMeanT0(AliTRDcalibDB::kNsector);
    TVectorD smRMST0(AliTRDcalibDB::kNsector);
    parName = "T0";
    if(getCalibrationInfo) {
      entry = manager->Get("TRD/Calib/ChamberT0");
      if(entry) {
	entry->SetOwner(kFALSE);
	chamberT0 = (AliTRDCalDet*)entry->GetObject();
      }
      entry = manager->Get("TRD/Calib/LocalT0");
      if(entry) {
	entry->SetOwner(kFALSE);
	padT0 = (AliTRDCalPad*)entry->GetObject();
      }
      ProcessTRDCalibArray(chamberT0, padT0, 
			   parName,
			   runMeanT0, runRMST0,
			   chamberMeanT0, chamberRMST0,
			   smMeanT0, smRMST0);
    }    

    // process pad and chamber status
    AliTRDCalChamberStatus* chamberStatus = 0;
    AliTRDCalPadStatus *padStatus = 0;
    Float_t runBadPadFraction=0.0;
    TVectorD chamberBadPadFraction(AliTRDcalibDB::kNdet);
    TVectorD chamberStatusValues(AliTRDcalibDB::kNdet);
    if(getCalibrationInfo) {
      entry = manager->Get("TRD/Calib/ChamberStatus");
      if(entry) {
	entry->SetOwner(kFALSE);
	chamberStatus = (AliTRDCalChamberStatus*)entry->GetObject();
      }
      entry = manager->Get("TRD/Calib/PadStatus");
      if(entry) {
	entry->SetOwner(kFALSE);
	padStatus = (AliTRDCalPadStatus*)entry->GetObject();
      }
      ProcessTRDstatus(chamberStatus, padStatus, 
		       runBadPadFraction, chamberBadPadFraction,
		       chamberStatusValues);
    }

    // get Goofie information
    AliTRDSensorArray *goofieGainSensors = 0x0;
    AliTRDSensorArray *goofieHvSensors = 0x0;
    AliTRDSensorArray *goofiePressureSensors = 0x0;
    AliTRDSensorArray *goofieTempSensors = 0x0;
    AliTRDSensorArray *goofieVelocitySensors = 0x0;
    AliTRDSensorArray *goofieCO2Sensors = 0x0;
    AliTRDSensorArray *goofieN2Sensors = 0x0;

    if(getGoofieInfo) {
      // goofie gain
      entry = manager->Get("TRD/Calib/trd_goofieGain");
      if(entry) {
	entry->SetOwner(kTRUE);
	goofieGainSensors = (AliTRDSensorArray*)entry->GetObject();
      }
      // goofie HV
      entry = manager->Get("TRD/Calib/trd_goofieHv");
      if(entry) {
	entry->SetOwner(kTRUE);
	goofieHvSensors = (AliTRDSensorArray*)entry->GetObject();
      }
      // goofie pressure
      entry = manager->Get("TRD/Calib/trd_goofiePressure");
      if(entry) {
	entry->SetOwner(kTRUE);
	goofiePressureSensors = (AliTRDSensorArray*)entry->GetObject();
      }
      // goofie temperature
      entry = manager->Get("TRD/Calib/trd_goofieTemp");
      if(entry) {
	entry->SetOwner(kTRUE);
	goofieTempSensors = (AliTRDSensorArray*)entry->GetObject();
      }
      // goofie drift velocity
      entry = manager->Get("TRD/Calib/trd_goofieVelocity");
      if(entry) {
	entry->SetOwner(kTRUE);
	goofieVelocitySensors = (AliTRDSensorArray*)entry->GetObject();
      }
      // goofie CO2
      entry = manager->Get("TRD/Calib/trd_goofieCO2");
      if(entry) {
	entry->SetOwner(kTRUE);
	goofieCO2Sensors = (AliTRDSensorArray*)entry->GetObject();
      }
      // goofie N2
      entry = manager->Get("TRD/Calib/trd_goofieN2");
      if(entry) {
	entry->SetOwner(kTRUE);
	goofieN2Sensors = (AliTRDSensorArray*)entry->GetObject();
      }
    }   // end if getGoofieInfo

    // process the DCS FEE arrays
    Int_t nSB1 = 0; Int_t nSB2 = 0; Int_t nSB3 = 0; Int_t nSB4 = 0; Int_t nSB5 = 0; 
    Int_t nChanged = 0;
    Bool_t sorAndEor = kFALSE;
    TVectorD statusArraySOR(AliTRDcalibDB::kNdet);
    TVectorD statusArrayEOR(AliTRDcalibDB::kNdet);
    Int_t dcsFeeGlobalNTimeBins = -1;
    Int_t dcsFeeGlobalConfigTag = -1;
    Int_t dcsFeeGlobalSingleHitThres = -1;
    Int_t dcsFeeGlobalThreePadClustThres = -1;
    Int_t dcsFeeGlobalSelectiveNoSZ = -1;
    Int_t dcsFeeGlobalTCFilterWeight = -1;
    Int_t dcsFeeGlobalTCFilterShortDecPar = -1;
    Int_t dcsFeeGlobalTCFilterLongDecPar = -1;
    Int_t dcsFeeGlobalModeFastStatNoise = -1;
    TObjString dcsFeeGlobalConfigVersion("");
    TObjString dcsFeeGlobalConfigName("");
    TObjString dcsFeeGlobalFilterType("");
    TObjString dcsFeeGlobalReadoutParam("");
    TObjString dcsFeeGlobalTestPattern("");
    TObjString dcsFeeGlobalTrackletMode("");
    TObjString dcsFeeGlobalTrackletDef("");
    TObjString dcsFeeGlobalTriggerSetup("");
    TObjString dcsFeeGlobalAddOptions("");
    if(getDCSInfo) {
      TObjArray *objArrayCDB = 0;
      AliTRDCalDCS* calDCSsor = 0x0;
      AliTRDCalDCS* calDCSeor = 0x0;
      entry = manager->Get("TRD/Calib/DCS");
      if(entry) {
	entry->SetOwner(kTRUE);
	objArrayCDB = (TObjArray*)entry->GetObject();
	if(objArrayCDB) {
	  objArrayCDB->SetOwner(kTRUE);
	  calDCSsor = (AliTRDCalDCS*)objArrayCDB->At(0);
	  calDCSeor = (AliTRDCalDCS*)objArrayCDB->At(1);
	}
      }
      ProcessTRDCalDCSFEE(calDCSsor, calDCSeor, 
			  nSB1, nSB2, nSB3, nSB4, nSB5, 
			  nChanged, sorAndEor, statusArraySOR, statusArrayEOR);
      
      if(calDCSsor || calDCSeor) {
	AliTRDCalDCS *caldcs = 0;
	if(calDCSsor) caldcs = calDCSsor;
	else caldcs = calDCSeor;
	dcsFeeGlobalNTimeBins = caldcs->GetGlobalNumberOfTimeBins();
	dcsFeeGlobalConfigTag = caldcs->GetGlobalConfigTag();
	dcsFeeGlobalSingleHitThres = caldcs->GetGlobalSingleHitThres();
	dcsFeeGlobalThreePadClustThres = caldcs->GetGlobalThreePadClustThres();
	dcsFeeGlobalSelectiveNoSZ = caldcs->GetGlobalSelectiveNoZS();
	dcsFeeGlobalTCFilterWeight = caldcs->GetGlobalTCFilterWeight();
	dcsFeeGlobalTCFilterShortDecPar = caldcs->GetGlobalTCFilterShortDecPar();
	dcsFeeGlobalTCFilterLongDecPar = caldcs->GetGlobalTCFilterLongDecPar();
	dcsFeeGlobalModeFastStatNoise = caldcs->GetGlobalModeFastStatNoise();
	dcsFeeGlobalConfigVersion = caldcs->GetGlobalConfigVersion().Data();
	dcsFeeGlobalConfigName = caldcs->GetGlobalConfigName().Data();
	dcsFeeGlobalFilterType = caldcs->GetGlobalFilterType().Data();
	dcsFeeGlobalReadoutParam = caldcs->GetGlobalReadoutParam().Data();
	dcsFeeGlobalTestPattern = caldcs->GetGlobalTestPattern().Data();
	dcsFeeGlobalTrackletMode = caldcs->GetGlobalTrackletMode().Data();
	dcsFeeGlobalTrackletDef = caldcs->GetGlobalTrackletDef().Data();
	dcsFeeGlobalTriggerSetup = caldcs->GetGlobalTriggerSetup().Data();
	dcsFeeGlobalAddOptions = caldcs->GetGlobalAddOptions().Data();
      }
      if(objArrayCDB) objArrayCDB->RemoveAll();
    }   // end if(getDCSInfo)
       

    // loop over time steps
    for(UInt_t iTime = startTime; iTime<=endTime; iTime += dTime) {
      // time stamp
      TTimeStamp iStamp(iTime);
      cout << "time step  " << iStamp.GetDate()/10000 << "/" 
	   << (iStamp.GetDate()/100)-(iStamp.GetDate()/10000)*100 << "/" 
	   << iStamp.GetDate()%100 << "   "
	   << iStamp.GetTime()/10000 << ":"
	   << (iStamp.GetTime()/100)-(iStamp.GetTime()/10000)*100 << ":" 
	   << iStamp.GetTime()%100 << endl;
      
      // cavern pressure
      Float_t pressure = -99.;
      if(cavern_pressure) 
	pressure = cavern_pressure->Eval(iStamp);
            
      // surface pressure
      Float_t surfacePressure = -99.;
      if(surface_pressure) 
	surfacePressure = surface_pressure->Eval(iStamp);
            
      // anode I sensors
      TVectorD anodeIValues(AliTRDcalibDB::kNdet);
      if(anodeISensors) 
	ProcessTRDSensorArray(anodeISensors, iStamp, anodeIValues);
            
      // anode U sensors
      TVectorD anodeUValues(AliTRDcalibDB::kNdet);
      if(anodeUSensors) 
	ProcessTRDSensorArray(anodeUSensors, iStamp, anodeUValues);
            
      // drift I sensors
      TVectorD driftIValues(AliTRDcalibDB::kNdet);
      if(driftISensors) 
	ProcessTRDSensorArray(driftISensors, iStamp, driftIValues);
      
      // drift U sensors
      TVectorD driftUValues(AliTRDcalibDB::kNdet);
      if(driftUSensors) 
	ProcessTRDSensorArray(driftUSensors, iStamp, driftUValues);
      
      // chamber temperatures
      TVectorD envTempValues(AliTRDcalibDB::kNdet);
      if(temperatureSensors) 
	ProcessTRDSensorArray(temperatureSensors, iStamp, envTempValues);

      // chamber status sensors
      TVectorD statusValues(AliTRDcalibDB::kNdet);
      if(chamberStatusSensors)
	ProcessTRDSensorArray(chamberStatusSensors, iStamp, statusValues);

      // gas overpressure
      TVectorD overpressureValues(overpressureSensors ? overpressureSensors->NumSensors() : 0);
      if(overpressureSensors)
	ProcessTRDSensorArray(overpressureSensors, iStamp, overpressureValues);
      
      // gas CO2
      TVectorD gasCO2Values(gasCO2Sensors ? gasCO2Sensors->NumSensors() : 0);
      if(gasCO2Sensors)
	ProcessTRDSensorArray(gasCO2Sensors, iStamp, gasCO2Values);
      
      // gas H2O
      TVectorD gasH2OValues(gasH2OSensors ? gasH2OSensors->NumSensors() : 0);
      if(gasH2OSensors)
	ProcessTRDSensorArray(gasH2OSensors, iStamp, gasH2OValues);
      
      // gas O2
      TVectorD gasO2Values(gasO2Sensors ? gasO2Sensors->NumSensors() : 0);
      if(gasO2Sensors)
	ProcessTRDSensorArray(gasO2Sensors, iStamp, gasO2Values);
      
      // ADC Clk phase
      //TVectorD adcClkPhaseValues(adcClkPhaseSensors ? adcClkPhaseSensors->NumSensors() : 0);
      //if(adcClkPhaseSensors)
//	ProcessTRDSensorArray(adcClkPhaseSensors, iStamp, adcClkPhaseValues);

      // goofie gain
      TVectorD goofieGainValues(goofieGainSensors ? goofieGainSensors->NumSensors() : 0);
      if(goofieGainSensors)
	ProcessTRDSensorArray(goofieGainSensors, iStamp, goofieGainValues);

      // goofie HV
      TVectorD goofieHvValues(goofieHvSensors ? goofieHvSensors->NumSensors() : 0);
      if(goofieHvSensors)
	ProcessTRDSensorArray(goofieHvSensors, iStamp, goofieHvValues);

      // goofie pressure
      TVectorD goofiePressureValues(goofiePressureSensors ? goofiePressureSensors->NumSensors() : 0);
      if(goofiePressureSensors)
	ProcessTRDSensorArray(goofiePressureSensors, iStamp, goofiePressureValues);
      
      // goofie temperature
      TVectorD goofieTempValues(goofieTempSensors ? goofieTempSensors->NumSensors() : 0);
      if(goofieTempSensors) 
	ProcessTRDSensorArray(goofieTempSensors, iStamp, goofieTempValues);
      
      // goofie drift velocity
      TVectorD goofieVelocityValues(goofieVelocitySensors ? goofieVelocitySensors->NumSensors() : 0);
      if(goofieVelocitySensors) 
	ProcessTRDSensorArray(goofieVelocitySensors, iStamp, goofieVelocityValues);
      
      // goofie CO2
      TVectorD goofieCO2Values(goofieCO2Sensors ? goofieCO2Sensors->NumSensors() : 0);
      if(goofieCO2Sensors) 
	ProcessTRDSensorArray(goofieCO2Sensors, iStamp, goofieCO2Values);
      
      // goofie N2
      TVectorD goofieN2Values(goofieN2Sensors ? goofieN2Sensors->NumSensors() : 0);
      if(goofieN2Sensors) 
	ProcessTRDSensorArray(goofieN2Sensors, iStamp, goofieN2Values);
      
            
      // fill the tree
      (*treeStreamer)<< "trdTree"
		     << "run=" << currRun
		     << "time=" << iTime
		     << "startTimeGRP=" << startTime
		     << "endTimeGRP=" << endTime
		     << "runType.=" << &runType
		     << "cavernPressure=" << pressure
		     << "surfacePressure=" << surfacePressure
		     << "detectorMask=" << detectorMask;
      if(getMonitoringInfo) {
	(*treeStreamer)<< "trdTree"
		       << "hvAnodeI.=" << &anodeIValues
		       << "hvAnodeU.=" << &anodeUValues
		       << "hvDriftI.=" << &driftIValues
		       << "hvDriftU.=" << &driftUValues
		       << "envTemp.=" << &envTempValues
		       << "sensorStatusValues.=" << &statusValues
		       << "gasOverPressure.=" << &overpressureValues
		       << "gasCO2.=" << &gasCO2Values
		       << "gasH2O.=" << &gasH2OValues
		       << "gasO2.=" << &gasO2Values;
		       //<< "adcClkPhase.=" << &adcClkPhaseValues;
      }
      if(getGoofieInfo) {
	(*treeStreamer)<< "trdTree"
		       << "goofieGain.=" << &goofieGainValues
		       << "goofieHV.=" << &goofieHvValues
		       << "goofiePressure.=" << &goofiePressureValues
		       << "goofieTemp.=" << &goofieTempValues
		       << "goofieVelocity.=" << &goofieVelocityValues
		       << "goofieCO2.=" << &goofieCO2Values
		       << "goofieN2.=" << &goofieN2Values;
      }
      if(getCalibrationInfo) {
	(*treeStreamer)<< "trdTree"
		       << "runMeanGain=" << runMeanGain
		       << "runRMSGain=" << runRMSGain
		       << "smMeanGain.=" << &smMeanGain
		       << "smRMSGain.=" << &smRMSGain
		       << "chamberMeanGain.=" << &chamberMeanGain
		       << "chamberRMSGain.=" << &chamberRMSGain
		       << "runMeanNoise=" << runMeanNoise
		       << "runRMSNoise=" << runRMSNoise
		       << "smMeanNoise.=" << &smMeanNoise
		       << "smRMSNoise.=" << &smRMSNoise
		       << "chamberMeanNoise.=" << &chamberMeanNoise
		       << "chamberRMSNoise.=" << &chamberRMSNoise
		       << "runMeanVdrift=" << runMeanVdrift
		       << "runRMSVdrift=" << runRMSVdrift
		       << "smMeanVdrift.=" << &smMeanVdrift
		       << "smRMSVdrift.=" << &smRMSVdrift
		       << "chamberMeanVdrift.=" << &chamberMeanVdrift
		       << "chamberRMSVdrift.=" << &chamberRMSVdrift
		       << "runMeanT0=" << runMeanT0
		       << "runRMST0=" << runRMST0
		       << "smMeanT0.=" << &smMeanT0
		       << "smRMST0.=" << &smRMST0
		       << "chamberMeanT0.=" << &chamberMeanT0
		       << "chamberRMST0.=" << &chamberRMST0
		       << "runBadPadFraction=" << runBadPadFraction
		       << "chamberBadPadFraction.=" << &chamberBadPadFraction
		       << "chamberStatusValues.=" << &chamberStatusValues;
      }
      if(getDCSInfo) {
	(*treeStreamer)<< "trdTree" 
		       << "dcsFeeGlobalNTimeBins=" << dcsFeeGlobalNTimeBins
		       << "dcsFeeGlobalConfigTag=" << dcsFeeGlobalConfigTag
		       << "dcsFeeGlobalSingleHitThres=" << dcsFeeGlobalSingleHitThres
		       << "dcsFeeGlobalThreePadClustThres=" << dcsFeeGlobalThreePadClustThres
		       << "dcsFeeGlobalSelectiveNoSZ=" << dcsFeeGlobalSelectiveNoSZ
		       << "dcsFeeGlobalTCFilterWeight=" << dcsFeeGlobalTCFilterWeight
		       << "dcsFeeGlobalTCFilterShortDecPar=" << dcsFeeGlobalTCFilterShortDecPar
		       << "dcsFeeGlobalTCFilterLongDecPar=" << dcsFeeGlobalTCFilterLongDecPar
		       << "dcsFeeGlobalModeFastStatNoise=" << dcsFeeGlobalModeFastStatNoise
	  //		     << "dcsFeeGlobalConfigVersion.=" << &dcsFeeGlobalConfigVersion
	  //		     << "dcsFeeGlobalConfigName.=" << &dcsFeeGlobalConfigName
	  //		     << "dcsFeeGlobalFilterType.=" << &dcsFeeGlobalFilterType
	  //		     << "dcsFeeGlobalReadoutParam.=" << &dcsFeeGlobalReadoutParam
	  //		     << "dcsFeeGlobalTestPattern.=" << &dcsFeeGlobalTestPattern
	  //		     << "dcsFeeGlobalTrackletMode.=" << &dcsFeeGlobalTrackletMode
	  //		     << "dcsFeeGlobalTrackletDef.=" << &dcsFeeGlobalTrackletDef
	  //		     << "dcsFeeGlobalTriggerSetup.=" << &dcsFeeGlobalTriggerSetup
	  //		     << "dcsFeeGlobalAddOptions.=" << &dcsFeeGlobalAddOptions
		       << "statusDCSFEESOR.=" << &statusArraySOR
		       << "statusDCSFEEEOR.=" << &statusArrayEOR
		       << "SORandEOR=" << sorAndEor
		       << "nChanged=" << nChanged
		       << "nSB1=" << nSB1
		       << "nSB2=" << nSB2
		       << "nSB3=" << nSB3
		       << "nSB4=" << nSB4
		       << "nSB5=" << nSB5;
      }
      (*treeStreamer)<< "trdTree"
		     << "\n";
    }  // end loop over time steps

    // add the run number to the list of runs
    runs.ResizeTo(runs.GetNoElements()+1);
    runs[runs.GetNoElements()-1] = currRun;

    // do some cleaning
    if(grpObject) delete grpObject;
    if(anodeISensors) anodeISensors->Clear();
    if(anodeUSensors) anodeUSensors->Clear();
    if(driftISensors) driftISensors->Clear();
    if(driftUSensors) driftUSensors->Clear();
    if(temperatureSensors) temperatureSensors->Clear();
    if(overpressureSensors) overpressureSensors->Clear();
    if(gasCO2Sensors) gasCO2Sensors->Clear();
    if(gasH2OSensors) gasH2OSensors->Clear();
    if(gasO2Sensors) gasO2Sensors->Clear();
    //if(adcClkPhaseSensors) adcClkPhaseSensors->Clear();
    if(goofieGainSensors) goofieGainSensors->Clear();
    if(goofieHvSensors) goofieHvSensors->Clear();
    if(goofiePressureSensors) goofiePressureSensors->Clear();
    if(goofieTempSensors) goofieTempSensors->Clear();
    if(goofieVelocitySensors) goofieVelocitySensors->Clear();
    if(goofieCO2Sensors) goofieCO2Sensors->Clear();
    if(goofieN2Sensors) goofieN2Sensors->Clear();
    if(chamberGainFactor) delete chamberGainFactor;
    if(padGainFactor) delete padGainFactor;
    if(chamberNoise) delete chamberNoise;
    if(padNoise) delete padNoise;
    if(chamberVdrift) delete chamberVdrift;
    if(padVdrift) delete padVdrift;
    if(chamberT0) delete chamberT0;
    if(padT0) delete padT0;
    if(chamberStatus) delete chamberStatus;
    if(padStatus) delete padStatus;
  } // end while (loop over runs)
  TTimeStamp loopEndTime;

  treeStreamer->GetFile()->cd();
  runs.Write("runs");
  delete treeStreamer;

  // output some job informations
  TTimeStamp jobEndTime;
  cout << "=============================================" << endl;
  cout << "Job launched at          :  " << jobStartTime.AsString() << endl;
  cout << "Loop over runs started at:  " << loopStartTime.AsString() << endl;
  cout << "Loop over runs ended at  :  " << loopEndTime.AsString() << endl;
  cout << "Job ended at             :  " << jobEndTime.AsString() << endl;
  cout << "Initialization duration  :  " 
       << loopStartTime.GetSec() - jobStartTime.GetSec() << " seconds" << endl;
  cout << "Loop over runs duration  :  " 
       << loopEndTime.GetSec() - loopStartTime.GetSec() << " seconds" << endl;
  cout << "Post loop                :  "
       << jobEndTime.GetSec() - loopEndTime.GetSec() << " seconds" << endl;
  cout << "Running time per processed run:  "
       << (loopEndTime.GetSec()-loopStartTime.GetSec())/(runs.GetNoElements()>0 ? Double_t(runs.GetNoElements()) : 1.0)
       << " sec./run" << endl;
  cout << "Running time per input run:  "
       << (loopEndTime.GetSec()-loopStartTime.GetSec())/((rejectedRuns.GetNoElements()+runs.GetNoElements())>0 ? Double_t(runs.GetNoElements()+rejectedRuns.GetNoElements()) : 1.0)
       << " sec./run" << endl;

  // print the runs that had problems
  cout << "number of rejected runs: " << rejectedRuns.GetNoElements() << endl;
  cout << "rejected run numbers" << endl;
  cout << "********************" << endl;
  for(Int_t iRun=0; iRun<rejectedRuns.GetNoElements(); iRun++) {
    cout << rejectedRuns[iRun] << "    ";
    if(iRun%10==0) cout << endl;
  }
  cout << "=============================================" << endl;
  return;
}

//__________________________________________________________________________________________
void ProcessTRDSensorArray(AliTRDSensorArray *sensorArray, TTimeStamp timeStamp, TVectorD &values) {
  // Fill a vector with sensor values for a given time stamp
  // The sensor->Eval() method makes interpolation inside the covered time interval
  // and returns the value at the closest time limit (start or end) outside the covered time range
  AliDCSSensor *sensor;
  for(Int_t i=0; i<sensorArray->NumSensors(); i++) {
    sensor = sensorArray->GetSensorNum(i);
    if(sensor && sensor->GetGraph()) 
      values[i] = sensor->Eval(timeStamp);
    else
      values[i] = -99.;
  }
  return;
}

//__________________________________________________________________________________________
void ProcessTRDCalibArray(AliTRDCalDet* chamberCalib, AliTRDCalPad *padCalib,
			  TString parName,
			  Double_t &runValue, Double_t &runRMS,
			  TVectorD &chamberValues, TVectorD &chamberValuesRMS,
			  TVectorD &superModuleValues, TVectorD &superModuleValuesRMS) {
  // Process the calibrations for a given run.
  // Calculates the run and chamber wise averages
  //

  // check if the calibration parameter is multiplicative or additive
  Bool_t multiplicative = kTRUE;
  if(!parName.CompareTo("T0")) multiplicative = kFALSE;

  // first iteration (calculate all averages and RMS without discrimination on the SM average)
  gRunWiseHisto->Reset();
  for(Int_t iSM = 0; iSM<AliTRDcalibDB::kNsector; iSM++) {   // loop over supermodules
    // reset the super module histogram
    gSuperModuleWiseHisto->Reset();
    // check if SM is installed
    if(!gkSuperModuleStatus[iSM]) continue;
    for(Int_t iChamber=iSM*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber < (iSM+1)*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber++) {  // loop over chambers in this supermodule
      // get the chamber value
      Float_t chamberValue = chamberCalib->GetValue(iChamber);
      // get the ROC object
      AliTRDCalROC *chamberROC = padCalib->GetCalROC(iChamber);
      if(!chamberROC) 
	continue;
      gChamberWiseHisto->Reset();
      for(Int_t iChannel = 0; iChannel < chamberROC->GetNchannels(); iChannel++){ // loop over channels
	// calculate the calibration parameter for this pad
	Float_t padValue;
	if(multiplicative)
	  padValue = chamberValue * chamberROC->GetValue(iChannel);
	else
	  padValue = chamberValue + chamberROC->GetValue(iChannel);
	// fill the run, SM and chamber wise histograms
	gChamberWiseHisto->Fill(padValue);
	// if the parameter is Noise then check if the pad value is not a default one
	// Default value is now 1.2!!!! Check with Raphaelle for more informations
	if(parName.Contains("Noise") &&
	   TMath::Abs(padValue-1.2)<0.00001) continue;
	gSuperModuleWiseHisto->Fill(padValue);
	gRunWiseHisto->Fill(padValue);
      }  // end loop over channels
      // get the chamber wise mean and RMS
      chamberValues[iChamber] = gChamberWiseHisto->GetMean();
      chamberValuesRMS[iChamber] = gChamberWiseHisto->GetRMS();
    }  // end loop over chambers
    // SM wise mean and RMS
    superModuleValues[iSM] = gSuperModuleWiseHisto->GetMean();
    superModuleValuesRMS[iSM] = gSuperModuleWiseHisto->GetRMS();
  }  // end loop over supermodules
  // run wise mean and RMS
  runValue = gRunWiseHisto->GetMean();
  runRMS = gRunWiseHisto->GetRMS();

  // Noise and Gain are finished processing
  if(parName.Contains("Noise") || parName.Contains("Gain"))
    return;
  // second iteration (calculate SM and run wise averages and RMS for Vdrift and T0)
  // The pads with calib parameter equal to the SM average are discarded (default value)
  gRunWiseHisto->Reset();
  for(Int_t iSM = 0; iSM<AliTRDcalibDB::kNsector; iSM++) {   // loop over supermodules
    gSuperModuleWiseHisto->Reset();
    // eliminate the uninstalled super modules
    if(!gkSuperModuleStatus[iSM]) continue;
    for(Int_t iChamber=iSM*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber < (iSM+1)*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber++) {  // loop over chambers
      // the chamber value
      Float_t chamberValue = chamberCalib->GetValue(iChamber);
      AliTRDCalROC *chamberROC = padCalib->GetCalROC(iChamber);
      if(!chamberROC) 
	continue;
      
      for(Int_t iChannel = 0; iChannel < chamberROC->GetNchannels(); iChannel++){ // loop over channels in a chamber
	// get the pad value
	Float_t padValue;
	if(multiplicative)
	  padValue = chamberValue * chamberROC->GetValue(iChannel);
	else
	  padValue = chamberValue + chamberROC->GetValue(iChannel);
	// eliminate from the average and RMS calculation all pads which
	// have the calib parameter equal with the SM average
	if((parName.Contains("Vdrift") || parName.Contains("T0")) && 
	   TMath::Abs(padValue-superModuleValues[iSM])<0.00001) continue;
	gSuperModuleWiseHisto->Fill(padValue);
	gRunWiseHisto->Fill(padValue);
      }   // end loop over channels
    }   // end loop over chambers 
    superModuleValues[iSM] = gSuperModuleWiseHisto->GetMean();
    superModuleValuesRMS[iSM] = gSuperModuleWiseHisto->GetRMS();
  }   // end loop over super modules
  runValue = gRunWiseHisto->GetMean();
  runRMS = gRunWiseHisto->GetRMS();
  return;
}

//__________________________________________________________________________________________
void ProcessTRDstatus(AliTRDCalChamberStatus* chamberStatus, AliTRDCalPadStatus* padStatus,
		      Float_t &runBadPadFraction, TVectorD &chamberBadPadFraction,
		      TVectorD &chamberStatusValues) {
  // Process the pad status. Calculates the fraction of pads with non 0 status
  // run and chamber wise
  //
  Int_t runPadStatusNot0 = 0;
  Int_t runPadStatusAll = 0;

  Int_t superModuleStatus[18] = {1, 1, 0, 0, 0, 0, 0, 1, 1,
				 1, 1, 0, 0, 0, 0, 0, 0, 1};  

  // loop over chambers
  for(Int_t iChamber=0; iChamber < AliTRDcalibDB::kNdet; iChamber++) {
    // check if the chamber is in an installed sector;
    Int_t sm = AliTRDgeometry::GetSector(iChamber);
    if(!superModuleStatus[sm]) continue;

    chamberStatusValues[iChamber] = chamberStatus->GetStatus(iChamber);
    AliTRDCalSingleChamberStatus *singleChamberStatus = padStatus->GetCalROC(iChamber);
    if(!singleChamberStatus)
      continue;
    Int_t chamberPadStatusNot0 = 0;
    Int_t chamberPadStatusAll = 0;
    // loop over channels in a chamber
    for(Int_t iChannel = 0; iChannel < singleChamberStatus->GetNchannels(); iChannel++) {
      if(singleChamberStatus->GetStatus(iChannel) > 0) {
	chamberPadStatusNot0++;
	runPadStatusNot0++;
      }
      chamberPadStatusAll++;
      runPadStatusAll++;
    }
    chamberBadPadFraction[iChamber] = (chamberPadStatusAll>0 ? 
				       Float_t(chamberPadStatusNot0)/Float_t(chamberPadStatusAll) : -99.);
  }
  runBadPadFraction = (runPadStatusAll>0 ? Float_t(runPadStatusNot0)/Float_t(runPadStatusAll) : -99.);
  return;
}

//__________________________________________________________________________________________
void ProcessTRDCalDCSFEE(AliTRDCalDCS *caldcsSOR, AliTRDCalDCS *caldcsEOR, 
			 Int_t &nsb1, Int_t &nsb2, Int_t &nsb3, Int_t &nsb4, Int_t &nsb5,
			 Int_t &nChanged, Bool_t &sorAndEor, 
			 TVectorD &statusArraySOR, TVectorD &statusArrayEOR) {
  //
  // Process the DCS information
  //
  sorAndEor = kTRUE;
  if(!caldcsSOR && !caldcsEOR) {
    sorAndEor = kFALSE;
    return;
  }
  else if(caldcsSOR && !caldcsEOR) {
    sorAndEor = kFALSE;
  }
  else if(!caldcsSOR && caldcsEOR) {
    caldcsSOR = caldcsEOR;
    sorAndEor = kFALSE;
  }

  nsb1 = 0; nsb2 = 0; nsb3 = 0; nsb4 = 0; nsb5 = 0; nChanged = 0;
  for(Int_t iROC=0; iROC<AliTRDcalibDB::kNdet && iROC<caldcsSOR->GetFEEArr()->GetSize(); iROC++) {
    AliTRDCalDCSFEE *dcsSorFee = caldcsSOR->GetCalDCSFEEObj(iROC);
    AliTRDCalDCSFEE *dcsEorFee = caldcsEOR->GetCalDCSFEEObj(iROC);
    if(dcsSorFee) {
      statusArraySOR[iROC] = dcsSorFee->GetStatusBit();
      if(statusArraySOR[iROC] == 1) nsb1++;
      if(statusArraySOR[iROC] == 2) nsb2++;
      if(statusArraySOR[iROC] == 3) nsb3++;
      if(statusArraySOR[iROC] == 4) nsb4++;
      if(statusArraySOR[iROC] == 5) nsb5++;
    }
    if(dcsEorFee) {
      statusArrayEOR[iROC] = dcsEorFee->GetStatusBit();
    }
    if(sorAndEor) {
      if((statusArraySOR[iROC]-statusArrayEOR[iROC]) != 0) nChanged++;
    } 
  }
  return;
}

//__________________________________________________________________________________________
void MakeRunListFromOCDB(const Char_t* directory, const Char_t* outfile, Bool_t fromAlien) {
  //
  // For a given OCDB path dump the list of available run numbers
  //
  if(fromAlien)
    gSystem->Exec(Form("alien_ls %s > temp.txt", directory));
  else
    gSystem->Exec(Form("ls %s > temp.txt", directory));

  ifstream inBuffer("temp.txt");
  if(!inBuffer.is_open()) {
    cout << "File temp.txt not opened! Exiting" << endl;
    return;
  }
  ofstream outBuffer(outfile);
  if(!outBuffer.is_open()) {
    cout << "File runList.txt cannot be opened! Exiting" << endl;
    return;
  }

  while(!inBuffer.eof()) {
    char inputLine[128];
    inBuffer.getline(inputLine, 128, '\n');
    int runLow, runHigh;
    const char* tempLine = inputLine;
    sscanf(tempLine, "Run%d_%d_v1_s0.root", &runLow, &runHigh);
    outBuffer << runLow << endl;
  }

  inBuffer.close();
  outBuffer.close();
  gSystem->Exec("rm temp.txt");
  return;
}
