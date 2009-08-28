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

//-------------------------------------------------------------------------
//                          Class AliGRPPreprocessor
//                  Global Run Parameters (GRP) preprocessor
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//    Modified: Ernesto.Lopez.Torres@cern.ch  CEADEN-CERN
//    Modified: Chiara.Zampolli@cern.ch  CERN
//-------------------------------------------------------------------------

#include <TChain.h>
#include <TList.h>
#include <TMap.h>
#include <TObjString.h>
#include <TGraph.h>

#include <float.h>

#include "AliGRPPreprocessor.h"
#include "AliGRPObject.h"
#include "AliDCSSensor.h"
#include "AliSplineFit.h"
#include "AliDCSSensorArray.h"
//#include "AliRawEventHeaderVersions.h"

#include "AliTriggerConfiguration.h"
#include "AliTriggerRunScalers.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

class AliDCSValue;
class AliShuttleInterface;

// needed for ReceivePromptRecoParameters

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <AliCDBManager.h>
#include <AliCDBMetaData.h>
#include <AliCDBId.h>
#include <AliTriggerConfiguration.h>

const Double_t kFitFraction = 0.7;                 // Fraction of DCS sensor fits required

ClassImp(AliGRPPreprocessor)

//_______________________________________________________________

  const Int_t AliGRPPreprocessor::fgknDAQLbPar = 8; // num parameters in the logbook for PHYSICS runs, when beamType from DAQ logbook != Cosmics
  const Int_t AliGRPPreprocessor::fgknDAQLbParReduced = 7; // num parameters in the logbook for the other cases
  const Int_t AliGRPPreprocessor::fgknDCSDP = 50;   // number of dcs dps
  const Int_t AliGRPPreprocessor::fgknDCSDPHallProbes = 40;   // number of dcs dps
  const char* AliGRPPreprocessor::fgkDCSDataPoints[AliGRPPreprocessor::fgknDCSDP] = {
                   "LHCState",              // missing in DCS
                   "L3Polarity",
                   "DipolePolarity",
                   "LHCLuminosity",         // missing in DCS
                   "BeamIntensity",         // missing in DCS
                   "L3Current",
                   "DipoleCurrent",
		   "L3_BSF17_H1",
		   "L3_BSF17_H2",
		   "L3_BSF17_H3",
		   "L3_BSF17_Temperature",
		   "L3_BSF4_H1",
		   "L3_BSF4_H2",
		   "L3_BSF4_H3",
		   "L3_BSF4_Temperature",
		   "L3_BKF17_H1",
		   "L3_BKF17_H2",
		   "L3_BKF17_H3",
		   "L3_BKF17_Temperature",
		   "L3_BKF4_H1",
		   "L3_BKF4_H2",
		   "L3_BKF4_H3",
		   "L3_BKF4_Temperature",
		   "L3_BSF13_H1",
		   "L3_BSF13_H2",
		   "L3_BSF13_H3",
		   "L3_BSF13_Temperature",
		   "L3_BSF8_H1",
		   "L3_BSF8_H2",
		   "L3_BSF8_H3",
		   "L3_BSF8_Temperature",
		   "L3_BKF13_H1",
		   "L3_BKF13_H2",
		   "L3_BKF13_H3",
		   "L3_BKF13_Temperature",
		   "L3_BKF8_H1",
		   "L3_BKF8_H2",
		   "L3_BKF8_H3",
		   "L3_BKF8_Temperature",
		   "Dipole_Inside_H1",
		   "Dipole_Inside_H2",
		   "Dipole_Inside_H3",
		   "Dipole_Inside_Temperature",
		   "Dipole_Outside_H1",
		   "Dipole_Outside_H2",
		   "Dipole_Outside_H3",
		   "Dipole_Outside_Temperature",
                   "CavernTemperature",
                   "CavernAtmosPressure",
                   "SurfaceAtmosPressure"
                 };

  const char* AliGRPPreprocessor::fgkDCSDataPointsHallProbes[AliGRPPreprocessor::fgknDCSDPHallProbes] = {
		   "L3_BSF17_H1",
		   "L3_BSF17_H2",
		   "L3_BSF17_H3",
		   "L3_BSF17_Temperature",
		   "L3_BSF4_H1",
		   "L3_BSF4_H2",
		   "L3_BSF4_H3",
		   "L3_BSF4_Temperature",
		   "L3_BKF17_H1",
		   "L3_BKF17_H2",
		   "L3_BKF17_H3",
		   "L3_BKF17_Temperature",
		   "L3_BKF4_H1",
		   "L3_BKF4_H2",
		   "L3_BKF4_H3",
		   "L3_BKF4_Temperature",
		   "L3_BSF13_H1",
		   "L3_BSF13_H2",
		   "L3_BSF13_H3",
		   "L3_BSF13_Temperature",
		   "L3_BSF8_H1",
		   "L3_BSF8_H2",
		   "L3_BSF8_H3",
		   "L3_BSF8_Temperature",
		   "L3_BKF13_H1",
		   "L3_BKF13_H2",
		   "L3_BKF13_H3",
		   "L3_BKF13_Temperature",
		   "L3_BKF8_H1",
		   "L3_BKF8_H2",
		   "L3_BKF8_H3",
		   "L3_BKF8_Temperature",
		   "Dipole_Inside_H1",
		   "Dipole_Inside_H2",
		   "Dipole_Inside_H3",
		   "Dipole_Inside_Temperature",
		   "Dipole_Outside_H1",
		   "Dipole_Outside_H2",
		   "Dipole_Outside_H3",
		   "Dipole_Outside_Temperature",
                 };
                 
  const Short_t kSensors = 48; // start index position of sensor in DCS DPs
  const Short_t kNumSensors = 2; // Number of sensors in DCS DPs

  const char* AliGRPPreprocessor::fgkLHCState[20] = {
                   "P", "PREPARE",
                   "J", "PREINJECTION",
                   "I", "INJECTION",
                   "F", "FILLING",
                   "A", "ADJUST",
                   "U", "UNSTABLE BEAMS",
                   "S", "STABLE BEAMS",
                   "D", "BEAM DUMP",
                   "R", "RECOVER",
                   "C", "PRECYCLE"
                 };

  const char* kppError[] = {
                   "",
                   "(DAQ logbook ERROR)",
                   "(DAQ FXS ERROR)",
                   "(Trigger Scalers not found in DCS FXS - ERROR)",
                   "(DCS data points ERROR)",
                   "(Trigger Configuration ERROR)",
                   "(DAQ logbook ERROR determining partition of the run)"
  };

//_______________________________________________________________

AliGRPPreprocessor::AliGRPPreprocessor(AliShuttleInterface* shuttle):
	AliPreprocessor("GRP",shuttle),  fPressure(0), fmaxFloat(0), fminFloat(0),fmaxDouble(0), fminDouble(0), fmaxInt(0), fminInt(0), fmaxUInt(0), fminUInt(0)
{
	// constructor - shuttle must be instantiated!

	AddRunType("COSMIC");
	AddRunType("LASER");
	AddRunType("PHYSICS");
	AddRunType("CALIBRATION_BC");
	AddRunType("CALIBRATION_CENTRAL");
	AddRunType("CALIBRATION_EMD");
	AddRunType("CALIBRATION_MB");
	AddRunType("CALIBRATION_SEMICENTRAL");
	AddRunType("CALIBRATION");
	AddRunType("PEDESTAL");
	AddRunType("STANDALONE");
	AddRunType("GAIN");
	AddRunType("NOISE");
	AddRunType("PULSER");
        AddRunType("STANDALONE_PULSER");

	fmaxFloat = FLT_MAX;
	fminFloat = -FLT_MAX;
	fmaxDouble = DBL_MAX;
	fminDouble = -DBL_MAX;
	fmaxInt = kMaxInt;
	fminInt = kMinInt;
	fmaxUInt = kMaxUInt;
	fminUInt = 0;

	AliInfo(Form("Max allowed float = %6.5e",fmaxFloat));
	AliInfo(Form("Min allowed float = %6.5e",fminFloat));
	AliInfo(Form("Max allowed double = %6.5e",fmaxDouble));
	AliInfo(Form("Min allowed double = %6.5e",fminDouble));
	AliInfo(Form("Max allowed integer = %d",fmaxInt));
	AliInfo(Form("Min allowed integer = %d",fminInt));
	AliInfo(Form("Max allowed unsigned integer = %u",(Int_t)fmaxUInt));
	AliInfo(Form("Min allowed unsigned integer = %u",(Int_t)fminUInt));

}

//_______________________________________________________________

AliGRPPreprocessor::~AliGRPPreprocessor()
{
	//destructor
	
	delete fPressure;
}

//_______________________________________________________________

void AliGRPPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  // Initialize preprocessor

  AliPreprocessor::Initialize(run, startTime, endTime);

  AliInfo("Initialization of the GRP preprocessor.");
  AliInfo(Form("Start Time DCS = %d",GetStartTimeDCSQuery()));
  AliInfo(Form("End Time DCS = %d",GetEndTimeDCSQuery()));
  TClonesArray * array = new TClonesArray("AliDCSSensor",kNumSensors); 
  for(Int_t j = 0; j < kNumSensors; j++) {
    AliDCSSensor * sens = new ((*array)[j])AliDCSSensor;
    sens->SetStringID(fgkDCSDataPoints[j+kSensors]);
  }
  AliInfo(Form("Pressure Entries: %d",array->GetEntries()));

  fPressure = new AliDCSSensorArray(GetStartTimeDCSQuery(), GetEndTimeDCSQuery(), array);
}

//_______________________________________________________________

UInt_t AliGRPPreprocessor::Process(TMap* valueMap)
{
	// process data retrieved by the Shuttle
	
	// retrieving "partition" and "detector" fields from DAQ logbook to 
	// determine the partition in which the run was taken
	// the partition is used to decide how to react in case of errors for CTP

	TString partition = (TString)GetRunParameter("partition");  
	TString detector = (TString)GetRunParameter("detector");   

	AliGRPObject *grpobj = new AliGRPObject();  // object to store data

	//=================//
	// DAQ logbook     //
	//=================//
	UInt_t error = 0;
	
	Int_t iDaqLB = ProcessDaqLB(grpobj);
	TString runType = (TString)GetRunType();
	TString beamType = (TString)GetRunParameter("beamType");
	if((runType == "PHYSICS" && iDaqLB == fgknDAQLbPar && beamType!="Cosmics") ||  (runType == "PHYSICS" && iDaqLB == fgknDAQLbParReduced && beamType=="Cosmics") || (runType != "PHYSICS" && iDaqLB == fgknDAQLbParReduced)) {
		Log(Form("DAQ Logbook, successful!"));
	} else {
		Log(Form("DAQ Logbook, could not get all expected entries!!!"));
		error |= 1;
	}

	//=================//
	// DAQ FXS         //
	//=================//
	UInt_t iDaqFxs = ProcessDaqFxs();
	if( iDaqFxs == 0 ) {
		Log(Form("DAQ FXS, successful!"));
	} else {
		Log(Form("DAQ FXS, could not store run raw tag file!!!"));
		error |= 2;
	}
	
	//=================//
	// DCS FXS         //
	//=================//
	UInt_t iDcsFxs = ProcessDcsFxs(partition, detector);
	if( iDcsFxs == 0 ) {
		Log(Form("DCS FXS, successful!"));
	} else  if (iDcsFxs ==1) {
		Log(Form("DCS FXS, Could not store CTP scalers!!!"));
		error |= 4;
	} else{
		Log(Form("Incorrect field in DAQ logbook for partition = %s and detector = %s, going into error without CTP scalers...",partition.Data(),detector.Data()));
		error |= 32;
	}
	
	//=================//
	// DCS data points //
	//=================//
	Log(Form("Starting DCS Query at %d and finishing at %d",GetStartTimeDCSQuery(),GetEndTimeDCSQuery()));
	Int_t entries = ProcessDcsDPs( valueMap, grpobj );
	Log(Form("entries found = %d (should be %d)",entries, fgknDCSDP-4));
	if( entries < fgknDCSDP-4 ) { // FIXME (!= ) LHState, LHCLuminosity, BeamIntensity, LH3_BSF4_H3 are not working yet...  
		Log(Form("Problem with the DCS data points!!! Only %d/%d entries found",entries,fgknDCSDP-4));
		error |= 8;
	} else  Log(Form("DCS data points, successful!"));
	
	//=======================//
	// Trigger Configuration //
	//=======================//

	const char * triggerConf = GetTriggerConfiguration();

	if (partition.IsNull() && !detector.IsNull()){ // standalone partition
		Log("STANDALONE partition for current run, using Trigger Configuration dummy value");
		AliCDBEntry *cdbEntry = GetFromOCDB("CTP","DummyConfig");
		if (!cdbEntry) {
			Log(Form("No dummy CTP configuration entry found, going into error..."));
			error |= 16;
		}
		else{
			AliTriggerConfiguration *runcfg = (AliTriggerConfiguration*)cdbEntry->GetObject();
			if (!runcfg){
				Log(Form("dummy CTP config not found in OCDB entry, going into error..."));
				error |= 16;
			}
			else {
				TString titleCTPcfg = Form("CTP cfg for run %i from Dummy entry in OCDB",fRun);
				runcfg->SetTitle(titleCTPcfg);
				AliCDBMetaData metaData;
				metaData.SetResponsible("Roman Lietava");
				metaData.SetComment("CTP run configuration from dummy entry in OCDB");
				if (!Store("CTP","Config", runcfg, &metaData, 0, 0)) {
					Log("Unable to store the dummy CTP run configuration object to OCDB!");
					error |= 16;
				}
			}
		}
	}

	else if (!partition.IsNull() && detector.IsNull()){ // global partition
		Log("GLOBAL partition for current run, using Trigger Configuration from DAQ Logbook");
		if (triggerConf!= NULL) {
			Log("Found trigger configuration in DAQ logbook");
			AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfigurationFromString(triggerConf);	  
			if (!runcfg) {
				Log("Bad CTP run configuration file from DAQ logbook! The corresponding CDB entry will not be filled!");
				error |= 16;
			}
			else {
				TString titleCTPcfg = Form("CTP cfg for run %i from DAQ",fRun);
				runcfg->SetTitle(titleCTPcfg);
				AliCDBMetaData metaData;
				metaData.SetBeamPeriod(0);
				metaData.SetResponsible("Roman Lietava");
				metaData.SetComment("CTP run configuration from DAQ logbook");
				if (!Store("CTP","Config", runcfg, &metaData, 0, 0)) {
					Log("Unable to store the CTP run configuration object to OCDB!");
					error |= 16;
				}
			}
		}

		else {
			Log("Trigger configuration NULL in DAQ logbook");
			error |= 16;
		}
	}

	else {
		Log(Form("Incorrect field in DAQ logbook for partition = %s and detector = %s, going into error without trigger configuration...",partition.Data(),detector.Data()));
		error |= 32;
	}

	// storing AliGRPObject in OCDB

	AliCDBMetaData md;
	md.SetResponsible("Chiara Zampolli");
	md.SetComment("Output parameters from the GRP preprocessor.");
	
	Bool_t result = kTRUE;
	result = Store("GRP", "Data", grpobj, &md); 
	delete grpobj;
	
	if (result && !error ) {
		Log("GRP Preprocessor Success");
		return 0;
	} else {
		Log( Form("GRP Preprocessor FAILS!!! %s%s%s%s%s%s",
			  kppError[(error&1)?1:0],
			  kppError[(error&2)?2:0],
			  kppError[(error&4)?3:0],
			  kppError[(error&8)?4:0],
			  kppError[(error&16)?5:0],
			  kppError[(error&32)?6:0]
			  ));
		return error;
	}
}

//_______________________________________________________________

Int_t AliGRPPreprocessor::ProcessDaqLB(AliGRPObject* grpObj)
{
	//Getting the DAQ lb information
	
	time_t timeStart = (time_t)(((TString)GetRunParameter("DAQ_time_start")).Atoi());
	time_t timeEnd = (time_t)(((TString)GetRunParameter("DAQ_time_end")).Atoi());
	Float_t beamEnergy = (Float_t)(((TString)GetRunParameter("beamEnergy")).Atof());
	TString beamType = (TString)GetRunParameter("beamType");
	Char_t numberOfDetectors = (Char_t)(((TString)GetRunParameter("numberOfDetectors")).Atoi());
	UInt_t  detectorMask = (UInt_t)(((TString)GetRunParameter("detectorMask")).Atoi());
	TString lhcPeriod = (TString)GetRunParameter("LHCperiod");
	TString runType = (TString)GetRunType();

	UInt_t nparameter = 0;
	if (timeStart != 0){
		grpObj->SetTimeStart(timeStart);
		Log(Form("Start time for run %d: %d",fRun, (Int_t)timeStart));
		nparameter++;
	} 
	else {
		Log(Form("Start time not put in logbook, setting to invalid in GRP entry!"));
	}

	if (timeEnd != 0){
		grpObj->SetTimeEnd(timeEnd);
		Log(Form("End time for run %d: %i",fRun, (Int_t)timeEnd));
		nparameter++;
	} 
	else {
		Log(Form("End time not put in logbook, setting to invalid in GRP entry!"));
	}

	if (beamEnergy != 0){
		grpObj->SetBeamEnergy(beamEnergy);
		Log(Form("Beam Energy for run %d: %f",fRun, beamEnergy));
		if ((runType == "PHYSICS" && beamType!="Cosmics")){
			nparameter++; // increasing nparameters only in case we're in PHYSICS runs with beamType != Cosmics
		}
	} 
	else {
		if ((runType == "PHYSICS" && beamType!="Cosmics")){
			Log(Form("Beam Energy not put in logbook, setting to invalid in GRP entry, and producing an error (beamType = %s, runType = %s)",beamType.Data(), runType.Data()));
		}
		else{
			Log(Form("Beam Energy not put in logbook, setting to invalid in GRP entry, but not producing any error (beamType = %s, runType = %s)",beamType.Data(), runType.Data()));
		}
	}

		
	if (beamType.Length() != 0){
		grpObj->SetBeamType(beamType);
		Log(Form("Beam Type for run %d: %s",fRun, beamType.Data()));
		nparameter++; 
	} 
	else {
		Log(Form("Beam Type not put in logbook, setting to invalid in GRP entry!"));
	}
		
	if (numberOfDetectors != 0){
		grpObj->SetNumberOfDetectors(numberOfDetectors);
		Log(Form("Number Of Detectors for run %d: %d",fRun, (Int_t)numberOfDetectors));
		nparameter++;
	} 
	else {
		Log(Form("Number Of Detectors not put in logbook, setting to invalid in GRP entry!"));
	}

	if (detectorMask != 0){
		grpObj->SetDetectorMask(detectorMask);
		Log(Form("Detector Mask for run %d: %d",fRun, detectorMask));
		nparameter++;
	} 
	else {
		Log(Form("Detector Mask not put in logbook, setting to invalid in GRP entry!"));
	}

	if (lhcPeriod.Length() != 0) {
		grpObj->SetLHCPeriod(lhcPeriod);
		Log(Form("LHC period (DAQ) for run %d: %s",fRun, lhcPeriod.Data()));
		nparameter++;
	} 
	else {
		Log(Form("LHCperiod not put in logbook, setting to invalid in GRP entry!"));
	}
	if (runType.Length() != 0) {
		grpObj->SetRunType(runType);
		Log(Form("Run Type (DAQ) for run %d: %s",fRun, runType.Data()));
		nparameter++;
	} 
	else {
		Log(Form("Run Type not put in logbook, setting to invalid in GRP entry!"));
	}

	return nparameter;
}

//_______________________________________________________________

UInt_t AliGRPPreprocessor::ProcessDaqFxs()
{
	//======DAQ FXS======//
	
	//	AliRawEventHeaderV3_9::Class()->IgnoreTObjectStreamer(); // to avoid trying reading TObject store in AliRawEventHeaderV3_9 - temporary fix 
	TList* list = GetFileSources(kDAQ);  
	if (!list) {
		Log("No raw data tag list: connection problems with DAQ FXS logbook!");
		return 1;
	}
	
	if (list->GetEntries() == 0) {
		Log("no raw data tags in this run: nothing to merge!");
		delete  list; list=0;
		return 0;
	}
	
	TChain *fRawTagChain = new TChain("T");
	Int_t nFiles=0;
	TIterator* iter = list->MakeIterator();
	TObject* obj = 0;
	while ((obj = iter->Next())) {
		TObjString* objStr = dynamic_cast<TObjString*> (obj);
		if (objStr) {
			Log(Form("Found source %s", objStr->String().Data()));
			TList* list2 = GetFileIDs(kDAQ, objStr->String());
			if (!list2) {
				Log("No list with ids from DAQ was found: connection problems with DAQ FXS logbook!");
				delete fRawTagChain; fRawTagChain=0;
				return 1;
			}
			Log(Form("Number of ids: %d",list2->GetEntries()));
			for(Int_t i = 0; i < list2->GetEntries(); i++) {
				TObjString *idStr = (TObjString *)list2->At(i);
				TString fileName = GetFile(kDAQ,idStr->String().Data(),objStr->String().Data());
				if (fileName.Length() > 0) {
					Log(Form("Adding file in the chain: %s",fileName.Data()));
					fRawTagChain->Add(fileName.Data());
					nFiles++;
				} else {
					Log(Form("Could not retrieve file with id %s from source %s: "
						 "connection problems with DAQ FXS!",
						 idStr->String().Data(), objStr->String().Data()));
					delete list; list=0;
					delete list2; list2=0;
					delete fRawTagChain; fRawTagChain=0;
					return 2;
				}
			}
			delete list2;
		}
	}
	
	TString fRawDataFileName = "GRP_Merged.tag.root";
	Log(Form("Merging %d raw data tags into file: %s", nFiles, fRawDataFileName.Data()));
	
	if( fRawTagChain->Merge(fRawDataFileName) < 1 ) {
		Log("Error merging raw data files!!!");
		return 3;
	}
	
	TString outputfile = Form("Run%d.Merged.RAW.tag.root", fRun);
	Bool_t result = StoreRunMetadataFile(fRawDataFileName.Data(),outputfile.Data());
	
	if (!result) {
		Log("Problem storing raw data tags in local file!!!");
	} else {
		Log("Raw data tags merged successfully!!");
	}
	
	delete iter;
	delete list;
	delete fRawTagChain; fRawTagChain=0;
	
	if (result == kFALSE) {
		return 4;
	}
	
	return 0;
	
}

//_______________________________________________________________
UInt_t AliGRPPreprocessor::ProcessDcsFxs(TString partition, TString detector)
{

	// processing the info
	// stored in the DCS FXS
	// coming from the trigger

	// Get the CTP counters information

	if (partition.IsNull() && !detector.IsNull()){ // standalone partition
		Log("STANDALONE partition for current run, using Trigger Scalers dummy value");
		AliCDBEntry *cdbEntry = GetFromOCDB("CTP","DummyScalers");
		if (!cdbEntry) {
			Log(Form("No dummy CTP scalers entry found, going into error..."));
			return 1;
		}
		else{
			AliTriggerRunScalers *scalers = (AliTriggerRunScalers*)cdbEntry->GetObject();
			if (!scalers){
				Log(Form("CTP dummy scalers not found in OCDB entry, going into error..."));
				return 1;
			}
			else {
				AliCDBMetaData metaData;
				metaData.SetResponsible("Roman Lietava");
				metaData.SetComment("CTP scalers from dummy entry in OCDB");
				if (!Store("CTP","Scalers", scalers, &metaData, 0, 0)) {
					Log("Unable to store the dummy CTP scalers object to OCDB!");
					return 1;
				}
			}
		}
	}

	else if (!partition.IsNull() && detector.IsNull()){ // global partition
		Log("GLOBAL partition for current run, using CTP scalers from DCS FXS");
		TString countersfile = GetFile(kDCS, "CTP_xcounters","");
		if (countersfile.IsNull()) {
			Log("No CTP counters files has been found: empty source!");
			return 1;
		}
		else {
			Log(Form("File with Id CTP_xcounters found in DCS FXS! Copied to %s",countersfile.Data()));
			AliTriggerRunScalers *scalers = AliTriggerRunScalers::ReadScalers(countersfile);
			if (!scalers) {
				Log("Bad CTP counters file! The corresponding CDB entry will not be filled!");
				return 1;
			}
			else {
				AliCDBMetaData metaData;
				metaData.SetBeamPeriod(0);
				metaData.SetResponsible("Roman Lietava");
				metaData.SetComment("CTP scalers");
				if (!Store("CTP","Scalers", scalers, &metaData, 0, 0)) {
					Log("Unable to store the CTP scalers object to OCDB!");
					return 1;
				}
			}
		}
	}
	

	else{	
		Log(Form("Incorrect field in DAQ logbook for partition = %s and detector = %s, going into error...",partition.Data(),detector.Data()));
		return 2;
	}

	return 0;

}
//_______________________________________________________________

Int_t AliGRPPreprocessor::ProcessDcsDPs(TMap* valueMap, AliGRPObject* grpObj)
{

	//
	// processing DCS DPs
	//

	Int_t entries = 0;  // counting the entries that are in the DCS DB, not taking care whether they have values or not
	Int_t nLHCEntries = 0;
	Int_t nL3Entries = 0;
	Int_t nDipoleEntries = 0;
	Int_t nEnvEntries = 0;
	Int_t nHallProbesEntries = 0;
	nLHCEntries = ProcessLHCDPs(valueMap, grpObj);
	nL3Entries = ProcessL3DPs(valueMap, grpObj);
	nDipoleEntries = ProcessDipoleDPs(valueMap, grpObj);
	nEnvEntries = ProcessEnvDPs(valueMap, grpObj);
	nHallProbesEntries = ProcessHPDPs(valueMap, grpObj);

	entries = nLHCEntries + nL3Entries + nDipoleEntries + nEnvEntries + nHallProbesEntries;
	return entries;

}

//_______________________________________________________________

Int_t AliGRPPreprocessor::ProcessL3DPs(const TMap* valueMap, AliGRPObject* grpObj)
{

	// processing DPs
	// related to 
	// L3 info

	Int_t nL3Entries = 0;
	TObjArray *array = 0x0;
	Int_t indexDP = -1;
	Bool_t isZero = kTRUE; // flag to monitor L3Current. If set to true, the magnet is OFF, and the polarity can change

	AliInfo(Form("==========L3Current==========="));
	Bool_t outOfRange = kFALSE;  // flag to monitor if any value collected by DCS is out of range
	indexDP = kL3Current;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			Float_t *floatDCS = ProcessFloatAllMagnet(array, indexDP, isZero);
			if (floatDCS != NULL){
				grpObj->SetL3Current(floatDCS);
			}
			else{
				outOfRange = kTRUE;
			}	
			if (floatDCS){
				delete[] floatDCS;
				floatDCS = 0x0;
			}
		}
		if (!outOfRange) nL3Entries++;
	}

	if (array) array = 0x0;

	AliInfo(Form("==========L3Polarity==========="));
	indexDP = kL3Polarity;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s Polarity to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			Bool_t change = kFALSE;
			Char_t charDCS = ProcessBool(array,change);
			if (change == kFALSE){
				grpObj->SetL3Polarity(charDCS);
				AliInfo(Form("%s set to %d",fgkDCSDataPoints[indexDP],(Int_t)(grpObj->GetL3Polarity())));
				nL3Entries++;
			}
			else if (isZero){
				AliInfo(Form("%s set to invalid, but magnet was OFF (according to the current), DP not considered wrong",fgkDCSDataPoints[indexDP]));
				nL3Entries++;
			}
			else {
				AliError(Form("%s value changed within the run, while the magnet was ON (according to the current), setting it to invalid and considering the DP as wrong",fgkDCSDataPoints[indexDP]));
			}
		}
	}

	return nL3Entries;

}
//_______________________________________________________________

Int_t AliGRPPreprocessor::ProcessDipoleDPs(const TMap* valueMap, AliGRPObject* grpObj)
{
	// processing DPs
	// related to 
	// the Dipole info

	Int_t nDipoleEntries = 0;
	TObjArray *array = 0x0;
	Int_t indexDP = -1;
	Bool_t isZero = kTRUE; // flag to monitor L3Current. If set to true, the magnet is OFF, and the polarity can change

	AliInfo(Form("==========DipoleCurrent==========="));
	Bool_t outOfRange = kFALSE;  // flag to monitor if any value collected by DCS is out of range
	indexDP = kDipoleCurrent;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			Float_t *floatDCS = ProcessFloatAllMagnet(array, indexDP, isZero);
			if (floatDCS != NULL){
				grpObj->SetDipoleCurrent(floatDCS);
			} 
			else{
				outOfRange=kTRUE;
			}
			if (floatDCS){
				delete[] floatDCS;
				floatDCS = 0x0;
			}
		}
		if (!outOfRange) nDipoleEntries++;
	}

	if (array) array = 0x0;

	AliInfo(Form("==========DipolePolarity==========="));
	indexDP = kDipolePolarity;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			Bool_t change = kFALSE;
			Char_t charDCS = ProcessBool(array,change);
			if (!change){
				grpObj->SetDipolePolarity(charDCS);
				AliInfo(Form("%s set to %d",fgkDCSDataPoints[indexDP],(Int_t)(grpObj->GetDipolePolarity())));
				nDipoleEntries++;
			}
			else if (isZero){
				AliInfo(Form("%s set to invalid, but magnet was OFF (according to the current), DP not considered wrong",fgkDCSDataPoints[indexDP]));
				nDipoleEntries++;
			}
			else{
				AliError(Form("%s value changed within the run while the magnet was ON (according to the current), setting it to invalid and considering the DP as wrong", fgkDCSDataPoints[indexDP]));
			}
		}
	}

	return nDipoleEntries;

}
//_______________________________________________________________

Int_t AliGRPPreprocessor::ProcessEnvDPs(TMap* valueMap, AliGRPObject* grpObj)
{
	// processing DPs
	// related to 
	// evironment conditions (temperature, pressure) info

	Int_t nEnvEntries = 0;
	TObjArray *array = 0x0;
	Int_t indexDP = -1;

	AliInfo(Form("==========CavernTemperature==========="));
	Bool_t outOfRange = kFALSE;  // flag to monitor if any value collected by DCS is out of range
	indexDP = kCavernTemperature;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			Float_t *floatDCS = ProcessFloatAll(array);
			if (floatDCS != NULL){
				grpObj->SetCavernTemperature(floatDCS);
			}
			else{
				outOfRange = kTRUE;
			}
			if (floatDCS){
				delete[] floatDCS;
				floatDCS = 0x0;
			}
		}
		if (!outOfRange) nEnvEntries++;
	}

	if (array) array = 0x0;

	AliInfo(Form("==========AtmosPressures (Cavern + Surface)==========="));
	AliDCSSensorArray *dcsSensorArray = GetPressureMap(valueMap);
	dcsSensorArray->Print();
	if( fPressure->NumFits()==0 ) {
		Log("Problem with the pressure sensor values!!!");
	} 
	else {
		AliInfo(Form("==========CavernAtmosPressure==========="));
		indexDP = kCavernAtmosPressure;
		AliDCSSensor* sensorCavernP2 = dcsSensorArray->GetSensor(fgkDCSDataPoints[indexDP]);
		AliDebug(2,Form("sensorCavernP2 = %p", sensorCavernP2));
		if( sensorCavernP2->GetFit() ) {
			Log(Form("<%s> for run %d: Sensor Fit found",fgkDCSDataPoints[indexDP], fRun));
			grpObj->SetCavernAtmosPressure(sensorCavernP2);
			nEnvEntries++;
		} 
		//if (sensorP2) delete sensorP2;
		else {
			Log(Form("ERROR Sensor Fit for %s not found ", fgkDCSDataPoints[indexDP] ));
		}
		AliInfo(Form("==========SurfaceAtmosPressure==========="));
		indexDP = kSurfaceAtmosPressure;
		AliDCSSensor* sensorP2 = dcsSensorArray->GetSensor(fgkDCSDataPoints[indexDP]);
		AliDebug(2,Form("sensorP2 = %p", sensorP2));
		if( sensorP2->GetFit() ) {
			Log(Form("<%s> for run %d: Sensor Fit found",fgkDCSDataPoints[indexDP], fRun));
			grpObj->SetSurfaceAtmosPressure(sensorP2);
			nEnvEntries++;
		} 
		//if (sensorP2) delete sensorP2;
		else {
			Log(Form("ERROR Sensor Fit for %s not found ", fgkDCSDataPoints[indexDP] ));
		}
		
	}

	return nEnvEntries;
}
//_______________________________________________________________

Int_t AliGRPPreprocessor::ProcessHPDPs(const TMap* valueMap, AliGRPObject* grpObj)
{
	// processing DPs
	// related to 
	// Hall Probes info

	Int_t nHPEntries = 0;
	TObjArray *array = 0x0;
	Int_t indexDP = -1;
	Bool_t outOfRange; // flag to monitor if any value collected by DCS is out of range

	if (fgknDCSDPHallProbes != AliGRPObject::GetNumberOfHP()){
		AliError(Form("Number of Hall probes expected in GRP Preprocessor (i.e. %d) different from number of Hall Probes foreseen in GRP object (i.e. %d). Looping on entries from GRP object anyway.", fgknDCSDPHallProbes, AliGRPObject::GetNumberOfHP()));
	}
	for (indexDP = 0; indexDP < AliGRPObject::GetNumberOfHP(); indexDP++){
		outOfRange = kFALSE; // resetting outOfRange flag at each HP
		AliInfo(Form("==========%s===========",AliGRPObject::GetHPDP(indexDP)));
		array = (TObjArray *)valueMap->GetValue(AliGRPObject::GetHPDP(indexDP));
		if(!array) {
			Log(Form("%s not found in the map!!!",AliGRPObject::GetHPDP(indexDP)));
		} 
		else {
			if (array->GetEntries() == 0){
				AliError(Form("No entries found in array! setting %s to invalid...",AliGRPObject::GetHPDP(indexDP)));
			}
			else {
				Float_t *floatDCS = ProcessFloatAll(array);
				if (floatDCS != NULL){
					AliDebug(2,Form("value[0] = %f, value[1] = %f, value[2] = %f, value[3] = %f, value[4] = %f",floatDCS[0],floatDCS[1],floatDCS[2],floatDCS[3],floatDCS[4])); 
					grpObj->SetHallProbes((AliGRPObject::DP_HallProbes)indexDP,floatDCS);
					for (Int_t kk = 0 ; kk< 5; kk++){
						AliDebug(2,Form("HallProbe[%d][%d]=%f",indexDP,kk,grpObj->GetHallProbes((AliGRPObject::DP_HallProbes)indexDP,(AliGRPObject::Stats)kk)));
					}
				}
				else{
					outOfRange = kTRUE;
				}
				if (floatDCS){
					delete[] floatDCS;
					floatDCS = 0x0;
				}
			}
			if (!outOfRange) nHPEntries++;
		}
	}
		
	Log(Form("Hall Probes = %d ", nHPEntries));
	return nHPEntries;
}

//_______________________________________________________________

Int_t AliGRPPreprocessor::ProcessLHCDPs(const TMap* valueMap, AliGRPObject* grpObj)
{

	//
	// processing of LHC related DCS DPs, i.e.:
	// LHCState
	// LHCLuminosity
	// BeamIntensity
	//

	Int_t nLHCEntries = 0;
	TObjArray *array = 0x0;
	Int_t indexDP = -1;

	AliInfo(Form("==========LHCState==========="));
	indexDP = kLHCState;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			TString stringDCS = ProcessChar(array);
			if (stringDCS.Length()!=0) {
				Bool_t found = kFALSE;
				for( Int_t i=0; i<20; i+=2 ) {
					if( stringDCS.CompareTo(fgkLHCState[i]) == 0 ) {
						stringDCS = fgkLHCState[i+1];
						found = kTRUE;
						break;
					}
				}
				if (found){
					Log(Form("<%s> for run %d: %s",fgkDCSDataPoints[indexDP],fRun, stringDCS.Data()));
					grpObj->SetLHCState(stringDCS);
				}
				else{
					Log(Form("%s values found not valid!",fgkDCSDataPoints[indexDP]));
					grpObj->SetLHCState(AliGRPObject::GetInvalidString());
				} 
			}
			else {
				Log(Form("%s not valid (null length), string set as invalid!",fgkDCSDataPoints[indexDP]));
				grpObj->SetLHCState(AliGRPObject::GetInvalidString());
			}	  
		}
		nLHCEntries++;
	}
	
	if (array) array = 0x0;

	AliInfo(Form("==========LHCLuminosity==========="));
	Bool_t outOfRange = kFALSE; // flag to monitor if any value collected by DCS is out of range
	indexDP = kLHCLuminosity;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s and its Spline Fit to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			Float_t *floatDCS = ProcessFloatAll(array);
			if (floatDCS != NULL){
				grpObj->SetLHCLuminosity(floatDCS);
				AliSplineFit* splfit = GetSplineFit(array,fgkDCSDataPoints[indexDP]);
				grpObj->SetLHCLuminositySplineFit(splfit);
			//		delete splfit;
			}
			else {
				outOfRange = kTRUE;
			}
			if (floatDCS){
				delete[] floatDCS;
				floatDCS = 0x0;
			}
		}
		if (!outOfRange) nLHCEntries++;
	}

	if (array) array = 0x0;

	AliInfo(Form("==========BeamIntensity==========="));
	if (outOfRange) outOfRange = kFALSE;  // resetting outOfRange if needed
	indexDP = kBeamIntensity;
	array = (TObjArray *)valueMap->GetValue(fgkDCSDataPoints[indexDP]);
	if(!array) {
		Log(Form("%s not found in the map!!!",fgkDCSDataPoints[indexDP]));
	} 
	else {
		if (array->GetEntries() == 0){
			AliError(Form("No entries found in array! setting %s and its Spline Fit to invalid...",fgkDCSDataPoints[indexDP]));
		}
		else {
			Float_t *floatDCS = ProcessFloatAll(array);
			if (floatDCS != NULL){
				grpObj->SetBeamIntensity(floatDCS);
				AliSplineFit* splfit1 = GetSplineFit(array,fgkDCSDataPoints[indexDP]);
				grpObj->SetBeamIntensitySplineFit(splfit1);
				//delete splfit;
			}
			else{
				outOfRange = kTRUE;
			}
			if (floatDCS){
				delete[] floatDCS;
				floatDCS = 0x0;
			}
		}
		if (!outOfRange) nLHCEntries++;
	}

	return nLHCEntries;
}
//_________________________________________________________________________

AliSplineFit* AliGRPPreprocessor::GetSplineFit(const TObjArray *array, const TString& stringID){


	// 
	// returning Spline Fit 
	// 

	Int_t entriesarray = array->GetEntries();
	Float_t* value = new Float_t[entriesarray];
	Float_t* time = new Float_t[entriesarray];
	AliDCSValue* v = 0x0;
	for (Int_t iarray = 0; iarray < entriesarray; iarray++){
		v = (AliDCSValue*)array->At(iarray);
		value[iarray] = v->GetFloat();
		time[iarray] = v->GetTimeStamp();
		AliDebug(2,Form("iarray = %d, value = %f, time = %f",iarray,value[iarray],time[iarray]));
	}
	TGraph* gr = new TGraph(entriesarray,value,time);
	if (!gr ) {
		AliWarning(Form("%s: no input graph to compute SplineFit",stringID.Data()));
		return NULL;
	}
	AliSplineFit *fit = new AliSplineFit();
	fit->SetMinPoints(10);
	fit->InitKnots(gr,10,10,0.0);
	fit->SplineFit(2);
	fit->Cleanup();
	if (!fit) {
		AliWarning(Form("%s: no fit performed",stringID.Data()));
		return NULL;
	} 
	return fit;
}

//_________________________________________________________________________

TString AliGRPPreprocessor::ProcessChar(const TObjArray *array)
{

	// 
	// processing char
	//

	TString aDCSString="";
	
	AliDCSValue *v = 0x0;
	for(Int_t iCount = 0; iCount < array->GetEntries(); iCount++) {
		v = (AliDCSValue *)array->At(iCount);
		if (((Int_t)(v->GetTimeStamp()) < (Int_t)GetStartTimeDCSQuery()) || ((Int_t)(v->GetTimeStamp()) > (Int_t)GetEndTimeDCSQuery())) {
			AliError(Form("DCS values for the parameter outside the queried interval"));
			continue;
		}
		if (iCount > 0) {
			if (aDCSString != v->GetChar())
			AliError(Form("DCS values for the parameter changed from %s to %c within the queried interval", aDCSString.Data(), (Char_t)v->GetChar()));
		}
		aDCSString = (TString)v->GetChar();  // keeping always last value in the array
	}
	return aDCSString;
}

//__________________________________________________________________________________________________________________

Float_t* AliGRPPreprocessor::ProcessFloatAll(const TObjArray* array)
{
	// 
	// processing Float values using Mean, Median, Standard Deviation wrt Mean, Standar Deviation wrt Median 
	//
	// parameters[0] = mean
	// parameters[1] = truncated mean (calculated excluding points outside +/- 3RMS from mean
	// parameters[2] = median
	// parameters[3] = standard deviation wrt mean
	// parameters[4] = standard deviation wrt median
	//

	Float_t* parameters = new Float_t[5];
	Double_t aDCSArrayMean = 0;     // Mean
	Double_t aDCSArrayTruncMean = 0;// Truncated Mean
	Double_t aDCSArrayMedian = 0;   // Median
	Double_t aDCSArraySDMean = 0;   // Standard Deviation wrt Mean
	Double_t aDCSArraySDMedian = 0; // Standard Deviation wrt Median
	Float_t aDCSArraySum = 0.0;
	Int_t iCounts = 0;
	Int_t iCounts1 = 0;
	Float_t temp = 0;
	Float_t temp1 = 0;
	Int_t nCounts = array->GetEntries();
	Float_t *tempArray = new Float_t[nCounts];
	for(Int_t i = 0; i < nCounts; i++) {
		AliDCSValue *v = (AliDCSValue *)array->At(i);
		if ((v->GetFloat() <= fminFloat) || (v->GetFloat() >= fmaxFloat)) {
			AliError(Form("Error! Float value found in DCS map at %d-th entry is OUT OF RANGE: value = %6.5e",i,v->GetFloat()));
			if (v->GetFloat() < fminFloat) AliInfo(Form("The value is smaller than %6.5e",fminFloat));
			if (v->GetFloat() > fmaxFloat) AliInfo(Form("The value is greater than %6.5e",fmaxFloat));
			return NULL;
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			aDCSArraySum += v->GetFloat();
			tempArray[i] = v->GetFloat();
			AliDebug(2,Form("%d-th entry = %f",i,tempArray[i]));
			iCounts += 1;
		}
		else {
			AliError(Form("DCS values for the parameter outside the queried interval"));
		}
	}

	AliDebug(2,Form("Using %i entries, starting from %i entries",iCounts,nCounts));
    	if(iCounts != 0) {
		aDCSArrayMean = TMath::Mean(iCounts,tempArray);
		aDCSArrayMedian = TMath::Median(iCounts,tempArray);
		aDCSArraySDMean = TMath::RMS(iCounts,tempArray);
		AliDebug(2,Form("SD = %f",aDCSArraySDMean));
		// computing standard deviation wrt median
		AliDebug(2,Form("maximum = %f, minimum = %f", aDCSArrayMean+3*aDCSArraySDMean, aDCSArrayMean-3*aDCSArraySDMean));
		for (Int_t i = 0; i < iCounts; i++){
			AliDCSValue *v = (AliDCSValue *)array->At(i);
			AliDebug(3,Form("maximum = %f, minimum = %f", aDCSArrayMean+3*aDCSArraySDMean, aDCSArrayMean-3*aDCSArraySDMean));
			AliDebug(3,Form("%i-th entry = %f",i, v->GetFloat())); 
			if ((v->GetFloat()<=aDCSArrayMean+3*aDCSArraySDMean) && (v->GetFloat()>=aDCSArrayMean-3*aDCSArraySDMean)){
				temp1+=v->GetFloat();
				iCounts1++;
				AliDebug(3,Form("temp1 = %f, iCounts1 = %i",temp1,iCounts1));
			}
    			temp += (v->GetFloat()-aDCSArrayMedian)*(v->GetFloat()-aDCSArrayMedian);
		}
		AliDebug(3,Form("temp before the ratio = %f, with %d counts", temp, iCounts));
		temp/=iCounts;
		AliDebug(3,Form("temp after the ratio = %f", temp));
    		if (temp>0) {
			aDCSArraySDMedian = TMath::Sqrt(temp);
		}
		else if (temp==0) {
			AliInfo(Form("Radical = 0 in computing standard deviation wrt median! Setting it to zero...."));
			aDCSArraySDMedian = 0;
		}
		else{
			AliError(Form("Radical < 0 in computing standard deviation! Setting it to invalid...."));
			aDCSArraySDMedian = AliGRPObject::GetInvalidFloat();
		}
	}
	else {
		aDCSArrayMean = AliGRPObject::GetInvalidFloat();
		aDCSArrayMedian = AliGRPObject::GetInvalidFloat();
		aDCSArraySDMean = AliGRPObject::GetInvalidFloat();
	}
	AliDebug(3,Form("iCounts1 = %d and temp1 = %f",iCounts1, temp1));
	if (iCounts1 > 0) {
		aDCSArrayTruncMean = temp1/iCounts1;
	}
	else{
		aDCSArrayTruncMean = AliGRPObject::GetInvalidFloat();
	}
	
	
	
	AliDebug(2,Form("mean within %d counts = %f ",iCounts,aDCSArrayMean));
	AliDebug(2,Form("truncated mean within %d counts = %f (%i values used)",iCounts,aDCSArrayTruncMean,iCounts1));
	AliDebug(2,Form("median within %d counts = %f ",iCounts,aDCSArrayMedian));
	AliDebug(2,Form("standard deviation with mean within %d counts = %f ",iCounts,aDCSArraySDMean));
	AliDebug(2,Form("standard deviation with median within %d counts = %f ",iCounts,aDCSArraySDMedian));
	
    	parameters[0] = aDCSArrayMean;
	parameters[1] = aDCSArrayTruncMean;
	parameters[2] = aDCSArrayMedian;
	parameters[3] = aDCSArraySDMean;
	parameters[4] = aDCSArraySDMedian;

	AliDebug(2,Form("mean = %f, truncated mean = %f, median = %f, SD wrt mean = %f, SD wrt median = %f ",parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]));
    	//AliInfo(Form("mean = %f, truncated mean = %f, median = %f, SD wrt mean = %f, SD wrt median = %f ",parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]));

	return parameters;
}



//__________________________________________________________________________________________________________________

Float_t* AliGRPPreprocessor::ProcessFloatAllMagnet(const TObjArray* array, Int_t indexDP, Bool_t &isZero)
{
	// 
	// processing Float values using Mean, Median, Standard Deviation wrt Mean, Standar Deviation wrt Median 
	// used for L3 and Dipole magnets, using isZero flag to decide whther the magnet was OFF/ON
	// the flag is set according to the L3/Dipole current value
	// current threshold for L3 = 350 A (value provided by DCS)
	// current threshold for Dipole = 450 A (value provided by DCS)
	//
	// parameters[0] = mean
	// parameters[1] = truncated mean (calculated excluding points outside +/- 3RMS from mean
	// parameters[2] = median
	// parameters[3] = standard deviation wrt mean
	// parameters[4] = standard deviation wrt median
	//

	AliInfo(Form("indexDP = %d",indexDP)); 
	Float_t* parameters = new Float_t[5];
	Double_t aDCSArrayMean = 0;     // Mean
	Double_t aDCSArrayTruncMean = 0;// Truncated Mean
	Double_t aDCSArrayMedian = 0;   // Median
	Double_t aDCSArraySDMean = 0;   // Standard Deviation wrt Mean
	Double_t aDCSArraySDMedian = 0; // Standard Deviation wrt Median
	Float_t aDCSArraySum = 0.0;
	Int_t iCounts = 0;
	Int_t iCounts1 = 0;
	Float_t temp = 0;
	Float_t temp1 = 0;
	Int_t nCounts = array->GetEntries();
	Float_t *tempArray = new Float_t[nCounts];
	for(Int_t i = 0; i < nCounts; i++) {
		AliDCSValue *v = (AliDCSValue *)array->At(i);
		if ((v->GetFloat() <= fminFloat) || (v->GetFloat() >= fmaxFloat)) {
			AliError(Form("Error! Float value found in DCS map at %d-th entry is OUT OF RANGE: value = %6.5e",i,v->GetFloat()));
			if (v->GetFloat() < fminFloat) AliInfo(Form("The value is smaller than %6.5e",fminFloat));
			if (v->GetFloat() > fmaxFloat) AliInfo(Form("The value is greater than %6.5e",fmaxFloat));
			return NULL;
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			aDCSArraySum += v->GetFloat();
			tempArray[i] = v->GetFloat();
			AliDebug(2,Form("%d-th entry = %f",i,tempArray[i]));
			iCounts += 1;
			if (indexDP == kL3Current && v->GetFloat() > 350 && isZero == kTRUE) isZero=kFALSE; 
			if (indexDP == kDipoleCurrent && v->GetFloat() > 450 && isZero == kTRUE) isZero=kFALSE; 
		}
		else {
			AliError(Form("DCS values for the parameter outside the queried interval"));
		}
	}

	AliDebug(2,Form("Using %i entries, starting from %i entries",iCounts,nCounts));
    	if(iCounts != 0) {
		aDCSArrayMean = TMath::Mean(iCounts,tempArray);
		aDCSArrayMedian = TMath::Median(iCounts,tempArray);
		aDCSArraySDMean = TMath::RMS(iCounts,tempArray);
		AliDebug(2,Form("SD = %f",aDCSArraySDMean));
		// computing standard deviation wrt median
		AliDebug(2,Form("maximum = %f, minimum = %f", aDCSArrayMean+3*aDCSArraySDMean, aDCSArrayMean-3*aDCSArraySDMean));
		for (Int_t i = 0; i < iCounts; i++){
			AliDCSValue *v = (AliDCSValue *)array->At(i);
			AliDebug(3,Form("maximum = %f, minimum = %f", aDCSArrayMean+3*aDCSArraySDMean, aDCSArrayMean-3*aDCSArraySDMean));
			AliDebug(3,Form("%i-th entry = %f",i, v->GetFloat())); 
			if ((v->GetFloat()<=aDCSArrayMean+3*aDCSArraySDMean) && (v->GetFloat()>=aDCSArrayMean-3*aDCSArraySDMean)){
				temp1+=v->GetFloat();
				iCounts1++;
				AliDebug(3,Form("temp1 = %f, iCounts1 = %i",temp1,iCounts1));
			}
    			temp += (v->GetFloat()-aDCSArrayMedian)*(v->GetFloat()-aDCSArrayMedian);
		}
		AliDebug(3,Form("temp before the ratio = %f, with %d counts", temp, iCounts));
		temp/=iCounts;
		AliDebug(3,Form("temp after the ratio = %f", temp));
    		if (temp>0) {
			aDCSArraySDMedian = TMath::Sqrt(temp);
		}
		else if (temp==0) {
			AliInfo(Form("Radical = 0 in computing standard deviation wrt median! Setting it to zero...."));
			aDCSArraySDMedian = 0;
		}
		else{
			AliError(Form("Radical < 0 in computing standard deviation! Setting it to invalid...."));
			aDCSArraySDMedian = AliGRPObject::GetInvalidFloat();
		}
	}
	else {
		aDCSArrayMean = AliGRPObject::GetInvalidFloat();
		aDCSArrayMedian = AliGRPObject::GetInvalidFloat();
		aDCSArraySDMean = AliGRPObject::GetInvalidFloat();
	}
	AliDebug(3,Form("iCounts1 = %d and temp1 = %f",iCounts1, temp1));
	if (iCounts1 > 0) {
		aDCSArrayTruncMean = temp1/iCounts1;
	}
	else{
		aDCSArrayTruncMean = AliGRPObject::GetInvalidFloat();
	}
	
	
	
	AliDebug(2,Form("mean within %d counts = %f ",iCounts,aDCSArrayMean));
	AliDebug(2,Form("truncated mean within %d counts = %f (%i values used)",iCounts,aDCSArrayTruncMean,iCounts1));
	AliDebug(2,Form("median within %d counts = %f ",iCounts,aDCSArrayMedian));
	AliDebug(2,Form("standard deviation with mean within %d counts = %f ",iCounts,aDCSArraySDMean));
	AliDebug(2,Form("standard deviation with median within %d counts = %f ",iCounts,aDCSArraySDMedian));
	
    	parameters[0] = aDCSArrayMean;
	parameters[1] = aDCSArrayTruncMean;
	parameters[2] = aDCSArrayMedian;
	parameters[3] = aDCSArraySDMean;
	parameters[4] = aDCSArraySDMedian;

	AliDebug(2,Form("mean = %f, truncated mean = %f, median = %f, SD wrt mean = %f, SD wrt median = %f ",parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]));
    	//AliInfo(Form("mean = %f, truncated mean = %f, median = %f, SD wrt mean = %f, SD wrt median = %f ",parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]));

	return parameters;
}


//_______________________________________________________________

Char_t AliGRPPreprocessor::ProcessBool(const TObjArray* array, Bool_t &change)
{
	// 
	// processing Boolean values
	//

	Bool_t aDCSBool = kTRUE;

	AliDCSValue *v = 0x0;

	for(Int_t iCount = 0; iCount < array->GetEntries(); iCount++) {
		v = (AliDCSValue *)array->At(iCount);
		if (((Int_t)(v->GetTimeStamp()) < (Int_t)GetStartTimeDCSQuery()) || ((Int_t)(v->GetTimeStamp()) > (Int_t)GetEndTimeDCSQuery())) {
			AliError(Form("DCS values for the parameter outside the queried interval"));
			continue;
		}
		if (iCount > 0) {
			if (aDCSBool != v->GetBool()) {
				AliError(Form("DCS values for the parameter changed from %d to %d within the queried interval", (UInt_t)aDCSBool, (UInt_t)v->GetBool()));
				change = kTRUE;
			}
		}
		aDCSBool = v->GetBool(); // always keeping last value
		AliDebug(2,Form("Bool = %d",(Int_t)aDCSBool));
	}
	
	Char_t caDCSBool = (Char_t) aDCSBool;
	return caDCSBool;
	
}

//_______________________________________________________________

Float_t AliGRPPreprocessor::ProcessInt(const TObjArray* array)
{
	// 
	// processing Int values, returning mean
	// AliGRPObject::GetInvalidFloat() is returned if any of the DCS values
	// are outside the queried time interval or their value is out of range
	//

	Float_t aDCSArraySum = 0.0;
	Float_t aDCSArrayMean = 0.0;
	Int_t iCounts = 0;
	AliDCSValue* v = 0x0;

	for(Int_t iCount = 0; iCount < array->GetEntries(); iCount++) {
		v = (AliDCSValue *)array->At(iCount);
		if ((v->GetInt() < fminInt) || (v->GetInt() > fmaxInt)) {
			AliError(Form("Error! Int value found in DCS map at %d-th entry is OUT OF RANGE: value = %d",iCount, v->GetInt()));
			return AliGRPObject::GetInvalidFloat();
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			aDCSArraySum += v->GetInt();
			iCounts += 1;
		}
	}

	if(iCounts != 0) aDCSArrayMean = aDCSArraySum/iCounts;
	else aDCSArrayMean = AliGRPObject::GetInvalidFloat();
	
	return aDCSArrayMean;

}
//_______________________________________________________________

Float_t AliGRPPreprocessor::ProcessUInt(const TObjArray* array)
{
	// 
	// processing Int values, returning mean 
	// AliGRPObject::GetInvalidFloat() is returned if any of the DCS values
	// are outside the queried time interval or their value is out of range
	//

	Float_t aDCSArraySum = 0.0;
	Float_t aDCSArrayMean = 0.0;
	Int_t iCounts = 0;
	AliDCSValue* v = 0x0;

	for(Int_t iCount = 0; iCount < array->GetEntries(); iCount++) {
		v = (AliDCSValue *)array->At(iCount);
		if ((v->GetUInt() < fminUInt) || (v->GetUInt() > fmaxUInt)) {
			AliError(Form("Error! UInt value found in DCS map at %d-th entry is OUT OF RANGE: value = %u",iCount,v->GetUInt()));
			return AliGRPObject::GetInvalidFloat();
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			aDCSArraySum += v->GetUInt();
			iCounts += 1;
		}
	}

	if(iCounts != 0) aDCSArrayMean = aDCSArraySum/iCounts;
	else aDCSArrayMean = AliGRPObject::GetInvalidFloat();
	
	return aDCSArrayMean;

}


//_______________________________________________________________

AliDCSSensorArray *AliGRPPreprocessor::GetPressureMap(TMap* dcsAliasMap)
{
	// extract DCS pressure maps. Perform fits to save space
	
	TMap *map = fPressure->ExtractDCS(dcsAliasMap);
	if (map) {
		fPressure->MakeSplineFit(map);
		Double_t fitFraction = fPressure->NumFits()/fPressure->NumSensors(); 
		if (fitFraction > kFitFraction ) {
			AliInfo(Form("Pressure values extracted, %d fits performed.", fPressure->NumFits()));
		} else { 
			AliInfo("Too few pressure maps fitted!!!");
		}
	} else {
		AliInfo("no atmospheric pressure map extracted!!!");
	}
	delete map;
	
	return fPressure;
}


  
//_______________________________________________________________
Int_t AliGRPPreprocessor::ReceivePromptRecoParameters(UInt_t run, const char* dbHost, Int_t dbPort, const char* dbName, const char* user, const char* password, const char *cdbRoot, TString &gdc)
{
	//
	// Retrieves logbook and trigger information from the online logbook
	// This information is needed for prompt reconstruction
	//
	// Parameters are:
	// Run number
	// DAQ params: dbHost, dbPort, dbName, user, password, logbookTable, triggerTable
	// cdbRoot
	//
	// returns:
	//         positive on success: the return code is the run number of last run processed of the same run type already processed by the SHUTTLE
	//         0 on success and no run was found
	//         negative on error
	//
	// This function is NOT called during the preprocessor run in the Shuttle!
	//
	
	// defaults
	if (dbPort == 0)
		dbPort = 3306;
	
	// CDB connection
	AliCDBManager* cdb = AliCDBManager::Instance();
	cdb->SetDefaultStorage(cdbRoot);
	
	// SQL connection
	TSQLServer* server = TSQLServer::Connect(Form("mysql://%s:%d/%s", dbHost, dbPort, dbName), user, password);
	
	if (!server)
		{
			Printf("ERROR: Could not connect to DAQ LB");
			return -1;
		}
	
	// main logbook
	TString sqlQuery;
	sqlQuery.Form("SELECT DAQ_time_start, run_type, detectorMask, L3_magnetCurrent, Dipole_magnetCurrent FROM logbook WHERE run = %d", run);
	TSQLResult* result = server->Query(sqlQuery);
	if (!result)
		{
			Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
			return -2;
		}
	
	if (result->GetRowCount() == 0)
		{
			Printf("ERROR: Run %d not found", run);
			delete result;
			return -3;
		}
	
	TSQLRow* row = result->Next();
	if (!row)
		{
			Printf("ERROR: Could not receive data from run %d", run);
			delete result;
			return -4;
		}
	
	TString timeStartString(row->GetField(0));
	TString runType(row->GetField(1));
	TString detectorMaskString(row->GetField(2));
	TString l3CurrentString(row->GetField(3));
	TString dipoleCurrentString(row->GetField(4));
	time_t timeStart = (time_t)(timeStartString.Atoi());
	UInt_t detectorMask = (UInt_t)(detectorMaskString.Atoi());
	Float_t l3Current = (Float_t)(TMath::Abs(l3CurrentString.Atof()));
	Float_t dipoleCurrent = (Float_t)(TMath::Abs(dipoleCurrentString.Atof()));
	Char_t l3Polarity = (l3CurrentString.Atof() < 0) ? 1 : 0;
	Char_t dipolePolarity = (dipoleCurrentString.Atof() < 0) ? 1 : 0;
	
	AliGRPObject * grpObj = new AliGRPObject();
	grpObj->SetTimeStart(timeStart); 
	grpObj->SetRunType((TString)(row->GetField(1)));
	grpObj->SetDetectorMask(detectorMask);
	grpObj->SetL3Current(l3Current,(AliGRPObject::Stats)0);
	grpObj->SetDipoleCurrent(dipoleCurrent,(AliGRPObject::Stats)0);
	grpObj->SetL3Polarity(l3Polarity);
	grpObj->SetDipolePolarity(dipolePolarity);

	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Printf("Storing GRP/GRP/Data object with the following content");
	grpObj->Dump();
	
	AliCDBMetaData metadata;
	metadata.SetResponsible("Jan Fiete Grosse-Oetringhaus & Chiara Zampolli");
	metadata.SetComment("GRP Output parameters received during online running");
	
	AliCDBId id("GRP/GRP/Data", run, run);
	Bool_t success = cdb->Put(grpObj, id, &metadata);
	
	delete grpObj;
	
	if (!success)
		{
			Printf("ERROR: Could not store GRP/GRP/Data into OCDB");
			return -5;
		}
	
	// Receive trigger information
	sqlQuery.Form("SELECT configFile FROM logbook_trigger_config WHERE run = %d", run);
	result = server->Query(sqlQuery);
	if (!result)
		{
			Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
			return -11;
		}
	
	if (result->GetRowCount() == 0)
		{
			Printf("ERROR: Run %d not found in logbook_trigger_config", run);
			delete result;
			return -12;
		}
	
	row = result->Next();
	if (!row)
		{
			Printf("ERROR: Could not receive logbook_trigger_config data from run %d", run);
			delete result;
			return -13;
		}
	
	TString triggerConfig(row->GetField(0));
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Printf("Found trigger configuration: %s", triggerConfig.Data());
	
	AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfigurationFromString(triggerConfig);
	if (!runcfg)
		{
			Printf("ERROR: Could not create CTP configuration object");
			return -14;
		}
	
	metadata.SetComment("CTP run configuration received during online running");
	
	AliCDBId id2("GRP/CTP/Config", run, run);
	success = cdb->Put(runcfg, id2, &metadata);
	
	delete runcfg;
	runcfg = 0;
	
	if (!success)
		{
			Printf("ERROR: Could not store GRP/CTP/Config into OCDB");
			return -15;
		}
	

	// Receive list of GDCs for this run
	sqlQuery.Form("SELECT GDC FROM logbook_stats_GDC WHERE run = %d", run);
	result = server->Query(sqlQuery);
	if (!result)
		{
			Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
			return -24;
		}
	
	if (result->GetRowCount() == 0)
		{
			Printf("ERROR: Run %d not found in logbook_stats_GDC", run);
			delete result;
			return -25;
		}
	
	row = result->Next();
	if (!row)
		{
			Printf("ERROR: Could not receive logbook_stats_GDC data from run %d", run);
			delete result;
			return -26;
		}
	
	// For the moment take the first GDC in the list
	gdc = row->GetField(0);
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	Printf("Found GDC: %s", gdc.Data());

	// get last run with same run type that was already processed by the SHUTTLE
	
	sqlQuery.Form("SELECT max(logbook.run) FROM logbook LEFT JOIN logbook_shuttle ON logbook_shuttle.run = logbook.run WHERE run_type = '%s' AND shuttle_done = 1", runType.Data());
	result = server->Query(sqlQuery);
	if (!result)
		{
			Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
			return -21;
		}
	
	if (result->GetRowCount() == 0)
		{
			Printf("ERROR: No result with query <%s>", sqlQuery.Data());
			delete result;
			return -22;
		}
	
	row = result->Next();
	if (!row)
		{
			Printf("ERROR: Could not receive data for query <%s>", sqlQuery.Data());
			delete result;
			return -23;
		}
	
	TString lastRunStr(row->GetField(0));
	Int_t lastRun = lastRunStr.Atoi();
	
	Printf("Last run with same run type %s is %d", runType.Data(), lastRun);
	
	delete row;
	row = 0;
	
	delete result;
	result = 0;
	
	server->Close();
	delete server;
	server = 0;
	
	return lastRun;
}
