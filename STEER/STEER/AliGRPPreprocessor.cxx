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
#include <TObjArray.h>
#include <TGraph.h>
#include <TString.h>
#include <TFile.h>

#include <float.h>

#include "AliGRPPreprocessor.h"
#include "AliGRPObject.h"
#include "AliDCSSensor.h"
#include "AliSplineFit.h"
#include "AliDCSSensorArray.h"
#include "AliRawEventHeaderVersions.h"

#include "AliTriggerConfiguration.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerInput.h"

#include "AliCDBMetaData.h"
#include "AliESDVertex.h"
#include "AliLHCReader.h"
#include "AliLHCData.h"
#include "AliDCSArray.h"
#include "AliDAQ.h"
#include "AliLTUConfig.h"

class AliLog;
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
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"

const Double_t kFitFraction = -1.;                 // Fraction of DCS sensor fits required

ClassImp(AliGRPPreprocessor)

//_______________________________________________________________

  const Int_t AliGRPPreprocessor::fgknDAQLbPar = 6; // num parameters in the logbook used to fill the GRP object
  const Int_t AliGRPPreprocessor::fgknDCSDP = 48;   // number of dcs dps
  const Int_t AliGRPPreprocessor::fgknDCSDPHallProbes = 40;   // number of dcs dps
  const Int_t AliGRPPreprocessor::fgknLHCDP = 9;   // number of dcs dps from LHC data
  const Int_t AliGRPPreprocessor::fgkDCSDPHallTopShift = 4;   // shift from the top to get tp the Hall Probes names in the list of DCS DPs
  const Int_t AliGRPPreprocessor::fgkDCSDPNonWorking = 5; // number of non working DCS DPs
  const char* AliGRPPreprocessor::fgkDCSDataPoints[AliGRPPreprocessor::fgknDCSDP] = {
                   "L3Polarity",
                   "DipolePolarity",
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
                   "SurfaceAtmosPressure",
                   "CavernAtmosPressure2"
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
		   "Dipole_Outside_Temperature"
                 };
                 
  const Short_t kSensors = 45; // start index position of sensor in DCS DPs
  const Short_t kNumSensors = 3; // Number of sensors in DCS DPs (CavernAtmosPressure, SurfaceAtmosPressure, CavernAtmosPressure2)


  const char* AliGRPPreprocessor::fgkLHCDataPoints[AliGRPPreprocessor::fgknLHCDP] = {
	  "LHC_Beam_Energy",
	  "LHC_MachineMode",
	  "LHC_BeamMode",
          "LHC_Beams_Particle_Type",
	  "BPTX_Phase_Shift_B1",
	  "BPTX_Phase_Shift_B2",
	  "LHC_Particle_Type_B1",
	  "LHC_Particle_Type_B2",
	  "LHC_Data_Quality_Flag"
  };

  const char* kppError[] = {
                   "",
                   "(DAQ logbook ERROR)",
                   "(DAQ FXS ERROR)",
                   "(Trigger Scalers not found in DCS FXS - ERROR)",
                   "(DCS data points ERROR)",
                   "(Trigger Configuration ERROR)",
                   "(DAQ logbook ERROR determining partition of the run)",
                   "(CTP timing ERROR)",
		   "(SPD Mean Vertex ERROR)",
		   "(DCS FXS Error for LHC Data)",
		   "(LHC Data Error)",
		   "(LHC Clock Phase Error (from LHC Data))",
		   "(LTU Configuration Error)",
		   "(DQM Failure)"
  };

//_______________________________________________________________

AliGRPPreprocessor::AliGRPPreprocessor(AliShuttleInterface* shuttle):
	AliPreprocessor("GRP",shuttle),  fPressure(0), fmaxFloat(0), fminFloat(0),fmaxDouble(0), fminDouble(0), fmaxInt(0), fminInt(0), fmaxUInt(0), fminUInt(0),fdaqStartEndTimeOk(kTRUE),ffailedDPs(new TObjArray(fgknDCSDP))
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
        AddRunType("STANDALONE_BC");

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

	ffailedDPs->SetOwner(kTRUE);
}

//_______________________________________________________________

AliGRPPreprocessor::~AliGRPPreprocessor()
{
	//destructor
	
	delete fPressure;
	delete ffailedDPs;

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

	ffailedDPs->Clear(); // cleaning ffailedDPs for current run
	for (Int_t iDP=0; iDP < fgknDCSDP; iDP++){
		TObjString* dp = new TObjString(fgkDCSDataPoints[iDP]);
		ffailedDPs->AddAt(dp,iDP);
	}

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
	grpobj->SetBeamEnergyIsSqrtSHalfGeV(); // new format

	//=================//
	// DAQ logbook     //
	//=================//

	Log("*************** Processing DAQ logbook");

	UInt_t error = 0;
	
	Int_t iDaqLB = ProcessDaqLB(grpobj);
	TString runType = (TString)GetRunType();
	TString beamType = (TString)GetRunParameter("beamType");
	if(iDaqLB == fgknDAQLbPar) {
		Log(Form("DAQ Logbook, successful! Retrieved %d/%d entries",iDaqLB,fgknDAQLbPar));
	} else {
		Log(Form("DAQ Logbook, could not get all expected entries!!! Retrieved only %d/%d entries",iDaqLB,fgknDAQLbPar));
		error |= 1;
	}

	//=================//
	// DAQ FXS         //
	//=================//

	Log("*************** Processing DAQ FXS");

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

	Log("*************** Processing DCS FXS");

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

	Log("*************** Processing DCS DPs");

	Log(Form("Starting DCS Query at %d and finishing at %d",GetStartTimeDCSQuery(),GetEndTimeDCSQuery()));
	Int_t entries = ProcessDcsDPs( valueMap, grpobj );
	Log(Form("entries found = %d (should be %d)",entries, fgknDCSDP-fgkDCSDPNonWorking));
	if (fdaqStartEndTimeOk){
		if( entries < fgknDCSDP - fgkDCSDPNonWorking ) { // L3_BSF4_H3, L3_BSF17_H1, L3_BSF17_H2, L3_BSF17_H3, L3_BSF17_Temperature are not working yet...  
			Log(Form("Possible problem with the DCS data points!!! Only %d/%d entries found - Please read further for more details",entries,fgknDCSDP-fgkDCSDPNonWorking));
			Log(Form("The DPs giving problems were:"));
			for (Int_t iDP = 0; iDP < fgknDCSDP; iDP++){
				TObjString *dpString = (TObjString*)ffailedDPs->At(iDP);
				if (dpString){
					TString name = dpString->String();
					if (name != "L3_BSF4_H3" && name != "L3_BSF17_H1" && name != "L3_BSF17_H2" && name != "L3_BSF17_H3" && name != "L3_BSF17_Temperature" ){
						Log(Form("******** %s ******** not present, but foreseen --> causing an ERROR",name.Data()));
					}
					else {
						Log(Form(" %s is not present, but was not generating any error since it is not ready in DCS - check the other DPs in this list!",name.Data()));
					}
				}
			}
			error |= 8;
		}
		else  Log(Form("DCS data points, successful!"));
	}
	else Log(Form("Statistical values for DCS DPs could not be computed due to missing DAQ_time_start and DAQ_time_end fields in DAQ logbook")); 
	
	//=======================//
	// Trigger Configuration //
	//=======================//

	Log("*************** Processing Trigger Configuration");

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

        //===========================//
	// Trigger Timing Parameters //
        //===========================//

	Log("*************** Processing Trigger Time Params");
	
	const char * triggerCTPtiming = GetCTPTimeParams();

	if (partition.IsNull() && !detector.IsNull()){ // standalone partition
		Log("STANDALONE partition for current run, using CTP timing params dummy value");
		AliCDBEntry *cdbEntry = GetFromOCDB("CTP","DummyCTPtime");
		if (!cdbEntry) {
			Log(Form("No dummy CTP timing parameters entry found, going into error..."));
			error |= 64;
		}
		else{
			AliCTPTimeParams *runCTPtiming = (AliCTPTimeParams*)cdbEntry->GetObject();
			if (!runCTPtiming){
				Log(Form("dummy CTP timing parameters not found in OCDB entry, going into error..."));
				error |= 64;
			}
			else {
				TString titleCTPtiming = Form("CTP timing params for run %i from Dummy entry in OCDB",fRun);
				runCTPtiming->SetTitle(titleCTPtiming);
				AliCDBMetaData metadata;
				metadata.SetResponsible("Roman Lietava");
				metadata.SetComment("CTP run timing parameters from dummy entry in OCDB");
				if (!Store("CTP","CTPtiming", runCTPtiming, &metadata, 0, 0)) {
					Log("Unable to store the dummy CTP timing params object to OCDB!");
					error |= 64;
				}
			}
		}
	}

	else if (!partition.IsNull() && detector.IsNull()){ // global partition
		Log("GLOBAL partition for current run, using Trigger Timing Parameters from DAQ Logbook");
		if (triggerCTPtiming!= NULL) {
			Log("Found trigger timing params in DAQ logbook");
			AliDebug(2,Form("%s",triggerCTPtiming));
			AliCTPTimeParams *runCTPtiming = AliCTPTimeParams::LoadCTPTimeParamsFromString(triggerCTPtiming);	  
			if (!runCTPtiming) {
				Log("Bad CTP trigger timing params file from DAQ logbook! The corresponding CDB entry will not be filled!");
				error |= 64;
			}
			else {
				TString titleCTPtiming = Form("CTP timing params for run %i from DAQ",fRun);
				runCTPtiming->SetTitle(titleCTPtiming);
				AliCDBMetaData metadata;
				metadata.SetBeamPeriod(0);
				metadata.SetResponsible("Roman Lietava");
				metadata.SetComment("CTP timing params from DAQ logbook");
				if (!Store("CTP","CTPtiming", runCTPtiming, &metadata, 0, 0)) {
					Log("Unable to store the CTP timing params object to OCDB!");
					error |= 64;
				}
			}
		}

		else {
			Log("Trigger timing params NULL in DAQ logbook");
			error |= 64;
		}
	}

	else {
		Log(Form("Incorrect field in DAQ logbook for partition = %s and detector = %s, going into error without trigger timing parameters...",partition.Data(),detector.Data()));
		error |= 32;
	}

        //===========================//
	// LTU Configuration         //
        //===========================//

	Log("*************** Processing LTU Configuration");
	
	if (partition.IsNull() && !detector.IsNull()){ // standalone partition
		Log("STANDALONE partition for current run, using LTU configuration dummy value");
		AliCDBEntry *cdbEntry = GetFromOCDB("CTP","DummyLTUConfig");
		if (!cdbEntry) {
			Log(Form("No dummy LTU Config entry found, going into error..."));
			error |= 2048;
		}
		else{
			TObjArray *ltuConfig = (TObjArray*)cdbEntry->GetObject();
			if (!ltuConfig){
				Log(Form("dummy LTU Config not found in OCDB entry, going into error..."));
				error |= 2048;
			}
			else {
				AliCDBMetaData metadata;
				metadata.SetResponsible("Roman Lietava");
				metadata.SetComment("LTU Config from dummy entry in OCDB");
				if (!Store("CTP","LTUConfig", ltuConfig, &metadata, 0, 0)) {
					Log("Unable to store the dummy LTU Config object to OCDB!");
					error |= 2048;
				}
			}
		}
	}

	else if (!partition.IsNull() && detector.IsNull()){ // global partition
	
		Log("GLOBAL partition for current run, getting LTU Config from DAQ Logbook (logbook_detectors table)");
		UInt_t  detectorMask = (UInt_t)(((TString)GetRunParameter("detectorMask")).Atoi());
		Printf ("detectormask = %d",detectorMask);
		TObjArray * ltuarray = new TObjArray();
		ltuarray->SetOwner(1);
		Bool_t isLTUok = kTRUE;
		for(Int_t i = 0; i<AliDAQ::kNDetectors-2; i++){
			if ((detectorMask >> i) & 0x1) {
				TString det = AliDAQ::OfflineModuleName(i);
				TString detCTPName = AliTriggerInput::fgkCTPDetectorName[i];
				if (detCTPName == "CTP") {
					detCTPName="TRG"; // converting according to what is found in DAQ logbook_detectors					
					Printf("Processing CTP (CTP Detector name %s) --> SKIPPING, CTP does not have any LTU!!!!!!",detCTPName.Data());
					continue;
				}				
				Printf("Processing detector %s (CTP Detector name %s)",det.Data(),detCTPName.Data());
				TString* ltu = GetLTUConfig(detCTPName.Data());
				if (!ltu){
					Log(Form("No LTU Configuration from DAQ logbook for detector %s (BUT it was expected)! The corresponding CDB entry will not be filled!",detCTPName.Data()));
					error |= 2048;
					isLTUok = kFALSE;
					break;
				}
				else{
					Float_t ltuFineDelay1 = ltu[0].Atof();
					Float_t ltuFineDelay2 = ltu[1].Atof();
					Float_t ltuBCDelayAdd = ltu[2].Atof();
					const char* name = AliDAQ::DetectorName(i);
					AliLTUConfig* ltuConfig = new AliLTUConfig((UChar_t)AliDAQ::DetectorID(name),ltuFineDelay1,ltuFineDelay2,ltuBCDelayAdd);
					ltuarray->AddAtAndExpand(ltuConfig,i);
				}				
			}
		}
		if (isLTUok){
			AliCDBMetaData metadata;
			metadata.SetBeamPeriod(0);
			metadata.SetResponsible("Roman Lietava");
			metadata.SetComment("LTU Configuration for current run");
			if (!Store("CTP","LTUConfig", ltuarray, &metadata, 0, 0)) {
				Log("Unable to store the LTU Config object to OCDB!");
				error |= 2048;
			}		
		}
		delete ltuarray;
	}

	else {
		Log(Form("Incorrect field in DAQ logbook for partition = %s and detector = %s, going into error without trigger timing parameters...",partition.Data(),detector.Data()));
		error |= 32;
	}


	//=================//
	// LHC Data        //
	//=================//

	if (runType == "PHYSICS"){  // processing the LHC file only in PHYSICS runs
		Log("*************** Processing LHC Data");

		UInt_t iLHCData = ProcessLHCData(grpobj);
		
		if( iLHCData == 0 ) {
			Log(Form("LHC Data from DCS FXS, successful!"));
		} else  if (iLHCData == 1) {
			Log(Form("LHC Data, problems with DCS FXS!"));
			error |= 256;
		} else  if (iLHCData == 2) {
			Log(Form("LHC Data, problems with DAQ_time_start/DAQ_time_end!"));
			error |= 512;
		} else if (iLHCData ==3){
			Log(Form("Problems in storing LHC Phase - going into Error"));
			error |= 1024;
		} else if (iLHCData ==4){
			Log(Form("Problems with LHC Phase - going into Error"));
			error |= 1024;
		} else{
			Log(Form("LHC Data problems"));
			error |= 512;
		}
	
	}

	//==================//
	// SPD Mean Vertex  //
	//==================//

	Log("*************** Processing SPD Mean Vertex");

	if (runType == "PHYSICS"){
		UInt_t iSPDMeanVertex = ProcessSPDMeanVertex();
		if( iSPDMeanVertex == 1 ) {
			Log(Form("SPD Mean Vertex, successful!"));
		} else {
			Log(Form("SPD Mean Vertex failed!!!"));
			error |= 128; 
		}
	}
	else {
		Log("SPD Mean Vertex not processed since runType != PHYSICS");
	}

	//=================//
	// DQM FXS         //
	//=================//

	Log("*************** Processing DQM FXS");

	UInt_t iDqmFxs = ProcessDqmFxs();
	if( iDqmFxs == 1 ) {
		Log(Form("DQM FXS, successful!"));
	} else {
		Log(Form("DQM FXS failed!!!"));
		error |= 4096;
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
		Log( Form("GRP Preprocessor FAILS!!! %s%s%s%s%s%s%s%s%s%s%s%s%s",
			  kppError[(error&1)?1:0],
			  kppError[(error&2)?2:0],
			  kppError[(error&4)?3:0],
			  kppError[(error&8)?4:0],
			  kppError[(error&16)?5:0],
			  kppError[(error&32)?6:0],
			  kppError[(error&64)?7:0],
			  kppError[(error&128)?8:0],
			  kppError[(error&256)?9:0],
			  kppError[(error&512)?10:0],
			  kppError[(error&1024)?11:0],
			  kppError[(error&2048)?12:0],
			  kppError[(error&4096)?12:0]
			  ));
		return error;
	}


}

//_______________________________________________________________

UInt_t AliGRPPreprocessor::ProcessLHCData(AliGRPObject *grpobj)
{
	//
	//Getting the LHC Data from DCS FXS
	//

	TString timeStartString = (TString)GetRunParameter("DAQ_time_start");
	TString timeEndString = (TString)GetRunParameter("DAQ_time_end");
	if (timeStartString.IsNull() || timeEndString.IsNull()){
		if (timeStartString.IsNull()){ 
			AliError("DAQ_time_start not set in logbook! Setting statistical values for current DP to invalid");
		}
		else if (timeEndString.IsNull()){
			AliError("DAQ_time_end not set in logbook! Setting statistical values for current DP to invalid");
		}
		return 2;
	}  

	Double_t timeStart = timeStartString.Atof();
	Double_t timeEnd = timeEndString.Atof();

	TString fileName = GetFile(kDCS, "LHCData","");
	if (fileName.Length()>0){
		AliInfo("Got The LHC Data file");
		AliLHCReader lhcReader;

		// Processing data to be put in AliGRPObject

		// Energy
		Log("*************Energy ");
		TObjArray* energyArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[0]);
		if (energyArray){			
			Float_t energy = ProcessEnergy(energyArray,timeStart);
			if (energy != -1.) {
				grpobj->SetBeamEnergy(energy);
				grpobj->SetBeamEnergyIsSqrtSHalfGeV(kTRUE);
			}
			delete energyArray;
		}
		else {
			AliError("Energy not found in LHC Data file!!!");
		}	

		Double_t timeBeamModeEnd = timeEnd;        // max validity for Beam Mode 
		Double_t timeMachineModeEnd = timeEnd;     // max validity for Machine Mode
 		Double_t timeBeamEnd = timeEnd;            // max validity for Beam Type
 		Double_t timeBeamTypeEnd[2] = {timeEnd, timeEnd}; // max validity for Beam Type1,2
		Double_t timeBeamModeStart = -1;    // min validity for Beam Mode
		Double_t timeMachineModeStart = -1; // min validity for Machine Mode
		Double_t timeBeamStart = -1;        // min validity for Beam Type
		Double_t timeBeamTypeStart[2] = {-1,-1};        // min validity for Beam Type1,2
		Int_t indexBeamMode = -1;                  // index of measurement used to set Beam Mode
	        Int_t indexMachineMode = -1;               // index of measurement used to set Machine Mode
		Int_t indexBeam = -1;                      // index of measurement used to set Beam Type
		Int_t indexBeamType[2] = {-1, -1};                      // index of measurement used to set Beam Type1,2
		Bool_t foundBeamModeStart = kFALSE;        // flag to be set in case an entry for the Beam Mode is found before (or at) SOR
		Bool_t foundMachineModeStart = kFALSE;     // flag to be set in case an entry for the Machine Mode is found before (or at) SOR
		Bool_t foundBeamStart = kFALSE;            // flag to be set in case an entry for the Beam Type is found before (or at) SOR
		Bool_t foundBeamTypeStart[2] = {kFALSE, kFALSE};            // flag to be set in case an entry for the Beam Type1,2 is found before (or at) SOR
		Bool_t flagBeamMode = kFALSE;  //flag set true if a changed occurred in BeamMode
		Bool_t flagMachineMode = kFALSE;  //flag set true if a changed occurred in MachineMode
		Bool_t flagBeam = kFALSE;  //flag set true if a changed occurred in BeamType
		Bool_t flagBeamType[2] = {kFALSE, kFALSE};  //flag set true if a changed occurred in BeamType1,2
		
		Double_t arrayTimes[5]={2.E9, 2.E9, 2.E9, 2.E9, 2.E9}; // array to keep track of the times of the possible changes of the LHC DPs; each entry set to Wed May 18 2033, 03:33:20 GMT (ALICE should not be running anymore...)
		                                                       // arrayTimes elements order correspond to the one used in the array of the strings fgkLHCDataPoints, i.e.:
		                                                       // arrayTimes[0] --> MachineMode
		                                                       // arrayTimes[1] --> BeamMode
		                                                       // arrayTimes[2] --> BeamType (when written together)
		                                                       // arrayTimes[3] --> BeamType1 (when written separate)
		                                                       // arrayTimes[4] --> BeamType2 (when written separate)

		// BeamMode
		Log("*************BeamMode (LHCState) ");
		TObjArray* beamModeArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[2]);
		Int_t nBeamMode = -1;
		if (beamModeArray){	
			nBeamMode = beamModeArray->GetEntries();	
			if (nBeamMode==0){
				AliInfo("Found zero entries for the Beam Mode, leaving it empty");
			}
			else{
				for (Int_t iBeamMode = 0; iBeamMode<nBeamMode; iBeamMode++){
					AliDCSArray* beamMode = (AliDCSArray*)beamModeArray->At(iBeamMode);
					if (beamMode){
						if (beamMode->GetTimeStamp()<=timeStart && beamMode->GetTimeStamp()>=timeBeamModeStart){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
							timeBeamModeStart = beamMode->GetTimeStamp();
							indexBeamMode = iBeamMode;
							foundBeamModeStart = kTRUE;
						}
						else {
							break;

						}
					}
				}
				if (!foundBeamModeStart){
					AliInfo("No value for the Beam Mode found before start of run, the Beam Mode will remain empty");
				}
				else {
					AliDCSArray* beamMode = (AliDCSArray*)beamModeArray->At(indexBeamMode);
					TObjString* beamModeString = beamMode->GetStringArray(0);
					AliInfo(Form("LHC State (corresponding to BeamMode) = %s (set at %f)",(beamModeString->String()).Data(),beamMode->GetTimeStamp()));
					grpobj->SetLHCState(beamModeString->String());
					if (indexBeamMode < nBeamMode-1){
						AliDCSArray* beamMode1 = (AliDCSArray*)beamModeArray->At(indexBeamMode+1);
						if (beamMode1){
							if (beamMode1->GetTimeStamp()<=timeStart){
								AliError("you did not choose the correct value! there is still something before (or at) SOR, but later than this!");
							}
							else if (beamMode1->GetTimeStamp()>timeStart && beamMode1->GetTimeStamp()<=timeEnd){
								timeBeamModeEnd = beamMode1->GetTimeStamp();
								TObjString* beamModeString1 = beamMode1->GetStringArray(0);
							        TString bmString0 = beamModeString->String();
							        TString bmString1 = beamModeString1->String();
								if (bmString0.CompareTo(bmString1.Data(),TString::kIgnoreCase) == -1){
									AliWarning(Form("The beam mode changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",bmString0.Data(), bmString1.Data(), timeBeamModeEnd, bmString0.Data()));
									flagBeamMode = kTRUE;
									arrayTimes[1]=timeBeamModeEnd;

								}
							}
						}
						else {
							AliInfo("Invalid pointer for the first entry for Beam Mode after the first valid one, not considering anything after what has already been found");
						}
					}
				}
			}
			delete beamModeArray;
		}
		else{
			AliError("Beam mode array not found in LHC Data file!!!");
		}
		
		// MachineMode
		Log("*************MachineMode ");
		TObjArray* machineModeArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[1]);
		Int_t nMachineMode = -1;
		if (machineModeArray){
			nMachineMode = machineModeArray->GetEntries();
			if (nMachineMode==0){
				AliInfo("No Machine Mode found, leaving it empty");
			}
			else{
				for (Int_t iMachineMode = 0; iMachineMode<nMachineMode; iMachineMode++){
					AliDCSArray* machineMode = (AliDCSArray*)machineModeArray->At(iMachineMode);
					if (machineMode){
						if (machineMode->GetTimeStamp()<=timeStart && machineMode->GetTimeStamp()>=timeMachineModeStart){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
							timeMachineModeStart = machineMode->GetTimeStamp();
							indexMachineMode = iMachineMode;
							foundMachineModeStart = kTRUE;
						}
						else{
							break;
						}
					}
				}
				if (!foundMachineModeStart){
					AliInfo("No value for the Machine Mode found before start of run, the Machine Mode will remain empty");
				}
				else {
					AliDCSArray* machineMode = (AliDCSArray*)machineModeArray->At(indexMachineMode);
					TObjString* machineModeString = machineMode->GetStringArray(0);
					AliInfo(Form("MachineMode = %s (set at %f)",(machineModeString->String()).Data(),machineMode->GetTimeStamp()));
					grpobj->SetMachineMode(machineModeString->String());
					if (indexMachineMode < nMachineMode-1){
						AliDCSArray* machineMode1 = (AliDCSArray*)machineModeArray->At(indexMachineMode+1);
						if (machineMode1){
							if (machineMode1->GetTimeStamp()>timeStart && machineMode1->GetTimeStamp()<=timeEnd){
								timeMachineModeEnd = machineMode1->GetTimeStamp();
								TObjString* machineModeString1 = machineMode1->GetStringArray(0);
							        TString mmString0 = machineModeString->String();
							        TString mmString1 = machineModeString1->String();
								if (mmString0.CompareTo(mmString1.Data(),TString::kIgnoreCase) == -1){
									AliWarning(Form("The machine mode changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",mmString0.Data(),mmString1.Data(),timeMachineModeEnd,mmString0.Data()));
									flagMachineMode = kTRUE;
									arrayTimes[0]=timeMachineModeEnd;
								}
							}
						}
						else {
							AliInfo("Invalid pointer for the first entry for Machine Mode after the first valid one, not considering anything after what has already been found");
						}
					}
				}
			}
			delete machineModeArray;
		}
		else{
			AliError("Machine mode array not found in LHC Data file!!!");
		}
		
		// BeamType1 and BeamType2 - both put in the same string
		Log("*************BeamType ");
		TObjArray* beamArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[3]);
		if (beamArray){			
			Int_t nBeam = beamArray->GetEntries();
			if (nBeam==0){
				AliInfo("No Beam Type found, leaving it empty");
			}
			else{
				for (Int_t iBeam = 0; iBeam<nBeam; iBeam++){
					AliDCSArray* beam = (AliDCSArray*)beamArray->At(iBeam);
					if (beam){
						if (beam->GetTimeStamp()<=timeStart && beam->GetTimeStamp()>=timeBeamStart){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
							timeBeamStart = beam->GetTimeStamp();
							indexBeam = iBeam;
							foundBeamStart = kTRUE;
						}
						else{
							break;
						}
					}
				}
				if (!foundBeamStart){
					AliInfo("No value for the Beam Type found before start of run, the (common) Beam Type will remain empty");
				}
				else {
					AliDCSArray* beam = (AliDCSArray*)beamArray->At(indexBeam);
					TObjString* beamString = beam->GetStringArray(0);
					TString beamType = beamString->String();
					AliInfo(Form("Beam Type = %s",beamType.Data()));	
					if (beamType.CompareTo("PROTON",TString::kIgnoreCase) == 0){
						AliInfo("Setting beam type to p-p");
						grpobj->SetBeamType("p-p");
					}
					else { // if there is no PROTON beam, we suppose it is Pb, and we put A-A
						AliInfo("Setting beam type to A-A");
						grpobj->SetBeamType("A-A");
					}
					/*
					  else if (beamType.CompareTo("LEAD82",TString::kIgnoreCase) == 0){
						AliInfo("Setting beam type to Pb-Pb");
						grpobj->SetBeamType("Pb-Pb");
					}
					else{
						AliError("Beam Type not known, leaving it empty");
					}
					*/
					if (indexBeam < nBeam-1){
						AliDCSArray* beam1 = (AliDCSArray*)beamArray->At(indexBeam+1);
						if (beam1){
							if (beam1->GetTimeStamp()>timeStart && beam1->GetTimeStamp()<=timeEnd){
								timeBeamEnd = beam1->GetTimeStamp();
								TObjString* beamString1 = beam1->GetStringArray(0);
								TString beamType1 = beamString1->String();
								if (beamType.CompareTo(beamType1.Data(),TString::kIgnoreCase) == -1){
									AliWarning(Form("The Beam Type changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",beamType.Data(),(beamString1->String()).Data(),timeBeamEnd,beamType.Data()));
									flagBeam = kTRUE;
									arrayTimes[2] = timeBeamEnd;
								}
							}
						}
						else {
							AliInfo("Invalid pointer for the first entry for Beam Type after the first valid one, not considering anything after what has already been found");
						}
					}
				}
			}
			delete beamArray;
		}
		else{
			AliError("Beam Type array not found in LHC Data file!!!");
		}		

		// BeamType1 and BeamType2 - in separete string
		Log("*************BeamType, 1 and 2 ");
		Int_t indexBeamTypeString = 6;  // index of the string with the alias of BeanType1 in the array fgkLHCDataPoints
		TString combinedBeamType = "-";  // combined beam tyope, built from beam type 1 and beam type 2
		for (Int_t ibeamType = 0; ibeamType<2; ibeamType++){
			beamArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[indexBeamTypeString+ibeamType]);
			if (beamArray){			
				Int_t nBeam = beamArray->GetEntries();
				if (nBeam==0){
					AliInfo(Form("No Beam Type %s found, leaving it empty",fgkLHCDataPoints[indexBeamTypeString+ibeamType]));
				}
				else{
					for (Int_t iBeam = 0; iBeam<nBeam; iBeam++){
						AliDCSArray* beam = (AliDCSArray*)beamArray->At(iBeam);
						if (beam){
							if (beam->GetTimeStamp()<=timeStart && beam->GetTimeStamp()>=timeBeamTypeStart[ibeamType]){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
								timeBeamTypeStart[ibeamType] = beam->GetTimeStamp();
								indexBeamType[ibeamType] = iBeam;
								foundBeamTypeStart[ibeamType] = kTRUE;
							}
							else{
								break;
							}
						}
					}
					if (!foundBeamTypeStart[ibeamType]){
						AliInfo(Form("No value for the Beam Type %s found before start of run, the Beam Type %d will remain empty", fgkLHCDataPoints[indexBeamTypeString+ibeamType], ibeamType));
					}
					else {
						AliDCSArray* beam = (AliDCSArray*)beamArray->At(indexBeam);
						TObjString* beamString = beam->GetStringArray(0);
						TString beamType = beamString->String();
						AliInfo(Form("Beam Type (for %s) = %s", fgkLHCDataPoints[indexBeamTypeString+ibeamType], beamType.Data()));	
						if (beamType.CompareTo("PROTON",TString::kIgnoreCase) == 0){
							AliInfo(Form("Setting beam type %s to p", fgkLHCDataPoints[indexBeamTypeString+ibeamType]));
							grpobj->SetSingleBeamType(ibeamType,"p");
							if (ibeamType == 0) combinedBeamType.Prepend("p");
							else combinedBeamType.Append("p");
						}
						else { // if there is no PROTON beam, we suppose it is Pb, and we put A-A
							AliInfo("Setting beam type to A");
							grpobj->SetSingleBeamType(ibeamType,"A");
							if (ibeamType == 0) combinedBeamType.Prepend("A");
							else combinedBeamType.Append("A");
						}
						/*
						  else if (beamType.CompareTo("LEAD82",TString::kIgnoreCase) == 0){
						  AliInfo("Setting beam type to Pb-Pb");
						  grpobj->SetSingleBeamType(ibeamType, "Pb-Pb");
						  }
						  else{
						  AliError("Beam Type not known, leaving it empty");
						  }
						*/
						if (indexBeamType[ibeamType] < nBeam-1){
							AliDCSArray* beam1 = (AliDCSArray*)beamArray->At(indexBeam+1);
							if (beam1){
								if (beam1->GetTimeStamp()>timeStart && beam1->GetTimeStamp()<=timeEnd){
									timeBeamTypeEnd[ibeamType] = beam1->GetTimeStamp();
									TObjString* beamString1 = beam1->GetStringArray(0);
									TString beamType1 = beamString1->String();
									if (beamType.CompareTo(beamType1.Data(),TString::kIgnoreCase) == -1){
										AliWarning(Form("The Beam Type for %s changed from %s to %s during the run at timestamp %f! Setting it to %s and keeping track of the time of the change to set MaxTimeLHCValidity afterward",fgkLHCDataPoints[indexBeamTypeString+ibeamType],beamType.Data(),(beamString1->String()).Data(),timeBeamEnd,beamType.Data()));
										flagBeamType[ibeamType] = kTRUE;
										arrayTimes[3+ibeamType] = timeBeamTypeEnd[ibeamType];
									}
								}
							}
							else {
								AliInfo(Form("Invalid pointer for the first entry for Beam Type %s after the first valid one, not considering anything after what has already been found",fgkLHCDataPoints[indexBeamTypeString+ibeamType]));
							}
						}
					}
				}
				delete beamArray;
			}
			else{
				AliError(Form("Beam Type %s array not found in LHC Data file!!!",fgkLHCDataPoints[indexBeamTypeString+ibeamType]));
			}		
		}
		AliInfo(Form("Setting combined beam type to %s",combinedBeamType.Data()));
		grpobj->SetBeamType(combinedBeamType);
		
		// Setting minTimeLHCValidity
		if (flagBeamMode == kTRUE || flagMachineMode == kTRUE || flagBeam == kTRUE || flagBeamType[0] == kTRUE || flagBeamType[1] == kTRUE){ 
			Double_t minTimeLHCValidity= TMath::MinElement(5,arrayTimes);
			AliWarning(Form("Setting MaxTimeLHCValidity to %f",minTimeLHCValidity));
			grpobj->SetMaxTimeLHCValidity(minTimeLHCValidity);
		}
		/* 
		   // Old way to determine the Maximum Time during which the LHC info is valid
		if (timeBeamModeEnd!=0 || timeMachineModeEnd!=0 || timeBeamEnd !=0){
			Double_t minTimeLHCValidity;
			if (flagBeamMode == kFALSE && flagMachineMode == kFALSE && flagBeam == kTRUE){ // flagBeam only true --> it is the only one that changed
				minTimeLHCValidity = timeBeamEnd;
			}
			else if (flagBeamMode == kFALSE && flagMachineMode == kTRUE && flagBeam == kFALSE){ // flagMachineMode only true
				minTimeLHCValidity = timeMachineModeEnd;
			}
			else if (flagBeamMode == kTRUE && flagMachineMode == kFALSE && flagBeam == kFALSE){ // flagBeamMode only true
				minTimeLHCValidity = timeBeamModeEnd;
			}
			else if (flagBeamMode == kFALSE && flagMachineMode == kTRUE && flagBeam == kTRUE){ // flagBeam and flagMachineMode only true
				minTimeLHCValidity= TMath::Min(timeBeamEnd,timeMachineModeEnd);
			}
			else if (flagBeamMode == kTRUE && flagMachineMode == kFALSE && flagBeam == kTRUE){ // flagBeam and flagBeamMode only true
				minTimeLHCValidity= TMath::Min(timeBeamEnd,timeBeamModeEnd);
			}
			else if (flagBeamMode == kTRUE && flagMachineMode == kTRUE && flagBeam == kFALSE){ // flagMachineMode and flagBeamMode only true
				minTimeLHCValidity= TMath::Min(timeMachineModeEnd,timeBeamModeEnd);
			}
			else {
				Double_t arrayTimes[3] = {timeBeamModeEnd,timeMachineModeEnd,timeBeamEnd};// flagMachineMode and flagBeamMode and flagBeam 
				minTimeLHCValidity= TMath::MinElement(3,arrayTimes);
			}
			AliWarning(Form("Setting MaxTimeLHCValidity to %f",minTimeLHCValidity));
			grpobj->SetMaxTimeLHCValidity(minTimeLHCValidity);
		}
		*/
		
		// Data Quality Flag --> storing start and end values of periods within the run during which the value was found to be FALSE
		Log("*************Data Quality Flag ");
		TObjArray* dataQualityArray = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[8]);
		Int_t nDataQuality = -1;
		Double_t timeDataQualityStart = -1; // min validity for Data Quality Flag
	        Int_t indexDataQuality = -1;               // index of first measurement used to set Data Quality Flag
		Bool_t foundDataQualityStart = kFALSE;     // flag to be set in case an entry for the Data Quality Flag is found before (or at) SOR

		if (dataQualityArray){
			nDataQuality = dataQualityArray->GetEntries();
			if (nDataQuality==0){
				AliInfo("No Data Quality Flag found, leaving it empty");
			}
			else{
				for (Int_t iDataQuality = 0; iDataQuality<nDataQuality; iDataQuality++){
					AliDCSArray* dataQuality = (AliDCSArray*)dataQualityArray->At(iDataQuality);
					if (dataQuality){
						if (dataQuality->GetTimeStamp()<=timeStart && dataQuality->GetTimeStamp()>=timeDataQualityStart){// taking always the very last entry: if two measurements have the same timestamp, the last one is taken
							timeDataQualityStart = dataQuality->GetTimeStamp();
							indexDataQuality = iDataQuality;
							foundDataQualityStart = kTRUE;
						}
						else{
							// we suppose here that if the first measurement is not before SOR, then none will be (they MUST be in chronological order!!!) 
							break;
						}
					}
				}
				if (!foundDataQualityStart){
					// The Data Quality Flag should be found and TRUE at the start of the run. For the time being, if it is not found, don't do anything, but it means there is a problem..
					AliInfo("No value for the Data Quality Flag found before start of run, the Data Quality Flag will remain empty");
				}
				else {
					// counting how many FALSE values there are
					Bool_t foundEndOfFalse = kFALSE;
					Int_t nFalse = 0;
					for (Int_t iDataQuality = indexDataQuality; iDataQuality < nDataQuality; iDataQuality ++){
						AliDCSArray* dataQuality = (AliDCSArray*)dataQualityArray->At(iDataQuality);
						AliDebug(4,Form("dataQuality->GetTimeStamp() = %f, timeDataQualityStart = %f, timeEnd = %f", dataQuality->GetTimeStamp(), timeDataQualityStart, timeEnd ));
						if (dataQuality->GetTimeStamp()>=timeDataQualityStart && dataQuality->GetTimeStamp()<=timeEnd){ // considering only values between the first valid and the end of the run
							Bool_t dataQualityFlag = dataQuality->GetBool(0);
							AliDebug(3,Form("DataQuality = %d (set at %f)",(Int_t)dataQualityFlag,dataQuality->GetTimeStamp()));
							if (dataQualityFlag != kTRUE){
								if (iDataQuality == indexDataQuality) {  // the first Data Quality value should be TRUE, but ignoring the problem now...
									AliError("The first value for the Data Quality MUST be TRUE! Ignoring for now...");
								}
								nFalse++;
							}
						}
					}

					AliInfo(Form("Found %d FALSE values for the Data Quality Flag",nFalse));
					Double_t falses[nFalse*2];  // dimensioning this to the maximum possible, as if each false value was followed by a true one --> the false periods correspond to the number of falses

					Int_t iDataQuality = indexDataQuality;
					if (nFalse > 0){
						Int_t iFalse = 0;
						// filling the info about the periods when the flag was set to FALSE
						// starting, like for the other DPS, from the measurement closest to SOR (the index of which is iDataQuality)
						while (iDataQuality < nDataQuality){
							AliDebug(3,Form("iDataQuality = %d",iDataQuality));
							AliDCSArray* dataQuality = (AliDCSArray*)dataQualityArray->At(iDataQuality);
							if (dataQuality->GetTimeStamp()>=timeDataQualityStart && dataQuality->GetTimeStamp()<=timeEnd){ // considering only values between the first valid and the end of the run
								Bool_t dataQualityFlag = dataQuality->GetBool(0);
								AliDebug(3,Form("DataQuality = %d (set at %f)",(Int_t)dataQualityFlag,dataQuality->GetTimeStamp()));
								if (dataQualityFlag == kTRUE){
									// found TRUE value, continuing
									iDataQuality++;
									continue;
								}
								else{
									/*
									// the check was already done before
									if (iDataQuality == indexDataQuality) {  // the first Data Quality value should be TRUE, but ignoring the problem now...
									AliError("The first value for the Data Quality MUST be TRUE! Ignoring for now...");
									}
									*/
									falses[iFalse*2] = dataQuality->GetTimeStamp();
									foundEndOfFalse = kFALSE;
									Int_t iDataQualityNext = iDataQuality+1;
									while (iDataQualityNext < nDataQuality){
										AliDCSArray* dataQualityNext = (AliDCSArray*)dataQualityArray->At(iDataQualityNext);
										if (dataQualityNext->GetTimeStamp()>timeDataQualityStart && dataQualityNext->GetTimeStamp()<=timeEnd && dataQualityNext->GetTimeStamp() > dataQuality->GetTimeStamp()){ // considering only values between the first valid and the end of the run, and subsequent to the current value
											Bool_t dataQualityFlagNext = dataQualityNext->GetBool(0);
											AliDebug(3,Form("DataQualityNext = %d (set at %f)",(Int_t)dataQualityFlagNext,dataQualityNext->GetTimeStamp()));
											if (dataQualityFlagNext == kTRUE){
												// found TRUE value, first FALSE period completed
												foundEndOfFalse = kTRUE;
												falses[iFalse*2+1] = dataQualityNext->GetTimeStamp();
												iFalse++;
												break;
											}
											iDataQualityNext++;
										}
									}
									if (!foundEndOfFalse) {
										AliInfo("Please, note that the last FALSE value lasted until the end of the run");
										falses[iFalse*2+1] = timeEnd;
										iFalse++;
										break;
									}
									iDataQuality = iDataQualityNext+1;
								}
							}
						}
						grpobj->SetNFalseDataQualityFlag(iFalse);
						grpobj->SetFalseDataQualityFlagPeriods(falses);
					}
				}
			}
			delete dataQualityArray;
		}
		else{
			AliError("Data Quality Flag array not found in LHC Data file!!!");
		}

		// Processing data to go to AliLHCData object
		AliLHCData* dt = new AliLHCData(fileName.Data(),timeStart,timeEnd);
		// storing AliLHCData in OCDB
		if (dt){
		  AliInfo(Form("Filled %d records to AliLHCData object",dt->GetData().GetEntriesFast()));
		  AliCDBMetaData md;
		  md.SetResponsible("Ruben Shahoyan");
		  md.SetComment("LHC data from the GRP preprocessor.");
		  Bool_t result = kTRUE;
		  result = Store("GRP", "LHCData", dt, &md); 
		  delete dt;
		  if (!result){
			Log(Form("Problems in storing LHC Data - but not going into Error"));
		  }
		}

		// processing LHC Phase

		TObjArray *beam1phase = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[4]);
		TObjArray *beam2phase = lhcReader.ReadSingleLHCDP(fileName.Data(),fgkLHCDataPoints[5]);
		if (beam1phase == 0x0 || beam2phase == 0x0){
			Log(Form("Problems in retrieving LHC Clock data from LHC file"));
			return 4;
		}			
		AliLHCClockPhase *phaseObj = ProcessLHCClockPhase(beam1phase,beam2phase,timeEnd);
		delete beam1phase;
		delete beam2phase;
		if (phaseObj){
			AliInfo(Form("LHC Phase found"));
			AliCDBMetaData mdPhase;
			mdPhase.SetResponsible("Cvetan Cheshkov");
			mdPhase.SetComment("LHC Clock Phase");
			Bool_t result = kTRUE;
			result = Store("Calib", "LHCClockPhase", phaseObj, &mdPhase); 
			delete phaseObj;
			if (!result) return 3;
		}
		else return 4;
	}
	
	else {
		AliError("No LHCData file found in DCS FXS");
		return 1;
	}

	return 0;
}

//_______________________________________________________________

UInt_t AliGRPPreprocessor::ProcessSPDMeanVertex()
{
	//
	//Getting the SPD Mean Vertex
	//

	TList* list = GetForeignFileSources("SPD", kDAQ, "VertexDiamond");
	Bool_t storeResult = kTRUE;
	if (list !=0x0 && list->GetEntries()!=0)
		{
			AliInfo("The following sources produced files with the id VertexDiamond from SPD");
			list->Print();
			for (Int_t jj=0;jj<list->GetEntries();jj++){
				TObjString * str = dynamic_cast<TObjString*> (list->At(jj)); 
				if (!str){
					AliError(Form("Expecting a TObjString in the list for the %d-th source, but something else was found.",jj));
					continue;
				}
				AliInfo(Form("found source %s", str->String().Data()));
				TString fileNameRun = GetForeignFile("SPD", kDAQ, "VertexDiamond", str->GetName());
				if (fileNameRun.Length()>0){
					AliInfo(Form("Got the file %s", fileNameRun.Data()));
					TFile daqFile(fileNameRun.Data(),"READ");
					if (daqFile.IsOpen()) {
						AliESDVertex* meanVtx = dynamic_cast<AliESDVertex*>(daqFile.Get("MeanVertexPos"));
						if (meanVtx){
							meanVtx->Print();	
							// storing in the OCDB 
							AliCDBMetaData md;
							md.SetResponsible("Cvetan Cheshkov");
							md.SetComment("SPD Mean Vertex");					
							storeResult = Store("Calib", "MeanVertexSPD", meanVtx, &md, 0, kTRUE); 
						}
						else{
							AliWarning("No SPD Mean Vertex object found in file");
						} 
					}
					else {
						AliError("Can't open file");
						storeResult = kFALSE;
					}
				}
				else{
					AliWarning("No file found for current source for SPD Mean Vertex");
				}
			}
		}
	else {
		AliWarning("No list found for SPD Mean Vertex");
	}

	if (list) delete list;

	return storeResult;
}
//_______________________________________________________________

UInt_t AliGRPPreprocessor::ProcessDqmFxs()
{
	//
	// Processing DQM fxs information
	//

	TList* list = GetFileSources(kDQM, "TriggerClassesAndHistosToClone");
	Bool_t storeResult = kTRUE;
	if (list !=0x0 && list->GetEntries()!=0){
		AliInfo("The following sources produced files with the id TriggerClassesAndHistosToClone for GRP");
		list->Print();
		for (Int_t jj=0;jj<list->GetEntries();jj++){
			TObjString * str = dynamic_cast<TObjString*> (list->At(jj)); 
			if (!str){
				AliError(Form("Expecting a TObjString in the list for the %d-th source, but something else was found.",jj));
				continue;
			}
			AliInfo(Form("found source %s", str->String().Data()));
			TString fileNameRun = GetFile(kDQM, "TriggerClassesAndHistosToClone", str->GetName());
			if (fileNameRun.Length()>0){
				AliInfo(Form("Got the file %s", fileNameRun.Data()));
				TFile dqmFile(fileNameRun.Data(),"READ");
				if (dqmFile.IsOpen()) {
					dqmFile.ls();
				}
				else {
					AliError("Can't open file");
					storeResult = kFALSE;
				}
			}
			else{
				AliWarning("No file found for current source for DQM TriggerClassesAndHistosToClone");
			}
		}
	}
	else {
		AliWarning("No list found for DQM TriggerClassesAndHistosToClone");
	}
	
	if (list) delete list;
	
	//	return storeResult;
	return kTRUE;  // temporary!!
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

	if (timeEnd >= 2.E9) AliFatal("ALICE run finshed later than Wed May 18 2033, 03:33:20 GMT, maximum time allowed for LHC data --> fix the GRP preprocessor!!!");

	UInt_t nparameter = 0;
	if (timeStart != 0){
		grpObj->SetTimeStart(timeStart);
		Log(Form("Start time for run %d: %d",fRun, (Int_t)timeStart));
		nparameter++;
	} 
	else {
		Log(Form("Start time not put in logbook, setting to invalid in GRP entry, and causing an error!"));
	}

	if (timeEnd != 0){
		grpObj->SetTimeEnd(timeEnd);
		Log(Form("End time for run %d: %i",fRun, (Int_t)timeEnd));
		nparameter++;
	} 
	else {
		Log(Form("End time not put in logbook, setting to invalid in GRP entry, and causing an error!"));
	}

	if (beamEnergy != 0){
		Log(Form("Beam Energy for run %d: %f (NOT USING IT TO FILL THE GRP OBJECT, taking it from the LHC file)",fRun, beamEnergy));
	} 
	else {
		Log(Form("Beam Energy not put in logbook, but not using it anyway for the GRP object (taking it from the LHC file)"));
	}

		
	if (beamType.Length() != 0){
		Log(Form("Beam Type for run %d: %s (NOT USING IT TO FILL THE GRP OBJECT, taking it from the LHC file)",fRun, beamType.Data()));
	} 
	else {
		Log(Form("Beam Type not put in logbook, but not using it anyway for the GRP entry (taking it from the LHC file)"));
	}
		
	if (numberOfDetectors != 0){
		grpObj->SetNumberOfDetectors(numberOfDetectors);
		Log(Form("Number Of Detectors for run %d: %d",fRun, (Int_t)numberOfDetectors));
		nparameter++;
	} 
	else {
		Log(Form("Number Of Detectors not put in logbook, setting to invalid in GRP entry, and causing an error!"));
	}

	if (detectorMask != 0){
		grpObj->SetDetectorMask(detectorMask);
		Log(Form("Detector Mask for run %d: %d",fRun, detectorMask));
		nparameter++;
	} 
	else {
		Log(Form("Detector Mask not put in logbook, setting to invalid in GRP entry, and causing an error!"));
	}

	if (lhcPeriod.Length() != 0) {
		grpObj->SetLHCPeriod(lhcPeriod);
		Log(Form("LHC period (DAQ) for run %d: %s",fRun, lhcPeriod.Data()));
		nparameter++;
	} 
	else {
		Log(Form("LHCperiod not put in logbook, setting to invalid in GRP entry, and causing an error!"));
	}
	if (runType.Length() != 0) {
		grpObj->SetRunType(runType);
		Log(Form("Run Type (DAQ) for run %d: %s",fRun, runType.Data()));
		nparameter++;
	} 
	else {
		Log(Form("Run Type not put in logbook, setting to invalid in GRP entry, and causing an error!"));
	}

	return nparameter;
}
//_______________________________________________________________

UInt_t AliGRPPreprocessor::ProcessDaqFxs()
{
	//======DAQ FXS======//
	
	AliRawEventHeaderV3_9::Class()->IgnoreTObjectStreamer(); // to avoid trying reading TObject store in AliRawEventHeaderV3_9 - temporary fix 
	AliRawEventHeaderV3_11::Class()->IgnoreTObjectStreamer(); // to avoid trying reading TObject store in AliRawEventHeaderV3_11 - temporary fix 
	AliRawEventHeaderV3_12::Class()->IgnoreTObjectStreamer(); // to avoid trying reading TObject store in AliRawEventHeaderV3_12 - temporary fix 
	Log("Processing DAQ FXS");
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
				if (idStr->String().CompareTo("QAThreshold") == 0 || idStr->String().CompareTo("TriggerClassesAndHistosToClone") == 0) {
					Log(Form("Skipping file with Id %s",idStr->String().Data()));
					continue; 
				}
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

	if (nFiles == 0){
		Log("no raw data tags in this run: it could be that one or more files were found in the DAQ FXS, but they were ignored, since not interesting for the raw data tag: nothing to merge!");
		delete iter;
		iter = 0;
		delete list;
		list = 0;
		delete fRawTagChain; 
		fRawTagChain=0;
		return 0;
	}
	
	TString fRawDataFileName = "GRP_Merged.tag.root";
	Log(Form("Merging %d raw data tags into file: %s", nFiles, fRawDataFileName.Data()));
	
	if (fRawTagChain->Merge(fRawDataFileName) < 1 ) {
		Log(Form("Error merging %d raw data files!!!",nFiles));
		delete iter;
		iter = 0;
		delete list;
		list = 0;
		delete fRawTagChain; 
		fRawTagChain=0;
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
	iter = 0;
	delete list;
	list = 0;
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
					delete scalers;					
					return 1;
				}
			}
			delete scalers;
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
	Int_t nL3Entries = 0;
	Int_t nDipoleEntries = 0;
	Int_t nEnvEntries = 0;
	Int_t nHallProbesEntries = 0;
	nL3Entries = ProcessL3DPs(valueMap, grpObj);
	nDipoleEntries = ProcessDipoleDPs(valueMap, grpObj);
	nEnvEntries = ProcessEnvDPs(valueMap, grpObj);
	nHallProbesEntries = ProcessHPDPs(valueMap, grpObj);
	grpObj->SetPolarityConventionLHC();  // after the dipole cables swap we comply with LHC convention
	Log(Form("L3Entries = %d, nDipoleEntries =%d, nEnvEntries = %d, nHallProbesEntries = %d", nL3Entries, nDipoleEntries, nEnvEntries, nHallProbesEntries));
	entries = nL3Entries + nDipoleEntries + nEnvEntries + nHallProbesEntries;
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
		if (!outOfRange) {
			nL3Entries++;
			ffailedDPs->RemoveAt(indexDP);
		}
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
				ffailedDPs->RemoveAt(indexDP);
				nL3Entries++;
			}
			else if (isZero){
				AliInfo(Form("%s set to invalid, but magnet was OFF (according to the current), DP not considered wrong",fgkDCSDataPoints[indexDP]));
				ffailedDPs->RemoveAt(indexDP);
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
		if (!outOfRange) {
			nDipoleEntries++;
			ffailedDPs->RemoveAt(indexDP);
		}
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
				ffailedDPs->RemoveAt(indexDP);
				nDipoleEntries++;
			}
			else if (isZero){
				AliInfo(Form("%s set to invalid, but magnet was OFF (according to the current), DP not considered wrong",fgkDCSDataPoints[indexDP]));
				ffailedDPs->RemoveAt(indexDP);
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
		if (!outOfRange) {
			ffailedDPs->RemoveAt(indexDP);
			nEnvEntries++;
		}
	}

	if (array) array = 0x0;

	AliInfo(Form("========== AtmosPressures (Cavern + Surface + Cavern2) ==========="));
	AliDCSSensorArray *dcsSensorArray = GetPressureMap(valueMap);
	//dcsSensorArray->Print();
	if( fPressure->NumFits()<kNumSensors ) {
		Log(Form("Check the pressure sensor values! Not all the %d pressure sensors have been fit",kNumSensors));
	} 
	Log(Form("Number of fits performed = %d",fPressure->NumFits()));

	AliInfo(Form("==========CavernAtmosPressure==========="));
	indexDP = kCavernAtmosPressure;
	AliDCSSensor* sensorCavernP2 = dcsSensorArray->GetSensor(fgkDCSDataPoints[indexDP]);
	TGraph* graph = sensorCavernP2->GetGraph();
	AliDebug(3,Form("index = %d",indexDP));
	AliDebug(3,Form("name = %s",fgkDCSDataPoints[indexDP]));
	AliDebug(2,Form("graph = %p",graph));
	AliDebug(3,Form("sensorCavernP2 = %p", sensorCavernP2));
	if(sensorCavernP2->GetFit() || graph) {
		if (sensorCavernP2->GetFit()){
			Log(Form("Fit for sensor %s found",fgkDCSDataPoints[indexDP]));
		}
		else {
			Log(Form("Fit for sensor %s not found, but the graph is there - NOT going into error",fgkDCSDataPoints[indexDP]));
		}
		grpObj->SetCavernAtmosPressure(sensorCavernP2);
		ffailedDPs->RemoveAt(indexDP);
		nEnvEntries++;
	} 
	//if (sensorP2) delete sensorP2;
	else {
		Log(Form("ERROR!!! Neither graph nor fit found for sensor %s - this will not increase the number of found DCS DPs and will cause an error", fgkDCSDataPoints[indexDP] ));
	}
	
	AliInfo(Form("==========SurfaceAtmosPressure==========="));
	indexDP = kSurfaceAtmosPressure;
	AliDCSSensor* sensorP2 = dcsSensorArray->GetSensor(fgkDCSDataPoints[indexDP]);
	graph = sensorP2->GetGraph();
	AliDebug(3,Form("index = %d",indexDP));
	AliDebug(3,Form("name = %s",fgkDCSDataPoints[indexDP]));
	AliDebug(2,Form("graph = %p",graph));	
	AliDebug(3,Form("sensorP2 = %p", sensorP2));
	if(sensorP2->GetFit() || graph) {
		if (sensorP2->GetFit()){
			Log(Form("Fit for sensor %s found",fgkDCSDataPoints[indexDP]));
		}
		else {
			Log(Form("Fit for sensor %s not found, but the graph is there - NOT going into error",fgkDCSDataPoints[indexDP]));
		}
		grpObj->SetSurfaceAtmosPressure(sensorP2);
		ffailedDPs->RemoveAt(indexDP);
		nEnvEntries++;
	} 
	//if (sensorP2) delete sensorP2;
	else {
		Log(Form("ERROR!!! Neither graph nor fit found for sensor %s - this will not increase the number of found DCS DPs and will cause an error", fgkDCSDataPoints[indexDP] ));
	}

	AliInfo(Form("==========CavernAtmosPressure2==========="));
	indexDP = kCavernAtmosPressure2;
	AliDCSSensor* sensorCavernP22 = dcsSensorArray->GetSensor(fgkDCSDataPoints[indexDP]);
	graph = sensorCavernP22->GetGraph();
	AliDebug(3,Form("index = %d",indexDP));
	AliDebug(3,Form("name = %s",fgkDCSDataPoints[indexDP]));
	AliDebug(2,Form("graph = %p",graph));	
	AliDebug(3,Form("sensorCavernP2_2 = %p", sensorCavernP22));
	if(sensorCavernP22->GetFit() || graph) {
		if (sensorCavernP22->GetFit()){
			Log(Form("Fit for sensor %s found",fgkDCSDataPoints[indexDP]));
		}
		else {
			Log(Form("Fit for sensor %s not found, but the graph is there - NOT going into error",fgkDCSDataPoints[indexDP]));
		}
		grpObj->SetCavernAtmosPressure2(sensorCavernP22);
		ffailedDPs->RemoveAt(indexDP);
		nEnvEntries++;
	} 
	//if (sensorP2) delete sensorP2;
	else {
		Log(Form("ERROR!!! Neither graph nor fit found for sensor %s - this will not increase the number of found DCS DPs and will cause an error", fgkDCSDataPoints[indexDP] ));
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
			if (!outOfRange) {
				ffailedDPs->RemoveAt(indexDP + fgkDCSDPHallTopShift);  // 7 = shift in the complete list of DPs to get to the Hall Probes
				nHPEntries++;
			}
		}
	}
		
	Log(Form("Hall Probes = %d ", nHPEntries));
	return nHPEntries;
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

	TString timeStartString = (TString)GetRunParameter("DAQ_time_start");
	TString timeEndString = (TString)GetRunParameter("DAQ_time_end");
	if (timeStartString.IsNull() || timeStartString.IsNull()){
		if (timeStartString.IsNull()){ 
			AliError("DAQ_time_start not set in logbook! Setting statistical values for current DP to invalid");
		}
		else if (timeStartString.IsNull()){
			AliError("DAQ_time_end not set in logbook! Setting statistical values for current DP to invalid");
		}
		fdaqStartEndTimeOk = kFALSE;
		return 0;
	}  

	Int_t timeStart = (Int_t)(timeStartString.Atoi());
	Int_t timeEnd = (Int_t)(timeEndString.Atoi());
	Float_t* parameters = new Float_t[5];
	Int_t iCounts = 0;
	Int_t iCountsRun = 0;
	Int_t nCounts = array->GetEntries();
	Float_t valueBeforeSOR = 0;
	Float_t valueAfterEOR = 0;
	Int_t timestampBeforeSOR = -1;
	Int_t timestampAfterEOR = -1;
	Int_t ientrySOR = -1;
	Int_t ientryEOR = -1;
	Float_t* arrayValues = 0x0; 
	Double_t* arrayWeights = 0x0; 
	Bool_t truncMeanFlag = kTRUE;  // flag to indicate whether Truncated Mean should be calculated or not
	Bool_t sdFlag = kTRUE;  // flag to indicate whether SD (wrt Mean/Median) should be calculated or not

	for(Int_t i = 0; i < nCounts; i++) {
		AliDCSValue *v = (AliDCSValue *)array->At(i);
		if ((v->GetFloat() <= fminFloat) || (v->GetFloat() >= fmaxFloat)) {
			AliError(Form("Error! Float value found in DCS map at %d-th entry is OUT OF RANGE: value = %6.5e",i,v->GetFloat()));
			if (v->GetFloat() < fminFloat) AliInfo(Form("The value is smaller than %6.5e",fminFloat));
			if (v->GetFloat() > fmaxFloat) AliInfo(Form("The value is greater than %6.5e",fmaxFloat));
			delete [] parameters;
			return NULL;
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			AliDebug(2,Form("%d-th entry = %f at timestamp %i",i,v->GetFloat(),v->GetTimeStamp()));
			iCounts += 1;
			// look for the last value before SOR and the first value before EOR
			if (((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) && (Int_t)(v->GetTimeStamp()) < timeStart) {
				timestampBeforeSOR = (Int_t)(v->GetTimeStamp());
				AliDebug(2,Form("timestamp of last value before SOR = %d, with DAQ_time_start = %d",timestampBeforeSOR,timeStart));
				valueBeforeSOR = v->GetFloat();
			}
			else if ((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery() && (Int_t)(v->GetTimeStamp()) > timeEnd && timestampAfterEOR == -1){
				timestampAfterEOR = (Int_t)(v->GetTimeStamp());
				valueAfterEOR = v->GetFloat();
				AliDebug(2,Form("timestamp of first value after EOR = %d, with DAQ_time_end = %d",timestampAfterEOR,timeEnd));
			}
			// check if there are DPs between DAQ_time_start and DAQ_time_end
			if(((Int_t)(v->GetTimeStamp()) >= timeStart) &&((Int_t)(v->GetTimeStamp()) <= timeEnd)) {
				if (ientrySOR == -1) ientrySOR = i;  // first entry after SOR
				if (ientryEOR < i) ientryEOR = i;  // last entry before EOR
				AliDebug(2,Form("entry between SOR and EOR"));
				iCountsRun += 1;
			}
		}
		else {
			AliError(Form("DCS values for the parameter outside the queried interval: timestamp = %d",v->GetTimeStamp()));
		}
	}

	if (timestampBeforeSOR == -1){
		AliWarning("No value found before SOR");
	}
	if (timestampAfterEOR == -1){
		AliWarning("No value found after EOR");
	}

	AliDebug(2,Form("Number of valid entries (within DCS query interval) = %i, from a total amount of %i entries",iCounts,nCounts));
	AliDebug(2,Form("Last value before DAQ_time_start (SOR) = %f at timestamp = %d",valueBeforeSOR,timestampBeforeSOR));
	AliDebug(2,Form("First value after DAQ_time_end (EOR)   = %f at timestamp = %d",valueAfterEOR,timestampAfterEOR));
	AliInfo(Form("Found %d entries between DAQ_time_start (SOR) and DAQ_time_end (EOR)",iCountsRun));
	AliDebug(2,Form("Index of first entry after DAQ_time_start (SOR) = %d ",ientrySOR));
	AliDebug(2,Form("Index of first entry before DAQ_time_end (EOR) = %d ",ientryEOR));

	Int_t nentriesUsed = 0;
	if (iCountsRun > 1){
		AliInfo("Using entries between DAQ_time_start (SOR) and DAQ_time_end (EOR)");
		AliDebug(2,"Calculating (weighted) Mean and Median");
		arrayValues = new Float_t[iCountsRun]; 
		arrayWeights = new Double_t[iCountsRun]; 
		nentriesUsed = iCountsRun;
		for (Int_t i = ientrySOR; i <= ientryEOR; i++){
			AliDCSValue *v = (AliDCSValue *)array->At(i);
			Int_t timestamp2 = 0;
			if (i < ientryEOR){
				AliDCSValue *v1 = (AliDCSValue *)array->At(i+1);
				timestamp2 = (Int_t)v1->GetTimeStamp();
			}
			else {
				timestamp2 = timeEnd+1;
			}
			arrayWeights[i-ientrySOR] = (Double_t)(timestamp2 - (Int_t)v->GetTimeStamp());
			arrayValues[i-ientrySOR] = v->GetFloat();
		}
		parameters[0] = TMath::Mean(iCountsRun,arrayValues,arrayWeights);
		parameters[2] = TMath::Median(iCountsRun,arrayValues,arrayWeights);
	}
	else if (iCountsRun == 1){
		AliDCSValue* v = (AliDCSValue *)array->At(ientrySOR);
		nentriesUsed = 2;
		if (timestampBeforeSOR != -1 && timestampBeforeSOR != (Int_t)v->GetTimeStamp()){
			AliWarning("Using single entry between DAQ_time_start (SOR) and DAQ_time_end (EOR) and last entry before SOR. Truncated mean won't be calculated.");
			arrayValues = new Float_t[2];
			arrayWeights = new Double_t[2];
			arrayValues[0] = valueBeforeSOR;
			arrayWeights[0] = (Double_t)((Int_t)v->GetTimeStamp()-timestampBeforeSOR);
			arrayValues[1] = v->GetFloat();
			arrayWeights[1] = (Double_t)(timeEnd+1-(Int_t)v->GetTimeStamp());
			AliDebug(2, Form("value0 = %f, with weight = %f",arrayValues[0],arrayWeights[0])); 
			AliDebug(2, Form("value1 = %f, with weight = %f",arrayValues[1],arrayWeights[1])); 
			parameters[0] = TMath::Mean(2,arrayValues,arrayWeights);
			parameters[2] = TMath::Median(2,arrayValues,arrayWeights);
			truncMeanFlag = kFALSE;
		}
		else{
			AliError("Cannot calculate mean, truncated mean, median, SD wrt mean, SD wrt median for current DP - only one value collected during the run, but no value before with which to calculate the statistical quantities");
			parameters[0] = AliGRPObject::GetInvalidFloat();
			parameters[1] = AliGRPObject::GetInvalidFloat();
			parameters[2] = AliGRPObject::GetInvalidFloat();
			parameters[3] = AliGRPObject::GetInvalidFloat();
			parameters[4] = AliGRPObject::GetInvalidFloat();
			return parameters;
		}
	}
	else { // iCountsRun == 0, using only the point immediately before SOR
		if (timestampBeforeSOR == -1){
			AliError("Cannot set mean, truncated mean, median, SD wrt mean, SD wrt median for current DP - no points during the run collected, and point before SOR missing");
			parameters[0] = AliGRPObject::GetInvalidFloat();
			parameters[1] = AliGRPObject::GetInvalidFloat();
			parameters[2] = AliGRPObject::GetInvalidFloat();
			parameters[3] = AliGRPObject::GetInvalidFloat();
			parameters[4] = AliGRPObject::GetInvalidFloat();
			return parameters;
		}
		else {
			AliWarning("Using only last entry before SOR. Truncated mean and Standard deviations (wrt mean/median) won't be calculated.");
			AliDebug(2,Form("value = %f",valueBeforeSOR)); 
			parameters[0] = valueBeforeSOR;
			parameters[2] = valueBeforeSOR;
			truncMeanFlag = kFALSE;
			sdFlag = kFALSE;
		}
	}

	Float_t temp = 0;
	Float_t temp1 = 0;
	Float_t sumweights = 0; 
	Int_t entriesTruncMean = 0;
	Float_t* arrayValuesTruncMean = new Float_t[nentriesUsed]; 
	Double_t* arrayWeightsTruncMean = new Double_t[nentriesUsed]; 

	// calculating SD wrt Mean and Median
	AliDebug(2,"Calculating SD wrt Mean and SD wrt Median");
	if (sdFlag){
		for (Int_t i =0; i< nentriesUsed; i++){
			AliDebug(2,Form("Entry %d: value = %f, weight = %f",i,arrayValues[i],arrayWeights[i]));
			temp += (arrayValues[i]-parameters[2])*(arrayValues[i]-parameters[2]);
			temp1 += arrayWeights[i]*(arrayValues[i]-parameters[0])*(arrayValues[i]-parameters[0]);
			sumweights += arrayWeights[i];
		}
		// setting SD wrt Mean 
		if (sumweights != 0 ){
			parameters[3] = TMath::Sqrt(temp1/sumweights);
		}
		else {
			AliError("Sum of weights to calculate Standard Deviation (wrt mean) <= 0, setting the SD to invalid");
			parameters[3] = AliGRPObject::GetInvalidFloat();
		}
		// setting SD wrt Median
		if (nentriesUsed != 0){
			parameters[4] = TMath::Sqrt(temp/nentriesUsed);
		}
		else{
			AliError("Number of entries used to calculate Standard Deviation (wrt median) <= 0, setting the SD to invalid");
			parameters[4] = AliGRPObject::GetInvalidFloat();
		}
	}
	else {
		parameters[3] = AliGRPObject::GetInvalidFloat();
		parameters[4] = AliGRPObject::GetInvalidFloat();
	}		

	// calculating truncated mean (this comes afterwards since you need the SD wrt Mean)
	if (truncMeanFlag){
		AliDebug(2,"Calculating Truncated Mean");
		for (Int_t i =0; i< nentriesUsed; i++){
			AliDebug(2,Form("Entry %d: value = %f, weight = %f",i,arrayValues[i],arrayWeights[i]));
			if ((arrayValues[i]<=parameters[0]+3*parameters[3]) && (arrayValues[i]>=parameters[0]-3*parameters[3])){
				arrayValuesTruncMean[entriesTruncMean]=arrayValues[i];
				arrayWeightsTruncMean[entriesTruncMean]=arrayWeights[i];
				AliDebug(2,Form("For Truncated Mean: Entry %d: value = %f, weight = %f",entriesTruncMean,arrayValuesTruncMean[entriesTruncMean],arrayWeightsTruncMean[entriesTruncMean]));
				entriesTruncMean++;			
			}
			else{
				AliDebug(2,"Discarding entry");
			}
		}
		// setting truncated mean 
		if (entriesTruncMean >1){
			AliDebug(2,Form("%d entries used for truncated mean",entriesTruncMean));
			parameters[1] = TMath::Mean(entriesTruncMean,arrayValuesTruncMean,arrayWeightsTruncMean);
		}
		else{	
			AliDebug(2,Form("Too few entries (%d) to calculate truncated mean",entriesTruncMean));
			parameters[1] = AliGRPObject::GetInvalidFloat();
		}
	}
	else{
			parameters[1] = AliGRPObject::GetInvalidFloat();
	}

	if (arrayValues){
		delete [] arrayValues;
	} 
	if (arrayWeights){
		delete [] arrayWeights;
	} 
	delete [] arrayValuesTruncMean;
	delete [] arrayWeightsTruncMean;

	AliInfo(Form("(weighted) mean = %f ",parameters[0]));
	AliInfo(Form("(weighted) truncated mean = %f ",parameters[1]));
	AliInfo(Form("median = %f ",parameters[2]));
	AliInfo(Form("(weighted) standard deviation with (weighted) mean = %f ",parameters[3]));
	AliInfo(Form("standard deviation with median = %f ",parameters[4]));
	
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

	Int_t nCounts = array->GetEntries();
	for(Int_t i = 0; i < nCounts; i++) {
		AliDCSValue *v = (AliDCSValue *)array->At(i);
		if ((v->GetFloat() <= fminFloat) || (v->GetFloat() >= fmaxFloat)) {
			AliError(Form("Error! Float value found in DCS map at %d-th entry is OUT OF RANGE: value = %6.5e",i,v->GetFloat()));
			if (v->GetFloat() < fminFloat) AliInfo(Form("The value is smaller than %6.5e",fminFloat));
			if (v->GetFloat() > fmaxFloat) AliInfo(Form("The value is greater than %6.5e",fmaxFloat));
			return NULL;
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			AliDebug(2,Form("%d-th entry = %f",i,v->GetFloat()));
			if (indexDP == kL3Current && v->GetFloat() > 350 && isZero == kTRUE) isZero=kFALSE; 
			if (indexDP == kDipoleCurrent && v->GetFloat() > 450 && isZero == kTRUE) isZero=kFALSE; 
		}
		else {
			AliError(Form("DCS values for the parameter outside the queried interval"));
		}
	}

	return ProcessFloatAll(array);
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

	TString timeStartString = (TString)GetRunParameter("DAQ_time_start");
	TString timeEndString = (TString)GetRunParameter("DAQ_time_end");
	if (timeStartString.IsNull() || timeStartString.IsNull()){
		if (timeStartString.IsNull()){ 
			AliError("DAQ_time_start not set in logbook! Setting statistical values for current DP to invalid");
		}
		else if (timeStartString.IsNull()){
			AliError("DAQ_time_end not set in logbook! Setting statistical values for current DP to invalid");
		}
		return 0;
	}  

	Int_t timeStart = (Int_t)(timeStartString.Atoi());
	Int_t timeEnd = (Int_t)(timeEndString.Atoi());
	Float_t aDCSArrayMean = 0.0;
	Int_t iCounts = 0;
	Float_t valueBeforeSOR = 0;
	Float_t valueAfterEOR = 0;
	Int_t timestampBeforeSOR = -1;
	Int_t timestampAfterEOR = -1;
	Int_t ientrySOR = -1;
	Int_t ientryEOR = -1;
	Float_t* arrayValues = 0x0; 
	Double_t* arrayWeights = 0x0; 
	Int_t iCountsRun = 0;
	Int_t nCounts = array->GetEntries();

	for(Int_t i = 0; i < nCounts; i++) {
		AliDCSValue* v = (AliDCSValue *)array->At(i);
		if ((v->GetInt() < fminInt) || (v->GetInt() > fmaxInt)) {
			AliError(Form("Error! Int value found in DCS map at %d-th entry is OUT OF RANGE: value = %d",i, v->GetInt()));
			return AliGRPObject::GetInvalidFloat();
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			AliDebug(2,Form("%d-th entry = %d at timestamp %i",i,v->GetInt(),v->GetTimeStamp()));
			iCounts += 1;
			// look for the last value before SOR and the first value before EOR
			if (((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) && (Int_t)(v->GetTimeStamp()) < timeStart) {
				timestampBeforeSOR = (Int_t)(v->GetTimeStamp());
				AliDebug(2,Form("timestamp of last entry before SOR = %d, with DAQ_time_start = %d",timestampBeforeSOR,timeStart));
				valueBeforeSOR = (Float_t) v->GetInt();
			}
			else if ((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery() && (Int_t)(v->GetTimeStamp()) > timeEnd && timestampAfterEOR == -1){
				timestampAfterEOR = (Int_t)(v->GetTimeStamp());
				valueAfterEOR = (Float_t) v->GetInt();
				AliDebug(2,Form("timestamp of first entry after EOR = %d, with DAQ_time_end = %d",timestampAfterEOR,timeEnd));
			}
			// check if there are DPs between DAQ_time_start and DAQ_time_end
			if(((Int_t)(v->GetTimeStamp()) >= timeStart) &&((Int_t)(v->GetTimeStamp()) <= timeEnd)) {
				if (ientrySOR == -1) ientrySOR = i;  // first entry after SOR
				if (ientryEOR < i) ientryEOR = i;  // last entry before EOR
				AliDebug(2,Form("entry between SOR and EOR"));
				iCountsRun += 1;
			}
		}
		else {
			AliError(Form("DCS values for the parameter outside the queried interval: timestamp = %d",v->GetTimeStamp()));
		}
	}

	if (timestampBeforeSOR == -1){
		AliWarning("No value found before SOR!");
	}
	if (timestampAfterEOR == -1){
		AliWarning("No value found after EOR!");
	}

	AliDebug(2,Form("Number of valid entries (within query interval) = %i, starting from %i entries",iCounts,nCounts));
	AliDebug(2,Form("Last value before DAQ_time_start (SOR) = %f at timestamp = %d",valueBeforeSOR,timestampBeforeSOR));
	AliDebug(2,Form("First value after DAQ_time_end (EOR)   = %f at timestamp = %d",valueAfterEOR,timestampAfterEOR));
	AliInfo(Form("Found %d entries between DAQ_time_start (SOR) and DAQ_time_end (EOR)",iCountsRun));
	AliDebug(2,Form("Index of first entry after DAQ_time_start (SOR) = %d ",ientrySOR));
	AliDebug(2,Form("Index of first entry before DAQ_time_end (EOR) = %d ",ientryEOR));

	Int_t nentriesUsed = 0;
	if (iCountsRun > 1){
		AliInfo("Using entries between DAQ_time_start (SOR) and DAQ_time_end (EOR)");
		AliDebug(2,"Calculating (weighted) Mean");
		arrayValues = new Float_t[iCountsRun]; 
		arrayWeights = new Double_t[iCountsRun]; 
		nentriesUsed = iCountsRun;
		for (Int_t i = ientrySOR; i <= ientryEOR; i++){
			AliDCSValue *v = (AliDCSValue *)array->At(i);
			Int_t timestamp2 = 0;
			if (i < ientryEOR){
				AliDCSValue *v1 = (AliDCSValue *)array->At(i+1);
				timestamp2 = (Int_t)v1->GetTimeStamp();
			}
			else {
				timestamp2 = timeEnd+1;
			}
			arrayWeights[i-ientrySOR] = (Double_t)(timestamp2 - (Int_t)v->GetTimeStamp());
			arrayValues[i-ientrySOR] = (Float_t)v->GetInt();
		}
		aDCSArrayMean = TMath::Mean(iCountsRun,arrayValues,arrayWeights);
		delete [] arrayValues;
		delete [] arrayWeights;
	}
	else if (iCountsRun == 1){
		AliDCSValue* v = (AliDCSValue *)array->At(ientrySOR);
		nentriesUsed = 2;
		if (timestampBeforeSOR != -1 && timestampBeforeSOR != (Int_t)v->GetTimeStamp()){
			AliWarning("Using single entry between DAQ_time_start (SOR) and DAQ_time_end (EOR) and last entry before SOR.");
			arrayValues = new Float_t[2];
			arrayWeights = new Double_t[2];
			arrayValues[0] = valueBeforeSOR;
			arrayWeights[0] = (Double_t)((Int_t)v->GetTimeStamp()-timestampBeforeSOR);
			arrayValues[1] = (Float_t)v->GetInt();
			arrayWeights[1] = (Double_t)(timeEnd+1-(Int_t)v->GetTimeStamp());
			AliDebug(2,Form("value0 = %f, with weight = %f",arrayValues[0],arrayWeights[0])); 
			AliDebug(2,Form("value1 = %f, with weight = %f",arrayValues[1],arrayWeights[1])); 
			aDCSArrayMean = TMath::Mean(2,arrayValues,arrayWeights);
			delete [] arrayValues;
			delete [] arrayWeights;
		}
		else{
			AliError("Cannot calculate mean - only one value collected during the run, but no value before with which to calculate the statistical quantities");
			return AliGRPObject::GetInvalidFloat();
		}
	}
	else { // iCountsRun == 0, using the point immediately before SOR and the one immediately after EOR
		if (timestampBeforeSOR == -1 || timestampAfterEOR == -1){
			if (timestampBeforeSOR == -1){
				AliError("Cannot calculate mean - no points during the run collected, and point before SOR missing");
			}
			if (timestampAfterEOR == -1){
				AliError("Cannot calculate maen - no points during the run collected, and point after EOR missing");
			}
			return AliGRPObject::GetInvalidFloat();
		}
		else {
			AliWarning("Using last entry before SOR and first entry after EOR.");
			nentriesUsed = 2;
			arrayValues = new Float_t[2];
			arrayWeights = new Double_t[2];
			arrayValues[0] = valueBeforeSOR;
			arrayWeights[0] = (Double_t)(timestampAfterEOR - timestampBeforeSOR);
			arrayValues[1] = valueAfterEOR;
			arrayWeights[1] = 1.;
			AliDebug(2,Form("value0 = %f, with weight = %f",arrayValues[0],arrayWeights[0])); 
			AliDebug(2,Form("value1 = %f, with weight = %f",arrayValues[1],arrayWeights[1])); 
			aDCSArrayMean = TMath::Mean(1,arrayValues,arrayWeights);
			delete [] arrayValues;
			delete [] arrayWeights;
		}
	}

	AliInfo(Form("mean = %f ", aDCSArrayMean));
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

	TString timeStartString = (TString)GetRunParameter("DAQ_time_start");
	TString timeEndString = (TString)GetRunParameter("DAQ_time_end");
	if (timeStartString.IsNull() || timeStartString.IsNull()){
		if (timeStartString.IsNull()){ 
			AliError("DAQ_time_start not set in logbook! Setting statistical values for current DP to invalid");
		}
		else if (timeStartString.IsNull()){
			AliError("DAQ_time_end not set in logbook! Setting statistical values for current DP to invalid");
		}
		return 0;
	}  

	Int_t timeStart = (Int_t)(timeStartString.Atoi());
	Int_t timeEnd = (Int_t)(timeEndString.Atoi());
	Float_t aDCSArrayMean = 0.0;
	Int_t iCounts = 0;
	Float_t valueBeforeSOR = 0;
	Float_t valueAfterEOR = 0;
	Int_t timestampBeforeSOR = -1;
	Int_t timestampAfterEOR = -1;
	Int_t ientrySOR = -1;
	Int_t ientryEOR = -1;
	Float_t* arrayValues = 0x0; 
	Double_t* arrayWeights = 0x0; 
	Int_t iCountsRun = 0;
	Int_t nCounts = array->GetEntries();

	for(Int_t i = 0; i < nCounts; i++) {
		AliDCSValue* v = (AliDCSValue *)array->At(i);
		if ((v->GetUInt() < fminUInt) || (v->GetUInt() > fmaxUInt)) {
			AliError(Form("Error! UInt value found in DCS map at %d-th entry is OUT OF RANGE: value = %u",i,v->GetUInt()));
			return AliGRPObject::GetInvalidFloat();
		}
		if(((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) &&((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery())) {
			AliDebug(2,Form("%d-th entry = %d at timestamp %i",i,v->GetUInt(),v->GetTimeStamp()));
			iCounts += 1;
			// look for the last value before SOR and the first value before EOR
			if (((Int_t)(v->GetTimeStamp()) >= (Int_t)GetStartTimeDCSQuery()) && (Int_t)(v->GetTimeStamp()) < timeStart) {
				timestampBeforeSOR = (Int_t)(v->GetTimeStamp());
				AliDebug(2,Form("timestamp of last entry before SOR = %d, with DAQ_time_start = %d",timestampBeforeSOR,timeStart));
				valueBeforeSOR = (Float_t)v->GetUInt();
			}
			else if ((Int_t)(v->GetTimeStamp()) <= (Int_t)GetEndTimeDCSQuery() && (Int_t)(v->GetTimeStamp()) > timeEnd && timestampAfterEOR == -1){
				timestampAfterEOR = (Int_t)(v->GetTimeStamp());
				valueAfterEOR = (Float_t)v->GetUInt();
				AliDebug(2,Form("timestamp of first entry after EOR = %d, with DAQ_time_end = %d",timestampAfterEOR,timeEnd));
			}
			// check if there are DPs between DAQ_time_start and DAQ_time_end
			if(((Int_t)(v->GetTimeStamp()) >= timeStart) &&((Int_t)(v->GetTimeStamp()) <= timeEnd)) {
				if (ientrySOR == -1) ientrySOR = i;  // first entry after SOR
				if (ientryEOR < i) ientryEOR = i;  // last entry before EOR
				AliDebug(2,Form("entry between SOR and EOR"));
				iCountsRun += 1;
			}
		}
		else {
			AliError(Form("DCS values for the parameter outside the queried interval: timestamp = %d",v->GetTimeStamp()));
		}
	}

	if (timestampBeforeSOR == -1){
		AliWarning("No value found before SOR!");
	}
	if (timestampAfterEOR == -1){
		AliWarning("No value found after EOR!");
	}

	AliDebug(2,Form("Number of valid entries (within query interval) = %i, starting from %i entries",iCounts,nCounts));
	AliDebug(2,Form("Last value before DAQ_time_start (SOR) = %f at timestamp = %d",valueBeforeSOR,timestampBeforeSOR));
	AliDebug(2,Form("First value after DAQ_time_end (EOR)   = %f at timestamp = %d",valueAfterEOR,timestampAfterEOR));
	AliInfo(Form("Found %d entries between DAQ_time_start (SOR) and DAQ_time_end (EOR)",iCountsRun));
	AliDebug(2,Form("Index of first entry after DAQ_time_start (SOR) = %d ",ientrySOR));
	AliDebug(2,Form("Index of first entry before DAQ_time_end (EOR) = %d ",ientryEOR));

	Int_t nentriesUsed = 0;
	if (iCountsRun > 1){
		AliInfo("Using entries between DAQ_time_start (SOR) and DAQ_time_end (EOR)");
		AliDebug(2,"Calculating (weighted) Mean");
		arrayValues = new Float_t[iCountsRun]; 
		arrayWeights = new Double_t[iCountsRun]; 
		nentriesUsed = iCountsRun;
		for (Int_t i = ientrySOR; i <= ientryEOR; i++){
			AliDCSValue *v = (AliDCSValue *)array->At(i);
			Int_t timestamp2 = 0;
			if (i < ientryEOR){
				AliDCSValue *v1 = (AliDCSValue *)array->At(i+1);
				timestamp2 = (Int_t)v1->GetTimeStamp();
			}
			else {
				timestamp2 = timeEnd+1;
			}
			arrayWeights[i-ientrySOR] = (Double_t)(timestamp2 - (Int_t)v->GetTimeStamp());
			arrayValues[i-ientrySOR] = (Float_t)v->GetUInt();
		}
		aDCSArrayMean = TMath::Mean(iCountsRun,arrayValues,arrayWeights);
		delete [] arrayValues;
		delete [] arrayWeights;
	}
	else if (iCountsRun == 1){
		AliDCSValue* v = (AliDCSValue *)array->At(ientrySOR);
		nentriesUsed = 2;
		if (timestampBeforeSOR != -1 && timestampBeforeSOR != (Int_t)v->GetTimeStamp()){
			AliWarning("Using single entry between DAQ_time_start (SOR) and DAQ_time_end (EOR) and last entry before SOR.");
			arrayValues = new Float_t[2];
			arrayWeights = new Double_t[2];
			arrayValues[0] = valueBeforeSOR;
			arrayWeights[0] = (Double_t)((Int_t)v->GetTimeStamp()-timestampBeforeSOR);
			arrayValues[1] = (Float_t)v->GetUInt();
			arrayWeights[1] = (Double_t)(timeEnd+1-(Int_t)v->GetTimeStamp());
			AliDebug(2,Form("value0 = %f, with weight = %f",arrayValues[0],arrayWeights[0])); 
			AliDebug(2,Form("value1 = %f, with weight = %f",arrayValues[1],arrayWeights[1])); 
			aDCSArrayMean = TMath::Mean(2,arrayValues,arrayWeights);
			delete [] arrayValues;
			delete [] arrayWeights;
		}
		else{
			AliError("Cannot calculate mean - only one value collected during the run, but no value before with which to calculate the statistical quantities");
			return AliGRPObject::GetInvalidFloat();
		}
	}
	else { // iCountsRun == 0, using the point immediately before SOR and the one immediately after EOR
		if (timestampBeforeSOR == -1 || timestampAfterEOR == -1){
			if (timestampBeforeSOR == -1){
				AliError("Cannot calculate mean - no points during the run collected, and point before SOR missing");
			}
			if (timestampAfterEOR == -1){
				AliError("Cannot calculate maen - no points during the run collected, and point after EOR missing");
			}
			return AliGRPObject::GetInvalidFloat();
		}
		else {
			AliWarning("Using last entry before SOR and first entry after EOR.");
			nentriesUsed = 2;
			arrayValues = new Float_t[2];
			arrayWeights = new Double_t[2];
			arrayValues[0] = valueBeforeSOR;
			arrayWeights[0] = (Double_t)(timestampAfterEOR - timestampBeforeSOR);
			arrayValues[1] = valueAfterEOR;
			arrayWeights[1] = 1.;
			AliDebug(2,Form("value0 = %f, with weight = %f",arrayValues[0],arrayWeights[0])); 
			AliDebug(2,Form("value1 = %f, with weight = %f",arrayValues[1],arrayWeights[1])); 
			aDCSArrayMean = TMath::Mean(1,arrayValues,arrayWeights);
			delete [] arrayValues;
			delete [] arrayWeights;
		}
	}

	AliInfo(Form("mean = %f ",aDCSArrayMean));
	return aDCSArrayMean;

}


//_______________________________________________________________

AliDCSSensorArray *AliGRPPreprocessor::GetPressureMap(TMap* dcsAliasMap)
{
	// extract DCS pressure maps. Perform fits to save space
	
	TMap *map = fPressure->ExtractDCS(dcsAliasMap);
	if (map) {
		AliDebug(2,Form("Map has %d entries",map->GetEntries()));
		fPressure->MakeSplineFit(map);
		Double_t fitFraction = fPressure->NumFits()/fPressure->NumSensors(); 
		if (fitFraction > kFitFraction ) {
			AliInfo(Form("Pressure values extracted, %d fits performed for %d sensors.", fPressure->NumFits(),fPressure->NumSensors()));
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
	sqlQuery.Form("SELECT DAQ_time_start, run_type, detectorMask, L3_magnetCurrent, Dipole_magnetCurrent,beamType FROM logbook WHERE run = %d", run);
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
	TString beamTypeString(row->GetField(5));
	time_t timeStart = (time_t)(timeStartString.Atoi());
	UInt_t detectorMask = (UInt_t)(detectorMaskString.Atoi());
	Float_t l3Current = (Float_t)(TMath::Abs(l3CurrentString.Atof()));
	Float_t dipoleCurrent = (Float_t)(TMath::Abs(dipoleCurrentString.Atof()));
	Char_t l3Polarity = (l3CurrentString.Atof() < 0) ? 1 : 0;
	Char_t dipolePolarity = (dipoleCurrentString.Atof() < 0) ? 1 : 0;
	if (beamTypeString.CompareTo("Pb-Pb",TString::kIgnoreCase) == 0){
		beamTypeString="A-A";
	}
	
	AliGRPObject * grpObj = new AliGRPObject();
	grpObj->SetTimeStart(timeStart); 
	grpObj->SetRunType((TString)(row->GetField(1)));
	grpObj->SetDetectorMask(detectorMask);
	grpObj->SetL3Current(l3Current,(AliGRPObject::Stats)0);
	grpObj->SetDipoleCurrent(dipoleCurrent,(AliGRPObject::Stats)0);
	grpObj->SetL3Polarity(l3Polarity);
	grpObj->SetDipolePolarity(dipolePolarity);
	grpObj->SetPolarityConventionLHC();  // after the dipole cables swap we comply with LHC convention
	grpObj->SetBeamType(beamTypeString);

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

	gdc = "";
	for (Int_t iGDC = 0; iGDC < result->GetRowCount(); iGDC++) {
	  row = result->Next();
	  if (!row)
	    {
	      Printf("ERROR: Could not receive logbook_stats_GDC data from run %d", run);
	      delete result;
	      return -26;
	    }
	  gdc += row->GetField(0);
	  gdc += " ";
	}

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
//------------------------------------------------------------------------------------------------------
Float_t AliGRPPreprocessor::ProcessEnergy(TObjArray* const array, Double_t timeStart){

	//
	// Method to processo LHC Energy information
	// Only the first value is returned, provided that it is withing DAQ_time_start and DAQ_time_end
	//

	Int_t nCounts = array->GetEntries();
	Float_t energy = -1;
	Double_t timeEnergy = -1;
	Int_t indexEnergy = -1;
	Bool_t foundEnergy = kFALSE;

	AliDebug(2,Form("Energy measurements = %d\n",nCounts));
	if (nCounts ==0){
		AliWarning("No Energy values found! Beam Energy remaining invalid!");
	}
	else{
		for (Int_t i = 0; i < nCounts; i++){
			AliDCSArray *dcs = (AliDCSArray*)array->At(i);
			if (dcs){
				if (dcs->GetTimeStamp()<=timeStart && dcs->GetTimeStamp()>=timeEnergy){// taking always the very last entry: of two measurements have the same timestamp, the last one is taken
					timeEnergy = dcs->GetTimeStamp();
					indexEnergy = i;
					foundEnergy = kTRUE;
				}
				else{
					break;
				}
			}
		}
		if (!foundEnergy){
			AliInfo("No value for the Energy found before start of run, the Energy will remain invalid");
		}
		else {
			AliDCSArray* dcs = (AliDCSArray*)array->At(indexEnergy);
			energy = (Float_t)(TMath::Nint(((Double_t)(dcs->GetInt(0)))*120/1000)); // sqrt(s)/2 energy in GeV
			AliInfo(Form("Energy value found = %d (at %f), converting --> sqrt(s)/2 = %f (GeV)", dcs->GetInt(0),dcs->GetTimeStamp(),energy));
		}
	}

	return energy;
}
//------------------------------------------------------------------------------------------------------
AliLHCClockPhase* AliGRPPreprocessor::ProcessLHCClockPhase(TObjArray *beam1phase,TObjArray *beam2phase, Double_t timeEnd)
{
  //
  // Method to process LHC-Clock Phase data
  // Only the values between DAQ_time_start and DAQ_time_end are kept
  //
  AliLHCClockPhase *phaseObj = new AliLHCClockPhase;

  Bool_t foundBeam1Phase = kFALSE, foundBeam2Phase = kFALSE;
  const Float_t threshold = 0.050; // we store the measurement only in case they differ with more 50ps from the previous one 

  TString timeCreatedStr = GetRunParameter("time_created");
  Double_t timeCreated = timeCreatedStr.Atof();

  Int_t nCounts = beam1phase->GetEntries();
  AliDebug(2,Form("Beam1 phase measurements = %d\n",nCounts));
  if (nCounts ==0){
    AliWarning("No beam1 LHC clock phase values found!");
    delete phaseObj;
    return NULL;
  }
  else{
    Double_t prevPhase = 0;
    for (Int_t i = 0; i < nCounts; i++){
      AliDCSArray *dcs = (AliDCSArray*)beam1phase->At(i);
      if (dcs){
	      //if (dcs->GetTimeStamp()>=timeStart && dcs->GetTimeStamp()<=timeEnd) {
	      if (dcs->GetTimeStamp()>=timeCreated && dcs->GetTimeStamp()<=timeEnd) {
	  if ((i == 0) || (i == (nCounts-1)) ||
	      !foundBeam1Phase ||
	      (TMath::Abs(dcs->GetDouble(0)-prevPhase) > threshold)) {
	    prevPhase = dcs->GetDouble(0);
	    foundBeam1Phase = kTRUE;
	    AliInfo(Form("B1 Clk Phase = %f at TS = %f",
			 (Float_t)dcs->GetDouble(0),dcs->GetTimeStamp()));  
	    phaseObj->AddPhaseB1DP((UInt_t)dcs->GetTimeStamp(),(Float_t)dcs->GetDouble(0));
	  }
	}
      }
    }
    if (!foundBeam1Phase){
      AliError("No beam1 LHC clock phase values found within the run!");
      delete phaseObj;
      return NULL;
    }
  }

  nCounts = beam2phase->GetEntries();
  AliDebug(2,Form("Beam2 phase measurements = %d\n",nCounts));
  if (nCounts ==0){
    AliWarning("No beam2 LHC clock phase values found!");
    delete phaseObj;
    return NULL;
  }
  else{
    Double_t prevPhase = 0;
    for (Int_t i = 0; i < nCounts; i++){
      AliDCSArray *dcs = (AliDCSArray*)beam2phase->At(i);
      if (dcs){
	if (dcs->GetTimeStamp()>=timeCreated && dcs->GetTimeStamp()<=timeEnd) {
	  if ((i == 0) || (i == (nCounts-1)) ||
	      !foundBeam2Phase ||
	      (TMath::Abs(dcs->GetDouble(0)-prevPhase) > threshold)) {
	    prevPhase = dcs->GetDouble(0);
	    foundBeam2Phase = kTRUE;
	    AliInfo(Form("B2 Clk Phase = %f at TS = %f",
			 (Float_t)dcs->GetDouble(0),dcs->GetTimeStamp()));  
	    phaseObj->AddPhaseB2DP((UInt_t)dcs->GetTimeStamp(),(Float_t)dcs->GetDouble(0));
	  }
	}
      }
    }
    if (!foundBeam2Phase){
      AliError("No beam2 LHC clock phase values found within the run!");
      delete phaseObj;
      return NULL;
    }
  }

  return phaseObj;
}
