// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: David Rohr <drohr@kip.uni-heidelberg.de>                *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include "AliHLTCreateGRP.h"
#include <dic.hxx>
#include <stdio.h>
#include <ctime>
#include <iostream>
#include <memory>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <AliGRPObject.h>
#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliDCSSensor.h>
#include <AliMagF.h>
#include <AliDAQ.h>
using namespace std;

ClassImp(AliHLTCreateGRP)

/** Static string to define a local storage for the OCDB contact. */
static const TString kLOCAL_STORAGE_DEFINE = "local://";
/** Static define of T-HCDB */
static const char* THCDB_BASE_FOLDER = getenv("ALIHLT_T_HCDBDIR");

int AliHLTCreateGRP::CreateGRP(Int_t runNumber, TString detectorList, TString beamType, TString runType, UInt_t startShift, UInt_t endShift, Int_t defaults)
{
	int l3Polarity;
	float l3Current;
	int dipolePolarity;
	float dipoleCurrent;
	float surfaceAtmosPressure;
	float cavernAtmosPressure;
	float cavernAtmosPressure2;
	int beamEnergy;
	
	if (defaults)
	{
		printf("WARNING: Defaults mode enabled, creating GRP with some default values for tests!\n");
	
		dipoleCurrent = 5.999e+03; //Valid settings for running
		l3Current = 3.000e+04;
		dipolePolarity = 0;
		l3Polarity = 0;
		surfaceAtmosPressure = 971.127;
		cavernAtmosPressure = 974.75;
		cavernAtmosPressure2 = 975.55;
		beamEnergy = 6501;
	}
	else
	{
		//Dim does not give a real error, but we have to define some values for the case of error, so let's take something unprobable. THIS SHOULD BE IMPROVED!!!
		DimCurrentInfo solenoidPolarityDim("DCS_GRP_L3MAGNET_POLARITY", -12345678);
		DimCurrentInfo solenoidCurrentDim("DCS_GRP_L3MAGNET_CURRENT", -123456.f);
		DimCurrentInfo dipolePolarityDim("DCS_GRP_DIPOLE_POLARITY", -12345678);
		DimCurrentInfo dipoleCurrentDim("DCS_GRP_DIPOLE_CURRENT", -123456.f);
		DimCurrentInfo pressureSurfaceDim("DCS_GRP_PRESSURE_SURFACE", -123456.f);
		DimCurrentInfo pressureCavernDim("DCS_GRP_PRESSURE_CAVERN", -123456.f);
		DimCurrentInfo pressureCavernDim2("DCS_GRP_PRESSURE_CAVERN2", -123456.f);
		//These 2 are not used yet, but it might make sense to use them perhaps?
		//DimCurrentInfo LHCBeamType1Dim("DCS_GRP_LHC_BEAM_TYPE_1", -1);
		//DimCurrentInfo LHCBeamType2Dim("DCS_GRP_LHC_BEAM_TYPE_2", -1);
		//LHC Energy has not been used before, but it was seto to constant ((cmsEnergy = 14000) / 0.12), I think it is better to use dim value
		DimCurrentInfo LHCBeamEnergyDim("DCS_GRP_LHC_BEAM_ENERGY", -12345678);
		l3Polarity = solenoidPolarityDim.getInt();
		l3Current = solenoidCurrentDim.getFloat();
		dipolePolarity = dipolePolarityDim.getInt();
		dipoleCurrent = dipoleCurrentDim.getFloat();
		surfaceAtmosPressure = pressureSurfaceDim.getFloat();
		cavernAtmosPressure = pressureCavernDim.getFloat();
		cavernAtmosPressure2 = pressureCavernDim2.getFloat();
		beamEnergy = LHCBeamEnergyDim.getInt();
	}
	
	if (l3Polarity == -12345678) {printf("Error obtaining L3 Polarity from DIM\n"); return(1);}
	if (l3Current == -123456.f) {printf("Error obtaining L3 Current from DIM\n"); return(1);}
	if (dipolePolarity == -12345678) {printf("Error obtaining Dipole Polarity from DIM\n"); return(1);}
	if (dipoleCurrent == -123456.f) {printf("Error obtaining Dipole Current from DIM\n"); return(1);}
	if (surfaceAtmosPressure == -123456.f) {printf("Error obtaining Surface Pressure from DIM\n"); return(1);}
	if (cavernAtmosPressure == -123456.f) {printf("Error obtaining Cavern Pressure from DIM\n"); return(1);}
	if (cavernAtmosPressure2 == -123456.f) {printf("Error obtaining Cavern Pressure 2 from DIM\n"); return(1);}
	if (beamEnergy == -123456.f) {printf("Error obtaining Beam Energy from DIM\n"); return(1);}
	
	AliMagF* field = AliMagF::CreateFieldMap(l3Current, dipoleCurrent, 0, kFALSE, beamEnergy, beamType, "$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root", true);
	if (field == NULL)
	{
		cout << "Cannot create fild map with current magnet currents: L3 " << l3Current << ", Diploe " << dipoleCurrent;
		return(4);
	}
	
	cout << "Values from DIM: L3 Polarity " << l3Polarity << ", L3 Current " << l3Current << ", Dipole Polarity " << dipolePolarity <<
		", Dipole Current " << dipoleCurrent << ", Surface Pressure " << surfaceAtmosPressure << ", Cavern Pressure " << cavernAtmosPressure << " / " << cavernAtmosPressure2 <<
		", Beam Energy " << beamEnergy << endl;

	if (THCDB_BASE_FOLDER == NULL)
	{
		cout << "ALIHLT_T_HCDBDIR env variable must be set!" << endl;
		return 5;
	}

	TString THCDBPath = kLOCAL_STORAGE_DEFINE + THCDB_BASE_FOLDER;

	AliGRPObject grpObject;

	// set GRP values to object
	int curtime = time(NULL);
	grpObject.SetTimeStart(curtime);   // ?? not the time issued by ECS
	grpObject.SetTimeEnd(curtime + 24 * 3600);
	grpObject.SetBeamType(beamType);
	grpObject.SetRunType(runType);
	grpObject.SetLHCPeriod(getenv("LHC_PERIOD")); // ? This variable is wrong that way!

	TObjArray* setOfDetectors = detectorList.Tokenize(",");
	Char_t numOfDetectors = setOfDetectors->GetEntries();
	grpObject.SetNumberOfDetectors(numOfDetectors);
	UInt_t detectMask = createDetectorMask(setOfDetectors);
	grpObject.SetDetectorMask(detectMask);

	grpObject.SetBeamEnergy(beamEnergy);

	grpObject.SetL3Current(l3Current, (AliGRPObject::Stats) 0);
	grpObject.SetDipoleCurrent(dipoleCurrent, (AliGRPObject::Stats) 0);  
	grpObject.SetL3Polarity(l3Polarity);  
	grpObject.SetDipolePolarity(dipolePolarity);
	grpObject.SetPolarityConventionLHC();                    // LHC convention +/+ current -> -/- field main components
	grpObject.SetCavernAtmosPressure(CreateSensor("CavernAtmosPressure", cavernAtmosPressure, startShift, endShift));
	//This sensor is not yet present in DIM
	grpObject.SetCavernAtmosPressure2(CreateSensor("CavernAtmosPressure2", cavernAtmosPressure2, startShift, endShift));
	grpObject.SetSurfaceAtmosPressure(CreateSensor("SurfaceAtmosPressure", surfaceAtmosPressure, startShift, endShift));

	cout << "   ### GRP entry prepared with: " << endl;
	cout << "       start time: " << grpObject.GetTimeStart() << " (Sensor shiftTime " << startShift << " - " << endShift << ")" << endl; 
	cout << "       beam type: " << grpObject.GetBeamType().Data() << ", run type: " << grpObject.GetRunType().Data() << endl; 
	cout << "       LHC period: " << grpObject.GetLHCPeriod() << endl;
	cout << "       num of detectors: " << (Int_t) (grpObject.GetNumberOfDetectors()) << ", detector mask: " << grpObject.GetDetectorMask() << endl;
	printf("       Detector mask hex: %x\n", grpObject.GetDetectorMask());

	// Prepare storage of entry
	// init connection to HCDB
	AliCDBManager *man = AliCDBManager::Instance();
	if (man)
	{
		cout << "    --- Got a CDB Manager reference. " << endl;
	}
	else
	{
		cout << "    *** ERROR cannot obtain a CDB Manager reference." << endl;
		cout << "    *** Exiting now " << endl << endl;
		return 2;
	}
	AliCDBStorage *thcdb_storage = man->GetStorage(THCDBPath.Data());
	if (thcdb_storage)
	{
		cout << "    --- Contacted HCDB storage: " << THCDBPath.Data() << endl;
	}
	else
	{
		cout << "    *** ERROR contacting HCDB storage:" << THCDBPath.Data() << endl;
		cout << "    *** Exiting now " << endl << endl;
		return 3;
	}
	
	const char* aliroot_version = getenv("ALIROOT_VERSION");
	if (aliroot_version == NULL) aliroot_version = "unknownAliRoot";
	
	AliCDBMetaData meta("HLT", 0, aliroot_version, "GRP entry created online for HLT");
	AliCDBPath path("GRP", "GRP", "Data");
	int runmin, runmax;
	if (defaults == 2)
	{
		runmin = 0;
		runmax = 999999999;
		printf("Creating Default Object\n");
	}
	else
	{
		runmin = runmax = runNumber;
	}
	AliCDBEntry THCDBEntry(&grpObject, path, runmin, runmax, &meta);
	if (thcdb_storage->Put(&THCDBEntry))
	{
		cout << "   +++ Successfully stored GRP entry." << endl;
	}
	else
	{
		cout << "   *** ERROR storing GRP entry to HCDB." << endl << endl;
		return 6;
	}

	return(0);
}

//COPIED FROM setGrpVal.C
/**
 * Map detector list to detector mask
 *
 * @param listOfDetectors array containing the participating detectors
 *
 * @returns the detector mask encoded in an UInt_t
 */
UInt_t AliHLTCreateGRP::createDetectorMask(TObjArray* listOfDetectors)
{
	UInt_t mask = 0x00000000;

	for (Int_t i = 0; i < listOfDetectors->GetEntries(); i++) {
		TObjString* detectorName = (TObjString*) listOfDetectors->At(i);
		TString name(detectorName->String());
		if (name.CompareTo("SPD") == 0) { // detector name
			mask |= AliDAQ::kSPD;
		} else if (name.CompareTo("SDD") == 0) {
			mask |= AliDAQ::kSDD;
		} else if (name.CompareTo("SSD") == 0) {
			mask |= AliDAQ::kSSD;
		} else if (name.CompareTo("TPC") == 0) {
			mask |= AliDAQ::kTPC;
		} else if (name.CompareTo("TRD") == 0) {
			mask |= AliDAQ::kTRD;
		} else if (name.CompareTo("TOF") == 0) {
			mask |= AliDAQ::kTOF;
		} else if (name.CompareTo("HMPID") == 0) {
			mask |= AliDAQ::kHMPID;
		} else if (name.CompareTo("PHOS") == 0) {
			mask |= AliDAQ::kPHOS;
		} else if (name.CompareTo("CPV") == 0) {
			mask |= AliDAQ::kCPV;
		} else if (name.CompareTo("PMD") == 0) {
			mask |= AliDAQ::kPMD;
		} else if (name.CompareTo("MUON_TRK") == 0) {
			mask |= AliDAQ::kMUONTRK;
		} else if (name.CompareTo("MUON_TRG") == 0) {
			mask |= AliDAQ::kMUONTRG;
		} else if (name.CompareTo("FMD") == 0) {
			mask |= AliDAQ::kFMD;
		} else if (name.CompareTo("T0") == 0) {
			mask |= AliDAQ::kT0;
		} else if (name.CompareTo("V0") == 0) {
			mask |= AliDAQ::kVZERO;
		} else if (name.CompareTo("ZDC") == 0) {
			mask |= AliDAQ::kZDC;
		} else if (name.CompareTo("ACORDE") == 0) {
			mask |= AliDAQ::kACORDE;
		} else if (name.CompareTo("TRG") == 0) {
			mask |= AliDAQ::kTRG;
		} else if (name.CompareTo("EMCAL") == 0) {
			mask |= AliDAQ::kEMCAL;
		} else if (name.CompareTo("DAQ_TEST") == 0) {
			mask |= AliDAQ::kDAQTEST;
		} else if (name.CompareTo("SHUTTLE") == 0) {
			mask |= 0x20000000;
		} else if (name.CompareTo("AD") == 0) {
			mask |= AliDAQ::kAD;
		} else if (name.CompareTo("MFT") == 0) {
			mask |= AliDAQ::kMFT;
		} else if (name.CompareTo("FIT") == 0) {
			mask |= AliDAQ::kFIT;
		} else if (name.CompareTo("HLT") == 0) {
			mask |= AliDAQ::kHLT;
		} else {
			// Unknown detector names
			cout << "   *** Detector list contains unknown detector name, skipping ..." << endl;
		}
	}
	
	//Always enable HLT and CTP flags
	mask |= AliDAQ::kHLT;
	mask |= AliDAQ::kTRG;

	return mask;
}

//COPIED FROM AliHLTPredictionProcessorInterface.h
template <typename T> AliDCSSensor* AliHLTCreateGRP::CreateSensor(const char* id, T value, UInt_t starttime, UInt_t endtime) 
{
	// create an AliDCSSensor object with specified id and a linear graph
	// There is no extrapolation in the online reconstruction, only the last value
	// counts

	if (!id) return NULL;

	std::auto_ptr<AliDCSSensor> pSensor(new AliDCSSensor);
	if (!pSensor.get())
	{
		//HLTFatal("memory allocation failed");
		return NULL;
	}

	// AliDCSSensor allows two types of value representation: a spline fit and
	// a graph. The online system uses a linear graph with const values between
	// start and end time
	// Note: AliDCSSensor::GetValue returns -99 if the requested time is before
	// the start time and the last value if the time is after end time
	// The measurements are stored in fractions of hours (see definition in
	// class AliDCSSensor
	const int points=2;
	const Double_t kSecInHour = 3600.; // seconds in one hour
	T x[points];
	T y[points];
	x[0]=0;
	x[1]=(endtime-starttime)/kSecInHour;
	y[0]=value;
	y[1]=value;
	std::auto_ptr<TGraph> pGraph(new TGraph(2,x,y));
	if (!pGraph.get())
	{
		//HLTFatal("can not create graph for id %s", id);
		return NULL;
	}

	AliSplineFit *fit = new AliSplineFit();
	if (!fit) return NULL;

	fit->SetMinPoints(10);
	fit->InitKnots(pGraph.get(),10, 10, 0.0);
	fit->SplineFit(2);

	pSensor->SetStringID(id);
	pSensor->SetStartTime(starttime);
	pSensor->SetEndTime(endtime);
	//note: AliDCSSensor has no correct cleanup in the destructor
	// so the fit object is lost if the sensor is deleted
	pSensor->SetFit(fit);

	return pSensor.release();
}
