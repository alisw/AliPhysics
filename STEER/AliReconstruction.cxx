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
//                                                                           //
// class for running the reconstruction                                      //
//                                                                           //
// Clusters and tracks are created for all detectors and all events by       //
// typing:                                                                   //
//                                                                           //
//   AliReconstruction rec;                                                  //
//   rec.Run();                                                              //
//                                                                           //
// The Run method returns kTRUE in case of successful execution.             //
//                                                                           //
// If the input to the reconstruction are not simulated digits but raw data, //
// this can be specified by an argument of the Run method or by the method   //
//                                                                           //
//   rec.SetInput("...");                                                    //
//                                                                           //
// The input formats and the corresponding argument are:                     //
// - DDL raw data files: directory name, ends with "/"                       //
// - raw data root file: root file name, extension ".root"                   //
// - raw data DATE file: DATE file name, any other non-empty string          //
// - MC root files     : empty string, default                               //
//                                                                           //
// By default all events are reconstructed. The reconstruction can be        //
// limited to a range of events by giving the index of the first and the     //
// last event as an argument to the Run method or by calling                 //
//                                                                           //
//   rec.SetEventRange(..., ...);                                            //
//                                                                           //
// The index -1 (default) can be used for the last event to indicate no      //
// upper limit of the event range.                                           //
//                                                                           //
// In case of raw-data reconstruction the user can modify the default        //
// number of events per digits/clusters/tracks file. In case the option      //
// is not used the number is set 1. In case the user provides 0, than        //
// the number of events is equal to the number of events inside the          //
// raw-data file (i.e. one digits/clusters/tracks file):                     //
//                                                                           //
//   rec.SetNumberOfEventsPerFile(...);                                      //
//                                                                           //
//                                                                           //
// The name of the galice file can be changed from the default               //
// "galice.root" by passing it as argument to the AliReconstruction          //
// constructor or by                                                         //
//                                                                           //
//   rec.SetGAliceFile("...");                                               //
//                                                                           //
// The local reconstruction can be switched on or off for individual         //
// detectors by                                                              //
//                                                                           //
//   rec.SetRunLocalReconstruction("...");                                   //
//                                                                           //
// The argument is a (case sensitive) string with the names of the           //
// detectors separated by a space. The special string "ALL" selects all      //
// available detectors. This is the default.                                 //
//                                                                           //
// The reconstruction of the primary vertex position can be switched off by  //
//                                                                           //
//   rec.SetRunVertexFinder(kFALSE);                                         //
//                                                                           //
// The tracking and the creation of ESD tracks can be switched on for        //
// selected detectors by                                                     //
//                                                                           //
//   rec.SetRunTracking("...");                                              //
//                                                                           //
// Uniform/nonuniform field tracking switches (default: uniform field)       //
//                                                                           //
//   rec.SetUniformFieldTracking(); ( rec.SetUniformFieldTracking(kFALSE); ) //
//                                                                           //
// The filling of additional ESD information can be steered by               //
//                                                                           //
//   rec.SetFillESD("...");                                                  //
//                                                                           //
// Again, for both methods the string specifies the list of detectors.       //
// The default is "ALL".                                                     //
//                                                                           //
// The call of the shortcut method                                           //
//                                                                           //
//   rec.SetRunReconstruction("...");                                        //
//                                                                           //
// is equivalent to calling SetRunLocalReconstruction, SetRunTracking and    //
// SetFillESD with the same detector selecting string as argument.           //
//                                                                           //
// The reconstruction requires digits or raw data as input. For the creation //
// of digits and raw data have a look at the class AliSimulation.            //
//                                                                           //
// The input data of a detector can be replaced by the corresponding HLT     //
// data by calling (usual detector string)                                   //
// SetUseHLTData("...");                                                     //
//                                                                           //
// For debug purposes the method SetCheckPointLevel can be used. If the      //
// argument is greater than 0, files with ESD events will be written after   //
// selected steps of the reconstruction for each event:                      //
//   level 1: after tracking and after filling of ESD (final)                //
//   level 2: in addition after each tracking step                           //
//   level 3: in addition after the filling of ESD for each detector         //
// If a final check point file exists for an event, this event will be       //
// skipped in the reconstruction. The tracking and the filling of ESD for    //
// a detector will be skipped as well, if the corresponding check point      //
// file exists. The ESD event will then be loaded from the file instead.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TArrayF.h>
#include <TFile.h>
#include <TList.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TGeoManager.h>
#include <TLorentzVector.h>
#include <TArrayS.h>
#include <TArrayD.h>
#include <TObjArray.h>
#include <TMap.h>

#include "AliReconstruction.h"
#include "AliCodeTimer.h"
#include "AliReconstructor.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliRawEventHeaderBase.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDfriend.h"
#include "AliESDVertex.h"
#include "AliESDcascade.h"
#include "AliESDkink.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliMultiplicity.h"
#include "AliTracker.h"
#include "AliVertexer.h"
#include "AliVertexerTracks.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDPmdTrack.h"

#include "AliESDTagCreator.h"
#include "AliAODTagCreator.h"

#include "AliGeomManager.h"
#include "AliTrackPointArray.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliAlignObj.h"

#include "AliCentralTrigger.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliCTPRawStream.h"

#include "AliQADataMakerRec.h" 
#include "AliGlobalQADataMaker.h" 
#include "AliQA.h"
#include "AliQADataMakerSteer.h"

#include "AliPlaneEff.h"

#include "AliSysInfo.h" // memory snapshots
#include "AliRawHLTManager.h"

#include "AliMagWrapCheb.h"

ClassImp(AliReconstruction)


//_____________________________________________________________________________
const char* AliReconstruction::fgkDetectorName[AliReconstruction::fgkNDetectors] = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "ACORDE", "HLT"};

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const char* gAliceFilename,
				     const char* name, const char* title) :
  TNamed(name, title),

  fUniformField(kFALSE),
  fForcedFieldMap(0x0),
  fRunVertexFinder(kTRUE),
  fRunVertexFinderTracks(kTRUE),
  fRunHLTTracking(kFALSE),
  fRunMuonTracking(kFALSE),
  fRunV0Finder(kTRUE),
  fRunCascadeFinder(kTRUE),
  fStopOnError(kFALSE),
  fWriteAlignmentData(kFALSE),
  fWriteESDfriend(kFALSE),
  fWriteAOD(kFALSE),
  fFillTriggerESD(kTRUE),

  fCleanESD(kTRUE),
  fV0DCAmax(3.),
  fV0CsPmin(0.),
  fDmax(50.),
  fZmax(50.),

  fRunLocalReconstruction("ALL"),
  fRunTracking("ALL"),
  fFillESD("ALL"),
  fUseTrackingErrorsForAlignment(""),
  fGAliceFileName(gAliceFilename),
  fInput(""),
  fEquipIdMap(""),
  fFirstEvent(0),
  fLastEvent(-1),
  fNumberOfEventsPerFile(1),
  fCheckPointLevel(0),
  fOptions(),
  fLoadAlignFromCDB(kTRUE),
  fLoadAlignData("ALL"),
  fESDPar(""),
  fUseHLTData(),

  fRunLoader(NULL),
  fRawReader(NULL),
  fParentRawReader(NULL),

  fVertexer(NULL),
  fDiamondProfile(NULL),
  fDiamondProfileTPC(NULL),
  fMeanVertexConstraint(kTRUE),

  fGRPData(NULL),

  fAlignObjArray(NULL),
  fCDBUri(),
  fSpecCDBUri(), 
  fInitCDBCalled(kFALSE),
  fSetRunNumberFromDataCalled(kFALSE),
  fQADetectors("ALL"), 
  fQATasks("ALL"), 
  fRunQA(kTRUE),  
  fRunGlobalQA(kTRUE),
  fInLoopQA(kFALSE),
  fSameQACycle(kFALSE),

  fRunPlaneEff(kFALSE),

  fesd(NULL),
  fhltesd(NULL),
  fesdf(NULL),
  ffile(NULL),
  ftree(NULL),
  fhlttree(NULL),
  ffileOld(NULL),
  ftreeOld(NULL),
  fhlttreeOld(NULL),
  ftVertexer(NULL),
  fIsNewRunLoader(kFALSE),
  fRunAliEVE(kFALSE)
{
// create reconstruction object with default parameters
  
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
    fQADataMaker[iDet] = NULL;
	fQACycles[iDet] = 999999;	
  }
  fQADataMaker[fgkNDetectors]=NULL;  //Global QA
  AliPID pid;
}

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const AliReconstruction& rec) :
  TNamed(rec),

  fUniformField(rec.fUniformField),
  fForcedFieldMap(0x0),
  fRunVertexFinder(rec.fRunVertexFinder),
  fRunVertexFinderTracks(rec.fRunVertexFinderTracks),
  fRunHLTTracking(rec.fRunHLTTracking),
  fRunMuonTracking(rec.fRunMuonTracking),
  fRunV0Finder(rec.fRunV0Finder),
  fRunCascadeFinder(rec.fRunCascadeFinder),
  fStopOnError(rec.fStopOnError),
  fWriteAlignmentData(rec.fWriteAlignmentData),
  fWriteESDfriend(rec.fWriteESDfriend),
  fWriteAOD(rec.fWriteAOD),
  fFillTriggerESD(rec.fFillTriggerESD),

  fCleanESD(rec.fCleanESD),
  fV0DCAmax(rec.fV0DCAmax),
  fV0CsPmin(rec.fV0CsPmin),
  fDmax(rec.fDmax),
  fZmax(rec.fZmax),

  fRunLocalReconstruction(rec.fRunLocalReconstruction),
  fRunTracking(rec.fRunTracking),
  fFillESD(rec.fFillESD),
  fUseTrackingErrorsForAlignment(rec.fUseTrackingErrorsForAlignment),
  fGAliceFileName(rec.fGAliceFileName),
  fInput(rec.fInput),
  fEquipIdMap(rec.fEquipIdMap),
  fFirstEvent(rec.fFirstEvent),
  fLastEvent(rec.fLastEvent),
  fNumberOfEventsPerFile(rec.fNumberOfEventsPerFile),
  fCheckPointLevel(0),
  fOptions(),
  fLoadAlignFromCDB(rec.fLoadAlignFromCDB),
  fLoadAlignData(rec.fLoadAlignData),
  fESDPar(rec.fESDPar),
  fUseHLTData(rec.fUseHLTData),

  fRunLoader(NULL),
  fRawReader(NULL),
  fParentRawReader(NULL),

  fVertexer(NULL),
  fDiamondProfile(NULL),
  fDiamondProfileTPC(NULL),
  fMeanVertexConstraint(rec.fMeanVertexConstraint),

  fGRPData(NULL),

  fAlignObjArray(rec.fAlignObjArray),
  fCDBUri(rec.fCDBUri),
  fSpecCDBUri(), 
  fInitCDBCalled(rec.fInitCDBCalled),
  fSetRunNumberFromDataCalled(rec.fSetRunNumberFromDataCalled),
  fQADetectors(rec.fQADetectors), 
  fQATasks(rec.fQATasks), 
  fRunQA(rec.fRunQA),  
  fRunGlobalQA(rec.fRunGlobalQA),
  fInLoopQA(rec.fInLoopQA),
  fSameQACycle(rec.fSameQACycle),
  fRunPlaneEff(rec.fRunPlaneEff),

  fesd(NULL),
  fhltesd(NULL),
  fesdf(NULL),
  ffile(NULL),
  ftree(NULL),
  fhlttree(NULL),
  ffileOld(NULL),
  ftreeOld(NULL),
  fhlttreeOld(NULL),
  ftVertexer(NULL),
  fIsNewRunLoader(rec.fIsNewRunLoader),
  fRunAliEVE(kFALSE)
{
// copy constructor

  for (Int_t i = 0; i < rec.fOptions.GetEntriesFast(); i++) {
    if (rec.fOptions[i]) fOptions.Add(rec.fOptions[i]->Clone());
  }
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
    fQADataMaker[iDet] = NULL;
	fQACycles[iDet] = rec.fQACycles[iDet];	
  }
  fQADataMaker[fgkNDetectors]=NULL;  //Global QA
  for (Int_t i = 0; i < rec.fSpecCDBUri.GetEntriesFast(); i++) {
    if (rec.fSpecCDBUri[i]) fSpecCDBUri.Add(rec.fSpecCDBUri[i]->Clone());
  }

  fForcedFieldMap=new AliMagWrapCheb(*((AliMagWrapCheb*)rec.fForcedFieldMap));
}

//_____________________________________________________________________________
AliReconstruction& AliReconstruction::operator = (const AliReconstruction& rec)
{
// assignment operator

  this->~AliReconstruction();
  new(this) AliReconstruction(rec);
  return *this;
}

//_____________________________________________________________________________
AliReconstruction::~AliReconstruction()
{
// clean up

  CleanUp();
  fOptions.Delete();
  fSpecCDBUri.Delete();
  delete fForcedFieldMap;

  AliCodeTimer::Instance()->Print();
}

//_____________________________________________________________________________
void AliReconstruction::InitCDB()
{
// activate a default CDB storage
// First check if we have any CDB storage set, because it is used 
// to retrieve the calibration and alignment constants

  if (fInitCDBCalled) return;
  fInitCDBCalled = kTRUE;

  AliCDBManager* man = AliCDBManager::Instance();
  if (man->IsDefaultStorageSet())
  {
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("Default CDB storage has been already set !");
    AliWarning(Form("Ignoring the default storage declared in AliReconstruction: %s",fCDBUri.Data()));
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    fCDBUri = man->GetDefaultStorage()->GetURI();
  }
  else {
    if (fCDBUri.Length() > 0) 
    {
    	AliDebug(2,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    	AliDebug(2, Form("Default CDB storage is set to: %s", fCDBUri.Data()));
    	AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    } else {
    	fCDBUri="local://$ALICE_ROOT";
    	AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    	AliWarning("Default CDB storage not yet set !!!!");
    	AliWarning(Form("Setting it now to: %s", fCDBUri.Data()));
    	AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    		
    }
    man->SetDefaultStorage(fCDBUri);
  }

  // Now activate the detector specific CDB storage locations
  for (Int_t i = 0; i < fSpecCDBUri.GetEntriesFast(); i++) {
    TObject* obj = fSpecCDBUri[i];
    if (!obj) continue;
    AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliDebug(2, Form("Specific CDB storage for %s is set to: %s",obj->GetName(),obj->GetTitle()));
    AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    man->SetSpecificStorage(obj->GetName(), obj->GetTitle());
  }
  
}

//_____________________________________________________________________________
void AliReconstruction::SetDefaultStorage(const char* uri) {
// Store the desired default CDB storage location
// Activate it later within the Run() method

  fCDBUri = uri;

}

//_____________________________________________________________________________
void AliReconstruction::SetSpecificStorage(const char* calibType, const char* uri) {
// Store a detector-specific CDB storage location
// Activate it later within the Run() method

  AliCDBPath aPath(calibType);
  if(!aPath.IsValid()){
	// if calibType is not wildcard but it is a valid detector, add "/*" to make it a valid path
	for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
		if(!strcmp(calibType, fgkDetectorName[iDet])) {
			aPath.SetPath(Form("%s/*", calibType));
			AliInfo(Form("Path for specific storage set to %s", aPath.GetPath().Data()));
			break;
		}
        }
	if(!aPath.IsValid()){
  		AliError(Form("Not a valid path or detector: %s", calibType));
  		return;
	}
  }

//  // check that calibType refers to a "valid" detector name
//  Bool_t isDetector = kFALSE;
//  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
//    TString detName = fgkDetectorName[iDet];
//    if(aPath.GetLevel0() == detName) {
//    	isDetector = kTRUE;
//	break;
//    }
//  }
//
//  if(!isDetector) {
//	AliError(Form("Not a valid detector: %s", aPath.GetLevel0().Data()));
//	return;
//  }

  TObject* obj = fSpecCDBUri.FindObject(aPath.GetPath().Data());
  if (obj) fSpecCDBUri.Remove(obj);
  fSpecCDBUri.Add(new TNamed(aPath.GetPath().Data(), uri));

}

//_____________________________________________________________________________
Bool_t AliReconstruction::SetRunNumberFromData()
{
  // The method is called in Run() in order
  // to set a correct run number.
  // In case of raw data reconstruction the
  // run number is taken from the raw data header

  if (fSetRunNumberFromDataCalled) return kTRUE;
  fSetRunNumberFromDataCalled = kTRUE;
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  if(man->GetRun() > 0) {
  	AliWarning("Run number is taken from raw-event header! Ignoring settings in AliCDBManager!");
  } 
  
  if (!fRunLoader) {
      AliError("No run loader is found !"); 
      return kFALSE;
    }
    // read run number from gAlice
    if(fRunLoader->GetAliRun())
      AliCDBManager::Instance()->SetRun(fRunLoader->GetHeader()->GetRun());
    else {
      if(fRawReader) {
	if(fRawReader->NextEvent()) {
	  AliCDBManager::Instance()->SetRun(fRawReader->GetRunNumber());
	  fRawReader->RewindEvents();
	}
	else {
	  if(man->GetRun() > 0) {
	    AliWarning("No raw events is found ! Using settings in AliCDBManager !");
	    man->Print();  
	    return kTRUE;
	  }
	  else {
	    AliWarning("Neither raw events nor settings in AliCDBManager are found !");
	    return kFALSE;
	  }
	}
      }
      else {
	AliError("Neither gAlice nor RawReader objects are found !");
	return kFALSE;
      }
  }

  man->Print();  
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::SetCDBLock() {
  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  
  AliCDBManager::Instance()->SetLock(1);
}

//_____________________________________________________________________________
Bool_t AliReconstruction::MisalignGeometry(const TString& detectors)
{
  // Read the alignment objects from CDB.
  // Each detector is supposed to have the
  // alignment objects in DET/Align/Data CDB path.
  // All the detector objects are then collected,
  // sorted by geometry level (starting from ALIC) and
  // then applied to the TGeo geometry.
  // Finally an overlaps check is performed.

  // Load alignment data from CDB and fill fAlignObjArray 
  if(fLoadAlignFromCDB){
  	
    TString detStr = detectors;
    TString loadAlObjsListOfDets = "";
    
    for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
      if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
      loadAlObjsListOfDets += fgkDetectorName[iDet];
      loadAlObjsListOfDets += " ";
    } // end loop over detectors
    loadAlObjsListOfDets.Prepend("GRP "); //add alignment objects for non-sensitive modules
    AliGeomManager::ApplyAlignObjsFromCDB(loadAlObjsListOfDets.Data());
  }else{
    // Check if the array with alignment objects was
    // provided by the user. If yes, apply the objects
    // to the present TGeo geometry
    if (fAlignObjArray) {
      if (gGeoManager && gGeoManager->IsClosed()) {
	if (AliGeomManager::ApplyAlignObjsToGeom(*fAlignObjArray) == kFALSE) {
	  AliError("The misalignment of one or more volumes failed!"
		   "Compare the list of simulated detectors and the list of detector alignment data!");
	  return kFALSE;
	}
      }
      else {
	AliError("Can't apply the misalignment! gGeoManager doesn't exist or it is still opened!");
	return kFALSE;
      }
    }
  }
  
  delete fAlignObjArray; fAlignObjArray=0;

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::SetGAliceFile(const char* fileName)
{
// set the name of the galice file

  fGAliceFileName = fileName;
}

//_____________________________________________________________________________
void AliReconstruction::SetInput(const char* input) 
{
  // In case the input string starts with 'mem://', we run in an online mode
  // and AliRawReaderDateOnline object is created. In all other cases a raw-data
  // file is assumed. One can give as an input:
  // mem://: - events taken from DAQ monitoring libs online
  //  or
  // mem://<filename> - emulation of the above mode (via DATE monitoring libs)
  fInput = input;
}

//_____________________________________________________________________________
void AliReconstruction::SetOption(const char* detector, const char* option)
{
// set options for the reconstruction of a detector

  TObject* obj = fOptions.FindObject(detector);
  if (obj) fOptions.Remove(obj);
  fOptions.Add(new TNamed(detector, option));
}

//_____________________________________________________________________________
Bool_t AliReconstruction::SetFieldMap(Float_t l3Current, Float_t diCurrent, Float_t factor, const char *path) {
  //------------------------------------------------
  // The magnetic field map, defined externally...
  // L3 current 30000 A  -> 0.5 T
  // L3 current 12000 A  -> 0.2 T
  // dipole current 6000 A
  // The polarities must be the same
  //------------------------------------------------
  const Float_t l3NominalCurrent1=30000.; // (A)
  const Float_t l3NominalCurrent2=12000.; // (A)
  const Float_t diNominalCurrent =6000. ; // (A)

  const Float_t tolerance=0.03; // relative current tolerance
  const Float_t zero=77.;       // "zero" current (A)

  Int_t map=0;
  Bool_t dipoleON=kFALSE;

  TString s=(factor < 0) ? "L3: -" : "L3: +";

  if (TMath::Abs(l3Current-l3NominalCurrent1)/l3NominalCurrent1 < tolerance) {
    map=AliMagWrapCheb::k5kG;
    s+="0.5 T;  ";
  } else
  if (TMath::Abs(l3Current-l3NominalCurrent2)/l3NominalCurrent2 < tolerance) {
    map=AliMagWrapCheb::k2kG;
    s+="0.2 T;  ";
  } else
  if (TMath::Abs(l3Current) < zero) {
    map=AliMagWrapCheb::k2kG;
    s+="0.0 T;  ";
    factor=0.;                  // in fact, this is a global factor...
  } else {
    AliError("Wrong L3 current !");
    return kFALSE;
  }

  if (TMath::Abs(diCurrent-diNominalCurrent)/diNominalCurrent < tolerance) {
    // 3% current tolerance...
    dipoleON=kTRUE;
    s+="Dipole ON";
  } else
  if (TMath::Abs(diCurrent) < zero) { // some small current..
    dipoleON=kFALSE;
    s+="Dipole OFF";
  } else {
    AliError("Wrong dipole current !");
    return kFALSE;
  }

  delete fForcedFieldMap;
  fForcedFieldMap=
    new AliMagWrapCheb("B field map  ",s,2,factor,10.,map,dipoleON,path);

  fForcedFieldMap->Print();

  AliTracker::SetFieldMap(fForcedFieldMap,fUniformField);    

  return kTRUE;
}


Bool_t AliReconstruction::InitGRP() {
  //------------------------------------
  // Initialization of the GRP entry 
  //------------------------------------
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  if (entry) fGRPData = dynamic_cast<TMap*>(entry->GetObject());  

  if (!fGRPData) {
     AliError("No GRP entry found in OCDB!");
     return kFALSE;
  }


  //*** Dealing with the magnetic field map
  if (AliTracker::GetFieldMap()) {
    AliInfo("Running with the externally set B field !");
  } else {
    // Construct the field map out of the information retrieved from GRP.

    Bool_t ok = kTRUE;

    // L3
    TObjString *l3Current=
       dynamic_cast<TObjString*>(fGRPData->GetValue("fL3Current"));
    if (!l3Current) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    TObjString *l3Polarity=
       dynamic_cast<TObjString*>(fGRPData->GetValue("fL3Polarity"));
    if (!l3Polarity) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    TObjString *diCurrent=
       dynamic_cast<TObjString*>(fGRPData->GetValue("fDipoleCurrent"));
    if (!diCurrent) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }
    TObjString *diPolarity=
       dynamic_cast<TObjString*>(fGRPData->GetValue("fDipolePolarity"));
    if (!diPolarity) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    if (ok) { 
       Float_t l3Cur=TMath::Abs(atof(l3Current->GetName()));
       Float_t diCur=TMath::Abs(atof(diCurrent->GetName()));
       Float_t l3Pol=atof(l3Polarity->GetName());
       Float_t factor=1.;
       if (l3Pol != 0.) factor=-1.;
    

      if (!SetFieldMap(l3Cur, diCur, factor)) {
         AliFatal("Failed to creat a B field map ! Exiting...");
      }
      AliInfo("Running with the B field constructed out of GRP !");
    }
    else {
      AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
    }

  }


  //*** Get the diamond profile from OCDB
  entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertex");
  if (entry) {
     fDiamondProfile = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  } else {
     AliError("No diamond profile found in OCDB!");
  }

  entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexTPC");
  if (entry) {
     fDiamondProfileTPC = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  } else {
     AliError("No diamond profile found in OCDB!");
  }

  return kTRUE;
} 

//_____________________________________________________________________________
Bool_t AliReconstruction::Run(const char* input)
{
  // Run Run Run
  AliCodeTimerAuto("");

  if (!InitRun(input)) return kFALSE;
  //******* The loop over events
  Int_t iEvent = 0;
  while ((iEvent < fRunLoader->GetNumberOfEvents()) ||
	 (fRawReader && fRawReader->NextEvent())) {
    if (!RunEvent(iEvent)) return kFALSE;
    iEvent++;
  }

  if (!FinishRun()) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::InitRun(const char* input)
{
  // Initialize all the stuff before
  // going into the event loop
  // If the second argument is given, the first one is ignored and
  // the reconstruction works in an online mode
  AliCodeTimerAuto("");

  // Overwrite the previous setting
  if (input) fInput = input;

  // set the input in case of raw data
  fRawReader = AliRawReader::Create(fInput.Data());
  if (!fRawReader)
    AliInfo("Reconstruction will run over digits");

  if (!fEquipIdMap.IsNull() && fRawReader)
    fRawReader->LoadEquipmentIdsMap(fEquipIdMap);

  if (!fUseHLTData.IsNull()) {
    // create the RawReaderHLT which performs redirection of HLT input data for
    // the specified detectors
    AliRawReader* pRawReader=AliRawHLTManager::CreateRawReaderHLT(fRawReader, fUseHLTData.Data());
    if (pRawReader) {
      fParentRawReader=fRawReader;
      fRawReader=pRawReader;
    } else {
      AliError(Form("can not create Raw Reader for HLT input %s", fUseHLTData.Data()));
    }
  }

   AliSysInfo::AddStamp("Start");
  // get the run loader
  if (!InitRunLoader()) return kFALSE;
   AliSysInfo::AddStamp("LoadLoader");

  // Initialize the CDB storage
  InitCDB();
  
  AliSysInfo::AddStamp("LoadCDB");

  // Set run number in CDBManager (if it is not already set by the user)
  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  
  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  SetCDBLock();
  
  // Import ideal TGeo geometry and apply misalignment
  if (!gGeoManager) {
    TString geom(gSystem->DirName(fGAliceFileName));
    geom += "/geometry.root";
    AliGeomManager::LoadGeometry(geom.Data());

    TString detsToCheck=fRunLocalReconstruction;
    if(!AliGeomManager::CheckSymNamesLUT(detsToCheck.Data()))
	  AliFatalClass("Current loaded geometry differs in the definition of symbolic names!");
    if (!gGeoManager) if (fStopOnError) return kFALSE;
  }

  if (!MisalignGeometry(fLoadAlignData)) if (fStopOnError) return kFALSE;
   AliSysInfo::AddStamp("LoadGeom");


  if (!InitGRP()) return kFALSE;


  ftVertexer = new AliVertexerTracks(AliTracker::GetBz());
  if(fDiamondProfile && fMeanVertexConstraint) ftVertexer->SetVtxStart(fDiamondProfile);

  // get vertexer
  if (fRunVertexFinder && !CreateVertexer()) {
    if (fStopOnError) {
      CleanUp(); 
      return kFALSE;
    }
  }
   AliSysInfo::AddStamp("Vertexer");

  // get trackers
  if (!fRunTracking.IsNull() && !CreateTrackers(fRunTracking)) {
    if (fStopOnError) {
      CleanUp(); 
      return kFALSE;
    }      
  }
   AliSysInfo::AddStamp("LoadTrackers");

  // get the possibly already existing ESD file and tree
  fesd = new AliESDEvent(); fhltesd = new AliESDEvent();
  if (!gSystem->AccessPathName("AliESDs.root")){
    gSystem->CopyFile("AliESDs.root", "AliESDs.old.root", kTRUE);
    ffileOld = TFile::Open("AliESDs.old.root");
    if (ffileOld && ffileOld->IsOpen()) {
      ftreeOld = (TTree*) ffileOld->Get("esdTree");
      if (ftreeOld)fesd->ReadFromTree(ftreeOld);
      fhlttreeOld = (TTree*) ffileOld->Get("HLTesdTree");
      if (fhlttreeOld)	fhltesd->ReadFromTree(fhlttreeOld);
    }
  }

  // create the ESD output file and tree
  ffile = TFile::Open("AliESDs.root", "RECREATE");
  ffile->SetCompressionLevel(2);
  if (!ffile->IsOpen()) {
    AliError("opening AliESDs.root failed");
    if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}    
  }

  ftree = new TTree("esdTree", "Tree with ESD objects");
  fesd = new AliESDEvent();
  fesd->CreateStdContent();
  fesd->WriteToTree(ftree);

  fhlttree = new TTree("HLTesdTree", "Tree with HLT ESD objects");
  fhltesd = new AliESDEvent();
  fhltesd->CreateStdContent();
  fhltesd->WriteToTree(fhlttree);


  if (fWriteESDfriend) {
    fesdf = new AliESDfriend();
    TBranch *br=ftree->Branch("ESDfriend.","AliESDfriend", &fesdf);
    br->SetFile("AliESDfriends.root");
    fesd->AddObject(fesdf);
  }



  if (fRawReader) fRawReader->RewindEvents();

  ProcInfo_t ProcInfo;
  gSystem->GetProcInfo(&ProcInfo);
  AliInfo(Form("Current memory usage %d %d", ProcInfo.fMemResident, ProcInfo.fMemVirtual));
  
  //QA
  if (fRunQA && fRawReader && fQATasks.Contains(Form("%d", AliQA::kRAWS))) { 
		AliQADataMakerSteer qas ; 
		qas.Run(fRunLocalReconstruction, fRawReader) ; 
		fSameQACycle = kTRUE ; 
  }
  //Initialize the QA and start of cycle for out-of-cycle QA
  if (fRunQA) {
	  TString detStr(fQADetectors) ; 
      for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
         if (!IsSelected(fgkDetectorName[iDet], detStr)) 
			 continue;
         AliQADataMakerRec *qadm = GetQADataMaker(iDet);  
         if (!qadm) 
			 continue;
         AliInfo(Form("Initializing the QA data maker for %s", 
                fgkDetectorName[iDet]));
		 if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) 
			 qadm->Init(AliQA::kRECPOINTS, AliCDBManager::Instance()->GetRun());
		 if (fQATasks.Contains(Form("%d", AliQA::kESDS)))  
			 qadm->Init(AliQA::kESDS, AliCDBManager::Instance()->GetRun());
         if (!fInLoopQA) {
			 if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) {
				 qadm->StartOfCycle(AliQA::kRECPOINTS, fSameQACycle);
				 fSameQACycle = kTRUE;
			 }
			 if (fQATasks.Contains(Form("%d", AliQA::kESDS))) {  
				 qadm->StartOfCycle(AliQA::kESDS, fSameQACycle);
				 fSameQACycle = kTRUE;
			 }
		 }
      }
	  if (fRunGlobalQA) {
		  AliQADataMakerRec *qadm = GetQADataMaker(AliQA::kGLOBAL);
		  AliInfo(Form("Initializing the global QA data maker"));
		  if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) {
			  TObjArray *arr=
				qadm->Init(AliQA::kRECPOINTS, AliCDBManager::Instance()->GetRun());
			  AliTracker::SetResidualsArray(arr);
		  }
		  if (fQATasks.Contains(Form("%d", AliQA::kESDS))) {
			  qadm->Init(AliQA::kESDS, AliCDBManager::Instance()->GetRun());
		  }
		  if (!fInLoopQA) {
			  fSameQACycle = kFALSE;
			  if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) { 				  
				  qadm->StartOfCycle(AliQA::kRECPOINTS, fSameQACycle);
				  fSameQACycle = kTRUE;
			  }
			  if (fQATasks.Contains(Form("%d", AliQA::kESDS))) { 
				  qadm->StartOfCycle(AliQA::kESDS, fSameQACycle);
				  fSameQACycle = kTRUE;	
			  }
		  }
	  }
  }

  //Initialize the Plane Efficiency framework
  if (fRunPlaneEff && !InitPlaneEff()) {
    if(fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
  }

  if (strcmp(gProgName,"alieve") == 0)
    fRunAliEVE = InitAliEVE();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunEvent(Int_t iEvent)
{
  // run the reconstruction over a single event
  // The event loop is steered in Run method

  AliCodeTimerAuto("");

  if (iEvent >= fRunLoader->GetNumberOfEvents()) {
    fRunLoader->SetEventNumber(iEvent);
    fRunLoader->GetHeader()->Reset(fRawReader->GetRunNumber(), 
				   iEvent, iEvent);
    //??      fRunLoader->MakeTree("H");
    fRunLoader->TreeE()->Fill();
  }

  if ((iEvent < fFirstEvent) || ((fLastEvent >= 0) && (iEvent > fLastEvent))) {
    // copy old ESD to the new one
    if (ftreeOld) {
      fesd->ReadFromTree(ftreeOld);
      ftreeOld->GetEntry(iEvent);
      ftree->Fill();
    }
    if (fhlttreeOld) {
      fhltesd->ReadFromTree(fhlttreeOld);
      fhlttreeOld->GetEntry(iEvent);
      fhlttree->Fill();
    }
    return kTRUE;
  }

  AliInfo(Form("processing event %d", iEvent));

    //Start of cycle for the in-loop QA
    if (fInLoopQA) {
       if (fRunQA) {
		  fSameQACycle = kFALSE ;
          TString detStr(fQADetectors); 
          for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
             if (!IsSelected(fgkDetectorName[iDet], detStr)) 
				 continue;
             AliQADataMakerRec *qadm = GetQADataMaker(iDet);  
             if (!qadm) 
				 continue;
			  if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) {
				  qadm->StartOfCycle(AliQA::kRECPOINTS, fSameQACycle);
				  fSameQACycle = kTRUE;
			  }
			  if (fQATasks.Contains(Form("%d", AliQA::kESDS))) {
				  qadm->StartOfCycle(AliQA::kESDS, fSameQACycle) ;
				  fSameQACycle = kTRUE;
			  }
          }
		   if (fRunGlobalQA) {
			   fSameQACycle = kFALSE;
			   AliQADataMakerRec *qadm = GetQADataMaker(AliQA::kGLOBAL);
			   if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) {
				   qadm->StartOfCycle(AliQA::kRECPOINTS, fSameQACycle);
				   fSameQACycle = kTRUE;
			   }
			   if (fQATasks.Contains(Form("%d", AliQA::kESDS))) {
				   qadm->StartOfCycle(AliQA::kESDS, fSameQACycle);
				   fSameQACycle = kTRUE;
			   }
		   }		   
	   }
    }

    fRunLoader->GetEvent(iEvent);

    char aFileName[256];
    sprintf(aFileName, "ESD_%d.%d_final.root", 
	    fRunLoader->GetHeader()->GetRun(), 
	    fRunLoader->GetHeader()->GetEventNrInRun());
    if (!gSystem->AccessPathName(aFileName)) return kTRUE;

    // local single event reconstruction
    if (!fRunLocalReconstruction.IsNull()) {
      TString detectors=fRunLocalReconstruction;
      // run HLT event reconstruction first
      // ;-( IsSelected changes the string
      if (IsSelected("HLT", detectors) &&
	  !RunLocalEventReconstruction("HLT")) {
	if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
      }
      detectors=fRunLocalReconstruction;
      detectors.ReplaceAll("HLT", "");
      if (!RunLocalEventReconstruction(detectors)) {
	if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
      }
    }

    fesd->SetRunNumber(fRunLoader->GetHeader()->GetRun());
    fhltesd->SetRunNumber(fRunLoader->GetHeader()->GetRun());
    fesd->SetEventNumberInFile(fRunLoader->GetHeader()->GetEventNrInRun());
    fhltesd->SetEventNumberInFile(fRunLoader->GetHeader()->GetEventNrInRun());
    
    // Set magnetic field from the tracker
    fesd->SetMagneticField(AliTracker::GetBz());
    fhltesd->SetMagneticField(AliTracker::GetBz());

    
    
    // Fill raw-data error log into the ESD
    if (fRawReader) FillRawDataErrorLog(iEvent,fesd);

    // vertex finder
    if (fRunVertexFinder) {
      if (!ReadESD(fesd, "vertex")) {
	if (!RunVertexFinder(fesd)) {
	  if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
	}
	if (fCheckPointLevel > 0) WriteESD(fesd, "vertex");
      }
    }

    // Muon tracking
    if (!fRunTracking.IsNull()) {
      if (fRunMuonTracking) {
	if (!RunMuonTracking(fesd)) {
	  if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
	}
      }
    }

    // barrel tracking
    if (!fRunTracking.IsNull()) {
      if (!ReadESD(fesd, "tracking")) {
	if (!RunTracking(fesd)) {
	  if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
	}
	if (fCheckPointLevel > 0) WriteESD(fesd, "tracking");
      }
    }

    // fill ESD
    if (!fFillESD.IsNull()) {
      TString detectors=fFillESD;
      // run HLT first and on hltesd
      // ;-( IsSelected changes the string
      if (IsSelected("HLT", detectors) &&
	  !FillESD(fhltesd, "HLT")) {
	if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
      }
      detectors=fFillESD;
      // Temporary fix to avoid problems with HLT that overwrites the offline ESDs
      if (detectors.Contains("ALL")) {
	detectors="";
	for (Int_t idet=0; idet<fgkNDetectors; ++idet){
	  detectors += fgkDetectorName[idet];
	  detectors += " ";
	}
      }
      detectors.ReplaceAll("HLT", "");
      if (!FillESD(fesd, detectors)) {
	if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
      }
    }
  
    // fill Event header information from the RawEventHeader
    if (fRawReader){FillRawEventHeaderESD(fesd);}

    // combined PID
    AliESDpid::MakePID(fesd);
    if (fCheckPointLevel > 1) WriteESD(fesd, "PID");

    if (fFillTriggerESD) {
      if (!ReadESD(fesd, "trigger")) {
	if (!FillTriggerESD(fesd)) {
	  if (fStopOnError) {CleanUp(ffile, ffileOld); return kFALSE;}
	}
	if (fCheckPointLevel > 1) WriteESD(fesd, "trigger");
      }
    }

    ffile->cd();

    //
    // Propagate track to the beam pipe  (if not already done by ITS)
    //
    const Int_t ntracks = fesd->GetNumberOfTracks();
    const Double_t kBz = fesd->GetMagneticField();
    const Double_t kRadius  = 2.8; //something less than the beam pipe radius

    TObjArray trkArray;
    UShort_t *selectedIdx=new UShort_t[ntracks];

    for (Int_t itrack=0; itrack<ntracks; itrack++){
      const Double_t kMaxStep = 5;   //max step over the material
      Bool_t ok;

      AliESDtrack *track = fesd->GetTrack(itrack);
      if (!track) continue;

      AliExternalTrackParam *tpcTrack =
           (AliExternalTrackParam *)track->GetTPCInnerParam();
      ok = kFALSE;
      if (tpcTrack)
	ok = AliTracker::
	  PropagateTrackTo(tpcTrack,kRadius,track->GetMass(),kMaxStep,kTRUE);

      if (ok) {
	Int_t n=trkArray.GetEntriesFast();
        selectedIdx[n]=track->GetID();
        trkArray.AddLast(tpcTrack);
      }

      //Tracks refitted by ITS should already be at the SPD vertex
      if (track->IsOn(AliESDtrack::kITSrefit)) continue;

      AliTracker::
         PropagateTrackTo(track,kRadius,track->GetMass(),kMaxStep,kTRUE);
      track->RelateToVertex(fesd->GetPrimaryVertexSPD(), kBz, kVeryBig);

    }

    //
    // Improve the reconstructed primary vertex position using the tracks
    //
    TObject *obj = fOptions.FindObject("ITS");
    if (obj) {
      TString optITS = obj->GetTitle();
      if (optITS.Contains("cosmics") || optITS.Contains("COSMICS")) 
	fRunVertexFinderTracks=kFALSE;
    }
    if (fRunVertexFinderTracks) {
       // TPC + ITS primary vertex
       ftVertexer->SetITSrefitRequired();
       if(fDiamondProfile && fMeanVertexConstraint) {
	 ftVertexer->SetVtxStart(fDiamondProfile);
       } else {
	 ftVertexer->SetConstraintOff();
       }
       AliESDVertex *pvtx=ftVertexer->FindPrimaryVertex(fesd);
       if (pvtx) {
          if (pvtx->GetStatus()) {
             fesd->SetPrimaryVertex(pvtx);
             for (Int_t i=0; i<ntracks; i++) {
	         AliESDtrack *t = fesd->GetTrack(i);
                 t->RelateToVertex(pvtx, kBz, kVeryBig);
             } 
          }
       }

       // TPC-only primary vertex
       ftVertexer->SetITSrefitNotRequired();
       if(fDiamondProfileTPC && fMeanVertexConstraint) {
	 ftVertexer->SetVtxStart(fDiamondProfileTPC);
       } else {
	 ftVertexer->SetConstraintOff();
       }
       pvtx=ftVertexer->FindPrimaryVertex(&trkArray,selectedIdx);
       if (pvtx) {
          if (pvtx->GetStatus()) {
             fesd->SetPrimaryVertexTPC(pvtx);
             for (Int_t i=0; i<ntracks; i++) {
	         AliESDtrack *t = fesd->GetTrack(i);
                 t->RelateToVertexTPC(pvtx, kBz, kVeryBig);
             } 
          }
       }

    }
    delete[] selectedIdx;

    if(fDiamondProfile) fesd->SetDiamond(fDiamondProfile);
    

    if (fRunV0Finder) {
       // V0 finding
       AliV0vertexer vtxer;
       vtxer.Tracks2V0vertices(fesd);

       if (fRunCascadeFinder) {
          // Cascade finding
          AliCascadeVertexer cvtxer;
          cvtxer.V0sTracks2CascadeVertices(fesd);
       }
    }
 
    // write ESD
    if (fCleanESD) CleanESD(fesd);

	if (fRunQA) {
		if (fRunGlobalQA) {
			AliQADataMakerRec *qadm = GetQADataMaker(AliQA::kGLOBAL);
			if (qadm && fQATasks.Contains(Form("%d", AliQA::kESDS)))
				qadm->Exec(AliQA::kESDS, fesd);
		}
	}

    if (fWriteESDfriend) {
      fesdf->~AliESDfriend();
      new (fesdf) AliESDfriend(); // Reset...
      fesd->GetESDfriend(fesdf);
    }
    ftree->Fill();

    // write HLT ESD
    fhlttree->Fill();

    // call AliEVE
    if (fRunAliEVE) RunAliEVE();

    if (fCheckPointLevel > 0)  WriteESD(fesd, "final"); 
    fesd->Reset();
    fhltesd->Reset();
    if (fWriteESDfriend) {
      fesdf->~AliESDfriend();
      new (fesdf) AliESDfriend(); // Reset...
    }
 
    ProcInfo_t ProcInfo;
    gSystem->GetProcInfo(&ProcInfo);
    AliInfo(Form("Event %d -> Current memory usage %d %d",iEvent, ProcInfo.fMemResident, ProcInfo.fMemVirtual));
  

  // End of cycle for the in-loop  
     if (fInLoopQA) {
        if (fRunQA) {
           RunQA(fesd);
           for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
			   if (!IsSelected(fgkDetectorName[iDet], fQADetectors)) 
				   continue;
			   AliQADataMakerRec * qadm = GetQADataMaker(iDet);
			   if (!qadm)
				   continue;
			   if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) 
				   qadm->EndOfCycle(AliQA::kRECPOINTS);
			   if (fQATasks.Contains(Form("%d", AliQA::kESDS))) 
				   qadm->EndOfCycle(AliQA::kESDS);
			   qadm->Finish();
		   }
        }
        if (fRunGlobalQA) {
           AliQADataMakerRec *qadm = GetQADataMaker(AliQA::kGLOBAL);
           if (qadm) {
			   if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) 
				   qadm->EndOfCycle(AliQA::kRECPOINTS);
			   if (fQATasks.Contains(Form("%d", AliQA::kESDS))) 
				   qadm->EndOfCycle(AliQA::kESDS);
			   qadm->Finish();
		   }
        }
     }

     return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FinishRun()
{
  // Finalize the run
  // Called after the exit
  // from the event loop
  AliCodeTimerAuto("");

  if (fIsNewRunLoader) { // galice.root didn't exist
    fRunLoader->WriteHeader("OVERWRITE");
    fRunLoader->CdGAFile();
    fRunLoader->Write(0, TObject::kOverwrite);
  }

  ftree->GetUserInfo()->Add(fesd);
  fhlttree->GetUserInfo()->Add(fhltesd);
  
  const TMap *cdbMap = AliCDBManager::Instance()->GetStorageMap();	 
  const TList *cdbList = AliCDBManager::Instance()->GetRetrievedIds();	 
	 	 
   TMap *cdbMapCopy = new TMap(cdbMap->GetEntries());	 
   cdbMapCopy->SetOwner(1);	 
   cdbMapCopy->SetName("cdbMap");	 
   TIter iter(cdbMap->GetTable());	 
 	 
   TPair* pair = 0;	 
   while((pair = dynamic_cast<TPair*> (iter.Next()))){	 
         TObjString* keyStr = dynamic_cast<TObjString*> (pair->Key());	 
         TObjString* valStr = dynamic_cast<TObjString*> (pair->Value());	 
         cdbMapCopy->Add(new TObjString(keyStr->GetName()), new TObjString(valStr->GetName()));	 
   }	 
 	 
   TList *cdbListCopy = new TList();	 
   cdbListCopy->SetOwner(1);	 
   cdbListCopy->SetName("cdbList");	 
 	 
   TIter iter2(cdbList);	 
 	 
   AliCDBId* id=0;	 
   while((id = dynamic_cast<AliCDBId*> (iter2.Next()))){	 
         cdbListCopy->Add(new TObjString(id->ToString().Data()));	 
   }	 
 	 
   ftree->GetUserInfo()->Add(cdbMapCopy);	 
   ftree->GetUserInfo()->Add(cdbListCopy);


  if(fESDPar.Contains("ESD.par")){
    AliInfo("Attaching ESD.par to Tree");
    TNamed *fn = CopyFileToTNamed(fESDPar.Data(),"ESD.par");
    ftree->GetUserInfo()->Add(fn);
  }


  ffile->cd();

  if (fWriteESDfriend)
    ftree->SetBranchStatus("ESDfriend*",0);
  // we want to have only one tree version number
  ftree->Write(ftree->GetName(),TObject::kOverwrite);
  fhlttree->Write();

// Finish with Plane Efficiency evaluation: before of CleanUp !!!
  if (fRunPlaneEff && !FinishPlaneEff()) {
   AliWarning("Finish PlaneEff evaluation failed");
  }

  gROOT->cd();
  CleanUp(ffile, ffileOld);
    
  if (fWriteAOD) {
    AliWarning("AOD creation not supported anymore during reconstruction. See ANALYSIS/AliAnalysisTaskESDfilter.cxx instead.");
  }

  // Create tags for the events in the ESD tree (the ESD tree is always present)
  // In case of empty events the tags will contain dummy values
  AliESDTagCreator *esdtagCreator = new AliESDTagCreator();
  esdtagCreator->CreateESDTags(fFirstEvent,fLastEvent,fGRPData);
  if (fWriteAOD) {
    AliWarning("AOD tag creation not supported anymore during reconstruction.");
  }

  //Finish QA and end of cycle for out-of-loop QA
  if (!fInLoopQA) {
	  if (fRunQA) {
		  AliQADataMakerSteer qas;
		  if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) 
			  qas.Run(fRunLocalReconstruction.Data(), AliQA::kRECPOINTS, fSameQACycle);
		  //qas.Reset() ;
		  if (fQATasks.Contains(Form("%d", AliQA::kESDS))) 
			  qas.Run(fRunLocalReconstruction.Data(), AliQA::kESDS, fSameQACycle);
		  if (fRunGlobalQA) {
			 AliQADataMakerRec *qadm = GetQADataMaker(AliQA::kGLOBAL);
			  if (qadm) {
				  if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) 
					  qadm->EndOfCycle(AliQA::kRECPOINTS);
				  if (fQATasks.Contains(Form("%d", AliQA::kESDS))) 
					  qadm->EndOfCycle(AliQA::kESDS);
				  qadm->Finish();
			  }
		  }
	  }
  }
  
  // Cleanup of CDB manager: cache and active storages!
  AliCDBManager::Instance()->ClearCache();
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::RunLocalReconstruction(const TString& /*detectors*/)
{
// run the local reconstruction
  static Int_t eventNr=0;
  AliCodeTimerAuto("")

 //  AliCDBManager* man = AliCDBManager::Instance();
//   Bool_t origCache = man->GetCacheFlag();

//   TString detStr = detectors;
//   for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
//     if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
//     AliReconstructor* reconstructor = GetReconstructor(iDet);
//     if (!reconstructor) continue;
//     if (reconstructor->HasLocalReconstruction()) continue;

//     AliCodeTimerStart(Form("running reconstruction for %s", fgkDetectorName[iDet]));
//     AliInfo(Form("running reconstruction for %s", fgkDetectorName[iDet]));
    
//     AliCodeTimerStart(Form("Loading calibration data from OCDB for %s", fgkDetectorName[iDet]));                          
//     AliInfo(Form("Loading calibration data from OCDB for %s", fgkDetectorName[iDet]));

//     man->SetCacheFlag(kTRUE);
//     TString calibPath = Form("%s/Calib/*", fgkDetectorName[iDet]);
//     man->GetAll(calibPath); // entries are cached!

//     AliCodeTimerStop(Form("Loading calibration data from OCDB for %s", fgkDetectorName[iDet]));
     
//     if (fRawReader) {
//       fRawReader->RewindEvents();
//       reconstructor->Reconstruct(fRunLoader, fRawReader);
//     } else {
//       reconstructor->Reconstruct(fRunLoader);
//     }
     
//      AliCodeTimerStop(Form("running reconstruction for %s", fgkDetectorName[iDet]));
    // AliSysInfo::AddStamp(Form("LRec%s_%d",fgkDetectorName[iDet],eventNr));

//     // unload calibration data
//     man->UnloadFromCache(calibPath);
//     //man->ClearCache();
//   }

//   man->SetCacheFlag(origCache);

//   if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
//     AliError(Form("the following detectors were not found: %s",
//                   detStr.Data()));
//     if (fStopOnError) return kFALSE;
//   }

	  eventNr++;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunLocalEventReconstruction(const TString& detectors)
{
// run the local reconstruction

  static Int_t eventNr=0;
  AliCodeTimerAuto("")

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
    AliLoader* loader = fLoader[iDet];
    // Matthias April 2008: temporary fix to run HLT reconstruction
    // although the HLT loader is missing
    if (strcmp(fgkDetectorName[iDet], "HLT")==0) {
      if (fRawReader) {
	reconstructor->Reconstruct(fRawReader, NULL);
      } else {
	TTree* dummy=NULL;
	reconstructor->Reconstruct(dummy, NULL);
      }
      continue;
    }
    if (!loader) {
      AliWarning(Form("No loader is defined for %s!",fgkDetectorName[iDet]));
      continue;
    }
    // conversion of digits
    if (fRawReader && reconstructor->HasDigitConversion()) {
      AliInfo(Form("converting raw data digits into root objects for %s", 
		   fgkDetectorName[iDet]));
      AliCodeTimerAuto(Form("converting raw data digits into root objects for %s", 
                            fgkDetectorName[iDet]));
      loader->LoadDigits("update");
      loader->CleanDigits();
      loader->MakeDigitsContainer();
      TTree* digitsTree = loader->TreeD();
      reconstructor->ConvertDigits(fRawReader, digitsTree);
      loader->WriteDigits("OVERWRITE");
      loader->UnloadDigits();
    }
    // local reconstruction
    AliInfo(Form("running reconstruction for %s", fgkDetectorName[iDet]));
    AliCodeTimerAuto(Form("running reconstruction for %s", fgkDetectorName[iDet]));
    loader->LoadRecPoints("update");
    loader->CleanRecPoints();
    loader->MakeRecPointsContainer();
    TTree* clustersTree = loader->TreeR();
    if (fRawReader && !reconstructor->HasDigitConversion()) {
      reconstructor->Reconstruct(fRawReader, clustersTree);
    } else {
      loader->LoadDigits("read");
      TTree* digitsTree = loader->TreeD();
      if (!digitsTree) {
	AliError(Form("Can't get the %s digits tree", fgkDetectorName[iDet]));
	if (fStopOnError) return kFALSE;
      } else {
	reconstructor->Reconstruct(digitsTree, clustersTree);
      }
      loader->UnloadDigits();
    }

    // In-loop QA for local reconstrucion 
    if (fRunQA && fInLoopQA) {
       AliQADataMakerRec * qadm = GetQADataMaker(iDet);
       if (qadm) {
	  //AliCodeTimerStart
	  //(Form("Running QA data maker for %s", fgkDetectorName[iDet]));
	  //AliInfo
          //(Form("Running QA data maker for %s", fgkDetectorName[iDet]));

		   if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) 
			   qadm->Exec(AliQA::kRECPOINTS, clustersTree) ;
	  //AliCodeTimerStop
          //(Form("Running QA data maker for %s", fgkDetectorName[iDet]));
       }
    }

    loader->WriteRecPoints("OVERWRITE");
    loader->UnloadRecPoints();
    AliSysInfo::AddStamp(Form("LRec%s_%d",fgkDetectorName[iDet],eventNr), iDet,1,eventNr);
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s",
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }
  eventNr++;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunVertexFinder(AliESDEvent*& esd)
{
// run the barrel tracking

  AliCodeTimerAuto("")

  AliESDVertex* vertex = NULL;
  Double_t vtxPos[3] = {0, 0, 0};
  Double_t vtxErr[3] = {0.07, 0.07, 0.1};
  TArrayF mcVertex(3); 
  if (fRunLoader->GetHeader() && fRunLoader->GetHeader()->GenEventHeader()) {
    fRunLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
    for (Int_t i = 0; i < 3; i++) vtxPos[i] = mcVertex[i];
  }

  if (fVertexer) {
    AliInfo("running the ITS vertex finder");
    if (fLoader[0]) {
      fLoader[0]->LoadRecPoints();
      TTree* cltree = fLoader[0]->TreeR();
      if (cltree) {
	if(fDiamondProfile) fVertexer->SetVtxStart(fDiamondProfile);
	vertex = fVertexer->FindVertexForCurrentEvent(cltree);
      }
      else {
	AliError("Can't get the ITS cluster tree");
      }
      fLoader[0]->UnloadRecPoints();
    }
    else {
      AliError("Can't get the ITS loader");
    }
    if(!vertex){
      AliWarning("Vertex not found");
      vertex = new AliESDVertex();
      vertex->SetName("default");
    }
    else {
      vertex->SetName("reconstructed");
    }

  } else {
    AliInfo("getting the primary vertex from MC");
    vertex = new AliESDVertex(vtxPos, vtxErr);
  }

  if (vertex) {
    vertex->GetXYZ(vtxPos);
    vertex->GetSigmaXYZ(vtxErr);
  } else {
    AliWarning("no vertex reconstructed");
    vertex = new AliESDVertex(vtxPos, vtxErr);
  }
  esd->SetPrimaryVertexSPD(vertex);
  // if SPD multiplicity has been determined, it is stored in the ESD
  AliMultiplicity *mult = fVertexer->GetMultiplicity();
  if(mult)esd->SetMultiplicity(mult);

  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (fTracker[iDet]) fTracker[iDet]->SetVertex(vtxPos, vtxErr);
  }  
  delete vertex;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunHLTTracking(AliESDEvent*& esd)
{
// run the HLT barrel tracking

  AliCodeTimerAuto("")

  if (!fRunLoader) {
    AliError("Missing runLoader!");
    return kFALSE;
  }

  AliInfo("running HLT tracking");

  // Get a pointer to the HLT reconstructor
  AliReconstructor *reconstructor = GetReconstructor(fgkNDetectors-1);
  if (!reconstructor) return kFALSE;

  // TPC + ITS
  for (Int_t iDet = 1; iDet >= 0; iDet--) {
    TString detName = fgkDetectorName[iDet];
    AliDebug(1, Form("%s HLT tracking", detName.Data()));
    reconstructor->SetOption(detName.Data());
    AliTracker *tracker = reconstructor->CreateTracker();
    if (!tracker) {
      AliWarning(Form("couldn't create a HLT tracker for %s", detName.Data()));
      if (fStopOnError) return kFALSE;
      continue;
    }
    Double_t vtxPos[3];
    Double_t vtxErr[3]={0.005,0.005,0.010};
    const AliESDVertex *vertex = esd->GetVertex();
    vertex->GetXYZ(vtxPos);
    tracker->SetVertex(vtxPos,vtxErr);
    if(iDet != 1) {
      fLoader[iDet]->LoadRecPoints("read");
      TTree* tree = fLoader[iDet]->TreeR();
      if (!tree) {
	AliError(Form("Can't get the %s cluster tree", detName.Data()));
	return kFALSE;
      }
      tracker->LoadClusters(tree);
    }
    if (tracker->Clusters2Tracks(esd) != 0) {
      AliError(Form("HLT %s Clusters2Tracks failed", fgkDetectorName[iDet]));
      return kFALSE;
    }
    if(iDet != 1) {
      tracker->UnloadClusters();
    }
    delete tracker;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunMuonTracking(AliESDEvent*& esd)
{
// run the muon spectrometer tracking

  AliCodeTimerAuto("")

  if (!fRunLoader) {
    AliError("Missing runLoader!");
    return kFALSE;
  }
  Int_t iDet = 7; // for MUON

  AliInfo("is running...");

  // Get a pointer to the MUON reconstructor
  AliReconstructor *reconstructor = GetReconstructor(iDet);
  if (!reconstructor) return kFALSE;

  
  TString detName = fgkDetectorName[iDet];
  AliDebug(1, Form("%s tracking", detName.Data()));
  AliTracker *tracker =  reconstructor->CreateTracker();
  if (!tracker) {
    AliWarning(Form("couldn't create a tracker for %s", detName.Data()));
    return kFALSE;
  }
     
  // read RecPoints
  fLoader[iDet]->LoadRecPoints("read");  

  tracker->LoadClusters(fLoader[iDet]->TreeR());
  
  Int_t rv = tracker->Clusters2Tracks(esd);
  
  if ( rv )
  {
    AliError(Form("%s Clusters2Tracks failed", fgkDetectorName[iDet]));
    return kFALSE;
  }
  
  fLoader[iDet]->UnloadRecPoints();

  tracker->UnloadClusters();
  
  delete tracker;
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::RunTracking(AliESDEvent*& esd)
{
// run the barrel tracking
  static Int_t eventNr=0;
  AliCodeTimerAuto("")

  AliInfo("running tracking");

  //Fill the ESD with the T0 info (will be used by the TOF) 
  if (fReconstructor[11] && fLoader[11]) {
    fLoader[11]->LoadRecPoints("READ");
    TTree *treeR = fLoader[11]->TreeR();
    GetReconstructor(11)->FillESD((TTree *)NULL,treeR,esd);
  }

  // pass 1: TPC + ITS inwards
  for (Int_t iDet = 1; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s tracking", fgkDetectorName[iDet]));

    // load clusters
    fLoader[iDet]->LoadRecPoints("read");
    AliSysInfo::AddStamp(Form("RLoadCluster%s_%d",fgkDetectorName[iDet],eventNr),iDet,1, eventNr);
    TTree* tree = fLoader[iDet]->TreeR();
    if (!tree) {
      AliError(Form("Can't get the %s cluster tree", fgkDetectorName[iDet]));
      return kFALSE;
    }
    fTracker[iDet]->LoadClusters(tree);
    AliSysInfo::AddStamp(Form("TLoadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,2, eventNr);
    // run tracking
    if (fTracker[iDet]->Clusters2Tracks(esd) != 0) {
      AliError(Form("%s Clusters2Tracks failed", fgkDetectorName[iDet]));
      return kFALSE;
    }
    if (fCheckPointLevel > 1) {
      WriteESD(esd, Form("%s.tracking", fgkDetectorName[iDet]));
    }
    // preliminary PID in TPC needed by the ITS tracker
    if (iDet == 1) {
      GetReconstructor(1)->FillESD((TTree*)NULL, (TTree*)NULL, esd);
      AliESDpid::MakePID(esd);
    } 
    AliSysInfo::AddStamp(Form("Tracking0%s_%d",fgkDetectorName[iDet],eventNr), iDet,3,eventNr);
  }

  // pass 2: ALL backwards

  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s back propagation", fgkDetectorName[iDet]));

    // load clusters
    if (iDet > 1) {     // all except ITS, TPC
      TTree* tree = NULL;
      fLoader[iDet]->LoadRecPoints("read");
      AliSysInfo::AddStamp(Form("RLoadCluster0%s_%d",fgkDetectorName[iDet],eventNr), iDet,1, eventNr);
      tree = fLoader[iDet]->TreeR();
      if (!tree) {
	AliError(Form("Can't get the %s cluster tree", fgkDetectorName[iDet]));
	return kFALSE;
      }
      fTracker[iDet]->LoadClusters(tree); 
      AliSysInfo::AddStamp(Form("TLoadCluster0%s_%d",fgkDetectorName[iDet],eventNr), iDet,2, eventNr);
    }

    // run tracking
    if (iDet>1) // start filling residuals for the "outer" detectors
    if (fRunGlobalQA) AliTracker::SetFillResiduals(kTRUE);     

    if (fTracker[iDet]->PropagateBack(esd) != 0) {
      AliError(Form("%s backward propagation failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    if (fCheckPointLevel > 1) {
      WriteESD(esd, Form("%s.back", fgkDetectorName[iDet]));
    }

    // unload clusters
    if (iDet > 2) {     // all except ITS, TPC, TRD
      fTracker[iDet]->UnloadClusters();
      fLoader[iDet]->UnloadRecPoints();
    }
    // updated PID in TPC needed by the ITS tracker -MI
    if (iDet == 1) {
      GetReconstructor(1)->FillESD((TTree*)NULL, (TTree*)NULL, esd);
      AliESDpid::MakePID(esd);
    }
    AliSysInfo::AddStamp(Form("Tracking1%s_%d",fgkDetectorName[iDet],eventNr), iDet,3, eventNr);
  }
  //stop filling residuals for the "outer" detectors
  if (fRunGlobalQA) AliTracker::SetFillResiduals(kFALSE);     

  // write space-points to the ESD in case alignment data output
  // is switched on
  if (fWriteAlignmentData)
    WriteAlignmentData(esd);

  // pass 3: TRD + TPC + ITS refit inwards

  for (Int_t iDet = 2; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s inward refit", fgkDetectorName[iDet]));

    // run tracking
    if (iDet<2) // start filling residuals for TPC and ITS
    if (fRunGlobalQA) AliTracker::SetFillResiduals(kTRUE);     

    if (fTracker[iDet]->RefitInward(esd) != 0) {
      AliError(Form("%s inward refit failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    // run postprocessing
    if (fTracker[iDet]->PostProcess(esd) != 0) {
      AliError(Form("%s postprocessing failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    if (fCheckPointLevel > 1) {
      WriteESD(esd, Form("%s.refit", fgkDetectorName[iDet]));
    }
    AliSysInfo::AddStamp(Form("Tracking2%s_%d",fgkDetectorName[iDet],eventNr), iDet,3, eventNr);
    // unload clusters
    fTracker[iDet]->UnloadClusters();
    AliSysInfo::AddStamp(Form("TUnloadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,4, eventNr);
    fLoader[iDet]->UnloadRecPoints();
    AliSysInfo::AddStamp(Form("RUnloadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,5, eventNr);
  }
  // stop filling residuals for TPC and ITS
  if (fRunGlobalQA) AliTracker::SetFillResiduals(kFALSE);     

  eventNr++;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CleanESD(AliESDEvent *esd){
  //
  // Remove the data which are not needed for the physics analysis.
  //

  Int_t nTracks=esd->GetNumberOfTracks();
  Int_t nV0s=esd->GetNumberOfV0s();
  AliInfo
  (Form("Number of ESD tracks and V0s before cleaning: %d %d",nTracks,nV0s));

  Float_t cleanPars[]={fV0DCAmax,fV0CsPmin,fDmax,fZmax};
  Bool_t rc=esd->Clean(cleanPars);

  nTracks=esd->GetNumberOfTracks();
  nV0s=esd->GetNumberOfV0s();
  AliInfo
  (Form("Number of ESD tracks and V0s after cleaning %d %d",nTracks,nV0s));

  return rc;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillESD(AliESDEvent*& esd, const TString& detectors)
{
// fill the event summary data

  AliCodeTimerAuto("")
    static Int_t eventNr=0; 
  TString detStr = detectors;
  
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
  if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
    if (!ReadESD(esd, fgkDetectorName[iDet])) {
      AliDebug(1, Form("filling ESD for %s", fgkDetectorName[iDet]));
      TTree* clustersTree = NULL;
      if (fLoader[iDet]) {
	fLoader[iDet]->LoadRecPoints("read");
	clustersTree = fLoader[iDet]->TreeR();
	if (!clustersTree) {
	  AliError(Form("Can't get the %s clusters tree", 
			fgkDetectorName[iDet]));
	  if (fStopOnError) return kFALSE;
	}
      }
      if (fRawReader && !reconstructor->HasDigitConversion()) {
        reconstructor->FillESD(fRawReader, clustersTree, esd);
      } else {
	TTree* digitsTree = NULL;
	if (fLoader[iDet]) {
	  fLoader[iDet]->LoadDigits("read");
	  digitsTree = fLoader[iDet]->TreeD();
	  if (!digitsTree) {
	    AliError(Form("Can't get the %s digits tree", 
			  fgkDetectorName[iDet]));
	    if (fStopOnError) return kFALSE;
	  }
	}
	reconstructor->FillESD(digitsTree, clustersTree, esd);
	if (fLoader[iDet]) fLoader[iDet]->UnloadDigits();
      }
      if (fLoader[iDet]) {
	fLoader[iDet]->UnloadRecPoints();
      }

      if (fCheckPointLevel > 2) WriteESD(esd, fgkDetectorName[iDet]);
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s", 
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }
  AliSysInfo::AddStamp(Form("FillESD%d",eventNr), 0,1, eventNr);
  eventNr++;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillTriggerESD(AliESDEvent*& esd)
{
  // Reads the trigger decision which is
  // stored in Trigger.root file and fills
  // the corresponding esd entries

  AliCodeTimerAuto("")
  
  AliInfo("Filling trigger information into the ESD");

  AliCentralTrigger *aCTP = NULL;

  if (fRawReader) {
    AliCTPRawStream input(fRawReader);
    if (!input.Next()) {
      AliWarning("No valid CTP (trigger) DDL raw data is found ! The trigger mask will be taken from the event header, trigger cluster mask will be empty !");
      ULong64_t mask = (((ULong64_t)fRawReader->GetTriggerPattern()[1]) << 32) +
	fRawReader->GetTriggerPattern()[0];
      esd->SetTriggerMask(mask);
      esd->SetTriggerCluster(0);
    }
    else {
      esd->SetTriggerMask(input.GetClassMask());
      esd->SetTriggerCluster(input.GetClusterMask());
    }

    aCTP = new AliCentralTrigger();
    TString configstr("");
    if (!aCTP->LoadConfiguration(configstr)) { // Load CTP config from OCDB
      AliError("No trigger configuration found in OCDB! The trigger classes information will no be stored in ESD!");
      delete aCTP;
      return kFALSE;
    }
  }
  else {
    AliRunLoader *runloader = AliRunLoader::GetRunLoader();
    if (runloader) {
      if (!runloader->LoadTrigger()) {
	aCTP = runloader->GetTrigger();
	esd->SetTriggerMask(aCTP->GetClassMask());
	esd->SetTriggerCluster(aCTP->GetClusterMask());
      }
      else {
	AliWarning("No trigger can be loaded! The trigger information is not stored in the ESD !");
	return kFALSE;
      }
    }
    else {
      AliError("No run loader is available! The trigger information is not stored in the ESD !");
      return kFALSE;
    }
  }

  // Now fill the trigger class names into AliESDRun object
  AliTriggerConfiguration *config = aCTP->GetConfiguration();
  if (!config) {
    AliError("No trigger configuration has been found! The trigger classes information will not be stored in ESD!");
    if (fRawReader) delete aCTP;
    return kFALSE;
  }

  const TObjArray& classesArray = config->GetClasses();
  Int_t nclasses = classesArray.GetEntriesFast();
  for( Int_t j=0; j<nclasses; j++ ) {
    AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At( j );
    Int_t trindex = (Int_t)TMath::Log2(trclass->GetMask());
    esd->SetTriggerClass(trclass->GetName(),trindex);
  }

  if (fRawReader) delete aCTP;
  return kTRUE;
}





//_____________________________________________________________________________
Bool_t AliReconstruction::FillRawEventHeaderESD(AliESDEvent*& esd)
{
  // 
  // Filling information from RawReader Header
  // 

  AliInfo("Filling information from RawReader Header");
  esd->SetBunchCrossNumber(0);
  esd->SetOrbitNumber(0);
  esd->SetPeriodNumber(0);
  esd->SetTimeStamp(0);
  esd->SetEventType(0);
  const AliRawEventHeaderBase * eventHeader = fRawReader->GetEventHeader();
  if (eventHeader){

    const UInt_t *id = eventHeader->GetP("Id");
    esd->SetBunchCrossNumber((id)[1]&0x00000fff);
    esd->SetOrbitNumber((((id)[0]<<20)&0xf00000)|(((id)[1]>>12)&0xfffff));
    esd->SetPeriodNumber(((id)[0]>>4)&0x0fffffff);

    esd->SetTimeStamp((eventHeader->Get("Timestamp")));  
    esd->SetEventType((eventHeader->Get("Type")));
  }

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::IsSelected(TString detName, TString& detectors) const
{
// check whether detName is contained in detectors
// if yes, it is removed from detectors

  // check if all detectors are selected
  if ((detectors.CompareTo("ALL") == 0) ||
      detectors.BeginsWith("ALL ") ||
      detectors.EndsWith(" ALL") ||
      detectors.Contains(" ALL ")) {
    detectors = "ALL";
    return kTRUE;
  }

  // search for the given detector
  Bool_t result = kFALSE;
  if ((detectors.CompareTo(detName) == 0) ||
      detectors.BeginsWith(detName+" ") ||
      detectors.EndsWith(" "+detName) ||
      detectors.Contains(" "+detName+" ")) {
    detectors.ReplaceAll(detName, "");
    result = kTRUE;
  }

  // clean up the detectors string
  while (detectors.Contains("  ")) detectors.ReplaceAll("  ", " ");
  while (detectors.BeginsWith(" ")) detectors.Remove(0, 1);
  while (detectors.EndsWith(" ")) detectors.Remove(detectors.Length()-1, 1);

  return result;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::InitRunLoader()
{
// get or create the run loader

  if (gAlice) delete gAlice;
  gAlice = NULL;

  if (!gSystem->AccessPathName(fGAliceFileName.Data())) { // galice.root exists
    // load all base libraries to get the loader classes
    TString libs = gSystem->GetLibraries();
    for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
      TString detName = fgkDetectorName[iDet];
      if (detName == "HLT") continue;
      if (libs.Contains("lib" + detName + "base.so")) continue;
      gSystem->Load("lib" + detName + "base.so");
    }
    fRunLoader = AliRunLoader::Open(fGAliceFileName.Data());
    if (!fRunLoader) {
      AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }

    fRunLoader->CdGAFile();
    fRunLoader->LoadgAlice();

    //PH This is a temporary fix to give access to the kinematics
    //PH that is needed for the labels of ITS clusters
    fRunLoader->LoadHeader();
    fRunLoader->LoadKinematics();

  } else {               // galice.root does not exist
    if (!fRawReader) {
      AliError(Form("the file %s does not exist", fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }
    fRunLoader = AliRunLoader::Open(fGAliceFileName.Data(),
				    AliConfig::GetDefaultEventFolderName(),
				    "recreate");
    if (!fRunLoader) {
      AliError(Form("could not create run loader in file %s", 
		    fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }
    fIsNewRunLoader = kTRUE;
    fRunLoader->MakeTree("E");

    if (fNumberOfEventsPerFile > 0)
      fRunLoader->SetNumberOfEventsPerFile(fNumberOfEventsPerFile);
    else
      fRunLoader->SetNumberOfEventsPerFile((UInt_t)-1);
  }

  return kTRUE;
}

//_____________________________________________________________________________
AliReconstructor* AliReconstruction::GetReconstructor(Int_t iDet)
{
// get the reconstructor object and the loader for a detector

  if (fReconstructor[iDet]) return fReconstructor[iDet];

  // load the reconstructor object
  TPluginManager* pluginManager = gROOT->GetPluginManager();
  TString detName = fgkDetectorName[iDet];
  TString recName = "Ali" + detName + "Reconstructor";

  if (!fIsNewRunLoader && !fRunLoader->GetLoader(detName+"Loader") && (detName != "HLT")) return NULL;

  AliReconstructor* reconstructor = NULL;
  // first check if a plugin is defined for the reconstructor
  TPluginHandler* pluginHandler = 
    pluginManager->FindHandler("AliReconstructor", detName);
  // if not, add a plugin for it
  if (!pluginHandler) {
    AliDebug(1, Form("defining plugin for %s", recName.Data()));
    TString libs = gSystem->GetLibraries();
    if (libs.Contains("lib" + detName + "base.so") ||
	(gSystem->Load("lib" + detName + "base.so") >= 0)) {
      pluginManager->AddHandler("AliReconstructor", detName, 
				recName, detName + "rec", recName + "()");
    } else {
      pluginManager->AddHandler("AliReconstructor", detName, 
				recName, detName, recName + "()");
    }
    pluginHandler = pluginManager->FindHandler("AliReconstructor", detName);
  }
  if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
    reconstructor = (AliReconstructor*) pluginHandler->ExecPlugin(0);
  }
  if (reconstructor) {
    TObject* obj = fOptions.FindObject(detName.Data());
    if (obj) reconstructor->SetOption(obj->GetTitle());
    reconstructor->Init();
    fReconstructor[iDet] = reconstructor;
  }

  // get or create the loader
  if (detName != "HLT") {
    fLoader[iDet] = fRunLoader->GetLoader(detName + "Loader");
    if (!fLoader[iDet]) {
      AliConfig::Instance()
	->CreateDetectorFolders(fRunLoader->GetEventFolder(), 
				detName, detName);
      // first check if a plugin is defined for the loader
      pluginHandler = 
	pluginManager->FindHandler("AliLoader", detName);
      // if not, add a plugin for it
      if (!pluginHandler) {
	TString loaderName = "Ali" + detName + "Loader";
	AliDebug(1, Form("defining plugin for %s", loaderName.Data()));
	pluginManager->AddHandler("AliLoader", detName, 
				  loaderName, detName + "base", 
				  loaderName + "(const char*, TFolder*)");
	pluginHandler = pluginManager->FindHandler("AliLoader", detName);
      }
      if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
	fLoader[iDet] = 
	  (AliLoader*) pluginHandler->ExecPlugin(2, detName.Data(), 
						 fRunLoader->GetEventFolder());
      }
      if (!fLoader[iDet]) {   // use default loader
	fLoader[iDet] = new AliLoader(detName, fRunLoader->GetEventFolder());
      }
      if (!fLoader[iDet]) {
	AliWarning(Form("couldn't get loader for %s", detName.Data()));
	if (fStopOnError) return NULL;
      } else {
	fRunLoader->AddLoader(fLoader[iDet]);
	fRunLoader->CdGAFile();
	if (gFile && !gFile->IsWritable()) gFile->ReOpen("UPDATE");
	fRunLoader->Write(0, TObject::kOverwrite);
      }
    }
  }
      
  return reconstructor;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateVertexer()
{
// create the vertexer

  fVertexer = NULL;
  AliReconstructor* itsReconstructor = GetReconstructor(0);
  if (itsReconstructor) {
    fVertexer = itsReconstructor->CreateVertexer();
  }
  if (!fVertexer) {
    AliWarning("couldn't create a vertexer for ITS");
    if (fStopOnError) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateTrackers(const TString& detectors)
{
// create the trackers

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
    TString detName = fgkDetectorName[iDet];
    if (detName == "HLT") {
      fRunHLTTracking = kTRUE;
      continue;
    }
    if (detName == "MUON") {
      fRunMuonTracking = kTRUE;
      continue;
    }


    fTracker[iDet] = reconstructor->CreateTracker();
    if (!fTracker[iDet] && (iDet < 7)) {
      AliWarning(Form("couldn't create a tracker for %s", detName.Data()));
      if (fStopOnError) return kFALSE;
    }
    AliSysInfo::AddStamp(Form("LTracker%s",fgkDetectorName[iDet]), iDet,0);
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::CleanUp(TFile* file, TFile* fileOld)
{
// delete trackers and the run loader and close and delete the file

  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    delete fReconstructor[iDet];
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    delete fTracker[iDet];
    fTracker[iDet] = NULL;
//    delete fQADataMaker[iDet];
//    fQADataMaker[iDet] = NULL;
  }
  delete fVertexer;
  fVertexer = NULL;

  if (ftVertexer) delete ftVertexer;
  ftVertexer = NULL;
  
  if(!(AliCDBManager::Instance()->GetCacheFlag())) {
  	delete fDiamondProfile;
  	fDiamondProfile = NULL;
  	delete fDiamondProfileTPC;
  	fDiamondProfileTPC = NULL;
	delete fGRPData;
	fGRPData = NULL;
  }


  delete fRunLoader;
  fRunLoader = NULL;
  delete fRawReader;
  fRawReader = NULL;
  if (fParentRawReader) delete fParentRawReader;
  fParentRawReader=NULL;

  if (file) {
    file->Close();
    delete file;
  }

  if (fileOld) {
    fileOld->Close();
    delete fileOld;
    gSystem->Unlink("AliESDs.old.root");
  }

}

//_____________________________________________________________________________

Bool_t AliReconstruction::ReadESD(AliESDEvent*& esd, const char* recStep) const
{
// read the ESD event from a file

  if (!esd) return kFALSE;
  char fileName[256];
  sprintf(fileName, "ESD_%d.%d_%s.root", 
	  esd->GetRunNumber(), esd->GetEventNumberInFile(), recStep);
  if (gSystem->AccessPathName(fileName)) return kFALSE;

  AliInfo(Form("reading ESD from file %s", fileName));
  AliDebug(1, Form("reading ESD from file %s", fileName));
  TFile* file = TFile::Open(fileName);
  if (!file || !file->IsOpen()) {
    AliError(Form("opening %s failed", fileName));
    delete file;
    return kFALSE;
  }

  gROOT->cd();
  delete esd;
  esd = (AliESDEvent*) file->Get("ESD");
  file->Close();
  delete file;
  return kTRUE;

}



//_____________________________________________________________________________
void AliReconstruction::WriteESD(AliESDEvent* esd, const char* recStep) const
{
// write the ESD event to a file

  if (!esd) return;
  char fileName[256];
  sprintf(fileName, "ESD_%d.%d_%s.root", 
	  esd->GetRunNumber(), esd->GetEventNumberInFile(), recStep);

  AliDebug(1, Form("writing ESD to file %s", fileName));
  TFile* file = TFile::Open(fileName, "recreate");
  if (!file || !file->IsOpen()) {
    AliError(Form("opening %s failed", fileName));
  } else {
    esd->Write("ESD");
    file->Close();
  }
  delete file;
}


void AliReconstruction::WriteAlignmentData(AliESDEvent* esd)
{
  // Write space-points which are then used in the alignment procedures
  // For the moment only ITS, TRD and TPC

  // Load TOF clusters
  if (fTracker[3]){
    fLoader[3]->LoadRecPoints("read");
    TTree* tree = fLoader[3]->TreeR();
    if (!tree) {
      AliError(Form("Can't get the %s cluster tree", fgkDetectorName[3]));
      return;
    }
    fTracker[3]->LoadClusters(tree);
  }
  Int_t ntracks = esd->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < ntracks; itrack++)
    {
      AliESDtrack *track = esd->GetTrack(itrack);
      Int_t nsp = 0;
      Int_t idx[200];
      for (Int_t iDet = 3; iDet >= 0; iDet--)
	nsp += track->GetNcls(iDet);
      if (nsp) {
	AliTrackPointArray *sp = new AliTrackPointArray(nsp);
	track->SetTrackPointArray(sp);
	Int_t isptrack = 0;
	for (Int_t iDet = 3; iDet >= 0; iDet--) {
	  AliTracker *tracker = fTracker[iDet];
	  if (!tracker) continue;
	  Int_t nspdet = track->GetNcls(iDet);
	  if (nspdet <= 0) continue;
	  track->GetClusters(iDet,idx);
	  AliTrackPoint p;
	  Int_t isp = 0;
	  Int_t isp2 = 0;
	  while (isp2 < nspdet) {
	    Bool_t isvalid;
            TString dets = fgkDetectorName[iDet];
            if ((fUseTrackingErrorsForAlignment.CompareTo(dets) == 0) ||
            fUseTrackingErrorsForAlignment.BeginsWith(dets+" ") ||
            fUseTrackingErrorsForAlignment.EndsWith(" "+dets) ||
            fUseTrackingErrorsForAlignment.Contains(" "+dets+" ")) {
              isvalid = tracker->GetTrackPointTrackingError(idx[isp2],p,track);
	    } else {
	      isvalid = tracker->GetTrackPoint(idx[isp2],p); 
	    } 
	    isp2++;
	    const Int_t kNTPCmax = 159;
	    if (iDet==1 && isp2>kNTPCmax) break;   // to be fixed
	    if (!isvalid) continue;
	    sp->AddPoint(isptrack,&p); isptrack++; isp++;
	  }
	}	
      }
    }
  if (fTracker[3]){
    fTracker[3]->UnloadClusters();
    fLoader[3]->UnloadRecPoints();
  }
}

//_____________________________________________________________________________
void AliReconstruction::FillRawDataErrorLog(Int_t iEvent, AliESDEvent* esd)
{
  // The method reads the raw-data error log
  // accumulated within the rawReader.
  // It extracts the raw-data errors related to
  // the current event and stores them into
  // a TClonesArray inside the esd object.

  if (!fRawReader) return;

  for(Int_t i = 0; i < fRawReader->GetNumberOfErrorLogs(); i++) {

    AliRawDataErrorLog *log = fRawReader->GetErrorLog(i);
    if (!log) continue;
    if (iEvent != log->GetEventNumber()) continue;

    esd->AddRawDataErrorLog(log);
  }

}

TNamed* AliReconstruction::CopyFileToTNamed(TString fPath,TString pName){
  // Dump a file content into a char in TNamed
  ifstream in;
  in.open(fPath.Data(),ios::in | ios::binary|ios::ate);
  Int_t kBytes = (Int_t)in.tellg();
  printf("Size: %d \n",kBytes);
  TNamed *fn = 0;
  if(in.good()){
    char* memblock = new char [kBytes];
    in.seekg (0, ios::beg);
    in.read (memblock, kBytes);
    in.close();
    TString fData(memblock,kBytes);
    fn = new TNamed(pName,fData);
    printf("fData Size: %d \n",fData.Sizeof());
    printf("pName Size: %d \n",pName.Sizeof());
    printf("fn    Size: %d \n",fn->Sizeof());
    delete[] memblock;
  }
  else{
    AliInfo(Form("Could not Open %s\n",fPath.Data()));
  }

  return fn;
}

void AliReconstruction::TNamedToFile(TTree* fTree, TString pName){
  // This is not really needed in AliReconstruction at the moment
  // but can serve as a template

  TList *fList = fTree->GetUserInfo();
  TNamed *fn = (TNamed*)fList->FindObject(pName.Data());
  printf("fn Size: %d \n",fn->Sizeof());

  TString fTmp(fn->GetName()); // to be 100% sure in principle pName also works
  const char* cdata = fn->GetTitle();
  printf("fTmp Size %d\n",fTmp.Sizeof());

  int size = fn->Sizeof()-fTmp.Sizeof()-sizeof(UChar_t)-sizeof(Int_t); // see dfinition of TString::SizeOf()...
  printf("calculated size %d\n",size);
  ofstream out(pName.Data(),ios::out | ios::binary);
  out.write(cdata,size);
  out.close();

}
  
//_____________________________________________________________________________
AliQADataMakerRec * AliReconstruction::GetQADataMaker(Int_t iDet)
{
 // get the quality assurance data maker object and the loader for a detector

  if (fQADataMaker[iDet]) 
    return fQADataMaker[iDet];

  AliQADataMakerRec * qadm = NULL;
  if (iDet == fgkNDetectors) { //Global QA
     qadm = new AliGlobalQADataMaker();
     fQADataMaker[iDet] = qadm;
     return qadm;
  }

  // load the QA data maker object
  TPluginManager* pluginManager = gROOT->GetPluginManager();
  TString detName = fgkDetectorName[iDet];
  TString qadmName = "Ali" + detName + "QADataMakerRec";
  if (!fIsNewRunLoader && !fRunLoader->GetLoader(detName+"Loader") && (detName != "HLT")) 
    return NULL;

  // first check if a plugin is defined for the quality assurance data maker
  TPluginHandler* pluginHandler = pluginManager->FindHandler("AliQADataMakerRec", detName);
  // if not, add a plugin for it
  if (!pluginHandler) {
    AliDebug(1, Form("defining plugin for %s", qadmName.Data()));
    TString libs = gSystem->GetLibraries();
    if (libs.Contains("lib" + detName + "base.so") ||
	(gSystem->Load("lib" + detName + "base.so") >= 0)) {
      pluginManager->AddHandler("AliQADataMakerRec", detName, 
				qadmName, detName + "qadm", qadmName + "()");
    } else {
      pluginManager->AddHandler("AliQADataMakerRec", detName, 
				qadmName, detName, qadmName + "()");
    }
    pluginHandler = pluginManager->FindHandler("AliQADataMakerRec", detName);
  }
  if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
    qadm = (AliQADataMakerRec *) pluginHandler->ExecPlugin(0);
  }

  fQADataMaker[iDet] = qadm;

  return qadm;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunQA(AliESDEvent *& esd)
{
  // run the Quality Assurance data producer

  AliCodeTimerAuto("")
  TString detStr = fQADetectors ;
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
   if (!IsSelected(fgkDetectorName[iDet], detStr)) 
     continue;
   AliQADataMakerRec * qadm = GetQADataMaker(iDet);
   if (!qadm) 
     continue;
   AliCodeTimerStart(Form("running quality assurance data maker for %s", fgkDetectorName[iDet]));
   AliInfo(Form("running quality assurance data maker for %s", fgkDetectorName[iDet]));
    
   if (fQATasks.Contains(Form("%d", AliQA::kESDS))) {
	   qadm->Exec(AliQA::kESDS, esd) ; 
	   qadm->Increment() ; 
   }
   AliCodeTimerStop(Form("running quality assurance data maker for %s", fgkDetectorName[iDet]));
 }
 if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
   AliError(Form("the following detectors were not found: %s",
		 detStr.Data()));
   if (fStopOnError) 
     return kFALSE;
 }
 
 return kTRUE;
  
}

//_____________________________________________________________________________
void AliReconstruction::CheckQA()
{
// check the QA of SIM for this run and remove the detectors 
// with status Fatal
  
	TString newRunLocalReconstruction ; 
	TString newRunTracking ;
	TString newFillESD ;
	 
	for (Int_t iDet = 0; iDet < AliQA::kNDET; iDet++) {
		TString detName(AliQA::GetDetName(iDet)) ;
		AliQA * qa = AliQA::Instance(AliQA::DETECTORINDEX_t(iDet)) ; 
		if ( qa->IsSet(AliQA::DETECTORINDEX_t(iDet), AliQA::kSIM, AliQA::kFATAL)) {
				AliInfo(Form("QA status for %s in Hits and/or SDIGITS  and/or Digits was Fatal; No reconstruction performed", detName.Data())) ;
		} else {
			if ( fRunLocalReconstruction.Contains(AliQA::GetDetName(iDet)) || 
					fRunLocalReconstruction.Contains("ALL") )  {
				newRunLocalReconstruction += detName ; 
				newRunLocalReconstruction += " " ; 			
			}
			if ( fRunTracking.Contains(AliQA::GetDetName(iDet)) || 
					fRunTracking.Contains("ALL") )  {
				newRunTracking += detName ; 
				newRunTracking += " " ; 			
			}
			if ( fFillESD.Contains(AliQA::GetDetName(iDet)) || 
					fFillESD.Contains("ALL") )  {
				newFillESD += detName ; 
				newFillESD += " " ; 			
			}
		}
	}
	fRunLocalReconstruction = newRunLocalReconstruction ; 
	fRunTracking            = newRunTracking ; 
	fFillESD                = newFillESD ; 
}

//_____________________________________________________________________________
Int_t AliReconstruction::GetDetIndex(const char* detector)
{
  // return the detector index corresponding to detector
  Int_t index = -1 ; 
  for (index = 0; index < fgkNDetectors ; index++) {
    if ( strcmp(detector, fgkDetectorName[index]) == 0 )
	break ; 
  }	
  return index ; 
}
//_____________________________________________________________________________
Bool_t AliReconstruction::FinishPlaneEff() {
 //
 // Here execute all the necessary operationis, at the end of the tracking phase,
 // in case that evaluation of PlaneEfficiencies was required for some detector. 
 // E.g., write into a DataBase file the PlaneEfficiency which have been evaluated. 
 //
 // This Preliminary version works only FOR ITS !!!!!
 // other detectors (TOF,TRD, etc. have to develop their specific codes)
 //
 //  Input: none
 //  Return: kTRUE if all operations have been done properly, kFALSE otherwise 
 //
 Bool_t ret=kFALSE;
 //for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
 for (Int_t iDet = 0; iDet < 1; iDet++) { // for the time being only ITS  
   //if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
   if(fTracker[iDet]) {
      AliPlaneEff *planeeff=fTracker[iDet]->GetPlaneEff(); 
      ret=planeeff->WriteIntoCDB();
      if(planeeff->GetCreateHistos()) {
        TString name="PlaneEffHisto";
        name+=fgkDetectorName[iDet];
        name+=".root";
        ret*=planeeff->WriteHistosToFile(name,"RECREATE");
      }
   }
 }
 return ret;
}
//_____________________________________________________________________________
Bool_t AliReconstruction::InitPlaneEff() {
//
 // Here execute all the necessary operations, before of the tracking phase,
 // for the evaluation of PlaneEfficiencies, in case required for some detectors.
 // E.g., read from a DataBase file a first evaluation of the PlaneEfficiency 
 // which should be updated/recalculated.
 //
 // This Preliminary version will work only FOR ITS !!!!!
 // other detectors (TOF,TRD, etc. have to develop their specific codes)
 //
 //  Input: none
 //  Return: kTRUE if all operations have been done properly, kFALSE otherwise
 //
 AliWarning(Form("Implementation of this method not yet done !! Method return kTRUE"));
 return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::InitAliEVE()
{
  // This method should be called only in case 
  // AliReconstruction is run
  // within the alieve environment.
  // It will initialize AliEVE in a way
  // so that it can visualize event processed
  // by AliReconstruction.
  // The return flag shows whenever the
  // AliEVE initialization was successful or not.

  TString macroStr;
  macroStr.Form("%s/EVE/macros/alieve_online.C",gSystem->ExpandPathName("$ALICE_ROOT"));
  AliInfo(Form("Loading AliEVE macro: %s",macroStr.Data()));
  if (gROOT->LoadMacro(macroStr.Data()) != 0) return kFALSE;

  gROOT->ProcessLine("if (!gAliEveEvent) {gAliEveEvent = new AliEveEventManager();gAliEveEvent->SetAutoLoad(kTRUE);gAliEveEvent->AddNewEventCommand(\"alieve_online_on_new_event()\");gEve->AddEvent(gAliEveEvent);};");
  gROOT->ProcessLine("alieve_online_init()");

  return kTRUE;
}
  
//_____________________________________________________________________________
void AliReconstruction::RunAliEVE()
{
  // Runs AliEVE visualisation of
  // the current event.
  // Should be executed only after
  // successful initialization of AliEVE.

  AliInfo("Running AliEVE...");
  gROOT->ProcessLine(Form("gAliEveEvent->SetEvent((AliRunLoader*)%p,(AliRawReader*)%p,(AliESDEvent*)%p);",fRunLoader,fRawReader,fesd));
  gROOT->ProcessLine("gAliEveEvent->StartStopAutoLoadTimer();");
  gSystem->Run();
}

//_____________________________________________________________________________
Bool_t AliReconstruction::SetRunQA(TString detAndAction) 
{
	// Allows to run QA for a selected set of detectors
	// and a selected set of tasks among RAWS, RECPOINTS and ESDS
	// all selected detectors run the same selected tasks
	
	if (!detAndAction.Contains(":")) {
		AliError( Form("%s is a wrong syntax, use \"DetectorList:ActionList\" \n", detAndAction.Data()) ) ;
		fRunQA = kFALSE ;
		return kFALSE ; 		
	}
	Int_t colon = detAndAction.Index(":") ; 
	fQADetectors = detAndAction(0, colon) ; 
	if (fQADetectors.Contains("ALL") )
		fQADetectors = fFillESD ; 
	fQATasks   = detAndAction(colon+1, detAndAction.Sizeof() ) ; 
	if (fQATasks.Contains("ALL") ) {
		fQATasks = Form("%d %d %d", AliQA::kRAWS, AliQA::kRECPOINTS, AliQA::kESDS) ; 
	} else {
		fQATasks.ToUpper() ; 
		TString tempo("") ; 
		if ( fQATasks.Contains("RAW") ) 
			tempo = Form("%d ", AliQA::kRAWS) ; 
		if ( fQATasks.Contains("RECPOINT") ) 
			tempo += Form("%d ", AliQA::kRECPOINTS) ; 
		if ( fQATasks.Contains("ESD") ) 
			tempo += Form("%d ", AliQA::kESDS) ; 
		fQATasks = tempo ; 
		if (fQATasks.IsNull()) {
			AliInfo("No QA requested\n")  ;
			fRunQA = kFALSE ;
			return kTRUE ; 
		}
	}	
	TString tempo(fQATasks) ; 
    tempo.ReplaceAll(Form("%d", AliQA::kRAWS), AliQA::GetTaskName(AliQA::kRAWS)) 	;
    tempo.ReplaceAll(Form("%d", AliQA::kRECPOINTS), AliQA::GetTaskName(AliQA::kRECPOINTS)) ;	
    tempo.ReplaceAll(Form("%d", AliQA::kESDS), AliQA::GetTaskName(AliQA::kESDS)) ; 	
	AliInfo( Form("QA will be done on \"%s\" for \"%s\"\n", fQADetectors.Data(), tempo.Data()) ) ;  
	fRunQA = kTRUE ;
	return kTRUE; 
} 

	
