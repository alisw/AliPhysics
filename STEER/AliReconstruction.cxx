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
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TArrayD.h>
#include <TArrayF.h>
#include <TArrayS.h>
#include <TChain.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TPRegexp.h>
#include <TParameter.h>
#include <TPluginManager.h>
#include <TProof.h>
#include <TProofOutputFile.h>
#include <TROOT.h>
#include <TSystem.h>

#include "AliAlignObj.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCTPRawStream.h"
#include "AliCascadeVertexer.h"
#include "AliCentralTrigger.h"
#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliDetectorRecoParam.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDPmdTrack.h"
#include "AliESDTagCreator.h"
#include "AliESDVertex.h"
#include "AliESDcascade.h"
#include "AliESDfriend.h"
#include "AliESDkink.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrack.h"
#include "AliEventInfo.h"
#include "AliGRPObject.h"
#include "AliGRPRecoParam.h"
#include "AliGenEventHeader.h"
#include "AliGeomManager.h"
#include "AliGlobalQADataMaker.h" 
#include "AliHeader.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliPlaneEff.h"
#include "AliQA.h"
#include "AliQADataMakerRec.h" 
#include "AliQADataMakerSteer.h"
#include "AliRawEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawHLTManager.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"
#include "AliReconstruction.h"
#include "AliReconstructor.h"
#include "AliRun.h"
#include "AliRunInfo.h"
#include "AliRunLoader.h"
#include "AliSysInfo.h" // memory snapshots
#include "AliTrackPointArray.h"
#include "AliTracker.h"
#include "AliTriggerClass.h"
#include "AliTriggerCluster.h"
#include "AliTriggerConfiguration.h"
#include "AliV0vertexer.h"
#include "AliVertexer.h"
#include "AliVertexerTracks.h"

ClassImp(AliReconstruction)

//_____________________________________________________________________________
const char* AliReconstruction::fgkDetectorName[AliReconstruction::kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "ACORDE", "HLT"};

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const char* gAliceFilename) :
  TSelector(),
  fUniformField(kFALSE),
  fRunVertexFinder(kTRUE),
  fRunVertexFinderTracks(kTRUE),
  fRunHLTTracking(kFALSE),
  fRunMuonTracking(kFALSE),
  fRunV0Finder(kTRUE),
  fRunCascadeFinder(kTRUE),
  fStopOnError(kFALSE),
  fWriteAlignmentData(kFALSE),
  fWriteESDfriend(kFALSE),
  fFillTriggerESD(kTRUE),

  fCleanESD(kTRUE),
  fV0DCAmax(3.),
  fV0CsPmin(0.),
  fDmax(50.),
  fZmax(50.),

  fRunLocalReconstruction("ALL"),
  fRunTracking("ALL"),
  fFillESD("ALL"),
  fLoadCDB(""),
  fUseTrackingErrorsForAlignment(""),
  fGAliceFileName(gAliceFilename),
  fRawInput(""),
  fEquipIdMap(""),
  fFirstEvent(0),
  fLastEvent(-1),
  fNumberOfEventsPerFile((UInt_t)-1),
  fOptions(),
  fLoadAlignFromCDB(kTRUE),
  fLoadAlignData("ALL"),
  fUseHLTData(),
  fRunInfo(NULL),
  fEventInfo(),

  fRunLoader(NULL),
  fRawReader(NULL),
  fParentRawReader(NULL),

  fRecoParam(),

  fDiamondProfileSPD(NULL),
  fDiamondProfile(NULL),
  fDiamondProfileTPC(NULL),
  
  fGRPData(NULL),

  fAlignObjArray(NULL),
  fCDBUri(),
  fSpecCDBUri(), 
  fInitCDBCalled(kFALSE),
  fSetRunNumberFromDataCalled(kFALSE),
  fQADetectors("ALL"), 
  fQASteer(NULL),  
  fQATasks("ALL"), 
  fRunQA(kTRUE),  
  fRunGlobalQA(kTRUE),
  fSameQACycle(kFALSE),

  fRunPlaneEff(kFALSE),

  fesd(NULL),
  fhltesd(NULL),
  fesdf(NULL),
  ffile(NULL),
  ftree(NULL),
  fhlttree(NULL),
  ftVertexer(NULL),
  fIsNewRunLoader(kFALSE),
  fRunAliEVE(kFALSE),
  fChain(NULL)
{
// create reconstruction object with default parameters
  gGeoManager = NULL;
  
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
  }
  for (Int_t iDet = 0; iDet < AliQA::kNDET; iDet++) {
    fQACycles[iDet] = 999999 ;
    fQAWriteExpert[iDet] = kFALSE ; 
  }
    
  AliPID pid;
}

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const AliReconstruction& rec) :
  TSelector(),
  fUniformField(rec.fUniformField),
  fRunVertexFinder(rec.fRunVertexFinder),
  fRunVertexFinderTracks(rec.fRunVertexFinderTracks),
  fRunHLTTracking(rec.fRunHLTTracking),
  fRunMuonTracking(rec.fRunMuonTracking),
  fRunV0Finder(rec.fRunV0Finder),
  fRunCascadeFinder(rec.fRunCascadeFinder),
  fStopOnError(rec.fStopOnError),
  fWriteAlignmentData(rec.fWriteAlignmentData),
  fWriteESDfriend(rec.fWriteESDfriend),
  fFillTriggerESD(rec.fFillTriggerESD),

  fCleanESD(rec.fCleanESD),
  fV0DCAmax(rec.fV0DCAmax),
  fV0CsPmin(rec.fV0CsPmin),
  fDmax(rec.fDmax),
  fZmax(rec.fZmax),

  fRunLocalReconstruction(rec.fRunLocalReconstruction),
  fRunTracking(rec.fRunTracking),
  fFillESD(rec.fFillESD),
  fLoadCDB(rec.fLoadCDB),
  fUseTrackingErrorsForAlignment(rec.fUseTrackingErrorsForAlignment),
  fGAliceFileName(rec.fGAliceFileName),
  fRawInput(rec.fRawInput),
  fEquipIdMap(rec.fEquipIdMap),
  fFirstEvent(rec.fFirstEvent),
  fLastEvent(rec.fLastEvent),
  fNumberOfEventsPerFile(rec.fNumberOfEventsPerFile),
  fOptions(),
  fLoadAlignFromCDB(rec.fLoadAlignFromCDB),
  fLoadAlignData(rec.fLoadAlignData),
  fUseHLTData(rec.fUseHLTData),
  fRunInfo(NULL),
  fEventInfo(),

  fRunLoader(NULL),
  fRawReader(NULL),
  fParentRawReader(NULL),

  fRecoParam(rec.fRecoParam),

  fDiamondProfileSPD(rec.fDiamondProfileSPD),
  fDiamondProfile(rec.fDiamondProfile),
  fDiamondProfileTPC(rec.fDiamondProfileTPC),
  
  fGRPData(NULL),

  fAlignObjArray(rec.fAlignObjArray),
  fCDBUri(rec.fCDBUri),
  fSpecCDBUri(), 
  fInitCDBCalled(rec.fInitCDBCalled),
  fSetRunNumberFromDataCalled(rec.fSetRunNumberFromDataCalled),
  fQADetectors(rec.fQADetectors), 
  fQASteer(NULL),  
  fQATasks(rec.fQATasks), 
  fRunQA(rec.fRunQA),  
  fRunGlobalQA(rec.fRunGlobalQA),
  fSameQACycle(rec.fSameQACycle),
  fRunPlaneEff(rec.fRunPlaneEff),

  fesd(NULL),
  fhltesd(NULL),
  fesdf(NULL),
  ffile(NULL),
  ftree(NULL),
  fhlttree(NULL),
  ftVertexer(NULL),
  fIsNewRunLoader(rec.fIsNewRunLoader),
  fRunAliEVE(kFALSE),
  fChain(NULL)
{
// copy constructor

  for (Int_t i = 0; i < rec.fOptions.GetEntriesFast(); i++) {
    if (rec.fOptions[i]) fOptions.Add(rec.fOptions[i]->Clone());
  }
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
  }  
  
  for (Int_t iDet = 0; iDet < AliQA::kNDET; iDet++) {
    fQACycles[iDet] = rec.fQACycles[iDet];
    fQAWriteExpert[iDet] = rec.fQAWriteExpert[iDet] ; 
  }

  for (Int_t i = 0; i < rec.fSpecCDBUri.GetEntriesFast(); i++) {
    if (rec.fSpecCDBUri[i]) fSpecCDBUri.Add(rec.fSpecCDBUri[i]->Clone());
  }
}

//_____________________________________________________________________________
AliReconstruction& AliReconstruction::operator = (const AliReconstruction& rec)
{
// assignment operator
// Used in PROOF mode
// Be very careful while modifing it!
// Simple rules to follow:
// for persistent data members - use their assignment operators
// for non-persistent ones - do nothing or take the default values from constructor
// TSelector members should not be touched
  if(&rec == this) return *this;

  fUniformField          = rec.fUniformField;
  fRunVertexFinder       = rec.fRunVertexFinder;
  fRunVertexFinderTracks = rec.fRunVertexFinderTracks;
  fRunHLTTracking        = rec.fRunHLTTracking;
  fRunMuonTracking       = rec.fRunMuonTracking;
  fRunV0Finder           = rec.fRunV0Finder;
  fRunCascadeFinder      = rec.fRunCascadeFinder;
  fStopOnError           = rec.fStopOnError;
  fWriteAlignmentData    = rec.fWriteAlignmentData;
  fWriteESDfriend        = rec.fWriteESDfriend;
  fFillTriggerESD        = rec.fFillTriggerESD;

  fCleanESD  = rec.fCleanESD;
  fV0DCAmax  = rec.fV0DCAmax;
  fV0CsPmin  = rec.fV0CsPmin;
  fDmax      = rec.fDmax;
  fZmax      = rec.fZmax;

  fRunLocalReconstruction        = rec.fRunLocalReconstruction;
  fRunTracking                   = rec.fRunTracking;
  fFillESD                       = rec.fFillESD;
  fLoadCDB                       = rec.fLoadCDB;
  fUseTrackingErrorsForAlignment = rec.fUseTrackingErrorsForAlignment;
  fGAliceFileName                = rec.fGAliceFileName;
  fRawInput                      = rec.fRawInput;
  fEquipIdMap                    = rec.fEquipIdMap;
  fFirstEvent                    = rec.fFirstEvent;
  fLastEvent                     = rec.fLastEvent;
  fNumberOfEventsPerFile         = rec.fNumberOfEventsPerFile;

  for (Int_t i = 0; i < rec.fOptions.GetEntriesFast(); i++) {
    if (rec.fOptions[i]) fOptions.Add(rec.fOptions[i]->Clone());
  }

  fLoadAlignFromCDB              = rec.fLoadAlignFromCDB;
  fLoadAlignData                 = rec.fLoadAlignData;
  fUseHLTData                    = rec.fUseHLTData;

  delete fRunInfo; fRunInfo = NULL;
  if (rec.fRunInfo) fRunInfo = new AliRunInfo(*rec.fRunInfo);

  fEventInfo                     = rec.fEventInfo;

  fRunLoader       = NULL;
  fRawReader       = NULL;
  fParentRawReader = NULL;

  fRecoParam = rec.fRecoParam;

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    delete fReconstructor[iDet]; fReconstructor[iDet] = NULL;
    delete fLoader[iDet]; fLoader[iDet] = NULL;
    delete fTracker[iDet]; fTracker[iDet] = NULL;
  }
  
  for (Int_t iDet = 0; iDet < AliQA::kNDET; iDet++) {
    fQACycles[iDet] = rec.fQACycles[iDet];
    fQAWriteExpert[iDet] = rec.fQAWriteExpert[iDet] ;
  } 
    
  delete fDiamondProfileSPD; fDiamondProfileSPD = NULL;
  if (rec.fDiamondProfileSPD) fDiamondProfileSPD = new AliESDVertex(*rec.fDiamondProfileSPD);
  delete fDiamondProfile; fDiamondProfile = NULL;
  if (rec.fDiamondProfile) fDiamondProfile = new AliESDVertex(*rec.fDiamondProfile);
  delete fDiamondProfileTPC; fDiamondProfileTPC = NULL;
  if (rec.fDiamondProfileTPC) fDiamondProfileTPC = new AliESDVertex(*rec.fDiamondProfileTPC);

  delete fGRPData; fGRPData = NULL;
  //  if (rec.fGRPData) fGRPData = (TMap*)((rec.fGRPData)->Clone());
  if (rec.fGRPData) fGRPData = (AliGRPObject*)((rec.fGRPData)->Clone());

  delete fAlignObjArray; fAlignObjArray = NULL;

  fCDBUri        = "";
  fSpecCDBUri.Delete();
  fInitCDBCalled               = rec.fInitCDBCalled;
  fSetRunNumberFromDataCalled  = rec.fSetRunNumberFromDataCalled;
  fQADetectors                 = rec.fQADetectors;
  fQASteer                     = NULL;  
  fQATasks                     = rec.fQATasks; 
  fRunQA                       = rec.fRunQA;  
  fRunGlobalQA                 = rec.fRunGlobalQA;
  fSameQACycle                 = rec.fSameQACycle;
  fRunPlaneEff                 = rec.fRunPlaneEff;

  fesd     = NULL;
  fhltesd  = NULL;
  fesdf    = NULL;
  ffile    = NULL;
  ftree    = NULL;
  fhlttree = NULL;
  ftVertexer = NULL;
  fIsNewRunLoader = rec.fIsNewRunLoader;
  fRunAliEVE = kFALSE;
  fChain = NULL;

  return *this;
}

//_____________________________________________________________________________
AliReconstruction::~AliReconstruction()
{
// clean up

  CleanUp();
  delete fGRPData;
  fOptions.Delete();
  if (fAlignObjArray) {
    fAlignObjArray->Delete();
    delete fAlignObjArray;
  }
  fSpecCDBUri.Delete();
  delete fQASteer;
  AliCodeTimer::Instance()->Print();
}

//_____________________________________________________________________________
void AliReconstruction::InitCDB()
{
// activate a default CDB storage
// First check if we have any CDB storage set, because it is used 
// to retrieve the calibration and alignment constants
  AliCodeTimerAuto("");

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
  AliSysInfo::AddStamp("InitCDB");
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
	for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
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
//  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
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
 
  if(fRawReader) {
    if(fRawReader->NextEvent()) {
      if(man->GetRun() > 0) {
  	AliWarning("Run number is taken from raw-event header! Ignoring settings in AliCDBManager!");
      } 
      man->SetRun(fRawReader->GetRunNumber());
      fRawReader->RewindEvents();
    }
    else {
      if(man->GetRun() > 0) {
	AliWarning("No raw-data events are found ! Using settings in AliCDBManager !");
      }
      else {
	AliWarning("Neither raw events nor settings in AliCDBManager are found !");
	return kFALSE;
      }
    }
  }
  else {
    AliRunLoader *rl = AliRunLoader::Open(fGAliceFileName.Data());
    if (!rl) {
      AliError(Form("No run loader found in file %s", fGAliceFileName.Data()));
      return kFALSE;
    }
    else {
      rl->LoadHeader();
      // read run number from gAlice
      if(rl->GetHeader()) {
	man->SetRun(rl->GetHeader()->GetRun());
	rl->UnloadHeader();
	delete rl;
      }
      else {
	AliError("Neither run-loader header nor RawReader objects are found !");
	delete rl;
	return kFALSE;
      }
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
    
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if(!IsSelected(fgkDetectorName[iDet], detStr)) continue;
      if(!strcmp(fgkDetectorName[iDet],"HLT")) continue;
      
      if(AliGeomManager::GetNalignable(fgkDetectorName[iDet]) != 0)
      {
	loadAlObjsListOfDets += fgkDetectorName[iDet];
	loadAlObjsListOfDets += " ";
      }
    } // end loop over detectors
    
    if(AliGeomManager::GetNalignable("GRP") != 0)
      loadAlObjsListOfDets.Prepend("GRP "); //add alignment objects for non-sensitive modules
    AliGeomManager::ApplyAlignObjsFromCDB(loadAlObjsListOfDets.Data());
    AliCDBManager::Instance()->UnloadFromCache("*/Align/*");
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
  
  if (fAlignObjArray) {
    fAlignObjArray->Delete();
    delete fAlignObjArray; fAlignObjArray=NULL;
  }

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
  if (input) fRawInput = input;
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
void AliReconstruction::SetRecoParam(const char* detector, AliDetectorRecoParam *par)
{
  // Set custom reconstruction parameters for a given detector
  // Single set of parameters for all the events

  // First check if the reco-params are global
  if(!strcmp(detector, "GRP")) {
    par->SetAsDefault();
    fRecoParam.AddDetRecoParam(kNDetectors,par);
    return;
  }

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if(!strcmp(detector, fgkDetectorName[iDet])) {
      par->SetAsDefault();
      fRecoParam.AddDetRecoParam(iDet,par);
      break;
    }
  }

}

//_____________________________________________________________________________
Bool_t AliReconstruction::SetFieldMap(Float_t l3Cur, Float_t diCur, Float_t l3Pol, 
				      Float_t diPol, Float_t beamenergy, 
				      const Char_t *beamtype, const Char_t *path) 
{
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
  //
  TString s=(l3Pol < 0) ? "L3: -" : "L3: +";
  //
  AliMagF::BMap_t map = AliMagF::k5kG;
  //
  double fcL3,fcDip;
  //
  l3Cur = TMath::Abs(l3Cur);
  if (TMath::Abs(l3Cur-l3NominalCurrent1)/l3NominalCurrent1 < tolerance) {
    fcL3 = l3Cur/l3NominalCurrent1;
    map  = AliMagF::k5kG;
    s   += "0.5 T;  ";
  } else if (TMath::Abs(l3Cur-l3NominalCurrent2)/l3NominalCurrent2 < tolerance) {
    fcL3 = l3Cur/l3NominalCurrent2;
    map  = AliMagF::k2kG;
    s   += "0.2 T;  ";
  } else if (l3Cur <= zero) {
    fcL3 = 0;
    map  = AliMagF::k5kGUniform;
    s   += "0.0 T;  ";
    fUniformField=kTRUE;        // track with the uniform (zero) B field
  } else {
    AliError(Form("Wrong L3 current (%f A)!",l3Cur));
    return kFALSE;
  }
  //
  diCur = TMath::Abs(diCur);
  if (TMath::Abs(diCur-diNominalCurrent)/diNominalCurrent < tolerance) {
    // 3% current tolerance...
    fcDip = diCur/diNominalCurrent;
    s    += "Dipole ON";
  } else if (diCur <= zero) { // some small current..
    fcDip = 0.;
    s    += "Dipole OFF";
  } else {
    AliError(Form("Wrong dipole current (%f A)!",diCur));
    return kFALSE;
  }
  //
  if (l3Pol!=diPol && (map==AliMagF::k5kG || map==AliMagF::k2kG) && fcDip!=0) {
    AliError("L3 and Dipole polarities must be the same");
    return kFALSE;
  }
  //
  if (l3Pol<0) fcL3  = -fcL3;
  if (diPol<0) fcDip = -fcDip;
  //
  AliMagF::BeamType_t btype = AliMagF::kNoBeamField;
  TString btypestr = beamtype;
  btypestr.ToLower();
  TPRegexp protonBeam("(proton|p)\\s*-?\\s*\\1");
  TPRegexp ionBeam("(lead|pb|ion|a)\\s*-?\\s*\\1");
  if (btypestr.Contains(ionBeam)) btype = AliMagF::kBeamTypeAA;
  else if (btypestr.Contains(protonBeam)) btype = AliMagF::kBeamTypepp;
  else {
    AliInfo(Form("Cannot determine the beam type from %s, assume no LHC magnet field",beamtype));
  }
  
  AliMagF* fld = new AliMagF("MagneticFieldMap", s.Data(), 2, fcL3, fcDip, 10., map, path, 
			     btype,beamenergy,kTRUE);
  TGeoGlobalMagField::Instance()->SetField( fld );
  TGeoGlobalMagField::Instance()->Lock();
  //
  return kTRUE;
}


Bool_t AliReconstruction::InitGRP() {
  //------------------------------------
  // Initialization of the GRP entry 
  //------------------------------------
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
       AliInfo("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       fGRPData = new AliGRPObject();
       fGRPData->ReadValuesFromMap(m);
    }

    else {
       AliInfo("Found an AliGRPObject in GRP/GRP/Data, reading it");
       fGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }

    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if (!fGRPData) {
     AliError("No GRP entry found in OCDB!");
     return kFALSE;
  }

  TString lhcState = fGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = fGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = fGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }
  // energy is provided in MeV*120
  beamEnergy /= 120E3;

  TString runType = fGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = fGRPData->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidUInt()) {
    AliError("GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }

  fRunInfo = new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);
  printf("qqqqqqqqqqqqqqqqqqqqqqq %s %s %f %s %d\n", lhcState.Data(), beamType.Data(), beamEnergy, runType.Data(), activeDetectors);
  fRunInfo->Dump();


  // Process the list of active detectors
  if (activeDetectors) {
    UInt_t detMask = activeDetectors;
    fRunLocalReconstruction = MatchDetectorList(fRunLocalReconstruction,detMask);
    fRunTracking = MatchDetectorList(fRunTracking,detMask);
    fFillESD = MatchDetectorList(fFillESD,detMask);
    fQADetectors = MatchDetectorList(fQADetectors,detMask);
    fLoadCDB.Form("%s %s %s %s",
		  fRunLocalReconstruction.Data(),
		  fRunTracking.Data(),
		  fFillESD.Data(),
		  fQADetectors.Data());
    fLoadCDB = MatchDetectorList(fLoadCDB,detMask);
    if (!((detMask >> AliDAQ::DetectorID("ITSSPD")) & 0x1)) {
      // switch off the vertexer
      AliInfo("SPD is not in the list of active detectors. Vertexer switched off.");
      fRunVertexFinder = kFALSE;
    }
    if (!((detMask >> AliDAQ::DetectorID("TRG")) & 0x1)) {
      // switch off the reading of CTP raw-data payload
      if (fFillTriggerESD) {
	AliInfo("CTP is not in the list of active detectors. CTP data reading switched off.");
	fFillTriggerESD = kFALSE;
      }
    }
  }

  AliInfo("===================================================================================");
  AliInfo(Form("Running local reconstruction for detectors: %s",fRunLocalReconstruction.Data()));
  AliInfo(Form("Running tracking for detectors: %s",fRunTracking.Data()));
  AliInfo(Form("Filling ESD for detectors: %s",fFillESD.Data()));
  AliInfo(Form("Quality assurance is active for detectors: %s",fQADetectors.Data()));
  AliInfo(Form("CDB and reconstruction parameters are loaded for detectors: %s",fLoadCDB.Data()));
  AliInfo("===================================================================================");

  //*** Dealing with the magnetic field map
  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {AliInfo("Running with the externally locked B field !");}
  else {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = fGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = fGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = fGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = fGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    /*
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
    */

    if (ok) { 
      if ( !SetFieldMap(l3Current, diCurrent, l3Polarity ? -1:1, diPolarity ? -1:1) )
	AliFatal("Failed to creat a B field map ! Exiting...");
      AliInfo("Running with the B field constructed out of GRP !");
    }
    else AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
    
  }
  
  //*** Get the diamond profiles from OCDB
  entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexSPD");
  if (entry) {
    fDiamondProfileSPD = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  } else {
     AliError("No SPD diamond profile found in OCDB!");
  }

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
     AliError("No TPC diamond profile found in OCDB!");
  }

  return kTRUE;
} 

//_____________________________________________________________________________
Bool_t AliReconstruction::LoadCDB()
{
  AliCodeTimerAuto("");

  AliCDBManager::Instance()->Get("GRP/CTP/Config");

  TString detStr = fLoadCDB;
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliCDBManager::Instance()->GetAll(Form("%s/Calib/*",fgkDetectorName[iDet]));
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::Run(const char* input)
{
  // Run Run Run
  AliCodeTimerAuto("");

  InitRun(input);
  if (GetAbort() != TSelector::kContinue) return kFALSE;

  TChain *chain = NULL;
  if (fRawReader && (chain = fRawReader->GetChain())) {
    // Proof mode
    if (gProof) {
      gProof->AddInput(this);
      TUrl outputFile;
      outputFile.SetProtocol("root",kTRUE);
      outputFile.SetHost(gSystem->HostName());
      outputFile.SetFile(Form("%s/AliESDs.root",gSystem->pwd()));
      AliInfo(Form("Output file with ESDs is %s",outputFile.GetUrl()));
      gProof->AddInput(new TNamed("PROOF_OUTPUTFILE",outputFile.GetUrl()));
      chain->SetProof();
      chain->Process("AliReconstruction");
    }
    else {
      chain->Process(this);
    }
  }
  else {
    Begin(NULL);
    if (GetAbort() != TSelector::kContinue) return kFALSE;
    SlaveBegin(NULL);
    if (GetAbort() != TSelector::kContinue) return kFALSE;
    //******* The loop over events
    AliInfo("Starting looping over events");
    Int_t iEvent = 0;
    while ((iEvent < fRunLoader->GetNumberOfEvents()) ||
	   (fRawReader && fRawReader->NextEvent())) {
      if (!ProcessEvent(iEvent)) {
        Abort("ProcessEvent",TSelector::kAbortFile);
        return kFALSE;
      }
      iEvent++;
    }
    SlaveTerminate();
    if (GetAbort() != TSelector::kContinue) return kFALSE;
    Terminate();
    if (GetAbort() != TSelector::kContinue) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::InitRawReader(const char* input)
{
  AliCodeTimerAuto("");

  // Init raw-reader and
  // set the input in case of raw data
  if (input) fRawInput = input;
  fRawReader = AliRawReader::Create(fRawInput.Data());
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
  AliSysInfo::AddStamp("CreateRawReader");
}

//_____________________________________________________________________________
void AliReconstruction::InitRun(const char* input)
{
  // Initialization of raw-reader,
  // run number, CDB etc.
  AliCodeTimerAuto("");
  AliSysInfo::AddStamp("Start");

  // Initialize raw-reader if any
  InitRawReader(input);

  // Initialize the CDB storage
  InitCDB();

  // Set run number in CDBManager (if it is not already set by the user)
  if (!SetRunNumberFromData()) {
    Abort("SetRunNumberFromData", TSelector::kAbortProcess);
    return;
  }

  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  SetCDBLock();
  
}

//_____________________________________________________________________________
void AliReconstruction::Begin(TTree *)
{
  // Initialize AlReconstruction before
  // going into the event loop
  // Should follow the TSelector convention
  // i.e. initialize only the object on the client side
  AliCodeTimerAuto("");

  AliReconstruction *reco = NULL;
  if (fInput) {
    if ((reco = (AliReconstruction*)fInput->FindObject("AliReconstruction"))) {
      *this = *reco;
    }
    AliSysInfo::AddStamp("ReadInputInBegin");
  }

  // Import ideal TGeo geometry and apply misalignment
  if (!gGeoManager) {
    TString geom(gSystem->DirName(fGAliceFileName));
    geom += "/geometry.root";
    AliGeomManager::LoadGeometry(geom.Data());
    if (!gGeoManager) {
      Abort("LoadGeometry", TSelector::kAbortProcess);
      return;
    }
    AliSysInfo::AddStamp("LoadGeom");
    TString detsToCheck=fRunLocalReconstruction;
    if(!AliGeomManager::CheckSymNamesLUT(detsToCheck.Data())) {
      Abort("CheckSymNamesLUT", TSelector::kAbortProcess);
      return;
    }
    AliSysInfo::AddStamp("CheckGeom");
  }

  if (!MisalignGeometry(fLoadAlignData)) {
    Abort("MisalignGeometry", TSelector::kAbortProcess);
    return;
  }
  AliCDBManager::Instance()->UnloadFromCache("GRP/Geometry/Data");
  AliSysInfo::AddStamp("MisalignGeom");

  if (!InitGRP()) {
    Abort("InitGRP", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("InitGRP");

  if (!LoadCDB()) {
    Abort("LoadCDB", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("LoadCDB");

  // Read the reconstruction parameters from OCDB
  if (!InitRecoParams()) {
    AliWarning("Not all detectors have correct RecoParam objects initialized");
  }
  AliSysInfo::AddStamp("InitRecoParams");

  if (fInput) {
    if (reco) *reco = *this;
    fInput->Add(gGeoManager);
    gGeoManager = NULL;
    fInput->Add(const_cast<TMap*>(AliCDBManager::Instance()->GetEntryCache()));
    fInput->Add(new TParameter<Int_t>("RunNumber",AliCDBManager::Instance()->GetRun()));
    AliMagF *magFieldMap = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    magFieldMap->SetName("MagneticFieldMap");
    fInput->Add(magFieldMap);
  }

}

//_____________________________________________________________________________
void AliReconstruction::SlaveBegin(TTree*)
{
  // Initialization related to run-loader,
  // vertexer, trackers, recontructors
  // In proof mode it is executed on the slave
  AliCodeTimerAuto("");

  TProofOutputFile *outProofFile = NULL;
  if (fInput) { 
    if (AliReconstruction *reco = (AliReconstruction*)fInput->FindObject("AliReconstruction")) {
      *this = *reco;
    }
    if (TGeoManager *tgeo = (TGeoManager*)fInput->FindObject("Geometry")) {
      gGeoManager = tgeo;
      AliGeomManager::SetGeometry(tgeo);
    }
    if (TMap *entryCache = (TMap*)fInput->FindObject("CDBEntryCache")) {
      Int_t runNumber = -1;
      if (TProof::GetParameter(fInput,"RunNumber",runNumber) == 0) {
	AliCDBManager *man = AliCDBManager::Instance(entryCache,runNumber);
	man->SetCacheFlag(kTRUE);
	man->SetLock(kTRUE);
	man->Print();
      }
    }
    if (AliMagF *map = (AliMagF*)fInput->FindObject("MagneticFieldMap")) {
      TGeoGlobalMagField::Instance()->SetField(map);
    }
    if (TNamed *outputFileName = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE")) {
      outProofFile = new TProofOutputFile(gSystem->BaseName(TUrl(outputFileName->GetTitle()).GetFile()));
      outProofFile->SetOutputFileName(outputFileName->GetTitle());
      fOutput->Add(outProofFile);
    }
    AliSysInfo::AddStamp("ReadInputInSlaveBegin");
  }

  // get the run loader
  if (!InitRunLoader()) {
    Abort("InitRunLoader", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("LoadLoader");
 
  ftVertexer = new AliVertexerTracks(AliTracker::GetBz());

  // get trackers
  if (!fRunTracking.IsNull() && !CreateTrackers(fRunTracking)) {
    Abort("CreateTrackers", TSelector::kAbortProcess);
    return;
  }      
  AliSysInfo::AddStamp("CreateTrackers");

  // create the ESD output file and tree
  if (!outProofFile) {
    ffile = TFile::Open("AliESDs.root", "RECREATE");
    ffile->SetCompressionLevel(2);
    if (!ffile->IsOpen()) {
      Abort("OpenESDFile", TSelector::kAbortProcess);
      return;
    }
  }
  else {
    if (!(ffile = outProofFile->OpenFile("RECREATE"))) {
      Abort(Form("Problems opening output PROOF file: %s/%s",
		 outProofFile->GetDir(), outProofFile->GetFileName()),
	    TSelector::kAbortProcess);
      return;
    }
  }

  ftree = new TTree("esdTree", "Tree with ESD objects");
  fesd = new AliESDEvent();
  fesd->CreateStdContent();

  fesd->WriteToTree(ftree);
  if (fWriteESDfriend) {
    // careful:
    // Since we add the branch manually we must 
    // book and add it after WriteToTree
    // otherwise it is created twice,
    // once via writetotree and once here.
    // The case for AliESDfriend is now 
    // caught also in AlIESDEvent::WriteToTree but 
    // be careful when changing the name (AliESDfriend is not 
    // a TNamed so we had to hardwire it)
    fesdf = new AliESDfriend();
    TBranch *br=ftree->Branch("ESDfriend.","AliESDfriend", &fesdf);
    br->SetFile("AliESDfriends.root");
    fesd->AddObject(fesdf);
  }
  ftree->GetUserInfo()->Add(fesd);

  fhlttree = new TTree("HLTesdTree", "Tree with HLT ESD objects");
  fhltesd = new AliESDEvent();
  fhltesd->CreateStdContent();

  // read the ESD template from CDB
  // HLT is allowed to put non-std content to its ESD, the non-std
  // objects need to be created before invocation of WriteToTree in
  // order to create all branches. Initialization is done from an
  // ESD layout template in CDB
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBPath hltESDConfigPath("HLT/ConfigHLT/esdLayout");
  AliCDBEntry* hltESDConfig=NULL;
  if (man->GetId(hltESDConfigPath)!=NULL &&
      (hltESDConfig=man->Get(hltESDConfigPath))!=NULL) {
    AliESDEvent* pESDLayout=dynamic_cast<AliESDEvent*>(hltESDConfig->GetObject());
    if (pESDLayout) {
      // init all internal variables from the list of objects
      pESDLayout->GetStdContent();

      // copy content and create non-std objects
      *fhltesd=*pESDLayout;
      fhltesd->Reset();
    } else {
      AliError(Form("error setting hltEsd layout from %s: invalid object type",
		    hltESDConfigPath.GetPath().Data()));
    }
  }

  fhltesd->WriteToTree(fhlttree);
  fhlttree->GetUserInfo()->Add(fhltesd);

  ProcInfo_t procInfo;
  gSystem->GetProcInfo(&procInfo);
  AliInfo(Form("Current memory usage %d %d", procInfo.fMemResident, procInfo.fMemVirtual));
  
  //QA
  //Initialize the QA and start of cycle 
  if (fRunQA) {
    fQASteer = new AliQADataMakerSteer("rec") ; 
    fQASteer->SetActiveDetectors(fQADetectors) ; 
    for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
      fQASteer->SetCycleLength(AliQA::DETECTORINDEX_t(det), fQACycles[det]) ;  
      fQASteer->SetWriteExpert(AliQA::DETECTORINDEX_t(det)) ;
    }
    if (!fRawReader && fQATasks.Contains(AliQA::kRAWS))
      fQATasks.ReplaceAll(Form("%d",AliQA::kRAWS), "") ;
    fQASteer->SetTasks(fQATasks) ; 
    fQASteer->InitQADataMaker(AliCDBManager::Instance()->GetRun()) ; 
  }
  
  if (fRunGlobalQA) {
    Bool_t sameCycle = kFALSE ;
    if (!fQASteer) fQASteer = new AliQADataMakerSteer("rec") ; 
    AliQADataMaker *qadm = fQASteer->GetQADataMaker(AliQA::kGLOBAL);
    AliInfo(Form("Initializing the global QA data maker"));
    if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) {
      qadm->StartOfCycle(AliQA::kRECPOINTS, AliCDBManager::Instance()->GetRun(), sameCycle) ; 
      TObjArray **arr=qadm->Init(AliQA::kRECPOINTS);
      AliTracker::SetResidualsArray(arr);
      sameCycle = kTRUE ; 
    }
    if (fQATasks.Contains(Form("%d", AliQA::kESDS))) {
      qadm->StartOfCycle(AliQA::kESDS, AliCDBManager::Instance()->GetRun(), sameCycle) ; 
      qadm->Init(AliQA::kESDS);
    }
  }

  //Initialize the Plane Efficiency framework
  if (fRunPlaneEff && !InitPlaneEff()) {
    Abort("InitPlaneEff", TSelector::kAbortProcess);
    return;
  }

  if (strcmp(gProgName,"alieve") == 0)
    fRunAliEVE = InitAliEVE();

  return;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::Process(Long64_t entry)
{
  // run the reconstruction over a single entry
  // from the chain with raw data
  AliCodeTimerAuto("");

  TTree *currTree = fChain->GetTree();
  AliRawEvent *event = new AliRawEvent;
  currTree->SetBranchAddress("rawevent",&event);
  currTree->GetEntry(entry);
  fRawReader = new AliRawReaderRoot(event);
  fStatus = ProcessEvent(fRunLoader->GetNumberOfEvents());  
  delete fRawReader;
  fRawReader = NULL;
  delete event;

  return fStatus;
}

//_____________________________________________________________________________
void AliReconstruction::Init(TTree *tree)
{
  if (tree == 0) {
    AliError("The input tree is not found!");
    return;
  }
  fChain = tree;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::ProcessEvent(Int_t iEvent)
{
  // run the reconstruction over a single event
  // The event loop is steered in Run method

  AliCodeTimerAuto("");

  if (iEvent >= fRunLoader->GetNumberOfEvents()) {
    fRunLoader->SetEventNumber(iEvent);
    fRunLoader->GetHeader()->Reset(fRawReader->GetRunNumber(), 
				   iEvent, iEvent);
    fRunLoader->TreeE()->Fill();
    if (fRawReader && fRawReader->UseAutoSaveESD())
      fRunLoader->TreeE()->AutoSave("SaveSelf");
  }

  if ((iEvent < fFirstEvent) || ((fLastEvent >= 0) && (iEvent > fLastEvent))) {
    return kTRUE;
  }

  AliInfo(Form("processing event %d", iEvent));

  fRunLoader->GetEvent(iEvent);

  // Fill Event-info object
  GetEventInfo();
  fRecoParam.SetEventSpecie(fRunInfo,fEventInfo);
  AliInfo(Form("Current event specie: %s",fRecoParam.PrintEventSpecie()));

  // Set the reco-params
  {
    TString detStr = fLoadCDB;
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
      AliReconstructor *reconstructor = GetReconstructor(iDet);
      if (reconstructor && fRecoParam.GetDetRecoParamArray(iDet)) {
        const AliDetectorRecoParam *par = fRecoParam.GetDetRecoParam(iDet);
        reconstructor->SetRecoParam(par);
        if (fRunQA) {
          fQASteer->SetRecoParam(iDet, par) ; 
        }
      }
    }
  }

    // QA on single raw 
  if (fRunQA) {
    fQASteer->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
    fQASteer->RunOneEvent(fRawReader) ;  
  }
    // local single event reconstruction
    if (!fRunLocalReconstruction.IsNull()) {
      TString detectors=fRunLocalReconstruction;
      // run HLT event reconstruction first
      // ;-( IsSelected changes the string
      if (IsSelected("HLT", detectors) &&
	  !RunLocalEventReconstruction("HLT")) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
      detectors=fRunLocalReconstruction;
      detectors.ReplaceAll("HLT", "");
      if (!RunLocalEventReconstruction(detectors)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }

    fesd->SetRunNumber(fRunLoader->GetHeader()->GetRun());
    fhltesd->SetRunNumber(fRunLoader->GetHeader()->GetRun());
    fesd->SetEventNumberInFile(fRunLoader->GetHeader()->GetEventNrInRun());
    fhltesd->SetEventNumberInFile(fRunLoader->GetHeader()->GetEventNrInRun());
    
    // Set magnetic field from the tracker
    fesd->SetMagneticField(AliTracker::GetBz());
    fhltesd->SetMagneticField(AliTracker::GetBz());

    // Set most probable pt, for B=0 tracking
    // Get the global reco-params. They are atposition 16 inside the array of detectors in fRecoParam
    const AliGRPRecoParam *grpRecoParam = dynamic_cast<const AliGRPRecoParam*>(fRecoParam.GetDetRecoParam(kNDetectors));
    if (grpRecoParam) AliExternalTrackParam::SetMostProbablePt(grpRecoParam->GetMostProbablePt());
    
    // Fill raw-data error log into the ESD
    if (fRawReader) FillRawDataErrorLog(iEvent,fesd);

    // vertex finder
    if (fRunVertexFinder) {
      if (!RunVertexFinder(fesd)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }

    // Muon tracking
    if (!fRunTracking.IsNull()) {
      if (fRunMuonTracking) {
	if (!RunMuonTracking(fesd)) {
	  if (fStopOnError) {CleanUp(); return kFALSE;}
	}
      }
    }

    // barrel tracking
    if (!fRunTracking.IsNull()) {
      if (!RunTracking(fesd)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }

    // fill ESD
    if (!fFillESD.IsNull()) {
      TString detectors=fFillESD;
      // run HLT first and on hltesd
      // ;-( IsSelected changes the string
      if (IsSelected("HLT", detectors) &&
	  !FillESD(fhltesd, "HLT")) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
      detectors=fFillESD;
      // Temporary fix to avoid problems with HLT that overwrites the offline ESDs
      if (detectors.Contains("ALL")) {
	detectors="";
	for (Int_t idet=0; idet<kNDetectors; ++idet){
	  detectors += fgkDetectorName[idet];
	  detectors += " ";
	}
      }
      detectors.ReplaceAll("HLT", "");
      if (!FillESD(fesd, detectors)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }
  
    // fill Event header information from the RawEventHeader
    if (fRawReader){FillRawEventHeaderESD(fesd);}

    // combined PID
    AliESDpid::MakePID(fesd);

    if (fFillTriggerESD) {
      if (!FillTriggerESD(fesd)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
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
      const Double_t kMaxStep = 1;   //max step over the material
      Bool_t ok;

      AliESDtrack *track = fesd->GetTrack(itrack);
      if (!track) continue;

      AliExternalTrackParam *tpcTrack =
           (AliExternalTrackParam *)track->GetTPCInnerParam();
      ok = kFALSE;
      if (tpcTrack)
	ok = AliTracker::
	  PropagateTrackTo(tpcTrack,kRadius,track->GetMass(),kMaxStep,kFALSE);

      if (ok) {
	Int_t n=trkArray.GetEntriesFast();
        selectedIdx[n]=track->GetID();
        trkArray.AddLast(tpcTrack);
      }

      //Tracks refitted by ITS should already be at the SPD vertex
      if (track->IsOn(AliESDtrack::kITSrefit)) continue;

      AliTracker::
         PropagateTrackTo(track,kRadius,track->GetMass(),kMaxStep,kFALSE);
      track->RelateToVertex(fesd->GetPrimaryVertexSPD(), kBz, kVeryBig);

    }

    //
    // Improve the reconstructed primary vertex position using the tracks
    //
    Bool_t runVertexFinderTracks = fRunVertexFinderTracks;
    if(fesd->GetPrimaryVertexSPD()) {
      TString vtitle = fesd->GetPrimaryVertexSPD()->GetTitle();
      if(vtitle.Contains("cosmics")) {
	runVertexFinderTracks=kFALSE;
      }
    }

    if (runVertexFinderTracks) {
       // TPC + ITS primary vertex
       ftVertexer->SetITSMode();
       ftVertexer->SetConstraintOff();
       // get cuts for vertexer from AliGRPRecoParam
       if (grpRecoParam) {
	 Int_t nCutsVertexer = grpRecoParam->GetVertexerTracksNCuts();
	 Double_t *cutsVertexer = new Double_t[nCutsVertexer];
	 grpRecoParam->GetVertexerTracksCutsITS(cutsVertexer);
	 ftVertexer->SetCuts(cutsVertexer);
	 delete [] cutsVertexer; cutsVertexer = NULL; 
	 if(fDiamondProfile && grpRecoParam->GetVertexerTracksConstraintITS())
	   ftVertexer->SetVtxStart(fDiamondProfile);
       }
       AliESDVertex *pvtx=ftVertexer->FindPrimaryVertex(fesd);
       if (pvtx) {
          if (pvtx->GetStatus()) {
             fesd->SetPrimaryVertexTracks(pvtx);
             for (Int_t i=0; i<ntracks; i++) {
	         AliESDtrack *t = fesd->GetTrack(i);
                 t->RelateToVertex(pvtx, kBz, kVeryBig);
             } 
          }
       }

       // TPC-only primary vertex
       ftVertexer->SetTPCMode();
       ftVertexer->SetConstraintOff();
       // get cuts for vertexer from AliGRPRecoParam
       if (grpRecoParam) {
	 Int_t nCutsVertexer = grpRecoParam->GetVertexerTracksNCuts();
	 Double_t *cutsVertexer = new Double_t[nCutsVertexer];
	 grpRecoParam->GetVertexerTracksCutsTPC(cutsVertexer);
	 ftVertexer->SetCuts(cutsVertexer);
	 delete [] cutsVertexer; cutsVertexer = NULL; 
	 if(fDiamondProfileTPC && grpRecoParam->GetVertexerTracksConstraintTPC())
	   ftVertexer->SetVtxStart(fDiamondProfileTPC);
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
    fQASteer->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
    fQASteer->RunOneEvent(fesd) ; 
  }
  if (fRunGlobalQA) {
      AliQADataMaker *qadm = fQASteer->GetQADataMaker(AliQA::kGLOBAL);
      qadm->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
      if (qadm && fQATasks.Contains(Form("%d", AliQA::kESDS)))
        qadm->Exec(AliQA::kESDS, fesd);
    }

    if (fWriteESDfriend) {
      //      fesdf->~AliESDfriend();
      //  new (fesdf) AliESDfriend(); // Reset...
      fesd->GetESDfriend(fesdf);
    }
    ftree->Fill();

    // Auto-save the ESD tree in case of prompt reco @P2
    if (fRawReader && fRawReader->UseAutoSaveESD()) {
      ftree->AutoSave("SaveSelf");
      TFile *friendfile = (TFile *)(gROOT->GetListOfFiles()->FindObject("AliESDfriends.root"));
      if (friendfile) friendfile->Save();
    }

    // write HLT ESD
    fhlttree->Fill();

    // call AliEVE
    if (fRunAliEVE) RunAliEVE();

    fesd->Reset();
    fhltesd->Reset();
    if (fWriteESDfriend) {
      fesdf->~AliESDfriend();
      new (fesdf) AliESDfriend(); // Reset...
    }
 
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);
    AliInfo(Form("Event %d -> Current memory usage %d %d",iEvent, procInfo.fMemResident, procInfo.fMemVirtual));
  
    fEventInfo.Reset();
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (fReconstructor[iDet])
	fReconstructor[iDet]->SetRecoParam(NULL);
    }
	
  if (fRunQA || fRunGlobalQA) 
      fQASteer->Increment() ; 
  
    return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::SlaveTerminate()
{
  // Finalize the run on the slave side
  // Called after the exit
  // from the event loop
  AliCodeTimerAuto("");

  if (fIsNewRunLoader) { // galice.root didn't exist
    fRunLoader->WriteHeader("OVERWRITE");
    fRunLoader->CdGAFile();
    fRunLoader->Write(0, TObject::kOverwrite);
  }

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

  // End of cycle for the in-loop  
  if (fRunQA) 
    fQASteer->EndOfCycle() ;
  
  if (fRunGlobalQA) {
    AliQADataMaker *qadm = fQASteer->GetQADataMaker(AliQA::kGLOBAL);
    if (qadm) {
      if (fQATasks.Contains(Form("%d", AliQA::kRECPOINTS))) 
        qadm->EndOfCycle(AliQA::kRECPOINTS);
      if (fQATasks.Contains(Form("%d", AliQA::kESDS))) 
        qadm->EndOfCycle(AliQA::kESDS);
      qadm->Finish();
    }
  }
  gROOT->cd();
  CleanUp();
}
    
//_____________________________________________________________________________
void AliReconstruction::Terminate()
{
  // Create tags for the events in the ESD tree (the ESD tree is always present)
  // In case of empty events the tags will contain dummy values
  AliCodeTimerAuto("");

  AliESDTagCreator *esdtagCreator = new AliESDTagCreator();
  esdtagCreator->CreateESDTags(fFirstEvent,fLastEvent,fGRPData, AliQA::Instance()->GetQA(), AliQA::Instance()->GetEventSpecies(), AliQA::kNDET, AliRecoParam::kNSpecies);

  // Cleanup of CDB manager: cache and active storages!
  AliCDBManager::Instance()->ClearCache();
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunLocalEventReconstruction(const TString& detectors)
{
// run the local reconstruction

  static Int_t eventNr=0;
  AliCodeTimerAuto("")

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
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
//      AliCodeTimerAuto(Form("converting raw data digits into root objects for %s", 
//                            fgkDetectorName[iDet]));
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
    //AliCodeTimerAuto(Form("running reconstruction for %s", fgkDetectorName[iDet]));
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

		TString detQAStr(fQADetectors) ; 
		if (fRunQA) {
      fQASteer->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
			fQASteer->RunOneEventInOneDetector(iDet, clustersTree) ; 
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

  AliVertexer *vertexer = CreateVertexer();
  if (!vertexer) return kFALSE;

  AliInfo("running the ITS vertex finder");
  AliESDVertex* vertex = NULL;
  if (fLoader[0]) {
    fLoader[0]->LoadRecPoints();
    TTree* cltree = fLoader[0]->TreeR();
    if (cltree) {
      if(fDiamondProfileSPD) vertexer->SetVtxStart(fDiamondProfileSPD);
      vertex = vertexer->FindVertexForCurrentEvent(cltree);
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

  Double_t vtxPos[3];
  Double_t vtxErr[3];
  vertex->GetXYZ(vtxPos);
  vertex->GetSigmaXYZ(vtxErr);

  esd->SetPrimaryVertexSPD(vertex);
  // if SPD multiplicity has been determined, it is stored in the ESD
  AliMultiplicity *mult = vertexer->GetMultiplicity();
  if(mult)esd->SetMultiplicity(mult);

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (fTracker[iDet]) fTracker[iDet]->SetVertex(vtxPos, vtxErr);
  }  
  delete vertex;

  delete vertexer;

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
  AliReconstructor *reconstructor = GetReconstructor(kNDetectors-1);
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
    if (treeR) {
      GetReconstructor(11)->FillESD((TTree *)NULL,treeR,esd);
    }
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
    // preliminary PID in TPC needed by the ITS tracker
    if (iDet == 1) {
      GetReconstructor(1)->FillESD((TTree*)NULL, (TTree*)NULL, esd);
      AliESDpid::MakePID(esd);
    } 
    AliSysInfo::AddStamp(Form("Tracking0%s_%d",fgkDetectorName[iDet],eventNr), iDet,3,eventNr);
  }

  // pass 2: ALL backwards

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
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
    if (fRunGlobalQA) AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kTRUE);     

    if (fTracker[iDet]->PropagateBack(esd) != 0) {
      AliError(Form("%s backward propagation failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }

    // unload clusters
    if (iDet > 3) {     // all except ITS, TPC, TRD and TOF
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
  if (fRunGlobalQA) AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kFALSE);     

  // pass 3: TRD + TPC + ITS refit inwards

  for (Int_t iDet = 2; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s inward refit", fgkDetectorName[iDet]));

    // run tracking
    if (iDet<2) // start filling residuals for TPC and ITS
    if (fRunGlobalQA) AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kTRUE);     

    if (fTracker[iDet]->RefitInward(esd) != 0) {
      AliError(Form("%s inward refit failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    // run postprocessing
    if (fTracker[iDet]->PostProcess(esd) != 0) {
      AliError(Form("%s postprocessing failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    AliSysInfo::AddStamp(Form("Tracking2%s_%d",fgkDetectorName[iDet],eventNr), iDet,3, eventNr);
  }

  // write space-points to the ESD in case alignment data output
  // is switched on
  if (fWriteAlignmentData)
    WriteAlignmentData(esd);

  for (Int_t iDet = 3; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    // unload clusters
    fTracker[iDet]->UnloadClusters();
    AliSysInfo::AddStamp(Form("TUnloadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,4, eventNr);
    fLoader[iDet]->UnloadRecPoints();
    AliSysInfo::AddStamp(Form("RUnloadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,5, eventNr);
  }
  // stop filling residuals for TPC and ITS
  if (fRunGlobalQA) AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kFALSE);     

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
  
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
  if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
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

  if (fRawReader) {
    AliCTPRawStream input(fRawReader);
    if (!input.Next()) {
      AliWarning("No valid CTP (trigger) DDL raw data is found ! The trigger info is taken from the event header!");
    }
    else {
      if (esd->GetTriggerMask() != input.GetClassMask())
	AliError(Form("Invalid trigger pattern found in CTP raw-data: %llx %llx",
		      input.GetClassMask(),esd->GetTriggerMask()));
      if (esd->GetOrbitNumber() != input.GetOrbitID())
	AliError(Form("Invalid orbit id found in CTP raw-data: %x %x",
		      input.GetOrbitID(),esd->GetOrbitNumber()));
      if (esd->GetBunchCrossNumber() != input.GetBCID())
	AliError(Form("Invalid bunch-crossing id found in CTP raw-data: %x %x",
		      input.GetBCID(),esd->GetBunchCrossNumber()));
    }

  // Here one has to add the filling of trigger inputs and
  // interaction records
  // ...
  }
  return kTRUE;
}





//_____________________________________________________________________________
Bool_t AliReconstruction::FillRawEventHeaderESD(AliESDEvent*& esd)
{
  // 
  // Filling information from RawReader Header
  // 

  if (!fRawReader) return kFALSE;

  AliInfo("Filling information from RawReader Header");

  esd->SetBunchCrossNumber(fRawReader->GetBCID());
  esd->SetOrbitNumber(fRawReader->GetOrbitID());
  esd->SetPeriodNumber(fRawReader->GetPeriod());

  esd->SetTimeStamp(fRawReader->GetTimestamp());  
  esd->SetEventType(fRawReader->GetType());

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
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
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

  if (fReconstructor[iDet]) {
    if (fRecoParam.GetDetRecoParamArray(iDet) && !AliReconstructor::GetRecoParam(iDet)) {
      const AliDetectorRecoParam *par = fRecoParam.GetDetRecoParam(iDet);
      fReconstructor[iDet]->SetRecoParam(par);
    }
    return fReconstructor[iDet];
  }

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
      
  if (fRecoParam.GetDetRecoParamArray(iDet) && !AliReconstructor::GetRecoParam(iDet)) {
    const AliDetectorRecoParam *par = fRecoParam.GetDetRecoParam(iDet);
    reconstructor->SetRecoParam(par);
  }
  return reconstructor;
}

//_____________________________________________________________________________
AliVertexer* AliReconstruction::CreateVertexer()
{
// create the vertexer
// Please note that the caller is the owner of the
// vertexer

  AliVertexer* vertexer = NULL;
  AliReconstructor* itsReconstructor = GetReconstructor(0);
  if (itsReconstructor) {
    vertexer = itsReconstructor->CreateVertexer();
  }
  if (!vertexer) {
    AliWarning("couldn't create a vertexer for ITS");
  }

  return vertexer;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateTrackers(const TString& detectors)
{
// create the trackers
	AliInfo("Creating trackers");

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
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
void AliReconstruction::CleanUp()
{
// delete trackers and the run loader and close and delete the file

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    delete fReconstructor[iDet];
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    delete fTracker[iDet];
    fTracker[iDet] = NULL;
  }
  delete fRunInfo;
  fRunInfo = NULL;

  delete ftVertexer;
  ftVertexer = NULL;
  
  if(!(AliCDBManager::Instance()->GetCacheFlag())) {
    delete fDiamondProfileSPD;
    fDiamondProfileSPD = NULL;
    delete fDiamondProfile;
    fDiamondProfile = NULL;
    delete fDiamondProfileTPC;
    fDiamondProfileTPC = NULL;
  }
  delete fRunLoader;
  fRunLoader = NULL;
  delete fRawReader;
  fRawReader = NULL;
  delete fParentRawReader;
  fParentRawReader=NULL;

  if (ffile) {
    ffile->Close();
    delete ffile;
    ffile = NULL;
  }
}

void AliReconstruction::WriteAlignmentData(AliESDEvent* esd)
{
  // Write space-points which are then used in the alignment procedures
  // For the moment only ITS, TPC, TRD and TOF

  Int_t ntracks = esd->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < ntracks; itrack++)
    {
      AliESDtrack *track = esd->GetTrack(itrack);
      Int_t nsp = 0;
      Int_t idx[200];
      for (Int_t iDet = 5; iDet >= 0; iDet--) {// TOF, TRD, TPC, ITS clusters
          nsp += track->GetNcls(iDet);

          if (iDet==0) { // ITS "extra" clusters
             track->GetClusters(iDet,idx);
             for (Int_t i=6; i<12; i++) if(idx[i] >= 0) nsp++;
          }  
      }

      if (nsp) {
	AliTrackPointArray *sp = new AliTrackPointArray(nsp);
	track->SetTrackPointArray(sp);
	Int_t isptrack = 0;
	for (Int_t iDet = 5; iDet >= 0; iDet--) {
	  AliTracker *tracker = fTracker[iDet];
	  if (!tracker) continue;
	  Int_t nspdet = track->GetClusters(iDet,idx);

	  if (iDet==0) // ITS "extra" clusters             
             for (Int_t i=6; i<12; i++) if(idx[i] >= 0) nspdet++;

	  if (nspdet <= 0) continue;
	  AliTrackPoint p;
	  Int_t isp = 0;
	  Int_t isp2 = 0;
	  while (isp2 < nspdet) {
	    Bool_t isvalid=kTRUE;

            Int_t index=idx[isp++];
            if (index < 0) continue;

            TString dets = fgkDetectorName[iDet];
            if ((fUseTrackingErrorsForAlignment.CompareTo(dets) == 0) ||
            fUseTrackingErrorsForAlignment.BeginsWith(dets+" ") ||
            fUseTrackingErrorsForAlignment.EndsWith(" "+dets) ||
            fUseTrackingErrorsForAlignment.Contains(" "+dets+" ")) {
              isvalid = tracker->GetTrackPointTrackingError(index,p,track);
	    } else {
	      isvalid = tracker->GetTrackPoint(index,p); 
	    } 
	    isp2++;
	    if (!isvalid) continue;
	    if (iDet==0 && (isp-1)>=6) p.SetExtra();
	    sp->AddPoint(isptrack,&p); isptrack++;
	  }
	}	
      }
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

//_____________________________________________________________________________
void AliReconstruction::CheckQA()
{
// check the QA of SIM for this run and remove the detectors 
// with status Fatal
  
//	TString newRunLocalReconstruction ; 
//	TString newRunTracking ;
//	TString newFillESD ;
//	 
//	for (Int_t iDet = 0; iDet < AliQA::kNDET; iDet++) {
//		TString detName(AliQA::GetDetName(iDet)) ;
//		AliQA * qa = AliQA::Instance(AliQA::DETECTORINDEX_t(iDet)) ;       
//      if ( qa->IsSet(AliQA::DETECTORINDEX_t(iDet), AliQA::kSIM, specie, AliQA::kFATAL)) {
//        AliInfo(Form("QA status for %s %s in Hits and/or SDIGITS  and/or Digits was Fatal; No reconstruction performed", 
//                   detName.Data(), AliRecoParam::GetEventSpecieName(es))) ;
//			} else {
//			if ( fRunLocalReconstruction.Contains(AliQA::GetDetName(iDet)) || 
//					fRunLocalReconstruction.Contains("ALL") )  {
//				newRunLocalReconstruction += detName ; 
//				newRunLocalReconstruction += " " ; 			
//			}
//			if ( fRunTracking.Contains(AliQA::GetDetName(iDet)) || 
//					fRunTracking.Contains("ALL") )  {
//				newRunTracking += detName ; 
//				newRunTracking += " " ; 			
//			}
//			if ( fFillESD.Contains(AliQA::GetDetName(iDet)) || 
//					fFillESD.Contains("ALL") )  {
//				newFillESD += detName ; 
//				newFillESD += " " ; 			
//			}
//		}
//	}
//	fRunLocalReconstruction = newRunLocalReconstruction ; 
//	fRunTracking            = newRunTracking ; 
//	fFillESD                = newFillESD ; 
}

//_____________________________________________________________________________
Int_t AliReconstruction::GetDetIndex(const char* detector)
{
  // return the detector index corresponding to detector
  Int_t index = -1 ; 
  for (index = 0; index < kNDetectors ; index++) {
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
 //for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
 for (Int_t iDet = 0; iDet < 1; iDet++) { // for the time being only ITS
   //if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
   if(fTracker[iDet]) {
      AliPlaneEff *planeeff=fTracker[iDet]->GetPlaneEff();
      TString name=planeeff->GetName();
      name+=".root";
      TFile* pefile = TFile::Open(name, "RECREATE");
      ret=(Bool_t)planeeff->Write();
      pefile->Close();
      if(planeeff->GetCreateHistos()) {
        TString hname=planeeff->GetName();
        hname+="Histo.root";
        ret*=planeeff->WriteHistosToFile(hname,"RECREATE");
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

  gROOT->ProcessLine("if (!AliEveEventManager::GetMaster()){new AliEveEventManager();AliEveEventManager::GetMaster()->AddNewEventCommand(\"alieve_online_on_new_event()\");gEve->AddEvent(AliEveEventManager::GetMaster());};");
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
  gROOT->ProcessLine(Form("AliEveEventManager::GetMaster()->SetEvent((AliRunLoader*)0x%lx,(AliRawReader*)0x%lx,(AliESDEvent*)0x%lx,(AliESDfriend*)0x%lx);",fRunLoader,fRawReader,fesd,fesdf));
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

//_____________________________________________________________________________
Bool_t AliReconstruction::InitRecoParams() 
{
  // The method accesses OCDB and retrieves all
  // the available reco-param objects from there.

  Bool_t isOK = kTRUE;

  TString detStr = fLoadCDB;
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {

    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;

    if (fRecoParam.GetDetRecoParamArray(iDet)) {
      AliInfo(Form("Using custom reconstruction parameters for detector %s",fgkDetectorName[iDet]));
      continue;
    }

    AliInfo(Form("Loading reconstruction parameter objects for detector %s",fgkDetectorName[iDet]));
  
    AliCDBPath path(fgkDetectorName[iDet],"Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(!entry){ 
      AliWarning(Form("Couldn't find RecoParam entry in OCDB for detector %s",fgkDetectorName[iDet]));
      isOK = kFALSE;
    }
    else {
      TObject *recoParamObj = entry->GetObject();
      if (dynamic_cast<TObjArray*>(recoParamObj)) {
	// The detector has a normal TobjArray of AliDetectorRecoParam objects
	// Registering them in AliRecoParam
	fRecoParam.AddDetRecoParamArray(iDet,dynamic_cast<TObjArray*>(recoParamObj));
      }
      else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
	// The detector has only onse set of reco parameters
	// Registering it in AliRecoParam
	AliInfo(Form("Single set of reconstruction parameters found for detector %s",fgkDetectorName[iDet]));
	dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
	fRecoParam.AddDetRecoParam(iDet,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
      }
      else {
	AliError(Form("No valid RecoParam object found in the OCDB for detector %s",fgkDetectorName[iDet]));
	isOK = kFALSE;
      }
      entry->SetOwner(0);
      AliCDBManager::Instance()->UnloadFromCache(path.GetPath());
    }
  }

  if (AliDebugLevel() > 0) fRecoParam.Print();

  return isOK;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::GetEventInfo() 
{
  // Fill the event info object
  // ...
  AliCodeTimerAuto("")

  AliCentralTrigger *aCTP = NULL;
  if (fRawReader) {
    fEventInfo.SetEventType(fRawReader->GetType());

    ULong64_t mask = fRawReader->GetClassMask();
    fEventInfo.SetTriggerMask(mask);
    UInt_t clmask = fRawReader->GetDetectorPattern()[0];
    fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(clmask));

    aCTP = new AliCentralTrigger();
    TString configstr("");
    if (!aCTP->LoadConfiguration(configstr)) { // Load CTP config from OCDB
      AliError("No trigger configuration found in OCDB! The trigger configuration information will not be used!");
      delete aCTP;
      return kFALSE;
    }
    aCTP->SetClassMask(mask);
    aCTP->SetClusterMask(clmask);
  }
  else {
    fEventInfo.SetEventType(AliRawEventHeaderBase::kPhysicsEvent);

    if (fRunLoader && (!fRunLoader->LoadTrigger())) {
      aCTP = fRunLoader->GetTrigger();
      fEventInfo.SetTriggerMask(aCTP->GetClassMask());
      fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(aCTP->GetClusterMask()));
    }
    else {
      AliWarning("No trigger can be loaded! The trigger information will not be used!");
      return kFALSE;
    }
  }

  AliTriggerConfiguration *config = aCTP->GetConfiguration();
  if (!config) {
    AliError("No trigger configuration has been found! The trigger configuration information will not be used!");
    if (fRawReader) delete aCTP;
    return kFALSE;
  }

  UChar_t clustmask = 0;
  TString trclasses;
  ULong64_t trmask = fEventInfo.GetTriggerMask();
  const TObjArray& classesArray = config->GetClasses();
  Int_t nclasses = classesArray.GetEntriesFast();
  for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
    AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
    if (trclass) {
      Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
      fesd->SetTriggerClass(trclass->GetName(),trindex);
      if (fRawReader) fRawReader->LoadTriggerClass(trclass->GetName(),trindex);
      if (trmask & (1 << trindex)) {
	trclasses += " ";
	trclasses += trclass->GetName();
	trclasses += " ";
	clustmask |= trclass->GetCluster()->GetClusterMask();
      }
    }
  }
  fEventInfo.SetTriggerClasses(trclasses);

  // Set the information in ESD
  fesd->SetTriggerMask(trmask);
  fesd->SetTriggerCluster(clustmask);

  if (!aCTP->CheckTriggeredDetectors()) {
    if (fRawReader) delete aCTP;
    return kFALSE;
  }    

  if (fRawReader) delete aCTP;

  // We have to fill also the HLT decision here!!
  // ...

  return kTRUE;
}

const char *AliReconstruction::MatchDetectorList(const char *detectorList, UInt_t detectorMask)
{
  // Match the detector list found in the rec.C or the default 'ALL'
  // to the list found in the GRP (stored there by the shuttle PP which
  // gets the information from ECS)
  static TString resultList;
  TString detList = detectorList;

  resultList = "";

  for(Int_t iDet = 0; iDet < (AliDAQ::kNDetectors-1); iDet++) {
    if ((detectorMask >> iDet) & 0x1) {
      TString det = AliDAQ::OfflineModuleName(iDet);
      if ((detList.CompareTo("ALL") == 0) ||
	  ((detList.BeginsWith("ALL ") ||
	    detList.EndsWith(" ALL") ||
	    detList.Contains(" ALL ")) &&
	   !(detList.BeginsWith("-"+det+" ") ||
	     detList.EndsWith(" -"+det) ||
	     detList.Contains(" -"+det+" "))) ||
	  (detList.CompareTo(det) == 0) ||
	  detList.BeginsWith(det+" ") ||
	  detList.EndsWith(" "+det) ||
	  detList.Contains( " "+det+" " )) {
	if (!resultList.EndsWith(det + " ")) {
	  resultList += det;
	  resultList += " ";
	}
      }	       
    }
  }

  // HLT
  if ((detectorMask >> AliDAQ::kHLTId) & 0x1) {
    TString hltDet = AliDAQ::OfflineModuleName(AliDAQ::kNDetectors-1);
    if ((detList.CompareTo("ALL") == 0) ||
	((detList.BeginsWith("ALL ") ||
	  detList.EndsWith(" ALL") ||
	  detList.Contains(" ALL ")) &&
	 !(detList.BeginsWith("-"+hltDet+" ") ||
	   detList.EndsWith(" -"+hltDet) ||
	   detList.Contains(" -"+hltDet+" "))) ||
	(detList.CompareTo(hltDet) == 0) ||
	detList.BeginsWith(hltDet+" ") ||
	detList.EndsWith(" "+hltDet) ||
	detList.Contains( " "+hltDet+" " )) {
      resultList += hltDet;
    }
  }

  return resultList.Data();

}

//______________________________________________________________________________
void AliReconstruction::Abort(const char *method, EAbort what)
{
  // Abort processing. If what = kAbortProcess, the Process() loop will be
  // aborted. If what = kAbortFile, the current file in a chain will be
  // aborted and the processing will continue with the next file, if there
  // is no next file then Process() will be aborted. Abort() can also  be
  // called from Begin(), SlaveBegin(), Init() and Notify(). After abort
  // the SlaveTerminate() and Terminate() are always called. The abort flag
  // can be checked in these methods using GetAbort().
  //
  // The method is overwritten in AliReconstruction for better handling of
  // reco specific errors 

  if (!fStopOnError) return;

  CleanUp();

  TString whyMess = method;
  whyMess += " failed! Aborting...";

  AliError(whyMess.Data());

  fAbort = what;
  TString mess = "Abort";
  if (fAbort == kAbortProcess)
    mess = "AbortProcess";
  else if (fAbort == kAbortFile)
    mess = "AbortFile";

  Info(mess, whyMess.Data());
}

