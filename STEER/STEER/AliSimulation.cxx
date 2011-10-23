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
// class for running generation, simulation and digitization                 //
//                                                                           //
// Hits, sdigits and digits are created for all detectors by typing:         //
//                                                                           //
//   AliSimulation sim;                                                      //
//   sim.Run();                                                              //
//                                                                           //
// The Run method returns kTRUE in case of successful execution.             //
// The number of events can be given as argument to the Run method or it     //
// can be set by                                                             //
//                                                                           //
//   sim.SetNumberOfEvents(n);                                               //
//                                                                           //
// The name of the configuration file can be passed as argument to the       //
// AliSimulation constructor or can be specified by                          //
//                                                                           //
//   sim.SetConfigFile("...");                                               //
//                                                                           //
// The generation of particles and the simulation of detector hits can be    //
// switched on or off by                                                     //
//                                                                           //
//   sim.SetRunGeneration(kTRUE);   // generation of primary particles       //
//   sim.SetRunSimulation(kFALSE);  // but no tracking                       //
//                                                                           //
// For which detectors sdigits and digits will be created, can be steered    //
// by                                                                        //
//                                                                           //
//   sim.SetMakeSDigits("ALL");     // make sdigits for all detectors        //
//   sim.SetMakeDigits("ITS TPC");  // make digits only for ITS and TPC      //
//                                                                           //
// The argument is a (case sensitive) string with the names of the           //
// detectors separated by a space. An empty string ("") can be used to       //
// disable the creation of sdigits or digits. The special string "ALL"       //
// selects all available detectors. This is the default.                     //
//                                                                           //
// The creation of digits from hits instead of from sdigits can be selected  //
// by                                                                        //
//                                                                           //
//   sim.SetMakeDigitsFromHits("TRD");                                       //
//                                                                           //
// The argument is again a string with the selected detectors. Be aware that //
// this feature is not available for all detectors and that merging is not   //
// possible, when digits are created directly from hits.                     //
//                                                                           //
// Background events can be merged by calling                                //
//                                                                           //
//   sim.MergeWith("background/galice.root", 2);                             //
//                                                                           //
// The first argument is the file name of the background galice file. The    //
// second argument is the number of signal events per background event.      //
// By default this number is calculated from the number of available         //
// background events. MergeWith can be called several times to merge more    //
// than two event streams. It is assumed that the sdigits were already       //
// produced for the background events.                                       //
//                                                                           //
// The output of raw data can be switched on by calling                      //
//                                                                           //
//   sim.SetWriteRawData("MUON");   // write raw data for MUON               //
//                                                                           //
// The default output format of the raw data are DDL files. They are         //
// converted to a DATE file, if a file name is given as second argument.     //
// For this conversion the program "dateStream" is required. If the file     //
// name has the extension ".root", the DATE file is converted to a root      //
// file. The program "alimdc" is used for this purpose. For the conversion   //
// to DATE and root format the two conversion programs have to be installed. //
// Only the raw data in the final format is kept if the third argument is    //
// kTRUE.                                                                    //
//                                                                           //
// The methods RunSimulation, RunSDigitization, RunDigitization,             //
// RunHitsDigitization and WriteRawData can be used to run only parts of     //
// the full simulation chain. The creation of raw data DDL files and their   //
// conversion to the DATE or root format can be run directly by calling      //
// the methods WriteRawFiles, ConvertRawFilesToDate and ConvertDateToRoot.   //
//                                                                           //
// The default number of events per file, which is usually set in the        //
// config file, can be changed for individual detectors and data types       //
// by calling                                                                //
//                                                                           //
//   sim.SetEventsPerFile("PHOS", "Reconstructed Points", 3);                //
//                                                                           //
// The first argument is the detector, the second one the data type and the  //
// last one the number of events per file. Valid data types are "Hits",      //
// "Summable Digits", "Digits", "Reconstructed Points" and "Tracks".         //
// The number of events per file has to be set before the simulation of      //
// hits. Otherwise it has no effect.                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TCint.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TVirtualMCApplication.h>
#include <TDatime.h>

#include "AliAlignObj.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliCDBStorage.h"
#include "AliCTPRawData.h"
#include "AliCentralTrigger.h"
#include "AliCentralTrigger.h"
#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliDigitizer.h"
#include "AliESDEvent.h"
#include "AliGRPObject.h"
#include "AliGenEventHeader.h"
#include "AliGenerator.h"
#include "AliGeomManager.h"
#include "AliHLTSimulation.h"
#include "AliHeader.h"
#include "AliLego.h"
#include "AliLegoGenerator.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliModule.h"
#include "AliPDG.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"
#include "AliRun.h"
#include "AliDigitizationInput.h"
#include "AliRunLoader.h"
#include "AliSimulation.h"
#include "AliSysInfo.h"
#include "AliVertexGenFile.h"

ClassImp(AliSimulation)

AliSimulation *AliSimulation::fgInstance = 0;
const char* AliSimulation::fgkDetectorName[AliSimulation::fgkNDetectors] = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "ACORDE"
// #ifdef MFT_UPGRADE
//                                                                             ,"MFT"
// #endif 
                                                                            ,"MFT"    // AU
									    ,"HLT"
};

//_____________________________________________________________________________
AliSimulation::AliSimulation(const char* configFileName,
			     const char* name, const char* title) :
  TNamed(name, title),

  fRunGeneration(kTRUE),
  fRunSimulation(kTRUE),
  fLoadAlignFromCDB(kTRUE),
  fLoadAlObjsListOfDets("ALL"),
  fMakeSDigits("ALL"),
  fMakeDigits("ALL"),
  fTriggerConfig(""),
  fMakeDigitsFromHits(""),
  fWriteRawData(""),
  fRawDataFileName(""),
  fDeleteIntermediateFiles(kFALSE),
  fWriteSelRawData(kFALSE),
  fStopOnError(kFALSE),
  fNEvents(1),
  fConfigFileName(configFileName),
  fGAliceFileName("galice.root"),
  fEventsPerFile(),
  fBkgrdFileNames(NULL),
  fAlignObjArray(NULL),
  fUseBkgrdVertex(kTRUE),
  fRegionOfInterest(kFALSE),
  fCDBUri(""),
  fQARefUri(""), 
  fSpecCDBUri(),
  fRun(-1),
  fSeed(0),
  fInitCDBCalled(kFALSE),
  fInitRunNumberCalled(kFALSE),
  fSetRunNumberFromDataCalled(kFALSE),
  fEmbeddingFlag(kFALSE),
  fLego(NULL),
  fKey(0),
  fUseVertexFromCDB(0),
  fUseMagFieldFromGRP(0),
  fGRPWriteLocation(Form("local://%s", gSystem->pwd())),
  fUseTimeStampFromCDB(0),
  fTimeStart(0),
  fTimeEnd(0),
  fQADetectors("ALL"),                  
  fQATasks("ALL"),	
  fRunQA(kTRUE), 
  fEventSpecie(AliRecoParam::kDefault),
  fWriteQAExpertData(kTRUE), 
  fGeometryFile(),
  fRunHLT("default"),
  fpHLT(NULL),
  fWriteGRPEntry(kTRUE)
{
// create simulation object with default parameters
  fgInstance = this;
  SetGAliceFile("galice.root");
  
// for QA
	AliQAManager * qam = AliQAManager::QAManager(AliQAv1::kSIMMODE) ; 
	qam->SetActiveDetectors(fQADetectors) ; 
	fQATasks = Form("%d %d %d", AliQAv1::kHITS, AliQAv1::kSDIGITS, AliQAv1::kDIGITS) ; 
	qam->SetTasks(fQATasks) ; 	
}

//_____________________________________________________________________________
AliSimulation::~AliSimulation()
{
// clean up

  fEventsPerFile.Delete();
//  if(fAlignObjArray) fAlignObjArray->Delete(); // fAlignObjArray->RemoveAll() ???
//  delete fAlignObjArray; fAlignObjArray=0;

  if (fBkgrdFileNames) {
    fBkgrdFileNames->Delete();
    delete fBkgrdFileNames;
  }

  fSpecCDBUri.Delete();
  if (fgInstance==this) fgInstance = 0;

  AliQAManager::QAManager()->ShowQA() ; 
  AliQAManager::Destroy() ; 	
  AliCodeTimer::Instance()->Print();
}


//_____________________________________________________________________________
void AliSimulation::SetNumberOfEvents(Int_t nEvents)
{
// set the number of events for one run

  fNEvents = nEvents;
}

//_____________________________________________________________________________
void AliSimulation::InitQA()
{
  // activate a default CDB storage
  // First check if we have any CDB storage set, because it is used 
  // to retrieve the calibration and alignment constants
  
  if (fInitCDBCalled) return;
  fInitCDBCalled = kTRUE;

  AliQAManager * qam = AliQAManager::QAManager(AliQAv1::kSIMMODE) ; 
  qam->SetActiveDetectors(fQADetectors) ; 
  fQATasks = Form("%d %d %d", AliQAv1::kHITS, AliQAv1::kSDIGITS, AliQAv1::kDIGITS) ; 
  qam->SetTasks(fQATasks) ;
 	if (fWriteQAExpertData)
    qam->SetWriteExpert() ; 
  
  if (qam->IsDefaultStorageSet()) {
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("Default QA reference storage has been already set !");
    AliWarning(Form("Ignoring the default storage declared in AliSimulation: %s",fQARefUri.Data()));
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    fQARefUri = qam->GetDefaultStorage()->GetURI();
  } else {
      if (fQARefUri.Length() > 0) {
        AliDebug(2,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        AliDebug(2, Form("Default QA reference storage is set to: %s", fQARefUri.Data()));
        AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      } else {
        fQARefUri="local://$ALICE_ROOT/QARef";
        AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        AliWarning("Default QA reference storage not yet set !!!!");
        AliWarning(Form("Setting it now to: %s", fQARefUri.Data()));
        AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      }
    qam->SetDefaultStorage(fQARefUri);
  }
}

//_____________________________________________________________________________
void AliSimulation::InitCDB()
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
    AliWarning(Form("Ignoring the default storage declared in AliSimulation: %s",fCDBUri.Data()));
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
    	fCDBUri="local://$ALICE_ROOT/OCDB";
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
void AliSimulation::InitRunNumber(){
// check run number. If not set, set it to 0 !!!!
  
  if (fInitRunNumberCalled) return;
  fInitRunNumberCalled = kTRUE;
  
  AliCDBManager* man = AliCDBManager::Instance();
  if (man->GetRun() >= 0)
  {
    	AliFatal(Form("Run number cannot be set in AliCDBManager before start of simulation: "
			"Use external variable DC_RUN or AliSimulation::SetRun()!"));
  }
    
  if(fRun >= 0) {
    	AliDebug(2, Form("Setting CDB run number to: %d",fRun));
  } else {
    	fRun=0;
    	AliWarning(Form("Run number not yet set !!!! Setting it now to: %d",
			fRun));
  }
  man->SetRun(fRun);

  man->Print();

}

//_____________________________________________________________________________
void AliSimulation::SetCDBLock() {
  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  
  ULong_t key = AliCDBManager::Instance()->SetLock(1);
  if (key) fKey = key;
}

//_____________________________________________________________________________
void AliSimulation::SetDefaultStorage(const char* uri) {
  // Store the desired default CDB storage location
  // Activate it later within the Run() method
  
  fCDBUri = uri;
  
}

//_____________________________________________________________________________
void AliSimulation::SetQARefDefaultStorage(const char* uri) {
  // Store the desired default CDB storage location
  // Activate it later within the Run() method
  
  fQARefUri = uri;
  AliQAv1::SetQARefStorage(fQARefUri.Data()) ;
}

//_____________________________________________________________________________
void AliSimulation::SetSpecificStorage(const char* calibType, const char* uri) {
// Store a detector-specific CDB storage location
// Activate it later within the Run() method

  AliCDBPath aPath(calibType);
  if(!aPath.IsValid()){
  	AliError(Form("Not a valid path: %s", calibType));
  	return;
  }

  TObject* obj = fSpecCDBUri.FindObject(calibType);
  if (obj) fSpecCDBUri.Remove(obj);
  fSpecCDBUri.Add(new TNamed(calibType, uri));

}

//_____________________________________________________________________________
void AliSimulation::SetRunNumber(Int_t run)
{
// sets run number
// Activate it later within the Run() method

	fRun = run;
}

//_____________________________________________________________________________
void AliSimulation::SetSeed(Int_t seed)
{
// sets seed number
// Activate it later within the Run() method

	fSeed = seed;
}

//_____________________________________________________________________________
Bool_t AliSimulation::SetRunNumberFromData()
{
  // Set the CDB manager run number
  // The run number is retrieved from gAlice

    if (fSetRunNumberFromDataCalled) return kTRUE;
    fSetRunNumberFromDataCalled = kTRUE;    
  
    AliCDBManager* man = AliCDBManager::Instance();
    Int_t runData = -1, runCDB = -1;
  
    AliRunLoader* runLoader = LoadRun("READ");
    if (!runLoader) return kFALSE;
    else {
    	runData = runLoader->GetHeader()->GetRun();
	delete runLoader;
    }
  
    runCDB = man->GetRun();
    if(runCDB >= 0) {
	if (runCDB != runData) {
    		AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    		AliWarning(Form("A run number was previously set in AliCDBManager: %d !", runCDB));
    		AliWarning(Form("It will be replaced with the run number got from run header: %d !", runData));
    		AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");	
	}
   	
    }
      
    man->SetRun(runData);
    fRun = runData;
    
    if(man->GetRun() < 0) {
    	AliError("Run number not properly initalized!");
	return kFALSE;
    }
  
    man->Print();
    
    return kTRUE;
}

//_____________________________________________________________________________
void AliSimulation::SetConfigFile(const char* fileName)
{
// set the name of the config file

  fConfigFileName = fileName;
}

//_____________________________________________________________________________
void AliSimulation::SetGAliceFile(const char* fileName)
{
// set the name of the galice file
// the path is converted to an absolute one if it is relative

  fGAliceFileName = fileName;
  if (!gSystem->IsAbsoluteFileName(fGAliceFileName)) {
    char* absFileName = gSystem->ConcatFileName(gSystem->WorkingDirectory(),
						fGAliceFileName);
    fGAliceFileName = absFileName;
    delete[] absFileName;
  }

  AliDebug(2, Form("galice file name set to %s", fileName));
}

//_____________________________________________________________________________
void AliSimulation::SetEventsPerFile(const char* detector, const char* type, 
				     Int_t nEvents)
{
// set the number of events per file for the given detector and data type
// ("Hits", "Summable Digits", "Digits", "Reconstructed Points" or "Tracks")

  TNamed* obj = new TNamed(detector, type);
  obj->SetUniqueID(nEvents);
  fEventsPerFile.Add(obj);
}

//_____________________________________________________________________________
Bool_t AliSimulation::MisalignGeometry(AliRunLoader *runLoader)
{
  // Read the alignment objects from CDB.
  // Each detector is supposed to have the
  // alignment objects in DET/Align/Data CDB path.
  // All the detector objects are then collected,
  // sorted by geometry level (starting from ALIC) and
  // then applied to the TGeo geometry.
  // Finally an overlaps check is performed.

  if (!AliGeomManager::GetGeometry() || !AliGeomManager::GetGeometry()->IsClosed()) {
    AliError("Can't apply the misalignment! Geometry is not loaded or it is still opened!");
    return kFALSE;
  }  
  
  // initialize CDB storage, run number, set CDB lock
  InitCDB();
//  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  SetCDBLock();
    
  Bool_t delRunLoader = kFALSE;
  if (!runLoader) {
    runLoader = LoadRun("READ");
    if (!runLoader) return kFALSE;
    delRunLoader = kTRUE;
  }
  
  // Export ideal geometry 
  if(!IsGeometryFromFile()) AliGeomManager::GetGeometry()->Export("geometry.root");

  // Load alignment data from CDB and apply to geometry through AliGeomManager
  if(fLoadAlignFromCDB){
    
    TString detStr = fLoadAlObjsListOfDets;
    TString loadAlObjsListOfDets = "";
    
    TObjArray* detArray = runLoader->GetAliRun()->Detectors();
    for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
      AliModule* det = (AliModule*) detArray->At(iDet);
      if (!det || !det->IsActive()) continue;
      if (IsSelected(det->GetName(), detStr)) {
        //add det to list of dets to be aligned from CDB
        loadAlObjsListOfDets += det->GetName();
        loadAlObjsListOfDets += " ";
      }
    } // end loop over detectors
    loadAlObjsListOfDets.Prepend("GRP "); //add alignment objects for non-sensitive modules
    AliGeomManager::ApplyAlignObjsFromCDB(loadAlObjsListOfDets.Data());
  }else{
    // Check if the array with alignment objects was
    // provided by the user. If yes, apply the objects
    // to the present TGeo geometry
    if (fAlignObjArray) {
      if (AliGeomManager::ApplyAlignObjsToGeom(*fAlignObjArray) == kFALSE) {
        AliError("The misalignment of one or more volumes failed!"
                 "Compare the list of simulated detectors and the list of detector alignment data!");
        if (delRunLoader) delete runLoader;
        return kFALSE;
      }
    }
  }

  // Update the internal geometry of modules (ITS needs it)
  TString detStr = fLoadAlObjsListOfDets;
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {

    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      det->UpdateInternalGeometry();
    }
  } // end loop over detectors


  if (delRunLoader) delete runLoader;

  return kTRUE;
}

//_____________________________________________________________________________
void AliSimulation::MergeWith(const char* fileName, Int_t nSignalPerBkgrd)
{
// add a file with background events for merging

  TObjString* fileNameStr = new TObjString(fileName);
  fileNameStr->SetUniqueID(nSignalPerBkgrd);
  if (!fBkgrdFileNames) fBkgrdFileNames = new TObjArray;
  fBkgrdFileNames->Add(fileNameStr);
}

void AliSimulation::EmbedInto(const char* fileName, Int_t nSignalPerBkgrd)
{
// add a file with background events for embeddin
  MergeWith(fileName, nSignalPerBkgrd);
  fEmbeddingFlag = kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::Run(Int_t nEvents)
{
// run the generation, simulation and digitization

 
  AliCodeTimerAuto("",0)
  AliSysInfo::AddStamp("Start_Run");
  
  // Load run number and seed from environmental vars
  ProcessEnvironmentVars();
  AliSysInfo::AddStamp("ProcessEnvironmentVars");

  gRandom->SetSeed(fSeed);
   
  if (nEvents > 0) fNEvents = nEvents;

  // create and setup the HLT instance
  if (!fRunHLT.IsNull() && !CreateHLT()) {
    if (fStopOnError) return kFALSE;
    // disable HLT
    fRunHLT="";
  }
  
  // generation and simulation -> hits
  if (fRunGeneration) {
    if (!RunSimulation()) if (fStopOnError) return kFALSE;
  }
  AliSysInfo::AddStamp("RunSimulation");
           
  // initialize CDB storage from external environment
  // (either CDB manager or AliSimulation setters),
  // if not already done in RunSimulation()
  InitCDB();
  AliSysInfo::AddStamp("InitCDB");
  
  // Set run number in CDBManager from data 
  // From this point on the run number must be always loaded from data!
  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  
  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  SetCDBLock();

  // If RunSimulation was not called, load the geometry and misalign it
  if (!AliGeomManager::GetGeometry()) {
    // Initialize the geometry manager
    AliGeomManager::LoadGeometry("geometry.root");
    AliSysInfo::AddStamp("GetGeometry");
//    // Check that the consistency of symbolic names for the activated subdetectors
//    // in the geometry loaded by AliGeomManager
//    AliRunLoader* runLoader = LoadRun("READ");
//    if (!runLoader) return kFALSE;
//
//    TString detsToBeChecked = "";
//    TObjArray* detArray = runLoader->GetAliRun()->Detectors();
//    for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
//      AliModule* det = (AliModule*) detArray->At(iDet);
//      if (!det || !det->IsActive()) continue;
//      detsToBeChecked += det->GetName();
//      detsToBeChecked += " ";
//    } // end loop over detectors
//    if(!AliGeomManager::CheckSymNamesLUT(detsToBeChecked.Data()))
    if(!AliGeomManager::CheckSymNamesLUT("ALL"))
	AliFatalClass("Current loaded geometry differs in the definition of symbolic names!");
	
    if (!AliGeomManager::GetGeometry()) if (fStopOnError) return kFALSE;
    // Misalign geometry
    if(!MisalignGeometry()) if (fStopOnError) return kFALSE;
  }
  AliSysInfo::AddStamp("MissalignGeometry");


  // hits -> summable digits
  AliSysInfo::AddStamp("Start_sdigitization");
  if (!fMakeSDigits.IsNull()) {
    if (!RunSDigitization(fMakeSDigits)) if (fStopOnError) return kFALSE;
 
  }
  AliSysInfo::AddStamp("Stop_sdigitization");
  
  AliSysInfo::AddStamp("Start_digitization");  
  // summable digits -> digits  
  if (!fMakeDigits.IsNull()) {
    if (!RunDigitization(fMakeDigits, fMakeDigitsFromHits)) {
      if (fStopOnError) return kFALSE;
    }
   }
  AliSysInfo::AddStamp("Stop_digitization");

  
  
  // hits -> digits
  if (!fMakeDigitsFromHits.IsNull()) {
    if (fBkgrdFileNames && (fBkgrdFileNames->GetEntriesFast() > 0)) {
      AliWarning(Form("Merging and direct creation of digits from hits " 
                 "was selected for some detectors. "
                 "No merging will be done for the following detectors: %s",
                 fMakeDigitsFromHits.Data()));
    }
    if (!RunHitsDigitization(fMakeDigitsFromHits)) {
      if (fStopOnError) return kFALSE;
    }
  }

  AliSysInfo::AddStamp("Hits2Digits");
  
  
  // digits -> trigger
  if (!fTriggerConfig.IsNull() && !RunTrigger(fTriggerConfig,fMakeDigits)) {
    if (fStopOnError) return kFALSE;
  }

  AliSysInfo::AddStamp("RunTrigger");
  
  
  // digits -> raw data
  if (!fWriteRawData.IsNull()) {
    if (!WriteRawData(fWriteRawData, fRawDataFileName, 
		      fDeleteIntermediateFiles,fWriteSelRawData)) {
      if (fStopOnError) return kFALSE;
    }
  }

  AliSysInfo::AddStamp("WriteRaw");
  
  // run HLT simulation on simulated digit data if raw data is not
  // simulated, otherwise its called as part of WriteRawData
  if (!fRunHLT.IsNull() && fWriteRawData.IsNull()) {
    if (!RunHLT()) {
      if (fStopOnError) return kFALSE;
    }
  }

  AliSysInfo::AddStamp("RunHLT");
  
  //QA
	if (fRunQA) {
		Bool_t rv = RunQA() ; 
		if (!rv)
			if (fStopOnError) 
				return kFALSE ;   	
	}

  AliSysInfo::AddStamp("RunQA");

  // Cleanup of CDB manager: cache and active storages!
  AliCDBManager::Instance()->ClearCache();

  return kTRUE;
}

//_______________________________________________________________________
Bool_t AliSimulation::RunLego(const char *setup, Int_t nc1, Float_t c1min,
		     Float_t c1max,Int_t nc2,Float_t c2min,Float_t c2max,
		     Float_t rmin,Float_t rmax,Float_t zmax, AliLegoGenerator* gener, Int_t nev)
{
  //
  // Generates lego plots of:
  //    - radiation length map phi vs theta
  //    - radiation length map phi vs eta
  //    - interaction length map
  //    - g/cm2 length map
  //
  //  ntheta    bins in theta, eta
  //  themin    minimum angle in theta (degrees)
  //  themax    maximum angle in theta (degrees)
  //  nphi      bins in phi
  //  phimin    minimum angle in phi (degrees)
  //  phimax    maximum angle in phi (degrees)
  //  rmin      minimum radius
  //  rmax      maximum radius
  //  
  //
  //  The number of events generated = ntheta*nphi
  //  run input parameters in macro setup (default="Config.C")
  //
  //  Use macro "lego.C" to visualize the 3 lego plots in spherical coordinates
  //Begin_Html
  /*
    <img src="picts/AliRunLego1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliRunLego2.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliRunLego3.gif">
  */
  //End_Html
  //

// run the generation and simulation

  AliCodeTimerAuto("",0)

  // initialize CDB storage and run number from external environment
  // (either CDB manager or AliSimulation setters)
  InitCDB();
  InitRunNumber();
  SetCDBLock();
  
  if (!gAlice) {
    AliError("no gAlice object. Restart aliroot and try again.");
    return kFALSE;
  }
  if (gAlice->Modules()->GetEntries() > 0) {
    AliError("gAlice was already run. Restart aliroot and try again.");
    return kFALSE;
  }

  AliInfo(Form("initializing gAlice with config file %s",
          fConfigFileName.Data()));

  // Number of events 
    if (nev == -1) nev  = nc1 * nc2;
    
  // check if initialisation has been done
  // If runloader has been initialized, set the number of events per file to nc1 * nc2
    
  // Set new generator
  if (!gener) gener  = new AliLegoGenerator();
  //
  // Configure Generator

  gener->SetRadiusRange(rmin, rmax);
  gener->SetZMax(zmax);
  gener->SetCoor1Range(nc1, c1min, c1max);
  gener->SetCoor2Range(nc2, c2min, c2max);
  
  
  //Create Lego object  
  fLego = new AliLego("lego",gener);

  //__________________________________________________________________________

  gAlice->Announce();

  gROOT->LoadMacro(setup);
  gInterpreter->ProcessLine(gAlice->GetConfigFunction());

  if(AliCDBManager::Instance()->GetRun() >= 0) { 
    SetRunNumber(AliCDBManager::Instance()->GetRun());
  } else {
    AliWarning("Run number not initialized!!");
  }

  AliRunLoader::Instance()->CdGAFile();
  
  AliPDG::AddParticlesToPdgDataBase();  
  
  gMC->SetMagField(TGeoGlobalMagField::Instance()->GetField());

  gAlice->GetMCApp()->Init();
  
  
  //Must be here because some MCs (G4) adds detectors here and not in Config.C
  gAlice->InitLoaders();
  AliRunLoader::Instance()->MakeTree("E");
  
  //
  // Save stuff at the beginning of the file to avoid file corruption
  AliRunLoader::Instance()->CdGAFile();
  gAlice->Write();

  //Save current generator
  AliGenerator *gen=gAlice->GetMCApp()->Generator();
  gAlice->GetMCApp()->ResetGenerator(gener);
  //Prepare MC for Lego Run
  gMC->InitLego();
  
  //Run Lego Object
  
  
  AliRunLoader::Instance()->SetNumberOfEventsPerFile(nev);
  gMC->ProcessRun(nev);
  
  // End of this run, close files
  FinishRun();
  // Restore current generator
  gAlice->GetMCApp()->ResetGenerator(gen);
  // Delete Lego Object
  delete fLego;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunTrigger(const char* config, const char* detectors)
{
  // run the trigger

  AliCodeTimerAuto("",0)

  // initialize CDB storage from external environment
  // (either CDB manager or AliSimulation setters),
  // if not already done in RunSimulation()
  InitCDB();
  
  // Set run number in CDBManager from data 
  // From this point on the run number must be always loaded from data!
  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  
  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  SetCDBLock();
   
   AliRunLoader* runLoader = LoadRun("READ");
   if (!runLoader) return kFALSE;
   TString trconfiguration = config;

   if (trconfiguration.IsNull()) {
     if(!fTriggerConfig.IsNull()) {
       trconfiguration = fTriggerConfig;
     }
     else
       AliWarning("No trigger descriptor is specified. Loading the one that is in the CDB.");
   }

   runLoader->MakeTree( "GG" );
   AliCentralTrigger* aCTP = runLoader->GetTrigger();
   // Load Configuration
   if (!aCTP->LoadConfiguration( trconfiguration ))
     return kFALSE;

   // digits -> trigger
   if( !aCTP->RunTrigger( runLoader , detectors ) ) {
      if (fStopOnError) {
	//  delete aCTP;
	return kFALSE;
      }
   }

   delete runLoader;

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::WriteTriggerRawData()
{
  // Writes the CTP (trigger) DDL raw data
  // Details of the format are given in the
  // trigger TDR - pages 134 and 135.
  AliCTPRawData writer;
  writer.RawData();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunSimulation(Int_t nEvents)
{
// run the generation and simulation

  AliCodeTimerAuto("",0)

  // initialize CDB storage and run number from external environment
  // (either CDB manager or AliSimulation setters)
  AliSysInfo::AddStamp("RunSimulation_Begin");
  InitCDB();
  AliSysInfo::AddStamp("RunSimulation_InitCDB");
  InitRunNumber();
  SetCDBLock();
  AliSysInfo::AddStamp("RunSimulation_SetCDBLock");
  
  if (!gAlice) {
    AliError("no gAlice object. Restart aliroot and try again.");
    return kFALSE;
  }
  if (gAlice->Modules()->GetEntries() > 0) {
    AliError("gAlice was already run. Restart aliroot and try again.");
    return kFALSE;
  }

  AliInfo(Form("initializing gAlice with config file %s",
          fConfigFileName.Data()));

  //
  // Initialize ALICE Simulation run
  //
  gAlice->Announce();

  //
  // If requested set the mag. field from the GRP entry.
  // After this the field is loccked and cannot be changed by Config.C
  if (fUseMagFieldFromGRP) {
      AliGRPManager grpM;
      grpM.ReadGRPEntry();
      grpM.SetMagField();
      AliInfo("Field is locked now. It cannot be changed in Config.C");
  }
//
// Execute Config.C
  gROOT->LoadMacro(fConfigFileName.Data());
  gInterpreter->ProcessLine(gAlice->GetConfigFunction());
  AliSysInfo::AddStamp("RunSimulation_Config");

//
// If requested obtain the vertex position and vertex sigma_z from the CDB
// This overwrites the settings from the Config.C  
  if (fUseVertexFromCDB) {
      Double_t vtxPos[3] = {0., 0., 0.}; 
      Double_t vtxSig[3] = {0., 0., 0.};
      AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertex");
      if (entry) {
	  AliESDVertex* vertex = dynamic_cast<AliESDVertex*> (entry->GetObject());
	  if (vertex) {
	      if(vertex->GetXRes()>2.8) { // > pipe radius --> it's a dummy object, don't use it 
		  entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexSPD");
		  if (entry) vertex = dynamic_cast<AliESDVertex*> (entry->GetObject());
	      }
	  }
	  if (vertex) {
	      vertex->GetXYZ(vtxPos);
	      vertex->GetSigmaXYZ(vtxSig);
	      AliInfo("Overwriting Config.C vertex settings !");
	      AliInfo(Form("Vertex position from OCDB entry: x = %13.3f, y = %13.3f, z = %13.3f (sigma = %13.3f)\n",
			   vtxPos[0], vtxPos[1], vtxPos[2], vtxSig[2]));
	      
	      AliGenerator *gen = gAlice->GetMCApp()->Generator();
	      gen->SetOrigin(vtxPos[0], vtxPos[1], vtxPos[2]);   // vertex position
	      gen->SetSigmaZ(vtxSig[2]);
	  }
      }
  }

  // If requested we take the SOR and EOR time-stamps from the GRP and use them
  // in order to generate the event time-stamps
  if (fUseTimeStampFromCDB) {
    AliGRPManager grpM;
    grpM.ReadGRPEntry();
    const AliGRPObject *grpObj = grpM.GetGRPData();
    if (!grpObj || (grpObj->GetTimeEnd() <= grpObj->GetTimeStart())) {
      AliError("Missing GRP or bad SOR/EOR time-stamps! Switching off the time-stamp generation from GRP!");
      fTimeStart = fTimeEnd = 0;
      fUseTimeStampFromCDB = kFALSE;
    }
    else {
      fTimeStart = grpObj->GetTimeStart();
      fTimeEnd = grpObj->GetTimeEnd();
    }
  }
  
  if(AliCDBManager::Instance()->GetRun() >= 0) { 
    AliRunLoader::Instance()->SetRunNumber(AliCDBManager::Instance()->GetRun());
    AliRunLoader::Instance()->SetNumberOfEventsPerRun(fNEvents);
  } else {
  	AliWarning("Run number not initialized!!");
  }
  
   AliRunLoader::Instance()->CdGAFile();
   

   AliPDG::AddParticlesToPdgDataBase();  

   gMC->SetMagField(TGeoGlobalMagField::Instance()->GetField());
   AliSysInfo::AddStamp("RunSimulation_GetField");
   
   gAlice->GetMCApp()->Init();
   AliSysInfo::AddStamp("RunSimulation_InitMCApp");

   //Must be here because some MCs (G4) adds detectors here and not in Config.C
   gAlice->InitLoaders();
   AliRunLoader::Instance()->MakeTree("E");
   AliRunLoader::Instance()->LoadKinematics("RECREATE");
   AliRunLoader::Instance()->LoadTrackRefs("RECREATE");
   AliRunLoader::Instance()->LoadHits("all","RECREATE");
   //
   // Save stuff at the beginning of the file to avoid file corruption
   AliRunLoader::Instance()->CdGAFile();
   gAlice->Write();
   gAlice->SetEventNrInRun(-1); //important - we start Begin event from increasing current number in run
   AliSysInfo::AddStamp("RunSimulation_InitLoaders");
  //___________________________________________________________________________________________
  
  AliSysInfo::AddStamp("RunSimulation_TriggerDescriptor");

  // Set run number in CDBManager
  AliInfo(Form("Run number: %d",AliCDBManager::Instance()->GetRun()));

  AliRunLoader* runLoader = AliRunLoader::Instance();
  if (!runLoader) {
             AliError(Form("gAlice has no run loader object. "
        		     "Check your config file: %s", fConfigFileName.Data()));
             return kFALSE;
  }
  SetGAliceFile(runLoader->GetFileName());
      
  // Misalign geometry
#if ROOT_VERSION_CODE < 331527
  AliGeomManager::SetGeometry(gGeoManager);
  
  // Check that the consistency of symbolic names for the activated subdetectors
  // in the geometry loaded by AliGeomManager
  TString detsToBeChecked = "";
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    detsToBeChecked += det->GetName();
    detsToBeChecked += " ";
  } // end loop over detectors
  if(!AliGeomManager::CheckSymNamesLUT(detsToBeChecked.Data()))
    AliFatalClass("Current loaded geometry differs in the definition of symbolic names!");
  MisalignGeometry(runLoader);
  AliSysInfo::AddStamp("RunSimulation_MisalignGeometry");
#endif

//   AliRunLoader* runLoader = AliRunLoader::Instance();
//   if (!runLoader) {
//     AliError(Form("gAlice has no run loader object. "
//                   "Check your config file: %s", fConfigFileName.Data()));
//     return kFALSE;
//   }
//   SetGAliceFile(runLoader->GetFileName());

  if (!gAlice->GetMCApp()->Generator()) {
    AliError(Form("gAlice has no generator object. "
                  "Check your config file: %s", fConfigFileName.Data()));
    return kFALSE;
  }

  // Write GRP entry corresponding to the setting found in Cofig.C
  if (fWriteGRPEntry)
    WriteGRPEntry();
  AliSysInfo::AddStamp("RunSimulation_WriteGRP");

  if (nEvents <= 0) nEvents = fNEvents;

  // get vertex from background file in case of merging
  if (fUseBkgrdVertex &&
      fBkgrdFileNames && (fBkgrdFileNames->GetEntriesFast() > 0)) {
    Int_t signalPerBkgrd = GetNSignalPerBkgrd(nEvents);
    const char* fileName = ((TObjString*)
			    (fBkgrdFileNames->At(0)))->GetName();
    AliInfo(Form("The vertex will be taken from the background "
                 "file %s with nSignalPerBackground = %d", 
                 fileName, signalPerBkgrd));
    AliVertexGenFile* vtxGen = new AliVertexGenFile(fileName, signalPerBkgrd);
    gAlice->GetMCApp()->Generator()->SetVertexGenerator(vtxGen);
  }

  if (!fRunSimulation) {
    gAlice->GetMCApp()->Generator()->SetTrackingFlag(0);
  }

  // set the number of events per file for given detectors and data types
  for (Int_t i = 0; i < fEventsPerFile.GetEntriesFast(); i++) {
    if (!fEventsPerFile[i]) continue;
    const char* detName = fEventsPerFile[i]->GetName();
    const char* typeName = fEventsPerFile[i]->GetTitle();
    TString loaderName(detName);
    loaderName += "Loader";
    AliLoader* loader = runLoader->GetLoader(loaderName);
    if (!loader) {
      AliError(Form("RunSimulation no loader for %s found\n Number of events per file not set for %s %s", 
                    detName, typeName, detName));
      continue;
    }
    AliDataLoader* dataLoader = 
      loader->GetDataLoader(typeName);
    if (!dataLoader) {
      AliError(Form("no data loader for %s found\n"
                    "Number of events per file not set for %s %s", 
                    typeName, detName, typeName));
      continue;
    }
    dataLoader->SetNumberOfEventsPerFile(fEventsPerFile[i]->GetUniqueID());
    AliDebug(1, Form("number of events per file set to %d for %s %s",
                     fEventsPerFile[i]->GetUniqueID(), detName, typeName));
  }

  AliInfo("running gAlice");
  AliSysInfo::AddStamp("Start_ProcessRun");

  // Create the Root Tree with one branch per detector
  //Hits moved to begin event -> now we are crating separate tree for each event

  gMC->ProcessRun(nEvents);

  // End of this run, close files
  if(nEvents>0) FinishRun();

  AliSysInfo::AddStamp("Stop_ProcessRun");
  delete runLoader;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunSDigitization(const char* detectors)
{
// run the digitization and produce summable digits
  static Int_t eventNr=0;
  AliCodeTimerAuto("",0) ;

  // initialize CDB storage, run number, set CDB lock
  InitCDB();
  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  SetCDBLock();
  
  AliRunLoader* runLoader = LoadRun();
  if (!runLoader) return kFALSE;

  TString detStr = detectors;
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      AliInfo(Form("creating summable digits for %s", det->GetName()));
      AliCodeTimerStart(Form("creating summable digits for %s", det->GetName()));
      det->Hits2SDigits();
      AliCodeTimerStop(Form("creating summable digits for %s", det->GetName()));
      AliSysInfo::AddStamp(Form("Digit_%s_%d",det->GetName(),eventNr), 0,1, eventNr);
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s",
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }
  eventNr++;
  delete runLoader;

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliSimulation::RunDigitization(const char* detectors, 
				      const char* excludeDetectors)
{
// run the digitization and produce digits from sdigits

  AliCodeTimerAuto("",0)

  // initialize CDB storage, run number, set CDB lock
  InitCDB();
  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  SetCDBLock();
  
  delete AliRunLoader::Instance();
  delete gAlice;
  gAlice = NULL;

  Int_t nStreams = 1;
  if (fBkgrdFileNames) nStreams = fBkgrdFileNames->GetEntriesFast() + 1;
  Int_t signalPerBkgrd = GetNSignalPerBkgrd();
  AliDigitizationInput digInp(nStreams, signalPerBkgrd);
  // digInp.SetEmbeddingFlag(fEmbeddingFlag);
  digInp.SetRegionOfInterest(fRegionOfInterest);
  digInp.SetInputStream(0, fGAliceFileName.Data());
  for (Int_t iStream = 1; iStream < nStreams; iStream++) {
    const char* fileName = ((TObjString*)(fBkgrdFileNames->At(iStream-1)))->GetName();
    digInp.SetInputStream(iStream, fileName);
  }
  TObjArray detArr;
  detArr.SetOwner(kTRUE);
  TString detStr = detectors;
  TString detExcl = excludeDetectors;
  if (!static_cast<AliStream*>(digInp.GetInputStream(0))->ImportgAlice()) {
    AliError("Error occured while getting gAlice from Input 0");
    return kFALSE;
  }
  AliRunLoader* runLoader = AliRunLoader::GetRunLoader(digInp.GetInputStream(0)->GetFolderName());
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (!IsSelected(det->GetName(), detStr) || IsSelected(det->GetName(), detExcl)) continue;
    AliDigitizer* digitizer = det->CreateDigitizer(&digInp);
    if (!digitizer || !digitizer->Init()) {
      AliError(Form("no digitizer for %s", det->GetName()));
      if (fStopOnError) return kFALSE;
      else continue;
    }
    detArr.AddLast(digitizer);    
    AliInfo(Form("Created digitizer from SDigits -> Digits for %s", det->GetName()));    
  }
  //
  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s", detStr.Data()));
    if (fStopOnError) return kFALSE;
  }
  //
  Int_t ndigs = detArr.GetEntriesFast();
  Int_t eventsCreated = 0;
  AliRunLoader* outRl =  digInp.GetOutRunLoader();
  while ((eventsCreated++ < fNEvents) || (fNEvents < 0)) {
    if (!digInp.ConnectInputTrees()) break;
    digInp.InitEvent(); //this must be after call of Connect Input tress.
    if (outRl) outRl->SetEventNumber(eventsCreated-1);
    static_cast<AliStream*>(digInp.GetInputStream(0))->ImportgAlice(); // use gAlice of the first input stream
    for (int id=0;id<ndigs;id++) ((AliDigitizer*)detArr[id])->Digitize("");
    digInp.FinishEvent();
  };
  digInp.FinishGlobal();
  //
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunHitsDigitization(const char* detectors)
{
// run the digitization and produce digits from hits

  AliCodeTimerAuto("",0)

  // initialize CDB storage, run number, set CDB lock
  InitCDB();
  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  SetCDBLock();
  
  AliRunLoader* runLoader = LoadRun("READ");
  if (!runLoader) return kFALSE;

  TString detStr = detectors;
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      AliInfo(Form("creating digits from hits for %s", det->GetName()));
      det->Hits2Digits();
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s", 
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::WriteRawData(const char* detectors, 
				   const char* fileName,
				   Bool_t deleteIntermediateFiles,
				   Bool_t selrawdata)
{
// convert the digits to raw data
// First DDL raw data files for the given detectors are created.
// If a file name is given, the DDL files are then converted to a DATE file.
// If deleteIntermediateFiles is true, the DDL raw files are deleted 
// afterwards.
// If the file name has the extension ".root", the DATE file is converted
// to a root file.
// If deleteIntermediateFiles is true, the DATE file is deleted afterwards.
// 'selrawdata' flag can be used to enable writing of detectors raw data
// accoring to the trigger cluster.

  AliCodeTimerAuto("",0)
  AliSysInfo::AddStamp("WriteRawData_Start");
  
  TString detStr = detectors;
  if (!WriteRawFiles(detStr.Data())) {
    if (fStopOnError) return kFALSE;
  }
  AliSysInfo::AddStamp("WriteRawFiles");

  // run HLT simulation on simulated DDL raw files
  // and produce HLT ddl raw files to be included in date/root file
  // bugfix 2009-06-26: the decision whether to write HLT raw data
  // is taken in RunHLT. Here HLT always needs to be run in order to
  // create HLT digits, unless its switched off. This is due to the
  // special placement of the HLT between the generation of DDL files
  // and conversion to DATE/Root file.
  detStr.ReplaceAll("HLT", "");
  if (!fRunHLT.IsNull()) {
    if (!RunHLT()) {
      if (fStopOnError) return kFALSE;
    }
  }
  AliSysInfo::AddStamp("WriteRawData_RunHLT");

  TString dateFileName(fileName);
  if (!dateFileName.IsNull()) {
    Bool_t rootOutput = dateFileName.EndsWith(".root");
    if (rootOutput) dateFileName += ".date";
    TString selDateFileName;
    if (selrawdata) {
      selDateFileName = "selected.";
      selDateFileName+= dateFileName;
    }
    if (!ConvertRawFilesToDate(dateFileName,selDateFileName)) {
      if (fStopOnError) return kFALSE;
    }
    AliSysInfo::AddStamp("ConvertRawFilesToDate");
    if (deleteIntermediateFiles) {
      AliRunLoader* runLoader = LoadRun("READ");
      if (runLoader) for (Int_t iEvent = 0; 
			  iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
	char command[256];
	snprintf(command, 256, "rm -r raw%d", iEvent);
	gSystem->Exec(command);
      }
      delete runLoader;
    }

    if (rootOutput) {
      if (!ConvertDateToRoot(dateFileName, fileName)) {
	if (fStopOnError) return kFALSE;
      }
      AliSysInfo::AddStamp("ConvertDateToRoot");
      if (deleteIntermediateFiles) {
	gSystem->Unlink(dateFileName);
      }
      if (selrawdata) {
	TString selFileName = "selected.";
	selFileName        += fileName;
	if (!ConvertDateToRoot(selDateFileName, selFileName)) {
	  if (fStopOnError) return kFALSE;
	}
	if (deleteIntermediateFiles) {
	  gSystem->Unlink(selDateFileName);
	}
      }
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::WriteRawFiles(const char* detectors)
{
// convert the digits to raw data DDL files

  AliCodeTimerAuto("",0)
  
  AliRunLoader* runLoader = LoadRun("READ");
  if (!runLoader) return kFALSE;

  // write raw data to DDL files
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    AliInfo(Form("processing event %d", iEvent));
    runLoader->GetEvent(iEvent);
    TString baseDir = gSystem->WorkingDirectory();
    char dirName[256];
    snprintf(dirName, 256, "raw%d", iEvent);
    gSystem->MakeDirectory(dirName);
    if (!gSystem->ChangeDirectory(dirName)) {
      AliError(Form("couldn't change to directory %s", dirName));
      if (fStopOnError) return kFALSE; else continue;
    }

    ofstream runNbFile(Form("run%u",runLoader->GetHeader()->GetRun()));
    runNbFile.close();

    TString detStr = detectors;
    if (IsSelected("HLT", detStr)) {
      // Do nothing. "HLT" will be removed from detStr and HLT raw
      // data files are generated in RunHLT.
    }

    TObjArray* detArray = runLoader->GetAliRun()->Detectors();
    for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
      AliModule* det = (AliModule*) detArray->At(iDet);
      if (!det || !det->IsActive()) continue;
      if (IsSelected(det->GetName(), detStr)) {
	AliInfo(Form("creating raw data from digits for %s", det->GetName()));
	det->Digits2Raw();
      }
    }

    if (!WriteTriggerRawData())
      if (fStopOnError) return kFALSE;

    gSystem->ChangeDirectory(baseDir);
    if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
      AliError(Form("the following detectors were not found: %s", 
                    detStr.Data()));
      if (fStopOnError) return kFALSE;
    }
  }

  delete runLoader;
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::ConvertRawFilesToDate(const char* dateFileName,
					    const char* selDateFileName)
{
// convert raw data DDL files to a DATE file with the program "dateStream"
// The second argument is not empty when the user decides to write
// the detectors raw data according to the trigger cluster.

  AliCodeTimerAuto("",0)
  
  char* path = gSystem->Which(gSystem->Getenv("PATH"), "dateStream");
  if (!path) {
    AliError("the program dateStream was not found");
    if (fStopOnError) return kFALSE;
  } else {
    delete[] path;
  }

  AliRunLoader* runLoader = LoadRun("READ");
  if (!runLoader) return kFALSE;

  AliInfo(Form("converting raw data DDL files to DATE file %s", dateFileName));
  Bool_t selrawdata = kFALSE;
  if (strcmp(selDateFileName,"") != 0) selrawdata = kTRUE;

  char command[256];
  // Note the option -s. It is used in order to avoid
  // the generation of SOR/EOR events.
  snprintf(command, 256, "dateStream -c -s -D -o %s -# %d -C -run %d", 
	  dateFileName, runLoader->GetNumberOfEvents(),runLoader->GetHeader()->GetRun());
  FILE* pipe = gSystem->OpenPipe(command, "w");

  if (!pipe) {
    AliError(Form("Cannot execute command: %s",command));
    return kFALSE;
  }

  Int_t selEvents = 0;
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {

    UInt_t detectorPattern = 0;
    runLoader->GetEvent(iEvent);
    if (!runLoader->LoadTrigger()) {
      AliCentralTrigger *aCTP = runLoader->GetTrigger();
      detectorPattern = aCTP->GetClusterMask();
      // Check if the event was triggered by CTP
      if (selrawdata) {
	if (aCTP->GetClassMask()) selEvents++;
      }
    }
    else {
      AliWarning("No trigger can be loaded! Some fields in the event header will be empty !");
      if (selrawdata) {
	AliWarning("No trigger can be loaded! Writing of selected raw data is abandoned !");
	selrawdata = kFALSE;
      }
    }

    fprintf(pipe, "GDC DetectorPattern %u Timestamp %ld\n", detectorPattern, runLoader->GetHeader()->GetTimeStamp());
    Float_t ldc = 0;
    Int_t prevLDC = -1;

    // loop over detectors and DDLs
    for (Int_t iDet = 0; iDet < AliDAQ::kNDetectors; iDet++) {
      for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls(iDet); iDDL++) {

        Int_t ddlID = AliDAQ::DdlID(iDet,iDDL);
        Int_t ldcID = Int_t(ldc + 0.0001);
        ldc += AliDAQ::NumberOfLdcs(iDet) / AliDAQ::NumberOfDdls(iDet);

        char rawFileName[256];
        snprintf(rawFileName, 256, "raw%d/%s", 
                iEvent, AliDAQ::DdlFileName(iDet,iDDL));

	// check existence and size of raw data file
        FILE* file = fopen(rawFileName, "rb");
        if (!file) continue;
        fseek(file, 0, SEEK_END);
        unsigned long size = ftell(file);
	fclose(file);
        if (!size) continue;

        if (ldcID != prevLDC) {
          fprintf(pipe, " LDC Id %d\n", ldcID);
          prevLDC = ldcID;
        }
        fprintf(pipe, "  Equipment Id %d Payload %s\n", ddlID, rawFileName);
      }
    }
  }

  Int_t result = gSystem->ClosePipe(pipe);

  if (!(selrawdata && selEvents > 0)) {
    delete runLoader;
    return (result == 0);
  }

  AliInfo(Form("converting selected by trigger cluster raw data DDL files to DATE file %s", selDateFileName));
  
  snprintf(command, 256, "dateStream -c -s -D -o %s -# %d -C -run %d", 
	  selDateFileName,selEvents,runLoader->GetHeader()->GetRun());
  FILE* pipe2 = gSystem->OpenPipe(command, "w");

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {

    // Get the trigger decision and cluster
    UInt_t detectorPattern = 0;
    TString detClust;
    runLoader->GetEvent(iEvent);
    if (!runLoader->LoadTrigger()) {
      AliCentralTrigger *aCTP = runLoader->GetTrigger();
      if (aCTP->GetClassMask() == 0) continue;
      detectorPattern = aCTP->GetClusterMask();
      detClust = AliDAQ::ListOfTriggeredDetectors(detectorPattern);
      AliInfo(Form("List of detectors to be read out: %s",detClust.Data()));
    }

    fprintf(pipe2, "GDC DetectorPattern %u Timestamp %ld\n", detectorPattern, runLoader->GetHeader()->GetTimeStamp());
    Float_t ldc = 0;
    Int_t prevLDC = -1;

    // loop over detectors and DDLs
    for (Int_t iDet = 0; iDet < AliDAQ::kNDetectors; iDet++) {
      // Write only raw data from detectors that
      // are contained in the trigger cluster(s)
      if (!IsSelected(AliDAQ::DetectorName(iDet),detClust)) continue;

      for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls(iDet); iDDL++) {

        Int_t ddlID = AliDAQ::DdlID(iDet,iDDL);
        Int_t ldcID = Int_t(ldc + 0.0001);
        ldc += AliDAQ::NumberOfLdcs(iDet) / AliDAQ::NumberOfDdls(iDet);

        char rawFileName[256];
        snprintf(rawFileName, 256, "raw%d/%s", 
                iEvent, AliDAQ::DdlFileName(iDet,iDDL));

	// check existence and size of raw data file
        FILE* file = fopen(rawFileName, "rb");
        if (!file) continue;
        fseek(file, 0, SEEK_END);
        unsigned long size = ftell(file);
	fclose(file);
        if (!size) continue;

        if (ldcID != prevLDC) {
          fprintf(pipe2, " LDC Id %d\n", ldcID);
          prevLDC = ldcID;
        }
        fprintf(pipe2, "  Equipment Id %d Payload %s\n", ddlID, rawFileName);
      }
    }
  }

  Int_t result2 = gSystem->ClosePipe(pipe2);

  delete runLoader;
  return ((result == 0) && (result2 == 0));
}

//_____________________________________________________________________________
Bool_t AliSimulation::ConvertDateToRoot(const char* dateFileName,
					const char* rootFileName)
{
// convert a DATE file to a root file with the program "alimdc"

  // ALIMDC setup
  const Int_t kDBSize = 2000000000;
  const Int_t kTagDBSize = 1000000000;
  const Bool_t kFilter = kFALSE;
  const Int_t kCompression = 1;

  char* path = gSystem->Which(gSystem->Getenv("PATH"), "alimdc");
  if (!path) {
    AliError("the program alimdc was not found");
    if (fStopOnError) return kFALSE;
  } else {
    delete[] path;
  }

  AliInfo(Form("converting DATE file %s to root file %s", 
               dateFileName, rootFileName));

  const char* rawDBFS[2] = { "/tmp/mdc1", "/tmp/mdc2" };
  const char* tagDBFS    = "/tmp/mdc1/tags";

  // User defined file system locations
  if (gSystem->Getenv("ALIMDC_RAWDB1")) 
    rawDBFS[0] = gSystem->Getenv("ALIMDC_RAWDB1");
  if (gSystem->Getenv("ALIMDC_RAWDB2")) 
    rawDBFS[1] = gSystem->Getenv("ALIMDC_RAWDB2");
  if (gSystem->Getenv("ALIMDC_TAGDB")) 
    tagDBFS = gSystem->Getenv("ALIMDC_TAGDB");

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));

  gSystem->Exec(Form("mkdir %s",rawDBFS[0]));
  gSystem->Exec(Form("mkdir %s",rawDBFS[1]));
  gSystem->Exec(Form("mkdir %s",tagDBFS));

  Int_t result = gSystem->Exec(Form("alimdc %d %d %d %d %s", 
				    kDBSize, kTagDBSize, kFilter, kCompression, dateFileName));
  gSystem->Exec(Form("mv %s/*.root %s", rawDBFS[0], rootFileName));

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));

  return (result == 0);
}


//_____________________________________________________________________________
AliRunLoader* AliSimulation::LoadRun(const char* mode) const
{
// delete existing run loaders, open a new one and load gAlice

  delete AliRunLoader::Instance();
  AliRunLoader* runLoader = 
    AliRunLoader::Open(fGAliceFileName.Data(), 
		       AliConfig::GetDefaultEventFolderName(), mode);
  if (!runLoader) {
    AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
    return NULL;
  }
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    AliError(Form("no gAlice object found in file %s", 
                  fGAliceFileName.Data()));
    return NULL;
  }
  return runLoader;
}

//_____________________________________________________________________________
Int_t AliSimulation::GetNSignalPerBkgrd(Int_t nEvents) const
{
// get or calculate the number of signal events per background event

  if (!fBkgrdFileNames) return 1;
  Int_t nBkgrdFiles = fBkgrdFileNames->GetEntriesFast();
  if (nBkgrdFiles == 0) return 1;

  // get the number of signal events
  if (nEvents <= 0) {
    AliRunLoader* runLoader = 
	AliRunLoader::Open(fGAliceFileName.Data(), "SIGNAL");
    if (!runLoader) return 1;
    
    nEvents = runLoader->GetNumberOfEvents();
    delete runLoader;
  }

  Int_t result = 0;
  for (Int_t iBkgrdFile = 0; iBkgrdFile < nBkgrdFiles; iBkgrdFile++) {
    // get the number of background events
    const char* fileName = ((TObjString*)
			    (fBkgrdFileNames->At(iBkgrdFile)))->GetName();
    AliRunLoader* runLoader =
      AliRunLoader::Open(fileName, "BKGRD");
    if (!runLoader) continue;
    Int_t nBkgrdEvents = runLoader->GetNumberOfEvents();
    delete runLoader;
  
    // get or calculate the number of signal per background events
    Int_t nSignalPerBkgrd = fBkgrdFileNames->At(iBkgrdFile)->GetUniqueID();
    if (nSignalPerBkgrd <= 0) {
      nSignalPerBkgrd = (nEvents-1) / nBkgrdEvents + 1;
    } else if (result && (result != nSignalPerBkgrd)) {
      AliInfo(Form("the number of signal events per background event "
                   "will be changed from %d to %d for stream %d", 
                   nSignalPerBkgrd, result, iBkgrdFile+1));
      nSignalPerBkgrd = result;
    }

    if (!result) result = nSignalPerBkgrd;
    if (nSignalPerBkgrd * nBkgrdEvents < nEvents) {
      AliWarning(Form("not enough background events (%d) for %d signal events "
                      "using %d signal per background events for stream %d",
                      nBkgrdEvents, nEvents, nSignalPerBkgrd, iBkgrdFile+1));
    }
  }

  return result;
}

//_____________________________________________________________________________
Bool_t AliSimulation::IsSelected(TString detName, TString& detectors) const
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
Bool_t AliSimulation::ConvertRaw2SDigits(const char* rawDirectory, const char* esdFileName, Int_t N) 
{
//
// Steering routine  to convert raw data in directory rawDirectory/ to fake SDigits. 
// These can be used for embedding of MC tracks into RAW data using the standard 
// merging procedure.
//
// If an ESD file is given the reconstructed vertex is taken from it and stored in the event header.
//
  if (!gAlice) {
    AliError("no gAlice object. Restart aliroot and try again.");
    return kFALSE;
  }
  if (gAlice->Modules()->GetEntries() > 0) {
    AliError("gAlice was already run. Restart aliroot and try again.");
    return kFALSE;
  }
  
  AliInfo(Form("initializing gAlice with config file %s",fConfigFileName.Data()));
  
  gAlice->Announce();
  
  gROOT->LoadMacro(fConfigFileName.Data());
  gInterpreter->ProcessLine(gAlice->GetConfigFunction());
  
  if(AliCDBManager::Instance()->GetRun() >= 0) { 
  	SetRunNumber(AliCDBManager::Instance()->GetRun());
  } else {
  	AliWarning("Run number not initialized!!");
  }
  
   AliRunLoader::Instance()->CdGAFile();
    
   AliPDG::AddParticlesToPdgDataBase();  

   gMC->SetMagField(TGeoGlobalMagField::Instance()->GetField());
   
   gAlice->GetMCApp()->Init();

   //Must be here because some MCs (G4) adds detectors here and not in Config.C
   gAlice->InitLoaders();
   AliRunLoader::Instance()->MakeTree("E");
   AliRunLoader::Instance()->LoadKinematics("RECREATE");
   AliRunLoader::Instance()->LoadTrackRefs("RECREATE");
   AliRunLoader::Instance()->LoadHits("all","RECREATE");

   //
   // Save stuff at the beginning of the file to avoid file corruption
   AliRunLoader::Instance()->CdGAFile();
   gAlice->Write();
//
//  Initialize CDB     
    InitCDB();
    //AliCDBManager* man = AliCDBManager::Instance();
    //man->SetRun(0); // Should this come from rawdata header ?
    
    Int_t iDet;
    //
    // Get the runloader
    AliRunLoader* runLoader = AliRunLoader::Instance();
    //
    // Open esd file if available
    TFile* esdFile = 0;
    TTree* treeESD = 0;
    AliESDEvent* esd = 0;
    if (esdFileName && (strlen(esdFileName)>0)) {
      esdFile = TFile::Open(esdFileName);
      if (esdFile) {
        esd = new AliESDEvent();
        esdFile->GetObject("esdTree", treeESD);
        if (treeESD) esd->ReadFromTree(treeESD);
      }
    }

    //
    // Create the RawReader
    TString fileName(rawDirectory);
    AliRawReader* rawReader = 0x0;
    if (fileName.EndsWith("/")) {
      rawReader = new AliRawReaderFile(fileName);
    } else if (fileName.EndsWith(".root")) {
      rawReader = new AliRawReaderRoot(fileName);
    } else if (!fileName.IsNull()) {
      rawReader = new AliRawReaderDate(fileName);
    }
    if (!rawReader) return (kFALSE);
    
//     if (!fEquipIdMap.IsNull() && fRawReader)
//       fRawReader->LoadEquipmentIdsMap(fEquipIdMap);
    //
    // Get list of detectors
    TObjArray* detArray = runLoader->GetAliRun()->Detectors();
    //
    // Get Header
    AliHeader* header = runLoader->GetHeader();
    // Event loop
    Int_t nev = 0;
    while(kTRUE) {
	if (!(rawReader->NextEvent())) break;
	runLoader->SetEventNumber(nev);
        runLoader->GetHeader()->Reset(rawReader->GetRunNumber(), 
                                      nev, nev);
        runLoader->GetEvent(nev);
        AliInfo(Form("We are at event %d",nev));
	//
	// Detector loop
        TString detStr = fMakeSDigits;
	for (iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
	    AliModule* det = (AliModule*) detArray->At(iDet);
	    if (!det || !det->IsActive()) continue;
	    if (IsSelected(det->GetName(), detStr)) {
	      AliInfo(Form("Calling Raw2SDigits for %s", det->GetName()));
	      det->Raw2SDigits(rawReader);
	      rawReader->Reset();
	    }
	} // detectors
        

	//
	//  If ESD information available obtain reconstructed vertex and store in header.
	if (treeESD) {
	    treeESD->GetEvent(nev);
	    const AliESDVertex* esdVertex = esd->GetPrimaryVertex();
	    Double_t position[3];
	    esdVertex->GetXYZ(position);
	    AliGenEventHeader* mcHeader = new  AliGenEventHeader("ESD");
	    TArrayF mcV;
	    mcV.Set(3);
	    for (Int_t i = 0; i < 3; i++) mcV[i] = position[i];
	    mcHeader->SetPrimaryVertex(mcV);
	    header->Reset(0,nev);
	    header->SetGenEventHeader(mcHeader);
	    AliInfo(Form("***** Saved vertex %f %f %f \n", position[0], position[1], position[2]));
	}
//
//      Finish the event
	runLoader->TreeE()->Fill();
        AliInfo(Form("Finished event %d",nev));
	nev++;
        if (N>0&&nev>=N)
          break;
    } // events
 
    delete rawReader;
//
//  Finish the run 
    runLoader->CdGAFile();
    runLoader->WriteHeader("OVERWRITE");
    runLoader->WriteRunLoader();

    return kTRUE;
}

//_____________________________________________________________________________
void AliSimulation::FinishRun()
{
  //
  // Called at the end of the run.
  //

  if(IsLegoRun()) 
   {
    AliDebug(1, "Finish Lego");
    AliRunLoader::Instance()->CdGAFile();
    fLego->FinishRun();
   }
  
  // Clean detector information
  TIter next(gAlice->Modules());
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    AliDebug(2, Form("%s->FinishRun()", detector->GetName()));
    detector->FinishRun();
  }
  
  AliDebug(1, "AliRunLoader::Instance()->WriteHeader(OVERWRITE)");
  AliRunLoader::Instance()->WriteHeader("OVERWRITE");

  // Write AliRun info and all detectors parameters
  AliRunLoader::Instance()->CdGAFile();
  gAlice->Write(0,TObject::kOverwrite);//write AliRun
  AliRunLoader::Instance()->Write(0,TObject::kOverwrite);//write RunLoader itself
  
  if(gAlice->GetMCApp()) gAlice->GetMCApp()->FinishRun();  
  AliRunLoader::Instance()->Synchronize();
}

//_____________________________________________________________________________
Int_t AliSimulation::GetDetIndex(const char* detector)
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
Bool_t AliSimulation::CreateHLT()
{
  // Init the HLT simulation.
  // The function  loads the library and creates the instance of AliHLTSimulation.
  // the main reason for the decoupled creation is to set the transient OCDB
  // objects before the OCDB is locked

  // load the library dynamically
  gSystem->Load(ALIHLTSIMULATION_LIBRARY);

  // check for the library version
  AliHLTSimulationGetLibraryVersion_t fctVersion=(AliHLTSimulationGetLibraryVersion_t)(gSystem->DynFindSymbol(ALIHLTSIMULATION_LIBRARY, ALIHLTSIMULATION_GET_LIBRARY_VERSION));
  if (!fctVersion) {
    AliError(Form("can not load library %s", ALIHLTSIMULATION_LIBRARY));
    return kFALSE;
  }
  if (fctVersion()!= ALIHLTSIMULATION_LIBRARY_VERSION) {
    AliWarning(Form("%s version does not match: compiled for version %d, loaded %d", ALIHLTSIMULATION_LIBRARY, ALIHLTSIMULATION_LIBRARY_VERSION, fctVersion()));
  }

  // print compile info
  typedef void (*CompileInfo)( const char*& date, const char*& time);
  CompileInfo fctInfo=(CompileInfo)gSystem->DynFindSymbol(ALIHLTSIMULATION_LIBRARY, "CompileInfo");
  if (fctInfo) {
    const char* date="";
    const char* time="";
    (*fctInfo)(date, time);
    if (!date) date="unknown";
    if (!time) time="unknown";
    AliInfo(Form("%s build on %s (%s)", ALIHLTSIMULATION_LIBRARY, date, time));
  } else {
    AliInfo(Form("no build info available for %s", ALIHLTSIMULATION_LIBRARY));
  }

  // create instance of the HLT simulation
  AliHLTSimulationCreateInstance_t fctCreate=(AliHLTSimulationCreateInstance_t)(gSystem->DynFindSymbol(ALIHLTSIMULATION_LIBRARY, ALIHLTSIMULATION_CREATE_INSTANCE));
  if (fctCreate==NULL || (fpHLT=(fctCreate()))==NULL) {
    AliError(Form("can not create instance of HLT simulation (creator %p)", fctCreate));
    return kFALSE;    
  }

  TString specObjects;
  for (Int_t i = 0; i < fSpecCDBUri.GetEntriesFast(); i++) {
    if (specObjects.Length()>0) specObjects+=" ";
    specObjects+=fSpecCDBUri[i]->GetName();
  }

  AliHLTSimulationSetup_t fctSetup=(AliHLTSimulationSetup_t)(gSystem->DynFindSymbol(ALIHLTSIMULATION_LIBRARY, ALIHLTSIMULATION_SETUP));
  if (fctSetup==NULL || fctSetup(fpHLT, this, specObjects.Data())<0) {
    AliWarning(Form("failed to setup HLT simulation (function %p)", fctSetup));
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunHLT()
{
  // Run the HLT simulation
  // HLT simulation is implemented in HLT/sim/AliHLTSimulation
  // Disabled if fRunHLT is empty, default vaule is "default".
  // AliSimulation::SetRunHLT can be used to set the options for HLT simulation
  // The default simulation depends on the HLT component libraries and their
  // corresponding agents which define components and chains to run. See
  // http://web.ift.uib.no/~kjeks/doc/alice-hlt-current/
  // http://web.ift.uib.no/~kjeks/doc/alice-hlt-current/classAliHLTModuleAgent.html
  //
  // The libraries to be loaded can be specified as an option.
  // <pre>
  // AliSimulation sim;
  // sim.SetRunHLT("libAliHLTSample.so");
  // </pre>
  // will only load <tt>libAliHLTSample.so</tt>

  // Other available options:
  // \li loglevel=<i>level</i> <br>
  //     logging level for this processing
  // \li alilog=off
  //     disable redirection of log messages to AliLog class
  // \li config=<i>macro</i>
  //     configuration macro
  // \li chains=<i>configuration</i>
  //     comma separated list of configurations to be run during simulation
  // \li rawfile=<i>file</i>
  //     source for the RawReader to be created, the default is <i>./</i> if
  //     raw data is simulated

  int iResult=0;

  if (!fpHLT && !CreateHLT()) {
    return kFALSE;
  }
  AliHLTSimulation* pHLT=fpHLT;

  AliRunLoader* pRunLoader = LoadRun("READ");
  if (!pRunLoader) return kFALSE;

  // initialize CDB storage, run number, set CDB lock
  // thats for the case of running HLT simulation without all the other steps
  // multiple calls are handled by the function, so we can just call
  InitCDB();
  if (!SetRunNumberFromData()) if (fStopOnError) return kFALSE;
  SetCDBLock();
  
  // init the HLT simulation
  TString options;
  if (fRunHLT.CompareTo("default")!=0) options=fRunHLT;
  TString detStr = fWriteRawData;
  if (!IsSelected("HLT", detStr)) {
    options+=" writerawfiles=";
  } else {
    options+=" writerawfiles=HLT";
  }

  if (!detStr.IsNull() && !options.Contains("rawfile=")) {
    // as a matter of fact, HLT will run reconstruction and needs the RawReader
    // in order to get detector data. By default, RawReaderFile is used to read
    // the already simulated ddl files. Date and Root files from the raw data
    // are generated after the HLT simulation.
    options+=" rawfile=./";
  }

  AliHLTSimulationInit_t fctInit=(AliHLTSimulationInit_t)(gSystem->DynFindSymbol(ALIHLTSIMULATION_LIBRARY, ALIHLTSIMULATION_INIT));
  if (fctInit==NULL || (iResult=(fctInit(pHLT, pRunLoader, options.Data())))<0) {
    AliError(Form("can not init HLT simulation: error %d (init %p)", iResult, fctInit));
  } else {
    // run the HLT simulation
    AliHLTSimulationRun_t fctRun=(AliHLTSimulationRun_t)(gSystem->DynFindSymbol(ALIHLTSIMULATION_LIBRARY, ALIHLTSIMULATION_RUN));
    if (fctRun==NULL || (iResult=(fctRun(pHLT, pRunLoader)))<0) {
      AliError(Form("can not run HLT simulation: error %d (run %p)", iResult, fctRun));
    }
  }

  // delete the instance
  AliHLTSimulationDeleteInstance_t fctDelete=(AliHLTSimulationDeleteInstance_t)(gSystem->DynFindSymbol(ALIHLTSIMULATION_LIBRARY, ALIHLTSIMULATION_DELETE_INSTANCE));
  if (fctDelete==NULL || fctDelete(pHLT)<0) {
    AliError(Form("can not delete instance of HLT simulation (creator %p)", fctDelete));
  }
  pHLT=NULL;

  return iResult>=0?kTRUE:kFALSE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunQA()
{
	// run the QA on summable hits, digits or digits
	
    //if(!gAlice) return kFALSE;
	AliQAManager::QAManager()->SetRunLoader(AliRunLoader::Instance()) ;

	TString detectorsw("") ;  
	Bool_t rv = kTRUE ; 
  AliQAManager::QAManager()->SetEventSpecie(fEventSpecie) ;
	detectorsw = AliQAManager::QAManager()->Run(fQADetectors.Data()) ; 
	if ( detectorsw.IsNull() ) 
		rv = kFALSE ; 
	return rv ; 
}

//_____________________________________________________________________________
Bool_t AliSimulation::SetRunQA(TString detAndAction) 
{
	// Allows to run QA for a selected set of detectors
	// and a selected set of tasks among HITS, SDIGITS and DIGITS
	// all selected detectors run the same selected tasks
	
	if (!detAndAction.Contains(":")) {
		AliError( Form("%s is a wrong syntax, use \"DetectorList:ActionList\" \n", detAndAction.Data()) ) ;
		fRunQA = kFALSE ;
		return kFALSE ; 		
	}
	Int_t colon = detAndAction.Index(":") ; 
	fQADetectors = detAndAction(0, colon) ; 
	if (fQADetectors.Contains("ALL") ){
    TString tmp = Form("%s %s", fMakeDigits.Data(), fMakeDigitsFromHits.Data()) ; 
    Int_t minus = fQADetectors.Last('-') ; 
    TString toKeep = Form("%s %s", fMakeDigits.Data(), fMakeDigitsFromHits.Data()) ; 
    TString toRemove("") ;
    while (minus >= 0) {
      toRemove = fQADetectors(minus+1, fQADetectors.Length()) ; 
      toRemove = toRemove.Strip() ; 
      toKeep.ReplaceAll(toRemove, "") ; 
      fQADetectors.ReplaceAll(Form("-%s", toRemove.Data()), "") ; 
      minus = fQADetectors.Last('-') ; 
    }
    fQADetectors = toKeep ; 
  }
  fQATasks   = detAndAction(colon+1, detAndAction.Sizeof() ) ; 
	if (fQATasks.Contains("ALL") ) {
		fQATasks = Form("%d %d %d", AliQAv1::kHITS, AliQAv1::kSDIGITS, AliQAv1::kDIGITS) ; 
	} else {
		fQATasks.ToUpper() ; 
		TString tempo("") ; 
		if ( fQATasks.Contains("HIT") ) 
			tempo = Form("%d ", AliQAv1::kHITS) ; 
		if ( fQATasks.Contains("SDIGIT") ) 
			tempo += Form("%d ", AliQAv1::kSDIGITS) ; 
		if ( fQATasks.Contains("DIGIT") ) 
			tempo += Form("%d ", AliQAv1::kDIGITS) ; 
		fQATasks = tempo ; 
		if (fQATasks.IsNull()) {
			AliInfo("No QA requested\n")  ;
			fRunQA = kFALSE ;
			return kTRUE ; 
		}
	}	
	TString tempo(fQATasks) ; 
    tempo.ReplaceAll(Form("%d", AliQAv1::kHITS), AliQAv1::GetTaskName(AliQAv1::kHITS)) 	;
    tempo.ReplaceAll(Form("%d", AliQAv1::kSDIGITS), AliQAv1::GetTaskName(AliQAv1::kSDIGITS)) ;	
    tempo.ReplaceAll(Form("%d", AliQAv1::kDIGITS), AliQAv1::GetTaskName(AliQAv1::kDIGITS)) ; 	
	AliInfo( Form("QA will be done on \"%s\" for \"%s\"\n", fQADetectors.Data(), tempo.Data()) ) ;  
	fRunQA = kTRUE ;
	AliQAManager::QAManager()->SetActiveDetectors(fQADetectors) ; 
	AliQAManager::QAManager()->SetTasks(fQATasks) ; 
  for (Int_t det = 0 ; det < AliQAv1::kNDET ; det++) 
    AliQAManager::QAManager()->SetWriteExpert(AliQAv1::DETECTORINDEX_t(det)) ;
  
	return kTRUE; 
} 

//_____________________________________________________________________________
void AliSimulation::ProcessEnvironmentVars()
{
// Extract run number and random generator seed from env variables

    AliInfo("Processing environment variables");
    
    // Random Number seed
    
    // first check that seed is not already set
    if (fSeed == 0) {
    	if (gSystem->Getenv("CONFIG_SEED")) {
     	 	fSeed = atoi(gSystem->Getenv("CONFIG_SEED"));
    	}
    } else {
    	if (gSystem->Getenv("CONFIG_SEED")) {
    		AliInfo(Form("Seed for random number generation already set (%d)"
			     ": CONFIG_SEED variable ignored!", fSeed));
    	}
    }
   
    AliInfo(Form("Seed for random number generation = %d ", fSeed)); 

    // Run Number
    
    // first check that run number is not already set
    if(fRun < 0) {    
    	if (gSystem->Getenv("DC_RUN")) {
		fRun = atoi(gSystem->Getenv("DC_RUN"));
    	}
    } else {
    	if (gSystem->Getenv("DC_RUN")) {
    		AliInfo(Form("Run number already set (%d): DC_RUN variable ignored!", fRun));
    	}
    }
    
    AliInfo(Form("Run number = %d", fRun)); 
}

//---------------------------------------------------------------------
void AliSimulation::WriteGRPEntry()
{
  // Get the necessary information from galice (generator, trigger etc) and
  // write a GRP entry corresponding to the settings in the Config.C used
  // note that Hall probes and Cavern and Surface Atmos pressures are not simulated.


  AliInfo("Writing global run parameters entry into the OCDB");

  AliGRPObject* grpObj = new AliGRPObject();

  grpObj->SetRunType("PHYSICS");
  grpObj->SetTimeStart(fTimeStart);
  grpObj->SetTimeEnd(fTimeEnd); 
  grpObj->SetBeamEnergyIsSqrtSHalfGeV(); // new format of GRP: store sqrt(s)/2 in GeV

  const AliGenerator *gen = gAlice->GetMCApp()->Generator();
  Int_t a = 0;
  Int_t z = 0;

  if (gen) {
    TString projectile;
    gen->GetProjectile(projectile,a,z);
    TString target;
    gen->GetTarget(target,a,z);
    TString beamType = projectile + "-" + target;
    beamType.ReplaceAll(" ","");
    if (!beamType.CompareTo("-")) {
      grpObj->SetBeamType("UNKNOWN");
      grpObj->SetBeamEnergy(gen->GetEnergyCMS()/2);
    }
    else {
      grpObj->SetBeamType(beamType);
      if (z != 0) {
	  grpObj->SetBeamEnergy(gen->GetEnergyCMS()/2 * a / z);
      } else {
	  grpObj->SetBeamEnergy(gen->GetEnergyCMS()/2 );
      }
      // Heavy ion run, the event specie is set to kHighMult
      fEventSpecie = AliRecoParam::kHighMult;
      if ((strcmp(beamType,"p-p") == 0) ||
          (strcmp(beamType,"p-")  == 0) ||
          (strcmp(beamType,"-p")  == 0) ||
          (strcmp(beamType,"P-P") == 0) ||
          (strcmp(beamType,"P-")  == 0) ||
          (strcmp(beamType,"-P")  == 0)) {
        // Proton run, the event specie is set to kLowMult
        fEventSpecie = AliRecoParam::kLowMult;
      } 
    }
  } else {
    AliWarning("Unknown beam type and energy! Setting energy to 0");
    grpObj->SetBeamEnergy(0);
    grpObj->SetBeamType("UNKNOWN");
  }

  UInt_t detectorPattern  = 0;
  Int_t nDets = 0;
  TObjArray *detArray = gAlice->Detectors();
  for (Int_t iDet = 0; iDet < AliDAQ::kNDetectors-1; iDet++) {
    if (detArray->FindObject(AliDAQ::OfflineModuleName(iDet))) {
      AliDebug(1, Form("Detector #%d found: %s", iDet, AliDAQ::OfflineModuleName(iDet)));
      detectorPattern |= (1 << iDet);
      nDets++;
    }
  }
  // CTP
  if (!fTriggerConfig.IsNull())
    detectorPattern |= (1 << AliDAQ::DetectorID("TRG"));

  // HLT
  if (!fRunHLT.IsNull())
    detectorPattern |= (1 << AliDAQ::kHLTId);

  grpObj->SetNumberOfDetectors((Char_t)nDets);
  grpObj->SetDetectorMask((Int_t)detectorPattern);
  grpObj->SetLHCPeriod("LHC08c");
  grpObj->SetLHCState("STABLE_BEAMS");
  //
  AliMagF *field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  Float_t solenoidField = field ? TMath::Abs(field->SolenoidField()) : 0;

  Float_t factorSol     = field ? field->GetFactorSol() : 0;
  Float_t currentSol    = TMath::Abs(factorSol)>1E-6 ? 
    TMath::Nint(TMath::Abs(solenoidField/factorSol))/5.*30000.*TMath::Abs(factorSol) : 0;
  //
  Float_t factorDip     = field ? field->GetFactorDip() : 0;
  Float_t currentDip    = 6000.*TMath::Abs(factorDip);
  //
  grpObj->SetL3Current(currentSol,(AliGRPObject::Stats)0);
  grpObj->SetDipoleCurrent(currentDip,(AliGRPObject::Stats)0);  
  grpObj->SetL3Polarity(factorSol>0 ? 0:1);  
  grpObj->SetDipolePolarity(factorDip>0 ? 0:1);
  if (field) grpObj->SetUniformBMap(field->IsUniform()); // for special MC with k5kGUniform map
  grpObj->SetPolarityConventionLHC();                    // LHC convention +/+ current -> -/- field main components
  //
  grpObj->SetCavernTemperature(0,(AliGRPObject::Stats)0);
  
  //grpMap->Add(new TObjString("fCavernPressure"),new TObjString("0")); ---> not inserted in simulation with the new object, since it is now an AliDCSSensor

  // Now store the entry in OCDB
  AliCDBManager* man = AliCDBManager::Instance();
  
  man->SetLock(0, fKey);
  
  AliCDBStorage* sto = man->GetStorage(fGRPWriteLocation.Data());
  

  AliCDBId id("GRP/GRP/Data", man->GetRun(), man->GetRun(), 1, 1);
  AliCDBMetaData *metadata= new AliCDBMetaData();

  metadata->SetResponsible("alice-off@cern.ch");
  metadata->SetComment("Automatically produced GRP entry for Monte Carlo");
 
  sto->Put(grpObj,id,metadata);
  man->SetLock(1, fKey);
}

//_____________________________________________________________________________
time_t AliSimulation::GenerateTimeStamp() const
{
  // Generate event time-stamp according to
  // SOR/EOR time from GRP
  if (fUseTimeStampFromCDB)
    return fTimeStart + gRandom->Integer(fTimeEnd-fTimeStart);
  else
    return 0;
}
