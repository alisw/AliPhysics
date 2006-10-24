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

#include <TGeoManager.h>
#include <TObjString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TFile.h>

#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliAlignObj.h"
#include "AliCentralTrigger.h"
#include "AliDAQ.h"
#include "AliDigitizer.h"
#include "AliGenerator.h"
#include "AliLog.h"
#include "AliModule.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliSimulation.h"
#include "AliVertexGenFile.h"
#include "AliCentralTrigger.h"
#include "AliCTPRawData.h"
#include "AliRawReaderFile.h"
#include "AliESD.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

ClassImp(AliSimulation)


//_____________________________________________________________________________
AliSimulation::AliSimulation(const char* configFileName, const char* cdbUri,
			     const char* name, const char* title) :
  TNamed(name, title),

  fRunGeneration(kTRUE),
  fRunSimulation(kTRUE),
  fLoadAlignFromCDB(kTRUE),
  fLoadAlignData("ALL"),
  fMakeSDigits("ALL"),
  fMakeDigits("ALL"),
  fMakeTrigger(""),
  fMakeDigitsFromHits(""),
  fWriteRawData(""),
  fRawDataFileName(""),
  fDeleteIntermediateFiles(kFALSE),
  fStopOnError(kFALSE),

  fNEvents(1),
  fConfigFileName(configFileName),
  fGAliceFileName("galice.root"),
  fEventsPerFile(),
  fBkgrdFileNames(NULL),
  fAlignObjArray(NULL),
  fUseBkgrdVertex(kTRUE),
  fRegionOfInterest(kFALSE),
  fCDBUri(cdbUri),
  fSpecCDBUri(),
  fEmbeddingFlag(kFALSE)
{
// create simulation object with default parameters

  SetGAliceFile("galice.root");
}

//_____________________________________________________________________________
AliSimulation::AliSimulation(const AliSimulation& sim) :
  TNamed(sim),

  fRunGeneration(sim.fRunGeneration),
  fRunSimulation(sim.fRunSimulation),
  fLoadAlignFromCDB(sim.fLoadAlignFromCDB),
  fLoadAlignData(sim.fLoadAlignData),
  fMakeSDigits(sim.fMakeSDigits),
  fMakeDigits(sim.fMakeDigits),
  fMakeTrigger(sim.fMakeTrigger),
  fMakeDigitsFromHits(sim.fMakeDigitsFromHits),
  fWriteRawData(sim.fWriteRawData),
  fRawDataFileName(""),
  fDeleteIntermediateFiles(kFALSE),
  fStopOnError(sim.fStopOnError),

  fNEvents(sim.fNEvents),
  fConfigFileName(sim.fConfigFileName),
  fGAliceFileName(sim.fGAliceFileName),
  fEventsPerFile(),
  fBkgrdFileNames(NULL),
  fAlignObjArray(NULL),
  fUseBkgrdVertex(sim.fUseBkgrdVertex),
  fRegionOfInterest(sim.fRegionOfInterest),
  fCDBUri(sim.fCDBUri),
  fSpecCDBUri(),
  fEmbeddingFlag(sim.fEmbeddingFlag)
{
// copy constructor

  for (Int_t i = 0; i < sim.fEventsPerFile.GetEntriesFast(); i++) {
    if (!sim.fEventsPerFile[i]) continue;
    fEventsPerFile.Add(sim.fEventsPerFile[i]->Clone());
  }

  fBkgrdFileNames = new TObjArray;
  for (Int_t i = 0; i < sim.fBkgrdFileNames->GetEntriesFast(); i++) {
    if (!sim.fBkgrdFileNames->At(i)) continue;
    fBkgrdFileNames->Add(sim.fBkgrdFileNames->At(i)->Clone());
  }

  for (Int_t i = 0; i < sim.fSpecCDBUri.GetEntriesFast(); i++) {
    if (sim.fSpecCDBUri[i]) fSpecCDBUri.Add(sim.fSpecCDBUri[i]->Clone());
  }

}

//_____________________________________________________________________________
AliSimulation& AliSimulation::operator = (const AliSimulation& sim)
{
// assignment operator

  this->~AliSimulation();
  new(this) AliSimulation(sim);
  return *this;
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
}


//_____________________________________________________________________________
void AliSimulation::SetNumberOfEvents(Int_t nEvents)
{
// set the number of events for one run

  fNEvents = nEvents;
}

//_____________________________________________________________________________
void AliSimulation::InitCDBStorage()
{
// activate a default CDB storage
// First check if we have any CDB storage set, because it is used 
// to retrieve the calibration and alignment constants

  AliCDBManager* man = AliCDBManager::Instance();
  if (man->IsDefaultStorageSet())
  {
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("Default CDB storage has been already set !");
    AliWarning(Form("Ignoring the default storage declared in AliReconstruction: %s",fCDBUri.Data()));
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    fCDBUri = "";
  }
  else {
    AliDebug(2,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliDebug(2, Form("Default CDB storage is set to: %s",fCDBUri.Data()));
    AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
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
  man->Print();
}

//_____________________________________________________________________________
void AliSimulation::SetDefaultStorage(const char* uri) {
// Store the desired default CDB storage location
// Activate it later within the Run() method

  fCDBUri = uri;

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
Bool_t AliSimulation::ApplyAlignObjsToGeom(TObjArray* alObjArray)
{
  // Read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName and apply them to the geometry
  // manager singleton.
  //
  alObjArray->Sort();
  Int_t nvols = alObjArray->GetEntriesFast();

  Bool_t flag = kTRUE;

  for(Int_t j=0; j<nvols; j++)
    {
      AliAlignObj* alobj = (AliAlignObj*) alObjArray->UncheckedAt(j);
      if (alobj->ApplyToGeometry() == kFALSE) flag = kFALSE;
    }

  if (AliDebugLevelClass() >= 1) {
    gGeoManager->GetTopNode()->CheckOverlaps(1);
    TObjArray* ovexlist = gGeoManager->GetListOfOverlaps();
    if(ovexlist->GetEntriesFast()){  
      AliError("The application of alignment objects to the geometry caused huge overlaps/extrusions!");
   }
  }

  return flag;

}

//_____________________________________________________________________________
Bool_t AliSimulation::ApplyAlignObjsToGeom(const char* fileName, const char* clArrayName)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the file fileName and apply
  // them to the TGeo geometry passed as argument
  //

  TFile* inFile = TFile::Open(fileName,"READ");
  if (!inFile || !inFile->IsOpen()) {
    AliErrorClass(Form("Could not open file %s !",fileName));
    return kFALSE;
  }

  TClonesArray* alObjArray = ((TClonesArray*) inFile->Get(clArrayName));
  inFile->Close();
  if (!alObjArray) {
    AliErrorClass(Form("Could not get array (%s) from file (%s) !",clArrayName,fileName));
    return kFALSE;
  }

  return ApplyAlignObjsToGeom(alObjArray);

}

//_____________________________________________________________________________
Bool_t AliSimulation::ApplyAlignObjsToGeom(AliCDBParam* param, AliCDBId& Id)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the AliCDBEntry identified by
  // param (to get the AliCDBStorage) and Id; apply the alignment objects
  // to the TGeo geometry passed as argument
  //

  AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage(param);
  AliCDBEntry* entry = storage->Get(Id);
  TClonesArray* AlObjArray = ((TClonesArray*) entry->GetObject());

  return ApplyAlignObjsToGeom(AlObjArray);

}

//_____________________________________________________________________________
Bool_t AliSimulation::ApplyAlignObjsToGeom(const char* uri, const char* path, Int_t runnum, Int_t version, Int_t sversion)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the AliCDBEntry identified by
  // param (to get the AliCDBStorage) and Id; apply the alignment objects
  // to the TGeo geometry passed as argument
  //

  AliCDBParam* param = AliCDBManager::Instance()->CreateParameter(uri);
  AliCDBId id(path, runnum, runnum, version, sversion);

  return ApplyAlignObjsToGeom(param, id);

}

//_____________________________________________________________________________
Bool_t AliSimulation::ApplyAlignObjsToGeom(const char* detName, Int_t runnum, Int_t version, Int_t sversion)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the AliCDBEntry identified by
  // param (to get the AliCDBStorage) and Id; apply the alignment objects
  // to the TGeo geometry passed as argument
  //

  AliCDBPath path(detName,"Align","Data");
  AliCDBEntry* entry = AliCDBManager::Instance()->Get(path.GetPath(),runnum,version,sversion);

  if(!entry) return kFALSE;
  TClonesArray* AlObjArray = ((TClonesArray*) entry->GetObject());

  return ApplyAlignObjsToGeom(AlObjArray);
}

//_____________________________________________________________________________
Bool_t AliSimulation::SetAlignObjArraySingleDet(const char* detName)
{
  // Fills array of single detector's alignable objects from CDB
  
  AliDebug(2, Form("Loading alignment data for detector: %s",detName));
  
  AliCDBEntry *entry;
  	
  AliCDBPath path(detName,"Align","Data");
	
  entry=AliCDBManager::Instance()->Get(path.GetPath());
  if(!entry){ 
  	AliDebug(2,Form("Couldn't load alignment data for detector %s",detName));
	return kFALSE;
  }
  entry->SetOwner(1);
  TClonesArray *alignArray = (TClonesArray*) entry->GetObject();	
  alignArray->SetOwner(0);
  AliDebug(2,Form("Found %d alignment objects for %s",
			alignArray->GetEntries(),detName));

  AliAlignObj *alignObj=0;
  TIter iter(alignArray);
	
  // loop over align objects in detector
  while( ( alignObj=(AliAlignObj *) iter.Next() ) ){
  	fAlignObjArray->Add(alignObj);
  }
  // delete entry --- Don't delete, it is cached!
	
  AliDebug(2, Form("fAlignObjArray entries: %d",fAlignObjArray->GetEntries() ));
  return kTRUE;

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

  Bool_t delRunLoader = kFALSE;
  if (!runLoader) {
    runLoader = LoadRun("READ");
    if (!runLoader) return kFALSE;
    delRunLoader = kTRUE;
  }

  // Load alignment data from CDB and fill fAlignObjArray 
  if(fLoadAlignFromCDB){
  	if(!fAlignObjArray) fAlignObjArray = new TObjArray();
  	
	//fAlignObjArray->RemoveAll(); 
 	fAlignObjArray->Clear();  	
	fAlignObjArray->SetOwner(0);
 
  	TString detStr = fLoadAlignData;
  	TString dataNotLoaded="";
  	TString dataLoaded="";
  
  	TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  	for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
  		AliModule* det = (AliModule*) detArray->At(iDet);
  		if (!det || !det->IsActive()) continue;
  		if (IsSelected(det->GetName(), detStr)) {
  			if(!SetAlignObjArraySingleDet(det->GetName())){
  				dataNotLoaded += det->GetName();
  				dataNotLoaded += " ";
  			} else {
				dataLoaded += det->GetName();
				dataLoaded += " ";
			}
  		}
  	} // end loop over detectors
  
  	if ((detStr.CompareTo("ALL") == 0)) detStr = "";
  	dataNotLoaded += detStr;
  	if(!dataLoaded.IsNull()) AliInfo(Form("Alignment data loaded for: %s",
  			  dataLoaded.Data()));
  	if(!dataNotLoaded.IsNull()) AliInfo(Form("Didn't/couldn't load alignment data for: %s",
  			  dataNotLoaded.Data()));
  } // fLoadAlignFromCDB flag
 
  // Check if the array with alignment objects was
  // provided by the user. If yes, apply the objects
  // to the present TGeo geometry
  if (fAlignObjArray) {
    if (gGeoManager && gGeoManager->IsClosed()) {
      if (ApplyAlignObjsToGeom(fAlignObjArray) == kFALSE) {
	AliError("The misalignment of one or more volumes failed!"
		 "Compare the list of simulated detectors and the list of detector alignment data!");
	if (delRunLoader) delete runLoader;
	return kFALSE;
      }
    }
    else {
      AliError("Can't apply the misalignment! gGeoManager doesn't exist or it is still opened!");
      if (delRunLoader) delete runLoader;
      return kFALSE;
    }
  }


  // Update the internal geometry of modules (ITS needs it)
  TString detStr = fLoadAlignData;
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
Bool_t AliSimulation::SetRunNumber()
{
  // Set the CDB manager run number
  // The run number is retrieved from gAlice

  if(AliCDBManager::Instance()->GetRun() < 0) {
    AliRunLoader* runLoader = LoadRun("READ");
    if (!runLoader) return kFALSE;
    else {
      AliCDBManager::Instance()->SetRun(runLoader->GetAliRun()->GetRunNumber());
      AliInfo(Form("Run number: %d",AliCDBManager::Instance()->GetRun()));
      delete runLoader;
    }
  }
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

  InitCDBStorage();

  if (nEvents > 0) fNEvents = nEvents;

  // generation and simulation -> hits
  if (fRunGeneration) {
    if (!RunSimulation()) if (fStopOnError) return kFALSE;
  }

  // Set run number in CDBManager (if it is not already set in RunSimulation)
  if (!SetRunNumber()) if (fStopOnError) return kFALSE;

  // Load and misalign the geometry
  if (!gGeoManager) {
    TGeoManager::Import("geometry.root");
    if (!gGeoManager) if (fStopOnError) return kFALSE;
    if (!MisalignGeometry()) if (fStopOnError) return kFALSE;
  }
  
  // hits -> summable digits
  if (!fMakeSDigits.IsNull()) {
    if (!RunSDigitization(fMakeSDigits)) if (fStopOnError) return kFALSE;
  }

  // summable digits -> digits
  if (!fMakeDigits.IsNull()) {
    if (!RunDigitization(fMakeDigits, fMakeDigitsFromHits)) {
      if (fStopOnError) return kFALSE;
    }
  }

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

  // digits -> trigger
  if (!RunTrigger(fMakeTrigger)) {
    if (fStopOnError) return kFALSE;
  }

  // digits -> raw data
  if (!fWriteRawData.IsNull()) {
    if (!WriteRawData(fWriteRawData, fRawDataFileName, 
		      fDeleteIntermediateFiles)) {
      if (fStopOnError) return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunTrigger(const char* descriptors)
{
  // run the trigger

   TStopwatch stopwatch;
   stopwatch.Start();

   AliRunLoader* runLoader = LoadRun("READ");
   if (!runLoader) return kFALSE;
   TString des = descriptors;

   if (des.IsNull()) {
     if (gAlice->GetTriggerDescriptor() != "") {
       des = gAlice->GetTriggerDescriptor();
     }
     else {
       AliWarning("No trigger descriptor is specified. Skipping the trigger simulation...");
       return kTRUE;
     }
   }

   runLoader->MakeTree( "CT" );
   AliCentralTrigger* aCTP = runLoader->GetTrigger();
  // Load Descriptors
   aCTP->LoadDescriptor( des );

  // digits -> trigger
   if( !aCTP->RunTrigger( runLoader ) ) {
      if (fStopOnError) {
    //  delete aCTP;
         return kFALSE;
      }
   }

   AliInfo(Form("Execution time: R:%.2fs C:%.2fs",
           stopwatch.RealTime(),stopwatch.CpuTime()));

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

  TStopwatch stopwatch;
  stopwatch.Start();

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
  StdoutToAliInfo(StderrToAliError(
    gAlice->Init(fConfigFileName.Data());
  ););

  // Get the trigger descriptor string
  // Either from AliSimulation or from
  // gAlice
  if (fMakeTrigger.IsNull()) {
    if (gAlice->GetTriggerDescriptor() != "")
      fMakeTrigger = gAlice->GetTriggerDescriptor();
  }
  else
    gAlice->SetTriggerDescriptor(fMakeTrigger.Data());

  // Set run number in CDBManager
  AliCDBManager::Instance()->SetRun(gAlice->GetRunNumber());
  AliInfo(Form("Run number: %d",AliCDBManager::Instance()->GetRun()));

  AliRunLoader* runLoader = gAlice->GetRunLoader();
  if (!runLoader) {
             AliError(Form("gAlice has no run loader object. "
        		     "Check your config file: %s", fConfigFileName.Data()));
             return kFALSE;
  }
  SetGAliceFile(runLoader->GetFileName());
 
  // Export ideal geometry 
  if (gGeoManager) gGeoManager->Export("geometry.root");

  // Misalign geometry
//   if (!MisalignGeometry(runLoader)) {
//     delete runLoader;
//     return kFALSE;
//   }
  MisalignGeometry(runLoader);

  // Export (mis)aligned geometry 
  if (gGeoManager) gGeoManager->Export("misaligned_geometry.root");

//   AliRunLoader* runLoader = gAlice->GetRunLoader();
//   if (!runLoader) {
//     AliError(Form("gAlice has no run loader object. "
//                   "Check your config file: %s", fConfigFileName.Data()));
//     return kFALSE;
//   }
//   SetGAliceFile(runLoader->GetFileName());

  if (!gAlice->Generator()) {
    AliError(Form("gAlice has no generator object. "
                  "Check your config file: %s", fConfigFileName.Data()));
    return kFALSE;
  }
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
    gAlice->Generator()->SetVertexGenerator(vtxGen);
  }

  if (!fRunSimulation) {
    gAlice->Generator()->SetTrackingFlag(0);
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
      AliError(Form("RunSimulation", "no loader for %s found\n"
                    "Number of events per file not set for %s %s", 
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
  StdoutToAliInfo(StderrToAliError(
    gAlice->Run(nEvents);
  ););

  delete runLoader;

  AliInfo(Form("Execution time: R:%.2fs C:%.2fs",
	       stopwatch.RealTime(),stopwatch.CpuTime()));

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunSDigitization(const char* detectors)
{
// run the digitization and produce summable digits

  TStopwatch stopwatch;
  stopwatch.Start();

  AliRunLoader* runLoader = LoadRun();
  if (!runLoader) return kFALSE;

  TString detStr = detectors;
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      AliInfo(Form("creating summable digits for %s", det->GetName()));
      TStopwatch stopwatchDet;
      stopwatchDet.Start();
      det->Hits2SDigits();
      AliInfo(Form("Execution time for %s: R:%.2fs C:%.2fs",
	   det->GetName(),stopwatchDet.RealTime(),stopwatchDet.CpuTime()));
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s",
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }

  delete runLoader;

  AliInfo(Form("Execution time: R:%.2fs C:%.2fs",
	   stopwatch.RealTime(),stopwatch.CpuTime()));

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliSimulation::RunDigitization(const char* detectors, 
				      const char* excludeDetectors)
{
// run the digitization and produce digits from sdigits

  TStopwatch stopwatch;
  stopwatch.Start();

  while (AliRunLoader::GetRunLoader()) delete AliRunLoader::GetRunLoader();
  if (gAlice) delete gAlice;
  gAlice = NULL;

  Int_t nStreams = 1;
  if (fBkgrdFileNames) nStreams = fBkgrdFileNames->GetEntriesFast() + 1;
  Int_t signalPerBkgrd = GetNSignalPerBkgrd();
  AliRunDigitizer* manager = new AliRunDigitizer(nStreams, signalPerBkgrd);
  // manager->SetEmbeddingFlag(fEmbeddingFlag);
  manager->SetInputStream(0, fGAliceFileName.Data());
  for (Int_t iStream = 1; iStream < nStreams; iStream++) {
    const char* fileName = ((TObjString*)
			    (fBkgrdFileNames->At(iStream-1)))->GetName();
    manager->SetInputStream(iStream, fileName);
  }

  TString detStr = detectors;
  TString detExcl = excludeDetectors;
  manager->GetInputStream(0)->ImportgAlice();
  AliRunLoader* runLoader = 
    AliRunLoader::GetRunLoader(manager->GetInputStream(0)->GetFolderName());
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr) && 
	!IsSelected(det->GetName(), detExcl)) {
      AliDigitizer* digitizer = det->CreateDigitizer(manager);
      
      if (!digitizer) {
	AliError(Form("no digitizer for %s", det->GetName()));
	if (fStopOnError) return kFALSE;
      } else {
	digitizer->SetRegionOfInterest(fRegionOfInterest);
      }
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s", 
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }

  if (!manager->GetListOfTasks()->IsEmpty()) {
    AliInfo("executing digitization");
    manager->Exec("");
  }

  delete manager;

  AliInfo(Form("Execution time: R:%.2fs C:%.2fs",
	       stopwatch.RealTime(),stopwatch.CpuTime()));
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunHitsDigitization(const char* detectors)
{
// run the digitization and produce digits from hits

  TStopwatch stopwatch;
  stopwatch.Start();

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

  delete runLoader;
  //PH Temporary fix to avoid interference with the PHOS loder/getter
  //PH The problem has to be solved in more general way 09/06/05

  AliInfo(Form("Execution time: R:%.2fs C:%.2fs",
	       stopwatch.RealTime(),stopwatch.CpuTime()));

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::WriteRawData(const char* detectors, 
				   const char* fileName,
				   Bool_t deleteIntermediateFiles)
{
// convert the digits to raw data
// First DDL raw data files for the given detectors are created.
// If a file name is given, the DDL files are then converted to a DATE file.
// If deleteIntermediateFiles is true, the DDL raw files are deleted 
// afterwards.
// If the file name has the extension ".root", the DATE file is converted
// to a root file.
// If deleteIntermediateFiles is true, the DATE file is deleted afterwards.

  TStopwatch stopwatch;
  stopwatch.Start();

  if (!WriteRawFiles(detectors)) {
    if (fStopOnError) return kFALSE;
  }

  TString dateFileName(fileName);
  if (!dateFileName.IsNull()) {
    Bool_t rootOutput = dateFileName.EndsWith(".root");
    if (rootOutput) dateFileName += ".date";
    if (!ConvertRawFilesToDate(dateFileName)) {
      if (fStopOnError) return kFALSE;
    }
    if (deleteIntermediateFiles) {
      AliRunLoader* runLoader = LoadRun("READ");
      if (runLoader) for (Int_t iEvent = 0; 
			  iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
	char command[256];
	sprintf(command, "rm -r raw%d", iEvent);
	gSystem->Exec(command);
      }
    }

    if (rootOutput) {
      if (!ConvertDateToRoot(dateFileName, fileName)) {
	if (fStopOnError) return kFALSE;
      }
      if (deleteIntermediateFiles) {
	gSystem->Unlink(dateFileName);
      }
    }
  }

  AliInfo(Form("Execution time: R:%.2fs C:%.2fs",
	       stopwatch.RealTime(),stopwatch.CpuTime()));

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::WriteRawFiles(const char* detectors)
{
// convert the digits to raw data DDL files

  AliRunLoader* runLoader = LoadRun("READ");
  if (!runLoader) return kFALSE;

  // write raw data to DDL files
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    AliInfo(Form("processing event %d", iEvent));
    runLoader->GetEvent(iEvent);
    TString baseDir = gSystem->WorkingDirectory();
    char dirName[256];
    sprintf(dirName, "raw%d", iEvent);
    gSystem->MakeDirectory(dirName);
    if (!gSystem->ChangeDirectory(dirName)) {
      AliError(Form("couldn't change to directory %s", dirName));
      if (fStopOnError) return kFALSE; else continue;
    }

    TString detStr = detectors;
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
Bool_t AliSimulation::ConvertRawFilesToDate(const char* dateFileName)
{
// convert raw data DDL files to a DATE file with the program "dateStream"

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
  char command[256];
  // Note the option -s. It is used in order to avoid
  // the generation of SOR/EOR events.
  sprintf(command, "dateStream -s -D -o %s -# %d -C", 
	  dateFileName, runLoader->GetNumberOfEvents());
  FILE* pipe = gSystem->OpenPipe(command, "w");

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    fprintf(pipe, "GDC\n");
    Float_t ldc = 0;
    Int_t prevLDC = -1;

    // loop over detectors and DDLs
    for (Int_t iDet = 0; iDet < AliDAQ::kNDetectors; iDet++) {
      for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls(iDet); iDDL++) {

        Int_t ddlID = AliDAQ::DdlID(iDet,iDDL);
        Int_t ldcID = Int_t(ldc + 0.0001);
        ldc += AliDAQ::NumberOfLdcs(iDet) / AliDAQ::NumberOfDdls(iDet);

        char rawFileName[256];
        sprintf(rawFileName, "raw%d/%s", 
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

  delete runLoader;
  return (result == 0);
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
  const Int_t kCompression = 0;

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
  const char* runDBFS    = "/tmp/mdc1/meta";

  // User defined file system locations
  if (gSystem->Getenv("ALIMDC_RAWDB1")) 
    rawDBFS[0] = gSystem->Getenv("ALIMDC_RAWDB1");
  if (gSystem->Getenv("ALIMDC_RAWDB2")) 
    rawDBFS[1] = gSystem->Getenv("ALIMDC_RAWDB2");
  if (gSystem->Getenv("ALIMDC_TAGDB")) 
    tagDBFS = gSystem->Getenv("ALIMDC_TAGDB");
  if (gSystem->Getenv("ALIMDC_RUNDB")) 
    runDBFS = gSystem->Getenv("ALIMDC_RUNDB");

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));
  gSystem->Exec(Form("rm -rf %s",runDBFS));

  gSystem->Exec(Form("mkdir %s",rawDBFS[0]));
  gSystem->Exec(Form("mkdir %s",rawDBFS[1]));
  gSystem->Exec(Form("mkdir %s",tagDBFS));
  gSystem->Exec(Form("mkdir %s",runDBFS));

  Int_t result = gSystem->Exec(Form("alimdc %d %d %d %d %s", 
				    kDBSize, kTagDBSize, kFilter, kCompression, dateFileName));
  gSystem->Exec(Form("mv %s/*.root %s", rawDBFS[0], rootFileName));

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));
  gSystem->Exec(Form("rm -rf %s",runDBFS));

  return (result == 0);
}


//_____________________________________________________________________________
AliRunLoader* AliSimulation::LoadRun(const char* mode) const
{
// delete existing run loaders, open a new one and load gAlice

  while (AliRunLoader::GetRunLoader()) delete AliRunLoader::GetRunLoader();
  AliRunLoader* runLoader = 
    AliRunLoader::Open(fGAliceFileName.Data(), 
		       AliConfig::GetDefaultEventFolderName(), mode);
  if (!runLoader) {
    AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
    return NULL;
  }
  runLoader->LoadgAlice();
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

Bool_t AliSimulation::ConvertRaw2SDigits(const char* rawDirectory, const char* esdFileName) 
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
    StdoutToAliInfo(StderrToAliError(gAlice->Init(fConfigFileName.Data());););
//
//  Initialize CDB     
    InitCDBStorage();
    AliCDBManager* man = AliCDBManager::Instance();
    man->SetRun(0); // Should this come from rawdata header ?
    
    Int_t iDet;
    //
    // Get the runloader
    AliRunLoader* runLoader = gAlice->GetRunLoader();
    //
    // Open esd file if available
    TFile* esdFile = TFile::Open(esdFileName);
    Bool_t esdOK = (esdFile != 0);
    AliESD* esd = new AliESD;
    TTree* treeESD = 0;
    if (esdOK) {
	treeESD = (TTree*) esdFile->Get("esdTree");
	if (!treeESD) {
	    AliWarning("No ESD tree found");
	    esdOK = kFALSE;
	} else {
	    treeESD->SetBranchAddress("ESD", &esd);
	}
    }
    //
    // Create the RawReader
    AliRawReaderFile* rawReader = new AliRawReaderFile(rawDirectory);
    //
    // Get list of detectors
    TObjArray* detArray = runLoader->GetAliRun()->Detectors();
    //
    // Get Header
    AliHeader* header = runLoader->GetHeader();
    //
    // Event loop
    Int_t nev = 0;
    while(kTRUE) {
	if (!(rawReader->NextEvent())) break;
	//
	// Detector loop
	for (iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
	    AliModule* det = (AliModule*) detArray->At(iDet);
	    AliInfo(Form("Calling Raw2SDigits for %s\n", det->GetName()));
	    det->Raw2SDigits(rawReader);
	    rawReader->Reset();
	} // detectors

	//
	//  If ESD information available obtain reconstructed vertex and store in header.
	if (esdOK) {
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
	    printf("***** Saved vertex %f %f %f \n", position[0], position[1], position[2]);
	}
	nev++;
//
//      Finish the event
	runLoader->TreeE()->Fill();
	runLoader->SetNextEvent();
    } // events
 
    delete rawReader;
//
//  Finish the run 
    runLoader->CdGAFile();
    runLoader->WriteHeader("OVERWRITE");
    runLoader->WriteRunLoader();

    return kTRUE;
}
