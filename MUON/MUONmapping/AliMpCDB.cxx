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

// $Id$
// $MpId: $
// Category: management

//-----------------------------------------------------------------------------
// Class AliMpCDB
// -----------------------
// Manager class for mapping CDB IO
// Author: Ivana Hrivnacova, IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpCDB.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"
#include "AliMpDataMap.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataStreams.h"
#include "AliMpManuStore.h"
#include "AliMpSegmentation.h"
#include <Riostream.h>
#include <TClass.h>
#include <TSystem.h>

/// \cond CLASSIMP
ClassImp(AliMpCDB)
/// \endcond

Bool_t AliMpCDB::fgLoadFromData = kTRUE;                                       

//
// private static methods
//

//______________________________________________________________________________
TObject*  AliMpCDB::GetCDBEntryObject(const char* dataPath)
{
/// Load CDB entry object with checks

  AliCDBManager* cdbManager = AliCDBManager::Instance();

  Int_t run = cdbManager->GetRun();
  if ( run < 0 ) {
    AliErrorClassStream() << "Cannot get run number from CDB manager." << endl; 
    return 0;
  }  

  AliCDBEntry* cdbEntry = cdbManager->Get(dataPath, run);
  if ( ! cdbEntry ) {
    AliErrorClassStream() << "Cannot get cdbEntry." << endl; 
    return 0;
  }
  
  TObject* object = cdbEntry->GetObject();
  if ( ! object ) {
    AliErrorClassStream() << "Cannot get object from cdbEntry." << endl; 
    return 0;
  }  

  return object;
}    
  

//______________________________________________________________________________
TObject*  AliMpCDB::GetCDBEntryObject(const char* dataPath, 
                                      const char* cdbpath,
                                      Int_t runNumber )
{
/// Load CDB entry from CDB and run specified in arguments

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage(cdbpath);

  AliCDBEntry* cdbEntry = cdbManager->Get(dataPath, runNumber);
  if ( ! cdbEntry ) {
    AliErrorClassStream() << "Cannot get cdbEntry." << endl; 
    return 0;
  }  
    
  TObject* object = cdbEntry->GetObject();
  if ( ! object ) {
    AliErrorClassStream() << "Cannot get object from cdbEntry." << endl; 
    return 0;
  }  

  return object;
}    
                                      
//
// public static methods
//

   
//______________________________________________________________________________
Bool_t AliMpCDB::LoadMpSegmentation(Bool_t warn)
{
/// Load the sementation from the mapping data from OCDB,
///  if it does not yet exist;
/// return false only in case loading from CDB failed

  if ( AliMpSegmentation::Instance(false) ) {
    if ( warn )  
      AliWarningClass("Segmentation has been already loaded."); 
    return true;
  }  
  
  if ( fgLoadFromData ) {
    AliDebugClassStream(1)
      << "Loading segmentation from MUON/Calib/MappingData" << endl;
  
    TObject* cdbEntryObject = GetCDBEntryObject("MUON/Calib/MappingData");
    if ( ! cdbEntryObject ) return kFALSE;
  
    // Pass the map to the streams and then read segmentation
    // from data map
    AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
    AliMpDataStreams dataStreams(dataMap);
    AliMpSegmentation::ReadData(dataStreams);
    return kTRUE;
  }
  else {
    AliDebugClassStream(1)
      << "Loading segmentation from MUON/Calib/Mapping" << endl;
  
    TObject* cdbEntryObject = GetCDBEntryObject("MUON/Calib/Mapping");
    return cdbEntryObject != 0x0;
  }  
}    

//______________________________________________________________________________
Bool_t AliMpCDB::LoadDDLStore(Bool_t warn)
{
/// Load the DDL store from the mapping data from OCDB,
///  if it does not yet exist;
/// return false only in case loading from CDB failed

  if ( AliMpDDLStore::Instance(false) ) {
    if ( warn )  
      AliWarningClass("DDL Store has been already loaded."); 
    return true;
  }  
  
  if ( fgLoadFromData ) {
    AliDebugClassStream(1)
      << "Loading DDL store from MUON/Calib/MappingData" << endl;
  
    TObject* cdbEntryObject = GetCDBEntryObject("MUON/Calib/MappingData");
    if ( ! cdbEntryObject ) return kFALSE;
  
    AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
    AliMpDataStreams dataStreams(dataMap);
    AliMpDDLStore::ReadData(dataStreams);
    return kTRUE;
  }
  else {
    AliDebugClassStream(1)
      << "Loading DDL store from MUON/Calib/DDLStore" << endl;
  
    // Load segmentation
    LoadMpSegmentation(warn); 
  
    // Load DDL store
    TObject* cdbEntryObject =  GetCDBEntryObject("MUON/Calib/DDLStore");
    return cdbEntryObject != 0x0;
  }     
}    

//______________________________________________________________________________
Bool_t AliMpCDB::LoadAll2(const char* cdbpath, Int_t runNumber, Bool_t warn)
{
  /// Load everything in one shot 
  return 
  LoadDDLStore2(cdbpath,runNumber,warn) && 
  LoadManuStore2(cdbpath,runNumber,warn);
}

//______________________________________________________________________________
Bool_t AliMpCDB::LoadAll(Bool_t warn)
{
  /// Load everything in one shot 
  return LoadDDLStore(warn) && LoadManuStore(warn);
}

//______________________________________________________________________________
Bool_t AliMpCDB::LoadManuStore(Bool_t warn)
{
/// Load the DDL store from the mapping data from OCDB,
///  if it does not yet exist;
/// return false only in case loading from CDB failed

  if ( AliMpManuStore::Instance(false) ) {
    if ( warn )  
      AliWarningClass("Manu Store has been already loaded."); 
    return true;
  }  
  
  if ( fgLoadFromData ) {
    AliDebugClassStream(1)
      << "Loading Manu store from MUON/Calib/MappingRunData" << endl;
  
    // Load segmentation
    if ( ! AliMpSegmentation::Instance(false) ) 
      LoadMpSegmentation(warn); 

    TObject* cdbEntryObject = GetCDBEntryObject("MUON/Calib/MappingRunData");
    if ( ! cdbEntryObject ) return kFALSE;
  
    AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
    AliMpDataStreams dataStreams(dataMap);
    AliMpManuStore::ReadData(dataStreams);
    return kTRUE;
  }
  else {
    AliDebugClassStream(1)
      << "Loading Manu store from MUON/Calib/ManuStore" << endl;
  
    // Load Manu store
    TObject* cdbEntryObject =  GetCDBEntryObject("MUON/Calib/ManuStore");
    return cdbEntryObject != 0x0;
  }     
}    
//______________________________________________________________________________
Bool_t AliMpCDB::LoadMpSegmentation2(const char* cdbpath, Int_t runNumber, 
                                     Bool_t warn)
{
/// Load the sementation from the CDB if it does not yet exist;
/// return false only in case loading from CDB failed.
/// In difference from LoadMpSegmentation(), in this method the CDB path
/// and run is set directly via arguments.


  if ( AliMpSegmentation::Instance(false) ) {
    if ( warn )  
      AliWarningClass("Segmentation has been already loaded."); 
    return true;
  }  
  
  if ( fgLoadFromData ) {
    AliDebugClassStream(1)
      << "Loading segmentation from MUON/Calib/MappingData" << endl;
  
    TObject* cdbEntryObject 
      = GetCDBEntryObject("MUON/Calib/MappingData", cdbpath, runNumber);
    if ( ! cdbEntryObject ) return kFALSE;
  
    // Pass the map to the streams and then read segmentation
    // from data map
    AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
    AliMpDataStreams dataStreams(dataMap);
    AliMpSegmentation::ReadData(dataStreams);
    return kTRUE;
  }
  else {
    AliDebugClassStream(1)
      << "Loading segmentation from MUON/Calib/Mapping" << endl;
  
    TObject* cdbEntryObject 
      = GetCDBEntryObject("MUON/Calib/Mapping", cdbpath, runNumber);
    return cdbEntryObject != 0x0;
  }  
}    

//______________________________________________________________________________
Bool_t AliMpCDB::LoadDDLStore2(const char* cdbpath, Int_t runNumber, 
                               Bool_t warn)
{
/// Load the DDL store from the CDB if it does not yet exist
/// return false only in case loading from CDB failed
/// In difference from LoadDDLStore(), in this method the CDB path
/// and run is set directly via arguments.

  if ( AliMpDDLStore::Instance(false) ) {
    if ( warn )  
      AliWarningClass("DDL Store has been already loaded."); 
    return true;
  }  
  
  if ( fgLoadFromData ) {
    AliDebugClassStream(1)
      << "Loading DDL store from MUON/Calib/MappingData" << endl;
  
    TObject* cdbEntryObject 
      = GetCDBEntryObject("MUON/Calib/MappingData", cdbpath, runNumber);
    if ( ! cdbEntryObject ) return kFALSE;
  
    AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
    AliMpDataStreams dataStreams(dataMap);
    AliMpDDLStore::ReadData(dataStreams);
    return kTRUE;
  }
  else {
    AliDebugClassStream(1)
      << "Loading DDL store from MUON/Calib/DDLStore" << endl;
  
    // Load segmentation
    LoadMpSegmentation2(cdbpath, runNumber, warn); 
  
    // Load DDL store
    TObject* cdbEntryObject =  GetCDBEntryObject("MUON/Calib/DDLStore");
    return cdbEntryObject != 0x0;
  }  
}    

//______________________________________________________________________________
Bool_t AliMpCDB::LoadManuStore2(const char* cdbpath, Int_t runNumber, 
                                Bool_t warn)
{
/// Load the DDL store from the CDB if it does not yet exist
/// return false only in case loading from CDB failed
/// In difference from LoadDDLStore(), in this method the CDB path
/// and run is set directly via arguments.

  if ( AliMpManuStore::Instance(false) ) {
    if ( warn )  
      AliWarningClass("Manu Store has been already loaded."); 
    return true;
  }  
  
  if ( fgLoadFromData ) {
    AliDebugClassStream(1)
      << "Loading Manu store from MUON/Calib/MappingRunData" << endl;
  
    // Load segmentation
    LoadMpSegmentation2(cdbpath, runNumber, warn); 

    TObject* cdbEntryObject 
      = GetCDBEntryObject("MUON/Calib/MappingRunData", cdbpath, runNumber);
    if ( ! cdbEntryObject ) return kFALSE;
  
    AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
    AliMpDataStreams dataStreams(dataMap);
    AliMpManuStore::ReadData(dataStreams);
    return kTRUE;
  }
  else {
    AliDebugClassStream(1)
      << "Loading Manu store from MUON/Calib/ManuStore" << endl;
  
    // Load Manu store
    TObject* cdbEntryObject =  GetCDBEntryObject("MUON/Calib/ManuStore");
    return cdbEntryObject != 0x0;
  }  
}    

//______________________________________________________________________________
Bool_t AliMpCDB::WriteMpData()
{
/// Write mapping data in OCDB

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON mapping");
  cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("MUON/Calib/MappingData", 0, AliCDBRunRange::Infinity()); 

  AliMpDataProcessor mp;
  AliMpDataMap* map = mp.CreateDataMap("data");
  return cdbManager->Put(map, id, cdbData);
}

//______________________________________________________________________________
Bool_t AliMpCDB::WriteMpRunData()
{
/// Write mapping data in OCDB

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON run-dependent mapping");
  cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("MUON/Calib/MappingRunData", 0, AliCDBRunRange::Infinity()); 

  AliMpDataProcessor mp;
  AliMpDataMap* map = mp.CreateDataMap("data_run");
  return cdbManager->Put(map, id, cdbData);
}

//______________________________________________________________________________
Bool_t AliMpCDB::WriteMpSegmentation(Bool_t readData)
{
/// Write mapping segmentation in OCDB

  if ( ! readData && ! AliMpSegmentation::Instance() ) return false;

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON mapping");
  cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("MUON/Calib/Mapping", 0, AliCDBRunRange::Infinity()); 

  if ( readData ) {
    AliMpDataStreams dataStreams;
    AliMpSegmentation::ReadData(dataStreams, false);
    AliMpDDLStore::ReadData(dataStreams, false);
  }
  
  return cdbManager->Put(AliMpSegmentation::Instance(), id, cdbData);
}

//______________________________________________________________________________
Bool_t AliMpCDB::WriteDDLStore(Bool_t readData)
{
/// Write mapping DDL store in OCDB

  if ( ! readData && ! AliMpDDLStore::Instance() ) return false;

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON DDL store");
  cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("MUON/Calib/DDLStore", 0, AliCDBRunRange::Infinity()); 

  if ( readData ) {
    AliMpDataStreams dataStreams;
    AliMpSegmentation::ReadData(dataStreams, false);
    AliMpDDLStore::ReadData(dataStreams, false);
  }
  return cdbManager->Put(AliMpDDLStore::Instance(), id, cdbData);
}   

//______________________________________________________________________________
Bool_t AliMpCDB::WriteManuStore(Bool_t readData)
{
/// Write mapping Manu store in OCDB

  if ( ! readData && ! AliMpManuStore::Instance() ) return false;

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON Manu store");
  cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("MUON/Calib/ManuStore", 0, AliCDBRunRange::Infinity()); 

  if ( readData ) {
    AliMpDataStreams dataStreams;
    AliMpSegmentation::ReadData(dataStreams, false);
    AliMpManuStore::ReadData(dataStreams, false);
  }
  return cdbManager->Put(AliMpManuStore::Instance(), id, cdbData);
}   

//______________________________________________________________________________
Bool_t  AliMpCDB::GenerateMpData(const char* cdbpath, Int_t runNumber)
{
/// Generate mapping data ASCII files from OCDB

  TObject* cdbEntryObject 
    = GetCDBEntryObject("MUON/Calib/MappingData", cdbpath, runNumber);
  if ( ! cdbEntryObject ) return kFALSE;
  
  AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
  AliMpDataProcessor mp;
  return mp.GenerateData(dataMap);
} 

//______________________________________________________________________________
Bool_t  AliMpCDB::GenerateMpRunData(const char* cdbpath, Int_t runNumber)
{
/// Generate mapping data ASCII files from OCDB

  TObject* cdbEntryObject 
    = GetCDBEntryObject("MUON/Calib/MappingRunData", cdbpath, runNumber);
  if ( ! cdbEntryObject ) return kFALSE;
  
  AliMpDataMap* dataMap = (AliMpDataMap*)cdbEntryObject;
  AliMpDataProcessor mp;
  return mp.GenerateData(dataMap);
} 

//______________________________________________________________________________
void AliMpCDB::UnloadAll()
{
  /// Unload all the mapping from the memory
  delete AliMpDDLStore::Instance(false);
  delete AliMpSegmentation::Instance(false);
  delete AliMpDEStore::Instance(false);
  delete AliMpManuStore::Instance(false);
}

