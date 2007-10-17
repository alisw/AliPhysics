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

#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"

#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include <TSystem.h>
#include <TClass.h>

/// \cond CLASSIMP
ClassImp(AliMpCDB)
/// \endcond

//
// static methods
//

//______________________________________________________________________________
Bool_t AliMpCDB::LoadMpSegmentation(Bool_t warn)
{
/// Load the sementation from the CDB if it does not yet exist;
/// return false only in case loading from CDB failed

  AliDebugClass(1,"");

  if ( AliMpSegmentation::Instance(false) ) {
    if ( warn )  
      AliWarningClass("Segmentation has been already loaded."); 
    return true;
  }  
  
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT");

  Int_t run = cdbManager->GetRun();
  if ( run < 0 ) {
    run = 0;
    cdbManager->SetRun(run);
  }  

  AliCDBEntry* cdbEntry = cdbManager->Get("MUON/Calib/Mapping", run);
  
  if ( cdbEntry ) 
  {
    return (AliMpSegmentation*)cdbEntry->GetObject() != 0x0;
  }
  else
  {
    return kFALSE;
  }
}    

//______________________________________________________________________________
Bool_t AliMpCDB::LoadDDLStore(Bool_t warn)
{
/// Load the DDL store from the CDB if it does not yet exist
/// return false only in case loading from CDB failed

  AliDebugClass(1,"");

  if ( AliMpDDLStore::Instance(false) ) {
    if ( warn )  
      AliWarningClass("DDL Store has been already loaded."); 
    return true;
  }  
  
  // Load segmentation
  //
  LoadMpSegmentation(warn); 
  
  // Load DDL store
  //
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT");

  Int_t run = cdbManager->GetRun();
  if ( run < 0 ) {
    run = 0;
    cdbManager->SetRun(run);
  }  

  AliCDBEntry* cdbEntry =  cdbManager->Get("MUON/Calib/DDLStore", run);
  
  if ( cdbEntry ) 
  {
    return (AliMpDDLStore*)cdbEntry->GetObject() != 0x0;
  }
  else
  {
    return kFALSE;
  }
}    

//______________________________________________________________________________
Bool_t AliMpCDB::WriteMpSegmentation()
{
/// Write mapping segmentation in OCDB

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON mapping");
  cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("MUON/Calib/Mapping", 0, 9999999); 

  AliMpSegmentation::ReadData(false);
  return cdbManager->Put(AliMpSegmentation::Instance(), id, cdbData);
}

//______________________________________________________________________________
Bool_t AliMpCDB::WriteDDLStore()
{
/// Write mapping DDL store in OCDB

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if ( ! cdbManager->GetDefaultStorage() )
    cdbManager->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON DDL store");
  cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("MUON/Calib/DDLStore", 0, 9999999); 

  AliMpSegmentation::ReadData(false);
  AliMpDDLStore::ReadData(false);
  return cdbManager->Put(AliMpDDLStore::Instance(), id, cdbData);
}   
