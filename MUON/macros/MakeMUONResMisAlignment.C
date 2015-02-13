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

/// \ingroup macros
/// \file MakeMUONResMisAlignment.C
/// \brief Macro for generating the residual misalignment data.
///
/// The macro is triggered from AliRoot/macros/MakeAllDETsResMisAlignment.C
///
/// \author: I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryMisAligner.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"

#include <TSystem.h>
#include <TError.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TFile.h>
#include <Riostream.h>

#endif

void MakeMUONResMisAlignment()
{
  const char* macroname = "MakeMUONResMisAlignment.C";
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  AliCDBStorage* storage = 0;
  
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
  }else{
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage
  }    

  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData();
  
  AliMUONGeometryMisAligner misAligner(0.0, 0.004, 0.0, 0.003, 0.0, 0.0023);
  AliMUONGeometryTransformer* newTransform 
    = misAligner.MisAlign(&transformer, true);
  const TClonesArray* array = newTransform->GetMisAlignmentData();

  // 100 mum residual resolution for chamber misalignments?
  misAligner.SetAlignmentResolution(array,-1,0.01,0.01);

  if ( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ) {
    // Save in file
    const char* filename = "MUONresidualMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f.IsOpen()){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"MUONAlignObjs","kSingleKey");
    f.Close();
  } else {
    // Save in CDB storage
    AliCDBMetaData* cdbData = new AliCDBMetaData();
    cdbData->SetResponsible("Dimuon Offline project");
    cdbData->SetComment("MUON alignment objects with residual misalignment");
    cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("MUON/Align/Data", 0, AliCDBRunRange::Infinity()); 
    storage->Put(const_cast<TClonesArray*>(array), id, cdbData);
  }
  delete newTransform;
 }   

