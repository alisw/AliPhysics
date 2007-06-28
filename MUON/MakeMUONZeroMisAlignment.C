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

// Macro for generating the zero misalignment data.
//
// Author: I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONGeometryTransformer.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"

#include <TSystem.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TFile.h>
#include <Riostream.h>

#endif

void MakeMUONZeroMisAlignment()
{
  // Load geometry, if not yet loaded,
  if ( ! AliGeomManager::GetGeometry() )
    AliGeomManager::LoadGeometry("geometry.root");

  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData();
  TClonesArray* array = transformer.CreateZeroAlignmentData();;

  if ( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ) {
    // save in file
    cout << "Generating zero misalignment data in a file ..." << endl;

    // Create a file to store the alignment data
    TFile f("MUONzeroMisalignment.root", "RECREATE");
    if( !f.IsOpen() ) {
      cerr<<"cannot open file for output\n";
    }
    
    f.cd();
    f.WriteObject(array,"MUONAlignObjs ","kSingleKey");
    f.Close();
  }
  else {
    cout << "Generating zero misalignment data in CDB ..." << endl;

    // save in CDB storage
    const char* Storage = gSystem->Getenv("STORAGE");
    AliCDBManager* cdbManager = AliCDBManager::Instance();
    AliCDBStorage* storage = cdbManager->GetStorage(Storage);
    AliCDBMetaData* cdbData = new AliCDBMetaData();
    cdbData->SetResponsible("Dimuon Offline project");
    cdbData->SetComment("MUON alignment objects with zero misalignment");
    cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("MUON/Align/Data", 0, 9999999); 
    storage->Put(array, id, cdbData);
  }
}   

