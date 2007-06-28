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

// Macro for generating the full misalignment data.
// The macro is trigger from AliRoot/macros/MakeAllDETsFullMisAlignment.C
//
// Author: I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryMisAligner.h"

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


void MakeMUONFullMisAlignment()
{
  // Load geometry, if not yet loaded,
  if ( ! AliGeomManager::GetGeometry() )
    AliGeomManager::LoadGeometry("geometry.root");

  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData();
  
  AliMUONGeometryMisAligner misAligner(0.0, 0.03, 0.0, 0.03, 0.0, 0.03);
  AliMUONGeometryTransformer* newTransform 
    = misAligner.MisAlign(&transformer, true);
  const TClonesArray* array = newTransform->GetMisAlignmentData();
  
  if ( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ) {
    cout << "Generating full misalignment data in a file" << endl;

    // Create a File to store the alignment data
    TFile f("MUONfullMisalignment.root","RECREATE");
    if ( !f.IsOpen() ) {
      cerr << "cannot open file for output" << endl;
    }
    
    f.cd();
    f.WriteObject(array,"MUONAlignObjs ","kSingleKey");
    f.Close();
  }
  else {
    cout << "Generating full misalignment data in CDB" << endl;

    // save in CDB storage
    const char* Storage = gSystem->Getenv("STORAGE");
    
    AliCDBManager* cdbManager = AliCDBManager::Instance();
    AliCDBStorage* storage = cdbManager->GetStorage(Storage);
    AliCDBMetaData* cdbData = new AliCDBMetaData();
    cdbData->SetResponsible("Dimuon Offline project");
    cdbData->SetComment("MUON alignment objects with full misalignment");
    cdbData->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("MUON/Align/Data", 0, 9999999); 
    storage->Put(const_cast<TClonesArray*>(array), id, cdbData);
  }   
  delete newTransform;
}
