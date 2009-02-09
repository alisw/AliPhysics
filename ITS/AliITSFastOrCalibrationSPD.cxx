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

#include "AliITSFastOrCalibrationSPD.h"
///////////////////////////////////////////////////////////////////////////
//
//  Calibration class for the SPD FastOr configuration info
//  
//  C. Di Giglio Carmelo.Digiglio@ba.infn.it
//  D. Elia      Domenico.Elia@ba.infn.it 
//
///////////////////////////////////////////////////////////////////////////

ClassImp(AliITSFastOrCalibrationSPD)

//-----------------------------------------------
//Deafault constructor
AliITSFastOrCalibrationSPD::AliITSFastOrCalibrationSPD():
TObject(),
fFastOrConfiguredChips(1200)
{
// constructor
}

//Default destructor
AliITSFastOrCalibrationSPD::~AliITSFastOrCalibrationSPD() {}
//____________________________________________________________________________________________
Bool_t AliITSFastOrCalibrationSPD::WriteFOConfToDB(Int_t runNrStart, Int_t runNrEnd)  {

  AliCDBManager* man = AliCDBManager::Instance();

  if(!man->IsDefaultStorageSet()) {
     man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }

  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Domenico Elia");
  metaData->SetComment("Created by storeFastOrConfToDB.C");
  AliCDBId idCalSPD("ITS/Calib/SPDFastOr",runNrStart,runNrEnd);
  AliCDBEntry* cdbEntry = new AliCDBEntry(this,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete cdbEntry;
  delete metaData;

  return kTRUE;
}
