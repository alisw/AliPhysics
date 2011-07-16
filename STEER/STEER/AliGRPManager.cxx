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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliGRPManager class                                                    //
// The class can be used in order to access and read the Global Run       //
// Parameters entry from OCDB.                                            //
// It has a methods to set the magnetic field instanton and return        //
// the run and event info objects.                                        //
//                                                                        //
// cvetan.cheshkov@cern.ch 15/06/2009                                     //
//                                                                        //
// Usage:                                                                 //
// AliGRPManager grpMan;                                                  //
// Bool_t status = kTRUE;                                                 //
// status = grpMan.ReadGRPEntry(); // Read the corresponding OCDB entry   //
// status = grpMan.SetMagField();  // Set global field instanton          //
// AliRunInfo *runInfo = grpMan.GetRunInfo();// Get instance of run info  //
//                                                                        //
// Note: CDB manager should be initialized beforehand                     //
////////////////////////////////////////////////////////////////////////////

#include <TGeoGlobalMagField.h>

#include "AliGRPManager.h"
#include "AliLog.h"
#include "AliRunInfo.h"
#include "AliGRPObject.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliMagF.h"

ClassImp(AliGRPManager)

//_____________________________________________________________________________
AliGRPManager::AliGRPManager() :
  TObject(),
  fGRPData(NULL)
{
  // Default constructor
}

//_____________________________________________________________________________
AliGRPManager::~AliGRPManager()
{
  // Destructor
  if (fGRPData) delete fGRPData;
}

//_____________________________________________________________________________
Bool_t AliGRPManager::ReadGRPEntry()
{
  //------------------------------------
  // Initialization of the GRP entry. 
  // Returns kTRUE in case of success.
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

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliGRPManager::SetMagField()
{
  // Dealing with the magnetic field map
  // Construct the mag field map from the data in GRP
  // Set the global mag field instance

  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      AliInfo("ExpertMode!!! GRP information will be ignored !");
      AliInfo("ExpertMode!!! Running with the externally locked B field !");
      return kTRUE;
    }
    else {
      AliInfo("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }
  }
  //
  if (!fGRPData) {
    AliError("GRP Data is not loaded");
    return kFALSE;
  }
  //
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
  
  TString beamType = fGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
    //ok = kFALSE;  // temprorary suppressed to make read cosmics data
  }
  
  Float_t beamEnergy = fGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
    //ok = kFALSE;  // temprorary suppressed to make read cosmics data
  }
  
  // read special bits for the polarity convention and map type
  Int_t  polConvention = fGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
  Bool_t uniformB = fGRPData->IsUniformBMap();
  
  if (ok) { 
    AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1), 
					   TMath::Abs(diCurrent) * (diPolarity ? -1:1), 
					   polConvention,uniformB,beamEnergy, beamType.Data());
    if (fld) {
      TGeoGlobalMagField::Instance()->SetField( fld );
      TGeoGlobalMagField::Instance()->Lock();
      AliInfo("Running with the B field constructed out of GRP !");
    }
    else {
      AliError("Failed to create a B field map !");
      ok = kFALSE;
    }
  }
  else {
    AliError("B field is neither set nor constructed from GRP ! Exitig...");
  }
  
  return ok;
}

//_____________________________________________________________________________
AliRunInfo* AliGRPManager::GetRunInfo()
{
  // Constructs and returns an object
  // containing the run information
  // The user code is the owner of the object

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

  return new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);
}

//_____________________________________________________________________________
void AliGRPManager::SetGRPEntry(AliGRPObject* source)
{
  // Create a GRP entry from the extrnaly provide GRP object
  // To be used by HLT to create an online GRP instance
  if (!source) return;
  if (fGRPData) delete fGRPData;
  fGRPData = new AliGRPObject(*source);
  AliInfo("Created GRP Data from external object");
  //
}
 
