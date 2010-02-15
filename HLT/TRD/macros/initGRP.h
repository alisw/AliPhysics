#if !defined (__CINT__) || defined (__MAKECINT__)

#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#endif

Bool_t InitGRP(TString OCDBpath, TString GRPpath="") {
  //------------------------------------
  // Initialization of the GRP entry 
  //------------------------------------

  AliCDBManager::Instance()->SetDefaultStorage(OCDBpath.Data());
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data", GRPpath.Data());
  AliCDBManager::Instance()->SetRun(0);

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* GRPData;

  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
      printf("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
      m->Print();
      GRPData = new AliGRPObject();
      GRPData->ReadValuesFromMap(m);
    }

    else {
      printf("Found an AliGRPObject in GRP/GRP/Data, reading it");
      GRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
      entry->SetOwner(0);
    }

    //    FIX ME: The unloading of GRP entry is temporarily disabled
    //    because ZDC and VZERO are using it in order to initialize
    //    their reconstructor objects. In the future one has to think
    //    of propagating AliRunInfo to the reconstructors.
    //    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if (!GRPData) {
    printf("No GRP entry found in OCDB!");
    return kFALSE;
  }

  TString lhcState = GRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    printf("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = GRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    printf("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = GRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    printf("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }
  // LHC: "multiply by 120 to get the energy in MeV"
  beamEnergy *= 0.120;

  TString runType = GRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    printf("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = GRPData->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidUInt()) {
    printf("GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }

  //*** Dealing with the magnetic field map
  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      printf("ExpertMode!!! GRP information will be ignored !");
      printf("ExpertMode!!! Running with the externally locked B field !");
    }
    else {
      printf("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }    
  }
  if ( !TGeoGlobalMagField::Instance()->IsLocked() ) {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = GRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      printf("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = GRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      printf("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = GRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      printf("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = GRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      printf("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    // read special bits for the polarity convention and map type
    Int_t  polConvention = GRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    Bool_t uniformB = GRPData->IsUniformBMap();

    if (ok) { 
      AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1), 
					     TMath::Abs(diCurrent) * (diPolarity ? -1:1), 
					     polConvention,uniformB,beamEnergy, beamType.Data());
      if (fld) {
	TGeoGlobalMagField::Instance()->SetField( fld );
	TGeoGlobalMagField::Instance()->Lock();
	printf("Running with the B field constructed out of GRP !");
      }
      else printf("Failed to create a B field map !");
    }
    else printf("B field is neither set nor constructed from GRP ! Exitig...");
  }
  
  //*** Get the diamond profiles from OCDB
  // entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexSPD");
  // if (entry) {
  //   fDiamondProfileSPD = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  // } else {
  //   printf("No SPD diamond profile found in OCDB!");
  // }

  // entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertex");
  // if (entry) {
  //   fDiamondProfile = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  // } else {
  //   printf("No diamond profile found in OCDB!");
  // }

  // entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexTPC");
  // if (entry) {
  //   fDiamondProfileTPC = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  // } else {
  //   printf("No TPC diamond profile found in OCDB!");
  // }

  // entry = AliCDBManager::Instance()->Get("GRP/Calib/CosmicTriggers");
  // if (entry) {
  //   fListOfCosmicTriggers = dynamic_cast<THashTable*>(entry->GetObject());
  //   entry->SetOwner(0);
  //   AliCDBManager::Instance()->UnloadFromCache("GRP/Calib/CosmicTriggers");
  // }

  // if (!fListOfCosmicTriggers) {
  //   AliWarning("Can not get list of cosmic triggers from OCDB! Cosmic event specie will be effectively disabled!");
  // }

  return kTRUE;
} 
