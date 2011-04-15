//
// Macro to initialize: 
// - the OCDB (run number required as input argument)
// - the geometry (expected to be in the current directory)
// to run the Calibration train.
// 

void ConfigCalibTrain(Int_t run, const char *ocdb="raw://"){

  // OCDB
 
  printf("setting run to %d\n",run);
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(run); 

  // geometry
  AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC");



  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = (AliGRPObject*)entry->GetObject();

  Bool_t ok=kTRUE;
  Float_t l3Current = grpData->GetL3Current((AliGRPObject::Stats)0);
  if (l3Current == AliGRPObject::GetInvalidFloat()) {
    printf("GRP/GRP/Data entry:  missing value for the L3 current !");
    ok = kFALSE;
  }
  
  Char_t l3Polarity = grpData->GetL3Polarity();
  if (l3Polarity == AliGRPObject::GetInvalidChar()) {
    printf("GRP/GRP/Data entry:  missing value for the L3 polarity !");
    ok = kFALSE;
  }
  
  // Dipole
  Float_t diCurrent = grpData->GetDipoleCurrent((AliGRPObject::Stats)0);
  if (diCurrent == AliGRPObject::GetInvalidFloat()) {
    printf("GRP/GRP/Data entry:  missing value for the dipole current !");
    ok = kFALSE;
  }
  
  Char_t diPolarity = grpData->GetDipolePolarity();
  if (diPolarity == AliGRPObject::GetInvalidChar()) {
    printf("GRP/GRP/Data entry:  missing value for the dipole polarity !");
    ok = kFALSE;
  }

  TString beamType = grpData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    printf("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = grpData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    printf("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  // read special bits for the polarity convention and map type
  //Int_t  polConvention = grpData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
  Int_t  polConvention = grpData->IsPolarityConventionLHC() ? 0 : 1;
  Bool_t uniformB = grpData->IsUniformBMap();
  
  if (ok) {
    AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1),
					   TMath::Abs(diCurrent) * (diPolarity ? -1:1),
					   polConvention,uniformB,beamEnergy, beamType.Data());
    if (fld) {
      TGeoGlobalMagField::Instance()->SetField( fld );
      TGeoGlobalMagField::Instance()->Lock();
      printf("Running with the B field constructed out of GRP !");
    }
  }
  printf("Problem with magnetic field setup\n");
}
