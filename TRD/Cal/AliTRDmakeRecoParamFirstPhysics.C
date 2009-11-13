//____________________________________________________
void AliTRDmakeRecoParamFirstPhysics()
{
  AliCDBMetaData *metaData= new AliCDBMetaData(); 
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexandru Bercuci");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-21-01"); //root version
  metaData->SetComment("First Physics reconstruction parameters for low, high and cosmic runs");
  
  AliCDBId id("TRD/Calib/RecoParam", 95352, AliCDBRunRange::Infinity()); 
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc) {
    return;
  }
  gStorLoc->Put(CreateRecoParamObject(), id, metaData); 

  return;
}


//____________________________________________________
TObjArray* CreateRecoParamObject()
{
  TObjArray *recos = new TObjArray(5);

  AliTRDrecoParam *rec = 0x0;
  recos->AddLast(rec = AliTRDrecoParam::GetLowFluxParam());
  rec->SetAsDefault();
  rec->SetNameTitle("Default", "TRD Default Reco Param");
  rec->SetRawStreamVersion("FAST");
  rec->SetXenon();
  rec->SetVertexConstrained();
  Double_t cov[3] = {2.,2.,0}
  rec->SetSysCovMatrix(cov);
  rec->SetChi2YSlope(0.11853);
  rec->SetChi2ZSlope(0.04527);
  rec->SetChi2YCut(1.);
  rec->SetPhiSlope(10.); //3.17954;
  rec->SetMaxTheta(2.1445);
  rec->SetMaxPhi(2.7475);
  rec->SetNMeanClusters(12.89);
  rec->SetNSigmaClusters(2.095);
  rec->SetRoadzMultiplicator = 3.;
  rec->SetStreamLevel(AliTRDrecoParam::kTracker, 1);

  recos->AddLast(rec = AliTRDrecoParam::GetLowFluxParam());
  rec->SetEventSpecie(AliRecoParam::kLowMult);
  rec->SetNameTitle("LOW", "TRD Low Flux Reco Param");
  rec->SetRawStreamVersion("FAST");
  rec->SetXenon();
  rec->SetVertexConstrained();
  rec->SetSysCovMatrix(cov);
  rec->SetChi2YSlope(0.11853);
  rec->SetChi2ZSlope(0.04527);
  rec->SetChi2YCut(1.);
  rec->SetPhiSlope(10.); //3.17954;
  rec->SetMaxTheta(2.1445);
  rec->SetMaxPhi(2.7475);
  rec->SetNMeanClusters(12.89);
  rec->SetNSigmaClusters(2.095);
  rec->SetRoadzMultiplicator(3.);
  rec->SetStreamLevel(AliTRDrecoParam::kTracker, 1);

  recos->AddLast(rec = AliTRDrecoParam::GetHighFluxParam());
  rec->SetEventSpecie(AliRecoParam::kHighMult);
  rec->SetNameTitle("HIGH", "TRD High Flux Reco Param");
  rec->SetRawStreamVersion("FAST");
  rec->SetXenon();
  rec->SetVertexConstrained();
  rec->SetSysCovMatrix(cov);
  rec->SetChi2YSlope(0.11853);
  rec->SetChi2ZSlope(0.04527);
  rec->SetChi2YCut(1.);
  rec->SetPhiSlope(10.); //3.17954;
  rec->SetMaxTheta(2.1445);
  rec->SetMaxPhi(2.7475);
  rec->SetNMeanClusters(12.89);
  rec->SetNSigmaClusters(2.095);
 
  recos->AddLast(rec = AliTRDrecoParam::GetCosmicTestParam());
  rec->SetEventSpecie(AliRecoParam::kCosmic);
  rec->SetNameTitle("COSMIC", "TRD Cosmic Reco Param");
  rec->SetRawStreamVersion("FAST");
  rec->SetXenon();

  recos->AddLast(rec = AliTRDrecoParam::GetCosmicTestParam());
  rec->SetEventSpecie(AliRecoParam::kCalib);
  rec->SetNameTitle("CALIBRATION", "TRD Calibration Reco Param");
  rec->SetRawStreamVersion("FAST");
  rec->SetXenon();

//  recos->AddLast(rec = AliTRDrecoParam::GetLowFluxParam());
//  rec->SetNameTitle("HLT", "TRD HLT Reco Param");
//  rec->SetChi2Y(.1);
//  rec->SetChi2Z(5.);

  return recos;
}
