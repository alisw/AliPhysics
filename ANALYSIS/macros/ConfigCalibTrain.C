//
// Macro to initialize: 
// - the OCDB (run number required as input argument)
// - the geometry (expected to be in the current directory)
// to run the Calibration train.
// 

void ConfigCalibTrain(Int_t run, const char *ocdb="raw://"){

  // OCDB

  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(run); 

  // geometry
  TGeoManager::Import("./geometry.root");
  AliGeomManager::LoadGeometry("./geometry.root");
  AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC");
}
