void CreateRecPars_CDB(){
  // Create TOF Calibration Object for Ideal calibration and 
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliTOFcalib *tofcalib = new AliTOFcalib();

  TObjArray * array = new TObjArray(2);
  array->Clear();
  AliTOFRecoParam *param = new AliTOFRecoParam();
  param->SetTimeResolution(100.);
  AliTOFRecoParam *paramPbPb = param->GetPbPbparam();
  paramPbPb->SetTimeResolution(100.);
  AliTOFRecoParam *parampp = param->GetPPparam();
  parampp->SetTimeResolution(100.);

  array->AddLast(parampp);
  array->AddLast(paramPbPb);

  tofcalib->WriteRecParOnCDB("TOF/Calib",0,999999999,array);
}
