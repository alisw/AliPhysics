void rec(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2009 RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////

  
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  
  // Reconstruction settings
  AliReconstruction rec;

  // Set protection against too many events in a chunk (should not happen)
  //  rec.SetEventRange(0,30000);

  // Set reconstruction flags (skip detectors here if neded
  rec.SetRunReconstruction("ITS TPC TRD TOF");
  //rec.SetFillESD("ITS TPC TRD");
  // QA options
  //rec.SetRunQA("Global:ESDs") ;
  
  //rec.SetRunQA(":") ;
  AliQAManager *qam = AliQAManager::QAManager(AliQAv1::kRECMODE) ;
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);
  //ITS QA Off (https://savannah.cern.ch/bugs/?60187)   
  //rec.SetRunQA("ALL -HLT:ALL") ;

  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  rec.SetUseTrackingErrorsForAlignment("ITS");

  // Magnetic field hack for L3 off, dipole on

  AliMagF* fld = new AliMagF("map","map",0,-1, AliMagF::k5kG,AliMagF::kBeamTypepp, 450);
  fld->SetBit(AliMagF::kOverrideGRP);
  TGeoGlobalMagField::Instance()->SetField(fld);
  TGeoGlobalMagField::Instance()->Lock();
  printf(" ATTENTION: Using external field with WRONG CURRENTS COMBINATION\n");

  // Specially for ITS (https://savannah.cern.ch/bugs/?59368)

  rec.SetRunPlaneEff(kTRUE); 

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

  AliLog::Flush();
  rec.Run();

}
