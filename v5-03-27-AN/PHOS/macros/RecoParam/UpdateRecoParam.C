UpdateRecoParam(const Int_t runNumber)
{
  // Read the array of PHOS recoparam objects from OCDB and update
  // EMC fitter version to "v4".
  // Write the updated object to OCDB.
  // Yuri Kharlov. 9.12.2011
  //
  /* $Id$ */

  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetRun(runNumber);

  AliCDBEntry* cdbEntry = AliCDBManager::Instance()->Get("PHOS/Calib/RecoParam");
  AliCDBMetaData *md = cdbEntry->GetMetaData();
  cout << "Responsible: " << md->GetResponsible() << endl;
  cout << "MD Comment : " << md->GetComment() << endl;
  TObjArray* arrayRecoParam = (TObjArray*)cdbEntry->GetObject();

  cout << "N recoparam = " << arrayRecoParam->GetEntries() << endl;

  AliPHOSRecoParam *rp = 0;
  for (Int_t i=0; i<arrayRecoParam->GetEntries(); i++) {
    rp = (AliPHOSRecoParam*)arrayRecoParam->At(i);
    printf("RP %d: event specie = %d, fitter version = %s\n",
	   i,rp->GetEventSpecie(),rp->EMCFitterVersion());
    rp->SetEMCFitterVersion("v4");
  }

  // Writing new recoparam to OCDB

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://OCDB");

  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("PHOS recoparameters: EMC fitter version is updated to v4");
  AliCDBId id("PHOS/Calib/RecoParam",167690,AliCDBRunRange::Infinity());
  cdb->Put(arrayRecoParam,id, md);

}
//-----------------------------------------------------------------------------

ReadRecoParam(const Int_t runNumber)
{
  // Read the array of PHOS recoparam objects from OCDB and
  // print its content to stdout

  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetRun(runNumber);

  AliCDBEntry* cdbEntry = AliCDBManager::Instance()->Get("PHOS/Calib/RecoParam");
  AliCDBMetaData *md = cdbEntry->GetMetaData();
  printf("Responsible: %s\n",md->GetResponsible());
  printf("MD Comment : %s\n",md->GetComment());
  TObjArray* arrayRecoParam = (TObjArray*)cdbEntry->GetObject();

  AliPHOSRecoParam *rp = 0;
  for (Int_t i=0; i<arrayRecoParam->GetEntries(); i++) {
    rp = (AliPHOSRecoParam*)arrayRecoParam->At(i);
    printf("Recoparam %d: event specie = %d\n",
	   i,rp->GetEventSpecie());

    printf("\tEMCClusteringThreshold = %g\n",rp->GetEMCClusteringThreshold());
    printf("\tEMCLocalMaxCut         = %g\n",rp->GetEMCLocalMaxCut());
    printf("\tEMCRawDigitThreshold   = %g\n",rp->GetEMCRawDigitThreshold());
    printf("\tEMCMinE                = %g\n",rp->GetEMCMinE());
    printf("\tEMCLogWeight           = %g\n",rp->GetEMCLogWeight());
    printf("\tEMCSampleQualityCut    = %g\n",rp->GetEMCSampleQualityCut());
    printf("\tEMCEcoreRadius         = %g\n",rp->GetEMCEcoreRadius());
    printf("\tEMCEcore2ESD           = %d\n",rp->EMCEcore2ESD());
    printf("\tEMCSubtractPedestals   = %d\n",rp->EMCSubtractPedestals());
    printf("\tEMCToUnfold            = %d\n",rp->EMCToUnfold());
    printf("\tEMCfitter version      = %s\n",rp->EMCFitterVersion());
    printf("\tEMCEnergyCorrectionOn  = %d\n",rp->GetEMCEnergyCorrectionOn());
    printf("\tGlobalAltroOffset      = %f\n",rp->GetGlobalAltroOffset());
    printf("\tGlobalAltroThreshold   = %d\n",rp->GetGlobalAltroThreshold());
    printf("\tTimeGateAmpThresh      = %g\n",rp->GetTimeGateAmpThresh());
    printf("\tTimeGateLow            = %g\n",rp->GetTimeGateLow());
    printf("\tTimeGateHigh           = %g\n",rp->GetTimeGateHigh());
    printf("\tNonlinearityCorrectionVersion = %s\n",rp->GetNonlinearityCorrectionVersion());
  }
}
