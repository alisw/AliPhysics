void TestMatrices(Bool_t test2011=1)
{
  AliLog::SetClassDebugLevel("AliCDBManager",1);

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  const char *geoname = "EMCAL_COMPLETEV1";
  if (test2011) {
    man->SetSpecificStorage("EMCAL/Align/Data",             "alien://folder=/alice/data/2011/OCDB");
    man->SetRun(146805);
  } else {
    man->SetSpecificStorage("EMCAL/Align/Data",             "alien://folder=/alice/data/2010/OCDB");
    man->SetRun(137366);
    geoname = "EMCAL_FIRSTYEARV1";
  }
  AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("EMCAL");
  AliEMCALGeometry *geo =  AliEMCALGeometry::GetInstance(geoname);
  for (Int_t i=0;i<(geo->GetEMCGeometry())->GetNumberOfSuperModules();++i)
    geo->GetMatrixForSuperModule(i)->Print();

  AliEMCALEMCGeometry *emc = geo->GetEMCGeometry();
  Double_t phimin = emc->GetArm1PhiMin();
  Double_t phimax = emc->GetArm1PhiMax();
  cout << phimin << " " << phimax << endl;
  emc->PrintGeometry();
}
