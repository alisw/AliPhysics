void 
SPDAcceptance(Int_t runNo=118560, Int_t year=2010)
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(Form("alien://Folder=/alice/data/%d/OCDB",year));
  cdb->SetSpecificStorage("ITS/Calib/CalibSSD", 
			  Form("alien://Folder=/alice/data/%d/OCDB",2012));
  cdb->SetRun(runNo);
  AliGeomManager::LoadGeometry();

  AliITSInitGeometry*  itsInit   = new AliITSInitGeometry();
  AliITSgeom*          itsGeom   = itsInit->CreateAliITSgeom();
  AliITSChannelStatus* itsStatus = new AliITSChannelStatus(cdb);

  Int_t spdFirst = itsGeom->GetStartSPD();
  Int_t spdLast  = 79; // itsGeom->GetLastSPD();

  Info("", "SPD range: %d - %d", spdFirst, spdLast);
  for (Int_t i = spdFirst; i <= spdLast; i++) {
    Double_t loc[] = { 0, 0, 0 };
    Double_t glb[] = { 0, 0, 0 };
    if (!AliITSgeomTGeo::LocalToGlobal(i, loc, glb)) continue;
    Double_t r = TMath::Sqrt(glb[0]*glb[0] + glb[1]*glb[1]);
    Double_t p = TMath::ATan2(glb[1], glb[0]) * TMath::RadToDeg();
    if (p < 0) p += 360;
    Printf("Module %3d: (x,y,z)=%10.6f,%10.6f,%10.6f, (r,phi)=%10.6f,%10.6f", 
	   i, glb[0], glb[1], glb[2], r, p);
  }
}
