AliESDVertex *
ReadMeanVertexFromOCDB(Int_t runNb, const Char_t* type = "MeanVertex")
{

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);
  AliCDBEntry *cdbe = cdb->Get(Form("GRP/Calib/%s", type));
  AliESDVertex *vertex = (AliESDVertex *)cdbe->GetObject();
  Double_t vmean[3], vsigma[3];
  vertex->GetXYZ(vmean);
  vertex->GetSigmaXYZ(vsigma);
  printf("vertex in run %d: z_mean = %f, z_sigma = %f\n", runNb, vmean[2], vsigma[2]);
  return vertex;
}
