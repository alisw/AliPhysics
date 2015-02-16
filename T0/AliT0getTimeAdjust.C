//------------------------------------------------------------------------
void AliT0getTimeAdjust(Int_t run)
{
  // Read calibration coefficients into the Calibration DB
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("T0/Calib/TimeAdjust");
 
  AliT0CalibSeasonTimeShift *clb = (AliT0CalibSeasonTimeShift*)entry->GetObject();
  clb->Print();
  Float_t *means= clb-> GetT0Means();
  Float_t *sigmas= clb-> GetT0Sigmas();
    for (Int_t i=0; i<4; i++) 
      cout<<means[i]<<" "<<sigmas[i]<<endl;

}
