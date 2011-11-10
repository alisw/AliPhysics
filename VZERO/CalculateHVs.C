void CalculateHVs(Float_t requiredADCperMIP = 2.6, Int_t run = 137366, Bool_t specificStorage = kTRUE)
{
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  if (specificStorage) {
    man->SetSpecificStorage("VZERO/Calib/PMGains","local://$ALICE_ROOT/OCDB");
    man->SetSpecificStorage("VZERO/Calib/LightYields","local://$ALICE_ROOT/OCDB");
  }
  man->SetRun(run);

  AliCDBEntry *ent = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calData = (AliVZEROCalibData*)ent->GetObject();
  for(Int_t i = 0; i < 64; ++i) {
    printf("%d %.0f (%.0f) S%d R%d     Delta=%.0f\n",
	   i,
	   calData->GetHV(i,requiredADCperMIP),
	   calData->GetMeanHV(i),
	   i%8,
	   (i<32) ? i/8 : (i-32)/8,
	   calData->GetHV(i,requiredADCperMIP)-calData->GetMeanHV(i));
  }
}
