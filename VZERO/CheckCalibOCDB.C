void CheckCalibOCDB(Int_t run)
{
  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");
  man->SetRun(run);

  AliCDBEntry *ent = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calData = (AliVZEROCalibData*)ent->GetObject();

  for(Int_t pmNumber = 0; pmNumber < 64; ++pmNumber) {
    printf("Ch=%d HV=%.1f MIP=%.3f Ped0=%.1f Ped1=%.1f Thr=%.1f CalThr=%.2f Dead=%s\n",
	   pmNumber,
	   calData->GetMeanHV(pmNumber),
	   1./calData->GetMIPperADC(pmNumber),
	   calData->GetPedestal(pmNumber),calData->GetPedestal(pmNumber+64),
	   calData->GetDiscriThr(pmNumber),calData->GetCalibDiscriThr(pmNumber,kTRUE,run),
	   calData->IsChannelDead(pmNumber)?"yes":"no");
  }
}
