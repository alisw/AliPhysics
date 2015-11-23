void CalculateThr(Float_t requiredThr = 4.5, Int_t run = 243398, Bool_t specificStorage = kFALSE)
{
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  if (specificStorage) {
    man->SetSpecificStorage("VZERO/Calib/Thresholds","local://$ALICE_ROOT/OCDB");
  }
  man->SetRun(run);

 AliCDBEntry *ent = man->Get("VZERO/Calib/Thresholds");
 TObjArray *arr = (TObjArray*)ent->GetObject();
 
 for(Int_t channel = 0; channel < 64; ++channel) {
   TF1 *fThr = (TF1*)arr->UncheckedAt(channel);
   Double_t newThr = fThr->GetX(requiredThr);
   printf("Ch=%d S%d R%d     Thr=%.1f\n",
	  channel,channel%8,(channel<32) ? channel/8 : (channel-32)/8,
	  newThr);
 }
}
