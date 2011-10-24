void CheckTriggerOCDB(Int_t run = 164744)
{
  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");
  man->SetRun(run);

  AliCDBEntry *ent1 = man->Get("VZERO/Trigger/Data");
  AliVZEROTriggerData *trigData = (AliVZEROTriggerData*)ent1->GetObject();
  AliCDBEntry *ent2 = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calData = (AliVZEROCalibData*)ent2->GetObject();

  printf("%d <= MTA <= %d    %d <= MTC <= %d\n",
	 trigData->GetMultV0AThrLow(),
	 trigData->GetMultV0AThrHigh(),
	 trigData->GetMultV0CThrLow(),
	 trigData->GetMultV0CThrHigh());
  printf("CTA1 = %hd   CTC1 = %hd\n",
	 trigData->GetCentralityV0AThrLow(),
	 trigData->GetCentralityV0CThrLow());
  printf("CTA2 = %hd   CTC2 = %hd\n",
	 trigData->GetCentralityV0AThrHigh(),
	 trigData->GetCentralityV0CThrHigh());

  UShort_t fPed[64][2];
  UShort_t fPedCut[64][2];

  for(Int_t pmNumber = 0; pmNumber < 64; ++pmNumber) {
    Int_t board   = AliVZEROCalibData::GetBoardNumber(pmNumber);
    Int_t channel = AliVZEROCalibData::GetFEEChannelNumber(pmNumber);
			
    //    if(trigData->GetEnableCharge(board,channel)) {
    //    if(trigData->GetPedestalSubtraction(board)) {
    fPedCut[pmNumber][0] = trigData->GetPedestalCut(0,board,channel);
    fPedCut[pmNumber][1] = trigData->GetPedestalCut(1,board,channel);
    fPed[pmNumber][0] = trigData->GetPedestal(0,board,channel);
    fPed[pmNumber][1] = trigData->GetPedestal(1,board,channel);
    printf("Ch=%d Q=%s PedSub=%s Ped0=%.2f TrPed0=%hd TrPedCut0=%hd   Ped1=%.2f TrPed1=%hd TrPedCut1=%hd   delta(ped)=%.2f %f nsigma=%.2f %.2f\n",
	   pmNumber,
	   trigData->GetEnableCharge(board,channel) ? "On" : "Off",
	   trigData->GetPedestalSubtraction(board)  ? "On" : "Off",
	   calData->GetPedestal(pmNumber),fPed[pmNumber][0],fPedCut[pmNumber][0],
	   calData->GetPedestal(pmNumber+64),fPed[pmNumber][1],fPedCut[pmNumber][1],
	   (Float_t)fPed[pmNumber][0]-calData->GetPedestal(pmNumber),(Float_t)fPed[pmNumber][1]-calData->GetPedestal(pmNumber+64),
	   ((Float_t)fPedCut[pmNumber][0]-(Float_t)fPed[pmNumber][0])/calData->GetSigma(pmNumber),
	   ((Float_t)fPedCut[pmNumber][1]-(Float_t)fPed[pmNumber][1])/calData->GetSigma(pmNumber+64));
  }
   
}
