void CheckDelaysOCDB(Int_t run)
{
  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");
  man->SetRun(run);

  AliCDBEntry *ent = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calData = (AliVZEROCalibData*)ent->GetObject();

  AliCDBEntry *ent1 = man->Get("VZERO/Calib/TimeDelays");
  TH1F *delays = (TH1F*)ent1->GetObject();

  TH1F *hitdelays = new TH1F("hitdelays","HitDelay values",64,-0.5,63.5);
  hitdelays->SetLineColor(kRed);
  TH1F *sumdelays = new TH1F("sumdelays","Sum of the delays and HitDelay values",64,-0.5,63.5);
  sumdelays->SetLineColor(kGreen);

  for(Int_t i = 0; i < 64; ++i) {
    printf("Ch=%d delay=%.3f hitdelay=%.3f sum=%.3f\n",
	   i,
	   delays->GetBinContent(i+1),
	   calData->GetTimeOffset(i),
	   delays->GetBinContent(i+1)+calData->GetTimeOffset(i));
    hitdelays->SetBinContent(i+1,calData->GetTimeOffset(i));
    sumdelays->SetBinContent(i+1,delays->GetBinContent(i+1)+calData->GetTimeOffset(i));
  }

  new TCanvas;
  delays->GetYaxis()->SetRangeUser(-10,10);
  delays->Draw();
  hitdelays->Draw("same");
  sumdelays->Draw("same");
}
