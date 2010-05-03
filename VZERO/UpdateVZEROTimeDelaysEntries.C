void UpdateVZEROTimeDelaysEntries(Int_t year = 2010, Int_t period = 0, Bool_t local = kFALSE)
{
  Int_t runRange[2];
  if (year == 2009) {
    runRange[0] = 0;
    runRange[1] = 105268;
  }
  else if (year == 2010) {
    if (period == 0) {
      runRange[0] = 106031;
      runRange[1] = 116353;
    }
    else if (period == 1) {
      runRange[0] = 116354;
      runRange[1] = 118285;
    }
    else if (period == 2) {
      runRange[0] = 118486;
      runRange[1] = 118556;
    }
    else {
      printf("Invalid run period...\n");
      return;
    }
  }
  else {
    printf("Invalid year...\n");
    return;
  }

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(Form("alien://?folder=/alice/data/%d/OCDB",year));

  man->SetRun(runRange[1]);

  AliCDBEntry *entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibda = (AliVZEROCalibData*)entry->GetObject();
  printf("Year %d Period %d:\n",year,period);
  printf("Calib:\n");
  for(Int_t i = 0; i < 64; ++i) printf("%.2f ",calibda->GetTimeOffset(i));
  printf("\n");
  entry = man->Get("VZERO/Calib/TimeDelays");
  TH1F *delays = (TH1F*)entry->GetObject();
  TH1F *delaysNew = new TH1F(*delays);
  printf("Delay:\n");
  for(Int_t i = 0; i < 64; ++i) printf("%.2f ",delaysNew->GetBinContent(i+1));
  printf("\n");
  for(Int_t i = 0; i < 64; ++i) {
    if (year == 2009) delaysNew->SetBinContent(i+1,delays->GetBinContent(i+1)+5.0);
    if (year == 2010 && period == 0) delaysNew->SetBinContent(i+1,5.0);
    if (year == 2010 && (period == 1 || period == 2)) {
      Int_t board = i / 8;
      Int_t channel = i % 8;
      Int_t j = AliVZEROCalibData::GetOfflineChannelNumber(board,channel);
      delaysNew->SetBinContent(j+1,
			       delays->GetBinContent(j+1)+
			       calibda->GetTimeOffset(i)-
			       calibda->GetTimeOffset(j));
    }
  }
  printf("CorrDelay:\n");
  for(Int_t i = 0; i < 64; ++i) printf("%.2f ",delaysNew->GetBinContent(i+1));
  printf("\n");

  if (local) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  {
    AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
    md->SetResponsible("Brigitte Cheynis");
    md->SetBeamPeriod(0);
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    md->SetComment(Form("Time delays channel by channel (corrected values for year %d, period %d)",year,period));
    md->PrintMetaData();

    AliCDBStorage *storLoc = man->GetDefaultStorage();
    AliCDBId id("VZERO/Calib/TimeDelays",runRange[0],runRange[1]);

    storLoc->Put(delaysNew, id, md);

    storLoc->Delete();
    delete md;
  }

}
