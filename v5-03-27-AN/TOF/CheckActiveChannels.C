CheckActiveChannels(const Char_t *);
CheckActiveChannels(Int_t nrun=-1);
CheckActiveChannelsFromCDBEntry(AliCDBEntry *);

//____________________________________________________________________

CheckActiveChannels(const Char_t *fileName)
{

  TFile *file = TFile::Open(fileName);
  AliCDBEntry *cdbe = (AliCDBEntry *)file->Get("AliCDBEntry");
  CheckActiveChannelsFromCDBEntry(cdbe);
  
}

//____________________________________________________________________

CheckActiveChannels(Int_t nrun)
{

  TGrid *alien = TGrid::Connect("alien://");
  if (!alien || !alien->IsConnected()) return;
  AliCDBManager *cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage("raw://");
  if (nrun==-1)
    cdbm->SetRun(80015);
  else
    cdbm->SetRun(nrun);
  AliCDBEntry *cdbe = cdbm->Get("TOF/Calib/Status");
  CheckActiveChannelsFromCDBEntry(cdbe);

}

//____________________________________________________________________

CheckActiveChannelsFromCDBEntry(AliCDBEntry *cdbe)
{

  AliTOFChannelOnlineStatusArray *array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();
  TH1F *hStatus_hw_ok = new TH1F("hStatus_hw_ok", "HW status;index;status ok", array->GetSize(), 0., array->GetSize(););
  TH1F *hStatus_pulser_ok = new TH1F("hStatus_pulser_ok", "pulser status;index;status ok", array->GetSize(), 0., array->GetSize(););
  TH1F *hStatus_noise_ok = new TH1F("hStatus_noise_ok", "noise status;index;status ok", array->GetSize(), 0., array->GetSize(););
  TH1F *hStatus_active = new TH1F("hStatus_active", "active status;index;status ok", array->GetSize(), 0., array->GetSize(););
  Bool_t hw_ok, pulser_ok, noise_ok;
  for (Int_t i = 0; i <  array->GetSize(); i++) {
    hw_ok = array->GetHWStatus(i) & AliTOFChannelOnlineStatusArray::kTOFHWOk;
    pulser_ok = !(array->GetPulserStatus(i) & AliTOFChannelOnlineStatusArray::kTOFPulserBad);
    noise_ok = array->GetNoiseStatus(i) & AliTOFChannelOnlineStatusArray::kTOFNoiseOk;
    hStatus_hw_ok->SetBinContent(i + 1, hw_ok);
    hStatus_pulser_ok->SetBinContent(i + 1, pulser_ok);
    hStatus_noise_ok->SetBinContent(i + 1, noise_ok);
    hStatus_active->SetBinContent(i + 1, hw_ok && pulser_ok && noise_ok);
  }
  printf("%d active channels\n", hStatus_active->Integral());

  TFile *fileout = TFile::Open("ChannelStatus.root", "RECREATE");
  hStatus_hw_ok->Write();
  hStatus_pulser_ok->Write();
  hStatus_noise_ok->Write();
  hStatus_active->Write();
  fileout->Close();

}


