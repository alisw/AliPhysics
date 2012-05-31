CheckTimeRawData(const Char_t *fileName, Int_t maxEv = kMaxInt)
{

  AliTOFRawStream::ApplyBCCorrections(kTRUE);

  gStyle->SetPalette(1);

  TH2F *hCrateTime = new TH2F("hCrateTime", ";crate;time (ns)", 72, 0, 72, 2000, 0, 2000);
  Float_t t;
  Int_t nPhysEv = 0;

  AliRawReaderRoot reader(fileName);
  AliTOFRawStream tofs(&reader);
  while (reader.NextEvent() && nPhysEv < maxEv) {
    if (reader.GetType() != 7) continue;
    if (nPhysEv % 100 == 0) printf("nPhysEv = %d\n", nPhysEv);
    nPhysEv++;
    for (Int_t i = 0; i < 72; i++) {
      tofs.LoadRawData(i);
      TClonesArray *array = tofs.GetRawData();
      for (Int_t j = 0; j < array->GetEntries(); j++) {
	AliTOFrawData *tofraw = (AliTOFrawData *)array->At(j);
	//	tofraw->Dump();
	t = tofraw->GetLeading() * 24.4e-3; /* ns */
	hCrateTime->Fill(i, t);
      }
    }
  }

  hCrateTime->Draw("colz");

}
