Char_t *partname[3] = {"pion", "kaon", "proton"};
Char_t *chargename[2] = {"positive", "negative"};
Char_t *chargename2[2] = {"plus", "minus"};

mcSpectra(Char_t *filename)
{

  TFile *fin = TFile::Open(filename);
  TH1D *hspectrumMB[3][2];
  TH1D *hspectrumCent[3][2][6];
  TH1 *h;
  THnSparseF *hsparse;
  TFile *fout = TFile::Open("mcSpectra.root", "RECREATE");

  for (Int_t ipart = 0; ipart < 3; ipart++)
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      hsparse = (THnSparseF *)fin->Get(Form("hHisto_hPrimaryTracks_%s_%s", partname[ipart], chargename[icharge]));
      h = hsparse->Projection(1);
      h->SetName(Form("mb_%s_%s", partname[ipart], chargename2[icharge]));
      h->Sumw2();
      for (Int_t ipt = 0; ipt < h->GetNbinsX(); ipt++) {
	h->SetBinContent(ipt + 1, h->GetBinContent(ipt + 1) / h->GetBinWidth(ipt + 1));
	h->SetBinError(ipt + 1, h->GetBinError(ipt + 1) / h->GetBinWidth(ipt + 1));
      }
      fout->cd();
      h->Write();
      for (Int_t icent = 0; icent < 6; icent++) {
	hsparse->GetAxis(0)->SetRange(icent + 1, icent + 1);
	h = hsparse->Projection(1);
	h->SetName(Form("cent%d_%s_%s", icent, partname[ipart], chargename2[icharge]));
	h->Sumw2();
	for (Int_t ipt = 0; ipt < h->GetNbinsX(); ipt++) {
	  h->SetBinContent(ipt + 1, h->GetBinContent(ipt + 1) / h->GetBinWidth(ipt + 1));
	  h->SetBinError(ipt + 1, h->GetBinError(ipt + 1) / h->GetBinWidth(ipt + 1));
	}
	fout->cd();
	h->Write();
      }
    }

}
