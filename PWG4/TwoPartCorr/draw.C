// sum.root made by doing
// cd output/
// hadd sum.root histos_*.root

void draw(const char* histFile = "output/sum.root")
{
  TFile* inFile = new TFile(histFile, "read");
  int nMultBinsRead = hMultBins->GetNbinsX();
  const int nMultBins = nMultBinsRead;
  int nPtBinsRead = hPtBins->GetNbinsX();
  const int nPtBins = nPtBinsRead;

  TH2F* hSig[nMultBins];
  TH2F* hBkg[nMultBins];
  TH2F* hSB[nMultBins];

  TH2F* hSigPt[nMultBins][nPtBins];
  TH2F* hBkgPt[nMultBins][nPtBins];
  TH2F* hSBPt[nMultBins][nPtBins];
  

  for (int iM=0; iM<nMultBins; iM++) {
    hSig[iM] = (TH2F*)inFile->Get(Form("hSig_%i", iM));
    hBkg[iM]  = (TH2F*)inFile->Get(Form("hBkg_%i", iM));
    hSB[iM] = (TH2F*)inFile->Get(Form("hSB_%i", iM));

    for (int iPt=0; iPt<nPtBins; iPt++) {
      hSigPt[iM][iPt] = (TH2F*)inFile->Get(Form("hSigPt_%i_%i",iM,iPt));
      hBkgPt[iM][iPt] = (TH2F*)inFile->Get(Form("hBkgPt_%i_%i",iM,iPt));
      hSBPt[iM][iPt]  = (TH2F*)inFile->Get(Form("hSBPt_%i_%i", iM,iPt));
    }

  }

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  hDEta->Draw();
  //  hEta->Draw("same");

  TCanvas* c2 = new TCanvas("c2", "c2", 1);
  c2->cd();
  hDPhi->Draw();

  TCanvas* c3 = new TCanvas("c3", "c3", 1);
  c3->cd();
  gPad->SetLogy();
  hPtAll->Draw("ep");
  hPtSel->SetLineColor(kRed);
  hPtSel->Draw("epsame");

  TCanvas* c4 = new TCanvas("c4", "c4", 1);
  c4->Divide(2, 1, 0.001, 0.001);
  c4->cd(1);
  hDR->Draw("surf2");
  c4->cd(2);
  hDRcut->Draw("colz");

  TCanvas* cEff = new TCanvas("cEff", "cEff", 1);
  cEff->cd();
  hMixEff->Draw("colz");

  TCanvas* cSig = new TCanvas("cSig", "cSig", 1);
  cSig->Divide(nMultBins, 1, 0.001, 0.001);
  for (int iM=0; iM<nMultBins; iM++) {
    cSig->cd(iM+1);
    hSig[iM]->Draw("surf1");
  }
  TCanvas* cBkg = new TCanvas("cBkg", "cBkg", 1);
  cBkg->Divide(nMultBins, 1, 0.001, 0.001);
  for (int iM=0; iM<nMultBins; iM++) {
    cBkg->cd(iM+1);
    hBkg[iM]->Draw("surf1");
  }
  
  TCanvas* cSB = new TCanvas("cSB", "cSB", 1);
  cSB->Divide(nMultBins, 1, 0.001, 0.001);
  for (int iM=0; iM<nMultBins; iM++) {
    cSB->cd(iM+1);
    hSB[iM]->Draw("surf1");
  }

  TCanvas* cSBpt[nPtBins];
  for (int iPt=0; iPt<nPtBins; iPt++) {
    cSBpt[iPt]= new TCanvas(Form("cSBpt_%i",iPt), Form("cSBpt_%i",iPt), 1);
    cSBpt[iPt]->Divide(nMultBins, 1, 0.001, 0.001);
    for (int iM=0; iM<nMultBins; iM++) {
      cSBpt[iPt]->cd(iM+1);
      hSBPt[iM][iPt]->Draw("surf1");
    }
  }
  
  TCanvas* cMult = new TCanvas("cMult", "cMult", 1);
  //  cMult->Divide(nMultBins, 1, 0.001, 0.001);
  cMult->cd();
  hMultAll->Draw();
  hMultSel->SetLineColor(kRed);
  hMultSel->Draw("same");

  return;
}
