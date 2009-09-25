{
 TFile *fin = TFile::Open("histosscaled.root");
 TCanvas *myc1 = new TCanvas("myc1","myc1",1);
// myc1->SetLogx(1);
// myc1->SetLogy(1);
 histosscaled.Print();
 histosscaled->FindObject("AnaElectron_hPtElectronScaled")->Draw("pe");
//AnaElectron_hPtElectronScaled->GetXaxis()->SetRangeUser(1,200);
// AnaElectron_hPtElectronScaled->GetYaxis()->SetRangeUser(1,4000);
// AnaElectron_hPtElectronScaled->SetTitle("pT of Electron");
// AnaElectron_hPtElectronScaled->SetStats(kFALSE);
// AnaElectron_hPtElectronScaled->Draw("pe");

// myc1->Print("electron.pdf");
}
