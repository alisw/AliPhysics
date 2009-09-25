{
 TFile *fin = TFile::Open("histos.root");
 TCanvas *myc1 = new TCanvas("myc1","myc1",1);
// myc1->SetLogx(1);
// myc1->SetLogy(1);
 histos.Print();
 histos->FindObject("AnaElectron_hPtElectron")->Draw("pe");
//AnaElectron_hPtElectron->GetXaxis()->SetRangeUser(1,200);
//AnaElectron_hPtElectron->GetYaxis()->SetRangeUser(1,4000);
//AnaElectron_hPtElectron->SetTitle("pT of Electron");
//AnaElectron_hPtElectron->SetStats(kFALSE);
//AnaElectron_hPtElectron->Draw("pe");

// myc1->Print("electron.pdf");
}
