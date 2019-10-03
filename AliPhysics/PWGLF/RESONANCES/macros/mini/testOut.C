void testOut(const char *file = "test.root", Int_t p1 = 0, Int_t p2 = -1)
{
   TFile *f = TFile::Open(file);
   TList *l = (TList*)f->Get("RsnOut");
   
   l->Print();
   
   TH2F *hhPM = (TH2F*)RsnOut->FindObject("RSN_phi_Unlike");
   TH2F *hhPP = (TH2F*)RsnOut->FindObject("RSN_phi_LikePP");
   TH2F *hhMM = (TH2F*)RsnOut->FindObject("RSN_phi_LikeMM");
   TH2F *hhMX = (TH2F*)RsnOut->FindObject("RSN_phi_Mixing");
   
   TH1D *hPM = 0x0; if (hhPM) hPM = hhPM->ProjectionX(Form("px1_%d_%d", p1, p2), p1, p2);
   TH1D *hPP = 0x0; if (hhPP) hPP = hhPP->ProjectionX(Form("px2_%d_%d", p1, p2), p1, p2);
   TH1D *hMM = 0x0; if (hhMM) hMM = hhMM->ProjectionX(Form("px3_%d_%d", p1, p2), p1, p2);
   TH1D *hMX = 0x0; if (hhMX) hMX = hhMX->ProjectionX(Form("px4_%d_%d", p1, p2), p1, p2);
   
   hPM->SetLineColor(kBlack);
   
   if (hPP) {
      hPP->SetLineColor(kGreen);
      hPP->Add(hMM);
   }
   
   if (hMX) {
      hMX->SetLineColor(kRed);
      Double_t intS = hPM->Integral(hPM->GetXaxis()->FindBin(1.1), hPM->GetXaxis()->FindBin(1.2));
      Double_t intM = hMX->Integral(hMX->GetXaxis()->FindBin(1.1), hMX->GetXaxis()->FindBin(1.2));
      hMX->Scale(intS / intM);
   }
   
   TCanvas *c1 = new TCanvas("c1", "", 0, 0, 1100, 600);
   hPM->Draw();
   if (hPP) hPP->Draw("same");
   if (hMX) hMX->Draw("same");
   
   TH3F *hhRes = (TH3F*)RsnOut->FindObject("RSN_phi_Res");
   TH3F *hhNum = (TH3F*)RsnOut->FindObject("RSN_phi_Trues");
   TH3F *hhDen = (TH3F*)RsnOut->FindObject("RSN_phi_TrueMC");
   if (hhNum && hhDen) {
      TCanvas *c2 = new TCanvas("c2", "EFF", 0, 650, 800, 600);
      Double_t pt[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0};
      Int_t    npt  = sizeof(pt) / sizeof(pt[0]);
      
      TH1D *hNum  = hhNum->ProjectionY("num");
      TH1D *hDen  = hhDen->ProjectionY("den");
      TH1D *hrNum = hNum->Rebin(npt - 1, "neff", pt);
      TH1D *hrDen = hDen->Rebin(npt - 1, "deff", pt);
      hrNum->Divide(hrDen);
      hrNum->Draw();
      
      TCanvas *c3 = new TCanvas("c3", "", 500, 500, 600, 600);
      TH1D *hRes = hhRes->ProjectionX();
      hRes->Draw();
   }
}
