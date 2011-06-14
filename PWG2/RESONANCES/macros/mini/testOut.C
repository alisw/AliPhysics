void testOut(const char *file = "test.root", Double_t scale = 0.1, Int_t p1 = 0, Int_t p2 = -1)
{
   TFile *f = TFile::Open(file);
   TList *l = (TList*)f->Get("RsnOut");
   
   l->Print();
   
   TH2F *hhPM = (TH2F*)RsnOut->FindObject("hist_RSN_phi_Unlike_default");
   TH2F *hhPP = (TH2F*)RsnOut->FindObject("hist_RSN_phi_LikePP_default");
   TH2F *hhMM = (TH2F*)RsnOut->FindObject("hist_RSN_phi_LikeMM_default");
   TH2F *hhMX = (TH2F*)RsnOut->FindObject("hist_RSN_phi_Mixing_default");
   
   TH1D *hPM = 0x0; if (hhPM) hPM = hhPM->ProjectionX(Form("px1_%d_%d", p1, p2), p1, p2);
   TH1D *hPP = 0x0; if (hhPP) hPP = hhPP->ProjectionX(Form("px2_%d_%d", p1, p2), p1, p2);
   TH1D *hMM = 0x0; if (hhMM) hMM = hhMM->ProjectionX(Form("px3_%d_%d", p1, p2), p1, p2);
   TH1D *hMX = 0x0; if (hhMX) hMX = hhMX->ProjectionX(Form("px4_%d_%d", p1, p2), p1, p2);
   
   hPM->SetLineColor(kBlack);
   
   hPP->SetLineColor(kGreen);
   hPP->Add(hMM);
   
   hMX->SetLineColor(kRed);
   hMX->Scale(scale);
   
   TCanvas *c1 = new TCanvas("c1");
   hPM->Draw();
   hPP->Draw("same");
   hMX->Draw("same");
}
