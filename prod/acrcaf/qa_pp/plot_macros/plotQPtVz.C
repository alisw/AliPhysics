{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetPalette(1);

//open files and get lists
TList *l;
TFile f0("QAsym.proof.root");
l = (TList *)f0.Get("QAsymHists");


TH2F* fHpos;
TH2F* fHpneg;

fHpos=  (TH2F*) l->FindObject("fRecQPtPosEtaVz");
fHneg=  (TH2F*) l->FindObject("fRecQPtNegEtaVz");






TCanvas * dca = new TCanvas("pt", "pt", 100, 100, 1020, 520);
dca->Divide(2,1);
dca->cd(1);
//fHpos->Divide(fHneg);
gPad->SetRightMargin(0.15);
gPad->SetLeftMargin(0.15);
fHpos->SetTitleOffset(1.8,"Y");
//fHpos->SetXTitle("(Charge/p_{T})_{pos} / (Charge/p_{T})_{neg}");
fHpos->SetXTitle("Charge/p_{T}");
fHpos->Draw("colz");
fHpos->SetTitle("positive #eta");
fHpos->SetYTitle("esdtrack->Zv()");

dca->cd(2);
fHneg->Draw("colz");
gPad->SetRightMargin(0.15);
gPad->SetLeftMargin(0.15);
fHneg->SetTitleOffset(1.8,"Y");
fHneg->SetXTitle("Charge/p_{T}");
fHneg->SetMaximum(fHpos->GetMaximum());
fHneg->SetTitle("negative #eta");
fHneg->SetYTitle("esdtrack->Zv()");



}
  


}
