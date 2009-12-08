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

fHpos=  (TH2F*) l->FindObject("fRecPtPosVz");
fHneg=  (TH2F*) l->FindObject("fRecPtNegVz");






TCanvas * dca = new TCanvas("pt", "pt", 100, 100, 1020, 520);
dca->Divide(2,1);
dca->cd(1);
fHpos->Draw("colz");
fHpos->SetTitle("positive tracks");
fHpos->SetYTitle("esdtrack->Zv()");

dca->cd(2);
fHneg->Draw("colz");
fHneg->SetMaximum(fHpos->GetMaximum());
fHneg->SetTitle("negative tracks");
fHneg->SetYTitle("esdtrack->Zv()");



}
  


}
