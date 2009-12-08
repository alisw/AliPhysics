{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetPalette(1);

//open files and get lists
TList *l;
TFile f0("QAsym.proof.root");
l = (TList *)f0.Get("QAsymHists");


TH1F* fHpos;
TH1F* fHpneg;

fHpos=  (TH1F*) l->FindObject("fRecEtaPosVz");
fHneg=  (TH1F*) l->FindObject("fRecEtaNegVz");






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
