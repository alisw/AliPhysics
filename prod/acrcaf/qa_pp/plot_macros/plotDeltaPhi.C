{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetPalette(1);

//open files and get lists
TList *l;
TFile f0("QAsym.proof.root");
l = (TList *)f0.Get("QAsymHists");


TH2F* fHpos;
fHpos=  (TH2F*) l->FindObject("fDeltaPhiLeading");

TH1D* px;
TH1D* py;

//pgood= fHpos->ProjectionY();
//pgood->Draw();

TCanvas * c1 = new TCanvas("c1", "c1", 100, 100, 1020, 480);
c1->Divide(2,1);
// c1->cd(1);
// fHpos->Draw("colz");


c1->cd(1);
px= fHpos->ProjectionX();
px->Draw();
px->SetTitle("");

c1->cd(2);
py= fHpos->ProjectionY();
py->Draw();
py->SetMinimum(0);
py->SetTitle("");

TCanvas * c2 = new TCanvas("c2", "c2", 100, 100, 620, 480);
c2->cd();
gPad->SetLeftMargin(0.13);
fHpos->SetTitleOffset(1.2, "Y");
fHpos->Draw("colz");
fHpos->SetTitle("");



}
