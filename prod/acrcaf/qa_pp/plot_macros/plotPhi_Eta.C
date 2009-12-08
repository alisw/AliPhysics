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

fHpos=  (TH1F*) l->FindObject("fRecPhiPosEta");
fHpos->SetLineWidth(2);

fHneg=  (TH1F*) l->FindObject("fRecPhiNegEta");
fHneg->SetLineColor(kRed);
fHneg->SetLineWidth(2);

TLegend *legp;
legp= new TLegend(0.9,0.65,0.65,0.9);
legp->SetFillColor(kWhite);

TF1 *fun0, *fun1;
fun0= new TF1("fun0","gaus",-5.0,5.0);
fun1= new TF1("fun1","gaus",-5.0,5.0);
fun0->SetLineColor(kBlack);
fun1->SetLineColor(kRed);

legp->AddEntry(fun0,"positive #eta","l");
legp->AddEntry(fun1,"negative #eta","l");

TCanvas * dca = new TCanvas("pt", "pt", 100, 100, 1020, 480);
dca->Divide(2,1);
dca->cd(1);
fHpos->DrawClone("");
fHpos->SetTitle("");
fHpos->SetYTitle("");
fHneg->DrawClone("same");
//gPad->SetLogy();
legp->Draw();
dca->cd(2);
fHpos->Divide(fHneg);
fHpos->Draw("");
fHpos->SetLineColor(kBlue);
fHpos->SetTitle("Yield_{positive #eta}/Yield_{negative #eta}");
fHpos->SetYTitle("");


}
  


}
