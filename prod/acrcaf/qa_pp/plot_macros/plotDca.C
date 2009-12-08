{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetPalette(1);

//open files and get lists
TList *l;
TFile f0("QAsym.proof.root");
l = (TList *)f0.Get("QAsymHists");


TH1F* fHpos;
TH1F* fHneg;
TH1F* fHnegInv;

fHpos=  (TH1F*) l->FindObject("fRecDcaPos");
fHpos->SetLineWidth(2);

fHneg=  (TH1F*) l->FindObject("fRecDcaNeg");
fHneg->SetLineColor(kRed);
fHneg->SetLineWidth(2);

fHnegInv=  (TH1F*) l->FindObject("fRecDcaNegInv");
fHnegInv->SetLineColor(kRed);
fHnegInv->SetLineWidth(2);

TLegend *legp;
legp= new TLegend(0.9,0.65,0.65,0.9);
legp->SetFillColor(kWhite);

TF1 *fun0, *fun1;
fun0= new TF1("fun0","gaus",-5.0,5.0);
fun1= new TF1("fun1","gaus",-5.0,5.0);
fun0->SetLineColor(kBlack);
fun1->SetLineColor(kRed);

legp->AddEntry(fun0,"positive","l");
legp->AddEntry(fun1,"negative","l");

TCanvas * dca = new TCanvas("pt", "pt", 100, 100, 1020, 480);
dca->Divide(2,1);
dca->cd(1);
fHpos->DrawClone("");
fHpos->SetTitle("");
fHpos->SetYTitle("");
fHneg->DrawClone("same");
gPad->SetLogy();
legp->Draw();
dca->cd(2);
fHpos->Divide(fHnegInv);
fHpos->Draw("");
fHpos->SetLineColor(kBlue);
fHpos->SetTitle("Yield_{positive}/Inverse Yield_{negative}");
fHpos->SetYTitle("");


}
  


}
