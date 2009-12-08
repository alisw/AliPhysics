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

fHpos=  (TH1F*) l->FindObject("fRecDcaPos");
fHpos->SetLineWidth(2);

fHneg=  (TH1F*) l->FindObject("fRecDcaNeg");
fHneg->SetLineColor(kRed);
fHneg->SetLineWidth(2);

TH1F* fHposD;
TH1F* fHnegD;

fHposD=  (TH1F*) l->FindObject("fRecDPos");
fHposD->SetLineWidth(2);

fHnegD=  (TH1F*) l->FindObject("fRecDNeg");
fHnegD->SetLineColor(kRed);
fHnegD->SetLineWidth(2);


fH=  (TH2F*) l->FindObject("fDiffDcaD");


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
dca->Divide(2,2);
dca->cd(1);
fHpos->DrawClone("");
fHpos->SetTitle("");
fHpos->SetYTitle("");
fHneg->DrawClone("same");
gPad->SetLogy();
legp->Draw();

dca->cd(2);
fHposD->DrawClone();
fHnegD->DrawClone("same");
legp->Draw();
gPad->SetLogy();


dca->cd(3);
fH->Draw("colz");

}
  


}
