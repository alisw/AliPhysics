{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetPalette(1);

//open files and get lists
TList *l;
TFile f0("QAsym.proof.root");
l = (TList *)f0.Get("QAsymHists");


TH1F* fHpos[7];
TH1F* fHneg[7];


fHpos[0]=  (TH1F*) l->FindObject("fRecPtPos");
fHneg[0]=  (TH1F*) l->FindObject("fRecPtNeg");
fHpos[1]=  (TH1F*) l->FindObject("fRecPhiPos");
fHneg[1]=  (TH1F*) l->FindObject("fRecPhiNeg");
fHpos[2]=  (TH1F*) l->FindObject("fRecEtaPos");
fHneg[2]=  (TH1F*) l->FindObject("fRecEtaNeg");
fHpos[3]=  (TH1F*) l->FindObject("fRecEtaPtPos");
fHneg[3]=  (TH1F*) l->FindObject("fRecEtaPtNeg");

fHpos[4]=  (TH1F*) l->FindObject("fRecPtPosEta");
fHneg[4]=  (TH1F*) l->FindObject("fRecPtNegEta");
fHpos[5]=  (TH1F*) l->FindObject("fRecPhiPosEta");
fHneg[5]=  (TH1F*) l->FindObject("fRecPhiNegEta");
fHpos[6]=  (TH1F*) l->FindObject("fRecQPtPosEta");
fHneg[6]=  (TH1F*) l->FindObject("fRecQPtNegEta");


TLegend *legp1;
legp1= new TLegend(0.9,0.65,0.65,0.9);
legp1->SetFillColor(kWhite);

TLegend *legp2;
legp2= new TLegend(0.9,0.65,0.65,0.9);
legp2->SetFillColor(kWhite);

TF1 *fun0, *fun1;
fun0= new TF1("fun0","gaus",-5.0,5.0);
fun1= new TF1("fun1","gaus",-5.0,5.0);
fun0->SetLineColor(kBlack);
fun1->SetLineColor(kRed);

legp1->AddEntry(fun0,"pos. charge","l");
legp1->AddEntry(fun1,"neg. charge","l");

legp2->AddEntry(fun0,"pos. #eta","l");
legp2->AddEntry(fun1,"neg. #eta","l");



for(Int_t i=0; i<7;i++){
  fHpos[i]->SetLineColor(kBlack);
  fHpos[i]->SetLineWidth(2);
  fHneg[i]->SetLineColor(kRed);
  fHneg[i]->SetLineWidth(2);
}


TCanvas * dca = new TCanvas("pt", "pt", 100, 100, 1020, 820);
dca->Divide(3,3);
for(Int_t i=0; i<7;i++){
  dca->cd(i+1);
  fHpos[i]->Draw();
  fHneg[i]->Draw("same");
  if(i==0||i==4)gPad->SetLogy();
  if (i < 4) 
    legp1->Draw();
  else
    legp2->Draw();
}
  


}
