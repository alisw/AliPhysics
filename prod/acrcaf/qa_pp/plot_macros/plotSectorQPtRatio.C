{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetPalette(1);

//open files and get lists
TList *l;
TFile f0("QAsym.proof.root");
l = (TList *)f0.Get("QAsymHists");


TH1F* fHpos[18];

for(Int_t sector=0;sector<18;sector++){
 fHpos[sector]=  (TH1F*) l->FindObject(Form("fRecQPtTpcSector%02d", sector));
} 

TH1F* fHposref =(TH1F*) l->FindObject("fQPt");

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



for(Int_t i=0; i<18;i++){
  fHpos[i]->SetLineColor(kRed);
  // fHpos[i]->SetLineWidth(2);
}
      

TCanvas * dca = new TCanvas("pt", "pt", 100, 100, 1020, 820);
dca->Divide(6,3);
for(Int_t i=0; i<18;i++){
  fHpos[i]->Divide(fHposref);
  dca->cd(i+1);

  fHpos[i]->Draw();
 
  fHpos[i]->SetMaximum(0.3);
  fHpos[i]->SetMinimum(-0.1);
//   fHpos[i]->SetMaximum(100);
//   fHpos[i]->SetMinimum(1)
    ;
 
  fHpos[i]->SetTitle("");
  fHpos[i]->SetXTitle(Form("(Q/p_{T})_{%d}/(Q/p_{T})_{all}",i));
  //fHpos[i]->SetXTitle(Form("Q/p_{T} sector %2d",i));

  gPad->SetLeftMargin(0.20);
  gPad->SetBottomMargin(0.17);

  fHpos[i]->SetTitleOffset(1., "X");
  fHpos[i]->SetTitleOffset(.8, "Y");
  fHpos[i]->SetTitleSize(0.08, "X");
  fHpos[i]->SetTitleSize(0.05, "Y");
  fHpos[i]->SetLabelSize(0.1, "X");
  fHpos[i]->SetLabelSize(0.1, "Y");
  fHpos[i]->GetYaxis()->SetNdivisions(5);
  fHpos[i]->GetXaxis()->SetNdivisions(5);


  //gPad->SetLogy();

}
  


}
