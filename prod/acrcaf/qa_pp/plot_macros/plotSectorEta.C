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
  fHpos[sector]=  (TH1F*) l->FindObject(Form("fRecEtaTpcSector%02d", sector));

}
//TH1F* fHposref =(TH1F*) l->FindObject("")

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
  // fHpos[i]->Divide(fHpos[0]);
  dca->cd(i+1);

  fHpos[i]->Draw();

  fHpos[i]->SetMaximum(5000);
  fHpos[i]->SetMinimum(0.1);
//   fHpos[i]->SetMaximum(3);
//   fHpos[i]->SetMinimum(-1);

  fHpos[i]->GetXaxis()->SetRangeUser(-1.5,1.5);
  fHpos[i]->SetTitle("");
  //fHpos[i]->SetXTitle(Form("#eta_{%d} / #eta_{0} ",i));
  fHpos[i]->SetXTitle(Form("#eta TPC sector %d ",i));

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.2);

  fHpos[i]->SetTitleOffset(1., "X");
  fHpos[i]->SetTitleOffset(1., "Y");
  fHpos[i]->SetTitleSize(0.1, "X");
  fHpos[i]->SetTitleSize(0.05, "Y");
  fHpos[i]->SetLabelSize(0.1, "X");
  fHpos[i]->SetLabelSize(0.1, "Y");
  //fHpos[i]->GetXaxis()->SetNdivision(5);
  fHpos[i]->GetXaxis()->SetNdivisions(5);
  fHpos[i]->GetYaxis()->SetNdivisions(7);

  gPad->SetLogy();

}
  


}
