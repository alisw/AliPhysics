{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetPalette(1);

//open files and get lists
TList *l;
TFile f0("QAsym.proof.root");
l = (TList *)f0.Get("QAsymHists");


TH1F* fHpos[6];
TH1F* fHneg[6];

for(Int_t sector=0;sector<6;sector++){
 fHpos[sector]=  (TH1F*) l->FindObject(Form("fRecPhiPosLadder%d", sector));
 fHneg[sector]=  (TH1F*) l->FindObject(Form("fRecPhiNegLadder%d", sector));
} 

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



for(Int_t i=0; i<6;i++){
  fHpos[i]->SetLineColor(kBlack);
  fHpos[i]->SetLineWidth(2);
  fHneg[i]->SetLineColor(kRed);
  fHneg[i]->SetLineWidth(2);
}
      

TCanvas * dca = new TCanvas("pt", "pt", 100, 100, 920, 720);
dca->Divide(3,2);
for(Int_t i=0; i<6;i++){
  //  fHpos[i]->Divide(fHpos[0]);
  dca->cd(i+1);

  fHpos[i]->Draw();
  fHpos[i]->SetMaximum(10000);
  fHpos[i]->SetMinimum(0.1);
 

  fHpos[i]->SetTitle("");
  fHpos[i]->SetXTitle(Form("#phi of case %d",i));

//   gPad->SetLeftMargin(0.15);
//   gPad->SetBottomMargin(0.15);
  fHneg[i]->Draw("same");

//   fHpos[i]->SetTitleOffset(1., "X");
//   fHpos[i]->SetTitleOffset(1., "Y");
//   fHpos[i]->SetTitleSize(0.08, "X");
//   fHpos[i]->SetTitleSize(0.05, "Y");
//   fHpos[i]->SetLabelSize(0.1, "X");
//   fHpos[i]->SetLabelSize(0.1, "Y");
//   //fHpos[i]->GetXaxis()->SetNdivision(5);
//   fHpos[i]->GetXaxis()->SetNdivisions(5);


  gPad->SetLogy();

}
  


}
