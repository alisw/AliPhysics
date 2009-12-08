//void posNegInv(TString fname)
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  //open files and get lists
  TList *l0;
  // TFile f0(fname);
  TFile f0("QAsym.proof.root");
  l0 = (TList *)f0.Get("QAsymHists");

  // TList *l1;
  // TFile f1("QAsym.proof.fp.all.ITSrefit.root");
  // l1 = (TList *)f1.Get("QAsymHists");

  TH1F* fHpos[6];
  TH1F* fHneg[6];

  for(int i=0;i<6;i++){
    fHpos[i]=  (TH1F*) l0->FindObject(Form("fSignDcaPos%d", i));
    fHpos[i]->SetLineColor(kBlack);
    fHpos[i]->SetLineWidth(2);
    fHneg[i]=  (TH1F*) l0->FindObject(Form("fSignDcaNegInv%d", i));
    fHneg[i]->SetLineColor(kRed);
    fHneg[i]->SetLineWidth(2);
  }
 


     
  TLegend *legp[6];
  for(int i=0;i<6;i++){
    legp[i]= new TLegend(0.9,0.65,0.55,0.9);
    legp[i]->SetFillColor(kWhite);
    TF1 *fun[2];
    fun[0]= new TF1(Form("fun%d", i),"gaus",-5.0,5.0);
    fun[0]->SetLineColor(kBlack);
    fun[1]= new TF1(Form("fun%d", i),"gaus",-5.0,5.0);
    fun[1]->SetLineColor(kRed);
    legp[i]->AddEntry(fun[0],"Positive Particles  f( x)","l");  
    legp[i]->AddEntry(fun[1],"Negative partivcles f(-x)","l");  
  }

 



  Double_t range= 4.;

  TCanvas * dca = new TCanvas("dca", "dca", 100, 100, 880, 680);
  dca->Divide(2,3);
  for(Int_t i=0;i<6;i++){

    dca->cd(i+1);
    fHpos[i]->SetTitle(Form("Case %d",i));
    fHpos[i]->SetMaximum(1000000);
    fHpos[i]->SetMinimum(0.1);
    fHpos[i]->Draw();
    fHneg[i]->Draw("same");
    gPad->SetLogy();
    legp[i]->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // fHpos[i]->GetXaxis()->SetRangeUser(-range,range);
 
  }


}
  


