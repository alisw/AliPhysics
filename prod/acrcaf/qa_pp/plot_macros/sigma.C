//void sigma( TString fname)
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  //open files and get lists
  TList *l;
 
  
  //TFile f0(fname);
  TFile f0("QAsym.proof.root");
  l = (TList *)f0.Get("QAsymHists");
  
  
  TH1F* fHpos[7];
  TH1F* fHneg[7];
  
  for(int i=0;i<7;i++){
    fHpos[i]=  (TH1F*) l->FindObject(Form("fPtSigmaPos%d", i));
    fHpos[i]->SetLineColor(i+1);
    fHpos[i]->SetLineWidth(2);
    fHpos[i]->GetXaxis()->SetLabelSize(0.06);
    fHpos[i]->GetYaxis()->SetLabelSize(0.06);
    fHpos[i]->GetXaxis()->SetTitleSize(0.06);
    fHpos[i]->GetYaxis()->SetTitleSize(0.06);
  
    fHneg[i]=  (TH1F*) l->FindObject(Form("fPtSigmaNeg%d", i));
    fHneg[i]->SetLineColor(i+1);
    fHneg[i]->SetLineWidth(2);
    fHneg[i]->GetXaxis()->SetLabelSize(0.06);
    fHneg[i]->GetYaxis()->SetLabelSize(0.06);
    fHneg[i]->GetXaxis()->SetTitleSize(0.06);
    fHneg[i]->GetYaxis()->SetTitleSize(0.06);
  }
      
  TLegend *legp;
  legp= new TLegend(0.9,0.45,0.83,0.9);
  legp->SetFillColor(kWhite);
  TF1 *fun[7];
  for(int i=0;i<7;i++){
    fun[i]= new TF1(Form("fun%d", i),"gaus",-5.0,5.0);
    fun[i]->SetLineColor(i+1);
    legp->AddEntry(fun[i],Form("case %d",i),"l");  
  }

 





  TCanvas * dca = new TCanvas("sigma", "sigma", 100, 100, 1020, 720);
  dca->Divide(1,2);

  for(int i=0;i<7;i++){

    dca->cd(1);
    gPad->SetBottomMargin(0.18);
    // fHpos[i]->GetXaxis()->SetRangeUser(-3.,3.);
    fHpos[i]->SetTitle("positive particles");
    fHpos[i]->SetMaximum(1000);
    fHpos[i]->SetMinimum(0.1);
    if(i==0)fHpos[i]->Draw();
    else fHpos[i]->Draw("same");
    gPad->SetLogy();
    //gPad->SetLogx();
    legp->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    //fHpos[i]->GetXaxis()->SetRangeUser(-4,0);

    dca->cd(2);
    gPad->SetBottomMargin(0.18);
    // fHneg[i]->GetXaxis()->SetRangeUser(-3.,3.);
    fHneg[i]->SetTitle("negative particles");
    fHneg[i]->SetMaximum(1000);
    fHneg[i]->SetMinimum(0.1);
    if(i==0)fHneg[i]->Draw();
    else fHneg[i]->Draw("same");

    gPad->SetLogy();
    //gPad->SetLogx();
    legp->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    //fHneg[i]->GetXaxis()->SetRangeUser(-4,0);
 


  }
  


}
