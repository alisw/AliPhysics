void makefig(const int trig=0, const int icut=0, const int icent=0){

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetFillColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderSize(0);
  //gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadBottomMargin(0.125);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetTitleYOffset(1.3);
  //gStyle->SetPadLeftMargin(0.1);
  gStyle->SetTitleW(0.7);
  gStyle->SetTitleH(0.1);
  cout<<"physics is always fun! "<<endl; 

  TH2D *hmasspt[7][11];

  ////// input root file ///////////////
  char name[100];
  sprintf(name,"result_trig%d_cut%d.root", trig, icut);
  TFile *fin = new TFile(name,"read");

  for(int i=0;i<7;i++){
    sprintf(name,"hmasspt_mb_%d",i);
    hmasspt[i][10] = new TH2D(name, name, 500, 0, 5, 500, 0, 5);
  }    

  ////// add into the histograms //////////////////
  for(int i=0;i<7;i++){
    for(int j=0;j<10;j++){
      sprintf(name,"hmasspt_cent%d_pair%d", j, i);
      hmasspt[i][j] = (TH2D*)fin->Get(name);
      hmasspt[i][10]->Add(hmasspt[i][j]);
    }
  }


  TH1D *hpy = hmasspt[0][0]->ProjectionY("hpy");
  
  ////// range of pT /////////////////////////////////
  int npt=3;
  float low[4] ={0, 1, 2};
  float high[4]={1, 2, 4};
  
  //// [pairtype][pt] //////////////////////////////
  TH1D *hmass[7][4]; ///sliced mass spectrum
  for(int i=0;i<7;i++){
    for(int ipt=0;ipt<npt;ipt++){
      int rangelow = hpy->FindBin(low[ipt]);  
      int rangehigh = hpy->FindBin(high[ipt]);
      sprintf(name,"hmass_tmp_pt%d_cent%d_%d",ipt, icent, i);
      hmass[i][ipt] = (TH1D*)hmasspt[i][icent]->ProjectionX(name, rangelow, rangehigh-1);
    }
  }

  //////// another historgram to merge like sign/mixed like sign
  //////// hmass_add will be scaled.
  ///////  hmass_org will be unchanged.
  TH1D *hmass_add[4][4];
  TH1D *hmass_org[4][4];
  for(int i=0;i<4;i++){
    for(int ipt=0;ipt<npt;ipt++){
      sprintf(name,"hmass_pt%d_cent%d_%d",ipt, j, i);
      hmass_add[i][ipt] = new TH1D(name, name, 500, 0, 5);
      hmass_add[i][ipt]->Sumw2();
      sprintf(name,"hmass_pt%d_cent%d_%d_org",ipt, j, i);
      hmass_org[i][ipt] = new TH1D(name, name, 500, 0, 5);
      hmass_org[i][ipt]->Sumw2();
    }
  }
  for(int ipt=0;ipt<npt;ipt++){
    hmass_add[0][ipt]->Add(hmass[0][ipt]);
    hmass_add[1][ipt]->Add(hmass[1][ipt]);
    hmass_add[1][ipt]->Add(hmass[2][ipt]);
    hmass_add[2][ipt]->Add(hmass[3][ipt]);
    hmass_add[2][ipt]->Add(hmass[4][ipt]);
    hmass_add[3][ipt]->Add(hmass[5][ipt]);
    hmass_add[3][ipt]->Add(hmass[6][ipt]);

    hmass_org[0][ipt]->Add(hmass[0][ipt]);
    hmass_org[1][ipt]->Add(hmass[1][ipt]);
    hmass_org[1][ipt]->Add(hmass[2][ipt]);
    hmass_org[2][ipt]->Add(hmass[3][ipt]);
    hmass_org[2][ipt]->Add(hmass[4][ipt]);
    hmass_org[3][ipt]->Add(hmass[5][ipt]);
    hmass_org[3][ipt]->Add(hmass[6][ipt]);


    ///// scaling factor for like sign and mixed like sign
    float scale1 = hmass_add[0][ipt]->GetEntries()/(2*sqrt(hmass[1][ipt]->GetEntries()*hmass[2][ipt]->GetEntries()));
    float scale2 = hmass_add[2][ipt]->GetEntries()/(2*sqrt(hmass[5][ipt]->GetEntries()*hmass[6][ipt]->GetEntries()));
    
    hmass_add[1][ipt]->Scale(scale1);
    hmass_add[3][ipt]->Scale(scale2);
    hmass_org[1][ipt]->Scale(scale1);
    hmass_org[3][ipt]->Scale(scale2);
  }


  ///// then take ratio of mixed unlike / mixed like
  TH1F *hdiv[4]; 
  for(int ipt=0;ipt<npt;ipt++){
    sprintf(name,"hdiv_pt%d",ipt);
    hdiv[ipt] = new TH1F(name, name, 500, 0, 5);
    hdiv[ipt]->Sumw2();
    hdiv[ipt]->Divide(hmass_add[2][ipt], hmass_add[3][ipt]);
    hdiv[ipt]->SetMaximum(2);
    hdiv[ipt]->SetMinimum(0);

    //// multiply to like sign 
    hmass_add[1][ipt]->Multiply(hdiv[ipt]);

  }

  TCanvas *c0 = new TCanvas("c0","c0",1200,800);
  c0->Divide(npt,2);
  for(int ipt=0;ipt<npt;ipt++){
    c0->cd(ipt+1);
    gPad->SetLogy();
    gPad->SetGrid(1);
    hmass_add[0][ipt]->SetXTitle("M_{ee} [GeV]");
    hmass_add[0][ipt]->SetYTitle("counts");
    hmass_add[0][ipt]->SetAxisRange(0,4);
    hmass_add[0][ipt]->SetLineColor(1);
    hmass_org[1][ipt]->SetLineColor(2);
    hmass_add[0][ipt]->SetMinimum(1);
    hmass_add[0][ipt]->Draw();
    hmass_org[1][ipt]->Draw("same");

    c0->cd(ipt+1+npt);
    gPad->SetLogy();
    gPad->SetGrid(1);
    hmass_add[2][ipt]->SetXTitle("M_{ee} [GeV]");
    hmass_add[2][ipt]->SetYTitle("counts");
    hmass_add[2][ipt]->SetAxisRange(0,4);
    hmass_add[2][ipt]->SetLineColor(4);
    hmass_org[3][ipt]->SetLineColor(6);
    hmass_add[2][ipt]->SetMinimum(1);
    hmass_add[2][ipt]->Draw();
    hmass_org[3][ipt]->Draw("same");
  }
    



  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(npt,2);
  for(int ipt=0;ipt<npt;ipt++){
    c1->cd(ipt+1);

    gPad->SetGrid(1);
    hdiv[ipt]->SetXTitle("M_{ee}");
    hdiv[ipt]->SetYTitle("R(Mixed unlike/Mixed like)");
    hdiv[ipt]->Draw();

    c1->cd(ipt+1+npt);

    gPad->SetGrid(1);
    gPad->SetLogy();
    hmass_add[0][ipt]->SetXTitle("M_{ee} [GeV]");
    hmass_add[0][ipt]->SetYTitle("counts");
    hmass_add[0][ipt]->SetAxisRange(0,4);
    hmass_add[0][ipt]->SetLineColor(1);
    hmass_add[1][ipt]->SetLineColor(2);

    hmass_add[0][ipt]->SetMinimum(1);
    hmass_add[0][ipt]->Draw();
    hmass_add[1][ipt]->Draw("same");
    hmass_add[0][ipt]->Draw("same");
  }
}
