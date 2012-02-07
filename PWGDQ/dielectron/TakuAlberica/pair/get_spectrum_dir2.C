void get_spectrum_dir2(int icent=1){
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  //gStyle->SetOptFit(0);
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

  const int ncentbin=9;
  //  int cent_low[ncentbin] ={  0,  0, 20,  40, 40,  0,  0,  0,  30};
  //  int cent_high[ncentbin]={100, 10, 40, 100, 80, 80, 40, 30, 100};

  //  int cent_low_bin[ncentbin]={0, 0, 2, 4, 4, 0, 0, 0, 3 };
  //  int cent_high_bin[ncentbin]={10, 1, 4, 10, 8, 8, 4, 3, 10};

  int cent_low[ncentbin] ={ 0, 10, 20, 30};
  int cent_high[ncentbin]={10, 20, 30, 50};

  int cent_low_bin[ncentbin] ={0, 1, 2, 3};
  int cent_high_bin[ncentbin]={1, 2, 3, 5};

  
  TH2D *hmasspt_inc[11][7];
  TH1D *hmasspt[4][5];
  TH1D *hmasspt_ls[4][5];
  TH1D *hdiv[5];

  char outname[100];
  //sprintf(outname,"spec_histo_loose_eid_pt04_2_cent%d_no_gst.root",icent);
  sprintf(outname,"spec_histo_pt04_2_cent%d_no_gst.root",icent);
  //sprintf(outname,"output/spec_histo_loose_eid_2_pt04_cent%d_6.root",icent);
  //sprintf(outname,"spec_histo_loose_eid_pt04_cent%d_conv2.root",icent);
  char name[100];
  //  TFile *fin = new TFile("result_histo_nmix40_pt04.root");
  // TFile *fin = new TFile("result_histo_loose_eid_nmix10_conv2_pt04.root");
  //TFile *fin = new TFile("output/result_histo_loose_eid_nmix40_rev1_6.root");
  //TFile *fin = new TFile("result_histo_loose_eid_nmix40_pt04_no_gst.root");
  TFile *fin = new TFile("result_trig1_cut1.root");
  //TFile *fin = new TFile("result_histo_nmix20_pt04.root");
  for(int i=0;i<11;i++){
    for(int j=0;j<7;j++){
      sprintf(name,"hmasspt_weight_cent%d_pair%d",i,j);
      //sprintf(name,"hmasspt_cent%d_pair%d",i,j);
      hmasspt_inc[i][j] = (TH2D*)fin->Get(name);
      //cout<<i<<" "<<j<<" "<<hmasspt_inc[i][j]->GetEntries()<<endl;
    }
  }

  TH1D *hpt = (TH1D*)hmasspt_inc[0][0]->ProjectionY("hpt");

  TH1F *hstat = (TH1F*)fin->Get("hCentrality");
  float nevents = hstat->Integral(cent_low[icent], cent_high[icent]-1);

  cout<<nevents<<endl;

  const int nptbin=4;
  float pt_low[nptbin] =  {0.5, 1.0, 1.5, 2.0};
  float pt_high[nptbin] = {1.0, 1.5, 2.0, 5.0};
  double mass_r_ptbin[nptbin+2]={0, 0.5, 1.0, 1.5, 2.0, 4.0};
  
  for(int i=0;i<4;i++){
    for(int j=0;j<nptbin;j++){
      sprintf(name,"hmass2_w_pair%d_pt%d", i, j);
      //sprintf(name,"hmass2_pair%d_pt%d", i, j);
      hmasspt[i][j] = new TH1D(name, name, 500, 0, 5);
      hmasspt[i][j]->Sumw2();

      sprintf(name,"hmass2_ls_w_pair%d_pt%d", i, j);
      //sprintf(name,"hmass2_pair%d_pt%d", i, j);
      hmasspt_ls[i][j] = new TH1D(name, name, 500, 0, 5);
      hmasspt_ls[i][j]->Sumw2();
    }
  }

  /// first []
  /// 0-->unlike
  /// 1-->like
  /// 2-->mixed unlike
  /// 3-->mixed like 

  /// second [] --> pt
  /// third [] --> ee or pp
  float npairs[4][nptbin][2];

  for(int j=0;j<nptbin;j++){
    int binl = hpt->FindBin(pt_low[j]);
    int binh = hpt->FindBin(pt_high[j]);

    npairs[0][j][0]=0; npairs[0][j][1]=0;
    npairs[1][j][0]=0; npairs[1][j][1]=0;
    npairs[2][j][0]=0; npairs[2][j][1]=0;
    npairs[3][j][0]=0; npairs[3][j][1]=0;



    for(int i=cent_low_bin[icent]; i<cent_high_bin[icent]; i++){

      ///////////////////////////////////////////
      ///// pair0
      ////////////////////////////////////////////
      sprintf(name,"hmass_%d_%d_0", i, j);
      hmasspt[0][j]->Add((TH1D*)hmasspt_inc[i][0]->ProjectionX(name,binl, binh-1));
      //      npairs[0][j][0] += hmasspt_inc[i][0]->ProjectionX(name,binl, binh-1)->Integral();
      npairs[0][j][0] += hmasspt_inc[i][0]->Integral();
      
      /////////////////////////////////////////
      //// pair2 (like sign)
      ////////////////////////////////////////
      sprintf(name,"hmass_%d_%d_1", i, j);
      hmasspt[2][j]->Add((TH1D*)hmasspt_inc[i][1]->ProjectionX(name,binl, binh-1));
      hmasspt_ls[0][j]->Add((TH1D*)hmasspt_inc[i][1]->ProjectionX(name,binl, binh-1));
      //npairs[2][j][0] += hmasspt_inc[i][1]->ProjectionX(name,binl, binh-1)->Integral();
      npairs[2][j][0] += hmasspt_inc[i][1]->Integral();

      sprintf(name,"hmass_%d_%d_2", i, j);
      hmasspt[2][j]->Add((TH1D*)hmasspt_inc[i][2]->ProjectionX(name,binl, binh-1));
      //npairs[2][j][1] += hmasspt_inc[i][2]->ProjectionX(name,binl, binh-1)->Integral();
      hmasspt_ls[1][j]->Add((TH1D*)hmasspt_inc[i][2]->ProjectionX(name,binl, binh-1));
      npairs[2][j][1] += hmasspt_inc[i][2]->Integral();

      /////////////////////////////////////////
      //// pair1 (mixed unlike sign)
      /////////////////////////////////////////
      sprintf(name,"hmass_%d_%d_1", i, j);
      hmasspt[1][j]->Add((TH1D*)hmasspt_inc[i][3]->ProjectionX(name,binl, binh-1));
      //npairs[1][j][0] += hmasspt_inc[i][3]->ProjectionX(name,binl, binh-1)->Integral();
      npairs[1][j][0] += hmasspt_inc[i][3]->Integral();

      sprintf(name,"hmass_%d_%d_2", i, j);
      hmasspt[1][j]->Add((TH1D*)hmasspt_inc[i][4]->ProjectionX(name,binl, binh-1));
      //      npairs[1][j][0] += hmasspt_inc[i][4]->ProjectionX(name,binl, binh-1)->Integral();
      npairs[1][j][0] += hmasspt_inc[i][4]->Integral();

      /////////////////////////////////////////
      //// pair3 (mixed like sign)
      ////////////////////////////////////////
      sprintf(name,"hmass_%d_%d_1", i, j);
      //      hmasspt[3][j]->Add((TH1D*)hmasspt_inc[i][5]->ProjectionX(name,binl, binh-1));
      hmasspt_ls[2][j]->Add((TH1D*)hmasspt_inc[i][5]->ProjectionX(name,binl, binh-1));
      //npairs[3][j][0] += hmasspt_inc[i][5]->ProjectionX(name,binl, binh-1)->Integral();
      npairs[3][j][0] += hmasspt_inc[i][5]->Integral();

      sprintf(name,"hmass_%d_%d_2", i, j);
      //hmasspt[3][j]->Add((TH1D*)hmasspt_inc[i][6]->ProjectionX(name,binl, binh-1));
      hmasspt_ls[3][j]->Add((TH1D*)hmasspt_inc[i][6]->ProjectionX(name,binl, binh-1));
      //npairs[3][j][1] += hmasspt_inc[i][6]->ProjectionX(name,binl, binh-1)->Integral();
      npairs[3][j][1] += hmasspt_inc[i][6]->Integral();

    }
  }

  /// for mixed like sign pairs 
  /// entry should be same as like sign pairs 
  /// re-add the hist for mixed like sign pairs using weights  
  for(int j=0;j<nptbin;j++){
    float frac1 = (npairs[2][j][0]*(npairs[3][j][0]+npairs[3][j][1]))/(npairs[3][j][0]*(npairs[2][j][0]+npairs[2][j][1]));
    float frac2 = (npairs[2][j][1]*(npairs[3][j][0]+npairs[3][j][1]))/(npairs[3][j][1]*(npairs[2][j][0]+npairs[2][j][1]));
    hmasspt[3][j]->Add( hmasspt_ls[2][j], frac1); ///mixed like sign 
    hmasspt[3][j]->Add( hmasspt_ls[3][j], frac2); // mixed like sign 
    //hmasspt[3][j]->Add( hmasspt_ls[2][j], 1);
    //hmasspt[3][j]->Add( hmasspt_ls[3][j], 1);
  }


  for(int i=0;i<4;i++){
    for(int j=0;j<nptbin;j++){
      hmasspt[i][j]->Scale(1.0/nevents);
      hmasspt_ls[i][j]->Scale(1.0/nevents);
      cout<<hmasspt[i][j]->Integral()<<endl;
      cout<<npairs[i][j][0]<<" "<<npairs[i][j][1]<<endl;
    }
  }

  TCanvas *c0 = new TCanvas("c0","c0",1500,500);
  c0->Divide(3,1);
  for(int iptbin=1;iptbin<nptbin;iptbin++){
    c0->cd(iptbin);
    gPad->SetGrid(1);
    gPad->SetLogy();
    hmasspt[0][iptbin]->SetXTitle("M_{ee}[GeV/c^{2}]");
    hmasspt[0][iptbin]->SetYTitle("N_{ee}/N_{evt}");
    hmasspt[0][iptbin]->Draw();
  }


  TCanvas *c1 = new TCanvas("c1","c1",1500,1000);
  c1->Divide(3,3);

  TF1 *f1 = new TF1("f1","pol0",0.05,2);

  double r[5];

  for(int iptbin=1;iptbin<nptbin;iptbin++){
    sprintf(name,"hdiv_pt%d",iptbin);
    hdiv[iptbin] = new TH1D(name, name, 500,0,5);
    hdiv[iptbin]->SetXTitle("M_{ee}[GeV/c^{2}]");
    hdiv[iptbin]->SetYTitle("mixed unlike/mixed like");
    hdiv[iptbin]->Sumw2();
    //hmasspt[3][iptbin]->Scale(npairs[1][iptbin][0]/(npairs[3][iptbin][0]+npairs[3][iptbin][1]));
    hdiv[iptbin]->Divide(hmasspt[1][iptbin], hmasspt[3][iptbin]); ///mixed unlike/mixed like
    /*
    for(int ii=0;ii<hdiv[iptbin]->GetNbinsX(); ii++){
      double vent1 = hmasspt[1][iptbin]->GetBinContent(ii+1);
      double e_vent1 = hmasspt[1][iptbin]->GetBinError(ii+1);


      double vent2 = hmasspt_ls[2][iptbin]->GetBinContent(ii+1);
      double e_vent2 = hmasspt_ls[2][iptbin]->GetBinError(ii+1);
      double vent3 = hmasspt_ls[3][iptbin]->GetBinContent(ii+1);
      double e_vent3 = hmasspt_ls[3][iptbin]->GetBinError(ii+1);



      double val1 = vent1; 
      double val2 = 2*sqrt(vent2*vent3);

      if(val2>0){
	hdiv[iptbin]->SetBinContent(ii+1, val1/val2);
	double e_val1 = e_vent1;
	double e_val2 = val2/2*sqrt(pow(e_vent2/vent2,2)+pow(e_vent3/vent3,2));
	hdiv[iptbin]->SetBinError(ii+1, val1/val2*sqrt(pow(e_val1/val1,2)+pow(e_val2/val2,2)));
      }
    }
    */   

    c1->cd(iptbin);
    gPad->SetGrid(1);
    hmasspt[1][iptbin]->SetXTitle("M_{ee}[GeV/c^{2}]");
    hmasspt[0][iptbin]->SetYTitle("N_{ee}/N_{evt}");
    hmasspt[1][iptbin]->Draw();
    hmasspt[3][iptbin]->SetLineColor(2);
    hmasspt[3][iptbin]->Draw("same");

    c1->cd(iptbin+3);
    gPad->SetGrid(1);
    //hdiv[iptbin]->Rebin(4); hdiv[iptbin]->Scale(1.0/4); 

    hdiv[iptbin]->SetAxisRange(0,3);
    hdiv[iptbin]->SetMinimum(0);
    hdiv[iptbin]->SetMaximum(2);
    hdiv[iptbin]->Draw();
    ///hdiv[iptbin]->Fit(f1,"R","");

    //r[iptbin] = f1->GetParameter(0);
    //int bin1 =  hmasspt[2][iptbin]->FindBin(0.1);
    //int bin2 =  hmasspt[2][iptbin]->FindBin(1.0);
    //r[iptbin] = hmasspt[2][iptbin]->Integral(bin1, bin2)/hmasspt[3][iptbin]->Integral(bin1, bin2);
    //cout<<" normalization  factor "<<r[iptbin]<<endl;

    /*
    for(int k=0;k<hmasspt[1][iptbin]->GetNbinsX();k++){
      float w = f1->Eval(hmasspt[1][iptbin]->GetBinCenter(k+1));
      float rr = hmasspt[1][iptbin]->GetBinContent(k+1);
      float err = hmasspt[1][iptbin]->GetBinError(k+1);
      hmasspt[1][iptbin]->SetBinContent(k+1, rr*w);
      hmasspt[1][iptbin]->SetBinError(k+1, err*w);
    }
    */
    //hmasspt[1][iptbin]->Scale(r[iptbin]);
    hmasspt[2][iptbin]->Multiply(hdiv[iptbin]);
    //c/out<<npairs[1][iptbin][0]/(2*sqrt(npairs[3][iptbin][0]*npairs[3][iptbin][1]))<<endl;
    //cout<<" aaa "<<(2*sqrt(npairs[2][iptbin][0]*npairs[2][iptbin][1]))<<" "<<npairs[2][iptbin][0]+npairs[2][iptbin][1]<<" "<<hmasspt[2][iptbin]->Integral()*nevents<<endl;
    //hmasspt[2][iptbin]->Scale((2*sqrt(npairs[2][iptbin][0]*npairs[2][iptbin][1])/(npairs[2][iptbin][0]+npairs[2][iptbin][1])));
    //hmasspt[2][iptbin]->Scale(npairs[0][iptbin][0]/(2*sqrt(npairs[2][iptbin][0]*npairs[2][iptbin][1])));

   
    for(int ii=0;ii<hdiv[iptbin]->GetNbinsX(); ii++){    
      double d_val1 = 2*sqrt(hmasspt_ls[0][iptbin]->GetBinContent(ii+1)*hmasspt_ls[1][iptbin]->GetBinContent(ii+1));
      if(d_val1>0){
	double d_eval1 = d_val1/2*sqrt(pow(hmasspt_ls[0][iptbin]->GetBinError(ii+1)/hmasspt_ls[0][iptbin]->GetBinContent(ii+1),2)+
					pow(hmasspt_ls[1][iptbin]->GetBinError(ii+1)/hmasspt_ls[1][iptbin]->GetBinContent(ii+1),2));
	//	double accep = hdiv[iptbin]->GetBinContent(ii+1)-hdiv[iptbin]->GetBinError(ii+1);
	double accep = hdiv[iptbin]->GetBinContent(ii+1);
	hmasspt[2][iptbin]->SetBinContent(ii+1, d_val1*accep);
	hmasspt[2][iptbin]->SetBinError(ii+1, d_eval1*accep);
      }
    }
   
    /*
    double norm1=0;
    double norm2=0;
    for(int ii=0;ii<hmasspt[2][iptbin]->GetNbinsX(); ii++){    
      double d_bg1 = hmasspt[1][iptbin]->GetBinContent(ii+1);
      double d_bg_pp = hmasspt_ls[2][iptbin]->GetBinContent(ii+1);
      double d_bg_mm = hmasspt_ls[3][iptbin]->GetBinContent(ii+1);
      double d_fg_pp = hmasspt_ls[0][iptbin]->GetBinContent(ii+1);
      double d_fg_mm = hmasspt_ls[1][iptbin]->GetBinContent(ii+1);
      double d_e_fg_pp = hmasspt_ls[0][iptbin]->GetBinError(ii+1);
      double d_e_fg_mm = hmasspt_ls[1][iptbin]->GetBinError(ii+1);

      if(d_fg_pp>0 && d_fg_mm>0 && d_bg_pp>0 && d_bg_mm>0){
	//double d_val = 2*sqrt(d_fg_pp*d_fg_mm)*d_bg1/(2*sqrt(d_bg_pp*d_bg_mm));
	//double d_e_val = sqrt(d_fg_pp*d_fg_mm)*sqrt(pow(d_e_fg_pp/d_fg_pp,2)+pow(d_e_fg_mm/d_fg_mm,2))*d_bg1/(2*sqrt(d_bg_pp*d_bg_mm));
	double d_val = 2*(d_fg_pp*(0.5*d_bg1/d_bg_pp)*sqrt(npairs[2][iptbin][0])+d_fg_mm*(0.5*d_bg1/d_bg_mm)*sqrt(npairs[2][iptbin][1]))/(sqrt(npairs[2][iptbin][0])+sqrt(npairs[2][iptbin][1]));
	double d_e_val = d_e_fg_pp*(0.5*d_bg1/d_bg_pp)+d_e_fg_mm*(0.5*d_bg1/d_bg_mm);

	hmasspt[2][iptbin]->SetBinContent(ii+1, d_val);
	hmasspt[2][iptbin]->SetBinError(ii+1, d_e_val);

      }
    }
    */
    //hmasspt[2][iptbin]->Scale(hmasspt[0][iptbin]->Integral()/(2*sqrt(norm1*norm2)));
    cout<<" entry ="<<hmasspt[0][iptbin]->Integral()<<" "<<hmasspt[2][iptbin]->Integral()<<" "<<hmasspt_ls[0][iptbin]->Integral()<<" "<<hmasspt_ls[1][iptbin]->Integral()<<endl;
    c1->cd(iptbin+6);
    hmasspt[0][iptbin]->Draw();
    hmasspt[2][iptbin]->SetLineColor(2);
    hmasspt[2][iptbin]->Draw("same");
  }


  const int nmassbin=8;
  float mass_low[nmassbin] ={0,    0.05, 0.09, 0.14, 0.20, 0.30, 0.1, 0.09};
  float mass_high[nmassbin]={0.03, 0.09, 0.14, 0.20, 0.30, 0.40, 0.3, 0.30};
  float nsig[5][nmassbin];
  float nbg[5][nmassbin];
  float ensig[5][nmassbin];
  float enbg[5][nmassbin];

  for(int ii=0;ii<nptbin; ii++){
    for(int j=0;j<nmassbin; j++){
      int binlow = hmasspt[0][ii]->FindBin(mass_low[j]);
      int binhigh = hmasspt[0][ii]->FindBin(mass_high[j]);
      nsig[ii][j] = hmasspt[0][ii]->Integral(binlow, binhigh-1);      
      nbg[ii][j] = hmasspt[2][ii]->Integral(binlow, binhigh-1);
      ensig[ii][j]=0;
      enbg[ii][j]=0;
      for(int ib=binlow; ib<binhigh;ib++){
	ensig[ii][j] += pow(hmasspt[0][ii]->GetBinError(ib),2);
	enbg[ii][j] += pow(hmasspt[2][ii]->GetBinError(ib),2);
      }

      ensig[ii][j] = sqrt(ensig[ii][j]);
      enbg[ii][j] = sqrt(enbg[ii][j]);


      float r1 = nsig[ii][j]-nbg[ii][j];
      float r0 =  nsig[ii][0]-nbg[ii][0];
      float er1 = sqrt(pow(ensig[ii][j],2)+pow(enbg[ii][j],2));
      float er0 = sqrt(pow(ensig[ii][0],2)+pow(enbg[ii][0],2));


      if(r0>0 && r1>0){
	cout<<" pt : "<<ii<<" mass "<<mass_low[j]<<"-"<<mass_high[j]<<" --> Nsig = "<<nsig[ii][j]<<", Nbg = "<<nbg[ii][j]
	    <<" : Nsig-Nbg "<<nsig[ii][j]-nbg[ii][j]<<" +- "<<sqrt(pow(ensig[ii][j],2)+pow(enbg[ii][j],2))<<" ("<<ensig[ii][j]<<" "<<enbg[ii][j]<<")"
	    <<" : Ratio = "<<r1/r0<<" +- "<<r1/r0*sqrt(pow(er1/r1,2)+pow(er0/r0,2))
	    <<endl;
      }else{
	cout<<" pt : "<<ii<<" mass "<<mass_low[j]<<"-"<<mass_high[j]<<" --> Nsig = "<<nsig[ii][j]<<", Nbg = "<<nbg[ii][j]
	    <<" : Nsig-Nbg "<<nsig[ii][j]-nbg[ii][j]<<" +- "<<sqrt(pow(ensig[ii][j],2)+pow(enbg[ii][j],2))<<" ("<<ensig[ii][j]<<" "<<enbg[ii][j]<<")"
	    <<endl;
      }
    }
  }


  /// for histogram
  TH1D *hsub[5];
  const int nbins = 55;
  double x[nbins+1]={0, 0.01, 0.02, 0.03, 0.05, 0.09, 0.14, 0.20, 0.30, 0.40,
		     0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 
		     1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 
		     2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 
		     3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 
		     4.5, 4.6, 4.7, 4.8, 4.9, 5.0};

//  const int nbins = 53;
//  double x[nbins+1]={0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.20, 0.40, 
//		     0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 
//		     1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 
//		     2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 
//		     3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 
//		     4.5, 4.6, 4.7, 4.8, 4.9, 5.0};

  int narrows[nptbin];
  TArrow *Arrow[nptbin][10];

  for(int i=0; i<nptbin; i++){
    narrows[i]=0;


    sprintf(name,"hsub_mass_pt%d", i);
    hsub[i] = new TH1D(name, name, nbins, x);

    for(int j=0; j<nbins; j++){
      int binl = hmasspt[0][i]->FindBin(x[j]);
      int binh = hmasspt[0][i]->FindBin(x[j+1]);
      float dw = x[j+1]-x[j];

      float ent = 0;
      float e_ent = 0;
      float sig = 0;
      float bg = 0;
      for(int k=binl; k<binh; k++){
	sig += hmasspt[0][i]->GetBinContent(k);
	bg += hmasspt[2][i]->GetBinContent(k);
	ent += (hmasspt[0][i]->GetBinContent(k)-hmasspt[2][i]->GetBinContent(k));
	e_ent += (pow(hmasspt[0][i]->GetBinError(k),2)+
		  pow(hmasspt[2][i]->GetBinError(k),2));
      }
      e_ent = sqrt(e_ent);
      hsub[i]->SetBinContent(j+1, ent/dw);
      hsub[i]->SetBinError(j+1, e_ent/dw);

      if(sig<bg && sig>0 && x[j]<0.5){ // calculate 90% confidence level
	double conf90=0;
	sig = sig*nevents;
	bg = bg*nevents;
	TH1F *htmp = new TH1F("htmp","htmp",1000,-2*bg, 2*bg);
	for(int kk=0;kk<100000;kk++){
	  double a1 = gRandom->PoissonD(sig);
	  double a2 = gRandom->PoissonD(bg);
	  htmp->Fill(a1-a2); 
	}
	int bin0 = htmp->FindBin(0);
	float ent0 = htmp->Integral(bin0, 1000);
	htmp->Scale(1.0/ent0);
	//htmp->Scale(1.0/100000);
	int bincent=0;
	for(int ib=bin0;ib<1000;ib++){
	  if(htmp->Integral(bin0, ib)>0.9){
	    bincent = ib;
	    break;
	  }
	}
	conf90 = htmp->GetBinCenter(bincent);
	cout<<i<<" "<<hsub[i]->GetBinCenter(j+1)<<" : 90% conf: "<<conf90<<" from sig :"<<sig<<" bg:"<<bg<<endl;
	Arrow[i][narrows[i]] = new TArrow(hsub[i]->GetBinCenter(j+1), 1.0e-07, hsub[i]->GetBinCenter(j+1), conf90/(dw*nevents),0.01,"<|");
	narrows[i]++;
	htmp->Clear();
      }// calculate 90% confidence level

    }
  }

  TH1D *hsub_clone[5];

  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  c3->Divide(4,1);

  c3->cd(1);
  gPad->SetGrid(1);  gPad->SetLogy();
  TH2F *hw = new TH2F("hw","hw",10,0,0.5,10,1.0e-07, 10);
  hw->Draw();
  hw->SetXTitle("M_{ee} [GeV]");
  hw->SetYTitle("dN_{ee}/dM_{ee}");
  c3->cd(2);
  gPad->SetGrid(1);  gPad->SetLogy();
  hw->Draw();
  c3->cd(3);
  gPad->SetGrid(1);  gPad->SetLogy();
  hw->Draw();
  c3->cd(4);
  gPad->SetGrid(1);  gPad->SetLogy();
  hw->Draw();

  for(int i=0; i<nptbin; i++){
    hsub_clone[i] = (TH1D*)hsub[i]->Clone();
    sprintf(name,"%s_cln",hsub[i]->GetName());
    hsub_clone[i]->SetName(name);
    hsub_clone[i]->SetLineColor(i+1);
    hsub_clone[i]->SetMarkerColor(i+1);
    hsub_clone[i]->SetMarkerStyle(20);
  }


  //  int bin = hsub_clone[1]->FindBin(0.03);
  //  float ent = hsub_clone[1]->Integral(0,bin);
  for(int i=0; i<nptbin; i++){
    //float ent1 = hsub_clone[i]->Integral(0,bin);
    //hsub_clone[i]->Scale(ent/ent1);
    c3->cd(i+1);
    hsub_clone[i]->Draw("same");
    for(int k=0;k<narrows[i];k++){
      Arrow[i][k]->SetLineColor(i+1);
      Arrow[i][k]->Draw();
    }
  }

  TCanvas *c5 = new TCanvas("c5","c5",1200,600);
  c5->Divide(3,1);

  c5->cd(1);
  gPad->SetGrid(1);  gPad->SetLogy();
  TH2F *hw2 = new TH2F("hw2","hw2",10,0,0.5,10,1.0e-07, 10);
  hw2->Draw();
  hw2->SetXTitle("M_{ee} [GeV]");
  hw2->SetYTitle("1/N_{evt}dN_{ee}/dM_{ee}");
  c5->cd(2);
  gPad->SetGrid(1);  gPad->SetLogy();
  hw2->Draw();
  c5->cd(3);
  gPad->SetGrid(1);  gPad->SetLogy();
  hw2->Draw();

  TH1D *hmass_bin[2][nptbin];
  for(int i=1; i<nptbin; i++){
    c5->cd(i);
    hmass_bin[0][i] = (TH1D*)hmasspt[0][i]->Clone();
    hmass_bin[1][i] = (TH1D*)hmasspt[2][i]->Clone();


    float binw = hmass_bin[0][i]->GetBinWidth(10);
    hmass_bin[0][i]->Scale(1/binw);
    hmass_bin[1][i]->Scale(1/binw);

    hmass_bin[1][i]->SetLineColor(2);
    hmass_bin[0][i]->Draw("same");    
    hmass_bin[1][i]->Draw("same");
    hsub_clone[i]->Draw("same");
    for(int k=0;k<narrows[i];k++){
      Arrow[i][k]->SetLineColor(i+1);
      Arrow[i][k]->Draw();
    }
  }






  //  TLegend *leg = new TLegend(0.5,0.7,0.875,0.875);
  //  leg->SetFillColor(10);
  //  if(nptbin==4){
  //    leg->AddEntry(hsub_clone[0],"0.5<p_{T}<1.0 GeV","pl");
  //    leg->AddEntry(hsub_clone[1],"1.0<p_{T}<1.5 GeV","pl");
  //    leg->AddEntry(hsub_clone[2],"1.5<p_{T}<2.0 GeV","pl");
  //    leg->AddEntry(hsub_clone[3],"2.0<p_{T}","pl");
  //  }else{
  //  }
  //  leg->Draw("same");

  /////draw the results 
  TH1D *hmass_cocktail[4];
  TH1D *hmass_cocktail_comp[4][6];
  TFile *fin1 = new TFile("../../ana2/cocktail/cocktail_allpt_pt04_no_smearing.root");
  fin1->cd();
  hmass_cocktail[0] = (TH1D*)fin1->Get("hmass_all");
  hmass_cocktail[0]->SetName("hmass_all_pt0");
  hmass_cocktail_comp[0][0] =  (TH1D*)fin1->Get("hmass_3");
  hmass_cocktail_comp[0][0]->SetName("hmass_pt0_eta_prime");
  hmass_cocktail_comp[0][1] =  (TH1D*)fin1->Get("hmass_4");
  hmass_cocktail_comp[0][1]->SetName("hmass_pt0_rho");
  hmass_cocktail_comp[0][2] =  (TH1D*)fin1->Get("hmass_7");
  hmass_cocktail_comp[0][2]->SetName("hmass_pt0_pi0");
  hmass_cocktail_comp[0][3] =  (TH1D*)fin1->Get("hmass_8");
  hmass_cocktail_comp[0][3]->SetName("hmass_pt0_eta");
  hmass_cocktail_comp[0][4] =  (TH1D*)fin1->Get("hmass_9");
  hmass_cocktail_comp[0][4]->SetName("hmass_pt0_omega");
  hmass_cocktail_comp[0][5] =  (TH1D*)fin1->Get("hmass_10");
  hmass_cocktail_comp[0][5]->SetName("hmass_pt0_phi");
  

  //TFile *fin2 = new TFile("./cocktail//cocktail_pt10_15.root");
  TFile *fin2 = new TFile("../../ana2/cocktail//cocktail_pt04_pairpt_10_15_no_smearing.root");
  fin2->cd();
  hmass_cocktail[1] = (TH1D*)fin2->Get("hmass_all");
  hmass_cocktail[1]->SetName("hmass_all_pt1");
  hmass_cocktail_comp[1][0] =  (TH1D*)fin1->Get("hmass_3");
  hmass_cocktail_comp[1][0]->SetName("hmass_pt1_eta_prime");
  hmass_cocktail_comp[1][1] =  (TH1D*)fin1->Get("hmass_4");
  hmass_cocktail_comp[1][1]->SetName("hmass_pt1_rho");
  hmass_cocktail_comp[1][2] =  (TH1D*)fin1->Get("hmass_7");
  hmass_cocktail_comp[1][2]->SetName("hmass_pt1_pi0");
  hmass_cocktail_comp[1][3] =  (TH1D*)fin1->Get("hmass_8");
  hmass_cocktail_comp[1][3]->SetName("hmass_pt1_eta");
  hmass_cocktail_comp[1][4] =  (TH1D*)fin1->Get("hmass_9");
  hmass_cocktail_comp[1][4]->SetName("hmass_pt1_omega");
  hmass_cocktail_comp[1][5] =  (TH1D*)fin1->Get("hmass_10");
  hmass_cocktail_comp[1][5]->SetName("hmass_pt1_phi");



  //TFile *fin3 = new TFile("./cocktail//cocktail_pt15_20.root");
  //TFile *fin3 = new TFile("./cocktail//cocktail_pt15_20_no_smearing.root");
  TFile *fin3 = new TFile("../../ana2/cocktail//cocktail_pt04_pairpt_15_20_no_smearing.root");
  fin3->cd();
  hmass_cocktail[2] = (TH1D*)fin3->Get("hmass_all");
  hmass_cocktail[2]->SetName("hmass_all_pt2");
  hmass_cocktail_comp[2][0] =  (TH1D*)fin1->Get("hmass_3");
  hmass_cocktail_comp[2][0]->SetName("hmass_pt2_eta_prime");
  hmass_cocktail_comp[2][1] =  (TH1D*)fin1->Get("hmass_4");
  hmass_cocktail_comp[2][1]->SetName("hmass_pt2_rho");
  hmass_cocktail_comp[2][2] =  (TH1D*)fin1->Get("hmass_7");
  hmass_cocktail_comp[2][2]->SetName("hmass_pt2_pi0");
  hmass_cocktail_comp[2][3] =  (TH1D*)fin1->Get("hmass_8");
  hmass_cocktail_comp[2][3]->SetName("hmass_pt2_eta");
  hmass_cocktail_comp[2][4] =  (TH1D*)fin1->Get("hmass_9");
  hmass_cocktail_comp[2][4]->SetName("hmass_pt2_omega");
  hmass_cocktail_comp[2][5] =  (TH1D*)fin1->Get("hmass_10");
  hmass_cocktail_comp[2][5]->SetName("hmass_pt2_phi");



  //TFile *fin4 = new TFile("./cocktail//cocktail_pt20_no_smearing.root");
  //TFile *fin4 = new TFile("./cocktail//cocktail_pt20.root");
TFile *fin4 = new TFile("../../ana2/cocktail//cocktail_pt04_pairpt_20_no_smearing.root");
  fin4->cd();
  hmass_cocktail[3] = (TH1D*)fin4->Get("hmass_all");
  hmass_cocktail[3]->SetName("hmass_all_pt3");
  hmass_cocktail_comp[3][0] =  (TH1D*)fin1->Get("hmass_3");
  hmass_cocktail_comp[3][0]->SetName("hmass_pt3_eta_prime");
  hmass_cocktail_comp[3][1] =  (TH1D*)fin1->Get("hmass_4");
  hmass_cocktail_comp[3][1]->SetName("hmass_pt3_rho");
  hmass_cocktail_comp[3][2] =  (TH1D*)fin1->Get("hmass_7");
  hmass_cocktail_comp[3][2]->SetName("hmass_pt3_pi0");
  hmass_cocktail_comp[3][3] =  (TH1D*)fin1->Get("hmass_8");
  hmass_cocktail_comp[3][3]->SetName("hmass_pt3_eta");
  hmass_cocktail_comp[3][4] =  (TH1D*)fin1->Get("hmass_9");
  hmass_cocktail_comp[3][4]->SetName("hmass_pt3_omega");
  hmass_cocktail_comp[3][5] =  (TH1D*)fin1->Get("hmass_10");
  hmass_cocktail_comp[3][5]->SetName("hmass_pt3_phi");



  for(int i=0;i<4;i++){
    hmass_cocktail[i]->Rebin(4);
    float ent = 0;
    for(int j=0;j<hsub_clone[i]->GetNbinsX();j++){
      if(hsub_clone[i]->GetBinCenter(j+1)<0.03){
	ent += hsub_clone[i]->GetBinWidth(j+1)*hsub_clone[i]->GetBinContent(j+1);
      }
    }
    float ent2=0;
    for(int j=0;j<hmass_cocktail[i]->GetNbinsX();j++){
      if(hmass_cocktail[i]->GetBinCenter(j+1)<0.03){
	ent2 += hmass_cocktail[i]->GetBinWidth(j+1)*hmass_cocktail[i]->GetBinContent(j+1);
      }
    }
    hmass_cocktail[i]->Scale(ent/ent2);
    hmass_cocktail[i]->SetLineColor(6);
    hmass_cocktail[i]->SetLineWidth(2);
    c3->cd(i+1);
    hmass_cocktail[i]->Draw("same,l");
    if(i>=1){
      c5->cd(i);
      hmass_cocktail[i]->Draw("same,l");
    }

    for(int j=0;j<6;j++){
      hmass_cocktail_comp[i][j]->Rebin(4);
      hmass_cocktail_comp[i][j]->Scale(ent/ent2*0.2);
      hmass_cocktail_comp[i][j]->SetLineColor(1);
      //hmass_cocktail_comp[i][j]->Draw("same");
    }
  }

  float r_cocktail[nptbin];
  cout<<" ***************** "<<endl;
  cout<<" mass ratio ("<<mass_low[nmassbin-1]<<"-"<<mass_high[nmassbin-1]<<"/0-30) in cocktail"<<endl;
  TH1F *hmass_r_ck = new TH1F("hmass_r_cocktail","mass ratio cocktail", nptbin+1, mass_r_ptbin);
  for(int i=1;i<nptbin;i++){
    int bin1=hmass_cocktail[i]->FindBin(0);
    int bin2=hmass_cocktail[i]->FindBin(0.03);
    float ent1 = hmass_cocktail[i]->Integral(bin1, bin2-1);
    bin1=hmass_cocktail[i]->FindBin(mass_low[nmassbin-1]);
    bin2=hmass_cocktail[i]->FindBin(mass_high[nmassbin-1]);
    float ent2 = hmass_cocktail[i]->Integral(bin1, bin2-1);
    cout<<" pt : "<<i<<" "
	<<ent2/ent1
	<<endl;
    int bin = hmass_r_ck->FindBin((0.5*(pt_low[i]+pt_high[i])));
    hmass_r_ck->SetBinContent(bin, ent2/ent1); 
    hmass_r_ck->SetBinError(bin,  ent2/ent1*0.05); 
    r_cocktail[i] = ent2/ent1;
  }

  float r_dir = 0.37; // this is from 90-300MeV

  TH1F *hmass_r = new TH1F("hmass_r","mass ratio", nptbin+1, mass_r_ptbin);
  TH1F *hrgamma = new TH1F("hrgamma","r gmamma", nptbin+1, mass_r_ptbin);
  cout<<" ***************** "<<endl;
  cout<<" mass ratio (100-300/0-30)"<<endl;
  for(int i=1; i<nptbin; i++){
    cout<<" pt : "<<i<<" "
	<<nsig[i][nmassbin-1]-nbg[i][nmassbin-1]<<" +- "<<sqrt(pow(ensig[i][nmassbin-1],2)+pow(enbg[i][nmassbin-1],2))<<" / "
	<<(nsig[i][0]-nbg[i][0])<<" ratio = "
	<<(nsig[i][nmassbin-1]-nbg[i][nmassbin-1])/(nsig[i][0]-nbg[i][0])
	<<" +- "<<sqrt(pow(ensig[i][nmassbin-1],2)+pow(enbg[i][nmassbin-1],2))/(nsig[i][0]-nbg[i][0])
	<<endl;

    int bin = hmass_r->FindBin((0.5*(pt_low[i]+pt_high[i])));
    float r_data = (nsig[i][nmassbin-1]-nbg[i][nmassbin-1])/(nsig[i][0]-nbg[i][0]);
    float er_data = sqrt(pow(ensig[i][nmassbin-1],2)+pow(enbg[i][nmassbin-1],2))/(nsig[i][0]-nbg[i][0]);

    hmass_r->SetBinContent(bin, r_data);
    hmass_r->SetBinError(bin, er_data);

    float r_gamma = (r_data-r_cocktail[i])/(r_dir - r_cocktail[i]);
    float e_r_gamma = (er_data)/(r_dir - r_cocktail[i]);
    hrgamma->SetBinContent(bin, r_gamma);
    hrgamma->SetBinError(bin, e_r_gamma);

  }

  hrgamma->SetXTitle("p_{T} [GeV]");
  hrgamma->SetYTitle("R_{#gamma}^{0-30}");
  hrgamma->SetMarkerColor(icent+1);
  hrgamma->SetLineColor(icent+1);
  hrgamma->SetMarkerStyle(20);
  hrgamma->SetMinimum(0);
  hrgamma->SetMaximum(1);

  TCanvas *c4 = new TCanvas("c4","c4");
  gPad->SetGrid(1);
  hmass_r->SetXTitle("p_{T} [GeV]");
  hmass_r->SetYTitle("R_{data}^{90-300}");
  hmass_r->SetMarkerColor(icent+1);
  hmass_r->SetLineColor(icent+1);
  hmass_r->SetMarkerStyle(20);
  hmass_r->SetMinimum(0);
  hmass_r->SetMaximum(0.5);

  hmass_r_ck->SetMarkerStyle(24);
  hmass_r->Draw(); hmass_r_ck->Draw("same");


  



  TFile *fout = new TFile(outname,"recreate");
  hw->Write();
  for(int i=0; i<nptbin; i++){
    hsub[i]->Write();
    hsub_clone[i]->Write();
  }
  for(int j=0;j<4;j++){
    hmass_cocktail[j]->Write();
    for(int k=0;k<6;k++){
      hmass_cocktail_comp[j][k]->Write();
    }
  }
  hmass_r->Write();
  hmass_r_ck->Write();
  hrgamma->Write();

}
      


  

