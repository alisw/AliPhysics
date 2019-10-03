void PlotJpsiCrsec7TeV(const Int_t SPD=0, const Char_t* met="Arith") {

// plot dsigma/dy and dsigma/dpt, 7 TeV data
//

    //gROOT->SetStyle("Plain");
    //gStyle->SetOptStat(0);

    char *loc = 0x0;

    Float_t BR=0.0594; // branching ratio ...as used to extract the ALICE cross sect.

  Float_t errL=0.08; //luminosity error

  const Int_t Npt=6;
  Float_t Pt[Npt],PtE[Npt],PtEs[Npt];
  Float_t PtErrDo[Npt];
  Float_t PtErrUp[Npt];
  Float_t y1[Npt],y2[Npt],y3[Npt],y1e[Npt],y2e[Npt],y3e[Npt]; //for graphs
  Float_t dsdpt[Npt],dsdptE[Npt],dsdptEs[Npt],dsdptEs1[Npt],dsdptEs2[Npt],dsdptEs3[Npt];
  char filename[200];
  Float_t IntData=0.;
 
  sprintf(filename, Form("dsigdpt_spd%1d_%s.txt",SPD,met));
  ifstream vfile(filename);
  int ipt=0;
  if(loc = gSystem->FindFile(".",filename)){
      while (vfile >> Pt[ipt] >> PtE[ipt] >> dsdpt[ipt] >> dsdptE[ipt] >> dsdptEs[ipt]){ 
	  PtEs[ipt]=0.12; //box width for syst err. plot
	  dsdptEs1[ipt]=TMath::Sqrt(dsdptE[ipt]*dsdptE[ipt]+dsdptEs[ipt]*dsdptEs[ipt]); //all err. (-lumi)
	  dsdptEs2[ipt]=TMath::Sqrt(dsdptE[ipt]*dsdptE[ipt]+dsdptEs[ipt]*dsdptEs[ipt]+errL*errL*dsdpt[ipt]*dsdpt[ipt]); //all err. added (incl. lumi 7%)
	  dsdptEs3[ipt]=errL*dsdpt[ipt]; //Luminosity error
	  ipt++; 
      }
      vfile.close();
  } else {
      fprintf(stderr, "Cannot open %s\n", filename);
  }

  const Int_t Np=ipt;
  const Int_t NpPl=Np-1; //ommit last point
  cout << Np << " points from file " << filename << " ...used for plot: " <<NpPl << endl;

  TVirtualPad *pad = 0x0;

  TGraphErrors *gSpect;
  gSpect = new TGraphErrors(NpPl,Pt,dsdpt,PtE,dsdptE);
  
  TGraphAsymmErrors *gSpect2; //syst. err.
  gSpect2 = new TGraphAsymmErrors(NpPl,Pt,dsdpt,PtEs,PtEs,dsdptEs,dsdptEs);

  TGraphAsymmErrors *gSpect3; //all err. (incl. lumi)
  gSpect3 = new TGraphAsymmErrors(NpPl,Pt,dsdpt,PtE,PtE,dsdptEs2,dsdptEs2);

  TGraphErrors *gSpect1=new TGraphErrors(NpPl,Pt,dsdpt,PtE,dsdptEs1);//all err-lumi
  TGraphErrors *gSpect4=new TGraphErrors(NpPl,Pt,dsdpt,PtEs,dsdptEs3);//lumi err

  Float_t Y1M=1.5*y1[1]; Float_t Y2M=1.5*y2[3]; Float_t Y3M=1.4*y3[3];
  Float_t Y1m=0.5*y1[0]; Float_t Y2m=0.5*y2[0]; Float_t Y3m=0.7*y3[0];

  const Int_t kMarkTyp=20; //marker type
  const Int_t kMarkCol=2; //...and color
  const Float_t kTitSize=0.055; //axis title size
  const Float_t kAsize=0.85*kTitSize; //...and label size
  const Float_t kToffset=0.8;
  const Int_t kFont=42;
  
  TCanvas *c4=new TCanvas("dsigdpt","dsigdpt",20,20,620,620);
  c4->SetTopMargin(0.03);  c4->SetRightMargin(0.03);
  c4->SetLeftMargin(0.15);  c4->SetBottomMargin(0.145);
  pad = c4->cd(1); pad->SetLogy();
  
  TH2D *ho1 = new TH2D("ho1", "ho1", 10, 0, 8, 10, 0.03, 3.); //...just for frame
  ho1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ho1->GetYaxis()->SetTitle("d^{2}#sigma_{J/#psi} /dp_{T}dy (#mub/GeV/c)");
  ho1->GetXaxis()->SetTitleOffset(1.1);
  ho1->GetYaxis()->SetTitleOffset(1.1);
  ho1->GetYaxis()->SetLabelOffset(.01);
  ho1->SetTitleSize(0.06,"XY");ho1->SetTitleFont(42,"XY");
  ho1->SetLabelSize(0.05,"XY");ho1->SetLabelFont(42,"XY");
  ho1->SetTitle("");
  ho1->Draw("");

  gSpect->SetMarkerStyle(kMarkTyp); gSpect->SetMarkerSize(1.5);
  gSpect->SetMarkerColor(kMarkCol); gSpect->SetLineColor(kMarkCol);
  gSpect->SetLineWidth(1.4);
  gSpect->Draw("psame");

  gSpect2->SetMarkerStyle(0);gSpect2->SetMarkerColor(2); gSpect2->SetMarkerSize(0.1);
  gSpect2->SetLineStyle(1); gSpect2->SetLineColor(2); gSpect2->SetLineWidth(1.4);
  gSpect2->SetFillColor(0); gSpect2->SetFillStyle(0);
  gSpect2->Draw("spe2");

  leg1 = new TLegend(.18,0.24,.54,0.35);
  leg1->SetTextFont(kFont);    leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);    leg1->SetMargin(0.16); //separation symbol-text
  leg1->SetEntrySeparation(0.2);    
  leg1->AddEntry(gSpect, "e^{+}e^{-}, |y|<0.9", "p");

  float ptMu[10],pteMu[10],sPtMu[10],sePtMu[10],se1PtMu[10],se2PtMu[10];
  float ScaleMu=1.; 
  float sesPtMu[10],ses1PtMu[10],ses2PtMu[10],ses3PtMu[10], PtEsMu[10];

  sprintf(filename, "dsigdpt_y2.5_4.0.txt");
  ifstream vfile(filename);
  if(loc = gSystem->FindFile(".",filename)) {
      int i0=0;
      while ( vfile >> ptMu[i0] >> pteMu[i0] >> sPtMu[i0] >> sePtMu[i0] >> se1PtMu[i0]>> se2PtMu[i0] ) { 
	  printf("%5.1f %7.2f \n",ptMu[i0],sPtMu[i0]);
	  sPtMu[i0]*=ScaleMu;
	  sePtMu[i0]*=ScaleMu;
	  se1PtMu[i0]*=ScaleMu;
	  se2PtMu[i0]*=ScaleMu;
	  sesPtMu[i0]=TMath::Sqrt(se1PtMu[i0]*se1PtMu[i0]+se2PtMu[i0]*se2PtMu[i0]-errL*errL*sPtMu[i0]*sPtMu[i0]); //add correl. and uncorrel. syst err. and subtract lumi
	  ses2PtMu[i0]=TMath::Sqrt(sePtMu[i0]*sePtMu[i0]+se1PtMu[i0]*se1PtMu[i0]+se2PtMu[i0]*se2PtMu[i0]); //add correl. and uncorrel. syst err. (incl. lumi)
	  ses1PtMu[i0]=TMath::Sqrt(sePtMu[i0]*sePtMu[i0]+sesPtMu[i0]*sesPtMu[i0]); //stst+syst. (-lumi)
	  ses3PtMu[i0]=errL*sPtMu[i0]; //syst err. lumi
	  PtEsMu[i0]=0.12; //fixed width for syst err. plot
	  i0++;	  
      }
      cout <<i0 << " points from file " << filename << endl;
      vfile.close();
      const int Nmu=i0;
      TGraphErrors *gMu = new TGraphErrors(Nmu,ptMu,sPtMu,pteMu,sePtMu); //stat
      TGraphErrors *gMu2 = new TGraphErrors(Nmu,ptMu,sPtMu,PtEsMu,sesPtMu); //syst
      TGraphErrors *gMu3 = new TGraphErrors(Nmu,ptMu,sPtMu,pteMu,ses2PtMu); //all

      TGraphErrors *gMu1 = new TGraphErrors(Nmu,ptMu,sPtMu,pteMu,ses1PtMu); //all-lumi
      TGraphErrors *gMu4 = new TGraphErrors(Nmu,ptMu,sPtMu,PtEsMu,ses3PtMu); //lumi

      gMu->SetMarkerStyle(22);gMu->SetMarkerColor(4); gMu->SetMarkerSize(1.4);
      gMu->SetLineStyle(1); gMu->SetLineColor(4); gMu->SetLineWidth(1.4);
      gMu->Draw("pSAME");
      gMu2->SetMarkerStyle(22);gMu2->SetMarkerColor(4); gMu2->SetMarkerSize(.1);
      gMu2->SetLineStyle(1); gMu2->SetLineColor(4); gMu2->SetLineWidth(1.4);
      gMu2->SetFillColor(0); gMu2->SetFillStyle(0);
      gMu2->Draw("spe2");
      leg1->AddEntry(gMu,"#mu^{+}#mu^{-}, 2.5<y<4.0", "p");
//prepare graph with all err.
      gMu3->SetMarkerStyle(22);gMu3->SetMarkerColor(4); gMu3->SetMarkerSize(1.4);
      gMu3->SetLineStyle(1); gMu3->SetLineColor(4); gMu3->SetLineWidth(1.4);

      gMu1->SetMarkerStyle(22);gMu1->SetMarkerColor(4); gMu1->SetMarkerSize(1.2);
      gMu1->SetLineStyle(1); gMu1->SetLineColor(4); gMu1->SetLineWidth(1.4);
      gMu4->SetMarkerStyle(22);gMu4->SetMarkerColor(4); gMu4->SetMarkerSize(.1);
      gMu4->SetLineStyle(1); gMu4->SetLineColor(4); gMu4->SetLineWidth(1.4);
      gMu4->SetFillColor(0); gMu4->SetFillStyle(0);
  }
  else{
      fprintf(stderr, "Cannot open %s\n", filename);
  }

  leg1->Draw();

  TLatex *lat=new TLatex();  lat->SetNDC(kTRUE);
  lat->SetTextColor(1); lat->SetTextFont(42); lat->SetTextSize(.042);
  lat->DrawLatex(0.56, 0.9, Form("ALICE  pp   #sqrt{s}=7 TeV"));
  lat->SetTextSize(.04);
  lat->DrawLatex(0.2, 0.18, "#pm7% scale uncertainty (luminosity)");
  //lat->DrawLatex(0.69, 0.77, "(luminosity)");

  c4->SaveAs(Form("jpsi_dsdpt_SPD%1d.eps",SPD));
  
// ...and dsigma/dy
  Float_t sj2e[1]={5.95};
  Float_t sej2e[1]={0.65}; //stat err.
  Float_t se2j2e[1]={0.815}; //13.7% syst. err. without Lumi; 0.94 total syst. err.
  Float_t se3j2e[1]; se3j2e[0]=TMath::Sqrt(sej2e[0]*sej2e[0]+se2j2e[0]*se2j2e[0]+errL*errL*sj2e[0]*sj2e[0]); //all err.
  Float_t se4j2e[1]; se4j2e[0]=TMath::Sqrt(sej2e[0]*sej2e[0]+se2j2e[0]*se2j2e[0]); //all err.-lumi
  Float_t se5j2e[1]; se5j2e[0]=errL*sj2e[0]; //only lumi

  Float_t yj2e[1]={0.};
  Float_t yej2e[1]={0.9}; 
  Float_t yej2e2[1]={0.2}; //for box
  
  TGraph *gj2e = new TGraphErrors(1,yj2e,sj2e,yej2e,sej2e); 
  TGraph *gj2e2 = new TGraphErrors(1,yj2e,sj2e,yej2e2,se2j2e); //syst. err.
  TGraph *gj2e3 = new TGraphErrors(1,yj2e,sj2e,yej2e,se3j2e); //all err.
  TGraph *gj2e4 = new TGraphErrors(1,yj2e,sj2e,yej2e,se4j2e); //all err.-lumi
  TGraph *gj2e5 = new TGraphErrors(1,yj2e,sj2e,yej2e2,se5j2e); //lumi

  TCanvas *c1 = new TCanvas("dsigdy","dsigdy",10,10,610,610);
  c1->SetTopMargin(0.04);    c1->SetBottomMargin(.1);
  c1->SetRightMargin(0.04);    c1->SetLeftMargin(0.11);
  c1->SetTicky(); //set ticks on right border (2) adds values too...
  pad = c1->cd(1);
  pad->SetLeftMargin(0.16); pad->SetRightMargin(0.04);
  pad->SetTopMargin(0.04); pad->SetBottomMargin(0.16);
  pad->SetFrameLineWidth(0.1); //histo frame ...bad in eps...

  TH2D *ho = new TH2D("ho", "ho", 10, -5, 5, 10, 0, 10.); //...just for frame
  ho->GetXaxis()->SetTitle("y");
  ho->GetYaxis()->SetTitle("d#sigma_{J/#psi} /dy (#mub)");
  ho->GetXaxis()->SetTitleOffset(1.);
  ho->GetYaxis()->SetTitleOffset(1.1);
  ho->GetYaxis()->SetLabelOffset(.01);
  ho->SetTitleSize(0.055,"XY");ho->SetTitleFont(42,"XY");
  ho->SetLabelSize(0.043,"XY");ho->SetLabelFont(42,"XY");
  ho->SetTitle("");
  ho->Draw(""); 
  
  gj2e->SetMarkerStyle(20);gj2e->SetMarkerColor(2); gj2e->SetMarkerSize(1.5);
  gj2e->SetLineStyle(1); gj2e->SetLineColor(2); gj2e->SetLineWidth(1);
  gj2e->Draw("psame");

  gj2e2->SetMarkerStyle(0);gj2e2->SetMarkerColor(2); gj2e2->SetMarkerSize(0.1);
  gj2e2->SetLineStyle(1); gj2e2->SetLineColor(2); //gj2e2->SetLineWidth(5);
  gj2e2->SetFillColor(0); gj2e2->SetFillStyle(0);
  gj2e2->Draw("spe2");

  //gj2e2->Draw("pSAME");

  leg2 = new TLegend(.2,0.78,.56,0.93);
  leg2->SetTextFont(42);    leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);    leg2->SetMargin(0.15); //separation symbol-text
  leg2->SetEntrySeparation(0.14);   
  //leg2->SetHeader("ALICE data"); // preliminary");
  leg2->AddEntry(gj2e,"e^{+}e^{-}","p");
  
//ALICE muons  
  float yMu[5],yrMu[5],yeMu[5],ye2Mu[5],sMu[5],seMu[5],se1Mu[5],se2Mu[5],sesMu[5],ses1Mu[5],ses2Mu[5],ses3Mu[5]; //
  char filename[200];
  //sprintf(filename, "Mu/jpsi2m_dsdy.dat");
  sprintf(filename, "dsigdy_y2.5_4.0.txt");
  ifstream mfile(filename);
  if(loc = gSystem->FindFile(".",filename)){
      int i0=0;
      while ( mfile >> yMu[i0] >> yeMu[i0] >> sMu[i0]  >> seMu[i0]  >> se1Mu[i0] >> se2Mu[i0] ) { 
	  se2Mu[i0]=TMath::Sqrt(se1Mu[i0]*se1Mu[i0]+se2Mu[i0]*se2Mu[i0]); //add correl. and uncorr. syst err
	  sesMu[i0]=TMath::Sqrt(se2Mu[i0]*se2Mu[i0]-errL*errL*sMu[i0]*sMu[i0]); //add correl. and uncorrel. syst err. and subtract lumi
	  printf("%4.1f %4.1f %4.1f  %4.1f %4.1f %4.1f \n", yMu[i0],sMu[i0],seMu[i0],se1Mu[i0],se2Mu[i0],sesMu[i0]);
	  ses2Mu[i0]=TMath::Sqrt(se2Mu[i0]*se2Mu[i0]+seMu[i0]*seMu[i0]); //add correl. and uncorrel. syst err. incl. lumi
	  ses1Mu[i0]=TMath::Sqrt(ses2Mu[i0]*ses2Mu[i0]-errL*errL*sMu[i0]*sMu[i0]); //all err. - 7% lumi
	  ses3Mu[i0]=errL*sMu[i0]; // lumi
	  yrMu[i0]=-yMu[i0];
	  i0++; 
      }
      cout <<i0 << " points from file " << filename << endl;
      mfile.close();
      const int Nmu=i0;
      TGraphErrors *gj2m = new TGraphErrors(Nmu,yrMu,sMu,yeMu,seMu); //stat. 
      TGraphErrors *gj2mr = new TGraphErrors(Nmu,yMu,sMu,yeMu,seMu); 
      TGraphErrors *gj2m2 = new TGraphErrors(Nmu,yrMu,sMu,yeMu,sesMu); //syst-7%
      TGraphErrors *gj2m2r = new TGraphErrors(Nmu,yMu,sMu,yeMu,sesMu);
      TGraphErrors *gj2m3 = new TGraphErrors(Nmu,yrMu,sMu,yeMu,ses2Mu); //all err.  
      TGraphErrors *gj2m3r = new TGraphErrors(Nmu,yMu,sMu,yeMu,ses2Mu); 
      TGraphErrors *gj2m4 = new TGraphErrors(Nmu,yrMu,sMu,yeMu,ses1Mu); //all err. -lumi 
      TGraphErrors *gj2m4r = new TGraphErrors(Nmu,yMu,sMu,yeMu,ses1Mu); 
      TGraphErrors *gj2m5 = new TGraphErrors(Nmu,yrMu,sMu,yeMu,ses3Mu); //lumi alone
      TGraphErrors *gj2m5r = new TGraphErrors(Nmu,yMu,sMu,yeMu,ses3Mu); 
      gj2m2->SetMarkerStyle(0);gj2m2->SetMarkerColor(3); gj2m2->SetMarkerSize(0.1);
      gj2m2->SetLineStyle(1); gj2m2->SetLineColor(4); gj2e2->SetLineWidth(1.2);
      gj2m2->SetFillColor(0); gj2m2->SetFillStyle(0);
      gj2m2->Draw("spe2");
      //gj2m2->Draw("pSAME");
      
      gj2m->SetMarkerStyle(22);gj2m->SetMarkerColor(4); gj2m->SetMarkerSize(1.4);
      gj2m->SetLineStyle(1); gj2m->SetLineColor(4); gj2m->SetLineWidth(1.2);
      gj2m->Draw("pSAME");
  
      gj2mr->SetMarkerStyle(26);gj2mr->SetMarkerColor(4); gj2mr->SetMarkerSize(1.4);
      gj2mr->SetLineStyle(1); gj2mr->SetLineColor(4); gj2mr->SetLineWidth(1.2);
      gj2mr->Draw("pSAME");
      
      gj2m2r->SetMarkerStyle(0);gj2m2r->SetMarkerColor(3); gj2m2r->SetMarkerSize(0.1);
      gj2m2r->SetLineStyle(1); gj2m2r->SetLineColor(4); gj2e2->SetLineWidth(1.2);
      gj2m2r->SetFillColor(0); gj2m2r->SetFillStyle(0);
      gj2m2r->Draw("spe2");
      
      //leg2->AddEntry(gj2m,"#mu^{+}#mu^{-}, 2.5<y<4.0","p"); 
      leg2->AddEntry(gj2m,"#mu^{+}#mu^{-}","p"); 
      leg2->AddEntry(gj2mr,"#mu^{+}#mu^{-}, reflected","p");
      
  }
  else{
	fprintf(stderr, "Cannot open muon file %s\n", filename);
  }
  
  leg2->Draw();
  
  lat->DrawLatex(0.30, 0.2, "#pm7% scale uncertainty (luminosity)");

  lat->SetTextSize(.042);
  lat->DrawLatex(0.56, 0.88, Form("ALICE  pp   #sqrt{s}=7 TeV"));

  c1->SaveAs(Form("jpsi_dsdy_SPD%1d.eps",SPD));

// ***************************************************
// ...and now with CMS http://arxiv.org/abs/1011.4193
// ***************************************************

  Double_t pMax=12.; //pt max for plot
  
  TCanvas *c4b=new TCanvas("dsigdpt 1","dsigdpt 1",25,25,625,625);
  c4b->SetTopMargin(0.03);  c4b->SetRightMargin(0.03);
  c4b->SetLeftMargin(0.158);  c4b->SetBottomMargin(0.145);
  pad = c4b->cd(1); pad->SetLogy();
  
  TH2D *ho2 = new TH2D("ho2", "ho2", 12, 0, pMax, 10, 0.003, 3.); //...just for frame
  ho2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ho2->GetYaxis()->SetTitle("d^{2}#sigma_{J/#psi} /dp_{T}dy (#mub/GeV/c)");
  ho2->GetXaxis()->SetTitleOffset(1.1);
  ho2->GetYaxis()->SetTitleOffset(1.3);
  ho2->GetYaxis()->SetLabelOffset(.01);
  ho2->SetTitleSize(0.055,"XY");ho2->SetTitleFont(42,"XY");
  ho2->SetLabelSize(0.043,"XY");ho2->SetLabelFont(42,"XY");
  ho2->SetTitle("");
  ho2->Draw("");
  //ho1->Draw("");

  gSpect3->SetMarkerStyle(20); gSpect3->SetMarkerSize(1.3);
  gSpect3->SetMarkerColor(kMarkCol); gSpect3->SetLineColor(kMarkCol);
  gSpect3->SetLineWidth(1.4);
  //gSpect3->Draw("psame"); //all err.

  gSpect1->SetMarkerStyle(20); gSpect1->SetMarkerSize(1.1);
  gSpect1->SetMarkerColor(kMarkCol); gSpect1->SetLineColor(kMarkCol);
  gSpect1->SetLineWidth(1.4);
  gSpect1->Draw("psame");
  gSpect4->SetMarkerStyle(0);gSpect4->SetMarkerColor(2); gSpect4->SetMarkerSize(0.1);
  gSpect4->SetLineStyle(1); gSpect4->SetLineColor(2); gSpect4->SetLineWidth(1.4);
  gSpect4->SetFillColor(0); gSpect4->SetFillStyle(0);
  gSpect4->Draw("spe2");

  //gMu3->Draw("pSAME");
  gMu1->SetMarkerStyle(22); gMu1->SetMarkerSize(1.1);
  gMu1->Draw("pSAME");
  gMu4->Draw("spe2");

  //leg1b = new TLegend(.25,0.18,.61,0.38);
  leg1b = new TLegend(.2,0.18,.61,0.42); //with Atlas
  leg1b->SetTextFont(kFont);    leg1b->SetBorderSize(0);
  leg1b->SetFillStyle(0);    leg1b->SetMargin(0.16); //separation symbol-text
  leg1b->SetEntrySeparation(0.2);    
  leg1b->AddEntry(gSpect1, "ALICE e^{+}e^{-}, |y|<0.9", "p"); 
  leg1b->AddEntry(gMu1,"ALICE  #mu^{+}#mu^{-}, 2.5<y<4.0", "p");

  float ptCMS1[10],pteCMS1[10],pte2CMS1[10],sCMS1[10],seCMS1[10],se1CMS1[10],se2CMS1[10],se3CMS1[10];
  float ptCMS[90],pteCMS[90],pte2CMS[90],sCMS[90],seCMS[90],se1CMS[90],se2CMS[90],se3CMS[90];
  Float_t ScaleCMS=.001/BR; // nb->ub
  Float_t sigCMS=0.,sigeCMS=0.,sige2CMS=0.;
  
  char filename[200];
  sprintf(filename, "dsigdpt-cms_y1.2.txt"); 
  ifstream vfile(filename);
  if(loc = gSystem->FindFile(".",filename)){
      int i0=0;
      while ( vfile >> ptCMS1[i0] >> pteCMS1[i0] >> sCMS1[i0] >> seCMS1[i0] >> se2CMS1[i0] ) { 
	  sCMS1[i0]*=ScaleCMS;
	  seCMS1[i0]*=ScaleCMS;
	  se2CMS1[i0]*=ScaleCMS;
	  //se2CMS1[i0]=TMath::Sqrt(seCMS1[i0]*seCMS1[i0]+se2CMS1[i0]*se2CMS1[i0]); //no lumi?
	  se2CMS1[i0]=TMath::Sqrt(seCMS1[i0]*seCMS1[i0]+se2CMS1[i0]*se2CMS1[i0]-0.11*0.11*sCMS1[i0]*sCMS1[i0]); //take out 11% lumi
	  se1CMS1[i0]=TMath::Sqrt(se2CMS1[i0]*se2CMS1[i0]+0.11*0.11*sCMS1[i0]*sCMS1[i0]);
	  se3CMS1[i0]=0.11*sCMS1[i0]; //11% lumi err.
	  pte2CMS1[i0]=0.12;
	  i0++; 
      }
      cout <<i0 << " points from file " << filename << endl;
      vfile.close();
      const int Ncms1=i0;
      TGraphErrors *gCMS1 = new TGraphErrors(Ncms1,ptCMS1,sCMS1,pteCMS1,se2CMS1); 
      gCMS1->SetMarkerStyle(24);gCMS1->SetMarkerColor(6); gCMS1->SetMarkerSize(1.1);
      gCMS1->SetLineStyle(1); gCMS1->SetLineColor(6); gCMS1->SetLineWidth(1);
      gCMS1->Draw("pSAME");
      TGraphErrors *gCMS1b = new TGraphErrors(Ncms1,ptCMS1,sCMS1,pte2CMS1,se3CMS1); 
      gCMS1b->SetMarkerStyle(24);gCMS1b->SetMarkerColor(6); gCMS1b->SetMarkerSize(0.1);
      gCMS1b->SetLineStyle(1); gCMS1b->SetLineColor(6); gCMS1b->SetLineWidth(1);
      gCMS1b->SetFillColor(0); gCMS1b->SetFillStyle(0);
      gCMS1b->Draw("spe2");

      TGraphErrors *gCMS1c = new TGraphErrors(Ncms1,ptCMS1,sCMS1,pteCMS1,se1CMS1); 
      gCMS1c->SetMarkerStyle(27);gCMS1c->SetMarkerColor(6); gCMS1c->SetMarkerSize(1.2);
      gCMS1c->SetLineStyle(1); gCMS1c->SetLineColor(6); gCMS1c->SetLineWidth(1);

      leg1b->AddEntry(gCMS1,"CMS, |y|<1.2", "p");
  }
  else{
      fprintf(stderr, "!!! Cannot open %s\n", filename);
  }

  //sprintf(filename, "CMS/dsigdpt_y1.6-2.4.dat");
  sprintf(filename, "dsigdpt-cms_y1.6-2.4.txt");
  ifstream vfile(filename);
  if(loc = gSystem->FindFile(".",filename)){
      int i0=0;
      while ( vfile >> ptCMS[i0] >> pteCMS[i0] >> sCMS[i0]>> seCMS[i0]>> se2CMS[i0] ) {
	  //printf("%5.1f %5.1f %7.2f \n",ptCMS[i0],pteCMS[i0],sCMS[i0]);
	  sCMS[i0]*=ScaleCMS;
	  seCMS[i0]*=ScaleCMS;
	  se2CMS[i0]*=ScaleCMS;
	  //se2CMS[i0]=TMath::Sqrt(seCMS[i0]*seCMS[i0]+se2CMS[i0]*se2CMS[i0]); //-no lumi ?
	  se2CMS[i0]=TMath::Sqrt(seCMS[i0]*seCMS[i0]+se2CMS[i0]*se2CMS[i0]-0.11*0.11*sCMS[i0]*sCMS[i0]); //take lumi out
	  se1CMS[i0]=TMath::Sqrt(se2CMS[i0]*se2CMS[i0]+0.11*0.11*sCMS[i0]*sCMS[i0]); //all err (+11% lumi)
	  se3CMS[i0]=0.11*sCMS[i0]; // 11% lumi
	  ptCMS[i0]=(pteCMS[i0]+ptCMS[i0])/2.;
	  pteCMS[i0]=(pteCMS[i0]-ptCMS[i0]);
	  pte2CMS[i0]=0.1*pMax;
	  sigCMS+=sCMS[i0]*pteCMS[i0]*2.;
	  //sigeCMS+=(seCMS[i0]*pteCMS[i0]*2.)*(seCMS[i0]*pteCMS[i0]*2.);
	  //sige2CMS+=(se2CMS[i0]*pteCMS[i0]*2.)*(se2CMS[i0]*pteCMS[i0]*2.);
	  sigeCMS+=(se2CMS[i0]*pteCMS[i0]*2.)*(se2CMS[i0]*pteCMS[i0]*2.); //stat. + syst.
	  sige2CMS+=(se3CMS[i0]*pteCMS[i0]*2.)*(se3CMS[i0]*pteCMS[i0]*2.);
	  i0++;
      }
      cout <<i0 << " points from file " << filename << endl;
      vfile.close();
      sigeCMS=TMath::Sqrt(sigeCMS);
      sige2CMS=TMath::Sqrt(sige2CMS);
      printf("CMS(): %5.2f +- %5.2f +- %5.2f \n",sigCMS,sigeCMS,sige2CMS);
      
      const int Ncms=i0;
      //TGraphErrors *gCMS = new TGraphErrors(Ncms,ptCMS,sCMS,pteCMS,seCMS); 
      TGraphErrors *gCMS = new TGraphErrors(Ncms,ptCMS,sCMS,pteCMS,se2CMS); 
      gCMS->SetMarkerStyle(26);gCMS->SetMarkerColor(6); gCMS->SetMarkerSize(1.3);
      gCMS->SetLineStyle(1); gCMS->SetLineColor(6); gCMS->SetLineWidth(1);
      //gCMS->Draw("pSAME");
      TGraphErrors *gCMS2 = new TGraphErrors(Ncms,ptCMS,sCMS,pte2CMS,se3CMS); 
      gCMS2->SetMarkerStyle(26);gCMS2->SetMarkerColor(6); gCMS2->SetMarkerSize(1.3);
      gCMS2->SetLineStyle(1); gCMS2->SetLineColor(6); gCMS2->SetLineWidth(1);
      gCMS2->SetFillColor(0); gCMS2->SetFillStyle(0);
      //gCMS2->Draw("spe2");
      //leg1b->AddEntry(gCMS,"CMS, 1.6<|y|<2.4", "p");
  }
  else{
      fprintf(stderr, "Cannot open %s\n", filename);
  }

// ATLAS http://arxiv.org/abs/1104.3038
  float ptA[90],pteA[90],pte2A[90],sA[90],seA[90],se1A[90],se2A[90],se3A[90];
  Float_t ScaleA=.000001/BR; // pb->ub
  
  sprintf(filename, "dsigdpt_atlas_y0.75.txt"); 
  ifstream vfile(filename);
  if(loc = gSystem->FindFile(".",filename)){
      int i0=0;
      while ( vfile >> ptA[i0] >> pteA[i0] >> sA[i0] >> seA[i0] >> se2A[i0] ){
	  ptA[i0]=(pteA[i0]+ptA[i0])/2.; //center of bin
	  pteA[i0]=(pteA[i0]-ptA[i0]);
	  sA[i0]*=ScaleA;
	  seA[i0]*=ScaleA;
	  se2A[i0]*=ScaleA;
	  //se2CMS1[i0]=TMath::Sqrt(seCMS1[i0]*seCMS1[i0]+se2CMS1[i0]*se2CMS1[i0]); //no lumi?
	  se2CMS1[i0]=TMath::Sqrt(seA[i0]*seA[i0]+se2A[i0]*se2A[i0]-0.034*0.034*sA[i0]*sA[i0]); //take out 3.4% lumi
	  se1A[i0]=TMath::Sqrt(se2A[i0]*se2A[i0]+0.034*0.034*sA[i0]*sA[i0]);
	  se3A[i0]=0.034*sA[i0]; //3.4% lumi err.
	  pte2A[i0]=0.12;
	  i0++; 
      }
      cout <<i0 << " points from file " << filename << endl;
      vfile.close();

      const int Natl=i0;
      TGraphErrors *gATL = new TGraphErrors(Natl,ptA,sA,pteA,se2A); 
      gATL->SetMarkerStyle(27);gATL->SetMarkerColor(1);gATL->SetMarkerSize(1.1);
      gATL->SetLineStyle(1); gATL->SetLineColor(1); gATL->SetLineWidth(1);
      gATL->Draw("pSAME");
      TGraphErrors *gATLb = new TGraphErrors(Natl,ptA,sA,pte2A,se3A); 
      gATLb->SetMarkerStyle(24);gATLb->SetMarkerColor(1); gATLb->SetMarkerSize(0.1);
      gATLb->SetLineStyle(1); gATLb->SetLineColor(1); gATLb->SetLineWidth(1);
      gATLb->SetFillColor(0); gATLb->SetFillStyle(0);
      gATLb->Draw("spe2");

      TGraphErrors *gATLc = new TGraphErrors(Natl,ptA,sA,pteA,se1A); 
      gATLc->SetMarkerStyle(27);gATLc->SetMarkerColor(1); gATLc->SetMarkerSize(1.2);
      gATLc->SetLineStyle(1); gATLc->SetLineColor(1); gATLc->SetLineWidth(1);

      leg1b->AddEntry(gATL,"ATLAS, |y|<0.75", "p");
  }
  else{
      fprintf(stderr, "!!! Cannot open %s\n", filename);
  }


//LHCb pt, http://arxiv.org/abs/1103.0423 ...output of LHCb_Jpsi_pp7TeV.C

  float ptL[90],pteL[90],pte2L[90],sLpt[90],seLpt[90],se1Lpt[90],se2Lpt[90],se3Lpt[90];
  float ScaleLHCbPt=1./1.5; //2.5-4.0 -> dsig/dptdy

  sprintf(filename, "dsigdpt_lhcb_y2.5-4.0.txt");
  ifstream vfile(filename);
  if(loc = gSystem->FindFile(".",filename)){
      int i0=0;
      while ( vfile >> ptL[i0] >> pteL[i0] >> sLpt[i0]>> seLpt[i0]>> se1Lpt[i0] ) {
	  //printf("%5.1f %5.1f %7.2f \n",ptL[i0],pteL[i0],sLpt[i0]);
	  se1Lpt[i0]=TMath::Sqrt(seLpt[i0]*seLpt[i0]+se1Lpt[i0]*se1Lpt[i0]); //all err.
	  se2Lpt[i0]=TMath::Sqrt(se1Lpt[i0]*se1Lpt[i0]-0.10*.1*sLpt[i0]*sLpt[i0]); //-lumi (10%)
	  sLpt[i0]*=ScaleLHCbPt;
	  se1Lpt[i0]*=ScaleLHCbPt;
	  se2Lpt[i0]*=ScaleLHCbPt;
	  se3Lpt[i0]=0.10*sLpt[i0]; // 10% lumi
	  pte2L[i0]=0.01*pMax;
	  i0++;
      }
      cout <<i0 << " points from file " << filename << endl;
      vfile.close();
      //sigeL=TMath::Sqrt(sigeL);
      //printf("LHCb: %5.2f +- %5.2f +- %5.2f \n",sigCMS,sigeCMS,sige2CMS);
      
      const int Nlhcb=i0;
      TGraphErrors *gLHCBpt1 = new TGraphErrors(Nlhcb,ptL,sLpt,pteL,se1Lpt);//all err.
      gLHCBpt1->SetMarkerStyle(26);gLHCBpt1->SetMarkerColor(3);gLHCBpt1->SetMarkerSize(1.);
      gLHCBpt1->SetLineStyle(1); gLHCBpt1->SetLineColor(3); gLHCBpt1->SetLineWidth(1);

      TGraphErrors *gLHCBpt = new TGraphErrors(Nlhcb,ptL,sLpt,pteL,se2Lpt); 
      gLHCBpt->SetMarkerStyle(26);gLHCBpt->SetMarkerColor(3);gLHCBpt->SetMarkerSize(1.);
      gLHCBpt->SetLineStyle(1); gLHCBpt->SetLineColor(3); gLHCBpt->SetLineWidth(1);
      gLHCBpt->Draw("pSAME");
      TGraphErrors *gLHCBpt2 = new TGraphErrors(Nlhcb,ptL,sLpt,pte2L,se3Lpt); 
      gLHCBpt2->SetMarkerStyle(1);gLHCBpt2->SetMarkerColor(3);gLHCBpt2->SetMarkerSize(.1);
      gLHCBpt2->SetLineStyle(1); gLHCBpt2->SetLineColor(3); gLHCBpt2->SetLineWidth(1);
      gLHCBpt2->SetFillColor(0); gLHCBpt2->SetFillStyle(0);
      gLHCBpt2->Draw("spe2");
      leg1b->AddEntry(gLHCBpt,"LHCb, 2.5<|y|<4.0", "p");
  }
  else{
      fprintf(stderr, "Cannot open %s\n", filename);
  }

  leg1b->Draw();

  TLatex *lat=new TLatex();  lat->SetNDC(kTRUE);
  lat->SetTextColor(1); lat->SetTextFont(42); lat->SetTextSize(.042);
  lat->DrawLatex(0.65, 0.9, Form("pp   #sqrt{s}=7 TeV"));
  lat->SetTextSize(.04);
  //lat->DrawLatex(0.2, 0.18, "#pm7% scale uncertainty (luminosity)");
  //lat->DrawLatex(0.69, 0.77, "(luminosity)");

  c4b->SaveAs(Form("jpsi_dsdpt1a_SPD%1d.eps",SPD));

// NOW ALL ERR.
  TCanvas *c4c=new TCanvas("dsigdpt 2","dsigdpt 2",35,35,635,635);
  c4c->SetTopMargin(0.03);  c4c->SetRightMargin(0.03);
  c4c->SetLeftMargin(0.15);  c4c->SetBottomMargin(0.145);
  pad = c4c->cd(1); pad->SetLogy();

  ho2->Draw(""); //->8 or 12 GeV
  //ho1->Draw(""); //->8 GeV
  gSpect3->Draw("psame"); //all err.
  gMu3->Draw("pSAME"); //all err
  gCMS1c->Draw("pSAME");
  gATLc->Draw("pSAME");
  gLHCBpt1->Draw("pSAME");
  leg1b->Draw();
  lat->DrawLatex(0.65, 0.9, Form("pp   #sqrt{s}=7 TeV"));

  c4c->SaveAs(Form("jpsi_dsdpt2a_SPD%1d.eps",SPD));

////////////////////////////////////////////////////////////////////////

  TCanvas *c1b = new TCanvas("dsigdy 1","dsigdy 1",15,15,615,615);
  c1b->SetTopMargin(0.04);    c1b->SetBottomMargin(.1);
  c1b->SetRightMargin(0.04);    c1b->SetLeftMargin(0.11);
  c1b->SetTicky(); //set ticks on right border (2) adds values too...
  pad = c1b->cd(1);
  pad->SetLeftMargin(0.16); pad->SetRightMargin(0.04);
  pad->SetTopMargin(0.04); pad->SetBottomMargin(0.16);
  pad->SetFrameLineWidth(0.1); //histo frame ...bad in eps...

  ho->Draw("");
  
  //gj2e->Draw("psame");
  //gj2e2->Draw("spe2");
  //gj2m->Draw("pSAME"); gj2mr->Draw("pSAME");  
  //gj2m2->Draw("spe2"); gj2m2r->Draw("spe2");
 
  gj2e4->SetMarkerStyle(20);gj2e4->SetMarkerColor(2); gj2e4->SetMarkerSize(1.5);
  gj2e4->SetLineStyle(1); gj2e4->SetLineColor(2); gj2e4->SetLineWidth(1);
  gj2e4->Draw("psame");

  gj2e5->SetMarkerStyle(0);gj2e5->SetMarkerColor(2); gj2e5->SetMarkerSize(0.1);
  gj2e5->SetLineStyle(1); gj2e5->SetLineColor(2); 
  gj2e5->SetFillColor(0); gj2e5->SetFillStyle(0);
  gj2e5->Draw("spe2");

  gj2m4->SetMarkerStyle(22);gj2m4->SetMarkerColor(4); gj2m4->SetMarkerSize(1.3);
  gj2m4->SetLineStyle(1); gj2m4->SetLineColor(4); gj2m4->SetLineWidth(1.4);
  gj2m4r->SetMarkerStyle(26);gj2m4r->SetMarkerColor(4); gj2m4r->SetMarkerSize(1.3);
  gj2m4r->SetLineStyle(1); gj2m4r->SetLineColor(4); gj2m4r->SetLineWidth(1.4);
  gj2m4->Draw("pSAME"); gj2m4r->Draw("pSAME");  
  gj2m5->SetMarkerSize(.01); gj2m5r->SetMarkerSize(.01); 
  gj2m5->SetFillColor(0); gj2m5->SetFillStyle(0); gj2m5->SetLineColor(4);
  gj2m5r->SetFillColor(0); gj2m5r->SetFillStyle(0); gj2m5r->SetLineColor(4);
  gj2m5->Draw("spe2"); gj2m5r->Draw("spe2");


  leg2b = new TLegend(.18,0.77,.6,0.94);
  leg2b->SetTextFont(42);    leg2b->SetBorderSize(0);
  leg2b->SetFillStyle(0);    leg2b->SetMargin(0.15); //separation symbol-text
  leg2b->SetEntrySeparation(0.14);   
  leg2b->AddEntry(gj2e4,"ALICE, e^{+}e^{-}","p");
  leg2b->AddEntry(gj2m4,"ALICE, #mu^{+}#mu^{-}","p"); 
  //leg2b->AddEntry(gj2mr,"ALICE, #mu^{+}#mu^{-}, reflected","p");

  Float_t sgCMS[1]; sgCMS[0]=sigCMS;
  Float_t sgeCMS[1]; sgeCMS[0]=sigeCMS; //stat.+syst. (-lumi)
  Float_t sge1CMS[1]; sge1CMS[0]=TMath::Sqrt(sigeCMS*sigeCMS+0.11*0.11*sigCMS*sigCMS); //+ 11% lumi
  Float_t sge2CMS[1]; sge2CMS[0]=0.11*sigCMS; //lumi
  Float_t yCMS[1]={2.};  Float_t yCMSr[1]={-2.};
  Float_t yeCMS[1]={0.4}; 
  Float_t ye2CMS[1]={0.15}; //for box

  TGraphErrors *gCMSi = new TGraphErrors(1,yCMS,sgCMS,yeCMS,sgeCMS); 
  TGraphErrors *gCMSir = new TGraphErrors(1,yCMSr,sgCMS,yeCMS,sgeCMS); 
  TGraphErrors *gCMSi2 = new TGraphErrors(1,yCMS,sgCMS,ye2CMS,sge2CMS); 
  TGraphErrors *gCMSi2r = new TGraphErrors(1,yCMSr,sgCMS,ye2CMS,sge2CMS); 
  TGraphErrors *gCMSi1 = new TGraphErrors(1,yCMS,sgCMS,yeCMS,sge1CMS); 
  TGraphErrors *gCMSi1r = new TGraphErrors(1,yCMSr,sgCMS,yeCMS,sge1CMS); 

  gCMSi->SetMarkerStyle(21);gCMSi->SetMarkerColor(6); gCMSi->SetMarkerSize(1.4);
  gCMSi->SetLineStyle(1); gCMSi->SetLineColor(6); gCMSi->SetLineWidth(1.2);
  gCMSi->Draw("pSAME");//
  gCMSir->SetMarkerStyle(25);gCMSir->SetMarkerColor(6); gCMSir->SetMarkerSize(1.4);
  gCMSir->SetLineStyle(1); gCMSir->SetLineColor(6); gCMSir->SetLineWidth(1.2);
  gCMSir->Draw("pSAME");//reflected

  //TGraphErrors *gCMSi2 = new TGraphErrors(1,yCMS,sgCMS,ye2CMS,sge2CMS); 
  gCMSi2->SetMarkerStyle(21);gCMSi2->SetMarkerColor(6); gCMSi2->SetMarkerSize(1.3);
  gCMSi2->SetLineStyle(1); gCMSi2->SetLineColor(6); gCMSi2->SetLineWidth(1);
  gCMSi2->SetFillColor(0); gCMSi2->SetFillStyle(0);
  gCMSi2->Draw("spe2");
  gCMSi2r->SetMarkerStyle(25);gCMSi2r->SetMarkerColor(6); gCMSi2r->SetMarkerSize(.1);
  gCMSi2r->SetLineStyle(1); gCMSi2r->SetLineColor(6); gCMSi2r->SetLineWidth(1);
  gCMSi2r->SetFillColor(0); gCMSi2r->SetFillStyle(0);
  gCMSi2r->Draw("spe2");

  gCMSi1->SetMarkerStyle(21);gCMSi1->SetMarkerColor(6); gCMSi1->SetMarkerSize(1.4);
  gCMSi1->SetLineStyle(1); gCMSi1->SetLineColor(6); gCMSi1->SetLineWidth(1.2);
  gCMSi1r->SetMarkerStyle(25);gCMSi1r->SetMarkerColor(6); gCMSi1r->SetMarkerSize(1.4);
  gCMSi1r->SetLineStyle(1); gCMSi1r->SetLineColor(6); gCMSi1r->SetLineWidth(1.2);

  leg2b->AddEntry(gCMSi,"CMS", "p");
  

// LHCb http://arxiv.org/abs/1103.0423

  float yL[10],yrL[10],yeL[10],sL[10],se1L[10],se2L[10],se3L[10],sLb[10],sb1L[10],sb2L[10],sb3L[10];
  float se4L[10],se5L[10],ye2L[10];
  float ScaleLHCb=.001; // nb->ub

  sprintf(filename, "dsigdy_lhcb.txt");
  ifstream vfile(filename);
  if(loc = gSystem->FindFile(".",filename)) {
      int i0=0;
      while ( vfile >> yL[i0] >> yeL[i0] >> sL[i0] >> se1L[i0] >> se2L[i0]>> se3L[i0]>> sLb[i0] >> sb1L[i0] >> sb2L[i0]>> sb3L[i0] ) { 
	  //printf("%5.1f %7.2f \n",ptMu[i0],sPtMu[i0]);
	  yrL[i0]=-yL[i0];
	  se1L[i0]=TMath::Sqrt(se1L[i0]*se1L[i0]+se2L[i0]*se2L[i0]+se3L[i0]*se3L[i0]);
	  sb1L[i0]=TMath::Sqrt(sb1L[i0]*sb1L[i0]+sb2L[i0]*sb2L[i0]+sb3L[i0]*sb3L[i0]);
	  sL[i0]=sL[i0]+sLb[i0]; //prompt + non-prompt
	  se1L[i0]=TMath::Sqrt(se1L[i0]*se1L[i0]+sb1L[i0]*sb1L[i0]); //all err.
	  sL[i0]*=ScaleLHCb;
	  se1L[i0]*=ScaleLHCb;
	  se4L[i0]=TMath::Sqrt(se1L[i0]*se1L[i0]-0.1*0.1*sL[i0]*sL[i0]); // -10% lumi
	  se5L[i0]=0.1*sL[i0];
	  ye2L[i0]=0.13;
	  i0++;	  
      }
      cout <<i0 << " points from file " << filename << endl;
      vfile.close();
      const int NL=i0;
      TGraphErrors *gLHCb = new TGraphErrors(NL,yL,sL,yeL,se4L); 
      gLHCb->SetMarkerStyle(29);gLHCb->SetMarkerColor(3); gLHCb->SetMarkerSize(1.4);
      gLHCb->SetLineStyle(1); gLHCb->SetLineColor(3); gLHCb->SetLineWidth(1.4);
      gLHCb->Draw("pSAME");
      leg2b->AddEntry(gLHCb,"LHCb", "p");

      TGraphErrors *gLHCbr = new TGraphErrors(NL,yrL,sL,yeL,se4L); 
      gLHCbr->SetMarkerStyle(30);gLHCbr->SetMarkerColor(3); gLHCbr->SetMarkerSize(1.4);
      gLHCbr->SetLineStyle(1); gLHCbr->SetLineColor(3); gLHCbr->SetLineWidth(1.4);
      gLHCbr->Draw("pSAME");

      TGraphErrors *gLHCb2 = new TGraphErrors(NL,yL,sL,ye2L,se5L); 
      TGraphErrors *gLHCb2r = new TGraphErrors(NL,yrL,sL,ye2L,se5L); 

      TGraphErrors *gLHCb1 = new TGraphErrors(NL,yL,sL,ye2L,se1L);  //all err.
      TGraphErrors *gLHCb1r = new TGraphErrors(NL,yrL,sL,ye2L,se1L); 
      gLHCb1->SetMarkerStyle(29);gLHCb1->SetMarkerColor(3); gLHCb1->SetMarkerSize(1.4);
      gLHCb1->SetLineStyle(1); gLHCb1->SetLineColor(3); gLHCb1->SetLineWidth(1.4);
      gLHCb1r->SetMarkerStyle(30);gLHCb1r->SetMarkerColor(3); gLHCb1r->SetMarkerSize(1.4);
      gLHCb1r->SetLineStyle(1); gLHCb1r->SetLineColor(3); gLHCb1r->SetLineWidth(1.4);

      gLHCb2->SetMarkerStyle(29);gLHCb2->SetMarkerColor(3); gLHCb2->SetMarkerSize(.1);
// ps err if gLHCb2->SetMarkerSize(.01); !!!
      gLHCb2->SetLineStyle(1); gLHCb2->SetLineColor(3); gLHCb2->SetLineWidth(1.4);
      gLHCb2->SetFillColor(0); gLHCb2->SetFillStyle(0);
      gLHCb2->Draw("spe2");

      gLHCb2r->SetMarkerStyle(29);gLHCb2r->SetMarkerColor(3);gLHCb2r->SetMarkerSize(.1);
      gLHCb2r->SetLineStyle(1); gLHCb2r->SetLineColor(3); gLHCb2r->SetLineWidth(1.4);
      gLHCb2r->SetFillColor(0); gLHCb2r->SetFillStyle(0);
      gLHCb2r->Draw("spe2");
  }
  else{
      fprintf(stderr, "Cannot open %s\n", filename);
  }

  leg2b->Draw();
  
  //lat->DrawLatex(0.30, 0.2, "#pm7% scale uncertainty (luminosity)");
  lat->SetTextSize(.040); lat->DrawLatex(0.44, 0.2, "open: reflected");
  lat->SetTextSize(.042); lat->DrawLatex(0.68, 0.89, Form("pp   #sqrt{s}=7 TeV"));

  c1b->SaveAs(Form("jpsi_dsdy1_SPD%1d.eps",SPD));

  TCanvas *c1c = new TCanvas("dsigdy 2","dsigdy 2",35,35,635,635);
  c1c->SetTopMargin(0.04);    c1c->SetBottomMargin(.1);
  c1c->SetRightMargin(0.04);    c1c->SetLeftMargin(0.11);
  c1c->SetTicky(); //set ticks on right border (2) adds values too...
  pad = c1c->cd(1);
  pad->SetLeftMargin(0.16); pad->SetRightMargin(0.04);
  pad->SetTopMargin(0.04); pad->SetBottomMargin(0.16);
  pad->SetFrameLineWidth(0.1); //histo frame ...bad in eps...

  ho->Draw("");
  
  gj2e3->SetMarkerStyle(20);gj2e3->SetMarkerColor(2); gj2e3->SetMarkerSize(1.5);
  gj2e3->SetLineStyle(1); gj2e3->SetLineColor(2); gj2e3->SetLineWidth(1.2);
  gj2e3->Draw("psame"); //all err.

  gj2m3->SetMarkerStyle(22);gj2m3->SetMarkerColor(4); gj2m3->SetMarkerSize(1.4);
  gj2m3->SetLineStyle(1); gj2m3->SetLineColor(4); gj2m3->SetLineWidth(1.2);
  gj2m3r->SetMarkerStyle(26);gj2m3r->SetMarkerColor(4); gj2m3r->SetMarkerSize(1.4);
  gj2m3r->SetLineStyle(1); gj2m3r->SetLineColor(4); gj2m3r->SetLineWidth(1.2);
  gj2m3->Draw("pSAME"); gj2m3r->Draw("pSAME"); //all err.

  gCMSi1->Draw("pSAME"); gCMSi1r->Draw("pSAME");//reflected
  gLHCb1->Draw("pSAME"); gLHCb1r->Draw("pSAME"); //all err.  
  leg2b->Draw();
  
  lat->SetTextSize(.040); lat->DrawLatex(0.44, 0.2, "open: reflected");
  lat->SetTextSize(.042); lat->DrawLatex(0.68, 0.89, Form("pp   #sqrt{s}=7 TeV"));

  c1c->SaveAs(Form("jpsi_dsdy2_SPD%1d.eps",SPD));

//*** sqrt(s) dependence ***//

  Float_t Scale_ccbar=0.0075; 
    const Int_t Ns=8;
    Float_t e[Ns]={0.2,0.9,2.75,3.94,5.5,7,10,14};
    //Float_t s_c[Ns],c[Ns],cy0[Ns]={0.112,0.300,0.683,0.822,0.940,1.147,1.386};
    Float_t s_c[Ns],c[Ns],cy0[Ns]={0.100,0.274,0.515,0.629,0.757,0.866,1.058,1.28};
    Float_t s_b[Ns],b[Ns],by0[Ns]={0.001,0.006,0.027,0.035,0.043,0.056,0.072};
    Float_t csca=1000*Scale_ccbar, bsca=1000*0.0116; //mb->ub * fraction to J/psi
    for (int i=0; i<Ns ;i++) { //there must be a clever way to do this...
      cy0[i]*=8.; //csca; //*0.83; //?
      by0[i]*=bsca;
    }

    TGraph *f1 = new TGraph(Ns,e,cy0); float f1M=20, f1m=0.4;
    TGraph *f2 = new TGraph(Ns,e,by0); 

    TCanvas *b1 = new TCanvas("b1","cross section",10,10,510,510);
    b1->SetTopMargin(0.03);    b1->SetBottomMargin(.12);
    b1->SetRightMargin(0.03);    b1->SetLeftMargin(0.12);
    b1->SetTicky(); //set ticks on right border (2) adds values too...
    //b1->Divide(1,2,0.0,0.0);
    b1->SetLogx(); b1->SetLogy();

    f1->SetMarkerStyle(25); f1->SetMarkerColor(1); f1->SetMarkerSize(1.);
    f1->SetLineStyle(1); f1->SetLineColor(1); f1->SetLineWidth(1.2);
    f1->SetMaximum(f1M);    f1->SetMinimum(f1m);
    f1->GetXaxis()->SetTitle("#sqrt{s} (TeV)");
    f1->GetYaxis()->SetTitle("d#sigma_{J/#psi} /dy (#mub)");
    f1->GetXaxis()->SetTitleOffset(1.1);
    f1->GetYaxis()->SetTitleOffset(1.1);
    f1->GetHistogram()->SetTitleSize(0.05,"XY");f1->GetHistogram()->SetTitleFont(42,"XY");
    f1->GetHistogram()->SetLabelSize(0.042,"XY");f1->GetHistogram()->SetLabelFont(42,"XY");
    f1->SetTitle("");
    f1->GetHistogram()->SetAxisRange(0.1,10000,"X"); //no effect!
    //f1->GetHistogram()->GetXaxis()->SetRange(0.1,10000.);
    f1->Draw("AL");

    f2->SetMarkerStyle(24); f2->SetMarkerColor(1); f2->SetMarkerSize(1.);
    f2->SetLineStyle(2); f2->SetLineColor(4); f2->SetLineWidth(1.3);
    f2->SetTitle("");
    //f2->Draw("lSAME");
    leg1 = new TLegend(.3,0.15,.9,0.3);
    leg1->SetHeader("line: NLO (MNR), m_{c}=1.2 GeV,  #mu_{F} = #mu_{R} = 2m_{c}");
    //leg1->AddEntry(f1,Form("d#sigma_{c#bar{c}} /dy*%6.4f",Scale_ccbar),"l");
    leg1->AddEntry(f1,"d#sigma_{c#bar{c}} /dy, scaled to CDF data point","l");
    leg1->SetTextFont(42);    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);    leg1->SetMargin(0.35); //separation symbol-text
    leg1->SetEntrySeparation(0.3);   
    //leg1->Draw();

    TLatex *lat3=new TLatex();  lat3->SetNDC(kTRUE);
    lat3->SetTextColor(1);lat3->SetTextFont(42);lat3->SetTextSize(.03);
    lat3->DrawLatex(0.38, 0.25,"line: d#sigma_{c#bar{c}} /dy, scaled to CDF data point");
    lat3->DrawLatex(0.38, 0.2,"NLO (MNR), m_{c}=1.2 GeV,  #mu_{F} = #mu_{R} = 2m_{c}");

    float sphe[1]={0.2}; //PHENIX
    float brjphe[1]={0.0443}; //http://arxiv.org/abs/hep-ex/0611020 44.3 +-1.4 +-5.10 +-4.50 (BRxdsigma/dy)
    float brjephe[1]={0.0069}; //3 errors in quadrature
    float jphe[1]={0}; float jephe[1]={0};
    jphe[0]=brjphe[0]/BR;
    jephe[0]=brjephe[0]/BR;

    float scdf[1]={1.96}; //CDF
    float jcdf[1]={3.4}; //4.08/1.2 (y) 4.08 \pm 0.02 (stat)^{+0.36}_{-0.33} (syst) \mu {\rm b}$
    float jecdf[1]={0.3};//0.36/1.2
 
//ALICE electrons
    //printf(" *** Total Syst. Err.: %4.2f (%4.1f%)\n",ae3_yall[0],100.*ae3_yall[0]/a_yall[0]);
    //printf("dsigma/dy(|y|<0.9): %6.2f+-%4.2f+-%4.2f\n",a_yall[0],ae_yall[0],ae3_yall[0]);

    Float_t sali[1]={7.}; 
    //Float_t yej2e2[1]={.5}; //for box

    TGraph *gphe = new TGraphErrors(1,sphe,jphe,0,jephe); 
    gphe->SetMarkerStyle(22);gphe->SetMarkerColor(4); gphe->SetMarkerSize(1.7);
    gphe->SetLineStyle(1); gphe->SetLineColor(4); gphe->SetLineWidth(1);
    gphe->Draw("pSAME");
    
    TGraph *gcdf = new TGraphErrors(1,scdf,jcdf,0,jecdf); 
    gcdf->SetMarkerStyle(21);gcdf->SetMarkerColor(6); gcdf->SetMarkerSize(1.7);
    gcdf->SetLineStyle(1); gcdf->SetLineColor(6); gcdf->SetLineWidth(1);
    gcdf->Draw("pSAME");

    TGraph *gali = new TGraphErrors(1,sali,sj2e,0,se3j2e);  //all err.
    gali->SetMarkerStyle(20);gali->SetMarkerColor(2); gali->SetMarkerSize(1.5);
    gali->SetLineStyle(1); gali->SetLineColor(2); gali->SetLineWidth(1);
    gali->Draw("pSAME");

    leg = new TLegend(.13,0.77,.7,0.94);
    //leg->AddEntry(gali,"ALICE, |y|<0.88, preliminary","p"); //PRELIM.
    leg->AddEntry(gali,"ALICE, |y|<0.9","p");
    leg->AddEntry(gcdf,"CDF, |y|<0.6","p");
    leg->AddEntry(gphe,"PHENIX, |y|<0.35","p");
    leg->SetTextFont(42);    leg->SetBorderSize(0);
    leg->SetFillStyle(0);    leg->SetMargin(0.15); //separation symbol-text
    leg->SetEntrySeparation(0.1);    leg->Draw();

    b1->SaveAs(Form("jpsi_s_spd%1d.eps",SPD));

}

