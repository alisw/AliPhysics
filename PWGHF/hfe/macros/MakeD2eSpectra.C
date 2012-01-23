#include <TF1.h>
#include <TFormula.h>
void MakeD2eSpectra(){

  // load libs
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("/alidata10/alice_u/minjung/032.heavy/hfe/bg4conv/util/hfe/libPWG3hfeDevel.so");

  setGeneralStyle();

  // read input files
  //char filename0[]="/alidata10/alice_u/minjung/train/GSI/V006.MC_pp/2011-03-08_2238.5742/LHC10f6a/HFEtask.root"; // pwg3 talk(linear binning)
  //char filename1[]="../TableOldnorm/HFPtSpectrum_D0Kpi_rebinnedth_combinedFD_021210_DecayBRSubtracted.root";     // pwg3 talk(use of 71.4mb)
  //char filename2[]="../TableOldnorm/HFPtSpectrum_DplusKpipi_combinedFD_021210_DecayBRSubtracted.root";           // pwg3 talk(use of 71.4mb)
  //char filename0[]="/alidata10/alice_u/minjung/train/GSI/V006.MC_pp/2011-03-15_2121.5954/LHC10f6a/HFEtask.root";   // preliminary plot (log binning)
  char filename0[]="./HFEtaskLHC10f6awExtendedDpt.root";   // QM plot 
  char filename1[]="../TableOldnorm/HFPtSpectrum_D0Kpi_combinedFD_rebinnedth_150311_newsigma_woDecayBR.root";      // preliminary plot(use of 62.1mb)
  char filename2[]="../TableOldnorm/HFPtSpectrum_DplusKpipi_combinedFD_150311_newsigma_woDecayBR.root";            // preliminary plot(use of 62.1mb) 

  TFile *_file[3]; 
  _file[0] = TFile::Open(filename0);
  _file[1] = TFile::Open(filename1);
  _file[2] = TFile::Open(filename2);

  TList *tl = (TList *)_file[0]->Get("TRD_HFE_QA2")->FindObject("MCqa");

  //count # of events
  TList *tl_result = (TList *)_file[0]->Get("TRD_HFE_Results2");
  AliCFContainer *containermc = (AliCFContainer *) tl_result->FindObject("eventContainer");
  AliCFDataGrid *mcGrid1 = new AliCFDataGrid("mcgrid1","",*containermc,AliHFEcuts::kEventStepGenerated);
  AliCFDataGrid *mcGrid2 = new AliCFDataGrid("mcgrid2","",*containermc,AliHFEcuts::kEventStepRecNoCut);
  printf("mc # of entries kEventStepGenerated= %lf kEventStepRecNoCut= %lf\n",mcGrid1->GetEntries(),mcGrid2->GetEntries());
 
  Int_t nEvt = (Int_t) mcGrid1->GetEntries();
  Double_t nevtnorm = 1./float(nEvt); //normalize with total number of event 


  // cross section for the old D meson data
  //Double_t Xsectrion4e = 71400000000.;
  //Double_t Xsectrion4D = 71400.; // consider the data point is represented as umb
  Double_t Xsectrion4e = 62100000000.;
  Double_t Xsectrion4D = 62100.; // consider the data point is represented as umb
  

  TString kDType[7];
  kDType[0]="421";
  kDType[1]="411";
  kDType[2]="431";
  kDType[3]="4122";
  kDType[4]="4132";
  kDType[5]="4232";
  kDType[6]="4332";
  
  TString kDType2[4];
  kDType2[0]="";
  kDType2[1]="D0";
  kDType2[2]="Dp";
  kDType2[3]="Drest";

  TString hname;

  // 2D histo having D, e pt correlation
  TH2D *hDeCorr[4]; // 4 groups.  "", "D0", "Dp", "Drest"
  for(int i=0; i<4; i++){
    hname = "mcqa_barrel_PtCorr"+kDType2[i]+"_c";
    hDeCorr[i]=(TH2D*)tl->FindObject(hname);  
  }

  // D mesons' pt, y 2D histograms from PYTHIA
  TH2D *hDPYTHIA[7]; // 7 particles. 421:D0, 411:D+, 431:Ds, 4122, 4132, 4232, 4332
  for(int i=0; i<7; i++){
    hname = "Dmeson"+kDType[i];
    hDPYTHIA[i]=(TH2D*)tl->FindObject(hname);
  }

  // PYTHIA pt distribution in |y|<0.5 of 7 D mesons
  TH1D *hDPYTHIApt[7];
  for(int i=0; i<7; i++){
    hname= "hD"+kDType[i]+"PYTHIApt";
    hDPYTHIApt[i]=(TH1D*)hDPYTHIA[i]->ProjectionX(hname,5,5); // |rapidity|<0.5 

    CorrectFromTheWidth(hDPYTHIApt[i]); // bin width correction
    hDPYTHIApt[i]->Scale(nevtnorm);     // # of event normalization
    hDPYTHIApt[i]->Scale(0.5);          // (particle + antiparticle)/2
    hDPYTHIApt[i]->Scale(Xsectrion4D);  // convert to cross section
    hDPYTHIApt[i]->SetMarkerStyle(20+i);
    hDPYTHIApt[i]->GetYaxis()->SetTitle("d#sigma/dp_{t} |_{|y|<0.5} [#mu b/GeV/c]");
  }

  // D measured (be careful with binnings)
  TH1D *hD0MeasuredStat= (TH1D *)_file[1]->Get("histoSigmaCorr");                        // data & statistic errors, start from bin number 4
  TGraphAsymmErrors *gD0measuredSys= (TGraphAsymmErrors *)_file[1]->Get("gSigmaCorr");   // data & systematic errors, start from bin number 4

  TH1D *hDpMeasuredStat = (TH1D *)_file[2]->Get("histoSigmaCorr");                       // data & statistic errors, start from bin number 1
  TGraphAsymmErrors *gDpMeasuredSys = (TGraphAsymmErrors *)_file[2]->Get("gSigmaCorr");  // data & statistic errors, start from bin number 1

  // drawing -------------------------------------------------------
  cPtDmeson = new TCanvas("cPtDmeson","pT of D mesons",0,0,800,500);
  cPtDmeson->Divide(2,1);

  cPtDmeson->cd(1);
  hDPYTHIApt[0]->Draw(""); //0: D0
  hD0MeasuredStat->Draw("same");
  TLegend *legD0 = new TLegend(0.60,0.75,0.75,0.85);
  legD0->AddEntry(hDPYTHIApt[0],"D^{0} PYTHIA","p");
  legD0->AddEntry(hD0MeasuredStat,"D^{0} Measured","p");
  setLegendStyle(*legD0,0);
  legD0->Draw("same");

  cPtDmeson->cd(2);
  hDPYTHIApt[1]->Draw(""); //1: Dp
  hDpMeasuredStat->Draw("same");
  TLegend *legDp = new TLegend(0.60,0.75,0.75,0.85);
  legDp->AddEntry(hDPYTHIApt[1],"D^{+} PYTHIA","p");
  legDp->AddEntry(hDpMeasuredStat,"D^{+} Measured","p");
  setLegendStyle(*legDp,0);
  legDp->Draw("same");
  //----------------------------------------------------------------


  Double_t ptval[6];
  Double_t XsecD0[6], XsecD0ErrStatFrac[6], XsecD0ErrSysFracLow[6], XsecD0ErrSysFracHigh[6];
  Double_t XsecDp[6], XsecDpErrStatFrac[6], XsecDpErrSysFracLow[6], XsecDpErrSysFracHigh[6];
  Double_t iM2PD0[6], iM2PDp[6];

  Double_t xbins[12]={0.,0.5,1.,2.,3.,4.,5.,6.,8.,12.,16.,20.};
  TH1D *hM2PD0= new TH1D("hM2PD0","D0: Measured/Pythia",11,xbins);
  TH1D *hM2PDp= new TH1D("hM2PDp","Dp: Measured/Pythia",11,xbins);

  for(int i=0; i<6; i++){ // 0 -> first measured point, there are 6 points for pp for the moment
     // D0
     iM2PD0[i]=0.; 
     ptval[i]=hD0MeasuredStat->GetBinCenter(i+4); // +4: starting from the first measured point
     XsecD0[i]=hD0MeasuredStat->GetBinContent(i+4);
     if(hDPYTHIApt[0]->GetBinContent(i+4)>0) iM2PD0[i]=hD0MeasuredStat->GetBinContent(i+4)/hDPYTHIApt[0]->GetBinContent(i+4);
     if(hD0MeasuredStat->GetBinContent(i+4)>0){
       XsecD0ErrStatFrac[i]    = hD0MeasuredStat->GetBinError(i+4)/hD0MeasuredStat->GetBinContent(i+4);
       XsecD0ErrSysFracLow[i]  = gD0measuredSys->GetErrorYlow(i+4)/hD0MeasuredStat->GetBinContent(i+4);
       XsecD0ErrSysFracHigh[i] = gD0measuredSys->GetErrorYhigh(i+4)/hD0MeasuredStat->GetBinContent(i+4);
     }

     // D+
     iM2PDp[i]=0.; 
     XsecDp[i]=hDpMeasuredStat->GetBinContent(i+1); // be careful with the fact that D+ measured starting from bin 1
     if(hDPYTHIApt[1]->GetBinContent(i+4)>0) iM2PDp[i]=hDpMeasuredStat->GetBinContent(i+1)/hDPYTHIApt[1]->GetBinContent(i+4); //D+
     if(hDpMeasuredStat->GetBinContent(i+1)>0) {
       XsecDpErrStatFrac[i]    = hDpMeasuredStat->GetBinError(i+1)/hDpMeasuredStat->GetBinContent(i+1);
       XsecDpErrSysFracLow[i]  = gDpMeasuredSys->GetErrorYlow(i+1)/hDpMeasuredStat->GetBinContent(i+1);
       XsecDpErrSysFracHigh[i] = gDpMeasuredSys->GetErrorYhigh(i+1)/hDpMeasuredStat->GetBinContent(i+1);
     }

     hM2PD0->Fill(ptval[i],iM2PD0[i]);
     hM2PDp->Fill(ptval[i],iM2PDp[i]);
  }
  // drawing -------------------------------------------------------
  cM2PDmeson = new TCanvas("cM2PDmeson","Measured/Pythia",0,0,800,500);	
  cM2PDmeson->Divide(2,1);
  cM2PDmeson->cd(1);
  hM2PD0->SetMarkerStyle(20);
  hM2PD0->SetXTitle("p_{t} [GeV/c]");
  hM2PD0->SetYTitle("ratio");
  hM2PD0->Draw("p");
  TLegend *legM2PD0 = new TLegend(0.45,0.75,0.75,0.85);
  legM2PD0->AddEntry(hM2PD0,"D0: Measured/PYTHIA","p");
  setLegendStyle(*legM2PD0,0);
  legM2PD0->Draw("same");

  cM2PDmeson->cd(2);
  hM2PDp->SetMarkerStyle(20);
  hM2PDp->SetXTitle("p_{t} [GeV/c]");
  hM2PDp->SetYTitle("ratio");
  hM2PDp->Draw("p");
  TLegend *legM2PDp = new TLegend(0.45,0.75,0.75,0.85);
  legM2PDp->AddEntry(hM2PDp,"D+: Measured/PYTHIA","p");
  setLegendStyle(*legM2PDp,0);
  legM2PDp->Draw("same");
  //----------------------------------------------------------------


  // electrons pt spectra from given D meson pt bins (before any scaling)
  TString kDPtbin[14];
  kDPtbin[0]="Dpt0";
  kDPtbin[1]="Dpt05";
  kDPtbin[2]="Dpt1";
  kDPtbin[3]="Dpt2";
  kDPtbin[4]="Dpt3";
  kDPtbin[5]="Dpt4";
  kDPtbin[6]="Dpt5";
  kDPtbin[7]="Dpt6";
  kDPtbin[8]="Dpt8";
  kDPtbin[9]="Dpt12";
  kDPtbin[10]="Dpt16";
  kDPtbin[11]="Dpt20";
  kDPtbin[12]="Dpt30";
  kDPtbin[13]="Dpt40";

  Int_t kbinl[14]={1, 6,11,21,31,41,51,61, 81,121,161,201,301,401};
  Int_t kbinh[14]={5,10,20,30,40,50,60,80,120,160,200,300,400,500};

  TH1D *hD2e[4][14];
  TH1D *hD2eNorm[4][14];
  TH1D *hD2eNormWm[2][14];
  Int_t nRebin = 1;
  for(int iDtype=0; iDtype<2; iDtype++){ // 4 groups.  0:"", 1:"D0", 2:"Dp", 3:"Drest"
    for(int iDptbin=0; iDptbin<14; iDptbin++){
      hname = "h"+kDType2[iDtype+1]+kDPtbin[iDptbin];
      hD2e[iDtype][iDptbin] = (TH1D*)hDeCorr[iDtype+1]->ProjectionY(hname,kbinl[iDptbin],kbinh[iDptbin]); 
      hD2e[iDtype][iDptbin]->Rebin(nRebin);
      CorrectFromTheWidth(hD2e[iDtype][iDptbin]);
      hname = "hNorm"+kDType2[iDtype]+kDPtbin[iDptbin];
      hD2eNorm[iDtype][iDptbin] = (TH1D*)hD2e[iDtype][iDptbin]->Clone(hname);
      hD2eNormWm[iDtype][iDptbin] = (TH1D*)hD2e[iDtype][iDptbin]->Clone(hname);
      if(iDtype==0) hD2eNorm[2][iDptbin] = (TH1D*)hD2e[iDtype][iDptbin]->Clone(hname); // for Ds, use D0 electron spectra
      if(iDtype==0) hD2eNorm[3][iDptbin] = (TH1D*)hD2e[iDtype][iDptbin]->Clone(hname); // for Lc, use D0 electron spectra
    }
  }

  hD0pt0to2=(TH1D*)hD2e[0][0]->Clone("hD0pt0to2");
  hD0pt2to12=(TH1D*)hD2e[0][3]->Clone("hD0pt2to12");
  hD0pt12to50=(TH1D*)hD2e[0][9]->Clone("hD0pt12to50");
  hD0pt0to2->Add(hD2e[0][1]);
  hD0pt0to2->Add(hD2e[0][2]);
  hD0pt2to12->Add(hD2e[0][4]);
  hD0pt2to12->Add(hD2e[0][5]);
  hD0pt2to12->Add(hD2e[0][6]);
  hD0pt2to12->Add(hD2e[0][7]);
  hD0pt2to12->Add(hD2e[0][8]);
  hD0pt12to50->Add(hD2e[0][10]);
  hD0pt12to50->Add(hD2e[0][11]);
  hD0pt12to50->Add(hD2e[0][12]);
  hD0pt12to50->Add(hD2e[0][13]);
  hD0ptall=(TH1D*)hD0pt0to2->Clone("hD0ptall");
  hD0ptall->Add(hD0pt2to12);
  hD0ptall->Add(hD0pt12to50);

  hD0lowptcontrib=(TH1D*)hD0pt0to2->Clone("hD0lowptcontrib");
  hD0lowptcontrib->Divide(hD0pt2to12);

  hD0pt0to2->Divide(hD0ptall);
  hD0pt2to12->Divide(hD0ptall);
  hD0pt12to50->Divide(hD0ptall);

  new TCanvas;
  hD0pt0to2->SetMarkerStyle(24);
  hD0pt2to12->SetMarkerStyle(24);
  hD0pt12to50->SetMarkerStyle(24);
  hD0pt0to2->SetMarkerColor(1);
  hD0pt0to2->SetLineColor(1);
  hD0pt2to12->SetMarkerColor(4);
  hD0pt2to12->SetLineColor(4);
  hD0pt12to50->SetMarkerColor(2);
  hD0pt12to50->SetLineColor(2);
  hD0pt0to2->Draw(); 
  hD0pt2to12->Draw("same");
  hD0pt12to50->Draw("same");
  gPad->SetLogy();
  gPad->SetGrid();

  // drawing -------------------------------------------------------
  cPtD2e = new TCanvas("cPtD2e","pT of e from certain D pt bin",0,0,1000,700);	
  cPtD2e->Divide(2,1);
  Int_t colorcode=1;
  TLegend *legD2e[4];
  for(int iDtype=0; iDtype<2; iDtype++){ 
    legD2e[iDtype] = new TLegend(0.75,0.45,0.88,0.88);
    legD2e[iDtype]->AddEntry(hD2e[0][0],kDType2[iDtype+1],"");
    for(int iDptbin=0; iDptbin<14; iDptbin++){
      cPtD2e->cd(iDtype+1); 
      colorcode=iDptbin+1;
      if(colorcode==10) colorcode=38;
      hD2e[iDtype][iDptbin]->SetLineColor(colorcode);
      hD2e[iDtype][iDptbin]->SetMarkerColor(colorcode);
      hD2e[iDtype][iDptbin]->SetMarkerStyle(20);
      if(iDptbin==0) hD2e[iDtype][iDptbin]->Draw();
      else hD2e[iDtype][iDptbin]->Draw("samep");
      gPad->SetLogy();
      gPad->SetLogx();
      legD2e[iDtype]->AddEntry(hD2e[iDtype][iDptbin],kDPtbin[iDptbin],"p");
      setLegendStyle(*legD2e[iDtype],0);
    }
    legD2e[iDtype]->Draw("same");
  }
  //----------------------------------------------------------------

  // electrons pt spectra from given D meson pt bins "with scaling based on measured X section"
  Double_t accfactorD[4][6] = { {1.73018, 1.72687, 1.76262, 1.94291, 1.82095, 1.69734,}, // acceptance factor of D0->e
                                {1.77726, 1.72664, 1.79331, 1.73913, 1.71084, 1.74054,}, // acceptance factor of Dp->e
                                {1.77726, 1.72664, 1.79331, 1.73913, 1.71084, 1.74054,}, // acceptance factor of Ds->e(copy of Dp)
                                {1.77726, 1.72664, 1.79331, 1.73913, 1.71084, 1.74054} }; // acceptance factor of Lc->e(copy of Dp)
  Double_t ptbwdth[6] = {1., 1., 1., 1., 2., 4.}; // pt bin width for binwidth correction

  // normalize PYTHIA electron spectra with the measured D meson cross section * branching ratio
  Double_t brachratioD[4];
  brachratioD[0] = 0.0653; // D0->e branching ratio
  brachratioD[1] = 0.16;   // Dp->e branching ratio
  brachratioD[2] = 0.08;   // Ds->e branching ratio
  brachratioD[3] = 0.045;  // Lc->e branching ratio
  //Double_t DiffMthdScaleD0= 1./1.07;  // considering difference between weighting method and X section method for D0 : new
  //Double_t DiffMthdScaleDp= 1./1.19;  // considering difference between weighting method and X section method for Dp : new
  Double_t DiffMthdScaleD0= 1./1.12;  // considering difference between weighting method and X section method for D0
  Double_t DiffMthdScaleDp= 1./1.25;  // considering difference between weighting method and X section method for Dp


  Double_t norm, xsec;
  Double_t xsecD[4][6];
  for(int iptb=0; iptb<6; iptb++){
    xsecD[0][iptb] = XsecD0[iptb]; //D0
    xsecD[1][iptb] = XsecDp[iptb]; //Dp
    xsecD[2][iptb] = (XsecD0[iptb]+XsecDp[iptb])*2.37/(8.49+2.65+5.07); //Ds: use of ZEUS ratio
    xsecD[3][iptb] = (XsecD0[iptb]+XsecDp[iptb])*3.59/(8.49+2.65+5.07); //Lc: use of ZEUS ratio 
  }

  for(int iDtype=0; iDtype<4; iDtype++){
    for(int iptb=0; iptb<6; iptb++){
      xsec = xsecD[iDtype][iptb]*1000000.; // change scale from ub to pb
      norm = hD2eNorm[iDtype][iptb+3]->Integral("width")/(xsec*brachratioD[iDtype]*ptbwdth[iptb]); //D0
      norm = accfactorD[iDtype][iptb]/norm;
      hD2eNorm[iDtype][iptb+3]->Scale(norm);
    }
  }

  // not measured pt data points for D0
  hD2eNorm[0][9]->Scale(nevtnorm);
  hD2eNorm[0][9]->Scale(0.44);  // measured/pythia for this pt bin
  hD2eNorm[0][9]->Scale(0.5);
  hD2eNorm[0][9]->Scale(DiffMthdScaleD0);
  hD2eNorm[0][9]->Scale(Xsectrion4e);
  hD2eNorm[0][10]->Scale(nevtnorm);
  hD2eNorm[0][10]->Scale(0.44); // measured/pythia for this pt bin
  hD2eNorm[0][10]->Scale(0.5);
  hD2eNorm[0][10]->Scale(DiffMthdScaleD0);
  hD2eNorm[0][10]->Scale(Xsectrion4e);
  hD2eNorm[0][11]->Scale(nevtnorm);
  hD2eNorm[0][11]->Scale(0.44); // measured/pythia for this pt bin
  hD2eNorm[0][11]->Scale(0.5);
  hD2eNorm[0][11]->Scale(DiffMthdScaleD0);
  hD2eNorm[0][11]->Scale(Xsectrion4e);
  hD2eNorm[0][12]->Scale(nevtnorm);
  hD2eNorm[0][12]->Scale(0.44); // measured/pythia for this pt bin
  hD2eNorm[0][12]->Scale(0.5);
  hD2eNorm[0][12]->Scale(DiffMthdScaleD0);
  hD2eNorm[0][12]->Scale(Xsectrion4e);
  hD2eNorm[0][13]->Scale(nevtnorm);
  hD2eNorm[0][13]->Scale(0.44); // measured/pythia for this pt bin
  hD2eNorm[0][13]->Scale(0.5);
  hD2eNorm[0][13]->Scale(DiffMthdScaleD0);
  hD2eNorm[0][13]->Scale(Xsectrion4e);


  // not measured pt data points for D+
  hD2eNorm[1][9]->Scale(nevtnorm); 
  hD2eNorm[1][9]->Scale(0.41);  // measured/pythia for this pt bin
  hD2eNorm[1][9]->Scale(0.5);
  hD2eNorm[1][9]->Scale(DiffMthdScaleDp);
  hD2eNorm[1][9]->Scale(Xsectrion4e);
  hD2eNorm[1][10]->Scale(nevtnorm);
  hD2eNorm[1][10]->Scale(0.41); // measured/pythia for this pt bin
  hD2eNorm[1][10]->Scale(0.5);
  hD2eNorm[1][10]->Scale(DiffMthdScaleDp);
  hD2eNorm[1][10]->Scale(Xsectrion4e);
  hD2eNorm[1][11]->Scale(nevtnorm);
  hD2eNorm[1][11]->Scale(0.41); // measured/pythia for this pt bin
  hD2eNorm[1][11]->Scale(0.5);
  hD2eNorm[1][11]->Scale(DiffMthdScaleDp);
  hD2eNorm[1][11]->Scale(Xsectrion4e);
  hD2eNorm[1][12]->Scale(nevtnorm);
  hD2eNorm[1][12]->Scale(0.41); // measured/pythia for this pt bin
  hD2eNorm[1][12]->Scale(0.5);
  hD2eNorm[1][12]->Scale(DiffMthdScaleDp);
  hD2eNorm[1][12]->Scale(Xsectrion4e);
  hD2eNorm[1][13]->Scale(nevtnorm);
  hD2eNorm[1][13]->Scale(0.41); // measured/pythia for this pt bin
  hD2eNorm[1][13]->Scale(0.5);
  hD2eNorm[1][13]->Scale(DiffMthdScaleDp);
  hD2eNorm[1][13]->Scale(Xsectrion4e);


  // estimation for Ds and Lc based on D0, D+ at high pt
  Double_t fragRatio=2.37/(8.49+2.65+5.07);
  Double_t brachRatio[2];
  Double_t brachRatio[0]=2*brachratioD[2]/(brachratioD[0]+brachratioD[1]);
  Double_t brachRatio[1]=2*brachratioD[3]/(brachratioD[0]+brachratioD[1]);
  for(int iptb=9; iptb<14; iptb++){
    hD2eNorm[2][iptb]->Add(hD2eNorm[0][iptb]);
    hD2eNorm[2][iptb]->Add(hD2eNorm[1][iptb]);
    hD2eNorm[2][iptb]->Scale(fragRatio);
    hD2eNorm[2][iptb]->Scale(brachRatio[0]);
    hD2eNorm[3][iptb]->Add(hD2eNorm[0][iptb]);
    hD2eNorm[3][iptb]->Add(hD2eNorm[1][iptb]);
    hD2eNorm[3][iptb]->Scale(fragRatio);
    hD2eNorm[3][iptb]->Scale(brachRatio[1]);
  }


  int kptmaxb = 11; // D0, D+
  TH1D *hD2eNormSum[4];
  for(int iDtype=0; iDtype<4; iDtype++){
    hname="hD2eNormSum" + iDtype;
    hD2eNormSum[iDtype]=(TH1D*)hD2eNorm[iDtype][3]->Clone(hname);
//    if(iDtype==2 || iDtype==3) kptmaxb=6; //Ds, Lc    //mmjj
    for(int iptb=1; iptb<kptmaxb; iptb++){
      hD2eNormSum[iDtype]->Add(hD2eNorm[iDtype][iptb+3]);
    }
  }

  for(int i=1; i<hD2eNormSum[0]->FindBin(2.0); i++){ // consider low pt contribution
     for(int iDtype=0; iDtype<4; iDtype++){
       Double_t icontrib =hD0lowptcontrib->GetBinContent(i)*hD2eNormSum[iDtype]->GetBinContent(i);
       Double_t ipt =hD2eNormSum[iDtype]->GetBinCenter(i);
       hD2eNormSum[iDtype]->Fill(ipt, icontrib);
     }
  }


  TH1D *hD2eNormSumStat[4]; // copy of hD2eNormSum[4] to get mc statistical errors
  for(int iDtype=0; iDtype<4; iDtype++){
    hname="hD2eNormSumStat" + iDtype;
    hD2eNormSumStat[iDtype]=(TH1D*)hD2eNormSum[iDtype]->Clone(hname);
  }

  // sum
  hD2eNormTotalStat=(TH1D*)hD2eNormSumStat[0]->Clone("hD2eNormTotalStat");
  hD2eNormTotalStat->Add(hD2eNormSumStat[1]);
  hD2eNormTotalStat->Add(hD2eNormSumStat[2]);
  hD2eNormTotalStat->Add(hD2eNormSumStat[3]);

  // drawing, here the error bars are only for MC statistics ---------------------------
  cPtD2eSum = new TCanvas("cPtD2eSum","pT of e from D",0,0,1000,600);
  cPtD2eSum->Divide(2,1);
  cPtD2eSum->cd(1);
  hD2eNormSumStat[0]->GetYaxis()->SetTitle("BRxd#sigma/dp_{t} | |#eta|<0.9  [pb/GeV/c]");
  hD2eNormSumStat[0]->SetMarkerStyle(20);
  hD2eNormSumStat[0]->Draw();
  hD2eNormSumStat[1]->SetMarkerColor(2);
  hD2eNormSumStat[1]->SetMarkerStyle(20);
  hD2eNormSumStat[1]->Draw("same");
  hD2eNormSumStat[2]->SetMarkerColor(4);
  hD2eNormSumStat[2]->SetMarkerStyle(20);
  hD2eNormSumStat[2]->Draw("same");
  hD2eNormSumStat[3]->SetMarkerColor(3);
  hD2eNormSumStat[3]->SetMarkerStyle(20);
  hD2eNormSumStat[3]->Draw("same");
  gPad->SetLogy();
  TLegend *legsum = new TLegend(0.25,0.70,0.75,0.85);
  legsum->AddEntry(hD2eNormSumStat[0],"D^{0}#rightarrow e 2.0<mother p_{t}<50.0 GeV/c","p");
  legsum->AddEntry(hD2eNormSumStat[1],"D^{+}#rightarrow e 2.0<mother p_{t}<50.0 GeV/c","p");
  legsum->AddEntry(hD2eNormSumStat[2],"D_{s}#rightarrow e 2.0<mother p_{t}<50.0 GeV/c","p");
  legsum->AddEntry(hD2eNormSumStat[3],"#Lambda_{c}#rightarrow e 2.0<mother p_{t}<50.0 GeV/c","p");
  setLegendStyle(*legsum,0);  legsum->Draw("same");

  cPtD2eSum->cd(2);
  hD2eNormTotalStat->SetMarkerStyle(20);
  hD2eNormTotalStat->GetYaxis()->SetTitle("BRxd#sigma/dp_{t} | |#eta|<0.9  [pb/GeV/c]");
  hD2eNormTotalStat->Draw();
  TLegend *legtotal = new TLegend(0.30,0.75,0.75,0.85);
  legtotal->AddEntry(hD2eNormTotalStat,"D#rightarrow e 2.0<mother p_{t}<50.0 GeV/c","p");
  setLegendStyle(*legtotal,0);
  legtotal->Draw("same");
  gPad->SetLogy();
  //----------------------------------------------------------------


  // calculating propagated error of D meson's statistic and systematic errors
  TGraphAsymmErrors *gD2eNormSum[4]; // store propagated error from D systematic errors
  for(iDtype=0; iDtype<4; iDtype++){
    gD2eNormSum[iDtype]= new TGraphAsymmErrors(hD2eNormSum[iDtype]);
  }

  Double_t statErr[11]; //kptmaxb
  Double_t sysErrLow[11]; //kptmaxb
  Double_t sysErrHigh[11]; //kptmaxb
  Double_t statErrSum[4]; // D0, D+, Ds. Lc
  Double_t sysErrLowSum[4]; // D0, D+, Ds. Lc
  Double_t sysErrHighSum[4]; // D0, D+, Ds. Lc
  Int_t ib4err = 0;
  for(int i=0; i<hD2eNormSum[0]->GetNbinsX(); i++){
     for(int iDtype=0; iDtype<4; iDtype++){ 
       statErrSum[iDtype]=0.;
       sysErrLowSum[iDtype]=0.;
       sysErrHighSum[iDtype]=0.;
     }

     for(int iptb=0; iptb<11; iptb++){
       // calculate propagated error from D statistical errors 
       if(iptb<6) ib4err = iptb;
       else ib4err = 5; // last two bins from pythia, so take the error from the last measured bin
       // statistical error of D propagated to systematic error of elec
       statErr[iptb]=hD2eNorm[0][iptb+3]->GetBinContent(i+1)*XsecD0ErrStatFrac[ib4err]; //D0
       if(!(iptb<6)) statErr[iptb] = statErr[iptb]*1.5; //1.5 put conservative error
       statErrSum[0] += statErr[iptb]*statErr[iptb];
       statErr[iptb]=hD2eNorm[1][iptb+3]->GetBinContent(i+1)*XsecDpErrStatFrac[ib4err]; //Dp
       if(!(iptb<6)) statErr[iptb] = statErr[iptb]*1.5;
       statErrSum[1] += statErr[iptb]*statErr[iptb];
       statErr[iptb]=hD2eNorm[2][iptb+3]->GetBinContent(i+1)*XsecDpErrStatFrac[ib4err]*1.5; //Ds
       if(!(iptb<6)) statErr[iptb] = statErr[iptb]*1.5;
       statErrSum[2] += statErr[iptb]*statErr[iptb];
       statErr[iptb]=hD2eNorm[3][iptb+3]->GetBinContent(i+1)*XsecDpErrStatFrac[ib4err]*1.5; //Lc
       if(!(iptb<6)) statErr[iptb] = statErr[iptb]*1.5;
       statErrSum[3] += statErr[iptb]*statErr[iptb];

       // systematic error of D propagated to systematic error of elec
       sysErrLow[iptb]=hD2eNorm[0][iptb+3]->GetBinContent(i+1)*XsecD0ErrSysFracLow[ib4err]; //D0
       if(!(iptb<6)) sysErrLow[iptb] = sysErrLow[iptb]*1.5;
       sysErrLowSum[0] += sysErrLow[iptb]*sysErrLow[iptb];
       sysErrLow[iptb]=hD2eNorm[1][iptb+3]->GetBinContent(i+1)*XsecDpErrSysFracLow[ib4err]; //Dp
       if(!(iptb<6)) sysErrLow[iptb] = sysErrLow[iptb]*1.5;
       sysErrLowSum[1] += sysErrLow[iptb]*sysErrLow[iptb];
       sysErrLow[iptb]=hD2eNorm[2][iptb+3]->GetBinContent(i+1)*XsecDpErrSysFracLow[ib4err]*1.5; //Ds
       if(!(iptb<6)) sysErrLow[iptb] = sysErrLow[iptb]*1.5;
       sysErrLowSum[2] += sysErrLow[iptb]*sysErrLow[iptb];
       sysErrLow[iptb]=hD2eNorm[3][iptb+3]->GetBinContent(i+1)*XsecDpErrSysFracLow[ib4err]*1.5; //Lc
       if(!(iptb<6)) sysErrLow[iptb] = sysErrLow[iptb]*1.5;
       sysErrLowSum[3] += sysErrLow[iptb]*sysErrLow[iptb];

       sysErrHigh[iptb]=hD2eNorm[0][iptb+3]->GetBinContent(i+1)*XsecD0ErrSysFracHigh[ib4err]; //D0
       if(!(iptb<6)) sysErrHigh[iptb] = sysErrHigh[iptb]*1.5;
       sysErrHighSum[0] += sysErrHigh[iptb]*sysErrHigh[iptb];
       sysErrHigh[iptb]=hD2eNorm[1][iptb+3]->GetBinContent(i+1)*XsecDpErrSysFracHigh[ib4err]; //Dp
       if(!(iptb<6)) sysErrHigh[iptb] = sysErrHigh[iptb]*1.5;
       sysErrHighSum[1] += sysErrHigh[iptb]*sysErrHigh[iptb];
       sysErrHigh[iptb]=hD2eNorm[2][iptb+3]->GetBinContent(i+1)*XsecDpErrSysFracHigh[ib4err]*1.5; //Ds
       if(!(iptb<6)) sysErrHigh[iptb] = sysErrHigh[iptb]*1.5;
       sysErrHighSum[2] += sysErrHigh[iptb]*sysErrHigh[iptb];
       sysErrHigh[iptb]=hD2eNorm[3][iptb+3]->GetBinContent(i+1)*XsecDpErrSysFracHigh[ib4err]*1.5; //Lc
       if(!(iptb<6)) sysErrHigh[iptb] = sysErrHigh[iptb]*1.5;
       sysErrHighSum[3] += sysErrHigh[iptb]*sysErrHigh[iptb];
     }
     
     for(iDtype=0; iDtype<4; iDtype++){
       hD2eNormSum[iDtype]->SetBinError(i+1,TMath::Sqrt(statErrSum[iDtype])); 
       gD2eNormSum[iDtype]->SetPointError(i,0,0,TMath::Sqrt(sysErrLowSum[iDtype]),TMath::Sqrt(sysErrHighSum[iDtype]));
     }
  }

  // sum -----------------------------------------------------------
  hD2eNormTotal=(TH1D*)hD2eNormSum[0]->Clone("hD2eNormTotal");
  hD2eNormTotal->Add(hD2eNormSum[1]);
  hD2eNormTotal->Add(hD2eNormSum[2]);
  hD2eNormTotal->Add(hD2eNormSum[3]);

  // sum up errors -------------------------------------------------
  TGraphAsymmErrors *gD2eNormTotal = new TGraphAsymmErrors(hD2eNormTotal);
  TGraphAsymmErrors *gD2eNormTotalwMCstat = new TGraphAsymmErrors(hD2eNormTotal);
  gD2eNormTotal->SetName("gD2eSys");
  gD2eNormTotal->SetTitle("electron spectra based on D measurement with systematic errors"); 
  gD2eNormTotalwMCstat->SetName("gD2eSysStat");
  gD2eNormTotalwMCstat->SetTitle("electron spectra based on D measurement with systematic and statistic errors"); 

  Double_t errLow[4], errHigh[4], errDStat[4], errMCStat[4];
  Double_t errLowSum, errHighSum, errDStatSum, errMCStatSum;
  for(int i=0; i<hD2eNormTotal->GetNbinsX(); i++){
    errLowSum = 0.;
    errHighSum = 0.;
    errDStatSum = 0.;
    errMCStatSum = 0.;
    
    for(int iDtype=0; iDtype<2; iDtype++){ 
      errLow[iDtype] = gD2eNormSum[iDtype]->GetErrorYlow(i); // error from D systematic error
      errLowSum += errLow[iDtype]; 

      errHigh[iDtype]= gD2eNormSum[iDtype]->GetErrorYhigh(i); // error from D systematic error
      errHighSum += errHigh[iDtype]; 

      errDStat[iDtype] = hD2eNormSum[iDtype]->GetBinError(i+1); //error from D statistical error
      errDStatSum += errDStat[iDtype]; 

      errMCStat[iDtype] = hD2eNormSumStat[iDtype]->GetBinError(i+1); //error from MC statistical error
      errMCStatSum += errMCStat[iDtype]; 
    }
    if(i<hD2eNormTotal->FindBin(2.0)){
      errDStatSum = (hD0lowptcontrib->GetBinContent(i+1)+1.)*errDStatSum; //put conservative error based on the yield fraction of low pt contribution
    } 

    gD2eNormTotal->SetPointError(i,0,0,TMath::Sqrt(errLowSum*errLowSum + errDStatSum*errDStatSum),TMath::Sqrt(errHighSum*errHighSum + errDStatSum*errDStatSum));
    gD2eNormTotalwMCstat->SetPointError(i,0,0,TMath::Sqrt(errLowSum*errLowSum + errDStatSum*errDStatSum + errMCStatSum*errMCStatSum),TMath::Sqrt(errHighSum*errHighSum + errDStatSum*errDStatSum + errMCStatSum*errMCStatSum));

  }

  // drawing ----------------------------------------------------------
  cPtD2eSumwErr = new TCanvas("cPtD2eSumwErr","pT of e from D",0,0,1000,600);
  cPtD2eSumwErr->cd(1);
  gD2eNormTotalwMCstat->SetLineColor(2);
  gD2eNormTotalwMCstat->SetMarkerColor(2);
  gD2eNormTotalwMCstat->SetMarkerStyle(20);
  gD2eNormTotalwMCstat->Draw("ACP");
  gD2eNormTotalwMCstat->GetXaxis()->SetTitle("p_{t} [GeV/c]");
  gD2eNormTotalwMCstat->GetYaxis()->SetTitle("BRxd#sigma/dp_{t} | |#eta|<0.9  [pb/GeV/c]");
  gD2eNormTotal->SetMarkerStyle(20);
  gD2eNormTotal->SetLineStyle(2);
  gD2eNormTotal->GetYaxis()->SetTitle("BRxd#sigma/dp_{t} | |#eta|<0.9  [pb/GeV/c]");
  gD2eNormTotal->Draw("Psame");
  TLegend *legtotalerr = new TLegend(0.50,0.75,0.75,0.85);
  legtotalerr->AddEntry(gD2eNormTotal,"D#rightarrow e 2.0<mother p_{t}<50.0 GeV/c","lp");  
  legtotalerr->AddEntry(gD2eNormTotalwMCstat,"(include MC stat error)","l");  
  setLegendStyle(*legtotalerr,0);
  legtotalerr->Draw("same");
  gPad->SetLogy();
  // ------------------------------------------------------------------

  // weighting method
  Double_t iM2PD[6];
  for(int iDtype=0; iDtype<2; iDtype++){
    for(int iptb=0; iptb<6; iptb++){
      if(iDtype==0) iM2PD[iptb]=iM2PD0[iptb];
      if(iDtype==1) iM2PD[iptb]=iM2PDp[iptb];
      hD2eNormWm[iDtype][iptb+3]->Scale(nevtnorm);
      hD2eNormWm[iDtype][iptb+3]->Scale(iM2PD[iptb]);  // measured/pythia for this pt bin
      hD2eNormWm[iDtype][iptb+3]->Scale(0.5);
      hD2eNormWm[iDtype][iptb+3]->Scale(Xsectrion4e);
    }
  }

  // sum up to 12 GeV/c for D0 and D+ to compare 2 different methods
  TH1D *hD2eNormSumM1[2];
  TH1D *hD2eNormSumM2[2];
  for(int iDtype=0; iDtype<2; iDtype++){
    hname="hD2eNormSumM1" + iDtype;
    hD2eNormSumM1[iDtype]=(TH1D*)hD2eNorm[iDtype][3]->Clone(hname);
    hname="hD2eNormSumM2" + iDtype;
    hD2eNormSumM2[iDtype]=(TH1D*)hD2eNormWm[iDtype][3]->Clone(hname);
    for(int iptb=1; iptb<6; iptb++){
      hD2eNormSumM1[iDtype]->Add(hD2eNorm[iDtype][iptb+3]);
      hD2eNormSumM2[iDtype]->Add(hD2eNormWm[iDtype][iptb+3]);
    }
  }
  // drawing -------------------------------------------------------
  c2Methods = new TCanvas("c2Methods","e spectra from 2 different methods",0,0,1000,700);
  c2Methods->Divide(3,1);
  TH1D* hRatio2mtd[2];
  TLegend *legM[3];
  for(int i=0; i<2; i++){
    c2Methods->cd(i+1);
    legM[i]= new TLegend(0.25,0.70,0.75,0.85);
    hD2eNormSumM1[i]->SetMarkerStyle(20);
    hD2eNormSumM1[i]->Draw();
    hD2eNormSumM1[i]->GetYaxis()->SetTitle("BRxd#sigma/dp_{t} | |#eta|<0.9  [pb/GeV/c]");
    hD2eNormSumM2[i]->SetMarkerColor(2);
    hD2eNormSumM2[i]->SetMarkerStyle(20);
    hD2eNormSumM2[i]->Draw("same");
    hname="hRatio2mtdD" + i;
    hRatio2mtd[i]=(TH1D*)hD2eNormSumM2[i]->Clone(hname);
    hRatio2mtd[i]->Divide(hD2eNormSumM1[i]);
    gPad->SetLogy();
  }
  legM[2]= new TLegend(0.25,0.70,0.75,0.85);
  legM[0]->AddEntry(hD2eNormSumM1[0],"2.0<mother p_{t}<12.0 GeV/c","");
  legM[0]->AddEntry(hD2eNormSumM1[0],"D^{0}#rightarrow e(use of Xsection)","p");
  legM[0]->AddEntry(hD2eNormSumM2[0],"D^{0}#rightarrow e(use of weighting)","p");
  legM[1]->AddEntry(hD2eNormSumM1[1],"2.0<mother p_{t}<12.0 GeV/c","");
  legM[1]->AddEntry(hD2eNormSumM1[1],"D^{+}#rightarrow e(use of Xsection)","p");
  legM[1]->AddEntry(hD2eNormSumM2[1],"D^{+}#rightarrow e(use of weighting)","p");

  legM[2]->AddEntry(hRatio2mtd[0],"2.0<mother p_{t}<12.0 GeV/c","");
  legM[2]->AddEntry(hRatio2mtd[0],"D^{0}#rightarrow e: (use of weighting)/(use of Xsection)","p");
  legM[2]->AddEntry(hRatio2mtd[1],"D^{+}#rightarrow e: (use of weighting)/(use of Xsection)","p");
  c2Methods->cd(1);
  legM[0]->Draw("same");
  c2Methods->cd(2);
  legM[1]->Draw("same");
  
  setLegendStyle(*legM[0],0);
  setLegendStyle(*legM[1],0);
  setLegendStyle(*legM[2],0);

  c2Methods->cd(3);
  hRatio2mtd[0]->Draw();
  hRatio2mtd[1]->SetLineColor(1);
  hRatio2mtd[1]->SetMarkerColor(1);
  hRatio2mtd[1]->Draw("same");
  legM[2]->Draw("same");
  //----------------------------------------------------------------


  // produce output file
  TFile *f = new TFile("Testoutput.root","recreate");
  f->cd();
  hD2eNormTotal->Write();
  gD2eNormTotal->Write();
  gD2eNormTotalwMCstat->Write();
  f->Close();


  // compare with FONLL ----------------------------------------------
  h4FONLL= new TH1F("h4FONLL","h4FONLL",60,0.25,30.25);
  ifstream stream1("./e_from_D-FONLL-abseta.lt.0.9.txt");

  double b;
  if(!stream1)
    cout << "While opening a file an error is encountered" << endl;

  int i=0;
  int ipt2=0;
  int k=0;
  Double_t val[8][60];
  Double_t exl[60];
  Double_t exh[60];
  Double_t eyl[60];
  Double_t eyh[60];
  while(!stream1.eof())
     {
         stream1 >> b;

         i=k%8;
         ipt2=int(double(k)/8.);
         val[i][ipt2] = b;
         k++;
         if(ipt2==59 && i==7) break;
     }

  for(int i=0; i<60; i++){
    h4FONLL->Fill(val[0][i],val[1][i]);
  }

  for(int m=0; m<60; m++){
    exl[m]=0.0;
    exh[m]=0.0;
    eyl[m]=val[1][m]-val[2][m];
    eyh[m]=val[3][m]-val[1][m];
  }

  TGraphAsymmErrors* gr=new  TGraphAsymmErrors(60,val[0],val[1],exl,exh,eyl,eyh);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  gr->SetMarkerSize(0.5);
  gr->SetLineColor(3);
  gr->SetFillColor(5);
  gr->GetXaxis()->SetTitle("p_{t} [GeV/c]");
  gr->GetYaxis()->SetTitle("BRxd#sigma/dp_{t} | |#eta|<0.9  [pb/GeV/c]");

  //gr->Draw("ZACPe3");
  //gr->Draw("ZCe3");
  //gr->Draw("Psame");
  cD2eFONLL = new TCanvas("cD2eFONLL","pT of e from D and FONLL",0,0,1000,600);
  cD2eFONLL->cd();
  gr->Draw("ACP");
  gD2eNormTotalwMCstat->Draw("psame");
  TLegend *legwFON = new TLegend(0.50,0.75,0.75,0.85);
  legwFON->AddEntry(gr,"D#rightarrow e FONLL","p");
  legwFON->AddEntry(gD2eNormTotalwMCstat,"D#rightarrow e 2.0<mother p_{t}<50.0 GeV/c","p");
  setLegendStyle(*legwFON,0);
  legwFON->Draw("same");
  gPad->SetLogy();

}

//--------------------------
void setDataStyle(TH1D &h, Int_t lc, Int_t lw, Int_t ls){

  h.SetLineColor(lc);
  h.SetLineWidth(lw);
  h.SetLineStyle(ls);
  h.SetMarkerColor(lc);

}

//--------------------------
void setLegendStyle(TLegend &legend, Int_t bs){

  legend.SetBorderSize(bs);
  legend.SetFillColor(0);
  legend.SetTextFont(132);
  legend.SetTextSize(0.04);
  legend.SetMargin(0.2);

}

//--------------------------
void setPadStyle(Int_t lvl, TPad *pad){

  pad->SetLogy();
  if(lvl>0) gPad->SetGridy();
  if(lvl>1) gPad->SetGridx();

}

//--------------------------
void setGeneralStyle(){

  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderSize(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleFont(132,"X");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(0.9);
  gStyle->SetLabelFont(132,"X");
  gStyle->SetLabelFont(132,"Y");
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetLineWidth(2);

}

void CorrectFromTheWidth(TH1D *h1) const {
  //
  // Correct from the width of the bins --> dN/dp_{T} (GeV/c)^{-1}
  //

  TAxis *axis = h1->GetXaxis();
  Int_t nbinX = h1->GetNbinsX();
  
  for(Int_t i = 1; i <= nbinX; i++) {
  
    Double_t width = axis->GetBinWidth(i);
    Double_t content = h1->GetBinContent(i);
    Double_t error = h1->GetBinError(i);
    h1->SetBinContent(i,content/width);
    h1->SetBinError(i,error/width);
    Double_t pt = h1->GetBinCenter(i);
  }
  
}
