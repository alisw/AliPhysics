#include "/Users/rashmiraniwala/ToALICE/ANALYSIS/SetStyle.C"
#include "/Users/rashmiraniwala/ToALICE/NumericalSim/readSimData.C"
ifstream fdata[8][5];
ifstream fTrk[8], fV0M[8], fVz5[8], fVz3[8];
Char_t fname1[256],fname2[256],fname[256],fnameData[256];
Float_t fpol1_TrkV01_Par1[8][5], fpol3_TrkV01_Par1[8][5], fpol3_TrkV01_Par2[8][5];
Float_t fpol1_TrkV01_dPar1[8][5], fpol3_TrkV01_dPar1[8][5], fpol3_TrkV01_dPar2[8][5];
Float_t fpol1_V01_Par1[8][5], fpol3_V01_Par1[8][5], fpol3_V01_Par2[8][5];
Float_t fpol1_V01_dPar1[8][5], fpol3_V01_dPar1[8][5], fpol3_V01_dPar2[8][5];
Float_t fpol1_Trk1_Par1[8][5], fpol3_Trk1_Par1[8][5], fpol3_Trk1_Par2[8][5];
Float_t fpol1_Trk1_dPar1[8][5], fpol3_Trk1_dPar1[8][5], fpol3_Trk1_dPar2[8][5];

Float_t fpol1_TrkV02_Par1[8][5], fpol3_TrkV02_Par1[8][5], fpol3_TrkV02_Par2[8][5];
Float_t fpol1_TrkV02_dPar1[8][5], fpol3_TrkV02_dPar1[8][5], fpol3_TrkV02_dPar2[8][5];
Float_t fpol1_V02_Par1[8][5], fpol3_V02_Par1[8][5], fpol3_V02_Par2[8][5];
Float_t fpol1_V02_dPar1[8][5], fpol3_V02_dPar1[8][5], fpol3_V02_dPar2[8][5];
Float_t fpol1_Trk2_Par1[8][5], fpol3_Trk2_Par1[8][5], fpol3_Trk2_Par2[8][5];
Float_t fpol1_Trk2_dPar1[8][5], fpol3_Trk2_dPar1[8][5], fpol3_Trk2_dPar2[8][5];
// To hold the systematic error due to different centrality i.e. TrkCentrality and V0Centrality
Float_t fpol1Par1Sys[8],fpol3Par1Sys[8], fpol3Par2Sys[8];
// To hold the systematic error due to different Vz cuts. 
Float_t fpol1Par1SysVz[8],fpol3Par1SysVz[8], fpol3Par2SysVz[8];
Float_t SysErrVzCentWeight[8]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
//Float_t SysErrVzCentWeight[8]={0.0};
Float_t SysErrfpol1fpol3[8]={0};

Int_t Reg1MinBin, Reg1MaxBin,Reg2MinBin, Reg2MaxBin,Reg3MinBin, Reg3MaxBin;
Int_t nevent1, nevent2, nevent3;
Float_t Reg12Limit, Reg3Limit;

Float_t etasigma = 3.6;
Float_t SlopeGl[8][5],dSlopeGl[8][5];

Float_t y0FromG[8]={0.0141288, 0.0293864, 0.0393269, 0.0465225, 0.0523195,0.0571578, 0.060533, 0.0637892};
Float_t dy0FromG[8]={4.31088e-05,8.13055e-05, 0.000121852, 0.000165566,0.000213216, 0.00026604,0.000321359,0.000391445};
Int_t imarker[10]={20,21,22,23,34,24,25,26,32,28};
Float_t Meanfpol3_TrkV0_Reg1_Par2[8], RMSfpol3_TrkV0_Reg1_Par2[8], dMeanfpol3_TrkV0_Reg1_Par2[8];
Float_t Meanfpol1_TrkV0_Reg1_Par1[8], RMSfpol1_TrkV0_Reg1_Par1[8], dMeanfpol1_TrkV0_Reg1_Par1[8];
Float_t Meanfpol3_TrkV0_Reg1_Par1[8], RMSfpol3_TrkV0_Reg1_Par1[8], dMeanfpol3_TrkV0_Reg1_Par1[8];
Float_t Meanfpol3_TrkV0_Reg2_Par2[8], RMSfpol3_TrkV0_Reg2_Par2[8], dMeanfpol3_TrkV0_Reg2_Par2[8];
Float_t Meanfpol1_TrkV0_Reg2_Par1[8], RMSfpol1_TrkV0_Reg2_Par1[8], dMeanfpol1_TrkV0_Reg2_Par1[8];
Float_t Meanfpol3_TrkV0_Reg2_Par1[8], RMSfpol3_TrkV0_Reg2_Par1[8], dMeanfpol3_TrkV0_Reg2_Par1[8];
TH1F * h_fpol1_TrkV0_Reg1[5];
TH1F * h_fpol3_TrkV0_Reg1_par1[5];
TH1F * h_fpol3_TrkV0_Reg1_par2[5];
TH1F * h_fpol1_TrkV0_Reg2[5];
TH1F * h_fpol3_TrkV0_Reg2_par1[5];
TH1F * h_fpol3_TrkV0_Reg2_par2[5];
TH1F *hMeanSlopeReg1, *hMeanSlopeReg2, *hMeanSlopeReg12WithStat, *hMeanSlopeReg12WithSys; 
Float_t meanfpol3par1[8], statfpol3par1[8];
TGraph * tg20, *tg07, *tg815;
TH1F * hc1High, *hc1Low;
Float_t fdummy;
Int_t dummy;

Float_t Final_fpol3par1_Reg1[8], Final_fpol3par1_Reg2[8];
Float_t Stat_fpol3par1_Reg1[8], Stat_fpol3par1_Reg2[8];

Float_t SysErrVz3_Reg1[8], SysErrVz3_Reg2[8]; 
Float_t SysErrTrkCent1[8], SysErrTrkCent2[8];
Float_t SysErrV0MCent1[8], SysErrV0MCent2[8], SysErrCent[8]={0};
Float_t SysErrReg3Width_Reg1[8], SysErrReg3Width_Reg2[8];
Float_t SysErrReg12[8];
Float_t Sys_fpol3par1_Reg1[8], Sys_fpol3par1_Reg2[8];
  
void DrawFitResultsv4(){
  
  SetStyle();
  //  readSimData();
  
  // Reads data files and fills information on fir parameters for data
  ReadData();
  cout<<" Read Data"<<endl;
  FillHistosData();
  
  // Draws plots for different reg ( 1 & 2) and different width of reg3 ( 0.04 to 0.2)
  DrawFitCoeff2();
  
  for(Int_t icent = 0;icent<8;icent++){
    // Fill Region 1 reference values
    Final_fpol3par1_Reg1[icent] = fpol3_TrkV01_Par1[icent][4];
    Stat_fpol3par1_Reg1[icent] = fpol3_TrkV01_dPar1[icent][4];
    // Fill Region 2 reference values
    Final_fpol3par1_Reg2[icent] = fpol3_TrkV02_Par1[icent][4];
    Stat_fpol3par1_Reg2[icent] = fpol3_TrkV02_dPar1[icent][4];
  }
    
  GetSysReg3Width();
  GetTrkCent_SysErr();
  GetV0Cent_SysErr();
  GetVzDiff_SystematicError();

  GetSysReg12();
  GetSysErrfpol1fpol3();
  // only symmetric errors are calculated here
  GetFinalSysErrors();
  // The asymmetric errors are added here. 
  DrawMeanStatSys();
  
  TFile * fout = new TFile("DataSlopeWithStatSys.root","RECREATE");
  fout->cd();
  hMeanSlopeReg12WithStat->Write();
  hMeanSlopeReg12WithSys->Write();
  tg20->Write();
  tg07->Write();
  tg815->Write();
  hc1High->Write();
  hc1Low->Write();
  fout->Close();
  
  //DrawMeanRMS_Reg12Mean();
  
}

void GetSysErrfpol1fpol3(){
  cout<<" Sys Error fpol1 fpol3 "<<endl;
  for(Int_t icent = 0;icent<8;icent++){
    SysErrfpol1fpol3[icent] = fpol1_TrkV01_Par1[icent][4] - fpol3_TrkV01_Par1[icent][4];
    cout<<icent<<"  "<<fpol1_TrkV01_Par1[icent][4]<<" "<<fpol3_TrkV01_Par1[icent][4]<<" "<<SysErrfpol1fpol3[icent]<<"     % "<<SysErrfpol1fpol3[icent]*100.0/fpol3_TrkV01_Par1[icent][4]<<endl;
  }

}


void DrawMeanStatSys(){
  Float_t cent[8]={2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5};
  Float_t ex[8]={2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5};
  Float_t SysLow[8], SysHigh[8];
  hMeanSlopeReg12WithStat = new TH1F("hMeanSlopeReg12WithStat","hMeanSlopeReg12WithStat",8,0.0,40.0);
  hMeanSlopeReg12WithSys = new TH1F("hMeanSlopeReg12WithSys","hMeanSlopeReg12WithSys",8,0.0,40.0);
  hc1High = new TH1F("hc1High","hc1High",8,0.0,40.0);
  hc1Low = new TH1F("hc1Low","hc1Low",8,0.0,40.0);

  tg20 = new TGraph(20);
  tg20->SetNameTitle("tg20ForSys","tg20ForSys");
  Float_t fpol3High[8];
  Float_t fpol3Low[8];
  for(Int_t icent=0;icent<8;icent++){
    hMeanSlopeReg12WithStat->SetBinContent(icent+1,TMath::Abs(Final_fpol3par1_Reg1[icent]));
    hMeanSlopeReg12WithStat->SetBinError(icent+1,Stat_fpol3par1_Reg1[icent]);
    //    hMeanSlopeReg12WithSys->SetBinContent(icent+1,TMath::Abs(Final_fpol3par1_Reg1[icent]));
    //    hMeanSlopeReg12WithSys->SetBinError(icent+1,Sys_fpol3par1_Reg1[icent]);
    SysLow[icent] = Sys_fpol3par1_Reg1[icent];
    SysHigh[icent] = Sys_fpol3par1_Reg1[icent]*Sys_fpol3par1_Reg1[icent] + SysErrfpol1fpol3[icent]*SysErrfpol1fpol3[icent];
    SysHigh[icent] = TMath::Sqrt(SysHigh[icent]);
    //cout<<SysHigh[icent]+Final_fpol3par1_Reg1[icent]<<" "<<Final_fpol3par1_Reg1[icent]-SysLow[icent]<<endl;
    
    fpol3High[icent] = SysHigh[icent]+Final_fpol3par1_Reg1[icent];
    fpol3Low[icent] = Final_fpol3par1_Reg1[icent]-SysLow[icent];
    hc1High->SetBinContent(icent+1,SysHigh[icent]+Final_fpol3par1_Reg1[icent]);
    hc1Low->SetBinContent(icent+1,Final_fpol3par1_Reg1[icent]-SysLow[icent]);
    //cout<<SysHigh[icent]<<" "<<SysLow[icent]<<" "<<SysHigh[icent]*100/Final_fpol3par1_Reg1[icent]<<" "<<SysLow[icent]*100/Final_fpol3par1_Reg1[icent]<<endl;
    cout<<"Higher % "<<SysHigh[icent]*100/Final_fpol3par1_Reg1[icent]<<" lower% "<<SysLow[icent]*100/Final_fpol3par1_Reg1[icent]<<endl;
}
  tg07 = new TGraph(8,cent,fpol3High);
  tg07->SetNameTitle("tg07","tg07");
  tg815 = new TGraph(8,cent,fpol3Low);
  tg815->SetNameTitle("tg815","tg815");
  
  hc1High->SetFillColor(kBlue-9);
  hc1High->SetLineColor(kBlue-9);
  hc1Low->SetFillColor(0);
  TCanvas * cMeanSys = new TCanvas("cMeanSys","cMeanSys",800,600);
  hMeanSlopeReg12WithStat->SetMarkerStyle(8);
  hMeanSlopeReg12WithStat->SetMaximum(0.006);
  hMeanSlopeReg12WithStat->SetMinimum(0.000);
  hMeanSlopeReg12WithStat->SetMarkerColor(4);
  hMeanSlopeReg12WithStat->SetLineColor(4);
  hMeanSlopeReg12WithStat->SetXTitle("Centrality %");
  hMeanSlopeReg12WithStat->SetYTitle("Slope of Ratio of dN/d#eta(c_{1})");
  hMeanSlopeReg12WithStat->Draw();
  /*
  hMeanSlopeReg12WithSys->SetMarkerStyle(8);
  hMeanSlopeReg12WithSys->SetMarkerColor(4);
  hMeanSlopeReg12WithSys->SetLineColor(4);
  hMeanSlopeReg12WithSys->SetFillColor(kBlue-9);
  hMeanSlopeReg12WithSys->SetFillStyle(3001);
  hMeanSlopeReg12WithSys->Draw("same E3");
  */
  tg20->SetFillStyle(3001);
  tg20->SetFillColor(kBlue-10);
  tg07->Fit("pol3");
  tg815->Fit("pol3");
  TF1 * f1 = tg07->GetFunction("pol3");
  TF1 * f2 = tg815->GetFunction("pol3");
  f1->SetLineColor(kBlue-9);
  f2->SetLineColor(kBlue-9);
  for(Int_t icent = 0;icent<8;icent++){
    Float_t y = f1(cent[icent]);
    tg20->SetPoint(icent+1,cent[icent],y);
  }
  for(Int_t icent = 0;icent<8;icent++){
    Float_t y = f2(cent[icent]);
    tg20->SetPoint(18-icent,cent[icent],y);
  }
  tg20->SetPoint(0,0.0,f1(0.0));
  tg20->SetPoint(9,40.0,f1(40.0));
  tg20->SetPoint(10,40.0,f2(40.0));
  tg20->SetPoint(19,0.0,f2(0.0));
  tg20->Draw("f");  
  TLatex * l = new TLatex(5,0.005,"Pb-Pb 2.76TeV");
  l->SetTextColor(4);
  l->Draw();
  l = new TLatex(20,0.005,"ALICE PRELIMINARY");
  l->SetTextColor(2);
  l->Draw();
  hMeanSlopeReg12WithStat->Draw("same");
  
}

void GetSysReg12(){

  for(Int_t icent =0;icent<8;icent++){
    SysErrReg12[icent] = TMath::Abs(Final_fpol3par1_Reg1[icent])-TMath::Abs(Final_fpol3par1_Reg2[icent]);
    SysErrReg12[icent]=TMath::Abs(SysErrReg12[icent]);
    cout<<" Reg12 Diff = "<<SysErrReg12[icent]<<" %= "<<SysErrReg12[icent]*100.0/Final_fpol3par1_Reg1[icent]<<endl;
  }
}

void GetSysReg3Width(){
  
  Float_t mean1,mean2;
  cout<<" Difference With fpol3Par1 for Reg3 = +/- 0.1"<<endl;
  for(Int_t icent = 0;icent<8;icent++){
    mean1 = Final_fpol3par1_Reg1[icent];
    mean2 = Final_fpol3par1_Reg2[icent];
    Float_t maxDiff = 0, diff;
    SysErrReg3Width_Reg1[icent]=0;
    SysErrReg3Width_Reg2[icent]=0;
    
    for(Int_t wreg = 0;wreg<4;wreg++){
      diff = fpol3_TrkV01_Par1[icent][wreg]-mean1;
      diff = TMath::Abs(diff);
      //    cout<<" Diff = "<<diff<<" MaxDiff="<<maxDiff<<endl;
      if(maxDiff<diff) maxDiff=diff;
      /*
      SysErrReg3Width_Reg1[icent] += diff*diff;
      diff = fpol3_TrkV02_Par1[icent][wreg]-mean2;
      SysErrReg3Width_Reg2[icent] += diff*diff;
      */
    }
    //    SysErrReg3Width_Reg1[icent]=TMath::Sqrt(SysErrReg3Width_Reg1[icent]/4);
    //    SysErrReg3Width_Reg2[icent]=TMath::Sqrt(SysErrReg3Width_Reg2[icent]/4);
    // This was overwritten to maximum difference instead of RMS of diff reg3 values about that for Reg3 = +/-0.1
    SysErrReg3Width_Reg1[icent]=maxDiff;
    cout<<"Re3=0.1 is "<< mean1 <<" Reg1  "<<TMath::Abs(SysErrReg3Width_Reg1[icent]*100/mean1)<<endl;
    //"%   Reg2 Reg3-0.1 is  "<<mean2<<"  "<<TMath::Abs(SysErrReg3Width_Reg2[icent]*100/mean2)<<" %   "<<endl; ; 
  }
  cout<<endl;
}

void GetFinalSysErrors(){

  cout<<" Sys Error on Reg1 = "<<endl;
  for(Int_t icent = 0;icent<8;icent++){
    SysErrCent[icent]=SysErrV0MCent1[icent];
    if(SysErrTrkCent1[icent]>SysErrCent[icent]) SysErrCent[icent]=SysErrTrkCent1[icent];
    Sys_fpol3par1_Reg1[icent]=0;
    Sys_fpol3par1_Reg1[icent]+=SysErrVz3_Reg1[icent]*SysErrVz3_Reg1[icent];
    //    Sys_fpol3par1_Reg1[icent]+=SysErrTrkCent1[icent]*SysErrTrkCent1[icent];
    //    Sys_fpol3par1_Reg1[icent]+=SysErrV0MCent1[icent]*SysErrV0MCent1[icent];
    Sys_fpol3par1_Reg1[icent]+=SysErrCent[icent]*SysErrCent[icent];
    //  Sys_fpol3par1_Reg1[icent]+=SysErrReg3Width_Reg1[icent]*SysErrReg3Width_Reg1[icent];
    Sys_fpol3par1_Reg1[icent]+=SysErrReg12[icent]*SysErrReg12[icent];
    Sys_fpol3par1_Reg1[icent]+=TMath::Power(SysErrVzCentWeight[icent]*Final_fpol3par1_Reg1[icent],2);
    Sys_fpol3par1_Reg1[icent] = TMath::Sqrt(Sys_fpol3par1_Reg1[icent]);
    cout<<icent<<" "<<Final_fpol3par1_Reg1[icent]<<" stat "<<Stat_fpol3par1_Reg1[icent]*100/Final_fpol3par1_Reg1[icent]<<"  Lower% "<<Sys_fpol3par1_Reg1[icent]*100/Final_fpol3par1_Reg1[icent]<<
      " Cent in% "<<SysErrCent[icent]*100/Final_fpol3par1_Reg1[icent]<<
      " Vz3% "<<SysErrVz3_Reg1[icent]*100/Final_fpol3par1_Reg1[icent]<<
      "  REg12Diff% "<<SysErrReg12[icent]*100/Final_fpol3par1_Reg1[icent]<<
      endl;
  }
  /*
  for(Int_t icent = 0;icent<8;icent++){
    Sys_fpol3par1_Reg2[icent]=0;
    Sys_fpol3par1_Reg2[icent]+=SysErrVz3_Reg2[icent]*SysErrVz3_Reg2[icent];
    Sys_fpol3par1_Reg2[icent]+=SysErrTrkCent2[icent]*SysErrTrkCent2[icent];
    Sys_fpol3par1_Reg2[icent]+=SysErrV0MCent2[icent]*SysErrV0MCent2[icent];
    Sys_fpol3par1_Reg2[icent]+=SysErrReg3Width_Reg2[icent]*SysErrReg3Width_Reg2[icent];
    Sys_fpol3par1_Reg2[icent] = TMath::Sqrt(Sys_fpol3par1_Reg2[icent]);
    cout<<icent<<" "<<Final_fpol3par1_Reg2[icent]<<" stat "<<Stat_fpol3par1_Reg2[icent]*100/Final_fpol3par1_Reg2[icent]<<"  % "<<Sys_fpol3par1_Reg2[icent]*100/Final_fpol3par1_Reg2[icent]<<endl;
  }
  */
}

void GetVzDiff_SystematicError(){

  Int_t dummy;

  Float_t fVz3pol11_Par1, fVz3pol12_Par1, fVz3pol31_Par1, fVz3pol32_Par1,fVz3pol31_Par2, fVz3pol32_Par2;
  Float_t fdummy;
  for(Int_t icent = 0;icent<8;icent++){

    sprintf(fname2,"/Users/rashmiraniwala/ToALICE/ANALYSIS/dNdEtaRatios/SystematicsVertexZ/dNdetaRatioFitData_Cent%3.3dto%3.3d_V0MTrk_TrkV0MCentrality_Reg3_%2.2f.dat",icent*5,(icent+1)*5,0.10);
    //  cout<<" Opening file "<<fname2<<endl;
    fVz3[icent].open(fname2,ios::in);
    fVz3[icent]>>fdummy>>fdummy;
    fVz3[icent]>>dummy>>dummy;
    fVz3[icent]>>dummy>>dummy;
    fVz3[icent]>>dummy>>dummy;
    fVz3[icent]>>dummy>>dummy>>dummy;
    // cout<<dummy<<endl;
    fVz3[icent]>>fVz3pol11_Par1>>fdummy;
    fVz3[icent]>>fVz3pol12_Par1>>fdummy;
    fVz3[icent]>>fVz3pol31_Par1>>fdummy>>fVz3pol31_Par2>>fdummy;
    fVz3[icent]>>fVz3pol32_Par1>>fdummy>>fVz3pol32_Par2>>fdummy;

    SysErrVz3_Reg1[icent]=TMath::Abs(TMath::Abs(Final_fpol3par1_Reg1[icent])-TMath::Abs(fVz3pol31_Par1));
    SysErrVz3_Reg2[icent]=TMath::Abs(TMath::Abs(Final_fpol3par1_Reg2[icent])-TMath::Abs(fVz3pol32_Par1));
    
  }
  cout<<"fpol3Par1 Diff between Vz5 and Vz3 for Region 1"<<endl;
  for(Int_t icent=0;icent<8;icent++){
    cout<<icent<<" "<<SysErrVz3_Reg1[icent]<<" %= "<<SysErrVz3_Reg1[icent]*100.0/Final_fpol3par1_Reg1[icent]<<endl;
  }
  /*
 cout<<"fpol3Par1 Diff between Vz5 and Vz3 for Region 2"<<endl;
  for(Int_t icent=0;icent<8;icent++){
    cout<<icent<<" "<<SysErrVz3_Reg2[icent]<<" %= "<<SysErrVz3_Reg2[icent]*100.0/Final_fpol3par1_Reg2[icent]<<endl;
  }
  */
}


void GetTrkCent_SysErr(){
  cout<<" Getting TrkCent_SysErr"<<endl;
  Int_t dummy;
  Float_t fpol1_Reg1_Par1, fpol1_Reg2_Par1, fpol3_Reg1_Par1, fpol3_Reg2_Par1,fpol3_Reg1_Par2, fpol3_Reg2_Par2;
  Float_t fdummy;
  

  
  for(Int_t icent = 0;icent<8;icent++){
    ifstream infile;
    sprintf(fname1,"/Users/rashmiraniwala/ToALICE/ANALYSIS/dNdEtaRatios/TrkCentrality/dNdetaRatioFitData_Cent%3.3dto%3.3d_Trk_TrkCentrality_Reg3_%2.2f.dat",icent*5,(icent+1)*5,0.10);
    // cout<<" Opening file "<<fname1<<endl;
    infile.open(fname1,ios::in);
    infile>>fdummy>>fdummy;
    infile>>dummy>>dummy;
    infile>>dummy>>dummy;
    infile>>dummy>>dummy;
    infile>>dummy>>dummy>>dummy;
    //   cout<<dummy<<endl;
    infile>>fpol1_Reg1_Par1>>fdummy;
    infile>>fpol1_Reg2_Par1>>fdummy;
    infile>>fpol3_Reg1_Par1>>fdummy>>fpol3_Reg1_Par2>>fdummy;
    infile>>fpol3_Reg2_Par1>>fdummy>>fpol3_Reg2_Par2>>fdummy;
    //    cout<<fpol3_Reg1_Par1<<" "<<fpol3_Reg2_Par1<<endl;
    SysErrTrkCent1[icent] =  TMath::Abs(TMath::Abs(Final_fpol3par1_Reg1[icent]) - TMath::Abs(fpol3_Reg1_Par1));
    
    SysErrTrkCent2[icent] =  TMath::Abs(TMath::Abs(Final_fpol3par1_Reg2[icent]) - TMath::Abs(fpol3_Reg2_Par1));
  }
  for(Int_t icent = 0;icent<8;icent++){
    cout<<icent<<" "<< SysErrTrkCent1[icent]<<" % "<< SysErrTrkCent1[icent]*100/Final_fpol3par1_Reg1[icent]<<endl;
  }
  /*
  for(Int_t icent = 0;icent<8;icent++){
    cout<<icent<<" "<< SysErrTrkCent2[icent]<<" % "<< SysErrTrkCent2[icent]*100/Final_fpol3par1_Reg2[icent]<<endl;
  }
  cout<<endl;
  */  
}


void GetV0Cent_SysErr(){

  cout<<" Getting V0Cent_SysErr"<<endl;
  Int_t dummy;
  Float_t fpol1_Reg1_Par1, fpol1_Reg2_Par1, fpol3_Reg1_Par1, fpol3_Reg2_Par1,fpol3_Reg1_Par2, fpol3_Reg2_Par2;
  Float_t fdummy;
  
  for(Int_t icent = 0;icent<8;icent++){
    ifstream infile;
    sprintf(fname1,"/Users/rashmiraniwala/ToALICE/ANALYSIS/dNdEtaRatios/V0MCentrality/dNdetaRatioFitData_Cent%3.3dto%3.3d_V0M_V0MCentrality_Reg3_%2.2f.dat",icent*5,(icent+1)*5,0.10);
    // cout<<" Opening file "<<fname1<<endl;
    infile.open(fname1,ios::in);
    infile>>fdummy>>fdummy;
    infile>>dummy>>dummy;
    infile>>dummy>>dummy;
    infile>>dummy>>dummy;
    infile>>dummy>>dummy>>dummy;
    //   cout<<dummy<<endl;
    infile>>fpol1_Reg1_Par1>>fdummy;
    infile>>fpol1_Reg2_Par1>>fdummy;
    infile>>fpol3_Reg1_Par1>>fdummy>>fpol3_Reg1_Par2>>fdummy;
    infile>>fpol3_Reg2_Par1>>fdummy>>fpol3_Reg2_Par2>>fdummy;
    //    cout<<fpol3_Reg1_Par1<<" "<<fpol3_Reg2_Par1<<endl;
    SysErrV0MCent1[icent] =  TMath::Abs(TMath::Abs(Final_fpol3par1_Reg1[icent]) - TMath::Abs(fpol3_Reg1_Par1));
    
    SysErrV0MCent2[icent] =  TMath::Abs(TMath::Abs(Final_fpol3par1_Reg2[icent]) - TMath::Abs(fpol3_Reg2_Par1));
  }
  for(Int_t icent = 0;icent<8;icent++){
    cout<<icent<<" "<< SysErrV0MCent1[icent]<<" % "<< SysErrV0MCent1[icent]*100/Final_fpol3par1_Reg1[icent]<<endl;
  }
  /*
  for(Int_t icent = 0;icent<8;icent++){
    cout<<icent<<" "<< SysErrV0MCent2[icent]<<" % "<< SysErrV0MCent2[icent]*100/Final_fpol3par1_Reg2[icent]<<endl;
  }
  cout<<endl;
  */
}


 
void DrawMeanRMS_Reg12Mean(){
  
  Float_t y3[5]={0.0018,0.00155,0.0013,0.00105,0.0008};
  Char_t mystr[48];
  TMarker *m;
  TLatex * l;
  
  TH1F *hMeanSlopeReg1 = new TH1F("hMeanSlopeReg1","hMeanSlopeReg1",8,0.0,40.0);
  TH1F *hMeanSlopeReg2 = new TH1F("hMeanSlopeReg2","hMeanSlopeReg2",8,0.0,40.0);
  TH1F *hMeanSlopeReg12WithStat = new TH1F("hMeanSlopeReg12WithStat","hMeanSlopeReg12WithStat",8,0.0,40.0);
  TH1F * hMeanSlopeReg12WithSys = new TH1F("hMeanSlopeReg12WithSys","hMeanSlopeReg12WithSys",8,0.0,40.0);
  
  for(Int_t icent=0;icent<8;icent++){
    
    hMeanSlopeReg1->SetBinContent(icent+1,Final_fpol3par1_Reg1[icent]);
    hMeanSlopeReg1->SetBinError(icent+1,Stat_fpol3par1_Reg1[icent]);
    hMeanSlopeReg2->SetBinContent(icent+1,-1.0*Final_fpol3par1_Reg2[icent]);
    hMeanSlopeReg2->SetBinError(icent+1,Stat_fpol3par1_Reg2[icent]);

  }
  TCanvas * cMean12 = new TCanvas("cMean12","cMean12",600,600);

  hMeanSlopeReg1->SetMarkerStyle(8);
  hMeanSlopeReg1->SetMaximum(0.007);
  hMeanSlopeReg1->SetMinimum(0.0005);
  hMeanSlopeReg1->SetMarkerColor(4);
  hMeanSlopeReg1->SetLineColor(4);
  hMeanSlopeReg1->SetXTitle("Centrality %");
  hMeanSlopeReg1->SetYTitle("c_{1}");
  hMeanSlopeReg1->Draw();
  hMeanSlopeReg2->SetMarkerStyle(8);
  hMeanSlopeReg2->SetMarkerColor(2);
  hMeanSlopeReg2->SetLineColor(2);
  hMeanSlopeReg2->Draw("same");
  h_c1Reg2_Fromfpol3par1_Sim->Draw(" same");

  TLatex * l = new TLatex(5,0.005,"Data");
  l->SetTextColor(4);
  l->Draw();
  l = new TLatex(5,0.0056,"MC Simulation");
  l->SetTextColor(6);
  l->Draw();
}
 
void Drawy0(){
  TCanvas * cy0 = new TCanvas("cy0","cy0",600,600);
  TH1F *hMeany0Reg12WithStat = (TH1F*)hMeanSlopeReg12WithStat->Clone("hMeany0Reg12WithStat");
  TH1F *hMeany0Reg12WithSys = (TH1F*)hMeanSlopeReg12WithSys->Clone("hMeany0Reg12WithSys");
  hMeany0Reg12WithStat->Reset();
  hMeany0Reg12WithSys->Scale(3.7*3.7);
  for(Int_t icent = 0;icent<8;icent++){
    Float_t y0 = 3.7*3.7*hMeanSlopeReg12WithStat->GetBinContent(icent+1);
    Float_t dy0 = y0*0.07;
    hMeany0Reg12WithStat->SetBinContent(icent+1,y0);
    hMeany0Reg12WithStat->SetBinError(icent+1, dy0);
    cout<<" y0 = "<<y0<<"+/-"<<dy0<<endl;
  }
  hMeany0Reg12WithStat->SetYTitle("Shift in Rapidity y_{0}");
  hMeany0Reg12WithStat->SetMinimum(0);
  hMeany0Reg12WithStat->SetMaximum(0.07);
  hMeany0Reg12WithStat->Draw();
  hMeany0Reg12WithSys->SetFillColor(4);
  hMeany0Reg12WithSys->SetFillStyle(3001);
  hMeany0Reg12WithSys->Draw("same E3");
  
  TH1F *h_y0Reg1_Fromfpol3par1_Sim = (TH1F*)h_c1Reg1_Fromfpol3par1_Sim->Clone("h_y0Reg1_Fromfpol3par1_Sim");
  TH1F *h_y0Reg1_Fromfpol3par1_SimWithSys = (TH1F*)h_c1Reg1_Fromfpol3par1_SimWithSys->Clone("h_y0Reg1_Fromfpol3par1_SimWithSys");
  for(Int_t icent = 0;icent<8;icent++){
    Float_t y0 = etaWidth[icent]*etaWidth[icent]*h_c1Reg1_Fromfpol3par1_Sim->GetBinContent(icent+1);
    Float_t dy01 = detaWidth[icent]/etaWidth[icent];
    Float_t dy02 = etaWidth[icent]*etaWidth[icent]*h_c1Reg1_Fromfpol3par1_Sim->GetBinError(icent+1);
    Float_t dy0 = TMath::Sqrt(dy01*dy01+dy02*dy02);
    
    // Added new by RR on 22.08.2015
    /*
      Float_t dy01 = detaWidth[icent]/etaWidth[icent];
      Float_t dy02 = h_c1Reg1_Fromfpol3par1_Sim->GetBinError(icent+1)/h_c1Reg1_Fromfpol3par1_Sim->GetBinContent(icent+1);
      Float_t dy0 = y0*TMath::Sqrt(2*dy01*dy01+dy02*dy02);
      //    cout<<" Sys dy01% = "<<1.414*dy01*100<<" dy02% = "<<dy02*100<<endl;
      //
      */ 
    h_y0Reg1_Fromfpol3par1_Sim->SetBinContent(icent+1,y0);
    h_y0Reg1_Fromfpol3par1_Sim->SetBinError(icent+1, dy0);
    cout<<" y0 = "<<y0<<"+/-"<<dy0<<" stat%="<<dy0*100/y0;
    
    y0 = etaWidth[icent]*etaWidth[icent]*h_c1Reg1_Fromfpol3par1_SimWithSys->GetBinContent(icent+1);
    dy0 = etaWidth[icent]*etaWidth[icent]*h_c1Reg1_Fromfpol3par1_SimWithSys->GetBinError(icent+1);
    h_y0Reg1_Fromfpol3par1_SimWithSys->SetBinContent(icent+1,y0);
    h_y0Reg1_Fromfpol3par1_SimWithSys->SetBinError(icent+1,dy0);
    cout<<" sys%="<<dy0*100/y0<<endl;
  }
  
  h_y0Reg1_Fromfpol3par1_Sim->Draw("same");
  h_y0Reg1_Fromfpol3par1_SimWithSys->SetFillStyle(3001);
  h_y0Reg1_Fromfpol3par1_SimWithSys->SetFillColor(6);
  h_y0Reg1_Fromfpol3par1_SimWithSys->Draw("same E3");
  
  TLatex * l = new TLatex(5,0.05,"Data");
  l->SetTextColor(4);
  l->Draw();
  l = new TLatex(5,0.056,"MC Simulation");
  l->SetTextColor(6);
  l->Draw();
}

void DrawFitCoeff1(){
  Float_t y3[5]={0.0018,0.00155,0.0013,0.00105,0.0008};
  Char_t mystr[48];
  TMarker *m;
  TLatex * l;
  
  TCanvas * cpol11 = new TCanvas("cpol11","cpol11",600,600);
  cpol11->cd();
  for(Int_t wreg=0;wreg<5;wreg++){
    h_fpol1_TrkV0_Reg1[wreg]->SetMinimum(0.0005);
    h_fpol1_TrkV0_Reg1[wreg]->SetMaximum(0.004);
    h_fpol1_TrkV0_Reg1[wreg]->SetMarkerStyle(imarker[wreg]);
    h_fpol1_TrkV0_Reg1[wreg]->SetMarkerColor(4);
    h_fpol1_TrkV0_Reg1[wreg]->SetLineColor(4);
    h_fpol1_TrkV0_Reg2[wreg]->SetMarkerStyle(imarker[wreg]);
    h_fpol1_TrkV0_Reg2[wreg]->SetMarkerColor(2);
    h_fpol1_TrkV0_Reg2[wreg]->SetLineColor(2);
    if(wreg==0) h_fpol1_TrkV0_Reg1[wreg]->Draw();
    h_fpol1_TrkV0_Reg1[wreg]->Draw("same");
    h_fpol1_TrkV0_Reg2[wreg]->Draw("same");
    sprintf(mystr,"%3.2f",(wreg+1)*0.04);
    l = new TLatex(17,y3[wreg]-0.00005,mystr);
    l->SetTextSize(0.04);
    l->Draw();
    m = new TMarker(15.5,y3[wreg],imarker[wreg]);
    m->Draw();
  }
  l = new TLatex(4.0,0.0035,"fpol1 = 1.0+c_{1}#eta");
  l->Draw();
  l = new TLatex(4,0.0032,"Reg 1");
  l->SetTextColor(4);
  l->Draw();
  l = new TLatex(4,0.0029,"Reg 2");
  l->SetTextColor(2);
  l->Draw();
  l = new TLatex(15,.00205,"Reg 3 Width");
  l->Draw();
}

void DrawFitCoeff2(){
  Float_t y3[5]={0.0018,0.00155,0.0013,0.00105,0.0008};
  Char_t mystr[48];
  TMarker *m;
  TLatex * l;
  TCanvas * cpol3 = new TCanvas("cpol3","cpol3",600,600);
  cpol3->cd();
  for(Int_t wreg=0;wreg<5;wreg++){
    h_fpol3_TrkV0_Reg1_par1[wreg]->SetMinimum(0.0005);
    h_fpol3_TrkV0_Reg1_par1[wreg]->SetMaximum(0.004);
    h_fpol3_TrkV0_Reg1_par1[wreg]->SetMarkerStyle(imarker[wreg]);
    h_fpol3_TrkV0_Reg1_par1[wreg]->SetLineColor(4);
    h_fpol3_TrkV0_Reg1_par1[wreg]->SetMarkerColor(4);
    h_fpol3_TrkV0_Reg2_par1[wreg]->SetMarkerStyle(imarker[wreg]);
    h_fpol3_TrkV0_Reg2_par1[wreg]->SetLineColor(2);
    h_fpol3_TrkV0_Reg2_par1[wreg]->SetMarkerColor(2);
    if(wreg==0) h_fpol3_TrkV0_Reg1_par1[wreg]->Draw();
    h_fpol3_TrkV0_Reg1_par1[wreg]->Draw("same");
    h_fpol3_TrkV0_Reg2_par1[wreg]->Draw("same");
    sprintf(mystr,"%3.2f",(wreg+1)*0.04);
    l = new TLatex(17,y3[wreg]-0.00005,mystr);
    l->SetTextSize(0.04);
    l->Draw();
    m = new TMarker(15.5,y3[wreg],imarker[wreg]);
    m->Draw();
  }
  l = new TLatex(4.0,0.0035,"fpol3 = 1.0+c_{1}#eta+c_{2}#eta^{3}");
  l->Draw();
  l = new TLatex(4,0.0032,"Reg 1");
  l->SetTextColor(4);
  l->Draw();
  l = new TLatex(4,0.0029,"Reg 2");
  l->SetTextColor(2);
  l->Draw();
  l = new TLatex(15,.00205,"Reg 3 Width");
  l->Draw();
  
  TCanvas * cpol32 = new TCanvas("cpol32","cpol32",600,600);
  cpol32->cd();
  for(Int_t wreg=0;wreg<5;wreg++){
    h_fpol3_TrkV0_Reg1_par2[wreg]->SetMinimum(0.000);
    h_fpol3_TrkV0_Reg1_par2[wreg]->SetMaximum(0.00007);
    h_fpol3_TrkV0_Reg1_par2[wreg]->SetMarkerStyle(imarker[wreg]);
    h_fpol3_TrkV0_Reg1_par2[wreg]->SetLineColor(4);
    h_fpol3_TrkV0_Reg1_par2[wreg]->SetMarkerColor(4);
    h_fpol3_TrkV0_Reg2_par2[wreg]->SetMarkerStyle(imarker[wreg]);
    h_fpol3_TrkV0_Reg2_par2[wreg]->SetLineColor(2);
    h_fpol3_TrkV0_Reg2_par2[wreg]->SetMarkerColor(2);
    if(wreg==0) h_fpol3_TrkV0_Reg1_par2[wreg]->Draw();
    h_fpol3_TrkV0_Reg1_par2[wreg]->Draw("same");
    h_fpol3_TrkV0_Reg2_par2[wreg]->Draw("same");
  }
  l = new TLatex(5.0,0.00015,"fpol3 = 1.0+Par1*#eta+Par2*#eta^{3}");
  l->Draw();
  
  cout<<" Par1 for fpol1"<<endl;

  for(Int_t icent =0;icent<8;icent++){
    Float_t inSlope[5]={0};
    Float_t indSlope[5]={0};
    for(Int_t wreg=0;wreg<5;wreg++){
      inSlope[wreg]=fpol1_TrkV01_Par1[icent][wreg];
      indSlope[wreg]=fpol1_TrkV01_dPar1[icent][wreg];
    }
    GetMeanRMS(5,inSlope,indSlope, Meanfpol1_TrkV0_Reg1_Par1[icent], dMeanfpol1_TrkV0_Reg1_Par1[icent], RMSfpol1_TrkV0_Reg1_Par1[icent]);
    cout<<"CentBin "<<icent<<" "<<Meanfpol1_TrkV0_Reg1_Par1[icent]<<" +/- stat "<<dMeanfpol1_TrkV0_Reg1_Par1[icent]<<" +/- sys "<<RMSfpol1_TrkV0_Reg1_Par1[icent]<<endl;
    Float_t inSlope[5]={0};
    Float_t indSlope[5]={0};
    for(Int_t wreg=0;wreg<5;wreg++){
      inSlope[wreg]=-1.0*fpol1_TrkV02_Par1[icent][wreg];
      indSlope[wreg]=fpol1_TrkV02_dPar1[icent][wreg];
    }
    GetMeanRMS(5,inSlope,indSlope, Meanfpol1_TrkV0_Reg2_Par1[icent], dMeanfpol1_TrkV0_Reg2_Par1[icent], RMSfpol1_TrkV0_Reg2_Par1[icent]);
    cout<<"CentBin "<<icent<<" "<<Meanfpol1_TrkV0_Reg2_Par1[icent]<<" +/- stat "<<dMeanfpol1_TrkV0_Reg2_Par1[icent]<<" +/- sys "<<RMSfpol1_TrkV0_Reg2_Par1[icent]<<endl;
    
  }
  
  
  cout<<" Par1 of fpol3"<<endl;

  for(Int_t icent =0;icent<8;icent++){
    Float_t inSlope[10]={0};
    Float_t indSlope[10]={0};
    for(Int_t wreg=0;wreg<5;wreg++){
      inSlope[wreg]=fpol3_TrkV01_Par1[icent][wreg];
      indSlope[wreg]=fpol3_TrkV01_dPar1[icent][wreg];
      //      inSlope[wreg+5]=-1.0*fpol3_TrkV02_Par1[icent][wreg];
    }
    GetMeanRMS(5,inSlope,indSlope,Meanfpol3_TrkV0_Reg1_Par1[icent],dMeanfpol3_TrkV0_Reg1_Par1[icent],RMSfpol3_TrkV0_Reg1_Par1[icent]);
    cout<<"CentBin "<<icent<<" "<<Meanfpol3_TrkV0_Reg1_Par1[icent]<<" +/-stat "<<dMeanfpol3_TrkV0_Reg1_Par1[icent]<<" +/- sys "<<RMSfpol3_TrkV0_Reg1_Par1[icent]<<endl;
    Float_t inSlope[5]={0};
    Float_t indSlope[5]={0};
    for(Int_t wreg=0;wreg<5;wreg++){
      inSlope[wreg]=-1.0*fpol3_TrkV02_Par1[icent][wreg];
      indSlope[wreg]=fpol3_TrkV02_dPar1[icent][wreg];
      //      inSlope[wreg+5]=-1.0*fpol3_TrkV02_Par1[icent][wreg];
    }
    GetMeanRMS(5,inSlope,indSlope,Meanfpol3_TrkV0_Reg2_Par1[icent],dMeanfpol3_TrkV0_Reg2_Par1[icent],RMSfpol3_TrkV0_Reg2_Par1[icent]);
    cout<<"CentBin "<<icent<<" "<<Meanfpol3_TrkV0_Reg2_Par1[icent]<<" +/-stat "<<dMeanfpol3_TrkV0_Reg2_Par1[icent]<<" +/- sys "<<RMSfpol3_TrkV0_Reg2_Par1[icent]<<endl;
  }
  
  cout<<" Par2 of fpol3"<<endl;
  for(Int_t icent =0;icent<8;icent++){
    Float_t inSlope[5]={0};
    Float_t indSlope[5]={0};
    for(Int_t wreg=0;wreg<5;wreg++){
      inSlope[wreg]=fpol3_TrkV01_Par2[icent][wreg];
      indSlope[wreg]=fpol3_TrkV01_dPar2[icent][wreg];
      //      inSlope[wreg+5]=-1.0*fpol3_TrkV02_Par2[icent][wreg];
    }
    GetMeanRMS(5,inSlope, indSlope,Meanfpol3_TrkV0_Reg1_Par2[icent],dMeanfpol3_TrkV0_Reg1_Par2[icent], RMSfpol3_TrkV0_Reg1_Par2[icent]);
    cout<<"CentBin "<<icent<<" "<<Meanfpol3_TrkV0_Reg1_Par2[icent]<<" +/- stat "<<dMeanfpol3_TrkV0_Reg1_Par2[icent]<<" +/- sys "<<RMSfpol3_TrkV0_Reg1_Par2[icent]<<endl;
  }
  for(Int_t icent =0;icent<8;icent++){
    Float_t inSlope[5]={0};
    Float_t indSlope[5]={0};
    for(Int_t wreg=0;wreg<5;wreg++){
      inSlope[wreg]=fpol3_TrkV02_Par2[icent][wreg];
      indSlope[wreg]=fpol3_TrkV02_dPar2[icent][wreg];
      //      inSlope[wreg+5]=-1.0*fpol3_TrkV02_Par2[icent][wreg];
    }
    GetMeanRMS(5,inSlope, indSlope,Meanfpol3_TrkV0_Reg2_Par2[icent],dMeanfpol3_TrkV0_Reg2_Par2[icent], RMSfpol3_TrkV0_Reg2_Par2[icent]);
    cout<<"CentBin "<<icent<<" "<<Meanfpol3_TrkV0_Reg2_Par2[icent]<<" +/- stat "<<dMeanfpol3_TrkV0_Reg2_Par2[icent]<<" +/- sys "<<RMSfpol3_TrkV0_Reg2_Par2[icent]<<endl;
  }
}


void ObtainRapShiftMS(){
  
  TH1F * hShiftGlauber = new TH1F("hShiftGlauber","hShiftGlauber",8.0,0.0,40.0);
  TH1F * hSlopeMS1 = new TH1F("hSlopeMS1","hSlopeMS1",8,0.0,40.0);
  TH1F * hSlopeMS2 = new TH1F("hSlopeMS2","hSlopeMS2",8,0.0,40.0);
  TH1F * hSlopefpol1MS1 = new TH1F("hSlopefpol1MS1","hSlopefpol1MS1",8,0.0,40.0);
  TH1F * hSlopefpol1MS2 = new TH1F("hSlopefpol1MS2","hSlopefpol1MS2",8,0.0,40.0);

  
  for(Int_t icent =0;icent<8;icent++){
    Float_t slope,dslope;
    cout<<fpol1_TrkV01_Par1[icent]<<" "<<fpol1_TrkV02_Par1[icent]<<endl;

    slope = fpol1_TrkV01_Par1[icent]*etasigma*etasigma;
    dslope = fpol1_TrkV01_dPar1[icent]*etasigma*etasigma;
    hSlopefpol1MS1->SetBinContent(icent+1,slope);
    hSlopefpol1MS1->SetBinError(icent+1,dslope);
    
    slope = TMath::Abs(fpol1_TrkV02_Par1[icent])*etasigma*etasigma;
    dslope = TMath::Abs(fpol1_TrkV02_dPar1[icent])*etasigma*etasigma;
    hSlopefpol1MS2->SetBinContent(icent+1,slope);
    hSlopefpol1MS2->SetBinError(icent+1,dslope);


    slope = fpol3_TrkV01_Par1[icent]*etasigma*etasigma;
    dslope = fpol3_TrkV01_dPar1[icent]*etasigma*etasigma;
    hSlopeMS1->SetBinContent(icent+1,slope);
    hSlopeMS1->SetBinError(icent+1,dslope);
    
    slope = TMath::Abs(fpol3_TrkV02_Par1[icent])*etasigma*etasigma;
    dslope = TMath::Abs(fpol3_TrkV02_dPar1[icent])*etasigma*etasigma;
    hSlopeMS2->SetBinContent(icent+1,slope);
    hSlopeMS2->SetBinError(icent+1,dslope);

    hShiftGlauber->SetBinContent(icent+1,SlopeGl[icent]);
    hShiftGlauber->SetBinError(icent+1,dSlopeGl[icent]);
  }
  TH1F * hy0FromGlauber = new TH1F("hy0FromGlauber","hy0FromGlauber",8,0.0,40.0);
  for(Int_t icent = 0;icent<8;icent++){
    hy0FromGlauber->SetBinContent(icent+1,y0FromG[icent]);
    hy0FromGlauber->SetBinError(icent+1,dy0FromG[icent]);
  }
  
  hShiftGlauber->SetXTitle("Centrality\%");
  hShiftGlauber->SetYTitle("y_{0} Shift ");
  hShiftGlauber->SetMarkerStyle(4);
  hShiftGlauber->SetMarkerColor(6);
  hShiftGlauber->SetLineColor(6);
  hShiftGlauber->SetMinimum(0.0);
  hShiftGlauber->SetMaximum(0.10);

  hSlopeMS1->SetXTitle("Centrality\%");
  hSlopeMS1->SetYTitle("y_{0} Shift from dN/d#eta");
  hSlopeMS1->SetMarkerStyle(8);
  hSlopeMS1->SetMarkerColor(4);
  hSlopeMS1->SetLineColor(4);
  hSlopeMS1->SetMinimum(0.0);

  hSlopeMS2->SetMarkerStyle(8);
  hSlopeMS2->SetMarkerColor(2);
  hSlopeMS2->SetLineColor(2);
  
 
  hSlopefpol1MS1->SetMarkerStyle(29);
  hSlopefpol1MS1->SetMarkerColor(4);
  hSlopefpol1MS1->SetLineColor(4);
  hSlopefpol1MS2->SetMarkerStyle(29);
  hSlopefpol1MS2->SetMarkerColor(2);
  hSlopefpol1MS2->SetLineColor(2);

  hy0FromGlauber->SetMarkerStyle(8);
  hy0FromGlauber->SetMarkerColor(6);
  

  
  TCanvas * cSlopeMS = new TCanvas("cSlopeMS","cSlopeMS",600,600);

  hShiftGlauber->Draw();
  hSlopefpol1MS1->Draw("same");
  hSlopefpol1MS2->Draw("same");
  //hSlopeMS1->Draw("same");
  //  hSlopeMS2->Draw("same");
  hy0FromGlauber->Draw("same");

  TMarker * m = new TMarker(3.0,0.091,8);
  m->SetMarkerColor(6);
  m->Draw();
  TLatex * l =  new TLatex(5,0.09,"GlauberModel");
  l->SetTextColor(6);
  l->SetTextSize(0.03);
  l->Draw();

  m = new TMarker(3.0,0.084,29);
  m->SetMarkerColor(4);
  m->Draw();
  //  l = new TLatex(5.0,0.083,"dN/d#eta slope (2 Par) Reg1");
  //  l = new TLatex(5.0,0.083,"Fit 1+c_{1}#eta+c_{3}#eta^{3} to dN/d#eta Reg1");
    l = new TLatex(4.0,0.083,"Fit 1+c_{1}#eta to dN/d#eta Reg1");
  l->SetTextColor(4);
  l->SetTextSize(0.03);
  l->Draw();
  
  m = new TMarker(3.0,0.077,29);
  m->SetMarkerColor(2);
  m->Draw();
  //  l = new TLatex(5.0,0.076,"Fit 1+c_{1}#eta+c_{3}#eta^{3} to dN/d#eta Reg2");
  l = new TLatex(4.0,0.076,"Fit 1+c_{1}#eta to dN/d#eta Reg2");
  //  l = new TLatex(5.0,0.076,"dN/d#eta slope (2 Par) Reg2");
  l->SetTextColor(2);
  l->SetTextSize(0.03);
  l->Draw();

  /*
  m = new TMarker(3.0,0.070,29);
  m->SetMarkerColor(4);
  m->Draw();
  l = new TLatex(5.0,0.069,"dN/d#eta slope (1 Par) Reg1");
  l->SetTextColor(4);
  l->SetTextSize(0.03);
  l->Draw();
  m = new TMarker(3.0,0.063,29);
  m->SetMarkerColor(2);
  m->Draw();
  l = new TLatex(5.0,0.062,"dN/d#eta slope (1 Par) Reg2");
  l->SetTextColor(2);
  l->SetTextSize(0.03);
  l->Draw();
  */

}

void ObtainRapShiftFromGlauber(){
  TH1F * hShiftgl[8];
  Char_t hname[64];
  TFile * fg[8];
  TNtuple * nt;
  Float_t Npart,Ncoll,NpA,NnA,NpB,NnB,SpA,SnA,SpB,SnB;
  Float_t B;

  
  for(Int_t icent = 0;icent<8;icent++){
    sprintf(hname,"hShiftGlauber_Cent%ito%i",icent*5,(icent+1)*5);
    hShiftgl[icent] = new TH1F(hname,hname,100,0,0.25);
    sprintf(fname,"/Users/rashmiraniwala/alice/tglaubermc_2.1/Glauber%ito%i.root",icent*5,(icent+1)*5);
    //    cout<<" Reading file "<<fname<<endl;
    fg[icent] = new TFile(fname);
    fg[icent]->GetObject("nt",nt);
    nt->SetBranchAddress("Npart",&Npart);
    nt->SetBranchAddress("Ncoll",&Ncoll);
    nt->SetBranchAddress("B",&B);
    nt->SetBranchAddress("NpA",&NpA);
    nt->SetBranchAddress("NpB",&NpB);
    nt->SetBranchAddress("NnA",&NnA);
    nt->SetBranchAddress("NnB",&NnB);
    nt->SetBranchAddress("SpA",&SpA);
    nt->SetBranchAddress("SpB",&SpB);
    nt->SetBranchAddress("SnA",&SnA);
    nt->SetBranchAddress("SnB",&SnB);
    //    cout<<" Entries in this file are "<<nt->GetEntries()<<endl;
    //    for(Int_t iev=0;iev<nt->GetEntries();iev++){
    for(Int_t iev=0;iev<100000;iev++){
      nt->GetEntry(iev);
      Float_t shift = 0.5*TMath::Log((NnA+NpA)/(NnB+NpB));
      hShiftgl[icent]->Fill(shift);
    }
    cout<<"icent bin  = "<<icent+1<<" Mean 0.5log(A/B) = "<<hShiftgl[icent]->GetMean()<<" +/- "<<hShiftgl[icent]->GetMeanError()<<endl;
    SlopeGl[icent]=hShiftgl[icent]->GetMean();
    dSlopeGl[icent]=hShiftgl[icent]->GetMeanError();
  }
}

void GetMeanRMS(Int_t imax, Float_t * val, Float_t *dval, Float_t& mean, Float_t& dmeanstat, Float_t& dmeansys){

  mean=0;
  dmeansys=0;
  dmeanstat = 0;
  for(Int_t i=0;i<imax;i++){
    //    cout<<val[i]<<" ";
    mean+=val[i];
    dmeansys+=val[i]*val[i];
    dmeanstat+=dval[i]*dval[i];
  }
  mean = mean/(imax);
  dmeansys = dmeansys/imax;
  dmeansys = TMath::Sqrt(dmeansys-mean*mean);
  dmeanstat = (1.0/imax)*TMath::Sqrt(dmeanstat);
  //  cout<<"="<<mean<<"+/-"<<dmean<<endl;

  return;
}

void ReadData(){
  for(Int_t icent = 0;icent<8;icent++){
    for(Int_t wreg3 = 0;wreg3<5;wreg3++){
      //      sprintf(fname,"/Users/rashmiraniwala/alice/InitFluct/ANALYSIS/Version4/NarrowAlphaBins/TrkCentrality/dNdetaRatioFitData_Cent%3.3dto%3.3d_V0MTrk_TrkV0MCentrality_Reg3_%2.2f.dat",icent*5,(icent+1)*5,(wreg3+1)*0.02);
      sprintf(fname,"/Users/rashmiraniwala/ToALICE/ANALYSIS/dNdEtaRatios/V0MTrkCentrality/dNdetaRatioFitData_Cent%3.3dto%3.3d_V0MTrk_TrkV0MCentrality_Reg3_%2.2f.dat",icent*5,(icent+1)*5,(wreg3+1)*0.02);
      // sprintf(fname,"dNdetaRatioFitData_Cent%3.3dto%3.3d_V0M_V0MCentrality_Reg3_%2.2f.dat",icent*5,(icent+1)*5,(wreg3+1)*0.02);
      //    sprintf(fname,"dNdetaRatioFitData_Cent%3.3dto%3.3d.dat",icent*5,(icent+1)*5);
      cout<<"Reading file "<<fname<<endl;
      fdata[icent][wreg3].open(fname,ios::in);
      fdata[icent][wreg3]>>Reg12Limit>>Reg3Limit;
      fdata[icent][wreg3]>>Reg1MinBin>>Reg1MaxBin;
      fdata[icent][wreg3]>>Reg2MinBin>>Reg2MaxBin;
      fdata[icent][wreg3]>>Reg3MinBin>>Reg3MaxBin;
      fdata[icent][wreg3]>>nevent1>>nevent2>>nevent3;
      
      fdata[icent][wreg3]>>fpol1_TrkV01_Par1[icent][wreg3]>>fpol1_TrkV01_dPar1[icent][wreg3];
      fdata[icent][wreg3]>>fpol1_TrkV02_Par1[icent][wreg3]>>fpol1_TrkV02_dPar1[icent][wreg3];
      fdata[icent][wreg3]>>fpol3_TrkV01_Par1[icent][wreg3]>>fpol3_TrkV01_dPar1[icent][wreg3]>>fpol3_TrkV01_Par2[icent][wreg3]>>fpol3_TrkV01_dPar2[icent][wreg3];
      fdata[icent][wreg3]>>fpol3_TrkV02_Par1[icent][wreg3]>>fpol3_TrkV02_dPar1[icent][wreg3]>>fpol3_TrkV02_Par2[icent][wreg3]>>fpol3_TrkV02_dPar2[icent][wreg3];
      fdata[icent][wreg3]>>fpol1_V01_Par1[icent][wreg3]>>fpol1_V01_dPar1[icent][wreg3];
      fdata[icent][wreg3]>>fpol1_V02_Par1[icent][wreg3]>>fpol1_V02_dPar1[icent][wreg3];
      fdata[icent][wreg3]>>fpol3_V01_Par1[icent][wreg3]>>fpol3_V01_dPar1[icent][wreg3]>>fpol3_V01_Par2[icent][wreg3]>>fpol3_V01_dPar2[icent][wreg3];
      fdata[icent][wreg3]>>fpol3_V02_Par1[icent][wreg3]>>fpol3_V02_dPar1[icent][wreg3]>>fpol3_V02_Par2[icent][wreg3]>>fpol3_V02_dPar2[icent][wreg3];
      fdata[icent][wreg3]>>fpol1_Trk1_Par1[icent][wreg3]>>fpol1_Trk1_dPar1[icent][wreg3];
      fdata[icent][wreg3]>>fpol1_Trk2_Par1[icent][wreg3]>>fpol1_Trk2_dPar1[icent][wreg3];
      fdata[icent][wreg3]>>fpol3_Trk1_Par1[icent][wreg3]>>fpol3_Trk1_dPar1[icent][wreg3]>>fpol3_Trk1_Par2[icent][wreg3]>>fpol3_Trk1_dPar2[icent][wreg3];
      fdata[icent][wreg3]>>fpol3_Trk2_Par1[icent][wreg3]>>fpol3_Trk2_dPar1[icent][wreg3]>>fpol3_Trk2_Par2[icent][wreg3]>>fpol3_Trk2_dPar2[icent][wreg3];
      
      cout<<fpol1_TrkV01_Par1[icent][wreg3]<<" "<<fpol1_TrkV01_dPar1[icent][wreg3]<<endl;
      cout<<fpol1_TrkV02_Par1[icent][wreg3]<<" "<<fpol1_TrkV02_dPar1[icent][wreg3]<<endl;
      cout<<fpol3_TrkV01_Par1[icent][wreg3]<<" "<<fpol3_TrkV01_dPar1[icent][wreg3]<<" "<<fpol3_TrkV01_Par2[icent][wreg3]<<" "<<fpol3_TrkV01_dPar2[icent][wreg3]<<endl;
      cout<<fpol3_TrkV02_Par1[icent][wreg3]<<" "<<fpol3_TrkV02_dPar1[icent][wreg3]<<" "<<fpol3_TrkV02_Par2[icent][wreg3]<<" "<<fpol3_TrkV02_dPar2[icent][wreg3]<<endl;
      cout<<fpol1_V01_Par1[icent][wreg3]<<" "<<fpol1_V01_dPar1[icent][wreg3]<<endl;
      cout<<fpol1_V02_Par1[icent][wreg3]<<" "<<fpol1_V02_dPar1[icent][wreg3]<<endl;
      cout<<fpol3_V01_Par1[icent][wreg3]<<" "<<fpol3_V01_dPar1[icent][wreg3]<<" "<<fpol3_V01_Par2[icent][wreg3]<<" "<<fpol3_V01_dPar2[icent][wreg3]<<endl;
      cout<<fpol3_V02_Par1[icent][wreg3]<<" "<<fpol3_V02_dPar1[icent][wreg3]<<" "<<fpol3_V02_Par2[icent][wreg3]<<" "<<fpol3_V02_dPar2[icent][wreg3]<<endl;
      cout<<fpol1_Trk1_Par1[icent][wreg3]<<" "<<fpol1_Trk1_dPar1[icent][wreg3]<<endl;
      cout<<fpol1_Trk2_Par1[icent][wreg3]<<" "<<fpol1_Trk2_dPar1[icent][wreg3]<<endl;
      cout<<fpol3_Trk1_Par1[icent][wreg3]<<" "<<fpol3_Trk1_dPar1[icent][wreg3]<<" "<<fpol3_Trk1_Par2[icent][wreg3]<<" "<<fpol3_Trk1_dPar2[icent][wreg3]<<endl;
      cout<<fpol3_Trk2_Par1[icent][wreg3]<<" "<<fpol3_Trk2_dPar1[icent][wreg3]<<" "<<fpol3_Trk2_Par2[icent][wreg3]<<" "<<fpol3_Trk2_dPar2[icent][wreg3]<<endl;
    }
  }
}

void FillHistosData(){
  Char_t hname[128];
  for(Int_t wreg=0;wreg<5;wreg++){
    sprintf(hname,"h_fpol1_TrkV0_Reg1_wreg3%d",wreg+1);
    h_fpol1_TrkV0_Reg1[wreg]= new TH1F(hname,hname,8.0,0.0,40.0);
    sprintf(hname,"h_fpol3_TrkV0_Reg1_par1_wreg3%d",wreg+1);
    h_fpol3_TrkV0_Reg1_par1[wreg]= new TH1F(hname,hname,8.0,0.0,40.0);
    sprintf(hname,"h_fpol1_TrkV0_Reg1_par2_wreg3%d",wreg+1);
    h_fpol3_TrkV0_Reg1_par2[wreg]= new TH1F(hname,hname,8.0,0.0,40.0);
    sprintf(hname,"h_fpol1_TrkV0_Reg2_wreg3%d",wreg+1);
    h_fpol1_TrkV0_Reg2[wreg]= new TH1F(hname,hname,8.0,0.0,40.0);
    sprintf(hname,"h_fpol3_TrkV0_Reg2_par1_wreg3%d",wreg+1);
    h_fpol3_TrkV0_Reg2_par1[wreg]= new TH1F(hname,hname,8.0,0.0,40.0);
    sprintf(hname,"h_fpol3_TrkV0_Reg2_par2_wreg3%d",wreg+1);
    h_fpol3_TrkV0_Reg2_par2[wreg]= new TH1F(hname,hname,8.0,0.0,40.0);
    
    h_fpol1_TrkV0_Reg1[wreg]->SetXTitle("Centrality");
    h_fpol1_TrkV0_Reg1[wreg]->SetYTitle("c_{1}");
    h_fpol3_TrkV0_Reg1_par1[wreg]->SetXTitle("Centrality");
    h_fpol3_TrkV0_Reg1_par1[wreg]->SetYTitle("c_{1}");
    h_fpol3_TrkV0_Reg1_par2[wreg]->SetXTitle("Centrality");
    h_fpol3_TrkV0_Reg1_par2[wreg]->SetYTitle("c_{2}");
  }
  
  
  for(Int_t icent = 0;icent<8;icent++){
    for(Int_t wreg=0;wreg<5;wreg++){
      h_fpol1_TrkV0_Reg1[wreg]->SetBinContent(icent+1,fpol1_TrkV01_Par1[icent][wreg]);
      h_fpol1_TrkV0_Reg1[wreg]->SetBinError(icent+1,fpol1_TrkV01_dPar1[icent][wreg]);
      h_fpol1_TrkV0_Reg2[wreg]->SetBinContent(icent+1,TMath::Abs(fpol1_TrkV02_Par1[icent][wreg]));
      h_fpol1_TrkV0_Reg2[wreg]->SetBinError(icent+1,TMath::Abs(fpol1_TrkV02_dPar1[icent][wreg]));
      
      h_fpol3_TrkV0_Reg1_par1[wreg]->SetBinContent(icent+1,fpol3_TrkV01_Par1[icent][wreg]);
      h_fpol3_TrkV0_Reg1_par1[wreg]->SetBinError(icent+1,fpol3_TrkV01_dPar1[icent][wreg]);
      h_fpol3_TrkV0_Reg2_par1[wreg]->SetBinContent(icent+1,TMath::Abs(fpol3_TrkV02_Par1[icent][wreg]));
      h_fpol3_TrkV0_Reg2_par1[wreg]->SetBinError(icent+1,TMath::Abs(fpol3_TrkV02_dPar1[icent][wreg]));
      
      h_fpol3_TrkV0_Reg1_par2[wreg]->SetBinContent(icent+1,fpol3_TrkV01_Par2[icent][wreg]);
      h_fpol3_TrkV0_Reg1_par2[wreg]->SetBinError(icent+1,fpol3_TrkV01_dPar2[icent][wreg]);
      h_fpol3_TrkV0_Reg2_par2[wreg]->SetBinContent(icent+1,TMath::Abs(fpol3_TrkV02_Par2[icent][wreg]));
      h_fpol3_TrkV0_Reg2_par2[wreg]->SetBinError(icent+1,TMath::Abs(fpol3_TrkV02_dPar2[icent][wreg]));
    }  
  }
}
