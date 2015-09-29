#include "/Users/rashmiraniwala/ToALICE/ANALYSIS/SetStyle.C"
Int_t centmin,centmax;
TFile *fin1, *fin2, *fout, *fin;
TH1F * hSum;
TH1F * hZDCAsym, *hZDCAsym1, *hZDCAsym2, *hZDCAsym3, *hZDCAsymAbs;
TH1F * hVzSum1, *hVzSum3, *hVzSum2;
TH1F * hVzSum, *hCentSum;
TH1F * hCentSum1, *hCentSum3, *hCentSum2;

TH1F * hTrkratio1, *hTrkratio2, *hTrkratio3;
TH1F * hV0Aratio1, *hV0Aratio2, *hV0Aratio3;
TH1F * hV0Cratio1, *hV0Cratio2, *hV0Cratio3;
TH1F *hV0ratio1, *hV0ratio2;
TH1F * hV0AOnlyRatio, *hV0COnlyRatio;

TH1F *hV0Ratio1,*hV0Ratio2, *hV0ACTrk1by3,*hV0ACTrk2by3, * hV0CReflect1, *hV0CREflect2;

TProfile * hpSum;
TProfile * hpTrkReg1, *hpTrkReg2, *hpTrkReg3;
TProfile * hpV0AReg1, *hpV0AReg2, *hpV0AReg3;
TProfile * hpV0CReg1, *hpV0CReg2, *hpV0CReg3;
TProfile *hpV0CReg1and2, *hpV0AReg1and2, *hpTrkReg1and2;
Float_t V0GlobalCorrReg1=1.0, V0GlobalCorrReg2=1.0;
Float_t TrkGlobalCorrReg1=1.0, TrkGlobalCorrReg2=1.0;

Float_t V0CEta[4]={-3.45,-2.95,-2.45,-1.95};
Float_t dV0CEta[4]={0.25,0.25,0.25,0.25};
Float_t V0AEta[4]={3.1,3.65,4.2,4.8};
Float_t dV0AEta[4]={0.3,0.25,0.3,0.3};
Float_t V0CEtaReflect[4]={3.45,2.95,2.45,1.95};

Float_t V0Aratio1[4], V0Aratio2[4], V0Cratio1[4], V0Cratio2[4];
Float_t dV0Aratio1[4], dV0Aratio2[4], dV0Cratio1[4], dV0Cratio2[4];
Float_t Trkratio1[20], Trkratio2[20], TrkEta[20];
Float_t dTrkratio1[20], dTrkratio2[20], dTrkEta[20];

Float_t zdcasym;
TF1 * fpol1, *fpol3, *fpol5;
Int_t OptDrawVzCentZDC = 0;
Char_t fname1[256],fname2[256],fnameroot[256], fnamedat[256];
ofstream fdata;
Int_t Reg1MinBin, Reg1MaxBin;
Int_t Reg2MinBin, Reg2MaxBin;
Int_t Reg3MinBin, Reg3MaxBin;
Double_t Reg12Limit = 0.1, Reg3Limit = 0.1;

void GetEtaRatioV0Trk(Int_t incentmin = 15, Int_t incentmax = 20, Double_t inReg3Limit = 0.1){

  Reg3Limit = inReg3Limit;
  centmin = incentmin;
  centmax = incentmax;
  SetStyle();
  gStyle->SetOptFit(0);

  // Input Files
  // Track points are taken from fin1 and V0M points are taken from fin2

  //sprintf(fname1,"/Users/rashmiraniwala/ToALICE/ANALYSIS/SystematicsVertexZ/V0MCentrality/AsymWithVzCentWeightsCent%3.3ito%3.3i.root",centmin,centmax);
  sprintf(fname1,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/AsymWithVzCentWeightsCent%3.3ito%3.3i.root",centmin,centmax);
  cout<<" Obtaining file for dNdeta of Tracks "<<fname1<<endl;
  fin1 = new TFile(fname1);
  
  //sprintf(fname2,"/Users/rashmiraniwala/ToALICE/ANALYSIS/SystematicsVertexZ/TrkCentrality/AsymWithVzCentWeightsCent%3.3ito%3.3i.root",centmin,centmax);
  sprintf(fname2,"/Users/rashmiraniwala/ToALICE/ANALYSIS/TrkCentrality/AsymWithVzCentWeightsCent%3.3ito%3.3i.root",centmin,centmax);
  cout<<" Obtaining file for dNdeta of V0M "<<fname2<<endl;
  fin2 = new TFile(fname2);
  
  // Output root file
  //sprintf(fnameroot,"SystematicsVertexZ/Ratio4Cent%3.3dto%3.3d_V0MTrk_TrkV0MCentrality_Reg3_%2.2f.root",centmin,centmax,Reg3Limit);
  sprintf(fnameroot,"V0MTrkCentrality/Ratio4Cent%3.3dto%3.3d_V0MTrk_TrkV0MCentrality_Reg3_%2.2f.root",centmin,centmax,Reg3Limit);
  cout<<" Creating a new output File for histograms"<<fnameroot<<endl;
  fout  = new TFile(fnameroot,"RECREATE");
  // Output dat file
  //  sprintf(fnamedat,"SystematicsVertexZ/dNdetaRatioFitData_Cent%3.3dto%3.3d_V0MTrk_TrkV0MCentrality_Reg3_%2.2f.dat",centmin,centmax,Reg3Limit);
  sprintf(fnamedat,"V0MTrkCentrality/dNdetaRatioFitData_Cent%3.3dto%3.3d_V0MTrk_TrkV0MCentrality_Reg3_%2.2f.dat",centmin,centmax,Reg3Limit);
  fdata.open(fnamedat,ios::out);


  Int_t NReg3Bins = Reg3Limit/0.02;
  //  cout<<"NReg3Bins = "<<2*NReg3Bins<<endl;
  Reg3MinBin = 50-NReg3Bins+1;
  Reg3MaxBin = 50+NReg3Bins;
  cout<<"Reg3MinBin="<<Reg3MinBin<<" MaxBin = "<<Reg3MaxBin<<endl;
  
  Int_t NReg1Bins = -1.0*Reg12Limit/0.02;
  //  cout<<"NReg1Bins = "<<NReg1Bins<<endl;
  Int_t NReg2Bins = Reg12Limit/0.02;
  //cout<<"NReg2Bins = "<<NReg2Bins<<endl;
  Reg1MinBin = 1; Reg1MaxBin =50+NReg1Bins;
  Reg2MinBin = 50+NReg2Bins+1; Reg2MaxBin =100;
  cout<<"Reg1MinBin = "<<Reg1MinBin<<" Reg1MaxBin = "<<Reg1MaxBin<<endl;
  cout<<"Reg2MinBin = "<<Reg2MinBin<<" Reg2MaxBin = "<<Reg2MaxBin<<endl;
  
  
  FillZDCHistoFromFile(fin1);

  FillVzFromFile(fin1);
  FillCentFromFile(fin1);
  fdata<<Reg12Limit<<" "<<Reg3Limit<<endl;
  fdata<<Reg1MinBin<<" "<<Reg1MaxBin<<endl;
  fdata<<Reg2MinBin<<" "<<Reg2MaxBin<<endl;
  fdata<<Reg3MinBin<<" "<<Reg3MaxBin<<endl;
  fdata<<hVzSum1->GetEntries()<<" "<<hVzSum2->GetEntries()<<" "<<hVzSum3->GetEntries()<<" "<<endl;

  FillTrkFromFile(fin1);
  FillV0AFromFile(fin2);
  FillV0CFromFile(fin2);

  fpol1 = new TF1("fpol1","1.0+[0]*x",-4.0,5.0);
  fpol3 = new TF1("fpol3","1.0+[0]*x+[1]*x*x*x",-4.0,5.0);
  //  fpol5 = new TF1("fpol5","1.0+[0]*x+[1]*x*x*x+[2]*x*x*x*x*x",-4.0,5.0);

  GetGlobalCorrection2V0();

  ObtainRatios();
  //DrawTrkRatio("V0MCentrality");
  // Give FitOption as "Trk" or "V0A+V0C" or "Trk+V0A+V0C"); 
  //DrawTrkV0Ratio("V0MCentrality","V0MCentrality","Trk","fpol1");
  // DrawTrkV0Ratio("V0MCentrality","TrkCentrality","Trk","fpol1");
  //DrawTrkV0Ratio("V0MCentrality","TrkCentrality","V0A+V0C","fpol1");
  DrawTrkV0Ratio("V0MCentrality","TrkCentrality","Trk+V0A+V0C","fpol1");
  //DrawTrkV0Ratio("V0MCentrality","TrkCentrality","Trk+V0A+V0C","fpol3");
  //   DrawTrkV0Ratio("TrkCentrality","TrkCentrality","Trk+V0A+V0C","fpol3");
  //  DrawReg1and2Plots("V0MCentrality","TrkCentrality");
  //WriteHistograms();
  cout<<"Events in Reg 1 2 3 are "<<hVzSum1->GetEntries()<<" "<<hVzSum2->GetEntries()<<" "<<hVzSum3->GetEntries()<<" "<<endl;

}

void ObtainRatios(){
  hV0Aratio1 = (TH1F*)((TProfile*)hpV0AReg1->Clone("hV0Aratio1"))->ProjectionX();
  hV0Aratio2 = (TH1F*)((TProfile*)hpV0AReg2->Clone("hV0Aratio2"))->ProjectionX();
  hV0Aratio3 = (TH1F*)((TProfile*)hpV0AReg3->Clone("hV0Aratio3"))->ProjectionX();
  //  cout<<" V0A Sum Wts for Reg 1,2,3 are "<<hV0Aratio1->GetSumOfWeights()<<" "<<hV0Aratio2->GetSumOfWeights()<<" "<<hV0Aratio3->GetSumOfWeights()<<endl;
  hV0Aratio1->Scale(V0GlobalCorrReg1);
  hV0Aratio2->Scale(V0GlobalCorrReg2);
  hV0Aratio1->Divide(hV0Aratio3);
  hV0Aratio2->Divide(hV0Aratio3);
  
  hV0Cratio1 = (TH1F*)((TProfile*)hpV0CReg1->Clone("hV0Cratio1"))->ProjectionX();
  hV0Cratio2 = (TH1F*)((TProfile*)hpV0CReg2->Clone("hV0Cratio2"))->ProjectionX();
  hV0Cratio3 = (TH1F*)((TProfile*)hpV0CReg3->Clone("hV0Cratio3"))->ProjectionX();
  //    cout<<" V0C Sum Wts for Reg 1,2,3 are "<<hV0Cratio1->GetSumOfWeights()<<" "<<hV0Cratio2->GetSumOfWeights()<<" "<<hV0Cratio3->GetSumOfWeights()<<endl;

  hV0Cratio1->Scale(V0GlobalCorrReg1);
  hV0Cratio2->Scale(V0GlobalCorrReg2);
  hV0Cratio1->Divide(hV0Cratio3);
  hV0Cratio2->Divide(hV0Cratio3);
  
  hTrkratio1 = (TH1F*)((TProfile*)hpTrkReg1->Clone("hTrkratio1"))->ProjectionX();
  hTrkratio2 = (TH1F*)((TProfile*)hpTrkReg2->Clone("hTrkratio2"))->ProjectionX();
  hTrkratio3 = (TH1F*)((TProfile*)hpTrkReg3->Clone("hTrkratio3"))->ProjectionX();
  //  cout<<" TrkSum Wts for Reg 1,2,3 are "<<hTrkratio1->GetSumOfWeights()<<" "<<hTrkratio2->GetSumOfWeights()<<" "<<hTrkratio3->GetSumOfWeights()<<endl;
  hTrkratio1->Scale(TrkGlobalCorrReg1);
  hTrkratio2->Scale(TrkGlobalCorrReg2);
  hTrkratio1->Divide(hTrkratio3);
  hTrkratio2->Divide(hTrkratio3);
  
  //  cout<<" Ratios && error for V0A different eta"<<endl;
  for(Int_t iv = 0;iv<4;iv++){
    V0Aratio1[iv] = hV0Aratio1->GetBinContent(iv+1);
    V0Aratio2[iv] = hV0Aratio2->GetBinContent(iv+1);
    dV0Aratio1[iv] = hV0Aratio1->GetBinError(iv+1);
    dV0Aratio2[iv] = hV0Aratio2->GetBinError(iv+1);
    //    cout<<V0AEta[iv]<<" "<<V0Aratio1[iv]<<" "<<dV0Aratio1[iv]<<" "<<V0Aratio2[iv]<<" "<<dV0Aratio2[iv]<<endl;
  }
  //  cout<<" Ratios && error for V0C different eta"<<endl;
  for(Int_t iv = 0;iv<4;iv++){
    V0Cratio1[iv] = hV0Cratio1->GetBinContent(iv+1);
    V0Cratio2[iv] = hV0Cratio2->GetBinContent(iv+1);
    dV0Cratio1[iv] = hV0Cratio1->GetBinError(iv+1);
    dV0Cratio2[iv] = hV0Cratio2->GetBinError(iv+1);
    //    cout<<V0CEta[iv]<<" "<<V0Cratio1[iv]<<" "<<dV0Cratio1[iv]<<" "<<V0Cratio2[iv]<<" "<<dV0Cratio2[iv]<<endl;
    
  }
  
  //  cout<<" Ratios and error for Trks different eta"<<endl;
  Int_t netabins = hTrkratio1->GetXaxis()->GetNbins();
  Int_t nbins = 0;
  for(Int_t ieta = 1;ieta<netabins-1;ieta++){
    Trkratio1[nbins] = hTrkratio1->GetBinContent(ieta+1);
    Trkratio2[nbins] = hTrkratio2->GetBinContent(ieta+1);
    dTrkratio1[nbins] = hTrkratio1->GetBinError(ieta+1);
    dTrkratio2[nbins] = hTrkratio2->GetBinError(ieta+1);
    TrkEta[nbins] = hTrkratio1->GetBinCenter(ieta+1);
    dTrkEta[nbins] = 0.5*hTrkratio1->GetBinWidth(ieta+1);
    //    cout<<TrkEta[nbins]<<" "<<Trkratio1[nbins]<<" "<<Trkratio2[nbins]<<" "<<dTrkratio1[nbins]<<" "<<dTrkratio2[nbins]<<endl;
    //    cout<<" ieta = "<<ieta<<" eta center "<<TrkEta[nbins]<<endl;
    nbins++;
    //    cout<<"nBins  = "<<nbins<<endl;
  }
  return;
  
}

void WriteHistograms(){
  
  fout->cd();
  
  hpTrkReg1->Write();
  hpTrkReg2->Write();
  hpTrkReg3->Write();

  hpV0AReg1->Write();
  hpV0AReg2->Write();
  hpV0AReg3->Write();
     
  hpV0CReg1->Write();
  hpV0CReg2->Write();
  hpV0CReg3->Write();

  hTrkratio1->Write();
  hTrkratio2->Write();

  hV0Aratio1->Write();
  hV0Aratio2->Write();
  hV0Cratio1->Write();
  hV0Cratio2->Write();

  hV0Ratio1->Write();
  hV0Ratio2->Write();

  hV0ACTrk1by3->Write();
  hV0ACTrk2by3->Write();

  hVzSum1->Write();
  hVzSum2->Write();
  hVzSum3->Write();

  hCentSum1->Write();
  hCentSum2->Write();
  hCentSum3->Write();

  hZDCAsym1->Write();
  hZDCAsym2->Write();
  hZDCAsym3->Write();
  
}

void DrawTrkRatio(Char_t *hname){
  
  hTrkratio1->Fit("fpol1","0","",-0.8,0.8);
  hTrkratio2->Fit("fpol1","0","",-0.8,0.8);
  
  TF1* mypol1 = hTrkratio1->GetFunction("fpol1");
  TF1* mypol2 = hTrkratio2->GetFunction("fpol1");
  cout<<mypol1->GetChisquare()<<"/"<<mypol1->GetNDF()<<endl; 
  cout<<mypol2->GetChisquare()<<"/"<<mypol2->GetNDF()<<endl; 
  Char_t chisqpdof1[128],chisqpdof2[128];
  sprintf(chisqpdof1,"#chi^{2}/ndf=%4.3f/%d",mypol1->GetChisquare(),mypol1->GetNDF());
  sprintf(chisqpdof2,"#chi^{2}/ndf=%4.3f/%d",mypol2->GetChisquare(),mypol2->GetNDF());

  
  mypol1->SetRange(-0.8,0.8);
  mypol2->SetRange(-0.8,0.8);
  TCanvas * cTrk = new TCanvas("cTrk","cTrk",1000,600);
  cTrk->cd();
  cTrk->SetGrid();
  hTrkratio1->GetXaxis()->SetRangeUser(-0.8,0.8);
  hTrkratio1->SetMarkerStyle(8);
  hTrkratio1->SetMarkerColor(kBlue+1);
  hTrkratio2->SetMarkerStyle(8);
  hTrkratio2->SetMarkerColor(kRed+1);

  hTrkratio1->SetMaximum(1.004);
  hTrkratio1->SetMinimum(0.997);
  hTrkratio1->SetXTitle("#eta");
  hTrkratio1->SetYTitle("(dN/d#eta) Ratios");
  hTrkratio1->Draw();
  hTrkratio2->Draw("same");
  mypol1->SetLineColor(kBlue+1);
  mypol1->SetLineStyle(2);
  mypol1->Draw("same");

  mypol2->SetLineColor(kRed+1);
  mypol2->SetLineStyle(2);
  mypol2->Draw("same");

  TLatex * l = new TLatex(-0.7,1.003,"Pb-Pb 2.76 TeV LHC10h");
  l->SetTextSize(0.04);
  l->Draw();
   
  Char_t sCent[40];
  sprintf(sCent,"%d%% < #sigma < %d%%",centmin,centmax);
  l = new TLatex(-0.7,1.0024,sCent);
  l->SetTextSize(0.04);
  l->Draw();

  sprintf(sCent,"Using %s",hname);
  l = new TLatex(-0.7,1.0018,sCent);
  l->SetTextSize(0.04);
  l->Draw();
  
  TMarker * m = new TMarker(0.07,1.003,8);
  m->SetMarkerColor(kRed+1);
  m->Draw();
  TLatex * l = new TLatex(0.1,1.003,"(dN/d#eta)_{#alpha_{ZDC}>0.1}/(dN/d#eta)_{|#alpha_{ZDC}|<0.1}");
  l->SetTextSize(0.04);
  l->SetTextColor(kRed+1);
  l->Draw();
  l = new TLatex(.1,1.0024,"(dN/d#eta)_{#alpha_{ZDC}<-0.1}/(dN/d#eta)_{|#alpha_{ZDC}|<0.1}");
  l->SetTextSize(0.04);
  l->SetTextColor(kBlue+1);
  l->Draw();
  m = new TMarker(0.07,1.0024,8);
  m->SetMarkerColor(kBlue+1);
  m->Draw();
  
  Char_t sSlope[20];
  sprintf(sSlope,"%f+/-%f",mypol1->GetParameter(0),mypol1->GetParError(0));
  TLatex* l = new TLatex(-.20, 1.00085,sSlope);
  l->SetTextColor(kBlue+1);
  l->Draw();
  sprintf(sSlope,"%f+/-%f",mypol2->GetParameter(0),mypol2->GetParError(0));
  TLatex* l = new TLatex(-.20, 0.9985,sSlope);
  l->SetTextColor(kRed+1);
  l->Draw();
  l = new TLatex(-0.6,0.9975,chisqpdof1);
  l->SetTextColor(kBlue+1);
  l->Draw();
  l = new TLatex(0.4,0.9975,chisqpdof2);
  l->SetTextColor(kRed+1);
  l->Draw();
}

void DrawTrkV0Ratio(Char_t *hname1, Char_t *hname2, Char_t* FitOption, Char_t *ffname){

  cout<<" Fitfunction name is "<<ffname<<endl;
  
  hV0CReflect1 = new TH1F("hV0CReflect1","hV0CReflect1",220,-5.0,6.0);
  hV0CReflect2 = new TH1F("hV0CReflect2","hV0CReflect2",220,-5.0,6.0);
  
  hV0ACTrk1by3 = new TH1F("hV0ACTrk1by3","hV0ACTrk1by3",220,-5.0,6.0);
  hV0ACTrk2by3 = new TH1F("hV0ACTrk2by3","hV0ACTrk2by3",220,-5.0,6.0);
  
  hV0Ratio1 = new TH1F("hV0Ratio1","hV0Ratio1",220,-5.0,6.0);
  hV0Ratio2 = new TH1F("hV0Ratio2","hV0Ratio2",220,-5.0,6.0);

  
  for(Int_t iv = 0;iv<4;iv++){
    Int_t ietabin = (V0CEta[iv]+5)*20;
    hV0ACTrk1by3->SetBinContent(ietabin+1,V0Cratio1[iv]);
    hV0ACTrk1by3->SetBinError(ietabin+1,dV0Cratio1[iv]);
    hV0Ratio1->SetBinContent(ietabin+1,V0Cratio1[iv]);
    hV0Ratio1->SetBinError(ietabin+1,dV0Cratio1[iv]);
    
    ietabin = (V0AEta[iv]+5)*20;
    hV0ACTrk1by3->SetBinContent(ietabin+1,V0Aratio1[iv]);
    hV0ACTrk1by3->SetBinError(ietabin+1,dV0Aratio1[iv]);
    hV0Ratio1->SetBinContent(ietabin+1,V0Aratio1[iv]);
    hV0Ratio1->SetBinError(ietabin+1,dV0Aratio1[iv]);
    //  hV0AOnlyRatio->SetBinContent(ietabin+1,V0Aratio1[iv]);
    // hV0AOnlyRatio->SetBinError(ietabin+1,dV0Aratio1[iv]);

    ietabin = (-1.0*V0CEta[iv]+5)*20;
    hV0CReflect1->SetBinContent(ietabin+1,V0Cratio1[iv]);
    hV0CReflect1->SetBinError(ietabin+1,dV0Cratio1[iv]);

    hV0CReflect2->SetBinContent(ietabin+1,V0Cratio2[iv]);
    hV0CReflect2->SetBinError(ietabin+1,dV0Cratio2[iv]);

    //  hV0COnlyRatio->SetBinContent(ietabin+1,V0Cratio1[iv]);
    //hV0COnlyRatio->SetBinError(ietabin+1,dV0Cratio1[iv]);

    ietabin = (-1.0*V0AEta[iv]+5)*20;
    // hV0AOnlyRatio->SetBinContent(ietabin+1,V0Aratio2[iv]);
    //hV0AOnlyRatio->SetBinError(ietabin+1,dV0Aratio2[iv]);
    
  }
  
  Int_t netabins = hTrkratio1->GetXaxis()->GetNbins();
  for(Int_t ieta = 0;ieta<netabins-2;ieta++){
    Int_t ietabin= (TrkEta[ieta]+5)*20;
    hV0ACTrk1by3->SetBinContent(ietabin+1,Trkratio1[ieta]);
    hV0ACTrk1by3->SetBinError(ietabin+1,dTrkratio1[ieta]);
  }
  
  
  for(Int_t iv = 0;iv<4;iv++){
    Int_t ietabin = (V0CEta[iv]+5)*20;
    hV0ACTrk2by3->SetBinContent(ietabin+1,V0Cratio2[iv]);
    hV0ACTrk2by3->SetBinError(ietabin+1,dV0Cratio2[iv]);
    hV0Ratio2->SetBinContent(ietabin+1,V0Cratio2[iv]);
    hV0Ratio2->SetBinError(ietabin+1,dV0Cratio2[iv]);
    //  hV0COnlyRatio->SetBinContent(ietabin+1,V0Cratio2[iv]);
    //  hV0COnlyRatio->SetBinError(ietabin+1,dV0Cratio2[iv]);

    ietabin = (V0AEta[iv]+5)*20;
    hV0ACTrk2by3->SetBinContent(ietabin+1,V0Aratio2[iv]);
    hV0ACTrk2by3->SetBinError(ietabin+1,dV0Aratio2[iv]);
    hV0Ratio2->SetBinContent(ietabin+1,V0Aratio2[iv]);
    hV0Ratio2->SetBinError(ietabin+1,dV0Aratio2[iv]);
 
  }
  
  for(Int_t ieta = 0;ieta<netabins-2;ieta++){
    Int_t ietabin= (TrkEta[ieta]+5)*20;
    hV0ACTrk2by3->SetBinContent(ietabin+1,Trkratio2[ieta]);
    hV0ACTrk2by3->SetBinError(ietabin+1,dTrkratio2[ieta]);
  }
  /*
  for(Int_t ieta = 0;ieta<4;ieta++){
    Int_t ietabin= (TrkEta[ieta]+5)*20;
    hV0COnlyRatio->SetBinContent(ietabin+1,Trkratio2[ieta]);
    ietabin = (-1.0*TrkEta[ieta]+5)*20;
    hV0COnlyRatio->SetBinContent(ietabin+1,Trkratio1[ieta]);
  }
  for(Int_t ieta = 4;ieta<8;ieta++){
    Int_t ietabin= (TrkEta[ieta]+5)*20;
    hV0AOnlyRatio->SetBinContent(ietabin+1,Trkratio1[ieta]);
    ietabin = (-1.0*TrkEta[ieta]+5)*20;
    hV0AOnlyRatio->SetBinContent(ietabin+1,Trkratio2[ieta]);
  }
  */
  //DrawV0AC();

  
  TCanvas * c1 = new TCanvas("c1","c1");
  c1->cd();
  hV0ACTrk1by3->SetMinimum(0.97);
  hV0ACTrk1by3->SetMaximum(1.03);
  
  TF1 * mypol1, *mypol2, *mypol3, *mypol4;

  
  //For Fitting on Trk + V0A + V0C points.
  hV0ACTrk1by3->Fit("fpol1","0");
  fdata<<fpol1->GetParameter(0)<<" "<<fpol1->GetParError(0)<<endl;
  
  hV0ACTrk2by3->Fit("fpol1","0");
  fdata<<fpol1->GetParameter(0)<<" "<<fpol1->GetParError(0)<<endl;
  
  hV0ACTrk1by3->Fit("fpol3","0+");
  fdata<<fpol3->GetParameter(0)<<" "<<fpol3->GetParError(0)<<" "<<fpol3->GetParameter(1)<<" "<<fpol3->GetParError(1)<<endl;
  hV0ACTrk2by3->Fit("fpol3","0+");
  fdata<<fpol3->GetParameter(0)<<" "<<fpol3->GetParError(0)<<" "<<fpol3->GetParameter(1)<<" "<<fpol3->GetParError(1)<<endl;

  /*
  hV0ACTrk1by3->Fit("fpol5","0+");
  fdata<<fpol5->GetParameter(0)<<" "<<fpol5->GetParError(0)<<" "<<fpol5->GetParameter(1)<<" "<<fpol5->GetParError(1)<<" "<<fpol5->GetParameter(2)<<" "<<fpol5->GetParError(2)<<endl;
  hV0ACTrk2by3->Fit("fpol5","0+");
  fdata<<fpol5->GetParameter(0)<<" "<<fpol5->GetParError(0)<<" "<<fpol5->GetParameter(1)<<" "<<fpol5->GetParError(1)<<" "<<fpol5->GetParameter(2)<<" "<<fpol5->GetParError(2)<<endl;
  */
  if(FitOption=="Trk+V0A+V0C"){
    mypol1 = hV0ACTrk1by3->GetFunction(ffname);
    mypol2 = hV0ACTrk2by3->GetFunction(ffname);
    mypol3 = hV0ACTrk1by3->GetFunction("fpol3");
    mypol4 = hV0ACTrk2by3->GetFunction("fpol3");

    /*
    hV0Ratio1_hFit = (TH1F*)hV0Ratio1->Clone("hV0Ratio1_hFit");
    hV0Ratio1_hFit->Reset();
    for(Int_t ibin = 1;ibin<=hV0Ratio1_hFit->GetNbinsX();ibin++){
      if(hV0Ratio1->GetBinContent(ibin)==0) continue;
      Float_t databyfit = hV0Ratio1->GetBinContent(ibin)/mypol1(hV0Ratio1->GetBinCenter(ibin));
      hV0Ratio1_hFit->SetBinContent(ibin,databyfit);
    }
    */
  }
  
  // For Fitting on V0A + V0C points.
  hV0Ratio1->SetMarkerStyle(8);
  hV0Ratio2->SetMarkerStyle(8);
  hV0Ratio1->Fit("fpol1","0");
  fdata<<fpol1->GetParameter(0)<<" "<<fpol1->GetParError(0)<<endl;
  hV0Ratio2->Fit("fpol1","0");
  fdata<<fpol1->GetParameter(0)<<" "<<fpol1->GetParError(0)<<endl;
  hV0Ratio1->Fit("fpol3","0+");
  fdata<<fpol3->GetParameter(0)<<" "<<fpol3->GetParError(0)<<" "<<fpol3->GetParameter(1)<<" "<<fpol3->GetParError(1)<<endl;
  hV0Ratio2->Fit("fpol3","0+");
  fdata<<fpol3->GetParameter(0)<<" "<<fpol3->GetParError(0)<<" "<<fpol3->GetParameter(1)<<" "<<fpol3->GetParError(1)<<endl;
  if(FitOption=="V0A+V0C"){
    mypol1 = hV0Ratio1->GetFunction(ffname);
    mypol2 = hV0Ratio2->GetFunction(ffname);
    /*
    hV0Ratio1_hFit = (TH1F*)hV0Ratio1->Clone("hV0Ratio1_hFit");
    hV0Ratio1_hFit->Reset();
    for(Int_t ibin = 1;ibin<=hV0Ratio1_hFit->GetNbinsX();ibin++){
      if(hV0Ratio1->GetBinContent(ibin)==0) continue;
      Float_t databyfit = hV0Ratio1->GetBinContent(ibin)/mypol1(hV0Ratio1->GetBinCenter(ibin));
      hV0Ratio1_hFit->SetBinContent(ibin,databyfit);
    }
    */
  }

  
  // For fittiing to only Trk Points.
  hTrkratio1->Fit("fpol1","0","",-0.8,0.8);
  hTrkratio2->Fit("fpol1","0","",-0.8,0.8);
  fdata<<fpol1->GetParameter(0)<<" "<<fpol1->GetParError(0)<<endl;
  fdata<<fpol1->GetParameter(0)<<" "<<fpol1->GetParError(0)<<endl;
  hTrkratio1->Fit("fpol3","0+","",-0.8,0.8);
  hTrkratio2->Fit("fpol3","0+","",-0.8,0.8);
  fdata<<fpol3->GetParameter(0)<<" "<<fpol3->GetParError(0)<<" "<<fpol3->GetParameter(1)<<" "<<fpol3->GetParError(1)<<endl;
  fdata<<fpol3->GetParameter(0)<<" "<<fpol3->GetParError(0)<<" "<<fpol3->GetParameter(1)<<" "<<fpol3->GetParError(1)<<endl;
  if(FitOption=="Trk"){
    mypol1 = hTrkratio1->GetFunction(ffname);
    mypol2 = hTrkratio2->GetFunction(ffname);
    mypol3 = hTrkratio1->GetFunction("fpol3");
    mypol4 = hTrkratio2->GetFunction("fpol3");
  }

  mypol1->SetRange(-5.0,5.0);
  mypol2->SetRange(-5.0,5.0);
  mypol3->SetRange(-5.0,5.0);
  mypol4->SetRange(-5.0,5.0);
  Char_t chisqpdof1[128],chisqpdof2[128];
  sprintf(chisqpdof1,"#chi^{2}/ndf=%4.3f/%d",mypol1->GetChisquare(),mypol1->GetNDF());
  sprintf(chisqpdof2,"#chi^{2}/ndf=%4.3f/%d",mypol2->GetChisquare(),mypol2->GetNDF());
  
  
  TCanvas * c = new TCanvas("c","c",1000,600);
  c->cd();
  c->SetGrid();
  
  hV0ACTrk1by3->SetMarkerStyle(21);
  hV0ACTrk1by3->SetMarkerColor(kBlue+1);
  hV0ACTrk1by3->SetLineColor(kBlue+1);
  hV0ACTrk1by3->SetXTitle("#eta");
  hV0ACTrk1by3->SetYTitle("(dN/d#eta) Ratios");

  hV0ACTrk1by3->GetXaxis()->SetRangeUser(-5.0,5.0);
  hV0ACTrk1by3->Draw();
  mypol1->SetLineColor(kBlue+1);
  mypol1->SetLineStyle(2);
  mypol1->Draw("same");
  
  hV0ACTrk2by3->SetMarkerStyle(20);
  hV0ACTrk2by3->SetMarkerColor(kRed+1);
  hV0ACTrk2by3->SetLineColor(kRed+1);
  //  hV0ACTrk2by3->GetFunction("fpol1")->SetLineColor(kRed+1);
  //  hV0ACTrk2by3->GetFunction("fpol1")->SetLineStyle(2);
  //  hTrkratio2->GetFunction("pol1")->SetLineColor(kRed+1);
  //  hTrkratio2->GetFunction("pol1")->SetLineStyle(2);
  //  hTrkratio2->GetFunction("pol1")->Draw("same");
  
  hV0ACTrk2by3->Draw("same");
  mypol2->SetLineColor(kRed+1);
  mypol2->SetLineStyle(2);
  mypol2->Draw("same");
  
  hV0CReflect1->SetMarkerStyle(25);
  hV0CReflect1->SetMarkerColor(kBlue+1);
  hV0CReflect1->Draw("same");
  hV0CReflect2->SetMarkerStyle(24);
  hV0CReflect2->SetMarkerColor(kRed+1);
  hV0CReflect2->Draw("same");

  mypol3->SetLineColor(kBlue+1);
  mypol3->SetLineStyle(3);
  mypol3->Draw("same");
  mypol4->SetLineColor(kRed+1);
  mypol4->SetLineStyle(3);
  mypol4->Draw("same");
  
  l = new TLatex(-4.2,1.025,"Pb-Pb 2.76 TeV ");
  l->SetTextSize(0.04);
  l->Draw();
  
  Char_t sCent[40];
  sprintf(sCent,"%d%% < #sigma < %d%%",centmin,centmax);
  l = new TLatex(-4.2,1.021,sCent);
  l->SetTextSize(0.04);
  l->Draw();

  sprintf(sCent,"Track Data Using %s",hname1);
  l = new TLatex(-4.2,1.017,sCent);
  l->SetTextSize(0.04);
  // l->Draw();

  sprintf(sCent,"V0 Data Using %s",hname2);
  l = new TLatex(-4.2,1.013,sCent);
  l->SetTextSize(0.04);
  // l->Draw();

  l =  new TLatex(2.0,1.025,"ALICE PRELIMINARY");
  l->SetTextSize(0.04);
  l->SetTextColor(1);
  l->Draw();
  
  TMarker * m = new TMarker(-1.5,1.022,20);
  m->SetMarkerColor(kRed+1);
  m->Draw();
  //  TLatex * l = new TLatex(1,1.026,"(dN/d#eta)_{#alpha_{ZDC}>0.1}/(dN/d#eta)_{|#alpha_{ZDC}|<0.1}");
  TLatex * l = new TLatex(-1.2,1.021,"Region2/Region3");
  l->SetTextSize(0.04);
  l->SetTextColor(kRed+1);
  l->Draw();
  // l = new TLatex(1,1.021,"(dN/d#eta)_{#alpha_{ZDC}<-0.1}/(dN/d#eta)_{|#alpha_{ZDC}|<0.1}");
  l = new TLatex(1.2,1.016,chisqpdof2);
  l->SetTextColor(kRed+1);
  l->SetTextSize(0.035);
  l->Draw();
  Char_t sSlope[20];
  sprintf(sSlope,"c_{1}= %5.5f+/-%5.5f",mypol2->GetParameter(0),mypol2->GetParError(0));
  l = new TLatex(1.2, 1.021,sSlope);
  l->SetTextColor(kRed+1);
  l->SetTextSize(0.035);
  l->Draw();

  
  l = new TLatex(-1.2,0.98,"Region1/Region3");
  l->SetTextSize(0.04);
  l->SetTextColor(kBlue+1);
  l->Draw();
  m = new TMarker(-1.5,0.981,21);
  m->SetMarkerColor(kBlue+1);
  m->Draw();
  sprintf(sSlope,"c_{1}=%5.5f+/-%5.5f",mypol1->GetParameter(0),mypol1->GetParError(0));
  l = new TLatex(1.2, 0.98,sSlope);
  l->SetTextColor(kBlue+1);
  l->SetTextSize(0.035);
  l->Draw();
  l = new TLatex(1.2,0.975,chisqpdof1);
  l->SetTextSize(0.035);
  l->SetTextColor(kBlue+1);
  l->Draw();
}

void FillZDCHistoFromFile(TFile * fin){
  hZDCAsymAbs = new TH1F("hZDCAsymAbs","hZDCAsymAbs",40,0.0,1.0);
  hZDCAsymAbs->SetXTitle("|#alpha_{ZDC}|");
  hZDCAsymAbs->SetYTitle("Number of Events");
  
  hZDCAsym = new TH1F("hZDCAsym","hZDCAsym",80,-1.0,1.0);
  hZDCAsym->SetXTitle("#alpha_{ZDC}");
  hZDCAsym->SetYTitle("Number of Events");
  
  hZDCAsym1 = new TH1F("hZDCAsym1","hZDCAsym1",80,-1.0,1.0);
  hZDCAsym1->SetXTitle("#alpha_{ZDC} -1.0 to -0.1");
  hZDCAsym1->SetYTitle("Number of Events");
  
  hZDCAsym3 = new TH1F("hZDCAsym3","hZDCAsym3",80,-1.0,1.0);
  hZDCAsym3->SetXTitle("#alpha_{ZDC} -0.1 to 0.1");
  hZDCAsym3->SetYTitle("Number of Events");
  
  hZDCAsym2 = new TH1F("hZDCAsym2","hZDCAsym2",80,-1.0,1.0);
  hZDCAsym2->SetXTitle("#alpha_{ZDC} 0.1 to 1.0");
  hZDCAsym2->SetYTitle("Number of Events");

  fin->cd();
  tree = (TTree*)fin->Get("asymTree");
  tree->SetBranchAddress("zdcasym",&zdcasym);
  Int_t nentries =tree->GetEntries();
  for(Int_t iev = 0;iev<nentries;iev++){
    tree->GetEntry(iev);
    hZDCAsym->Fill(zdcasym);
    hZDCAsymAbs->Fill(TMath::Abs(zdcasym));
    if(zdcasym < -1.0*Reg12Limit && zdcasym > -1.0) hZDCAsym1->Fill(zdcasym);
    if(zdcasym > -1.0*Reg3Limit && zdcasym < Reg3Limit) hZDCAsym3->Fill(zdcasym);
    if(zdcasym >  1.0*Reg12Limit && zdcasym < 1.0) hZDCAsym2->Fill(zdcasym);
  }
  
  if(OptDrawVzCentZDC==1) DrawZDCAsym();
}


TH1F * FillHistoSum4ZDCAsym(Int_t minbin,Int_t maxbin, Char_t * hname, TFile * f4sum){
  //  cout<<" minbin = "<<minbin<<" maxbin = "<<maxbin<<endl;
  TH1F * h=NULL;
  Int_t nevent = 0;
  Char_t hnamebin[40], hnamesum[40];
  //  cout<<"Inputfile  "<<f4sum->GetName()<<endl;
  for(Int_t ibin = minbin;ibin<=maxbin;ibin++){
    sprintf(hnamebin,"%s%i",hname,ibin);
    //    cout<<" Reading "<<hnamebin<<" from file "<<f4sum<<endl;
    h = (TH1F *)f4sum->Get(hnamebin);
    if(ibin==minbin){
      sprintf(hnamebin,"%s%ito%i",hname,minbin,maxbin);
      //      cout<<"hnamebin = "<<hnamebin<<endl;
      hSum = (TH1F*)h->Clone(hnamebin);
    }else{
      hSum->Add((TH1F*)f4sum->Get(hnamebin));
    }
    nevent+=h->GetEntries();
  }
  //  cout<<"Entries in "<<hname<< " is "<<hSum->GetEntries()<<endl;
  //  cout<<" nevent = "<<nevent<<endl;

  return hSum;
}

TProfile * FillEta4ZDCAsym(Int_t minbin, Int_t maxbin, Char_t* hname, TFile * f4sum ){
  Int_t SumWeight = 0;
  f4sum->cd();
  TProfile * hp=NULL;
  Char_t hnamebin[48];
  //  cout<<"Inputfile  "<<f4sum->GetName()<<endl;
  for(Int_t ibin = minbin;ibin<=maxbin;ibin++){
    sprintf(hnamebin,"%s%i",hname,ibin);
    //    cout<<" Reading histogram "<<hnamebin<<" From file "<<f4sum->GetName()<<endl;
    hp= (TProfile*)f4sum->Get(hnamebin);
    SumWeight+=hp->GetSumOfWeights();
    if(ibin==minbin) {
      sprintf(hnamebin,"%s%ito%i",hname,minbin,maxbin);
      hpSum = (TProfile*)hp->Clone(hnamebin);
    }else{
      hpSum->Add(hp);
    }
  }
  //  cout<<"Entries in "<<hpSum->GetName()<<" are "<<hpSum->GetEntries()<<endl;
  //  cout<<" SumWeights = "<<SumWeight<<" from hsum "<<hpSum->GetSumOfWeights()<<endl;
  //  hpTrkReg->Scale(1.0/SumWeight);
  return hpSum;
}


void DrawReg1and2Plots(Char_t *hname1, Char_t *hname2 ){
  TH1F * h1and2 = (TH1F*)hV0ACTrk1by3->Clone();
  h1and2->Reset();
  Int_t nbin = hV0ACTrk1by3->GetNbinsX();
  for(Int_t i = 1;i<=nbin;i++){
    if(hV0ACTrk1by3->GetBinContent(i)==0) continue;
    Float_t y1 = hV0ACTrk1by3->GetBinContent(i);
    Float_t dy1 = hV0ACTrk1by3->GetBinError(i);
    Float_t y2 = hV0ACTrk2by3->GetBinContent(i);
    y2 = 2.0-y2;
    Float_t dy2 =hV0ACTrk2by3->GetBinError(i);
    Float_t mean = 0.5*(y1+y2);
    Float_t Error = 0.5*TMath::Sqrt(dy1*dy1+dy2*dy2);
    cout<<y1<<" "<<dy1<<" "<<y2<<" "<<dy2<<"   "<<mean<<" "<<Error<<endl;
    
    h1and2->SetBinContent(i,mean);
    h1and2->SetBinError(i,Error);
  }
  //  cout<<"Entries in "<<hpSum->GetName()<<" are "<<hpSum->GetEntries()<<endl;
  //  cout<<" SumWeights = "<<SumWeight<<" from hsum "<<hpSum->GetSumOfWeights()<<endl;
  //  hpTrkReg->Scale(1.0/SumWeight);
  h1and2->Fit("fpol1","0");
  h1and2->Fit("fpol3","0+");

  TCanvas * c12 = new TCanvas("c12","c12");
  c12->cd();
  c12->SetGrid();
  h1and2->SetMinimum(0.97);
  h1and2->SetMaximum(1.03);
  h1and2->SetMarkerColor(1);
  h1and2->Draw();
  
  TF1 * mypol1, *mypol2, *mypol3;
  mypol1 = h1and2->GetFunction("fpol1");
  mypol2 = h1and2->GetFunction("fpol3");
  Char_t chisqpdof1[128],chisqpdof2[128];
  sprintf(chisqpdof1,"#chi^{2}/ndf=%4.3f/%d",mypol1->GetChisquare(),mypol1->GetNDF());
  sprintf(chisqpdof2,"#chi^{2}/ndf=%4.3f/%d",mypol2->GetChisquare(),mypol2->GetNDF());
  mypol2->SetLineColor(1);
  mypol2->SetLineStyle(2);
  mypol2->Draw("same");
  
  l = new TLatex(-4.2,1.025,"Pb-Pb 2.76 TeV LHC10h");
  l->SetTextSize(0.04);
  l->Draw();
  
  Char_t sCent[40];
  sprintf(sCent,"%d%% < #sigma < %d%%",centmin,centmax);
  l = new TLatex(-4.2,1.021,sCent);
  l->SetTextSize(0.04);
  l->Draw();

  sprintf(sCent,"Track Data Using %s",hname1);
  l = new TLatex(-4.2,1.017,sCent);
  l->SetTextSize(0.04);
  l->Draw();

  sprintf(sCent,"V0 Data Using %s",hname2);
  l = new TLatex(-4.2,1.013,sCent);
  l->SetTextSize(0.04);
  l->Draw();

  TMarker * m = new TMarker(0.7,1.025,8);
  m->SetMarkerColor(kRed+1);
  m->Draw();
  TLatex * l = new TLatex(1,1.025,"(dN/d#eta)_{#alpha_{ZDC}>0.1}/(dN/d#eta)_{|#alpha_{ZDC}|<0.1}");
  l->SetTextSize(0.04);
  l->SetTextColor(kRed+1);
  l->Draw();
  l = new TLatex(1,1.021,"(dN/d#eta)_{#alpha_{ZDC}<-0.1}/(dN/d#eta)_{|#alpha_{ZDC}|<0.1}");
  l->SetTextSize(0.04);
  l->SetTextColor(kBlue+1);
  l->Draw();
  m = new TMarker(0.7,1.021,8);
  m->SetMarkerColor(kBlue+1);
  m->Draw();

  Char_t sSlope[20];
  sprintf(sSlope,"%f+/-%f",mypol2->GetParameter(0),mypol2->GetParError(0));
  TLatex* l = new TLatex(-2.0, 1.0085,sSlope);
  l->SetTextColor(1);
  l->Draw();
  l = new TLatex(-4.0,0.975,chisqpdof2);
  l->SetTextColor(kBlue+1);
  l->Draw();

}


void DrawZDCAsym(){
  TCanvas * cZDCAsym = new TCanvas("cZDCAsym","cZDCAsym",600,600);
  cZDCAsym->cd();
  hZDCAsym->Draw();
  hZDCAsym3->SetFillStyle(2);
  hZDCAsym3->SetLineColor(1);
  hZDCAsym3->SetFillColor(1);
  hZDCAsym3->Draw("same");
  hZDCAsym3->SetMaximum(60000);
  hZDCAsym1->SetFillStyle(2);
  hZDCAsym1->SetLineColor(kBlue+1);
  hZDCAsym1->SetFillColor(kBlue+1);
  hZDCAsym1->Draw("same");
  hZDCAsym2->SetFillStyle(2);
  hZDCAsym2->SetLineColor(kRed+1);
  hZDCAsym2->SetFillColor(kRed+1);
  hZDCAsym2->Draw("same");
  hZDCAsym3->Draw("same");
  TLatex * l = new TLatex(-0.5,10000," Region1");
  l->SetTextColor(kBlue+1);
  l->Draw();
  l = new TLatex(0.15,10000," Region2");
  l->SetTextColor(kRed+1);
  l->Draw();
  l = new TLatex(-0.15,5000,"Region3");
  l->SetTextColor(1);
  l->Draw();
  
  l = new TLatex(-0.9,48000,"LHC10h Pb-Pb 2.76 TeV");
  l->SetTextSize(0.04);
  l->Draw();
  Char_t sCent[40];
  sprintf(sCent,"%d%% < #sigma < %d%%",centmin,centmax);
  l = new TLatex(-0.9,44000,sCent);
  l->SetTextSize(0.04);
  l->Draw();

  TCanvas * cZDCAsymAbs = new TCanvas("cZDCAsymAbs","cZDCAsymAbs",600,600);
  cZDCAsymAbs->cd();
  hZDCAsymAbs->Fit("gaus");

  return;
}

void FillVzFromFile(TFile * fin){
  // Getting stat errors on Vz plot

  Char_t hname[]="hVertexZ_ZDCAsymBin";
  //  cout<<" Reg1MinBin ="<<Reg1MinBin<<" Reg1MaxBin = "<<Reg1MaxBin<<endl;
  hVzSum1 = FillHistoSum4ZDCAsym(Reg1MinBin,Reg1MaxBin,hname,fin);
  //  cout<<" Number of events in Region1  are "<<hVzSum1->GetEntries()<<endl;
  hVzSum3 = FillHistoSum4ZDCAsym(Reg3MinBin,Reg3MaxBin,hname,fin);
  //  cout<<" Number of events in REgion3 are "<<hVzSum3->GetEntries()<<endl;
  hVzSum2 = FillHistoSum4ZDCAsym(Reg2MinBin,Reg2MaxBin,hname,fin);
  //  cout<<" Number of events in Region2 are "<<hVzSum2->GetEntries()<<endl;
  if(OptDrawVzCentZDC==1) {
    TCanvas * cvz = new TCanvas("cvz","cvz",600,600);
    hVzSum1->SetLineColor(kBlue+1);
    hVzSum1->SetMarkerColor(kBlue+1);
    hVzSum1->Draw();
    hVzSum2->SetLineColor(kRed+1);
    hVzSum2->SetMarkerColor(kRed+1);
    hVzSum2->Draw("same");
    hVzSum3->Draw("same");
  }
  return;
}

void FillCentFromFile(TFile * fin){
  // Getting stat errors on Cent plot
  Char_t hname[]="hCent_ZDCAsymBin";
  hCentSum1 = FillHistoSum4ZDCAsym(Reg1MinBin,Reg1MaxBin,hname,fin);
  //  cout<<" Number of events in Region1  are "<<hCentSum1->GetEntries()<<endl;
  hCentSum3 = FillHistoSum4ZDCAsym(Reg3MinBin,Reg3MaxBin,hname,fin);
  //  cout<<" Number of events in REgion3 are "<<hCentSum3->GetEntries()<<endl;
  hCentSum2 = FillHistoSum4ZDCAsym(Reg2MinBin,Reg2MaxBin,hname,fin);
  //  cout<<" Number of events in Region2 are "<<hCentSum2->GetEntries()<<endl;

  if(OptDrawVzCentZDC==1) {
    TCanvas * cCent = new TCanvas("cCent","cCent",600,600);
    hCentSum1->SetLineColor(kBlue+1);
    hCentSum1->SetMarkerColor(kBlue+1);
    hCentSum1->Draw();
    hCentSum2->SetLineColor(kRed+1);
    hCentSum2->SetMarkerColor(kRed+1);
    hCentSum2->Draw("same");
    hCentSum3->Draw("same");
  }
  return;
}

void FillTrkFromFile(TFile * fin){

  //  cout<<" Reading Track data from file "<<fin->GetName()<<endl;
  Char_t hname[]="hpTrackEta_ZDCAsymBin";
  hpTrkReg1 = FillEta4ZDCAsym(Reg1MinBin,Reg1MaxBin,hname,fin);
  hpTrkReg3 = FillEta4ZDCAsym(Reg3MinBin,Reg3MaxBin,hname,fin);
  hpTrkReg2 = FillEta4ZDCAsym(Reg2MinBin,Reg2MaxBin,hname,fin);
  hpTrkReg123 = FillEta4ZDCAsym(Reg1MinBin,Reg2MaxBin,hname,fin);
  //  hpTrkReg1and2 = FillReg1and2Eta4ZDCAsym(hpTrkReg1,hpTrkReg2,hname,fin);
  
  hpTrkReg1->Rebin(2);
  hpTrkReg2->Rebin(2);
  hpTrkReg3->Rebin(2);

  return;
}

void FillV0AFromFile(TFile * fin){

  //  cout<<" Reading V0A data from file "<<fin->GetName()<<endl;
  Char_t hname[]="hpV0AEta_ZDCAsymBin";
  hpV0AReg1 = FillEta4ZDCAsym(Reg1MinBin,Reg1MaxBin,hname,fin);
  hpV0AReg3 = FillEta4ZDCAsym(Reg3MinBin,Reg3MaxBin,hname,fin);
  hpV0AReg2 = FillEta4ZDCAsym(Reg2MinBin,Reg2MaxBin,hname,fin);
  hpV0AReg123 = FillEta4ZDCAsym(Reg1MinBin,Reg2MaxBin,hname,fin);
  //  hpV0AReg1and2 = FillReg1and2Eta4ZDCAsym(hpV0AReg1,hpV0AReg2,hname,fin);
  
  return;
}

void FillV0CFromFile(TFile * fin){

  //  cout<<" Reading V0C data from file "<<fin->GetName()<<endl;
  Char_t hname[]="hpV0CEta_ZDCAsymBin";
  hpV0CReg1 = FillEta4ZDCAsym(Reg1MinBin,Reg1MaxBin,hname,fin);
  hpV0CReg3 = FillEta4ZDCAsym(Reg3MinBin,Reg3MaxBin,hname,fin);
  hpV0CReg2 = FillEta4ZDCAsym(Reg2MinBin,Reg2MaxBin,hname,fin);
  hpV0CReg123 = FillEta4ZDCAsym(Reg1MinBin,Reg2MaxBin,hname,fin);
  //  hpV0CReg1and2 = FillReg1and2Eta4ZDCAsym(hpV0CReg1,hpV0CReg2,hname,fin);
    
  return;
}

void GetGlobalCorrection2V0(){

  V0GlobalCorrReg1 = (hpV0CReg3->GetSumOfWeights()+hpV0AReg3->GetSumOfWeights())/(hpV0CReg1->GetSumOfWeights()+hpV0AReg1->GetSumOfWeights());
  V0GlobalCorrReg2 = (hpV0CReg3->GetSumOfWeights()+hpV0AReg3->GetSumOfWeights())/(hpV0CReg2->GetSumOfWeights()+hpV0AReg2->GetSumOfWeights());

  TrkGlobalCorrReg1 = (hpTrkReg3->GetSumOfWeights())/(hpTrkReg1->GetSumOfWeights());
  TrkGlobalCorrReg2 = (hpTrkReg3->GetSumOfWeights())/(hpTrkReg2->GetSumOfWeights());
  
  return;
}

void DrawV0AC(){
  hV0AOnlyRatio = new TH1F("hV0AOnlyRatio","hV0AOnlyRatio",220,-5.0,6.0);
  hV0COnlyRatio = new TH1F("hV0COnlyRatio","hV0COnlyRatio",220,-5.0,6.0);
  
  for(Int_t iv = 0;iv<4;iv++){
    Int_t ietabin = (V0AEta[iv]+5)*20;
    hV0AOnlyRatio->SetBinContent(ietabin+1,V0Aratio1[iv]);
    hV0AOnlyRatio->SetBinError(ietabin+1,dV0Aratio1[iv]);
    ietabin = (-1.0*V0CEta[iv]+5)*20;
    hV0COnlyRatio->SetBinContent(ietabin+1,V0Cratio1[iv]);
    hV0COnlyRatio->SetBinError(ietabin+1,dV0Cratio1[iv]);
    ietabin = (-1.0*V0AEta[iv]+5)*20;
    hV0AOnlyRatio->SetBinContent(ietabin+1,V0Aratio2[iv]);
    hV0AOnlyRatio->SetBinError(ietabin+1,dV0Aratio2[iv]);
  }
  for(Int_t iv = 0;iv<4;iv++){
    Int_t ietabin = (V0CEta[iv]+5)*20;
    hV0COnlyRatio->SetBinContent(ietabin+1,V0Cratio2[iv]);
    hV0COnlyRatio->SetBinError(ietabin+1,dV0Cratio2[iv]);
  }
  for(Int_t ieta = 0;ieta<4;ieta++){
    Int_t ietabin= (TrkEta[ieta]+5)*20;
    hV0COnlyRatio->SetBinContent(ietabin+1,Trkratio2[ieta]);
    hV0COnlyRatio->SetBinError(ietabin+1,dTrkratio2[ieta]);
    ietabin = (-1.0*TrkEta[ieta]+5)*20;
    hV0COnlyRatio->SetBinContent(ietabin+1,Trkratio1[ieta]);
    hV0COnlyRatio->SetBinError(ietabin+1,dTrkratio1[ieta]);
  }
  for(Int_t ieta = 4;ieta<8;ieta++){
    Int_t ietabin= (TrkEta[ieta]+5)*20;
    hV0AOnlyRatio->SetBinContent(ietabin+1,Trkratio1[ieta]);
    hV0AOnlyRatio->SetBinError(ietabin+1,dTrkratio1[ieta]);
    ietabin = (-1.0*TrkEta[ieta]+5)*20;
    hV0AOnlyRatio->SetBinContent(ietabin+1,Trkratio2[ieta]);
    hV0AOnlyRatio->SetBinError(ietabin+1,dTrkratio2[ieta]);
  }
  
  TCanvas * cV0AC = new TCanvas("cV0AC","cV0AC",1000,600);
  cV0AC->SetGrid();
  hV0AOnlyRatio->SetMarkerStyle(21);
  hV0AOnlyRatio->SetMarkerColor(kBlue);
  hV0AOnlyRatio->SetMaximum(1.03);
  hV0AOnlyRatio->SetMinimum(0.97);
  hV0COnlyRatio->SetMarkerStyle(20);
  hV0COnlyRatio->SetMarkerColor(kRed);
  hV0AOnlyRatio->Fit("fpol3");
  hV0COnlyRatio->Fit("fpol3");
  hV0AOnlyRatio->Draw();
  hV0COnlyRatio->Draw("same");
  TF1* mypol7 = hV0AOnlyRatio->GetFunction("fpol3");
  TF1* mypol8 = hV0COnlyRatio->GetFunction("fpol3");
  cout<<mypol7->GetChisquare()<<"/"<<mypol7->GetNDF()<<endl; 
  cout<<mypol8->GetChisquare()<<"/"<<mypol8->GetNDF()<<endl;
  
  Char_t chisqpdof1[128],chisqpdof2[128];
  sprintf(chisqpdof1,"#chi^{2}/ndf=%4.3f/%d",mypol7->GetChisquare(),mypol7->GetNDF());
  sprintf(chisqpdof2,"#chi^{2}/ndf=%4.3f/%d",mypol8->GetChisquare(),mypol8->GetNDF());
  Char_t sSlope[20];
  sprintf(sSlope,"c_{1}= %5.5f+/-%5.5f",mypol7->GetParameter(0),mypol7->GetParError(0));
  l = new TLatex(0, 1.024,sSlope);
  l->SetTextColor(kBlue);
  l->SetTextSize(0.035);
  l->Draw();
  sprintf(sSlope,"c_{1}= %5.5f+/-%5.5f",mypol8->GetParameter(0),mypol8->GetParError(0));
  l = new TLatex(0, 0.984,sSlope);
  l->SetTextColor(kRed);
  l->SetTextSize(0.035);
  l->Draw();
  
  TLatex * l  = new TLatex(0.0,1.02,chisqpdof1);
  l->SetTextColor(kBlue);
  l->Draw();
  l  = new TLatex(0.0,0.98,chisqpdof2);
  l->SetTextColor(kRed);
  l->Draw();
  l = new TLatex(-2.0,1.02,"V0A");
  l->SetTextColor(kBlue);
  l->Draw();
  l = new TLatex(-2.0,0.98,"V0C");
  l->SetTextColor(kRed);
  l->Draw();
  
}
