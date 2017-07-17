TString strPtAss[3]={"0.3to1.0","1.0to99.0","0.3to99.0"};
TString strPtAssCanvas[3]={"0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}","#it{p}_{T}^{assoc} > 1 GeV/#it{c}","#it{p}_{T}^{assoc} > 0.3 GeV/#it{c}"};
TString strSystem[2]={"pp","pPb"};
Color_t colSystem[2]={kBlack,kRed};
Int_t markerStyle[2]={20,21};
Bool_t useLegendForData=kTRUE;
Bool_t plotv2unc=kTRUE;
TString strFitResultPP=""; //  "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";
TString strFitResultPPb=""; //  "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";
Double_t canvasheight=801;
Double_t resizeTextFactor=1.;//1.*canvasheight/1200.; // size was tuned for canvasheight =1800. 
Double_t referencePadHeight=0.48; // Do not touch unless the canvas division is changed from 2x3, change canvasheight and resizeTextFactor instead
Double_t innerPadHeight;// not touch, set internally
Double_t innerPadWidth;// not touch, set internally
Int_t style=1;
Double_t scaleHeightPads=1;// do not touch this, it is regulated automatically in the macro
Double_t scaleWidthPads=1;// do not touch this, it is regulated automatically in the macro
TString strFitResultMC[2]={"",""};
//TString strFitResultMC[2]={"/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/MCTemplates/Templates_pp_12May15/FitResults/",
//			   "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May19UseScriptPWGHF/MCTemplates/Templates_pPb_12May15/FitResults/"};

TString strquantityFile[5]={"NSYield","NSSigma","Pedestal","ASYield","ASSigma"};
Double_t maxRangePP[3][5]={{2.5,0.64,4,4.3,1.48},
			   {2.5,0.64,3,4.3,1.48},
			   {3.3,0.64,4.4,4.3,1.48}};

Double_t maxRangePPb[3][5]={{2.5,0.64,13,4.3,1.48},
			   {2.5,0.64,5,4.3,1.48},
			   {3.3,0.64,13,4.3,1.48}};
Double_t minptMC=0.,maxptMC=99.;
Double_t minptData=0.,maxptData=99.;
Int_t ncolumns=3;
Int_t nrows=2;
Double_t ytitleoffset=2.45,xtitleoffset=1.95;
Bool_t skip3to5=kTRUE;
Double_t markersize=1.5;
Double_t markersizeMC=1.2;
Bool_t drawSystMC=kTRUE;
Bool_t runonpPb2016=kFALSE;

void SetMinPtDisplayData(Double_t ptmin){minptMC=ptmin;}
void SetMaxPtDisplayData(Double_t ptmax){maxptMC=ptmax;}
void SetMinPtDisplayMC(Double_t ptmin){minptData=ptmin;}
void SetMaxPtDisplayMC(Double_t ptmax){maxptData=ptmax;}
void SetRunOn2016(Bool_t run){runonpPb2016=run;}

void AdaptRangeHist(TH1D *h,Double_t minpt,Double_t maxpt){
  for(Int_t i=1;i<=h->GetNbinsX();i++){
    if(h->GetXaxis()->GetBinUpEdge(i)<minpt*1.001 ||  h->GetXaxis()->GetBinLowEdge(i)>maxpt*0.999 ){
      h->SetBinContent(i,0);
      h->SetBinError(i,0);
    }
  }
}

void AdaptRangeTGraph(TGraph *gr,Double_t minpt,Double_t maxpt){
  Double_t *x=gr->GetX();
  Double_t *exl=gr->GetEXlow();
  Double_t *exh=gr->GetEXhigh();
  Int_t np=gr->GetN();
  for(Int_t i=0;i<np;i++){
    if(x[i]+exh[i]<minpt*1.0001){
      gr->RemovePoint(i);
      AdaptRangeTGraph(gr,minpt,maxpt);
      return;
    }
    if(x[i]-exl[i]>maxpt*0.999){
      gr->RemovePoint(i);
      AdaptRangeTGraph(gr,minpt,maxpt);
      return;
    }
  }
  return;
}

void SetDrawSystMC(Bool_t drawsystmc){
  drawSystMC=drawsystmc;
}
void Init3x3Settings(){
  referencePadHeight=0.35;
  resizeTextFactor=0.85;
  ytitleoffset=4.8;
  xtitleoffset=4;
  ncolumns=3;
  nrows=3;
  markersize=1.;
  markersizeMC=1.;
}

TString yaxisTitle[5]={"Associated yield","#sigma_{fit,NS} (rad)","Baseline (rad^{-1})","Associated yield","#sigma_{fit,AS} (rad)"};
Double_t leftMarginCanvas=0.17;
Double_t rightMarginCanvas=0.055;
Double_t bottomMarginCanvas=0.13;
Double_t topMarginCanvas=0.07;
const Int_t nmodels=8;
Bool_t includemodel[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE,kTRUE};
TString strModelDir[nmodels]={"Perugia0","Perugia2010","Perugia2011","PYTHIA8","HERWIG","POWHEG","POWHEG","EPOS3"};
TString strModelDirLeg[nmodels]={"PYTHIA6, Perugia 0","PYTHIA6, Perugia 2010","PYTHIA6, Perugia 2011","PYTHIA8, Tune 4C","HERWIG","POWHEG+PYTHIA6","POWHEG+PYTHIA6 EPS09","EPOS 3.117"};
Color_t modelColors[nmodels]={kRed+2,kCyan,kGreen+2,kMagenta+1,kViolet,kBlue,kBlue,kOrange+1};
Bool_t includeinlegend[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE,kTRUE,kFALSE};// this is also used to split the legend in 2!!
Int_t modelMarkerStyle[nmodels]={kOpenSquare,kOpenCircle,kOpenDiamond,28,26,3,3,33};
TH1D **hMC;
TGraphAsymmErrors **grMC;

void SetSkip3to5(Bool_t skip){
  skip3to5=skip;
}
void IncludePerugia0(Bool_t incl){
  includemodel[0]=incl;
  includeinlegend[0]=incl;
}
void IncludePerugia2010(Bool_t incl){
  includemodel[1]=incl;
  includeinlegend[1]=incl;
}
void IncludePerugia2011(Bool_t incl){
  includemodel[2]=incl;
  includeinlegend[2]=incl;
}
void IncludePythia8(Bool_t incl){
  includemodel[4]=incl;
  includeinlegend[4]=incl;
}
void IncludePowheg(Bool_t incl){
  includemodel[3]=incl;
  includeinlegend[3]=incl;
}
void IncludePowhegEPS09(Bool_t incl){
  includemodel[6]=incl;
  includeinlegend[6]=incl;
}
void IncludeHerwig(Bool_t incl){
  includemodel[5]=incl;
  includeinlegend[5]=incl;
}
void IncludeEPOS(Bool_t incl){
  includemodel[7]=incl;
  includeinlegend[7]=incl;
}
void IncludeModel(Int_t imod,Bool_t incl){
  includemodel[imod]=incl;
  includeinlegend[imod]=incl;
}
Int_t GetNPythia6Tunes(){
  Int_t ncount=0;
  for(Int_t j=0;j<3;j++){
    if(includemodel[j])ncount++;
  }
  return ncount;
}

Int_t CountNmodels(){
  Int_t ncount=0;
  for(Int_t j=0;j<nmodels;j++){
    if(includemodel[j])ncount++;
  }
  return ncount;
}

Int_t CountNmodelsInLegend(){
  Int_t ncount=0;
  for(Int_t j=0;j<nmodels;j++){
    if(includeinlegend[j])ncount++;
    if(j==6&&includeinlegend[j])ncount++;
  }
  return ncount;
}

void SetDirectoryFitResultPP(TString str){
  strFitResultPP=str;
}

void SetDirectoryFitResultPPb(TString str){
  strFitResultPPb=str;
}

void SetDirectoryFitResultsMCpp(TString str){
  strFitResultMC[0]=str;
}

void SetDirectoryFitResultsMCpPb(TString str){
  strFitResultMC[1]=str;
}

TCanvas *CreateCanvasWithDefaultStyle(TString name){
  gStyle->SetOptStat(0000);
  TString nameNoPoints=name;
  nameNoPoints.ReplaceAll(".","");
  TCanvas *cout=new TCanvas(nameNoPoints.Data(),name.Data(),800,800);
  cout->SetTicky();
  cout->SetTickx();
  cout->SetFrameBorderMode(0);
  cout->SetLeftMargin(leftMarginCanvas);
  cout->SetRightMargin(rightMarginCanvas);
  cout->SetBottomMargin(bottomMarginCanvas);
  cout->SetTopMargin(topMarginCanvas);

  return cout;
}

void ConvertTH1ToTGraphAsymmError(TH1D* h,TGraphAsymmErrors *&gr, Double_t shift) {
  const Int_t nbinsxx=2;
  Double_t x[nbinsxx], y[nbinsxx], ex1[nbinsxx], ex2[nbinsxx], ey1[nbinsxx], ey2[nbinsxx];

  for(int i=0; i<nbinsxx; i++) {
    x[i] = h->GetBinCenter(i+2)+shift;
    y[i] = h->GetBinContent(i+2);
    ex1[i] = h->GetBinCenter(i+2)-h->GetBinLowEdge(i+2)+shift;
    ex2[i] = ex1[i]-2*shift;
    ey1[i] = h->GetBinError(i+2);
    ey2[i] = h->GetBinError(i+2);
  }

  gr = new TGraphAsymmErrors(nbinsxx,x,y,ex1,ex2,ey1,ey2);
  return;

}

void ConvertTH1ToTGraphAsymmError2016(TH1D* h,TGraphAsymmErrors *&gr, Double_t shift) {
  const Int_t nbinsxx=3;
  Double_t x[nbinsxx], y[nbinsxx], ex1[nbinsxx], ex2[nbinsxx], ey1[nbinsxx], ey2[nbinsxx];

  for(int i=0; i<nbinsxx; i++) {
    x[i] = h->GetBinCenter(i+2)+shift;
    y[i] = h->GetBinContent(i+2);
    ex1[i] = h->GetBinCenter(i+2)-h->GetBinLowEdge(i+2)+shift;
    ex2[i] = ex1[i]-2*shift;
    ey1[i] = h->GetBinError(i+2);
    ey2[i] = h->GetBinError(i+2);
  }

  gr = new TGraphAsymmErrors(nbinsxx,x,y,ex1,ex2,ey1,ey2);
  return;

}

void SetPadStyle(TPad *p){
  p->SetTicky();
  p->SetTickx();
  p->SetFrameBorderMode(0);
  p->SetFillStyle(0);
  p->SetFrameFillStyle(4000);
  //  p->SetBottomMargin(0);
  //  p->SetTopMargin(0);
  //  p->Range(-2.6,-1.5,5.4,4.4);
  //  p->SetLeftMargin(leftMarginCanvas);
  //  p->SetRightMargin(rightMarginCanvas);
  //  p->SetBottomMargin(bottomMarginCanvas);
  //  p->SetTopMargin(topMarginCanvas);
}

void Set2x3PadPositions(TCanvas* c){
    
  TPad * pd1 = (TPad*)c->GetPad(1);
  TPad * pd2 =(TPad*) c->GetPad(2);
  TPad * pd3 =(TPad*) c->GetPad(3);
  TPad * pd4 =(TPad*) c->GetPad(4);
  TPad * pd5 =(TPad*) c->GetPad(5);
  TPad * pd6 =(TPad*) c->GetPad(6);
  
  SetPadStyle(pd1);
  SetPadStyle(pd2);
  SetPadStyle(pd3);
  SetPadStyle(pd4);
  SetPadStyle(pd5);
  SetPadStyle(pd6);
  

  Double_t xl,xu,yl,yu;
  Double_t marginLeft=0.08;
  Double_t marginRight=0.04;
  Double_t marginTop=0.04;
  Double_t marginBottom=0.09;
  Double_t marginLeftForXAxis=0.02;
  Double_t marginBottomForYAxis=0.02;
  innerPadWidth=(1-marginLeft-marginRight)/3.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/2.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!


    // Bottom row
    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd4->SetRightMargin(0.);
    pd4->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd4->SetTopMargin(0.);

    pd4->Modified();
    pd4->Update();

    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd5->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2.*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd5->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd5->SetRightMargin(0.);
    pd5->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd5->SetTopMargin(0.);

    pd5->Modified();
    pd5->Update();
    
    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd6->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd6->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd6->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd6->SetTopMargin(0.);

    pd6->Modified();
    pd6->Update();

    // Top Row
    pd1->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd1->SetPad(0,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,1.);
    pd1->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd1->SetRightMargin(0.);
    pd1->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd1->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd1->Modified();
    pd1->Update();

    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,2.*innerPadWidth+marginLeft,1);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    //    pd2->SetLeftMargin(0.);
    pd2->SetRightMargin(0);
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd2->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(2.*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    //    pd3->SetLeftMargin(0.);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd3->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd3->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();
    

}





void Set3x3PadPositions(TCanvas* c){
    
  TPad * pd1 = (TPad*)c->GetPad(1);
  TPad * pd2 =(TPad*) c->GetPad(2);
  TPad * pd3 =(TPad*) c->GetPad(3);
  TPad * pd4 =(TPad*) c->GetPad(4);
  TPad * pd5 =(TPad*) c->GetPad(5);
  TPad * pd6 =(TPad*) c->GetPad(6);
  TPad * pd7 =(TPad*) c->GetPad(7);
  TPad * pd8 =(TPad*) c->GetPad(8);
  TPad * pd9 =(TPad*) c->GetPad(9);
  
  SetPadStyle(pd1);
  SetPadStyle(pd2);
  SetPadStyle(pd3);
  SetPadStyle(pd4);
  SetPadStyle(pd5);
  SetPadStyle(pd6);
  SetPadStyle(pd7);
  SetPadStyle(pd8);
  SetPadStyle(pd9);

  Double_t xl,xu,yl,yu;
  Double_t marginLeft=0.08;
  Double_t marginRight=0.04;
  Double_t marginTop=0.04;
  Double_t marginBottom=0.09;
  Double_t marginLeftForXAxis=0.02;
  Double_t marginBottomForYAxis=0.02;
  innerPadWidth=(1-marginLeft-marginRight)/3.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/3.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!


    // Bottom row
    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd7->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd7->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd7->SetRightMargin(0.);
    pd7->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd7->SetTopMargin(0.);

    pd7->Modified();
    pd7->Update();

    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd8->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2.*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd8->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd8->SetRightMargin(0.);
    pd8->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd8->SetTopMargin(0.);

    pd8->Modified();
    pd8->Update();
    
    pd9->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd9->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd9->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd9->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd9->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd9->SetTopMargin(0.);

    pd9->Modified();
    pd9->Update();

// Middle row
    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2*innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd4->SetRightMargin(0.);
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd4->SetTopMargin(0.);

    pd4->Modified();
    pd4->Update();

    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd5->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,2.*innerPadWidth+marginLeft,2*innerPadHeight+marginBottom);
    pd5->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd5->SetRightMargin(0.);
    pd5->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd5->SetTopMargin(0.);

    pd5->Modified();
    pd5->Update();
    
    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd6->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1.,2.*innerPadHeight+marginBottom);
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd6->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd6->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd6->SetTopMargin(0.);

    pd6->Modified();
    pd6->Update();

    // Top Row
    pd1->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd1->SetPad(0,2.*innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,1.);
    pd1->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd1->SetRightMargin(0.);
    pd1->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd1->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd1->Modified();
    pd1->Update();

    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,2.*innerPadWidth+marginLeft,1);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    //    pd2->SetLeftMargin(0.);
    pd2->SetRightMargin(0);
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd2->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(2.*innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd3->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd3->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();
    

}




TCanvas* CompareNSyieldPPtoPPb(Int_t binass){
  return ComparePPtoPPb(binass,0);
}

TCanvas* CompareNSsigmaPPtoPPb(Int_t binass){
  return  ComparePPtoPPb(binass,1);

}

TCanvas* ComparePedestalPPtoPPb(Int_t binass){
  return ComparePPtoPPb(binass,2);
}

TLatex *GetDRapForSystem(Int_t collSystem,Int_t identifier,Int_t includeDEta=0){
  
  TLatex *tlrap;
  Double_t x=0.015,y=0.28;
  if(ncolumns==3){
    x=0.015;// this and the following numbers does not make too much sense, they come just from an optimization
    if(includeDEta)x=0.060; // was 0.065
    if(collSystem==1&&includeDEta)x=0.026;// was 0.048
  }
  if(nrows==3){// these are hard coded number from an optimization
    if(includeDEta)y=0.25;
    else y=0.18;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.3850;// was 0.390 at round 1    
  }
  TString str;
  if(collSystem==0)str="|#it{y}^{D}_{cms}| < 0.5";
  if(collSystem==1)str="-0.96 < #it{y}^{D}_{cms} < 0.04";
  if(includeDEta==1){
    str.Append(", |#Delta#eta| < 1");
  }
  if(style==-1){
    if(collSystem==0)tlrap=new TLatex(0.24,0.75,"|#it{y}^{D}_{cms}| < 0.5");
    else if(collSystem==1)tlrap=new TLatex(0.24,0.75,"-0.96 < #it{y}^{D}_{cms} < 0.04");
    else return;
    tlrap->SetNDC();
    tlrap->SetTextFont(42);
    tlrap->SetTextSize(0.03);
  }
  else{
    if(collSystem==0)tlrap=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),str.Data()); 
    else if(collSystem==1)tlrap=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),str.Data()); 
    tlrap->SetNDC();
    tlrap->SetTextFont(43);
    tlrap->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  }  
  tlrap->SetTextAlign(11);
  tlrap->SetName(Form("tlrap_%d",identifier));
  return tlrap;
}

TLatex *GetDEtaD(Int_t identifier){
  
  TLatex *tlDEta;
  Double_t x=0.015,y=0.28;
  if(ncolumns==3){
    x=0.19;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.25;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.28;
  }
  
  tlDEta=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"|#Delta#eta| < 1"); 

  tlDEta->SetNDC();
  tlDEta->SetTextFont(43);
  tlDEta->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  
  tlDEta->SetTextAlign(12);
  tlDEta->SetName(Form("tlDEta_%d",identifier));
  return tlDEta;
}


TLatex *GetCollSystem(Int_t collSystem,Int_t identifier){
  Double_t x=0.007,y=0.390;
  if(ncolumns==3){
    x=0.007;
  }
  if(nrows==3){// these are hard coded number from an optimization; has precedence on ncolumns
    x=0.15;
    y=0.25;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.3850; // was 0.390 at round 1   
    if(collSystem==0){
      x=0.15;
    }
    else if(collSystem==1){
      x=0.15;
    }
  }

  TLatex *tlsystem;
  if(style==-1){
    if(collSystem==0)tlsystem=new TLatex(0.24,0.8,"pp, #sqrt{#it{s}} = 7 TeV");
    else if(collSystem==1)tlsystem=new TLatex(0.24,0.8,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    else return;
    tlsystem->SetNDC();
    tlsystem->SetTextFont(42);
    tlsystem->SetTextSize(0.03);
    tlsystem->SetTextAlign(12);
  }
  else{
    if(collSystem==0)tlsystem=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"pp, #sqrt{#it{s}} = 7 TeV"); 
    else if(collSystem==1)tlsystem=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV"); 
    tlsystem->SetNDC();
    tlsystem->SetTextFont(43);
    tlsystem->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
    tlsystem->SetTextAlign(21);
  }
  

  tlsystem->SetName(Form("tlSystem_%d",identifier));
  return tlsystem;
}

TLatex *GetALICEtext(Int_t identifier){
  TLatex *alice;
Double_t x=0.21,y=0.390;
  if(ncolumns==3){
    x=0.21;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.25;
    x=0.20;// draft 2 was not present -> above value 0.21
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.385;// draft 2 was 0.39
    x=0.20;// draft 2 was not present -> above value 0.21
  }

  if(style==-1){
    alice=new TLatex(0.75,0.85,"ALICE");
    alice->SetNDC();
    alice->SetTextFont(42);
    alice->SetTextSize(0.03);
    alice->SetTextAlign(11);
  }
  else{
    alice= new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"ALICE"); 
    alice->SetNDC();
    alice->SetTextFont(43);
    alice->SetTextSize(32*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);// draft 2 was: 28 *...
    alice->SetTextAlign(11);
  }


  //  TPaveText *alice = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  alice->SetName(Form("tlALICE_%d",identifier));
  return alice;
}

TLatex* GetAssocPtText(Int_t binassoc,Int_t identifier,Int_t addDEta=1){
  TLatex *tlasspt;
  Double_t x=0.035,y=0.35;
  if(ncolumns==3){
    //    x=0.035;
    x=0.15;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.215;
    //    if(addDEta==0)x=0.055;
    x=0.15;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.34;// was 0.35 in draft 2
    x=0.15;
  }

  if(style==-1){
    tlasspt=new TLatex(0.25,0.78,Form("%s, |#Delta#eta| < 1",strPtAssCanvas[binassoc].Data()));
    tlasspt->SetNDC();
    tlasspt->SetTextFont(42);
    tlasspt->SetTextSize(0.03);
  }
  else{
    TString strTot=strPtAssCanvas[binassoc];
    if(addDEta==1)strTot.Append(", |#Delta#eta| < 1");
    tlasspt= new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),strTot.Data()); 
    tlasspt->SetNDC();
    tlasspt->SetTextFont(43);
    tlasspt->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//  if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor) but see notes on top
    Printf("Height pad: %f, scal Height = %f",gPad->GetHNDC(),scaleHeightPads);
  }
  //  TPaveText *tlasspt = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  tlasspt->SetTextAlign(22);
  tlasspt->SetName(Form("tlAssocPt_%d",identifier));
  return tlasspt;

}

TLatex* GetTextSide(Int_t variable,Int_t identifier){
  TLatex *tlSide=new TLatex();    
  Double_t x=0.015,y=0.390;
  if(ncolumns==3){
    x=0.015;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.250;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.385;// draft 2 was 0.39, round 1 was 0.378
  }
  
  if(variable==2)return 0x0;
  else if(variable==0 || variable ==1)
    {
      if(style==-1){
	tlSide=new TLatex(0.25,0.85,"Near side");
	tlSide->SetNDC();
	tlSide->SetTextFont(42);
	tlSide->SetTextAlign(12);
	tlSide->SetTextSize(0.03*resizeTextFactor);
      }
      else{
	tlSide=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"Near side");
	tlSide->SetNDC();
	tlSide->SetTextAlign(11);
	tlSide->SetTextFont(43);
	tlSide->SetTextSize(32*innerPadHeight/referencePadHeight*resizeTextFactor);// draft 2 was 28*...
      }
    }
  else {
    if(style==-1){
      tlSide=new TLatex(0.25,0.85,"Away side");
      tlSide->SetNDC();
      tlSide->SetTextFont(42);
	tlSide->SetTextAlign(12);
      tlSide->SetTextSize(0.03*resizeTextFactor);
    }
    else{
      tlSide=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"Away side");
      tlSide->SetNDC();
      tlSide->SetTextFont(43);
      tlSide->SetTextAlign(11);
      tlSide->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    }
    

  }
  
  tlSide->SetName(Form("tlSideName_%d",identifier));
  return tlSide;
}
TLegend *GetLegendMCDataPoints(TH1D *hpp,TH1D *hpPb,Int_t identifier,TString strlegendHeader=""){
  TLegend * legend;
  Double_t xl=0.002,xr=0.3,yl=0.2,yh=0.3;
  if(ncolumns==3){
    xl=0.0015;
    xr=0.3;
  }
 
 if(nrows==3){// these are hard coded number from an optimization
   Int_t neffmod=CountNmodelsInLegend();
   if(!(strlegendHeader.IsNull()))neffmod++;
   //    yl=0.11;
    yh=0.20;
    yl=0.12+(3-neffmod)*0.07/3.;// for neffmod>3 an optimization might be needed
    if(neffmod==4&& strlegendHeader.EqualTo(" ")){// is just from optimization to move up the legends on the rightmost panels in the paper 
      yl+=0.07/5.;
      yh+=0.07/5.;
    }
  }
  if(nrows==2){// these are hard coded number from an optimization
   Int_t neffmod=CountNmodelsInLegend();   
   if(!(strlegendHeader.IsNull()))neffmod++;
   //    yl=0.2;
    yl=0.215+(3-neffmod)*0.1/3.;// for neffmod>3 an optimization might be needed
    yh=0.315;
    if(neffmod==5 && strlegendHeader.EqualTo(" ")){// is just from optimization to move up the legends on the rightmost panels in the paper 
      yl+=0.1/5.;
      yh+=0.1/5.;
    }
  }

  if(style==-1){
    legend=new TLegend(0.21,0.63,0.47,0.75);  
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
  }
  else{
    legend = new TLegend(xl/gPad->GetWNDC()+gPad->GetLeftMargin(),yl/gPad->GetHNDC()+gPad->GetBottomMargin(),xr/gPad->GetWNDC()+gPad->GetLeftMargin(),yh/gPad->GetHNDC()+gPad->GetBottomMargin());
    legend->SetTextFont(43);
    // Double_t hsiz=gPad->GetHNDC();
    //    Double_t wsiz=gPad->GetWNDC();
    //    Double_t textSiz=1./hsiz;
    //    if(wsiz>hsiz)textSiz=1./wsiz;
    // DID NOT UNDERSTAND WHY DIFFERNT SIZE IF TEXT SIZE IS IN PIXEL
    legend->SetTextSize(20*innerPadHeight/referencePadHeight*resizeTextFactor);//*(gPad->GetWNDC()-gPad->GetLeftMargin()-gPad->GetRightMargin())*textSiz);//*innerPadHeight/referencePadHeight*resizeTextFactor);//gPad->UtoPixel(gPad->GetWNDC()-gPad->GetLeftMargin()-gPad->GetRightMargin()));
			//*innerPadHeight/referencePadHeight*resizeTextFactor);//if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor); but see the notes on top
    legend->SetTextAlign(12);
  }
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  if(!(strlegendHeader.IsNull()))legend->SetHeader(strlegendHeader.Data());
  
  if(hpp)legend->AddEntry(hpp,"pp, #sqrt{#it{s}} = 7 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
  if(hpPb)legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, -0.96<#it{y}^{D}_{cms} < 0.04","lep");
  legend->SetName(Form("LegendDataAndMCPPandpPb_%d",identifier));
  return legend;
}

TLegend *GetLegendDataPoints(TH1D *hpp,TH1D *hpPb,Int_t identifier){
  TLegend * legend;
  Double_t xl=0.002,xr=0.15,yl=0.23,yh=0.3;
  if(ncolumns==3){
    xl=0.002;
    xr=0.15;
  }
  if(nrows==3){// these are hard coded number from an optimization
    yl=0.16;
    yh=0.2;
  }
  if(nrows==2){// these are hard coded number from an optimization
    yl=0.23;
    yh=0.3;
  }

  if(style==-1){
    legend=new TLegend(0.21,0.63,0.47,0.75);  
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
  }
  else{
    legend = new TLegend(0.01/gPad->GetWNDC()+gPad->GetLeftMargin(),0.185/gPad->GetHNDC()+gPad->GetBottomMargin(),0.18/gPad->GetWNDC()+gPad->GetLeftMargin(),0.312/gPad->GetHNDC()+gPad->GetBottomMargin());// draft 2 (2 lines only, rapidity on the same line also for p-Pb): 0.002/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.30/gPad->GetHNDC()+gPad->GetBottomMargin()
    legend->SetTextFont(43);
    legend->SetTextAlign(12);
    legend->SetTextSize(27*innerPadHeight/referencePadHeight*resizeTextFactor);// draft 2 was 20*... //if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor); but see the notes on top
  }
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  
  if(hpp)legend->AddEntry(hpp,"pp, #sqrt{#it{s}} = 7 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
  if(hpPb){// draft 2 was only:  legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04","lep");
    legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV,","lep");
    legend->AddEntry((TObject*)0,"-0.96 < #it{y}^{D}_{cms} < 0.04","");
  }
  legend->SetName(Form("LegendDataAndMCPPandpPb_%d",identifier));
  return legend;
}

TLegend *GetLegendDataPointsFake(TH1D *hpp,TH1D *hpPb,Int_t identifier){
  TLegend * legend;
  Double_t xl=0.002,xr=0.15,yl=0.23,yh=0.3;
  if(ncolumns==3){
    xl=0.002;
    xr=0.15;
  }
  if(nrows==3){// these are hard coded number from an optimization
    yl=0.16;
    yh=0.2;
  }
  if(nrows==2){// these are hard coded number from an optimization
    yl=0.23;
    yh=0.3;
  }

  if(style==-1){
    legend=new TLegend(0.21,0.63,0.47,0.75);  
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
  }
  else{
    legend = new TLegend(0.01/gPad->GetWNDC()+gPad->GetLeftMargin(),0.185/gPad->GetHNDC()+gPad->GetBottomMargin(),0.18/gPad->GetWNDC()+gPad->GetLeftMargin(),0.312/gPad->GetHNDC()+gPad->GetBottomMargin());// draft 2 (2 lines only, rapidity on the same line also for p-Pb): 0.002/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.30/gPad->GetHNDC()+gPad->GetBottomMargin()
    legend->SetTextFont(43);
    legend->SetTextAlign(12);
    legend->SetTextSize(27*innerPadHeight/referencePadHeight*resizeTextFactor);// draft 2 was 20*... //if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor); but see the notes on top
  }
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  
  if(hpp)legend->AddEntry(hpp,"","lep");
  if(hpPb){// draft 2 was only:  legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04","lep");
    legend->AddEntry(hpPb,"","lep");
    legend->AddEntry((TObject*)0,"","");
  }
  legend->SetName(Form("LegendDataAndMCPPandpPb_%d",identifier));
  return legend;
}

TH1D *GetAndPreparePP(Int_t binass,Int_t quantity,TGraphAsymmErrors *&gr){
  Printf("Opening file: %s", Form("%s/Trends_pp/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()));
  TFile *f=TFile::Open(Form("%s/Trends_pp/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
  TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  gr=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  gr->SetName(Form("%sPP",gr->GetName()));
  for(Int_t iPoint=0;iPoint<3;iPoint++) gr->SetPointError(iPoint,0.7*gr->GetErrorXlow(iPoint),0.7*gr->GetErrorXhigh(iPoint),gr->GetErrorYlow(iPoint),gr->GetErrorYhigh(iPoint));
  TH1D *hPP=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  hPP->SetName(Form("%sPP",hPP->GetName()));


  hPP->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
  hPP->SetYTitle(yaxisTitle[quantity].Data());
  if(style==-1){
    hPP->GetYaxis()->SetTitleSize(0.04);
    hPP->GetYaxis()->SetTitleOffset(1.2);
    hPP->GetYaxis()->SetLabelSize(0.04);
    hPP->GetXaxis()->SetTitleSize(0.04);
    hPP->GetXaxis()->SetLabelSize(0.04);
  }
  else {
    hPP->GetYaxis()->SetTitle("");      

    hPP->GetXaxis()->SetRangeUser(0,16.9);
    hPP->GetYaxis()->SetTitleFont(43);
    hPP->GetYaxis()->SetLabelFont(43);
    hPP->GetXaxis()->SetTitleFont(43);
    hPP->GetXaxis()->CenterTitle();
    hPP->GetYaxis()->CenterTitle();
    hPP->GetXaxis()->SetLabelFont(43);
    hPP->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPP->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPP->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPP->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPP->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPP->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
  }
  if(style<=0){
    hPP->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[binass][quantity],maxRangePPb[binass][quantity]));
  }


  hPP->SetLineColor(colSystem[0]);
  hPP->SetLineWidth(2);
  hPP->SetMarkerColor(colSystem[0]);
  hPP->SetMarkerStyle(markerStyle[0]);
  hPP->SetMarkerSize(markersize);

  gr->SetMarkerColor(colSystem[0]);
  gr->SetLineColor(colSystem[0]);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(markerStyle[0]);
  gr->SetMarkerSize(markersize);

  AdaptRangeHist(hPP,minptData,maxptData);
  AdaptRangeTGraph(gr,minptData,maxptData);
  
  return hPP;
}


TH1D *GetAndPreparePPb(Int_t binass,Int_t quantity,TGraphAsymmErrors *&gr, TGraphAsymmErrors *&grV2=0){


  TFile *f=TFile::Open(Form("%s/Trends_pPb/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
  TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  gr=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  gr->SetName(Form("%sPPb",gr->GetName()));
  for(Int_t iPoint=0;iPoint<2;iPoint++) gr->SetPointError(iPoint,0.7*gr->GetErrorXlow(iPoint),0.7*gr->GetErrorXhigh(iPoint),gr->GetErrorYlow(iPoint),gr->GetErrorYhigh(iPoint));
  if(plotv2unc==kTRUE) {
      grV2=(TGraphAsymmErrors*)c->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));
      grV2->SetName(Form("%sPPb",grV2->GetName()));
  }

  TH1D *hPPb=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  hPPb->SetName(Form("%sPPb",hPPb->GetName()));

  hPPb->SetLineColor(colSystem[1]);
  hPPb->SetLineWidth(2);
  hPPb->SetMarkerColor(colSystem[1]);
  hPPb->SetMarkerStyle(markerStyle[1]);
  hPPb->SetMarkerSize(markersize);
  gr->SetMarkerColor(colSystem[1]);
  gr->SetLineColor(colSystem[1]);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(markerStyle[1]);
  gr->SetMarkerSize(markersize);
  if(plotv2unc==kTRUE) {
    grV2->SetMarkerColor(kGreen-2);
    grV2->SetLineColor(kGreen-2);
    grV2->SetFillStyle(3001);
  }

  hPPb->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
  hPPb->SetYTitle(yaxisTitle[quantity].Data());
  if(style==-1){
    hPPb->GetYaxis()->SetTitleSize(0.04);
    hPPb->GetYaxis()->SetTitleOffset(1.2);
    hPPb->GetYaxis()->SetLabelSize(0.04);
    hPPb->GetXaxis()->SetTitleSize(0.04);
    hPPb->GetXaxis()->SetLabelSize(0.04);
  }
  else {
    hPPb->GetYaxis()->SetTitle("");      

    hPPb->GetXaxis()->SetRangeUser(0,16.9);
    hPPb->GetYaxis()->SetTitleFont(43);
    hPPb->GetYaxis()->SetLabelFont(43);
    hPPb->GetXaxis()->SetTitleFont(43);
    hPPb->GetXaxis()->CenterTitle();
    hPPb->GetYaxis()->CenterTitle();
    hPPb->GetXaxis()->SetLabelFont(43);
    hPPb->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
  }
  if(style<=0){
    hPPb->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[binass][quantity],maxRangePPb[binass][quantity]));
  }

  AdaptRangeHist(hPPb,minptData,maxptData);
  AdaptRangeTGraph(gr,minptData,maxptData);
  if(plotv2unc==kTRUE) AdaptRangeTGraph(grV2,minptData,maxptData);

  return hPPb;
}


TCanvas* ComparePPtoPPb(Int_t binass,Int_t quantity,TPad *pd=0x0,Int_t textlegendOptions=0){
//   Printf("Opening file: %s", Form("%s/Trends_pp/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()));
//   TFile *f=TFile::Open(Form("%s/Trends_pp/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
//   TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
//   TGraphAsymmErrors *grPP=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
//   grPP->SetName(Form("%sPP",grPP->GetName()));
//   TH1D *hPP=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
//   hPP->SetName(Form("%sPP",hPP->GetName()));
  TGraphAsymmErrors *grPP;
  TH1D *hPP=GetAndPreparePP(binass,quantity,grPP);

//     f=TFile::Open(Form("%s/Trends_pPb/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
//   c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
//   TGraphAsymmErrors *grPPb=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
//   grPPb->SetName(Form("%sPPb",grPPb->GetName()));
//   TH1D *hPPb=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
//   hPPb->SetName(Form("%sPPb",hPPb->GetName()));
  TGraphAsymmErrors *grPPb, *grPPbV2;
  TH1D *hPPb=GetAndPreparePPb(binass,quantity,grPPb,grPPbV2);


  TCanvas *cout=0x0;
  if(!pd){
    cout=CreateCanvasWithDefaultStyle(Form("%sComparisonPPtoPPbBinAss%d",strquantityFile[quantity].Data(),binass));
    cout->cd();
  }
  else {
    pd->cd();
  }
  //new TCanvas(Form("NSyieldComparisonBinAss%d",binass),Form("NSyieldComparisonBinAss%d",binass),800,800);

  TH2D *hDraw;
  //  hPP->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
  //  hPP->SetYTitle(yaxisTitle[quantity].Data());
  if(style==-1){
//     hPP->GetYaxis()->SetTitleSize(0.04);
//     hPP->GetYaxis()->SetTitleOffset(1.2);
//     hPP->GetYaxis()->SetLabelSize(0.04);
//     hPP->GetXaxis()->SetTitleSize(0.04);
//     hPP->GetXaxis()->SetLabelSize(0.04);
    hPP->Draw();
  }
  else {
    hDraw=new TH2D(Form("hDraw%d",10*quantity+binass),"",100,0,25,200,0,10);
    hDraw->GetYaxis()->SetTitle("");      
    //    hPP->GetYaxis()->SetTitle("");      

    hDraw->GetXaxis()->SetRangeUser(0,17.5);
    hDraw->GetYaxis()->SetTitleFont(43);
    hDraw->GetYaxis()->SetLabelFont(43);
    hDraw->GetXaxis()->SetTitleFont(43);
    hDraw->GetXaxis()->CenterTitle();
    hDraw->GetYaxis()->CenterTitle();
    hDraw->GetXaxis()->SetLabelFont(43);
    hDraw->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hDraw->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hDraw->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);

    if(textlegendOptions%10==2||textlegendOptions%10==3){
      hDraw->SetYTitle(yaxisTitle[quantity].Data());
    }
    else {
      hDraw->GetYaxis()->SetLabelSize(0);
    }
    if(textlegendOptions%10==1||textlegendOptions%10==3){
      hDraw->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
    }
    else {
      hDraw->GetXaxis()->SetLabelSize(0);
    }

    // Set same drawing settings also to hPP: not needed, just to have it also in that histo
//     hPP->GetXaxis()->SetRangeUser(0,16.9);
//     hPP->GetYaxis()->SetTitleFont(43);
//     hPP->GetYaxis()->SetLabelFont(43);
//     hPP->GetXaxis()->SetTitleFont(43);
//     hPP->GetXaxis()->CenterTitle();
//     hPP->GetYaxis()->CenterTitle();
//     hPP->GetXaxis()->SetLabelFont(43);
//     hPP->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
//     hPP->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
//     hPP->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
//     hPP->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
//     hPP->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
//     hPP->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    
    hDraw->Draw();
    hPP->Draw("same");
  }

  //  hPP->Draw("E0X0");// to avoid plotting the error along x
  if(style>0){
    hDraw->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[2][quantity],maxRangePPb[2][quantity]));
  }
  else {
    hPP->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[binass][quantity],maxRangePPb[binass][quantity]));
  }
 
//   hPP->SetLineColor(colSystem[0]);
//   hPP->SetLineWidth(2);
//   hPP->SetMarkerColor(colSystem[0]);
//   hPP->SetMarkerStyle(20);
//   hPP->SetMarkerSize(markersize);

//   grPP->SetMarkerColor(colSystem[0]);
//   grPP->SetLineColor(colSystem[0]);
//   grPP->SetLineWidth(2);
//   grPP->SetMarkerStyle(20);
//   grPP->SetMarkerSize(markersize);
  grPP->Draw("E2");
  
//  hPPb->Draw("same");

// p-Pb: displace the syst errors along x axis (Fabio)
  Double_t shift=0.3;
  Double_t binsx=grPPb->GetN();
  for(int i=0;i<binsx;i++) {
    Double_t x,y,ex1,ex2,ey1,ey2;
    
    grPPb->GetPoint(i,x,y);
    ex1=grPPb->GetErrorXlow(i); ex2=grPPb->GetErrorXhigh(i); ey1=grPPb->GetErrorYlow(i); ey2=grPPb->GetErrorYhigh(i);
    grPPb->SetPoint(i,x+shift,y);
    grPPb->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPbV2->GetPoint(i,x,y);
    ex1=grPPbV2->GetErrorXlow(i); ex2=grPPbV2->GetErrorXhigh(i); ey1=grPPbV2->GetErrorYlow(i); ey2=grPPbV2->GetErrorYhigh(i);
    grPPbV2->SetPoint(i,x+shift,y);
    grPPbV2->SetPointError(i,ex1,ex2,ey1,ey2);
  }

  grPPb->Draw("E2");
  grPPbV2->Draw("E2");

  TH1D* hPPbSuperimp = hPPb->Clone();
  hPPbSuperimp->SetMarkerStyle(25);
  hPPbSuperimp->SetMarkerColor(kRed+1);
//  hPPbSuperimp->Draw("same");

//Conversion of pPb TH1F to TGraph to displace the points along x axis (Fabio)
  TGraphAsymmErrors *gr_points_PPb, *gr_pointCount_PPb;
  if(!runonpPb2016) {
    ConvertTH1ToTGraphAsymmError(hPPb,gr_points_PPb,shift);
    ConvertTH1ToTGraphAsymmError(hPPb,gr_pointCount_PPb,shift);
  } else {
    ConvertTH1ToTGraphAsymmError2016(hPPb,gr_points_PPb,shift);
    ConvertTH1ToTGraphAsymmError2016(hPPb,gr_pointCount_PPb,shift);    
  }
  gr_points_PPb->SetLineColor(colSystem[1]);
  gr_points_PPb->SetLineWidth(2);
  gr_points_PPb->SetMarkerColor(colSystem[1]);
  gr_points_PPb->SetMarkerStyle(markerStyle[1]);
  gr_points_PPb->SetMarkerSize(markersize);
  gr_points_PPb->Draw("samePZ");
  
  gr_pointCount_PPb->SetLineColor(colSystem[1]);
  gr_pointCount_PPb->SetLineWidth(2);
  gr_pointCount_PPb->SetMarkerStyle(25);
  gr_pointCount_PPb->SetMarkerColor(kRed+1);
  gr_pointCount_PPb->SetMarkerSize(markersize);
  gr_pointCount_PPb->Draw("samePZ");

  if(quantity==1 && binass==2) {
    TLatex *tlDispl=new TLatex(0.27,0.86,"p-Pb points and error boxes");
    TLatex *tlDispl2=new TLatex(0.27,0.80,"shifted by #Delta#it{p}_{T} = +0.3 GeV/#it{c}");
    tlDispl->SetNDC();
    tlDispl->SetTextFont(42);
    tlDispl->SetTextSize(0.048);
    tlDispl->Draw();
    tlDispl2->SetNDC();
    tlDispl2->SetTextFont(42);
    tlDispl2->SetTextSize(0.048);
    tlDispl2->Draw();
  }

//   hPPb->SetLineColor(colSystem[1]);
//   hPPb->SetLineWidth(2);
//   hPPb->SetMarkerColor(colSystem[1]);
//   hPPb->SetMarkerStyle(21);
//   hPPb->SetMarkerSize(markersize);
//   grPPb->SetMarkerColor(colSystem[1]);
//   grPPb->SetLineColor(colSystem[1]);
//   grPPb->SetLineWidth(2);
//   grPPb->SetMarkerStyle(21);

//   grPPb->SetMarkerSize(markersize);

  if(style==-1){
    TLegend *legend=GetLegendDataPoints(hPP,hPPb,10*quantity+binass);
    legend->Draw();

    TLegend *legendSuperimp=GetLegendDataPointsFake(hPP,hPPbSuperimp,10*quantity+binass);
    legendSuperimp->Draw("same");
        
    TLatex *tlAssYieldPt=GetAssocPtText(binass,10*quantity+binass);
    tlAssYieldPt->Draw();

    TLatex *tlALICE=GetALICEtext(10*quantity+binass);
    tlALICE->Draw();
    
    
    if(quantity!=2){
      TLatex *tlSide=GetTextSide(quantity,10*quantity+binass);
      tlSide->Draw();
    }
  }
  else{

    if(textlegendOptions%100>=10){
      TLatex *tlAssYieldPt=GetAssocPtText(binass,10*quantity+binass);
      tlAssYieldPt->Draw();
    }
    if(textlegendOptions%1000>=100){
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10*quantity+binass);
      legend->Draw();

      TLegend *legendSuperimp=GetLegendDataPointsFake(hPP,hPPbSuperimp,10*quantity+binass);
      legendSuperimp->Draw("same");
    }
    if(textlegendOptions%10000>=1000){
      TLatex *tlALICE=GetALICEtext(10*quantity+binass);
      tlALICE->Draw();
    }
    if(textlegendOptions%100000>=10000){
      if(quantity!=2){
	TLatex *tlSide=GetTextSide(quantity,10*quantity+binass);
	tlSide->Draw();
      }
    }
  }
  
  return cout;
  
}


void InitMCobjects(){
  hMC=new TH1D*[nmodels];
  grMC=new TGraphAsymmErrors*[nmodels];
}


TCanvas* CompareDatatoModels(Int_t collsystem,Int_t binass,Int_t quantity,TPad *pd=0x0,Int_t textlegendOptions=0,Int_t drawMCasLines=0,TString legendHeader=""){
  Int_t system=collsystem;
  if(system==0){
    Printf("Opening file: %s", Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()));
  }
  else if(system==1){
    Printf("Opening file: %s", Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()));
  }
  else if(system==-1){
    Printf("Opening file: %s", Form("%s/Trends_%s/CanvasFinalTrendXXX_pthad%s.root",strFitResultPP.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()));
  }
  TH1D **hData=new TH1D*[2];
  TGraphAsymmErrors **grData=new  TGraphAsymmErrors*[2];

  TCanvas *cout=0x0;
  if(system==-1){
    cout=CreateCanvasWithDefaultStyle(Form("%sComparisonMCtoDataBothSystemsBinAss%d",strquantityFile[quantity].Data(),binass));
  }
  else{
    if(!pd){
      cout=CreateCanvasWithDefaultStyle(Form("%sComparisonMCto%sDataBinAss%d",strquantityFile[quantity].Data(),strSystem[system].Data(),binass));
      pd=(TPad*)cout->cd();
    }
    else{
      pd->cd();
    }
  }
  Int_t neffmod=CountNmodels();//InLegend();
  TLegend * legend;
  if(style==-1){
    legend= new TLegend(0.61,0.54,0.95,0.54+neffmod*0.05);
    //    legend= new TLegend(0.61,0.54,0.95,0.54+neffmod*0.08);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
  }
  else{
    if(textlegendOptions%1000>=100){
      TLegend *legend=GetLegendMCDataPoints(0x0,0x0,10*quantity+binass,legendHeader);
      legend->Draw();	
    }
  }
  //new TCanvas(Form("NSyieldComparisonBinAss%d",binass),Form("NSyieldComparisonBinAss%d",binass),800,800);
  TLatex *tlCollSystem=0x0;
  TLatex *tlDrap=0x0;
  
  if(collsystem==-1){
    
    system=0;   
    TFile *f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
     
    TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
    grData[0]=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
    grData[0]->SetName(Form("%sPP",grData[0]->GetName()));
    for(Int_t iPoint=0;iPoint<3;iPoint++) grData[0]->SetPointError(iPoint,0.7*grData[0]->GetErrorXlow(iPoint),0.7*grData[0]->GetErrorXhigh(iPoint),grData[0]->GetErrorYlow(iPoint),grData[0]->GetErrorYhigh(iPoint));
    hData[0]=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
    hData[0]->SetName(Form("%sPP",hData[0]->GetName()));
    
    AdaptRangeHist(hData[0],minptData,maxptData);
    AdaptRangeTGraph(grData[0],minptData,maxptData);

    system=1;
    f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
    c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
    grData[1]=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
    grData[1]->SetName(Form("%sPPb",grData[1]->GetName()));
    for(Int_t iPoint=0;iPoint<2;iPoint++) grData[1]->SetPointError(iPoint,0.7*grData[1]->GetErrorXlow(iPoint),0.7*grData[1]->GetErrorXhigh(iPoint),grData[1]->GetErrorYlow(iPoint),grData[1]->GetErrorYhigh(iPoint));
    hData[1]=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
    hData[1]->SetName(Form("%sPPb",hData[1]->GetName()));

    AdaptRangeData(hData[1],minptData,maxptData);
    AdaptRangeData(grData[1],minptData,maxptData);

    system=0;   
  }
  else{
    TFile *f;
    if(system==0)f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
    else if(system==1)f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
    
    TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
    grData[0]=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
    grData[0]->SetName(Form("%sPP",grData[0]->GetName()));

    Int_t iPointEnd=3;
    if(system==1) iPointEnd=2;
    for(Int_t iPoint=0;iPoint<iPointEnd;iPoint++) grData[0]->SetPointError(iPoint,0.7*grData[0]->GetErrorXlow(iPoint),0.7*grData[0]->GetErrorXhigh(iPoint),grData[0]->GetErrorYlow(iPoint),grData[0]->GetErrorYhigh(iPoint));
    hData[0]=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
    hData[0]->SetName(Form("%sPP",hData[0]->GetName()));

    AdaptRangeHist(hData[0],minptData,maxptData);
    AdaptRangeTGraph(grData[0],minptData,maxptData);

  }
    
  pd->cd();
  TH2D *hDraw;
  if(style==-1){
    hData[0]->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
    hData[0]->SetYTitle(yaxisTitle[quantity].Data());
    hData[0]->GetYaxis()->SetTitleSize(0.04);
    hData[0]->GetYaxis()->SetTitleOffset(1.2);
    hData[0]->GetYaxis()->SetLabelSize(0.04);
    hData[0]->GetXaxis()->SetTitleSize(0.04);
    hData[0]->GetXaxis()->SetLabelSize(0.04);
    hData[0]->Draw();
    //  hData[0]->Draw("E0X0");// to avoid plotting the error along x
    if(system==0)hData[0]->GetYaxis()->SetRangeUser(0,maxRangePP[binass][quantity]);
    if(system==1)hData[0]->GetYaxis()->SetRangeUser(0,maxRangePPb[binass][quantity]);
    
    hData[0]->SetLineColor(colSystem[system]);
    hData[0]->SetLineWidth(2);
    hData[0]->SetMarkerColor(colSystem[system]);
    hData[0]->SetMarkerStyle(markerStyle[system]);
    hData[0]->SetMarkerSize(markersize);
    
  }
  else{
    hDraw=new TH2D(Form("hDraw%d",10*quantity+binass),"",100,0,25,200,0,10);
    hDraw->GetYaxis()->SetTitle("");      

    hDraw->GetXaxis()->SetRangeUser(0,17.5);
    if(system==0){
      hData[0]->GetYaxis()->SetRangeUser(0,maxRangePP[2][quantity]);
      hDraw->GetYaxis()->SetRangeUser(0,maxRangePP[2][quantity]);
    }
    if(system==1){
      hData[0]->GetYaxis()->SetRangeUser(0,maxRangePPb[2][quantity]);
      hDraw->GetYaxis()->SetRangeUser(0,maxRangePPb[2][quantity]);
    }

    //    hData[0]->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[0][quantity],maxRangePPb[0][quantity]));

    //    hDraw->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[0][quantity],maxRangePPb[0][quantity]));
    hDraw->GetYaxis()->SetTitleFont(43);
    hDraw->GetYaxis()->SetLabelFont(43);
    hDraw->GetXaxis()->SetTitleFont(43);
    hDraw->GetXaxis()->CenterTitle();
    hDraw->GetYaxis()->CenterTitle();
    hDraw->GetXaxis()->SetLabelFont(43);
    hDraw->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetYaxis()->SetTitleOffset(ytitleoffset*innerPadHeight/referencePadHeight*resizeTextFactor);//ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hDraw->GetXaxis()->SetTitleOffset(xtitleoffset*innerPadHeight/referencePadHeight*resizeTextFactor);//xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hDraw->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);

    if(textlegendOptions%10==2||textlegendOptions%10==3){
      hDraw->SetYTitle(yaxisTitle[quantity].Data());
    }
    else {
      hDraw->GetYaxis()->SetLabelSize(0);
    }
    if(textlegendOptions%10==1||textlegendOptions%10==3){
      hDraw->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
    }
    else {
      hDraw->GetXaxis()->SetLabelSize(0);
    }
    hDraw->Draw();
    hData[0]->Draw("same");
  }
  hData[0]->SetMarkerSize(markersize);
  hData[0]->SetMarkerStyle(markerStyle[system]);
  hData[0]->SetLineColor(colSystem[system]);
  hData[0]->SetMarkerColor(colSystem[system]);
    grData[0]->SetMarkerColor(colSystem[system]);
    grData[0]->SetLineColor(colSystem[system]);
    grData[0]->SetLineWidth(2);
    grData[0]->SetMarkerStyle(markerStyle[system]);
    grData[0]->SetMarkerSize(markersize);
    grData[0]->Draw("E2");
   
    if(collsystem==-1){    
      hData[1]->Draw("same");
      hData[1]->SetLineColor(colSystem[1]);
      hData[1]->SetLineWidth(2);
      hData[1]->SetMarkerColor(colSystem[1]);
      hData[1]->SetMarkerStyle(markerStyle[1]);
      hData[1]->SetMarkerSize(markersize);
      grData[1]->SetMarkerColor(colSystem[1]);
      grData[1]->SetLineColor(colSystem[1]);
      grData[1]->SetLineWidth(2);
      grData[1]->SetMarkerStyle(markerStyle[1]);
      grData[1]->SetMarkerSize(markersize);
      grData[1]->Draw("E2");
    }



    if(collsystem==-1){
      legend->AddEntry(hData[0],"pp, #sqrt{#it{s}} = 7 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
      legend->AddEntry(hData[1],"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, -0.96 < #it{y}^{D}_{cms} < 0.04","lep");      
    }
    else {    
      //      legend->AddEntry(hData[0],"data","lep");
      if(textlegendOptions%1000000>=100000)tlCollSystem=GetCollSystem(system,10*quantity+binass);
      if(textlegendOptions%10000000>=1000000)tlDrap=GetDRapForSystem(system,10*quantity+binass);
    }
    
    
    // INITIATE MC HISTOS AND GRAPHS
    InitMCobjects();
    // NOW LOOP OVER MODELS
    for(Int_t kmc=0;kmc<nmodels;kmc++){
      if(!includemodel[kmc])continue;
      if(quantity==2&&strModelDir[kmc].Contains("POWHEG"))continue;
      f=TFile::Open(Form("%s/Trends_%s/%s/CanvasFinalTrend%s_pthad%s.root",strFitResultMC[system].Data(),strSystem[system].Data(),strModelDir[kmc].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
      c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
      grMC[kmc]=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
      grMC[kmc]->SetName(Form("%s%s",grMC[kmc]->GetName(),strModelDir[kmc].Data()));
      hMC[kmc]=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
      hMC[kmc]->SetName(Form("%s%s",hMC[kmc]->GetName(),strModelDir[kmc].Data()));

      AdaptRangeHist(hMC[kmc],minptMC,maxptMC);
      AdaptRangeTGraph(grMC[kmc],minptMC,maxptMC);

      pd->cd();
      if(drawMCasLines==2){
	hMC[kmc]->Draw("C X0 same");
      }
      else{
	hMC[kmc]->Draw("same");
      }
      hMC[kmc]->SetLineColor(modelColors[kmc]);
      hMC[kmc]->SetLineWidth(2);
      if(drawMCasLines==1){
	hMC[kmc]->SetLineStyle(2);
      }
      hMC[kmc]->SetMarkerColor(modelColors[kmc]);
      hMC[kmc]->SetMarkerStyle(modelMarkerStyle[kmc]);
      hMC[kmc]->SetMarkerSize(markersizeMC);
      grMC[kmc]->SetMarkerColor(modelColors[kmc]);
      grMC[kmc]->SetLineColor(modelColors[kmc]);
      grMC[kmc]->SetLineWidth(2);
      grMC[kmc]->SetMarkerStyle(modelMarkerStyle[kmc]);
      grMC[kmc]->SetFillStyle(3001+kmc);
      grMC[kmc]->SetFillColor(modelColors[kmc]);
      grMC[kmc]->SetMarkerSize(markersizeMC);
      if(drawSystMC){
	if(drawMCasLines==2){
	  grMC[kmc]->Draw("C0");
	}
	else {
	  grMC[kmc]->Draw("E5");
	}
      }
      
      //      if(legend)legend->AddEntry(hMC[kmc],Form("%s",strModelDir[kmc].Data()),"lep");  
      if(legend&&includeinlegend[kmc]){
	if(kmc==6){
	  strModelDirLeg[6].ReplaceAll(" EPS09","");
	  legend->AddEntry(hMC[kmc],Form("%s",strModelDirLeg[kmc].Data()),"lep");  
	  legend->AddEntry((TObject*)0,"with EPS09 nPDF","");  
	}
	else {
	  legend->AddEntry(hMC[kmc],Form("%s",strModelDirLeg[kmc].Data()),"lep");  
	}
      }
      

    }
    
    if(legend){
      legend->Draw();
    }

    if(style==-1){    
      TLatex *tlALICE=new TLatex(0.68,0.86,"ALICE");
      tlALICE->SetNDC();
      tlALICE->SetTextSize(0.04);
      tlALICE->Draw();
      if(quantity!=2){
	
	TLatex *tlSide;
	if(quantity==0||quantity==1)tlSide=new TLatex(0.24,0.86,"Near side");
	else tlSide=new TLatex(0.75,0.9,"Away side");
	tlSide->SetNDC();
	tlSide->SetTextSize(0.04);
	tlSide->Draw();
	
      }
      TLatex *tlAssYieldPt=new TLatex(0.24,0.7,Form("%s, |#Delta#eta| < 1",strPtAssCanvas[binass].Data()));
      tlAssYieldPt->SetNDC();
      tlAssYieldPt->SetTextSize(0.03);
      tlAssYieldPt->Draw();
    }
    else{
      
      if(textlegendOptions%100>=10){
	TLatex *tlAssYieldPt=GetAssocPtText(binass,10*quantity+binass);
	tlAssYieldPt->Draw();
      }
      if(textlegendOptions%10000>=1000){
	TLatex *tlALICE=GetALICEtext(10*quantity+binass);
	tlALICE->Draw();
      }
      if(textlegendOptions%100000>=10000){
	if(quantity!=2){
	  TLatex *tlSide=GetTextSide(quantity,10*quantity+binass);
	  tlSide->Draw();
	}
      }
    }
    

    hData[0]->Draw("same");    
    grData[0]->Draw("E2");   

    TH1D* hSuperimp = (TH1D*)hData[0]->Clone();
    hSuperimp->SetMarkerColor(kRed+1);
    hSuperimp->SetMarkerStyle(25);
    if(system==1) hSuperimp->Draw("same");        


    if(collsystem!=-1){
      if(tlCollSystem!=0x0){
	if(useLegendForData){
	  TLegend *legData=new TLegend(gPad->GetLeftMargin()+0.02,tlCollSystem->GetY()-0.02,0.95-gPad->GetRightMargin(),tlCollSystem->GetY()+tlCollSystem->GetYsize()-0.02,"");//0.388/gPad->GetHNDC()+gPad->GetBottomMargin());
	  //0.378/gPad->GetHNDC()+gPad->GetBottomMargin()
	  //tlCollSystem->GetY()+tlCollSystem->GetYsize(),"");
	  //tlCollSystem->GetX(),tlCollSystem->GetY(),tlCollSystem->GetX()+tlCollSystem->GetXsize(),tlCollSystem->GetY()+tlCollSystem->GetYsize());
	  legData->AddEntry(hData[0],tlCollSystem->GetTitle(),"lp");
	  //gPad->GetWNDC()
	  legData->SetTextAlign(12);
	  legData->SetTextFont(43);
	  legData->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
	  legData->Draw();
	  if(system==1) {
	  TLegend *legDataSuperimp=new TLegend(gPad->GetLeftMargin()+0.02,tlCollSystem->GetY()-0.02,0.95-gPad->GetRightMargin(),tlCollSystem->GetY()+tlCollSystem->GetYsize()-0.02,"");
	    legDataSuperimp->AddEntry(hSuperimp,"","lp");
	    legDataSuperimp->SetFillStyle(0);
	    legDataSuperimp->SetTextAlign(12);
	    legDataSuperimp->SetTextFont(43);
	    legDataSuperimp->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
	    legDataSuperimp->Draw("same");
 	  }
	}
	else tlCollSystem->Draw();
      }
      if(tlDrap!=0x0){
	tlDrap->Draw();
      }
    }
    else{
      hData[1]->Draw("same");    
      grData[1]->Draw("E2");   
    }
 
    return cout;
    
}




void CompareFitResultsPPtoPPbUniqueCanvas(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPvsPPbFitResultsFinalPaperStyle","pp vs.pPb fit results ",1800./1200.*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(3,2,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set2x3PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);
  }
  else{
    Printf("NOT READY YET");
    return;
  }  

  Int_t orderAssoc[3]={2,0,1};  
  for(Int_t jp=0;jp<=2;jp++){// First loop for NS yield
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      ComparePPtoPPb(orderAssoc[jp],0,pd,100+needTitle);    
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      ComparePPtoPPb(orderAssoc[jp],0,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderAssoc[jp]);
      tlSide->Draw();
    }
    else{
      ComparePPtoPPb(orderAssoc[jp],0,pd,0+needTitle);    
    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp]);
    tlAssYieldPt->Draw();
    if(jp==2){
      TLatex *tlPythiaStatement=new TLatex(0.11*gPad->GetWNDC()+gPad->GetLeftMargin(),0.27/gPad->GetHNDC()+gPad->GetBottomMargin(),"<7% variation expected from energy and");
      TLatex *tlPythiaStatementTwo=new TLatex(0.05*gPad->GetWNDC()+gPad->GetLeftMargin(),0.24/gPad->GetHNDC()+gPad->GetBottomMargin(),"rapidity difference (PYTHIA6, Perugia 2011)");//"#frac{Assoc. yield (#sqrt{#it{s}}=7 TeV, |#it{y}^{D}_{cms}|<0.5)}{Assoc. yield(#sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04)}~1.07");
      tlPythiaStatement->SetTextAlign(11);
      tlPythiaStatement->SetTextFont(43);
      tlPythiaStatement->SetNDC();
      tlPythiaStatement->SetTextSize(20*innerPadHeight/referencePadHeight*resizeTextFactor);
      tlPythiaStatement->Draw();
      tlPythiaStatementTwo->SetTextAlign(11);
      tlPythiaStatementTwo->SetTextFont(43);
      tlPythiaStatementTwo->SetNDC();
      tlPythiaStatementTwo->SetTextSize(20*innerPadHeight/referencePadHeight*resizeTextFactor);
      tlPythiaStatementTwo->Draw();
    }
 
    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=2;jp++){// second loop --> sigma
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+4);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    ComparePPtoPPb(orderAssoc[jp],1,pd,needTitle); 
    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binass; Sigma --> 1 ; binass == orderAssoc[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderAssoc[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binass; NS:
      TLatex *tlALICE=GetALICEtext(10+orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderAssoc[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],10+orderAssoc[jp]);
    tlAssYieldPt->Draw();
    */

    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=6;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  cFinalPaperStyle->SaveAs("ComparePPtoPPbFitResults.root");
  cFinalPaperStyle->SaveAs("ComparePPtoPPbFitResults.eps");
  cFinalPaperStyle->SaveAs("ComparePPtoPPbFitResults.png");
  cFinalPaperStyle->SaveAs("ComparePPtoPPbFitResults.pdf");

  return;

}

void CompareFitResultsPPtoPPb(){
  
  TCanvas *c=CompareNSyieldPPtoPPb(0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareNSyieldPPtoPPb(1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareNSyieldPPtoPPb(2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

  c=CompareNSsigmaPPtoPPb(0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareNSsigmaPPtoPPb(1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareNSsigmaPPtoPPb(2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

  c=ComparePedestalPPtoPPb(0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=ComparePedestalPPtoPPb(1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=ComparePedestalPPtoPPb(2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

}


void CompareFitResultsPPDataToMC(){
  
  TCanvas *c=CompareDatatoModels(0,0,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(0,0,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(0,0,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));


  c=CompareDatatoModels(0,1,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(0,1,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(0,1,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));


  c=CompareDatatoModels(0,2,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(0,2,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(0,2,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

}



void CompareFitResultsPPbDataToMC(){
  
  TCanvas *c=CompareDatatoModels(1,0,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,0,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,0,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));


  c=CompareDatatoModels(1,1,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,1,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,1,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));


  c=CompareDatatoModels(1,2,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,2,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,2,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

}



void CompareFitResultsDataBothSystemToMCPP(){
  
  TCanvas *c=CompareDatatoModels(-1,0,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(-1,0,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(-1,0,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));


  c=CompareDatatoModels(-1,1,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(-1,1,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(-1,1,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));


  c=CompareDatatoModels(-1,2,0);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(-1,2,1);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(-1,2,2);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

}





void CompareFitResultsPPtoMCUniqueCanvas(){
  gStyle->SetOptStat(0000);
  Init3x3Settings();
  TCanvas *cFinalPaperStyle;
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPvsMCFitResultsFinalPaperStyle","pp vs. MC fit results ",canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(3,3,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set3x3PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);
  }
  else{
    Printf("NOT READY YET");
    return;
  }  

  Bool_t includeinlegendOrig[nmodels];
  for(Int_t k=0;k<nmodels;k++)includeinlegendOrig[k]=includeinlegend[k];
  Int_t orderAssoc[3]={2,0,1};  
  for(Int_t jp=0;jp<=2;jp++){

    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      includeinlegend[3]=kFALSE;
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(0,orderAssoc[jp],0,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 7 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[3]=includeinlegendOrig[3];
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      CompareDatatoModels(0,orderAssoc[jp],0,pd,needTitle);//+1000000);    
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderAssoc[jp]);
      tlSide->Draw();
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      CompareDatatoModels(0,orderAssoc[jp],0,pd,needTitle+100,0,"");// " "); the latter " " needed for counting lines properly    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];
    }
    else{
      CompareDatatoModels(0,orderAssoc[jp],0,pd,needTitle);    
    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp],0);
    tlAssYieldPt->Draw();
 
    if(jp==2){
      TLatex *tlDeta=GetDRapForSystem(0,orderAssoc[jp],1);
      tlDeta->Draw();
    }
      gStyle->SetOptStat(0000);   
  }
  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+4);
    SetPadStyle(pd);
    pd->cd();
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderAssoc[jp],1,pd,needTitle);    

    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binass; Sigma --> 1 ; binass == orderAssoc[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderAssoc[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binass; NS:
      TLatex *tlALICE=GetALICEtext(10+orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderAssoc[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],10+orderAssoc[jp]);
    tlAssYieldPt->Draw();
    */

    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+7);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderAssoc[jp],2,pd,needTitle);    
    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=9;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  if(drawSystMC){
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults.pdf");
  }
  else{
  cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults.pdf");
  }
  return;

}




void CompareFitResultsPPtoPpbAndMCUniqueCanvas(){
  gStyle->SetOptStat(0000);
  Init3x3Settings();
  TCanvas *cFinalPaperStyle;
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPvsMCFitResultsFinalPaperStyle","pp vs. MC fit results ",canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(3,3,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set3x3PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);
  }
  else{
    Printf("NOT READY YET");
    return;
  }  

  Int_t orderAssoc[3]={2,0,1};  
  for(Int_t jp=0;jp<=2;jp++){// first loop -> assoc yield
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      includeinlegend[3]=kFALSE;
      includeinlegend[4]=kFALSE;
      CompareDatatoModels(0,orderAssoc[jp],0,pd,needTitle,1+100);    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      CompareDatatoModels(0,orderAssoc[jp],0,pd,needTitle,1);//+1000000);          
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderAssoc[jp]);
      tlSide->Draw();
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      includeinlegend[3]=kTRUE;
      includeinlegend[4]=kFALSE;
      CompareDatatoModels(0,orderAssoc[jp],0,pd,needTitle+100,1);    
    }
    else{
      CompareDatatoModels(0,orderAssoc[jp],0,pd,needTitle,1);    
    }

    TGraphAsymmErrors *grPPb;
    TH1D *hPPb=GetAndPreparePPb(orderAssoc[jp],0,grPPb);
    Printf("PPb histo and graph: %p, %p",hPPb,grPPb);
    pd->cd();
    hPPb->Draw("same");
    grPPb->Draw("E2");

    if(jp==1){
      TH1D *hPP=(TH1D*)pd->FindObject("FinalTrendNSYieldPP");
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,orderAssoc[jp]);
      legend->SetX1(0.02/gPad->GetWNDC()+gPad->GetLeftMargin());
      legend->SetX2(0.16/gPad->GetWNDC()+gPad->GetLeftMargin());
      legend->SetY1(0.11/gPad->GetHNDC()+gPad->GetBottomMargin());
      legend->SetY2(0.2/gPad->GetHNDC()+gPad->GetBottomMargin());
      legend->Draw();
    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp],0);
    tlAssYieldPt->Draw();
 
    if(jp==2){
      TLatex *tlDeta=GetDRapForSystem(0,orderAssoc[jp],1);
      tlDeta->Draw();
    }
      gStyle->SetOptStat(0000);   
  }
  
  for(Int_t jp=0;jp<=2;jp++){// second loop -> sigma
    Int_t needTitle=0;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+4);
    SetPadStyle(pd);
    pd->cd();
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderAssoc[jp],1,pd,needTitle,1);    

  
    TGraphAsymmErrors *grsigmaPPb;
    TH1D *hPPbsigma=GetAndPreparePPb(orderAssoc[jp],1,grsigmaPPb);
    pd->cd();
    hPPbsigma->Draw("same");
    grsigmaPPb->Draw("E2");

    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binass; Sigma --> 1 ; binass == orderAssoc[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderAssoc[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binass; NS:
      TLatex *tlALICE=GetALICEtext(10+orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderAssoc[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],10+orderAssoc[jp]);
    tlAssYieldPt->Draw();
    */

    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=2;jp++){ // third loop -> baseline
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+7);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderAssoc[jp],2,pd,needTitle,1);    
    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=9;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  cFinalPaperStyle->SaveAs("ComparePPtoPPbandPPMCFitResults.root");
  cFinalPaperStyle->SaveAs("ComparePPtoPPbandPPMCFitResults.eps");
  cFinalPaperStyle->SaveAs("ComparePPtoPPbandPPMCFitResults.png");
  cFinalPaperStyle->SaveAs("ComparePPtoPPbandPPMCFitResults.pdf");

  return;

}








void CompareFitResultsPPbtoMCUniqueCanvas(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
  xtitleoffset=2.5;
  strModelDirLeg[7]="EPOS 3.117, p-Pb simulation";
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPbvsMCFitResultsFinalPaperStyle","pp vs. MC fit results ",1800./1200.*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(3,2,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set2x3PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);
  }
  else{
    Printf("NOT READY YET");
    return;
  }  

  Bool_t includeinlegendOrig[nmodels];
  for(Int_t k=0;k<nmodels;k++)includeinlegendOrig[k]=includeinlegend[k];
  Int_t orderAssoc[3]={2,0,1};  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      includeinlegend[3]=kFALSE;
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(1,orderAssoc[jp],0,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 5.02 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[3]=includeinlegendOrig[3];
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];

    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      CompareDatatoModels(1,orderAssoc[jp],0,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderAssoc[jp]);
      tlSide->Draw();
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      CompareDatatoModels(1,orderAssoc[jp],0,pd,needTitle+100,0,"");//" ");// the latter empty space needed for counting properly lines    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];
      TLatex *tlDeta=GetDRapForSystem(1,orderAssoc[jp],1);
      tlDeta->Draw();
    
    }
    else{
      CompareDatatoModels(1,orderAssoc[jp],0,pd,needTitle);    
    }

    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp],0);
    tlAssYieldPt->Draw();
 
    gStyle->SetOptStat(0000);   
  }
  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+4);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    CompareDatatoModels(1,orderAssoc[jp],1,pd,needTitle);    

    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binass; Sigma --> 1 ; binass == orderAssoc[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderAssoc[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binass; NS:
      TLatex *tlALICE=GetALICEtext(10+orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderAssoc[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],10+orderAssoc[jp]);
    tlAssYieldPt->Draw();
    */

    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=6;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  if(drawSystMC){
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults.pdf");
  }
  else {
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults.pdf");
  }
  return;

}



void CompareFitResultsPPtoMCUniqueCanvasAwaySide(){
  gStyle->SetOptStat(0000);
  Init3x3Settings();
  TCanvas *cFinalPaperStyle;
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPvsMCFitResultsFinalPaperStyleAS","pp vs. MC fit results Away Side ",canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(3,3,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set3x3PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);
  }
  else{
    Printf("NOT READY YET");
    return;
  }  

  Bool_t includeinlegendOrig[nmodels];
  for(Int_t k=0;k<nmodels;k++)includeinlegendOrig[k]=includeinlegend[k];
  Int_t orderAssoc[3]={2,0,1};  
  for(Int_t jp=0;jp<=2;jp++){

    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      includeinlegend[3]=kFALSE;
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(0,orderAssoc[jp],3,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 7 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[3]=includeinlegendOrig[3];
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      CompareDatatoModels(0,orderAssoc[jp],3,pd,needTitle);//+1000000);    
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(3,orderAssoc[jp]);
      tlSide->Draw();
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      CompareDatatoModels(0,orderAssoc[jp],3,pd,needTitle+100,0,"");// " "); the latter " " needed for counting lines properli    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];
    }
    else{
      CompareDatatoModels(0,orderAssoc[jp],3,pd,needTitle);    
    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp],0);
    tlAssYieldPt->Draw();
 
    if(jp==2){
      TLatex *tlDeta=GetDRapForSystem(0,orderAssoc[jp],1);
      tlDeta->Draw();
    }
      gStyle->SetOptStat(0000);   
  }
  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+4);
    SetPadStyle(pd);
    pd->cd();
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderAssoc[jp],4,pd,needTitle);    

    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binass; Sigma --> 1 ; binass == orderAssoc[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderAssoc[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binass; NS:
      TLatex *tlALICE=GetALICEtext(10+orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderAssoc[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],10+orderAssoc[jp]);
    tlAssYieldPt->Draw();
    */

    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+7);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderAssoc[jp],2,pd,needTitle);    
    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=9;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  if(drawSystMC){
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS.pdf");
  }
  else{
  cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS.pdf");
  }
  return;

}


void CompareFitResultsPPbtoMCUniqueCanvasAwaySide(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
  xtitleoffset=2.5;
  strModelDirLeg[7]="EPOS 3.117, p-Pb simulation";
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPbvsMCFitResultsFinalPaperStyleAS","pp vs. MC fit results Away Side ",1800./1200.*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(3,2,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set2x3PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);
  }
  else{
    Printf("NOT READY YET");
    return;
  }  

  Bool_t includeinlegendOrig[nmodels];
  for(Int_t k=0;k<nmodels;k++)includeinlegendOrig[k]=includeinlegend[k];
  Int_t orderAssoc[3]={2,0,1};  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      includeinlegend[3]=kFALSE;
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(1,orderAssoc[jp],3,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 5.02 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[3]=includeinlegendOrig[3];
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];

    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      CompareDatatoModels(1,orderAssoc[jp],3,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(3,orderAssoc[jp]);
      tlSide->Draw();
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      CompareDatatoModels(1,orderAssoc[jp],3,pd,needTitle+100,0,"");//" ");// the latter empty space needed for counting properly lines    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];
      TLatex *tlDeta=GetDRapForSystem(1,orderAssoc[jp],1);
      tlDeta->Draw();
    
    }
    else{
      CompareDatatoModels(1,orderAssoc[jp],3,pd,needTitle);    
    }

    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp],0);
    tlAssYieldPt->Draw();
 
    gStyle->SetOptStat(0000);   
  }
  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+4);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    CompareDatatoModels(1,orderAssoc[jp],4,pd,needTitle);    

    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binass; Sigma --> 1 ; binass == orderAssoc[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderAssoc[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binass; NS:
      TLatex *tlALICE=GetALICEtext(10+orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderAssoc[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],10+orderAssoc[jp]);
    tlAssYieldPt->Draw();
    */

    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=6;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  if(drawSystMC){
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS.pdf");
  }
  else {
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS.pdf");
  }
  return;

}


