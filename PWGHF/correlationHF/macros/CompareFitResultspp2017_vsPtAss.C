TString strPtAss[3]={"0.3to1.0","1.0to2.0","2.0to3.0"};
TString strPtDCanvas[4]={"3 < #it{p}_{T}^{D} < 5 GeV/#it{c}","5 < #it{p}_{T}^{D} < 8 GeV/#it{c}","8 < #it{p}_{T}^{D} < 16 GeV/#it{c}","16 < #it{p}_{T}^{D} < 24 GeV/#it{c}"};
TString strSystem[2]={"pp","pPb"};
Color_t colSystem[2]={kBlack,kRed};
Int_t markerStyle[2]={20,21};
Bool_t useLegendForData=kFALSE;
Bool_t plotv2unc=kFALSE;
TString strFitResultPP=""; //  "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";
TString strFitResultPPb=""; //  "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";
Double_t canvasheight=1000;
Double_t resizeTextFactor=1.2;//1.*canvasheight/1200.; // size was tuned for canvasheight =1800. 
Double_t referencePadHeight=0.48; // Do not touch unless the canvas division is changed from 2x3, change canvasheight and resizeTextFactor instead
Double_t innerPadHeight;// not touch, set internally
Double_t innerPadWidth;// not touch, set internally
Int_t style=1;
Double_t scaleHeightPads=1;// do not touch this, it is regulated automatically in the macro
Double_t scaleWidthPads=1;// do not touch this, it is regulated automatically in the macro
TString strFitResultMC[2]={"",""};
//TString strFitResultMC[2]={"/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/MCTemplates/Templates_pp_12May15/FitResults/",
//			   "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May19UseScriptPWGHF/MCTemplates/Templates_pPb_12May15/FitResults/"};

/********************************************************************/
// TO REMOVE/INCLUDE POWHEG BASELINE, COMMENT/UNCOMMENT THE LINE:  //
// "if(quantity==2&&strModelDir[kmc].Contains("POWHEG"))continue;" //
/********************************************************************/

TString strquantityFile[5]={"NSYield","NSSigma","Pedestal","ASYield","ASSigma"};
Double_t maxRangePP[8][5]={   //used in L1403. MODIFY THESE!!!
         {1.7,0.64,2.33,2.1,1.09},
         {1.7,0.64,2.33,2.1,1.09},
         {1.7,0.64,2.33,2.1,1.09},
         {1.7,0.64,2.33,2.1,1.09},
         {1.35,0.64,2.33,2,1.28},
         {1.35,0.64,2.33,2,1.28},
         {1.35,0.64,2.33,2,1.28},
         {1.35,0.64,2.33,2,1.28}};

Double_t maxRangePP_Ratio[4][5]={   //used in ratios [quantity: NSY,NSs,bas,ASy,ASs][binD]
         {1.7,2.1,2,2.8,1.7},
         {1.7,2.1,2,2.8,1.7},
         {1.7,2.1,2,2.8,1.7},
         {1.7,2.1,2,2.8,1.7}};

Double_t minRangePP_Ratio[4][5]={   //used in ratios [quantity][binass]
         {0.15,0.4,0.5,0.4,0.5},
         {0.15,0.4,0.5,0.4,0.5},
         {0.15,0.4,0.5,0.4,0.5},
         {0.15,0.4,0.5,0.4,0.5}};

Double_t maxRangePPb[6][5]={
         {3.3,0.64,13,4.3,1.28},
         {3.3,0.64,13,4.3,1.28},
         {3.3,0.64,13,4.3,1.28},
         {1.35,0.64,13,2,1.28},
         {1.35,0.64,13,2,1.28},
         {1.35,0.64,13,2,1.28}};

Double_t maxRangePPbSinglePanel[6][5]={
         {3.2,0.8,9.0,4.3,1.48},
         {1.6,0.8,6.0,2.5,1.48},
         {2.1,0.8,3.0,2.5,1.48},
         {1.5,0.8,2.0,1.5,1.48},
         {0.6,0.8,0.5,1.0,1.48},
         {0.6,0.8,0.3,1.0,1.48}};

Double_t minptMC=0.,maxptMC=99.;
Double_t minptData=0.5,maxptData=3.5; //you have to be tighter than real limits...
Int_t ncolumns=4;
Int_t nrows=3;
Double_t ytitleoffset=2.9,xtitleoffset=2.2;
Bool_t skip3to5=kTRUE;
Double_t markersize=1.2;
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
      h->SetBinContent(i,-1);
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

TString yaxisTitle[5]={"Associated yield","Peak width (rad)","Baseline (rad^{-1})","Associated yield","Peak width (rad)"};
Double_t leftMarginCanvas=0.17;
Double_t rightMarginCanvas=0.055;
Double_t bottomMarginCanvas=0.13;
Double_t topMarginCanvas=0.07;
const Int_t nmodels=8;
Bool_t includemodel[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
TString strModelDir[nmodels]={"Perugia0","Perugia2010","Perugia2011","PYTHIA8","HERWIG","POWHEG","POWHEG_LO","EPOS3"};
TString strModelDirLeg[nmodels]={"PYTHIA6, Perugia 0","PYTHIA6, Perugia 2010","PYTHIA6, Perugia 2011","PYTHIA8, Tune 4C","HERWIG 7","POWHEG+PYTHIA6","POWHEG LO+PYTHIA6","EPOS3"};
//Color_t modelColors[nmodels]={kRed+2,kCyan,kGreen+2,kMagenta+1,kOrange+1,kBlue,kViolet,kYellow+1};
Color_t modelColors[nmodels]={kCyan,kYellow+1,kGreen+2,kViolet,kOrange+1,kBlue,kRed,kMagenta+1};
Bool_t includeinlegend[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kFALSE};// this is also used to split the legend in 2!!
Int_t modelMarkerStyle[nmodels]={4,33,kFullSquare,kOpenDiamond,kFullDiamond,kOpenSquare,kOpenCircle,3};
Int_t modelMarkerStyleRatio[nmodels]={4,33,kFullSquare,kOpenDiamond,kFullDiamond,kOpenSquare,kOpenCircle,3};
TString strRefForRatios="POWHEG"; //**model for which doing the division of data and other theaory curves**

TH1D **hMC, **hMC1, **hMC2, **hMC3;
TGraphAsymmErrors **grMC, **grMC1, **grMC2, **grMC3;

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
  includemodel[3]=incl;
  includeinlegend[3]=incl;
}
void IncludePowheg(Bool_t incl){
  includemodel[4]=incl;
  includeinlegend[4]=incl;
}
void IncludeEPOS(Bool_t incl){
  includemodel[5]=incl;
  includeinlegend[5]=incl;
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

void SetDirectoryFitResultsMCPP(TString str){
  strFitResultMC[0]=str;
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
  const Int_t nbinsxx=4;
  Double_t x[nbinsxx], y[nbinsxx], ex1[nbinsxx], ex2[nbinsxx], ey1[nbinsxx], ey2[nbinsxx];

  for(int i=0 ; i<nbinsxx; i++) {
    x[i] = h->GetBinCenter(i+3)+shift;
    y[i] = h->GetBinContent(i+3);
    ex1[i] = h->GetBinCenter(i+3)-h->GetBinLowEdge(i+3)+shift;
    ex2[i] = ex1[i]-2*shift;
    ey1[i] = h->GetBinError(i+3);
    ey2[i] = h->GetBinError(i+3);
  }

  gr = new TGraphAsymmErrors(nbinsxx,x,y,ex1,ex2,ey1,ey2);
  return;

}
/*void ConvertTH1ToTGraphAsymmError2016(TH1D* h,TGraphAsymmErrors *&gr, Double_t shift) {
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
*/
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

void Set4x2PadPositions(TCanvas* c){
    
  TPad * pd1 =(TPad*) c->GetPad(1);
  TPad * pd2 =(TPad*) c->GetPad(2);
  TPad * pd3 =(TPad*) c->GetPad(3);
  TPad * pd4 =(TPad*) c->GetPad(4);
  TPad * pd5 =(TPad*) c->GetPad(5);
  TPad * pd6 =(TPad*) c->GetPad(6);
  TPad * pd7 =(TPad*) c->GetPad(7);
  TPad * pd8 =(TPad*) c->GetPad(8);

  SetPadStyle(pd1);
  SetPadStyle(pd2);
  SetPadStyle(pd3);
  SetPadStyle(pd4);
  SetPadStyle(pd5);
  SetPadStyle(pd6);
  SetPadStyle(pd7);
  SetPadStyle(pd8);

  Double_t xl,xu,yl,yu;
  Double_t marginLeft=0.06;
  Double_t marginRight=0.02;
  Double_t marginTop=0.04;
  Double_t marginBottom=0.1;
  innerPadWidth=(1-marginLeft-marginRight)/4.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/2.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!
  Printf("innerPadHeight: %f",innerPadHeight);
  Printf("innerPadWidth: %f",innerPadWidth);
  Double_t marginLeftForXAxis=0.2;
  Double_t marginBottomForYAxis=0.2;

 // Bottom row

    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 7 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd5->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd5->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd5->SetRightMargin(0);
    pd5->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd5->SetTopMargin(0.);
    pd5->Modified();
    pd5->SetFillStyle(0);

    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 10 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd6->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    //pd8->SetLeftMargin(0.004/(1.-innerPadWidth-marginLeft+0.004));
    pd6->SetRightMargin(0.);
    pd6->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd6->SetTopMargin(0);
    pd6->SetFillStyle(0);

    pd6->Modified();
    pd6->Update();

    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 11 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
  
    pd7->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,3*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd7->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd7->SetRightMargin(0);
    pd7->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd7->SetTopMargin(0.);
    pd7->SetFillStyle(0);

    pd7->Modified();
    pd7->Update();

    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 12 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd8->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd8->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd8->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd8->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd8->SetTopMargin(0.);
    pd8->SetFillStyle(0);

    pd8->Modified();
    pd8->Update();

    // Top Row
    pd1->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd1->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd1->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd1->SetRightMargin(0.);
    pd1->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd1->SetTopMargin(0.);
    pd1->SetFillStyle(0);

    pd1->Modified();
    pd1->Update();

    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,2*innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd2->SetRightMargin(0.);
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd2->SetTopMargin(0.);
    pd2->SetFillStyle(0);

    pd2->Modified();
    pd2->Update();

    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,3*innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd3->SetRightMargin(0.);
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd3->SetTopMargin(0.);

    pd3->Modified();
    pd3->Update();

    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1.,2.*innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd4->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd4->SetTopMargin(0.);
    
    pd4->Modified();
    pd4->Update();

    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();

}

TCanvas* CompareNSyieldPPtoPPb(Int_t binD){
  return ComparePPtoPPb(binD,0);
}

TCanvas* CompareNSsigmaPPtoPPb(Int_t binD){
  return  ComparePPtoPPb(binD,1);

}

TCanvas* ComparePedestalPPtoPPb(Int_t binD){
  return ComparePPtoPPb(binD,2);
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
    if(collSystem==0)tlrap=new TLatex(0.016/gPad->GetWNDC()+gPad->GetLeftMargin(),0.29/gPad->GetHNDC()+gPad->GetBottomMargin(),str.Data()); 
    else if(collSystem==1)tlrap=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),str.Data()); 
    tlrap->SetNDC();
    tlrap->SetTextFont(43);
    tlrap->SetTextSize(24*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
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
  
  tlDEta=new TLatex(0.016/gPad->GetWNDC()+gPad->GetLeftMargin(),0.21/gPad->GetHNDC()+gPad->GetBottomMargin(),"|#Delta#eta| < 1"); 

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
    if(collSystem==0)tlsystem=new TLatex(0.24,0.8,"pp, #sqrt{#it{s}} = 5.02 TeV");
    else if(collSystem==1)tlsystem=new TLatex(0.24,0.8,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    else return;
    tlsystem->SetNDC();
    tlsystem->SetTextFont(42);
    tlsystem->SetTextSize(0.03);
    tlsystem->SetTextAlign(12);
  }
  else{
    if(collSystem==0)tlsystem=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"pp, #sqrt{#it{s}} = 5.02 TeV"); 
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
Double_t x=0.018,y=0.390;
  if(ncolumns==3){
    x=0.018;
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
    alice=new TLatex(0.75,0.85,"ALICE Preliminary");
    alice->SetNDC();
    alice->SetTextFont(42);
    alice->SetTextSize(0.03);
    alice->SetTextAlign(11);
  }
  else{
    alice= new TLatex(0.016/gPad->GetWNDC()+gPad->GetLeftMargin(),0.38/gPad->GetHNDC()+gPad->GetBottomMargin(),"ALICE Preliminary"); 
    alice->SetNDC();
    alice->SetTextFont(43);
    alice->SetTextSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);// draft 2 was: 28 *...
    alice->SetTextAlign(11);
  }


  //  TPaveText *alice = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  alice->SetName(Form("tlALICE_%d",identifier));
  return alice;
}

TLatex* GetDPtText(Int_t binD,Int_t identifier,Int_t addDEta=1){
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
    tlasspt=new TLatex(0.25,0.78,Form("%s, |#Delta#eta| < 1",strPtSCanvas[binD].Data()));
    tlasspt->SetNDC();
    tlasspt->SetTextFont(42);
    tlasspt->SetTextSize(0.03);
  }
  else{
    TString strTot=strPtDCanvas[binD];
    if(addDEta==1)strTot.Append("|#Delta#eta| < 1");
    tlasspt= new TLatex(0.115/gPad->GetWNDC()+gPad->GetLeftMargin(),0.345/gPad->GetHNDC()+gPad->GetBottomMargin(),strTot.Data()); 
    tlasspt->SetNDC();
    tlasspt->SetTextFont(43);
    tlasspt->SetTextSize(24*innerPadHeight/referencePadHeight*resizeTextFactor);//  if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor) but see notes on top
    Printf("Height pad: %f, scal Height = %f",gPad->GetHNDC(),scaleHeightPads);
  }
  //  TPaveText *tlasspt = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  tlasspt->SetTextAlign(22);
  tlasspt->SetName(Form("tlAssocPt_%d",identifier));
  return tlasspt;

}

TLatex* GetTextSide(Int_t variable,Int_t identifier){
  TLatex *tlSide=new TLatex();    
  Double_t x=0.18,y=0.390;
  if(ncolumns==3){
    x=0.18;
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
        tlSide=new TLatex(0.155/gPad->GetWNDC()+gPad->GetLeftMargin(),0.38/gPad->GetHNDC()+gPad->GetBottomMargin(),"Near side");
        tlSide->SetNDC();
        tlSide->SetTextAlign(11);
        tlSide->SetTextFont(43);
        tlSide->SetTextSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);// draft 2 was 28*...
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
      tlSide=new TLatex(0.155/gPad->GetWNDC()+gPad->GetLeftMargin(),0.38/gPad->GetHNDC()+gPad->GetBottomMargin(),"Away side");
      tlSide->SetNDC();
      tlSide->SetTextFont(43);
      tlSide->SetTextAlign(11);
      tlSide->SetTextSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);
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
    yh=0.32;
    yl=0.22+(3-neffmod)*0.10/3.;// for neffmod>3 an optimization might be needed
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
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  if(!(strlegendHeader.IsNull()))legend->SetHeader(strlegendHeader.Data());
  
  if(hpp)legend->AddEntry(hpp,"pp, #sqrt{#it{s}} = 5.02 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
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
    legend = new TLegend(0.2,0.2,0.4,0.4);// draft 2 (2 lines only, rapidity on the same line also for p-Pb): 0.002/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.30/gPad->GetHNDC()+gPad->GetBottomMargin()
    legend->SetTextFont(43);
    legend->SetTextAlign(12);
    legend->SetTextSize(20*innerPadHeight/referencePadHeight*resizeTextFactor);// draft 2 was 20*... //if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor); but see the notes on top
  }
 /* legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);*/
  
  legend->AddEntry(hpp,"pp, #sqrt{#it{s}} = 5.02 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
  legend->AddEntry((TObject*)0,"EPJC 77 (2017) 245","");
  
/*  if(hpPb){// draft 2 was only:  legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04","lep");
    legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV,-0.96 < #it{y}^{D}_{cms} < 0.04","lep");
  }*/
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
    legend = new TLegend(0.005/gPad->GetWNDC()+gPad->GetLeftMargin(),0.235/gPad->GetHNDC()+gPad->GetBottomMargin(),0.18/gPad->GetWNDC()+gPad->GetLeftMargin(),0.312/gPad->GetHNDC()+gPad->GetBottomMargin());// draft 2 (2 lines only, rapidity on the same line also for p-Pb): 0.002/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.30/gPad->GetHNDC()+gPad->GetBottomMargin()
    legend->SetTextFont(43);
    legend->SetTextAlign(12);
    legend->SetTextSize(20*innerPadHeight/referencePadHeight*resizeTextFactor);// draft 2 was 20*... //if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor); but see the notes on top
  }
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  
  if(hpp) {
    legend->AddEntry(hpp,"","lep");
    legend->AddEntry((TObject*)0,"","");
  }
  if(hpPb){// draft 2 was only:  legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04","lep");
    legend->AddEntry(hpPb,"","lep");
  }
  legend->SetName(Form("LegendDataAndMCPPandpPb_%d",identifier));
  return legend;
}


void InitMCobjects(){
  hMC=new TH1D*[nmodels];
  hMC1=new TH1D*[nmodels];
  hMC2=new TH1D*[nmodels];
  hMC3=new TH1D*[nmodels];
  grMC=new TGraphAsymmErrors*[nmodels];
  grMC1=new TGraphAsymmErrors*[nmodels];
  grMC2=new TGraphAsymmErrors*[nmodels];
  grMC3=new TGraphAsymmErrors*[nmodels];  
}


TCanvas* CompareDatatoModels(Int_t collsystem,Int_t binD,Int_t quantity,TPad *pd=0x0,Int_t textlegendOptions=0,Int_t drawMCasLines=0,TString legendHeader=""){
  Int_t system=collsystem;

  TH1D **hData=new TH1D*[2];
  TGraphAsymmErrors **grData=new  TGraphAsymmErrors*[2];
  TGraphAsymmErrors **grDatav2=new  TGraphAsymmErrors*[1];

  TCanvas *cout=0x0;

  if(system==-1){
    cout=CreateCanvasWithDefaultStyle(Form("%sComparisonMCtoDataBothSystemsbinD%d",strquantityFile[quantity].Data(),binD));
  }
  else{
    if(!pd){
      cout=CreateCanvasWithDefaultStyle(Form("%sComparisonMCto%sDatabinD%d",strquantityFile[quantity].Data(),strSystem[system].Data(),binD));
      pd=(TPad*)cout->cd();
    }
    else{
      pd->cd();
    }
  }
  Int_t neffmod=CountNmodels();//InLegend();

  if(textlegendOptions!=999) {
    if(style==-1){
      TLegend *legend = new TLegend(0.61,0.54,0.95,0.54+neffmod*0.05);
      //    legend= new TLegend(0.61,0.54,0.95,0.54+neffmod*0.08);
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.03);
    }
    else{
      if(textlegendOptions%1000>=100){
        TLegend *legend=GetLegendMCDataPoints(0x0,0x0,10*quantity+binD,legendHeader);
        legend->Draw();   
      }
    }
  }

  //new TCanvas(Form("NSyieldComparisonbinD%d",binD),Form("NSyieldComparisonbinD%d",binD),800,800);
  TLatex *tlCollSystem=0x0;
  TLatex *tlDrap=0x0;
  
  TFile *f1,*f2,*f3, *f1mc, *f2mc, *f3mc;
  if(system==0){
    f1=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[0].Data()),"READ");
    f2=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[1].Data()),"READ");
    f3=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[2].Data()),"READ");
  } else if(system==1){
    f1=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[0].Data()),"READ");
    f2=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[1].Data()),"READ");
    f3=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[2].Data()),"READ");
  }
  TCanvas *c1=(TCanvas*)f1->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TCanvas *c2=(TCanvas*)f2->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TCanvas *c3=(TCanvas*)f3->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TH1D *hInput1=(TH1D*)c1->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  TH1D *hInput2=(TH1D*)c2->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  TH1D *hInput3=(TH1D*)c3->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));

  Double_t xaxis[6] = {0.,0.3,1.,2.,3.,4.};
  hData[0]=new TH1D(Form("%sPPb",hInput1->GetName()),"Hist_obs_vsPtAss",5,xaxis); //the real one!
  if(system==1) {
    hData[0]->SetBinContent(1,-1);
    hData[0]->SetBinContent(2,hInput1->GetBinContent(binD+2));
    hData[0]->SetBinContent(3,hInput2->GetBinContent(binD+2));
    hData[0]->SetBinContent(4,hInput3->GetBinContent(binD+2));
    hData[0]->SetBinContent(5,-1);
    hData[0]->SetBinError(1,0);
    hData[0]->SetBinError(2,hInput1->GetBinError(binD+2));
    hData[0]->SetBinError(3,hInput2->GetBinError(binD+2));
    hData[0]->SetBinError(4,hInput3->GetBinError(binD+2));
    hData[0]->SetBinError(5,0);
    //printf("BinD = %d\n, pPb in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+2));
  } else {  //pp has also additional bin 2-3 to be skipped, so bin 3-5 is #3 (0-2, 2-3, 3-5)
    hData[0]->SetBinContent(1,-1);
    hData[0]->SetBinContent(2,hInput1->GetBinContent(binD+3));
    hData[0]->SetBinContent(3,hInput2->GetBinContent(binD+3));
    hData[0]->SetBinContent(4,hInput3->GetBinContent(binD+3));
    hData[0]->SetBinContent(5,-1);
    hData[0]->SetBinError(1,0);
    hData[0]->SetBinError(2,hInput1->GetBinError(binD+3));
    hData[0]->SetBinError(3,hInput2->GetBinError(binD+3));
    hData[0]->SetBinError(4,hInput3->GetBinError(binD+3));
    hData[0]->SetBinError(5,0);
    //printf("BinD = %d\n, pp in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+3));
  }

  TGraphAsymmErrors *gr1=(TGraphAsymmErrors*)c1->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr2=(TGraphAsymmErrors*)c2->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr3=(TGraphAsymmErrors*)c3->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  
  Double_t x[3] = {0.65,1.5,2.5};
  Double_t ex[3] = {0.35,0.5,0.5};
  Double_t y[3];
  Double_t eyl[3];  
  Double_t eyh[3];  
  Double_t *y1 = gr1->GetY();
  Double_t *y2 = gr2->GetY();
  Double_t *y3 = gr3->GetY();
  if(system==0) {
    y[0] = y1[binD+1]; y[1] = y2[binD+1]; y[2] = y3[binD+1];
    eyl[0] = gr1->GetErrorYlow(binD+1); eyl[1] = gr2->GetErrorYlow(binD+1); eyl[2] = gr3->GetErrorYlow(binD+1);
    eyh[0] = gr1->GetErrorYhigh(binD+1); eyh[1] = gr2->GetErrorYhigh(binD+1); eyh[2] = gr3->GetErrorYhigh(binD+1);
  }    
  if(system==1) {
    y[0] = y1[binD]; y[1] = y2[binD]; y[2] = y3[binD];
    eyl[0] = gr1->GetErrorYlow(binD); eyl[1] = gr2->GetErrorYlow(binD); eyl[2] = gr3->GetErrorYlow(binD);
    eyh[0] = gr1->GetErrorYhigh(binD); eyh[1] = gr2->GetErrorYhigh(binD); eyh[2] = gr3->GetErrorYhigh(binD);
  }

  grData[0] = new TGraphAsymmErrors(3,x,y,ex,ex,eyl,eyh); //the real one!
  grData[0]->SetFillColor(kGreen+2);
  grData[0]->SetFillStyle(0);
  grData[0]->SetName(Form("%sPPb",gr1->GetName()));
  grData[0]->SetMarkerColor(colSystem[system]);
  grData[0]->SetLineColor(colSystem[system]);
  grData[0]->SetLineWidth(1);
  grData[0]->SetMarkerStyle(markerStyle[system]);
  grData[0]->SetMarkerSize(markersize);

  if(textlegendOptions==10){ //Fabiom, only dta legend
    TLegend *legend= new TLegend(0.24,0.67,0.94,0.75);
    legend->SetLineWidth(0);
    legend->SetLineColor(kWhite);
    legend->SetTextFont(43);
    legend->SetTextAlign(12);
    legend->SetTextSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);
    legend->AddEntry(hData[0],"pp, #sqrt{#it{s}} = 5.02 TeV","lep");
    legend->Draw();   
  }

  if(system==1 && plotv2unc==kTRUE) {
    TGraphAsymmErrors *gr1v2=(TGraphAsymmErrors*)c1->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));
    TGraphAsymmErrors *gr2v2=(TGraphAsymmErrors*)c2->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));
    TGraphAsymmErrors *gr3v2=(TGraphAsymmErrors*)c3->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));    
    Double_t yv2[3];
    Double_t eylv2[3];  
    Double_t eyhv2[3];  
    Double_t *y1v2 = gr1v2->GetY();
    Double_t *y2v2 = gr2v2->GetY();
    Double_t *y3v2 = gr3v2->GetY();
    yv2[0] = y1v2[binD]; yv2[1] = y2v2[binD]; yv22] = y3v2[binD];
    eylv2[0] = gr1v2->GetErrorYlow(binD); eylv2[1] = gr2v2->GetErrorYlow(binD); eylv2[2] = gr3v2->GetErrorYlow(binD);
    eyhv2[0] = gr1v2->GetErrorYhigh(binD); eyhv2[1] = gr2v2->GetErrorYhigh(binD); eyhv2[2] = gr3v2->GetErrorYhigh(binD);

    grDatav2[0]= new TGraphAsymmErrors(3,x,y,ex,ex,eyl,eyh); //the real one!
    grDatav2[0]->SetName(Form("%sv2",gr1v2->GetName()));
    grDatav2[0]->SetFillColor(kGreen+2);
    grDatav2[0]->SetFillStyle(3005);
  }

  for(Int_t iPoint=0;iPoint<3;iPoint++) {
    grData[0]->SetPointError(iPoint,0.3*grData[0]->GetErrorXlow(iPoint),0.3*grData[0]->GetErrorXhigh(iPoint),grData[0]->GetErrorYlow(iPoint),grData[0]->GetErrorYhigh(iPoint));
    if(plotv2unc==kTRUE) grDatav2[0]->SetPointError(iPoint,0.65*grDatav2[0]->GetErrorXlow(iPoint),0.65*grDatav2[0]->GetErrorXhigh(iPoint),grDatav2[0]->GetErrorYlow(iPoint),grDatav2[0]->GetErrorYhigh(iPoint));
  }

 // AdaptRangeHist(hData[0],minptData,maxptData);
 // AdaptRangeTGraph(grData[0],minptData,maxptData);
 // if(plotv2unc==kTRUE) AdaptRangeTGraph(grDatav2[0],minptData,maxptData);
    
  printf("binD= %d, quantity = %d\n",binD,quantity);
/*  if((binD==1 && quantity==3) || (binD==1 && quantity==4)) { //remove NSy and NSw of 0.3-1 in 3-5 and 16-24
     grData[0]->RemovePoint(3);
     grData[0]->RemovePoint(0);
     hData[0]->SetBinContent(3,0);
     hData[0]->SetBinError(3,0);
     hData[0]->SetBinContent(6,0);
     hData[0]->SetBinError(6,0);
     printf("Removing NSy and NSw of 0.3-1 in 3-5 and 16-24\n");
  }*/
  if(quantity==3 || quantity==4) { //remove NSy and NSw of 0.3-1 in 3-5 and 16-24
     if(binD==0) {hData[0]->SetBinContent(2,-1); hData[0]->SetBinError(2,0);}
     if(binD==3) {hData[0]->SetBinContent(2,-1); hData[0]->SetBinError(2,0);}
     if(binD==0) grData[0]->RemovePoint(0);
     if(binD==3) grData[0]->RemovePoint(0);
     printf("Removing NSy and NSw of 0.3-1 in 3-5 and 16-24\n");
   }  

  pd->cd();
  TH2D *hDraw;
  if(style==-1){
    /*hData[0]->SetXTitle("Assoc. track #it{p}_{T} (GeV/#it{c})");
    hData[0]->SetYTitle(yaxisTitle[quantity].Data());
    hData[0]->GetYaxis()->SetTitleSize(0.04);
    hData[0]->GetYaxis()->SetTitleOffset(1.2);
    hData[0]->GetYaxis()->SetLabelSize(0.04);
    hData[0]->GetXaxis()->SetTitleSize(0.04);
    hData[0]->GetXaxis()->SetLabelSize(0.04);
    hData[0]->Draw();
    //  hData[0]->Draw("E0X0");// to avoid plotting the error along x
    if(system==0)hData[0]->GetYaxis()->SetRangeUser(0,maxRangePP[2][quantity]);
    if(system==1)hData[0]->GetYaxis()->SetRangeUser(0,maxRangePPb[binD][quantity]);
    
    hData[0]->SetLineColor(colSystem[system]);
    hData[0]->SetLineWidth(1);
    hData[0]->SetMarkerColor(colSystem[system]);
    hData[0]->SetMarkerStyle(markerStyle[system]);
    hData[0]->SetMarkerSize(markersize);
    */
  }
  else{
    hDraw=new TH2D(Form("hDraw%d",10*quantity+binD),"",100,0,28,200,0,10);
    hDraw->GetYaxis()->SetTitle("");      

    hDraw->GetXaxis()->SetRangeUser(0,3.2);
    if(system==0){
      hData[0]->GetYaxis()->SetRangeUser(0,maxRangePP[2][quantity]);  //THIS IS THE LINE WHICH IS CALLED!
      hDraw->GetYaxis()->SetRangeUser(0,maxRangePP[binD][quantity]);
    }
    if(system==1){
      hData[0]->GetYaxis()->SetRangeUser(0,maxRangePPb[binD][quantity]);
      hDraw->GetYaxis()->SetRangeUser(0,maxRangePPb[binD][quantity]);
      if(textlegendOptions==999) hData[0]->GetYaxis()->SetRangeUser(0,maxRangePPbSinglePanel[binD][quantity]);
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
    hDraw->GetYaxis()->SetLabelOffset(0.012);

    if(textlegendOptions%10==2||textlegendOptions%10==3||textlegendOptions==10){
      hDraw->SetYTitle(yaxisTitle[quantity].Data());
    }
    else {
      hDraw->GetYaxis()->SetLabelSize(0);
    }
    if(textlegendOptions%10==1||textlegendOptions%10==3){
      hDraw->SetXTitle("Assoc. track #it{p}_{T} (GeV/#it{c})");
    }
    else {
      hDraw->GetXaxis()->SetLabelSize(0);
    }
   //the following is always done
    hData[0]->SetXTitle("Assoc. track #it{p}_{T} (GeV/#it{c})");
    hData[0]->SetYTitle(yaxisTitle[quantity].Data());
    hData[0]->GetYaxis()->SetTitleSize(0.04);
    hData[0]->GetYaxis()->SetTitleOffset(1.2);
    hData[0]->GetYaxis()->SetLabelSize(0.04);
    hData[0]->GetYaxis()->SetLabelOffset(0.012);
    hData[0]->GetXaxis()->SetTitleSize(0.04);
    hData[0]->GetXaxis()->SetLabelSize(0.04);

    if(textlegendOptions==999) {
      TPaveText *pv = new TPaveText(0.6,0.82,0.88,0.91,"NDC");
      pv->AddText(Form("%s",strPtAssCanvas[binD].Data()));
      pv->SetBorderSize(0);
      pv->SetFillColor(0);
      pv->SetTextFont(42);
      pv->SetTextSize(0.033);
      hData[0]->SetTitle("");
      hData[0]->Draw();
      pv->Draw();
    }
    else {
      hDraw->Draw(); //FABIO... WHAT WAS THIS???
      hData[0]->Draw("same");      
    }

  }
  hData[0]->SetMarkerSize(markersize);
  hData[0]->SetMarkerStyle(markerStyle[system]);
  hData[0]->SetLineColor(colSystem[system]);
  hData[0]->SetMarkerColor(colSystem[system]);
  grData[0]->SetMarkerColor(colSystem[system]);
  grData[0]->SetLineColor(colSystem[system]);
  grData[0]->SetLineWidth(1);
  grData[0]->SetMarkerStyle(markerStyle[system]);
  grData[0]->SetMarkerSize(markersize);
  grData[0]->Draw("E2");
  if(plotv2unc==kTRUE) {
    grDatav2[0]->SetMarkerColor(kGreen-2);
    grDatav2[0]->SetLineColor(kGreen-2);
    grDatav2[0]->SetFillStyle(3001);
    grDatav2[0]->Draw("E2");
  }

  TLegend *legend;
  if(textlegendOptions==999) {
    if(style==-1){
      legend = new TLegend(0.61,0.44,0.95,0.54+neffmod*0.05);
      //    legend= new TLegend(0.61,0.54,0.95,0.54+neffmod*0.08);
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.03);
    }
    else{
      if(textlegendOptions%1000>=100 && textlegendOptions!=10){
        legend=GetLegendMCDataPoints(0x0,0x0,10*quantity+binD,legendHeader);
        legend->AddEntry(hData[0],"pp, #sqrt{#it{s}} = 5.02 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
        legend->AddEntry(grData[0],"Total data syst unc","f");
        legend->AddEntry(grDatav2[0],"syst unc from v_{2}","f");
        legend->SetX1(0.2);
        legend->SetX2(0.55);
        legend->SetY1(0.65);
        legend->SetY2(0.9); 
        legend->Draw();   
      }
    }
  }

    if(collsystem==-1){
      legend->AddEntry(hData[0],"pp, #sqrt{#it{s}} = 5.02 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
      legend->AddEntry(hData[1],"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, -0.96 < #it{y}^{D}_{cms} < 0.04","lep");      
    }
    else {    
      //      legend->AddEntry(hData[0],"data","lep");
      if(textlegendOptions%1000000>=100000)tlCollSystem=GetCollSystem(system,10*quantity+binD);
      if(textlegendOptions%10000000>=1000000)tlDrap=GetDRapForSystem(system,10*quantity+binD);
    }
    
    
    // INITIATE MC HISTOS AND GRAPHS
    InitMCobjects();
    // NOW LOOP OVER MODELS
    for(Int_t kmc=0;kmc<nmodels;kmc++){
      if(!includemodel[kmc])continue;
      if(quantity==2&&strModelDir[kmc].Contains("POWHEG"))continue;
 
      f1mc=TFile::Open(Form("%s/Trends_%s/%s/CanvasFinalTrend%s_pthad%s.root",strFitResultMC[system].Data(),strSystem[system].Data(),strModelDir[kmc].Data(),strquantityFile[quantity].Data(),strPtAss[0].Data()),"READ");
      f2mc=TFile::Open(Form("%s/Trends_%s/%s/CanvasFinalTrend%s_pthad%s.root",strFitResultMC[system].Data(),strSystem[system].Data(),strModelDir[kmc].Data(),strquantityFile[quantity].Data(),strPtAss[1].Data()),"READ");
      f3mc=TFile::Open(Form("%s/Trends_%s/%s/CanvasFinalTrend%s_pthad%s.root",strFitResultMC[system].Data(),strSystem[system].Data(),strModelDir[kmc].Data(),strquantityFile[quantity].Data(),strPtAss[2].Data()),"READ");
      TCanvas* c1mc=(TCanvas*)f1mc->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
      TCanvas* c2mc=(TCanvas*)f2mc->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
      TCanvas* c3mc=(TCanvas*)f3mc->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
      grMC1[kmc]=(TGraphAsymmErrors*)c1mc->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
      grMC2[kmc]=(TGraphAsymmErrors*)c2mc->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
      grMC3[kmc]=(TGraphAsymmErrors*)c3mc->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
      Double_t xMC[3] = {0.65,1.5,2.5};
      Double_t exMC[3] = {0.35,0.5,0.5};
      Double_t yMC[3];
      Double_t eylMC[3];  
      Double_t eyhMC[3];  
      Double_t *y1MC = grMC1[kmc]->GetY();
      Double_t *y2MC = grMC2[kmc]->GetY();
      Double_t *y3MC = grMC3[kmc]->GetY();
      if(system==0) {
        yMC[0] = y1MC[binD+1]; yMC[1] = y2MC[binD+1]; yMC[2] = y3MC[binD+1];
        eylMC[0] = grMC1[kmc]->GetErrorYlow(binD+1); eylMC[1] = grMC2[kmc]->GetErrorYlow(binD+1); eylMC[2] = grMC3[kmc]->GetErrorYlow(binD+1);
        eyhMC[0] = grMC1[kmc]->GetErrorYhigh(binD+1); eyhMC[1] = grMC2[kmc]->GetErrorYhigh(binD+1); eyhMC[2] = grMC3[kmc]->GetErrorYhigh(binD+1);
      }    
      if(system==1) {
        yMC[0] = y1[binD]; yMC[1] = y2[binD]; yMC[2] = y3[binD];
        eylMC[0] = grMC1[kmc]->GetErrorYlow(binD); eylMC[1] = grMC2[kmc]->GetErrorYlow(binD); eylMC[2] = grMC3[kmc]>GetErrorYlow(binD);
        eyhMC[0] = grMC1[kmc]->GetErrorYhigh(binD); eyhMC[1] = grMC2[kmc]->GetErrorYhigh(binD); eyhMC[2] = grMC3[kmc]->GetErrorYhigh(binD);
      }

      grMC[kmc] = new TGraphAsymmErrors(3,xMC,yMC,exMC,exMC,eylMC,eyhMC); //the real one!
      grMC[kmc]->SetFillStyle(0);
      grMC[kmc]->SetName(Form("%s%s",grMC1[kmc]->GetName(),strModelDir[kmc].Data()));
      grMC[kmc]->SetMarkerColor(colSystem[system]);
      grMC[kmc]->SetLineColor(colSystem[system]);
      grMC[kmc]->SetLineWidth(1);
      grMC[kmc]->SetMarkerStyle(markerStyle[system]);
      grMC[kmc]->SetMarkerSize(markersize);

      hMC1[kmc]=(TH1D*)c1mc->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
      hMC2[kmc]=(TH1D*)c2mc->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
      hMC3[kmc]=(TH1D*)c3mc->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
      
      hMC[kmc]=new TH1D(Form("%sPPb",hMC1[kmc]->GetName()),"Hist_obs_vsPtAss",5,xaxis); //the real one!
      if(system==1) {
        hMC[kmc]->SetBinContent(1,-1);
        hMC[kmc]->SetBinContent(2,hMC1[kmc]->GetBinContent(binD+2));
        hMC[kmc]->SetBinContent(3,hMC2[kmc]->GetBinContent(binD+2));
        hMC[kmc]->SetBinContent(4,hMC3[kmc]->GetBinContent(binD+2));
        hMC[kmc]->SetBinContent(5,-1);
        hMC[kmc]->SetBinError(1,0);
        hMC[kmc]->SetBinError(2,hMC1[kmc]->GetBinError(binD+2));
        hMC[kmc]->SetBinError(3,hMC2[kmc]->GetBinError(binD+2));
        hMC[kmc]->SetBinError(4,hMC3[kmc]->GetBinError(binD+2));
        hMC[kmc]->SetBinError(5,0);
        //printf("BinD = %d\n, pPb in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+2));
      } else {  //pp has also additional bin 2-3 to be skipped, so bin 3-5 is #3 (0-2, 2-3, 3-5)
        hMC[kmc]->SetBinContent(1,-1);
        hMC[kmc]->SetBinContent(2,hMC1[kmc]->GetBinContent(binD+3));
        hMC[kmc]->SetBinContent(3,hMC2[kmc]->GetBinContent(binD+3));
        hMC[kmc]->SetBinContent(4,hMC3[kmc]->GetBinContent(binD+3));
        hMC[kmc]->SetBinContent(5,-1);
        hMC[kmc]->SetBinError(1,0);
        hMC[kmc]->SetBinError(2,hMC1[kmc]->GetBinError(binD+3));
        hMC[kmc]->SetBinError(3,hMC2[kmc]->GetBinError(binD+3));
        hMC[kmc]->SetBinError(4,hMC3[kmc]->GetBinError(binD+3));
        hMC[kmc]->SetBinError(5,0);
        //printf("BinD = %d\n, pp in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+3));
      }

      // AdaptRangeHist(hMC[kmc],minptMC,maxptMC);
      // AdaptRangeTGraph(grMC[kmc],minptMC,maxptMC);

      pd->cd();
      if(drawMCasLines==2){
	hMC[kmc]->Draw("C X0 same");
      }
      else{
	hMC[kmc]->Draw("same");
      }
      hMC[kmc]->SetLineColor(modelColors[kmc]);
      hMC[kmc]->SetLineWidth(1);
      if(drawMCasLines==1){
	hMC[kmc]->SetLineStyle(2);
      }
      hMC[kmc]->SetMarkerColor(modelColors[kmc]);
      hMC[kmc]->SetMarkerStyle(modelMarkerStyle[kmc]);
      hMC[kmc]->SetMarkerSize(markersizeMC);
      grMC[kmc]->SetMarkerColor(modelColors[kmc]);
      grMC[kmc]->SetLineColor(modelColors[kmc]);
      grMC[kmc]->SetLineWidth(1);
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
      if(legend&&includeinlegend[kmc] && textlegendOptions!=10){
      legend->AddEntry(hMC[kmc],Form("%s",strModelDirLeg[kmc].Data()),"lep");  
      }
      

    } //kmc loop
    
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
      TLatex *tlAssYieldPt=new TLatex(0.24,0.7,Form("%s, |#Delta#eta| < 1",strPtAssCanvas[binD].Data()));
      tlAssYieldPt->SetNDC();
      tlAssYieldPt->SetTextSize(0.03);
      tlAssYieldPt->Draw();
    }
    else{
      
      if(textlegendOptions%100>=10){
	TLatex *tlAssYieldPt=GetDPtText(binD,10*quantity+binD,0);
	tlAssYieldPt->Draw();
      }
      if(textlegendOptions%10000>=1000){
	TLatex *tlALICE=GetALICEtext(10*quantity+binD);
	tlALICE->Draw();
      }
      if(textlegendOptions%100000>=10000){
	if(quantity!=2){
	  TLatex *tlSide=GetTextSide(quantity,10*quantity+binD);
	  tlSide->Draw();
	}
      }
    }
    

    if(plotv2unc==kTRUE) grDatav2[0]->Draw("E2");
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
	//  legData->AddEntry(hData[0],tlCollSystem->GetTitle(),"lep");
	  //gPad->GetWNDC()
	  legData->SetTextAlign(12);
	  legData->SetTextFont(43);
          legData->SetLineColor(kWhite);
	  legData->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
	  legData->Draw();
	  if(system==1) {
	  TLegend *legDataSuperimp=new TLegend(gPad->GetLeftMargin()+0.02,tlCollSystem->GetY()-0.02,0.95-gPad->GetRightMargin(),tlCollSystem->GetY()+tlCollSystem->GetYsize()-0.02,"");
	    legDataSuperimp->AddEntry(hSuperimp,"","lep");
	    legDataSuperimp->SetFillStyle(0);
	    legDataSuperimp->SetTextAlign(12);
	    legDataSuperimp->SetTextFont(43);
	    legDataSuperimp->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
	    legDataSuperimp->Draw("same");
 	  }
	}
	//else tlCollSystem->Draw();
      }
      if(tlDrap!=0x0){
	tlDrap->Draw();
      }
    }
 
    return cout;
    
}


void CompareFitResultsPPtoMCUniqueCanvas_vsPtAss() {
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
  cFinalPaperStyle=new TCanvas("cPPvsPPbFitResultsFinalPaperStyle","pp vs.pPb fit results ",1800./1200.*canvasheight,canvasheight);
  //SetPadStyle(cFinalPaperStyle);
  //    cFinalPaperStyle->SetBottomMargin(0);
  //    cFinalPaperStyle->SetTopMargin(0);
  cFinalPaperStyle->Divide(4,2,0.0,0.0,0);
  cFinalPaperStyle->SetTicky();
  cFinalPaperStyle->SetTickx();
  cFinalPaperStyle->SetFrameBorderMode(0);

  Set4x2PadPositions(cFinalPaperStyle);
  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  //    Convert3x2Matrix(cFinalPaperStyle,kTRUE); 

  Bool_t includeinlegendOrig[nmodels];
  for(Int_t k=0;k<nmodels;k++)includeinlegendOrig[k]=includeinlegend[k];

  Int_t orderD[4]={0,1,2,3};  
  for(Int_t jp=0;jp<=3;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=10;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binD; NS --> 0 ; binD == jp (not orderD[jp])
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(0,orderD[jp],0,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 5.02 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];
    }
    else if(jp==0){// identifier set as: 10*quantity+binD; NS:
      CompareDatatoModels(0,orderD[jp],0,pd,needTitle);//+1000000);    
      TLatex *tlALICE=GetALICEtext(orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderD[jp]);
      tlSide->Draw();
      TLatex *tlDeta=GetDRapForSystem(0,orderD[jp],1);
      tlDeta->Draw();      
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      includeinlegend[3]=kFALSE;
      CompareDatatoModels(0,orderD[jp],0,pd,needTitle+100,0,"");// " "); the latter " " needed for counting lines properly    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];
      includeinlegend[3]=includeinlegendOrig[3];
    }
    else{
      CompareDatatoModels(0,orderD[jp],0,pd,needTitle);    
    }
    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],orderD[jp],0);
    tlAssYieldPt->Draw();

    gStyle->SetOptStat(0000);   
  }
  
  for(Int_t jp=0;jp<=3;jp++){
    Int_t needTitle=1;
    if(jp==0)needTitle=3;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+5);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderD[jp],1,pd,needTitle);    
    gStyle->SetOptStat(0000);   
  }

  for(Int_t j=8;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  if(drawSystMC){
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResults_vsPtAss.pdf");
  }
  else{
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResults_vsPtAss.pdf");
  }
  return;

}


void CompareFitResultsPPtoMCUniqueCanvasAwaySide_vsPtAss(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
  cFinalPaperStyle=new TCanvas("cPPvsPPbFitResultsFinalPaperStyle","pp vs.pPb fit results ",1800./1200.*canvasheight,canvasheight);
  //SetPadStyle(cFinalPaperStyle);
  //    cFinalPaperStyle->SetBottomMargin(0);
  //    cFinalPaperStyle->SetTopMargin(0);
  cFinalPaperStyle->Divide(4,2,0.0,0.0,0);
  cFinalPaperStyle->SetTicky();
  cFinalPaperStyle->SetTickx();
  cFinalPaperStyle->SetFrameBorderMode(0);

  Set4x2PadPositions(cFinalPaperStyle);
  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  //    Convert3x2Matrix(cFinalPaperStyle,kTRUE); 

  Bool_t includeinlegendOrig[nmodels];
  for(Int_t k=0;k<nmodels;k++)includeinlegendOrig[k]=includeinlegend[k];

  Int_t orderD[4]={0,1,2,3};  
  for(Int_t jp=0;jp<=3;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=10;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binD; NS --> 0 ; binD == jp (not orderD[jp])
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(0,orderD[jp],3,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 5.02 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];
    }
    else if(jp==0){// identifier set as: 10*quantity+binD; NS:
      CompareDatatoModels(0,orderD[jp],3,pd,needTitle);//+1000000);    
      TLatex *tlALICE=GetALICEtext(orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(3,orderD[jp]);
      tlSide->Draw();
      TLatex *tlDeta=GetDRapForSystem(0,orderD[jp],1);
      tlDeta->Draw();      
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      includeinlegend[3]=kFALSE;      
      CompareDatatoModels(0,orderD[jp],3,pd,needTitle+100,0,"");// " "); the latter " " needed for counting lines properly    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];
      includeinlegend[3]=includeinlegendOrig[3];
    }
    else{
      CompareDatatoModels(0,orderD[jp],3,pd,needTitle);    
    }
    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],orderD[jp],0);
    tlAssYieldPt->Draw();

    gStyle->SetOptStat(0000);   
  }
  
  for(Int_t jp=0;jp<=3;jp++){
    Int_t needTitle=1;
    if(jp==0)needTitle=3;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+5);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    CompareDatatoModels(0,orderD[jp],4,pd,needTitle);    

    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binD; Sigma --> 1 ; binD == orderD[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderD[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binD; NS:
      TLatex *tlALICE=GetALICEtext(10+orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderD[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],10+orderD[jp]);
    tlAssYieldPt->Draw();
    */

    gStyle->SetOptStat(0000);   
  }

  for(Int_t j=8;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  if(drawSystMC){
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss.pdf");
  }
  else{
  cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPtoMCnoSystFitResultsAS_vsPtAss.pdf");
  }
  return;

}

void CompareFitResults_Ratios_NS_vsPtAss() {
        
        TFile fIn("ComparePPtoMCFitResults_vsPtAss.root");
        TCanvas *cRat = (TCanvas*)fIn.Get("cPPvsPPbFitResultsFinalPaperStyle");
        cRat->Draw();

        Int_t nmodelsOn = CountNmodels();
        Int_t modRef=0;
        printf("nmodels on = %d\n",nmodelsOn);
        for(int k=0;k<nmodelsOn;k++) {
          while(includemodel[modRef]==kFALSE) modRef++; //because k index don't have to be updated is model is excluded, but name index yes!!                  
          printf("k = %d, modRef = %d, included? %d\n",k,modRef,includemodel[modRef]);
          modRef++;
          if(!strModelDir[modRef].CompareTo(strRefForRatios)) {
            modRef=k+1; //reconcile with k index!!!
            break;
          }
        }
        printf("at the end, modRef=%d\n",modRef);

        for(int i=0;i<8;i++) {

                TPad *pad = (TPad*)cRat->FindObject(Form("cPPvsPPbFitResultsFinalPaperStyle_%d",i+1));
                TList *l = (TList*)pad->GetListOfPrimitives();
                printf("i = %d\n",i);
                l->ls();

                Int_t nmodelsOn = CountNmodels();
                TH1D *hModRat[20];
                TGraphAsymmErrors *grModRat[20];

                TH1D *hData = (TH1D*)pad->GetListOfPrimitives()->At(2);
                hData->GetYaxis()->SetRangeUser(-5,5);
                TGraphAsymmErrors *grData = (TGraphAsymmErrors*)l->At(3);
                Double_t *xvalDa = grData->GetX();
                Double_t *yvalDa = grData->GetY();
                Double_t *eyvalDa = grData->GetEYlow();           

                TH1D *hMod = (TH1D*)pad->GetListOfPrimitives()->At(4+modRef*2);
                hModRat[modRef] = (TH1D*)hMod->Clone(Form("%s_Ratio",hMod->GetName()));
                TGraphAsymmErrors *grModRef = (TGraphAsymmErrors*)l->At(5+modRef*2);
                Double_t *xval = grModRef->GetX();
                Double_t *exval = grModRef->GetEXlow();
                Double_t *yvalMCref = grModRef->GetY();
                Double_t *eyvalMCref = grModRef->GetEYlow();

                for(int k=0;k<nmodelsOn;k++) {
                        if(k==modRef) continue;
                        TH1D *hMod = (TH1D*)pad->GetListOfPrimitives()->At(4+k*2);
                        hModRat[k] = (TH1D*)hMod->Clone(Form("%s_Ratio",hMod->GetName()));
                        hModRat[k]->Divide(hModRat[modRef]);
                        TGraphAsymmErrors *grMod = (TGraphAsymmErrors*)l->At(5+k*2);
                        grModRat[k] = (TGraphAsymmErrors*)grMod->Clone(Form("%s_Ratio",grMod->GetName()));
                        Double_t *yvalMC = grMod->GetY();
                        Double_t *eyvalMC = grMod->GetEYlow();
                        for(int ip=0;ip<grModRef->GetN();ip++) {
                          Double_t err = yvalMC[ip]/yvalMCref[ip]*(TMath::Sqrt((eyvalMC[ip]/yvalMC[ip])*(eyvalMC[ip]/yvalMC[ip])+(eyvalMCref[ip]/yvalMCref[ip])*(eyvalMCref[ip]/yvalMCref[ip])));
                          grModRat[k]->SetPoint(ip,xval[ip],yvalMC[ip]/yvalMCref[ip]);
                          grModRat[k]->SetPointError(ip,exval[ip],exval[ip],err,err);
                        }
                }

                TH1D* hDataRat = (TH1D*)hData->Clone(Form("%s_Ratio",hData->GetName()));
                hDataRat->Divide(hModRat[modRef]);
                grDataRat = (TGraphAsymmErrors*)grData->Clone(Form("%s_Ratio",grData->GetName()));

                for(int ip=0;ip<grModRef->GetN();ip++) {
                    Double_t err = yvalDa[ip]/yvalMCref[ip]*(TMath::Sqrt((eyvalDa[ip]/yvalDa[ip])*(eyvalDa[ip]/yvalDa[ip])+(eyvalMCref[ip]/yvalMCref[ip])*(eyvalMCref[ip]/yvalMCref[ip])));
                    grDataRat->SetPoint(ip,xval[ip],yvalDa[ip]/yvalMCref[ip]);
                    grDataRat->SetPointError(ip,exval[ip],exval[ip],err,err);
                }

                pad->cd();
                TH2D *hframe = (TH2D*)l->At(1);
                if(i==0) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[0][0],maxRangePP_Ratio[0][0]);
                if(i==1) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[1][0],maxRangePP_Ratio[1][0]);
                if(i==2) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[2][0],maxRangePP_Ratio[2][0]);
                if(i==3) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[3][0],maxRangePP_Ratio[3][0]);
                if(i==4) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[0][1],maxRangePP_Ratio[0][1]);
                if(i==5) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[1][1],maxRangePP_Ratio[1][1]);
                if(i==6) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[2][1],maxRangePP_Ratio[2][1]);
                if(i==7) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[3][1],maxRangePP_Ratio[3][1]);
                
                Int_t num = l->GetEntries()-1;
                while(num>1) {
                  TObject *obj=l->At(num);
                  TString strName=obj->ClassName();
                  if(strName.Contains("TH1D")||strName.Contains("TGraphAsymmErrors")) l->Remove(l->At(num));
                  if(strName.Contains("TLegend")) l->Remove(l->At(num));
                  num--;
                }

                for(int k=0;k<nmodelsOn;k++) {
                  if(k==modRef) continue;
                  hModRat[k]->SetBinContent(1,-99);
                  hModRat[k]->SetBinError(1,-99);
                  hModRat[k]->DrawCopy("same");
                  grModRat[k]->SetFillStyle(1);
                  grModRat[k]->SetMarkerStyle(modelMarkerStyleRatio[k]);
                  grModRat[k]->Draw("E2");
                }
                hDataRat->SetBinContent(1,-99);
                hDataRat->SetBinError(1,-99);
                hDataRat->Draw("same");
                grDataRat->Draw("E2");

                pad->Update(); 
        }

        cRat->SaveAs("ComparePPtoMCFitResults_vsPtAss_Ratio.root");
        cRat->SaveAs("ComparePPtoMCFitResults_vsPtAss_Ratio.eps");
        cRat->SaveAs("ComparePPtoMCFitResults_vsPtAss_Ratio.png");
        cRat->SaveAs("ComparePPtoMCFitResults_vsPtAss_Ratio.pdf");    
}

void CompareFitResults_Ratios_AS_vsPtAss() {
        
        TFile fIn("ComparePPtoMCFitResultsAS_vsPtAss.root");
        TCanvas *cRat = (TCanvas*)fIn.Get("cPPvsPPbFitResultsFinalPaperStyle");
        cRat->Draw();

        Int_t nmodelsOn = CountNmodels();
        Int_t modRef=0;
        printf("nmodels on = %d\n",nmodelsOn);
        for(int k=0;k<nmodelsOn;k++) {
          while(includemodel[modRef]==kFALSE) modRef++; //because k index don't have to be updated is model is excluded, but name index yes!!                  
          printf("k = %d, modRef = %d, included? %d\n",k,modRef,includemodel[modRef]);
          modRef++;
          if(!strModelDir[modRef].CompareTo(strRefForRatios)) {
            modRef=k+1; //reconcile with k index!!!
            break;
          }
        }
        printf("at the end, modRef=%d\n",modRef);

        for(int i=0;i<8;i++) {

                TPad *pad = (TPad*)cRat->FindObject(Form("cPPvsPPbFitResultsFinalPaperStyle_%d",i+1));
                TList *l = (TList*)pad->GetListOfPrimitives();
                printf("i = %d\n",i);
                l->ls();

                Int_t nmodelsOn = CountNmodels();
                TH1D *hModRat[20];
                TGraphAsymmErrors *grModRat[20];

                TH1D *hData = (TH1D*)pad->GetListOfPrimitives()->At(2);
                hData->GetYaxis()->SetRangeUser(-5,5);
                TGraphAsymmErrors *grData = (TGraphAsymmErrors*)l->At(3);
                Double_t *xvalDa = grData->GetX();
                Double_t *yvalDa = grData->GetY();
                Double_t *eyvalDa = grData->GetEYlow();           

                TH1D *hMod = (TH1D*)pad->GetListOfPrimitives()->At(4+modRef*2);
                hModRat[modRef] = (TH1D*)hMod->Clone(Form("%s_Ratio",hMod->GetName()));
                TGraphAsymmErrors *grModRef = (TGraphAsymmErrors*)l->At(5+modRef*2);
                Double_t *xval = grModRef->GetX();
                Double_t *exval = grModRef->GetEXlow();
                Double_t *yvalMCref = grModRef->GetY();
                Double_t *eyvalMCref = grModRef->GetEYlow();

                for(int k=0;k<nmodelsOn;k++) {
                        if(k==modRef) continue;
                        TH1D *hMod = (TH1D*)pad->GetListOfPrimitives()->At(4+k*2);
                        hModRat[k] = (TH1D*)hMod->Clone(Form("%s_Ratio",hMod->GetName()));
                        hModRat[k]->Divide(hModRat[modRef]);
                        TGraphAsymmErrors *grMod = (TGraphAsymmErrors*)l->At(5+k*2);
                        grModRat[k] = (TGraphAsymmErrors*)grMod->Clone(Form("%s_Ratio",grMod->GetName()));
                        Double_t *yvalMC = grMod->GetY();
                        Double_t *eyvalMC = grMod->GetEYlow();
                        for(int ip=0;ip<grModRef->GetN();ip++) {
                          Double_t err = yvalMC[ip]/yvalMCref[ip]*(TMath::Sqrt((eyvalMC[ip]/yvalMC[ip])*(eyvalMC[ip]/yvalMC[ip])+(eyvalMCref[ip]/yvalMCref[ip])*(eyvalMCref[ip]/yvalMCref[ip])));
                          grModRat[k]->SetPoint(ip,xval[ip],yvalMC[ip]/yvalMCref[ip]);
                          grModRat[k]->SetPointError(ip,exval[ip],exval[ip],err,err);
                        }
                }

                TH1D* hDataRat = (TH1D*)hData->Clone(Form("%s_Ratio",hData->GetName()));
                hDataRat->Divide(hModRat[modRef]);
                grDataRat = (TGraphAsymmErrors*)grData->Clone(Form("%s_Ratio",grData->GetName()));

                if(xvalDa[0]==0.65) {
                    for(int ip=0;ip<grModRef->GetN();ip++) {
                        Double_t err = yvalDa[ip]/yvalMCref[ip]*(TMath::Sqrt((eyvalDa[ip]/yvalDa[ip])*(eyvalDa[ip]/yvalDa[ip])+(eyvalMCref[ip]/yvalMCref[ip])*(eyvalMCref[ip]/yvalMCref[ip])));
                        grDataRat->SetPoint(ip,xval[ip],yvalDa[ip]/yvalMCref[ip]);
                        grDataRat->SetPointError(ip,exval[ip],exval[ip],err,err);
                    }                
                }
                else {
                    for(int ip=0;ip<grData->GetN();ip++) {
                        Double_t err = yvalDa[ip]/yvalMCref[ip+1]*(TMath::Sqrt((eyvalDa[ip]/yvalDa[ip])*(eyvalDa[ip]/yvalDa[ip])+(eyvalMCref[ip+1]/yvalMCref[ip+1])*(eyvalMCref[ip+1]/yvalMCref[ip+1])));
                        grDataRat->SetPoint(ip,xval[ip+1],yvalDa[ip]/yvalMCref[ip+1]);
                        grDataRat->SetPointError(ip,exval[ip+1],exval[ip+1],err,err);
                    }
                }

                pad->cd();
                TH2D *hframe = (TH2D*)l->At(1);
                if(i==0) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[0][3],maxRangePP_Ratio[0][3]);
                if(i==1) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[1][3],maxRangePP_Ratio[1][3]);
                if(i==2) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[2][3],maxRangePP_Ratio[2][3]);
                if(i==3) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[3][3],maxRangePP_Ratio[3][3]);
                if(i==4) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[0][4],maxRangePP_Ratio[0][4]);
                if(i==5) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[1][4],maxRangePP_Ratio[1][4]);
                if(i==6) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[2][4],maxRangePP_Ratio[2][4]);
                if(i==7) hframe->GetYaxis()->SetRangeUser(minRangePP_Ratio[3][4],maxRangePP_Ratio[3][4]);
                
                Int_t num = l->GetEntries()-1;
                while(num>1) {
                  TObject *obj=l->At(num);
                  TString strName=obj->ClassName();
                  if(strName.Contains("TH1D")||strName.Contains("TGraphAsymmErrors")) l->Remove(l->At(num));
                  if(strName.Contains("TLegend")) l->Remove(l->At(num));
                  num--;
                }

                for(int k=0;k<nmodelsOn;k++) {
                  if(k==modRef) continue;
                  hModRat[k]->SetBinContent(1,-99);
                  hModRat[k]->SetBinError(1,-99);
                  hModRat[k]->DrawCopy("same");
                  grModRat[k]->SetFillStyle(1);
                  grModRat[k]->SetMarkerStyle(modelMarkerStyleRatio[k]);
                  grModRat[k]->Draw("E2");
                }
                hDataRat->SetBinContent(1,-99);
                hDataRat->SetBinError(1,-99);
                hDataRat->Draw("same");
                grDataRat->Draw("E2");

                pad->Update(); 
        }

        cRat->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss_Ratio.root");
        cRat->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss_Ratio.eps");
        cRat->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss_Ratio.png");
        cRat->SaveAs("ComparePPtoMCFitResultsAS_vsPtAss_Ratio.pdf");    
}

void CompareFitResultsPPbtoMCUniqueCanvas(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
  xtitleoffset=2.5;
  strModelDirLeg[7]="EPOS 3.117, p-Pb simulation";
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPbvsMCFitResultsFinalPaperStyle","pp vs. MC fit results ",canvasheight,canvasheight);
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
  Int_t orderD[3]={0,1,2};  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binD; NS --> 0 ; binD == jp (not orderD[jp])
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(1,orderD[jp],0,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 5.02 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];

    }
    else if(jp==0){// identifier set as: 10*quantity+binD; NS:
      CompareDatatoModels(1,orderD[jp],0,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderD[jp]);
      tlSide->Draw();
      TLatex *tlDeta=GetDRapForSystem(1,orderD[jp],1);
      tlDeta->Draw();      
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      includeinlegend[3]=kFALSE;
      CompareDatatoModels(1,orderD[jp],0,pd,needTitle+100,0,"");//" ");// the latter empty space needed for counting properly lines    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];   
      includeinlegend[3]=includeinlegendOrig[3];      
    }
    else{
      CompareDatatoModels(1,orderD[jp],0,pd,needTitle);    
    }

    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],orderD[jp],0);
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
    CompareDatatoModels(1,orderD[jp],1,pd,needTitle);    

    /* NOTHING NEEDS TO BE WRITTEN IN BOTTOM ROW
       if(jp==1){// identifier set as: 10*quantity+binD; Sigma --> 1 ; binD == orderD[jp]
      TLegend *legend=GetLegendDataPoints(hPP,hPPb,10+orderD[jp]);
      legend->Draw();
    }
    if(jp==0){// identifier set as: 10*quantity+binD; NS:
      TLatex *tlALICE=GetALICEtext(10+orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,10+orderD[jp]);
      tlSide->Draw();

    }
    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],10+orderD[jp]);
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
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResults_vsPtAss.pdf");
  }
  else {
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResults_vsPtAss.pdf");
  }
  return;

}

void CompareFitResultsPPbtoMCUniqueCanvasAwaySide(){
  printf("AAA\n");
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
  xtitleoffset=2.5;
  strModelDirLeg[7]="EPOS 3.117, p-Pb simulation";
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cPPbvsMCFitResultsFinalPaperStyleAS","pp vs. MC fit results Away Side ",1550./1200.*canvasheight,canvasheight);
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
  Int_t orderD[3]={0,1,2};  
  for(Int_t jp=0;jp<=2;jp++){
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binD; NS --> 0 ; binD == jp (not orderD[jp])
      includeinlegend[4]=kFALSE;
      includeinlegend[5]=kFALSE;
      includeinlegend[6]=kFALSE;
      includeinlegend[7]=kFALSE;
      CompareDatatoModels(1,orderD[jp],3,pd,100000+needTitle+100,0,"Simulations, pp, #sqrt{#it{s}} = 5.02 TeV");    // title + 10*asspt+100*legendDataMC+1000*ALICE+10000*side+100000*collSyst+1000000*Drap
      includeinlegend[4]=includeinlegendOrig[4];
      includeinlegend[5]=includeinlegendOrig[5];
      includeinlegend[6]=includeinlegendOrig[6];
      includeinlegend[7]=includeinlegendOrig[7];
    }
    else if(jp==0){// identifier set as: 10*quantity+binD; NS:
      CompareDatatoModels(1,orderD[jp],3,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(3,orderD[jp]);
      tlSide->Draw();
    }
    else if(jp==2){
      includeinlegend[0]=kFALSE;
      includeinlegend[1]=kFALSE;
      includeinlegend[2]=kFALSE;
      includeinlegend[3]=kFALSE;
      CompareDatatoModels(1,orderD[jp],3,pd,needTitle+100,0,"");//" ");// the latter empty space needed for counting properly lines    
      includeinlegend[0]=includeinlegendOrig[0];
      includeinlegend[1]=includeinlegendOrig[1];
      includeinlegend[2]=includeinlegendOrig[2];
      includeinlegend[3]=includeinlegendOrig[3];
      TLatex *tlDeta=GetDRapForSystem(1,orderD[jp],1);
      tlDeta->Draw();
    
    }
    else{
      CompareDatatoModels(1,orderD[jp],3,pd,needTitle);    
    }

    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],orderD[jp],0);
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
    CompareDatatoModels(1,orderD[jp],4,pd,needTitle);    

    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=6;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  if(drawSystMC){
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCFitResultsAS_vsPtAss.pdf");
  }
  else {
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS_vsPtAss.root");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS_vsPtAss.eps");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS_vsPtAss.png");
    cFinalPaperStyle->SaveAs("ComparePPbtoMCnoSystFitResultsAS_vsPtAss.pdf");
  }
  return;

}
