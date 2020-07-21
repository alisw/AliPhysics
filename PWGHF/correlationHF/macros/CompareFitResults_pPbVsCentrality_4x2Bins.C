TString strPtAss[4]={"0.3to99.0","0.3to1.0","1.0to2.0","2.0to3.0"};
TString strPtAssCanvas[4]={"#it{p}_{T}^{assoc} > 0.3 GeV/#it{c}","0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}","1 < #it{p}_{T}^{assoc} < 2 GeV/#it{c}","2 < #it{p}_{T}^{assoc} < 3 GeV/#it{c}"};
TString strSystem[4]={"0-0.1%","0.1-10%","10-30%","30-100%"};
Color_t colSystem[4]={kBlue,kRed,kGreen+2,kBlack};
Int_t markerStyle[4]={20,21,22,23};
Bool_t useLegendForData=kTRUE;
Bool_t plotv2unc=kTRUE;
TString strFitResultPPb[4]={"./001","./0110","./1030","./30100"}; //  "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";
Double_t canvasheight=1000;
Double_t resizeTextFactor=1.3;//1.*canvasheight/1200.; // size was tuned for canvasheight =1800. 
Double_t referencePadHeight=0.48; // Do not touch unless the canvas division is changed from 2x3, change canvasheight and resizeTextFactor instead
Double_t innerPadHeight;// not touch, set internally
Double_t innerPadWidth;// not touch, set internally
Int_t style=1;
Double_t scaleHeightPads=1;// do not touch this, it is regulated automatically in the macro
Double_t scaleWidthPads=1;// do not touch this, it is regulated automatically in the macro
TString strFitResultMC[4]={"","","",""};
//TString strFitResultMC[2]={"/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/MCTemplates/Templates_pp_12May15/FitResults/",
//			   "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May19UseScriptPWGHF/MCTemplates/Templates_pPb_12May15/FitResults/"};

TString strquantityFile[5]={"NSYield","NSSigma","Pedestal","ASYield","ASSigma"};

Double_t maxRangePPb[6][5]={
         {4,0.83,13,3,1.28},
         {4,0.83,13,3,1.28},
         {4,0.83,13,3,1.28},
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
Double_t minptData=0.,maxptData=99.;
Int_t ncolumns=4;
Int_t nrows=2;
Double_t ytitleoffset=2.45,xtitleoffset=1.95;
Bool_t skip3to5=kFALSE;
Double_t markersize=1.2;
Double_t markersizeMC=1.2;
Bool_t drawSystMC=kTRUE;
Bool_t runonpPb2016=kTRUE;

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
Bool_t includemodel[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE,kTRUE,kFALSE};
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
  const Int_t nbinsxx=4;
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

void Set2x4PadPositions(TCanvas* c){
    
  TPad * pd1 = (TPad*)c->GetPad(1);
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
  Double_t marginLeft=0.08;
  Double_t marginRight=0.04;
  Double_t marginTop=0.04;
  Double_t marginBottom=0.09;
  Double_t marginLeftForXAxis=0.02;
  Double_t marginBottomForYAxis=0.02;
  innerPadWidth=(1-marginLeft-marginRight)/4.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/2.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!


    // Bottom row
    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd5->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd5->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd5->SetRightMargin(0.);
    pd5->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd5->SetTopMargin(0.);

    pd5->Modified();
    pd5->Update();

    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd6->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2.*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd6->SetRightMargin(0.);
    pd6->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd6->SetTopMargin(0.);

    pd6->Modified();
    pd6->Update();
    
    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd7->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,3.*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd7->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd7->SetRightMargin(0.);
    pd7->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd7->SetTopMargin(0.);

    pd7->Modified();
    pd7->Update();

    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd8->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd8->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd8->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd8->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd8->SetTopMargin(0.);

    pd8->Modified();
    pd8->Update();    

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
    pd3->SetPad(2.*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,3.*innerPadWidth+marginLeft,1);
    //    pd3->SetLeftMargin(0.);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd3->SetRightMargin(0.);
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd3->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(3.*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    //    pd3->SetLeftMargin(0.);
    pd4->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd4->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd4->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();
    

}


TCanvas* CompareNSyieldPPtoPPb(Int_t binass){
  return ComparePPbVsCent(binass,0);
}

TCanvas* CompareNSsigmaPPtoPPb(Int_t binass){
  return  ComparePPbVsCent(binass,1);

}

TCanvas* ComparePedestalPPtoPPb(Int_t binass){
  return ComparePPbVsCent(binass,2);
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


TLatex *GetSystemtext(Int_t identifier){
  
  TLatex *tlSyst;
  Double_t x=0.015,y=0.4;
  
  tlSyst=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"pp, #sqrt{#it{s}} = 13 TeV"); 

  tlSyst->SetNDC();
  tlSyst->SetTextFont(43);
  tlSyst->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  
  tlSyst->SetTextAlign(12);
  tlSyst->SetName(Form("tlSystem_%d",identifier));
  return tlSyst;
}


TLatex *GetRaptext(Int_t identifier){
  
  TLatex *tlRap;
  Double_t x=0.015,y=0.28;
  
  tlRap=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"|#it{y}^{D}_{cms}| < 0.5, |#Delta#eta| < 1"); 

  tlRap->SetNDC();
  tlRap->SetTextFont(43);
  tlRap->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  
  tlRap->SetTextAlign(12);
  tlRap->SetName(Form("tlRap_%d",identifier));
  return tlRap;
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
    x=0.17;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.25;
    x=0.20;// draft 2 was not present -> above value 0.21
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.385;// draft 2 was 0.39
    x=0.16;// draft 2 was not present -> above value 0.21
  }

  if(style==-1){
    alice=new TLatex(0.55,0.85,"ALICE");
    alice->SetNDC();
    alice->SetTextFont(42);
    alice->SetTextSize(0.03);
    alice->SetTextAlign(11);
  }
  else{
    alice= new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"ALICE"); 
    alice->SetNDC();
    alice->SetTextFont(43);
    alice->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);// draft 2 was: 28 *...
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
    tlasspt->SetTextSize(0.025);
  }
  else{
    TString strTot=strPtAssCanvas[binassoc];
    if(addDEta==1)strTot.Append(", |#Delta#eta| < 1");
    tlasspt= new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),strTot.Data()); 
    tlasspt->SetNDC();
    tlasspt->SetTextFont(43);
    tlasspt->SetTextSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);//  if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor) but see notes on top
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

TH1D *GetAndPreparePPb(Int_t binass,Int_t quantity,Int_t numsyst,TGraphAsymmErrors *&gr, TGraphAsymmErrors *&grV2=0){


  TFile *f=TFile::Open(Form("%s/Trends_pp/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
  TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  gr=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  gr->SetName(Form("%sPPb",gr->GetName()));
  for(Int_t iPoint=0;iPoint<4;iPoint++) gr->SetPointError(iPoint,0.7*gr->GetErrorXlow(iPoint),0.7*gr->GetErrorXhigh(iPoint),gr->GetErrorYlow(iPoint),gr->GetErrorYhigh(iPoint));
  if(plotv2unc==kTRUE) {
      grV2=(TGraphAsymmErrors*)c->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));
      grV2->SetName(Form("%sPPb",grV2->GetName()));
      for(Int_t iPoint=0;iPoint<4;iPoint++) grV2->SetPointError(iPoint,0.65*grV2->GetErrorXlow(iPoint),0.65*grV2->GetErrorXhigh(iPoint),grV2->GetErrorYlow(iPoint),grV2->GetErrorYhigh(iPoint));
  }

  TH1D *hPPb=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  hPPb->SetName(Form("%sPPb",hPPb->GetName()));
  hPPb->SetLineColor(colSystem[numsyst]);
  hPPb->SetLineWidth(2);
  hPPb->SetMarkerColor(colSystem[numsyst]);
  hPPb->SetMarkerStyle(markerStyle[numsyst]);
  hPPb->SetMarkerSize(markersize);
  gr->SetMarkerColor(colSystem[numsyst]);
  gr->SetLineColor(colSystem[numsyst]);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(markerStyle[numsyst]);
  gr->SetMarkerSize(markersize);

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
    hPPb->GetYaxis()->SetRangeUser(0,maxRangePPb[binass][quantity]);
  }

  AdaptRangeHist(hPPb,minptData,maxptData);
  AdaptRangeTGraph(gr,minptData,maxptData);

  return hPPb;
}


TCanvas* ComparePPbVsCent(Int_t binass,Int_t quantity,TPad *pd=0x0,Int_t textlegendOptions=0){

  TGraphAsymmErrors *grPPb_1, *grPPbV2_1;
  TH1D *hPPb_1=GetAndPreparePPb(binass,quantity,0,grPPb_1,grPPbV2_1);

  TGraphAsymmErrors *grPPb_2, *grPPbV2_2;
  TH1D *hPPb_2=GetAndPreparePPb(binass,quantity,1,grPPb_2,grPPbV2_2);

  TGraphAsymmErrors *grPPb_3, *grPPbV2_3;
  TH1D *hPPb_3=GetAndPreparePPb(binass,quantity,2,grPPb_3,grPPbV2_3);

  TGraphAsymmErrors *grPPb_4, *grPPbV2_4;
  TH1D *hPPb_4=GetAndPreparePPb(binass,quantity,3,grPPb_4,grPPbV2_4);

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
    hPPb_1->Draw();
  }
  else {
    hDraw=new TH2D(Form("hDraw%d",10*quantity+binass),"",100,0,28,200,0,10);
    hDraw->GetYaxis()->SetTitle("");      
    //    hPP->GetYaxis()->SetTitle("");      

    hDraw->GetXaxis()->SetRangeUser(0,27);
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
  //  hPPb_1->Draw("same");
  }

  //  hPP->Draw("E0X0");// to avoid plotting the error along x
  if(style>0){
    hDraw->GetYaxis()->SetRangeUser(0,maxRangePPb[2][quantity]);
  }
  else {
    hPPb_1->GetYaxis()->SetRangeUser(0,maxRangePPb[binass][quantity]);
  }
 
//  grPPb_1->Draw("E2");
//  grPPbV2_1->Draw("E2");

//  hPPb->Draw("same");

// p-Pb: displace the syst errors along x axis (Fabio)
  Double_t shift=0.2;
  Double_t binsx=grPPb_1->GetN();
  for(int i=0;i<binsx;i++) {
    Double_t x,y,ex1,ex2,ey1,ey2;
    
    grPPb_1->GetPoint(i,x,y);
    ex1=grPPb_1->GetErrorXlow(i); ex2=grPPb_1->GetErrorXhigh(i); ey1=grPPb_1->GetErrorYlow(i); ey2=grPPb_1->GetErrorYhigh(i);
    grPPb_1->SetPoint(i,x+shift*(-1),y);
    grPPb_1->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPbV2_1->GetPoint(i,x,y);
    ex1=grPPbV2_1->GetErrorXlow(i); ex2=grPPbV2_1->GetErrorXhigh(i); ey1=grPPbV2_1->GetErrorYlow(i); ey2=grPPbV2_1->GetErrorYhigh(i);
    grPPbV2_1->SetPoint(i,x+shift*(-1),y);
    grPPbV2_1->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPb_2->GetPoint(i,x,y);
    ex1=grPPb_2->GetErrorXlow(i); ex2=grPPb_2->GetErrorXhigh(i); ey1=grPPb_2->GetErrorYlow(i); ey2=grPPb_2->GetErrorYhigh(i);
    grPPb_2->SetPoint(i,x+shift*(0),y);
    grPPb_2->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPbV2_2->GetPoint(i,x,y);
    ex1=grPPbV2_2->GetErrorXlow(i); ex2=grPPbV2_2->GetErrorXhigh(i); ey1=grPPbV2_2->GetErrorYlow(i); ey2=grPPbV2_2->GetErrorYhigh(i);
    grPPbV2_2->SetPoint(i,x+shift*(0),y);
    grPPbV2_2->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPb_3->GetPoint(i,x,y);
    ex1=grPPb_3->GetErrorXlow(i); ex2=grPPb_3->GetErrorXhigh(i); ey1=grPPb_3->GetErrorYlow(i); ey2=grPPb_3->GetErrorYhigh(i);
    grPPb_3->SetPoint(i,x+shift*(1),y);
    grPPb_3->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPbV2_3->GetPoint(i,x,y);
    ex1=grPPbV2_3->GetErrorXlow(i); ex2=grPPbV2_3->GetErrorXhigh(i); ey1=grPPbV2_3->GetErrorYlow(i); ey2=grPPbV2_3->GetErrorYhigh(i);
    grPPbV2_3->SetPoint(i,x+shift*(1),y);
    grPPbV2_3->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPb_4->GetPoint(i,x,y);
    ex1=grPPb_4->GetErrorXlow(i); ex2=grPPb_4->GetErrorXhigh(i); ey1=grPPb_4->GetErrorYlow(i); ey2=grPPb_4->GetErrorYhigh(i);
    grPPb_4->SetPoint(i,x+shift*(2),y);
    grPPb_4->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPbV2_4->GetPoint(i,x,y);
    ex1=grPPbV2_4->GetErrorXlow(i); ex2=grPPbV2_4->GetErrorXhigh(i); ey1=grPPbV2_4->GetErrorYlow(i); ey2=grPPbV2_4->GetErrorYhigh(i);
    grPPbV2_4->SetPoint(i,x+shift*(2),y);
    grPPbV2_4->SetPointError(i,ex1,ex2,ey1,ey2);
  }

  grPPb_1->Draw("E2");
  grPPbV2_1->SetFillColor(kBlue);
  grPPbV2_1->Draw("E2");

  grPPb_2->Draw("E2");
  grPPbV2_2->SetFillColor(kRed);
  grPPbV2_2->Draw("E2");

  grPPb_3->Draw("E2");
  grPPbV2_3->SetFillColor(kGreen+2);
  grPPbV2_3->Draw("E2");

  grPPb_4->Draw("E2");
  grPPbV2_4->SetFillColor(kBlack);
  grPPbV2_4->Draw("E2");

 /*TH1D* hPPbSuperimp = hPPb->Clone();
  hPPbSuperimp->SetMarkerStyle(25);
  hPPbSuperimp->SetMarkerColor(kRed+1);
//  hPPbSuperimp->Draw("same");*/

//Conversion of pPb TH1F to TGraph to displace the points along x axis (Fabio)
  TGraphAsymmErrors *gr_points_PPb1, *gr_pointCount_PPb1, *gr_points_PPb2, *gr_pointCount_PPb2, *gr_points_PPb3, *gr_pointCount_PPb3, *gr_points_PPb4, *gr_pointCount_PPb4;
    ConvertTH1ToTGraphAsymmError2016(hPPb_1,gr_points_PPb1,shift*(-1));
    ConvertTH1ToTGraphAsymmError2016(hPPb_1,gr_pointCount_PPb1,shift*(-1));    
    ConvertTH1ToTGraphAsymmError2016(hPPb_2,gr_points_PPb2,shift*(0));
    ConvertTH1ToTGraphAsymmError2016(hPPb_2,gr_pointCount_PPb2,shift*(0));    
    ConvertTH1ToTGraphAsymmError2016(hPPb_3,gr_points_PPb3,shift*(1));
    ConvertTH1ToTGraphAsymmError2016(hPPb_3,gr_pointCount_PPb3,shift*(1));    
    ConvertTH1ToTGraphAsymmError2016(hPPb_4,gr_points_PPb4,shift*(2));
    ConvertTH1ToTGraphAsymmError2016(hPPb_4,gr_pointCount_PPb4,shift*(2));  
  
  gr_points_PPb1->SetLineColor(colSystem[0]);
  gr_points_PPb1->SetLineWidth(2);
  gr_points_PPb1->SetMarkerColor(colSystem[0]);
  gr_points_PPb1->SetMarkerStyle(markerStyle[0]);
  gr_points_PPb1->SetMarkerSize(markersize);
  gr_points_PPb1->Draw("samePZ");
  
  gr_pointCount_PPb1->SetLineColor(colSystem[0]);
  gr_pointCount_PPb1->SetLineWidth(2);
  gr_pointCount_PPb1->SetMarkerStyle(colSystem[0]);
  gr_pointCount_PPb1->SetMarkerColor(kBlue);
  gr_pointCount_PPb1->SetMarkerSize(markersize);
  gr_pointCount_PPb1->Draw("samePZ");

  gr_points_PPb2->SetLineColor(colSystem[1]);
  gr_points_PPb2->SetLineWidth(2);
  gr_points_PPb2->SetMarkerColor(colSystem[1]);
  gr_points_PPb2->SetMarkerStyle(markerStyle[1]);
  gr_points_PPb2->SetMarkerSize(markersize);
  gr_points_PPb2->Draw("samePZ");
  
  gr_pointCount_PPb2->SetLineColor(colSystem[1]);
  gr_pointCount_PPb2->SetLineWidth(2);
  gr_pointCount_PPb2->SetMarkerStyle(colSystem[1]);
  gr_pointCount_PPb2->SetMarkerColor(kRed);
  gr_pointCount_PPb2->SetMarkerSize(markersize);
  gr_pointCount_PPb2->Draw("samePZ");

  gr_points_PPb3->SetLineColor(colSystem[2]);
  gr_points_PPb3->SetLineWidth(2);
  gr_points_PPb3->SetMarkerColor(colSystem[2]);
  gr_points_PPb3->SetMarkerStyle(markerStyle[2]);
  gr_points_PPb3->SetMarkerSize(markersize);
  gr_points_PPb3->Draw("samePZ");
  
  gr_pointCount_PPb3->SetLineColor(colSystem[2]);
  gr_pointCount_PPb3->SetLineWidth(2);
  gr_pointCount_PPb3->SetMarkerStyle(colSystem[2]);
  gr_pointCount_PPb3->SetMarkerColor(kGreen+2);
  gr_pointCount_PPb3->SetMarkerSize(markersize);
  gr_pointCount_PPb3->Draw("samePZ");  

  gr_points_PPb4->SetLineColor(colSystem[3]);
  gr_points_PPb4->SetLineWidth(2);
  gr_points_PPb4->SetMarkerColor(colSystem[3]);
  gr_points_PPb4->SetMarkerStyle(markerStyle[3]);
  gr_points_PPb4->SetMarkerSize(markersize);
  gr_points_PPb4->Draw("samePZ");
  
  gr_pointCount_PPb4->SetLineColor(colSystem[3]);
  gr_pointCount_PPb4->SetLineWidth(2);
  gr_pointCount_PPb4->SetMarkerStyle(colSystem[3]);
  gr_pointCount_PPb4->SetMarkerColor(kBlack+1);
  gr_pointCount_PPb4->SetMarkerSize(markersize);
  gr_pointCount_PPb4->Draw("samePZ");   
/*
  if(quantity==1 && binass==2) {
    TLatex *tlDispl=new TLatex(0.20,0.86,"p-Pb points and error boxes");
    TLatex *tlDispl2=new TLatex(0.20,0.80,"shifted by #Delta#it{p}_{T} = +0.3 GeV/#it{c}");
    tlDispl->SetNDC();
    tlDispl->SetTextFont(42);
    tlDispl->SetTextSize(0.048);
    tlDispl->Draw();
    tlDispl2->SetNDC();
    tlDispl2->SetTextFont(42);
    tlDispl2->SetTextSize(0.048);
    tlDispl2->Draw();
  }
*/
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

  /*  TLegend *legendSuperimp=GetLegendDataPointsFake(hPP,hPPbSuperimp,10*quantity+binass);
    legendSuperimp->Draw("same");*/
        
    TLatex *tlAssYieldPt=GetAssocPtText(binass,10*quantity+binass,0);
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
      TLatex *tlAssYieldPt=GetAssocPtText(binass,10*quantity+binass,0);
      tlAssYieldPt->Draw();
    }
    if(textlegendOptions%1000>=100){
      TLegend * legend;
      legend = new TLegend(0.005/gPad->GetWNDC()+gPad->GetLeftMargin(),0.215/gPad->GetHNDC()+gPad->GetBottomMargin(),0.18/gPad->GetWNDC()+gPad->GetLeftMargin(),0.312/gPad->GetHNDC()+gPad->GetBottomMargin());// draft 2 (2 lines only, rapidity on the same line also for p-Pb): 0.002/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.30/gPad->GetHNDC()+gPad->GetBottomMargin()
      legend->SetTextFont(43);
      legend->SetTextAlign(12);
      legend->SetLineColor(kWhite);
      legend->SetTextSize(20*innerPadHeight/referencePadHeight*resizeTextFactor);
      legend->AddEntry(hPPb_1,"0-0.1% V0M","lep");
      legend->AddEntry(hPPb_2,"0.1-10% V0M","lep");
      legend->AddEntry(hPPb_3,"10-30% V0M","lep");
      legend->AddEntry(hPPb_4,"30-100% V0M","lep");
      legend->Draw();

     /*TLegend *legendSuperimp=GetLegendDataPointsFake(hPP,hPPbSuperimp,10*quantity+binass);
      legendSuperimp->Draw("same");*/
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

void CompareFitResults_pPbVsCentrality_UniqueCanvas(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
    cFinalPaperStyle=new TCanvas("cPPvsPPbFitResultsFinalPaperStyle","pp vs.pPb fit results ",1800./900.*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(4,2,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set2x4PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);

  Int_t orderAssoc[4]={0,1,2,3};  
  for(Int_t jp=0;jp<=3;jp++){// First loop for NS yield
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      ComparePPbVsCent(orderAssoc[jp],0,pd,100+needTitle);    
      TLatex *tlSyst=GetSystemtext(orderAssoc[jp]);
      tlSyst->Draw();
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      ComparePPbVsCent(orderAssoc[jp],0,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderAssoc[jp]);
      tlSide->Draw();
    }
    else if(jp==2){// identifier set as: 10*quantity+binass; NS:
      ComparePPbVsCent(orderAssoc[jp],0,pd,needTitle);    
      TLatex *tlRap=GetRaptext(orderAssoc[jp]);
      tlRap->Draw();
    }         
    else{
      ComparePPbVsCent(orderAssoc[jp],0,pd,0+needTitle);    
    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp],0);
    tlAssYieldPt->Draw();
 
    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=3;jp++){// second loop --> sigma
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+5);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    ComparePPbVsCent(orderAssoc[jp],1,pd,needTitle); 
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
  for(Int_t j=8;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_4x2Bins.root");
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_4x2Bins.eps");
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_4x2Bins.png");
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_4x2Bins.pdf");

  return;

}

void CompareFitResults_pPbVsCentrality_UniqueCanvas_AwaySide(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
    cFinalPaperStyle=new TCanvas("cPPvsPPbFitResultsFinalPaperStyle_AwaySide","pp vs.pPb fit results ",1800./900.*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(4,2,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set2x4PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);

  Int_t orderAssoc[4]={0,1,2,3};  
  for(Int_t jp=0;jp<=3;jp++){// First loop for NS yield
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      ComparePPbVsCent(orderAssoc[jp],3,pd,100+needTitle);    
      TLatex *tlSyst=GetSystemtext(orderAssoc[jp]);
      tlSyst->Draw();      
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      ComparePPbVsCent(orderAssoc[jp],3,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderAssoc[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(3,orderAssoc[jp]);
      tlSide->Draw();
    }
    else if(jp==2){// identifier set as: 10*quantity+binass; NS:
      ComparePPbVsCent(orderAssoc[jp],3,pd,needTitle);    
      TLatex *tlRap=GetRaptext(orderAssoc[jp]);
      tlRap->Draw();
    }     
    else{
      ComparePPbVsCent(orderAssoc[jp],3,pd,0+needTitle);    
    }
    TLatex *tlAssYieldPt=GetAssocPtText(orderAssoc[jp],orderAssoc[jp],0);
    tlAssYieldPt->Draw();
 
    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=3;jp++){// second loop --> sigma
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+5);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    ComparePPbVsCent(orderAssoc[jp],4,pd,needTitle); 
    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=8;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_AwaySide_4x2Bins.root");
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_AwaySide_4x2Bins.eps");
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_AwaySide_4x2Bins.png");
  cFinalPaperStyle->SaveAs("ComparePPbVsCentFitResults_AwaySide_4x2Bins.pdf");

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
  
  TCanvas *c=CompareDatatoModels(1,0,0,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,0,1,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,0,2,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,0,3,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,0,4,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));  


  c=CompareDatatoModels(1,1,0,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,1,1,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,1,2,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,1,3,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,1,4,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));


  c=CompareDatatoModels(1,2,0,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,2,1,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,2,2,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,2,3,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,2,4,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

  c=CompareDatatoModels(1,3,0,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,3,1,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,3,2,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,3,3,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,3,4,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

  c=CompareDatatoModels(1,4,0,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,4,1,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,4,2,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,4,3,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,4,4,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));

  c=CompareDatatoModels(1,5,0,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,5,1,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,5,2,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,5,3,0x0,999);
  c->SaveAs(Form("%s.root",c->GetName()));
  c->SaveAs(Form("%s.eps",c->GetName()));
  c->SaveAs(Form("%s.png",c->GetName()));
  c=CompareDatatoModels(1,5,4,0x0,999);
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
