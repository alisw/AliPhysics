//Int_t nhistos=6;
const Int_t ncollsyst=2;
const Int_t nbinDpt=3;
const Int_t nbinAssocpt=3;
Int_t firstDpt=0;
Double_t canvasheight=900;
Double_t ratioForCanvasWidth=2./3.*1.25;
Double_t resizeTextFactor=0.9*900./1000.;//canvasheight/1800.; // size was tuned for canvasheight =1800. 
Double_t markersize=1.*resizeTextFactor;
Int_t style=1;
// FOR ELENA: HERE THERE ARE JUST STRINGS AND NAMES/LEGENDS
TString collsyst[ncollsyst]={"pp","pPb"};
TString pthadron[nbinAssocpt]={"0.3to99.0","0.3to1.0","1.0to99.0"};
TString strmesonpt[nbinDpt]={"3to5","5to8","8to16"};
TString strPtAssocText[nbinAssocpt]={"#it{p}_{T}^{assoc} > 0.3 GeV/#it{c}","0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}","#it{p}_{T}^{assoc} > 1 GeV/#it{c}"};
TString strPtMesonText[nbinDpt]={ "3 < #it{p}_{T}^{D} < 5 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}"};
TString strDeltaEta="|#Delta#eta| < 1";
Double_t minYaxis[nbinAssocpt]={-0.6,-0.6,-0.29}; // or -0.6 for all  
Double_t maxYaxis[nbinAssocpt]={5.2,1.8,1.9};// or 4.1 for the first , 2.9 for the middle and 1.9 for the last
Double_t mesonptcenter[nbinDpt]={4.,6.5,12.};// needed for getting the baseline
Color_t colourSystem[ncollsyst]={kBlack,kRed}; 
Color_t colourSystemBaseUnc[ncollsyst]={kGray+1,kRed-9};//was {kBlack,kRed-7};
Int_t fillColourBaselineStyle=1001;
Int_t markerstyle[2]={20,21};// pp, p-Pb
Double_t innerPadHeight;// not touch, set internally
Double_t innerPadWidth;// not touch, set internally
Int_t rebin =1;
Double_t **ptbinsD;
Int_t ndesiredbins=3;
TString inputdatadirectory = "./";
TString baselinedirectory="./";
TString avType="Weighted";
TString filenames[ncollsyst][nbinAssocpt][nbinDpt];// [coll syst][ptassoc][ptmes]
TString pedestalfilenames[ncollsyst][nbinAssocpt];// [coll syst][ptassoc]
TH1D ****histo;
TGraphAsymmErrors ****err;
TLatex ****ltscale;
TH1D ****subtractedhisto;
TH1D ***histopPbSuperimp;
TGraphAsymmErrors ****suberr;
TGraphAsymmErrors ****grbase;
TGraphAsymmErrors ****grv2; 

Bool_t skip3to5=kTRUE;
Int_t ihskip=0;
TString fitplotmacrodir=gSystem->ExpandPathName("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/");
Bool_t isReflected=kFALSE;
//Double_t leftMarginCanvas=0.17;
//Double_t rightMarginCanvas=0.055;
//Double_t bottomMarginCanvas=0.1;
//Double_t topMarginCanvas=0.1;

Double_t scaleHeightPads=1;// do not touch this, it is regulated automatically in the macro
Double_t scaleWidthPads=1;// do not touch this, it is regulated automatically in the macro
Double_t ytitleoffset=7.5,xtitleoffset=6.;// was 4.8(y) and 3.45(x)
Double_t referencePadHeight=0.2933; // Do not touch unless the canvas division is changed from 2x3, change canvasheight and resizeTextFactor instead
Double_t referencePadHeight=0.44; // Do not touch unless the canvas division is changed from 2x3, change canvasheight and resizeTextFactor instead

void SetSkip3to5pPb(Bool_t skip){
  skip3to5=skip;
}

void SetIsReflected(Bool_t isrefl){
  isReflected=isrefl;
}

void SetBaselineDirectory(TString dirbase){
  baselinedirectory=dirbase;
}

void SetInputDataDirectory(TString inputdir){
  inputdatadirectory=inputdir;
}

void SetAverageMode(Int_t avmode){
  if(avmode==0)avType="Weighted";
  else if(avmode==1)avType="Arithmetic";
  else Printf("DO COMPARISON WITH MC: WRONGE AVERAGE METHOD SET");
}

void SetPaveStyle(TPaveText *pv){
  pv->SetFillColor(0);
  pv->SetFillStyle(0);
  //  pv->SetBorderStyle(0);
  pv->SetBorderSize(0);
  pv->SetTextFont(43);
  //  pv->SetTextSize(20);
  pv->SetTextAlign(12);

}

void SetPadStyle(TPad *p){
  p->SetTicky();
  p->SetTickx();
  //  p->SetFrameBorderMode(0);
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

void LoadFileNamesppVspPb(){

  for(Int_t icoll=0;icoll<ncollsyst;icoll++){
    for(Int_t kassoc=0;kassoc<nbinAssocpt;kassoc++){
      for(Int_t jmes=0;jmes<nbinDpt;jmes++){
	filenames[icoll][kassoc][jmes]=Form("%s/%sAverage%sDzeroDstarDplus%s_assoc%s.root",inputdatadirectory.Data(),avType.Data(),collsyst[icoll].Data(),strmesonpt[jmes].Data(),pthadron[kassoc].Data());
      }
      pedestalfilenames[icoll][kassoc]=Form("%s/Trends_%s/CanvasBaselineVariationTrendPedestal_pthad%s.root",baselinedirectory.Data(),collsyst[icoll].Data(),pthadron[kassoc].Data());
    }
  }
}

void Init(){
  if(skip3to5)firstDpt=1;
  
  histo=new TH1D***[ncollsyst];
  err=new TGraphAsymmErrors***[ncollsyst]; 
  ltscale = new TLatex ***[ncollsyst]; 
  subtractedhisto=new TH1D***[ncollsyst];
  suberr=new TGraphAsymmErrors***[ncollsyst]; 
  grbase=new TGraphAsymmErrors***[ncollsyst]; 
  grv2=new TGraphAsymmErrors***[ncollsyst]; 
  for(Int_t icoll=0;icoll<ncollsyst;icoll++){
    histo[icoll]=new TH1D**[nbinAssocpt];
    err[icoll]= new TGraphAsymmErrors**[nbinAssocpt];
    ltscale[icoll] = new TLatex **[nbinAssocpt];
    subtractedhisto[icoll]=new TH1D**[nbinAssocpt];
    histopPbSuperimp=new TH1D**[nbinAssocpt];
    suberr[icoll]=new TGraphAsymmErrors**[nbinAssocpt]; 
    grbase[icoll]=new TGraphAsymmErrors**[nbinAssocpt]; 
    grv2[icoll]=new TGraphAsymmErrors**[nbinAssocpt]; 
    for(Int_t kassoc=0;kassoc<nbinAssocpt;kassoc++){
      histo[icoll][kassoc]=new TH1D*[nbinDpt];
      subtractedhisto[icoll][kassoc]=new TH1D*[nbinDpt];
      histopPbSuperimp[kassoc]=new TH1D*[nbinDpt];
      err[icoll][kassoc]= new TGraphAsymmErrors*[nbinDpt]; 
      ltscale[icoll][kassoc] = new TLatex*[nbinDpt]; 
      suberr[icoll][kassoc]=new TGraphAsymmErrors*[nbinDpt]; 
      grbase[icoll][kassoc]=new TGraphAsymmErrors*[nbinDpt]; 
      grv2[icoll][kassoc]=new TGraphAsymmErrors*[nbinDpt]; 
      for(Int_t jmes=0;jmes<nbinDpt;jmes++){
      histo[icoll][kassoc][jmes]=0x0;
      subtractedhisto[icoll][kassoc][jmes]=0x0;
      histopPbSuperimp[kassoc][jmes]=0x0;
      err[icoll][kassoc][jmes]=0x0; 
      ltscale[icoll][kassoc][jmes] =0x0; 
      suberr[icoll][kassoc][jmes]=0x0; 
      grbase[icoll][kassoc][jmes]=0x0; 
      grv2[icoll][kassoc][jmes]=0x0; 
      }
    }
  }
}

Int_t rebin =1;

//_______________________________________________________________________
TH1D * GetHistoAndSyst(TString path, Int_t collsyst,TString hname, TString hnamesyst,TGraphAsymmErrors *&gr2, TLatex *&tUncertainty){
    
    
  Printf("Opening file: %s",path.Data());
  TFile *f=TFile::Open(path.Data(),"READ");
  TH1D *hFDsub=(TH1D*)f->Get(hname.Data());
  AliHFDhadronCorrSystUnc *syst=(AliHFDhadronCorrSystUnc*)f->Get(hnamesyst.Data());
  TH1D *hUncCorrMin=syst->GetHistoTotFlatMin();
  TH1D *hUncCorrMax=syst->GetHistoTotFlatMax();

  gr2=syst->GetTotNonFlatUncGraph();
  gr2->SetLineColor(kBlack);
  gr2->SetMarkerColor(kBlack);
  gr2->SetFillStyle(0);

  if(collsyst==0){
    if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.55,0.46,Form("%.0f#% scale uncertainty pp",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.55,0.46,Form("{}^{#plus%.0f%s}_{#minus%.0f%s} scale uncertainty (pp)","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
   
  }
  if(collsyst==1){
    if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.55,0.42,Form("#bf{%.0f#% scale uncertainty p-Pb}",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.55,0.42,Form("{}^{#plus%.0f%s}_{#minus%.0f%s} scale uncertainty (p-Pb)","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
  }

  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.025);
  tUncertainty->SetTextFont(62);
    
  return hFDsub;
    
}



//_______________________________________________________________________
TH1D * GetPedestalHistoAndSystAndSubtractPedpPb(Int_t binSystem, Int_t binAssoc,Int_t binMeson,TH1D *histo, TGraphAsymmErrors* gr,TGraphAsymmErrors *&grout, TString canvasname,TGraphAsymmErrors *&grbaseOut,TGraphAsymmErrors *&grv2Out){

  Double_t value = 0, pedestal=0;
  TGraphAsymmErrors *grBaseHelp,*grV2;
  grbaseOut=new TGraphAsymmErrors();
  grbaseOut->SetName(Form("grbaselineUncFull_%s_%s_%s",collsyst[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data()));

  Double_t xuncFull,errxuncFull;
  Double_t xuncv2,errxuncv2;
  Int_t bin,bingr;
  TString path = pedestalfilenames[binSystem][binAssoc];//pPb, proper bin assoc

  if(binSystem==1){
    // get pedestal from fit outputs
    cout << "pPb -->  Reading File from path: " << path << endl;
    
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    grBaseHelp=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fBaselineVariationSystematicsPedestal");    
    grv2Out=new TGraphAsymmErrors();
    grv2Out->SetName(Form("grbaselineUncV2_%s_%s_%s",collsyst[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data()));

    grV2=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fv2SystematicsPedestal");    

    if(isReflected){
      xuncFull=3.4;
      errxuncFull=0.075;
      xuncv2=-0.3;
      errxuncv2=0.075;
    }
    else{
      xuncFull=4.98;
      errxuncFull=0.075;
      xuncv2=-1.85;
      errxuncv2=0.075;
    }
    bin=h->FindBin(mesonptcenter[binMeson]);
    bingr=GetBinGraph(mesonptcenter[binMeson],grBaseHelp);
  
    pedestal=h->GetBinContent(bin);
    Double_t x,y,erryl,erryh;
    grBaseHelp->GetPoint(bingr,x,y);
    Printf("histo: x=%f, graph: %f",h->GetBinCenter(bin),x);
    erryl=grBaseHelp->GetErrorYlow(bingr);
    erryh=grBaseHelp->GetErrorYhigh(bingr);
    
    grbaseOut->SetPoint(0,xuncFull,0);
    grbaseOut->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    if(grV2){
      grV2->GetPoint(bingr,x,y);
      erryl=grV2->GetErrorYlow(bingr);
      erryh=grV2->GetErrorYhigh(bingr);	
      grv2Out->SetPoint(0,xuncFull,0);
      grv2Out->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    }
    
  }
    

  if(binSystem==0){
    grV2=0x0;
    cout << "pp -->  Reading File from path: " << path << endl;
    
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    grBaseHelp=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fBaselineVariationSystematicsPedestal");    

    if(isReflected){
      xuncFull=3.25;
      errxuncFull=0.075;
      xuncv2=-0.15; 
      errxuncv2=0.075;
    }
    else{
      xuncFull=4.83;
      errxuncFull=0.075;
      xuncv2=-1.7;
      errxuncv2=0.075;
    }

    bin=h->FindBin(mesonptcenter[binMeson]);
    bingr=GetBinGraph(mesonptcenter[binMeson],grBaseHelp);


    pedestal=h->GetBinContent(bin);
    Double_t x,y,erryl,erryh;
    Printf("histo: x=%f, graph: %f",h->GetBinCenter(bin),x);
    grBaseHelp->GetPoint(bingr,x,y);
    erryl=grBaseHelp->GetErrorYlow(bingr);
    erryh=grBaseHelp->GetErrorYhigh(bingr);
    
    grbaseOut->SetPoint(0,xuncFull,0);
    grbaseOut->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    if(grV2){
      grV2->GetPoint(bingr,x,y);
      erryl=grV2->GetErrorYlow(bingr);
      erryh=grV2->GetErrorYhigh(bingr);	
      grv2Out->SetPoint(0,xuncFull,0);
      grv2Out->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    }
  }

  cout<<"Baseline being subtracted for:"<<collsyst[binSystem].Data()<<strmesonpt[binMeson].Data()<<pthadron[binAssoc].Data()<<endl;
  
  grout=(TGraphAsymmErrors*)gr->Clone(Form("grSub_%s_%s_%s",collsyst[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data()));
  
  TString nameoutput = histo->GetName();
  nameoutput += "_subtr_";
  nameoutput += "pedestal";
  nameoutput += Form("_%s_%s_%s",collsyst[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data());
   
  TH1D * outputhisto = (TH1D*)histo->Clone(nameoutput.Data());
  outputhisto->Reset();
  outputhisto->SetStats(kFALSE);
    

  for(Int_t iBin = 1; iBin <= histo->GetNbinsX();iBin++){
        
        
    value = histo->GetBinContent(iBin);
    value -= pedestal;
      
    outputhisto->SetBinContent(iBin,value);
      
    outputhisto->SetBinError(iBin,histo->GetBinError(iBin));
    Double_t x,y,eyl,eyh;
    gr->GetPoint(iBin-1,x,y);
    eyl=gr->GetErrorYlow(iBin-1);
    eyh=gr->GetErrorYhigh(iBin-1);
    //cout<<x<<"  "<<y<<"  "<<eyl<<"  "<<eyh<<endl;
    grout->SetPoint(iBin-1,x,y-pedestal);

  }
  cout<<"sub -> "<<outputhisto->GetBinContent(5)<<endl;
  
  outputhisto->SetXTitle("#Delta#varphi (rad)");
  outputhisto->SetYTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} - baseline (rad^{-1})");
  //  outputhisto->GetYaxis()->SetTitleOffset(1.5);
  //  outputhisto->GetYaxis()->SetTitleFont(42);
  //  outputhisto->GetXaxis()->SetTitleFont(42);
  //  outputhisto->GetYaxis()->SetTitleSize(25);
  //  outputhisto->GetXaxis()->SetTitleSize(25);
  //  outputhisto->GetYaxis()->SetLabelSize(0.04);
  //  outputhisto->GetXaxis()->SetLabelSize(0.04);

  cout<<"now return histogram"<<endl;

  return outputhisto;
    
}


Int_t GetBinGraph(Double_t xbincenter, TGraph *gr,Double_t tolerance=0.1){// works only if the bin center is passed (or a value within tolerance). In other cases might not giving back what is expected
  Int_t size=gr->GetN();
  Double_t xg,yg;
  for(Int_t j=0;j<size;j++){
    gr->GetPoint(j,xg,yg);
    if(j==0){
      if(xbincenter<xg&&TMath::Abs(xg-xbincenter)>tolerance)return -1;      
    }

    if(TMath::Abs(xg-xbincenter)<tolerance)return j;
    
  }
  return -1;
}


void Set3x2PadPositions(TCanvas* c){
    
  TPad * pd1 = (TPad*)c->GetPad(1);
    TPad * pd2 =(TPad*) c->GetPad(2);
    TPad * pd3 =(TPad*) c->GetPad(3);
    TPad * pd4 =(TPad*) c->GetPad(4);
    TPad * pd5 =(TPad*) c->GetPad(5);
    TPad * pd6 =(TPad*) c->GetPad(6);
    
    SetPadStyle(pd1);

    Double_t xl,xu,yl,yu;
    Double_t marginLeft=0.1;
    Double_t marginRight=0.02;
    Double_t marginTop=0.04;
    Double_t marginBottom=0.08;
    innerPadWidth=(1-marginLeft-marginRight)/2.;// this is the width w/o margin, not the real pad width!!
    innerPadHeight=(1-marginTop-marginBottom)/3.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!
    Printf("innerPadHeight: %f",innerPadHeight);
    Printf("innerPadWidth: %f",innerPadWidth);
    Double_t marginLeftForXAxis=0.02;
    Double_t marginBottomForYAxis=0.02;

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
    pd6->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));//0.02/(1.-innerPadWidth-marginLeft));
    pd6->SetRightMargin(marginRight/(innerPadWidth+marginRight+marginLeftForXAxis));
    pd6->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd6->SetTopMargin(0.);

        pd6->Modified();
        pd6->Update();


    // Middle Row
    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd3->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd3->SetRightMargin(0.);
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd3->SetTopMargin(0.);

        pd3->Modified();
        pd3->Update();

    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1.,2.*innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));//0.02/(1.-innerPadWidth-marginLeft));
    pd4->SetRightMargin(marginRight/(innerPadWidth+marginRight+marginLeftForXAxis));
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd4->SetTopMargin(0.);

        pd4->Modified();
        pd4->Update();


    // Top Row
    pd1->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd1->SetPad(0,2.*innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,1.);
    pd1->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd1->SetRightMargin(0.);
    pd1->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd1->SetTopMargin(marginTop/(innerPadHeight+marginTop+marginBottomForYAxis));

    pd1->Modified();
    pd1->Update();
    
   
    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();

    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));//0.02/(1.-innerPadWidth-marginLeft));
    pd2->SetRightMargin(marginRight/(innerPadWidth+marginRight+marginLeftForXAxis));
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd2->SetTopMargin(marginTop/(innerPadHeight+marginTop+marginBottomForYAxis));
    

}

void SetColorGraphs(Int_t system,TGraph* gr,Color_t color=-1){
  if(color==-1){
    gr->SetLineColor(colourSystem[system]);
    gr->SetFillColor(colourSystem[system]);
    gr->SetMarkerColor(colourSystem[system]);
  }
  else{
    gr->SetLineColor(color);
    gr->SetFillColor(color);
    gr->SetMarkerColor(color);
  }

  
}

void SetScaleUncertaintyPositionAndSize(TLatex *tlpp,TLatex *tlpPb){
  tlpp->SetTextFont(43);
  tlpPb->SetTextFont(43);
  /* old settings
 tlpp->SetX(0.2/gPad->GetWNDC()+gPad->GetLeftMargin());
  tlpp->SetY(0.19/gPad->GetHNDC()+gPad->GetBottomMargin());
  tlpPb->SetX(0.2/gPad->GetWNDC()+gPad->GetLeftMargin());
  tlpPb->SetY(0.156/gPad->GetHNDC()+gPad->GetBottomMargin());
  */

  //  tlpp->SetX(0.13/gPad->GetWNDC()+gPad->GetLeftMargin()); //draft 1
  //   tlpp->SetY(0.19/gPad->GetHNDC()+gPad->GetBottomMargin());
  tlpp->SetY(0.21/gPad->GetHNDC()+gPad->GetBottomMargin());
  tlpp->SetX(0.083/gPad->GetWNDC()+gPad->GetLeftMargin());

  //  tlpPb->SetX(0.13/gPad->GetWNDC()+gPad->GetLeftMargin()); // draft 1
  //  tlpPb->SetY(0.156/gPad->GetHNDC()+gPad->GetBottomMargin());
  tlpPb->SetX(0.083/gPad->GetWNDC()+gPad->GetLeftMargin());
  tlpPb->SetY(0.176/gPad->GetHNDC()+gPad->GetBottomMargin());

  tlpp->SetTextSize(26*innerPadHeight/referencePadHeight*resizeTextFactor);// old settings with font 42: 0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  tlpPb->SetTextSize(26*innerPadHeight/referencePadHeight*resizeTextFactor);// old settings with font 42: 0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  tlpp->SetTextAlign(12);
  tlpPb->SetTextAlign(12);

  return;
}

TLegend *GetLegendDataPoints(TH1D *hpp,TH1D *hpPb,Int_t identifier){
    
  TLegend * legend = new TLegend(0.008/gPad->GetWNDC()+gPad->GetLeftMargin(),0.18/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.245/gPad->GetHNDC()+gPad->GetBottomMargin());
    //new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.14/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.2/gPad->GetHNDC()+gPad->GetBottomMargin());
  //previous settings: (0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.16/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.22/gPad->GetHNDC()+gPad->GetBottomMargin());
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(25*innerPadHeight/referencePadHeight*resizeTextFactor);
			// old settings with font 42, still good (0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
    legend->AddEntry(hpp,"pp, #sqrt{#it{s}} = 7 TeV, |#it{y}^{D}_{cms}| < 0.5","lep");
    legend->AddEntry((TObject*)0,"EPJC 77 (2017) 245","");
    legend->AddEntry(hpPb,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, -0.96 < #it{y}^{D}_{cms} < 0.04","lep");
    legend->SetName(Form("LegendDataPPandpPb_%d",identifier));
    return legend;
  }

TLegend *GetLegendDataPointsFake(TH1D *hpp,TH1D *hpPb,Int_t identifier){
    
  TLegend * legend = new TLegend(0.008/gPad->GetWNDC()+gPad->GetLeftMargin(),0.18/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.245/gPad->GetHNDC()+gPad->GetBottomMargin());
    //new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.14/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.2/gPad->GetHNDC()+gPad->GetBottomMargin());
  //previous settings: (0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.16/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.22/gPad->GetHNDC()+gPad->GetBottomMargin());
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(25*innerPadHeight/referencePadHeight*resizeTextFactor);
			// old settings with font 42, still good (0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
    legend->AddEntry(hpp,"","lep");
    legend->AddEntry((TObject*)0,"","");
    legend->AddEntry(hpPb,"","lep");
    legend->SetName(Form("LegendDataPPandpPb_%d",identifier));
    return legend;
  }

TLegend *GetLegendBaselines(TGraphAsymmErrors *gpp,TGraphAsymmErrors *gpPb,Int_t identifier,Int_t positionOption=1){
  TLegend * legend;
  if(positionOption==0){
    legend=new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.21/gPad->GetHNDC()+gPad->GetBottomMargin(),0.10/gPad->GetWNDC()+gPad->GetLeftMargin(),0.275/gPad->GetHNDC()+gPad->GetBottomMargin());
  }
  else if(positionOption==1){
    legend= new TLegend(0.20/gPad->GetWNDC()+gPad->GetLeftMargin(),0.21/gPad->GetHNDC()+gPad->GetBottomMargin(),0.31/gPad->GetWNDC()+gPad->GetLeftMargin(),0.275/gPad->GetHNDC()+gPad->GetBottomMargin());
  }
  else if(positionOption==2){// standard
    legend= new TLegend(0.13/gPad->GetWNDC()+gPad->GetLeftMargin(),0.17/gPad->GetHNDC()+gPad->GetBottomMargin(),0.23/gPad->GetWNDC()+gPad->GetLeftMargin(),0.235/gPad->GetHNDC()+gPad->GetBottomMargin());
  }
  else if(positionOption==3){
    legend= new TLegend(0.035/gPad->GetWNDC()+gPad->GetLeftMargin(),0.18/gPad->GetHNDC()+gPad->GetBottomMargin(),0.155/gPad->GetWNDC()+gPad->GetLeftMargin(),0.245/gPad->GetHNDC()+gPad->GetBottomMargin());
  }
  else return 0x0;
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextAlign(12);
  legend->SetTextSize(23*innerPadHeight/referencePadHeight*resizeTextFactor);// old settings with font 42, 0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  legend->AddEntry(gpp,"baseline-subtraction uncertainty (pp)","f");
  legend->AddEntry(gpPb,"baseline-subtraction uncertainty (p-Pb)","f");
  // legend->AddEntry(box2[4],"v2 subtr p-Pb","f");
  legend->SetName(Form("LegendBaselineUncPPandpPb_%d",identifier));
  return legend;
}


TPaveText *GetALICEpavetext(Int_t identifier){
  TPaveText *alice = new TPaveText(0.33/gPad->GetWNDC()+gPad->GetLeftMargin(),0.255/gPad->GetHNDC()+gPad->GetBottomMargin(),0.37/gPad->GetWNDC()+gPad->GetLeftMargin(),0.275/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  //  TPaveText *alice = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  SetPaveStyle(alice);
  alice->SetTextFont(43);
  alice->SetTextSize(32*innerPadHeight/referencePadHeight*resizeTextFactor);//old settings with font 42 : 0.08/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  alice->AddText("ALICE");//commented
  // fitvalueslow->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  alice->SetName(Form("paveALICE_%d",identifier));
  return alice;
}

TPaveText *GetAveragepavetext(Int_t identifier){
  TPaveText* pvAverage = new TPaveText(0.0107/gPad->GetWNDC()+gPad->GetLeftMargin(),0.255/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.275/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  //  TPaveText* pvAverage = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.25/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  //0.21,0.83,0.5,0.863,"NDC");
  SetPaveStyle(pvAverage);
  pvAverage->SetTextFont(43);
  pvAverage->SetTextSize(32*innerPadHeight/referencePadHeight*resizeTextFactor);//old settings with font 42: 0.08/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  pvAverage->AddText("Average D^{0}, D^{+}, D^{*+}");
  pvAverage->SetName(Form("paveAverage_%d",identifier));
  return pvAverage;
}


TPaveText *GetPaveKineInfo(TString strPtMeson,TString strPtAssoc,Int_t identifier,Int_t mode=10/* mode: alignment * 10 + 0,1,2,3 (two rows, two rows + include DeltaEta, single row, single row + included delta eta, with alignment=0 -> align left, 1->align right*/){
  TPaveText *pvKineInfo;
  Double_t extraSize=0.;
  if(mode%10==1){
    extraSize=0.02;
  }
  if(mode%10==1){
    extraSize=0.02;
  }
  if(mode%10==2){
    extraSize=0.218;
  }
  if(mode%10==3){
    extraSize=0.238;
  }
  
  if(mode/10==0){
    pvKineInfo= new TPaveText(0.009/gPad->GetWNDC()+gPad->GetLeftMargin(),0.215/gPad->GetHNDC()+gPad->GetBottomMargin(),0.23/gPad->GetWNDC()+gPad->GetLeftMargin()+extraSize,0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
    SetPaveStyle(pvKineInfo);
    pvKineInfo->SetTextAlign(12);
  }
  else if(mode/10==2){
    pvKineInfo= new TPaveText(0.025/gPad->GetWNDC()+gPad->GetLeftMargin(),0.215/gPad->GetHNDC()+gPad->GetBottomMargin(),0.42/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
    SetPaveStyle(pvKineInfo);
    pvKineInfo->SetTextAlign(22);
  }
  else if(mode/10==1){
    pvKineInfo= new TPaveText(0.32/gPad->GetWNDC()+gPad->GetLeftMargin()-extraSize,0.215/gPad->GetHNDC()+gPad->GetBottomMargin(),0.42/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
    SetPaveStyle(pvKineInfo);
    pvKineInfo->SetTextAlign(32);
  }
  else return 0x0;
  TString strname=Form("pvKineInfo");
  TString strall;
  pvKineInfo->SetTextFont(43);
  pvKineInfo->SetTextSize(26.*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  if(!strPtMeson.IsNull()){
    if(mode%10<2)pvKineInfo->AddText(strPtMeson.Data());
    else {
      strall=strPtMeson;
      if(strPtAssoc.IsNull()){
	pvKineInfo->AddText(strall.Data());
      }
    }
    strname.Append("Dpt");
  }
  if(!strPtAssoc.IsNull()){
    strname.Append("AssocPt");
    if(mode%10==1||mode%10==3){
      strPtAssoc.Append(Form(", %s",strDeltaEta.Data()));
      strname.Append("DeltaEta");
    }
    if(mode%10<2){
      pvKineInfo->AddText(strPtAssoc.Data());
    }
    else{
      strall.Append(Form(", %s",strPtAssoc.Data()));
      pvKineInfo->AddText(strall.Data());
    }
  }
  else if(strPtMeson.IsNull()&&(mode%10==1||mode%10==3)){
    strname.Append("DeltaEta");
    pvKineInfo->AddText(strDeltaEta.Data());
  }
  pvKineInfo->SetName(Form("%s_%d",strname.Data(),identifier));
  return pvKineInfo;
}

void DoComparison_ppVspPballPanels2016(){
  gStyle->SetOptStat(0000);
  Init();
  LoadFileNamesppVspPb();
  
  for(Int_t icoll=0;icoll<ncollsyst;icoll++){
    for(Int_t kassoc=0;kassoc<nbinAssocpt;kassoc++){
      for(Int_t jmes=firstDpt;jmes<nbinDpt;jmes++){
	//	if(skip3to5&&jmes==0)continue;
	//	TH1D * GetHistoAndSyst(TString path, Int_t collsyst,TString hname, TString hnamesyst,TGraphAsymmErrors *&gr2, TLatex *&tUncertainty){
	histo[icoll][kassoc][jmes]=GetHistoAndSyst(filenames[icoll][kassoc][jmes],icoll,"fhDaverage","AverageSystematicUncertainty",err[icoll][kassoc][jmes],ltscale[icoll][kassoc][jmes]); 
	

	if(icoll==0) { // pp style
	  histo[icoll][kassoc][jmes]->SetMarkerColor(colourSystem[icoll]); 
	  histo[icoll][kassoc][jmes]->SetMarkerStyle(markerstyle[icoll]); 
	  histo[icoll][kassoc][jmes]->SetMarkerSize(markersize); 
	  histo[icoll][kassoc][jmes]->SetLineColor(colourSystem[icoll]);  
	  

	  subtractedhisto[icoll][kassoc][jmes] = GetPedestalHistoAndSystAndSubtractPedpPb(icoll,kassoc,jmes,histo[icoll][kassoc][jmes],err[icoll][kassoc][jmes],suberr[icoll][kassoc][jmes],"CanvasBaselineVariationTrendPedestal",grbase[icoll][kassoc][jmes],grv2[icoll][kassoc][jmes]);
	  Printf("Histo subtrcated obtained");
	  cout<<"sub -> "<<subtractedhisto[icoll][kassoc][jmes]->GetBinContent(5)<<endl;	
	  Printf("Histo subtrcated: pointer is working");
	  grbase[icoll][kassoc][jmes]->SetFillStyle(fillColourBaselineStyle);
	  grbase[icoll][kassoc][jmes]->SetFillColor(colourSystemBaseUnc[icoll]);// was kRed-7
	  grbase[icoll][kassoc][jmes]->SetLineColor(colourSystemBaseUnc[icoll]);//was kRed-7

	  if(grv2[icoll][kassoc][jmes]){
	    Printf("SHOULD NOT ENTER HERE");
	    grv2[icoll][kassoc][jmes]->SetFillStyle(1001);
	    grv2[icoll][kassoc][jmes]->SetFillColor(kMagenta);  
	  }
	  suberr[icoll][kassoc][jmes]->SetLineColor(colourSystem[icoll]);
	  suberr[icoll][kassoc][jmes]->SetLineWidth(2);
	}
	
	if(icoll==1) { // p-Pb style
	  histo[icoll][kassoc][jmes]->SetMarkerColor(colourSystem[icoll]); 
	  histo[icoll][kassoc][jmes]->SetLineColor(colourSystem[icoll]); 
	  histo[icoll][kassoc][jmes]->SetMarkerStyle(markerstyle[icoll]); 
	  histo[icoll][kassoc][jmes]->SetMarkerSize(markersize); 

	    subtractedhisto[icoll][kassoc][jmes] = GetPedestalHistoAndSystAndSubtractPedpPb(icoll,kassoc,jmes,histo[icoll][kassoc][jmes],err[icoll][kassoc][jmes],suberr[icoll][kassoc][jmes],"CanvasBaselineVariationTrendPedestal",grbase[icoll][kassoc][jmes],grv2[icoll][kassoc][jmes]);
	    cout<<"sub -> "<<subtractedhisto[icoll][kassoc][jmes]->GetBinContent(5)<<endl;
	    subtractedhisto[icoll][kassoc][jmes]->GetYaxis()->SetRangeUser(-0.8,2*subtractedhisto[icoll][kassoc][jmes]->GetBinContent(subtractedhisto[icoll][kassoc][jmes]->GetMaximumBin()));    
	    grbase[icoll][kassoc][jmes]->SetFillStyle(fillColourBaselineStyle);
	    grbase[icoll][kassoc][jmes]->SetFillColor(colourSystemBaseUnc[icoll]);// was kRed-7
	    grbase[icoll][kassoc][jmes]->SetLineColor(colourSystemBaseUnc[icoll]);//was kRed-7
	    
	    if(grv2[icoll][kassoc][jmes]){
	      grv2[icoll][kassoc][jmes]->SetFillStyle(3002);
	      grv2[icoll][kassoc][jmes]->SetFillColor(kMagenta);  
	    }
	    suberr[icoll][kassoc][jmes]->SetLineColor(colourSystem[icoll]);
	    
	    
	}
	subtractedhisto[icoll][kassoc][jmes]->SetMinimum(minYaxis[kassoc]);
      	subtractedhisto[icoll][kassoc][jmes]->SetMaximum(maxYaxis[kassoc]);


      }      
    }
  }

  Printf("All Histos and Graphs created");
  // UP TO HERE SHOULD BE OK
  
  TCanvas *cFinalPaperStyle;
  if(skip3to5){
    cFinalPaperStyle=new TCanvas("cFinalPaperStyle","cFinalPaperStyle",ratioForCanvasWidth*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
	cFinalPaperStyle->Divide(2,3,0.0,0.0,0);
	Set3x2PadPositions(cFinalPaperStyle);
	cFinalPaperStyle->Modified();
	cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);
  }
  else{
    cFinalPaperStyle=new TCanvas("cFinalPaperStyle","cFinalPaperStyle",1000,1000);
    SetPadStyle(cFinalPaperStyle);
    cFinalPaperStyle->Divide(3,3,0,0);
    Set3x3PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
  }



  for(Int_t iassoc=0;iassoc<nbinAssocpt;iassoc++){
    for(Int_t jDpt=firstDpt;jDpt<nbinDpt;jDpt++){
      TPad *pd=(TPad*)cFinalPaperStyle->cd((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      SetPadStyle(pd);
      //    gPad->SetLeftMargin(0.3);
      //    gPad->SetRightMargin(0.);
      //    gPad->SetTopMargin(0.3);
      //    gPad->SetBottomMargin(0.3);
      
      Double_t startline,endline,wi;
      wi=subtractedhisto[1][iassoc][jDpt]->GetBinWidth(1);
      startline=subtractedhisto[1][iassoc][jDpt]->GetBinLowEdge(1);
      endline=subtractedhisto[1][iassoc][jDpt]->GetBinLowEdge(subtractedhisto[1][iassoc][jDpt]->GetNbinsX())+wi+0.3;// 0.15
      // cout<<"******************* "<<startline<<"  "<<endline<<endl;
      TLine* line=new TLine(startline, 0, endline, 0);
      line->SetLineStyle(2);
      
      TH1D* h=new TH1D(*subtractedhisto[1][iassoc][jDpt]);
      h->Reset();
      h->GetXaxis()->SetLimits(startline,endline);
      h->SetLineColor(0);
      h->GetYaxis()->SetNdivisions(505,kTRUE);
      h->GetXaxis()->CenterTitle();
      h->GetYaxis()->CenterTitle();
      /*
	// older settings, still valid
	h->GetYaxis()->SetTitleFont(42);
	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelFont(42);
	h->GetYaxis()->SetTitleSize(0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
	h->GetXaxis()->SetTitleSize(0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
	h->GetYaxis()->SetLabelSize(0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
	h->GetXaxis()->SetLabelSize(0.07/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);

*/

      h->GetYaxis()->SetTitleFont(43);
      h->GetXaxis()->SetTitleFont(43);
      h->GetYaxis()->SetLabelFont(43);
      h->GetXaxis()->SetLabelFont(43);
      h->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
      h->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
      h->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
      h->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);

      h->GetYaxis()->SetTitleOffset(ytitleoffset*innerPadHeight/referencePadHeight*resizeTextFactor);//0.95*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
      h->GetXaxis()->SetTitleOffset(xtitleoffset*innerPadHeight/referencePadHeight*resizeTextFactor);// not number before, no number was optimized before font 43 was introduced
      if(jDpt-firstDpt!=0){
	h->GetYaxis()->SetTitleSize(0);
	h->GetYaxis()->SetLabelSize(0);
      }

      if(iassoc!=nbinAssocpt-1){
	h->GetXaxis()->SetTitleSize(0);
	h->GetXaxis()->SetLabelSize(0);
      }
      
      //      h->GetXaxis()->Delete();
      //      h->GetYaxis()->Delete();
      h->Draw();


      //  box2[4]->Draw("same");
      SetColorGraphs(0,grbase[0][iassoc][jDpt],colourSystemBaseUnc[0]);// redundant...
      SetColorGraphs(1,grbase[1][iassoc][jDpt],colourSystemBaseUnc[1]);// redundant...
      grbase[0][iassoc][jDpt]->SetFillStyle(fillColourBaselineStyle);
      grbase[1][iassoc][jDpt]->SetFillStyle(fillColourBaselineStyle);
      grbase[1][iassoc][jDpt]->Draw("E2");
      grbase[0][iassoc][jDpt]->Draw("E2"); 
      
      line->Draw(); 
      subtractedhisto[1][iassoc][jDpt]->Draw("same");
      suberr[1][iassoc][jDpt]->Draw("E2");
      Printf("Setting scale unc text pos and size");
      SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],ltscale[1][iassoc][jDpt]);
      
      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1){
	ltscale[0][iassoc][jDpt]->SetY(0.12/gPad->GetHNDC()+gPad->GetBottomMargin());// older settings 0.135/gPad->GetHNDC()+gPad->GetBottomMargin());
	ltscale[1][iassoc][jDpt]->SetY(0.09/gPad->GetHNDC()+gPad->GetBottomMargin());// older settings 0.105/gPad->GetHNDC()+gPad->GetBottomMargin());
	ltscale[0][iassoc][jDpt]->SetX(0.11/gPad->GetWNDC()+gPad->GetLeftMargin());
	ltscale[1][iassoc][jDpt]->SetX(0.11/gPad->GetWNDC()+gPad->GetLeftMargin());
      }
      // not in older settings:
      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==2){
	ltscale[0][iassoc][jDpt]->SetY(0.15/gPad->GetHNDC()+gPad->GetBottomMargin());
	ltscale[1][iassoc][jDpt]->SetY(0.12/gPad->GetHNDC()+gPad->GetBottomMargin());	
      }

      histopPbSuperimp[iassoc][jDpt] = subtractedhisto[1][iassoc][jDpt]->Clone();
      histopPbSuperimp[iassoc][jDpt]->SetMarkerStyle(25);
      histopPbSuperimp[iassoc][jDpt]->SetMarkerColor(kRed+1);
      histopPbSuperimp[iassoc][jDpt]->Draw("same");
      //      ltscale[0][iassoc][jDpt]->SetTextSize(0.05/(gPad->GetHNDC())*scaleHeightPads);


      //      ltscale[1][iassoc][jDpt]->SetTextSize(0.05/(gPad->GetHNDC())*scaleHeightPads);


      ltscale[1][iassoc][jDpt]->Draw();//commented
      suberr[0][iassoc][jDpt]->Draw("E2");
      subtractedhisto[0][iassoc][jDpt]->Draw("same");
      ltscale[0][iassoc][jDpt]->Draw();//commented
      Printf("Getting data legend");
      TLegend *legend=GetLegendDataPoints(histo[0][iassoc][firstDpt],histo[1][iassoc][firstDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TLegend *legendSuperimp=GetLegendDataPointsFake(histo[0][iassoc][firstDpt],histopPbSuperimp[iassoc][firstDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      Printf("Getting baseline unc legend");
      TLegend *legend2;
      if(style==0)legend2=GetLegendBaselines(grbase[0][iassoc][firstDpt],grbase[1][iassoc][firstDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,0);
      else if(style==1)legend2=GetLegendBaselines(grbase[0][iassoc][firstDpt],grbase[1][iassoc][firstDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,3); // default option was 1 in older settings; 2 for draft 1 version;
      else legend2=GetLegendBaselines(grbase[0][iassoc][firstDpt],grbase[1][iassoc][firstDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,1);
      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1){legend->Draw(); legendSuperimp->Draw();}
      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==2)legend2->Draw();// was in pad 3 before
      Printf("Getting average pave text");
      TPaveText *pvAverage=GetAveragepavetext((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      Printf("Getting ALICE pave text");
      TPaveText *alice=GetALICEpavetext((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);  
      alice->SetX1(0.8);
      TPaveText *pvKineInfo;

      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1){
	Printf("Getting Pave Kine infos in bin main, aligning left");
	pvKineInfo=GetPaveKineInfo(strPtMesonText[jDpt],strPtAssocText[iassoc],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,22);
	//was: pvKineInfo=GetPaveKineInfo(strPtMesonText[jDpt],0x0,(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,2);		
	//	pvKineInfo->SetY1(0.225/gPad->GetHNDC()+gPad->GetBottomMargin());
	//	pvKineInfo->SetY2(0.24/gPad->GetHNDC()+gPad->GetBottomMargin());

	pvKineInfo->SetY1(0.15/gPad->GetHNDC()+gPad->GetBottomMargin());
	pvKineInfo->SetY2(0.175/gPad->GetHNDC()+gPad->GetBottomMargin());
	pvKineInfo->SetX1(0.008/gPad->GetWNDC()+gPad->GetLeftMargin());
	pvKineInfo->SetX2(0.23/gPad->GetWNDC()+gPad->GetLeftMargin());
	pvKineInfo->SetTextAlign(12);	
	pvKineInfo->Draw("same");
	//	Printf("Getting Pave Kine infos in bin main, aligning right");
	//	pvKineInfo=GetPaveKineInfo(0x0,strPtAssocText[iassoc],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,2);// was option 11
	//	pvKineInfo->SetY1(0.225/gPad->GetHNDC()+gPad->GetBottomMargin());
	//	pvKineInfo->SetY2(0.24/gPad->GetHNDC()+gPad->GetBottomMargin());
	//	pvKineInfo->Draw("same");	
	Printf("Getting Pave Kine infos (DeltaEta) in bin main, aligning left");
	pvKineInfo=GetPaveKineInfo(0x0,0x0,(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,1);// just a trick to get DeltaEta
	Printf("Delta Eta pave obtained");
	//	pvKineInfo->SetY1(0.20/gPad->GetHNDC()+gPad->GetBottomMargin());
	//	pvKineInfo->SetY2(0.215/gPad->GetHNDC()+gPad->GetBottomMargin());
	pvKineInfo->SetY1(0.127/gPad->GetHNDC()+gPad->GetBottomMargin());
	pvKineInfo->SetY2(0.142/gPad->GetHNDC()+gPad->GetBottomMargin());
	pvKineInfo->SetX1(0.008/gPad->GetWNDC()+gPad->GetLeftMargin());
	pvKineInfo->SetX2(0.23/gPad->GetWNDC()+gPad->GetLeftMargin());
	pvKineInfo->SetTextAlign(12);
	pvKineInfo->Draw("same");
      }
      else{
	Printf("Getting Pave Kine infos");
	if(style==0){
	  pvKineInfo=GetPaveKineInfo(strPtMesonText[jDpt],strPtAssocText[iassoc],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,10);
	}
	else if(style==1){
	  pvKineInfo=GetPaveKineInfo(strPtMesonText[jDpt],strPtAssocText[iassoc],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,22);
	  if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==2){
	    pvKineInfo->SetY2(0.255/gPad->GetHNDC()+gPad->GetBottomMargin());
	    pvKineInfo->SetY1(0.275/gPad->GetHNDC()+gPad->GetBottomMargin());
	  }
	}
	else {// very old default
	  pvKineInfo=GetPaveKineInfo(strPtMesonText[jDpt],strPtAssocText[iassoc],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1,0);
	}
	pvKineInfo->SetX1(0.2);
	pvKineInfo->SetX2(0.8);
        
	pvKineInfo->Draw("same");
      }
      if(jDpt==0)  {
        pvKineInfo->SetX1(0.4);
	pvKineInfo->SetX2(0.85);
      }

      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1) alice->Draw("same");
      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1) pvAverage->Draw("same");

      Printf("Bin D :%d, bin assoc %d done",jDpt,iassoc);
    }
  }
  

  for(Int_t j=6;j>=1;j--){
    Printf("IN THE LOOP WHERE PADS ARE DRAWN, JUST BEFORE SAVING");
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }
  cFinalPaperStyle->SaveAs(Form("plotComparison_%sAverage_pp_pPb_UniqueCanvas_Style%d.root",avType.Data(),style));
  cFinalPaperStyle->SaveAs(Form("plotComparison_%sAverage_pp_pPb_UniqueCanvas_Style%d.png",avType.Data(),style));
  cFinalPaperStyle->SaveAs(Form("plotComparison_%sAverage_pp_pPb_UniqueCanvas_Style%d.pdf",avType.Data(),style));
  cFinalPaperStyle->SaveAs(Form("plotComparison_%sAverage_pp_pPb_UniqueCanvas_Style%d.eps",avType.Data(),style));

}


void Set3x3PadPositions(TCanvas* c){
    
  TPad * pd1 =(TPad*) c->GetPad(1);
  TPad * pd2 =(TPad*) c->GetPad(2);
  TPad * pd3 =(TPad*) c->GetPad(3);
  TPad * pd4 =(TPad*) c->GetPad(4);
  TPad * pd5 =(TPad*) c->GetPad(5);
  TPad * pd6 =(TPad*) c->GetPad(6);
  TPad * pd7 =(TPad*) c->GetPad(7);
  TPad * pd8 =(TPad*) c->GetPad(8);
  TPad * pd9 =(TPad*) c->GetPad(9);

  SetPadStyle(pd1);

  Double_t xl,xu,yl,yu;
  Double_t marginLeft=0.1;
  Double_t marginRight=0.02;
  Double_t marginTop=0.04;
  Double_t marginBottom=0.08;
  innerPadWidth=(1-marginLeft-marginRight)/3.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/3.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!
  Printf("innerPadHeight: %f",innerPadHeight);
  Printf("innerPadWidth: %f",innerPadWidth);
  Double_t marginLeftForXAxis=0.02;
  Double_t marginBottomForYAxis=0.0;

 // Bottom row

    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 7 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd7->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd7->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd7->SetRightMargin(0);
    pd7->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd7->SetTopMargin(0.);
    pd7->Modified();
    pd7->SetFillStyle(0);

    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 8 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    Float_t xLowPad8=innerPadWidth+marginLeft-0.006;
    Float_t marginLeftPad8=0.015;
    pd8->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd8->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    //pd8->SetLeftMargin(0.004/(1.-innerPadWidth-marginLeft+0.004));
    pd8->SetRightMargin(0.);
    pd8->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd8->SetTopMargin(0);
    pd8->SetFillStyle(0);

    pd8->Modified();
    pd8->Update();

    pd9->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 9 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    Float_t xLowPad9=2*innerPadWidth+marginLeft-0.004;
    Float_t xHighPad9=3*innerPadWidth+marginLeft+0.02;
    Float_t marginLeftPad9=0.008/(1.-innerPadWidth-marginLeft+0.008);
  
    pd9->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd9->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd9->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd9->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd9->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd9->SetTopMargin(0.);
    pd9->SetFillStyle(0);

    pd9->Modified();
    pd9->Update();

    // Middle Row
    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd4->SetRightMargin(0.);
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd4->SetTopMargin(0.);
    pd4->SetFillStyle(0);

    pd4->Modified();
    pd4->Update();

    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd5->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,2*innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd5->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd5->SetRightMargin(0.);
    pd5->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd5->SetTopMargin(0.);
    pd5->SetFillStyle(0);

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
    pd1->SetFillStyle(0);
    pd1->Modified();
    pd1->Update();
    
    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,2*innerPadWidth+marginLeft,1);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd2->SetFillStyle(0);
    pd2->SetRightMargin(0);
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd2->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd3->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd3->SetTopMargin(marginTop/(innerPadHeight+marginTop));
    pd3->SetFillStyle(0);
 
    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();

}

