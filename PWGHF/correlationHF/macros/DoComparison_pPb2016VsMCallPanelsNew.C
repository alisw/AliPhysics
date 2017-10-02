TString inputdatadirectory = "/Users/elenabruna/Documents/ALICE/HFCJ/Dh/plots_June18/ReflectedPlots/StdRebin/AllPlots/Averages";//***elena****
TString inputtemplatedirecotry="/Users/elenabruna/Documents/ALICE/HFCJ/Dh/plots_June18/Templates_pp_12May15";//***elena****
TString baselinedirectory="/Users/elenabruna/Documents/ALICE/HFCJ/Dh/plots_June18/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";//***elena****
// TString inputdatadirectory = "./";
// TString inputtemplatedirecotry="./";
// TString baselinedirectory="./";
TString avType="Weighted";
Bool_t reflTempl=kTRUE;//***elena****
Bool_t skip3to5=kFALSE;//***elena****
const Int_t nbinDpt=4;
const Int_t nbinAssocpt=6;
const Int_t nSets=7;
Int_t firstDpt=0;
Bool_t isReflectedData=kTRUE;
Double_t canvasheight=900;
Double_t ratioForCanvasWidth=3./3.*1.25;
Double_t resizeTextFactor=1.1;//canvasheight/1800.; // size was tuned for canvasheight =1800. 
Double_t innerPadHeight;// not touch, set internally
Double_t innerPadWidth;// not touch, set internally
Double_t referencePadHeight=0.44; 
TString strsyst="pPb";
TString sets[nSets]={"pPb","Perugia0","Perugia2010","Perugia2011","PYTHIA8","POWHEG","EPOS3"};
TString setsBoost[nSets]={"pPb","Perugia0wBoost","Perugia2010wBoost","Perugia2011wBoost","PYTHIA8wBoost","POWHEG","EPOS3"};
Bool_t includeset[nSets]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kFALSE};
void SetIncludePerugia0(Bool_t incl){
  includeset[1]=incl;
}
void SetIncludePerugia2010(Bool_t incl){
  includeset[2]=incl;
}
void SetIncludePerugia2011(Bool_t incl){
  includeset[3]=incl;
}
void SetIncludePYTHIA8(Bool_t incl){
  includeset[4]=incl;
}
void SetIncludePOWHEG(Bool_t incl){
  includeset[5]=incl;
}
void SetIncludeEPOS(Bool_t incl){
  includeset[6]=incl;
}
void SetIncludeAllMCmodels(Bool_t incl=kTRUE){
  for(Int_t j=1;j<nSets;j++){
    includeset[j]=incl;
  }
}
const Int_t nmodels=7;// does not matter that they are 6 maximum now; this is used only to define the 2 arrays modelColors and modelMarkerStyle
Bool_t splitLegendMC=kFALSE;
Color_t modelColors[nmodels]={kRed+2,kCyan,kGreen+2,kMagenta+1,kBlue,kOrange+1,kViolet};
Int_t modelMarkerStyle[nmodels]={kOpenSquare,kOpenCircle,kOpenDiamond,3,28,26,33};
TString pthadron[nbinAssocpt]={"0.3to99.0","0.3to1.0","1.0to99.0","1.0to2.0","2.0to3.0","3.0to99.0"};
TString strmesonpt[nbinDpt]={"3to5","5to8","8to16","16to24"};
TString strmesonMCpt[nbinDpt]={"3To5","5To8","8To16","16To24"};
Double_t minYaxis[nbinAssocpt]={-0.7,-0.7,-0.45,-0.39,-0.45,-0.23}; // or -0.6 for all  
Double_t maxYaxis[nbinAssocpt]={5.0,1.9,3.2,1.5,1.1,1.1};// or 2.9 for the first and 1.9 for the last
//Double_t maxYaxis[nbinAssocpt]={5.7,2.9,3.2};// or 2.9 for the first and 1.9 for the last  //used so far
//Double_t maxYaxis[nbinAssocpt]={4.8,2.9,3.2};// or 2.9 for the first and 1.9 for the last
//Double_t maxYaxis[nbinAssocpt]={3.8,1.9,2.2};// or 2.9 for the first and 1.9 for the last
Double_t mesonptcenter[nbinDpt]={4.,6.5,12.,20.};// needed for getting the baseline
TString strPtAssocText[nbinAssocpt]={"#it{p}_{T}^{assoc} > 0.3 GeV/#it{c}, |#Delta#eta| < 1.0","0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0","#it{p}_{T}^{assoc} > 1 GeV/#it{c}, |#Delta#eta| < 1.0","1 < #it{p}_{T}^{assoc} < 2 GeV/#it{c}, |#Delta#eta| < 1.0","2 < #it{p}_{T}^{assoc} < 3 GeV/#it{c}, |#Delta#eta| < 1.0","#it{p}_{T}^{assoc} > 3 GeV/#it{c}, |#Delta#eta| < 1.0"};
TString strPtRangeText[nbinDpt*nbinAssocpt]={ "3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} > 0.3 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} > 0.3 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} > 0.3 GeV/#it{c}", "16 < #it{p}_{T}^{D} < 24 GeV/#it{c}, #it{p}_{T}^{assoc} > 0.3 GeV/#it{c}", "3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}", "16 < #it{p}_{T}^{D} < 24 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} < 1 GeV/#it{c}","3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} > 1 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} > 1 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} > 1 GeV/#it{c}", "16 < #it{p}_{T}^{D} < 24 GeV/#it{c}, #it{p}_{T}^{assoc} > 1 GeV/#it{c}", "3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, 1 < #it{p}_{T}^{assoc} < 2 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, 1 < #it{p}_{T}^{assoc} < 2 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, 1 < #it{p}_{T}^{assoc} < 2 GeV/#it{c}", "16 < #it{p}_{T}^{D} < 24 GeV/#it{c}, 1 < #it{p}_{T}^{assoc} < 2 GeV/#it{c}", "3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, 2 < #it{p}_{T}^{assoc} < 3 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, 2 < #it{p}_{T}^{assoc} < 3 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, 2 < #it{p}_{T}^{assoc} < 3 GeV/#it{c}", "16 < #it{p}_{T}^{D} < 24 GeV/#it{c}, 2 < #it{p}_{T}^{assoc} < 3 GeV/#it{c}","3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} > 3 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} > 3 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} > 3 GeV/#it{c}", "16 < #it{p}_{T}^{D} < 24 GeV/#it{c}, #it{p}_{T}^{assoc} > 3 GeV/#it{c}"};
TString strPtMesonText[nbinDpt]={ "3 < #it{p}_{T}^{D} < 5 GeV/#it{c}", "5 < #it{p}_{T}^{D} < 8 GeV/#it{c}", "8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, |#it{y}^{D}| < 0.5", "16 < #it{p}_{T}^{D} < 24 GeV/#it{c}, |#it{y}^{D}| < 0.5"};

TString strYText= "|#it{y}^{D}_{cms}| < 0.5, |#Delta#eta| < 1";

TString filenames[nSets][nbinAssocpt][nbinDpt];// [coll syst][ptassoc][ptmes]
TString pedestalfilenames[nSets][nbinAssocpt];// [coll syst][ptassoc]
Double_t scaleHeightPads=1;// do not touch this, it is regulated automatically in the macro
Double_t scaleWidthPads=1;// do not touch this, it is regulated automatically in the macro


TH1D ****histo;
TGraphAsymmErrors ****err;
TLatex ****ltscale;
TH1D ****subtractedhisto;
TGraphAsymmErrors ****suberr;
TGraphAsymmErrors ****grbase;
TGraphAsymmErrors ****grv2; 

TString fitplotmacrodir=gSystem->ExpandPathName("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/");

TString strSystemFDtempl="none";
void SetFDtemplateSystemString(TString str){
  strSystemFDtempl=str;
}

void SetFitPlotMacroPath(TString strdir){
  fitplotmacrodir=strdir;
}

Bool_t isReflectedData=kTRUE;
void SetSkip3to5pPb(Bool_t skip){
  skip3to5=skip;
}
void SetIsDataReflected(Bool_t isrefl){
  isReflectedData=isrefl;
}
void SetBaselineDirectory(TString dirbase){
  baselinedirectory=dirbase;
}
void SetReflectTemplate(Bool_t doreflTempl){
  reflTempl=doreflTempl;
}
void SetInputDataDirectory(TString inputdir){
  inputdatadirectory=inputdir;
}
void SetInputTemplateDirectory(TString inputdir){
  inputtemplatedirecotry=inputdir;
}
void SetAverageMode(Int_t avmode){
  if(avmode==0)avType="Weighted";
  else if(avmode==1)avType="Arithmetic";
  else Printf("DO COMPARISON WITH MC: WRONGE AVERAGE METHOD SET");
}
void SetSplitMClegendInTwoPanels(Bool_t split=kTRUE){splitLegendMC=split;}

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
  /// p->SetFrameBorderMode(0);
  p->SetFillStyle(0);
  p->SetFrameFillStyle(4000);
  p->SetFrameBorderMode(0);

  //p->SetBorderMode(0);
  //  p->SetBorderSize(0);
  //  p->SetBottomMargin(0);
  //  p->SetTopMargin(0);
  //  p->Range(-2.6,-1.5,5.4,4.4);
  //  p->SetLeftMargin(leftMarginCanvas);
  //  p->SetRightMargin(rightMarginCanvas);
  //  p->SetBottomMargin(bottomMarginCanvas);
  //  p->SetTopMargin(topMarginCanvas);
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

void Init(){

  
  histo=new TH1D***[nSets];
  err=new TGraphAsymmErrors***[nSets]; 
  ltscale = new TLatex ***[nSets]; 
  subtractedhisto=new TH1D***[nSets];
  suberr=new TGraphAsymmErrors***[nSets]; 
  grbase=new TGraphAsymmErrors***[nSets]; 
  grv2=new TGraphAsymmErrors***[nSets]; 
  for(Int_t iset=0;iset<nSets;iset++){
    histo[iset]=new TH1D**[nbinAssocpt];
    err[iset]= new TGraphAsymmErrors**[nbinAssocpt];
    ltscale[iset] = new TLatex **[nbinAssocpt];
    subtractedhisto[iset]=new TH1D**[nbinAssocpt];
    suberr[iset]=new TGraphAsymmErrors**[nbinAssocpt]; 
    grbase[iset]=new TGraphAsymmErrors**[nbinAssocpt]; 
    grv2[iset]=new TGraphAsymmErrors**[nbinAssocpt]; 
    for(Int_t kassoc=0;kassoc<nbinAssocpt;kassoc++){
      histo[iset][kassoc]=new TH1D*[nbinDpt];
      subtractedhisto[iset][kassoc]=new TH1D*[nbinDpt];
      err[iset][kassoc]= new TGraphAsymmErrors*[nbinDpt]; 
      ltscale[iset][kassoc] = new TLatex*[nbinDpt]; 
      suberr[iset][kassoc]=new TGraphAsymmErrors*[nbinDpt]; 
      grbase[iset][kassoc]=new TGraphAsymmErrors*[nbinDpt]; 
      grv2[iset][kassoc]=new TGraphAsymmErrors*[nbinDpt]; 
      for(Int_t jmes=0;jmes<nbinDpt;jmes++){
      histo[iset][kassoc][jmes]=0x0;
      subtractedhisto[iset][kassoc][jmes]=0x0;
      err[iset][kassoc][jmes]=0x0; 
      ltscale[iset][kassoc][jmes] =0x0; 
      suberr[iset][kassoc][jmes]=0x0; 
      grbase[iset][kassoc][jmes]=0x0; 
      grv2[iset][kassoc][jmes]=0x0; 
      }
    }
  }
}


//_______________________________________________________________________
void SaveCanvas(TCanvas * c, TString directory, TString name){
  //    
    
  if(directory != ""){outputDir += directory;
    TString exec = "mkdir -p ";
    exec += outputDir;
    cout << exec << endl;
    gSystem->Exec(exec.Data());
  }
    
    
  TString plotsout = c->GetName();"Canvas_pT_05_";
  plotsout += name;
    
  c->SaveAs(Form("%s/%s.root",outputDir.Data(),plotsout.Data()));
  c->SaveAs(Form("%s/%s.eps",outputDir.Data(),plotsout.Data()));
  c->SaveAs(Form("%s/%s.png",outputDir.Data(),plotsout.Data()));
}

//_______________________________________________________________________
TH1D * GetHistoAndSyst(TString path, Int_t iset,TString hname, TString hnamesyst,TGraphAsymmErrors *&gr2, TLatex *&tUncertainty){
    
    
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

  if(iset==0){//DATA
     if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.55,0.46,Form("#bf{%.0f#% scale uncertainty pp}",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.65,0.6,Form("{}^{#plus%.0f%s}_{#minus%.0f%s} scale uncertainty","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
   
  }

  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.025);
  tUncertainty->SetTextFont(62);
    
  return hFDsub;
    
}
//_______________________________________________________________________
TH1D * GetHisto(TString path, TString hname){
    
  Printf("Opening file: %s",path.Data());
  TFile *file=TFile::Open(path.Data(),"READ");

  TH1D * histo=0x0;
  histo=(TH1D*)file->Get(hname.Data());
  if(!histo){
    TCanvas * c = (TCanvas*)file->Get("cDeltaPhi");
    c->cd();
    histo = (TH1D*)c->FindObject(hname.Data());
  }
  return histo;
    
}
//_______________________________________________________________________
TH1D * GetPedestalHistoAndSystAndSubtractPed(Int_t binSystem, Int_t binAssoc,Int_t binMeson,TH1D *histo, TGraphAsymmErrors* gr,TGraphAsymmErrors *&grout, TString canvasname,TGraphAsymmErrors *&grbaseOut,TGraphAsymmErrors *&grv2Out){

  Double_t value = 0, pedestal=0;
  TGraphAsymmErrors *grBaseHelp,*grV2,*grV2Out;
  grbaseOut=new TGraphAsymmErrors();
  grv2Out=new TGraphAsymmErrors();
  grbaseOut->SetName(Form("grbaselineUncFull_%s_%s_%s",sets[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data()));
  grv2Out->SetName(Form("grbaselineUncFull_%s_%s_%s",sets[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data()));

  Double_t xuncFull,errxuncFull;
  Double_t xuncv2,errxuncv2;
  Int_t bin,bingr;
  TString path = pedestalfilenames[binSystem][binAssoc];


  if(binSystem==0){//pPb
    grV2=0x0;
    cout << "pp -->  Reading File from path: " << path << endl;
    
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    grBaseHelp=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fBaselineVariationSystematicsPedestal");    
    grV2=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fv2SystematicsPedestal");    

    if(isReflectedData){
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
    printf("PEDESTAL ERROR: Point has x %f, y %f, errlowy %f, errlowx %f\n",x,y,erryl,erryh);

    grbaseOut->SetPoint(0,xuncFull,0);
    grbaseOut->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    if(grV2){
      grV2->GetPoint(bingr,x,y);
      erryl=grV2->GetErrorYlow(bingr);
      erryh=grV2->GetErrorYhigh(bingr);	
      grv2Out->SetPoint(0,xuncFull,0);
      grv2Out->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
      printf("PEDESTAL ERROR FROM v2: Point has x %f, y %f, errlowy %f, errlowx %f\n",x,y,erryl,erryh);
    }
  }

  cout<<"Baseline being subtracted for:"<<sets[binSystem].Data()<<strmesonpt[binMeson].Data()<<pthadron[binAssoc].Data()<<endl;
  
  grout=(TGraphAsymmErrors*)gr->Clone(Form("grSub_%s_%s_%s",sets[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data()));
  
  TString nameoutput = histo->GetName();
  nameoutput += "_subtr_";
  nameoutput += "pedestal";
  nameoutput += Form("_%s_%s_%s",sets[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data());
   
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
  outputhisto->GetYaxis()->CenterTitle();
  outputhisto->SetMarkerColor(kBlack);
  outputhisto->SetLineColor(kBlack);
  outputhisto->SetMarkerStyle(20);
  
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
//_______________________________________________________________________
TH1D * GetPedestalHistoAndSystAndSubtractPedMC(Int_t binSystem, Int_t binAssoc,Int_t binMeson,TH1D *histo, TString canvasname){

  Double_t value = 0, pedestal=0;
 
  Double_t xuncFull,errxuncFull;
  Double_t xuncv2,errxuncv2;
  Int_t bin,bingr;
  TString path = pedestalfilenames[binSystem][binAssoc];


  cout << "pPb -->  Reading File from path: " << path << endl;
    
  TFile * file = TFile::Open(path.Data(),"READ");
  TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
  TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
  

  if(isReflectedData){
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
   
    pedestal=h->GetBinContent(bin);
    Double_t x,y,erryl,erryh;
    Printf("histo: x=%f, graph: %f",h->GetBinCenter(bin),x);
   

    cout<<"pedestal = "<<pedestal<<"   Baseline being subtracted for:"<<sets[binSystem].Data()<<strmesonpt[binMeson].Data()<<pthadron[binAssoc].Data()<<endl;
  
  
  TString nameoutput = histo->GetName();
  nameoutput += "_subtr_";
  nameoutput += "pedestal";
  nameoutput += Form("_%s_%s_%s",sets[binSystem].Data(),strmesonpt[binMeson].Data(),pthadron[binAssoc].Data());
   
  TH1D * outputhisto = (TH1D*)histo->Clone(nameoutput.Data());
  outputhisto->Reset();
  outputhisto->SetStats(kFALSE);
  outputhisto->GetYaxis()->CenterTitle();
    
  cout<<"*******"<<endl;
  for(Int_t iBin = 1; iBin <= histo->GetNbinsX();iBin++){
        
        
    value = histo->GetBinContent(iBin);
    value -= pedestal;
    cout<<iBin<<"  "<<value<<endl;
    outputhisto->SetBinContent(iBin,value);
      
    outputhisto->SetBinError(iBin,histo->GetBinError(iBin));

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

void DoComparison_pPb2016VsMCallPanels(){

  Init();
  LoadFileNamesAll();
  
  if(!includeset[0]){
    Printf("The pPb dataset is expected to be included! Cannot proceed"); return;
  }
  
  TH1D *h;
  // for(Int_t iset=0;iset<1;iset++){
  for(Int_t iset=0;iset<nSets;iset++){
    if(!includeset[iset])continue;
    for(Int_t kassoc=0;kassoc<nbinAssocpt;kassoc++){
      for(Int_t jmes=firstDpt;jmes<nbinDpt;jmes++){

	if(iset==0){
	  histo[iset][kassoc][jmes]=GetHistoAndSyst(filenames[iset][kassoc][jmes],iset,"fhDaverage","AverageSystematicUncertainty",err[iset][kassoc][jmes],ltscale[iset][kassoc][jmes]); 
	
	  
	// pp style
	  histo[iset][kassoc][jmes]->SetMarkerColor(4); 
	  histo[iset][kassoc][jmes]->SetMarkerStyle(21); 
	  histo[iset][kassoc][jmes]->SetLineColor(4);  
	  
	  subtractedhisto[iset][kassoc][jmes] = GetPedestalHistoAndSystAndSubtractPed(iset,kassoc,jmes,histo[iset][kassoc][jmes],err[iset][kassoc][jmes],suberr[iset][kassoc][jmes],"CanvasBaselineVariationTrendPedestal",grbase[iset][kassoc][jmes],grv2[iset][kassoc][jmes]);
	  Printf("Histo subtrcated obtained");
	  cout<<"sub -> "<<subtractedhisto[iset][kassoc][jmes]->GetBinContent(5)<<endl;	
	  Printf("Histo subtrcated: pointer is working");
	  if(grv2[iset][kassoc][jmes]){
	    Printf("SHOULD NOT ENTER HERE");
	    grv2[iset][kassoc][jmes]->SetFillStyle(3002);
	    grv2[iset][kassoc][jmes]->SetFillColor(kMagenta);  
	  }
	  suberr[iset][kassoc][jmes]->SetLineColor(kBlack);
	}
	

	else{
	  h=GetHisto(filenames[iset][kassoc][jmes],"hCorrDeltaPhi");
	  if(reflTempl)histo[iset][kassoc][jmes]=AliHFCorrelationUtils::ReflectHisto(h,0.5);
	  else histo[iset][kassoc][jmes]=h;

	    histo[iset][kassoc][jmes]->SetMarkerColor(modelColors[iset-1]);// -1 because first iset is data
	    histo[iset][kassoc][jmes]->SetLineColor(modelColors[iset-1]); 
	    histo[iset][kassoc][jmes]->SetLineWidth(2); 
	    histo[iset][kassoc][jmes]->SetMarkerStyle(kDot); 
	    //  histo[iset][kassoc][jmes]->SetMarkerStyle(kOpenSquare); 

	   subtractedhisto[iset][kassoc][jmes] = GetPedestalHistoAndSystAndSubtractPedMC(iset,kassoc,jmes,histo[iset][kassoc][jmes],"CanvasBaselineVariationTrendPedestal");
	   
	}

	
	subtractedhisto[iset][kassoc][jmes]->SetMinimum(minYaxis[kassoc]);
      	subtractedhisto[iset][kassoc][jmes]->SetMaximum(maxYaxis[kassoc]);


      }      
    }
  }

  Printf("All Histos and Graphs created");
  // UP TO HERE SHOULD BE OK
  TCanvas *cFinalPaperStyle;
  TCanvas *cFinalPaperStyle2;

  cFinalPaperStyle=new TCanvas("cFinalPaperStyle","cFinalPaperStyle",1200.,900.);
  cFinalPaperStyle->Divide(4,3,0.0,0.0,0);
  Set4x6PadPositions(cFinalPaperStyle);
  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();

  cFinalPaperStyle2=new TCanvas("cFinalPaperStyle2","cFinalPaperStyle2",1200.,900.);
  cFinalPaperStyle2->Divide(4,3,0.0,0.0,0);
  Set4x6PadPositions(cFinalPaperStyle2);
  cFinalPaperStyle2->Modified();
  cFinalPaperStyle2->Update();
  
  for(Int_t iassoc=0;iassoc<nbinAssocpt;iassoc++){

    if(iassoc>=3) cFinalPaperStyle2->cd();
    for(Int_t jDpt=firstDpt;jDpt<nbinDpt;jDpt++){
      cout<<"#################***************** "<<iassoc<<"  "<<jDpt<<endl;
      TPad *pd;
      if(iassoc<3) pd=(TPad*)cFinalPaperStyle->cd((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      else pd=(TPad*)cFinalPaperStyle2->cd((nbinDpt-firstDpt)*(iassoc-3)+jDpt-firstDpt+1);
      SetPadStyle(pd);
      pd->Modified();
      pd->Update();
      //    gPad->SetLeftMargin(0.3);
      //    gPad->SetRightMargin(0.);
      //    gPad->SetTopMargin(0.3);
      //    gPad->SetBottomMargin(0.3);
      
      Double_t startline,endline,wi;
      wi=subtractedhisto[0][iassoc][jDpt]->GetBinWidth(1);
      startline=subtractedhisto[0][iassoc][jDpt]->GetBinLowEdge(1);
      // endline=subtractedhisto[0][iassoc][jDpt]->GetBinLowEdge(subtractedhisto[0][iassoc][jDpt]->GetNbinsX())+wi+0.4;
      endline=subtractedhisto[0][iassoc][jDpt]->GetBinLowEdge(subtractedhisto[0][iassoc][jDpt]->GetNbinsX())+wi+0.15;
      // cout<<"******************* "<<startline<<"  "<<endline<<endl;
      TLine* line=new TLine(startline, 0, endline, 0);
      line->SetLineStyle(2);
      
      TH1D* h=new TH1D(*subtractedhisto[0][iassoc][jDpt]);
      h->Reset();
      h->GetXaxis()->SetLimits(startline,endline);
      h->SetLineColor(0);

      h->GetYaxis()->SetTitleFont(43);
      h->GetXaxis()->SetTitleFont(43);
      h->GetYaxis()->SetLabelFont(43);
      h->GetXaxis()->SetLabelFont(43);
      h->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.07/(gPad->GetHNDC())*scaleHeightPads);
      h->GetXaxis()->SetTitleSize(30*innerPadHeight/referencePadHeight*resizeTextFactor);//;0.07/(gPad->GetWNDC())*scaleHeightPads);
      h->GetYaxis()->SetLabelSize(30*innerPadHeight/referencePadHeight*resizeTextFactor);//0.07/(gPad->GetHNDC())*scaleHeightPads);
      h->GetXaxis()->SetLabelSize(30*innerPadHeight/referencePadHeight*resizeTextFactor);//0.07/(gPad->GetWNDC())*scaleHeightPads);
      h->GetYaxis()->SetTitleOffset(6.5*innerPadHeight/referencePadHeight*resizeTextFactor);//8*innerPadHeight/referencePadHeight*resizeTextFactor);//1.3*(gPad->GetHNDC())/scaleHeightPads);
      h->GetXaxis()->CenterTitle();
      h->GetYaxis()->CenterTitle();

      //      h->GetXaxis()->Delete();
      //      h->GetYaxis()->Delete();
      if(jDpt>0){
	     h->GetYaxis()->SetTitleSize(0);
	     h->GetYaxis()->SetLabelSize(0);
      }

      if(iassoc!=2 && iassoc!=5){
	     h->GetXaxis()->SetTitleSize(0);
	     h->GetXaxis()->SetLabelSize(0);
      }

      if((iassoc==2 && jDpt>=0) || (iassoc==5 && jDpt>=0)){
	     //	h->GetXaxis()->SetLabelOffset(-0.005);
	     h->GetXaxis()->SetTitleOffset(3.5*innerPadHeight/referencePadHeight*resizeTextFactor);//8*innerPadHeight/referencePadHeight*resizeTextFactor);

      }

      /*     if(iassoc==2 && jDpt>0==1){
	     h->GetXaxis()->SetLabelOffset(-0.005);
	     h->GetXaxis()->SetTitleOffset(0.8);

      }
      if(iassoc==2 && jDpt>0==2){
	     h->GetXaxis()->SetLabelOffset(-0.004);
	     h->GetXaxis()->SetTitleOffset(0.8);

      }
      */

      h->Draw();


      //  box2[4]->Draw("same");
      // SetColorGraphs(0,grbase[0][iassoc][jDpt]);

      //      grbase[0][iassoc][jDpt]->SetFillStyle(3002);
      grbase[0][iassoc][jDpt]->SetFillColor(kGray+2);
      grbase[0][iassoc][jDpt]->SetLineColor(kGray+2);
      grbase[0][iassoc][jDpt]->Draw("E2");
/* //DISABLED BY FABIO - COMPLETELY NEGLIGIBLE UNCERTAINTIES...
      grv2[0][iassoc][jDpt]->SetFillColor(kGreen+1);
      grv2[0][iassoc][jDpt]->SetLineColor(kGreen+1);
      grv2[0][iassoc][jDpt]->Draw("E2");
*/
      line->Draw(); 
      subtractedhisto[0][iassoc][jDpt]->Draw("same");

      suberr[0][iassoc][jDpt]->SetLineColor(kBlack);
      suberr[0][iassoc][jDpt]->Draw("E2");

      Float_t size=0.058;
      SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.104,22);//0.08
      //      if(iassoc==2)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.11,0.2,0.06);
      if(iassoc==0)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.104,22);//0.087
      if(iassoc==1)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.087
      //      if(iassoc==0 && jDpt==1)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.11,0.11,0.062);
      if(iassoc==2 && jDpt==0)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.07
      if(iassoc==2 && jDpt==1)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.085
      if(iassoc==2 && jDpt==2)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.078
      if(iassoc!=0 && jDpt==3)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.078
      if(iassoc>=4 && jDpt==0)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.07
      if(iassoc>=4 && jDpt==1)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.085
      if(iassoc>=4 && jDpt==2)SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.04,0.2,22);//0.078
      

      //     if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1){
      //ltscale[0][iassoc][jDpt]->SetY(0.11/gPad->GetHNDC()+gPad->GetBottomMargin());
      //}
      ltscale[0][iassoc][jDpt]->Draw();
      for(Int_t kmod=1;kmod<nSets;kmod++){
	if(includeset[kmod])subtractedhisto[kmod][iassoc][jDpt]->Draw("hist same c");
      }
//       subtractedhisto[1][iassoc][jDpt]->Draw("hist same c");//Perugia0
//       subtractedhisto[2][iassoc][jDpt]->Draw("hist same c");//Perugia2010
//       subtractedhisto[3][iassoc][jDpt]->Draw("hist same c");//Perugia2011
//       subtractedhisto[4][iassoc][jDpt]->Draw("hist same c");//PYTHIA8
//       subtractedhisto[5][iassoc][jDpt]->Draw("hist same c");//POWHEG
//       subtractedhisto[6][iassoc][jDpt]->Draw("hist same c");//EPOS3

     // subtractedhisto[1][iassoc][jDpt]->Draw("same");//Perugia0
     //  subtractedhisto[2][iassoc][jDpt]->Draw("same");//Perugia2010
     //  subtractedhisto[3][iassoc][jDpt]->Draw("same");//Perugia2011
     //  subtractedhisto[4][iassoc][jDpt]->Draw("same");//PYTHIA8
     //  subtractedhisto[5][iassoc][jDpt]->Draw("same");//POWHEG
      subtractedhisto[0][iassoc][jDpt]->Draw("same");//pp DATA 
      suberr[0][iassoc][jDpt]->Draw("E2");
    
      TPaveText *pvAverage=GetAveragepavetext((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TPaveText *alice=GetALICEpavetext((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);  
      //TPaveText *pvKineInfo = GetPaveKineInfo(strPtMesonText[jDpt],strPtAssocText[iassoc],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TPaveText *pvKineInfo = GetPaveKineInfo(jDpt,iassoc,(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TPaveText *pvKineInfo2 = GetPaveKineInfo2((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TLegend *legendbase=GetLegendBaselines(grbase[0][iassoc][firstDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      if(jDpt==0 && iassoc==0) legendbase->Draw();
      //   if(jDpt==1 && iassoc==0) legendbase->Draw();
      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1)alice->Draw();
      pvKineInfo->Draw("same");
      if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1)pvAverage->Draw("same");
      //if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1)pvKineInfo2->Draw("same");
      if(jDpt==0 && iassoc==0)pvKineInfo2->Draw("same");

      if(jDpt==0 && iassoc==3) legendbase->Draw();
      //   if(jDpt==1 && iassoc==0) legendbase->Draw();
      if((nbinDpt-firstDpt)*(iassoc-3)+jDpt-firstDpt+1==1)alice->Draw();
      pvKineInfo->Draw("same");
      if((nbinDpt-firstDpt)*(iassoc-3)+jDpt-firstDpt+1==1)pvAverage->Draw("same");
      //if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1)pvKineInfo2->Draw("same");
      if(jDpt==0 && iassoc==3)pvKineInfo2->Draw("same");
      //   if(jDpt==1 && iassoc==0)pvKineInfo2->Draw("same");

      // N.B. GetLegendData does not assume that the pointers to histograms that are passed are not null --> no need to check which model is included 
      if(jDpt==0 && iassoc==0){
	TLegend *legendData=GetLegendData(subtractedhisto[0][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
	legendData->Draw();
	
      }
      if(jDpt==1 && iassoc==0){
	TLegend *legendMC;
	if(splitLegendMC){
	  legendMC=GetLegendMC(subtractedhisto[1][iassoc][jDpt],subtractedhisto[2][iassoc][jDpt],subtractedhisto[3][iassoc][jDpt],0x0,0x0,0x0,(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
	}
	else{
	  legendMC=GetLegendMC(subtractedhisto[1][iassoc][jDpt],subtractedhisto[2][iassoc][jDpt],subtractedhisto[3][iassoc][jDpt],subtractedhisto[4][iassoc][jDpt],subtractedhisto[5][iassoc][jDpt],subtractedhisto[6][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
	}
	legendMC->Draw();
      }
      if(jDpt==2 && iassoc==0 && splitLegendMC){
	
	TLegend *legendMC=GetLegendMC(0x0,0x0,0x0,subtractedhisto[4][iassoc][jDpt],subtractedhisto[5][iassoc][jDpt],subtractedhisto[6][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
	
	legendMC->Draw();
      }

      if(jDpt==0 && iassoc==3){
  TLegend *legendData=GetLegendData(subtractedhisto[0][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
  legendData->Draw();
  
      }
      if(jDpt==1 && iassoc==3){
  TLegend *legendMC;
  if(splitLegendMC){
    legendMC=GetLegendMC(subtractedhisto[1][iassoc][jDpt],subtractedhisto[2][iassoc][jDpt],subtractedhisto[3][iassoc][jDpt],0x0,0x0,0x0,(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
  }
  else{
    legendMC=GetLegendMC(subtractedhisto[1][iassoc][jDpt],subtractedhisto[2][iassoc][jDpt],subtractedhisto[3][iassoc][jDpt],subtractedhisto[4][iassoc][jDpt],subtractedhisto[5][iassoc][jDpt],subtractedhisto[6][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
  }
  legendMC->Draw();
      }
      if(jDpt==2 && iassoc==3 && splitLegendMC){
  
  TLegend *legendMC=GetLegendMC(0x0,0x0,0x0,subtractedhisto[4][iassoc][jDpt],subtractedhisto[5][iassoc][jDpt],subtractedhisto[6][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
  
  legendMC->Draw();
      }

    }
  }


  // for(Int_t j=1;j<=9;j++){
  //   // for(Int_t j=6;j>=1;j--){
  //   TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
  //   pd->Draw();
  // }

  TString nameout="1";
  cFinalPaperStyle->SaveAs(Form("CorrelationppMC4x6_%sNew.pdf",nameout.Data()));
  cFinalPaperStyle->SaveAs(Form("CorrelationppMC4x6_%sNew.eps",nameout.Data()));
  cFinalPaperStyle->SaveAs(Form("CorrelationppMC4x6_%sNew.gif",nameout.Data()));
  cFinalPaperStyle->SaveAs(Form("CorrelationppMC4x6_%sNew.png",nameout.Data()));
  cFinalPaperStyle->SaveAs(Form("CorrelationppMC4x6_%sNew.root",nameout.Data()));

  nameout="2";
  cFinalPaperStyle2->SaveAs(Form("CorrelationppMC4x6_%sNew.pdf",nameout.Data()));
  cFinalPaperStyle2->SaveAs(Form("CorrelationppMC4x6_%sNew.eps",nameout.Data()));
  cFinalPaperStyle2->SaveAs(Form("CorrelationppMC4x6_%sNew.gif",nameout.Data()));
  cFinalPaperStyle2->SaveAs(Form("CorrelationppMC4x6_%sNew.png",nameout.Data()));
  cFinalPaperStyle2->SaveAs(Form("CorrelationppMC4x6_%sNew.root",nameout.Data()));
  
}

void Set4x6PadPositions(TCanvas* c){
    
  TPad * pd1 =(TPad*) c->GetPad(1);
  TPad * pd2 =(TPad*) c->GetPad(2);
  TPad * pd3 =(TPad*) c->GetPad(3);
  TPad * pd4 =(TPad*) c->GetPad(4);
  TPad * pd5 =(TPad*) c->GetPad(5);
  TPad * pd6 =(TPad*) c->GetPad(6);
  TPad * pd7 =(TPad*) c->GetPad(7);
  TPad * pd8 =(TPad*) c->GetPad(8);
  TPad * pd9 =(TPad*) c->GetPad(9);
  TPad * pd10 =(TPad*) c->GetPad(10);
  TPad * pd11 =(TPad*) c->GetPad(11);
  TPad * pd12 =(TPad*) c->GetPad(12);

  SetPadStyle(pd1);

  Double_t xl,xu,yl,yu;
  Double_t marginLeft=0.1;
  Double_t marginRight=0.02;
  Double_t marginTop=0.04;
  Double_t marginBottom=0.08;
  innerPadWidth=(1-marginLeft-marginRight)/4.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/3.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!
  Printf("innerPadHeight: %f",innerPadHeight);
  Printf("innerPadWidth: %f",innerPadWidth);
  Double_t marginLeftForXAxis=0.02;
  Double_t marginBottomForYAxis=0.0;

 // Bottom row

    pd9->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 7 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd9->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd9->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd9->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd9->SetRightMargin(0);
    pd9->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd9->SetTopMargin(0.);
    pd9->Modified();
    pd9->SetFillStyle(0);

    pd10->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 10 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd10->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd10->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd10->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    //pd8->SetLeftMargin(0.004/(1.-innerPadWidth-marginLeft+0.004));
    pd10->SetRightMargin(0.);
    pd10->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd10->SetTopMargin(0);
    pd10->SetFillStyle(0);

    pd10->Modified();
    pd10->Update();

    pd11->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 11 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
  
    pd11->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,3*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd11->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd11->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd11->SetRightMargin(0);
    pd11->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd11->SetTopMargin(0.);
    pd11->SetFillStyle(0);

    pd11->Modified();
    pd11->Update();

    pd12->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 12 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd12->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd12->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd12->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd12->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd12->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd12->SetTopMargin(0.);
    pd12->SetFillStyle(0);

    pd12->Modified();
    pd12->Update();

    // Middle Row
    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd5->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd5->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd5->SetRightMargin(0.);
    pd5->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd5->SetTopMargin(0.);
    pd5->SetFillStyle(0);

    pd5->Modified();
    pd5->Update();

    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd6->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,2*innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd6->SetRightMargin(0.);
    pd6->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd6->SetTopMargin(0.);
    pd6->SetFillStyle(0);

    pd6->Modified();
    pd6->Update();

    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd7->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,3*innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd7->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd7->SetRightMargin(0.);
    pd7->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd7->SetTopMargin(0.);

    pd7->Modified();
    pd7->Update();

    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd8->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1.,2.*innerPadHeight+marginBottom);
    pd8->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd8->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd8->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd8->SetTopMargin(0.);
    
    pd8->Modified();
    pd8->Update();

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
    pd3->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,3*innerPadWidth+marginLeft,1);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd3->SetRightMargin(0);
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd3->SetTopMargin(marginTop/(innerPadHeight+marginTop));
    pd3->SetFillStyle(0);
 
    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    pd4->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd4->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd4->SetTopMargin(marginTop/(innerPadHeight+marginTop));
    pd4->SetFillStyle(0);

    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();

}


//_______________________________________________________________________
void LoadFileNamesAll(){
    
  
    for(Int_t iset=0;iset<nSets;iset++){
      if(!includeset[iset])continue;
      for(Int_t kassoc=0;kassoc<nbinAssocpt;kassoc++){
	for(Int_t jmes=0;jmes<nbinDpt;jmes++){
	  
	  if(iset==0){
	    filenames[iset][kassoc][jmes]=Form("%s/%sAverage%sDzeroDstarDplus%s_assoc%s.root",inputdatadirectory.Data(),avType.Data(),sets[iset].Data(),strmesonpt[jmes].Data(),pthadron[kassoc].Data());//pPb data
	  }
          else if(iset>=1 && iset<=4){ //add "wBoost"
            filenames[iset][kassoc][jmes]=Form("%s/%sCorrelationPlots%sPtDzerofromC%s_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsyst.Data(),setsBoost[iset].Data(),strmesonMCpt[jmes].Data(),pthadron[kassoc].Data());//Pythia templates
          }
	  else if(iset==5){
	    filenames[iset][kassoc][jmes] = Form("%s/%sCorrelationPlots%sPtDzerofromC%s_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),"",sets[iset].Data(),strmesonMCpt[jmes].Data(),pthadron[kassoc].Data());//POWHEG
	  }
	  else{
	    filenames[iset][kassoc][jmes] = Form("%s/%sCorrelationPlots%sPtDzerofromC%s_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsyst.Data(),sets[iset].Data(),strmesonMCpt[jmes].Data(),pthadron[kassoc].Data());//MC
	  }
	  cout<<iset<<"  "<<kassoc<<"  "<<jmes<<endl;
	  cout<<filenames[iset][kassoc][jmes]<<endl;
	  
      }
	
	
	if(iset==0)  pedestalfilenames[iset][kassoc]=Form("%s/Trends_%s/CanvasBaselineVariationTrendPedestal_pthad%s.root",baselinedirectory.Data(),strsyst.Data(),pthadron[kassoc].Data());
	else  pedestalfilenames[iset][kassoc]=Form("%s/FitResults/Trends_%s/%s/CanvasBaselineVariationTrendPedestal_pthad%s.root",inputtemplatedirecotry.Data(),strsyst.Data(),sets[iset].Data(),pthadron[kassoc].Data());
	cout<<"pedestal -> "<<pedestalfilenames[iset][kassoc]<<endl<<endl<<endl;
	
      }
    }
    
}

void SetScaleUncertaintyPositionAndSize(TLatex *tlpp,Float_t xx,Float_t yy,Float_t size){
  // tlpp->SetX(0.6);
  //tlpp->SetY(0.6);

  tlpp->SetX(xx/gPad->GetWNDC()+gPad->GetLeftMargin());
  tlpp->SetY(yy/gPad->GetHNDC()+gPad->GetBottomMargin());
  // tlpp->SetX(xx/gPad->GetWNDC()+gPad->GetLeftMargin());
  //tlpp->SetY(yy/gPad->GetHNDC()+gPad->GetBottomMargin());
  //  tlpp->SetX(0.09/gPad->GetWNDC()+gPad->GetLeftMargin());
  //tlpp->SetY(0.15/gPad->GetHNDC()+gPad->GetBottomMargin());
  //cout<<"innerPadHeight="<<innerPadHeight<<"  referencePadHeight="<<referencePadHeight<<"  resizeTextFactor="<<resizeTextFactor<<"  size="<<size<<"  total="<<size*innerPadHeight/referencePadHeight*resizeTextFactor<<endl;
  tlpp->SetTextFont(43);
  //  TString str=tlpp->GetTitle();
  //  str.Prepend("#font[43]{");
  //  str.Append("}");
  tlpp->SetTextSize(size*innerPadHeight/referencePadHeight*resizeTextFactor);
  //   tlpp->SetTextSize(size/(gPad->GetHNDC())*scaleHeightPads);
  //tlpp->SetTextSize(0.055);

  return;
}

TPaveText *GetALICEpavetext(Int_t identifier){
   TPaveText *alice = new TPaveText(0.157/gPad->GetWNDC()+gPad->GetLeftMargin(),0.195/gPad->GetHNDC()+gPad->GetBottomMargin(),0.29/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
   cout<<"gPad->GetLeftMargin()="<<gPad->GetLeftMargin()<<"  gPad->GetBottomMargin()="<<gPad->GetBottomMargin()<<"   gPad->GetWNDC="<<gPad->GetWNDC()<<"  gPad->GetHNDC()="<<gPad->GetHNDC()<<"  totX="<<0.33/gPad->GetWNDC()+gPad->GetLeftMargin()<<"  totY="<<0.255/gPad->GetHNDC()+gPad->GetBottomMargin()<<endl;
   // TPaveText *alice = new TPaveText(0.78,0.77,0.9,0.83,"NDC");
  //TPaveText *alice = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
//0.72,0.78,0.85,0.83,"NDC");
  SetPaveStyle(alice);
  alice->SetTextFont(43);
  //  alice->SetTextSize(0.07/(gPad->GetHNDC())*scaleHeightPads);
  alice->SetTextSize(26*innerPadHeight/referencePadHeight*resizeTextFactor);
  alice->AddText("ALICE");//commented
  // fitvalueslow->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  alice->SetName(Form("paveALICE_%d",identifier));
  return alice;
}

TPaveText *GetAveragepavetext(Int_t identifier){
  //TPaveText* pvAverage = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.25/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  TPaveText *pvAverage = new TPaveText(0.005/gPad->GetWNDC()+gPad->GetLeftMargin(),0.195/gPad->GetHNDC()+gPad->GetBottomMargin(),0.17/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  // TPaveText * pvAverage= new TPaveText(0.22,0.77,0.64,0.83,"NDC");
  //0.21,0.83,0.5,0.863,"NDC");
  SetPaveStyle(pvAverage);
  pvAverage->SetTextFont(43);
  //  pvAverage->SetTextSize(0.07/(gPad->GetHNDC())*scaleHeightPads);
  pvAverage->SetTextSize(26*innerPadHeight/referencePadHeight*resizeTextFactor);
  pvAverage->AddText("Average D^{0},D^{+},D^{*+}");
  pvAverage->SetName(Form("paveAverage_%d",identifier));
  return pvAverage;
}


TPaveText *GetPaveKineInfo(Int_t iPtMeson,Int_t iPtAssoc,Int_t identifier){
  //TPaveText *GetPaveKineInfo(TString strPtMeson,TString strPtAssoc,Int_t identifier){
  // cout<<" ********* PAD info -> left m="<< gPad->GetLeftMargin()<<"  WNDC="<<gPad->GetWNDC()<<"  HNDC="<<gPad->GetHNDC()<<endl;
   Int_t binrange=4*iPtAssoc+iPtMeson;
   cout<<" ********* PAD info -> left m="<< gPad->GetLeftMargin()<<"  WNDC="<<gPad->GetWNDC()<<"  HNDC="<<gPad->GetHNDC()<<"  binrange = "<<binrange<<endl;
  
   //TPaveText *pvKineInfo = new TPaveText(0.01/gPad->GetWNDC()+gPad->GetLeftMargin(),0.2/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.25/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
   TPaveText *pvKineInfo;
   // pvKineInfo = new TPaveText(0.6,0.64,0.95,0.73,"NDC");

   // if(binrange==1)pvKineInfo = new TPaveText(0.62,0.75,0.9,0.82,"NDC");
   // if(binrange==2)pvKineInfo = new TPaveText(0.62,0.75,0.9,0.82,"NDC");
   
   // if(binrange==3)pvKineInfo = new TPaveText(0.62,0.87,0.98,0.97,"NDC");
   // if(binrange==4)pvKineInfo = new TPaveText(0.8,0.87,0.97,0.97,"NDC");
   // if(binrange==5)pvKineInfo = new TPaveText(0.83,0.87,0.9,0.97,"NDC");

   // if(binrange==6)pvKineInfo = new TPaveText(0.35,0.9,0.9,0.95,"NDC");
   // if(binrange==7)pvKineInfo = new TPaveText(0.42,0.9,0.88,0.95,"NDC");
   // if(binrange==8)pvKineInfo = new TPaveText(0.31,0.9,0.86,0.95,"NDC");

   // if(binrange==0) pvKineInfo = new TPaveText(0.02/gPad->GetWNDC()+gPad->GetLeftMargin(),0.22/gPad->GetHNDC()+gPad->GetBottomMargin(),0.19/gPad->GetWNDC()+gPad->GetLeftMargin(),0.24/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
   // if(binrange==1 || binrange==2)pvKineInfo = new TPaveText(0.02/gPad->GetWNDC()+gPad->GetLeftMargin(),0.255/gPad->GetHNDC()+gPad->GetBottomMargin(),0.19/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
   // if(binrange==3 || binrange==4 || binrange==5)pvKineInfo = new TPaveText(0.0035/gPad->GetWNDC()+gPad->GetLeftMargin(),0.255/gPad->GetHNDC()+gPad->GetBottomMargin(),0.19/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  
   // if(binrange>=6)pvKineInfo = new TPaveText(0.02/gPad->GetWNDC()+gPad->GetLeftMargin(),0.255/gPad->GetHNDC()+gPad->GetBottomMargin(),0.19/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");

   // if(binrange==0) pvKineInfo = new TPaveText(0.013/gPad->GetWNDC()+gPad->GetLeftMargin(),0.21/gPad->GetHNDC()+gPad->GetBottomMargin(),0.19/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
   // if(binrange==1 || binrange==2)pvKineInfo = new TPaveText(0.013/gPad->GetWNDC()+gPad->GetLeftMargin(),0.245/gPad->GetHNDC()+gPad->GetBottomMargin(),0.19/gPad->GetWNDC()+gPad->GetLeftMargin(),0.27/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
    pvKineInfo = new TPaveText(0.004/gPad->GetWNDC()+gPad->GetLeftMargin(),0.245/gPad->GetHNDC()+gPad->GetBottomMargin(),0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.27/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");

   // if(binrange==1)pvKineInfo = new TPaveText(0.62,0.75,0.9,0.82,"NDC");
   // if(binrange==2)pvKineInfo = new TPaveText(0.62,0.75,0.9,0.82,"NDC");
   
   // if(binrange==3)pvKineInfo = new TPaveText(0.62,0.87,0.98,0.97,"NDC");
   // if(binrange==4)pvKineInfo = new TPaveText(0.8,0.87,0.97,0.97,"NDC");
   // if(binrange==5)pvKineInfo = new TPaveText(0.83,0.87,0.9,0.97,"NDC");

   // if(binrange==6)pvKineInfo = new TPaveText(0.35,0.9,0.9,0.95,"NDC");
   // if(binrange==7)pvKineInfo = new TPaveText(0.42,0.9,0.88,0.95,"NDC");
   // if(binrange==8)pvKineInfo = new TPaveText(0.31,0.9,0.86,0.95,"NDC");

 
   SetPaveStyle(pvKineInfo);

   pvKineInfo->SetTextAlign(10);
   //  pvKineInfo->SetTextAlign(12);
   //  pvKineInfo->SetTextSize(20/(gPad->GetHNDC())*scaleHeightPads);
   pvKineInfo->SetTextFont(43);
   pvKineInfo->SetTextSize(18*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
   //  pvKineInfo->SetTextSize(21.9*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight


   // if(binrange==1)pvKineInfo->SetTextSize(0.056/(gPad->GetHNDC())*scaleHeightPads);
   //if(binrange==7)pvKineInfo->SetTextSize(0.065/(gPad->GetHNDC())*scaleHeightPads);
   //if(binrange==8)pvKineInfo->SetTextSize(0.062/(gPad->GetHNDC())*scaleHeightPads);
   pvKineInfo->AddText(strPtRangeText[binrange].Data());
   printf("Binrange %d, string %s\n",binrange,strPtRangeText[binrange]);

   pvKineInfo->SetName(Form("pvKineInf0_%d",identifier));
   return pvKineInfo;
}
TPaveText *GetPaveKineInfo2(Int_t identifier){
  //TPaveText *pvKineInfo = new TPaveText(0.01/gPad->GetWNDC()+gPad->GetLeftMargin(),0.18/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.22/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  //  TPaveText *pvKineInfo = new TPaveText(0.4,0.56,0.7,0.65,"NDC");//in pad 0
  //  TPaveText *pvKineInfo = new TPaveText(0.24,0.56,0.6,0.65,"NDC");//in pad 0
  // TPaveText *pvKineInfo = new TPaveText(0.01,0.6,0.5,0.8,"brNDC");//in pad 1
TPaveText *  pvKineInfo = new TPaveText(0.02/gPad->GetWNDC()+gPad->GetLeftMargin(),0.197/gPad->GetHNDC()+gPad->GetBottomMargin(),0.19/gPad->GetWNDC()+gPad->GetLeftMargin(),0.21/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  
  SetPaveStyle(pvKineInfo);
  pvKineInfo->SetTextAlign(12);
  //pvKineInfo->SetTextSize(0.06/(gPad->GetHNDC())*scaleHeightPads);
  pvKineInfo->SetTextFont(43);
  pvKineInfo->SetTextSize(18*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
  
  pvKineInfo->AddText(strYText.Data());
  //  pvKineInfo->AddText("|#it{y}^{D}|<0.5, |#Delta#it{#eta}| < 1.0 ");
  pvKineInfo->SetName(Form("pvKineInf0_%d",identifier));
  return pvKineInfo;
}
TLegend *GetLegendBaselines(TGraphAsymmErrors *gpp,Int_t identifier){
  TLegend * legend = new TLegend(0.34,0.39,0.55,0.47);//pad 0
  //  TLegend * legend = new TLegend(0.40,0.45,0.85,0.55);//pad 1
  // TLegend * legend = new TLegend(0.40,0.6,0.85,0.67);
  // TLegend * legend = new TLegend(0.40,0.64,0.85,0.72);
  //  TLegend * legend = new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.22/gPad->GetHNDC()+gPad->GetBottomMargin(),0.27/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin());
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextAlign(12);
  // legend->SetTextSize(0.055/(gPad->GetHNDC())*scaleHeightPads);
  legend->SetTextSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
 
  legend->AddEntry(gpp,"baseline-subtraction uncertainty","f");
  legend->SetName(Form("LegendBaselineUncPPandpPb_%d",identifier));
  return legend;
}
TLegend *GetLegendData(TH1D *hpp,Int_t identifier){
    
  // TLegend * legend = new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.13/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.18/gPad->GetHNDC()+gPad->GetBottomMargin());
   TLegend * legend = new TLegend(0.3,0.48,0.6,0.57,NULL,"brNDC");
  
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(22.1*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
 
    legend->AddEntry(hpp,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
    legend->SetName(Form("LegendDataPP_%d",identifier));
    return legend;
  }
TLegend *GetLegendMC(TH1D *hMC1,TH1D *hMC2,TH1D *hMC3,TH1D *hMC4=0x0,TH1D *hMC5=0x0,TH1D *hMC6=0x0,Int_t identifier=0){
    
  // TLegend * legend = new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.1/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.13/gPad->GetHNDC()+gPad->GetBottomMargin());
  Int_t nEffectiveModels=0;
  if(hMC1)nEffectiveModels++;
  if(hMC2)nEffectiveModels++;
  if(hMC3)nEffectiveModels++;
  if(hMC4)nEffectiveModels++;
  if(hMC5)nEffectiveModels++;
  if(hMC6)nEffectiveModels++;

  Double_t ylegMin=0.4,ylegMax=0.75;
  if(hMC1==0x0 || nEffectiveModels>=6)ylegMax=0.675; // no header in TLegend
  if(nEffectiveModels>=6){
    ylegMin=0.37;
  }
  else if(nEffectiveModels<=4){
    ylegMin=0.46;
  }
  if(nEffectiveModels==2){
    ylegMin=0.535;
  }
printf("in max %f %f\n",ylegMin,ylegMax);
  TLegend * legend= new TLegend(0.1,ylegMin,0.7,ylegMax,NULL,"brNDC");
//   if(nEffectiveModels<6){
//     if(hMC1!=0x0) legend= new TLegend(0.1,0.4,0.7,0.75,NULL,"brNDC");
//     else legend= new TLegend(0.1,0.4,0.7,0.68,NULL,"brNDC");
//   }
//   else legend= new TLegend(0.1,0.37,0.7,0.68,NULL,"brNDC");

   //  TLegend * legend = new TLegend(0.25,0.44,0.79,0.84,NULL,"brNDC");
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    //legend->SetTextSize(0.06/(gPad->GetHNDC())*scaleHeightPads);
    if(nEffectiveModels<6)legend->SetTextSize(22.1*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
    else legend->SetTextSize(18.1*innerPadHeight/referencePadHeight*resizeTextFactor);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
    // Draft 5,6: 22.1
    //   with EPOS: 21.5

    if(nEffectiveModels<6){if(hMC1!=0x0)legend->SetHeader("Simulations, pp, #sqrt{#it{s}} = 7 TeV");}
    else {
      TLatex *tlatHeader=new TLatex(0.115,0.71,"Simulations, pp, #sqrt{#it{s}} = 7 TeV");// the legend is 0.115, 71 with standard size
      tlatHeader->SetNDC();
      tlatHeader->SetTextFont(43);
      tlatHeader->SetTextSize(22.1*innerPadHeight/referencePadHeight*resizeTextFactor);
      tlatHeader->Draw();
    }
    if(hMC1)legend->AddEntry(hMC1,"PYTHIA6, Perugia 0","l");
    if(hMC2)legend->AddEntry(hMC2,"PYTHIA6, Perugia 2010","l");
    if(hMC3)legend->AddEntry(hMC3,"PYTHIA6, Perugia 2011","l");   
    if(hMC4)legend->AddEntry(hMC4,"PYTHIA8, Tune 4C","l");
    if(hMC5)legend->AddEntry(hMC5,"POWHEG+PYTHIA6","l");
    if(hMC6)legend->AddEntry(hMC6,"EPOS 3.117","l");
    

    //legend->AddEntry(hMC1,"PYTHIA6, Perugia0","lep");
    //legend->AddEntry(hMC2,"PYTHIA6, Perugia2010","lep");
    //legend->AddEntry(hMC3,"PYTHIA6, Perugia2011","lep");
    //legend->AddEntry(hMC4,"POWHEG+PYTHIA6","lep");

    legend->SetName(Form("LegendMC_%d",identifier));
    return legend;
  }


TLegend *GetLegendBaselinesSinglePanel(TGraphAsymmErrors *gpp,Int_t identifier){
  TLegend * legend = new TLegend(0.3,0.65,0.6,0.68,NULL,"brNDC");//pad 0
  //  TLegend * legend = new TLegend(0.40,0.45,0.85,0.55);//pad 1
  // TLegend * legend = new TLegend(0.40,0.6,0.85,0.67);
  // TLegend * legend = new TLegend(0.40,0.64,0.85,0.72);
  //  TLegend * legend = new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.22/gPad->GetHNDC()+gPad->GetBottomMargin(),0.27/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin());
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextAlign(12);
  // legend->SetTextSize(0.055/(gPad->GetHNDC())*scaleHeightPads);
  legend->SetTextSize(21);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
 
  legend->AddEntry(gpp,"baseline-subtraction uncertainty","f");
  legend->SetName(Form("LegendBaselineUncPPandpPb_%d",identifier));
  return legend;
}
TLegend *GetLegendDataSinglePanel(TH1D *hpp,Int_t identifier){
    
  // TLegend * legend = new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.13/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.18/gPad->GetHNDC()+gPad->GetBottomMargin());
   TLegend * legend = new TLegend(0.3,0.68,0.6,0.71,NULL,"brNDC");
  
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(21);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
 
    legend->AddEntry(hpp,"pPb, #sqrt{#it{s}_{NN}} = 5.02 TeV","lep");
    legend->SetName(Form("LegendDataPP_%d",identifier));
    return legend;
  }
TLegend *GetLegendMCSinglePanel(TH1D *hMC1,TH1D *hMC2,TH1D *hMC3,TH1D *hMC4=0x0,TH1D *hMC5=0x0,TH1D *hMC6=0x0,Int_t identifier=0){
    
  // TLegend * legend = new TLegend(0.011/gPad->GetWNDC()+gPad->GetLeftMargin(),0.1/gPad->GetHNDC()+gPad->GetBottomMargin(),0.2/gPad->GetWNDC()+gPad->GetLeftMargin(),0.13/gPad->GetHNDC()+gPad->GetBottomMargin());
  Int_t nEffectiveModels=0;
  if(hMC1)nEffectiveModels++;
  if(hMC2)nEffectiveModels++;
  if(hMC3)nEffectiveModels++;
  if(hMC4)nEffectiveModels++;
  if(hMC5)nEffectiveModels++;
  if(hMC6)nEffectiveModels++;

  Double_t ylegMin=0.4,ylegMax=0.65;
  if(hMC1==0x0 || nEffectiveModels>=6)ylegMax=0.675; // no header in TLegend
  if(nEffectiveModels>=6){
    ylegMin=0.37;
  }
  else if(nEffectiveModels<=4){
    ylegMin=0.46;
  }
  if(nEffectiveModels==2){
    ylegMin=0.535;
  }
printf("in max %f %f\n",ylegMin,ylegMax);
  TLegend * legend= new TLegend(0.3,ylegMin,0.6,ylegMax,NULL,"brNDC");
//   if(nEffectiveModels<6){
//     if(hMC1!=0x0) legend= new TLegend(0.1,0.4,0.7,0.75,NULL,"brNDC");
//     else legend= new TLegend(0.1,0.4,0.7,0.68,NULL,"brNDC");
//   }
//   else legend= new TLegend(0.1,0.37,0.7,0.68,NULL,"brNDC");

   //  TLegend * legend = new TLegend(0.25,0.44,0.79,0.84,NULL,"brNDC");
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    //legend->SetTextSize(0.06/(gPad->GetHNDC())*scaleHeightPads);
    if(nEffectiveModels<6)legend->SetTextSize(21);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
    else legend->SetTextSize(21);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
    // Draft 5,6: 22.1
    //   with EPOS: 21.5

    if(nEffectiveModels<6){if(hMC1!=0x0)legend->SetHeader("Simulations, pp, #sqrt{#it{s}} = 7 TeV");}
    else {
      TLatex *tlatHeader=new TLatex(0.115,0.71,"Simulations, pp, #sqrt{#it{s}} = 7 TeV");// the legend is 0.115, 71 with standard size
      tlatHeader->SetNDC();
      tlatHeader->SetTextFont(43);
      tlatHeader->SetTextSize(21);
      tlatHeader->Draw();
    }
    if(hMC1)legend->AddEntry(hMC1,"PYTHIA6, Perugia 0","l");
    if(hMC2)legend->AddEntry(hMC2,"PYTHIA6, Perugia 2010","l");
    if(hMC3)legend->AddEntry(hMC3,"PYTHIA6, Perugia 2011","l");   
    if(hMC4)legend->AddEntry(hMC4,"PYTHIA8, Tune 4C","l");
    if(hMC5)legend->AddEntry(hMC5,"POWHEG+PYTHIA6","l");
    if(hMC6)legend->AddEntry(hMC6,"EPOS 3.117","l");
    

    //legend->AddEntry(hMC1,"PYTHIA6, Perugia0","lep");
    //legend->AddEntry(hMC2,"PYTHIA6, Perugia2010","lep");
    //legend->AddEntry(hMC3,"PYTHIA6, Perugia2011","lep");
    //legend->AddEntry(hMC4,"POWHEG+PYTHIA6","lep");

    legend->SetName(Form("LegendMC_%d",identifier));
    return legend;
  }

TPaveText *GetALICEpavetextSinglePanel(Int_t identifier){
   TPaveText *alice = new TPaveText(0.2,0.81,0.4,0.87,"NDC");
   cout<<"gPad->GetLeftMargin()="<<gPad->GetLeftMargin()<<"  gPad->GetBottomMargin()="<<gPad->GetBottomMargin()<<"   gPad->GetWNDC="<<gPad->GetWNDC()<<"  gPad->GetHNDC()="<<gPad->GetHNDC()<<"  totX="<<0.33/gPad->GetWNDC()+gPad->GetLeftMargin()<<"  totY="<<0.255/gPad->GetHNDC()+gPad->GetBottomMargin()<<endl;
  SetPaveStyle(alice);
  alice->SetTextFont(43);
  alice->SetTextSize(30);
  alice->AddText("ALICE");//commented
  alice->SetName(Form("paveALICE_%d",identifier));
  return alice;
}

TPaveText *GetAveragepavetextSinglePanel(Int_t identifier){

  TPaveText *pvAverage = new TPaveText(0.2,0.75,0.4,0.82,"NDC");
  SetPaveStyle(pvAverage);
  pvAverage->SetTextFont(43);
  pvAverage->SetTextSize(30);
  pvAverage->AddText("Average D^{0},D^{+},D^{*+}");
  pvAverage->SetName(Form("paveAverage_%d",identifier));
  return pvAverage;
}


TPaveText *GetPaveKineInfoSinglePanel(Int_t iPtMeson,Int_t iPtAssoc,Int_t identifier){
   Int_t binrange=4*iPtAssoc+iPtMeson;
   cout<<" ********* PAD info -> left m="<< gPad->GetLeftMargin()<<"  WNDC="<<gPad->GetWNDC()<<"  HNDC="<<gPad->GetHNDC()<<"  binrange = "<<binrange<<endl;
   TPaveText *pvKineInfo = new TPaveText(0.4,0.81,0.9,0.86,"NDC");
 
   SetPaveStyle(pvKineInfo);
   pvKineInfo->SetTextFont(43);
   pvKineInfo->SetTextSize(30);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
   pvKineInfo->AddText(strPtRangeText[binrange].Data());
   printf("Binrange %d, string %s\n",binrange,strPtRangeText[binrange]);
   pvKineInfo->SetName(Form("pvKineInf0_%d",identifier));
   return pvKineInfo;
}

TPaveText *GetPaveKineInfo2SinglePanel(Int_t identifier){
  TPaveText *  pvKineInfo = new TPaveText(0.4,0.75,0.9,0.82,"NDC");
  SetPaveStyle(pvKineInfo);
  pvKineInfo->SetTextFont(43);
  pvKineInfo->SetTextSize(30);// settings for font 42: 0.07/(gPad->GetHNDC())*scaleHeight
  pvKineInfo->AddText(strYText.Data());
  pvKineInfo->SetName(Form("pvKineInf0_%d",identifier));
  return pvKineInfo;
}

void DoComparison_pPb2016VsMCSinglePanel(){

  Init();
  LoadFileNamesAll();
  
  if(!includeset[0]){
    Printf("The pPb dataset is expected to be included! Cannot proceed"); return;
  }
  
  TH1D *h;
  // for(Int_t iset=0;iset<1;iset++){
  for(Int_t iset=0;iset<nSets;iset++){
    if(!includeset[iset])continue;
    for(Int_t kassoc=0;kassoc<nbinAssocpt;kassoc++){
      for(Int_t jmes=firstDpt;jmes<nbinDpt;jmes++){

  if(iset==0){
    histo[iset][kassoc][jmes]=GetHistoAndSyst(filenames[iset][kassoc][jmes],iset,"fhDaverage","AverageSystematicUncertainty",err[iset][kassoc][jmes],ltscale[iset][kassoc][jmes]); 
  
    
  // pp style
    histo[iset][kassoc][jmes]->SetMarkerColor(4); 
    histo[iset][kassoc][jmes]->SetMarkerStyle(21); 
    histo[iset][kassoc][jmes]->SetLineColor(4);  
    
    
    subtractedhisto[iset][kassoc][jmes] = GetPedestalHistoAndSystAndSubtractPed(iset,kassoc,jmes,histo[iset][kassoc][jmes],err[iset][kassoc][jmes],suberr[iset][kassoc][jmes],"CanvasBaselineVariationTrendPedestal",grbase[iset][kassoc][jmes],grv2[iset][kassoc][jmes]);
    Printf("Histo subtrcated obtained");
    cout<<"sub -> "<<subtractedhisto[iset][kassoc][jmes]->GetBinContent(5)<<endl; 
    Printf("Histo subtrcated: pointer is working");
    if(grv2[iset][kassoc][jmes]){
      Printf("SHOULD NOT ENTER HERE");
      grv2[iset][kassoc][jmes]->SetFillStyle(3002);
      grv2[iset][kassoc][jmes]->SetFillColor(kMagenta);  
    }
    suberr[iset][kassoc][jmes]->SetLineColor(kBlack);
  }
  

  else{
    h=GetHisto(filenames[iset][kassoc][jmes],"hCorrDeltaPhi");
    if(reflTempl)histo[iset][kassoc][jmes]=AliHFCorrelationUtils::ReflectHisto(h,0.5);
    else histo[iset][kassoc][jmes]=h;

      histo[iset][kassoc][jmes]->SetMarkerColor(modelColors[iset-1]);// -1 because first iset is data
      histo[iset][kassoc][jmes]->SetLineColor(modelColors[iset-1]); 
      histo[iset][kassoc][jmes]->SetLineWidth(2); 
      histo[iset][kassoc][jmes]->SetMarkerStyle(kDot); 
      //  histo[iset][kassoc][jmes]->SetMarkerStyle(kOpenSquare); 

     subtractedhisto[iset][kassoc][jmes] = GetPedestalHistoAndSystAndSubtractPedMC(iset,kassoc,jmes,histo[iset][kassoc][jmes],"CanvasBaselineVariationTrendPedestal");
     
  }

  
        subtractedhisto[iset][kassoc][jmes]->SetMinimum(minYaxis[kassoc]);
        subtractedhisto[iset][kassoc][jmes]->SetMaximum(maxYaxis[kassoc]);


      }      
    }
  }

  Printf("All Histos and Graphs created");
  // UP TO HERE SHOULD BE OK

  for(Int_t iassoc=0;iassoc<nbinAssocpt;iassoc++){

    for(Int_t jDpt=firstDpt;jDpt<nbinDpt;jDpt++){

      TCanvas *cFinalPaperStyle;

      cFinalPaperStyle=new TCanvas("cFinalPaperStyle","cFinalPaperStyle",1200.,900.);
      cFinalPaperStyle->SetTicky();
      cFinalPaperStyle->SetTickx();
      cFinalPaperStyle->SetFillStyle(0);
      cFinalPaperStyle->SetFrameFillStyle(4000);
      cFinalPaperStyle->SetFrameBorderMode(0);
      cFinalPaperStyle->Modified();
      cFinalPaperStyle->Update();
    
      Double_t startline,endline,wi;
      wi=subtractedhisto[0][iassoc][jDpt]->GetBinWidth(1);
      startline=subtractedhisto[0][iassoc][jDpt]->GetBinLowEdge(1);
      // endline=subtractedhisto[0][iassoc][jDpt]->GetBinLowEdge(subtractedhisto[0][iassoc][jDpt]->GetNbinsX())+wi+0.4;
      endline=subtractedhisto[0][iassoc][jDpt]->GetBinLowEdge(subtractedhisto[0][iassoc][jDpt]->GetNbinsX())+wi+0.15;
      // cout<<"******************* "<<startline<<"  "<<endline<<endl;
      TLine* line=new TLine(startline, 0, endline, 0);
      line->SetLineStyle(2);
      
      TH1D* h=new TH1D(*subtractedhisto[0][iassoc][jDpt]);

      h->Draw();

      //      grbase[0][iassoc][jDpt]->SetFillStyle(3002);
      grbase[0][iassoc][jDpt]->SetFillColor(kGray+2);
      grbase[0][iassoc][jDpt]->SetLineColor(kGray+2);
      grbase[0][iassoc][jDpt]->Draw("E2");
/* //DISABLED BY FABIO - COMPLETELY NEGLIGIBLE UNCERTAINTIES...
      grv2[0][iassoc][jDpt]->SetFillColor(kGreen+1);
      grv2[0][iassoc][jDpt]->SetLineColor(kGreen+1);
      grv2[0][iassoc][jDpt]->Draw("E2");
*/
      line->Draw(); 
      subtractedhisto[0][iassoc][jDpt]->Draw("same");

      suberr[0][iassoc][jDpt]->SetLineColor(kBlack);
      suberr[0][iassoc][jDpt]->Draw("E2");

      Float_t size=0.058;
      SetScaleUncertaintyPositionAndSize(ltscale[0][iassoc][jDpt],0.08,0.104,22);//0.08
     

      ltscale[0][iassoc][jDpt]->SetX(0.65);
      ltscale[0][iassoc][jDpt]->SetY(0.67);
      ltscale[0][iassoc][jDpt]->SetTextFont(43);

      ltscale[0][iassoc][jDpt]->SetTextSize(25);


      //     if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1){
      //ltscale[0][iassoc][jDpt]->SetY(0.11/gPad->GetHNDC()+gPad->GetBottomMargin());
      //}
      ltscale[0][iassoc][jDpt]->Draw();
      for(Int_t kmod=1;kmod<nSets;kmod++){
        if(includeset[kmod])subtractedhisto[kmod][iassoc][jDpt]->Draw("hist same c");
      }
//       subtractedhisto[1][iassoc][jDpt]->Draw("hist same c");//Perugia0
//       subtractedhisto[2][iassoc][jDpt]->Draw("hist same c");//Perugia2010
//       subtractedhisto[3][iassoc][jDpt]->Draw("hist same c");//Perugia2011
//       subtractedhisto[4][iassoc][jDpt]->Draw("hist same c");//PYTHIA8
//       subtractedhisto[5][iassoc][jDpt]->Draw("hist same c");//POWHEG
//       subtractedhisto[6][iassoc][jDpt]->Draw("hist same c");//EPOS3

     // subtractedhisto[1][iassoc][jDpt]->Draw("same");//Perugia0
     //  subtractedhisto[2][iassoc][jDpt]->Draw("same");//Perugia2010
     //  subtractedhisto[3][iassoc][jDpt]->Draw("same");//Perugia2011
     //  subtractedhisto[4][iassoc][jDpt]->Draw("same");//PYTHIA8
     //  subtractedhisto[5][iassoc][jDpt]->Draw("same");//POWHEG
      subtractedhisto[0][iassoc][jDpt]->Draw("same");//pp DATA 
      suberr[0][iassoc][jDpt]->Draw("E2");
    
      TPaveText *pvAverage=GetAveragepavetextSinglePanel((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TPaveText *alice=GetALICEpavetextSinglePanel((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);  
      //TPaveText *pvKineInfo = GetPaveKineInfo(strPtMesonText[jDpt],strPtAssocText[iassoc],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TPaveText *pvKineInfo = GetPaveKineInfoSinglePanel(jDpt,iassoc,(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TPaveText *pvKineInfo2 = GetPaveKineInfo2SinglePanel((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      TLegend *legendbase=GetLegendBaselinesSinglePanel(grbase[0][iassoc][firstDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
      legendbase->Draw();
      alice->Draw();
      pvKineInfo->Draw("same");
      pvAverage->Draw("same");
      //if((nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1==1)pvKineInfo2->Draw("same");
      pvKineInfo2->Draw("same");
      legendbase->Draw();
      //   if(jDpt==1 && iassoc==0)pvKineInfo2->Draw("same");

      // N.B. GetLegendData does not assume that the pointers to histograms that are passed are not null --> no need to check which model is included 

        TLegend *legendData=GetLegendDataSinglePanel(subtractedhisto[0][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
        legendData->Draw();
  
        TLegend *legendMC;
        legendMC=GetLegendMCSinglePanel(subtractedhisto[1][iassoc][jDpt],subtractedhisto[2][iassoc][jDpt],subtractedhisto[3][iassoc][jDpt],subtractedhisto[4][iassoc][jDpt],subtractedhisto[5][iassoc][jDpt],subtractedhisto[6][iassoc][jDpt],(nbinDpt-firstDpt)*iassoc+jDpt-firstDpt+1);
        legendMC->Draw();


      cFinalPaperStyle->SaveAs(Form("cFinalPaperStyle_SingleKineRange_%s_%s.png",strmesonpt[jDpt].Data(),pthadron[iassoc].Data()));
      cFinalPaperStyle->SaveAs(Form("cFinalPaperStyle_SingleKineRange_%s_%s.root",strmesonpt[jDpt].Data(),pthadron[iassoc].Data()));

    }
  }
  
}
