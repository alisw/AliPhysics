//############################################################################
// 
//  Macros to perform and draw the comparison of D-haron correlation results
//  in pp and p-Pb collisions
//  elena.bruna@to.infn.it
//
//############################################################################
Int_t nhistos=6;
Int_t rebin =1;
Double_t **ptbinsD;
Int_t ndesiredbins=3;
TString inputdatadirectory = "./";
TString baselinedirectory="./";
TString avType="Weighted";
TString *filenames;
TString *pedestalfilenames;
TH1D ** histo;
TGraphAsymmErrors ** err;
TLatex** ltscale;
Bool_t skip3to5=kTRUE;
Int_t ihskip=0;
TString fitplotmacrodir=gSystem->ExpandPathName("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/");
Bool_t isReflected=kFALSE;
Double_t leftMarginCanvas=0.17;
Double_t rightMarginCanvas=0.055;
Double_t bottomMarginCanvas=0.1;
Double_t topMarginCanvas=0.1;
Color_t colpp=kBlack;
Color_t colpPb=kRed;

void SetIsReflected(Bool_t isrefl){
  isReflected=isrefl;
}
void SetDesiredPtBins(Int_t nbins,Double_t **bins){// maybe useful in future to replace hard-coded values, not used right now
  ndesiredbins=nbins;
  ptbinsD=new Double_t[ndesiredbins];
  for(Int_t j=0;j<nbins;j++){
    ptbinsD[j]=bins[j];
  }
}
void SetFitPlotMacroPath(TString strdir){
  fitplotmacrodir=strdir;
}

void SetSkip3to5pPb(Bool_t skip){
  skip3to5=skip;
  if(skip)ihskip=1;
  else ihskip=0;
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
  pv->SetTextFont(42);
  pv->SetTextAlign(12);

}
TCanvas *GetCanvasWithStyle(TString name,TString title){
  TCanvas *c=new TCanvas(name.Data(),title.Data(),800,800);
  c->SetTicky();
  c->SetTickx();

  c->SetFrameBorderMode(0);
  c->Range(-2.6,-1.5,5.4,4.4);
  c->SetLeftMargin(leftMarginCanvas);
  c->SetRightMargin(rightMarginCanvas);
  c->SetBottomMargin(bottomMarginCanvas);
  c->SetTopMargin(topMarginCanvas);

  return c;
}

void Init(){
  if(filenames)delete [] filenames;
  if(pedestalfilenames)delete [] pedestalfilenames;
  if(histo){
    for(Int_t k=0;k<nhistos;k++){
      delete histo[k];
    }
    delete [] histo;
  }
  if(err){
    for(Int_t k=0;k<nhistos;k++){
      delete err[k];
    }
    delete [] err;
  }
  if(ltscale){
    for(Int_t k=0;k<nhistos;k++){
      delete ltscale[k];
    }
    delete [] ltscale;
  }
  

  filenames = new TString[nhistos];
  pedestalfilenames = new TString[2];
  
  histo = new TH1D *[nhistos];
  err = new TGraphAsymmErrors *[nhistos];
  ltscale = new TLatex *[nhistos];
  
}

Int_t rebin =1;
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

void test(){
  TH1D* h1;
  TGraphAsymmErrors* gr;
  TLatex* lt;

  h1=(TH1D*)OpenOutputFileAndDraw("Average_D_pPb_Final_Plots/AverageDzeroDstarDplus8to16_assoc1.0.root",8,16,"D",1,1,1,10,&gr,&lt);

  TCanvas* c=new TCanvas();
  h1->Draw();
  gr->Draw("E2");
  lt->Draw();
}


TH1D* OpenOutputFileAndDraw(TString strfile,Double_t ptminD,Double_t ptmaxD,TString strMeson="D",Double_t ptminAss=0.3,Double_t deltaeta=1, Int_t system=0, Double_t max=10.,TGraphAsymmErrors** gr2, TLatex** lt){
  gStyle->SetOptStat(0000);
  gROOT->LoadMacro(Form("%s/FitPlots.C",fitplotmacrodir.Data()));
 
  TFile *f=TFile::Open(strfile.Data(),"READ");
  AliHFDhadronCorrSystUnc *syst=(AliHFDhadronCorrSystUnc*)f->Get("AverageSystematicUncertainty");
  TH1D *hUncCorrMin=syst->GetHistoTotFlatMin();
  TH1D *hUncCorrMax=syst->GetHistoTotFlatMax();
  TH1D *hFDsub=(TH1D*)f->Get("fhDaverage");
  TGraphAsymmErrors *gr=syst->GetTotNonFlatUncGraph();
  gStyle->SetOptStat(0000);
  (*gr2)=gr;

  TCanvas *cDraw=new TCanvas("cDraw","cDraw",700,700);
  cDraw->cd();
  cDraw->SetLeftMargin(0.15);
  cDraw->SetRightMargin(0.05);
  cDraw->SetTicks();

  hFDsub->SetLineColor(kBlack);
  hFDsub->SetMarkerColor(kBlack);
  hFDsub->SetXTitle("#Delta#varphi(D,h) (rad)");
  hFDsub->SetYTitle(Form("#frac{1}{N_{%s}}#frac{dN^{assoc}}{d#Delta#varphi} (rad^{-1})",strMeson.Data()));
  hFDsub->GetYaxis()->SetTitleOffset(1.3);
  hFDsub->GetYaxis()->SetRangeUser(0,max);
  hFDsub->SetTitle("");
  gr->SetLineColor(kBlack);
  gr->SetMarkerColor(kBlack);
  gr->SetFillStyle(0);
  hFDsub->Draw();
  gr->Draw("E2");

  TLatex *tlTitle=new TLatex(0.18,0.85,Form("#bf{%s meson-charged hadron azimuthal correlations}",strMeson.Data()));
  tlTitle->SetNDC();
  tlTitle->Draw();
  tlTitle->SetTextSize(0.03);

  TLatex *tSystem=new TLatex(0.18,0.80,"#bf{pp, #sqrt{s}=7 TeV, L_{int} = 5 nb^{-1}}");
  if(system==1) tSystem->SetTitle("#bf{p-Pb, #sqrt{s_{NN}}=5.02 TeV, L_{int} = 50 #mub^{-1}}");
  tSystem->SetNDC();
  tSystem->SetTextSize(0.03);
  tSystem->Draw();
  
  TLatex *tptD=new TLatex(0.18,0.75,Form("#bf{%.0f < p_{T}^{%s} < %.0f GeV/c,  p_{T}^{assoc} > %.1f GeV/c}",ptminD,ptmaxD,strMeson.Data(),ptminAss));
  tptD->SetNDC();
  tptD->SetTextSize(0.03);
  tptD->Draw();

  TLatex *tDEta=new TLatex(0.18,0.70,Form("#bf{|#Delta#eta|<%.1f}",deltaeta));
  tDEta->SetNDC();
  tDEta->SetTextSize(0.03);
  tDEta->Draw();


  TLatex *tlTitle=new TLatex(0.60,0.15,"ALICE Preliminary");
  tlTitle->SetNDC();
  tlTitle->Draw();
  tlTitle->SetTextFont(42);
  tlTitle->SetTextSize(0.04);

  TLatex *tUncertainty;
  if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.2,0.6,Form("#bf{%.0f#% scale uncertainty}",hUncCorrMin->GetBinContent(1)*100.));
  else tUncertainty=new TLatex(0.40,0.65,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.04);
  tUncertainty->Draw();
  tUncertainty->SetTextFont(62);
  (*lt)=tUncertainty;
  
  TCanvas *cVaryHisto=new TCanvas("cVaryHisto","cVaryHisto",700,700);
  cVaryHisto->cd();

  hFDsub->DrawCopy();

  TH1D *hVaryUp=syst->GetVariedHisto(hFDsub,gr,1);
  TH1D *hVaryDown=syst->GetVariedHisto(hFDsub,gr,0);
  hVaryUp->Draw("same");
  hVaryDown->Draw("same");

  TString strfileout="CanvaAndVariedHisto";
  strfileout.Append(Form("AverageDzeroDstarPt%.0fto%.0fassocPt%.1f.root",ptminD,ptmaxD,ptminAss));
  
  cDraw->Update();
  TFile *fout=new TFile(strfileout.Data(),"RECREATE");
  strfileout.ReplaceAll(".root",".png");
  fout->cd();
  cDraw->Write();
  cDraw->Print(strfileout.Data());
  hVaryUp->Write();
  hVaryDown->Write();
  fout->Close();
  return hFDsub;
}


//_______________________________________________________________________
TH1D * GetHisto(Int_t i){
    
  // AliDCorrelationPlotter * plotter = new AliDCorrelationPlotter();
    
  TString path = filenames[i];
    
  cout << "Reading File from path: " << path << endl;
    
  TFile * file = TFile::Open(path.Data(),"READ");
  TH1D * histo = (TH1D*)file->Get("hDataCorrectedTempl0CentrFprompt");
  if(!histo) {
    cout << "smtg is wrong" << endl;
    file->ls();
    return 0;
  }
    
  TString histoname = "DeltaPhiBfeedDown_";
  histoname += Form("%d",i);
  TH1D * outputhisto = (TH1D*)histo->Clone("blabla");
  outputhisto->SetName(histoname.Data());
    
    
    
    
  //  TH1D * histo2 = plotter->ReflectHistogram(histo);
    
  return histo;
    
}

//_______________________________________________________________________
TH1D * GetPedestalHisto(Int_t i, TString canvasname, TString hname){
    
  // get pedestal from fit outputs
  TString path = pedestalfilenames[i];
    
  cout << "Reading File from path: " << path << endl;
    
  TFile * file = TFile::Open(path.Data(),"WRITE");
  TCanvas * c = (TCanvas*)file->Get(canvasname.Data());
  c->cd();
  
  TH1D * histo = (TH1D*)c->FindObject(hname.Data());
    
  TString histoname = "Pedestal_";
  histoname += Form("%d",i);
  TH1D * outputhisto = (TH1D*)histo->Clone("blabla");
  outputhisto->SetName(histoname.Data());

  return outputhisto;
    
}
//_______________________________________________________________________
TH1D * GetPedestalHistoAndSystAndSubtractPedpPb(TString system="pPb",Int_t i,TH1D *histo, TGraphAsymmErrors* gr,TGraphAsymmErrors** grout, TString canvasname,TBox** bsyst, TBox** bsyst2){// OBSOLETE METHOD

  //i=D pt bin/
  //Double_t* arrbox1, Double_t* arrbox2
  TBox* bout1;
  TBox* bout2;

  Double_t value = 0, pedestal=0;
  if(system.Contains("pPb")){
    // get pedestal from fit outputs
    TString path = pedestalfilenames[1];//pPb ----------------
    
    cout << "pPb -->  Reading File from path: " << path << endl;
    
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    TBox* b1=(TBox*)c->GetListOfPrimitives()->At(2);
    TBox* b2=(TBox*)c->GetListOfPrimitives()->At(3);
    TBox* b3=(TBox*)c->GetListOfPrimitives()->At(4);
    TBox* b4=(TBox*)c->GetListOfPrimitives()->At(5);
    

    //example:
    //    TBox  X1= 5.600000 Y1=6.724682 X2=7.400000 Y2=6.802155     bin 5-8 : v2
    //    TBox  X1= 11.100000 Y1=6.563838 X2=12.900000 Y2=6.638230   bin 8-16 :v2
    //    TBox  X1= 5.800000 Y1=6.146429 X2=7.200000 Y2=7.380408     bin 5-8 : Total (no v2 modulation)
    //    TBox  X1= 11.300000 Y1=5.990975 X2=12.700000 Y2=7.211093   bin 8-16 : Total (no v2 modulation)
 
     
    if(i==3)pedestal=h->GetBinContent(1);
    if(i==4)pedestal=h->GetBinContent(2);
    if(i==5)pedestal=h->GetBinContent(3);
    cout<<i<<"  pedestal pPb="<<pedestal<<endl;

    bout1=(TBox*)b1->Clone("b1");//Total
    bout2=(TBox*)b1->Clone("b2");//v2

    if(i==3){
      bout2->SetX1(0);
      bout2->SetX2(0);
      bout2->SetY1(0);
      bout2->SetY2(0);

      bout1->SetX1(0);
      bout1->SetX2(0);
      bout1->SetY1(0);
      bout1->SetY2(0);
    }
    if(i==4){
      bout2->SetX1(-1.5);
      bout2->SetX2(-1.2);
      bout2->SetY1(b1->GetY1()-pedestal);
      bout2->SetY2(b1->GetY2()-pedestal);

      bout1->SetX1(4.9);
      bout1->SetX2(5.05);
      //      bout1->SetX1(-1.4);
      //bout1->SetX2(-1.3);
      bout1->SetY1(b3->GetY1()-pedestal);
      bout1->SetY2(b3->GetY2()-pedestal);
    }

    if(i==5){
      bout2->SetX1(-1.5);
      bout2->SetX2(-1.2);
      bout2->SetY1(b2->GetY1()-pedestal);
      bout2->SetY2(b2->GetY2()-pedestal);

      bout1->SetX1(4.9);
      bout1->SetX2(5.05);
      //      bout1->SetX1(-1.4);
      //bout1->SetX2(-1.3);
      bout1->SetY1(b4->GetY1()-pedestal);
      bout1->SetY2(b4->GetY2()-pedestal);
    }

    cout << "pPb --> subtracting pedestal for histo " << i << " = " << pedestal << endl;
    cout<<"BOX1 = "<<bout1->GetX1()<<"  "<<bout1->GetX2()<<"  "<<bout1->GetY1()<<"  "<<bout1->GetY2()<<"  "<<endl;
    cout<<"BOX2 = "<<bout2->GetX1()<<"  "<<bout2->GetX2()<<"  "<<bout2->GetY1()<<"  "<<bout2->GetY2()<<"  "<<endl;
  }
    
  if(system.Contains("pp")){
    // get pedestal from fit outputs
    TString path = pedestalfilenames[0];//pp ----------------
    
    cout << "pp -->  Reading File from path: " << path << endl;
      
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    TBox* b1=(TBox*)c->GetListOfPrimitives()->At(2);//3-5
    TBox* b2=(TBox*)c->GetListOfPrimitives()->At(3);//5-8
    TBox* b3=(TBox*)c->GetListOfPrimitives()->At(4);//8-16
    bout1=(TBox*)b1->Clone("b1");
    bout2=(TBox*)b1->Clone("b2");

    if(i==0)pedestal=h->GetBinContent(1);
    if(i==1)pedestal=h->GetBinContent(2);
    if(i==2)pedestal=h->GetBinContent(3);
    cout<<i<<"  pedestal pp="<<pedestal<<endl;

    if(i==0){
      bout1->SetX1(-1..2);
      bout1->SetX2(-1.1);
      bout1->SetY1(b1->GetY1()-pedestal);
      bout1->SetY2(b1->GetY2()-pedestal);

      bout2->SetX1(0);
      bout2->SetX2(0);
      bout2->SetY1(0);
      bout2->SetY2(0);
    }
    if(i==1){
      //  bout1->SetX1(4.2);
      //bout1->SetX2(4.5);
      bout1->SetX1(4.75);
      bout1->SetX2(4.9);
      //      bout1->SetX1(-1.2);
      //bout1->SetX2(-1.1);
      bout1->SetY1(b2->GetY1()-pedestal);
      bout1->SetY2(b2->GetY2()-pedestal);
	  
      bout2->SetX1(0);
      bout2->SetX2(0);
      bout2->SetY1(0);
      bout2->SetY2(0);
    }

    if(i==2){
      bout1->SetX1(4.75);
      bout1->SetX2(4.9);
      //      bout1->SetX1(-1.2);
      //bout1->SetX2(-1.1);
      bout1->SetY1(b3->GetY1()-pedestal);
      bout1->SetY2(b3->GetY2()-pedestal);
	  
      bout2->SetX1(0);
      bout2->SetX2(0);
      bout2->SetY1(0);
      bout2->SetY2(0);
    }

    cout << "pp --> subtracting pedestal for histo " << i << " = " << pedestal << endl;
    cout<<"BOX1 = "<<bout1->GetX1()<<"  "<<bout1->GetX2()<<"  "<<bout1->GetY1()<<"  "<<bout1->GetY2()<<"  "<<endl;
    cout<<"BOX2 = "<<bout2->GetX1()<<"  "<<bout2->GetX2()<<"  "<<bout2->GetY1()<<"  "<<bout2->GetY2()<<"  "<<endl;
  }
  
  (*bsyst)=bout1;
  (*bsyst2)=bout2;

  TGraphAsymmErrors* gr2=(TGraphAsymmErrors*)gr->Clone("gr_out");
  
  TString nameoutput = histo->GetName();
  nameoutput += "_subtr_";
  nameoutput += "pedestal";
  nameoutput += Form("%d",i);
   
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
    gr2->SetPoint(iBin-1,x,y-pedestal);

  }
  (*grout)=gr2;
  cout<<"sub -> "<<outputhisto->GetBinContent(5)<<endl;
  cout<<"now return histogram"<<endl;
  
  outputhisto->SetXTitle("#Delta#varphi (D,h) (rad)");
  outputhisto->SetYTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} - baseline (rad^{-1})");
  outputhisto->GetYaxis()->SetTitleOffset(1.5);
  return outputhisto;
    
}

//_______________________________________________________________________
TH1D * GetHisto(Int_t i, TString canvasname, TString hname){
    
  //load the histogram with
  TString path = filenames[i];
    
  cout << "Reading File from path: " << path << endl;
    
  TFile * file = TFile::Open(path.Data(),"WRITE");
  //if(i ==2){cout << "file content = " << endl; file->ls();}
  TCanvas * c = (TCanvas*)file->Get(canvasname.Data());
  c->cd();
  // if(i ==2) {
  //   cout << "canvas content = " << endl; c->ls(); }
  TH1D * histo = (TH1D*)c->FindObject(hname.Data());
    
  TString histoname = "DeltaPhi_";
  histoname += Form("%d",i);
  TH1D * outputhisto = (TH1D*)histo->Clone("blabla");
  outputhisto->SetName(histoname.Data());
  // cout << "Scaling By Bin Width = " << outputhisto->GetBinWidth(3) << endl;
  // outputhisto->Scale(1./outputhisto->GetBinWidth(3));
    
  //  TH1D * histo2 = plotter->ReflectHistogram(histo);
    
  return histo;
    
}

//_______________________________________________________________________
TH1D * GetHistoAndSyst(Int_t i, TString hname, TString hnamesyst,TGraphAsymmErrors *&gr2, TLatex *&tUncertainty){
    
   
  //load the histogram with
  TString path = filenames[i];
    
  cout << "Reading File from path: " << path << endl;
    
  //TFile * file = TFile::Open(path.Data(),"WRITE");
  //if(i ==2){cout << "file content = " << endl; file->ls();}
  //TCanvas * c = (TCanvas*)file->Get(canvasname.Data());
  //c->cd();
  // if(i ==2) {
  //   cout << "canvas content = " << endl; c->ls(); }
  //    TH1D * histo = (TH1D*)c->FindObject(hname.Data());

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

  //  (*gr2)=gr;

  if(i<3){
    if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.55,0.46,Form("#bf{%.0f#% scale uncertainty pp}",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.55,0.46,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty (pp)}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
   
  }
if(i>=3){
    if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.55,0.42,Form("#bf{%.0f#% scale uncertainty p-Pb}",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.55,0.42,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty (p-Pb)}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
  }

  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.025);
  tUncertainty->SetTextFont(62);
  //tUncertainty->Draw();

  

  // TString histoname = "DeltaPhi_";
  //histoname += Form("%d",i);
  //TH1D * outputhisto = (TH1D*)histo->Clone("blabla");
  //outputhisto->SetName(histoname.Data());
  // cout << "Scaling By Bin Width = " << outputhisto->GetBinWidth(3) << endl;
  // outputhisto->Scale(1./outputhisto->GetBinWidth(3));
    
  //  TH1D * histo2 = plotter->ReflectHistogram(histo);
    
  return hFDsub;
    
}


//_______________________________________________________________________
void SaveCanvas(TCanvas * c, TString directory, TString name){
    
  TString outputDir = "";//Plots/15_May/";
  gSystem->Exec(Form("mkdir %s",outputDir.Data()));
    
  if(directory != ""){outputDir += directory;
    TString exec = "mkdir -p ";
    exec += outputDir;
    cout << exec << endl;
    gSystem->Exec(exec.Data());
  }
    
    
  TString plotsout = c->GetName();//"Canvas_pT_05_";
  plotsout += name;
   
  c->SaveAs(Form("%s/%s.root",outputDir.Data(),plotsout.Data()));
  c->SaveAs(Form("%s/%s.png",outputDir.Data(),plotsout.Data()));
}

//_______________________________________________________________________
TH1D * subtractpedestal(Int_t i,TH1D *histo,Double_t pedestal){
  TString nameoutput = histo->GetName();
  nameoutput += "_subtr_";
  nameoutput += "pedestal";
  nameoutput += Form("%d",i);
  cout << "subtracting pedestal for histo " << i << " = " << pedestal << endl;
    
  TH1D * outputhisto = (TH1D*)histo->Clone(nameoutput.Data());
  outputhisto->Reset();
  Double_t value = 0;
    
  for(Int_t iBin = 1; iBin <= histo->GetNbinsX();iBin++){
        
        
    value = histo->GetBinContent(iBin);
    value -= pedestal;
        
    outputhisto->SetBinContent(iBin,value);
        
    outputhisto->SetBinError(iBin,histo->GetBinError(iBin));
        
  }
    
  //outputhisto->SetMarkerColor(color);
  //outputhisto->SetLineColor(color);
  //outputhisto->SetMarkerStyle(19+color);
    
  return outputhisto;
}

//_______________________________________________________________________
TH1D * ReflectHistogram(TH1D * h){
    
  TString   histname = h->GetName();
  histname += "_reflected";//Form("_reflected",i);
    
  TH1D *h2=new TH1D(histname.Data(),"; #Delta#varphi (rad); #frac{1}{N_{D}} #frac{dN}{d#Delta#varphi}^{ass} (rad^{-1})",h->GetNbinsX()/2.,0.,TMath::Pi());
  for(Int_t j=1;j<=h->GetNbinsX();j++){
    Double_t x=h->GetBinCenter(j);
    Double_t y0=h->GetBinContent(j);
    Double_t ey0=h->GetBinError(j);
    Int_t j2 = -1;
    if(x>0&&x<TMath::Pi())j2=h2->FindBin(x);
    else if(x<0)j2=h2->FindBin(-1.*x);
    else if(x>TMath::Pi())j2=h2->FindBin(2.*TMath::Pi()-x);
    else std::cout << "Point "<< j << " excluded " <<std::endl;
    Double_t y=h2->GetBinContent(j2);
    Double_t ey=h2->GetBinError(j2);
    h2->SetBinContent(j2,y+y0);
    h2->SetBinError(j2,TMath::Sqrt(ey0*ey0+ey*ey));
  }
    
  h2->GetYaxis()->SetRangeUser(0,1.4*h2->GetBinContent(h2->GetMaximumBin()));
  h2->SetLineColor(1);
  h2->SetMarkerColor(1);
  h2->SetMarkerStyle(20);
  return h2;
    
}
//_______________________________________________________________________
void DoComparison_ppVspPbTEST(TString pthad="0.3_1.0"){
    
  gStyle->SetOptStat(0000);
    
  // set all the paths for the averaged plots

  Init();
  LoadFileNamesppVspPb(pthad);       ///modified w proper files===========================
  histo = new TH1D *[nhistos];
  err = new TGraphAsymmErrors *[nhistos];
  ltscale = new TLatex *[nhistos];        
        
  // loop on histos
  for(Int_t k=0; k<nhistos; k++){ //GetHisto has to get graph/syst uncer===========================
    TLatex* lt;
    if(skip3to5){
      if(k%3==0)continue;
    }
    histo[k] = GetHistoAndSyst(k,"fhDaverage","AverageSystematicUncertainty",err[k],ltscale[k]); 
             
    if(k<3) { 
      histo[k]->SetMarkerColor(colpPb); 
      histo[k]->SetMarkerStyle(21); 
      histo[k]->SetLineColor(colpPb); 

	 
    } // get pp

    // if(k>=3) { 
    if(k>3) { //removed here case k=3 which corresponds to pPb 3-5, removed in the analysis
      //histo[k] = GetHisto(k,"cDraw","fhDaverage"); 
      histo[k]->SetMarkerColor(colpp); 
      histo[k]->SetLineColor(colpp); 
      histo[k]->SetMarkerStyle(20); 
    } // get p-pPb
    // 
    // histo[k]->Rebin(rebin); histo[k]->Scale(1./rebin);
  }

  //  return;
  Double_t BaselinePP, BaselineErrPP;// = 0;
  Double_t BaselinepPb, BaselineErrpPb;// = 0;

  TH1D ** subtractedhisto = new TH1D *[nhistos];
  TGraphAsymmErrors **suberr = new TGraphAsymmErrors *[nhistos];
  TGraphAsymmErrors **grbase=new TGraphAsymmErrors*[nhistos];
  TGraphAsymmErrors **grv2=new TGraphAsymmErrors*[nhistos]; 

  for(Int_t k=0; k<nhistos; k++){
    if(skip3to5){
      if(k%3==0)continue;
    }
    if(k<3) {
      subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pp",k,histo[k],err[k],suberr[k],"CanvasBaselineVariationTrendPedestal",grbase[k],grv2[k]);
      //      subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pp",k,histo[k],err[k],&gr,"CanvasBaselineVariationTrendPedestal",&b1,&b2);
      //    subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pPb",k,histo[k],err[k],&gr,"CanvasBaselineVariationTrendPedestal",&b1,&b2);
      //      cout<<"BOX1 = "<<b1->GetX1()<<"  "<<b1->GetX2()<<"  "<<b1->GetY1()<<"  "<<b1->GetY2()<<"  "<<endl;
      //      cout<<"BOX2 = "<<b2->GetX1()<<"  "<<b2->GetX2()<<"  "<<b2->GetY1()<<"  "<<b2->GetY2()<<"  "<<endl;
      cout<<"sub -> "<<subtractedhisto[k]->GetBinContent(5)<<endl;
      subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.7,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));    
      grbase[k]->SetFillStyle(3002);
      grbase[k]->SetFillColor(kBlue-7);
      grbase[k]->SetLineColor(kBlue-7);

      if(grv2[k]){
	grv2[k]->SetFillStyle(3002);
	grv2[k]->SetFillColor(kMagenta);  
      }
      suberr[k]->SetLineColor(kBlue);
    }

    if(k>=3) { 
      subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pPb",k,histo[k],err[k],suberr[k],"CanvasBaselineVariationTrendPedestal",grbase[k],grv2[k]);
      //      subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pp",k,histo[k],err[k],&gr,"CanvasBaselineVariationTrendPedestal",&b1,&b2);
      //    subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pPb",k,histo[k],err[k],&gr,"CanvasBaselineVariationTrendPedestal",&b1,&b2);
      //      cout<<"BOX1 = "<<b1->GetX1()<<"  "<<b1->GetX2()<<"  "<<b1->GetY1()<<"  "<<b1->GetY2()<<"  "<<endl;
      //      cout<<"BOX2 = "<<b2->GetX1()<<"  "<<b2->GetX2()<<"  "<<b2->GetY1()<<"  "<<b2->GetY2()<<"  "<<endl;
      cout<<"sub -> "<<subtractedhisto[k]->GetBinContent(5)<<endl;
      subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.8,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));    
      grbase[k]->SetFillStyle(3002);
      grbase[k]->SetFillColor(kRed-7);
      grbase[k]->SetLineColor(kRed-7);

      if(grv2[k]){
	grv2[k]->SetFillStyle(3002);
	grv2[k]->SetFillColor(kMagenta);  
      }
      suberr[k]->SetLineColor(kRed);
    } 
  }
  subtractedhisto[5]->SetMaximum(3.7);

  if(pthad.Contains("0.3")){
    if(pthad.Contains("1.0")){
      subtractedhisto[1]->SetMinimum(-0.9);
      subtractedhisto[2]->SetMinimum(-0.9);
      subtractedhisto[4]->SetMinimum(-0.9);
      subtractedhisto[5]->SetMinimum(-0.9);
      subtractedhisto[4]->SetMaximum(3);
      subtractedhisto[5]->SetMaximum(4);///QUI
    }
    else { // >0.3 case
      subtractedhisto[1]->SetMinimum(-0.9);
      subtractedhisto[2]->SetMinimum(-0.9);
      subtractedhisto[4]->SetMinimum(-0.9);
      subtractedhisto[5]->SetMinimum(-0.9);
      subtractedhisto[4]->SetMaximum(5);
      subtractedhisto[5]->SetMaximum(5);///QUI
    }
  }
  else if(pthad.Contains("1.0")){
    subtractedhisto[4]->SetMinimum(-0.5);
    subtractedhisto[5]->SetMinimum(-0.5);
    subtractedhisto[4]->SetMaximum(2.5);
    subtractedhisto[5]->SetMaximum(3);
 
  }
 
    
  TLegend * legend = new TLegend(0.21,0.63,0.47,0.72);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.027);
  legend->AddEntry(histo[ihskip],"pp #sqrt{#it{s}}=7 TeV, |#it{y}^{D}_{cms}|<0.5","lep");
  legend->AddEntry(histo[4],"p-Pb #sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04","lep");
  
  TLegend * legend2 = new TLegend(0.55,0.50,0.74,0.60);
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetBorderSize(0);
  legend2->SetTextSize(0.025);
  legend2->AddEntry(grbase[1],"baseline uncertainty pp","f");
  legend2->AddEntry(grbase[4],"baseline uncertainty p-Pb","f");
  // legend2->AddEntry(box2[4],"v2 subtr p-Pb","f");
  

  TPaveText *alice = new TPaveText(0.72,0.78,0.85,0.83,"NDC");
  SetPaveStyle(alice);
  alice->SetTextSize(0.04);
  alice->AddText("ALICE");//commented

  TPaveText *fitvalueslow = new TPaveText(0.21,0.835,0.75,0.873,"NDC"); // not used
  SetPaveStyle(fitvalueslow);
  //  TPaveText *fitvalueslow = new TPaveText(0.54,0.72,0.89,0.82,"NDC");
  fitvalueslow->SetTextSize(0.032);
  fitvalueslow->AddText("D meson - charged particle correlation");

  // fitvalueslow->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvalueslow1 = new TPaveText(0.21,0.83,0.5,0.863,"NDC");
  SetPaveStyle(fitvalueslow1);
  fitvalueslow1->SetTextSize(0.032);
  fitvalueslow1->AddText("Average D^{0},D^{+},D^{*+}");
  //  TPaveText *fitvalueslow2 = new TPaveText(0.13,0.774,0.81,0.81,"NDC");
  //  SetPaveStyle(fitvalueslow2);
  //  fitvalueslow2->AddText(Form("3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} > %s GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad.Data()));
    

  TPaveText *fitvaluesmid = new TPaveText(0.136,0.85,0.88,0.893,"NDC");// not used
  SetPaveStyle(fitvaluesmid);
  fitvaluesmid->SetTextFont(42);
  fitvaluesmid->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  
  TPaveText *fitvaluesmid2 = new TPaveText(0.21,0.73,0.64,0.82,"NDC");
  SetPaveStyle(fitvaluesmid2);
  //  TPaveText *fitvaluesmid2 = new TPaveText(0.15,0.78,0.84,0.825,"NDC");
  // TPaveText *fitvaluesmid2 = new TPaveText(0.18,0.79,0.84,0.843,"NDC");
  fitvaluesmid2->SetTextSize(0.03);
  fitvaluesmid2->AddText("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}");
  if(pthad.Contains("0.3")){
    if(pthad.Contains("1.0"))fitvaluesmid2->AddText("0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
    else fitvaluesmid2->AddText("#it{p}_{T}^{assoc} >0.3 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  }
  else if(pthad.Contains("1.0"))  fitvaluesmid2->AddText("#it{p}_{T}^{assoc} >1 GeV/#it{c}, |#Delta#eta| < 1.0 ");


  TPaveText *fitvalueshigh = new TPaveText(0.136,0.85,0.88,0.893,"NDC"); // not used
  SetPaveStyle(fitvalueshigh);
  fitvalueshigh->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvalueshigh2 = new TPaveText(0.21,0.73,0.64,0.82,"NDC");
  SetPaveStyle(fitvalueshigh2);
 //TPaveText *fitvalueshigh2 = new TPaveText(0.15,0.78,0.84,0.825,"NDC");
  //  TPaveText *fitvalueshigh2 = new TPaveText(0.18,0.79,0.84,0.843,"NDC");
  fitvalueshigh2->SetTextSize(0.03);
  fitvalueshigh2->AddText("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, |#it{y}^{D}|<0.5");
  //  fitvalueshigh2->AddText(Form("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} > %s GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad.Data()));
  if(pthad.Contains("0.3"))  {
    if(pthad.Contains("1.0"))fitvalueshigh2->AddText("0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
    else fitvalueshigh2->AddText("#it{p}_{T}^{assoc} >0.3 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  }
  else if(pthad.Contains("1.0"))  {
    fitvalueshigh2->AddText("#it{p}_{T}^{assoc} >1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  }
 
  Double_t startline,endline,wi;
  wi=subtractedhisto[4]->GetBinWidth(1);
  startline=subtractedhisto[4]->GetBinLowEdge(1);
  endline=subtractedhisto[4]->GetBinLowEdge(subtractedhisto[4]->GetNbinsX())+wi+0.4;
  // cout<<"******************* "<<startline<<"  "<<endline<<endl;
  TLine* line=new TLine(startline, 0, endline, 0);
  line->SetLineStyle(2);

  TCanvas * ptmidsub = GetCanvasWithStyle("ptmidsub","mid pt subtracted correlations");
  ptmidsub->cd();
  
  TH1D* h=new TH1D(*subtractedhisto[4]);
  h->Reset();
  h->GetXaxis()->SetLimits(startline,endline);
  h->SetLineColor(0);
  h->GetXaxis()->Delete();
  h->GetYaxis()->Delete();
  h->Draw();
  
  //  box2[4]->Draw("same");
  grbase[4]->Draw("E2");
  grbase[1]->Draw("E2"); 
  line->Draw(); 
  subtractedhisto[4]->Draw("same");
  suberr[4]->Draw("E2");
  ltscale[4]->Draw();//commented
  subtractedhisto[1]->Draw("same");
  suberr[1]->Draw("E2");
  ltscale[1]->Draw();//commented
  legend->Draw();
  legend2->Draw();
  //  fitvalueslow->Draw();
  //  fitvaluesmid->Draw();
  alice->Draw();
  fitvaluesmid2->Draw("same");
  fitvalueslow1->Draw("same");


  TCanvas * pthighsub = GetCanvasWithStyle("pthighsub","high pt subtracted correlations");
  pthighsub->cd();
  TH1D* h2=new TH1D(*subtractedhisto[5]);
  h2->Reset();
  h2->GetXaxis()->SetLimits(startline,endline);
  h2->SetLineColor(0);
  h2->GetXaxis()->Delete();
  h2->GetYaxis()->Delete();
  h2->Draw();

  subtractedhisto[5]->Draw("epsame");
  //box2[5]->Draw("same");
  grbase[5]->Draw("E2");
  grbase[2]->Draw("E2"); 
  line->Draw(); 
  subtractedhisto[5]->Draw("same");
  suberr[5]->Draw("E2");
  ltscale[5]->Draw();//commented
  subtractedhisto[2]->Draw("same");
  subtractedhisto[2]->Draw("same");
  suberr[2]->Draw("E2");
  ltscale[2]->Draw();//commented
  legend->Draw();
  legend2->Draw();
  //  fitvalueslow->Draw();
  // fitvalueshigh->Draw();
  alice->Draw();
  fitvalueshigh2->Draw("same");
  fitvalueslow1->Draw("same");

  cout<<pthighsub->GetFrame()->GetX1()<<"  "<<pthighsub->GetFrame()->GetY1()<<endl;

  Printf("Now saving canvases");

  ptmidsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_58_%s.root",avType.Data(),pthad.Data()));
  ptmidsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_58_%s.pdf",avType.Data(),pthad.Data()));
  ptmidsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_58_%s.eps",avType.Data(),pthad.Data()));
  //ptmidsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_58_%s.gif",avType.Data(),pthad.Data()));

  pthighsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_816_%s.root",avType.Data(),pthad.Data()));
  pthighsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_816_%s.pdf",avType.Data(),pthad.Data()));
  pthighsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_816_%s.eps",avType.Data(),pthad.Data()));
  //  pthighsub->SaveAs(Form("plotComparison_%sAverage_pp_pPb_816_%s.gif",avType.Data(),pthad.Data()));

 
}


void LoadFileNamesppVspPb(TString pthadron){
  //load pp
  filenames[0] = Form("%s/%sAverageppDzeroDstarDplus3to5_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  filenames[1] = Form("%s/%sAverageppDzeroDstarDplus5to8_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  filenames[2] = Form("%s/%sAverageppDzeroDstarDplus8to16_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  //load pPb
  filenames[3] = Form("%s/%sAveragepPbDzeroDstarDplus3to5_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  filenames[4] = Form("%s/%sAveragepPbDzeroDstarDplus5to8_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  filenames[5] = Form("%s/%sAveragepPbDzeroDstarDplus8to16_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());

    //load pedestals from fit (STILL TEMPORARY)
    for(Int_t k = 0; k<2; k++){
      pedestalfilenames[k] = baselinedirectory;
    }
    //    pedestalfilenames[0] += Form("/finalCanvas_DAverage_Baseline%s.root",pthadron.Data());//pp
    //    pedestalfilenames[1] += Form("/p_Pb_finalCanvas_DAverage_Baseline%s.root",pthadron.Data());//pPb
    pedestalfilenames[0] += Form("/Trends_pp/CanvasBaselineVariationTrendPedestal_pthad%s.root",pthadron.Data());//pp
    pedestalfilenames[1] += Form("/Trends_pPb/CanvasBaselineVariationTrendPedestal_pthad%s.root",pthadron.Data());//pPb

}

void GetBaseline(Double_t &av=0.,Double_t &errAv=0.,TH1D * h){
    
  Int_t nbins = 0;
  Double_t sum = 0;
  for(Int_t binPhi =0; binPhi<h->GetNbinsX();binPhi++){
        
    if(h->GetBinLowEdge(binPhi)>=-0.5*TMath::Pi() && h->GetBinLowEdge(binPhi+1)<=-0.25*TMath::Pi()){
      cout << "iBin = " << binPhi << endl;
      av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
      errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));
    }
        
    if(h->GetBinLowEdge(binPhi)>=0.25*TMath::Pi() && h->GetBinLowEdge(binPhi+1)<=0.5*TMath::Pi()){
      cout << "iBin = " << binPhi << endl;
      av+=h->GetBinContent(binPhi)/(h->GetBinError(binPhi)*h->GetBinError(binPhi));
      errAv+=1./(h->GetBinError(binPhi)*h->GetBinError(binPhi));
    }
  }
    
  av/=errAv;
  errAv=TMath::Sqrt(1./errAv);
  //  av/=2;
  //  errAv/=2;
  printf("Average baseline: %f +- %f \n",av,errAv);
  //fitFunction->FixParameter(0,av);
  //baseline=av;
  //errbaseline=errAv;
    
}

//_______________________________________________________________________
TH1D * GetPedestalHistoAndSystAndSubtractPedpPb(TString system="pPb",Int_t i,TH1D *histo, TGraphAsymmErrors* gr,TGraphAsymmErrors *&grout, TString canvasname,TGraphAsymmErrors *&grbaseOut,TGraphAsymmErrors *&grv2Out){
  //i=D pt bin
  //Double_t* arrbox1, Double_t* arrbox2
  Double_t value = 0, pedestal=0;
  TGraphAsymmErrors *grBase,*grV2;
  grbaseOut=new TGraphAsymmErrors();
  grbaseOut->SetName(Form("grbaselineUncFull_%sBin%d",system.Data(),i));
  grv2Out=new TGraphAsymmErrors();
  grv2Out->SetName(Form("grbaselineUncV2_%sBin%d",system.Data(),i));

  Double_t xuncFull,errxuncFull;
  Double_t xuncv2,errxuncv2;
  Int_t bin,bingr;
  if(system.Contains("pPb")){
    // get pedestal from fit outputs
    TString path = pedestalfilenames[1];//pPb ----------------
    cout << "pPb -->  Reading File from path: " << path << endl;
    
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    grBase=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fBaselineVariationSystematicsPedestal");    
    grV2=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fv2SystematicsPedestal");    
    //example:
    //    TBox  X1= 5.600000 Y1=6.724682 X2=7.400000 Y2=6.802155     bin 5-8 : v2
    //    TBox  X1= 11.100000 Y1=6.563838 X2=12.900000 Y2=6.638230   bin 8-16 :v2
    //    TBox  X1= 5.800000 Y1=6.146429 X2=7.200000 Y2=7.380408     bin 5-8 : Total (no v2 modulation)
    //    TBox  X1= 11.300000 Y1=5.990975 X2=12.700000 Y2=7.211093   bin 8-16 : Total (no v2 modulation)
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
    if(i==3){
      bin=h->FindBin(4.);
      bingr=GetBinGraph(4.,grBase);
    }
    else if(i==4){
      bin=h->FindBin(6.5);
      bingr=GetBinGraph(6.5,grBase);
    }
    else if(i==5){
      bin=h->FindBin(12.);
      bingr=GetBinGraph(12.,grBase);
    }
      
    pedestal=h->GetBinContent(bin);
    Double_t x,y,erryl,erryh;
    grBase->GetPoint(bingr,x,y);
    Printf("histo: x=%f, graph: %f",h->GetBinCenter(bin),x);
    erryl=grBase->GetErrorYlow(bingr);
    erryh=grBase->GetErrorYhigh(bingr);
    
    grbaseOut->SetPoint(0,xuncFull,0);
    grbaseOut->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    if(grV2){
      grV2->GetPoint(bingr,x,y);
      erryl=grV2->GetErrorYlow(bingr);
      erryh=grV2->GetErrorYhigh(bingr);	
      grv2Out->SetPoint(0,xuncFull,0);
      grv2Out->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    }
    
    cout<<i<<"  pedestal pPb="<<pedestal<<endl;

   
   
    cout << "pPb --> subtracting pedestal for histo " << i << " = " << pedestal << endl;
    //    cout<<"GR1 = "<<bout1->GetX1()<<"  "<<bout1->GetX2()<<"  "<<bout1->GetY1()<<"  "<<bout1->GetY2()<<"  "<<endl;
    //    if(grv2){
    //      cout<<"GR2 = "<<bout2->GetX1()<<"  "<<bout2->GetX2()<<"  "<<bout2->GetY1()<<"  "<<bout2->GetY2()<<"  "<<endl;
    //    }
    //    else{
    //      Printf("v2 uncertainty on baseline graph not present, won't be used");
    //    }
  }
    
  if(system.Contains("pp")){
    // get pedestal from fit outputs
    TString path = pedestalfilenames[0];//pp ----------------
    grV2=0x0;
    cout << "pp -->  Reading File from path: " << path << endl;
    
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    grBase=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fBaselineVariationSystematicsPedestal");    

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

    if(i==0){
      bin=h->FindBin(4.);
      bingr=GetBinGraph(4.,grBase);
    }
    else if(i==1){
      bin=h->FindBin(6.5);
      bingr=GetBinGraph(6.5,grBase);
    }
    else if(i==2){
      bin=h->FindBin(12.);
      bingr=GetBinGraph(12.,grBase);
    }
    pedestal=h->GetBinContent(bin);
    Double_t x,y,erryl,erryh;
    Printf("histo: x=%f, graph: %f",h->GetBinCenter(bin),x);
    grBase->GetPoint(bingr,x,y);
    erryl=grBase->GetErrorYlow(bingr);
    erryh=grBase->GetErrorYhigh(bingr);
    
    grbaseOut->SetPoint(0,xuncFull,0);
    grbaseOut->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    if(grV2){
      grV2->GetPoint(bingr,x,y);
      erryl=grV2->GetErrorYlow(bingr);
      erryh=grV2->GetErrorYhigh(bingr);	
      grv2Out->SetPoint(0,xuncFull,0);
      grv2Out->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
    }
    cout<<i<<"  pedestal pp="<<pedestal<<endl;


    cout << "pp --> subtracting pedestal for histo " << i << " = " << pedestal << endl;
    //    cout<<"BOX1 = "<<bout1->GetX1()<<"  "<<bout1->GetX2()<<"  "<<bout1->GetY1()<<"  "<<bout1->GetY2()<<"  "<<endl;
    //    cout<<"BOX2 = "<<bout2->GetX1()<<"  "<<bout2->GetX2()<<"  "<<bout2->GetY1()<<"  "<<bout2->GetY2()<<"  "<<endl;
  }
  
  grout=(TGraphAsymmErrors*)gr->Clone(Form("grSub_%sBin%d",system.Data(),i));
  
  TString nameoutput = histo->GetName();
  nameoutput += "_subtr_";
  nameoutput += "pedestal";
  nameoutput += Form("%d",i);
   
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
  cout<<"now return histogram"<<endl;
  
  outputhisto->SetXTitle("#Delta#varphi (D,h) (rad)");
  outputhisto->SetYTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} - baseline (rad^{-1})");
  outputhisto->GetYaxis()->SetTitleOffset(1.5);
  outputhisto->GetYaxis()->SetTitleSize(0.04);
  outputhisto->GetXaxis()->SetTitleSize(0.04);
  outputhisto->GetYaxis()->SetLabelSize(0.04);
  outputhisto->GetXaxis()->SetLabelSize(0.04);

  return outputhisto;
    
}
