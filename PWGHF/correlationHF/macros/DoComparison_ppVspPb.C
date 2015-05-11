//############################################################################
// 
//  Macros to perform and draw the comparison of D-haron correlation results
//  in pp and p-Pb collisions
//  elena.bruna@to.infn.it
//
//############################################################################
Int_t nhistos=6;
Int_t rebin =1;
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
  hFDsub->SetXTitle("#Delta#varphi (rad)");
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
TH1D * GetPedestalHistoAndSystAndSubtractPedpPb(TString system="pPb",Int_t i,TH1D *histo, TGraphAsymmErrors* gr,TGraphAsymmErrors** grout, TString canvasname,TBox** bsyst, TBox** bsyst2){
  //i=D pt bin
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
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->At(1);
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
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->At(1);
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
  
  outputhisto->SetXTitle("#Delta#varphi (rad)");
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
TH1D * GetHistoAndSyst(Int_t i, TString hname, TString hnamesyst,TGraphAsymmErrors** gr2, TLatex** lt){
    
   
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

  TGraphAsymmErrors *gr=syst->GetTotNonFlatUncGraph();
  gr->SetLineColor(kBlack);
  gr->SetMarkerColor(kBlack);
  gr->SetFillStyle(0);

  (*gr2)=gr;

  TLatex *tUncertainty;
  if(i<3){
    if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.2,0.6,Form("#bf{%.0f#% scale uncertainty pp}",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.55,0.59,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty pp}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
   
  }
if(i>=3){
    if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.2,0.6,Form("#bf{%.0f#% scale uncertainty p-Pb}",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.55,0.55,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty p-Pb}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
  }

  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.025);
  //tUncertainty->Draw();
  (*lt)=tUncertainty;
  

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
    
  TH1D *h2=new TH1D(histname.Data(),"; #Delta#varphi(D,had) (rad); #frac{1}{N_{D}} #frac{dN}{d#Delta#varphi}^{ass} (rad^{-1})",h->GetNbinsX()/2.,0.,TMath::Pi());
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
    TGraphAsymmErrors* gr;
    TLatex* lt;
    if(skip3to5){
      if(k%3==0)continue;
    }
    histo[k] = GetHistoAndSyst(k,"fhDaverage","AverageSystematicUncertainty",&gr,&lt); 
       
    err[k]=gr;
    ltscale[k]=lt;      
      
    if(k<3) { 
      histo[k]->SetMarkerColor(4); 
      histo[k]->SetMarkerStyle(21); 
      histo[k]->SetLineColor(4); 

	 
    } // get pp

    // if(k>=3) { 
    if(k>3) { //removed here case k=3 which corresponds to pPb 3-5, removed in the analysis
      //histo[k] = GetHisto(k,"cDraw","fhDaverage"); 
      histo[k]->SetMarkerColor(2); 
      histo[k]->SetLineColor(2); 
      histo[k]->SetMarkerStyle(20); 
    } // get p-pPb
    // 
    // histo[k]->Rebin(rebin); histo[k]->Scale(1./rebin);
  }

  //  return;
  Double_t BaselinePP, BaselineErrPP;// = 0;
  Double_t BaselinepPb, BaselineErrpPb;// = 0;

  TH1D ** subtractedhisto = new TH1D *[nhistos];
  TGraphAsymmErrors ** suberr = new TGraphAsymmErrors *[nhistos];
  TBox** box1=new TBox*[nhistos];
  TBox** box2=new TBox*[nhistos];
 

  for(Int_t k=0; k<nhistos; k++){
    if(skip3to5){
      if(k%3==0)continue;
    }
    if(k<3) {
      TGraphAsymmErrors* gr;
      TBox* b1;
      TBox* b2;
      subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pp",k,histo[k],err[k],&gr,"finalCanvas",&b1,&b2);
      //    subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pPb",k,histo[k],err[k],&gr,"finalCanvas",&b1,&b2);
      suberr[k]=gr;
      box1[k]=b1;
      box2[k]=b2;
      cout<<"BOX1 = "<<b1->GetX1()<<"  "<<b1->GetX2()<<"  "<<b1->GetY1()<<"  "<<b1->GetY2()<<"  "<<endl;
      cout<<"BOX2 = "<<b2->GetX1()<<"  "<<b2->GetX2()<<"  "<<b2->GetY1()<<"  "<<b2->GetY2()<<"  "<<endl;
      cout<<"sub -> "<<subtractedhisto[k]->GetBinContent(5)<<endl;
      subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.7,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));    
      box1[k]->SetFillStyle(3002);
      box1[k]->SetFillColor(kBlue-7);
      box1[k]->SetLineColor(kBlue-7);

      box2[k]->SetFillStyle(3002);
      box2[k]->SetFillColor(kMagenta);  
      suberr[k]->SetLineColor(kBlue);
    }

    if(k>3) { //removed k=3 (corresponding to pPb 3-5)
      //if(k>=3) {
      TGraphAsymmErrors* gr;
      TBox* b1;
      TBox* b2;
      //subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pp",k,histo[k],err[k],&gr,"finalCanvas",&b1,&b2);
      subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pPb",k,histo[k],err[k],&gr,"finalCanvas",&b1,&b2);
      suberr[k]=gr;
      box1[k]=b1;
      box2[k]=b2;
      cout<<"BOX1 = "<<b1->GetX1()<<"  "<<b1->GetX2()<<"  "<<b1->GetY1()<<"  "<<b1->GetY2()<<"  "<<endl;
      cout<<"BOX2 = "<<b2->GetX1()<<"  "<<b2->GetX2()<<"  "<<b2->GetY1()<<"  "<<b2->GetY2()<<"  "<<endl;
      cout<<"sub -> "<<subtractedhisto[k]->GetBinContent(5)<<endl;
      subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.8,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));    
      box1[k]->SetFillStyle(3002);
      box1[k]->SetFillColor(kRed-7);
      box1[k]->SetLineColor(kRed-7);
  
      box2[k]->SetFillStyle(3002);
      box2[k]->SetFillColor(kMagenta); 
      suberr[k]->SetLineColor(kRed); 
    }
    
  }
  subtractedhisto[5]->SetMaximum(3.7);
 
  if(pthad.Contains("0.3")){
    subtractedhisto[1]->SetMinimum(-0.9);
    subtractedhisto[2]->SetMinimum(-0.9);
    subtractedhisto[4]->SetMinimum(-0.9);
    subtractedhisto[5]->SetMinimum(-0.9);
    subtractedhisto[4]->SetMaximum(4);
    subtractedhisto[5]->SetMaximum(4);///QUI
  }
  else if(pthad.Contains("1.0")){
    subtractedhisto[4]->SetMinimum(-0.9);
    subtractedhisto[5]->SetMinimum(-0.9);
    subtractedhisto[4]->SetMaximum(4);
    subtractedhisto[5]->SetMaximum(4);
 
  }
 
    
  TLegend * legend = new TLegend(0.14,0.67,0.47,0.77);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.027);
  legend->AddEntry(histo[ihskip],"pp #sqrt{#it{s}}=7 TeV","lep");
  legend->AddEntry(histo[4],"p-Pb #sqrt{#it{s}_{NN}}=5.02 TeV","lep");
  
  TLegend * legend2 = new TLegend(0.55,0.62,0.74,0.72);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);
  legend2->SetTextSize(0.025);
  legend2->AddEntry(box1[1],"baseline uncertainty pp","f");
  legend2->AddEntry(box1[4],"baseline uncertainty p-Pb","f");
  // legend2->AddEntry(box2[4],"v2 subtr p-Pb","f");
  

  TPaveText *alice = new TPaveText(0.58,0.73,0.82,0.76,"NDC");
  alice->SetBorderSize(0);
  alice->SetFillColor(0);
  alice->SetTextFont(42);
  alice->SetTextSize(0.04);
  //  alice->AddText("ALICE Preliminary");//commented

  TPaveText *fitvalueslow = new TPaveText(0.17,0.855,0.7,0.893,"NDC");
  //  TPaveText *fitvalueslow = new TPaveText(0.54,0.72,0.89,0.82,"NDC");
  fitvalueslow->SetBorderSize(0);
  fitvalueslow->SetFillColor(0);
  fitvalueslow->SetTextFont(42);
  fitvalueslow->AddText("D meson  - charged particle correlation");
  // fitvalueslow->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvalueslow1 = new TPaveText(0.134,0.82,0.45,0.855,"NDC");
  fitvalueslow1->SetBorderSize(0);
  fitvalueslow1->SetFillColor(0);
  fitvalueslow1->SetTextFont(42);
  fitvalueslow1->AddText("Average D^{0},D^{+},D^{*+}");
  TPaveText *fitvalueslow2 = new TPaveText(0.13,0.774,0.81,0.81,"NDC");
  fitvalueslow2->SetBorderSize(0);
  fitvalueslow2->SetFillColor(0);
  fitvalueslow2->SetTextFont(42);
  fitvalueslow2->AddText(Form("3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} > %s GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad.Data()));
    

TPaveText *fitvaluesmid = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  fitvaluesmid->SetBorderSize(0);
  fitvaluesmid->SetFillColor(0);
  fitvaluesmid->SetTextFont(42);
  fitvaluesmid->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvaluesmid2 = new TPaveText(0.15,0.774,0.84,0.81,"NDC");
  //  TPaveText *fitvaluesmid2 = new TPaveText(0.15,0.78,0.84,0.825,"NDC");
  // TPaveText *fitvaluesmid2 = new TPaveText(0.18,0.79,0.84,0.843,"NDC");
  fitvaluesmid2->SetBorderSize(0);
  fitvaluesmid2->SetFillColor(0);
  fitvaluesmid2->SetTextFont(42);
  if(pthad.Contains("0.3"))  fitvaluesmid2->AddText("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  else if(pthad.Contains("1.0"))  fitvaluesmid2->AddText("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} >1 GeV/#it{c}, |#Delta#eta| < 1.0 ");


 TPaveText *fitvalueshigh = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  fitvalueshigh->SetBorderSize(0);
  fitvalueshigh->SetFillColor(0);
  fitvalueshigh->SetTextFont(42);
  fitvalueshigh->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvalueshigh2 = new TPaveText(0.15,0.774,0.84,0.81,"NDC");
 //TPaveText *fitvalueshigh2 = new TPaveText(0.15,0.78,0.84,0.825,"NDC");
  //  TPaveText *fitvalueshigh2 = new TPaveText(0.18,0.79,0.84,0.843,"NDC");
  fitvalueshigh2->SetBorderSize(0);
  fitvalueshigh2->SetFillColor(0);
  fitvalueshigh2->SetTextFont(42);
  //  fitvalueshigh2->AddText(Form("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} > %s GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad.Data()));
  if(pthad.Contains("0.3"))  fitvalueshigh2->AddText("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  else if(pthad.Contains("1.0"))  fitvalueshigh2->AddText("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} >1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
 
  Double_t startline,endline,wi;
  wi=subtractedhisto[4]->GetBinWidth(1);
  startline=subtractedhisto[4]->GetBinLowEdge(1);
  endline=subtractedhisto[4]->GetBinLowEdge(subtractedhisto[4]->GetNbinsX())+wi+0.5;
  // cout<<"******************* "<<startline<<"  "<<endline<<endl;
  TLine* line=new TLine(startline, 0, endline, 0);
  line->SetLineStyle(2);

  TCanvas * ptmidsub = new TCanvas("cc","mid pt subtracted correlations",800,800);
  ptmidsub->cd();
  ptmidsub->SetTicky();
  ptmidsub->SetTickx();

  ptmidsub->Range(-2.6,-1.5,5.4,4.4);
  ptmidsub->SetLeftMargin(0.125);
  
  TH1D* h=new TH1D(*subtractedhisto[4]);
  h->Reset();
  h->GetXaxis()->SetLimits(startline,endline);
  h->SetLineColor(0);
  h->GetXaxis()->Delete();
  h->GetYaxis()->Delete();
  h->Draw();
  
  //  box2[4]->Draw("same");
  box1[4]->Draw("lsame");
  box1[1]->Draw("lsame"); 
  line->Draw(); 
  subtractedhisto[4]->Draw("same");
  suberr[4]->Draw("E2");
  ltscale[4]->Draw();//commented
  subtractedhisto[1]->Draw("same");
  subtractedhisto[1]->Draw("same");
  suberr[1]->Draw("E2");
  ltscale[1]->Draw();//commented
  legend->Draw();
  legend2->Draw();
  fitvalueslow->Draw();
  //  fitvaluesmid->Draw();
  alice->Draw();
  fitvaluesmid2->Draw("same");
  fitvalueslow1->Draw("same");


  TCanvas * pthighsub = new TCanvas("pthighsub","high pt subtracted correlations",416,90,800,800);
  pthighsub->cd();
  pthighsub->SetTicky();
  pthighsub->SetTickx();

  pthighsub->Range(-2.6,-1.5,5.4,4.4);
  pthighsub->SetLeftMargin(0.125);
  TH1D* h2=new TH1D(*subtractedhisto[5]);
  h2->Reset();
  h2->GetXaxis()->SetLimits(startline,endline);
  h2->SetLineColor(0);
  h2->GetXaxis()->Delete();
  h2->GetYaxis()->Delete();
  h2->Draw();

  subtractedhisto[5]->Draw("epsame");
  //box2[5]->Draw("same");
  box1[5]->Draw("lsame");
  box1[2]->Draw("lsame"); 
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
  fitvalueslow->Draw();
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
//_______________________________________________________________________
void DoComparison_ppVspPb(Double_t pthad=0.3){
 
    
  gStyle->SetOptStat(0);
    
  // set all the paths for the averaged plots
  for(Int_t k = 0; k<nhistos; k++){
    filenames[k] = directory;
  }
  for(Int_t k = 0; k<2; k++){
    pedestalfilenames[k] = directory;
  }
        
  LoadFileNamesppVspPb(pthad);       ///modified w proper files===========================

   
  TH1D ** histo = new TH1D *[nhistos];
  // loop on histos
  //   for(Int_t k=0; k<nhistos; k++){ //GetHisto has to get graph/syst uncer===========================

  //         if(k<3) { histo[k] = GetHisto(k,"cDraw","fhDaverage"); histo[k]->SetMarkerColor(2); histo[k]->SetMarkerStyle(20); } // get p-Pb
  //         if(k>=3) { histo[k] = GetHisto(k,"cDraw","fhDaverage"); histo[k]->SetMarkerColor(4); histo[k]->SetMarkerStyle(21); } // get p-p
  //           histo[k]->Rebin(rebin); histo[k]->Scale(1./rebin);
  //     }

    
    
  Double_t BaselinePP, BaselineErrPP;// = 0;
  Double_t BaselinepPb, BaselineErrpPb;// = 0;

  TH1D ** subtractedhisto = new TH1D *[nhistos];
  for(Int_t k=0; k<nhistos; k++){
    if(k<3) {//baseline with syst uncert?===========================
      GetBaselne(BaselinePP,BaselineErrPP, histo[k]);
      cout <<"------pp-------> " <<BaselinePP << " +/- " <<BaselineErrPP<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinePP);
        
    } // get p-Pb
    if(k>=3) {
      GetBaselne(BaselinepPb,BaselineErrpPb, histo[k]);
      cout <<"------pPb-------> " <<BaselinepPb << " +/- " <<BaselineErrpPb<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb);
        
    }  // subtract pp
        
    subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.5,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));
    //===========================subtract baseline also form point syst uncer.

  }

  // reflect histos
  TH1D ** reflectedhisto = new TH1D *[nhistos];
  cout << "test 3 "<< endl;
  for(Int_t k=0; k<nhistos; k++){
    reflectedhisto[k] = ReflectHistogram(subtractedhisto[k]);
    if(k<3) {reflectedhisto[k]->SetMarkerColor(2); reflectedhisto[k]->SetMarkerStyle(20);}
    if(k>=3) {reflectedhisto[k]->SetMarkerColor(4); reflectedhisto[k]->SetMarkerStyle(21);}
    reflectedhisto[k]->GetYaxis()->SetRangeUser(-0.5,2*reflectedhisto[k]->GetBinContent(reflectedhisto[k]->GetMaximumBin()));
  }
  cout << "test 4 "<< endl;
    

  TLegend *legend = new TLegend(0.1,0.8,0.5,0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.025);
  legend->AddEntry(histo[0+ihskip],"pp #sqrt(s) 7 TeV Data","lep");
  legend->AddEntry(histo[3+ihskip],"p-Pb #sqrt(s_{NN}) 5.02 TeV Data","lep");
    
    
  TPaveText *fitvalueslow = new TPaveText(0.51,0.78,0.85,0.88,"NDC");
  fitvalueslow->SetBorderSize(0);
  fitvalueslow->SetFillColor(0);
  fitvalueslow->AddText("3 GeV/c < p_{T}^{D} < 5 GeV/c");
  fitvalueslow->AddText(Form("p_{T}^{Ass} > %.1f GeV/c ",pthad));
  fitvalueslow->AddText(Form("|#Delta#eta| < %.1f ",1.0));
    
  TPaveText *fitvaluesmid = new TPaveText(0.51,0.78,0.85,0.88,"NDC");
  fitvaluesmid->SetBorderSize(0);
  fitvaluesmid->SetFillColor(0);
  fitvaluesmid->AddText("5 GeV/c < p_{T}^{D} < 8 GeV/c");
  fitvaluesmid->AddText(Form("p_{T}^{Ass} > %.1f GeV/c ",pthad));
  fitvaluesmid->AddText(Form("|#Delta#eta| < %.1f ",1.0));
    
  TPaveText *fitvalueshigh = new TPaveText(0.51,0.78,0.85,0.88,"NDC");
  fitvalueshigh->SetBorderSize(0);
  fitvalueshigh->SetFillColor(0);
  fitvalueshigh->AddText("8 GeV/c < p_{T}^{D} < 16 GeV/c");
  fitvalueshigh->AddText(Form("p_{T}^{Ass} > %.1f GeV/c ",pthad));
  fitvalueshigh->AddText(Form("|#Delta#eta| < %.1f ",1.0));
    
  TCanvas * ptlow = new TCanvas("Correlation_ptlow","pT 3-5 correlations",0,0,800,800);
  TCanvas * ptmid = new TCanvas("Correlation_ptmid","pT 5-8 correlations",0,0,800,800);
  TCanvas * pthigh = new TCanvas("Correlation_pthigh","pT 8-16 correlations",0,0,800,800);
    
  TCanvas * ptlowRef = new TCanvas("CorrelationReflected_ptlow","pT 3-5 correlations",0,0,800,800);
  TCanvas * ptmidRef = new TCanvas("CorrelationReflected_ptmid","pT 5-8 correlations",0,0,800,800);
  TCanvas * pthighRef = new TCanvas("CorrelationReflected_pthigh","pT 8-16 correlations",0,0,800,800);
    
  ptlow->cd();
  subtractedhisto[0+ihskip]->Draw("ep");
  subtractedhisto[3+ihskip]->Draw("sameep");
  legend->Draw("same");
  fitvalueslow->Draw("same");
    
  ptlowRef->cd();
  reflectedhisto[1]->Draw("ep");
  reflectedhisto[4]->Draw("sameep");
  legend->Draw("same");
  fitvalueslow->Draw("same");
    
    
  ptmid->cd();
  subtractedhisto[2]->Draw("ep");
  subtractedhisto[5]->Draw("sameep");
  legend->Draw("same");
  fitvaluesmid->Draw("same");
    
  ptmidRef->cd();
  reflectedhisto[0+ihskip]->Draw("ep");
  reflectedhisto[3+ihskip]->Draw("sameep");
  legend->Draw("same");
  fitvaluesmid->Draw("same");
    
    
  pthigh->cd();
  subtractedhisto[1]->Draw("ep");
  subtractedhisto[4]->Draw("sameep");
  legend->Draw("same");
  fitvalueshigh->Draw("same");
    
  pthighRef->cd();
  reflectedhisto[2]->Draw("ep");
  reflectedhisto[5]->Draw("sameep");
  legend->Draw("same");
  fitvalueshigh->Draw("same");
    
    
    
    
  // saving the canvases in .root and .png
  TString ptoutput="",ptoutputRef="", direcname="", direcnameRef="";
  ptoutput += Form("%sAverage_ppVspPb",avType.Data());
  direcname += "Output_Plots/ppVspPbPlots";
    
    
    
  // changes name of final canvas/histo
  if(pthad == 0.3) ptoutput += "_pt03";
  if(pthad == 0.5) ptoutput += "_pt05";
  if(pthad == 1.0) ptoutput += "_pt1";
    
    
  SaveCanvas(ptlow, direcname, ptoutput);
  SaveCanvas(ptmid, direcname, ptoutput);
  SaveCanvas(pthigh, direcname, ptoutput);
    
  //ptlow->Close();
  //ptmid->Close();
  //pthigh->Close();
    
    
  ptoutputRef += Form("%sAverage_ppVspPb",avType.Data());
  direcnameRef += "Output_PlotsRef/ppVspPbPlots";
    
  // changes name of final canvas/histo
  if(pthad == 0.3) ptoutputRef += "_Reflected_pt03";
  if(pthad == 0.5) ptoutputRef += "_Reflected_pt05";
  if(pthad == 1.0) ptoutputRef += "_Reflected_pt1";
    
  SaveCanvas(ptlowRef, direcnameRef, ptoutputRef);
  SaveCanvas(ptmidRef, direcnameRef, ptoutputRef);
  SaveCanvas(pthighRef, direcnameRef, ptoutputRef);
    
  //ptlowRef->Close();
  //ptmidRef->Close();
  //pthighRef->Close();
    
  return;
    
    
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
    pedestalfilenames[0] += Form("/finalCanvas_DAverage_Baseline%s.root",pthadron.Data());//pp
    pedestalfilenames[1] += Form("/p_Pb_finalCanvas_DAverage_Baseline%s.root",pthadron.Data());//pPb

}

void GetBaselne(Double_t &av=0.,Double_t &errAv=0.,TH1D * h){
    
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


void DoAllComparison(){
    
    
  DoComparison_ppVspPb(0.3);
  DoComparison_ppVspPb(0.5);
  DoComparison_ppVspPb(1.0);
  cout << "Done ! Check yours plots in folder ...  " << endl;
    
    
}


