//############################################################################
// 
//  Macros to perform and draw the comparison of D-haron correlation results
//  in pPb collisions and MC simulations
//  elena.bruna@to.infn.it
//
//############################################################################
Int_t nhistos=12;
Int_t rebin =1;
TString inputdatadirectory = "./";
TString inputtemplatedirecotry="./";
TString baselinedirectory="./";
TString avType="Weighted";
TString *filenames;
TString *pedestalfilenames;
Bool_t reflTempl=kFALSE;
Bool_t skip3to5=kFALSE;
TH1D ** histo;
TGraphAsymmErrors **err;
TLatex** ltscale;
TString fitplotmacrodir=gSystem->ExpandPathName("$ALICE_PHYSICS/../src/PWGHF/correlationHF/macros/");
Bool_t isReflectedData=kFALSE;
TString strSystemFDtempl="none";
void SetFDtemplateSystemString(TString str){
  strSystemFDtempl=str;
}

void SetIsDataReflected(Bool_t isrefl){
  isReflectedData=isrefl;
}
void SetFitPlotMacroPath(TString strdir){
  fitplotmacrodir=strdir;
}

void SetSkip3to5pPb(Bool_t skip){
  skip3to5=skip;
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
TH1D * GetHistoAndSyst(Int_t i, TString hname, TString hnamesyst,TGraphAsymmErrors *&gr, TLatex *&tUncertainty){
    
   
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

  TFile *f=TFile::Open(path.Data(),"READ");
  TH1D *hFDsub=(TH1D*)f->Get(hname.Data());
  AliHFDhadronCorrSystUnc *syst=(AliHFDhadronCorrSystUnc*)f->Get(hnamesyst.Data());
  TH1D *hUncCorrMin=syst->GetHistoTotFlatMin();
  TH1D *hUncCorrMax=syst->GetHistoTotFlatMax();

  gr=syst->GetTotNonFlatUncGraph();
  gr->SetLineColor(kBlack);
  gr->SetMarkerColor(kBlack);
  gr->SetFillStyle(0);


  if(i<3){

    if(TMath::Abs(hUncCorrMin->GetBinContent(1)-hUncCorrMax->GetBinContent(1))<0.001)tUncertainty=new TLatex(0.2,0.6,Form("#bf{%.0f#% scale uncertainty p-Pb}",hUncCorrMin->GetBinContent(1)*100.));
    else tUncertainty=new TLatex(0.63,0.61,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
   
    //  else tUncertainty=new TLatex(0.63,0.64,Form("#bf{{}^{+%.0f%s}_{-%.0f%s} scale uncertainty}","%","%",TMath::Abs(hUncCorrMax->GetBinContent(1))*100.,TMath::Abs(hUncCorrMin->GetBinContent(1)*100.)));
   
  }
 

  tUncertainty->SetNDC();
  tUncertainty->SetTextSize(0.0255);
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
TH1D * GetPedestalHisto(Int_t i, TString canvasname, TString hname){
    
  // get pedestal from fit outputs
  TString path = pedestalfilenames[i];
    
  cout << "Reading File from path: " << path << endl;
    
  TFile * file = TFile::Open(path.Data(),"READ");
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
TH1D * GetPedestalHistoAndSystAndSubtractPedpPb(Int_t i,TH1D *histo, TGraphAsymmErrors* gr,TGraphAsymmErrors** grout, TString canvasname,TBox** bsyst, TBox** bsyst2){
  //i=D pt bin
  //Double_t* arrbox1, Double_t* arrbox2
  TBox* bout1;
  TBox* bout2;

  Double_t value = 0, pedestal=0;
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
  //    TBox  X1= 11.300000 Y1=5.990975 X2=12.700000 Y2=7.211093   bin 8-18 : Total (no v2 modulation)
 
     
  if(i==0)pedestal=h->GetBinContent(1);
  if(i==1)pedestal=h->GetBinContent(2);
  if(i==2)pedestal=h->GetBinContent(3);
  bout1=(TBox*)b1->Clone("b1");//Total
  bout2=(TBox*)b1->Clone("b2");//v2

  if(i==0){
    bout2->SetX1(0);
    bout2->SetX2(0);
    bout2->SetY1(0);
    bout2->SetY2(0);

    bout1->SetX1(0);
    bout1->SetX2(0);
    bout1->SetY1(0);
    bout1->SetY2(0);
  }
  if(i==1){
    bout2->SetX1(-1.5);
    bout2->SetX2(-1.2);
    bout2->SetY1(b1->GetY1()-pedestal);
    bout2->SetY2(b1->GetY2()-pedestal);

    //    bout1->SetX1(-1.4);
    //     bout1->SetX2(-1.3);
    bout1->SetX1(4.75);
    bout1->SetX2(4.9);
    bout1->SetY1(b3->GetY1()-pedestal);
    bout1->SetY2(b3->GetY2()-pedestal);
  }

  if(i==2){
    bout2->SetX1(-1.5);
    bout2->SetX2(-1.2);
    bout2->SetY1(b2->GetY1()-pedestal);
    bout2->SetY2(b2->GetY2()-pedestal);

    bout1->SetX1(4.75);
    bout1->SetX2(4.9);
    //bout1->SetX1(-1.4);
    //bout1->SetX2(-1.3);
    bout1->SetY1(b4->GetY1()-pedestal);
    bout1->SetY2(b4->GetY2()-pedestal);
  }

  cout << "pPb --> subtracting pedestal for histo " << i << " = " << pedestal << endl;
  cout<<"BOX1 = "<<bout1->GetX1()<<"  "<<bout1->GetX2()<<"  "<<bout1->GetY1()<<"  "<<bout1->GetY2()<<"  "<<endl;
  cout<<"BOX2 = "<<bout2->GetX1()<<"  "<<bout2->GetX2()<<"  "<<bout2->GetY1()<<"  "<<bout2->GetY2()<<"  "<<endl;
 
  
  
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
  //  TCanvas * c = (TCanvas*)file->Get(canvasname.Data());
  //     c->cd();
  
  //     TH1D * histo = (TH1D*)c->FindObject(hname.Data());
    
  //     TString histoname = "Pedestal_";
  //     histoname += Form("%d",i);
  //     TH1D * outputhisto = (TH1D*)histo->Clone("blabla");
  //     outputhisto->SetName(histoname.Data());
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
  TH1D *histo=0x0;    
  TFile * file = TFile::Open(path.Data(),"WRITE");
  histo=(TH1D*)file->Get(hname.Data());
  if(!histo){
    //if(i ==2){cout << "file content = " << endl; file->ls();}
    TCanvas * c = (TCanvas*)file->Get(canvasname.Data());
    c->cd();
  // if(i ==2) {
  //   cout << "canvas content = " << endl; c->ls(); }
   histo = (TH1D*)c->FindObject(hname.Data());
  }
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
  c->SaveAs(Form("%s/%s.eps",outputDir.Data(),plotsout.Data()));
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
    
  TH1D *h2=new TH1D(histname.Data(),"; #Delta#phi(D,had) (rad); #frac{1}{N_{D}} #frac{dN}{d#Delta#phi}^{ass} (rad^{-1})",h->GetNbinsX()/2.,0.,TMath::Pi());
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
void DoComparison_pPbVsMCTEST(TString pthad = "0.3_1.0"){
 
  gROOT->LoadMacro(Form("%s/FitPlots.C",fitplotmacrodir.Data()));
  gStyle->SetOptStat(0);

        
  //set all the paths for the averaged plots

  Init();
  LoadFileNamespPbVsMCtemplates(pthad);
   histo = new TH1D *[nhistos];
  err = new TGraphAsymmErrors *[nhistos];
  ltscale = new TLatex *[nhistos];        

  cout <<"------> You are comparing pPb Vs MC pPb-templates" << endl;        
  // loop on histos
  for(Int_t k=0; k<nhistos; k++){ //GetHisto has to get graph/syst uncer===========================
  cout <<"---------> Looping on histos... " << k << endl;                
    if(k<3) { 
      if(skip3to5){
	if(k%3==0)continue;
      }

      histo[k] = GetHistoAndSyst(k,"fhDaverage","AverageSystematicUncertainty",err[k], ltscale[k]); 
	
      histo[k]->SetMarkerColor(1); 
      histo[k]->SetLineColor(1); 
      histo[k]->SetMarkerStyle(20); 
	
    } // get p-Pb
      
    if(k>=3 && k <6) { 
    if(skip3to5){
      if(k%3==0)continue;
    }

      h = GetHisto(k,"cDeltaPhi","hCorrDeltaPhi"); 
      if(reflTempl)histo[k]=AliHFCorrelationUtils::ReflectHisto(h,0.5);
      else histo[k]=h;
      histo[k]->SetMarkerColor(kMagenta+1);
      histo[k]->SetLineColor(kMagenta+1); 
      //histo[k]->SetMarkerColor(kOrange+2);
      //histo[k]->SetLineColor(kOrange+2); 
      histo[k]->SetLineWidth(1); 
      histo[k]->SetMarkerStyle(kOpenSquare); 
     
    } // get p-p
    if(k>=6 && k <9) { 
    if(skip3to5){
      if(k%3==0)continue;
    }

      h = GetHisto(k,"cDeltaPhi","hCorrDeltaPhi"); 
      if(reflTempl)histo[k]=AliHFCorrelationUtils::ReflectHisto(h,0.5);
      else histo[k]=h;
      histo[k]->SetMarkerColor(kGreen+2); 
      histo[k]->SetLineColor(kGreen+2);
      histo[k]->SetLineWidth(1); 
      histo[k]->SetMarkerStyle(kOpenCircle); 
    } // get p-p
    if(k>=9) { 
    if(skip3to5){
      if(k%3==0)continue;
    }

      h = GetHisto(k,"cDeltaPhi","hCorrDeltaPhi"); 
      if(reflTempl)histo[k]=AliHFCorrelationUtils::ReflectHisto(h,0.5);
      else histo[k]=h;    
      histo[k]->SetMarkerColor(kBlue);
      histo[k]->SetLineColor(kBlue); 
      histo[k]->SetLineWidth(1); 
      histo[k]->SetMarkerStyle(kOpenDiamond); 
      histo[k]->SetMarkerSize(1.5); 
    } // get p-p
  }

  Double_t BaselinepPb_Per0000, BaselineErrpPb_Per0000;// = 0;
  Double_t BaselinepPb_Per2010, BaselineErrpPb_Per2010;// = 0;
  Double_t BaselinepPb_Per2011, BaselineErrpPb_Per2011;// = 0;
  
  TH1D ** subtractedhisto = new TH1D *[nhistos];
  TGraphAsymmErrors ** suberr = new TGraphAsymmErrors *[nhistos];
  TGraphAsymmErrors **grbase=new TGraphAsymmErrors*[nhistos];
  TGraphAsymmErrors **grv2=new TGraphAsymmErrors*[nhistos]; 

 
  for(Int_t k=0; k<nhistos; k++){
    if(skip3to5){
      if(k%3==0)continue;
    }
    if(k<3) {
      Printf("Checking pointers, bin %d %p, %p",k,histo[k],err[k]);
      subtractedhisto[k] = GetPedestalHistoAndSystAndSubtractPedpPb("pPb",k,histo[k],err[k],suberr[k],"CanvasBaselineVariationTrendPedestal",grbase[k],grv2[k]);
      grbase[k]->SetFillStyle(3002);
      grbase[k]->SetFillColor(kBlue-7);
      grbase[k]->SetLineColor(kBlue-7);
      
      if(grv2[k]){
	grv2[k]->SetFillStyle(3002);
	grv2[k]->SetFillColor(kMagenta);  
      }
      suberr[k]->SetLineColor(kBlue);
      cout<<"sub -> "<<subtractedhisto[k]->GetBinContent(5)<<endl;
      //subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.7,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));    
    }
    if(k>=3 && k <6) {
      GetBaselne(BaselinepPb_Per0000,BaselineErrpPb_Per0000, histo[k]);
      cout <<"------Perugia8pPb-------> " <<BaselinepPb_Per0000 << " +/- " <<BaselineErrpPb_Per0000<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb_Per0000);
      
    }  // subtract Perugia8
    
    if(k>=6 && k <9) {
      GetBaselne(BaselinepPb_Per2010,BaselineErrpPb_Per2010, histo[k]);
      cout <<"------Perugia2010pPb-------> " <<BaselinepPb_Per2010 << " +/- " <<BaselineErrpPb_Per2010<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb_Per2010);
      
    }  // subtract Perugia2010
    
    if(k>=9) {
      GetBaselne(BaselinepPb_Per2011,BaselineErrpPb_Per2011, histo[k]);
      cout <<"------Perugia2011pPb-------> " <<BaselinepPb_Per2011 << " +/- " <<BaselineErrpPb_Per2011<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb_Per2011);
      
    }  // subtract Perugia2011
    
    
    
    subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.8,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));
  }
  // if(pthad>=1){
 if(pthad.Contains("99")){
 
    subtractedhisto[1]->SetMaximum(2.);
    //    subtractedhisto[2]->SetMaximum(3.2);
    subtractedhisto[2]->SetMaximum(4.2);
  }
  else if(pthad.Contains("0.5")) subtractedhisto[2]->SetMaximum(4.2);
  else if(pthad.Contains("0.3")) subtractedhisto[2]->SetMaximum(4.6);

  Int_t ihelpskip=0;
  if(skip3to5){
    ihelpskip=1;
  }

  TLegend * legend = new TLegend(0.14,0.748,0.51,0.748);
  //  TLegend * legend = new TLegend(0.14,0.765,0.51,0.765);
  //  TLegend * legend = new TLegend(0.15,0.75,0.5,0.88);
  legend->SetFillColor(0);
  legend->SetTextSize(0.028);
  legend->SetBorderSize(0);
  legend->AddEntry(histo[0+ihelpskip],"p-Pb #sqrt{#it{s}_{NN}}=5.02 TeV Data ","lep");
 
  TLegend * legendb = new TLegend(0.14,0.6,0.51,0.72);
  // TLegend * legendb = new TLegend(0.14,0.63,0.51,0.74);
  legendb->SetFillColor(0);
  legendb->SetTextSize(0.028);
  legendb->SetBorderSize(0);
  legendb->SetHeader("Simulations, pp #sqrt{#it{s}}=5.02 TeV");
  // legendb->AddEntry(histo[3],"Pythia8","lep");
  legendb->AddEntry(histo[3+ihelpskip],"Pythia6, Perugia 0","lep");
  legendb->AddEntry(histo[6+ihelpskip],"Pythia6, Perugia2010","lep");
  legendb->AddEntry(histo[9+ihelpskip],"Pythia6, Perugia2011","lep");

  //   legend->AddEntry(histo[0],"p-Pb #sqrt{s_{NN}}=5.02 TeV Data ","lep");
  //   legend->AddEntry(histo[3],"p-Pb #sqrt{s_{NN}}=5.02 TeV MC-Pythia8","lep");
  //   legend->AddEntry(histo[6],"p-Pb #sqrt{s_{NN}}=5.02 TeV MC-Perugia2010","lep");
  //   legend->AddEntry(histo[9],"p-Pb #sqrt{s_{NN}}=5.02 TeV MC-Perugia2011","lep");
  //  legend->AddEntry(histo[0],"p-Pb #sqrt{s}=5.02 TeV Data ","lep");
  //   legend->AddEntry(histo[3],"Pythia Pythia8 pp #sqrt{s}=5.02 TeV","lep");
  //   legend->AddEntry(histo[6],"Pythia Perugia2010 pp #sqrt{s}=5.02 TeV","lep");
  //   legend->AddEntry(histo[9],"Pythia Perugia2011 pp #sqrt{s}=5.02 TeV","lep");

  TLegend * legend2 = new TLegend(0.62,0.64,0.82,0.72);
  //  TLegend * legend2 = new TLegend(0.62,0.67,0.82,0.74);
  //  TLegend * legend2 = new TLegend(0.15,0.67,0.37,0.747);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);
  legend2->SetTextSize(0.02555);
  //legend2->AddEntry(box1[1],"correlated syst. pp","f");
  legend2->AddEntry(grbase[1],"baseline uncertainty","f");
  // legend2->AddEntry(box2[4],"v2 subtr p-Pb","f");

  // TPaveText *alice = new TPaveText(0.62,0.84,0.82,0.89,"NDC");
  //   alice->SetBorderSize(0);
  //   alice->SetFillColor(0);
  //   alice->SetTextFont(42);
  //   alice->SetTextSize(0.04);
  //   alice->AddText("ALICE Preliminary");

  TPaveText *alice = new TPaveText(0.62,0.725,0.82,0.76,"NDC");
  // TPaveText *alice = new TPaveText(0.62,0.74,0.82,0.79,"NDC");
  alice->SetBorderSize(0);
  alice->SetFillColor(0);
  alice->SetTextFont(42);
  alice->SetTextSize(0.04);
  //  alice->AddText("ALICE Preliminary");//commented

  TPaveText *fitvalueslow = new TPaveText(0.136,0.855,0.7,0.893,"NDC");
  //TPaveText *fitvalueslow = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  //TPaveText *fitvalueslow = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  fitvalueslow->SetBorderSize(0);
  fitvalueslow->SetFillColor(0);
  fitvalueslow->SetTextFont(42);
  fitvalueslow->AddText("D meson  - charged particle correlation");
  //  fitvalueslow->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvalueslow1 = new TPaveText(0.134,0.82,0.42,0.855,"NDC");
  fitvalueslow1->SetBorderSize(0);
  fitvalueslow1->SetFillColor(0);
  fitvalueslow1->SetTextFont(42);
  fitvalueslow1->AddText("Average D^{0},D^{+},D^{*+}");
  TPaveText *fitvalueslow2 = new TPaveText(0.15,0.774,0.84,0.81,"NDC");
  //  TPaveText *fitvalueslow2 = new TPaveText(0.18,0.79,0.84,0.843,"NDC");
  fitvalueslow2->SetBorderSize(0);
  fitvalueslow2->SetFillColor(0);
  fitvalueslow2->SetTextFont(42);
  //  fitvalueslow2->AddText(Form("3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} > %.1f GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad));
 if(pthad.Contains("0.3"))  fitvalueslow2->AddText("3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  else if(pthad.Contains("1.0"))  fitvalueslow2->AddText("3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} >1 GeV/#it{c}, |#Delta#eta| < 1.0 ");


  TPaveText *fitvaluesmid = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  //  TPaveText *fitvaluesmid = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  fitvaluesmid->SetBorderSize(0);
  fitvaluesmid->SetFillColor(0);
  fitvaluesmid->SetTextFont(42);
  fitvaluesmid->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvaluesmid2 = new TPaveText(0.15,0.774,0.84,0.81,"NDC"); 
  //TPaveText *fitvaluesmid2 = new TPaveText(0.18,0.79,0.84,0.843,"NDC");
  fitvaluesmid2->SetBorderSize(0);
  fitvaluesmid2->SetFillColor(0);
  fitvaluesmid2->SetTextFont(42);
  //  fitvaluesmid2->AddText(Form("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} > %.1f GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad));
  if(pthad.Contains("0.3"))  fitvaluesmid2->AddText("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  else if(pthad.Contains("1.0"))  fitvaluesmid2->AddText("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} >1 GeV/#it{c}, |#Delta#eta| < 1.0 ");


  TPaveText *fitvalueshigh = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  //TPaveText *fitvalueshigh = new TPaveText(0.136,0.85,0.88,0.893,"NDC");
  fitvalueshigh->SetBorderSize(0);
  fitvalueshigh->SetFillColor(0);
  fitvalueshigh->SetTextFont(42);
  fitvalueshigh->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  TPaveText *fitvalueshigh2 = new TPaveText(0.15,0.774,0.84,0.81,"NDC");
  //  TPaveText *fitvalueshigh2 = new TPaveText(0.18,0.79,0.84,0.843,"NDC");
  fitvalueshigh2->SetBorderSize(0);
  fitvalueshigh2->SetFillColor(0);
  fitvalueshigh2->SetTextFont(42);
  //  fitvalueshigh2->AddText(Form("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} > %.1f GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad));
  if(pthad.Contains("0.3"))  fitvalueshigh2->AddText("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, 0.3 < #it{p}_{T}^{assoc} <1 GeV/#it{c}, |#Delta#eta| < 1.0 ");
  else if(pthad.Contains("1.0"))  fitvalueshigh2->AddText("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} >1 GeV/#it{c}, |#Delta#eta| < 1.0 ");

  // TPaveText *fitvalueslow = new TPaveText(0.136,0.77,0.75,0.89,"NDC");
  //   fitvalueslow->SetBorderSize(0);
  //   fitvalueslow->SetFillColor(0);
  //   fitvalueslow->SetTextFont(42);
  //   fitvalueslow->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  //   fitvalueslow->AddText(Form("3 < #it{p}_{T}^{D} < 5 GeV/#it{c}, #it{p}_{T}^{assoc} > %.1f GeV/#it{c} ",pthad));
  //    fitvalueslow->AddText(Form("|#Delta#eta| < %.1f ",1.0));
    
  //   TPaveText *fitvaluesmid = new TPaveText(0.136,0.8,0.85,0.89,"NDC"); 
  //    fitvaluesmid->SetBorderSize(0);
  //   fitvaluesmid->SetFillColor(0);
  //   fitvaluesmid->SetTextFont(42);
  //   fitvaluesmid->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  //   fitvaluesmid->AddText(Form("5 < #it{p}_{T}^{D} < 8 GeV/#it{c}, #it{p}_{T}^{assoc} > %.1f GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad));
     
  //   TPaveText *fitvalueshigh = new TPaveText(0.136,0.8,0.85,0.89,"NDC");
  //    fitvalueshigh->SetBorderSize(0);
  //   fitvalueshigh->SetFillColor(0);
  //   fitvalueshigh->SetTextFont(42);
  //   fitvalueshigh->AddText("D meson (average D^{0},D^{+},D^{*+}) - charged particle correlation");
  //   fitvalueshigh->AddText(Form("8 < #it{p}_{T}^{D} < 16 GeV/#it{c}, #it{p}_{T}^{assoc} > %.1f GeV/#it{c}, |#Delta#eta| < 1.0 ",pthad));

  Double_t startline,endline,wi;
  wi=subtractedhisto[4]->GetBinWidth(1);
  startline=subtractedhisto[4]->GetBinLowEdge(1);
  endline=subtractedhisto[4]->GetBinLowEdge(subtractedhisto[4]->GetNbinsX())+wi+0.3;
  // cout<<"******************* "<<startline<<"  "<<endline<<endl;
  TLine* line=new TLine(startline, 0, endline, 0);
  line->SetLineStyle(2);

  //TCanvas * ptlow = new TCanvas("Correlation_ptlow","pT 3-5 correlations",0,0,800,800);
  TCanvas * ptmid = new TCanvas("Correlation_ptmid","pT 5-8 correlations",0,0,800,800);
  TCanvas * pthigh = new TCanvas("Correlation_pthigh","pT 8-16 correlations",0,0,800,800);
   

  ptmid->cd();
  ptmid->SetTicky();
  ptmid->SetTickx();
  ptmid->SetLeftMargin(0.125);
  TH1D* h=new TH1D(*subtractedhisto[2]);
  h->Reset();
  h->GetXaxis()->SetLimits(startline,endline);
  h->SetLineColor(0);
  h->GetXaxis()->Delete();
  h->GetYaxis()->Delete();
  h->Draw();
  subtractedhisto[1]->Draw("epsame");
  //box2[1]->Draw("same");
  grbase[1]->Draw("E2");
  line->Draw(); 
  subtractedhisto[4]->Draw("sameep");
  subtractedhisto[7]->Draw("sameep");
  subtractedhisto[10]->Draw("sameep");
  TH1D* hpy1=new TH1D(*subtractedhisto[4]);
  TH1D* hpy2=new TH1D(*subtractedhisto[7]);
  TH1D* hpy3=new TH1D(*subtractedhisto[10]);
  hpy1->GetXaxis()->Delete();
  hpy1->GetYaxis()->Delete();
  hpy2->GetXaxis()->Delete();
  hpy2->GetYaxis()->Delete();
  hpy3->GetXaxis()->Delete();
  hpy3->GetYaxis()->Delete();
  hpy1->SetMarkerSize(1.25);
  hpy2->SetMarkerSize(1.25);
  hpy3->SetMarkerSize(1.7);
  hpy1->Draw("samep");
  hpy2->Draw("samep");
  hpy3->Draw("samep");
  subtractedhisto[1]->Draw("sameep");
 
  suberr[1]->Draw("E2");
  ltscale[1]->Draw();//commented

  legend->Draw("same");
  legendb->Draw("same");
  legend2->Draw();
  fitvalueslow->Draw("same");
  //  fitvaluesmid->Draw("same");
  alice->Draw();
  fitvaluesmid2->Draw("same");
  fitvalueslow1->Draw("same");
 
  pthigh->cd();
  pthigh->SetTicky();
  pthigh->SetTickx();
  pthigh->SetLeftMargin(0.125);
  h->Draw();
  subtractedhisto[2]->Draw("epsame");
  //box2[2]->Draw("same");
  grbase[2]->Draw("E2");
  line->Draw(); 
  subtractedhisto[5]->Draw("sameep");
  subtractedhisto[8]->Draw("sameep");
  subtractedhisto[11]->Draw("sameep");
  TH1D* hpy1=new TH1D(*subtractedhisto[5]);
  TH1D* hpy2=new TH1D(*subtractedhisto[8]);
  TH1D* hpy3=new TH1D(*subtractedhisto[11]);
  hpy1->GetXaxis()->Delete();
  hpy1->GetYaxis()->Delete();
  hpy2->GetXaxis()->Delete();
  hpy2->GetYaxis()->Delete();
  hpy3->GetXaxis()->Delete();
  hpy3->GetYaxis()->Delete();
  hpy1->SetMarkerSize(1.25);
  hpy2->SetMarkerSize(1.25);
  hpy3->SetMarkerSize(1.7);
 
  hpy1->Draw("samep");
  hpy2->Draw("samep");
  hpy3->Draw("samep");
  subtractedhisto[2]->Draw("sameep");
  suberr[2]->Draw("E2");
  ltscale[2]->Draw();//commented
  legend->Draw("same");
  legendb->Draw("same");
  legend2->Draw();
  fitvalueslow->Draw("same");
  //  fitvalueshigh->Draw("same");
  alice->Draw();
  fitvalueshigh2->Draw("same");
  fitvalueslow1->Draw("same");

  ptmid->SaveAs(Form("plotComparison_%sAverage_pPb-MC_58_%s.root",avType.Data(),pthad.Data()));
  ptmid->SaveAs(Form("plotComparison_%sAverage_pPb-MC_58_%s.pdf",avType.Data(),pthad.Data()));
  ptmid->SaveAs(Form("plotComparison_%sAverage_pPb-MC_58_%s.eps",avType.Data(),pthad.Data()));
  ptmid->SaveAs(Form("plotComparison_%sAverage_pPb-MC_58_%s.gif",avType.Data(),pthad.Data()));

  pthigh->SaveAs(Form("plotComparison_%sAverage_pPb-MC_816_%s.root",avType.Data(),pthad.Data()));
  pthigh->SaveAs(Form("plotComparison_%sAverage_pPb-MC_816_%s.pdf",avType.Data(),pthad.Data()));
  pthigh->SaveAs(Form("plotComparison_%sAverage_pPb-MC_816_%s.eps",avType.Data(),pthad.Data()));
  pthigh->SaveAs(Form("plotComparison_%sAverage_pPb-MC_816_%s.gif",avType.Data(),pthad.Data()));
 
}
//_______________________________________________________________________
void DoComparison_pPbVsMC(Double_t pthad = 0.3){
 

    
  gStyle->SetOptStat(0);

    
  // set all the paths for the averaged plots
  for(Int_t k = 0; k<nhistos; k++){
    filenames[k] = directory;
  }
  for(Int_t k = 0; k<2; k++){
    pedestalfilenames[k] = directory;
  }
        

  LoadFileNamespPbVsMCtemplates(pthad);
  cout <<"------> You are comparing pPb Vs MC pPb-templates" << endl;        
  TH1D ** histo = new TH1D *[nhistos];

    
  // loop on histos
  for(Int_t k=0; k<nhistos; k++){
    if(skip3to5){
      if(k%3==0)continue;
    }
    if(k<3) { histo[k] = GetHisto(k,"cDraw","fhDaverage"); histo[k]->SetMarkerColor(2); histo[k]->SetMarkerStyle(20); } // get p-Pb
    if(k>=3 && k <6) { histo[k] = GetHisto(k,"cDeltaPhi","hCorrDeltaPhi"); histo[k]->SetMarkerColor(3); histo[k]->SetMarkerStyle(21); } // get p-p
    if(k>=6 && k <9) { histo[k] = GetHisto(k,"cDeltaPhi","hCorrDeltaPhi"); histo[k]->SetMarkerColor(4); histo[k]->SetMarkerStyle(22); } // get p-p
    if(k>=9) { histo[k] = GetHisto(k,"cDeltaPhi","hCorrDeltaPhi"); histo[k]->SetMarkerColor(6); histo[k]->SetMarkerStyle(29); } // get p-p

    //if(k>=3) { histo[k] = GetHisto(k,"cDeltaPhi","hCorrDeltaPhi"); histo[k]->SetMarkerColor(4); histo[k]->SetMarkerStyle(21); } // get MC
    histo[k]->Rebin(rebin); histo[k]->Scale(1./rebin);
  }

    
    

  Double_t BaselinepPb, BaselineErrpPb;// = 0;
  Double_t BaselinepPb_Per0000, BaselineErrpPb_Per0000;// = 0;
  Double_t BaselinepPb_Per2010, BaselineErrpPb_Per2010;// = 0;
  Double_t BaselinepPb_Per2011, BaselineErrpPb_Per2011;// = 0;
    
  TH1D ** subtractedhisto = new TH1D *[nhistos];
  for(Int_t k=0; k<nhistos; k++){
    if(skip3to5){
      if(k%3==0)continue;
    }
    if(k<3) {
      GetBaselne(BaselinepPb,BaselineErrpPb, histo[k]);
      cout <<"------pPb-------> " <<BaselinepPb << " +/- " <<BaselineErrpPb<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb);
            
    } // get p-Pb
    if(k>=3 && k <6) {
      GetBaselne(BaselinepPb_Per0000,BaselineErrpPb_Per0000, histo[k]);
      cout <<"------Perugia8pPb-------> " <<BaselinepPb_Per0000 << " +/- " <<BaselineErrpPb_Per0000<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb_Per0000);
            
    }  // subtract Perugia8
        
    if(k>=6 && k <9) {
      GetBaselne(BaselinepPb_Per2010,BaselineErrpPb_Per2010, histo[k]);
      cout <<"------Perugia2010pPb-------> " <<BaselinepPb_Per2010 << " +/- " <<BaselineErrpPb_Per2010<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb_Per2010);
            
    }  // subtract Perugia2010
        
    if(k>=9) {
      GetBaselne(BaselinepPb_Per2011,BaselineErrpPb_Per2011, histo[k]);
      cout <<"------Perugia2011pPb-------> " <<BaselinepPb_Per2011 << " +/- " <<BaselineErrpPb_Per2011<<endl;
      subtractedhisto[k] = subtractpedestal(k,histo[k],BaselinepPb_Per2011);
            
    }  // subtract Perugia2011
        
        
        
    subtractedhisto[k]->GetYaxis()->SetRangeUser(-0.5,2*subtractedhisto[k]->GetBinContent(subtractedhisto[k]->GetMaximumBin()));
  }
    

  // reflect histos
  TH1D ** reflectedhisto = new TH1D *[nhistos];
  cout << "test 3 "<< endl;
  for(Int_t k=0; k<nhistos; k++){
    if(skip3to5){
      if(k%3==0)continue;
    }
    reflectedhisto[k] = ReflectHistogram(subtractedhisto[k]);
    if(k<3) {reflectedhisto[k]->SetMarkerColor(2); reflectedhisto[k]->SetMarkerStyle(20);}
    if(k>=3 && k <6) {reflectedhisto[k]->SetMarkerColor(3); reflectedhisto[k]->SetMarkerStyle(21);}
    if(k>=6 && k <9) {reflectedhisto[k]->SetMarkerColor(4); reflectedhisto[k]->SetMarkerStyle(22);}
    if(k>=9) {reflectedhisto[k]->SetMarkerColor(6); reflectedhisto[k]->SetMarkerStyle(29);}

    reflectedhisto[k]->GetYaxis()->SetRangeUser(-0.5,2*reflectedhisto[k]->GetBinContent(reflectedhisto[k]->GetMaximumBin()));
  }
  cout << "test 4 "<< endl;
    
    
    
  TLegend * legend = new TLegend(0.1,0.8,0.60,0.9);
  legend->SetFillColor(0);
  legend->SetTextSize(0.025);
  Int_t ihelpskip=0;
  if(skip3to5){
    ihelpskip=1;
  }
  legend->AddEntry(histo[0+ihelpskip],"p-Pb @ 5.02 TeV Data ","lep");
  legend->AddEntry(histo[3+ihelpskip],"p-Pb @ 5.02 TeV MC-Pythia8","lep");
  legend->AddEntry(histo[6+ihelpskip],"p-Pb @ 5.02 TeV MC-Perugia2010","lep");
  legend->AddEntry(histo[9+ihelpskip],"p-Pb @ 5.02 TeV MC-Perugia2011","lep");
    
    
  TPaveText *fitvalueslow = new TPaveText(0.515,0.66,0.85,0.78,"NDC");
  fitvalueslow->SetBorderSize(0);
  fitvalueslow->SetFillColor(0);
  fitvalueslow->AddText("3 GeV/c < p_{T}^{D} < 5 GeV/c");
  fitvalueslow->AddText(Form("p_{T}^{Ass} > %.1f GeV/c ",pthad));
  fitvalueslow->AddText(Form("|#Delta#eta| < %.1f ",1.0));
    
  TPaveText *fitvaluesmid = new TPaveText(0.515,0.66,0.85,0.78,"NDC");
  fitvaluesmid->SetBorderSize(0);
  fitvaluesmid->SetFillColor(0);
  fitvaluesmid->AddText("5 GeV/c < p_{T}^{D} < 8 GeV/c");
  fitvaluesmid->AddText(Form("p_{T}^{Ass} > %.1f GeV/c ",pthad));
  fitvaluesmid->AddText(Form("|#Delta#eta| < %.1f ",1.0));
    
  TPaveText *fitvalueshigh = new TPaveText(0.515,0.66,0.85,0.78,"NDC");
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
  
  if(!skip3to5){
    ptlow->cd();
    subtractedhisto[0]->Draw("ep");
    subtractedhisto[3]->Draw("sameep");
    subtractedhisto[6]->Draw("sameep");
    subtractedhisto[9]->Draw("sameep");
    legend->Draw("same");
    fitvalueslow->Draw("same");
    
    ptlowRef->cd();
    reflectedhisto[0]->Draw("ep");
    reflectedhisto[3]->Draw("sameep");
    reflectedhisto[6]->Draw("sameep");
    reflectedhisto[9]->Draw("sameep");
    legend->Draw("same");
    fitvalueslow->Draw("same");
  }
    
  ptmid->cd();
  subtractedhisto[1]->Draw("ep");
  subtractedhisto[4]->Draw("sameep");
  subtractedhisto[7]->Draw("sameep");
  subtractedhisto[10]->Draw("sameep");
  legend->Draw("same");
  fitvaluesmid->Draw("same");
    
  ptmidRef->cd();
  reflectedhisto[1]->Draw("ep");
  reflectedhisto[4]->Draw("sameep");
  reflectedhisto[7]->Draw("sameep");
  reflectedhisto[10]->Draw("sameep");
  legend->Draw("same");
  fitvaluesmid->Draw("same");
    
    
  pthigh->cd();
  subtractedhisto[2]->Draw("ep");
  subtractedhisto[5]->Draw("sameep");
  subtractedhisto[8]->Draw("sameep");
  subtractedhisto[11]->Draw("sameep");
  legend->Draw("same");
  fitvalueshigh->Draw("same");
    
  pthighRef->cd();
  reflectedhisto[2]->Draw("ep");
  reflectedhisto[5]->Draw("sameep");
  reflectedhisto[8]->Draw("sameep");
  reflectedhisto[11]->Draw("sameep");
  legend->Draw("same");
  fitvalueshigh->Draw("same");
    
    
    
    
  // saving the canvases in .root and .png
  TString ptoutput="",ptoutputRef="", direcname="", direcnameRef="";
  ptoutput += Form("%sAveragepPbVsMCtemplates",avType.Data());
  direcname += "Output_Plots/pPbVsMCtemplatesPlots";
    
    
    
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
    
    
  ptoutputRef += Form("%sAveragepPbVsMCtemplates",avType.Data());
  direcnameRef += "Output_PlotsRef/pPbVsMCtemplatesPlots";
    
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

void LoadFileNamespPbVsMCtemplates(TString pthadron){
  filenames[0] = Form("%s/%sAveragepPbDzeroDstarDplus3to5_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  filenames[1] = Form("%s/%sAveragepPbDzeroDstarDplus5to8_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  filenames[2] = Form("%s/%sAveragepPbDzeroDstarDplus8to16_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
 // filenames[3] = Form("%s/%sAveragepPbDzeroDstarDplus16to24_assoc%s.root",inputdatadirectory.Data(),avType.Data(),pthadron.Data());
  
    TString strsystHelp=strSystemFDtempl;
    if(strSystemFDtempl.EqualTo("none")){
      strsystHelp="pPb";
    }

  //load MC Templates1
  filenames[3] = Form("%s/%sCorrelationPlotsPerugia0wBoostPtDzerofromC3To5_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());//in the preliminary it was Puthia8
  filenames[4] = Form("%s/%sCorrelationPlotsPerugia0wBoostPtDzerofromC5To8_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());
  filenames[5] = Form("%s/%sCorrelationPlotsPerugia0wBoostPtDzerofromC8To16_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());
  
  
  filenames[6] = Form("%s/%sCorrelationPlotsPerugia2010wBoostPtDzerofromC3To5_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());//in the preliminary it was Puthia8
  filenames[7] = Form("%s/%sCorrelationPlotsPerugia2010wBoostPtDzerofromC5To8_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());
  filenames[8] = Form("%s/%sCorrelationPlotsPerugia2010wBoostPtDzerofromC8To16_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());
  
  
  filenames[9] = Form("%s/%sCorrelationPlotsPerugia2011wBoostPtDzerofromC3To5_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());//in the preliminary it was Puthia8
  filenames[10] = Form("%s/%sCorrelationPlotsPerugia2011wBoostPtDzerofromC5To8_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());
  filenames[11] = Form("%s/%sCorrelationPlotsPerugia2011wBoostPtDzerofromC8To16_ptAssall%s_DeltaEta10.root",inputtemplatedirecotry.Data(),strsystHelp.Data(),pthadron.Data());
  
  //load pedestals from fit (STILL TEMPORARY)
  for(Int_t k = 0; k<2; k++){
    pedestalfilenames[k] = baselinedirectory;
  }
  pedestalfilenames[0] += Form("/Trends_pp/CanvasBaselineVariationTrendPedestal_pthad%s.root",pthadron.Data());//pp
  pedestalfilenames[1] += Form("/Trends_pPb/CanvasBaselineVariationTrendPedestal_pthad%s.root",pthadron.Data());//pPb
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
    

  DoComparison_pPbVsMCTEST("0.3to1.0");
  DoComparison_pPbVsMCTEST("0.3to99.0");
  DoComparison_pPbVsMCTEST("1.0to99.0");
  cout << "Done ! Check yours plots in folder ...  " << endl;
    
    
    
}

//_______________________________________________________________________
TH1D * GetPedestalHistoAndSystAndSubtractPedpPb(TString system,Int_t i,TH1D *histo, TGraphAsymmErrors* gr,TGraphAsymmErrors *&grout, TString canvasname,TGraphAsymmErrors *&grbaseOut,TGraphAsymmErrors *&grv2Out){
  //i=D pt bin
  //Double_t* arrbox1, Double_t* arrbox2
  Double_t value = 0, pedestal=0;
  TGraphAsymmErrors *grBase,*grV2;
  grbaseOut=new TGraphAsymmErrors();
  grbaseOut->SetName(Form("grbaselineUncFull_%sBin%d",system.Data(),i));
  grv2Out=new TGraphAsymmErrors();
  grv2Out->SetName(Form("grbaselineUncV2_%sBin%d",system.Data(),i));
  TH1D* h;  
  Double_t xuncFull,errxuncFull;
  Double_t xuncv2,errxuncv2;
  Int_t bin,bingr;

  if(system.Contains("pPb")){
    // get pedestal from fit outputs
    TString path = pedestalfilenames[1];//pPb ----------------
    cout << "pPb -->  Reading File from path: " << path << endl;
    
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    grBase=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fBaselineVariationSystematicsPedestal");    
    grV2=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fv2SystematicsPedestal");    
  }
  else if(system.Contains("pp")){
    // get pedestal from fit outputs
    TString path = pedestalfilenames[0];//pp ----------------
    grV2=0x0;
    cout << "pp -->  Reading File from path: " << path << endl;
    TFile * file = TFile::Open(path.Data(),"READ");
    TCanvas* c=(TCanvas*)file->Get(canvasname.Data());
    TH1D* h = (TH1D*)c->GetListOfPrimitives()->FindObject("FinalTrendPedestal");
    grBase=(TGraphAsymmErrors*)c->GetListOfPrimitives()->FindObject("fBaselineVariationSystematicsPedestal");    
  }
    
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
  grBase->GetPoint(bingr,x,y);
  Printf("histo: x=%f, graph: %f",h->GetBinCenter(bin),x);
  erryl=grBase->GetErrorYlow(bingr);
  erryh=grBase->GetErrorYhigh(bingr);
  
  grbaseOut->SetPoint(0,xuncFull,0);
  grbaseOut->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
  Printf("Point set to grbaseOut");
  if(grV2){
    Printf("Setting ve gr prop");
    grV2->GetPoint(bingr,x,y);
    erryl=grV2->GetErrorYlow(bingr);
    erryh=grV2->GetErrorYhigh(bingr);	
    grv2Out->SetPoint(0,xuncFull,0);
    grv2Out->SetPointError(0,errxuncFull,errxuncFull,erryl,erryh);
  }
  grout=(TGraphAsymmErrors*)gr->Clone(Form("grSub_%sBin%d",system.Data(),i));
  
  Printf("gr out clone");
  TString nameoutput = histo->GetName();
  nameoutput += "_subtr_";
  nameoutput += "pedestal";
  nameoutput += Form("%d",i);
  
  TH1D *outputhisto = (TH1D*)histo->Clone(nameoutput.Data());
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
  
  outputhisto->SetXTitle("#Delta#varphi (rad)");
  outputhisto->SetYTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} - baseline (rad^{-1})");
  outputhisto->GetYaxis()->SetTitleOffset(1.5);
  return outputhisto;
    
}

