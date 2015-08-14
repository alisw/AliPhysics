TString strPtAss[3]={"0.3to1.0","1.0to99.0","0.3to99.0"};
TString strPtAssCanvas[3]={"0.3<#it{p}_{T}^{assoc}<1 GeV/#it{c}","#it{p}_{T}^{assoc}>1 GeV/#it{c}","#it{p}_{T}^{assoc}>0.3 GeV/#it{c}"};
TString strSystem[2]={"pp","pPb"};
TString strFitResultPP="";//"/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";
TString strFitResultPPb="";//"/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";

TString strFitResultMC[2]={"",""};
//TString strFitResultMC[2]={"/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/MCTemplates/Templates_pp_12May15/FitResults/",
//			   "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May19UseScriptPWGHF/MCTemplates/Templates_pPb_12May15/FitResults/"};

TString strquantityFile[5]={"NSYield","NSSigma","Pedestal","ASYield","ASSigma"};
Double_t maxRangePP[3][5]={{4,1,5,4,1},
			   {4,1,3,4,1},
			   {6,1,5,4,1}};

Double_t maxRangePPb[3][5]={{4,1,13,4,1},
			   {4,1,5,4,1},
			   {6,1,13,4,1}};


TString yaxisTitle[5]={"Associated yield","#sigma_{fit} (rad)","Baseline (rad^{-1})","Associated yield","#sigma_{fit} (rad)"};
Double_t leftMarginCanvas=0.17;
Double_t rightMarginCanvas=0.055;
Double_t bottomMarginCanvas=0.13;
Double_t topMarginCanvas=0.07;
const Int_t nmodels=6;
Bool_t includemodel[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE};
TString strModelDir[nmodels]={"Perugia0","Perugia2010","Perugia2011","POWHEG","PYTHIA8","HERWIG"};
Color_t modelColors[nmodels]={kMagenta+1,kGreen+2,kBlue,kRed+2,kViolet,kCyan};
Int_t modelMarkerStyle[nmodels]={kOpenSquare,kOpenCircle,kOpenDiamond,kOpenDiamond,26,28};
TH1D **hMC;
TGraphAsymmErrors **grMC;

void IncludePerugia0(Bool_t incl){
  includemodel[0]=incl;
}
void IncludePerugia2010(Bool_t incl){
  includemodel[1]=incl;
}
void IncludePerugia2011(Bool_t incl){
  includemodel[2]=incl;
}
void IncludePythia8(Bool_t incl){
  includemodel[4]=incl;
}
void IncludePowheg(Bool_t incl){
  includemodel[3]=incl;
}
void IncludeHerwig(Bool_t incl){
  includemodel[5]=incl;
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

TCanvas* CompareNSyieldPPtoPPb(Int_t binass){
  ComparePPtoPPb(binass,0);
}

TCanvas* CompareNSsigmaPPtoPPb(Int_t binass){
  ComparePPtoPPb(binass,1);
}

TCanvas* ComparePedestalPPtoPPb(Int_t binass){
  ComparePPtoPPb(binass,2);
}


TCanvas* ComparePPtoPPb(Int_t binass,Int_t quantity){
  Printf("Opening file: %s", Form("%s/Trends_pp/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()));
  TFile *f=TFile::Open(Form("%s/Trends_pp/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
  TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *grPP=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  grPP->SetName(Form("%sPP",grPP->GetName()));
  TH1D *hPP=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  hPP->SetName(Form("%sPP",hPP->GetName()));
  
  f=TFile::Open(Form("%s/Trends_pPb/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
  c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *grPPb=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  grPPb->SetName(Form("%sPPb",grPPb->GetName()));
  TH1D *hPPb=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  hPPb->SetName(Form("%sPPb",hPPb->GetName()));
  
  TCanvas *cout=CreateCanvasWithDefaultStyle(Form("%sComparisonPPtoPPbBinAss%d",strquantityFile[quantity].Data(),binass));
  //new TCanvas(Form("NSyieldComparisonBinAss%d",binass),Form("NSyieldComparisonBinAss%d",binass),800,800);


  cout->cd();
  hPP->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
  hPP->SetYTitle(yaxisTitle[quantity].Data());
  hPP->GetYaxis()->SetTitleSize(0.04);
  hPP->GetYaxis()->SetTitleOffset(1.2);
  hPP->GetYaxis()->SetLabelSize(0.04);
  hPP->GetXaxis()->SetTitleSize(0.04);
  hPP->GetXaxis()->SetLabelSize(0.04);
  hPP->Draw();
  //  hPP->Draw("E0X0");// to avoid plotting the error along x
  hPP->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[binass][quantity],maxRangePPb[binass][quantity]));
 
  hPP->SetLineColor(kBlack);
  hPP->SetLineWidth(2);
  hPP->SetMarkerColor(kBlack);
  hPP->SetMarkerStyle(20);

  grPP->SetMarkerColor(kBlack);
  grPP->SetLineColor(kBlack);
  grPP->SetLineWidth(2);
  grPP->SetMarkerStyle(20);
  grPP->Draw("E2");
  
  

  hPPb->Draw("same");
  hPPb->SetLineColor(kRed);
  hPPb->SetLineWidth(2);
  hPPb->SetMarkerColor(kRed);
  hPPb->SetMarkerStyle(21);
  grPPb->SetMarkerColor(kRed);
  grPPb->SetLineColor(kRed);
  grPPb->SetLineWidth(2);
  grPPb->SetMarkerStyle(21);
  grPPb->Draw("E2");


  TLegend * legend = new TLegend(0.21,0.63,0.47,0.75);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.035);
  legend->AddEntry(hPP,"pp #sqrt{#it{s}}=7 TeV, |#it{y}^{D}_{cms}|<0.5","lep");
  legend->AddEntry(hPPb,"p-Pb #sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04","lep");
  
  legend->Draw();

  TLatex *tlALICE=new TLatex(0.75,0.85,"ALICE");
  tlALICE->SetNDC();
  tlALICE->SetTextSize(0.03);
  tlALICE->Draw();
  if(quantity!=2){
    
    TLatex *tlSide;
    if(quantity==0||quantity==1)tlSide=new TLatex(0.25,0.85,"Near side");
    else tlSide=new TLatex(0.75,0.85,"Away side");
    tlSide->SetNDC();
    tlSide->SetTextSize(0.03);
    tlSide->Draw();

  }
  TLatex *tlAssYieldPt=new TLatex(0.25,0.78,Form("%s, |#Delta#eta|<1",strPtAssCanvas[binass].Data()));
  tlAssYieldPt->SetNDC();
  tlAssYieldPt->SetTextSize(0.03);
  tlAssYieldPt->Draw();

  return cout;
  
}


void InitMCobjects(){
  hMC=new TH1D*[nmodels];
  grMC=new TGraphAsymmErrors*[nmodels];
}


TCanvas* CompareDatatoModels(Int_t collsystem,Int_t binass,Int_t quantity){
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

  TCanvas *cout;
  if(system==-1){
    cout=CreateCanvasWithDefaultStyle(Form("%sComparisonMCtoDataBothSystemsBinAss%d",strquantityFile[quantity].Data(),binass));
  }
  else{
    cout=CreateCanvasWithDefaultStyle(Form("%sComparisonMCto%sDataBinAss%d",strquantityFile[quantity].Data(),strSystem[system].Data(),binass));
  }
  Int_t neffmod=CountNmodels();
  TLegend * legend = new TLegend(0.61,0.54,0.95,0.54+neffmod*0.05);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.03);
  //new TCanvas(Form("NSyieldComparisonBinAss%d",binass),Form("NSyieldComparisonBinAss%d",binass),800,800);
  TLatex *tlCollSystem;
  TLatex *tlDrap;
  
  if(collsystem==-1){
    
    system=0;   
    TFile *f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
     
    TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
    grData[0]=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
    grData[0]->SetName(Form("%sPP",grData[0]->GetName()));
    hData[0]=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
    hData[0]->SetName(Form("%sPP",hData[0]->GetName()));
    
    system=1;
    f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
    c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
    grData[1]=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
    grData[1]->SetName(Form("%sPPb",grData[1]->GetName()));
    hData[1]=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
    hData[1]->SetName(Form("%sPPb",hData[1]->GetName()));
    system=0;   
  }
  else{
    TFile *f;
    if(system==0)f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPP.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
    else if(system==1)f=TFile::Open(Form("%s/Trends_%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb.Data(),strSystem[system].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
    
    TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
    grData[0]=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
    grData[0]->SetName(Form("%sPP",grData[0]->GetName()));
    hData[0]=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
    hData[0]->SetName(Form("%sPP",hData[0]->GetName()));
  }
    
    cout->cd();
    hData[0]->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
    hData[0]->SetYTitle(yaxisTitle[quantity].Data());
    hData[0]->GetYaxis()->SetTitleSize(0.04);
    hData[0]->GetYaxis()->SetTitleOffset(1.2);
    hData[0]->GetYaxis()->SetLabelSize(0.04);
    hData[0]->GetXaxis()->SetTitleSize(0.04);
    hData[0]->GetXaxis()->SetLabelSize(0.04);
    hData[0]->Draw();
    //  hData[0]->Draw("E0X0");// to avoid plotting the error along x
    hData[0]->GetYaxis()->SetRangeUser(0,TMath::Max(maxRangePP[binass][quantity],maxRangePPb[binass][quantity]));
    
    hData[0]->SetLineColor(kBlack);
    hData[0]->SetLineWidth(2);
    hData[0]->SetMarkerColor(kBlack);
    hData[0]->SetMarkerStyle(20);
    
    
    grData[0]->SetMarkerColor(kBlack);
    grData[0]->SetLineColor(kBlack);
    grData[0]->SetLineWidth(2);
    grData[0]->SetMarkerStyle(20);
    grData[0]->Draw("E2");
    

    
    if(collsystem==-1){    
      hData[1]->Draw("same");
      hData[1]->SetLineColor(kRed);
      hData[1]->SetLineWidth(2);
      hData[1]->SetMarkerColor(kRed);
      hData[1]->SetMarkerStyle(21);
      grData[1]->SetMarkerColor(kRed);
      grData[1]->SetLineColor(kRed);
      grData[1]->SetLineWidth(2);
      grData[1]->SetMarkerStyle(21);
      grData[1]->Draw("E2");
    }



    if(collsystem==-1){
      legend->AddEntry(hData[0],"pp #sqrt{#it{s}}=7 TeV, |#it{y}^{D}_{cms}|<0.5","lep");
      legend->AddEntry(hData[1],"p-Pb #sqrt{#it{s}_{NN}}=5.02 TeV, -0.96<#it{y}^{D}_{cms}<0.04","lep");      
    }
    else {    
      //      legend->AddEntry(hData[0],"data","lep");
      if(system==0){
	tlCollSystem=new TLatex(0.24,0.8,Form("pp #sqrt{#it{s}}=7 TeV"));
	tlDrap=new TLatex(0.24,0.75,Form("|#it{y}^{D}_{cms}|<0.5"));
      }
      else if(system==1) {
	tlCollSystem=new TLatex(0.24,0.8,Form("p-Pb #sqrt{#it{s}_{NN}}=5.02 TeV"));
	tlDrap=new TLatex(0.24,0.75,Form("-0.96<#it{y}^{D}_{cms}<0.04"));
      }
      tlDrap->SetNDC();
      tlDrap->SetTextSize(0.03);      
      
      tlCollSystem->SetNDC();
      tlCollSystem->SetTextSize(0.03);      
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
      
      cout->cd();
      hMC[kmc]->Draw("same");
      hMC[kmc]->SetLineColor(modelColors[kmc]);
      hMC[kmc]->SetLineWidth(2);
      hMC[kmc]->SetMarkerColor(modelColors[kmc]);
      hMC[kmc]->SetMarkerStyle(modelMarkerStyle[kmc]);
      grMC[kmc]->SetMarkerColor(modelColors[kmc]);
      grMC[kmc]->SetLineColor(modelColors[kmc]);
      grMC[kmc]->SetLineWidth(2);
      grMC[kmc]->SetMarkerStyle(21);
      grMC[kmc]->SetFillStyle(3001+kmc);
      grMC[kmc]->SetFillColor(modelColors[kmc]);
      grMC[kmc]->Draw("E5");
      
      legend->AddEntry(hMC[kmc],Form("%s",strModelDir[kmc].Data()),"lep");  
    }
    
    legend->Draw();
    
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
    TLatex *tlAssYieldPt=new TLatex(0.24,0.7,Form("%s, |#Delta#eta|<1",strPtAssCanvas[binass].Data()));
    tlAssYieldPt->SetNDC();
    tlAssYieldPt->SetTextSize(0.03);
    tlAssYieldPt->Draw();

    hData[0]->Draw("same");    
    grData[0]->Draw("E2");   
    
    if(collsystem!=-1){
      tlCollSystem->Draw();
      tlDrap->Draw();
    }
    else{
      hData[1]->Draw("same");    
      grData[1]->Draw("E2");   
    }
 
    return cout;
    
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
