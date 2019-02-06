///
/// \file AddChi2JJDecLowHighGraphs
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Add graphs with Chi2 dsitributions of different MC productions in same plot
///
/// Add graphs with Chi2 dsitributions of different MC productions in same plot. 
/// Root files with graphs produced with MakeDataMCComparisonPerSMClusterEbinAndChi2.C executed
/// with options titleMC "JJDecLow" and "JJDecHigh"
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TObject.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphErrors.h>
#endif

//------------------------------------------------------------------------------
/// Main method
///
/// \param quantity   :  Name of histogram with the parameter to study = "SMM02NoCut", "SMM02NoCut","SMM02",
///                      "SMM20LowM02NoCut","SMM20LowM02","SMM20HighM02NoCut","SMM20HighM02","SMNCell"
/// \param titleData  : Simplified string acronym of input data
/// \param opt        : Additional string in input file name containing analysis options
/// \param fileFormat : String with file format starting with .
//------------------------------------------------------------------------------
void AddChi2JJDecLowHighGraphs
(
 TString quantity = "SMM02NoCut",//"SMM02NoCut","SMM02","SMM20LowM02NoCut","SMM02LowM02","SMM20HighM02NoCut","SMM02HighM02","SMNCell"
 TString titleData = "LHC11cd_EMC7",
 TString opt       = "", //"TMDep_Neutral_V1_..." 
 TString fileFormat = ".eps" 
 )
{  
  // SM
  const Int_t nparam = 10;
  Double_t pbins[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Int_t ncol = 4;
  Int_t nrow = 3;
  
//  Float_t eBinsL[] = {2.5,3.5,5.5,7.5,9.5,11.5,13.5,16.5,18.5};
//  Float_t eBinsH[] = {16.5,18.5,20.5,22.5,24.5,26.5,50.5};

  Float_t eBinsL[] = { 7.5, 9.5,11.5,13.5,15.5,17.5,19.5};
  Float_t eBinsH[] = {17.5,19.5,24.5,29.5,34.5,39.5,49.5};
  
  const Int_t nEBins = 8;
  
  Int_t colorL[] = {1,4,2,6,8,9,kYellow-6,kCyan,kOrange+2,kViolet,kOrange-2};
  Int_t colorH[] = {9,kOrange+2,kViolet,kOrange-2,kRed-3,kRed+3,kBlue-3,kBlue+3};
  
  TFile* fileL = TFile::Open(Form("figures/%s/Chi2Distributions_%s_JJDecLow%s_PerSM_Ebin.root",
                                  quantity.Data(),titleData.Data(),opt.Data()));
  TFile* fileH = TFile::Open(Form("figures/%s/Chi2Distributions_%s_JJDecHigh%s_PerSM_Ebin.root",
                                  quantity.Data(),titleData.Data(),opt.Data()));
  
  printf("File JJ Low  %p %s\n",fileL,fileL->GetName());
  printf("File JJ High %p %s\n",fileH,fileH->GetName());
  
  TGraphErrors* gChi2L    [nEBins][nparam];
  TGraphErrors* gChi2LMin [nparam];
  TGraphErrors* gChi2H    [nEBins][nparam];
  TGraphErrors* gChi2HMin [nparam];
  
  for(Int_t ipbin = 0; ipbin < nparam; ipbin++)
  {
    gChi2LMin[ipbin] = (TGraphErrors*) fileL->Get(Form("Chi2Min_SM%d",ipbin));
    if( gChi2LMin[ipbin] )
    {
      gChi2LMin[ipbin]->SetMarkerStyle(20);
      gChi2LMin[ipbin]->SetMarkerColor(4);
      gChi2LMin[ipbin]->SetLineColor  (4);
    }
    
    gChi2HMin[ipbin] = (TGraphErrors*) fileH->Get(Form("Chi2Min_SM%d",ipbin));
    if( gChi2HMin[ipbin] )
    {
      gChi2HMin[ipbin]->SetMarkerStyle(24);
      gChi2HMin[ipbin]->SetMarkerColor(4);   
      gChi2HMin[ipbin]->SetLineColor  (4);   
    }
    //printf("\t SM %d, %p %p\n",ipbin,gChi2LMin[ipbin],gChi2LMin[ipbin]);
    
    for(Int_t iebin = 0; iebin < nEBins; iebin++)
    {
      gChi2L[iebin][ipbin] = (TGraphErrors*) fileL->Get(Form("Chi2_SM%d_EBin%d",ipbin,iebin));
      if(gChi2L[iebin][ipbin])
      {
        gChi2L[iebin][ipbin]->SetMarkerStyle(20);
        gChi2L[iebin][ipbin]->SetMarkerSize(4);
        gChi2L[iebin][ipbin]->SetMarkerColor(colorL[iebin]);
        gChi2L[iebin][ipbin]->SetLineColor(colorL[iebin]);
      }
      
      gChi2H[iebin][ipbin] = (TGraphErrors*)fileH->Get(Form("Chi2_SM%d_EBin%d",ipbin,iebin));
      if(gChi2H[iebin][ipbin])
      {
        gChi2H[iebin][ipbin]->SetMarkerStyle(24);
        gChi2H[iebin][ipbin]->SetMarkerSize(4);
        gChi2H[iebin][ipbin]->SetMarkerColor(colorH[iebin]);
        gChi2H[iebin][ipbin]->SetLineColor(colorH[iebin]);
      }
      //printf("\t ie %d, %p %p\n",iebin,gChi2L[iebin][ipbin],gChi2H[iebin][ipbin]);
    }
  }
  
  //
  TCanvas * gChi2graph = new TCanvas(Form("cChi2_graph_%s",quantity.Data()),
                                     Form("Chi2, %s",quantity.Data()),
                                     ncol*1000,nrow*1000);
  gChi2graph->Divide(ncol,nrow);
  
  TLegend *lgchi2 = new TLegend(0,0,1,1);
  lgchi2->SetFillColor(0);
  lgchi2->SetFillStyle(0);
  lgchi2->SetLineColor(0);
  lgchi2->SetBorderSize(0);
  lgchi2->SetTextSize(0.07);
  
  for(Int_t ipbin = 0; ipbin < nparam; ipbin++)
  {
    gChi2graph->cd(ipbin+1);
    
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    
    // LOW
    Bool_t set = kFALSE;
    for(Int_t iebin = 0; iebin < 5; iebin++)
    {
      if(!gChi2L[iebin][ipbin]) continue;
      
      gChi2L[iebin][ipbin]->GetHistogram()->SetTitleOffset(1.5,"Y");
      //gChi2L[iebin][ipbin]->GetHistogram()->SetYTitle("#chi^{2} / #nu");
      // gChi2L[iebin][ipbin]->SetAxisRange(0,0.8,"X");
      
      //printf("iebin %d, ipbin %d 0 %d %p\n",iebin,ipbin,0,h[0][iebin][ipbin]);
      
      if(!set) { 
        set = kTRUE;
        gChi2L[iebin][ipbin]->Draw("APL");
      }
      else
        gChi2L[iebin][ipbin]->Draw("PL");
      
      //      gChi2L[iebin][ipbin]->SetMaximum(200);
      //      gChi2L[iebin][ipbin]->SetMinimum(0.25);
      
      if(ipbin==0) 
      {
        lgchi2->AddEntry(gChi2L[iebin][ipbin],
                         gChi2L[iebin][ipbin]->GetTitle(),
                         "PL");
      }
      if(quantity.Contains("SM")) gChi2L[iebin][ipbin]->SetTitle(Form("SM %d",ipbin));
      
    } // iebin
    
    // HIGH
    for(Int_t iebin = 0; iebin < 2; iebin++)
    {
      if(!gChi2H[iebin][ipbin]) continue;
      
      gChi2H[iebin][ipbin]->Draw("PL");
      
      if(ipbin==0) 
      {
        lgchi2->AddEntry(gChi2H[iebin][ipbin],
                         gChi2H[iebin][ipbin]->GetTitle(),
                         "PL");
      }
    } // iebin
      //printf("\t end\n");
  } // param
  
  gChi2graph->cd(ncol*nrow-1);
  lgchi2->Draw();
  
  gChi2graph->cd(ncol*nrow);
  TLegend *legend2 = new TLegend(0.,0.,0.9,0.9);
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetLineColor(0);
  legend2->SetBorderSize(0);
  legend2->SetTextSize(0.07);
  legend2->SetHeader("#chi^{2}/#nu of Data-MC");
  legend2->AddEntry("","Data: LHC11c+d, EMC7","");
  legend2->AddEntry("","MC:","");
  legend2->AddEntry(gChi2L[0][0],"#gammaJ+JJ #it{p}^{#gamma,EMC}_{T}>3.5 GeV/#it{c}","PL");
  legend2->AddEntry(gChi2H[0][0],"#gammaJ+JJ #it{p}^{#gamma,EMC}_{T}>7 GeV/#it{c}","PL");
  
  legend2->Draw();
  
  TString fileName = Form("figures/%s/Comparison_Chi2_%s_JJDecLowHigh%s_PerSM_Ebin",
                          quantity.Data(),titleData.Data(),opt.Data());
  fileName+=fileFormat;
  gChi2graph->Print(fileName);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadRightMargin(0.01);

  TCanvas * gChi2Mingraph = new TCanvas(Form("cChi2Min_graph_%s",
                                             quantity.Data()),
                                        Form("Diff, %s",
                                             quantity.Data()),
                                        ncol*500,nrow*500);
  gChi2Mingraph->Divide(ncol,nrow);
  
  for(Int_t ipbin = 0; ipbin < nparam; ipbin++)
  {
    gChi2Mingraph->cd(ipbin+1);
    
    if(!gChi2LMin[ipbin]) continue;
    
    //gChi2LMin[ipbin]->GetHistogram()->SetTitleOffset(1.8,"Y");
   
    //gChi2LMin[ipbin]->GetHistogram()->SetYTitle("#mu_{1} (%)");
    //if(quantity.Contains("SM"   )) gChi2LMin[ipbin]->SetTitle(Form("SM %2.0f"            ,pbins[ipbin]));
    
    TH1F* haxis= new TH1F(Form("haxis%d",ipbin),Form("SM %d",ipbin),30,0.5,45.5);
    haxis->SetYTitle(gChi2LMin[ipbin]->GetHistogram()->GetYaxis()->GetTitle());
    haxis->SetXTitle(gChi2LMin[ipbin]->GetHistogram()->GetXaxis()->GetTitle());
    haxis->SetTitle(Form("SM %d",ipbin));
    haxis->SetMaximum(1.8);
    haxis->SetMinimum(0.0);
    
    //gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    
    haxis->Draw("axis");
    //printf("iebin %d, ipbin %d 0 %d %p\n",iebin,ipbin,0,h[0][ipbin]);
    gChi2LMin[ipbin]->Draw("PL");
    gChi2HMin[ipbin]->Draw("PL");
    
    //gChi2LMin[ipbin]->SetMaximum(1.7);
    //gChi2LMin[ipbin]->SetMinimum(-0.1);
    TLegend *lgChi2Min = 0;
    if(ipbin==3 ||ipbin==7) lgChi2Min= new TLegend(0.15,0.20,0.98,0.43);
    else                    lgChi2Min= new TLegend(0.15,0.75,0.98,0.98);
    lgChi2Min->SetFillColor(0);
    lgChi2Min->SetFillStyle(0);
    lgChi2Min->SetLineColor(0);
    lgChi2Min->SetBorderSize(0);
    lgChi2Min->SetTextSize(0.06);
    lgChi2Min->SetHeader(Form("SM %d",ipbin));
    lgChi2Min->Draw();

    Bool_t doFit = kFALSE;
    if(doFit)
    {
      gChi2HMin[ipbin]->Fit("pol0","QR","",16,30);   
      TF1* fun0 =  gChi2HMin[ipbin]->GetFunction("pol0");
      Float_t par0 = fun0->GetParameter(0);
      printf("sm %d parameter0 %2.2f\n",ipbin,par0);

      TF1* fun1 =  new TF1("plin1","[0]-x*[1]",2.5,16.);
      fun1->SetParameters(par0,-1.1);
      fun1->SetParLimits(0,par0*0.001,par0*5.1);
      fun1->SetParLimits(1,0.07,0.1);
      gChi2LMin[ipbin]->Fit("plin1","","",7.5,16.);   
      fun1->SetRange(0.5,16.0);
      fun1->Draw("same");
    }
    else
    {
      if(ipbin == 3 || ipbin == 7)
      {
        TF1* fun1 =  new TF1("plin1","2.2-x*0.09",0.5,16.);
        TF1* fun0 =  new TF1("plin0","0.8",16,45);
        fun1->Draw("same");
        fun0->Draw("same");
        fun0->SetLineColor(kViolet);
        lgChi2Min->AddEntry(fun0,"0.8%","L");
        lgChi2Min->AddEntry(fun1,"2.2% - #it{E}_{cluster} #times 0.09%","L");
      }

      else if(ipbin == 1 || ipbin == 2)
      {
        TF1* fun1 =  new TF1("plin1","1.5-x*0.08",0.5,13.);
        TF1* fun0 =  new TF1("plin0","0.5",13,45);
        fun1->Draw("same");
        fun0->Draw("same");
        fun0->SetLineColor(kViolet);
        lgChi2Min->AddEntry(fun0,"0.5%","L");
        lgChi2Min->AddEntry(fun1,"1.5% - #it{E}_{cluster} #times 0.08%","L");
      }
      else //if(ipbin == 1 || ipbin == 2)
      {
        TF1* fun1 =  new TF1("plin1","1.1-x*0.08",0.5,9);
        TF1* fun0 =  new TF1("plin0","0.4",9,45);
        fun1->Draw("same");
        fun0->Draw("same");
        fun0->SetLineColor(kViolet);
        lgChi2Min->AddEntry(fun0,"0.4%","L");
        lgChi2Min->AddEntry(fun1,"1.1% - #it{E}_{cluster} #times 0.08%","L");
      }
      
    }
    
  } // param
  
  gChi2Mingraph->cd(ncol*nrow);
  
  TLegend *legend = new TLegend(0.,0.,0.9,0.9);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.07);
  legend->SetHeader("Min of #chi^{2}/#nu of Data-MC");
  legend->AddEntry("","Data: LHC11c+d, EMC7","");
  legend->AddEntry("","MC:","");
  legend->AddEntry(gChi2LMin[0],"#gammaJ+JJ #it{p}^{#gamma,EMC}_{T}>3.5 GeV/#it{c}","PL");
  legend->AddEntry(gChi2HMin[0],"#gammaJ+JJ #it{p}^{#gamma,EMC}_{T}>7 GeV/#it{c}","PL");
  
  legend->Draw();
  
  fileName = Form("figures/%s/Comparison_Chi2Min_%s_JJDecLowHigh%s_PerSM_Ebin",
                  quantity.Data(),titleData.Data(),opt.Data());
  fileName+=fileFormat;
  gChi2Mingraph->Print(fileName);
}

