#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TImage.h>
#include <TPaveText.h>
#include <TLine.h>
#endif

void RatioPlot(TH1F **num,TH1F **den,TString t1,TString t2,TString opt,Double_t Low[],Double_t Up[]);
void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1);
void ResetOutOfRange(TH1F *histo, Int_t ipart, Double_t lowRange[3], Double_t upRange[3]);

void MakeCorrectedSpectraITSsaNsigma( TString period="LHC10d",
				      TString MCname="LHC10f6a",
				      Float_t DCAcut=7,
				      Int_t lowmult=-1,
				      Int_t upmult=-1
				      )    
{

  // parameters

  Bool_t MakeMCPlots=1;
  Bool_t MakeSystErr=0;
  Bool_t NormalizeData=1;
  Bool_t NormalizeToINEL=0;
  Bool_t SaveOnlyAsymm=0;
  Bool_t CorrectOnlyForEfficiency=0;
  Bool_t MarekNorm=0;
  Bool_t DCACorrection=1;
  Bool_t GFCorrection=1;

  //MC Histos
  TH1F *fHistPrimMCposBefEvSel[3];
  TH1F *fHistPrimMCnegBefEvSel[3];
  TH1F *fHistPrimMCpos[3];
  TH1F *fHistPrimMCneg[3];
  TH1F *hHistPosNSigma[3];
  TH1F *hHistNegNSigma[3];
  TH1F *hHistPosNSigmaMean[3];
  TH1F *hHistNegNSigmaMean[3];
  TH1F *hHistPosNSigmaPrim[3];
  TH1F *hHistNegNSigmaPrim[3];
  TH1F *hHistPosNSigmaPrimMean[3];
  TH1F *hHistNegNSigmaPrimMean[3];
  TH1F *hHistPosNSigmaPrimMC[3];
  TH1F *hHistNegNSigmaPrimMC[3];
  TH1F *hHistPosNSigmaPrimMCMean[3];
  TH1F *hHistNegNSigmaPrimMCMean[3];
  TH1F *hHistPosNSigmaMC[3];
  TH1F *hHistNegNSigmaMC[3];
  TH1F *hHistPosNSigmaMCMean[3];
  TH1F *hHistNegNSigmaMCMean[3];
  TH1F *hCorrFactPos[3];
  TH1F *hCorrFactNeg[3];
  TH1F *hCorrFactPosMean[3];
  TH1F *hCorrFactNegMean[3];
  TH1F *hEffPos[3];
  TH1F *hEffNeg[3];
  TH1F *hEffPosMean[3];
  TH1F *hEffNegMean[3];
  
  //DATA Histos
  TH1F *hHistPosNSigmaDATA[3];
  TH1F *hHistNegNSigmaDATA[3];
  TH1F *hHistPosNSigmaMeanDATA[3];
  TH1F *hHistNegNSigmaMeanDATA[3];
  TH1F *hITSsaPos[3];
  TH1F *hITSsaNeg[3];
  TH1F *hITSsaPosMean[3];
  TH1F *hITSsaNegMean[3];
  TH1F *hITSsaRawPos[3];
  TH1F *hITSsaRawNeg[3];
  TH1F *hITSsaRawPosMean[3];
  TH1F *hITSsaRawNegMean[3];
  
  //  Float_t nmin=1.5;
  Double_t Low[3]={0.10,0.2,0.3};
  Double_t Up[3]={0.4,0.5,0.7};
  Double_t lowRange[3]={0.10,0.2,0.3};
  Double_t upRange[3]={0.6,0.6,0.6};

  

  const Int_t nbinspt=22;
  Double_t xbins[nbinspt+1]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
    
  TString fname_MC=Form("./%s/AnalysisResults.root",MCname.Data());
  cout<<fname_MC.Data()<<endl;

  TString foutMC=Form("outSpectraMC_%s_%s_%.0fDCA.root",period.Data(),MCname.Data(),DCAcut);
  TFile *outMC=new TFile(foutMC.Data(),"recreate");
  outMC->Close();
  delete outMC;

  TString foutDATA;
  if(CorrectOnlyForEfficiency){
    foutDATA=Form("outSpectraData_%s_%s_%.0fDCA.root",period.Data(),MCname.Data(),DCAcut);
  }else{
    foutDATA=Form("outCorrSpectraData_%s_%s_%.0fDCA.root",period.Data(),MCname.Data(),DCAcut);
  }

  TString fnameDATA=Form("%s/AnalysisResults.root",period.Data());
  cout<<fnameDATA.Data()<<endl;
  TFile *outDATA=new TFile(foutDATA.Data(),"recreate");
  outDATA->Close();
  delete outDATA;
  

  TString lnameMC=Form("clistITSsaMult%ito%i",lowmult,upmult);
  Printf("\n------------------ READING MC  %s -------------------\n",lnameMC.Data());
  TFile *finMC = new TFile(fname_MC.Data());
  TDirectory *dirFileMC=(TDirectory*)finMC->Get("PWG2SpectraITSsa");
  TList * cOutputMC = (TList*)dirFileMC->Get(lnameMC.Data());
  //cOutputMC->Print();
  
  TH1F *fHistNEventsMC = (TH1F*)cOutputMC->FindObject("fHistNEvents");
  Float_t normMC=(Float_t)1/fHistNEventsMC->GetBinContent(fHistNEventsMC->FindBin(1.));
  printf("- NormalizationMC (NEvents): %.0f \n",fHistNEventsMC->GetBinContent(fHistNEventsMC->FindBin(1.))); 
  
  for(Int_t ipart=0;ipart<=2;ipart++){
    fHistPrimMCposBefEvSel[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistPrimMCposBefEvSel%i",ipart));//primaries generated (before event selection)
    fHistPrimMCnegBefEvSel[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistPrimMCnegBefEvSel%i",ipart));
    fHistPrimMCpos[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistPrimMCpos%i",ipart)); //primaries generated (after event selection)
    fHistPrimMCneg[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistPrimMCneg%i",ipart));
    hHistPosNSigma[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigma%i",ipart));//NSigma recontructed, no MCtruth
    hHistNegNSigma[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigma%i",ipart));
    hHistPosNSigmaPrimMC[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigmaPrimMC%i",ipart));//NSigma recontructed, MCtruth
    hHistNegNSigmaPrimMC[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigmaPrimMC%i",ipart));
    hHistPosNSigmaMC[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigmaMC%i",ipart));//NSigma recontructed, MCtruth
    hHistNegNSigmaMC[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigmaMC%i",ipart));
    hHistPosNSigmaPrim[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigmaPrim%i",ipart));//NSigma recontructed, truth
    hHistNegNSigmaPrim[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigmaPrim%i",ipart));
    hHistPosNSigmaMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigmaMean%i",ipart));//NSigmaMean recontructed, no MCtruth
    hHistNegNSigmaMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigmaMean%i",ipart));
    hHistPosNSigmaPrimMCMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigmaPrimMCMean%i",ipart));//NSigmaMean recontructed, no MCtruth
    hHistNegNSigmaPrimMCMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigmaPrimMCMean%i",ipart));
    hHistPosNSigmaMCMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigmaMCMean%i",ipart));//NSigmaMean recontructed, no MCtruth
    hHistNegNSigmaMCMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigmaMCMean%i",ipart));
    hHistPosNSigmaPrimMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistPosNSigmaPrimMean%i",ipart));//NSigmaMean recontructed, no MCtruth
    hHistNegNSigmaPrimMean[ipart]=(TH1F*)cOutputMC->FindObject(Form("hHistNegNSigmaPrimMean%i",ipart));
  }
    
  for(Int_t ipart=0;ipart<=2;ipart++){
    hCorrFactPos[ipart]=new TH1F(Form("hCorrFactPos%i",ipart),Form("hCorrFactPos%i",ipart),nbinspt,xbins);
    hCorrFactNeg[ipart]=new TH1F(Form("hCorrFactNeg%i",ipart),Form("hCorrFactNeg%i",ipart),nbinspt,xbins);
    hCorrFactPosMean[ipart]=new TH1F(Form("hCorrFactPosMean%i",ipart),Form("hCorrFactPosMean%i",ipart),nbinspt,xbins);
    hCorrFactNegMean[ipart]=new TH1F(Form("hCorrFactNegMean%i",ipart),Form("hCorrFactNegMean%i",ipart),nbinspt,xbins);
    hEffPos[ipart]=new TH1F(Form("hEffPos%i",ipart),Form("hEffPos%i",ipart),nbinspt,xbins);
    hEffNeg[ipart]=new TH1F(Form("hEffNeg%i",ipart),Form("hEffNeg%i",ipart),nbinspt,xbins);
    hEffPosMean[ipart]=new TH1F(Form("hEffPosMean%i",ipart),Form("hEffPosMean%i",ipart),nbinspt,xbins);
    hEffNegMean[ipart]=new TH1F(Form("hEffNegMean%i",ipart),Form("hEffNegMean%i",ipart),nbinspt,xbins);
  }
    
  for(Int_t ipart=0;ipart<=2;ipart++){     
    hHistPosNSigma[ipart]->Sumw2();
    hHistNegNSigma[ipart]->Sumw2();
    hHistPosNSigmaMean[ipart]->Sumw2();
    hHistNegNSigmaMean[ipart]->Sumw2();
    hHistPosNSigmaPrimMC[ipart]->Sumw2();
    hHistNegNSigmaPrimMC[ipart]->Sumw2();
    hHistPosNSigmaMC[ipart]->Sumw2();
    hHistNegNSigmaMC[ipart]->Sumw2();
    hHistPosNSigmaPrimMCMean[ipart]->Sumw2();
    hHistNegNSigmaPrimMCMean[ipart]->Sumw2();
    hHistPosNSigmaPrim[ipart]->Sumw2();
    hHistNegNSigmaPrim[ipart]->Sumw2();
    hHistPosNSigmaPrimMean[ipart]->Sumw2();
    hHistNegNSigmaPrimMean[ipart]->Sumw2();
    hHistPosNSigmaMCMean[ipart]->Sumw2();
    hHistNegNSigmaMCMean[ipart]->Sumw2();
    
    hCorrFactPos[ipart]->Divide(hHistPosNSigma[ipart],fHistPrimMCposBefEvSel[ipart],1,1,"B");
    hCorrFactNeg[ipart]->Divide(hHistNegNSigma[ipart],fHistPrimMCnegBefEvSel[ipart],1,1,"B");
    hCorrFactPosMean[ipart]->Divide(hHistPosNSigmaMean[ipart],fHistPrimMCposBefEvSel[ipart],1,1,"B");
    hCorrFactNegMean[ipart]->Divide(hHistNegNSigmaMean[ipart],fHistPrimMCnegBefEvSel[ipart],1,1,"B");
    
    // hCorrFactPos[ipart]->Divide(hHistPosNSigmaPrim[ipart],fHistPrimMCposBefEvSel[ipart],1,1,"B");
    // hCorrFactNeg[ipart]->Divide(hHistNegNSigmaPrim[ipart],fHistPrimMCnegBefEvSel[ipart],1,1,"B");
    // hCorrFactPosMean[ipart]->Divide(hHistPosNSigmaPrimMean[ipart],fHistPrimMCposBefEvSel[ipart],1,1,"B");
    // hCorrFactNegMean[ipart]->Divide(hHistNegNSigmaPrimMean[ipart],fHistPrimMCnegBefEvSel[ipart],1,1,"B");
    
    //hCorrFactPos[ipart]->Divide(hHistPosNSigmaPrimMC[ipart],fHistPrimMCposBefEvSel[ipart],1,1,"B");
    //hCorrFactNeg[ipart]->Divide(hHistNegNSigmaPrimMC[ipart],fHistPrimMCnegBefEvSel[ipart],1,1,"B");
    //hCorrFactPosMean[ipart]->Divide(hHistPosNSigmaPrimMCMean[ipart],fHistPrimMCposBefEvSel[ipart],1,1,"B");
    //hCorrFactNegMean[ipart]->Divide(hHistNegNSigmaPrimMCMean[ipart],fHistPrimMCnegBefEvSel[ipart],1,1,"B");
    
    hCorrFactPos[ipart]->SetTitle(Form("Correction Factor NSigma Pos%i",ipart));
    hCorrFactNeg[ipart]->SetTitle(Form("Correction Factor NSigma Neg%i",ipart));
    hCorrFactPosMean[ipart]->SetTitle(Form("Corr Fact NSigma Pos%i",ipart));
    hCorrFactNegMean[ipart]->SetTitle(Form("Corr Fact NSigma Neg%i",ipart));
    hCorrFactPos[ipart]->SetName(Form("hCorrFactPos%i",ipart));
    hCorrFactNeg[ipart]->SetName(Form("hCorrFactNeg%i",ipart));
    hCorrFactPosMean[ipart]->SetName(Form("hCorrFactPosMean%i",ipart));
    hCorrFactNegMean[ipart]->SetName(Form("hCorrFactNegMean%i",ipart));
    SetDrawAtt(21+ipart,ipart+4,1,ipart+4,1,hCorrFactPos[ipart]);
    SetDrawAtt(25+ipart,ipart+4,1,ipart+4,1,hCorrFactNeg[ipart]);
    SetDrawAtt(21+ipart,ipart+1,1,ipart+1,1,hCorrFactPosMean[ipart]);
    SetDrawAtt(25+ipart,ipart+1,1,ipart+1,1,hCorrFactNegMean[ipart]); 
    ResetOutOfRange(hCorrFactPos[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hCorrFactNeg[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hCorrFactPosMean[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hCorrFactNegMean[ipart],ipart,lowRange,upRange);
    
    hHistPosNSigmaPrimMC[ipart]->Sumw2();
    hHistNegNSigmaPrimMC[ipart]->Sumw2();
    hHistPosNSigmaPrimMCMean[ipart]->Sumw2();
    hHistNegNSigmaPrimMCMean[ipart]->Sumw2();
    hEffPos[ipart]->Divide(hHistPosNSigmaPrimMC[ipart],fHistPrimMCpos[ipart],1,1,"B");
    hEffNeg[ipart]->Divide(hHistNegNSigmaPrimMC[ipart],fHistPrimMCneg[ipart],1,1,"B");
    hEffPosMean[ipart]->Divide(hHistPosNSigmaPrimMCMean[ipart],fHistPrimMCpos[ipart],1,1,"B");
    hEffNegMean[ipart]->Divide(hHistNegNSigmaPrimMCMean[ipart],fHistPrimMCneg[ipart],1,1,"B");
    hEffPos[ipart]->SetTitle(Form("Efficiency NSigma Pos%i",ipart));
    hEffNeg[ipart]->SetTitle(Form("Efficiency NSigma Neg%i",ipart));
    hEffPosMean[ipart]->SetTitle(Form("Efficiency NSigma Asymm Pos%i",ipart));
    hEffNegMean[ipart]->SetTitle(Form("Efficiency NSigma Asymm Neg%i",ipart));
    hEffPos[ipart]->SetName(Form("hEffPos%i",ipart));
    hEffNeg[ipart]->SetName(Form("hEffNeg%i",ipart));
    hEffPosMean[ipart]->SetName(Form("hEffPosMean%i",ipart));
    hEffNegMean[ipart]->SetName(Form("hEffNegMean%i",ipart));
    
    SetDrawAtt(21+ipart,ipart+4,1,ipart+4,1,hEffPos[ipart]);
    SetDrawAtt(25+ipart,ipart+4,1,ipart+4,1,hEffNeg[ipart]);
    SetDrawAtt(21+ipart,ipart+1,1,ipart+1,1,hEffPosMean[ipart]);
    SetDrawAtt(25+ipart,ipart+1,1,ipart+1,1,hEffNegMean[ipart]);
    ResetOutOfRange(hEffPos[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hEffNeg[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hEffPosMean[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hEffNegMean[ipart],ipart,lowRange,upRange);
      
      
    hCorrFactPos[ipart]->SetMinimum(0.);
    hCorrFactNeg[ipart]->SetMinimum(0.);
    hCorrFactPosMean[ipart]->SetMinimum(0.);
    hCorrFactNegMean[ipart]->SetMinimum(0.);
    hEffPos[ipart]->SetMinimum(0.);
    hEffNeg[ipart]->SetMinimum(0.);
    hEffPosMean[ipart]->SetMinimum(0.);
    hEffNegMean[ipart]->SetMinimum(0.);
    hCorrFactPos[ipart]->SetMaximum(1.);
    hCorrFactNeg[ipart]->SetMaximum(1.);
    hCorrFactPosMean[ipart]->SetMaximum(1.);
    hCorrFactNegMean[ipart]->SetMaximum(1.);
    hEffPos[ipart]->SetMaximum(1.);
    hEffNeg[ipart]->SetMaximum(1.);
    hEffPosMean[ipart]->SetMaximum(1.);
    hEffNegMean[ipart]->SetMaximum(1.);
  }
    
  if(MakeMCPlots){
    TCanvas *cCorr=new TCanvas("cCorr","cCorr");
    cCorr->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0)hCorrFactPos[ipart]->Draw("P");
      else hCorrFactPos[ipart]->Draw("PSAME");
      hCorrFactNeg[ipart]->Draw("PSAME");
    }
    cCorr->BuildLegend();
    
    TCanvas *cCorrMean=new TCanvas("cCorrMean","cCorrMean");
    cCorrMean->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0){
	TH1F * hempty=(TH1F*)hCorrFactPosMean[ipart]->Clone();
	hempty->SetXTitle("P_{T} (GeV/c)");
	hempty->SetYTitle("Acc x #epsilon x Contamination");
	hempty->Reset("all");
	hempty->Draw();
	
      }
      hCorrFactPosMean[ipart]->Draw("PSAME");
      hCorrFactNegMean[ipart]->Draw("PSAME");
    }
    cCorrMean->BuildLegend();
    
    TCanvas *cCorrcfPos=new TCanvas("cCorrcfPos","cCorrcfPos");
    cCorrcfPos->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0)hCorrFactPos[ipart]->Draw("P");
      else hCorrFactPos[ipart]->Draw("PSAME");
      hCorrFactPosMean[ipart]->Draw("PSAME");
    }
    cCorrcfPos->BuildLegend();
    
    TCanvas *cCorrcfNeg=new TCanvas("cCorrcfNeg","cCorrcfNeg");
    cCorrcfNeg->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0)hCorrFactNeg[ipart]->Draw("P");
      else hCorrFactNeg[ipart]->Draw("PSAME");
      hCorrFactNegMean[ipart]->Draw("PSAME");
    }
    cCorrcfNeg->BuildLegend();
    
    TCanvas *cEff=new TCanvas("cEff","cEff");
    cEff->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0)hEffPos[ipart]->Draw("P");
      else hEffPos[ipart]->Draw("PSAME");
      hEffNeg[ipart]->Draw("PSAME");
    }
    cEff->BuildLegend();
    
    TCanvas *cEffMean=new TCanvas("cEffMean","cEffMean");
    cEffMean->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0)hEffPosMean[ipart]->Draw("P");
      else hEffPosMean[ipart]->Draw("PSAME");
      hEffNegMean[ipart]->Draw("PSAME");
    }
    cEffMean->BuildLegend();
    
    RatioPlot(hEffNegMean,hEffPosMean,"Neg","Pos","Efficiency",Low,Up);
    
    
    TCanvas *cEffcfPos=new TCanvas("cEffcfPos","cEffcfPos");
    cEffcfPos->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0)hEffPos[ipart]->Draw("P");
      else hEffPos[ipart]->Draw("PSAME");
      hEffPosMean[ipart]->Draw("PSAME");
    }
    cEffcfPos->BuildLegend();
    
    TCanvas *cEffcfNeg=new TCanvas("cEffcfNeg","cEffcfNeg");
    cEffcfNeg->SetGridy();
    for(Int_t ipart=0;ipart<=2;ipart++){
      if(ipart==0)hEffNeg[ipart]->Draw("P");
      else hEffNeg[ipart]->Draw("PSAME");
      hEffNegMean[ipart]->Draw("PSAME");
    }
    cEffcfNeg->BuildLegend();
  }
    
    
    
    
  for(Int_t ipart=0;ipart<=2;ipart++) {
    hHistPosNSigma[ipart]->Scale(normMC,"width");
    hHistNegNSigma[ipart]->Scale(normMC,"width");
    hHistPosNSigmaMean[ipart]->Scale(normMC,"width");
    hHistNegNSigmaMean[ipart]->Scale(normMC,"width");
  }
    

  //////////////////////////////Prim/All Sec/All SecStr/All
  TH1F *hPrimPos[3];
  TH1F *hPrimNeg[3];
  TH1F *hSecPos[3];
  TH1F *hSecNeg[3];
  TH1F *hSecStrPos[3];
  TH1F *hSecStrNeg[3];
  TH1F *hAllPos[3];
  TH1F *hAllNeg[3];
  TH1F *hAllSecPos[3];
  TH1F *hAllSecNeg[3];
    
  for(Int_t ipart=0;ipart<=2;ipart++){
    hPrimPos[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistPrimMCposReco%i",ipart));
    hPrimNeg[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistPrimMCnegReco%i",ipart));
    hSecPos[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistSecMatMCposReco%i",ipart));
    hSecNeg[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistSecMatMCnegReco%i",ipart));
    hSecStrPos[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistSecStrMCposReco%i",ipart));
    hSecStrNeg[ipart]=(TH1F*)cOutputMC->FindObject(Form("fHistSecStrMCnegReco%i",ipart));
    
    hAllSecPos[ipart]=(TH1F*)hSecPos[ipart]->Clone(Form("hAllSecPos%i",ipart));
    hAllSecPos[ipart]->Add(hSecStrPos[ipart]);
    hAllSecNeg[ipart]=(TH1F*)hSecNeg[ipart]->Clone(Form("hAllSecNeg%i",ipart));
    hAllSecNeg[ipart]->Add(hSecStrNeg[ipart]);
    hAllPos[ipart]=(TH1F*)hPrimPos[ipart]->Clone(Form("hAllPos%i",ipart));
    hAllPos[ipart]->Add(hSecPos[ipart]);
    hAllPos[ipart]->Add(hSecStrPos[ipart]);
    hAllNeg[ipart]=(TH1F*)hPrimNeg[ipart]->Clone(Form("hAllNeg%i",ipart));
    hAllNeg[ipart]->Add(hSecNeg[ipart]);
    hAllNeg[ipart]->Add(hSecStrNeg[ipart]);
    
    hPrimPos[ipart]->Divide(hAllPos[ipart]);
    hSecPos[ipart]->Divide(hAllPos[ipart]);
    hSecStrPos[ipart]->Divide(hAllPos[ipart]);
    hPrimNeg[ipart]->Divide(hAllNeg[ipart]);
    hSecNeg[ipart]->Divide(hAllNeg[ipart]);
    hSecStrNeg[ipart]->Divide(hAllNeg[ipart]);
    
    hAllSecPos[ipart]->Divide(hAllPos[ipart]);
    hAllSecNeg[ipart]->Divide(hAllNeg[ipart]);
    for(Int_t binnum=1;binnum<=hAllSecNeg[ipart]->GetNbinsX();binnum++){
      hAllSecPos[ipart]->SetBinContent(binnum,1-hAllSecPos[ipart]->GetBinContent(binnum));
      hAllSecNeg[ipart]->SetBinContent(binnum,1-hAllSecNeg[ipart]->GetBinContent(binnum));
    }
  }
  
  /////////////////////////////////////////////
  TFile *outMC2=new TFile(foutMC.Data(),"update");
  TList *lMC=new TList();
  lMC->SetOwner();
  lMC->SetName(Form("MC_Mult%ito%i",lowmult,upmult));
  for(Int_t ipart=0;ipart<=2;ipart++){
    lMC->Add(hCorrFactPos[ipart]);
    lMC->Add(hCorrFactNeg[ipart]);
    lMC->Add(hCorrFactPosMean[ipart]);
    lMC->Add(hCorrFactNegMean[ipart]);
    lMC->Add(hEffPos[ipart]);
    lMC->Add(hEffNeg[ipart]);
    lMC->Add(hEffPosMean[ipart]);
    lMC->Add(hEffNegMean[ipart]);
    lMC->Add(hHistPosNSigma[ipart]);
    lMC->Add(hHistNegNSigma[ipart]);
    lMC->Add(hHistPosNSigmaMean[ipart]);
    lMC->Add(hHistNegNSigmaMean[ipart]);
    lMC->Add(hPrimPos[ipart]);
    lMC->Add(hPrimNeg[ipart]);
    lMC->Add(hSecPos[ipart]);
    lMC->Add(hSecNeg[ipart]);
    lMC->Add(hSecStrPos[ipart]);
    lMC->Add(hSecStrNeg[ipart]);
    lMC->Add(hAllSecPos[ipart]);
    lMC->Add(hAllSecNeg[ipart]);
  }
  lMC->Write(Form("MC_Mult%ito%i",lowmult,upmult),1);
  outMC2->Close();
  delete outMC2;
  
    
  TString lnameDATA=Form("clistITSsaMult%ito%i",lowmult,upmult);
  Printf("\n\n------------------ READING DATA  %s -------------------\n",lnameDATA.Data());
  TFile *finDATA = new TFile(fnameDATA.Data());
  TDirectory *dirFileDATA=(TDirectory*)finDATA->Get("PWG2SpectraITSsa");
  TList *cOutputDATA = (TList*)dirFileDATA->Get(lnameDATA.Data());
  
  
  for(Int_t ipart=0;ipart<=2;ipart++){
    hHistPosNSigmaDATA[ipart]=(TH1F*)cOutputDATA->FindObject(Form("hHistPosNSigma%i",ipart));
    hHistNegNSigmaDATA[ipart]=(TH1F*)cOutputDATA->FindObject(Form("hHistNegNSigma%i",ipart));
    hHistPosNSigmaMeanDATA[ipart]=(TH1F*)cOutputDATA->FindObject(Form("hHistPosNSigmaMean%i",ipart));
    hHistNegNSigmaMeanDATA[ipart]=(TH1F*)cOutputDATA->FindObject(Form("hHistNegNSigmaMean%i",ipart)); 
    
    hITSsaPos[ipart]=(TH1F*)hHistPosNSigmaDATA[ipart]->Clone(Form("hITSsaPos%i",ipart));
    hITSsaNeg[ipart]=(TH1F*)hHistNegNSigmaDATA[ipart]->Clone(Form("hITSsaNeg%i",ipart));
    hITSsaPosMean[ipart]=(TH1F*)hHistPosNSigmaMeanDATA[ipart]->Clone(Form("hITSsaPosAsymm%i",ipart));
    hITSsaNegMean[ipart]=(TH1F*)hHistNegNSigmaMeanDATA[ipart]->Clone(Form("hITSsaNegAsymm%i",ipart));
    hITSsaRawPos[ipart]=(TH1F*)hHistPosNSigmaDATA[ipart]->Clone(Form("hITSsaRawPos%i",ipart));
    hITSsaRawNeg[ipart]=(TH1F*)hHistNegNSigmaDATA[ipart]->Clone(Form("hITSsaRawNeg%i",ipart));
    hITSsaRawPosMean[ipart]=(TH1F*)hHistPosNSigmaMeanDATA[ipart]->Clone(Form("hITSsaRawPosAsymm%i",ipart));
    hITSsaRawNegMean[ipart]=(TH1F*)hHistNegNSigmaMeanDATA[ipart]->Clone(Form("hITSsaRawNegAsymm%i",ipart));
    
    hITSsaPos[ipart]->Sumw2();
    hITSsaNeg[ipart]->Sumw2();
    hITSsaPosMean[ipart]->Sumw2();
    hITSsaNegMean[ipart]->Sumw2();
    hITSsaRawPos[ipart]->Sumw2();
    hITSsaRawNeg[ipart]->Sumw2();
    hITSsaRawPosMean[ipart]->Sumw2();
    hITSsaRawNegMean[ipart]->Sumw2();
    
    hITSsaPos[ipart]->Divide(hCorrFactPos[ipart]);
    hITSsaNeg[ipart]->Divide(hCorrFactNeg[ipart]);
    hITSsaPosMean[ipart]->Divide(hCorrFactPosMean[ipart]);
    hITSsaNegMean[ipart]->Divide(hCorrFactNegMean[ipart]);
    
    SetDrawAtt(21+ipart,ipart+4,1,ipart+4,1,hITSsaPos[ipart]);
    SetDrawAtt(25+ipart,ipart+4,1,ipart+4,1,hITSsaNeg[ipart]);
    SetDrawAtt(21+ipart,ipart+1,1,ipart+1,1,hITSsaPosMean[ipart]);
    SetDrawAtt(25+ipart,ipart+1,1,ipart+1,1,hITSsaNegMean[ipart]); 
    SetDrawAtt(21+ipart,ipart+4,1,ipart+4,1,hITSsaRawPos[ipart]);
    SetDrawAtt(25+ipart,ipart+4,1,ipart+4,1,hITSsaRawNeg[ipart]);
    SetDrawAtt(21+ipart,ipart+1,1,ipart+1,1,hITSsaRawPosMean[ipart]);
    SetDrawAtt(25+ipart,ipart+1,1,ipart+1,1,hITSsaRawNegMean[ipart]); 
  }
  
  if(!CorrectOnlyForEfficiency){
    
    
    if(DCACorrection){
      printf("- Correction from DCA\n");
      
      // Correction factor based on fit to DCA distr for P and Pibar
      TString fnameDCAPosP=Form("./DCACorr%s_%s_%.0fDCA_PosP_TFraction.root",period.Data(),MCname.Data(),DCAcut);
      TFile *fDCAPosP= new TFile(fnameDCAPosP.Data());
      TH1F * hDCAPosP= (TH1F*)fDCAPosP->Get("hPrimAllDATAMC");
      //TH1F * hDCAPosP= (TH1F*)fDCAPosP->Get("hPrimAllDATA");
      TString fnameDCANegP=Form("./DCACorr%s_%s_%.0fDCA_NegP_TFraction.root",period.Data(),MCname.Data(),DCAcut);
      TFile *fDCANegP= new TFile(fnameDCANegP.Data());
      TH1F * hDCANegP= (TH1F*)fDCANegP->Get("hPrimAllDATAMC");
      //TH1F * hDCANegP= (TH1F*)fDCANegP->Get("hPrimAllDATA");
      
      for(Int_t pbin=0;pbin<=hDCAPosP->GetNbinsX();pbin++)hDCAPosP->SetBinError(pbin,0);
      for(Int_t pbin=0;pbin<=hDCANegP->GetNbinsX();pbin++)hDCANegP->SetBinError(pbin,0);
      
      hITSsaPos[2]->Multiply(hDCAPosP);
      hITSsaNeg[2]->Multiply(hDCANegP);
      hITSsaPosMean[2]->Multiply(hDCAPosP);
      hITSsaNegMean[2]->Multiply(hDCANegP);
      
      hDCAPosP->Draw("");
      hDCANegP->Draw("same"); 

    }
    else{
      Printf("NO DCA correction for protons");
    }
    
    if(GFCorrection){
      printf("- Geant/FLUKA correction\n"); 
      //GeantFluka correction for Pi and K
      TString fnameGeanFlukaPi="$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/ITSsa/RootFilesGeantFlukaCorrection/correctionForCrossSection.211.root";
      TFile *fGeanFlukaPi= new TFile(fnameGeanFlukaPi.Data());
      TH1F *hGeantFlukaPiPos=(TH1F*)fGeanFlukaPi->Get("gHistCorrectionForCrossSectionParticles");
      TH1F *hGeantFlukaPiNeg=(TH1F*)fGeanFlukaPi->Get("gHistCorrectionForCrossSectionAntiParticles");
      for(Int_t binPi=0;binPi<=hITSsaPos[0]->GetNbinsX();binPi++){
	Float_t FlukaCorrPiPos=hGeantFlukaPiPos->GetBinContent(hGeantFlukaPiPos->FindBin(hITSsaPos[0]->GetBinCenter(binPi)));
	Float_t FlukaCorrPiNeg=hGeantFlukaPiNeg->GetBinContent(hGeantFlukaPiNeg->FindBin(hITSsaNeg[0]->GetBinCenter(binPi)));
	hITSsaPos[0]->SetBinContent(binPi,hITSsaPos[0]->GetBinContent(binPi)*FlukaCorrPiPos);
	hITSsaPos[0]->SetBinError(binPi,hITSsaPos[0]->GetBinError(binPi)*FlukaCorrPiPos);
	hITSsaNeg[0]->SetBinContent(binPi,hITSsaNeg[0]->GetBinContent(binPi)*FlukaCorrPiNeg);
	hITSsaNeg[0]->SetBinError(binPi,hITSsaNeg[0]->GetBinError(binPi)*FlukaCorrPiNeg);
	hITSsaPosMean[0]->SetBinContent(binPi,hITSsaPosMean[0]->GetBinContent(binPi)*FlukaCorrPiPos);
	hITSsaPosMean[0]->SetBinError(binPi,hITSsaPosMean[0]->GetBinError(binPi)*FlukaCorrPiPos);
	hITSsaNegMean[0]->SetBinContent(binPi,hITSsaNegMean[0]->GetBinContent(binPi)*FlukaCorrPiNeg);
	hITSsaNegMean[0]->SetBinError(binPi,hITSsaNegMean[0]->GetBinError(binPi)*FlukaCorrPiNeg);
      }
      
      TString fnameGeanFlukaK="$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/ITSsa/RootFilesGeantFlukaCorrection/correctionForCrossSection.321.root";
      TFile *fGeanFlukaK= new TFile(fnameGeanFlukaK.Data());
      TH1F *hGeantFlukaKPos=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionParticles");
      TH1F *hGeantFlukaKNeg=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionAntiParticles");
      for(Int_t binK=0;binK<=hITSsaPos[1]->GetNbinsX();binK++){
	Float_t FlukaCorrKPos=hGeantFlukaKPos->GetBinContent(hGeantFlukaKPos->FindBin(hITSsaPos[1]->GetBinCenter(binK)));
	Float_t FlukaCorrKNeg=hGeantFlukaKNeg->GetBinContent(hGeantFlukaKNeg->FindBin(hITSsaNeg[1]->GetBinCenter(binK)));
	hITSsaPos[1]->SetBinContent(binK,hITSsaPos[1]->GetBinContent(binK)*FlukaCorrKPos);
	hITSsaNeg[1]->SetBinContent(binK,hITSsaNeg[1]->GetBinContent(binK)*FlukaCorrKNeg);
	hITSsaPos[1]->SetBinError(binK,hITSsaPos[1]->GetBinError(binK)*FlukaCorrKPos);
	hITSsaNeg[1]->SetBinError(binK,hITSsaNeg[1]->GetBinError(binK)*FlukaCorrKNeg);
	hITSsaPosMean[1]->SetBinContent(binK,hITSsaPosMean[1]->GetBinContent(binK)*FlukaCorrKPos);
	hITSsaNegMean[1]->SetBinContent(binK,hITSsaNegMean[1]->GetBinContent(binK)*FlukaCorrKNeg);
	hITSsaPosMean[1]->SetBinError(binK,hITSsaPosMean[1]->GetBinError(binK)*FlukaCorrKPos);
	hITSsaNegMean[1]->SetBinError(binK,hITSsaNegMean[1]->GetBinError(binK)*FlukaCorrKNeg);
      }
      
      //Geant Fluka for P
      // ITS specific file for protons/antiprotons
      Int_t kNCharge=2;
      Int_t kPos=0;
      Int_t kNeg=1;
      TFile* fITS = new TFile ("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/ITSsa/RootFilesGeantFlukaCorrection/correctionForCrossSectionITS_20100719.root");
      // TH2D * hCorrFlukaITS[kNCharge];
      TH2D * hCorrFlukaITS[2];
      hCorrFlukaITS[kPos] = (TH2D*)fITS->Get("gHistCorrectionForCrossSectionProtons");
      hCorrFlukaITS[kNeg] = (TH2D*)fITS->Get("gHistCorrectionForCrossSectionAntiProtons");
      for(Int_t icharge = 0; icharge < kNCharge; icharge++){
	Int_t nbins = hITSsaPos[2]->GetNbinsX();
	Int_t nbinsy=hCorrFlukaITS[icharge]->GetNbinsY();
	//for(Int_t binP=0;binP<=hITSsaPos[2]->GetNbinsX();binP++){
	for(Int_t ibin = 0; ibin < nbins; ibin++){
	  Float_t pt = hITSsaPos[2]->GetBinCenter(ibin);
	  Float_t minPtCorrection = hCorrFlukaITS[icharge]->GetYaxis()->GetBinLowEdge(1);
	  Float_t maxPtCorrection = hCorrFlukaITS[icharge]->GetYaxis()->GetBinLowEdge(nbinsy+1);
	  if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
	  if (pt > maxPtCorrection) pt = maxPtCorrection;
	  Float_t correction = hCorrFlukaITS[icharge]->GetBinContent(1,hCorrFlukaITS[icharge]->GetYaxis()->FindBin(pt));
	  cout<<correction<<"     charge "<<icharge<<endl;
	  if(icharge==0){
	    if (correction != 0) {// If the bin is empty this is a  0
	      hITSsaPos[2]->SetBinContent(ibin,hITSsaPos[2]->GetBinContent(ibin)*correction);
	      hITSsaPos[2]->SetBinError(ibin,hITSsaPos[2]->GetBinError  (ibin)*correction);
	      hITSsaPosMean[2]->SetBinContent(ibin,hITSsaPosMean[2]->GetBinContent(ibin)*correction);
	      hITSsaPosMean[2]->SetBinError(ibin,hITSsaPosMean[2]->GetBinError  (ibin)*correction);
	      
	    }else if (hITSsaPos[2]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
	      cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries, ITS, " << endl;
	      cout << " Bin content: " << hITSsaPos[2]->GetBinContent(ibin)  << endl;
	    }
	  }
	  if(icharge==1){
	    if (correction != 0) {// If the bin is empty this is a  0
	      hITSsaNeg[2]->SetBinContent(ibin,hITSsaNeg[2]->GetBinContent(ibin)*correction);
	      hITSsaNeg[2]->SetBinError(ibin,hITSsaNeg[2]->GetBinError  (ibin)*correction);
	      hITSsaNegMean[2]->SetBinContent(ibin,hITSsaNegMean[2]->GetBinContent(ibin)*correction);
	      hITSsaNegMean[2]->SetBinError(ibin,hITSsaNegMean[2]->GetBinError  (ibin)*correction);
	    }else if (hITSsaNeg[2]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
	      cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for Antiprotons secondaries, ITS, " << endl;
	      cout << " Bin content: " << hITSsaNeg[2]->GetBinContent(ibin)  << endl;
	    }
	  }
	  
	}	
      }
    }else{
      Printf("SPECTRA ARE NOT CORRECTED FOR GEANT/FLUKA");
      
    }
    
  }else{
    Printf("SPECTRA ARE NOT CORRECTED FOR DCA FIT + GEANT/FLUKA");
  }

  Float_t norm;
  TH1F *fHistNEvents =0x0;
  if(MarekNorm){
    Printf("Normalization to the number of events with good vertex");
    TFile *finMAREK = new TFile(fnameDATA.Data());
    TDirectory *dirFileMAREK=(TDirectory*)finMAREK->Get("PWG2SpectraITSTPC");
    if(!dirFileMAREK)Printf("!dirFileMAREK");
    TString lnameMAREK=Form("outputlow%dup%dHI0",lowmult,upmult);
    TList *cOutputMAREK = (TList*)dirFileMAREK->Get(lnameMAREK.Data());
    if(!cOutputMAREK)Printf("!cOutputMAREK");
    TH1F *StatsHistMAREK=(TH1F*)cOutputMAREK->FindObject("StatsHist");
    norm=(Float_t)1/StatsHistMAREK->GetBinContent(3);
  }
  else{
    Printf("Normalization to the number of events after phys sel");
    fHistNEvents = (TH1F*)cOutputDATA->FindObject("fHistNEvents");
    norm=(Float_t)1/fHistNEvents->GetBinContent(fHistNEvents->FindBin(1.));
  }
  if(NormalizeData)printf("- Normalization (NEvents): %.0f \n",1/norm); 
  else printf("- Normalization (NEvents): %.0f !!!!!!!SPECTRA ARE NOT NORMALIZED!!!!!!\n",fHistNEvents->GetBinContent(fHistNEvents->FindBin(1.))); 
  for(Int_t ipart=0;ipart<=2;ipart++){
    if(NormalizeData){
      hITSsaPos[ipart]->Scale(norm,"width");
      hITSsaNeg[ipart]->Scale(norm,"width");
      hITSsaPosMean[ipart]->Scale(norm,"width");
      hITSsaNegMean[ipart]->Scale(norm,"width");
      hITSsaRawPos[ipart]->Scale(norm,"width");
      hITSsaRawNeg[ipart]->Scale(norm,"width");
      hITSsaRawPosMean[ipart]->Scale(norm,"width");
      hITSsaRawNegMean[ipart]->Scale(norm,"width");
    }
    else{
      hITSsaPos[ipart]->Scale(1,"width");
      hITSsaNeg[ipart]->Scale(1,"width");
      hITSsaPosMean[ipart]->Scale(1,"width");
      hITSsaNegMean[ipart]->Scale(1,"width");
      hITSsaRawPos[ipart]->Scale(1,"width");
      hITSsaRawNeg[ipart]->Scale(1,"width");
      hITSsaRawPosMean[ipart]->Scale(1,"width");
      hITSsaRawNegMean[ipart]->Scale(1,"width");
    }
    if(NormalizeToINEL){
      Float_t INEL=0.86;
      printf("Normalization to INEL : %.2f \n",INEL);
      hITSsaPos[ipart]->Scale(INEL,"");
      hITSsaNeg[ipart]->Scale(INEL,"");
      hITSsaPosMean[ipart]->Scale(INEL,"");
      hITSsaNegMean[ipart]->Scale(INEL,"");
      hITSsaRawPos[ipart]->Scale(INEL,"");
      hITSsaRawNeg[ipart]->Scale(INEL,"");
      hITSsaRawPosMean[ipart]->Scale(INEL,"");
      hITSsaRawNegMean[ipart]->Scale(INEL,"");
    }
    else printf("Normalization to INEL ignored \n");
    
    hITSsaPos[ipart]->SetMinimum(0.00001);
    hITSsaNeg[ipart]->SetMinimum(0.00001);
    hITSsaPosMean[ipart]->SetMinimum(0.00001);
    hITSsaNegMean[ipart]->SetMinimum(0.00001);
    hITSsaRawPos[ipart]->SetMinimum(0.00001);
    hITSsaRawNeg[ipart]->SetMinimum(0.00001);
    hITSsaRawPosMean[ipart]->SetMinimum(0.00001);
    hITSsaRawNegMean[ipart]->SetMinimum(0.00001);
    hITSsaPos[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaNeg[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaPosMean[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaNegMean[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaRawPos[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaRawNeg[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaRawPosMean[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaRawNegMean[ipart]->GetXaxis()->SetTitle("Pt [GeV/c]");
    hITSsaPos[ipart]->SetTitle(Form("ITSsa NSigma Symm Pos%i",ipart));
    hITSsaNeg[ipart]->SetTitle(Form("ITSsa NSigma Symm Neg%i",ipart));
    hITSsaPosMean[ipart]->SetTitle(Form("ITSsa NSigma Asymm Pos%i",ipart));
    hITSsaNegMean[ipart]->SetTitle(Form("ITSsa NSigma Asymm Neg%i",ipart));
    hITSsaRawPos[ipart]->SetTitle(Form("ITSsaRaw NSigma Symm Pos%i",ipart));
    hITSsaRawNeg[ipart]->SetTitle(Form("ITSsaRaw NSigma Symm Neg%i",ipart));
    hITSsaRawPosMean[ipart]->SetTitle(Form("ITSsaRaw NSigma Asymm Pos%i",ipart));
    hITSsaRawNegMean[ipart]->SetTitle(Form("ITSsaRaw NSigma Asymm Neg%i",ipart));
      

    ResetOutOfRange(hITSsaPos[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hITSsaNeg[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hITSsaPosMean[ipart],ipart,lowRange,upRange);
    ResetOutOfRange(hITSsaNegMean[ipart],ipart,lowRange,upRange);

  }
  TH1F *hSystematicsITSsaPos[3];
  TH1F *hSystematicsITSsaNeg[3];
  TH1F *hSystematicsITSsaPosMean[3];
  TH1F *hSystematicsITSsaNegMean[3];
  
  if(MakeSystErr){
    printf("- making Systematic Error\n");
    //preliminary systematic error
    TFile *fsys = new TFile("RootFilesSystError/ITSsa-systematics-20101014.root");
    //fsys->Print("all");
    //cout<<fsys->GetSize()<<endl;
    TH1F* hSystPos[3];
    TH1F* hSystNeg[3];
    hSystPos[0]=(TH1F*)fsys->Get("hSystTotPosPion");
    hSystNeg[0]=(TH1F*)fsys->Get("hSystTotNegPion");
    hSystPos[1]=(TH1F*)fsys->Get("hSystTotPosKaon");
    hSystNeg[1]=(TH1F*)fsys->Get("hSystTotNegKaon");
    hSystPos[2]=(TH1F*)fsys->Get("hSystTotPosProton");
    hSystNeg[2]=(TH1F*)fsys->Get("hSystTotNegProton");

    for(Int_t ipart=0;ipart<=2;ipart++){
      hSystematicsITSsaPos[ipart]=(TH1F*)hITSsaPos[ipart]->Clone(Form("hSystematicsITSsaPos%i",ipart));  
      hSystematicsITSsaNeg[ipart]=(TH1F*)hITSsaNeg[ipart]->Clone(Form("hSystematicsITSsaNeg%i",ipart));  
      hSystematicsITSsaPosMean[ipart]=(TH1F*)hITSsaPosMean[ipart]->Clone(Form("hSystematicsITSsaPosAsymm%i",ipart));  
      hSystematicsITSsaNegMean[ipart]=(TH1F*)hITSsaNegMean[ipart]->Clone(Form("hSystematicsITSsaNegAsymm%i",ipart));  
    
      for(Int_t ibin=0;ibin<=hSystematicsITSsaPos[ipart]->GetNbinsX();ibin++){
	Float_t syserr=hSystPos[ipart]->GetBinContent(hSystPos[ipart]->FindBin(hSystematicsITSsaPos[ipart]->GetBinCenter(ibin)));
	hSystematicsITSsaPos[ipart]->SetBinError(ibin,syserr*hSystematicsITSsaPos[ipart]->GetBinContent(ibin));
      	hSystematicsITSsaPosMean[ipart]->SetBinError(ibin,syserr*hSystematicsITSsaPosMean[ipart]->GetBinContent(ibin));
      }
      for(Int_t ibin=0;ibin<=hSystematicsITSsaNeg[ipart]->GetNbinsX();ibin++){
	Float_t syserr=hSystNeg[ipart]->GetBinContent(hSystNeg[ipart]->FindBin(hSystematicsITSsaNeg[ipart]->GetBinCenter(ibin)));
	hSystematicsITSsaNeg[ipart]->SetBinError(ibin,syserr*hSystematicsITSsaNeg[ipart]->GetBinContent(ibin));
	hSystematicsITSsaNegMean[ipart]->SetBinError(ibin,syserr*hSystematicsITSsaNegMean[ipart]->GetBinContent(ibin));
      }
    }
  }
    
  TCanvas *cSpectraSymm=new TCanvas("cSpectraSymm","cSpectraSymm");
  cSpectraSymm->SetGridy();
  for(Int_t ipart=0;ipart<=2;ipart++){
    if(ipart==0)hITSsaPos[ipart]->DrawCopy("P");
    else hITSsaPos[ipart]->DrawCopy("PSAME");
    hITSsaNeg[ipart]->DrawCopy("PSAME");
  }
  cSpectraSymm->BuildLegend();
  
  TCanvas *cSpectraAsymm=new TCanvas("cSpectraAsymm","cSpectraAsymm");
  cSpectraAsymm->SetGridy();
  for(Int_t ipart=0;ipart<=2;ipart++){
    if(ipart==0)hITSsaPosMean[ipart]->DrawCopy("P");
    else hITSsaPosMean[ipart]->DrawCopy("PSAME");
    hITSsaNegMean[ipart]->DrawCopy("PSAME");
  }
  cSpectraAsymm->BuildLegend();
  

  
  ///////////////////////////  SCALING SPECTRA FOR TRACKING EFFICIENCY   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Float_t ScalingFactor=1.015;
  Printf("SCALING FACTOR: %f",ScalingFactor);
  for(Int_t ipart=0;ipart<=2;ipart++){
    hITSsaPosMean[ipart]->Scale(ScalingFactor);
    Printf("Scaling %s by %f",hITSsaPosMean[ipart]->GetName(),ScalingFactor);
    hITSsaNegMean[ipart]->Scale(ScalingFactor);
    Printf("Scaling %s by %f",hITSsaNegMean[ipart]->GetName(),ScalingFactor);
    hITSsaPos[ipart]->Scale(ScalingFactor);
    Printf("Scaling %s by %f",hITSsaPos[ipart]->GetName(),ScalingFactor);
    hITSsaNeg[ipart]->Scale(ScalingFactor);
    Printf("Scaling %s by %f",hITSsaNeg[ipart]->GetName(),ScalingFactor);
  }
  
  TFile *outDATA2=new TFile(foutDATA.Data(),"update");
  TList *lDATA=new TList();
  lDATA->SetOwner();
  lDATA->SetName(Form("DATA_Mult%ito%i",lowmult,upmult));
  if(SaveOnlyAsymm){
    for(Int_t ipart=0;ipart<=2;ipart++)  lDATA->Add(hITSsaPosMean[ipart]);
    for(Int_t ipart=0;ipart<=2;ipart++)  lDATA->Add(hITSsaNegMean[ipart]);
  }else{
    for(Int_t ipart=0;ipart<=2;ipart++){
      lDATA->Add(fHistPrimMCposBefEvSel[ipart]);
      lDATA->Add(fHistPrimMCnegBefEvSel[ipart]);
      lDATA->Add(hHistPosNSigmaPrimMean[ipart]);
      lDATA->Add(hHistNegNSigmaPrimMean[ipart]);
      lDATA->Add(hHistPosNSigmaPrimMCMean[ipart]);
      lDATA->Add(hHistNegNSigmaPrimMCMean[ipart]);
      lDATA->Add(hITSsaPos[ipart]);
      lDATA->Add(hITSsaNeg[ipart]);
      lDATA->Add(hITSsaPosMean[ipart]);
      lDATA->Add(hITSsaNegMean[ipart]);
      lDATA->Add(hITSsaRawPos[ipart]);
      lDATA->Add(hITSsaRawNeg[ipart]);
      lDATA->Add(hITSsaRawPosMean[ipart]);
      lDATA->Add(hITSsaRawNegMean[ipart]);
    }
    if(MakeSystErr){
      for(Int_t ipart=0;ipart<=2;ipart++){
	lDATA->Add(hSystematicsITSsaPos[ipart]);
	lDATA->Add(hSystematicsITSsaNeg[ipart]);
	lDATA->Add( hSystematicsITSsaPosMean[ipart]);
	lDATA->Add( hSystematicsITSsaNegMean[ipart]);
      }
    }
  }
  lDATA->Write(Form("DATA_Mult%ito%i",lowmult,upmult),1);
  outDATA2->Close();
  delete outDATA2;
  
  
  RatioPlot(hITSsaNeg,hITSsaPos,"Neg","Pos","Symm NSigma",Low,Up);
  RatioPlot(hITSsaNegMean,hITSsaPosMean,"Neg","Pos","Assymm NSigma",Low,Up);
  RatioPlot(hITSsaPos,hITSsaPosMean,"Symm","Asymm","Positive",Low,Up);
  RatioPlot(hITSsaNeg,hITSsaNegMean,"Symm","Asymm","Negative",Low,Up);
  
  
} //end main


void ResetOutOfRange(TH1F *histo, Int_t ipart, Double_t lowRange[3], Double_t upRange[3]){
  // set to -1 the bin contents out of the selcted pt range
  for(Int_t ibin=0;ibin<=histo->GetNbinsX();ibin++){
    if(histo->GetBinCenter(ibin)<lowRange[ipart] || histo->GetBinCenter(ibin)>upRange[ipart]) histo->SetBinContent(ibin,-1);
  }
}

void RatioPlot(TH1F **num,TH1F **den,TString t1,TString t2,TString opt,Double_t Low[],Double_t Up[])
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  TH1F *Num[3];
  TString title=Form("%s/%s - %s",t1.Data(),t2.Data(),opt.Data());
  TCanvas *c=new TCanvas(title.Data(),title.Data());
  c->Divide(3,2);
  for(Int_t ipart=0;ipart<=2;ipart++){
    num[ipart]->GetXaxis()->SetRangeUser(Low[ipart],Up[ipart]-0.00001);
    num[ipart]->SetMarkerColor(2);
    num[ipart]->SetLineColor(2);
    num[ipart]->SetStats(kFALSE);
    den[ipart]->GetXaxis()->SetRangeUser(Low[ipart],Up[ipart]-0.00001);
    den[ipart]->SetMarkerColor(1);
    den[ipart]->SetLineColor(1);
    den[ipart]->SetStats(kFALSE);
    c->cd(ipart+1);
    gPad->SetGridy();
    //gPad->SetLogy();
    num[ipart]->Draw("");
    den[ipart]->Draw("same");
    TPaveText *tpave=new TPaveText(0.1,0.2,0.5,0.29,"brNDC");
    tpave->SetBorderSize(0);
    tpave->SetFillStyle(0);
    tpave->SetFillColor(0);
    tpave->SetTextColor(2);
    tpave->SetTextSize(.04);
    TText *txt1=tpave->AddText(t1.Data());
    txt1->SetTextFont(62);
    txt1->SetTextColor(2);
    TText *txt2=tpave->AddText(t2.Data());
    txt2->SetTextFont(62);
    txt2->SetTextColor(1);
    tpave->Draw();
    Num[ipart]=(TH1F*)num[ipart]->Clone(title.Data());
    Num[ipart]->SetTitle(title.Data());
    Num[ipart]->Divide(den[ipart]);
    Num[ipart]->GetXaxis()->SetRangeUser(Low[ipart],Up[ipart]-0.00001);
    Num[ipart]->SetMarkerColor(4);
    Num[ipart]->SetLineColor(4);
    Num[ipart]->SetMaximum(1.2);
    Num[ipart]->SetMinimum(.8);
    c->cd(ipart+4);
    gPad->SetGridy();
    Num[ipart]->SetStats(kFALSE);
    Num[ipart]->Draw();
    TPaveText *tpave_1=new TPaveText(0.3,0.7,0.7,0.99,"brNDC");
    tpave_1->SetBorderSize(0);
    tpave_1->SetFillStyle(0);
    tpave_1->SetFillColor(0);
    tpave_1->SetTextColor(2);
    tpave_1->SetTextSize(.04);
    TText *txt1_1=tpave_1->AddText(title.Data());
    txt1_1->SetTextFont(62);
    tpave_1->Draw();
    // TLine line = TLine(Low[ipart],1,Up[ipart],1);
    // line.SetLineColor(2);
    // line.SetLineStyle(2);
    // line.SetLineWidth(2);
    // line.DrawClone("same");
  }
}




 void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1)
{ 
  h1->SetMarkerStyle(markerstyle);
  h1->SetMarkerColor(markercolor);
  h1->SetMarkerSize(markersize);
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
}


