const Int_t numberOfCentralityBins = 8;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};

const Int_t gRebin = 1;

void drawCorrelationFunctionPsiSummary(TString histoName = "PN",
				       const char* lhcPeriod = "LHC11h",
				       Int_t gTrainID = 222,			      
				       Int_t gCentrality = 1,
				       Double_t psiMin = -0.5, 
				       Double_t psiMax = 3.5) {
  // Macro that draws the fit results for the 
  // correlation functions from the balance function analysis
  // Author: m.weber@cern.ch

  gROOT->LoadMacro("~/SetPlotStyle.C");
  SetPlotStyle();
  gStyle->SetPalette(1,0);

  //Load the PWG2ebye library
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGCFebye.so");

  // pt bins
  // this could also be retrieved directly from AliBalancePsi
  //const Int_t kNPtBins = 16;
  //Double_t ptBins[kNPtBins+1] = {0.2,0.6,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.,12.,15.,20.};
  //const Int_t kNPtBins = 5;
  //Double_t ptBins[kNPtBins+1] = {0.6,1.0,1.5,2.0,4.0,20.0};
  const Int_t kNPtBins = 4;
  Double_t ptBins[kNPtBins+1] = {1.0,2.0,3.0,4.0,8.0};
  //const Int_t kNPtBins = 1;
  //Double_t ptBins[kNPtBins+1] = {1.0,2.0};

  Double_t pt[kNPtBins*kNPtBins];
  Double_t ptE[kNPtBins*kNPtBins];
  for(Int_t i = 0; i < kNPtBins; i++){
    for(Int_t j = 0; j < kNPtBins; j++){
      pt[i*kNPtBins+j] = 10*i + (ptBins[j]+ptBins[j+1])/2.;
      ptE[i*kNPtBins+j] = 0.2;
    }
  }


  // Canvases
  TCanvas *cQA[kNPtBins][kNPtBins];
  for(Int_t i = 0; i < kNPtBins; i++){
    for(Int_t j = 0; j <= i; j++){
      cQA[i][j] = new TCanvas(Form("cQA%d%d",i,j),Form("Fitting QA for bin %d %d",i,j),1200,900);
      cQA[i][j]->Divide(3,3);
    }
  }


  // Loop over pt bins
  Double_t ptTriggerMin = 0.0;
  Double_t ptTriggerMax = 0.;
  Double_t ptAssociatedMin = 0.0;
  Double_t ptAssociatedMax = 0.0;
  TString inFileName = "";
  
  //Fit Parameters
  Double_t p[17][kNPtBins*kNPtBins];
  Double_t pE[17][kNPtBins*kNPtBins];
  TString pNames[17] = {
    "Normalization",
    "NearSideN",
    "NearSideSigmaDeltaEta",
    "NearSideSigmaDeltaPhi",
    "NearSideSigmaExponent",
    "AwaySideN",
    "AwaySideSigmaDeltaPhi",
    "AwaySideSigmaExponent",
    "LongRidgeN",
    "LongRidgeSigma",
    "LongRidgeExponent",
    "Wing",
    "FlowN",
    "FlowV1",
    "FlowV2",
    "FlowV3",
    "FlowV4"
  }
  
  for(Int_t iPar = 0; iPar < 17; iPar++){
    for(Int_t i = 0; i < kNPtBins; i++){
      for(Int_t j = 0; j < kNPtBins; j++){
	p[iPar][i*kNPtBins+j] = -1.;
	pE[iPar][i*kNPtBins+j] = 0.;
      }
    }
  }

  TFile *inFile  = NULL; 
  TH2D *hTMPData = NULL;
  TH2D *hTMPRes  = NULL;
  TH2D *hTMPFit  = NULL;
  TF2 *fFit  = NULL;

  
  for(Int_t i = 0; i < kNPtBins; i++){
    for(Int_t j = 0; j <= i; j++){

      cout<<" PROCESSING PT BIN "<<i<<" "<<j<<endl;

      ptTriggerMin = ptBins[i];
      ptTriggerMax = ptBins[i+1];
      ptAssociatedMin = ptBins[j];
      ptAssociatedMax = ptBins[j+1];


      //Latex
      TString centralityLatex = "Centrality: ";
      centralityLatex += centralityArray[gCentrality-1]; 
      centralityLatex += "%";
      
      TString psiLatex;
      if((psiMin == -0.5)&&(psiMax == 0.5))
	psiLatex = " -7.5^{o} < #varphi - #Psi_{2} < 7.5^{o}"; 
      else if((psiMin == 0.5)&&(psiMax == 1.5))
	psiLatex = " 37.5^{o} < #varphi - #Psi_{2} < 52.5^{o}"; 
      else if((psiMin == 1.5)&&(psiMax == 2.5))
	psiLatex = " 82.5^{o} < #varphi - #Psi_{2} < 97.5^{o}"; 
      else 
	psiLatex = " 0^{o} < #varphi - #Psi_{2} < 180^{o}"; 
      
      TString pttLatex = Form("%.1f",ptTriggerMin);
      pttLatex += " < p_{T,trig} < "; pttLatex += Form("%.1f",ptTriggerMax);
      pttLatex += " GeV/c";
      
      TString ptaLatex = Form("%.1f",ptAssociatedMin);
      ptaLatex += " < p_{T,assoc} < "; ptaLatex += Form("%.1f",ptAssociatedMax);
      ptaLatex += " GeV/c";
      
      TLatex *latexInfo1 = new TLatex();
      latexInfo1->SetNDC();
      latexInfo1->SetTextSize(0.045);
      latexInfo1->SetTextColor(1);
      
      // Open input file
      inFileName = Form("PbPb/%s/Train%d/Fits/correlationFunctionFit",lhcPeriod,gTrainID);
      if(histoName.Contains("PN")) inFileName += "PN";
      else if(histoName.Contains("NP")) inFileName += "NP";
      else if(histoName.Contains("PP")) inFileName += "PP";
      else if(histoName.Contains("NN")) inFileName += "NN";
      inFileName += ".Centrality";  
      inFileName += gCentrality; inFileName += ".Psi";
      if((psiMin == -0.5)&&(psiMax == 0.5)) inFileName += "InPlane.Ptt";
      else if((psiMin == 0.5)&&(psiMax == 1.5)) inFileName += "Intermediate.Ptt";
      else if((psiMin == 1.5)&&(psiMax == 2.5)) inFileName += "OutOfPlane.Ptt";
      else if((psiMin == 2.5)&&(psiMax == 3.5)) inFileName += "Rest.PttFrom";
      else inFileName += "All.PttFrom";
      inFileName += Form("%.1f",ptTriggerMin); inFileName += "To"; 
      inFileName += Form("%.1f",ptTriggerMax); inFileName += "PtaFrom";
      inFileName += Form("%.1f",ptAssociatedMin); inFileName += "To"; 
      inFileName += Form("%.1f",ptAssociatedMax); 
      inFileName += ".root";
      inFile = TFile::Open(inFileName.Data(),"read");
      inFile->ls();
      hTMPData = (TH2D*)inFile->Get(Form("gHist%sCorrelationFunctions",histoName.Data()));
      hTMPRes  = (TH2D*)inFile->Get("gHistResidual");
      hTMPFit  = (TH2D*)inFile->Get(Form("gHist%sCorrelationFunctionsFit",histoName.Data()));

      cQA[i][j]->cd(1);
      hTMPData->DrawCopy("surf1fb");
      latexInfo1->DrawLatex(0.2,0.95,"Data");

      latexInfo1->DrawLatex(0.44,0.88,centralityLatex.Data());
      //latexInfo1->DrawLatex(0.44,0.82,psiLatex.Data());
      latexInfo1->DrawLatex(0.44,0.82,pttLatex.Data());
      latexInfo1->DrawLatex(0.44,0.76,ptaLatex.Data());

      cQA[i][j]->cd(2);
      hTMPFit->DrawCopy("surf1fb");
      latexInfo1->DrawLatex(0.2,0.95,"Fit");

      cQA[i][j]->cd(3);
      hTMPRes->DrawCopy("surf1fb");
      latexInfo1->DrawLatex(0.2,0.95,"Residual");

      cQA[i][j]->cd(4);
      hTMPData->ProjectionX()->DrawCopy("");
      latexInfo1->DrawLatex(0.2,0.95,"Data");

      cQA[i][j]->cd(5);
      hTMPFit->ProjectionX()->DrawCopy("");
      latexInfo1->DrawLatex(0.2,0.95,"Fit");

      cQA[i][j]->cd(6);
      hTMPRes->ProjectionX()->DrawCopy("");
      latexInfo1->DrawLatex(0.2,0.95,"Residual");

      cQA[i][j]->cd(7);
      hTMPData->ProjectionY()->DrawCopy("");
      latexInfo1->DrawLatex(0.2,0.95,"Data");

      cQA[i][j]->cd(8);
      hTMPFit->ProjectionY()->DrawCopy("");
      latexInfo1->DrawLatex(0.2,0.95,"Fit");

      cQA[i][j]->cd(9);
      hTMPRes->ProjectionY()->DrawCopy("");
      latexInfo1->DrawLatex(0.2,0.95,"Residual");


      cQA[i][j]->SaveAs(Form(Form("PbPb/%s/Train%d/figs/correlationFunctionFit_%s_PttFrom%.1fTo%.1fPtaFrom%.1fTo%.1f.eps",lhcPeriod,gTrainID,histoName.Data(),ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax)));
      cQA[i][j]->SaveAs(Form(Form("PbPb/%s/Train%d/figs/correlationFunctionFit_%s_PttFrom%.1fTo%.1fPtaFrom%.1fTo%.1f.pdf",lhcPeriod,gTrainID,histoName.Data(),ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax)));

      // fit parameters
      fFit = (TF2*)inFile->Get("gFitFunction");
      for(Int_t iPar = 0; iPar < 17; iPar++){
	p[iPar][i*kNPtBins+j] = fFit->GetParameter(iPar);
	pE[iPar][i*kNPtBins+j] = fFit->GetParError(iPar);
      }

      inFile->Close();
    }
  }

  TGraphErrors *gPar[17];
  for(Int_t iPar = 0; iPar < 17; iPar++){
    gPar[iPar]  = new TGraphErrors(kNPtBins*kNPtBins,pt,p[iPar],ptE,pE[iPar]);
    gPar[iPar]->SetTitle(pNames[iPar].Data());
    gPar[iPar]->GetXaxis()->SetTitle("p_{T}");
    gPar[iPar]->GetYaxis()->SetTitle(pNames[iPar].Data());
    gPar[iPar]->SetMinimum(0.01);
    gPar[iPar]->SetMaximum(2);
    gPar[iPar]->SetMarkerStyle(20);
    gPar[iPar]->SetMarkerColor(2);
    gPar[iPar]->SetLineColor(2);
  }

  TLatex *latexInfo2 = new TLatex();
  latexInfo2->SetNDC();
  latexInfo2->SetTextSize(0.045);
  latexInfo2->SetTextColor(1);

  TCanvas *cPar = new TCanvas("cPar","parameters",1200,900);
  cPar->Divide(3,3);

  cPar->cd(1);
  gPar[1]->Draw("AP");

  cPar->cd(2);
  gPar[2]->Draw("AP");

  cPar->cd(3);
  gPar[3]->Draw("AP");

  cPar->cd(4);
  gPar[6]->Draw("AP");

  cPar->cd(5);
  gPar[9]->Draw("AP");

  cPar->cd(6);
  gPar[12]->Draw("AP");

  cPar->cd(7);
  gPar[13]->Draw("AP");

  cPar->cd(8);
  gPar[14]->Draw("AP");

  cPar->cd(9);
  gPar[15]->Draw("AP");

  cPar->SaveAs(Form("PbPb/%s/Train%d/figs/correlationFunctionFit_%s_FitParameters.eps",lhcPeriod,gTrainID,histoName.Data()));
  cPar->SaveAs(Form("PbPb/%s/Train%d/figs/correlationFunctionFit_%s_FitParameters.pdf",lhcPeriod,gTrainID,histoName.Data()));

  // delete canvases
  for(Int_t i = 0; i < kNPtBins; i++){
    for(Int_t j = 0; j <= i; j++){
      //delete cQA[i][j];
    }
  }

}
