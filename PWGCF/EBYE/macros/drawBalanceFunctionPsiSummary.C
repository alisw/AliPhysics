const Int_t numberOfCentralityBins = 8;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};

void drawBalanceFunctionPsiSummarySummary(const char* lhcPeriod = "LHC11h",
					      Int_t gTrainID = 250,			      
					      Double_t psiMin = -0.5, 
					      Double_t psiMax = 3.5){
  TFile        *fPar[3];
  TGraphErrors *gPar[3][18];

  Int_t iCentrality[3] = {1,2,3};
  
  for(Int_t iCent = 0 ; iCent < 3; iCent++){

    // open file
    fPar[iCent] = TFile::Open(Form("PbPb/%s/Train%d/figs/balanceFunctionFit_Cent%d_FitParameters.root",lhcPeriod,gTrainID,iCentrality[iCent]));
    if(!fPar[iCent]){
      cerr<<"FILE "<<Form("PbPb/%s/Train%d/figs/correlationFunctionFit_Cent%d_FitParameters.root",lhcPeriod,gTrainID,iCentrality[iCent])<<" not found!"<<endl;
      return;
    } 
    
    // open graph
    for(Int_t iPar = 0 ; iPar < 18; iPar++){
      gPar[iCent][iPar] = (TGraphErrors*)fPar[iCent]->Get(Form("gPar%d",iPar));
      if(!gPar[iCent][iPar]){
	cerr<<"Graph for parameter "<<iPar<<" not found!"<<endl;
	return;
      } 
    }
  }
  

  TCanvas *cSummary[18]; 
  for(Int_t iPar = 0 ; iPar < 18; iPar++){
    
    cSummary[iPar]  = new TCanvas(Form("cSummary%d",iPar),Form("Summary %d",iPar));
    
    // compare centralities
    cSummary[iPar]->cd();
    gPar[0][iPar]->SetMarkerColor(1);
    gPar[0][iPar]->SetLineColor(1);
    gPar[0][iPar]->Draw("AP");
    gPar[1][iPar]->SetMarkerColor(2);
    gPar[1][iPar]->SetLineColor(2);
    gPar[1][iPar]->Draw("P");
    gPar[2][iPar]->SetMarkerColor(4);
    gPar[2][iPar]->SetLineColor(4);
    gPar[2][iPar]->Draw("P");
    
    
    TLegend *legend2 = new TLegend(0.2,0.6,0.85,0.88,"","brNDC");
    setupLegend(legend2,0.065);
    for(Int_t iCent = 0 ; iCent < 3; iCent++){
      legend2->AddEntry(gPar[iCent][iPar],Form("%s \%",centralityArray[iCentrality[iCent]-1].Data()),"lp");
    }
    legend2->Draw();

    cSummary[iPar]->SaveAs(Form("PbPb/%s/Train%d/figs/balanceFunctionFit_FitParameter%d_Summary.eps",lhcPeriod,gTrainID,iPar));
  }
  
  
}
  

void drawBalanceFunctionPsiSummary(const char* lhcPeriod = "LHC11h",
				   Int_t gTrainID = 208,			      
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
  //Double_t ptBins[kNPtBins+1] = {0.2,1.0,2.0,3.0,4.0,8.0};
  const Int_t kNPtBins = 4;
  Double_t ptBins[kNPtBins+1] = {1.0,2.0,3.0,4.0,8.0};
  //const Int_t kNPtBins = 1;
  //Double_t ptBins[kNPtBins+1] = {0.2,2.0};

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
  Double_t p[18][kNPtBins*kNPtBins];
  Double_t pE[18][kNPtBins*kNPtBins];

  for(Int_t iPar = 0; iPar < 18; iPar++){
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
      inFileName = Form("PbPb/%s/Train%d/Fits/balanceFunctionFit2D",lhcPeriod,gTrainID);
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
      inFileName += "_2pMethod.root";
      inFile = TFile::Open(inFileName.Data(),"read");
      inFile->ls();
      hTMPData = (TH2D*)inFile->Get(Form("gHistBalanceFunctionSubtracted"));
      hTMPRes  = (TH2D*)inFile->Get("gHistResidual");
      hTMPFit  = (TH2D*)inFile->Get(Form("gHistBalanceFunctionSubtractedFit"));

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


      cQA[i][j]->SaveAs(Form(Form("PbPb/%s/Train%d/figs/balanceFunctionFit_PttFrom%.1fTo%.1fPtaFrom%.1fTo%.1f.eps",lhcPeriod,gTrainID,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax)));
      cQA[i][j]->SaveAs(Form(Form("PbPb/%s/Train%d/figs/balanceFunctionFit_PttFrom%.1fTo%.1fPtaFrom%.1fTo%.1f.pdf",lhcPeriod,gTrainID,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax)));

      // fit parameters
      fFit = (TF2*)inFile->Get("gFitFunction");
      for(Int_t iPar = 0; iPar < 18; iPar++){
	p[iPar][i*kNPtBins+j] = fFit->GetParameter(iPar);
	pE[iPar][i*kNPtBins+j] = fFit->GetParError(iPar);
	cout<<iPar<<") Parameter "<<fFit->GetParName(iPar)<<" : "<<p[iPar][i*kNPtBins+j]<<" +- "<<pE[iPar][i*kNPtBins+j]<<endl;
      }

      inFile->Close();
    }
  }

  TGraphErrors *gPar[18];
  for(Int_t iPar = 0; iPar < 18; iPar++){
    gPar[iPar]  = new TGraphErrors(kNPtBins*kNPtBins,pt,p[iPar],ptE,pE[iPar]);
    gPar[iPar]->SetName(Form("gPar%d",iPar));
    gPar[iPar]->SetTitle(fFit->GetParName(iPar));
    gPar[iPar]->GetXaxis()->SetTitle("p_{T}");
    gPar[iPar]->GetYaxis()->SetTitle(fFit->GetParName(iPar));
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
  gPar[8]->Draw("AP");

  cPar->cd(5);
  gPar[9]->Draw("AP");

  cPar->cd(6);
  gPar[17]->Draw("AP");

  cPar->cd(7);
  gPar[12]->Draw("AP");

  cPar->cd(8);
  gPar[14]->Draw("AP");

  cPar->cd(9);
  gPar[15]->Draw("AP");

  cPar->SaveAs(Form("PbPb/%s/Train%d/figs/balanceFunctionFit_Cent%d_FitParameters.eps",lhcPeriod,gTrainID,gCentrality));
  cPar->SaveAs(Form("PbPb/%s/Train%d/figs/balanceFunctionFit_Cent%d_FitParameters.pdf",lhcPeriod,gTrainID,gCentrality));

  TFile *fOut = TFile::Open(Form("PbPb/%s/Train%d/figs/balanceFunctionFit_Cent%d_FitParameters.root",lhcPeriod,gTrainID,gCentrality),"RECREATE");
  for(Int_t iPar = 0; iPar < 18; iPar++){
    gPar[iPar]->Write();
  }

  // delete canvases
  for(Int_t i = 0; i < kNPtBins; i++){
    for(Int_t j = 0; j <= i; j++){
      //delete cQA[i][j];
    }
  }

}


