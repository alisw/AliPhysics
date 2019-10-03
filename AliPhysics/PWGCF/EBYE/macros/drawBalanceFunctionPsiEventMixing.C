const Int_t numberOfCentralityBins = 8;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};


void drawBalanceFunctionPsiEventMixing(const char* lhcPeriod = "LHC11h",
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
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libEventMixing");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGCFebye");

  const Int_t kNPtBins = 3;
  Double_t ptBins[kNPtBins+1] = {1.0,2.0,3.0,4.0};

  TString type[3] = {"PN","NN","PP"}

  TCanvas *cEM = new TCanvas("cEM","",1200,900);
  cEM->Divide(3,3);

  TFile *inFile  = NULL; 
  TH2D *hTMPData = NULL;
  TH2D *hTMPEM = NULL;
  TH1D *hTMPData1D = NULL;
  TH1D *hTMPEM1D = NULL;
  TH1D *hTMPRatio1D = NULL;

  // Loop over pt bins
  Double_t ptTriggerMin = 0.0;
  Double_t ptTriggerMax = 0.;
  Double_t ptAssociatedMin = 0.0;
  Double_t ptAssociatedMax = 0.0;
  TString inFileName = "";
  
  for(Int_t i = 0; i < kNPtBins; i++){

      cout<<" PROCESSING PT BIN "<<i<<" "<<endl;

      ptTriggerMin = ptBins[i];
      ptTriggerMax = ptBins[i+1];
      ptAssociatedMin = ptBins[i];
      ptAssociatedMax = ptBins[i+1];


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
      inFileName = Form("PbPb/%s/Train%d/Centrality%d/correlationFunction",lhcPeriod,gTrainID,gCentrality);
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

      for(Int_t j = 0; j < 3; j++){

	hTMPData = (TH2D*)inFile->Get(Form("gHist%sRaw",type[j].Data()));
	hTMPEM   = (TH2D*)inFile->Get(Form("gHist%sMixed",type[j].Data()));
	
	cEM->cd(3*j+i+1);
	hTMPData1D = (TH1D*)hTMPData->ProjectionX(Form("hTMP%d%d",i,j),33,39);
	hTMPEM1D = (TH1D*)hTMPEM->ProjectionX(Form("hTMPEM%d%d",i,j),33,39);
	hTMPData1D->Fit("pol1","0","0",0,1.6);
	hTMPEM1D->Fit("pol1","0","0",0,1.6);
	hTMPRatio1D = (TH1D*)hTMPData1D->Clone(Form("hTMPRatio%d%d",i,j));
	hTMPRatio1D->Divide(hTMPEM1D);
	hTMPData1D->Divide(hTMPData1D->GetFunction("pol1"));
	hTMPEM1D->Divide(hTMPEM1D->GetFunction("pol1"));
	hTMPData1D->SetMinimum(0.8);
	hTMPData1D->SetMaximum(1.2);
	hTMPData1D->GetXaxis()->SetRangeUser(0,1.6);
	hTMPData1D->Draw();
	hTMPEM1D->SetMarkerColor(2);
	hTMPEM1D->SetLineColor(2);
	hTMPEM1D->Draw("same");
	hTMPRatio1D->SetMarkerColor(4);
	hTMPRatio1D->SetLineColor(4);
	hTMPRatio1D->Draw("same");


	TLegend *legend1 = new TLegend(0.24,0.17,0.5,0.4,"","brNDC");
	setupLegend(legend1,0.065);
	legend1->AddEntry(hTMPRatio1D,"Correlation Function","lp");
	legend1->AddEntry(hTMPEM1D,"Event Mixing/Pol1(EM)","lp");
	legend1->AddEntry(hTMPData1D,"Raw Data/Pol1(Raw)","lp");	  
	legend1->Draw();


	latexInfo1->DrawLatex(0.24,0.82,centralityLatex.Data());
	latexInfo1->DrawLatex(0.24,0.76,pttLatex.Data());
	latexInfo1->DrawLatex(0.24,0.70,ptaLatex.Data());
      }
      
  }

  cEM->SaveAs(Form(Form("PbPb/%s/Train%d/figs/eventMixingDivided_Cent%d.eps",lhcPeriod,gTrainID,gCentrality)));

}

//____________________________________________________________//
void setupLegend(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}
