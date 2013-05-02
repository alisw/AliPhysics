static const Int_t nCentralityBins = 9;
TString strCentrality[nCentralityBins] = {"0-5","5-10","10-20",
                                          "20-30","30-40","40-50",
                                          "50-60","60-70","70-80"};
TString strCentralityBinsLabel[nCentralityBins] = {"Centrality: 0-5%",
						   "Centrality: 5-10%",
						   "Centrality: 10-20%",
						   "Centrality: 20-30%",
						   "Centrality: 30-40%",
						   "Centrality: 40-50%",
						   "Centrality: 50-60%",
						   "Centrality: 60-70%",
						   "Centrality: 70-80%"};

//void drawCorrection(const char* filenameEffCont = "AnalysisResults_HIJING.root"){  

void drawCorrection(const char* filenameEffCont = "mergedAnalysisResults_proofPbPb_128.root"){

  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  //gSystem->Load("libANALYSIS.so");
  //gSystem->Load("libANALYSISalice.so");
  Int_t markerStyle = 20;
  Int_t markerColor = 1;
  Int_t fillColor   = 9;

  //_______________________________________________________________//
  //Open the input file
  TFile *f = TFile::Open(filenameEffCont);
  if(!f->IsOpen()) {
    Printf("File not found!!!");
    break;
  }
  
  //_______________________________________________________________//
  //Get the TDirectoryFile
  TDirectoryFile *dirEffCont = dynamic_cast<TDirectoryFile *>(f->Get("PWGCFEbyE.outputBalanceFunctionEffContAnalysis"));
  if(!dirEffCont) {
    Printf("TDirectoryFile not found!!!");
    break;
  }

  //_______________________________________________________________//
  TList *listEffCont[nCentralityBins];
  TString listName;
  TString listName1;
  TString histName;
  TString histName1;

  TH3D* h1d[nCentralityBins];
  TH3D* h1n[nCentralityBins];
  TH3D* h2d[nCentralityBins];
  TH3D* h2n[nCentralityBins];
  TH3F* h3d[nCentralityBins];
  TH3F* h4d[nCentralityBins];
  TH3F* fHistCorrectionPlus[nCentralityBins];
  TH3F* fHistCorrectionMinus[nCentralityBins];

  TCanvas *correctionMatrix[nCentralityBins]; 
  TString canvasName;

  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    Printf("================Centrality: %s================",  strCentrality[iBin].Data());
    
    listName = "listEffContBF_V0M_"; listName += strCentrality[iBin].Data(); listName +="_10";
    listName1 = "listQA_V0M_"; listName1 += strCentrality[iBin].Data(); listName1 +="_10";
    listEffCont[iBin] = dynamic_cast<TList *>(dirEffCont->Get(listName.Data()));
    listEffCont[iBin]->ls();

    //____________Efficiency Plus
    //correction->cd(1);
    h1d[iBin] = dynamic_cast<TH3D*>(listEffCont[iBin]->FindObject("fHistGeneratedEtaPtPhiPlus"));
    h1d[iBin]->Sumw2();
    h1n[iBin] = dynamic_cast<TH3D*>(listEffCont[iBin]->FindObject("fHistSurvivedEtaPtPhiPlus"));
    h1n[iBin]->Sumw2();
    h1d[iBin]->Divide(h1n[iBin]);
    //h1d->GetYaxis()->SetTitleOffset(1.5);
    //h1d->SetTitle("Efficiency (+)");
    //h1d->SetName("fHistEfficiencyPlus");
    //h1d->Draw("P");
    
    //____________Efficiency Minus
    //correction->cd(2);
    h2d[iBin] = dynamic_cast<TH3D*>(listEffCont[iBin]->FindObject("fHistGeneratedEtaPtPhiMinus"));
    h2d[iBin]->Sumw2();
    h2n[iBin] = dynamic_cast<TH3D*>(listEffCont[iBin]->FindObject("fHistSurvivedEtaPtPhiMinus"));
    h2n[iBin]->Sumw2();
    h2d[iBin]->Divide(h2n[iBin]);
    //h2d->GetYaxis()->SetTitleOffset(1.5);
    //h2d->SetTitle("Efficiency (-)");
    //h2d->SetName("fHistEfficiencyMinus");
    //h2d->Draw("P");
    
    //__________Contamination
    //correction->cd(3);
    h4d[iBin] = dynamic_cast<TH3F*>(listEffCont[iBin]->FindObject("fHistContaminationPrimaries"));
    h4d[iBin]->SetName("fHistContaminationPrimaries");
    h4d[iBin]->Sumw2();
  
    h3d[iBin] = dynamic_cast<TH3F*>(listEffCont[iBin]->FindObject("fHistContaminationSecondaries"));
    h3d[iBin]->SetName("fHistContaminationSecondaries");
    h3d[iBin]->Sumw2();

    h4d[iBin]->Add(h4d[iBin]);
    h3d[iBin]->Divide(h4d[iBin]);
    //h3d->GetYaxis()->SetTitleOffset(2.0);
    //h3d->GetXaxis()->SetTitleOffset(1.5);
    //h3d->SetTitle("Contamination Secondaries");
    //h3d->Draw("P");
 
    //_____________________________________________________________
    //CORRECTION

    Int_t binsX = h3d[0]->GetNbinsX();
    Int_t binsY = h3d[0]->GetNbinsY();
    Int_t binsZ = h3d[0]->GetNbinsZ();

    for (Int_t iHistBinsX = 1; iHistBinsX <binsX+1 ; iHistBinsX++) {
      for (Int_t iHistBinsY = 1; iHistBinsY <binsY+1 ; iHistBinsY++) {
	for (Int_t iHistBinsZ = 1; iHistBinsZ <binsZ+1 ; iHistBinsZ++) {
	  h3d[iBin]->SetBinContent(iHistBinsX,iHistBinsY,iHistBinsZ,1 - h3d[iBin]->GetBinContent(iHistBinsX,iHistBinsY,iHistBinsZ));
	}
      }
    }
        
    correctionMatrix[iBin] = new TCanvas(canvasName.Data(),canvasName.Data(),0,0+iBin*50,1400,1000);
    correctionMatrix[iBin]->Divide(2,1);
    
    //Correction Maps Plus 
    histName = "fHistCorrectionPlus"; 
    histName += strCentrality[iBin].Data();
    correctionMatrix[iBin]->cd(1);
    fHistCorrectionPlus[iBin] = dynamic_cast<TH3F *>(h3d[iBin]->Clone());
    fHistCorrectionPlus[iBin]->Divide(h1d[iBin]);  
    fHistCorrectionPlus[iBin]->SetName(histName.Data());
    fHistCorrectionPlus[iBin]->GetYaxis()->SetTitleOffset(2.0);
    fHistCorrectionPlus[iBin]->GetXaxis()->SetTitleOffset(1.5);
    fHistCorrectionPlus[iBin]->GetZaxis()->SetTitleOffset(2.0);
    fHistCorrectionPlus[iBin]->SetTitle("Correction Plus");
    fHistCorrectionPlus[iBin]->DrawCopy("");
    
    //Correction Maps Minus 
    histName1 = "fHistCorrectionMinus"; 
    histName1 += strCentrality[iBin].Data();
    correctionMatrix[iBin]->cd(2);
    fHistCorrectionMinus[iBin] = dynamic_cast<TH3F *>(h3d[iBin]->Clone());
    fHistCorrectionMinus[iBin]->Divide(h2d[iBin]);
    fHistCorrectionMinus[iBin]->SetName(histName1.Data());
    fHistCorrectionMinus[iBin]->GetYaxis()->SetTitleOffset(2.0);
    fHistCorrectionMinus[iBin]->GetXaxis()->SetTitleOffset(1.5);
    fHistCorrectionMinus[iBin]->GetZaxis()->SetTitleOffset(2.0);
    fHistCorrectionMinus[iBin]->SetTitle("Correction Minus");
    fHistCorrectionMinus[iBin]->DrawCopy("");
  }
  
  //____________________________________________________________
  //Output files
  TFile *fCorrectionMaps = TFile::Open("CorrectionMaps.root",
					  "recreate");
  
  for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
    fHistCorrectionPlus[iBin]->Write();
    fHistCorrectionMinus[iBin]->Write();
  }
  fCorrectionMaps->Close();  
}
