const Int_t numberOfCentralityBins = 9;
TString centralityArray[numberOfCentralityBins] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
Double_t gMinCentrality[numberOfCentralityBins] = {0.,5.,10.,20.,30.,40.,50.,60.,70.};
Double_t gMaxCentrality[numberOfCentralityBins] = {5.,10.,20.,30.,40.,50.,60.,70.,80.};
TString gAnalysisType[7] = {"y","eta","qlong","qout","qside","qinv","phi"};

const Int_t gRebin = 1;
void drawCorrelationFunctionPsi(const char* filename = "AnalysisResults.root", 
				Double_t psiMin = -0.5, 
				Double_t psiMax = 0.5,
				Double_t ptTriggerMin = -1.,
				Double_t ptTriggerMax = -1.,
				Double_t ptAssociatedMin = -1.,
				Double_t ptAssociatedMax = -1.) {
  //Macro that draws the correlation functions from the balance function
  //analysis vs the reaction plane
  //Author: Panos.Christakoglou@nikhef.nl
  //Load the PWG2ebye library
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGCFebye.so");

  //Prepare the objects and return them
  TList *list = GetListOfObjects(filename);
  if(!list) {
    Printf("The TList object was not created");
    return;
  }
  else 
    draw(list,psiMin,psiMax,
	 ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
}

//______________________________________________________//
TList *GetListOfObjects(const char* filename) {
  //Get the TList objects (QA, bf, bf shuffled)
  TList *listQA = 0x0;
  TList *listBF = 0x0;
  TList *listBFShuffling = 0x0;
  
  //Open the file
  TFile *f = TFile::Open(filename);
  if((!f)||(!f->IsOpen())) {
    Printf("The file %s is not found. Aborting...",filename);
    return listBF;
  }
  //f->ls();
  
  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("PWGCFEbyE.outputBalanceFunctionPsiAnalysis"));
  if(!dir) {   
    Printf("The TDirectoryFile is not found. Aborting...",filename);
    return listBF;
  }
  //dir->ls();
  
  TString listBFName = "listBFPsi";
  listBF = dynamic_cast<TList *>(dir->Get(listBFName.Data()));
  //listBF->ls();

  //Get the histograms
  TString histoName = "fHistPV0M";
  AliTHn *fHistP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));  
  if(!fHistP) {
    Printf("fHistP %s not found!!!",histoName.Data());
    break;
  }
  fHistP->FillParent(); fHistP->DeleteContainers();

  histoName = "fHistNV0M"; 
  AliTHn *fHistN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistN) {
    Printf("fHistN %s not found!!!",histoName.Data());
    break;
  }
  fHistN->FillParent(); fHistN->DeleteContainers();
    
  histoName = "fHistPNV0M"; 
  AliTHn *fHistPN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistPN) {
    Printf("fHistPN %s not found!!!",histoName.Data());
    break;
  }
  fHistPN->FillParent(); fHistPN->DeleteContainers();
  
  histoName = "fHistNPV0M";
  AliTHn *fHistNP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistNP) {
    Printf("fHistNP %s not found!!!",histoName.Data());
    break;
  }
  fHistNP->FillParent(); fHistNP->DeleteContainers();

  histoName = "fHistPPV0M";
  AliTHn *fHistPP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistPP) {
    Printf("fHistPP %s not found!!!",histoName.Data());
    break;
  }
  fHistPP->FillParent(); fHistPP->DeleteContainers();

  histoName = "fHistNNV0M";
  AliTHn *fHistNN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistNN) {
    Printf("fHistNN %s not found!!!",histoName.Data());
    break;
  }
  fHistNN->FillParent(); fHistNN->DeleteContainers();
  
  return listBF;
}

//______________________________________________________//
void draw(TList *list, Double_t psiMin, Double_t psiMax,
	  Double_t ptTriggerMin, Double_t ptTriggerMax,
	  Double_t ptAssociatedMin, Double_t ptAssociatedMax) {
  //Draws the correlation functions for every centrality bin
  //(+-), (-+), (++), (--)
  gROOT->LoadMacro("~/SetPlotStyle.C");
  SetPlotStyle();
  gStyle->SetPalette(1,0);
  
  AliTHn *hP = NULL;
  AliTHn *hN = NULL;
  AliTHn *hPN = NULL;
  AliTHn *hNP = NULL;
  AliTHn *hPP = NULL;
  AliTHn *hNN = NULL;
  
  hP = (AliTHn*) list->FindObject("fHistPV0M");
  hN = (AliTHn*) list->FindObject("fHistNV0M");
  hPN = (AliTHn*) list->FindObject("fHistPNV0M");
  hNP = (AliTHn*) list->FindObject("fHistNPV0M");
  hPP = (AliTHn*) list->FindObject("fHistPPV0M");
  hNN = (AliTHn*) list->FindObject("fHistNNV0M");

  //Create the AliBalancePsi object and fill it with the AliTHn objects
  AliBalancePsi *b = new AliBalancePsi();
  b->SetHistNp(hP);
  b->SetHistNn(hN);
  b->SetHistNpn(hPN);
  b->SetHistNnp(hNP);
  b->SetHistNpp(hPP);
  b->SetHistNnn(hNN);
  TH2D *gHistPN;
  TH2D *gHistNP;
  TH2D *gHistPP;
  TH2D *gHistNN;
  
  TCanvas *cPN;
  TCanvas *cNP;
  TCanvas *cPP;
  TCanvas *cNN;
  TString histoTitle, pngName;
  //loop over the centrality bins
  //for(Int_t iCentralityBin = 0; iCentralityBin < numberOfCentralityBins; iCentralityBin++) {
  
  histoTitle = "(+-) | Centrality: "; 
  histoTitle += centralityArray[6]; 
  histoTitle += "%";
  histoTitle += " | "; histoTitle += psiMin; 
  histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
  gHistPN = b->GetCorrelationFunctionPN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistPN->GetYaxis()->SetTitleOffset(1.5);
  gHistPN->SetTitle(histoTitle.Data());
  cPN = new TCanvas(histoTitle.Data(),"",0,0,400,400);
  cPN->SetFillColor(10); 
  cPN->SetHighLightColor(10);
  gHistPN->DrawCopy("lego2");
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[6]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".PositiveNegative.png";
  cPN->SaveAs(pngName.Data());
  
  histoTitle = "(-+) | Centrality: "; 
  histoTitle += centralityArray[6]; 
  histoTitle += "%";
  histoTitle += " | "; histoTitle += psiMin; 
  histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
  gHistNP = b->GetCorrelationFunctionNP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistNP->GetYaxis()->SetTitleOffset(1.5);
  gHistNP->SetTitle(histoTitle.Data());
  cNP = new TCanvas(histoTitle.Data(),"",400,0,400,400);
  cNP->SetFillColor(10); 
  cNP->SetHighLightColor(10);
  gHistNP->DrawCopy("lego2");
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[6]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".NegativePositive.png";
  cNP->SaveAs(pngName.Data());
  
  histoTitle = "(++) | Centrality: "; 
  histoTitle += centralityArray[6]; 
  histoTitle += "%";
  histoTitle += " | "; histoTitle += psiMin; 
  histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
  gHistPP = b->GetCorrelationFunctionPP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistPP->GetYaxis()->SetTitleOffset(1.5);
  gHistPP->SetTitle(histoTitle.Data());
  cPP = new TCanvas(histoTitle.Data(),"",0,400,400,400);
  cPP->SetFillColor(10); 
  cPP->SetHighLightColor(10);
  gHistPP->DrawCopy("lego2");
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[6]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".PositivePositive.png";
  cPP->SaveAs(pngName.Data());
  
  histoTitle = "(--) | Centrality: "; 
  histoTitle += centralityArray[6]; 
  histoTitle += "%";
  histoTitle += " | "; histoTitle += psiMin; 
  histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
  gHistNN = b->GetCorrelationFunctionNN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistNN->GetYaxis()->SetTitleOffset(1.5);
  gHistNN->SetTitle(histoTitle.Data());
  cNN = new TCanvas(histoTitle.Data(),"",400,400,400,400);
  cNN->SetFillColor(10); 
  cNN->SetHighLightColor(10);
  gHistNN->DrawCopy("lego2");
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[6]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".NegativeNegative.png";
  cNN->SaveAs(pngName.Data());
  //}//end of loop over centralities
}

