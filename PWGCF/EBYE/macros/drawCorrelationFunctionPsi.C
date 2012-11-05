const Int_t numberOfCentralityBins = 8;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};

const Int_t gRebin = 1;
void drawCorrelationFunctionPsi(const char* filename = "AnalysisResults.root", 
				Int_t gCentrality = 1,
				Int_t gBit = -1,
				const char* gCentralityEstimator = 0x0,
				Bool_t kShowShuffled = kFALSE, 
				Bool_t kShowMixed = kTRUE, 
				Double_t psiMin = -0.5, 
				Double_t psiMax = 3.5,
				Double_t ptTriggerMin = -1.,
				Double_t ptTriggerMax = -1.,
				Double_t ptAssociatedMin = -1.,
				Double_t ptAssociatedMax = -1.) {
  //Macro that draws the correlation functions from the balance function
  //analysis vs the reaction plane
  //Author: Panos.Christakoglou@nikhef.nl
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

  //Prepare the objects and return them
  TList *list = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,0);
  TList *listShuffled = NULL;
  if(kShowShuffled) listShuffled = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,1);
  TList *listMixed = NULL;
  if(kShowMixed) listMixed = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,2);

  if(!list) {
    Printf("The TList object was not created");
    return;
  }
  else 
    draw(list,listShuffled,listMixed,gCentrality,psiMin,psiMax,
	 ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
}

//______________________________________________________//
TList *GetListOfObjects(const char* filename,
			Int_t gCentrality,
			Int_t gBit,
			const char *gCentralityEstimator,
			Int_t kData = 1) {
  //Get the TList objects (QA, bf, bf shuffled)
  TList *listBF = 0x0;
  
  //Open the file
  TFile *f = TFile::Open(filename,"UPDATE");
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
  
  TString listBFName;
  if(kData == 0) {
    //cout<<"no shuffling - no mixing"<<endl;
    listBFName = "listBFPsi_";
  }
  else if(kData == 1) {
    //cout<<"shuffling - no mixing"<<endl;
    listBFName = "listBFPsiShuffled_";
  }
  else if(kData == 2) {
    //cout<<"no shuffling - mixing"<<endl;
    listBFName = "listBFPsiMixed_";
  }
  listBFName += centralityArray[gCentrality-1];
  if(gBit > -1) {
    listBFName += "_Bit"; listBFName += gBit; }
  if(gCentralityEstimator) {
    listBFName += "_"; listBFName += gCentralityEstimator;}

  // histograms were already retrieved (in first iteration)
  if(dir->Get(Form("%s_histograms",listBFName.Data()))){
    listBF = dynamic_cast<TList *>(dir->Get(Form("%s_histograms",listBFName.Data())));
  }

  // histograms were not yet retrieved (this is the first iteration)
  else{

    listBF = dynamic_cast<TList *>(dir->Get(listBFName.Data()));
    cout<<"======================================================="<<endl;
    cout<<"List name: "<<listBF->GetName()<<endl;
    //listBF->ls();
    
    //Get the histograms
    TString histoName;
    if(kData == 0)
      histoName = "fHistPV0M";
    else if(kData == 1)
      histoName = "fHistP_shuffleV0M";
    else if(kData == 2)
      histoName = "fHistPV0M";
    AliTHn *fHistP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));  
    if(!fHistP) {
      Printf("fHistP %s not found!!!",histoName.Data());
      break;
    }
    fHistP->FillParent(); fHistP->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistNV0M";
    if(kData == 1)
      histoName = "fHistN_shuffleV0M";
    if(kData == 2)
      histoName = "fHistNV0M";
    AliTHn *fHistN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistN) {
      Printf("fHistN %s not found!!!",histoName.Data());
      break;
    }
    fHistN->FillParent(); fHistN->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistPNV0M";
    if(kData == 1)
      histoName = "fHistPN_shuffleV0M";
    if(kData == 2)
      histoName = "fHistPNV0M";
    AliTHn *fHistPN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistPN) {
      Printf("fHistPN %s not found!!!",histoName.Data());
      break;
    }
    fHistPN->FillParent(); fHistPN->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistNPV0M";
    if(kData == 1)
      histoName = "fHistNP_shuffleV0M";
    if(kData == 2)
      histoName = "fHistNPV0M";
    AliTHn *fHistNP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistNP) {
      Printf("fHistNP %s not found!!!",histoName.Data());
      break;
    }
    fHistNP->FillParent(); fHistNP->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistPPV0M";
    if(kData == 1)
      histoName = "fHistPP_shuffleV0M";
    if(kData == 2)
      histoName = "fHistPPV0M";
    AliTHn *fHistPP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistPP) {
      Printf("fHistPP %s not found!!!",histoName.Data());
      break;
    }
    fHistPP->FillParent(); fHistPP->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistNNV0M";
    if(kData == 1)
      histoName = "fHistNN_shuffleV0M";
    if(kData == 2)
      histoName = "fHistNNV0M";
    AliTHn *fHistNN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistNN) {
      Printf("fHistNN %s not found!!!",histoName.Data());
      break;
    }
    fHistNN->FillParent(); fHistNN->DeleteContainers();
    
    dir->cd();
    listBF->Write(Form("%s_histograms",listBFName.Data()), TObject::kSingleKey);
    
  }// first iteration
  
  f->Close();
  
  return listBF;
}

//______________________________________________________//
void draw(TList *list, TList *listBFShuffled, TList *listBFMixed, 
	  Int_t gCentrality, Double_t psiMin, Double_t psiMax,
	  Double_t ptTriggerMin, Double_t ptTriggerMax,
	  Double_t ptAssociatedMin, Double_t ptAssociatedMax) {
  //Draws the correlation functions for every centrality bin
  //(+-), (-+), (++), (--)  
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

  //balance function shuffling
  AliTHn *hPShuffled = NULL;
  AliTHn *hNShuffled = NULL;
  AliTHn *hPNShuffled = NULL;
  AliTHn *hNPShuffled = NULL;
  AliTHn *hPPShuffled = NULL;
  AliTHn *hNNShuffled = NULL;
  if(listBFShuffled) {
    //listBFShuffled->ls();
    
    hPShuffled = (AliTHn*) listBFShuffled->FindObject("fHistP_shuffleV0M");
    hPShuffled->SetName("gHistPShuffled");
    hNShuffled = (AliTHn*) listBFShuffled->FindObject("fHistN_shuffleV0M");
    hNShuffled->SetName("gHistNShuffled");
    hPNShuffled = (AliTHn*) listBFShuffled->FindObject("fHistPN_shuffleV0M");
    hPNShuffled->SetName("gHistPNShuffled");
    hNPShuffled = (AliTHn*) listBFShuffled->FindObject("fHistNP_shuffleV0M");
    hNPShuffled->SetName("gHistNPShuffled");
    hPPShuffled = (AliTHn*) listBFShuffled->FindObject("fHistPP_shuffleV0M");
    hPPShuffled->SetName("gHistPPShuffled");
    hNNShuffled = (AliTHn*) listBFShuffled->FindObject("fHistNN_shuffleV0M");
    hNNShuffled->SetName("gHistNNShuffled");
    
    AliBalancePsi *bShuffled = new AliBalancePsi();
    bShuffled->SetHistNp(hPShuffled);
    bShuffled->SetHistNn(hNShuffled);
    bShuffled->SetHistNpn(hPNShuffled);
    bShuffled->SetHistNnp(hNPShuffled);
    bShuffled->SetHistNpp(hPPShuffled);
    bShuffled->SetHistNnn(hNNShuffled);
  }

  //balance function mixing
  AliTHn *hPMixed = NULL;
  AliTHn *hNMixed = NULL;
  AliTHn *hPNMixed = NULL;
  AliTHn *hNPMixed = NULL;
  AliTHn *hPPMixed = NULL;
  AliTHn *hNNMixed = NULL;

  if(listBFMixed) {
    //listBFMixed->ls();

    hPMixed = (AliTHn*) listBFMixed->FindObject("fHistPV0M");
    hPMixed->SetName("gHistPMixed");
    hNMixed = (AliTHn*) listBFMixed->FindObject("fHistNV0M");
    hNMixed->SetName("gHistNMixed");
    hPNMixed = (AliTHn*) listBFMixed->FindObject("fHistPNV0M");
    hPNMixed->SetName("gHistPNMixed");
    hNPMixed = (AliTHn*) listBFMixed->FindObject("fHistNPV0M");
    hNPMixed->SetName("gHistNPMixed");
    hPPMixed = (AliTHn*) listBFMixed->FindObject("fHistPPV0M");
    hPPMixed->SetName("gHistPPMixed");
    hNNMixed = (AliTHn*) listBFMixed->FindObject("fHistNNV0M");
    hNNMixed->SetName("gHistNNMixed");
    
    AliBalancePsi *bMixed = new AliBalancePsi();
    bMixed->SetHistNp(hPMixed);
    bMixed->SetHistNn(hNMixed);
    bMixed->SetHistNpn(hPNMixed);
    bMixed->SetHistNnp(hNPMixed);
    bMixed->SetHistNpp(hPPMixed);
    bMixed->SetHistNnn(hNNMixed);
  }

  TH2D *gHistPN[4];
  TH2D *gHistNP[4];
  TH2D *gHistPP[4];
  TH2D *gHistNN[4];
  
  TCanvas *cPN[4];
  TCanvas *cNP[4];
  TCanvas *cPP[4];
  TCanvas *cNN[4];
  TString histoTitle, pngName;
  
  //(+-)
  histoTitle = "(+-) | Centrality: ";
  histoTitle += centralityArray[gCentrality-1]; 
  histoTitle += "%";
  if((psiMin == -0.5)&&(psiMax == 0.5))
    histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
  else if((psiMin == 0.5)&&(psiMax == 1.5))
    histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
  else if((psiMin == 1.5)&&(psiMax == 2.5))
    histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
  else 
    histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 

  gHistPN[0] = b->GetCorrelationFunctionPN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistPN[0]->GetYaxis()->SetTitleOffset(1.5);
  gHistPN[0]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
  gHistPN[0]->SetTitle(histoTitle.Data());
  cPN[0] = new TCanvas("cPN0","",0,0,600,500);
  cPN[0]->SetFillColor(10); cPN[0]->SetHighLightColor(10);
  gHistPN[0]->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  //gPad->SetPhi(130); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".PositiveNegative.png";
  cPN[0]->SaveAs(pngName.Data());
  
  if(listBFShuffled) {
    histoTitle = "(+-) shuffled | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistPN[1] = bShuffled->GetCorrelationFunctionPN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistPN[1]->GetYaxis()->SetTitleOffset(1.5);
    gHistPN[1]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistPN[1]->SetTitle(histoTitle.Data());
    cPN[1] = new TCanvas("cPN1","",0,100,600,500);
    cPN[1]->SetFillColor(10); 
    cPN[1]->SetHighLightColor(10);
    gHistPN[1]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "DeltaPhiDeltaEtaShuffled.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositiveNegative.png";
    cPN[1]->SaveAs(pngName.Data());
  }

  if(listBFMixed) {
    histoTitle = "(+-) mixed | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistPN[2] = bMixed->GetCorrelationFunctionPN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistPN[2]->GetYaxis()->SetTitleOffset(1.5);
    gHistPN[2]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistPN[2]->SetTitle(histoTitle.Data());
    cPN[2] = new TCanvas("cPN2","",0,200,600,500);
    cPN[2]->SetFillColor(10); 
    cPN[2]->SetHighLightColor(10);
    gHistPN[2]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "DeltaPhiDeltaEtaMixed.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositiveNegative.png";
    cPN[2]->SaveAs(pngName.Data());

    //Correlation function (+-)
    gHistPN[3] = dynamic_cast<TH2D *>(gHistPN[0]->Clone());
    gHistPN[3]->Divide(gHistPN[2]);
    gHistPN[3]->GetXaxis()->SetRangeUser(-1.5,1.5);
    gHistPN[3]->GetZaxis()->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
    cPN[3] = new TCanvas("cPN3","",0,300,600,500);
    cPN[3]->SetFillColor(10); 
    cPN[3]->SetHighLightColor(10);
    gHistPN[3]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "CorrelationFunction.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositiveNegative.png";
    cPN[3]->SaveAs(pngName.Data());
  }

  //(-+)
  histoTitle = "(-+) | Centrality: "; 
  histoTitle += centralityArray[gCentrality-1]; 
  histoTitle += "%";
  if((psiMin == -0.5)&&(psiMax == 0.5))
    histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
  else if((psiMin == 0.5)&&(psiMax == 1.5))
    histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
  else if((psiMin == 1.5)&&(psiMax == 2.5))
    histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
  else 
    histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 

  gHistNP[0] = b->GetCorrelationFunctionNP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistNP[0]->GetYaxis()->SetTitleOffset(1.5);
  gHistNP[0]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
  gHistNP[0]->SetTitle(histoTitle.Data());
  cNP[0] = new TCanvas("cNP0","",100,0,600,500);
  cNP[0]->SetFillColor(10); 
  cNP[0]->SetHighLightColor(10);
  gHistNP[0]->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  //gPad->SetPhi(130); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".NegativePositive.png";
  cNP[0]->SaveAs(pngName.Data());

  if(listBFShuffled) {
    histoTitle = "(-+) shuffled | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistNP[1] = bShuffled->GetCorrelationFunctionNP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistNP[1]->GetYaxis()->SetTitleOffset(1.5);
    gHistNP[1]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistNP[1]->SetTitle(histoTitle.Data());
    cNP[1] = new TCanvas("cNP1","",100,100,600,500);
    cNP[1]->SetFillColor(10); 
    cNP[1]->SetHighLightColor(10);
    gHistNP[1]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();
    pngName = "DeltaPhiDeltaEtaShuffled.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativePositive.png";
    cNP[1]->SaveAs(pngName.Data());
  }

  if(listBFMixed) {
    histoTitle = "(-+) mixed | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistNP[2] = bMixed->GetCorrelationFunctionNP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistNP[2]->GetYaxis()->SetTitleOffset(1.5);
    gHistNP[2]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistNP[2]->SetTitle(histoTitle.Data());
    cNP[2] = new TCanvas("cNP2","",100,200,600,500);
    cNP[2]->SetFillColor(10); 
    cNP[2]->SetHighLightColor(10);
    gHistNP[2]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();
    pngName = "DeltaPhiDeltaEtaMixed.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativePositive.png";
    cNP[2]->SaveAs(pngName.Data());

    //Correlation function (-+)
    gHistNP[3] = dynamic_cast<TH2D *>(gHistNP[0]->Clone());
    gHistNP[3]->Divide(gHistNP[2]);
    gHistNP[3]->GetXaxis()->SetRangeUser(-1.5,1.5);
    gHistNP[3]->GetZaxis()->SetTitle("C_{-+}(#Delta#eta,#Delta#varphi)");
    cNP[3] = new TCanvas("cNP3","",100,300,600,500);
    cNP[3]->SetFillColor(10); 
    cNP[3]->SetHighLightColor(10);
    gHistNP[3]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "CorrelationFunction.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativePositive.png";
    cNP[3]->SaveAs(pngName.Data());
  }
  
  //(++)
  histoTitle = "(++) | Centrality: "; 
  histoTitle += centralityArray[gCentrality-1]; 
  histoTitle += "%";
  if((psiMin == -0.5)&&(psiMax == 0.5))
    histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
  else if((psiMin == 0.5)&&(psiMax == 1.5))
    histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
  else if((psiMin == 1.5)&&(psiMax == 2.5))
    histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
  else 
    histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 

  gHistPP[0] = b->GetCorrelationFunctionPP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistPP[0]->GetYaxis()->SetTitleOffset(1.5);
  gHistPP[0]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
  gHistPP[0]->SetTitle(histoTitle.Data());
  cPP[0] = new TCanvas("cPP0","",200,0,600,500);
  cPP[0]->SetFillColor(10); 
  cPP[0]->SetHighLightColor(10);
  gHistPP[0]->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  //gPad->SetPhi(130); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".PositivePositive.png";
  cPP[0]->SaveAs(pngName.Data());
  
  if(listBFShuffled) {
    histoTitle = "(++) shuffled | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistPP[1] = bShuffled->GetCorrelationFunctionPP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistPP[1]->GetYaxis()->SetTitleOffset(1.5);
    gHistPP[1]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistPP[1]->SetTitle(histoTitle.Data());
    cPP[1] = new TCanvas("cPP1","",200,100,600,500);
    cPP[1]->SetFillColor(10); 
    cPP[1]->SetHighLightColor(10);
    gHistPP[1]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();
    pngName = "DeltaPhiDeltaEtaShuffled.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositivePositive.png";
    cPP[1]->SaveAs(pngName.Data());
  }

  if(listBFMixed) {
    histoTitle = "(++) mixed | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistPP[2] = bMixed->GetCorrelationFunctionPP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistPP[2]->GetYaxis()->SetTitleOffset(1.5);
    gHistPP[2]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistPP[2]->SetTitle(histoTitle.Data());
    cPP[2] = new TCanvas("cPP2","",200,200,600,500);
    cPP[2]->SetFillColor(10); 
    cPP[2]->SetHighLightColor(10);
    gHistPP[2]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();
    pngName = "DeltaPhiDeltaEtaMixed.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositivePositive.png";
    cPP[2]->SaveAs(pngName.Data());

    //Correlation function (++)
    gHistPP[3] = dynamic_cast<TH2D *>(gHistPP[0]->Clone());
    gHistPP[3]->Divide(gHistPP[2]);
    gHistPP[3]->GetXaxis()->SetRangeUser(-1.5,1.5);
    gHistPP[3]->GetZaxis()->SetTitle("C_{++}(#Delta#eta,#Delta#varphi)");
    cPP[3] = new TCanvas("cPP3","",200,300,600,500);
    cPP[3]->SetFillColor(10); 
    cPP[3]->SetHighLightColor(10);
    gHistPP[3]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "CorrelationFunction.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositivePositive.png";
    cPP[3]->SaveAs(pngName.Data());
  }

  //(--)
  histoTitle = "(--) | Centrality: "; 
  histoTitle += centralityArray[gCentrality-1]; 
  histoTitle += "%";
  if((psiMin == -0.5)&&(psiMax == 0.5))
    histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
  else if((psiMin == 0.5)&&(psiMax == 1.5))
    histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
  else if((psiMin == 1.5)&&(psiMax == 2.5))
    histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
  else 
    histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 

  gHistNN[0] = b->GetCorrelationFunctionNN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistNN[0]->GetYaxis()->SetTitleOffset(1.5);
  gHistNN[0]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
  gHistNN[0]->SetTitle(histoTitle.Data());
  cNN[0] = new TCanvas("cNN0","",300,0,600,500);
  cNN[0]->SetFillColor(10); 
  cNN[0]->SetHighLightColor(10);
  gHistNN[0]->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  //gPad->SetPhi(-60); // default is 30
  gPad->Update();
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".NegativeNegative.png";
  cNN[0]->SaveAs(pngName.Data());

  if(listBFShuffled) {
    histoTitle = "(--) shuffled | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistNN[1] = bShuffled->GetCorrelationFunctionNN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistNN[1]->GetYaxis()->SetTitleOffset(1.5);
    gHistNN[1]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistNN[1]->SetTitle(histoTitle.Data());
    cNN[1] = new TCanvas("cNN1","",300,100,600,500);
    cNN[1]->SetFillColor(10); 
    cNN[1]->SetHighLightColor(10);
    gHistNN[1]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();
    pngName = "DeltaPhiDeltaEtaShuffled.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativeNegative.png";
    cNN[1]->SaveAs(pngName.Data());
  }

  if(listBFMixed) {
    histoTitle = "(--) mixed | Centrality: "; 
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
    
    gHistNN[2] = bMixed->GetCorrelationFunctionNN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistNN[2]->GetYaxis()->SetTitleOffset(1.5);
    gHistNN[2]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistNN[2]->SetTitle(histoTitle.Data());
    cNN[2] = new TCanvas("cNN2","",300,200,600,500);
    cNN[2]->SetFillColor(10); 
    cNN[2]->SetHighLightColor(10);
    gHistNN[2]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();
    pngName = "DeltaPhiDeltaEtaMixed.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativeNegative.png";
    cNN[2]->SaveAs(pngName.Data());

    //Correlation function (--)
    gHistNN[3] = dynamic_cast<TH2D *>(gHistNN[0]->Clone());
    gHistNN[3]->Divide(gHistNN[2]);
    gHistNN[3]->GetXaxis()->SetRangeUser(-1.5,1.5);
    gHistNN[3]->GetZaxis()->SetTitle("C_{--}(#Delta#eta,#Delta#varphi)");
    cNN[3] = new TCanvas("cNN3","",300,300,600,500);
    cNN[3]->SetFillColor(10); 
    cNN[3]->SetHighLightColor(10);
    gHistNN[3]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "CorrelationFunction.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativeNegative.png";
    cNN[3]->SaveAs(pngName.Data());
  }
}

