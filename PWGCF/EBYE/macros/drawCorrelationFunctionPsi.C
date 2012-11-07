const Int_t numberOfCentralityBins = 8;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};

const Int_t gRebin = 1;

//____________________________________________________________//
void drawCorrelationFunctionPsiAllPtCombinations(const char* filename = "AnalysisResults.root", 
						 Int_t gCentrality = 1,
						 Int_t gBit = -1,
						 const char* gCentralityEstimator = 0x0,
						 Bool_t kShowShuffled = kFALSE, 
						 Bool_t kShowMixed = kTRUE, 
						 Double_t psiMin = -0.5, 
						 Double_t psiMax = 3.5){

  // this could also be retrieved directly from AliBalancePsi
  const Int_t kNPtBins = 16;
  Double_t ptBins[kNPtBins+1] = {0.2,0.6,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.,12.,15.,20.};

  cout<<"You have chosen to do all pT combinations --> this could take some time."<<endl;

  for(Int_t iTrig = 0; iTrig < 2/*kNPtBins*/; iTrig++){
    for(Int_t iAssoc = 0; iAssoc < 2/*kNPtBins*/; iAssoc++){
      cout<<"================================================================="<<endl;
      cout<<"PROCESS NOW: "<<endl; 
      cout<<" -> "<< ptBins[iTrig]<<" < pTtrig < "<<ptBins[iTrig+1]<<"   "<<ptBins[iAssoc]<<" < pTassoc < "<<ptBins[iAssoc+1]<<endl;
      drawCorrelationFunctionPsi(filename,gCentrality,gBit,gCentralityEstimator,kShowShuffled,kShowMixed,psiMin,psiMax,ptBins[iTrig],ptBins[iTrig+1],ptBins[iAssoc],ptBins[iAssoc+1]);
      cout<<"================================================================="<<endl;
      cout<<endl;
    }
  }

}

//____________________________________________________________//
void drawCorrelationFunctionsAllPtCombinations(const char* lhcPeriod = "LHC11h",
					       Int_t gTrainID = 171,			      
					       Int_t gCentrality = 1,
					       Double_t psiMin = -0.5, Double_t psiMax = 3.5) 
{

 // this could also be retrieved directly from AliBalancePsi
  const Int_t kNPtBins = 16;
  Double_t ptBins[kNPtBins+1] = {0.2,0.6,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.,12.,15.,20.};

  cout<<"You have chosen to do all pT combinations --> this could take some time."<<endl;

  for(Int_t iTrig = 0; iTrig < 2/*kNPtBins*/; iTrig++){
    for(Int_t iAssoc = 0; iAssoc < 2/*kNPtBins*/; iAssoc++){
      cout<<"================================================================="<<endl;
      cout<<"FIT NOW: "<<endl; 
      cout<<" -> "<< ptBins[iTrig]<<" < pTtrig < "<<ptBins[iTrig+1]<<"   "<<ptBins[iAssoc]<<" < pTassoc < "<<ptBins[iAssoc+1]<<endl;
      drawCorrelationFunctions(lhcPeriod,gTrainID,gCentrality,psiMin,psiMax,ptBins[iTrig],ptBins[iTrig+1],ptBins[iAssoc],ptBins[iAssoc+1]);
      cout<<"================================================================="<<endl;
      cout<<endl;
    }
  }  

}


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
  //cPN[0]->SaveAs(pngName.Data());
  
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
    //cPN[1]->SaveAs(pngName.Data());
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
    //cPN[2]->SaveAs(pngName.Data());

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
    //cPN[3]->SaveAs(pngName.Data());
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
  //cNP[0]->SaveAs(pngName.Data());

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
    //cNP[1]->SaveAs(pngName.Data());
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
    //cNP[2]->SaveAs(pngName.Data());

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
    //cNP[3]->SaveAs(pngName.Data());
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
  //cPP[0]->SaveAs(pngName.Data());
  
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
    //cPP[1]->SaveAs(pngName.Data());
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
    //cPP[2]->SaveAs(pngName.Data());

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
    //cPP[3]->SaveAs(pngName.Data());
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
  //cNN[0]->SaveAs(pngName.Data());

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
    //cNN[1]->SaveAs(pngName.Data());
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
    //cNN[2]->SaveAs(pngName.Data());

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
    //cNN[3]->SaveAs(pngName.Data());
  }

  //Write to output file
  TString newFileName = "correlationFunction.Centrality";  
  newFileName += gCentrality; newFileName += ".Psi";
  if((psiMin == -0.5)&&(psiMax == 0.5)) newFileName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) newFileName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) newFileName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) newFileName += "Rest.PttFrom";
  else newFileName += "All.PttFrom";
  newFileName += Form("%.1f",ptTriggerMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptTriggerMax); newFileName += "PtaFrom";
  newFileName += Form("%.1f",ptAssociatedMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptAssociatedMax); 
  newFileName += ".root";
  TFile *newFile = TFile::Open(newFileName.Data(),"recreate");
  gHistPN[0]->SetName("gHistPNRaw"); gHistPN[0]->Write();
  gHistNP[0]->SetName("gHistNPRaw"); gHistNP[0]->Write();
  gHistPP[0]->SetName("gHistPPRaw"); gHistPP[0]->Write();
  gHistNN[0]->SetName("gHistNNRaw"); gHistNN[0]->Write();
  if(listBFShuffled) {
    gHistPN[1]->SetName("gHistPNShuffled"); gHistPN[1]->Write();
    gHistNP[1]->SetName("gHistNPShuffled"); gHistNP[1]->Write();
    gHistPP[1]->SetName("gHistPPShuffled"); gHistPP[1]->Write();
    gHistNN[1]->SetName("gHistNNShuffled"); gHistNN[1]->Write();
  }
  if(listBFMixed) {
    gHistPN[2]->SetName("gHistPNMixed"); gHistPN[2]->Write();
    gHistNP[2]->SetName("gHistNPMixed"); gHistNP[2]->Write();
    gHistPP[2]->SetName("gHistPPMixed"); gHistPP[2]->Write();
    gHistNN[2]->SetName("gHistNNMixed"); gHistNN[2]->Write();

    gHistPN[3]->SetName("gHistPNCorrelationFunctions"); gHistPN[3]->Write();
    gHistNP[3]->SetName("gHistNPCorrelationFunctions"); gHistNP[3]->Write();
    gHistPP[3]->SetName("gHistPPCorrelationFunctions"); gHistPP[3]->Write();
    gHistNN[3]->SetName("gHistNNCorrelationFunctions"); gHistNN[3]->Write();
  }
  newFile->Close();

  // some cleaning
  for(Int_t i = 0; i < 4; i++){

    if(!listBFShuffled && i == 1) continue;
    if(!listBFMixed && (i == 2 || i == 3)) continue;

    if(gHistPP[i]) delete gHistPP[i];
    if(gHistPN[i]) delete gHistPN[i];
    if(gHistNP[i]) delete gHistNP[i];
    if(gHistNN[i]) delete gHistNN[i];
    
    if(cPN[i]) delete cPN[i];
    if(cNP[i]) delete cNP[i];
    if(cPP[i]) delete cPP[i];
    if(cNN[i]) delete cNN[i];
  }

  delete hP;
  delete hN;
  delete hPP;
  delete hPN;
  delete hNP;
  delete hNN;

  delete hPMixed;
  delete hNMixed;
  delete hPPMixed;
  delete hPNMixed;
  delete hNPMixed;
  delete hNNMixed;

  delete hPShuffled;
  delete hNShuffled;
  delete hPPShuffled;
  delete hPNShuffled;
  delete hNPShuffled;
  delete hNNShuffled;

}

//____________________________________________________________//
void drawCorrelationFunctions(const char* lhcPeriod = "LHC11h",
			      Int_t gTrainID = 171,			      
			      Int_t gCentrality = 1,
			      Double_t psiMin = -0.5, Double_t psiMax = 3.5,
			      Double_t ptTriggerMin = -1.,
			      Double_t ptTriggerMax = -1.,
			      Double_t ptAssociatedMin = -1.,
			      Double_t ptAssociatedMax = -1.) {
  //Macro that draws the charge dependent correlation functions
  //for each centrality bin for the different pT of trigger and 
  //associated particles
  //Author: Panos.Christakoglou@nikhef.nl
  TGaxis::SetMaxDigits(3);

  //Get the input file
  TString filename = "PbPb/"; filename += lhcPeriod; 
  filename +="/Train"; filename += gTrainID;
  filename +="/Centrality"; filename += gCentrality;
  filename += "/correlationFunction.Centrality";
  filename += gCentrality; filename += ".Psi";
  if((psiMin == -0.5)&&(psiMax == 0.5)) filename += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) filename += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) filename += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) filename += "Rest.Ptt";
  else filename += "All.PttFrom";
  filename += Form("%.1f",ptTriggerMin); filename += "To"; 
  filename += Form("%.1f",ptTriggerMax); filename += "PtaFrom";
  filename += Form("%.1f",ptAssociatedMin); filename += "To"; 
  filename += Form("%.1f",ptAssociatedMax); 
  filename += ".root";  

  //Open the file
  TFile *f = TFile::Open(filename.Data());
  if((!f)||(!f->IsOpen())) {
    Printf("The file %s is not found. Aborting...",filename);
    return listBF;
  }
  //f->ls();
  
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

  TString pngName;

  //============================================================//
  //Get the +- correlation function
  TH2D *gHistPN = dynamic_cast<TH2D *>(f->Get("gHistPNCorrelationFunctions"));
  gHistPN->SetStats(kFALSE);
  gHistPN->SetTitle("");
  gHistPN->GetXaxis()->SetRangeUser(-1.45,1.45);
  gHistPN->GetXaxis()->CenterTitle();
  gHistPN->GetXaxis()->SetTitleOffset(1.2);
  gHistPN->GetYaxis()->CenterTitle();
  gHistPN->GetYaxis()->SetTitleOffset(1.2);
  gHistPN->GetZaxis()->SetTitleOffset(1.2);
  TCanvas *cPN = new TCanvas("cPN","",0,0,600,500);
  cPN->SetFillColor(10); cPN->SetHighLightColor(10);
  cPN->SetLeftMargin(0.15);
  gHistPN->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.44,0.88,centralityLatex.Data());
  //latexInfo1->DrawLatex(0.44,0.82,psiLatex.Data());
  latexInfo1->DrawLatex(0.44,0.82,pttLatex.Data());
  latexInfo1->DrawLatex(0.44,0.76,ptaLatex.Data());

  pngName = "CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if((psiMin == -0.5)&&(psiMax == 0.5)) pngName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) pngName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) pngName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) pngName += "Rest.Ptt";
  else pngName += "All.PttFrom";
  pngName += Form("%.1f",ptTriggerMin); pngName += "To"; 
  pngName += Form("%.1f",ptTriggerMax); pngName += "PtaFrom";
  pngName += Form("%.1f",ptAssociatedMin); pngName += "To"; 
  pngName += Form("%.1f",ptAssociatedMax); 
  pngName += ".PositiveNegative.png";
  cPN->SaveAs(pngName.Data());
  fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			  ptTriggerMin,ptTriggerMax,
			  ptAssociatedMin, ptAssociatedMax,gHistPN);
  //============================================================//
  //Get the -+ correlation function
  TH2D *gHistNP = dynamic_cast<TH2D *>(f->Get("gHistNPCorrelationFunctions"));
  gHistNP->SetStats(kFALSE);
  gHistNP->SetTitle("");
  gHistNP->GetXaxis()->SetRangeUser(-1.45,1.45);
  gHistNP->GetXaxis()->CenterTitle();
  gHistNP->GetXaxis()->SetTitleOffset(1.2);
  gHistNP->GetYaxis()->CenterTitle();
  gHistNP->GetYaxis()->SetTitleOffset(1.2);
  gHistNP->GetZaxis()->SetTitleOffset(1.2);
  TCanvas *cNP = new TCanvas("cNP","",50,50,600,500);
  cNP->SetFillColor(10); cNP->SetHighLightColor(10);
  cNP->SetLeftMargin(0.15);
  gHistNP->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.44,0.88,centralityLatex.Data());
  //latexInfo1->DrawLatex(0.44,0.82,psiLatex.Data());
  latexInfo1->DrawLatex(0.44,0.82,pttLatex.Data());
  latexInfo1->DrawLatex(0.44,0.76,ptaLatex.Data());

  pngName = "CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if((psiMin == -0.5)&&(psiMax == 0.5)) pngName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) pngName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) pngName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) pngName += "Rest.Ptt";
  else pngName += "All.PttFrom";
  pngName += Form("%.1f",ptTriggerMin); pngName += "To"; 
  pngName += Form("%.1f",ptTriggerMax); pngName += "PtaFrom";
  pngName += Form("%.1f",ptAssociatedMin); pngName += "To"; 
  pngName += Form("%.1f",ptAssociatedMax); 
  pngName += ".NegativePositive.png";
  cNP->SaveAs(pngName.Data());

  //============================================================//
  //Get the ++ correlation function
  TH2D *gHistPP = dynamic_cast<TH2D *>(f->Get("gHistPPCorrelationFunctions"));
  gHistPP->SetStats(kFALSE);
  gHistPP->SetTitle("");
  gHistPP->GetXaxis()->SetRangeUser(-1.45,1.45);
  gHistPP->GetXaxis()->CenterTitle();
  gHistPP->GetXaxis()->SetTitleOffset(1.2);
  gHistPP->GetYaxis()->CenterTitle();
  gHistPP->GetYaxis()->SetTitleOffset(1.2);
  gHistPP->GetZaxis()->SetTitleOffset(1.2);
  TCanvas *cPP = new TCanvas("cPP","",100,100,600,500);
  cPP->SetFillColor(10); cPP->SetHighLightColor(10);
  cPP->SetLeftMargin(0.15);
  gHistPP->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.44,0.88,centralityLatex.Data());
  //latexInfo1->DrawLatex(0.44,0.82,psiLatex.Data());
  latexInfo1->DrawLatex(0.44,0.82,pttLatex.Data());
  latexInfo1->DrawLatex(0.44,0.76,ptaLatex.Data());

  pngName = "CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if((psiMin == -0.5)&&(psiMax == 0.5)) pngName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) pngName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) pngName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) pngName += "Rest.Ptt";
  else pngName += "All.PttFrom";
  pngName += Form("%.1f",ptTriggerMin); pngName += "To"; 
  pngName += Form("%.1f",ptTriggerMax); pngName += "PtaFrom";
  pngName += Form("%.1f",ptAssociatedMin); pngName += "To"; 
  pngName += Form("%.1f",ptAssociatedMax); 
  pngName += ".PositivePositive.png";
  cPP->SaveAs(pngName.Data());

  //============================================================//
  //Get the -- correlation function
  TH2D *gHistNN = dynamic_cast<TH2D *>(f->Get("gHistNNCorrelationFunctions"));
  gHistNN->SetStats(kFALSE);
  gHistNN->SetTitle("");
  gHistNN->GetXaxis()->SetRangeUser(-1.45,1.45);
  gHistNN->GetXaxis()->CenterTitle();
  gHistNN->GetXaxis()->SetTitleOffset(1.2);
  gHistNN->GetYaxis()->CenterTitle();
  gHistNN->GetYaxis()->SetTitleOffset(1.2);
  gHistNN->GetZaxis()->SetTitleOffset(1.2);
  TCanvas *cNN = new TCanvas("cNN","",150,150,600,500);
  cNN->SetFillColor(10); cNN->SetHighLightColor(10);
  cNN->SetLeftMargin(0.15);
  gHistNN->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.44,0.88,centralityLatex.Data());
  //latexInfo1->DrawLatex(0.44,0.82,psiLatex.Data());
  latexInfo1->DrawLatex(0.44,0.82,pttLatex.Data());
  latexInfo1->DrawLatex(0.44,0.76,ptaLatex.Data());

  pngName = "CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if((psiMin == -0.5)&&(psiMax == 0.5)) pngName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) pngName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) pngName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) pngName += "Rest.Ptt";
  else pngName += "All.PttFrom";
  pngName += Form("%.1f",ptTriggerMin); pngName += "To"; 
  pngName += Form("%.1f",ptTriggerMax); pngName += "PtaFrom";
  pngName += Form("%.1f",ptAssociatedMin); pngName += "To"; 
  pngName += Form("%.1f",ptAssociatedMax); 
  pngName += ".NegativeNegative.png";
  cNN->SaveAs(pngName.Data());
}

//____________________________________________________________//
void fitCorrelationFunctions(Int_t gCentrality = 1,
			     Double_t psiMin = -0.5, Double_t psiMax = 3.5,
			     Double_t ptTriggerMin = -1.,
			     Double_t ptTriggerMax = -1.,
			     Double_t ptAssociatedMin = -1.,
			     Double_t ptAssociatedMax = -1.,
			     TH2D *gHist) {

  cout<<"FITTING FUNCTION"<<endl;

  //near side peak: [1]*TMath::Exp(-TMath::Power((0.5*TMath::Power((x/[2]),2)+0.5*TMath::Power((y/[3]),2)),[4]))
  //away side ridge: [5]*TMath::Exp(-TMath::Power((0.5*TMath::Power(((y-TMath::Pi())/[6]),2)),[7]))
  //longitudinal ridge: [8]*TMath::Exp(-TMath::Power((0.5*TMath::Power((x/[9]),2)),[10]))
  //wing structures: [11]*TMath::Power(x,2)
  //flow contribution (v1 up to v4): 2.*([12]*TMath::Cos(y) + [13]*TMath::Cos(2.*y) + [14]*TMath::Cos(3.*y) + [15]*TMath::Cos(4.*y))
  TF2 *gFitFunction = new TF2("gFitFunction","[0]+[1]*TMath::Exp(-TMath::Power((0.5*TMath::Power((x/[2]),2)+0.5*TMath::Power((y/[3]),2)),[4]))+[5]*TMath::Exp(-TMath::Power((0.5*TMath::Power(((y-TMath::Pi())/[6]),2)),[7]))+[8]*TMath::Exp(-TMath::Power((0.5*TMath::Power((x/[9]),2)),[10]))+[11]*TMath::Power(x,2)+2.*[12]*([13]*TMath::Cos(y) + [14]*TMath::Cos(2.*y) + [15]*TMath::Cos(3.*y) + [16]*TMath::Cos(4.*y))",-2.0,2.0,-TMath::Pi()/2.,3.*TMath::Pi()/2.); 
  gFitFunction->SetName("gFitFunction");
  //Normalization
  gFitFunction->SetParName(0,"N1"); gFitFunction->SetParameter(0,1.0);
  //near side peak
  gFitFunction->SetParName(1,"N_{near side}"); gFitFunction->SetParameter(1,0.3);
  gFitFunction->SetParName(2,"Sigma_{near side}(delta eta)"); gFitFunction->SetParameter(2,0.3);
  gFitFunction->SetParName(3,"Sigma_{near side}(delta phi)"); gFitFunction->SetParameter(3,0.1);
  gFitFunction->SetParName(4,"Exponent_{near side}"); gFitFunction->SetParameter(4,1.1);
  //away side ridge
  gFitFunction->SetParName(5,"N_{away side}"); gFitFunction->SetParameter(5,0.1);
  gFitFunction->SetParName(6,"Sigma_{away side}(delta phi)"); gFitFunction->SetParameter(6,1.1);
  gFitFunction->SetParName(7,"Exponent_{away side}"); gFitFunction->SetParameter(7,1.0);
  //longitudianl ridge
  gFitFunction->SetParName(8,"N_{long. ridge}"); gFitFunction->SetParameter(8,0.05);
  gFitFunction->SetParName(9,"Sigma_{long. ridge}(delta eta)"); gFitFunction->SetParameter(9,0.6);
  gFitFunction->SetParName(10,"Exponent_{long. ridge}"); gFitFunction->SetParameter(10,1.0);
  //wing structures
  gFitFunction->SetParName(11,"N_{wing}"); gFitFunction->SetParameter(11,0.01);
  //flow contribution
  gFitFunction->SetParName(12,"N_{flow}"); gFitFunction->SetParameter(12,0.2);
  gFitFunction->SetParName(13,"V1"); gFitFunction->SetParameter(13,0.005);
  gFitFunction->SetParName(14,"V2"); gFitFunction->SetParameter(14,0.1);
  gFitFunction->SetParName(15,"V3"); gFitFunction->SetParameter(15,0.05);
  gFitFunction->SetParName(16,"V4"); gFitFunction->SetParameter(16,0.005);  

  //Fitting the correlation function
  gHist->Fit("gFitFunction","nm");

  //Cloning the histogram
  TString histoName = gHist->GetName(); histoName += "Fit"; 
  TH2D *gHistFit = new TH2D(histoName.Data(),";#Delta#eta;#Delta#varphi (rad);C(#Delta#eta,#Delta#varphi)",gHist->GetNbinsX(),gHist->GetXaxis()->GetXmin(),gHist->GetXaxis()->GetXmax(),gHist->GetNbinsY(),gHist->GetYaxis()->GetXmin(),gHist->GetYaxis()->GetXmax());
  TH2D *gHistResidual = dynamic_cast<TH2D *>(gHist->Clone());
  gHistResidual->SetName("gHistResidual");
  gHistResidual->Sumw2();

  for(Int_t iBinDeltaEta = 1; iBinDeltaEta <= gHist->GetNbinsX(); iBinDeltaEta++) {
    for(Int_t iBinDeltaPhi = 1; iBinDeltaPhi <= gHist->GetNbinsY(); iBinDeltaPhi++) {
      gHistFit->SetBinContent(iBinDeltaEta,iBinDeltaPhi,gFitFunction->Eval(gHist->GetXaxis()->GetBinCenter(iBinDeltaEta),gHist->GetYaxis()->GetBinCenter(iBinDeltaPhi)));
    }
  }
  gHistResidual->Add(gHistFit,-1);

  //Write to output file
  TString newFileName = "correlationFunctionFit";
  if(histoName.Contains("PN")) newFileName += "PN";
  else if(histoName.Contains("NP")) newFileName += "NP";
  else if(histoName.Contains("PP")) newFileName += "PP";
  else if(histoName.Contains("NN")) newFileName += "NN";
  newFileName += ".Centrality";  
  newFileName += gCentrality; newFileName += ".Psi";
  if((psiMin == -0.5)&&(psiMax == 0.5)) newFileName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) newFileName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) newFileName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) newFileName += "Rest.PttFrom";
  else newFileName += "All.PttFrom";
  newFileName += Form("%.1f",ptTriggerMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptTriggerMax); newFileName += "PtaFrom";
  newFileName += Form("%.1f",ptAssociatedMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptAssociatedMax); 
  newFileName += ".root";
  TFile *newFile = TFile::Open(newFileName.Data(),"recreate");
  gHistFit->Write();
  gHistResidual->Write();
  gFitFunction->Write();
  newFile->Close();
  

}
