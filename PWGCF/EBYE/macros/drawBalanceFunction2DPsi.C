const Int_t numberOfCentralityBins = 12;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","0-100","0-1","1-2","2-3"};


const Int_t gRebin = 1;
void drawBalanceFunction2DPsi(const char* filename = "AnalysisResultsPsi.root", 
			      Int_t gCentrality = 1,
			      Int_t gBit = -1,
			      const char* gCentralityEstimator = 0x0,
			      Bool_t kShowShuffled = kFALSE, 
			      Bool_t kShowMixed = kTRUE, 
			      Double_t psiMin = -0.5, Double_t psiMax = 0.5,
			      Double_t vertexZMin = -10.,
			      Double_t vertexZMax = 10.,
			      Double_t ptTriggerMin = -1.,
			      Double_t ptTriggerMax = -1.,
			      Double_t ptAssociatedMin = -1.,
			      Double_t ptAssociatedMax = -1.,
			      Bool_t kUseVzBinning = kTRUE,
			      Bool_t k2pMethod = kTRUE,
			      TString eventClass = "EventPlane") //Can be "EventPlane", "Centrality", "Multiplicity"
{
  //Macro that draws the BF distributions for each centrality bin
  //for reaction plane dependent analysis
  //Author: Panos.Christakoglou@nikhef.nl
  //Load the PWG2ebye library
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGCFebye.so");

  //gROOT->LoadMacro("~/SetPlotStyle.C");
  //SetPlotStyle();
  gStyle->SetPalette(1,0);

  //Prepare the objects and return them
  TList *listBF = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,0);
  TList *listBFShuffled = NULL;
  if(kShowShuffled) listBFShuffled = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,1);
  TList *listBFMixed = NULL;
  if(kShowMixed) listBFMixed = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,2);
  if(!listBF) {
    Printf("The TList object was not created");
    return;
  }
  else 
    draw(listBF,listBFShuffled,listBFMixed,gCentrality,gCentralityEstimator,
	 psiMin,psiMax,vertexZMin,vertexZMax,
	 ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,
	 kUseVzBinning,k2pMethod,eventClass);  
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
    //cout<<"second iteration"<<endl;
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
      histoName = "fHistP";
    else if(kData == 1)
      histoName = "fHistP_shuffle";
    else if(kData == 2)
      histoName = "fHistP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    AliTHn *fHistP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));  
    if(!fHistP) {
      Printf("fHistP %s not found!!!",histoName.Data());
      break;
    }
    fHistP->FillParent(); fHistP->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistN";
    if(kData == 1)
      histoName = "fHistN_shuffle";
    if(kData == 2)
      histoName = "fHistN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    AliTHn *fHistN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistN) {
      Printf("fHistN %s not found!!!",histoName.Data());
      break;
    }
    fHistN->FillParent(); fHistN->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistPN";
    if(kData == 1)
      histoName = "fHistPN_shuffle";
    if(kData == 2)
      histoName = "fHistPN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    AliTHn *fHistPN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistPN) {
      Printf("fHistPN %s not found!!!",histoName.Data());
      break;
    }
    fHistPN->FillParent(); fHistPN->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistNP";
    if(kData == 1)
      histoName = "fHistNP_shuffle";
    if(kData == 2)
      histoName = "fHistNP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    AliTHn *fHistNP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistNP) {
      Printf("fHistNP %s not found!!!",histoName.Data());
      break;
    }
    fHistNP->FillParent(); fHistNP->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistPP";
    if(kData == 1)
      histoName = "fHistPP_shuffle";
    if(kData == 2)
      histoName = "fHistPP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    AliTHn *fHistPP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
    if(!fHistPP) {
      Printf("fHistPP %s not found!!!",histoName.Data());
      break;
    }
    fHistPP->FillParent(); fHistPP->DeleteContainers();
    
    if(kData == 0)
      histoName = "fHistNN";
    if(kData == 1)
      histoName = "fHistNN_shuffle";
    if(kData == 2)
      histoName = "fHistNN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
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
void draw(TList *listBF, TList *listBFShuffled, TList *listBFMixed,
	  Int_t gCentrality, const char* gCentralityEstimator,
	  Double_t psiMin, Double_t psiMax,
	  Double_t vertexZMin,
	  Double_t vertexZMax,
	  Double_t ptTriggerMin, Double_t ptTriggerMax,
	  Double_t ptAssociatedMin, Double_t ptAssociatedMax,
	  Bool_t kUseVzBinning=kFALSE,
	  Bool_t k2pMethod = kFALSE, TString eventClass) {  
  //balance function
  AliTHn *hP = NULL;
  AliTHn *hN = NULL;
  AliTHn *hPN = NULL;
  AliTHn *hNP = NULL;
  AliTHn *hPP = NULL;
  AliTHn *hNN = NULL;
  //listBF->ls();
  //Printf("=================");
  TString histoName = "fHistP";
  if(gCentralityEstimator) histoName += gCentralityEstimator;
  hP = (AliTHn*) listBF->FindObject(histoName.Data());
  hP->SetName("gHistP");
  histoName = "fHistN";
  if(gCentralityEstimator) histoName += gCentralityEstimator;
  hN = (AliTHn*) listBF->FindObject(histoName.Data());
  hN->SetName("gHistN");
  histoName = "fHistPN";
  if(gCentralityEstimator) histoName += gCentralityEstimator;
  hPN = (AliTHn*) listBF->FindObject(histoName.Data());
  hPN->SetName("gHistPN");
  histoName = "fHistNP";
  if(gCentralityEstimator) histoName += gCentralityEstimator;
  hNP = (AliTHn*) listBF->FindObject(histoName.Data());
  hNP->SetName("gHistNP");
  histoName = "fHistPP";
  if(gCentralityEstimator) histoName += gCentralityEstimator;
  hPP = (AliTHn*) listBF->FindObject(histoName.Data());
  hPP->SetName("gHistPP");
  histoName = "fHistNN";
  if(gCentralityEstimator) histoName += gCentralityEstimator;
  hNN = (AliTHn*) listBF->FindObject(histoName.Data());
  hNN->SetName("gHistNN");

  AliBalancePsi *b = new AliBalancePsi();
  b->SetEventClass(eventClass);
  b->SetHistNp(hP);
  b->SetHistNn(hN);
  b->SetHistNpn(hPN);
  b->SetHistNnp(hNP);
  b->SetHistNpp(hPP);
  b->SetHistNnn(hNN);
  if(kUseVzBinning) b->SetVertexZBinning(kTRUE);


  //balance function shuffling
  AliTHn *hPShuffled = NULL;
  AliTHn *hNShuffled = NULL;
  AliTHn *hPNShuffled = NULL;
  AliTHn *hNPShuffled = NULL;
  AliTHn *hPPShuffled = NULL;
  AliTHn *hNNShuffled = NULL;
  if(listBFShuffled) {
    //listBFShuffled->ls();
    histoName = "fHistP_shuffle";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hPShuffled->SetName("gHistPShuffled");
    histoName = "fHistN_shuffle";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hNShuffled->SetName("gHistNShuffled");
    histoName = "fHistPN_shuffle";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPNShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hPNShuffled->SetName("gHistPNShuffled");
    histoName = "fHistNP_shuffle";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNPShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hNPShuffled->SetName("gHistNPShuffled");
    histoName = "fHistPP_shuffle";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPPShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hPPShuffled->SetName("gHistPPShuffled");
    histoName = "fHistNN_shuffle";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNNShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hNNShuffled->SetName("gHistNNShuffled");
    
    AliBalancePsi *bShuffled = new AliBalancePsi();
    bShuffled->SetEventClass(eventClass);
    bShuffled->SetHistNp(hPShuffled);
    bShuffled->SetHistNn(hNShuffled);
    bShuffled->SetHistNpn(hPNShuffled);
    bShuffled->SetHistNnp(hNPShuffled);
    bShuffled->SetHistNpp(hPPShuffled);
    bShuffled->SetHistNnn(hNNShuffled);
  if(kUseVzBinning) bShuffled->SetVertexZBinning(kTRUE);

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
    histoName = "fHistP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPMixed = (AliTHn*) listBFMixed->FindObject(histoName.Data());
    hPMixed->SetName("gHistPMixed");
    histoName = "fHistN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNMixed = (AliTHn*) listBFMixed->FindObject(histoName.Data());
    hNMixed->SetName("gHistNMixed");
    histoName = "fHistPN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPNMixed = (AliTHn*) listBFMixed->FindObject(histoName.Data());
    hPNMixed->SetName("gHistPNMixed");
    histoName = "fHistNP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNPMixed = (AliTHn*) listBFMixed->FindObject(histoName.Data());
    histoName = "fHistNP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNPMixed->SetName("gHistNPMixed");
    histoName = "fHistPP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPPMixed = (AliTHn*) listBFMixed->FindObject(histoName.Data());
    hPPMixed->SetName("gHistPPMixed");
    histoName = "fHistNN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNNMixed = (AliTHn*) listBFMixed->FindObject(histoName.Data());
    hNNMixed->SetName("gHistNNMixed");
    
    AliBalancePsi *bMixed = new AliBalancePsi();
    bMixed->SetEventClass(eventClass);
    bMixed->SetHistNp(hPMixed);
    bMixed->SetHistNn(hNMixed);
    bMixed->SetHistNpn(hPNMixed);
    bMixed->SetHistNnp(hNPMixed);
    bMixed->SetHistNpp(hPPMixed);
    bMixed->SetHistNnn(hNNMixed);
    if(kUseVzBinning) bMixed->SetVertexZBinning(kTRUE);
  
  }

  TH2D *gHistBalanceFunction;
  TH2D *gHistBalanceFunctionSubtracted;
  TH2D *gHistBalanceFunctionShuffled;
  TH2D *gHistBalanceFunctionMixed;
  TString histoTitle, pngName;
  
  if(eventClass == "Centrality"){
    histoTitle = "Centrality: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " % ";
    histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
  }
  else if(eventClass == "Multiplicity"){
    histoTitle = "Multiplicity: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " tracks";
    histoTitle += " (0^{o} < #varphi - #Psi_{2} < 180^{o})"; 
  }
  else{ // "EventPlane" (default)
    histoTitle = "Centrality: ";
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
  }

  if(k2pMethod) 
    if(bMixed)
      gHistBalanceFunction = b->GetBalanceFunctionDeltaEtaDeltaPhi2pMethod(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
    else{
      cerr<<"NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
      return;
    }
  else
    gHistBalanceFunction = b->GetBalanceFunctionDeltaEtaDeltaPhi(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistBalanceFunction->SetTitle(histoTitle.Data());
  gHistBalanceFunction->GetYaxis()->SetTitleOffset(1.3);
  gHistBalanceFunction->SetName("gHistBalanceFunction");

  if(listBFShuffled) {
    
    if(k2pMethod) 
      if(bMixed)
	gHistBalanceFunctionShuffled = bShuffled->GetBalanceFunctionDeltaEtaDeltaPhi2pMethod(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
      else{
	cerr<<"NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
	return;
      }
    else
      gHistBalanceFunctionShuffled = bShuffled->GetBalanceFunctionDeltaEtaDeltaPhi(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistBalanceFunctionShuffled->SetTitle(histoTitle.Data());
    gHistBalanceFunctionShuffled->GetYaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionShuffled->SetName("gHistBalanceFunctionShuffled");
  }

  if(listBFMixed) {
    if(k2pMethod) 
      if(bMixed)
	gHistBalanceFunctionMixed = bMixed->GetBalanceFunctionDeltaEtaDeltaPhi2pMethod(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
      else{
	cerr<<"NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
	return;
      }
    else
      gHistBalanceFunctionMixed = bMixed->GetBalanceFunctionDeltaEtaDeltaPhi(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    gHistBalanceFunctionMixed->SetTitle(histoTitle.Data());
    gHistBalanceFunctionMixed->GetYaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionMixed->SetName("gHistBalanceFunctionMixed");
  
    gHistBalanceFunctionSubtracted = dynamic_cast<TH2D *>(gHistBalanceFunction->Clone());
    gHistBalanceFunctionSubtracted->Add(gHistBalanceFunctionMixed,-1);
    gHistBalanceFunctionSubtracted->SetTitle(histoTitle.Data());
    gHistBalanceFunctionSubtracted->GetYaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionSubtracted->SetName("gHistBalanceFunctionSubtracted");
  }

  //Draw the results
  TCanvas *c1 = new TCanvas("c1","",0,0,600,500);
  c1->SetFillColor(10); 
  c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.15);
  gHistBalanceFunction->DrawCopy("lego2");
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();  
  TCanvas *c1a = new TCanvas("c1a","",600,0,600,500);
  c1a->SetFillColor(10); 
  c1a->SetHighLightColor(10);
  c1a->SetLeftMargin(0.15);
  gHistBalanceFunction->DrawCopy("colz");

  if(listBFShuffled) {
    TCanvas *c2 = new TCanvas("c2","",100,100,600,500);
    c2->SetFillColor(10); 
    c2->SetHighLightColor(10);
    c2->SetLeftMargin(0.15);
    gHistBalanceFunctionShuffled->DrawCopy("lego2");
    gPad->SetTheta(30); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();  
    TCanvas *c2a = new TCanvas("c2a","",700,100,600,500);
    c2a->SetFillColor(10); 
    c2a->SetHighLightColor(10);
    c2a->SetLeftMargin(0.15);
    gHistBalanceFunctionShuffled->DrawCopy("colz");
  }

  if(listBFMixed) {
    TCanvas *c3 = new TCanvas("c3","",200,200,600,500);
    c3->SetFillColor(10); 
    c3->SetHighLightColor(10);
    c3->SetLeftMargin(0.15);
    gHistBalanceFunctionMixed->DrawCopy("lego2");
    gPad->SetTheta(30); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();  
    TCanvas *c3a = new TCanvas("c3a","",800,200,600,500);
    c3a->SetFillColor(10); 
    c3a->SetHighLightColor(10);
    c3a->SetLeftMargin(0.15);
    gHistBalanceFunctionMixed->DrawCopy("colz");

    TCanvas *c4 = new TCanvas("c4","",300,300,600,500);
    c4->SetFillColor(10); 
    c4->SetHighLightColor(10);
    c4->SetLeftMargin(0.15);
    gHistBalanceFunctionSubtracted->DrawCopy("lego2");
    gPad->SetTheta(30); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();  
    TCanvas *c4a = new TCanvas("c4a","",900,300,600,500);
    c4a->SetFillColor(10); 
    c4a->SetHighLightColor(10);
    c4a->SetLeftMargin(0.15);
    gHistBalanceFunctionSubtracted->DrawCopy("colz");

    fitbalanceFunction(gCentrality, psiMin , psiMax,
		       ptTriggerMin, ptTriggerMax,
		       ptAssociatedMin, ptAssociatedMax,
		       gHistBalanceFunctionSubtracted,k2pMethod, eventClass);
  }

  TString newFileName = "balanceFunction2D."; 
  if(eventClass == "Centrality"){
    newFileName += Form("Centrality%.1fTo%.1f",psiMin,psiMax);
    newFileName += ".PsiAll.PttFrom";
  }
  else if(eventClass == "Multiplicity"){
    newFileName += Form("Multiplicity%.0fTo%.0f",psiMin,psiMax);
    newFileName += ".PsiAll.PttFrom";
  }
  else{ // "EventPlane" (default)
    newFileName += "Centrality";
    newFileName += gCentrality; newFileName += ".Psi";
    if((psiMin == -0.5)&&(psiMax == 0.5)) newFileName += "InPlane.Ptt";
    else if((psiMin == 0.5)&&(psiMax == 1.5)) newFileName += "Intermediate.Ptt";
    else if((psiMin == 1.5)&&(psiMax == 2.5)) newFileName += "OutOfPlane.Ptt";
    else if((psiMin == 2.5)&&(psiMax == 3.5)) newFileName += "Rest.PttFrom";
    else newFileName += "All.PttFrom";
  }  
  newFileName += Form("%.1f",ptTriggerMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptTriggerMax); newFileName += "PtaFrom";
  newFileName += Form("%.1f",ptAssociatedMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptAssociatedMax); 
  if(k2pMethod) newFileName += "_2pMethod";
  newFileName += ".root";

  TFile *fOutput = new TFile(newFileName.Data(),"recreate");
  fOutput->cd();
  /*hP->Write(); hN->Write();
  hPN->Write(); hNP->Write();
  hPP->Write(); hNN->Write();
  hPShuffled->Write(); hNShuffled->Write();
  hPNShuffled->Write(); hNPShuffled->Write();
  hPPShuffled->Write(); hNNShuffled->Write();
  hPMixed->Write(); hNMixed->Write();
  hPNMixed->Write(); hNPMixed->Write();
  hPPMixed->Write(); hNNMixed->Write();*/
  gHistBalanceFunction->Write();
  if(listBFShuffled) gHistBalanceFunctionShuffled->Write();
  if(listBFMixed) {
    gHistBalanceFunctionMixed->Write();
    gHistBalanceFunctionSubtracted->Write();
  }
  fOutput->Close();
}

//____________________________________________________________//
void fitbalanceFunction(Int_t gCentrality = 1,
			Double_t psiMin = -0.5, Double_t psiMax = 3.5,
			Double_t ptTriggerMin = -1.,
			Double_t ptTriggerMax = -1.,
			Double_t ptAssociatedMin = -1.,
			Double_t ptAssociatedMax = -1.,
			TH2D *gHist,
			Bool_t k2pMethod = kFALSE, 
			TString eventClass="EventPlane") {
  //balancing charges: [1]*TMath::Exp(-0.5*TMath::Power(((x - [3])/[2]),2)-0.5*TMath::Power(((y - [5])/[4]),2)) 
  //short range correlations: [6]*TMath::Exp(-0.5*TMath::Power(((x - [8])/[7]),2)-0.5*TMath::Power(((y - [10])/[9]),2))
  cout<<"FITTING FUNCTION"<<endl;

  TF2 *gFitFunction = new TF2("gFitFunction","[0] + [1]*TMath::Exp(-0.5*TMath::Power(((x - [3])/[2]),2)-0.5*TMath::Power(((y - [5])/[4]),2)) + [6]*TMath::Exp(-0.5*TMath::Power(((x - [8])/[7]),2)-0.5*TMath::Power(((y - [10])/[9]),2))",-1.2,1.2,-TMath::Pi()/2.,3.*TMath::Pi()/2.); 
  gFitFunction->SetName("gFitFunction");

  //Normalization
  gFitFunction->SetParName(0,"N1"); 
  gFitFunction->SetParameter(0,1.0);

  //2D balance function
  gFitFunction->SetParName(1,"N_{BF}");
  gFitFunction->SetParameter(1,1.0);
  gFitFunction->SetParLimits(1, 0., 100.);
  gFitFunction->SetParName(2,"Sigma_{BF}(delta eta)"); 
  gFitFunction->SetParameter(2,0.6);
  gFitFunction->SetParLimits(2, 0., 1.);
  gFitFunction->SetParName(3,"Mean_{BF}(delta eta)"); 
  gFitFunction->SetParameter(3,0.0);
  gFitFunction->SetParLimits(3, -0.2, 0.2);
  gFitFunction->SetParName(4,"Sigma_{BF}(delta phi)"); 
  gFitFunction->SetParameter(4,0.6);
  gFitFunction->SetParLimits(4, 0., 1.);
  gFitFunction->SetParName(5,"Mean_{BF}(delta phi)"); 
  gFitFunction->SetParameter(5,0.0);
  gFitFunction->SetParLimits(5, -0.2, 0.2);

  //short range structure
  gFitFunction->SetParName(6,"N_{SR}");
  gFitFunction->SetParameter(6,5.0);
  gFitFunction->SetParLimits(6, 0., 100.);
  gFitFunction->SetParName(7,"Sigma_{SR}(delta eta)"); 
  gFitFunction->SetParameter(7,0.01);
  gFitFunction->SetParLimits(7, 0.0, 0.1);
  gFitFunction->SetParName(8,"Mean_{SR}(delta eta)"); 
  gFitFunction->SetParameter(8,0.0);
  gFitFunction->SetParLimits(8, -0.01, 0.01);
  gFitFunction->SetParName(9,"Sigma_{SR}(delta phi)"); 
  gFitFunction->SetParameter(9,0.01);
  gFitFunction->SetParLimits(9, 0.0, 0.1);
  gFitFunction->SetParName(10,"Mean_{SR}(delta phi)"); 
  gFitFunction->SetParameter(10,0.0);
  gFitFunction->SetParLimits(10, -0.01, 0.01);


  //Cloning the histogram
  TH2D *gHistResidual = dynamic_cast<TH2D *>(gHist->Clone());
  gHistResidual->SetName("gHistResidual");
  gHistResidual->Sumw2();

  //Fitting the 2D bf
  for(Int_t iAttempt = 0; iAttempt < 10; iAttempt++) {
    gHist->Fit("gFitFunction","nm");
    for(Int_t iParam = 0; iParam < 11; iParam++) 
      gFitFunction->SetParameter(iParam,gFitFunction->GetParameter(iParam));
  }
  cout<<"======================================================"<<endl;
  cout<<"Fit chi2/ndf: "<<gFitFunction->GetChisquare()/gFitFunction->GetNDF()<<" - chi2: "<<gFitFunction->GetChisquare()<<" - ndf: "<<gFitFunction->GetNDF()<<endl;
  cout<<"======================================================"<<endl;

  //Getting the residual
  gHistResidual->Add(gFitFunction,-1);

  //Write to output file
  TString newFileName = "balanceFunctionFit2D.";
  if(eventClass == "Centrality"){
    newFileName += Form("Centrality%.1fTo%.1f",psiMin,psiMax);
    newFileName += ".PsiAll.PttFrom";
  }
  else if(eventClass == "Multiplicity"){
    newFileName += Form("Multiplicity%.0fTo%.0f",psiMin,psiMax);
    newFileName += ".PsiAll.PttFrom";
  }
  else{ // "EventPlane" (default)
    newFileName += "Centrality";
    newFileName += gCentrality; newFileName += ".Psi";
    if((psiMin == -0.5)&&(psiMax == 0.5)) newFileName += "InPlane.Ptt";
    else if((psiMin == 0.5)&&(psiMax == 1.5)) newFileName += "Intermediate.Ptt";
    else if((psiMin == 1.5)&&(psiMax == 2.5)) newFileName += "OutOfPlane.Ptt";
    else if((psiMin == 2.5)&&(psiMax == 3.5)) newFileName += "Rest.PttFrom";
    else newFileName += "All.PttFrom";
  }  
  newFileName += Form("%.1f",ptTriggerMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptTriggerMax); newFileName += "PtaFrom";
  newFileName += Form("%.1f",ptAssociatedMin); newFileName += "To"; 
  newFileName += Form("%.1f",ptAssociatedMax); 
  if(k2pMethod) newFileName += "_2pMethod";
  newFileName += ".root";
  TFile *newFile = TFile::Open(newFileName.Data(),"recreate");
  gHist->Write();
  gHistResidual->Write();
  gFitFunction->Write();
  newFile->Close();
}

//____________________________________________________________//
void drawBFPsi2D(const char* lhcPeriod = "LHC11h",
		 const char* gCentralityEstimator = "V0M",
		 Int_t gBit = 128,
		 const char* gEventPlaneEstimator = "VZERO",
		 Int_t gCentrality = 1,
		 Bool_t kShowShuffled = kFALSE, 
		 Bool_t kShowMixed = kFALSE, 
		 Double_t psiMin = -0.5, Double_t psiMax = 0.5,
		 Double_t ptTriggerMin = -1.,
		 Double_t ptTriggerMax = -1.,
		 Double_t ptAssociatedMin = -1.,
		 Double_t ptAssociatedMax = -1.,
		 Bool_t k2pMethod = kTRUE) {
  //Macro that draws the BF distributions for each centrality bin
  //for reaction plane dependent analysis
  //Author: Panos.Christakoglou@nikhef.nl
  TGaxis::SetMaxDigits(3);

  //Get the input file
  TString filename = lhcPeriod; 
  filename += "/Centrality"; filename += gCentralityEstimator;
  filename += "_Bit"; filename += gBit;
  filename += "_"; filename += gEventPlaneEstimator;
  filename +="/PttFrom";
  filename += Form("%.1f",ptTriggerMin); filename += "To"; 
  filename += Form("%.1f",ptTriggerMax); filename += "PtaFrom";
  filename += Form("%.1f",ptAssociatedMin); filename += "To"; 
  filename += Form("%.1f",ptAssociatedMax); 
  filename += "/balanceFunction2D.Centrality"; 
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
  if(k2pMethod) filename += "_2pMethod";
  filename += ".root";  

  //Open the file
  TFile *f = TFile::Open(filename.Data());
  if((!f)||(!f->IsOpen())) {
    Printf("The file %s is not found. Aborting...",filename);
    return listBF;
  }
  //f->ls();
  
  //Raw balance function
  TH1D *gHistBalanceFunction = dynamic_cast<TH1D *>(f->Get("gHistBalanceFunction"));
  gHistBalanceFunction->SetStats(kFALSE);
  gHistBalanceFunction->GetXaxis()->SetNdivisions(10);
  gHistBalanceFunction->GetYaxis()->SetNdivisions(10);
  gHistBalanceFunction->GetZaxis()->SetNdivisions(10);
  gHistBalanceFunction->GetXaxis()->SetTitleOffset(1.3);
  gHistBalanceFunction->GetYaxis()->SetTitleOffset(1.3);
  gHistBalanceFunction->GetZaxis()->SetTitleOffset(1.3);
  gHistBalanceFunction->GetXaxis()->SetTitle("#Delta #eta");
  gHistBalanceFunction->GetYaxis()->SetTitle("#Delta #varphi (rad)");
  gHistBalanceFunction->GetZaxis()->SetTitle("B(#Delta #eta, #Delta #varphi)");

  //Shuffled balance function
  if(kShowShuffled) {
    TH1D *gHistBalanceFunctionShuffled = dynamic_cast<TH1D *>(f->Get("gHistBalanceFunctionShuffled"));
    gHistBalanceFunctionShuffled->SetStats(kFALSE);
    gHistBalanceFunctionShuffled->GetXaxis()->SetNdivisions(10);
    gHistBalanceFunctionShuffled->GetYaxis()->SetNdivisions(10);
    gHistBalanceFunctionShuffled->GetZaxis()->SetNdivisions(10);
    gHistBalanceFunctionShuffled->GetXaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionShuffled->GetYaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionShuffled->GetZaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionShuffled->GetXaxis()->SetTitle("#Delta #eta");
    gHistBalanceFunctionShuffled->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistBalanceFunctionShuffled->GetZaxis()->SetTitle("B(#Delta #eta, #Delta #varphi)");
  }

  //Mixed balance function
  if(kShowMixed) {
    TH1D *gHistBalanceFunctionMixed = dynamic_cast<TH1D *>(f->Get("gHistBalanceFunctionMixed"));
    gHistBalanceFunctionMixed->SetStats(kFALSE);
    gHistBalanceFunctionMixed->GetXaxis()->SetNdivisions(10);
    gHistBalanceFunctionMixed->GetYaxis()->SetNdivisions(10);
    gHistBalanceFunctionMixed->GetZaxis()->SetNdivisions(10);
    gHistBalanceFunctionMixed->GetXaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionMixed->GetYaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionMixed->GetZaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionMixed->GetXaxis()->SetTitle("#Delta #eta");
    gHistBalanceFunctionMixed->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistBalanceFunctionMixed->GetZaxis()->SetTitle("B(#Delta #eta, #Delta #varphi)");
  }

  //Subtracted balance function
  if(kShowMixed) {
    TH1D *gHistBalanceFunctionSubtracted = dynamic_cast<TH1D *>(f->Get("gHistBalanceFunctionSubtracted"));
    gHistBalanceFunctionSubtracted->SetStats(kFALSE);
    gHistBalanceFunctionSubtracted->GetXaxis()->SetNdivisions(10);
    gHistBalanceFunctionSubtracted->GetYaxis()->SetNdivisions(10);
    gHistBalanceFunctionSubtracted->GetZaxis()->SetNdivisions(10);
    gHistBalanceFunctionSubtracted->GetXaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionSubtracted->GetYaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionSubtracted->GetZaxis()->SetTitleOffset(1.3);
    gHistBalanceFunctionSubtracted->GetXaxis()->SetTitle("#Delta #eta");
    gHistBalanceFunctionSubtracted->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHistBalanceFunctionSubtracted->GetZaxis()->SetTitle("B(#Delta #eta, #Delta #varphi)");
  }

  TString pngName;
  
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
  pttLatex += " < p_{T}^{t} < "; pttLatex += Form("%.1f",ptTriggerMax);
  pttLatex += " GeV/c";

  TString ptaLatex = Form("%.1f",ptAssociatedMin);
  ptaLatex += " < p_{T}^{a} < "; ptaLatex += Form("%.1f",ptAssociatedMax);
  ptaLatex += " GeV/c";

  TLatex *latexInfo1 = new TLatex();
  latexInfo1->SetNDC();
  latexInfo1->SetTextSize(0.045);
  latexInfo1->SetTextColor(1);

  //Draw the results
  TCanvas *c1 = new TCanvas("c1","Raw balance function 2D",0,0,600,500);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.17); c1->SetTopMargin(0.05);
  gHistBalanceFunction->SetTitle("");
  gHistBalanceFunction->GetYaxis()->SetTitleOffset(1.4);
  gHistBalanceFunction->GetYaxis()->SetNdivisions(10);
  gHistBalanceFunction->GetXaxis()->SetRangeUser(-1.4,1.4); 
  gHistBalanceFunction->GetXaxis()->SetNdivisions(10);
  gHistBalanceFunction->GetYaxis()->SetTitle("#Delta #varphi (rad)");
  gHistBalanceFunction->DrawCopy("lego2");
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();  

  latexInfo1->DrawLatex(0.64,0.88,centralityLatex.Data());
  latexInfo1->DrawLatex(0.64,0.82,psiLatex.Data());
  latexInfo1->DrawLatex(0.64,0.76,pttLatex.Data());
  latexInfo1->DrawLatex(0.64,0.70,ptaLatex.Data());

  TString pngName = "BalanceFunction2D."; 
  pngName += "Centrality";
  pngName += gCentrality; 
  if((psiMin == -0.5)&&(psiMax == 0.5)) pngName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) pngName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) pngName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) pngName += "Rest.PttFrom";
  else pngName += "All.PttFrom";  
  pngName += Form("%.1f",ptTriggerMin); pngName += "To"; 
  pngName += Form("%.1f",ptTriggerMax); pngName += "PtaFrom";
  pngName += Form("%.1f",ptAssociatedMin); pngName += "To"; 
  pngName += Form("%.1f",ptAssociatedMax); 
  if(k2pMethod) pngName += "_2pMethod";
  pngName += ".png";

  c1->SaveAs(pngName.Data());

  if(kShowShuffled) {
    TCanvas *c2 = new TCanvas("c2","Shuffled balance function 2D",100,100,600,500);
    c2->SetFillColor(10); c2->SetHighLightColor(10);
    c2->SetLeftMargin(0.17); c2->SetTopMargin(0.05);
    gHistBalanceFunctionShuffled->SetTitle("Shuffled events");
    gHistBalanceFunctionShuffled->GetYaxis()->SetTitleOffset(1.4);
    gHistBalanceFunctionShuffled->GetYaxis()->SetNdivisions(10);
    gHistBalanceFunctionShuffled->GetXaxis()->SetNdivisions(10);
    gHistBalanceFunctionShuffled->DrawCopy("lego2");
    gPad->SetTheta(30); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();  

    latexInfo1->DrawLatex(0.64,0.88,centralityLatex.Data());
    latexInfo1->DrawLatex(0.64,0.82,psiLatex.Data());
    latexInfo1->DrawLatex(0.64,0.76,pttLatex.Data());
    latexInfo1->DrawLatex(0.64,0.70,ptaLatex.Data());
  }

  if(kShowMixed) {
    TCanvas *c3 = new TCanvas("c3","Mixed balance function 2D",200,200,600,500);
    c3->SetFillColor(10); c3->SetHighLightColor(10);
    c3->SetLeftMargin(0.17); c3->SetTopMargin(0.05);
    gHistBalanceFunctionMixed->SetTitle("Mixed events");
    gHistBalanceFunctionMixed->GetYaxis()->SetTitleOffset(1.4);
    gHistBalanceFunctionMixed->GetYaxis()->SetNdivisions(10);
    gHistBalanceFunctionMixed->GetXaxis()->SetNdivisions(10);
    gHistBalanceFunctionMixed->DrawCopy("lego2");
    gPad->SetTheta(30); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();  

    latexInfo1->DrawLatex(0.64,0.88,centralityLatex.Data());
    latexInfo1->DrawLatex(0.64,0.82,psiLatex.Data());
    latexInfo1->DrawLatex(0.64,0.76,pttLatex.Data());
    latexInfo1->DrawLatex(0.64,0.70,ptaLatex.Data());
  }

  if(kShowMixed) {
    TCanvas *c4 = new TCanvas("c4","Subtracted balance function 2D",300,300,600,500);
    c4->SetFillColor(10); c4->SetHighLightColor(10);
    c4->SetLeftMargin(0.17); c4->SetTopMargin(0.05);
    gHistBalanceFunctionSubtracted->SetTitle("Subtracted balance function");
    gHistBalanceFunctionSubtracted->GetYaxis()->SetTitleOffset(1.4);
    gHistBalanceFunctionSubtracted->GetYaxis()->SetNdivisions(10);
    gHistBalanceFunctionSubtracted->GetXaxis()->SetNdivisions(10);
    gHistBalanceFunctionSubtracted->DrawCopy("lego2");
    gPad->SetTheta(30); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();  

    latexInfo1->DrawLatex(0.64,0.88,centralityLatex.Data());
    latexInfo1->DrawLatex(0.64,0.82,psiLatex.Data());
    latexInfo1->DrawLatex(0.64,0.76,pttLatex.Data());
    latexInfo1->DrawLatex(0.64,0.70,ptaLatex.Data());
  }
}

//____________________________________________________________//
void drawProjections(const char* lhcPeriod = "LHC10h",
		     const char* gCentralityEstimator = "V0M",
		     Int_t gBit = 128,
		     const char* gEventPlaneEstimator = "VZERO",
		     Bool_t kProjectInEta = kFALSE,
		     Int_t binMin = 1,
		     Int_t binMax = 80,
		     Int_t gCentrality = 1,
		     Double_t psiMin = -0.5, 
		     Double_t psiMax = 3.5,
		     Double_t vertexZMin = -10., 
		     Double_t vertexZMax = 10.,
		     Double_t ptTriggerMin = -1.,
		     Double_t ptTriggerMax = -1.,
		     Double_t ptAssociatedMin = -1.,
		     Double_t ptAssociatedMax = -1.,
		     Bool_t kUseZYAM = kFALSE,
		     Bool_t k2pMethod = kTRUE,
		     TString eventClass = "Centrality",
		     Bool_t bRootMoments = kFALSE) {
  //Macro that draws the charge dependent correlation functions PROJECTIONS 
  //for each centrality bin for the different pT of trigger and 
  //associated particles
  TGaxis::SetMaxDigits(3);

  //first we need some libraries
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGCFebye.so");

  //Get the input file
  TString filename = "balanceFunction2D."; 
  if(eventClass == "Centrality"){
    filename += Form("Centrality%.1fTo%.1f",psiMin,psiMax);
    filename += ".PsiAll.PttFrom";
  }
  else if(eventClass == "Multiplicity"){
    filename += Form("Multiplicity%.0fTo%.0f",psiMin,psiMax);
    filename += ".PsiAll.PttFrom";
  }
  else{ // "EventPlane" (default)
    filename += "Centrality";
    filename += gCentrality; filename += ".Psi";
    if((psiMin == -0.5)&&(psiMax == 0.5)) filename += "InPlane.Ptt";
    else if((psiMin == 0.5)&&(psiMax == 1.5)) filename += "Intermediate.Ptt";
    else if((psiMin == 1.5)&&(psiMax == 2.5)) filename += "OutOfPlane.Ptt";
    else if((psiMin == 2.5)&&(psiMax == 3.5)) filename += "Rest.PttFrom";
    else filename += "All.PttFrom";
  }  
  filename += Form("%.1f",ptTriggerMin); filename += "To"; 
  filename += Form("%.1f",ptTriggerMax); filename += "PtaFrom";
  filename += Form("%.1f",ptAssociatedMin); filename += "To"; 
  filename += Form("%.1f",ptAssociatedMax); 
  if(k2pMethod) filename += "_2pMethod";
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
    psiLatex = " -7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o}"; 
  else if((psiMin == 0.5)&&(psiMax == 1.5))
    psiLatex = " 37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o}"; 
  else if((psiMin == 1.5)&&(psiMax == 2.5))
    psiLatex = " 82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o}"; 
  else 
    psiLatex = " 0^{o} < #varphi^{t} - #Psi_{2} < 180^{o}"; 
 
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
  //Get subtracted and mixed balance function

  TH2D *gHistBalanceFunctionSubtracted2D = (TH2D*)f->Get("gHistBalanceFunctionSubtracted");
  TH2D *gHistBalanceFunctionMixed2D      = (TH2D*)f->Get("gHistBalanceFunctionMixed");

  TH1D *gHistBalanceFunctionSubtracted = NULL;
  TH1D *gHistBalanceFunctionMixed      = NULL;

  if(kProjectInEta){
    gHistBalanceFunctionSubtracted = dynamic_cast<TH1D *>(gHistBalanceFunctionSubtracted2D->ProjectionX());
    gHistBalanceFunctionMixed      = dynamic_cast<TH1D *>(gHistBalanceFunctionMixed2D->ProjectionX());
    gHistBalanceFunctionSubtracted->SetTitle("B(#Delta#eta)");
    gHistBalanceFunctionMixed->SetTitle("B_{mix}(#Delta#eta)");  
  }
  else{
    gHistBalanceFunctionSubtracted = dynamic_cast<TH1D *>(gHistBalanceFunctionSubtracted2D->ProjectionY());
    gHistBalanceFunctionMixed      = dynamic_cast<TH1D *>(gHistBalanceFunctionMixed2D->ProjectionY());
    gHistBalanceFunctionSubtracted->SetTitle("B(#Delta#varphi)");
    gHistBalanceFunctionMixed->SetTitle("B_{mix}(#Delta#varphi)");  
  }

  gHistBalanceFunctionSubtracted->SetMarkerStyle(20);
  gHistBalanceFunctionSubtracted->GetYaxis()->SetTitleOffset(1.3);
  gHistBalanceFunctionSubtracted->SetName("gHistBalanceFunctionSubtracted");

  gHistBalanceFunctionMixed->SetMarkerStyle(25);
  gHistBalanceFunctionMixed->SetName("gHistBalanceFunctionMixed");

  TCanvas *c1 = new TCanvas("c1","",0,0,600,500);
  c1->SetFillColor(10); 
  c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.15);
  gHistBalanceFunctionSubtracted->DrawCopy("E");
  gHistBalanceFunctionMixed->DrawCopy("E, SAME");
  
  legend = new TLegend(0.18,0.62,0.45,0.82,"","brNDC");
  legend->SetTextSize(0.045); 
  legend->SetTextFont(42); 
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); 
  legend->SetFillColor(10);
  legend->SetMargin(0.25); 
  legend->SetShadowColor(10);
  legend->AddEntry(gHistBalanceFunctionSubtracted,"Data","lp");
  legend->AddEntry(gHistBalanceFunctionMixed,"Mixed data","lp");
  legend->Draw();
  
  pngName = "BalanceFunction."; 
  if(eventClass == "Centrality"){
    pngName += Form("Centrality%.1fTo%.1f",psiMin,psiMax);
    if(kProjectInEta) pngName += ".InDeltaEta.PsiAll.PttFrom";
    else pngName += ".InDeltaPhi.PsiAll.PttFrom";
  }
  else if(eventClass == "Multiplicity"){
    pngName += Form("Multiplicity%.0fTo%.0f",psiMin,psiMax);
    if(kProjectInEta) pngName += ".InDeltaEta.PsiAll.PttFrom";
    else pngName += ".InDeltaPhi.PsiAll.PttFrom";  
  }
  else{ // "EventPlane" (default)
    pngName += "Centrality";
    pngName += gCentrality; 
    if(kProjectInEta) pngName += ".InDeltaEta.Psi";
    else pngName += ".InDeltaPhi.Psi";
    if((psiMin == -0.5)&&(psiMax == 0.5)) pngName += "InPlane.Ptt";
    else if((psiMin == 0.5)&&(psiMax == 1.5)) pngName += "Intermediate.Ptt";
    else if((psiMin == 1.5)&&(psiMax == 2.5)) pngName += "OutOfPlane.Ptt";
    else if((psiMin == 2.5)&&(psiMax == 3.5)) pngName += "Rest.PttFrom";
    else pngName += "All.PttFrom";
  }  
  pngName += Form("%.1f",ptTriggerMin); pngName += "To"; 
  pngName += Form("%.1f",ptTriggerMax); pngName += "PtaFrom";
  pngName += Form("%.1f",ptAssociatedMin); pngName += "To"; 
  pngName += Form("%.1f",ptAssociatedMax); 
  if(k2pMethod) pngName += "_2pMethod";
  pngName += ".png";

  c1->SaveAs(pngName.Data());
    
  TString meanLatex, rmsLatex, skewnessLatex, kurtosisLatex;

  if(bRootMoments){
    meanLatex = "#mu = "; 
    meanLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetMean());
    meanLatex += " #pm "; 
    meanLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetMeanError());
    
    rmsLatex = "#sigma = "; 
    rmsLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetRMS());
    rmsLatex += " #pm "; 
    rmsLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetRMSError());
    
    skewnessLatex = "S = "; 
    skewnessLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetSkewness(1));
    skewnessLatex += " #pm "; 
    skewnessLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetSkewness(11));
    
    kurtosisLatex = "K = "; 
    kurtosisLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetKurtosis(1));
    kurtosisLatex += " #pm "; 
    kurtosisLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetKurtosis(11));
    Printf("Mean: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetMean(),gHistBalanceFunctionSubtracted->GetMeanError());
    Printf("RMS: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetRMS(),gHistBalanceFunctionSubtracted->GetRMSError());
    Printf("Skeweness: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetSkewness(1),gHistBalanceFunctionSubtracted->GetSkewness(11));
    Printf("Kurtosis: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetKurtosis(1),gHistBalanceFunctionSubtracted->GetKurtosis(11));
  }
  // calculate the moments by hand
  else{

    Double_t meanAnalytical, meanAnalyticalError;
    Double_t sigmaAnalytical, sigmaAnalyticalError;
    Double_t skewnessAnalytical, skewnessAnalyticalError;
    Double_t kurtosisAnalytical, kurtosisAnalyticalError;

    Int_t gDeltaEtaPhi = 2;
    if(kProjectInEta) gDeltaEtaPhi = 1;

    AliBalancePsi *bHelper = new AliBalancePsi;
    bHelper->GetMomentsAnalytical(gDeltaEtaPhi,gHistBalanceFunctionSubtracted,meanAnalytical,meanAnalyticalError,sigmaAnalytical,sigmaAnalyticalError,skewnessAnalytical,skewnessAnalyticalError,kurtosisAnalytical,kurtosisAnalyticalError);

    meanLatex = "#mu = "; 
    meanLatex += Form("%.3f",meanAnalytical);
    meanLatex += " #pm "; 
    meanLatex += Form("%.3f",meanAnalyticalError);
    
    rmsLatex = "#sigma = "; 
    rmsLatex += Form("%.3f",sigmaAnalytical);
    rmsLatex += " #pm "; 
    rmsLatex += Form("%.3f",sigmaAnalyticalError);
    
    skewnessLatex = "S = "; 
    skewnessLatex += Form("%.3f",skewnessAnalytical);
    skewnessLatex += " #pm "; 
    skewnessLatex += Form("%.3f",skewnessAnalyticalError);
    
    kurtosisLatex = "K = "; 
    kurtosisLatex += Form("%.3f",kurtosisAnalytical);
    kurtosisLatex += " #pm "; 
    kurtosisLatex += Form("%.3f",kurtosisAnalyticalError);
    Printf("Mean: %lf - Error: %lf",meanAnalytical, meanAnalyticalError);
    Printf("Sigma: %lf - Error: %lf",sigmaAnalytical, sigmaAnalyticalError);
    Printf("Skeweness: %lf - Error: %lf",skewnessAnalytical, skewnessAnalyticalError);
    Printf("Kurtosis: %lf - Error: %lf",kurtosisAnalytical, kurtosisAnalyticalError);
  }

  TCanvas *c2 = new TCanvas("c2","",600,0,600,500);
  c2->SetFillColor(10); 
  c2->SetHighLightColor(10);
  c2->SetLeftMargin(0.15);
  gHistBalanceFunctionSubtracted->DrawCopy("E");
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.035);
  latex->SetTextColor(1);
  latex->DrawLatex(0.64,0.85,meanLatex.Data());
  latex->DrawLatex(0.64,0.81,rmsLatex.Data());
  latex->DrawLatex(0.64,0.77,skewnessLatex.Data());
  latex->DrawLatex(0.64,0.73,kurtosisLatex.Data());

  TString outFileName = filename;
  if(kProjectInEta) outFileName.ReplaceAll(".root","_DeltaEtaProjection.root");
  else              outFileName.ReplaceAll(".root","_DeltaPhiProjection.root");
  TFile *fProjection = TFile::Open(outFileName.Data(),"recreate");  
  gHistBalanceFunctionSubtracted->Write();
  gHistBalanceFunctionMixed->Write();
  fProjection->Close();
}
