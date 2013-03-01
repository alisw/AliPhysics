const Int_t numberOfCentralityBins = 12;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","0-100","0-1","1-2","2-3"};

void drawBalanceFunctionPsi(const char* filename = "AnalysisResultsPsi.root", 
			    Int_t gCentrality = 1,
			    Int_t gDeltaEtaDeltaPhi = 2,
			    Int_t gBit = -1,
			    const char* gCentralityEstimator = 0x0,
			    Double_t psiMin = -0.5, Double_t psiMax = 0.5,
			    Double_t ptTriggerMin = -1.,
			    Double_t ptTriggerMax = -1.,
			    Double_t ptAssociatedMin = -1.,
			    Double_t ptAssociatedMax = -1.,
			    Bool_t k2pMethod = kFALSE,
			    Bool_t k2pMethod2D = kFALSE,
			    TString eventClass = "EventPlane", //Can be "EventPlane", "Centrality", "Multiplicity"
			    Bool_t bRootMoments = kTRUE)
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

  //correction method check
  if(k2pMethod2D&&!k2pMethod){
    Printf("Chosen 2D 2particle correction method w/o 2particle correction --> not possible");
    return;
  }

  //Prepare the objects and return them
  TList *listBF = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,0);
  TList *listBFShuffled = 0x0;//GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,1);
  TList *listBFMixed = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,2);
  if(!listBF) {
    Printf("The TList object was not created");
    return;
  }
  else 
    draw(listBF,listBFShuffled,listBFMixed,
	 gCentrality,gDeltaEtaDeltaPhi,
	 gCentralityEstimator,
	 psiMin,psiMax,
	 ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,
	 k2pMethod,k2pMethod2D,eventClass,bRootMoments);  
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
	  Int_t gCentrality, Int_t gDeltaEtaDeltaPhi, 
	  const char* gCentralityEstimator,
	  Double_t psiMin, Double_t psiMax,
	  Double_t ptTriggerMin, Double_t ptTriggerMax,
	  Double_t ptAssociatedMin, Double_t ptAssociatedMax,
	  Bool_t k2pMethod = kFALSE,Bool_t k2pMethod2D = kFALSE, TString eventClass="EventPlane",Bool_t bRootMoments=kTRUE) {
  gROOT->LoadMacro("~/SetPlotStyle.C");
  SetPlotStyle();
  gStyle->SetPalette(1,0);
  
  const Int_t gRebin = gDeltaEtaDeltaPhi; //rebin by 2 the Delta phi projection

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

  //balance function shuffling
  AliTHn *hPShuffled = NULL;
  AliTHn *hNShuffled = NULL;
  AliTHn *hPNShuffled = NULL;
  AliTHn *hNPShuffled = NULL;
  AliTHn *hPPShuffled = NULL;
  AliTHn *hNNShuffled = NULL;
  if(listBFShuffled) {
    //listBFShuffled->ls();
    histoName = "fHistP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hPShuffled->SetName("gHistPShuffled");
    histoName = "fHistN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hNShuffled->SetName("gHistNShuffled");
    histoName = "fHistPN";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPNShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hPNShuffled->SetName("gHistPNShuffled");
    histoName = "fHistNP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hNPShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hNPShuffled->SetName("gHistNPShuffled");
    histoName = "fHistPP";
    if(gCentralityEstimator) histoName += gCentralityEstimator;
    hPPShuffled = (AliTHn*) listBFShuffled->FindObject(histoName.Data());
    hPPShuffled->SetName("gHistPPShuffled");
    histoName = "fHistNN";
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
  }

  TH1D *gHistBalanceFunction;
  TH1D *gHistBalanceFunctionShuffled;
  TH1D *gHistBalanceFunctionMixed;
  TH1D *gHistBalanceFunctionSubtracted;
  TString histoTitle, pngName;
  TLegend *legend;
  
  if(eventClass == "Centrality"){
    histoTitle = "Centrality: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " % ";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else if(eventClass == "Multiplicity"){
    histoTitle = "Multiplicity: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " tracks";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else{ // "EventPlane" (default)
    histoTitle = "Centrality: ";
    histoTitle += centralityArray[gCentrality-1]; 
    histoTitle += "%";
    if((psiMin == -0.5)&&(psiMax == 0.5))
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }

  //Raw balance function
  if(k2pMethod && !k2pMethod2D){ 
    if(bMixed){
      gHistBalanceFunction = b->GetBalanceFunctionHistogram2pMethod(0,gDeltaEtaDeltaPhi,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
    }
    else{
      cerr<<"RAW: NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
      return;
    }
  }
  else if(k2pMethod && k2pMethod2D){ 
    if(bMixed){
      if(gDeltaEtaDeltaPhi==1) //Delta eta
	gHistBalanceFunction = b->GetBalanceFunction1DFrom2D2pMethod(0,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
      else //Delta phi
	gHistBalanceFunction = b->GetBalanceFunction1DFrom2D2pMethod(1,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
    }
    else{
      cerr<<"RAW: NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
      return;
    }
  }
  else
    gHistBalanceFunction = b->GetBalanceFunctionHistogram(0,gDeltaEtaDeltaPhi,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistBalanceFunction->SetMarkerStyle(20);
  gHistBalanceFunction->SetTitle(histoTitle.Data());
  gHistBalanceFunction->GetYaxis()->SetTitleOffset(1.3);
  gHistBalanceFunction->SetName("gHistBalanceFunction");
  
  //Shuffled balance function
  //if(k2pMethod){ 
  //if(bMixed)
      //gHistBalanceFunctionShuffled = bShuffled->GetBalanceFunctionHistogram2pMethod(0,gDeltaEtaDeltaPhi,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
  //else{
  //cerr<<"SHUFFLE: NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
  //return;
  //}
  //}
  //else if(k2pMethod2D){ 
  //if(bMixed){
  //  if(gDeltaEtaDeltaPhi==1) //Delta eta
  //gHistBalanceFunctionShuffled = bShuffled->GetBalanceFunction1DFrom2D2pMethod(0,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
  //  else //Delta phi
  //gHistBalanceFunctionShuffled = bShuffled->GetBalanceFunction1DFrom2D2pMethod(1,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
  //}
  //else{
  //  cerr<<"SHUFFLE: NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
  //  return;
  //}
  //}
  //else
  //gHistBalanceFunctionShuffled = bShuffled->GetBalanceFunctionHistogram(0,gDeltaEtaDeltaPhi,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  //gHistBalanceFunctionShuffled->SetMarkerStyle(24);
  //gHistBalanceFunctionShuffled->SetName("gHistBalanceFunctionShuffled");
  
  //Mixed balance function
  if(k2pMethod && !k2pMethod2D){ 
    if(bMixed)
      gHistBalanceFunctionMixed = bMixed->GetBalanceFunctionHistogram2pMethod(0,gDeltaEtaDeltaPhi,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
    else{
      cerr<<"MIXED: NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
      return;
    }
  }
  else if(k2pMethod && k2pMethod2D){ 
    if(bMixed){
      if(gDeltaEtaDeltaPhi==1) //Delta eta
	gHistBalanceFunctionMixed = bMixed->GetBalanceFunction1DFrom2D2pMethod(0,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
      else //Delta phi
	gHistBalanceFunctionMixed = bMixed->GetBalanceFunction1DFrom2D2pMethod(1,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
    }
    else{
      cerr<<"MIXED: NO MIXED BF BUT REQUESTED CORRECTING WITH IT! --> FAIL"<<endl;
      return;
    }
  }
  else
    gHistBalanceFunctionMixed = bMixed->GetBalanceFunctionHistogram(0,gDeltaEtaDeltaPhi,psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  gHistBalanceFunctionMixed->SetMarkerStyle(25);
  gHistBalanceFunctionMixed->SetName("gHistBalanceFunctionMixed");
  
  //Subtracted balance function
  gHistBalanceFunctionSubtracted = dynamic_cast<TH1D *>(gHistBalanceFunction->Clone());
  gHistBalanceFunctionSubtracted->Add(gHistBalanceFunctionMixed,-1);
  gHistBalanceFunctionSubtracted->Rebin(gRebin);
  gHistBalanceFunctionSubtracted->Scale(1./(Double_t)(gRebin));    
  gHistBalanceFunctionSubtracted->SetMarkerStyle(20);
  gHistBalanceFunctionSubtracted->SetTitle(histoTitle.Data());
  gHistBalanceFunctionSubtracted->GetYaxis()->SetTitleOffset(1.3);
  gHistBalanceFunctionSubtracted->SetName("gHistBalanceFunctionSubtracted");

  TCanvas *c1 = new TCanvas("c1","",0,0,600,500);
  c1->SetFillColor(10); 
  c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.15);
  gHistBalanceFunction->DrawCopy("E");
  //gHistBalanceFunctionShuffled->DrawCopy("ESAME");
  gHistBalanceFunctionMixed->DrawCopy("ESAME");
  
  legend = new TLegend(0.18,0.62,0.45,0.82,"","brNDC");
  legend->SetTextSize(0.045); 
  legend->SetTextFont(42); 
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); 
  legend->SetFillColor(10);
  legend->SetMargin(0.25); 
  legend->SetShadowColor(10);
  legend->AddEntry(gHistBalanceFunction,"Data","lp");
  //legend->AddEntry(gHistBalanceFunctionShuffled,"Shuffled data","lp");
  legend->AddEntry(gHistBalanceFunctionMixed,"Mixed data","lp");
  legend->Draw();
  
  pngName = "BalanceFunction."; 
  if(eventClass == "Centrality"){
    pngName += Form("Centrality%.1fTo%.1f",psiMin,psiMax);
    if(gDeltaEtaDeltaPhi == 1) pngName += ".InDeltaEta.PsiAll.PttFrom";
    else if(gDeltaEtaDeltaPhi == 2) pngName += ".InDeltaPhi.PsiAll.PttFrom";
  }
  else if(eventClass == "Multiplicity"){
    pngName += Form("Multiplicity%.0fTo%.0f",psiMin,psiMax);
    if(gDeltaEtaDeltaPhi == 1) pngName += ".InDeltaEta.PsiAll.PttFrom";
    else if(gDeltaEtaDeltaPhi == 2) pngName += ".InDeltaPhi.PsiAll.PttFrom";  
  }
  else{ // "EventPlane" (default)
    pngName += "Centrality";
    pngName += gCentrality; 
    if(gDeltaEtaDeltaPhi == 1) pngName += ".InDeltaEta.Psi";
    else if(gDeltaEtaDeltaPhi == 2) pngName += ".InDeltaPhi.Psi";
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
  if(k2pMethod2D) pngName += "_2pMethod2D";
  else if(k2pMethod) pngName += "_2pMethod";
  pngName += ".png";

  c1->SaveAs(pngName.Data());
  
  GetWeightedMean(gHistBalanceFunction);
  //GetWeightedMean(gHistBalanceFunctionShuffled);
  
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

    b->GetMomentsAnalytical(gHistBalanceFunctionSubtracted,meanAnalytical,meanAnalyticalError,sigmaAnalytical,sigmaAnalyticalError,skewnessAnalytical,skewnessAnalyticalError,kurtosisAnalytical,kurtosisAnalyticalError);

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

  TString newFileName = "balanceFunction."; 
  if(eventClass == "Centrality"){
    newFileName += Form("Centrality%.1fTo%.1f",psiMin,psiMax);
    if(gDeltaEtaDeltaPhi == 1) newFileName += ".InDeltaEta.PsiAll.PttFrom";
    else if(gDeltaEtaDeltaPhi == 2) newFileName += ".InDeltaPhi.PsiAll.PttFrom";
  }
  else if(eventClass == "Multiplicity"){
    newFileName += Form("Multiplicity%.0fTo%.0f",psiMin,psiMax);
    if(gDeltaEtaDeltaPhi == 1) newFileName += ".InDeltaEta.PsiAll.PttFrom";
    else if(gDeltaEtaDeltaPhi == 2) newFileName += ".InDeltaPhi.PsiAll.PttFrom";  
  }
  else{ // "EventPlane" (default)
    newFileName += "Centrality";
    newFileName += gCentrality; 
    if(gDeltaEtaDeltaPhi == 1) newFileName += ".InDeltaEta.Psi";
    else if(gDeltaEtaDeltaPhi == 2) newFileName += ".InDeltaPhi.Psi";
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
  if(k2pMethod2D) newFileName += "_2pMethod2D";
  else if(k2pMethod) newFileName += "_2pMethod";
  newFileName += ".root";

  TFile *fOutput = new TFile(newFileName.Data(),"recreate");
  fOutput->cd();
  gHistBalanceFunction->Write();
  //gHistBalanceFunctionShuffled->Write();
  gHistBalanceFunctionMixed->Write();
  gHistBalanceFunctionSubtracted->Write();
  fOutput->Close();
}

//____________________________________________________________________//
void GetWeightedMean(TH1D *gHistBalance, Int_t fStartBin = 1) {
  //Prints the calculated width of the BF and its error
  Double_t gSumXi = 0.0, gSumBi = 0.0, gSumBiXi = 0.0;
  Double_t gSumBiXi2 = 0.0, gSumBi2Xi2 = 0.0;
  Double_t gSumDeltaBi2 = 0.0, gSumXi2DeltaBi2 = 0.0;
  Double_t deltaBalP2 = 0.0, integral = 0.0;
  Double_t deltaErrorNew = 0.0;

  //Retrieve this variables from Histogram
  Int_t fNumberOfBins = gHistBalance->GetNbinsX();
  Double_t fP2Step    = gHistBalance->GetBinWidth(1); // assume equal binning!
  
  //cout<<"=================================================="<<endl;
  //cout<<"RECALCULATION OF BF WIDTH (StartBin = "<<fStartBin<<")"<<endl;
  //cout<<"HISTOGRAM has "<<fNumberOfBins<<" bins with bin size of "<<fP2Step<<endl;
  //cout<<"=================================================="<<endl;
  for(Int_t i = 1; i <= fNumberOfBins; i++) {
    // this is to simulate |Delta eta| or |Delta phi|
    if(fNumberOfBins/2 - fStartBin + 1 < i && i < fNumberOfBins/2 + fStartBin ) continue;

    //cout<<"B: "<<gHistBalance->GetBinContent(i)<<"\t Error: "<<gHistBalance->GetBinError(i)<<"\t bin: "<<TMath::Abs(gHistBalance->GetBinCenter(i))<<endl;

    gSumXi += TMath::Abs(gHistBalance->GetBinCenter(i)); // this is to simulate |Delta eta| or |Delta phi|
    gSumBi += gHistBalance->GetBinContent(i); 
    gSumBiXi += gHistBalance->GetBinContent(i)*TMath::Abs(gHistBalance->GetBinCenter(i));
    gSumBiXi2 += gHistBalance->GetBinContent(i)*TMath::Power(TMath::Abs(gHistBalance->GetBinCenter(i)),2);
    gSumBi2Xi2 += TMath::Power(gHistBalance->GetBinContent(i),2)*TMath::Power(TMath::Abs(gHistBalance->GetBinCenter(i)),2);
    gSumDeltaBi2 +=  TMath::Power(gHistBalance->GetBinError(i),2);
    gSumXi2DeltaBi2 += TMath::Power(TMath::Abs(gHistBalance->GetBinCenter(i)),2) * TMath::Power(gHistBalance->GetBinError(i),2);
    
    deltaBalP2 += fP2Step*TMath::Power(gHistBalance->GetBinError(i),2);
    integral += fP2Step*gHistBalance->GetBinContent(i);
  }
  for(Int_t i = fStartBin; i < fNumberOfBins; i++)
    deltaErrorNew += gHistBalance->GetBinError(i)*(TMath::Abs(gHistBalance->GetBinCenter(i))*gSumBi - gSumBiXi)/TMath::Power(gSumBi,2);
  
  Double_t integralError = TMath::Sqrt(deltaBalP2);
  
  Double_t delta = gSumBiXi / gSumBi;
  Double_t deltaError = (gSumBiXi / gSumBi) * TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumXi2DeltaBi2)/gSumBiXi),2) + TMath::Power((gSumDeltaBi2/gSumBi),2) );
  cout<<"=================================================="<<endl;
  cout<<"Width: "<<delta<<"\t Error: "<<deltaError<<endl;
  cout<<"New error: "<<deltaErrorNew<<endl;
  cout<<"Integral: "<<integral<<"\t Error: "<<integralError<<endl;
  cout<<"=================================================="<<endl;
  cout<<endl;
}

//______________________________________________________//
void drawBFPsi(const char* lhcPeriod = "LHC10h",
	       const char* gCentralityEstimator = "V0M",
	       Int_t gBit = 128,
	       const char* gEventPlaneEstimator = "VZERO",
	       Bool_t kShowShuffled = kFALSE, 
	       Bool_t kShowMixed = kFALSE, 
	       Int_t gCentrality = 1,
	       Int_t gDeltaEtaDeltaPhi = 2,
	       Double_t psiMin = -0.5, Double_t psiMax = 0.5,
	       Double_t ptTriggerMin = -1.,
	       Double_t ptTriggerMax = -1.,
	       Double_t ptAssociatedMin = -1.,
	       Double_t ptAssociatedMax = -1., 
	       Bool_t k2pMethod = kTRUE) {
  //Macro that draws the BF distributions for each centrality bin
  //for reaction plane dependent analysis
  //Author: Panos.Christakoglou@nikhef.nl
  gROOT->LoadMacro("~/SetPlotStyle.C");
  SetPlotStyle();

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
  filename += "/balanceFunction.Centrality"; 
  filename += gCentrality; filename += ".In";
  if(gDeltaEtaDeltaPhi == 1) filename += "DeltaEta.Psi";
  else if(gDeltaEtaDeltaPhi == 2) filename += "DeltaPhi.Psi";
  if((psiMin == -0.5)&&(psiMax == 0.5)) filename += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) filename += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) filename += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) filename += "Rest.Ptt";
  else filename += "All.PttFrom";
  filename += Form("%.1f",ptTriggerMin); filename += "To"; 
  filename += Form("%.1f",ptTriggerMax); filename += "PtaFrom";
  filename += Form("%.1f",ptAssociatedMin); filename += "To"; 
  filename += Form("%.1f",ptAssociatedMax);    
  if(k2pMethod) filename += "_2pMethod2D";
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
  gHistBalanceFunction->SetMarkerStyle(20);
  gHistBalanceFunction->GetYaxis()->SetTitleOffset(1.3);
  if(gDeltaEtaDeltaPhi == 2) {
    gHistBalanceFunction->GetYaxis()->SetTitle("B(#Delta #varphi)");
    gHistBalanceFunction->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  }

  //Shuffled balance function
  TH1D *gHistBalanceFunctionShuffled = 0x0;
  if(kShowShuffled) {
    gHistBalanceFunctionShuffled = dynamic_cast<TH1D *>(f->Get("gHistBalanceFunctionShuffled"));
    gHistBalanceFunction->SetStats(kFALSE);
    gHistBalanceFunctionShuffled->SetMarkerStyle(24);
  }

  //Mixed balance function
  TH1D *gHistBalanceFunctionMixed = 0x0;
  TH1D *gHistBalanceFunctionSubtracted = 0x0;
  if(kShowMixed) {
    gHistBalanceFunctionMixed = dynamic_cast<TH1D *>(f->Get("gHistBalanceFunctionMixed"));
    gHistBalanceFunction->SetStats(kFALSE);
    gHistBalanceFunctionMixed->SetMarkerStyle(25);

    //Subtracted balance function
    gHistBalanceFunctionSubtracted = dynamic_cast<TH1D *>(f->Get("gHistBalanceFunctionSubtracted"));
    gHistBalanceFunctionSubtracted->SetStats(kFALSE);
    gHistBalanceFunctionSubtracted->SetMarkerStyle(20);
    gHistBalanceFunctionSubtracted->GetYaxis()->SetTitleOffset(1.3);
    if(gDeltaEtaDeltaPhi == 2) {
      gHistBalanceFunctionSubtracted->GetYaxis()->SetTitle("B(#Delta #varphi)");
      gHistBalanceFunctionSubtracted->GetXaxis()->SetTitle("#Delta#varphi (rad)");
    }
  }
  
  TString pngName;
  TLegend *legend;
  
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
    psiLatex = " 0^{o} < #varphi - #Psi_{2} < 180^{o}"; 
 
  TString pttLatex = Form("%.1f",ptTriggerMin);
  pttLatex += " < p_{T}^{t} < "; pttLatex += Form("%.1f",ptTriggerMax);
  pttLatex += " GeV/c";

  TString ptaLatex = Form("%.1f",ptAssociatedMin);
  ptaLatex += " < p_{T}^{a} < "; ptaLatex += Form("%.1f",ptAssociatedMax);
  ptaLatex += " GeV/c";

  //Draw the results
  TCanvas *c1 = new TCanvas("c1","",0,0,600,500);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.17); c1->SetTopMargin(0.05);
  gHistBalanceFunction->SetTitle("");
  gHistBalanceFunction->GetYaxis()->SetTitleOffset(1.4);
  gHistBalanceFunction->GetYaxis()->SetNdivisions(10);
  gHistBalanceFunction->GetXaxis()->SetNdivisions(10);
  gHistBalanceFunction->DrawCopy("E");
  if(kShowShuffled) gHistBalanceFunctionShuffled->DrawCopy("ESAME");
  if(kShowMixed) gHistBalanceFunctionMixed->DrawCopy("ESAME");
  
  legend = new TLegend(0.2,0.72,0.45,0.92,"","brNDC");
  legend->SetTextSize(0.05); 
  legend->SetTextFont(42); 
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); 
  legend->SetFillColor(10);
  legend->SetMargin(0.25); 
  legend->SetShadowColor(10);
  legend->AddEntry(gHistBalanceFunction,"Data","lp");
  if(kShowShuffled) 
    legend->AddEntry(gHistBalanceFunctionShuffled,"Shuffled data","lp");
  if(kShowMixed) 
    legend->AddEntry(gHistBalanceFunctionMixed,"Mixed data","lp");
  legend->Draw();
  
  TLatex *latexInfo1 = new TLatex();
  latexInfo1->SetNDC();
  latexInfo1->SetTextSize(0.04);
  latexInfo1->SetTextColor(1);
  latexInfo1->DrawLatex(0.58,0.88,centralityLatex.Data());
  latexInfo1->DrawLatex(0.58,0.82,psiLatex.Data());
  latexInfo1->DrawLatex(0.58,0.76,pttLatex.Data());
  latexInfo1->DrawLatex(0.58,0.70,ptaLatex.Data());

  //pngName = "BalanceFunctionDeltaEta.Centrality"; 
  //pngName += centralityArray[gCentrality-1]; 
  //pngName += ".Psi"; //pngName += psiMin; pngName += "To"; pngName += psiMax;
  //pngName += ".png";
  //c1->SaveAs(pngName.Data());
    
  TCanvas *c2 = new TCanvas("c2","",600,0,600,500);
  c2->SetFillColor(10); c2->SetHighLightColor(10);
  c2->SetLeftMargin(0.17); c2->SetTopMargin(0.05);
  gHistBalanceFunctionSubtracted->SetTitle("");
  gHistBalanceFunctionSubtracted->GetYaxis()->SetTitleOffset(1.4);
  gHistBalanceFunctionSubtracted->GetYaxis()->SetNdivisions(5);
  gHistBalanceFunctionSubtracted->GetXaxis()->SetNdivisions(10);
  gHistBalanceFunctionSubtracted->DrawCopy("E");
  
  //Opening the output ascii files
  TString filenameMean = "meanIn"; 
  TString filenameSigma = "sigmaIn"; 
  TString filenameSkewness = "skewnessIn"; 
  TString filenameKurtosis = "kurtosisIn"; 
  if(gDeltaEtaDeltaPhi == 1) {
    filenameMean += "DeltaEta.Psi"; filenameSigma += "DeltaEta.Psi";
    filenameSkewness += "DeltaEta.Psi"; filenameKurtosis += "DeltaEta.Psi";
  }
  else if(gDeltaEtaDeltaPhi == 2) {
    filenameMean += "DeltaPhi.Psi"; filenameSigma += "DeltaPhi.Psi";
    filenameSkewness += "DeltaPhi.Psi"; filenameKurtosis += "DeltaPhi.Psi";
  }
  if((psiMin == -0.5)&&(psiMax == 0.5)) {
    filenameMean += "InPlane.Ptt"; filenameSigma += "InPlane.Ptt";
    filenameSkewness += "InPlane.Ptt"; filenameKurtosis += "InPlane.Ptt";
  }
  else if((psiMin == 0.5)&&(psiMax == 1.5)) {
    filenameMean += "Intermediate.Ptt"; filenameSigma += "Intermediate.Ptt";
    filenameSkewness += "Intermediate.Ptt"; 
    filenameKurtosis += "Intermediate.Ptt";
  }
  else if((psiMin == 1.5)&&(psiMax == 2.5)) {
    filenameMean += "OutOfPlane.Ptt"; filenameSigma += "OutOfPlane.Ptt";
    filenameSkewness += "OutOfPlane.Ptt"; 
    filenameKurtosis += "OutOfPlane.Ptt";    
  }
  else if((psiMin == 2.5)&&(psiMax == 3.5)) {
    filenameMean += "Rest.Ptt"; filenameSigma += "Rest.Ptt";
    filenameSkewness += "Rest.Ptt"; filenameKurtosis += "Rest.Ptt";
  }
  else {
    filenameMean += "All.Ptt"; filenameSigma += "All.Ptt";
    filenameSkewness += "All.Ptt"; filenameKurtosis += "All.Ptt";
  }
  filenameMean += Form("%.1f",ptTriggerMin); filenameMean += "To"; 
  filenameMean += Form("%.1f",ptTriggerMax); filenameMean += "PtaFrom";
  filenameMean += Form("%.1f",ptAssociatedMin); filenameMean += "To"; 
  filenameMean += Form("%.1f",ptAssociatedMax);  filenameMean += ".txt";  
  filenameSigma += Form("%.1f",ptTriggerMin); filenameSigma += "To"; 
  filenameSigma += Form("%.1f",ptTriggerMax); filenameSigma += "PtaFrom";
  filenameSigma += Form("%.1f",ptAssociatedMin); filenameSigma += "To"; 
  filenameSigma += Form("%.1f",ptAssociatedMax);  filenameSigma += ".txt";  
  filenameSkewness += Form("%.1f",ptTriggerMin); filenameSkewness += "To"; 
  filenameSkewness += Form("%.1f",ptTriggerMax); filenameSkewness += "PtaFrom";
  filenameSkewness += Form("%.1f",ptAssociatedMin); filenameSkewness += "To"; 
  filenameSkewness += Form("%.1f",ptAssociatedMax);  
  filenameSkewness += ".txt";  
  filenameKurtosis += Form("%.1f",ptTriggerMin); filenameKurtosis += "To"; 
  filenameKurtosis += Form("%.1f",ptTriggerMax); filenameKurtosis += "PtaFrom";
  filenameKurtosis += Form("%.1f",ptAssociatedMin); filenameKurtosis += "To"; 
  filenameKurtosis += Form("%.1f",ptAssociatedMax);  
  filenameKurtosis += ".txt";  
 
  //==================================================================//
  //In Delta phi we calculate the moments from the near side structure//
  if(gDeltaEtaDeltaPhi == 2)
    gHistBalanceFunctionSubtracted->GetXaxis()->SetRangeUser(-TMath::Pi()/2.,TMath::Pi()/2.);
  //==================================================================//

  TString meanLatex, rmsLatex, skewnessLatex, kurtosisLatex;
  meanLatex = "#mu"; 
  if(gDeltaEtaDeltaPhi == 1) meanLatex += "_{#Delta#eta} = "; 
  else if(gDeltaEtaDeltaPhi == 2) meanLatex += "_{#Delta#varphi} = "; 
  meanLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetMean());
  meanLatex += " #pm "; 
  meanLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetMeanError());
  ofstream fileMean(filenameMean.Data(),ios::app);
  fileMean << gCentrality << " " << gHistBalanceFunctionSubtracted->GetMean() << " " <<gHistBalanceFunctionSubtracted->GetMeanError()<<endl;
  fileMean.close();

  rmsLatex = "#sigma"; 
  if(gDeltaEtaDeltaPhi == 1) rmsLatex += "_{#Delta#eta} = "; 
  else if(gDeltaEtaDeltaPhi == 2) rmsLatex += "_{#Delta#varphi} = "; 
  rmsLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetRMS());
  rmsLatex += " #pm "; 
  rmsLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetRMSError());
  ofstream fileSigma(filenameSigma.Data(),ios::app);
  fileSigma << gCentrality << " " << gHistBalanceFunctionSubtracted->GetRMS() << " " <<gHistBalanceFunctionSubtracted->GetRMSError()<<endl;
  fileSigma.close();
  
  skewnessLatex = "S"; 
  if(gDeltaEtaDeltaPhi == 1) skewnessLatex += "_{#Delta#eta} = "; 
  else if(gDeltaEtaDeltaPhi == 2) skewnessLatex += "_{#Delta#varphi} = "; 
  skewnessLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetSkewness(1));
  skewnessLatex += " #pm "; 
  skewnessLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetSkewness(11));
  ofstream fileSkewness(filenameSkewness.Data(),ios::app);
  fileSkewness << gCentrality << " " << gHistBalanceFunctionSubtracted->GetSkewness(1) << " " <<gHistBalanceFunctionSubtracted->GetSkewness(11)<<endl;
  fileSkewness.close();

  kurtosisLatex = "K"; 
  if(gDeltaEtaDeltaPhi == 1) kurtosisLatex += "_{#Delta#eta} = "; 
  else if(gDeltaEtaDeltaPhi == 2) kurtosisLatex += "_{#Delta#varphi} = "; 
  kurtosisLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetKurtosis(1));
  kurtosisLatex += " #pm "; 
  kurtosisLatex += Form("%.3f",gHistBalanceFunctionSubtracted->GetKurtosis(11));
  ofstream fileKurtosis(filenameKurtosis.Data(),ios::app);
  fileKurtosis << gCentrality << " " << gHistBalanceFunctionSubtracted->GetKurtosis(1) << " " <<gHistBalanceFunctionSubtracted->GetKurtosis(11)<<endl;
  fileKurtosis.close();

  Printf("Mean: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetMean(),gHistBalanceFunctionSubtracted->GetMeanError());
  Printf("RMS: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetRMS(),gHistBalanceFunctionSubtracted->GetRMSError());
  Printf("Skeweness: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetSkewness(1),gHistBalanceFunctionSubtracted->GetSkewness(11));
  Printf("Kurtosis: %lf - Error: %lf",gHistBalanceFunctionSubtracted->GetKurtosis(1),gHistBalanceFunctionSubtracted->GetKurtosis(11));

  //latexInfo1->DrawLatex(0.18,0.88,centralityLatex.Data());
  //latexInfo1->DrawLatex(0.18,0.82,psiLatex.Data());
  //latexInfo1->DrawLatex(0.18,0.76,pttLatex.Data());
  //latexInfo1->DrawLatex(0.18,0.70,ptaLatex.Data());
  latexInfo1->DrawLatex(0.59,0.88,centralityLatex.Data());
  latexInfo1->DrawLatex(0.58,0.82,psiLatex.Data());
  latexInfo1->DrawLatex(0.59,0.76,pttLatex.Data());
  latexInfo1->DrawLatex(0.59,0.70,ptaLatex.Data());


  TLatex *latexResults = new TLatex();
  latexResults->SetNDC();
  latexResults->SetTextSize(0.04);
  latexResults->SetTextColor(1);

  if(gDeltaEtaDeltaPhi == 1) {
    latexResults->DrawLatex(0.18,0.88,meanLatex.Data());
    latexResults->DrawLatex(0.18,0.82,rmsLatex.Data());
    latexResults->DrawLatex(0.18,0.76,skewnessLatex.Data());
    latexResults->DrawLatex(0.18,0.70,kurtosisLatex.Data());
  }
  else if(gDeltaEtaDeltaPhi == 2) {
    latexResults->DrawLatex(0.59,0.60,meanLatex.Data());
    latexResults->DrawLatex(0.59,0.54,rmsLatex.Data());
    latexResults->DrawLatex(0.59,0.48,skewnessLatex.Data());
    latexResults->DrawLatex(0.59,0.42,kurtosisLatex.Data());
  }

  TString pngName = "BalanceFunction."; 
  pngName += "Centrality";
  pngName += gCentrality; 
  if(gDeltaEtaDeltaPhi == 1) pngName += ".InDeltaEta.Psi";
  else if(gDeltaEtaDeltaPhi == 2) pngName += ".InDeltaPhi.Psi";
  if((psiMin == -0.5)&&(psiMax == 0.5)) pngName += "InPlane.Ptt";
  else if((psiMin == 0.5)&&(psiMax == 1.5)) pngName += "Intermediate.Ptt";
  else if((psiMin == 1.5)&&(psiMax == 2.5)) pngName += "OutOfPlane.Ptt";
  else if((psiMin == 2.5)&&(psiMax == 3.5)) pngName += "Rest.PttFrom";
  else pngName += "All.PttFrom";  
  pngName += Form("%.1f",ptTriggerMin); pngName += "To"; 
  pngName += Form("%.1f",ptTriggerMax); pngName += "PtaFrom";
  pngName += Form("%.1f",ptAssociatedMin); pngName += "To"; 
  pngName += Form("%.1f",ptAssociatedMax); 
  pngName += ".png";

  c2->SaveAs(pngName.Data());
}
