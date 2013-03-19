const Int_t numberOfCentralityBins = 12;
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","0-100","0-1","1-2","2-3"};

const Int_t gRebin = 1;

void drawCorrelationFunctionPsiChargeIndependent(const char* filename = "AnalysisResults.root", 
						 Int_t gCentrality = 1,
						 Int_t gBit = -1,
						 const char* gCentralityEstimator = 0x0,
						 Bool_t kShowShuffled = kFALSE, 
						 Bool_t kShowMixed = kTRUE, 
						 Double_t psiMin = -0.5, 
						 Double_t psiMax = 3.5,
						 Double_t vertexZMin = -10.,
						 Double_t vertexZMax = 10.,
						 Double_t ptTriggerMin = -1.,
						 Double_t ptTriggerMax = -1.,
						 Double_t ptAssociatedMin = -1.,
						 Double_t ptAssociatedMax = -1.,
						 Bool_t normToTrig = kFALSE,
						 Int_t rebinEta = 1,
						 Int_t rebinPhi = 1) {
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
    draw(list,listShuffled,listMixed,
	 gCentralityEstimator,gCentrality,psiMin,psiMax,vertexZMin,vertexZMax,
	 ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,normToTrig,rebinEta,rebinPhi);
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
    if(gCentralityEstimator) 
      histoName += gCentralityEstimator;   
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
    if(gCentralityEstimator)
      histoName += gCentralityEstimator;
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
    if(gCentralityEstimator)
      histoName += gCentralityEstimator;
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
    if(gCentralityEstimator) 
      histoName += gCentralityEstimator;    
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
    if(gCentralityEstimator)
      histoName += gCentralityEstimator;   
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
    if(gCentralityEstimator) 
      histoName += gCentralityEstimator;
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
	  const char *gCentralityEstimator,
	  Int_t gCentrality, Double_t psiMin, Double_t psiMax,
	  Double_t vertexZMin,
	  Double_t vertexZMax,
	  Double_t ptTriggerMin, Double_t ptTriggerMax,
	  Double_t ptAssociatedMin, Double_t ptAssociatedMax,
	  Bool_t normToTrig, Int_t rebinEta, Int_t rebinPhi) {
  //Draws the correlation functions for every centrality bin
  //(+-), (-+), (++), (--)  
  AliTHn *hP = NULL;
  AliTHn *hN = NULL;
  AliTHn *hPN = NULL;
  AliTHn *hNP = NULL;
  AliTHn *hPP = NULL;
  AliTHn *hNN = NULL;
  
  TString gHistPname = "fHistP"; 
  if(gCentralityEstimator) gHistPname += gCentralityEstimator;
  TString gHistNname = "fHistN";
  if(gCentralityEstimator) gHistNname += gCentralityEstimator;
  TString gHistPNname = "fHistPN"; 
  if(gCentralityEstimator) gHistPNname += gCentralityEstimator;
  TString gHistNPname = "fHistNP";
  if(gCentralityEstimator) gHistNPname += gCentralityEstimator;
  TString gHistPPname = "fHistPP";
  if(gCentralityEstimator) gHistPPname += gCentralityEstimator;
  TString gHistNNname = "fHistNN";
  if(gCentralityEstimator) gHistNNname += gCentralityEstimator;

  hP = (AliTHn*) list->FindObject(gHistPname.Data());
  hN = (AliTHn*) list->FindObject(gHistNname.Data());
  hPN = (AliTHn*) list->FindObject(gHistPNname.Data());
  hNP = (AliTHn*) list->FindObject(gHistNPname.Data());
  hPP = (AliTHn*) list->FindObject(gHistPPname.Data());
  hNN = (AliTHn*) list->FindObject(gHistNNname.Data());

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
    
    gHistPname = "fHistP_shuffle"; 
    if(gCentralityEstimator) gHistPname += gCentralityEstimator;
    gHistNname = "fHistN_shuffle";
    if(gCentralityEstimator) gHistNname += gCentralityEstimator;
    gHistPNname = "fHistPN_shuffle"; 
    if(gCentralityEstimator) gHistPNname += gCentralityEstimator;
    gHistNPname = "fHistNP_shuffle";
    if(gCentralityEstimator) gHistNPname += gCentralityEstimator;
    gHistPPname = "fHistPP_shuffle";
    if(gCentralityEstimator) gHistPPname += gCentralityEstimator;
    gHistNNname = "fHistNN_shuffle";
    if(gCentralityEstimator) gHistNNname += gCentralityEstimator;

    hPShuffled = (AliTHn*) listBFShuffled->FindObject(gHistPname.Data());
    hPShuffled->SetName("gHistPShuffled");
    hNShuffled = (AliTHn*) listBFShuffled->FindObject(gHistNname.Data());
    hNShuffled->SetName("gHistNShuffled");
    hPNShuffled = (AliTHn*) listBFShuffled->FindObject(gHistPNname.Data());
    hPNShuffled->SetName("gHistPNShuffled");
    hNPShuffled = (AliTHn*) listBFShuffled->FindObject(gHistNPname.Data());
    hNPShuffled->SetName("gHistNPShuffled");
    hPPShuffled = (AliTHn*) listBFShuffled->FindObject(gHistPPname.Data());
    hPPShuffled->SetName("gHistPPShuffled");
    hNNShuffled = (AliTHn*) listBFShuffled->FindObject(gHistNNname.Data());
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

    gHistPname = "fHistP"; 
    if(gCentralityEstimator) gHistPname += gCentralityEstimator;
    gHistNname = "fHistN";
    if(gCentralityEstimator) gHistNname += gCentralityEstimator;
    gHistPNname = "fHistPN"; 
    if(gCentralityEstimator) gHistPNname += gCentralityEstimator;
    gHistNPname = "fHistNP";
    if(gCentralityEstimator) gHistNPname += gCentralityEstimator;
    gHistPPname = "fHistPP";
    if(gCentralityEstimator) gHistPPname += gCentralityEstimator;
    gHistNNname = "fHistNN";
    if(gCentralityEstimator) gHistNNname += gCentralityEstimator;
    hPMixed = (AliTHn*) listBFMixed->FindObject(gHistPname.Data());
    hPMixed->SetName("gHistPMixed");
    hNMixed = (AliTHn*) listBFMixed->FindObject(gHistNname.Data());
    hNMixed->SetName("gHistNMixed");
    hPNMixed = (AliTHn*) listBFMixed->FindObject(gHistPNname.Data());
    hPNMixed->SetName("gHistPNMixed");
    hNPMixed = (AliTHn*) listBFMixed->FindObject(gHistNPname.Data());
    hNPMixed->SetName("gHistNPMixed");
    hPPMixed = (AliTHn*) listBFMixed->FindObject(gHistPPname.Data());
    hPPMixed->SetName("gHistPPMixed");
    hNNMixed = (AliTHn*) listBFMixed->FindObject(gHistNNname.Data());
    hNNMixed->SetName("gHistNNMixed");
    
    AliBalancePsi *bMixed = new AliBalancePsi();
    bMixed->SetHistNp(hPMixed);
    bMixed->SetHistNn(hNMixed);
    bMixed->SetHistNpn(hPNMixed);
    bMixed->SetHistNnp(hNPMixed);
    bMixed->SetHistNpp(hPPMixed);
    bMixed->SetHistNnn(hNNMixed);
  }

  TH2D *gHist[6];
  
  TCanvas *c[6];
  TString histoTitle, pngName;
  
  // all charges together
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

  gHist[0] = b->GetCorrelationFunctionChargeIndependent(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  if(rebinEta > 1 || rebinPhi > 1){
    gHist[0]->Rebin2D(rebinEta,rebinPhi);
    gHist[0]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
  }
  gHist[0]->GetYaxis()->SetTitleOffset(1.5);
  gHist[0]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
  gHist[0]->SetTitle(histoTitle.Data());
  c[0] = new TCanvas("c0","",0,0,600,500);
  c[0]->SetFillColor(10); c[0]->SetHighLightColor(10);
  gHist[0]->DrawCopy("surf1fb");
  gPad->SetTheta(30); // default is 30
  //gPad->SetPhi(130); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();
  pngName = "DeltaPhiDeltaEta.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
  pngName += ".PositiveNegative.png";
  //c[0]->SaveAs(pngName.Data());
  
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
    
    gHist[1] = bShuffled->GetCorrelationFunctionChargeIndependent(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHist[1]->Rebin2D(rebinEta,rebinPhi);
      gHist[1]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
    gHist[1]->GetYaxis()->SetTitleOffset(1.5);
    gHist[1]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHist[1]->SetTitle(histoTitle.Data());
    c[1] = new TCanvas("c1","",0,100,600,500);
    c[1]->SetFillColor(10); 
    c[1]->SetHighLightColor(10);
    gHist[1]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "DeltaPhiDeltaEtaShuffled.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositiveNegative.png";
    //c[1]->SaveAs(pngName.Data());
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
    
    // if normalization to trigger then do not divide Event mixing by number of trigger particles
    gHist[2] = bMixed->GetCorrelationFunctionChargeIndependent(psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHist[2]->Rebin2D(rebinEta,rebinPhi);
      gHist[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHist[2]->Integral(gHist[2]->GetXaxis()->FindBin(0-10e-5),gHist[2]->GetXaxis()->FindBin(0+10e-5),1,gHist[2]->GetNbinsX());
      mixedNorm /= gHist[2]->GetNbinsY()*(gHist[2]->GetXaxis()->FindBin(0.01) - gHist[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHist[2]->Scale(1./mixedNorm);
    } 

    gHist[2]->GetYaxis()->SetTitleOffset(1.5);
    gHist[2]->GetYaxis()->SetTitle("#Delta #varphi (rad)");
    gHist[2]->SetTitle(histoTitle.Data());
    c[2] = new TCanvas("c2","",0,200,600,500);
    c[2]->SetFillColor(10); 
    c[2]->SetHighLightColor(10);
    gHist[2]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "DeltaPhiDeltaEtaMixed.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositiveNegative.png";
    //c[2]->SaveAs(pngName.Data());

    //Correlation function (+-)
    gHist[3] = b->GetCorrelationFunction("ALL",psiMin,psiMax,vertexZMin,vertexZMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,bMixed);
    gHist[3]->GetXaxis()->SetRangeUser(-1.5,1.5);
    gHist[3]->GetZaxis()->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
    c[3] = new TCanvas("c3","",0,300,600,500);
    c[3]->SetFillColor(10); 
    c[3]->SetHighLightColor(10);
    gHist[3]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
    pngName = "CorrelationFunction.Centrality"; 
    pngName += centralityArray[gCentrality-1]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositiveNegative.png";
    //c[3]->SaveAs(pngName.Data());

    //Correlation function subtracted awayside (+-)
    gHist[4] = dynamic_cast<TH2D *>(gHist[3]->Clone()); // this will be used to imitate twice the away-side
    gHist[4]->GetXaxis()->SetRangeUser(-1.5,1.5);
    gHist[4]->GetZaxis()->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
    gHist[5] = dynamic_cast<TH2D *>(gHist[3]->Clone()); // this will be the subtracted one
    gHist[5]->GetXaxis()->SetRangeUser(-1.5,1.5);
    gHist[5]->GetZaxis()->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");

    //prepare the double away side histo
    for(Int_t ix = 0; ix < gHist[4]->GetNbinsX(); ix++ ){
      for(Int_t iy = 0; iy < gHist[4]->GetNbinsX(); iy++ ){
	if(iy<gHist[4]->GetNbinsY()/2) gHist[4]->SetBinContent(ix+1,iy+1,gHist[4]->GetBinContent(ix+1,iy+1+gHist[4]->GetNbinsY()/2));
      }
    }
    gHist[5]->Add(gHist[4],-1);

    c[4] = new TCanvas("c3","",0,300,600,500);
    c[4]->SetFillColor(10); 
    c[4]->SetHighLightColor(10);
    gHist[5]->DrawCopy("surf1fb");
    gPad->SetTheta(30); // default is 30
    //gPad->SetPhi(130); // default is 30
    gPad->SetPhi(-60); // default is 30
    gPad->Update();    
  }

  //Write to output file
  TString newFileName = "correlationFunctionChargeIndependent.Centrality";  
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
  gHist[0]->SetName("gHistRaw"); gHist[0]->Write();
  if(listBFShuffled) {
    gHist[1]->SetName("gHistShuffled"); gHist[1]->Write();
  }
  if(listBFMixed) {
    gHist[2]->SetName("gHistMixed"); gHist[2]->Write();

    gHist[3]->SetName("gHistCorrelationFunctions");gHist[3]->Write();
    gHist[4]->SetName("gHistCorrelationFunctionsAwaySide"); gHist[4]->Write();
    gHist[5]->SetName("gHistCorrelationFunctionsSubtracted"); gHist[5]->Write();
  }
  newFile->Close();

  // // some cleaning
  // for(Int_t i = 0; i < 6; i++){

  //   if(!listBFShuffled && i == 1) continue;
  //   if(!listBFMixed && (i == 2 || i == 3 || i == 4 || i == 5)) continue;

  //   if(gHist[i]) delete gHist[i];
    
  //   if(c[i]) delete c[i];
  // }

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
			      Int_t gTrainID = 208,			      
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

  pngName = "CorrelationFunctionChargeIndependent.Centrality"; 
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

  fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			  ptTriggerMin,ptTriggerMax,
			  ptAssociatedMin, ptAssociatedMax,gHistNP);
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

  fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			  ptTriggerMin,ptTriggerMax,
			  ptAssociatedMin, ptAssociatedMax,gHistPP);
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

  fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			  ptTriggerMin,ptTriggerMax,
			  ptAssociatedMin, ptAssociatedMax,gHistNN);
}

// //____________________________________________________________//
// void fitCorrelationFunctions(Int_t gCentrality = 1,
// 			     Double_t psiMin = -0.5, Double_t psiMax = 3.5,
// 			     Double_t ptTriggerMin = -1.,
// 			     Double_t ptTriggerMax = -1.,
// 			     Double_t ptAssociatedMin = -1.,
// 			     Double_t ptAssociatedMax = -1.,
// 			     TH2D *gHist) {

//   cout<<"FITTING FUNCTION (MW style)"<<endl;

//   //near side peak(s): [1]*(TMath::Exp(-TMath::Power((0.5*TMath::Power(((x-[5])/[2]),2)+0.5*TMath::Power((y/[3]),2)),[4])) + 
//   //                        TMath::Exp(-TMath::Power((0.5*TMath::Power(((x+[5])/[6]),2)+0.5*TMath::Power((y/[3]),2)),[4])))
//   //away side peak(s): [7]*(TMath::Exp(-TMath::Power((0.5*TMath::Power(((x-[11])/[8]),2)+0.5*TMath::Power((y/[9]),2)),[10])) + 
//   //                        TMath::Exp(-TMath::Power((0.5*TMath::Power(((x+[11])/[12]),2)+0.5*TMath::Power((y/[9]),2)),[10])))
//   //flow contribution (v1 up to v4): 2.*[13]*([14]*TMath::Cos(y) + [15]*TMath::Cos(2.*y) + [16]*TMath::Cos(3.*y) + [17]*TMath::Cos(4.*y))



//   TF2 *gFitFunction = new TF2("gFitFunction","[0]+[1]*(TMath::Exp(-TMath::Power((0.5*TMath::Power(((x-[5])/[2]),2)+0.5*TMath::Power((y/[3]),2)),[4])) + TMath::Exp(-TMath::Power((0.5*TMath::Power(((x+[5])/[6]),2)+0.5*TMath::Power((y/[3]),2)),[4]))) + [7]*(TMath::Exp(-TMath::Power((0.5*TMath::Power(((x-[11])/[8]),2)+0.5*TMath::Power(((y-TMath::Pi())/[9]),2)),[10])) + TMath::Exp(-TMath::Power((0.5*TMath::Power(((x+[11])/[12]),2)+0.5*TMath::Power(((y-TMath::Pi())/[9]),2)),[10]))) + 2.*[13]*([14]*TMath::Cos(y) + [15]*TMath::Cos(2.*y) + [16]*TMath::Cos(3.*y) + [17]*TMath::Cos(4.*y))",-2.0,2.0,-TMath::Pi()/2.,3.*TMath::Pi()/2.); 
//   gFitFunction->SetName("gFitFunction");


//   //Normalization
//   gFitFunction->SetParName(0,"N1"); 
//   //near side peak(s)
//   gFitFunction->SetParName(1,"N_{near side}");gFitFunction->FixParameter(1,0.0);
//   gFitFunction->SetParName(2,"Sigma_{near side}(delta eta 1)"); gFitFunction->FixParameter(2,0.0);
//   gFitFunction->SetParName(3,"Sigma_{near side}(delta phi)"); gFitFunction->FixParameter(3,0.0);
//   gFitFunction->SetParName(4,"Exponent_{near side}"); gFitFunction->FixParameter(4,1.0);
//   gFitFunction->SetParName(5,"Offset_{near side}"); gFitFunction->FixParameter(5,0.0);
//   gFitFunction->SetParName(6,"Sigma_{near side}(delta eta 2)"); gFitFunction->FixParameter(6,0.0);

//   //away side peak(s)
//   gFitFunction->SetParName(7,"N_{near side}");gFitFunction->FixParameter(7,0.0);
//   gFitFunction->SetParName(8,"Sigma_{near side}(delta eta 1)"); gFitFunction->FixParameter(8,0.0);
//   gFitFunction->SetParName(9,"Sigma_{near side}(delta phi)"); gFitFunction->FixParameter(9,0.0);
//   gFitFunction->SetParName(10,"Exponent_{near side}"); gFitFunction->FixParameter(10,1.0);
//   gFitFunction->SetParName(11,"Offset_{near side}"); gFitFunction->FixParameter(11,0.0);
//   gFitFunction->SetParName(12,"Sigma_{near side}(delta eta 2)"); gFitFunction->FixParameter(12,0.0);

//   //flow contribution
//   gFitFunction->SetParName(13,"N_{flow}"); gFitFunction->SetParameter(13,0.2);
//   gFitFunction->SetParName(14,"V1"); gFitFunction->SetParameter(14,0.005);
//   gFitFunction->SetParName(15,"V2"); gFitFunction->SetParameter(15,0.1);
//   gFitFunction->SetParName(16,"V3"); gFitFunction->SetParameter(16,0.05);
//   gFitFunction->SetParName(17,"V4"); gFitFunction->SetParameter(17,0.005);

//   // flow parameters
//   Double_t fNV = 0.;
//   Double_t fV1 = 0.;
//   Double_t fV2 = 0.;
//   Double_t fV3 = 0.;
//   Double_t fV4 = 0.;
 
//   //Fitting the correlation function (first the edges to extract flow)
//   gHist->Fit("gFitFunction","nm","",1.0,1.6);

//   fNV += gFitFunction->GetParameter(13);
//   fV1 += gFitFunction->GetParameter(14);
//   fV2 += gFitFunction->GetParameter(15);
//   fV3 += gFitFunction->GetParameter(16);
//   fV4 += gFitFunction->GetParameter(17);

//   gHist->Fit("gFitFunction","nm","",-1.6,-1.0);

//   fNV += gFitFunction->GetParameter(13);
//   fV1 += gFitFunction->GetParameter(14);
//   fV2 += gFitFunction->GetParameter(15);
//   fV3 += gFitFunction->GetParameter(16);
//   fV4 += gFitFunction->GetParameter(17);

//   // Now fit the whole with fixed flow
//   gFitFunction->FixParameter(13,fNV/2.);
//   gFitFunction->FixParameter(14,fV1/2.);
//   gFitFunction->FixParameter(15,fV2/2.);
//   gFitFunction->FixParameter(16,fV3/2.);
//   gFitFunction->FixParameter(17,fV4/2.);
  
//   gFitFunction->ReleaseParameter(1);gFitFunction->SetParameter(1,0.3);
//   gFitFunction->ReleaseParameter(2);gFitFunction->SetParameter(2,0.3);gFitFunction->SetParLimits(2,0.05,0.7);
//   gFitFunction->ReleaseParameter(3);gFitFunction->SetParameter(3,0.3);gFitFunction->SetParLimits(3,0.05,1.7);
//   gFitFunction->ReleaseParameter(5);gFitFunction->SetParameter(5,0.7);gFitFunction->SetParLimits(5,0.0,0.9);
//   gFitFunction->ReleaseParameter(6);gFitFunction->SetParameter(6,0.3);gFitFunction->SetParLimits(6,0.01,1.7);

//   gFitFunction->ReleaseParameter(7);gFitFunction->SetParameter(1,0.3);
//   gFitFunction->ReleaseParameter(8);gFitFunction->SetParameter(2,0.3);gFitFunction->SetParLimits(2,0.05,0.7);
//   gFitFunction->ReleaseParameter(9);gFitFunction->SetParameter(3,0.3);gFitFunction->SetParLimits(3,0.05,1.7);
//   gFitFunction->ReleaseParameter(11);gFitFunction->SetParameter(5,0.7);gFitFunction->SetParLimits(5,0.0,0.9);
//   gFitFunction->ReleaseParameter(12);gFitFunction->SetParameter(6,0.3);gFitFunction->SetParLimits(6,0.01,1.7);

//   gHist->Fit("gFitFunction","nm");


//   //Cloning the histogram
//   TString histoName = gHist->GetName(); histoName += "Fit"; 
//   TH2D *gHistFit = new TH2D(histoName.Data(),";#Delta#eta;#Delta#varphi (rad);C(#Delta#eta,#Delta#varphi)",gHist->GetNbinsX(),gHist->GetXaxis()->GetXmin(),gHist->GetXaxis()->GetXmax(),gHist->GetNbinsY(),gHist->GetYaxis()->GetXmin(),gHist->GetYaxis()->GetXmax());
//   TH2D *gHistResidual = dynamic_cast<TH2D *>(gHist->Clone());
//   gHistResidual->SetName("gHistResidual");
//   gHistResidual->Sumw2();

//   for(Int_t iBinDeltaEta = 1; iBinDeltaEta <= gHist->GetNbinsX(); iBinDeltaEta++) {
//     for(Int_t iBinDeltaPhi = 1; iBinDeltaPhi <= gHist->GetNbinsY(); iBinDeltaPhi++) {
//       gHistFit->SetBinContent(iBinDeltaEta,iBinDeltaPhi,gFitFunction->Eval(gHist->GetXaxis()->GetBinCenter(iBinDeltaEta),gHist->GetYaxis()->GetBinCenter(iBinDeltaPhi)));
//     }
//   }
//   gHistResidual->Add(gHistFit,-1);

//   //Write to output file
//   TString newFileName = "correlationFunctionFit";
//   if(histoName.Contains("PN")) newFileName += "PN";
//   else if(histoName.Contains("NP")) newFileName += "NP";
//   else if(histoName.Contains("PP")) newFileName += "PP";
//   else if(histoName.Contains("NN")) newFileName += "NN";
//   newFileName += ".Centrality";  
//   newFileName += gCentrality; newFileName += ".Psi";
//   if((psiMin == -0.5)&&(psiMax == 0.5)) newFileName += "InPlane.Ptt";
//   else if((psiMin == 0.5)&&(psiMax == 1.5)) newFileName += "Intermediate.Ptt";
//   else if((psiMin == 1.5)&&(psiMax == 2.5)) newFileName += "OutOfPlane.Ptt";
//   else if((psiMin == 2.5)&&(psiMax == 3.5)) newFileName += "Rest.PttFrom";
//   else newFileName += "All.PttFrom";
//   newFileName += Form("%.1f",ptTriggerMin); newFileName += "To"; 
//   newFileName += Form("%.1f",ptTriggerMax); newFileName += "PtaFrom";
//   newFileName += Form("%.1f",ptAssociatedMin); newFileName += "To"; 
//   newFileName += Form("%.1f",ptAssociatedMax); 
//   newFileName += ".root";
//   TFile *newFile = TFile::Open(newFileName.Data(),"recreate");
//   gHist->Write();
//   gHistFit->Write();
//   gHistResidual->Write();
//   gFitFunction->Write();
//   newFile->Close();
  

// }

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
  TF2 *gFitFunction = new TF2("gFitFunction","[0]+[1]*TMath::Exp(-TMath::Power((0.5*TMath::Power((x/[2]),2)+0.5*TMath::Power((y/[3]),2)),[4]))+[5]*TMath::Exp(-TMath::Power((0.5*TMath::Power(((y-TMath::Pi())/[6]),2)),[7]))+[8]*TMath::Exp(-TMath::Power((0.5*TMath::Power(((x+[17])/[9]),2)),[10]))+[8]*TMath::Exp(-TMath::Power((0.5*TMath::Power(((x-[17])/[9]),2)),[10]))+[11]*TMath::Power(x,2)+2.*[12]*([13]*TMath::Cos(y) + [14]*TMath::Cos(2.*y) + [15]*TMath::Cos(3.*y) + [16]*TMath::Cos(4.*y))",-2.0,2.0,-TMath::Pi()/2.,3.*TMath::Pi()/2.); 
  gFitFunction->SetName("gFitFunction");


  //Normalization
  gFitFunction->SetParName(0,"N1"); 
  //near side peak
  gFitFunction->SetParName(1,"N_{near side}");gFitFunction->FixParameter(1,0.0);
  gFitFunction->SetParName(2,"Sigma_{near side}(delta eta)"); gFitFunction->FixParameter(2,0.0);
  gFitFunction->SetParName(3,"Sigma_{near side}(delta phi)"); gFitFunction->FixParameter(3,0.0);
  gFitFunction->SetParName(4,"Exponent_{near side}"); gFitFunction->FixParameter(4,1.0);
  //away side ridge
  gFitFunction->SetParName(5,"N_{away side}"); gFitFunction->FixParameter(5,0.0);
  gFitFunction->SetParName(6,"Sigma_{away side}(delta phi)"); gFitFunction->FixParameter(6,0.0);
  gFitFunction->SetParName(7,"Exponent_{away side}"); gFitFunction->FixParameter(7,1.0);
  //longitudinal ridge
  gFitFunction->SetParName(8,"N_{long. ridge}"); gFitFunction->SetParameter(8,0.05);//
  gFitFunction->FixParameter(8,0.0);
  gFitFunction->SetParName(9,"Sigma_{long. ridge}(delta eta)"); gFitFunction->FixParameter(9,0.0);
  gFitFunction->SetParName(10,"Exponent_{long. ridge}"); gFitFunction->FixParameter(10,1.0);
  //wing structures
  gFitFunction->SetParName(11,"N_{wing}"); gFitFunction->FixParameter(11,0.0);
  //flow contribution
  gFitFunction->SetParName(12,"N_{flow}"); gFitFunction->SetParameter(12,0.2);gFitFunction->SetParLimits(12,0.0,10);
  gFitFunction->SetParName(13,"V1"); gFitFunction->SetParameter(13,0.005);gFitFunction->SetParLimits(13,0.0,10);
  gFitFunction->SetParName(14,"V2"); gFitFunction->SetParameter(14,0.1);gFitFunction->SetParLimits(14,0.0,10);
  gFitFunction->SetParName(15,"V3"); gFitFunction->SetParameter(15,0.05);gFitFunction->SetParLimits(15,0.0,10);
  gFitFunction->SetParName(16,"V4"); gFitFunction->SetParameter(16,0.005);gFitFunction->SetParLimits(16,0.0,10);
  gFitFunction->SetParName(17,"Offset"); gFitFunction->FixParameter(17,0.0);

  // flow parameters
  Double_t fNV = 0.;
  Double_t fV1 = 0.;
  Double_t fV2 = 0.;
  Double_t fV3 = 0.;
  Double_t fV4 = 0.;
 
  //Fitting the correlation function (first the edges to extract flow)
  gHist->Fit("gFitFunction","nm","",1.0,1.6);

  fNV += gFitFunction->GetParameter(12);
  fV1 += gFitFunction->GetParameter(13);
  fV2 += gFitFunction->GetParameter(14);
  fV3 += gFitFunction->GetParameter(15);
  fV4 += gFitFunction->GetParameter(16);

  gHist->Fit("gFitFunction","nm","",-1.6,-1.0);

  fNV += gFitFunction->GetParameter(12);
  fV1 += gFitFunction->GetParameter(13);
  fV2 += gFitFunction->GetParameter(14);
  fV3 += gFitFunction->GetParameter(15);
  fV4 += gFitFunction->GetParameter(16);

  // Now fit the whole with fixed flow
  gFitFunction->FixParameter(12,fNV/2.);
  gFitFunction->FixParameter(13,fV1/2.);
  gFitFunction->FixParameter(14,fV2/2.);
  gFitFunction->FixParameter(15,fV3/2.);
  gFitFunction->FixParameter(16,fV4/2.);
  
  gFitFunction->ReleaseParameter(0);gFitFunction->SetParameter(0,1.0);
  gFitFunction->ReleaseParameter(1);gFitFunction->SetParameter(1,0.3);
  gFitFunction->ReleaseParameter(2);gFitFunction->SetParameter(2,0.3);gFitFunction->SetParLimits(2,0.05,0.7);
  gFitFunction->ReleaseParameter(3);gFitFunction->SetParameter(3,0.3);gFitFunction->SetParLimits(3,0.05,1.7);
  //gFitFunction->ReleaseParameter(4);gFitFunction->SetParameter(4,1.0);gFitFunction->SetParLimits(4,0.0,2.0);
  gFitFunction->ReleaseParameter(5);gFitFunction->SetParameter(5,1.0);//gFitFunction->SetParLimits(5,0.0,10);
  gFitFunction->ReleaseParameter(6);gFitFunction->SetParameter(6,0.5);//gFitFunction->SetParLimits(6,0.0,10);
  //gFitFunction->ReleaseParameter(7);gFitFunction->SetParameter(7,1.0);gFitFunction->SetParLimits(7,0.0,2.0);
  gFitFunction->ReleaseParameter(8);gFitFunction->SetParameter(8,0.05);
  gFitFunction->ReleaseParameter(9);gFitFunction->SetParameter(9,0.6);gFitFunction->SetParLimits(9,0.1,10.0);
  //gFitFunction->ReleaseParameter(10);gFitFunction->SetParameter(10,1.0);gFitFunction->SetParLimits(10,0.0,2.0);
  gFitFunction->ReleaseParameter(17);gFitFunction->SetParameter(17,0.7);gFitFunction->SetParLimits(17,0.0,0.9);

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
  gHist->Write();
  gHistFit->Write();
  gHistResidual->Write();
  gFitFunction->Write();
  newFile->Close();
  

}



// //____________________________________________________________//
// void fitCorrelationFunctions(Int_t gCentrality = 1,
// 			     Double_t psiMin = -0.5, Double_t psiMax = 3.5,
// 			     Double_t ptTriggerMin = -1.,
// 			     Double_t ptTriggerMax = -1.,
// 			     Double_t ptAssociatedMin = -1.,
// 			     Double_t ptAssociatedMax = -1.,
// 			     TH2D *gHist) {

//   cout<<"FITTING FUNCTION (HOUSTON)"<<endl;

//   // Fit Function
//    //x axis = delta_eta
//    //y axis = delta_phi
//                                                                                                                                                                                         //  [9]*exp(-1*pow(((x/[10])^2 + (y/[10])^2),0.5)) Hyper expo
//   TF2 *fit1 = new TF2("fit1","[0] + [1]*cos(y) + [2]*cos(2*y) + [3]*cos(3*y) + [4]*cos(4*y) +[5]*cos(5*y)+ [6]*exp(-0.5*pow(((x/[7])^2 + (y/[8])^2),[11])) + [6]*exp(-0.5*pow(((x/[7])^2 + ((y-6.283)/[8])^2),[11]))+ [9]*exp(-1*((x/[10])^2 + (y/[10])^2)) ",-2.0,2.0,-TMath::Pi()/2.,3.*TMath::Pi()/2.);
 
  
   
//   Double_t Parameters[] = {0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.3,0.3,1.0,0.1,0.1};
  	
//   fit1->SetParameters(Parameters);  // input pars from macro arguments
	
//   fit1->SetParName(0,"offset");
//   fit1->SetParName(1,"v1");
//   fit1->SetParName(2,"v2");
//   fit1->SetParName(3,"v3");
//   fit1->SetParName(4,"v4");
//   fit1->SetParName(5,"v5");
//   fit1->SetParName(6,"Ridgeamp");
//   fit1->SetParName(7,"Ridgesigx");
//   fit1->SetParName(8,"Ridgesigy");
//   fit1->SetParName(9,"Expoamp");
//   fit1->SetParName(10,"Exposig");
//   fit1->SetParName(11,"Gausspara");
    

//   //Fit Parameter Ranges
//   fit1->SetParLimits(0,-2.0,2.0);   //offset
//   fit1->SetParLimits(1,-1.0,0.1);   //v1
//   fit1->SetParLimits(2,-1.6,0.9);   //v2
//   fit1->SetParLimits(3,0.0,0.5);    //v3
//   fit1->SetParLimits(4,0.0,0.9);    //v4
//   fit1->SetParLimits(5,0.0,0.9);    //v5
//   fit1->SetParLimits(6,0.0,1.5);    //Ridgeamp
//   fit1->SetParLimits(7,0.1,3.0);    //Ridgesigx
//   fit1->SetParLimits(8,0.1,2.0);    //Ridgesigy
//   fit1->SetParLimits(9,0.0,6.0);      //Expoamp
//   fit1->SetParLimits(10,0.05,0.5); //Exposig
//   fit1->SetParLimits(11,0.0,2.0);    //Gausspara

//   //Fitting the correlation function
//   gHist->Fit("fit1","nm");

//   //Cloning the histogram
//   TString histoName = gHist->GetName(); histoName += "Fit"; 
//   TH2D *gHistFit = new TH2D(histoName.Data(),";#Delta#eta;#Delta#varphi (rad);C(#Delta#eta,#Delta#varphi)",gHist->GetNbinsX(),gHist->GetXaxis()->GetXmin(),gHist->GetXaxis()->GetXmax(),gHist->GetNbinsY(),gHist->GetYaxis()->GetXmin(),gHist->GetYaxis()->GetXmax());
//   TH2D *gHistResidual = dynamic_cast<TH2D *>(gHist->Clone());
//   gHistResidual->SetName("gHistResidual");
//   gHistResidual->Sumw2();

//   for(Int_t iBinDeltaEta = 1; iBinDeltaEta <= gHist->GetNbinsX(); iBinDeltaEta++) {
//     for(Int_t iBinDeltaPhi = 1; iBinDeltaPhi <= gHist->GetNbinsY(); iBinDeltaPhi++) {
//       gHistFit->SetBinContent(iBinDeltaEta,iBinDeltaPhi,fit1->Eval(gHist->GetXaxis()->GetBinCenter(iBinDeltaEta),gHist->GetYaxis()->GetBinCenter(iBinDeltaPhi)));
//     }
//   }
//   gHistResidual->Add(gHistFit,-1);

//   //Write to output file
//   TString newFileName = "correlationFunctionFit";
//   if(histoName.Contains("PN")) newFileName += "PN";
//   else if(histoName.Contains("NP")) newFileName += "NP";
//   else if(histoName.Contains("PP")) newFileName += "PP";
//   else if(histoName.Contains("NN")) newFileName += "NN";
//   newFileName += ".Centrality";  
//   newFileName += gCentrality; newFileName += ".Psi";
//   if((psiMin == -0.5)&&(psiMax == 0.5)) newFileName += "InPlane.Ptt";
//   else if((psiMin == 0.5)&&(psiMax == 1.5)) newFileName += "Intermediate.Ptt";
//   else if((psiMin == 1.5)&&(psiMax == 2.5)) newFileName += "OutOfPlane.Ptt";
//   else if((psiMin == 2.5)&&(psiMax == 3.5)) newFileName += "Rest.PttFrom";
//   else newFileName += "All.PttFrom";
//   newFileName += Form("%.1f",ptTriggerMin); newFileName += "To"; 
//   newFileName += Form("%.1f",ptTriggerMax); newFileName += "PtaFrom";
//   newFileName += Form("%.1f",ptAssociatedMin); newFileName += "To"; 
//   newFileName += Form("%.1f",ptAssociatedMax); 
//   newFileName += ".root";
//   TFile *newFile = TFile::Open(newFileName.Data(),"recreate");
//   gHist->Write();
//   gHistFit->Write();
//   gHistResidual->Write();
//   fit1->Write();
//   newFile->Close();
  

// }
