const Int_t numberOfCentralityBins = 12;
//TString centralityArray[numberOfCentralityBins] = {"0-100","10-20","20-30","30-40","40-50","50-60","60-70","70-80","0-100","0-1","1-2","2-3"};
TString centralityArray[numberOfCentralityBins] = {"0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","0-100","0-1","1-2","2-3"};

const Int_t gRebin = 1;

void drawCorrelationFunctionPsi(const char* filename = "AnalysisResultsPsi_train_06feb.root", 
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
				Double_t ptAssociatedMax = -1.,
				Bool_t normToTrig = kTRUE,
				Int_t rebinEta = 1,
				Int_t rebinPhi = 1,
				TString eventClass = "EventPlane") //Can be "EventPlane", "Centrality", "Multiplicity"
{ 
  //Macro that draws the correlation functions from the balance function
  //analysis vs the reaction plane
  //Author: Panos.Christakoglou@nikhef.nl
  //gROOT->LoadMacro("~/SetPlotStyle.C");
  //SetPlotStyle();
  gStyle->SetPalette(1,0);

  //Load the PWG2ebye library
  gSystem->Load("libCore.so");        
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGCFebye.so");

  //Prepare the objects and return them
  TList *listQA = NULL;
  TList *list = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,0);
  TList *listShuffled = NULL;
  if(kShowShuffled) listShuffled = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,1);
  TList *listMixed = NULL;
  if(kShowMixed) listMixed = GetListOfObjects(filename,gCentrality,gBit,gCentralityEstimator,2);

  // else  get the QA histograms (for convolution)
  else{    
    //Open the file again
    TFile *f = TFile::Open(filename,"UPDATE");
    if((!f)||(!f->IsOpen())) {
      Printf("The file %s is not found.",filename);
    }
    
    
    TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("PWGCFEbyE.outputBalanceFunctionPsiAnalysis"));
    if(!dir) {   
      Printf("The TDirectoryFile is not found.",filename);
    }

    TString listQAName = "listQAPsi_";
    
    listQAName += centralityArray[gCentrality-1];
    if(gBit > -1) {
      listQAName += "_Bit"; listQAName += gBit; }
    if(gCentralityEstimator) {
      listQAName += "_"; listQAName += gCentralityEstimator;}

    listQA = (TList*)dir->Get(Form("%s",listQAName.Data()));
    if(!listQA) {
      Printf("TList QA not found!!!");
    }
  }

  if(!list) {
    Printf("The TList object was not created");
    return;
  }
  else 
    draw(list,listShuffled,listMixed,listQA,
	 gCentralityEstimator,gCentrality,psiMin,psiMax,
	 ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax,normToTrig,rebinEta,rebinPhi,eventClass);
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
    cout<<"List name (control): "<<listBFName.Data()<<endl;
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
void draw(TList *list, TList *listBFShuffled, TList *listBFMixed, TList *listQA,
	  const char *gCentralityEstimator,
	  Int_t gCentrality, Double_t psiMin, Double_t psiMax,
	  Double_t ptTriggerMin, Double_t ptTriggerMax,
	  Double_t ptAssociatedMin, Double_t ptAssociatedMax,
	  Bool_t normToTrig, Int_t rebinEta, Int_t rebinPhi,TString eventClass) {
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
  hNN->Print();


  //Create the AliBalancePsi object and fill it with the AliTHn objects
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
    bMixed->SetEventClass(eventClass);
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

  // if no mixing then divide by convoluted histograms
  if(!listBFMixed && listQA){

    if(!listQA->FindObject("fHistEtaPhiPos") || !listQA->FindObject("fHistEtaPhiNeg")){

      Printf("fHistEtaPhiPos or fHistEtaPhiNeg not found! --> no convolution");
      listQA = NULL;

    }
    else{ 

      TH2D* fHistPos = (TH2D*)((TH3D*)listQA->FindObject("fHistEtaPhiPos"))->Project3D("xy");
      fHistPos->GetYaxis()->SetRangeUser(-0.79,0.79);
      
      TH2D* fHistNeg = (TH2D*)((TH3D*)listQA->FindObject("fHistEtaPhiNeg"))->Project3D("xy");
      fHistNeg->GetYaxis()->SetRangeUser(-0.79,0.79);
      
      gHistPN[2] = convolute2D(fHistPos, fHistNeg, "hConvPN");
      gHistPN[2]->Scale(1./gHistPN[2]->GetBinContent(gHistPN[2]->FindBin(0,0)));
      
      gHistNP[2] = convolute2D(fHistNeg, fHistPos, "hConvNP");
      gHistNP[2]->Scale(1./gHistNP[2]->GetBinContent(gHistNP[2]->FindBin(0,0)));
      
      gHistNN[2] = convolute2D(fHistNeg, fHistNeg, "hConvNN");
      gHistNN[2]->Scale(1./gHistNN[2]->GetBinContent(gHistNN[2]->FindBin(0,0)));
      
      gHistPP[2] = convolute2D(fHistPos, fHistPos, "hConvPP");
      gHistPP[2]->Scale(1./gHistPP[2]->GetBinContent(gHistPP[2]->FindBin(0,0)));
    }
  }
  
  //(+-)
  if(eventClass == "Centrality"){
    histoTitle = "(+-) | Centrality: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " % ";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else if(eventClass == "Multiplicity"){
    histoTitle = "(+-) | Multiplicity: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " tracks";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else{ // "EventPlane" (default)
    histoTitle = "(+-) | Centrality: ";
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

  gHistPN[0] = b->GetCorrelationFunctionPN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  if(rebinEta > 1 || rebinPhi > 1){
    gHistPN[0]->Rebin2D(rebinEta,rebinPhi);
    gHistPN[0]->Scale(1./(Double_t)(rebinEta*rebinPhi));
  }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
    
    gHistPN[1] = bShuffled->GetCorrelationFunctionPN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistPN[1]->Rebin2D(rebinEta,rebinPhi);
      gHistPN[1]->Scale(1./(Double_t)(rebinEta*rebinPhi));
    }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
    
    // if normalization to trigger then do not divide Event mixing by number of trigger particles
    gHistPN[2] = bMixed->GetCorrelationFunctionPN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistPN[2]->Rebin2D(rebinEta,rebinPhi);
      gHistPN[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));
    }
    
    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistPN[2]->Integral(gHistPN[2]->GetXaxis()->FindBin(0-10e-5),gHistPN[2]->GetXaxis()->FindBin(0+10e-5),1,gHistPN[2]->GetNbinsX());
      mixedNorm /= gHistPN[2]->GetNbinsY()*(gHistPN[2]->GetXaxis()->FindBin(0.01) - gHistPN[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistPN[2]->Scale(1./mixedNorm);
    } 

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
    //gHistPN[3]->GetZaxis()->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
    gHistPN[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    //gHistPN[3]->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
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
  // if no mixing then divide by convoluted histograms
  else if(listQA){
    histoTitle = "(+-) mixed | Centrality: "; 
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
    
    if(rebinEta > 1 || rebinPhi > 1){
      gHistPN[2]->Rebin2D(rebinEta,rebinPhi);
      gHistPN[0]->Scale(1./(Double_t)(rebinEta*rebinPhi));
    }

    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistPN[2]->Integral(gHistPN[2]->GetXaxis()->FindBin(0-10e-5),gHistPN[2]->GetXaxis()->FindBin(0+10e-5),1,gHistPN[2]->GetNbinsX());
      mixedNorm /= gHistPN[2]->GetNbinsY()*(gHistPN[2]->GetXaxis()->FindBin(0.01) - gHistPN[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistPN[2]->Scale(1./mixedNorm);
    } 

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
    //gHistPN[3]->GetZaxis()->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
    gHistPN[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    //gHistPN[3]->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
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
  if(eventClass == "Centrality"){
    histoTitle = "(-+) | Centrality: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " % ";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else if(eventClass == "Multiplicity"){
    histoTitle = "(-+) | Multiplicity: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " tracks";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else{ // "EventPlane" (default)
    histoTitle = "(-+) | Centrality: ";
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

  gHistNP[0] = b->GetCorrelationFunctionNP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  if(rebinEta > 1 || rebinPhi > 1){
    gHistNP[0]->Rebin2D(rebinEta,rebinPhi);
    gHistNP[0]->Scale(1./(Double_t)(rebinEta*rebinPhi));
  }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
    
    gHistNP[1] = bShuffled->GetCorrelationFunctionNP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistNP[1]->Rebin2D(rebinEta,rebinPhi);
      gHistNP[1]->Scale(1./(Double_t)(rebinEta*rebinPhi));
    }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 

    // if normalization to trigger then do not divide Event mixing by number of trigger particles
    gHistNP[2] = bMixed->GetCorrelationFunctionNP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistNP[2]->Rebin2D(rebinEta,rebinPhi);
      gHistNP[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));
    }
    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistNP[2]->Integral(gHistNP[2]->GetXaxis()->FindBin(0-10e-5),gHistNP[2]->GetXaxis()->FindBin(0+10e-5),1,gHistNP[2]->GetNbinsX());
      mixedNorm /= gHistNP[2]->GetNbinsY()*(gHistNP[2]->GetXaxis()->FindBin(0.01) - gHistNP[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistNP[2]->Scale(1./mixedNorm);
    } 

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
    //gHistNP[3]->GetZaxis()->SetTitle("C_{-+}(#Delta#eta,#Delta#varphi)");
    gHistNP[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    //gHistNP[3]->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
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
  // if no mixing then divide by convoluted histograms
  else if(listQA){
    histoTitle = "(-+) mixed | Centrality: "; 
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

    if(rebinEta > 1 || rebinPhi > 1){
      gHistNP[2]->Rebin2D(rebinEta,rebinPhi);
      gHistNP[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));
    }

    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistNP[2]->Integral(gHistNP[2]->GetXaxis()->FindBin(0-10e-5),gHistNP[2]->GetXaxis()->FindBin(0+10e-5),1,gHistNP[2]->GetNbinsX());
      mixedNorm /= gHistNP[2]->GetNbinsY()*(gHistNP[2]->GetXaxis()->FindBin(0.01) - gHistNP[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistNP[2]->Scale(1./mixedNorm);
    } 

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
    //gHistNP[3]->GetZaxis()->SetTitle("C_{-+}(#Delta#eta,#Delta#varphi)");
    gHistNP[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    //gHistNP[3]->SetTitle("C_{-+}(#Delta#eta,#Delta#varphi)");
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
  if(eventClass == "Centrality"){
    histoTitle = "(++) | Centrality: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " % ";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else if(eventClass == "Multiplicity"){
    histoTitle = "(++) | Multiplicity: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " tracks";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else{ // "EventPlane" (default)
    histoTitle = "(++) | Centrality: ";
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

  gHistPP[0] = b->GetCorrelationFunctionPP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  if(rebinEta > 1 || rebinPhi > 1){
    gHistPP[0]->Rebin2D(rebinEta,rebinPhi);
    gHistPP[0]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
  }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
    
    gHistPP[1] = bShuffled->GetCorrelationFunctionPP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistPP[1]->Rebin2D(rebinEta,rebinPhi);
      gHistPP[1]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
    
    // if normalization to trigger then do not divide Event mixing by number of trigger particles
    gHistPP[2] = bMixed->GetCorrelationFunctionPP(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistPP[2]->Rebin2D(rebinEta,rebinPhi);
      gHistPP[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistPP[2]->Integral(gHistPP[2]->GetXaxis()->FindBin(0-10e-5),gHistPP[2]->GetXaxis()->FindBin(0+10e-5),1,gHistPP[2]->GetNbinsX());
      mixedNorm /= gHistPP[2]->GetNbinsY()*(gHistPP[2]->GetXaxis()->FindBin(0.01) - gHistPP[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistPP[2]->Scale(1./mixedNorm);
    } 

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
    //gHistPP[3]->GetZaxis()->SetTitle("C_{++}(#Delta#eta,#Delta#varphi)");
    gHistPP[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    // gHistPP[3]->SetTitle("C_{++}(#Delta#eta,#Delta#varphi)");
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
  // if no mixing then divide by convoluted histograms
  else if(listQA){
        histoTitle = "(++) mixed | Centrality: "; 
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
    
    if(rebinEta > 1 || rebinPhi > 1){
      gHistPP[2]->Rebin2D(rebinEta,rebinPhi);
      gHistPP[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistPP[2]->Integral(gHistPP[2]->GetXaxis()->FindBin(0-10e-5),gHistPP[2]->GetXaxis()->FindBin(0+10e-5),1,gHistPP[2]->GetNbinsX());
      mixedNorm /= gHistPP[2]->GetNbinsY()*(gHistPP[2]->GetXaxis()->FindBin(0.01) - gHistPP[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistPP[2]->Scale(1./mixedNorm);
    } 

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
    //gHistPP[3]->GetZaxis()->SetTitle("C_{++}(#Delta#eta,#Delta#varphi)");
    gHistPP[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    //gHistPP[3]->SetTitle("C_{++}(#Delta#eta,#Delta#varphi)");
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
  if(eventClass == "Centrality"){
    histoTitle = "(--) | Centrality: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " % ";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else if(eventClass == "Multiplicity"){
    histoTitle = "(--) | Multiplicity: ";
    histoTitle += psiMin;
    histoTitle += " - ";
    histoTitle += psiMax;
    histoTitle += " tracks";
    histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
  }
  else{ // "EventPlane" (default)
    histoTitle = "(--) | Centrality: ";
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

  gHistNN[0] = b->GetCorrelationFunctionNN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
  if(rebinEta > 1 || rebinPhi > 1){
    gHistNN[0]->Rebin2D(rebinEta,rebinPhi);
    gHistNN[0]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
  }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
    
    gHistNN[1] = bShuffled->GetCorrelationFunctionNN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistNN[1]->Rebin2D(rebinEta,rebinPhi);
      gHistNN[1]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
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
      histoTitle += " (-7.5^{o} < #varphi^{t} - #Psi_{2} < 7.5^{o})"; 
    else if((psiMin == 0.5)&&(psiMax == 1.5))
      histoTitle += " (37.5^{o} < #varphi^{t} - #Psi_{2} < 52.5^{o})"; 
    else if((psiMin == 1.5)&&(psiMax == 2.5))
      histoTitle += " (82.5^{o} < #varphi^{t} - #Psi_{2} < 97.5^{o})"; 
    else 
      histoTitle += " (0^{o} < #varphi^{t} - #Psi_{2} < 180^{o})"; 
    
    // if normalization to trigger then do not divide Event mixing by number of trigger particles
    gHistNN[2] = bMixed->GetCorrelationFunctionNN(psiMin,psiMax,ptTriggerMin,ptTriggerMax,ptAssociatedMin,ptAssociatedMax);
    if(rebinEta > 1 || rebinPhi > 1){
      gHistNN[2]->Rebin2D(rebinEta,rebinPhi);
      gHistNN[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistNN[2]->Integral(gHistNN[2]->GetXaxis()->FindBin(0-10e-5),gHistNN[2]->GetXaxis()->FindBin(0+10e-5),1,gHistNN[2]->GetNbinsX());
      mixedNorm /= gHistNN[2]->GetNbinsY()*(gHistNN[2]->GetXaxis()->FindBin(0.01) - gHistNN[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistNN[2]->Scale(1./mixedNorm);
    } 

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
    //gHistNN[3]->GetZaxis()->SetTitle("C_{--}(#Delta#eta,#Delta#varphi)");
    gHistNN[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    // gHistNN[3]->SetTitle("C_{--}(#Delta#eta,#Delta#varphi)");
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
  // if no mixing then divide by convoluted histograms
  else if(listQA){
        histoTitle = "(--) mixed | Centrality: "; 
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
    
    if(rebinEta > 1 || rebinPhi > 1){
      gHistNN[2]->Rebin2D(rebinEta,rebinPhi);
      gHistNN[2]->Scale(1./(Double_t)(rebinEta*rebinPhi));  
    }
    // normalization to 1 at (0,0) --> Jan Fietes method
    if(normToTrig){
      Double_t mixedNorm = gHistNN[2]->Integral(gHistNN[2]->GetXaxis()->FindBin(0-10e-5),gHistNN[2]->GetXaxis()->FindBin(0+10e-5),1,gHistNN[2]->GetNbinsX());
      mixedNorm /= gHistNN[2]->GetNbinsY()*(gHistNN[2]->GetXaxis()->FindBin(0.01) - gHistNN[2]->GetXaxis()->FindBin(-0.01) + 1);
      gHistNN[2]->Scale(1./mixedNorm);
    } 

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
    //gHistNN[3]->GetZaxis()->SetTitle("C_{--}(#Delta#eta,#Delta#varphi)");
    gHistNN[3]->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N_{assoc}}{d#Delta#eta#Delta#varphi} (rad^{-1})");
    //gHistNN[3]->SetTitle("C_{--}(#Delta#eta,#Delta#varphi)");
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
  TString newFileName = "correlationFunction."; 
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
  if(listBFMixed || (!listBFMixed&&listQA)) {
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
void drawCorrelationFunctions(const char* lhcPeriod = "LHC10h",
			      const char* gCentralityEstimator = "V0M",
			      Int_t gBit = 128,
			      const char* gEventPlaneEstimator = "VZERO",
			      Int_t gCentrality = 1,
			      Double_t psiMin = -0.5, Double_t psiMax = 3.5,
			      Double_t ptTriggerMin = -1.,
			      Double_t ptTriggerMax = -1.,
			      Double_t ptAssociatedMin = -1.,
			      Double_t ptAssociatedMax = -1.,
			      Bool_t kFit = kFALSE) {
  //Macro that draws the charge dependent correlation functions
  //for each centrality bin for the different pT of trigger and 
  //associated particles
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
    return;
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
  //Get the +- correlation function
  TH2D *gHistPN = dynamic_cast<TH2D *>(f->Get("gHistPNCorrelationFunctions"));
  gHistPN->SetStats(kFALSE);
  gHistPN->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
  gHistPN->GetXaxis()->SetRangeUser(-1.4,1.4);
  gHistPN->GetXaxis()->CenterTitle();
  gHistPN->GetXaxis()->SetTitleOffset(1.2);
  gHistPN->GetYaxis()->CenterTitle();
  gHistPN->GetYaxis()->SetTitleOffset(1.2);
  gHistPN->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cPN = new TCanvas("cPN","",0,0,600,500);
  cPN->SetFillColor(10); cPN->SetHighLightColor(10);
  cPN->SetLeftMargin(0.15);
  gHistPN->DrawCopy("surf1fb");

  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

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
  if(kFit)
    fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			    ptTriggerMin,ptTriggerMax,
			    ptAssociatedMin, ptAssociatedMax,gHistPN);

  //============================================================//
  //Get the -+ correlation function
  TH2D *gHistNP = dynamic_cast<TH2D *>(f->Get("gHistNPCorrelationFunctions"));
  gHistNP->SetStats(kFALSE);
  gHistNP->SetTitle("C_{-+}(#Delta#eta,#Delta#varphi)");
  gHistNP->GetXaxis()->SetRangeUser(-1.4,1);
  gHistNP->GetXaxis()->CenterTitle();
  gHistNP->GetXaxis()->SetTitleOffset(1.2);
  gHistNP->GetYaxis()->CenterTitle();
  gHistNP->GetYaxis()->SetTitleOffset(1.2);
  gHistNP->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cNP = new TCanvas("cNP","",50,50,600,500);
  cNP->SetFillColor(10); cNP->SetHighLightColor(10);
  cNP->SetLeftMargin(0.15);
  gHistNP->DrawCopy("surf1fb");

  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

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

  if(kFit)
    fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			    ptTriggerMin,ptTriggerMax,
			    ptAssociatedMin, ptAssociatedMax,gHistNP);

  //============================================================//
  //Get the ++ correlation function
  TH2D *gHistPP = dynamic_cast<TH2D *>(f->Get("gHistPPCorrelationFunctions"));
  gHistPP->SetStats(kFALSE);
  gHistPP->SetTitle("C_{++}(#Delta#eta,#Delta#varphi)");
  gHistPP->GetXaxis()->SetRangeUser(-1.4,1.4);
  gHistPP->GetXaxis()->CenterTitle();
  gHistPP->GetXaxis()->SetTitleOffset(1.2);
  gHistPP->GetYaxis()->CenterTitle();
  gHistPP->GetYaxis()->SetTitleOffset(1.2);
  gHistPP->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cPP = new TCanvas("cPP","",100,100,600,500);
  cPP->SetFillColor(10); cPP->SetHighLightColor(10);
  cPP->SetLeftMargin(0.15);
  gHistPP->DrawCopy("surf1fb");

  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

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

  if(kFit)
    fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			    ptTriggerMin,ptTriggerMax,
			    ptAssociatedMin, ptAssociatedMax,gHistPP);

  //============================================================//
  //Get the -- correlation function
  TH2D *gHistNN = dynamic_cast<TH2D *>(f->Get("gHistNNCorrelationFunctions"));
  gHistNN->SetStats(kFALSE);
  gHistNN->SetTitle("C_{--}(#Delta#eta,#Delta#varphi)");
  gHistNN->GetXaxis()->SetRangeUser(-1.4,1.4);
  gHistNN->GetXaxis()->CenterTitle();
  gHistNN->GetXaxis()->SetTitleOffset(1.2);
  gHistNN->GetYaxis()->CenterTitle();
  gHistNN->GetYaxis()->SetTitleOffset(1.2);
  gHistNN->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cNN = new TCanvas("cNN","",150,150,600,500);
  cNN->SetFillColor(10); cNN->SetHighLightColor(10);
  cNN->SetLeftMargin(0.15);
  gHistNN->DrawCopy("surf1fb");

  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

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

  if(kFit)
    fitCorrelationFunctions(gCentrality, psiMin, psiMax,
			    ptTriggerMin,ptTriggerMax,
			    ptAssociatedMin, ptAssociatedMax,gHistNN);
}

//____________________________________________________________//
void drawProjections(const char* lhcPeriod = "LHC10h",
		     const char* gCentralityEstimator = "V0M",
		     Int_t gBit = 128,
		     const char* gEventPlaneEstimator = "VZERO",
		     Bool_t kProjectInEta = kFALSE,
		     Int_t gCentrality = 1,
		     Double_t psiMin = -0.5, 
		     Double_t psiMax = 3.5,
		     Double_t ptTriggerMin = -1.,
		     Double_t ptTriggerMax = -1.,
		     Double_t ptAssociatedMin = -1.,
		     Double_t ptAssociatedMax = -1.) {
  //Macro that draws the charge dependent correlation functions PROJECTIONS 
  //for each centrality bin for the different pT of trigger and 
  //associated particles
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
  //Get the +- correlation function
  TH2D *gHistPN = dynamic_cast<TH2D *>(f->Get("gHistPNCorrelationFunctions"));
  gHistPN->SetStats(kFALSE);
  gHistPN->SetTitle("C_{+-}(#Delta#eta,#Delta#varphi)");
  gHistPN->GetXaxis()->SetRangeUser(-1.4,1.4);
  gHistPN->GetXaxis()->CenterTitle();
  gHistPN->GetXaxis()->SetTitleOffset(1.2);
  gHistPN->GetYaxis()->CenterTitle();
  gHistPN->GetYaxis()->SetTitleOffset(1.2);
  gHistPN->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cPN = new TCanvas("cPN","",0,0,600,500);
  cPN->SetFillColor(10); 
  cPN->SetHighLightColor(10);
  cPN->SetLeftMargin(0.15);

  //================//
  TH1D* gHistPNprojection = 0x0;
  Double_t sum = 0.0;
  Double_t gError = 0.0;
  Int_t nCounter = 0;
  if(kProjectInEta) {
    gHistPNprojection = new TH1D("gHistPNprojection","",gHistPN->GetNbinsX(),gHistPN->GetXaxis()->GetXmin(),gHistPN->GetXaxis()->GetXmax());
    for(Int_t iBinX = 1; iBinX <= gHistPN->GetNbinsX(); iBinX++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinY = 1; iBinY <= gHistPN->GetNbinsY(); iBinY++) {
	sum += gHistPN->GetBinContent(iBinX,iBinY);
	if(gHistPN->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistPN->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistPNprojection->SetBinContent(iBinX,sum);
      gHistPNprojection->SetBinError(iBinX,gError);
    }
    gHistPNprojection->GetXaxis()->SetRangeUser(-1.4,1.4);
    //gHistPNprojection = (TH1D*)gHistPN->ProjectionX("gHistPNprojection",1,-1);
    //gHistPNprojection->Scale(1./gHistPN->GetNbinsY());
    gHistPNprojection->Scale(2.*TMath::Pi());
    gHistPNprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#eta} - b_{ZYAM}");
    gHistPNprojection->GetXaxis()->SetTitle("#Delta#eta");
  }
  else {
    gHistPNprojection = new TH1D("gHistPNprojection","",gHistPN->GetNbinsY(),gHistPN->GetYaxis()->GetXmin(),gHistPN->GetYaxis()->GetXmax());
    for(Int_t iBinY = 1; iBinY <= gHistPN->GetNbinsY(); iBinY++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinX = 1; iBinX <= gHistPN->GetNbinsX(); iBinX++) {
	sum += gHistPN->GetBinContent(iBinX,iBinY);
	if(gHistPN->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistPN->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistPNprojection->SetBinContent(iBinY,sum);
      gHistPNprojection->SetBinError(iBinY,gError);
    }
    gHistPNprojection->Scale(3.2);
    gHistPNprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#varphi} - b_{ZYAM} (rad^{-1})");
    gHistPNprojection->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  }
  //ZYAM
  Double_t reference = gHistPNprojection->GetBinContent(gHistPNprojection->GetMinimumBin());
  for(Int_t iBinX = 1; iBinX <= gHistPNprojection->GetNbinsX(); iBinX++) 
    gHistPNprojection->SetBinContent(iBinX,gHistPNprojection->GetBinContent(iBinX) - reference);
  
  gHistPNprojection->GetYaxis()->SetTitleOffset(1.5);
  gHistPNprojection->SetMarkerStyle(20);
  gHistPNprojection->SetStats(kFALSE);
  gHistPNprojection->DrawCopy("E");
  //=================//
  
  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

  pngName = "Projection.CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if(kProjectInEta)
    pngName += ".ETAprojection.";
  else
    pngName += ".PHIprojection.";
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
  
  //============================================================//
  //Get the -+ correlation function
  TH2D *gHistNP = dynamic_cast<TH2D *>(f->Get("gHistNPCorrelationFunctions"));
  gHistNP->SetStats(kFALSE);
  gHistNP->SetTitle("C_{-+}(#Delta#eta,#Delta#varphi)");
  gHistNP->GetXaxis()->SetRangeUser(-1.4,1.4);
  gHistNP->GetXaxis()->CenterTitle();
  gHistNP->GetXaxis()->SetTitleOffset(1.2);
  gHistNP->GetYaxis()->CenterTitle();
  gHistNP->GetYaxis()->SetTitleOffset(1.2);
  gHistNP->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cNP = new TCanvas("cNP","",50,50,600,500);
  cNP->SetFillColor(10); cNP->SetHighLightColor(10);
  cNP->SetLeftMargin(0.15);

  //================//
  TH1D* gHistNPprojection = 0x0;
  Double_t sum = 0.0;
  Double_t gError = 0.0;
  Int_t nCounter = 0;
  if(kProjectInEta) {
    gHistNPprojection = new TH1D("gHistNPprojection","",gHistNP->GetNbinsX(),gHistNP->GetXaxis()->GetXmin(),gHistNP->GetXaxis()->GetXmax());
    for(Int_t iBinX = 1; iBinX <= gHistNP->GetNbinsX(); iBinX++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinY = 1; iBinY <= gHistNP->GetNbinsY(); iBinY++) {
	sum += gHistNP->GetBinContent(iBinX,iBinY);
	if(gHistNP->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistNP->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistNPprojection->SetBinContent(iBinX,sum);
      gHistNPprojection->SetBinError(iBinX,gError);
    }
    gHistNPprojection->GetXaxis()->SetRangeUser(-1.4,1.4);
    //gHistNPprojection = (TH1D*)gHistNP->ProjectionX("gHistNPprojection",1,-1);
    //gHistNPprojection->Scale(1./gHistNP->GetNbinsY());
    gHistNPprojection->Scale(2.*TMath::Pi());
    gHistNPprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#eta} - b_{ZYAM}");
    gHistNPprojection->GetXaxis()->SetTitle("#Delta#eta");
  }
  else {
    gHistNPprojection = new TH1D("gHistNPprojection","",gHistNP->GetNbinsY(),gHistNP->GetYaxis()->GetXmin(),gHistNP->GetYaxis()->GetXmax());
    for(Int_t iBinY = 1; iBinY <= gHistNP->GetNbinsY(); iBinY++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinX = 1; iBinX <= gHistNP->GetNbinsX(); iBinX++) {
	sum += gHistNP->GetBinContent(iBinX,iBinY);
	if(gHistNP->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistNP->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistNPprojection->SetBinContent(iBinY,sum);
      gHistNPprojection->SetBinError(iBinY,gError);
    }
    gHistNPprojection->Scale(3.2);
    gHistNPprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#varphi} - b_{ZYAM} (rad^{-1})");
    gHistNPprojection->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  }
  //ZYAM
  Double_t reference = gHistNPprojection->GetBinContent(gHistNPprojection->GetMinimumBin());
  for(Int_t iBinX = 1; iBinX <= gHistNPprojection->GetNbinsX(); iBinX++) 
    gHistNPprojection->SetBinContent(iBinX,gHistNPprojection->GetBinContent(iBinX) - reference);

  gHistNPprojection->GetYaxis()->SetTitleOffset(1.5);
  gHistNPprojection->SetMarkerStyle(20);
  gHistNPprojection->SetStats(kFALSE);
  gHistNPprojection->DrawCopy("E");
  //================//

  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

  pngName = "Projection.CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if(kProjectInEta)
    pngName += ".ETAprojection.";
  else
    pngName += ".PHIprojection.";
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
  gHistPP->SetTitle("C_{++}(#Delta#eta,#Delta#varphi)");
  gHistPP->GetXaxis()->SetRangeUser(-1.4,1.4);
  gHistPP->GetXaxis()->CenterTitle();
  gHistPP->GetXaxis()->SetTitleOffset(1.2);
  gHistPP->GetYaxis()->CenterTitle();
  gHistPP->GetYaxis()->SetTitleOffset(1.2);
  gHistPP->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cPP = new TCanvas("cPP","",100,100,600,500);
  cPP->SetFillColor(10); cPP->SetHighLightColor(10);
  cPP->SetLeftMargin(0.15);
 
  //=================//
  TH1D* gHistPPprojection = 0x0;
  Double_t sum = 0.0;
  Double_t gError = 0.0;
  Int_t nCounter = 0;
  if(kProjectInEta) {
    gHistPPprojection = new TH1D("gHistPPprojection","",gHistPP->GetNbinsX(),gHistPP->GetXaxis()->GetXmin(),gHistPP->GetXaxis()->GetXmax());
    for(Int_t iBinX = 1; iBinX <= gHistPP->GetNbinsX(); iBinX++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinY = 1; iBinY <= gHistPP->GetNbinsY(); iBinY++) {
	sum += gHistPP->GetBinContent(iBinX,iBinY);
	if(gHistPP->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistPP->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistPPprojection->SetBinContent(iBinX,sum);
      gHistPPprojection->SetBinError(iBinX,gError);
    }
    gHistPPprojection->GetXaxis()->SetRangeUser(-1.4,1.4);
    //gHistPPprojection = (TH1D*)gHistPP->ProjectionX("gHistPPprojection",1,-1);
    //gHistPPprojection->Scale(1./gHistPP->GetNbinsY());
    gHistPPprojection->Scale(2.*TMath::Pi());
    gHistPPprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#eta} - b_{ZYAM}");
    gHistPPprojection->GetXaxis()->SetTitle("#Delta#eta");
  }
  else {
    gHistPPprojection = new TH1D("gHistPPprojection","",gHistPP->GetNbinsY(),gHistPP->GetYaxis()->GetXmin(),gHistPP->GetYaxis()->GetXmax());
    for(Int_t iBinY = 1; iBinY <= gHistPP->GetNbinsY(); iBinY++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinX = 1; iBinX <= gHistPP->GetNbinsX(); iBinX++) {
	sum += gHistPP->GetBinContent(iBinX,iBinY);
	if(gHistPP->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistPP->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistPPprojection->SetBinContent(iBinY,sum);
      gHistPPprojection->SetBinError(iBinY,gError);
    }
    gHistPPprojection->Scale(3.2);
    gHistPPprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#varphi} - b_{ZYAM} (rad^{-1})");
    gHistPPprojection->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  }
  //ZYAM
  Double_t reference = gHistPPprojection->GetBinContent(gHistPPprojection->GetMinimumBin());
  for(Int_t iBinX = 1; iBinX <= gHistPPprojection->GetNbinsX(); iBinX++) 
    gHistPPprojection->SetBinContent(iBinX,gHistPPprojection->GetBinContent(iBinX) - reference);
  
  gHistPPprojection->GetYaxis()->SetTitleOffset(1.5);
  gHistPPprojection->SetMarkerStyle(20);
  gHistPPprojection->SetStats(kFALSE);
  gHistPPprojection->DrawCopy("E");
  //================//

  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

  pngName = "Projection.CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if(kProjectInEta) 
    pngName += ".ETAprojection.";
  else
    pngName += ".PHIprojection.";
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
  gHistNN->SetTitle("C_{--}(#Delta#eta,#Delta#varphi)");
  gHistNN->GetXaxis()->SetRangeUser(-1.4,1.4);
  gHistNN->GetXaxis()->CenterTitle();
  gHistNN->GetXaxis()->SetTitleOffset(1.2);
  gHistNN->GetYaxis()->CenterTitle();
  gHistNN->GetYaxis()->SetTitleOffset(1.2);
  gHistNN->GetZaxis()->SetTitleOffset(1.5);
  TCanvas *cNN = new TCanvas("cNN","",150,150,600,500);
  cNN->SetFillColor(10); cNN->SetHighLightColor(10);
  cNN->SetLeftMargin(0.15);

  //=================//
  TH1D* gHistNNprojection = 0x0;
  Double_t sum = 0.0;
  Double_t gError = 0.0;
  Int_t nCounter = 0;
  if(kProjectInEta) {
    gHistNNprojection = new TH1D("gHistNNprojection","",gHistNN->GetNbinsX(),gHistNN->GetXaxis()->GetXmin(),gHistNN->GetXaxis()->GetXmax());
    for(Int_t iBinX = 1; iBinX <= gHistNN->GetNbinsX(); iBinX++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinY = 1; iBinY <= gHistNN->GetNbinsY(); iBinY++) {
	sum += gHistNN->GetBinContent(iBinX,iBinY);
	if(gHistNN->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistNN->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistNNprojection->SetBinContent(iBinX,sum);
      gHistNNprojection->SetBinError(iBinX,gError);
    }
    gHistNNprojection->GetXaxis()->SetRangeUser(-1.4,1.4);
    //gHistNNprojection = (TH1D*)gHistNN->ProjectionX("gHistNNprojection",1,-1);
    //gHistNNprojection->Scale(1./gHistNN->GetNbinsY());
    gHistNNprojection->Scale(2.*TMath::Pi());
    gHistNNprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#eta} - b_{ZYAM}");
    gHistNNprojection->GetXaxis()->SetTitle("#Delta#eta");
  }
  else {
    gHistNNprojection = new TH1D("gHistNNprojection","",gHistNN->GetNbinsY(),gHistNN->GetYaxis()->GetXmin(),gHistNN->GetYaxis()->GetXmax());
    for(Int_t iBinY = 1; iBinY <= gHistNN->GetNbinsY(); iBinY++) {
      sum = 0.; gError = 0.0; nCounter = 0;
      for(Int_t iBinX = 1; iBinX <= gHistNN->GetNbinsX(); iBinX++) {
	sum += gHistNN->GetBinContent(iBinX,iBinY);
	if(gHistNN->GetBinContent(iBinX,iBinY) != 0.) nCounter += 1;
        Double_t exy = gHistNN->GetCellError(iBinX,iBinY);
	gError += exy*exy;
      }
      if(nCounter != 0) {
	sum /= nCounter;
	gError = TMath::Sqrt(gError)/nCounter;
      }
      gHistNNprojection->SetBinContent(iBinY,sum);
      gHistNNprojection->SetBinError(iBinY,gError);
    }
    gHistNNprojection->Scale(3.2);
    gHistNNprojection->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN_{assoc}}{#Delta#varphi} - b_{ZYAM} (rad^{-1})");
    gHistNNprojection->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  }
  //ZYAM
  Double_t reference = gHistNNprojection->GetBinContent(gHistNNprojection->GetMinimumBin());
  for(Int_t iBinX = 1; iBinX <= gHistNNprojection->GetNbinsX(); iBinX++) 
    gHistNNprojection->SetBinContent(iBinX,gHistNNprojection->GetBinContent(iBinX) - reference); 

  gHistNNprojection->GetYaxis()->SetTitleOffset(1.5);
  gHistNNprojection->SetMarkerStyle(20);
  gHistNNprojection->SetStats(kFALSE);
  gHistNNprojection->DrawCopy("E");
  //=================//

  gPad->SetTheta(30); // default is 30
  gPad->SetPhi(-60); // default is 30
  gPad->Update();

  latexInfo1->DrawLatex(0.6,0.95,centralityLatex.Data());
  latexInfo1->DrawLatex(0.6,0.89,psiLatex.Data());
  latexInfo1->DrawLatex(0.6,0.83,pttLatex.Data());
  latexInfo1->DrawLatex(0.6,0.77,ptaLatex.Data());

  pngName = "Projection.CorrelationFunction.Centrality"; 
  pngName += centralityArray[gCentrality-1]; 
  pngName += ".Psi"; 
  if(kProjectInEta) 
    pngName += ".ETAprojection.";
  else
    pngName += ".PHIprojection.";
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

  TFile *fProjection = TFile::Open("ProjectionCorrelationFunction.root",
					  "recreate");
  
  gHistNPprojection->Write();
  gHistPNprojection->Write();
  gHistNNprojection->Write();
  gHistPPprojection->Write();
  fProjection->Close();
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
  //gFitFunction->ReleaseParameter(8);gFitFunction->SetParameter(8,0.05);
  //gFitFunction->ReleaseParameter(9);gFitFunction->SetParameter(9,0.6);gFitFunction->SetParLimits(9,0.1,10.0);
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

//____________________________________________________________//
TH2D *convolute2D(TH2D* h1, TH2D* h2, TString hname) {
  //
  // this function does the convolution of 2 eta or phi "efficiencies" in a deta or dphi "efficiency"
  // and returns a new histogram which is normalized to the number of bin combinations

  // histogram definition
  TH2D *hConv = NULL;
  hConv = new TH2D(hname.Data(),hname.Data(), h1->GetNbinsY(),-2,2,h1->GetNbinsX(),-TMath::Pi()/2.,3*TMath::Pi()/2.);

  Double_t x1 = 0.;
  Double_t x2 = 0.;
  Double_t y1 = 0.;
  Double_t y2 = 0.;
  Double_t z1 = 0.;
  Double_t z2 = 0.;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;


  // convolution
  for(Int_t i = 0; i < h1->GetNbinsX(); i ++){
    cout<<"CONVOLUTING BIN = "<<i<<" of "<<h1->GetNbinsX()<<endl;
    for(Int_t j = 0; j < h1->GetNbinsY(); j ++){
      for(Int_t k = 0; k < h2->GetNbinsX(); k ++){
	for(Int_t l = 0; l < h2->GetNbinsY(); l ++){

	  x1 = (Double_t)h1->GetXaxis()->GetBinCenter(i+1);
	  y1 = (Double_t)h1->GetYaxis()->GetBinCenter(j+1);
	  x2 = (Double_t)h2->GetXaxis()->GetBinCenter(k+1);
	  y2 = (Double_t)h2->GetYaxis()->GetBinCenter(l+1);
      
	  z1 = (Double_t)h1->GetBinContent(i+1,j+1);
	  z2 = (Double_t)h2->GetBinContent(k+1,l+1);

	  // need the gymnastics to keep the same binning
	  dx = x1 - x2 - (h1->GetXaxis()->GetBinWidth(1)/2.);
	  if(gRandom->Gaus() > 0) dy = y1 - y2 + (h1->GetYaxis()->GetBinWidth(1)/2.);
	  else                    dy = y1 - y2 - (h1->GetYaxis()->GetBinWidth(1)/2.);

	  if(dx>3./2.*TMath::Pi())  dx = dx - 2.*TMath::Pi();  
	  if(dx<-1./2.*TMath::Pi()) dx = 2*TMath::Pi() + dx;  

	  dz = z1*z2;

	  hConv->Fill(dy,dx,dz);

	}
      }
    }
  }

  return hConv;

}  
