const Int_t numberOfCentralityBins = 9;
TString centralityArray[numberOfCentralityBins] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
Double_t gMinCentrality[numberOfCentralityBins] = {0.,5.,10.,20.,30.,40.,50.,60.,70.};
Double_t gMaxCentrality[numberOfCentralityBins] = {5.,10.,20.,30.,40.,50.,60.,70.,80.};
TString gAnalysisType[7] = {"y","eta","qlong","qout","qside","qinv","phi"};

const Int_t gRebin = 1;
void drawCorrelationFunctionPsi(const char* filename = "AnalysisResults.root", 
				Double_t psiMin = 0., Double_t psiMax = 7.5) {
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
    draw(list,psiMin,psiMax);
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
  
  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("PWGCFEbyE.outputBalanceFunctionAnalysis"));
  if(!dir) {   
    Printf("The TDirectoryFile is not found. Aborting...",filename);
    return listBF;
  }
  //dir->ls();
  
  TString listBFName = "listBF_0-100_V0M_vZ10.0_DCAxy-1.0_DCAz-1.0_Pt0.3-5.0_Eta-0.8-0.8_Chi-1.0_nClus-1_Bit1_withCentralTrigger";
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
void draw(TList *list, Double_t psiMin, Double_t psiMax) {
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
  TH2D *gHistPN[numberOfCentralityBins];
  TH2D *gHistNP[numberOfCentralityBins];
  TH2D *gHistPP[numberOfCentralityBins];
  TH2D *gHistNN[numberOfCentralityBins];
  
  TCanvas *cPN[numberOfCentralityBins];
  TCanvas *cNP[numberOfCentralityBins];
  TCanvas *cPP[numberOfCentralityBins];
  TCanvas *cNN[numberOfCentralityBins];
  TString histoTitle, pngName;
  //loop over the centrality bins
  for(Int_t iCentralityBin = 0; iCentralityBin < numberOfCentralityBins; iCentralityBin++) {
    histoTitle = "(+-) | Centrality: "; 
    histoTitle += centralityArray[iCentralityBin]; 
    histoTitle += "%";
    histoTitle += " | "; histoTitle += psiMin; 
    histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
    gHistPN[iCentralityBin] = b->GetCorrelationFunctionPN(gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    gHistPN[iCentralityBin]->GetYaxis()->SetTitleOffset(1.5);
    gHistPN[iCentralityBin]->SetTitle(histoTitle.Data());
    cPN[iCentralityBin] = new TCanvas(histoTitle.Data(),"",0,0,400,400);
    cPN[iCentralityBin]->SetFillColor(10); 
    cPN[iCentralityBin]->SetHighLightColor(10);
    gHistPN[iCentralityBin]->DrawCopy("lego2");
    pngName = "DeltaPhiDeltaEta.Centrality"; 
    pngName += centralityArray[iCentralityBin]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositiveNegative.png";
    cPN[iCentralityBin]->SaveAs(pngName.Data());

    histoTitle = "(-+) | Centrality: "; 
    histoTitle += centralityArray[iCentralityBin]; 
    histoTitle += "%";
    histoTitle += " | "; histoTitle += psiMin; 
    histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
   gHistNP[iCentralityBin] = b->GetCorrelationFunctionNP(gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    gHistNP[iCentralityBin]->GetYaxis()->SetTitleOffset(1.5);
    gHistNP[iCentralityBin]->SetTitle(histoTitle.Data());
    cNP[iCentralityBin] = new TCanvas(histoTitle.Data(),"",400,0,400,400);
    cNP[iCentralityBin]->SetFillColor(10); 
    cNP[iCentralityBin]->SetHighLightColor(10);
    gHistNP[iCentralityBin]->DrawCopy("lego2");
    pngName = "DeltaPhiDeltaEta.Centrality"; 
    pngName += centralityArray[iCentralityBin]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativePositive.png";
    cNP[iCentralityBin]->SaveAs(pngName.Data());

    histoTitle = "(++) | Centrality: "; 
    histoTitle += centralityArray[iCentralityBin]; 
    histoTitle += "%";
    histoTitle += " | "; histoTitle += psiMin; 
    histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
    gHistPP[iCentralityBin] = b->GetCorrelationFunctionPP(gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    gHistPP[iCentralityBin]->GetYaxis()->SetTitleOffset(1.5);
    gHistPP[iCentralityBin]->SetTitle(histoTitle.Data());
    cPP[iCentralityBin] = new TCanvas(histoTitle.Data(),"",0,400,400,400);
    cPP[iCentralityBin]->SetFillColor(10); 
    cPP[iCentralityBin]->SetHighLightColor(10);
    gHistPP[iCentralityBin]->DrawCopy("lego2");
    pngName = "DeltaPhiDeltaEta.Centrality"; 
    pngName += centralityArray[iCentralityBin]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".PositivePositive.png";
    cPP[iCentralityBin]->SaveAs(pngName.Data());

    histoTitle = "(--) | Centrality: "; 
    histoTitle += centralityArray[iCentralityBin]; 
    histoTitle += "%";
    histoTitle += " | "; histoTitle += psiMin; 
    histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;
    gHistNN[iCentralityBin] = b->GetCorrelationFunctionNN(gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    gHistNN[iCentralityBin]->GetYaxis()->SetTitleOffset(1.5);
    gHistNN[iCentralityBin]->SetTitle(histoTitle.Data());
    cNN[iCentralityBin] = new TCanvas(histoTitle.Data(),"",400,400,400,400);
    cNN[iCentralityBin]->SetFillColor(10); 
    cNN[iCentralityBin]->SetHighLightColor(10);
    gHistNN[iCentralityBin]->DrawCopy("lego2");
    pngName = "DeltaPhiDeltaEta.Centrality"; 
    pngName += centralityArray[iCentralityBin]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".NegativeNegative.png";
    cNN[iCentralityBin]->SaveAs(pngName.Data());
  }//end of loop over centralities
}

