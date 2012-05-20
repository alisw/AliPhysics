const Int_t numberOfCentralityBins = 9;
TString centralityArray[numberOfCentralityBins] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
Double_t gMinCentrality[numberOfCentralityBins] = {0.,5.,10.,20.,30.,40.,50.,60.,70.};
Double_t gMaxCentrality[numberOfCentralityBins] = {5.,10.,20.,30.,40.,50.,60.,70.,80.};
TString gAnalysisType[7] = {"y","eta","qlong","qout","qside","qinv","phi"};

const Int_t gRebin = 1;
void drawBalanceFunctionPsi(const char* filename = "AnalysisResultsPsi.root", 
			    Double_t psiMin = 0., Double_t psiMax = 7.5) {
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

  //Prepare the objects and return them
  TList *listBF = GetListOfObjects(filename);
  TList *listBFShuffled = GetListOfObjects(filename,kTRUE);
  if(!listBF) {
    Printf("The TList object was not created");
    return;
  }
  else 
    draw(listBF,listBFShuffled,psiMin,psiMax);  
}

//______________________________________________________//
TList *GetListOfObjects(const char* filename, Bool_t kShuffling = kFALSE) {
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
  
  TString listBFName;
  if(!kShuffling) 
    listBFName = "listBF_0-100_V0M_vZ10.0_DCAxy-1.0_DCAz-1.0_Pt0.3-5.0_Eta-0.8-0.8_Chi-1.0_nClus-1_Bit1_withCentralTrigger";
  else if(kShuffling)
    listBFName = "listBFShuffled_0-100_V0M_vZ10.0_DCAxy-1.0_DCAz-1.0_Pt0.3-5.0_Eta-0.8-0.8_Chi-1.0_nClus-1_Bit1_withCentralTrigger";
  listBF = dynamic_cast<TList *>(dir->Get(listBFName.Data()));
  //listBF->ls();

  //Get the histograms
  TString histoName;
  if(!kShuffling)
    histoName = "fHistPV0M";
  else if(kShuffling)
    histoName = "fHistP_shuffleV0M";
  AliTHn *fHistP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));  
  if(!fHistP) {
    Printf("fHistP %s not found!!!",histoName.Data());
    break;
  }
  fHistP->FillParent(); fHistP->DeleteContainers();

  if(!kShuffling)
    histoName = "fHistNV0M";
  else if(kShuffling)
    histoName = "fHistN_shuffleV0M";
  AliTHn *fHistN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistN) {
    Printf("fHistN %s not found!!!",histoName.Data());
    break;
  }
  fHistN->FillParent(); fHistN->DeleteContainers();
    
  if(!kShuffling)
    histoName = "fHistPNV0M";
  else if(kShuffling)
    histoName = "fHistPN_shuffleV0M";
  AliTHn *fHistPN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistPN) {
    Printf("fHistPN %s not found!!!",histoName.Data());
    break;
  }
  fHistPN->FillParent(); fHistPN->DeleteContainers();
  
  if(!kShuffling)
    histoName = "fHistNPV0M";
  else if(kShuffling)
    histoName = "fHistNP_shuffleV0M";
  AliTHn *fHistNP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistNP) {
    Printf("fHistNP %s not found!!!",histoName.Data());
    break;
  }
  fHistNP->FillParent(); fHistNP->DeleteContainers();

  if(!kShuffling)
    histoName = "fHistPPV0M";
  else if(kShuffling)
    histoName = "fHistPP_shuffleV0M";
  AliTHn *fHistPP = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistPP) {
    Printf("fHistPP %s not found!!!",histoName.Data());
    break;
  }
  fHistPP->FillParent(); fHistPP->DeleteContainers();

  if(!kShuffling)
    histoName = "fHistNNV0M";
  else if(kShuffling)
    histoName = "fHistNN_shuffleV0M";
  AliTHn *fHistNN = dynamic_cast<AliTHn *>(listBF->FindObject(histoName.Data()));
  if(!fHistNN) {
    Printf("fHistNN %s not found!!!",histoName.Data());
    break;
  }
  fHistNN->FillParent(); fHistNN->DeleteContainers();
  
  return listBF;
}

//______________________________________________________//
void draw(TList *listBF, TList *listBFShuffled,
	  Double_t psiMin, Double_t psiMax) {
  gROOT->LoadMacro("~/SetPlotStyle.C");
  SetPlotStyle();
  gStyle->SetPalette(1,0);
  
  //balance function
  AliTHn *hP = NULL;
  AliTHn *hN = NULL;
  AliTHn *hPN = NULL;
  AliTHn *hNP = NULL;
  AliTHn *hPP = NULL;
  AliTHn *hNN = NULL;
  //listBF->ls();
  //Printf("=================");
  hP = (AliTHn*) listBF->FindObject("fHistPV0M");
  hN = (AliTHn*) listBF->FindObject("fHistNV0M");
  hPN = (AliTHn*) listBF->FindObject("fHistPNV0M");
  hNP = (AliTHn*) listBF->FindObject("fHistNPV0M");
  hPP = (AliTHn*) listBF->FindObject("fHistPPV0M");
  hNN = (AliTHn*) listBF->FindObject("fHistNNV0M");

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
  //listBFShuffled->ls();
  hPShuffled = (AliTHn*) listBFShuffled->FindObject("fHistP_shuffleV0M");
  hNShuffled = (AliTHn*) listBFShuffled->FindObject("fHistN_shuffleV0M");
  hPNShuffled = (AliTHn*) listBFShuffled->FindObject("fHistPN_shuffleV0M");
  hNPShuffled = (AliTHn*) listBFShuffled->FindObject("fHistNP_shuffleV0M");
  hPPShuffled = (AliTHn*) listBFShuffled->FindObject("fHistPP_shuffleV0M");
  hNNShuffled = (AliTHn*) listBFShuffled->FindObject("fHistNN_shuffleV0M");

  AliBalancePsi *bShuffled = new AliBalancePsi();
  bShuffled->SetHistNp(hPShuffled);
  bShuffled->SetHistNn(hNShuffled);
  bShuffled->SetHistNpn(hPNShuffled);
  bShuffled->SetHistNnp(hNPShuffled);
  bShuffled->SetHistNpp(hPPShuffled);
  bShuffled->SetHistNnn(hNNShuffled);

  TH1D *gHistBalanceFunction[numberOfCentralityBins];
  TH1D *gHistBalanceFunctionShuffled[numberOfCentralityBins];
  TCanvas *c1[numberOfCentralityBins];
  TString histoTitle, pngName;
  TLegend *legend[numberOfCentralityBins];
  
  //loop over the centrality bins
  for(Int_t iCentralityBin = 0; iCentralityBin < numberOfCentralityBins; iCentralityBin++) {  
    histoTitle = "Centrality: "; 
    histoTitle += centralityArray[iCentralityBin]; 
    histoTitle += "%";
    histoTitle += " | "; histoTitle += psiMin; 
    histoTitle += " < #phi - #Psi_{2} < "; histoTitle += psiMax;

    gHistBalanceFunction[iCentralityBin] = b->GetBalanceFunctionHistogram(2,2,gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    //gHistBalanceFunction[iCentralityBin] = b->GetBalanceFunctionHistogram(3,3,gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    gHistBalanceFunction[iCentralityBin]->SetMarkerStyle(20);
    gHistBalanceFunction[iCentralityBin]->SetTitle(histoTitle.Data());
    gHistBalanceFunction[iCentralityBin]->GetYaxis()->SetTitleOffset(1.3);

    gHistBalanceFunctionShuffled[iCentralityBin] = bShuffled->GetBalanceFunctionHistogram(2,2,gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    //gHistBalanceFunctionShuffled[iCentralityBin] = b->GetBalanceFunctionHistogram(3,3,gMinCentrality[iCentralityBin],gMaxCentrality[iCentralityBin],psiMin,psiMax);
    gHistBalanceFunctionShuffled[iCentralityBin]->SetMarkerStyle(24);

    c1[iCentralityBin] = new TCanvas(histoTitle.Data(),"",0,0,600,500);
    c1[iCentralityBin]->SetFillColor(10); 
    c1[iCentralityBin]->SetHighLightColor(10);
    c1[iCentralityBin]->SetLeftMargin(0.15);
    gHistBalanceFunction[iCentralityBin]->Draw("E");
    gHistBalanceFunctionShuffled[iCentralityBin]->Draw("ESAME");

    legend[iCentralityBin] = new TLegend(0.18,0.6,0.45,0.82,"","brNDC");
    legend[iCentralityBin]->SetTextSize(0.045); 
    legend[iCentralityBin]->SetTextFont(42); 
    legend[iCentralityBin]->SetBorderSize(0);
    legend[iCentralityBin]->SetFillStyle(0); 
    legend[iCentralityBin]->SetFillColor(10);
    legend[iCentralityBin]->SetMargin(0.25); 
    legend[iCentralityBin]->SetShadowColor(10);
    legend[iCentralityBin]->AddEntry(gHistBalanceFunction[iCentralityBin],"Data","lp");
    legend[iCentralityBin]->AddEntry(gHistBalanceFunctionShuffled[iCentralityBin],"Shuffled data","lp");
    legend[iCentralityBin]->Draw();

    pngName = "BalanceFunctionDeltaEta.Centrality"; 
    pngName += centralityArray[iCentralityBin]; 
    pngName += ".Psi"; pngName += psiMin; pngName += "To"; pngName += psiMax;
    pngName += ".png";
    c1[iCentralityBin]->SaveAs(pngName.Data());

    GetWeightedMean(gHistBalanceFunction[iCentralityBin]);
    GetWeightedMean(gHistBalanceFunctionShuffled[iCentralityBin]);

    TString meanLatex, rmsLatex, skewnessLatex, kurtosisLatex;
    meanLatex = "#mu = "; 
    meanLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetMean());
    meanLatex += " #pm "; 
    meanLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetMeanError());

    rmsLatex = "#sigma = "; 
    rmsLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetRMS());
    rmsLatex += " #pm "; 
    rmsLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetRMSError());

    skewnessLatex = "S = "; 
    skewnessLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetSkewness(1));
    skewnessLatex += " #pm "; 
    skewnessLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetSkewness(11));
 
    kurtosisLatex = "K = "; 
    kurtosisLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetKurtosis(1));
    kurtosisLatex += " #pm "; 
    kurtosisLatex += Form("%.3f",gHistBalanceFunction[iCentralityBin]->GetKurtosis(11));

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.035);
    latex->SetTextColor(1);
    latex->DrawLatex(0.64,0.85,meanLatex.Data());
    latex->DrawLatex(0.64,0.81,rmsLatex.Data());
    latex->DrawLatex(0.64,0.77,skewnessLatex.Data());
    latex->DrawLatex(0.64,0.73,kurtosisLatex.Data());
    Printf("Mean: %lf - Error: %lf",gHistBalanceFunction[iCentralityBin]->GetMean(),gHistBalanceFunction[iCentralityBin]->GetMeanError());
    Printf("RMS: %lf - Error: %lf",gHistBalanceFunction[iCentralityBin]->GetRMS(),gHistBalanceFunction[iCentralityBin]->GetRMSError());
    Printf("Skeweness: %lf - Error: %lf",gHistBalanceFunction[iCentralityBin]->GetSkewness(1),gHistBalanceFunction[iCentralityBin]->GetSkewness(11));
    Printf("Kurtosis: %lf - Error: %lf",gHistBalanceFunction[iCentralityBin]->GetKurtosis(1),gHistBalanceFunction[iCentralityBin]->GetKurtosis(11));
  }
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
  
  cout<<"=================================================="<<endl;
  cout<<"RECALCULATION OF BF WIDTH (StartBin = "<<fStartBin<<")"<<endl;
  cout<<"HISTOGRAM has "<<fNumberOfBins<<" bins with bin size of "<<fP2Step<<endl;
  cout<<"=================================================="<<endl;
  for(Int_t i = 1; i <= fNumberOfBins; i++) {

    // this is to simulate |Delta eta| or |Delta phi|
    if(fNumberOfBins/2 - fStartBin + 1 < i && i < fNumberOfBins/2 + fStartBin ) continue;

    cout<<"B: "<<gHistBalance->GetBinContent(i)<<"\t Error: "<<gHistBalance->GetBinError(i)<<"\t bin: "<<TMath::Abs(gHistBalance->GetBinCenter(i))<<endl;

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
