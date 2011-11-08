const TString gBFAnalysisType[7] = {"y","eta","qlong","qout","qside","qinv","phi"};

const Double_t cent[9]  = {382.8,329.7,260.5,186.4,128.9,85.,52.8,30.,15.8};   // hard coded at the moment for centrality percentiles 
const Double_t centE[9] = {3.1,4.6,4.4,3.9,3.3,2.6,2.0,1.3,0.6};               // (0-5,5-10,10-20,20-30,...,70-80)

void readBalanceFunction(Bool_t bHistos = kTRUE, TString inFile = "AnalysisResults.root",Int_t fStartBinBFWidth = 1, Int_t fRebin = 1) {
  // Macro to read the output of the BF analysis:  MW: CHANGE THIS!!!!
  //i) Prints and draws the final BF output
  //ii) Plots the QA part of the analysis
  //iii) store BF in output file
  //Author: Panos.Christakoglou@cern.ch, m.weber@cern.ch
  //Loading the needed libraries
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG2ebye.so");

  //Draw BF              
  drawBF(bHistos,inFile, fStartBinBFWidth, fRebin);

  //Merge the output
  //mergeOutput("/alice/cern.ch/user/p/pchrist/Balance/pp/7TeV/LHC10b/output/");
}

//___________________________________________________________//
void drawBF(Bool_t bHistos = kTRUE, TString inFile = "AnalysisResults.root", Int_t fStartBinBFWidth = 1, Int_t fRebin = 1) {
  //Function to draw the BF objects and write them into the output file

  Int_t maximumCanvases = 10;
  Int_t iCanvas         = 0;
  TCanvas *cQA[10];
  TCanvas *cQAV0M = new TCanvas("cQAV0M","V0M multiplicities");
  cQAV0M->Divide(4,3);
  TCanvas *cQARef = new TCanvas("cQARef","reference track multiplicities");
  cQARef->Divide(4,3);
  TCanvas *cBF[10];
  TCanvas *cBFS[10];

  // get the file
  TFile *f = TFile::Open(inFile.Data());
  if(!f) {
    Printf("File not found!!!");
    break;
  }

  // get the BF output directory
  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("PWG2EbyE.outputBalanceFunctionAnalysis"));
  if(!dir) {
    Printf("Output directory not found!!!");
    break;
  }

  // loop over all lists and plot the BF and QA
  TList *list = NULL;
  TString listName;
  TIter nextkey( dir->GetListOfKeys() );
  TKey *key;

  AliBalance *bf[10][7];
  AliBalance *bfs[10][7];
  TH1D *gbf[10][7];
  TH1D *gbfs[10][7];

  for(Int_t i = 0; i < 10; i++){
    for(Int_t j = 0; j < 7; j++){
      gbf[i][j] = NULL;
      gbfs[i][j] = NULL;
    }
  }

  TH1D *fHistP[7]; //N+
  TH1D *fHistN[7]; //N-
  TH1D *fHistPN[7]; //N+-
  TH1D *fHistNP[7]; //N-+
  TH1D *fHistPP[7]; //N++
  TH1D *fHistNN[7]; //N--

  TH1D *fHistPS[7]; //N+
  TH1D *fHistNS[7]; //N-
  TH1D *fHistPNS[7]; //N+-
  TH1D *fHistNPS[7]; //N-+
  TH1D *fHistPPS[7]; //N++
  TH1D *fHistNNS[7]; //N--

  Double_t WM[10];     // weighted mean for eta (recalculated from fStartBin)
  Double_t WME[10];    // error
  Double_t WMS[10];     // weighted mean for eta (recalculated from fStartBin) (shuffled)
  Double_t WMSE[10];    // error (shuffled)

  while ( (key = (TKey*)nextkey())) {

    list = (TList*)key->ReadObj();
    listName = TString(list->GetName());

    cout<<"Processing list "<<listName<<endl;

    // ----------------------------------------------------
    // plot QA histograms
    if(listName.Contains("QA")){     

      iCanvas ++;
      if(iCanvas == 10) {cout<<"TOO MANY LISTS --> increase MAXIMUM"<<endl; return;}
      cQA[iCanvas] = new TCanvas(listName,listName);
      cQA[iCanvas]->Divide(4,3);

     cQA[iCanvas]->cd(1);
      TH1F* histEventStats = (TH1F*)list->FindObject("fHistEventStats");
      if(histEventStats){
	histEventStats->SetFillColor(9);
	histEventStats->Draw();
      }

      cQA[iCanvas]->cd(2);
      TH1F* histTriggerStats = (TH1F*)list->FindObject("fHistTriggerStats");
      if(histTriggerStats){
	histTriggerStats->SetFillColor(9);
	histTriggerStats->Draw();
      }

     cQA[iCanvas]->cd(3);
      TH1F* histTrackStats = (TH1F*)list->FindObject("fHistTrackStats");
      if(histTrackStats){
	histTrackStats->SetFillColor(9);
	histTrackStats->Draw();
      }

      cQA[iCanvas]->cd(4);
      TH1F* histCentStats = (TH1F*)list->FindObject("fHistCentStats");
      if(histCentStats){
	histCentStats->SetFillColor(9);
	histCentStats->Draw("colz");
      }



      cQA[iCanvas]->cd(5);
      TH1F* histVx = (TH1F*)list->FindObject("fHistVx");
      if(histVx){
	histVx->SetFillColor(9);
	histVx->Draw();
      }
      cQA[iCanvas]->cd(6);
      TH1F* histVy = (TH1F*)list->FindObject("fHistVy");
      if(histVy){
	histVy->SetFillColor(9);
	histVy->Draw();
      }
      
      cQA[iCanvas]->cd(7);
      TH1F* histVz = (TH1F*)list->FindObject("fHistVz");
      if(histVz){
	histVz->SetFillColor(9);
	histVz->Draw();
      }

      cQA[iCanvas]->cd(8);
      cQA[iCanvas]->cd(8)->SetLogz();
      TH2F* histDCA = (TH2F*)list->FindObject("fHistDCA");  
      if(histDCA) histDCA->Draw("colz");
      

      cQA[iCanvas]->cd(9);
      cQA[iCanvas]->cd(9)->SetLogz();
      TH2F* histClus = (TH2F*)list->FindObject("fHistClus");
      if(histClus) histClus->Draw("colz");
      
      cQA[iCanvas]->cd(10);
      TH1F* histPt = (TH1F*)list->FindObject("fHistPt");
      if(histPt){
	histPt->SetFillColor(9);
	histPt->Draw();
      }
      
      cQA[iCanvas]->cd(11);
      TH1F* histEta = (TH1F*)list->FindObject("fHistEta");
      if(histEta){
	histEta->SetFillColor(9);
	histEta->Draw();
      }
      
      cQA[iCanvas]->cd(12);
      TH1F* histPhi = (TH1F*)list->FindObject("fHistPhi");
      if(histPhi){
	histPhi->SetFillColor(9);
	histPhi->Draw();
      }

      // centrality estimator QA
      cQAV0M->cd(iCanvas);
      cQAV0M->cd(iCanvas)->SetLogz();
      TH1F* histV0M = (TH1F*)list->FindObject("fHistV0M");
      if(histV0M){
	histV0M->Draw("colz");
      }

      cQARef->cd(iCanvas);
      cQARef->cd(iCanvas)->SetLogz();
      TH1F* histRef = (TH1F*)list->FindObject("fHistRefTracks");
      if(histRef){
	histRef->Draw("colz");
      }
    }
    // ----------------------------------------------------

    // ----------------------------------------------------
    // calculate and plot BF 
    if(listName.Contains("BF_")){
      cBF[iCanvas] = new TCanvas(listName,listName);
      cBF[iCanvas]->Divide(3,3);

      for(Int_t a = 0; a < 7; a++){

	cout<<"ANALYSE "<<gBFAnalysisType[a]<<endl;

	// create the BF object
	bf[iCanvas][a]  = new AliBalance();

	fHistP[a] = (TH1D*)list->FindObject(Form("fHistP%s",gBFAnalysisType[a].Data()));
	fHistN[a] = (TH1D*)list->FindObject(Form("fHistN%s",gBFAnalysisType[a].Data()));
	fHistPP[a] = (TH1D*)list->FindObject(Form("fHistPP%s",gBFAnalysisType[a].Data()));
	fHistPN[a] = (TH1D*)list->FindObject(Form("fHistPN%s",gBFAnalysisType[a].Data()));
	fHistNP[a] = (TH1D*)list->FindObject(Form("fHistNP%s",gBFAnalysisType[a].Data()));
	fHistNN[a] = (TH1D*)list->FindObject(Form("fHistNN%s",gBFAnalysisType[a].Data()));

	// rebin histograms (be careful with divider!)
	fHistP[a]->Rebin(fRebin);
	fHistN[a]->Rebin(fRebin);
	fHistPP[a]->Rebin(fRebin);
	fHistPN[a]->Rebin(fRebin);
	fHistNP[a]->Rebin(fRebin);
	fHistNN[a]->Rebin(fRebin);

	// set histograms in AliBalance object
	bf[iCanvas][a]->SetHistNp(a, fHistP[a]);
	bf[iCanvas][a]->SetHistNn(a, fHistN[a]);
	bf[iCanvas][a]->SetHistNpp(a, fHistPP[a]);
	bf[iCanvas][a]->SetHistNpn(a, fHistPN[a]);
	bf[iCanvas][a]->SetHistNnp(a, fHistNP[a]);
	bf[iCanvas][a]->SetHistNnn(a, fHistNN[a]);

	gbf[iCanvas][a] = bf[iCanvas][a]->GetBalanceFunctionHistogram(a);
	gbf[iCanvas][a]->SetName(Form("%s_BF_%s",listName.Data(),gBFAnalysisType[a].Data()));

	cBF[iCanvas]->cd(a+1);
	gbf[iCanvas][a]->SetMarkerStyle(20);

	if(!bHistos){
	  gbf[iCanvas][a]->DrawCopy("AP");
	  if(a==1) GetWeightedMean(gbf[iCanvas][a],fStartBinBFWidth,WM[iCanvas-1],WME[iCanvas-1]); // for eta recalculate width (from 0.1 only!)
	}
	else{
	  fHistPN[a]->SetLineColor(2);
	  fHistPN[a]->Draw();
	  fHistPP[a]->SetLineColor(1);
	  fHistPP[a]->Draw("same");
	  fHistNP[a]->SetLineColor(4);
	  fHistNP[a]->Draw("same");
	  fHistNN[a]->SetLineColor(8);
	  fHistNN[a]->Draw("same");
	}
      }
    }
    // ----------------------------------------------------

    // ----------------------------------------------------
    // calculate and plot BF (shuffled)
    if(listName.Contains("BFShuffled")){
      cBFS[iCanvas] = new TCanvas(listName,listName);
      cBFS[iCanvas]->Divide(3,3);

      for(Int_t a = 0; a < 7; a++){

	// create the BF object
	bfs[iCanvas][a]  = new AliBalance();

	fHistPS[a] = (TH1D*)list->FindObject(Form("fHistP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistNS[a] = (TH1D*)list->FindObject(Form("fHistP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistPPS[a] = (TH1D*)list->FindObject(Form("fHistPP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistPNS[a] = (TH1D*)list->FindObject(Form("fHistPN%s_shuffle",gBFAnalysisType[a].Data()));
	fHistNPS[a] = (TH1D*)list->FindObject(Form("fHistNP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistNNS[a] = (TH1D*)list->FindObject(Form("fHistNN%s_shuffle",gBFAnalysisType[a].Data()));

	// rebin histograms (be careful with divider!)
	fHistPS[a]->Rebin(fRebin);
	fHistNS[a]->Rebin(fRebin);
	fHistPPS[a]->Rebin(fRebin);
	fHistPNS[a]->Rebin(fRebin);
	fHistNPS[a]->Rebin(fRebin);
	fHistNNS[a]->Rebin(fRebin);

	// set histograms in AliBalance object
	bfs[iCanvas][a]->SetHistNp(a, fHistPS[a]);
	bfs[iCanvas][a]->SetHistNn(a, fHistNS[a]);
	bfs[iCanvas][a]->SetHistNpp(a, fHistPPS[a]);
	bfs[iCanvas][a]->SetHistNpn(a, fHistPNS[a]);
	bfs[iCanvas][a]->SetHistNnp(a, fHistNPS[a]);
	bfs[iCanvas][a]->SetHistNnn(a, fHistNNS[a]);

	gbfs[iCanvas][a] = bfs[iCanvas][a]->GetBalanceFunctionHistogram(a);
	gbfs[iCanvas][a]->SetName(Form("%s_BF_%s",listName.Data(),gBFAnalysisType[a].Data()));

	cBFS[iCanvas]->cd(a+1);
	gbfs[iCanvas][a]->SetMarkerStyle(20);
	if(!bHistos){
	  gbfs[iCanvas][a]->DrawCopy("AP");
	  if(a==1) GetWeightedMean(gbfs[iCanvas][a],fStartBinBFWidth,WMS[iCanvas-1],WMSE[iCanvas-1]); // for eta recalculate width (from 0.1 only!)
	}
	else{
	  fHistPNS[a]->SetLineColor(2);
	  fHistPNS[a]->Draw();
	  fHistPPS[a]->SetLineColor(1);
	  fHistPPS[a]->Draw("same");
	  fHistNPS[a]->SetLineColor(4);
	  fHistNPS[a]->Draw("same");
	  fHistNNS[a]->SetLineColor(8);
	  fHistNNS[a]->Draw("same");
	}
      }
    }
    // ----------------------------------------------------
  }

  // for BF calculation create also graphs with weighted mean for eta
  TGraphErrors *gWM  = NULL;
  TGraphErrors *gWMS = NULL;
  if(!bHistos && gbf[1][1]){
    gWM = new TGraphErrors(iCanvas,cent,WM,centE,WME);
    gWM->SetName("gCentrality");    
  }
  if(!bHistos && gbfs[1][1]){
    gWMS = new TGraphErrors(iCanvas,cent,WMS,centE,WMSE);  
    gWMS->SetName("gCentralityS");
  }

  TFile *fOut = TFile::Open(Form("Histograms_WMstart%d_rebin%d_%s", fStartBinBFWidth, fRebin,inFile.Data()),"RECREATE");
  fOut->cd();
  for(Int_t i = 0; i <= iCanvas; i++){
    for(Int_t a = 0; a < 7; a++){
      if(gbf[i][a]){
	gbf[i][a]->Write();
	gbf[i][a]->Delete();
      }
      if(gbfs[i][a]){
	gbfs[i][a]->Write();
	gbfs[i][a]->Delete();	
      }
    }
  }

  if(gWM) gWM->Write();
  if(gWMS) gWMS->Write();

  fOut->Close();
  f->Close();
  gROOT->Reset();
}

//____________________________________________________________________//
void GetWeightedMean(TH1D *gHistBalance, Int_t fStartBin = 1, Double_t &WM, Double_t &WME) {

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
  for(Int_t i = fStartBin; i <= fNumberOfBins; i++) { 
    cout<<"B: "<<gHistBalance->GetBinContent(i)<<"\t Error: "<<gHistBalance->GetBinError(i)<<"\t bin: "<<gHistBalance->GetBinCenter(i)<<endl;
  } 
  cout<<"=================================================="<<endl;
  for(Int_t i = fStartBin; i <= fNumberOfBins; i++) {
    gSumXi += gHistBalance->GetBinCenter(i);
    gSumBi += gHistBalance->GetBinContent(i);
    gSumBiXi += gHistBalance->GetBinContent(i)*gHistBalance->GetBinCenter(i);
    gSumBiXi2 += gHistBalance->GetBinContent(i)*TMath::Power(gHistBalance->GetBinCenter(i),2);
    gSumBi2Xi2 += TMath::Power(gHistBalance->GetBinContent(i),2)*TMath::Power(gHistBalance->GetBinCenter(i),2);
    gSumDeltaBi2 +=  TMath::Power(gHistBalance->GetBinError(i),2);
    gSumXi2DeltaBi2 += TMath::Power(gHistBalance->GetBinCenter(i),2) * TMath::Power(gHistBalance->GetBinError(i),2);
    
    deltaBalP2 += fP2Step*TMath::Power(gHistBalance->GetBinError(i),2);
    integral += fP2Step*gHistBalance->GetBinContent(i);
  }
  for(Int_t i = 1; i < fNumberOfBins; i++)
    deltaErrorNew += gHistBalance->GetBinError(i)*(gHistBalance->GetBinCenter(i)*gSumBi - gSumBiXi)/TMath::Power(gSumBi,2);
  
  Double_t integralError = TMath::Sqrt(deltaBalP2);
  
  Double_t delta = gSumBiXi / gSumBi;
  Double_t deltaError = (gSumBiXi / gSumBi) * TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumXi2DeltaBi2)/gSumBiXi),2) + TMath::Power((gSumDeltaBi2/gSumBi),2) );
  
  cout<<"Width: "<<delta<<"\t Error: "<<deltaError<<endl;
  cout<<"New error: "<<deltaErrorNew<<endl;
  cout<<"Integral: "<<integral<<"\t Error: "<<integralError<<endl;
  cout<<"=================================================="<<endl;

  WM  = delta;
  WME = deltaError;
}

//___________________________________________________________//
// NOT USED any more
//___________________________________________________________//
void mergeOutput(const char* outputDir) {
  //Function to merge the output of the sub-jobs
  //Create a BF object
  AliBalance *bf = new AliBalance();

  //connect to AliEn's API services
  TGrid::Connect("alien://"); 

  //Getting the output dir from the env. variable 
  //(JDL field: JDLVariables={"OutputDir"};)
  TGridResult* result = gGrid->Query(outputDir,"*/root_archive.zip","","-l 1000");
  
  Int_t nEntries = result->GetEntries();

  TString alienUrl;
  TDirectoryFile *dirSubJob;

  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  TH1F *fHistEventStats = new TH1F("fHistEventStats",
				   "Event statistics;;N_{events}",
				   4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());

  AliESDtrackCuts *bfTrackCuts = new AliESDtrackCuts("bfTrackCuts");
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = result->GetKey(i,"turl");
    alienUrl += "#AnalysisResults.root";
    Printf("Opening file: %s",alienUrl.Data());
    TFile *file = TFile::Open(alienUrl.Data());
    dirSubJob = dynamic_cast<TDirectoryFile *>(file->Get("PWG2EbyE.outputBalanceFunctionAnalysis.root"));

    //merge BF
    AliBalance *bfSubJob = dynamic_cast<AliBalance *>(dirSubJob->Get("AliBalance"));
    //bfSubJob->PrintResults();
    bf->Merge(bfSubJob);
    //delete bfSubJob;

    //merge event stats
    TList *listSubJob = dynamic_cast<TList *>(dirSubJob->Get("listQA"));
    fHistEventStats->Add(dynamic_cast<TH1F *>(listSubJob->At(0)));

    bfTrackCuts = dynamic_cast<AliESDtrackCuts *>(listSubJob->At(1));
    delete listSubJob;
  }

  //Create the output file
  TString outputFile = "AnalysisResults.Merged.root";
  TFile *foutput = TFile::Open(outputFile.Data(),"recreate");
  TDirectoryFile *dirOutput = new TDirectoryFile();
  dirOutput->SetName("PWG2EbyE.outputBalanceFunctionAnalysis.root");
  //dirOutput->cd();
  dirOutput->Add(bf);
  TList *list = new TList();
  list->SetName("listQA");
  list->Add(fHistEventStats);
  list->Add(bfTrackCuts);
  dirOutput->Add(list);
  dirOutput->Write();
  bf->Write();
  list->Write();
  foutput->Close();

    //cout<<alienUrl.Data()<<endl;
}
