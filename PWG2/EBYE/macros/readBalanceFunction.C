const TString gBFAnalysisType[7] = {"y","eta","qlong","qout","qside","qinv","phi"};

const Int_t nrOfCentralities = 9;
const Double_t centralityArray[nrOfCentralities+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};  // in centrality percentile (0-5,5-10,10-20,20-30,...,70-80)
//const Double_t centralityArray[nrOfCentralities+1] = {0.,1.,2.,3.,4.,6.,10.,20.,30.,80.};  // in centrality percentile (0-5,5-10,10-20,20-30,...,70-80)
const Double_t cent[nrOfCentralities]  = {382.8,329.7,260.5,186.4,128.9,85.,52.8,30.,15.8};   // hard coded at the moment for centrality percentiles 
const Double_t centE[nrOfCentralities] = {3.1,4.6,4.4,3.9,3.3,2.6,2.0,1.3,0.6};               // (0-5,5-10,10-20,20-30,...,70-80)

void readBalanceFunction(Bool_t bHistos = kTRUE, TString inFile = "AnalysisResults.root",Int_t fStartBinBFWidth = 1, Int_t fRebin = 1,TString centEst = "V0M",TString extraString = "") {
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
  drawBF(bHistos,inFile, fStartBinBFWidth, fRebin,centEst,extraString);

  //Merge the output
  //mergeOutput("/alice/cern.ch/user/p/pchrist/Balance/pp/7TeV/LHC10b/output/");
}

//___________________________________________________________//
void drawBF(Bool_t bHistos = kTRUE, TString inFile = "AnalysisResults.root", Int_t fStartBinBFWidth = 1, Int_t fRebin = 1, TString centEst = "V0M",TString extraString = "") {
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

  AliBalance *bf[7];
  AliBalance *bfs[7];
  TH1D *gbf[10][7];
  TH1D *gbfs[10][7];

  for(Int_t i = 0; i < 10; i++){
    for(Int_t j = 0; j < 7; j++){
      gbf[i][j] = NULL;
      gbfs[i][j] = NULL;
    }
  }

  TH2D *fHistP[7]; //N+
  TH2D *fHistN[7]; //N-
  TH2D *fHistPN[7]; //N+-
  TH2D *fHistNP[7]; //N-+
  TH2D *fHistPP[7]; //N++
  TH2D *fHistNN[7]; //N--

  TH2D *fHistPS[7]; //N+
  TH2D *fHistNS[7]; //N-
  TH2D *fHistPNS[7]; //N+-
  TH2D *fHistNPS[7]; //N-+
  TH2D *fHistPPS[7]; //N++
  TH2D *fHistNNS[7]; //N--

  Double_t WM[10];     // weighted mean for eta (recalculated from fStartBin)
  Double_t WME[10];    // error
  Double_t WMS[10];     // weighted mean for eta (recalculated from fStartBin) (shuffled)
  Double_t WMSE[10];    // error (shuffled)

  for(iCanvas = 0; iCanvas < nrOfCentralities; iCanvas++){
    
    cBF[iCanvas] = new TCanvas(Form("cBF_%d",iCanvas),Form("cBF_%d",iCanvas));
    cBF[iCanvas]->Divide(3,3);
    
    cBFS[iCanvas] = new TCanvas(Form("Shuffled_%d",iCanvas),Form("Shuffled_%d",iCanvas));
    cBFS[iCanvas]->Divide(3,3);
    
  }

  while ( (key = (TKey*)nextkey())) {

    list = (TList*)key->ReadObj();
    listName = TString(list->GetName());

    cout<<"Processing list "<<listName<<endl;

    // ----------------------------------------------------
    // plot QA histograms
    if(listName.Contains("QA")&&listName.Contains(extraString.Data())){     

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
      cQAV0M->cd(iCanvas+1);
      cQAV0M->cd(iCanvas+1)->SetLogz();
      TH1F* histV0M = (TH1F*)list->FindObject("fHistV0M");
      if(histV0M){
	histV0M->Draw("colz");
      }

      cQARef->cd(iCanvas+1);
      cQARef->cd(iCanvas+1)->SetLogz();
      TH1F* histRef = (TH1F*)list->FindObject("fHistRefTracks");
      if(histRef){
	histRef->Draw("colz");
      }
    }
    // ----------------------------------------------------

    // ----------------------------------------------------
    // calculate and plot BF 
    if(listName.Contains("BF_")&&listName.Contains(extraString.Data())){


      for(Int_t a = 0; a < 7; a++){

	cout<<"ANALYSE "<<gBFAnalysisType[a]<<endl;

	// create the BF object
	bf[a]  = new AliBalance();

	fHistP[a] = (TH2D*)list->FindObject(Form("fHistP%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistN[a] = (TH2D*)list->FindObject(Form("fHistN%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistPP[a] = (TH2D*)list->FindObject(Form("fHistPP%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistPN[a] = (TH2D*)list->FindObject(Form("fHistPN%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistNP[a] = (TH2D*)list->FindObject(Form("fHistNP%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistNN[a] = (TH2D*)list->FindObject(Form("fHistNN%s%s",gBFAnalysisType[a].Data(),centEst.Data()));

	// rebin histograms (be careful with divider!)
	fHistP[a]->RebinY(fRebin);
	fHistN[a]->RebinY(fRebin);
	fHistPP[a]->RebinY(fRebin);
	fHistPN[a]->RebinY(fRebin);
	fHistNP[a]->RebinY(fRebin);
	fHistNN[a]->RebinY(fRebin);

	fHistP[a]->SetName(Form("%s_%d",fHistP[a]->GetName(),iCanvas));
	fHistN[a]->SetName(Form("%s_%d",fHistN[a]->GetName(),iCanvas));
	fHistPP[a]->SetName(Form("%s_%d",fHistPP[a]->GetName(),iCanvas));
	fHistPN[a]->SetName(Form("%s_%d",fHistPN[a]->GetName(),iCanvas));
	fHistNP[a]->SetName(Form("%s_%d",fHistNP[a]->GetName(),iCanvas));
	fHistNN[a]->SetName(Form("%s_%d",fHistNN[a]->GetName(),iCanvas));

	// set histograms in AliBalance object
	bf[a]->SetHistNp(a, fHistP[a]);
	bf[a]->SetHistNn(a, fHistN[a]);
	bf[a]->SetHistNpp(a, fHistPP[a]);
	bf[a]->SetHistNpn(a, fHistPN[a]);
	bf[a]->SetHistNnp(a, fHistNP[a]);
	bf[a]->SetHistNnn(a, fHistNN[a]);

  	for(iCanvas = 0; iCanvas < nrOfCentralities; iCanvas++){

	  gbf[iCanvas][a] = bf[a]->GetBalanceFunctionHistogram(a,centralityArray[iCanvas],centralityArray[iCanvas+1]);
	  gbf[iCanvas][a]->SetName(Form("BF_%s_Cent_%.0f_%.0f",gBFAnalysisType[a].Data(),centralityArray[iCanvas],centralityArray[iCanvas+1]));

	  cBF[iCanvas]->cd(a+1);
	  gbf[iCanvas][a]->SetMarkerStyle(20);

	  if(!bHistos){
	    gbf[iCanvas][a]->DrawCopy("AP");
	    if(a==1) GetWeightedMean(gbf[iCanvas][a],fStartBinBFWidth,WM[iCanvas],WME[iCanvas]); // for eta recalculate width (from 0.1 only!)
	  }
	  else{
	    fHistPN[a]->SetLineColor(2);
	    fHistPN[a]->ProjectionY(Form("pn%d",a))->DrawCopy();
	    fHistPP[a]->SetLineColor(1);
	    fHistPP[a]->ProjectionY(Form("pp%d",a))->DrawCopy("same");
	    fHistNP[a]->SetLineColor(4);
	    fHistNP[a]->ProjectionY(Form("np%d",a))->DrawCopy("same");
	    fHistNN[a]->SetLineColor(8);
	    fHistNN[a]->ProjectionY(Form("nn%d",a))->DrawCopy("same");
	  }
	}
      }
    }
    // ----------------------------------------------------

    // ----------------------------------------------------
    // calculate and plot BF (shuffled)
    if(listName.Contains("BFShuffled")&&listName.Contains(extraString.Data())){

      for(Int_t a = 0; a < 7; a++){

	// create the BF object
	bfs[a]  = new AliBalance();

	fHistPS[a] = (TH2D*)list->FindObject(Form("fHistP%s%s_shuffle",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistNS[a] = (TH2D*)list->FindObject(Form("fHistN%s%s_shuffle",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistPPS[a] = (TH2D*)list->FindObject(Form("fHistPP%s%s_shuffle",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistPNS[a] = (TH2D*)list->FindObject(Form("fHistPN%s%s_shuffle",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistNPS[a] = (TH2D*)list->FindObject(Form("fHistNP%s%s_shuffle",gBFAnalysisType[a].Data(),centEst.Data()));
	fHistNNS[a] = (TH2D*)list->FindObject(Form("fHistNN%s%s_shuffle",gBFAnalysisType[a].Data(),centEst.Data()));

	// rebin histograms (be careful with divider!)
	fHistPS[a]->RebinY(fRebin);
	fHistNS[a]->RebinY(fRebin);
	fHistPPS[a]->RebinY(fRebin);
	fHistPNS[a]->RebinY(fRebin);
	fHistNPS[a]->RebinY(fRebin);
	fHistNNS[a]->RebinY(fRebin);

	fHistPS[a]->SetName(Form("%s_%d",fHistPS[a]->GetName(),iCanvas));
	fHistNS[a]->SetName(Form("%s_%d",fHistNS[a]->GetName(),iCanvas));
	fHistPPS[a]->SetName(Form("%s_%d",fHistPPS[a]->GetName(),iCanvas));
	fHistPNS[a]->SetName(Form("%s_%d",fHistPNS[a]->GetName(),iCanvas));
	fHistNPS[a]->SetName(Form("%s_%d",fHistNPS[a]->GetName(),iCanvas));
	fHistNNS[a]->SetName(Form("%s_%d",fHistNNS[a]->GetName(),iCanvas));

	// set histograms in AliBalance object
	bfs[a]->SetHistNp(a, fHistPS[a]);
	bfs[a]->SetHistNn(a, fHistNS[a]);
	bfs[a]->SetHistNpp(a, fHistPPS[a]);
	bfs[a]->SetHistNpn(a, fHistPNS[a]);
	bfs[a]->SetHistNnp(a, fHistNPS[a]);
	bfs[a]->SetHistNnn(a, fHistNNS[a]);

	for(iCanvas = 0; iCanvas < nrOfCentralities; iCanvas++){
	  
	  gbfs[iCanvas][a] = bfs[a]->GetBalanceFunctionHistogram(a,centralityArray[iCanvas],centralityArray[iCanvas+1]);
	  gbf[iCanvas][a]->SetName(Form("BFS_%s_Cent_%.0f_%.0f",gBFAnalysisType[a].Data(),centralityArray[iCanvas],centralityArray[iCanvas+1]));
	  
	  cBFS[iCanvas]->cd(a+1);
	  gbfs[iCanvas][a]->SetMarkerStyle(20);
	  if(!bHistos){
	    gbfs[iCanvas][a]->DrawCopy("AP");
	    if(a==1) GetWeightedMean(gbfs[iCanvas][a],fStartBinBFWidth,WMS[iCanvas],WMSE[iCanvas]); // for eta recalculate width (from 0.1 only!)
	  }
	  else{
	    fHistPNS[a]->SetLineColor(2);
	    fHistPNS[a]->ProjectionY(Form("pns%d",a))->DrawCopy();
	    fHistPPS[a]->SetLineColor(1);
	    fHistPPS[a]->ProjectionY(Form("pps%d",a))->DrawCopy("same");
	    fHistNPS[a]->SetLineColor(4);
	    fHistNPS[a]->ProjectionY(Form("nps%d",a))->DrawCopy("same");
	    fHistNNS[a]->SetLineColor(8);
	    fHistNNS[a]->ProjectionY(Form("nns%d",a))->DrawCopy("same");
	  }
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
  for(Int_t a = 0; a < 7; a++){

    if(fHistPN[a]){
      (fHistPN[a]->ProjectionY(Form("hPN_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistPP[a]->ProjectionY(Form("hPP_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistNP[a]->ProjectionY(Form("hNP_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistNN[a]->ProjectionY(Form("hNN_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistP[a]->ProjectionY(Form("hP_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistN[a]->ProjectionY(Form("hN_%s",gBFAnalysisType[a].Data())))->Write();
    }

    if(fHistPNS[a]){
      (fHistPNS[a]->ProjectionY(Form("hPNS_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistPPS[a]->ProjectionY(Form("hPPS_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistNPS[a]->ProjectionY(Form("hNPS_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistNNS[a]->ProjectionY(Form("hNNS_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistPS[a]->ProjectionY(Form("hPS_%s",gBFAnalysisType[a].Data())))->Write();
      (fHistNS[a]->ProjectionY(Form("hNS_%s",gBFAnalysisType[a].Data())))->Write();
    }

    for(Int_t i = 0; i < iCanvas; i++){
      
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
