const TString gBFAnalysisType[7] = {"y","eta","qlong","qout","qside","qinv","phi"};

const Int_t nrOfCentralities = 9;
const Double_t centralityArray[nrOfCentralities+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};  // in centrality percentile (0-5,5-10,10-20,20-30,...,70-80)
//const Double_t centralityArray[nrOfCentralities+1] = {0.,1.,2.,3.,4.,6.,10.,20.,30.,80.};  // in centrality percentile (0-5,5-10,10-20,20-30,...,70-80)
const Double_t cent[nrOfCentralities]  = {382.8,329.7,260.5,186.4,128.9,85.,52.8,30.,15.8};   // hard coded at the moment for centrality percentiles 
const Double_t centE[nrOfCentralities] = {3.1,4.6,4.4,3.9,3.3,2.6,2.0,1.3,0.6};               // (0-5,5-10,10-20,20-30,...,70-80)

void readBalanceFunction(Bool_t bHistos = kFALSE, TString inFile = "AnalysisResults.root",Int_t fStartBinBFWidth = 3, Int_t fRebin = 2,Int_t fStartBinBFWidthPhi = 2, Int_t fRebinPhi = 2,TString centEst = "V0M",Double_t etaWindow = -1) {
  // Macro to read the output of the BF analysis:  MW: CHANGE THIS!!!!
  //i) Prints and draws the final BF output
  //ii) Plots the QA part of the analysis
  //iii) store BF in output file
  //Author: Panos.Christakoglou@cern.ch, m.weber@cern.ch
  //Loading the needed libraries
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGCFebye.so");

  //Draw BF       
  drawBF(bHistos,inFile, fStartBinBFWidth, fRebin,fStartBinBFWidthPhi, fRebinPhi,centEst, "",  etaWindow);    
  

  //Merge the output
  //mergeOutput("/alice/cern.ch/user/p/pchrist/Balance/pp/7TeV/LHC10b/output/");
}

//___________________________________________________________//
void drawBF(Bool_t bHistos = kFALSE, TString inFile = "AnalysisResults.root", Int_t fStartBinBFWidth = 1, Int_t fRebin = 1, Int_t fStartBinBFWidthPhi = 1, Int_t fRebinPhi = 1, TString centEst = "V0M",TString extraString = "", Double_t etaWindow = -1) {
  //Function to draw the BF objects and write them into the output file

  Int_t maximumCanvases = 13;
  Int_t iCanvas         = 0;
  Int_t iList           = -1;
  TCanvas *cQA[13][10];
  TCanvas *cQAV0M = new TCanvas("cQAV0M","V0M multiplicities");
  cQAV0M->Divide(4,3);
  TCanvas *cQARef = new TCanvas("cQARef","reference track multiplicities");
  cQARef->Divide(4,3);
  TCanvas *cBF[13][10];
  TCanvas *cBFS[13][10];

  // get the file
  TFile *f = TFile::Open(inFile.Data());
  if(!f) {
    Printf("File not found!!!");
    break;
  }

  // get the BF output directory
  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("PWGCFEbyE.outputBalanceFunctionAnalysis"));
  if(!dir) {
    Printf("Output directory not found!!!");
    break;
  }

  // loop over all lists and plot the BF and QA
  TList *list = NULL;
  TString listName;
  TIter nextkey( dir->GetListOfKeys() );
  TKey *key;

  AliBalance *bf[13][7];
  AliBalance *bfs[13][7];
  TH1D *gbf[13][10][7];
  TH1D *gbfs[13][10][7];

  for(Int_t i = 0; i < 13; i++){
    for(Int_t j = 0; j < 10; j++){
      for(Int_t k = 0; k < 7; k++){
	gbf[i][j][k]  = NULL;
	gbfs[i][j][k] = NULL;
      }
    }
  }

  TH2D *fHistP[13][7]; //N+
  TH2D *fHistN[13][7]; //N-
  TH2D *fHistPN[13][7]; //N+-
  TH2D *fHistNP[13][7]; //N-+
  TH2D *fHistPP[13][7]; //N++
  TH2D *fHistNN[13][7]; //N--

  TH2D *fHistPS[13][7]; //N+
  TH2D *fHistNS[13][7]; //N-
  TH2D *fHistPNS[13][7]; //N+-
  TH2D *fHistNPS[13][7]; //N-+
  TH2D *fHistPPS[13][7]; //N++
  TH2D *fHistNNS[13][7]; //N--

  Double_t WM[13][10];     // weighted mean for eta (recalculated from fStartBin)
  Double_t WME[13][10];    // error
  Double_t WMS[13][10];     // weighted mean for eta (recalculated from fStartBin) (shuffled)
  Double_t WMSE[13][10];    // error (shuffled)

  Double_t WMP[13][10];     // weighted mean for phi (recalculated from fStartBinPhi)
  Double_t WMPE[13][10];    // error
  Double_t WMPS[13][10];     // weighted mean for phi (recalculated from fStartBin) (shuffled)
  Double_t WMPSE[13][10];    // error (shuffled)

  Double_t integ[13][10];     // integral for eta (calculated from bin 1)
  Double_t integE[13][10];    // error
  Double_t integS[13][10];     // integral for eta (calculated from bin 1) (shuffled)
  Double_t integSE[13][10];    // error (shuffled)

  Double_t integP[13][10];     // integral for phi (calculated from bin 1)
  Double_t integPE[13][10];    // error
  Double_t integPS[13][10];     // integral for phi (calculated from bin 1) (shuffled)
  Double_t integPSE[13][10];    // error (shuffled)


  while ( (key = (TKey*)nextkey())) {

    list = (TList*)key->ReadObj();
    listName = TString(list->GetName());

    cout<<"Processing list "<<listName<<endl;

 
    // ----------------------------------------------------
    // plot QA histograms
    if(listName.Contains("QA")){   

      if(iList<13) iList++;
      else{
	cerr<<"TOO MANY LISTS!!!"<<endl;
	return;
      }  
      
      cQA[iList][iCanvas] = new TCanvas(listName,listName);
      cQA[iList][iCanvas]->Divide(4,3);

     cQA[iList][iCanvas]->cd(1);
      TH1F* histEventStats = (TH1F*)list->FindObject("fHistEventStats");
      if(histEventStats){
	histEventStats->SetFillColor(9);
	histEventStats->Draw();
      }

      cQA[iList][iCanvas]->cd(2);
      TH1F* histTriggerStats = (TH1F*)list->FindObject("fHistTriggerStats");
      if(histTriggerStats){
	histTriggerStats->SetFillColor(9);
	histTriggerStats->Draw();
      }

     cQA[iList][iCanvas]->cd(3);
      TH1F* histTrackStats = (TH1F*)list->FindObject("fHistTrackStats");
      if(histTrackStats){
	histTrackStats->SetFillColor(9);
	histTrackStats->Draw();
      }

      cQA[iList][iCanvas]->cd(4);
      TH1F* histCentStats = (TH1F*)list->FindObject("fHistCentStats");
      if(histCentStats){
	histCentStats->SetFillColor(9);
	histCentStats->Draw("colz");
      }



      cQA[iList][iCanvas]->cd(5);
      TH1F* histVx = (TH1F*)list->FindObject("fHistVx");
      if(histVx){
	histVx->SetFillColor(9);
	histVx->Draw();
      }
      cQA[iList][iCanvas]->cd(6);
      TH1F* histVy = (TH1F*)list->FindObject("fHistVy");
      if(histVy){
	histVy->SetFillColor(9);
	histVy->Draw();
      }
      
      cQA[iList][iCanvas]->cd(7);
      TH1F* histVz = (TH1F*)list->FindObject("fHistVz");
      if(histVz){
	histVz->SetFillColor(9);
	histVz->Draw();
      }

      cQA[iList][iCanvas]->cd(8);
      cQA[iList][iCanvas]->cd(8)->SetLogz();
      TH2F* histDCA = (TH2F*)list->FindObject("fHistDCA");  
      if(histDCA) histDCA->Draw("colz");
      

      cQA[iList][iCanvas]->cd(9);
      cQA[iList][iCanvas]->cd(9)->SetLogz();
      TH2F* histClus = (TH2F*)list->FindObject("fHistClus");
      if(histClus) histClus->Draw("colz");
      
      cQA[iList][iCanvas]->cd(10);
      TH1F* histPt = (TH1F*)list->FindObject("fHistPt");
      if(histPt){
	histPt->SetFillColor(9);
	histPt->Draw();
      }
      
      cQA[iList][iCanvas]->cd(11);
      TH1F* histEta = (TH1F*)list->FindObject("fHistEta");
      if(histEta){
	histEta->SetFillColor(9);
	histEta->Draw();
      }
      
      cQA[iList][iCanvas]->cd(12);
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
    if(listName.Contains("BF_")&&listName.Contains(centEst.Data())){

      for(iCanvas = 0; iCanvas < nrOfCentralities; iCanvas++){
	
	cBF[iList][iCanvas] = new TCanvas(Form("cBF_%d_%d",iList,iCanvas),Form("cBF_%d_%d",iList,iCanvas));
	cBF[iList][iCanvas]->Divide(3,3);
	
      }
      

      for(Int_t a = 0; a < 7; a++){

    	cout<<"ANALYSE "<<gBFAnalysisType[a]<<endl;

    	// create the BF object
    	bf[iList][a]  = new AliBalance();

    	fHistP[iList][a] = (TH2D*)list->FindObject(Form("fHistP%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
    	fHistN[iList][a] = (TH2D*)list->FindObject(Form("fHistN%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
    	fHistPP[iList][a] = (TH2D*)list->FindObject(Form("fHistPP%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
    	fHistPN[iList][a] = (TH2D*)list->FindObject(Form("fHistPN%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
    	fHistNP[iList][a] = (TH2D*)list->FindObject(Form("fHistNP%s%s",gBFAnalysisType[a].Data(),centEst.Data()));
    	fHistNN[iList][a] = (TH2D*)list->FindObject(Form("fHistNN%s%s",gBFAnalysisType[a].Data(),centEst.Data()));

    	// rebin histograms (be careful with divider!)
	if(a==6){
	  fHistP[iList][a]->RebinY(fRebinPhi);
	  fHistN[iList][a]->RebinY(fRebinPhi);
	  fHistPP[iList][a]->RebinY(fRebinPhi);
	  fHistPN[iList][a]->RebinY(fRebinPhi);
	  fHistNP[iList][a]->RebinY(fRebinPhi);
	  fHistNN[iList][a]->RebinY(fRebinPhi);
	}
	else{
	  fHistP[iList][a]->RebinY(fRebin);
	  fHistN[iList][a]->RebinY(fRebin);
	  fHistPP[iList][a]->RebinY(fRebin);
	  fHistPN[iList][a]->RebinY(fRebin);
	  fHistNP[iList][a]->RebinY(fRebin);
	  fHistNN[iList][a]->RebinY(fRebin);
	}

    	fHistP[iList][a]->SetName(Form("%s_%d_%d",fHistP[iList][a]->GetName(),iList,iCanvas));
    	fHistN[iList][a]->SetName(Form("%s_%d_%d",fHistN[iList][a]->GetName(),iList,iCanvas));
    	fHistPP[iList][a]->SetName(Form("%s_%d_%d",fHistPP[iList][a]->GetName(),iList,iCanvas));
    	fHistPN[iList][a]->SetName(Form("%s_%d_%d",fHistPN[iList][a]->GetName(),iList,iCanvas));
    	fHistNP[iList][a]->SetName(Form("%s_%d_%d",fHistNP[iList][a]->GetName(),iList,iCanvas));
    	fHistNN[iList][a]->SetName(Form("%s_%d_%d",fHistNN[iList][a]->GetName(),iList,iCanvas));

    	// set histograms in AliBalance object
    	bf[iList][a]->SetHistNp(a, fHistP[iList][a]);
    	bf[iList][a]->SetHistNn(a, fHistN[iList][a]);
    	bf[iList][a]->SetHistNpp(a, fHistPP[iList][a]);
    	bf[iList][a]->SetHistNpn(a, fHistPN[iList][a]);
    	bf[iList][a]->SetHistNnp(a, fHistNP[iList][a]);
    	bf[iList][a]->SetHistNnn(a, fHistNN[iList][a]);

    	for(iCanvas = 0; iCanvas < nrOfCentralities; iCanvas++){

    	  gbf[iList][iCanvas][a] = bf[iList][a]->GetBalanceFunctionHistogram(a,centralityArray[iCanvas],centralityArray[iCanvas+1],etaWindow);
    	  gbf[iList][iCanvas][a]->SetName(Form("BF_%s_Cent_%.0f_%.0f_%d",gBFAnalysisType[a].Data(),centralityArray[iCanvas],centralityArray[iCanvas+1],iList));

    	  cBF[iList][iCanvas]->cd(a+1);
    	  gbf[iList][iCanvas][a]->SetMarkerStyle(20);

    	  if(!bHistos){
    	    gbf[iList][iCanvas][a]->DrawCopy("AP");
    	    if(a==1){
	      GetWeightedMean(gbf[iList][iCanvas][a],fStartBinBFWidth,WM[iList][iCanvas],WME[iList][iCanvas]); // for eta recalculate width 
	      GetIntegral(gbf[iList][iCanvas][a],integ[iList][iCanvas],integE[iList][iCanvas]);
	    }
	    else if(a==6){
	      GetWeightedMean(gbf[iList][iCanvas][a],fStartBinBFWidthPhi,WMP[iList][iCanvas],WMPE[iList][iCanvas]); // for phi calculate width 
	      GetIntegral(gbf[iList][iCanvas][a],integP[iList][iCanvas],integPE[iList][iCanvas]);
	    }
	  }
    	  else{
    	    fHistPN[iList][a]->SetLineColor(2);
    	    fHistPN[iList][a]->ProjectionY(Form("pn%d",a))->DrawCopy();
    	    fHistPP[iList][a]->SetLineColor(1);
    	    fHistPP[iList][a]->ProjectionY(Form("pp%d",a))->DrawCopy("same");
    	    fHistNP[iList][a]->SetLineColor(4);
    	    fHistNP[iList][a]->ProjectionY(Form("np%d",a))->DrawCopy("same");
    	    fHistNN[iList][a]->SetLineColor(8);
    	    fHistNN[iList][a]->ProjectionY(Form("nn%d",a))->DrawCopy("same");
    	  }
    	}
      }
    }
    // ----------------------------------------------------

    // ----------------------------------------------------
    // calculate and plot BF (shuffled)
    if(listName.Contains("BFShuffled")&&listName.Contains(centEst.Data())&&listName.Contains(extraString.Data())){

      for(iCanvas = 0; iCanvas < nrOfCentralities; iCanvas++){
	
	cBFS[iList][iCanvas] = new TCanvas(Form("Shuffled_%d_%d",iList,iCanvas),Form("Shuffled_%d_%d",iList,iCanvas));
	cBFS[iList][iCanvas]->Divide(3,3);
	
      }

      for(Int_t a = 0; a < 7; a++){

  	// create the BF object
  	bfs[iList][a]  = new AliBalance();

  	fHistPS[iList][a] = (TH2D*)list->FindObject(Form("fHistP%s_shuffle%s",gBFAnalysisType[a].Data(),centEst.Data()));
  	fHistNS[iList][a] = (TH2D*)list->FindObject(Form("fHistN%s_shuffle%s",gBFAnalysisType[a].Data(),centEst.Data()));
  	fHistPPS[iList][a] = (TH2D*)list->FindObject(Form("fHistPP%s_shuffle%s",gBFAnalysisType[a].Data(),centEst.Data()));
  	fHistPNS[iList][a] = (TH2D*)list->FindObject(Form("fHistPN%s_shuffle%s",gBFAnalysisType[a].Data(),centEst.Data()));
  	fHistNPS[iList][a] = (TH2D*)list->FindObject(Form("fHistNP%s_shuffle%s",gBFAnalysisType[a].Data(),centEst.Data()));
  	fHistNNS[iList][a] = (TH2D*)list->FindObject(Form("fHistNN%s_shuffle%s",gBFAnalysisType[a].Data(),centEst.Data()));

  	// rebin histograms (be careful with divider!)
	if(a==6){
	  fHistPS[iList][a]->RebinY(fRebinPhi);
	  fHistNS[iList][a]->RebinY(fRebinPhi);
	  fHistPPS[iList][a]->RebinY(fRebinPhi);
	  fHistPNS[iList][a]->RebinY(fRebinPhi);
	  fHistNPS[iList][a]->RebinY(fRebinPhi);
	  fHistNNS[iList][a]->RebinY(fRebinPhi);
	}
	else{
	  fHistPS[iList][a]->RebinY(fRebin);
	  fHistNS[iList][a]->RebinY(fRebin);
	  fHistPPS[iList][a]->RebinY(fRebin);
	  fHistPNS[iList][a]->RebinY(fRebin);
	  fHistNPS[iList][a]->RebinY(fRebin);
	  fHistNNS[iList][a]->RebinY(fRebin);
	}

  	fHistPS[iList][a]->SetName(Form("%s_%d_%d",fHistPS[iList][a]->GetName(),iList,iCanvas));
  	fHistNS[iList][a]->SetName(Form("%s_%d_%d",fHistNS[iList][a]->GetName(),iList,iCanvas));
  	fHistPPS[iList][a]->SetName(Form("%s_%d_%d",fHistPPS[iList][a]->GetName(),iList,iCanvas));
  	fHistPNS[iList][a]->SetName(Form("%s_%d_%d",fHistPNS[iList][a]->GetName(),iList,iCanvas));
  	fHistNPS[iList][a]->SetName(Form("%s_%d_%d",fHistNPS[iList][a]->GetName(),iList,iCanvas));
  	fHistNNS[iList][a]->SetName(Form("%s_%d_%d",fHistNNS[iList][a]->GetName(),iList,iCanvas));

  	// set histograms in AliBalance object
  	bfs[iList][a]->SetHistNp(a, fHistPS[iList][a]);
  	bfs[iList][a]->SetHistNn(a, fHistNS[iList][a]);
  	bfs[iList][a]->SetHistNpp(a, fHistPPS[iList][a]);
  	bfs[iList][a]->SetHistNpn(a, fHistPNS[iList][a]);
  	bfs[iList][a]->SetHistNnp(a, fHistNPS[iList][a]);
  	bfs[iList][a]->SetHistNnn(a, fHistNNS[iList][a]);

  	for(iCanvas = 0; iCanvas < nrOfCentralities; iCanvas++){
	  
  	  gbfs[iList][iCanvas][a] = bfs[iList][a]->GetBalanceFunctionHistogram(a,centralityArray[iCanvas],centralityArray[iCanvas+1],etaWindow);
  	  gbfs[iList][iCanvas][a]->SetName(Form("BFS_%s_Cent_%.0f_%.0f_%d",gBFAnalysisType[a].Data(),centralityArray[iCanvas],centralityArray[iCanvas+1],iList));
	  
  	  cBFS[iList][iCanvas]->cd(a+1);
  	  gbfs[iList][iCanvas][a]->SetMarkerStyle(20);
  	  if(!bHistos){
  	    gbfs[iList][iCanvas][a]->DrawCopy("AP");
  	    if(a==1){
	      GetWeightedMean(gbfs[iList][iCanvas][a],fStartBinBFWidth,WMS[iList][iCanvas],WMSE[iList][iCanvas]); 
	      GetIntegral(gbfs[iList][iCanvas][a],integS[iList][iCanvas],integSE[iList][iCanvas]); 
	    }
	    else if(a==6){
	      GetWeightedMean(gbfs[iList][iCanvas][a],fStartBinBFWidthPhi,WMPS[iList][iCanvas],WMPSE[iList][iCanvas]); // for phi calculate width 
	      GetIntegral(gbfs[iList][iCanvas][a],integPS[iList][iCanvas],integPSE[iList][iCanvas]); 
	    }
  	  }
  	  else{
  	    fHistPNS[iList][a]->SetLineColor(2);
  	    fHistPNS[iList][a]->ProjectionY(Form("pns%d",a))->DrawCopy();
  	    fHistPPS[iList][a]->SetLineColor(1);
  	    fHistPPS[iList][a]->ProjectionY(Form("pps%d",a))->DrawCopy("same");
  	    fHistNPS[iList][a]->SetLineColor(4);
  	    fHistNPS[iList][a]->ProjectionY(Form("nps%d",a))->DrawCopy("same");
  	    fHistNNS[iList][a]->SetLineColor(8);
  	    fHistNNS[iList][a]->ProjectionY(Form("nns%d",a))->DrawCopy("same");
  	  }
  	}
      }
    }
    // ----------------------------------------------------
  }
  
  // for BF calculation create also graphs with weighted mean for eta
  if(iCanvas == 0) return; 

  TGraphErrors *gWM[13];
  TGraphErrors *gWMS[13];
  TGraphErrors *gWMP[13];
  TGraphErrors *gWMPS[13];
  TGraphErrors *ginteg[13];
  TGraphErrors *gintegS[13];
  TGraphErrors *gintegP[13];
  TGraphErrors *gintegPS[13];
  for(Int_t i = 0; i < 13; i++){
    gWM[i] = NULL;
    gWMS[i] = NULL;
    ginteg[i] = NULL;
    gintegS[i] = NULL;
    gWMP[i] = NULL;
    gWMPS[i] = NULL;
    gintegP[i] = NULL;
    gintegPS[i] = NULL;
  }


  if(!bHistos){
    for(Int_t i = 0; i < iList+1; i++){
      gWM[i] = new TGraphErrors(iCanvas,cent,WM[i],centE,WME[i]);
      if(iList==0) gWM[i]->SetName("gCentrality");  
      else gWM[i]->SetName(Form("gCentrality_%d",i));
      gWMS[i] = new TGraphErrors(iCanvas,cent,WMS[i],centE,WMSE[i]); 
      if(iList==0) gWMS[i]->SetName("gCentralityS");  
      else gWMS[i]->SetName(Form("gCentralityS_%d",i)); 

      gWMP[i] = new TGraphErrors(iCanvas,cent,WMP[i],centE,WMPE[i]);
      if(iList==0) gWMP[i]->SetName("gCentralityPhi");  
      else gWMP[i]->SetName(Form("gCentralityPhi_%d",i));
      gWMPS[i] = new TGraphErrors(iCanvas,cent,WMPS[i],centE,WMPSE[i]); 
      if(iList==0) gWMPS[i]->SetName("gCentralityPhiS");  
      else gWMPS[i]->SetName(Form("gCentralityPhiS_%d",i)); 

      ginteg[i] = new TGraphErrors(iCanvas,cent,integ[i],centE,integE[i]);
      if(iList==0) ginteg[i]->SetName("gIntegral");  
      else ginteg[i]->SetName(Form("gIntegral_%d",i));
      gintegS[i] = new TGraphErrors(iCanvas,cent,integS[i],centE,integSE[i]); 
      if(iList==0) gintegS[i]->SetName("gIntegralS");  
      else gintegS[i]->SetName(Form("gIntegralS_%d",i)); 

      gintegP[i] = new TGraphErrors(iCanvas,cent,integP[i],centE,integPE[i]);
      if(iList==0) gintegP[i]->SetName("gIntegraPhil");  
      else gintegP[i]->SetName(Form("gIntegralPhi_%d",i));
      gintegPS[i] = new TGraphErrors(iCanvas,cent,integPS[i],centE,integPSE[i]); 
      if(iList==0) gintegPS[i]->SetName("gIntegralPhiS");  
      else gintegPS[i]->SetName(Form("gIntegralPhiS_%d",i)); 
    }
  }

  TFile *fOut = TFile::Open(Form("Histograms_WMstart%d_rebin%d_WMstartPhi%d_rebinPhi%d_%s", fStartBinBFWidth, fRebin,fStartBinBFWidthPhi, fRebinPhi,inFile.Data()),"RECREATE");
  fOut->cd();
  for(Int_t i = 0; i < iList+1; i++){
    cout<<"PROCESS LIST "<<i<<" NOW!"<<endl;
    for(Int_t a = 0; a < 7; a++){
    cout<<"PROCESS VARIABLE "<<a<<" NOW!"<<endl;
      
      if(fHistPN[i][a]){
	(fHistPN[i][a]->ProjectionY(Form("hPN_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistPP[i][a]->ProjectionY(Form("hPP_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistNP[i][a]->ProjectionY(Form("hNP_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistNN[i][a]->ProjectionY(Form("hNN_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistP[i][a]->ProjectionY(Form("hP_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistN[i][a]->ProjectionY(Form("hN_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
      }

      if(fHistPNS[i][a]){
	(fHistPNS[i][a]->ProjectionY(Form("hPNS_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistPPS[i][a]->ProjectionY(Form("hPPS_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistNPS[i][a]->ProjectionY(Form("hNPS_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistNNS[i][a]->ProjectionY(Form("hNNS_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistPS[i][a]->ProjectionY(Form("hPS_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
	(fHistNS[i][a]->ProjectionY(Form("hNS_%s_%d",gBFAnalysisType[a].Data(),i)))->Write();
      }

      //printout in text format for delta eta      
      for(Int_t j = 0; j < iCanvas; j++){
	cout<<"//=========================Centrality "<<centralityArray[j]<<"-"<<centralityArray[j+1]<<"%==================//"<<endl;
	
	if(gbf[i][j][a]){
	  if(a==1){
	    cout<<"Double_t gALICEDataBalanceFunctionInDeltaEtaCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"[nBinsInDeltaEta] = {";
	    for(Int_t k = 0; k < gbf[i][j][a]->GetNbinsX()-1;k++) cout<<gbf[i][j][a]->GetBinContent(k+1)<<", ";
	    cout<<gbf[i][j][a]->GetBinContent(k+1)<<"};"<<endl;
	    cout<<"Double_t gALICEDataBalanceFunctionInDeltaEtaCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"Error[nBinsInDeltaEta] = {";
	    for(Int_t k = 0; k < gbf[i][j][a]->GetNbinsX()-1;k++) cout<<gbf[i][j][a]->GetBinError(k+1)<<", ";
	    cout<<gbf[i][j][a]->GetBinError(k+1)<<"};"<<endl;
	  } 
	  else if(a==6){
	    cout<<"Double_t gALICEDataBalanceFunctionInDeltaPhiCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"[nBinsInDeltaPhi] = {";
	    for(Int_t k = 0; k < gbf[i][j][a]->GetNbinsX()-1;k++) cout<<gbf[i][j][a]->GetBinContent(k+1)<<", ";
	    cout<<gbf[i][j][a]->GetBinContent(k+1)<<"};"<<endl;
	    cout<<"Double_t gALICEDataBalanceFunctionInDeltaPhiCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"Error[nBinsInDeltaPhi] = {";
	    for(Int_t k = 0; k < gbf[i][j][a]->GetNbinsX()-1;k++) cout<<gbf[i][j][a]->GetBinError(k+1)<<", ";
	    cout<<gbf[i][j][a]->GetBinError(k+1)<<"};"<<endl;
	  } 
	  gbf[i][j][a]->Write();
	  gbf[i][j][a]->Delete();
	}
	if(gbfs[i][j][a]){
	  if(a==1){
	    cout<<"Double_t gALICEShuffledDataBalanceFunctionInDeltaEtaCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"[nBinsInDeltaEta] = {";
	    for(Int_t k = 0; k < gbfs[i][j][a]->GetNbinsX()-1;k++) cout<<gbfs[i][j][a]->GetBinContent(k+1)<<", ";
	    cout<<gbfs[i][j][a]->GetBinContent(k+1)<<"};"<<endl;
	    cout<<"Double_t gALICEShuffledDataBalanceFunctionInDeltaEtaCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"Error[nBinsInDeltaEta] = {";
	    for(Int_t k = 0; k < gbfs[i][j][a]->GetNbinsX()-1;k++) cout<<gbfs[i][j][a]->GetBinError(k+1)<<", ";
	    cout<<gbfs[i][j][a]->GetBinError(k+1)<<"};"<<endl;
	  } 
	  else if(a==6){
	    cout<<"Double_t gALICEShuffledDataBalanceFunctionInDeltaPhiCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"[nBinsInDeltaPhi] = {";
	    for(Int_t k = 0; k < gbfs[i][j][a]->GetNbinsX()-1;k++) cout<<gbfs[i][j][a]->GetBinContent(k+1)<<", ";
	    cout<<gbfs[i][j][a]->GetBinContent(k+1)<<"};"<<endl;
	    cout<<"Double_t gALICEShuffledDataBalanceFunctionInDeltaPhiCentrality"<<centralityArray[j]<<"to"<<centralityArray[j+1]<<"Error[nBinsInDeltaPhi] = {";
	    for(Int_t k = 0; k < gbfs[i][j][a]->GetNbinsX()-1;k++) cout<<gbfs[i][j][a]->GetBinError(k+1)<<", ";
	    cout<<gbfs[i][j][a]->GetBinError(k+1)<<"};"<<endl;
	  } 
	  gbfs[i][j][a]->Write();
	  gbfs[i][j][a]->Delete();
	}
	cout<<"//=========================Centrality "<<centralityArray[j]<<"-"<<centralityArray[j+1]<<"%==================//"<<endl;
	cout<<endl;
      }
    }

    Double_t x,y;

    cout<<"//================================ALICE================================//"<<endl;
    if(gWM[i]){
      cout<<"Double_t gWeightedMeanInEtaAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gWM[i]->GetN()-1;k++){
	gWM[i]->GetPoint(k,x,y);
	cout<<y<<", ";      
      }
      gWM[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;    
      
      cout<<"Double_t gWeightedMeanInEtaAliceError[nCentralityBins] = {";
      for(Int_t k = 0; k < gWM[i]->GetN()-1;k++){
	cout<<gWM[i]->GetErrorY(k)<<", ";      
      }
      cout<<gWM[i]->GetErrorY(k)<<"};"<<endl;
      gWM[i]->Write();
    } 
    if(gWMS[i]){
      cout<<"Double_t gShuffledWeightedMeanInEtaAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gWMS[i]->GetN()-1;k++){
	gWMS[i]->GetPoint(k,x,y);
	cout<<y<<", ";      
      }
      gWMS[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;    
      
      cout<<"Double_t gShuffledWeightedMeanInEtaAliceError[nCentralityBins] = {";
      for(Int_t k = 0; k < gWMS[i]->GetN()-1;k++){
	cout<<gWMS[i]->GetErrorY(k)<<", ";      
      }
      cout<<gWMS[i]->GetErrorY(k)<<"};"<<endl;
      cout<<endl;
      gWMS[i]->Write();
    } 

    if(gWMP[i]){
      cout<<"Double_t gWeightedMeanInPhiAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gWMP[i]->GetN()-1;k++){
	gWMP[i]->GetPoint(k,x,y);
	cout<<y<<", ";      
      }
      gWMP[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;    
      
      cout<<"Double_t gWeightedMeanInPhiAliceError[nCentralityBins] = {";
      for(Int_t k = 0; k < gWMP[i]->GetN()-1;k++){
	cout<<gWMP[i]->GetErrorY(k)<<", ";      
      }
      cout<<gWMP[i]->GetErrorY(k)<<"};"<<endl;
      gWMP[i]->Write();
    } 
    if(gWMPS[i]){
      cout<<"Double_t gShuffledWeightedMeanInPhiAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gWMPS[i]->GetN()-1;k++){
	gWMPS[i]->GetPoint(k,x,y);
	cout<<y<<", ";      
      }
      gWMPS[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;    
      
      cout<<"Double_t gShuffledWeightedMeanInPhiAliceError[nCentralityBins] = {";
      for(Int_t k = 0; k < gWMPS[i]->GetN()-1;k++){
	cout<<gWMPS[i]->GetErrorY(k)<<", ";      
      }
      cout<<gWMPS[i]->GetErrorY(k)<<"};"<<endl;
      cout<<endl;
      gWMPS[i]->Write();
    }
    if(ginteg[i]){
      cout<<"Double_t gIntegralEtaAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < ginteg[i]->GetN()-1;k++){
	ginteg[i]->GetPoint(k,x,y);
	cout<<y<<", ";      
      }
      ginteg[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;
      cout<<endl;

      cout<<"Double_t gIntegralErrorEtaAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < ginteg[i]->GetN()-1;k++){
	cout<<ginteg[i]->GetErrorY(k)<<", ";      
      }
      cout<<ginteg[i]->GetErrorY(k)<<"};"<<endl;
      cout<<endl;

      ginteg[i]->Write();

      cout<<"Double_t gIntegralPhiAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gintegP[i]->GetN()-1;k++){
	gintegP[i]->GetPoint(k,x,y);
	cout<<y<<", ";      
      }
      gintegP[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;
      cout<<endl;

      cout<<"Double_t gIntegralErrorPhiAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gintegP[i]->GetN()-1;k++){
	cout<<gintegP[i]->GetErrorY(k)<<", ";      
      }
      cout<<gintegP[i]->GetErrorY(k)<<"};"<<endl;
      cout<<endl;

      gintegP[i]->Write();

      cout<<"Double_t gIntegralEtaShuffledAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gintegS[i]->GetN()-1;k++){
	gintegS[i]->GetPoint(k,x,y);
	cout<<y<<", ";      
      }
      gintegS[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;
      cout<<endl;

      cout<<"Double_t gIntegralErrorEtaShuffledAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gintegS[i]->GetN()-1;k++){
	cout<<gintegS[i]->GetErrorY(k)<<", ";      
      }
      cout<<gintegS[i]->GetErrorY(k)<<"};"<<endl;
      cout<<endl;

      gintegS[i]->Write();

      cout<<"Double_t gIntegralPhiShuffledAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gintegPS[i]->GetN()-1;k++){
	cout<<y<<", ";      
      }
      gintegPS[i]->GetPoint(k,x,y);
      cout<<y<<"};"<<endl;
      cout<<endl;

      cout<<"Double_t gIntegralErrorPhiShuffledAlice[nCentralityBins] = {";
      for(Int_t k = 0; k < gintegPS[i]->GetN()-1;k++){
	cout<<gintegPS[i]->GetErrorY(k)<<", ";      
      }
      cout<<gintegPS[i]->GetErrorY(k)<<"};"<<endl;
      cout<<endl;

      gintegPS[i]->Write();

    }
  }
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
  for(Int_t i = fStartBin; i < fNumberOfBins; i++)
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

//____________________________________________________________________//
void GetIntegral(TH1D *gHistBalance, Double_t &integ, Double_t &integE) {

  //always start with 1
  Int_t fStartBin = 1; 

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
  for(Int_t i = fStartBin; i < fNumberOfBins; i++)
    deltaErrorNew += gHistBalance->GetBinError(i)*(gHistBalance->GetBinCenter(i)*gSumBi - gSumBiXi)/TMath::Power(gSumBi,2);
  
  Double_t integralError = TMath::Sqrt(deltaBalP2);
  
  Double_t delta = gSumBiXi / gSumBi;
  Double_t deltaError = (gSumBiXi / gSumBi) * TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumXi2DeltaBi2)/gSumBiXi),2) + TMath::Power((gSumDeltaBi2/gSumBi),2) );
  
  cout<<"Width: "<<delta<<"\t Error: "<<deltaError<<endl;
  cout<<"New error: "<<deltaErrorNew<<endl;
  cout<<"Integral: "<<integral<<"\t Error: "<<integralError<<endl;
  cout<<"=================================================="<<endl;

  integ  = integral;
  integE = integralError;
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
    dirSubJob = dynamic_cast<TDirectoryFile *>(file->Get("PWGCFEbyE.outputBalanceFunctionAnalysis.root"));

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
  dirOutput->SetName("PWGCFEbyE.outputBalanceFunctionAnalysis.root");
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
