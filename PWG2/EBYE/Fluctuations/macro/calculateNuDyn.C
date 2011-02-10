const Int_t nCentralityBins = 20;

void calculateNuDyn(Int_t runNumber = 1) {
  //macro to read the output of the fluctuations task
  //calculate nu_dyn for each centrality bin 
  //and draw the centrality dependence
  TString filename = "mergedOutputCentrality.Pass1_4Plus."; filename += runNumber; filename += ".root";
  TString outputFileName = "output."; outputFileName += runNumber; outputFileName += ".txt";

  //===============================================================//
  Double_t nEvents[nCentralityBins] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
  Double_t nMeanPlus[nCentralityBins] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
  Double_t nMeanMinus[nCentralityBins] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
  Double_t gNuStat[nCentralityBins] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
  Double_t gNu[nCentralityBins] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
  Double_t gNuDyn[nCentralityBins] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
  //===============================================================//

  //_______________________________________________________________//
  //Open the input file
  TFile *f = TFile::Open(filename.Data());
  if(!f->IsOpen()) {
    Printf("File not found!!!");
    break;
  }

  //_______________________________________________________________//
  //Get the TDirectoryFile
  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("outputFluctuationAnalysis.root"));
  if(!dir) {
    Printf("TDirectoryFile not found!!!");
    break;
  }

  //_______________________________________________________________//
  //Get the TList
  TList *list = dynamic_cast<TList *>(dir->Get("fluctuationsOutput"));
  if(!list) {
    Printf("TList not found!!!");
    break;
  }

  //_______________________________________________________________//
  //Get the histograms
  TH1F *gHistEventStats = dynamic_cast<TH1F *>(list->FindObject("fHistEventStats"));
  PrintEventStats(gHistEventStats);
  TH1F *gHistCentrality = dynamic_cast<TH1F *>(list->FindObject("fHistCentrality"));
  GetCentralityStats(gHistCentrality,nEvents);

  TH2F *gHistNPlusNMinus[nCentralityBins];
  TH1F *gHistNetCharge[nCentralityBins];
  TCanvas *c[nCentralityBins];
  TString histName;

  ofstream outputAscii;
  outputAscii.open(outputFileName.Data());

  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    histName = "fHistNPlusNMinusCentrality"; histName += iBin;
    gHistNPlusNMinus[iBin-1] = dynamic_cast<TH2F *>(list->FindObject(histName.Data()));
    histName = "fHistNetChargeCentrality"; histName += iBin;
    gHistNetCharge[iBin-1] = new TH1F(histName.Data(),
				      "Net Charge;#Delta Q = N_{+} - N_{-};Entries",100,-100,100);
    gHistNetCharge[iBin-1]->SetStats(kFALSE);
    //Printf("%s",gHistNPlusNMinus[iBin-1]->GetName());
    if(nEvents[iBin-1] != 0) {
      gNuDyn[iBin-1] = GetNuDun(gHistNPlusNMinus[iBin-1],
				gHistNetCharge[iBin-1],
				nEvents[iBin-1]);
      cout<<"Centrality: "<<iBin<<" - nuDyn: "<<gNuDyn[iBin-1]<<endl;
      outputAscii << iBin<<" "<<gNuDyn[iBin-1]<<" "<<gHistNetCharge[iBin-1]->GetMean()<<" "<<gHistNetCharge[iBin-1]->GetRMS()<<" "<<nEvents[iBin-1]<<endl;

      //c[iBin-1] = new TCanvas(histName.Data(),histName.Data(),
      //(iBin-1)*50,(iBin-1)*50,600,500);
      //c[iBin-1]->SetFillColor(10); c[iBin-1]->SetHighLightColor(10);
      //gHistNetCharge[iBin-1]->Draw("E");
    }//centrality with events
  }//centrality loop
  outputAscii.close();
}

//_______________________________________________________________//
void PrintEventStats(TH1F *gHistEventStats) {
  Printf("================EVENT STATISTICS================");
  Printf("Total events: %d",(Int_t)(gHistEventStats->GetBinContent(1)));  
  Printf("Total triggered events: %d",(Int_t)(gHistEventStats->GetBinContent(2)));
  Printf("Total events with proper vertex: %d",(Int_t)(gHistEventStats->GetBinContent(3)));  
  Printf("Total events analyzed: %d",(Int_t)(gHistEventStats->GetBinContent(4)));
  Printf("================EVENT STATISTICS================\n\n\n");

  return;
}

//_______________________________________________________________//
void GetCentralityStats(TH1F *gHistCentrality,
			Double_t *nEvents) {
  Printf("================CENTRALITY STATISTICS================");
  for(Int_t iBin = 1; iBin <= gHistCentrality->GetNbinsX(); iBin++) {
    nEvents[iBin-1] = gHistCentrality->GetBinContent(iBin);  
    Printf("Centrality %d: %d",iBin,
	   (Int_t)(gHistCentrality->GetBinContent(iBin)));  
  }
  Printf("================CENTRALITY STATISTICS================");

  return;
}

//_______________________________________________________________//
Double_t GetNuDun(TH2F *gHistNPlusNMinus,
		  TH1F *gHistNetCharge,
		  Int_t numberOfEvents) {
  Double_t avgNPlus = gHistNPlusNMinus->GetMean(1);
  Double_t avgNMinus = gHistNPlusNMinus->GetMean(2);
  Double_t gNuStat = (1./avgNPlus) + (1./avgNMinus);
  Double_t gNuDyn = 0., gNu = 0.;

  for(Int_t iBinX = 1; iBinX <= gHistNPlusNMinus->GetNbinsX(); iBinX++) {
    for(Int_t iBinY = 1; iBinY <= gHistNPlusNMinus->GetNbinsY(); iBinY++) {
      if(gHistNPlusNMinus->GetBinContent(iBinX,iBinY) != 0) {
	for(Int_t nEvents = 0; nEvents < gHistNPlusNMinus->GetBinContent(iBinX,iBinY); nEvents++) {
	  Int_t nPlus = (Int_t)(gHistNPlusNMinus->GetXaxis()->GetBinCenter(iBinX));
	  Int_t nMinus = (Int_t)(gHistNPlusNMinus->GetYaxis()->GetBinCenter(iBinY));
	  if((avgNPlus != 0)&&(avgNMinus != 0)) {
	    gHistNetCharge->Fill(nPlus-nMinus);
	    gNu += TMath::Power(((nPlus/avgNPlus)-(nMinus/avgNMinus)),2); 
	  }
	}
      }
    }
  }

  gNu /= numberOfEvents;
  gNuDyn = gNu - gNuStat;

  return gNuDyn;
}
