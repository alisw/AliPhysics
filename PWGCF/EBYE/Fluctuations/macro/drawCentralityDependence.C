//=======================================================================//
//Macro to draw the main results of the charge fluctuation analysis:
//=======================================================================//

//+++++++++++++++++++++GLOBAL VARIABLES+++++++++++++++++++++//
const Int_t nCentralityBins = 20;
//Double_t gCentralityPercentile[nCentralityBins] = {2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5,82.5,87.5,92.5,97.5};
//Double_t gCentralityPercentileError[nCentralityBins] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
Double_t gCentralityPercentile[nCentralityBins] = {5.,15.,25.,35.,45.,55.,65.,75.,85.,95.};
Double_t gCentralityPercentileError[nCentralityBins] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};

//================================ALICE================================//
Double_t gNuDynALICEData[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNuDynALICEDataError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gNetChargeALICEData[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNetChargeALICEDataError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gNetChargeRmsALICEData[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNetChargeRmsALICEDataError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//================================ALICE================================//
//+++++++++++++++++++++END OF VARIABLES+++++++++++++++++++++//

//_____________________________________________________//
void drawCentralityDependence(const char*resultPath = "../LHC10h/Pass1_4Plus/20CentralityBins/") {
  //Draws the nu_dyn vs centrality percentile
  SetDataPoints(resultPath);
  //================================================//
  //ALICE PbPb @ 2.76 TeV
  TGraphErrors *grALICEDataNudyn = new TGraphErrors(nCentralityBins,
						    gCentralityPercentile,
						    gNuDynALICEData,
						    gCentralityPercentileError,
						    gNuDynALICEDataError);
  grALICEDataNudyn->SetMarkerStyle(20);
  grALICEDataNudyn->SetMarkerColor(2);

  TGraphErrors *grALICEDataNetCharge = new TGraphErrors(nCentralityBins,
							gCentralityPercentile,
							gNetChargeALICEData,
							gCentralityPercentileError,
							gNetChargeALICEDataError);
  grALICEDataNetCharge->SetMarkerStyle(20);
  grALICEDataNetCharge->SetMarkerColor(2);

  TGraphErrors *grALICEDataNetChargeRms = new TGraphErrors(nCentralityBins,
							   gCentralityPercentile,
							   gNetChargeRmsALICEData,
							   gCentralityPercentileError,
							   gNetChargeRmsALICEDataError);
  grALICEDataNetChargeRms->SetMarkerStyle(20);
  grALICEDataNetChargeRms->SetMarkerColor(2);
  
  //_____________________________________________________//
  //Draw the results
  //_____________________________________________________//
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  //====================================//
  //Results vs centrality
  TH2F *gEmpty1 = new TH2F("gEmpty1",
			   ";Centrality percentile;",
			   nCentralityBins,0,100,10000,-0.1,50.);
  gEmpty1->SetStats(kFALSE);
  gEmpty1->GetYaxis()->SetTitleOffset(1.5);
  gEmpty1->GetXaxis()->SetTitleOffset(1.5);
  gEmpty1->GetYaxis()->SetNdivisions(10);
  gEmpty1->GetXaxis()->SetNdivisions(10);

  TF1 *f1 = new TF1("f1","0",0,1000);
  f1->SetLineColor(1); f1->SetLineStyle(3); f1->SetLineWidth(2);

  //============================================================//
  //nu_{dyn.}
  TCanvas *c1 = new TCanvas("c1","Centrality dependence: nu_dyn",
			    0,0,500,500);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.15);
  c1->SetGridx(); c1->SetGridy();
  gEmpty1->GetYaxis()->SetRangeUser(-0.08,0.01);
  gEmpty1->GetYaxis()->SetTitle("#nu_{dyn.}");
  gEmpty1->DrawCopy();
  f1->Draw("same");
  grALICEDataNudyn->Draw("P");

  DrawMarker(47., 0.0043, 20, 1.6, 2);
  latex->DrawLatex(52.,0.0035,"PbPb @ #sqrt{s_{NN}} = 2.76 TeV");

  c1->SaveAs("centralityPercentileDependenceNuDyn.png");

  //============================================================//
  //net charge
  TCanvas *c2 = new TCanvas("c2","Centrality dependence: net charge",
			    100,100,500,500);
  c2->SetFillColor(10); c2->SetHighLightColor(10);
  c2->SetLeftMargin(0.15); c2->SetBottomMargin(0.15);
  c2->SetGridx(); c2->SetGridy();
  gEmpty1->GetYaxis()->SetRangeUser(0.0,10.0);
  gEmpty1->GetYaxis()->SetTitle("#LT #Delta Q = N_{+} - N_{-}#GT");
  gEmpty1->DrawCopy();
  grALICEDataNetCharge->Draw("P");

  DrawMarker(47., 8.3, 20, 1.6, 2);
  latex->DrawLatex(52.,8.1,"PbPb @ #sqrt{s_{NN}} = 2.76 TeV");

  c2->SaveAs("centralityPercentileDependenceNetCharge.png");

  //============================================================//
  //rms net charge
  TCanvas *c3 = new TCanvas("c3","Centrality dependence: net charge rms",
			    200,200,500,500);
  c3->SetFillColor(10); c3->SetHighLightColor(10);
  c3->SetLeftMargin(0.15); c3->SetBottomMargin(0.15);
  c3->SetGridx(); c3->SetGridy();
  gEmpty1->GetYaxis()->SetRangeUser(0.0,40.0);
  gEmpty1->GetYaxis()->SetTitle("#sigma_{#Delta Q}");
  gEmpty1->DrawCopy();
  grALICEDataNetChargeRms->Draw("P");

  DrawMarker(47., 38.3, 20, 1.6, 2);
  latex->DrawLatex(52.,38.1,"PbPb @ #sqrt{s_{NN}} = 2.76 TeV");

  c3->SaveAs("centralityPercentileDependenceNetChargeRms.png");
}

//_______________________________________________________________//
void SetDataPoints(const char* resultPath) {
  //Calculate the mean and the statistical error of the data points
  const Int_t nRuns = 7;
  Double_t nuDyn[nRuns][nCentralityBins];
  Double_t gNetCharge[nRuns][nCentralityBins];
  Double_t gNetChargeRms[nRuns][nCentralityBins];
  Double_t nuDyn[nRuns][nCentralityBins];
  Int_t nEvents[nRuns][nCentralityBins];
  Int_t nRunNumbers[nRuns] = {137161,137431,137549,
			      137595,137638,137639,137693};
  Int_t gCentrality;
  Double_t netChargeValues[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
					   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t netChargeValuesError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,
						    0.,0.,0.,0.,0.,0.,0.,0.,
						    0.,0.,0.,0.};
  Double_t netChargeRmsValues[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,
						  0.,0.,0.,0.,0.,0.,0.,
						  0.,0.,0.,0.,0.,0.};
  Double_t netChargeRmsValuesError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,
						       0.,0.,0.,0.,0.,0.,0.,0.,
						    0.,0.,0.,0.};
  Double_t nuDynValues[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
					   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t nuDynValuesError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
						0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t nEventsValues[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
					     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  TString inputFileName;
  for(Int_t iRun = 0; iRun < nRuns; iRun++) {
    ifstream inputAscii;
    Printf("Adding run %d",nRunNumbers[iRun]);
    inputFileName = resultPath;
    inputFileName += "/output."; inputFileName += nRunNumbers[iRun];
    inputFileName += ".txt";

    //Printf("Filename: %s",inputFileName.Data());
    inputAscii.open(inputFileName.Data());
    for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
      inputAscii>>gCentrality>>nuDyn[iRun][iCentrality]>>gNetCharge[iRun][iCentrality]>>gNetChargeRms[iRun][iCentrality]>>nEvents[iRun][iCentrality];
      //cout<<nuDyn[iRun][iCentrality]<<" "<<nEvents[iRun][iCentrality]<<endl;
    }
    inputAscii.close();
  }

  Int_t nRunCounter = 0;
  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
    for(Int_t iRun = 0; iRun < nRuns; iRun++) {
      nuDynValues[iCentrality] += nuDyn[iRun][iCentrality]*nEvents[iRun][iCentrality];
      netChargeValues[iCentrality] += gNetCharge[iRun][iCentrality]*nEvents[iRun][iCentrality];
      netChargeRmsValues[iCentrality] += gNetChargeRms[iRun][iCentrality]*nEvents[iRun][iCentrality];
      nEventsValues[iCentrality] += nEvents[iRun][iCentrality];
    }

    if(nEventsValues[iCentrality] != 0) {
      nRunCounter += 1;

      netChargeValues[iCentrality] /= nEventsValues[iCentrality];
      netChargeRmsValues[iCentrality] /= nEventsValues[iCentrality];
      nuDynValues[iCentrality] /= nEventsValues[iCentrality];
    }
    else {
      netChargeValues[iCentrality] = 999.;
      netChargeRmsValues[iCentrality] = 999.;
      nuDynValues[iCentrality] = 999.;
    }

    gNetChargeALICEData[iCentrality] = netChargeValues[iCentrality];
    gNetChargeRmsALICEData[iCentrality] = netChargeRmsValues[iCentrality];
    gNuDynALICEData[iCentrality] = nuDynValues[iCentrality];
  }

  ofstream outputAscii;
  outputAscii.open("output.txt");
  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
    for(Int_t iRun = 0; iRun < nRuns; iRun++) {
      netChargeValuesError[iCentrality] += TMath::Power((netChargeValues[iCentrality]-gNetCharge[iRun][iCentrality]),2);
      netChargeRmsValuesError[iCentrality] += TMath::Power((netChargeRmsValues[iCentrality]-gNetChargeRms[iRun][iCentrality]),2);
      nuDynValuesError[iCentrality] += TMath::Power((nuDynValues[iCentrality]-nuDyn[iRun][iCentrality]),2);
    }
    gNetChargeALICEDataError[iCentrality] = TMath::Sqrt(netChargeValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    gNetChargeRmsALICEDataError[iCentrality] = TMath::Sqrt(netChargeRmsValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    gNuDynALICEDataError[iCentrality] = TMath::Sqrt(nuDynValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    Printf("Centrality: %d - nu_dyn: %lf +- %lf - Net charge: %lf +- %lf - RMS: %lf - %lf",iCentrality+1,
	   gNuDynALICEData[iCentrality],
	   gNuDynALICEDataError[iCentrality],
	   gNetChargeALICEData[iCentrality],
	   gNetChargeALICEDataError[iCentrality],
	   gNetChargeRmsALICEData[iCentrality],
	   gNetChargeRmsALICEDataError[iCentrality]);

    outputAscii<<iCentrality+1<<" "<<gNuDynALICEData[iCentrality]<<" "<<gNuDynALICEDataError[iCentrality]<<" "<<gNetChargeALICEData[iCentrality]<<" "<<gNetChargeALICEDataError[iCentrality]<<" "<<gNetChargeRmsALICEData[iCentrality]<<" "<<gNetChargeRmsALICEDataError[iCentrality]<<endl;
  }
  outputAscii.close();

}


//_______________________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, 
		Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}
