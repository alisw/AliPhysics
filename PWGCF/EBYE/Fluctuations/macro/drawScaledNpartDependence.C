//=======================================================================//
//Macro to draw the main results of the charge fluctuation analysis:
//=======================================================================//

//+++++++++++++++++++++GLOBAL VARIABLES+++++++++++++++++++++//
const Int_t nCentralityBins = 20;

Double_t gNParticipants[nCentralityBins] = {382.646,329.213,280.911,238.617,201.81,169.339,141.067,116.208,94.2515,75.4558,59.4054,45.7565,34.4839,25.3127,18.0293,12.6535,8.84139,6.16348,4.37454,3.06182};
Double_t gNParticipantsError[nCentralityBins] = {16.9428,18.0505,17.0971,16.2076,15.5418,14.9458,14.3174,13.9067,13.2661,12.6134,11.8133,11.0495,10.0939,8.99737,7.7884,6.48725,5.21602,3.91988,2.78741,1.75066};
//Double_t gNParticipants[nCentralityBins] = {356.132,260.175,185.76,128.428,84.5666,52.3432,29.7072,15.2207,7.47192,3.71973};//10 bins
//Double_t gNParticipantsError[nCentralityBins] = {31.8228,26.9093,22.3767,18.8802,15.9735,13.3103,10.5806,7.62745,4.78538,2.41942};//10 bins

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
void drawScaledNpartDependence(const char*resultPath = "../LHC10h/Pass1_4Plus/20CentralityBins/") {
  //Draws the nu_dyn vs centrality percentile
  SetDataPoints(resultPath);
  //================================================//
  //ALICE PbPb @ 2.76 TeV
  TGraphErrors *grALICEDataNudyn = new TGraphErrors(nCentralityBins,
						    gNParticipants,
						    gNuDynALICEData,
						    gNParticipantsError,
						    gNuDynALICEDataError);
  grALICEDataNudyn->SetMarkerStyle(20);
  grALICEDataNudyn->SetMarkerColor(2);

  TGraphErrors *grALICEDataNetCharge = new TGraphErrors(nCentralityBins,
							gNParticipants,
							gNetChargeALICEData,
							gNParticipantsError,
							gNetChargeALICEDataError);
  grALICEDataNetCharge->SetMarkerStyle(20);
  grALICEDataNetCharge->SetMarkerColor(2);

  TGraphErrors *grALICEDataNetChargeRms = new TGraphErrors(nCentralityBins,
							   gNParticipants,
							   gNetChargeRmsALICEData,
							   gNParticipantsError,
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
			   ";N_{part.};",
			   100,0,400,10000,-2.,1.);
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
  gEmpty1->GetYaxis()->SetRangeUser(-2.,0.0);
  gEmpty1->GetYaxis()->SetTitle("#nu_{dyn.} #times N_{part.}");
  gEmpty1->DrawCopy();
  f1->Draw("same");
  grALICEDataNudyn->Draw("P");

  DrawMarker(40., 0.0043, 20, 1.6, 2);
  latex->DrawLatex(50.,0.0033,"PbPb @ #sqrt{s_{NN}} = 2.76 TeV");

  c1->SaveAs("scaledNparticipantsDependenceNuDyn.png");

  //============================================================//
  //net charge
  TCanvas *c2 = new TCanvas("c2","Centrality dependence: net charge",
			    100,100,500,500);
  c2->SetFillColor(10); c2->SetHighLightColor(10);
  c2->SetLeftMargin(0.15); c2->SetBottomMargin(0.15);
  c2->SetGridx(); c2->SetGridy();
  gEmpty1->GetYaxis()->SetRangeUser(0.0,0.05);
  gEmpty1->GetYaxis()->SetTitle("#LT #Delta Q = N_{+} - N_{-}#GT/N_{part.}");
  gEmpty1->DrawCopy();
  grALICEDataNetCharge->Draw("P");

  DrawMarker(40., 8.3, 20, 1.6, 2);
  latex->DrawLatex(50.,8.1,"PbPb @ #sqrt{s_{NN}} = 2.76 TeV");

  c2->SaveAs("scaledNparticipantsDependenceNetCharge.png");

  //============================================================//
  //rms net charge
  TCanvas *c3 = new TCanvas("c3","Centrality dependence: net charge rms",
			    200,200,500,500);
  c3->SetFillColor(10); c3->SetHighLightColor(10);
  c3->SetLeftMargin(0.15); c3->SetBottomMargin(0.15);
  c3->SetGridx(); c3->SetGridy();
  gEmpty1->GetYaxis()->SetRangeUser(0.0,0.6);
  gEmpty1->GetYaxis()->SetTitle("#sigma_{#Delta Q}/N_{part.}");
  gEmpty1->DrawCopy();
  grALICEDataNetChargeRms->Draw("P");

  DrawMarker(40., 38., 20, 1.6, 2);
  latex->DrawLatex(50.,37.,"PbPb @ #sqrt{s_{NN}} = 2.76 TeV");

  c3->SaveAs("scaledNparticipantsDependenceNetChargeRms.png");
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

    Printf("Filename: %s",inputFileName.Data());
    inputAscii.open(inputFileName.Data());
    for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
      inputAscii>>gCentrality>>nuDyn[iRun][iCentrality]>>gNetCharge[iRun][iCentrality]>>gNetChargeRms[iRun][iCentrality]>>nEvents[iRun][iCentrality];
      cout<<nuDyn[iRun][iCentrality]<<" "<<nEvents[iRun][iCentrality]<<endl;
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

    gNetChargeALICEData[iCentrality] = netChargeValues[iCentrality]/gNParticipants[iCentrality];
    gNetChargeRmsALICEData[iCentrality] = netChargeRmsValues[iCentrality]/gNParticipants[iCentrality];;
    gNuDynALICEData[iCentrality] = nuDynValues[iCentrality]*gNParticipants[iCentrality];
  }

  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
    for(Int_t iRun = 0; iRun < nRuns; iRun++) {
      netChargeValuesError[iCentrality] += TMath::Power((netChargeValues[iCentrality]-gNetCharge[iRun][iCentrality]),2);
      netChargeRmsValuesError[iCentrality] += TMath::Power((netChargeRmsValues[iCentrality]-gNetChargeRms[iRun][iCentrality]),2);
      nuDynValuesError[iCentrality] += TMath::Power((nuDynValues[iCentrality]-nuDyn[iRun][iCentrality]),2);
    }
    gNetChargeALICEDataError[iCentrality] = TMath::Sqrt(netChargeValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)))/gNParticipants[iCentrality];
    gNetChargeRmsALICEDataError[iCentrality] = TMath::Sqrt(netChargeRmsValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)))/gNParticipants[iCentrality];
    gNuDynALICEDataError[iCentrality] = TMath::Sqrt(nuDynValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    Printf("Centrality: %d - nu_dyn: %lf +- %lf - Net charge: %lf +- %lf - RMS: %lf - %lf",iCentrality+1,
	   gNuDynALICEData[iCentrality],
	   gNuDynALICEDataError[iCentrality],
	   gNetChargeALICEData[iCentrality],
	   gNetChargeALICEDataError[iCentrality],
	   gNetChargeRmsALICEData[iCentrality],
	   gNetChargeRmsALICEDataError[iCentrality]);
  }

}


//_______________________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, 
		Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}
