//=======================================================================//
//Macro to draw the main results of the charge fluctuation analysis:
//=======================================================================//

//+++++++++++++++++++++GLOBAL VARIABLES+++++++++++++++++++++//
const Int_t nCentralityBins = 20;

Double_t gNParticipants20[nCentralityBins] = {382.646,329.213,280.911,238.617,201.81,169.339,141.067,116.208,94.2515,75.4558,59.4054,45.7565,34.4839,25.3127,18.0293,12.6535,8.84139,6.16348,4.37454,3.06182};
Double_t gNParticipants20Error[nCentralityBins] = {16.9428,18.0505,17.0971,16.2076,15.5418,14.9458,14.3174,13.9067,13.2661,12.6134,11.8133,11.0495,10.0939,8.99737,7.7884,6.48725,5.21602,3.91988,2.78741,1.75066};
Double_t gNParticipants10[nCentralityBins] = {356.132,260.175,185.76,128.428,84.5666,52.3432,29.7072,15.2207,7.47192,3.71973};//10 bins
Double_t gNParticipants10Error[nCentralityBins] = {31.8228,26.9093,22.3767,18.8802,15.9735,13.3103,10.5806,7.62745,4.78538,2.41942};//10 bins

//================================ALICE================================//
Double_t gNuDynALICEData20[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNuDynALICEData20Error[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gNetChargeALICEData20[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNetChargeALICEData20Error[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gNetChargeRmsALICEData20[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNetChargeRmsALICEData20Error[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

Double_t gNuDynALICEData10[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNuDynALICEData10Error[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gNetChargeALICEData10[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNetChargeALICEData10Error[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
Double_t gNetChargeRmsALICEData10[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNetChargeRmsALICEData10Error[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//================================ALICE================================//
//+++++++++++++++++++++END OF VARIABLES+++++++++++++++++++++//

//_____________________________________________________//
void compareCentralitySets(const char*resultPath1 = "../LHC10h/Pass1_4Plus/20CentralityBins/", const char*resultPath2 = "../LHC10h/Pass1_4Plus/10CentralityBins/") {
  //Draws the nu_dyn vs centrality percentile
  SetDataPoints(resultPath1,resultPath2);
  //================================================//
  //ALICE PbPb @ 2.76 TeV
  TGraphErrors *grALICEData20Nudyn = new TGraphErrors(nCentralityBins,
						    gNParticipants20,
						    gNuDynALICEData20,
						    gNParticipants20Error,
						    gNuDynALICEData20Error);
  grALICEData20Nudyn->SetMarkerStyle(20);
  grALICEData20Nudyn->SetMarkerColor(2);

  TGraphErrors *grALICEData20NetCharge = new TGraphErrors(nCentralityBins,
							gNParticipants20,
							gNetChargeALICEData20,
							gNParticipants20Error,
							gNetChargeALICEData20Error);
  grALICEData20NetCharge->SetMarkerStyle(20);
  grALICEData20NetCharge->SetMarkerColor(2);

  TGraphErrors *grALICEData20NetChargeRms = new TGraphErrors(nCentralityBins,
							   gNParticipants20,
							   gNetChargeRmsALICEData20,
							   gNParticipants20Error,
							   gNetChargeRmsALICEData20Error);
  grALICEData20NetChargeRms->SetMarkerStyle(20);
  grALICEData20NetChargeRms->SetMarkerColor(2);

  //10 bins
  TGraphErrors *grALICEData10NetCharge = new TGraphErrors(nCentralityBins,
							gNParticipants10,
							gNetChargeALICEData10,
							gNParticipants10Error,
							gNetChargeALICEData10Error);
  grALICEData10NetCharge->SetMarkerStyle(24);
  grALICEData10NetCharge->SetMarkerColor(2);

  TGraphErrors *grALICEData10NetChargeRms = new TGraphErrors(nCentralityBins,
							   gNParticipants10,
							   gNetChargeRmsALICEData10,
							   gNParticipants10Error,
							   gNetChargeRmsALICEData10Error);
  grALICEData10NetChargeRms->SetMarkerStyle(24);
  grALICEData10NetChargeRms->SetMarkerColor(2);

  TGraphErrors *grALICEData10Nudyn = new TGraphErrors(nCentralityBins,
						    gNParticipants10,
						    gNuDynALICEData10,
						    gNParticipants10Error,
						    gNuDynALICEData10Error);
  grALICEData10Nudyn->SetMarkerStyle(24);
  grALICEData10Nudyn->SetMarkerColor(2);
  
  //_____________________________________________________//
  //Draw the results
  //_____________________________________________________//
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  //====================================//
  //Results vs centrality
  TH2F *gEmpty1 = new TH2F("gEmpty1",
			   ";N_{part.};",
			   100,0,400,10000,-0.1,50.);
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
  grALICEData20Nudyn->Draw("P");
  grALICEData10Nudyn->Draw("P");

  DrawMarker(140., -0.043, 20, 1.6, 2);
  latex->DrawLatex(150.,-0.0445,"20 centrality bins");
  DrawMarker(140., -0.05, 24, 1.6, 2);
  latex->DrawLatex(150.,-0.0515,"10 centrality bins");

  c1->SaveAs("comparissonCentralitySetsNparticipantsDependenceNuDyn.png");

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
  grALICEData20NetCharge->Draw("P");
  grALICEData10NetCharge->Draw("P");

  DrawMarker(40., 8.3, 20, 1.6, 2);
  latex->DrawLatex(50.,8.1,"20 centrality bins");
  DrawMarker(40., 7.3, 24, 1.6, 2);
  latex->DrawLatex(50.,7.1,"10 centrality bins");

  c2->SaveAs("comparissonCentralitySetsNparticipantsDependenceNetCharge.png");

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
  grALICEData20NetChargeRms->Draw("P");
  grALICEData10NetChargeRms->Draw("P");

  DrawMarker(40., 38., 20, 1.6, 2);
  latex->DrawLatex(50.,37.5,"20 centrality bins");
  DrawMarker(40., 34., 24, 1.6, 2);
  latex->DrawLatex(50.,33.5,"10 centrality bins");

  c3->SaveAs("comparissonCentralitySetsNparticipantsDependenceNetChargeRms.png");
}

//_______________________________________________________________//
void SetDataPoints(const char* resultPath1,
		   const char* resultPath2) {
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

  //20 centrality bins
  TString inputFileName;
  for(Int_t iRun = 0; iRun < nRuns; iRun++) {
    ifstream inputAscii20;
    Printf("Adding run %d",nRunNumbers[iRun]);
    inputFileName = resultPath1;
    inputFileName += "/output."; inputFileName += nRunNumbers[iRun];
    inputFileName += ".txt";

    Printf("Filename: %s",inputFileName.Data());
    inputAscii20.open(inputFileName.Data());
    for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
      inputAscii20>>gCentrality>>nuDyn[iRun][iCentrality]>>gNetCharge[iRun][iCentrality]>>gNetChargeRms[iRun][iCentrality]>>nEvents[iRun][iCentrality];
      cout<<nuDyn[iRun][iCentrality]<<" "<<nEvents[iRun][iCentrality]<<endl;
    }
    inputAscii20.close();
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

    gNetChargeALICEData20[iCentrality] = netChargeValues[iCentrality];
    gNetChargeRmsALICEData20[iCentrality] = netChargeRmsValues[iCentrality];
    gNuDynALICEData20[iCentrality] = nuDynValues[iCentrality];
  }

  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
    for(Int_t iRun = 0; iRun < nRuns; iRun++) {
      netChargeValuesError[iCentrality] += TMath::Power((netChargeValues[iCentrality]-gNetCharge[iRun][iCentrality]),2);
      netChargeRmsValuesError[iCentrality] += TMath::Power((netChargeRmsValues[iCentrality]-gNetChargeRms[iRun][iCentrality]),2);
      nuDynValuesError[iCentrality] += TMath::Power((nuDynValues[iCentrality]-nuDyn[iRun][iCentrality]),2);
    }
    gNetChargeALICEData20Error[iCentrality] = TMath::Sqrt(netChargeValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    gNetChargeRmsALICEData20Error[iCentrality] = TMath::Sqrt(netChargeRmsValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    gNuDynALICEData20Error[iCentrality] = TMath::Sqrt(nuDynValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    Printf("Centrality: %d - nu_dyn: %lf +- %lf - Net charge: %lf +- %lf - RMS: %lf - %lf",iCentrality+1,
	   gNuDynALICEData20[iCentrality],
	   gNuDynALICEData20Error[iCentrality],
	   gNetChargeALICEData20[iCentrality],
	   gNetChargeALICEData20Error[iCentrality],
	   gNetChargeRmsALICEData20[iCentrality],
	   gNetChargeRmsALICEData20Error[iCentrality]);
  }

  //10 centrality bins
  nRunCounter = 0;
  for(Int_t iRun = 0; iRun < nRuns; iRun++) {
    for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
      nuDyn[iRun][iCentrality] = 0.;
      gNetCharge[iRun][iCentrality] = 0.;
      gNetChargeRms[iRun][iCentrality] = 0.;
      nEvents[iRun][iCentrality] = 0.;
      nuDynValues[iCentrality] = 0.;
      netChargeValues[iCentrality] = 0.;
      netChargeRmsValues[iCentrality] = 0.;
      nEventsValues[iCentrality] = 0.;
      netChargeValuesError[iCentrality] = 0.;
      netChargeRmsValuesError[iCentrality] = 0.;
      nuDynValuesError[iCentrality] = 0.;
    }
  }

  for(Int_t iRun = 0; iRun < nRuns; iRun++) {
    ifstream inputAscii10;
    Printf("Adding run %d",nRunNumbers[iRun]);
    inputFileName = resultPath2;
    inputFileName += "/output."; inputFileName += nRunNumbers[iRun];
    inputFileName += ".txt";

    Printf("Filename: %s",inputFileName.Data());
    inputAscii10.open(inputFileName.Data());
    for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
      inputAscii10>>gCentrality>>nuDyn[iRun][iCentrality]>>gNetCharge[iRun][iCentrality]>>gNetChargeRms[iRun][iCentrality]>>nEvents[iRun][iCentrality];
      //cout<<nuDyn[iRun][iCentrality]<<" "<<nEvents[iRun][iCentrality]<<endl;
    }
    inputAscii10.close();
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

    gNetChargeALICEData10[iCentrality] = netChargeValues[iCentrality];
    gNetChargeRmsALICEData10[iCentrality] = netChargeRmsValues[iCentrality];
    gNuDynALICEData10[iCentrality] = nuDynValues[iCentrality];
  }

  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
    for(Int_t iRun = 0; iRun < nRuns; iRun++) {
      netChargeValuesError[iCentrality] += TMath::Power((netChargeValues[iCentrality]-gNetCharge[iRun][iCentrality]),2);
      netChargeRmsValuesError[iCentrality] += TMath::Power((netChargeRmsValues[iCentrality]-gNetChargeRms[iRun][iCentrality]),2);
      nuDynValuesError[iCentrality] += TMath::Power((nuDynValues[iCentrality]-nuDyn[iRun][iCentrality]),2);
    }
    gNetChargeALICEData10Error[iCentrality] = TMath::Sqrt(netChargeValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    gNetChargeRmsALICEData10Error[iCentrality] = TMath::Sqrt(netChargeRmsValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    gNuDynALICEData10Error[iCentrality] = TMath::Sqrt(nuDynValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    Printf("Centrality: %d - nu_dyn: %lf +- %lf - Net charge: %lf +- %lf - RMS: %lf - %lf",iCentrality+1,
	   gNuDynALICEData10[iCentrality],
	   gNuDynALICEData10Error[iCentrality],
	   gNetChargeALICEData10[iCentrality],
	   gNetChargeALICEData10Error[iCentrality],
	   gNetChargeRmsALICEData10[iCentrality],
	   gNetChargeRmsALICEData10Error[iCentrality]);
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
