//=======================================================================//
//Macro to draw the main results of the charge fluctuation analysis:
//=======================================================================//

//+++++++++++++++++++++GLOBAL VARIABLES+++++++++++++++++++++//
const Int_t nCentralityBins = 20;

//Double_t gNParticipants[nCentralityBins] = {382.646,329.213,280.911,238.617,201.81,169.339,141.067,116.208,94.2515,75.4558,59.4054,45.7565,34.4839,25.3127,18.0293,12.6535,8.84139,6.16348,4.37454,3.06182};
//Double_t gNParticipantsError[nCentralityBins] = {16.9428,18.0505,17.0971,16.2076,15.5418,14.9458,14.3174,13.9067,13.2661,12.6134,11.8133,11.0495,10.0939,8.99737,7.7884,6.48725,5.21602,3.91988,2.78741,1.75066};
Double_t gNParticipants[nCentralityBins] = {356.132,260.175,185.76,128.428,84.5666,52.3432,29.7072,15.2207,7.47192,3.71973};//10 bins
Double_t gNParticipantsError[nCentralityBins] = {31.8228,26.9093,22.3767,18.8802,15.9735,13.3103,10.5806,7.62745,4.78538,2.41942};//10 bins

//================================ALICE================================//
Double_t gNuDynALICEData[nCentralityBins] = {10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
Double_t gNuDynALICEDataError[nCentralityBins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//================================ALICE================================//

//================================ALICE================================//
Double_t gNParticipants200[nCentralityBins] = {352.,299.,235.,167.,116.,77.,47.,27.,14.};
Double_t gNParticipants200Error[nCentralityBins] = {4.,7.,9.,10.,10.,9.,8.,6.,4.};
Double_t gNuDynSTARData200[nCentralityBins] = {-0.002308,-0.002821,-0.003567,-0.004950,-0.007078,-0.010419,-0.010419,-0.016410,-0.028150,-0.056151};
Double_t gNuDynSTARData200Error[nCentralityBins] = {2e-06,2e-06,2e-06,2e-06,2e-06,2e-06,2e-06,2e-06,3e-06};

Double_t gNParticipants130[nCentralityBins] = {352.,298.,234.,166.,115.,76.,47.,27.,14.};
Double_t gNParticipants130Error[nCentralityBins] = {4.,7.,9.,10.,10.,9.,8.,5.,4.};
Double_t gNuDynSTARData130[nCentralityBins] = {-0.0023,-0.0029,-0.0036,-0.0051,-0.0070,-0.0103,-0.0103,-0.0154,-0.0274,-0.0547};
Double_t gNuDynSTARData130Error[nCentralityBins] = {0.0003,0.0002,0.0001,0.0001,0.0002,0.0002,0.0002,0.0002,0.0002};

Double_t gNParticipants62[nCentralityBins] = {347.,293.,229.,162.,112.,74.,46.,26.,13.};
Double_t gNParticipants62Error[nCentralityBins] = {4.,7.,9.,10.,10.,9.,7.,6.,3.};
Double_t gNuDynSTARData62[nCentralityBins] = {-0.002954,-0.003582,-0.004628,-0.006292,-0.008876,-0.013046,-0.020570,-0.035563,-0.071856};
Double_t gNuDynSTARData62Error[nCentralityBins] = {10e-06,8e-06,6e-06,6e-06,6e-06,6e-06,7e-06,7e-06,7e-06};
//================================ALICE================================//
//+++++++++++++++++++++END OF VARIABLES+++++++++++++++++++++//

//_____________________________________________________//
void compareRHIC(const char*resultPath = "../LHC10h/Pass1_4Plus/10CentralityBins/") {
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

  //================================================//
  //STAR AuAu @ 200 GeV
  TGraphErrors *grSTARData200Nudyn = new TGraphErrors(9,
						      gNParticipants,
						      gNuDynSTARData200,
						      gNParticipantsError,
						      gNuDynSTARData200Error);
  grSTARData200Nudyn->SetMarkerStyle(25);
  grSTARData200Nudyn->SetMarkerColor(1);

  //================================================//
  //STAR AuAu @ 130 GeV
  TGraphErrors *grSTARData130Nudyn = new TGraphErrors(9,
						      gNParticipants,
						      gNuDynSTARData130,
						      gNParticipantsError,
						      gNuDynSTARData130Error);
  grSTARData130Nudyn->SetMarkerStyle(26);
  grSTARData130Nudyn->SetMarkerColor(3);
  grSTARData130Nudyn->SetMarkerSize(1.4);

  //================================================//
  //STAR AuAu @ 62 GeV
  TGraphErrors *grSTARData62Nudyn = new TGraphErrors(9,
						     gNParticipants,
						     gNuDynSTARData62,
						     gNParticipantsError,
						     gNuDynSTARData62Error);
  grSTARData62Nudyn->SetMarkerStyle(30);
  grSTARData62Nudyn->SetMarkerColor(4);
  grSTARData62Nudyn->SetMarkerSize(1.6);
  
  //_____________________________________________________//
  //Draw the results
  //_____________________________________________________//
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  //====================================//
  //Results vs centrality
  TH2F *gEmpty1 = new TH2F("gEmpty1",
			   ";N_{part.};",
			   100,0,400,10000,-0.1,0.01);
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
  grSTARData200Nudyn->Draw("P");
  grSTARData130Nudyn->Draw("P");
  grSTARData62Nudyn->Draw("P");
  grALICEDataNudyn->Draw("P");

  DrawMarker(140., -0.043, 20, 1.6, 2);
  latex->DrawLatex(150.,-0.0445,"ALICE PbPb @ #sqrt{s_{NN}} = 2.76 TeV");
  DrawMarker(140., -0.05, 25, 1.6, 1);
  latex->DrawLatex(150.,-0.0515,"STAR AuAu @ #sqrt{s_{NN}} = 200 GeV");
  DrawMarker(140., -0.057, 26, 1.6, 3);
  latex->DrawLatex(150.,-0.0585,"STAR AuAu @ #sqrt{s_{NN}} = 200 GeV");
  DrawMarker(140., -0.064, 30, 1.6, 4);
  latex->DrawLatex(150.,-0.0655,"STAR AuAu @ #sqrt{s_{NN}} = 200 GeV");

  c1->SaveAs("comparisonRHICNparticipantsDependenceNuDyn.png");
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

    gNuDynALICEData[iCentrality] = nuDynValues[iCentrality];
  }

  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
    for(Int_t iRun = 0; iRun < nRuns; iRun++) {
      netChargeValuesError[iCentrality] += TMath::Power((netChargeValues[iCentrality]-gNetCharge[iRun][iCentrality]),2);
      netChargeRmsValuesError[iCentrality] += TMath::Power((netChargeRmsValues[iCentrality]-gNetChargeRms[iRun][iCentrality]),2);
      nuDynValuesError[iCentrality] += TMath::Power((nuDynValues[iCentrality]-nuDyn[iRun][iCentrality]),2);
    }
    gNuDynALICEDataError[iCentrality] = TMath::Sqrt(nuDynValuesError[iCentrality]/(nRunCounter*(nRunCounter-1)));
    Printf("Centrality: %d - nu_dyn: %lf +- %lf ",iCentrality+1,
	   gNuDynALICEData[iCentrality],
	   gNuDynALICEDataError[iCentrality]);
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
