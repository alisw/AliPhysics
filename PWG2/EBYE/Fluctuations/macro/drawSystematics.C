//=======================================================================//
//Macro to draw the main results of the charge fluctuation analysis:
//=======================================================================//

//+++++++++++++++++++++GLOBAL VARIABLES+++++++++++++++++++++//
const Int_t nCentralityBins = 10;
Double_t gCentralityPercentileSet1[nCentralityBins];
Double_t gCentralityPercentileSet2[nCentralityBins];
Double_t gCentralityPercentileSet3[nCentralityBins];
Double_t gCentralityPercentileError[nCentralityBins];

//================================ALICE================================//
Double_t gNuDynALICEDataSet1[nCentralityBins];
Double_t gNuDynALICEDataSet1Error[nCentralityBins];
Double_t gNuDynALICEDataSet2[nCentralityBins];
Double_t gNuDynALICEDataSet2Error[nCentralityBins];
Double_t gNuDynALICEDataSet3[nCentralityBins];
Double_t gNuDynALICEDataSet3Error[nCentralityBins];

Double_t gNetChargeALICEDataSet1[nCentralityBins];
Double_t gNetChargeALICEDataSet1Error[nCentralityBins];
Double_t gNetChargeALICEDataSet2[nCentralityBins];
Double_t gNetChargeALICEDataSet2Error[nCentralityBins];
Double_t gNetChargeALICEDataSet3[nCentralityBins];
Double_t gNetChargeALICEDataSet3Error[nCentralityBins];

Double_t gNetChargeRmsALICEDataSet1[nCentralityBins];
Double_t gNetChargeRmsALICEDataSet1Error[nCentralityBins];
Double_t gNetChargeRmsALICEDataSet2[nCentralityBins];
Double_t gNetChargeRmsALICEDataSet2Error[nCentralityBins];
Double_t gNetChargeRmsALICEDataSet3[nCentralityBins];
Double_t gNetChargeRmsALICEDataSet3Error[nCentralityBins];
//================================ALICE================================//
//+++++++++++++++++++++END OF VARIABLES+++++++++++++++++++++//

//_____________________________________________________//
void drawSystematics(const char* resultsPath) {
  //Draws the nu_dyn vs centrality percentile
  for(Int_t i = 0; i < nCentralityBins; i++) {
    if(nCentralityBins == 20) {
      gCentralityPercentileSet1[i] = 1.0 + i*5.0;
      gCentralityPercentileSet2[i] = 2.5 + i*5.0;
      gCentralityPercentileSet3[i] = 4.0 + i*5.0;
    }
    if(nCentralityBins == 10) {
      gCentralityPercentileSet1[i] = 3.0 + i*10.0;
      gCentralityPercentileSet2[i] = 5.0 + i*10.0;
      gCentralityPercentileSet3[i] = 7.0 + i*10.0;
    }
    gCentralityPercentileError[i] = 0.5;
  }

  SetDataPoints(resultsPath);

  //================================================//
  //ALICE PbPb @ 2.76 TeV
  //Nu dyn
  TGraphErrors *grNuDynALICEDataSet1 = new TGraphErrors(nCentralityBins,
							gCentralityPercentileSet1,
							gNuDynALICEDataSet1,
							gCentralityPercentileError,
							gNuDynALICEDataSet1Error);
  grNuDynALICEDataSet1->SetMarkerStyle(20);
  grNuDynALICEDataSet1->SetMarkerColor(4);
  
  TGraphErrors *grNuDynALICEDataSet2 = new TGraphErrors(nCentralityBins,
							gCentralityPercentileSet2,
							gNuDynALICEDataSet2,
							gCentralityPercentileError,
							gNuDynALICEDataSet2Error);
  grNuDynALICEDataSet2->SetMarkerStyle(25);
  grNuDynALICEDataSet2->SetMarkerColor(2);
  
  TGraphErrors *grNuDynALICEDataSet3 = new TGraphErrors(nCentralityBins,
							gCentralityPercentileSet3,
							gNuDynALICEDataSet3,
							gCentralityPercentileError,
							gNuDynALICEDataSet3Error);
  grNuDynALICEDataSet3->SetMarkerStyle(29);
  grNuDynALICEDataSet3->SetMarkerColor(3);
  
  //Net charge
  TGraphErrors *grNetChargeALICEDataSet1 = new TGraphErrors(nCentralityBins,
							    gCentralityPercentileSet1,
							    gNetChargeALICEDataSet1,
							    gCentralityPercentileError,
							    gNetChargeALICEDataSet1Error);
  grNetChargeALICEDataSet1->SetMarkerStyle(20);
  grNetChargeALICEDataSet1->SetMarkerColor(4);
  
  TGraphErrors *grNetChargeALICEDataSet2 = new TGraphErrors(nCentralityBins,
							    gCentralityPercentileSet2,
							    gNetChargeALICEDataSet2,
							    gCentralityPercentileError,
							    gNetChargeALICEDataSet2Error);
  grNetChargeALICEDataSet2->SetMarkerStyle(25);
  grNetChargeALICEDataSet2->SetMarkerColor(2);
  
  TGraphErrors *grNetChargeALICEDataSet3 = new TGraphErrors(nCentralityBins,
							    gCentralityPercentileSet3,
							    gNetChargeALICEDataSet3,
							    gCentralityPercentileError,
							    gNetChargeALICEDataSet3Error);
  grNetChargeALICEDataSet3->SetMarkerStyle(29);
  grNetChargeALICEDataSet3->SetMarkerColor(3);
  
  //Net charge rms
  TGraphErrors *grNetChargeRmsALICEDataSet1 = new TGraphErrors(nCentralityBins,
							       gCentralityPercentileSet1,
							       gNetChargeRmsALICEDataSet1,
							       gCentralityPercentileError,
							       gNetChargeRmsALICEDataSet1Error);
  grNetChargeRmsALICEDataSet1->SetMarkerStyle(20);
  grNetChargeRmsALICEDataSet1->SetMarkerColor(4);
  
  TGraphErrors *grNetChargeRmsALICEDataSet2 = new TGraphErrors(nCentralityBins,
							       gCentralityPercentileSet2,
							       gNetChargeRmsALICEDataSet2,
							       gCentralityPercentileError,
							       gNetChargeRmsALICEDataSet2Error);
  grNetChargeRmsALICEDataSet2->SetMarkerStyle(25);
  grNetChargeRmsALICEDataSet2->SetMarkerColor(2);
  
  TGraphErrors *grNetChargeRmsALICEDataSet3 = new TGraphErrors(nCentralityBins,
							       gCentralityPercentileSet3,
							       gNetChargeRmsALICEDataSet3,
							       gCentralityPercentileError,
							       gNetChargeRmsALICEDataSet3Error);
  grNetChargeRmsALICEDataSet3->SetMarkerStyle(29);
  grNetChargeRmsALICEDataSet3->SetMarkerColor(3);
  
  //_____________________________________________________//
  //Draw the results
  //_____________________________________________________//
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  //====================================//
  //results vs centrality
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
  grNuDynALICEDataSet1->Draw("P");
  grNuDynALICEDataSet2->Draw("P");
  grNuDynALICEDataSet3->Draw("P");

  DrawMarker(7., -0.053, 20, 1.6, 1);
  DrawMarker(7., -0.058, 25, 1.6, 2);
  DrawMarker(7., -0.063, 29, 1.6, 3);
  latex->DrawLatex(12.,-0.0545,"V0");
  latex->DrawLatex(12.,-0.0595,"TPC tracks");
  latex->DrawLatex(12.,-0.0645,"SPD clusters");

  c1->SaveAs("systematicsCentralityPercentileDependenceNuDyn.png");

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
  grNetChargeALICEDataSet1->Draw("P");
  grNetChargeALICEDataSet2->Draw("P");
  grNetChargeALICEDataSet3->Draw("P");

  DrawMarker(47., 9.0, 20, 1.6, 1);
  DrawMarker(47., 8.0, 25, 1.6, 2);
  DrawMarker(47., 7.0, 29, 1.6, 3);
  latex->DrawLatex(52.,8.8,"V0");
  latex->DrawLatex(52.,7.8,"TPC tracks");
  latex->DrawLatex(52.,6.8,"SPD clusters");

  c2->SaveAs("systematicsCentralityPercentileDependenceNetCharge.png");

  //============================================================//
  //net charge
  TCanvas *c3 = new TCanvas("c3","Centrality dependence: net charge rms",
			    200,200,500,500);
  c3->SetFillColor(10); c3->SetHighLightColor(10);
  c3->SetLeftMargin(0.15); c3->SetBottomMargin(0.15);
  c3->SetGridx(); c3->SetGridy();
  gEmpty1->GetYaxis()->SetRangeUser(0.0,40.0);
  gEmpty1->GetYaxis()->SetTitle("#sigma_{#Delta Q}");
  gEmpty1->DrawCopy();
  grNetChargeRmsALICEDataSet1->Draw("P");
  grNetChargeRmsALICEDataSet2->Draw("P");
  grNetChargeRmsALICEDataSet3->Draw("P");

  DrawMarker(47., 37., 20, 1.6, 1);
  DrawMarker(47., 35., 25, 1.6, 2);
  DrawMarker(47., 33., 29, 1.6, 3);
  latex->DrawLatex(52.,36.,"V0");
  latex->DrawLatex(52.,34.,"TPC tracks");
  latex->DrawLatex(52.,32,"SPD clusters");

  c3->SaveAs("systematicsCentralityPercentileDependenceNetChargeRms.png");
}

//_______________________________________________________________//
void SetDataPoints(const char* resultsPath) {
  //Calculate the mean and the statistical error of the data points

  TString inputFileName;
  Int_t gCentrality = 0;

  //Set 1
  ifstream inputAscii1;
  inputFileName = resultsPath; inputFileName += "V0/output.txt";     
  inputAscii1.open(inputFileName.Data());

  //Set 2
  ifstream inputAscii2;
  inputFileName = resultsPath; inputFileName += "TPC/output.txt";     
  inputAscii2.open(inputFileName.Data());

  //Set 3
  ifstream inputAscii3;
  inputFileName = resultsPath; inputFileName += "SPDClusters/output.txt";     
  inputAscii3.open(inputFileName.Data());

  for(Int_t iCentrality = 0; iCentrality < nCentralityBins; iCentrality++) {
    inputAscii1>>gCentrality>>gNuDynALICEDataSet1[iCentrality]>>gNuDynALICEDataSet1Error[iCentrality]>>gNetChargeALICEDataSet1[iCentrality]>>gNetChargeALICEDataSet1Error[iCentrality]>>gNetChargeRmsALICEDataSet1[iCentrality]>>gNetChargeRmsALICEDataSet1Error[iCentrality];
    inputAscii2>>gCentrality>>gNuDynALICEDataSet2[iCentrality]>>gNuDynALICEDataSet2Error[iCentrality]>>gNetChargeALICEDataSet2[iCentrality]>>gNetChargeALICEDataSet2Error[iCentrality]>>gNetChargeRmsALICEDataSet2[iCentrality]>>gNetChargeRmsALICEDataSet2Error[iCentrality];
    inputAscii3>>gCentrality>>gNuDynALICEDataSet3[iCentrality]>>gNuDynALICEDataSet3Error[iCentrality]>>gNetChargeALICEDataSet3[iCentrality]>>gNetChargeALICEDataSet3Error[iCentrality]>>gNetChargeRmsALICEDataSet3[iCentrality]>>gNetChargeRmsALICEDataSet3Error[iCentrality];

  }
  inputAscii1.close();

}


//_______________________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, 
		Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}
