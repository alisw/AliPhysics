void drawEnergyDependence() {
  //Macro to draw the energy dependence of the nu_dyn

  //================================ALICE================================//
  Double_t gNuDynALICEData[1] = {-0.000991451};
  Double_t gNuDynALICEDataStatError[1] = {3.06304e-05};

  Double_t gALICEEnergy[1] = {2760.};
  Double_t gALICEEnergyError[1] = {10.};

  TGraphErrors *grNuDynALICEData = new TGraphErrors(1,gALICEEnergy,
						    gNuDynALICEData,
						    gALICEEnergyError,
						    gNuDynALICEDataStatError);
  grNuDynALICEData->SetMarkerStyle(20);
  grNuDynALICEData->SetMarkerColor(2);
  grNuDynALICEData->SetMarkerSize(1.4);

  //================================CERES================================//
  Double_t gNuDynCERESData[3] = {-0.00006,-0.00069,-0.00108};
  Double_t gNuDynCERESDataStatError[3] = {5.0e-05,5.0e-05,5.0e-05};

  Double_t gCERESEnergy[3] = {7.8,12.3,17.3};
  Double_t gCERESEnergyError[3] = {0.5,0.5,0.5};

  TGraphErrors *grNuDynCERESData = new TGraphErrors(3,gCERESEnergy,
						    gNuDynCERESData,
						    gCERESEnergyError,
						    gNuDynCERESDataStatError);
  grNuDynCERESData->SetMarkerStyle(23);
  grNuDynCERESData->SetMarkerColor(1);
  grNuDynCERESData->SetMarkerSize(1.4);

  //================================STAR================================//
  Double_t gNuDynSTARData[4] = {-0.00113,-0.00147,-0.00122,-0.00163};
  Double_t gNuDynSTARDataStatError[4] = {0.00026,0.00018,0.00014,7e-05};

  Double_t gSTAREnergy[4] = {19.6,62.4,130.,200.};
  Double_t gSTAREnergyError[4] = {0.5,0.5,0.5,0.5};

  TGraphErrors *grNuDynSTARData = new TGraphErrors(4,gSTAREnergy,
						   gNuDynSTARData,
						   gSTAREnergyError,
						   gNuDynSTARDataStatError);
  grNuDynSTARData->SetMarkerStyle(29);
  grNuDynSTARData->SetMarkerColor(3);
  grNuDynSTARData->SetMarkerSize(1.6);

  //================================PHENIX================================//
  Double_t gNuDynPHENIXData[1] = {-0.00086};
  Double_t gNuDynPHENIXDataStatError[1] = {0.00038};

  Double_t gPHENIXEnergy[1] = {150.};
  Double_t gPHENIXEnergyError[4] = {0.5};

  TGraphErrors *grNuDynPHENIXData = new TGraphErrors(1,gPHENIXEnergy,
						   gNuDynPHENIXData,
						   gPHENIXEnergyError,
						   gNuDynPHENIXDataStatError);
  grNuDynPHENIXData->SetMarkerStyle(21);
  grNuDynPHENIXData->SetMarkerColor(4);
  grNuDynPHENIXData->SetMarkerSize(1.2);

  //_____________________________________________________//
  //Draw the results
  //_____________________________________________________//
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);
  //====================================//
  //results vs centrality
  TH2F *gEmpty1 = new TH2F("gEmpty1",
			   ";#sqrt{s_{NN}} (GeV);",
			   100,1,10000,10000,-0.003,0.5e-03);
  gEmpty1->SetStats(kFALSE);
  gEmpty1->GetYaxis()->SetTitleOffset(1.5);
  gEmpty1->GetXaxis()->SetTitleOffset(1.5);
  gEmpty1->GetYaxis()->SetNdivisions(10);
  gEmpty1->GetXaxis()->SetNdivisions(10);

  TF1 *f1 = new TF1("f1","0",0,10000);
  f1->SetLineColor(1); f1->SetLineStyle(3); f1->SetLineWidth(3);

  //============================================================//
  //nu_{dyn.}
  TCanvas *c1 = new TCanvas("c1","Energy dependence: nu_dyn",
			    0,0,500,500);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.15);
  c1->SetGridx(); c1->SetGridy(); c1->SetLogx();
  gEmpty1->GetYaxis()->SetTitle("#nu_{dyn.}");
  gEmpty1->DrawCopy();
  f1->Draw("same");
  grNuDynALICEData->Draw("P");
  grNuDynCERESData->Draw("P");
  grNuDynSTARData->Draw("P");
  grNuDynPHENIXData->Draw("P");

  DrawMarker(7., -0.0019, 20, 1.6, 2);
  DrawMarker(7., -0.0021, 29, 1.6, 3);
  DrawMarker(7., -0.0023, 21, 1.6, 4);
  DrawMarker(7., -0.0025, 23, 1.6, 1);
  latex->DrawLatex(10.,-0.00195,"ALICE");
  latex->DrawLatex(10.,-0.00215,"STAR");
  latex->DrawLatex(10.,-0.00235,"PHENIX");
  latex->DrawLatex(10.,-0.00255,"CERES");

  c1->SaveAs("energyDependenceNuDyn.png");
}

//_______________________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, 
		Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}
