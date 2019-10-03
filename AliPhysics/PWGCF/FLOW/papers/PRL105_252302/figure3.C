// Objects holding results are following TGraphErrors:
//  1.) Cumulant2ndTPCrefMultTPConlyAll = v2{2}, TPC ref mult, TPConly tracks, all charges
//  2.) Cumulant2ndTPCrefMultTPConlySame = v2{2}, TPC ref mult, TPConly tracks, same charges
//  3.) Cumulant4thTPCrefMultTPConlyAll = v2{4}, TPC ref mult, TPConly tracks, all charges
//  4.) Cumulant4thTPCrefMultTPConlySame = v2{4}, TPC ref mult, TPConly tracks, same charges
//  5.) FQDTPCrefMultTPConlyAll = v2{FQD}, TPC ref mult, TPConly tracks, all charges
//  6.) LYZSTPCrefMultTPConlyAll = v2{LYZ} (Lee-Yang zero, sum generating function), TPC ref mult, TPConly tracks, all charges
//  7.) EventPlaneSTAR - published STAR results for event plane method
//  8.) Cumulant4thSTAR - published STAR results for v2{4} (generating function formalism)
//  9.) LYZPSTAR - published STAR results for v2{LYZ} (Lee-Yang zero, product generating function)

void figure3()
{
 // Set style:
 SetFlowStyle();
 
 //=================================================================================================================================================
 // v2{2}, TPC ref mult, TPConly tracks, all charges
 Double_t xCumulant2ndTPCrefMultTPConlyAll[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
 Double_t yCumulant2ndTPCrefMultTPConlyAll[] = {0.028811,0.044020,0.064765,0.084340,0.096771,0.104257,0.105902,0.104897,0.104811};
 Double_t xErrCumulant2ndTPCrefMultTPConlyAll[9] = {0.};
 Double_t yErrCumulant2ndTPCrefMultTPConlyAll[] = {0.000470,0.000510,0.000433,0.000528,0.000621,0.000772,0.000955,0.001254,0.002317};
 Int_t nPointsCumulant2ndTPCrefMultTPConlyAll = sizeof(xCumulant2ndTPCrefMultTPConlyAll)/sizeof(Double_t);         
 TGraphErrors *Cumulant2ndTPCrefMultTPConlyAll = new TGraphErrors(nPointsCumulant2ndTPCrefMultTPConlyAll,xCumulant2ndTPCrefMultTPConlyAll,
                                                              yCumulant2ndTPCrefMultTPConlyAll,xErrCumulant2ndTPCrefMultTPConlyAll,
                                                              yErrCumulant2ndTPCrefMultTPConlyAll );
 Cumulant2ndTPCrefMultTPConlyAll->SetMarkerStyle(kFullCircle);
 Cumulant2ndTPCrefMultTPConlyAll->SetMarkerColor(kBlue);
 ShiftAlongXaxis(Cumulant2ndTPCrefMultTPConlyAll,0.0);
 
 // v2{2}, TPC ref mult, TPC only tracks, same charges
 Double_t xCumulant2ndTPCrefMultTPConlySame[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
 Double_t yCumulant2ndTPCrefMultTPConlySame[] = {0.027453,0.042691,0.063489,0.083387,0.095468,0.102824,0.104025,0.101708,0.100819};
 Double_t xErrCumulant2ndTPCrefMultTPConlySame[9] = {0.};
 Double_t yErrCumulant2ndTPCrefMultTPConlySame[] = {0.000463,0.000454,0.000369,0.000447,0.000541,0.000689,0.000916,0.001331,0.003072};
 Int_t nPointsCumulant2ndTPCrefMultTPConlySame = sizeof(xCumulant2ndTPCrefMultTPConlySame)/sizeof(Double_t);         
 TGraphErrors *Cumulant2ndTPCrefMultTPConlySame = new TGraphErrors(nPointsCumulant2ndTPCrefMultTPConlySame,xCumulant2ndTPCrefMultTPConlySame,
                                                              yCumulant2ndTPCrefMultTPConlySame,xErrCumulant2ndTPCrefMultTPConlySame,
                                                              yErrCumulant2ndTPCrefMultTPConlySame );
 Cumulant2ndTPCrefMultTPConlySame->SetMarkerStyle(kOpenCircle);
 Cumulant2ndTPCrefMultTPConlySame->SetMarkerColor(kBlue);
 ShiftAlongXaxis(Cumulant2ndTPCrefMultTPConlySame,-1.44);
 //=================================================================================================================================================

 //=================================================================================================================================================
 // v2{4}, TPC ref mult, TPConly tracks, all charges
 Double_t xCumulant4thTPCrefMultTPConlyAll[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
 Double_t yCumulant4thTPCrefMultTPConlyAll[] = {0.017855,0.032440,0.055818,0.073137,0.083898,0.086690,0.082040,0.077777,0.000000};
 Double_t xErrCumulant4thTPCrefMultTPConlyAll[9] = {0.};
 Double_t yErrCumulant4thTPCrefMultTPConlyAll[] = {0.004327,0.001690,0.001049,0.001150,0.001177,0.001785,0.002793,0.005153,0.000000};
 Int_t nPointsCumulant4thTPCrefMultTPConlyAll = sizeof(xCumulant4thTPCrefMultTPConlyAll)/sizeof(Double_t);         
 TGraphErrors *Cumulant4thTPCrefMultTPConlyAll = new TGraphErrors(nPointsCumulant4thTPCrefMultTPConlyAll,xCumulant4thTPCrefMultTPConlyAll,
                                                              yCumulant4thTPCrefMultTPConlyAll,xErrCumulant4thTPCrefMultTPConlyAll,
                                                              yErrCumulant4thTPCrefMultTPConlyAll );
 Cumulant4thTPCrefMultTPConlyAll->SetMarkerStyle(kFullSquare);
 Cumulant4thTPCrefMultTPConlyAll->SetMarkerColor(kRed);
 ShiftAlongXaxis(Cumulant4thTPCrefMultTPConlyAll,0.1);
 
 // v2{4}, TPC ref mult, TPC only tracks, same charges
 Double_t xCumulant4thTPCrefMultTPConlySame[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
 Double_t yCumulant4thTPCrefMultTPConlySame[] = {0.015724,0.031731,0.055660,0.072992,0.083332,0.087978,0.080705,0.077744,0.000000};
 Double_t xErrCumulant4thTPCrefMultTPConlySame[9] = {0.};
 Double_t yErrCumulant4thTPCrefMultTPConlySame[] = {0.006343,0.001764,0.000895,0.001011,0.001108,0.001685,0.003890,0.007947,0.000000};
 Int_t nPointsCumulant4thTPCrefMultTPConlySame = sizeof(xCumulant4thTPCrefMultTPConlySame)/sizeof(Double_t);         
 TGraphErrors *Cumulant4thTPCrefMultTPConlySame = new TGraphErrors(nPointsCumulant4thTPCrefMultTPConlySame,xCumulant4thTPCrefMultTPConlySame,
                                                              yCumulant4thTPCrefMultTPConlySame,xErrCumulant4thTPCrefMultTPConlySame,
                                                              yErrCumulant4thTPCrefMultTPConlySame );
 Cumulant4thTPCrefMultTPConlySame->SetMarkerStyle(kOpenSquare);
 Cumulant4thTPCrefMultTPConlySame->SetMarkerColor(kRed);
 ShiftAlongXaxis(Cumulant4thTPCrefMultTPConlySame,-1.44);
 //=================================================================================================================================================
 
 //=================================================================================================================================================
 // v2{FQD}, TPC ref mult, TPConly tracks, all charges
 Double_t xFQDTPCrefMultTPConlyAll[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
 Double_t yFQDTPCrefMultTPConlyAll[] = {0.018442,0.031239,0.054123,0.072266,0.081993,0.084699,0.080743,0.079986,0.069490};
 Double_t xErrFQDTPCrefMultTPConlyAll[9] = {0.};
 Double_t yErrFQDTPCrefMultTPConlyAll[] = {0.002872,0.001597,0.000636,0.000679,0.000887,0.001285,0.002367,0.003960,0.014574};
 Int_t nPointsFQDTPCrefMultTPConlyAll = sizeof(xFQDTPCrefMultTPConlyAll)/sizeof(Double_t);         
 TGraphErrors *FQDTPCrefMultTPConlyAll = new TGraphErrors(nPointsFQDTPCrefMultTPConlyAll,xFQDTPCrefMultTPConlyAll,
                                                              yFQDTPCrefMultTPConlyAll,xErrFQDTPCrefMultTPConlyAll,
                                                              yErrFQDTPCrefMultTPConlyAll );
 FQDTPCrefMultTPConlyAll->SetMarkerStyle(kOpenCross);
 FQDTPCrefMultTPConlyAll->SetMarkerColor(kBlack);
 ShiftAlongXaxis(FQDTPCrefMultTPConlyAll,1.25);
 //=================================================================================================================================================  

 //=================================================================================================================================================
 // v2{LYZ} (Lee-Yang zero, sum generating function), TPC ref mult, TPConly tracks, all charges
 Double_t xLYZSTPCrefMultTPConlyAll[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
 Double_t yLYZSTPCrefMultTPConlyAll[] = {0.000000,0.027581,0.053395,0.072002,0.082107,0.085281,0.081629,0.069773,0.086976};
 Double_t xErrLYZSTPCrefMultTPConlyAll[9] = {0.};
 Double_t yErrLYZSTPCrefMultTPConlyAll[] = {0.000000,0.004444,0.000708,0.000735,0.000898,0.001296,0.002460,0.019011,0.016716};
 Int_t nPointsLYZSTPCrefMultTPConlyAll = sizeof(xLYZSTPCrefMultTPConlyAll)/sizeof(Double_t);         
 TGraphErrors *LYZSTPCrefMultTPConlyAll = new TGraphErrors(nPointsLYZSTPCrefMultTPConlyAll,xLYZSTPCrefMultTPConlyAll,
                                                           yLYZSTPCrefMultTPConlyAll,xErrLYZSTPCrefMultTPConlyAll,
                                                           yErrLYZSTPCrefMultTPConlyAll );
 LYZSTPCrefMultTPConlyAll->SetMarkerStyle(kFullTriangleUp);
 LYZSTPCrefMultTPConlyAll->SetMarkerColor(kGreen+2);
 ShiftAlongXaxis(LYZSTPCrefMultTPConlyAll,2.25);
 //=================================================================================================================================================  
 
 //=================================================================================================================================================
 // Event Plane STAR - published in Phys. Rev. C 77 (2008) 54901, Figure 5a:
 Double_t xEventPlaneSTAR[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
 Double_t yEventPlaneSTAR[] = {0.0232,0.0339,0.0476,0.0618,0.0703,0.074,0.0744,0.0723,0.0696}; 
 Double_t xErrEventPlaneSTAR[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
 Double_t yErrEventPlaneSTAR[] = {0.00017,0.00014,0.0001,0.00011,0.00014,0.00019,0.00028,0.00049,0.00089};
 
 Int_t nPointsEventPlaneSTAR = sizeof(xEventPlaneSTAR)/sizeof(Double_t);                                     
 TGraphErrors *EventPlaneSTAR = new TGraphErrors(nPointsEventPlaneSTAR,xEventPlaneSTAR,
                                                 yEventPlaneSTAR,xErrEventPlaneSTAR,yErrEventPlaneSTAR);
 //EventPlaneSTAR->SetFillStyle(3001);
 //EventPlaneSTAR->SetFillColor(kBlack);
 EventPlaneSTAR->SetFillStyle(1001);
 EventPlaneSTAR->SetLineColor(kBlue);
 EventPlaneSTAR->SetFillColor(kBlue-10);
 //=================================================================================================================================================

 //=================================================================================================================================================
 // 4th order cumulant STAR - published pt-integrated v2{4} in Phys. Rev. C 77 (2008) 54901, Figure 5):
 Double_t xCumulant4thSTAR[] = {7.5,15.,25.,35.,45.,55.,65.};
 Double_t yCumulant4thSTAR[] = {0.0269707,0.0425772,0.0563609,0.0630763,0.0655523,0.0617547,0.0581859};
 Double_t xErrCumulant4thSTAR[] = {0.,0.,0.,0.,0.,0.,0.};
 Double_t yErrCumulant4thSTAR[] = {0.000241449,0.0000805643,0.0000674778,0.0000873136,0.000144613,0.000318834,0.000947424};
 Int_t nPointsCumulant4thSTAR = sizeof(xCumulant4thSTAR)/sizeof(Double_t);                                     
 TGraphErrors *Cumulant4thSTAR = new TGraphErrors(nPointsCumulant4thSTAR,xCumulant4thSTAR,
                                                          yCumulant4thSTAR,xErrCumulant4thSTAR,yErrCumulant4thSTAR);
 Cumulant4thSTAR->SetFillStyle(1001);
 Cumulant4thSTAR->SetFillColor(kRed-10);
 //=================================================================================================================================================
 
 //=================================================================================================================================================
 // LYZ product STAR - published pt-integrated LYZ product result in Phys. Rev. C 77 (2008) 54901, Figure 5a:
 Double_t xLYZPSTAR[] = {7.5,15.,25.,35.,45.,55.};
 Double_t yLYZPSTAR[] = {0.0292,0.0404,0.0541,0.0612,0.0627,0.063};
 Double_t xErrLYZPSTAR[] = {0.,0.,0.,0.,0.};
 Double_t yErrLYZPSTAR[] = {0.0023,0.00032,0.00024,0.0003,0.0017,0.0037};
 Int_t nPointsLYZPSTAR = sizeof(xLYZPSTAR)/sizeof(Double_t);                                     
 TGraphErrors *LYZPSTAR = new TGraphErrors(nPointsLYZPSTAR,xLYZPSTAR,
                                                   yLYZPSTAR,xErrLYZPSTAR,yErrLYZPSTAR);
 LYZPSTAR->SetMarkerStyle(kFullTriangleUp);
 LYZPSTAR->SetMarkerColor(kBlue); 
 //LYZPSTAR->SetFillStyle(3003);
 //LYZPSTAR->SetFillColor(kBlack);
 LYZPSTAR->SetFillStyle(1001);
 LYZPSTAR->SetLineColor(kRed);
 LYZPSTAR->SetFillColor(kRed-10);
 ShiftAlongXaxis(LYZPSTAR,2.0);
 //=================================================================================================================================================

 
 // Canvas:
 TCanvas *c1 = new TCanvas("c1","v_{2} reference flow"); 
 c1->cd(0);
 // Legend:
 TLegend *l1 = new TLegend(0.45,0.18,0.75,0.50);
 l1->SetFillStyle(0); // white legend background
 l1->AddEntry(Cumulant2ndTPCrefMultTPConlyAll,"v_{2}{2}","p"); 
 l1->AddEntry(Cumulant2ndTPCrefMultTPConlySame,"v_{2}{2} (same charge)","p"); 
 l1->AddEntry(Cumulant4thTPCrefMultTPConlyAll,"v_{2}{4}","p"); 
 l1->AddEntry(Cumulant4thTPCrefMultTPConlySame,"v_{2}{4} (same charge)","p"); 
 l1->AddEntry(FQDTPCrefMultTPConlyAll,"v_{2}{q-dist}","p"); 
 l1->AddEntry(LYZSTPCrefMultTPConlyAll,"v_{2}{LYZ}","p");  
 l1->AddEntry(EventPlaneSTAR,"v_{2}{EP} STAR","l"); 
 l1->AddEntry(LYZPSTAR,"v_{2}{LYZ} STAR","l"); 
 //l1->AddEntry(hydro,"hydro #eta/s=0.08","l"); 
 l1->SetTextSize(0.038);
 l1->SetTextFont(22); // 22 = Times New Roman (bold)
  
 // Final drawing (order is important!):
 StyleHistogram()->Draw(); 
 EventPlaneSTAR->Draw("same3"); // mash
 EventPlaneSTAR->Draw("lsameX"); // line
 LYZPSTAR->Draw("same3"); // mash
 LYZPSTAR->Draw("lsameX"); // line
 //Cumulant4thSTAR->Draw("same3");  
 Cumulant2ndTPCrefMultTPConlyAll->Draw("psame");
 Cumulant2ndTPCrefMultTPConlySame->Draw("psame");
 Cumulant4thTPCrefMultTPConlyAll->Draw("psame");
 Cumulant4thTPCrefMultTPConlySame->Draw("psame");
 LYZSTPCrefMultTPConlyAll->Draw("psame"); 
 FQDTPCrefMultTPConlyAll->Draw("psame"); 
 //hydro->Draw("lsame");
 l1->Draw("same"); 
 
//================================================================================================================================================= 
 
} // end of void figure3()

//=================================================================================================================================================

void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
 // Shift original TGraphErrors along x-axis by Double_t shift.
 
 if(!ge){cout<<"!ge"<<endl;exit(0);}
 
 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 Double_t y = 0.;
 for(Int_t p=0;p<nPoints;p++)
 { 
  ge->GetPoint(p,x,y);
  x+=shift;
  ge->SetPoint(p,x,y);
 } // end of for(Int_t p=0;p<nPoints;p++)
} // end of void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)

//=================================================================================================================================================

TH1D* StyleHistogram()
{
 // Style histogram:
 
 const Int_t nBins = 80; // to be improved - hardwired 10
 Double_t min = 0.; // in [GeV/c] // to be improved - hardwired 0
 Double_t max = 80.; // in [GeV/c] // to be improved - hardwired 5.75
 //TString xTitle  = "hadronic cross section #sigma/#sigma_{geo}";
 TString xTitle  = "centrality percentile";
 //TString xTitle  = "% most central"; 
 //TString binLabelsX[nBins] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%"};
 //TString binLabelsX[nBins] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%"};
 TString yTitle  = "v_{2}";
   
 TH1D* hist = new TH1D("","",nBins,min,max);
 hist->SetXTitle(xTitle.Data());
 //for(Int_t b=0;b<nBins;b++)
 //{
 // hist->GetXaxis()->SetBinLabel(b+1,binLabelsX[b].Data());
 //} 
 hist->SetYTitle(yTitle.Data());
 hist->SetMinimum(0.000001); // minimum y-value - to be improved
 hist->SetMaximum(0.12); // maximum y-value - to be improved
 
 return hist;
 
} // end of TH1D* StyleHistogramVsPt()
  
// =====================================================================================

void SetFlowStyle()
{
 // Set style which will affect all plots.
 
 gStyle->Reset();
 // gStyle->SetOptitle(0);
 // gStyle->SetOptStat(0);
 //gStyle->SetOptDate(1);
 // gStyle->SetPalette(8,0);  // (1,0)
 gStyle->SetPalette(1);  // (1,0)
 gStyle->SetDrawBorder(0);
 // gStyle->SetFillColor(0);  // kills palete ???
 gStyle->SetCanvasColor(0);
 gStyle->SetPadColor(0);
 // gStyle->SetFillColor(0); // otherwize it affects Fill colors later
 gStyle->SetFrameFillColor(0);
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetFrameLineWidth(2);
 // gStyle->SetFrameFillStyle(4000);
 gStyle->SetPadBorderMode(0);
 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);
 gStyle->SetPadBottomMargin(0.15);
 gStyle->SetPadLeftMargin(0.15);
 gStyle->SetHistLineWidth(2);
 gStyle->SetFuncWidth(2);
 gStyle->SetLineWidth(2);
 gStyle->SetLabelSize(0.05,"xyz");
 gStyle->SetLabelOffset(0.01,"y");
 gStyle->SetLabelColor(kBlack,"xyz");
 gStyle->SetTitleSize(0.06,"xyz");
 gStyle->SetTitleOffset(1.3,"y");
 gStyle->SetTitleFillColor(0);
 gStyle->SetLineWidth(2);  
 gStyle->SetHistLineColor(1);
 gStyle->SetTextColor(1);
 gStyle->SetTitleTextColor(1);
 TGaxis::SetMaxDigits(4);
 gStyle->SetOptStat(0); // removes stat. box from all histos
 gROOT->ForceStyle();

} // end of void SetFlowStyle()

