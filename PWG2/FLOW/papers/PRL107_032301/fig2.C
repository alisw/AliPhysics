// Plotting macro for v2{2}/e{2} and v3{2}/e{3} vs centrality (Figure 2b in "Other harmonics" paper).
// Objects holding measurements are following TGraphErrors:
// A.) v2:
//  1.) QC2_v2 = v2{2}, 2-particle Q-cumulant, TPConly tracks, all charges, reference flow;
// B.) v3:
//  1.) QC2_v3 = v3{2}, 2-particle Q-cumulant, TPConly tracks, all charges, reference flow;
// C.) Eccentricities:
//  1.) Epsilon2_Wound: e{2} for centralities 0-1%, 1-2%, 2-3%, 3-4% and 4-5%;
//  2.) Epsilon3_Wound: e{3} for centralities 0-1%, 1-2%, 2-3%, 3-4% and 4-5%;
//  3.) Epsilon2_Ncoll: e{2} for centralities 0-1%, 1-2%, 2-3%, 3-4% and 4-5%;
//  4.) Epsilon3_Ncoll: e{3} for centralities 0-1%, 1-2%, 2-3%, 3-4% and 4-5%.

static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

void fig2(Int_t rWrite = 0, Int_t rPerformance = 0)
{
 myOptions();
 gROOT->ForceStyle();   
 TDatime now;
 int iDate = now.GetDate();
 int iYear=iDate/10000;
 int iMonth=(iDate%10000)/100;
 int iDay=iDate%100;
 char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
 char cStamp1[25],cStamp2[25];
 sprintf(cStamp1,"%i %s %i",iDay,cMonth[iMonth-1],iYear);
 sprintf(cStamp2,"%i/%.2d/%i",iDay,iMonth,iYear);
    

    Double_t Scale_factor = 0.1;
    
    //===================================================================================================================
    // Eccentricities:    
    // values You Zhou
    
    // 64mb-r^2weight-e2 ,3,4,5W- 0-5%

    Double_t r2E2Cum2WCent[]= {0.0679313,0.0743901,0.0803903,0.0875735,0.0951992};
    Double_t r2E3Cum2WCent[]= {0.0678642,0.0744068,0.0792529,0.0854844,0.0902971};
    Double_t r2E4Cum2WCent[]= {0.0674344,0.0735389,0.0790758,0.084159,0.0894306};
    Double_t r2E5Cum2WCent[]= {0.067565,0.0726308,0.0773171,0.0831959,0.0865153};
    
    // 64mb-r^2weight-e2,3,4,5B -0-5%
    Double_t r2E2Cum2BCent[]= {0.106073,0.113651,0.122,0.133783,0.145521};
    Double_t r2E3Cum2BCent[]= {0.102277,0.106538,0.108152,0.11282,0.115907};
    Double_t r2E4Cum2BCent[]= {0.0990573,0.101967,0.105939,0.1086,0.110835};
    Double_t r2E5Cum2BCent[]= {0.0962965,0.0991143,0.100108,0.104052,0.10675};
    
    // 64mb-r^n-weight-e2,3,4,5 -Wound-0-5%
    Double_t rnE2Cum2WCent[]= {0.0698327,0.0758561,0.0759677,0.0904921,0.0949813}; //looks like the third point is not correct
    Double_t rnE3Cum2WCent[]= {0.0822765,0.0921005,0.0961503,0.10326,0.11043};
    Double_t rnE4Cum2WCent[]= {0.0935775,0.103275,0.113841,0.123022,0.13041};
    Double_t rnE5Cum2WCent[]= {0.107025,0.12147,0.131068,0.140556,0.14948};
    
    // 64mb-r^nweight-e2,3,4,5-Binary-0-5%
    Double_t rnE2Cum2BCent[]= {0.108189,0.116104,0.120313,0.137813,0.146756};
    Double_t rnE3Cum2BCent[]= {0.110244,0.119488,0.120747,0.121706,0.127783};
    Double_t rnE4Cum2BCent[]= {0.119342,0.122386,0.12686,0.134593,0.135264};
    Double_t rnE5Cum2BCent[]= {0.11959,0.129879,0.135822,0.140458,0.14582};
    
    // Constantin
    //Double_t CGCE2[] = {0.0689396, 0.0801229, 0.0941998, 0.109266, 0.126151};
    //Double_t CGCE3[] = {0.0641664, 0.0679419, 0.0715216, 0.0745122, 0.0780054};
    
    // CGC
    Double_t CGCE2[] = {0.0742864, 0.0841532, 0.0983202, 0.114496, 0.130138};
    Double_t CGCE3[] = {0.0703217, 0.0752441, 0.0792893, 0.0829627, 0.0868921};
    
    Double_t x_eps[] = {0.5,1.5,2.5,3.5,4.5};
    Int_t nPoints_e = sizeof(x_eps)/sizeof(Double_t);
    Double_t sr2E2Cum2BCent[5]; 
    Double_t sr2E3Cum2BCent[5]; 
    Double_t sr2E2Cum2WCent[5]; 
    Double_t sr2E3Cum2WCent[5]; 
    Double_t sCGCE2[5]; 
    Double_t sCGCE3[5]; 

    for(Int_t p=0;p<5;p++)
    {
        sr2E2Cum2BCent[p] = 0.23*r2E2Cum2BCent[p];
        sr2E3Cum2BCent[p] = 0.193*r2E3Cum2BCent[p];
        sr2E2Cum2WCent[p] = 0.35*r2E2Cum2WCent[p];
        sr2E3Cum2WCent[p] = 0.263*r2E3Cum2WCent[p];
        sCGCE2[p] = 0.283*CGCE2[p];
        sCGCE3[p] = 0.263*CGCE3[p];
    }
    
    TGraph * e2ncoll = new TGraph(nPoints_e,x_eps,sr2E2Cum2BCent);
    TGraph * e2nwound = new TGraph(nPoints_e,x_eps,sr2E2Cum2WCent);
    TGraph * e3ncoll = new TGraph(nPoints_e,x_eps,sr2E3Cum2BCent);
    TGraph * e3nwound = new TGraph(nPoints_e,x_eps,sr2E3Cum2WCent);
    TGraph * e2cgc = new TGraph(nPoints_e,x_eps,sCGCE2);
    TGraph * e3cgc = new TGraph(nPoints_e,x_eps,sCGCE3);
    
    e3cgc->SetMarkerStyle(kOpenCircle);  
    e3cgc->SetMarkerColor(kBlue);
    e3cgc->SetLineColor(kBlue);
    e3cgc->SetLineStyle(0);
    e3cgc->SetLineWidth(2);

    e2cgc->SetMarkerStyle(kOpenSquare);
    e2cgc->SetMarkerColor(kRed);
    e2cgc->SetLineColor(kRed);
    e2cgc->SetLineStyle(0);
    e2cgc->SetLineWidth(2);

    e3ncoll->SetMarkerStyle(kOpenCircle);  
    e3ncoll->SetMarkerColor(kBlue);
    e3ncoll->SetLineColor(kBlue);
    e3ncoll->SetLineStyle(3);
    e3ncoll->SetLineWidth(2);
    
    e2ncoll->SetMarkerStyle(kOpenSquare);
    e2ncoll->SetMarkerColor(kRed);
    e2ncoll->SetLineColor(kRed);
    e2ncoll->SetLineStyle(3);
    e2ncoll->SetLineWidth(2);

    e3nwound->SetMarkerStyle(kOpenCircle);  
    e3nwound->SetMarkerColor(kBlue);
    e3nwound->SetLineColor(kBlue);
    e3nwound->SetLineStyle(kDashed);
    e3nwound->SetLineWidth(2);
    
    e2nwound->SetMarkerStyle(kOpenSquare);
    e2nwound->SetMarkerColor(kRed);
    e2nwound->SetLineColor(kRed);
    e2nwound->SetLineStyle(kDashed);
    e2nwound->SetLineWidth(2);
    

    //===================================================================================================================
    //===============================================================================================================================
    // SP_v2 = v2{SP}, Scalar Product, TPConly tracks, all charges, reference flow:
 
    Double_t xSP_v2_etaGap10[] = {0.5.,1.5,2.5,3.5,4.5};
    Double_t ySP_v2_etaGap10[] = {0.02120716,0.02424403,0.02788042,0.03212382,0.03541490};                              
    Double_t xErrSP_v2_etaGap10[5] = {0.};
    Double_t yErrSP_v2_etaGap10[] = {0.00026893,0.00025887,0.00024884,0.00024094,0.00023628};
    Double_t ySystErrSP_v2_etaGap10[5] = {0.};
    Int_t nPointsSP_v2_etaGap10 = sizeof(xSP_v2_etaGap10)/sizeof(Double_t);         
    TGraphErrors *SP_v2_etaGap10 = new TGraphErrors(nPointsSP_v2_etaGap10,xSP_v2_etaGap10,ySP_v2_etaGap10,
                                                    xErrSP_v2_etaGap10,yErrSP_v2_etaGap10);
    SP_v2_etaGap10->SetMarkerStyle(kFullCircle);
    SP_v2_etaGap10->SetMarkerColor(kRed);
    
    for (Int_t p=0;p<5;p++) {
        Double_t estsys = 0.03*ySP_v2_etaGap10[p];
        ySystErrSP_v2_etaGap10[p] = TMath::Sqrt(estsys*estsys+yErrSP_v2_etaGap10[p]*yErrSP_v2_etaGap10[p]);
    }
    
    TGraphErrors *SP_v2_2 = new TGraphErrors(nPointsSP_v2_etaGap10,xSP_v2_etaGap10,ySP_v2_etaGap10,
                                                    xErrSP_v2_etaGap10,ySystErrSP_v2_etaGap10);
    SP_v2_2->SetMarkerStyle(kFullCircle);
    SP_v2_2->SetMarkerColor(kRed-9);
    SP_v2_2->SetLineColor(kRed-9); 
    SP_v2_2->SetFillColor(kRed-9); 
    SP_v2_2->SetFillStyle(1001); 

    //===============================================================================================================================
    // SP_v3 = v3{SP}, Scalar Product, TPConly tracks, all charges, reference flow:
      
    // SP_v3_etaGap10 = v3{SP}:
    Double_t xSP_v3_etaGap10[] = {0.5.,1.5,2.5,3.5,4.5};
    Double_t ySP_v3_etaGap10[] = {0.01890397,0.02056964,0.02079589,0.02170537,0.02297969};                              
    Double_t xErrSP_v3_etaGap10[5] = {0.};
    Double_t yErrSP_v3_etaGap10[] = {0.00027369,0.00027552,0.00028550,0.00028783,0.00028729};
    Double_t ySystErrSP_v3_etaGap10[5] = {0.};
    Int_t nPointsSP_v3_etaGap10 = sizeof(xSP_v3_etaGap10)/sizeof(Double_t);         
    TGraphErrors *SP_v3_etaGap10 = new TGraphErrors(nPointsSP_v3_etaGap10,xSP_v3_etaGap10,ySP_v3_etaGap10,
                                                    xErrSP_v3_etaGap10,yErrSP_v3_etaGap10);
    SP_v3_etaGap10->SetMarkerStyle(kFullSquare);
    SP_v3_etaGap10->SetMarkerColor(kBlue);
    
    
    for (Int_t p=0;p<5;p++) {
        Double_t estsys = 0.03*ySP_v3_etaGap10[p];
        ySystErrSP_v3_etaGap10[p] = TMath::Sqrt(estsys*estsys+yErrSP_v3_etaGap10[p]*yErrSP_v3_etaGap10[p]);
    }
    
    TGraphErrors *SP_v3_2 = new TGraphErrors(nPointsSP_v3_etaGap10,xSP_v3_etaGap10,ySP_v3_etaGap10,
                                             xErrSP_v3_etaGap10,ySystErrSP_v3_etaGap10);
    SP_v3_2->SetMarkerStyle(kFullCircle);
    SP_v3_2->SetMarkerColor(kBlue-9);
    SP_v3_2->SetLineColor(kBlue-9); 
    SP_v3_2->SetFillColor(kBlue-9); 
    SP_v3_2->SetFillStyle(1001);
                
 //=========================================================================================================================

 // Style histogramm:
 TH1F *myBlankHisto = new TH1F("myBlankHisto","Blank Histogram",5,0.,5.);
 myBlankHisto->SetXTitle("centrality percentile");
 myBlankHisto->SetNdivisions(505,"y");
 myBlankHisto->GetYaxis()->SetTitleOffset(0.6);
 //myBlankHisto->GetYaxis()->SetLabelSize(0.05);
 //myBlankHisto->GetYaxis()->SetTitleSize(0.07);
 myBlankHisto->GetYaxis()->SetRangeUser(0.015,0.045);
 //myBlankHisto->GetXaxis()->SetRangeUser(0.,81.);    

 // Main canvas:   
 TCanvas *myCan = new TCanvas("myCan",cStamp1);
 myCan->Draw();
 myCan->cd(0);   
 gPad->SetBottomMargin(0.15); 
 myBlankHisto->Draw();
 TLegend *myLegend = new TLegend(0.12,0.56,0.34,0.86);
 myLegendSetUp(myLegend,0.04); 
 myLegend->AddEntry(SP_v2_etaGap10,"v_{2}{2, |#Delta#eta| > 1}","p"); 
 myLegend->AddEntry(SP_v3_etaGap10,"v_{3}{2, |#Delta#eta| > 1}","p"); 
 myLegend->AddEntry(e2cgc,"k_{1} #times #varepsilon_{n}^{CGC}{2}","l"); 
 myLegend->AddEntry(e2nwound,"k_{2} #times #varepsilon_{n}^{W}{2}","l"); 
 //myLegend->AddEntry(e2ncoll,"k_{3} #times #varepsilon_{n}^{B}{2}","l");  
 myLegend->Draw("same");
 
    SP_v2_2->Draw("same3"); 
    SP_v3_2->Draw("same3"); 
    SP_v2_etaGap10->Draw("psame"); 
    SP_v3_etaGap10->Draw("psame"); 
    e2cgc->Draw("lsame");
    e3cgc->Draw("lsame");
    //e2ncoll->Draw("lsame");
    //e3ncoll->Draw("lsame");
    e2nwound->Draw("lsame");
    e3nwound->Draw("lsame");
    //e2ncoll->Draw("lsame");
    //e3ncoll->Draw("lsame"); 
 
 
  return;    


 
 
 // Final drawing (order is important!):
 //f5pCorrelator->Draw("psame");
 //QC2_v2->Draw("psame");
 //QC4_v2->Draw("psame");
 //QC2_v3->Draw("psame"); 
 //QC2_v3_same->Draw("psame"); 
 //SP_v3_etaGap02->Draw("psame");    
 SP_TPCrefMultTPConlyAll2->Draw("same3");
 SP_TPCrefMultTPConlyAll->Draw("lsame");
 SP_TPCrefMultTPConlyAll->Draw("psame");   
 //SP_v3_etaGap10->Draw("psame"); 
 SP_v4_etaGap10->Draw("psame"); 
 SP_v3_corr2->Draw("same3");
 SP_v3_corr->Draw("lsame");
 SP_v3_corr->Draw("psame");
 QC4_v3->Draw("psame");
 //QC4_v3_same->Draw("psame");
 QC5_v3_squared->Draw("psame");
 ZDC_v3->Draw("psame");
 //LuzOll->Draw("lsame"); 
 
 // Figure 1b: 
 myCan->cd(2);   
 gPad->SetTopMargin(0.); 
 gPad->SetBottomMargin(0.15); 
 gPad->SetRightMargin(0.001); 
 TH1D* styleHistPad2 = myBlankHisto->Clone("styleHistPad2");
 styleHistPad2->GetYaxis()->SetRangeUser(0.,0.4799);  
 styleHistPad2->GetXaxis()->SetTitle("centrality percentile");  
 styleHistPad2->GetYaxis()->SetTitle("v_{n}/#varepsilon_{n}");   
 styleHistPad2->Draw();
 TLatex *b = new TLatex(0.95,0.93,"(b)");
 b->SetNDC();
 b->SetTextSize(0.05);
 b->Draw();
 TLegend *myLegendPad2 = new TLegend(0.6,0.64,0.86,0.94);
 myLegendSetUp(myLegendPad2,0.04); 
 myLegendPad2->AddEntry(Ratio_v22_e2CGC,"v_{2}{2, #Delta#eta > 1}/#varepsilon_{2}^{CGC}{2}","p");
 myLegendPad2->AddEntry(Ratio_v32_e2CGC,"v_{3}{2, #Delta#eta > 1}/#varepsilon_{3}^{CGC}{2}","p");
 myLegendPad2->AddEntry(Ratio_v22_e2w,"v_{2}{2, #Delta#eta > 1}/#varepsilon_{2}^{W}{2}","p");
 myLegendPad2->AddEntry(Ratio_v32_e2w,"v_{3}{2, #Delta#eta > 1}/#varepsilon_{3}^{W}{2}","p");
 //myLegendPad2->AddEntry(Ratio_v32_e3b,"v_{3}{2}/#varepsilon_{3,Ncoll}{2}","l"); 
 //myLegendPad2->AddEntry(Ratio_v34_e2b,"v_{3}{4}/#varepsilon_{3,Ncoll}{4}","l"); 
 //myLegendPad2->AddEntry(Ratio_v34_e2w,"v_{3}{4}/#varepsilon_{3}{4}","p"); 
 //myLegendPad2->AddEntry(Ratio_v32_e3w,"v_{3}{2}/#varepsilon_{3}{2}","p"); 
 myLegendPad2->Draw("same");
 //Ratio_v34_e2b->Draw("same3"); 
 //Ratio_v34_e2b->Draw("lsameX"); 
 //Ratio_v32_e3b->Draw("same3"); 
 //Ratio_v32_e3b->Draw("lsameX"); 
 //Ratio_v34_e2w->Draw("psame"); 
 //Ratio_v32_e3w->Draw("psame"); 
 
 
 
 
 Ratio_v22_e2CGC2->Draw("same3");
 Ratio_v22_e2w2->Draw("same3");
 Ratio_v32_e2CGC2->Draw("same3");
 Ratio_v32_e2w2->Draw("same3");   
 Ratio_v32_e2CGC->Draw("lsameX");
 Ratio_v22_e2CGC->Draw("lsameX");
 Ratio_v32_e2w->Draw("lsameX");
 Ratio_v22_e2w->Draw("lsameX");   
 Ratio_v32_e2CGC->Draw("psame");
 Ratio_v22_e2CGC->Draw("psame");
 Ratio_v32_e2w->Draw("psame");
 Ratio_v22_e2w->Draw("psame"); 
 
 return;


    /*
    TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.47,0.37,0.62,0.60);
    //myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
    myPadSetUp(myPadLogo,0,0,0,0);
    myPadLogo->Draw();
    myPadLogo->cd();
    TASImage *myAliceLogo = new TASImage("alice_logo_transparent.png");
    //.TASImage *myAliceLogo = new TASImage("alice_logo.pdf");     
    myAliceLogo->Draw();
    */
    
    if (rPerformance){
        TLatex *alice = new TLatex(0.75,0.34,"Performance");
        alice->SetNDC();
        alice->SetTextColor(myDarkRed);
        //    alice->SetTextFont(42);
        alice->SetTextSize(0.05);
        alice->SetLineWidth(2);
        alice->Draw();
        
        TText *date = new TText(0.78,0.28,cStamp2);
        date->SetNDC();
        date->SetTextFont(42);
        date->SetTextSize(0.04);
        date->Draw();
        
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.77,0.37,0.92,0.60);
        myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
        myPadSetUp(myPadLogo,0,0,0,0);
        myPadLogo->Draw();
        myPadLogo->cd();
        TASImage *myAliceLogo = new TASImage("alice_logo_transparent.png");
        //.TASImage *myAliceLogo = new TASImage("alice_logo.pdf");
        
        myAliceLogo->Draw();
    }
    if (rWrite == 1)  myCan->SaveAs("fig_template.pdf");
    if (rWrite == 2)  myCan->SaveAs("fig_template.png");
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
    currentLegend->SetTextFont(42);
    currentLegend->SetBorderSize(0);
    currentLegend->SetFillStyle(0);
    currentLegend->SetFillColor(0);
    currentLegend->SetMargin(0.25);
    currentLegend->SetTextSize(currentTextSize);
    currentLegend->SetEntrySeparation(0.5);
    return;
}

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
    currentPad->SetLeftMargin(currentLeft);
    currentPad->SetTopMargin(currentTop);
    currentPad->SetRightMargin(currentRight);
    currentPad->SetBottomMargin(currentBottom);
    return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
                  int currentMarkerStyle=20, int currentMarkerColor=0,
                  int currentLineStyle=1, int currentLineColor=0, int currentFillColor=0, int currentFillStyle=0)
{
    currentGraph->SetMarkerSize(currentMarkerSize);
    currentGraph->SetMarkerStyle(currentMarkerStyle);
    currentGraph->SetMarkerColor(currentMarkerColor);
    currentGraph->SetLineStyle(currentLineStyle);
    currentGraph->SetLineWidth(2);
    currentGraph->SetLineColor(currentLineColor);
    currentGraph->SetFillColor(currentFillColor);
    currentGraph->SetFillStyle(currentFillStyle);
    return;
}

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

void myOptions(Int_t lStat=0){
    // Set gStyle
    int font = 42;
    // From plain
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    //
    gStyle->SetDrawBorder(0);
    gStyle->SetTextFont(font);
    gStyle->SetStatFont(font);
    gStyle->SetStatFontSize(0.05);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetTickLength(0.02,"y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelFont(font,"xyz"); 
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetTitleFont(font,"xyz");  
    gStyle->SetTitleOffset(1.0,"xyz");  
    gStyle->SetTitleSize(0.06,"xyz");  
    gStyle->SetMarkerSize(1); 
    gStyle->SetPalette(1,0); 
    if (lStat){
        gStyle->SetOptTitle(1);
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(1111);
    }
    else {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
    }
}
