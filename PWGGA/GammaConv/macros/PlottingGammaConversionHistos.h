/***********************************************************************************************
*** provided by Gamma Conversion Group, PWGGA, 									******
***	    Friederike Bock, fbock@physi.uni-heidelberg.de ***							******
************************************************************************************************/


/************************************************************************************************
/************************************************************************************************
/* This header was created to make things easier with making plots for the gamma conversion group.
	it offers you several functions for drawing and styling your histogramms.
/************************************************************************************************
/************************************************************************************************

  The functions are 
	- StyleSettingsThesis
	- StyleSettings
	- GammaScalingHistogramm

	- DrawAutoGammaHistos
	- DrawAutoGammaHisto
	- DrawAutoGammaHisto2D

	- DrawRatioGammaHisto

	- DrawCutGammaHisto
	- DrawCutGammaHistos

	- DrawGammaLines
*/

#include <Riostream.h>

//extern TRandom *kgRandom;
//extern TBenchmark *kgBenchmark;
//extern TSystem *kgSystem;


// ---------------------------- Function definiton --------------------------------------------------------------------------------------------


/* StyleSettingsThesis will make some standard settings for gStyle 
*/
void StyleSettingsThesis(){
	//gStyle->SetOptTitle(kFALSE);
	gStyle->SetOptDate(0);   //show day and time
	gStyle->SetOptStat(0);  //show statistic
	gStyle->SetPalette(1,0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTextSize(0.5);
	gStyle->SetLabelSize(0.03,"xyz");
	gStyle->SetLabelOffset(0.002,"xyz");
	gStyle->SetTitleFontSize(0.04);
	gStyle->SetTitleOffset(1,"y");
	gStyle->SetTitleOffset(0.7,"x");		
	gStyle->SetCanvasColor(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetLineWidth(0.01);

	gStyle->SetPadTopMargin(0.03);
	gStyle->SetPadBottomMargin(0.09);
	gStyle->SetPadRightMargin(0.03);
	gStyle->SetPadLeftMargin(0.13);

	
	TGaxis::SetMaxDigits(3);
}


/* StyleSettings will make some standard settings for gStyle 
*/
void StyleSettings(){
	//gStyle->SetOptTitle(kFALSE);
	gStyle->SetOptDate(0);   //show day and time
	gStyle->SetOptStat(0);  //show statistic
	gStyle->SetPalette(1,0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTextSize(0.5);
	gStyle->SetLabelSize(0.03,"xyz");
	gStyle->SetLabelOffset(0.002,"xyz");
	gStyle->SetTitleFontSize(0.04);
	gStyle->SetTitleOffset(1,"y");
	gStyle->SetTitleOffset(0.7,"x");		
	gStyle->SetCanvasColor(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetLineWidth(0.01);

	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadBottomMargin(0.1);
	gStyle->SetPadRightMargin(0.08);
	gStyle->SetPadLeftMargin(0.12);

	
	TGaxis::SetMaxDigits(3);
}


void SetPlotStyle() {
    const Int_t nRGBs = 5;
    const Int_t nCont = 255;

    Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
    gStyle->SetNumberContours(nCont);
}




// GammaScalingHistogram will scale the histogram by "Factor" 
void GammaScalingHistogramm(TH1 *histo, Float_t Factor){
	histo->Sumw2();
	histo->Scale(Factor);
} 

// GammaScalingHistogram will scale the histogram by "Factor" 
void GammaScalingHistogramm(TH2 *histo, Float_t Factor){
	histo->Sumw2();
	histo->Scale(Factor);
} 

void StylingSliceHistos(TH1 *histo, Float_t markersize){
   	histo->SetMarkerStyle(22);
	histo->SetMarkerSize(markersize);
}



/************************************ OLD VERSION ************************************

/*void DrawAutoGammaHistos( TH1* histo1, 
					 TH1*histo2, 
					 const char *Title, const char *XTitle, const char *YTitle, 
					 Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum, 
					 Bool_t YRange, Float_t YMin ,Float_t YMax, 
					 Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		if(maxRangeR < histo2->GetMaximum()){
			maxRangeR = histo2->GetMaximum();
			}
		Double_t minRangeR = histo1->GetMinimum();		
		if(minRangeR > histo2->GetMinimum()){
			minRangeR = histo2->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		if(maxRangeR < histo2->GetMaximum()){
			maxRangeR = histo2->GetMaximum();
			}
		Double_t minRangeR = histo1->GetMinimum();		
		if(minRangeR > histo2->GetMinimum()){
			minRangeR = histo2->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title != ""){
		histo1->SetTitle(Title);
	}
	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetTitleSize(0.025);	
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.8);
	histo1->GetXaxis()->SetTitleSize(0.025);	
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->DrawCopy("e,hist");
	histo2->SetLineColor(2);
	histo2->DrawCopy("e,hist,same");
	leg1 = new TLegend( 0.6,0.82,0.92,0.9);
	leg1->SetTextSize(0.04);			
	leg1->SetFillColor(0);
	leg1->AddEntry(histo1,("Data"));
	leg1->AddEntry(histo2,("MC"));
	leg1->Draw();

}
*/


/* DrawAutoGammaHistos is function used for styling the histograms of the gamma conversion group for two histos and standart settings
	* histo1 - first histogram (Data)
	* histo2 - second histogram (MC)
	* Title - histogram title
	* XTitle - X-axis title
	* YTitle - Y-axis title
	* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
	*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
	*YMinimum - this will be used if YRangeMax is set
	*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
		- will be set to kFAlSE if YRangeMax is set
	*YMin - minimum Y
	*YMax - maximum Y
	*XRange 	= kTRUE will Cut x-axis by XMin and XMax
	*XMin - minimum Y
	*XMax - maximum Y
*/ 

void DrawAutoGammaHistos( TH1* histo1, 
					 TH1*histo2, 
					 const char *Title, const char *XTitle, const char *YTitle, 
					 Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum, 
					 Bool_t YRange, Float_t YMin ,Float_t YMax, 
					 Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		if(maxRangeR < histo2->GetMaximum()){
			maxRangeR = histo2->GetMaximum();
			}
		Double_t minRangeR = histo1->GetMinimum();		
		if(minRangeR > histo2->GetMinimum()){
			minRangeR = histo2->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		if(maxRangeR < histo2->GetMaximum()){
			maxRangeR = histo2->GetMaximum();
			}
		Double_t minRangeR = histo1->GetMinimum();		
		if(minRangeR > histo2->GetMinimum()){
			minRangeR = histo2->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title != ""){
		histo1->SetTitle(Title);
	}else{histo1->SetTitle();
		 histo2->SetTitle();}
	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetLabelSize(0.035);
	histo1->GetYaxis()->SetTitleSize(0.04);	
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.8);
	histo1->GetXaxis()->SetTitleSize(0.04);	
	histo1->GetXaxis()->SetLabelSize(0.035);
	histo1->DrawCopy("e,hist");
	histo2->SetLineColor(2);
	histo2->DrawCopy("e,hist,same");
	leg1 = new TLegend( 0.7,0.87,0.97,0.97);
	leg1->SetTextSize(0.04);			
	leg1->SetFillColor(0);
	leg1->AddEntry(histo1,("Data"));
	leg1->AddEntry(histo2,("MC"));
	leg1->Draw();

}


/* DrawAutoGammaHisto is function used for styling a histograma of the gamma conversion group with standart settings
	* histo1 - first histogram (Data)
	* Title - histogram title
	* XTitle - X-axis title
	* YTitle - Y-axis title
	* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
	*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
	*YMinimum - this will be used if YRangeMax is set
	*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
		- will be set to kFAlSE if YRangeMax is set
	*YMin - minimum Y
	*YMax - maximum Y
	*XRange 	= kTRUE will Cut x-axis by XMin and XMax
	*XMin - minimum Y
	*XMax - maximum Y
*/ 
void DrawAutoGammaHisto( TH1* histo1, 
					const char *Title, const char *XTitle, const char *YTitle,
					Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum,
					Bool_t YRange, Float_t YMin ,Float_t YMax,  
					Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}

	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title != ""){
		histo1->SetTitle(Title);
	}
	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetTitleSize(0.025);	
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.8);
	histo1->GetXaxis()->SetTitleSize(0.025);
	histo1->GetXaxis()->SetLabelSize(0.02);	
	histo1->DrawCopy("e,hist");
}

/*DrawAutoGammaHisto2D is a function for drawing a 2D-histogram of the gamma conversion group
	* histo - histogramm which need to be drawn
	* Title - histogram title
	* XTitle - X- axis-title
	* YTitle - Y-axis-title
	* Input - Legend 
	* YRange - if kTRUE will scale by YMin and YMay
	* YMin  - Y minimum
	* YMax - Y maximum
	* XRange - if kTRUE will scale by XMin and XMax
	* XMin - X minimum
	* XMax - X maximum
*/
void DrawAutoGammaHisto2D(	TH2 *histo,  
						const char *Title, const char *XTitle, const char *YTitle, const char *Input,
						Bool_t YRange, Float_t YMin ,Float_t YMax, 
						Bool_t XRange, Float_t XMin, Float_t XMax) {


	if (YRange && XRange){
		histo->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if ( !YRange && XRange){
		histo->GetXaxis()->SetRangeUser(XMin, XMax);	
	}

	if (YRange && !XRange){
		histo->GetYaxis()->SetRangeUser(YMin, YMax);
	}
	
	if(Title != ""){
		histo->SetTitle(Title);
	}
	if(XTitle ! =""){
		histo->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo->SetYTitle(YTitle);
	}
	histo->GetYaxis()->SetTitleSize(0.025);	
	histo->GetYaxis()->SetLabelSize(0.02);
	histo->GetXaxis()->SetLabelSize(0.02);
	histo->GetYaxis()->SetDecimals();
	histo->GetYaxis()->SetTitleOffset(1.8);
	histo->GetXaxis()->SetTitleSize(0.025);	
	histo->DrawCopy("colz");
	if(Input!=""){
	leg2 = new TLegend(0.6,0.82,0.83,0.9);
	leg2->SetTextSize(0.04);			
	leg2->SetFillColor(0);
	leg2->AddEntry(histo,(Input));
	leg2->Draw("same");
	}
}


/* DrawRatioGammaHisto is function used for styling the ratio-histograms of the gamma conversion group
	* histo1 - histogram 
	* Title - histogram title
	* XTitle - X-axis title
	* YTitle - Y-axis title
	* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
	*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
	*YMinimum - this will be used if YRangeMax is set
	*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
		- will be set to kFAlSE if YRangeMax is set
	*YMin - minimum Y
	*YMax - maximum Y
	*XRange 	= kTRUE will Cut x-axis by XMin and XMax
	*XMin - minimum Y
	*XMax - maximum Y
*/ 
void DrawRatioGammaHisto( TH1* histo1, 
					const char *Title, const char *XTitle, const char *YTitle,
					Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum,
					Bool_t YRange, Float_t YMin ,Float_t YMax,  
					Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}

	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title != ""){	histo1->SetTitle(Title);
	}else{	histo1->SetTitle();}

	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetTitleSize(0.04);	
	histo1->GetYaxis()->SetLabelSize(0.03);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.3);
	histo1->GetXaxis()->SetTitleOffset(1.1);
	histo1->GetXaxis()->SetTitleSize(0.04);
	histo1->GetXaxis()->SetLabelSize(0.03);	
	histo1->DrawCopy("e");
}

/* DrawCutGammaHistos is function used for styling the Cut-histograms of the gamma conversion group for 4 histos combined
	* histo1 - histogram Data
	* histo2 - histogram Data Comparision
	* histo3 - histogram MC
	* histo4 - histogram MC Comparision
	* Title - histogram title
	* XTitle - X-axis title
	* YTitle - Y-axis title
	* Legend1 - additional Legend for histo2
	* Legend2 - additional Legend for histo4	
	* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
	*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
	*YMinimum - this will be used if YRangeMax is set
	*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
		- will be set to kFAlSE if YRangeMax is set
	*YMin - minimum Y
	*YMax - maximum Y
	*XRange 	= kTRUE will Cut x-axis by XMin and XMax
	*XMin - minimum Y
	*XMax - maximum Y
*/ 
void DrawCutGammaHistos( TH1* histo1, TH1* histo2, 
					TH1* histo3, TH1*histo4, 
					const char *Title, const char *XTitle, const char *YTitle, const char *Legend1, const char *Legend2,
					Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum,
					Bool_t YRange, Float_t YMin ,Float_t YMax,  
					Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo2->GetMaximum();
		if(maxRangeR < histo4->GetMaximum()){
			maxRangeR = histo4->GetMaximum();
			}
		Double_t minRangeR = histo2->GetMinimum();		
		if(minRangeR > histo4->GetMinimum()){
			minRangeR = histo4->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo2->GetMaximum();
		if(maxRangeR < histo4->GetMaximum()){
			maxRangeR = histo4->GetMaximum();
			}
		Double_t minRangeR = histo2->GetMinimum();		
		if(minRangeR > histo4->GetMinimum()){
			minRangeR = histo4->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title != ""){
		histo1->SetTitle(Title);
	}
	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetTitleSize(0.025);		
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.8);
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->GetXaxis()->SetTitleSize(0.025);	
	histo1->Draw("e,hist");
	
	histo2->SetLineColor(15);
	histo2->Draw("e,hist,same");

	histo3->SetLineColor(2);
	histo3->Draw("e,hist,same");
	
	histo4->SetLineColor(46);		
	histo4->Draw("e,hist,same");

	leg1 = new TLegend( 0.6,0.82,0.92,0.9);
	leg1->SetTextSize(0.02);			
	leg1->SetFillColor(0);
	leg1->AddEntry(histo1,("Data"));
	leg1->AddEntry(histo2,(Legend1));
	leg1->AddEntry(histo3,("MC"));
	leg1->AddEntry(histo4,(Legend2));

	leg1->Draw();

	
}

/* DrawCutGammaHisto is function used for styling the Cut-histograms of the gamma conversion group for 2 histos combined
	* histo1 - histogram Data
	* histo2 - histogram Data Comparision
	* Title - histogram title
	* XTitle - X-axis title
	* YTitle - Y-axis title
	* Legend - additional Legend for histo2
	* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
	*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
	*YMinimum - this will be used if YRangeMax is set
	*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
		- will be set to kFAlSE if YRangeMax is set
	*YMin - minimum Y
	*YMax - maximum Y
	*XRange 	= kTRUE will Cut x-axis by XMin and XMax
	*XMin - minimum Y
	*XMax - maximum Y
*/ 
void DrawCutGammaHisto( TH1* histo1, TH1* histo2, 
					const char *Title, const char *XTitle, const char *YTitle, const char *Legend,
					Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum,
					Bool_t YRange, Float_t YMin ,Float_t YMax,  
					Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo2->GetMaximum();
		Double_t minRangeR = histo2->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo2->GetMaximum();
		Double_t minRangeR = histo2->GetMinimum();				
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title != ""){
		histo1->SetTitle(Title);
	}
	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetTitleSize(0.025);	
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.8);
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->GetXaxis()->SetTitleSize(0.025);	
	histo1->Draw("e,hist");
	
	histo2->SetLineColor(15);
	histo2->Draw("e,hist,same");

	leg1 = new TLegend( 0.6,0.82,0.92,0.9);
	leg1->SetTextSize(0.04);			
	leg1->SetFillColor(0);
	leg1->AddEntry(histo1,("Data"));
	leg1->AddEntry(histo2,(Legend));
	leg1->Draw();

}

/* DrawRatioGammaHisto is function used for styling the ratio-histograms of the gamma conversion group
	* histo1 - histogram 
	* Title - histogram title
	* XTitle - X-axis title
	* YTitle - Y-axis title
	* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
	*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
	*YMinimum - this will be used if YRangeMax is set
	*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
		- will be set to kFAlSE if YRangeMax is set
	*YMin - minimum Y
	*YMax - maximum Y
	*XRange 	= kTRUE will Cut x-axis by XMin and XMax
	*XMin - minimum Y
	*XMax - maximum Y
*/ 
void DrawResolutionGammaHisto( TH1* histo1, 
					const char *Title, const char *XTitle, const char *YTitle,
					Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum,
					Bool_t YRange, Float_t YMin ,Float_t YMax,  
					Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}

	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title != ""){
		histo1->SetTitle(Title);
	}else { histo1->SetTitle();}
	
	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetTitleSize(0.04);	
	histo1->GetYaxis()->SetLabelSize(0.03);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.3);
	histo1->GetXaxis()->SetTitleOffset(1.1);
	histo1->GetXaxis()->SetTitleSize(0.04);
	histo1->GetXaxis()->SetLabelSize(0.03);	
	histo1->DrawCopy("e");
}

/*DrawAutoGammaHisto2Dres is a function for drawing a resolution 2D-histogram of the gamma conversion group
	* histo - histogramm which need to be drawn
	* Title - histogram title
	* XTitle - X- axis-title
	* YTitle - Y-axis-title
	* Input - Legend 
	* YRange - if kTRUE will scale by YMin and YMay
	* YMin  - Y minimum
	* YMax - Y maximum
	* XRange - if kTRUE will scale by XMin and XMax
	* XMin - X minimum
	* XMax - X maximum
*/
void DrawAutoGammaHisto2DRes(	TH2 *histo,  
						const char *Title, const char *XTitle, const char *YTitle, const char *Input,
						Bool_t YRange, Float_t YMin ,Float_t YMax, 
						Bool_t XRange, Float_t XMin, Float_t XMax) {


	if (YRange && XRange){
		histo->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if ( !YRange && XRange){
		histo->GetXaxis()->SetRangeUser(XMin, XMax);	
	}

	if (YRange && !XRange){
		histo->GetYaxis()->SetRangeUser(YMin, YMax);
	}
	
	if(Title != ""){
		histo->SetTitle(Title);
	}
	if(XTitle ! =""){
		histo->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo->SetYTitle(YTitle);
	}
	histo->GetYaxis()->SetTitleSize(0.045);	
	histo->GetYaxis()->SetLabelSize(0.03);
	histo->GetXaxis()->SetLabelSize(0.03);
	histo->GetYaxis()->SetDecimals();
	histo->GetYaxis()->SetTitleOffset(1.5);
	histo->GetXaxis()->SetTitleSize(0.045);	
	histo->GetYaxis()->SetTitleOffset(1.5);
	histo->DrawCopy("colz");
	if(Input!=""){
	leg2 = new TLegend(0.6,0.82,0.83,0.9);
	leg2->SetTextSize(0.04);			
	leg2->SetFillColor(0);
	leg2->AddEntry(histo,(Input));
	leg2->Draw("same");
	}
}


/* DrawAutoGammaMesonHistos is function used for styling the histograms of the gamma conversion group for two histos and standart settings
* histo1 - first histogram
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 

void DrawAutoGammaMesonHistos( TH1* histo1, 
					 const char *Title, const char *XTitle, const char *YTitle, 
					 Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum, 
					 Bool_t YRange, Float_t YMin ,Float_t YMax, 
					 Bool_t XRange, Float_t XMin, Float_t XMax) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();		
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}
	
	if(Title != ""){
		histo1->SetTitle(Title);
	}else{
		histo1->SetTitle();
	}
	if(XTitle ! =""){
		histo1->SetXTitle(XTitle);
	}
	if(YTitle ! =""){
		histo1->SetYTitle(YTitle);
	}
	histo1->GetYaxis()->SetLabelSize(0.03);
	histo1->GetYaxis()->SetTitleSize(0.04);	
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.2);
	histo1->GetXaxis()->SetTitleSize(0.04);	
	histo1->GetXaxis()->SetLabelSize(0.03);
	
	histo1->DrawCopy("e1,p");
	
}

void DrawGammaSetMarker( TH1* histo1, 
				Style_t MarkerStyle, Size_t MarkerSize, Color_t MarkerColor, Color_t LineColor ) {
	histo1->SetMarkerStyle(MarkerStyle);
	histo1->SetMarkerSize(MarkerSize);
	histo1->SetMarkerColor(MarkerColor);
	histo1->SetLineColor(LineColor);	
}

void DrawGammaCanvasSettings( TCanvas* c1, Double_t LeftMargin, Double_t RightMargin, Double_t TopMargin, Double_t BottomMargin){
	c1->SetTickx();
	c1->SetTicky();
	c1->SetGridx(0);
	c1->SetGridy(0);
	c1->SetLogy(0);	
	c1->SetLeftMargin(LeftMargin);
	c1->SetRightMargin(RightMargin);
	c1->SetTopMargin(TopMargin);				
	c1->SetBottomMargin(BottomMargin);				
	c1->SetFillColor(0);
}

void DrawGammaPadSettings( TPad* pad1, Double_t LeftMargin, Double_t RightMargin, Double_t TopMargin, Double_t BottomMargin){
	pad1->SetFillColor(0);
	pad1->GetFrame()->SetFillColor(0);
	pad1->SetBorderMode(0);
	pad1->SetLeftMargin(LeftMargin);
	pad1->SetBottomMargin(BottomMargin);
	pad1->SetRightMargin(RightMargin);
	pad1->SetTopMargin(TopMargin);
	pad1->SetTickx();
	pad1->SetTicky();
	
}

void DrawGammaSetMarkerTGraph( TGraph* graph, 
					Style_t MarkerStyle, Size_t MarkerSize, Color_t MarkerColor, Color_t LineColor ) {
	graph->SetMarkerStyle(MarkerStyle);
	graph->SetMarkerSize(MarkerSize);
	graph->SetMarkerColor(MarkerColor);
	graph->SetLineColor(LineColor);	
}

void DrawGammaSetMarkerTGraphErr( TGraphErrors* graph, 
						 Style_t MarkerStyle, Size_t MarkerSize, Color_t MarkerColor, Color_t LineColor ) {
	graph->SetMarkerStyle(MarkerStyle);
	graph->SetMarkerSize(MarkerSize);
	graph->SetMarkerColor(MarkerColor);
	graph->SetLineColor(LineColor);	
}


void DrawGammaSetMarkerTGraphAsym( TGraphAsymmErrors* graph, 
						 Style_t MarkerStyle, Size_t MarkerSize, Color_t MarkerColor, Color_t LineColor ) {
	graph->SetMarkerStyle(MarkerStyle);
	graph->SetMarkerSize(MarkerSize);
	graph->SetMarkerColor(MarkerColor);
	graph->SetLineColor(LineColor);	
}


void DrawGammaSetMarkerTF1( TF1* fit1, 
						 Style_t LineStyle, Size_t LineWidth, Color_t LineColor ) {
	fit1->SetLineColor(LineColor);	
	fit1->SetLineStyle(LineStyle);	
	fit1->SetLineWidth(LineWidth);	
}
