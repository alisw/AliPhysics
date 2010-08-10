//void StyleSettingsThesis();
void StyleSettings();

      class TGradientParFunction {

      public:

         TGradientParFunction(int ipar, TF1 * f)  :
            fPar(ipar),
            fFunc(f)
         {}

         double operator() (double * x, double *) const
         {
            // evaluate gradient vector of functions at point x
            return fFunc->GradientPar(fPar,x);
         }

      private:

         unsigned int fPar;
         mutable TF1 * fFunc;
      };


void DrawGammaHisto(TH1*,TString, TString, TString, TString, Bool_t, Float_t, Float_t, Bool_t, Float_t, Float_t, Bool_t, Float_t, Float_t,Int_t);
/*
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


	//	TGaxis::SetMaxDigits(3);

}

/* DrawAutoGammaHisto is function used for styling a histograma of the gamma conversion group with standart settings
        * histo1 - first histogram (Data)
        * Title - histogram title
        * XTitle - X-axis title
        * YTitle - Y-axis title
        * YRangeMax     = kTRUE will scale by Maximum and Minimum Range in Y
        *YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE
        *YMinimum - this will be used if YRangeMax is set
        *YRange         = kTRUE will Cut y-axis by YMin and YMax
                - will be set to kFAlSE if YRangeMax is set
        *YMin - minimum Y
        *YMax - maximum Y
        *XRange         = kTRUE will Cut x-axis by XMin and XMax
        *XMin - minimum Y
        *XMax - maximum Y
*/
void DrawGammaHisto( TH1* histo1,
		     TString CutSelection,TString Title, TString XTitle, TString YTitle,
		     Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum,
		     Bool_t YRange, Float_t YMin ,Float_t YMax,
		     Bool_t XRange, Float_t XMin, Float_t XMax,Int_t bck) {


  
  
  /*
        if (YRangeMax && !XRange){
                YRange = kFALSE;
                Double_t maxRange_R = histo1->GetMaximum();
                Double_t minRange_R = histo1->GetMinimum();
                if(YMinimum > minRange_R){minRange_R = YMinimum;}
                histo1->GetYaxis()->SetRangeUser(minRange_R, maxRange_R*YMaxFactor);
        }
        if (YRangeMax && XRange){
                YRange = kFALSE;
                Double_t maxRange_R = histo1->GetMaximum();
                Double_t minRange_R = histo1->GetMinimum();
                if(YMinimum > minRange_R){minRange_R = YMinimum;}
                histo1->GetYaxis()->SetRangeUser(minRange_R, maxRange_R*YMaxFactor);
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
  */
  //	histo1->GetYaxis()->SetRangeUser(YMin, YMax);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	if(Title.Length() > 0){
	  histo1->SetTitle(Form("%s_%s",Title.Data(),CutSelection.Data()));
	}
	if(XTitle.Length() > 0){
	  histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
	  histo1->SetYTitle(YTitle.Data());
	}
        histo1->GetYaxis()->SetLabelSize(0.02);
        histo1->GetYaxis()->SetTitleSize(0.025);
        histo1->GetYaxis()->SetDecimals();
        histo1->GetYaxis()->SetTitleOffset(1.8);
        histo1->GetXaxis()->SetTitleSize(0.025);
        histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerColor(1);
	histo1->SetLineColor(1);
	histo1->SetMarkerSize(0.8);
	histo1->SetTitleOffset(1.2,"xy");		
	histo1->SetTitleSize(0.05,"xy");		
	histo1->GetYaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetNdivisions(507,kTRUE);
	if(bck){
	  histo1->SetLineStyle(1);		
	  histo1->SetLineColor(4);
	  histo1->SetLineWidth(2);
	  histo1->DrawCopy("hist,same");
	}else{
	  histo1->DrawCopy("e1,p");
	}
}
