/***********************************************************************************************
*** provided by Gamma Conversion Group, PWG4, 									******
***	    Friederike Bock, fbock@physi.uni-heidelberg.de ***							******
************************************************************************************************/

/************************************************************************************************
/************************************************************************************************
	This header contains the functions to plot additional things, like logos and Lines 
	in your histograms
/************************************************************************************************
/************************************************************************************************

  The functions are 
	- DrawAliceLogo
	- DrawAliceLogo1D
	- DrawAliceLogoPerformance
	- DrawAliceLogoPerformance2D
	- DrawAliceText
	
	- DrawStructure
	- DrawArmenteros
	- DrawdEdxLabels
	
	- DrawGammaLines
// */

#include <Riostream.h>

/*********************************************************************************************************
	DrawAliceLogoPi0Performace 
		* will draw you the ALICE Logo as well the text "ALICE Performance", "pp @ 7 TeV" and the date
		  which you can hand over
		  Float_t start_text_x, Float_t start_text_y, Float_t start_pi0_text_x, Float_t difference_text,
		  Float_t start_logo_x, Float_t start_logo_y, Float_t width_logo, 
		  Float_t text_size, 
		  Bool_t LowEnergy, Bool_t MCfile
		* float_t start_x, float_t start_y - give starting Point of Logo
		* float_t width_logo - gives the width of the logo
	 	* float_t text_height - gives you the heigth of the text
		* char * date - date handed over
**********************************************************************************************************
**********************************************************************************************************/

void DrawAliceLogoPi0Performance(Float_t start_text_x, Float_t start_text_y, Float_t start_pi0_text_x, Float_t difference_text,
						   Float_t start_logo_x, Float_t start_logo_y, Float_t width_logo, 
						   Float_t text_size, 
						   Bool_t LowEnergy, Bool_t MCfile , Bool_t RawData){

		string alicetext = "ALICE Performance";
		if(RawData) alicetext = "RAW DATA";
		
		TLatex *alice = new TLatex(start_text_x,(start_text_y+(2*difference_text)),Form("%s",alicetext.c_str())); // Bo: this was modified
		TLatex *energy = new TLatex(start_text_x,(start_text_y+difference_text),"pp @ 7TeV"); // Bo: this was modified
		TLatex *process = new TLatex(start_pi0_text_x, (start_text_y + (4* difference_text)), "#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}");
		TLatex *events = new TLatex(start_text_x,start_text_y,"1.233 x 10^{8} MinBias events"); // Bo: this was modified
		TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",start_logo_x,start_logo_y,(start_logo_x + width_logo),(start_logo_y + (width_logo*0.975)));
		
		
		if(LowEnergy){
			TLatex *alice = new TLatex(start_text_x,0.8,Form("%s",alicetext.c_str())); // Bo: this was modified
			TLatex *energy = new TLatex(start_text_x,0.75,"pp @ 7TeV"); // Bo: this was modified
			TLatex *process = new TLatex(0.35, 0.9, "#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}");
			TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",start_logo_x,start_logo_y,(start_logo_x + width_logo),(start_logo_y + (width_logo*0.975)));
		}
		
		if(MCfile != "") {
			TLatex *alice = new TLatex(start_text_x,0.8,Form("%s",alicetext.c_str())); // Bo: this was modified
			TLatex *energy = new TLatex(start_text_x,0.75,"pp @ 7TeV"); // Bo: this was modified
			TLatex *process = new TLatex(0.35, 0.9, "#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}");
			TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",start_logo_x,start_logo_y,(start_logo_x + width_logo),(start_logo_y + (width_logo*0.975)));
		}
		
		
		alice->SetNDC();
		alice->SetTextColor(2);
		alice->SetTextSize(text_size);
		alice->Draw();
		
		energy->SetNDC();
		energy->SetTextColor(2);
		energy->SetTextSize(text_size);
		energy->Draw();
		
		process->SetNDC(kTRUE); // <- use NDC coordinate
		process->SetTextSize(text_size);
		process->Draw();
		
		events->SetNDC();
		events->SetTextColor(2);
		events->SetTextSize(text_size);
		events->Draw();
		
		//myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
		myPadLogo->SetBorderMode(0);
		myPadLogo->SetBorderSize(2);
		myPadLogo->SetFrameBorderMode(0);
		myPadLogo->SetLeftMargin(0.0);
		myPadLogo->SetTopMargin(0.0);
		myPadLogo->SetBottomMargin(0.0);
		myPadLogo->SetRightMargin(0.0);
		TASImage *myAliceLogo = new TASImage("/home/fredi/Dokumente/CERN/ALICE/Photon/ALICE_logo.eps");
		myPadLogo->Draw();  // to take out for not using a logo.
		myPadLogo->cd();
		myAliceLogo->Draw();
}


/*************************************************************************************************
 DrawAliceLogo draws you the Alice logo + "work in progress" and " ALICE Performance" 
   be careful you have to set the path of the alice-logo for 	your system
	* float_t start_x, float_t start_y - give starting Point of Logo
	* float_t width_logo - gives the width of the logo
 	* float_t text_height - gives you the heigth of the text
**************************************************************************************************
**************************************************************************************************/

void DrawAliceLogo(Float_t start_x, Float_t start_y, Float_t width_logo, Float_t text_height){
		
		Float_t aliceStartY = start_y - text_height * 1.1;  
		TLatex *alice = new TLatex(start_x-0.0225,aliceStartY,"ALICE Performance"); // Bo: this was modified
	            alice->SetNDC();
	            alice->SetTextColor(2);
	            alice->SetTextFont(42);
	            alice->SetTextSize(text_height);
	            alice->SetLineWidth(2);
 		  alice->Draw("same");
	         
		Float_t wipStartY = start_y - text_height *(1 + 1.1);           
	            TLatex *wip = new TLatex((start_x-0.011),wipStartY,"work in progress"); // Bo: this was modified
	            wip->SetNDC();
	            wip->SetTextColor(2);
	            wip->SetTextFont(51);
	            wip->SetTextSize(text_height);
	            wip->SetLineWidth(2);
		  wip->Draw("same");
		Float_t endX = start_x + width_logo;
		Float_t endY = start_y + width_logo;

	            TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",start_x,start_y,endX,endY);
	            //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
	            myPadLogo->SetBorderMode(0);
	            myPadLogo->SetBorderSize(2);
	            myPadLogo->SetFrameBorderMode(0);
	            myPadLogo->SetLeftMargin(0.0);
	            myPadLogo->SetTopMargin(0.0);
	            myPadLogo->SetBottomMargin(0.0);
	            myPadLogo->SetRightMargin(0.0);
	            TASImage *myAliceLogo = new TASImage("/home/fredi/Dokumente/CERN/ALICE/Photon/ALICE_logo.eps");
		  myPadLogo->Draw();  // to take out for not using a logo.
	            myPadLogo->cd();
	            myAliceLogo->Draw("same");

}

/*************************************************************************************************
 DrawAliceLogo1D draws you the Alice logo + "work in progress" and " ALICE Performance" 
   be careful you have to set the path of the alice-logo for 	your system
	* float_t start_x, float_t start_y - give starting Point of Logo
	* float_t width_logo - gives the width of the logo
 	* float_t text_height - gives you the heigth of the text
**************************************************************************************************
**************************************************************************************************/

void DrawAliceLogo1D(Float_t start_x, Float_t start_y, Float_t width_logo, Float_t text_height){
		
		Float_t aliceStartY = start_y - text_height * 1.1;  
		TLatex *alice = new TLatex(start_x+0.008,aliceStartY,"ALICE Performance"); // Bo: this was modified
	            alice->SetNDC();
	            alice->SetTextColor(2);
	            alice->SetTextFont(42);
	            alice->SetTextSize(text_height);
	            alice->SetLineWidth(2);
 		  alice->Draw("same");
	         
		Float_t wipStartY = start_y - text_height *(1 + 1.1);           
	            TLatex *wip = new TLatex((start_x+0.022),wipStartY,"work in progress"); // Bo: this was modified
	            wip->SetNDC();
	            wip->SetTextColor(2);
	            wip->SetTextFont(51);
	            wip->SetTextSize(text_height);
	            wip->SetLineWidth(2);
		  wip->Draw("same");
		Float_t endX = start_x + width_logo;
		Float_t endY = start_y + width_logo;

	            TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",start_x,start_y,endX,endY);
	            //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
	            myPadLogo->SetBorderMode(0);
	            myPadLogo->SetBorderSize(2);
	            myPadLogo->SetFrameBorderMode(0);
	            myPadLogo->SetLeftMargin(0.0);
	            myPadLogo->SetTopMargin(0.0);
	            myPadLogo->SetBottomMargin(0.0);
	            myPadLogo->SetRightMargin(0.0);
	            TASImage *myAliceLogo = new TASImage("/home/fredi/Dokumente/CERN/ALICE/Photon/ALICE_logo.eps");
		  myPadLogo->Draw();  // to take out for not using a logo.
	            myPadLogo->cd();
	            myAliceLogo->Draw("same");
}


/*********************************************************************************************************
	DrawAliceLogoPerformace 
		* will draw you the ALICE Logo as well the text "ALICE Performance", "pp @ 7 TeV" and the date
		  which you can hand over
		* float_t start_x, float_t start_y - give starting Point of Logo
		* float_t width_logo - gives the width of the logo
	 	* float_t text_height - gives you the heigth of the text
		* char * date - date handed over
**********************************************************************************************************
**********************************************************************************************************/

void DrawAliceLogoPerformance(Float_t start_x, Float_t start_y, Float_t width_logo, Float_t text_height, Float_t decrease, char *date){
		
		Float_t aliceStartY = start_y - text_height * 1.1;  
		TLatex *alice = new TLatex((start_x-decrease),aliceStartY,"ALICE Performance"); // Bo: this was modified
	            alice->SetNDC();
	            alice->SetTextColor(2);
	            alice->SetTextFont(62);
	            alice->SetTextSize(text_height);
	            alice->SetLineWidth(2);
 		  alice->Draw("same");
		TLatex *pp7 = new TLatex((start_x+decrease*0.3),(aliceStartY-text_height*1.1),"pp @ 7 TeV"); // Bo: this was modified
	            pp7->SetNDC();
	            pp7->SetTextColor(2);
	            pp7->SetTextFont(62);	
	            pp7->SetTextSize(text_height);
	            pp7->SetLineWidth(2);
 		  pp7->Draw("same");
		TLatex *today = new TLatex((start_x-decrease*0),(aliceStartY-2*text_height*1.1),date); // Bo: this was modified
	            today->SetNDC();
	            today->SetTextColor(2);
	            today->SetTextFont(62);
	            today->SetTextSize(text_height);
	            today->SetLineWidth(2);
 		  today->Draw("same");
		
	         
		Float_t endX = start_x + width_logo;
		Float_t endY = start_y + width_logo;

	            TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",start_x,start_y,endX,endY);
	            //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
	            myPadLogo->SetBorderMode(0);
	            myPadLogo->SetBorderSize(2);
	            myPadLogo->SetFrameBorderMode(0);
	            myPadLogo->SetLeftMargin(0.0);
	            myPadLogo->SetTopMargin(0.0);
	            myPadLogo->SetBottomMargin(0.0);
	            myPadLogo->SetRightMargin(0.0);
	            TASImage *myAliceLogo = new TASImage("/home/fredi/Dokumente/CERN/ALICE/Photon/ALICE_logo.eps");
		  myPadLogo->Draw();  // to take out for not using a logo.
	            myPadLogo->cd();
	            myAliceLogo->Draw("same");

}
/*********************************************************************************************************
	DrawAliceLogoPerformace2D 
		* will draw you the ALICE Logo as well the text "ALICE Performance", "pp @ 7 TeV" and the date
		  which you can hand over
		* float_t start_x, float_t start_y - give starting Point of Logo
		* float_t width_logo - gives the width of the logo
	 	* float_t text_height - gives you the heigth of the text
		* float_t decrease - gives percentage in the canvas to which the text should be decreased in y 
						compared to the logo
		* char * date - date handed over
**********************************************************************************************************
**********************************************************************************************************/

void DrawAliceLogoPerformance2D(Float_t start_x, Float_t start_y, Float_t width_logo, Float_t text_height, Float_t decrease, char *date){
		
		Float_t aliceStartY = start_y + 3* text_height * 1.1;  
		TLatex *alice = new TLatex((start_x-decrease),aliceStartY,"ALICE Performance"); // Bo: this was modified
	            alice->SetNDC();
	            alice->SetTextColor(2);
	            alice->SetTextFont(62);
	            alice->SetTextSize(text_height);
	            alice->SetLineWidth(2);
 		  alice->Draw("same");
		TLatex *pp7 = new TLatex((start_x+decrease*0.3),(aliceStartY-text_height*1.1),"pp @ 7 TeV"); // Bo: this was modified
	            pp7->SetNDC();
	            pp7->SetTextColor(2);
	            pp7->SetTextFont(62);	
	            pp7->SetTextSize(text_height);
	            pp7->SetLineWidth(2);
 		  pp7->Draw("same");
		TLatex *today = new TLatex((start_x-decrease*0),(aliceStartY-2*text_height*1.1),date); // Bo: this was modified
	            today->SetNDC();
	            today->SetTextColor(2);
	            today->SetTextFont(62);
	            today->SetTextSize(text_height);
	            today->SetLineWidth(2);
 		  today->Draw("same");
		
	         
		Float_t endX = start_x + width_logo;
		Float_t endY = start_y + width_logo;

	            TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",start_x , (start_y - width_logo), endX,start_y);
	            //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
	            myPadLogo->SetBorderMode(0);
	            myPadLogo->SetBorderSize(2);
	            myPadLogo->SetFrameBorderMode(0);
	            myPadLogo->SetLeftMargin(0.0);
	            myPadLogo->SetTopMargin(0.0);
	            myPadLogo->SetBottomMargin(0.0);
	            myPadLogo->SetRightMargin(0.0);
	            TASImage *myAliceLogo = new TASImage("/home/fredi/Dokumente/CERN/ALICE/Photon/ALICE_logo.eps");
		  myPadLogo->Draw();  // to take out for not using a logo.
	            myPadLogo->cd();
	            myAliceLogo->Draw("same");

}


//***********************************************************************************************************
/* DrawAliceText draws you the "work in progress" and " ALICE Performance" 
	* float_t start_x, float_t start_y - give starting Point of Logo
 	* float_t text_height - gives you the heigth of the text
//***********************************************************************************************************/
//***********************************************************************************************************/

void DrawAliceText(Float_t start_x, Float_t start_y, Float_t text_height){
		Float_t aliceStartY = start_y - text_height * 1.1;  
		TLatex *alice = new TLatex(start_x,aliceStartY,"ALICE Performance"); // Bo: this was modified
	            alice->SetNDC();
	            alice->SetTextColor(2);
	            alice->SetTextFont(42);
	            alice->SetTextSize(text_height);
	            alice->SetLineWidth(2);
 		  alice->Draw();
	         
		Float_t wipStartY = start_y - text_height *(1 + 1.1);           
	            TLatex *wip = new TLatex(start_x,wipStartY,"work in progress"); // Bo: this was modified
	            wip->SetNDC();
	            wip->SetTextColor(2);
	            wip->SetTextFont(51);
	            wip->SetTextSize(text_height);
	            wip->SetLineWidth(2);
		  wip->Draw();
}

/***************************************************************************************************** 
	DrawStructure() draws the structure of the Inner Alice Detectors labeled for a xy - Plot

******************************************************************************************************
******************************************************************************************************/

void DrawStructure(){
  TLatex *ssdText = new TLatex(0.16,0.78,"SSD");
  ssdText->SetNDC();
  ssdText->SetTextFont(72);
  ssdText->SetTextSize(0.03);
  ssdText->SetLineWidth(4);
  ssdText->Draw();

  TLatex *sddText = new TLatex(0.14,0.7,"SDD");
  sddText->SetNDC();
  sddText->SetTextFont(72);
  sddText->SetTextSize(0.03);
  sddText->SetLineWidth(4);
  sddText->Draw();

  TLatex *spdText = new TLatex(0.14,0.25,"SPD");
  spdText->SetNDC();
  spdText->SetTextFont(72);
  spdText->SetTextSize(0.03);
  spdText->SetLineWidth(4);
  spdText->Draw();

  TLatex *tpcRText = new TLatex(0.52,0.115,"TPC Rods");
  tpcRText->SetNDC();
  tpcRText->SetTextFont(72);
  tpcRText->SetTextSize(0.03);
  tpcRText->SetLineWidth(4);
  tpcRText->Draw();

  TLatex *tpcIText = new TLatex(0.14,0.20,"TPC inner");
  tpcIText->SetNDC();
  tpcIText->SetTextFont(72);
  tpcIText->SetTextSize(0.03);
  tpcIText->SetLineWidth(4);
  tpcIText->Draw();

  TLatex *tpcIFText = new TLatex(0.14,0.17,"field cage");
  tpcIFText->SetNDC();
  tpcIFText->SetTextFont(72);
  tpcIFText->SetTextSize(0.03);
  tpcIFText->SetLineWidth(4);
  tpcIFText->Draw();

  TLatex *tpcIFVText = new TLatex(0.16,0.14,"vessel");
  tpcIFVText->SetNDC();
  tpcIFVText->SetTextFont(72);
  tpcIFVText->SetTextSize(0.03);
  tpcIFVText->SetLineWidth(4);
  tpcIFVText->Draw();

  TLatex *tpc1IText = new TLatex(0.705,0.20,"TPC inner");
  tpc1IText->SetNDC();
  tpc1IText->SetTextFont(72);
  tpc1IText->SetTextSize(0.03);
  tpc1IText->SetLineWidth(4);
  tpc1IText->Draw();

  TLatex *tpcICText = new TLatex(0.69,0.17,"containment");
  tpcICText->SetNDC();
  tpcICText->SetTextFont(72);
  tpcICText->SetTextSize(0.03);
  tpcICText->SetLineWidth(4);
  tpcICText->Draw();

  TLatex *tpcICVText = new TLatex(0.72,0.14,"vessel");
  tpcICVText->SetNDC();
  tpcICVText->SetTextFont(72);
  tpcICVText->SetTextSize(0.03);
  tpcICVText->SetLineWidth(10);
  tpcICVText->Draw();

  TLatex *tpcGasText = new TLatex(0.78,0.28,"TPC");
  tpcGasText->SetNDC();
  tpcGasText->SetTextFont(72);
  tpcGasText->SetTextSize(0.03);
  tpcGasText->SetLineWidth(10);
  tpcGasText->Draw();

  TLatex *tpcGas2Text = new TLatex(0.76,0.25,"drift gas");
  tpcGas2Text->SetNDC();
  tpcGas2Text->SetTextFont(72);
  tpcGas2Text->SetTextSize(0.03);
  tpcGas2Text->SetLineWidth(10);
  tpcGas2Text->Draw();


   TArrow * arrow = new TArrow(-150.8049,145.,-11.99843,38.629599,0.02,">"); //SSD arrow
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineWidth(2.);
   arrow->Draw();

   TArrow * arrow1 = new TArrow(-160.,105.,-11.99843,25.,0.02,">"); //SDD arrow
   arrow1->SetFillColor(1);
   arrow1->SetFillStyle(1001);
   arrow1->SetLineWidth(2.);
   arrow1->Draw();


   TArrow * arrow2 = new TArrow(-150,-108.,-7.,2.,0.02,">");  //SPD arrow
   arrow2->SetFillColor(1);
   arrow2->SetFillStyle(1001);
   arrow2->SetLineWidth(2.);
  arrow2->Draw();

   TArrow * arrow3 = new TArrow(-105.,-160.,-30.,-75.,0.02,">"); //TPC field cage vessel arrow
   arrow3->SetFillColor(1);
   arrow3->SetFillStyle(1001);
   arrow3->SetLineWidth(2.);
  arrow3->Draw();


  TArrow * arrow4 = new TArrow(50.,-178.,15.,-88.,0.02,">"); //TPC rods arrow
   arrow4->SetFillColor(1);
   arrow4->SetFillStyle(1001);
 arrow4->SetLineWidth(2.);
  arrow4->Draw();

  TArrow * arrow5 = new TArrow(130.,-130.,50.,-38.,0.02,">");// TPC inner constainment vessel arrow
   arrow5->SetFillColor(1);
   arrow5->SetFillStyle(1001);
 arrow5->SetLineWidth(2.);
  arrow5->Draw();


  TArrow * arrow6 = new TArrow(160.,-90.,130.,-70.,0.02,">");// TPC gas
   arrow6->SetFillColor(1);
   arrow6->SetFillStyle(1001);
 arrow6->SetLineWidth(2.);
  arrow6->Draw();

}
/*************************************************************************************************************
	DrawArmenteros() draws the labels for the particles in the Armenteros plot
**************************************************************************************************************/

void DrawArmenteros(){
  TLatex *k0s = new TLatex(0.45,0.87,"K^{0}_{s}");
  k0s->SetNDC();
  k0s->SetTextFont(62);
  k0s->SetTextSize(0.05);
  k0s->SetLineWidth(4);
  k0s->Draw();

  TLatex *lambda = new TLatex(0.71,0.47,"#Lambda");
  lambda->SetNDC();
  lambda->SetTextFont(62);
  lambda->SetTextSize(0.05);
  lambda->SetLineWidth(4);
  lambda->Draw();

  TLatex *lambdabar = new TLatex(0.22,0.48,"#bar{#Lambda}");
  lambdabar->SetNDC();
  lambdabar->SetTextFont(62);
  lambdabar->SetTextSize(0.05);
  lambdabar->SetLineWidth(4);
  lambdabar->Draw();

  TLatex *gamma = new TLatex(0.49,0.18,"#gamma");
  gamma->SetNDC();
  gamma->SetTextFont(62);
  gamma->SetTextSize(0.05);
  gamma->SetLineWidth(4);
  gamma->Draw();

}


/************************************************************************************************
	DrawdEdxLabels()
		* Float_t linewidth - will give the final linewidth in the plot

*************************************************************************************************
*************************************************************************************************/
void DrawdEdxLabel(){
		TLatex *pion = new TLatex(0.23,0.7,"#pi"); //text at pion line
	             pion->SetNDC();
	             pion->SetTextColor(2);
	             pion->SetTextFont(62);
	             pion->SetTextSize(0.03);
	             pion->SetLineWidth(2);	

		TLatex *kaon = new TLatex(0.375,0.7,"K"); //text at kaon line
	             kaon->SetNDC();
	             kaon->SetTextColor(4);
	             kaon->SetTextFont(62);
	             kaon->SetTextSize(0.03);
	             kaon->SetLineWidth(2);	

		TLatex *proton = new TLatex(0.455,0.7,"p"); //text at proton line
	             proton->SetNDC();
	             proton->SetTextColor(kGreen+3);
	             proton->SetTextFont(62);
	             proton->SetTextSize(0.03);
	             proton->SetLineWidth(2);	
		
		TLatex *electron = new TLatex(0.72,0.50,"e"); //text at electron line
	             electron->SetNDC();
	             electron->SetTextColor(1);
	             electron->SetTextFont(62);
	             electron->SetTextSize(0.03);
	             electron->SetLineWidth(2);	
		
		electron->Draw("same");
		kaon->Draw("same");
		proton->Draw("same");
		pion->Draw("same");
}

/* // DrawGammaLines will draw the lines in the histogram for you
		* start_x - starting point of drawing in x
		* endX - end point of drawing in x
		* start_y -starting point of drawing in y
		* endY - end point of drawing in y
		* linew - line width
*/
void DrawGammaLines(Float_t start_x, Float_t endX,
				 Float_t start_y, Float_t endY,
				Float_t linew){
			TLine * l1 = new TLine (start_x,start_y,endX,endY);
			l1->SetLineColor(4);
			l1->SetLineWidth(linew);
			l1->Draw("same");
}

