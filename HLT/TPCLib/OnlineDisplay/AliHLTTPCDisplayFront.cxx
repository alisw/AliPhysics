// $Id$

/** \class AliHLTTPCDisplayFront
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayFront
//
// Display class for the HLT TPC-Pad events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TAttText.h>
#include <TAxis.h>
#include <TCanvas.h>

#ifdef use_aliroot
#include <TClonesArray.h>
#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPCParam.h>
#endif

#include "AliHLTStdIncludes.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTransform.h"

#include "AliHLTTPCDisplayMain.h"
#include "AliHLTTPCDisplayFront.h"
#include "AliHLTTPCDisplayPad.h"
#if __GNUC__ >= 3
using namespace std;
#endif

#define TESTCODE 0

ClassImp(AliHLTTPCDisplayFront)

//____________________________________________________________________________________________________
AliHLTTPCDisplayFront::AliHLTTPCDisplayFront(AliHLTTPCDisplayMain* display) {
    // constructor
    fDisplay = display;

    fNTimes = display->GetNTimeBins();
    
    fBinY[0] = 0;
    fBinY[1] = AliHLTTPCTransform::GetNRows() - 1;
#if TESTCODE    
    fBinX[0] = (-4) * AliHLTTPCTransform::GetNPads(fBinY[1]);      
    fBinX[1] = (4) * AliHLTTPCTransform::GetNPads(fBinY[1]);
    Int_t Bins =  (8 * AliHLTTPCTransform::GetNPads(fBinY[1]) ) + 1;
#else
    fBinX[0] = 0;      
    fBinX[1] = AliHLTTPCTransform::GetNPads(fBinY[1]);
#endif    
    fTmpEvent = 0;    
    
    Int_t fBinningFaktor = 4 ;
    



 
#if TESTCODE
    fHistfront = new TH2F("fHistfront","FrontView of selected slice;Pad #;Padrow #",Bins,fBinX[0],fBinX[1],fBinY[1]+1,fBinY[0],fBinY[1]);
#else
    fHistfront = new TH2F("fHistfront","FrontView of selected slice;Pad #;Padrow #",fBinX[1]+1,fBinX[0],fBinX[1],fBinY[1]+1,fBinY[0],fBinY[1]);
    fHistfrontcl = new TH1F("fHistfrontcl","cvcv;ddd;kkk",fBinX[1]+1,fBinX[0],fBinX[1]);
    gStyle->SetPalette(1);
#endif

    fHistfront->SetOption("COLZ");  
    fHistfront->SetTitleSize(0.03);
    fHistfront->GetXaxis()->SetLabelSize(0.03);
    fHistfront->GetXaxis()->SetTitleSize(0.03);
    fHistfront->GetYaxis()->SetLabelSize(0.03);
    fHistfront->GetYaxis()->SetTitleSize(0.03);
}

//____________________________________________________________________________________________________
AliHLTTPCDisplayFront::~AliHLTTPCDisplayFront() {
    // destructor   
    if ( fHistfront ){
        delete fHistfront;
        fHistfront = NULL;
    }
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayFront::Reset(){  
    fHistfront->Reset(); 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayFront::Save(){  
    fCanvas->SaveAs("HLT-FrontView.eps"); 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayFront::Fill(Int_t patch, ULong_t dataBlock, ULong_t dataLen){
  // Fill Pad Histogram

  Int_t timeSwitch = fDisplay->GetFrontDataSwitch();

  // --- TEST CODE beginn
  Int_t fBinning = 8; // == 1/0.125
  Int_t fBinningFaktor = 4; // binning / 2 because of width half

#if TESTCODE
    // use sum
    if (timeSwitch == 0) {
      for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){

	Int_t width_half = AliHLTTPCTransform::GetNPads(row) * AliHLTTPCTransform::GetPadLength(row) * fBinningFaktor;
	Int_t pad_corrected_loop = (Int_t) ( AliHLTTPCTransform::GetPadLength(row) *fBinning );

 	for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){

	  Int_t pad_corrected = ( pad * AliHLTTPCTransform::GetPadLength(row) * fBinning ) - width_half;

	  UInt_t timeSum = 0;
	  for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	    timeSum += fDisplay->fRawData[row][pad][timeBin];
	  } // end time
	  for (Int_t ii=0;ii < pad_corrected_loop; ii++){
	    pad_corrected++;
	    fHistfront->Fill(pad_corrected,row,(Int_t) timeSum);
	  }
	} // end pad
      }  // end row
    } // end use sum

    return;
    // --- TEST CODE end
#endif

  // !!
  // !!  DO unrolling because of cache effects (avoid cache trashing) !!
  // !!
  
  if ( fDisplay->GetZeroSuppression() ){
    // use sum
    if (timeSwitch == 0) {
      for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){
	for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){
	  UInt_t timeSum = 0;
	  for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	    timeSum += fDisplay->fRawDataZeroSuppressed[row][pad][timeBin];
	  } // end time
	  
	  fHistfront->Fill(pad,row,(Int_t) timeSum);
	} // end pad
      }  // end row
    } // end use sum
    
    // use average
    else if (timeSwitch == 1){
      for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){
	for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){
	  UInt_t timeSum = 0;
	  Int_t NTimeBins = 0;
	  Float_t timeAverage = 0.;
	  for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	    timeSum += fDisplay->fRawDataZeroSuppressed[row][pad][timeBin];
	    NTimeBins++;
	  } // end time
	  
	  if (NTimeBins <= 0) 
	    HLTFatal("Division by Zero - NTimeBins == 0");
	  else
	    timeAverage = ((Float_t) timeSum) / ((Float_t) NTimeBins);
	  
	  fHistfront->Fill(pad,row, timeAverage);
	} // end pad
      }  // end row
    }// end use average
    
    // use maximum
    else if (timeSwitch == 2){
      for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){
	for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){
	  UInt_t timeMax = 0;
	  for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	    UInt_t tmp = fDisplay->fRawDataZeroSuppressed[row][pad][timeBin];
	    if (tmp > timeMax) timeMax = tmp;
	  } // end time
	
	  fHistfront->Fill(pad,row,(Int_t) timeMax);
	} // end pad
      }  // end row
    }// end use maximum
  }  // end - if ( fDisplay->GetZeroSuppression() ){

  else {
    // use sum
    if (timeSwitch == 0) {
      for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){
	for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){
	  UInt_t timeSum = 0;
	  for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	    timeSum += fDisplay->fRawData[row][pad][timeBin];
	  } // end time
	  
	  fHistfront->Fill(pad,row,(Int_t) timeSum);
	} // end pad
      }  // end row
    } // end use sum
    
    // use average
    else if (timeSwitch == 1){
      for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){
	for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){
	  UInt_t timeSum = 0;
	  Int_t NTimeBins = 0;
	  Float_t timeAverage = 0.;
	  for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	    timeSum += fDisplay->fRawData[row][pad][timeBin];
	    NTimeBins++;
	  } // end time
	  
	  if (NTimeBins <= 0) 
	    HLTFatal("Division by Zero - NTimeBins == 0");
	  else
	    timeAverage = ((Float_t) timeSum) / ((Float_t) NTimeBins);
	  
	  fHistfront->Fill(pad,row, timeAverage);
	} // end pad
      }  // end row
    }// end use average
    
    // use maximum
    else if (timeSwitch == 2){
      for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){
	for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){
	  UInt_t timeMax = 0;
	  for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	    UInt_t tmp = fDisplay->fRawData[row][pad][timeBin];
	    if (tmp > timeMax) timeMax = tmp;
	  } // end time
	
	  fHistfront->Fill(pad,row,(Int_t) timeMax);
	} // end pad
      }  // end row
    }// end use maximum
  } // end - else of  if ( fDisplay->GetZeroSuppression() ){


  if (fDisplay->ExistsClusterData()){
    for (patch=0; patch < 6; patch++){
      AliHLTTPCSpacePointData *points = fDisplay->GetSpacePointDataPointer(fDisplay->GetSlicePadRow(),patch);
      if(!points) return;

      cout << "fill" << patch << endl;

      Float_t xyz[3];
      for(Int_t i=0; i< fDisplay->GetNumberSpacePoints(fDisplay->GetSlicePadRow(),patch); i++){
	xyz[0] = points[i].fX;
	xyz[1] = points[i].fY;
	xyz[2] = points[i].fZ;
	Int_t padRow = AliHLTTPCTransform::GetPadRow(xyz[0]);
	

	// select padrow to fill in histogramm
	//      if (padRow == AliHLTTPCTransform::GetPadRow(xyz[0])){
	AliHLTTPCTransform::LocHLT2Raw(xyz, 0, padRow);
	fHistfrontcl->Fill(xyz[1],padRow);
	//      }
      }
    }
  } // END if (fDisplay->ExistsClusterData()){


}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayFront::Draw(){
    fDisplay->GetCanvasFront()->cd();
    fDisplay->GetCanvasFront()->Clear();

    if (fDisplay->GetSplitFront()){
      fDisplay->GetCanvasFront()->Divide(1,2);
      fDisplay->GetCanvasFront()->cd(1);
    }

    Char_t title[256];
    sprintf(title,"FrontView of selected slice%d",fDisplay->GetSlicePadRow());

    // Keep Zoom
    if (!fDisplay->GetKeepView() ){
	fBinY[0] = 0;
	fBinY[1] = AliHLTTPCTransform::GetNRows() - 1;
	fBinX[0] = 0;      
	fBinX[1] = AliHLTTPCTransform::GetNPads(fBinY[1]);
    }

    fHistfront->SetAxisRange(fBinX[0],fBinX[1]);
    fHistfront->SetAxisRange(fBinY[0],fBinY[1],"Y");

    fHistfront->SetTitle(title);
    fHistfront->SetStats(kFALSE);
    fHistfront->Draw("COLZ");

    if ( fDisplay->ExistsClusterData() ){
     	fHistfrontcl->SetAxisRange(fBinX[0],fBinX[1]);
	fHistfrontcl->SetAxisRange(fBinY[0],fBinY[1],"Y");
	fHistfrontcl->SetStats(kFALSE);
	fHistfrontcl->SetMarkerStyle(28);
	fHistfrontcl->SetMarkerSize(2);
	fHistfrontcl->SetMarkerColor(1);
	fHistfrontcl->Draw("psame");

	cout << "draw" << endl;
    }

    if (fDisplay->GetSplitFront()){
      fDisplay->GetCanvasFront()->cd(2);
      fDisplay->GetPadPointer()->fHistpad2->Draw();
    }

    fDisplay->GetCanvasFront()->Modified();
    fDisplay->GetCanvasFront()->Update();

    // Select Pad
    fDisplay->GetCanvasFront()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliHLTTPCDisplayMain",(void*) fDisplay,"ExecPadEvent(Int_t,Int_t,Int_t,TObject*)");
    // Keep Zoom
    fDisplay->GetCanvasFront()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliHLTTPCDisplayFront",(void*) this,"ExecEvent(Int_t,Int_t,Int_t,TObject*)");
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayFront::ExecEvent(Int_t event, Int_t px, Int_t py, TObject *selected){
   // Saves the Zoom Position of the Histogram 

   // - Mouse down on Axis : StartPoint of Range
   if (event == 1 && selected->InheritsFrom("TAxis"))
       fTmpEvent = 1;

   // - Mouse pressed on Axis : Real Zoom process not only click
   else if (event == 21 && selected->InheritsFrom("TAxis") && fTmpEvent == 1)
       fTmpEvent = 21;

   // - Mouse pressed on Axis : Still pressed
   else if (event == 21 && selected->InheritsFrom("TAxis") && fTmpEvent == 21) 
       return;

    // - Mouse up on Axis : End Point of Rangex
   else if (event == 11 && selected->InheritsFrom("TAxis") && fTmpEvent == 21){
       TAxis * axis = (TAxis*) selected;
       
       if (selected == fHistfront->GetXaxis() || selected == fHistfront->GetXaxis()){ 
	   fBinX[0] = axis->GetFirst() -1;
	   fBinX[1] = axis->GetLast() -1;
       }
       else {
	   fBinY[0] = axis->GetFirst() -1;
	   fBinY[1] = axis->GetLast() -1;
       }

       fTmpEvent = 0; 
   }
   else fTmpEvent = 0;
}
