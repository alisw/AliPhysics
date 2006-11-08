// $Id$

/** \class AliHLTTPCDisplayPadRow
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayPadRow
//
// Display class for the HLT TPC-PadRow events.
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
#include <TPolyMarker3D.h>

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
#include "AliHLTTPCDigitReaderPacked.h"
#include "AliHLTTPCDigitReaderRaw.h"


#include "AliHLTTPCDisplayMain.h"
#include "AliHLTTPCDisplayPadRow.h"
#include "AliHLTTPCDisplayPad.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCDisplayPadRow)

//____________________________________________________________________________________________________
AliHLTTPCDisplayPadRow::AliHLTTPCDisplayPadRow(AliHLTTPCDisplayMain* display) {
    // constructor
    fDisplay = display;

    fNTimes = AliHLTTPCTransform::GetNTimeBins();
    
    fBinX[0] = 0;      
    fBinX[1] = fDisplay->GetNPads(); 
    fBinY[0] = 0;
    fBinY[1] = fNTimes -1;
    fTmpEvent = 0;         

    for (Int_t ii = 0;ii < 20; ii++)
      fpmarr[ii] = 0;

    Int_t maxpads = 150;
    Int_t padbinning = maxpads*10;

    fHistraw = new TH2F("fHistraw","Selected PadRow with found Clusters;Pad #;Timebin #",maxpads,0,maxpads-1,fNTimes,0,fNTimes-1);
    fHistrawcl = new TH1F("fHistrawcl","cvcv;ddd;kkk",padbinning,0,maxpads-1);

    fHistraw->SetOption("COLZ");    
    fHistraw->SetTitleSize(0.03);
    fHistraw->GetXaxis()->SetLabelSize(0.03);
    fHistraw->GetXaxis()->SetTitleSize(0.03);
    fHistraw->GetYaxis()->SetLabelSize(0.03);
    fHistraw->GetYaxis()->SetTitleSize(0.03);

    gStyle->SetPalette(1);
}

//____________________________________________________________________________________________________
AliHLTTPCDisplayPadRow::~AliHLTTPCDisplayPadRow() {
    // destructor     
    if ( fHistraw ){
	delete fHistraw;
	fHistraw = NULL;
    }
    if ( fHistrawcl ){
	delete fHistrawcl;
	fHistrawcl = NULL;
    }
    for (Int_t ii=19; ii <= 0; ii--){ // try fwd
      if ( fpmarr[ii] ){
	delete [] fpmarr[ii];
	fpmarr[ii] = NULL;
      }
    }
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPadRow::Reset(){  
  fHistraw->Reset();   
  fHistrawcl->Reset(); 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPadRow::Save(){  
  fDisplay->GetCanvasCharge()->SaveAs("HLT-PadRowView.eps"); 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPadRow::Fill(Int_t patch, ULong_t dataBlock, ULong_t dataLen){
    // Fill PadRow Histogram
    
#if defined(HAVE_TPC_MAPPING)
    AliHLTTPCDigitReaderRaw digitReader(0);

    bool readValue = true;
    Int_t rowOffset = 0;
    
    Int_t padRow = fDisplay->GetPadRow();
    Int_t slice = fDisplay->GetSlicePadRow();

    // Initialize RAW DATA
    Int_t firstRow = AliHLTTPCTransform::GetFirstRow(patch);
    Int_t lastRow = AliHLTTPCTransform::GetLastRow(patch);

    // Outer sector, patches 2, 3, 4, 5 -  start counting in patch 2 with row 0
    if ( patch >= 2 ) rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

    // Initialize block for reading packed data
    void* tmpdataBlock = (void*) dataBlock;
    digitReader.InitBlock(tmpdataBlock,dataLen,firstRow,lastRow,patch,slice);

    readValue = digitReader.Next();

    if (!readValue){	
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplayPadRow::Fill","Read first value") << "No value in data block" << ENDLOG;
	return;
    }

    if ( fDisplay->GetZeroSuppression() ){
      for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(padRow); pad++){
	//	for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){

	  fHistraw->Fill(pad,timeBin,fDisplay->fRawDataZeroSuppressed[padRow][pad][timeBin]);
	} // end time
      }  // end pad
      
    }  // end - if ( fDisplay->GetZeroSuppression() ){

    else {
      for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(padRow); pad++){
	//	for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
	for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){
	  fHistraw->Fill(pad,timeBin,fDisplay->fRawData[padRow][pad][timeBin]);
	} // end time
      }  // end pad
    }  // end - else of if ( fDisplay->GetZeroSuppression() ){


    // FILL PADROW 3D --- Initialize the colorbins
    if ( fDisplay->Get3DSwitchPadRow() ){
	for (UInt_t ii=0;ii < 20;ii++){
	    fbinct[ii] = 0;
	    fcolorbin[ii] = 0;
	}

	// read number of entries in colorbin
	while ( readValue ){ 
	    
	    Int_t row = digitReader.GetRow() + rowOffset;
	    
	    if (row == padRow){    
		UInt_t charge = digitReader.GetSignal();
		
		for (UInt_t ii=0;ii < 19;ii++){
		    if ( charge > (ii*15) && charge <= ((ii*15) + 15) )	fcolorbin[ii]++;
		}
                // larger than 19 * 15	
		if (charge > 285 ) fcolorbin[19]++;
	    }

	    // read next value
	    readValue = digitReader.Next();
      
	    if(!readValue) break; //No more value
	} 

	// Initialize fpmarr[color][3*colorbin[ii]] 
	for (Int_t ii = 0; ii < 20; ii++) {
	  if ( fpmarr[ii] ) delete[] fpmarr[ii];
	  fpmarr[ii] = new Float_t[fcolorbin[ii]*3]; 
	}

	// Rewind the raw reader and fill the polymarker3D
	digitReader.InitBlock(tmpdataBlock,dataLen,firstRow,lastRow,patch,slice);
	
	readValue = digitReader.Next();
    } // END if (fDisplay->Get3DSwitchPadRow())

    // -- Fill Raw Data
    while ( readValue ){ 

	Int_t row = digitReader.GetRow() + rowOffset;

	// select padrow to fill in histogramm
	if (row == padRow){    
	    UChar_t pad = digitReader.GetPad();
	    UShort_t time = digitReader.GetTime();
	    UInt_t charge = digitReader.GetSignal();
	    Float_t xyz[3];
	    //fHistraw->Fill(pad,time,charge);

	    if ( fDisplay->Get3DSwitchPadRow() ) {
		// Transform raw coordinates to local coordinates
		AliHLTTPCTransform::RawHLT2Global(xyz, slice, padRow, pad, time);

		for (UInt_t ii=0;ii < 19;ii++){
		    if ( charge > (ii*15) && charge <= ((ii*15) + 15) ){
			fpmarr[ii][fbinct[ii]] = xyz[0];
			fpmarr[ii][fbinct[ii]+1] = xyz[1];
			fpmarr[ii][fbinct[ii]+2] = xyz[2];
			fbinct[ii] += 3;
		    }
		}
		// larger than 19 * 15
		if (charge > 285 ) {
		    fpmarr[19][fbinct[19]] = xyz[0];
		    fpmarr[19][fbinct[19]+1] = xyz[1];
		    fpmarr[19][fbinct[19]+2] = xyz[2];
		    fbinct[19] += 3;
		}
	    } // END if (fSwitch3DPadRow)
	
	} // END if (row == padRow){ 
	
	// read next value
	readValue = digitReader.Next();
      
	//Check where to stop:
	if(!readValue) break; //No more value
    } 
    
     if (fDisplay->ExistsClusterData()){
	AliHLTTPCSpacePointData *points = fDisplay->GetSpacePointDataPointer(slice,patch);
	if(!points) return;
	
	Float_t xyz[3];
	for(Int_t i=0; i< fDisplay->GetNumberSpacePoints(slice,patch); i++){
	    xyz[0] = points[i].fX;
	    xyz[1] = points[i].fY;
	    xyz[2] = points[i].fZ;
	
	    // select padrow to fill in histogramm
	    if (padRow == AliHLTTPCTransform::GetPadRow(xyz[0])){
		AliHLTTPCTransform::LocHLT2Raw(xyz, slice, padRow);
		fHistrawcl->Fill(xyz[1],xyz[2]);
	    }
	}
     } // END if (fDisplay->ExistsClusterData()){
#else //! defined(HAVE_TPC_MAPPING)
      HLTFatal("DigitReaderRaw not available - check your build");
#endif //defined(HAVE_TPC_MAPPING)
 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPadRow::Draw(){

    fDisplay->GetCanvasPadRow()->cd();
    fDisplay->GetCanvasPadRow()->Clear();

    if (fDisplay->GetSplitPadRow()){
	fDisplay->GetCanvasPadRow()->Divide(1,2);
	fDisplay->GetCanvasPadRow()->cd(1);
    }
  
    Char_t title[256];
    sprintf(title,"Selected PadRow %d with found Clusters",fDisplay->GetPadRow());
  
    // Keep Zoom
    if (!fDisplay->GetKeepView() || (fBinX[1]>fDisplay->GetNPads()) || (fBinY[1]>fNTimes)){
	fBinX[0] = 0;
	fBinX[1] = fDisplay->GetNPads();
	fBinY[0] = 0;
	fBinY[1] = fNTimes - 1;
    }
    
    fHistraw->SetAxisRange(fBinX[0],fBinX[1]); 
    fHistraw->SetAxisRange(fBinY[0],fBinY[1],"Y");
    fHistraw->SetTitle(title);
    fHistraw->SetStats(kFALSE);
    fHistraw->Draw("COLZ");
    
    if ( fDisplay->ExistsClusterData() ){
      //cout << "XX" <<  endl;
	fHistrawcl->SetAxisRange(fBinX[0],fBinX[1]);
	fHistrawcl->SetAxisRange(fBinY[0],fBinY[1],"Y");
	fHistrawcl->SetStats(kFALSE);
	fHistrawcl->SetMarkerStyle(28);
	fHistrawcl->SetMarkerSize(2);
	fHistrawcl->SetMarkerColor(1);
	fHistrawcl->Draw("psame");
    }
    
    if (fDisplay->GetSplitPadRow()){
	fDisplay->GetCanvasPadRow()->cd(2);
	fDisplay->GetPadPointer()->fHistpad2->Draw();
    }
    
    fDisplay->GetCanvasPadRow()->Modified();
    fDisplay->GetCanvasPadRow()->Update();
    
    // Select Pad
    fDisplay->GetCanvasPadRow()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliHLTTPCDisplayMain",(void*) fDisplay,"ExecPadEvent(Int_t,Int_t,Int_t,TObject*)");
    // Keep Zoom
    fDisplay->GetCanvasPadRow()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliHLTTPCDisplayPadRow",(void*) this,"ExecEvent(Int_t,Int_t,Int_t,TObject*)");
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPadRow::Draw3D(){
    Int_t markercolor = 51;
    
    for (UInt_t ii=0;ii < 20;ii++){
	if (fcolorbin[ii]> 0){
	    
	    TPolyMarker3D* pm = new  TPolyMarker3D(fcolorbin[ii], fpmarr[ii], 7 );
	    pm->SetBit(kCanDelete);

	    pm->SetMarkerColor(markercolor); 
	    pm->Draw(""); 
	}

	// in order to have the SetPalette(1), so called "pretty"
	if (ii % 2 == 0 ) markercolor += 2;
	else  markercolor += 3;
    }
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPadRow::ExecEvent(Int_t event, Int_t px, Int_t py, TObject *selected){
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
       
       if (selected == fHistraw->GetXaxis() || selected == fHistrawcl->GetXaxis()){ 
	   fBinX[0] = axis->GetFirst() -1;
	   fBinX[1] = axis->GetLast() -1;
       }
       else {
	   fBinY[0] = axis->GetFirst() -1;
	   fBinY[1] = axis->GetLast() -1;
       }

       
       //    Reset();
//     Fill();
       //   Draw();

       fTmpEvent = 0; 
   }
   else fTmpEvent = 0;
}
