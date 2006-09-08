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
#include "AliHLTTPCDigitReaderPacked.h"
#include "AliHLTTPCDigitReaderRaw.h"

#include "AliHLTTPCDisplayMain.h"
#include "AliHLTTPCDisplayFront.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCDisplayFront)

//____________________________________________________________________________________________________
AliHLTTPCDisplayFront::AliHLTTPCDisplayFront(AliHLTTPCDisplayMain* display) {
    // constructor
    fDisplay = display;

    fNTimes = AliHLTTPCTransform::GetNTimeBins();

    fBinY[0] = 0;
    fBinY[1] = AliHLTTPCTransform::GetNRows() - 1;
    fBinX[0] = 0;      
    fBinX[1] = AliHLTTPCTransform::GetNPads(fBinY[1]);
    fTmpEvent = 0;    

    fHistfront = new TH2F("fHistfront","FrontView of selected slice;Pad #;Padrow #",fBinX[1]+1,fBinX[0],fBinX[1],fBinY[1]+1,fBinY[0],fBinY[1]);

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

#if defined(HAVE_TPC_MAPPING)
    AliHLTTPCDigitReaderRaw digitReader(0);

    bool readValue = true;
    Int_t rowOffset = 0;
    
    Int_t timebin = fDisplay->GetTimebin();
    Bool_t allTimebins = fDisplay->GetAllTimebins();
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
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplayFrontRow::Fill","Read first value") << "No value in data block" << ENDLOG;
	return;
    }

    // -- Fill Raw Data
    while ( readValue ){ 

	UShort_t time = digitReader.GetTime();

	// check if all timebins or just one
	if (allTimebins || time == timebin ) 
	    fHistfront->Fill(digitReader.GetPad(),(digitReader.GetRow() + rowOffset),digitReader.GetSignal());
	
	// read next value
	readValue = digitReader.Next();
      
	//Check where to stop:
	if(!readValue) break; //No more value
    } 
#else //! defined(HAVE_TPC_MAPPING)
      HLTFatal("DigitReaderRaw not available - check your build");
#endif //defined(HAVE_TPC_MAPPING)
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayFront::Draw(){
    fDisplay->GetCanvasFront()->cd();
    fDisplay->GetCanvasFront()->Clear();

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

    fDisplay->GetCanvasFront()->Modified();
    fDisplay->GetCanvasFront()->Update();

    // Select Pad
    // fCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliHLTTPCDisplayFront",(void*) fDisplay,"ExecPadEvent(Int_t,Int_t,Int_t,TObject*)");

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
