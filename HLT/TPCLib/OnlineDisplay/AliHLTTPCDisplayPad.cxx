// $Id$

/** \class AliHLTTPCDisplayPad
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayPad
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
#include "AliHLTTPCDisplayPad.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCDisplayPad)

//____________________________________________________________________________________________________
AliHLTTPCDisplayPad::AliHLTTPCDisplayPad(AliHLTTPCDisplayMain* display) {
    // constructor
    fDisplay = display;

    fNTimes = AliHLTTPCTransform::GetNTimeBins();

    fBinX[0] = 0;      
    fBinX[1] = fNTimes -1; 
    fTmpEvent = 0;    

    fHistpad1 = new TH1F ("fHistpad1","Selected Pad -1;Timebin #",fNTimes,0,fNTimes-1);
    fHistpad2 = new TH1F ("fHistpad2","Selected Pad;Timebin #",fNTimes,0,fNTimes-1); 
    fHistpad3 = new TH1F ("fHistpad3","Selected Pad +1;Timebin #",fNTimes,0,fNTimes-1);

    fHistpad1->SetTitleSize(0.03);
    fHistpad1->GetXaxis()->SetLabelSize(0.03);
    fHistpad1->GetXaxis()->SetTitleSize(0.03);
    fHistpad1->GetYaxis()->SetLabelSize(0.03);
    fHistpad1->GetYaxis()->SetTitleSize(0.03);

    fHistpad2->SetTitleSize(0.03);
    fHistpad2->GetXaxis()->SetLabelSize(0.03);
    fHistpad2->GetXaxis()->SetTitleSize(0.03);
    fHistpad2->GetYaxis()->SetLabelSize(0.03);
    fHistpad2->GetYaxis()->SetTitleSize(0.03);

    fHistpad3->SetTitleSize(0.03);
    fHistpad3->GetXaxis()->SetLabelSize(0.03);
    fHistpad3->GetXaxis()->SetTitleSize(0.03);
    fHistpad3->GetYaxis()->SetLabelSize(0.03);
    fHistpad3->GetYaxis()->SetTitleSize(0.03);
}

//____________________________________________________________________________________________________
AliHLTTPCDisplayPad::~AliHLTTPCDisplayPad() {
    // destructor   
    if ( fHistpad1 ){
        delete fHistpad1;
        fHistpad1 = NULL;
    }

    if ( fHistpad2 ){
        delete fHistpad2;
        fHistpad2 = NULL;
    }

    if ( fHistpad3 ){
        delete fHistpad3;
        fHistpad3 = NULL;
    }  
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPad::Reset(){  
    fHistpad1->Reset(); 
    fHistpad2->Reset();  
    fHistpad3->Reset(); 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPad::Save(){  
    fDisplay->GetCanvasPad()->SaveAs("HLT-PadView.eps"); 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPad::Fill(Int_t patch, ULong_t dataBlock, ULong_t dataLen){
  // Fill Pad Histogram

  Int_t padRow = fDisplay->GetPadRow();
  Int_t pad = fDisplay->GetPad();
  Int_t nextPad = pad + 1;
  Int_t prevPad = pad - 1;
  
  // !!
  // !!  DO unrolling because of cache effects (avoid cache trashing) !!
  // !!

  if ( fDisplay->GetZeroSuppression() ){
    for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){
	fHistpad1->Fill( timeBin, fDisplay->fRawDataZeroSuppressed[padRow][prevPad][timeBin] );
    } // end time
    
    for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){
      fHistpad2->Fill( timeBin, fDisplay->fRawDataZeroSuppressed[padRow][pad][timeBin] );
    } // end time
    
    for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){
      fHistpad3->Fill( timeBin, fDisplay->fRawDataZeroSuppressed[padRow][nextPad][timeBin] );
    } // end time
    
  }  // end - if ( fDisplay->GetZeroSuppression() ){
  
  else {
    for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){
	fHistpad1->Fill( timeBin, fDisplay->fRawData[padRow][prevPad][timeBin] );
    } // end time
    
    for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){
      fHistpad2->Fill( timeBin, fDisplay->fRawData[padRow][pad][timeBin] );
      } // end time
    
    for (Int_t timeBin=0; timeBin <= fDisplay->GetNTimeBins(); timeBin++){
      fHistpad3->Fill( timeBin, fDisplay->fRawData[padRow][nextPad][timeBin] );
    } // end time
    
  }  // end - else of if ( fDisplay->GetZeroSuppression() ){
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPad::Draw(){
    Char_t title[256];

    fDisplay->GetCanvasPad()->cd();
    fDisplay->GetCanvasPad()->Clear();
    fDisplay->GetCanvasPad()->Divide(1,3);
    
    // Keep Zoom
    if (!fDisplay->GetKeepView() || (fBinX[1]>fNTimes)){
	fBinX[0] = 0;
	fBinX[1] = fNTimes -1;
    }
    
    fHistpad1->SetAxisRange(fBinX[0],fBinX[1]);
    fHistpad2->SetAxisRange(fBinX[0],fBinX[1]);
    fHistpad3->SetAxisRange(fBinX[0],fBinX[1]);


    fDisplay->GetCanvasPad()->cd(1);
    sprintf(title,"Selected Pad %d",fDisplay->GetPad() -1);
    fHistpad1->SetStats(kFALSE);
    fHistpad1->SetTitle(title);
    fHistpad1->Draw();

    fDisplay->GetCanvasPad()->cd(2); 
    sprintf(title,"Selected Pad %d",fDisplay->GetPad() );
    fHistpad2->SetStats(kFALSE);
    fHistpad2->SetTitle(title);
    fHistpad2->Draw();

    fDisplay->GetCanvasPad()->cd(3);
    sprintf(title,"Selected Pad %d",fDisplay->GetPad() +1);
    fHistpad3->SetStats(kFALSE);
    fHistpad3->SetTitle(title);
    fHistpad3->Draw();

    fDisplay->GetCanvasPad()->Modified();
    fDisplay->GetCanvasPad()->Update();

    // Keep Zoom
    fDisplay->GetCanvasPad()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliHLTTPCDisplayPad",(void*) this,"ExecEvent(Int_t,Int_t,Int_t,TObject*)");
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayPad::ExecEvent(Int_t event, Int_t px, Int_t py, TObject *selected){
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
       
       if (selected == fHistpad1->GetXaxis() || selected == fHistpad2->GetXaxis() || selected == fHistpad3->GetXaxis()){ 
	   fBinX[0] = axis->GetFirst() -1;
	   fBinX[1] = axis->GetLast() -1;
       }

       fTmpEvent = 0; 
   }
   else fTmpEvent = 0;
}
