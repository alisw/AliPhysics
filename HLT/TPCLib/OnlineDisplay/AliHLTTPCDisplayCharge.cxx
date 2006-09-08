// $Id$

/** \class AliHLTTPCDisplayCharge
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayCharge
//
// Display class for the HLT TPC-Charge events.
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
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"

#include "AliHLTTPCDisplayMain.h"
#include "AliHLTTPCDisplayCharge.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCDisplayCharge)

//____________________________________________________________________________________________________
AliHLTTPCDisplayCharge::AliHLTTPCDisplayCharge(AliHLTTPCDisplayMain* display) {
    // constructor
    fDisplay = display;

    fBinX[0] = 0;      
    fBinX[1] = 1; 
    fTmpEvent = 0;         
    fMaxCharge = 0;

    fHistcharge = new TH1F ("fHistcharge","Cluster distribution per charge;charge;#cluster",5000,0,30000);
    fHistcharge->SetTitleSize(0.03);
    fHistcharge->GetXaxis()->SetLabelSize(0.03);
    fHistcharge->GetXaxis()->SetTitleSize(0.03);
    fHistcharge->GetYaxis()->SetLabelSize(0.03);
    fHistcharge->GetYaxis()->SetTitleSize(0.03);
}

//____________________________________________________________________________________________________
AliHLTTPCDisplayCharge::~AliHLTTPCDisplayCharge() {
    // destructor   
    if ( fHistcharge){
	delete fHistcharge;
	fHistcharge = NULL;
    }
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayCharge::Reset(){  
    fHistcharge->Reset();
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayCharge::Save(){  
  fDisplay->GetCanvasCharge()->SaveAs("HLT-ChargeView.eps"); 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayCharge::Fill(){  
    // Fill Charge Histogram

    Int_t maxCharge = 0;

    for (Int_t slice=0; slice <= 35; slice++){

	Int_t currenttrack = -1;

	if (fDisplay->ExistsTrackData() && fDisplay->GetSelectTrackSwitch()) currenttrack = fDisplay->GetGlobalTrack(slice);
		
	if (!fDisplay->GetDisplaySlice(slice)) continue;
	
	for(Int_t patch=0;patch<6;patch++){

	    AliHLTTPCSpacePointData *points = fDisplay->GetSpacePointDataPointer(slice,patch);
	    if(!points) continue;
	  
	    for(Int_t i=0; i< fDisplay->GetNumberSpacePoints(slice,patch); i++){
		// Used  cluster only
		if (fDisplay->GetSelectCluster() == 1  && points[i].fUsed == kFALSE) continue; 
		// Unused cluster only
		if (fDisplay->GetSelectCluster() == 2  && points[i].fUsed == kTRUE) continue; 

		// if single track is selcted draw only cluster for this track
		if (fDisplay->GetSelectCluster() == 1 && fDisplay->GetSelectTrackSwitch() && points[i].fTrackN != currenttrack) continue;
		
		// Fill Charge Histogram
		fHistcharge->Fill(points[i].fCharge);
		if ((Int_t)points[i].fCharge > maxCharge ) maxCharge = (Int_t) points[i].fCharge; 
	    }

	} // END - PATCH LOOP	    
    }  // END - SLICE LOOP

    fMaxCharge = maxCharge;     
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayCharge::Draw(){
    fDisplay->GetCanvasCharge()->cd();
    fDisplay->GetCanvasCharge()->Clear();

    // Keep Zoom
    if (!fDisplay->GetKeepView() || (fBinX[1]>fMaxCharge)){
	fBinX[0] = 0;
	fBinX[1] = fMaxCharge;
    }

    fHistcharge->SetAxisRange(fBinX[0],fBinX[1]);

    fHistcharge->SetStats(kFALSE);
    fHistcharge->Draw();       

    fDisplay->GetCanvasCharge()->Modified();
    fDisplay->GetCanvasCharge()->Update();

    // Keep Zoom
    fDisplay->GetCanvasCharge()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliHLTTPCDisplayCharge",(void*) this,"ExecEvent(Int_t,Int_t,Int_t,TObject*)");
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayCharge::ExecEvent(Int_t event, Int_t px, Int_t py, TObject *selected){
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
       
       if (selected == fHistcharge->GetXaxis() ){
	   fBinX[0] = axis->GetFirst() -1;
	   fBinX[1] = axis->GetLast() -1;
       }

       fTmpEvent = 0; 
   }
   else fTmpEvent = 0;
}
