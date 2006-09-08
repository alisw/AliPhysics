// $Id$

/** \class AliHLTTPCDisplayPadRow
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayResiduals
//
// Display class for the HLT TPC-Residuals events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include "AliHLTTPCDisplayResiduals.h"
#include "AliHLTTPCDisplayPadRow.h"

#include "AliHLTStdIncludes.h"
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAttText.h>
#include <TAxis.h>
#include <TCanvas.h>


#ifdef use_aliroot
#include <TClonesArray.h>
#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPCParam.h>
#endif

#include "AliHLTTPCDefinitions.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"


#include "AliHLTTPCDigitReader.h"
#include "AliHLT_C_Component_WrapperInterface.h"

#include "AliHLTTPCDisplayMain.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCDisplay.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDigitReaderPacked.h"


#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCDisplayResiduals)

//____________________________________________________________________________________________________
AliHLTTPCDisplayResiduals::AliHLTTPCDisplayResiduals(AliHLTTPCDisplayMain* display) {
    // constructor
    fDisplay = display;

    fGraphresidualsY = NULL;
    fGraphresidualsZ = NULL;
    fGraphresidualsYLength = NULL;

    fHistallresidualsY = new TH1F ("fHistallresidualsY","Y Residuals of all tracks in selected slices;residuals",5000,-100,100);
    fHistallresidualsZ = new TH1F ("fHistallresidualsZ","Z Residuals of all tracks in selected slices;residuals",5000,-100,100);
    fHistHits_S = new TH1F ("fHistHits_S","Number of cluster over track length;#Hits/s",20,0,20);
    fHistQ_Track = new TH1F ("fHistQ_Track","Average cluster charge per track;<q>/track",500,0,500);
    fHistQ_S   = new TH1F ("fHistQ_S","Total cluster charge over track length;Q/s",500,0,500);

    fHistallresidualsY->SetTitleSize(0.03);
    fHistallresidualsY->GetXaxis()->SetLabelSize(0.03);
    fHistallresidualsY->GetXaxis()->SetTitleSize(0.03);
    fHistallresidualsY->GetYaxis()->SetLabelSize(0.03);
    fHistallresidualsY->GetYaxis()->SetTitleSize(0.03);

    fHistallresidualsZ->SetTitleSize(0.03);
    fHistallresidualsZ->GetXaxis()->SetLabelSize(0.03);
    fHistallresidualsZ->GetXaxis()->SetTitleSize(0.03);
    fHistallresidualsZ->GetYaxis()->SetLabelSize(0.03);
    fHistallresidualsZ->GetYaxis()->SetTitleSize(0.03);

    fHistHits_S->SetTitleSize(0.03);
    fHistHits_S->GetXaxis()->SetLabelSize(0.03);
    fHistHits_S->GetXaxis()->SetTitleSize(0.03);
    fHistHits_S->GetYaxis()->SetLabelSize(0.03);
    fHistHits_S->GetYaxis()->SetTitleSize(0.03);

    fHistQ_Track->SetTitleSize(0.03);
    fHistQ_Track->GetXaxis()->SetLabelSize(0.03);
    fHistQ_Track->GetXaxis()->SetTitleSize(0.03);
    fHistQ_Track->GetYaxis()->SetLabelSize(0.03);
    fHistQ_Track->GetYaxis()->SetTitleSize(0.03);

    fHistQ_S->SetTitleSize(0.03);
    fHistQ_S->GetXaxis()->SetLabelSize(0.03);
    fHistQ_S->GetXaxis()->SetTitleSize(0.03);
    fHistQ_S->GetYaxis()->SetLabelSize(0.03);
    fHistQ_S->GetYaxis()->SetTitleSize(0.03);
}

//____________________________________________________________________________________________________
AliHLTTPCDisplayResiduals::~AliHLTTPCDisplayResiduals() {
    // destructor   
    if ( fHistallresidualsY){
	delete fHistallresidualsY;
	fHistallresidualsY = NULL;
    }
    fHistallresidualsY = NULL;
    if ( fHistallresidualsZ){
	delete fHistallresidualsZ;
	fHistallresidualsZ = NULL;
    }
    if ( fHistHits_S){
	delete fHistHits_S;
	fHistHits_S = NULL;
    }
    if ( fHistQ_Track){
	delete fHistQ_Track;
	fHistQ_Track = NULL;
    }
    if ( fHistQ_S){
	delete fHistQ_S;
	fHistQ_S = NULL;
    }
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayResiduals::Reset(){    
    fHistallresidualsY->Reset();
    fHistallresidualsZ->Reset();
    fHistHits_S->Reset();              
    fHistQ_Track->Reset();           
    fHistQ_S->Reset();

    fDisplay->GetCanvasResiduals()->Clear();
    fDisplay->GetCanvasHits_S()->Clear();
    fDisplay->GetCanvasQ_Track()->Clear();
    fDisplay->GetCanvasQ_S()->Clear();
}    

//____________________________________________________________________________________________________
void AliHLTTPCDisplayResiduals::Save(){  
    fDisplay->GetCanvasResiduals()->SaveAs("HLT-ResidualsView.eps"); 
    fDisplay->GetCanvasHits_S()->SaveAs("HLT-Hits_S.eps");
    fDisplay->GetCanvasQ_Track()->SaveAs("HLT-Q_Track.eps");
    fDisplay->GetCanvasQ_S()->SaveAs("HLT-Q_S.eps");
}


//____________________________________________________________________________________________________
void AliHLTTPCDisplayResiduals::Fill(){
    // Fill Resiudals Histograms / Graphs
   
    AliHLTTPCTrackArray* tracks = fDisplay->GetTrackArrayPointer();
    Int_t ntracks = tracks->GetNTracks();

    Double_t maxResidualY = 0.;
    Double_t maxResidualZ = 0.;

    fGraphresidualsY = NULL;
    fGraphresidualsYLength = NULL;
    fGraphresidualsZ = NULL;

    for(Int_t j=0; j<ntracks; j++) { 	
	
	AliHLTTPCTrack *gtrack = tracks->GetCheckedTrack(j); 
	if(!gtrack) continue;	

	// ---------------------------------------------------------------------
	// ++ Select if slice should be displayed or not ( selection in DISPLAY )

	Int_t slice = gtrack->GetSector();
	if (!fDisplay->GetDisplaySlice(slice)) continue;
	
	// -------------------
	// Get track parameter

	Int_t nHits = gtrack->GetNHits();  // Number of associated hits to track
	Double_t radius = gtrack->GetRadius();      // radius
	Double_t kappa = gtrack->GetKappa();        // curvature = 1/R , signed
	Double_t lambda = atan( gtrack->GetTgl() ); // dipAngle lambda
	
	// -----------------------------
	// ++ Check for straightr tracks

	if (kappa == 0 && AliHLTTPCTransform::GetBFieldValue() > 0.) {
	    printf("-#- KAPPA == 0 -#-");
	    // continue;
	}

	// ------------------------------------
        // ++ Get first/last point of the track

	Double_t xyzL[3];      // lastpoint of track
	Double_t xyzF[3];      // firstpoint of track
	
	xyzF[0] = gtrack->GetFirstPointX();
	xyzF[1] = gtrack->GetFirstPointY();
	xyzF[2] = gtrack->GetFirstPointZ();
	
	xyzL[0] = gtrack->GetLastPointX();
	xyzL[1] = gtrack->GetLastPointY();
	xyzL[2] = gtrack->GetLastPointZ();

	// --------------------------
	// ++ Calculate length of the track
	
	Double_t s = 0.;       // length of the track
	if (  AliHLTTPCTransform::GetBFieldValue() == 0. || kappa == 0 ) 
	  s = sqrt ( (xyzL[0] - xyzF[0])*(xyzL[0] - xyzF[0]) + (xyzL[1] - xyzF[1])*(xyzL[1] - xyzF[1]) ); 
	else {
	  // Calculate the length of the track. If it is to flat in in s,z plane use sxy, otherwise use sz
	  if (fabs(lambda) > 0.05){
	    // length of track calculated out of z
	    s = fabs( (xyzL[2] - xyzF[2]) / sin(lambda) ); // length of track calculated out of z
	  }
	  else {
	    Double_t d = (xyzL[0] - xyzF[0])*(xyzL[0] - xyzF[0]) + (xyzL[1] - xyzF[1])*(xyzL[1] - xyzF[1]); 
	    // length of track calculated out of xy
	    s = fabs ( acos( 0.5 * (2 - (d / (radius*radius)))) / ( kappa * cos(lambda) ) ); 		
	  }
	}
	
	// -----------------------
	// ++ Apply cuts on tracks

	if (nHits < fDisplay->GetCutHits() ) continue;
	if (s < fDisplay->GetCutS() ) continue;
	if (gtrack->GetPsi() < fDisplay->GetCutPsi() ) continue;
	if (lambda < fDisplay->GetCutLambda() ) continue; 
	if (gtrack->GetPt() < fDisplay->GetCutPt()  &&  AliHLTTPCTransform::GetBFieldValue() != 0. ) continue;
	if ( AliHLTTPCTransform::GetPadRow((Float_t)xyzF[0]) >   fDisplay->GetIncidentPadrow() ) continue;

	// ============================================
	// ## ROTATED Track to local coordinates BEGINN

	gtrack->Rotate(slice,kTRUE);
	
	Int_t nRes = 0;                                            // number of resiudals
	Double_t totalQ = 0.;                                      // total charge of track
	UInt_t *hitnum = gtrack->GetHitNumbers();                  // hist per track
	
	Double_t *resY= new Double_t[nHits];                       // Y residuals of every hit
	Double_t *resZ= new Double_t[nHits];                       // Z residuals of every hit
	
	Double_t *resYLength= new Double_t[2*nHits];               // length of pad in y direction
	
	Double_t *padrows = new Double_t[nHits];                   
	Double_t *padrowsLength = new Double_t[2*nHits];

	// ---------------------
	// ++ Loop over all hits
	
	for(Int_t h=0; h<nHits; h++){
	    UInt_t id=hitnum[h];
	    Int_t patch = (id>>22) & 0x7;
	    UInt_t pos = id&0x3fffff; 
	    
	    AliHLTTPCSpacePointData *points = fDisplay->GetSpacePointDataPointer(slice,patch);
	    if (!points) continue;

	    Float_t xyzCtmp[3];    // cluster tmp
	    Float_t xyzTtmp[3];    // track tmp
	    
	    // -------------------------
	    // ++ Get coordinates of Hit
	    xyzCtmp[0] = points[pos].fX;
	    xyzCtmp[1] = points[pos].fY;
	    xyzCtmp[2] = points[pos].fZ;
	    totalQ += points[pos].fCharge;

	    // ---------------------------------
	    // ++ Get X Coordinate of the padrow

	    Int_t padrow = AliHLTTPCTransform::GetPadRow(points[pos].fX);
	    xyzTtmp[0] = gtrack->GetFirstPointX();
	    
	    // --------------------------------------
	    // ++ Check for CrossingPoint with Padrow
	    if(gtrack->GetCrossingPoint(padrow,xyzTtmp)) {
		
		// ----------------------
		// ++ Calculate Residuals

		Float_t deltaY = ( xyzCtmp[1] - xyzTtmp[1] );
		Float_t deltaZ = ( xyzCtmp[2] - xyzTtmp[2] );
		
		padrows[nRes] = (Double_t) padrow;
		resY[nRes] = (Double_t) deltaY;
		resZ[nRes] = (Double_t) deltaZ;
		
		resYLength[(2*nRes)] = 0.5 * AliHLTTPCTransform::GetPadLength(padrow);
		resYLength[(2*nRes)+1] = -0.5 * AliHLTTPCTransform::GetPadLength(padrow);
		padrowsLength[(2*nRes)] = (Double_t) padrow;
		padrowsLength[(2*nRes)+1] = (Double_t) padrow;
		
		// ---------------------------
		// ++ Fill residuals histogram

		fHistallresidualsY->Fill(resY[nRes]);
		fHistallresidualsZ->Fill(resZ[nRes]);

		if (resY[nRes] > maxResidualY ) maxResidualY = resY[nRes];
		if (resZ[nRes] > maxResidualZ ) maxResidualZ = resZ[nRes];
		nRes++;

	    } // END CrossingPoint
	} // END cluster loop

	gtrack->Rotate(slice,kFALSE);

	// ## ROTATED Track to local coordinates END
 	// =========================================

	// -------------------------------------
	// ++ Fill Number Hits over track length

	Double_t hits_S = nHits / s;
	fHistHits_S->Fill(hits_S);

	// --------------------------------
	// ++ Fill Average charge per track

	Double_t q_Track = totalQ / nHits;
	fHistQ_Track->Fill(q_Track);

	// -------------------------------------
	// ++ Fill total charge per track length

	Double_t Q_S = totalQ / s;
	fHistQ_S->Fill(Q_S);

	// --------------------------
	// ++ Fill Graphs for 1 track

	if (fDisplay->GetSelectTrackSwitch() && fDisplay->GetGlobalTrack(slice) == j){
	    // FILL Y RESIDUALS GRAPH
	    fGraphresidualsY = new TGraph(nRes-1,padrows,resY);
	    fGraphresidualsYLength = new TGraph((2*nRes)-2,padrowsLength,resYLength);
	    // FILL Z RESIDUALS GRAPH
	    fGraphresidualsZ = new TGraph(nRes-1,padrows,resZ);
	}
	
	// --------------
	// ++ Free Memory
	
	if ( resY ){
	    delete [] resY;
	    resY = NULL;
	}
	if ( resZ ){
	    delete [] resZ;
	    resZ = NULL;
	}
	if ( resYLength ){
	    delete [] resYLength;
	    resYLength = NULL;
	}
	if ( padrows ){
	    delete [] padrows;
	    padrows = NULL;
	}
	if ( padrowsLength ){
	    delete [] padrowsLength;
	    padrowsLength = NULL;
	}
	
    } // END track loop

    // ----------------------------------------
    // ++ Set Axis Range of residual histograms

    fHistallresidualsY->SetAxisRange(-maxResidualY,maxResidualY);
    fHistallresidualsZ->SetAxisRange(-maxResidualZ,maxResidualZ);
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayResiduals::Draw(){

    fDisplay->GetCanvasResiduals()->cd();
    fDisplay->GetCanvasResiduals()->Clear();
    
    fDisplay->GetCanvasResiduals()->Divide(1,2);
   
    // ++ Y Residuals     
    // ==============
    fDisplay->GetCanvasResiduals()->cd(1);
   
    if (fDisplay->GetSelectTrackSwitch() && fGraphresidualsY){

	// -------------------------
	// ++ Y Graph for single track
	
	Char_t title[256];
	sprintf(title,"Y Residuals of Track %d in Slice %d",fDisplay->GetSelectTrack(), fDisplay->GetSelectTrackSlice() );

	TMultiGraph *mgY = new TMultiGraph() ;
	mgY->SetBit(kCanDelete);
	
	fGraphresidualsY->GetXaxis()->SetTitle("padrow");	
	fGraphresidualsY->GetYaxis()->SetTitle("residuals");
	fGraphresidualsY->GetXaxis()->SetLabelSize(0.02);
	fGraphresidualsY->GetXaxis()->SetTitleSize(0.02);
	fGraphresidualsY->GetYaxis()->SetLabelSize(0.02);
	fGraphresidualsY->GetYaxis()->SetTitleSize(0.02);

	fGraphresidualsYLength->SetMarkerColor(2);
	fGraphresidualsYLength->SetMarkerStyle(5);
	fGraphresidualsY->SetMarkerColor(1);
	fGraphresidualsY->SetMarkerStyle(3);

	mgY->Add(fGraphresidualsY);
	mgY->Add(fGraphresidualsYLength);
	mgY->SetTitle(title);
	mgY->Draw("AP");
    }

    else{

	// -------------------------
	// ++ Y Histogram  -> global
	
	fHistallresidualsY->SetStats(kFALSE);
	fHistallresidualsY->Draw();
    }
    
    // ++ Z Residuals     
    // ==============
    fDisplay->GetCanvasResiduals()->cd(2);

    // graph for single track
    if (fDisplay->GetSelectTrackSwitch() && fGraphresidualsZ){

	// -------------------------
	// ++ Z Graph for single track

	Char_t title[256];
	sprintf(title,"Z Residuals of Track %d in Slice %d",fDisplay->GetSelectTrack(), fDisplay->GetSelectTrackSlice() );
	
	TMultiGraph *mgZ = new TMultiGraph();
	mgZ->SetBit(kCanDelete);

	fGraphresidualsZ->SetTitle(title);	
	fGraphresidualsZ->GetXaxis()->SetTitle("padrow");	
	fGraphresidualsZ->GetYaxis()->SetTitle("residuals");
	fGraphresidualsZ->GetXaxis()->SetLabelSize(0.02);
	fGraphresidualsZ->GetXaxis()->SetTitleSize(0.02);
	fGraphresidualsZ->GetYaxis()->SetLabelSize(0.02);
	fGraphresidualsZ->GetYaxis()->SetTitleSize(0.02);
	
	mgZ->Add(fGraphresidualsZ);
	mgZ->SetTitle(title);
	mgZ->Draw("A*");
    }
    
	// -------------------------
	// ++ Z Histogram  -> global

    else{
	fHistallresidualsZ->SetStats(kFALSE);
	fHistallresidualsZ->Draw();
    }

    fDisplay->GetCanvasResiduals()->Modified();
    fDisplay->GetCanvasResiduals()->Update();

    // -------------------------------------
    // ++ Draw Number Hits over track length

    fDisplay->GetCanvasHits_S()->cd ();
    fDisplay->GetCanvasHits_S()->Clear();
    
    fHistHits_S->Draw();

    fDisplay->GetCanvasHits_S()->Modified();
    fDisplay->GetCanvasHits_S()->Update();

    // --------------------------------
    // ++ Draw Average charge per track

    fDisplay->GetCanvasQ_Track()->cd();
    fDisplay->GetCanvasQ_Track()->Clear();

    fHistQ_Track->Draw();

    fDisplay->GetCanvasQ_Track()->Modified();
    fDisplay->GetCanvasQ_Track()->Update();

    // -------------------------------------
    // ++ Draw total charge per track length

    fDisplay->GetCanvasQ_S()->cd();
    fDisplay->GetCanvasQ_S()->Clear();
    
    fHistQ_S->Draw();
  
    fDisplay->GetCanvasQ_S()->Modified();
    fDisplay->GetCanvasQ_S()->Update();
}

