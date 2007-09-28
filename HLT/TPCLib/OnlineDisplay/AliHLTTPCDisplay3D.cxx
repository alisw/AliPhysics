// $Id$

/**
<pre>
//_____________________________________________________________
// AliHLTTPCDisplay3D
//
// Display class for the HLT TPC-3D events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#define TRACKHELIX 0       // use THelix for tracks
#define TRACKPOLYMARKER 0  // use TPolymarker3D for tracks
#define FIRSTLASTPOINT 0   // show first / last point of tracks
#define DEBUG 1

#define DRAWSTEP 0.2

#define UNUSEDCLUSTERCOLOR 2
#define USEDCLUSTERCOLOR 3
#define TRACKCOLOR 4
#define TRACKPOLYMARKERCOLOR 5
#define TRACKHELIXCOLOR 6

#if defined(HAVE_HOMERREADER) 
#include HOMERREADER_HEADER
#endif // defined(HAVE_HOMERREADER) 

#include "AliHLTTPCDisplay3D.h"
#include "AliHLTTPCDisplayPadRow.h"

#include "AliHLTStdIncludes.h"
#if defined(HAVE_TVIEW3D_H)
#include <TView3D.h>
#else
#include <TView.h>
#endif
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TH2.h>
#include <TTree.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TShape.h>
#include <TParticle.h>
#include <TFile.h>
#include <THelix.h>
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
#include "AliHLTTPCDigitReaderRaw.h"
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

ClassImp(AliHLTTPCDisplay3D)

//____________________________________________________________________________________________________
AliHLTTPCDisplay3D::AliHLTTPCDisplay3D(AliHLTTPCDisplayMain* display, Char_t* gfile ) {
    // constructor
    fDisplay = display;
    
    fGeom = NULL;
    LoadGeometrie(gfile);
}

//____________________________________________________________________________________________________
AliHLTTPCDisplay3D::~AliHLTTPCDisplay3D() {
    // destructor   

}

//____________________________________________________________________________________________________
void AliHLTTPCDisplay3D::Save(){
  fDisplay->GetCanvas3D()->SaveAs("HLT-3D-View.eps");
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplay3D::Draw(){
    fDisplay->GetCanvas3D()->cd();
    fDisplay->GetCanvas3D()->Clear();

#if defined(HAVE_TVIEW3D_H)
    TView3D *v = new TView3D();
    if (v) v->SetSystem(1);
#else
    TView *v = new TView(1);
#endif
    if (v==NULL) {
      HLTFatal("can not create viewer");
      return;
    }
    v->SetRange(-800,-800,-800,800,800,800);

    Float_t* etaRange = NULL;   // ------  STILL TO FIX


    //--------------------------------------------------------------------------------------------
    // DRAW 3D GEOMETRY
    //--------------------------------------------------------------------------------------------
    if (fDisplay->Get3DSwitchGeometry()){

	TList* masterNodeList = fGeom->GetListOfNodes();
	TNode* masterNode=0;
	TIter next(masterNodeList);
	
   	while ((masterNode = static_cast<TNode*> (next()))) {  
	
	    TList* nodeList = masterNode->GetListOfNodes();
	    TNode* node=0;
	    TIter next(nodeList);
	    
	    while ((node = static_cast<TNode*> (next()))) {  
		
		ULong_t tmpslice = atol(node->GetName() + 2);

		if (fDisplay->GetDisplaySlice(tmpslice)) {
		    node->SetVisibility(1);
		    node->SetFillColor(0);
		    node->SetLineColor(fDisplay->GetLineColor());
		}
		else node->SetVisibility(0);
	    } // end while son Nodes
    	} // end while master Node

	fGeom->Draw("");
    }   // END - DRAW 3D GEOMETRY

    //--------------------------------------------------------------------------------------------
    // DRAW 3D PadRow
    //--------------------------------------------------------------------------------------------
    
    if (fDisplay->ExistsRawData() &&  fDisplay->Get3DSwitchPadRow() && fDisplay->GetDisplaySlice(fDisplay->GetSlicePadRow())){

      // -- only one padrow
      if ( fDisplay->Get3DSwitchRaw() == 1 ) {
      	fDisplay->GetPadRowPointer()->Draw3D();
      }
      // show all padrows 
      else {
	Int_t nRawdata = 0;
	TPolyMarker3D *pmRawdata = new TPolyMarker3D(1,6);
	pmRawdata->SetBit(kCanDelete);
	
	if ( fDisplay->GetZeroSuppression() ){
	  for (Int_t row=0; row < AliHLTTPCTransform::GetNRows(); row++){
	    for (Int_t pad=0; pad < AliHLTTPCTransform::GetNPads(row); pad++){
	      for (Int_t timeBin=fDisplay->GetTimeBinMin(); timeBin <= fDisplay->GetTimeBinMax(); timeBin++){
		
		Int_t time = (Int_t) fDisplay->fRawDataZeroSuppressed[row][pad][timeBin];

		if (time != 0) {
		  Float_t xyz[3];    
		  
		  AliHLTTPCTransform::RawHLT2Global(xyz,fDisplay->GetSlicePadRow(), row, pad, timeBin);
		  
		  pmRawdata->SetPoint(nRawdata,xyz[0],xyz[1],xyz[2]);
		  nRawdata++;
		}
	      } // end time
	    } // end pad
	  } // end row
	}

	pmRawdata->SetMarkerSize(1);
	pmRawdata->SetMarkerColor(6); 
	pmRawdata->Draw("same");

      } // end show all padrows

    }

    //--------------------------------------------------------------------------------------------
    // DRAW 3D CLUSTER
    //--------------------------------------------------------------------------------------------
    if (fDisplay->Get3DSwitchCluster() && fDisplay->ExistsClusterData()){

	for (Int_t slice=0; slice <= 35; slice++){

	    Int_t currenttrack = -1;

	    if (fDisplay->ExistsTrackData() && fDisplay->GetSelectTrackSwitch()) currenttrack = fDisplay->GetGlobalTrack(slice);

	    if (!fDisplay->GetDisplaySlice(slice)) continue;
	    
	    for(Int_t patch=0;patch<6;patch++){

		AliHLTTPCSpacePointData *points = fDisplay->GetSpacePointDataPointer(slice,patch);
		if(!points) continue;

		TPolyMarker3D *pmUsed = new TPolyMarker3D(1,6);
		TPolyMarker3D *pmUnused = new TPolyMarker3D(1,6);
		pmUnused->SetBit(kCanDelete);
		pmUsed->SetBit(kCanDelete);


		Int_t nUsedCluster = 0;
		Int_t nUnusedCluster = 0;

		Float_t xyz[3];    
		for(Int_t i=0; i< fDisplay->GetNumberSpacePoints(slice,patch); i++){
		    // Used  cluster only
		    if (fDisplay->GetSelectCluster() == 1  && points[i].fUsed == kFALSE) continue; 
		    // Unused cluster only
		    if (fDisplay->GetSelectCluster() == 2  && points[i].fUsed == kTRUE) continue; 
		    
		    // if single track is selcted draw only cluster for this track
		    if (fDisplay->GetSelectCluster() == 1 && fDisplay->GetSelectTrackSwitch() && points[i].fTrackN != currenttrack) continue;
		    
		    xyz[0] = points[i].fX;
		    xyz[1] = points[i].fY;
		    xyz[2] = points[i].fZ;
		    
		    if ( etaRange ){		  
			// Do this before the transform, because the tracker also uses
			// local coordinates when using this limit to determine 
			// which clusters to use for tracking
			Double_t pointEta = AliHLTTPCTransform::GetEta( xyz );
			if ( pointEta<etaRange[0] || pointEta>etaRange[1] )
			    continue;
		    }
		    
		    AliHLTTPCTransform::Local2Global(xyz,slice);
		    
		    if (points[i].fUsed == kTRUE){
			pmUsed->SetPoint(nUsedCluster,xyz[0],xyz[1],xyz[2]);
			nUsedCluster++;
		    }
		    else {
			pmUnused->SetPoint(nUnusedCluster,xyz[0],xyz[1],xyz[2]);
			nUnusedCluster++;
		    }
		    
		}
		pmUsed->SetMarkerSize(1);
		pmUsed->SetMarkerColor(USEDCLUSTERCOLOR); 
		pmUsed->Draw("same");

		pmUnused->SetMarkerSize(1);
		pmUnused->SetMarkerColor(UNUSEDCLUSTERCOLOR); 
		pmUnused->Draw("same");

		fDisplay->GetCanvas3D()->Modified();
		fDisplay->GetCanvas3D()->Update();

	    } // END - PATCH LOOP	    
	}  // END - SLICE LOOP	
    }   // END - DRAW 3D CLUSTER 

    //--------------------------------------------------------------------------------------------
    // DRAW 3D TRACKS
    //--------------------------------------------------------------------------------------------
    if (fDisplay->Get3DSwitchTracks() && fDisplay->ExistsTrackData()){

      AliHLTTPCTransform::SetBField( 0.5 );  // ++++++


	AliHLTTPCTrackArray* tracks = fDisplay->GetTrackArrayPointer();
	Int_t ntracks = tracks->GetNTracks();

	//	TPolyLine3D **line = new (TPolyLine3D*)[ntracks];
	//	for(Int_t j=0; j<ntracks; j++) line[j] = 0;
#if TRACKHELIX
	//	THelix **helix = new (THelix*)[ntracks];
	//	for(Int_t j=0; j<ntracks; j++) helix[j] = 0;

#endif
	for(Int_t j=0; j<ntracks; j++) { 	

	    AliHLTTPCTrack *gtrack = tracks->GetCheckedTrack(j); 
	    if(!gtrack) continue;

	    Int_t nHits = gtrack->GetNHits();  // Number of associated hits to track
	    Int_t slice = gtrack->GetSector();

            // --- CHECK if track is should be drawn
	    // select if slice should be displayed or not
	    if (!fDisplay->GetDisplaySlice(slice)) continue;

	    if (fDisplay->GetSelectTrackSwitch() && fDisplay->GetGlobalTrack(slice) != j) continue;

	    Double_t radius = gtrack->GetRadius();      // radius
	    Double_t kappa = gtrack->GetKappa();        // curvature = 1/R , signed
	    Double_t lambda = atan( gtrack->GetTgl() ); // dipAngle lambda
	    Double_t phi0 = gtrack->GetPsi() + (gtrack->GetCharge() * AliHLTTPCTransform::PiHalf() ); // azimuthal angle of startingpoint, with respect to helix axis

	    if (kappa == 0 && AliHLTTPCTransform::GetBFieldValue() > 0.) {
		printf("================================KAPPA == 0");
		continue;
	    }

	    Double_t xyzL[3];      // lastpoint of track
	    Double_t xyzF[3];      // firstpoint of track

	    xyzF[0] = gtrack->GetFirstPointX();
	    xyzF[1] = gtrack->GetFirstPointY();
	    xyzF[2] = gtrack->GetFirstPointZ();

	    xyzL[0] = gtrack->GetLastPointX();
	    xyzL[1] = gtrack->GetLastPointY();
	    xyzL[2] = gtrack->GetLastPointZ();

	    if ( xyzL[0] ==  xyzF[0] ) continue;
	    if ( xyzL[1] ==  xyzF[1] ) continue;
	    if ( xyzL[2] ==  xyzF[2] ) continue;


#if FIRSTLASTPOINT	    
            TPolyMarker3D *pmL = new TPolyMarker3D(1,2);
	    TPolyMarker3D *pmF = new TPolyMarker3D(1,2);
	    pmL->SetBit(kCanDelete);
	    pmF->SetBit(kCanDelete);

	    //TPolyMarker3D pmL(1,2);
	    //TPolyMarker3D pmF(1,2);


	    pmF->SetPoint(0,xyzF[0],xyzF[1],xyzF[2]);
	    pmL->SetPoint(0,xyzL[0],xyzL[1],xyzL[2]);
#endif

	    Double_t s = 0.;       // length of the track


	    //if (  AliHLTTPCTransform::GetBFieldValue() == 0.) 
            if ( radius == 999999.0){ // IVAN HACK

		s = sqrt ( (xyzL[0] - xyzF[0])*(xyzL[0] - xyzF[0]) + (xyzL[1] - xyzF[1])*(xyzL[1] - xyzF[1]) ); 
	    }
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


	    // CUTS on tracks
/*
	    if (nHits < fDisplay->GetCutHits() ) continue;
	    if (s < fDisplay->GetCutS() ) continue;
	    if (gtrack->GetPsi() < fDisplay->GetCutPsi() ) continue;
	    if (lambda < fDisplay->GetCutLambda() ) continue; 
	    if (gtrack->GetPt() < fDisplay->GetCutPt()  &&  AliHLTTPCTransform::GetBFieldValue() != 0. ) continue;
	    if ( AliHLTTPCTransform::GetPadRow((Float_t)xyzF[0]) >   fDisplay->GetIncidentPadrow() ) continue;
*/    
	    Int_t nTrackPoints = 2 + (Int_t) floor(s / DRAWSTEP);
			    
#if TRACKPOLYMARKER
	    //    TPolyMarker3D *pmT = new TPolyMarker3D(nTrackPoints,6);
	    TPolyMarker3D pmT(nTrackPoints,6);

#endif
	    Double_t *xT = new Double_t[nTrackPoints];
	    Double_t *yT = new Double_t[nTrackPoints];
	    Double_t *zT = new Double_t[nTrackPoints];

	    Int_t trackPointCounter = 0;

	    //Write Track Parameters for single track
	    if (fDisplay->GetSelectTrackSwitch() ){
		fDisplay->fTrackParam.id = j;
		fDisplay->fTrackParam.nHits = nHits;
		fDisplay->fTrackParam.charge = gtrack->GetCharge();
		fDisplay->fTrackParam.lambda = lambda;
		fDisplay->fTrackParam.kappa = kappa;
		fDisplay->fTrackParam.radius = radius;
		fDisplay->fTrackParam.slice = slice;
		fDisplay->fTrackParam.phi0 = phi0;
		fDisplay->fTrackParam.pt = gtrack->GetPt();
		fDisplay->fTrackParam.bfield = AliHLTTPCTransform::GetBFieldValue();
		fDisplay->fTrackParam.psi = gtrack->GetPsi();
		fDisplay->fTrackParam.s = s;
	    }


	    //  if (  AliHLTTPCTransform::GetBFieldValue() == 0.) {

		    if ( radius == 999999.0) { // IVAN HACK

		      //	      cout << "ENTER IVAN 2   radius="<< radius << endl;
		for (Double_t ds = 0.; ds < s; ds = ds + DRAWSTEP){
		    // FILL ARRAYS IN ORDER TO DRAW THE TRACKPOINTS, OUT OF THE PARAMETER
		    //		    xT[trackPointCounter] = xyzF[0] + ds * cos(phi0); 
		    //		    yT[trackPointCounter] = xyzF[1] + ds * sin(phi0);
		  xT[trackPointCounter] = xyzF[0] + ds * cos(phi0 - 1.57); 
		  yT[trackPointCounter] = xyzF[1] + ds * sin(phi0 - 1.57);
		  zT[trackPointCounter] = xyzF[2] + ds * sin(lambda);
#if TRACKPOLYMARKER
		    pmT.SetPoint(trackPointCounter,xT[trackPointCounter],yT[trackPointCounter],zT[trackPointCounter]);
#endif
		    trackPointCounter++;
		}
		
	        if (trackPointCounter > nTrackPoints) printf("N=%d  n=%d", nTrackPoints,trackPointCounter);
		else {
		  //xT[trackPointCounter] = xyzF[0] + s * cos(phi0); 
		  //yT[trackPointCounter] = xyzF[1] + s * sin(phi0);

		  xT[trackPointCounter] = xyzF[0] + s * cos(phi0 - 1.57); 
		  yT[trackPointCounter] = xyzF[1] + s * sin(phi0 - 1.57);
		  zT[trackPointCounter] = xyzF[2] + s * sin(lambda);
#if TRACKPOLYMARKER       
		    pmT.SetPoint(trackPointCounter,xT[trackPointCounter],yT[trackPointCounter],zT[trackPointCounter]);
#endif
		}

	    }
	    else {

		for (Double_t ds = 0.; ds < s; ds = ds + DRAWSTEP){
		    // FILL ARRAYS IN ORDER TO DRAW THE TRACKPOINTS, OUT OF THE PARAMETER
		    xT[trackPointCounter] = xyzF[0] + radius * ( cos( phi0 + (ds*kappa*cos(lambda)) ) - cos(phi0) );
		    yT[trackPointCounter] = xyzF[1] + radius * ( sin( phi0 + (ds*kappa*cos(lambda)) ) - sin(phi0) );
		    zT[trackPointCounter] = xyzF[2] + ds * sin(lambda);
#if TRACKPOLYMARKER
		    pmT.SetPoint(trackPointCounter,xT[trackPointCounter],yT[trackPointCounter],zT[trackPointCounter]);
#endif
		    trackPointCounter++;
		}
		
	        if (trackPointCounter > nTrackPoints) printf("N=%d  n=%d", nTrackPoints,trackPointCounter);
		else {
		    xT[trackPointCounter] = xyzF[0] + radius * ( cos( phi0 + (s*kappa*cos(lambda)) ) - cos(phi0) );
		    yT[trackPointCounter] = xyzF[1] + radius * ( sin( phi0 + (s*kappa*cos(lambda)) ) - sin(phi0) );
		    zT[trackPointCounter] = xyzF[2] + s * sin(lambda);
#if TRACKPOLYMARKER       
		    pmT.SetPoint(trackPointCounter,xT[trackPointCounter],yT[trackPointCounter],zT[trackPointCounter]);
#endif
		}
	    }

	    // Draw Track -- as line
	    //line[j] = new TPolyLine3D(nTrackPoints,xT,yT,zT,"");
	    //TPolyLine3D *currentline = line[j];
	    //* currentline = new TPolyLine3D(nTrackPoints,xT,yT,zT,"");
	    //	    TPolyLine3D currentline(nTrackPoints,xT,yT,zT,"");

	    TPolyLine3D *currentline = new TPolyLine3D(nTrackPoints,xT,yT,zT,"");
	    currentline->SetBit(kCanDelete);
	    currentline->SetLineColor(TRACKCOLOR);   
	    currentline->SetLineWidth(2);
	    currentline->Draw("same");



	    // ---  ADDITIONAL DRAW OPTIONS
#if FIRSTLASTPOINT
	    // Draw last point of Track
	    pmL->SetMarkerSize(3);
	    pmL->SetMarkerColor(4); 
	    pmL->Draw();

	    // Draw first point of Track
	    pmF->SetMarkerSize(3);
	    pmF->SetMarkerColor(5); 
	    pmF->Draw();
#endif
#if TRACKPOLYMARKER
	    // Draw Track -- as polymarker
	    pmT.SetMarkerSize(3);
	    pmT.SetMarkerColor(TRACKPOLYMARKERCOLOR); 
	    pmT.Draw();
#endif
#if TRACKHELIX
	    // Draw Track -- as helix
            // works ok, execpt for very small dipangles -> track almost horizontal
	    Double_t hrange[2];
	    Double_t v0[3];
	    Double_t omega;
	    hrange[0] = xyzF[2];
	    hrange[1] = xyzL[2];
	    v0[0] = gtrack->GetPx();
	    v0[1] = gtrack->GetPy();
	    v0[2] = gtrack->GetPz();
	    omega = AliHLTTPCTransform::GetBFieldValue() * gtrack->GetCharge();

	    //	    helix[j] = new THelix(xyzF,v0,omega,hrange,kHelixZ,0);
	    //	    THelix *currenthelix = helix[j];
	    //	    currenthelix = new THelix(xyzF,v0,omega,hrange,kHelixZ,0);
	    THelix currenthelix(xyzF,v0,omega,hrange,kHelixZ,0);
	    currenthelix.SetLineColor(TRACKHELIXCOLOR);   
	    currenthelix.SetLineWidth(1);
	    currenthelix.Draw("same");   	    
#endif


	    // delete[]
	    //Double_t *xT = new Double_t[nTrackPoints];
	    //	    Double_t *yT = new Double_t[nTrackPoints];
	    //	    Double_t *zT = new Double_t[nTrackPoints];
	    if (xT){ 
	      delete[] xT; 
	      xT = NULL;
	    }
	    if (yT){ 
	      delete[] yT; 
	      yT = NULL;
	    }
	    if (zT){ 
	      delete[] zT; 
	      zT = NULL;
	    }

      	} // END for track loop


	// NO !!! DELETE line #ifdef helix delete helix


    }   // END - DRAW 3D Tracks
 
    //--------------------------------------------------------------------------------------------
    // DRAW 3D 
    //--------------------------------------------------------------------------------------------
    // v->ZoomView(0,3);
    // v->Draw();   

    fDisplay->GetCanvas3D()->SetFillColor(fDisplay->GetBackColor());

    if ( !fDisplay->GetKeepView() ){ 
	fDisplay->GetCanvas3D()->SetTheta(fDisplay->GetTheta());
	fDisplay->GetCanvas3D()->SetPhi(fDisplay->GetPhi());
    }

    // fDisplay->GetCanvas3D()->Modified();

    fDisplay->GetCanvas3D()->Update();

}

//____________________________________________________________________________________________________
void AliHLTTPCDisplay3D::LoadGeometrie(Char_t *gfile) {
    if (gfile) {
	TFile *file = TFile::Open(gfile);
	if(!file) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay3D::AliHLTDisplay","File Open") 
	      <<"Geometry file " << gfile << " does not exist!"<<ENDLOG;
	    exit(-1);
	}
	
	fGeom = (TGeometry*)file->Get("AliceGeom");

	file->Close();
    }
}
