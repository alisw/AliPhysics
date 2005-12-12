// @(#) $Id$

/** \class AliHLTTPCDisplay
<pre>
//_____________________________________________________________
// AliHLTTPCDisplay
//
// Display class for the HLT TPC events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//         Anders Vestbo <mailto:vestbo@fi.uib.no>      
//*-- Copyright &copy ALICE HLT Group 


// Recent Changes:
// ==============
// - Rename and Merge of functions / Complete new arrangement in order to use the AliHLTGUI
// - 3D Geometry
//   - display padrows, cluster, tracks
//   - select single Tracks
//   - select used / unused cluster
// - Histogram
//   - display padrows
//   - display pads in padrows

#include "AliHLTTPCStandardIncludes.h"
#include <TView.h>
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
#include <TAttText.h>
#include <TAxis.h>


#ifdef use_aliroot
#include <TClonesArray.h>
#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPCParam.h>
#endif

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCDisplay.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDigitReaderPacked.h"

#if __GNUC__ == 3
using namespace std;
#endif

ClassImp(AliHLTTPCDisplay)

// #############################################################################
AliHLTTPCDisplay::AliHLTTPCDisplay(Char_t *gfile) {
    //constructor
    memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
    memset(fNcl, 0, 36*6*sizeof(UInt_t)); 

    fTracks = NULL;
    fHistrawcl = NULL;
    fHistraw = NULL;
    fHistpad1 = NULL;
    fHistpad2 = NULL;
    fHistpad3 = NULL;
    fHistallresiduals = NULL;
    fHistcharge = NULL;
    fGraphresiduals = NULL;

    fGeom = NULL;

    fNPads = 0;
    fNTimes = 0;
    fMinHits = 0;
    fPtThreshold = 0.;
    fPad = -1;
    fPadRow = 0;
    fSlicePadRow = 0; 
    fSelectTrack = -1;
    fSelectTrackSlice = 0;
    fSelectTrackSwitch = kFALSE;
    fSelectCluster = 0;

    fMinSlice = 0;
    fMaxSlice = 35;
    fSlicePair = kFALSE;

    SetSliceArray();

    fBackColor = 1; 
    fLineColor = 0;

    fSwitch3DCluster = kFALSE;
    fSwitch3DTracks = kFALSE;
    fSwitch3DPadRow = kFALSE;
    fSwitch3DGeometry = kFALSE;

    //ctor. Specify which slices you want to look at.
    LoadGeometrie(gfile);
}


// #############################################################################
AliHLTTPCDisplay::~AliHLTTPCDisplay() {
    //destructor
    if(fTracks) delete fTracks;
    fTracks = NULL;
}

// #############################################################################
Bool_t AliHLTTPCDisplay::LoadGeometrie(Char_t *gfile) {
    if (gfile) {
	TFile *file = TFile::Open(gfile);
	if(!file) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::AliHLTTPCDisplay","File Open") <<"Geometry file " << gfile << " does not exist!"<<ENDLOG;
	    return kFALSE;
	}
	
	fGeom = (TGeometry*)file->Get("AliceGeom");

	file->Close();
	delete file;
    }
    return kTRUE;
}

// #############################################################################
//                 SETTER
// #############################################################################
void AliHLTTPCDisplay::SetHistPadRowAxis() {
    // Set Axis range of Histogramm, due to variable NPads per padrow

    fNPads = AliHLTTPCTransform::GetNPads(fPadRow);
    fHistrawcl->SetAxisRange(0,fNPads);
    fHistraw->SetAxisRange(0,fNPads);
    fHistrawcl->SetAxisRange(0,fNTimes,"Y");
    fHistraw->SetAxisRange(0,fNTimes,"Y");
}

void AliHLTTPCDisplay::SetSliceArray() {
    Int_t slice=0;
    Int_t minSlice = fMinSlice; 
    Int_t maxSlice = fMaxSlice; 
    Int_t realslice = 0;

    for (slice=0;slice<=35;slice++){
	fSliceArray[slice] = kFALSE;
    }

    // Single Slice, or Range
    if (minSlice > maxSlice) maxSlice += 17;
	
    for (slice=minSlice;slice<=maxSlice;slice++){
	realslice = slice % 18;
	fSliceArray[realslice] = kTRUE;
	fSliceArray[realslice+18] = kTRUE;
    }

    // Pair of Slices
    if (fSlicePair) {
	minSlice = fMinSlice + 9;
	maxSlice = fMaxSlice + 9;
	
	if (minSlice > maxSlice) maxSlice += 17;
	
	for (slice=minSlice;slice<=maxSlice;slice++){
	    realslice = slice % 18;
	    fSliceArray[realslice] = kTRUE;	
	    fSliceArray[realslice+18] = kTRUE;
	}
    }
}

// #############################################################################
//                 SETUP
// #############################################################################
void AliHLTTPCDisplay::SetupCluster(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data)  {  

    if (data && slice>=0 && slice<36 && patch>=0 && patch<AliHLTTPCTransform::GetNPatches()) {
	if (fClusters[slice][patch]!=NULL) {
	    delete(fClusters[slice][patch]);
	    fClusters[slice][patch]=NULL;
	}
	Int_t arraysize=nofClusters*sizeof(AliHLTTPCSpacePointData);
	fClusters[slice][patch] = (AliHLTTPCSpacePointData*)new Byte_t[arraysize];
	if (fClusters[slice][patch]) {
	    memcpy(fClusters[slice][patch], data, arraysize);
	    fNcl[slice][patch]=nofClusters;
	} else {
	    fNcl[slice][patch]=nofClusters;
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupCluster","memory allocation") << "memory allocation failed "<<ENDLOG; 
	}
    } else LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupCluster","argument check") << "invalid argument "<<ENDLOG; 
}

// #############################################################################
void AliHLTTPCDisplay::SetupTracks(AliHLTTPCTrackArray *tracks) {
    fTracks=tracks;

    // Set USED cluster
    Int_t ntracks = fTracks->GetNTracks();

    for(Int_t j=0; j<ntracks; j++) { 	
	AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
	if(!gtrack) continue;

	Int_t nHits = gtrack->GetNHits();
	UInt_t *hitnum = gtrack->GetHitNumbers();

	for(Int_t h=0; h<nHits; h++){
	  
	    UInt_t id=hitnum[h];
	    Int_t slice = (id>>25) & 0x7f;
	    Int_t patch = (id>>22) & 0x7;
	    UInt_t pos = id&0x3fffff;	    
		
	    AliHLTTPCSpacePointData *points = fClusters[slice][patch];
	    
	    if(!points) {
		LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::Draw3D","Clusterarray") <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
		continue;
	    }
 
	    if(pos>=fNcl[slice][patch]) {
		LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::Draw3D","Clusterarray") <<"Pos is too large: pos "<<pos <<" ncl "<<fNcl[slice][patch]<<ENDLOG;
		continue;
	    }
	    points[pos].fUsed = kTRUE;
	}
    }
}

// #############################################################################
void AliHLTTPCDisplay::SetupHist(){

    Int_t maxpads = 150;
    fNTimes = AliHLTTPCTransform::GetNTimeBins();
    Float_t xyz[3];
    AliHLTTPCTransform::RawHLT2Global(xyz, 0, 0, 0, 0);
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::ccc","") << "t="<< fNTimes  <<"maxZ="<<xyz[2] << "|" << xyz[0]<< "|"<< xyz[1]<< ENDLOG;

    if ( fHistraw ){
	delete fHistraw;
	fHistraw = NULL;
    }
    if ( fHistrawcl ){
	delete fHistrawcl;
	fHistrawcl = NULL;
    }

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

    if ( fHistallresiduals){
	delete fHistallresiduals;
	fHistallresiduals = NULL;
    }

    if ( fHistcharge){
	delete fHistcharge;
	fHistcharge = NULL;
    }

    // Setup the histograms
    Int_t padbinning = maxpads*10;
    fHistraw = new TH2F("fHistraw","Selected PadRow with found Clusters;Pad #;Timebin #",maxpads,0,maxpads-1,fNTimes,0,fNTimes-1);
    fHistrawcl = new TH1F("fHistrawcl","",padbinning,0,maxpads-1);
    fHistpad1 = new TH1F ("fHistpad1","Selected Pad -1;Timebin #",fNTimes,0,fNTimes-1);
    fHistpad2 = new TH1F ("fHistpad2","Selected Pad;Timebin #",fNTimes,0,fNTimes-1); 
    fHistpad3 = new TH1F ("fHistpad3","Selected Pad +1;Timebin #",fNTimes,0,fNTimes-1);
    fHistallresiduals = new TH1F ("fHistallresiduals","Residuals of all Tracks in selected slices;residuals",5000,0,100);
    fHistcharge = new TH1F ("fHistcharge","Cluster distribution per charge;charge;#cluster",500,0,6000);

    fHistraw->SetOption("COLZ"); 

  

    fHistallresiduals->SetTitleSize(0.03);
    fHistallresiduals->GetXaxis()->SetLabelSize(0.03);
    fHistallresiduals->GetXaxis()->SetTitleSize(0.03);
    fHistallresiduals->GetYaxis()->SetLabelSize(0.03);
    fHistallresiduals->GetYaxis()->SetTitleSize(0.03);

    fHistcharge->SetTitleSize(0.03);
    fHistcharge->GetXaxis()->SetLabelSize(0.03);
    fHistcharge->GetXaxis()->SetTitleSize(0.03);
    fHistcharge->GetYaxis()->SetLabelSize(0.03);
    fHistcharge->GetYaxis()->SetTitleSize(0.03);

    fHistraw->SetTitleSize(0.03);
    fHistraw->GetXaxis()->SetLabelSize(0.03);
    fHistraw->GetXaxis()->SetTitleSize(0.03);
    fHistraw->GetYaxis()->SetLabelSize(0.03);
    fHistraw->GetYaxis()->SetTitleSize(0.03);

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

    gStyle->SetPalette(1);
    
    SetHistPadRowAxis();

}

// ####################################################################################################
void AliHLTTPCDisplay::FillPadRow(Int_t patch, ULong_t dataBlock, ULong_t dataLen){
    AliHLTTPCDigitReaderPacked* fDigitReader = new AliHLTTPCDigitReaderPacked();
    bool readValue = true;
    Int_t rowOffset = 0;

    // Initialize RAW DATA
    Int_t firstRow = AliHLTTPCTransform::GetFirstRow(patch);
    Int_t lastRow = AliHLTTPCTransform::GetLastRow(patch);

    // Outer sector, patches 2, 3, 4, 5 -  start counting in patch 2 with row 0
    if ( patch >= 2 ) rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

    // Initialize block for reading packed data
    void* tmpdataBlock = (void*) dataBlock;
    fDigitReader->InitBlock(tmpdataBlock,dataLen,firstRow,lastRow);

    readValue = fDigitReader->Next();

    if (!readValue){	
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::FillPadRow","Read first value") << "No value in data block" << ENDLOG;
	return;
    }
    
    // FILL PADROW 3D --- Initialize the colorbins
    if (fSwitch3DPadRow){
	for (UInt_t ii=0;ii < 20;ii++){
	    fbinct[ii] = 0;
	    fcolorbin[ii] = 0;
	}

	// read number of entries in colorbin
	while ( readValue ){ 

	    Int_t row = fDigitReader->GetRow() + rowOffset;
	    
	    if (row == fPadRow){    
		UInt_t charge = fDigitReader->GetSignal();
		
		for (UInt_t ii=0;ii < 19;ii++){
		    if ( charge > (ii*21) && charge <= ((ii*21) + 21) )	fcolorbin[ii]++;
		}
                // larger than 19 * 21	
		if (charge > 399 ) fcolorbin[19]++;
	    }

	    // read next value
	    readValue = fDigitReader->Next();
      
	    if(!readValue) break; //No more value
	} 
	//Initialize fpmarr[color][3*colorbin[ii]]  
	fpmarr[0] = new Float_t[fcolorbin[0]*3]; 
	fpmarr[1] = new Float_t[fcolorbin[1]*3]; 
	fpmarr[2] = new Float_t[fcolorbin[2]*3]; 
	fpmarr[3] = new Float_t[fcolorbin[3]*3]; 
	fpmarr[4] = new Float_t[fcolorbin[4]*3];  
	fpmarr[5] = new Float_t[fcolorbin[5]*3]; 
	fpmarr[6] = new Float_t[fcolorbin[6]*3]; 
	fpmarr[7] = new Float_t[fcolorbin[7]*3]; 
	fpmarr[8] = new Float_t[fcolorbin[8]*3]; 
	fpmarr[9] = new Float_t[fcolorbin[9]*3]; 
	fpmarr[10] = new Float_t[fcolorbin[10]*3]; 
	fpmarr[11] = new Float_t[fcolorbin[11]*3]; 
	fpmarr[12] = new Float_t[fcolorbin[12]*3]; 
	fpmarr[13] = new Float_t[fcolorbin[13]*3]; 
	fpmarr[14] = new Float_t[fcolorbin[14]*3]; 
	fpmarr[15] = new Float_t[fcolorbin[15]*3]; 
	fpmarr[16] = new Float_t[fcolorbin[16]*3]; 
	fpmarr[17] = new Float_t[fcolorbin[17]*3]; 
	fpmarr[18] = new Float_t[fcolorbin[18]*3]; 
	fpmarr[19] = new Float_t[fcolorbin[19]*3]; 
	
	// Rewind the raw reader and fill the polymarker3D
	fDigitReader->InitBlock(tmpdataBlock,dataLen,firstRow,lastRow);
	
	readValue = fDigitReader->Next();
    } // END if (fSwitch3DPadRow)

    // -- Fill Raw Data
    while ( readValue ){ 

	Int_t row = fDigitReader->GetRow() + rowOffset;

	// select padrow to fill in histogramm
	if (row == fPadRow){    
	    UChar_t pad = fDigitReader->GetPad();
	    UShort_t time = fDigitReader->GetTime();
	    UInt_t charge = fDigitReader->GetSignal();
	    Float_t xyz[3];
	    fHistraw->Fill(pad,time,charge);

	    if (pad == (fPad-1) ) fHistpad1->Fill(time,charge);
	    if (pad == fPad) fHistpad2->Fill(time,charge);
	    if (pad == (fPad+1) ) fHistpad3->Fill(time,charge);

	    if (fSwitch3DPadRow) {
		// Transform raw coordinates to local coordinates
		AliHLTTPCTransform::RawHLT2Global(xyz, fSlicePadRow, fPadRow, pad, time);

		for (UInt_t ii=0;ii < 19;ii++){
		    if ( charge > (ii*21) && charge <= ((ii*21) + 21) ){
			fpmarr[ii][fbinct[ii]] = xyz[0];
			fpmarr[ii][fbinct[ii]+1] = xyz[1];
			fpmarr[ii][fbinct[ii]+2] = xyz[2];
			fbinct[ii] += 3;
		    }
		}
		// larger than 19 * 21
		if (charge > 399 ) {
		    fpmarr[19][fbinct[19]] = xyz[0];
		    fpmarr[19][fbinct[19]+1] = xyz[1];
		    fpmarr[19][fbinct[19]+2] = xyz[2];
		    fbinct[19] += 3;
		}
	    } // END if (fSwitch3DPadRow)
	
	}
	
	// read next value
	readValue = fDigitReader->Next();
      
	//Check where to stop:
	if(!readValue) break; //No more value
    } 
    
    if ( fDigitReader )
	delete fDigitReader;
    fDigitReader = NULL;

    AliHLTTPCSpacePointData *points = fClusters[fSlicePadRow][patch];
    if(!points) return;
    Int_t npoints = fNcl[fSlicePadRow][patch];
    
    Float_t xyz[3];
    for(Int_t i=0; i<npoints; i++){
	xyz[0] = points[i].fX;
	xyz[1] = points[i].fY;
	xyz[2] = points[i].fZ;
	
	Int_t clrow = AliHLTTPCTransform::GetPadRow(xyz[0]);
	// select padrow to fill in histogramm
	if (clrow == fPadRow){
	    AliHLTTPCTransform::LocHLT2Raw(xyz, fSlicePadRow, fPadRow);
	    fHistrawcl->Fill(xyz[1],xyz[2]);
	}
    }
}


// #############################################################################
void AliHLTTPCDisplay::ResetHistPadRow(){  
    fHistraw->Reset();   
    fHistrawcl->Reset(); 
    fHistpad1->Reset(); 
    fHistpad2->Reset();  
    fHistpad3->Reset(); 
}

// #############################################################################
void AliHLTTPCDisplay::ResetHistResiduals(){  
    fHistallresiduals->Reset();
}

// #############################################################################
void AliHLTTPCDisplay::ResetHistCharge(){  
    fHistcharge->Reset();
}


// #############################################################################
//                 DRAWER
// #############################################################################
void AliHLTTPCDisplay::DrawGeomSector(Int_t sector) {  
  Char_t fname[256];
  Int_t realsector = sector;// % 18;
  
  if (realsector < 10){
    sprintf(fname,"LS0%d",realsector);
    fGeom->GetNode(fname)->SetLineColor(fLineColor);
    fGeom->GetNode(fname)->Draw("same");
    sprintf(fname,"US0%d",realsector);
    fGeom->GetNode(fname)->SetLineColor(fLineColor); 
    fGeom->GetNode(fname)->Draw("same");
  }
  else {
    sprintf(fname,"LS%d",realsector);
    fGeom->GetNode(fname)->SetLineColor(fLineColor);
    fGeom->GetNode(fname)->Draw("same");
    sprintf(fname,"US%d",realsector);
    fGeom->GetNode(fname)->SetLineColor(fLineColor); 
    fGeom->GetNode(fname)->Draw("same");
  }   
}
// #############################################################################
void AliHLTTPCDisplay::DrawHistPadRow(){  
    Char_t title[256];
    sprintf(title,"Selected PadRow %d with found Clusters",fPadRow);

    fHistraw->SetTitle(title);
    fHistraw->SetStats(kFALSE);
    fHistraw->Draw("COLZ");

    fHistrawcl->SetStats(kFALSE);
    fHistrawcl->SetMarkerStyle(28);
    fHistrawcl->SetMarkerSize(2);
    fHistrawcl->SetMarkerColor(1);
    fHistrawcl->Draw("psame");
}

// #############################################################################
void AliHLTTPCDisplay::DrawHistPad1(){  
    Char_t title[256];
    sprintf(title,"Selected Pad %d",fPad -1);
    fHistpad1->SetStats(kFALSE);
    fHistpad1->SetTitle(title);
    fHistpad1->Draw();
}

// #############################################################################
void AliHLTTPCDisplay::DrawHistPad2(){  
    Char_t title[256];
    sprintf(title,"Selected Pad %d",fPad);

    fHistpad2->SetStats(kFALSE);
    fHistpad2->SetTitle(title);
    fHistpad2->Draw();
}

// #############################################################################
void AliHLTTPCDisplay::DrawHistPad3(){  
    Char_t title[256];
    sprintf(title,"Selected Pad %d",fPad +1);

    fHistpad3->SetStats(kFALSE);
    fHistpad3->SetTitle(title);
    fHistpad3->Draw();
}

// #############################################################################
void AliHLTTPCDisplay::DrawHistResiduals(){  
    if (fSwitch3DTracks){
    
        // Residual histogram for 1 track
	if (fSelectTrackSwitch){
	    Char_t title[256];
	    sprintf(title,"Residuals of Track %d in Slice %d",fSelectTrack, fSelectTrackSlice );
	 
//	    fHistresiduals->SetMarkerStyle(2);
//	    fHistresiduals->SetStats(kFALSE);
//	    fHistresiduals->SetTitle(title);
//	    fHistresiduals->Draw("p");
	    fGraphresiduals->SetTitle(title);	
	    fGraphresiduals->GetXaxis()->SetTitle("z");	
	    fGraphresiduals->GetYaxis()->SetTitle("residuals");
	    fGraphresiduals->Draw("A*");
	    
	}
	// Global residuals histogram
	else{
	    fHistallresiduals->SetStats(kFALSE);
	    fHistallresiduals->Draw();
	}
    }
}

// #############################################################################
void AliHLTTPCDisplay::DrawHistCharge(){  
    if (fSwitch3DCluster){
//	fHistcharge->SetStats(kFALSE);
	fHistcharge->Draw();
    }
}

// #############################################################################
void AliHLTTPCDisplay::Draw3D(){    	
    
    TView *v = new TView(1);
    v->SetRange(-800,-800,-800,800,800,800);

    Float_t* etaRange = NULL;   // ------  STILL TO FIX
    

    //--------------------------------------------------------------------------------------------
    // DRAW 3D CLUSTER
    //--------------------------------------------------------------------------------------------
    if (fSwitch3DCluster){
	for (Int_t slice=0; slice <= 35; slice++){

	    UInt_t maxCharge = 0;
	    Float_t maxXYZ[3];
	    UChar_t padrow;

	    if (!fSliceArray[slice]) continue;
	    
	    for(Int_t p=0;p<6;p++){

		AliHLTTPCSpacePointData *points = fClusters[slice][p];
		if(!points) continue;
		Int_t npoints = fNcl[slice][p];
		TPolyMarker3D *pm = new TPolyMarker3D(npoints);
	
		Float_t xyz[3];
		for(Int_t i=0; i<npoints; i++){
		    // Used  cluster only
		    if (fSelectCluster == 1  && points[i].fUsed == kFALSE) continue; 
		    // Unused cluster only
		    if (fSelectCluster == 2  && points[i].fUsed == kTRUE) continue; 
		    
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
		    
		    pm->SetPoint(i,xyz[0],xyz[1],xyz[2]);
		    
		    // Fill Charge Histogram
		    fHistcharge->Fill(points[i].fCharge);
		    if (points[i].fCharge > maxCharge ){
			maxCharge = points[i].fCharge; 
			maxXYZ[0] = points[i].fX;
			maxXYZ[1] = points[i].fY;
			maxXYZ[2] = points[i].fZ;
			padrow = points[i].fPadRow;
		    }
		    
		}
		pm->SetMarkerSize(4);
		pm->SetMarkerColor(2); 
		pm->Draw("");
	    } 	    
	
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::CHARGE","")<< "MAX CHARGE =" <<  maxCharge <<  "  slice=" << slice 
								   << " z=" << maxXYZ[0]<< " y=" << maxXYZ[1]<< " z=" << maxXYZ[2]

								   << ENDLOG;


	}
    }   // END - DRAW 3D CLUSTER 

    //--------------------------------------------------------------------------------------------
    // DRAW 3D TRACKS
    //--------------------------------------------------------------------------------------------
    if (fSwitch3DTracks){

	Int_t ntracks = fTracks->GetNTracks();

	TPolyLine3D *line = new TPolyLine3D[ntracks];
	TPolyLine3D *lineT = new TPolyLine3D[ntracks];
	Float_t xCl[176];
	Float_t yCl[176];
	Float_t zCl[176];

	Float_t xT[176];
	Float_t yT[176];
	Float_t zT[176];

	Float_t res[176];

	Int_t trackcounter = 0;

	for(Int_t j=0; j<ntracks; j++) { 	

	    AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
	    if(!gtrack) continue;

	    Int_t nHits = gtrack->GetNHits();
	    UInt_t *hitnum = gtrack->GetHitNumbers();
	    Int_t hitcount=0;	

	    Bool_t nexttrack = kFALSE;

	    TPolyMarker3D *pm = new TPolyMarker3D(nHits,7);
	    TPolyMarker3D *pmT = new TPolyMarker3D(nHits,7);

	    TPolyMarker3D *pmL = new TPolyMarker3D(1,2);
	    TPolyMarker3D *pmF = new TPolyMarker3D(1,2);

	    Double_t lambda = 0.;  // dipAngle lambda
	    Double_t r = 0.;       // radius
	    Double_t kappa = 0.;   // curvature = 1/R , signed
	    Double_t xyz0[3];      // startingpoint of track
	    Double_t xyzT[3];      // point on track
	    Double_t xyzC[3];      // cluster
	    Double_t s = 0.;       // length of track
	    Double_t phi0 = 0.;    // azimuthal angle of startingpoint, with respect to helix axis
	    Double_t bfield = 0;   // BField

	    Double_t xyzL[3];      // lastpoint of track
	    Double_t xyzF[3];      // firstpoint of track

	    Double_t maxZ = 0;     // range of the histogram
	    Double_t minZ =99999.; // range of the histogram

	    for(Int_t h=0; h<nHits; h++){

		UInt_t id=hitnum[h];
		Int_t slice = (id>>25) & 0x7f;
		Int_t patch = (id>>22) & 0x7;
		UInt_t pos = id&0x3fffff; 

		// select if slice should be displayed or not
		if (!fSliceArray[slice]) {	
		    nexttrack = kTRUE;
		    break;	   
		}
		
		// select Single Track
		if (fSelectTrackSwitch){
		    if(slice != fSelectTrackSlice) {
			nexttrack = kTRUE;
			break;
		    }

		    if (trackcounter != fSelectTrack && h==0){
			trackcounter++;  
			nexttrack = kTRUE;
			break;
		    }

		    trackcounter++;
    		}

		// --> in the hit loop because of 'trackcounter++', otherwise wrong single track in slice will be selected
		if((fPtThreshold > 0) && (gtrack->GetPt()< fPtThreshold)) {	
		    nexttrack = kTRUE;
		    break;
		}
		
		// --> in the hit loop because of 'trackcounter++', otherwise wrong single track in slice will be selected
		if(nHits < fMinHits) {	
		    nexttrack = kTRUE;
		    break;
		}

		AliHLTTPCSpacePointData *points = fClusters[slice][patch];
		
		if(!points) {
		    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::Draw3D","Clusterarray") <<"No points at slice "<<slice<<" patch "<<patch
											<<" pos "<<pos<<ENDLOG;
		    continue;
		}
 
		if(pos>=fNcl[slice][patch]) {
		    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::Draw3D","Clusterarray") <<"Pos is too large: pos "<<pos <<" ncl "
											<<fNcl[slice][patch]<<ENDLOG;
		    continue;
		}

		// set the data for the residuals
		if (h == 0){
		    // -> Curvature / Radius / Phi0
		    //kappa = gtrack->GetKappa(); 
		    //r = gtrack->GetRadius();
		    //phi0 = gtrack->GetPhi0();
		    //bfield = AliHLTTPCTransform::GetBFieldValue();
		    
		    lambda = atan( gtrack->GetTgl() ); 
		    
		    bfield = 0.0029980 * 0.4; // KORRIGIERE B für =0.4 T Feld
		    
		    xyz0[0] = gtrack->GetFirstPointX();
		    xyz0[1] = gtrack->GetFirstPointY();
		    xyz0[2] = gtrack->GetFirstPointZ();
		    
		    xyzF[0] = gtrack->GetFirstPointX();
		    xyzF[1] = gtrack->GetFirstPointY();
		    xyzF[2] = gtrack->GetFirstPointZ();

		    xyzL[0] = gtrack->GetLastPointX();
		    xyzL[1] = gtrack->GetLastPointY();
		    xyzL[2] = gtrack->GetLastPointZ();

		    pmL->SetPoint(0,xyzL[0],xyzL[1],xyzL[2]);
		    pmF->SetPoint(0,xyzF[0],xyzF[1],xyzF[2]);

		    phi0 = gtrack->GetPsi() + (gtrack->GetCharge() * AliHLTTPCTransform::PiHalf() ); 
		    
		    if (bfield != 0.){
			r = gtrack->GetPt() / bfield; 
			kappa = - gtrack->GetCharge() * 1. / r;
		    }
		    else {
			r = 999999; // just infinity
			kappa = 0;
		    }

		    //Write Track Parameters for single track
		    if (fSelectTrackSwitch){
			fTrackParam.id = trackcounter - 1;
			fTrackParam.nHits = nHits;
			fTrackParam.charge = gtrack->GetCharge();
			fTrackParam.lambda = lambda;
			fTrackParam.kappa = kappa;
			fTrackParam.radius = r;
			fTrackParam.slice = slice;
			fTrackParam.phi0 = phi0;
			fTrackParam.pt = gtrack->GetPt();
			fTrackParam.bfield = bfield;
			fTrackParam.xyzF[0] = gtrack->GetFirstPointX();
			fTrackParam.xyzF[1] = gtrack->GetFirstPointY();
			fTrackParam.xyzF[2] = gtrack->GetFirstPointZ();
			fTrackParam.xyzL[0] = gtrack->GetLastPointX();
			fTrackParam.xyzL[1] = gtrack->GetLastPointY();
			fTrackParam.xyzL[2] = gtrack->GetLastPointZ();
			fTrackParam.psi = gtrack->GetPsi();
		    }
		}


		Float_t xyzCtmp[3];    // cluster tmp

		xyzCtmp[0] = points[pos].fX;
		xyzCtmp[1] = points[pos].fY;
		xyzCtmp[2] = points[pos].fZ;

		AliHLTTPCTransform::Local2Global(xyzCtmp,slice);
		    
		xCl[h] = xyzCtmp[0];
		yCl[h] = xyzCtmp[1];
		zCl[h] = xyzCtmp[2];

		// FILL POLYMARKER FOR THE ORIGINAL TRACKS
		pm->SetPoint(h,xCl[h],yCl[h],zCl[h]);  

		xyzC[0] = (Double_t) xyzCtmp[0];
		xyzC[1] = (Double_t) xyzCtmp[1];
		xyzC[2] = (Double_t) xyzCtmp[2];

		xyzT[2] = xyzC[2];

		//calculate length
		s = ( xyzT[2] - xyz0[2] ) / sin(lambda);

		// calculate the corresponding coordinates on the track
		xyzT[0] = xyz0[0] +  r * ( cos( phi0 + (s*kappa*cos(lambda)) ) - cos( phi0 ) );
		xyzT[1] = xyz0[1] +  r * ( sin( phi0 + (s*kappa*cos(lambda)) ) - sin( phi0 ) );

		Double_t deltaX = ( xyzC[0] - xyzT[0] );
		Double_t deltaY = ( xyzC[1] - xyzT[1] );
		
		Double_t residual = sqrt( deltaX*deltaX + deltaY*deltaY );
		res[h] = (Float_t) residual;

		if (maxZ < fabs(xyzT[2])) maxZ = fabs(xyzT[2]);
		if (minZ > fabs(xyzT[2])) minZ = fabs(xyzT[2]);

		// FILL RESIDUALS HISTOGRAM
		fHistallresiduals->Fill(residual);
		
		// FILL ARRAYS IN ORDER TO DRAW THE TRACKPOINTS, OUT OF THE PARAMETER
		xT[h] = xyzT[0];
		yT[h] = xyzT[1];
		zT[h] = xyzT[2];

		// FILL POLYMARKER FOR THE NEW TRACKS
		pmT->SetPoint(h,xT[h],yT[h],zT[h]);
       	  
		hitcount++;
	    }

	    if(nexttrack) continue;
	    if(hitcount==0) continue;
	

	    if ( fGraphresiduals){
		delete fGraphresiduals;
		fGraphresiduals = NULL;
	    }

	    // FILL RESIDUALS GRAPH
	    fGraphresiduals = new TGraph(nHits,zT,res);
	    fGraphresiduals->GetXaxis()->SetLabelSize(0.02);
	    fGraphresiduals->GetXaxis()->SetTitleSize(0.02);
	    fGraphresiduals->GetYaxis()->SetLabelSize(0.02);
	    fGraphresiduals->GetYaxis()->SetTitleSize(0.02);
	    
	    TPolyLine3D *currentline = &(line[j]);
	    currentline = new TPolyLine3D(nHits,xCl,yCl,zCl,"");
	    currentline->SetLineColor(4);   
	    currentline->SetLineWidth(2);
//	    currentline->Draw("same");

	    TPolyLine3D *currentlineT = &(lineT[j]);
	    currentlineT = new TPolyLine3D(nHits,xT,yT,zT,"");
	    currentlineT->SetLineColor(7);   
	    currentlineT->SetLineWidth(1);
//	    currentlineT->Draw("same");

	    
	    //Last Point of Track
	    pmL->SetMarkerSize(3);
	    pmL->SetMarkerColor(4); 
	    pmL->Draw();

	    //First Point of Track
	    pmF->SetMarkerSize(3);
	    pmF->SetMarkerColor(5); 
	    pmF->Draw();
	    
            //Original Track 
	    pm->SetMarkerSize(4);
	    pm->SetMarkerColor(6); 
	    pm->Draw();

	    //New Track
	    pmT->SetMarkerSize(4);
            pmT->SetMarkerColor(3); 
	    pmT->Draw();

	} // END for tracks

    }   // END - DRAW 3D Tracks

    //--------------------------------------------------------------------------------------------
    // DRAW 3D GEOMETRY
    //--------------------------------------------------------------------------------------------
    if (fSwitch3DGeometry){

	for (Int_t slice=0; slice <= 17; slice++){
	    if (!fSliceArray[slice]) continue;
	    DrawGeomSector(slice);
	}
    }   // END - DRAW 3D GEOMETRY
    
    //--------------------------------------------------------------------------------------------
    // DRAW 3D PadRow
    //--------------------------------------------------------------------------------------------
    if (fSwitch3DPadRow && fSliceArray[fSlicePadRow]){
	Int_t markercolor = 51;

	for (UInt_t ii=0;ii < 20;ii++){
	    if (fcolorbin[ii]> 0){
		
		TPolyMarker3D *pm = new TPolyMarker3D(fcolorbin[ii], fpmarr[ii], 7 );

		pm->SetMarkerColor(markercolor); 
		pm->Draw(""); 
	    }

	    // in order to have the SetPalette(1), so called "pretty"
	    if (ii % 2 == 0 ) markercolor += 2;
	    else  markercolor += 3;
	}
    }

    //--------------------------------------------------------------------------------------------
    // DRAW 3D 
    //--------------------------------------------------------------------------------------------
    v->ZoomView(0,4);
    v->Draw();   
}



 
