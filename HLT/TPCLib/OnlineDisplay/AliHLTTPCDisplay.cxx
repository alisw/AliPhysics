// @(#) $Id$
// Original: AliHLTDisplay.cxx,v 1.26 2005/06/14 10:55:21 cvetan 

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

#define TRACKHELIX 0
#define TRACKPOLYMARKER 0
#define BACKWARD 0
#define FIRSTLASTPOINT 0

#define TRACKCOLOR 
#define USEDCLUSTERCOLOR
#define UNUSEDCLUSTERCOLOR

#include "AliHLTStdIncludes.h"
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
#include <TMultiGraph.h>
#include <TAttText.h>
#include <TAxis.h>

#if TRACKHELIX
#include <THelix.h>
#endif

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
void AliHLTTPCDisplay::InitDisplay(Char_t *gfile) {
    //constructor
    memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
    memset(fNcl, 0, 36*6*sizeof(UInt_t)); 

    fTracks = NULL;
    fHistrawcl = NULL;
    fHistraw = NULL;
    fHistpad1 = NULL;
    fHistpad2 = NULL;
    fHistpad3 = NULL;
    fHistallresidualsY = NULL;   
    fHistallresidualsZ = NULL;
    fHistcharge = NULL;
    fGraphresidualsY = NULL;
    fGraphresidualsZ = NULL;
    fGraphresidualsYLength = NULL;
    fGraphresidualsZLength = NULL;


    fGeom = NULL;
// ---------------------------------------------------
// In order to be backward compatible
// ---------------------------------------------------
#if BACKWARD
    //fc1 = NULL;
#endif 
// ---------------------------------------------------
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
    fKeepView = kFALSE;

    fSwitch3DCluster = kFALSE;
    fSwitch3DTracks = kFALSE;
    fSwitch3DPadRow = kFALSE;
    fSwitch3DGeometry = kFALSE;

    AliHLTTPCTransform::SetBField(0.4);
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
//                 EXECUTER
// #############################################################################
void AliHLTTPCDisplay::ExecPadRow(){
   int event = gPad->GetEvent();
   if (event != 11) return;

   printf("TEST !!!!!!!!!!!!!!!");
/*   int px = gPad->GetEventX();
   TObject *select = gPad->GetSelected();
   if (!select) return;
   if (select->InheritsFrom("TH1")) {
      TH1 *h = (TH1*)select;
      Float_t xx = gPad->AbsPixeltoX(px);
      Float_t x  = gPad->PadtoX(xx);
      Int_t binx = h->GetXaxis()->FindBin(x);
      printf("event=%d, hist:%s, bin=%d, content=%f\n",event,h->GetName(),binx,h->GetBinContent(binx));
   }

*/

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
	    points[pos].fTrackN = j;
	}
    }
}

// #############################################################################
void AliHLTTPCDisplay::SetupHist(){

    Int_t maxpads = 150;
    fNTimes = AliHLTTPCTransform::GetNTimeBins();

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

    if ( fHistallresidualsY){
	delete fHistallresidualsY;
	fHistallresidualsY = NULL;
    }

    if ( fHistallresidualsZ){
	delete fHistallresidualsZ;
	fHistallresidualsZ = NULL;
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
    fHistallresidualsY = new TH1F ("fHistallresiduals","Y Residuals of all Tracks in selected slices;residuals",5000,-100,100);
    fHistallresidualsZ = new TH1F ("fHistallresiduals","Z Residuals of all Tracks in selected slices;residuals",5000,-100,100);
    fHistcharge = new TH1F ("fHistcharge","Cluster distribution per charge;charge;#cluster",5000,0,30000);

    fHistraw->SetOption("COLZ"); 

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
#if defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
    AliHLTTPCDigitReader* digitReader = new AliHLTTPCDigitReaderPacked();
    bool readValue = true;
    Int_t rowOffset = 0;

    // Initialize RAW DATA
    Int_t firstRow = AliHLTTPCTransform::GetFirstRow(patch);
    Int_t lastRow = AliHLTTPCTransform::GetLastRow(patch);

    // Outer sector, patches 2, 3, 4, 5 -  start counting in patch 2 with row 0
    if ( patch >= 2 ) rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

    // Initialize block for reading packed data
    void* tmpdataBlock = (void*) dataBlock;
    digitReader->InitBlock(tmpdataBlock,dataLen,firstRow,lastRow,patch,0);

    readValue = digitReader->Next();

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

	    Int_t row = digitReader->GetRow() + rowOffset;
	    
	    if (row == fPadRow){    
		UInt_t charge = digitReader->GetSignal();
		
		for (UInt_t ii=0;ii < 19;ii++){
		    if ( charge > (ii*15) && charge <= ((ii*15) + 15) )	fcolorbin[ii]++;
		}
                // larger than 19 * 15	
		if (charge > 285 ) fcolorbin[19]++;
	    }

	    // read next value
	    readValue = digitReader->Next();
      
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
	digitReader->InitBlock(tmpdataBlock,dataLen,firstRow,lastRow,patch,0);
	
	readValue = digitReader->Next();
    } // END if (fSwitch3DPadRow)

    // -- Fill Raw Data
    while ( readValue ){ 

	Int_t row = digitReader->GetRow() + rowOffset;

	// select padrow to fill in histogramm
	if (row == fPadRow){    
	    UChar_t pad = digitReader->GetPad();
	    UShort_t time = digitReader->GetTime();
	    UInt_t charge = digitReader->GetSignal();
	    Float_t xyz[3];
	    fHistraw->Fill(pad,time,charge);

	    if (pad == (fPad-1) ) fHistpad1->Fill(time,charge);
	    if (pad == fPad) fHistpad2->Fill(time,charge);
	    if (pad == (fPad+1) ) fHistpad3->Fill(time,charge);

	    if (fSwitch3DPadRow) {
		// Transform raw coordinates to local coordinates
		AliHLTTPCTransform::RawHLT2Global(xyz, fSlicePadRow, fPadRow, pad, time);

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
	
	}
	
	// read next value
	readValue = digitReader->Next();
      
	//Check where to stop:
	if(!readValue) break; //No more value
    } 
    
    if ( digitReader )
	delete digitReader;
    digitReader = NULL;

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
#else //! if defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
    HLTFatal("DigitReaderPacked not available - check your build");
#endif //defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
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
    fHistallresidualsY->Reset();
    fHistallresidualsZ->Reset();
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
void AliHLTTPCDisplay::DrawHistResiduals(Bool_t ySwitch){  
    if (fSwitch3DTracks){
	if (ySwitch){
	    // Y Residual histogram for 1 track

	    if (fSelectTrackSwitch){
		Char_t title[256];
		sprintf(title,"Y Residuals of Track %d in Slice %d",fSelectTrack, fSelectTrackSlice );

		TMultiGraph *mgY = new TMultiGraph();


//		fGraphresidualsY->SetTitle(title);	
		fGraphresidualsY->GetXaxis()->SetTitle("padrow");	
		fGraphresidualsY->GetYaxis()->SetTitle("residuals");
//		fGraphresidualsY->Draw("A*");
		fGraphresidualsY->GetXaxis()->SetLabelSize(0.02);
		fGraphresidualsY->GetXaxis()->SetTitleSize(0.02);
		fGraphresidualsY->GetYaxis()->SetLabelSize(0.02);
		fGraphresidualsY->GetYaxis()->SetTitleSize(0.02);
		fGraphresidualsYLength->SetMarkerColor(2);
		fGraphresidualsYLength->SetMarkerStyle(5);
		fGraphresidualsY->SetMarkerColor(1);
		fGraphresidualsY->SetMarkerStyle(3);

//		fGraphresidualsY->Draw("A*");
//		fGraphresidualsYLength->Draw("*");

		mgY->Add(fGraphresidualsY);
		mgY->Add(fGraphresidualsYLength);
		mgY->SetTitle(title);
//		mgY->GetXaxis()->SetTitle("padrow");	
//		mgY->GetYaxis()->SetTitle("residuals");
		mgY->Draw("AP");
	    }
	    // Global residuals histogram
	    else{
		fHistallresidualsY->SetStats(kFALSE);
		fHistallresidualsY->Draw();
	    }
	}
	else {
            // Z Residual histogram for 1 track
	    if (fSelectTrackSwitch){
		Char_t title[256];
		sprintf(title,"Z Residuals of Track %d in Slice %d",fSelectTrack, fSelectTrackSlice );

		TMultiGraph *mgZ = new TMultiGraph();

		fGraphresidualsZ->SetTitle(title);	
		fGraphresidualsZ->GetXaxis()->SetTitle("padrow");	
		fGraphresidualsZ->GetYaxis()->SetTitle("residuals");
		fGraphresidualsZ->GetXaxis()->SetLabelSize(0.02);
		fGraphresidualsZ->GetXaxis()->SetTitleSize(0.02);
		fGraphresidualsZ->GetYaxis()->SetLabelSize(0.02);
		fGraphresidualsZ->GetYaxis()->SetTitleSize(0.02);
//		fGraphresidualsZLength->Draw("F*");
//		fGraphresidualsZ->Draw("A*");
	
		mgZ->Add(fGraphresidualsZ);
//		mgZ->Add(fGraphresidualsZLength);
		mgZ->SetTitle(title);
		mgZ->Draw("A*");
	    }
	    // Global residuals histogram
	    else{
		fHistallresidualsZ->SetStats(kFALSE);
		fHistallresidualsZ->Draw();
	    }
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
	Int_t maxCharge = 0;

	for (Int_t slice=0; slice <= 35; slice++){

	    Int_t currenttrack = -1;

	    if (fSelectCluster == 1 && fSelectTrackSwitch && slice == fSelectTrackSlice ){

		Int_t trackcounter = 0;
		Int_t ntracks = fTracks->GetNTracks();
	
		for(Int_t j=0; j<ntracks; j++) { 	

		    AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
		    if(!gtrack) continue;

		    Int_t nHits = gtrack->GetNHits();  // Number of associated hits to track
		    Int_t tmpslice = gtrack->GetSector();

		    // --- CHECK if track is should be drawn
		    // select Single Track
		    if(tmpslice != fSelectTrackSlice) continue;
			
		    if (trackcounter != fSelectTrack){
			trackcounter++;  
			continue;
		    }
		    trackcounter++;

		    if((fPtThreshold > 0) && (gtrack->GetPt()< fPtThreshold)) continue;
		    if(nHits < fMinHits) continue;

		    currenttrack = j;
		    break;
		}
	    }
	    
	    if (!fSliceArray[slice]) continue;
	    
	    for(Int_t p=0;p<6;p++){

		AliHLTTPCSpacePointData *points = fClusters[slice][p];
		if(!points) continue;
		Int_t npoints = fNcl[slice][p];
		TPolyMarker3D *pmUsed = new TPolyMarker3D(1,6);
		TPolyMarker3D *pmUnused = new TPolyMarker3D(1,6);
		Int_t nUsedCluster = 0;
		Int_t nUnusedCluster = 0;

		Float_t xyz[3];
		for(Int_t i=0; i<npoints; i++){
		    // Used  cluster only
		    if (fSelectCluster == 1  && points[i].fUsed == kFALSE) continue; 
		    // Unused cluster only
		    if (fSelectCluster == 2  && points[i].fUsed == kTRUE) continue; 

		    // if single track is selcted draw only cluster for this track
		    if (fSelectCluster == 1 && fSelectTrackSwitch && points[i].fTrackN != currenttrack) continue;
		    
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

		    // Fill Charge Histogram
		    fHistcharge->Fill(points[i].fCharge);
		    if ((Int_t)points[i].fCharge > maxCharge ) maxCharge = (Int_t) points[i].fCharge; 
		}
		pmUsed->SetMarkerSize(1);
		pmUsed->SetMarkerColor(3); 
		pmUsed->Draw("");

		pmUnused->SetMarkerSize(1);
		pmUnused->SetMarkerColor(2); 
		pmUnused->Draw("");
	    } // END - PATCH LOOP	    
	}  // END - SLICE LOOP
	fHistcharge->SetAxisRange(0,maxCharge);
    }   // END - DRAW 3D CLUSTER 

    //--------------------------------------------------------------------------------------------
    // DRAW 3D TRACKS
    //--------------------------------------------------------------------------------------------
    if (fSwitch3DTracks){

	Int_t trackcounter = 0;
	Int_t ntracks = fTracks->GetNTracks();
	Double_t drawStep = 0.2;

	Double_t maxResidualY = 0.;
	Double_t maxResidualZ = 0.;

	TPolyLine3D *line = new TPolyLine3D[ntracks];
#if TRACKHELIX
	THelix *helix = new THelix[ntracks];
#endif
	for(Int_t j=0; j<ntracks; j++) { 	

	    AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
	    if(!gtrack) continue;

	    Int_t nHits = gtrack->GetNHits();  // Number of associated hits to track
	    Int_t slice = gtrack->GetSector();

            // --- CHECK if track is should be drawn
	    // select if slice should be displayed or not
	    if (!fSliceArray[slice]) continue;	
	    
	    // select Single Track
	    if (fSelectTrackSwitch){
		if(slice != fSelectTrackSlice) continue;
	
		if (trackcounter != fSelectTrack){
		    trackcounter++;  
		    continue;
		}
		trackcounter++;
	    }
    
	    if((fPtThreshold > 0) && (gtrack->GetPt()< fPtThreshold)) continue;
	    if(nHits < fMinHits) continue;
	    
	    TPolyMarker3D *pmL = new TPolyMarker3D(1,2);
	    TPolyMarker3D *pmF = new TPolyMarker3D(1,2);

	    Double_t radius = gtrack->GetRadius();      // radius
	    Double_t kappa = gtrack->GetKappa();        // curvature = 1/R , signed
	    Double_t lambda = atan( gtrack->GetTgl() ); // dipAngle lambda
	    Double_t phi0 = gtrack->GetPsi() + (gtrack->GetCharge() * AliHLTTPCTransform::PiHalf() ); // azimuthal angle of startingpoint, with respect to helix axis

	    Double_t xyzL[3];      // lastpoint of track
	    Double_t xyzF[3];      // firstpoint of track

	    xyzF[0] = gtrack->GetFirstPointX();
	    xyzF[1] = gtrack->GetFirstPointY();
	    xyzF[2] = gtrack->GetFirstPointZ();
	    pmF->SetPoint(0,xyzF[0],xyzF[1],xyzF[2]);

	    xyzL[0] = gtrack->GetLastPointX();
	    xyzL[1] = gtrack->GetLastPointY();
	    xyzL[2] = gtrack->GetLastPointZ();
	    pmL->SetPoint(0,xyzL[0],xyzL[1],xyzL[2]);

	    Double_t s = 0.;       // length of the track

	    // Calculate the length of the track. If it is to flat in in s,z plane use sxy, otherwise use sz
	    if (fabs(lambda) > 0.05){
                // length of track calculated out of z
		s = fabs( (xyzL[2] - xyzF[2]) / sin(lambda) ); // length of track calculated out of z
	    }
	    else {
		Double_t d =  (xyzL[0] - xyzF[0])*(xyzL[0] - xyzF[0]) + (xyzL[1] - xyzF[1])*(xyzL[1] - xyzF[1]);
		// length of track calculated out of xy
		s = fabs ( acos( 0.5 * (2 - (d / (radius*radius)))) / ( kappa * cos(lambda) ) ); 		
	    }
	    
	    Int_t nTrackPoints = 2 + (Int_t) floor(s / drawStep);

#if TRACKPOLYMARKER
	    TPolyMarker3D *pmT = new TPolyMarker3D(nTrackPoints,6);
#endif

	    Double_t *xT = new Double_t[nTrackPoints];
	    Double_t *yT = new Double_t[nTrackPoints];
	    Double_t *zT = new Double_t[nTrackPoints];

	    //Write Track Parameters for single track
	    if (fSelectTrackSwitch){
		fTrackParam.id = trackcounter - 1;
		fTrackParam.nHits = nHits;
		fTrackParam.charge = gtrack->GetCharge();
		fTrackParam.lambda = lambda;
		fTrackParam.kappa = kappa;
		fTrackParam.radius = radius;
		fTrackParam.slice = slice;
		fTrackParam.phi0 = phi0;
		fTrackParam.pt = gtrack->GetPt();
		fTrackParam.bfield = AliHLTTPCTransform::GetBFieldValue();
		fTrackParam.xyzF[0] = gtrack->GetFirstPointX();
		fTrackParam.xyzF[1] = gtrack->GetFirstPointY();
		fTrackParam.xyzF[2] = gtrack->GetFirstPointZ();
		fTrackParam.xyzL[0] = gtrack->GetLastPointX();
		fTrackParam.xyzL[1] = gtrack->GetLastPointY();
		fTrackParam.xyzL[2] = gtrack->GetLastPointZ();
		fTrackParam.psi = gtrack->GetPsi();
		fTrackParam.s = s;
	    }

	    Int_t trackPointCounter = 0;

	    for (Double_t ds = 0.; ds < s; ds = ds + drawStep){
		// FILL ARRAYS IN ORDER TO DRAW THE TRACKPOINTS, OUT OF THE PARAMETER
		xT[trackPointCounter] = xyzF[0] + radius * ( cos( phi0 + (ds*kappa*cos(lambda)) ) - cos(phi0) );
		yT[trackPointCounter] = xyzF[1] + radius * ( sin( phi0 + (ds*kappa*cos(lambda)) ) - sin(phi0) );
		zT[trackPointCounter] = xyzF[2] + ds * sin(lambda);
#if TRACKPOLYMARKER
		pmT->SetPoint(trackPointCounter,xT[trackPointCounter],yT[trackPointCounter],zT[trackPointCounter]);
#endif
		trackPointCounter++;
	    }

	    xT[trackPointCounter] = xyzF[0] + radius * ( cos( phi0 + (s*kappa*cos(lambda)) ) - cos(phi0) );
	    yT[trackPointCounter] = xyzF[1] + radius * ( sin( phi0 + (s*kappa*cos(lambda)) ) - sin(phi0) );
	    zT[trackPointCounter] = xyzF[2] + s * sin(lambda);
#if TRACKPOLYMARKER       
	    pmT->SetPoint(trackPointCounter,xT[trackPointCounter],yT[trackPointCounter],zT[trackPointCounter]);
#endif
	    // --- RESIDUALS ---
	    gtrack->Rotate(slice,kTRUE);
	    Int_t nRes = 0;  // number of resiudals
	    
	    UInt_t *hitnum = gtrack->GetHitNumbers();

	    Double_t *resY= new Double_t[nHits];
	    Double_t *resZ= new Double_t[nHits];

	    Double_t *resYLength= new Double_t[2*nHits];
	    Double_t *resZLength= new Double_t[2*nHits];

	    Double_t *padrows = new Double_t[nHits];
	    Double_t *padrowsLength = new Double_t[2*nHits];

	    for(Int_t h=0; h<nHits; h++){
		UInt_t id=hitnum[h];
		Int_t patch = (id>>22) & 0x7;
		UInt_t pos = id&0x3fffff; 

		AliHLTTPCSpacePointData *points = fClusters[slice][patch];

		Float_t xyzCtmp[3];    // cluster tmp
		Float_t xyzTtmp[3];    // track tmp

		xyzCtmp[0] = points[pos].fX;
		xyzCtmp[1] = points[pos].fY;
		xyzCtmp[2] = points[pos].fZ;

		Int_t padrow = AliHLTTPCTransform::GetPadRow(points[pos].fX);
		xyzTtmp[0] = gtrack->GetFirstPointX();
		if(gtrack->GetCrossingPoint(padrow,xyzTtmp)) {

		    Float_t deltaY = ( xyzCtmp[1] - xyzTtmp[1] );
		    Float_t deltaZ = ( xyzCtmp[2] - xyzTtmp[2] );
//		    Float_t residual = sqrt( deltaY*deltaY + deltaZ*deltaZ );
		    
		    padrows[nRes] = (Double_t) padrow;
		    resY[nRes] = (Double_t) deltaY;
		    resZ[nRes] = (Double_t) deltaZ;

		    resYLength[(2*nRes)] = 0.5 * AliHLTTPCTransform::GetPadLength(padrow);
		    resYLength[(2*nRes)+1] = -0.5 * AliHLTTPCTransform::GetPadLength(padrow);
		    resZLength[nRes] = AliHLTTPCTransform::GetZLength();
		    padrowsLength[(2*nRes)] = (Double_t) padrow;
		    padrowsLength[(2*nRes)+1] = (Double_t) padrow;

		    // FILL RESIDUALS HISTOGRAM
		    fHistallresidualsY->Fill(resY[nRes]);
		    fHistallresidualsZ->Fill(resZ[nRes]);
		    if (resY[nRes] > maxResidualY ) maxResidualY = resY[nRes];
		    if (resZ[nRes] > maxResidualZ ) maxResidualZ = resZ[nRes];
		    nRes++;
		}
	    }

	    gtrack->Rotate(slice,kFALSE);
	    // --- RESIDUALS ---

	    // Draw last point of Track
	    pmL->SetMarkerSize(3);
	    pmL->SetMarkerColor(4); 
//	    pmL->Draw();

	    // Draw first point of Track
	    pmF->SetMarkerSize(3);
	    pmF->SetMarkerColor(5); 
//	    pmF->Draw();

#if TRACKPOLYMARKER
	    // Draw Track -- as polymarker
	    pmT->SetMarkerSize(3);
	    pmT->SetMarkerColor(3); 
	    pmT->Draw();
#endif
	    // Draw Track -- as line
	    TPolyLine3D *currentline = &(line[j]);
	    currentline = new TPolyLine3D(nTrackPoints,xT,yT,zT,"");
	    currentline->SetLineColor(4);   
	    currentline->SetLineWidth(2);
	    currentline->Draw("same");

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

	    THelix *currenthelix = &(helix[j]);
	    currenthelix = new THelix(xyzF,v0,omega,hrange,kHelixZ,0);
	    currenthelix->SetLineColor(6);   
	    currenthelix->SetLineWidth(1);
	    currenthelix->Draw("same");   	    
#endif

	    //Residuals
	    if ( fGraphresidualsY){
		delete fGraphresidualsY;
		fGraphresidualsY = NULL;
	    }

	    if ( fGraphresidualsZ){
		delete fGraphresidualsZ;
		fGraphresidualsZ = NULL;
	    }
	    //Residuals
	    if ( fGraphresidualsYLength){
		delete fGraphresidualsYLength;
		fGraphresidualsYLength = NULL;
	    }

	    if ( fGraphresidualsZLength){
		delete fGraphresidualsZLength;
		fGraphresidualsZLength = NULL;
	    }



	    // FILL Y RESIDUALS GRAPH
	    fGraphresidualsY = new TGraph(nRes-1,padrows,resY);
	    fGraphresidualsYLength = new TGraph((2*nRes)-2,padrowsLength,resYLength);
	    // FILL Z RESIDUALS GRAPH
	    fGraphresidualsZ = new TGraph(nRes-1,padrows,resZ);
	    fGraphresidualsZLength = new TGraph(nRes-1,padrows,resZLength);

	    if (xT) delete xT;
	    if (yT) delete yT;
	    if (zT) delete zT;

	} // END for tracks

	fHistallresidualsY->SetAxisRange(-maxResidualY,maxResidualY);
	fHistallresidualsZ->SetAxisRange(-maxResidualZ,maxResidualZ);

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

// ---------------------------------------------------
// In order to be backward compatible
// ---------------------------------------------------
#if BACKWARD
void AliHLTTPCDisplay::DisplayClusters(Bool_t x3don,Float_t* etaRange) {
    if (!fc1){
	fc1 = new TCanvas("c1","",900,900);
	fc1->cd();
    }

    fSwitch3DTracks = kFALSE; 
    fSwitch3DCluster = kTRUE; 
    fSwitch3DPadRow = kFALSE; 
    fSwitch3DGeometry = kFALSE;

    Draw3D();
}
// ---------------------------------------------------
void AliHLTTPCDisplay::DisplayTracks(Int_t minhits,Bool_t x3don,Float_t thr) {
    if (!fc1){
	fc1 = new TCanvas("c1","",900,900);
	fc1->cd();
    }

    fMinHits = minhits; 
    fPtThreshold = thr;
    fSwitch3DTracks = kTRUE; 
    fSwitch3DCluster = kFALSE; 
    fSwitch3DPadRow = kFALSE; 
    fSwitch3DGeometry = kFALSE;

    Draw3D();
}
// ---------------------------------------------------
void AliHLTTPCDisplay::DisplayAll(Int_t minhits,Bool_t clusterswitch,Bool_t trackswitch,Bool_t x3don, Float_t thr, Float_t* etaRange){
    if (!fc1){
	fc1 = new TCanvas("c1","",900,900);
	fc1->cd();
    }

    fMinHits = minhits; 
    fPtThreshold = thr;
    fSwitch3DTracks = trackswitch; 
    fSwitch3DCluster = clusterswitch; 
    fSwitch3DPadRow = kFALSE; 
    fSwitch3DGeometry = kFALSE;

    Draw3D();
}
#endif
// ---------------------------------------------------
