// @(#) $Id$

/** \class AliHLTTPCDisplay
<pre>
//_____________________________________________________________
// AliHLTTPCDisplay
//
// Simple display class for the HLT tracks.
</pre>
*/
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//         Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>      
//*-- Copyright &copy ALICE HLT Group 


// Recent Changes:
// ==============
// - Displaying Padrows in Histograms and in the 3D Geometry, with the SetupPadRow, FillPadRow, DrawPadRow functions
// - Modification of the SetSlice, etc functions
// - Rewrite of the draw geometry functions





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

AliHLTTPCDisplay::AliHLTTPCDisplay()
{
  //constructor
  fGeom = NULL;
  fTracks = NULL;

  fc1 = new TCanvas("c1","",900,900);
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  memset(fNcl, 0, 36*6*sizeof(UInt_t));
  fBackColor = 1;
  fLineColor = 0;
  fSlicePair = -1;
  fTheta = 90.;

  fHistrawcl = NULL;
  fHistraw = NULL;
}

AliHLTTPCDisplay::AliHLTTPCDisplay(Int_t *slice,Char_t *gfile)
{
  //ctor. Specify which slices you want to look at.
  LoadGeometrie(gfile);
  if (slice) {
    SetSlices(slice[0], slice[1]);
  }

  fc1 = new TCanvas("c1","",900,900);
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  memset(fNcl, 0, 36*6*sizeof(UInt_t)); 
  fBackColor = 1;
  fLineColor = 0;
  fSlicePair = -1;
  fTheta = 90.;  

  fHistrawcl = NULL;
  fHistraw = NULL;
}

AliHLTTPCDisplay::~AliHLTTPCDisplay()
{
  //destructor
  if(fTracks)
    delete fTracks;
  if (fc1)
    delete fc1;
}

void AliHLTTPCDisplay::SetSlices(Int_t minslice, Int_t maxslice) {
  // set slice range
  fTheta = 90.;

  fMinSlice = minslice;
  fMaxSlice = maxslice;
  
  fSlicePair = -1;
}

void AliHLTTPCDisplay::SetSlices(Int_t slice) {
  // Set on slice
  fTheta = 0.;
  fMinSlice = slice;
  fMaxSlice = slice;
  fSlicePair = -1;
}

void AliHLTTPCDisplay::SetSlices() { 
  // Set all Slices
  fTheta = 90.;
  fMinSlice = 0;
  fMaxSlice = 35;
  fSlicePair = -1;
}

void AliHLTTPCDisplay::SetSlicesPair(Int_t slice) {
  // set pair of slices
  fTheta = 0.;
  fSlicePair = slice;
}

void AliHLTTPCDisplay::SetSlicesPair(Int_t minslice, Int_t maxslice) {
  // set range of pair of slices
  fTheta = 90.;
  fSlicePair = -2;
  fSlicePairMax = maxslice;
  fSlicePairMin = minslice;
}

void AliHLTTPCDisplay::SetInvert(Bool_t invert) {
  Int_t tmp;
  if (invert){ 
    tmp = fBackColor;
    fBackColor = fLineColor;
    fLineColor = tmp ;
  }
}


void AliHLTTPCDisplay::SetDrawGeo(Bool_t drawgeo) {
  if (drawgeo){
      if (fDrawGeo == kTRUE){
	  fDrawGeo = kFALSE;
      }
      else{
	  fDrawGeo = kTRUE;
      }
  }
}

// ########################################################################################################################################
void AliHLTTPCDisplay::DrawGeom(Int_t slice) {  
  Char_t fname[256];
  Int_t realslice = slice % 18;
  
  if (realslice < 10){
    sprintf(fname,"LS0%d",realslice);
    fGeom->GetNode(fname)->SetLineColor(fLineColor);
    fGeom->GetNode(fname)->Draw("same");
    sprintf(fname,"US0%d",realslice);
    fGeom->GetNode(fname)->SetLineColor(fLineColor); 
    fGeom->GetNode(fname)->Draw("same");
  }
  else {
    sprintf(fname,"LS%d",realslice);
    fGeom->GetNode(fname)->SetLineColor(fLineColor);
    fGeom->GetNode(fname)->Draw("same");
    sprintf(fname,"US%d",realslice);
    fGeom->GetNode(fname)->SetLineColor(fLineColor); 
    fGeom->GetNode(fname)->Draw("same");
  }   
}

// ########################################################################################################################################
Bool_t AliHLTTPCDisplay::LoadGeometrie(Char_t *gfile)
{
  if (gfile) {
  TFile *file = TFile::Open(gfile);
  if(!file)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::AliHLTTPCDisplay","File Open")
	<<"Geometry file " << gfile << " does not exist!"<<ENDLOG;
      return kFALSE;
    }
  
  fGeom = (TGeometry*)file->Get("AliceGeom");

  file->Close();
  delete file;
  }
  return kTRUE;
}

// ########################################################################################################################################
void AliHLTTPCDisplay::SetupClusterDataForPatch(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data)
{
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
      LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupClusterDataForPatch","memory allocation")
	<<"mmemory allocation failed "<<ENDLOG; 
    }
  } else {
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupClusterDataForPatch","argument check")
      <<"invalid argument "<<ENDLOG; 
  } 
}

// ########################################################################################################################################
void AliHLTTPCDisplay::Setup(Char_t *trackfile,Char_t *path,Int_t event,Bool_t sp)
{
  //Read in the hit and track information from produced files.
  
  Char_t fname[256];
  AliHLTTPCMemHandler *clusterfile[36][6];
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<AliHLTTPCTransform::GetNPatches(); p++)
	{
	  Int_t patch;
	  if(sp==kTRUE)
	    patch=-1;
	  else
	    patch=p;
	  clusterfile[s][p] = new AliHLTTPCMemHandler();
	  if(event<0)
	    sprintf(fname,"%s/points_%d_%d.raw",path,s,patch);
	  else
	    sprintf(fname,"%s/points_%d_%d_%d.raw",path,event,s,patch);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliHLTTPCLog::kError,"AliHLTTPCEvaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	      delete clusterfile[s][p];
              clusterfile[s][p] = 0; 
	      continue;
	    }
	  fClusters[s][p] = (AliHLTTPCSpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	  if(sp==kTRUE)
	    break;
	}
    }
//LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","SETUPTrack") <<" ========================= 1 " << ENDLOG;
  if(!trackfile) return;
//LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","SETUPTrack") <<" ========================= 2 " << ENDLOG;
  AliHLTTPCMemHandler *tfile = new AliHLTTPCMemHandler();
//LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","SETUPTrack") <<" ========================= 3 " << ENDLOG;
 if(!tfile->SetBinaryInput(trackfile)){
    
      LOG(AliHLTTPCLog::kError,"AliHLTTPCEvaluation::Setup","File Open")
	<<"Inputfile "<<trackfile<<" does not exist"<<ENDLOG; 
      return;
    }
 //LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","SETUPTrack") <<" ========================= 4 " << ENDLOG;
  fTracks = new AliHLTTPCTrackArray();
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
  //LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","SETUPTrack") <<" ========================= 5 " << ENDLOG;
}

// ########################################################################################################################################
void AliHLTTPCDisplay::DisplayAll(Int_t minhits,Bool_t clusterswitch,Bool_t trackswitch,Bool_t x3don, Float_t thr, Float_t* etaRange){
  //Display tracks & clusters

  // Set up canvas
  if (!fc1) return;
  fc1->cd();

  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  fc1->Clear();
  fc1->SetFillColor(fBackColor);
  fc1->SetTheta(fTheta);
  fc1->SetPhi(0.);
  
  if (clusterswitch){ // Display Clusters ---------------------------------------
      
      for(Int_t s=fMinSlice; s<=fMaxSlice; s++){
      
	  for(Int_t p=0;p<6;p++){
	  
	      AliHLTTPCSpacePointData *points = fClusters[s][p];
	      if(!points) continue;
	      Int_t npoints = fNcl[s][p];
	      TPolyMarker3D *pm = new TPolyMarker3D(npoints);
	      
	      Float_t xyz[3];
	      for(Int_t i=0; i<npoints; i++){
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
		  AliHLTTPCTransform::Local2Global(xyz,s);
		  
		  pm->SetPoint(i,xyz[0],xyz[1],xyz[2]); 
	      
	      }
	      pm->SetMarkerColor(2); 
	      pm->Draw("");
	  }
      }

      // ------------------- 
      Int_t minslice = fMinSlice + 18;
      Int_t maxslice = fMaxSlice + 18;
      if ( minslice < 36 && maxslice < 36) {
	  for(Int_t s=minslice; s<=maxslice; s++){
	      
	      for(Int_t p=0;p<6;p++){
		  
		  AliHLTTPCSpacePointData *points = fClusters[s][p];
		  if(!points) continue;
		  Int_t npoints = fNcl[s][p];
		  TPolyMarker3D *pm = new TPolyMarker3D(npoints);
		  
		  Float_t xyz[3];
		  for(Int_t i=0; i<npoints; i++){
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
		      AliHLTTPCTransform::Local2Global(xyz,s);
		      
		      pm->SetPoint(i,xyz[0],xyz[1],xyz[2]); 
		      
		  }
		  pm->SetMarkerColor(2); 
		  pm->Draw("");
	      }
	  }
      }

  } // END Display Clusters  ----------------------------------


  if (trackswitch){  // Display Tracks -----------------------------------------

    Int_t ntracks = fTracks->GetNTracks();

    TPolyLine3D *line = new TPolyLine3D[ntracks];
    TPolyLine3D *linetest = new TPolyLine3D[ntracks];
    //THelix *helix = new THelix[ntracks];

    Float_t xcl[176];
    Float_t ycl[176];
    Float_t zcl[176];

    for(Int_t j=0; j<ntracks; j++)
    { 

	AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
	if(!gtrack) continue;
	if((thr>0)&&(gtrack->GetPt()<thr)) continue;    

	// -------------------------- INSERT THELIX ----------------------------------------------------------
#if 0
/*	// TODO: check if IS local, then local to global

	Double_t xyz0[3];
	Double_t v0[3];


	// testmass is a parameter to fit the helix to the trackpoints --> has to be removed
	// but that is the problem right now 21 seems to be ok, for most of the tracks, but not for all
	// thats why i assume "mass"
	Double_t testmass = 21;

        xyz0[0] = gtrack->GetFirstPointX(); 
	xyz0[1] = gtrack->GetFirstPointY(); 
	xyz0[2] = gtrack->GetFirstPointZ();

      	v0[0] = testmass * gtrack->GetPx();
	v0[1] = testmass * gtrack->GetPy(); 
        v0[2] = testmass * gtrack->GetPz(); 

	Double_t xEnd = gtrack->GetLastPointX();       
	Double_t yEnd = gtrack->GetLastPointY();
	Double_t zEnd = gtrack->GetLastPointZ();

	Double_t tgl = gtrack->GetTgl();
	Double_t psi = gtrack->GetPsi();
	Double_t bField = AliHLTTPCTransform::GetBFieldValue();

        Int_t charge = gtrack->GetCharge();

	Double_t omega;	

	// used because bfield is not set yet in the Display class
	bField = 0.5;
	omega = charge * bField / testmass;
	
	Double_t hrange[2];
	hrange[0] = xyz0[2];
	hrange[1] = zEnd;

	Double_t rapidity = gtrack->GetRapidity();
	Int_t slice0 = gtrack->GetSector();

//if (j == 0){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","HELIX") << "local="<< gtrack->IsLocal() <<" rapidity="<<rapidity << " slice" << slice0 <<"\n"
	    << "x=" << xyz0[0] << " y=" << xyz0[1] << " z=" << xyz0[2] << " bfield=" << bField << "\n"
	    << "vx=" << v0[0] << " vy=" << v0[1] << " vz=" << v0[2] << " charge=" << charge << "\n"
	    << "x=" << xEnd  << " y=" << yEnd << " z=" << zEnd << " tgl=" << tgl << " psi=" << psi << "\n" <<  ENDLOG;


       	THelix *currenthelix = &(helix[j]);
	currenthelix = new THelix(xyz0,v0,omega,hrange,kHelixZ,0);
	currenthelix->SetLineColor(5);   
	currenthelix->SetLineWidth(3);
	currenthelix->Draw("same");
//	}
*/
#endif
	// --------------------------

	Int_t nHits = gtrack->GetNHits();
	UInt_t *hitnum = gtrack->GetHitNumbers();
	if(nHits < minhits) continue;
 	TPolyMarker3D *pm = new TPolyMarker3D(nHits);
	Int_t hitcount=0;

	Float_t xArr[2];
	Float_t yArr[2];
	Float_t zArr[2];

#if 0
/*	nHits = 50;
	Double_t dist =  xEnd - xyz0[0];
	Double_t deltax = dist / nHits;
	
	gtrack->CalculateHelix();
	Double_t pt = gtrack->GetPt();
	Double_t xcenter = gtrack->GetCenterX();       
	Double_t ycenter = gtrack->GetCenterY();
	Double_t radius = gtrack->GetRadius();


	for(Int_t h=0; h<5; h++) {
	    Double_t cx = xyz0[0] - xcenter;
	    Double_t cy = xyz0[1] - ycenter;

	    Double_t newx = cx + (h * deltax) ; 
	    Double_t newy =  cy - sqrt(radius+radius - h*h*deltax*deltax);
	    Double_t newz =  0.0; 
	    newx += xcenter;
	    newy += ycenter;
	    
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","HELIX") << "local="<< gtrack->IsLocal() << " slice" << slice0 <<"\n"
	    << "x=" << newx << " y=" << newy << " z=" << newz << " xcenter=" << xcenter << "\n"
	    << "dist=" << dist << " deltax=" << deltax << " radius=" << radius <<  ENDLOG;

	    pm->SetPoint(h,newx,newy,newz);
	    hitcount++;
	}
*/
#endif

	for(Int_t h=0; h<nHits; h++)
	  {
	    UInt_t id=hitnum[h];
	    Int_t slice = (id>>25) & 0x7f;
	    Int_t patch = (id>>22) & 0x7;
	    UInt_t pos = id&0x3fffff;	    
	    Int_t testslice;
	    if (slice > 17 &&  slice <36) testslice = slice - 17;
	    else testslice = slice;

	    if(testslice < fMinSlice || testslice > fMaxSlice)
	      continue;
	
	    AliHLTTPCSpacePointData *points = fClusters[slice][patch];

	    if(!points) {
	      LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","Clusterarray")
		<<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	      continue;
	    }
 
	    if(pos>=fNcl[slice][patch]) {
	      LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","Clusterarray")
		<<"Pos is too large: pos "<<pos <<" ncl "<<fNcl[slice][patch]<<ENDLOG;
	      continue;
	    }

	    Float_t xyztmp[3];
	    xyztmp[0] = points[pos].fX;
	    xyztmp[1] = points[pos].fY;
	    xyztmp[2] = points[pos].fZ;

	    AliHLTTPCTransform::Local2Global(xyztmp,slice);
	  	  
	    xcl[h] = xyztmp[0];
	    ycl[h] = xyztmp[1];
	    zcl[h] = xyztmp[2];

	    if (h == 0){
		xArr[0]=xyztmp[0];
		yArr[0]=xyztmp[1];
		zArr[0]=xyztmp[2];
	    }

	    Int_t maxH = nHits - 1; 

	    if (h == maxH){
		xArr[1]=xyztmp[0];
		yArr[1]=xyztmp[1];
		zArr[1]=xyztmp[2];
	    }
	    
	  
	    pm->SetPoint(h,xcl[h],ycl[h],zcl[h]);
	  
	    hitcount++;
	  }

	if(hitcount==0) continue;
	
	pm->SetMarkerColor(6); 
	//pm->Draw();

	TPolyLine3D *currentlinetest = &(linetest[j]);
	currentlinetest = new TPolyLine3D(2,xArr,yArr,zArr,"");
	currentlinetest->SetLineColor(3);   
	currentlinetest->SetLineWidth(3);
//	currentlinetest->Draw("same");

	TPolyLine3D *currentline = &(line[j]);
	currentline = new TPolyLine3D(nHits,xcl,ycl,zcl,"");
	currentline->SetLineColor(4);   
	currentline->SetLineWidth(2);
	currentline->Draw("same");
      }
  
LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","Track") <<" LEAVING TRACKS " << ENDLOG;
  } // END Display Tracks ------------------------------------------------

  DisplayGeom(x3don);
}

// ########################################################################################################################################
void AliHLTTPCDisplay::DisplayGeom(Bool_t x3don){
  // Draw GEO
  Int_t slice=0;
  Int_t slicepair[1];

  if (fSlicePair != -1){
    // Draw Pair of Slices
    slicepair[0] = fSlicePair;

    if (fSlicePair > 8) slicepair[1] = fSlicePair - 9;
    else slicepair[1] = fSlicePair + 9;

    DrawGeom(slicepair[0]);
    DrawGeom(slicepair[1]);

  }

  else if (fSlicePair != -1){
    for (slice=fSlicePairMin;slice<=fSlicePairMax;slice++){
      DrawGeom(slice);
    }
    // Draw range Pair of Slices
    if (fSlicePairMin > 8) {
      slicepair[0] = fSlicePairMin - 9;
      slicepair[1] = fSlicePairMax - 9;
    }
    else {
      slicepair[0] = fSlicePairMin + 9;
      slicepair[1] = fSlicePairMax + 9;
    }

    for (slice=slicepair[0];slice<=slicepair[1];slice++){
      DrawGeom(slice);
    }
    
  }
  else { 
    // Single Slice, or Range
      if (fMinSlice > fMaxSlice){
	  fMaxSlice += 17;
      }

      for (slice=fMinSlice;slice<=fMaxSlice;slice++){
	  if (fDrawGeo){
	      DrawGeom(slice);
	  }
      }
  }


  fc1->Draw();
  
  if(x3don) fc1->x3d(); 
  fc1->Modified();
  fc1->Update();


}

// ########################################################################################################################################
void AliHLTTPCDisplay::SetupPadRow(Int_t histSwitch, Int_t slice, Int_t padrow){

    fhistSwitch = histSwitch;
    fSlice = slice;
    fPadRow = padrow;

    fNPads = AliHLTTPCTransform::GetNPads(fPadRow);
    fNTimes = AliHLTTPCTransform::GetNTimeBins();

    if (!fc1) {
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupPadRow","Setup histograms") << "No valid canvas" << ENDLOG;
	return;
    }
    
    // ---  PADROW with HISTOGRAM "1"---
    if (fhistSwitch == 1){
	if ( fHistraw ){
	    delete fHistraw;
	    fHistraw = NULL;
	}
	if ( fHistrawcl ){
	    delete fHistrawcl;
	    fHistrawcl = NULL;
	}

 	// Setup the histograms
	Int_t padbinning = fNPads*10;
	fHistraw = new TH2F("fHistraw","Selected PadRow with found Clusters;Pad #;Timebin #",fNPads,0,fNPads-1,fNTimes,0,fNTimes-1);
	fHistrawcl = new TH1F("fHistrawcl","",padbinning,0,fNPads-1);

	fHistraw-> SetOption("COLZ"); 
    }



    gStyle->SetPalette(1);
}

// ########################################################################################################################################
void AliHLTTPCDisplay::FillPadRow(Int_t patch, ULong_t dataBlock, ULong_t dataLen, UInt_t nofClusters, AliHLTTPCSpacePointData* data){
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::FillPadRow","Enter") << ENDLOG;
    AliHLTTPCDigitReaderPacked* fDigitReader = new AliHLTTPCDigitReaderPacked();

    // reset histogram "1"
    if (fhistSwitch == 1){ 
	fHistraw->Reset();
	fHistrawcl->Reset();
    }

    bool readValue = true;
    Int_t rowOffset = 0;

    // Initialize RAW DATA
    Int_t firstRow = AliHLTTPCTransform::GetFirstRow(patch);
    Int_t lastRow = AliHLTTPCTransform::GetLastRow(patch);

    // Initialize block for reading packed data
    void* tmpdataBlock = (void*) dataBlock;
    fDigitReader->InitBlock(tmpdataBlock,dataLen,firstRow,lastRow);

    readValue = fDigitReader->Next();

    if (!readValue){	
	LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::FillPadRow","Read first value") << "No value in data block" << ENDLOG;
	return;
    }

    // Outer sector, patches 2, 3, 4, 5 -  start counting in patch 2 with row 0
    if ( patch >= 2 ) rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

    //  --- PADROW with GEOMETRY --- Initialize the colorbins
    if (fhistSwitch == 0){

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

    } // END --- PADROW with GEOMETRY --- Initialize the colorbins

    while ( readValue ){ 

	Int_t row = fDigitReader->GetRow() + rowOffset;

	// select padrow to fill in histogramm
	if (row == fPadRow){    

	    UChar_t pad = fDigitReader->GetPad();
	    UShort_t time = fDigitReader->GetTime();
	    UInt_t charge = fDigitReader->GetSignal();
	    Float_t xyz[3];
	   
            // --- PADROW with GEOMETRY ---
	    if (fhistSwitch == 0){
                // Transform raw coordinates to local coordinates
		AliHLTTPCTransform::RawHLT2Global(xyz, fSlice, fPadRow, pad, time);

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
	    }

	    // --- PADROW with HISTOGRAM ---
	    else {
		// Transform raw coordinates to local coordinates
		AliHLTTPCTransform::RawHLT2Local(xyz, fSlice, fPadRow, pad, time);

		fHistraw->Fill(pad,time,charge);
	    }
	}

	// read next value
	readValue = fDigitReader->Next();
      
	//Check where to stop:
	if(!readValue) break; //No more value
    } 
    
    if ( fDigitReader )
	delete fDigitReader;
    fDigitReader = NULL;

    //  --- PADROW with HISTOGRAM ---
    if (fhistSwitch != 0){
	// ========= FILL CLUSTER DATA 
	if (data && fSlice>=0 && fSlice<36 && patch>=0 && patch<AliHLTTPCTransform::GetNPatches()) {
	    if (fClusters[fSlice][patch]!=NULL) {
		delete(fClusters[fSlice][patch]);
		fClusters[fSlice][patch]=NULL;
	    }
	    
	    Int_t arraysize=nofClusters*sizeof(AliHLTTPCSpacePointData);
	    
	    fClusters[fSlice][patch] = (AliHLTTPCSpacePointData*)new Byte_t[arraysize];
	    
	    if (fClusters[fSlice][patch]) {
		memcpy(fClusters[fSlice][patch], data, arraysize);
		fNcl[fSlice][patch]=nofClusters;
	    } 
	    else {
		fNcl[fSlice][patch]=nofClusters;
		LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupClusterDataForPatch","memory allocation")
		    <<"mmemory allocation failed "<<ENDLOG; 
	    }
	} 

	else {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupClusterDataForPatch","argument check")
		<<"invalid argument "<<ENDLOG; 
	} 
	
	AliHLTTPCSpacePointData *points = fClusters[fSlice][patch];
	if(!points) return;
	Int_t npoints = fNcl[fSlice][patch];
	
	Float_t xyz[3];
	for(Int_t i=0; i<npoints; i++){
	    xyz[0] = points[i].fX;
	    xyz[1] = points[i].fY;
	    xyz[2] = points[i].fZ;
	    
	    Int_t clrow = AliHLTTPCTransform::GetPadRow(xyz[0]);
	    // select padrow to fill in histogramm
	    if (clrow == fPadRow){
		AliHLTTPCTransform::LocHLT2Raw(xyz, fSlice, fPadRow);
		fHistrawcl->Fill(xyz[1],xyz[2]);
	    }
	}
    } // END --- PADROW with HISTOGRAM ---
  LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::FillPadRow","LEAVE") << ENDLOG;
}

// ########################################################################################################################################
void AliHLTTPCDisplay::DrawPadRow(Bool_t x3don){  
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DrawPadRow","Enter") << ENDLOG;

//    fc1->Clear();

    //  --- PADROW with GEOMETRY ---
    if (fhistSwitch == 0){
	Int_t markercolor = 51;

	DisplayGeom(x3don);

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
    //  --- PADROW with HISTOGRAM ---
    else { 
	if (fBackColor == 1) SetInvert(kTRUE);
	fc1->SetFillColor(fBackColor);
	fHistraw->Draw("COLZ");

	fHistrawcl->SetMarkerStyle(28);
	fHistrawcl->SetMarkerSize(2);
	fHistrawcl->SetMarkerColor(1);
	
	fHistrawcl->Draw("psame");
    }


    fc1->Modified();
    fc1->Update();	
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DrawPadRow","LEAVE") << ENDLOG;
}
