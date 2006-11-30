// XEmacs -*-C++-*-
// @(#) $Id$
// Original: AliHLTDisplay.h,v 1.14 2004/06/30 08:30:14 loizides 

#ifndef ALIHLTTPCDISPLAY_H
#define ALIHLTTPCDISPLAY_H

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

#include <TGeometry.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraph.h>
#include "AliHLTLogging.h"

class AliHLTTPCSpacePointData;
class AliHLTTPCTrackArray;

class AliHLTTPCDisplay : public TObject, public AliHLTLogging {

 public:
    AliHLTTPCDisplay(Char_t *gfile="$(ALIHLT_BASEDIR)/geo/alice.geom") {InitDisplay(gfile);}
    virtual ~AliHLTTPCDisplay();

    void InitDisplay(Char_t *gfile);
  
    // SETUP
    void SetupHist();
    void SetupCluster(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data);
    void SetupTracks(AliHLTTPCTrackArray *tracks);
    
    // HISTOGRAM FILL/RESET FUNCTIONS
    void FillPadRow(Int_t patch, ULong_t dataBlock, ULong_t dataLen);
    void ResetHistPadRow();
    void ResetHistResiduals();
    void ResetHistCharge();

    // DRAWER
    void Draw3D();
    void DrawHistPadRow(); 
    void DrawGeomSector(Int_t sector);
    void DrawHistPad1();
    void DrawHistPad2();
    void DrawHistPad3();
    void DrawHistResiduals(Bool_t ySwitch);
    void DrawHistCharge();
    
    // SETTER  
    void SetSlices(){fMinSlice = 0; fMaxSlice = 35; fSlicePair = kFALSE; SetSliceArray();}
    void SetSlices(Int_t s){fMinSlice = s; fMaxSlice = s; fSlicePair = kFALSE; SetSliceArray();}
    void SetSlices(Int_t mins, Int_t maxs){fMinSlice = mins; fMaxSlice = maxs; fSlicePair = kFALSE; SetSliceArray();}
    void SetSlicesPair(Int_t s){fMinSlice = s; fMaxSlice = s; fSlicePair = kTRUE; SetSliceArray();}
    void SetSlicesPair(Int_t mins, Int_t maxs){fMinSlice = mins; fMaxSlice = maxs; fSlicePair = kTRUE; SetSliceArray();}

    void SetPad(Int_t f){fPad = f;}
    void SetPadRow(Int_t f){fPadRow = f;}
    void SetSlicePadRow(Int_t f){fSlicePadRow = f;}

    void SetMinHits(Int_t f){fMinHits = f;}
    void SetPtThreshold(Float_t f){fPtThreshold = f;}

    void SetSwitches(Bool_t f1, Bool_t f2, Bool_t f3, Bool_t f4) {fSwitch3DTracks = f1; fSwitch3DCluster = f2; fSwitch3DPadRow = f3; fSwitch3DGeometry = f4;} 
    void Set3DSwitchTracks(Bool_t f1) {fSwitch3DTracks = f1;}
    void Set3DSwitchCluster(Bool_t f1) {fSwitch3DCluster = f1;}
    void Set3DSwitchPadRow(Bool_t f1) {fSwitch3DPadRow = f1;}
    void Set3DSwitchGeometry(Bool_t f1) {fSwitch3DGeometry = f1;}

    void SetSelectTrack(Int_t f) {fSelectTrack = f;}
    void SetSelectTrackSlice(Int_t f) {fSelectTrackSlice = f;}
    void SetSelectTrackSwitch(Bool_t f) {fSelectTrackSwitch = f;}
    void SetSelectCluster(Int_t f) {fSelectCluster = f;}

    void SetHistPadRowAxis();

    void SetInvert() {Int_t tmp = fBackColor; fBackColor = fLineColor; fLineColor = tmp; }
    void SetKeepView(Bool_t f){fKeepView = f;} 

    // GETTER
    Int_t GetPadRow(){return fPadRow;}
    Int_t GetSlicePadRow(){return fSlicePadRow;}
    Int_t GetNPads(){return fNPads;}
    Int_t GetBackColor(){return fBackColor;}
    Bool_t GetKeepView(){return fKeepView;}
    
    Bool_t Get3DSwitchTracks() {return fSwitch3DTracks;}
    Bool_t Get3DSwitchCluster() {return fSwitch3DCluster;}
    Bool_t Get3DSwitchPadRow() {return fSwitch3DPadRow;}

    Bool_t GetSelectTrackSwitch() {return fSelectTrackSwitch;}

    // EXECUTER
    void ExecPadRow();

    struct AliHLTTPCTrackParameter{
	Int_t nHits;
	Int_t charge;
	Double_t kappa;
	Double_t radius;
	Double_t xyzF[3];
	Double_t xyzL[3];
	Int_t slice;
	Double_t phi0;
	Double_t psi;
	Double_t lambda;
	Double_t pt;
	Int_t id;
	Double_t bfield;
	Double_t s;
    };
    
    AliHLTTPCTrackParameter fTrackParam;

// ---------------------------------------------------
// In order to be backward compatible
// ---------------------------------------------------
/*
    AliHLTTPCDisplay(Int_t *slice, Char_t *gfile ) { InitDisplay(gfile);}

    // SETUP
    void SetupClusterDataForPatch(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data) {SetupCluster(slice, patch, nofClusters, data);}
    void SetTracks(AliHLTTPCTrackArray *tracks) { SetupTracks(tracks);}
    void DisplayClusters(Bool_t x3don=kTRUE,Float_t* etaRange=NULL);
    void DisplayTracks(Int_t minhits=10,Bool_t x3don=kTRUE,Float_t thr=0.);
    void DisplayAll(Int_t minhits=10,Bool_t clusterswitch=kTRUE,Bool_t trackswitch=kTRUE,Bool_t x3don=kTRUE, Float_t thr=0., Float_t* etaRange=NULL);
*/
// ---------------------------------------------------
 private:
// ---------------------------------------------------
// In order to be backward compatible
// ---------------------------------------------------
/*    TCanvas *fc1; */
// ---------------------------------------------------
    Bool_t LoadGeometrie(Char_t *gfile);
    void SetSliceArray();

    AliHLTTPCDisplay(const AliHLTTPCDisplay &/*d*/):TObject(){;}
    AliHLTTPCDisplay& operator=(const AliHLTTPCDisplay &/*d*/){return *this;}
    
    AliHLTTPCSpacePointData *fClusters[36][6]; 
    AliHLTTPCTrackArray *fTracks; 
    
    UInt_t fNcl[36][6];//number of cluster

    TH1F *fHistrawcl;  // histogram for cluster in padrow
    TH2F *fHistraw;    // histogram for signals in padrow
    TH1F *fHistpad1;   // histogram for pad in padrow
    TH1F *fHistpad2;   // histogram for pad in padrow
    TH1F *fHistpad3;   // histogram for pad in padrow
    TH1F *fHistallresidualsY;//histogram for all Y residuals
    TH1F *fHistallresidualsZ;//histogram for all Z residuals
    TH1F *fHistcharge; // histogram for clustercharge
    
    TGraph *fGraphresidualsY; // graph of the Y residuals for one track
    TGraph *fGraphresidualsZ; // graph of the Z residuals for one track
    TGraph *fGraphresidualsYLength; // graph of the Y residuals for one track
    TGraph *fGraphresidualsZLength; // graph of the Z residuals for one track


    TGeometry *fGeom;  // geometry
    Int_t fBackColor;  // Background color
    Int_t fLineColor;  // Line color
    Bool_t fKeepView;  // Keep View when redisplaying
    
    Int_t fPad;        // pad
    Int_t fPadRow;     // padrow
    Int_t fSlicePadRow;// slice where padrow is in
    Int_t fNPads;      // number of pads in padrow
    Int_t fNTimes;     // number of timebins
    Int_t fMinHits;    // minimum cluster per track
    Float_t fPtThreshold;// pt threshold for tracks
    
    Bool_t fSelectTrackSwitch;// switch ti single track mode
    Int_t fSelectTrack;// select single track
    Int_t fSelectTrackSlice; // select slice for single track
    
    Int_t fSelectCluster; // select all=0, used=1, unused=2 cluster
    
    Int_t fMinSlice;   //min slice
    Int_t fMaxSlice;   //max slice
    Bool_t fSlicePair; //draw pair of slices;
    Bool_t fSliceArray[36];//Array if slice should be drawn or not
    
    Bool_t fDrawGeo;
    Int_t fcolorbin[20]; // number of entries per colorbin
    Int_t fbinct[20];    // index of colorbin
    Float_t *fpmarr[20]; // contains point data
    
    Bool_t fSwitch3DCluster;
    Bool_t fSwitch3DTracks;
    Bool_t fSwitch3DPadRow;
    Bool_t fSwitch3DGeometry;

    ClassDef(AliHLTTPCDisplay,1) //Display class
};

#endif
