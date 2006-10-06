// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAYMAIN_H
#define ALIHLTTPCDISPLAYMAIN_H
/** \class AliHLTTPCDisplayMain
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayMain
//
// Display class for the HLT events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TGeometry.h>
#include <TObject.h>
#include <TCanvas.h>

#include "AliHLTTPCTransform.h"
#include "AliHLTLogging.h"

class AliHLTTPCSpacePointData;
class AliHLTTPCTrackArray;
class AliHLTTPCDisplayCharge;
class AliHLTTPCDisplayPadRow;
class AliHLTTPCDisplayPad;
class AliHLTTPCDisplay3D;
class AliHLTTPCDisplayResiduals;
class AliHLTTPCDisplayFront;

class AliHLTTPCDisplayMain : public TObject , public AliHLTLogging {

 public:
    AliHLTTPCDisplayMain(void* pt2GUI, void (*pt2Function)(void*, Int_t));
    virtual ~AliHLTTPCDisplayMain();

    Int_t Connect( unsigned int cnt, const char** hostnames, unsigned short* ports,Char_t *gfile="$(ALIHLT_BASEDIR)/geo/alice.geom" );
    
    void ReadData(Bool_t nextSwitch = kTRUE);
    void DisplayEvent();
    void SaveHistograms();

    // SETUP
    void SetupCluster(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data);
    void SetupTracks();

    // SETTER   
    void SetConnectionStatus(Bool_t f) {fConnect = f;}

    // canvas ----
    void SetCanvasCharge(TCanvas *f){fCanvasCharge = f;}
    void SetCanvasPadRow(TCanvas *f){fCanvasPadRow = f;}
    void SetCanvasPad(TCanvas *f){fCanvasPad = f;}
    void SetCanvas3D(TCanvas *f){fCanvas3D = f;}
    void SetCanvasResiduals(TCanvas *f){fCanvasResiduals = f;}
    void SetCanvasFront(TCanvas *f){fCanvasFront = f;}
    void SetCanvasHits_S(TCanvas *f){fCanvasHits_S = f;}
    void SetCanvasQ_Track(TCanvas *f){fCanvasQ_Track = f;}
    void SetCanvasQ_S(TCanvas *f){fCanvasQ_S = f;}

    // slices ----
    void SetSlices(){fMinSlice = 0; fMaxSlice = 35; fSlicePair = kFALSE; SetSliceArray();}
    void SetSlices(Int_t s){fMinSlice = s; fMaxSlice = s; fSlicePair = kFALSE; SetSliceArray();}
    void SetSlices(Int_t mins, Int_t maxs){fMinSlice = mins; fMaxSlice = maxs; fSlicePair = kFALSE; SetSliceArray();}
    void SetSlicesPair(Int_t s){fMinSlice = s; fMaxSlice = s; fSlicePair = kTRUE; SetSliceArray();}
    void SetSlicesPair(Int_t mins, Int_t maxs){fMinSlice = mins; fMaxSlice = maxs; fSlicePair = kTRUE; SetSliceArray();}

    // cuts ----
    //    void SetCutHits(Int_t f){fMinHits = f;}
    //    void SetCutPt(Float_t f){fPtThreshold = f;}
    void SetCutHits(Int_t f){fCutHits = f;}
    void SetCutPt(Float_t f){fCutPt = f;}
    void SetCutPsi(Float_t f){fCutPsi = f;}
    void SetCutLambda(Float_t f){fCutLambda = f;}
    void SetCutS(Float_t f){fCut_S = f;}
    void SetIncidentPadrow(Int_t f){fCutPadrow = f;}

    void SetInvert() {Int_t tmp = fBackColor; fBackColor = fLineColor; fLineColor = tmp; }
    void SetKeepView(Bool_t f){fKeepView = f;} 

    // raw ----
    void SetPad(Int_t f){fPad = f;}
    void SetPadRow(Int_t f){fPadRow = f; fNPads = AliHLTTPCTransform::GetNPads(f); }
    void SetSlicePadRow(Int_t f){fSlicePadRow = f;}
    void SetTimebin(Int_t f){fTimebin = f;}

    void SetAllTimebins(Bool_t f){fAllTimebins = f;}
    void SetSplitPadRow(Bool_t f){fSplitPadRow = f;}

    void SetSelectTrack(Int_t f) {fSelectTrack = f;}
    void SetSelectTrackSlice(Int_t f) {fSelectTrackSlice = f;}
    void SetSelectTrackSwitch(Bool_t f) {fSelectTrackSwitch = f;}
    void SetSelectCluster(Int_t f) {fSelectCluster = f;}

    void SetSwitches(Bool_t f1, Bool_t f2, Bool_t f3, Bool_t f4) {fSwitch3DTracks = f1; fSwitch3DCluster = f2; fSwitch3DPadRow = f3; fSwitch3DGeometry = f4;} 

    void SetTheta(Float_t f) {fTheta = f;}


    // GETTER
    Bool_t GetConnectionStatus() {return fConnect;}

    AliHLTTPCSpacePointData* GetSpacePointDataPointer(Int_t slice,Int_t patch){return fClusters[slice][patch];}
    AliHLTTPCTrackArray* GetTrackArrayPointer() {return fTracks;}
    AliHLTTPCDisplayPadRow * GetPadRowPointer() {return fDisplayPadRow;}
    AliHLTTPCDisplayPad * GetPadPointer() {return fDisplayPad;}

    // canvas ----
    TCanvas * GetCanvasHits_S(){return fCanvasHits_S;}
    TCanvas * GetCanvasQ_Track(){return fCanvasQ_Track;}
    TCanvas * GetCanvasQ_S(){return fCanvasQ_S;}
    TCanvas * GetCanvasResiduals(){return fCanvasResiduals;}
    TCanvas * GetCanvasCharge(){return fCanvasCharge;}
    TCanvas * GetCanvasPadRow(){return fCanvasPadRow;}
    TCanvas * GetCanvasPad(){return fCanvasPad;}
    TCanvas * GetCanvas3D(){return fCanvas3D;}
    TCanvas * GetCanvasFront(){return fCanvasFront;}

    // cuts ----
    Int_t GetCutHits(){return fCutHits;}
    Float_t GetCutPt(){return fCutPt;}
    Float_t GetCutPsi(){return fCutPsi;}
    Float_t GetCutLambda(){return fCutLambda;}
    Float_t GetCutS(){return fCut_S;}
    Int_t GetIncidentPadrow(){return fCutPadrow;}

    Int_t GetNumberSpacePoints(Int_t slice,Int_t patch){return fNcl[slice][patch];}
    Int_t GetTracksPerSlice(Int_t slice){return fTracksPerSlice[slice];}
    Bool_t GetDisplaySlice(Int_t slice){ return fSliceArray[slice];}

    Int_t GetSelectCluster(){return fSelectCluster;}
    Bool_t GetSelectTrackSwitch() {return fSelectTrackSwitch;}
    Int_t GetSelectTrack() {return fSelectTrack;}
    Int_t GetSelectTrackSlice() {return fSelectTrackSlice;}

    Int_t GetGlobalTrack(Int_t slice);

    Int_t GetPad(){return fPad;}
    Int_t GetPadRow(){return fPadRow;}
    Int_t GetSlicePadRow(){return fSlicePadRow;}
    Bool_t GetSplitPadRow(){return fSplitPadRow;}

    Int_t GetNPads(){return fNPads;}
    Int_t GetTimebin(){return fTimebin;}
    Bool_t GetAllTimebins(){return fAllTimebins;}

    Int_t GetLineColor(){return fLineColor;}
    Int_t GetBackColor(){return fBackColor;}
    Bool_t GetKeepView(){return fKeepView;}
    Float_t GetTheta(){return fTheta;}
    Float_t GetPhi(){return fPhi;}

    Bool_t Get3DSwitchTracks() {return fSwitch3DTracks;}
    Bool_t Get3DSwitchCluster() {return fSwitch3DCluster;}
    Bool_t Get3DSwitchPadRow() {return fSwitch3DPadRow;}
    Bool_t Get3DSwitchGeometry() {return fSwitch3DGeometry;}

    Int_t GetTrackParamNHits(){ return fTrackParam.nHits;}

#if 0
    Char_t GetTrackParamSlice(){ return Char_t GetTrackParamId();
    Char_t GetTrackParamKappa());
    Char_t GetTrackParamPt());
Char_t GetTrackParamNHits());
Char_t GetTrackParamCharge());
Char_t GetTrackParamRadius());
Char_t GetTrackParamPhi0());
Char_t GetTrackParamPsi());
Char_t GetTrackParamLambda());
Char_t GetTrackParamBfield());
#endif
    // EXISTS
    Bool_t ExistsRawData() {return fExistsRawData;}   
    Bool_t ExistsClusterData() {return fExistsClusterData;}   
    Bool_t ExistsTrackData() {return fExistsTrackData;}   

    // EVENTS
    void ExecPadEvent(Int_t event, Int_t x, Int_t y, TObject *selected);

    // Callback Handler
    void * fPt2Gui;
    void (*fPadCallback)(void*, Int_t);

    struct AliHLTTPCTrackParameter{
	Int_t nHits;
	Int_t charge;
	Double_t kappa;
	Double_t radius;
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
 private:
    AliHLTTPCDisplayMain(const AliHLTTPCDisplayMain &/*d*/):TObject(), AliHLTLogging(){;}
    AliHLTTPCDisplayMain& operator=(const AliHLTTPCDisplayMain &/*d*/){return *this;}
   
    void SetSliceArray();

    Bool_t fConnect;               // Connection status
  
    Bool_t fExistsRawData;         // Raw data present
    Bool_t fExistsClusterData;     // Cluster data present
    Bool_t fExistsTrackData;       // Track data present

    void* fReader;                 // really HOMERReader*

    Int_t fTracksPerSlice[36];     // TrackCount per slice

    AliHLTTPCSpacePointData *fClusters[36][6]; 
    AliHLTTPCTrackArray *fTracks; 
    
    UInt_t fNcl[36][6];            // Number of cluster

    Int_t fMinSlice;               // Min slice
    Int_t fMaxSlice;               // Max slice
    Bool_t fSlicePair;             // Pair of slices;
    Bool_t fSliceArray[36];        // Array if slice should be drawn or not

    Int_t fSelectCluster;          // select all=0, used=1, unused=2 cluster

    // cuts ----
//    Int_t fMinHits;                // minimum cluster per track
//    Float_t fPtThreshold;          // pt threshold for tracks
    Int_t fCutHits;
    Float_t fCutPt;
    Float_t fCutPsi;
    Float_t fCutLambda;
    Float_t fCut_S;
    Int_t fCutPadrow;

    Int_t fBackColor;              // Background color
    Int_t fLineColor;              // Line color
    Bool_t fKeepView;              // Keep View when redisplaying

    Bool_t fSwitch3DCluster;
    Bool_t fSwitch3DTracks;
    Bool_t fSwitch3DPadRow;
    Bool_t fSwitch3DGeometry;

    Int_t fPad;                    // pad
    Int_t fPadRow;                 // padrow
    Int_t fSlicePadRow;            // slice where padrow is in
    Int_t fNPads;                  // number of pads in padrow
    Int_t fTimebin;

    Bool_t fAllTimebins;
    Bool_t fSplitPadRow;

    Bool_t fSelectTrackSwitch;     // switch ti single track mode
    Int_t fSelectTrack;            // select single track
    Int_t fSelectTrackSlice;       // select slice for single track

    Float_t fTheta;
    Float_t fPhi;

    AliHLTTPCDisplayCharge* fDisplayCharge; 
    AliHLTTPCDisplayPadRow* fDisplayPadRow;
    AliHLTTPCDisplayPad* fDisplayPad;
    AliHLTTPCDisplay3D* fDisplay3D;
    AliHLTTPCDisplayResiduals* fDisplayResiduals;
    AliHLTTPCDisplayFront * fDisplayFront;

    TCanvas * fCanvasCharge;
    TCanvas * fCanvasPadRow;
    TCanvas * fCanvasPad;
    TCanvas * fCanvas3D;
    TCanvas * fCanvasResiduals;
    TCanvas * fCanvasFront;
    TCanvas * fCanvasHits_S;
    TCanvas * fCanvasQ_Track;
    TCanvas * fCanvasQ_S;

    ClassDef(AliHLTTPCDisplayMain,0) //Display class
};
#endif
