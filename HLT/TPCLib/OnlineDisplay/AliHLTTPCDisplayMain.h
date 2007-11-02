// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAYMAIN_H
#define ALIHLTTPCDISPLAYMAIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCDisplayMain.h
    @author Jochen Thaeder
    @date   
    @brief  Interface class for ALICE HLT Online Display
*/

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

/**
 * @class AliHLTTPCDisplayMain
 * The class handels the interface between the graphical user interface (AliLHLTGUI) 
 * and the worker classes for the ALICE HLT online display. It also handles the interface
 * to the AnalysisChain via HOMER reader class. Furthermore all relavant global variables 
 * are stored here via get and set functions.
 * @ingroup alihlt_tpc
 */

class AliHLTTPCDisplayMain : public TObject , public AliHLTLogging {

 public:  
  /** standard constructor */
  //  AliHLTTPCDisplayMain();

  /**
   * Constructor
   * @param pt2GUI          Pointer to class AliHLTGUI
   * @param pt2Function     Pointer to callback Function class AliHLTGUI
   */
  AliHLTTPCDisplayMain(void* pt2GUI, void (*pt2Function)(void*, Int_t));

  /** standard destructor */
  virtual ~AliHLTTPCDisplayMain();

  /** Methods for connecting to data sources, reading data and writeing data */
  /** ---------------------------------------------------------------------- */

  //-->   TODO check this file location
  /** Connect to Analysis Chain via hosts and Ports specified in AliHLTGUI 
   *  @param cnt          Number hosts
   *  @param hostnames    Array of hostnames of the TDS
   *  @param ports        Array of ports of the TDS
   *  @param gfile        Path to the TPC geometry file
   */
  Int_t Connect( unsigned int cnt, const char** hostnames, unsigned short* ports,Char_t *gfile="$(ALIHLT_TOPDIR)/build/share/alice-hlt/alice.geom" );

  /** Disconnet from the Analysis Chain */
  Int_t Disconnect();
  
  /** Reads next event and call reader for raw, clsuter and track data 
   *  @param nextSwitch  Switch for reading the next event.<br> 
                         Standard is KTRUE.
   */
  Int_t ReadData(Bool_t nextSwitch = kTRUE);

  /** Read raw data,do zero suppression and check for occupancy limit. Afterwards fill in two arrays (raw and zerosuppressed). */
  void ReadRawData();

  /** Read cluster data and fill cluster structure. */
  void ReadClusterData();

  /** Read track data and fill in track structure. */
  void ReadTrackData();

  /** Main working function. Calls worker classes and handles displaying 
   *  @param newRawSlice   If set to kTRUE, raw data from a new slice is read and zero suppression is applied.<br>
   *                       Standard is kFALSE.
   */
  void DisplayEvent(Bool_t newRawSlice = kFALSE);

  /** Save all histograms as eps files **/
  void SaveHistograms();



  /** Setup methods */
  /** ------------- */

  /** Fill cluster data in structure
   *  @param slice         Slice ID of datablock where data is in.
   *  @param patch         Patch ID of datablock where data is in.
   *  @param nofClusters   Number of clusters.in data block.
   *  @param data          Pointer to structure where cluster data is written.
   */
   void SetupCluster(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data);

  /** Fill track data in structure */
  void SetupTracks();

  /** Register methods */
  /** ---------------- */

  /** Canvas types
   */
  enum AliHLTTPCDisplayCanvasType_t {
    Ccharge,
    Cpadrow,
    Cpad,
    CthreeD,
    Cresiduals,
    Cfront,
    Chits_s,
    Cq_track,
    Cq_s,
    Cpadrow_pad,
    nCanvasTypes
  };

  /** Worker types
   */
  enum AliHLTTPCDisplayWorkerType_t {
    Wcharge,
    Wpadrow,
    Wpad,
    WthreeD,
    Wresiduals,
    Wfront,
    nWorkerTypes
  };


  /** Register canvas 
   *  @param canvasType    @ref AliHLTTPCDisplayCanvasType_t  
   *  @param canvas        Pointer to TCanvas.
   */
  void RegisterCanvas(AliHLTTPCDisplayCanvasType_t canvasType, TCanvas * canvas ){
    fCanvasArray[canvasType] = canvas;
  }

  /** Register worker
   *  @param workerType    @ref AliHLTTPCDisplayWorkerType_t  
   *  @param worker        Pointer to worker class.
   */
  void RegisterWorker(AliHLTTPCDisplayWorkerType_t workerType, void * worker ){
    fWorkerArray[workerType] = worker;
  }

  /** Get pointer to canvas 
   *  @param  canvasType    @ref AliHLTTPCDisplayCanvasType_t 
   *  @return Returns pointer to "canvasType" (@ref AliHLTTPCDisplayCanvasType_t) or NULL if fails.
   */
  TCanvas * GetCanvas(AliHLTTPCDisplayCanvasType_t canvasType){
    if ( !fCanvasArray[canvasType] )
      return NULL;
    else
      return fCanvasArray[canvasType];
  }

  /** Get pointer to worker 
   *  @param  workerType    @ref AliHLTTPCDisplayWorkerType_t 
   *  @return Returns pointer to "workerType" (@ref AliHLTTPCDisplayWorkerType_t) or NULL if fails.
   */
  void * GetWorker(AliHLTTPCDisplayWorkerType_t workerType){
    if ( !fWorkerArray[workerType] )
      return NULL;
    else
      return fWorkerArray[workerType];
  }



  void SetRawReaderMode(Int_t f) {fRawReaderMode = f;}
  void SetZeroSuppressionThreshold(Int_t f) {fZeroSuppressionThreshold = f;}
  void SetOccupancyLimit(Float_t f) {fOccupancyLimit = f;}
  void SetBField(Float_t f) {fBField = f;}
  void SetNTimeBins(Int_t f) {fNTimeBins = f;}

  Int_t GetRawReaderMode(){return fRawReaderMode;}
  Int_t GetZeroSuppressionThreshold(){return fZeroSuppressionThreshold;}
  Float_t GetOccupancyLimit() {return fOccupancyLimit;}
  Float_t GetBField(){return fBField;}
  Int_t GetNTimeBins() {return fNTimeBins;}

  /** Set methods */
  /** ----------- */


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
  void SetCanvasPadRow_Pad(TCanvas *f){fCanvasPadRow_Pad = f;}


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
  
  // 3D display ----
  void SetInvert() {Int_t tmp = fBackColor; fBackColor = fLineColor; fLineColor = tmp; }
  void SetKeepView(Bool_t f){fKeepView = f;} 
  void SetSwitches(Bool_t f1, Bool_t f2, Bool_t f3, Bool_t f4) {fSwitch3DTracks = f1; fSwitch3DCluster = f2; fSwitch3DPadRow = f3; fSwitch3DGeometry = f4;} 
  void Set3DRawSwitch(Int_t f) {fSwitch3DRaw = f;}
  void SetTheta(Float_t f) {fTheta = f;}

  // raw ----
  void SetZeroSuppression(Bool_t f){fZeroSuppression = f;}
  //     -- pad/padrow
  void SetPad(Int_t f){fPad = f;}
  void SetPadRow(Int_t f){fPadRow = f; fNPads = AliHLTTPCTransform::GetNPads(f); }
  void SetSlicePadRow(Int_t f){fSlicePadRow = f;}
  void SetSplitPadRow(Bool_t f){fSplitPadRow = f;}
  //     -- front
  void SetFrontDataSwitch(Int_t f){fFrontDataSwitch = f;}
  void SetTimeBinMinMax(Int_t f1, Int_t f2){fTimeBinMin = f1; fTimeBinMax = f2; }
  void SetSplitFront(Bool_t f){fSplitFront = f;}

  // track ----
  void SetSelectTrack(Int_t f) {fSelectTrack = f;}
  void SetSelectTrackSlice(Int_t f) {fSelectTrackSlice = f;}
  void SetSelectTrackSwitch(Bool_t f) {fSelectTrackSwitch = f;}
  void SetSelectCluster(Int_t f) {fSelectCluster = f;}

  /*
   * **********************************
   *              GETTER
   * **********************************
   */

    
    Bool_t GetConnectionStatus() {return fConnect;}
    ULong64_t GetEventID() {return fEventID;}

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
  
    // raw ----
    Bool_t GetZeroSuppression(){return fZeroSuppression;}
    //     -- pad/padrow
    Int_t GetPad(){return fPad;}
    Int_t GetNPads(){return fNPads;}
    Int_t GetPadRow(){return fPadRow;}
    Int_t GetSlicePadRow(){return fSlicePadRow;}
    Bool_t GetSplitPadRow(){return fSplitPadRow;}
    //     -- front
  //    Int_t GetNTimeBins(){return fgNTimeBins;}
    Int_t GetFrontDataSwitch(){return fFrontDataSwitch;}
    Int_t GetTimeBinMin(){return fTimeBinMin;}
    Int_t GetTimeBinMax(){return fTimeBinMax;}
    Bool_t GetSplitFront(){return fSplitFront;}
  


    Int_t GetLineColor(){return fLineColor;}
    Int_t GetBackColor(){return fBackColor;}
    Bool_t GetKeepView(){return fKeepView;}
    Float_t GetTheta(){return fTheta;}
    Float_t GetPhi(){return fPhi;}

    Bool_t Get3DSwitchTracks() {return fSwitch3DTracks;}
    Bool_t Get3DSwitchCluster() {return fSwitch3DCluster;}
    Bool_t Get3DSwitchPadRow() {return fSwitch3DPadRow;}
    Bool_t Get3DSwitchGeometry() {return fSwitch3DGeometry;}
    Int_t  Get3DSwitchRaw() {return fSwitch3DRaw;}

    Int_t GetTrackParamNHits(){ return fTrackParam.nHits;}
    /*    
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
    */

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
    
    void* fReader;                 // really HOMERReader*

    UInt_t fRawData[159][140][1024];                    // Raw Data of one Slice
    UInt_t fRawDataZeroSuppressed[159][140][1024];      // Raw Data of one Slice zero suppressed

    AliHLTTPCTrackParameter fTrackParam;

  TCanvas **fCanvasArray;
  void **fWorkerArray;


// ---------------------------------------------------
 private:
  /** copy constructor prohibited */
  AliHLTTPCDisplayMain(const AliHLTTPCDisplayMain&);
  /** assignment operator prohibited */
  AliHLTTPCDisplayMain& operator=(const AliHLTTPCDisplayMain&);

  //    AliHLTTPCDisplayMain(const AliHLTTPCDisplayMain &/*d*/):TObject(){;}
  //    AliHLTTPCDisplayMain& operator=(const AliHLTTPCDisplayMain &/*d*/){return *this;}
      
    void SetSliceArray();          // Fill Array with slices which 

  // ** global constants **
  Int_t fgNTimeBins;             // Number of TimeBins
  Int_t fRawReaderMode;           // raw reader mode
  Int_t fZeroSuppressionThreshold;
  Float_t fOccupancyLimit;
  Float_t fBField;
  Int_t fNTimeBins;


    // **  HOMER parameter / connection / data exist **
    ULong64_t fEventID;            // Event ID

    Bool_t fConnect;               // Connection status
  
    Bool_t fExistsRawData;         // Raw data present
    Bool_t fExistsClusterData;     // Cluster data present
    Bool_t fExistsTrackData;       // Track data present

    // ** pointer to display classes **
    AliHLTTPCDisplayCharge* fDisplayCharge; 
    AliHLTTPCDisplayPadRow* fDisplayPadRow;
    AliHLTTPCDisplayPad* fDisplayPad;
    AliHLTTPCDisplay3D* fDisplay3D;
    AliHLTTPCDisplayResiduals* fDisplayResiduals;
    AliHLTTPCDisplayFront * fDisplayFront;

    // ** pointer to canvases in GUI ** 
    TCanvas * fCanvasCharge;
    TCanvas * fCanvasPadRow;
    TCanvas * fCanvasPad;
    TCanvas * fCanvas3D;
    TCanvas * fCanvasResiduals;
    TCanvas * fCanvasFront;
    TCanvas * fCanvasHits_S;
    TCanvas * fCanvasQ_Track;
    TCanvas * fCanvasQ_S;
    TCanvas * fCanvasPadRow_Pad;

    Int_t fTracksPerSlice[36];     // TrackCount per slice

    // ** cluster / tarck container **
    AliHLTTPCSpacePointData *fClusters[36][6]; 
    AliHLTTPCTrackArray *fTracks; 
    
    UInt_t fNcl[36][6];            // Number of cluster

    // ** selected selected slice(s) **
    Int_t fMinSlice;               // Min slice
    Int_t fMaxSlice;               // Max slice
    Bool_t fSlicePair;             // Pair of slices;
    Bool_t fSliceArray[36];        // Array if slice should be drawn or not

    // ** select type of clusters **
    Int_t fSelectCluster;          // select all=0, used=1, unused=2 cluster

    // ** raw data variables **
    Bool_t fZeroSuppression;       // enable zero suppression
    //     -- pad/padrow
    Int_t fPad;                    // pad
    Int_t fNPads;                  // number of pads in row
    Int_t fPadRow;                 // padrow
    Int_t fSlicePadRow;            // slice where padrow is in
    Bool_t fSplitPadRow;           // Split PadRow Canvas
    //     -- front
    Int_t fFrontDataSwitch;        // select average/sum/maximum
    Int_t fTimeBinMin;             // min TimeBin
    Int_t fTimeBinMax;             // max TimeBin
    Bool_t fSplitFront;            // Split Front Canvas

  
 
    // ** select tracks **
    Bool_t fSelectTrackSwitch;     // switch ti single track mode
    Int_t fSelectTrack;            // select single track
    Int_t fSelectTrackSlice;       // select slice for single track


    // ** cuts on tracks **
    Int_t fCutHits;
    Float_t fCutPt;
    Float_t fCutPsi;
    Float_t fCutLambda;
    Float_t fCut_S;
    Int_t fCutPadrow;

    // ** "keepview", angles and colors for 3D view **
    Bool_t fKeepView;              // Keep View when redisplaying

    Float_t fTheta;
    Float_t fPhi;

    Int_t fBackColor;              // Background color
    Int_t fLineColor;              // Line color

    // ** 3D switches **
    Bool_t fSwitch3DCluster;
    Bool_t fSwitch3DTracks;
    Bool_t fSwitch3DPadRow;
    Bool_t fSwitch3DGeometry;
    Int_t fSwitch3DRaw;    

    ClassDef(AliHLTTPCDisplayMain,0) //Display class
};
#endif
