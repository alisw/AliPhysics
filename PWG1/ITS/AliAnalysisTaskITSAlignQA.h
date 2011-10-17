#ifndef ALIANALYSISTASKITSALIGNQA
#define ALIANALYSISTASKITSALIGNQA

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysiTaskITSAlignQA
// AliAnalysisTaskSE to extract from ESD + ESDfriends 
// the track-to-point residuals and dE/dx vs, time for SDD modules
//
// Author: F. Prino, prino@to.infn.it
//*************************************************************************

class TList;
class TH1F;
class TH2F;
class TProfile;
class TTree;
class TString;
class AliESDEvent;
class AliESDfriend;
class AliITSTPArrayFit;
class AliTrackPointArray;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskITSAlignQA : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskITSAlignQA();
  virtual ~AliAnalysisTaskITSAlignQA();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetDoSPDResiduals(Bool_t opt){
    fDoSPDResiduals=opt;
  }
  void SetDoSDDResiduals(Bool_t opt){
    fDoSDDResiduals=opt;
  }
  void SetDoSSDResiduals(Bool_t opt){
    fDoSSDResiduals=opt;
  }
  void SetDoSDDdEdxCalib(Bool_t opt){
    fDoSDDdEdxCalib=opt;
  }
  void SetDoSDDVDriftCalib(Bool_t opt){
    fDoSDDVDriftCalib=opt;
  }
  void SetDoSDDDriftTime(Bool_t opt){
    fDoSDDDriftTime=opt;
  }
  void SetDoAllResiduals(){
    fDoSPDResiduals=kTRUE;
    fDoSDDResiduals=kTRUE;
    fDoSSDResiduals=kTRUE;
  }
  void SetDoAll(){
    SetDoAllResiduals();
    fDoSDDdEdxCalib=kTRUE;    
  }
  void SetUseITSstandaloneTracks(Bool_t use){
    fUseITSsaTracks=use;
  }
  void SetLoadGeometryFromOCDB(Bool_t opt){
    fLoadGeometry=opt;
  }

  void SetMinITSPoints(Int_t minp=3){
    fMinITSpts=minp;
  }
  void SetMinTPCPoints(Int_t minp=70){
    fMinTPCpts=minp;
  }
  void SetMinPt(Float_t minpt=1.0){
    fMinPt=minpt;
  }
  void SetMinVtxContributors(Int_t n=5)     { fMinVtxContributors = n; }
  void SetUseVertex(Bool_t v=kTRUE)         { fUseVertex = v; }
  void SetUseVertexForZOnly(Bool_t v=kTRUE) { fUseVertexForZOnly = v; } // Use the vertex for SDD Z residuals only
  void SetRemovePileupWithSPD(Bool_t opt=kTRUE) { fRemovePileupWithSPD = opt; }
  
  void     SetOCDBInfo(UInt_t runNb, const char *location) {
    fRunNb=runNb; 
    fOCDBLocation=location;
  }

  Bool_t   AcceptTrack(const AliESDtrack * track);
  Bool_t   AcceptVertex(const AliESDVertex * vtx, const AliESDVertex * vtxSPD);
  void     CreateSPDHistos();
  void     CreateSDDHistos();
  void     CreateSSDHistos();

  void     FitAndFillSPD(Int_t iLayer, const AliTrackPointArray *array, Int_t npts, AliESDtrack * track);
  void     FitAndFillSDDrphi(const AliTrackPointArray *array, Int_t npts, AliESDtrack * track);
  void     FitAndFillSDDz(Int_t iLayer, const AliTrackPointArray *array, Int_t npts, AliESDtrack * track);
  void     FitAndFillSSD(Int_t iLayer, const AliTrackPointArray *array, Int_t npts, AliESDtrack * track);
  void     SetPtBinLimits(Int_t nBins, Double_t* xbins){
    fNPtBins=nBins;
    if(nBins>kMaxPtBins) fNPtBins=kMaxPtBins;
    for(Int_t iBin=0; iBin<=fNPtBins; iBin++) fPtBinLimits[iBin]=xbins[iBin];
  }
  void     LoadGeometryFromOCDB();
  AliTrackPointArray* PrepareTrack(const AliTrackPointArray* inp, const AliESDVertex* vtx=0);
  void                PrepareVertexConstraint(const AliESDVertex* vtx, AliTrackPoint &point);
 private:
  AliAnalysisTaskITSAlignQA(const AliAnalysisTaskITSAlignQA &source);
  AliAnalysisTaskITSAlignQA& operator=(const AliAnalysisTaskITSAlignQA &source);
  
  enum {kNSPDmods = 240};
  enum {kNSDDmods = 260};
  enum {kNSSDmods = 1698};
  enum {kMaxPtBins = 12};
  enum {kVtxSensVID=14371};    // dummy VID for "vertex" point

  TList* fOutput;              //! Histos with residuals
  TH1F*  fHistNEvents;         //! histo with N of events  
  TH1F*  fHistPtAccept;        //! histo of pt distribution of accepted tracks 

  TH2F*  fHistSPDResidX[kNSPDmods];       //! histos of SPD residuals along Xloc vs. Pt
  TH2F*  fHistSPDResidZ[kNSPDmods];       //! histos of SPD residuals along Zloc vs. Pt
  TH2F*  fHistSDDResidX[kNSDDmods];       //! histos of SDD residuals along Xloc vs. Pt
  TH2F*  fHistSDDResidZ[kNSDDmods];       //! histos of SDD residuals along Zloc vs. Pt
  TH2F*  fHistSSDResidX[kNSSDmods];       //! histos of SSD residuals along Xloc vs. Pt
  TH2F*  fHistSSDResidZ[kNSSDmods];       //! histos of SSD residuals along Zloc vs. Pt

  TH2F*  fHistSDDResidXvsX[kNSDDmods];    //! histos of SDD residuals along Xloc vs. Xloc
  TH2F*  fHistSDDResidXvsZ[kNSDDmods];    //! histos of SDD residuals along Xloc vs. Zloc
  TH2F*  fHistSDDResidZvsX[kNSDDmods];    //! histos of SDD residuals along Zloc vs. Xloc
  TH2F*  fHistSDDResidZvsZ[kNSDDmods];    //! histos of SDD residuals along Zloc vs. Zloc
  TH2F*  fHistSDDdEdxvsDrTime[kNSDDmods]; //! histos of SDD dE/dx vs. drift time
  TH1F*  fHistSDDDrTimeAll[kNSDDmods];    //! histos of SDD drift time (all clusters)
  TH1F*  fHistSDDDrTimeExtra[kNSDDmods];  //! histos of SDD drift time (extra clusters)
  TH1F*  fHistSDDDrTimeAttac[kNSDDmods];  //! histos of SDD drift time (attached clusters)
  //
  // RS
  TProfile* fHProfSDDResidXvsXD[kNSDDmods][2]; // ! profile histos of SDD residuals along Xloc vs. Drift distance, each side separately
  TProfile* fHProfSDDDrTimevsXD[kNSDDmods][2]; // ! profile histos of SDD drift time vs. Drift distance, each side separately
  TProfile* fHProfSDDResidXvsZ[kNSDDmods][2];  // ! profile histos of SDD residuals along Xloc vs. Z (anode), each side separately
  TProfile* fHProfSDDDrTimevsZ[kNSDDmods][2];  // ! profile histos of SDD drift time vs. Z (anode), each side separately
  //
  Bool_t   fDoSPDResiduals;   // Flag to enable histos of SPD residuals
  Bool_t   fDoSDDResiduals;   // Flag to enable histos of SDD residuals
  Bool_t   fDoSSDResiduals;   // Flag to enable histos of SSD residuals
  Bool_t   fDoSDDdEdxCalib;   // Flag to enable histos for SDD dE/dx calibration
  Bool_t   fDoSDDVDriftCalib; // Flag to enable histos for SDD VDrift calibration
  Bool_t   fDoSDDDriftTime;   // Flag to enable histos for SDD Drift times
  Bool_t   fUseITSsaTracks;   // Flag for using standalone ITS tracks
  Bool_t   fLoadGeometry;     // Flag to control the loading of geometry from OCDB
  Bool_t   fUseVertex;        // Use the vertex as an extra point
  Bool_t   fUseVertexForZOnly; // Use the vertex for SDD Z residuals only
  Int_t    fMinVtxContributors; // min N contributors to accept vertex if fUseVertex is on
  Bool_t   fRemovePileupWithSPD; // Use/not use pileup rejection with SPD
  Int_t    fMinITSpts;        // Minimum number of ITS points per track
  Int_t    fMinTPCpts;        // Minimum number of TPC points per track
  Float_t  fMinPt;            // Minimum pt to accept tracks
  Int_t    fNPtBins;          // number of pt bins
  Double_t fPtBinLimits[kMaxPtBins+1];  // limits of Pt bins

  AliITSTPArrayFit* fFitter;  // Track Point fitter
  Int_t fRunNb;               // Run number
  TString fOCDBLocation;      // OCDB location

  ClassDef(AliAnalysisTaskITSAlignQA,4);
};


#endif


