#ifndef AliAnalysisTaskVertexESD_cxx
#define AliAnalysisTaskVertexESD_cxx

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskVertexESD
// AliAnalysisTask to extract from ESD the information for the analysis
// of the primary vertex reconstruction efficiency and resolution
// (for MC events) and distributions (for real data). Three vertices:
// - SPD tracklets
// - ITS+TPC tracks
// - TPC-only tracks
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
//*************************************************************************

class TNtuple;
class TH1F;
class TH2F;
class AliESDEvent;
class AliESDVertex;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskVertexESD : public AliAnalysisTaskSE 
{
 public:

  AliAnalysisTaskVertexESD(const char *name = "AliAnalysisTaskVertexESD");
  virtual ~AliAnalysisTaskVertexESD(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void           SetCheckEventType(Bool_t check=kTRUE) {fCheckEventType=check;}
  Bool_t         GetReadMC() const { return fReadMC; }
  void           SetReadMC(Bool_t flag=kTRUE) { fReadMC=flag; if(flag) fCheckEventType=kFALSE;}
  void           SetRerecoVertexTPC(Bool_t flag=kTRUE) { fRecoVtxTPC=flag; }
  void           SetRerecoVertexITSTPC(Bool_t flag=kTRUE) { fRecoVtxITSTPC=flag; }
  void           SetRerecoVertexITSTPCHalfEvent(Bool_t flag=kTRUE) { fRecoVtxITSTPCHalfEvent=flag; }
  void           SetOnlyITSTPCTracks() {fOnlyITSTPCTracks=kTRUE;}
  void           SetOnlyITSSATracks() {fOnlyITSSATracks=kTRUE;}
  void           SetFillNtuple(Bool_t fill=kTRUE) {fFillNtuple=fill;}  
  void           SetFillNtupleBeamSpot(Bool_t fillBeamSpot=kFALSE){fFillNtupleBeamSpot = fillBeamSpot;}

 protected:
  Bool_t       fCheckEventType; // read only events of type 7
  Bool_t       fReadMC;         // read Monte Carlo
  Bool_t       fRecoVtxTPC;     // reco TPC vertex on the flight
  Bool_t       fRecoVtxITSTPC;  // reco ITS+TPC vertex on the flight
  Bool_t       fRecoVtxITSTPCHalfEvent;  // reco ITS+TPC vertex with even and odd tracks
  Bool_t       fOnlyITSTPCTracks; // only ITS-TPC tracks to redo ITSTPC vertex
  Bool_t       fOnlyITSSATracks;  // only ITS-SA tracks to redo ITSTPC vertex
  Bool_t       fFillNtuple;      // fill ntuple 
  Bool_t       fFillNtupleBeamSpot; //beam spot info 
  AliESDEvent *fESD;            // ESD object
  TList       *fOutput;         //! list send on output slot 0
  TNtuple     *fNtupleVertexESD;//! output ntuple
  TH1F        *fhSPDVertexX; //! output histo
  TH1F        *fhSPDVertexY; //! output histo
  TH1F        *fhSPDVertexZ; //! output histo
  TH1F        *fhTRKVertexX; //! output histo
  TH1F        *fhTRKVertexY; //! output histo
  TH1F        *fhTRKVertexZ; //! output histo
  TH1F        *fhTPCVertexX; //! output histo
  TH1F        *fhTPCVertexY; //! output histo
  TH1F        *fhTPCVertexZ; //! output histo
  TH2F        *fhTrackRefs;     //! output histo
  TNtuple     *fNtupleBeamSpot; //! output ntuple beam spot  

 private:    

  AliAnalysisTaskVertexESD(const AliAnalysisTaskVertexESD&); // not implemented
  AliAnalysisTaskVertexESD& operator=(const AliAnalysisTaskVertexESD&); // not implemented
  AliESDVertex* ReconstructPrimaryVertexTPC(Bool_t constr=kFALSE) const;
  AliESDVertex* ReconstructPrimaryVertexITSTPC(Bool_t constr=kFALSE,Int_t mode=0) const;
  
  ClassDef(AliAnalysisTaskVertexESD,8); // primary vertex analysis
};

#endif
