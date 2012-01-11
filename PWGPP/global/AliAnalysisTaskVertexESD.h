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
class AliVEvent;

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
  void           SetFillNtupleBeamSpot(Bool_t fillBeamSpot=kFALSE){fFillTreeBeamSpot = fillBeamSpot;}
  void           SetTriggerType(AliVEvent::EOfflineTriggerTypes triggerType) {fTriggerType = triggerType;}

 protected:
  Bool_t       fCheckEventType; // read only events of type 7
  Bool_t       fReadMC;         // read Monte Carlo
  Bool_t       fRecoVtxTPC;     // reco TPC vertex on the flight
  Bool_t       fRecoVtxITSTPC;  // reco ITS+TPC vertex on the flight
  Bool_t       fRecoVtxITSTPCHalfEvent;  // reco ITS+TPC vertex with even and odd tracks
  Bool_t       fOnlyITSTPCTracks; // only ITS-TPC tracks to redo ITSTPC vertex
  Bool_t       fOnlyITSSATracks;  // only ITS-SA tracks to redo ITSTPC vertex
  Bool_t       fFillNtuple;      // fill ntuple  
  Bool_t       fFillTreeBeamSpot; //beam spot info in a tree
  AliESDEvent *fESD;            // ESD object
  TList       *fOutput;         //! list send on output slot 0
  TNtuple     *fNtupleVertexESD;//! output ntuple
  TH1F        *fhSPDVertexX; //! output histo
  TH1F        *fhSPDVertexY; //! output histo
  TH1F        *fhSPDVertexZ; //! output histo
  TH1F        *fhSPDVertexZonly; //! output histo
  TH1F        *fhTRKVertexX; //! output histo
  TH1F        *fhTRKVertexY; //! output histo
  TH1F        *fhTRKVertexZ; //! output histo
  TH1F        *fhTPCVertexX; //! output histo
  TH1F        *fhTPCVertexY; //! output histo
  TH1F        *fhTPCVertexZ; //! output histo
  TH2F        *fhTrackRefs;     //! output histo
  TTree       *fTreeBeamSpot;  //! output tree beam spot

  TH1F        *fhTriggeredTrklets; //! output histo
  TH1F        *fhSPD3DTrklets; //! output histo
  TH1F        *fhSPDZTrklets; //! output histo
  TH1F        *fhTRKTrklets; //! output histo
  TH1F        *fhTRKcTrklets; //! output histo
  TH1F        *fhTRKncTrklets; //! output histo
  TH1F        *fhTPCTrklets; //! output histo
  TH1F        *fhTPCcTrklets; //! output histo
  TH1F        *fhTPCncTrklets; //! output histo
  TH1F        *fhSPD3DZreco; //! output histo
  TH1F        *fhSPDZZreco; //! output histo


  TH1F        *fhSPDVertexXPile; //! output histo
  TH1F        *fhSPDVertexYPile; //! output histo
  TH1F        *fhSPDVertexZPile; //! output histo
  TH1F        *fhSPDVertexDiffZPileContr2; //! output histo
  TH1F        *fhSPDVertexDiffZPileContr3; //! output histo
  TH1F        *fhSPDVertexDiffZPileContr4; //! output histo
  TH1F        *fhSPDVertexDiffZPileContr5; //! output histo
  TH1F        *fhSPDVertexDiffZPileDefault; //! output histo
  TH1F        *fhSPDContributorsPile; //! output histo
  TH2F        *fhSPDDispContributors; //! output histo
  AliVEvent::EOfflineTriggerTypes    fTriggerType; //flag to set trigger type
	
  TH2F        *fhntrksSPDvsSPDcls;//!  output histo correlation contributors vs number of clustrers spd 
  TH2F        *fhntrksZvsSPDcls; //! output histo correlation contributors vs number of clustrers spd

 private:    

  AliAnalysisTaskVertexESD(const AliAnalysisTaskVertexESD&); // not implemented
  AliAnalysisTaskVertexESD& operator=(const AliAnalysisTaskVertexESD&); // not implemented
  AliESDVertex* ReconstructPrimaryVertexTPC(Bool_t constr=kFALSE) const;
  AliESDVertex* ReconstructPrimaryVertexITSTPC(Bool_t constr=kFALSE,Int_t mode=0) const;
  
  ClassDef(AliAnalysisTaskVertexESD,11); // primary vertex analysis
};

#endif
