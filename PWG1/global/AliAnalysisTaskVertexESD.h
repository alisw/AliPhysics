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

#include "AliAnalysisTask.h"

class AliAnalysisTaskVertexESD : public AliAnalysisTask 
{
 public:

  AliAnalysisTaskVertexESD(const char *name = "AliAnalysisTaskVertexESD");
  virtual ~AliAnalysisTaskVertexESD(); 
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t         GetReadMC() const { return fReadMC; }
  void           SetReadMC(Bool_t flag=kTRUE) { fReadMC=flag; }
  void           SetRerecoVertexTPC(Bool_t flag=kTRUE) { fRecoVtxTPC=flag; }
  void           SetRerecoVertexITSTPC(Bool_t flag=kTRUE) { fRecoVtxITSTPC=flag; }
  
 protected:
  Bool_t       fReadMC;         // read Monte Carlo
  Bool_t       fRecoVtxTPC;     // reco TPC vertex on the flight
  Bool_t       fRecoVtxITSTPC;  // reco ITS+TPC vertex on the flight
  AliESDEvent *fESD;            // ESD object
  TList       *fOutput;         //! list send on output slot 0
  TNtuple     *fNtupleVertexESD;//! output ntuple
  TH2F        *fhTrackRefs;     //! output histo
  TH1F *fhVtxTPCx;  //! output histo
  TH1F *fhVtxTPCy;  //! output histo
  TH1F *fhVtxTPCz;  //! output histo
  TH1F *fhVtxTRKx;  //! output histo
  TH1F *fhVtxTRKy;  //! output histo
  TH1F *fhVtxTRKz;  //! output histo
  TH1F *fhVtxSPDx;  //! output histo
  TH1F *fhVtxSPDy;  //! output histo
  TH1F *fhVtxSPDz3D;  //! output histo
  TH1F *fhVtxSPDzZ;  //! output histo
  TH2F *fhVtxTPCvsSPDx;  //! output histo
  TH2F *fhVtxTPCvsSPDy;  //! output histo
  TH2F *fhVtxTPCvsSPDz;  //! output histo
  TH2F *fhVtxTRKvsSPDx;  //! output histo
  TH2F *fhVtxTRKvsSPDy;  //! output histo
  TH2F *fhVtxTRKvsSPDz;  //! output histo
  TH2F *fhVtxSPDContrvsMult;  //! output histo
  TH2F *fhVtxTRKContrvsTrks56;  //! output histo
  TH2F *fhVtxTPCContrvsTrks;  //! output histo

 private:    

  AliAnalysisTaskVertexESD(const AliAnalysisTaskVertexESD&); // not implemented
  AliAnalysisTaskVertexESD& operator=(const AliAnalysisTaskVertexESD&); // not implemented
  AliESDVertex* ReconstructPrimaryVertexTPC() const;
  AliESDVertex* ReconstructPrimaryVertexITSTPC() const;
  
  ClassDef(AliAnalysisTaskVertexESD,3); // primary vertex analysis
};

#endif
