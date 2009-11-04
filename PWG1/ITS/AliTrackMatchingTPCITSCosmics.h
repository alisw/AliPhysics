#ifndef ALITRACKMATCHINGTPCITSCOSMICS_H
#define ALITRACKMATCHINGTPCITSCOSMICS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliTrackMatchingTPCITSCosmics
// AliAnalysisTask to study the matching efficiency of ITS and TPC
// Author: A.Dainese, andrea.dainese@pd.infn.it
//*************************************************************************

class TNtuple;
class TList;
class TH1F;

#include "AliAnalysisTask.h"

class AliTrackMatchingTPCITSCosmics : public AliAnalysisTask
{
 public:

  AliTrackMatchingTPCITSCosmics(const char *name="matching");
  virtual ~AliTrackMatchingTPCITSCosmics();


  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetOnlySPDFO(Bool_t set=kTRUE) {fOnlySPDFO=set;}  
  void SetReadHLTESD(Bool_t set=kTRUE) {fReadHLTESD=set;}  
  void SetGeometryFileName(TString name="geometry.root") {fGeometryFileName=name;}

 private:


  AliTrackMatchingTPCITSCosmics(const AliTrackMatchingTPCITSCosmics &source);
  AliTrackMatchingTPCITSCosmics& operator=(const AliTrackMatchingTPCITSCosmics& source); 

  Bool_t fOnlySPDFO;        // only fastOR events
  Bool_t fReadHLTESD;        // read the ESD from the HLT tree
  TString fGeometryFileName; // where to find the geometry.root
  AliESDEvent  *fESD;        // ESD object
  TList   *fList;    //! list of histos and ntuples: output slot 0
  TH1F    *fHistEvCount;     //! output histogram
  TNtuple *fntTrks;         //! output ntuple  
  TNtuple *fntPairs;        //! output ntuple  
  TNtuple *fntITSPairs;        //! output ntuple  

  ClassDef(AliTrackMatchingTPCITSCosmics,2); // AliAnalysisTask to check ITS-TPC matching
};

#endif

