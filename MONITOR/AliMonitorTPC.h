#ifndef ALIMONITORTPC_H
#define ALIMONITORTPC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitor.h"
#include "AliMonitorHisto.h"
#include "AliMonitorTrend.h"
#include "AliTPCParam.h"


class AliMonitorDataTPC : public TObject {
public:
  AliMonitorDataTPC();
  AliMonitorDataTPC(Int_t size);
  virtual ~AliMonitorDataTPC();
  void     SetSize(Int_t size);

  Int_t    fNTracks;   // number of TPC tracks
  Float_t* fPt;        //[fNTracks]
  Float_t* fEta;       //[fNTracks]
  Float_t* fPhi;       //[fNTracks]

private:
  Int_t    fSize;      //! size of the arrays

  ClassDef(AliMonitorDataTPC, 1)   // data structure for the TPC monitor tree branch
};


class AliMonitorTPC : public AliMonitor {
public:
  AliMonitorTPC(AliTPCParam* param);
  virtual ~AliMonitorTPC();

  virtual void     CreateHistos(TFolder* folder);
  virtual void     CreateBranches(TTree* tree);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader);

private:
  AliTPCParam*     fParam;              // TPC parameters

  AliMonitorHisto* fPadsCharge;         // charge distribution of TPC pads
  AliMonitorHisto* fClustersCharge;     // charge distribution of TPC clusters
  AliMonitorHisto* fNClustersVsRow;     // mean number of TPC clusters per pad row
  AliMonitorHisto* fNClustersVsSector;  // mean number of TPC clusters per sector
  AliMonitorTrend* fNTracks;            // number of TPC tracks per event
  AliMonitorHisto* fTrackPt;            // pt distribution of TPC tracks
  AliMonitorHisto* fTrackEta;           // eta distribution of TPC tracks
  AliMonitorHisto* fTrackPhi;           // phi distribution of TPC tracks
  AliMonitorHisto* fTrackDEdxVsP;       // dE/dx vs momentum distribution of TPC tracks

  AliMonitorDataTPC* fData;             // data for the monitor tree

  ClassDef(AliMonitorTPC, 0)   // creation and filling of monitor histograms for TPC
};
 

#endif









