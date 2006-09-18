#ifndef ALIMONITORTPC_H
#define ALIMONITORTPC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitor.h"

class AliTPCParam;
class AliMonitorDataTPC;


class AliMonitorTPC : public AliMonitor {
public:
  AliMonitorTPC(AliTPCParam* param);
  virtual ~AliMonitorTPC();

  virtual void     CreateHistos(TFolder* folder);
  virtual void     CreateBranches(TTree* tree);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader, AliESD* esd);

private:
  AliMonitorTPC(const AliMonitorTPC& monitor);
  AliMonitorTPC& operator = (const AliMonitorTPC& monitor);

  AliTPCParam*     fParam;              // TPC parameters

  AliMonitorHisto* fPadsCharge;         // charge distribution of TPC pads
  AliMonitorHisto* fClustersCharge;     // charge distribution of TPC clusters
  AliMonitorHisto* fNClustersVsRow;     // mean number of TPC clusters per pad row
  AliMonitorHisto* fNClustersVsSector;  // mean number of TPC clusters per sector
  AliMonitorTrend* fNTracks;            // number of TPC tracks per event
  AliMonitorHisto* fTrackPt;            // pt distribution of TPC tracks
  AliMonitorHisto* fTrackEta;           // eta distribution of TPC tracks
  AliMonitorHisto* fTrackPhi;           // phi distribution of TPC tracks
  AliMonitorHisto* fTrackNCl;           // number of clusters per track
  AliMonitorHisto* fTrackDEdxVsP;       // dE/dx vs momentum distribution of TPC tracks
  AliMonitorHisto* fTrackDEdx;          // dE/dx distribution of TPC tracks for a given momentum region
  AliMonitorHisto* fTrackEtaVsPhi;      // phi vs eta for TPC tracks
  AliMonitorHisto* fPtEtaVsPhi;      // phi vs eta for TPC tracks

  AliMonitorDataTPC* fData;             // data for the monitor tree

  ClassDef(AliMonitorTPC, 0)   // creation and filling of monitor histograms for TPC
};
 

#endif









