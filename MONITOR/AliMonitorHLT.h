#ifndef ALIMONITORHLT_H
#define ALIMONITORHLT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitor.h"
#include "AliMonitorHisto.h"
#include "AliMonitorTrend.h"
#include "AliTPCParam.h"


class AliMonitorHLT : public AliMonitor {
public:
  AliMonitorHLT(AliTPCParam* param);
  virtual ~AliMonitorHLT();

  virtual void     CreateHistos(TFolder* folder);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader);

private:
  AliTPCParam*     fParam;              // TPC parameters

  AliMonitorHisto* fClustersCharge;     // charge distribution of HLT clusters
  AliMonitorHisto* fNClustersVsRow;     // mean number of HLT clusters per pad row
  AliMonitorHisto* fNClustersVsSector;  // mean number of HLT clusters per sector
  AliMonitorTrend* fNTracks;            // number of HLT tracks per event
  AliMonitorHisto* fTrackPt;            // pt distribution of HLT tracks
  AliMonitorHisto* fTrackEta;           // eta distribution of HLT tracks
  AliMonitorHisto* fTrackPhi;           // phi distribution of HLT tracks

  ClassDef(AliMonitorHLT, 0)   // creation and filling of monitor histograms for HLT
};
 

#endif









