#ifndef ALIMONITORHLT_H
#define ALIMONITORHLT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitor.h"

class AliTPCParam;


class AliMonitorHLT : public AliMonitor {
public:
  AliMonitorHLT(AliTPCParam* param);
  AliMonitorHLT(const AliMonitorHLT& monitor);
  AliMonitorHLT& operator = (const AliMonitorHLT& monitor);
  virtual ~AliMonitorHLT() {};

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
  AliMonitorHisto* fTrackNHits;         // number of hits per HLT track
  AliMonitorHisto* fTrackDEdxVsP;       // dedx distribution of HLT tracks
  AliMonitorHisto* fTrackDz0;           // dz0 distribution of HLT tracks
  AliMonitorHisto* fTrackDr0;           // dr0 distribution of HLT tracks
  AliMonitorHisto* fTrackAngle;         // azimutal distribution of HLT tracks

  ClassDef(AliMonitorHLT, 0)   // creation and filling of monitor histograms for HLT
};
 

#endif









