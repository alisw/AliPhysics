#ifndef ALIMONITORHLTHOUGH_H
#define ALIMONITORHLTHOUGH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliMonitor.h"

class AliTPCParam;


class AliMonitorHLTHough : public AliMonitor {
public:
  AliMonitorHLTHough(AliTPCParam* param);
  virtual ~AliMonitorHLTHough() {};

  virtual void     CreateHistos(TFolder* folder);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader, AliESD* esd);

private:
  AliMonitorHLTHough(const AliMonitorHLTHough& monitor);
  AliMonitorHLTHough& operator = (const AliMonitorHLTHough& monitor);

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
  AliMonitorHisto* fTrackDEdx;          // dedx distribution of HLT tracks for a given momentum region
  AliMonitorHisto* fTrackEtaVsPhi;      // phi vs eta for HLT tracks
  AliMonitorHisto* fPtEtaVsPhi;         // phi vs eta for HLT tracks

  ClassDef(AliMonitorHLTHough, 0)   // creation and filling of monitor histograms for HLT Hough transform
};
 

#endif









