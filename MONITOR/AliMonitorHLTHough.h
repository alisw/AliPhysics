#ifndef ALIMONITORHLTHOUGH_H
#define ALIMONITORHLTHOUGH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliMonitor.h"

class AliTPCParam;


class AliMonitorHLTHough : public AliMonitor {
public:
  AliMonitorHLTHough(AliTPCParam* param);
  AliMonitorHLTHough(const AliMonitorHLTHough& monitor);
  AliMonitorHLTHough& operator = (const AliMonitorHLTHough& monitor);
  virtual ~AliMonitorHLTHough();

  virtual void     CreateHistos(TFolder* folder);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader);

private:
  AliTPCParam*     fParam;              // TPC parameters

  AliMonitorTrend* fNTracks;            // number of HLT tracks per event
  AliMonitorHisto* fTrackPt;            // pt distribution of HLT tracks
  AliMonitorHisto* fTrackEta;           // eta distribution of HLT tracks
  AliMonitorHisto* fTrackPhi;           // phi distribution of HLT tracks
  AliMonitorHisto* fTrackEtaVsPhi;      // phi vs eta for HLT tracks
  AliMonitorHisto* fPtEtaVsPhi;         // phi vs eta for HLT tracks

  ClassDef(AliMonitorHLTHough, 0)   // creation and filling of monitor histograms for HLT Hough transform
};
 

#endif









