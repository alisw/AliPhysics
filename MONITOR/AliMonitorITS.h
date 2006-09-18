#ifndef ALIMONITORITS_H
#define ALIMONITORITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMonitor.h"

class AliITSgeom;


class AliMonitorITS : public AliMonitor {
public:
  AliMonitorITS(AliITSgeom* param);
  virtual ~AliMonitorITS() {};

  virtual void     CreateHistos(TFolder* folder);
  virtual void     FillHistos(AliRunLoader* runLoader, 
			      AliRawReader* rawReader, AliESD* esd);

private:
  AliMonitorITS(const AliMonitorITS& monitor);
  AliMonitorITS& operator = (const AliMonitorITS& monitor);

  AliITSgeom*      fGeom;               // ITS geometry

  AliMonitorHisto* fSDDDigitsCharge;    // charge distribution of ITS-SDD digits
  AliMonitorHisto* fSSDDigitsCharge;    // charge distribution of ITS-SSD digits
  AliMonitorHisto* fSDDClustersCharge;  // charge distribution of ITS-SDD clusters
  AliMonitorHisto* fSSDClustersCharge;  // charge distribution of ITS-SSD clusters
  AliMonitorHisto* fSPDNClustersVsModule;  // mean number of ITS-SPD clusters per module
  AliMonitorHisto* fSDDNClustersVsModule;  // mean number of ITS-SDD clusters per module
  AliMonitorHisto* fSSDNClustersVsModule;  // mean number of ITS-SSD clusters per module
  AliMonitorHisto* fNClustersVsLayer;   // mean number of ITS clusters per layer
  AliMonitorHisto* fNTracks;            // number of ITS tracks per event
  AliMonitorHisto* fNTracksITSTPC;      // correlation of number of ITS and TPC tracks per event
  AliMonitorHisto* fTrackPt;            // pt distribution of ITS tracks
  AliMonitorHisto* fTrackEta;           // eta distribution of ITS tracks
  AliMonitorHisto* fTrackPhi;           // phi distribution of ITS tracks
  AliMonitorHisto* fTrackDEdxVsP;       // dE/dx vs momentum distribution of TPC tracks

  ClassDef(AliMonitorITS, 0)   // creation and filling of monitor histograms for ITS
};
 

#endif









