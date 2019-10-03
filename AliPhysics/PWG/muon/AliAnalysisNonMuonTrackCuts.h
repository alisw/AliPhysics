#ifndef ALIANALYSISNONMUONTRACKCUTS_H
#define ALIANALYSISNONMUONTRACKCUTS_H

/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */

#include "AliAnalysisCuts.h"

/// \brief Class to filter out tracks which are not (forward) muon tracks
/// \author L. Aphecetche (Subatech)
///

class AliAnalysisNonMuonTrackCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNonMuonTrackCuts();
  virtual ~AliAnalysisNonMuonTrackCuts() {}
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  
  ClassDef(AliAnalysisNonMuonTrackCuts,1); // Select muon spectrometer tracks
};

#endif