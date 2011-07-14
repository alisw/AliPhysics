#ifndef ALIANALYSISMUMUFROMESD_H
#define ALIANALYSISMUMUFROMESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// 
/// AliAnalysisMuMuFromESD : implementation of AliAnalysisMuMu 
/// when reading data from ESD
/// 
/// author : Laurent Aphecetche (Subatech)

#ifndef ALIANALYSISMUMU_H
#  include "AliAnalysisMuMu.h"
#endif
#ifndef ALIESDVZERO_H
#  include "AliESDVZERO.h"
#endif

class AliESDEvent;
class TH1;
class TH2;
class AliESDMuonTrack;
class AliESDVertex;

class AliAnalysisMuMuFromESD : public AliAnalysisMuMu
{
public:
  AliAnalysisMuMuFromESD();
  AliAnalysisMuMuFromESD(TList* triggerClassesToConsider);
  AliAnalysisMuMuFromESD(Bool_t aa);
  virtual ~AliAnalysisMuMuFromESD();

protected:
  virtual void MuUserExec(Option_t *option);

private:
  
  AliAnalysisMuMuFromESD(const AliAnalysisMuMuFromESD&);
  AliAnalysisMuMuFromESD& operator=(const AliAnalysisMuMuFromESD&);

  UInt_t GetTrackMask(const AliESDMuonTrack& track) const;

  void FillHistogramCollection(const char* physics, const char* triggerClassName);
  
  void FillHistosForTrack(const char* physics, const char* triggerClassName, const char* centrality, const AliESDMuonTrack& track, const char* runNumber);
  
  void FillHistos(const char* physics, const char* triggerClassName, const char* centrality, const AliESDEvent& esd);

  Double_t CorrectedDCA(const AliESDMuonTrack& track) const;

  Double_t PDCACutValue(const AliESDMuonTrack& track) const;

private:
  void Ctor();
  const char* RunNumber(const AliESDEvent& esd) const;
  Bool_t TrackMatchCut(const AliESDMuonTrack& track) const;
  Bool_t TrackMatchLowCut(const AliESDMuonTrack& track) const;
  Bool_t TrackMatchHighCut(const AliESDMuonTrack& track) const;
  Bool_t TrackRabsCut(const AliESDMuonTrack& track) const;
  Bool_t TrackPtCut(const AliESDMuonTrack& track) const;
  Bool_t TrackChi2(const AliESDMuonTrack& track) const;
  Bool_t TrackEtaCut(const AliESDMuonTrack& track) const;
  Bool_t TrackDCACut(const AliESDMuonTrack& track) const;
    
  AliESDVertex* fVertex; //! current event vertex
  
  ClassDef(AliAnalysisMuMuFromESD,5) // Analysis of mu-mu pairs from ESD
};

#endif
