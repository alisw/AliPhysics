#ifndef ALICFTREEMAPPING_H
#define ALICFTREEMAPPING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Class for mapping the CF trees of lightweight events
//
// Author: Igor Lakomov <Igor.Lakomov@cern.ch>
//

#include "TObject.h"

class AliCFTreeMapping : public TObject {
 public:
  AliCFTreeMapping();
  virtual ~AliCFTreeMapping();

  enum {kStatus=0, kITSClusterMap, kTPCNCrossedRows, kTPCNclsF, kTPCnclsS,  kTPCncls, kChi2perNDF,kX, kY, kZ, kTPCsignalN, kdEdx, kbeta, kLabel, kDCAxy, kDCAz,
        kNsigmaTPCe, kNsigmaTPCpi, kNsigmaTPCk, kNsigmaTPCp, kNsigmaTOFe, kNsigmaTOFpi, kNsigmaTOFk, kNsigmaTOFp, kMappingTracks};
  enum {kTklptMC=0, kTkletaMC, kTklphiMC, kTklpdg, kMappingTracklets};
  enum {kMuonDCA=0, kMuonChi2perNDF, kMuonRabs, kMuonpDCA, kMuonMCPt, kMuonMCEta, kMuonMCPhi, kMuonMCPdg, kMuonMCprimPt, kMuonMCprimEta, kMuonMCprimPhi, kMuonMCprimPdg, kMuonMCoriginpt, kMuonMCoriginEta, kMuonMCoriginPhi, kMuonMCoriginPdg, kMappingMuons};
  enum {kMCindex=0, kMCLabel, kMCIsPhysPrim, kMCmotherpdg, kMappingMCTracks};

  Int_t* MappingTracks() { return fMappingTracks; }
  Int_t* MappingTracklets() { return fMappingTracklets; }
  Int_t* MappingMuons() { return fMappingMuons; }
  Int_t* MappingMCTracks() { return fMappingMCTracks; }

 protected:
  Int_t fMappingTracks[kMappingTracks];		//  mapping array of positions of additional parameters of the AliCFParticle tracks
  Int_t fMappingTracklets[kMappingTracklets];	//  mapping array of positions of additional parameters of the AliCFParticle tracklets
  Int_t fMappingMuons [kMappingMuons];		//  mapping array of positions of additional parameters of the AliCFParticle muons
  Int_t fMappingMCTracks[kMappingMCTracks];	//  mapping array of positions of additional parameters of the AliCFParticle MC generated tracks

  ClassDef(AliCFTreeMapping,2);
};
#endif
