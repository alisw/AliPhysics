#ifndef __CORRELDEFS_H__
#define __CORRELDEFS_H__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

// CorrelDefs.h is at the top of preprocessor includes,
// hence we add the general headers here
//-- Author: Paul Constantin

// C++ headers:
#include <iostream>

// ROOT headers:
#include <TROOT.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TNtuple.h>
#include <TString.h>
#include <TRandom2.h>

// AliRoot headers:
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"

namespace JetCorrelHD {

  enum FillType_t {real, mixed};
  enum PoolType_t {triggs, assocs};
  enum BinType_t  {centr, zvert}; 
  enum PartType_t {unknown, hadron, proton, kaon, pion, photon, electron, jet, 
		   dihadron, diphoton, dielectron, dijet};

  const Float_t kEPS = 1.e-6;
  const Float_t pi = TMath::Pi();
  const Bool_t kUseAliKF = kFALSE; // reconstruction with AliKFParticle or TLorentzVector
  const Float_t kZ0MassMean = 91.;
  const Float_t kZ0MassSig  = 3.;
  const Float_t kPi0MassMean = 0.140;  // Pi0 mass range for diphoton selection
  const Float_t kPi0MassSig  = 0.010;

  const UInt_t kDPhiNumBins = 62;       // number of bins in DeltaPhi histos (2pi-bin=0.1)
  const UInt_t kDEtaNumBins = 40;       // number of bins in DeltaEta histos (3.6-bin=0.1)
  
  const UInt_t kMAXNUMCORREL = 3;     // Maximum no of correlations
  const UInt_t kMAXVERTBIN   = 10;    // Maximum no of vertex bins
  const UInt_t kMAXCENTBIN   = 3;     // Maximum no of centrality bins

} // namespace declaration

#endif
