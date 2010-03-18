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
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"

enum FillType_t {real, mixed};
enum PoolType_t {triggs, assocs};
enum BinType_t  {centr, zvert, trigg, assoc}; 
enum PartType_t {unknown, hadron, proton, kaon, pion, photon, electron, jet, 
		 dihadron, diphoton, dielectron, dijet};

const UInt_t kMAXNUMCORREL = 1;     // Maximum no of correlations
const UInt_t kMAXVERTBIN   = 10;    // Maximum no of vertex bins
const UInt_t kMAXCENTBIN   = 3;     // Maximum no of centrality bins
const UInt_t kMAXTRIGBIN   = 10;    // Maximum no of trigger bins
const UInt_t kMAXASSOBIN   = 10;    // Maximum no of associated bins

#endif
