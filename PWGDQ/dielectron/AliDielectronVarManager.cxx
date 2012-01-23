/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                Dielectron Variables Manager class                     //
//                                                                       //
/*

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliDielectronVarManager.h"

ClassImp(AliDielectronVarManager)

const char* AliDielectronVarManager::fgkParticleNames[AliDielectronVarManager::kNMaxValues] = {
  "Px",
  "Py",
  "Pz",
  "Pt",
  "P",
  "Xv",
  "Yv",
  "Zv",
  "OneOverPt",
  "Phi",
  "Theta",
  "Eta",
  "Y",
  "E",
  "M",
  "Charge",
  "NclsITS",
  "NclsTPC",
  "NclsTPCiter1",
  "NFclsTPC",
  "NFclsTPCrobust",
  "NFclsTPCrobustFraction",
  "TPCsignalN",
  "TPCsignalNfrac",
  "TPCchi2PerCluster",
  "TrackStatus",
    
  "NclsTRD",
  "TRDntracklets",
  "TRDpidQuality",
  "TRDpidProb_Electrons",
  "TRDpidProb_Pions",
  "ImpactParXY",
  "ImpactParZ",
  "TrackLength",
  "PdgCode",

  "PdgCodeMother",
  "PdgCodeGrandMother",
    
  "NumberOfDaughters",
  "HaveSameMother",
  "IsJpsiPrimary",
  "ITS_signal",
  "SSD1_signal",
  "SSD2_signal",
  "SDD1_signal",
  "SDD2_signal",
  "ITS_clusterMap",
  "ITS_nSigma_Electrons",
  "ITS_nSigma_Pions",
  "ITS_nSigma_Muons",
  "ITS_nSigma_Kaons",
  "ITS_nSigma_Protons",

  "P_InnerParam",
  "TPC_signal",
  "TOF_signal",
  "TOF_beta",
  
  "TPC_nSigma_Electrons",
  "TPC_nSigma_Pions",
  "TPC_nSigma_Muons",
  "TPC_nSigma_Kaons",
  "TPC_nSigma_Protons",

  "TOF_nSigma_Electrons",
  "TOF_nSigma_Pions",
  "TOF_nSigma_Muons",
  "TOF_nSigma_Kaons",
  "TOF_nSigma_Protons",

  "KinkIndex0",
  //
  "Chi2NDF",
  "DecayLength",
  "R",
  "OpeningAngle",
  "ThetaHE",
  "PhiHE",
  "ThetaCS",
  "PhiCS",
  "LegDistance",
  "LegDistanceXY",
  "DeltaEta",
  "DeltaPhi",
  "Merr",
  "DCA",
  "PairType",
  "PseudoProperTime",
  //
  "X",
  "Y",
  "Z",
  "XRes",
  "YRes",
  "ZRes",
  "NTrk",
  "Tracks",
  "Nacc",
  "kNaccTrcklts",
  "kNch",
  "Centrality",
  "Nevents"
};

AliPIDResponse* AliDielectronVarManager::fgPIDResponse = 0x0;
AliVEvent*      AliDielectronVarManager::fgEvent       = 0x0;
AliKFVertex*    AliDielectronVarManager::fgKFVertex    = 0x0;
Double_t        AliDielectronVarManager::fgData[AliDielectronVarManager::kNMaxValues] = {};
//________________________________________________________________
AliDielectronVarManager::AliDielectronVarManager() :
  TNamed("AliDielectronVarManager","AliDielectronVarManager")
{
  //
  // Default constructor
  //

}

//________________________________________________________________
AliDielectronVarManager::AliDielectronVarManager(const char* name, const char* title) :
  TNamed(name,title)
{
  //
  // Named constructor
  //
  
}

//________________________________________________________________
AliDielectronVarManager::~AliDielectronVarManager()
{
  //
  // Default destructor
  //
}
