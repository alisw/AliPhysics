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
//                Dielectron Variables Container Class                   //
//                                                                       //
/*

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliDielectronVarContainer.h"

ClassImp(AliDielectronVarContainer)

const char* AliDielectronVarContainer::fgkParticleNames[AliDielectronVarContainer::kNMaxValues] = {
  // Leg specific variables
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
  // Pair specific variables
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
  // Event specific variables
  "X",
  "Y",
  "Z",
  "XRes",
  "YRes",
  "ZRes",
  "NContributors",
  "BzkG",
  "NTrk",
  "Nacc",
  "kNaccTrcklts",
  "Centrality",
  "Nevents"
};


const char* AliDielectronVarContainer::fgkParticleNamesMC[AliDielectronVarContainer::kNMaxValues_MC] = {
  // Leg specific variables
  "Px_MC",
  "Py_MC",
  "Pz_MC",
  "Pt_MC",
  "P_MC",
  "Xv_MC",
  "Yv_MC",
  "Zv_MC",
  "OneOverPt_MC",
  "Phi_MC",
  "Theta_MC",
  "Eta_MC",
  "Y_MC",
  "E_MC",
  "M_MC",
  "Charge_MC",
  "ImpactParXY_MC",
  "ImpactParZ_MC",
  "PdgCode",
  "PdgCodeMother",
  "PdgCodeGrandMother",
  "NumberOfDaughters",
  "HaveSameMother",
  "IsJpsiPrimary",
  // Pair specific variables
  "DecayLength_MC",
  "R_MC",
  "OpeningAngle_MC",
  "ThetaHE_MC",
  "PhiHE_MC",
  "ThetaCS_MC",
  "PhiCS_MC",
  "LegDistance_MC",
  "LegDistanceXY_MC",
  "DeltaEta_MC",
  "DeltaPhi_MC",
  "DCA_MC",
  "PairType_MC",
  "PseudoProperTime_MC",
  // Event specific variables
  "X_MC",
  "Y_MC",
  "Z_MC",
  "kNch",
  "Centrality_MC",
  "Nevents"
};


AliPIDResponse* AliDielectronVarContainer::fgPIDResponse = 0x0;
AliKFVertex*    AliDielectronVarContainer::fgKFVertex    = 0x0;
AliAODVertex*   AliDielectronVarContainer::fgAODVertex   = 0x0;
Double_t        AliDielectronVarContainer::fgData[AliDielectronVarContainer::kNMaxValues] = {};
Double_t        AliDielectronVarContainer::fgDataMC[AliDielectronVarContainer::kNMaxValues_MC] = {};


//________________________________________________________________
AliDielectronVarContainer::AliDielectronVarContainer() :
  TNamed("AliDielectronVarContainer","AliDielectronVarContainer")
{
  //
  // Default constructor
  //
}


//________________________________________________________________
AliDielectronVarContainer::AliDielectronVarContainer(const char* name, const char* title) :
  TNamed(name,title)
{
  //
  // Named constructor
  //
}


//________________________________________________________________
AliDielectronVarContainer::~AliDielectronVarContainer()
{
  //
  // Default destructor
  //

  if (fgKFVertex) delete fgKFVertex;
  fgKFVertex = 0x0;
  if (fgAODVertex) delete fgAODVertex;
  fgAODVertex = 0x0;
  if (fgPIDResponse) delete fgPIDResponse;
  fgPIDResponse = 0x0;
}

