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
  "NFclsTPC",
  "TPCsignalN",
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

  "NumberOfDaughters",
  "HaveSameMother",
  "ITS_signal",
  "SSD1_signal",
  "SSD2_signal",
  "SDD1_signal",
  "SDD2_signal",

  "P_InnerParam",
  "TPC_signal",
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
  "Merr",
  "DCA",
  "PairType",
  //
  "X",
  "Y",
  "Z",
  "XRes",
  "YRes",
  "ZRes",
  "NTrk",
  "Tracks",
  "Nevents"
};

AliESDpid* AliDielectronVarManager::fgESDpid = 0x0;
AliVEvent* AliDielectronVarManager::fgEvent  = 0x0;
AliKFVertex* AliDielectronVarManager::fgKFVertex  = 0x0;
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

