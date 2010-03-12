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
  "ImpactParXY",
  "ImpactParZ",
  "TrackLength",
  "PdgCode",
  "P_InnerParam",
  "TPC_signal",
  "TPC_nSigma_Electrons",
  //
  "Chi2NDF",
  "DecayLength",
  "R",
  "OpeningAngle",
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

AliESDpid* AliDielectronVarManager::fgESDpid = new AliESDpid;
AliVEvent* AliDielectronVarManager::fgEvent  = 0x0;
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

