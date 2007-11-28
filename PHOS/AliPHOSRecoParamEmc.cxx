/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

/* $Id$ */

// This class contains the PHOS EMC reconstruction parameters.
// To get default parameters from any place of your code, use:
//   AliPHOSRecoParam* defPar = AliPHOSRecoParamEmc::GetEmcDefaultParameters();
//   Float_t cluth = defPar->GetClusteringThreshold();
//   ...

// --- ROOT header files ---
#include "TObjArray.h"
#include "TString.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

// --- AliRoot header files ---
#include "AliPHOSRecoParamEmc.h"

ClassImp(AliPHOSRecoParamEmc)

TObjArray* AliPHOSRecoParamEmc::fgkMaps =0; //ALTRO mappings

//-----------------------------------------------------------------------------
AliPHOSRecoParamEmc::AliPHOSRecoParamEmc() : AliPHOSRecoParam()
{
  //Default constructor.

  SetClusteringThreshold(0.2);
  SetLocalMaxCut(0.03);
  SetMinE(0.01);
  SetLogWeight(4.5);
}

//-----------------------------------------------------------------------------
AliPHOSRecoParam* AliPHOSRecoParamEmc::GetEmcDefaultParameters()
{
  //Default parameters for the reconstruction in EMC.

  AliPHOSRecoParam* params = new AliPHOSRecoParamEmc();
  return params;
}

//-----------------------------------------------------------------------------
const TObjArray* AliPHOSRecoParamEmc::GetMappings()
{
  //Returns array of AliAltroMappings for RCU0..RCU3.
  //If not found, read it from OCDB.

  //Quick check as follows:
  //  root [0] AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  //  root [1] AliCDBManager::Instance()->SetRun(1);
  //  root [2] TObjArray* maps = AliPHOSRecoParamEmc::GetMappings();
  //  root [3] maps->Print();
  
  if(fgkMaps) return fgkMaps;
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("PHOS/Calib/Mapping");
  if(entry)
    fgkMaps = (TObjArray*)entry->GetObject();
  
  return fgkMaps;
  
}
