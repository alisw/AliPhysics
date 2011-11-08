/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//  
// Object used in the electron identification
// Combines reconstructed tracks with additional information (like apriori PID, centrality)
// which are missing in the reconstructed track itself
//
// Authors: 
//   Markus Fasel <M.Fasel@gsi.de> 
// 
#include "AliHFEpidObject.h"
#include "AliHFEtools.h"

//___________________________________________________________________
AliHFEpidObject &AliHFEpidObject::operator=(const AliHFEpidObject &ref){
  //
  // Assignment operator
  //
  if(&ref != this){
    fkRecTrack = ref.fkRecTrack;
    fAnalysisType = ref.fAnalysisType;
    fAbInitioPID = ref.fAbInitioPID;
    fCentrality = ref.fCentrality;
  }
  return *this;
}

//___________________________________________________________________
void AliHFEpidObject::SetMCTrack(const AliVParticle *mctrack){
  //
  // Set the aprioriPID information coming from the MC truth
  //
  if(mctrack) fAbInitioPID = AliHFEtools::PDG2AliPID(AliHFEtools::GetPdg(mctrack));
}

