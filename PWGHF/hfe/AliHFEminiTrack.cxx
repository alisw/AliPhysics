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
// Originator:  M.Fasel <M.Fasel@gsi.de>
// Base class for AliHFEminiEventCreator.cxx/h
// Contributors: Nirbhay K. Behera, Jiyeon Kwon, Jonghan Park

#include <cstring>

#include "AliHFEminiTrack.h"

ClassImp(AliHFEminiTrack)


//_______________________________________
AliHFEminiTrack::AliHFEminiTrack():
  TObject(),
  fSignedPt(0.),
  fMCSignedPt(0.),
  fMCPDG(0),
  fMCMotherPdg(0),
  fMCGrandMotherPdg(0),
  fMCSource(5),
  fMCEleSource(0),
  fMCEleSourcePt(0)
{
  // 
  // Default Constuctor
  //
  fHFEImpactParam[0] = -999.;
  fHFEImpactParam[1] = -999.;
}

//_______________________________________
AliHFEminiTrack::AliHFEminiTrack(const AliHFEminiTrack &ref):
  TObject(ref),
  fSignedPt(ref.fSignedPt),
  fMCSignedPt(ref.fMCSignedPt),
  fMCPDG(ref.fMCPDG),
  fMCMotherPdg(ref.fMCMotherPdg),
  fMCGrandMotherPdg(ref.fMCGrandMotherPdg),
  fMCSource(ref.fMCSource),
  fMCEleSource(ref.fMCEleSource),
  fMCEleSourcePt(ref.fMCEleSourcePt)
{
  // 
  // Copy Constuctor
  //
  memcpy(fHFEImpactParam, ref.fHFEImpactParam, sizeof(Double_t) * 2);
}

//_______________________________________
AliHFEminiTrack &AliHFEminiTrack::operator=(const AliHFEminiTrack &ref){
  //
  // Assignment Operator
  //
  if(&ref != this){
    TObject::operator=(ref);
    fSignedPt = ref.fSignedPt;
    fMCSignedPt = ref.fMCSignedPt;
    fMCPDG = ref.fMCPDG;
    fMCMotherPdg = ref.fMCMotherPdg;
    fMCGrandMotherPdg = ref.fMCGrandMotherPdg;
    fMCSource = ref.fMCSource;
    fMCEleSource = ref.fMCEleSource;
    fMCEleSourcePt = ref.fMCEleSourcePt;
    memcpy(fHFEImpactParam, ref.fHFEImpactParam, sizeof(Double_t) * 2);
  }
  return *this;
}
