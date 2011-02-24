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

/* $Id$ */

//
// Abstract PID base class for Detector PID classes 
// Supplies detector PID classes with basic informations (i.e. Debug 
// Level)
//  
// Authors: 
//   Markus Fasel <M.Fasel@gsi.de> 
// 

#include "AliAODpidUtil.h"
#include "AliESDpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEtools.h"

ClassImp(AliHFEpidBase)

//___________________________________________________________________
AliHFEpidBase::AliHFEpidBase():
  TNamed(),
  fESDpid(NULL),
  fAODpid(NULL)
{
  //
  // Default constructor
  //
}

//___________________________________________________________________
AliHFEpidBase::AliHFEpidBase(const Char_t *name):
  TNamed(name, ""),
  fESDpid(NULL),
  fAODpid(NULL)
{
  //
  // Default constructor
  //
}

//___________________________________________________________________
AliHFEpidBase::AliHFEpidBase(const AliHFEpidBase &c):
  TNamed(),
  fESDpid(NULL),
  fAODpid(NULL)
{
  //
  //Copy constructor
  //
  c.Copy(*this);
}

//___________________________________________________________________
AliHFEpidBase &AliHFEpidBase::operator=(const AliHFEpidBase &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  }

  return *this;
}

//___________________________________________________________________
void AliHFEpidBase::Copy(TObject &ref) const {
  //
  // Copy function for assignment operator
  //
  AliHFEpidBase &target = dynamic_cast<AliHFEpidBase &>(ref);

  target.fESDpid = fESDpid;
  target.fAODpid = fAODpid;

  TNamed::Copy(ref);
}

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

