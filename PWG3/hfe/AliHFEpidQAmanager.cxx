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
// Class AliHFEpidQAmanager
// Does steering of the PID QA. The PID QA manager is initialized according
// to the configuration used for PID. It contains PID objects inheriting 
// from AliHFEdetPIDqa, which are
//     AliHFEtpcPIDqa
//     AliHFEtrdPIDqaV1
//     AliHFEtofPIDqa
// PID QA objects are filled for every detector before PID decision and 
// after PID decision for tracks which are selected by the given detector
// as electron candidates.
//
// Author
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TList.h>

#include "AliAODpidUtil.h"
#include "AliESDpid.h"
#include "AliVParticle.h"

#include "AliHFEtpcPIDqa.h"
#include "AliHFEtrdPIDqaV1.h"
#include "AliHFEtofPIDqa.h"
#include "AliHFEemcalPIDqa.h"  //s.s
#include "AliHFEbayesPIDqa.h"
#include "AliHFEpidQAmanager.h"

ClassImp(AliHFEpidQAmanager)

//____________________________________________________________
AliHFEpidQAmanager::AliHFEpidQAmanager():
  TObject()
{
  //
  // Dummy constructor
  //
  memset(fDetPIDqa, 0, sizeof(AliHFEdetPIDqa *) * AliHFEpid::kNdetectorPID);
  memset(fDetPID, 0, sizeof(AliHFEdetPIDqa *) * AliHFEpid::kNdetectorPID);
  SetOwner(); 
}

//____________________________________________________________
AliHFEpidQAmanager::AliHFEpidQAmanager(const AliHFEpidQAmanager &ref):
  TObject(ref)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
  SetOwner();
}

//____________________________________________________________
AliHFEpidQAmanager &AliHFEpidQAmanager::operator=(const AliHFEpidQAmanager &ref){
  //
  // Assignment operator
  //
  if(this != &ref) ref.Copy(*this);
  SetOwner();
  return *this;
}

//____________________________________________________________
void AliHFEpidQAmanager::Copy(TObject &o) const{
  //
  // Make copy
  //
  TObject::Copy(o);
  AliHFEpidQAmanager &target = dynamic_cast<AliHFEpidQAmanager &>(o);

  for(Int_t idet = 0; idet < AliHFEpid::kNdetectorPID; idet++){
    target.fDetPID[idet] = fDetPID[idet]; 
    if(target.fDetPIDqa[idet]) delete target.fDetPIDqa[idet];
    if(fDetPIDqa[idet]) target.CreateDetPIDqa(static_cast<AliHFEpid::EDETtype_t>(idet));
  }
}

//____________________________________________________________
AliHFEpidQAmanager::~AliHFEpidQAmanager(){
  //
  // Destructor
  //
  if(IsOwner())
    for(Int_t idet = 0; idet < AliHFEpid::kNdetectorPID; idet++){
      if(fDetPIDqa[idet]) delete fDetPIDqa[idet];
    }
}

//____________________________________________________________
void AliHFEpidQAmanager::Initialize(AliHFEpid *pid){
  //
  // Initialize PID QA manager according to the detector
  // configuration used in the PID
  //
  for(Int_t idet = 0; idet < AliHFEpid::kNdetectorPID; idet++){
    // Initialize Array for detector PID for all detectors
    fDetPID[idet] = pid->GetDetPID(static_cast<AliHFEpid::EDETtype_t>(idet));
    if(pid->HasDetector(static_cast<AliHFEpid::EDETtype_t>(idet))){
      CreateDetPIDqa(static_cast<AliHFEpid::EDETtype_t>(idet));
      if(fDetPIDqa[idet]){
        fDetPIDqa[idet]->SetPIDqaManager(this);
        fDetPIDqa[idet]->Initialize();
      }
    }
  }
}

//____________________________________________________________
void AliHFEpidQAmanager::CreateDetPIDqa(AliHFEpid::EDETtype_t idet){
  //
  // Create new PID QA object
  //
  switch(idet){
    case AliHFEpid::kTPCpid: 
          fDetPIDqa[idet] = new AliHFEtpcPIDqa("TPCQA"); 
          break;
    case AliHFEpid::kTRDpid: 
          fDetPIDqa[idet] = new AliHFEtrdPIDqaV1("TRDQA"); 
          break;
    case AliHFEpid::kTOFpid: 
          fDetPIDqa[idet] = new AliHFEtofPIDqa("TOFQA"); 
          break;
    case AliHFEpid::kEMCALpid: //s.s (name) 
          fDetPIDqa[idet] = new AliHFEemcalPIDqa("EMCALQA"); // s.s 
	  break;  //s.s
    case AliHFEpid::kBAYESpid:
          fDetPIDqa[idet] = new AliHFEbayesPIDqa("BAYESQA");
          break;

    default:
          break;
  };
}

//____________________________________________________________
void AliHFEpidQAmanager::ProcessTrack(const AliHFEpidObject *track, AliHFEpid::EDETtype_t det, AliHFEdetPIDqa::EStep_t step){
  //
  // Process single Track
  //
  if(!fDetPIDqa[det]){
    AliDebug(1, Form("QA for detector %d not available", det));
    return;
  }
  AliDebug(1, Form("Doing QA for detector %d\n", det));
  fDetPIDqa[det]->ProcessTrack(track, step);
}

//____________________________________________________________
TList *AliHFEpidQAmanager::MakeList(const Char_t *name){
  //
  // Make List of PID QA objects
  //
  TList *list = new TList;
  list->SetName(name);
  list->SetOwner();
  for(Int_t idet = 0; idet < AliHFEpid::kNdetectorPID; idet++){
    if(fDetPIDqa[idet]) list->Add(fDetPIDqa[idet]);
  }
  ReleaseOwnerShip();
  return list;
}

