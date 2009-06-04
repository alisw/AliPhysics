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
/************************************************************************
 *                                                                      *
 * PID Steering Class                                                   *
 * Interface to the user task                                           *
 *                                                                      *
 * Authors:                                                             *
 *   Markus Fasel <M.Fasel@gsi.de>                                      *
 *                                                                      *
 ************************************************************************/
#include <TClass.h>
#include <TIterator.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include "AliESDtrack.h"

#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidTRD.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidMC.h"

ClassImp(AliHFEpid)

//____________________________________________________________
AliHFEpid::AliHFEpid():
  fEnabledDetectors(0),
  fQAlist(0x0)
{
  //
  // Default constructor
  //
  memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
}

//____________________________________________________________
AliHFEpid::AliHFEpid(const AliHFEpid &c):
  TObject(c),
  fEnabledDetectors(c.fEnabledDetectors),
  fQAlist(0x0)
{
  //
  // Copy Constructor
  //
  memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
  if(c.fDetectorPID[kMCpid])
    fDetectorPID[kMCpid] = new AliHFEpidMC(*(dynamic_cast<AliHFEpidMC *>(c.fDetectorPID[kMCpid])));
  if(c.fDetectorPID[kTPCpid])
    fDetectorPID[kTPCpid] = new AliHFEpidTPC(*(dynamic_cast<AliHFEpidTPC *>(c.fDetectorPID[kTPCpid])));
  if(c.fDetectorPID[kTRDpid])
    fDetectorPID[kTRDpid] = new AliHFEpidTRD(*(dynamic_cast<AliHFEpidTRD *>(c.fDetectorPID[kTOFpid])));
  if(c.fDetectorPID[kTOFpid])
    fDetectorPID[kTOFpid] = new AliHFEpidTOF(*(dynamic_cast<AliHFEpidTOF *>(c.fDetectorPID[kTOFpid])));
  if(c.IsQAOn()) SetQAOn();
  if(c.HasMCData()) SetHasMCData(kTRUE);
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(c.IsQAOn() && fDetectorPID[idet]) fDetectorPID[idet]->SetQAOn(fQAlist);
    if(c.HasMCData() && fDetectorPID[idet]) fDetectorPID[idet]->SetHasMCData(kTRUE);
  }
}

//____________________________________________________________
AliHFEpid& AliHFEpid::operator=(const AliHFEpid &c){
  //
  // Assignment operator
  //
  TObject::operator=(c);

  if(this != &c){
    fEnabledDetectors = c.fEnabledDetectors;
    fQAlist = 0x0;
  
    memset(fDetectorPID, 0, sizeof(AliHFEpidBase *) * kNdetectorPID);
    if(c.fDetectorPID[kMCpid])
      fDetectorPID[kMCpid] = new AliHFEpidMC(*(dynamic_cast<AliHFEpidMC *>(c.fDetectorPID[kMCpid])));
    if(c.fDetectorPID[kTPCpid])
      fDetectorPID[kTPCpid] = new AliHFEpidTPC(*(dynamic_cast<AliHFEpidTPC *>(c.fDetectorPID[kTPCpid])));
    if(c.fDetectorPID[kTRDpid])
      fDetectorPID[kTRDpid] = new AliHFEpidTRD(*(dynamic_cast<AliHFEpidTRD *>(c.fDetectorPID[kTOFpid])));
    if(c.fDetectorPID[kTOFpid])
      fDetectorPID[kTOFpid] = new AliHFEpidTOF(*(dynamic_cast<AliHFEpidTOF *>(c.fDetectorPID[kTOFpid])));
    if(c.IsQAOn()) SetQAOn();
    if(c.HasMCData()) SetHasMCData(kTRUE);
    for(Int_t idet = 0; idet < kNdetectorPID; idet++){
      if(c.IsQAOn() && fDetectorPID[idet]) fDetectorPID[idet]->SetQAOn(fQAlist);
      if(c.HasMCData() && fDetectorPID[idet]) fDetectorPID[idet]->SetHasMCData();
    }
  }
  return *this; 
}

//____________________________________________________________
AliHFEpid::~AliHFEpid(){
  //
  // Destructor
  //
  if(fQAlist) delete fQAlist; fQAlist = 0x0;  // Each detector has to care about its Histograms
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(fDetectorPID[idet]) delete fDetectorPID[idet];
  }
}

//____________________________________________________________
Bool_t AliHFEpid::InitializePID(TString detectors){
  //
  // Initializes PID Object:
  // + Defines which detectors to use
  // + Initializes Detector PID objects
  // + Handles QA
  //
  fDetectorPID[kMCpid] = new AliHFEpidMC("Monte Carlo PID"); // Always there
  SETBIT(fEnabledDetectors, kMCpid);
  
  TObjArray *detsEnabled = detectors.Tokenize(":");
  TIterator *detIterator = detsEnabled->MakeIterator();
  TObjString *det = 0x0;
  while((det = dynamic_cast<TObjString *>(detIterator->Next()))){
    if(det->String().CompareTo("TPC") == 0){
      fDetectorPID[kTPCpid] = new AliHFEpidTPC("TPC PID");
      SETBIT(fEnabledDetectors, kTPCpid);
    } else if(det->String().CompareTo("TRD") == 0){
      fDetectorPID[kTRDpid] = new AliHFEpidTRD("TRD PID");
      SETBIT(fEnabledDetectors, kTRDpid);
    } else if(det->String().CompareTo("TOF") == 0){
      fDetectorPID[kTOFpid] = new AliHFEpidTOF("TOF PID");
      SETBIT(fEnabledDetectors, kTOFpid);
    }
    // Here is still space for ESD PID
  }
  // Initialize PID Objects
  Bool_t status = kTRUE;
  for(Int_t idet = 0; idet < kNdetectorPID; idet++){
    if(fDetectorPID[idet]){ 
      status &= fDetectorPID[idet]->InitializePID();
      if(IsQAOn() && status) fDetectorPID[idet]->SetQAOn(fQAlist);
      if(HasMCData() && status) fDetectorPID[idet]->SetHasMCData();
    }
  }
  return status;
}

//____________________________________________________________
Bool_t AliHFEpid::IsSelected(AliVParticle *track){
  //
  // Steers PID decision for single detectors respectively combined
  // PID decision
  //
  
  if(TString(track->IsA()->GetName()).CompareTo("AliMCparticle") == 0){
    return (TMath::Abs(fDetectorPID[kMCpid]->IsSelected(track)) == 11);
  }
  if(TString(track->IsA()->GetName()).CompareTo("AliESDtrack") == 0){
    if(TESTBIT(fEnabledDetectors, kTPCpid) && TESTBIT(fEnabledDetectors, kTOFpid)){
      // case TPC-TOF
      return MakePID_TPC_TOF(dynamic_cast<AliESDtrack *>(track));
    } else if(TESTBIT(fEnabledDetectors, kTPCpid)){
      return (TMath::Abs(fDetectorPID[kTPCpid]->IsSelected(track)) ==11);
    } else if(TESTBIT(fEnabledDetectors, kTRDpid)){
      return (TMath::Abs(fDetectorPID[kTRDpid]->IsSelected(track)) ==11);
    } else if(TESTBIT(fEnabledDetectors, kTOFpid)){
      return (TMath::Abs(fDetectorPID[kTOFpid]->IsSelected(track)) ==11);
    }
    
  }
  return kFALSE;
}

//____________________________________________________________
Bool_t AliHFEpid::MakePID_TPC_TOF(AliESDtrack *track){
  //
  // Combines TPC and TOF PID decision
  //
  if(fDetectorPID[kTOFpid]->IsSelected(track)) return fDetectorPID[kTPCpid]->IsSelected(track);
  return kFALSE;
}

//____________________________________________________________
void AliHFEpid::SetQAOn(){
  //
  // Switch on QA
  //
  SetBit(kIsQAOn);
  fQAlist = new TList;
  fQAlist->SetName("PIDqa");
}

void AliHFEpid::SetMCEvent(AliMCEvent *event){
  for(Int_t idet = 0; idet < kNdetectorPID; idet++)
    if(fDetectorPID[idet]) fDetectorPID[idet]->SetMCEvent(event);
}
