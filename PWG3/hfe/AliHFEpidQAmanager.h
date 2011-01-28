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
// Steering PID QA
// More information can be found inside the implementation file
//
#ifndef ALIHFEPIDQAMANAGER_H
#define ALIHFEPIDQAMANAGER_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ALIHFEPID_H
#include "AliHFEpid.h"
#endif

#ifndef ALIHFEDETPIDQA_H
#include "AliHFEdetPIDqa.h"
#endif

class TList; 
class AliVParticle;
class AliESDpid;
class AliAODpidUtil;
class AliHFEpidObject;

class AliHFEpidQAmanager : public TObject{
  public:
    AliHFEpidQAmanager();
    AliHFEpidQAmanager(const AliHFEpidQAmanager &ref);
    AliHFEpidQAmanager &operator=(const AliHFEpidQAmanager &ref);
    void Copy(TObject &o) const;
    ~AliHFEpidQAmanager();

    void Initialize(AliHFEpid *pid);
    void ProcessTrack(const AliHFEpidObject *track, AliHFEpid::EDETtype_t det, AliHFEdetPIDqa::EStep_t step);
    AliHFEdetPIDqa *GetDetectorPIDqa(AliHFEpid::EDETtype_t detector) const { return fDetPIDqa[detector]; }
    AliHFEpidBase *GetDetectorPID(AliHFEpid::EDETtype_t detector) const { return fDetPID[detector]; }
    TList *MakeList(const Char_t *name);
    void SetHighResolutionHistos() { SetBit(kHighResolutionHistos, kTRUE); };
    Bool_t HasHighResolutionHistos() const { return TestBit(kHighResolutionHistos); }

  protected:
    enum{
      kIsOwner = BIT(14),
      kHighResolutionHistos = BIT(15)
    };
    Bool_t IsOwner() const { return TestBit(kIsOwner); }
    void SetOwner() { SetBit(kIsOwner, kTRUE); }
    void ReleaseOwnerShip() { SetBit(kIsOwner, kFALSE); }
    void CreateDetPIDqa(AliHFEpid::EDETtype_t detector);

  private:
    AliHFEdetPIDqa *fDetPIDqa[AliHFEpid::kNdetectorPID]; //! Detector PID QA objects
    AliHFEpidBase *fDetPID[AliHFEpid::kNdetectorPID];    //  Detector PID objects

  ClassDef(AliHFEpidQAmanager, 0)
};
#endif
