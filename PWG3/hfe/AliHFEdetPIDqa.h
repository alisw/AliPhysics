#ifndef ALIHFEDETPIDQA_H
#define ALIHFEDETPIDQA_H

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
// Class AliHFEdetPIDqa
// Base class for detector PID QA describing the interface to the PID QA
// manager, keeping also commom functionality
// More information can be found inside the implementation file
//
#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class AliHFEpidObject;
class AliHFEpidQAmanager;

class AliHFEdetPIDqa : public TNamed{
  public:
    enum EStep_t{
      kBeforePID = 0,
      kAfterPID = 1
    };
    AliHFEdetPIDqa();
    AliHFEdetPIDqa(const Char_t *name, const Char_t *title);
    AliHFEdetPIDqa(const AliHFEdetPIDqa &c);
    AliHFEdetPIDqa &operator=(const AliHFEdetPIDqa &o);
    ~AliHFEdetPIDqa(){}

    virtual void Initialize() = 0;
    virtual void ProcessTrack(const AliHFEpidObject *track, EStep_t step)= 0;

    void SetPIDqaManager(AliHFEpidQAmanager *qamgr) { fQAmanager = qamgr; }
    AliHFEpidQAmanager *GetPIDqaManager() const { return fQAmanager; }

  protected:
    AliHFEpidQAmanager *fQAmanager; // PID QA manager
  
    ClassDef(AliHFEdetPIDqa, 1)     // Base class for detector PID QA
};

#endif
