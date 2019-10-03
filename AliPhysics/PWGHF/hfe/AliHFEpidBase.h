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
// Base Class for Detector PID Objects
// For more information see the implementation file
//
#ifndef ALIHFEPIDBASE_H
#define ALIHFEPIDBASE_H
 
#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ALIHFEPIDOBJECT_H
#include "AliHFEpidObject.h"
#endif

class TList;
class AliPIDResponse;
class AliVParticle;
class AliMCParticle;
class AliHFEpidQAmanager;

class AliHFEpidBase : public TNamed{
  public:
    AliHFEpidBase();
    AliHFEpidBase(const Char_t *name);
    AliHFEpidBase(const AliHFEpidBase &c);
    AliHFEpidBase &operator=(const AliHFEpidBase &c);
    virtual ~AliHFEpidBase() {};
    // Framework functions that have to be implemented by the detector PID classes
    virtual Bool_t InitializePID(Int_t run) = 0;
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa = NULL) const = 0;

    Bool_t HasMCData() const { return TestBit(kHasMCData); };

    void SetPIDResponse(const AliPIDResponse * const pid) { fkPIDResponse = pid; }
    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };

    const AliPIDResponse *GetPIDResponse() const { return fkPIDResponse; }; 

  protected:
    const AliPIDResponse *fkPIDResponse;        //! PID Response
    void Copy(TObject &ref) const;

  private:
    enum{
      kHasMCData = BIT(14)
    };

    ClassDef(AliHFEpidBase, 2)      // Base class for detector Electron ID
};
#endif
