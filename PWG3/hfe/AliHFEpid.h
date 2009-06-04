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
#ifndef __ALIHFEPID_H__
#define __ALIHFEPID_H__

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class AliHFEpidBase;
class AliESDtrack;
class AliVParticle;
class AliMCEvent;

class TList;

class AliHFEpid : public TObject{
  enum{
    kIsQAOn = BIT(14),
    kHasMCData = BIT(15)
  };
  enum{
    kMCpid = 0,
    kESDpid = 1,
    kTPCpid = 2,
    kTRDpid = 3,
    kTOFpid = 4,
    kNdetectorPID = 5
  };
  public:
    AliHFEpid();
    AliHFEpid(const AliHFEpid &c);
    AliHFEpid &operator=(const AliHFEpid &c);
    ~AliHFEpid();
    
    Bool_t InitializePID(TString detectors);
    Bool_t IsSelected(AliVParticle *track);
    void SetMCEvent(AliMCEvent *mc);

    Bool_t IsQAOn() const { return TestBit(kIsQAOn); };
    Bool_t HasMCData() const { return TestBit(kHasMCData); };
    void SetQAOn();
    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData, hasMCdata); };
    TList *GetQAhistograms() const { return fQAlist; };

  protected:
    Bool_t MakePID_TPC_TOF(AliESDtrack *track);
  private:
    AliHFEpidBase *fDetectorPID[kNdetectorPID];    //! Detector PID classes
    UInt_t fEnabledDetectors;             // Enabled Detectors
    TList *fQAlist;                       //! QA histograms

  ClassDef(AliHFEpid, 1)      // Steering class for Electron ID
};

#endif
