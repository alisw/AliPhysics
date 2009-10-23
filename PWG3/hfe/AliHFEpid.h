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
// Steering class for electron identification
// Combines detector PID objects
// For more information please check the implementation file
//
#ifndef ALIHFEPID_H
#define ALIHFEPID_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

class AliESDtrack;
class AliVParticle;
class AliMCParticle;

class TList;

class AliHFEpid : public TObject{
 public:
    AliHFEpid();
    AliHFEpid(const AliHFEpid &c);
    AliHFEpid &operator=(const AliHFEpid &c);
    ~AliHFEpid();
    
    Bool_t InitializePID(TString argument);
    Bool_t IsSelected(AliHFEpidObject *track);

    Bool_t IsQAOn() const { return TestBit(kIsQAOn); };
    Bool_t HasMCData() const { return TestBit(kHasMCData); };
    void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel; }
    void SetQAOn();
    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData, hasMCdata); };
    TList *GetQAhistograms() const { return fQAlist; };

  protected:
    Bool_t MakePidTpcTof(AliHFEpidObject *track);
    Bool_t MakePidTpcTrd(AliHFEpidObject *track);
    void MakePlotsItsTpc(AliHFEpidObject *track);

    // Stratgies
    void InitStrategy1();
    void InitStrategy2();
    void InitStrategy3();
    void InitStrategy4();
    void InitStrategy5();
    Bool_t IdentifyStrategy1(AliHFEpidObject *track);
    Bool_t IdentifyStrategy2(AliHFEpidObject *track);
    Bool_t IdentifyStrategy3(AliHFEpidObject *track);
    Bool_t IdentifyStrategy4(AliHFEpidObject *track);
    Bool_t IdentifyStrategy5(AliHFEpidObject *track);
  private:
    enum{
      kIsQAOn = BIT(14),
      kHasMCData = BIT(15)
    };
    enum{
      kMCpid = 0,
      kESDpid = 1,
      kITSpid = 2,
      kTPCpid = 3,
      kTRDpid = 4,
      kTOFpid = 5,
      kNdetectorPID = 6
    };
    enum{
      kCombinedTPCTRD=0
    };
    enum{
      kTRDSignal = 0,
      kITSSignal = 1
    };
    AliHFEpidBase *fDetectorPID[kNdetectorPID];     //! Detector PID classes
    UInt_t fEnabledDetectors;                       //  Enabled Detectors
    UInt_t fPIDstrategy;                            //  PID Strategy
    TList *fQAlist;                                 //! QA histograms
    Int_t fDebugLevel;                              //  Debug Level

  ClassDef(AliHFEpid, 1)      // Steering class for Electron ID
};

#endif
