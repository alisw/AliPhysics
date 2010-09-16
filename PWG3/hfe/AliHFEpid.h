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

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

#include <climits>

class AliESDpid;
class AliESDtrack;
class AliHFEpidBase;
class AliVParticle;
class AliMCParticle;

class TList;

class AliHFEpid : public TNamed{
 public:
    enum{
      kUndefined = UINT_MAX 
    };
    enum DETtype_t {
      kMCpid = 0,
      kESDpid = 1,
      kITSpid = 2,
      kTPCpid = 3,
      kTRDpid = 4,
      kTOFpid = 5,
      kNdetectorPID = 6
    };
    AliHFEpid();
    AliHFEpid(const Char_t *name);
    AliHFEpid(const AliHFEpid &c);
    AliHFEpid &operator=(const AliHFEpid &c);
    void Copy(TObject &o) const;
    ~AliHFEpid();
    
    Bool_t InitializePID(TString argument);
    Bool_t IsSelected(AliHFEpidObject *track);

    Bool_t IsQAOn() const { return TestBit(kIsQAOn); };
    Bool_t HasMCData() const { return TestBit(kHasMCData); };
    void SetESDpid(AliESDpid *pid);
    void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel; }
    void SetQAOn();
    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData, hasMCdata); };
    TList *GetQAhistograms() const { return fQAlist; };
    AliHFEpidBase *GetDetPID(DETtype_t det) const { return det < kNdetectorPID ? fDetectorPID[det] : NULL; }
    void PrintStatus() const;

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
    void InitStrategy6();
    void InitStrategy7();
    void InitStrategy8();
    Bool_t IdentifyStrategy0(AliHFEpidObject *track);
    Bool_t IdentifyStrategy1(AliHFEpidObject *track);
    Bool_t IdentifyStrategy2(AliHFEpidObject *track);
    Bool_t IdentifyStrategy3(AliHFEpidObject *track);
    Bool_t IdentifyStrategy4(AliHFEpidObject *track);
    Bool_t IdentifyStrategy5(AliHFEpidObject *track);
    Bool_t IdentifyStrategy6(AliHFEpidObject *track);
    Bool_t IdentifyStrategy7(AliHFEpidObject *track);
    Bool_t IdentifyStrategy8(AliHFEpidObject *track);
  private:
    enum{
      kIsQAOn = BIT(14),
      kHasMCData = BIT(15)
    };
    enum{
      kCombinedTPCTRD=0
    };
    enum{
      kTRDSignal = 0,
      kITSSignal = 1
    };

    void AddCommonObject(TObject * const o);
    void ClearCommonObjects();
    //-----Switch on/off detectors in PID sequence------
    void SwitchOnDetector(UInt_t det){ 
      if(det < kNdetectorPID) SETBIT(fEnabledDetectors, det);
    }
    void SwitchOffDetector(UInt_t det){
      if(det < kNdetectorPID) CLRBIT(fEnabledDetectors, det);
    }
    Bool_t IsDetectorOn(UInt_t det) const {
      return det < kNdetectorPID ? TESTBIT(fEnabledDetectors, det): kFALSE;
    }
    //--------------------------------------------------

    AliHFEpidBase *fDetectorPID[kNdetectorPID];     //! Detector PID classes
    UInt_t fEnabledDetectors;                       //  Enabled Detectors
    UInt_t fPIDstrategy;                            //  PID Strategy
    TList *fQAlist;                                 //! QA histograms
    Int_t fDebugLevel;                              //  Debug Level
    TObjArray *fCommonObjects;                       // Garbage Collector

  ClassDef(AliHFEpid, 1)      // Steering class for Electron ID
};

#endif
