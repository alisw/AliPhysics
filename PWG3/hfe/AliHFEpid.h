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
//#include "AliPIDResponse.h"

class AliHFEcontainer;
class AliHFEvarManager;
class AliPIDResponse;
class AliHFEpidBase;
class AliVParticle;
class AliMCParticle;

class TList;

class AliHFEpid : public TNamed{
 public:
    enum{
      kUndefined = UINT_MAX 
    };
    enum EDETtype_t {
      kMCpid = 0,
      kBAYESpid = 1,
      kITSpid = 2,
      kTPCpid = 3,
      kTRDpid = 4,
      kTOFpid = 5,
      kEMCALpid = 6,
      kNdetectorPID = 7
    };
    AliHFEpid();
    AliHFEpid(const Char_t *name);
    AliHFEpid(const AliHFEpid &c);
    AliHFEpid &operator=(const AliHFEpid &c);
    void Copy(TObject &o) const;
    ~AliHFEpid();
    
    Bool_t InitializePID(Int_t run = 0);
    Bool_t IsSelected(const AliHFEpidObject * const track, AliHFEcontainer *cont = NULL, const Char_t *contname = "trackContainer", AliHFEpidQAmanager *qa = NULL);

    Bool_t HasMCData() const { return TestBit(kHasMCData); };

    void AddDetector(TString detector, UInt_t position);
    void SetPIDResponse(const AliPIDResponse * const pid);
    void SetVarManager(AliHFEvarManager *vm) { fVarManager = vm; }
    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData, hasMCdata); };

    const AliPIDResponse *GetPIDResponse() const;
    UInt_t GetNumberOfPIDdetectors() const { return fNPIDdetectors; }
    Bool_t HasDetector(EDETtype_t det) const { return IsDetectorOn(det); }
    Bool_t IsInitialized() const { return TestBit(kIsInit); }
    void SortDetectors();
    AliHFEpidBase *GetDetPID(EDETtype_t det) const { return det < kNdetectorPID ? fDetectorPID[det] : NULL; }

    void PrintStatus() const;
    const Char_t *SortedDetectorName(Int_t det) {
      if(!TestBit(kDetectorsSorted)) SortDetectors();
      if(det < kNdetectorPID) return fgkDetectorName[fSortedOrder[det]]; 
      else return fgkDetectorName[kNdetectorPID];
    }    
    //-----Configure PID detectors with predefined stettings------
    void ConfigureTOF(Float_t TOFcut = 3.);
    void ConfigureTPCasymmetric(Double_t pmin = 0.1, Double_t pmax = 20., Double_t sigmamin = -0.2, Double_t sigmamax = 5.);
    void ConfigureTPCrejectionSimple();
    void ConfigureTPCcentralityCut(Int_t centralityBin, const char *lowerCutParam = NULL, const Double_t * const params = NULL, Float_t upperTPCCut=3.0);
    void ConfigureTPCdefaultCut(const char *lowerCutParam = NULL, const Double_t * const params = NULL, Float_t upperTPCCut=3.0);
    void ConfigureBayesDetectorMask(Int_t detmask = 10);
    void ConfigureBayesPIDThreshold(Float_t pidthres = 0.9);
    //------------------------------------------------------------

  protected:
    Bool_t MakePidTpcTof(AliHFEpidObject *track);
  private:
    enum{
      kHasMCData = BIT(14),
      kIsInit = BIT(15),
      kDetectorsSorted = BIT(16)
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
    void ConfigureTPCcut(Int_t centralityBin, const char *lowerCutParam, const Double_t * const params, Float_t upperTPCCut);

    static const Char_t *fgkDetectorName[kNdetectorPID + 1]; // PID Detector Names
    AliHFEpidBase *fDetectorPID[kNdetectorPID];     //   Detector PID classes
    UInt_t fDetectorOrder[kNdetectorPID];           //   Position requested by the user
    UInt_t fSortedOrder[kNdetectorPID];             //   Sorted array of detectorIDs
    UInt_t fEnabledDetectors;                       //   Enabled Detectors
    UInt_t fNPIDdetectors;                          //   Number of PID detectors
    AliHFEvarManager *fVarManager;                  //!  HFE Var Manager
    TObjArray *fCommonObjects;                      //   Garbage Collector

  ClassDef(AliHFEpid, 1)      // Steering class for Electron ID
};

#endif
