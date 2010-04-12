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
// Class for PID QA
// Several studies done on clean samples of electrons, pions and kaons
// coming from V0 PID
// More information can be found in the implementation file
//
#ifndef ALIHFEPIDQA_H
#define ALIHFEPIDQA_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TList;
class TObjArray;

class AliAODTrack;
class AliESDtrack;
class AliHFEcollection;
class AliHFEV0pid;
class AliMCEvent;
//class AliTRDPIDResponseLQ1D;
class AliVEvent;
class AliESDpid;
class AliHFEV0pidMC;

class AliHFEpidQA : public TObject{
  public:
    AliHFEpidQA();
    ~AliHFEpidQA();

    void Init();
    void Process(AliVEvent *inputEvent);

    TList *GetOutput();
    TList *GetV0pidQA();
    TList *GetV0pidMC();

    Bool_t   HasV0pidQA() const { return TestBit(kV0pidQA); };
    Bool_t   HasRecalculateTRDpid() const { return TestBit(kRecalculateTRDpid); };

    void     SetMCEvent(AliMCEvent * const mc) { fMC = mc; };
    void     SetV0pidQA(Bool_t v0pidQA = kTRUE) { SetBit(kV0pidQA, v0pidQA); };
    void     SetRecalculateTRDpid(Bool_t recal = kTRUE) { SetBit(kRecalculateTRDpid, recal); };

    void     SetRun(Int_t run) { fRun = run; };
    // temporary solutions for correction the T0 for pass4 & pass5
    void     CorrectT0();
    void     SetT0(Float_t t0) { fT0 = t0; };
    Float_t  TOFbeta(AliESDtrack *const track) const;

  protected:
    enum{
      kV0pidQA = BIT(14),
      kRecalculateTRDpid = BIT(15)
    };
    enum{  // detectors for histogram names
      kITS = 0,
      kTPC = 1,
      kTRD = 2,
      kTOF = 4
    };

    void MakePurity(TObjArray *tracks, Int_t species);
    void FillTRDelectronLikelihoods(TObjArray * const particles, Int_t species);
    void FillPIDresponse(TObjArray * const particles, Int_t species);
    void RecalculateTRDpid(AliESDtrack *track, Double_t *pidProbs) const;
    void RecalculateTRDpid(AliAODTrack *track, Double_t *pidProbs) const;

  private:
    AliHFEpidQA(const AliHFEpidQA &ref);
    AliHFEpidQA &operator=(const AliHFEpidQA &ref);

    AliMCEvent        *fMC;           // MC Event
    AliHFEV0pid       *fV0pid;        // V0 PID 
    AliHFEV0pidMC     *fV0pidMC;      // V0 MC PID
    //AliTRDPIDResponseLQ1D *fTRDpidResponse;   // TRD PID
    AliHFEcollection  *fOutput;       // Output container
    Float_t            fT0;           // corrected T0 for pass4 & pass5
    Int_t              fRun;          // Run Number
    AliESDpid *fESDpid;               // ESD PID object
  
  ClassDef(AliHFEpidQA, 1)            // PID QA tool
};
#endif
