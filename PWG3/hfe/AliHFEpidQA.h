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
class TPDGCode;
class TParticle;
class TFile;
class TMultiLayerPerceptron;

class AliAODTrack;
class AliESDtrack;
class AliMCEvent;
class AliVEvent;
class AliPIDResponse;
class AliExternalTrackParam;
class AliLog;

class AliHFEV0pidMC;
class AliHFEcollection;
class AliHFEtrdPIDqa;
class AliHFEV0pid;
class AliHFEpidTRD;


class AliHFEpidQA : public TObject{
  public:
    AliHFEpidQA();
    ~AliHFEpidQA();
    AliHFEpidQA(const AliHFEpidQA &ref);
    AliHFEpidQA &operator=(const AliHFEpidQA &ref);
    virtual void Copy(TObject &o) const;

    void Init();
    void Process();

    TList *GetOutput();
    TList *GetV0pidQA();
    TList *GetV0pidMC();

    Bool_t   HasV0pidQA() const { return TestBit(kV0pidQA); };
    Bool_t   HasRecalculateTRDpid() const { return TestBit(kRecalculateTRDpid); };

    void     SetEvent(AliVEvent* const ev) { fEvent = ev; };
    void     SetMCEvent(AliMCEvent * const mc) { fMC = mc; };
    void     SetV0pidQA(Bool_t v0pidQA = kTRUE) { SetBit(kV0pidQA, v0pidQA); };
    void     SetRecalculateTRDpid(Bool_t recal = kTRUE) { SetBit(kRecalculateTRDpid, recal); };
    void     SetTRDTotalChargeInSlice0() { fTRDTotalChargeInSlice0 = kTRUE; }

    void     SetPIDResponse(AliPIDResponse* const pid) { fESDpid = pid; }
    Float_t  TOFbeta(const AliESDtrack* const track) const;

    void     CheckEvent();
    void     SetNNref(TFile *f) { fNNref = f; };
    
    AliHFEtrdPIDqa *GetTRDQA() const { return fTRDpidQA; }

  protected:
    enum{
      kV0pidQA = BIT(14),
      kRecalculateTRDpid = BIT(15)
    };
    enum{  // detectors for histogram names
      kITS = 0,
      kTPC = 1,
      kTRD = 2,
      kTOF = 3
    };

    TObjArray *MakeTrackList(const TObjArray *tracks) const;
    TObjArray *MakeCleanListElectrons(const TObjArray *tracks) const;

    void MakePurity(const TObjArray *tracks, Int_t species);
    TObjArray *MakeCleanListForTRD(const TObjArray * const track, Int_t species);
    void FillElectronLikelihoods(const TObjArray * const particles, Int_t species);
    void FillPIDresponse(const TObjArray * const particles, Int_t species);
    void FillIllumination(const TObjArray *const particles, Int_t species);
    void FillTPCinfo(AliESDtrack * const track, Int_t species);
    void TestTRDResponse(const TObjArray * const tracks, Int_t species);
    void RecalculateTRDpid(AliESDtrack *track, Double_t *pidProbs) const;
    void RecalculateTRDpid(AliAODTrack *track, Double_t *pidProbs) const;
    void CheckTenderV0pid(const TObjArray * const particles, Int_t species);
    Int_t GetTenderV0pid(AliESDtrack * const track);
    
    Double_t TRDlikeTracklet(Int_t layer, AliESDtrack * const track, Double_t *likelihood);
    Int_t TRDmomBin(Double_t p) const;

 protected:
    Int_t GetMaxPID(const Double_t *pidProbs) const;
    Int_t GetPDG(Int_t index);
    
  private:
    AliVEvent         *fEvent;        // event pointer
    AliMCEvent        *fMC;           // MC Event
    AliHFEV0pid       *fV0pid;        // V0 PID 
    AliHFEV0pidMC     *fV0pidMC;      // V0 MC PID
    AliHFEtrdPIDqa    *fTRDpidQA;     //! TRD PID QA object
    AliHFEcollection  *fOutput;       // Output container
    AliPIDResponse    *fESDpid;       // ESD PID object
 private:
    TFile             *fNNref;        // reference file for NN pid 
    TMultiLayerPerceptron *fNet[11];  //  reference networks
    Bool_t fTRDTotalChargeInSlice0;     // Fix for Foreward/Backward compatibility
  
  ClassDef(AliHFEpidQA, 1)            // PID QA tool
};
#endif
