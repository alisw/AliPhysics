#ifndef ALIANALYSISTASKPIDV0BASE_H
#define ALIANALYSISTASKPIDV0BASE_H

/*
This class is a base class for all other
analysis tasks that use V0's.
It provides basics for V0 identification.
In addition, some further basic functions are provided.

Class written by Benjamin Hess.
Contact: bhess@cern.ch
*/

class TF1;
class TRandom3;
class TObjArray;
class AliVEvent;
class AliESDEvent;
class AliMCEvent;
class AliPIDResponse;
class AliESDv0KineCuts;
class AliPID;
class AliAnalysisFilter;
class AliVTrack;

#include "AliTOFPIDResponse.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPIDV0base : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPIDV0base();
  AliAnalysisTaskPIDV0base(const char *name);
  virtual ~AliAnalysisTaskPIDV0base();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  
  virtual Bool_t GetVertexIsOk(AliVEvent* event) const;
  
  virtual Bool_t GetIsPbpOrpPb() const { return fIsPbpOrpPb; };
  virtual void SetIsPbpOrpPb(Bool_t newValue) { fIsPbpOrpPb = newValue; };
  
  virtual Double_t GetZvtxCutEvent() const { return fZvtxCutEvent; };
  virtual void SetZvtxCutEvent(Double_t newValue) { fZvtxCutEvent = newValue; };
  
  virtual Bool_t GetUsePhiCut() const { return fUsePhiCut; };
  virtual void SetUsePhiCut(Bool_t newValue) { fUsePhiCut = newValue; };
  
  virtual Double_t GetEtaCut() const { return fEtaCut; };     
  virtual void  SetEtaCut(Double_t etaCut){ fEtaCut = etaCut; };
  
  virtual const AliAnalysisFilter* GetTrackFilter() const { return fTrackFilter; };
  virtual void  SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
  
  virtual Char_t GetV0tag(Int_t trackIndex) const;
  
  virtual Bool_t GetStoreMotherIndex() const { return fStoreMotherIndex; };
  virtual void SetStoreMotherIndex(Bool_t newValue) { fStoreMotherIndex = newValue; };
  
  virtual Int_t GetV0motherIndex(Int_t trackIndex) const;
  
  virtual Double_t GetPhiPrime(Double_t phi, Double_t magField, Int_t charge) const;
  virtual Bool_t PhiPrimeCut(const AliVTrack* track, Double_t magField) const;
  virtual Bool_t PhiPrimeCut(Double_t trackPt, Double_t trackPhi, Short_t trackCharge, Double_t magField) const;
  virtual Float_t GetDeltaTOF(const AliVTrack *track, const AliTOFPIDResponse* tofPIDresponse, const Double_t* times, 
                              AliPID::EParticleType type) const;
  
 protected:
  
  AliVEvent   *fEvent;    //! VEvent object
  AliESDEvent *fESD;      //! ESDEvent object, if ESD
  AliMCEvent  *fMC;       //! MC object

  AliPIDResponse *fPIDResponse;    //! PID response Handler
  AliESDv0KineCuts *fV0KineCuts;       //! ESD V0 kine cuts
  
  Bool_t fIsPbpOrpPb;       // Pbp/pPb collision or something else?
  Bool_t fUsePhiCut;        // Use cut on phi (useful for TPC)
  Double_t fZvtxCutEvent;   // Vertex z cut for the event (cm)
  Double_t fEtaCut;         // Eta cut
  
  TF1* fPhiCutLow;          // phi prime cut, low
  TF1* fPhiCutHigh;         // phi prime cut, high
  
  TRandom3* fRandom;        //! Can be used to statistically determine the shape in the pt bins e.g.
  
  AliAnalysisFilter* fTrackFilter; //! Track Filter
  

  Int_t fNumTagsStored;     // Number of entries of fV0tags
  Char_t* fV0tags;         //! Pointer to array with tags for identified particles from V0 decays
  
  Bool_t fStoreMotherIndex; // Switch on/off storing the mother indices of V0 daughters
  Int_t* fV0motherIndex; //! Pointer to array with index of the mother V0
  
 private:
  void FillV0PIDlist(AliESDEvent* esdEvent = 0x0);
  void ClearV0PIDlist();
  
  AliAnalysisTaskPIDV0base(const AliAnalysisTaskPIDV0base&); // not implemented
  AliAnalysisTaskPIDV0base& operator=(const AliAnalysisTaskPIDV0base&); // not implemented
  
  ClassDef(AliAnalysisTaskPIDV0base, 1); // example of analysis
};



inline Float_t AliAnalysisTaskPIDV0base::GetDeltaTOF(const AliVTrack *track, const AliTOFPIDResponse* tofPIDresponse,
                                                     const Double_t* times, AliPID::EParticleType type) const
{
  return (track->GetTOFsignal() - tofPIDresponse->GetStartTime(track->P()) - times[type]);
}

#endif
