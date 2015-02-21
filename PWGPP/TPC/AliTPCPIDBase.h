#ifndef ALITPCPIDBASE_H
#define ALITPCPIDBASE_H

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

#include <TTreeStream.h>
#include "AliInputEventHandler.h"
#include "AliTOFPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"

class AliTPCPIDBase : public AliAnalysisTaskSE {
 public:
  enum PileUpRejectionType { kPileUpRejectionOff = 0, kPileUpRejectionSPD = 1, kPileUpRejectionMV = 2 };
  enum TPCcutType { kNoCut = 0, kTPCCutMIGeo = 1, kTPCnclCut = 2 };
  AliTPCPIDBase();
  AliTPCPIDBase(const char *name);
  virtual ~AliTPCPIDBase();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  
  virtual Bool_t GetVertexIsOk(AliVEvent* event, Bool_t doVtxZcut = kTRUE) const;
 
  virtual Bool_t GetIsPileUp(AliVEvent* event, PileUpRejectionType pileUpRejectionType) const;
  
  virtual Bool_t GetIsPbpOrpPb() const { return fIsPbpOrpPb; };
  virtual void SetIsPbpOrpPb(Bool_t newValue) { fIsPbpOrpPb = newValue; };
  
  virtual Double_t GetZvtxCutEvent() const { return fZvtxCutEvent; };
  virtual void SetZvtxCutEvent(Double_t newValue) { fZvtxCutEvent = newValue; if (fAnaUtils) fAnaUtils->SetMaxVtxZ(fZvtxCutEvent);};
  
  virtual Bool_t GetUsePhiCut() const { return fUsePhiCut; };
  virtual void SetUsePhiCut(Bool_t newValue) { fUsePhiCut = newValue; };
  
  virtual TPCcutType GetTPCcutType() const { return fTPCcutType; };
  virtual Bool_t GetUseTPCCutMIGeo() const { return (fTPCcutType == kTPCCutMIGeo); };
  virtual Bool_t GetUseTPCnclCut() const { return (fTPCcutType == kTPCnclCut); };
  
  virtual void SetTPCcutType(TPCcutType newType) { fTPCcutType = newType; };
  
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
  
  static Double_t GetCutGeo() { return fgCutGeo; };
  static Double_t GetCutNcr() { return fgCutNcr; };
  static Double_t GetCutNcl() { return fgCutNcl; };
  
  static void SetCutGeo(Double_t value) { fgCutGeo = value; };
  static void SetCutNcr(Double_t value) { fgCutNcr = value; };
  static void SetCutNcl(Double_t value) { fgCutNcl = value; };
  
  static Bool_t TPCCutMIGeo(const AliVTrack* track, const AliVEvent* evt, TTreeStream* streamer = 0x0);
  static Bool_t TPCCutMIGeo(const AliVTrack* track, const AliInputEventHandler* evtHandler, TTreeStream* streamer = 0x0)
    { if (!evtHandler) return kFALSE; return TPCCutMIGeo(track, evtHandler->GetEvent(), streamer); };
  
  static UShort_t GetCutPureNcl() { return fgCutPureNcl; };
  static void SetCutPureNcl(UShort_t value) { fgCutPureNcl = value; };
  
  static Bool_t TPCnclCut(const AliVTrack* track);
  
 protected:
  void FillV0PIDlist(AliESDEvent* esdEvent = 0x0);
  void ClearV0PIDlist();
  
  static Double_t fgCutGeo;   // Cut variable for TPCCutMIGeo concerning geometry
  static Double_t fgCutNcr; // Cut variable for TPCCutMIGeo concerning num crossed rows
  static Double_t fgCutNcl;  // Cut variable for TPCCutMIGeo concerning num clusters
  
  static UShort_t fgCutPureNcl; // Cut variable for TPCnclCut
  
  AliVEvent   *fEvent;    //! VEvent object
  AliESDEvent *fESD;      //! ESDEvent object, if ESD
  AliMCEvent  *fMC;       //! MC object

  AliPIDResponse *fPIDResponse;    //! PID response Handler
  AliESDv0KineCuts *fV0KineCuts;       //! ESD V0 kine cuts
  
  AliAnalysisUtils *fAnaUtils; //! Object to use analysis utils like pile-up rejection
  
  Bool_t fIsPbpOrpPb;       // Pbp/pPb collision or something else?
  Bool_t fUsePhiCut;        // Use cut on phi (useful for TPC)
  TPCcutType fTPCcutType;   // Type of TPC cut to be used
  Double_t fZvtxCutEvent;   // Vertex z cut for the event (cm)
  Double_t fEtaCut;         // Eta cut
  
  TF1* fPhiCutLow;          // phi prime cut, low
  TF1* fPhiCutHigh;         // phi prime cut, high
  
  TRandom3* fRandom;        //! Can be used to statistically determine the shape in the pt bins e.g.
  
  AliAnalysisFilter* fTrackFilter; // Track Filter
  

  Int_t fNumTagsStored;     // Number of entries of fV0tags
  Char_t* fV0tags;         //! Pointer to array with tags for identified particles from V0 decays
  
  Bool_t fStoreMotherIndex; // Switch on/off storing the mother indices of V0 daughters
  Int_t* fV0motherIndex; //! Pointer to array with index of the mother V0
  
 private:
  
  AliTPCPIDBase(const AliTPCPIDBase&); // not implemented
  AliTPCPIDBase& operator=(const AliTPCPIDBase&); // not implemented
  
  ClassDef(AliTPCPIDBase, 2);
};



inline Float_t AliTPCPIDBase::GetDeltaTOF(const AliVTrack *track, const AliTOFPIDResponse* tofPIDresponse,
                                                     const Double_t* times, AliPID::EParticleType type) const
{
  return (track->GetTOFsignal() - tofPIDresponse->GetStartTime(track->P()) - times[type]);
}

#endif
