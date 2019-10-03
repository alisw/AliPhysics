#ifndef ALISPECTRAAODTRACKCUTS_H
#define ALISPECTRAAODTRACKCUTS_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraAODTrackCuts
//
//
//
//
// Authors: Michele Floris, CERN, Philip Versteeg, UU, Redmer Bertens, UU
//-------------------------------------------------------------------------

class TH1I;
class AliAODEvent;
class AliPIDResponse;  
class AliAODMCParticle;
class AliAODTrack;
class TH2F;
class TH1F;

#include "TNamed.h"

class AliSpectraAODTrackCuts : public TNamed
{
 public:
  
  enum { kTrkBit = 0, kTrkCuts, kTrkEta, kTrkDCA, kTrkP, kTrkPt,kTrkPtTOF,kTOFMatching,kTrTOFout,kTrTIME,kTrTOFpid,kAccepted,kNTrkCuts};
  
  
 AliSpectraAODTrackCuts() : TNamed(), fIsSelected(0), fTrackBits(0), fMinTPCcls(0), fRequestSPDcls(0), fEtaCutMin(0), fEtaCutMax(0),fDCACut(0), fPCut(0), fPtCut(0),fYCut(0), fPtCutTOFMatching(0), fHistoCuts(0), fHistoNSelectedPos(0), fHistoNSelectedNeg(0), fHistoNMatchedPos(0), fHistoNMatchedNeg(0), fHistoEtaPhiHighPt(0), fTrack(0), fPIDResponse(0) {}
  
  AliSpectraAODTrackCuts(const char *name);
  virtual  ~AliSpectraAODTrackCuts() {} // To be implemented
  
  Bool_t IsSelected(AliAODTrack * track,Bool_t FillHistStat);
  
  void SetEta(Float_t etamin,Float_t etamax)   { fEtaCutMin = etamin;fEtaCutMax = etamax; }
  void SetDCA(Float_t dca)   { fDCACut = dca; }
  void SetP(Float_t p)       { fPCut = p; }
  void SetPt(Float_t pt)     { fPtCut = pt; }
  void SetY(Float_t y) { fYCut = y;}
  void SetPtTOFMatching(Float_t pt)     { fPtCutTOFMatching = pt; }
  void SetTrackType(UInt_t bit);
  void SetTrackBits(UInt_t TrackBits) {fTrackBits=TrackBits;}
  void SetMinTPCcls(UInt_t MinTPCcls) {fMinTPCcls=MinTPCcls;}
  void SetRequestSPDcls(Bool_t RequestSPDcls) {fRequestSPDcls=RequestSPDcls;}
  
  UInt_t GetTrackType()  const    { return fTrackBits;}
  TH1I * GetHistoCuts()      { return fHistoCuts; }
  TH1F * GetHistoNSelectedPos()      { return fHistoNSelectedPos; } 
  TH1F * GetHistoNSelectedNeg()      { return fHistoNSelectedNeg; }
  TH1F * GetHistoNMatchedPos()      { return fHistoNMatchedPos; }
  TH1F * GetHistoNMatchedNeg()      { return fHistoNMatchedNeg; }
  TH2F * GetHistoEtaPhiHighPt()      { return fHistoEtaPhiHighPt; }
  Float_t GetEtaMin()       const    { return fEtaCutMin; }
  Float_t GetEtaMax()       const    { return fEtaCutMax; }
  Float_t GetY()         const    { return fYCut; }
  Float_t GetDCA()       const    { return fDCACut; }
  Float_t GetP()         const    { return fPCut; }
  Float_t GetPt()        const    { return fPtCut; }
  Float_t GetPtTOFMatching()        const    { return fPtCutTOFMatching; }
  
  Bool_t CheckTrackType();
  Bool_t CheckTrackCuts();
  Bool_t CheckEtaCut();
  Bool_t CheckYCut(Double_t mass); // not included in standard cuts
  Bool_t CheckDCACut();
  Bool_t CheckPCut();
  Bool_t CheckPtCut();
  Bool_t CheckTOFMatching(Bool_t FillHistStat);
  void PrintCuts() const;
  
  Long64_t Merge(TCollection* list);
   
   
 private:
  
  Bool_t           fIsSelected;      // True if cuts are selected
  UInt_t           fTrackBits;       // Type of track to be used
  UInt_t           fMinTPCcls;       // min number of clusters in the TPC
  Bool_t           fRequestSPDcls;         // request a hit in the SPD
  Float_t          fEtaCutMin;          // Allowed absolute maximum value of Eta
  Float_t          fEtaCutMax;          // Allowed absolute maximum value of Eta
  Float_t          fDCACut;          // Maximum value of DCA
  Float_t          fPCut;            // Maximum value of P
  Float_t          fPtCut;           // Maximum value of Pt
  Float_t          fYCut;           // Maximum value of Y
  Float_t          fPtCutTOFMatching;           // TOF Matching
  TH1I             *fHistoCuts;       // Cuts statistics
  TH1F             *fHistoNSelectedPos;       // Selected positive tracks
  TH1F             *fHistoNSelectedNeg;       // Selected negative tracks
  TH1F             *fHistoNMatchedPos;       // Matched positive tracks
  TH1F             *fHistoNMatchedNeg;       // Matched negative tracks
  TH2F             *fHistoEtaPhiHighPt;       // EtaPhi distr at high pt (>1.5 GeV/c)
  AliAODTrack      *fTrack;           //! Track pointer
  AliPIDResponse   *fPIDResponse;     // ! PID response object
  static const char * kBinLabel[]; // labels of stat histo

   
  AliSpectraAODTrackCuts(const AliSpectraAODTrackCuts&);
  AliSpectraAODTrackCuts& operator=(const AliSpectraAODTrackCuts&);
   
  ClassDef(AliSpectraAODTrackCuts, 3);
};
#endif

