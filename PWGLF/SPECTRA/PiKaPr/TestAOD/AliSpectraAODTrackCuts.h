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

class AliAODEvent;
class AliSpectraAODHistoManager;
#include "AliAODMCParticle.h"
#include "AliSpectraAODHistoManager.h"
#include "TNamed.h"

using namespace AliSpectraNameSpace;
 
const char * kBinLabel[] ={"TrkBit",
			   "TrkEta",
			   "TrkDCA",
			   "TrkP",
			   "TrkPt",
			   "TrkPtTOF",
			   "TOFMatching",
			   "kTOFout",
			   "kTIME",
			   "kTOFpid",
			   "Accepted"};

class AliSpectraAODTrackCuts : public TNamed
{
 public:
  
  enum { kTrkBit = 0, kTrkEta, kTrkDCA, kTrkP, kTrkPt,kTrkPtTOF,kTOFMatching,kTrTOFout,kTrTIME,kTrTOFpid,kAccepted,kNTrkCuts};
  
  
 AliSpectraAODTrackCuts() : TNamed(), fIsSelected(0), fTrackBits(0), fEtaCut(0), fPCut(0), fPtCut(0), fPtCutTOFMatching(0), fQvecCutMin(0), fQvecCutMax(0), fHistoCuts(0), fHistoNSelectedPos(0), fHistoNSelectedNeg(0), fHistoNMatchedPos(0), fHistoNMatchedNeg(0), fTrack(0) {}
  
  AliSpectraAODTrackCuts(const char *name);
  virtual  ~AliSpectraAODTrackCuts() {} // To be implemented
  
  Bool_t IsSelected(AliAODTrack * track);
  
  void SetTrackType(UInt_t bit);
  Bool_t CheckTrackType();
  Bool_t CheckEtaCut();
  Bool_t CheckYCut(AODParticleSpecies_t specie); // not included in standard cuts
  Bool_t CheckDCACut();
  Bool_t CheckPCut();
   Bool_t CheckPtCut();
   Bool_t CheckTOFMatching();
   void PrintCuts() const;

   UInt_t GetTrackType()  const    { return fTrackBits; }
   TH1I * GetHistoCuts()      { return fHistoCuts; }
   TH1F * GetHistoNSelectedPos()      { return fHistoNSelectedPos; } 
   TH1F * GetHistoNSelectedNeg()      { return fHistoNSelectedNeg; }
   TH1F * GetHistoNMatchedPos()      { return fHistoNMatchedPos; }
   TH1F * GetHistoNMatchedNeg()      { return fHistoNMatchedNeg; }
   void SetEta(Float_t eta)   { fEtaCut = eta; }
   void SetDCA(Float_t dca)   { fDCACut = dca; }
   void SetP(Float_t p)       { fPCut = p; }
   void SetPt(Float_t pt)     { fPtCut = pt; }
   void SetY(Float_t y) { fYCut = y;}
   void SetPtTOFMatching(Float_t pt)     { fPtCutTOFMatching = pt; }
   void SetQvecMin(Float_t qvecmin)     { fQvecCutMin = qvecmin; }
   void SetQvecMax(Float_t qvecmax)     { fQvecCutMax = qvecmax; }
   Float_t GetEta()       const    { return fEtaCut; }
   Float_t GetY()         const    { return fYCut; }
   Float_t GetDCA()       const    { return fDCACut; }
   Float_t GetP()         const    { return fPCut; }
   Float_t GetPt()        const    { return fPtCut; }
   Float_t GetPtTOFMatching()        const    { return fPtCutTOFMatching; }
   Float_t GetQvecMin()        const    { return fQvecCutMin; }
   Float_t GetQvecMax()        const    { return fQvecCutMax; }
    
   Long64_t Merge(TCollection* list);
   
   
 private:
   
   Bool_t           fIsSelected;      // True if cuts are selected
   UInt_t           fTrackBits;       // Type of track to be used
   Float_t          fEtaCut;          // Allowed absolute maximum value of Eta
   Float_t          fDCACut;          // Maximum value of DCA
   Float_t          fPCut;            // Maximum value of P
   Float_t          fPtCut;           // Maximum value of Pt
   Float_t          fYCut;           // Maximum value of Y
   Float_t          fPtCutTOFMatching;           // TOF Matching
   Float_t          fQvecCutMin;           // Minimum value of Qvec
   Float_t          fQvecCutMax;           // Minimum value of Qvec
   TH1I             *fHistoCuts;       // Cuts statistics
   TH1F             *fHistoNSelectedPos;       // Selected positive tracks
   TH1F             *fHistoNSelectedNeg;       // Selected negative tracks
   TH1F             *fHistoNMatchedPos;       // Matched positive tracks
   TH1F             *fHistoNMatchedNeg;       // Matched negative tracks
   AliAODTrack      *fTrack;           //! Track pointer
   
   AliSpectraAODTrackCuts(const AliSpectraAODTrackCuts&);
   AliSpectraAODTrackCuts& operator=(const AliSpectraAODTrackCuts&);
   
   ClassDef(AliSpectraAODTrackCuts, 1);
};
#endif

