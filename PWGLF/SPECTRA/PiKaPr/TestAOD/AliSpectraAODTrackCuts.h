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

#include "TNamed.h"

class AliSpectraAODTrackCuts : public TNamed
{
public:

   enum { kTrkBit = 0, kTrkEta, kTrkDCA, kTrkP, kTrkPt, kNTrkCuts};


   AliSpectraAODTrackCuts() : TNamed(), fIsSelected(0), fTrackBits(0), fEtaCut(0), fPCut(0), fPtCut(0), fHistoCuts(0), fTrack(0) {}

   AliSpectraAODTrackCuts(const char *name);
   virtual  ~AliSpectraAODTrackCuts() {} // To be implemented

   Bool_t IsSelected(AliAODTrack * track);

   void SetTrackType(UInt_t bit);
   Bool_t CheckTrackType();
   Bool_t CheckEtaCut();
   Bool_t CheckDCACut();
   Bool_t CheckPCut();
   Bool_t CheckPtCut();
   void PrintCuts() const;

   UInt_t GetTrackType()  const    { return fTrackBits; }
   TH1I * GetHistoCuts()      { return fHistoCuts; }
   void SetEta(Float_t eta)   { fEtaCut = eta; }
   void SetDCA(Float_t dca)   { fDCACut = dca; }
   void SetP(Float_t p)       { fPCut = p; }
   void SetPt(Float_t pt)     { fPtCut = pt; }
   Float_t GetEta()       const    { return fEtaCut; }
   Float_t GetDCA()       const    { return fDCACut; }
   Float_t GetP()         const    { return fPCut; }
   Float_t GetPt()        const    { return fPtCut; }

  Long64_t Merge(TCollection* list);


private:

   Bool_t         fIsSelected;      // True if cuts are selected
   UInt_t         fTrackBits;       // Type of track to be used
   Float_t        fEtaCut;          // Allowed absolute maximum value of Eta
   Float_t        fDCACut;          // Maximum value of DCA
   Float_t        fPCut;            // Maximum value of P
   Float_t        fPtCut;           // Maximum value of Pt

   TH1I *         fHistoCuts;       // Cuts statistics
   AliAODTrack *  fTrack;           //! Track pointer
   
   AliSpectraAODTrackCuts(const AliSpectraAODTrackCuts&);
   AliSpectraAODTrackCuts& operator=(const AliSpectraAODTrackCuts&);

   ClassDef(AliSpectraAODTrackCuts, 1);
};
#endif

