#ifndef ALIHFEPIDTOF_H
#define ALIHFEPIDTOF_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice   */   

//
// Class for TOF PID
// Rejects protons and kaons at the TPC dE/dx line crossings
// For more information please check the implementation file
//
#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

#ifndef ROOT_TARRAY
#include <TArrayD.h>
#endif

class AliVParticle;
class AliVTrack;
class AliPID;

class AliHFEpidQAmanager;

class AliHFEpidTOF : public AliHFEpidBase{
  public:
    AliHFEpidTOF();
    AliHFEpidTOF(const Char_t *name);
    virtual ~AliHFEpidTOF();
    AliHFEpidTOF(const AliHFEpidTOF &c);
    AliHFEpidTOF &operator=(const AliHFEpidTOF &c);
  
    virtual Bool_t    InitializePID(Int_t /*run*/);
    virtual Int_t     IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *piqa) const;
  
    void SetTOFnSigma(Float_t nSigma) { fNsigmaTOF = nSigma; };
    void SetTOFnSigmaBand(Float_t lower, Float_t upper);
    void SetTOFnSigmaBandCentrality(Float_t lower, Float_t upper, Int_t centralityBin); 
    void SetGenerateTOFmismatch(Bool_t gen = kTRUE, Int_t ntrk = 10) { fGenerateTOFmismatch = gen; fNmismatchTracks = ntrk; }
    Bool_t IsGenerateTOFmismatch() const { return fGenerateTOFmismatch; }
    Int_t GetNmismatchTracks() const { return fNmismatchTracks; }
    void UseTOFonlyIfAvailable() { fUseOnlyIfAvailable = kTRUE; }
    void SetRejectTOFmismatch() { fRejectMismatch = kTRUE; }
    void GenerateTOFmismatch(const AliVTrack * const trk, int ntrk, TArrayD &sigmaEl);

  protected:
    void Copy(TObject &ref) const;
    Bool_t IsMismatch(const AliVTrack *const track) const;

  private:
    enum {
      kSigmaBand = BIT(15)
    };
    Float_t    fNsigmaTOF;          // TOF sigma band
    Float_t    fSigmaBordersTOFLower[12]; // Min.  sigma cut
    Float_t    fSigmaBordersTOFUpper[12]; // Max.  sigma cut
    Bool_t     fUseOnlyIfAvailable;       // Use TOF obly if available
    Bool_t     fRejectMismatch;           // Reject TOF mismatch
    Bool_t     fGenerateTOFmismatch;      // Generate TOF mismatch
    Int_t      fNmismatchTracks;          // Number of mismatch tracks to generate

    ClassDef(AliHFEpidTOF, 1)
};

#endif
