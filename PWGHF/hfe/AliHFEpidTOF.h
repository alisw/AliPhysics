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

class AliVParticle;
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
    void SetTOFnSigmaBand(Float_t lower, Float_t upper) { fSigmaBordersTOF[0] = lower; fSigmaBordersTOF[1] = upper; SetBit(kSigmaBand, kTRUE); }

  protected:
    void Copy(TObject &ref) const;
  private:
    enum {
      kSigmaBand = BIT(15)
    };
    Float_t    fNsigmaTOF;          // TOF sigma band
    Float_t    fSigmaBordersTOF[2]; // Min. and Max TOF sigma cut

    ClassDef(AliHFEpidTOF, 1)
};

#endif
