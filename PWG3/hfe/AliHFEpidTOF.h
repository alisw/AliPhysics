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
  
    virtual Bool_t    InitializePID();
    virtual Int_t     IsSelected(AliHFEpidObject *track, AliHFEpidQAmanager *piqa);
  
    void SetTOFnSigma(Short_t nSigma) { fNsigmaTOF = nSigma; };

  protected:
    void Copy(TObject &ref) const;
    Double_t NumberOfSigmas(const AliVParticle *track, AliPID::EParticleType species, AliHFEpidObject::AnalysisType_t anaType);
  private:
    AliPID        *fPID;           //! PID Object

    Short_t    fNsigmaTOF;         // TOF sigma band

    ClassDef(AliHFEpidTOF, 1)
};

#endif
