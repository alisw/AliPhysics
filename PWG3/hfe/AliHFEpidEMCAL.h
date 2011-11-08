#ifndef ALIHFEPIDEMCAL_H
#define ALIHFEPIDEMCAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice   */   

//
// Class for EMCAL PID
// electron selection with energy-momnetum matching (e/p)
// For more information please check the implementation file
//
#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

class AliVParticle;
class AliPID;

class AliHFEpidQAmanager;

class AliHFEpidEMCAL : public AliHFEpidBase{
  public:
    AliHFEpidEMCAL();
    AliHFEpidEMCAL(const Char_t *name);
    virtual ~AliHFEpidEMCAL();
    AliHFEpidEMCAL(const AliHFEpidEMCAL &c);
    AliHFEpidEMCAL &operator=(const AliHFEpidEMCAL &c);
  
    virtual Bool_t    InitializePID(Int_t /*run*/);
    virtual Int_t     IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *piqa) const;
      
    //Double_t MomentumEnergyMatchV1(const AliVParticle *track) const;
    Double_t MomentumEnergyMatchV2(const AliVParticle *track) const;
    Double_t CalEopCutMax(const AliVParticle *const track, Int_t flageop) const;
    Double_t CalEopCutMim(const AliVParticle *const track, Int_t flageop) const;

    void SetEoPMax(Float_t eopmax) {feopMax = eopmax;}
    void SetEoPMim(Float_t eopmim) {feopMim = eopmim;}


  protected:
    void Copy(TObject &ref) const;
  private:
    AliPID        *fPID;           //! PID Object
    Float_t    feopMim;         // EMCAL eop mim. cut
    Float_t    feopMax;         // EMCAL eop max. cut

    ClassDef(AliHFEpidEMCAL, 1)
};

#endif
