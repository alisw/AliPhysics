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
//class AliEMCALGeoUtils;
//class AliEMCALGeometry;
//class AliEMCALRecoUtils;
/*
class AliEMCALGeoUtils;
class AliEMCALGeometry;
#include "AliEMCALGeoParams.h"
class AliEMCALRecoUtils;
class AliEMCALPIDUtils;
*/


class AliHFEpidEMCAL : public AliHFEpidBase{
  public:
    AliHFEpidEMCAL();
    AliHFEpidEMCAL(const Char_t *name);
    virtual ~AliHFEpidEMCAL();
    AliHFEpidEMCAL(const AliHFEpidEMCAL &c);
    AliHFEpidEMCAL &operator=(const AliHFEpidEMCAL &c);
  
    virtual Bool_t    InitializePID();
    virtual Int_t     IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *piqa) const;
      
    //Double_t MomentumEnergyMatch(const AliVParticle *track) const;

  protected:
    void Copy(TObject &ref) const;
  private:
    AliPID        *fPID;           //! PID Object
    Float_t    feopMim;         // EMCAL eop mim. cut
    Float_t    feopMax;         // EMCAL eop max. cut

    ClassDef(AliHFEpidEMCAL, 1)
};

#endif
