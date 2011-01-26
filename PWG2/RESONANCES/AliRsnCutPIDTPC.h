//
// *** Class AliRsnCutPIDTPC ***
//
// This class implements all cuts which have to be used for the 2010 runs
// for phi and generic resonance analysis.
// It contains an AliESDtrackCuts object for track quality selection
// and some criteria for particle identification with ITS, TPC and TOF.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTPIDTPC_H
#define ALIRSNCUTPIDTPC_H

#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include "AliRsnCut.h"

class AliRsnCutPIDTPC : public AliRsnCut
{
  public:

    AliRsnCutPIDTPC(const char *name = "cutPIDTPC", AliPID::EParticleType type = AliPID::kKaon, Double_t momLimit = 0.350, Double_t cut1 = 5.0, Double_t cut2 = 3.0);
    AliRsnCutPIDTPC(const AliRsnCutPIDTPC& copy);
    AliRsnCutPIDTPC& operator=(const AliRsnCutPIDTPC& copy);
    virtual ~AliRsnCutPIDTPC() { }

    AliESDpid*       GetESDpid()                         {return &fESDpid;}    
    void             SetPIDType(AliPID::EParticleType t) {fPIDtype = t;}
        
    void             SetMomentumLimit(Double_t v)        {fMomentumLimit = v;}
    void             SetLargeCut(Double_t v)             {fLargeCut = v;}
    void             SetSmallCut(Double_t v)             {fSmallCut = v;}
    void             SetBBParam(Double_t *p)             {SetBBParam(p[0], p[1], p[2], p[3], p[4]);}
    void             SetBBParam(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4);
    
    virtual Bool_t   IsSelected(TObject *object);
    virtual void     Print(const Option_t *option = "") const;

  protected:
  
    AliPID::EParticleType   fPIDtype;          //  PID reference type used for checks
    
    Double_t                fMomentumLimit;    //  below this value, large cut is used; above, small one is used (value in GeV/c)
    Double_t                fLargeCut;         //  range for TPC de/dx large cut
    Double_t                fSmallCut;         //  range for TPC de/dx small cut
    
    AliESDpid               fESDpid;           //  ESD PID object
    AliAODpidUtil           fAODpid;           //  AOD PID object

    ClassDef(AliRsnCutPIDTPC, 1)
};

#endif
