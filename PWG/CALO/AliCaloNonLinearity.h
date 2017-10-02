#ifndef ALICALONONLINEARITY_H
#define ALICALONONLINEARITY_H
// Class handling Non linearity correction for Calo analysis
// Authors: Nicolas Schmidt, Friederike Bock, Daniel Muehlheim

#include "AliAnalysisCuts.h"
#include "AliVCluster.h"
#include "TF1.h"
#include <vector>


class AliCaloNonLinearity : public AliAnalysisCuts {

  public:

    enum MCSetEnum {
      // MC data sets
      kNoMC=0,
      k14e2a,
      k14e2b,
      k14e2c,
      k12f1a,
      k12f1b,
      k12i3,
      k15g1a,
      k15g1b,
      k15g2,
      k15a3a,
      k15a3a_plus,
      k15a3b,
      k13b2_efix,
      k13e7,
      k15h1,
      k15h2,
      k14j4,
      k16c2,
      k16c2_plus,
      k16c3a,
      k16c3b,
      k16c3c,
      k16h3,
      k16h3b,
      k16h8a,
      k16h8b,
      k16k3a,
      k16k3b,
      k16k5a,
      k16k5b,
      // Data starts here
      k10pp7TeV,
      k10pp900GeV,
      k10PbPb2760GeV,
      k11pp2760GeV,
      k11pp7TeV,
      k11PbPb2760GeV,
      k12pp8TeV,
      k13pPb5023GeV,
      k13pp2760GeV,
      k15pp13TeV,
      k15pp5TeV,
      k15PbPb5TeV,
      k16pp13TeV
    };

    //Constructors
    AliCaloNonLinearity(const char *name="NonLinearity", const char *title="NonLinearity");

    //virtual destructor
    virtual     ~AliCaloNonLinearity();


    //correct NonLinearity
    MCSetEnum   FindEnumForMCSetString(TString namePeriod);
    Float_t     GetCorrectedEnergy(Float_t clusterEnergy, Int_t isMC, Int_t switchNonLin, Int_t clusterType, MCSetEnum periodEnum);

    //predefined functions
    Float_t     FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6);
    Float_t     FunctionNL_PHOS(Float_t e, Float_t p0, Float_t p1, Float_t p2);
    Float_t     FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2);
    Float_t     FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5);
    Float_t     FunctionM02(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4);
    Float_t     FunctionNL_kPi0MCv1(Float_t e);
    Float_t     FunctionNL_kPi0MCv2(Float_t e);
    Float_t     FunctionNL_kPi0MCv3(Float_t e);
    Float_t     FunctionNL_kPi0MCv5(Float_t e);
    Float_t     FunctionNL_kPi0MCv6(Float_t e);
    Float_t     FunctionNL_kSDMv5(Float_t e);
    Float_t     FunctionNL_kSDMv6(Float_t e);
    Float_t     FunctionNL_kTestBeamv2(Float_t e);
    Float_t     FunctionNL_kTestBeamv3(Float_t e);
    Bool_t      IsSelected(TList*   /* list */ ) { return kTRUE; }

  protected:

    MCSetEnum   fCurrentMC;                               // enum for current MC set being processed
    Int_t       fClusterType;                             // which cluster do we have

  private:

    ClassDef(AliCaloNonLinearity,1)
};

#endif
