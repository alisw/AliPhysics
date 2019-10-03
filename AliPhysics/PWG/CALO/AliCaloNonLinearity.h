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
      // pp 7 TeV 2010
      k14j4,
      // pp 2.76 TeV 2011
      k12f1a,
      k12f1b,
      k12i3,
      kPP2T11P4JJ,
      k15g1b,
      // PbPb 2.76 TeV 2011
      k14a1,
      // pp 8 TeV 2012
      k14e2b,
      kPP8T12P2Pyt8,
      kPP8T12P2Pho,
      kPP8T12P2JJ,
      // pPb 5 TeV 2013
      kPPb5T13P2DPMJet,
      kPPb5T13P2HIJAdd,
      k16c3a,
      k16c3b,
      k16c3c,
      // pp 2.76TeV 2013
      k15g2,
      kPP2T13P1JJ,
      k15a3b,
      // pp 13 TeV 2015
      kPP13T15P2Pyt8,
      kPP13T15P2EPOS,
      k15k5,
      // pp 5 TeV 2015
      k16h8a,
      k16h8b,
      k16k3a,
      k16k5a,
      k16k5b,
      k17e2,
      k16h3,
      // pp 5 TeV 2017
      k17l4b,
      k17l3b,
      // PbPb 5 TeV 2015
      kPbPb5T15HIJING,
      k16k3b,
      // pp 13 TeV 2016
      kPP13T16P1Pyt8,
      kPP13T16P1Pyt8LowB,
      kPP13T16P1EPOS,
      kPP13T16P1JJ,
      kPP13T16P1JJLowB,
      k17h8a,
      k17h8b,
      k17h8c,
      k17c3b1,
      k17c3a1,
      k17c3b2,
      k17c3a2,
      // pPb 5 TeV 2016
      kPPb5T16EPOS,
      kPPb5T16DPMJet,
      k17g8a,
      k17d2a,
      k17d2b,
      // pPb 8 TeV 2016
      k17f3a,
      k17f3b,
      k17f4a,
      k17f4b,
      k17g8b,
      k17g8c,
      // pp 13 TeV 2017
      k17k1,
      kPP13T17P1Pyt8,
      kPP13T17P1Pho,
      kPP13T17P1Pyt6,
      kPP13T17P1Pyt8Str,
      kPP13T17P1Pyt8LowB,
      // Xe-Xe MC
      kXeXe5T17HIJING,

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
      k16pp13TeV,
      k16pp13TeVLow,
      k16pPb5023GeV,
      k16pPb8TeV,
      k17pp13TeV,
      k17pp13TeVLow,
      k17pp13TeVNo,
      k17XeXe5440GeV,
      k17pp5TeV
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
    Float_t     FunctionNL_PHOSOnlyMC(Float_t e, Float_t p0, Float_t p1, Float_t p2);

    Float_t     FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2);
    Float_t     FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5);
    Float_t     FunctionNL_DExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5);
    //predefined functions
    Float_t     FunctionNL_kPi0MCv1(Float_t e);
    Float_t     FunctionNL_kPi0MCv2(Float_t e);
    Float_t     FunctionNL_kPi0MCv3(Float_t e);
    Float_t     FunctionNL_kPi0MCv5(Float_t e);
    Float_t     FunctionNL_kPi0MCv6(Float_t e);
    Float_t     FunctionNL_kSDMv5(Float_t e);
    Float_t     FunctionNL_kSDMv6(Float_t e);
    Float_t     FunctionNL_kTestBeamv2(Float_t e);
    Float_t     FunctionNL_kTestBeamv3(Float_t e);
    Float_t     FunctionM02 (Float_t E, Float_t a, Float_t b, Float_t c, Float_t d, Float_t e);
    Bool_t      IsSelected(TList*   /* list */ ) { return kTRUE; }

  protected:

    MCSetEnum   fCurrentMC;                               // enum for current MC set being processed
    Int_t       fClusterType;                             // which cluster do we have

  private:

    ClassDef(AliCaloNonLinearity,2)
};

#endif
