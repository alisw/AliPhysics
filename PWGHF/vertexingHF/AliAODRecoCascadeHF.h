#ifndef ALIAODRECOCASCADEHF_H
#define ALIAODRECOCASCADEHF_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliAODRecoCascadeHF
// base class for AOD reconstructed heavy-flavour cascade decays
// (D*->D0pi, ...)
// The convention is: prong 0 is the bachelor, prong 1 is the "V0"
//
// Author: X-M. Zhang, zhangxm@iopp.ccnu.edu.cn
//***********************************************************

#include <TRef.h>
#include <TRefArray.h>
#include <TClonesArray.h>
#include <TClass.h>
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliAODRecoDecayHF2Prong.h"

class AliAODRecoCascadeHF : public AliAODRecoDecayHF2Prong {

 public:

  AliAODRecoCascadeHF();
  AliAODRecoCascadeHF(AliAODVertex *vtx2, Short_t charge,
		      Double_t *px, Double_t *py, Double_t *pz,
		      Double_t *d0, Double_t *d0err, Double_t dca);
  AliAODRecoCascadeHF(AliAODVertex *vtx2, Short_t charge,
		      Double_t *d0, Double_t *d0err, Double_t dca);
  virtual ~AliAODRecoCascadeHF();

  AliAODRecoCascadeHF(const AliAODRecoCascadeHF& source);
  AliAODRecoCascadeHF& operator=(const AliAODRecoCascadeHF& source);

  // 2prong (D0 for Dstar)
  AliAODRecoDecayHF2Prong* Get2Prong() const {
    if(!GetDaughter(1)) return 0;
    if ( ! ((AliAODRecoDecay*)GetDaughter(1))->IsA()->InheritsFrom("AliAODRecoDecayHF2Prong") ){
      AliWarning("Object is not of type AliAODRecoDecayHF2Prong");
      return 0;
    }
    return (AliAODRecoDecayHF2Prong*)GetDaughter(1);
  }

  // Bachelor (soft pion for Dstar)
  AliAODTrack* GetBachelor() const {return (AliAODTrack*)GetDaughter(0);}

  // v0 (Ks or Lambda for Lambda_c)
  AliAODv0* Getv0() const {
    if ( ! ((AliAODRecoDecay*)GetDaughter(1))->IsA()->InheritsFrom("AliAODv0") ){
       AliWarning("Object is not of type v0");
       return 0;
      }
    return (AliAODv0*)GetDaughter(1);
    }

  // Get v0 positive track
  AliAODTrack* Getv0PositiveTrack() const { return  (AliAODTrack*)Getv0()->GetDaughter(0);  }
  // Get v0 negative track
  AliAODTrack* Getv0NegativeTrack() const { return  (AliAODTrack*)Getv0()->GetDaughter(1);  }

  // D*->D0pi, D0->Kpi
  Double_t EDstar() const {return E(413);} 
  Double_t YDstar() const {return Y(413);} 
  Bool_t   SelectDstar(const Double_t *cutsDstar,const Double_t *cutsD0,Bool_t testD0=kTRUE) const;
  Double_t InvMassD0() const {return (Charge()>0 ? Get2Prong()->InvMassD0() : Get2Prong()->InvMassD0bar());}
  Double_t InvMassDstarKpipi() const;
  Double_t DeltaInvMass() const {return (InvMassDstarKpipi()-InvMassD0());}
  Double_t AngleD0dkpPisoft() const;
  Bool_t   TrigonometricalCut() const;

  // Lc invariant mass
  Double_t InvMassLctoK0sP() const {
    UInt_t pdg[2]={2212,310}; return InvMass(2,pdg);
  }
  Double_t InvMassLctoLambdaPi() const {
    UInt_t pdg[2]={211,3122}; return InvMass(2,pdg);
  }
  Bool_t SelectLctoV0(const Double_t *cutsLctoV0, Bool_t okLck0sp, Bool_t okLcLpi, Bool_t okLcLbarpi) const;

  Int_t MatchToMC(Int_t pdgabs,Int_t pdgabs2prong,
                  Int_t *pdgDg,Int_t *pdgDg2prong,
                  TClonesArray *mcArray, Bool_t isV0=kFALSE) const;

  Double_t CosV0PointingAngle() const;
  Double_t DecayLengthV0() const;
  Double_t DecayLengthXYV0() const;

 protected:

  ClassDef(AliAODRecoCascadeHF, 2); // heavy-flavour cascade class
};

#endif
