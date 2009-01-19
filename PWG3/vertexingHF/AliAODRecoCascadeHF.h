#ifndef ALIAODRECOCASCADEHF_H
#define ALIAODRECOCASCADEHF_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
#include "AliAODVertex.h"
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
  AliAODRecoDecayHF2Prong* Get2Prong() const {return (AliAODRecoDecayHF2Prong*)GetDaughter(1);}

  // Bachelor (soft pion for Dstar)
  AliAODTrack* GetBachelor() const {return (AliAODTrack*)GetDaughter(0);}

  // D*->D0pi, D0->Kpi
  Double_t EDstar() const {return E(413);} 
  Double_t YDstar() const {return Y(413);} 
  Bool_t   SelectDstar(const Double_t *cutsDstar,const Double_t *cutsD0,Bool_t testD0=kTRUE) const;
  Double_t InvMassD0() const {return (Charge()>0 ? Get2Prong()->InvMassD0() : Get2Prong()->InvMassD0bar());}
  Double_t InvMassDstarKpipi() const;
  Double_t DeltaInvMass() const {return (InvMassDstarKpipi()-InvMassD0());}

 protected:

  ClassDef(AliAODRecoCascadeHF, 2); // heavy-flavour cascade class
};

#endif
