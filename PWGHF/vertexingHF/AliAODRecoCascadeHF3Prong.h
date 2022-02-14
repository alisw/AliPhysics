#ifndef ALIAODRECOCASCADEHF3PRONG_H
#define ALIAODRECOCASCADEHF3PRONG_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

///***********************************************************
/// \class Class AliAODRecoCascadeHF3Prong
/// \brief base class for AOD reconstructed heavy-flavour cascade 3Prong decays (Xic+->pi Xi pi, ...)
/// The convention is: prong 0 is bachelor, prong 1 is cascade
/// prong 2 is bachelor
///
/// \author Author: Y.S. Watanabe, wyosuke@cns.s.u-tokyo.ac.jp
///***********************************************************

#include <TRef.h>
#include <TRefArray.h>
#include <TClonesArray.h>
#include <TClass.h>
#include "AliAODVertex.h"
#include "AliAODcascade.h"
#include "AliAODv0.h"
#include "AliAODRecoDecayHF3Prong.h"

class AliAODRecoCascadeHF3Prong : public AliAODRecoDecayHF3Prong {

 public:

  AliAODRecoCascadeHF3Prong();
  AliAODRecoCascadeHF3Prong(AliAODVertex *vtx2, Short_t charge,
			    Double_t *px, Double_t *py, Double_t *pz,
			    Double_t *d0, Double_t *d0err, 
			    Double_t *dca, Double_t sigvert,
			    Double_t dist12,Double_t dist23);
  virtual ~AliAODRecoCascadeHF3Prong();

  AliAODRecoCascadeHF3Prong(const AliAODRecoCascadeHF3Prong& source);
  AliAODRecoCascadeHF3Prong& operator=(const AliAODRecoCascadeHF3Prong& source);

  AliAODTrack* GetBachelor1() const {return (AliAODTrack*)GetDaughter(0);}
  AliAODTrack* GetBachelor2() const {return (AliAODTrack*)GetDaughter(2);}
  AliAODcascade* GetCascade() const {
    if ( ! ((AliAODRecoDecay*)GetDaughter(1))->IsA()->InheritsFrom("AliAODcascade") ){
      AliWarning("Object is not of type cascade");
      return 0;
    }
    return (AliAODcascade*)GetDaughter(1);
  }

  AliAODTrack* GetCascadePositiveTrack() const { return  (AliAODTrack*)GetCascade()->GetDaughter(0);  }
  AliAODTrack* GetCascadeNegativeTrack() const { return  (AliAODTrack*)GetCascade()->GetDaughter(1);  }
  AliAODTrack* GetCascadeBachelorTrack() const { return  (AliAODTrack*)GetCascade()->GetDecayVertexXi()->GetDaughter(0);  }

  /// Xic invariant mass
  Double_t InvMassPiXiPi() const {
    UInt_t pdg[3]={211,3312,211}; return InvMass(3,pdg);
  }
  Double_t YXicPlus() const {return Y(4232);}
  
  Int_t MatchToMC(Int_t pdgabs,Int_t pdgabscasc, Int_t *pdgDg,Int_t *pdgDgcasc, Int_t *pdgDgv0,TClonesArray *mcArray) const;
  Int_t MatchToMCCascade(AliAODcascade *casc, Int_t pdgabscasc, Int_t *pdgDgcasc,Int_t *pdgDgv0, TClonesArray *mcArray) const;
	Int_t MatchToMCXicPlus(Int_t pdgabs,TClonesArray *mcArray,Int_t dgLabels[10],Int_t ndg,Int_t ndgCk, const Int_t *pdgDg) const;

  Double_t CascDcaXiDaughters() const;
  Double_t CascDcaV0Daughters() const;
  Double_t CascDecayLength() const;
  Double_t CascDecayLengthV0() const;
  Double_t CascCosPointingAngle() const;
  Double_t CascCosPointingAngleV0() const;
  Double_t CascDcaV0ToPrimVertex() const;
  Double_t CascDcaPosToPrimVertex() const;
  Double_t CascDcaNegToPrimVertex() const;
  Double_t CascDcaBachToPrimVertex() const;
  Double_t CascMassXi() const;
  Double_t CascMassLambda() const;
  Double_t CascMassAntiLambda() const;

  Double_t XicCosPointingAngle() const;
  Double_t BachelorsCosPointingAngle() const;

 protected:

  /// \cond CLASSIMP
  ClassDef(AliAODRecoCascadeHF3Prong, 2); /// heavy-flavour cascade 3prong class
  /// \endcond
};

#endif
