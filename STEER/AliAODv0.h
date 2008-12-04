#ifndef ALIAODV0_H
#define ALIAODV0_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Analysis Oriented Data (AOD) V0 vertex class
//     Authors: B.Hippolyte, IPHC, hippolyt@in2p3.fr 
//              G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//-------------------------------------------------------------------------

#include "AliAODRecoDecay.h"

class AliAODv0 : public AliAODRecoDecay {

public:

  AliAODv0();
  AliAODv0(AliAODVertex *rAODVertex, Double_t rDcaV0Daughters, Double_t rDcaV0ToPrimVertex,
	   const Double_t *rMomPos, const Double_t *rMomNeg, Double_t *rDcaDaughterToPrimVertex);
  virtual ~AliAODv0();

  AliAODv0(const AliAODv0& rAliAODv0);
  AliAODv0& operator=(const AliAODv0& rAliAODv0);

  void     Fill(AliAODVertex *rAODVertex, Double_t rDcaV0Daughters, Double_t rDcaV0ToPrimVertex,
		const Double_t *rMomPos, const Double_t *rMomNeg, const Double_t *rDcaDaughterToPrimVertex);
  void     ResetV0();
  void     Print(Option_t* option = "") const;

  void     SetOnFlyStatus(Bool_t status){fOnFlyStatus=status;}
  Bool_t   GetOnFlyStatus() const {return fOnFlyStatus;}

  Double_t DecayVertexV0X() const;
  Double_t DecayVertexV0Y() const;
  Double_t DecayVertexV0Z() const;

  Double_t DecayLengthV0(const Double_t *) const;
                     
  Double_t DcaV0Daughters()     const;
  Double_t DcaV0ToPrimVertex()  const;
  Double_t DcaPosToPrimVertex() const; 
  Double_t DcaNegToPrimVertex() const; 
  Double_t RadiusV0()           const;
  Double_t OpenAngleV0()        const;

  Double_t MomPosX() const;
  Double_t MomPosY() const;
  Double_t MomPosZ() const;
  Double_t MomNegX() const;
  Double_t MomNegY() const;
  Double_t MomNegZ() const;

  Double_t Chi2V0()  const;

  Double_t MomV0X()  const;
  Double_t MomV0Y()  const;
  Double_t MomV0Z()  const;

  Double_t Ptot2Pos() const;
  Double_t Ptot2Neg() const;
  Double_t Ptot2V0()  const;
  Double_t Pt2V0()    const;
  Double_t MomPosAlongV0() const;
  Double_t MomNegAlongV0() const;
  Double_t AlphaV0() const;
  Double_t PtArmV0() const;
  Double_t EPosProton() const;
  Double_t ENegProton() const;
  Double_t EPosPion() const;
  Double_t ENegPion() const;
  Double_t ELambda() const;
  Double_t EK0Short() const;
  Double_t MassLambda() const;
  Double_t MassAntiLambda() const;
  Double_t MassK0Short() const;
  Double_t RapK0Short() const;
  Double_t RapLambda() const;
  Double_t PseudoRapV0()    const;
  Double_t PseudoRapPos()   const;
  Double_t PseudoRapNeg()   const;
  
  Short_t  GetPosID()       const;
  Short_t  GetNegID()       const;
  Int_t    GetLabel() const {return -1;} // Dummy

protected:
  Double32_t fDcaV0ToPrimVertex;    // dca of V0 to primary vertex 
  Bool_t     fOnFlyStatus;          // if kTRUE, then this V0 is recontructed
                                    // "on fly" during the tracking
  ClassDef(AliAODv0,2)
};

inline Double_t AliAODv0::DecayVertexV0X() const {return this->GetSecVtxX();}
inline Double_t AliAODv0::DecayVertexV0Y() const {return this->GetSecVtxY();}
inline Double_t AliAODv0::DecayVertexV0Z() const {return this->GetSecVtxZ();}

inline Double_t AliAODv0::DecayLengthV0(const Double_t *tParentVertexPosition) const {
  return ::sqrt(::pow(DecayVertexV0X() - tParentVertexPosition[0],2) +
		::pow(DecayVertexV0Y() - tParentVertexPosition[1],2) +
		::pow(DecayVertexV0Z() - tParentVertexPosition[2],2));
}

inline Double_t AliAODv0::DcaV0Daughters() const {return fDCA[0];}
inline Double_t AliAODv0::DcaV0ToPrimVertex() const {return fDcaV0ToPrimVertex;}
inline Double_t AliAODv0::DcaPosToPrimVertex() const {return fd0[0];}
inline Double_t AliAODv0::DcaNegToPrimVertex() const {return fd0[1];}

inline Double_t AliAODv0::RadiusV0() const {
  return RadiusSecVtx();
}

inline Double_t AliAODv0::OpenAngleV0() const {
  Double_t lScalPtot1Ptot2 = PxProng(0)*PxProng(1)+PyProng(0)*PyProng(1)+PzProng(0)*PzProng(1);
  Double_t lPtot1xPtot2 = Ptot2Pos()*Ptot2Neg();
  return ::acos(lScalPtot1Ptot2/::sqrt(lPtot1xPtot2) );
}

inline Double_t AliAODv0::MomPosX() const {return fPx[0];}
inline Double_t AliAODv0::MomPosY() const {return fPy[0];}
inline Double_t AliAODv0::MomPosZ() const {return fPz[0];}
inline Double_t AliAODv0::MomNegX() const {return fPx[1];}
inline Double_t AliAODv0::MomNegY() const {return fPy[1];}
inline Double_t AliAODv0::MomNegZ() const {return fPz[1];}

inline Double_t AliAODv0::Chi2V0() const {return GetSecondaryVtx()->GetChi2perNDF();}

// Compare eventually AliAODv0::MomV0X() and AliAODRecoDecay::Px()
inline Double_t AliAODv0::MomV0X() const {return MomPosX()+MomNegX();}
inline Double_t AliAODv0::MomV0Y() const {return MomPosY()+MomNegY();}
inline Double_t AliAODv0::MomV0Z() const {return MomPosZ()+MomNegZ();}

inline Double_t AliAODv0::Ptot2Pos() const {
  return (::pow(MomPosX(),2) + ::pow(MomPosY(),2) + ::pow(MomPosZ(),2) );
}
inline Double_t AliAODv0::Ptot2Neg() const {
  return (::pow(MomNegX(),2) + ::pow(MomNegY(),2) + ::pow(MomNegZ(),2) );
}
inline Double_t AliAODv0::Ptot2V0() const {return ( Pt2V0() + ::pow(MomV0Z(),2) );}
inline Double_t AliAODv0::Pt2V0() const {
  return (::pow(MomV0X(),2) + ::pow(MomV0Y(),2) );
}

inline Double_t AliAODv0::MomPosAlongV0() const {
  Double_t lPtot2V0 = Ptot2V0();
  if (lPtot2V0)
    return (MomPosX()*MomV0X() +
	    MomPosY()*MomV0Y() +
	    MomPosZ()*MomV0Z()) / ::sqrt(lPtot2V0);
  return 0.;
}

inline Double_t AliAODv0::MomNegAlongV0() const {
  Double_t lPtot2V0 = Ptot2V0();
  if (lPtot2V0)
    return (MomNegX()*MomV0X() +
	    MomNegY()*MomV0Y() +
	    MomNegZ()*MomV0Z()) / ::sqrt(lPtot2V0);
  return 0.;
}

inline Double_t AliAODv0::AlphaV0() const {
  return Alpha();
}
inline Double_t AliAODv0::PtArmV0() const {
  return QtProng(0);
}

inline Double_t AliAODv0::EPosProton() const {
  return EProng(0,2212);
}

inline Double_t AliAODv0::ENegProton() const {
  return EProng(1,2212);
}

inline Double_t AliAODv0::EPosPion() const {
  return EProng(0,211);
}

inline Double_t AliAODv0::ENegPion() const {
  return EProng(1,211);
}

inline Double_t AliAODv0::ELambda() const {
  return E(3122);
}

inline Double_t AliAODv0::EK0Short() const {
  return E(310);
}

inline Double_t AliAODv0::MassLambda() const {
  return InvMass2Prongs(0,1,2212,211);
}

inline Double_t AliAODv0::MassAntiLambda() const {
  return InvMass2Prongs(0,1,211,2212);
}

inline Double_t AliAODv0::MassK0Short() const {
  return InvMass2Prongs(0,1,211,211);
}

inline Double_t AliAODv0::RapK0Short() const {
  return Y(310);
}

inline Double_t AliAODv0::RapLambda() const {
  return Y(3122);
}

inline Double_t AliAODv0::PseudoRapV0() const {
  return Eta();
}

inline Double_t AliAODv0::PseudoRapPos()   const {
  return EtaProng(0);
}

inline Double_t AliAODv0::PseudoRapNeg()   const {
  return EtaProng(1);
}
//----------------------------------------------------------------------------

#endif
