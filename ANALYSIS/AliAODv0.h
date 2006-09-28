#ifndef ALIAODV0_H
#define ALIAODV0_H

/* Copyright(c) 2004-2005, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) V0 vertex class
//     Origin: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//     Purpose: Having observables for physics available for V0s
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TDatabasePDG.h>

#define MASS(PID)  TDatabasePDG::Instance()->GetParticle((PID))->Mass()
#define MASS2(PID) MASS((PID))*MASS((PID))

class AliESD;
class AliESDVertex;
class AliESDv0;
class AliESDtrack;


class AliAODv0 : public TObject {

public:
  AliAODv0();
  AliAODv0(AliESDv0 *rV0Vertex, AliESD *rEvent);
  AliAODv0(const AliAODv0& rAliAODv0);
  virtual ~AliAODv0();

  AliAODv0& operator=(const AliAODv0& rAliAODv0);

  void     Fill(AliESDv0 *rV0Vertex, AliESD *rEvent);
  void     ResetV0();

  Double_t DecayVertexV0X() const;
  Double_t DecayVertexV0Y() const;
  Double_t DecayVertexV0Z() const;

  Double_t DecayLengthV0(double*) const;
                     
  Double_t DcaV0Daughters() const;
  Double_t DcaV0ToPrimVertex() const;
  Double_t DcaPosToPrimVertex() const; 
  Double_t DcaNegToPrimVertex() const; 
  Double_t CosPointAngle(Double_t&, Double_t&, Double_t&) const;
  Double_t RadiusV0()           const;
  Double_t OpenAngleV0()        const;

  Double_t MomPosX() const;
  Double_t MomPosY() const;
  Double_t MomPosZ() const;
  Double_t MomNegX() const;
  Double_t MomNegY() const;
  Double_t MomNegZ() const;

  Int_t    KeyPos()  const;
  Int_t    KeyNeg()  const;

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


protected:
  Double_t fDecayVertexV0X;       // decay vertex of V0 along X
  Double_t fDecayVertexV0Y;       // decay vertex of V0 along Y
  Double_t fDecayVertexV0Z;       // decay vertex of V0 along Z
  Double_t fDcaV0Daughters;       // dca between V0 daughters
  Double_t fDcaV0ToPrimVertex;    // dca of V0 to primary vertex 
  Double_t fDcaPosToPrimVertex;   // dca of pos daughter to primary vertex 
  Double_t fDcaNegToPrimVertex;   // dca of pos daughter to primary vertex 
  Double_t fMomPosX;              // momemtum of pos daughter along X
  Double_t fMomPosY;              // momemtum of pos daughter along Y
  Double_t fMomPosZ;              // momemtum of pos daughter along Z
  Double_t fMomNegX;              // momemtum of neg daughter along X
  Double_t fMomNegY;              // momemtum of neg daughter along Y
  Double_t fMomNegZ;              // momemtum of neg daughter along Z

  UInt_t   fKeyPos;               // track key/index to pos daughter 
  UInt_t   fKeyNeg;               // track key/index to neg daughter 

  Double_t fChi2;                 // main quality variable of V0

  ClassDef(AliAODv0,1)
};

inline Double_t AliAODv0::DecayVertexV0X() const {return fDecayVertexV0X;}
inline Double_t AliAODv0::DecayVertexV0Y() const {return fDecayVertexV0Y;}
inline Double_t AliAODv0::DecayVertexV0Z() const {return fDecayVertexV0Z;}

inline Double_t AliAODv0::DecayLengthV0(double *tParentVertexPosition) const {
  return ::sqrt(::pow(DecayVertexV0X() - tParentVertexPosition[0],2) +
		::pow(DecayVertexV0Y() - tParentVertexPosition[1],2) +
		::pow(DecayVertexV0Z() - tParentVertexPosition[2],2));
}

inline Double_t AliAODv0::DcaV0Daughters() const {return fDcaV0Daughters;}
inline Double_t AliAODv0::DcaV0ToPrimVertex() const {return fDcaV0ToPrimVertex;}
inline Double_t AliAODv0::DcaPosToPrimVertex() const {return fDcaPosToPrimVertex;}
inline Double_t AliAODv0::DcaNegToPrimVertex() const {return fDcaNegToPrimVertex;}

inline Double_t AliAODv0::MomPosX() const {return fMomPosX;}
inline Double_t AliAODv0::MomPosY() const {return fMomPosY;}
inline Double_t AliAODv0::MomPosZ() const {return fMomPosZ;}
inline Double_t AliAODv0::MomNegX() const {return fMomNegX;}
inline Double_t AliAODv0::MomNegY() const {return fMomNegY;}
inline Double_t AliAODv0::MomNegZ() const {return fMomNegZ;}

inline Int_t AliAODv0::KeyPos() const {return fKeyPos;}
inline Int_t AliAODv0::KeyNeg() const {return fKeyNeg;}

inline Double_t AliAODv0::Chi2V0() const {return fChi2;}

inline Double_t AliAODv0::CosPointAngle(Double_t& refPointX, Double_t& refPointY, Double_t& refPointZ) const {
  
  Double_t deltaPos[3]; //vector between the reference point and the V0 vertex
  deltaPos[0] = fDecayVertexV0X - refPointX;
  deltaPos[1] = fDecayVertexV0Y - refPointY;
  deltaPos[2] = fDecayVertexV0Z - refPointZ;

  Double_t deltaPos2 = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];
  
  Double_t cosinePointingAngle = (deltaPos[0]*MomV0X() +
				  deltaPos[1]*MomV0Y() +
				  deltaPos[2]*MomV0Z() ) /
    TMath::Sqrt(Ptot2V0() * deltaPos2);
  return cosinePointingAngle;
}

inline Double_t AliAODv0::MomV0X() const {return MomPosX()+MomNegX();}
inline Double_t AliAODv0::MomV0Y() const {return MomPosY()+MomNegY();}
inline Double_t AliAODv0::MomV0Z() const {return MomPosZ()+MomNegZ();}

inline Double_t AliAODv0::Ptot2Pos() const {
  return (::pow(MomPosX(),2) + ::pow(MomPosY(),2) + ::pow(MomPosZ(),2) );
}
inline Double_t AliAODv0::Ptot2Neg() const {
  return (::pow(MomNegX(),2) + ::pow(MomNegY(),2) + ::pow(MomNegZ(),2) );
}
inline Double_t AliAODv0::Pt2V0() const {
  return (::pow(MomV0X(),2) + ::pow(MomV0Y(),2) );
}

inline Double_t AliAODv0::Ptot2V0() const {return ( Pt2V0() + ::pow(MomV0Z(),2) );}

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
  return 1.-(2./(1.+(MomPosAlongV0()/MomNegAlongV0())));
}
inline Double_t AliAODv0::PtArmV0() const {
  return ::sqrt(Ptot2Pos()-MomPosAlongV0()*MomPosAlongV0());
}

inline Double_t AliAODv0::EPosProton() const {
  return ::sqrt(Ptot2Pos()+MASS2("proton"));
}

inline Double_t AliAODv0::ENegProton() const {
  return ::sqrt(Ptot2Neg()+MASS2("antiproton"));
}

inline Double_t AliAODv0::EPosPion() const {
  return ::sqrt(Ptot2Pos()+MASS2("pi+"));
}

inline Double_t AliAODv0::ENegPion() const {
  return ::sqrt(Ptot2Neg()+MASS2("pi-"));
}

inline Double_t AliAODv0::ELambda() const {
  return ::sqrt(Ptot2V0()+MASS2("Lambda0"));
}

inline Double_t AliAODv0::EK0Short() const {
  return ::sqrt(Ptot2V0()+MASS2("K_S0"));
}

inline Double_t AliAODv0::MassLambda() const {
  return ::sqrt(::pow(EPosProton()+ENegPion(),2)-Ptot2V0());
}

inline Double_t AliAODv0::MassAntiLambda() const {
  return ::sqrt(::pow(ENegProton()+EPosPion(),2)-Ptot2V0());
}

inline Double_t AliAODv0::MassK0Short() const {
  return ::sqrt(::pow(EPosPion()+ENegPion(),2)-Ptot2V0());
}

inline Double_t AliAODv0::RapK0Short() const {
  Double_t ek0 = EK0Short();
  Double_t mMomV0Z = MomV0Z();
  return 0.5*::log((ek0+mMomV0Z)/(ek0-mMomV0Z));
}

inline Double_t AliAODv0::RapLambda() const {
  Double_t eLambda = ELambda();
  Double_t mMomV0Z = MomV0Z();
  return 0.5*::log((eLambda+mMomV0Z)/(eLambda-mMomV0Z));
}

inline Double_t AliAODv0::RadiusV0() const {
  return ::sqrt( fDecayVertexV0X*fDecayVertexV0X
		 + fDecayVertexV0Y*fDecayVertexV0Y );
}

inline Double_t AliAODv0::PseudoRapV0() const {
  Double_t lTheta = ::acos( MomV0Z()/::sqrt(Ptot2V0()) );
  return ( -::log(TMath::Tan(lTheta/2.)) );
}

inline Double_t AliAODv0::PseudoRapPos()   const {
  Double_t lTheta = ::acos( MomPosZ()/::sqrt(Ptot2Pos()) );
  return ( -::log(TMath::Tan(lTheta/2.)) );
}

inline Double_t AliAODv0::PseudoRapNeg()   const {
  Double_t lTheta = ::acos( MomNegZ()/::sqrt(Ptot2Neg()) );
  return ( -::log(TMath::Tan(lTheta/2.)) );
}

inline Double_t AliAODv0::OpenAngleV0() const {
  Double_t lPtot1xPtot2 = fMomPosX*fMomNegX+fMomPosY*fMomNegY+fMomPosZ*fMomNegZ;
  Double_t lPtot1Ptot2_2 = Ptot2Pos()*Ptot2Neg();
  return ::acos(lPtot1xPtot2/::sqrt(lPtot1Ptot2_2) );
}

#endif
