#ifndef ALIAODXI_H
#define ALIAODXI_H

/* Copyright(c) 2004-2005, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) Xi vertex class
//     Origin: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//     Purpose: Having observables for physics available for Xis
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TDatabasePDG.h>

#define MASS(PID)  TDatabasePDG::Instance()->GetParticle((PID))->Mass()
#define MASS2(PID) MASS((PID))*MASS((PID))

class AliESD;
class AliESDVertex;
class AliESDv0;
class AliESDcascade;
class AliESDtrack;

class AliAODxi : public AliAODv0 {

public:
  AliAODxi();
  AliAODxi(AliESDcascade *rXiVertex, AliESD *rEvent);
  AliAODxi(const AliAODxi& rAliAODxi);
  virtual ~AliAODxi();

  AliAODxi& operator=(const AliAODxi& rAliAODxi);

  void     Fill(AliESDcascade *rXiVertex, AliESD *rEvent);
  void     ResetXi();


  Int_t    Charge()                  const;
  Double_t DecayVertexXiX()          const;
  Double_t DecayVertexXiY()          const;
  Double_t DecayVertexXiZ()          const;
  virtual Double_t DecayLengthV0()   const;
  Double_t DcaXiDaughters()          const;
  Double_t DcaBachelorToPrimVertex() const;
  Double_t DcaV0ToPrimVertex(Double_t&, Double_t&, Double_t&) const;
  Double_t DcaXiToPrimVertex()       const;
  Double_t CosPointAngle(Double_t&, Double_t&, Double_t&) const;

  Double_t DecayLengthXi(Double_t*)  const;
                     
  Double_t MomBachelorX()       const;
  Double_t MomBachelorY()       const;
  Double_t MomBachelorZ()       const;
  UInt_t   KeyBachelor()        const;
  Double_t Chi2Xi()             const;
  Double_t MomXiX()             const;
  Double_t MomXiY()             const;
  Double_t MomXiZ()             const;

  Double_t Ptot2Bachelor()      const;
  Double_t Ptot2Xi()            const;
  Double_t Pt2Xi()              const;
  Double_t MomBachelorAlongXi() const;
  Double_t MomV0AlongXi()       const;
  Double_t AlphaXi()            const;
  Double_t PtArmXi()            const;
  Double_t EBachelorPion()      const;
  Double_t EBachelorKaon()      const;
  Double_t EXi()                const;
  Double_t EOmega()             const;
  Double_t MassXi()             const;
  Double_t MassOmega()          const;
  Double_t RapXi()              const;
  Double_t RapOmega()           const;


protected:

  Int_t    fCharge;                  // charge of Xi
  Double_t fDecayVertexXiX;          // decay vertex of Xi along X
  Double_t fDecayVertexXiY;          // decay vertex of Xi along Y
  Double_t fDecayVertexXiZ;          // decay vertex of Xi along Z
  Double_t fDcaXiDaughters;          // dca between Xi daughters
  Double_t fDcaXiToPrimVertex;       // dca of Xi to primary vertex 
  Double_t fDcaBachelorToPrimVertex; // dca of bachelor to primary vertex 
  Double_t fMomBachelorX;            // momemtum of bachelor along X
  Double_t fMomBachelorY;            // momemtum of bachelor along Y
  Double_t fMomBachelorZ;            // momemtum of bachelor along Z

  UInt_t   fKeyBachelor;             // track key/index to bachelor 

  Double_t fChi2Xi;                  // main quality variable of Xi

  ClassDef(AliAODxi,1)
};

inline Int_t    AliAODxi::Charge() const {return fCharge;}

inline Double_t AliAODxi::DecayVertexXiX() const {return fDecayVertexXiX;}
inline Double_t AliAODxi::DecayVertexXiY() const {return fDecayVertexXiY;}
inline Double_t AliAODxi::DecayVertexXiZ() const {return fDecayVertexXiZ;}

inline Double_t AliAODxi::DecayLengthV0() const {
    return ::sqrt(::pow(fDecayVertexV0X - fDecayVertexXiX,2) +
		  ::pow(fDecayVertexV0Y - fDecayVertexXiY,2) +
		  ::pow(fDecayVertexV0Z - fDecayVertexXiZ,2));
}

inline Double_t AliAODxi::DcaV0ToPrimVertex(Double_t& primVertexX, Double_t& primVertexY, Double_t& primVertexZ) const {
  Double_t momV0X=MomV0X();
  Double_t momV0Y=MomV0Y();
  Double_t momV0Z=MomV0Z();
  Double_t dx=(primVertexY-fDecayVertexV0Y)*momV0Z - (primVertexZ-fDecayVertexV0Z)*momV0Y; 
  Double_t dy=(primVertexX-fDecayVertexV0X)*momV0Z - (primVertexZ-fDecayVertexV0Z)*momV0X;
  Double_t dz=(primVertexX-fDecayVertexV0X)*momV0Y - (primVertexY-fDecayVertexV0Y)*momV0X;
  return TMath::Sqrt((dx*dx+dy*dy+dz*dz)/(momV0X*momV0X+momV0Y*momV0Y+momV0Z*momV0Z));
}

inline Double_t AliAODxi::DcaXiDaughters() const {return fDcaXiDaughters;}
inline Double_t AliAODxi::DcaBachelorToPrimVertex() const {return fDcaBachelorToPrimVertex;}
inline Double_t AliAODxi::DcaXiToPrimVertex() const {return fDcaXiToPrimVertex;}

inline Double_t AliAODxi::CosPointAngle(Double_t& refPointX, Double_t& refPointY, Double_t& refPointZ) const {
  
  Double_t deltaPos[3]; //vector between the reference point and the cascade vertex
  deltaPos[0] = fDecayVertexXiX - refPointX;
  deltaPos[1] = fDecayVertexXiY - refPointY;
  deltaPos[2] = fDecayVertexXiZ - refPointZ;

  Double_t deltaPos2 = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];
  
  Double_t cosinePointingAngle = (deltaPos[0]*MomXiX() +
				  deltaPos[1]*MomXiY() +
				  deltaPos[2]*MomXiZ() ) /
    TMath::Sqrt(Ptot2Xi() * deltaPos2);
  return cosinePointingAngle;
}

inline Double_t AliAODxi::DecayLengthXi(Double_t *tPrimaryVertexPosition) const {
  return ::sqrt(::pow(fDecayVertexXiX - tPrimaryVertexPosition[0],2) +
		::pow(fDecayVertexXiY - tPrimaryVertexPosition[1],2) +
		::pow(fDecayVertexXiZ - tPrimaryVertexPosition[2],2));
}

inline Double_t AliAODxi::MomBachelorX() const {return fMomBachelorX;}
inline Double_t AliAODxi::MomBachelorY() const {return fMomBachelorY;}
inline Double_t AliAODxi::MomBachelorZ() const {return fMomBachelorZ;}
inline UInt_t   AliAODxi::KeyBachelor() const {return fKeyBachelor;}
inline Double_t AliAODxi::Chi2Xi() const {return fChi2Xi;}
inline Double_t AliAODxi::MomXiX() const {return MomV0X()+fMomBachelorX;}
inline Double_t AliAODxi::MomXiY() const {return MomV0Y()+fMomBachelorY;}
inline Double_t AliAODxi::MomXiZ() const {return MomV0Z()+fMomBachelorZ;}

inline Double_t AliAODxi::Ptot2Bachelor() const {
  return (::pow(fMomBachelorX,2) + ::pow(fMomBachelorY,2) + ::pow(fMomBachelorZ,2) );
}
inline Double_t AliAODxi::Ptot2Xi() const {return ( Pt2Xi() + ::pow(MomXiZ(),2) );}
inline Double_t AliAODxi::Pt2Xi() const {
  return (::pow(MomXiX(),2) + ::pow(MomXiY(),2) );
}

inline Double_t AliAODxi::MomBachelorAlongXi() const {
  Double_t lPtot2Xi = Ptot2Xi();
  if (lPtot2Xi)
    return (MomBachelorX()*MomXiX() +
	    MomBachelorY()*MomXiY() +
	    MomBachelorZ()*MomXiZ()) / ::sqrt(lPtot2Xi);
  return 0.;
}

inline Double_t AliAODxi::MomV0AlongXi() const {
  Double_t lPtot2Xi = Ptot2Xi();
  if (lPtot2Xi)
    return (MomV0X()*MomXiX() +
	    MomV0Y()*MomXiY() +
	    MomV0Z()*MomXiZ()) / ::sqrt(lPtot2Xi);
  return 0.;
}

inline Double_t AliAODxi::AlphaXi() const {
  Double_t lMomV0AlongXi       = MomV0AlongXi();
  Double_t lMomBachelorAlongXi = MomBachelorAlongXi();

  return (((Float_t) Charge()) * (lMomBachelorAlongXi-lMomV0AlongXi)/
                                 (lMomBachelorAlongXi+lMomV0AlongXi));
}

inline Double_t AliAODxi::PtArmXi() const {
  return ::sqrt(Ptot2V0()-MomBachelorAlongXi()*MomBachelorAlongXi());
}

inline Double_t AliAODxi::EBachelorPion() const {
  return ::sqrt(Ptot2Bachelor()+MASS2("pi-"));
}

inline Double_t AliAODxi::EBachelorKaon() const {
  return ::sqrt(Ptot2Bachelor()+MASS2("K-"));
}

inline Double_t AliAODxi::EXi() const {
  return ::sqrt(Ptot2Xi()+MASS2("Xi-"));
}

inline Double_t AliAODxi::EOmega() const {
  return ::sqrt(Ptot2Xi()+MASS2("Omega-"));
}

inline Double_t AliAODxi::MassXi() const {
  return ::sqrt(::pow(ELambda()+EBachelorPion(),2)-Ptot2Xi());
}

inline Double_t AliAODxi::MassOmega() const {
  return ::sqrt(::pow(ELambda()+EBachelorKaon(),2)-Ptot2Xi());
}

inline Double_t AliAODxi::RapXi() const {
  Double_t exi = EXi();
  Double_t lMomXiZ = MomXiZ();
  return 0.5*::log((exi+lMomXiZ)/(exi-lMomXiZ));
}

inline Double_t AliAODxi::RapOmega() const {
  Double_t eom = EOmega();
  Double_t lMomXiZ = MomXiZ();
  return 0.5*::log((eom+lMomXiZ)/(eom-lMomXiZ));
}

#endif
