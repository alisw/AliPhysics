#ifndef ALIAODCASCADE_H
#define ALIAODCASCADE_H

/* Copyright(c) 2004-2005, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) Xi vertex class
//     Origin: A.Maire, IReS, antonin.maire@ires.in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//
//     Purpose: Having physics observables available for Xis
//-------------------------------------------------------------------------


#include "AliAODv0.h"

class AliAODTrack;
class AliAODVertex;

class AliAODcascade : public AliAODv0 {

public:
  AliAODcascade();
  
  AliAODcascade(const AliAODcascade& rSource);
  
  AliAODcascade( AliAODVertex* rAODVertexXi,  // No "const" param, constructor "TRef(const TObject*)" doesn't exist.
                      Int_t         rChargeXi,
		      Double_t      rDcaXiDaughters,
		      Double_t      rDcaXiToPrimVertex,
		      Double_t      rDcaBachToPrimVertex,
		const Double_t*     rMomBach,
		
		 AliAODVertex* rAODVertexV0,  // No "const" param, see above.
		      Double_t      rDcaV0Daughters,    
		      Double_t      rDcaV0ToPrimVertex, 
		const Double_t*     rMomPos,            
		const Double_t*     rMomNeg,
		      Double_t*     rDcaDaughterToPrimVertex ); // const? -> Need agreement at AliAODRecoDecay level

  		
  AliAODcascade( AliAODVertex* rAODVertexXi,  // No "const" param, see above.
                      Int_t         rChargeXi,
		      Double_t      rDcaXiDaughters,
		      Double_t      rDcaXiToPrimVertex,
		      Double_t      rDcaBachToPrimVertex,
		const Double_t*     rMomBach,
		const AliAODv0&     rAODv0 );
		
		
		
  virtual ~AliAODcascade();

  AliAODcascade& operator=(const AliAODcascade& rSource);
  
  void  Fill(AliAODVertex* rAODVertexXi,  // No "const" param, see above.
                      Int_t     rChargeXi,
		      Double_t  rDcaXiDaughters,
		      Double_t  rDcaXiToPrimVertex,
		      Double_t  rDcaBachToPrimVertex,
		const Double_t* rMomBach,
			
		AliAODVertex*   rAODVertexV0, // No "const" param, see above.
		      Double_t  rDcaV0Daughters,
		      Double_t  rDcaV0ToPrimVertex,
		const Double_t* rMomPos,
		const Double_t* rMomNeg,
		      Double_t* rDcaDaughterToPrimVertex ); // const? -> Need agreement at AliAODRecoDecay level
  
//   void  Fill(   AliAODVertex*   rAODVertexXi, 
//                       Int_t     rChargeXi,
// 		      Double_t  rDcaXiDaughters,
// 		      Double_t  rDcaXiToPrimVertex,
// 		      Double_t  rDcaBachToPrimVertex,
// 		const Double_t* rMomBach,
// 		const AliAODv0& rAODv0  );          // -> To be implemented ...       
		      
		      
  void  ResetXi();                                       
  void  PrintXi(const Double_t& rPrimVtxX, 
                const Double_t& rPrimVtxY, 
                const Double_t& rPrimVtxZ) const; 

// ----
  Int_t    ChargeXi()                const;
  Int_t	   GetBachID()               const;
  Int_t    GetLabel()                const {return -1;}
  
  AliAODVertex* GetDecayVertexXi()   const; 
  Double_t DecayVertexXiX()          const;
  Double_t DecayVertexXiY()          const;
  Double_t DecayVertexXiZ()          const;
  Double_t Chi2Xi()                  const;
    
  Double_t DcaBachToPrimVertex() const;
  Double_t DcaXiDaughters()          const;
  Double_t DcaXiToPrimVertex()       const;
  Double_t DcaXiToPrimVertex(const Double_t& rPrimVtxX, // hopefully, temporary method ...
                             const Double_t& rPrimVtxY, 
                             const Double_t& rPrimVtxZ) const;
  Double_t CosPointingAngleXi(const Double_t& rPrimVtxX, 
                              const Double_t& rPrimVtxY, 
                              const Double_t& rPrimVtxZ) const;  
  
  Double_t DecayLengthV0()           const;
  Double_t DecayLengthXi(const Double_t& rPrimVtxX, 
                         const Double_t& rPrimVtxY, 
                         const Double_t& rPrimVtxZ) const;
                        
  Double_t MomBachX()       const;
  Double_t MomBachY()       const;
  Double_t MomBachZ()       const;
  
  Double_t MomXiX()         const;
  Double_t MomXiY()         const;
  Double_t MomXiZ()         const;

// ---- 
  Double_t Ptot2Bach()      const;
  Double_t Ptot2Xi()        const;
  Double_t Pt2Xi()          const;
  Double_t MomBachAlongXi() const;
  Double_t MomV0AlongXi()   const;
  Double_t AlphaXi()        const;
  Double_t PtArmXi()        const;
  Double_t EBachPion()      const;
  Double_t EBachKaon()      const;
  Double_t EXi()            const;
  Double_t EOmega()         const;
  Double_t MassXi()         const;
  Double_t MassOmega()      const;
  Double_t RapXi()          const;
  Double_t RapOmega()       const;

protected:

  TRef          fDecayVertexXi;           // ref to decay vertex of the cascade (Xi vertex)
  Short_t       fChargeXi;                // charge of Xi
    
  Double32_t    fDcaXiDaughters;          // dca between Xi daughters
  Double32_t    fDcaXiToPrimVertex;       // dca of Xi to primary vertex 
  Double32_t    fDcaBachToPrimVertex; // dca of bachelor to primary vertex 
  
  Double32_t    fMomBachX;            // momemtum of bachelor along X
  Double32_t    fMomBachY;            // momemtum of bachelor along Y
  Double32_t    fMomBachZ;            // momemtum of bachelor along Z
  
  ClassDef(AliAODcascade,1)   
};

//-----------------------------------------------------------

inline Int_t    AliAODcascade::ChargeXi()   const     {return fChargeXi; }

inline AliAODVertex* AliAODcascade::GetDecayVertexXi() const { return  (AliAODVertex*)fDecayVertexXi.GetObject(); }
inline Double_t AliAODcascade::DecayVertexXiX() const {return GetDecayVertexXi()->GetX(); }
inline Double_t AliAODcascade::DecayVertexXiY() const {return GetDecayVertexXi()->GetY(); }
inline Double_t AliAODcascade::DecayVertexXiZ() const {return GetDecayVertexXi()->GetZ(); }

inline Double_t AliAODcascade::Chi2Xi()         const {return GetDecayVertexXi()->GetChi2(); }

inline Double_t AliAODcascade::DcaBachToPrimVertex() const {return fDcaBachToPrimVertex;}
inline Double_t AliAODcascade::DcaXiDaughters()          const {return fDcaXiDaughters;}
inline Double_t AliAODcascade::DcaXiToPrimVertex()       const {return fDcaXiToPrimVertex;}

inline Double_t AliAODcascade::DecayLengthV0() const {
    return ::sqrt(::pow(DecayVertexV0X() - DecayVertexXiX(),2) +
		  ::pow(DecayVertexV0Y() - DecayVertexXiY(),2) +
		  ::pow(DecayVertexV0Z() - DecayVertexXiZ(),2));
}

inline Double_t AliAODcascade::DecayLengthXi(const Double_t& rPrimVtxX, 
                                             const Double_t& rPrimVtxY, 
                                             const Double_t& rPrimVtxZ) const {
  return ::sqrt(::pow(DecayVertexXiX() - rPrimVtxX,2) +
		::pow(DecayVertexXiY() - rPrimVtxY,2) +
		::pow(DecayVertexXiZ() - rPrimVtxZ,2));
}

inline Double_t AliAODcascade::MomBachX() const {return fMomBachX;}
inline Double_t AliAODcascade::MomBachY() const {return fMomBachY;}
inline Double_t AliAODcascade::MomBachZ() const {return fMomBachZ;}

inline Double_t AliAODcascade::MomXiX() const {return MomV0X()+fMomBachX;}
inline Double_t AliAODcascade::MomXiY() const {return MomV0Y()+fMomBachY;}
inline Double_t AliAODcascade::MomXiZ() const {return MomV0Z()+fMomBachZ;}

inline Double_t AliAODcascade::Ptot2Bach() const {
  return (::pow(fMomBachX,2) + ::pow(fMomBachY,2) + ::pow(fMomBachZ,2) );
}
inline Double_t AliAODcascade::Ptot2Xi() const {return ( Pt2Xi() + ::pow(MomXiZ(),2) );}
inline Double_t AliAODcascade::Pt2Xi() const {
  return (::pow(MomXiX(),2) + ::pow(MomXiY(),2) );
}

inline Double_t AliAODcascade::MomBachAlongXi() const {
  Double_t rPtot2Xi = Ptot2Xi();
  if (rPtot2Xi)
    return (MomBachX()*MomXiX() +
	    MomBachY()*MomXiY() +
	    MomBachZ()*MomXiZ()) / ::sqrt(rPtot2Xi);
  else return 0.;
}

inline Double_t AliAODcascade::MomV0AlongXi() const {
  Double_t rPtot2Xi = Ptot2Xi();
  if (rPtot2Xi)
    return (MomV0X()*MomXiX() +
	    MomV0Y()*MomXiY() +
	    MomV0Z()*MomXiZ()) / ::sqrt(rPtot2Xi);
  return 0.;
}

inline Double_t AliAODcascade::AlphaXi() const {
  Double_t rMomV0AlongXi   = MomV0AlongXi();
  Double_t rMomBachAlongXi = MomBachAlongXi();

  return (((Float_t) ChargeXi()) * (rMomBachAlongXi - rMomV0AlongXi)/
                                   (rMomBachAlongXi + rMomV0AlongXi));
}

inline Double_t AliAODcascade::PtArmXi() const {
  return ::sqrt(Ptot2Bach()-MomBachAlongXi()*MomBachAlongXi());
}

inline Double_t AliAODcascade::MassXi() const {
  return ::sqrt(::pow(ELambda()+EBachPion(),2)-Ptot2Xi());
}

inline Double_t AliAODcascade::MassOmega() const {
  return ::sqrt(::pow(ELambda()+EBachKaon(),2)-Ptot2Xi());
}

inline Double_t AliAODcascade::RapXi() const {
  Double_t exi = EXi();
  Double_t rMomXiZ = MomXiZ();
  return 0.5*::log((exi+rMomXiZ)/(exi-rMomXiZ));
}

inline Double_t AliAODcascade::RapOmega() const {
  Double_t eom = EOmega();
  Double_t rMomXiZ = MomXiZ();
  return 0.5*::log((eom+rMomXiZ)/(eom-rMomXiZ));
}

#endif
