#ifndef ALIAODV0_H
#define ALIAODV0_H

/* Copyright(c) 2004-2005, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) V0 vertex class
//
//     Origin: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//
//     Purpose: Having observables for physics available for V0s
//-------------------------------------------------------------------------

#include <TObject.h>

class AliESD;
class AliESDVertex;
class AliESDv0;
class AliESDtrack;


class AliAODv0 : public TObject {

public:
  AliAODv0();
  AliAODv0(AliESDv0 *rV0Vertex, AliESD *rEvent);
  void     Fill(AliESDv0 *rV0Vertex, AliESD *rEvent);
  void     ResetV0();

  Double_t DecayVertexV0X() const;
  Double_t DecayVertexV0Y() const;
  Double_t DecayVertexV0Z() const;
                     
  Double_t DcaV0Daughters() const;
  Double_t DcaV0ToPrimVertex() const;
  Double_t DcaPosToPrimVertex() const; 
  Double_t DcaNegToPrimVertex() const; 

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

  Int_t   fKeyPos;                // track key/index to pos daughter 
  Int_t   fKeyNeg;                // track key/index to neg daughter 

  Double_t fChi2;                 // main quality variable of V0
  AliESD  *fEvent;                // pointer to current event

  ClassDef(AliAODv0,1)
};

inline Double_t AliAODv0::DecayVertexV0X() const {return fDecayVertexV0X;}
inline Double_t AliAODv0::DecayVertexV0Y() const {return fDecayVertexV0Y;}
inline Double_t AliAODv0::DecayVertexV0Z() const {return fDecayVertexV0Z;}

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

#endif


