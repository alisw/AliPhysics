#ifndef ALIAODV0_H
#define ALIAODV0_H

/* Copyright(c) 2004-2005, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) V0 vertex class
//
//     Origin: B.Hippolyte, IReS, hippolyt@in2p3.fr 
//             G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TPDGCode.h>

class AliESD;
class AliESDVertex;
class AliESDv0;
class AliESDtrack;


class AliAODv0 : public TObject {

public:
  AliAODv0();
  AliAODv0(AliESDv0*,AliESD*);
  void  Fill(AliESDv0*,AliESD*);
  void  ResetV0();

  Double_t decayVertexV0X() const;
  Double_t decayVertexV0Y() const;
  Double_t decayVertexV0Z() const;
                     
  Double_t dcaV0Daughters() const;
  Double_t dcaV0ToPrimVertex() const;
  Double_t dcaPosToPrimVertex() const; 
  Double_t dcaNegToPrimVertex() const; 

  Double_t momPosX() const;
  Double_t momPosY() const;
  Double_t momPosZ() const;
  Double_t momNegX() const;
  Double_t momNegY() const;
  Double_t momNegZ() const;

  Int_t keyPos() const;
  Int_t keyNeg() const;

  Double_t chi2V0() const;

  // Following Need to be moved to Base Class
  Double_t momV0X();
  Double_t momV0Y();
  Double_t momV0Z();

  Double_t Ptot2Pos();
  Double_t Ptot2Neg();
  Double_t Ptot2V0();
  Double_t Pt2V0();
  Double_t MomPosAlongV0();
  Double_t MomNegAlongV0();
  Double_t alphaV0();
  Double_t ptArmV0();
  // Above  Need to be moved to Base Class


protected:
  Double_t fDecayVertexV0X;
  Double_t fDecayVertexV0Y;
  Double_t fDecayVertexV0Z;
  Double_t fDcaV0Daughters;
  Double_t fDcaV0ToPrimVertex;
  Double_t fDcaPosToPrimVertex;
  Double_t fDcaNegToPrimVertex;
  Double_t fMomPosX;
  Double_t fMomPosY;
  Double_t fMomPosZ;
  Double_t fMomNegX;
  Double_t fMomNegY;
  Double_t fMomNegZ;

  Int_t   fKeyPos;
  Int_t   fKeyNeg;

  Double_t fChi2;
  AliESD  *fEvent;

  ClassDef(AliAODv0,1)    // AOD V0 vertex
};

inline Double_t AliAODv0::decayVertexV0X() const {return fDecayVertexV0X;}
inline Double_t AliAODv0::decayVertexV0Y() const {return fDecayVertexV0Y;}
inline Double_t AliAODv0::decayVertexV0Z() const {return fDecayVertexV0Z;}

inline Double_t AliAODv0::dcaV0Daughters() const {return fDcaV0Daughters;}
inline Double_t AliAODv0::dcaV0ToPrimVertex() const {return fDcaV0ToPrimVertex;}
inline Double_t AliAODv0::dcaPosToPrimVertex() const {return fDcaPosToPrimVertex;}
inline Double_t AliAODv0::dcaNegToPrimVertex() const {return fDcaNegToPrimVertex;}

inline Double_t AliAODv0::momPosX() const {return fMomPosX;}
inline Double_t AliAODv0::momPosY() const {return fMomPosY;}
inline Double_t AliAODv0::momPosZ() const {return fMomPosZ;}
inline Double_t AliAODv0::momNegX() const {return fMomNegX;}
inline Double_t AliAODv0::momNegY() const {return fMomNegY;}
inline Double_t AliAODv0::momNegZ() const {return fMomNegZ;}

inline Int_t AliAODv0::keyPos() const {return fKeyPos;}
inline Int_t AliAODv0::keyNeg() const {return fKeyNeg;}

inline Double_t AliAODv0::chi2V0() const {return fChi2;}

// Following Need to be moved to Base Class
inline Double_t AliAODv0::momV0X(){return momPosX()+momNegX();}
inline Double_t AliAODv0::momV0Y(){return momPosY()+momNegY();}
inline Double_t AliAODv0::momV0Z(){return momPosZ()+momNegZ();}

inline Double_t AliAODv0::Ptot2Pos(){
  return (::pow(momPosX(),2) + ::pow(momPosY(),2) + ::pow(momPosZ(),2) );
}
inline Double_t AliAODv0::Ptot2Neg(){
  return (::pow(momNegX(),2) + ::pow(momNegY(),2) + ::pow(momNegZ(),2) );
}
inline Double_t AliAODv0::Pt2V0(){
  return (::pow(momV0X(),2) + ::pow(momV0Y(),2) );
}

inline Double_t AliAODv0::Ptot2V0(){return ( Pt2V0() + ::pow(momV0Z(),2) );}

inline Double_t AliAODv0::MomPosAlongV0(){
  Double_t mPtot2V0 = Ptot2V0();
  if (mPtot2V0)
    return (momPosX()*momV0X() +
	    momPosY()*momV0Y() +
	    momPosZ()*momV0Z()) / ::sqrt(mPtot2V0);
  return 0.;
}

inline Double_t AliAODv0::MomNegAlongV0(){
  Double_t mPtot2V0 = Ptot2V0();
  if (mPtot2V0)
    return (momNegX()*momV0X() +
	    momNegY()*momV0Y() +
	    momNegZ()*momV0Z()) / ::sqrt(mPtot2V0);
  return 0.;
}

inline Double_t AliAODv0::alphaV0(){
  return 1.-(2./(1.+(MomPosAlongV0()/MomNegAlongV0())));
}
inline Double_t AliAODv0::ptArmV0(){
  return ::sqrt(Ptot2Pos()-MomPosAlongV0()*MomPosAlongV0());
}

// Above Need to be moved to Base Class

#endif


