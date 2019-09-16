#ifndef ALIAODKINK_H
#define ALIAODKINK_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliAODkink
/// \brief Analysis Oriented Data (AOD) V0 vertex class
///
/// G.Van Buren, BNL,  gene@bnl.gov      (original STAR MuDsts)
///
/// \author B.Hippolyte, IPHC, hippolyt@in2p3.fr

#include "TVector3.h"
#include "AliAODRecoDecay.h"


class AliAODkink : public AliAODRecoDecay {

public:

  AliAODkink();
  AliAODkink(AliAODVertex *rAODVertex, Float_t rkinkAngle, Float_t rradius, Float_t qT, const TVector3& rmotherMfromKink, const TVector3& rdaughterMKink);
  virtual ~AliAODkink();

  AliAODkink(const AliAODkink& rAliAODkink);
  AliAODkink& operator=(const AliAODkink& rAliAODkink);

  void     Fill(AliAODVertex *rAODVertex, Float_t rkinkAngle, Float_t rradius, Float_t qT, const TVector3& rmotherMfromKink, const TVector3& rdaughterMKink);
  void     ResetKink();
  void     Print(Option_t* option = "") const;

  Double_t DecayVertexKinkX() const;
  Double_t DecayVertexKinkY() const;
  Double_t DecayVertexKinkZ() const;

  Double_t RadiusKink() const;  
  Float_t  KinkAngle()const; 
  Float_t  KinkQt() const;

  Double_t MomMotherX() const;
  Double_t MomMotherY() const;
  Double_t MomMotherZ() const;
  Double_t MomDaughterX() const;
  Double_t MomDaughterY() const;
  Double_t MomDaughterZ() const;

  Double_t PseudoRapMKink()    const;
  
  Int_t    GetLabel()       const {return -1;} // Dummy
  Int_t    PdgCode()        const {return  0;} // Dummy
  
  virtual void     SetID(Short_t /*id*/) {;}
  
protected:

  Double_t fRadius;
  Float_t  fkinkAngle;
  Float_t  fQt;
  TVector3 fMotherMfromKink; 
   
  ClassDef(AliAODkink,1)
};

inline Double_t AliAODkink::DecayVertexKinkX() const {return this->GetSecVtxX();}
inline Double_t AliAODkink::DecayVertexKinkY() const {return this->GetSecVtxY();}
inline Double_t AliAODkink::DecayVertexKinkZ() const {return this->GetSecVtxZ();}

inline Double_t AliAODkink::MomMotherX() const {return fMotherMfromKink[0];}
inline Double_t AliAODkink::MomMotherY() const {return fMotherMfromKink[1];}
inline Double_t AliAODkink::MomMotherZ() const {return fMotherMfromKink[2];}
inline Double_t AliAODkink::MomDaughterX() const {return fPx[0];}
inline Double_t AliAODkink::MomDaughterY() const {return fPy[0];}
inline Double_t AliAODkink::MomDaughterZ() const {return fPz[0];}

inline Double_t AliAODkink::RadiusKink() const {return fRadius;}  
inline Float_t  AliAODkink::KinkAngle() const {return fkinkAngle;}
inline Float_t  AliAODkink::KinkQt() const {return fQt;}

inline Double_t AliAODkink::PseudoRapMKink() const {
  return Eta();
}

//----------------------------------------------------------------------------

#endif
