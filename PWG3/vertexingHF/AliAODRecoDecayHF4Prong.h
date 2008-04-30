#ifndef ALIAODRECODECAYHF4PRONG_H
#define ALIAODRECODECAYHF4PRONG_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliAODRecoDecayHF4Prong
// base class for AOD reconstructed 4-prong heavy-flavour decays
// (D0->Kpipipi, etc...)
// Authors: G.E.Bruno Giuseppe.Bruno@ba.infn.it, R.Romita Rossella.Romita@ba.infn.it
//***********************************************************

#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"

class AliAODRecoDecayHF4Prong : public AliAODRecoDecayHF {

 public:

  AliAODRecoDecayHF4Prong();
  AliAODRecoDecayHF4Prong(AliAODVertex *vtx2,
			   Double_t *px,Double_t *py,Double_t *pz,
			   Double_t *d0,Double_t *d0err,
			   Double_t *dca, //Double_t sigvert,
			   Double_t dist12,Double_t dist23, 
			   Double_t dist14,Double_t dist34, 
                           Short_t charge);
   AliAODRecoDecayHF4Prong(AliAODVertex *vtx2,
			   Double_t *d0,Double_t *d0err,
			   Double_t *dca, //Double_t sigvert,
			   Double_t dist12,Double_t dist23, 
			   Double_t dist14,Double_t dist34, 
			   Short_t charge);

  AliAODRecoDecayHF4Prong(const AliAODRecoDecayHF4Prong& source);
  AliAODRecoDecayHF4Prong& operator=(const AliAODRecoDecayHF4Prong& source); 

  virtual ~AliAODRecoDecayHF4Prong() {}  
 
  void GetDCAs(Double_t dca[6]) const 
    {for(Int_t i=0;i<6;i++) dca[i]=GetDCA(i);} // convention:fDCA[0]=p0p1,fDCA[1]=p0p2,fDCA[2]=p0p3,fDCA[3]=p1p2...
  Double_t GetDist12toPrim() const {return fDist12toPrim;}
  Double_t GetDist23toPrim() const {return fDist12toPrim;}
  Double_t GetDist14toPrim() const {return fDist12toPrim;}
  Double_t GetDist34toPrim() const {return fDist12toPrim;}

  // D0->pi+K- pipi and D0bar->K+pi- pipi (in the order given) 
  Double_t ED0() const {return E(421);}
  Double_t YD0() const {return Y(421);} 
  Double_t CtD0() const {return Ct(421);} 
  Double_t CtD0(Double_t point[3]) const {return AliAODRecoDecay::Ct(421,point);}
  Double_t CtD0(AliAODVertex *vtx1) const {return AliAODRecoDecay::Ct(421,vtx1);}
  Double_t InvMassRho() const {return InvMass2Prongs(2,3,211,211);} 
  Double_t InvMassD0() const {UInt_t pdg[4]={211,321,211,211};return InvMass(4,pdg);}
  Double_t InvMassD0bar() const {UInt_t pdg[4]={321,211,211,211};return InvMass(4,pdg);}
  void InvMassD0(Double_t &mD0,Double_t &mD0bar) const
    {mD0=InvMassD0();mD0bar=InvMassD0bar();return;}


  Bool_t   SelectD0(const Double_t* cuts,Int_t &okD0,Int_t &okD0bar) const;


 private:

  //Double_t fSigmaVert; // track dispersion around the secondary vertex
  Double_t fDist12toPrim; //distance prim vert - 2 opposite sign track vertex 
  Double_t fDist23toPrim; //distance prim vert - 2 opposite sign track vertex 
  Double_t fDist14toPrim; //distance prim vert - 2 opposite sign track vertex 
  Double_t fDist34toPrim; //distance prim vert - 2 opposite sign track vertex 
  //Double_t fDist123toPrim;  //distance 

  ClassDef(AliAODRecoDecayHF4Prong,1)  // base class for AOD reconstructed 
                                       // heavy-flavour 3-prong decays
};

#endif
