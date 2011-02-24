#ifndef ALIAODRECODECAYHF2PRONG_H
#define ALIAODRECODECAYHF2PRONG_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliAODRecoDecayHF2Prong
// base class for AOD reconstructed 2-prong heavy-flavour decays
// (D0->Kpi, J/psi->ee, ...)
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//         G.E.Bruno, giuseppe.bruno@ba.infn.it 
//***********************************************************

#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"

class AliAODRecoDecayHF2Prong : public AliAODRecoDecayHF {

 public:

  AliAODRecoDecayHF2Prong();
  AliAODRecoDecayHF2Prong(AliAODVertex *vtx2,
			  Double_t *px,Double_t *py,Double_t *pz,
			  Double_t *d0,Double_t *d0err,Float_t dca);
  AliAODRecoDecayHF2Prong(AliAODVertex *vtx2,
			  Double_t *d0,Double_t *d0err,Float_t dca);    
  AliAODRecoDecayHF2Prong(const AliAODRecoDecayHF2Prong& source);
  AliAODRecoDecayHF2Prong& operator=(const AliAODRecoDecayHF2Prong& source); 

  virtual ~AliAODRecoDecayHF2Prong() {}  
 
  Double_t Prodd0d0() const {return AliAODRecoDecay::Prodd0d0(0,1);} 

  // D0->Kpi
  Double_t ED0() const {return E(421);} 
  Double_t YD0() const {return Y(421);} 

  Double_t CtD0() const {return Ct(421);} 
  Double_t CtD0(Double_t point[3]) const {return AliAODRecoDecay::Ct(421,point);}
  Double_t CtD0(AliAODVertex *vtx1) const {return AliAODRecoDecay::Ct(421,vtx1);}

  Double_t CosThetaStarD0() const {return CosThetaStar(1,421,211,321);} // angle of K
  Double_t CosThetaStarD0bar() const {return CosThetaStar(0,421,321,211);} // angle of K
  void CosThetaStarD0(Double_t &ctsD0,Double_t &ctsD0bar) const 
    {ctsD0=CosThetaStarD0();ctsD0bar=CosThetaStarD0bar();return;}

  Double_t InvMassD0() const {UInt_t pdg[2]={211,321};return InvMass(2,pdg);}
  Double_t InvMassD0bar() const {UInt_t pdg[2]={321,211};return InvMass(2,pdg);}
  void InvMassD0(Double_t &mD0,Double_t &mD0bar) const 
    {mD0=InvMassD0();mD0bar=InvMassD0bar();return;}

  Bool_t   SelectD0(const Double_t* cuts,Int_t &okD0,Int_t &okD0bar) const;

  // Jpsi (from B) -> ee
  Double_t EJPSI() const {return E(443);}
  Double_t YJPSI() const {return Y(443);}

  Double_t CtJPSI() const {return Ct(443);}
  Double_t CtJPSI(Double_t point[3]) const {return AliAODRecoDecay::Ct(443,point);}
  Double_t CtJPSI(AliAODVertex *vtx1) const {return AliAODRecoDecay::Ct(443,vtx1);}

  Double_t CosThetaStarJPSI() const {return CosThetaStar(1,443,11,11);} // angle of e-

  Double_t InvMassJPSIee() const {UInt_t pdg[2]={11,11};return InvMass(2,pdg);}

  Bool_t   SelectBtoJPSI(const Double_t* cuts,Int_t &okB) const;
  

 private:

  ClassDef(AliAODRecoDecayHF2Prong,1)  // base class for AOD reconstructed 
                                       // heavy-flavour 2-prong decays
};

#endif
