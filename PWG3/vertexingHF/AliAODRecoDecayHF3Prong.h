#ifndef ALIAODRECODECAYHF3PRONG_H
#define ALIAODRECODECAYHF3PRONG_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliAODRecoDecayHF3Prong
// base class for AOD reconstructed 3-prong heavy-flavour decays
// (D+->Kpipi, Ds->KKpi ...)
// Author: E.Bruna bruna@to.infn.it, F.Prino prino@to.infn.it
//***********************************************************

#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODEvent.h"

class AliAODRecoDecayHF3Prong : public AliAODRecoDecayHF {

 public:
  
  AliAODRecoDecayHF3Prong();
   AliAODRecoDecayHF3Prong(AliAODVertex *vtx2,
			   Double_t *px,Double_t *py,Double_t *pz,
			   Double_t *d0,Double_t *d0err,
			   Double_t *dca, Double_t sigvert,
			   Double_t dist12,Double_t dist23,Short_t charge);
   AliAODRecoDecayHF3Prong(AliAODVertex *vtx2,
			   Double_t *d0,Double_t *d0err,
			   Double_t *dca, Double_t sigvert,
			   Double_t dist12,Double_t dist23, Short_t charge);

  AliAODRecoDecayHF3Prong(const AliAODRecoDecayHF3Prong& source);
  AliAODRecoDecayHF3Prong& operator=(const AliAODRecoDecayHF3Prong& source); 
  void GetDCAs(Double_t dca[3]) const 
    {for(Int_t i=0;i<3;i++) dca[i]=GetDCA(i);}
  Double_t GetSigmaVert(AliAODEvent* aod=0x0) { 
    if(fSigmaVert>0.00001) return fSigmaVert; 
    if(aod) fSigmaVert=ComputeSigmaVert(aod);
    return fSigmaVert;
  }
  Double_t ComputeSigmaVert(AliAODEvent* aod) const;
  Double_t GetDist12toPrim() const { return fDist12toPrim; }
  Double_t GetDist23toPrim() const { return fDist23toPrim; }
  void SetDist12toPrim(Double_t d) { fDist12toPrim=d; }
  void SetDist23toPrim(Double_t d) { fDist23toPrim=d; }


  // D+->Kpipi
  Double_t EDplus() const {return E(411);} 
  Double_t YDplus() const {return Y(411);} 
  Double_t CtDplus() const {return Ct(411);} 
  Double_t CtDplus(Double_t point[3]) const {return AliAODRecoDecay::Ct(411,point);}
  Double_t CtDplus(AliAODVertex *vtx1) const {return AliAODRecoDecay::Ct(411,vtx1);}
  Double_t InvMassDplus() const {UInt_t pdg[3]={211,321,211};return InvMass(3,pdg);}
  Bool_t   SelectDplus(const Double_t* cuts) const;

  // Ds+->KKpi
  Double_t EDs() const {return E(431);} 
  Double_t YDs() const {return Y(431);} 
  Double_t CtDs() const {return Ct(431);} 
  Double_t CtDs(Double_t point[3]) const {return AliAODRecoDecay::Ct(431,point);}
  Double_t CtDs(AliAODVertex *vtx1) const {return AliAODRecoDecay::Ct(431,vtx1);}
  Double_t InvMassDsKKpi() const {UInt_t pdg[3]={321,321,211};return InvMass(3,pdg);}
  Double_t InvMassDspiKK() const {UInt_t pdg[3]={211,321,321};return InvMass(3,pdg);}
  
  Double_t CosPiKPhiRFrameKKpi() const {return CosPiKPhiRFrame(0);}
  Double_t CosPiKPhiRFramepiKK() const {return CosPiKPhiRFrame(1);}
  Double_t CosPiDsLabFrameKKpi() const {return CosPiDsLabFrame(0);}
  Double_t CosPiDsLabFramepiKK() const {return CosPiDsLabFrame(1);}
  
  Double_t CosPiKPhiRFrame(Int_t option) const;
  Double_t CosPiDsLabFrame(Int_t option) const;
  Bool_t   SelectDs(const Double_t* cuts,Int_t &okDsKKpi,Int_t &okDspiKK, Int_t &okMassPhi, Int_t &okMassK0star) 
    const; // same variables as D+, for now

  // Lambdac+->pKpi
  Double_t ELc() const {return E(4122);} 
  Double_t YLc() const {return Y(4122);} 
  Double_t CtLc() const {return Ct(4122);} 
  Double_t CtLc(Double_t point[3]) const {return AliAODRecoDecay::Ct(4122,point);}
  Double_t CtLc(AliAODVertex *vtx1) const {return AliAODRecoDecay::Ct(4122,vtx1);}
  Double_t InvMassLcpKpi() const {UInt_t pdg[3]={2212,321,211};return InvMass(3,pdg);}
  Double_t InvMassLcpiKp() const {UInt_t pdg[3]={211,321,2212};return InvMass(3,pdg);}
  Bool_t   SelectLc(const Double_t* cuts,Int_t &okLcpKpi,Int_t &okLcpiKp) 
    const; // same variables as D+, for now

 private:

  Double_t fSigmaVert; // track dispersion around the secondary vertex
  Double_t fDist12toPrim; //distance prim vert - 2 opposite sign track vertex 
  Double_t fDist23toPrim; //distance prim vert - 2 opposite sign track vertex



  ClassDef(AliAODRecoDecayHF3Prong,1)  // base class for AOD reconstructed 
                                       // heavy-flavour 3-prong decays
};

#endif
