/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnDaughter
//
//           A simple object which describes a reconstructed track
//           with some references to its related generated particle
//           and some facilities which could help in composing its
//           4-momentum, for resonance study.
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNDAUGHTER_H
#define ALIRSNDAUGHTER_H

#include "AliPID.h"
#include <TVector3.h>
#include <TLorentzVector.h>

class TParticle;
class AliESDtrack;

class AliRsnDaughter : public TObject
{
public:
	
	AliRsnDaughter();
	
	virtual       ~AliRsnDaughter() { }
	
	// 4-momentum
	Double_t       GetEnergy() const {return TMath::Sqrt(fMass*fMass + GetP2());}
	Double_t       GetMass() const {return fMass;}
	Double_t       GetP2() const {return fP[0]*fP[0] + fP[1]*fP[1] + fP[2]*fP[2];}
	Double_t       GetP() const  {return TMath::Sqrt(GetP2());}
	Double_t       GetPx() const {return fP[0];}
	Double_t       GetPy() const {return fP[1];}
	Double_t       GetPz() const {return fP[2];}
	Double_t       GetPt() const {return TMath::Sqrt(fP[0]*fP[0] + fP[1]*fP[1]);}
	void           SetPxPyPz(Double_t px, Double_t py, Double_t pz) {fP[0] = px; fP[1] = py; fP[2] = pz;}
	void           SetMass(Double_t m) {fMass = m;}
	
	// DCA vertex
	Double_t       GetVx() const {return fV[0];}
	Double_t       GetVy() const {return fV[1];}
	Double_t       GetVz() const {return fV[2];}
	Double_t       GetVt() const {return TMath::Sqrt(fV[0]*fV[0] + fV[1]*fV[1]);}
	void           SetVxVyVz(Double_t vx,Double_t vy,Double_t vz) {fV[0]=vx;fV[1]=vy;fV[2]=vz;}
	
	// charge sign
	Char_t         GetSign() const {return fSign;}
	void           SetSign(Char_t value) {fSign = value;}
	void           SetSign(Int_t value) {fSign = (Char_t)value;}
	
	// PID
	UShort_t       GetPDG() const {return fPDG;}
	Double_t       GetPIDweight(Int_t i) const {return ( (i>=0&&i<AliPID::kSPECIES)?fPIDwgt[i]:-1.0 );}
	void           SetPDG(Int_t pdg) {fPDG = (UShort_t)TMath::Abs(pdg);}
	void           SetPDG(UShort_t pdg) {fPDG = TMath::Abs(pdg);}
	void           SetPIDweights(const Double_t *pid) {Int_t i;for(i=0;i<AliPID::kSPECIES;i++)fPIDwgt[i]=pid[i];}
	
	// position in stack/ESD array
	UShort_t       GetIndex() const {return fIndex;}
	Int_t          GetLabel() const {return fLabel;}
	void           SetIndex(UShort_t value) {fIndex = value;}
	void           SetIndex(Int_t value) {fIndex = (UShort_t)value;}
	void           SetLabel(Int_t l) {fLabel = l;}
	
	// MC info (when available)
	Int_t          GetMother() const {return fMother;}
	Short_t        GetMotherPDG() const {return fMotherPDG;}
	Short_t        GetTruePDG() const {return fTruePDG;}
	void           SetMother(Int_t l) {fMother = l;}
	void           SetMotherPDG(Short_t pdg) {fMotherPDG = pdg;}
	void           SetTruePDG(Short_t pdg) {fTruePDG = pdg;}
	
	// print info
	void           Print(Option_t *option = "") const;
	
	// static methods
	static AliRsnDaughter * Adopt(AliESDtrack* track, Int_t index);
	static AliRsnDaughter * Adopt(TParticle* particle, Int_t label);
	static AliRsnDaughter   Sum(AliRsnDaughter t1, AliRsnDaughter t2);
	
private:
	
	AliRsnDaughter(Int_t label, UShort_t index, Double_t *p, Double_t *v, Char_t sign);
	
	Char_t     fSign;                     // charge sign
	UShort_t   fPDG;                      // assigned PDG code from PID (0=undefined)
	UShort_t   fIndex;                    // reference index in AliESD container

	Double_t   fP[3];                     // vector momentum
	Double_t   fV[3];                     // production vertex
	Double_t   fMass;                     // mass
	
	Double_t   fPIDwgt[AliPID::kSPECIES]; // particle PID weights
	
	// The following data are useful for simulated events only
	// and can be left blank when using real or 'realistic' data
	
	Int_t      fLabel;        // GEANT label of corresponding particle
	Short_t    fTruePDG;      // PDG code of corresponding particle
	Int_t      fMother;       // GEANT label of mother particle
	Short_t    fMotherPDG;    // PDG code of mother particle
	
	ClassDef(AliRsnDaughter,1)
};

#endif
