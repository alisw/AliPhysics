/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
//-------------------------------------------------------------------------
//                      Class AliRsnDaughter
//                     ----------------------
//           A simple object which describes a reconstructed track
//           with some references to its related generated particle
//           and some facilities which could help in composing its
//           4-momentum, for resonance study.
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#include <TParticle.h>
#include <TString.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

//--------------------------------------------------------------------------------------------------------
AliRsnDaughter::AliRsnDaughter() :
  TObject(),
  fSign((Char_t)0),
  fPDG((UShort_t)0),
  fIndex((UShort_t)0),
  fMass(0.0),
  fLabel(-1),
  fTruePDG((Short_t)0),
  fMother(-1),
  fMotherPDG((Short_t)0)
{
// 
// Default constructor.
// Its unique argument defines how many PID weights are allowed. 
// Actually, it should be the same as AliESDtrack::kSPECIES (=5).
// 
	
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		if (i < 3) {
			fP[i] = 0.0;
			fV[i] = 0.0;
		}
		fPIDwgt[i] = 0.0;
	}
}
//--------------------------------------------------------------------------------------------------------
AliRsnDaughter::AliRsnDaughter(Int_t label, UShort_t index, Double_t *p, Double_t *v, Char_t sign) :
  TObject(),
  fSign(sign),
  fPDG((UShort_t)0),
  fIndex(index),
  fMass(0.0),
  fLabel(label),
  fTruePDG((Short_t)0),
  fMother(-1),
  fMotherPDG((Short_t)0)
{
//
// [PRIVATE]
// Constructor with arguments.
// To create an object with all its standard data members filled (the ones coming from ESD)
//
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		if (i < 3) {
			fP[i] = p[i];
			fV[i] = v[i];
		}
		fPIDwgt[i] = 0.0;
	}
}
//--------------------------------------------------------------------------------------------------------
AliRsnDaughter * AliRsnDaughter::Adopt(AliESDtrack* esdTrack, Int_t index)
{
//
// Copies reconstructed data from an AliESDtrack:
//
// - charge sign
// - momentum
// - point of closest approach to primary vertex
// - ESD pid weights
// - track label (AliESDtrack::GetLabel())
//
// Given that it happened that tracks have exactly zero momentum and DCA point,
// these quantities are checked and skipped when caught.
//
	Double_t p[3], v[3], pid[AliPID::kSPECIES];
	
	esdTrack->GetPxPyPz(p);
	esdTrack->GetXYZ(v);
	esdTrack->GetESDpid(pid);
	
	if (p[0] == 0.0 || p[1] == 0.0 || p[2] == 0.0) return 0x0;
	if (v[0] == 0.0 || v[1] == 0.0 || v[2] == 0.0) return 0x0;
		
	AliRsnDaughter *out = new AliRsnDaughter(esdTrack->GetLabel(), (UShort_t)index, p, v, (Char_t)esdTrack->GetSign());
	
	// store PID weights
	out->SetPIDweights(pid);
	
	return out;
}
//--------------------------------------------------------------------------------------------------------
AliRsnDaughter * AliRsnDaughter::Adopt(TParticle* particle, Int_t label)
{
//
// Copies data from a generated particle:
// 
// - PDG code
// - charge sign
// - momentum
// - production vertex
// - GEANT label of mother track
//
	Int_t    pdg;
	Char_t   sign;
	Double_t p[3], v[3];
	
	// get particle sign form the sign of PDG code
	pdg = particle->GetPdgCode();
	if (TMath::Abs(pdg) < 20) {
		if (pdg > 0) sign = -1; else sign = 1;
	}
	else if (TMath::Abs(pdg) < 3000) {
		if (pdg > 0) sign = 1; else sign = -1;
	}
	else {
		return 0x0;
	}
	
	p[0] = particle->Px();
	p[1] = particle->Py();
	p[2] = particle->Pz();
	
	v[0] = particle->Vx();
	v[1] = particle->Vy();
	v[2] = particle->Vz();
	
	AliRsnDaughter *out = new AliRsnDaughter(label, (UShort_t)TMath::Abs(label), p, v, sign);
	
	// set simulation data
	out->SetPDG(TMath::Abs(pdg));
	out->SetTruePDG(TMath::Abs(pdg));
	out->SetMother(particle->GetFirstMother());
	
	return out;
}
//--------------------------------------------------------------------------------------------------------
void AliRsnDaughter::Print(Option_t *option) const
{
// 
// Print informations about track.
// All options are used to add some specifical information as required:
// "S" --> charge sign
// "L" --> track label
// "I" --> assigned PID (PDG code)
// "T" --> true PDG code
// "V" --> coordinates of track vertex
// "P" --> coordinates of track momentum
// "M" --> mother label
// "N" --> mother PDG code
// "W" --> ESD PID weights
//

	TString output("Track info: ");
	TString opt(option);
	opt.ToUpper();
	
	if (opt.Contains("S")) {
		output.Append(Form("sign = %d -- ", (Int_t)fSign));
	}
	if (opt.Contains("L")) {
		output.Append(Form("label = %d -- ", fLabel));
	}
	if (opt.Contains("I")) {
		output.Append(Form("PDG = %d -- ", fPDG));
	}
	if (opt.Contains("T")) {
		output.Append(Form("true PDG = %d -- ", fTruePDG));
	}
	if (opt.Contains("V")) {
		output.Append(Form("v = %f, %f, %f -- ", fV[0], fV[1], fV[2]));
	}
	if (opt.Contains("P")) {
		output.Append(Form("p = %f, %f, %f -- ", fP[0], fP[1], fP[2]));
	}
	if (opt.Contains("M")) {
		output.Append(Form("mum = %d -- ", fMother));
	}
	if (opt.Contains("N")) {
		output.Append(Form("mum PDG = %d -- ", fMotherPDG));
	}
	if (opt.Contains("W")) {
		output.Append(Form("PID wgts (e, mu, pi, K, p) = %f, %f, %f, %f, %f", fPIDwgt[0], fPIDwgt[1], fPIDwgt[2], fPIDwgt[3], fPIDwgt[4]));
	}
	
	AliInfo(output.Data());
}
//--------------------------------------------------------------------------------------------------------
AliRsnDaughter AliRsnDaughter::Sum(AliRsnDaughter t1, AliRsnDaughter t2)
{
//
// Builds a new AliRsnDaughter object with the sum of momenta of two particles.
//
	// create new AliRsnDaughter with default useless values for datamembers
	AliRsnDaughter out;
	
	// if the summed particles are daughters of the same resonance
	// their common mother label becomes the label of the sum
	Int_t mum1 = t1.GetMother();
	Int_t mum2 = t2.GetMother();
	if (mum1 == mum2) {
		out.SetLabel(mum1);
		out.SetMotherPDG(t1.GetMotherPDG());
	}
	else {
		out.SetLabel(-1);
		out.SetMotherPDG(0);
	}

	// compute total 4-momentum
	Double_t etot  = t1.GetEnergy() + t2.GetEnergy();
	Double_t pxTot = t1.GetPx() + t2.GetPx();
	Double_t pyTot = t1.GetPy() + t2.GetPy();
	Double_t pzTot = t1.GetPz() + t2.GetPz();
	Double_t mass  = TMath::Sqrt(etot*etot - pxTot*pxTot - pyTot*pyTot - pzTot*pzTot);
	
	out.SetPxPyPz(pxTot, pyTot, pzTot);
	out.SetMass(mass);
	
	return out;
}
