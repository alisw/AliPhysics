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

#include <Riostream.h>

#include <TParticle.h>

#include "AliESDtrack.h"
#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

//--------------------------------------------------------------------------------------------------------
AliRsnDaughter::AliRsnDaughter()
{
// 
// Default constructor.
// Its unique argument defines how many PID weights are allowed. 
// Actually, it should be the same as AliESDtrack::kSPECIES (=5).
// 
	fSign = (Char_t)0;
	fPDG = (UShort_t)0;
	fIndex = (UShort_t)0;
	
	fP[0] = fP[1] = fP[2] = 0.0;
	fV[0] = fV[1] = fV[2] = 0.0;
	fMass = 0.0;
	
	Int_t i;
	for (i = 0; i < AliPID::kSPECIES; i++) {
		fPIDwgt[i] = 0.0;
	}
			
	fLabel = -1;
	fTruePDG = 0;
	fMother = -1;
	fMotherPDG = 0;
}
//--------------------------------------------------------------------------------------------------------
AliRsnDaughter::AliRsnDaughter(const AliRsnDaughter &copy) : TObject(copy)
{
//
// Copy constructor
//
	fSign = copy.fSign;
	fPDG = copy.fPDG;
	fIndex = copy.fIndex;
	
	Int_t i;
	for (i = 0; i < 3; i++) {
		fP[i] = copy.fP[i];
		fV[i] = copy.fV[i];
	}
	fMass = copy.fMass;
	
	for (i = 0; i < AliPID::kSPECIES; i++) {
		fPIDwgt[i] = copy.fPIDwgt[i];
	}
		
	fLabel = copy.fLabel;
	fTruePDG = copy.fTruePDG;
	fMother = copy.fMother;
	fMotherPDG = copy.fMotherPDG;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnDaughter::Adopt(const AliESDtrack* esdTrack, Bool_t checkITSRefit)
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
// Makes the following checks:
//
// - if argument 'checkITSRefit' is TRUE and the "ITS refit" flag is FALSE
//   in the track, the "fIsOK" flag of (this) is set to "false" (track should be rejected)
// - if the passed label is negative, track is considered as "fake", and 
//   this info is kept to allow fake tracks exclusion when doing analysis
//
	// check for refit in the ITS (if requested)
	if (checkITSRefit) {
		if ( !(esdTrack->GetStatus() & AliESDtrack::kITSrefit) ) {
			return kFALSE;
		}
	}
	
	// get sign and number of species allowed for PID
	fSign = (Char_t)esdTrack->GetSign();
	
	// get (and check) momentum
	esdTrack->GetPxPyPz(fP);
	if (fP[0] == 0.0 || fP[1] == 0.0 || fP[2] == 0.0) {
		return kFALSE;
	}
	
	// get (and check) vertex
	esdTrack->GetXYZ(fV);
	if (fV[0] == 0.0 || fV[1] == 0.0 || fV[2] == 0.0) {
		return kFALSE;
	}
	
	// get label 
	// (other kinematics informations are set to default and meaningless values)
	fLabel = esdTrack->GetLabel();
	
	// get PID weights
	esdTrack->GetESDpid(fPIDwgt);
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnDaughter::Adopt(TParticle* particle)
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
// When an AliRsnDaughter is copied from a TParticle, it is 
// considered always good for analysis and never fake.
//
	// get particle sign form the sign of PDG code
	Int_t pdg = particle->GetPdgCode();
	if (TMath::Abs(pdg) < 20) {
		if (pdg > 0) fSign = -1; else fSign = 1;
	}
	else if (TMath::Abs(pdg) < 3000) {
		if (pdg > 0) fSign = 1; else fSign = -1;
	}
	
	// get momentum
	fP[0] = particle->Px();
	fP[1] = particle->Py();
	fP[2] = particle->Pz();
	
	// get vertex
	fV[0] = particle->Vx();
	fV[1] = particle->Vy();
	fV[2] = particle->Vz();
	
	// set simulation data
	fPDG = fTruePDG = (Short_t)particle->GetPdgCode();
	fMother = particle->GetFirstMother();
	fMotherPDG = 0;
	
	return kTRUE;
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
	
	//TLorentzVector v1 = track1.Get4Momentum();
	//TLorentzVector v2 = track2.Get4Momentum();
	//TLorentzVector sm = v1 + v2;
	//Double_t pxTot = sum.X();
	//Double_t pyTot = sum.Y();
	//Double_t pzTot = sum.Z();
	//Double_t mass = sm.M();
	
	out.SetPxPyPz(pxTot, pyTot, pzTot);
	out.SetMass(mass);
	
	return out;
}
