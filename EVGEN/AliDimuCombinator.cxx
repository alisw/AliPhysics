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

/*
$Log$
Revision 1.9  2001/03/08 13:30:43  morsch
Make it work with particle stack of V3.05.

Revision 1.8  2000/12/21 16:24:06  morsch
Coding convention clean-up

Revision 1.7  2000/10/02 15:16:37  morsch
Correct coding rule violation for member data names of type fi -> fI.

Revision 1.6  2000/06/14 15:19:47  morsch
Include clean-up (IH)

Revision 1.5  2000/06/09 20:35:32  morsch
All coding rule violations except RS3 corrected

Revision 1.4  2000/03/20 18:03:24  morsch
Change muon particle code to PDG code.

Revision 1.3  1999/09/29 09:24:08  fca
Introduction of the Copyright and cvs Log

*/

/*  
 Class for dimuon analysis and fast dimuon simulation.
 It provides single and dimuon iterators, cuts, weighting, kinematic
 It uses the AliRun particle tree.
 Comments and suggestions to 
 andreas.morsch@cern.ch
*/


#include "AliDimuCombinator.h" 
#include "AliPDG.h" 
#include "AliRun.h" 
#include <TRandom.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TTree.h>

//
ClassImp(AliDimuCombinator)
    AliDimuCombinator::AliDimuCombinator() 
{
// Constructor
    fNParticle = (Int_t) (gAlice->TreeK())->GetEntries();
    fImuon1 = 0;
    fImuon2 = 0;
    fMuon1  = 0;
    fMuon2  = 0;
    fImin1  = 0;
    fImin2  = 0;
    fImax1  = fNParticle;
    fImax2  = fNParticle;
    fPtMin  = 0;
    fEtaMin = -10;
    fEtaMax = -10;
    fRate1  = 1.;
    fRate2  = 1.;
}

AliDimuCombinator::AliDimuCombinator(const AliDimuCombinator & combinator)
{
// Dummy copy constructor
}


//
//                       Iterators
// 
TParticle* AliDimuCombinator::Particle(Int_t i)
{
    return gAlice->Particle(i);
}

TParticle* AliDimuCombinator::FirstMuon()
{
// Single muon iterator: initialisation
    fImuon1 = fImin1;
    fMuon1  = Particle(fImuon1);
    while(Type(fMuon1) != kMuonPlus && Type(fMuon1) != kMuonMinus) {
	fImuon1++;
	if (fImuon1 >= fImax1) {fMuon1 = 0; break;}
	fMuon1 = Particle(fImuon1);
    }
    return fMuon1;
}

TParticle* AliDimuCombinator::FirstMuonSelected()
{
// Single selected muon iterator: initialisation
    TParticle* muon = FirstMuon();
    while(muon != 0 && !Selected(muon)) {muon = NextMuon();}
    return muon;
}


TParticle* AliDimuCombinator::NextMuon()
{
// Single muon iterator: increment
    fImuon1++;
    if (fImuon1 >= fNParticle) {fMuon1 = 0; return fMuon1;}
    
    fMuon1 = Particle(fImuon1);
    while(Type(fMuon1) != kMuonPlus && Type(fMuon1) != kMuonMinus) {
	fImuon1++;
	if (fImuon1 >= fImax1) {fMuon1 = 0; break;}
	fMuon1 = Particle(fImuon1);
    }
    return fMuon1;
}

TParticle* AliDimuCombinator::NextMuonSelected()
{
// Single selected muon iterator: increment
    TParticle * muon = NextMuon();
    while(muon !=0 && !Selected(muon)) {muon = NextMuon();}
    return muon;
}


void AliDimuCombinator::FirstPartner()
{
// Helper for  dimuon iterator: initialisation
    if (fImin1 == fImin2) {
	fImuon2 = fImuon1+1;
    } else {
	fImuon2 = fImin2;
    }
    if (fImuon2 >= fImax2) {fMuon2 = 0; return;}
    fMuon2 = Particle(fImuon2);
    while(Type(fMuon2) != kMuonPlus && Type(fMuon2) != kMuonMinus) {
	fImuon2++;
	if (fImuon2 >= fImax2) {fMuon2 = 0; break;}
	fMuon2 = Particle(fImuon2);
    }
}

void AliDimuCombinator::FirstPartnerSelected()
{
// Helper for selected dimuon iterator: initialisation
    FirstPartner();
    while(fMuon2 !=0 && !Selected(fMuon2)) {NextPartner();}
}


void AliDimuCombinator::NextPartner()
{
// Helper for dimuon iterator: increment    
    fImuon2++;
    if (fImuon2 >= fImax2) {fMuon2 = 0; return;}
    
    
    fMuon2 = Particle(fImuon2);
    
    while(Type(fMuon2) != kMuonPlus && Type(fMuon2) != kMuonMinus) {
	fImuon2++;
	if (fImuon2 >= fImax2) {fMuon2 = 0; break;}
	fMuon2 = Particle(fImuon2);
    }
}

void AliDimuCombinator::NextPartnerSelected()
{
// Helper for selected dimuon iterator: increment    
    NextPartner();
    while(fMuon2 !=0 && !Selected(fMuon2)) {NextPartner();}
}


TParticle*  AliDimuCombinator::Partner()
{
// Returns current partner for muon to form a dimuon
    return fMuon2;
}

void AliDimuCombinator::FirstMuonPair(TParticle* & muon1, TParticle* & muon2)
{
// Dimuon iterator: initialisation
    FirstMuon();
    FirstPartner();
    muon1 = fMuon1;
    muon2 = fMuon2;	 
}

void AliDimuCombinator::NextMuonPair(TParticle* & muon1, TParticle* & muon2)
{
// Dimuon iterator: increment    
    NextPartner();
    if (!Partner()) {
	NextMuon();
	FirstPartner();
    }
    muon1 = fMuon1;
    muon2 = fMuon2;	 
}
void AliDimuCombinator::FirstMuonPairSelected(TParticle* & muon1, 
					      TParticle* & muon2)
{
// Selected dimuon iterator: initialisation    
    FirstMuonSelected();
    FirstPartnerSelected();
    muon1 = fMuon1;
    muon2 = fMuon2;	 
}

void AliDimuCombinator::NextMuonPairSelected(TParticle* & muon1, 
					     TParticle* & muon2)
{
// Selected dimuon iterator: increment    
    NextPartnerSelected();
    if (!Partner()) {
	NextMuonSelected();
	FirstPartnerSelected();
    }
    muon1 = fMuon1;
    muon2 = fMuon2;	 
}

void AliDimuCombinator::ResetRange()
{
// Reset index ranges for single muons
    fImin1 = fImin2 = 0;
    fImax1 = fImax2 = fNParticle;
}

void AliDimuCombinator::SetFirstRange(Int_t from, Int_t to)
{
// Reset index range for first muon
    fImin1 = from;
    fImax1 = to;
    if (fImax1 > fNParticle) fImax1 = fNParticle;
}

void AliDimuCombinator::SetSecondRange(Int_t from, Int_t to)
{
// Reset index range for second muon
    fImin2 = from;
    fImax2 = to;
    if (fImax2 > fNParticle) fImax2 = fNParticle;
}
//
//                       Selection
//

Bool_t AliDimuCombinator::Selected(TParticle* part)
{
// Selection cut for single muon 
//
    if (part == 0) {return 0;}
    
    if (part->Pt() > fPtMin && part->Eta() > fEtaMin && part->Eta() < fEtaMax) {
	return 1;
    } else {
	return 0;
    }
}

Bool_t AliDimuCombinator::Selected(TParticle* part1, TParticle* part2)
{
// Selection cut for dimuons
//
     return Selected(part1)*Selected(part2);
}
//
//                       Kinematics
//
Float_t AliDimuCombinator::Mass(TParticle* part1, TParticle* part2)
{
// Invariant mass
//
    Float_t px,py,pz,e;
    px = part1->Px()+part2->Px();
    py = part1->Py()+part2->Py();
    pz = part1->Pz()+part2->Pz();    
    e  = part1->Energy()+part2->Energy();
    Float_t p = px*px+py*py+pz*pz;
    if (e*e < p) {
	return -1; 
    } else {
	return TMath::Sqrt(e*e-p);
    }
}

Float_t AliDimuCombinator::PT(TParticle* part1, TParticle* part2)
{
// Transverse momentum of dimuons
//
    Float_t px,py;
    px = part1->Px()+part2->Px();
    py = part1->Py()+part2->Py();
    return TMath::Sqrt(px*px+py*py);
}

Float_t AliDimuCombinator::Pz(TParticle* part1, TParticle* part2)
{
// Pz of dimuon system
//
    return part1->Pz()+part2->Pz();
}

Float_t AliDimuCombinator::Y(TParticle* part1, TParticle* part2)
{
// Rapidity of dimuon system
//
    Float_t pz,e;
    pz = part1->Pz()+part2->Pz();
    e  = part1->Energy()+part2->Energy();
    return 0.5*TMath::Log((e+pz)/(e-pz));
}
//                  Response
//
void AliDimuCombinator::SmearGauss(Float_t width, Float_t & value)
{
// Apply gaussian smearing
//
    value+=gRandom->Gaus(0, width);
}
//              Weighting
// 

Float_t AliDimuCombinator::DecayProbability(TParticle* part)
{
// Calculate decay probability for muons from pion and kaon decays
// 
    Float_t d, h, theta, cTau;
    TParticle* parent = Parent(part);
    Int_t ipar = Type(parent);
    if (ipar == kPiPlus || ipar == kPiMinus) {
	cTau=780.4;
    } else if (ipar == kKPlus || ipar == kKMinus) {
	cTau = 370.9;
    } else {
	cTau = 0;
    }
    
    
    Float_t gammaBeta=(parent->P())/(parent->GetMass());
//
// this part is still very ALICE muon-arm specific
//
    theta = parent->Theta();
    h = 90*TMath::Tan(theta);
    
    if (h<4) {
	d=4/TMath::Sin(theta);
    } else {
	d=90/TMath::Cos(theta);
    }
    
    if (cTau > 0) {
	return 1-TMath::Exp(-d/cTau/gammaBeta);
    } else {
	return 1;
    }
}

Float_t AliDimuCombinator::Weight(TParticle* part1, TParticle* part2)
{
// Dimuon weight

    Float_t wgt = (part1->GetWeight())*(part2->GetWeight());
    
    if (Correlated(part1, part2)) {
	if ( part1->GetFirstMother() == part2->GetFirstMother()) {
	    return part1->GetWeight()*fRate1;
	} else {
	    return wgt/(Parent(part1)->GetWeight())*fRate1;
	}
    } else {
	return wgt*fRate1*fRate2;
    }
} 


Float_t AliDimuCombinator::Weight(TParticle* part)
{
// Single muon weight
    return (part->GetWeight())*(Parent(part)->GetWeight())*fRate1;
}

Bool_t  AliDimuCombinator::Correlated(TParticle* part1, TParticle* part2)
{
// Check if muons are correlated
//
    if ((Origin(part1) >= 0) && Origin(part1) == Origin(part2)) {
/*
	printf("\n origin %d %d ", 
	       Type(Particle(Origin(part1))),
	       Type(Particle(Origin(part2))));
	printf("\n parent %d %d \n \n ", 
	       Type(Parent(part1)),
	       Type(Parent(part2)));
*/	
	return kTRUE;
    } else {
	return kFALSE;
    }
}

TParticle* AliDimuCombinator::Parent(TParticle* part)
{
// Return pointer to parent
//
    return Particle(part->GetFirstMother());
}

Int_t AliDimuCombinator::Origin(TParticle* part)
{
// Return pointer to primary particle
//
    Int_t iparent= part->GetFirstMother();
    if (iparent < 0) return iparent;
    Int_t ip;
    while(1) {
	ip = (Particle(iparent))->GetFirstMother();
	if (ip < 0) {
	    break;
	} else {
	    iparent = ip;
	}
    }
    return iparent;
}

Int_t AliDimuCombinator::Type(TParticle *part) 
{
// Return particle type for 
return part->GetPdgCode();
}

AliDimuCombinator& AliDimuCombinator::operator=(const  AliDimuCombinator& rhs)
{
// Assignment operator
    return *this;
}


void AliDimuCombinator::Copy(AliDimuCombinator &combi) const
{
  //
  // Copy *this onto lego -- not implemented
  //
  Fatal("Copy","Not implemented!\n");
}





