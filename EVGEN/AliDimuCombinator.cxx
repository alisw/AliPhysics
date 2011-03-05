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

/* $Id$ */

// 
// Class for dimuon analysis and fast dimuon simulation.
// It provides single and dimuon iterators, cuts, weighting, kinematic
// It uses the AliRun particle tree.
// Comments and suggestions to 
// andreas.morsch@cern.ch


#include <TClonesArray.h>
#include <TParticle.h>
#include <TPDGCode.h> 
#include <TRandom.h>
#include <TTree.h>

#include "AliDimuCombinator.h" 
#include "AliRun.h" 
#include "AliMC.h"

//
ClassImp(AliDimuCombinator)

AliDimuCombinator::AliDimuCombinator():
    fNParticle((Int_t) (AliRunLoader::Instance()->TreeK())->GetEntries()),
    fImuon1(0),
    fImuon2(0),
    fImin1(0),
    fImin2(0),
    fImax1(fNParticle),
    fImax2(fNParticle),
    fRate1(1.),
    fRate2(1.),
    fMuon1(0),
    fMuon2(0),
    fPtMin(0.),
    fEtaMin(-10.),
    fEtaMax(10.)
{
// Constructor
    fNParticle = (Int_t) (AliRunLoader::Instance()->TreeK())->GetEntries();
}

AliDimuCombinator::AliDimuCombinator(const AliDimuCombinator & combinator)
    :TObject(combinator),
     fNParticle(0),
     fImuon1(0),
     fImuon2(0),
     fImin1(0),
     fImin2(0),
     fImax1(0),
     fImax2(0),
     fRate1(0),
     fRate2(0),
     fMuon1(0),
     fMuon2(0),
     fPtMin(0.),
     fEtaMin(0.),
     fEtaMax(0.)
{
// Dummy copy constructor
    combinator.Copy(*this);
}


//
//                       Iterators
// 
TParticle* AliDimuCombinator::Particle(Int_t i) const
{
//  Return next particle
//
    return gAlice->GetMCApp()->Particle(i);
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


TParticle*  AliDimuCombinator::Partner() const
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

Bool_t AliDimuCombinator::Selected(const TParticle* part) const
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

Bool_t AliDimuCombinator::Selected(const TParticle* part1, const TParticle* part2) const
{
// Selection cut for dimuons
//
     return Selected(part1)*Selected(part2);
}
//
//                       Kinematics
//
Float_t AliDimuCombinator::Mass(const TParticle* part1, const TParticle* part2) const
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

Float_t AliDimuCombinator::PT(const TParticle* part1, const TParticle* part2) const
{
// Transverse momentum of dimuons
//
    Float_t px,py;
    px = part1->Px()+part2->Px();
    py = part1->Py()+part2->Py();
    return TMath::Sqrt(px*px+py*py);
}

Float_t AliDimuCombinator::Pz(const TParticle* part1, const TParticle* part2) const
{
// Pz of dimuon system
//
    return part1->Pz()+part2->Pz();
}

Float_t AliDimuCombinator::Y(const TParticle* part1, const TParticle* part2) const
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
void AliDimuCombinator::SmearGauss(Float_t width, Float_t & value) const
{
// Apply gaussian smearing
//
    value+=gRandom->Gaus(0, width);
}
//              Weighting
// 

Float_t AliDimuCombinator::DecayProbability(const TParticle* part) const
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

//Begin_Html
/*
<p> In the the code above :
<P>If h is less than 4 cm, pions or kaons go in the beam pipe and can have a long way
<BR>If h is greater than 4 cm, pions or kaons crash into the front absorber
<P><IMG SRC="absorbeur.jpg" HEIGHT=292 WIDTH=819>
*/
//End_Html


Float_t AliDimuCombinator::Weight(const TParticle* part1, const TParticle* part2) const
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

//Begin_Html
/*
<p>Some clarifications on the calculation of the dimuons weight :
<P>We must keep in mind that if we force the meson decay in muons and we put
lot of mesons (J/psi, upsilon, ...) to have a good statistic we are
obliged to calculate different weights to correct the number
of muons
<BR>&nbsp;
<P>First -->
<BR>The particle weight is given by w=R*M*Br
<BR>&nbsp;with&nbsp; :
<UL>R&nbsp;&nbsp; =&nbsp; the rate by event. This number gives the number
of produced J/psi, upsilon, pion ... in a collision.
<BR>It corresponds of the weight 0.06 given for example in&nbsp; gener->AddGenerator(jpsi,"J/Psi",
0.06); from the config.C macro.
<BR>In this example R=0.06

<P>M&nbsp; = the rate of the mother production. This number depend on :
<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - the number of generated events --> fParentWeight=1./Float_t(fNpart) in AliGenPythia.cxx . This
is a normalization to 1 of the number of generated particles.
<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - the kinematic bias coming
from the y and Pt cuts.&nbsp; Method&nbsp; AliGenPythia::AdjustWeights() in AliGenPythia.cxx
<BR>(in AliGenParam.cxx this 2 things are taken into account in fParentWeight
= fYWgt*fPtWgt*phiWgt/fNpart )

<P>Br = the branching ratio in muon from the mother decay</UL>

<P><BR>In this method, part->GetWeight() = M*Br
<UL>&nbsp;</UL>
Next -->
<BR>The weight of the dimuon depends on the correlation between muons
<BR>&nbsp;
<UL>If the muons are correlated and come from a resonance (for example
J/psi -> mu+ mu-) , the weight of the dimuon is the weight of one muon then
<BR>w12= R*M*Br = w1* R1 (in this method this gives part1->GetWeight()*fRate1)

<P>If the muons are correlated and come from a charm or a bottom pair then
w12 = M*R*Br1*Br2 = w1*w2*R1/M1
<BR>(in this method this gives wgt/(Parent(part1)->GetWeight())*fRate1).
Indeed the 2 muons come from the same mother so the
<BR>weight of a DD~ or BB~ is M*Br and they are no correlation in the decay
(Br1*Br2)

<P>If the muons are not correlated w12 = M1*M2*R1*R2*Br1*Br2 = w1*w2*R1*R2
(in this method this gives wgt*fRate1*fRate2)
<BR>&nbsp;</UL>
*/
//End_Html


Float_t AliDimuCombinator::Weight(const TParticle* part) const
{
// Single muon weight
    return (part->GetWeight())*(Parent(part)->GetWeight())*fRate1;
}

Bool_t  AliDimuCombinator::Correlated(const TParticle* part1, const TParticle* part2) const
{
// Check if muons are correlated
//
    if ((Origin(part1) >= 0) && Origin(part1) == Origin(part2)) {

	return kTRUE;
    } else {
	return kFALSE;
    }
}

TParticle* AliDimuCombinator::Parent(const TParticle* part) const
{
// Return pointer to parent
//
    return Particle(part->GetFirstMother());
}

Int_t AliDimuCombinator::Origin(const TParticle* part) const
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

Int_t AliDimuCombinator::Type(const TParticle *part)  const
{
// Return particle type for 
return part->GetPdgCode();
}

AliDimuCombinator& AliDimuCombinator::operator=(const  AliDimuCombinator& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}


void AliDimuCombinator::Copy(TObject&) const
{
  //
  // Copy 
  //
  Fatal("Copy","Not implemented!\n");
}





