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
Revision 1.3  1999/09/29 09:24:08  fca
Introduction of the Copyright and cvs Log

*/

//
//
//
//
#include "AliDimuCombinator.h" 
#include "AliRun.h" 
#include "AliPDG.h" 
#include "TRandom.h" 
//
ClassImp(AliDimuCombinator)
//
//                       Iterators
// 
 TParticle* AliDimuCombinator::FirstMuon()
     {
	 fimuon1=fimin1;
	 fmuon1 = (TParticle*) fPartArray->UncheckedAt(fimuon1);
	 while(Type(fmuon1)!=kMuonPlus && Type(fmuon1)!=kMuonMinus) {
	     fimuon1++;
	     if (fimuon1 >= fimax1) {fmuon1=0; break;}
	     fmuon1 = (TParticle*) fPartArray->UncheckedAt(fimuon1);
	 }
	 return fmuon1;
     }
	     
 TParticle* AliDimuCombinator::FirstMuonSelected()
     {
	 TParticle * muon=FirstMuon();
	 while(muon!=0 && !Selected(muon)) {muon=NextMuon();}
	 return muon;
     }
	     

 TParticle* AliDimuCombinator::NextMuon()
     {
	 fimuon1++;
	 if (fimuon1>=fNParticle) {fmuon1 = 0; return fmuon1;}
	 
	 fmuon1 = (TParticle*) fPartArray->UncheckedAt(fimuon1);
	 while(Type(fmuon1)!=kMuonPlus && Type(fmuon1)!=kMuonMinus) {
	     fimuon1++;
	     if (fimuon1>=fimax1) {fmuon1 = 0; break;}
	     fmuon1 = (TParticle*) fPartArray->UncheckedAt(fimuon1);
	 }
	 return fmuon1;
     }

TParticle* AliDimuCombinator::NextMuonSelected()
     {
	 TParticle * muon=NextMuon();
	 while(muon !=0 && !Selected(muon)) {muon=NextMuon();}
	 return muon;
     }


 void AliDimuCombinator::FirstPartner()
     {
	 if (fimin1==fimin2) {
	     fimuon2=fimuon1+1;
	 } else {
	     fimuon2=fimin2;
	 }
	 if (fimuon2 >= fimax2) {fmuon2=0; return;}
	 fmuon2 = (TParticle*) fPartArray->UncheckedAt(fimuon2);
	 while(Type(fmuon2)!=kMuonPlus && Type(fmuon2)!=kMuonMinus) {
	     fimuon2++;
	     if (fimuon2 >= fimax2) {fmuon2=0; break;}
	     fmuon2 = (TParticle*) fPartArray->UncheckedAt(fimuon2);
	 }
     }
void AliDimuCombinator::FirstPartnerSelected()
{
	 FirstPartner();
	 while(fmuon2 !=0 && !Selected(fmuon2)) {NextPartner();}
}


 void AliDimuCombinator::NextPartner()
     {
	 fimuon2++;
	 if (fimuon2>=fimax2) {fmuon2 = 0; return;}

	 
	 fmuon2 = (TParticle*) fPartArray->UncheckedAt(fimuon2);

	 while(Type(fmuon2)!=kMuonPlus && Type(fmuon2)!=kMuonMinus) {
	     fimuon2++;
	     if (fimuon2>=fimax2) {fmuon2 = 0; break;}
	     fmuon2 = (TParticle*) fPartArray->UncheckedAt(fimuon2);
	 }

     }

void AliDimuCombinator::NextPartnerSelected()
{
	 NextPartner();
	 while(fmuon2 !=0 && !Selected(fmuon2)) {NextPartner();}
}


 TParticle*  AliDimuCombinator::Partner()
     {
	 return fmuon2;
     }

void AliDimuCombinator::FirstMuonPair(TParticle* & muon1, TParticle* & muon2)
     {
	 FirstMuon();
	 FirstPartner();
	 muon1=fmuon1;
	 muon2=fmuon2;	 
     }
void AliDimuCombinator::NextMuonPair(TParticle* & muon1, TParticle* & muon2)
     {
	 NextPartner();
	 if (!Partner()) {
	     NextMuon();
	     FirstPartner();
	 }
	 muon1=fmuon1;
	 muon2=fmuon2;	 
     }
void AliDimuCombinator::FirstMuonPairSelected(TParticle* & muon1, TParticle* & muon2)
     {
	 FirstMuonSelected();
	 FirstPartnerSelected();
	 muon1=fmuon1;
	 muon2=fmuon2;	 
     }
void AliDimuCombinator::NextMuonPairSelected(TParticle* & muon1, TParticle* & muon2)
     {
	 NextPartnerSelected();
	 if (!Partner()) {
	     NextMuonSelected();
	     FirstPartnerSelected();
	 }
	 muon1=fmuon1;
	 muon2=fmuon2;	 
     }
void AliDimuCombinator::ResetRange()
{
    fimin1=fimin2=0;
    fimax1=fimax2=fNParticle;
}

void AliDimuCombinator::SetFirstRange(Int_t from, Int_t to)
{
    fimin1=from;
    fimax1=to;
    if (fimax1 > fNParticle) fimax1=fNParticle;
}

void AliDimuCombinator::SetSecondRange(Int_t from, Int_t to)
{
    fimin2=from;
    fimax2=to;
    if (fimax2 > fNParticle) fimax2=fNParticle;
}
//
//                       Selection
//

Bool_t AliDimuCombinator::Selected(TParticle* part)
{
// 
//
    if (part==0) {return 0;}
    
    if (part->Pt() > fPtMin && part->Eta()>fEtaMin && part->Eta()<fEtaMax) {
	return 1;
    } else {
	return 0;
    }
    
    
}

Bool_t AliDimuCombinator::Selected(TParticle* part1, TParticle* part2)
{
     return Selected(part1)*Selected(part2);
}
//
//                       Kinematics
//
Float_t AliDimuCombinator::Mass(TParticle* part1, TParticle* part2)
{
    Float_t px,py,pz,e;
    px=part1->Px()+part2->Px();
    py=part1->Py()+part2->Py();
    pz=part1->Pz()+part2->Pz();    
    e =part1->Energy()+part2->Energy();
    Float_t p=px*px+py*py+pz*pz;
    if (e*e < p) {
	return -1; 
    } else {
	return TMath::Sqrt(e*e-p);
    }
}

Float_t AliDimuCombinator::PT(TParticle* part1, TParticle* part2)
{
    Float_t px,py;
    px=part1->Px()+part2->Px();
    py=part1->Py()+part2->Py();
    return TMath::Sqrt(px*px+py*py);
}

Float_t AliDimuCombinator::Pz(TParticle* part1, TParticle* part2)
{
    return part1->Pz()+part2->Pz();
}

Float_t AliDimuCombinator::Y(TParticle* part1, TParticle* part2)
{
    Float_t pz,e;
    pz=part1->Pz()+part2->Pz();
    e =part1->Energy()+part2->Energy();
    return 0.5*TMath::Log((e+pz)/(e-pz));
}
//                  Response
//
void AliDimuCombinator::SmearGauss(Float_t width, Float_t & value)
{
    value+=gRandom->Gaus(0, width);
}
//              Weighting
// 

Float_t AliDimuCombinator::Decay_Prob(TParticle* part)
{
    Float_t d, h, theta, CTau;
    TParticle* parent = Parent(part);
    Int_t ipar=Type(parent);
    if (ipar==8 || ipar==9) {
	CTau=780.4;
    } else if (ipar==11 || ipar==12) {
	CTau=370.9;
    } else {
	CTau=0;
    }
    
    
    Float_t GammaBeta=(parent->P())/(parent->GetMass());
//
// this part is still very ALICE muon-arm specific
//
    theta=parent->Theta();
    h=90*TMath::Tan(theta);
    
    if (h<4) {
	d=4/TMath::Sin(theta);
    } else {
	d=90/TMath::Cos(theta);
    }
    
    if (CTau > 0) {
	return 1-TMath::Exp(-d/CTau/GammaBeta);
    } else {
	return 1;
    }
}

Float_t AliDimuCombinator::Weight(TParticle* part1, TParticle* part2)
{
    Float_t wgt=(part1->GetWeight())*(part2->GetWeight());
    
    if (Correlated(part1, part2)) {
	return wgt/(Parent(part1)->GetWeight())*fRate1;
    } else {
	return wgt*fRate1*fRate2;
    }
} 


Float_t AliDimuCombinator::Weight(TParticle* part)
{
    return (part->GetWeight())*(Parent(part)->GetWeight())*fRate1;
}
Bool_t  AliDimuCombinator::Correlated(TParticle* part1, TParticle* part2)
{
    if (Origin(part1) == Origin(part2)) {
	return kTRUE;
    } else {
	return kFALSE;
    }
}

TParticle* AliDimuCombinator::Parent(TParticle* part)
{
    return (TParticle*) (fPartArray->UncheckedAt(part->GetFirstMother()));
}

Int_t AliDimuCombinator::Origin(TParticle* part)
{
    Int_t iparent= part->GetFirstMother();
    if (iparent < 0) return iparent;
    Int_t ip;
    while(1) {
	ip=((TParticle*) fPartArray->UncheckedAt(iparent))->GetFirstMother();
	if (ip < 0) {
	    break;
	} else {
	    iparent=ip;
	}
    }
    return iparent;
}

