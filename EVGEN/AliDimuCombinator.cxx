//
//
//
//
#include "AliDimuCombinator.h" 
#include "AliRun.h" 
#include "TRandom.h" 
//
ClassImp(AliDimuCombinator)
//
//                       Iterators
// 
 GParticle* AliDimuCombinator::FirstMuon()
     {
	 fimuon1=fimin1;
	 fmuon1 = (GParticle*) fPartArray->UncheckedAt(fimuon1);
	 while(Type(fmuon1)!=5 && Type(fmuon1)!=6) {
	     fimuon1++;
	     if (fimuon1 >= fimax1) {fmuon1=0; break;}
	     fmuon1 = (GParticle*) fPartArray->UncheckedAt(fimuon1);
	 }
	 return fmuon1;
     }
	     
 GParticle* AliDimuCombinator::FirstMuonSelected()
     {
	 GParticle * muon=FirstMuon();
	 while(muon!=0 && !Selected(muon)) {muon=NextMuon();}
	 return muon;
     }
	     

 GParticle* AliDimuCombinator::NextMuon()
     {
	 fimuon1++;
	 if (fimuon1>=fNParticle) {fmuon1 = 0; return fmuon1;}
	 
	 fmuon1 = (GParticle*) fPartArray->UncheckedAt(fimuon1);
	 while(Type(fmuon1)!=5 && Type(fmuon1)!=6) {
	     fimuon1++;
	     if (fimuon1>=fimax1) {fmuon1 = 0; break;}
	     fmuon1 = (GParticle*) fPartArray->UncheckedAt(fimuon1);
	 }
	 return fmuon1;
     }

GParticle* AliDimuCombinator::NextMuonSelected()
     {
	 GParticle * muon=NextMuon();
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
	 fmuon2 = (GParticle*) fPartArray->UncheckedAt(fimuon2);
	 while(Type(fmuon2)!=5 && Type(fmuon2)!=6) {
	     fimuon2++;
	     if (fimuon2 >= fimax2) {fmuon2=0; break;}
	     fmuon2 = (GParticle*) fPartArray->UncheckedAt(fimuon2);
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

	 
	 fmuon2 = (GParticle*) fPartArray->UncheckedAt(fimuon2);

	 while(Type(fmuon2)!=5 && Type(fmuon2)!=6) {
	     fimuon2++;
	     if (fimuon2>=fimax2) {fmuon2 = 0; break;}
	     fmuon2 = (GParticle*) fPartArray->UncheckedAt(fimuon2);
	 }

     }

void AliDimuCombinator::NextPartnerSelected()
{
	 NextPartner();
	 while(fmuon2 !=0 && !Selected(fmuon2)) {NextPartner();}
}


 GParticle*  AliDimuCombinator::Partner()
     {
	 return fmuon2;
     }

void AliDimuCombinator::FirstMuonPair(GParticle* & muon1, GParticle* & muon2)
     {
	 FirstMuon();
	 FirstPartner();
	 muon1=fmuon1;
	 muon2=fmuon2;	 
     }
void AliDimuCombinator::NextMuonPair(GParticle* & muon1, GParticle* & muon2)
     {
	 NextPartner();
	 if (!Partner()) {
	     NextMuon();
	     FirstPartner();
	 }
	 muon1=fmuon1;
	 muon2=fmuon2;	 
     }
void AliDimuCombinator::FirstMuonPairSelected(GParticle* & muon1, GParticle* & muon2)
     {
	 FirstMuonSelected();
	 FirstPartnerSelected();
	 muon1=fmuon1;
	 muon2=fmuon2;	 
     }
void AliDimuCombinator::NextMuonPairSelected(GParticle* & muon1, GParticle* & muon2)
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

Bool_t AliDimuCombinator::Selected(GParticle* part)
{
// 
//
    if (part==0) {return 0;}
    
    if (part->GetPT() > fPtMin && part->GetEta()>fEtaMin && part->GetEta()<fEtaMax) {
	return 1;
    } else {
	return 0;
    }
    
    
}

Bool_t AliDimuCombinator::Selected(GParticle* part1, GParticle* part2)
{
     return Selected(part1)*Selected(part2);
}
//
//                       Kinematics
//
Float_t AliDimuCombinator::Mass(GParticle* part1, GParticle* part2)
{
    Float_t px,py,pz,e;
    px=part1->GetPx()+part2->GetPx();
    py=part1->GetPy()+part2->GetPy();
    pz=part1->GetPz()+part2->GetPz();    
    e =part1->GetEnergy()+part2->GetEnergy();
    Float_t p=px*px+py*py+pz*pz;
    if (e*e < p) {
	return -1; 
    } else {
	return TMath::Sqrt(e*e-p);
    }
}

Float_t AliDimuCombinator::PT(GParticle* part1, GParticle* part2)
{
    Float_t px,py;
    px=part1->GetPx()+part2->GetPx();
    py=part1->GetPy()+part2->GetPy();
    return TMath::Sqrt(px*px+py*py);
}

Float_t AliDimuCombinator::Pz(GParticle* part1, GParticle* part2)
{
    return part1->GetPz()+part2->GetPz();
}

Float_t AliDimuCombinator::Y(GParticle* part1, GParticle* part2)
{
    Float_t pz,e;
    pz=part1->GetPz()+part2->GetPz();
    e =part1->GetEnergy()+part2->GetEnergy();
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

Float_t AliDimuCombinator::Decay_Prob(GParticle* part)
{
    Float_t d, h, theta, CTau;
    GParticle* parent = Parent(part);
    Int_t ipar=Type(parent);
    if (ipar==8 || ipar==9) {
	CTau=780.4;
    } else if (ipar==11 || ipar==12) {
	CTau=370.9;
    } else {
	CTau=0;
    }
    
    
    Float_t GammaBeta=(parent->GetMomentum())/(parent->GetMass());
//
// this part is still very ALICE muon-arm specific
//
    theta=parent->GetTheta();
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

Float_t AliDimuCombinator::Weight(GParticle* part1, GParticle* part2)
{
    Float_t wgt=(part1->GetWgt())*(part2->GetWgt());
    
    if (Correlated(part1, part2)) {
	return wgt/(Parent(part1)->GetWgt())*fRate1;
    } else {
	return wgt*fRate1*fRate2;
    }
} 


Float_t AliDimuCombinator::Weight(GParticle* part)
{
    return (part->GetWgt())*(Parent(part)->GetWgt())*fRate1;
}
Bool_t  AliDimuCombinator::Correlated(GParticle* part1, GParticle* part2)
{
    if (Origin(part1) == Origin(part2)) {
	return kTRUE;
    } else {
	return kFALSE;
    }
}
GParticle* AliDimuCombinator::Parent(GParticle* part)
{
    return (GParticle*) (fPartArray->UncheckedAt(part->GetParent()));
}

Int_t AliDimuCombinator::Origin(GParticle* part)
{
    Int_t iparent= part->GetParent();
    if (iparent < 0) return iparent;
    Int_t ip;
    while(1) {
	ip=((GParticle*) fPartArray->UncheckedAt(iparent))->GetParent();
	if (ip < 0) {
	    break;
	} else {
	    iparent=ip;
	}
    }
    return iparent;
}

