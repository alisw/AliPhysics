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
*/


#include "AliPythia.h"
#include "AliMC.h"
ClassImp(AliPythia)

#ifndef WIN32
# define lu1ent lu1ent_
# define type_of_call
#else
# define lu1ent LU1ENT
# define type_of_call _stdcall
#endif

extern "C" void type_of_call 
          lu1ent(Int_t&, Int_t&, Float_t&, Float_t&, Float_t&);



//_____________________________________________________________________________

Int_t AliPythia::fgInit=0;

AliPythia::AliPythia()
{
    for (Int_t i=0; i<501; i++) {
	fGPCode[i][0]=0; 
	fGPCode[i][1]=0;
    }
}

void  AliPythia::Lu1Ent(Int_t flag, Int_t idpart, 
		      Float_t mom, Float_t theta,Float_t phi)
{
  printf("%d %d %f %f %f\n",flag, idpart, mom, theta, phi);
  lu1ent(flag, idpart, mom, theta, phi);

}
void AliPythia::DecayParticle(Int_t idpart, 
			      Float_t mom, Float_t theta,Float_t phi)
{
    Lu1Ent(0, idpart, mom, theta, phi);
    GetPrimaries();
}

void AliPythia::ProcInit(Process_t process, Float_t energy, StrucFunc_t strucfunc)
{
    fProcess = process;
    fEcms = energy;
    fStrucFunc = strucfunc;
//  don't decay p0
    SetMDCY(LuComp(111),1,0);
//  select structure function 
    SetMSTP(52,2);
    SetMSTP(51,strucfunc);
//
// Pythia initialisation for selected processes//
//
// Make MSEL clean
//
    for (Int_t i=1; i<= 200; i++) {
	SetMSUB(i,0);
    }
//  select charm production
    switch (process) 
    {
    case charm:
	SetMSEL(4);
//
//  heavy quark masses

	SetPMAS(4,1,1.2);

//
//    primordial pT
	SetMSTP(91,1);
	SetPARP(91,1);
	SetPARP(93,3);
//
	break;
    case beauty:
	SetMSEL(5);
	SetPMAS(5,1,4.75);
	break;
    case jpsi:
	SetMSEL(0);
// gg->J/Psi g
	SetMSUB(86,1);
	break;
    case jpsi_chi:
	SetMSEL(0);
// gg->J/Psi g
	SetMSUB(86,1);
// gg-> chi_0c g
	SetMSUB(87,1);
// gg-> chi_1c g
	SetMSUB(88,1);
// gg-> chi_2c g
	SetMSUB(89,1);	
    case charm_unforced:
	SetMSEL(0);
// gq->qg   
	SetMSUB(28,1);
// gg->qq
	SetMSUB(53,1);
// gg->gg
	SetMSUB(68,1);
    case beauty_unforced:
	SetMSEL(0);
// gq->qg   
	SetMSUB(28,1);
// gg->qq
	SetMSUB(53,1);
// gg->gg
	SetMSUB(68,1);
	break;
    case mb:
// Minimum Bias pp-Collisions
//
// Tuning of parameters descibed in G. Ciapetti and A. Di Ciaccio
// Proc. of the LHC Workshop, Aachen 1990, Vol. II p. 155
//   
//      select Pythia min. bias model
	SetMSEL(2);
	SetMSUB(92,1);
	SetMSUB(93,1);
	SetMSUB(94,1);
	SetMSUB(95,1);	
//      Multiple interactions switched on
	SetMSTP(81,1);
	SetMSTP(82,1);
//      Low-pT cut-off for hard scattering
	SetPARP(81,1.9);
//      model for subsequent non-hardest interaction
//      90% gg->gg 10% gg->qq
	SetPARP(86,0.9);
//      90% of gluon interactions have minimum string length
	SetPARP(85,0.9);
    }
//
//  Initialize PYTHIA
    Initialize("CMS","p","p",fEcms);
}

Int_t AliPythia::CountProducts(Int_t channel, Int_t particle)
{
    Int_t np=0;
    for (Int_t i=1; i<=5; i++) {
	if (TMath::Abs(GetKFDP(channel,i)) == particle) np++;
    }
    return np;
}

void AliPythia::AllowAllDecays()
{
    Int_t i;
    for (i=1; i<= 2000; i++) {
	SetMDME(i,1,1);
    }
//
    for (i=0; i<501; i++){
	fBraPart[i]=1;
    }
}

void AliPythia::ForceParticleDecay(Int_t particle, Int_t product, Int_t mult)
{
//
//  force decay of particle into products with multiplicity mult

    Int_t kc=LuComp(particle);
    SetMDCY(kc,1,1);
    Int_t ifirst=GetMDCY(kc,2);
    Int_t ilast=ifirst+GetMDCY(kc,3)-1;
    fBraPart[kc] = 1;
//
//  Loop over decay channels
    for (Int_t channel=ifirst; channel<=ilast;channel++) {
	if (CountProducts(channel,product) >= mult) {
	    SetMDME(channel,1,1);
	} else {
	    SetMDME(channel,1,0);
	    fBraPart[kc]-=GetBRAT(channel);
	}
    }
}

void AliPythia::ForceDecay(Decay_t decay)
{
    fDecay=decay;
//
// Make clean
// AllowAllDecays();
//
// select mode    

    switch (decay) 
    {
    case semimuonic:
	if (fProcess==charm || fProcess == charm_unforced) {
	    ForceParticleDecay(  411,13,1); // D+/-     
	    ForceParticleDecay(  421,13,1); // D0     
	    ForceParticleDecay(  431,13,1); // D_s     
	    ForceParticleDecay( 4122,13,1); // Lambda_c     
	}
	if (fProcess==beauty || fProcess == beauty_unforced) {
	    ForceParticleDecay(  511,13,1); // B0     
	    ForceParticleDecay(  521,13,1); // B+/-     
	    ForceParticleDecay(  531,13,1); // B_s     
	    ForceParticleDecay( 5122,13,1); // Lambda_b    
	}
    break;
    case dimuon:
	ForceParticleDecay(   41,13,2); // phi
	ForceParticleDecay(  443,13,2); // J/Psi
	ForceParticleDecay(30443,13,2); // Psi'
	ForceParticleDecay(  553,13,2); // Upsilon
	ForceParticleDecay(30553,13,2); // Upsilon'
	break;
    case semielectronic:
	
	ForceParticleDecay(  411,11,1); // D+/-     
	ForceParticleDecay(  421,11,1); // D0     
	ForceParticleDecay(  431,11,1); // D_s     
	ForceParticleDecay( 4122,11,1); // Lambda_c     
	
	ForceParticleDecay(  511,11,1); // B0     
	ForceParticleDecay(  521,11,1); // B+/-     
	ForceParticleDecay(  531,11,1); // B_s     
	ForceParticleDecay( 5122,11,1); // Lambda_b     
	break;
    case dielectron:

	ForceParticleDecay(   41,11,2); // phi
	ForceParticleDecay(  443,11,2); // J/Psi
	ForceParticleDecay(30443,11,2); // Psi'
	ForceParticleDecay(  553,11,2); // Upsilon
	ForceParticleDecay(30553,11,2); // Upsilon'
	break;
    case b_jpsi_dimuon:
	ForceParticleDecay(  511,443,1); // B0     
	ForceParticleDecay(  521,443,1); // B+/-     
	ForceParticleDecay(  531,443,1); // B_s     
	ForceParticleDecay( 5122,443,1); // Lambda_b
	ForceParticleDecay(  443,13,2);  // J/Psi    
	break;
    case b_psip_dimuon:
	ForceParticleDecay(  511,30443,1); // B0     
	ForceParticleDecay(  521,30443,1); // B+/-     
	ForceParticleDecay(  531,30443,1); // B_s     
	ForceParticleDecay( 5122,30443,1); // Lambda_b 
	ForceParticleDecay(30443,13,2);    // Psi'   
	break;
    case b_jpsi_dielectron:
	ForceParticleDecay(  511,443,1); // B0     
	ForceParticleDecay(  521,443,1); // B+/-     
	ForceParticleDecay(  531,443,1); // B_s     
	ForceParticleDecay( 5122,443,1); // Lambda_b
	ForceParticleDecay(  443,11,2);  // J/Psi    
	break;
    case b_psip_dielectron:
	ForceParticleDecay(  511,30443,1); // B0     
	ForceParticleDecay(  521,30443,1); // B+/-     
	ForceParticleDecay(  531,30443,1); // B_s     
	ForceParticleDecay( 5122,30443,1); // Lambda_b 
	ForceParticleDecay(30443,11,2);    // Psi'   
	break;
    case pitomu:
	ForceParticleDecay(211,13,1); // pi->mu     
	break;
    case katomu:
	ForceParticleDecay(321,13,1); // K->mu     
	break;
    }
}


    void AliPythia::DefineParticles()
{
    if (fgInit) return;
    fgInit=1;
    
    Float_t mass;
    Float_t tlife;
    Int_t kc, nkc, i;
//
//
// Some particles cloned for rare decays     
//
//  phi-> mu+mu- and phi -> e+e-
//  clone the original phi
    kc  = LuComp(333);
    nkc = 41;
    
    for (i=1;i<=3;i++) {
	SetKCHG(nkc,i,GetKCHG(kc,i));
    }
    
    for (i=1;i<=4;i++) {
    	SetPMAS(nkc,i,GetPMAS(kc,i));
    }
    SetCHAF(nkc,GetCHAF(kc));
    fBraPart[kc]=1;
//    
//  decay
    SetMDCY(nkc,1,1);
    SetMDCY(nkc,2,993);
    SetMDCY(nkc,3,2);
//
//  phi-> e+e-
    SetMDME(993,1,1);
    SetMDME(993,2,0);
    SetBRAT(993,2.99e-4);
    SetKFDP(993,1,+11);
    SetKFDP(993,2,-11);
    SetKFDP(993,3,0);
    SetKFDP(993,4,0);
    SetKFDP(993,5,0);
//
//  phi-> mu+mu-
    SetMDME(994,1,1);
    SetMDME(994,2,0);
    SetBRAT(994,2.5e-4);
    SetKFDP(994,1,+13);
    SetKFDP(994,2,-13);
    SetKFDP(994,3,0);
    SetKFDP(994,4,0);
    SetKFDP(994,5,0);
//
//          Vector mesons
//
// phi clone for dilepton decay-channel
    kc=LuComp(41);	    
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(113,"Phi",3,mass,0,tlife);
    fGPCode[kc][0]=113;
    fGPCode[kc][1]=113;
    // J/Psi  
    kc=LuComp(443);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(114,"J/Psi",3,mass,0,tlife);
    fGPCode[kc][0]=114;
    fGPCode[kc][1]=114;
    // psi prime
    kc=LuComp(30443);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(115,"Psi'",3,mass,0,tlife);
    fGPCode[kc][0]=115;
    fGPCode[kc][1]=115;
    // upsilon(1s) 
    kc=LuComp(553); 	
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(116,"Upsilon",3,mass,0,tlife);
    fGPCode[kc][0]=116;
    fGPCode[kc][1]=116;
    // upsilon(2s)	
    kc=LuComp(30553);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(117,"Upsilon'",3,mass,0,tlife);
    fGPCode[kc][0]=117;
    fGPCode[kc][1]=117;
    // upsilon(3s)	
    kc=LuComp(30553);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(118,"Upsilon''",3,mass,0,tlife);
    fGPCode[kc][0]=118;
    fGPCode[kc][1]=118;
//
// charmed mesons
//
    //  D^+/-
    kc=LuComp(411);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(119,"D^+",3,mass, 1,tlife);
    gMC->Gspart(120,"D^-",3,mass,-1,tlife);
    fGPCode[kc][0]=119;
    fGPCode[kc][1]=120;
    // D^0
    kc=LuComp(421);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(121,"D^0",3,mass,0,tlife);
    gMC->Gspart(122,"D^0bar",3,mass,0,tlife);
    fGPCode[kc][0]=121;
    fGPCode[kc][1]=122;
    // D_s
    kc=LuComp(431);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(123,"D_s^+",3,mass, 1,tlife);
    gMC->Gspart(124,"D_s^-",3,mass,-1,tlife);
    fGPCode[kc][0]=123;
    fGPCode[kc][1]=124;
    // Lambda_c
    kc=LuComp(4122);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(125,"Lambda_c+",3,mass, 1,tlife);
    gMC->Gspart(126,"Lambda_c-",3,mass,-1,tlife);
    fGPCode[kc][0]=125;
    fGPCode[kc][1]=126;
    //
    //  beauty mesons
    // B_0
    kc=LuComp(511);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(127,"B^0",3,mass, 0,tlife);
    gMC->Gspart(128,"B^0bar",3,mass, 0,tlife);
    fGPCode[kc][0]=127;
    fGPCode[kc][1]=128;
    // B^+-
    kc=LuComp(521);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(129,"B^+",3,mass, 1,tlife);
    gMC->Gspart(130,"B^-",3,mass,-1,tlife);
    fGPCode[kc][0]=129;
    fGPCode[kc][1]=130;
    // B_s
    kc=LuComp(531);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(131,"B_s",3,mass, 0,tlife);
    gMC->Gspart(132,"B_s^bar",3,mass,0,tlife);
    fGPCode[kc][0]=131;
    fGPCode[kc][1]=132;
    // Lambda_b
    kc=LuComp(5122);
    mass =GetPMAS(kc,1);
    tlife=GetPMAS(kc,4);
    gMC->Gspart(133,"Lambda_b",3,mass, 0,tlife);
    gMC->Gspart(134,"Lambda_b^bar",3,mass,0,tlife);
    fGPCode[kc][0]=133;
    fGPCode[kc][1]=134;
//
//          set up correspondance between standard GEANT particle codes
//          and PYTHIA kf

    kc=LuComp(22);  // gamma
    fGPCode[kc][0]=1;
    fGPCode[kc][1]=1;
    
    kc=LuComp(11);  // positron
    fGPCode[kc][0]=2;
    fGPCode[kc][1]=3;
    
    kc=LuComp(12);  // neutrino
    fGPCode[kc][0]=4;
    fGPCode[kc][1]=4;

    kc=LuComp(13);  // muon
    fGPCode[kc][0]=5;
    fGPCode[kc][1]=6;
    
    kc=LuComp(111); // pi0
    fGPCode[kc][0]=7;
    fGPCode[kc][1]=7;

    kc=LuComp(211); // pi+
    fGPCode[kc][0]=8;
    fGPCode[kc][1]=9;

    kc=LuComp(130); // K0 short
    fGPCode[kc][0]=10;
    fGPCode[kc][1]=10;

    kc=LuComp(321); // K+/-
    fGPCode[kc][0]=11;
    fGPCode[kc][1]=12;

    kc=LuComp(2112); // neutron/anti-neutron
    fGPCode[kc][0]=13;
    fGPCode[kc][1]=25;
    
    kc=LuComp(2212); // proton/anti-proton
    fGPCode[kc][0]=14;
    fGPCode[kc][1]=15;
    
    kc=LuComp(310);  // K0 short
    fGPCode[kc][0]=16; 
    fGPCode[kc][1]=16;

    kc=LuComp(221);  // eta
    fGPCode[kc][0]=17;
    fGPCode[kc][1]=17;

    kc=LuComp(3122); // lambda
    fGPCode[kc][0]=18;
    fGPCode[kc][1]=18;

    kc=LuComp(3222); // sigma+/antisigma+
    fGPCode[kc][0]=19;
    fGPCode[kc][1]=29;

    kc=LuComp(3212); // sigma0/antisigma0
    fGPCode[kc][0]=20;
    fGPCode[kc][1]=28;

    kc=LuComp(3112); // sigma-/antisigma-
    fGPCode[kc][0]=21;
    fGPCode[kc][1]=27;

    kc=LuComp(3322); // xsi0-/antixsi0
    fGPCode[kc][0]=22;
    fGPCode[kc][1]=30;

    kc=LuComp(3312); // xsi-/antixsi+
    fGPCode[kc][0]=23;
    fGPCode[kc][1]=31;

    kc=LuComp(3334); // omega/antiomega
    fGPCode[kc][0]=24;
    fGPCode[kc][1]=32;
}


    
Int_t  AliPythia::GetGeantCode(Int_t kf)
{
    Int_t kc=LuComp(TMath::Abs(kf));
    return (kf > 0) ? fGPCode[kc][0] : fGPCode[kc][1];
}
    
Float_t  AliPythia::GetBraPart(Int_t kf)
{
    Int_t kc=LuComp(TMath::Abs(kf));
    return fBraPart[kc];
}

    






