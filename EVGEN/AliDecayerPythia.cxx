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
Revision 1.1  2000/09/06 14:23:43  morsch
Realisation of AliDecayer using Pythia6

*/

#include "AliDecayerPythia.h"
#include <TLorentzVector.h>

ClassImp(AliDecayerPythia)

#ifndef WIN32
# define py1ent py1ent_
# define type_of_call
#else
# define lu1ent PY1ENT
# define type_of_call _stdcall
#endif

extern "C" void type_of_call 
          py1ent(Int_t&, Int_t&, Double_t&, Double_t&, Double_t&);


AliDecayerPythia::AliDecayerPythia()
{
// Default Constructor
    fPythia=AliPythia::Instance();
    for (Int_t i=0; i< 501; i++) fBraPart[i]=1;
}

void AliDecayerPythia::Init()
{
// Initialisation
//
// Switch on heavy flavor decays
    Int_t kc, i, j;
    Int_t heavy[8] = {411, 421, 431, 4122, 511, 521, 531, 5122};
    
    for (j=0; j < 8; j++) {
	kc=fPythia->Pycomp(heavy[j]);
	fPythia->SetMDCY(kc,1,1);
	for (i=fPythia->GetMDCY(kc,2);i<fPythia->GetMDCY(kc,2)+fPythia->GetMDCY(kc,3); i++) {
	    fPythia->SetMDME(i,1,1);
	}
    }
    ForceDecay();
}

void AliDecayerPythia::Decay(Int_t idpart, TLorentzVector* p)
{
    Float_t energy = p->Energy();
    Float_t theta  = p->Theta();
    Float_t phi    = p->Phi();
    
    Lu1Ent(0, idpart, energy, theta, phi);
    fPythia->GetPrimaries();
}

void AliDecayerPythia::ForceDecay()
{
// Force a particle decay mode
    Decay_t decay=fDecay;
    
//
// Make clean
// AllowAllDecays();
//
// select mode    

    switch (decay) 
    {
    case semimuonic:
	ForceParticleDecay(  411,13,1); // D+/-     
	ForceParticleDecay(  421,13,1); // D0     
	ForceParticleDecay(  431,13,1); // D_s     
	ForceParticleDecay( 4122,13,1); // Lambda_c     
	ForceParticleDecay(  511,13,1); // B0     
	ForceParticleDecay(  521,13,1); // B+/-     
	ForceParticleDecay(  531,13,1); // B_s     
	ForceParticleDecay( 5122,13,1); // Lambda_b    
    break;
    case dimuon:
//	ForceParticleDecay(   41,13,2); // phi
	ForceParticleDecay(  443,13,2); // J/Psi
	ForceParticleDecay(20443,13,2); // Psi'
	ForceParticleDecay(  553,13,2); // Upsilon
	ForceParticleDecay(20553,13,2); // Upsilon'
	ForceParticleDecay(30553,13,2); // Upsilon''
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
    case hadronicD:
	ForceHadronicD();
    case all:
	break;
    case nodecay:
	break;
    }
}

void  AliDecayerPythia::Lu1Ent(Int_t flag, Int_t idpart, 
		      Double_t mom, Double_t theta, Double_t phi)
{
// Wrap of Pythia lu1ent subroutine
//  printf("%d %d %f %f %f\n",flag, idpart, mom, theta, phi);
  py1ent(flag, idpart, mom, theta, phi);
  
}



Int_t AliDecayerPythia::CountProducts(Int_t channel, Int_t particle)
{
// Count number of decay products
    Int_t np=0;
    for (Int_t i=1; i<=5; i++) {
	if (TMath::Abs(fPythia->GetKFDP(channel,i)) == particle) np++;
    }
    return np;
}


void AliDecayerPythia::ForceHadronicD()
{
// Force golden D decay modes
//

    Int_t channel;
    
// D+ -> K- pi+ pi+ 
    Int_t kc=fPythia->Pycomp(411);
    fPythia->SetMDCY(kc,1,1);
    Int_t ifirst=fPythia->GetMDCY(kc,2);
    Int_t ilast=ifirst+fPythia->GetMDCY(kc,3)-1;
    for (channel=ifirst; channel<=ilast;channel++) {
	if (channel==837) {
	    fPythia->SetMDME(channel,1,1);
	} else {
	    fPythia->SetMDME(channel,1,0);
	    fBraPart[kc]-=fPythia->GetBRAT(channel);
	}
    }

// D0 -> K- pi+
    kc=fPythia->Pycomp(421);
    fPythia->SetMDCY(kc,1,1);
    ifirst=fPythia->GetMDCY(kc,2);
    ilast=ifirst+fPythia->GetMDCY(kc,3)-1;
    for (channel=ifirst; channel<=ilast;channel++) {
	if (channel==881) {
	    fPythia->SetMDME(channel,1,1);
	} else {
	    fPythia->SetMDME(channel,1,0);
	    fBraPart[kc]-=fPythia->GetBRAT(channel);
	}
    }

// no D_s decays
    kc=fPythia->Pycomp(431);
    fPythia->SetMDCY(kc,1,1);
    ifirst=fPythia->GetMDCY(kc,2);
    ilast=ifirst+fPythia->GetMDCY(kc,3)-1;
    for (channel=ifirst; channel<=ilast;channel++) {
	fPythia->SetMDME(channel,1,0);
	fBraPart[kc]-=fPythia->GetBRAT(channel);
    }

// no Lambda_c decays
    kc=fPythia->Pycomp(4122);
    fPythia->SetMDCY(kc,1,1);
    ifirst=fPythia->GetMDCY(kc,2);
    ilast=ifirst+fPythia->GetMDCY(kc,3)-1;
    for (channel=ifirst; channel<=ilast;channel++) {
	fPythia->SetMDME(channel,1,0);
	fBraPart[kc]-=fPythia->GetBRAT(channel);
    }
}

void AliDecayerPythia::ForceParticleDecay(Int_t particle, Int_t product, Int_t mult)
{
//
//  force decay of particle into products with multiplicity mult

    Int_t kc=fPythia->Pycomp(particle);
    fPythia->SetMDCY(kc,1,1);
    Int_t ifirst=fPythia->GetMDCY(kc,2);
    Int_t ilast=ifirst+fPythia->GetMDCY(kc,3)-1;
    fBraPart[kc] = 1;
//
//  Loop over decay channels
    for (Int_t channel=ifirst; channel<=ilast;channel++) {
	if (CountProducts(channel,product) >= mult) {
	    fPythia->SetMDME(channel,1,1);
	} else {
	    fPythia->SetMDME(channel,1,0);
	    fBraPart[kc]-=fPythia->GetBRAT(channel);
	}
    }
}


void AliDecayerPythia::AllowAllDecays()
{
// Reset decay flags
    Int_t i;
    for (i=1; i<= 2000; i++) {
	fPythia->SetMDME(i,1,1);
    }
//
    for (i=0; i<501; i++){
	fBraPart[i]=1;
    }
}


void AliDecayerPythia::DefineParticles()
{
    Float_t mass;
    Float_t tlife;
    Int_t kc, nkc, i;
//
//
// Some particles cloned for rare decays     
//
//  phi-> mu+mu- and phi -> e+e-
//  clone the original phi
    kc  = fPythia->Pycomp(333);
    nkc = 41;
    
    for (i=1;i<=3;i++) {
	 fPythia->SetKCHG(nkc,i,fPythia->GetKCHG(kc,i));
    }
    
    for (i=1;i<=4;i++) {
    	 fPythia->SetPMAS(nkc,i,fPythia->GetPMAS(kc,i));
    }
//    fPythia->SetCHAF(nkc, fPythia->GetCHAF(kc));
    fBraPart[kc]=1;
//    
//  decay
    fPythia-> SetMDCY(nkc,1,1);
    fPythia-> SetMDCY(nkc,2,993);
    fPythia-> SetMDCY(nkc,3,2);
//
//  phi-> e+e-
    fPythia->SetMDME(993,1,1);
    fPythia->SetMDME(993,2,0);
    fPythia->SetBRAT(993,2.99e-4);
    fPythia->SetKFDP(993,1,+11);
    fPythia->SetKFDP(993,2,-11);
    fPythia->SetKFDP(993,3,0);
    fPythia->SetKFDP(993,4,0);
    fPythia->SetKFDP(993,5,0);
//
//  phi-> mu+mu-
    fPythia->SetMDME(994,1,1);
    fPythia->SetMDME(994,2,0);
    fPythia->SetBRAT(994,2.5e-4);
    fPythia->SetKFDP(994,1,+13);
    fPythia->SetKFDP(994,2,-13);
    fPythia->SetKFDP(994,3,0);
    fPythia->SetKFDP(994,4,0);
    fPythia->SetKFDP(994,5,0);
//
//          Vector mesons
//
// phi clone for dilepton decay-channel
    kc =   fPythia->Pycomp(41);	    
    mass = fPythia->GetPMAS(kc,1);
    tlife= fPythia->GetPMAS(kc,4);
//    gMC->Gspart(113,"Phi",3,mass,0,tlife);
}

Float_t  AliDecayerPythia::GetPartialBranchingRatio(Int_t kf)
{
// Get branching ratio
    Int_t kc=fPythia->Pycomp(TMath::Abs(kf));
    return fBraPart[kc];
}

void AliDecayerPythia::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliDecayerPythia.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDecayer::Streamer(R__b);
      (AliPythia::Instance())->Streamer(R__b);
      R__b >> (Int_t&)fDecay;
      R__b.ReadStaticArray(fBraPart);
   } else {
      R__b.WriteVersion(AliDecayerPythia::IsA());
      AliDecayer::Streamer(R__b);
      R__b << fPythia;
      R__b << (Int_t)fDecay;
      R__b.WriteArray(fBraPart, 501);
   }
}
/*

                              Particle/parton data table

        KF     KC    particle        antiparticle      chg  col  anti        mass       width       w-cut     lifetime decay
           IDC on/off ME   Br.rat.    decay products

         1      1    d               dbar               -1    1    1      0.33000     0.00000     0.00000   0.00000E+00    0
             1    1  102    0.000000    g               d                                                               
             2    1  102    0.000000    gamma           d                                                               
             3    1  102    0.000000    Z0              d                                                               
             4    1  102    0.000000    W-              u                                                               
             5    1  102    0.000000    W-              c                                                               
             6    1  102    0.000000    W-              t                                                               
             7   -1  102    0.000000    W-              t'                                                              
             8    1  102    0.000000    h0              d                                                               

         2      2    u               ubar                2    1    1      0.33000     0.00000     0.00000   0.00000E+00    0
             9    1  102    0.000000    g               u                                                               
            10    1  102    0.000000    gamma           u                                                               
            11    1  102    0.000000    Z0              u                                                               
            12    1  102    0.000000    W+              d                                                               
            13    1  102    0.000000    W+              s                                                               
            14    1  102    0.000000    W+              b                                                               
            15   -1  102    0.000000    W+              b'                                                              
            16    1  102    0.000000    h0              u                                                               

         3      3    s               sbar               -1    1    1      0.50000     0.00000     0.00000   0.00000E+00    0
            17    1  102    0.000000    g               s                                                               
            18    1  102    0.000000    gamma           s                                                               
            19    1  102    0.000000    Z0              s                                                               
            20    1  102    0.000000    W-              u                                                               
            21    1  102    0.000000    W-              c                                                               
            22    1  102    0.000000    W-              t                                                               
            23   -1  102    0.000000    W-              t'                                                              
            24    1  102    0.000000    h0              s                                                               

         4      4    c               cbar                2    1    1      1.20000     0.00000     0.00000   0.00000E+00    0
            25    1  102    0.000000    g               c                                                               
            26    1  102    0.000000    gamma           c                                                               
            27    1  102    0.000000    Z0              c                                                               
            28    1  102    0.000000    W+              d                                                               
            29    1  102    0.000000    W+              s                                                               
            30    1  102    0.000000    W+              b                                                               
            31   -1  102    0.000000    W+              b'                                                              
            32    1  102    0.000000    h0              c                                                               

         5      5    b               bbar               -1    1    1      4.80000     0.00000     0.00000   0.00000E+00    0
            33    1  102    0.000000    g               b                                                               
            34    1  102    0.000000    gamma           b                                                               
            35    1  102    0.000000    Z0              b                                                               
            36    1  102    0.000000    W-              u                                                               
            37    1  102    0.000000    W-              c                                                               
            38    1  102    0.000000    W-              t                                                               
            39   -1  102    0.000000    W-              t'                                                              
            40    1  102    0.000000    h0              b                                                               

         6      6    t               tbar                2    1    1    175.00000     1.40244    14.02444   0.00000E+00    1
            41    1  102    0.000000    g               t                                                               
            42    1  102    0.000000    gamma           t                                                               
            43    1  102    0.000000    Z0              t                                                               
            44    1    0    0.000030    W+              d                                                               
            45    1    0    0.001765    W+              s                                                               
            46    1    0    0.998205    W+              b                                                               
            47   -1    0    0.000000    W+              b'                                                              
            48    1  102    0.000000    h0              t                                                               
            49   -1    0    0.000000    H+              b                                                               
            50   -1   53    0.000000    ~chi_10         ~t_1                                                            
            51   -1   53    0.000000    ~chi_20         ~t_1                                                            
            52   -1   53    0.000000    ~chi_30         ~t_1                                                            
            53   -1   53    0.000000    ~chi_40         ~t_1                                                            
            54   -1   53    0.000000    ~g              ~t_1                                                            
            55   -1   53    0.000000    ~gravitino      ~t_1                                                            

         7      7    b'              b'bar              -1    1    1    400.00000     0.00000     0.00000   0.00000E+00    1
            56    1  102    0.000000    g               b'                                                              
            57    1  102    0.000000    gamma           b'                                                              
            58    1  102    0.000000    Z0              b'                                                              
            59    1    0    0.000000    W-              u                                                               
            60    1    0    0.000000    W-              c                                                               
            61    1    0    0.000000    W-              t                                                               
            62    1    0    0.000000    W-              t'                                                              
            63    1  102    0.000000    h0              b'                                                              
            64   -1    0    0.000000    H-              c                                                               
            65   -1    0    0.000000    H-              t                                                               

         8      8    t'              t'bar               2    1    1    400.00000     0.00000     0.00000   0.00000E+00    1
            66    1  102    0.000000    g               t'                                                              
            67    1  102    0.000000    gamma           t'                                                              
            68    1  102    0.000000    Z0              t'                                                              
            69    1    0    0.000000    W+              d                                                               
            70    1    0    0.000000    W+              s                                                               
            71    1    0    0.000000    W+              b                                                               
            72    1    0    0.000000    W+              b'                                                              
            73    1  102    0.000000    h0              t'                                                              
            74   -1    0    0.000000    H+              b                                                               
            75   -1    0    0.000000    H+              b'                                                              

        11     11    e-              e+                 -3    0    1      0.00051     0.00000     0.00000   0.00000E+00    0
            76    1  102    0.000000    gamma           e-                                                              
            77    1  102    0.000000    Z0              e-                                                              
            78    1  102    0.000000    W-              nu_e                                                            
            79    1  102    0.000000    h0              e-                                                              

        12     12    nu_e            nu_ebar             0    0    1      0.00000     0.00000     0.00000   0.00000E+00    0
            80    1  102    0.000000    Z0              nu_e                                                            
            81    1  102    0.000000    W+              e-                                                              

        13     13    mu-             mu+                -3    0    1      0.10566     0.00000     0.00000   6.58654E+05    0
            82    1   42    1.000000    nu_ebar         e-              nu_mu                                           
            83    1  102    0.000000    gamma           mu-                                                             
            84    1  102    0.000000    Z0              mu-                                                             
            85    1  102    0.000000    W-              nu_mu                                                           
            86    1  102    0.000000    h0              mu-                                                             

        14     14    nu_mu           nu_mubar            0    0    1      0.00000     0.00000     0.00000   0.00000E+00    0
            87    1  102    0.000000    Z0              nu_mu                                                           
            88    1  102    0.000000    W+              mu-                                                             

        15     15    tau-            tau+               -3    0    1      1.77700     0.00000     0.00000   8.72000E-02    1
            89    1   42    0.178300    nu_ebar         e-              nu_tau                                          
            90    1   42    0.173500    nu_mubar        mu-             nu_tau                                          
            91    1    0    0.113100    nu_tau          pi-                                                             
            92    1    0    0.249400    nu_tau          rho-                                                            
            93    1   41    0.003000    nu_tau          pi-             pi0                                             
            94    1   41    0.090000    nu_tau          rho-            pi0                                             
            95    1   41    0.002700    nu_tau          pi-             pi0             pi0                             
            96    1   41    0.010000    nu_tau          rho-            pi0             pi0                             
            97    1   41    0.001400    nu_tau          pi-             pi0             pi0             pi0             
            98    1   41    0.001200    nu_tau          rho-            pi0             pi0             pi0             
            99    1   41    0.000250    nu_tau          pi-             K_S0                                            
           100    1   41    0.000250    nu_tau          pi-             K_L0                                            
           101    1    0    0.007100    nu_tau          K-                                                              
           102    1    0    0.012000    nu_tau          K*-                                                             
           103    1   41    0.000400    nu_tau          K-              pi0                                             
           104    1   41    0.000750    nu_tau          K*-             pi0                                             
           105    1   41    0.000060    nu_tau          K*-             pi0             pi0                             
           106    1   41    0.000780    nu_tau          K-              K_S0                                            
           107    1   41    0.000780    nu_tau          K-              K_L0                                            
           108    1   41    0.003400    nu_tau          K-              K+              pi-                             
           109    1   41    0.080000    nu_tau          pi-             rho0                                            
           110    1   41    0.011000    nu_tau          pi-             pi+             pi-                             
           111    1   41    0.019100    nu_tau          pi-             omega                                           
           112    1   41    0.000060    nu_tau          pi-             eta                                             
           113    1   41    0.005000    nu_tau          rho-            rho0                                            
           114    1   41    0.013300    nu_tau          pi-             rho0            pi0                             
           115    1   41    0.006700    nu_tau          rho-            pi+             pi-                             
           116    1   41    0.000500    nu_tau          pi-             pi+             pi-             pi0             
           117    1   41    0.003500    nu_tau          rho-            omega                                           
           118    1   41    0.000600    nu_tau          pi-             omega           pi0                             
           119    1   41    0.001500    nu_tau          rho-            eta                                             
           120    1   41    0.000210    nu_tau          pi-             eta             pi0                             
           121    1   41    0.000200    nu_tau          rho-            rho0            pi0                             
           122    1   41    0.000750    nu_tau          pi-             rho0            rho0                            
           123    1   41    0.000100    nu_tau          pi-             eta             eta                             
           124    1   41    0.000200    nu_tau          pi-             rho0            pi0             pi0             
           125    1   41    0.001100    nu_tau          rho-            rho0            pi0             pi0             
           126    1   41    0.000200    nu_tau          pi-             rho+            rho-                            
           127    1   41    0.000200    nu_tau          pi-             rho+            pi-             pi0             
           128    1   41    0.000200    nu_tau          pi-             rho-            pi+             pi0             
           129    1   41    0.000220    nu_tau          pi-             rho0            rho0            pi0             
           130    1   41    0.000400    nu_tau          K*-             pi0             pi0                             
           131    1   41    0.000100    nu_tau          K-              pi0             pi0             pi0             
           132    1   41    0.002050    nu_tau          pi-             K_S0            pi0                             
           133    1   41    0.002050    nu_tau          pi-             K_L0            pi0                             
           134    1   41    0.000690    nu_tau          K-              K_S0            pi0                             
           135    1   41    0.000690    nu_tau          K-              K_L0            pi0                             
           136    1   41    0.000250    nu_tau          pi-             K_S0            K_S0                            
           137    1   41    0.000510    nu_tau          pi-             K_S0            K_L0                            
           138    1   41    0.000250    nu_tau          pi-             K_L0            K_L0                            
           139    1  102    0.000000    gamma           tau-                                                            
           140    1  102    0.000000    Z0              tau-                                                            
           141    1  102    0.000000    W-              nu_tau                                                          
           142    1  102    0.000000    h0              tau-                                                            

        16     16    nu_tau          nu_taubar           0    0    1      0.00000     0.00000     0.00000   0.00000E+00    0
           143    1  102    0.000000    Z0              nu_tau                                                          
           144    1  102    0.000000    W+              tau-                                                            

        17     17    tau'-           tau'+              -3    0    1    400.00000     0.00000     0.00000   0.00000E+00    1
           145    1  102    0.000000    gamma           tau'-                                                           
           146    1  102    0.000000    Z0              tau'-                                                           
           147    1    0    0.000000    W-              nu'_tau                                                         
           148    1  102    0.000000    h0              tau'-                                                           
           149   -1    0    0.000000    H-              nu'_tau                                                         

        18     18    nu'_tau         nu'_taubar          0    0    1      0.00000     0.00000     0.00000   0.00000E+00    0
           150    1  102    0.000000    Z0              nu'_tau                                                         
           151    1    0    0.000000    W+              tau'-                                                           
           152   -1    0    0.000000    H+              tau'-                                                           

        21     21    g                                   0    2    0      0.00000     0.00000     0.00000   0.00000E+00    0
           153    0  102    0.000000    d               dbar                                                            
           154    0  102    0.000000    u               ubar                                                            
           155    0  102    0.000000    s               sbar                                                            
           156    1  102    0.000000    c               cbar                                                            
           157    0  102    0.000000    b               bbar                                                            
           158    0  102    0.000000    t               tbar                                                            
           159    0  102    0.000000    b'              b'bar                                                           
           160    0  102    0.000000    t'              t'bar                                                           
           161    1  102    0.000000    g               g                                                               

        22     22    gamma                               0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0
           162    0  102    0.000000    d               dbar                                                            
           163    0  102    0.000000    u               ubar                                                            
           164    0  102    0.000000    s               sbar                                                            
           165    1  102    0.000000    c               cbar                                                            
           166    0  102    0.000000    b               bbar                                                            
           167    0  102    0.000000    t               tbar                                                            
           168    0  102    0.000000    b'              b'bar                                                           
           169    0  102    0.000000    t'              t'bar                                                           
           170    0  102    0.000000    e-              e+                                                              
           171    0  102    0.000000    mu-             mu+                                                             
           172    0  102    0.000000    tau-            tau+                                                            
           173    0  102    0.000000    tau'-           tau'+                                                           

        23     23    Z0                                  0    0    0     91.18700     2.47872    24.78720   0.00000E+00    1
           174    1   32    0.153998    d               dbar                                                            
           175    1   32    0.119422    u               ubar                                                            
           176    1   32    0.153988    s               sbar                                                            
           177    1   32    0.119322    c               cbar                                                            
           178    1   32    0.152275    b               bbar                                                            
           179    1   32    0.000000    t               tbar                                                            
           180   -1   32    0.000000    b'              b'bar                                                           
           181   -1   32    0.000000    t'              t'bar                                                           
           182    1    0    0.033568    e-              e+                                                              
           183    1    0    0.066789    nu_e            nu_ebar                                                         
           184    1    0    0.033567    mu-             mu+                                                             
           185    1    0    0.066789    nu_mu           nu_mubar                                                        
           186    1    0    0.033492    tau-            tau+                                                            
           187    1    0    0.066789    nu_tau          nu_taubar                                                       
           188   -1    0    0.000000    tau'-           tau'+                                                           
           189   -1    0    0.000000    nu'_tau         nu'_taubar                                                      

        24     24    W+              W-                  3    0    1     80.33000     2.06856    20.68560   0.00000E+00    1
           190    1   32    0.321379    dbar            u                                                               
           191    1   32    0.016498    dbar            c                                                               
           192    1   32    0.000000    dbar            t                                                               
           193   -1   32    0.000000    dbar            t'                                                              
           194    1   32    0.016502    sbar            u                                                               
           195    1   32    0.320685    sbar            c                                                               
           196    1   32    0.000000    sbar            t                                                               
           197   -1   32    0.000000    sbar            t'                                                              
           198    1   32    0.000010    bbar            u                                                               
           199    1   32    0.000591    bbar            c                                                               
           200    1   32    0.000000    bbar            t                                                               
           201   -1   32    0.000000    bbar            t'                                                              
           202   -1   32    0.000000    b'bar           u                                                               
           203   -1   32    0.000000    b'bar           c                                                               
           204   -1   32    0.000000    b'bar           t                                                               
           205   -1   32    0.000000    b'bar           t'                                                              
           206    1    0    0.108138    e+              nu_e                                                            
           207    1    0    0.108138    mu+             nu_mu                                                           
           208    1    0    0.108059    tau+            nu_tau                                                          
           209   -1    0    0.000000    tau'+           nu'_tau                                                         

        25     25    h0                                  0    0    0     80.00000     0.00240     0.02402   0.00000E+00    1
           210    1   32    0.000001    d               dbar                                                            
           211    1   32    0.000000    u               ubar                                                            
           212    1   32    0.000378    s               sbar                                                            
           213    1   32    0.054441    c               cbar                                                            
           214    1   32    0.853399    b               bbar                                                            
           215    1   32    0.000000    t               tbar                                                            
           216   -1   32    0.000000    b'              b'bar                                                           
           217   -1   32    0.000000    t'              t'bar                                                           
           218    1    0    0.000000    e-              e+                                                              
           219    1    0    0.000241    mu-             mu+                                                             
           220    1    0    0.067867    tau-            tau+                                                            
           221   -1    0    0.000000    tau'-           tau'+                                                           
           222    1    0    0.022178    g               g                                                               
           223    1    0    0.000867    gamma           gamma                                                           
           224    1    0    0.000000    gamma           Z0                                                              
           225    1    0    0.000134    Z0              Z0                                                              
           226    1    0    0.000494    W+              W-                                                              
           227   -1   53    0.000000    ~chi_10         ~chi_10                                                         
           228   -1   53    0.000000    ~chi_20         ~chi_10                                                         
           229   -1   53    0.000000    ~chi_20         ~chi_20                                                         
           230   -1   53    0.000000    ~chi_30         ~chi_10                                                         
           231   -1   53    0.000000    ~chi_30         ~chi_20                                                         
           232   -1   53    0.000000    ~chi_30         ~chi_30                                                         
           233   -1   53    0.000000    ~chi_40         ~chi_10                                                         
           234   -1   53    0.000000    ~chi_40         ~chi_20                                                         
           235   -1   53    0.000000    ~chi_40         ~chi_30                                                         
           236   -1   53    0.000000    ~chi_40         ~chi_40                                                         
           237   -1   53    0.000000    ~chi_1+         ~chi_1-                                                         
           238   -1   53    0.000000    ~chi_1+         ~chi_2-                                                         
           239   -1   53    0.000000    ~chi_2+         ~chi_1-                                                         
           240   -1   53    0.000000    ~chi_2+         ~chi_2-                                                         
           241   -1   53    0.000000    ~d_L            ~d_Lbar                                                         
           242   -1   53    0.000000    ~d_R            ~d_Rbar                                                         
           243   -1   53    0.000000    ~d_L            ~d_Rbar                                                         
           244   -1   53    0.000000    ~d_Lbar         ~d_R                                                            
           245   -1   53    0.000000    ~u_L            ~u_Lbar                                                         
           246   -1   53    0.000000    ~u_R            ~u_Rbar                                                         
           247   -1   53    0.000000    ~u_L            ~u_Rbar                                                         
           248   -1   53    0.000000    ~u_Lbar         ~u_R                                                            
           249   -1   53    0.000000    ~s_L            ~s_Lbar                                                         
           250   -1   53    0.000000    ~s_R            ~s_Rbar                                                         
           251   -1   53    0.000000    ~s_L            ~s_Rbar                                                         
           252   -1   53    0.000000    ~s_Lbar         ~s_R                                                            
           253   -1   53    0.000000    ~c_L            ~c_Lbar                                                         
           254   -1   53    0.000000    ~c_R            ~c_Rbar                                                         
           255   -1   53    0.000000    ~c_L            ~c_Rbar                                                         
           256   -1   53    0.000000    ~c_Lbar         ~c_R                                                            
           257   -1   53    0.000000    ~b_1            ~b_1bar                                                         
           258   -1   53    0.000000    ~b_2            ~b_2bar                                                         
           259   -1   53    0.000000    ~b_1            ~b_2bar                                                         
           260   -1   53    0.000000    ~b_1bar         ~b_2                                                            
           261   -1   53    0.000000    ~t_1            ~t_1bar                                                         
           262   -1   53    0.000000    ~t_2            ~t_2bar                                                         
           263   -1   53    0.000000    ~t_1            ~t_2bar                                                         
           264   -1   53    0.000000    ~t_1bar         ~t_2                                                            
           265   -1   53    0.000000    ~e_L-           ~e_L+                                                           
           266   -1   53    0.000000    ~e_R-           ~e_R+                                                           
           267   -1   53    0.000000    ~e_L-           ~e_R+                                                           
           268   -1   53    0.000000    ~e_L+           ~e_R-                                                           
           269   -1   53    0.000000    ~nu_eL          ~nu_eLbar                                                       
           270   -1   53    0.000000    ~nu_eR          ~nu_eRbar                                                       
           271   -1   53    0.000000    ~nu_eL          ~nu_eRbar                                                       
           272   -1   53    0.000000    ~nu_eLbar       ~nu_eR                                                          
           273   -1   53    0.000000    ~mu_L-          ~mu_L+                                                          
           274   -1   53    0.000000    ~mu_R-          ~mu_R+                                                          
           275   -1   53    0.000000    ~mu_L-          ~mu_R+                                                          
           276   -1   53    0.000000    ~mu_L+          ~mu_R-                                                          
           277   -1   53    0.000000    ~nu_muL         ~nu_muLbar                                                      
           278   -1   53    0.000000    ~nu_muR         ~nu_muRbar                                                      
           279   -1   53    0.000000    ~nu_muL         ~nu_muRbar                                                      
           280   -1   53    0.000000    ~nu_muLbar      ~nu_muR                                                         
           281   -1   53    0.000000    ~tau_1-         ~tau_1+                                                         
           282   -1   53    0.000000    ~tau_2-         ~tau_2+                                                         
           283   -1   53    0.000000    ~tau_1-         ~tau_2+                                                         
           284   -1   53    0.000000    ~tau_1+         ~tau_2-                                                         
           285   -1   53    0.000000    ~nu_tauL        ~nu_tauLbar                                                     
           286   -1   53    0.000000    ~nu_tauR        ~nu_tauRbar                                                     
           287   -1   53    0.000000    ~nu_tauL        ~nu_tauRbar                                                     
           288   -1   53    0.000000    ~nu_tauLbar     ~nu_tauR                                                        

        28     28    reggeon                             0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        29     29    pomeron                             0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        32     32    Z'0                                 0    0    0    500.00000    14.54208   145.42084   0.00000E+00    1
           289    1   32    0.145842    d               dbar                                                            
           290    1   32    0.113282    u               ubar                                                            
           291    1   32    0.145842    s               sbar                                                            
           292    1   32    0.113278    c               cbar                                                            
           293    1   32    0.145788    b               bbar                                                            
           294    1   32    0.049004    t               tbar                                                            
           295   -1   32    0.000000    b'              b'bar                                                           
           296   -1   32    0.000000    t'              t'bar                                                           
           297    1    0    0.032021    e-              e+                                                              
           298    1    0    0.063634    nu_e            nu_ebar                                                         
           299    1    0    0.032021    mu-             mu+                                                             
           300    1    0    0.063634    nu_mu           nu_mubar                                                        
           301    1    0    0.032018    tau-            tau+                                                            
           302    1    0    0.063634    nu_tau          nu_taubar                                                       
           303   -1    0    0.000000    tau'-           tau'+                                                           
           304   -1    0    0.000000    nu'_tau         nu'_taubar                                                      
           305   -1    0    0.000000    W+              W-                                                              
           306   -1    0    0.000000    H+              H-                                                              
           307   -1    0    0.000000    Z0              gamma                                                           
           308   -1    0    0.000000    Z0              h0                                                              
           309   -1    0    0.000000    h0              A0                                                              
           310   -1    0    0.000000    H0              A0                                                              

        33     33    Z"0                                 0    0    0    900.00000     0.00000     0.00000   0.00000E+00    0

        34     34    W'+             W'-                 3    0    1    500.00000    16.66312   166.63122   0.00000E+00    1
           311    1   32    0.251235    dbar            u                                                               
           312    1   32    0.012901    dbar            c                                                               
           313    1   32    0.000006    dbar            t                                                               
           314   -1   32    0.000000    dbar            t'                                                              
           315    1   32    0.012901    sbar            u                                                               
           316    1   32    0.250776    sbar            c                                                               
           317    1   32    0.000380    sbar            t                                                               
           318   -1   32    0.000000    sbar            t'                                                              
           319    1   32    0.000008    bbar            u                                                               
           320    1   32    0.000465    bbar            c                                                               
           321    1   32    0.215427    bbar            t                                                               
           322   -1   32    0.000000    bbar            t'                                                              
           323   -1   32    0.000000    b'bar           u                                                               
           324   -1   32    0.000000    b'bar           c                                                               
           325   -1   32    0.000000    b'bar           t                                                               
           326   -1   32    0.000000    b'bar           t'                                                              
           327    1    0    0.085301    e+              nu_e                                                            
           328    1    0    0.085301    mu+             nu_mu                                                           
           329    1    0    0.085299    tau+            nu_tau                                                          
           330   -1    0    0.000000    tau'+           nu'_tau                                                         
           331   -1    0    0.000000    W+              Z0                                                              
           332   -1    0    0.000000    W+              gamma                                                           
           333   -1    0    0.000000    W+              h0                                                              

        35     35    H0                                  0    0    0    300.00000     8.42840    84.28402   0.00000E+00    1
           334    1   32    0.000000    d               dbar                                                            
           335    1   32    0.000000    u               ubar                                                            
           336    1   32    0.000000    s               sbar                                                            
           337    1   32    0.000048    c               cbar                                                            
           338    1   32    0.000768    b               bbar                                                            
           339    1   32    0.000000    t               tbar                                                            
           340   -1   32    0.000000    b'              b'bar                                                           
           341   -1   32    0.000000    t'              t'bar                                                           
           342    1    0    0.000000    e-              e+                                                              
           343    1    0    0.000000    mu-             mu+                                                             
           344    1    0    0.000074    tau-            tau+                                                            
           345   -1    0    0.000000    tau'-           tau'+                                                           
           346    1    0    0.000422    g               g                                                               
           347    1    0    0.000015    gamma           gamma                                                           
           348    1    0    0.000061    gamma           Z0                                                              
           349    1    0    0.306171    Z0              Z0                                                              
           350    1    0    0.688641    W+              W-                                                              
           351    1    0    0.000000    Z0              h0                                                              
           352    1    0    0.003799    h0              h0                                                              
           353    1    0    0.000000    A0              A0                                                              
           354   -1   53    0.000000    ~chi_10         ~chi_10                                                         
           355   -1   53    0.000000    ~chi_20         ~chi_10                                                         
           356   -1   53    0.000000    ~chi_20         ~chi_20                                                         
           357   -1   53    0.000000    ~chi_30         ~chi_10                                                         
           358   -1   53    0.000000    ~chi_30         ~chi_20                                                         
           359   -1   53    0.000000    ~chi_30         ~chi_30                                                         
           360   -1   53    0.000000    ~chi_40         ~chi_10                                                         
           361   -1   53    0.000000    ~chi_40         ~chi_20                                                         
           362   -1   53    0.000000    ~chi_40         ~chi_30                                                         
           363   -1   53    0.000000    ~chi_40         ~chi_40                                                         
           364   -1   53    0.000000    ~chi_1+         ~chi_1-                                                         
           365   -1   53    0.000000    ~chi_1+         ~chi_2-                                                         
           366   -1   53    0.000000    ~chi_2+         ~chi_1-                                                         
           367   -1   53    0.000000    ~chi_2+         ~chi_2-                                                         
           368   -1   53    0.000000    ~d_L            ~d_Lbar                                                         
           369   -1   53    0.000000    ~d_R            ~d_Rbar                                                         
           370   -1   53    0.000000    ~d_L            ~d_Rbar                                                         
           371   -1   53    0.000000    ~d_Lbar         ~d_R                                                            
           372   -1   53    0.000000    ~u_L            ~u_Lbar                                                         
           373   -1   53    0.000000    ~u_R            ~u_Rbar                                                         
           374   -1   53    0.000000    ~u_L            ~u_Rbar                                                         
           375   -1   53    0.000000    ~u_Lbar         ~u_R                                                            
           376   -1   53    0.000000    ~s_L            ~s_Lbar                                                         
           377   -1   53    0.000000    ~s_R            ~s_Rbar                                                         
           378   -1   53    0.000000    ~s_L            ~s_Rbar                                                         
           379   -1   53    0.000000    ~s_Lbar         ~s_R                                                            
           380   -1   53    0.000000    ~c_L            ~c_Lbar                                                         
           381   -1   53    0.000000    ~c_R            ~c_Rbar                                                         
           382   -1   53    0.000000    ~c_L            ~c_Rbar                                                         
           383   -1   53    0.000000    ~c_Lbar         ~c_R                                                            
           384   -1   53    0.000000    ~b_1            ~b_1bar                                                         
           385   -1   53    0.000000    ~b_2            ~b_2bar                                                         
           386   -1   53    0.000000    ~b_1            ~b_2bar                                                         
           387   -1   53    0.000000    ~b_1bar         ~b_2                                                            
           388   -1   53    0.000000    ~t_1            ~t_1bar                                                         
           389   -1   53    0.000000    ~t_2            ~t_2bar                                                         
           390   -1   53    0.000000    ~t_1            ~t_2bar                                                         
           391   -1   53    0.000000    ~t_1bar         ~t_2                                                            
           392   -1   53    0.000000    ~e_L-           ~e_L+                                                           
           393   -1   53    0.000000    ~e_R-           ~e_R+                                                           
           394   -1   53    0.000000    ~e_L-           ~e_R+                                                           
           395   -1   53    0.000000    ~e_L+           ~e_R-                                                           
           396   -1   53    0.000000    ~nu_eL          ~nu_eLbar                                                       
           397   -1   53    0.000000    ~nu_eR          ~nu_eRbar                                                       
           398   -1   53    0.000000    ~nu_eL          ~nu_eRbar                                                       
           399   -1   53    0.000000    ~nu_eLbar       ~nu_eR                                                          
           400   -1   53    0.000000    ~mu_L-          ~mu_L+                                                          
           401   -1   53    0.000000    ~mu_R-          ~mu_R+                                                          
           402   -1   53    0.000000    ~mu_L-          ~mu_R+                                                          
           403   -1   53    0.000000    ~mu_L+          ~mu_R-                                                          
           404   -1   53    0.000000    ~nu_muL         ~nu_muLbar                                                      
           405   -1   53    0.000000    ~nu_muR         ~nu_muRbar                                                      
           406   -1   53    0.000000    ~nu_muL         ~nu_muRbar                                                      
           407   -1   53    0.000000    ~nu_muLbar      ~nu_muR                                                         
           408   -1   53    0.000000    ~tau_1-         ~tau_1+                                                         
           409   -1   53    0.000000    ~tau_2-         ~tau_2+                                                         
           410   -1   53    0.000000    ~tau_1-         ~tau_2+                                                         
           411   -1   53    0.000000    ~tau_1+         ~tau_2-                                                         
           412   -1   53    0.000000    ~nu_tauL        ~nu_tauLbar                                                     
           413   -1   53    0.000000    ~nu_tauR        ~nu_tauRbar                                                     
           414   -1   53    0.000000    ~nu_tauL        ~nu_tauRbar                                                     
           415   -1   53    0.000000    ~nu_tauLbar     ~nu_tauR                                                        

        36     36    A0                                  0    0    0    300.00000     4.91995    49.19946   0.00000E+00    1
           416    1   32    0.000000    d               dbar                                                            
           417    1   32    0.000000    u               ubar                                                            
           418    1   32    0.000001    s               sbar                                                            
           419    1   32    0.000082    c               cbar                                                            
           420    1   32    0.001318    b               bbar                                                            
           421    1   32    0.000000    t               tbar                                                            
           422   -1   32    0.000000    b'              b'bar                                                           
           423   -1   32    0.000000    t'              t'bar                                                           
           424    1    0    0.000000    e-              e+                                                              
           425    1    0    0.000000    mu-             mu+                                                             
           426    1    0    0.000126    tau-            tau+                                                            
           427   -1    0    0.000000    tau'-           tau'+                                                           
           428    1    0    0.002164    g               g                                                               
           429    1    0    0.000010    gamma           gamma                                                           
           430    1    0    0.000002    gamma           Z0                                                              
           431    1    0    0.000000    Z0              Z0                                                              
           432    1    0    0.000000    W+              W-                                                              
           433    1    0    0.996296    Z0              h0                                                              
           434   -1   53    0.000000    ~chi_10         ~chi_10                                                         
           435   -1   53    0.000000    ~chi_20         ~chi_10                                                         
           436   -1   53    0.000000    ~chi_20         ~chi_20                                                         
           437   -1   53    0.000000    ~chi_30         ~chi_10                                                         
           438   -1   53    0.000000    ~chi_30         ~chi_20                                                         
           439   -1   53    0.000000    ~chi_30         ~chi_30                                                         
           440   -1   53    0.000000    ~chi_40         ~chi_10                                                         
           441   -1   53    0.000000    ~chi_40         ~chi_20                                                         
           442   -1   53    0.000000    ~chi_40         ~chi_30                                                         
           443   -1   53    0.000000    ~chi_40         ~chi_40                                                         
           444   -1   53    0.000000    ~chi_1+         ~chi_1-                                                         
           445   -1   53    0.000000    ~chi_1+         ~chi_2-                                                         
           446   -1   53    0.000000    ~chi_2+         ~chi_1-                                                         
           447   -1   53    0.000000    ~chi_2+         ~chi_2-                                                         
           448   -1   53    0.000000    ~d_L            ~d_Lbar                                                         
           449   -1   53    0.000000    ~d_R            ~d_Rbar                                                         
           450   -1   53    0.000000    ~d_L            ~d_Rbar                                                         
           451   -1   53    0.000000    ~d_Lbar         ~d_R                                                            
           452   -1   53    0.000000    ~u_L            ~u_Lbar                                                         
           453   -1   53    0.000000    ~u_R            ~u_Rbar                                                         
           454   -1   53    0.000000    ~u_L            ~u_Rbar                                                         
           455   -1   53    0.000000    ~u_Lbar         ~u_R                                                            
           456   -1   53    0.000000    ~s_L            ~s_Lbar                                                         
           457   -1   53    0.000000    ~s_R            ~s_Rbar                                                         
           458   -1   53    0.000000    ~s_L            ~s_Rbar                                                         
           459   -1   53    0.000000    ~s_Lbar         ~s_R                                                            
           460   -1   53    0.000000    ~c_L            ~c_Lbar                                                         
           461   -1   53    0.000000    ~c_R            ~c_Rbar                                                         
           462   -1   53    0.000000    ~c_L            ~c_Rbar                                                         
           463   -1   53    0.000000    ~c_Lbar         ~c_R                                                            
           464   -1   53    0.000000    ~b_1            ~b_1bar                                                         
           465   -1   53    0.000000    ~b_2            ~b_2bar                                                         
           466   -1   53    0.000000    ~b_1            ~b_2bar                                                         
           467   -1   53    0.000000    ~b_1bar         ~b_2                                                            
           468   -1   53    0.000000    ~t_1            ~t_1bar                                                         
           469   -1   53    0.000000    ~t_2            ~t_2bar                                                         
           470   -1   53    0.000000    ~t_1            ~t_2bar                                                         
           471   -1   53    0.000000    ~t_1bar         ~t_2                                                            
           472   -1   53    0.000000    ~e_L-           ~e_L+                                                           
           473   -1   53    0.000000    ~e_R-           ~e_R+                                                           
           474   -1   53    0.000000    ~e_L-           ~e_R+                                                           
           475   -1   53    0.000000    ~e_L+           ~e_R-                                                           
           476   -1   53    0.000000    ~nu_eL          ~nu_eLbar                                                       
           477   -1   53    0.000000    ~nu_eR          ~nu_eRbar                                                       
           478   -1   53    0.000000    ~nu_eL          ~nu_eRbar                                                       
           479   -1   53    0.000000    ~nu_eLbar       ~nu_eR                                                          
           480   -1   53    0.000000    ~mu_L-          ~mu_L+                                                          
           481   -1   53    0.000000    ~mu_R-          ~mu_R+                                                          
           482   -1   53    0.000000    ~mu_L-          ~mu_R+                                                          
           483   -1   53    0.000000    ~mu_L+          ~mu_R-                                                          
           484   -1   53    0.000000    ~nu_muL         ~nu_muLbar                                                      
           485   -1   53    0.000000    ~nu_muR         ~nu_muRbar                                                      
           486   -1   53    0.000000    ~nu_muL         ~nu_muRbar                                                      
           487   -1   53    0.000000    ~nu_muLbar      ~nu_muR                                                         
           488   -1   53    0.000000    ~tau_1-         ~tau_1+                                                         
           489   -1   53    0.000000    ~tau_2-         ~tau_2+                                                         
           490   -1   53    0.000000    ~tau_1-         ~tau_2+                                                         
           491   -1   53    0.000000    ~tau_1+         ~tau_2-                                                         
           492   -1   53    0.000000    ~nu_tauL        ~nu_tauLbar                                                     
           493   -1   53    0.000000    ~nu_tauR        ~nu_tauRbar                                                     
           494   -1   53    0.000000    ~nu_tauL        ~nu_tauRbar                                                     
           495   -1   53    0.000000    ~nu_tauLbar     ~nu_tauR                                                        

        37     37    H+              H-                  3    0    1    300.00000     5.76067    57.60673   0.00000E+00    1
           496    1   32    0.000000    dbar            u                                                               
           497    1   32    0.000015    sbar            c                                                               
           498    1   32    0.067644    bbar            t                                                               
           499   -1   32    0.000000    b'bar           t'                                                              
           500    1    0    0.000000    e+              nu_e                                                            
           501    1    0    0.000010    mu+             nu_mu                                                           
           502    1    0    0.002701    tau+            nu_tau                                                          
           503   -1    0    0.000000    tau'+           nu'_tau                                                         
           504    1    0    0.929631    W+              h0                                                              
           505   -1   53    0.000000    ~chi_10         ~chi_1+                                                         
           506   -1   53    0.000000    ~chi_10         ~chi_2+                                                         
           507   -1   53    0.000000    ~chi_20         ~chi_1+                                                         
           508   -1   53    0.000000    ~chi_20         ~chi_2+                                                         
           509   -1   53    0.000000    ~chi_30         ~chi_1+                                                         
           510   -1   53    0.000000    ~chi_30         ~chi_2+                                                         
           511   -1   53    0.000000    ~chi_40         ~chi_1+                                                         
           512   -1   53    0.000000    ~chi_40         ~chi_2+                                                         
           513   -1   53    0.000000    ~t_1            ~b_1bar                                                         
           514   -1   53    0.000000    ~t_2            ~b_1bar                                                         
           515   -1   53    0.000000    ~t_1            ~b_2bar                                                         
           516   -1   53    0.000000    ~t_2            ~b_2bar                                                         
           517   -1   53    0.000000    ~d_Lbar         ~u_L                                                            
           518   -1   53    0.000000    ~s_Lbar         ~c_L                                                            
           519   -1   53    0.000000    ~e_L+           ~nu_eL                                                          
           520   -1   53    0.000000    ~mu_L+          ~nu_muL                                                         
           521   -1   53    0.000000    ~tau_1+         ~nu_tauL                                                        
           522   -1   53    0.000000    ~tau_2+         ~nu_tauL                                                        

        38     38    eta_tech0                           0    2    0    350.00000     0.09572     0.95720   0.00000E+00    1
           523    1   32    0.442959    b               bbar                                                            
           524    1   32    0.000000    t               tbar                                                            
           525    1   32    0.557041    g               g                                                               

        39     39    LQ_ue           LQ_uebar           -1    1    1    200.00000     0.39162     3.91621   0.00000E+00    1
           526    1    0    1.000000    u               e-                                                              

        40     40    R0              Rbar0               0    0    1   5000.00000   417.32877  4173.28769   0.00000E+00    1
           527    1   32    0.215122    d               sbar                                                            
           528    1   32    0.215122    u               cbar                                                            
           529    1   32    0.215122    s               bbar                                                            
           530    1   32    0.214727    c               tbar                                                            
           531   -1   32    0.000000    b               b'bar                                                           
           532   -1   32    0.000000    t               t'bar                                                           
           533    1    0    0.069953    e-              mu+                                                             
           534    1    0    0.069953    mu-             tau+                                                            
           535   -1    0    0.000000    tau-            tau'+                                                           

        51     51    pi_tech0                            0    0    0    110.00000     0.02899     0.28994   0.00000E+00    1
           536    1   32    0.017504    s               sbar                                                            
           537    1   32    0.053796    c               cbar                                                            
           538    1   32    0.857596    b               bbar                                                            
           539    1   32    0.000000    t               tbar                                                            
           540    1    0    0.000000    e-              e+                                                              
           541    1    0    0.000251    mu-             mu+                                                             
           542    1    0    0.070854    tau-            tau+                                                            
           543    1   32    0.000000    g               g                                                               

        52     52    pi_tech+        pi_tech-            3    0    1    110.00000     0.01070     0.10704   0.00000E+00    1
           544    1   32    0.042758    c               sbar                                                            
           545    1   32    0.909078    c               bbar                                                            
           546    1   32    0.000000    W+              b               bbar                                            
           547    1    0    0.000000    e+              nu_e                                                            
           548    1    0    0.000170    mu+             nu_mu                                                           
           549    1    0    0.047994    tau+            nu_tau                                                          

        53     53    pi'_tech0                           0    0    0    110.00000     0.04547     0.45469   0.00000E+00    1
           550    1   32    0.011162    s               sbar                                                            
           551    1   32    0.034304    c               cbar                                                            
           552    1   32    0.546865    b               bbar                                                            
           553    1   32    0.000000    t               tbar                                                            
           554    1    0    0.000000    e-              e+                                                              
           555    1    0    0.000160    mu-             mu+                                                             
           556    1    0    0.045181    tau-            tau+                                                            
           557    1   32    0.362328    g               g                                                               

        54     54    rho_tech0                           0    0    0    210.00000     0.87415     8.74152   0.00000E+00    1
           558    1    0    0.144048    W+              W-                                                              
           559    1    0    0.352384    W+              pi_tech-                                                        
           560    1    0    0.352384    pi_tech+        W-                                                              
           561    1    0    0.000000    pi_tech+        pi_tech-                                                        
           562    1    0    0.081586    gamma           pi_tech0                                                        
           563    1    0    0.029378    gamma           pi'_tech0                                                       
           564    1    0    0.001501    Z0              pi_tech0                                                        
           565    1    0    0.000721    Z0              pi'_tech0                                                       
           566    1   32    0.004490    d               dbar                                                            
           567    1   32    0.006482    u               ubar                                                            
           568    1   32    0.004490    s               sbar                                                            
           569    1   32    0.006482    c               cbar                                                            
           570    1   32    0.004485    b               bbar                                                            
           571    1   32    0.000000    t               tbar                                                            
           572   -1   32    0.000000    b'              b'bar                                                           
           573   -1   32    0.000000    t'              t'bar                                                           
           574    1    0    0.002889    e-              e+                                                              
           575    1    0    0.000967    nu_e            nu_ebar                                                         
           576    1    0    0.002889    mu-             mu+                                                             
           577    1    0    0.000967    nu_mu           nu_mubar                                                        
           578    1    0    0.002889    tau-            tau+                                                            
           579    1    0    0.000967    nu_tau          nu_taubar                                                       
           580   -1    0    0.000000    tau'-           tau'+                                                           
           581   -1    0    0.000000    nu'_tau         nu'_taubar                                                      

        55     55    rho_tech+       rho_tech-           3    0    1    210.00000     0.62673     6.26729   0.00000E+00    1
           582    1    0    0.143941    W+              Z0                                                              
           583    1    0    0.491500    W+              pi_tech0                                                        
           584    1    0    0.194259    pi_tech+        Z0                                                              
           585    1    0    0.000000    pi_tech+        pi_tech0                                                        
           586    1    0    0.113795    pi_tech+        gamma                                                           
           587    1    0    0.008460    W+              pi'_tech0                                                       
           588    1   32    0.014790    dbar            u                                                               
           589    1   32    0.000759    dbar            c                                                               
           590    1   32    0.000000    dbar            t                                                               
           591   -1   32    0.000000    dbar            t'                                                              
           592    1   32    0.000759    sbar            u                                                               
           593    1   32    0.014762    sbar            c                                                               
           594    1   32    0.000003    sbar            t                                                               
           595   -1   32    0.000000    sbar            t'                                                              
           596    1   32    0.000000    bbar            u                                                               
           597    1   32    0.000027    bbar            c                                                               
           598    1   32    0.001934    bbar            t                                                               
           599   -1   32    0.000000    bbar            t'                                                              
           600   -1   32    0.000000    b'bar           u                                                               
           601   -1   32    0.000000    b'bar           c                                                               
           602   -1   32    0.000000    b'bar           t                                                               
           603   -1   32    0.000000    b'bar           t'                                                              
           604    1    0    0.005003    e+              nu_e                                                            
           605    1    0    0.005003    mu+             nu_mu                                                           
           606    1    0    0.005002    tau+            nu_tau                                                          
           607   -1    0    0.000000    tau'+           nu'_tau                                                         

        56     56    omega_tech                          0    0    0    210.00000     0.19204     1.92039   0.00000E+00    1
           608    1    0    0.133696    gamma           pi_tech0                                                        
           609    1    0    0.003283    Z0              pi_tech0                                                        
           610    1    0    0.371467    gamma           pi'_tech0                                                       
           611    1    0    0.006835    Z0              pi'_tech0                                                       
           612    1    0    0.031199    W+              pi_tech-                                                        
           613    1    0    0.031199    pi_tech+        W-                                                              
           614    1    0    0.001639    W+              W-                                                              
           615    1    0    0.000000    pi_tech+        pi_tech-                                                        
           616    1   32    0.047205    d               dbar                                                            
           617    1   32    0.073708    u               ubar                                                            
           618    1   32    0.047205    s               sbar                                                            
           619    1   32    0.073705    c               cbar                                                            
           620    1   32    0.047161    b               bbar                                                            
           621    1   32    0.000000    t               tbar                                                            
           622   -1   32    0.000000    b'              b'bar                                                           
           623   -1   32    0.000000    t'              t'bar                                                           
           624    1    0    0.034740    e-              e+                                                              
           625    1    0    0.009160    nu_e            nu_ebar                                                         
           626    1    0    0.034740    mu-             mu+                                                             
           627    1    0    0.009160    nu_mu           nu_mubar                                                        
           628    1    0    0.034738    tau-            tau+                                                            
           629    1    0    0.009160    nu_tau          nu_taubar                                                       
           630   -1    0    0.000000    tau'-           tau'+                                                           
           631   -1    0    0.000000    nu'_tau         nu'_taubar                                                      

        61     61    H_L++           H_L--               6    0    1    200.00000     0.88161     8.81606   0.00000E+00    1
           632    1    0    0.090264    e+              e+                                                              
           633    1    0    0.001805    e+              mu+                                                             
           634    1    0    0.001805    e+              tau+                                                            
           635    1    0    0.090264    mu+             mu+                                                             
           636    1    0    0.001805    mu+             tau+                                                            
           637    1    0    0.812250    tau+            tau+                                                            
           638    1    0    0.001806    W+              W+                                                              

        62     62    H_R++           H_R--               6    0    1    200.00000     0.88001     8.80013   0.00000E+00    1
           639    1    0    0.090428    e+              e+                                                              
           640    1    0    0.001809    e+              mu+                                                             
           641    1    0    0.001808    e+              tau+                                                            
           642    1    0    0.090428    mu+             mu+                                                             
           643    1    0    0.001808    mu+             tau+                                                            
           644    1    0    0.813720    tau+            tau+                                                            
           645    1    0    0.000000    W_R+            W_R+                                                            

        63     63    W_R+            W_R-                3    0    1    750.00000    19.32815   193.28147   0.00000E+00    1
           646    1   32    0.325914    dbar            u                                                               
           647    1   32    0.016735    dbar            c                                                               
           648    1   32    0.000009    dbar            t                                                               
           649    1   32    0.016735    sbar            u                                                               
           650    1   32    0.325320    sbar            c                                                               
           651    1   32    0.000554    sbar            t                                                               
           652    1   32    0.000010    bbar            u                                                               
           653    1   32    0.000603    bbar            c                                                               
           654    1   32    0.314119    bbar            t                                                               
           655    1    0    0.000000    e+              nu_Re                                                           
           656    1    0    0.000000    mu+             nu_Rmu                                                          
           657    1    0    0.000000    tau+            nu_Rtau                                                         

        64     64    nu_Re           nu_Rebar            0    0    1    750.00000     0.00000     0.00000   0.00000E+00    0

        65     65    nu_Rmu          nu_Rmubar           0    0    1    750.00000     0.00000     0.00000   0.00000E+00    0

        66     66    nu_Rtau         nu_Rtaubar          0    0    1    750.00000     0.00000     0.00000   0.00000E+00    0

        81     81    specflav                            0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        82     82    rndmflav        rndmflavbar         0    0    1      0.00000     0.00000     0.00000   0.00000E+00    0

        83     83    phasespa                            0    0    0      1.00000     0.00000     0.00000   0.00000E+00    1
           658    1   12    1.000000    rndmflav        rndmflavbar                                                     

        84     84    c-hadron        c-hadronbar         2    0    1      2.00000     0.00000     0.00000   1.00000E-01    1
           659    1   42    0.080000    e+              nu_e            s               specflav                        
           660    1   42    0.080000    mu+             nu_mu           s               specflav                        
           661    1   11    0.760000    u               dbar            s               specflav                        
           662    1   11    0.080000    u               sbar            s               specflav                        

        85     85    b-hadron        b-hadronbar        -1    0    1      5.00000     0.00000     0.00000   3.87000E-01    1
           663    1   42    0.105000    nu_ebar         e-              c               specflav                        
           664    1   42    0.105000    nu_mubar        mu-             c               specflav                        
           665    1   42    0.040000    nu_taubar       tau-            c               specflav                        
           666    1   42    0.500000    ubar            d               c               specflav                        
           667    1   42    0.080000    ubar            c               d               specflav                        
           668    1   42    0.140000    cbar            s               c               specflav                        
           669    1   42    0.010000    cbar            c               s               specflav                        
           670    1   42    0.015000    ubar            d               u               specflav                        
           671    1   42    0.005000    cbar            s               u               specflav                        

        91     91    cluster                             0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        92     92    string                              0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        93     93    indep.                              0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        94     94    CMshower                            0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        95     95    SPHEaxis                            0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        96     96    THRUaxis                            0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        97     97    CLUSjet                             0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        98     98    CELLjet                             0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

        99     99    table                               0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

       110    101    rho_diff0                           0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

       111    102    pi0                                 0    0    0      0.13498     0.00000     0.00000   3.00000E-05    0
           672    1    0    0.988000    gamma           gamma                                                           
           673    1    2    0.012000    gamma           e-              e+                                              

       113    103    rho0                                0    0    0      0.76850     0.15100     0.40000   0.00000E+00    1
           674    1    3    0.998739    pi+             pi-                                                             
           675    1    0    0.000790    pi0             gamma                                                           
           676    1    0    0.000380    eta             gamma                                                           
           677    1    0    0.000046    mu-             mu+                                                             
           678    1    0    0.000045    e-              e+                                                              

       115    104    a_20                                0    0    0      1.31800     0.10700     0.25000   0.00000E+00    1
           679    1    0    0.347250    rho+            pi-                                                             
           680    1    0    0.347250    rho-            pi+                                                             
           681    1    0    0.144000    eta             pi0                                                             
           682    1    0    0.104000    omega           pi+             pi-                                             
           683    1    0    0.024500    K+              K-                                                              
           684    1    0    0.012250    K_L0            K_L0                                                            
           685    1    0    0.012250    K_S0            K_S0                                                            
           686    1    0    0.002800    pi0             gamma                                                           
           687    1    0    0.005700    eta'            pi0                                                             

       130    105    K_L0                                0    0    0      0.49767     0.00000     0.00000   1.55000E+04    0
           688    1    0    0.211200    pi0             pi0             pi0                                             
           689    1    0    0.125600    pi+             pi-             pi0                                             
           690    1   42    0.193900    nu_ebar         e-              pi+                                             
           691    1   42    0.193900    nu_e            e+              pi-                                             
           692    1   42    0.135900    nu_mubar        mu-             pi+                                             
           693    1   42    0.135900    nu_mu           mu+             pi-                                             
           694    1    0    0.002000    pi+             pi-                                                             
           695    1    0    0.001000    pi0             pi0                                                             
           696    1    0    0.000600    gamma           gamma                                                           

       210    106    pi_diffr+       pi_diffr-           3    0    1      0.00000     0.00000     0.00000   0.00000E+00    0

       211    107    pi+             pi-                 3    0    1      0.13957     0.00000     0.00000   7.80450E+03    0
           697    1    0    0.999877    mu+             nu_mu                                                           
           698    1    0    0.000123    e+              nu_e                                                            

       213    108    rho+            rho-                3    0    1      0.76690     0.14900     0.40000   0.00000E+00    1
           699    1    3    0.999550    pi+             pi0                                                             
           700    1    0    0.000450    pi+             gamma                                                           

       215    109    a_2+            a_2-                3    0    1      1.31800     0.10700     0.25000   0.00000E+00    1
           701    1    0    0.347250    rho+            pi0                                                             
           702    1    0    0.347250    rho0            pi+                                                             
           703    1    0    0.144000    eta             pi+                                                             
           704    1    0    0.104000    omega           pi+             pi0                                             
           705    1    0    0.049000    K+              Kbar0                                                           
           706    1    0    0.002800    pi+             gamma                                                           
           707    1    0    0.005700    eta'            pi+                                                             

       220    110    omega_di                            0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

       221    111    eta                                 0    0    0      0.54745     0.00000     0.00000   0.00000E+00    1
           708    1    0    0.392300    gamma           gamma                                                           
           709    1    0    0.321000    pi0             pi0             pi0                                             
           710    1    0    0.231700    pi+             pi-             pi0                                             
           711    1    0    0.047800    gamma           pi+             pi-                                             
           712    1    2    0.004900    gamma           e-              e+                                              
           713    1    0    0.001300    pi+             pi-             e-              e+                              
           714    1    0    0.000300    gamma           mu-             mu+                                             
           715    1    0    0.000700    pi0             gamma           gamma                                           

       223    112    omega                               0    0    0      0.78194     0.00843     0.10000   0.00000E+00    1
           716    1    1    0.890000    pi+             pi-             pi0                                             
           717    1    0    0.086930    gamma           pi0                                                             
           718    1    3    0.022100    pi+             pi-                                                             
           719    1    0    0.000830    eta             gamma                                                           
           720    1    0    0.000070    pi0             pi0             gamma                                           
           721    1    0    0.000070    e-              e+                                                              

       225    113    f_2                                 0    0    0      1.27500     0.18500     0.17000   0.00000E+00    1
           722    1    0    0.564000    pi+             pi-                                                             
           723    1    0    0.282000    pi0             pi0                                                             
           724    1    0    0.072000    pi+             pi-             pi0             pi0                             
           725    1    0    0.028000    pi+             pi-             pi+             pi-                             
           726    1    0    0.023000    K+              K-                                                              
           727    1    0    0.011500    K_L0            K_L0                                                            
           728    1    0    0.011500    K_S0            K_S0                                                            
           729    1    0    0.005000    eta             eta                                                             
           730    1    0    0.003000    pi0             pi0             pi0             pi0                             

       310    114    K_S0                                0    0    0      0.49767     0.00000     0.00000   2.67620E+01    1
           731    1    0    0.686100    pi+             pi-                                                             
           732    1    0    0.313900    pi0             pi0                                                             

       311    115    K0              Kbar0               0    0    1      0.49767     0.00000     0.00000   0.00000E+00    1
           733    1    0    0.500000    K_L0                                                                            
           734    1    0    0.500000    K_S0                                                                            

       313    116    K*0             K*bar0              0    0    1      0.89610     0.05050     0.20000   0.00000E+00    1
           735    1    3    0.665000    K+              pi-                                                             
           736    1    3    0.333000    K0              pi0                                                             
           737    1    0    0.002000    K0              gamma                                                           

       315    117    K*_20           K*_2bar0            0    0    1      1.43200     0.10900     0.12000   0.00000E+00    1
           738    1    0    0.333000    K+              pi-                                                             
           739    1    0    0.166000    K0              pi0                                                             
           740    1    0    0.168000    K*+             pi-                                                             
           741    1    0    0.084000    K*0             pi0                                                             
           742    1    0    0.087000    K*+             pi-             pi0                                             
           743    1    0    0.043000    K*0             pi+             pi-                                             
           744    1    0    0.059000    K+              rho-                                                            
           745    1    0    0.029000    K0              rho0                                                            
           746    1    0    0.029000    K0              omega                                                           
           747    1    0    0.002000    K0              eta                                                             

       321    118    K+              K-                  3    0    1      0.49360     0.00000     0.00000   3.70900E+03    0
           748    1    0    0.635200    mu+             nu_mu                                                           
           749    1    0    0.211600    pi+             pi0                                                             
           750    1    0    0.055900    pi+             pi+             pi-                                             
           751    1    0    0.017300    pi+             pi0             pi0                                             
           752    1   42    0.048200    nu_e            e+              pi0                                             
           753    1   42    0.031800    nu_mu           mu+             pi0                                             

       323    119    K*+             K*-                 3    0    1      0.89160     0.04980     0.20000   0.00000E+00    1
           754    1    3    0.666000    K0              pi+                                                             
           755    1    3    0.333000    K+              pi0                                                             
           756    1    0    0.001000    K+              gamma                                                           

       325    120    K*_2+           K*_2-               3    0    1      1.42500     0.09800     0.12000   0.00000E+00    1
           757    1    0    0.332000    K0              pi+                                                             
           758    1    0    0.166000    K+              pi0                                                             
           759    1    0    0.168000    K*0             pi+                                                             
           760    1    0    0.084000    K*+             pi0                                                             
           761    1    0    0.086000    K*0             pi+             pi0                                             
           762    1    0    0.043000    K*+             pi+             pi-                                             
           763    1    0    0.059000    K0              rho+                                                            
           764    1    0    0.029000    K+              rho0                                                            
           765    1    0    0.029000    K+              omega                                                           
           766    1    0    0.002000    K+              eta                                                             
           767    1    0    0.002000    K+              gamma                                                           

       330    121    phi_diff                            0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

       331    122    eta'                                0    0    0      0.95777     0.00020     0.00200   0.00000E+00    1
           768    1    0    0.437000    pi+             pi-             eta                                             
           769    1    0    0.208000    pi0             pi0             eta                                             
           770    1    0    0.302000    gamma           rho0                                                            
           771    1    0    0.030200    gamma           omega                                                           
           772    1    0    0.021200    gamma           gamma                                                           
           773    1    0    0.001600    pi0             pi0             pi0                                             

       333    123    phi                                 0    0    0      1.01940     0.00443     0.01500   0.00000E+00    1
           774    1    3    0.489470    K+              K-                                                              
           775    1    3    0.340000    K_L0            K_S0                                                            
           776    1    0    0.043000    rho-            pi+                                                             
           777    1    0    0.043000    rho0            pi0                                                             
           778    1    0    0.043000    rho+            pi-                                                             
           779    1    1    0.027000    pi+             pi-             pi0                                             
           780    1    0    0.012600    gamma           eta                                                             
           781    1    0    0.001300    pi0             gamma                                                           
           782    1    0    0.000300    e-              e+                                                              
           783    1    0    0.000250    mu-             mu+                                                             
           784    1    0    0.000080    pi+             pi-                                                             

       335    124    f'_2                                0    0    0      1.52500     0.07600     0.20000   0.00000E+00    1
           785    1    0    0.444000    K+              K-                                                              
           786    1    0    0.222000    K_L0            K_L0                                                            
           787    1    0    0.222000    K_S0            K_S0                                                            
           788    1    0    0.104000    eta             eta                                                             
           789    1    0    0.004000    pi+             pi-                                                             
           790    1    0    0.004000    pi0             pi0                                                             

       411    125    D+              D-                  3    0    1      1.86930     0.00000     0.00000   3.17000E-01    1
           791    0   42    0.070000    e+              nu_e            Kbar0                                           
           792    0   42    0.065000    e+              nu_e            K*bar0                                          
           793    0   42    0.005000    e+              nu_e            Kbar0           pi0                             
           794    0   42    0.005000    e+              nu_e            K-              pi+                             
           795    0   42    0.011000    e+              nu_e            K*bar0          pi0                             
           796    0   42    0.011000    e+              nu_e            K*-             pi+                             
           797    0   42    0.001000    e+              nu_e            pi0                                             
           798    0   42    0.001000    e+              nu_e            eta                                             
           799    0   42    0.001000    e+              nu_e            eta'                                            
           800    0   42    0.001000    e+              nu_e            rho0                                            
           801    0   42    0.001000    e+              nu_e            omega                                           
           802    1   42    0.070000    mu+             nu_mu           Kbar0                                           
           803    1   42    0.065000    mu+             nu_mu           K*bar0                                          
           804    1   42    0.005000    mu+             nu_mu           Kbar0           pi0                             
           805    1   42    0.005000    mu+             nu_mu           K-              pi+                             
           806    1   42    0.011000    mu+             nu_mu           K*bar0          pi0                             
           807    1   42    0.011000    mu+             nu_mu           K*-             pi+                             
           808    1   42    0.001000    mu+             nu_mu           pi0                                             
           809    1   42    0.001000    mu+             nu_mu           eta                                             
           810    1   42    0.001000    mu+             nu_mu           eta'                                            
           811    1   42    0.001000    mu+             nu_mu           rho0                                            
           812    1   42    0.001000    mu+             nu_mu           omega                                           
           813    0    0    0.026000    Kbar0           pi+                                                             
           814    0    0    0.019000    K*bar0          pi+                                                             
           815    0    0    0.066000    Kbar0           rho+                                                            
           816    0    0    0.041000    K*bar0          rho+                                                            
           817    0    0    0.045000    K*_1bar0        pi+                                                             
           818    0    0    0.076000    Kbar0           a_1+                                                            
           819    0    0    0.007300    Kbar0           K+                                                              
           820    0    0    0.004700    K*bar0          K+                                                              
           821    0    0    0.004700    Kbar0           K*+                                                             
           822    0    0    0.026000    K*bar0          K*+                                                             
           823    0    0    0.001000    pi0             pi+                                                             
           824    0    0    0.000600    pi0             rho+                                                            
           825    0    0    0.006600    eta             pi+                                                             
           826    0    0    0.005000    eta             rho+                                                            
           827    0    0    0.003000    eta'            pi+                                                             
           828    0    0    0.003000    eta'            rho+                                                            
           829    0    0    0.000600    rho0            pi+                                                             
           830    0    0    0.000600    rho0            rho+                                                            
           831    0    0    0.001000    omega           pi+                                                             
           832    0    0    0.001000    omega           rho+                                                            
           833    0    0    0.006000    phi             pi+                                                             
           834    0    0    0.005000    phi             rho+                                                            
           835    0    0    0.012000    Kbar0           pi+             pi0                                             
           836    0    0    0.005700    K*bar0          pi+             rho0                                            
           837    0    0    0.067000    K-              pi+             pi+                                             
           838    0    0    0.008000    K-              rho+            pi+                                             
           839    0    0    0.002200    pi+             pi+             pi-                                             
           840    0    0    0.027000    Kbar0           K+              Kbar0                                           
           841    0    0    0.004000    K-              K+              pi+                                             
           842    0    0    0.019000    phi             pi+             pi0                                             
           843    0    0    0.012000    Kbar0           pi+             pi+             pi-                             
           844    0    0    0.002000    K*bar0          pi+             pi+             pi-                             
           845    0    0    0.009000    K-              pi+             pi+             pi0                             
           846    0    0    0.021800    pi+             pi+             pi-             pi0                             
           847    0    0    0.001000    K-              pi+             pi+             pi+             pi-             
           848    0    0    0.022000    K-              pi+             pi+             pi0             pi0             
           849    0    0    0.087000    Kbar0           pi+             pi+             pi-             pi0             
           850    0    0    0.001000    Kbar0           rho0            pi+             pi+             pi-             
           851    0    0    0.001900    K-              rho0            pi+             pi+             pi0             
           852    0    0    0.001500    pi+             pi+             pi+             pi-             pi-             
           853    0    0    0.002800    rho0            pi+             pi+             pi-             pi0             

       413    126    D*+             D*-                 3    0    1      2.01000     0.00000     0.00000   0.00000E+00    1
           854    1    3    0.683000    D0              pi+                                                             
           855    1    3    0.306000    D+              pi0                                                             
           856    1    0    0.011000    D+              gamma                                                           

       415    127    D*_2+           D*_2-               3    0    1      2.46000     0.02300     0.12000   0.00000E+00    1
           857    1    0    0.300000    D0              pi+                                                             
           858    1    0    0.150000    D+              pi0                                                             
           859    1    0    0.160000    D*0             pi+                                                             
           860    1    0    0.080000    D*+             pi0                                                             
           861    1    0    0.130000    D*0             pi+             pi0                                             
           862    1    0    0.060000    D*+             pi+             pi-                                             
           863    1    0    0.080000    D0              pi+             pi0                                             
           864    1    0    0.040000    D+              pi+             pi-                                             

       421    128    D0              Dbar0               0    0    1      1.86450     0.00000     0.00000   1.24400E-01    1
           865    0   42    0.034000    e+              nu_e            K-                                              
           866    0   42    0.027000    e+              nu_e            K*-                                             
           867    0   42    0.002000    e+              nu_e            Kbar0           pi-                             
           868    0   42    0.002000    e+              nu_e            K-              pi0                             
           869    0   42    0.004000    e+              nu_e            K*bar0          pi-                             
           870    0   42    0.004000    e+              nu_e            K*-             pi0                             
           871    0   42    0.002000    e+              nu_e            pi-                                             
           872    0   42    0.002000    e+              nu_e            rho-                                            
           873    1   42    0.034000    mu+             nu_mu           K-                                              
           874    1   42    0.027000    mu+             nu_mu           K*-                                             
           875    1   42    0.002000    mu+             nu_mu           Kbar0           pi-                             
           876    1   42    0.002000    mu+             nu_mu           K-              pi0                             
           877    1   42    0.004000    mu+             nu_mu           K*bar0          pi-                             
           878    1   42    0.004000    mu+             nu_mu           K*-             pi0                             
           879    1   42    0.002000    mu+             nu_mu           pi-                                             
           880    1   42    0.002000    mu+             nu_mu           rho-                                            
           881    0    0    0.036500    K-              pi+                                                             
           882    0    0    0.045000    K*-             pi+                                                             
           883    0    0    0.073000    K-              rho+                                                            
           884    0    0    0.062000    K*-             rho+                                                            
           885    0    0    0.021000    Kbar0           pi0                                                             
           886    0    0    0.021000    K*bar0          pi0                                                             
           887    0    0    0.021000    K*bar0          eta                                                             
           888    0    0    0.006100    Kbar0           rho0                                                            
           889    0    0    0.015000    K*bar0          rho0                                                            
           890    0    0    0.025000    Kbar0           omega                                                           
           891    0    0    0.008800    Kbar0           phi                                                             
           892    0    0    0.074000    K-              a_1+                                                            
           893    0    0    0.010900    K_1-            pi+                                                             
           894    0    0    0.004100    K-              K+                                                              
           895    0    0    0.002000    K*-             K+                                                              
           896    0    0    0.003500    K-              K*+                                                             
           897    0    0    0.001100    Kbar0           K0                                                              
           898    0    0    0.001000    K*bar0          K0                                                              
           899    0    0    0.002700    K*bar0          K*0                                                             
           900    0    0    0.001600    pi+             pi-                                                             
           901    0    0    0.001600    pi0             pi0                                                             
           902    0    0    0.001800    phi             rho0                                                            
           903    0    0    0.011000    K-              pi+             pi0                                             
           904    0    0    0.006300    K-              pi+             rho0                                            
           905    0    0    0.005200    K-              K+              Kbar0                                           
           906    0    0    0.018000    Kbar0           pi+             pi-                                             
           907    0    0    0.016000    K*bar0          pi+             pi-                                             
           908    0    0    0.003400    K-              K0              pi+                                             
           909    0    0    0.003600    K*bar0          K+              pi-                                             
           910    0    0    0.000900    K_S0            K_S0            K_S0                                            
           911    0    0    0.000600    phi             pi+             pi-                                             
           912    0    0    0.015000    pi+             pi-             pi0                                             
           913    0    0    0.092300    K-              pi+             pi0             pi0                             
           914    0    0    0.018000    K-              pi+             pi+             pi-                             
           915    0    0    0.022000    Kbar0           pi+             pi-             pi0                             
           916    0    0    0.007700    K*bar0          pi+             pi-             pi0                             
           917    0    0    0.009000    Kbar0           K+              K-              pi0                             
           918    0    0    0.007500    pi+             pi+             pi-             pi-                             
           919    0    0    0.024000    K-              pi+             pi+             pi-             pi0             
           920    0    0    0.008500    Kbar0           pi+             pi+             pi-             pi-             
           921    0    0    0.067000    Kbar0           pi+             pi-             pi0             pi0             
           922    0    0    0.051100    Kbar0           rho0            pi0             pi0             pi0             
           923    0    0    0.017000    pi+             pi+             pi-             pi-             pi0             
           924    0    0    0.000400    rho0            pi+             pi+             pi-             pi-             
           925    0    0    0.002800    K+              K-              pi+             pi-             pi0             

       423    129    D*0             D*bar0              0    0    1      2.00670     0.00000     0.00000   0.00000E+00    1
           926    1    3    0.619000    D0              pi0                                                             
           927    1    0    0.381000    D0              gamma                                                           

       425    130    D*_20           D*_2bar0            0    0    1      2.46000     0.02300     0.12000   0.00000E+00    1
           928    1    0    0.300000    D+              pi-                                                             
           929    1    0    0.150000    D0              pi0                                                             
           930    1    0    0.160000    D*+             pi-                                                             
           931    1    0    0.080000    D*0             pi0                                                             
           932    1    0    0.130000    D*+             pi-             pi0                                             
           933    1    0    0.060000    D*0             pi+             pi-                                             
           934    1    0    0.080000    D+              pi-             pi0                                             
           935    1    0    0.040000    D0              pi+             pi-                                             

       431    131    D_s+            D_s-                3    0    1      1.96850     0.00000     0.00000   1.40000E-01    1
           936    0    0    0.010000    tau+            nu_tau                                                          
           937    0   42    0.020000    e+              nu_e            eta                                             
           938    0   42    0.020000    e+              nu_e            eta'                                            
           939    0   42    0.030000    e+              nu_e            phi                                             
           940    0   42    0.005000    e+              nu_e            K+              K-                              
           941    0   42    0.005000    e+              nu_e            K0              Kbar0                           
           942    1   42    0.020000    mu+             nu_mu           eta                                             
           943    1   42    0.020000    mu+             nu_mu           eta'                                            
           944    1   42    0.030000    mu+             nu_mu           phi                                             
           945    1   42    0.005000    mu+             nu_mu           K+              K-                              
           946    1   42    0.005000    mu+             nu_mu           K0              Kbar0                           
           947    0    0    0.015000    eta             pi+                                                             
           948    0    0    0.037000    eta'            pi+                                                             
           949    0    0    0.028000    phi             pi+                                                             
           950    0    0    0.079000    eta             rho+                                                            
           951    0    0    0.095000    eta'            rho+                                                            
           952    0    0    0.052000    phi             rho+                                                            
           953    0    0    0.007800    f_0             pi+                                                             
           954    0    0    0.001000    pi+             pi0                                                             
           955    0    0    0.001000    rho+            pi0                                                             
           956    0    0    0.001000    pi+             rho0                                                            
           957    0    0    0.001000    rho+            rho0                                                            
           958    0    0    0.028000    K+              Kbar0                                                           
           959    0    0    0.033000    K*+             Kbar0                                                           
           960    0    0    0.026000    K+              K*bar0                                                          
           961    0    0    0.050000    K*+             K*bar0                                                          
           962    0    0    0.010000    p+              nbar0                                                           
           963    0    0    0.005000    eta             K+                                                              
           964    0    0    0.005000    eta'            K+                                                              
           965    0    0    0.005000    phi             K+                                                              
           966    0    0    0.005000    eta             K*+                                                             
           967    0   13    0.250000    u               dbar            s               sbar                            
           968    0   13    0.095200    u               dbar                                                            

       433    132    D*_s+           D*_s-               3    0    1      2.11240     0.00000     0.00000   0.00000E+00    1
           969    1    0    0.940000    D_s+            gamma                                                           
           970    1    0    0.060000    D_s+            pi0                                                             

       435    133    D*_2s+          D*_2s-              3    0    1      2.57350     0.01500     0.05000   0.00000E+00    1
           971    1    0    0.400000    D0              K+                                                              
           972    1    0    0.400000    D+              K0                                                              
           973    1    0    0.100000    D*0             K+                                                              
           974    1    0    0.100000    D*+             K0                                                              

       440    134    J/psi_di                            0    0    0      0.00000     0.00000     0.00000   0.00000E+00    0

       441    135    eta_c                               0    0    0      2.97980     0.00130     0.00500   0.00000E+00    1
           975    1   12    1.000000    rndmflav        rndmflavbar                                                     

       443    136    J/psi                               0    0    0      3.09688     0.00000     0.00000   0.00000E+00    1
           976    1    0    0.060200    e-              e+                                                              
           977    1    0    0.060100    mu-             mu+                                                             
           978    1   12    0.879700    rndmflav        rndmflavbar                                                     

       445    137    chi_2c                              0    0    0      3.55620     0.00200     0.01000   0.00000E+00    1
           979    1    0    0.135000    J/psi           gamma                                                           
           980    1   12    0.865000    rndmflav        rndmflavbar                                                     

       511    138    B0              Bbar0               0    0    1      5.27920     0.00000     0.00000   4.68000E-01    1
           981    0   42    0.020000    nu_e            e+              D-                                              
           982    0   42    0.055000    nu_e            e+              D*-                                             
           983    0   42    0.005000    nu_e            e+              D_1-                                            
           984    0   42    0.005000    nu_e            e+              D*_0-                                           
           985    0   42    0.008000    nu_e            e+              D*_1-                                           
           986    0   42    0.012000    nu_e            e+              D*_2-                                           
           987    1   42    0.020000    nu_mu           mu+             D-                                              
           988    1   42    0.055000    nu_mu           mu+             D*-                                             
           989    1   42    0.005000    nu_mu           mu+             D_1-                                            
           990    1   42    0.005000    nu_mu           mu+             D*_0-                                           
           991    1   42    0.008000    nu_mu           mu+             D*_1-                                           
           992    1   42    0.012000    nu_mu           mu+             D*_2-                                           
           993    0   42    0.010000    nu_tau          tau+            D-                                              
           994    0   42    0.030000    nu_tau          tau+            D*-                                             
           995    0    0    0.003500    D-              pi+                                                             
           996    0    0    0.011000    D-              rho+                                                            
           997    0    0    0.005500    D-              a_1+                                                            
           998    0    0    0.004200    D*-             pi+                                                             
           999    0    0    0.009000    D*-             rho+                                                            
          1000    0    0    0.018000    D*-             a_1+                                                            
          1001    0    0    0.015000    D-              D_s+                                                            
          1002    0    0    0.018500    D-              D*_s+                                                           
          1003    0    0    0.013500    D*-             D_s+                                                            
          1004    0    0    0.025000    D*-             D*_s+                                                           
          1005    0    0    0.000400    eta_c           K0                                                              
          1006    0    0    0.000700    eta_c           K*0                                                             
          1007    0    0    0.000800    J/psi           K0                                                              
          1008    0    0    0.001400    J/psi           K*0                                                             
          1009    0    0    0.001900    chi_1c          K0                                                              
          1010    0    0    0.002500    chi_1c          K*0                                                             
          1011    0   48    0.429100    u               dbar            cbar            d                               
          1012    0   13    0.080000    u               cbar            dbar            d                               
          1013    0   13    0.070000    c               sbar            cbar            d                               
          1014    0   13    0.020000    c               cbar            sbar            d                               
          1015    0   42    0.015000    u               dbar            ubar            d                               
          1016    0   42    0.005000    c               sbar            ubar            d                               

       513    139    B*0             B*bar0              0    0    1      5.32480     0.00000     0.00000   0.00000E+00    1
          1017    1    0    1.000000    B0              gamma                                                           

       515    140    B*_20           B*_2bar0            0    0    1      5.83000     0.02000     0.05000   0.00000E+00    1
          1018    1    0    0.300000    B+              pi-                                                             
          1019    1    0    0.150000    B0              pi0                                                             
          1020    1    0    0.160000    B*+             pi-                                                             
          1021    1    0    0.080000    B*0             pi0                                                             
          1022    1    0    0.130000    B*+             pi-             pi0                                             
          1023    1    0    0.060000    B*0             pi+             pi-                                             
          1024    1    0    0.080000    B+              pi-             pi0                                             
          1025    1    0    0.040000    B0              pi+             pi-                                             

       521    141    B+              B-                  3    0    1      5.27890     0.00000     0.00000   4.62000E-01    1
          1026    0   42    0.020000    nu_e            e+              Dbar0                                           
          1027    0   42    0.055000    nu_e            e+              D*bar0                                          
          1028    0   42    0.005000    nu_e            e+              D_1bar0                                         
          1029    0   42    0.005000    nu_e            e+              D*_0bar0                                        
          1030    0   42    0.008000    nu_e            e+              D*_1bar0                                        
          1031    0   42    0.012000    nu_e            e+              D*_2bar0                                        
          1032    1   42    0.020000    nu_mu           mu+             Dbar0                                           
          1033    1   42    0.055000    nu_mu           mu+             D*bar0                                          
          1034    1   42    0.005000    nu_mu           mu+             D_1bar0                                         
          1035    1   42    0.005000    nu_mu           mu+             D*_0bar0                                        
          1036    1   42    0.008000    nu_mu           mu+             D*_1bar0                                        
          1037    1   42    0.012000    nu_mu           mu+             D*_2bar0                                        
          1038    0   42    0.010000    nu_tau          tau+            Dbar0                                           
          1039    0   42    0.030000    nu_tau          tau+            D*bar0                                          
          1040    0    0    0.003500    Dbar0           pi+                                                             
          1041    0    0    0.011000    Dbar0           rho+                                                            
          1042    0    0    0.005500    Dbar0           a_1+                                                            
          1043    0    0    0.004200    D*bar0          pi+                                                             
          1044    0    0    0.009000    D*bar0          rho+                                                            
          1045    0    0    0.018000    D*bar0          a_1+                                                            
          1046    0    0    0.015000    Dbar0           D_s+                                                            
          1047    0    0    0.018500    Dbar0           D*_s+                                                           
          1048    0    0    0.013500    D*bar0          D_s+                                                            
          1049    0    0    0.025000    D*bar0          D*_s+                                                           
          1050    0    0    0.000400    eta_c           K+                                                              
          1051    0    0    0.000700    eta_c           K*+                                                             
          1052    0    0    0.000800    J/psi           K+                                                              
          1053    0    0    0.001400    J/psi           K*+                                                             
          1054    0    0    0.001900    chi_1c          K+                                                              
          1055    0    0    0.002500    chi_1c          K*+                                                             
          1056    0   48    0.429100    u               dbar            cbar            u                               
          1057    0   13    0.080000    u               cbar            dbar            u                               
          1058    0   13    0.070000    c               sbar            cbar            u                               
          1059    0   13    0.020000    c               cbar            sbar            u                               
          1060    0   42    0.015000    u               dbar            ubar            u                               
          1061    0   42    0.005000    c               sbar            ubar            u                               

       523    142    B*+             B*-                 3    0    1      5.32480     0.00000     0.00000   0.00000E+00    1
          1062    1    0    1.000000    B+              gamma                                                           

       525    143    B*_2+           B*_2-               3    0    1      5.83000     0.02000     0.05000   0.00000E+00    1
          1063    1    0    0.300000    B0              pi+                                                             
          1064    1    0    0.150000    B+              pi0                                                             
          1065    1    0    0.160000    B*0             pi+                                                             
          1066    1    0    0.080000    B*+             pi0                                                             
          1067    1    0    0.130000    B*0             pi+             pi0                                             
          1068    1    0    0.060000    B*+             pi+             pi-                                             
          1069    1    0    0.080000    B0              pi+             pi0                                             
          1070    1    0    0.040000    B+              pi+             pi-                                             

       531    144    B_s0            B_sbar0             0    0    1      5.36930     0.00000     0.00000   4.83000E-01    1
          1071    0   42    0.020000    nu_e            e+              D_s-                                            
          1072    0   42    0.055000    nu_e            e+              D*_s-                                           
          1073    0   42    0.005000    nu_e            e+              D_1s-                                           
          1074    0   42    0.005000    nu_e            e+              D*_0s-                                          
          1075    0   42    0.008000    nu_e            e+              D*_1s-                                          
          1076    0   42    0.012000    nu_e            e+              D*_2s-                                          
          1077    1   42    0.020000    nu_mu           mu+             D_s-                                            
          1078    1   42    0.055000    nu_mu           mu+             D*_s-                                           
          1079    1   42    0.005000    nu_mu           mu+             D_1s-                                           
          1080    1   42    0.005000    nu_mu           mu+             D*_0s-                                          
          1081    1   42    0.008000    nu_mu           mu+             D*_1s-                                          
          1082    1   42    0.012000    nu_mu           mu+             D*_2s-                                          
          1083    0   42    0.010000    nu_tau          tau+            D_s-                                            
          1084    0   42    0.030000    nu_tau          tau+            D*_s-                                           
          1085    0    0    0.003500    D_s-            pi+                                                             
          1086    0    0    0.011000    D_s-            rho+                                                            
          1087    0    0    0.005500    D_s-            a_1+                                                            
          1088    0    0    0.004200    D*_s-           pi+                                                             
          1089    0    0    0.009000    D*_s-           rho+                                                            
          1090    0    0    0.018000    D*_s-           a_1+                                                            
          1091    0    0    0.015000    D_s-            D_s+                                                            
          1092    0    0    0.018500    D_s-            D*_s+                                                           
          1093    0    0    0.013500    D*_s-           D_s+                                                            
          1094    0    0    0.025000    D*_s-           D*_s+                                                           
          1095    0    0    0.000200    eta_c           eta                                                             
          1096    0    0    0.000200    eta_c           eta'                                                            
          1097    0    0    0.000700    eta_c           phi                                                             
          1098    0    0    0.000400    J/psi           eta                                                             
          1099    0    0    0.000400    J/psi           eta'                                                            
          1100    0    0    0.001400    J/psi           phi                                                             
          1101    0    0    0.001000    chi_1c          eta                                                             
          1102    0    0    0.000900    chi_1c          eta'                                                            
          1103    0    0    0.002500    chi_1c          phi                                                             
          1104    0   48    0.429100    u               dbar            cbar            s                               
          1105    0   13    0.080000    u               cbar            dbar            s                               
          1106    0   13    0.070000    c               sbar            cbar            s                               
          1107    0   13    0.020000    c               cbar            sbar            s                               
          1108    0   42    0.015000    u               dbar            ubar            s                               
          1109    0   42    0.005000    c               sbar            ubar            s                               

       533    145    B*_s0           B*_sbar0            0    0    1      5.41630     0.00000     0.00000   0.00000E+00    1
          1110    1    0    1.000000    B_s0            gamma                                                           

       535    146    B*_2s0          B*_2sbar0           0    0    1      6.07000     0.02000     0.05000   0.00000E+00    1
          1111    1    0    0.300000    B+              K-                                                              
          1112    1    0    0.300000    B0              Kbar0                                                           
          1113    1    0    0.200000    B*+             K-                                                              
          1114    1    0    0.200000    B*0             Kbar0                                                           

       541    147    B_c+            B_c-                3    0    1      6.59400     0.00000     0.00000   1.50000E-01    1
          1115    1    0    0.047000    nu_tau          tau+                                                            
          1116    1   11    0.122000    c               sbar                                                            
          1117    1   11    0.006000    c               dbar                                                            
          1118    1   42    0.012000    nu_e            e+              eta_c                                           
          1119    1   42    0.035000    nu_e            e+              J/psi                                           
          1120    1   42    0.012000    nu_mu           mu+             eta_c                                           
          1121    1   42    0.035000    nu_mu           mu+             J/psi                                           
          1122    1   42    0.003000    nu_tau          tau+            eta_c                                           
          1123    1   42    0.007000    nu_tau          tau+            J/psi                                           
          1124    1   42    0.150000    u               dbar            cbar            c                               
          1125    1   42    0.037000    u               cbar            dbar            c                               
          1126    1   42    0.008000    u               sbar            cbar            c                               
          1127    1   42    0.002000    u               cbar            sbar            c                               
          1128    1   42    0.050000    c               sbar            cbar            c                               
          1129    1   42    0.015000    c               cbar            sbar            c                               
          1130    1   42    0.003000    c               dbar            cbar            c                               
          1131    1   42    0.001000    c               cbar            dbar            c                               
          1132    1   42    0.014000    e+              nu_e            B_s0                                            
          1133    1   42    0.042000    e+              nu_e            B*_s0                                           
          1134    1   42    0.014000    mu+             nu_mu           B_s0                                            
          1135    1   42    0.042000    mu+             nu_mu           B*_s0                                           
          1136    1   42    0.240000    dbar            u               s               bbar                            
          1137    1   42    0.065000    dbar            s               u               bbar                            
          1138    1   42    0.012000    sbar            u               s               bbar                            
          1139    1   42    0.003000    sbar            s               u               bbar                            
          1140    1   42    0.001000    e+              nu_e            B0                                              
          1141    1   42    0.002000    e+              nu_e            B*0                                             
          1142    1   42    0.001000    mu+             nu_mu           B0                                              
          1143    1   42    0.002000    mu+             nu_mu           B*0                                             
          1144    1   42    0.014000    dbar            u               d               bbar                            
          1145    1   42    0.003000    dbar            d               u               bbar                            

       543    148    B*_c+           B*_c-               3    0    1      6.60200     0.00000     0.00000   0.00000E+00    1
          1146    1    0    1.000000    B_c+            gamma                                                           

       545    149    B*_2c+          B*_2c-              3    0    1      7.35000     0.02000     0.05000   0.00000E+00    1
          1147    1    0    0.300000    B0              D+                                                              
          1148    1    0    0.300000    B+              D0                                                              
          1149    1    0    0.200000    B*0             D+                                                              
          1150    1    0    0.200000    B*+             D0                                                              

       551    150    eta_b                               0    0    0      9.40000     0.00000     0.00000   0.00000E+00    1
          1151    1   32    1.000000    g               g                                                               

       553    151    Upsilon                             0    0    0      9.46030     0.00000     0.00000   0.00000E+00    1
          1152    1    0    0.025200    e-              e+                                                              
          1153    1    0    0.024800    mu-             mu+                                                             
          1154    1    0    0.026700    tau-            tau+                                                            
          1155    1   32    0.015000    d               dbar                                                            
          1156    1   32    0.045000    u               ubar                                                            
          1157    1   32    0.015000    s               sbar                                                            
          1158    1   32    0.045000    c               cbar                                                            
          1159    1    4    0.774300    g               g               g                                               
          1160    1    4    0.029000    gamma           g               g                                               

       555    152    chi_2b                              0    0    0      9.91320     0.00000     0.00000   0.00000E+00    1
          1161    1    0    0.220000    Upsilon         gamma                                                           
          1162    1   32    0.780000    g               g                                                               

      1103    153    dd_1            dd_1bar            -2   -1    1      0.77133     0.00000     0.00000   0.00000E+00    0

      1114    154    Delta-          Deltabar+          -3    0    1      1.23400     0.12000     0.14000   0.00000E+00    1
          1163    1    0    1.000000    n0              pi-                                                             

      2101    155    ud_0            ud_0bar             1   -1    1      0.57933     0.00000     0.00000   0.00000E+00    0

      2103    156    ud_1            ud_1bar             1   -1    1      0.77133     0.00000     0.00000   0.00000E+00    0

      2110    157    n_diffr0        n_diffrbar0         0    0    1      0.00000     0.00000     0.00000   0.00000E+00    0

      2112    158    n0              nbar0               0    0    1      0.93957     0.00000     0.00000   0.00000E+00    0

      2114    159    Delta0          Deltabar0           0    0    1      1.23300     0.12000     0.14000   0.00000E+00    1
          1164    1    0    0.331000    p+              pi-                                                             
          1165    1    0    0.663000    n0              pi0                                                             
          1166    1    0    0.006000    n0              gamma                                                           

      2203    160    uu_1            uu_1bar             4   -1    1      0.77133     0.00000     0.00000   0.00000E+00    0

      2210    161    p_diffr+        p_diffrbar-         3    0    1      0.00000     0.00000     0.00000   0.00000E+00    0

      2212    162    p+              pbar-               3    0    1      0.93827     0.00000     0.00000   0.00000E+00    0

      2214    163    Delta+          Deltabar-           3    0    1      1.23200     0.12000     0.14000   0.00000E+00    1
          1167    1    0    0.663000    p+              pi0                                                             
          1168    1    0    0.331000    n0              pi+                                                             
          1169    1    0    0.006000    p+              gamma                                                           

      2224    164    Delta++         Deltabar--          6    0    1      1.23100     0.12000     0.14000   0.00000E+00    1
          1170    1    0    1.000000    p+              pi+                                                             

      3101    165    sd_0            sd_0bar            -2   -1    1      0.80473     0.00000     0.00000   0.00000E+00    0

      3103    166    sd_1            sd_1bar            -2   -1    1      0.92953     0.00000     0.00000   0.00000E+00    0

      3112    167    Sigma-          Sigmabar+          -3    0    1      1.19744     0.00000     0.00000   4.43400E+01    1
          1171    1    0    0.999000    n0              pi-                                                             
          1172    1    0    0.001000    nu_ebar         e-              n0                                              

      3114    168    Sigma*-         Sigma*bar+         -3    0    1      1.38720     0.03940     0.04000   0.00000E+00    1
          1173    1    0    0.880000    Lambda0         pi-                                                             
          1174    1    0    0.060000    Sigma0          pi-                                                             
          1175    1    0    0.060000    Sigma-          pi0                                                             

      3122    169    Lambda0         Lambdabar0          0    0    1      1.11568     0.00000     0.00000   7.88800E+01    1
          1176    1    0    0.639000    p+              pi-                                                             
          1177    1    0    0.358000    n0              pi0                                                             
          1178    1    0    0.002000    n0              gamma                                                           
          1179    1    0    0.001000    nu_ebar         e-              p+                                              

      3201    170    su_0            su_0bar             1   -1    1      0.80473     0.00000     0.00000   0.00000E+00    0

      3203    171    su_1            su_1bar             1   -1    1      0.92953     0.00000     0.00000   0.00000E+00    0

      3212    172    Sigma0          Sigmabar0           0    0    1      1.19255     0.00000     0.00000   0.00000E+00    1
          1180    1    0    1.000000    Lambda0         gamma                                                           

      3214    173    Sigma*0         Sigma*bar0          0    0    1      1.38370     0.03600     0.03500   0.00000E+00    1
          1181    1    0    0.880000    Lambda0         pi0                                                             
          1182    1    0    0.060000    Sigma+          pi-                                                             
          1183    1    0    0.060000    Sigma-          pi+                                                             

      3222    174    Sigma+          Sigmabar-           3    0    1      1.18937     0.00000     0.00000   2.39600E+01    1
          1184    1    0    0.516000    p+              pi0                                                             
          1185    1    0    0.483000    n0              pi+                                                             
          1186    1    0    0.001000    p+              gamma                                                           

      3224    175    Sigma*+         Sigma*bar-          3    0    1      1.38280     0.03580     0.03500   0.00000E+00    1
          1187    1    0    0.880000    Lambda0         pi+                                                             
          1188    1    0    0.060000    Sigma+          pi0                                                             
          1189    1    0    0.060000    Sigma0          pi+                                                             

      3303    176    ss_1            ss_1bar            -2   -1    1      1.09361     0.00000     0.00000   0.00000E+00    0

      3312    177    Xi-             Xibar+             -3    0    1      1.32130     0.00000     0.00000   4.91000E+01    1
          1190    1    0    0.998800    Lambda0         pi-                                                             
          1191    1    0    0.000100    Sigma-          gamma                                                           
          1192    1    0    0.000600    nu_ebar         e-              Lambda0                                         
          1193    1    0    0.000400    nu_mubar        mu-             Lambda0                                         
          1194    1    0    0.000100    nu_ebar         e-              Sigma0                                          

      3314    178    Xi*-            Xi*bar+            -3    0    1      1.53500     0.00990     0.05000   0.00000E+00    1
          1195    1    0    0.667000    Xi0             pi-                                                             
          1196    1    0    0.333000    Xi-             pi0                                                             

      3322    179    Xi0             Xibar0              0    0    1      1.31490     0.00000     0.00000   8.71000E+01    1
          1197    1    0    0.995400    Lambda0         pi0                                                             
          1198    1    0    0.001100    Lambda0         gamma                                                           
          1199    1    0    0.003500    Sigma0          gamma                                                           

      3324    180    Xi*0            Xi*bar0             0    0    1      1.53180     0.00910     0.05000   0.00000E+00    1
          1200    1    0    0.333000    Xi0             pi0                                                             
          1201    1    0    0.667000    Xi-             pi+                                                             

      3334    181    Omega-          Omegabar+          -3    0    1      1.67245     0.00000     0.00000   2.46000E+01    1
          1202    1    0    0.676000    Lambda0         K-                                                              
          1203    1    0    0.234000    Xi0             pi-                                                             
          1204    1    0    0.085000    Xi-             pi0                                                             
          1205    1    0    0.005000    nu_ebar         e-              Xi0                                             

      4101    182    cd_0            cd_0bar             1   -1    1      1.96908     0.00000     0.00000   0.00000E+00    0

      4103    183    cd_1            cd_1bar             1   -1    1      2.00808     0.00000     0.00000   0.00000E+00    0

      4112    184    Sigma_c0        Sigma_cbar0         0    0    1      2.45210     0.00000     0.00000   0.00000E+00    1
          1206    1    0    1.000000    Lambda_c+       pi-                                                             

      4114    185    Sigma*_c0       Sigma*_cbar0        0    0    1      2.50000     0.00000     0.00000   0.00000E+00    1
          1207    1    0    1.000000    Lambda_c+       pi-                                                             

      4122    186    Lambda_c+       Lambda_cbar-        3    0    1      2.28490     0.00000     0.00000   6.18000E-02    1
          1208    0   42    0.018000    e+              nu_e            Lambda0                                         
          1209    0   42    0.005000    e+              nu_e            Sigma0                                          
          1210    0   42    0.005000    e+              nu_e            Sigma*0                                         
          1211    0   42    0.003000    e+              nu_e            n0                                              
          1212    0   42    0.002000    e+              nu_e            Delta0                                          
          1213    0   42    0.006000    e+              nu_e            p+              pi-                             
          1214    0   42    0.006000    e+              nu_e            n0              pi0                             
          1215    1   42    0.018000    mu+             nu_mu           Lambda0                                         
          1216    1   42    0.005000    mu+             nu_mu           Sigma0                                          
          1217    1   42    0.005000    mu+             nu_mu           Sigma*0                                         
          1218    1   42    0.003000    mu+             nu_mu           n0                                              
          1219    1   42    0.002000    mu+             nu_mu           Delta0                                          
          1220    1   42    0.006000    mu+             nu_mu           p+              pi-                             
          1221    1   42    0.006000    mu+             nu_mu           n0              pi0                             
          1222    0    0    0.006600    Delta++         K-                                                              
          1223    0    0    0.025000    Delta++         K*-                                                             
          1224    0    0    0.016000    p+              Kbar0                                                           
          1225    0    0    0.008800    p+              K*bar0                                                          
          1226    0    0    0.005000    Delta+          Kbar0                                                           
          1227    0    0    0.005000    Delta+          K*bar0                                                          
          1228    0    0    0.005800    Lambda0         pi+                                                             
          1229    0    0    0.005000    Lambda0         rho+                                                            
          1230    0    0    0.005500    Sigma0          pi+                                                             
          1231    0    0    0.004000    Sigma0          rho+                                                            
          1232    0    0    0.004000    Sigma*0         pi+                                                             
          1233    0    0    0.004000    Sigma*0         rho+                                                            
          1234    0    0    0.004000    Sigma+          pi0                                                             
          1235    0    0    0.002000    Sigma+          eta                                                             
          1236    0    0    0.002000    Sigma+          eta'                                                            
          1237    0    0    0.004000    Sigma+          rho0                                                            
          1238    0    0    0.004000    Sigma+          omega                                                           
          1239    0    0    0.003000    Sigma*+         pi0                                                             
          1240    0    0    0.002000    Sigma*+         eta                                                             
          1241    0    0    0.003000    Sigma*+         rho0                                                            
          1242    0    0    0.003000    Sigma*+         omega                                                           
          1243    0    0    0.002000    Xi0             K+                                                              
          1244    0    0    0.002000    Xi0             K*+                                                             
          1245    0    0    0.002000    Xi*0            K+                                                              
          1246    0    0    0.001000    Delta++         pi-                                                             
          1247    0    0    0.001000    Delta++         rho-                                                            
          1248    0    0    0.002000    p+              pi0                                                             
          1249    0    0    0.001000    p+              eta                                                             
          1250    0    0    0.001000    p+              eta'                                                            
          1251    0    0    0.002000    p+              rho0                                                            
          1252    0    0    0.002000    p+              omega                                                           
          1253    0    0    0.001300    p+              phi                                                             
          1254    0    0    0.001800    p+              f_0                                                             
          1255    0    0    0.001000    Delta+          pi0                                                             
          1256    0    0    0.001000    Delta+          eta                                                             
          1257    0    0    0.001000    Delta+          eta'                                                            
          1258    0    0    0.001000    Delta+          rho0                                                            
          1259    0    0    0.001000    Delta+          omega                                                           
          1260    0    0    0.003000    n0              pi+                                                             
          1261    0    0    0.003000    n0              rho+                                                            
          1262    0    0    0.003000    Delta0          pi+                                                             
          1263    0    0    0.003000    Delta0          rho+                                                            
          1264    0    0    0.005000    Lambda0         K+                                                              
          1265    0    0    0.005000    Lambda0         K*+                                                             
          1266    0    0    0.002000    Sigma0          K+                                                              
          1267    0    0    0.002000    Sigma0          K*+                                                             
          1268    0    0    0.001000    Sigma*0         K+                                                              
          1269    0    0    0.001000    Sigma*0         K*+                                                             
          1270    0    0    0.002000    Sigma+          K0                                                              
          1271    0    0    0.002000    Sigma+          K*0                                                             
          1272    0    0    0.001000    Sigma*+         K0                                                              
          1273    0    0    0.001000    Sigma*+         K*0                                                             
          1274    0   13    0.243200    u               dbar            s               ud_0                            
          1275    0   13    0.057000    u               dbar            s               ud_1                            
          1276    0   13    0.035000    u               sbar            s               ud_0                            
          1277    0   13    0.035000    u               dbar            d               ud_0                            
          1278    0   13    0.150000    s               uu_1                                                            
          1279    0   13    0.075000    u               su_0                                                            
          1280    0   13    0.075000    u               su_1                                                            
          1281    0   13    0.030000    d               uu_1                                                            
          1282    0   13    0.015000    u               ud_0                                                            
          1283    0   13    0.015000    u               ud_1                                                            

      4132    187    Xi_c0           Xi_cbar0            0    0    1      2.47030     0.00000     0.00000   2.90000E-02    1
          1284    1   42    0.080000    e+              nu_e            s               specflav                        
          1285    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1286    1   11    0.760000    u               dbar            s               specflav                        
          1287    1   11    0.080000    u               sbar            s               specflav                        

      4201    188    cu_0            cu_0bar             4   -1    1      1.96908     0.00000     0.00000   0.00000E+00    0

      4203    189    cu_1            cu_1bar             4   -1    1      2.00808     0.00000     0.00000   0.00000E+00    0

      4212    190    Sigma_c+        Sigma_cbar-         3    0    1      2.45350     0.00000     0.00000   0.00000E+00    1
          1288    1    0    1.000000    Lambda_c+       pi0                                                             

      4214    191    Sigma*_c+       Sigma*_cbar-        3    0    1      2.50000     0.00000     0.00000   0.00000E+00    1
          1289    1    0    1.000000    Lambda_c+       pi0                                                             

      4222    192    Sigma_c++       Sigma_cbar--        6    0    1      2.45290     0.00000     0.00000   0.00000E+00    1
          1290    1    0    1.000000    Lambda_c+       pi+                                                             

      4224    193    Sigma*_c++      Sigma*_cbar--       6    0    1      2.50000     0.00000     0.00000   0.00000E+00    1
          1291    1    0    1.000000    Lambda_c+       pi+                                                             

      4232    194    Xi_c+           Xi_cbar-            3    0    1      2.46560     0.00000     0.00000   1.06000E-01    1
          1292    1   42    0.080000    e+              nu_e            s               specflav                        
          1293    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1294    1   11    0.760000    u               dbar            s               specflav                        
          1295    1   11    0.080000    u               sbar            s               specflav                        

      4301    195    cs_0            cs_0bar             1   -1    1      2.15432     0.00000     0.00000   0.00000E+00    0

      4303    196    cs_1            cs_1bar             1   -1    1      2.17967     0.00000     0.00000   0.00000E+00    0

      4312    197    Xi'_c0          Xi'_cbar0           0    0    1      2.55000     0.00000     0.00000   0.00000E+00    1
          1296    1    0    1.000000    Xi_c0           gamma                                                           

      4314    198    Xi*_c0          Xi*_cbar0           0    0    1      2.63000     0.00000     0.00000   0.00000E+00    1
          1297    1    0    0.500000    Xi_c0           pi0                                                             
          1298    1    0    0.500000    Xi_c0           gamma                                                           

      4322    199    Xi'_c+          Xi'_cbar-           3    0    1      2.55000     0.00000     0.00000   0.00000E+00    1
          1299    1    0    1.000000    Xi_c+           gamma                                                           

      4324    200    Xi*_c+          Xi*_cbar-           3    0    1      2.63000     0.00000     0.00000   0.00000E+00    1
          1300    1    0    0.500000    Xi_c+           pi0                                                             
          1301    1    0    0.500000    Xi_c+           gamma                                                           

      4332    201    Omega_c0        Omega_cbar0         0    0    1      2.70400     0.00000     0.00000   1.90000E-02    1
          1302    1   42    0.080000    e+              nu_e            s               specflav                        
          1303    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1304    1   11    0.760000    u               dbar            s               specflav                        
          1305    1   11    0.080000    u               sbar            s               specflav                        

      4334    202    Omega*_c0       Omega*_cbar0        0    0    1      2.80000     0.00000     0.00000   0.00000E+00    1
          1306    1    0    1.000000    Omega_c0        gamma                                                           

      4403    203    cc_1            cc_1bar             4   -1    1      3.27531     0.00000     0.00000   0.00000E+00    0

      4412    204    Xi_cc+          Xi_ccbar-           3    0    1      3.59798     0.00000     0.00000   1.00000E-01    1
          1307    1   42    0.080000    e+              nu_e            s               specflav                        
          1308    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1309    1   11    0.760000    u               dbar            s               specflav                        
          1310    1   11    0.080000    u               sbar            s               specflav                        

      4414    205    Xi*_cc+         Xi*_ccbar-          3    0    1      3.65648     0.00000     0.00000   1.00000E-01    1
          1311    1   42    0.080000    e+              nu_e            s               specflav                        
          1312    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1313    1   11    0.760000    u               dbar            s               specflav                        
          1314    1   11    0.080000    u               sbar            s               specflav                        

      4422    206    Xi_cc++         Xi_ccbar--          6    0    1      3.59798     0.00000     0.00000   1.00000E-01    1
          1315    1   42    0.080000    e+              nu_e            s               specflav                        
          1316    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1317    1   11    0.760000    u               dbar            s               specflav                        
          1318    1   11    0.080000    u               sbar            s               specflav                        

      4424    207    Xi*_cc++        Xi*_ccbar--         6    0    1      3.65648     0.00000     0.00000   1.00000E-01    1
          1319    1   42    0.080000    e+              nu_e            s               specflav                        
          1320    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1321    1   11    0.760000    u               dbar            s               specflav                        
          1322    1   11    0.080000    u               sbar            s               specflav                        

      4432    208    Omega_cc+       Omega_ccbar-        3    0    1      3.78663     0.00000     0.00000   1.00000E-01    1
          1323    1   42    0.080000    e+              nu_e            s               specflav                        
          1324    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1325    1   11    0.760000    u               dbar            s               specflav                        
          1326    1   11    0.080000    u               sbar            s               specflav                        

      4434    209    Omega*_cc+      Omega*_ccbar-       3    0    1      3.82466     0.00000     0.00000   1.00000E-01    1
          1327    1   42    0.080000    e+              nu_e            s               specflav                        
          1328    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1329    1   11    0.760000    u               dbar            s               specflav                        
          1330    1   11    0.080000    u               sbar            s               specflav                        

      4444    210    Omega*_ccc++    Omega*_cccbar-      6    0    1      4.91594     0.00000     0.00000   1.00000E-01    1
          1331    1   42    0.080000    e+              nu_e            s               specflav                        
          1332    1   42    0.080000    mu+             nu_mu           s               specflav                        
          1333    1   11    0.760000    u               dbar            s               specflav                        
          1334    1   11    0.080000    u               sbar            s               specflav                        

      5101    211    bd_0            bd_0bar            -2   -1    1      5.38897     0.00000     0.00000   0.00000E+00    0

      5103    212    bd_1            bd_1bar            -2   -1    1      5.40145     0.00000     0.00000   0.00000E+00    0

      5112    213    Sigma_b-        Sigma_bbar+        -3    0    1      5.80000     0.00000     0.00000   0.00000E+00    1
          1335    1    0    1.000000    Lambda_b0       pi-                                                             

      5114    214    Sigma*_b-       Sigma*_bbar+       -3    0    1      5.81000     0.00000     0.00000   0.00000E+00    1
          1336    1    0    1.000000    Lambda_b0       pi-                                                             

      5122    215    Lambda_b0       Lambda_bbar0        0    0    1      5.64100     0.00000     0.00000   3.42000E-01    1
          1337    0   42    0.105000    nu_ebar         e-              Lambda_c+                                       
          1338    1   42    0.105000    nu_mubar        mu-             Lambda_c+                                       
          1339    0   42    0.040000    nu_taubar       tau-            Lambda_c+                                       
          1340    0    0    0.007700    Lambda_c+       pi-                                                             
          1341    0    0    0.020000    Lambda_c+       rho-                                                            
          1342    0    0    0.023500    Lambda_c+       a_1-                                                            
          1343    0    0    0.028500    Lambda_c+       D_s-                                                            
          1344    0    0    0.043500    Lambda_c+       D*_s-                                                           
          1345    0    0    0.001100    eta_c           Lambda0                                                         
          1346    0    0    0.002200    J/psi           Lambda0                                                         
          1347    0    0    0.004400    chi_1c          Lambda0                                                         
          1348    0   48    0.429100    ubar            d               c               ud_0                            
          1349    0   13    0.080000    ubar            c               d               ud_0                            
          1350    0   13    0.070000    cbar            s               c               ud_0                            
          1351    0   13    0.020000    cbar            c               s               ud_0                            
          1352    0   42    0.015000    ubar            d               u               ud_0                            
          1353    0   42    0.005000    cbar            s               u               ud_0                            

      5132    216    Xi_b-           Xi_bbar+           -3    0    1      5.84000     0.00000     0.00000   3.87000E-01    1
          1354    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1355    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1356    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1357    1   42    0.500000    ubar            d               c               specflav                        
          1358    1   42    0.080000    ubar            c               d               specflav                        
          1359    1   42    0.140000    cbar            s               c               specflav                        
          1360    1   42    0.010000    cbar            c               s               specflav                        
          1361    1   42    0.015000    ubar            d               u               specflav                        
          1362    1   42    0.005000    cbar            s               u               specflav                        

      5142    217    Xi_bc0          Xi_bcbar0           0    0    1      7.00575     0.00000     0.00000   3.87000E-01    1
          1363    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1364    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1365    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1366    1   42    0.500000    ubar            d               c               specflav                        
          1367    1   42    0.080000    ubar            c               d               specflav                        
          1368    1   42    0.140000    cbar            s               c               specflav                        
          1369    1   42    0.010000    cbar            c               s               specflav                        
          1370    1   42    0.015000    ubar            d               u               specflav                        
          1371    1   42    0.005000    cbar            s               u               specflav                        

      5201    218    bu_0            bu_0bar             1   -1    1      5.38897     0.00000     0.00000   0.00000E+00    0

      5203    219    bu_1            bu_1bar             1   -1    1      5.40145     0.00000     0.00000   0.00000E+00    0

      5212    220    Sigma_b0        Sigma_bbar0         0    0    1      5.80000     0.00000     0.00000   0.00000E+00    1
          1372    1    0    1.000000    Lambda_b0       pi0                                                             

      5214    221    Sigma*_b0       Sigma*_bbar0        0    0    1      5.81000     0.00000     0.00000   0.00000E+00    1
          1373    1    0    1.000000    Lambda_b0       pi0                                                             

      5222    222    Sigma_b+        Sigma_bbar-         3    0    1      5.80000     0.00000     0.00000   0.00000E+00    1
          1374    1    0    1.000000    Lambda_b0       pi+                                                             

      5224    223    Sigma*_b+       Sigma*_bbar-        3    0    1      5.81000     0.00000     0.00000   0.00000E+00    1
          1375    1    0    1.000000    Lambda_b0       pi+                                                             

      5232    224    Xi_b0           Xi_bbar0            0    0    1      5.84000     0.00000     0.00000   3.87000E-01    1
          1376    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1377    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1378    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1379    1   42    0.500000    ubar            d               c               specflav                        
          1380    1   42    0.080000    ubar            c               d               specflav                        
          1381    1   42    0.140000    cbar            s               c               specflav                        
          1382    1   42    0.010000    cbar            c               s               specflav                        
          1383    1   42    0.015000    ubar            d               u               specflav                        
          1384    1   42    0.005000    cbar            s               u               specflav                        

      5242    225    Xi_bc+          Xi_bcbar-           3    0    1      7.00575     0.00000     0.00000   3.87000E-01    1
          1385    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1386    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1387    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1388    1   42    0.500000    ubar            d               c               specflav                        
          1389    1   42    0.080000    ubar            c               d               specflav                        
          1390    1   42    0.140000    cbar            s               c               specflav                        
          1391    1   42    0.010000    cbar            c               s               specflav                        
          1392    1   42    0.015000    ubar            d               u               specflav                        
          1393    1   42    0.005000    cbar            s               u               specflav                        

      5301    226    bs_0            bs_0bar            -2   -1    1      5.56725     0.00000     0.00000   0.00000E+00    0

      5303    227    bs_1            bs_1bar            -2   -1    1      5.57536     0.00000     0.00000   0.00000E+00    0

      5312    228    Xi'_b-          Xi'_bbar+          -3    0    1      5.96000     0.00000     0.00000   0.00000E+00    1
          1394    1    0    1.000000    Xi_b-           gamma                                                           

      5314    229    Xi*_b-          Xi*_bbar+          -3    0    1      5.97000     0.00000     0.00000   0.00000E+00    1
          1395    1    0    1.000000    Xi_b-           gamma                                                           

      5322    230    Xi'_b0          Xi'_bbar0           0    0    1      5.96000     0.00000     0.00000   0.00000E+00    1
          1396    1    0    1.000000    Xi_b0           gamma                                                           

      5324    231    Xi*_b0          Xi*_bbar0           0    0    1      5.97000     0.00000     0.00000   0.00000E+00    1
          1397    1    0    1.000000    Xi_b0           gamma                                                           

      5332    232    Omega_b-        Omega_bbar+        -3    0    1      6.12000     0.00000     0.00000   3.87000E-01    1
          1398    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1399    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1400    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1401    1   42    0.500000    ubar            d               c               specflav                        
          1402    1   42    0.080000    ubar            c               d               specflav                        
          1403    1   42    0.140000    cbar            s               c               specflav                        
          1404    1   42    0.010000    cbar            c               s               specflav                        
          1405    1   42    0.015000    ubar            d               u               specflav                        
          1406    1   42    0.005000    cbar            s               u               specflav                        

      5334    233    Omega*_b-       Omega*_bbar+       -3    0    1      6.13000     0.00000     0.00000   0.00000E+00    1
          1407    1    0    1.000000    Omega_b-        gamma                                                           

      5342    234    Omega_bc0       Omega_bcbar0        0    0    1      7.19099     0.00000     0.00000   3.87000E-01    1
          1408    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1409    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1410    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1411    1   42    0.500000    ubar            d               c               specflav                        
          1412    1   42    0.080000    ubar            c               d               specflav                        
          1413    1   42    0.140000    cbar            s               c               specflav                        
          1414    1   42    0.010000    cbar            c               s               specflav                        
          1415    1   42    0.015000    ubar            d               u               specflav                        
          1416    1   42    0.005000    cbar            s               u               specflav                        

      5401    235    bc_0            bc_0bar             1   -1    1      6.67143     0.00000     0.00000   0.00000E+00    0

      5403    236    bc_1            bc_1bar             1   -1    1      6.67397     0.00000     0.00000   0.00000E+00    0

      5412    237    Xi'_bc0         Xi'_bcbar0          0    0    1      7.03724     0.00000     0.00000   3.87000E-01    1
          1417    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1418    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1419    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1420    1   42    0.500000    ubar            d               c               specflav                        
          1421    1   42    0.080000    ubar            c               d               specflav                        
          1422    1   42    0.140000    cbar            s               c               specflav                        
          1423    1   42    0.010000    cbar            c               s               specflav                        
          1424    1   42    0.015000    ubar            d               u               specflav                        
          1425    1   42    0.005000    cbar            s               u               specflav                        

      5414    238    Xi*_bc0         Xi*_bcbar0          0    0    1      7.04850     0.00000     0.00000   3.87000E-01    1
          1426    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1427    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1428    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1429    1   42    0.500000    ubar            d               c               specflav                        
          1430    1   42    0.080000    ubar            c               d               specflav                        
          1431    1   42    0.140000    cbar            s               c               specflav                        
          1432    1   42    0.010000    cbar            c               s               specflav                        
          1433    1   42    0.015000    ubar            d               u               specflav                        
          1434    1   42    0.005000    cbar            s               u               specflav                        

      5422    239    Xi'_bc+         Xi'_bcbar-          3    0    1      7.03724     0.00000     0.00000   3.87000E-01    1
          1435    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1436    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1437    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1438    1   42    0.500000    ubar            d               c               specflav                        
          1439    1   42    0.080000    ubar            c               d               specflav                        
          1440    1   42    0.140000    cbar            s               c               specflav                        
          1441    1   42    0.010000    cbar            c               s               specflav                        
          1442    1   42    0.015000    ubar            d               u               specflav                        
          1443    1   42    0.005000    cbar            s               u               specflav                        

      5424    240    Xi*_bc+         Xi*_bcbar-          3    0    1      7.04850     0.00000     0.00000   3.87000E-01    1
          1444    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1445    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1446    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1447    1   42    0.500000    ubar            d               c               specflav                        
          1448    1   42    0.080000    ubar            c               d               specflav                        
          1449    1   42    0.140000    cbar            s               c               specflav                        
          1450    1   42    0.010000    cbar            c               s               specflav                        
          1451    1   42    0.015000    ubar            d               u               specflav                        
          1452    1   42    0.005000    cbar            s               u               specflav                        

      5432    241    Omega'_bc0      Omega'_bcba         0    0    1      7.21101     0.00000     0.00000   3.87000E-01    1
          1453    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1454    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1455    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1456    1   42    0.500000    ubar            d               c               specflav                        
          1457    1   42    0.080000    ubar            c               d               specflav                        
          1458    1   42    0.140000    cbar            s               c               specflav                        
          1459    1   42    0.010000    cbar            c               s               specflav                        
          1460    1   42    0.015000    ubar            d               u               specflav                        
          1461    1   42    0.005000    cbar            s               u               specflav                        

      5434    242    Omega*_bc0      Omega*_bcbar0       0    0    1      7.21900     0.00000     0.00000   3.87000E-01    1
          1462    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1463    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1464    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1465    1   42    0.500000    ubar            d               c               specflav                        
          1466    1   42    0.080000    ubar            c               d               specflav                        
          1467    1   42    0.140000    cbar            s               c               specflav                        
          1468    1   42    0.010000    cbar            c               s               specflav                        
          1469    1   42    0.015000    ubar            d               u               specflav                        
          1470    1   42    0.005000    cbar            s               u               specflav                        

      5442    243    Omega_bcc+      Omega_bccbar-       3    0    1      8.30945     0.00000     0.00000   3.87000E-01    1
          1471    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1472    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1473    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1474    1   42    0.500000    ubar            d               c               specflav                        
          1475    1   42    0.080000    ubar            c               d               specflav                        
          1476    1   42    0.140000    cbar            s               c               specflav                        
          1477    1   42    0.010000    cbar            c               s               specflav                        
          1478    1   42    0.015000    ubar            d               u               specflav                        
          1479    1   42    0.005000    cbar            s               u               specflav                        

      5444    244    Omega*_bcc+     Omega*_bccbar-      3    0    1      8.31325     0.00000     0.00000   3.87000E-01    1
          1480    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1481    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1482    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1483    1   42    0.500000    ubar            d               c               specflav                        
          1484    1   42    0.080000    ubar            c               d               specflav                        
          1485    1   42    0.140000    cbar            s               c               specflav                        
          1486    1   42    0.010000    cbar            c               s               specflav                        
          1487    1   42    0.015000    ubar            d               u               specflav                        
          1488    1   42    0.005000    cbar            s               u               specflav                        

      5503    245    bb_1            bb_1bar            -2   -1    1     10.07354     0.00000     0.00000   0.00000E+00    0

      5512    246    Xi_bb-          Xi_bbbar+          -3    0    1     10.42272     0.00000     0.00000   3.87000E-01    1
          1489    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1490    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1491    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1492    1   42    0.500000    ubar            d               c               specflav                        
          1493    1   42    0.080000    ubar            c               d               specflav                        
          1494    1   42    0.140000    cbar            s               c               specflav                        
          1495    1   42    0.010000    cbar            c               s               specflav                        
          1496    1   42    0.015000    ubar            d               u               specflav                        
          1497    1   42    0.005000    cbar            s               u               specflav                        

      5514    247    Xi*_bb-         Xi*_bbbar+         -3    0    1     10.44144     0.00000     0.00000   3.87000E-01    1
          1498    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1499    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1500    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1501    1   42    0.500000    ubar            d               c               specflav                        
          1502    1   42    0.080000    ubar            c               d               specflav                        
          1503    1   42    0.140000    cbar            s               c               specflav                        
          1504    1   42    0.010000    cbar            c               s               specflav                        
          1505    1   42    0.015000    ubar            d               u               specflav                        
          1506    1   42    0.005000    cbar            s               u               specflav                        

      5522    248    Xi_bb0          Xi_bbbar0           0    0    1     10.42272     0.00000     0.00000   3.87000E-01    1
          1507    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1508    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1509    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1510    1   42    0.500000    ubar            d               c               specflav                        
          1511    1   42    0.080000    ubar            c               d               specflav                        
          1512    1   42    0.140000    cbar            s               c               specflav                        
          1513    1   42    0.010000    cbar            c               s               specflav                        
          1514    1   42    0.015000    ubar            d               u               specflav                        
          1515    1   42    0.005000    cbar            s               u               specflav                        

      5524    249    Xi*_bb0         Xi*_bbbar0          0    0    1     10.44144     0.00000     0.00000   3.87000E-01    1
          1516    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1517    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1518    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1519    1   42    0.500000    ubar            d               c               specflav                        
          1520    1   42    0.080000    ubar            c               d               specflav                        
          1521    1   42    0.140000    cbar            s               c               specflav                        
          1522    1   42    0.010000    cbar            c               s               specflav                        
          1523    1   42    0.015000    ubar            d               u               specflav                        
          1524    1   42    0.005000    cbar            s               u               specflav                        

      5532    250    Omega_bb-       Omega_bbbar+       -3    0    1     10.60209     0.00000     0.00000   3.87000E-01    1
          1525    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1526    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1527    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1528    1   42    0.500000    ubar            d               c               specflav                        
          1529    1   42    0.080000    ubar            c               d               specflav                        
          1530    1   42    0.140000    cbar            s               c               specflav                        
          1531    1   42    0.010000    cbar            c               s               specflav                        
          1532    1   42    0.015000    ubar            d               u               specflav                        
          1533    1   42    0.005000    cbar            s               u               specflav                        

      5534    251    Omega*_bb-      Omega*_bbbar+      -3    0    1     10.61426     0.00000     0.00000   3.87000E-01    1
          1534    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1535    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1536    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1537    1   42    0.500000    ubar            d               c               specflav                        
          1538    1   42    0.080000    ubar            c               d               specflav                        
          1539    1   42    0.140000    cbar            s               c               specflav                        
          1540    1   42    0.010000    cbar            c               s               specflav                        
          1541    1   42    0.015000    ubar            d               u               specflav                        
          1542    1   42    0.005000    cbar            s               u               specflav                        

      5542    252    Omega_bbc0      Omega_bbcbar0       0    0    1     11.70767     0.00000     0.00000   3.87000E-01    1
          1543    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1544    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1545    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1546    1   42    0.500000    ubar            d               c               specflav                        
          1547    1   42    0.080000    ubar            c               d               specflav                        
          1548    1   42    0.140000    cbar            s               c               specflav                        
          1549    1   42    0.010000    cbar            c               s               specflav                        
          1550    1   42    0.015000    ubar            d               u               specflav                        
          1551    1   42    0.005000    cbar            s               u               specflav                        

      5544    253    Omega*_bbc0     Omega*_bbcbar0      0    0    1     11.71147     0.00000     0.00000   3.87000E-01    1
          1552    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1553    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1554    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1555    1   42    0.500000    ubar            d               c               specflav                        
          1556    1   42    0.080000    ubar            c               d               specflav                        
          1557    1   42    0.140000    cbar            s               c               specflav                        
          1558    1   42    0.010000    cbar            c               s               specflav                        
          1559    1   42    0.015000    ubar            d               u               specflav                        
          1560    1   42    0.005000    cbar            s               u               specflav                        

      5554    254    Omega*_bbb-     Omega*_bbbbar+     -3    0    1     15.11061     0.00000     0.00000   3.87000E-01    1
          1561    1   42    0.105000    nu_ebar         e-              c               specflav                        
          1562    1   42    0.105000    nu_mubar        mu-             c               specflav                        
          1563    1   42    0.040000    nu_taubar       tau-            c               specflav                        
          1564    1   42    0.500000    ubar            d               c               specflav                        
          1565    1   42    0.080000    ubar            c               d               specflav                        
          1566    1   42    0.140000    cbar            s               c               specflav                        
          1567    1   42    0.010000    cbar            c               s               specflav                        
          1568    1   42    0.015000    ubar            d               u               specflav                        
          1569    1   42    0.005000    cbar            s               u               specflav                        

     10111    255    a_00                                0    0    0      0.98350     0.06000     0.05000   0.00000E+00    1
          1570    1    0    1.000000    eta             pi0                                                             

     10113    256    b_10                                0    0    0      1.23100     0.14200     0.25000   0.00000E+00    1
          1571    1    0    1.000000    omega           pi0                                                             

     10211    257    a_0+            a_0-                3    0    1      0.98350     0.06000     0.05000   0.00000E+00    1
          1572    1    0    1.000000    eta             pi+                                                             

     10213    258    b_1+            b_1-                3    0    1      1.23100     0.14200     0.25000   0.00000E+00    1
          1573    1    0    1.000000    omega           pi+                                                             

     10221    259    f_0                                 0    0    0      1.00000     0.00000     0.00000   0.00000E+00    1
          1574    1    0    0.520000    pi+             pi-                                                             
          1575    1    0    0.260000    pi0             pi0                                                             
          1576    1    0    0.110000    K+              K-                                                              
          1577    1    0    0.055000    K_L0            K_L0                                                            
          1578    1    0    0.055000    K_S0            K_S0                                                            

     10223    260    h_1                                 0    0    0      1.17000     0.36000     0.20000   0.00000E+00    1
          1579    1    0    0.333000    rho+            pi-                                                             
          1580    1    0    0.334000    rho0            pi0                                                             
          1581    1    0    0.333000    rho-            pi+                                                             

     10311    261    K*_00           K*_0bar0            0    0    1      1.42900     0.28700     0.40000   0.00000E+00    1
          1582    1    0    0.667000    K+              pi-                                                             
          1583    1    0    0.333000    K0              pi0                                                             

     10313    262    K_10            K_1bar0             0    0    1      1.29000     0.09000     0.00500   0.00000E+00    1
          1584    1    0    0.280000    K+              rho-                                                            
          1585    1    0    0.140000    K0              rho0                                                            
          1586    1    0    0.313000    K*+             pi-                                                             
          1587    1    0    0.157000    K*0             pi0                                                             
          1588    1    0    0.110000    K0              omega                                                           

     10321    263    K*_0+           K*_0-               3    0    1      1.42900     0.28700     0.40000   0.00000E+00    1
          1589    1    0    0.667000    K0              pi+                                                             
          1590    1    0    0.333000    K+              pi0                                                             

     10323    264    K_1+            K_1-                3    0    1      1.29000     0.09000     0.01000   0.00000E+00    1
          1591    1    0    0.280000    K0              rho+                                                            
          1592    1    0    0.140000    K+              rho0                                                            
          1593    1    0    0.313000    K*0             pi+                                                             
          1594    1    0    0.157000    K*+             pi0                                                             
          1595    1    0    0.110000    K+              omega                                                           

     10331    265    f'_0                                0    0    0      1.40000     0.25000     0.35000   0.00000E+00    1
          1596    1    0    0.360000    pi+             pi-                                                             
          1597    1    0    0.180000    pi0             pi0                                                             
          1598    1    0    0.030000    K+              K-                                                              
          1599    1    0    0.015000    K_L0            K_L0                                                            
          1600    1    0    0.015000    K_S0            K_S0                                                            
          1601    1    0    0.200000    pi+             pi-             pi+             pi-                             
          1602    1    0    0.200000    pi+             pi-             pi0             pi0                             

     10333    266    h'_1                                0    0    0      1.40000     0.08000     0.00100   0.00000E+00    1
          1603    1    0    0.250000    K*0             Kbar0                                                           
          1604    1    0    0.250000    K*bar0          K0                                                              
          1605    1    0    0.250000    K*+             K-                                                              
          1606    1    0    0.250000    K*-             K+                                                              

     10411    267    D*_0+           D*_0-               3    0    1      2.27200     0.05000     0.10000   0.00000E+00    1
          1607    1    0    0.667000    D0              pi+                                                             
          1608    1    0    0.333000    D+              pi0                                                             

     10413    268    D_1+            D_1-                3    0    1      2.42400     0.02000     0.08000   0.00000E+00    1
          1609    1    0    0.667000    D*0             pi+                                                             
          1610    1    0    0.333000    D*+             pi0                                                             

     10421    269    D*_00           D*_0bar0            0    0    1      2.27200     0.05000     0.10000   0.00000E+00    1
          1611    1    0    0.667000    D+              pi-                                                             
          1612    1    0    0.333000    D0              pi0                                                             

     10423    270    D_10            D_1bar0             0    0    1      2.42400     0.02000     0.08000   0.00000E+00    1
          1613    1    0    0.667000    D*+             pi-                                                             
          1614    1    0    0.333000    D*0             pi0                                                             

     10431    271    D*_0s+          D*_0s-              3    0    1      2.50000     0.05000     0.10000   0.00000E+00    1
          1615    1    0    0.500000    D+              K0                                                              
          1616    1    0    0.500000    D0              K+                                                              

     10433    272    D_1s+           D_1s-               3    0    1      2.53600     0.00000     0.00000   0.00000E+00    1
          1617    1    0    0.500000    D*0             K+                                                              
          1618    1    0    0.500000    D*+             K0                                                              

     10441    273    chi_0c                              0    0    0      3.41510     0.01400     0.05000   0.00000E+00    1
          1619    1    0    0.007000    J/psi           gamma                                                           
          1620    1   12    0.993000    rndmflav        rndmflavbar                                                     

     10443    274    h_1c                                0    0    0      3.46000     0.01000     0.02000   0.00000E+00    1
          1621    1   12    1.000000    rndmflav        rndmflavbar                                                     

     10511    275    B*_00           B*_0bar0            0    0    1      5.68000     0.05000     0.10000   0.00000E+00    1
          1622    1    0    0.667000    B+              pi-                                                             
          1623    1    0    0.333000    B0              pi0                                                             

     10513    276    B_10            B_1bar0             0    0    1      5.73000     0.05000     0.10000   0.00000E+00    1
          1624    1    0    0.667000    B*+             pi-                                                             
          1625    1    0    0.333000    B*0             pi0                                                             

     10521    277    B*_0+           B*_0-               3    0    1      5.68000     0.05000     0.10000   0.00000E+00    1
          1626    1    0    0.667000    B0              pi+                                                             
          1627    1    0    0.333000    B+              pi0                                                             

     10523    278    B_1+            B_1-                3    0    1      5.73000     0.05000     0.10000   0.00000E+00    1
          1628    1    0    0.667000    B*0             pi+                                                             
          1629    1    0    0.333000    B*+             pi0                                                             

     10531    279    B*_0s0          B*_0sbar0           0    0    1      5.92000     0.05000     0.10000   0.00000E+00    1
          1630    1    0    0.500000    B+              K-                                                              
          1631    1    0    0.500000    B0              Kbar0                                                           

     10533    280    B_1s0           B_1sbar0            0    0    1      5.97000     0.05000     0.10000   0.00000E+00    1
          1632    1    0    0.500000    B*+             K-                                                              
          1633    1    0    0.500000    B*0             Kbar0                                                           

     10541    281    B*_0c+          B*_0c-              3    0    1      7.25000     0.05000     0.05000   0.00000E+00    1
          1634    1    0    0.500000    B0              D+                                                              
          1635    1    0    0.500000    B+              D0                                                              

     10543    282    B_1c+           B_1c-               3    0    1      7.30000     0.05000     0.10000   0.00000E+00    1
          1636    1    0    0.500000    B*0             D+                                                              
          1637    1    0    0.500000    B*+             D0                                                              

     10551    283    chi_0b                              0    0    0      9.85980     0.00000     0.00000   0.00000E+00    1
          1638    1    0    0.020000    Upsilon         gamma                                                           
          1639    1   32    0.980000    g               g                                                               

     10553    284    h_1b                                0    0    0      9.87500     0.01000     0.02000   0.00000E+00    1
          1640    1   32    1.000000    g               g                                                               

     20113    285    a_10                                0    0    0      1.23000     0.40000     0.30000   0.00000E+00    1
          1641    1    0    0.500000    rho+            pi-                                                             
          1642    1    0    0.500000    rho-            pi+                                                             

     20213    286    a_1+            a_1-                3    0    1      1.23000     0.40000     0.30000   0.00000E+00    1
          1643    1    0    0.500000    rho0            pi+                                                             
          1644    1    0    0.500000    rho+            pi0                                                             

     20223    287    f_1                                 0    0    0      1.28200     0.02500     0.05000   0.00000E+00    1
          1645    1    0    0.146000    a_0+            pi-                                                             
          1646    1    0    0.146000    a_00            pi0                                                             
          1647    1    0    0.146000    a_0-            pi+                                                             
          1648    1    0    0.050000    eta             pi+             pi-                                             
          1649    1    0    0.050000    eta             pi0             pi0                                             
          1650    1    0    0.050000    rho+            pi-             pi0                                             
          1651    1    0    0.150000    rho0            pi+             pi-                                             
          1652    1    0    0.050000    rho0            pi0             pi0                                             
          1653    1    0    0.050000    rho-            pi+             pi0                                             
          1654    1    0    0.024000    K+              K-              pi0                                             
          1655    1    0    0.024000    K+              Kbar0           pi-                                             
          1656    1    0    0.024000    K0              Kbar0           pi0                                             
          1657    1    0    0.024000    K0              K-              pi+                                             
          1658    1    0    0.066000    rho0            gamma                                                           

     20313    288    K*_10           K*_1bar0            0    0    1      1.40200     0.17400     0.30000   0.00000E+00    1
          1659    1    0    0.667000    K*+             pi-                                                             
          1660    1    0    0.333000    K*0             pi0                                                             

     20323    289    K*_1+           K*_1-               3    0    1      1.40200     0.17400     0.30000   0.00000E+00    1
          1661    1    0    0.667000    K*0             pi+                                                             
          1662    1    0    0.333000    K*+             pi0                                                             

     20333    290    f'_1                                0    0    0      1.42700     0.05300     0.02000   0.00000E+00    1
          1663    1    0    0.250000    K*0             Kbar0                                                           
          1664    1    0    0.250000    K*bar0          K0                                                              
          1665    1    0    0.250000    K*+             K-                                                              
          1666    1    0    0.250000    K*-             K+                                                              

     20413    291    D*_1+           D*_1-               3    0    1      2.37200     0.05000     0.10000   0.00000E+00    1
          1667    1    0    0.667000    D*0             pi+                                                             
          1668    1    0    0.333000    D*+             pi0                                                             

     20423    292    D*_10           D*_1bar0            0    0    1      2.37200     0.05000     0.10000   0.00000E+00    1
          1669    1    0    0.667000    D*+             pi-                                                             
          1670    1    0    0.333000    D*0             pi0                                                             

     20433    293    D*_1s+          D*_1s-              3    0    1      2.56000     0.05000     0.03000   0.00000E+00    1
          1671    1    0    0.500000    D*0             K+                                                              
          1672    1    0    0.500000    D*+             K0                                                              

     20443    294    chi_1c                              0    0    0      3.51060     0.00090     0.00100   0.00000E+00    1
          1673    1    0    0.273000    J/psi           gamma                                                           
          1674    1   12    0.727000    rndmflav        rndmflavbar                                                     

     20513    295    B*_10           B*_1bar0            0    0    1      5.78000     0.05000     0.10000   0.00000E+00    1
          1675    1    0    0.667000    B*+             pi-                                                             
          1676    1    0    0.333000    B*0             pi0                                                             

     20523    296    B*_1+           B*_1-               3    0    1      5.78000     0.05000     0.10000   0.00000E+00    1
          1677    1    0    0.667000    B*0             pi+                                                             
          1678    1    0    0.333000    B*+             pi0                                                             

     20533    297    B*_1s0          B*_1sbar0           0    0    1      6.02000     0.05000     0.10000   0.00000E+00    1
          1679    1    0    0.500000    B*+             K-                                                              
          1680    1    0    0.500000    B*0             Kbar0                                                           

     20543    298    B*_1c+          B*_1c-              3    0    1      7.30000     0.05000     0.10000   0.00000E+00    1
          1681    1    0    0.500000    B*0             D+                                                              
          1682    1    0    0.500000    B*+             D0                                                              

     20553    299    chi_1b                              0    0    0      9.89190     0.00000     0.00000   0.00000E+00    1
          1683    1    0    0.350000    Upsilon         gamma                                                           
          1684    1   32    0.650000    g               g                                                               

    100443    300    psi'                                0    0    0      3.68600     0.00000     0.00000   0.00000E+00    1
          1685    1    0    0.008300    e-              e+                                                              
          1686    1    0    0.008300    mu-             mu+                                                             
          1687    1   12    0.186600    rndmflav        rndmflavbar                                                     
          1688    1    0    0.324000    J/psi           pi+             pi-                                             
          1689    1    0    0.184000    J/psi           pi0             pi0                                             
          1690    1    0    0.027000    J/psi           eta                                                             
          1691    1    0    0.001000    J/psi           pi0                                                             
          1692    1    0    0.093000    chi_0c          gamma                                                           
          1693    1    0    0.087000    chi_1c          gamma                                                           
          1694    1    0    0.078000    chi_2c          gamma                                                           
          1695    1    0    0.002800    eta_c           gamma                                                           

    100553    301    Upsilon'                            0    0    0     10.02330     0.00000     0.00000   0.00000E+00    1
          1696    1    0    0.014000    e-              e+                                                              
          1697    1    0    0.014000    mu-             mu+                                                             
          1698    1    0    0.014000    tau-            tau+                                                            
          1699    1   32    0.008000    d               dbar                                                            
          1700    1   32    0.024000    u               ubar                                                            
          1701    1   32    0.008000    s               sbar                                                            
          1702    1   32    0.024000    c               cbar                                                            
          1703    1    4    0.425000    g               g               g                                               
          1704    1    4    0.020000    gamma           g               g                                               
          1705    1    0    0.185000    Upsilon         pi+             pi-                                             
          1706    1    0    0.088000    Upsilon         pi0             pi0                                             
          1707    1    0    0.043000    chi_0b          gamma                                                           
          1708    1    0    0.067000    chi_1b          gamma                                                           
          1709    1    0    0.066000    chi_2b          gamma                                                           

   1000001    302    ~d_L            ~d_Lbar            -1    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1710    1   53    0.000000    ~gravitino      d                                                               
          1711    1   53    0.000000    ~chi_1-         u                                                               
          1712    1   53    0.000000    ~chi_2-         u                                                               
          1713    1   53    0.000000    ~chi_10         d                                                               
          1714    1   53    0.000000    ~chi_20         d                                                               
          1715    1   53    0.000000    ~chi_30         d                                                               
          1716    1   53    0.000000    ~chi_40         d                                                               
          1717    1   53    0.000000    ~u_L            W-                                                              
          1718    1   53    0.000000    ~u_R            W-                                                              
          1719    1   53    0.000000    ~u_L            H-                                                              
          1720    1   53    0.000000    ~u_R            H-                                                              
          1721    1   53    0.000000    ~g              d                                                               

   1000002    303    ~u_L            ~u_Lbar             2    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1722    1   53    0.000000    ~gravitino      u                                                               
          1723    1   53    0.000000    ~chi_1+         d                                                               
          1724    1   53    0.000000    ~chi_2+         d                                                               
          1725    1   53    0.000000    ~chi_10         u                                                               
          1726    1   53    0.000000    ~chi_20         u                                                               
          1727    1   53    0.000000    ~chi_30         u                                                               
          1728    1   53    0.000000    ~chi_40         u                                                               
          1729    1   53    0.000000    ~d_L            W+                                                              
          1730    1   53    0.000000    ~d_R            W+                                                              
          1731    1   53    0.000000    ~d_L            H+                                                              
          1732    1   53    0.000000    ~d_R            H+                                                              
          1733    1   53    0.000000    ~g              u                                                               

   1000003    304    ~s_L            ~s_Lbar            -1    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1734    1   53    0.000000    ~gravitino      s                                                               
          1735    1   53    0.000000    ~chi_1-         c                                                               
          1736    1   53    0.000000    ~chi_2-         c                                                               
          1737    1   53    0.000000    ~chi_10         s                                                               
          1738    1   53    0.000000    ~chi_20         s                                                               
          1739    1   53    0.000000    ~chi_30         s                                                               
          1740    1   53    0.000000    ~chi_40         s                                                               
          1741    1   53    0.000000    ~c_L            W-                                                              
          1742    1   53    0.000000    ~c_R            W-                                                              
          1743    1   53    0.000000    ~c_L            H-                                                              
          1744    1   53    0.000000    ~c_R            H-                                                              
          1745    1   53    0.000000    ~g              s                                                               

   1000004    305    ~c_L            ~c_Lbar             2    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1746    1   53    0.000000    ~gravitino      c                                                               
          1747    1   53    0.000000    ~chi_1+         s                                                               
          1748    1   53    0.000000    ~chi_2+         s                                                               
          1749    1   53    0.000000    ~chi_10         c                                                               
          1750    1   53    0.000000    ~chi_20         c                                                               
          1751    1   53    0.000000    ~chi_30         c                                                               
          1752    1   53    0.000000    ~chi_40         c                                                               
          1753    1   53    0.000000    ~s_L            W+                                                              
          1754    1   53    0.000000    ~s_R            W+                                                              
          1755    1   53    0.000000    ~s_L            H+                                                              
          1756    1   53    0.000000    ~s_R            H+                                                              
          1757    1   53    0.000000    ~g              c                                                               

   1000005    306    ~b_1            ~b_1bar            -1    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1758    1   53    0.000000    ~gravitino      b                                                               
          1759    1   53    0.000000    ~chi_1-         t                                                               
          1760    1   53    0.000000    ~chi_2-         t                                                               
          1761    1   53    0.000000    ~chi_10         b                                                               
          1762    1   53    0.000000    ~chi_20         b                                                               
          1763    1   53    0.000000    ~chi_30         b                                                               
          1764    1   53    0.000000    ~chi_40         b                                                               
          1765    1   53    0.000000    ~t_1            W-                                                              
          1766    1   53    0.000000    ~t_2            W-                                                              
          1767    1   53    0.000000    ~t_1            H-                                                              
          1768    1   53    0.000000    ~t_2            H-                                                              
          1769    1   53    0.000000    ~g              b                                                               

   1000006    307    ~t_1            ~t_1bar             2    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1770    1   53    0.000000    ~gravitino      t                                                               
          1771    1   53    0.000000    ~chi_1+         b                                                               
          1772    1   53    0.000000    ~chi_2+         b                                                               
          1773    1   53    0.000000    ~chi_10         t                                                               
          1774    1   53    0.000000    ~chi_20         t                                                               
          1775    1   53    0.000000    ~chi_30         t                                                               
          1776    1   53    0.000000    ~chi_40         t                                                               
          1777    1   53    0.000000    ~b_1            W+                                                              
          1778    1   53    0.000000    ~b_2            W+                                                              
          1779    1   53    0.000000    ~b_1            H+                                                              
          1780    1   53    0.000000    ~b_2            H+                                                              
          1781    1   53    0.000000    ~g              t                                                               
          1782    1   53    0.000000    ~chi_10         c                                                               
          1783   -1   53    0.000000    ~nu_tauL        tau+            b                                               
          1784   -1   53    0.000000    ~tau_1+         nu_tau          b                                               

   1000011    308    ~e_L-           ~e_L+              -3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1785    1   53    0.000000    ~gravitino      e-                                                              
          1786    1   53    0.000000    ~chi_1-         nu_e                                                            
          1787    1   53    0.000000    ~chi_2-         nu_e                                                            
          1788    1   53    0.000000    ~chi_10         e-                                                              
          1789    1   53    0.000000    ~chi_20         e-                                                              
          1790    1   53    0.000000    ~chi_30         e-                                                              
          1791    1   53    0.000000    ~chi_40         e-                                                              
          1792    1   53    0.000000    ~nu_eL          W-                                                              
          1793    1   53    0.000000    ~nu_eR          W-                                                              
          1794    1   53    0.000000    ~nu_eL          H-                                                              
          1795    1   53    0.000000    ~nu_eR          H-                                                              

   1000012    309    ~nu_eL          ~nu_eLbar           0    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1796    1   53    0.000000    ~gravitino      nu_e                                                            
          1797    1   53    0.000000    ~chi_1+         e-                                                              
          1798    1   53    0.000000    ~chi_2+         e-                                                              
          1799    1   53    0.000000    ~chi_10         nu_e                                                            
          1800    1   53    0.000000    ~chi_20         nu_e                                                            
          1801    1   53    0.000000    ~chi_30         nu_e                                                            
          1802    1   53    0.000000    ~chi_40         nu_e                                                            
          1803    1   53    0.000000    ~e_L-           W+                                                              
          1804    1   53    0.000000    ~e_R-           W+                                                              
          1805    1   53    0.000000    ~e_L-           H+                                                              
          1806    1   53    0.000000    ~e_R-           H+                                                              

   1000013    310    ~mu_L-          ~mu_L+             -3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1807    1   53    0.000000    ~gravitino      mu-                                                             
          1808    1   53    0.000000    ~chi_1-         nu_mu                                                           
          1809    1   53    0.000000    ~chi_2-         nu_mu                                                           
          1810    1   53    0.000000    ~chi_10         mu-                                                             
          1811    1   53    0.000000    ~chi_20         mu-                                                             
          1812    1   53    0.000000    ~chi_30         mu-                                                             
          1813    1   53    0.000000    ~chi_40         mu-                                                             
          1814    1   53    0.000000    ~nu_muL         W-                                                              
          1815    1   53    0.000000    ~nu_muR         W-                                                              
          1816    1   53    0.000000    ~nu_muL         H-                                                              
          1817    1   53    0.000000    ~nu_muR         H-                                                              

   1000014    311    ~nu_muL         ~nu_muLbar          0    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1818    1   53    0.000000    ~gravitino      nu_mu                                                           
          1819    1   53    0.000000    ~chi_1+         mu-                                                             
          1820    1   53    0.000000    ~chi_2+         mu-                                                             
          1821    1   53    0.000000    ~chi_10         nu_mu                                                           
          1822    1   53    0.000000    ~chi_20         nu_mu                                                           
          1823    1   53    0.000000    ~chi_30         nu_mu                                                           
          1824    1   53    0.000000    ~chi_40         nu_mu                                                           
          1825    1   53    0.000000    ~mu_L-          W+                                                              
          1826    1   53    0.000000    ~mu_R-          W+                                                              
          1827    1   53    0.000000    ~mu_L-          H+                                                              
          1828    1   53    0.000000    ~mu_R-          H+                                                              

   1000015    312    ~tau_1-         ~tau_1+            -3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1829    1   53    0.000000    ~gravitino      tau-                                                            
          1830    1   53    0.000000    ~chi_1-         nu_tau                                                          
          1831    1   53    0.000000    ~chi_2-         nu_tau                                                          
          1832    1   53    0.000000    ~chi_10         tau-                                                            
          1833    1   53    0.000000    ~chi_20         tau-                                                            
          1834    1   53    0.000000    ~chi_30         tau-                                                            
          1835    1   53    0.000000    ~chi_40         tau-                                                            
          1836    1   53    0.000000    ~nu_tauL        W-                                                              
          1837    1   53    0.000000    ~nu_tauR        W-                                                              
          1838    1   53    0.000000    ~nu_tauL        H-                                                              
          1839    1   53    0.000000    ~nu_tauR        H-                                                              

   1000016    313    ~nu_tauL        ~nu_tauLbar         0    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          1840    1   53    0.000000    ~gravitino      nu_tau                                                          
          1841    1   53    0.000000    ~chi_1+         tau-                                                            
          1842    1   53    0.000000    ~chi_2+         tau-                                                            
          1843    1   53    0.000000    ~chi_10         nu_tau                                                          
          1844    1   53    0.000000    ~chi_20         nu_tau                                                          
          1845    1   53    0.000000    ~chi_30         nu_tau                                                          
          1846    1   53    0.000000    ~chi_40         nu_tau                                                          
          1847    1   53    0.000000    ~tau_1-         W+                                                              
          1848    1   53    0.000000    ~tau_2-         W+                                                              
          1849    1   53    0.000000    ~tau_1-         H+                                                              
          1850    1   53    0.000000    ~tau_2-         H+                                                              

   1000021    314    ~g                                  0    2    0    500.00000     1.00000    10.00000   0.00000E+00    1
          1851    1   53    0.000000    ~gravitino      g                                                               
          1852    1   53    0.000000    ~d_L            dbar                                                            
          1853    1   53    0.000000    ~d_Lbar         d                                                               
          1854    1   53    0.000000    ~d_R            dbar                                                            
          1855    1   53    0.000000    ~d_Rbar         d                                                               
          1856    1   53    0.000000    ~u_L            ubar                                                            
          1857    1   53    0.000000    ~u_Lbar         u                                                               
          1858    1   53    0.000000    ~u_R            ubar                                                            
          1859    1   53    0.000000    ~u_Rbar         u                                                               
          1860    1   53    0.000000    ~s_L            sbar                                                            
          1861    1   53    0.000000    ~s_Lbar         s                                                               
          1862    1   53    0.000000    ~s_R            sbar                                                            
          1863    1   53    0.000000    ~s_Rbar         s                                                               
          1864    1   53    0.000000    ~c_L            cbar                                                            
          1865    1   53    0.000000    ~c_Lbar         c                                                               
          1866    1   53    0.000000    ~c_R            cbar                                                            
          1867    1   53    0.000000    ~c_Rbar         c                                                               
          1868    1   53    0.000000    ~b_1            bbar                                                            
          1869    1   53    0.000000    ~b_1bar         b                                                               
          1870    1   53    0.000000    ~b_2            bbar                                                            
          1871    1   53    0.000000    ~b_2bar         b                                                               
          1872    1   53    0.000000    ~t_1            tbar                                                            
          1873    1   53    0.000000    ~t_1bar         t                                                               
          1874    1   53    0.000000    ~t_2            tbar                                                            
          1875    1   53    0.000000    ~t_2bar         t                                                               
          1876    1   53    0.000000    ~chi_10         d               dbar                                            
          1877    1   53    0.000000    ~chi_10         s               sbar                                            
          1878    1   53    0.000000    ~chi_10         b               bbar                                            
          1879    1   53    0.000000    ~chi_10         u               ubar                                            
          1880    1   53    0.000000    ~chi_10         c               cbar                                            
          1881    1   53    0.000000    ~chi_10         t               tbar                                            
          1882    1   53    0.000000    ~chi_20         d               dbar                                            
          1883    1   53    0.000000    ~chi_20         s               sbar                                            
          1884    1   53    0.000000    ~chi_20         b               bbar                                            
          1885    1   53    0.000000    ~chi_20         u               ubar                                            
          1886    1   53    0.000000    ~chi_20         c               cbar                                            
          1887    1   53    0.000000    ~chi_20         t               tbar                                            
          1888    1   53    0.000000    ~chi_30         d               dbar                                            
          1889    1   53    0.000000    ~chi_30         s               sbar                                            
          1890    1   53    0.000000    ~chi_30         b               bbar                                            
          1891    1   53    0.000000    ~chi_30         u               ubar                                            
          1892    1   53    0.000000    ~chi_30         c               cbar                                            
          1893    1   53    0.000000    ~chi_30         t               tbar                                            
          1894    1   53    0.000000    ~chi_40         d               dbar                                            
          1895    1   53    0.000000    ~chi_40         s               sbar                                            
          1896    1   53    0.000000    ~chi_40         b               bbar                                            
          1897    1   53    0.000000    ~chi_40         u               ubar                                            
          1898    1   53    0.000000    ~chi_40         c               cbar                                            
          1899    1   53    0.000000    ~chi_40         t               tbar                                            
          1900    1   53    0.000000    ~chi_1+         d               ubar                                            
          1901    1   53    0.000000    ~chi_1-         dbar            u                                               
          1902    1   53    0.000000    ~chi_1+         s               cbar                                            
          1903    1   53    0.000000    ~chi_1-         sbar            c                                               
          1904    1   53    0.000000    ~chi_1+         b               tbar                                            
          1905    1   53    0.000000    ~chi_1-         bbar            t                                               
          1906    1   53    0.000000    ~chi_2+         d               ubar                                            
          1907    1   53    0.000000    ~chi_2-         dbar            u                                               
          1908    1   53    0.000000    ~chi_2+         s               cbar                                            
          1909    1   53    0.000000    ~chi_2-         sbar            c                                               
          1910    1   53    0.000000    ~chi_2+         b               tbar                                            
          1911    1   53    0.000000    ~chi_2-         bbar            t                                               

   1000022    315    ~chi_10                             0    0    0    500.00000     1.00000    10.00000   0.00000E+00    1
          1912    1   53    0.000000    ~gravitino      gamma                                                           
          1913    1   53    0.000000    ~gravitino      Z0                                                              
          1914    1   53    0.000000    ~gravitino      h0                                                              
          1915    1   53    0.000000    ~gravitino      H0                                                              
          1916    1   53    0.000000    ~gravitino      A0                                                              
          1917   -1   53    0.000000    c               dbar            e-                                              
          1918   -1   53    0.000000    d               sbar            nu_e                                            

   1000023    316    ~chi_20                             0    0    0    500.00000     1.00000    10.00000   0.00000E+00    1
          1919    1   53    0.000000    ~gravitino      gamma                                                           
          1920    1   53    0.000000    ~gravitino      Z0                                                              
          1921    1   53    0.000000    ~gravitino      h0                                                              
          1922    1   53    0.000000    ~gravitino      H0                                                              
          1923    1   53    0.000000    ~gravitino      A0                                                              
          1924    1   53    0.000000    ~chi_10         gamma                                                           
          1925    1   53    0.000000    ~chi_10         Z0                                                              
          1926    1   53    0.000000    ~chi_10         e-              e+                                              
          1927    1   53    0.000000    ~chi_10         mu-             mu+                                             
          1928    1   53    0.000000    ~chi_10         tau-            tau+                                            
          1929    1   53    0.000000    ~chi_10         nu_e            nu_ebar                                         
          1930    1   53    0.000000    ~chi_10         nu_mu           nu_mubar                                        
          1931    1   53    0.000000    ~chi_10         nu_tau          nu_taubar                                       
          1932    1   53    0.000000    ~chi_10         d               dbar                                            
          1933    1   53    0.000000    ~chi_10         s               sbar                                            
          1934    1   53    0.000000    ~chi_10         b               bbar                                            
          1935    1   53    0.000000    ~chi_10         u               ubar                                            
          1936    1   53    0.000000    ~chi_10         c               cbar                                            
          1937    1   53    0.000000    ~chi_10         h0                                                              
          1938    1   53    0.000000    ~chi_10         H0                                                              
          1939    1   53    0.000000    ~chi_10         A0                                                              
          1940    1   53    0.000000    ~chi_1+         W-                                                              
          1941    1   53    0.000000    ~chi_1-         W+                                                              
          1942    1   53    0.000000    ~chi_1+         e-              nu_ebar                                         
          1943    1   53    0.000000    ~chi_1-         e+              nu_e                                            
          1944    1   53    0.000000    ~chi_1+         mu-             nu_mubar                                        
          1945    1   53    0.000000    ~chi_1-         mu+             nu_mu                                           
          1946    1   53    0.000000    ~chi_1+         tau-            nu_taubar                                       
          1947    1   53    0.000000    ~chi_1-         tau+            nu_tau                                          
          1948    1   53    0.000000    ~chi_1+         d               ubar                                            
          1949    1   53    0.000000    ~chi_1-         dbar            u                                               
          1950    1   53    0.000000    ~chi_1+         s               cbar                                            
          1951    1   53    0.000000    ~chi_1-         sbar            c                                               
          1952    1   53    0.000000    ~chi_2+         W-                                                              
          1953    1   53    0.000000    ~chi_2-         W+                                                              
          1954    1   53    0.000000    ~chi_2+         e-              nu_ebar                                         
          1955    1   53    0.000000    ~chi_2-         e+              nu_e                                            
          1956    1   53    0.000000    ~chi_2+         mu-             nu_mubar                                        
          1957    1   53    0.000000    ~chi_2-         mu+             nu_mu                                           
          1958    1   53    0.000000    ~chi_2+         tau-            nu_taubar                                       
          1959    1   53    0.000000    ~chi_2-         tau+            nu_tau                                          
          1960    1   53    0.000000    ~chi_2+         d               ubar                                            
          1961    1   53    0.000000    ~chi_2-         dbar            u                                               
          1962    1   53    0.000000    ~chi_2+         s               cbar                                            
          1963    1   53    0.000000    ~chi_2-         sbar            c                                               
          1964    1   53    0.000000    ~chi_1+         H-                                                              
          1965    1   53    0.000000    ~chi_1-         H+                                                              
          1966    1   53    0.000000    ~chi_2+         H-                                                              
          1967    1   53    0.000000    ~chi_2-         H+                                                              
          1968    1   53    0.000000    ~d_L            dbar                                                            
          1969    1   53    0.000000    ~d_Lbar         d                                                               
          1970    1   53    0.000000    ~d_R            dbar                                                            
          1971    1   53    0.000000    ~d_Rbar         d                                                               
          1972    1   53    0.000000    ~u_L            ubar                                                            
          1973    1   53    0.000000    ~u_Lbar         u                                                               
          1974    1   53    0.000000    ~u_R            ubar                                                            
          1975    1   53    0.000000    ~u_Rbar         u                                                               
          1976    1   53    0.000000    ~s_L            sbar                                                            
          1977    1   53    0.000000    ~s_Lbar         s                                                               
          1978    1   53    0.000000    ~s_R            sbar                                                            
          1979    1   53    0.000000    ~s_Rbar         s                                                               
          1980    1   53    0.000000    ~c_L            cbar                                                            
          1981    1   53    0.000000    ~c_Lbar         c                                                               
          1982    1   53    0.000000    ~c_R            cbar                                                            
          1983    1   53    0.000000    ~c_Rbar         c                                                               
          1984    1   53    0.000000    ~b_1            bbar                                                            
          1985    1   53    0.000000    ~b_1bar         b                                                               
          1986    1   53    0.000000    ~b_2            bbar                                                            
          1987    1   53    0.000000    ~b_2bar         b                                                               
          1988    1   53    0.000000    ~t_1            tbar                                                            
          1989    1   53    0.000000    ~t_1bar         t                                                               
          1990    1   53    0.000000    ~t_2            tbar                                                            
          1991    1   53    0.000000    ~t_2bar         t                                                               
          1992    1   53    0.000000    ~e_L-           e+                                                              
          1993    1   53    0.000000    ~e_L+           e-                                                              
          1994    1   53    0.000000    ~e_R-           e+                                                              
          1995    1   53    0.000000    ~e_R+           e-                                                              
          1996    1   53    0.000000    ~nu_eL          nu_ebar                                                         
          1997    1   53    0.000000    ~nu_eLbar       nu_e                                                            
          1998    1   53    0.000000    ~nu_eR          nu_ebar                                                         
          1999    1   53    0.000000    ~nu_eRbar       nu_e                                                            
          2000    1   53    0.000000    ~mu_L-          mu+                                                             
          2001    1   53    0.000000    ~mu_L+          mu-                                                             
          2002    1   53    0.000000    ~mu_R-          mu+                                                             
          2003    1   53    0.000000    ~mu_R+          mu-                                                             
          2004    1   53    0.000000    ~nu_muL         nu_mubar                                                        
          2005    1   53    0.000000    ~nu_muLbar      nu_mu                                                           
          2006    1   53    0.000000    ~nu_muR         nu_mubar                                                        
          2007    1   53    0.000000    ~nu_muRbar      nu_mu                                                           
          2008    1   53    0.000000    ~tau_1-         tau+                                                            
          2009    1   53    0.000000    ~tau_1+         tau-                                                            
          2010    1   53    0.000000    ~tau_2-         tau+                                                            
          2011    1   53    0.000000    ~tau_2+         tau-                                                            
          2012    1   53    0.000000    ~nu_tauL        nu_taubar                                                       
          2013    1   53    0.000000    ~nu_tauLbar     nu_tau                                                          
          2014    1   53    0.000000    ~nu_tauR        nu_taubar                                                       
          2015    1   53    0.000000    ~nu_tauRbar     nu_tau                                                          
          2016    1   53    0.000000    ~g              d               dbar                                            
          2017    1   53    0.000000    ~g              s               sbar                                            
          2018    1   53    0.000000    ~g              b               bbar                                            
          2019    1   53    0.000000    ~g              u               ubar                                            
          2020    1   53    0.000000    ~g              c               cbar                                            

   1000024    317    ~chi_1+         ~chi_1-             3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2021    1   53    0.000000    ~gravitino      W+                                                              
          2022    1   53    0.000000    ~gravitino      H+                                                              
          2023    1   53    0.000000    ~chi_10         W+                                                              
          2024    1   53    0.000000    ~chi_10         e+              nu_e                                            
          2025    1   53    0.000000    ~chi_10         mu+             nu_mu                                           
          2026    1   53    0.000000    ~chi_10         tau+            nu_tau                                          
          2027    1   53    0.000000    ~chi_10         dbar            u                                               
          2028    1   53    0.000000    ~chi_10         sbar            c                                               
          2029    1   53    0.000000    ~chi_20         W+                                                              
          2030    1   53    0.000000    ~chi_20         e+              nu_e                                            
          2031    1   53    0.000000    ~chi_20         mu+             nu_mu                                           
          2032    1   53    0.000000    ~chi_20         tau+            nu_tau                                          
          2033    1   53    0.000000    ~chi_20         dbar            u                                               
          2034    1   53    0.000000    ~chi_20         sbar            c                                               
          2035    1   53    0.000000    ~chi_30         W+                                                              
          2036    1   53    0.000000    ~chi_30         e+              nu_e                                            
          2037    1   53    0.000000    ~chi_30         mu+             nu_mu                                           
          2038    1   53    0.000000    ~chi_30         tau+            nu_tau                                          
          2039    1   53    0.000000    ~chi_30         dbar            u                                               
          2040    1   53    0.000000    ~chi_30         sbar            c                                               
          2041    1   53    0.000000    ~chi_40         W+                                                              
          2042    1   53    0.000000    ~chi_40         e+              nu_e                                            
          2043    1   53    0.000000    ~chi_40         mu+             nu_mu                                           
          2044    1   53    0.000000    ~chi_40         tau+            nu_tau                                          
          2045    1   53    0.000000    ~chi_40         dbar            u                                               
          2046    1   53    0.000000    ~chi_40         sbar            c                                               
          2047    1   53    0.000000    ~chi_10         H+                                                              
          2048    1   53    0.000000    ~chi_20         H+                                                              
          2049    1   53    0.000000    ~chi_30         H+                                                              
          2050    1   53    0.000000    ~chi_40         H+                                                              
          2051    1   53    0.000000    ~u_L            dbar                                                            
          2052    1   53    0.000000    ~u_R            dbar                                                            
          2053    1   53    0.000000    ~d_Lbar         u                                                               
          2054    1   53    0.000000    ~d_Rbar         u                                                               
          2055    1   53    0.000000    ~c_L            sbar                                                            
          2056    1   53    0.000000    ~c_R            sbar                                                            
          2057    1   53    0.000000    ~s_Lbar         c                                                               
          2058    1   53    0.000000    ~s_Rbar         c                                                               
          2059    1   53    0.000000    ~t_1            bbar                                                            
          2060    1   53    0.000000    ~t_2            bbar                                                            
          2061    1   53    0.000000    ~b_1bar         t                                                               
          2062    1   53    0.000000    ~b_2bar         t                                                               
          2063    1   53    0.000000    ~nu_eL          e+                                                              
          2064    1   53    0.000000    ~nu_eR          e+                                                              
          2065    1   53    0.000000    ~e_L+           nu_e                                                            
          2066    1   53    0.000000    ~e_R+           nu_e                                                            
          2067    1   53    0.000000    ~nu_muL         mu+                                                             
          2068    1   53    0.000000    ~nu_muR         mu+                                                             
          2069    1   53    0.000000    ~mu_L+          nu_mu                                                           
          2070    1   53    0.000000    ~mu_R+          nu_mu                                                           
          2071    1   53    0.000000    ~nu_tauL        tau+                                                            
          2072    1   53    0.000000    ~nu_tauR        tau+                                                            
          2073    1   53    0.000000    ~tau_1+         nu_tau                                                          
          2074    1   53    0.000000    ~tau_2+         nu_tau                                                          
          2075    1   53    0.000000    ~g              dbar            u                                               
          2076    1   53    0.000000    ~g              sbar            c                                               

   1000025    318    ~chi_30                             0    0    0    500.00000     1.00000    10.00000   0.00000E+00    1
          2077    1   53    0.000000    ~gravitino      gamma                                                           
          2078    1   53    0.000000    ~gravitino      Z0                                                              
          2079    1   53    0.000000    ~gravitino      h0                                                              
          2080    1   53    0.000000    ~gravitino      H0                                                              
          2081    1   53    0.000000    ~gravitino      A0                                                              
          2082    1   53    0.000000    ~chi_10         gamma                                                           
          2083    1   53    0.000000    ~chi_10         Z0                                                              
          2084    1   53    0.000000    ~chi_10         e-              e+                                              
          2085    1   53    0.000000    ~chi_10         mu-             mu+                                             
          2086    1   53    0.000000    ~chi_10         tau-            tau+                                            
          2087    1   53    0.000000    ~chi_10         nu_e            nu_ebar                                         
          2088    1   53    0.000000    ~chi_10         nu_mu           nu_mubar                                        
          2089    1   53    0.000000    ~chi_10         nu_tau          nu_taubar                                       
          2090    1   53    0.000000    ~chi_10         d               dbar                                            
          2091    1   53    0.000000    ~chi_10         s               sbar                                            
          2092    1   53    0.000000    ~chi_10         b               bbar                                            
          2093    1   53    0.000000    ~chi_10         u               ubar                                            
          2094    1   53    0.000000    ~chi_10         c               cbar                                            
          2095    1   53    0.000000    ~chi_10         h0                                                              
          2096    1   53    0.000000    ~chi_10         H0                                                              
          2097    1   53    0.000000    ~chi_10         A0                                                              
          2098    1   53    0.000000    ~chi_20         gamma                                                           
          2099    1   53    0.000000    ~chi_20         Z0                                                              
          2100    1   53    0.000000    ~chi_20         e-              e+                                              
          2101    1   53    0.000000    ~chi_20         mu-             mu+                                             
          2102    1   53    0.000000    ~chi_20         tau-            tau+                                            
          2103    1   53    0.000000    ~chi_20         nu_e            nu_ebar                                         
          2104    1   53    0.000000    ~chi_20         nu_mu           nu_mubar                                        
          2105    1   53    0.000000    ~chi_20         nu_tau          nu_taubar                                       
          2106    1   53    0.000000    ~chi_20         d               dbar                                            
          2107    1   53    0.000000    ~chi_20         s               sbar                                            
          2108    1   53    0.000000    ~chi_20         b               bbar                                            
          2109    1   53    0.000000    ~chi_20         u               ubar                                            
          2110    1   53    0.000000    ~chi_20         c               cbar                                            
          2111    1   53    0.000000    ~chi_20         h0                                                              
          2112    1   53    0.000000    ~chi_20         H0                                                              
          2113    1   53    0.000000    ~chi_20         A0                                                              
          2114    1   53    0.000000    ~chi_1+         W-                                                              
          2115    1   53    0.000000    ~chi_1-         W+                                                              
          2116    1   53    0.000000    ~chi_1+         e-              nu_ebar                                         
          2117    1   53    0.000000    ~chi_1-         e+              nu_e                                            
          2118    1   53    0.000000    ~chi_1+         mu-             nu_mubar                                        
          2119    1   53    0.000000    ~chi_1-         mu+             nu_mu                                           
          2120    1   53    0.000000    ~chi_1+         tau-            nu_taubar                                       
          2121    1   53    0.000000    ~chi_1-         tau+            nu_tau                                          
          2122    1   53    0.000000    ~chi_1+         d               ubar                                            
          2123    1   53    0.000000    ~chi_1-         dbar            u                                               
          2124    1   53    0.000000    ~chi_1+         s               cbar                                            
          2125    1   53    0.000000    ~chi_1-         sbar            c                                               
          2126    1   53    0.000000    ~chi_2+         W-                                                              
          2127    1   53    0.000000    ~chi_2-         W+                                                              
          2128    1   53    0.000000    ~chi_2+         e-              nu_ebar                                         
          2129    1   53    0.000000    ~chi_2-         e+              nu_e                                            
          2130    1   53    0.000000    ~chi_2+         mu-             nu_mubar                                        
          2131    1   53    0.000000    ~chi_2-         mu+             nu_mu                                           
          2132    1   53    0.000000    ~chi_2+         tau-            nu_taubar                                       
          2133    1   53    0.000000    ~chi_2-         tau+            nu_tau                                          
          2134    1   53    0.000000    ~chi_2+         d               ubar                                            
          2135    1   53    0.000000    ~chi_2-         dbar            u                                               
          2136    1   53    0.000000    ~chi_2+         s               cbar                                            
          2137    1   53    0.000000    ~chi_2-         sbar            c                                               
          2138    1   53    0.000000    ~chi_1+         H-                                                              
          2139    1   53    0.000000    ~chi_1-         H+                                                              
          2140    1   53    0.000000    ~chi_2+         H-                                                              
          2141    1   53    0.000000    ~chi_2-         H+                                                              
          2142    1   53    0.000000    ~d_L            dbar                                                            
          2143    1   53    0.000000    ~d_Lbar         d                                                               
          2144    1   53    0.000000    ~d_R            dbar                                                            
          2145    1   53    0.000000    ~d_Rbar         d                                                               
          2146    1   53    0.000000    ~u_L            ubar                                                            
          2147    1   53    0.000000    ~u_Lbar         u                                                               
          2148    1   53    0.000000    ~u_R            ubar                                                            
          2149    1   53    0.000000    ~u_Rbar         u                                                               
          2150    1   53    0.000000    ~s_L            sbar                                                            
          2151    1   53    0.000000    ~s_Lbar         s                                                               
          2152    1   53    0.000000    ~s_R            sbar                                                            
          2153    1   53    0.000000    ~s_Rbar         s                                                               
          2154    1   53    0.000000    ~c_L            cbar                                                            
          2155    1   53    0.000000    ~c_Lbar         c                                                               
          2156    1   53    0.000000    ~c_R            cbar                                                            
          2157    1   53    0.000000    ~c_Rbar         c                                                               
          2158    1   53    0.000000    ~b_1            bbar                                                            
          2159    1   53    0.000000    ~b_1bar         b                                                               
          2160    1   53    0.000000    ~b_2            bbar                                                            
          2161    1   53    0.000000    ~b_2bar         b                                                               
          2162    1   53    0.000000    ~t_1            tbar                                                            
          2163    1   53    0.000000    ~t_1bar         t                                                               
          2164    1   53    0.000000    ~t_2            tbar                                                            
          2165    1   53    0.000000    ~t_2bar         t                                                               
          2166    1   53    0.000000    ~e_L-           e+                                                              
          2167    1   53    0.000000    ~e_L+           e-                                                              
          2168    1   53    0.000000    ~e_R-           e+                                                              
          2169    1   53    0.000000    ~e_R+           e-                                                              
          2170    1   53    0.000000    ~nu_eL          nu_ebar                                                         
          2171    1   53    0.000000    ~nu_eLbar       nu_e                                                            
          2172    1   53    0.000000    ~nu_eR          nu_ebar                                                         
          2173    1   53    0.000000    ~nu_eRbar       nu_e                                                            
          2174    1   53    0.000000    ~mu_L-          mu+                                                             
          2175    1   53    0.000000    ~mu_L+          mu-                                                             
          2176    1   53    0.000000    ~mu_R-          mu+                                                             
          2177    1   53    0.000000    ~mu_R+          mu-                                                             
          2178    1   53    0.000000    ~nu_muL         nu_mubar                                                        
          2179    1   53    0.000000    ~nu_muLbar      nu_mu                                                           
          2180    1   53    0.000000    ~nu_muR         nu_mubar                                                        
          2181    1   53    0.000000    ~nu_muRbar      nu_mu                                                           
          2182    1   53    0.000000    ~tau_1-         tau+                                                            
          2183    1   53    0.000000    ~tau_1+         tau-                                                            
          2184    1   53    0.000000    ~tau_2-         tau+                                                            
          2185    1   53    0.000000    ~tau_2+         tau-                                                            
          2186    1   53    0.000000    ~nu_tauL        nu_taubar                                                       
          2187    1   53    0.000000    ~nu_tauLbar     nu_tau                                                          
          2188    1   53    0.000000    ~nu_tauR        nu_taubar                                                       
          2189    1   53    0.000000    ~nu_tauRbar     nu_tau                                                          
          2190    1   53    0.000000    ~g              d               dbar                                            
          2191    1   53    0.000000    ~g              s               sbar                                            
          2192    1   53    0.000000    ~g              b               bbar                                            
          2193    1   53    0.000000    ~g              u               ubar                                            
          2194    1   53    0.000000    ~g              c               cbar                                            

   1000035    319    ~chi_40                             0    0    0    500.00000     1.00000    10.00000   0.00000E+00    1
          2195    1   53    0.000000    ~gravitino      gamma                                                           
          2196    1   53    0.000000    ~gravitino      Z0                                                              
          2197    1   53    0.000000    ~gravitino      h0                                                              
          2198    1   53    0.000000    ~gravitino      H0                                                              
          2199    1   53    0.000000    ~gravitino      A0                                                              
          2200    1   53    0.000000    ~chi_10         gamma                                                           
          2201    1   53    0.000000    ~chi_10         Z0                                                              
          2202    1   53    0.000000    ~chi_10         e-              e+                                              
          2203    1   53    0.000000    ~chi_10         mu-             mu+                                             
          2204    1   53    0.000000    ~chi_10         tau-            tau+                                            
          2205    1   53    0.000000    ~chi_10         nu_e            nu_ebar                                         
          2206    1   53    0.000000    ~chi_10         nu_mu           nu_mubar                                        
          2207    1   53    0.000000    ~chi_10         nu_tau          nu_taubar                                       
          2208    1   53    0.000000    ~chi_10         d               dbar                                            
          2209    1   53    0.000000    ~chi_10         s               sbar                                            
          2210    1   53    0.000000    ~chi_10         b               bbar                                            
          2211    1   53    0.000000    ~chi_10         u               ubar                                            
          2212    1   53    0.000000    ~chi_10         c               cbar                                            
          2213    1   53    0.000000    ~chi_10         h0                                                              
          2214    1   53    0.000000    ~chi_10         H0                                                              
          2215    1   53    0.000000    ~chi_10         A0                                                              
          2216    1   53    0.000000    ~chi_20         gamma                                                           
          2217    1   53    0.000000    ~chi_20         Z0                                                              
          2218    1   53    0.000000    ~chi_20         e-              e+                                              
          2219    1   53    0.000000    ~chi_20         mu-             mu+                                             
          2220    1   53    0.000000    ~chi_20         tau-            tau+                                            
          2221    1   53    0.000000    ~chi_20         nu_e            nu_ebar                                         
          2222    1   53    0.000000    ~chi_20         nu_mu           nu_mubar                                        
          2223    1   53    0.000000    ~chi_20         nu_tau          nu_taubar                                       
          2224    1   53    0.000000    ~chi_20         d               dbar                                            
          2225    1   53    0.000000    ~chi_20         s               sbar                                            
          2226    1   53    0.000000    ~chi_20         b               bbar                                            
          2227    1   53    0.000000    ~chi_20         u               ubar                                            
          2228    1   53    0.000000    ~chi_20         c               cbar                                            
          2229    1   53    0.000000    ~chi_20         h0                                                              
          2230    1   53    0.000000    ~chi_20         H0                                                              
          2231    1   53    0.000000    ~chi_20         A0                                                              
          2232    1   53    0.000000    ~chi_30         gamma                                                           
          2233    1   53    0.000000    ~chi_30         Z0                                                              
          2234    1   53    0.000000    ~chi_30         e-              e+                                              
          2235    1   53    0.000000    ~chi_30         mu-             mu+                                             
          2236    1   53    0.000000    ~chi_30         tau-            tau+                                            
          2237    1   53    0.000000    ~chi_30         nu_e            nu_ebar                                         
          2238    1   53    0.000000    ~chi_30         nu_mu           nu_mubar                                        
          2239    1   53    0.000000    ~chi_30         nu_tau          nu_taubar                                       
          2240    1   53    0.000000    ~chi_30         d               dbar                                            
          2241    1   53    0.000000    ~chi_30         s               sbar                                            
          2242    1   53    0.000000    ~chi_30         b               bbar                                            
          2243    1   53    0.000000    ~chi_30         u               ubar                                            
          2244    1   53    0.000000    ~chi_30         c               cbar                                            
          2245    1   53    0.000000    ~chi_30         h0                                                              
          2246    1   53    0.000000    ~chi_30         H0                                                              
          2247    1   53    0.000000    ~chi_30         A0                                                              
          2248    1   53    0.000000    ~chi_1+         W-                                                              
          2249    1   53    0.000000    ~chi_1-         W+                                                              
          2250    1   53    0.000000    ~chi_1+         e-              nu_ebar                                         
          2251    1   53    0.000000    ~chi_1-         e+              nu_e                                            
          2252    1   53    0.000000    ~chi_1+         mu-             nu_mubar                                        
          2253    1   53    0.000000    ~chi_1-         mu+             nu_mu                                           
          2254    1   53    0.000000    ~chi_1+         tau-            nu_taubar                                       
          2255    1   53    0.000000    ~chi_1-         tau+            nu_tau                                          
          2256    1   53    0.000000    ~chi_1+         d               ubar                                            
          2257    1   53    0.000000    ~chi_1-         dbar            u                                               
          2258    1   53    0.000000    ~chi_1+         s               cbar                                            
          2259    1   53    0.000000    ~chi_1-         sbar            c                                               
          2260    1   53    0.000000    ~chi_2+         W-                                                              
          2261    1   53    0.000000    ~chi_2-         W+                                                              
          2262    1   53    0.000000    ~chi_2+         e-              nu_ebar                                         
          2263    1   53    0.000000    ~chi_2-         e+              nu_e                                            
          2264    1   53    0.000000    ~chi_2+         mu-             nu_mubar                                        
          2265    1   53    0.000000    ~chi_2-         mu+             nu_mu                                           
          2266    1   53    0.000000    ~chi_2+         tau-            nu_taubar                                       
          2267    1   53    0.000000    ~chi_2-         tau+            nu_tau                                          
          2268    1   53    0.000000    ~chi_2+         d               ubar                                            
          2269    1   53    0.000000    ~chi_2-         dbar            u                                               
          2270    1   53    0.000000    ~chi_2+         s               cbar                                            
          2271    1   53    0.000000    ~chi_2-         sbar            c                                               
          2272    1   53    0.000000    ~chi_1+         H-                                                              
          2273    1   53    0.000000    ~chi_1-         H+                                                              
          2274    1   53    0.000000    ~chi_2+         H-                                                              
          2275    1   53    0.000000    ~chi_2-         H+                                                              
          2276    1   53    0.000000    ~d_L            dbar                                                            
          2277    1   53    0.000000    ~d_Lbar         d                                                               
          2278    1   53    0.000000    ~d_R            dbar                                                            
          2279    1   53    0.000000    ~d_Rbar         d                                                               
          2280    1   53    0.000000    ~u_L            ubar                                                            
          2281    1   53    0.000000    ~u_Lbar         u                                                               
          2282    1   53    0.000000    ~u_R            ubar                                                            
          2283    1   53    0.000000    ~u_Rbar         u                                                               
          2284    1   53    0.000000    ~s_L            sbar                                                            
          2285    1   53    0.000000    ~s_Lbar         s                                                               
          2286    1   53    0.000000    ~s_R            sbar                                                            
          2287    1   53    0.000000    ~s_Rbar         s                                                               
          2288    1   53    0.000000    ~c_L            cbar                                                            
          2289    1   53    0.000000    ~c_Lbar         c                                                               
          2290    1   53    0.000000    ~c_R            cbar                                                            
          2291    1   53    0.000000    ~c_Rbar         c                                                               
          2292    1   53    0.000000    ~b_1            bbar                                                            
          2293    1   53    0.000000    ~b_1bar         b                                                               
          2294    1   53    0.000000    ~b_2            bbar                                                            
          2295    1   53    0.000000    ~b_2bar         b                                                               
          2296    1   53    0.000000    ~t_1            tbar                                                            
          2297    1   53    0.000000    ~t_1bar         t                                                               
          2298    1   53    0.000000    ~t_2            tbar                                                            
          2299    1   53    0.000000    ~t_2bar         t                                                               
          2300    1   53    0.000000    ~e_L-           e+                                                              
          2301    1   53    0.000000    ~e_L+           e-                                                              
          2302    1   53    0.000000    ~e_R-           e+                                                              
          2303    1   53    0.000000    ~e_R+           e-                                                              
          2304    1   53    0.000000    ~nu_eL          nu_ebar                                                         
          2305    1   53    0.000000    ~nu_eLbar       nu_e                                                            
          2306    1   53    0.000000    ~nu_eR          nu_ebar                                                         
          2307    1   53    0.000000    ~nu_eRbar       nu_e                                                            
          2308    1   53    0.000000    ~mu_L-          mu+                                                             
          2309    1   53    0.000000    ~mu_L+          mu-                                                             
          2310    1   53    0.000000    ~mu_R-          mu+                                                             
          2311    1   53    0.000000    ~mu_R+          mu-                                                             
          2312    1   53    0.000000    ~nu_muL         nu_mubar                                                        
          2313    1   53    0.000000    ~nu_muLbar      nu_mu                                                           
          2314    1   53    0.000000    ~nu_muR         nu_mubar                                                        
          2315    1   53    0.000000    ~nu_muRbar      nu_mu                                                           
          2316    1   53    0.000000    ~tau_1-         tau+                                                            
          2317    1   53    0.000000    ~tau_1+         tau-                                                            
          2318    1   53    0.000000    ~tau_2-         tau+                                                            
          2319    1   53    0.000000    ~tau_2+         tau-                                                            
          2320    1   53    0.000000    ~nu_tauL        nu_taubar                                                       
          2321    1   53    0.000000    ~nu_tauLbar     nu_tau                                                          
          2322    1   53    0.000000    ~nu_tauR        nu_taubar                                                       
          2323    1   53    0.000000    ~nu_tauRbar     nu_tau                                                          
          2324    1   53    0.000000    ~g              d               dbar                                            
          2325    1   53    0.000000    ~g              s               sbar                                            
          2326    1   53    0.000000    ~g              b               bbar                                            
          2327    1   53    0.000000    ~g              u               ubar                                            
          2328    1   53    0.000000    ~g              c               cbar                                            

   1000037    320    ~chi_2+         ~chi_2-             3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2329    1   53    0.000000    ~gravitino      W+                                                              
          2330    1   53    0.000000    ~gravitino      H+                                                              
          2331    1   53    0.000000    ~chi_1+         Z0                                                              
          2332    1   53    0.000000    ~chi_1+         e-              e+                                              
          2333    1   53    0.000000    ~chi_1+         mu-             mu+                                             
          2334    1   53    0.000000    ~chi_1+         tau-            tau+                                            
          2335    1   53    0.000000    ~chi_1+         nu_e            nu_ebar                                         
          2336    1   53    0.000000    ~chi_1+         nu_mu           nu_mubar                                        
          2337    1   53    0.000000    ~chi_1+         nu_tau          nu_taubar                                       
          2338    1   53    0.000000    ~chi_1+         d               dbar                                            
          2339    1   53    0.000000    ~chi_1+         s               sbar                                            
          2340    1   53    0.000000    ~chi_1+         b               bbar                                            
          2341    1   53    0.000000    ~chi_1+         u               ubar                                            
          2342    1   53    0.000000    ~chi_1+         c               cbar                                            
          2343    1   53    0.000000    ~chi_1+         h0                                                              
          2344    1   53    0.000000    ~chi_1+         H0                                                              
          2345    1   53    0.000000    ~chi_1+         A0                                                              
          2346    1   53    0.000000    ~chi_10         W+                                                              
          2347    1   53    0.000000    ~chi_10         e+              nu_e                                            
          2348    1   53    0.000000    ~chi_10         mu+             nu_mu                                           
          2349    1   53    0.000000    ~chi_10         tau+            nu_tau                                          
          2350    1   53    0.000000    ~chi_10         dbar            u                                               
          2351    1   53    0.000000    ~chi_10         sbar            c                                               
          2352    1   53    0.000000    ~chi_20         W+                                                              
          2353    1   53    0.000000    ~chi_20         e+              nu_e                                            
          2354    1   53    0.000000    ~chi_20         mu+             nu_mu                                           
          2355    1   53    0.000000    ~chi_20         tau+            nu_tau                                          
          2356    1   53    0.000000    ~chi_20         dbar            u                                               
          2357    1   53    0.000000    ~chi_20         sbar            c                                               
          2358    1   53    0.000000    ~chi_30         W+                                                              
          2359    1   53    0.000000    ~chi_30         e+              nu_e                                            
          2360    1   53    0.000000    ~chi_30         mu+             nu_mu                                           
          2361    1   53    0.000000    ~chi_30         tau+            nu_tau                                          
          2362    1   53    0.000000    ~chi_30         dbar            u                                               
          2363    1   53    0.000000    ~chi_30         sbar            c                                               
          2364    1   53    0.000000    ~chi_40         W+                                                              
          2365    1   53    0.000000    ~chi_40         e+              nu_e                                            
          2366    1   53    0.000000    ~chi_40         mu+             nu_mu                                           
          2367    1   53    0.000000    ~chi_40         tau+            nu_tau                                          
          2368    1   53    0.000000    ~chi_40         dbar            u                                               
          2369    1   53    0.000000    ~chi_40         sbar            c                                               
          2370    1   53    0.000000    ~chi_10         H+                                                              
          2371    1   53    0.000000    ~chi_20         H+                                                              
          2372    1   53    0.000000    ~chi_30         H+                                                              
          2373    1   53    0.000000    ~chi_40         H+                                                              
          2374    1   53    0.000000    ~u_L            dbar                                                            
          2375    1   53    0.000000    ~u_R            dbar                                                            
          2376    1   53    0.000000    ~d_Lbar         u                                                               
          2377    1   53    0.000000    ~d_Rbar         u                                                               
          2378    1   53    0.000000    ~c_L            sbar                                                            
          2379    1   53    0.000000    ~c_R            sbar                                                            
          2380    1   53    0.000000    ~s_Lbar         c                                                               
          2381    1   53    0.000000    ~s_Rbar         c                                                               
          2382    1   53    0.000000    ~t_1            bbar                                                            
          2383    1   53    0.000000    ~t_2            bbar                                                            
          2384    1   53    0.000000    ~b_1bar         t                                                               
          2385    1   53    0.000000    ~b_2bar         t                                                               
          2386    1   53    0.000000    ~nu_eL          e+                                                              
          2387    1   53    0.000000    ~nu_eR          e+                                                              
          2388    1   53    0.000000    ~e_L+           nu_e                                                            
          2389    1   53    0.000000    ~e_R+           nu_e                                                            
          2390    1   53    0.000000    ~nu_muL         mu+                                                             
          2391    1   53    0.000000    ~nu_muR         mu+                                                             
          2392    1   53    0.000000    ~mu_L+          nu_mu                                                           
          2393    1   53    0.000000    ~mu_R+          nu_mu                                                           
          2394    1   53    0.000000    ~nu_tauL        tau+                                                            
          2395    1   53    0.000000    ~nu_tauR        tau+                                                            
          2396    1   53    0.000000    ~tau_1+         nu_tau                                                          
          2397    1   53    0.000000    ~tau_2+         nu_tau                                                          
          2398    1   53    0.000000    ~g              dbar            u                                               
          2399    1   53    0.000000    ~g              sbar            c                                               

   1000039    321    ~gravitino                          0    0    0    500.00000     0.00000     0.00001   0.00000E+00    0

   2000001    322    ~d_R            ~d_Rbar            -1    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2400    1   53    0.000000    ~gravitino      d                                                               
          2401    1   53    0.000000    ~chi_1-         u                                                               
          2402    1   53    0.000000    ~chi_2-         u                                                               
          2403    1   53    0.000000    ~chi_10         d                                                               
          2404    1   53    0.000000    ~chi_20         d                                                               
          2405    1   53    0.000000    ~chi_30         d                                                               
          2406    1   53    0.000000    ~chi_40         d                                                               
          2407    1   53    0.000000    ~d_L            Z0                                                              
          2408    1   53    0.000000    ~d_L            h0                                                              
          2409    1   53    0.000000    ~d_L            H0                                                              
          2410    1   53    0.000000    ~d_L            A0                                                              
          2411    1   53    0.000000    ~u_L            W-                                                              
          2412    1   53    0.000000    ~u_R            W-                                                              
          2413    1   53    0.000000    ~u_L            H-                                                              
          2414    1   53    0.000000    ~u_R            H-                                                              
          2415    1   53    0.000000    ~g              d                                                               

   2000002    323    ~u_R            ~u_Rbar             2    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2416    1   53    0.000000    ~gravitino      u                                                               
          2417    1   53    0.000000    ~chi_1+         d                                                               
          2418    1   53    0.000000    ~chi_2+         d                                                               
          2419    1   53    0.000000    ~chi_10         u                                                               
          2420    1   53    0.000000    ~chi_20         u                                                               
          2421    1   53    0.000000    ~chi_30         u                                                               
          2422    1   53    0.000000    ~chi_40         u                                                               
          2423    1   53    0.000000    ~u_L            Z0                                                              
          2424    1   53    0.000000    ~u_L            h0                                                              
          2425    1   53    0.000000    ~u_L            H0                                                              
          2426    1   53    0.000000    ~u_L            A0                                                              
          2427    1   53    0.000000    ~d_L            W+                                                              
          2428    1   53    0.000000    ~d_R            W+                                                              
          2429    1   53    0.000000    ~d_L            H+                                                              
          2430    1   53    0.000000    ~d_R            H+                                                              
          2431    1   53    0.000000    ~g              u                                                               

   2000003    324    ~s_R            ~s_Rbar            -1    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2432    1   53    0.000000    ~gravitino      s                                                               
          2433    1   53    0.000000    ~chi_1-         c                                                               
          2434    1   53    0.000000    ~chi_2-         c                                                               
          2435    1   53    0.000000    ~chi_10         s                                                               
          2436    1   53    0.000000    ~chi_20         s                                                               
          2437    1   53    0.000000    ~chi_30         s                                                               
          2438    1   53    0.000000    ~chi_40         s                                                               
          2439    1   53    0.000000    ~s_L            Z0                                                              
          2440    1   53    0.000000    ~s_L            h0                                                              
          2441    1   53    0.000000    ~s_L            H0                                                              
          2442    1   53    0.000000    ~s_L            A0                                                              
          2443    1   53    0.000000    ~c_L            W-                                                              
          2444    1   53    0.000000    ~c_R            W-                                                              
          2445    1   53    0.000000    ~c_L            H-                                                              
          2446    1   53    0.000000    ~c_R            H-                                                              
          2447    1   53    0.000000    ~g              s                                                               

   2000004    325    ~c_R            ~c_Rbar             2    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2448    1   53    0.000000    ~gravitino      c                                                               
          2449    1   53    0.000000    ~chi_1+         s                                                               
          2450    1   53    0.000000    ~chi_2+         s                                                               
          2451    1   53    0.000000    ~chi_10         c                                                               
          2452    1   53    0.000000    ~chi_20         c                                                               
          2453    1   53    0.000000    ~chi_30         c                                                               
          2454    1   53    0.000000    ~chi_40         c                                                               
          2455    1   53    0.000000    ~c_L            Z0                                                              
          2456    1   53    0.000000    ~c_L            h0                                                              
          2457    1   53    0.000000    ~c_L            H0                                                              
          2458    1   53    0.000000    ~c_L            A0                                                              
          2459    1   53    0.000000    ~s_L            W+                                                              
          2460    1   53    0.000000    ~s_R            W+                                                              
          2461    1   53    0.000000    ~s_L            H+                                                              
          2462    1   53    0.000000    ~s_R            H+                                                              
          2463    1   53    0.000000    ~g              c                                                               

   2000005    326    ~b_2            ~b_2bar            -1    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2464    1   53    0.000000    ~gravitino      b                                                               
          2465    1   53    0.000000    ~chi_1-         t                                                               
          2466    1   53    0.000000    ~chi_2-         t                                                               
          2467    1   53    0.000000    ~chi_10         b                                                               
          2468    1   53    0.000000    ~chi_20         b                                                               
          2469    1   53    0.000000    ~chi_30         b                                                               
          2470    1   53    0.000000    ~chi_40         b                                                               
          2471    1   53    0.000000    ~b_1            Z0                                                              
          2472    1   53    0.000000    ~b_1            h0                                                              
          2473    1   53    0.000000    ~b_1            H0                                                              
          2474    1   53    0.000000    ~b_1            A0                                                              
          2475    1   53    0.000000    ~t_1            W-                                                              
          2476    1   53    0.000000    ~t_2            W-                                                              
          2477    1   53    0.000000    ~t_1            H-                                                              
          2478    1   53    0.000000    ~t_2            H-                                                              
          2479    1   53    0.000000    ~g              b                                                               

   2000006    327    ~t_2            ~t_2bar             2    1    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2480    1   53    0.000000    ~gravitino      t                                                               
          2481    1   53    0.000000    ~chi_1+         b                                                               
          2482    1   53    0.000000    ~chi_2+         b                                                               
          2483    1   53    0.000000    ~chi_10         t                                                               
          2484    1   53    0.000000    ~chi_20         t                                                               
          2485    1   53    0.000000    ~chi_30         t                                                               
          2486    1   53    0.000000    ~chi_40         t                                                               
          2487    1   53    0.000000    ~t_1            Z0                                                              
          2488    1   53    0.000000    ~t_1            h0                                                              
          2489    1   53    0.000000    ~t_1            H0                                                              
          2490    1   53    0.000000    ~t_1            A0                                                              
          2491    1   53    0.000000    ~b_1            W+                                                              
          2492    1   53    0.000000    ~b_2            W+                                                              
          2493    1   53    0.000000    ~b_1            H+                                                              
          2494    1   53    0.000000    ~b_2            H+                                                              
          2495    1   53    0.000000    ~g              t                                                               

   2000011    328    ~e_R-           ~e_R+              -3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2496    1   53    0.000000    ~gravitino      e-                                                              
          2497    1   53    0.000000    ~chi_1-         nu_e                                                            
          2498    1   53    0.000000    ~chi_2-         nu_e                                                            
          2499    1   53    0.000000    ~chi_10         e-                                                              
          2500    1   53    0.000000    ~chi_20         e-                                                              
          2501    1   53    0.000000    ~chi_30         e-                                                              
          2502    1   53    0.000000    ~chi_40         e-                                                              
          2503    1   53    0.000000    ~e_L-           Z0                                                              
          2504    1   53    0.000000    ~e_L-           h0                                                              
          2505    1   53    0.000000    ~e_L-           H0                                                              
          2506    1   53    0.000000    ~e_L-           A0                                                              
          2507    1   53    0.000000    ~nu_eL          W-                                                              
          2508    1   53    0.000000    ~nu_eR          W-                                                              
          2509    1   53    0.000000    ~nu_eL          H-                                                              
          2510    1   53    0.000000    ~nu_eR          H-                                                              

   2000012    329    ~nu_eR          ~nu_eRbar           0    0    1    500.00000     0.00000     0.00001   0.00000E+00    0

   2000013    330    ~mu_R-          ~mu_R+             -3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2511    1   53    0.000000    ~gravitino      mu-                                                             
          2512    1   53    0.000000    ~chi_1-         nu_mu                                                           
          2513    1   53    0.000000    ~chi_2-         nu_mu                                                           
          2514    1   53    0.000000    ~chi_10         mu-                                                             
          2515    1   53    0.000000    ~chi_20         mu-                                                             
          2516    1   53    0.000000    ~chi_30         mu-                                                             
          2517    1   53    0.000000    ~chi_40         mu-                                                             
          2518    1   53    0.000000    ~mu_L-          Z0                                                              
          2519    1   53    0.000000    ~mu_L-          h0                                                              
          2520    1   53    0.000000    ~mu_L-          H0                                                              
          2521    1   53    0.000000    ~mu_L-          A0                                                              
          2522    1   53    0.000000    ~nu_muL         W-                                                              
          2523    1   53    0.000000    ~nu_muR         W-                                                              
          2524    1   53    0.000000    ~nu_muL         H-                                                              
          2525    1   53    0.000000    ~nu_muR         H-                                                              

   2000014    331    ~nu_muR         ~nu_muRbar          0    0    1    500.00000     0.00000     0.00001   0.00000E+00    0

   2000015    332    ~tau_2-         ~tau_2+            -3    0    1    500.00000     1.00000    10.00000   0.00000E+00    1
          2526    1   53    0.000000    ~gravitino      tau-                                                            
          2527    1   53    0.000000    ~chi_1-         nu_tau                                                          
          2528    1   53    0.000000    ~chi_2-         nu_tau                                                          
          2529    1   53    0.000000    ~chi_10         tau-                                                            
          2530    1   53    0.000000    ~chi_20         tau-                                                            
          2531    1   53    0.000000    ~chi_30         tau-                                                            
          2532    1   53    0.000000    ~chi_40         tau-                                                            
          2533    1   53    0.000000    ~tau_1-         Z0                                                              
          2534    1   53    0.000000    ~tau_1-         h0                                                              
          2535    1   53    0.000000    ~tau_1-         H0                                                              
          2536    1   53    0.000000    ~tau_1-         A0                                                              
          2537    1   53    0.000000    ~nu_tauL        W-                                                              
          2538    1   53    0.000000    ~nu_tauR        W-                                                              
          2539    1   53    0.000000    ~nu_tauL        H-                                                              
          2540    1   53    0.000000    ~nu_tauR        H-                                                              

   2000016    333    ~nu_tauR        ~nu_tauRbar         0    0    1    500.00000     0.00000     0.00001   0.00000E+00    0

   4000001    334    d*              d*bar              -1    1    1    400.00000     2.60608    26.06076   0.00000E+00    1
          2541    1   53    0.851667    g               d                                                               
          2542    1    0    0.005385    gamma           d                                                               
          2543    1    0    0.044810    Z0              d                                                               
          2544    1    0    0.098138    W-              u                                                               

   4000002    335    u*              u*bar               2    1    1    400.00000     2.60935    26.09355   0.00000E+00    1
          2545    1    0    0.850597    g               u                                                               
          2546    1    0    0.021514    gamma           u                                                               
          2547    1    0    0.029875    Z0              u                                                               
          2548    1    0    0.098014    W+              d                                                               

   4000011    336    e*-             e*bar+             -3    0    1    400.00000     0.42901     4.29011   0.00000E+00    1
          2549    1    0    0.294414    gamma           e-                                                              
          2550    1    0    0.109437    Z0              e-                                                              
          2551    1    0    0.596149    W-              nu_e                                                            

   4000012    337    nu*_e0          nu*_ebar0           0    0    1    400.00000     0.41917     4.19173   0.00000E+00    1
          2552    1    0    0.389861    Z0              nu_e                                                            
          2553    1    0    0.610139    W+              e-                                                              

*/
