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
Revision 1.25  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.20.6.1  2002/06/10 14:57:41  hristov
Merged with v3-08-02

Revision 1.24  2002/05/22 13:22:53  morsch
Process kPyMbNonDiffr added.

Revision 1.23  2002/05/06 07:17:29  morsch
Pyr gives random number r in interval 0 < r < 1.

Revision 1.22  2002/04/26 10:28:48  morsch
Option kPyBeautyPbMNR added (N. Carrer).

Revision 1.21  2002/03/25 14:46:16  morsch
Case  kPyD0PbMNR added (N. Carrer).

Revision 1.20  2002/03/03 13:48:50  morsch
Option  kPyCharmPbMNR added. Produce charm pairs in agreement with MNR
NLO calculations (Nicola Carrer).

Revision 1.19  2002/02/20 08:52:20  morsch
Correct documentation of SetNuclei method.

Revision 1.18  2002/02/07 10:43:06  morsch
Tuned pp-min.bias settings (M.Monteno, R.Ugoccioni and N.Carrer)

Revision 1.17  2001/12/19 15:40:43  morsch
For kPyJets enforce simple jet topology, i.e no initial or final state
gluon radiation and no primordial pT.

Revision 1.16  2001/10/12 11:13:59  morsch
Missing break statements added (thanks to Nicola Carrer)

Revision 1.15  2001/03/27 10:54:50  morsch
Add ResetDecayTable() and SsetDecayTable() methods.

Revision 1.14  2001/03/09 13:03:40  morsch
Process_t and Struc_Func_t moved to AliPythia.h

Revision 1.13  2000/12/18 08:55:35  morsch
Make AliPythia dependent generartors work with new scheme of random number generation

Revision 1.12  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.11  2000/10/20 06:30:06  fca
Use version 0 to avoid streamer generation

Revision 1.10  2000/10/06 14:18:44  morsch
Upper cut of prim. pT distribution set to 5. GeV

Revision 1.9  2000/09/18 10:41:35  morsch
Add possibility to use nuclear structure functions from PDF library V8.

Revision 1.8  2000/09/06 14:26:24  morsch
Decayer functionality of AliPythia has been moved to AliDecayerPythia.
Class is now a singleton.

Revision 1.7  2000/06/09 20:34:50  morsch
All coding rule violations except RS3 corrected

Revision 1.6  1999/11/09 07:38:48  fca
Changes for compatibility with version 2.23 of ROOT

Revision 1.5  1999/11/03 17:43:20  fca
New version from G.Martinez & A.Morsch

Revision 1.4  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/


#include "AliPythia.h"

ClassImp(AliPythia)

//_____________________________________________________________________________

AliPythia* AliPythia::fgAliPythia=NULL;

AliPythia::AliPythia()
{
// Default Constructor
//
//  Set random number
    if (!sRandom) sRandom=fRandom;

}

void AliPythia::ProcInit(Process_t process, Float_t energy, StrucFunc_t strucfunc)
{
// Initialise the process to generate 
    fProcess = process;
    fEcms = energy;
    fStrucFunc = strucfunc;
//  don't decay p0
    SetMDCY(Pycomp(111),1,0);
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
    case kPyCharm:
	SetMSEL(4);
//
//  heavy quark masses

	SetPMAS(4,1,1.2);
	SetMSTU(16,2);
//
//    primordial pT
	SetMSTP(91,1);
	SetPARP(91,1.);
	SetPARP(93,5.);
//
	break;
    case kPyBeauty:
	SetMSEL(5);
	SetPMAS(5,1,4.75);
	SetMSTU(16,2);
	break;
    case kPyJpsi:
	SetMSEL(0);
// gg->J/Psi g
	SetMSUB(86,1);
	break;
    case kPyJpsiChi:
	SetMSEL(0);
// gg->J/Psi g
	SetMSUB(86,1);
// gg-> chi_0c g
	SetMSUB(87,1);
// gg-> chi_1c g
	SetMSUB(88,1);
// gg-> chi_2c g
	SetMSUB(89,1);	
	break;
    case kPyCharmUnforced:
	SetMSEL(0);
// gq->qg   
	SetMSUB(28,1);
// gg->qq
	SetMSUB(53,1);
// gg->gg
	SetMSUB(68,1);
	break;
    case kPyBeautyUnforced:
	SetMSEL(0);
// gq->qg   
	SetMSUB(28,1);
// gg->qq
	SetMSUB(53,1);
// gg->gg
	SetMSUB(68,1);
	break;
    case kPyMb:
// Minimum Bias pp-Collisions
//
//   
//      select Pythia min. bias model
	SetMSEL(0);
	SetMSUB(92,1);      // single diffraction AB-->XB
	SetMSUB(93,1);      // single diffraction AB-->AX
	SetMSUB(94,1);      // double diffraction
	SetMSUB(95,1);	    // low pt production
	SetMSTP(81,1);      // multiple interactions switched on
	SetMSTP(82,3);      // model with varying impact param. & a single Gaussian
	SetPARP(82,3.47);   // set value pT_0  for turn-off of the cross section of                  
                            // multiple interaction at a reference energy = 14000 GeV
	SetPARP(89,14000.); // reference energy for the above parameter
	SetPARP(90,0.174);  // set exponent for energy dependence of pT_0
    case kPyMbNonDiffr:
// Minimum Bias pp-Collisions
//
//   
//      select Pythia min. bias model
	SetMSEL(0);
	SetMSUB(95,1);	    // low pt production
	SetMSTP(81,1);      // multiple interactions switched on
	SetMSTP(82,3);      // model with varying impact param. & a single Gaussian
	SetPARP(82,3.47);   // set value pT_0  for turn-off of the cross section of                  
                            // multiple interaction at a reference energy = 14000 GeV
	SetPARP(89,14000.); // reference energy for the above parameter
	SetPARP(90,0.174);  // set exponent for energy dependence of pT_0
 
	break;
    case kPyJets:
//
	printf("\n*************************************************\n");
	printf("\nWARNING !\n");
	printf("The kPyJet option uses simplified jet-production\n");
	printf("without gluon radiation \n");
	printf("\n*************************************************\n");
//
	SetMSEL(1);
// no initial state radiation   
	SetMSTP(61,0);
// no final state radiation
	SetMSTP(71,0);
// no primordial pT
	SetMSTP(91,0);
//	SetMSTP(111,0);	
	SetMSTU(16,1);	
	SetMSTJ(1,1);
	
	break;
    case kPyDirectGamma:
	SetMSEL(10);
	break;
    case kPyCharmPbMNR:
    case kPyD0PbMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // c-cbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with Pb-Pb collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ_4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.1GeV. Example in ConfigCharmPPR.C.

      // All QCD processes
      SetMSEL(1);

      // No multiple interactions
      SetMSTP(81,0);
      SetPARP(81,0.0);
      SetPARP(82,0.0);

      // Initial/final parton shower on (Pythia default)
      SetMSTP(61,1);
      SetMSTP(71,1);

      // 2nd order alpha_s
      SetMSTP(2,2);

      // QCD scales
      SetMSTP(32,2);
      SetPARP(34,1.0);

      // Intrinsic <kT^2>
      SetMSTP(91,1);
      SetPARP(91,1.304);
      SetPARP(93,6.52);

      // Set c-quark mass
      SetPMAS(4,1,1.2);

      break;
    case kPyBeautyPbMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // b-bbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with Pb-Pb collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.75GeV. Example in ConfigBeautyPPR.C.

      // All QCD processes
      SetMSEL(1);

      // No multiple interactions
      SetMSTP(81,0);
      SetPARP(81,0.0);
      SetPARP(82,0.0);

      // Initial/final parton shower on (Pythia default)
      SetMSTP(61,1);
      SetMSTP(71,1);

      // 2nd order alpha_s
      SetMSTP(2,2);

      // QCD scales
      SetMSTP(32,2);
      SetPARP(34,1.0);
      SetPARP(67,1.0);
      SetPARP(71,1.0);

      // Intrinsic <kT^2>
      SetMSTP(91,1);
      SetPARP(91,2.035);
      SetPARP(93,10.17);

      // Set b-quark mass
      SetPMAS(5,1,4.75);

      break;
    }
//
//  Initialize PYTHIA
    SetMSTP(41,1);   // all resonance decays switched on

    Initialize("CMS","p","p",fEcms);

}

Int_t AliPythia::CheckedLuComp(Int_t kf)
{
// Check Lund particle code (for debugging)
    Int_t kc=Pycomp(kf);
    printf("\n Lucomp kf,kc %d %d",kf,kc);
    return kc;
}

void AliPythia::SetNuclei(Int_t a1, Int_t a2)
{
// Treat protons as inside nuclei with mass numbers a1 and a2  
//    The MSTP array in the PYPARS common block is used to enable and 
//    select the nuclear structure functions. 
//    MSTP(52)  : (D=1) choice of proton and nuclear structure-function library
//            =1: internal PYTHIA acording to MSTP(51) 
//            =2: PDFLIB proton  s.f., with MSTP(51)  = 1000xNGROUP+NSET
//    If the following mass number both not equal zero, nuclear corrections of the stf are used.
//    MSTP(192) : Mass number of nucleus side 1
//    MSTP(193) : Mass number of nucleus side 2
    SetMSTP(52,2);
    SetMSTP(192, a1);
    SetMSTP(193, a2);  
}
	

AliPythia* AliPythia::Instance()
{ 
// Set random number generator 
    if (fgAliPythia) {
	return fgAliPythia;
    } else {
	fgAliPythia = new AliPythia();
	return fgAliPythia;
    }
}

void AliPythia::PrintParticles()
{ 
// Print list of particl properties
    Int_t np = 0;
    
    for (Int_t kf=0; kf<1000000; kf++) {
	for (Int_t c = 1;  c > -2; c-=2) {
	    
	    Int_t kc = Pycomp(c*kf);
	    if (kc) {
		Float_t mass  = GetPMAS(kc,1);
		Float_t width = GetPMAS(kc,2);	
		Float_t tau   = GetPMAS(kc,4);
		
		char*   name = new char[8];
		Pyname(kf,name);
	
		np++;
		
		printf("\n mass, width, tau: %6d %s %10.3f %10.3e %10.3e", 
		       c*kf, name, mass, width, tau);
	    }
	}
    }
    printf("\n Number of particles %d \n \n", np);
}

void  AliPythia::ResetDecayTable()
{
//  Set default values for pythia decay switches
    Int_t i;
    for (i = 1; i <  501; i++) SetMDCY(i,1,fDefMDCY[i]);
    for (i = 1; i < 2001; i++) SetMDME(i,1,fDefMDME[i]);
}

void  AliPythia::SetDecayTable()
{
//  Set default values for pythia decay switches
//
    Int_t i;
    for (i = 1; i <  501; i++) fDefMDCY[i] = GetMDCY(i,1);
    for (i = 1; i < 2001; i++) fDefMDME[i] = GetMDME(i,1);
}


#ifndef WIN32
#define pyr    pyr_
#define pyrset pyrset_
#define pyrget pyrget_
#else
#define pyr    PYR
#define pyrset PYRSET
#define pyrget PYRGET
#endif

extern "C" {
  Double_t pyr(Int_t*) 
{
      Float_t r;
      do r=sRandom->Rndm(); while(0 >= r || r >= 1);
      return r;
}
  void pyrset(Int_t*,Int_t*) {}
  void pyrget(Int_t*,Int_t*) {}
}




