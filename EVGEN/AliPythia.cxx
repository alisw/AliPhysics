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
#include "AliRun.h"

ClassImp(AliPythia)

//_____________________________________________________________________________

AliPythia* AliPythia::fgAliPythia=NULL;

AliPythia::AliPythia()
{
// Default Constructor
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
    SetMSTP(41,1);

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
//            =3: PDFLIB proton  s.f. with nuclar correction:
//                MSTP( 51)  = 1000xNPGROUP+NPSET
//                MSTP(151)  = 1000xNAGROUP+NASET
//    MSTP(192) : Mass number of nucleus side 1
//    MSTP(193) : Mass number of nucleus side 2

    SetMSTP(52,3);
    SetMSTP(191, 1001);
    SetMSTP(192, a1);
    SetMSTP(193, a2);  
}
	

AliPythia* AliPythia::Instance()
{
    if (fgAliPythia) {
	return fgAliPythia;
    } else {
	fgAliPythia = new AliPythia();
	return fgAliPythia;
    }
}
void AliPythia::Streamer(TBuffer &R__b) {} 







