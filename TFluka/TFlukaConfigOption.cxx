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

/* $Id$*/

#include "TFlukaConfigOption.h"
#include "TFlukaMCGeometry.h"
#include "TFluka.h"
#include "TFlukaCerenkov.h"

#include <TString.h>
#include <TObjArray.h>
#include <TVirtualMC.h>
#include <TGeoMaterial.h>

Float_t            TFlukaConfigOption::fgMatMin(-1.);
Float_t            TFlukaConfigOption::fgMatMax(-1.);
FILE*              TFlukaConfigOption::fgFile(0x0);
TFlukaMCGeometry*  TFlukaConfigOption::fgGeom(0x0);

Double_t TFlukaConfigOption::fgDCutValue[11];
Int_t    TFlukaConfigOption::fgDProcessFlag[15];


ClassImp(TFlukaConfigOption)


TFlukaConfigOption::TFlukaConfigOption()
{
    // Default constructor
    fMedium  = -1;
    fCMatMin = -1.;
    fCMatMax = -1.;    
    Int_t i;
    for (i = 0; i < 11; i++) fCutValue[i]    = -1.;
    for (i = 0; i < 15; i++) fProcessFlag[i] = -1;
}


TFlukaConfigOption::TFlukaConfigOption(Int_t medium)
{
    // Constructor
    fMedium = medium;
    fCMatMin = -1.;
    fCMatMax = -1.;    
    Int_t i;
    for (i = 0; i < 11; i++) fCutValue[i]    = -1.;
    for (i = 0; i < 15; i++) fProcessFlag[i] = -1;
}

void TFlukaConfigOption::SetCut(const char* flagname, Double_t val)
{
    // Set a cut value
    const TString cuts[11] = 
	{"CUTGAM", "CUTELE", "CUTNEU", "CUTHAD", "CUTMUO", "BCUTE", "BCUTM", "DCUTE", "DCUTM", "PPCUTM", "TOFMAX"};
    Int_t i;
    for (i = 0; i < 11; i++) {
	if (cuts[i].CompareTo(flagname) == 0) {
	    fCutValue[i] = val;
	    if (fMedium == -1) fgDCutValue[i] = val;
	    break;
	}
    }
}

void TFlukaConfigOption::SetProcess(const char* flagname, Int_t flag)
{
    // Set a process flag
    const TString process[15] = 
	{"DCAY", "PAIR", "COMP", "PHOT", "PFIS", "DRAY", "ANNI", "BREM", "MUNU", "CKOV", 
	 "HADR", "LOSS", "MULS", "RAYL", "STRA"};
    Int_t i;
    for (i = 0; i < 15; i++) {
	if (process[i].CompareTo(flagname) == 0) {
	    fProcessFlag[i] = flag;
	    if (fMedium == -1) fgDProcessFlag[i] = flag;
	    break;
	}
    }
}

void TFlukaConfigOption::WriteFlukaInputCards()
{
    // Write the FLUKA input cards for the set of process flags and cuts
    //
    //
    if (fMedium > -1) {
	TFluka* fluka = (TFluka*) gMC;
	TObjArray *matList = fluka->GetFlukaMaterials();
	Int_t nmaterial =  matList->GetEntriesFast();
	TGeoMaterial* material = 0;
	for (Int_t im = 0; im < nmaterial; im++)
	{
	    material = dynamic_cast<TGeoMaterial*> (matList->At(im));
	    Int_t idmat = material->GetIndex();
	    if (idmat == fMedium) break;	    
	}
	
	
//
// Check if global option

	fprintf(fgFile,"*\n*Material specific process and cut settings for #%8d %s\n", fMedium, material->GetName());
	fCMatMin = fMedium;
	fCMatMax = fMedium;
    } else {
	fprintf(fgFile,"*\n*Global process and cut settings \n");
	fCMatMin = fgMatMin;
	fCMatMax = fgMatMax;
    }

//
// Handle Process Flags 
//    
    if (fProcessFlag[kDCAY] != -1) ProcessDCAY();
    if (fProcessFlag[kPAIR] != -1) ProcessPAIR();
    if (fProcessFlag[kBREM] != -1) ProcessBREM();
    if (fProcessFlag[kCOMP] != -1) ProcessCOMP();
    if (fProcessFlag[kPHOT] != -1) ProcessPHOT();
    if (fProcessFlag[kPFIS] != -1) ProcessPFIS();
    if (fProcessFlag[kANNI] != -1) ProcessANNI();
    if (fProcessFlag[kMUNU] != -1) ProcessMUNU();
    if (fProcessFlag[kHADR] != -1) ProcessHADR();
    if (fProcessFlag[kMULS] != -1) ProcessMULS();
    if (fProcessFlag[kRAYL] != -1) ProcessRAYL();

    if (fProcessFlag[kLOSS] != -1 || fProcessFlag[kDRAY] != -1)                                    ProcessLOSS();
    if ((fMedium == -1 && fProcessFlag[kCKOV] > 0) || (fMedium > -1 && fProcessFlag[kCKOV] != -1)) ProcessCKOV();
    
//
// Handle Cuts
//
    if (fCutValue[kCUTGAM] >= 0.) ProcessCUTGAM();
    if (fCutValue[kCUTELE] >= 0.) ProcessCUTELE();
    if (fCutValue[kCUTNEU] >= 0.) ProcessCUTNEU();
    if (fCutValue[kCUTHAD] >= 0.) ProcessCUTHAD();
    if (fCutValue[kCUTMUO] >= 0.) ProcessCUTMUO();

    if (fCutValue[kTOFMAX] >= 0.) ProcessTOFMAX();
}

void TFlukaConfigOption::ProcessDCAY()
{
    // Process DCAY option
    fprintf(fgFile,"*\n* --- DCAY --- Decays. Flag = %5d\n", fProcessFlag[kDCAY]);
    if (fProcessFlag[kDCAY] == 0) {
	printf("Decays cannot be switched off \n");
    } else {
	fprintf(fgFile, "* Decays are on by default\n");
    }
}


void TFlukaConfigOption::ProcessPAIR()
{
    // Process PAIR option
    fprintf(fgFile,"*\n* --- PAIR --- Pair production by gammas, muons and hadrons. Flag = %5d, PPCUTM = %13.4g \n", 
	    fProcessFlag[kPAIR], fCutValue[kPPCUTM]);
    //
    // gamma -> e+ e-
    //
    if (fProcessFlag[kPAIR] > 0) {
	fprintf(fgFile,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fPHOT-THR\n",0., 0., 0.,    fCMatMin, fCMatMax, 1.);
    } else {
	fprintf(fgFile,"EMFCUT    %10.1f%10.1f%10.4g%10.1f%10.1f%10.1fPHOT-THR\n",0., 0., 1e10,  fCMatMin, fCMatMax, 1.);
    }
    
    //
    // Direct pair production by Muons and Hadrons
    //
    //
    // Attention ! This card interferes with BREM
    //

    if (fProcessFlag[kBREM] == -1 ) fProcessFlag[kBREM] = fgDProcessFlag[kBREM];
    if (fCutValue[kBCUTM]   == -1.) fCutValue[kBCUTM]   = fgDCutValue[kBCUTM];	


    Float_t flag = -3.;    
    if (fProcessFlag[kPAIR] >  0 && fProcessFlag[kBREM] == 0) flag =  1.;
    if (fProcessFlag[kPAIR] == 0 && fProcessFlag[kBREM]  > 0) flag =  2.;
    if (fProcessFlag[kPAIR] >  0 && fProcessFlag[kBREM]  > 0) flag =  3.;    
    if (fProcessFlag[kPAIR] == 0 && fProcessFlag[kBREM] == 0) flag = -3.;
    // Flag BREM card as handled
    fProcessFlag[kBREM] = -1;
    
    //
    // Energy cut for pair prodution
    //
    Float_t cutP = fCutValue[kPPCUTM];
    if (fCutValue[kPPCUTM]   == -1.) cutP = fgDCutValue[kPPCUTM];	
    // In G3 this is the cut on the total energy of the e+e- pair
    // In FLUKA the cut is on the kinetic energy of the electron and poistron
    cutP = cutP / 2. - 0.51099906e-3;
    if (cutP < 0.) cutP = 0.;
    // No explicite generation of e+/e-
    if (fProcessFlag[kPAIR] == 2) cutP = -1.;
    //
    // Energy cut for bremsstrahlung
    //
    Float_t cutB = 0.;
    if (flag > 1.) {
	fprintf(fgFile,"*\n* +++ BREM --- Bremsstrahlung by muons/hadrons. Flag = %5d, BCUTM = %13.4g \n",
	    fProcessFlag[kBREM], fCutValue[kBCUTM]);

	cutB =  fCutValue[kBCUTM];
	// No explicite production of gammas
	if (fProcessFlag[kBREM] == 2) cutB = -1.;
    }

    fprintf(fgFile,"PAIRBREM  %10.1f%10.4g%10.4g%10.1f%10.1f\n",flag, cutP, cutB, fCMatMin, fCMatMax);
}


void TFlukaConfigOption::ProcessBREM()
{
    // Process BREM option
    fprintf(fgFile,"*\n* --- BREM --- Bremsstrahlung by e+/- and muons/hadrons. Flag = %5d, BCUTE = %13.4g, BCUTM = %13.4g \n",
	    fProcessFlag[kBREM], fCutValue[kBCUTE], fCutValue[kBCUTM]);

    //
    // e+/- -> e+/- gamma
    //
    Float_t cutB = fCutValue[kBCUTE];
    if (fCutValue[kBCUTE]   == -1.) cutB = fgDCutValue[kBCUTE];	
    
    
    if (fProcessFlag[kBREM] > 0) {
	fprintf(fgFile,"EMFCUT    %10.4g%10.1f%10.1f%10.1f%10.1f%10.1fELPO-THR\n",cutB,  0., 0.,  fCMatMin, fCMatMax, 1.);
    } else {
	fprintf(fgFile,"EMFCUT    %10.4g%10.1f%10.1f%10.1f%10.1f%10.1fELPO-THR\n",1.e10, 0., 0.,  fCMatMin, fCMatMax, 1.);
    }

    //
    // Bremsstrahlung by muons and hadrons
    //
    cutB = fCutValue[kBCUTM];
    if (fCutValue[kBCUTM]   == -1.) cutB = fgDCutValue[kBCUTM];	
    if (fProcessFlag[kBREM] == 2) cutB = -1.;
    Float_t flag = 2.;
    if (fProcessFlag[kBREM] == 0) flag = -2.;
    
    fprintf(fgFile,"PAIRBREM  %10.1f%10.4g%10.4g%10.1f%10.1f\n", flag, 0., cutB, fCMatMin, fCMatMax);
}

void TFlukaConfigOption::ProcessCOMP()
{
    // Process COMP option
    fprintf(fgFile,"*\n* --- COMP --- Compton scattering  Flag = %5d \n", fProcessFlag[kCOMP]);

    //
    // Compton scattering
    //

    if (fProcessFlag[kCOMP] > 0) {
	fprintf(fgFile,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fPHOT-THR\n",0.   , 0., 0.,  fCMatMin, fCMatMax, 1.);
    } else {
	fprintf(fgFile,"EMFCUT    %10.4g%10.1f%10.1f%10.1f%10.1f%10.1fPHOT-THR\n",1.e10, 0., 0.,  fCMatMin, fCMatMax, 1.);
    }
}

void TFlukaConfigOption::ProcessPHOT()
{
    // Process PHOS option
    fprintf(fgFile,"*\n* --- PHOT --- Photoelectric effect. Flag = %5d\n", fProcessFlag[kPHOT]);

    //
    // Photoelectric effect
    //

    if (fProcessFlag[kPHOT] > 0) {
	fprintf(fgFile,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fPHOT-THR\n",0.   , 0., 0.,  fCMatMin, fCMatMax, 1.);
    } else {
	fprintf(fgFile,"EMFCUT    %10.1f%10.4g%10.1f%10.1f%10.1f%10.1fPHOT-THR\n",0., 1.e10, 0.,  fCMatMin, fCMatMax, 1.);
    }
}

void TFlukaConfigOption::ProcessANNI()
{
    // Process ANNI option
    fprintf(fgFile,"*\n* --- ANNI --- Positron annihilation. Flag = %5d \n", fProcessFlag[kANNI]);

    //
    // Positron annihilation
    //

    if (fProcessFlag[kANNI] > 0) {
	fprintf(fgFile,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fANNH-THR\n",0.   , 0., 0.,  fCMatMin, fCMatMax, 1.);
    } else {
	fprintf(fgFile,"EMFCUT    %10.4g%10.1f%10.1f%10.1f%10.1f%10.1fANNH-THR\n",1.e10, 0., 0.,  fCMatMin, fCMatMax, 1.);
    }
}


void TFlukaConfigOption::ProcessPFIS()
{
    // Process PFIS option
    fprintf(fgFile,"*\n* --- PFIS --- Photonuclear interaction  Flag = %5d\n", fProcessFlag[kPFIS]);

    //
    // Photonuclear interactions
    //

    if (fProcessFlag[kPFIS] > 0) {
	fprintf(fgFile,"PHOTONUC  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",(Float_t) fProcessFlag[kPFIS], 0., 0.,  fCMatMin, fCMatMax, 1.);
    } else {
	fprintf(fgFile,"PHOTONUC  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",-1.                          , 0., 0.,  fCMatMin, fCMatMax, 1.);
    }
}

void TFlukaConfigOption::ProcessMUNU()
{
    // Process MUNU option
    fprintf(fgFile,"*\n* --- MUNU --- Muon nuclear interaction. Flag = %5d\n", fProcessFlag[kMUNU]);
    
    //
    // Muon nuclear interactions
    //
    if (fProcessFlag[kMUNU] > 0) {
	fprintf(fgFile,"MUPHOTON  %10.1f%10.3f%10.3f%10.1f%10.1f%10.1f\n",(Float_t )fProcessFlag[kMUNU], 0.25, 0.75,  fCMatMin, fCMatMax, 1.);
    } else {
	fprintf(fgFile,"MUPHOTON  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",-1.                          , 0.,   0.,    fCMatMin, fCMatMax, 1.);
    }
}

void TFlukaConfigOption::ProcessRAYL()
{
    // Process RAYL option
    fprintf(fgFile,"*\n* --- RAYL --- Rayleigh Scattering. Flag = %5d\n", fProcessFlag[kRAYL]);

    //
    // Rayleigh scattering
    //
    Int_t nreg;
    Int_t* reglist = fgGeom->GetMaterialList(fMedium, nreg);
    //
    // Loop over regions of a given material
    for (Int_t k = 0; k < nreg; k++) {
	Float_t ireg = reglist[k];
	if (fProcessFlag[kRAYL] > 0) {
	    fprintf(fgFile,"EMFRAY    %10.1f%10.1f%10.1f%10.1f\n", 1., ireg, ireg, 1.);
	} else {
	    fprintf(fgFile,"EMFRAY    %10.1f%10.1f%10.1f%10.1f\n", 3., ireg, ireg, 1.);
	}
    }
}

void TFlukaConfigOption::ProcessCKOV()
{
    // Process CKOV option
    fprintf(fgFile,"*\n* --- CKOV --- Cerenkov Photon production.  %5d\n", fProcessFlag[kCKOV]);

    //
    // Cerenkov photon production
    //

    TFluka* fluka = (TFluka*) gMC;
    TObjArray *matList = fluka->GetFlukaMaterials();
    Int_t nmaterial =  matList->GetEntriesFast();
    for (Int_t im = 0; im < nmaterial; im++)
    {
	TGeoMaterial* material = dynamic_cast<TGeoMaterial*> (matList->At(im));
	Int_t idmat = material->GetIndex();
//
// Check if global option
	if (fMedium != -1 && idmat != fMedium) continue;
	
	TFlukaCerenkov* cerenkovProp;
	if (!(cerenkovProp = dynamic_cast<TFlukaCerenkov*>(material->GetCerenkovProperties()))) continue;
	//
	// This medium has Cerenkov properties 
	//
	//
	if (fMedium == -1 || (fMedium != -1 && fProcessFlag[kCKOV] > 0)) {
	    // Write OPT-PROD card for each medium 
	    Float_t  emin  = cerenkovProp->GetMinimumEnergy();
	    Float_t  emax  = cerenkovProp->GetMaximumEnergy();	      
	    fprintf(fgFile, "OPT-PROD  %10.4g%10.4g%10.4g%10.4g%10.4g%10.4gCERENKOV\n", emin, emax, 0., 
		    Float_t(idmat), Float_t(idmat), 0.); 
	    //
	    // Write OPT-PROP card for each medium 
	    // Forcing FLUKA to call user routines (queffc.cxx, rflctv.cxx, rfrndx.cxx)
	    //
	    fprintf(fgFile, "OPT-PROP  %10.4g%10.4g%10.4g%10.1f%10.1f%10.1fWV-LIMIT\n",  
		    cerenkovProp->GetMinimumWavelength(), cerenkovProp->GetMaximumWavelength(), cerenkovProp->GetMaximumWavelength(), 
		    Float_t(idmat), Float_t(idmat), 0.0);
	    

	    fprintf(fgFile, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n", -100., -100., -100., 
		    Float_t(idmat), Float_t(idmat), 0.0);
	    
	    for (Int_t j = 0; j < 3; j++) {
		fprintf(fgFile, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f&\n", -100., -100., -100., 
			Float_t(idmat), Float_t(idmat), 0.0);
	    }


	    // Photon detection efficiency user defined	    
	    if (cerenkovProp->IsSensitive())
		fprintf(fgFile, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fSENSITIV\n", -100., -100., -100., 
			Float_t(idmat), Float_t(idmat), 0.0);
	    // Material has a reflective surface
	    if (cerenkovProp->IsMetal())
		fprintf(fgFile, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fMETAL\n", -100., -100., -100., 
			Float_t(idmat), Float_t(idmat), 0.0);

	} else {
	    fprintf(fgFile,"OPT-PROD  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fCERE-OFF\n",0., 0., 0., fCMatMin, fCMatMax, 1.);
	}
    }
}


void TFlukaConfigOption::ProcessHADR()
{
    // Process HADR option
    fprintf(fgFile,"*\n* --- HADR --- Hadronic interactions. Flag = %5d\n", fProcessFlag[kHADR]);
    
    if (fProcessFlag[kHADR] > 0) {
	fprintf(fgFile,"*\n*Hadronic interaction is ON by default in FLUKA\n");
    } else {
	if (fMedium != -1) printf("Hadronic interactions cannot be switched off material by material !\n");
	fprintf(fgFile,"THRESHOL  %10.1f%10.1f%10.1f%10.1e%10.1f\n",0., 0., 0., 1.e10, 0.);
    }
}



void TFlukaConfigOption::ProcessMULS()
{
    // Process MULS option
    fprintf(fgFile,"*\n* --- MULS --- Muliple Scattering. Flag = %5d\n", fProcessFlag[kMULS]);
    //
    // Multiple scattering
    //
    if (fProcessFlag[kMULS] > 0) {
	fprintf(fgFile,"*\n*Multiple scattering is  ON by default in FLUKA\n");
    } else {
	fprintf(fgFile,"MULSOPT   %10.1f%10.1f%10.1f%10.1f%10.1f\n",0., 3., 3., fCMatMin, fCMatMax);
    }
}

void TFlukaConfigOption::ProcessLOSS()
{
    // Process LOSS option
    fprintf(fgFile,"*\n* --- LOSS --- Ionisation energy loss. Flags: LOSS = %5d, DRAY = %5d, STRA = %5d; Cuts: DCUTE = %13.4g, DCUTM = %13.4g \n",
	    fProcessFlag[kLOSS], fProcessFlag[kDRAY], fProcessFlag[kSTRA], fCutValue[kDCUTE], fCutValue[kDCUTM]);
    //
    // Ionisation energy loss
    //
    //
    // Impose consistency
    
    if (fProcessFlag[kLOSS] == 1 || fProcessFlag[kLOSS] == 3) {
	fProcessFlag[kDRAY] = 1;
    } else if (fProcessFlag[kLOSS] == 2) {
	fProcessFlag[kDRAY] = 0;
	fCutValue[kDCUTE] = 1.e10;
	fCutValue[kDCUTM] = 1.e10;	
    } else {
	if (fProcessFlag[kDRAY] == 1) {
	    fProcessFlag[kLOSS] = 1;
	} else if (fProcessFlag[kDRAY] == 0) {
	    fProcessFlag[kLOSS] = 2;
	    fCutValue[kDCUTE] = 1.e10;
	    fCutValue[kDCUTM] = 1.e10;	
	}
    }
    
    if (fCutValue[kDCUTE] == -1.) fCutValue[kDCUTE] = fgDCutValue[kDCUTE];
    if (fCutValue[kDCUTM] == -1.) fCutValue[kDCUTM] = fgDCutValue[kDCUTM];    
    
    Float_t cutM =  fCutValue[kDCUTM];    
	

    if (fProcessFlag[kSTRA] == -1) fProcessFlag[kSTRA] = fgDProcessFlag[kSTRA];
    if (fProcessFlag[kSTRA] ==  0) fProcessFlag[kSTRA] = 1;
    Float_t stra = (Float_t) fProcessFlag[kSTRA];
    

    if (fProcessFlag[kLOSS] == 1 || fProcessFlag[kLOSS] == 3) {
//
// Restricted energy loss fluctuations 
//
	fprintf(fgFile,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f\n", 1., 1., stra, fCMatMin, fCMatMax);
	fprintf(fgFile,"DELTARAY  %10.4g%10.1f%10.1f%10.1f%10.1f%10.1f\n", cutM, 0., 0., fCMatMin, fCMatMax, 1.);
    } else if (fProcessFlag[kLOSS] == 4) {
//
// No fluctuations
//
	fprintf(fgFile,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f\n",-1., -1., stra, fCMatMin, fCMatMax);	
	fprintf(fgFile,"DELTARAY  %10.4g%10.1f%10.1f%10.1f%10.1f%10.1f\n", 1.e10, 0., 0., fCMatMin, fCMatMax, 1.);	
    }  else {
//
// Full fluctuations
//
	fprintf(fgFile,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f\n",1., 1., stra, fCMatMin, fCMatMax);	
	fprintf(fgFile,"DELTARAY  %10.4g%10.1f%10.1f%10.1f%10.1f%10.1f\n", 1.e10, 0., 0., fCMatMin, fCMatMax, 1.);	
    }
}


void TFlukaConfigOption::ProcessCUTGAM()
{
// Cut on gammas
//
    fprintf(fgFile,"*\n*Cut for Gammas. CUTGAM = %13.4g\n", fCutValue[kCUTGAM]);
    if (fMedium == -1) {
	fprintf(fgFile,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n", 
		0., fCutValue[kCUTGAM], 0., 0., Float_t(fgGeom->NofVolumes()), 1.);
    } else {
	Int_t nreg, *reglist;
	Float_t ireg;
	reglist = fgGeom->GetMaterialList(fMedium, nreg);
	// Loop over regions of a given material
	for (Int_t k = 0; k < nreg; k++) {
	    ireg = reglist[k];
	    fprintf(fgFile,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n", 0.,fCutValue[kCUTGAM], 0., ireg, ireg, 1.);
	}
    }
}

void TFlukaConfigOption::ProcessCUTELE()
{
// Cut on e+/e-
//
    fprintf(fgFile,"*\n*Cut for e+/e-. CUTELE = %13.4g\n", fCutValue[kCUTELE]);
    if (fMedium == -1) {
	fprintf(fgFile,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n", 
		-fCutValue[kCUTELE], 0., 0., 0., Float_t(fgGeom->NofVolumes()), 1.);
    } else {
	Int_t nreg, *reglist;
	Float_t ireg;
	reglist = fgGeom->GetMaterialList(fMedium, nreg);
	// Loop over regions of a given material
	for (Int_t k = 0; k < nreg; k++) {
	    ireg = reglist[k];
	    fprintf(fgFile,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n", -fCutValue[kCUTELE], 0., 0., ireg, ireg, 1.);
	}
    }
}

void TFlukaConfigOption::ProcessCUTNEU()
{
    // Cut on neutral hadrons
    fprintf(fgFile,"*\n*Cut for neutral hadrons. CUTNEU = %13.4g\n", fCutValue[kCUTNEU]);
    if (fMedium == -1) {
	Float_t cut = fCutValue[kCUTNEU];
	//
	// 8.0 = Neutron
	// 9.0 = Antineutron
	//
	// If the cut is > 19.6 MeV it is assumed the low energy neutron transport is requested.
	// In this case the cut has to coincide with the upper  limit of the first energy group.
	//
	Float_t neutronCut = cut;
	if (neutronCut < 0.0196) {
	    neutronCut = 0.0196;
	    printf("Cut on neutron lower than upper limit of first energy group.\n");
	    printf("Cut reset to 19.6 MeV !\n");
	}
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -neutronCut,  8.0,  9.0);
	//
	// 12.0 = Kaon zero long
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 12.0, 12.0);
	// 17.0 = Lambda, 18.0 = Antilambda
	// 19.0 = Kaon zero short
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 17.0, 19.0);
	// 22.0 = Sigma zero, Pion zero, Kaon zero
	// 25.0 = Antikaon zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 22.0, 25.0);
	// 32.0 = Antisigma zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 32.0, 32.0);
	// 34.0 = Xi zero
	// 35.0 = AntiXi zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 34.0, 35.0);
	// 47.0 = D zero
	// 48.0 = AntiD zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 47.0, 48.0);
	// 53.0 = Xi_c zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 53.0, 53.0);
	// 55.0 = Xi'_c zero
	// 56.0 = Omega_c zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 55.0, 56.0);
	// 59.0 = AntiXi_c zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 59.0, 59.0);
	// 61.0 = AntiXi'_c zero
	// 62.0 = AntiOmega_c zero
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 61.0, 62.0);
    } else {
	printf("Cuts on neutral hadrons material by material not yet implemented !\n");
    }
}

void TFlukaConfigOption::ProcessCUTHAD()
{
    // Cut on charged hadrons
    fprintf(fgFile,"*\n*Cut for charge hadrons. CUTHAD = %13.4g\n", fCutValue[kCUTHAD]);
    if (fMedium == -1) {
	Float_t cut = fCutValue[kCUTHAD];
	// 1.0 = Proton
	// 2.0 = Antiproton
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut,  1.0,  2.0);
	// 13.0 = Positive Pion, Negative Pion, Positive Kaon
	// 16.0 = Negative Kaon
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 13.0, 16.0);
	// 20.0 = Negative Sigma
	// 21.0 = Positive Sigma
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 20.0, 21.0);
	// 31.0 = Antisigma minus
	// 33.0 = Antisigma plus
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 31.0, 31.0);
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 33.0, 33.0);
	// 36.0 = Negative Xi, Positive Xi, Omega minus
	// 39.0 = Antiomega
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 36.0, 39.0);
	// 45.0 = D plus
	// 46.0 = D minus
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 45.0, 46.0);
	// 49.0 = D_s plus, D_s minus, Lambda_c plus
	// 52.0 = Xi_c plus
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 49.0, 52.0);
	// 54.0 = Xi'_c plus
	// 60.0 = AntiXi'_c minus
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 54.0, 54.0);
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 60.0, 60.0);
	// 57.0 = Antilambda_c minus
	// 58.0 = AntiXi_c minus
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -cut, 57.0, 58.0);
    } else {
	printf("Cuts on charged hadrons material by material not yet implemented !\n");
    }
}

void TFlukaConfigOption::ProcessCUTMUO()
{
    // Cut on muons
    fprintf(fgFile,"*\n*Cut for muons. CUTMUO = %13.4g\n", fCutValue[kCUTMUO]);
    Float_t cut = fCutValue[kCUTMUO];
    if (fMedium == -1) {
	fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n",-cut, 10.0, 11.0);
    } else {
	printf("Cuts on muons material by material not yet implemented !\n");
    }
    
    
}

void TFlukaConfigOption::ProcessTOFMAX()
{
    // Cut time of flight
    Float_t cut = fCutValue[kTOFMAX];
    fprintf(fgFile,"*\n*Cut on time of flight. TOFMAX = %13.4g\n", fCutValue[kTOFMAX]);
    fprintf(fgFile,"TIME-CUT  %10.4g%10.1f%10.1f%10.1f%10.1f\n",cut*1.e9,0.,0.,-6.0,64.0);
}
