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
#include <TList.h>
#include <TObjArray.h>
#include <TVirtualMC.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoManager.h>
#include <TGeoMedium.h>

Float_t            TFlukaConfigOption::fgMatMin(-1.);
Float_t            TFlukaConfigOption::fgMatMax(-1.);
FILE*              TFlukaConfigOption::fgFile(0x0);
TFlukaMCGeometry*  TFlukaConfigOption::fgGeom(0x0);

Double_t TFlukaConfigOption::fgDCutValue[11];
Int_t    TFlukaConfigOption::fgDProcessFlag[15];


ClassImp(TFlukaConfigOption)


TFlukaConfigOption::TFlukaConfigOption()
  :fMedium(-1),
   fCMatMin(-1),
   fCMatMax(-1),
   fCMaterial(0)
{
    // Default constructor
//    fMedium  = -1;
//    fCMatMin = -1.;
//    fCMatMax = -1.;    
    Int_t i;
    for (i = 0; i < 11; i++) fCutValue[i]    = -1.;
    for (i = 0; i < 15; i++) fProcessFlag[i] = -1;
}


TFlukaConfigOption::TFlukaConfigOption(Int_t medium)
  :fMedium(medium),
   fCMatMin(-1),
   fCMatMax(-1),
   fCMaterial(0)
{
    // Constructor
//    fMedium = medium;
//    fCMatMin = -1.;
//    fCMatMax = -1.;    
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

void TFlukaConfigOption::SetModelParameter(const char* flagname, Double_t val)
{
    // Set a model parameter value
    const TString parms[2] = {"PRIMIO_N", "PRIMIO_E"};
    Int_t i;
    for (i = 0; i < 2; i++) {
       if (parms[i].CompareTo(flagname) == 0) {
           fModelParameter[i] = val;
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
    //
    // Check if global option or medium specific

    Bool_t mediumIsSensitive = kFALSE;
    TGeoMedium*   med    = 0x0;
    TGeoMedium*   medium = 0x0;    
    TGeoMaterial* mat    = 0x0;

    if (fMedium != -1) {
       TFluka* fluka = (TFluka*) gMC;
       fMedium = fgGeom->GetFlukaMaterial(fMedium);
       //
       // Check if material is actually used
       Int_t* reglist;
       Int_t nreg;
       reglist = fgGeom->GetMaterialList(fMedium, nreg);
       if (nreg == 0) {
           // Material not used -- return
           return;
       }
       //
       // Find material
       TObjArray *matList = fluka->GetFlukaMaterials();
       Int_t nmaterial =  matList->GetEntriesFast();
       fCMaterial = 0;
       for (Int_t im = 0; im < nmaterial; im++)
       {
           fCMaterial = dynamic_cast<TGeoMaterial*> (matList->At(im));
           Int_t idmat = fCMaterial->GetIndex();
           if (idmat == fMedium) break;
       } // materials
        //
       // Find medium
       TList *medlist = gGeoManager->GetListOfMedia();
       TIter next(medlist);
       while((med = (TGeoMedium*)next()))
       {
           mat = med->GetMaterial();
           if (mat->GetIndex() == fMedium) {
              medium = med;
              break;
           }
       } // media
       //
       // Check if sensitive
       if (medium->GetParam(0) != 0.) mediumIsSensitive = kTRUE;


       fprintf(fgFile,"*\n*Material specific process and cut settings for #%8d %s\n", fMedium, fCMaterial->GetName());
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
//
//  First make sure that all cuts are taken into account
    if (DefaultProcessFlag(kPAIR) > 0 && fProcessFlag[kPAIR] == -1 && (fCutValue[kCUTELE] >= 0. || fCutValue[kPPCUTM] >= 0.)) 
       fProcessFlag[kPAIR] = DefaultProcessFlag(kPAIR);
    if (DefaultProcessFlag(kBREM) > 0 && fProcessFlag[kBREM] == -1 && (fCutValue[kBCUTE]  >= 0. || fCutValue[kBCUTM] >= 0.)) 
       fProcessFlag[kBREM] = DefaultProcessFlag(kBREM);
    if (DefaultProcessFlag(kDRAY) > 0 && fProcessFlag[kDRAY] == -1 && (fCutValue[kDCUTE]  >= 0. || fCutValue[kDCUTM] >= 0.)) 
       fProcessFlag[kDRAY] = DefaultProcessFlag(kDRAY);
//    
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
//
//  Time of flight 
    if (fCutValue[kTOFMAX] >= 0.) ProcessTOFMAX();

//
//  Tracking precission
    if (mediumIsSensitive) ProcessSensitiveMedium();
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
    fprintf(fgFile,"*\n* --- PAIR --- Pair production by gammas, muons and hadrons. Flag = %5d, PPCUTM = %13.4g, PPCUTE = %13.4g\n", 
           fProcessFlag[kPAIR], fCutValue[kCUTELE], fCutValue[kPPCUTM]);
    //
    // gamma -> e+ e-
    //
    if (fProcessFlag[kPAIR] > 0) {
       fprintf(fgFile,"EMFCUT    %10.1f%10.1f%10.4g%10.1f%10.1f%10.1fPHOT-THR\n",0., 0., 0.0, fCMatMin, fCMatMax, 1.);
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
    fProcessFlag[kBREM] = - fProcessFlag[kBREM];
    
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
    
    
    if (TMath::Abs(fProcessFlag[kBREM]) > 0) {
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
       fprintf(fgFile,"EMFCUT    %10.4g%10.4g%10.4g%10.1f%10.1f%10.1fPHOT-THR\n",0., 0., 0.,  fCMatMin, fCMatMax, 1.);
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
    
    if (fProcessFlag[kLOSS] == 1 || fProcessFlag[kLOSS] == 3 || fProcessFlag[kLOSS] > 10) {
    // Restricted fluctuations
       fProcessFlag[kDRAY] = 1;
    } else if (fProcessFlag[kLOSS] == 2) {
    // Full fluctuations
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
    } else if (fProcessFlag[kLOSS] > 10) {
//
// Primary ionisation electron generation
//
       // Ionisation model
       Float_t ioModel = Float_t (fProcessFlag[kLOSS]-10);
       //  Effective 1st ionisation potential
       Float_t ePot    = ModelParameter(kPRIMIOE);
       // Number of primary ionisations per cm for a mip
       Float_t nPrim   = ModelParameter(kPRIMION);
       
       fprintf(fgFile,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f\n", 1., 1., stra, fCMatMin, fCMatMax);
       fprintf(fgFile,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fPRIM-ION\n", ePot, nPrim, ioModel, fCMatMin, fCMatMax, 1.);
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
    
    // Transport production cut used for pemf
    //
    //  FUDGEM paramter. The parameter takes into account th contribution of atomic electrons to multiple scattering.
    //  For production and transport cut-offs larger than 100 keV it must be set = 1.0, while in the keV region it must be
    Float_t parFudgem = (fCutValue[kCUTGAM] > 1.e-4)? 1.0 : 0.0 ;
    fprintf(fgFile,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1fPROD-CUT\n", 
           0., fCutValue[kCUTGAM], parFudgem, fCMatMin, fCMatMax, 1.);
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
    // Transport production cut used for pemf
    //
    //  FUDGEM paramter. The parameter takes into account th contribution of atomic electrons to multiple scattering.
    //  For production and transport cut-offs larger than 100 keV it must be set = 1.0, while in the keV region it must be
    Float_t parFudgem = (fCutValue[kCUTELE] > 1.e-4)? 1.0 : 0.0;
    fprintf(fgFile,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1fPROD-CUT\n", 
           -fCutValue[kCUTELE], 0., parFudgem, fCMatMin, fCMatMax, 1.);
}

void TFlukaConfigOption::ProcessCUTNEU()
{
  // Cut on neutral hadrons
  fprintf(fgFile,"*\n*Cut for neutral hadrons. CUTNEU = %13.4g\n", fCutValue[kCUTNEU]);
  
  // Cut on neutral hadrons
  fprintf(fgFile,"*\n*Cut for neutral hadrons. CUTNEU = %13.4g\n", fCutValue[kCUTNEU]);
  
  // Energy group structure of the 72-neutron ENEA data set:
  const Float_t neutronGroupUpLimit[72] = {
    1.9600E-02, 1.7500E-02, 1.4918E-02, 1.3499E-02, 1.2214E-02, 1.1052E-02, 1.0000E-02, 9.0484E-03,
    8.1873E-03, 7.4082E-03, 6.7032E-03, 6.0653E-03, 5.4881E-03, 4.9659E-03, 4.4933E-03, 4.0657E-03,
    3.6788E-03, 3.3287E-03, 3.0119E-03, 2.7253E-03, 2.4660E-03, 2.2313E-03, 2.0190E-03, 1.8268E-03,
    1.6530E-03, 1.4957E-03, 1.3534E-03, 1.2246E-03, 1.1080E-03, 1.0026E-03, 9.0718E-04, 8.2085E-04,
    7.4274E-04, 6.0810E-04, 4.9787E-04, 4.0762E-04, 3.3373E-04, 2.7324E-04, 2.2371E-04, 1.8316E-04,
    1.4996E-04, 1.2277E-04, 8.6517E-05, 5.2475E-05, 3.1828E-05, 2.1852E-05, 1.5034E-05, 1.0332E-05,
    7.1018E-06, 4.8809E-06, 3.3546E-06, 2.3054E-06, 1.5846E-06, 1.0446E-06, 6.8871E-07, 4.5400E-07,
    2.7537E-07, 1.6702E-07, 1.0130E-07, 6.1442E-08, 3.7267E-08, 2.2603E-08, 1.5535E-08, 1.0677E-08,
    7.3375E-09, 5.0435E-09, 3.4662E-09, 2.3824E-09, 1.6374E-09, 1.1254E-09, 6.8257E-10, 4.1400E-10
  };

  Float_t cut = fCutValue[kCUTNEU];
  //
  // 8.0 = Neutron
  // 9.0 = Antineutron
  // Find the FLUKA neutron group corresponding to the cut
  //
  Float_t neutronCut = cut;
  Int_t groupCut = 1; // if cut is > 19.6 MeV not low energy neutron transport is done
  if (neutronCut < 0.0196) {
    neutronCut = 0.0196;
    // Search the group cutoff for the energy cut
    Int_t i;
    for( i=0; i<72; i++ ) {
      if (cut > neutronGroupUpLimit[i]) break;
    }
    groupCut = i+1;
  } 
  
  if (fMedium == -1) {
        Float_t cut = fCutValue[kCUTNEU];
        // 8.0 = Neutron
        // 9.0 = Antineutron
        fprintf(fgFile,"PART-THR  %10.4g%10.1f%10.1f\n", -neutronCut,  8.0,  9.0);
        fprintf(fgFile,"LOW-BIAS  %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n",
              Float_t(groupCut), 73.0, 0.95, 2., Float_t(fgGeom->NofVolumes()), 1.);
        //
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
        Int_t nreg, *reglist;
        Float_t ireg;
        reglist = fgGeom->GetMaterialList(fMedium, nreg);
        // Loop over regions of a given material
        for (Int_t k = 0; k < nreg; k++) {
         ireg = reglist[k];
         fprintf(fgFile,"LOW-BIAS  %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n",
                Float_t(groupCut), 73.0, 0.95, ireg, ireg, 1.);
       }

       Warning("ProcessCUTNEU",
              "Material #%4d %s: Cut on neutral hadrons (Ekin > %9.3e) material by material only implemented for low-energy neutrons !\n",
              fMedium, fCMaterial->GetName(), cut);
    }
}

void TFlukaConfigOption::ProcessCUTHAD()
{
    // Cut on charged hadrons
    fprintf(fgFile,"*\n*Cut for charge hadrons. CUTHAD = %13.4g\n", fCutValue[kCUTHAD]);
    Float_t cut = fCutValue[kCUTHAD];
    if (fMedium == -1) {
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
      Warning("ProcessCUTHAD", 
              "Material #%4d %s: Cut on charged hadrons (Ekin > %9.3e) material by material not yet implemented !\n",
             fMedium, fCMaterial->GetName(), cut);
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
       Warning("ProcessCUTMUO", "Material #%4d %s: Cut on muons (Ekin > %9.3e) material by material not yet implemented !\n",
              fMedium, fCMaterial->GetName(), cut);
    }
    
    
}

void TFlukaConfigOption::ProcessTOFMAX()
{
    // Cut time of flight
    Float_t cut = fCutValue[kTOFMAX];
    fprintf(fgFile,"*\n*Cut on time of flight. TOFMAX = %13.4g\n", fCutValue[kTOFMAX]);
    fprintf(fgFile,"TIME-CUT  %10.4g%10.1f%10.1f%10.1f%10.1f\n",cut*1.e9,0.,0.,-6.0,64.0);
}

void  TFlukaConfigOption::ProcessSensitiveMedium()
{
    //
    // Special options for sensitive media
    //

    fprintf(fgFile,"*\n*Options for sensitive medium \n");
    //
    // EMFFIX
    fprintf(fgFile,"EMFFIX    %10.1f%10.3f%10.1f%10.1f%10.1f%10.1f\n", fCMatMin, 0.05, 0., 0., 0., 0.);
    //
    // FLUKAFIX
    fprintf(fgFile,"FLUKAFIX  %10.3f                    %10.3f\n", 0.05, fCMatMin);
}
