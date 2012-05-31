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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Cosmic Rays ALICE Trigger                                                //
//  This class contains the basic functions for the Cosmic Ray ALICE         //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//
// Begin_Html
/*
<img src="picts/AliACORDEClass.gif">
</pre>
<p>The responsible person for this module is
<a href="mailto:Enrique.Gamez.Flores@cern.ch">Enrique Gamez Flores</a>.
</font>
<pre>
*/
//End_Html
//             
//
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliACORDE.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliACORDERawData.h"
#include "AliACORDERawStream.h"

ClassImp(AliACORDE)

//_____________________________________________________________________________
AliACORDE::AliACORDE()
  : AliDetector(),
    fCreateCavern(0),
    f4CentralModulesGeometry(0)
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliACORDE::AliACORDE(const char *name, const char *title)
  : AliDetector(name, title),
    fCreateCavern(kFALSE),
    f4CentralModulesGeometry(kTRUE)

{
  //
  // Standard constructor
}

//_____________________________________________________________________________
AliACORDE::~AliACORDE()
{
  //
  // Default destructor
  //
}

//_____________________________________________________________________________
void AliACORDE::CreateMaterials()
{
  // Magnatic field inside the pit
  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();

  //Magnetic field above the Magnet.
  Int_t xfield = 0;   // no Magnetic field.
  Float_t xfieldm = 0;
  Float_t xepsil = 0.1; // Tracking precission in cm. obove the pit

  // --- Define the various materials for GEANT --- 
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  //
  //     Aluminum 
  AliMaterial(9,  "ALUMINIUM0$", 26.98, 13., 2.7, 8.9, 37.2);
  //  AliMaterial(29, "ALUMINIUM1$", 26.98, 13., 2.7, 8.9, 37.2);
  //AliMaterial(49, "ALUMINIUM2$", 26.98, 13., 2.7, 8.9, 37.2);
  //
  //     Iron 
  //AliMaterial(10, "IRON0$    ", 55.85, 26., 7.87, 1.76, 17.1);
  //AliMaterial(30, "IRON1$    ", 55.85, 26., 7.87, 1.76, 17.1);
  //AliMaterial(50, "IRON2$    ", 55.85, 26., 7.87, 1.76, 17.1);
  //
  //     Air 
  Float_t as[] = { 12.0107, 14.0067,   15.9994,  39.948 };
  Float_t zs[] = {  6.,      7.,       8.,       18. };
  Float_t ws[] = { 0.000124, 0.755267, 0.231781, 0.012827 }; 
  Double_t density      = .00120479;
  AliMixture(15, "AIR0$", as, zs, density, 4, ws);

  //AliMaterial(15, "AIR0$     ", 14.61, 7.3, .001205, 30423.24, 67500.);
  //AliMaterial(35, "AIR1$     ", 14.61, 7.3, .001205, 30423.24, 67500.);
  //AliMaterial(55, "AIR2$     ", 14.61, 7.3, .001205, 30423.24, 67500.);
  //AliMaterial(75, "AIR3$     ", 14.61, 7.3, .001205, 30423.24, 67500.);
  //AliMaterial(95, "AIR4$     ", 14.61, 7.3, .001205, 30423.24, 67500.);


  // Scintillator material polystyrene 
  Float_t aP[2] = {12.011, 1.00794};
  Float_t zP[2] = {6.0, 1.0};
  Float_t wP[2] = {1.0, 1.0};
  Float_t dP = 1.032;
  AliMixture(13, "Polystyrene$", aP, zP, dP, -2, wP);
  // Subalpine Molasse over the ALICE hall. 
  Float_t aMolasse[10] = { 1., 12.01, 15.994, 22.99, 24.305, 26.98, 28.086, 39.1, 40.08, 55.85 };
  Float_t zMolasse[10] = {1., 6., 8., 11., 12., 13., 14., 19., 20., 26.};
  Float_t wMolasse[10] = {0.008, 0.043, 0.485, 0.007, 0.042, 0.037, 0.215, 0.023, 0.1, 0.04};
  Float_t dMolasse = 2.40;
  AliMixture(24, "Molasse$", aMolasse, zMolasse, 2*dMolasse, 10, wMolasse); // correction to density of molasse alpine newDen=2*dmolasse

  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001;  // Tracking precision, Inside the pit
  stemax = -1.;   // Maximum displacement for multiple scattering 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 

  Float_t atmaxfd = 10.;
  Float_t adeemax = -0.1;
  Float_t aepsil = 0.1;
  Float_t astmin = -10.;

  //
  //    Aluminum 
  AliMedium(9,  "ALU_C0          ",  9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  // AliMedium(29, "ALU_C1          ", 29, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //AliMedium(49, "ALU_C2          ", 49, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Iron 
  //AliMedium(10, "FE_C0           ", 10, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //AliMedium(30, "FE_C1           ", 30, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //AliMedium(50, "FE_C2           ", 50, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Air 
  AliMedium(15, "AIR_C0          ", 15, 0, isxfld, sxmgmx, atmaxfd, stemax, adeemax, aepsil, astmin);
  // AliMedium(35, "AIR_C1          ", 35, 0, isxfld, sxmgmx, atmaxfd, stemax, adeemax, aepsil, astmin);
  //AliMedium(55, "AIR_C2          ", 55, 0, isxfld, sxmgmx, atmaxfd, stemax, adeemax, aepsil, astmin);
  //AliMedium(75, "AIR_C4          ", 75, 0, isxfld, sxmgmx, atmaxfd, stemax, adeemax, aepsil, astmin);
  //AliMedium(95, "AIR_C5          ", 95, 0, isxfld, sxmgmx, atmaxfd, stemax, adeemax, aepsil, astmin);



  // The scintillator of the CPV made of Polystyrene 
  // scintillator -> idtmed[1112]
  //AliMedium(12 , "CPV scint.0     ", 13, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(13 , "CPV scint.1     ", 13, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  //AliMedium(14 , "CPV scint.2     ", 13, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);

  //     Molasse -> idtmed[1123]
  AliMedium(24 , "Molasse         ", 24, 0, xfield, xfieldm, tmaxfd, stemax, deemax, xepsil, stmin);

  // Concrete, in case if we need to put hte shafts by ourselves.

  Float_t aconc[10] = { 1.,12.01,15.994,22.99,24.305,26.98,28.086,39.1,40.08,55.85 };
  Float_t zconc[10] = { 1.,6.,8.,11.,12.,13.,14.,19.,20.,26. };
  Float_t wconc[10] = { .01,.001,.529107,.016,.002,.033872,.337021,.013,.044,.014 };

  AliMixture(17, "CONCRETE$", aconc, zconc, 2.35, 10, wconc);
  //    Concrete 
  AliMedium(17, "CC_C0            ", 17, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //  AliMedium(27, "CC_C1            ", 17, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin); // MX24
  //AliMedium(37, "CC_C2            ", 17, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin); // PM25
  //AliMedium(47, "CC_C3            ", 17, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin); // PGC2

}

//_____________________________________________________________________________
void AliACORDE::SetTreeAddress()
{

  TBranch *branch;
  char branchname[20];
  snprintf(branchname,19,"%s",GetName());
  // Branch address for hit tree
  TTree *treeH = fLoader->TreeH();
  if (treeH ) {
    branch = treeH->GetBranch(branchname);
    if (branch) branch->SetAddress(&fHits);
  }
}

//_____________________________________________________________________________
void AliACORDE::MakeBranch(Option_t* opt)
{
  //
  // Initializes the branches of the ACORDE inside the trees written
  // for each event.
  //
  const char* oH = strstr(opt, "H");
  if ( fLoader->TreeH() && oH && (fHits == 0x0) ) {
    fHits = new TClonesArray("AliACORDEhit", 1000);
    fNhits = 0;
  }
  AliDetector::MakeBranch(opt);
}

AliLoader* AliACORDE::MakeLoader(const char* topfoldername)
{ 
 
  AliDebug(1,Form("Creating AliACORDELoader, Top folder is %s ",
		  topfoldername));
  fLoader = new AliACORDELoader(GetName(),topfoldername);
  return fLoader;
}


AliDigitizer* AliACORDE::CreateDigitizer(AliDigitizationInput* digInput) const
{
  //
  //
  return new AliACORDEDigitizer(digInput);
}

void AliACORDE::Digits2Raw()
{
  // Produce Raw data starting from digits
  // 1. Get digits
  // 2. From digits get an array with the state of the modules
  // 3. Unload digits
  // 4. Write raw data

  // 1. Get digits

  // 1.1 Get detector, load digits and set branch
  AliACORDE* acorde = (AliACORDE*)gAlice->GetDetector("ACORDE");
  fLoader->LoadDigits("READ");
  TTree* treeD = fLoader->TreeD();
  if (!treeD) {
    Error("Digits2Raw", "no digits tree");
    return;
  }
  TClonesArray *adigits = new TClonesArray ("AliACORDEdigit", 1000);
  treeD->GetBranch("ACORDEdigit")->SetAddress(&adigits);
  // 1.2 Get first entry (there is always only one)
  acorde->ResetDigits();
  treeD->GetEvent(0);
  
  // 2. From digits get an array with the state of the modules
  // 2.1 Define and initialize the array
  Bool_t Modules[60];
  for (Int_t i=0;i<60;i++) Modules[i]= kFALSE;
  // 2.2 Loop over all digits
  Int_t ndig = adigits->GetEntriesFast();
  for (Int_t idig=0;idig<ndig;idig++) {
    // 2.3 set the array entry for each digit
    AliACORDEdigit* digit = (AliACORDEdigit*) adigits->At(idig);
    Int_t mod = digit->GetModule();
    Modules[mod]=kTRUE;
  } 
  // 3. Unload digits
  fLoader->UnloadDigits();

  // 4. Write raw data
  AliACORDERawData rawdata;
  rawdata.WriteACORDERawData(Modules,(ndig > 1));
}

//_____________________________________________________________________________
Bool_t AliACORDE::Raw2SDigits(AliRawReader* rawReader)
{
  //
  // Reads the raw data stream and exracts the digits
  //
  // Input:
  //         rawReader : pointer to the current AliRawReader
  // Output:
  //
  // Created:      31 Jan 2008  Mario Sitta
  //
  TStopwatch timer;
  timer.Start();

  if(!fLoader) {
    AliError("no ACORDE loader found");
    return kFALSE;
  }

  TTree* treeD  = fLoader->TreeD();
  if(!treeD) {
      fLoader->MakeTree("D");
      treeD = fLoader->TreeD();
  }
        
  AliACORDEdigit  digit;
  AliACORDEdigit* pdigit = &digit;
  const Int_t kBufferSize = 4000;
   
  treeD->Branch("ACORDE", "AliACORDEdigit",  &pdigit, kBufferSize);

  //  rawReader->Reset();
  AliACORDERawStream* rawStream  = new AliACORDERawStream(rawReader);    
     
  if (!rawStream->Next()) return kFALSE; // No ACORDE data found
  /*  
  for(Int_t i=0; i<64; i++) {
      new(pdigit) AliACORDEdigit(i, (Int_t)rawStream->GetADC(i), (Int_t)rawStream->GetTime(i)); 
      treeD->Fill();
  }
  */ 
  fLoader->WriteDigits("OVERWRITE");
  fLoader->UnloadDigits();	
	
  delete rawStream;

  timer.Stop();
  timer.Print();

  return kTRUE;
}

//_____________________________________________________________________________
