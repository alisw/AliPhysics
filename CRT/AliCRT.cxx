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
Revision 1.1  2002/06/16 17:08:19  hristov
First version of CRT


*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Cosmic Rays ALICE Trigger                                                //
//  This class contains the basic functions for the Cosmic Ray ALICE         //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//
// Begin_Html
/*
<img src="picts/AliCRTClass.gif">
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


#include <iostream.h>
#include <fstream.h>

#include <TTree.h>


#include "AliMC.h"
#include "AliRun.h"
#include "AliMagF.h"

#include "AliCRT.h"
#include "AliCRTConstants.h"
#include "AliCRTdigit.h"
#include "AliCRThit.h"

 
ClassImp(AliCRT)


static  AliCRTdigit *digits;
 
//_____________________________________________________________________________
AliCRT::AliCRT()
{
  //
  // Default constructor for L3 magnet
  //

  fIshunt   = 1;
  fHits     = 0;
  fDigits   = 0;
}
 
//_____________________________________________________________________________
AliCRT::AliCRT(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for the CRT module
  //

  fIshunt       =  1; // All hits are associated with primary particles  

  fHits         =  new TClonesArray("AliCRThit",400) ; 
  fDigits       =  new TClonesArray("AliCRTdigit",400) ; 

  gAlice->AddHitList(fHits);

  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);

}

//_____________________________________________________________________________
AliCRT::~AliCRT()
{
  //
  // Standar destructor.
  //
  if (fHits) {
    fHits->Delete();
    delete fHits;
  }
}

//_____________________________________________________________________________
void AliCRT::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a CRT hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliCRThit(fIshunt,track,vol,hits);
}

//_____________________________________________________________________________
void AliCRT::AddDigit(Int_t *tracks,Int_t *digits)
{
  
  //  Add a CRT digit to the list. Dummy function.
  
}

//_____________________________________________________________________________
void AliCRT::Init() const
{
  //
  // Initialise ...
  //

  Int_t i;
  //
  if(fDebug) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" CRT_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the CRT initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

//_____________________________________________________________________________
void AliCRT::ResetHits ()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits ();
}

//_____________________________________________________________________________
void AliCRT::ResetDigits ()
{
  //
  // Reset number of digits and the digits array for this detector
  AliDetector::ResetDigits ();
  //
} 

//____________________________________________________________________________
void AliCRT::FinishEvent()
{
// do nothing
}

//_____________________________________________________________________________
void AliCRT::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //
}

//_____________________________________________________________________________
void AliCRT::CreateGeometry()
{
  //
  // Build simple ROOT TNode geometry for GEANT simulations
  //
}

//_____________________________________________________________________________
void AliCRT::CreateMaterials()
{
  // Magnatic field inside the pit
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();

  //Magnetic field above the Magnet.
  Int_t xfield = 0;   // no Magnetic field.
  Float_t xfieldm = 0;
  Float_t xepsil = 0.1; // Tracking precission in cm. obove the pit

  // --- Define the various materials for GEANT --- 
  Float_t aconc[10] = { 1.,12.01,15.994,22.99,24.305,26.98,28.086,39.1,40.08,55.85 };
  Float_t zconc[10] = { 1.,6.,8.,11.,12.,13.,14.,19.,20.,26. };
  Float_t wconc[10] = { .01,.001,.529107,.016,.002,.033872,.337021,.013,.044,.014 };
  
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  //
  //     Aluminum 
  AliMaterial(9,  "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  //
  //     Iron 
  AliMaterial(10, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  //
  //     Air 
  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  AliMaterial(35, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  AliMaterial(55, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  //
  //     Vacuum 
  AliMaterial(16, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(36, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(56, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  //
  //     Concrete
  AliMixture(17, "CONCRETE$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(37, "CONCRETE$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(57, "CONCRETE$", aconc, zconc, 2.35, 10, wconc);


  // Scintillator material polystyrene 
  Float_t aP[2] = {12.011, 1.00794};
  Float_t zP[2] = {6.0, 1.0};
  Float_t wP[2] = {1.0, 1.0};
  Float_t dP = 1.032;
  AliMixture(3, "Polystyrene$", aP, zP, dP, -2, wP);
  // Standard ROCK, ROCK1 and molasse. 
  Float_t aRock[1] = {22.};
  Float_t zRock[1] = {11.};
  Float_t wRock[1] = {1.};
  Float_t dRock = 2.65;
  AliMixture(4, "Standard Rock$", aRock, zRock, dRock, 1, wRock);
  Float_t dRock1 = 2.2;
  AliMixture(5, "Rock1$", aRock, zRock, dRock1, 1, wRock);
  Float_t dMolasse = 2.50;
  AliMixture(6, "Molasse$", aRock, zRock, dMolasse, 1, wRock);
  
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

  //
  //    Aluminum 
  AliMedium(9,  "ALU_C0          ",  9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(29, "ALU_C1          ", 29, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(49, "ALU_C2          ", 49, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Iron 
  AliMedium(10, "FE_C0           ", 10, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1           ", 30, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(50, "FE_C2           ", 50, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Air 
  AliMedium(15, "AIR_C0          ", 15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_C1          ", 35, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(55, "AIR_C2          ", 55, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Vacuum 
  AliMedium(16, "VA_C0           ", 16, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(36, "VA_C1           ", 36, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(56, "VA_C2           ", 56, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Concrete 
  AliMedium(17, "CC_C0           ", 17, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(37, "CC_C1           ", 37, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(57, "CC_C2           ", 57, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);


  // The scintillator of the CPV made of Polystyrene 
  // scintillator -> idtmed[1102]
  AliMedium(3, "CPV scint.          ", 3, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  
  //     Standard ROCK  -> idtmed[1103]
  AliMedium(4, "Standard Rock ",  4, 0, xfield, xfieldm, tmaxfd, stemax, deemax, xepsil, stmin);
  
  //     ROCK1   -> idtmed[1104]
  AliMedium(5, "Rock 1                 ", 5, 0, xfield, xfieldm, tmaxfd, stemax, deemax, xepsil, stmin);


  //     Molasse -> idtmed[1105]
  AliMedium(6, "Molasse          ", 6, 0, xfield, xfieldm, tmaxfd, stemax, deemax, xepsil, stmin);


}

//_____________________________________________________________________________
void AliCRT::MakeBranch(Option_t* option, const char *file)
{
  //
  // Specific CRT branches
  //
  // Create Tree branches for the CRT.
  Int_t buffersize = 400;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option,file);

  const char *cD = strstr(option,"D");
  
  if (cD) {
    digits = new AliCRTdigit();
    MakeBranchInTree(gAlice->TreeD(), branchname, "AliCRTdigit", 
		     digits, buffersize, 1, file);
  } 

}
//_____________________________________________________________________________
