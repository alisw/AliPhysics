///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  L3 Magnet                                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliMAGClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include "AliMAG.h"
#include "AliRun.h"
 
ClassImp(AliMAG)
 
//_____________________________________________________________________________
AliMAG::AliMAG()
{
  //
  // Default constructor for L3 magnet
  //
}
 
//_____________________________________________________________________________
AliMAG::AliMAG(const char *name, const char *title)
  : AliModule(name,title)
{
  //
  // Standard constructor for L3 magnet
  //
  //Begin_Html
  /*
    <img src="picts/aliMAG.gif">
  */
  //End_Html
  
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}

//_____________________________________________________________________________
void AliMAG::CreateGeometry()
{
  //
  // Create geometry for L3 magnet
  //
  //Begin_Html
  /*
    <img src="picts/mag.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/tree_mag.gif">
  */
  //End_Html
  
  AliMC* pMC = AliMC::GetMC();
  
  Int_t *idtmed = gAlice->Idtmed();
  
  Float_t dpar[13];
  Int_t idrotm[399];
  Float_t par[10];
  
  // ANGLE POLAIRE MAXIMUM 
  
  //  Define Mother 
  
  par[0] = 22.5;
  par[1] = 360.;
  par[2] = 8.;
  par[3] = 2.;
  par[4] = -600.;
  par[5] = 580.;
  par[6] = 790.;
  par[7] = 600.;
  par[8] = 580.;
  par[9] = 790.;
  pMC->Gsvolu("L3MO", "PGON", idtmed[334], par, 10);
  pMC->Gspos("L3MO", 1, "ALIC", 0., -30., 0., 0, "ONLY");
  
  // Define coils 
  
  par[5] = 585.;
  par[6] = 690.;
  par[8] = 585.;
  par[9] = 690.;
  pMC->Gsvolu("L3CO", "PGON", idtmed[328], par, 10);
  pMC->Gspos("L3CO", 1, "L3MO", 0., 0., 0., 0, "ONLY");
  
  par[5] = 580.;
  par[6] = 585.;
  par[8] = 580.;
  par[9] = 585.;
  pMC->Gsvolu("L3C1", "PGON", idtmed[308], par, 10);
  pMC->Gspos("L3C1", 1, "L3MO", 0., 0., 0., 0, "ONLY");
  
  // Define yoke 
  
  par[5] = 690.;
  par[6] = 790.;
  par[8] = 690.;
  par[9] = 790.;
  pMC->Gsvolu("L3YO", "PGON", idtmed[329], par, 10);
  pMC->Gspos("L3YO", 1, "L3MO", 0., 0., 0., 0, "ONLY");
  
  // Define the return yoke of L3 (DOOR) 
  
  par[4] = 600.;
  par[5] = 0.;
  par[6] = 790.;
  par[7] = 700.;
  par[8] = par[5];
  par[9] = par[6];
  pMC->Gsvolu("L3DO", "PGON", idtmed[334], par, 10);
  
  par[4] = 610.;
  par[5] = 0.;
  par[6] = 790.;
  par[7] = 700.;
  par[8] = par[5];
  par[9] = par[6];
  pMC->Gsvolu("L3FR", "PGON", idtmed[329], par, 10);
  
  // INNER LAYER 
  
  par[4] = 600.;
  par[7] = 610.;
  pMC->Gsvolu("L3IR", "PGON", idtmed[309], par, 10);
  
  //     DOOR OPENING 
  
  dpar[0] = 22.5;
  dpar[1] = 360.;
  dpar[2] = 8.;
  dpar[3] = 3.;
  dpar[4] = 610.;
  dpar[5] = 0.;
  dpar[6] = 163.5;
  dpar[7] = 670.;
  dpar[8] = dpar[5];
  dpar[9] = dpar[6];
  dpar[10] = 700.;
  dpar[11] = dpar[5];
  dpar[12] = dpar[6] + 50.;
  pMC->Gsvolu("L3O1", "PGON", idtmed[314], dpar, 13);
  par[4] = 600.;
  par[5] = 0.;
  par[6] = 163.5;
  par[7] = 610.;
  par[8] = 0.;
  par[9] = 163.5;
  pMC->Gsvolu("L3O2", "PGON", idtmed[314], par, 10);
  
  //     THE DOOR OPENING HAS TO BE PLACED WITH 'MANY' SINCE THE REGION 
  //     WILL CONTAIN A MUON CHAMBER, BEAM PIPE AND BEAM SHIELD 
  //     PLACED WITH 'ONLY'. 
  
  pMC->Gspos("L3O1", 1, "L3FR", 0., 30., 0., 0, "MANY");
  pMC->Gspos("L3O2", 1, "L3IR", 0., 30., 0., 0, "MANY");
  
  pMC->Gspos("L3FR", 1, "L3DO", 0., 0., 0., 0, "MANY");
  pMC->Gspos("L3IR", 1, "L3DO", 0., 0., 0., 0, "MANY");
  
  pMC->Gspos("L3DO", 1, "ALIC", 0., -30., 0., 0, "MANY");
  AliMatrix(idrotm[300], 90., 0., 90., 90., 180., 0.);
  pMC->Gspos("L3DO", 2, "ALIC", 0., -30., 0., idrotm[300], "MANY");
}

//_____________________________________________________________________________
void AliMAG::CreateMaterials()
{
  //
  // Create materials for L3 magnet
  //
  
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Float_t epsil, stmin, deemax, tmaxfd, stemax;
  // --- Define the various materials for GEANT --- 
  
  //     Aluminum 
  AliMaterial(9, "Al$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "Al$", 26.98, 13., 2.7, 8.9, 37.2);
  
  //     Iron 
  AliMaterial(10, "Fe$", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "Fe$", 55.85, 26., 7.87, 1.76, 17.1);
  
  //     Air 
  AliMaterial(15, "Air$", 14.61, 7.3, .001205, 30420., 67500);
  AliMaterial(35, "Air$", 14.61, 7.3, .001205, 30420., 67500);
  
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  
  //    IRON 
  
  AliMedium(310, "FE_C0             ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(330, "FE_C1             ", 30, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //     ALUMINUM 

  AliMedium(309, "ALU_C0            ",  9, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(329, "ALU_C1            ", 29, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //     AIR 
  
  AliMedium(315, "AIR_C0            ", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(335, "AIR_C1            ", 35, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
void AliMAG::DrawModule()
{
  //
  // Draw a shaded view of the L3 magnet
  //
}

//_____________________________________________________________________________
void AliMAG::Init()
{
  //
  // Initialise L3 magnet after it has been built
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" MAG_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the MAG initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

