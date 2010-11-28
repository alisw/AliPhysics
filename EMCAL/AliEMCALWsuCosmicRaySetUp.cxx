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
//  Wsu Cosmic Ray SetUp                                                     //
//  This class contains the description of the  Wsu Cosmic Ray SetUp         //
//  external volume                                                          //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliEMCALWsuCosmicRaySetUpClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:pavlinov@physics.wayne.edu">Aleksei Pavlino, WSU</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TVirtualMC.h>

#include "AliEMCALWsuCosmicRaySetUp.h"
//#include "AliMagF.h"
#include "AliRun.h"

ClassImp(AliEMCALWsuCosmicRaySetUp)
 
//_____________________________________________________________________________
AliEMCALWsuCosmicRaySetUp::AliEMCALWsuCosmicRaySetUp()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliEMCALWsuCosmicRaySetUp::AliEMCALWsuCosmicRaySetUp(const char *name, const char *title)
       : AliModule(name,title)
{
  //
  // Standard constructor of the  Wsu Cosmic Ray SetUp external volume
  //
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliEMCALWsuCosmicRaySetUp::CreateGeometry()
{
  //
  // Create the geometry of the Alice external body
  //
  //Begin_Html
  /*
    <img src="picts/AliEMCALWsuCosmicRaySetUpTree.gif">
  */
  //End_Html

  // Master Volume
  fMasterVolume[0] = fMasterVolume[1] = 25.0;
  fMasterVolume[2] = 300.;

  Int_t *idtmed = fIdtmed->GetArray()+1;
  int idAir = idtmed[0];
  gMC->Gsvolu(GetName(),"BOX",idAir, fMasterVolume,3); // Master volume
  //
  // Sc counters
  //
  Float_t sc[3]; // tube
  sc[0] = 0.0;
  sc[1] = 5.0;
  sc[2] = 0.5;
  Float_t zsc[3]={10.,110., 310.};
  int idSC = idtmed[1];
  gMC->Gsvolu("SCOU","TUBE",idSC, sc,3); // Master volume
  Int_t idRot=0; // no rotation
  for(Int_t i=0; i<3; i++) {
    Float_t zpos = zsc[i] - fMasterVolume[2];
    gMC->Gspos("SCOU", i+1, "WSUC", 0.0, 0.0, zpos, idRot, "ONLY"); 
  }
}
 
//_____________________________________________________________________________
void AliEMCALWsuCosmicRaySetUp::CreateMaterials()
{
// Create materials and media
  Int_t   isxfld = 0;
  Float_t sxmgmx = 0., deemax = 0.1;  
  // AIR
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  AliMixture(1,"Air     $",aAir,zAir,dAir,4,wAir);

  // --- The polysterene scintillator (CH) ---
  Float_t aP[2] = {12.011, 1.00794} ;
  Float_t zP[2] = {6.0, 1.0} ;
  Float_t wP[2] = {1.0, 1.0} ;
  Float_t dP = 1.032 ;
  AliMixture(2, "Polystyrene$", aP, zP, dP, -2, wP) ;

  //
  AliMedium(1,"Air     $",1,0,isxfld,sxmgmx,10,-1,-0.1,0.1 ,-10);
  AliMedium(2, "Scintillator$", 2, 1,
            isxfld, sxmgmx, 10.0, 0.001, deemax, 0.001, 0.001, 0, 0) ;
  
}
 
//_____________________________________________________________________________
void AliEMCALWsuCosmicRaySetUp::DrawWSUC(float cxy) const
{
  //
  // Draw a view of the Wsu Cosmic Ray SetUp 
  //
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set WSUC mother visible
  gMC->Gsatt("WSUC","SEEN",1);
  //
  // Set the volumes visible
  //
  gMC->Gdopt("hide","off");

  gMC->Gdraw("WSUC", 40, 30, 0, 10, 9, cxy, cxy);
  gMC->Gdhead(1111, "WSU Cosmic Ray Setup ");

  gMC->Gdman(18, 4, "MAN");
}
 

