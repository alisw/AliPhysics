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

  Float_t dASUC[3];
  Int_t *idtmed = fIdtmed->GetArray()+1;
  int idSC = idtmed[0];
  //
  dASUC[0]=50;
  dASUC[1]=50;
  dASUC[2]=50;
  //  TString tmp(GetTitle());
  gMC->Gsvolu(GetName(),"BOX",idSC, dASUC,3); // WSUC - Wsu Cosmic Ray SetUp
}
 
//_____________________________________________________________________________
void AliEMCALWsuCosmicRaySetUp::CreateMaterials()
{
// Create materials and media
  Int_t   isxfld = 0;
  Float_t sxmgmx = 0.;
  
  // AIR
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  //  Float_t dAir1 = 1.20479E-10;
  //
  AliMixture(1,"Air     $",aAir,zAir,dAir,4,wAir);
  //
  AliMedium(1,"Air     $",1,0,isxfld,sxmgmx,10,-1,-0.1,0.1 ,-10);
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
 

