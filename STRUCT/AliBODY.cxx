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
//  Alice external volume                                                    //
//  This class contains the description of the Alice external volume         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliBODYClass.gif">
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
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TGeoGlobalMagField.h>
#include <TVirtualMC.h>
#include <TArrayI.h>

#include "AliBODY.h"
#include "AliMagF.h"
#include "AliRun.h"

ClassImp(AliBODY)
 
//_____________________________________________________________________________
AliBODY::AliBODY()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliBODY::AliBODY(const char *name, const char *title)
       : AliModule(name,title)
{
  //
  // Standard constructor of the Alice external volume
  //
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliBODY::CreateGeometry()
{
  //
  // Create the geometry of the Alice external body
  //
  //Begin_Html
  /*
    <img src="picts/AliBODYTree.gif">
  */
  //End_Html
  //
  // If the ZDC is present we have an asymmetric box
  // made by a four sides polygone
  //  
  //Begin_Html
  /*
    <img src="picts/AliBODYLarge.gif">
  */
  //End_Html
  //
  // If the ZDC is not present make just a BOX
  //
  //Begin_Html
  /*
    <img src="picts/AliBODYSmall.gif">
  */
  //End_Html

  Float_t dALIC[10];
  Int_t *idtmed = fIdtmed->GetArray()+1;
  //
  if(gAlice->GetModule("ZDC")) {
    //
    dALIC[0]=2500.;
    dALIC[1]=2500.;
    dALIC[2]=15000.;
    gMC->Gsvolu("ALIC","BOX ",idtmed[1],dALIC,3);
  } else if ( gAlice->GetModule("ACORDE")) {
    //
    // If the Cosmic Ray Trigger  is present we need a large box
    // 
    //
    dALIC[0]=13000.;
    dALIC[1]=5000.;
    dALIC[2]=13000.;
    gMC->Gsvolu("ALIC","BOX ",idtmed[1],dALIC,3);
      
  } else {
    //
    // If the ZDC and ACORDE are not present make just a BOX
    //
    dALIC[0]=2000;
    dALIC[1]=2000;
    dALIC[2]=3000;
    gMC->Gsvolu("ALIC","BOX ",idtmed[1],dALIC,3);
  }
}
 
//_____________________________________________________________________________
void AliBODY::CreateMaterials()
{
// Create materials and media
  Int_t isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  
  // AIR

  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3 * 960./1014.;
  Float_t dAir1 = 1.20479E-10;
  //
  AliMixture(1,"Vacuum  $",aAir,zAir,dAir1,4,wAir);
  AliMixture(2,"Air     $",aAir,zAir,dAir,4,wAir);
  AliMaterial(3,"Be      $", 9.01,4 ,1.848   ,35.30,36.70);
  //
  AliMedium(1,"Vacuum  $",1,0,isxfld,sxmgmx,10,1,0.1,0.1,10);
  AliMedium(2,"Air     $",2,0,isxfld,sxmgmx,10,-1,-0.1,0.1 ,-10);
  AliMedium(3,"Be pipe $",3,0,isxfld,sxmgmx,10,0.1,0.1,0.01,0.01);
}
 
