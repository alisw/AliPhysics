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
//  Experimental Hall                                                        //
//  This class contains the description of the experimental hall             //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliHALLClass.gif">
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

#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliHALL.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliLog.h"
 
ClassImp(AliHALL)
 
//_____________________________________________________________________________
AliHALL::AliHALL()
{
  //
  // Default constructor for the experimental Hall
  //
}
 
//_____________________________________________________________________________
AliHALL::AliHALL(const char *name, const char *title)
       : AliModule(name,title)
{
  //
  // Standard constructor for the experimental Hall
  //
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliHALL::CreateGeometry()
{
  //
  // Create the geometry of the exprimental hall
  //
  //Begin_Html
  /*
    <img src="picts/AliHALLTree.gif">
  */
  //End_Html
  //
  // If ZDC is not present the experimental hall includes a short
  // section of the accelerator tunnel
  //
  //Begin_Html
  /*
    <img src="picts/AliHALLSmall.gif">
  */
  //End_Html
  //
  // If ZDC is present the experimental hall includes the accelerator
  // tunnel beyond the ZDC
  //
  //Begin_Html
  /*
    <img src="picts/AliHALLLarge.gif">
  */
  //End_Html

  
  Float_t r2;
  Float_t phid, phim, pbox[3], h, r, tspar[5];
  Float_t w1, dh, am, bm, dl,cm, hm, dr, dx, xl;
  Int_t idrotm[1999];
  Float_t trdpar[4], trapar[11], hullen;
  Float_t phi;
  
  Int_t *idtmed = fIdtmed->GetArray()-1899;
  
  //     RB24/26 TUNNEL FLOOR 
  
  r   = 220.;
  h   = 140.;
  phi = TMath::ACos(h / r);
  xl  = r * TMath::Sin(phi);
  dr  = 1600.;
  dh  = dr * TMath::Cos(phi);
  dl  = dr * TMath::Sin(phi);
  if (gAlice->GetModule("ZDC") == 0) {
    
    //     No ZDC 
    hullen = 370.;
  } else {
    
    //     ZDC is present 
    hullen = 6520.;
  }
  trdpar[0] = xl + dl;
  trdpar[1] = xl;
  trdpar[2] = hullen;
  trdpar[3] = dh / 2.;
  AliMatrix(idrotm[1900], 90., 0., 0., 0., 90., 90.);
  AliMatrix(idrotm[1901], 270., 0., 90., 90., 0., 0.);
  gMC->Gsvolu("HUFL", "TRD1", idtmed[1956], trdpar, 4);
  r2 = hullen + 1900.;
  gMC->Gspos("HUFL", 1, "ALIC", 70.,-100-trdpar[3] , -r2, idrotm[1900], "ONLY");
  
  //     RB24/26 wall 
  
  phid     = phi * 57.296;
  tspar[0] = r;
  tspar[1] = r + dr;
  tspar[2] = hullen;
  tspar[3] = phid - 90.;
  tspar[4] = 270. - phid;
  gMC->Gsvolu("HUWA", "TUBS", idtmed[1956], tspar, 5);
  gMC->Gspos("HUWA", 1, "ALIC", 70., 40., -1900 - hullen , 0, "ONLY");
  
  //     END WALL 
  //gMC->Gsvolu("HEW1", "BOX ", idtmed[1956], pbox, 0);
  //pbox[0] = 600.;
  //pbox[1] = 418.;
  //pbox[2] = 60.;
  //gMC->Gsposp("HEW1", 1, "ALIC", 0., -pbox[1]-60., -1960, 0, "ONLY", pbox, 3);
  //pbox[1] = 822.;
  //gMC->Gsposp("HEW1", 2, "ALIC", 0.,  pbox[1]+60., -1960, 0, "ONLY", pbox, 3);
  //pbox[0] = 270.;
  //pbox[1] =  60.;
  //gMC->Gsposp("HEW1", 3, "ALIC",  pbox[0]+60.,  0. , -1960, 0, "ONLY", pbox, 3);
  //gMC->Gsposp("HEW1", 4, "ALIC", -pbox[0]-60.,  0. , -1960, 0, "ONLY", pbox, 3);

  //     hall floor 
  
  phid      = 16.197;
  trdpar[0] = 700.;
  trdpar[1] = TMath::Tan(phid * kDegrad) * 190. + 700.;
  trdpar[2] = 550.;
  trdpar[3] = 95.;
  gMC->Gsvolu("HHF1", "TRD1", idtmed[1956], trdpar, 4);
  gMC->Gspos("HHF1", 1, "ALIC", 0., -801., 1350., idrotm[1900], "ONLY");
  gMC->Gspos("HHF1", 2, "ALIC", 0., -801.,-1350., idrotm[1900], "ONLY");
  
  //     hall side walls 
  
  trapar[0] = 550.;
  trapar[1] = 0.;
  trapar[2] = 0.;
  trapar[3] = 1273.78/2;
  trapar[4] = 235.;
  trapar[5] = 50.;
  trapar[6] = TMath::ATan((trapar[4] - trapar[5]) / 2. / trapar[3]) * kRaddeg;
  trapar[7] = trapar[3];
  trapar[8] = trapar[4];
  trapar[9] = trapar[5];
  trapar[10] = trapar[6];
  dx = trapar[4] * 1.5 + 700. - trapar[5] * .5;
  gMC->Gsvolu("HHW1", "TRAP", idtmed[1956], trapar, 11);
  gMC->Gspos("HHW1", 1, "ALIC", dx, -896+trapar[3],  1350., 0, "ONLY");
  gMC->Gspos("HHW1", 2, "ALIC",-dx, -896+trapar[3],  1350., idrotm[1901], "ONLY");
  gMC->Gspos("HHW1", 3, "ALIC", dx, -896+trapar[3], -1350., 0, "ONLY");
  gMC->Gspos("HHW1", 4, "ALIC",-dx, -896+trapar[3], -1350., idrotm[1901], "ONLY");
  pbox[0] = 50.;
  pbox[1] = (500. - (trapar[3] * 2. - 896.)) / 2.;
  pbox[2] = 1900.;
  gMC->Gsvolu("HBW1", "BOX ", idtmed[1956], pbox, 3);
  gMC->Gspos("HBW1", 1, "ALIC",  1120., 500-pbox[1], 0., 0, "ONLY");
  gMC->Gspos("HBW1", 2, "ALIC", -1120., 500-pbox[1], 0., 0, "ONLY");
  
  //     slanted wall close to L3 magnet 
  
  phim = 45.;
  hm   = 790.;
  //rm   = hm / TMath::Cos(phim / 2. * kDegrad);
  am   = hm * TMath::Tan(phim / 2. * kDegrad);
  bm   = (hm + 76.) / hm * am;
  cm   = bm * 2. / TMath::Sqrt(2.);
  trapar[0] = 800.;
  trapar[1] = 0.;
  trapar[2] = 0.;
  trapar[3] = (1273.78 - cm) / 2.;
  trapar[4] = 235. - cm * TMath::Tan(phid * kDegrad) / 2.;
  trapar[5] = 50.;
  trapar[6] = TMath::ATan((trapar[4] - trapar[5]) / 2. / trapar[3]) * kRaddeg;
  trapar[7] = trapar[3];
  trapar[8] = trapar[4];
  trapar[9] = trapar[5];
  trapar[10] = trapar[6];
  w1 = trapar[4];
  dx = cm*TMath::Tan(phid * kDegrad) + 700. + trapar[4] * 1.5 - trapar[5] * .5;
  gMC->Gsvolu("HHW2", "TRAP", idtmed[1956], trapar, 11);
  r2 = cm - 896. + trapar[3];
  gMC->Gspos("HHW2", 1, "ALIC", dx, r2, 0., 0, "ONLY");
  gMC->Gspos("HHW2", 2, "ALIC",-dx, r2, 0., idrotm[1901], "ONLY");
  trapar[3]  = cm / 2.;
  trapar[4]  = w1 + cm / 2.;
  trapar[5]  = w1;
  trapar[6]  = TMath::ATan(.5) * kRaddeg;
  trapar[7]  = trapar[3];
  trapar[8]  = trapar[4];
  trapar[9]  = trapar[5];
  trapar[10] = trapar[6];
  dx = 1170. - trapar[4] * .5 - trapar[5] * .5;
  gMC->Gsvolu("HHW3", "TRAP", idtmed[1956], trapar, 11);
  r2 = trapar[3] - 896.;
  gMC->Gspos("HHW3", 1, "ALIC", dx, r2, 0., 0, "ONLY");
  gMC->Gspos("HHW3", 2, "ALIC",-dx, r2, 0., idrotm[1901], "ONLY");
  

  tspar[0] = 1070.;
  tspar[1] = 1170.;
  tspar[2] = 1900.;
  tspar[3] = 0.;
  tspar[4] = 180.;
  gMC->Gsvolu("HHC1", "TUBS", idtmed[1956], tspar, 5);
  gMC->Gspos("HHC1", 1, "ALIC", 0., 500., 0., 0, "ONLY");
  trdpar[0] = 1170 - trapar[4] * 2.;
  trdpar[1] = trdpar[0] + TMath::Tan(phim * kDegrad) * 76.;
  trdpar[2] = 800.;
  trdpar[3] = 38.;
  gMC->Gsvolu("HHF2", "TRD1", idtmed[1956], trdpar, 4);
  gMC->Gspos("HHF2", 1, "ALIC", 0., -858., 0., idrotm[1900], "ONLY");
  
  //     pillars for working platform 
  
  pbox[0] = 40.;
  pbox[1] = 96.;
  pbox[2] = 550.;
  gMC->Gsvolu("HPIL", "BOX ", idtmed[1956], pbox, 3);
  gMC->Gspos("HPIL", 1, "ALIC", 165.,-706+pbox[1] , -1350., 0, "ONLY");
  gMC->Gspos("HPIL", 2, "ALIC",-165.,-706+pbox[1] , -1350., 0, "ONLY");
  
  //     simple concrete beam shield 
  
  Float_t ppgon[10];
  ppgon[0] =    45.;
  ppgon[1] =   360.;
  ppgon[2] =     4.;
  ppgon[3] =     2.;
  ppgon[7] =  1800.;
  ppgon[8] =   150.;
  ppgon[9] =   250.;
  ppgon[4] =   720.;
  ppgon[5] =   150.;
  ppgon[6] =   250.;

  gMC->Gsvolu("HMBS", "PGON", idtmed[1956], ppgon, 10);
  gMC->Gspos("HMBS", 1, "ALIC", 0., 70., 0., 0, "ONLY");
}

//_____________________________________________________________________________
void AliHALL::CreateMaterials()
{
  //
  // Create materials for the experimental hall
  //
    
  Int_t   isxfld = 0;
  Float_t sxmgmx = gAlice->Field()->Max();
  
  Float_t aconc[10] = { 1.,12.01,15.994,22.99,24.305,26.98,28.086,39.1,40.08,55.85 };
  Float_t zconc[10] = { 1.,6.,8.,11.,12.,13.,14.,19.,20.,26. };
  Float_t wconc[10] = { .01,.001,.529107,.016,.002,.033872,.337021,.013,.044,.014 };
  
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;

  Float_t epsil, stmin, deemax, tmaxfd, stemax;

  //     FOR CONCRETE 
  
  AliMaterial(10, "IRON0$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON1$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON2$     ", 55.85, 26., 7.87, 1.76, 17.1); 

  AliMixture(15, "AIR0$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(35, "AIR1$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(55, "AIR2$      ", aAir, zAir, dAir, 4, wAir);

  AliMixture(17, "CONCRETE0$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(37, "CONCRETE1$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(57, "CONCRETE2$", aconc, zconc, 2.35, 10, wconc);
  
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
  
  //     IRON 
  
  AliMedium(10, "FE_C0             ", 10, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1             ", 30, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(50, "FE_C2             ", 50, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Air 
  
  AliMedium(15, "AIR_C0           ", 15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_C1           ", 35, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(55, "AIR_C2           ", 55, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Concrete 
  
  AliMedium(17, "CC_C0            ", 17, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(37, "CC_C1            ", 37, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(57, "CC_C2            ", 57, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
void AliHALL::Init()
{
  //
  // Initialise the HALL after it has been built
  //
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" HALL_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the HALL initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

//_____________________________________________________________________________
void AliHALL::DrawModule() const
{
  //
  // Draw a shaded view of Experimental Hall
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("HUFL","seen",1);
  gMC->Gsatt("HUWA","seen",1);
  gMC->Gsatt("HUP2","seen",1);
  gMC->Gsatt("HEW1","seen",1);
  gMC->Gsatt("HHF1","seen",1);
  gMC->Gsatt("HHW1","seen",1);
  gMC->Gsatt("HBW1","seen",1);
  gMC->Gsatt("HHW2","seen",1);
  gMC->Gsatt("HHW3","seen",1);
  gMC->Gsatt("HHC1","seen",1);
  gMC->Gsatt("HHF2","seen",1);
  gMC->Gsatt("HPIL","seen",1);
  gMC->Gsatt("HMBS","seen",1);
  gMC->Gsatt("HBBS","seen",1);
  gMC->Gsatt("HPBS","seen",1);
  gMC->Gsatt("HXFI","seen",1);
  gMC->Gsatt("HXII","seen",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  if (gAlice->GetModule("ZDC") == 0) {
    //
    // ZDC is not present
    //
    gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
    gMC->DefaultRange();
    gMC->Gdraw("alic", 40, 30, 0, 12, 7.5, .005, .005);
  } else {
    //
    // ZDC is present
    //
    gMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 16000);
    gMC->DefaultRange();
    gMC->Gdraw("alic", 40, 30, 0, 17.5, 10, .0019, .0019);
  }
  gMC->Gdhead(1111, "Experimental Hall");
  gMC->Gdman(18, 2, "MAN");
}
 
