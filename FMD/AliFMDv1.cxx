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
Revision 1.11  2000/12/22 11:17:39  hristov
New FMD code from Alla + code cleaning

Revision 1.10  2000/10/02 21:28:07  fca
Removal of useless dependecies via forward declarations

Revision 1.9  2000/05/10 21:56:07  fca
Avoid clashes with ITS and add supports

Revision 1.7  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Forward Multiplicity Detector version 1                                  //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliFMDv1Class.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:Valeri.Kondratiev@cern.ch">Valeri Kondratiev</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"
#include "AliFMDv1.h"
#include "AliMC.h"
#include "AliMagF.h"
#include <stdlib.h>
 
ClassImp(AliFMDv1)
 
//_____________________________________________________________________________
AliFMDv1::AliFMDv1()
{
  //
  // Default constructor for FMD version 1
  //
}
 
//_____________________________________________________________________________
AliFMDv1::AliFMDv1(const char *name, const char *title)
  : AliFMD(name,title)
{
  //
  // Standard constructor for FMD version 1
  //
  AliModule *start = gAlice->GetModule("START");
  if(start) {
    Error("ctor","This version of FMD is incompatible with START\n");
    exit(1);
  }
}
 
//___________________________________________
void AliFMDv1::CreateGeometry()
{
  //
  // Creation of the geometry of the FMD version 1
  //
  //Begin_Html
  /*
    <img src="picts/AliFMDv1Tree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliFMDv1.gif">
  */
  //End_Html

  
  Float_t rout, z;
  Float_t par[3], rin;

  Int_t irotm[100];
  
  TArrayI &idtmed = *fIdtmed;    
  
  // ******************************************************** 
  //       DEFINE DISK#3  OF FMD 
  // ******************************************************** 
  
  gMC->Gsvolu("BR3_", "BOX ", idtmed[4], par, 0);
  gMC->Gsvolu("CB3_", "BOX ", idtmed[5], par, 0);
  gMC->Gsvolu("BR5_", "BOX ", idtmed[4], par, 0);
  gMC->Gsvolu("CB5_", "BOX ", idtmed[5], par, 0);

  //       Define parameters for the right disk 
  
  rin  = 4.5;
  rout = 10.5;
  z    = 77.;

  //       Ring #1 

  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1R3", "TUBE", idtmed[1], par, 3);
  gMC->Gspos("R1R3", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");

  //       Ring #2 

  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2R3", "TUBE", idtmed[2], par, 3);
  gMC->Gspos("R2R3", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");

  //       Ring #3 

  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3R3", "TUBE", idtmed[3], par, 3);
  gMC->Gspos("R3R3", 1, "ALIC", 0., 0., z + 1.525, 0, "ONLY");

  //       Bracket #1

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[11], 90., 0., 161.2, 90., 71.2, 90.);
  gMC->Gsposp("BR3_",1,"ALIC", 0., 25.25, 85.0, irotm[11],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",1,"ALIC", 0., 25.41, 84.53,irotm[11],"ONLY",par,3); 


  //       Bracket #2

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[12], 90., 180., 161.2, 270., 71.2, 270.);
  gMC->Gsposp("BR3_",2,"ALIC", 0.,-25.25, 85.0, irotm[12],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",2,"ALIC", 0.,-25.41, 84.53,irotm[12],"ONLY",par,3); 


  //       Bracket #3

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[13], 90., 270., 161.2, 0., 71.2, 0.);
  gMC->Gsposp("BR3_",3,"ALIC", 25.25, 0., 85.0, irotm[13],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",3,"ALIC", 25.41, 0., 84.53,irotm[13],"ONLY",par,3); 


  //       Bracket #4

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[14], 90., 90., 161.2, 180., 71.2, 180.);
  gMC->Gsposp("BR3_",4,"ALIC", -25.25, 0., 85.0, irotm[14],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",4,"ALIC", -25.41, 0., 84.53,irotm[14],"ONLY",par,3); 

  //       Right support ring 

  par[0] = 39.;
  par[1] = 41.;
  par[2] = .5;
  gMC->Gsvolu("R1SP", "TUBE", idtmed[4], par, 3);
  gMC->Gspos("R1SP", 1, "ALIC", 0., 0., 89.5, 0, "ONLY");

  //       Define parameters for the left disk
  
  rin  = 4.5;
  rout = 10.5;
  z    = -77.;

  //       Ring #1 

  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L3", "TUBE", idtmed[1], par, 3);
  gMC->Gspos("R1L3", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #2 

  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L3", "TUBE", idtmed[2], par, 3);
  gMC->Gspos("R2L3", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #3 

  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L3", "TUBE", idtmed[3], par, 3);
  gMC->Gspos("R3L3", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");

  //       Bracket #1

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[61], 90., 0., 180.-161.2, 90., 180.-71.2, 90.);
  gMC->Gsposp("BR3_",5,"ALIC", 0., 25.25, -85.0, irotm[61],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",5,"ALIC", 0., 25.41, -84.53,irotm[61],"ONLY",par,3); 

  gMC->Matrix(irotm[91], 90., 0., 90.-10.6, 90., 180.-10.6, 90.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",1,"ALIC", 0., 55.0, -170.0,irotm[91],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",1,"ALIC", 0., 55.5, -169.5,irotm[91],"ONLY",par,3); 

  //       Bracket #2

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[62], 90., 180., 180.-161.2, 270., 180.-71.2, 270.);
  gMC->Gsposp("BR3_",6,"ALIC", 0.,-25.25, -85.0,irotm[62],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",6,"ALIC", 0.,-25.41, -84.53,irotm[62],"ONLY",par,3); 

  gMC->Matrix(irotm[92], 90., 180., 90.-10.6, 270., 180.-10.6, 270.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",2,"ALIC", 0., -55.0, -170.0,irotm[92],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",2,"ALIC", 0., -55.5, -169.5,irotm[92],"ONLY",par,3); 

  //       Bracket #3

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[63], 90., 270., 180.-161.2, 0., 180.-71.2, 0.);
  gMC->Gsposp("BR3_",7,"ALIC", 25.25, 0., -85.0, irotm[63],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",7,"ALIC", 25.41, 0., -84.53,irotm[63],"ONLY",par,3); 

  gMC->Matrix(irotm[93], 90., 270., 90.-10.6, 0., 180.-10.6, 0.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",3,"ALIC", 55., 0., -170.,irotm[93],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",3,"ALIC", 55.5, 0., -169.5,irotm[93],"ONLY",par,3); 

  //       Bracket #4

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 15.35;
  gMC->Matrix(irotm[64], 90., 90., 180.-161.2, 180., 180.-71.2, 180.);
  gMC->Gsposp("BR3_",8,"ALIC", -25.25, 0., -85., irotm[64],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB3_",8,"ALIC", -25.41, 0., -84.53,irotm[64],"ONLY",par,3); 

  gMC->Matrix(irotm[94], 90., 90., 90.-10.6, 180., 180.-10.6, 180.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",4,"ALIC", -55., 0., -170.,irotm[94],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",4,"ALIC", -55.5, 0., -169.5,irotm[94],"ONLY",par,3); 

  //       Central support ring 

  par[0] = 39.;
  par[1] = 41.;
  par[2] = .5;
  gMC->Gsvolu("R2SP", "TUBE", idtmed[4], par, 3);
  gMC->Gspos("R2SP", 1, "ALIC", 0., 0., -89.5, 0, "ONLY");

  //        Left support ring 

  par[0] = 69.;
  par[1] = 71.;
  par[2] = .5;
  gMC->Gsvolu("R3SP", "TUBE", idtmed[4], par, 3);
  gMC->Gspos("R3SP", 1, "ALIC", 0., 0., -249.5, 0, "ONLY");

  // ******************************************************** 
  //       DEFINE  DISK#2  OF FMD 
  // ******************************************************** 
  
  gMC->Gsvolu("BR2_", "BOX ", idtmed[4], par, 0);
  gMC->Gsvolu("CB2_", "BOX ", idtmed[5], par, 0);

  //       Define parameters for right disk #2
  
  rin  = 7.7;
  rout = 13.7;
  z    = 64.4;

  //       Ring #1 

  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1R2", "TUBE", idtmed[1], par, 3);
  gMC->Gspos("R1R2", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");

  //       Ring #2 

  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2R2", "TUBE", idtmed[2], par, 3);
  gMC->Gspos("R2R2", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");

  //       Ring #3 

  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3R2", "TUBE", idtmed[3], par, 3);
  gMC->Gspos("R3R2", 1, "ALIC", 0., 0., z + 1.525, 0, "ONLY");

  //       Bracket #1

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[21], 90., 30., 139.3, 120., 49.3, 120.);
  gMC->Gsposp("BR2_",1,"ALIC", -13.4, 23.3, 78.7, irotm[21],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",1,"ALIC", -13.56, 23.6, 78.4, irotm[21],"ONLY",par,3); 

  //       Bracket #2

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[22], 90., 210., 139.3, 300., 49.3, 300.);
  gMC->Gsposp("BR2_",2,"ALIC", 13.4,-23.3, 78.7, irotm[22],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",2,"ALIC", 13.5,-23.6, 78.4,irotm[22],"ONLY",par,3); 


  //       Bracket #3

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[23], 90., 300., 139.3, 30., 49.3, 30.);
  gMC->Gsposp("BR2_",3,"ALIC", 23.3, 13.4, 78.7, irotm[23],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",3,"ALIC", 23.6, 13.56, 78.4,irotm[23],"ONLY",par,3); 

  //       Bracket #4

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[24], 90., 120., 139.3, 210., 49.3, 210.);
  gMC->Gsposp("BR2_",4,"ALIC", -23.3, -13.4, 78.7,irotm[24],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",4,"ALIC", -23.6, -13.56, 78.4,irotm[24],"ONLY",par,3); 

  //       Define parameters for left disk
  
  rin  = 7.7;
  rout = 13.7;
  z    = -64.4;

  //       Ring #1 

  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L2", "TUBE", idtmed[1], par, 3);
  gMC->Gspos("R1L2", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #2 

  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L2", "TUBE", idtmed[2], par, 3);
  gMC->Gspos("R2L2", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #3 

  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L2", "TUBE", idtmed[3], par, 3);
  gMC->Gspos("R3L2", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");

  //       Bracket #1

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[51], 90., 30., 180.-139.3, 120., 180.-49.3, 120.);
  gMC->Gsposp("BR2_",5,"ALIC", -13.4, 23.3, -78.7, irotm[51],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",5,"ALIC", -13.56, 23.6, -78.4,irotm[51],"ONLY",par,3); 

  gMC->Matrix(irotm[81], 90., 30., 90.-10.6, 120., 180.-10.6, 120.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",5,"ALIC", -27.5, 47.6, -170.0,irotm[81],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",5,"ALIC", -27.85, 48., -169.5,irotm[81],"ONLY",par,3); 

  //       Bracket #2

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[52], 90., 210., 180.-139.3, 300., 180.-49.3, 300.);
  gMC->Gsposp("BR2_",6,"ALIC", 13.4, -23.3, -78.7,irotm[52],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",6,"ALIC", 13.56, -23.6, -78.4,irotm[52],"ONLY",par,3); 

  gMC->Matrix(irotm[82], 90., 210., 90.-10.6, 300., 180.-10.6, 300.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",6,"ALIC", 27.5, -47.6, -170.0,irotm[82],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",6,"ALIC", 27.85, -48., -169.5,irotm[82],"ONLY",par,3); 

  //       Bracket #3

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[53], 90., 300., 180.-139.3, 30., 180.-49.3, 30.);
  gMC->Gsposp("BR2_",7,"ALIC", 23.3, 13.4, -78.7, irotm[53],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",7,"ALIC", 23.6, 13.56, -78.4,irotm[53],"ONLY",par,3); 

  gMC->Matrix(irotm[83], 90., 300., 90.-10.6, 30., 180.-10.6, 30.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",7,"ALIC", 47.6, 27.5, -170.,irotm[83],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",7,"ALIC", 48., 27.85, -169.5,irotm[83],"ONLY",par,3); 

  //       Bracket #4

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.35;
  gMC->Matrix(irotm[54], 90., 120., 180.-139.3, 210., 180.-49.3, 210.);
  gMC->Gsposp("BR2_",8,"ALIC", -23.3, -13.4, -78.7, irotm[54],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB2_",8,"ALIC", -23.6, -13.56, -78.4,irotm[54],"ONLY",par,3); 

  gMC->Matrix(irotm[84], 90., 120., 90.-10.6, 210., 180.-10.6, 210.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",8,"ALIC", -47.6, -27.5, -170.,irotm[84],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",8,"ALIC", -48., -27.85, -169.5,irotm[84],"ONLY",par,3); 

  // ******************************************************** 
  //       DEFINE  DISK#1  OF FMD 
  // ******************************************************** 
  
  gMC->Gsvolu("BR1_", "BOX ", idtmed[4], par, 0);
  gMC->Gsvolu("CB1_", "BOX ", idtmed[5], par, 0);
  
  //       Define parameters for right disk #1
  
  rin  = 12.4;
  rout = 18.4;
  z    = 59.4;

  //       Ring #1 

  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1R1", "TUBE", idtmed[1], par, 3);
  gMC->Gspos("R1R1", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");

  //       Ring #2 

  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2R1", "TUBE", idtmed[2], par, 3);
  gMC->Gspos("R2R1", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");

  //       Ring #3 

  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3R1", "TUBE", idtmed[3], par, 3);
  gMC->Gspos("R3R1", 1, "ALIC", 0., 0., z + 1.525, 0, "ONLY");

  //       Bracket #1

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[31], 90., 60., 128., 150., 38., 150.);
  gMC->Gsposp("BR1_",1,"ALIC", -25.3, 14.6, 76.2, irotm[31],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",1,"ALIC", -25.35, 14.8, 75.9, irotm[31],"ONLY",par,3); 

  //       Bracket #2

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[32], 90., 240., 128., 330., 38., 330.);
  gMC->Gsposp("BR1_",2,"ALIC", 25.3, -14.6, 76.2, irotm[32],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",2,"ALIC", 25.35, -14.8, 75.9,irotm[32],"ONLY",par,3); 

  //       Bracket #3

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[33], 90., 330., 128., 60., 38., 60.);
  gMC->Gsposp("BR1_",3,"ALIC", 14.6, 25.3, 76.2, irotm[33],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",3,"ALIC", 14.8, 25.35, 75.9,irotm[33],"ONLY",par,3); 

  //       Bracket #4

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[34], 90., 150., 128., 240., 38., 240.);
  gMC->Gsposp("BR1_",4,"ALIC", -14.6, -25.3, 76.2,irotm[34],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",4,"ALIC", -14.8, -25.35, 75.9,irotm[34],"ONLY",par,3); 

  //       Define parameters for left disk #1
  
  rin  = 12.4;
  rout = 18.4;
  z    = -59.4;

  //       Ring #1 

  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L1", "TUBE", idtmed[1], par, 3);
  gMC->Gspos("R1L1", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #2 

  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L1", "TUBE", idtmed[2], par, 3);
  gMC->Gspos("R2L1", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #3 

  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L1", "TUBE", idtmed[3], par, 3);
  gMC->Gspos("R3L1", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");

  //       Bracket #1

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[41], 90., 60., 180.-128., 150., 180.-38.0, 150.);
  gMC->Gsposp("BR1_",5,"ALIC", -25.3, 14.6, -76.2, irotm[41],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",5,"ALIC", -25.35, 14.8, -75.9,irotm[41],"ONLY",par,3); 

  gMC->Matrix(irotm[1], 90., 60., 90.-10.6, 150., 180.-10.6, 150.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",9,"ALIC", -47.6, 27.5, -170.,irotm[1],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",9,"ALIC", -48., 27.85, -169.5,irotm[1],"ONLY",par,3); 

  //       Bracket #2

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[42], 90., 240., 180.-128., 330., 180.-38., 330.);
  gMC->Gsposp("BR1_",6,"ALIC", 25.3,-14.6, -76.2,irotm[42],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",6,"ALIC", 25.35,-14.8, -75.9,irotm[42],"ONLY",par,3); 

  gMC->Matrix(irotm[2], 90., 240., 90.-10.6, 330., 180.-10.6, 330.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",10,"ALIC", 47.6, -27.5, -170.0,irotm[2],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",10,"ALIC", 48., -27.85, -169.5,irotm[2],"ONLY",par,3); 

  //       Bracket #3

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[43], 90., 330., 180.-128., 60., 180.-38., 60.);
  gMC->Gsposp("BR1_",7,"ALIC", 14.6, 25.3, -76.2, irotm[43],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",7,"ALIC", 14.8, 25.35, -75.9,irotm[43],"ONLY",par,3); 

  gMC->Matrix(irotm[3], 90., 330., 90.-10.6, 60., 180.-10.6, 60.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",11,"ALIC", 27.5, 47.6, -170.,irotm[3],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",11,"ALIC", 27.85, 48., -169.5,irotm[3],"ONLY",par,3); 

  //       Bracket #4

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 17.5;
  gMC->Matrix(irotm[44], 90., 150., 180.-128., 240., 180.-38., 240.);
  gMC->Gsposp("BR1_",8,"ALIC", -14.6, -25.3, -76.2, irotm[44],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB1_",8,"ALIC", -14.8, -25.35, -75.9,irotm[44],"ONLY",par,3); 

  gMC->Matrix(irotm[4], 90., 150., 90.-10.6, 240., 180.-10.6, 240.);
  par[1] = 0.1;
  par[2] = 81.4;
  gMC->Gsposp("BR5_",12,"ALIC", -27.5, -47.6, -170.,irotm[4],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB5_",12,"ALIC", -27.85, -48., -169.5,irotm[4],"ONLY",par,3); 

  // *********************************************************** 
  //       DEFINE LEFT DISK#4 OF FMD 
  // *********************************************************** 

  gMC->Gsvolu("BR4_", "BOX ", idtmed[4], par, 0);
  gMC->Gsvolu("CB4_", "BOX ", idtmed[5], par, 0);

  //       Define parameters 

  rin  = 4.5;
  rout = 10.5;
  z    = -229.5;

  //       Ring #1 

  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L4", "TUBE", idtmed[1], par, 3);
  gMC->Gspos("R1L4", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #2 

  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L4", "TUBE", idtmed[2], par, 3);
  gMC->Gspos("R2L4", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");

  //       Ring #3 

  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L4", "TUBE", idtmed[3], par, 3);
  gMC->Gspos("R3L4", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");

  //       Bracket #1

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 31.25;
  gMC->Matrix(irotm[71], 90., 0., 90.-71., 90., 180.-71., 90.);
  gMC->Gsposp("BR4_",1,"ALIC", 0., 40.25, -240., irotm[71],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB4_",1,"ALIC", 0., 40.75, -239.5, irotm[71],"ONLY",par,3); 

  //       Bracket #2

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 31.25;
  gMC->Matrix(irotm[72], 90., 180., 90.-71., 270., 180.-71., 270.);
  gMC->Gsposp("BR4_",2,"ALIC", 0., -40.25, -240., irotm[72],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB4_",2,"ALIC", 0., -40.75, -239.5,irotm[72],"ONLY",par,3); 

  //       Bracket #3

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 31.25;
  gMC->Matrix(irotm[73], 90., 270., 90.-71., 0., 180.-71., 0.);
  gMC->Gsposp("BR4_",3,"ALIC", 40.25, 0., -240., irotm[73],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB4_",3,"ALIC", 40.75, 0., -239.5,irotm[73],"ONLY",par,3); 

  //       Bracket #4

  par[0] = 1.5;
  par[1] = 0.1;
  par[2] = 31.25;
  gMC->Matrix(irotm[74], 90., 90., 90.-71., 180., 180.-71., 180.);
  gMC->Gsposp("BR4_",4,"ALIC", -40.25, 0., -240.,irotm[74],"ONLY",par,3); 
  par[1] = 0.5;
  gMC->Gsposp("CB4_",4,"ALIC", -40.75, 0., -239.5,irotm[74],"ONLY",par,3); 
}

//_____________________________________________________________________________
void AliFMDv1::DrawModule()
{
  //
  // Draw a shaded view of the FMD version 1
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("R1R3","SEEN",1);
  gMC->Gsatt("R2R3","SEEN",1);
  gMC->Gsatt("R3R3","SEEN",1);
  gMC->Gsatt("R1L3","SEEN",1);
  gMC->Gsatt("R2L3","SEEN",1);
  gMC->Gsatt("R3L3","SEEN",1);
  gMC->Gsatt("R1R2","SEEN",1);
  gMC->Gsatt("R2R2","SEEN",1);
  gMC->Gsatt("R3R2","SEEN",1);
  gMC->Gsatt("R1L2","SEEN",1);
  gMC->Gsatt("R2L2","SEEN",1);
  gMC->Gsatt("R3L2","SEEN",1);
  gMC->Gsatt("R1R1","SEEN",1);
  gMC->Gsatt("R2R1","SEEN",1);
  gMC->Gsatt("R3R1","SEEN",1);
  gMC->Gsatt("R1L1","SEEN",1);
  gMC->Gsatt("R2L1","SEEN",1);
  gMC->Gsatt("R3L1","SEEN",1);
  gMC->Gsatt("R1L4","SEEN",1);
  gMC->Gsatt("R2L4","SEEN",1);
  gMC->Gsatt("R3L4","SEEN",1);
  gMC->Gsatt("BR1_","SEEN",1);
  gMC->Gsatt("BR2_","SEEN",1);
  gMC->Gsatt("BR3_","SEEN",1);
  gMC->Gsatt("BR4_","SEEN",1);
  gMC->Gsatt("BR5_","SEEN",1);
  gMC->Gsatt("CB1_","SEEN",1);
  gMC->Gsatt("CB2_","SEEN",1);
  gMC->Gsatt("CB3_","SEEN",1);
  gMC->Gsatt("CB4_","SEEN",1);
  gMC->Gsatt("CB5_","SEEN",1);
  gMC->Gsatt("R1SP","SEEN",1);
  gMC->Gsatt("R2SP","SEEN",1);
  gMC->Gsatt("R3SP","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 6, 9, .08, .08);
  gMC->Gdhead(1111, "Forward Multiplicity Detector version 1");
  gMC->Gdman(13, 9, "MAN");
}

//_____________________________________________________________________________
void AliFMDv1::CreateMaterials()
{
  //
  // Create Materials for version 1 of FMD
  //

  //     Material for ring #1 
  Float_t ar1[8] = { 55.8,58.7,52.,47.9,16.,28.,207.2,27. };
  Float_t zr1[8] = { 26.,28.,24.,22.,8.,14.,82.,13. };
  Float_t wr1[8] = { .27,.081,.054,.045,.18,.25,.06,.06 };
  //     Material for ring #2 
  Float_t ar2[3] = { 55.8,27.,16. };
  Float_t zr2[3] = { 26.,13.,8. };
  Float_t wr2[3] = { .35,.34,.31 };
  //     Material for ring #3 
  Float_t ar3[3] = { 28.,27.,16. };
  Float_t zr3[3] = { 14.,13.,8. };
  Float_t wr3[3] = { .37,.33,.3 };
  //     Material for brackets 
  Float_t abr[2] = { 1.,12. };
  Float_t zbr[2] = { 1.,6. };
  Float_t wbr[2] = { .1,.9 };
  //     Material for cables 
  Float_t acb[4] = { 1.,12.,37.,63.54 };
  Float_t zcb[4] = { 1.,6.,17.,2.9 };
  Float_t wcb[4] = { .02,.14,.12,.72 };
  
  Float_t epsil, stmin, deemax, tmaxfd, stemax;
  
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  //     Ring #1 
  
  AliMixture(1, "FMD_R1$", ar1, zr1, 2.69, 8, wr1);
  
  //     Ring #2 
  
  AliMixture(2, "FMD_R2$", ar2, zr2, 2.63, 3, wr2);
  
  //     Ring #3 
  
  AliMixture(3, "FMD_R3$", ar3, zr3, 3.15, 3, wr3);

  //     Brackets 
  
  AliMixture(4, "FMD_BR$", abr, zbr, 1.8, 2, wbr);

  //     Cables
  
  AliMixture(5, "FMD_CB$", acb, zcb, 3.11, 4, wcb);

  // ******************************************************* 
  //     Defines tracking media parameters. 
  // ******************************************************* 
  epsil  = .001; // Tracking precision, DLS 
  stemax = -1.;  // Maximum displacement for multiple scattering 
  tmaxfd = -20.; // Maximum angle due to field deflection 
  deemax = -.3;  // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // ******************************************************** 
  AliMedium(1, "FMD_R1_L3        ", 1, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "FMD_R2_L3        ", 2, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(3, "FMD_R3_L3        ", 3, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(4, "FMD_BR_L3        ", 4, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(5, "FMD_CB_L3        ", 5, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
}


 
