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
Revision 1.21  2000/11/30 11:48:50  kowal2
TLorentzVector.h adde to satisfy the latest changes by Federico

Revision 1.20  2000/11/14 10:48:57  kowal2
Correct material used for TSA4. Thanks to J. Barbosa.

Revision 1.19  2000/11/06 17:24:10  kowal2
Corrected bug in the outer containment vessel and
the outer field cage geometry.
Thanks to J. Barbosa.

Revision 1.18  2000/11/02 16:55:24  kowal2
Corrected bug in the inner containment vessel geometry.
Thanks to J. Belikov

Revision 1.17  2000/11/02 07:24:11  kowal2
Correction in the TPC geometry.
Changes due to the new hit structure.

Revision 1.16  2000/10/02 21:28:18  fca
Removal of useless dependecies via forward declarations

Revision 1.15  2000/07/10 20:57:39  hristov
Update of TPC code and macros by M.Kowalski

Revision 1.14  2000/06/30 12:07:50  kowal2
Updated from the TPC-PreRelease branch

Revision 1.13.2.4  2000/06/26 07:39:42  kowal2
Changes to obey the coding rules

Revision 1.13.2.3  2000/06/25 08:38:41  kowal2
Splitted from AliTPCtracking

Revision 1.13.2.2  2000/06/16 12:58:13  kowal2
Changed parameter settings

Revision 1.13.2.1  2000/06/09 07:15:07  kowal2

Defaults loaded automatically (hard-wired)
Optional parameters can be set via macro called in the constructor

Revision 1.13  2000/05/15 10:00:30  kowal2
Corrected bug in the TPC geometry, thanks to Ivana Hrivnacova

Revision 1.12  2000/04/17 09:37:33  kowal2
removed obsolete AliTPCDigitsDisplay.C

Revision 1.11.8.2  2000/04/10 08:36:12  kowal2

Updated readout chambers
Some modifications to StepManager by M. Kowalski

Revision 1.11.8.1  2000/04/10 07:56:53  kowal2
Not used anymore - removed

Revision 1.11  1999/11/04 17:28:07  fca
Correct barrel part of HV Degrader

Revision 1.10  1999/10/14 16:52:08  fca
Only use PDG codes and not GEANT ones

Revision 1.9  1999/10/08 06:27:23  fca
Corrected bug in the HV degrader geometry, thanks to G.Tabary

Revision 1.8  1999/10/04 13:39:55  fca
Correct array index problem

Revision 1.7  1999/09/29 09:24:34  fca
Introduction of the Copyright and cvs Log

*/

//
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber version 3 -- detailed TPC and slow simulation    //
//                                                                     
//Begin_Html
/*
<img src="picts/AliTPCv3Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <TMath.h>
#include "AliTPCv3.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"
#include "AliTPCDigitsArray.h"
#include"AliTPCParam.h"
#include"AliTPCParamSR.h"
#include "AliPDG.h"
#include <TInterpreter.h>
#include "TLorentzVector.h"

ClassImp(AliTPCv3)
 
//_____________________________________________________________________________
AliTPCv3::AliTPCv3(const char *name, const char *title) :
  AliTPC(name, title) 
{
  //
  // Standard constructor for Time Projection Chamber version 3
  //

  SetBufferSize(128000);
  fIdSens=0;

  SetGasMixt(2,20,10,-1,0.9,0.1,0.);

  if (fTPCParam)
     fTPCParam->Write(fTPCParam->GetTitle());
}
 
//_____________________________________________________________________________
void AliTPCv3::CreateGeometry()
{
  //
  // Creation of the TPC coarse geometry (version 3)
  // The whole drift volume is sensitive
  // Origin Marek Kowalski Cracow
  //
  //Begin_Html
  /*
    <img src="picts/AliTPCgif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTPCv3Tree.gif">
  */
  //End_Html
  Float_t dm[50];

  Int_t *idtmed = fIdtmed->GetArray(); // TPC media

  Int_t idrotm[120]; // rotation matrices

  Int_t nRotMat = 0; // actual number of rotation matrices


  //
  //  Mother volume (Air) - all volumes will be positioned in it
  //
  
  dm[0]=0.;
  dm[1]=360.;
  dm[2]=12.;

  //
 
  dm[3]= -283.7;
  dm[4]= 66.2;
  dm[5]= 277.95;

  //

  dm[6]= -255.6;
  dm[7]= 66.2;
  dm[8]= 277.95;

  //

  dm[9]= -73.3;
  dm[10]= 59.0;
  dm[11]= 277.95;

  //

  dm[12]= -73.3;
  dm[13]= 56.9;
  dm[14]= 277.95;

  //

  dm[15]= -72.1;
  dm[16]= 56.9;
  dm[17]= 277.95;

  //

  dm[18]= -72.1;
  dm[19]= 60.65;
  dm[20]= 277.95;

  //

  dm[21]= 72.1;
  dm[22]= 60.65;
  dm[23]= 277.95;  

  //

  dm[24]= 72.1;
  dm[25]= 56.9;
  dm[26]= 277.95;

  //

  dm[27]= 73.3;
  dm[28]= 56.9;
  dm[29]= 277.95;

  //

  dm[30]= 73.3;
  dm[31]= 59.;
  dm[32]= 277.95;

  //

  dm[33]= 250.4;
  dm[34]= 66.0;
  dm[35]= 277.95;

  //

  dm[36]= 283.7;
  dm[37]= 66.0;
  dm[38]= 277.95;


  gMC->Gsvolu("TPC ","PCON",idtmed[0],dm,39);


  //-------------------------------------------------------------------
  //  Tpc Outer INsulator (CO2)
  //-------------------------------------------------------------------

  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]= 6.;

  //
 
  dm[3]= -253.6;
  dm[4]= 258.;
  dm[5]= 266.65;

  //

  dm[6]= -253.;
  dm[7]= 258.;
  dm[8]= 266.65;

  dm[9]= -253.;
  dm[10]= 258.;
  dm[11]= 277.97;

  dm[12]= 253.6;
  dm[13]= 258.;
  dm[14]= 277.97;

  //

  dm[15]= 253.6;
  dm[16]= 265.2;
  dm[17]= 277.95;

  //

  dm[18]= 255.6;
  dm[19]= 265.2;
  dm[20]= 277.95;


  gMC->Gsvolu("TOIN","PCON",idtmed[3],dm,21);

  //---------------------------------------------------------------
  // shreds (G10) - TPC Rings
  //---------------------------------------------------------------

  gMC->Gsvolu("TPCR","TUBE",idtmed[12],dm,0);

  dm[0]= 258.;
  dm[1]= 266.65;
  dm[2]= 0.3;

  gMC->Gsposp("TPCR",1,"TOIN",0.,0.,-253.3,0,"ONLY",dm,3); // left bottom

  //

  dm[0]= 258.;
  dm[1]= 270.9;
  dm[2]= 0.3;

  gMC->Gsposp("TPCR",2,"TOIN",0.,0.,253.3,0,"ONLY",dm,3); // right

  //

  dm[0]= 272.2;
  dm[1]= 277.95;
  dm[2]= 0.3;

  gMC->Gsposp("TPCR",3,"TOIN",0.,0.,-250.7,0,"ONLY",dm,3); // left top

  //----------------------------------------------------------------
  // Tpc Outer Contaiment Vessel  
  //  mother volume - Al, daughters - composite (sandwich)
  //----------------------------------------------------------------
  
  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]=6.;

  //

  dm[3]= -250.4;
  dm[4]= 272.2;
  dm[5]= 277.95;

  //

  dm[6]= -248.4;
  dm[7]= 272.2;
  dm[8]= 277.95;

  //

  dm[9]= -248.4;
  dm[10]= 274.81;
  dm[11]= 277.95;

  //

  dm[12]= 253.6;
  dm[13]= 274.81;
  dm[14]= 277.95;

  //

  dm[15]= 253.6;
  dm[16]= 265.2;
  dm[17]= 277.95;

  // 

  dm[18]= 255.6;
  dm[19]= 265.2;
  dm[20]= 277.95;

  gMC->Gsvolu("TOCV","PCON",idtmed[4],dm,21);

  // Daughter volumes

  // Tpc SAndwich 1 - Al
 
  dm[0]= 274.81;
  dm[1]= 277.95;
  dm[2]= 251.7;

  gMC->Gsvolu("TSA1","TUBE",idtmed[4],dm,3);

  // Tpc SAndwich 2 - Tedlar

  dm[0] += 5.e-3;
  dm[1] -= 5.e-3;
  
  gMC->Gsvolu("TSA2","TUBE",idtmed[9],dm,3);

  // Tpc SAndwich 3 - Kevlar

  dm[0] += 5e-3;
  dm[1] -= 5.e-3;

  gMC->Gsvolu("TSA3","TUBE",idtmed[5],dm,3);

  // Tpc SAndwich 4 - NOMEX honeycomb

  dm[0] += 0.06;
  dm[1] -= 0.06;

  gMC->Gsvolu("TSA4","TUBE",idtmed[6],dm,3);  
  
  // 4->3->2->1->TOCV

  gMC->Gspos("TSA4",1,"TSA3",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TSA3",1,"TSA2",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TSA2",1,"TSA1",0.,0.,0.,0,"ONLY");

  gMC->Gspos("TSA1",1,"TOCV",0.,0.,2.6,0,"ONLY");

  // TCOV-> TOIN

  gMC->Gspos("TOCV",1,"TOIN",0.,0.,0.,0,"ONLY");

  //-------------------------------------------------------
  //  Tpc Outer Field Cage
  //  mother volume - Al, daughters - composite (sandwich)
  //-------------------------------------------------------

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=6.;

  dm[3]= -253.;
  dm[4]= 258.;
  dm[5]= 277.95;

  //

  dm[6]= -251.;
  dm[7]= 258.;
  dm[8]= 277.95;

  //
 
  dm[9]= -251.;
  dm[10]= 258.;
  dm[11]= 260.05;

  //

  dm[12]= 251.;
  dm[13]= 258.;
  dm[14]= 260.05;

  //

  dm[15]= 251.;
  dm[16]= 258.;
  dm[17]= 270.9;

  //

  dm[18]= 253.;
  dm[19]= 258.;
  dm[20]= 270.9;

  gMC->Gsvolu("TOFC","PCON",idtmed[4],dm,21);  

  // Daughter volumes 

  // Tpc SAndwich 5 - Tedlar

  dm[0]= 258.;
  dm[1]= 260.05;
  dm[2]= 251.7;

  gMC->Gsvolu("TSA5","TUBE",idtmed[9],dm,3);

  // Tpc SAndwich 6 - Kevlar

  dm[0] += 5.e-3;
  dm[1] -= 5.e-3;

  gMC->Gsvolu("TSA6","TUBE",idtmed[5],dm,3);


  // Tpc SAndwich 7 - NOMEX

  dm[0] += 0.02;
  dm[1] -= 0.02;

  gMC->Gsvolu("TSA7","TUBE",idtmed[6],dm,3);    

  // 7->6->5->TOFC

  gMC->Gspos("TSA7",1,"TSA6",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TSA6",1,"TSA5",0.,0.,0.,0,"ONLY"); 

  gMC->Gspos("TSA5",1,"TOFC",0.,0.,0.,0,"ONLY");

  // TOFC->TOIN

  gMC->Gspos("TOFC",1,"TOIN",0.,0.,0.,0,"ONLY"); 

  // TOIN->TPC

  gMC->Gspos("TOIN",1,"TPC ",0.,0.,0.,0,"ONLY");

  //--------------------------------------------------------------------
  // Tpc Inner INsulator (CO2)
  //--------------------------------------------------------------------


  dm[0]=0.;
  dm[1]= 360.;
  dm[2]= 15.;

  //

  dm[3]= -255.6;
  dm[4]= 66.2;
  dm[5]= 74.8;

  //

  Float_t tanL = (66.2-59.0)/(255.6-73.3); // tangent of the left cone part

  dm[6]= -253.6;
  dm[7]= 59.0+ (253.6-73.3)*tanL;
  dm[8]= 74.8;

  //

  dm[9]= -253.6;
  dm[10]= dm[7];
  dm[11]= 79.2;

  //

  dm[12]= -73.3;
  dm[13]= 59.0;
  dm[14]= 79.2;

  //

  dm[15]= -73.3;
  dm[16]= 56.9;
  dm[17]= 79.2;

  //

  dm[18]= -72.1;
  dm[19]= 59.6;
  dm[20]= 79.2;

  //

  dm[21]= -72.1;
  dm[22]= 60.65;
  dm[23]= 79.2;

  //

  dm[24]= 72.1;
  dm[25]= 60.65;
  dm[26]= 79.2;  

  //

  dm[27]= 72.1;
  dm[28]= 59.6;
  dm[29]= 79.2;

  //

  dm[30]= 73.3;
  dm[31]= 56.9;
  dm[32]= 79.2;

  //

  dm[33]= 73.3;
  dm[34]= 59.0;
  dm[35]= 79.2;

  //

  dm[36]= 250.4;
  dm[37]= 66.0;
  dm[38]= 79.2;

  //
  
  dm[39]= 253.0;
  dm[40]= 66.0;
  dm[41]= 79.2;

  //

  dm[42]= 253.0;
  dm[43]= 75.3;
  dm[44]= 79.2;

  //

  dm[45]= 253.6;
  dm[46]= 75.3;
  dm[47]= 79.2;

  gMC->Gsvolu("TIIN","PCON",idtmed[3],dm,48);


  //--------------------------------------------------------------------
  //  Tpc Inner Containment vessel, Left part
  //  mother volume - Al, daughter - composite (sandwich)
  //--------------------------------------------------------------------

  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]= 8.;

  //

  dm[3]= -255.6;
  dm[4]= 66.2;
  dm[5]= 74.8;

  //

  Float_t cosL = 1./TMath::Sqrt(1.+tanL*tanL);
  Float_t sandThick = 2.14; // cone composite thickness


  //

  dm[6]= -253.6;
  dm[7]= 59.0 + (253.6-73.3)*tanL;
  dm[8]= 74.8;

  //

  dm[9]= -253.6;
  dm[10]= dm[7];
  dm[11]= dm[7]+sandThick/cosL;

  //

  dm[12]= -75.6;
  dm[13]= 59.0+(75.6-73.3)*tanL;
  dm[14]= dm[13]+sandThick/cosL;

  //

  dm[15]= -75.6;
  dm[16]= dm[13];
  dm[17]= 60.65;

  //

  dm[18]= -73.3;
  dm[19]= 59.0;
  dm[20]= 60.65;

  //

  dm[21]= -73.3;
  dm[22]= 56.9;
  dm[23]= 60.65;

  //

  dm[24]= -72.1;
  dm[25]= 56.9;
  dm[26]= 60.65;

  gMC->Gsvolu("TICL","PCON",idtmed[4],dm,27);

  // Daughter volumes 

  // Tpc SAndwich 9 - Al

  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]= 2.;

  //

  dm[3]= - 254.3;
  dm[4]= 59.0+(254.3-73.3)*tanL;
  dm[5]= dm[4]+sandThick/cosL;

  //

  dm[6]= -78.3;
  dm[7]= 59.0+(78.3-73.3)*tanL;
  dm[8]= dm[7]+sandThick/cosL;

  //

  gMC->Gsvolu("TSA9","PCON",idtmed[4],dm,9);

  // Tpc SAndwich 10 - Tedlar

  dm[4]+= 5.e-3/cosL;
  dm[5]-= 5.e-3/cosL;

  //

  dm[7]+= 5.e-3/cosL;
  dm[8]+= 5.e-3/cosL;

  gMC->Gsvolu("TS10","PCON",idtmed[9],dm,9); 

  // Tpc SAndwich 11 - Kevlar

  dm[4]+= 5.e-3/cosL;
  dm[5]-= 5.e-3/cosL;

  //

  dm[7]+= 5.e-3/cosL;
  dm[8]+= 5.e-3/cosL;

  gMC->Gsvolu("TS11","PCON",idtmed[5],dm,9); 

  // Tpc SAndwich 12 - NOMEX

  dm[4]+= 0.06/cosL;
  dm[5]-= 0.06/cosL;

  //

  dm[7]+= 0.06/cosL;
  dm[8]+= 0.06/cosL;

  gMC->Gsvolu("TS12","PCON",idtmed[6],dm,9);   

  // 12->11->10->9

  gMC->Gspos("TS12",1,"TS11",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS11",1,"TS10",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS10",1,"TSA9",0.,0.,0.,0,"ONLY");

  // TSA9->TICL

  gMC->Gspos("TSA9",1,"TICL",0.,0.,0.,0,"ONLY");
 
  //--------------------------------------------------------------------
  //  Tpc Inner Containment vessel, Right part
  //  mother volume - Al, daughter - composite (sandwich)
  //--------------------------------------------------------------------

  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]=8.;

  //

  dm[3]= 72.1;
  dm[4]= 56.9;
  dm[5]= 60.65;

  //

  dm[6]= 73.3;
  dm[7]= 56.9;
  dm[8]= 60.65;

  //

  dm[9]= 73.3;
  dm[10]= 59.0;
  dm[11]= 60.65;

  //

  Float_t tanR = (66.0-59.0)/(250.5-73.3); // to avoid accuracy problems
  Float_t cosR = 1./TMath::Sqrt(1.+tanR*tanR); //as above

  //

  dm[12]= 75.6;
  dm[13]= 59.0+(75.6-73.3)*tanR;
  dm[14]= 60.65;

  //

  dm[15]= 75.6;
  dm[16]= dm[13];
  dm[17]= dm[16]+sandThick/cosR;

  //

  dm[18]= 248.4;
  dm[19]= 59.0+(248.4-73.3)*tanR;
  dm[20]= dm[19]+sandThick/cosR;

  //

  dm[21]= 248.4;
  dm[22]= dm[19];
  dm[23]= 70.2;

  //

  dm[24]= 250.4;
  dm[25]= 66.0;
  dm[26]= 70.2;

  gMC->Gsvolu("TICR","PCON",idtmed[4],dm,27);



  // Daughter volumes 

  // Tpc SAndwich 13 - Al

  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]= 2.;

  //

  dm[3]= 78.3;
  dm[4]= 59.0+(78.3-73.3)*tanR;
  dm[5]= dm[4]+sandThick/cosR;

  //

  dm[6]= 249.1;
  dm[7]= 59.0+(249.1-73.3)*tanR;
  dm[8]= dm[7]+sandThick/cosR;

  //

  gMC->Gsvolu("TS13","PCON",idtmed[4],dm,9);

  // Tpc SAndwich 14 - Tedlar

  dm[4]+= 5.e-3/cosR;
  dm[5]-= 5.e-3/cosR;

  //

  dm[7]+= 5.e-3/cosR;
  dm[8]+= 5.e-3/cosR;

  gMC->Gsvolu("TS14","PCON",idtmed[9],dm,9); 

  // Tpc SAndwich 15 - Kevlar

  dm[4]+= 5.e-3/cosR;
  dm[5]-= 5.e-3/cosR;

  //

  dm[7]+= 5.e-3/cosR;
  dm[8]+= 5.e-3/cosR;

  gMC->Gsvolu("TS15","PCON",idtmed[5],dm,9); 

  // Tpc SAndwich 16 - NOMEX

  dm[4]+= 0.06/cosR;
  dm[5]-= 0.06/cosR;

  //

  dm[7]+= 0.06/cosR;
  dm[8]+= 0.06/cosR;

  gMC->Gsvolu("TS16","PCON",idtmed[6],dm,9);   

  // 16->15->14->13

  gMC->Gspos("TS16",1,"TS15",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS15",1,"TS14",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS14",1,"TS13",0.,0.,0.,0,"ONLY");

  // TS12->TICR

  gMC->Gspos("TS13",1,"TICR",0.,0.,0.,0,"ONLY");

  //------------------------------------------------------
  // Tpc Inner Field Cage
  // mother volume - Al, daughters - composite (sandwich)
  //------------------------------------------------------

  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]=6.;

  //

  dm[3]= -253.0;
  dm[4]= 70.7;
  dm[5]= 79.2;

  //

  dm[6]= -251.0;
  dm[7]= 70.7;
  dm[8]= 79.2;

  //

  dm[9]= -251.0;
  dm[10]= 77.15;
  dm[11]= 79.2;

  //
  
  dm[12]= 251.0;
  dm[13]= 77.15;
  dm[14]= 79.2;

  //

  dm[15]= 251.0;
  dm[16]= 66.0;
  dm[17]= 79.2;

  //

  dm[18]= 253.0;
  dm[19]= 66.0;
  dm[20]= 79.2;

  gMC->Gsvolu("TIFC","PCON",idtmed[4],dm,21);

  // Daughter volumes

  // Tpc Sandwich 17 - Tedlar

  dm[0]= 77.15;
  dm[1]= 79.2;
  dm[2]= 251.7;

  gMC->Gsvolu("TS17","TUBE",idtmed[9],dm,3);

  // Tpc Sandwich 18 - Kevlar

  dm[0]+= 5.e-3;
  dm[1]-= 5.e-3;

  gMC->Gsvolu("TS18","TUBE",idtmed[5],dm,3);


  // Tpc Sandwich 19 - NOMEX

  dm[0]+= 0.02;
  dm[1]-= 0.02;

  gMC->Gsvolu("TS19","TUBE",idtmed[6],dm,3);

  // 19->18->17

  gMC->Gspos("TS19",1,"TS18",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS18",1,"TS17",0.,0.,0.,0,"ONLY");

  // TS17->TIFC

  gMC->Gspos("TS17",1,"TIFC",0.,0.,0.,0,"ONLY");

  // TPC Rings

  dm[0]= 70.7;
  dm[1]= 79.2;
  dm[2]= 0.3;

  gMC->Gsposp("TPCR",4,"TIIN",0.,0.,-253.3,0,"ONLY",dm,3);

  dm[0]= 66.0;
  dm[1]= 70.2;  

  gMC->Gsposp("TPCR",5,"TIIN",0.,0.,250.7,0,"ONLY",dm,3);

  dm[0]= 75.3;
  dm[1]= 79.2;  

  gMC->Gsposp("TPCR",6,"TIIN",0.,0.,253.3,0,"ONLY",dm,3);  

  // TICL->TIIN

  gMC->Gspos("TICL",1,"TIIN",0.,0.,0.,0,"ONLY");

  // TICR->TIIN

  gMC->Gspos("TICR",1,"TIIN",0.,0.,0.,0,"ONLY");

  // TIFC->TIIN

  gMC->Gspos("TIFC",1,"TIIN",0.,0.,0.,0,"ONLY");

  // Tpc Sandwich 21 - Al (central barrel)

  dm[0]= 60.65;
  dm[1]= 61.21;
  dm[2]= 75.2;

  gMC->Gsvolu("TS21","TUBE",idtmed[4],dm,3);

  // Tpc Sandwich 22 - Tedlar (central barrel) 

  dm[0]+= 5.e-3;
  dm[1]-= 5.e-3;

  gMC->Gsvolu("TS22","TUBE",idtmed[9],dm,3); 

  // Tpc Sandwich 23 - Kevlar (central barrel) 

  dm[0]+= 5.e-3;
  dm[1]-= 5.e-3;

  gMC->Gsvolu("TS23","TUBE",idtmed[5],dm,3); 

  // Tpc Sandwich 24 - NOMEX (central barrel) 

  dm[0]+= 0.02;
  dm[1]-= 0.02;

  gMC->Gsvolu("TS24","TUBE",idtmed[6],dm,3); 

  // 24->23->22->21

  gMC->Gspos("TS24",1,"TS23",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS23",1,"TS22",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS22",1,"TS21",0.,0.,0.,0,"ONLY");

  gMC->Gspos("TS21",1,"TIIN",0.,0.,0.,0,"ONLY");

  // put everything into the TPC 

  gMC->Gspos("TIIN",1,"TPC ",0.,0.,0.,0,"ONLY");


  //---------------------------------------------------------
  //  Tpc Dift Gas volume Nonsensitive (Ne-CO2 90/10)
  //  and its daughters (HV membrane, rods, readout chambers)
  //---------------------------------------------------------

  dm[0]= 79.2;
  dm[1]= 258.0;
  dm[2]= 253.6;

  gMC->Gsvolu("TDGS","TUBE",idtmed[2],dm,3);  

  // sector opening angles

  Float_t innerOpenAngle = fTPCParam->GetInnerAngle();

  // sector angle shift

  Float_t innerAngleShift = fTPCParam->GetInnerAngleShift();

  // number of sectors

  Int_t nInnerSector = fTPCParam->GetNInnerSector()/2;
  Int_t nOuterSector = fTPCParam->GetNOuterSector()/2;

  // All above parameters are identical for inner and outer
  // sectors. The distinction is kept for the historical reasons
  // and eventually will disappear.

  Float_t tanAlpha = TMath::Tan(0.5*innerOpenAngle);
  Float_t cosAlpha = TMath::Sqrt(1.+tanAlpha*tanAlpha);
  Float_t space;

  //-------------------------------------------------------------------------
  //   Tpc Inner Readout Chambers 
  //-------------------------------------------------------------------------

  dm[0]= 14.483;
  dm[1]= 23.3345; 
  dm[2]= 1.6; // thickness
  dm[3]= 25.1;

  gMC->Gsvolu("TIRC","TRD1",idtmed[4],dm,4);

  // this volume will be positioned in the empty space
  // of the end-cap to avoid overlaps

  dm[0]= 13.7305;
  dm[1]= 21.1895;
  dm[2]= 2.25;
  dm[3]= 21.15;

  gMC->Gsvolu("TIC1","TRD1",idtmed[4],dm,4);


  //------------------------------------------------
  // Tpc Inner readout chamber Pad Plane
  //------------------------------------------------

  dm[0]= 14.483;
  dm[1]= 23.3345;
  dm[2]= 0.5;
  dm[3]= 25.1;

  gMC->Gsvolu("TIPP","TRD1",idtmed[12],dm,4);

  // 

  dm[0] -= 1.218511934;
  dm[1] -= 1.218511934;
  dm[2] = 0.35;

  gMC->Gsvolu("TIC3","TRD1",idtmed[2],dm,4);

  gMC->Gspos("TIC3",1,"TIPP",0.,0.15,0.,0,"ONLY");

  gMC->Gspos("TIPP",1,"TIRC",0.,1.1,0.,0,"ONLY");


  //----------------------------------------------
  // Tpc Readout Chambers Empty spaces - for both
  // inner and outer sectors
  //----------------------------------------------

  gMC->Gsvolu("TRCE","TRD1",idtmed[0],dm,0);

  // Inner sector - 4 spaces


  dm[3] = 4.7625;
  dm[0] = 12.472;

  Float_t rr = 90.52;
  Float_t zz;

  zz= -12.7875;
  
  space = rr*tanAlpha-dm[0];

  for(Int_t nsLow=0;nsLow<4;nsLow++){

    rr += 9.525;
    dm[1]= rr*tanAlpha - space;  

    dm[2]=0.8;

    gMC->Gsposp("TRCE",nsLow+1,"TIRC",0.,-0.8,zz,0,"ONLY",dm,4);

    //

    dm[2]= 1.2;

    gMC->Gsposp("TRCE",nsLow+5,"TIC1",0.,1.05,zz-2.1,0,"ONLY",dm,4);

    rr += 0.4;
    dm[0] = rr*tanAlpha - space;
    zz += (0.4+9.525); 

  }

  dm[0]= 12.472;
  // dm[1] - this is the dm[1] from the previous TRCE
  dm[2]= 1.05;
  dm[3]= 19.65;

  gMC->Gsposp("TRCE",9,"TIC1",0.,-1.,0.,0,"ONLY",dm,4);   

  //
  // TPc Space for Connectors
  //

  dm[0]= .3;
  dm[1]= .3;
  dm[2]= 4.5;

  gMC->Gsvolu("TPSC","BOX ",idtmed[0],dm,3);

  // TPC Connectors

  dm[0]= .25;
  dm[1]= .15;
  dm[2]= 3.75;

  gMC->Gsvolu("TPCC","BOX ",idtmed[13],dm,3); 

  gMC->Gspos("TPCC",1,"TPSC",0.,0.15,0.,0,"ONLY");

  zz = -12.7875;


  Float_t alpha;
  Float_t astep;

  Float_t phi1,phi2,phi3,theta1,theta2,theta3; // rotation angles

  // inner part of the inner sector - 2 x 20 holes
  
  astep = 20.00096874/19.;

  alpha = 10.00048437-astep;

  Float_t x1,x2;

    x1 = 13.31175725;
    x1 -= 0.996357832; 

    x2 = 15.06180253;
    x2 -= 1.163028812;

  Int_t ncon;

  for(ncon=0;ncon<20;ncon++){

    phi1 = 0.;
    theta1 = 90.+alpha;
    phi2=90.;
    theta2 = 90.;
    phi3 = (alpha>0) ? 0. : 180.;
    theta3 = TMath::Abs(alpha);

    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);

 

    gMC->Gspos("TPSC",ncon+1,"TIRC",x1,0.3,-12.7875,idrotm[nRotMat],"ONLY");
    gMC->Gspos("TPSC",ncon+21,"TIRC",x2,0.3,-2.8625,idrotm[nRotMat],"ONLY");


    x1 -= 1.296357833;
    x2 -= 1.463028812;

    alpha -= astep;   
    nRotMat++; 

  }

  // outer part of the inner sector - 2 x 25 holes

   astep = 20.00096874/24.; 
   alpha = 10.00048437-astep;

   x1 = 16.81184781;
   x1 -= 1.016295986;

   x2 = 18.5618931;
   x2 -= 1.150914854;

  for(ncon=0;ncon<25;ncon++){

    phi1 = 0.;
    theta1 = 90.+alpha;
    phi2=90.;
    theta2 = 90.;
    phi3 = (alpha>0) ? 0. : 180.;
    theta3 = TMath::Abs(alpha);

    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);

 

    gMC->Gspos("TPSC",ncon+41,"TIRC",x1,0.3,7.0625,idrotm[nRotMat],"ONLY");
    gMC->Gspos("TPSC",ncon+66,"TIRC",x2,0.3,16.9875,idrotm[nRotMat],"ONLY");


    x1 -= 1.316295986;
    x2 -= 1.450914854;

    alpha -= astep;   
    nRotMat++; 

  }  

  //--------------------------------------------------------------------------
  //  TPC Outer Readout Chambers
  //  this is NOT a final design
  //--------------------------------------------------------------------------

  dm[0]= 23.3875;
  dm[1]= 43.524;
  dm[2]= 1.5; //thickness
  dm[3]= 57.1;

  gMC->Gsvolu("TORC","TRD1",idtmed[4],dm,4);

  //------------------------------------------------
  // Tpc Outerr readout chamber Pad Plane
  //------------------------------------------------

  dm[2]= 0.5;

  gMC->Gsvolu("TOPP","TRD1",idtmed[12],dm,4);

  dm[0] -= 1.218511934;
  dm[1] -= 1.218511934;
  dm[2] = 0.35;

  gMC->Gsvolu("TOC3","TRD1",idtmed[2],dm,4);

  gMC->Gspos("TOC3",1,"TOPP",0.,0.15,0.,0,"ONLY");

  gMC->Gspos("TOPP",1,"TORC",0.,1.0,0.,0,"ONLY");

  // empty space


  dm[0]= 21.035;
  dm[1]= 38.7205;
  dm[2]= 0.7; 
  dm[3]= 50.15;

  gMC->Gsposp("TRCE",10,"TORC",0.,-0.8,-2.15,0,"ONLY",dm,4);  

  dm[0]= 22.2935;
  dm[1]= 40.5085;
  dm[2]= 2.25;
  dm[3]= 51.65;

  gMC->Gsvolu("TOC1","TRD1",idtmed[4],dm,4);


  dm[0]= 21.35;
  dm[1]= 38.7205;
  dm[2]= 2.25;
  dm[3]= 50.15;

  gMC->Gsposp("TRCE",11,"TOC1",0.,0.,0.,0,"ONLY",dm,4); 


  //-----------------------------------------------
  // Tpc Services Support Wheel
  //-----------------------------------------------

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=18.;
  dm[3]=2.;

  dm[4]= -5.;
  dm[5]= 77.017;
  dm[6]= 255.267;

  dm[7]= 5.;
  dm[8]= dm[5];
  dm[9]= dm[6];

  gMC->Gsvolu("TSSW","PGON",idtmed[4],dm,10);

  // Tpc Services Wheel Cover

  dm[4]= -0.5;
  dm[7]= 0.5;

  gMC->Gsvolu("TSWC","PGON",idtmed[4],dm,10);

  // Tpc Service wheel Cover Empty space
   
  dm[0]= 10.99;
  dm[1]= 39.599;
  dm[2]= .5;
  dm[3]= 81.125;

  gMC->Gsvolu("TSCE","TRD1",idtmed[0],dm,4);

  // Tpc services Wheel Empty Spaces

  dm[0]= 13.18017507;
  dm[1]= 44.61045938;
  dm[2]= 4.;
  dm[3]= 89.125;

  gMC->Gsvolu("TWES","TRD1",idtmed[0],dm,4);

  // Tpc Services Wheel Bars

  gMC->Gsvolu("TSWB","TRD1",idtmed[4],dm,0);

  // bars-> TWES

  dm[2]= 4.;
  dm[3]= .4;

  dm[0]= 13.8149522;
  dm[1]= 13.95601379;
  
  gMC->Gsposp("TSWB",1,"TWES",0.,0.,-85.125,0,"ONLY",dm,4);

  dm[0]= 43.83462067; 
  dm[1]= 43.97568225;

  gMC->Gsposp("TSWB",2,"TWES",0.,0.,85.125,0,"ONLY",dm,4);

  // TPc ELectronics - right now 30% X0 Si

  dm[0]= 14.03813696;
  dm[1]= 43.3524075;
  dm[2]= 1.404;
  dm[3]= 83.125;

  gMC->Gsvolu("TPEL","TRD1",idtmed[11],dm,4);
  gMC->Gspos("TPEL",1,"TWES",0.,0.,0.,0,"ONLY");


  //--------------------------------------------------------------------------
  //  End caps
  //--------------------------------------------------------------------------

  // TPc Main Wheel - Al

  dm[0]= 75.3;
  dm[1]= 264.8;
  dm[2]= 3.0;

  gMC->Gsvolu("TPMW","TUBE",idtmed[4],dm,3);

  // TPc Extra Wheel (to avoid overlapping) - Al

  dm[0]= 264.8;
  dm[1]= 277.0;
  dm[2]= 1.95;

  gMC->Gsvolu("TPEW","TUBE",idtmed[4],dm,3);

  //--------------------------------------------------------------------------
  //  Tpc Empty Space for the Readout chambers
  //--------------------------------------------------------------------------  

  Float_t rLow= 86.2;
  Float_t rUp= 243.5;
  Float_t dR = 0.5*(rUp-rLow);

  space= 1.4/cosAlpha; // wheel ribs are 2.8 cm wide

  dm[0]= rLow*tanAlpha-space;
  dm[1]= rUp*tanAlpha-space;
  dm[2]= 3.0;
  dm[3]= dR;

  gMC->Gsvolu("TESR","TRD1",idtmed[0],dm,4);

  // TIC1->TESR


  gMC->Gspos("TIC1",1,"TESR",0.,0.75,-dR+23.97,0,"ONLY");


  // TOC1->TESR

  gMC->Gspos("TOC1",1,"TESR",0.,0.75,dR-55.02,0,"ONLY");

  // Tpc Empty Space Bars - Al (daughters of TESR)

  Float_t zBar;

  gMC->Gsvolu("TESB","TRD1",idtmed[4],dm,0);

  // lower bar

  dm[0]= rLow*tanAlpha-space;
  dm[1]= 88.7*tanAlpha-space;
  dm[2]= 0.95;
  dm[3]= 1.275;

  zBar = -dR+dm[3];

  gMC->Gsposp("TESB",1,"TESR",0.,2.05,zBar,0,"ONLY",dm,4);

  // middle bar

  dm[0]= 131.65*tanAlpha-space;
  dm[1]= 136.5*tanAlpha-space;
  dm[3]= 2.425;

  zBar = -dR +131.65+dm[3]-rLow;

  gMC->Gsposp("TESB",2,"TESR",0.,2.05,zBar,0,"ONLY",dm,4);

  // upper bar

  dm[0]= 240.4*tanAlpha-space;
  dm[1]= rUp*tanAlpha-space;
  dm[3]= 1.55;

  zBar = dR-dm[3];

  gMC->Gsposp("TESB",3,"TESR",0.,2.05,zBar,0,"ONLY",dm,4);

  //  positioning of the empty spaces into the main wheel

  Float_t rCenter,xc,yc;
  Float_t rInner,rOuter; // center of the inner and outer chamber

  rCenter = rLow+dR;

  rInner = 108.07;
  rOuter = 190.68;

  for(Int_t ns=0; ns<nInnerSector;ns++){

    phi1 = ns * innerOpenAngle + innerAngleShift;
    phi1 *= kRaddeg; // in degrees

    phi1 = (Float_t)TMath::Nint(phi1) + 270.;

    if (phi1 > 360.) phi1 -= 360.;

    theta1 = 90.;
    phi2   = 90.;
    theta2 = 180.;
    phi3   = ns * innerOpenAngle + innerAngleShift;
    phi3 *= kRaddeg; // in degrees

    phi3 = (Float_t)TMath::Nint(phi3);
      
    if(phi3 > 360.) phi3 -= 360.;

    theta3 = 90.;

    // "holes"->End plate

    xc = rCenter*TMath::Cos(phi3*kDegrad);
    yc = rCenter*TMath::Sin(phi3*kDegrad);

    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);

    gMC->Gspos("TESR",ns+1,"TPMW",xc,yc,0.,idrotm[nRotMat],"ONLY");

    // TSCE->TSWC (services wheel volumes)

    xc = 166.142*TMath::Cos(phi3*kDegrad);
    yc = 166.142*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TSCE",ns+1,"TSWC",xc,yc,0.,idrotm[nRotMat],"ONLY");
    gMC->Gspos("TWES",ns+1,"TSSW",xc,yc,0.,idrotm[nRotMat],"ONLY");


    // readout chambers->TDGS (drift gas)

    xc = rInner*TMath::Cos(phi3*kDegrad);
    yc = rInner*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TIRC",ns+1,"TDGS",xc,yc,252.,idrotm[nRotMat],"ONLY");
     
    xc = rOuter*TMath::Cos(phi3*kDegrad);
    yc = rOuter*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TORC",ns+1,"TDGS",xc,yc,252.1,idrotm[nRotMat],"ONLY");

    nRotMat++;

    theta2 = 0.; // reflection

    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);

    xc = rInner*TMath::Cos(phi3*kDegrad);
    yc = rInner*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TIRC",ns+nInnerSector+1,"TDGS",xc,yc,-252.,idrotm[nRotMat],"ONLY");

    xc = rOuter*TMath::Cos(phi3*kDegrad);
    yc = rOuter*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TORC",ns+nOuterSector+1,"TDGS",xc,yc,-252.1,idrotm[nRotMat],"ONLY");

    nRotMat++;


  } 


  // reflection matrix

  theta1 = 90.;
  phi1   = 0.;
  theta2 = 90.;
  phi2   = 270.;
  theta3 = 180.;
  phi3   = 0.;

  AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);


  // TPMW->TPC

  gMC->Gspos("TPMW",1,"TPC ",0.,0.,256.6,0,"ONLY");
  gMC->Gspos("TPMW",2,"TPC ",0.,0.,-256.6,idrotm[nRotMat],"ONLY");
  gMC->Gspos("TPEW",1,"TPC ",0.,0.,257.65,0,"ONLY");
  gMC->Gspos("TPEW",2,"TPC ",0.,0.,-257.65,0,"ONLY");



  //-------------------------------------------------------
  // Tpc High Voltage Membrane - NOMEX honeycomb
  //-------------------------------------------------------

  dm[0]=0.,
  dm[1]=360.;
  dm[2]=18.;
  dm[3]=2.;

  //

  dm[4]= -0.3;
  dm[5]= 81.156;
  dm[6]= 253.386;

  //

  dm[7]= 0.3;
  dm[8]= dm[5];
  dm[9]= dm[6];

  gMC->Gsvolu("THVM","PGON",idtmed[6],dm,10);

  gMC->Gspos("THVM",1,"TDGS",0.,0.,0.,0,"ONLY");

  //----------------------------------------------------------
  // TPc Support Rods - MAKROLON
  //----------------------------------------------------------

  dm[0]= 0.9;
  dm[1]= 1.2;
  dm[2]= 126.65;

  gMC->Gsvolu("TPSR","TUBE",idtmed[7],dm,3);

  for(Int_t nrod=1;nrod<18;nrod++){
    Float_t angle=innerOpenAngle*(Float_t)nrod;

    xc=82.4*TMath::Cos(angle);
    yc=82.4*TMath::Sin(angle);

    gMC->Gspos("TPSR",nrod,"TDGS",xc,yc,126.95,0,"ONLY");
    gMC->Gspos("TPSR",nrod+17,"TDGS",xc,yc,-126.95,0,"ONLY");

    xc=254.2*TMath::Cos(angle);
    yc=254.2*TMath::Sin(angle);

    gMC->Gspos("TPSR",nrod+34,"TDGS",xc,yc,126.95,0,"ONLY");
    gMC->Gspos("TPSR",nrod+51,"TDGS",xc,yc,-126.95,0,"ONLY");    

  }

  //----------------------------------------------------------
  // Tpc High Voltage rod - MAKROLON + Copper cable
  //----------------------------------------------------------

  // rod with cable (Left)

  dm[0]=0.;
  dm[1]=2.25;
  dm[2]=126.65;

  gMC->Gsvolu("THVL","TUBE",idtmed[7],dm,3);

  // HV cable
 
  dm[0]=0.;
  dm[1]=0.3;
  dm[2]=126.65;

  gMC->Gsvolu("THVC","TUBE",idtmed[10],dm,3);

  // empty space

  dm[0]=0.3;
  dm[1]=1.;
  dm[2]=126.65;

  gMC->Gsvolu("THVE","TUBE",idtmed[1],dm,3);

  gMC->Gspos("THVC",1,"THVL",0.,0.,0.,0,"ONLY");
  gMC->Gspos("THVE",1,"THVL",0.,0.,0.,0,"ONLY");

  // rod without cable

  dm[0]=1.8;
  dm[1]=2.25;
  dm[2]=126.65;

  gMC->Gsvolu("THVR","TUBE",idtmed[7],dm,3);

  
  
  gMC->Gspos("THVL",1,"TDGS",82.4,0.,-126.95,0,"ONLY");
  gMC->Gspos("THVL",2,"TDGS",254.2,0.,-126.95,0,"ONLY");

  gMC->Gspos("THVR",1,"TDGS",82.4,0.,126.95,0,"ONLY");
  gMC->Gspos("THVR",2,"TDGS",254.2,0.,126.95,0,"ONLY");  
  


  gMC->Gspos("TDGS",1,"TPC ",0.,0.,0.,0,"ONLY"); 

  // services wheel cover -> wheel


  gMC->Gspos("TSWC",1,"TSSW",0.,0.,4.5,0,"ONLY");
  gMC->Gspos("TSWC",2,"TSSW",0.,0.,-4.5,0,"ONLY");


  // put the wheel into the TPC

  gMC->Gspos("TSSW",1,"TPC ",0.,0.,278.7,0,"ONLY");
  gMC->Gspos("TSSW",2,"TPC ",0.,0.,-278.7,0,"ONLY");

  gMC->Gsord("TPMW",6);
  gMC->Gsord("TSSW",6);
  gMC->Gsord("TSWC",6);

  // put the TPC into ALIC (main mother volume)

  gMC->Gspos("TPC ",1,"ALIC",0.,0.,0.,0,"ONLY");

 
} // end of function
 
//_____________________________________________________________________________
void AliTPCv3::DrawDetector()
{
  //
  // Draw a shaded view of the Time Projection Chamber version 1
  //


  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  //

  gMC->Gsatt("TPC ","SEEN",0);
  gMC->Gsatt("TOIN","SEEN",1);
  gMC->Gsatt("TOIN","COLO",7);
  gMC->Gsatt("TPCR","SEEN",0);
  gMC->Gsatt("TOCV","SEEN",1);
  gMC->Gsatt("TOCV","COLO",4);
  gMC->Gsatt("TSA1","SEEN",0);
  gMC->Gsatt("TSA2","SEEN",0);
  gMC->Gsatt("TSA3","SEEN",0);
  gMC->Gsatt("TSA4","SEEN",0);
  gMC->Gsatt("TOFC","SEEN",1);
  gMC->Gsatt("TOFC","COLO",4);
  gMC->Gsatt("TSA5","SEEN",0);
  gMC->Gsatt("TSA6","SEEN",0);
  gMC->Gsatt("TSA7","SEEN",0);
  gMC->Gsatt("TIIN","COLO",7);
  gMC->Gsatt("TIIN","SEEN",1);
  gMC->Gsatt("TICL","SEEN",0);
  gMC->Gsatt("TSA9","SEEN",0);
  gMC->Gsatt("TS10","SEEN",0);
  gMC->Gsatt("TS11","SEEN",0);
  gMC->Gsatt("TS12","SEEN",0);
  gMC->Gsatt("TICR","SEEN",0); 
  gMC->Gsatt("TS13","SEEN",0);
  gMC->Gsatt("TS14","SEEN",0);
  gMC->Gsatt("TS15","SEEN",0);
  gMC->Gsatt("TS16","SEEN",0);
  gMC->Gsatt("TIFC","SEEN",1);
  gMC->Gsatt("TIFC","COLO",4); 
  gMC->Gsatt("TS17","SEEN",0);
  gMC->Gsatt("TS18","SEEN",0);
  gMC->Gsatt("TS19","SEEN",0);
  gMC->Gsatt("TS21","SEEN",0);
  gMC->Gsatt("TS22","SEEN",0);
  gMC->Gsatt("TS23","SEEN",0);
  gMC->Gsatt("TS24","SEEN",0);
  gMC->Gsatt("TDGS","SEEN",0);
  gMC->Gsatt("TIRC","SEEN",0);
  gMC->Gsatt("TIC1","SEEN",1);
  gMC->Gsatt("TIPP","SEEN",0);
  gMC->Gsatt("TIC3","SEEN",0);
  gMC->Gsatt("TRCE","SEEN",0);
  gMC->Gsatt("TPSC","SEEN",0);
  gMC->Gsatt("TPCC","SEEN",0);
  gMC->Gsatt("TORC","SEEN",0);
  gMC->Gsatt("TOPP","SEEN",0);
  gMC->Gsatt("TOC3","SEEN",0);
  gMC->Gsatt("TOC1","SEEN",1);
  gMC->Gsatt("TSSW","SEEN",1);
  gMC->Gsatt("TSWC","SEEN",1);
  gMC->Gsatt("TSCE","SEEN",1); 
  gMC->Gsatt("TSSW","COLO",3);
  gMC->Gsatt("TSWC","COLO",3);
  gMC->Gsatt("TSCE","COLO",6);
  gMC->Gsatt("TWES","SEEN",0);
  gMC->Gsatt("TSWB","SEEN",0);
  gMC->Gsatt("TPEL","SEEN",0);
  gMC->Gsatt("TPMW","SEEN",1);
  gMC->Gsatt("TPEW","SEEN",1);
  gMC->Gsatt("TESR","SEEN",1);
  gMC->Gsatt("TPMW","COLO",12);
  gMC->Gsatt("TPEW","COLO",12);
  gMC->Gsatt("TWES","COLO",5);
  gMC->Gsatt("TIC1","COLO",5);
  gMC->Gsatt("TOC1","COLO",5);  
  gMC->Gsatt("TESB","SEEN",0);
  gMC->Gsatt("TPLS","SEEN",0);
  gMC->Gsatt("TPUS","SEEN",0);
  gMC->Gsatt("TPSS","SEEN",0);
  gMC->Gsatt("THVM","SEEN",1);
  gMC->Gsatt("THVM","COLO",11);
  gMC->Gsatt("TPSR","SEEN",0);
  gMC->Gsatt("THVL","SEEN",0);
  gMC->Gsatt("THVC","SEEN",0);
  gMC->Gsatt("THVE","SEEN",0);
  gMC->Gsatt("THVR","SEEN",0);

  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("TPMW",-300,300,-300,300,254.,270.);
  gMC->SetClipBox("TESR",-300,300,-300,300,254.,270.);
  gMC->SetClipBox("TSSW",-300,300,-300,300,283.,284.);
  gMC->SetClipBox("TSWC",-300,300,-300,300,283.,284.);
  gMC->SetClipBox("*", 0, 300, -300, 300, -290, 290);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .025, .025);
  gMC->Gdhead(1111, "Time Projection Chamber");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTPCv3::CreateMaterials()
{
  //
  // Define materials for version 2 of the Time Projection Chamber
  //


  //
  AliTPC::CreateMaterials();
}

//_____________________________________________________________________________
void AliTPCv3::Init()
{
  //
  // Initialises version 3 of the TPC after that it has been built
  //
  Int_t *idtmed = fIdtmed->GetArray();

  AliTPC::Init();

  fIdSens=gMC->VolId("TDGS"); // drift gas as a sensitive volume

  gMC->SetMaxNStep(30000); // max. number of steps increased

  gMC->Gstpar(idtmed[2],"LOSS",5);

  printf("*** TPC version 3 initialized ***\n");
  printf("Maximum number of steps = %d\n",gMC->GetMaxNStep());

  //
  
}

//_____________________________________________________________________________
void AliTPCv3::StepManager()
{
  //
  // Called for every step in the Time Projection Chamber
  //

  //
  // parameters used for the energy loss calculations
  //
  const Float_t kprim = 14.35; // number of primary collisions per 1 cm
  const Float_t kpoti = 20.77e-9; // first ionization potential for Ne/CO2
  const Float_t kwIon = 35.97e-9; // energy for the ion-electron pair creation 
 
 
  const Float_t kbig = 1.e10;

  Int_t id,copy;
  TLorentzVector pos;
  Float_t hits[4];
  Int_t vol[2];  
  TClonesArray &lhits = *fHits;
  
  vol[1]=0;
  vol[0]=0;

  //

  gMC->SetMaxStep(kbig);
  
  if(!gMC->IsTrackAlive()) return; // particle has disappeared
  
  Float_t charge = gMC->TrackCharge();
  
  if(TMath::Abs(charge)<=0.) return; // take only charged particles
  
  
  id=gMC->CurrentVolID(copy);
  
  // Check the sensitive volume
  
  if (id != fIdSens) return;
  
  //
  //  charged particle is in the sensitive volume
  //
  
  if(gMC->TrackStep() > 0) {

    
    Int_t nel = (Int_t)(((gMC->Edep())-kpoti)/kwIon) + 1;
    nel=TMath::Min(nel,300); // 300 electrons corresponds to 10 keV
    
    gMC->TrackPosition(pos);
    hits[0]=pos[0];
    hits[1]=pos[1];
    hits[2]=pos[2];

    //
    // check the selected side of the TPC
    //
 
    if(fSide && fSide*hits[2]<=0.) return;

    hits[3]=(Float_t)nel;
    
    // Add this hit
   
    new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
    
  } 
  
  // Stemax calculation for the next step
  
  Float_t pp;
  TLorentzVector mom;
  gMC->TrackMomentum(mom);
  Float_t ptot=mom.Rho();
  Float_t betaGamma = ptot/gMC->TrackMass();
  
  Int_t pid=gMC->TrackPid();
  if((pid==kElectron || pid==kPositron) && ptot > 0.002)
    { 
      pp = kprim*1.58; // electrons above 20 MeV/c are on the plateau!
    }
  else
    {
      pp=kprim*BetheBloch(betaGamma);    
      if(TMath::Abs(charge) > 1.) pp *= (charge*charge);
    }
  
  Float_t random[1];
  gMC->Rndm(random,1); // good, old GRNDM from Geant3
  
  Double_t rnd = (Double_t)random[0];
  
  gMC->SetMaxStep(-TMath::Log(rnd)/pp);
  
}

//_____________________________________________________________________________
Float_t AliTPCv3::BetheBloch(Float_t bg)
{
  //
  // Bethe-Bloch energy loss formula
  //
  const Double_t kp1=0.76176e-1;
  const Double_t kp2=10.632;
  const Double_t kp3=0.13279e-4;
  const Double_t kp4=1.8631;
  const Double_t kp5=1.9479;

  Double_t dbg = (Double_t) bg;

  Double_t beta = dbg/TMath::Sqrt(1.+dbg*dbg);

  Double_t aa = TMath::Power(beta,kp4);
  Double_t bb = TMath::Power(1./dbg,kp5);

  bb=TMath::Log(kp3+bb);
  
  return ((Float_t)((kp2-aa-bb)*kp1/aa));
}
