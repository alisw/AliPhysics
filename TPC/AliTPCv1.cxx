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

//
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber version 1 -- detailed TPC and fast simulation    //
//                                                                           // 
//Begin_Html
/*
<img src="picts/AliTPCv1Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include <TBRIK.h>
#include <TGeometry.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliRun.h"
#include "AliTPCDigitsArray.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCv1.h"
#include "TLorentzVector.h"
#include "AliMC.h"

ClassImp(AliTPCv1)

//_____________________________________________________________________________

  AliTPCv1::AliTPCv1(){

  fHitType = 1;

}
 
//_____________________________________________________________________________
AliTPCv1::AliTPCv1(const char *name, const char *title) 
  :AliTPC(name, title) {
  //
  // Standard constructor for Time Projection Chamber
  //
  fIdSens=fIdLSec=fIdUSec=0;

  SetGasMixt(2,20,10,-1,0.9,0.1,0.);

  if (fTPCParam)
     fTPCParam->Write(fTPCParam->GetTitle());

  fHitType=1;
}

//_____________________________________________________________________________
void AliTPCv1::CreateGeometry()
{
  //-------------------------------------------------------------
  // Creates geometry for Time Projection Chamber version 1
  // Detailed geometry -- Fast simulation - space points
  // at the pad-rows centers
  // Origin M.Kowalski, INP Cracow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------

  //Begin_Html
  /*
    <img src="picts/AliTPC.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTPCv1Tree.gif">
  */
  //End_Html
  Float_t dm[50];

  Int_t *idtmed = fIdtmed->GetArray(); // TPC media

  Int_t idrotm[120]; // rotation matrices

  Int_t nRotMat = 0; // actual number of rotation matrices

  //
  //  Mother volume TPC (Air) - all volumes will be positioned in it
  //

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=8.;

  //

  dm[3]=-283.7;
  dm[4]=77.017;
  dm[5]=278.;

  //

  dm[6]=-253.6;
  dm[7]=65.6;
  dm[8]=278.;

  //

  dm[9]=-73.3;
  dm[10]=60.9;
  dm[11]=278.;  

  //

  dm[12]=-73.3;
  dm[13]=56.9;
  dm[14]=278.;

  //

  dm[15]=73.3;
  dm[16]=56.9;
  dm[17]=278.;

  //

  dm[18]=73.3;
  dm[19]=60.9;
  dm[20]=278.;

  //

  dm[21]=253.6;
  dm[22]=65.6;
  dm[23]=278.;

  //

  dm[24]=283.7;
  dm[25]=77.017;
  dm[26]=278.;
  
  gMC->Gsvolu("TPC ","PCON",idtmed[0],dm,27);

  // outer part

  //-------------------------------------------------------------------
  //   Tpc Outer INsulator (CO2) - contains cont. vessel and field cage
  //-------------------------------------------------------------------

  dm[0]= 0.;
  dm[1]= 360.;
  dm[2]= 6.;

  //

  dm[3]=-253.6;
  dm[4]=258.;
  dm[5]=275.5;

  //

  dm[6]=-250.6;
  dm[7]=258.;
  dm[8]=275.5; 

  //

  dm[9]=-250.6;
  dm[10]=258.;
  dm[11]=278.;

  //

  dm[12]=253.6;
  dm[13]=258.;
  dm[14]=278.; 

  //

  dm[15]=253.6;
  dm[16]=264.8;
  dm[17]=278.;  

  //

  dm[18]=256.6;
  dm[19]=264.8;
  dm[20]=278.;

  gMC->Gsvolu("TOIN","PCON",idtmed[3],dm,21);

  //----------------------------------------------------------------
  // Tpc Outer Contaiment Vessel  
  //  mother volume - Al, daughters - composite (sandwich)
  //----------------------------------------------------------------

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=6.;

  //

  dm[3]=-250.6;
  dm[4]=270.4;
  dm[5]=278.;

  //

  dm[6]=-247.6;
  dm[7]=270.4;
  dm[8]=278.; 

  //

  dm[9]=-247.6;
  dm[10]=274.8124;
  dm[11]=278.;

  //

  dm[12]=253.6;
  dm[13]=274.8124;
  dm[14]=278.;

  //

  dm[15]=253.6;
  dm[16]=264.8;
  dm[17]=278.;

  //

  dm[18]=256.6;
  dm[19]=264.8;
  dm[20]=278.;

  gMC->Gsvolu("TOCV","PCON",idtmed[4],dm,21);

  // Daughter volumes - sandwich

  // Tpc SAndwich 1 - Al

  dm[0]=274.8124;
  dm[1]=278.;
  dm[2]=252.1;

  gMC->Gsvolu("TSA1","TUBE",idtmed[4],dm,3);

  // Tpc SAndwich 2 - epoxy glue (I use Lexan)

  dm[0] += 5.e-3;
  dm[1] -= 5.e-3;

  gMC->Gsvolu("TSA2","TUBE",idtmed[14],dm,3);

  // Tpc SAndwich 3 - Tedlar

  dm[0] += 0.01;
  dm[1] -= 0.01;
  
  gMC->Gsvolu("TSA3","TUBE",idtmed[9],dm,3);


  // Tpc SAndwich 4 - fiber glass (G10)

  dm[0] += 3.8e-3;
  dm[1] -= 3.8e-3;

  gMC->Gsvolu("TSA4","TUBE",idtmed[12],dm,3);  

  // Tpc SAndwich 5 - NOMEX honeycomb

  dm[0] += 0.075;
  dm[1] -= 0.075;   
  
  gMC->Gsvolu("TSA5","TUBE",idtmed[6],dm,3);

  // 5->4->3->2->1->TCOV


  gMC->Gspos("TSA5",1,"TSA4",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TSA4",1,"TSA3",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TSA3",1,"TSA2",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TSA2",1,"TSA1",0.,0.,0.,0,"ONLY");  

  gMC->Gspos("TSA1",1,"TOCV",0.,0.,3.,0,"ONLY");

  // TCOV-> TOIN

  gMC->Gspos("TOCV",1,"TOIN",0.,0.,0.,0,"ONLY");

  //-------------------------------------------------------
  //  Tpc Outer Field Cage
  //  mother volume - Al, daughters - composite (sandwich)
  //-------------------------------------------------------

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=6.;

  //

  dm[3]=-253.6;
  dm[4]=258.;
  dm[5]=275.5;

  //

  dm[6]=-250.6;
  dm[7]=258.;
  dm[8]=275.5;

  //

  dm[9]=-250.6;
  dm[10]=258.;
  dm[11]=260.0476;

  //

  dm[12]=250.6;
  dm[13]=258.;
  dm[14]=260.0476;

  //

  dm[15]=250.6;
  dm[16]=258.;
  dm[17]=269.6;

  //

  dm[18]=253.6;
  dm[19]=258.;
  dm[20]=269.6;

  gMC->Gsvolu("TOFC","PCON",idtmed[4],dm,21);

  // Daughter volumes 

  // Tpc SAndwich 6 - Tedlar

  dm[0]= 258.;
  dm[1]= 260.0476;
  dm[2]= 252.1;

  gMC->Gsvolu("TSA6","TUBE",idtmed[9],dm,3);

  // Tpc SAndwich 7 - fiber glass

  dm[0] += 3.8e-3;
  dm[1] -= 3.8e-3;

  gMC->Gsvolu("TSA7","TUBE",idtmed[12],dm,3);


  // Tpc SAndwich 8 - NOMEX

  dm[0] += 0.02;
  dm[1] -= 0.02;

  gMC->Gsvolu("TSA8","TUBE",idtmed[6],dm,3);    

  // 8->7->6->TOFC

  gMC->Gspos("TSA8",1,"TSA7",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TSA7",1,"TSA6",0.,0.,0.,0,"ONLY"); 
  gMC->Gspos("TSA6",1,"TOFC",0.,0.,0.,0,"ONLY");

  // TOFC->TOIN

  gMC->Gspos("TOFC",1,"TOIN",0.,0.,0.,0,"ONLY");

  // TOIN->TPC

  gMC->Gspos("TOIN",1,"TPC ",0.,0.,0.,0,"ONLY");

  // inner part

  //--------------------------------------------------------------------
  // Tpc Inner INsulator (CO2) - inner f.c. will be placed there
  // Inner containment vessel will be placed directly in the TPC
  //-------------------------------------------------------------------- 

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=4.;

  // 

  dm[3]=-253.6;
  dm[4]=65.9;
  dm[5]=79.2;

  //

  dm[6]=-73.3;
  dm[7]=61.2;
  dm[8]=79.2;  

  //

  dm[9]=73.3;
  dm[10]=61.2;
  dm[11]=79.2;

  //

  dm[12]=253.6;
  dm[13]=65.9;
  dm[14]=79.2;

  gMC->Gsvolu("TIIN","PCON",idtmed[3],dm,15);

  // the middle part of the F.C. is thinner - carve out the strip - Ne-CO2

  dm[0]=79.16;
  dm[1]=79.2;
  dm[2]=88.;

  gMC->Gsvolu("TII1","TUBE",idtmed[1],dm,3);

  gMC->Gspos("TII1",1,"TIIN",0.,0.,0.,0,"ONLY");

  //-----------------------------------------------------
  // Tpc Inner Field Cage
  // mother volume - Al, daughters - composite (sandwich)
  //------------------------------------------------------

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=10.;

  //

  dm[3]=-253.6;
  dm[4]=70.3;
  dm[5]=79.2;

  //

  dm[6]=-250.6;
  dm[7]=70.3;
  dm[8]=79.2;

  //

  dm[9]=-250.6;
  dm[10]=77.0524;
  dm[11]=79.2;

  //

  dm[12]=-88.;
  dm[13]=77.0524;
  dm[14]=79.2;

  //

  dm[15]=-88.;
  dm[16]=77.0924;
  dm[17]=79.16;

  //

  dm[18]=88.;
  dm[19]=77.0924;
  dm[20]=79.16;

  //

  dm[21]=88.;
  dm[22]=77.0524;
  dm[23]=79.2;

  //

  dm[24]=250.6;
  dm[25]=77.0524;
  dm[26]=79.2;

  //

  dm[27]=250.6;
  dm[28]=70.3;
  dm[29]=79.2;

  //

  dm[30]=253.6;
  dm[31]=70.3;
  dm[32]=79.2;

  gMC->Gsvolu("TIFC","PCON",idtmed[4],dm,33);

  // daughter volumes - central part

  // Tpc Sandwich 9 -Tedlar

  dm[0]=77.0924;
  dm[1]=79.16;
  dm[2]=88.;

  gMC->Gsvolu("TSA9","TUBE",idtmed[9],dm,3); 

  // Tpc Sandwich 10 - fiber glass (G10) 

  dm[0] += 3.8e-3;
  dm[1] -= 3.8e-3;

  gMC->Gsvolu("TS10","TUBE",idtmed[12],dm,3);

  // Tpc Sandwich 11 - NOMEX

  dm[0] += 0.03;
  dm[1] -= 0.03; 

  gMC->Gsvolu("TS11","TUBE",idtmed[6],dm,3);

  // 11->10->9->TIFC

  gMC->Gspos("TS11",1,"TS10",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS10",1,"TSA9",0.,0.,0.,0,"ONLY");

  gMC->Gspos("TSA9",1,"TIFC",0.,0.,0.,0,"ONLY");

  // daughter volumes - outer parts (reinforced)

  // Tpc Sandwich 12 -Tedlar

  dm[0]=77.0524;
  dm[1]=79.2;
  dm[2]=82.05;

  gMC->Gsvolu("TS12","TUBE",idtmed[9],dm,3);

  // Tpc Sandwich 13 - fiber glass (G10) 

  dm[0] += 3.8e-3;
  dm[1] -= 3.8e-3;

  gMC->Gsvolu("TS13","TUBE",idtmed[12],dm,3);

  // Tpc Sandwich 14 - NOMEX

  dm[0] += 0.07;
  dm[1] -= 0.07;  

  gMC->Gsvolu("TS14","TUBE",idtmed[6],dm,3);

  // 14->13->12->TIFC

  gMC->Gspos("TS14",1,"TS13",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS13",1,"TS12",0.,0.,0.,0,"ONLY"); 

  gMC->Gspos("TS12",1,"TIFC",0.,0.,170.05,0,"ONLY");
  gMC->Gspos("TS12",2,"TIFC",0.,0.,-170.05,0,"ONLY"); 

  // place this inside the inner insulator

  gMC->Gspos("TIFC",1,"TIIN",0.,0.,0.,0,"ONLY");

  // and now in the TPC...

  gMC->Gspos("TIIN",1,"TPC ",0.,0.,0.,0,"ONLY");

  //---------------------------------------------------------
  // Tpc Inner Containment vessel - Cones
  //---------------------------------------------------------

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=8.;

  //

  dm[3]=71.8;
  dm[4]=56.9;
  dm[5]=59.4;

  //

  dm[6]=73.;
  dm[7]=56.9;
  dm[8]=59.4;

  //

  dm[9]=73.;
  dm[10]=56.9;
  dm[11]=61.2;

  //

  dm[12]=73.3;
  dm[13]=56.9;
  dm[14]=61.2;

  //
   
  dm[15]=73.3;
  dm[16]=60.9;
  dm[17]=61.2;

  // 

  dm[18]=253.6;
  dm[19]=65.6;
  dm[20]=65.9; 

  //

  dm[21]=253.6;
  dm[22]=65.6;
  dm[23]=74.6;

  //

  dm[24]=256.6;
  dm[25]=65.6;
  dm[26]=74.6;

  gMC->Gsvolu("TICC","PCON",idtmed[4],dm,27);

  Float_t phi1,phi2,phi3,theta1,theta2,theta3; // rotation angles

  // reflection matrix
  
  theta1 = 90.;
  phi1   = 0.;
  theta2 = 90.;
  phi2   = 270.;
  theta3 = 180.;
  phi3   = 0.;

  AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);

  gMC->Gspos("TICC",1,"TPC ",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TICC",2,"TPC ",0.,0.,0.,idrotm[nRotMat],"ONLY");


  //---------------------------------------------------------
  // Tpc Inner Containment vessel - Middle part -Al
  //---------------------------------------------------------

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=6.;

  //

  dm[3]=-71.6;
  dm[4]=60.2;
  dm[5]=61.2;

  //

  dm[6]=-69.1;
  dm[7]=60.2;
  dm[8]=61.2;

  //

  dm[9]=-69.1;
  dm[10]=60.6224;
  dm[11]=61.2;  

  //

  dm[12]=69.1;
  dm[13]=60.6224;
  dm[14]=61.2;

  //

  dm[15]=69.1;
  dm[16]=60.2;
  dm[17]=61.2;

  //

  dm[18]=71.6;
  dm[19]=60.2;
  dm[20]=61.2;

  gMC->Gsvolu("TICM","PCON",idtmed[4],dm,21);

  // Tpc Sandwich 15 - Al

  dm[0]=60.6224;
  dm[1]=61.2;
  dm[2]=70.1;

  gMC->Gsvolu("TS15","TUBE",idtmed[4],dm,3);

  // Tpc Sandwich 16 -  epoxy glue

  dm[0] += 5.e-3;
  dm[1] -= 5.e-3;

  gMC->Gsvolu("TS16","TUBE",idtmed[14],dm,3);

  // Tpc Sandwich 17 - Tedlar

  dm[0] += 0.01;
  dm[1] -= 0.01;

  gMC->Gsvolu("TS17","TUBE",idtmed[9],dm,3);

  // Tpc Sandwich 18 - carbon fiber

  dm[0] += 3.8e-3;
  dm[1] -= 3.8e-3;

  gMC->Gsvolu("TS18","TUBE",idtmed[15],dm,3);  

  // Tpc Sandwich 19 - Nomex

  dm[0] += 0.02;
  dm[1] -= 0.02;

  gMC->Gsvolu("TS19","TUBE",idtmed[6],dm,3); 

  // 19->18->17->16->15-> TICM

  gMC->Gspos("TS19",1,"TS18",0.,0.,0.,0,"ONLY"); 
  gMC->Gspos("TS18",1,"TS17",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS17",1,"TS16",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TS16",1,"TS15",0.,0.,0.,0,"ONLY");

  gMC->Gspos("TS15",1,"TICM ",0.,0.,0.,0,"ONLY");
 

  // TPc inner cont. vessel Joints

  dm[0]=60.2;
  dm[1]=61.2;
  dm[2]=0.5;

  gMC->Gsvolu("TPJ1","TUBE",idtmed[4],dm,3);

  gMC->Gspos("TPJ1",1,"TPC ",0.,0.,72.3,0,"ONLY");
  gMC->Gspos("TPJ1",2,"TPC ",0.,0.,-72.3,0,"ONLY");

  //

  dm[0]=0.;
  dm[1]=360.;
  dm[2]=4.;

  //

  dm[3]=70.8;
  dm[4]=58.4;
  dm[5]=60.1;

  //

  dm[6]=71.2;
  dm[7]=58.4;
  dm[8]=60.1;

  //

  dm[9]=71.2;
  dm[10]=58.4;
  dm[11]=59.4;

  //

  dm[12]=71.6;
  dm[13]=58.4;
  dm[14]=59.4;

  gMC->Gsvolu("TPJ2","PCON",idtmed[4],dm,15);

  gMC->Gspos("TPJ2",1,"TPC ",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TPJ2",2,"TPC ",0.,0.,0.,idrotm[nRotMat],"ONLY");



  // Tpc Inner Containment vessel Seal (Viton, I use Lexan for a time being)

  dm[0]=58.4;
  dm[1]=61.2;
  dm[2]=0.1;

  gMC->Gsvolu("TICS","TUBE",idtmed[14],dm,3);

  gMC->Gspos("TICS",1,"TPC ",0.,0.,71.7,0,"ONLY");
  gMC->Gspos("TICS",2,"TPC ",0.,0.,-71.7,0,"ONLY"); 

  // TICM -> TPC

  gMC->Gspos("TICM",1,"TPC ",0.,0.,0.,0,"ONLY");

  //

  nRotMat++; // prepare for the next rotation matrix 

  //---------------------------------------------------------
  //  Tpc Dift Gas volume Nonsensitive (Ne-CO2 90/10)
  //  and its daughters (HV membrane, rods, readout chambers)
  //---------------------------------------------------------

  dm[0]= 79.2;
  dm[1]= 258.0;
  dm[2]= 253.6;

  gMC->Gsvolu("TDGN","TUBE",idtmed[1],dm,3); 

  // sector opening angles

  Float_t innerOpenAngle = fTPCParam->GetInnerAngle();
  Float_t outerOpenAngle = fTPCParam->GetOuterAngle();

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

  gMC->Gsvolu("TIC3","TRD1",idtmed[1],dm,4);

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

    dm[2] = 1.2;

    gMC->Gsposp("TRCE",nsLow+5,"TIC1",0.,1.05,zz-2.1,0,"ONLY",dm,4);

    rr += 0.4;
    dm[0] = rr*tanAlpha - space;
    zz += (0.4+9.525); 

  }

  dm[0]= 12.472;
  // dm[1] - this is the dm[1] from the previous TRCE
  dm[2]= 1.05;
  dm[3]= 19.65;

  gMC->Gsposp("TRCE",9,"TIC1",0.,-1.2,0.,0,"ONLY",dm,4);   

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
  // Tpc Outer readout chamber Pad Plane
  //------------------------------------------------

  dm[2]= 0.5;

  gMC->Gsvolu("TOPP","TRD1",idtmed[12],dm,4);

  dm[0] -= 1.218511934;
  dm[1] -= 1.218511934;
  dm[2] = 0.35;

  gMC->Gsvolu("TOC3","TRD1",idtmed[1],dm,4);

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

  dm[0]= 74.9;
  dm[1]= 264.4;
  dm[2]= 3.0;

  gMC->Gsvolu("TPMW","TUBE",idtmed[4],dm,3);

  //--------------------------------------------------------------------------
  //  Tpc Empty Space for the Readout chambers
  //--------------------------------------------------------------------------  

  Float_t rLow= 86.2;
  Float_t rUp= 243.5;
  Float_t dR = 0.5*(rUp-rLow);

  space= 1.5/cosAlpha; // wheel ribs are 3.0 cm wide

  dm[0]= rLow*tanAlpha-space;
  dm[1]= rUp*tanAlpha-space;
  dm[2] = 3.0;
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
  dm[2]= 2.25;
  dm[3]= 1.275;

  zBar = -dR+dm[3];

  gMC->Gsposp("TESB",1,"TESR",0.,0.75,zBar,0,"ONLY",dm,4);

  // middle bar

  dm[0]= 131.65*tanAlpha-space;
  dm[1]= 136.5*tanAlpha-space;
  dm[3]= 2.425;

  zBar = -dR +131.65+dm[3]-rLow;

  gMC->Gsposp("TESB",2,"TESR",0.,0.75,zBar,0,"ONLY",dm,4);  

  // upper bar

  dm[0]= 240.4*tanAlpha-space;
  dm[1]= rUp*tanAlpha-space;
  dm[3]= 1.55;

  zBar = dR-dm[3];

  gMC->Gsposp("TESB",3,"TESR",0.,0.75,zBar,0,"ONLY",dm,4);


  //------------------------------------------------------
  //  TPc Lower "S" Sectors 
  //------------------------------------------------------
  
  Float_t deadSpace;


  Float_t inSecLowEdge = fTPCParam->GetInnerRadiusLow();
  Float_t inSecUpEdge =  fTPCParam->GetInnerRadiusUp();

  dm[0] = inSecLowEdge*TMath::Tan(0.5*innerOpenAngle)-0.01;
  dm[1] = inSecUpEdge*TMath::Tan(0.5*innerOpenAngle)-0.01;
  dm[2] = 0.5*(250. - 5.e-3);
  dm[3] = 0.5*(inSecUpEdge-inSecLowEdge);  

  gMC->Gsvolu("TPLS", "TRD1", idtmed[1], dm, 4); // nonsensitive 

  //----------------------------------------------------------
  //  TPc Upper Sectors
  //----------------------------------------------------------

  Float_t ouSecLowEdge = fTPCParam->GetOuterRadiusLow();
  Float_t ouSecUpEdge = fTPCParam->GetOuterRadiusUp();

  dm[0] = ouSecLowEdge*TMath::Tan(0.5*outerOpenAngle)-0.01;  
  dm[1] = ouSecUpEdge*TMath::Tan(0.5*outerOpenAngle)-0.01;  
  dm[2] = 0.5*(250. - 5.e-3);
  dm[3] = 0.5*(ouSecUpEdge-ouSecLowEdge);

  gMC->Gsvolu("TPUS", "TRD1", idtmed[1], dm, 4); // nonsensitive 

  //--------------------------------------------------------
  // TPc Sensitive Strips, the same for S and L sectors
  //--------------------------------------------------------

  gMC->Gsvolu("TPSS","TRD1",idtmed[2],dm,0); // sensitive

  Int_t nofStrips,nstr;
  Float_t r1,r2,zs;
  Float_t stripThick = 0.01; // 100 microns

  // inner sector

  nofStrips = fTPCParam->GetNRowLow();
  deadSpace = fTPCParam->GetInnerWireMount();

  dm[2] = 0.5*(250. - 5.e-3);
  dm[3] = 0.5 * stripThick;

  for(nstr=0;nstr<nofStrips;nstr++){

    r1 = fTPCParam->GetPadRowRadiiLow(nstr);
    r2 = r1 + stripThick;     
    dm[0] = r1 * TMath::Tan(0.5*innerOpenAngle) - deadSpace;
    dm[1] = r2 * TMath::Tan(0.5*innerOpenAngle) - deadSpace;
    zs = -inSecLowEdge -0.5*(inSecUpEdge-inSecLowEdge);
    zs += r1;
    zs += dm[3];
    
    gMC->Gsposp("TPSS", nstr+1, "TPLS", 0., 0., zs, 0, "ONLY", dm, 4);

  }

  Int_t nsSave = nofStrips;

  // outer sector

  nofStrips = fTPCParam->GetNRowUp();
  deadSpace = fTPCParam->GetOuterWireMount();

  dm[2] = 0.5*(250. - 5.e-3);
  dm[3] = 0.5 * stripThick;
  
   for(nstr=0;nstr<nofStrips;nstr++){
    
    r1 = fTPCParam->GetPadRowRadiiUp(nstr); 
    r2 = r1 + stripThick;
    dm[0] = r1 * TMath::Tan(0.5*outerOpenAngle) - deadSpace;
    dm[1] = r2 * TMath::Tan(0.5*outerOpenAngle) - deadSpace;
    zs = -ouSecLowEdge -0.5*(ouSecUpEdge-ouSecLowEdge);
    zs += r1;
    zs += dm[3];

    gMC->Gsposp("TPSS", nstr+1+nsSave, "TPUS", 0., 0., zs, 0, "ONLY", dm, 4);

   } 

  //-------------------------------------------------------
  //  positioning of the empty spaces into the main wheel
  //  and readout chambers and sectors into the drift gas
  //-------------------------------------------------------


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


    // readout chambers->TDGN (drift gas)

    xc = rInner*TMath::Cos(phi3*kDegrad);
    yc = rInner*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TIRC",ns+1,"TDGN",xc,yc,252.,idrotm[nRotMat],"ONLY");

    // here lower sectors 

    gMC->Gspos("TPLS",ns+1,"TDGN",xc,yc,125.0025,idrotm[nRotMat],"ONLY");
    gMC->Gspos("TPLS",ns+nInnerSector+1,"TDGN",xc,yc,-125.0025,idrotm[nRotMat],"ONLY");


     
    xc = rOuter*TMath::Cos(phi3*kDegrad);
    yc = rOuter*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TORC",ns+1,"TDGN",xc,yc,252.1,idrotm[nRotMat],"ONLY");

    // here upper sectors 

    gMC->Gspos("TPUS",ns+1,"TDGN",xc,yc,125.0025,idrotm[nRotMat],"ONLY");
    gMC->Gspos("TPUS",ns+nOuterSector+1,"TDGN",xc,yc,-125.0025,idrotm[nRotMat],"ONLY");


    nRotMat++;

    theta2 = 0.; // reflection

    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);

    xc = rInner*TMath::Cos(phi3*kDegrad);
    yc = rInner*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TIRC",ns+nInnerSector+1,"TDGN",xc,yc,-252.,idrotm[nRotMat],"ONLY");

    xc = rOuter*TMath::Cos(phi3*kDegrad);
    yc = rOuter*TMath::Sin(phi3*kDegrad);

    gMC->Gspos("TORC",ns+nOuterSector+1,"TDGN",xc,yc,-252.1,idrotm[nRotMat],"ONLY");

    nRotMat++;

  } 
  // TPMW->TPC

  gMC->Gspos("TPMW",1,"TPC ",0.,0.,256.6,0,"ONLY");
  gMC->Gspos("TPMW",2,"TPC ",0.,0.,-256.6,idrotm[0],"ONLY");

  //---------------------------------------------------------
  //  Tpc High Voltage Membrane - 100 microns of mylar
  //---------------------------------------------------------

  dm[0]=82.8;
  dm[1]=252.;
  dm[2]=0.005;

  gMC->Gsvolu("THVM","TUBE",idtmed[8],dm,3);

  gMC->Gspos("THVM",1,"TDGN",0.,0.,0.,0,"ONLY");

  // Tpc High Voltage membrane Holders

  gMC->Gsvolu("THVH","TUBE",idtmed[4],dm,0);

  

  // inner

  dm[0]=79.3;
  dm[1]=82.8;
  dm[2]=0.2;

  gMC->Gsposp("THVH",1,"TDGN",0.,0.,0.,0,"ONLY",dm,3);
  
  // outer

  dm[0]= 252.;
  dm[1]= 257.9;
  dm[2]= 0.4;

  gMC->Gsposp("THVH",2,"TDGN",0.,0.,0.,0,"ONLY",dm,3);

  //----------------------------------------------------------
  // TPc Support Rods - MAKROLON
  //----------------------------------------------------------

  dm[0]= 0.9;
  dm[1]= 1.2;

  gMC->Gsvolu("TPSR","TUBE",idtmed[7],dm,0); // inner and outer rods differ


  for(Int_t nrod=0;nrod<18;nrod++){
    Float_t angle=innerOpenAngle*(Float_t)nrod;

    xc=81.5*TMath::Cos(angle);
    yc=81.5*TMath::Sin(angle); 

    dm[2]=126.7;

    gMC->Gsposp("TPSR",nrod+1,"TDGN",xc,yc,126.9,0,"ONLY",dm,3); 
    gMC->Gsposp("TPSR",nrod+19,"TDGN",xc,yc,-126.9,0,"ONLY",dm,3);

    dm[2]=126.6;

    xc=254.25*TMath::Cos(angle);
    yc=254.25*TMath::Sin(angle);   
      
    // rod number 54 contans the HV cable

    if(nrod<17) {
      gMC->Gsposp("TPSR",nrod+37,"TDGN",xc,yc,127.,0,"ONLY",dm,3);
      gMC->Gsposp("TPSR",nrod+54,"TDGN",xc,yc,-127.,0,"ONLY",dm,3);
    }
    
  }

  //----------------------------------------------------------
  // Tpc High Voltage Rod - MAKROLON + Copper cable
  //----------------------------------------------------------

  // rod with cable (Left)

  dm[0]=0.;
  dm[1]=2.25;
  dm[2]=126.6;

  gMC->Gsvolu("THVL","TUBE",idtmed[7],dm,3);

  // HV cable
 
  dm[0]=0.;
  dm[1]=0.3;
  dm[2]=126.6;

  gMC->Gsvolu("THVC","TUBE",idtmed[10],dm,3);  

  // empty space

  dm[0]=0.3;
  dm[1]=1.;
  dm[2]=126.6;

  gMC->Gsvolu("THVE","TUBE",idtmed[1],dm,3);

  gMC->Gspos("THVC",1,"THVL",0.,0.,0.,0,"ONLY");
  gMC->Gspos("THVE",1,"THVL",0.,0.,0.,0,"ONLY");

  // rod without cable

  dm[0]=1.8;
  dm[1]=2.25;
  dm[2]=126.6;

  gMC->Gsvolu("THVR","TUBE",idtmed[7],dm,3);

  gMC->Gspos("THVL",1,"TDGN",xc,yc,-127.,0,"ONLY");  
  gMC->Gspos("THVR",1,"TDGN",xc,yc,127.,0,"ONLY");

  gMC->Gspos("TDGN",1,"TPC ",0.,0.,0.,0,"ONLY"); 
 
  // services wheel cover -> wheel


  gMC->Gspos("TSWC",1,"TSSW",0.,0.,4.5,0,"ONLY");
  gMC->Gspos("TSWC",2,"TSSW",0.,0.,-4.5,0,"ONLY");


  // put the wheel into the TPC

  gMC->Gspos("TSSW",1,"TPC ",0.,0.,278.7,0,"ONLY");
  gMC->Gspos("TSSW",2,"TPC ",0.,0.,-278.7,0,"ONLY");

  //

  gMC->Gsord("TPMW",6);
  gMC->Gsord("TSSW",6);
  gMC->Gsord("TSWC",6);
  gMC->Gsord("TPLS",3);
  gMC->Gsord("TPUS",3);
  gMC->Gsord("TDGN",6);


  // put the TPC into ALIC (main mother volume)

  gMC->Gspos("TPC ",1,"ALIC",0.,0.,0.,0,"ONLY");



} // end of function

//_____________________________________________________________________________
void AliTPCv1::DrawDetector()
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
  gMC->Gsatt("TPC ","SEEN",0);
  gMC->Gsatt("TOIN","SEEN",1);
  gMC->Gsatt("TOIN","COLO",7);
  gMC->Gsatt("TOCV","SEEN",1);
  gMC->Gsatt("TOCV","COLO",4);
  gMC->Gsatt("TSA1","SEEN",0);
  gMC->Gsatt("TSA2","SEEN",0);
  gMC->Gsatt("TSA3","SEEN",0);
  gMC->Gsatt("TSA4","SEEN",0);  
  gMC->Gsatt("TSA5","SEEN",0);
  gMC->Gsatt("TOFC","SEEN",1);
  gMC->Gsatt("TOFC","COLO",4);
  gMC->Gsatt("TSA6","SEEN",0);
  gMC->Gsatt("TSA7","SEEN",0);
  gMC->Gsatt("TSA8","SEEN",0);    
  gMC->Gsatt("TIIN","SEEN",1);
  gMC->Gsatt("TIIN","COLO",7);
  gMC->Gsatt("TII1","SEEN",0);
  gMC->Gsatt("TIFC","SEEN",1);
  gMC->Gsatt("TIFC","COLO",4);
  gMC->Gsatt("TSA9","SEEN",0); 
  gMC->Gsatt("TS10","SEEN",0);
  gMC->Gsatt("TS11","SEEN",0);
  gMC->Gsatt("TS12","SEEN",0);
  gMC->Gsatt("TS13","SEEN",0);
  gMC->Gsatt("TS14","SEEN",0);
  gMC->Gsatt("TICC","SEEN",0);
  gMC->Gsatt("TICM","SEEN",0);
  gMC->Gsatt("TS15","SEEN",0);
  gMC->Gsatt("TS16","SEEN",0);
  gMC->Gsatt("TS17","SEEN",0);
  gMC->Gsatt("TS18","SEEN",0);  
  gMC->Gsatt("TS19","SEEN",0); 
  gMC->Gsatt("TPJ1","SEEN",0);
  gMC->Gsatt("TPJ2","SEEN",0);
  gMC->Gsatt("TICS","SEEN",0);
  gMC->Gsatt("TDGN","SEEN",0); 
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
  gMC->Gsatt("TSSW","COLO",3);
  gMC->Gsatt("TSWC","COLO",3);
  gMC->Gsatt("TSCE","COLO",6);
  gMC->Gsatt("TSCE","SEEN",1);
  gMC->Gsatt("TWES","SEEN",0);
  gMC->Gsatt("TSWB","SEEN",0);
  gMC->Gsatt("TPEL","SEEN",0);
  gMC->Gsatt("TPMW","SEEN",1);
  gMC->Gsatt("TESR","SEEN",1);
  gMC->Gsatt("TPMW","COLO",12);
  gMC->Gsatt("TIC1","COLO",5);
  gMC->Gsatt("TOC1","COLO",5);
  gMC->Gsatt("TESB","SEEN",0);
  gMC->Gsatt("THVM","SEEN",1);
  gMC->Gsatt("THVM","COLO",11);
  gMC->Gsatt("THVH","SEEN",0);
  gMC->Gsatt("TPSR","SEEN",0); 
  gMC->Gsatt("THVL","SEEN",0);
  gMC->Gsatt("THVC","SEEN",0);  
  gMC->Gsatt("THVE","SEEN",0);
  gMC->Gsatt("THVR","SEEN",0);
  gMC->Gsatt("TPSS","SEEN",0);
  gMC->Gsatt("TPUS","SEEN",0);
  gMC->Gsatt("TPLS","SEEN",0);

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
void AliTPCv1::CreateMaterials()
{
  //
  // Define materials for Time Projection Chamber
  //
  AliTPC::CreateMaterials();
}

//_____________________________________________________________________________
void AliTPCv1::Init()
{
  //
  // Initialises TPC detector after it has been created
  //

  fIdSens=gMC->VolId("TPSS"); // sensitive strip
  fIdLSec=gMC->VolId("TPLS"); // lower sector
  fIdUSec=gMC->VolId("TPUS"); // upper sector

  printf("%s: *** TPC version 1 initialized***\n",ClassName());

  // printf("TPC version 1 initialized\n");
}

//_____________________________________________________________________________
void AliTPCv1::StepManager()
{
  //
  // Called at every step in the Time Projection Chamber
  //
  Int_t         copy, id, i;
  Float_t       hits[5];
  Int_t         vol[2];
  TLorentzVector p;
 //

  if(gMC->TrackCharge() && gMC->IsTrackEntering()) {

    // 
    // Only charged tracks entering the sensitive volume
    //

    id = gMC->CurrentVolID(copy);

    if(id == fIdSens){

      // sensitive volume (strip)

      vol[1]= copy-1; // row number (absolute)
    
      // sector type
      
      id = gMC->CurrentVolOffID(1,copy); 
      
      if(id == fIdLSec){

        // lower sector
     
        vol[0] = copy-1; // sector number

      }
      else {
   
        // upper sector

        vol[0] = copy-1+fTPCParam->GetNInnerSector(); // sector number
        vol[1] -= fTPCParam->GetNRowLow(); // row number  

      } 
    
      if(vol[1] == 0){
  
        // because Jouri wants to have this

        gMC->TrackMomentum(p);
        hits[0]=p[0];
        hits[1]=p[1];
        hits[2]=p[2];
        hits[3]=0.; // this hit has no energy loss
	// new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->GetCurrentTrackNumber(),vol,hits);

	// Get also the track time for pileup simulation
	hits[4]=gMC->TrackTime();

        AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits); // M.I. 

        gMC->TrackPosition(p);
        hits[0]=p[0];
        hits[1]=p[1];
        hits[2]=p[2];
        hits[3]=0.; // this hit has no energy loss
        // new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->GetCurrentTrackNumber(),vol,hits);

	// Get also the track time for pileup simulation
	hits[4]=gMC->TrackTime();

	AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits); // M.I. 

      } 

      gMC->TrackPosition(p);
      for(i=0;i<3;++i) hits[i]=p[i];
      hits[3]=1; //I'd like to have something positive here (I.Belikov)
      // new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->GetCurrentTrackNumber(),vol,hits);

      // Get also the track time for pileup simulation
      hits[4]=gMC->TrackTime();

      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits); // M.I. 

    }

  }  
   
}
//_____________________________________________________________________________



