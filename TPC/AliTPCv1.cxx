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
Revision 1.17.2.4  2000/06/26 07:39:42  kowal2
Changes to obey the coding rules

Revision 1.17.2.3  2000/06/25 08:38:41  kowal2
Splitted from AliTPCtracking

Revision 1.17.2.2  2000/06/16 12:58:13  kowal2
Changed parameter settings

Revision 1.17.2.1  2000/06/09 07:15:07  kowal2

Defaults loaded automatically (hard-wired)
Optional parameters can be set via macro called in the constructor

Revision 1.17  2000/05/15 10:00:30  kowal2
Corrected bug in the TPC geometry, thanks to Ivana Hrivnacova

Revision 1.16  2000/04/17 09:37:33  kowal2
removed obsolete AliTPCDigitsDisplay.C

Revision 1.15.8.2  2000/04/10 08:29:48  kowal2

Different pad geometry for different sectors
Updated readouch chambers geometry

Revision 1.15.8.1  2000/04/10 07:56:53  kowal2
Not used anymore - removed

Revision 1.15  1999/11/04 17:28:07  fca
Correct barrel part of HV Degrader

Revision 1.14  1999/10/08 06:27:23  fca
Corrected bug in the HV degrader geometry, thanks to G.Tabary

Revision 1.13  1999/10/04 13:39:54  fca
Correct array index problem

Revision 1.12  1999/09/29 09:24:34  fca
Introduction of the Copyright and cvs Log

*/

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

#include <TMath.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TBRIK.h>
#include <TTUBE.h>
#include "AliTPCv1.h"
#include "AliRun.h"
#include <iostream.h>
#include <fstream.h>

#include "AliMC.h"
#include "AliConst.h"

#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCDigitsArray.h"
#include <TInterpreter.h>

ClassImp(AliTPCv1)
 
//_____________________________________________________________________________
AliTPCv1::AliTPCv1(const char *name, const char *title) 
  :AliTPC(name, title) {
  //
  // Standard constructor for Time Projection Chamber
  //
  fIdSens1=fIdSens2=0;

  SetGasMixt(2,20,10,-1,0.9,0.1,0.);

  if (fTPCParam)
     fTPCParam->Write(fTPCParam->GetTitle());
}

//_____________________________________________________________________________
void AliTPCv1::CreateGeometry()
{
  //
  // Creates geometry for Time Projection Chamber version 1
  // Detailed geometry -- Fast simulation
  // Origin M.Kowalski 
  //
  //Begin_Html
  /*
    <img src="picts/AliTPCv1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTPCv1Tree.gif">
  */
  //End_Html


  Int_t *idtmed = fIdtmed->GetArray();

  Float_t dm[21];
  Int_t idrotm[120];

  Int_t nRotMat = 0;


  // ---------------------------------------------------- 
  //          FIELD CAGE WITH ENDCAPS - G10
  //          THIS IS ALSO A TPC MOTHER VOLUME 
  // ---------------------------------------------------- 

  dm[0] = 76.;
  dm[1] = 278.;
  dm[2] = 275.;

  gMC->Gsvolu("TPC ", "TUBE", idtmed[8], dm, 3); 

  //-----------------------------------------------------
  //  Endcap cover c-fibre 0.86% X0
  //-----------------------------------------------------

  dm[0] = 78.;
  dm[1] = 258.;
  dm[2] = 0.95;

  gMC->Gsvolu("TPEC","TUBE",idtmed[10],dm,3);

  //-----------------------------------------------------
  // Drift gas , leave 2 cm at the outer radius
  // and inner raddius
  //-----------------------------------------------------

  dm[0] = 78.;
  dm[1] = 258.;
  dm[2] = 250.;

  gMC->Gsvolu("TGAS", "TUBE", idtmed[3], dm, 3);


  //------------------------------------------------------
  //  membrane holder - carbon fiber
  //------------------------------------------------------


  gMC->Gsvolu("TPMH","TUBE",idtmed[6],dm,0);

  dm[0] = 252.;
  dm[1] = 258.;
  dm[2] = 0.2;

  gMC->Gsposp("TPMH",1,"TGAS",0.,0.,0.,0,"ONLY",dm,3);
 
  dm[0] = 78.;
  dm[1] = 82.;
  dm[2] = 0.1;

  gMC->Gsposp("TPMH",2,"TGAS",0.,0.,0.,0,"ONLY",dm,3);

  //----------------------------------------------------------
  //  HV membrane - 25 microns of mylar
  //----------------------------------------------------------

  dm[0] = 82.;
  dm[1] = 252.;
  dm[2] = 0.00125;

  gMC->Gsvolu("TPHV","TUBE",idtmed[5],dm,3);

  gMC->Gspos("TPHV",1,"TGAS",0.,0.,0.,0,"ONLY");

  gMC->Gspos("TGAS",1,"TPC ",0.,0.,0.,0,"ONLY");

  //----------------------------------------------------------
  // "side" gas volume, the same as the drift gas
  // the readout chambers are placed there.  
  //----------------------------------------------------------

  dm[0] = 78.;
  dm[1] = 258.;
  dm[2] = 0.5*(275. - 250.);
   
  gMC->Gsvolu("TPSG", "TUBE", idtmed[2], dm, 3);

  Float_t zSide = dm[2]; // 1/2 of the side gas thickness

  //-----------------------------------------------------------
  //   Readout chambers , 25% of X0, I use Al as the material
  //-----------------------------------------------------------

  Float_t innerOpenAngle = fTPCParam->GetInnerAngle();
  Float_t outerOpenAngle = fTPCParam->GetOuterAngle();

  Float_t innerAngleShift = fTPCParam->GetInnerAngleShift();
  Float_t outerAngleShift = fTPCParam->GetOuterAngleShift();


  Int_t nInnerSector = fTPCParam->GetNInnerSector()/2;
  Int_t nOuterSector = fTPCParam->GetNOuterSector()/2;


  Float_t inSecLowEdge = fTPCParam->GetInnerRadiusLow();
  Float_t inSecUpEdge =  fTPCParam->GetInnerRadiusUp();

  Float_t ouSecLowEdge = fTPCParam->GetOuterRadiusLow();
  Float_t ouSecUpEdge = fTPCParam->GetOuterRadiusUp();

  Float_t secThick = 2.225; // Al

  //  S (Inner) sectors

  Float_t lowEdge = fTPCParam->GetInnerFrameSpace();

  dm[0] = inSecLowEdge*TMath::Tan(0.5*innerOpenAngle)-lowEdge;
  dm[1] = inSecUpEdge*TMath::Tan(0.5*innerOpenAngle)-lowEdge;
  dm[2] = 0.5*secThick;
  dm[3] = 0.5*(inSecUpEdge-inSecLowEdge);

  Float_t xCenterS = inSecLowEdge+dm[3];

  gMC->Gsvolu("TRCS", "TRD1", idtmed[0], dm, 4); 

  //  L (Outer) sectors

  Float_t upEdge = fTPCParam->GetOuterFrameSpace();

  dm[0] = ouSecLowEdge*TMath::Tan(0.5*outerOpenAngle)-upEdge;
  dm[1] = ouSecUpEdge*TMath::Tan(0.5*outerOpenAngle)-upEdge;
  dm[2] = 0.5*secThick;
  dm[3] = 0.5*(ouSecUpEdge-ouSecLowEdge);

  Float_t xCenterL = ouSecLowEdge+dm[3];  

  gMC->Gsvolu("TRCL", "TRD1", idtmed[0], dm, 4);

  Float_t z1 = -zSide + secThick*0.5;

  //------------------------------------------------------------------
  // S sectors - "gas sectors" (TRD1)
  //------------------------------------------------------------------

  dm[0] = inSecLowEdge*TMath::Tan(0.5*innerOpenAngle)-0.01;
  dm[1] = inSecUpEdge*TMath::Tan(0.5*innerOpenAngle)-0.01;
  dm[2] = 0.5*(250. - 0.001);
  dm[3] = 0.5*(inSecUpEdge-inSecLowEdge);  

  gMC->Gsvolu("TSGA", "TRD1", idtmed[3], dm, 4); // nonsensitive

  //------------------------------------------------------------------
  // sensitive strips - "pad rows"
  //------------------------------------------------------------------

  Int_t ns;
  Float_t r1,r2,zz;

  Float_t stripThick = 0.01; // 100 microns
  Float_t dead;

  gMC->Gsvolu("TSST", "TRD1", idtmed[4], dm, 0);

  dm[2] = 0.5*(250. - 0.002);
  dm[3] = 0.5 * stripThick;


  // S-sector

  dead = fTPCParam->GetInnerWireMount();

  for (ns = 0; ns < fTPCParam->GetNRowLow(); ns++) {

    r1 = fTPCParam->GetPadRowRadiiLow(ns);
    r2 = r1 + stripThick;     
    dm[0] = r1 * TMath::Tan(0.5*innerOpenAngle) - dead;
    dm[1] = r2 * TMath::Tan(0.5*innerOpenAngle) - dead;

    zz = -inSecLowEdge -0.5*(inSecUpEdge-inSecLowEdge);
    zz += r1;
    zz += dm[3];

    gMC->Gsposp("TSST", ns+1, "TSGA", 0., 0., zz, 0, "ONLY", dm, 4);
    

  }

  //-----------------------------------------------------------------
  //  L sectors - "gas sectors" (PGON to avoid overlaps)
  //-----------------------------------------------------------------

  dm[0] = 360.*kDegrad - 0.5*outerOpenAngle;
  dm[0] *= kRaddeg;
  dm[0] = (Float_t)TMath::Nint(dm[0]);

  dm[1] = outerOpenAngle*kRaddeg;
  dm[1] = (Float_t)TMath::Nint(dm[1]);

  dm[2] = 1.;
  dm[3] = 4.;

  dm[4] = 0.002;
  dm[5] = ouSecLowEdge;
  dm[6] = 252.*TMath::Cos(0.5*outerOpenAngle)-0.002;

  dm[7] = dm[4]+0.2;
  dm[8] = dm[5];
  dm[9] = dm[6];

  dm[10] = dm[7];
  dm[11] = ouSecLowEdge;
  dm[12] = ouSecUpEdge;

  dm[13] = 250.;
  dm[14] = dm[11];
  dm[15] = dm[12];

  gMC->Gsvolu("TLGA","PGON",idtmed[3],dm,16);

  //------------------------------------------------------------------
  // sensitive strips - "pad rows"
  //------------------------------------------------------------------

  Float_t rmax = dm[6];

  // L-sectors

  gMC->Gsvolu("TLST", "PGON", idtmed[4], dm, 0);

  dm[0] = 360.*kDegrad - 0.5*outerOpenAngle;
  dm[0] *= kRaddeg;
  dm[0] = (Float_t)TMath::Nint(dm[0]);

  dm[1] = outerOpenAngle*kRaddeg;
  dm[1] = (Float_t)TMath::Nint(dm[1]);

  dm[2] = 1.;
  dm[3] = 2.;

  dm[7] = 250.;

  dead = fTPCParam->GetOuterWireMount();

  Float_t xx = dead/TMath::Tan(0.5*outerOpenAngle);

  for(ns=0;ns<fTPCParam->GetNRowUp();ns++){

    r1 = fTPCParam->GetPadRowRadiiUp(ns)-xx;
    r2 = r1 + stripThick;

    dm[5] = r1;
    dm[6] = r2;

    dm[8] = r1;
    dm[9] = r2;

    if(r2+xx < rmax){
      dm[4] = 0.002;
    }
    else{
      dm[4] = 0.202;
    }

    gMC->Gsposp("TLST",ns+1,"TLGA",xx,0.,0.,0,"ONLY",dm,10);

  }
  
  //------------------------------------------------------------------
  // Positioning of the S-sector readout chambers
  //------------------------------------------------------------------

  Float_t zs = 0.5*(250.+0.002);

  Float_t theta1,theta2,theta3;
  Float_t phi1,phi2,phi3;
  Float_t alpha;
  Float_t x,y;

  for(ns=0;ns<nInnerSector;ns++){
    
    phi1 = ns * innerOpenAngle + 270.*kDegrad + innerAngleShift;
    phi1 *= kRaddeg; // in degrees

    phi1 = (Float_t)TMath::Nint(phi1);

    if (phi1 > 360.) phi1 -= 360.;

      
    theta1 = 90.;
    phi2   = 90.;
    theta2 = 180.;
    phi3   = ns * innerOpenAngle + innerAngleShift;
    phi3 *= kRaddeg; // in degrees

    phi3 = (Float_t)TMath::Nint(phi3);
      
    if(phi3 > 360.) phi3 -= 360.;

    theta3 = 90.;

    alpha = phi3*kDegrad;

    x = xCenterS * TMath::Cos(alpha);
    y = xCenterS * TMath::Sin(alpha); 
 
    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);  
     
    gMC->Gspos("TRCS", ns+1, "TPSG", x, y, z1, idrotm[nRotMat], "ONLY");

    gMC->Gspos("TSGA",ns+1,"TGAS",x,y,zs,idrotm[nRotMat], "ONLY");
    gMC->Gspos("TSGA",ns+1+nInnerSector,"TGAS",x,y,-zs,idrotm[nRotMat], "ONLY");

    nRotMat++;     

  }
    
  //-------------------------------------------------------------------
  //  Positioning of the L-sectors readout chambers
  //-------------------------------------------------------------------
    
  for(ns=0;ns<nOuterSector;ns++){
    phi1 = ns * outerOpenAngle + 270.*kDegrad + outerAngleShift;
    phi1 *= kRaddeg; // in degrees

    phi1 = (Float_t)TMath::Nint(phi1);
    

    if (phi1 > 360.) phi1 -= 360.;
      
    theta1 = 90.;
    phi2   = 90.;
    theta2 = 180.;
    phi3   = ns * outerOpenAngle+outerAngleShift;
    phi3 *= kRaddeg; // in degrees

    phi3 = (Float_t)TMath::Nint(phi3);

      
    if(phi3 > 360.) phi3 -= 360.;

    theta3 = 90.;

    alpha = phi3*kDegrad;

    x = xCenterL * TMath::Cos(alpha);
    y = xCenterL * TMath::Sin(alpha); 
 
    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);  
     

    gMC->Gspos("TRCL", ns+1, "TPSG", x, y, z1, idrotm[nRotMat], "ONLY"); 

    nRotMat++;   

  }

  //-------------------------------------------------------------------
  // Positioning of the L-sectors (gas sectors)
  //-------------------------------------------------------------------

  for(ns=0;ns<nOuterSector;ns++){

     phi1 = ns*outerOpenAngle + outerAngleShift;
     phi1 *= kRaddeg;
    
     phi1 = (Float_t)TMath::Nint(phi1);
     if(phi1>360.) phi1 -= 360.;

     theta1 = 90.;

     phi2 = 90. + phi1;
     if(phi2>360.) phi2 -= 360.;

     theta2 = 90.; 

     phi3 = 0.;
     theta3 = 0.;

     alpha = phi1*kDegrad;


     AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3); 


     gMC->Gspos("TLGA",ns+1,"TGAS" ,0.,0.,0.,idrotm[nRotMat],"ONLY");

     nRotMat++;
   
     // reflection !!

     phi3 = 0.;
     theta3 = 180.;
     
     AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);
  
     gMC->Gspos("TLGA",ns+1+nOuterSector,"TGAS" ,0.,0.,0.,idrotm[nRotMat],"ONLY");
          
    nRotMat++;
  }
  
  Float_t z0 = zSide - 0.95;

  gMC->Gspos("TPEC",1,"TPSG",0.,0.,z0,0,"ONLY");

  // ========================================================== 
  //                  wheels 
  // ========================================================== 

  //
  //  auxilary structures
  //


  gMC->Gsvolu("TPWI","TUBE",idtmed[24],dm,0); // "air" 

  // ---------------------------------------------------------- 
  //       Large wheel -> positioned in the TPC 
  // ---------------------------------------------------------- 
  

  z0 = 263.5; // TPC length - 1/2 spoke wheel width

  dm[0] = 258.;
  dm[1] = 278.;
  dm[2] = 11.5;
  
  gMC->Gsvolu("TPWL", "TUBE", idtmed[0], dm, 3); 

  dm[0] = dm[0]+2.;
  dm[1] = 278.;
  dm[2] = dm[2]-2.;

  gMC->Gsposp("TPWI",1,"TPWL",0.,0.,0.,0,"ONLY",dm,3);

  gMC->Gspos("TPWL", 1, "TPC ", 0, 0, z0, 0, "ONLY");
  gMC->Gspos("TPWL", 2, "TPC ", 0, 0, -z0, 0, "ONLY");

  //
  //  Outer vessel + CO2 HV degrader
  //

  dm[0] = 260.;
  dm[1] = 278.;
  dm[2] = 252.;

  gMC->Gsvolu("TPCO","TUBE",idtmed[12],dm,3);

  dm[0] = 275.;
  dm[1] = 278.;
  
  gMC->Gsvolu("TPOV","TUBE",idtmed[10],dm,3);

  gMC->Gspos("TPOV",1,"TPCO",0.,0.,0.,0,"ONLY");


 

  gMC->Gspos("TPCO",1,"TPC ",0.,0.,0.,0,"ONLY");


  //----------------------------------------------------------
  //  Small wheel -> positioned in "side gas
  //----------------------------------------------------------

  dm[0] = 78.;
  dm[1] = 82.;
  dm[2] = 11.5;

  gMC->Gsvolu("TPWS", "TUBE", idtmed[0], dm, 3);

  dm[0] = 78.;
  dm[1] = dm[1]-2;
  dm[2] = dm[2]-2.;

  gMC->Gsvolu("TPW1", "TUBE", idtmed[2], dm, 3);
  
  gMC->Gspos("TPW1", 1, "TPWS", 0., 0., 0., 0, "ONLY");

  z0 = 1.; // spoke wheel is shifted w.r.t. center of the "side gas"

  gMC->Gspos("TPWS", 1, "TPSG", 0, 0, z0, 0, "ONLY");


  // to avoid overlaps

  dm[0] = 76.;
  dm[1] = 78.;
  dm[2] = 11.5;

  gMC->Gsvolu("TPS1","TUBE",idtmed[0],dm,3);

  dm[2] = 9.5;

  gMC->Gsvolu("TPS2","TUBE",idtmed[24],dm,3);

  gMC->Gspos("TPS2",1,"TPS1",0.,0.,0.,0,"ONLY");

  z0= 263.5;
  
  gMC->Gspos("TPS1",1,"TPC ",0.,0.,z0,0,"ONLY");
  gMC->Gspos("TPS1",2,"TPC ",0.,0.,-z0,0,"ONLY");

  // G10 plug

  dm[0] = 76.;
  dm[1] = 78.;
  dm[2] = 1.;

  gMC->Gsvolu("TPG2","TUBE",idtmed[8],dm,3);

  z0 = 251.;

  gMC->Gspos("TPG2",1,"TPC ",0.,0.,z0,0,"ONLY");
  gMC->Gspos("TPG2",2,"TPC ",0.,0.,-z0,0,"ONLY");


  //---------------------------------------------------------
  //  central wheel  6 (radial direction) x 4 (along z) cm2
  //---------------------------------------------------------

  dm[0] = 140.;
  dm[1] = 146.;
  dm[2] = 2.;

  gMC->Gsvolu("TPWC","TUBE",idtmed[0],dm,3);

  dm[0] = dm[0] + 2.;
  dm[1] = dm[1] - 2.;
  dm[2] = dm[2] - 1.;

  gMC->Gsposp("TPWI",2,"TPWC",0.,0.,0.,0,"ONLY",dm,3);

  z0 = zSide - 1.9 - 2.;

  gMC->Gspos("TPWC",1,"TPSG",0.,0.,z0,0,"ONLY");

  //

  gMC->Gsvolu("TPSE","BOX ",idtmed[24],dm,0); // "empty" part of the spoke 

 
  //---------------------------------------------------------
  //  inner spokes (nSectorInner)
  //---------------------------------------------------------

  dm[0] = 0.5*(139.9-82.1);
  dm[1] = 3.;
  dm[2] = 2.;

  Float_t x1 = dm[0]+82.;

  gMC->Gsvolu("TPSI","BOX",idtmed[0],dm,3);

  dm[1] = dm[1]-1.;
  dm[2] = dm[2]-1.;

  gMC->Gsposp("TPSE",1,"TPSI",0.,0.,0.,0,"ONLY",dm,3);

  for(ns=0;ns<nInnerSector;ns++){

    phi1 = 0.5*innerOpenAngle + ns*innerOpenAngle + innerAngleShift;
    theta1=90.;
    phi1 *=kRaddeg;

    phi1 = (Float_t)TMath::Nint(phi1);
    if(phi1>360.) phi1 -= 360.;    

    phi2 = phi1+90.;
    if(phi2>360.) phi2 -= 360.;
    theta2=90.;
    phi3=0.;
    theta3=0.;

    alpha = phi1 * kDegrad;
    x     = x1 * TMath::Cos(alpha);
    y     = x1 * TMath::Sin(alpha);    

   AliMatrix(idrotm[nRotMat],theta1,phi1,theta2,phi2,theta3,phi3);

   gMC->Gspos("TPSI",ns+1,"TPSG",x,y,z0,idrotm[nRotMat],"ONLY");  

   nRotMat++;

  }

  //-------------------------------------------------------------
  // outer spokes (nSectorOuter)
  //-------------------------------------------------------------

  dm[0] = 0.5*(257.9-146.1);
  dm[1] = 3.;
  dm[2] = 2.;

  x1 = dm[0] + 146.;

  gMC->Gsvolu("TPSO","BOX ",idtmed[0],dm,3);

  dm[1] = dm[1] - 1.;
  dm[2] = dm[2] - 1.;

  gMC->Gsposp("TPSE",2,"TPSO",0.,0.,0.,0,"ONLY",dm,3);

  for(ns=0;ns<nOuterSector;ns++){

    phi1 = 0.5*outerOpenAngle + ns*outerOpenAngle + outerAngleShift;
    theta1=90.;
    phi1 *=kRaddeg;

    phi1 = (Float_t)TMath::Nint(phi1);
    if(phi1>360.) phi1 -= 360.;

    phi2 = phi1+90.;
    if(phi2>360.) phi2 -= 360.;
    theta2=90.;
    phi3=0.;
    theta3=0.;

    alpha = phi1 * kDegrad;
    x     = x1 * TMath::Cos(alpha);
    y     = x1 * TMath::Sin(alpha);    

   AliMatrix(idrotm[nRotMat],theta1,phi1,theta2,phi2,theta3,phi3);

   gMC->Gspos("TPSO",ns+1,"TPSG",x,y,z0,idrotm[nRotMat],"ONLY");  

   nRotMat++;

  }  
  

  
  // -------------------------------------------------------- 
  //         put the readout chambers into the TPC 
  // -------------------------------------------------------- 

  theta1 = 90.;
  phi1   = 0.;
  theta2 = 90.;
  phi2   = 270.;
  theta3 = 180.;
  phi3   = 0.;
  
  AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);
  
  z0 = zSide + 250.;
  
  gMC->Gspos("TPSG", 1, "TPC ", 0, 0, z0, 0, "ONLY");
  gMC->Gspos("TPSG", 2, "TPC ", 0, 0, -z0, idrotm[nRotMat], "ONLY");
  
  gMC->Gspos("TPC ", 1, "ALIC", 0, 0, 0, 0, "ONLY");

  //----------------------------------------------------
  //  Inner vessel and HV degrader
  //----------------------------------------------------

  dm[0] = 0.;
  dm[1] = 360.;
  dm[2] = 4.;
  
  dm[3] = -250.;
  dm[4] = 74.4;
  dm[5] = 76.;

  dm[6] = -64.5;
  dm[7] = 50.;
  dm[8] = 76.;

  dm[9] = 64.5;
  dm[10] = 50.;
  dm[11] = 76.;

  dm[12] = 250.;
  dm[13] = 74.4;
  dm[14] = 76.;

  gMC->Gsvolu("TPVD", "PCON", idtmed[12], dm, 15); // CO2

  // cone parts

  dm[0] = 0.;
  dm[1] = 360.;
  dm[2] = 2.;

  dm[3] = 64.5;
  dm[4] = 50.;
  dm[5] = 51.6;
 
  dm[6] = 250.;
  dm[7] = 74.4;
  dm[8] = 76.;


  gMC->Gsvolu("TIVC","PCON",idtmed[11],dm,9); // C-fibre

  gMC->Gspos("TIVC",1,"TPVD",0.,0.,0.,0,"ONLY");
  gMC->Gspos("TIVC",2,"TPVD",0.,0.,0.,idrotm[nRotMat],"ONLY");

  // barrel part

  dm[0] = 50.;
  dm[1] = 50.5;
  dm[2] = 64.5;

  gMC->Gsvolu("TIVB","TUBE",idtmed[9],dm,3);

  gMC->Gspos("TIVB",1,"TPVD",0.,0.,0.,0,"ONLY");

  gMC->Gspos("TPVD",1,"ALIC",0.,0.,0.,0,"ONLY");

  

  

  // --------------------------------------------------- 
  //               volumes ordering 
  // --------------------------------------------------- 
  gMC->Gsord("TGAS", 6);
  gMC->Gsord("TPSG", 6);
  gMC->Gsord("TSGA", 3);
  gMC->Gsord("TLGA", 4);

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
  gMC->Gsatt("TPC","SEEN",0);
  gMC->Gsatt("TGAS","SEEN",0);
  gMC->Gsatt("TPSG","SEEN",0);
  gMC->Gsatt("TPHV","SEEN",1);
  gMC->Gsatt("TPMH","SEEN",1);
  gMC->Gsatt("TPEC","SEEN",0);
  gMC->Gsatt("TRCS","SEEN",1);
  gMC->Gsatt("TRCL","SEEN",1);
  gMC->Gsatt("TPWL","SEEN",1);
  gMC->Gsatt("TPWI","SEEN",1);
  gMC->Gsatt("TPWS","SEEN",1);
  gMC->Gsatt("TPW1","SEEN",1);
  gMC->Gsatt("TPS1","SEEN",1);
  gMC->Gsatt("TPS2","SEEN",1);
  gMC->Gsatt("TPG1","SEEN",1);
  gMC->Gsatt("TPG2","SEEN",1);
  gMC->Gsatt("TPWC","SEEN",1);
  gMC->Gsatt("TPSI","SEEN",1); 
  gMC->Gsatt("TPSO","SEEN",1);
  gMC->Gsatt("TPCO","SEEN",1);
  gMC->Gsatt("TPOV","SEEN",1);
  gMC->Gsatt("TPVD","SEEN",1);
  
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
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

  fIdSens1=gMC->VolId("TLST"); // L-sector
  fIdSens2=gMC->VolId("TSST"); // S-sector

  printf("TPC version 1 initialized\n");
}

//_____________________________________________________________________________
void AliTPCv1::StepManager()
{
  //
  // Called at every step in the Time Projection Chamber
  //
  Int_t         copy, id, i;
  Float_t       hits[4];
  Int_t         vol[2];
  TLorentzVector p;
  TClonesArray &lhits = *fHits;
  
  //

  if(gMC->TrackCharge() && gMC->IsTrackEntering()) {
    //
    // Only entering charged tracks
    //
    if((id=gMC->CurrentVolID(copy))==fIdSens1) {

      // L

      vol[1]=copy-1; // row
      id=gMC->CurrentVolOffID(1,copy);
      vol[0]=copy+fTPCParam->GetNInnerSector()-1; // sector

      //I.Belikov's modification
      if (vol[1]==0) {
       gMC->TrackMomentum(p);
       hits[0]=p[0];
       hits[1]=p[1];
       hits[2]=p[2];
       hits[3]=0.; // this hit has no energy loss
       new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);

       gMC->TrackPosition(p);
       hits[0]=p[0];
       hits[1]=p[1];
       hits[2]=p[2];
       hits[3]=0.; // this hit has no energy loss
       new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
      }
      //end of I.Belikov's modification
    } else if(id==fIdSens2) {

      // S

      vol[1]=copy-1; // row
      id=gMC->CurrentVolOffID(1,copy); // sector
      vol[0]=copy-1;

      //I.Belikov's modification
      if (vol[1]==0) {
       gMC->TrackMomentum(p);
       hits[0]=p[0];
       hits[1]=p[1];
       hits[2]=p[2];
       hits[3]=0.; // this hit has no energy loss
       new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);

       gMC->TrackPosition(p);
       hits[0]=p[0];
       hits[1]=p[1];
       hits[2]=p[2];
       hits[3]=0.; // this hit has no energy loss
       new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
      }
      //end of I.Belikov's modification

    } else return;

    gMC->TrackPosition(p);
    for(i=0;i<3;++i) hits[i]=p[i];
    hits[3]=1; //I'd like to have something positive here (I.Belikov)
    new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }
}



