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
Revision 1.14  2000/04/17 09:37:33  kowal2
removed obsolete AliTPCDigitsDisplay.C

Revision 1.13.8.2  2000/04/10 08:33:44  kowal2

Updated readout chambers

Revision 1.13.8.1  2000/04/10 07:56:53  kowal2
Not used anymore - removed

Revision 1.13  1999/11/04 17:28:06  fca
Correct barrel part of HV Degrader

Revision 1.12  1999/10/08 06:27:23  fca
Corrected bug in the HV degrader geometry, thanks to G.Tabary

Revision 1.11  1999/10/04 13:39:54  fca
Correct array index problem

Revision 1.10  1999/09/29 09:24:34  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber version 0 -- "coarse" TPC                        //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTPCv0Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include "AliTPCv0.h"
#include "AliRun.h"
#include <iostream.h>
#include <fstream.h>
#include "AliMC.h"
#include "AliConst.h"

#include "AliTPCParam.h"
#include "AliTPCDigitsArray.h"

ClassImp(AliTPCv0)
 
//_____________________________________________________________________________
AliTPCv0::AliTPCv0(const char *name, const char *title) 
         :AliTPC(name, title)
{
  //
  // Standard creator for TPC version 0
  //
}

//_____________________________________________________________________________
void AliTPCv0::CreateGeometry()
{
  //
  // Creation of the TPC coarse geometry (version 0)
  // Origin Marek Kowalski Cracow
  //
  //Begin_Html
  /*
    <img src="picts/AliTPCv0.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTPCv0Tree.gif">
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

  Float_t z_side = dm[2]; // 1/2 of the side gas thickness

  //-----------------------------------------------------------
  //   Readout chambers , 25% of X0, I use Al as the material
  //-----------------------------------------------------------

  Float_t InnerOpenAngle = fTPCParam->GetInnerAngle();
  Float_t OuterOpenAngle = fTPCParam->GetOuterAngle();

  Float_t InnerAngleShift = fTPCParam->GetInnerAngleShift();
  Float_t OuterAngleShift = fTPCParam->GetOuterAngleShift();


  Int_t nInnerSector = fTPCParam->GetNInnerSector()/2;
  Int_t nOuterSector = fTPCParam->GetNOuterSector()/2;


  Float_t InSecLowEdge = fTPCParam->GetInnerRadiusLow();
  Float_t InSecUpEdge =  fTPCParam->GetInnerRadiusUp();

  Float_t OuSecLowEdge = fTPCParam->GetOuterRadiusLow();
  Float_t OuSecUpEdge = fTPCParam->GetOuterRadiusUp();

  Float_t SecThick = 2.225; // Al


  //  S (Inner) sectors

  Float_t LowEdge = fTPCParam->GetInnerFrameSpace();

  dm[0] = InSecLowEdge*TMath::Tan(0.5*InnerOpenAngle)-LowEdge;
  dm[1] = InSecUpEdge*TMath::Tan(0.5*InnerOpenAngle)-LowEdge;
  dm[2] = SecThick;
  dm[3] = 0.5*(InSecUpEdge-InSecLowEdge);

  Float_t xCenterS = InSecLowEdge+dm[3];

  gMC->Gsvolu("TRCS", "TRD1", idtmed[0], dm, 4); 

  //  L (Outer) sectors

  Float_t UpEdge = fTPCParam->GetOuterFrameSpace();

  dm[0] = OuSecLowEdge*TMath::Tan(0.5*OuterOpenAngle)-UpEdge;
  dm[1] = OuSecUpEdge*TMath::Tan(0.5*OuterOpenAngle)-UpEdge;
  dm[2] = SecThick;
  dm[3] = 0.5*(OuSecUpEdge-OuSecLowEdge);

  Float_t xCenterL = OuSecLowEdge+dm[3];  

  gMC->Gsvolu("TRCL", "TRD1", idtmed[0], dm, 4);

  Float_t z1 = -z_side + SecThick*0.5;

  //------------------------------------------------------------------
  // Positioning of the S-sector readout chambers
  //------------------------------------------------------------------

  Int_t ns;
  Float_t theta1,theta2,theta3;
  Float_t phi1,phi2,phi3;
  Float_t alpha;
  Float_t x,y;

  for(ns=0;ns<nInnerSector;ns++){
    
    phi1 = ns * InnerOpenAngle + 270.*kDegrad + InnerAngleShift;
    phi1 *= kRaddeg; // in degrees

    phi1 = (Float_t)TMath::Nint(phi1);

    if (phi1 > 360.) phi1 -= 360.;
      
    theta1 = 90.;
    phi2   = 90.;
    theta2 = 180.;
    phi3   = ns * InnerOpenAngle + InnerAngleShift;
    phi3 *= kRaddeg; // in degrees

    phi3 = (Float_t)TMath::Nint(phi3);
      
    if(phi3 > 360.) phi3 -= 360.;

    theta3 = 90.;

    alpha = phi3*kDegrad;

    x = xCenterS * TMath::Cos(alpha);
    y = xCenterS * TMath::Sin(alpha); 
 
    AliMatrix(idrotm[nRotMat], theta1, phi1, theta2, phi2, theta3, phi3);  
     
    gMC->Gspos("TRCS", ns+1, "TPSG", x, y, z1, idrotm[nRotMat], "ONLY");

    nRotMat++;     

  }
    
  //-------------------------------------------------------------------
  //  Positioning of the L-sectors readout chambers
  //-------------------------------------------------------------------
    
  for(ns=0;ns<nOuterSector;ns++){
    phi1 = ns * OuterOpenAngle + 270.*kDegrad + OuterAngleShift;
    phi1 *= kRaddeg; // in degrees

    phi1 = (Float_t)TMath::Nint(phi1);
    

    if (phi1 > 360.) phi1 -= 360.;
      
    theta1 = 90.;
    phi2   = 90.;
    theta2 = 180.;
    phi3   = ns * OuterOpenAngle+OuterAngleShift;
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

  Float_t z0 = z_side - 0.95;

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

  z0 = z_side - 1.9 - 2.;

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

    phi1 = 0.5*InnerOpenAngle + ns*InnerOpenAngle + InnerAngleShift;
    theta1=90.;
    phi1 *=kRaddeg;

    phi1 = (Float_t)TMath::Nint(phi1);

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

    phi1 = 0.5*OuterOpenAngle + ns*OuterOpenAngle + OuterAngleShift;
    theta1=90.;
    phi1 *=kRaddeg;

    phi1 = (Float_t)TMath::Nint(phi1);

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
  
  z0 = z_side + 250.;
  
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
  gMC->Gsord("TPSG", 6);
 
} // end of function


//_____________________________________________________________________________
void AliTPCv0::CreateMaterials()
{
  //
  // Define materials for the TPC
  //
  AliTPC::CreateMaterials();
}

//_____________________________________________________________________________
void AliTPCv0::DrawDetector()
{
  //
  // Draw a shaded view of the Time Projection Chamber version 0
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
  //
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
void AliTPCv0::Init()
{
  //
  // Initialise Time Projection Chamber version 0
  //
  printf("TPC version 0 initialized\n");
}

//_____________________________________________________________________________
void AliTPCv0::StepManager()
{
  //
  // Procedure called at each step in the TPC
  //
}
