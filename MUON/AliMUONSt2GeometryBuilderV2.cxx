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

// $Id$
//
// Class AliMUONSt2GeometryBuilderV2
// -------------------------------
// MUON Station2 coarse geometry construction class.
//********************************************************************
// Author: SANJOY PAL ,Prof. SUKALYAN CHATTOPADHAYAY  [SINP, KOLKATA]
//         &  Dr.SHAKEEL AHMAD (AMU), INDIA
//********************************************************************


#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <Riostream.h>

#include "AliRun.h"
#include "AliLog.h"

#include "AliMUONSt2GeometryBuilderV2.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONConstants.h"

#define PI 3.14159
ClassImp(AliMUONSt2GeometryBuilderV2)

//______________________________________________________________________________
AliMUONSt2GeometryBuilderV2::AliMUONSt2GeometryBuilderV2(AliMUON* muon)
 : AliMUONVGeometryBuilder("st2V2.dat",
                           muon->Chamber(2).GetGeometry(),
			   muon->Chamber(3).GetGeometry()),
   fMUON(muon)
{
// Standard constructor

}

//______________________________________________________________________________
AliMUONSt2GeometryBuilderV2::AliMUONSt2GeometryBuilderV2()
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONSt2GeometryBuilderV2::AliMUONSt2GeometryBuilderV2(const AliMUONSt2GeometryBuilderV2& rhs)
  : AliMUONVGeometryBuilder(rhs)
{
// Protected copy constructor

  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONSt2GeometryBuilderV2::~AliMUONSt2GeometryBuilderV2() {
//
}

//______________________________________________________________________________
AliMUONSt2GeometryBuilderV2&
AliMUONSt2GeometryBuilderV2::operator = (const AliMUONSt2GeometryBuilderV2& rhs)
{
// Protected assignement operator

  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");

  return *this;
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSt2GeometryBuilderV2::CreateGeometry()
{
// From AliMUONv1::CreateGeometry()

//
//********************************************************************
//                            Station 2                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers

//     AliMUONChamber* iChamber = &fMUON->Chamber(2);
//     AliMUONChamber* iChamber1 = iChamber;
//     AliMUONChamber* iChamber2 = &fMUON->Chamber(3);


     // Get tracking medias Ids
     Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;
     Int_t idAir  = idtmed[1100]; // medium 1
     Int_t idGas  = idtmed[1108]; // medium Ar-CO2 gas (80%+20%)
     Int_t idPCB  = idtmed[1122]; // medium FR4
     Int_t idCU   = idtmed[1110]; // medium copper
     Int_t idRoha = idtmed[1113]; // medium roha cell
     Int_t idPGF30= idtmed[1123]; // medium for Frame Eq.to Bakelite
     Int_t idScru = idtmed[1128]; // screw material - Stainless Steel(18%Cr,9%Ni,Fe)


     Int_t irot1, irot2, irot3, irot4, irot5, irot6;

     fMUON->AliMatrix(irot1,  90.,  13., 90., 103.,  0., 0.);  //+13deg in x-y Plane
     fMUON->AliMatrix(irot2,  90.,  34., 90., 124.,  0., 0.); // +34deg in x-y plane
     fMUON->AliMatrix(irot3,  90.,  56., 90., 146.,  0., 0.); // +56 deg in x-y plane
     fMUON->AliMatrix(irot4,  90.,  76., 90., 166.,  0., 0.); // +76 deg in x-y plane
     fMUON->AliMatrix(irot5,  90.,  90., 90., 180.,  0., 0.); // +90 deg in x-y plane
     fMUON->AliMatrix(irot6,  90., 22.5, 90., 112.5, 0., 0.); //22.5 deg in x-y Plane

/*########################################################################################
    Create volume for one Quadrant connsist two plane
##########################################################################################*/
     Float_t tpar1[5];
     tpar1[0] = 20.6;
     tpar1[1] = 122.0;
     tpar1[2] = 6.6/2;
     tpar1[3] = -10.0;
     tpar1[4] = 100.0;


     gMC->Gsvolu("SQM3","TUBS", idAir, tpar1, 5);
     gMC->Gsvolu("SQM4","TUBS", idAir, tpar1, 5);



//==================================================================================
//                                 Plane - 1       L is used for one Plane while
//                                                   R for the other
//==================================================================================

//Thickness of variour parts
       Float_t zcbb  = 0.04;       //cathode pcb
       Float_t zcu   = 0.004;      // eff. cu in cathode pcb
       Float_t zRoha = 2.5;        // Rhocell
       Float_t zmeb = 0.08;      //Mech. exit board
       Float_t zcu2 = 0.02;     //Effective electronic exit board
    //Z-positions of various parts--- in Plane-1

       Float_t zpos_cbb = 0.25;   // 2.5 mm => gap between anode & chatode plane
       Float_t zpos_cu = zcbb;
       Float_t zpos_Roha  = zcu;
       Float_t zpos_meb = zmeb;
       Float_t zpos_cu2 = zcu2;

       Float_t zpos_cbb_bar = zpos_cbb + zcbb/2.;  //for segment 0 & 6
       Float_t zpos_cubar = zcbb/2. + zcu/2.;
       Float_t zpos_Roha_bar = zcu/2. + zRoha/2.;
       Float_t zpos_meb_bar =  zRoha/2. + zmeb/2.;
       Float_t zpos_eeb_bar = zmeb/2. + zcu2/2.;


 //Cathode PCB + Copper sheet + Rohacell + mech exit board + eff. electronic exit board

 //Segment-0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Float_t bpar_h[3];
       bpar_h[0] = 95.5/2.;
       bpar_h[1] = 1.6/2.;
       bpar_h[2] = zcbb/2.;
       gMC->Gsvolu("CB0L", "BOX", idPCB, bpar_h, 3);

       bpar_h[2] = zcu/2.;     //Thickness of Copper sheet
       gMC->Gsvolu("CU0L", "BOX", idCU, bpar_h, 3);

       bpar_h[2] = zRoha/2.;     //Thickness of Roha cell
       gMC->Gsvolu("RH0L", "BOX", idRoha, bpar_h, 3);

       bpar_h[0] = (100.6/2)-(0.9/2);
       bpar_h[1] = 2.8/2.;
       bpar_h[2] = zmeb/2;
       gMC->Gsvolu("MB0L", "BOX", idPCB, bpar_h, 3);

       bpar_h[2] = zcu2/2;
       gMC->Gsvolu("EB0L", "BOX", idCU, bpar_h, 3);

 //Segment-1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Float_t pgpar[10];
        pgpar[0] = 0.;
        pgpar[1] = 13.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = 0.;
        pgpar[5] = 22.5;
        pgpar[6] = 120.2;
        pgpar[7] = pgpar[4] + zcbb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB1L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU1L", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] + zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH1L", "PGON", idRoha, pgpar, 10);

        pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 121.2;
        pgpar[7] = pgpar[4] + zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB1L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu2;
        gMC->Gsvolu("EB1L", "PGON", idCU, pgpar, 10);

//Segment-2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 21.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = 0.;
        pgpar[5] = 22.5;
        pgpar[6] = 121.5;
        pgpar[7] = pgpar[4] +zcbb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB2L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU2L", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] + zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH2L", "PGON", idRoha, pgpar, 10);

        pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 122.5;
        pgpar[7] = pgpar[4] + zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB2L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu2;
        gMC->Gsvolu("EB2L", "PGON", idCU, pgpar, 10);

//Segment-3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 22.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = 0.;
        pgpar[5] = 22.5;
        pgpar[6] = 119.9;
        pgpar[7] = pgpar[4] + zcbb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB3L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU3L", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] + zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH3L", "PGON", idRoha, pgpar, 10);

        pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 120.9;
        pgpar[7] = pgpar[4] + zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB3L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu2;
        gMC->Gsvolu("EB3L", "PGON", idCU, pgpar, 10);

//Segment-4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 20.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = 0.;
        pgpar[5] = 22.5;
        pgpar[6] = 119.9;
        pgpar[7] = pgpar[4] + zcbb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB4L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU4L", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] + zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH4L", "PGON", idRoha, pgpar, 10);

        pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 120.9;
        pgpar[7] = pgpar[4] + zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB4L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu2;
        gMC->Gsvolu("EB4L", "PGON", idCU, pgpar, 10);

//Segment-5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 14.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = 0.;
        pgpar[5] = 22.5;
        pgpar[6] = 117.5;
        pgpar[7] = pgpar[4] + zcbb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB5L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU5L", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] + zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH5L", "PGON", idRoha, pgpar, 10);

        pgpar[5] = 21.5;// radius of inner rib
        pgpar[6] = 118.5;
        pgpar[7] = pgpar[4] + zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB5L", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] + zcu2;
        gMC->Gsvolu("EB5L", "PGON", idCU, pgpar, 10);

//Segment-6 - vertical box ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	Float_t bpar_v[3];

        bpar_v[0] = 1.6/2.;
	bpar_v[1] = 95.5/2.;
	bpar_v[2] = zcbb/2.;
        gMC->Gsvolu("CB6L", "BOX", idPCB, bpar_v, 3);

	bpar_v[2] = zcu/2.;
        gMC->Gsvolu("CU6L", "BOX", idCU, bpar_v, 3);

	bpar_v[2] = zRoha/2.;
        gMC->Gsvolu("RH6L", "BOX", idRoha, bpar_v, 3);

        bpar_v[0] = 2.8/2.;
	bpar_v[1] = (97.9/2)-(0.9/2.);
        bpar_v[2] = zmeb/2;
        gMC->Gsvolu("MB6L", "BOX", idPCB, bpar_v, 3);

        bpar_v[2] = zcu2/2;
        gMC->Gsvolu("EB6L", "BOX", idCU, bpar_v, 3);

//...........................................................................................
//Positioning of Electronic exit board
     gMC->Gspos("EB0L",1, "MB0L", 0.,0.,zpos_eeb_bar,0, "only");
     gMC->Gspos("EB1L",1, "MB1L", 0.,0.,zpos_cu2,0, "only");
     gMC->Gspos("EB2L",1, "MB2L", 0.,0.,zpos_cu2,0, "only");
     gMC->Gspos("EB3L",1, "MB3L", 0.,0.,zpos_cu2,0, "only");
     gMC->Gspos("EB4L",1, "MB4L", 0.,0.,zpos_cu2,0, "only");
     gMC->Gspos("EB5L",1, "MB5L", 0.,0.,zpos_cu2,0, "only");
     gMC->Gspos("EB6L",1, "MB6L", 0.,0.,zpos_eeb_bar,0, "only");

//Positioning of Mech. exit board
     Float_t xpos_meb_hor =  1.10;
     Float_t ypos_meb_hor = -0.60;

     Float_t xpos_meb_ver = -0.60;
     Float_t ypos_meb_ver = -0.35;

     gMC->Gspos("MB0L",1, "RH0L", xpos_meb_hor,ypos_meb_hor,zpos_meb_bar,0, "only");
     gMC->Gspos("MB1L",1, "RH1L", 0.,0.,zpos_meb,0, "only");
     gMC->Gspos("MB2L",1, "RH2L", 0.,0.,zpos_meb,0, "only");
     gMC->Gspos("MB3L",1, "RH3L", 0.,0.,zpos_meb,0, "only");
     gMC->Gspos("MB4L",1, "RH4L", 0.,0.,zpos_meb,0, "only");
     gMC->Gspos("MB5L",1, "RH5L", 0.,0.,zpos_meb,0, "only");
     gMC->Gspos("MB6L",1, "RH6L", xpos_meb_ver,ypos_meb_ver,zpos_meb_bar,0, "only");

//Positioning Roha cell over copper sheet
     gMC->Gspos("RH0L",1, "CU0L", 0.,0.,zpos_Roha_bar,0, "only");// box horizontal
     gMC->Gspos("RH1L",1, "CU1L", 0.,0.,zpos_Roha,0, "only");
     gMC->Gspos("RH2L",1, "CU2L", 0.,0.,zpos_Roha,0, "only");
     gMC->Gspos("RH3L",1, "CU3L", 0.,0.,zpos_Roha,0, "only");
     gMC->Gspos("RH4L",1, "CU4L", 0.,0.,zpos_Roha,0, "only");
     gMC->Gspos("RH5L",1, "CU5L", 0.,0.,zpos_Roha,0, "only");
     gMC->Gspos("RH6L",1, "CU6L", 0.,0.,zpos_Roha_bar,0, "only");

//Positioning Copper sheet over PCB
     gMC->Gspos("CU0L",1, "CB0L", 0.,0.,zpos_cubar,0, "only");// box horizontal
     gMC->Gspos("CU1L",1, "CB1L", 0.,0.,zpos_cu,0, "only");
     gMC->Gspos("CU2L",1, "CB2L", 0.,0.,zpos_cu,0, "only");
     gMC->Gspos("CU3L",1, "CB3L", 0.,0.,zpos_cu,0, "only");
     gMC->Gspos("CU4L",1, "CB4L", 0.,0.,zpos_cu,0, "only");
     gMC->Gspos("CU5L",1, "CB5L", 0.,0.,zpos_cu,0, "only");
     gMC->Gspos("CU6L",1, "CB6L", 0.,0.,zpos_cubar,0, "only");


//Positioning the PCB
      Float_t x_hor_pos = 95.5/2.+ 22.6; // 22.6 is inner radius of PCB
      Float_t y_hor_pos = -(1.6/2); //

      Float_t x_ver_pos = -(1.6/2); //
      Float_t y_ver_pos = 95.5/2.+ 22.6; //

      // chamber 3
     gMC->Gspos("CB0L",1, "SQM3", x_hor_pos,y_hor_pos,zpos_cbb_bar,0, "only");// box horizontal
     gMC->Gspos("CB1L",1, "SQM3", 0.,0.,zpos_cbb,0, "only");
     gMC->Gspos("CB2L",1, "SQM3", 0.,0.,zpos_cbb,irot1,"only");
     gMC->Gspos("CB3L",1, "SQM3", 0.,0.,zpos_cbb,irot2, "only");
     gMC->Gspos("CB4L",1, "SQM3", 0.,0.,zpos_cbb,irot3, "only");
     gMC->Gspos("CB5L",1, "SQM3", 0.,0.,zpos_cbb,irot4, "only");
     gMC->Gspos("CB6L",1, "SQM3", x_ver_pos,y_ver_pos,zpos_cbb_bar,0, "only");// box vertical

      // chamber 4
     gMC->Gspos("CB0L",1, "SQM4", x_hor_pos,y_hor_pos,zpos_cbb_bar,0, "only");// box horizontal
     gMC->Gspos("CB1L",1, "SQM4", 0.,0.,zpos_cbb,0, "only");
     gMC->Gspos("CB2L",1, "SQM4", 0.,0.,zpos_cbb,irot1,"only");
     gMC->Gspos("CB3L",1, "SQM4", 0.,0.,zpos_cbb,irot2, "only");
     gMC->Gspos("CB4L",1, "SQM4", 0.,0.,zpos_cbb,irot3, "only");
     gMC->Gspos("CB5L",1, "SQM4", 0.,0.,zpos_cbb,irot4, "only");
     gMC->Gspos("CB6L",1, "SQM4", x_ver_pos,y_ver_pos,zpos_cbb,0, "only");// box vertical

//----------------------------------------------------------------------
//                          Frames
//----------------------------------------------------------------------
//Frame-1
      Float_t frame1[3] ;                //Index-1 is used for horizontal frame bar
      frame1[0] = 100.6/2.;              //and 2,3,4..... is used for the next frames
      frame1[1] = 2.5/2.;                //from hor. to the vertical
      frame1[2] = .95/2.;

      Float_t rib1[3];
      rib1[0] = frame1[0] - 0.9/2.;
      rib1[1] = 1./2.;
      rib1[2] = 1.84/2.;

      Float_t xpos_fr1 = frame1[0] + 20.6;
      Float_t ypos_fr1 = -3.7 + frame1[1] ;
      Float_t zpos_fr1 = frame1[2];

      Float_t xpos_rb1 = 0.9/2;
      Float_t ypos_rb1 = -frame1[1] + rib1[1]+0.9;
      Float_t zpos_rb1 = frame1[2] + rib1[2];

      gMC->Gsvolu("FRHL", "BOX", idPGF30, frame1, 3); //Frame - 1
      gMC->Gsvolu("RBHL", "BOX", idPGF30, rib1, 3); //Rib - 1

//Fixing Screws ...........................................
//---------screw parameters (on frame) M4 ------------------
   //in this plane screw-length is taken eq. to thickness of frame/frame+rib
   // the remaining of the actual size is adjusted in the other plane

       Float_t spar[3];
       spar[0] = 0.0;
       spar[1] = 0.2;
       spar[2] = 0.95/2.;

       Float_t s_head[3];
       s_head[0] = 0.0; // screw-head
       s_head[1] = 0.4;
       s_head[2] = 0.4/2.;

       Float_t xpos_h = 0.0;
       Float_t ypos_h = 0.0;
       Float_t zpos_h = spar[2] + s_head[2];

       gMC->Gsvolu("SCHL","TUBE",idScru, spar,3);    //screw-vertical part for Frame M4
       gMC->Gsvolu("HDFL","TUBE",idScru, s_head,3);    //screw-head
   //positioning head over screws
       gMC->Gspos("HDFL",1,"SCHL",xpos_h,ypos_h,zpos_h,0,"ONLY");//positioning Screw-head

   //Screws on Frame
       Float_t xpos_s = frame1[0] - 0.4;
       Float_t ypos_s = -frame1[1] + 0.4;
       Float_t zpos_s =  0.0;

       Int_t Tot_Scru1 = 21;   //Total no. of Scru
       for(Int_t nos1 = 0; nos1 < Tot_Scru1; nos1++)
          {
           gMC->Gspos("SCHL",nos1+1,"FRHL",xpos_s,ypos_s,zpos_s,0,"ONLY");
           xpos_s -= 5.;
	  } // nos1 is no. of scru

       gMC->Gspos("RBHL",1, "FRHL", xpos_rb1, ypos_rb1, zpos_rb1,0, "only");// Rib-1
//       gMC->Gspos("FRHL",1, "SQM3", xpos_fr1, ypos_fr1, zpos_fr1,0, "only");// frame -1
//       gMC->Gspos("FRHL",1, "SQM4", xpos_fr1, ypos_fr1, zpos_fr1,0, "only");// frame -1

//......................................................................................
//Frame-2

       Float_t frame2[3] ;
       frame2[0] = 4./2.;
       frame2[1] = 28.9/2.;
       frame2[2] = frame1[2];

       Float_t xpos_fr2 = frame1[0]  -frame2[0];
       Float_t ypos_fr2 = frame2[1]+frame1[1] ;
       Float_t zpos_fr2 = 0.;

       Float_t rib2[3];
       rib2[0] = 1./2.;
       rib2[1] = frame2[1]+0.6;
       rib2[2] = 1.84/2.;

       Float_t xpos_rb2 = frame2[0] - rib2[0];
       Float_t ypos_rb2 = - 0.6/2;
       Float_t zpos_rb2 = frame2[2] + rib2[2];

       gMC->Gsvolu("FR1L", "BOX", idPGF30, frame2, 3); //Frame - 2
       gMC->Gsvolu("RB1L", "BOX", idPGF30, rib2, 3); //Rib - 2

//Fixing Screws ...........................................
   //screw parameters (on Ribs) M4 x 65
   //Head size & positions remianing the same-----------
       Float_t spar2[3];
       spar2[0] = 0.0;
       spar2[1] = 0.2;
       spar2[2] = 2.9/2.;
       gMC->Gsvolu("SCRL","TUBE",idScru, spar2,3);    //screw-vertical part for Rib m4 x 65
  //positioning head over screws
       Float_t xpos2_h = 0.0;
       Float_t ypos2_h = 0.0;
       Float_t zpos2_h = spar2[2] + s_head[2];
       gMC->Gspos("HDFL",2,"SCRL",xpos2_h,ypos2_h,zpos2_h,0,"ONLY");//positioning Screw-head

       Float_t xpos2_s = -rib2[0] + 0.4;
       Float_t ypos2_s = -rib2[1] + 1.7;
       Float_t zpos2_s = -frame2[2];

       Int_t Tot_Scru2 = 6;   // 6 screws
       for(Int_t nos2 = 0; nos2 < Tot_Scru2; nos2++)
          {
           gMC->Gspos("SCRL",nos2+1,"RB1L",xpos2_s,ypos2_s,zpos2_s,0,"ONLY");
           ypos2_s += 5.;
          }

       gMC->Gspos("RB1L",1, "FR1L", xpos_rb2, ypos_rb2, zpos_rb2,0, "only");
       gMC->Gspos("FR1L",1, "FRHL", xpos_fr2, ypos_fr2, zpos_fr2,0, "only");

//......................................................................................
//Frame-3

      Float_t frame3[3] ;
      frame3[0] = 4./2.;
      frame3[1] = 40.67/2.;
      frame3[2] = frame1[2];
      Float_t bend_ang = (22.5*PI/180.);  //bending angle of frame-3 w.r.t frame-2

      Float_t xpos_fr3 =-frame3[1]* sin(bend_ang);
      Float_t ypos_fr3 = frame2[1]+(frame3[1]*cos(bend_ang)-frame3[0]* sin(bend_ang));
      Float_t zpos_fr3 = 0.;

  //Frame part extended inside
      Float_t fr_ex[3];
      fr_ex[0] = 3.53/2.;
      fr_ex[1] = 56.28/2.;
      fr_ex[2] = frame1[2];

      Float_t xpos_ex = -frame3[0] - fr_ex[0];
      Float_t ypos_ex = 0.0;
      Float_t zpos_ex = 0.0;

      Float_t rib3[3];
      rib3[0] = 1./2.;
      rib3[1] = frame3[1];
      rib3[2] = 1.84/2.;

      Float_t xpos_rb3 = frame3[0] - rib3[0];
      Float_t ypos_rb3 = 0.0;
      Float_t zpos_rb3 = frame3[2] + rib3[2];

      gMC->Gsvolu("FR2L", "BOX", idPGF30, frame3, 3); //Frame - 3
      gMC->Gsvolu("FEXL", "BOX", idPGF30, fr_ex, 3); //frame- extended part
      gMC->Gsvolu("RB2L", "BOX", idPGF30, rib3, 3); //Rib - 3

      Float_t xpos3_s = -rib3[0] + 0.4;
      Float_t ypos3_s = -rib3[1] + 1.1;
      Float_t zpos3_s = -frame3[2];

      Int_t Tot_Scru3 = Tot_Scru2 + 9;      //Toal screw on this rib is 9--index continues
      for(Int_t nos3 = Tot_Scru2; nos3 < Tot_Scru3; nos3++)
         {
          gMC->Gspos("SCRL",nos3+1,"RB2L",xpos3_s,ypos3_s,zpos3_s,0,"ONLY");
          ypos3_s += 5.;
         }

      gMC->Gspos("FEXL",1, "FR2L", xpos_ex, ypos_ex, zpos_ex,0, "only");
      gMC->Gspos("RB2L",1, "FR2L", xpos_rb3, ypos_rb3, zpos_rb3,0, "only");
      gMC->Gspos("FR2L",1, "FR1L", xpos_fr3, ypos_fr3, zpos_fr3,irot6, "only");
//......................................................................................
//Frame-4

      Float_t frame4[3] ;
      frame4[0] = 4./2;
      frame4[1] = 52.54/2.;
      frame4[2] = frame1[2];

      Float_t xpos_fr4 =-frame4[1]* sin(bend_ang);
      Float_t ypos_fr4 = frame3[1]+(frame4[1]*cos(bend_ang)-frame4[0]* sin(bend_ang));
      Float_t zpos_fr4 =0.;

      Float_t rib4[3];
      rib4[0] = 1./2.;
      rib4[1] = frame4[1];
      rib4[2] = 1.84/2.;

      Float_t xpos_rb4 = frame4[0] - rib4[0];
      Float_t ypos_rb4 = 0.0;
      Float_t zpos_rb4 = frame4[2] + rib4[2];

      gMC->Gsvolu("FR3L", "BOX", idPGF30, frame4, 3); //Frame - 4
      gMC->Gsvolu("RB3L", "BOX", idPGF30, rib4, 3); //Rib - 4

      Float_t xpos4_s = -rib4[0] + 0.4;
      Float_t ypos4_s = -rib4[1] + 4.33;
      Float_t zpos4_s = -frame4[2];

      Int_t Tot_Scru4 = Tot_Scru3 + 10;  //Toal screw on this rib is 10--index continues
      for(Int_t nos4 = Tot_Scru3; nos4 < Tot_Scru4; nos4++)
         {
          gMC->Gspos("SCRL",nos4+1,"RB3L",xpos4_s,ypos4_s,zpos4_s,0,"ONLY");
          ypos4_s += 5.;
         }

     gMC->Gspos("RB3L",1, "FR3L", xpos_rb4, ypos_rb4, zpos_rb4,0, "only");
     gMC->Gspos("FR3L",1, "FR2L", xpos_fr4, ypos_fr4, zpos_fr4,irot6, "only");
//......................................................................................
//Frame-5

      Float_t frame5[3] ;
      frame5[0] = 4./2;
      frame5[1] = 41.83/2;
      frame5[2] = frame1[2];

      Float_t xpos_fr5 =-frame5[1]*sin(bend_ang);
      Float_t ypos_fr5 = frame4[1]+(frame5[1]*cos(bend_ang)-frame5[0]* sin(bend_ang));
      Float_t zpos_fr5 = 0.;

      Float_t rib5[3];
      rib5[0] = 1./2.;
      rib5[1] = frame5[1];
      rib5[2] = 1.84/2.;

      Float_t xpos_rb5 = frame5[0] - rib5[0];
      Float_t ypos_rb5 = 0.0;
      Float_t zpos_rb5 = frame5[2] + rib5[2];

      gMC->Gsvolu("FR4L", "BOX", idPGF30, frame5, 3); //Frame - 5
      gMC->Gsvolu("RB4L", "BOX", idPGF30, rib5, 3); //Rib - 5

      Float_t xpos5_s = -rib5[0] + 0.4;
      Float_t ypos5_s = -rib5[1] + 2.79;
      Float_t zpos5_s = -frame5[2];

      Int_t Tot_Scru5 = Tot_Scru4 + 8; //Toal screw on this rib is 8--index continues
      for(Int_t nos5 = Tot_Scru4; nos5 < Tot_Scru5; nos5++)
         {
          gMC->Gspos("SCRL",nos5+1,"RB4L",xpos5_s,ypos5_s,zpos5_s,0,"ONLY");
          ypos5_s += 5.;
         }

      gMC->Gspos("RB4L",1, "FR4L", xpos_rb5, ypos_rb5, zpos_rb5,0, "only");
      gMC->Gspos("FR4L",1, "FR3L", xpos_fr5, ypos_fr5, zpos_fr5,irot6, "only");
//......................................................................................
//Frame-6

      Float_t frame6[3] ;
      frame6[0] = 3./2.;
      frame6[1] = 30.84/2.;
      frame6[2] = frame1[2];

      Float_t xpos_fr6 = (-frame6[1]*sin(bend_ang))+0.5;
      Float_t ypos_fr6 = frame5[1]+(frame6[1]*cos(bend_ang)-frame6[0]* sin(bend_ang));
      Float_t zpos_fr6 = 0;

      Float_t rib6[3];
      rib6[0] = 1./2.;
      rib6[1] = frame6[1]+0.8;
      rib6[2] = 1.84/2.;

      Float_t xpos_rb6 = frame6[0] - rib6[0];
      Float_t ypos_rb6 = 0.4;
      Float_t zpos_rb6 = frame6[2] + rib6[2];

      gMC->Gsvolu("FR5L", "BOX", idPGF30, frame6, 3); //Frame - 6
      gMC->Gsvolu("RB5L", "BOX", idPGF30, rib6, 3); //Rib - 6

      Float_t xpos6_s = -rib6[0] + 0.4;
      Float_t ypos6_s = -rib6[1] + 1.36;
      Float_t zpos6_s = -frame6[2];

      Int_t Tot_Scru6 = Tot_Scru5 + 6;      //Toal screw on this rib is 7--index continues
      for(Int_t nos6 = Tot_Scru5; nos6 < Tot_Scru6; nos6++)
         {
          gMC->Gspos("SCRL",nos6+1,"RB5L",xpos6_s,ypos6_s,zpos6_s,0,"ONLY");
          ypos6_s += 5.;
         }

      gMC->Gspos("RB5L",1, "FR5L", xpos_rb6, ypos_rb6, zpos_rb6,0, "only");
      gMC->Gspos("FR5L",1, "FR4L", xpos_fr6, ypos_fr6, zpos_fr6,irot6, "only");
//......................................................................................
//Frame-7 Vertical frame

      Float_t frame7[3] ;
      frame7[0] = 2.7/2.;
      frame7[1] = 97.9/2.;
      frame7[2] = .95/2.;

      Float_t rib7[3];
      rib7[0] = 1./2.;
      rib7[1] = frame7[1] - .9/2.;
      rib7[2] = 1.84/2.;

      Float_t xpos_fr7 =  frame7[0] -3.7;
      Float_t ypos_fr7 =  frame7[1] +20.6;
      Float_t zpos_fr7 = frame7[2];

      Float_t xpos_rb7 = -frame7[0] + rib7[0]+0.9;
      Float_t ypos_rb7 = .9/2.;
      Float_t zpos_rb7 =  frame7[2] + rib7[2];

      gMC->Gsvolu("FRVL", "BOX", idPGF30, frame7, 3); //Frame - vertical
      gMC->Gsvolu("RBVL", "BOX", idPGF30, rib7, 3); //Rib

  //Fixing Screws-- screw parameter and screw-head are taken from horizontal frame bar
      gMC->Gsvolu("SCVL","TUBE",idScru, spar,3);    //screw-vertical part for Frame M4 x 25
  //positioning head over screws
      gMC->Gspos("HDFL",3,"SCVL",xpos_h,ypos_h,zpos_h,0,"ONLY");//positioning Screw-head
 //Screws on Frame
      Float_t xposv_s = -frame7[0] + 0.4;
      Float_t yposv_s = -frame7[1] + 0.4;
      Float_t zposv_s =  0.0;

      Int_t Tot_Scru_v = 20;   //Total no. of Screws
      for(Int_t nos_v = 0; nos_v < Tot_Scru_v; nos_v++)
         {
          gMC->Gspos("SCVL",nos_v+1,"FRVL",xposv_s,yposv_s,zposv_s,0,"ONLY");
          yposv_s += 5.;
	 } // nos2 is no. of scru

      gMC->Gspos("RBVL",1, "FRVL", xpos_rb7, ypos_rb7, zpos_rb7,0, "only");
//    gMC->Gspos("FRVL",1, "SQM3", xpos_fr7, ypos_fr7, zpos_fr7,0, "only");// frame vertical
//    gMC->Gspos("FRVL",1, "SQM4", xpos_fr7, ypos_fr7, zpos_fr7,0, "only");// frame vertical

//......................................................................................
//Frame - inner
      Float_t fr_c[5];   //semi-circular frame
      fr_c[0] = 20.6;
      fr_c[1] = 23.1;
      fr_c[2] = 0.95/2.;
      fr_c[3] = -3.5;
      fr_c[4] = 93.5;

      Float_t xpos_frc = 0.0;
      Float_t ypos_frc = 0.0;
      Float_t zpos_frc = fr_c[2];

      Float_t rib_c[5];
      rib_c[0] = 21.5;
      rib_c[1] = 22.5;
      rib_c[2] = 1.84/2.;
      rib_c[3] = fr_c[3]-2.0;
      rib_c[4] = fr_c[4]+2.0;

      Float_t xpos_rc = 0.0;
      Float_t ypos_rc = 0.0;
      Float_t zpos_rc = fr_c[2] + rib_c[2];

      gMC->Gsvolu("FRCL", "TUBS", idPGF30, fr_c, 5); //Frame - semi circular
      gMC->Gsvolu("RBCL", "TUBS", idPGF30, rib_c, 5); //Rib

//Screws
      gMC->Gsvolu("SCYL","TUBE",idScru, spar,3);    //screw-vertical part for extended part in -Y
      gMC->Gsvolu("SCIL","TUBE",idScru, spar,3);    //screw-vertical part
      gMC->Gsvolu("SCXL","TUBE",idScru, spar,3);    //screw-vertical part for extended part in -X
      gMC->Gspos("HDFL",3,"SCIL",xpos_h,ypos_h,zpos_h,0,"ONLY");//positioning Screw-head

    // on circular part
      Float_t zpos_is2 = 0.0;
      Float_t theta[7];
      Float_t radius = fr_c[0] + 0.4 ;   //inner radius + 0.4
      Float_t arc = 3.667;  // for 10-degree angle
      for(Int_t i = 0; i<8; i++)
         {
	  theta[i] = arc/radius;
	  Float_t xpos_is2 = radius * cos(theta[i]);
	  Float_t ypos_is2 = radius * sin(theta[i]);
          gMC->Gspos("SCIL",i+1,"FRCL",xpos_is2, ypos_is2, zpos_is2,0,"ONLY");
          arc +=3.667;
	 }



      gMC->Gspos("RBCL",1, "FRCL", xpos_rc, ypos_rc, zpos_rc,0, "only");  //Rib

       gMC->Gspos("FRHL",1, "SQM3", xpos_fr1, ypos_fr1, zpos_fr1,0, "only");// frame -1
       gMC->Gspos("FRHL",1, "SQM4", xpos_fr1, ypos_fr1, zpos_fr1,0, "only");// frame -1

       gMC->Gspos("FRVL",1, "SQM3", xpos_fr7, ypos_fr7, zpos_fr7,0, "only");// frame vertical
       gMC->Gspos("FRVL",1, "SQM4", xpos_fr7, ypos_fr7, zpos_fr7,0, "only");// frame vertical

       gMC->Gspos("FRCL",1, "SQM3", xpos_frc, ypos_frc, zpos_frc,0, "only");// frame semi circular
       gMC->Gspos("FRCL",1, "SQM4", xpos_frc, ypos_frc, zpos_frc,0, "only");// frame semi circular

//=============================================================================================

//                                   Plane - 2

//=============================================================================================

 //Cathode PCB + Copper sheet over PCB + Roha cell over copper sheet

//Segment - 0
       bpar_h[0] = 95.5/2.;
       bpar_h[1] = 1.6/2.;
       bpar_h[2] = zcbb/2.;
       gMC->Gsvolu("CB0R", "BOX", idPCB, bpar_h, 3);

       bpar_h[2] = zcu/2.;     //Thickness of Copper sheet
       gMC->Gsvolu("CU0R", "BOX", idCU, bpar_h, 3);

       bpar_h[2] = zRoha/2.;     //Thickness of Roha cell
       gMC->Gsvolu("RH0R", "BOX", idRoha, bpar_h, 3);

       bpar_h[0] = (100.6/2)-(0.9/2);
       bpar_h[1] = 2.8/2.;
       bpar_h[2] = zmeb/2.;
       gMC->Gsvolu("MB0R", "BOX", idPCB, bpar_h, 3);

       bpar_h[2] = zcu2/2;
       gMC->Gsvolu("EB0R", "BOX", idCU, bpar_h, 3);

//Segment-1
	pgpar[1] = 13.;
	pgpar[5] = 22.5;
	pgpar[6] = 120.2;
        pgpar[7] = pgpar[4] - zcbb;
	pgpar[8] = pgpar[5];
	pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB1R", "PGON", idPCB, pgpar, 10);

	pgpar[7] = pgpar[4] - zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU1R", "PGON", idCU, pgpar, 10);

	pgpar[7] = pgpar[4] - zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH1R", "PGON", idRoha, pgpar, 10);

	pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 121.2;
        pgpar[7] = pgpar[4] - zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB1R", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] - zcu2;
        gMC->Gsvolu("EB1R", "PGON", idCU, pgpar, 10);

//Segment-2

	pgpar[1] = 21.;
	pgpar[5] = 22.5;
	pgpar[6] = 121.5;
        pgpar[7] = pgpar[4] - zcbb;
	pgpar[8] = pgpar[5];
	pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB2R", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] - zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU2R", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] - zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH2R", "PGON", idRoha, pgpar, 10);

	pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 122.5;
        pgpar[7] = pgpar[4] - zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB2R", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] - zcu2;
        gMC->Gsvolu("EB2R", "PGON", idCU, pgpar, 10);

//Segment-3

	pgpar[1] = 22.;
	pgpar[5] = 22.5;
	pgpar[6] = 119.9;
        pgpar[7] = pgpar[4] - zcbb;
	pgpar[8] = pgpar[5];
	pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB3R", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] - zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU3R", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] - zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH3R", "PGON", idRoha, pgpar, 10);

        pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 120.9;
	pgpar[7] = pgpar[4] - zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB3R", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] - zcu2;
        gMC->Gsvolu("EB3R", "PGON", idCU, pgpar, 10);

//Segment-4

	pgpar[1] = 20.;
	pgpar[5] = 22.5;
	pgpar[6] = 119.9;
        pgpar[7] = pgpar[4] - zcbb;
	pgpar[8] = pgpar[5];
	pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB4R", "PGON", idPCB, pgpar, 10);

	pgpar[7] = pgpar[4] - zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU4R", "PGON", idCU, pgpar, 10);

	pgpar[7] = pgpar[4] - zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH4R", "PGON", idRoha, pgpar, 10);

	pgpar[5] = 21.5; // radius of inner rib
        pgpar[6] = 120.9;
	pgpar[7] = pgpar[4] - zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB4R", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] - zcu2;
        gMC->Gsvolu("EB4R", "PGON", idCU, pgpar, 10);

//Segment-5

	pgpar[1] = 14.;
	pgpar[5] = 22.5;
	pgpar[6] = 117.5;
        pgpar[7] = pgpar[4] - zcbb;
	pgpar[8] = pgpar[5];
	pgpar[9] = pgpar[6];
        gMC->Gsvolu("CB5R", "PGON", idPCB, pgpar, 10);

	pgpar[7] = pgpar[4] - zcu;  // Thickness of copper-sheet
        gMC->Gsvolu("CU5R", "PGON", idCU, pgpar, 10);

        pgpar[7] = pgpar[4] - zRoha;  // Thickness of Roha cell
        gMC->Gsvolu("RH5R", "PGON", idRoha, pgpar, 10);

        pgpar[5] = 21.5;// radius of inner rib
        pgpar[6] = 118.5;
        pgpar[7] = pgpar[4] - zmeb;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("MB5R", "PGON", idPCB, pgpar, 10);

        pgpar[7] = pgpar[4] - zcu2;
        gMC->Gsvolu("EB5R", "PGON", idCU, pgpar, 10);

//Segment-6 - vertical box

        bpar_v[0] = 1.6/2.;
        bpar_v[1] = 95.5/2.;
        bpar_v[2] = zcbb/2.;
        gMC->Gsvolu("CB6R", "BOX", idPCB, bpar_v, 3);
        bpar_v[2] = zcu/2.;
        gMC->Gsvolu("CU6R", "BOX", idCU, bpar_v, 3);

        bpar_v[2] = zRoha/2.;
        gMC->Gsvolu("RH6R", "BOX", idRoha, bpar_v, 3);

        bpar_v[0] = 2.8/2.;
        bpar_v[1] = (97.9/2)-(0.9/2.);
        bpar_v[2] = zmeb/2;
        gMC->Gsvolu("MB6R", "BOX", idPCB, bpar_v, 3);

        bpar_v[2] = zcu2/2;
        gMC->Gsvolu("EB6R", "BOX", idCU, bpar_v, 3);


//...........................................................................................
//Positioning of Electronic exit board
     gMC->Gspos("EB0R",1, "MB0R", 0.,0.,-zpos_eeb_bar,0, "only");
     gMC->Gspos("EB1R",1, "MB1R", 0.,0.,-zpos_cu2,0, "only");
     gMC->Gspos("EB2R",1, "MB2R", 0.,0.,-zpos_cu2,0, "only");
     gMC->Gspos("EB3R",1, "MB3R", 0.,0.,-zpos_cu2,0, "only");
     gMC->Gspos("EB4R",1, "MB4R", 0.,0.,-zpos_cu2,0, "only");
     gMC->Gspos("EB5R",1, "MB5R", 0.,0.,-zpos_cu2,0, "only");
     gMC->Gspos("EB6R",1, "MB6R", 0.,0.,-zpos_eeb_bar,0, "only");

//Positioning of Mech. exit board
     xpos_meb_hor =  1.1;
     ypos_meb_hor = -0.6;
     xpos_meb_ver = -0.6;
     ypos_meb_ver = -0.35;

     gMC->Gspos("MB0R",1, "RH0R", xpos_meb_hor,ypos_meb_hor,-zpos_meb_bar,0, "only");
     gMC->Gspos("MB1R",1, "RH1R", 0.,0.,-zpos_meb,0, "only");
     gMC->Gspos("MB2R",1, "RH2R", 0.,0.,-zpos_meb,0, "only");
     gMC->Gspos("MB3R",1, "RH3R", 0.,0.,-zpos_meb,0, "only");
     gMC->Gspos("MB4R",1, "RH4R", 0.,0.,-zpos_meb,0, "only");
     gMC->Gspos("MB5R",1, "RH5R", 0.,0.,-zpos_meb,0, "only");
     gMC->Gspos("MB6R",1, "RH6R", xpos_meb_ver,ypos_meb_ver,-zpos_meb_bar,0, "only");

//Positioning Roha cell over copper sheet
     gMC->Gspos("RH0R",1, "CU0R", 0.,0.,-zpos_Roha_bar,0, "only");// box horizontal
     gMC->Gspos("RH1R",1, "CU1R", 0.,0.,-zpos_Roha,0, "only");
     gMC->Gspos("RH2R",1, "CU2R", 0.,0.,-zpos_Roha,0, "only");
     gMC->Gspos("RH3R",1, "CU3R", 0.,0.,-zpos_Roha,0, "only");
     gMC->Gspos("RH4R",1, "CU4R", 0.,0.,-zpos_Roha,0, "only");
     gMC->Gspos("RH5R",1, "CU5R", 0.,0.,-zpos_Roha,0, "only");
     gMC->Gspos("RH6R",1, "CU6R", 0.,0.,-zpos_Roha_bar,0, "only");

//Positioning Copper sheet over PCB
     gMC->Gspos("CU0R",1, "CB0R", 0.,0.,-zpos_cubar,0, "only");// box horizontal
     gMC->Gspos("CU1R",1, "CB1R", 0.,0.,-zpos_cu,0, "only");
     gMC->Gspos("CU2R",1, "CB2R", 0.,0.,-zpos_cu,0, "only");
     gMC->Gspos("CU3R",1, "CB3R", 0.,0.,-zpos_cu,0, "only");
     gMC->Gspos("CU4R",1, "CB4R", 0.,0.,-zpos_cu,0, "only");
     gMC->Gspos("CU5R",1, "CB5R", 0.,0.,-zpos_cu,0, "only");
     gMC->Gspos("CU6R",1, "CB6R", 0.,0.,-zpos_cubar,0, "only");


//Positioning the PCB
     x_hor_pos = 95.5/2.+ 22.6; // 22.6 is inner radius of PCB
     y_hor_pos = -(1.6/2); //

     x_ver_pos = -(1.6/2); //
     y_ver_pos = 95.5/2.+ 22.6;

     gMC->Gspos("CB0R",1, "SQM3", x_hor_pos,-bpar_h[1],-zpos_cbb_bar,0, "only");// box horizontal
     gMC->Gspos("CB1R",1, "SQM3", 0.,0.,-zpos_cbb,0, "only");
     gMC->Gspos("CB2R",1, "SQM3", 0.,0.,-zpos_cbb,irot1,"only");
     gMC->Gspos("CB3R",1, "SQM3", 0.,0.,-zpos_cbb,irot2, "only");
     gMC->Gspos("CB4R",1, "SQM3", 0.,0.,-zpos_cbb,irot3, "only");
     gMC->Gspos("CB5R",1, "SQM3", 0.,0.,-zpos_cbb,irot4, "only");
     gMC->Gspos("CB6R",1, "SQM3", x_ver_pos,y_ver_pos,-zpos_cbb_bar,0, "only");// box vertical

     gMC->Gspos("CB0R",1, "SQM4", x_hor_pos,-bpar_h[1],-zpos_cbb_bar,0, "only");// box horizontal
     gMC->Gspos("CB1R",1, "SQM4", 0.,0.,-zpos_cbb,0, "only");
     gMC->Gspos("CB2R",1, "SQM4", 0.,0.,-zpos_cbb,irot1,"only");
     gMC->Gspos("CB3R",1, "SQM4", 0.,0.,-zpos_cbb,irot2, "only");
     gMC->Gspos("CB4R",1, "SQM4", 0.,0.,-zpos_cbb,irot3, "only");
     gMC->Gspos("CB5R",1, "SQM4", 0.,0.,-zpos_cbb,irot4, "only");
     gMC->Gspos("CB6R",1, "SQM4", x_ver_pos,y_ver_pos,-zpos_cbb_bar,0, "only");// box vertical

//----------------------------------------------------------------------
//                          Frames P2
//----------------------------------------------------------------------
//Frame-1 P2

       ypos_s = -frame1[1] + 0.4;
       zpos_fr1 = -frame1[2];

       zpos_rb1 = -(frame1[2] + rib1[2]);

       gMC->Gsvolu("FRHR", "BOX", idPGF30, frame1, 3); //Frame - 1 P2
       gMC->Gsvolu("RBHR", "BOX", idPGF30, rib1, 3);   //Rib - 1 P2

   //Fixing Screws ...........................................

       zpos_h = -(spar[2] + s_head[2]);

       gMC->Gsvolu("SCHR","TUBE",idScru, spar,3);    //screw-vertical part for Frame M4 x 25
       gMC->Gsvolu("HDFR","TUBE",idScru, s_head,3);    //screw-head
   //positioning head over screws
       gMC->Gspos("HDFR",1,"SCHR",xpos_h,ypos_h,zpos_h,0,"ONLY");//positioning Screw-head

   //Screws on Frame

       xpos_s = frame1[0] - 0.4;

       Tot_Scru1 = 21;   //Total no. of Scru
       for(Int_t nos1 = 0; nos1 < Tot_Scru1; nos1++)
          {
           gMC->Gspos("SCHR",nos1+1,"FRHR",xpos_s,ypos_s,zpos_s,irot1,"ONLY");
           xpos_s -= 5.;
          } // nos1 is no. of scru

       gMC->Gspos("RBHR",1,"FRHR",xpos_rb1,ypos_rb1,zpos_rb1,0, "only");
       gMC->Gspos("FRHR",1,"SQM3",xpos_fr1,ypos_fr1,zpos_fr1,0, "only");
       gMC->Gspos("FRHR",1,"SQM4",xpos_fr1,ypos_fr1,zpos_fr1,0, "only");

//......................................................................................

//Frame-2

      zpos_rb2 = -(frame2[2] + rib2[2]);

      gMC->Gsvolu("FR1R", "BOX", idPGF30, frame2, 3); //Frame - 2
      gMC->Gsvolu("RB1R", "BOX", idPGF30, rib2, 3);   //Rib - 2

//Fixing Screws ...........................................
//---------screw parameters (on Ribs) M4

       spar2[2] = 2.9/2.;
       gMC->Gsvolu("SCRR","TUBE",idScru, spar2,3); //screw-vertical part for Rib m4
   //positioning head over screws

       zpos2_h = -(spar2[2] + s_head[2]);
       gMC->Gspos("HDFR",2,"SCRR",xpos2_h,ypos2_h,zpos2_h,0,"ONLY");//positioning Screw-head

       ypos2_s = -rib2[1] + 1.7;
       zpos2_s = frame2[2];

       Tot_Scru2 = 6;   // 6 screws
       for(Int_t nos2 = 0; nos2 < Tot_Scru2; nos2++)
          {
           gMC->Gspos("SCRR",nos2+1,"RB1R",xpos2_s,ypos2_s,zpos2_s,0,"ONLY");
           ypos2_s += 5.;
          }

       gMC->Gspos("RB1R",1, "FR1R", xpos_rb2, ypos_rb2, zpos_rb2,0, "only");
       gMC->Gspos("FR1R",1, "FRHR", xpos_fr2, ypos_fr2, zpos_fr2,0, "only");

//......................................................................................

//Frame-3 P2

       zpos_rb3 = -(frame3[2] + rib3[2]);

       gMC->Gsvolu("FR2R", "BOX", idPGF30, frame3, 3); //Frame - 3
       gMC->Gsvolu("FEXR", "BOX", idPGF30, fr_ex, 3); //frame- extended part
       gMC->Gsvolu("RB2R", "BOX", idPGF30, rib3, 3); //Rib - 3

   //Fixing Screws ...........................................
        ypos3_s = -rib3[1] + 1.1;
        zpos3_s = frame3[2];

        Tot_Scru3 = Tot_Scru2 + 9;      //Toal screw on this rib is 9--index continues
        for(Int_t nos3 = Tot_Scru2; nos3 < Tot_Scru3; nos3++)
           {
            gMC->Gspos("SCRR",nos3+1,"RB2R",xpos3_s,ypos3_s,zpos3_s,0,"ONLY");
            ypos3_s += 5.;
           }

        gMC->Gspos("FEXR",1, "FR2R", xpos_ex, ypos_ex, zpos_ex,0, "only");
        gMC->Gspos("RB2R",1, "FR2R", xpos_rb3, ypos_rb3, zpos_rb3,0, "only");
        gMC->Gspos("FR2R",1, "FR1R", xpos_fr3, ypos_fr3, zpos_fr3,irot6, "only");
//......................................................................................

//Frame-4 P2

        zpos_rb4 = -(frame4[2] + rib4[2]);

        gMC->Gsvolu("FR3R", "BOX", idPGF30, frame4, 3); //Frame - 4
        gMC->Gsvolu("RB3R", "BOX", idPGF30, rib4, 3); //Rib - 4

   //Fixing Screws ...........................................

        ypos4_s = -rib4[1] + 4.33;
        zpos4_s = frame4[2];

        Tot_Scru4 = Tot_Scru3 + 10;      //Toal screw on this rib is 10--index continues
        for(Int_t nos4 = Tot_Scru3; nos4 < Tot_Scru4; nos4++)
           {
            gMC->Gspos("SCRR",nos4+1,"RB3R",xpos4_s,ypos4_s,zpos4_s,0,"ONLY");
            ypos4_s += 5.;
           }

        gMC->Gspos("RB3R",1, "FR3R", xpos_rb4, ypos_rb4, zpos_rb4,0, "only");
        gMC->Gspos("FR3R",1, "FR2R", xpos_fr4, ypos_fr4, zpos_fr4,irot6, "only");
//......................................................................................
//Frame-5 P2

        zpos_rb5 = -(frame5[2] + rib5[2]);

        gMC->Gsvolu("FR4R", "BOX", idPGF30, frame5, 3); //Frame - 5
        gMC->Gsvolu("RB4R", "BOX", idPGF30, rib5, 3); //Rib - 5

     //Fixing Screws ...........................................

        ypos5_s = -rib5[1] + 2.79;
        zpos5_s = frame5[2];

        Tot_Scru5 = Tot_Scru4 + 8;      //Toal screw on this rib is 8--index continues
        for(Int_t nos5 = Tot_Scru4; nos5 < Tot_Scru5; nos5++)
           {
            gMC->Gspos("SCRR",nos5+1,"RB4R",xpos5_s,ypos5_s,zpos5_s,0,"ONLY");
            ypos5_s += 5.;
           }

        gMC->Gspos("RB4R",1, "FR4R", xpos_rb5, ypos_rb5, zpos_rb5,0, "only");
        gMC->Gspos("FR4R",1, "FR3R", xpos_fr5, ypos_fr5, zpos_fr5,irot6, "only");
//......................................................................................
//Frame-6 P2

        zpos_rb6 = -(frame6[2] + rib6[2]);

        gMC->Gsvolu("FR5R", "BOX", idPGF30, frame6, 3); //Frame - 6
        gMC->Gsvolu("RB5R", "BOX", idPGF30, rib6, 3); //Rib - 6

     //Fixing Screws ...........................................

        ypos6_s = -rib6[1] + 1.36;
        zpos6_s = frame6[2];

        Tot_Scru6 = Tot_Scru5 + 6;      //Toal screw on this rib is 7--index continues
        for(Int_t nos6 = Tot_Scru5; nos6 < Tot_Scru6; nos6++)
           {
            gMC->Gspos("SCRR",nos6+1,"RB5R",xpos6_s,ypos6_s,zpos6_s,0,"ONLY");
            ypos6_s += 5.;
           }

        gMC->Gspos("RB5R",1, "FR5R", xpos_rb6, ypos_rb6, zpos_rb6,0, "only");
        gMC->Gspos("FR5R",1, "FR4R", xpos_fr6, ypos_fr6, zpos_fr6,irot6, "only");
//......................................................................................

//Frame-7 Vertical frame P2

        zpos_fr7 = -frame7[2];
        zpos_rb7 =  -(frame7[2] + rib7[2]);

        gMC->Gsvolu("FRVR", "BOX", idPGF30, frame7, 3); //Frame - vertical
        gMC->Gsvolu("RBVR", "BOX", idPGF30, rib7, 3); //Rib

   //Fixing Screws-- screw parameter and screw-head are taken from horizontal frame bar
        gMC->Gsvolu("SCVR","TUBE",idScru, spar,3);    //screw-vertical part for Frame M4
   //positioning head over screws
        gMC->Gspos("HDFR",3,"SCVR",xpos_h,ypos_h,zpos_h,0,"ONLY");//positioning Screw-head

       //Screws on Frame
       yposv_s = -frame7[1] + 0.4;
       zposv_s =  0.0;

       Tot_Scru_v = 20;   //Total no. of Screws
       for(Int_t nos_v = 0; nos_v < Tot_Scru_v; nos_v++)
          {
           gMC->Gspos("SCVR",nos_v+1,"FRVR",xposv_s,yposv_s,zposv_s,0,"ONLY");
           yposv_s += 5.;
          }

       gMC->Gspos("RBVR",1, "FRVR", xpos_rb7, ypos_rb7, zpos_rb7,0, "only");
       gMC->Gspos("FRVR",1, "SQM3", xpos_fr7, ypos_fr7, zpos_fr7,0, "only");// frame vertical
       gMC->Gspos("FRVR",1, "SQM4", xpos_fr7, ypos_fr7, zpos_fr7,0, "only");// frame vertical

//......................................................................................
//Frame - inner

       zpos_frc = - fr_c[2];
       zpos_rc = -(fr_c[2] + rib_c[2]);

       gMC->Gsvolu("FRCR", "TUBS", idPGF30, fr_c, 5); //Frame - semi circular
       gMC->Gsvolu("RBCR", "TUBS", idPGF30, rib_c, 5); //Rib

   //Screws -------------------------------------------------------
       gMC->Gsvolu("SCYR","TUBE",idScru, spar,3);    //screw-vertical part for extended part in -Y
       gMC->Gsvolu("SCIR","TUBE",idScru, spar,3);    //screw-vertical part
       gMC->Gsvolu("SCXR","TUBE",idScru, spar,3);    //screw-vertical part for extended part in -X
       gMC->Gspos("HDFR",3,"SCIR",xpos_h,ypos_h,zpos_h,0,"ONLY");//positioning Screw-head

   // on circular part
      radius = fr_c[0] + 0.4 ;   //inner radius + 0.4
      zpos_is2 = 0.0;
      arc = 3.667;  // for 10-degree angle
      for(Int_t i = 0; i<8; i++)
         {
	  theta[i] = arc/radius;
	  Float_t xpos_is2 = radius * cos(theta[i]);
	  Float_t ypos_is2 = radius * sin(theta[i]);
          gMC->Gspos("SCIR",i+1,"FRCR",xpos_is2, ypos_is2, zpos_is2,0,"ONLY");
          arc +=3.667;
         }

      gMC->Gspos("RBCR",1, "FRCR", xpos_rc, ypos_rc, zpos_rc,0, "only");  //Rib
      gMC->Gspos("FRCR",1, "SQM3", xpos_frc, ypos_frc, zpos_frc,0, "only");// frame semi circular
      gMC->Gspos("FRCR",1, "SQM4", xpos_frc, ypos_frc, zpos_frc,0, "only");// frame semi circular


//Plane 2 -----------------------------------------------------------------------------------

//^^^^^^^^^^^^^^^^^^^^^^^^^ Sensitive volumes ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       Float_t zsenv = 0.5; // distance between two cathode plane

 //Segment-0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       bpar_h[0] = 95.5/2.;
       bpar_h[1] = 1.6/2.;
       bpar_h[2] = zsenv/2.;
       gMC->Gsvolu("C3G0", "BOX", idGas, bpar_h, 3);
       gMC->Gsvolu("C4G0", "BOX", idGas, bpar_h, 3);

 //Segment-1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 13.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = -zsenv/2.;
        pgpar[5] = 22.5;
        pgpar[6] = 117.2;
        pgpar[7] = zsenv/2.;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("C3G1", "PGON", idGas, pgpar, 10);
        gMC->Gsvolu("C4G1", "PGON", idGas, pgpar, 10);

//Segment-2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 21.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = -zsenv/2.;
        pgpar[5] = 22.5;
        pgpar[6] = 115.0;
        pgpar[7] = zsenv/2.;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("C3G2", "PGON", idGas, pgpar, 10);
        gMC->Gsvolu("C4G2", "PGON", idGas, pgpar, 10);

//Segment-3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 22.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = -zsenv/2.;
        pgpar[5] = 22.5;
        pgpar[6] = 116.9;
        pgpar[7] = zsenv/2.;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("C3G3", "PGON", idGas, pgpar, 10);
        gMC->Gsvolu("C4G3", "PGON", idGas, pgpar, 10);

//Segment-4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 20.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = -zsenv/2.;
        pgpar[5] = 22.5;
        pgpar[6] = 116.9;
        pgpar[7] = zsenv/2.;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("C3G4", "PGON", idGas, pgpar, 10);
        gMC->Gsvolu("C4G4", "PGON", idGas, pgpar, 10);

//Segment-5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 14.;
        pgpar[2] = 1.;
        pgpar[3] = 2.;
        pgpar[4] = -zsenv/2.;
        pgpar[5] = 22.5;
        pgpar[6] = 115.5;
        pgpar[7] = zsenv/2.;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("C3G5", "PGON", idGas, pgpar, 10);
        gMC->Gsvolu("C4G5", "PGON", idGas, pgpar, 10);

//Segment-6 - vertical box ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        bpar_v[0] = 2.16/2.;
	bpar_v[1] = 95.5/2.;
	bpar_v[2] = zsenv/2.;
        gMC->Gsvolu("C3G6", "BOX", idGas, bpar_v, 3);
        gMC->Gsvolu("C4G6", "BOX", idGas, bpar_v, 3);

//...........................................................................................

//Positioning the PCB
      x_hor_pos = 95.5/2.+ 22.6; // 22.6 is inner radius of PCB
      y_hor_pos = -(1.6/2); //

      x_ver_pos = -(1.6/2); //
      y_ver_pos = 95.5/2.+ 22.6;

     gMC->Gspos("C3G0",1, "SQM3", x_hor_pos,y_hor_pos,0.,0, "only");// box horizontal
     gMC->Gspos("C3G1",1, "SQM3", 0.,0.,0.,0, "only");
     gMC->Gspos("C3G2",1, "SQM3", 0.,0.,0.,irot1,"only");
     gMC->Gspos("C3G3",1, "SQM3", 0.,0.,0.,irot2, "only");
     gMC->Gspos("C3G4",1, "SQM3", 0.,0.,0.,irot3, "only");
     gMC->Gspos("C3G5",1, "SQM3", 0.,0.,0.,irot4, "only");
     gMC->Gspos("C3G6",1, "SQM3", x_ver_pos,y_ver_pos,0.,0, "only");// box vertical

     gMC->Gspos("C4G0",1, "SQM4", x_hor_pos,y_hor_pos,0.,0, "only");// box horizontal
     gMC->Gspos("C4G1",1, "SQM4", 0.,0.,0.,0, "only");
     gMC->Gspos("C4G2",1, "SQM4", 0.,0.,0.,irot1,"only");
     gMC->Gspos("C4G3",1, "SQM4", 0.,0.,0.,irot2, "only");
     gMC->Gspos("C4G4",1, "SQM4", 0.,0.,0.,irot3, "only");
     gMC->Gspos("C4G5",1, "SQM4", 0.,0.,0.,irot4, "only");
     gMC->Gspos("C4G6",1, "SQM4", x_ver_pos,y_ver_pos,0.,0, "only");// box vertical

//^^^^^^^^^^^^^^^^^^^^^^^^^ Sensitive volumes ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


//##################################################################################################
//   Positioning Quadrant  in chamber#3 and chamber#4
//##################################################################################################
/******Transformations for  Quadrant**********************************************
      ||        I  => Quadrant I:   no rotation
      ||
  II. || I.    II  => Quadrant II:  Reflaction of Quadrant I in XZ plane
      ||           => TGeoRotation("Qrot2",90.,180.,90.,90.,180.,0.);
=============
      ||       III => Quadrant III: 180 degree rotation of Quadrant I in XY plane
 III. || IV.       => TGeoRotation("Qrot3",90.,180.,90.,90.,180.,0.);
      ||        IV => Quadrant IV:-180 degree rotation of Quadrant II in XY plane
                   => TGeoRotation("Qrot4",90.,0.,90.,-90.,180.,0.);
**********************************************************************************************/

 Int_t detElemId1 =  1;  // quadrant I
 Int_t detElemId2 =  0;  // quadrant II
 Int_t detElemId3 =  3;  // quadrant III
 Int_t detElemId4 =  2;  // quadrant IV

 Float_t half_chamber =6.6/2;
 //Float_t half_chamber = zcbb + zcu + zRoha + zmeb + zcu2 + zsenv/2;

// ------------------------------St2 Chamber3------------------------------------------------

 //    GetEnvelopes(2)->AddEnvelope("S3M0", 300, true,TGeoTranslation(0.,0.,0.));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId1, 1, TGeoTranslation( 0., 0., - half_chamber));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId2, 2, TGeoTranslation( 0., 0., + half_chamber),
                                 TGeoRotation("Qrot3",90.,180.,90.,90.,180.,0.));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId3, 3, TGeoTranslation( 0., 0., - half_chamber),
                                  TGeoRotation("Qrot3",90.,180.,90.,270.,0.,0.));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId4, 4, TGeoTranslation( 0., 0., + half_chamber),
                                  TGeoRotation("Qrot4",90.,0.,90.,-90.,180.,0.));

//--------------------------------St2 Chamber4-------------------------------------------------

    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId1, 1, TGeoTranslation( 0., 0., - half_chamber));
    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId2, 2, TGeoTranslation( 0., 0., + half_chamber),
                                 TGeoRotation("Qrot2",90.,180.,90.,90.,180.,0.));
    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId3, 3, TGeoTranslation( 0., 0., - half_chamber),
                                  TGeoRotation("Qrot3",90.,180.,90.,270.,0.,0.));
    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId4, 4, TGeoTranslation( 0., 0., + half_chamber),
                                  TGeoRotation("Qrot4",90.,0.,90.,-90.,180.,0.));

//**********************************************************************************************

   }

//______________________________________________________________________________
void AliMUONSt2GeometryBuilderV2::SetTransformations()
{
// Defines the transformations for the station2 chambers.
// ---

  AliMUONChamber* iChamber1 = &fMUON->Chamber(2);
  Double_t zpos1 = - iChamber1->Z();
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  AliMUONChamber* iChamber2 = &fMUON->Chamber(3);
  Double_t zpos2 = - iChamber2->Z();
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));
}

//______________________________________________________________________________
void AliMUONSt2GeometryBuilderV2::SetSensitiveVolumes()
{
// Defines the sensitive volumes for station2 chambers.
// ---

  GetGeometry(2)->SetSensitiveVolume("C3G0");
  GetGeometry(2)->SetSensitiveVolume("C3G1");
  GetGeometry(2)->SetSensitiveVolume("C3G2");
  GetGeometry(2)->SetSensitiveVolume("C3G3");
  GetGeometry(2)->SetSensitiveVolume("C3G4");
  GetGeometry(2)->SetSensitiveVolume("C3G5");
  GetGeometry(2)->SetSensitiveVolume("C3G6");

  GetGeometry(3)->SetSensitiveVolume("C4G0");
  GetGeometry(3)->SetSensitiveVolume("C4G1");
  GetGeometry(3)->SetSensitiveVolume("C4G2");
  GetGeometry(3)->SetSensitiveVolume("C4G3");
  GetGeometry(3)->SetSensitiveVolume("C4G4");
  GetGeometry(3)->SetSensitiveVolume("C4G5");
  GetGeometry(3)->SetSensitiveVolume("C4G6");

}
