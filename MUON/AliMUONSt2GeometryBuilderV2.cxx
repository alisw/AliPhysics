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
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONConstants.h"

#define PI 3.14159

/// \cond CLASSIMP
ClassImp(AliMUONSt2GeometryBuilderV2)
/// \endcond

//______________________________________________________________________________
AliMUONSt2GeometryBuilderV2::AliMUONSt2GeometryBuilderV2(AliMUON* muon)
 : AliMUONVGeometryBuilder(2, 2),
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
AliMUONSt2GeometryBuilderV2::~AliMUONSt2GeometryBuilderV2() {
//
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSt2GeometryBuilderV2::CreateGeometry()
{

//
//********************************************************************
//                            Station 2                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers

     // Get tracking medias Ids
     Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;
     Int_t idAir  = idtmed[1100]; // medium 1
     Int_t idGas  = idtmed[1108]; // medium Ar-CO2 gas (80%+20%)
     Int_t idPCB  = idtmed[1122]; // medium FR4
     Int_t idCU   = idtmed[1110]; // medium copper
     Int_t idRoha = idtmed[1113]; // medium roha cell
     Int_t idPGF30= idtmed[1123]; // medium for Frame Eq.to Bakelite
     //Int_t idScru = idtmed[1128]; // screw material - Stainless Steel(18%Cr,9%Ni,Fe)



/*########################################################################################
    Create volume for one Quadrant
##########################################################################################*/
     Float_t tpar1[5];
     tpar1[0] = 20.55;
     tpar1[1] = 123.5;
     tpar1[2] = 6.8/2;
     tpar1[3] = -12.0;
     tpar1[4] = 102.0;


     gMC->Gsvolu("SQM3","TUBS", idAir, tpar1, 5);
     gMC->Gsvolu("SQM4","TUBS", idAir, tpar1, 5);



//==================================================================================
//                                 Plane      
//==================================================================================

//Thickness of variour parts
       Float_t zCbb  = 0.04;       //cathode pcb
       Float_t zCu   = 0.004;      // eff. cu in cathode pcb
       Float_t zRoha = 2.5;        // Rhocell
       Float_t zMeb = 0.08;      //Mech. exit board
       Float_t zEeb = 0.02;     //Effective electronic exit board
    //Z-positions of various parts--- in Plane-1

       Float_t zposCbb   = 0.25;   // 2.5 mm => gap between anode & chatode plane
       Float_t zposCu    =zposCbb + zCbb;
       Float_t zposRoha  =zposCu + zCu;
       Float_t zposMeb   =zposRoha + zRoha;
       Float_t zposEeb   =zposMeb + zMeb;

       Float_t zposCbbBar  = zposCbb + zCbb/2.;  //for segment 0 & 6
       Float_t zposCuBar   = zCbb/2. + zCu/2.;
       Float_t zposRohaBar = zCu/2. + zRoha/2.;
       Float_t zposMebBar  =  zRoha/2. + zMeb/2.;
       Float_t zposEebBar  = zMeb/2. + zEeb/2.;


 //Cathode PCB + Copper sheet + Rohacell + mech exit board + eff. electronic exit board

 //Segment-0 ~~~Horizantal box  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Float_t bparH[3];
       bparH[0] = 94.5/2.; // extension beyond 0 deg in x direction
       bparH[1] = 1.2/2.;  // extension beyond 0 deg in y direction (3.7[total extn] - 2.5[frame dim])
       bparH[2] = zCbb/2.; // thickness of cathode sheet in z direction
       gMC->Gsvolu("CB0L", "BOX", idPCB, bparH, 3);
       
       bparH[2] = zCu/2.;     //Thickness of Copper sheet
       gMC->Gsvolu("CU0L", "BOX", idCU, bparH, 3);

       bparH[2] = zRoha/2.;     //Thickness of Roha cell
       gMC->Gsvolu("RH0L", "BOX", idRoha, bparH, 3);

       bparH[2] = zMeb/2;       //Thickness of mechanical exit board
       gMC->Gsvolu("MB0L", "BOX", idPCB, bparH, 3);

       bparH[2] = zEeb/2;          //Thickness of effective electronic  exit board
       gMC->Gsvolu("EB0L", "BOX", idCU, bparH, 3);

 //Segment-1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Float_t pgpar[10]; // polygon
       pgpar[0] = 0.;  // initial angle
       pgpar[1] = 90.; // increment in angle starting from initial angle 
       pgpar[2] = 5;   // number of side
       pgpar[3] = 2.; // number of plane
       pgpar[4] = 0.; // z-position of the first plane
       pgpar[5] = 23.1; // innner radius first plane
       pgpar[6] = 117.6;  // outer radious first plane
       pgpar[7] = pgpar[4] + zCbb; // z-position of the second plane
       pgpar[8] = pgpar[5];  // innner radius of second plane
       pgpar[9] = pgpar[6];  // outer radious of second plane
       gMC->Gsvolu("CB1L", "PGON", idPCB, pgpar, 10);

       pgpar[7] = pgpar[4] + zCu;  // Thickness of copper-sheet
       gMC->Gsvolu("CU1L", "PGON", idCU, pgpar, 10);
       
       pgpar[7] = pgpar[4] + zRoha;  // Thickness of Roha cell
       gMC->Gsvolu("RH1L", "PGON", idRoha, pgpar, 10);
       
       pgpar[7] = pgpar[4] + zMeb;  // Thickness of mechanical exit board
       gMC->Gsvolu("MB1L", "PGON", idPCB, pgpar, 10);
       
       pgpar[7] = pgpar[4] + zEeb;  // Thickness of electronic exit board
       gMC->Gsvolu("EB1L", "PGON", idCU, pgpar, 10);
       

//Segment-2 - vertical box (simalar to horizontal bar as in Segment 0)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
       Float_t bparV[3];
       
       bparV[0] = 1.0/2.;
       bparV[1] = 94.5/2.;
       bparV[2] = zCbb/2.;
       gMC->Gsvolu("CB2L", "BOX", idPCB, bparV, 3);
       
       bparV[2] = zCu/2.;
       gMC->Gsvolu("CU2L", "BOX", idCU, bparV, 3);

       bparV[2] = zRoha/2.;
       gMC->Gsvolu("RH2L", "BOX", idRoha, bparV, 3);
       
       bparV[2] = zMeb/2;
       gMC->Gsvolu("MB2L", "BOX", idPCB, bparV, 3);
       
       bparV[2] = zEeb/2;
       gMC->Gsvolu("EB2L", "BOX", idCU, bparV, 3);
       
//....(Setting posion of Segment 0,1,2)..................................................................
       
       Float_t xposHorBox =  bparH[0] + 23.1; // 23.1 = 20.6(inner radius of qrd) + 2.5 (width of frame)
       Float_t yposHorBox = -bparH[1];  
       
       Float_t xposVerBox = -bparV[0];
       Float_t yposVerBox =  bparV[1] + 23.1;



//Positioning the PCB
//

      // chamber 3
     gMC->Gspos("CB0L",1, "SQM3",xposHorBox,yposHorBox,zposCbbBar,0, "ONLY");// box horizontal
     gMC->Gspos("CB1L",1, "SQM3", 0.0,0.0,zposCbb,0, "ONLY");
     gMC->Gspos("CB2L",1, "SQM3",xposVerBox,yposVerBox,zposCbbBar,0, "ONLY");// box vertical

     gMC->Gspos("CB0L",2, "SQM3",xposHorBox,yposHorBox,-zposCbbBar,0, "ONLY");// box horizontal
     gMC->Gspos("CB1L",2, "SQM3", 0.0,0.0,-(zposCbb+ zCbb),0, "ONLY");
     gMC->Gspos("CB2L",2, "SQM3",xposVerBox,yposVerBox,zposCbbBar,0, "ONLY");// box vertical

      // chamber 4
     gMC->Gspos("CB0L",3, "SQM4",xposHorBox,yposHorBox,zposCbbBar,0, "ONLY");
     gMC->Gspos("CB1L",3, "SQM4", 0.0,0.0,zposCbb,0, "ONLY");
     gMC->Gspos("CB2L",3, "SQM4",xposVerBox,yposVerBox,zposCbbBar,0, "ONLY");

     gMC->Gspos("CB0L",4, "SQM4",xposHorBox,yposHorBox,-zposCbbBar,0, "ONLY");
     gMC->Gspos("CB1L",4, "SQM4", 0.0,0.0,-(zposCbb+ zCbb),0, "ONLY");
     gMC->Gspos("CB2L",4, "SQM4",xposVerBox,yposVerBox,zposCbbBar,0, "ONLY");


//Positioning Copper sheet

      // chamber 3

     gMC->Gspos("CU0L",1, "SQM3",xposHorBox,yposHorBox,zposCuBar,0, "ONLY");
     gMC->Gspos("CU1L",1, "SQM3", 0.0,0.0,zposCu,0, "ONLY");
     gMC->Gspos("CU2L",1, "SQM3",xposVerBox,yposVerBox,zposCuBar,0, "ONLY");

     gMC->Gspos("CU0L",2, "SQM3",xposHorBox,yposHorBox,-zposCuBar,0, "ONLY");
     gMC->Gspos("CU1L",2, "SQM3", 0.0,0.0,-(zposCu+ zCu),0, "ONLY");
     gMC->Gspos("CU2L",2, "SQM3",xposVerBox,yposVerBox,zposCuBar,0, "ONLY");
      // chamber 4

     gMC->Gspos("CU0L",3, "SQM4",xposHorBox,yposHorBox,zposCuBar,0, "ONLY");
     gMC->Gspos("CU1L",3, "SQM4", 0.0,0.0,zposCu,0, "ONLY");
     gMC->Gspos("CU2L",3, "SQM4",xposVerBox,yposVerBox,zposCuBar,0, "ONLY");

     gMC->Gspos("CU0L",4, "SQM4",xposHorBox,yposHorBox,-zposCuBar,0, "ONLY");
     gMC->Gspos("CU1L",4, "SQM4", 0.0,0.0,-(zposCu+ zCu),0, "ONLY");
     gMC->Gspos("CU2L",4, "SQM4",xposVerBox,yposVerBox,zposCuBar,0, "ONLY");

//Positioning Roha cell 

      // chamber 3

     gMC->Gspos("RH0L",1, "SQM3",xposHorBox,yposHorBox,zposRohaBar,0, "ONLY");
     gMC->Gspos("RH1L",1, "SQM3", 0.0,0.0,zposRoha,0, "ONLY");
     gMC->Gspos("RH2L",1, "SQM3",xposVerBox,yposVerBox,zposRohaBar,0, "ONLY");

     gMC->Gspos("RH0L",2, "SQM3",xposHorBox,yposHorBox,-zposRohaBar,0, "ONLY");
     gMC->Gspos("RH1L",2, "SQM3", 0.0,0.0,-(zposRoha+ zRoha),0, "ONLY");
     gMC->Gspos("RH2L",2, "SQM3",xposVerBox,yposVerBox,zposRohaBar,0, "ONLY");

      // chamber 4

     gMC->Gspos("RH0L",3, "SQM4",xposHorBox,yposHorBox,zposRohaBar,0, "ONLY");
     gMC->Gspos("RH1L",3, "SQM4", 0.0,0.0,zposRoha,0, "ONLY");
     gMC->Gspos("RH2L",3, "SQM4",xposVerBox,yposVerBox,zposRohaBar,0, "ONLY");

     gMC->Gspos("RH0L",4, "SQM4",xposHorBox,yposHorBox,-zposRohaBar,0, "ONLY");
     gMC->Gspos("RH1L",4, "SQM4", 0.0,0.0,-(zposRoha+ zRoha),0, "ONLY");
     gMC->Gspos("RH2L",4, "SQM4",xposVerBox,yposVerBox,zposRohaBar,0, "ONLY");

//Positioning of Mech. exit board

      // chamber 3

     gMC->Gspos("MB0L",1, "SQM3",xposHorBox,yposHorBox,zposMebBar,0, "ONLY");
     gMC->Gspos("MB1L",1, "SQM3", 0.0,0.0,zposMeb,0, "ONLY");
     gMC->Gspos("MB2L",1, "SQM3",xposVerBox,yposVerBox,zposMebBar,0, "ONLY");

     gMC->Gspos("MB0L",2, "SQM3",xposHorBox,yposHorBox,-zposMebBar,0, "ONLY");
     gMC->Gspos("MB1L",2, "SQM3", 0.0,0.0,-(zposMeb+ zMeb),0, "ONLY");
     gMC->Gspos("MB2L",2, "SQM3",xposVerBox,yposVerBox,zposMebBar,0, "ONLY");
      // chamber 4

     gMC->Gspos("MB0L",3, "SQM4",xposHorBox,yposHorBox,zposMebBar,0, "ONLY");
     gMC->Gspos("MB1L",3, "SQM4", 0.0,0.0,zposMeb,0, "ONLY");
     gMC->Gspos("MB2L",3, "SQM4",xposVerBox,yposVerBox,zposMebBar,0, "ONLY");

     gMC->Gspos("MB0L",4, "SQM4",xposHorBox,yposHorBox,-zposMebBar,0, "ONLY");
     gMC->Gspos("MB1L",4, "SQM4", 0.0,0.0,-(zposMeb+ zMeb),0, "ONLY");
     gMC->Gspos("MB2L",4, "SQM4",xposVerBox,yposVerBox,zposMebBar,0, "ONLY");

//Positioning of Electronic exit board


      // chamber 3

     gMC->Gspos("EB0L",1, "SQM3",xposHorBox,yposHorBox,zposEebBar,0, "ONLY");
     gMC->Gspos("EB1L",1, "SQM3", 0.0,0.0,zposEeb,0, "ONLY");
     gMC->Gspos("EB2L",1, "SQM3",xposVerBox,yposVerBox,zposEebBar,0, "ONLY");

     gMC->Gspos("EB0L",2, "SQM3",xposHorBox,yposHorBox,-zposEebBar,0, "ONLY");
     gMC->Gspos("EB1L",2, "SQM3", 0.0,0.0,-(zposEeb+ zEeb),0, "ONLY");
     gMC->Gspos("EB2L",2, "SQM3",xposVerBox,yposVerBox,zposEebBar,0, "ONLY");
      // chamber 4

     gMC->Gspos("EB0L",3, "SQM4",xposHorBox,yposHorBox,zposEebBar,0, "ONLY");
     gMC->Gspos("EB1L",3, "SQM4", 0.0,0.0,zposEeb,0, "ONLY");
     gMC->Gspos("EB2L",3, "SQM4",xposVerBox,yposVerBox,zposEebBar,0, "ONLY");

     gMC->Gspos("EB0L",4, "SQM4",xposHorBox,yposHorBox,-zposEebBar,0, "ONLY");
     gMC->Gspos("EB1L",4, "SQM4", 0.0,0.0,-(zposEeb+ zEeb),0, "ONLY");
     gMC->Gspos("EB2L",4, "SQM4",xposVerBox,yposVerBox,zposEebBar,0, "ONLY");

 
//----------------------------------------------------------------------
//                          Frames
//----------------------------------------------------------------------
//Frame-1
     Float_t frame1[3] ;                
     frame1[0] = 101.0/2.;             //100.6 = 94.5 + 2.5 + 3.6 
     frame1[1] = 2.5/2.;               
     frame1[2] = 5.0/2.;

     gMC->Gsvolu("FRM1", "BOX", idPGF30, frame1, 3); //Frame - 1 // fill with pkk GF30
     
     Float_t arib1[3];
     arib1[0] = frame1[0];
     arib1[1] = 0.9/2.;
     arib1[2] =(frame1[2]-0.95)/2.0;
     
     gMC->Gsvolu("FRA1", "BOX", idAir, arib1, 3); // fill with air
     
     Float_t xposarib1 = 0;
     Float_t yposarib1 = -frame1[1] + arib1[1];
     Float_t zposarib1 = frame1[2] - arib1[2];
     
     gMC->Gspos("FRA1",1, "FRM1", xposarib1, yposarib1, zposarib1,0, "ONLY");  //replace pkk GF30 with air(b)
     gMC->Gspos("FRA1",2, "FRM1", xposarib1, yposarib1, -zposarib1,0, "ONLY"); //replace pkk GF30 with air(nb)
     
     Float_t rrib1[3];
     rrib1[0] = frame1[0];
     rrib1[1] = 0.6/2.;
     rrib1[2] =(frame1[2]-0.95)/2.0;
     
     gMC->Gsvolu("FRR1", "BOX", idRoha, rrib1, 3); // fill with rohacell
     
     Float_t xposrrib1 = 0.0;
     Float_t yposrrib1 = frame1[1] - rrib1[1];
     Float_t zposrrib1 = frame1[2] - rrib1[2];
     
     gMC->Gspos("FRR1",1, "FRM1", xposrrib1, yposrrib1, zposrrib1,0, "ONLY");//replace pkk GF30 with rohacell
     gMC->Gspos("FRR1",2, "FRM1", xposrrib1, yposrrib1, -zposrrib1,0, "ONLY");//replace pkk GF30 with rohacell

     
     Float_t xposFr1 = frame1[0] + 20.6;
     Float_t yposFr1 = -3.7 + frame1[1] ;
     Float_t zposFr1 = 0.0;
     
     gMC->Gspos("FRM1",1, "SQM3", xposFr1, yposFr1, zposFr1,0, "ONLY");// frame -1
     gMC->Gspos("FRM1",2, "SQM4", xposFr1, yposFr1, zposFr1,0, "ONLY");// frame -1
     

//......................................................................................
//Frame-2

      Float_t frame2[3]; 
      frame2[0] = 4.0/2.;
      frame2[1] = 1.2/2.;
      frame2[2] = 5.0/2;

      gMC->Gsvolu("FRM2", "BOX", idPGF30, frame2, 3); //Frame - 2

      Float_t rrib2[3];
      rrib2[0] = frame2[0]-1.0/2.0;
      rrib2[1] = frame2[1];
      rrib2[2] =(frame2[2]-0.95)/2.0;

      gMC->Gsvolu("FRR2", "BOX", idRoha, rrib2, 3);
      
      Float_t xposrrib2 = -1.0/2.0;
      Float_t yposrrib2 = 0.0;
      Float_t zposrrib2 = frame2[2] - rrib2[2];

      gMC->Gspos("FRR2",1, "FRM2", xposrrib2, yposrrib2, zposrrib2,0, "ONLY");//replace pkk GF30 with rohacell
      gMC->Gspos("FRR2",2, "FRM2", xposrrib2, yposrrib2, -zposrrib2,0, "ONLY");//replace pkk GF30 with roha



      Float_t xposFr2 = frame2[0] + 117.6;
      Float_t yposFr2 = -frame2[1];
      Float_t zposFr2 = 0.0;

      gMC->Gspos("FRM2",1, "SQM3", xposFr2, yposFr2, zposFr2,0, "ONLY");//global positing of frame in SQM3
      gMC->Gspos("FRM2",2, "SQM4", xposFr2, yposFr2, zposFr2,0, "ONLY");//global positing of frame in SQM4
      
  
//......................................................................................
//Frame-3

      Float_t pgparFr3[10];
      pgparFr3[0] = 0.;
      pgparFr3[1] = 90.;
      pgparFr3[2] = 5;
      pgparFr3[3] = 2.;
      pgparFr3[4] = 0.;
      pgparFr3[5] = 117.6;
      pgparFr3[6] = 121.6;
      pgparFr3[7] = pgparFr3[4] + 5.0;
      pgparFr3[8] = pgparFr3[5];
      pgparFr3[9] = pgparFr3[6];
      
      gMC->Gsvolu("FRM3", "PGON", idPGF30, pgparFr3, 10);
      
      Float_t pgparRrib3[10];
      pgparRrib3[0] = 0.;
      pgparRrib3[1] = 90.;
      pgparRrib3[2] = 5;
      pgparRrib3[3] = 2.;
      pgparRrib3[4] = 0.;
      pgparRrib3[5] = 117.6;
      pgparRrib3[6] = 120.6;
      pgparRrib3[7] = pgparRrib3[4] +1.55 ;
      pgparRrib3[8] = pgparRrib3[5];
      pgparRrib3[9] = pgparRrib3[6];
      
      gMC->Gsvolu("FRR3", "PGON", idRoha, pgparRrib3, 10);
      
      Float_t xposrrib3 = 0.0;
      Float_t yposrrib3 = 0.0;
      Float_t zposrrib3 = 0.0;

      gMC->Gspos("FRR3",1, "FRM3", xposrrib3, yposrrib3, zposrrib3,0, "ONLY");

      zposrrib3 = 3.45;
      
      gMC->Gspos("FRR3",2, "FRM3", xposrrib3, yposrrib3, zposrrib3,0, "ONLY");
      
      
      
      Float_t xposFr3 = 0.0;
      Float_t yposFr3 = 0.0;
      Float_t zposFr3 = -frame1[2];
      
      gMC->Gspos("FRM3",1, "SQM3", xposFr3, yposFr3, zposFr3,0, "ONLY");// frame -1
      gMC->Gspos("FRM3",2, "SQM4", xposFr3, yposFr3, zposFr3,0, "ONLY");// frame -1
      

//......................................................................................
//Frame-4
 
      Float_t frame4[3]; 
      frame4[0] = 1.0/2.;
      frame4[1] = 4.0/2.;
      frame4[2] = frame1[2];

      gMC->Gsvolu("FRM4", "BOX", idPGF30, frame4, 3); 

      Float_t rrib4[3];
      rrib4[0] = frame4[0];
      rrib4[1] = frame4[1]-1.0/2;
      rrib4[2] =(frame4[2]-0.95)/2.0;

      gMC->Gsvolu("FRR4", "BOX", idRoha, rrib4, 3);
      
      Float_t xposrrib4 = 0.0;
      Float_t yposrrib4 = -1.0/2;
      Float_t zposrrib4 = frame4[2] - rrib4[2];

      gMC->Gspos("FRR4",1, "FRM4", xposrrib4, yposrrib4, zposrrib4,0, "ONLY");
      gMC->Gspos("FRR4",2, "FRM4", xposrrib4, yposrrib4, -zposrrib4,0, "ONLY");



      Float_t xposFr4 = -frame4[0];
      Float_t yposFr4 = -frame4[1] + 117.6;
      Float_t zposFr4 = 0.0;

      gMC->Gspos("FRM4",1, "SQM3", xposFr4, yposFr4, zposFr4,0, "ONLY");
      gMC->Gspos("FRM4",2, "SQM4", xposFr4, yposFr4, zposFr4,0, "ONLY");
      

//......................................................................................
//Frame-5
      Float_t frame5[3] ;               
      frame5[0] = 2.7/2.;            
      frame5[1] = 101.0/2.;              
      frame5[2] = 5.0/2.;

      gMC->Gsvolu("FRM5", "BOX", idPGF30, frame5, 3); //Frame - 1
      
      Float_t arib5[3];
      arib5[0] = 0.9/2.0;
      arib5[1] = frame5[1];
      arib5[2] = (frame5[2]-0.95)/2.0;

      gMC->Gsvolu("FRA5", "BOX", idAir, arib5, 3);

      Float_t xposarib5 = -frame5[0] + arib5[0];
      Float_t yposarib5 = 0.0;
      Float_t zposarib5 = frame5[2] - arib5[2];

      gMC->Gspos("FRA5",1, "FRM5", xposarib5, yposarib5, zposarib5,0, "ONLY");
      gMC->Gspos("FRA5",2, "FRM5", xposarib5, yposarib5, -zposarib5,0, "ONLY");

      Float_t rrib5[3];
      rrib5[0] = 0.8/2.0;
      rrib5[1] = frame5[1];
      rrib5[2] = (frame5[2]-0.95)/2.0;

      gMC->Gsvolu("FRR5", "BOX", idRoha, rrib5, 3);

      Float_t xposrrib5 = frame5[0] - rrib5[0];
      Float_t yposrrib5 = 0.0;
      Float_t zposrrib5 = frame5[2] - rrib5[2];

      gMC->Gspos("FRR5",1, "FRM5", xposrrib5, yposrrib5, zposrrib5,0, "ONLY");
      gMC->Gspos("FRR5",2, "FRM5", xposrrib5, yposrrib5, -zposrrib5,0, "ONLY");

      Float_t xposFr5 = -3.7 + frame5[0];
      Float_t yposFr5 = frame5[1] + 20.6;
      Float_t zposFr5 = 0.0;

      gMC->Gspos("FRM5",1, "SQM3", xposFr5, yposFr5, zposFr5,0, "ONLY");
      gMC->Gspos("FRM5",2, "SQM4", xposFr5, yposFr5, zposFr5,0, "ONLY");

//......................................................................................
//Frame -6 
 
      Float_t frame6[3]; 
      frame6[0] = 1.0/2.;
      frame6[1] = 2.5/2.;
      frame6[2] = frame1[2];

      gMC->Gsvolu("FRM6", "BOX", idPGF30, frame6, 3); 

      Float_t rrib6[3];
      rrib6[0] = frame6[0];
      rrib6[1] = 1.5/2.;
      rrib6[2] =(frame2[2]-0.95)/2.0;

      gMC->Gsvolu("FRR6", "BOX", idRoha, rrib6, 3);
      
      Float_t xposrrib6 = 0.0;
      Float_t yposrrib6 = 1.0/2.0;
      Float_t zposrrib6 = frame6[2] - rrib6[2];

      gMC->Gspos("FRR6",1, "FRM6", xposrrib6, yposrrib6, zposrrib6,0, "ONLY");
      gMC->Gspos("FRR6",2, "FRM6", xposrrib6, yposrrib6, -zposrrib6,0, "ONLY");



      Float_t xposFr6 = -frame6[0];
      Float_t yposFr6 = frame6[1] + 20.6;
      Float_t zposFr6 = 0.0;

      gMC->Gspos("FRM6",1, "SQM3", xposFr6, yposFr6, zposFr6,0, "ONLY");
      gMC->Gspos("FRM6",2, "SQM4", xposFr6, yposFr6, zposFr6,0, "ONLY");
      

//......................................................................................
//Frame - 7 inner pgon

        Float_t pgparFr7[10];
        pgparFr7[0] = 0.;
        pgparFr7[1] = 90.;
        pgparFr7[2] = 5;
        pgparFr7[3] = 2.;
        pgparFr7[4] = 0.;
        pgparFr7[5] = 20.6;
        pgparFr7[6] = 23.1;
        pgparFr7[7] = pgparFr7[4] + 5.0;
        pgparFr7[8] = pgparFr7[5];
        pgparFr7[9] = pgparFr7[6];

        gMC->Gsvolu("FRM7", "PGON", idPGF30, pgparFr7, 10);

        Float_t pgparRrib7[10];
        pgparRrib7[0] = 0.;
        pgparRrib7[1] = 90.;
        pgparRrib7[2] = 5;
        pgparRrib7[3] = 2.;
        pgparRrib7[4] = 0.;
        pgparRrib7[5] = 21.6;
        pgparRrib7[6] = 23.1;
        pgparRrib7[7] = pgparRrib7[4] + 1.55;
        pgparRrib7[8] = pgparRrib7[5];
        pgparRrib7[9] = pgparRrib7[6];

        gMC->Gsvolu("FRR7", "PGON", idRoha, pgparRrib7, 10);
      
	Float_t xposrrib7 = 0.0;
	Float_t yposrrib7 = 0.0;
	Float_t zposrrib7 = 0.0;

	gMC->Gspos("FRR7",1, "FRM7", xposrrib7, yposrrib7, zposrrib7,0, "ONLY");

	zposrrib7 = 3.45;

	gMC->Gspos("FRR7",2, "FRM7", xposrrib7, yposrrib7, zposrrib7,0, "ONLY");

	Float_t xposFr7 = 0.0;
	Float_t yposFr7 = 0.0;
	Float_t zposFr7 = -frame1[2];

	gMC->Gspos("FRM7",1, "SQM3", xposFr7, yposFr7, zposFr7,0, "ONLY");
	gMC->Gspos("FRM7",2, "SQM4", xposFr7, yposFr7, zposFr7,0, "ONLY");
      

//......................................................................................
//Frame - 8 


      Float_t frame8[3] ;
      frame8[0] = 2.5/2.0;
      frame8[1] = 1.2/2.0;
      frame8[2] = frame1[2];

      gMC->Gsvolu("FRM8", "BOX", idPGF30, frame8, 3); //Frame - 2

      Float_t rrib8[3];
      rrib8[0] = frame8[0]-1.0/2;
      rrib8[1] = frame8[1];
      rrib8[2] =(frame8[2]-0.95)/2.0;

      gMC->Gsvolu("FRR8", "BOX", idRoha, rrib8, 3);
      
      Float_t xposrrib8 = -1.0/2;
      Float_t yposrrib8 = 0.0;
      Float_t zposrrib8 = frame8[2] - rrib8[2];

      gMC->Gspos("FRR8",1, "FRM8", xposrrib8, yposrrib8, zposrrib8,0, "ONLY");
      gMC->Gspos("FRR8",2, "FRM8", xposrrib8, yposrrib8, -zposrrib8,0, "ONLY");



      Float_t xposFr8 = frame8[0] + 20.6;
      Float_t yposFr8 = -frame8[1];
      Float_t zposFr8 = 0.0;

      gMC->Gspos("FRM8",1, "SQM3", xposFr8, yposFr8, zposFr8,0, "ONLY");
      gMC->Gspos("FRM8",2, "SQM4", xposFr8, yposFr8, zposFr8,0, "ONLY");
      
 


//^^^^^^^^^^^^^^^^^^^^^^^^^ Sensitive volumes ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       Float_t zsenv = 0.5; // distance between two cathode plane

 //Segment-0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       bparH[0] = 94.5/2.;
       bparH[1] = 1.2/2.;
       bparH[2] = zsenv/2.;
       gMC->Gsvolu("C3G0", "BOX", idGas, bparH, 3);
       gMC->Gsvolu("C4G0", "BOX", idGas, bparH, 3);

 //Segment-1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pgpar[0] = 0.;
        pgpar[1] = 90.;
        pgpar[2] = 5;
        pgpar[3] = 2.;
        pgpar[4] = -zsenv/2.;
        pgpar[5] = 23.1;
        pgpar[6] = 117.6;
        pgpar[7] = zsenv/2.;
        pgpar[8] = pgpar[5];
        pgpar[9] = pgpar[6];
        gMC->Gsvolu("C3G1", "PGON", idGas, pgpar, 10);
        gMC->Gsvolu("C4G1", "PGON", idGas, pgpar, 10);


//Segment-2 - vertical box ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        bparV[0] = 1.0/2.;
	bparV[1] = 95.5/2.;
	bparV[2] = zsenv/2.;
        gMC->Gsvolu("C3G2", "BOX", idGas, bparV, 3);
        gMC->Gsvolu("C4G2", "BOX", idGas, bparV, 3);

//...........................................................................................

     xposHorBox =  bparH[0] + 23.1;
     yposHorBox = -bparH[1];

     xposVerBox = -bparV[0];
     yposVerBox = bparV[1] + 23.1;

     gMC->Gspos("C3G0",1, "SQM3", xposHorBox,yposHorBox,0.,0, "ONLY");
     gMC->Gspos("C3G1",1, "SQM3", 0.,0.,0.,0, "ONLY");
     gMC->Gspos("C3G2",1, "SQM3", xposVerBox,yposVerBox,0.,0, "ONLY");


     gMC->Gspos("C4G0",1, "SQM4", xposHorBox,yposHorBox,0.,0, "ONLY");
     gMC->Gspos("C4G1",1, "SQM4", 0.,0.,0.,0, "ONLY");
     gMC->Gspos("C4G2",1, "SQM4", xposVerBox,yposVerBox,0.,0, "ONLY");


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

 //Float_t halfChamber = zCbb + zCu + zRoha + zMeb + zEeb + zsenv/2;
//   cout<<  "\n half_chamber \t" << half_chamber << endl;

 Float_t halfChamber = tpar1[2];// tpar1[2] = 6.8/2;

// ------------------------------St2 Chamber3------------------------------------------------

 //    GetEnvelopes(2)->AddEnvelope("S3M0", 300, true,TGeoTranslation(0.,0.,0.));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId1, 1, TGeoTranslation( 0., 0., - halfChamber));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId2, 2, TGeoTranslation( 0., 0., + halfChamber),
                                 TGeoRotation("Qrot3",90.,180.,90.,90.,180.,0.));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId3, 3, TGeoTranslation( 0., 0., - halfChamber),
                                  TGeoRotation("Qrot3",90.,180.,90.,270.,0.,0.));
    GetEnvelopes(2)->AddEnvelope("SQM3", 300+detElemId4, 4, TGeoTranslation( 0., 0., + halfChamber),
                                  TGeoRotation("Qrot4",90.,0.,90.,-90.,180.,0.));

//--------------------------------St2 Chamber4-------------------------------------------------

    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId1, 1, TGeoTranslation( 0., 0., - halfChamber));
    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId2, 2, TGeoTranslation( 0., 0., + halfChamber),
                                 TGeoRotation("Qrot2",90.,180.,90.,90.,180.,0.));
    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId3, 3, TGeoTranslation( 0., 0., - halfChamber),
                                  TGeoRotation("Qrot3",90.,180.,90.,270.,0.,0.));
    GetEnvelopes(3)->AddEnvelope("SQM4", 400+detElemId4, 4, TGeoTranslation( 0., 0., + halfChamber),
                                  TGeoRotation("Qrot4",90.,0.,90.,-90.,180.,0.));

//**********************************************************************************************

}
//______________________________________________________________________________
void AliMUONSt2GeometryBuilderV2::SetTransformations()
{
// Defines the transformations for the station2 chambers.
// ---

  if (gAlice->GetModule("SHIL")) {
    SetMotherVolume(2, "YOUT1");
    SetMotherVolume(3, "YOUT1");
  }  

  // Define chamber volumes as virtual
  SetVolume(2, "SC03", true);
  SetVolume(3, "SC04", true);

  Double_t zpos1 = - AliMUONConstants::DefaultChamberZ(2); 
  SetTranslation(2, TGeoTranslation(0., 0., zpos1));

  Double_t zpos2 = - AliMUONConstants::DefaultChamberZ(3); 
  SetTranslation(3, TGeoTranslation(0., 0., zpos2));
}

//______________________________________________________________________________
void AliMUONSt2GeometryBuilderV2::SetSensitiveVolumes()
{
// Defines the sensitive volumes for station2 chambers.
// ---

  GetGeometry(2)->SetSensitiveVolume("C3G0");
  GetGeometry(2)->SetSensitiveVolume("C3G1");
  GetGeometry(2)->SetSensitiveVolume("C3G2");

  GetGeometry(3)->SetSensitiveVolume("C4G0");
  GetGeometry(3)->SetSensitiveVolume("C4G1");
  GetGeometry(3)->SetSensitiveVolume("C4G2");
}

