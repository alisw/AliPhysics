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
Revision 1.1  1999/10/15 15:35:20  fca
New version for frame1099 with and without holes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight: design of C.Williams                                //
//  This class contains the functions for version 1 of the Time Of Flight    //
//  detector.                                                                //
//
//  VERSION WITH 5 MODULES AND FLAT PLATES
//  
//   WITH HOLES FOR PHOS AND HMPID 
//   INSIDE THE FULL COVERAGE SPACE FRAME
//
//
//   Authors:
//
//   Alessio Seganti
//   Domenico Vicinanza
//
//   University of Salerno - Italy
//
//Begin_Html
/*
<img src="picts/AliTOFv6Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFv6.h"
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliTOFv6)
 
//_____________________________________________________________________________
AliTOFv6::AliTOFv6()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv6::AliTOFv6(const char *name, const char *title)
       : AliTOF(name,title)
{
  //
  // Standard constructor
  //
}
 
//_____________________________________________________________________________
void AliTOFv6::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv6.gif">
  */
  //End_Html
  //
  // Creates common geometry
  //
  AliTOF::CreateGeometry();
}
 
//_____________________________________________________________________________
void AliTOFv6::TOFpc(Float_t xtof, Float_t ytof, Float_t zlen1,
		     Float_t zlen2, Float_t zlen3, Float_t ztof0)
{
  //
  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xFLT, yFLT, zFLT - sizes of TOF modules (large)
  
  Float_t ycoor;
  Float_t par[10];
  Int_t idrotm[100];

  Int_t *idtmed = fIdtmed->GetArray()-499;
  

  par[0] =  xtof / 2.;
  par[1] =  ytof / 2.;
  par[2] = zlen1 / 2.;
  gMC->Gsvolu("FTO1", "BOX ", idtmed[506], par, 3);
  par[2] = zlen2 / 2.;
  gMC->Gsvolu("FTO2", "BOX ", idtmed[506], par, 3);
  par[2] = zlen3 / 2.;
  gMC->Gsvolu("FTO3", "BOX ", idtmed[506], par, 3);


// Positioning of modules

   Float_t zcor1 = ztof0 - zlen1/2;
   Float_t zcor2 = ztof0 - zlen1 - zlen2/2.;
   Float_t zcor3 = 0.;

   AliMatrix(idrotm[0], 90., 0., 0., 0., 90, -90.);
   AliMatrix(idrotm[1], 90., 180., 0., 0., 90, 90.);
   gMC->Gspos("FTO1", 1, "BTO1", 0,  zcor1, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTO1", 2, "BTO1", 0, -zcor1, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTO1", 1, "BTO2", 0,  zcor1, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTO1", 2, "BTO2", 0, -zcor1, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTO1", 1, "BTO3", 0,  zcor1, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTO1", 2, "BTO3", 0, -zcor1, 0, idrotm[1], "ONLY");

   gMC->Gspos("FTO2", 1, "BTO1", 0,  zcor2, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTO2", 2, "BTO1", 0, -zcor2, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTO2", 1, "BTO2", 0,  zcor2, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTO2", 2, "BTO2", 0, -zcor2, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTO2", 1, "BTO3", 0,  zcor2, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTO2", 2, "BTO3", 0, -zcor2, 0, idrotm[1], "ONLY");

   gMC->Gspos("FTO3", 0, "BTO1", 0, zcor3,  0, idrotm[0], "ONLY");
   gMC->Gspos("FTO3", 0, "BTO2", 0, zcor3,  0, idrotm[0], "ONLY");
   gMC->Gspos("FTO3", 0, "BTO3", 0, zcor3,  0, idrotm[0], "ONLY");

// Subtraction the distance to TOF module boundaries 

  Float_t db = 7.;
  Float_t xFLT, yFLT, zFLT1, zFLT2, zFLT3;

  xFLT = xtof -(.5 +.5)*2;
  yFLT = ytof;
  zFLT1 = zlen1 - db;
  zFLT2 = zlen2 - db;
  zFLT3 = zlen3 - db;

  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xFLT, yFLT, zFLT - sizes of TOF modules (large)
  
  Float_t yFREON, xp, yp, zp;
  
// fron gaps in MRPC chamber 
  yFREON = .11; //cm

// Sizes of MRPC pads

  xp = 3.0; 
  yp = 12.3*0.05; // 5% X0 of glass 
  zp = 3.0;

//  Subtraction of dead boundaries in X=2 cm and Z=7/2 cm 

cout <<"************************* TOF geometry **************************"<<endl;

  Int_t nz1, nz2, nz3, nx; //- numbers of pixels
  nx = Int_t (xFLT/xp);

  printf("Number of pixel along x axis = %i",nx);

  par[0] = xFLT/2;
  par[1] = yFLT/2;
  par[2] = (zFLT1 / 2.);
  nz1 = Int_t (par[2]*2/zp);
  gMC->Gsvolu("FLT1", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT1", 0, "FTO1", 0., 0., 0., 0, "ONLY");
  printf("Number of pixel along z axis (module 1) = %i",nz1);

  par[2] = (zFLT2 / 2.);
  nz2 = Int_t (par[2]*2/zp);
  gMC->Gsvolu("FLT2", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT2", 0, "FTO2", 0., 0., 0., 0, "ONLY");
  printf("Number of pixel along z axis (module 2) = %i",nz2);

  par[2] = (zFLT3 / 2.); 
  nz3 = Int_t (par[2]*2/zp);
  gMC->Gsvolu("FLT3", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT3", 0, "FTO3", 0., 0., 0., 0, "ONLY");
  printf("Number of pixel along z axis (module 3) = %i",nz3);

////////// Layers before detector ////////////////////

// Alluminium layer in front 1.0 mm thick at the beginning
  par[0] = -1;
  par[1] = 0.1;
  par[2] = -1;
  ycoor = -yFLT/2 + par[1];
  gMC->Gsvolu("FMY1", "BOX ", idtmed[508], par, 3); // Alluminium
  gMC->Gspos("FMY1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FMY2", "BOX ", idtmed[508], par, 3); // Alluminium
  gMC->Gspos("FMY2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FMY3", "BOX ", idtmed[508], par, 3); // Alluminium 
  gMC->Gspos("FMY3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

// Honeycomb layer (1cm of special polyethilene)
  ycoor = ycoor + par[1];
  par[0] = -1;
  par[1] = 0.5;
  par[2] = -1;
  ycoor = ycoor + par[1];
  gMC->Gsvolu("FPL1", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPL1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPL2", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPL2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPL3", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPL3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

///////////////// Detector itself //////////////////////

  const Float_t SpaceBefore = 2.;

  par[0] = -1;
  par[1] = yp/2; // 5 %X0 thick of glass  
  par[2] = -1;
  ycoor = -yFLT/2 + SpaceBefore;
  gMC->Gsvolu("FLD1", "BOX ", idtmed[514], par, 3); // Glass
  gMC->Gspos("FLD1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FLD2", "BOX ", idtmed[514], par, 3); // Glass
  gMC->Gspos("FLD2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FLD3", "BOX ", idtmed[514], par, 3); // Glass
  gMC->Gspos("FLD3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

  gMC->Gsdvn("FLZ1", "FLD1", nz1, 3); //pixel size xp=zp=3
  gMC->Gsdvn("FLZ2", "FLD2", nz2, 3); 
  gMC->Gsdvn("FLZ3", "FLD3", nz3, 3); 
  gMC->Gsdvn("FLX1", "FLZ1", nx, 1);
  gMC->Gsdvn("FLX2", "FLZ2", nx, 1); 
  gMC->Gsdvn("FLX3", "FLZ3", nx, 1); 

// MRPC pixel itself 
  par[0] = -1;//xp/2;
  par[1] = -1;//yp/2; // 5 %X0 thick of glass  
  par[2] = -1;//zp/2;
  gMC->Gsvolu("FPA0", "BOX ", idtmed[514], par, 3);// Glass
  gMC->Gspos("FPA0", 1, "FLX1", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("FPA0", 2, "FLX2", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("FPA0", 3, "FLX3", 0., 0., 0., 0, "ONLY");

// Freon gas sencitive vol.ume
  par[0] = -1;
  par[1] = yFREON/2;
  par[2] = -1;
  gMC->Gsvolu("FPAD", "BOX ", idtmed[513], par, 3);// Freon 
  gMC->Gspos("FPAD", 0, "FPA0", 0., 0., 0., 0, "ONLY");

////////// Layers after detector ////////////////////

  const Float_t SpaceAfter=6.;

// Honeycomb layer after (3cm)
  par[0] = -1;
  par[1] = 0.6;
  par[2] = -1;
  ycoor = -yFLT/2 + SpaceAfter - par[1];
  gMC->Gsvolu("FPE1", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPE2", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPE3", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

// Electronics (Cu) after
  par[0] = -1;
  par[1] = 1.43*0.05 / 2.; // 5% of X0
  par[2] = -1;
  ycoor = -yFLT/2 + SpaceAfter +par[1];
  gMC->Gsvolu("FEC1", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEC2", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEC3", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

// Cooling water after
  ycoor = ycoor+par[1];
  par[0] = -1;
  par[1] = 36.1*0.02 / 2.; // 2% of X0
  par[2] = -1;
  ycoor = ycoor+par[1];
  gMC->Gsvolu("FWA1", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos("FWA1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FWA2", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos("FWA2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FWA3", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos("FWA3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

//back plate honycomb (2cm)
  par[0] = -1;
  par[1] = 2 / 2.;
  par[2] = -1;
  ycoor = yFLT/2 - par[1];
  gMC->Gsvolu("FEG1", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FEG1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEG2", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FEG2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEG3", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FEG3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");


}

//_____________________________________________________________________________
void AliTOFv6::DrawModule()
{
  //
  // Draw a shaded view of the Time Of Flight version 1
  //
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("ALIC","SEEN",0);
  gMC->Gsatt("FBAR","SEEN",1);
  gMC->Gsatt("FTO1","SEEN",1);
  gMC->Gsatt("FTO2","SEEN",1);
  gMC->Gsatt("FTO3","SEEN",1);
  gMC->Gsatt("FBT1","SEEN",1);
  gMC->Gsatt("FBT2","SEEN",1);
  gMC->Gsatt("FBT3","SEEN",1);
  gMC->Gsatt("FDT1","SEEN",1);
  gMC->Gsatt("FDT2","SEEN",1);
  gMC->Gsatt("FDT3","SEEN",1);
  gMC->Gsatt("FLT1","SEEN",1);
  gMC->Gsatt("FLT2","SEEN",1);
  gMC->Gsatt("FLT3","SEEN",1);
  gMC->Gsatt("FPL1","SEEN",1);
  gMC->Gsatt("FPL2","SEEN",1);
  gMC->Gsatt("FPL3","SEEN",1);
  gMC->Gsatt("FLD1","SEEN",1);
  gMC->Gsatt("FLD2","SEEN",1);
  gMC->Gsatt("FLD3","SEEN",1);
  gMC->Gsatt("FLZ1","SEEN",1);
  gMC->Gsatt("FLZ2","SEEN",1);
  gMC->Gsatt("FLZ3","SEEN",1);
  gMC->Gsatt("FLX1","SEEN",1);
  gMC->Gsatt("FLX2","SEEN",1);
  gMC->Gsatt("FLX3","SEEN",1);
  gMC->Gsatt("FPA0","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .02, .02);
  gMC->Gdhead(1111, "Time Of Flight");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTOFv6::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //
  AliTOF::CreateMaterials();
}
 
//_____________________________________________________________________________
void AliTOFv6::Init()
{
  //
  // Initialise the detector after the geometry has been defined
  //
  AliTOF::Init();
  fIdFTO2=gMC->VolId("FTO2");
  fIdFTO3=gMC->VolId("FTO3");
  fIdFLT1=gMC->VolId("FLT1");
  fIdFLT2=gMC->VolId("FLT2");
  fIdFLT3=gMC->VolId("FLT3");
}
 
//_____________________________________________________________________________
void AliTOFv6::StepManager()
{
  //
  // Procedure called at each step in the Time Of Flight
  //
  TLorentzVector mom, pos;
  Float_t hits[8];
  Int_t vol[3];
  Int_t copy, id, i;
  Int_t *idtmed = fIdtmed->GetArray()-499;
  if(gMC->GetMedium()==idtmed[514-1] && 
     gMC->IsTrackEntering() && gMC->TrackCharge()
     && gMC->CurrentVolID(copy)==fIdSens) {
    TClonesArray &lhits = *fHits;
    //
    // Record only charged tracks at entrance
    gMC->CurrentVolOffID(1,copy);
    vol[2]=copy;
    gMC->CurrentVolOffID(3,copy);
    vol[1]=copy;
    id=gMC->CurrentVolOffID(8,copy);
    vol[0]=copy;
    if(id==fIdFTO3) {
      vol[0]+=22;
      id=gMC->CurrentVolOffID(5,copy);
      if(id==fIdFLT3) vol[1]+=6;
    } else if (id==fIdFTO2) {
      vol[0]+=20;
      id=gMC->CurrentVolOffID(5,copy);
      if(id==fIdFLT2) vol[1]+=8;
    } else {
      id=gMC->CurrentVolOffID(5,copy);
      if(id==fIdFLT1) vol[1]+=14;
    }
    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);
    //
    Double_t ptot=mom.Rho();
    Double_t norm=1/ptot;
    for(i=0;i<3;++i) {
      hits[i]=pos[i];
      hits[i+3]=mom[i]*norm;
    }
    hits[6]=ptot;
    hits[7]=pos[3];
    new(lhits[fNhits++]) AliTOFhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }
}

