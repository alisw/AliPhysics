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
Revision 1.12  1999/10/22 08:04:14  fca
Correct improper use of negative parameters

Revision 1.11  1999/10/16 19:30:05  fca
Corrected Rotation Matrix and CVS log

Revision 1.10  1999/10/15 15:35:20  fca
New version for frame1099 with and without holes

Revision 1.9  1999/09/29 09:24:33  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight: design of C.Williams                FCA                  //
//  This class contains the functions for version 1 of the Time Of Flight    //
//  detector.                                                                //
//
//  VERSION WITH 5 MODULES AND TILTED STRIPS 
//  
//   WITH HOLES FOR PHOS AND HMPID 
//   INSIDE A FULL COVERAGE SPACE FRAME
//
//
//   Authors: 
//
//   Alessio Seganti
//   Domenico Vicinanza
//
//   University of Salerno - Italy
//
//
//
//Begin_Html
/*
<img src="picts/AliTOFv2Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include "AliTOFv2.h"
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliTOFv2)
 
//_____________________________________________________________________________
AliTOFv2::AliTOFv2()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv2::AliTOFv2(const char *name, const char *title)
       : AliTOF(name,title)
{
  //
  // Standard constructor
  //
}
 
//_____________________________________________________________________________
void AliTOFv2::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv2.gif">
  */
  //End_Html
  //
  // Creates common geometry
  //
  AliTOF::CreateGeometry();
}
 
//_____________________________________________________________________________
void AliTOFv2::TOFpc(Float_t xtof, Float_t ytof, Float_t zlen1,
		     Float_t zlen2, Float_t zlen3, Float_t ztof0)
{
  //
  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xFLT, yFLT, zFLT - sizes of TOF modules (large)
  
  Int_t idrotm[100];
  Int_t nrot = 0;
  Float_t  ycoor, zcoor;
  Float_t par[10];
  
  Int_t *idtmed = fIdtmed->GetArray()-499;


  par[0] =  xtof / 2.;
  par[1] =  ytof / 2.;
  par[2] = zlen1 / 2.;
  gMC->Gsvolu("FTO1", "BOX ", idtmed[506], par, 3);
  par[2] = zlen2 / 2.;
  gMC->Gsvolu("FTO2", "BOX ", idtmed[506], par, 3);
  par[2] = zlen3 / 2.;
  gMC->Gsvolu("FTO3", "BOX ", idtmed[506], par, 3);


// Position of modules
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

   gMC->Gspos("FTO3", 0, "BTO1", 0, zcor3,  0, idrotm[0], "ONLY");

// Subtraction the distance to TOF module boundaries 

  Float_t db = 7.;
  Float_t xFLT, yFLT, zFLT1, zFLT2, zFLT3; 

  xFLT = xtof -(.5 +.5)*2;
  yFLT = ytof;
  zFLT1 = zlen1 - db;
  zFLT2 = zlen2 - db;
  zFLT3 = zlen3 - db;

  
// Sizes of MRPC pads

  Float_t yPad = 0.505; 
  
// Large not sensitive volumes with CO2 
  par[0] = xFLT/2;
  par[1] = yFLT/2;

  cout <<"************************* TOF geometry **************************"<<endl;

  par[2] = (zFLT1 / 2.);
  gMC->Gsvolu("FLT1", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT1", 0, "FTO1", 0., 0., 0., 0, "ONLY");

  par[2] = (zFLT2 / 2.);
  gMC->Gsvolu("FLT2", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT2", 0, "FTO2", 0., 0., 0., 0, "ONLY");

  par[2] = (zFLT3 / 2.); 
  gMC->Gsvolu("FLT3", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT3", 0, "FTO3", 0., 0., 0., 0, "ONLY");

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

  const Float_t StripWidth = 7.81;//cm
  const Float_t DeadBound = 1.;//cm non-sensitive between the pad edge and the boundary of the strip
  const Int_t nx = 40; // number of pads along x
  const Int_t nz = 2;  // number of pads along z
  const Float_t Gap=4.; //cm  distance between the strip axis
  const Float_t Space = 5.5; //cm distance from the front plate of the box

  Float_t zSenStrip;
  zSenStrip = StripWidth-2*DeadBound;//cm

  par[0] = xFLT/2;
  par[1] = yPad/2; 
  par[2] = StripWidth/2.;
  
  // Glass Layer of detector
  gMC->Gsvolu("FSTR","BOX",idtmed[514],par,3);

  // Freon for non-sesitive boundaries
  par[0] = xFLT/2;
  par[1] = 0.110/2;
  par[2] = -1;
  gMC->Gsvolu("FNSF","BOX",idtmed[512],par,3);
  gMC->Gspos("FNSF",0,"FSTR",0.,0.,0.,0,"ONLY");
  // Mylar for non-sesitive boundaries
  par[1] = 0.025;
  gMC->Gsvolu("FMYI","BOX",idtmed[510],par,3); 
  gMC->Gspos("FMYI",0,"FNSF",0.,0.,0.,0,"ONLY");

  // Mylar for outer layers
  par[1] = 0.035/2;
  ycoor = -yPad/2.+par[1];
  gMC->Gsvolu("FMYX","BOX",idtmed[510],par,3);
  gMC->Gspos("FMYX",1,"FSTR",0.,ycoor,0.,0,"ONLY");
  gMC->Gspos("FMYX",2,"FSTR",0.,-ycoor,0.,0,"ONLY");
  ycoor += par[1];
 
  // Graphyte layers
  par[1] = 0.003/2;
  ycoor += par[1];
  gMC->Gsvolu("FGRL","BOX",idtmed[502],par,3);
  gMC->Gspos("FGRL",1,"FSTR",0.,ycoor,0.,0,"ONLY");
  gMC->Gspos("FGRL",2,"FSTR",0.,-ycoor,0.,0,"ONLY");

  // Freon sensitive layer
  par[0] = -1;
  par[1] = 0.110/2.;
  par[2] = zSenStrip/2.;
  gMC->Gsvolu("FCFC","BOX",idtmed[513],par,3);
  gMC->Gspos("FCFC",0,"FNSF",0.,0.,0.,0,"ONLY");
  
  // Pad definition x & z
  gMC->Gsdvn("FLZ","FCFC", nz, 3); 
  gMC->Gsdvn("FLX","FLZ" , nx, 1); 

  // MRPC pixel itself 
  par[0] = -1;
  par[1] = -1; 
  par[2] = -1;
  gMC->Gsvolu("FPAD", "BOX ", idtmed[513], par, 3);
  gMC->Gspos("FPAD", 0, "FLX", 0., 0., 0., 0, "ONLY");


////  Positioning the Strips  (FSTR) in the FLT volumes  /////

 
  // 3 (Central) Plate 
  Float_t t = zFLT1+zFLT2+zFLT3/2.+7.*2.5;//Half Width of Barrel
  Float_t zpos = 0;
  Float_t ang;
  Float_t Offset;  
  Float_t last;
  nrot = 0;
  Int_t i=1,j=1;
  zcoor=0;
  Int_t UpDown=-1; // UpDown=-1 -> Upper strip, UpDown=+1 -> Lower strip
 
  do{
     ang = atan(zcoor/t);
     ang = ang*kRaddeg;
     AliMatrix (idrotm[nrot]  ,90.,  0.,90.-ang,90.,-ang,90.);
     AliMatrix (idrotm[nrot+1],90.,180.,90.+ang,90.,ang,90.);
     ycoor = -29./2.+ Space; //2 cm over front plate
     ycoor += (1-(UpDown+1)/2)*Gap;
     gMC->Gspos("FSTR",j,"FLT3",0.,ycoor,zcoor,idrotm[nrot],"ONLY");
     gMC->Gspos("FSTR",j+1,"FLT3",0.,ycoor,-zcoor,idrotm[nrot+1],"ONLY");
     ang  = ang/kRaddeg;
     
     zcoor=zcoor-(zSenStrip/2)/TMath::Cos(ang)+UpDown*Gap*TMath::Tan(ang)-(zSenStrip/2)/TMath::Cos(ang);
     UpDown*= -1; // Alternate strips 
     i++;
     j+=2;
  } while (zcoor-(StripWidth/2)*TMath::Cos(ang)>-t+zFLT1+zFLT2+7*2.5);
  
  ycoor = -29./2.+ Space; //2 cm over front plate

  // Plate  2
  zpos = -zFLT3/2-7;
  ang = atan(zpos/sqrt(2*t*t-zpos*zpos));
  Offset = StripWidth*TMath::Cos(ang)/2;
  zpos -= Offset;
  nrot = 0;
  i=1;
  // UpDown has not to be reinitialized, so that the arrangement of the strips can continue coherently

  do {
     ang = atan(zpos/sqrt(2*t*t-zpos*zpos));
     ang = ang*kRaddeg;
     AliMatrix (idrotm[nrot], 90., 0., 90.-ang,90.,ang, 270.);
     ycoor = -29./2.+ Space ; //2 cm over front plate
     ycoor += (1-(UpDown+1)/2)*Gap;
     zcoor = zpos+(zFLT3/2.+7+zFLT2/2); // Moves to the system of the centre of the modulus FLT2
     gMC->Gspos("FSTR",i, "FLT2", 0., ycoor, zcoor,idrotm[nrot], "ONLY");
     ang  = ang/kRaddeg;
     zpos = zpos - (zSenStrip/2)/TMath::Cos(ang)+UpDown*Gap*TMath::Tan(ang)-(zSenStrip/2)/TMath::Cos(ang);
     last = StripWidth*TMath::Cos(ang)/2;
     UpDown*=-1;
     i++; 
  } while (zpos-(StripWidth/2)*TMath::Cos(ang)>-t+zFLT1+7);

  // Plate  1
  zpos = -t+zFLT1+3.5;
  ang = atan(zpos/sqrt(2*t*t-zpos*zpos));
  Offset = StripWidth*TMath::Cos(ang)/2.;
  zpos -= Offset;
  nrot = 0;
  i=0;
  ycoor= -29./2.+Space+Gap/2;

 do {
     ang = atan(zpos/sqrt(2*t*t-zpos*zpos));
     ang = ang*kRaddeg;
     AliMatrix (idrotm[nrot], 90., 0., 90.-ang,90.,ang, 270.);
     i++;
     zcoor = zpos+(zFLT1/2+zFLT2+zFLT3/2+7.*2.);
     gMC->Gspos("FSTR",i, "FLT1", 0., ycoor, zcoor,idrotm[nrot], "ONLY");
     ang  = ang /kRaddeg;
     zpos = zpos - zSenStrip/TMath::Cos(ang);
     last = StripWidth*TMath::Cos(ang)/2.;
  }  while (zpos>-t+7.+last);

printf("#######################################################\n");
printf("     Distance from the bound of the FLT3: zFLT3- %f cm \n", t+zpos-(zSenStrip/2)/TMath::Cos(ang));
     ang = atan(zpos/sqrt(2*t*t-zpos*zpos));
     zpos = zpos - zSenStrip/TMath::Cos(ang);
printf("NEXT Distance from the bound of the FLT3: zFLT3- %f cm \n", t+zpos-(zSenStrip/2)/TMath::Cos(ang));
printf("#######################################################\n");

////////// Layers after detector /////////////////

// Honeycomb layer after (3cm)

  Float_t OverSpace = Space + 7.3;
///  StripWidth*TMath::Sin(ang) + 1.3;

  par[0] = -1;
  par[1] = 0.6;
  par[2] = -1;
  ycoor = -yFLT/2 + OverSpace + par[1];
  gMC->Gsvolu("FPE1", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPE2", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPE3", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

// Electronics (Cu) after
  ycoor += par[1];
  par[0] = -1;
  par[1] = 1.43*0.05 / 2.; // 5% of X0
  par[2] = -1;
  ycoor += par[1];
  gMC->Gsvolu("FEC1", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEC2", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEC3", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");

// Cooling water after
  ycoor += par[1];
  par[0] = -1;
  par[1] = 36.1*0.02 / 2.; // 2% of X0
  par[2] = -1;
  ycoor += par[1];
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
void AliTOFv2::DrawModule()
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
void AliTOFv2::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //
  AliTOF::CreateMaterials();
}
 
//_____________________________________________________________________________
void AliTOFv2::Init()
{
  //
  // Initialise the detector after the geometry has been defined
  //
  printf("**************************************"
	 "  TOF  "
	 "**************************************\n");
  printf("\n     Version 2 of TOF initialing, "
	 "with openings for PHOS and RICH in symmetric frame\n\n");

  AliTOF::Init();

  //
  // Check that FRAME is there otherwise we have no place where to
  // put TOF
  AliModule* FRAME=gAlice->GetModule("FRAME");
  if(!FRAME) {
    Error("Ctor","TOF needs FRAME to be present\n");
    exit(1);
  } else 
    if(FRAME->IsVersion()!=1) {
      Error("Ctor","FRAME version 1 needed with this version of TOF\n");
      exit(1);
    }

  fIdFTO2=gMC->VolId("FTO2");
  fIdFTO3=gMC->VolId("FTO3");
  fIdFLT1=gMC->VolId("FLT1");
  fIdFLT2=gMC->VolId("FLT2");
  fIdFLT3=gMC->VolId("FLT3");
  printf("**************************************"
	 "  TOF  "
	 "**************************************\n");
}
 
//_____________________________________________________________________________
void AliTOFv2::StepManager()
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

