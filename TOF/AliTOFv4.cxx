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
Revision 1.6  2000/05/10 16:52:18  vicinanz
New TOF version with holes for PHOS/RICH

Revision 1.4.2.1  2000/05/10 09:37:16  vicinanz
New version with Holes for PHOS/RICH

Revision 1.14  1999/11/05 22:39:06  fca
New hits structure

Revision 1.13  1999/11/02 11:26:39  fca
added stdlib.h for exit

Revision 1.12  1999/11/01 20:41:57  fca
Added protections against using the wrong version of FRAME

Revision 1.11  1999/10/22 08:04:14  fca
Correct improper use of negative parameters

Revision 1.10  1999/10/16 19:30:06  fca
Corrected Rotation Matrix and CVS log

Revision 1.9  1999/10/15 15:35:20  fca
New version for frame1099 with and without holes

Revision 1.8  1999/09/29 09:24:33  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight: design of C.Williams                  
//
//  This class contains the functions for version 1 of the Time Of Flight    //
//  detector.                                                                //
//
//  VERSION WITH 5 MODULES AND TILTED STRIPS 
//  
//   FULL COVERAGE VERSION
//
//   Authors:
//
//   Alessio Seganti
//   Domenico Vicinanza
//
//   University of Salerno - Italy
//
//
//Begin_Html
/*
<img src="picts/AliTOFv4Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <stdlib.h>

#include "AliTOFv4.h"
#include "TBRIK.h"
#include "TGeometry.h"
#include "TNode.h"
#include "TObject.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

 
ClassImp(AliTOFv4)
 
//_____________________________________________________________________________
AliTOFv4::AliTOFv4()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv4::AliTOFv4(const char *name, const char *title)
        : AliTOF(name,title)
{
  //
  // Standard constructor
  //
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
 
}

//_____________________________________________________________________________
void AliTOFv4::BuildGeometry()
{
  //
  // Build TOF ROOT geometry for the ALICE event display
  //
  TNode *Node, *Top;
  const int kColorTOF  = 27;

  // Find top TNODE
  Top = gAlice->GetGeometry()->GetNode("alice");

  // Position the different copies
  const Float_t rTof  =(fRmax+fRmin)/2;
  const Float_t hTof  = fRmax-fRmin;
  const Int_t   fNTof = 18;
  const Float_t kPi   = TMath::Pi();
  const Float_t angle = 2*kPi/fNTof;
  Float_t ang;

  // Define TOF basic volume
  
  char NodeName0[6], NodeName1[6], NodeName2[6]; 
  char NodeName3[6], NodeName4[6], RotMatNum[6];

  new TBRIK("S_TOF_C","TOF box","void",
            120*0.5,hTof*0.5,fZlenC*0.5);
  new TBRIK("S_TOF_B","TOF box","void",
            120*0.5,hTof*0.5,fZlenB*0.5);
  new TBRIK("S_TOF_A","TOF box","void",
            120*0.5,hTof*0.5,fZlenA*0.5);

  for (Int_t NodeNum=1;NodeNum<19;NodeNum++){
     
      if (NodeNum<10) {
           sprintf(RotMatNum,"rot50%i",NodeNum);
           sprintf(NodeName0,"FTO00%i",NodeNum);
           sprintf(NodeName1,"FTO10%i",NodeNum);
           sprintf(NodeName2,"FTO20%i",NodeNum);
           sprintf(NodeName3,"FTO30%i",NodeNum);
           sprintf(NodeName4,"FTO40%i",NodeNum);
      }
      if (NodeNum>9) {
           sprintf(RotMatNum,"rot5%i",NodeNum);
           sprintf(NodeName0,"FTO0%i",NodeNum);
           sprintf(NodeName1,"FTO1%i",NodeNum);
           sprintf(NodeName2,"FTO2%i",NodeNum);
           sprintf(NodeName3,"FTO3%i",NodeNum);
           sprintf(NodeName4,"FTO4%i",NodeNum);
      }
 
      new TRotMatrix(RotMatNum,RotMatNum,90,-20*NodeNum,90,90-20*NodeNum,0,0);
      ang = (4.5-NodeNum) * angle;

      Top->cd();
      Node = new TNode(NodeName0,NodeName0,"S_TOF_C",rTof*TMath::Cos(ang),rTof*TMath::Sin(ang),299.15,RotMatNum);
      Node->SetLineColor(kColorTOF);
      fNodes->Add(Node); 

      Top->cd(); 
      Node = new TNode(NodeName1,NodeName1,"S_TOF_C",rTof*TMath::Cos(ang),rTof*TMath::Sin(ang),-299.15,RotMatNum);
      Node->SetLineColor(kColorTOF);
      fNodes->Add(Node); 

      Top->cd();
      Node = new TNode(NodeName2,NodeName2,"S_TOF_B",rTof*TMath::Cos(ang),rTof*TMath::Sin(ang),146.45,RotMatNum);
      Node->SetLineColor(kColorTOF);
      fNodes->Add(Node); 

      Top->cd();
      Node = new TNode(NodeName3,NodeName3,"S_TOF_B",rTof*TMath::Cos(ang),rTof*TMath::Sin(ang),-146.45,RotMatNum);
      Node->SetLineColor(kColorTOF);
      fNodes->Add(Node); 

      Top->cd();
      Node = new TNode(NodeName4,NodeName4,"S_TOF_A",rTof*TMath::Cos(ang),rTof*TMath::Sin(ang),0.,RotMatNum);
      Node->SetLineColor(kColorTOF);
      fNodes->Add(Node); 
  }
}


 
//_____________________________________________________________________________
void AliTOFv4::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv4.gif">
  */
  //End_Html
  //
  // Creates common geometry
  //
  AliTOF::CreateGeometry();
}
 
//_____________________________________________________________________________
void AliTOFv4::TOFpc(Float_t xtof, Float_t ytof, Float_t zlenC,
		     Float_t zlenB, Float_t zlenA, Float_t ztof0)
{
  //
  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xFLT, yFLT, zFLT - sizes of TOF modules (large)
  
  Float_t  ycoor, zcoor;
  Float_t  par[10];
  Int_t    *idtmed = fIdtmed->GetArray()-499;
  Int_t    idrotm[100];
  Int_t    nrot = 0;
  Float_t  hTof = fRmax-fRmin;
  
  Float_t Radius = fRmin+2.;//cm

  par[0] =  xtof * 0.5;
  par[1] =  ytof * 0.5;
  par[2] = zlenC * 0.5;
  gMC->Gsvolu("FTOC", "BOX ", idtmed[506], par, 3);
  par[2] = zlenB * 0.5;
  gMC->Gsvolu("FTOB", "BOX ", idtmed[506], par, 3);
  par[2] = zlenA * 0.5;
  gMC->Gsvolu("FTOA", "BOX ", idtmed[506], par, 3);


// Positioning of modules

   Float_t zcor1 = ztof0 - zlenC*0.5;
   Float_t zcor2 = ztof0 - zlenC - zlenB*0.5;
   Float_t zcor3 = 0.;

   AliMatrix(idrotm[0], 90.,  0., 0., 0., 90,-90.);
   AliMatrix(idrotm[1], 90.,180., 0., 0., 90, 90.);
   gMC->Gspos("FTOC", 1, "BTO1", 0,  zcor1, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTOC", 2, "BTO1", 0, -zcor1, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTOC", 1, "BTO2", 0,  zcor1, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTOC", 2, "BTO2", 0, -zcor1, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTOC", 1, "BTO3", 0,  zcor1, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTOC", 2, "BTO3", 0, -zcor1, 0, idrotm[1], "ONLY");

   gMC->Gspos("FTOB", 1, "BTO1", 0,  zcor2, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTOB", 2, "BTO1", 0, -zcor2, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTOB", 1, "BTO2", 0,  zcor2, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTOB", 2, "BTO2", 0, -zcor2, 0, idrotm[1], "ONLY");
   gMC->Gspos("FTOB", 1, "BTO3", 0,  zcor2, 0, idrotm[0], "ONLY");
   gMC->Gspos("FTOB", 2, "BTO3", 0, -zcor2, 0, idrotm[1], "ONLY");

   gMC->Gspos("FTOA", 0, "BTO1", 0, zcor3,  0, idrotm[0], "ONLY");
   gMC->Gspos("FTOA", 0, "BTO2", 0, zcor3,  0, idrotm[0], "ONLY");
   gMC->Gspos("FTOA", 0, "BTO3", 0, zcor3,  0, idrotm[0], "ONLY");

  Float_t db = 0.5;//cm
  Float_t xFLT, xFST, yFLT, zFLTA, zFLTB, zFLTC;

  xFLT = fStripLn;
  yFLT = ytof;
  zFLTA = zlenA;
  zFLTB = zlenB;
  zFLTC = zlenC;

  xFST = xFLT-fDeadBndX*2;//cm

// Sizes of MRPC pads

  Float_t yPad = 0.505;//cm 
  
// Large not sensitive volumes with CO2 
  par[0] = xFLT*0.5;
  par[1] = yFLT*0.5;

  cout <<"************************* TOF geometry **************************"<<endl;

  par[2] = (zFLTA *0.5);
  gMC->Gsvolu("FLTA", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos ("FLTA", 0, "FTOA", 0., 0., 0., 0, "ONLY");

  par[2] = (zFLTB * 0.5);
  gMC->Gsvolu("FLTB", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos ("FLTB", 0, "FTOB", 0., 0., 0., 0, "ONLY");

  par[2] = (zFLTC * 0.5); 
  gMC->Gsvolu("FLTC", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos ("FLTC", 0, "FTOC", 0., 0., 0., 0, "ONLY");

////////// Layers before detector ////////////////////

// MYlar layer in front 1.0 mm thick at the beginning
  par[0] = -1;
  par[1] = 0.1;//cm
  par[2] = -1;
  ycoor = -yFLT/2 + par[1];
  gMC->Gsvolu("FMYA", "BOX ", idtmed[508], par, 3); // Alluminium
  gMC->Gspos ("FMYA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FMYB", "BOX ", idtmed[508], par, 3); // Alluminium
  gMC->Gspos ("FMYB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FMYC", "BOX ", idtmed[508], par, 3); // Alluminium 
  gMC->Gspos ("FMYC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");

// honeycomb (special Polyethilene Layer of 1cm)
  ycoor = ycoor + par[1];
  par[0] = -1;
  par[1] = 0.5;//cm
  par[2] = -1;
  ycoor = ycoor + par[1];
  gMC->Gsvolu("FPLA", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FPLA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPLB", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FPLB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPLC", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FPLC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");

///////////////// Detector itself //////////////////////

  const Float_t  DeadBound  =  fDeadBndZ; //cm non-sensitive between the pad edge 
                                          //and the boundary of the strip
  const Int_t    nx    = fNpadX;          // number of pads along x
  const Int_t    nz    = fNpadZ;          // number of pads along z
  const Float_t  Space = fSpace;            //cm distance from the front plate of the box

  Float_t zSenStrip  = fZpad*fNpadZ;//cm
  Float_t StripWidth = zSenStrip + 2*DeadBound;

  par[0] = xFLT*0.5;
  par[1] = yPad*0.5; 
  par[2] = StripWidth*0.5;
  
  // glass layer of detector STRip
  gMC->Gsvolu("FSTR","BOX",idtmed[514],par,3);

  // Non-Sesitive Freon boundaries
  par[0] =  xFLT*0.5;
  par[1] =  0.110*0.5;//cm
  par[2] = -1;
  gMC->Gsvolu("FNSF","BOX",idtmed[512],par,3);
  gMC->Gspos ("FNSF",0,"FSTR",0.,0.,0.,0,"ONLY");

  // MYlar for Internal non-sesitive boundaries
//  par[1] = 0.025;//cm
//  gMC->Gsvolu("FMYI","BOX",idtmed[510],par,3); 
//  gMC->Gspos ("FMYI",0,"FNSF",0.,0.,0.,0,"MANY");

  // MYlar eXternal layers
  par[1] = 0.035*0.5;//cm
  ycoor = -yPad*0.5+par[1];
  gMC->Gsvolu("FMYX","BOX",idtmed[510],par,3);
  gMC->Gspos ("FMYX",1,"FSTR",0.,ycoor,0.,0,"ONLY");
  gMC->Gspos ("FMYX",2,"FSTR",0.,-ycoor,0.,0,"ONLY");
  ycoor += par[1];
 
  // GRaphyte Layers
  par[1] = 0.003*0.5;
  ycoor += par[1];
  gMC->Gsvolu("FGRL","BOX",idtmed[502],par,3);
  gMC->Gspos ("FGRL",1,"FSTR",0.,ycoor,0.,0,"ONLY");
  gMC->Gspos ("FGRL",2,"FSTR",0.,-ycoor,0.,0,"ONLY");

  // freon sensitive layer (Chlorine-Fluorine-Carbon)
  par[0] = xFST*0.5;
  par[1] =  0.110*0.5;
  par[2] = zSenStrip*0.5;
  gMC->Gsvolu("FCFC","BOX",idtmed[513],par,3);
  gMC->Gspos ("FCFC",0,"FNSF",0.,0.,0.,0,"ONLY");
  
  // Pad definition x & z
  gMC->Gsdvn("FLZ","FCFC", nz, 3); 
  gMC->Gsdvn("FLX","FLZ" , nx, 1); 

  // MRPC PAD itself 
  par[0] = -1;
  par[1] = -1; 
  par[2] = -1;
  gMC->Gsvolu("FPAD", "BOX ", idtmed[513], par, 3);
  gMC->Gspos ("FPAD", 0, "FLX", 0., 0., 0., 0, "ONLY");

////  Positioning the Strips  (FSTR) in the FLT volumes  /////

  // Plate A (Central) 
  
  Float_t t = zFLTC+zFLTB+zFLTA*0.5+ 2*db;//Half Width of Barrel

  Float_t Gap  = fGapA; //cm  distance between the strip axis
  Float_t zpos = 0;
  Float_t ang  = 0;
  Int_t i=1,j=1;
  nrot  = 0;
  zcoor = 0;
  ycoor = -14.5 + Space ; //2 cm over front plate

  AliMatrix (idrotm[0],  90.,  0.,90.,90.,0., 90.);   
  gMC->Gspos("FSTR",j,"FLTA",0.,ycoor, 0.,idrotm[0],"ONLY");

     printf("%f,  St. %2i, Pl.3 ",ang*kRaddeg,i); 
     printf("y = %f,  z = %f, zpos = %f \n",ycoor,zcoor,zpos);

  zcoor -= zSenStrip;
  j++;
  Int_t UpDown = -1; // UpDown=-1 -> Upper strip
                     // UpDown=+1 -> Lower strip
  do{
     ang = atan(zcoor/Radius);
     ang *= kRaddeg;
     AliMatrix (idrotm[nrot],  90.,  0.,90.-ang,90.,-ang, 90.);   
     AliMatrix (idrotm[nrot+1],90.,180.,90.+ang,90., ang, 90.);
     ang /= kRaddeg;
     ycoor = -14.5+ Space; //2 cm over front plate
     ycoor += (1-(UpDown+1)/2)*Gap;
     gMC->Gspos("FSTR",j  ,"FLTA",0.,ycoor, zcoor,idrotm[nrot],  "ONLY");
     gMC->Gspos("FSTR",j+1,"FLTA",0.,ycoor,-zcoor,idrotm[nrot+1],"ONLY");

     printf("%f,  St. %2i, Pl.3 ",ang*kRaddeg,i); 
     printf("y = %f,  z = %f, zpos = %f \n",ycoor,zcoor,zpos);

     j += 2;
     UpDown*= -1; // Alternate strips 
     zcoor = zcoor-(zSenStrip/2)/TMath::Cos(ang)-
             UpDown*Gap*TMath::Tan(ang)-
	     (zSenStrip/2)/TMath::Cos(ang);
  } while (zcoor-(StripWidth/2)*TMath::Cos(ang)>-t+zFLTC+zFLTB+db*2);
  
  zcoor = zcoor+(zSenStrip/2)/TMath::Cos(ang)+
          UpDown*Gap*TMath::Tan(ang)+
          (zSenStrip/2)/TMath::Cos(ang);

  Gap = fGapB;
  zcoor = zcoor-(zSenStrip/2)/TMath::Cos(ang)-
          UpDown*Gap*TMath::Tan(ang)-
          (zSenStrip/2)/TMath::Cos(ang);

  ang = atan(zcoor/Radius);
  ang *= kRaddeg;
  AliMatrix (idrotm[nrot],  90.,  0.,90.-ang,90.,-ang, 90.);   
  AliMatrix (idrotm[nrot+1],90.,180.,90.+ang,90., ang, 90.);
  ang /= kRaddeg;
	  
  ycoor = -14.5+ Space; //2 cm over front plate
  ycoor += (1-(UpDown+1)/2)*Gap;
  gMC->Gspos("FSTR",j  ,"FLTA",0.,ycoor, zcoor,idrotm[nrot],  "ONLY");
  gMC->Gspos("FSTR",j+1,"FLTA",0.,ycoor,-zcoor,idrotm[nrot+1],"ONLY");

     printf("%f,  St. %2i, Pl.3 ",ang*kRaddeg,i); 
     printf("y = %f,  z = %f, zpos = %f \n",ycoor,zcoor,zpos);

  ycoor = -hTof/2.+ Space;//2 cm over front plate

  // Plate  B

  nrot = 0;
  i=1;
  UpDown = 1;
  Float_t DeadRegion = 1.0;//cm
  
  zpos = zcoor - (zSenStrip/2)/TMath::Cos(ang)-
         UpDown*Gap*TMath::Tan(ang)-
	 (zSenStrip/2)/TMath::Cos(ang)-
	 DeadRegion/TMath::Cos(ang);

  ang = atan(zpos/Radius);
  ang *= kRaddeg;
  AliMatrix (idrotm[nrot], 90., 0., 90.-ang,90.,ang, 270.);
  ang /= kRaddeg;
  ycoor = -hTof*0.5+ Space ; //2 cm over front plate
  ycoor += (1-(UpDown+1)/2)*Gap;
  zcoor = zpos+(zFLTA*0.5+zFLTB*0.5+db); // Moves to the system of the modulus FLTB
  gMC->Gspos("FSTR",i, "FLTB", 0., ycoor, zcoor,idrotm[nrot], "ONLY");

     printf("%f,  St. %2i, Pl.4 ",ang*kRaddeg,i); 
     printf("y = %f,  z = %f, zpos = %f \n",ycoor,zcoor,zpos);

  i++;
  UpDown*=-1;

  do {
     zpos = zpos - (zSenStrip/2)/TMath::Cos(ang)-
            UpDown*Gap*TMath::Tan(ang)-
	    (zSenStrip/2)/TMath::Cos(ang);
     ang = atan(zpos/Radius);
     ang *= kRaddeg;
     AliMatrix (idrotm[nrot], 90., 0., 90.-ang,90.,ang, 270.);
     ang /= kRaddeg;
     ycoor = -hTof*0.5+ Space ; //2 cm over front plate
     ycoor += (1-(UpDown+1)/2)*Gap;
     zcoor = zpos+(zFLTA*0.5+zFLTB*0.5+db); // Moves to the system of the modulus FLTB
     gMC->Gspos("FSTR",i, "FLTB", 0., ycoor, zcoor,idrotm[nrot], "ONLY");

     printf("%f,  St. %2i, Pl.4 ",ang*kRaddeg,i); 
     printf("y = %f,  z = %f, zpos = %f \n",ycoor,zcoor,zpos);

     UpDown*=-1;
     i++;
  } while (TMath::Abs(ang*kRaddeg)<22.5);
  //till we reach a tilting angle of 22.5 degrees

  ycoor = -hTof*0.5+ Space ; //2 cm over front plate
  zpos = zpos - zSenStrip/TMath::Cos(ang);

  do {
     ang = atan(zpos/Radius);
     ang *= kRaddeg;
     AliMatrix (idrotm[nrot], 90., 0., 90.-ang,90.,ang, 270.);
     ang /= kRaddeg;
     zcoor = zpos+(zFLTB/2+zFLTA/2+db);
     gMC->Gspos("FSTR",i, "FLTB", 0., ycoor, zcoor,idrotm[nrot], "ONLY");
     zpos = zpos - zSenStrip/TMath::Cos(ang);
     printf("%f,  St. %2i, Pl.4 ",ang*kRaddeg,i); 
     printf("y = %f,  z = %f, zpos = %f \n",ycoor,zcoor,zpos);
     i++;

  }  while (zpos-StripWidth*0.5/TMath::Cos(ang)>-t+zFLTC+db);

  // Plate  C
  
  zpos = zpos + zSenStrip/TMath::Cos(ang);

  zpos = zpos - (zSenStrip/2)/TMath::Cos(ang)+
         Gap*TMath::Tan(ang)-
	 (zSenStrip/2)/TMath::Cos(ang);

  nrot = 0;
  i=0;
  ycoor= -hTof*0.5+Space+Gap;

  do {
     i++;
     ang = atan(zpos/Radius);
     ang *= kRaddeg;
     AliMatrix (idrotm[nrot], 90., 0., 90.-ang,90.,ang, 270.);
     ang /= kRaddeg;
     zcoor = zpos+(zFLTC*0.5+zFLTB+zFLTA*0.5+db*2);
     gMC->Gspos("FSTR",i, "FLTC", 0., ycoor, zcoor,idrotm[nrot], "ONLY");

     printf("%f,  St. %2i, Pl.5 ",ang*kRaddeg,i); 
     printf("y = %f,  z = %f, zpos = %f \n",ycoor,zcoor,zpos);

     zpos = zpos - zSenStrip/TMath::Cos(ang);
  }  while (zpos-StripWidth*TMath::Cos(ang)*0.5>-t);


////////// Layers after detector /////////////////

// honeycomb (Polyethilene) Layer after (3cm)

  Float_t OverSpace = fOverSpc;//cm

  par[0] = -1;
  par[1] = 0.6;
  par[2] = -1;
  ycoor = -yFLT/2 + OverSpace + par[1];
  gMC->Gsvolu("FPEA", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FPEA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPEB", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FPEB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPEC", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FPEC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");

// Electronics (Cu) after
  ycoor += par[1];
  par[0] = -1;
  par[1] = 1.43*0.05*0.5; // 5% of X0
  par[2] = -1;
  ycoor += par[1];
  gMC->Gsvolu("FECA", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos ("FECA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FECB", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos ("FECB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FECC", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos ("FECC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");

// cooling WAter after
  ycoor += par[1];
  par[0] = -1;
  par[1] = 36.1*0.02*0.5; // 2% of X0
  par[2] = -1;
  ycoor += par[1];
  gMC->Gsvolu("FWAA", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos ("FWAA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FWAB", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos ("FWAB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FWAC", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos ("FWAC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");

//Back Plate honycomb (2cm)
  par[0] = -1;
  par[1] = 2 *0.5;
  par[2] = -1;
  ycoor = yFLT/2 - par[1];
  gMC->Gsvolu("FBPA", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FBPA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FBPB", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FBPB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FBPC", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos ("FBPC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");
}

//_____________________________________________________________________________
void AliTOFv4::DrawModule()
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

  gMC->Gsatt("FTOA","SEEN",1);
  gMC->Gsatt("FTOB","SEEN",1);
  gMC->Gsatt("FTOC","SEEN",1);
  gMC->Gsatt("FLTA","SEEN",1);
  gMC->Gsatt("FLTB","SEEN",1);
  gMC->Gsatt("FLTC","SEEN",1);
  gMC->Gsatt("FPLA","SEEN",1);
  gMC->Gsatt("FPLB","SEEN",1);
  gMC->Gsatt("FPLC","SEEN",1);
  gMC->Gsatt("FSTR","SEEN",1);
  gMC->Gsatt("FPEA","SEEN",1);
  gMC->Gsatt("FPEB","SEEN",1);
  gMC->Gsatt("FPEC","SEEN",1);
  
  gMC->Gsatt("FLZ1","SEEN",0);
  gMC->Gsatt("FLZ2","SEEN",0);
  gMC->Gsatt("FLZ3","SEEN",0);
  gMC->Gsatt("FLX1","SEEN",0);
  gMC->Gsatt("FLX2","SEEN",0);
  gMC->Gsatt("FLX3","SEEN",0);
  gMC->Gsatt("FPAD","SEEN",0);

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
void AliTOFv4::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //
  AliTOF::CreateMaterials();
}
 
//_____________________________________________________________________________
void AliTOFv4::Init()
{
  //
  // Initialise the detector after the geometry has been defined
  //
  printf("**************************************"
	 "  TOF  "
	 "**************************************\n");
  printf("\n   Version 4 of TOF initialing, "
	      "symmetric TOF - Full Coverage version\n");

  AliTOF::Init();

  fIdFTOA = gMC->VolId("FTOA");
  fIdFTOB = gMC->VolId("FTOB");
  fIdFTOC = gMC->VolId("FTOC");
  fIdFLTA = gMC->VolId("FLTA");
  fIdFLTB = gMC->VolId("FLTB");
  fIdFLTC = gMC->VolId("FLTC");

  printf("**************************************"
	 "  TOF  "
	 "**************************************\n");
}
 
//_____________________________________________________________________________
void AliTOFv4::StepManager()
{
  //
  // Procedure called at each step in the Time Of Flight
  //
  TLorentzVector mom, pos;
  Float_t xm[3],pm[3],xpad[3],ppad[3];
  Float_t hits[13],phi,phid,z;
  Int_t   vol[5];
  Int_t   sector, plate, pad_x, pad_z, strip;
  Int_t   copy, pad_z_id, pad_x_id, strip_id, i;
  Int_t   *idtmed = fIdtmed->GetArray()-499;
  Float_t IncidenceAngle;
  
  if(gMC->GetMedium()==idtmed[513] && 
     gMC->IsTrackEntering() && gMC->TrackCharge()
     && gMC->CurrentVolID(copy)==fIdSens) 
  {    
    // getting information about hit volumes
    
    pad_z_id=gMC->CurrentVolOffID(2,copy);
    pad_z=copy;  
    
    pad_x_id=gMC->CurrentVolOffID(1,copy);
    pad_x=copy;  
    
    strip_id=gMC->CurrentVolOffID(5,copy);
    strip=copy;  

    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);

//    Double_t NormPos=1./pos.Rho();
    Double_t NormMom=1./mom.Rho();

//  getting the cohordinates in pad ref system
    xm[0] = (Float_t)pos.X();
    xm[1] = (Float_t)pos.Y();
    xm[2] = (Float_t)pos.Z();

    pm[0] = (Float_t)mom.X()*NormMom;
    pm[1] = (Float_t)mom.Y()*NormMom;
    pm[2] = (Float_t)mom.Z()*NormMom;
 
    gMC->Gmtod(xm,xpad,1);
    gMC->Gmtod(pm,ppad,2);

    IncidenceAngle = TMath::ACos(ppad[1])*kRaddeg;

    z = pos[2];

    plate = 0;   
    if (TMath::Abs(z) <=  fZlenA*0.5)  plate = 3;
    if (z < (fZlenA*0.5+fZlenB) && 
        z >  fZlenA*0.5)               plate = 4;
    if (z >-(fZlenA*0.5+fZlenB) &&
        z < -fZlenA*0.5)               plate = 2;
    if (z > (fZlenA*0.5+fZlenB))       plate = 5;
    if (z <-(fZlenA*0.5+fZlenB))       plate = 1;

    phi = pos.Phi();
    phid = phi*kRaddeg+180.;
    sector = Int_t (phid/20.);
    sector++;

    for(i=0;i<3;++i) {
      hits[i]   = pos[i];
      hits[i+3] = pm[i];
    }

    hits[6] = mom.Rho();
    hits[7] = pos[3];
    hits[8] = xpad[0];
    hits[9] = xpad[1];
    hits[10]= xpad[2];
    hits[11]= IncidenceAngle;
    hits[12]= gMC->Edep();
    
    vol[0]= sector;
    vol[1]= plate;
    vol[2]= strip;
    vol[3]= pad_x;
    vol[4]= pad_z;
    
    AddHit(gAlice->CurrentTrack(),vol, hits);
  }
}
