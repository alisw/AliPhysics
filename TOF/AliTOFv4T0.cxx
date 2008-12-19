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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the functions for version 4 of the Time Of Flight    //
//  detector.                                                                //
//                                                                           //
//  VERSION WITH 5 MODULES AND TILTED STRIPS                                 //
//                                                                           //
//   FULL COVERAGE VERSION +OPTION for PHOS holes                            //
//                                                                           //
//   Author:                                                                 //
//   Fabrizio Pierella                                                       //
//   University of Bologna - Italy                                           //
//                                                                           //
//                                                                           //
//Begin_Html                                                                 //
/*                                                                           //
<img src="picts/AliTOFv4T0Class.gif">                                        //
*/                                                                           //
//End_Html                                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TBRIK.h"
#include "TGeometry.h"
#include "TLorentzVector.h"
#include "TNode.h"
#include "TVirtualMC.h"

#include "AliConst.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliTrackReference.h"

#include "AliTOFGeometry.h"
#include "AliTOFGeometryV4.h"
#include "AliTOFv4T0.h"

extern TDirectory *gDirectory;
extern TVirtualMC *gMC;

extern AliRun *gAlice;

ClassImp(AliTOFv4T0)

//_____________________________________________________________________________
  AliTOFv4T0::AliTOFv4T0():
  fIdFTOA(-1),
  fIdFTOB(-1),
  fIdFTOC(-1),
  fIdFLTA(-1),
  fIdFLTB(-1),
  fIdFLTC(-1),
  fTOFHoles(kFALSE)
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv4T0::AliTOFv4T0(const char *name, const char *title):
  AliTOF(name,title,"tzero"),
  fIdFTOA(-1),
  fIdFTOB(-1),
  fIdFTOC(-1),
  fIdFLTA(-1),
  fIdFLTB(-1),
  fIdFLTC(-1),
  fTOFHoles(kFALSE)
{
  //
  // Standard constructor
  //
  //
  // Check that FRAME is there otherwise we have no place where to
  // put TOF


  AliModule* frame = (AliModule*)gAlice->GetModule("FRAME");
  if(!frame) {
    AliFatal("TOF needs FRAME to be present");
  } else{
    
    if (fTOFGeometry) delete fTOFGeometry;
    fTOFGeometry = new AliTOFGeometryV4();

    if(frame->IsVersion()==1) {
      AliInfo(Form("Frame version %d", frame->IsVersion())); 
      AliInfo("Full Coverage for TOF");
      fTOFHoles=false;}    
    else {
      AliInfo(Form("Frame version %d", frame->IsVersion())); 
      AliInfo("TOF with Holes for PHOS");
      fTOFHoles=true;}      
  }
  fTOFGeometry->SetHoles(fTOFHoles);

  // Save the geometry
  TDirectory* saveDir = gDirectory;
  AliRunLoader::GetRunLoader()->CdGAFile();
  fTOFGeometry->Write("TOFgeometry");
  saveDir->cd();

} 
 
//_____________________________________________________________________________
void AliTOFv4T0::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv4T0.gif">
  */
  //End_Html
  //
  // Creates common geometry
  //
  AliTOF::CreateGeometry();
}
 

//_____________________________________________________________________________
void AliTOFv4T0::TOFpc(Float_t xtof, Float_t ytof, Float_t zlenC,
		     Float_t zlenB, Float_t zlenA, Float_t ztof0)
{
  //
  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xFLT, yFLT, zFLT - sizes of TOF modules (large)

  Float_t  ycoor;
  Float_t  par[3];
  Int_t    *idtmed = fIdtmed->GetArray()-499;
  Int_t    idrotm[100];
  Int_t    nrot = 0;

  Float_t radius = fTOFGeometry->Rmin()+2.;//cm

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
  if(!fTOFHoles)gMC->Gspos("FTOA", 0, "BTO2", 0, zcor3,  0, idrotm[0], "ONLY");
  gMC->Gspos("FTOA", 0, "BTO3", 0, zcor3,  0, idrotm[0], "ONLY");
  

  Float_t db = 0.5; // cm
  Float_t xFLT, xFST, yFLT, zFLTA, zFLTB, zFLTC;
  
  xFLT = fTOFGeometry->StripLength();
  yFLT = ytof;
  zFLTA = zlenA;
  zFLTB = zlenB;
  zFLTC = zlenC;
  
  xFST = xFLT - dynamic_cast<AliTOFGeometryV4*>(fTOFGeometry)->DeadBndX()*2.; // cm
  
  // Sizes of MRPC pads
  
  Float_t yPad = 0.505; //cm
  
  // Large not sensitive volumes with Insensitive Freon
  par[0] = xFLT*0.5;
  par[1] = yFLT*0.5;
  
  AliDebug(1, "************************* TOF geometry **************************");
  
  par[2] = (zFLTA *0.5);
  gMC->Gsvolu("FLTA", "BOX ", idtmed[512], par, 3); // Insensitive Freon
  gMC->Gspos ("FLTA", 0, "FTOA", 0., 0., 0., 0, "ONLY");
  
  par[2] = (zFLTB * 0.5);
  gMC->Gsvolu("FLTB", "BOX ", idtmed[512], par, 3); // Insensitive Freon
  gMC->Gspos ("FLTB", 0, "FTOB", 0., 0., 0., 0, "ONLY");
  
  par[2] = (zFLTC * 0.5);
  gMC->Gsvolu("FLTC", "BOX ", idtmed[512], par, 3); // Insensitive Freon
  gMC->Gspos ("FLTC", 0, "FTOC", 0., 0., 0., 0, "ONLY");
  
  ///// Layers of Aluminum before and after detector /////
  ///// Aluminum Box for Modules (1.8 mm thickness)  /////
  ///// lateral walls not simulated for the time being
  //    const Float_t khAlWall = 0.18;
  // fp to be checked
  const Float_t khAlWall = 0.11;
  par[0] = xFLT*0.5;
  par[1] = khAlWall/2.; // cm
  ycoor = -yFLT/2 + par[1];
  par[2] = (zFLTA *0.5);
  gMC->Gsvolu("FALA", "BOX ", idtmed[508], par, 3); // Alluminium
  gMC->Gspos ("FALA", 1, "FLTA", 0., ycoor, 0., 0, "ONLY");
  gMC->Gspos ("FALA", 2, "FLTA", 0.,-ycoor, 0., 0, "ONLY");
  par[2] = (zFLTB *0.5);
  gMC->Gsvolu("FALB", "BOX ", idtmed[508], par, 3); // Alluminium 
  gMC->Gspos ("FALB", 1, "FLTB", 0., ycoor, 0., 0, "ONLY");
  gMC->Gspos ("FALB", 2, "FLTB", 0.,-ycoor, 0., 0, "ONLY");
  par[2] = (zFLTC *0.5);
  gMC->Gsvolu("FALC", "BOX ", idtmed[508], par, 3); // Alluminium
  gMC->Gspos ("FALC", 1, "FLTC", 0., ycoor, 0., 0, "ONLY");
  gMC->Gspos ("FALC", 2, "FLTC", 0.,-ycoor, 0., 0, "ONLY");
  
  ///////////////// Detector itself //////////////////////
  
  const Float_t  kdeadBound  =  dynamic_cast<AliTOFGeometryV4*>(fTOFGeometry)->DeadBndZ(); //cm non-sensitive between the pad edge 
  //and the boundary of the strip
  const Int_t    knx    = fTOFGeometry->NpadX();  // number of pads along x
  const Int_t    knz    = fTOFGeometry->NpadZ();  // number of pads along z
  
  Float_t zSenStrip  = fTOFGeometry->ZPad() * fTOFGeometry->NpadZ(); // cm
  Float_t stripWidth = zSenStrip + 2*kdeadBound;
  
  par[0] = xFLT*0.5;
  par[1] = yPad*0.5;
  par[2] = stripWidth*0.5;
  
  // new description for strip volume -double stack strip-
  // -- all constants are expressed in cm
  // heigth of different layers
  const Float_t khhony = 0.8     ;   // heigth of HONY  Layer
  const Float_t khpcby = 0.08    ;   // heigth of PCB   Layer
  const Float_t khmyly = 0.035   ;   // heigth of MYLAR Layer
  const Float_t khgraphy = 0.02  ;   // heigth of GRAPHITE Layer
  const Float_t khglasseiy = 0.135;  // 0.6 Ext. Glass + 1.1 i.e. (Int. Glass/2) (mm)
  const Float_t khsensmy = 0.11  ;   // heigth of Sensitive Freon Mixture
  const Float_t kwsensmz = 2*3.5 ;   // cm
  const Float_t klsensmx = 48*2.5;   // cm
  const Float_t kwpadz = 3.5;   // cm z dimension of the FPAD volume
  const Float_t klpadx = 2.5;   // cm x dimension of the FPAD volume
  
  // heigth of the FSTR Volume (the strip volume)
  const Float_t khstripy = 2*khhony+3*khpcby+4*(khmyly+khgraphy+khglasseiy)+2*khsensmy;
  // width  of the FSTR Volume (the strip volume)
  const Float_t kwstripz = 10.;
  // length of the FSTR Volume (the strip volume)
  const Float_t klstripx = 122.;
  
  Float_t parfp[3]={klstripx*0.5,khstripy*0.5,kwstripz*0.5};
  // Coordinates of the strip center in the strip reference frame;
  // used for positioninG internal strip volumes
  Float_t posfp[3]={0.,0.,0.};  
  
  
  // FSTR volume definition-filling this volume with non sensitive Gas Mixture
  gMC->Gsvolu("FSTR","BOX",idtmed[512],parfp,3);
  //-- HONY Layer definition
  //  parfp[0] = -1;
  parfp[1] = khhony*0.5;
  //  parfp[2] = -1;
  gMC->Gsvolu("FHON","BOX",idtmed[503],parfp,3);
  // positioning 2 HONY Layers on FSTR volume
  
  posfp[1]=-khstripy*0.5+parfp[1];
  gMC->Gspos("FHON",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FHON",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  
  //-- PCB Layer definition 
 
  parfp[1] = khpcby*0.5;
  gMC->Gsvolu("FPCB","BOX",idtmed[504],parfp,3);
  // positioning 2 PCB Layers on FSTR volume
  posfp[1]=-khstripy*0.5+khhony+parfp[1];
  gMC->Gspos("FPCB",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FPCB",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  // positioning the central PCB layer
  gMC->Gspos("FPCB",3,"FSTR",0.,0.,0.,0,"ONLY");
  
  
  
  //-- MYLAR Layer definition

  parfp[1] = khmyly*0.5;
  gMC->Gsvolu("FMYL","BOX",idtmed[511],parfp,3);
  // positioning 2 MYLAR Layers on FSTR volume
  posfp[1] = -khstripy*0.5+khhony+khpcby+parfp[1];
  gMC->Gspos("FMYL",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FMYL",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  // adding further 2 MYLAR Layers on FSTR volume
  posfp[1] = khpcby*0.5+parfp[1];
  gMC->Gspos("FMYL",3,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FMYL",4,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  
  
  //-- Graphite Layer definition
 
  parfp[1] = khgraphy*0.5;
  gMC->Gsvolu("FGRP","BOX",idtmed[502],parfp,3);
  // positioning 2 Graphite Layers on FSTR volume
  posfp[1] = -khstripy*0.5+khhony+khpcby+khmyly+parfp[1];
  gMC->Gspos("FGRP",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGRP",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  // adding further 2 Graphite Layers on FSTR volume
  posfp[1] = khpcby*0.5+khmyly+parfp[1];
  gMC->Gspos("FGRP",3,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGRP",4,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  
  
  //-- Glass (EXT. +Semi INT.) Layer definition

  parfp[1] = khglasseiy*0.5;
  gMC->Gsvolu("FGLA","BOX",idtmed[514],parfp,3);
  // positioning 2 Glass Layers on FSTR volume
  posfp[1] = -khstripy*0.5+khhony+khpcby+khmyly+khgraphy+parfp[1];
  gMC->Gspos("FGLA",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGLA",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  // adding further 2 Glass Layers on FSTR volume
  posfp[1] = khpcby*0.5+khmyly+khgraphy+parfp[1];
  gMC->Gspos("FGLA",3,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGLA",4,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  
  
  //-- Sensitive Mixture Layer definition
 
  parfp[0] = klsensmx*0.5;
  parfp[1] = khsensmy*0.5;
  parfp[2] = kwsensmz*0.5;
  gMC->Gsvolu("FSEN","BOX",idtmed[513],parfp,3);
  gMC->Gsvolu("FNSE","BOX",idtmed[512],parfp,3);
  // positioning 2 gas Layers on FSTR volume
  // the upper is insensitive freon
  // while the remaining is sensitive
  posfp[1] = khpcby*0.5+khmyly+khgraphy+khglasseiy+parfp[1];
  gMC->Gspos("FNSE",0,"FSTR", 0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FSEN",0,"FSTR", 0.,-posfp[1],0.,0,"ONLY");
  
  // dividing FSEN along z in knz=2 and along x in knx=48

  gMC->Gsdvn("FSEZ","FSEN",knz,3);
  gMC->Gsdvn("FSEX","FSEZ",knx,1);
  
  // FPAD volume definition

  parfp[0] = klpadx*0.5;   
  parfp[1] = khsensmy*0.5;
  parfp[2] = kwpadz*0.5;
  gMC->Gsvolu("FPAD","BOX",idtmed[513],parfp,3);
  // positioning the FPAD volumes on previous divisions
  gMC->Gspos("FPAD",0,"FSEX",0.,0.,0.,0,"ONLY");
  

  ///////////////////Positioning A module//////////////////////////


  for(Int_t istrip =0; istrip < fTOFGeometry->NStripA(); istrip++){

    Float_t ang = fTOFGeometry->GetAngles(2,istrip);
    AliMatrix (idrotm[0],90.,0.,90.-ang,90.,-ang, 90.);  
    ang /= kRaddeg;
    Float_t zpos = tan(ang)*radius;
    Float_t ypos= fTOFGeometry->GetHeights(2,istrip);
    gMC->Gspos("FSTR",fTOFGeometry->NStripA()-istrip,"FLTA",0.,ypos, zpos,idrotm[0],  "ONLY");
    AliDebug(1, Form("y = %f,  z = %f, , z coord = %f, Rot ang = %f, St. %2i",ypos,zpos,tan(ang)*radius ,ang*kRaddeg,istrip));
  }

  
  ///////////////////Positioning B module//////////////////////////

  for(Int_t istrip =0; istrip < fTOFGeometry->NStripB(); istrip++){

    Float_t ang = fTOFGeometry->GetAngles(3,istrip);
    AliMatrix (idrotm[0],90.,0.,90.-ang,90.,-ang, 90.);  
    ang /= kRaddeg;
    Float_t zpos = tan(ang)*radius+(zFLTA*0.5+zFLTB*0.5+db);
    Float_t ypos= fTOFGeometry->GetHeights(3,istrip);
    gMC->Gspos("FSTR",istrip+1,"FLTB",0.,ypos, zpos,idrotm[nrot],  "ONLY");
    AliDebug(1, Form("y = %f,  z = %f, , z coord = %f, Rot ang = %f, St. %2i",ypos,zpos,tan(ang)*radius,ang*kRaddeg,istrip));
  }

  
  ///////////////////Positioning C module//////////////////////////

  for(Int_t istrip =0; istrip < fTOFGeometry->NStripC(); istrip++){

    Float_t ang = fTOFGeometry->GetAngles(4,istrip);
    AliMatrix (idrotm[0],90.,0.,90.-ang,90.,-ang, 90.);  
    ang /= kRaddeg;
    Float_t zpos = tan(ang)*radius+(zFLTC*0.5+zFLTB+zFLTA*0.5+db*2);
    Float_t ypos= fTOFGeometry->GetHeights(4,istrip);
    gMC->Gspos("FSTR",istrip+1,"FLTC",0.,ypos, zpos,idrotm[nrot],  "ONLY");
    AliDebug(1, Form("y = %f,  z = %f, z coord = %f, Rot ang = %f, St. %2i",ypos,zpos,tan(ang)*radius,ang*kRaddeg,istrip));
  }
   
  ////////// Layers after strips /////////////////
  // Al Layer thickness (2.3mm) factor 0.7
  
  Float_t overSpace = dynamic_cast<AliTOFGeometryV4*>(fTOFGeometry)->OverSpc();//cm
  
  par[0] = xFLT*0.5;
  par[1] = 0.115*0.7; // factor 0.7
  par[2] = (zFLTA *0.5);
  ycoor = -yFLT/2 + overSpace + par[1];
  gMC->Gsvolu("FPEA", "BOX ", idtmed[508], par, 3); // Al
  gMC->Gspos ("FPEA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  par[2] = (zFLTB *0.5);
  gMC->Gsvolu("FPEB", "BOX ", idtmed[508], par, 3); // Al
  gMC->Gspos ("FPEB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  par[2] = (zFLTC *0.5);
  gMC->Gsvolu("FPEC", "BOX ", idtmed[508], par, 3); // Al
  gMC->Gspos ("FPEC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");
  
  
  // plexiglass thickness: 1.5 mm ; factor 0.3

  ycoor += par[1];
  par[0] = xFLT*0.5;
  par[1] = 0.075*0.3; // factor 0.3 
  par[2] = (zFLTA *0.5);
  ycoor += par[1];
  gMC->Gsvolu("FECA", "BOX ", idtmed[505], par, 3); // Plexigl.
  gMC->Gspos ("FECA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  par[2] = (zFLTB *0.5);
  gMC->Gsvolu("FECB", "BOX ", idtmed[505], par, 3); // Plexigl.
  gMC->Gspos ("FECB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  par[2] = (zFLTC *0.5);
  gMC->Gsvolu("FECC", "BOX ", idtmed[505], par, 3); // Plexigl.
  gMC->Gspos ("FECC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");
  
  // frame of Air

  ycoor += par[1];
  par[0] = xFLT*0.5;
  par[1] = (yFLT/2-ycoor-khAlWall)*0.5; // Aluminum layer considered (0.18 cm)
  par[2] = (zFLTA *0.5);
  ycoor += par[1];
  gMC->Gsvolu("FAIA", "BOX ", idtmed[500], par, 3); // Air
  gMC->Gspos ("FAIA", 0, "FLTA", 0., ycoor, 0., 0, "ONLY");
  par[2] = (zFLTB *0.5);
  gMC->Gsvolu("FAIB", "BOX ", idtmed[500], par, 3); // Air
  gMC->Gspos ("FAIB", 0, "FLTB", 0., ycoor, 0., 0, "ONLY");
  par[2] = (zFLTC *0.5);
  gMC->Gsvolu("FAIC", "BOX ", idtmed[500], par, 3); // Air
  gMC->Gspos ("FAIC", 0, "FLTC", 0., ycoor, 0., 0, "ONLY");
  
  
  // start with cards and cooling tubes
  // finally, cards, cooling tubes and layer for thermal dispersion
  // 3 volumes
  // card volume definition
  
  // see GEOM200 in GEANT manual


  AliMatrix(idrotm[98], 90., 0., 90., 90., 0., 0.); // 0 deg
  
  Float_t cardpar[3];
  cardpar[0]= 61.;
  cardpar[1]= 5.;
  cardpar[2]= 0.1;
  gMC->Gsvolu("FCAR", "BOX ", idtmed[504], cardpar, 3); // PCB Card 
  //alu plate volume definition
  cardpar[1]= 3.5;
  cardpar[2]= 0.05;
  gMC->Gsvolu("FALP", "BOX ", idtmed[508], cardpar, 3); // Alu Plate
  
  
  // central module positioning (FAIA)
  Float_t cardpos[3], aplpos2, stepforcardA=6.625;
  cardpos[0]= 0.;
  cardpos[1]= -0.5;
  cardpos[2]= -53.;
  Float_t aplpos1 = -2.;
  Int_t icard;
  for (icard=0; icard < fTOFGeometry->NStripA(); ++icard) {
    cardpos[2]= cardpos[2]+stepforcardA;
    aplpos2 = cardpos[2]+0.15;
    gMC->Gspos("FCAR",icard,"FAIA",cardpos[0],cardpos[1],cardpos[2],idrotm[98],"ONLY");
    gMC->Gspos("FALP",icard,"FAIA",cardpos[0],aplpos1,aplpos2,idrotm[98],"ONLY");
    
  }
  
  
  // intermediate module positioning (FAIB)
  Float_t stepforcardB= 7.05;
  cardpos[2]= -70.5;
  for (icard=0; icard < fTOFGeometry->NStripB(); ++icard) {
    cardpos[2]= cardpos[2]+stepforcardB;
    aplpos2 = cardpos[2]+0.15; 
    gMC->Gspos("FCAR",icard,"FAIB",cardpos[0],cardpos[1],cardpos[2],idrotm[98],"ONLY");
    gMC->Gspos("FALP",icard,"FAIB",cardpos[0],aplpos1,aplpos2,idrotm[98],"ONLY");
  }
  
  
  // outer module positioning (FAIC)
  Float_t stepforcardC= 8.45238;
  cardpos[2]= -88.75;
  for (icard=0; icard < fTOFGeometry->NStripC(); ++icard) {
    cardpos[2]= cardpos[2]+stepforcardC;
    aplpos2 = cardpos[2]+0.15;
    gMC->Gspos("FCAR",icard,"FAIC",cardpos[0],cardpos[1],cardpos[2],idrotm[98],"ONLY");
    gMC->Gspos("FALP",icard,"FAIC",cardpos[0],aplpos1,aplpos2,idrotm[98],"ONLY");
  }
  
  // tube volume definition

  Float_t tubepar[3];
  tubepar[0]= 0.;
  tubepar[1]= 0.4;
  tubepar[2]= 61.;
  gMC->Gsvolu("FTUB", "TUBE", idtmed[516], tubepar, 3); // cooling tubes (steel)
  tubepar[0]= 0.;
  tubepar[1]= 0.35;
  tubepar[2]= 61.;
  gMC->Gsvolu("FITU", "TUBE", idtmed[515], tubepar, 3); // cooling water
  // positioning water tube into the steel one
  gMC->Gspos("FITU",1,"FTUB",0.,0.,0.,0,"ONLY");
  
  
  // rotation matrix
  AliMatrix(idrotm[99], 180., 90., 90., 90., 90., 0.);
  // central module positioning (FAIA)
  Float_t tubepos[3], tdis=0.6;
  tubepos[0]= 0.;
  tubepos[1]= cardpos[1];
  tubepos[2]= -53.+tdis;
  //  tub1pos = 5.;
  Int_t itub;
  for (itub=0; itub < fTOFGeometry->NStripA(); ++itub) {
    tubepos[2]= tubepos[2]+stepforcardA;
    gMC->Gspos("FTUB",itub,"FAIA",tubepos[0],tubepos[1],tubepos[2],idrotm[99],
	       "ONLY");
  }
  
  
  // intermediate module positioning (FAIB)
  tubepos[2]= -70.5+tdis;
  for (itub=0; itub < fTOFGeometry->NStripB(); ++itub) {
    tubepos[2]= tubepos[2]+stepforcardB;
    gMC->Gspos("FTUB",itub,"FAIB",tubepos[0],tubepos[1],tubepos[2],idrotm[99],
	       "ONLY");
  }
  
  // outer module positioning (FAIC)
  tubepos[2]= -88.75+tdis;
  for (itub=0; itub < fTOFGeometry->NStripC(); ++itub) {
    tubepos[2]= tubepos[2]+stepforcardC;
    gMC->Gspos("FTUB",itub,"FAIC",tubepos[0],tubepos[1],tubepos[2],idrotm[99],
	       "ONLY");
  }

}
//_____________________________________________________________________________
void AliTOFv4T0::DrawModule() const
{
  //
  // Draw a shaded view of the Time Of Flight version 4
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
void AliTOFv4T0::DrawDetectorModules() const
{
//
// Draw a shaded view of the TOF detector version 4
//
 
 
//Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);

//
//Set volumes visible
// 
//=====> Level 1
  // Level 1 for TOF volumes
  gMC->Gsatt("B077","seen",0);
 
 
//==========> Level 2
  // Level 2
  gMC->Gsatt("B076","seen",-1); // all B076 sub-levels skipped -
  gMC->Gsatt("B071","seen",0);
  gMC->Gsatt("B074","seen",0);
  gMC->Gsatt("B075","seen",0);
  gMC->Gsatt("B080","seen",0); // B080 does not has sub-level                


  // Level 2 of B071
  gMC->Gsatt("B063","seen",-1); // all B063 sub-levels skipped   -
  gMC->Gsatt("B065","seen",-1); // all B065 sub-levels skipped   -
  gMC->Gsatt("B067","seen",-1); // all B067 sub-levels skipped   -
  gMC->Gsatt("B069","seen",-1); // all B069 sub-levels skipped   -
  gMC->Gsatt("B056","seen",0);  // B056 does not has sub-levels  -
  gMC->Gsatt("B059","seen",-1); // all B059 sub-levels skipped   -
  gMC->Gsatt("B072","seen",-1); // all B072 sub-levels skipped   -
  gMC->Gsatt("BTR1","seen",0);  // BTR1 do not have sub-levels   -
  gMC->Gsatt("BTO1","seen",0);

 
  // Level 2 of B074
  gMC->Gsatt("BTR2","seen",0); // BTR2 does not has sub-levels -
  gMC->Gsatt("BTO2","seen",0);

  // Level 2 of B075
  gMC->Gsatt("BTR3","seen",0); // BTR3 do not have sub-levels -
  gMC->Gsatt("BTO3","seen",0);

// ==================> Level 3
  // Level 3 of B071 / Level 2 of BTO1
  gMC->Gsatt("FTOC","seen",-2);
  gMC->Gsatt("FTOB","seen",-2);
  gMC->Gsatt("FTOA","seen",-2);
 
  // Level 3 of B074 / Level 2 of BTO2
  // -> cfr previous settings
 
  // Level 3 of B075 / Level 2 of BTO3
  // -> cfr previous settings

  gMC->Gdopt("hide","on");
  gMC->Gdopt("shad","on");
  gMC->Gsatt("*", "fill", 5);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, 0, 1000, 0, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 45, 40, 0, 10, 10, .015, .015);
  gMC->Gdhead(1111,"TOF detector V1");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}                                 

//_____________________________________________________________________________
void AliTOFv4T0::DrawDetectorStrips() const
{
  //
  // Draw a shaded view of the TOF strips for version 4
  //
  
  //Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  
  //
  //Set volumes visible 
  //=====> Level 1
  // Level 1 for TOF volumes
  gMC->Gsatt("B077","seen",0);
  
  //==========> Level 2
  // Level 2
  gMC->Gsatt("B076","seen",-1); // all B076 sub-levels skipped -
  gMC->Gsatt("B071","seen",0);
  gMC->Gsatt("B074","seen",0);
  gMC->Gsatt("B075","seen",0);
  gMC->Gsatt("B080","seen",0); // B080 does not has sub-level
  
  // Level 2 of B071
  gMC->Gsatt("B063","seen",-1); // all B063 sub-levels skipped   -
  gMC->Gsatt("B065","seen",-1); // all B065 sub-levels skipped   -
  gMC->Gsatt("B067","seen",-1); // all B067 sub-levels skipped   -
  gMC->Gsatt("B069","seen",-1); // all B069 sub-levels skipped   -
  gMC->Gsatt("B056","seen",0);  // B056 does not has sub-levels  -
  gMC->Gsatt("B059","seen",-1); // all B059 sub-levels skipped   -
  gMC->Gsatt("B072","seen",-1); // all B072 sub-levels skipped   -
  gMC->Gsatt("BTR1","seen",0);  // BTR1 do not have sub-levels   -
  gMC->Gsatt("BTO1","seen",0);
  
  // ==================> Level 3
  // Level 3 of B071 / Level 2 of BTO1
  gMC->Gsatt("FTOC","seen",0);
  gMC->Gsatt("FTOB","seen",0);
  gMC->Gsatt("FTOA","seen",0);
  
  // Level 3 of B074 / Level 2 of BTO2
  // -> cfr previous settings
  
  // Level 3 of B075 / Level 2 of BTO3
  // -> cfr previous settings
  
  
  // ==========================> Level 4
  // Level 4 of B071 / Level 3 of BTO1 / Level 2 of FTOC
  gMC->Gsatt("FLTC","seen",0);
  // Level 4 of B071 / Level 3 of BTO1 / Level 2 of FTOB
  gMC->Gsatt("FLTB","seen",0);
  // Level 4 of B071 / Level 3 of BTO1 / Level 2 of FTOA
  gMC->Gsatt("FLTA","seen",0);
  
  // Level 4 of B074 / Level 3 of BTO2 / Level 2 of FTOC
  // -> cfr previous settings
  // Level 4 of B074 / Level 3 of BTO2 / Level 2 of FTOB
  // -> cfr previous settings
  
  // Level 4 of B075 / Level 3 of BTO3 / Level 2 of FTOC
  // -> cfr previous settings
  
  //======================================> Level 5
  // Level 5 of B071 / Level 4 of BTO1 / Level 3 of FTOC / Level 2 of FLTC
  gMC->Gsatt("FALC","seen",0); // no children for FALC
  gMC->Gsatt("FSTR","seen",-2);
  gMC->Gsatt("FPEC","seen",0); // no children for FPEC
  gMC->Gsatt("FECC","seen",0); // no children for FECC
  gMC->Gsatt("FWAC","seen",0); // no children for FWAC
  gMC->Gsatt("FAIC","seen",0); // no children for FAIC
  
  // Level 5 of B071 / Level 4 of BTO1 / Level 3 of FTOB / Level 2 of FLTB
  gMC->Gsatt("FALB","seen",0); // no children for FALB
  //-->  gMC->Gsatt("FSTR","seen",-2);
  
  
  // -> cfr previous settings
  gMC->Gsatt("FPEB","seen",0); // no children for FPEB
  gMC->Gsatt("FECB","seen",0); // no children for FECB
  gMC->Gsatt("FWAB","seen",0); // no children for FWAB
  gMC->Gsatt("FAIB","seen",0); // no children for FAIB
  
  // Level 5 of B071 / Level 4 of BTO1 / Level 3 of FTOA / Level 2 of FLTA
  gMC->Gsatt("FALA","seen",0); // no children for FALB
  //-->  gMC->Gsatt("FSTR","seen",-2);
  // -> cfr previous settings
  gMC->Gsatt("FPEA","seen",0); // no children for FPEA
  gMC->Gsatt("FECA","seen",0); // no children for FECA
  gMC->Gsatt("FWAA","seen",0); // no children for FWAA
  gMC->Gsatt("FAIA","seen",0); // no children for FAIA
  
  // Level 2 of B074
  gMC->Gsatt("BTR2","seen",0); // BTR2 does not has sub-levels -
  gMC->Gsatt("BTO2","seen",0);
  
  // Level 2 of B075
  gMC->Gsatt("BTR3","seen",0); // BTR3 do not have sub-levels -
  gMC->Gsatt("BTO3","seen",0);
  
  // for others Level 5, cfr. previous settings

  gMC->Gdopt("hide","on");
  gMC->Gdopt("shad","on");
  gMC->Gsatt("*", "fill", 5);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, 0, 1000, 0, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 45, 40, 0, 10, 10, .015, .015);
  gMC->Gdhead(1111,"TOF Strips V1");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTOFv4T0::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //
  //AliTOF::CreateMaterials();

  //
  // Defines TOF materials for all versions
  // Revision: F. Pierella 18-VI-2002
  //

  AliMagF *magneticField = (AliMagF*)gAlice->Field();

  Int_t   isxfld = magneticField->Integ();
  Float_t sxmgmx = magneticField->Max();

  //--- Quartz (SiO2) to simulate float glass
  //    density tuned to have correct float glass 
  //    radiation length
  Float_t   aq[2] = { 28.0855,15.9994 };
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  Float_t   dq = 2.55; // std value: 2.2
  Int_t nq = -2;

  // --- Freon C2F4H2 (TOF-TDR pagg.)
  // Geant Manual CONS110-1, pag. 43 (Geant, Detector Description and Simulation Tool)
  Float_t afre[3]  = {12.011,18.998,1.007};
  Float_t zfre[3]  = { 6., 9., 1.}; 
  Float_t wfre[3]  = { 2., 4., 2.};
  Float_t densfre  = 0.00375;   
// http://www.fi.infn.it/sezione/prevprot/gas/freon.html
  Int_t nfre = -3; 
/*
  //-- Isobutane quencher C4H10 (5% in the sensitive mixture)
  Float_t aiso[2]  = {12.011,1.007};
  Float_t ziso[2]  = { 6.,  1.};
  Float_t wiso[2]  = { 4., 10.};
  Float_t densiso  = .......;  // (g/cm3) density
  Int_t nfre = -2; // < 0 i.e. proportion by number of atoms of each kind
  //-- SF6 (5% in the sensitive mixture)
  Float_t asf[3]  = {32.066,18.998};
  Float_t zsf[3]  = { 16., 9.};
  Float_t wsf[3]  = {  1., 6.}; 
  Float_t denssf  = .....;   // (g/cm3) density
  Int_t nfre = -2; // < 0 i.e. proportion by number of atoms of each kind
*/
  // --- CO2 
  Float_t ac[2]   = {12.,16.};
  Float_t zc[2]   = { 6., 8.};
  Float_t wc[2]   = { 1., 2.};
  Float_t dc = .001977;
  Int_t nc = -2;
   // For mylar (C5H4O2) 
  Float_t amy[3] = { 12., 1., 16. };
  Float_t zmy[3] = {  6., 1.,  8. };
  Float_t wmy[3] = {  5., 4.,  2. };
  Float_t dmy    = 1.39;
  Int_t nmy = -3;
 // For polyethilene (CH2) - honeycomb -
  Float_t ape[2] = { 12., 1. };
  Float_t zpe[2] = {  6., 1. };
  Float_t wpe[2] = {  1., 2. };
  Float_t dpe    = 0.935*0.479; //To have 1%X0 for 1cm as for honeycomb
  Int_t npe = -2;
  // --- G10 
  Float_t ag10[4] = { 12.,1.,16.,28. };
  Float_t zg10[4] = {  6.,1., 8.,14. };
  Float_t wmatg10[4] = { .259,.288,.248,.205 };
  Float_t densg10  = 1.7;
  Int_t nlmatg10 = -4;

  // plexiglass CH2=C(CH3)CO2CH3
  Float_t aplex[3] = { 12.,1.,16.};
  Float_t zplex[3] = {  6.,1., 8.};
  Float_t wmatplex[3] = {5.,8.,2.};
  Float_t densplex  =1.16;
  Int_t nplex = -3;

  // ---- ALUMINA (AL203) 
  Float_t aal[2] = { 27.,16.};
  Float_t zal[2] = { 13., 8.};
  Float_t wmatal[2] = { 2.,3. };
  Float_t densal  = 2.3;
  Int_t nlmatal = -2;
  // -- Water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;
  Int_t nwa = -2;

// stainless steel
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };

  //AliMaterial(0, "Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);

  // AIR
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;

  AliMixture( 1, "Air$", aAir, zAir, dAir, 4, wAir);

  AliMaterial( 2, "Cu $",  63.54, 29.0, 8.96, 1.43, 14.8);
  AliMaterial( 3, "C  $",  12.01,  6.0, 2.265,18.8, 74.4);
  AliMixture ( 4, "Polyethilene$", ape, zpe, dpe, npe, wpe);
  AliMixture ( 5, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  AliMixture ( 6, "PLE$", aplex, zplex, densplex, nplex, wmatplex);
  AliMixture ( 7, "CO2$", ac, zc, dc, nc, wc);
  AliMixture ( 8, "ALUMINA$", aal, zal, densal, nlmatal, wmatal);
  AliMaterial( 9, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(10, "C-TRD$", 12.01, 6., 2.265*18.8/69.282*15./100, 18.8, 74.4); // for 15%
  AliMixture (11, "Mylar$",  amy, zmy, dmy, nmy, wmy);
  AliMixture (12, "Freon$",  afre, zfre, densfre, nfre, wfre);
  AliMixture (13, "Glass$", aq, zq, dq, nq, wq);
  AliMixture (14, "Water$",  awa, zwa, dwa, nwa, wwa);
  AliMixture (15, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);

  Float_t epsil, stmin, deemax, stemax;
 
  //   Previous data
  //       EPSIL  =  0.1   ! Tracking precision, 
  //       STEMAX = 0.1      ! Maximum displacement for multiple scattering
  //       DEEMAX = 0.1    ! Maximum fractional energy loss, DLS 
  //       STMIN  = 0.1 
  //
  //   New data  
  epsil  = .001;  // Tracking precision,
  stemax = -1.;   // Maximum displacement for multiple scattering
  deemax = -.3;   // Maximum fractional energy loss, DLS
  stmin  = -.8;

  AliMedium( 1, "Air$"  ,  1, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 2, "Cu $"  ,  2, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 3, "C  $"  ,  3, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 4, "Pol$"  ,  4, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 5, "G10$"  ,  5, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 6, "PLE$"  ,  6, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 7, "CO2$"  ,  7, 0, isxfld, sxmgmx, 10., -.01, -.1, .01, -.01);
  AliMedium( 8,"ALUMINA$", 8, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 9,"Al Frame$",9, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(10, "DME-S$",  6, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(11, "C-TRD$", 10, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(12, "Myl$"  , 11, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(13, "Fre$"  , 12, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(14, "Fre-S$", 12, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(15, "Glass$", 13, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(16, "Water$", 14, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(17, "STEEL$", 15, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);

}

//_____________________________________________________________________________
void AliTOFv4T0::Init()
{
  //
  // Initialise the detector after the geometry has been defined
  //
  AliDebug(1, "**************************************"
           "  TOF  "
           "**************************************");
  AliDebug(1, "  Version 4 of TOF initialing, "
	   "symmetric TOF - Full Coverage version");
  
  AliTOF::Init();
  
  fIdFTOA = gMC->VolId("FTOA");
  fIdFTOB = gMC->VolId("FTOB");
  fIdFTOC = gMC->VolId("FTOC");
  fIdFLTA = gMC->VolId("FLTA");
  fIdFLTB = gMC->VolId("FLTB");
  fIdFLTC = gMC->VolId("FLTC");

  AliDebug(1, "**************************************"
           "  TOF  "
           "**************************************");
}
 
//_____________________________________________________________________________
void AliTOFv4T0::StepManager()
{

  //
  // Procedure called at each step in the Time Of Flight
  //

  TLorentzVector mom, pos;
  Float_t xm[3],pm[3],xpad[3],ppad[3];
  Float_t hits[14];
  Int_t   vol[5];
  Int_t   sector, plate, padx, padz, strip;
  Int_t   copy, padzid, padxid, stripid, i;
  Int_t   *idtmed = fIdtmed->GetArray()-499;
  Float_t incidenceAngle;
      
  if(
     gMC->IsTrackEntering()
     && gMC->TrackCharge()
     //&& gMC->GetMedium()==idtmed[513]
     && gMC->CurrentMedium()==idtmed[513]
     && gMC->CurrentVolID(copy)==fIdSens
     )
  {

    AliMC *mcApplication = (AliMC*)gAlice->GetMCApp();

    AddTrackReference(mcApplication->GetCurrentTrackNumber(), AliTrackReference::kTOF);
    //AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());

    // getting information about hit volumes
    
    padzid=gMC->CurrentVolOffID(2,copy);
    padz=copy; 

    
    padxid=gMC->CurrentVolOffID(1,copy);
    padx=copy; 
    
    stripid=gMC->CurrentVolOffID(4,copy);
    strip=copy; 

    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);


    //    Double_t NormPos=1./pos.Rho();

    Double_t normMom=1./mom.Rho();

    //  getting the cohordinates in pad ref system

    xm[0] = (Float_t)pos.X();
    xm[1] = (Float_t)pos.Y();
    xm[2] = (Float_t)pos.Z();

    pm[0] = (Float_t)mom.X()*normMom;
    pm[1] = (Float_t)mom.Y()*normMom;
    pm[2] = (Float_t)mom.Z()*normMom;
 
    gMC->Gmtod(xm,xpad,1);
    gMC->Gmtod(pm,ppad,2);

    
    if (TMath::Abs(ppad[1])>1) {
      AliWarning("Abs(ppad) > 1");
      ppad[1]=TMath::Sign((Float_t)1,ppad[1]);
    }
    incidenceAngle = TMath::ACos(ppad[1])*kRaddeg;


    const char * pathA="FTOA";
    const char * pathB="FTOB";
    const char * pathC="FTOC";
    const char * path71="B071";
    const char * path75="B075";
    const char * path74="B074";
    const char* volpath;    

    Int_t index=0;
    volpath=gMC->CurrentVolOffName(6);
    index=gMC->CurrentVolOffID(6,copy);
    index=copy;

    
    plate=-1;
    if(strcmp(pathC,volpath)==0 && index==1)plate=0;
    if(strcmp(pathB,volpath)==0 && index==1)plate=1;
    if(strcmp(pathA,volpath)==0 && index==0)plate=2;
    if(strcmp(pathB,volpath)==0 && index==2)plate=3;
    if(strcmp(pathC,volpath)==0 && index==2)plate=4;



    if (plate==0) strip=fTOFGeometry->NStripC()-strip;
    else if (plate==1) strip=fTOFGeometry->NStripB()-strip;
    else strip--;
 
    //Apply ALICE conventions for volume numbering increasing with theta, phi

    if (plate==3 || plate==4){
      padx=fTOFGeometry->NpadX()-padx;
      padz=fTOFGeometry->NpadZ()-padz;
      xpad[0]=-xpad[0];      
      xpad[2]=-xpad[2];      
    }
    else {
     padx--;
     padz--;
    }



    volpath=gMC->CurrentVolOffName(8);
    index=gMC->CurrentVolOffID(8,copy);
    index=copy;

    sector=-1;
    if(strcmp(path71,volpath)==0 && index <6) sector=12+index;
    if(strcmp(path71,volpath)==0 && index >=6) sector=index-3;
    if(strcmp(path75,volpath)==0) sector=index-1;
    if(strcmp(path74,volpath)==0) sector=10+index;

    for(i=0;i<3;++i) {
      hits[i]   = pos[i];
      hits[i+3] = pm[i];
    }

    hits[6] = mom.Rho();
    hits[7] = pos[3];
    hits[8] = xpad[0];
    hits[9] = xpad[1];
    hits[10]= xpad[2];
    hits[11]= incidenceAngle;
    hits[12]= gMC->Edep();
    hits[13]= gMC->TrackLength();
    
    vol[0]= sector;
    vol[1]= plate;
    vol[2]= strip;
    vol[3]= padx;
    vol[4]= padz;    

    AddT0Hit(mcApplication->GetCurrentTrackNumber(),vol, hits);
    //AddT0Hit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol, hits);
  }
}
