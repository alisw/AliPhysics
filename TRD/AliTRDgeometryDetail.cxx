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
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry for the spaceframe without holes                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliMC.h"

#include "AliTRDgeometryDetail.h"

ClassImp(AliTRDgeometryDetail)

//_____________________________________________________________________________
AliTRDgeometryDetail::AliTRDgeometryDetail():AliTRDgeometryFull()
{
  //
  // AliTRDgeometryDetail default constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDgeometryDetail::~AliTRDgeometryDetail()
{
  //
  // AliTRDgeometryDetail destructor
  //

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::Init()
{
  //
  // Initializes the geometry parameter
  //

  AliTRDgeometryFull::Init();

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::CreateGeometry(Int_t *idtmed)
{
  //
  // Create the detailed TRD geometry without hole
  //

  Int_t iplan;

  const Int_t kNparTrd = 4;
  const Int_t kNparCha = 3;
  const Int_t kNplan   = fgkNplan;

  Float_t parDum[3];
  Float_t parTrd[kNparTrd];
  Float_t parCha[kNparCha];

  Float_t xpos, ypos, zpos;

  // The aluminum frames - readout + electronics (Al)
  // The inner chambers
  gMC->Gsvolu("UAFI","BOX ",idtmed[1301-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UAFM","BOX ",idtmed[1301-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UAFO","BOX ",idtmed[1301-1],parDum,0);

  // The inner part of the aluminum frames (Air)
  // The inner chambers
  gMC->Gsvolu("UAII","BOX ",idtmed[1302-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UAIM","BOX ",idtmed[1302-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UAIO","BOX ",idtmed[1302-1],parDum,0);

  // The carbon frames - radiator + driftchamber (C)
  // The inner chambers
  gMC->Gsvolu("UCFI","BOX ",idtmed[1307-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UCFM","BOX ",idtmed[1307-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UCFO","BOX ",idtmed[1307-1],parDum,0);

  // The inner part of the carbon frames (Air)
  // The inner chambers
  gMC->Gsvolu("UCII","BOX ",idtmed[1302-1],parDum,0);
  // The middle chambers
  gMC->Gsvolu("UCIM","BOX ",idtmed[1302-1],parDum,0);
  // The outer chambers
  gMC->Gsvolu("UCIO","BOX ",idtmed[1302-1],parDum,0);

  // The material layers inside the chambers
  parCha[0] = -1.;
  parCha[1] = -1.;
  // Rohacell layer (radiator)
  parCha[2] = fgkRaThick/2;
  gMC->Gsvolu("UL03","BOX ",idtmed[1315-1],parCha,kNparCha);
  // Mylar layer (entrance window + HV cathode) 
  parCha[2] = fgkMyThick/2;
  gMC->Gsvolu("UL04","BOX ",idtmed[1308-1],parCha,kNparCha);
  // Xe/Isobutane layer (drift volume) 
  parCha[2] = fgkDrThick/2.;
  gMC->Gsvolu("UL05","BOX ",idtmed[1309-1],parCha,kNparCha);
  // Xe/Isobutane layer (amplification volume)
  parCha[2] = fgkAmThick/2.;
  gMC->Gsvolu("UL06","BOX ",idtmed[1309-1],parCha,kNparCha);
  
  // Cu layer (pad plane)
  parCha[2] = fgkCuThick/2;
  gMC->Gsvolu("UL07","BOX ",idtmed[1305-1],parCha,kNparCha);
  // G10 layer (support structure)
  parCha[2] = fgkSuThick/2;
  gMC->Gsvolu("UL08","BOX ",idtmed[1313-1],parCha,kNparCha);

  // Create the readout volumina
  CreateReadout(idtmed);

  // Create the volumina for the cooling
  CreateCooling(idtmed);

  // Position the layers in the chambers
  xpos = 0;
  ypos = 0;

  // Rohacell layer (radiator)
  zpos = fgkRaZpos;
  gMC->Gspos("UL03",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL03",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL03",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Mylar layer (entrance window + HV cathode)   
  zpos = fgkMyZpos;
  gMC->Gspos("UL04",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (drift volume) 
  zpos = fgkDrZpos;
  gMC->Gspos("UL05",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (amplification volume)
  zpos = fgkAmZpos;
  gMC->Gspos("UL06",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Cu layer (pad plane)
  zpos = fgkCuZpos;
  gMC->Gspos("UL07",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // G10 layer (support structure)
  zpos = fgkSuZpos;
  gMC->Gspos("UL08",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",3,"UAIO",xpos,ypos,zpos,0,"ONLY");

  // The TRD mother volume for one sector (Air), full length in z-direction
  parTrd[0] = fgkSwidth1/2.;
  parTrd[1] = fgkSwidth2/2.;
  parTrd[2] = fgkSlenTR1/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("TRD1","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  // Make the aluminum frame of the chamber 0.5cm shorter to acommodate
  // the volumes for the detailed readout electronics
  const Float_t kcaframeOff = 0.5;
  Float_t caframe = fgkCaframe - kcaframeOff;

  // Position the chambers in the TRD mother volume
  for (iplan = 1; iplan <= kNplan; iplan++) {

    // The inner chambers ---------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthI[iplan-1]/2.;
    parCha[2] = caframe/2.;
    xpos      = 0.;
    ypos      = 0.;
    zpos      = fgkCheight    - fgkSheight/2. - kcaframeOff - caframe/2. 
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFI",iplan         ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.   - fgkCathick;
    parCha[1] = fClengthI[iplan-1]/2. - fgkCathick;
    parCha[2] = caframe/2.;
    xpos      = 0.;
    ypos      = 0.;
    zpos      = fgkCheight    - fgkSheight/2. - kcaframeOff - caframe/2.
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAII",iplan         ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthI[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos      = 0.;
    ypos      = 0.;
    zpos      = fgkCcframe/2. - fgkSheight/2. 
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFI",iplan         ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.   - fgkCcthick;
    parCha[1] = fClengthI[iplan-1]/2. - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos      = 0.;
    ypos      = 0.;
    zpos      = fgkCcframe/2. - fgkSheight/2. 
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCII",iplan         ,"TRD1",xpos,ypos,zpos,0,"ONLY",parCha,kNparCha);

    PositionReadout(iplan-1,2);
    PositionCooling(iplan-1,2);

    // The middle chambers --------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM1[iplan-1]/2.;
    parCha[2] = caframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos      = fgkCheight    - fgkSheight/2. - kcaframeOff - caframe/2.
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.      - fgkCathick;
    parCha[1] = fClengthM1[iplan-1]/2.   - fgkCathick;
    parCha[2] = caframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos      = fgkCheight    - fgkSheight/2. - kcaframeOff - caframe/2.
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthM1[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos      = fgkCcframe/2. - fgkSheight/2. 
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.      - fgkCcthick;
    parCha[1] = fClengthM1[iplan-1]/2.   - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos      = fgkCcframe/2. - fgkSheight/2. 
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    PositionReadout(iplan-1,1);
    PositionReadout(iplan-1,3);
    PositionCooling(iplan-1,1);
    PositionCooling(iplan-1,3);

    // The outer chambers ---------------------------------------------------------------

    // the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO1[iplan-1]/2.;
    parCha[2] = caframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos      = fgkCheight    - fgkSheight/2. - kcaframeOff - caframe/2.
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAFO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the inner part of the aluminum frame
    parCha[0] = fCwidth[iplan-1]/2.      - fgkCathick;
    parCha[1] = fClengthO1[iplan-1]/2.   - fgkCathick;
    parCha[2] = caframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos      = fgkCheight    - fgkSheight/2. - kcaframeOff - caframe/2.
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UAIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UAIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.;
    parCha[1] = fClengthO1[iplan-1]/2.;
    parCha[2] = fgkCcframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos      = fgkCcframe/2. - fgkSheight/2. 
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCFO",iplan,         "TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    // the inner part of the carbon frame
    parCha[0] = fCwidth[iplan-1]/2.      - fgkCcthick;
    parCha[1] = fClengthO1[iplan-1]/2.   - fgkCcthick;
    parCha[2] = fgkCcframe/2.;
    xpos      = 0.;
    ypos      = fClengthI[iplan-1]/2. + fClengthM1[iplan-1] + fClengthO1[iplan-1]/2.;
    zpos      = fgkCcframe/2. - fgkSheight/2. 
              + (iplan-1) * (fgkCheight + fgkCspace);
    gMC->Gsposp("UCIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",parCha,kNparCha);
    gMC->Gsposp("UCIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",parCha,kNparCha);

    PositionReadout(iplan-1,0);
    PositionReadout(iplan-1,4);
    PositionCooling(iplan-1,0);
    PositionCooling(iplan-1,4);

  }

  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("TRD1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD1",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD1",3,"BTR3",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::CreateReadout(Int_t *idtmed)
{
  //
  // Create the volumina of the readout electronics
  //

  const Int_t   kNparBox    = 3;

  Float_t parBox[kNparBox];
  Float_t xpos = 0.0;
  Float_t ypos = 0.0;
  Float_t zpos = 0.0;

  // The mother volume for the MCMs + connectors (air)
  parBox[0] = 3.0/2.;
  parBox[1] = 3.4/2.;
  parBox[2] = 0.5/2.;
  gMC->Gsvolu("UMCM","BOX",idtmed[1302-1],parBox,kNparBox);

  // The MCM carrier G10 layer
  parBox[0] = 3.0/2.;
  parBox[1] = 3.0/2.;
  parBox[2] = 0.1/2.;
  gMC->Gsvolu("UMC1","BOX",idtmed[1319-1],parBox,kNparBox);
  // The MCM carrier Cu layer
  parBox[0] = 3.0/2.;
  parBox[1] = 3.0/2.;
  parBox[2] = 0.0034/2.;
  gMC->Gsvolu("UMC2","BOX",idtmed[1318-1],parBox,kNparBox);
  // The MCM carrier Sn layer
  parBox[0] = 3.0/2.;
  parBox[1] = 3.0/2.;
  parBox[2] = 0.004/2.;
  gMC->Gsvolu("UMC3","BOX",idtmed[1317-1],parBox,kNparBox);
  // The MCM carrier Al layer
  parBox[0] = 3.0/2.;
  parBox[1] = 3.0/2.;
  parBox[2] = 0.05/2.;
  gMC->Gsvolu("UMC4","BOX",idtmed[1316-1],parBox,kNparBox);

  // The epoxy of chip no.1
  parBox[0] = 0.548/2.;
  parBox[1] = 0.548/2.;
  parBox[2] = 0.1/2.;
  gMC->Gsvolu("UCE1","BOX",idtmed[1321-1],parBox,kNparBox);
  // The silicon of chip no.1
  parBox[0] = 0.316/2.;
  parBox[1] = 0.316/2.;
  parBox[2] = 0.03/2.;
  gMC->Gsvolu("UCS1","BOX",idtmed[1320-1],parBox,kNparBox);

  // The epoxy of chip no.2
  parBox[0] = 1.549/2.;
  parBox[1] = 1.549/2.;
  parBox[2] = 0.1/2.;
  gMC->Gsvolu("UCE2","BOX",idtmed[1321-1],parBox,kNparBox);
  // The silicon of chip no.2
  parBox[0] = 0.894/2.;
  parBox[1] = 0.894/2.;
  parBox[2] = 0.03/2.;
  gMC->Gsvolu("UCS2","BOX",idtmed[1320-1],parBox,kNparBox);

  // The PE of the connector
  parBox[0] = 2.25/2.;
  parBox[1] = 0.4/2.;
  parBox[2] = 0.3/2.;
  gMC->Gsvolu("UCN1","BOX",idtmed[1322-1],parBox,kNparBox);
  // The Cu of the connector
  parBox[0] = 2.25/2.;
  parBox[1] = 0.4/2.;
  parBox[2] = 0.005/2.;
  gMC->Gsvolu("UCN2","BOX",idtmed[1323-1],parBox,kNparBox);

  xpos  =  0.0;
  ypos  = -0.4/2.;
  zpos  = -0.25      + 0.1/2.;
  gMC->Gspos("UMC1",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  0.1/2.   + 0.0034/2.;
  gMC->Gspos("UMC2",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  0.0034/2 + 0.004/2.;
  gMC->Gspos("UMC3",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  0.004/2  + 0.05/2.;
  gMC->Gspos("UMC4",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  0.05/2.  + 0.1/2.;
  xpos  =  1.0;
  gMC->Gspos("UCE1",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  xpos  = -0.5;
  gMC->Gspos("UCE2",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  zpos +=  0.1/2.   + 0.03/2.;
  xpos  =  1.0;
  gMC->Gspos("UCS1",1,"UMCM",xpos,ypos,zpos,0,"ONLY");
  xpos  = -0.5;
  gMC->Gspos("UCS2",1,"UMCM",xpos,ypos,zpos,0,"ONLY");  
  xpos  =  0.0;
  ypos  =  3.4/2.   - 0.4/2.;
  zpos  = -0.25     + 0.3/2.;
  gMC->Gspos("UCN1",1,"UMCM",xpos,ypos,zpos,0,"ONLY");  
  zpos +=  0.3/2.   + 0.005/2.;
  gMC->Gspos("UCN2",1,"UMCM",xpos,ypos,zpos,0,"ONLY");  

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::PositionReadout(Int_t ipla, Int_t icha)
{
  //
  // Position the volumina inside the readout mother volume
  //

  const Float_t kcaframeOff = 0.5;

  Int_t nMCMrow = GetRowMax(ipla,icha,0);
  Int_t nMCMcol = 6;

  Float_t xSize = GetChamberWidth(ipla) / ((Float_t) nMCMcol);
  Float_t ySize = 0.0;
  Float_t x0    = - GetChamberWidth(ipla) /2. + fgkCcthick;
  Float_t y0    = 0.0;
  switch (icha) {
  case 0:
    ySize =   GetChamberLengthO(ipla) / ((Float_t) nMCMrow);
    y0    =   fClengthI[ipla]/2. + fClengthM1[ipla] + fClengthO1[ipla]/2.
            - GetChamberLengthO(ipla) / 2. + fgkCcthick;
    break;
  case 1:
    ySize =   GetChamberLengthM(ipla) / ((Float_t) nMCMrow);
    y0    =   fClengthI[ipla]/2. + fClengthM1[ipla]/2.
            - GetChamberLengthM(ipla) / 2. + fgkCcthick;
    break;
  case 2:
    ySize =   GetChamberLengthI(ipla) / ((Float_t) nMCMrow);
    y0    = - GetChamberLengthI(ipla) / 2. + fgkCcthick;
    break;
  case 3:
    ySize =   GetChamberLengthM(ipla) / ((Float_t) nMCMrow);
    y0    = - fClengthI[ipla]/2. - fClengthM1[ipla]/2.
            - GetChamberLengthM(ipla) / 2. + fgkCcthick;
    break;
  case 4:
    ySize =   GetChamberLengthO(ipla) / ((Float_t) nMCMrow);
    y0    = - fClengthI[ipla]/2. - fClengthM1[ipla] - fClengthO1[ipla]/2.
            - GetChamberLengthO(ipla) / 2. + fgkCcthick;
    break;
  };

  Int_t iCopy = GetDetector(ipla,icha,0) * 1000;
  for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
    for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) {
      iCopy++;
      Float_t xpos = (0.5 + iMCMcol) * xSize + x0; 
      Float_t ypos = (0.5 + iMCMrow) * ySize + y0;
      Float_t zpos = fgkCheight - fgkSheight/2. - kcaframeOff/2. 
                   + ipla * (fgkCheight + fgkCspace);
      gMC->Gspos("UMCM",iCopy,"TRD1",xpos,ypos,zpos,0,"ONLY");    
    }
  }

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::CreateCooling(Int_t *idtmed)
{
  //
  // Create the volumina of the cooling
  //

  const Int_t kNparBox = 3;

  Float_t parBox[kNparBox];

  // The aluminum pipe for the cooling
  parBox[0] = 0.0;
  parBox[1] = 0.0;
  parBox[2] = 0.0;
  gMC->Gsvolu("UCOA","BOX",idtmed[1324-1],parBox,0);

  // The cooling water
  parBox[0] = -1.;
  parBox[1] =  0.2/2.;
  parBox[2] =  0.2/2.;
  gMC->Gsvolu("UCOW","BOX",idtmed[1314-1],parBox,kNparBox);

  // Water inside he cooling pipe
  Float_t xpos = 0.0;
  Float_t ypos = 0.0;
  Float_t zpos = 0.0;
  gMC->Gspos("UCOW",1,"UCOA",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::PositionCooling(Int_t ipla, Int_t icha)
{
  //
  // Position the volumina of the cooling
  //

  const Int_t   kNpar       = 3;
  
  const Float_t kcaframeOff = 0.5;

  Float_t par[kNpar];
  Float_t xpos    = 0.0;
  Float_t ypos    = 0.0;
  Float_t zpos    = 0.0;

  Int_t   iCopy   = GetDetector(ipla,icha,0) * 100;
  Int_t   nMCMrow = GetRowMax(ipla,icha,0);

  Float_t xSize   = 0.0;
  Float_t x0      = 0.0;
  switch (icha) {
  case 0:
    xSize =   GetChamberLengthO(ipla) / ((Float_t) nMCMrow);
    x0    =   fClengthI[ipla]/2. + fClengthM1[ipla] + fClengthO1[ipla]/2.
            - GetChamberLengthO(ipla) / 2. + fgkCcthick;
    break;
  case 1:
    xSize =   GetChamberLengthM(ipla) / ((Float_t) nMCMrow);
    x0    =   fClengthI[ipla]/2. + fClengthM1[ipla]/2.
            - GetChamberLengthM(ipla) / 2. + fgkCcthick;
    break;
  case 2:
    xSize =   GetChamberLengthI(ipla) / ((Float_t) nMCMrow);
    x0    = - GetChamberLengthI(ipla) / 2. + fgkCcthick;
    break;
  case 3:
    xSize =   GetChamberLengthM(ipla) / ((Float_t) nMCMrow);
    x0    = - fClengthI[ipla]/2. - fClengthM1[ipla]/2.
            - GetChamberLengthM(ipla) / 2. + fgkCcthick;
    break;
  case 4:
    xSize =   GetChamberLengthO(ipla) / ((Float_t) nMCMrow);
    x0    = - fClengthI[ipla]/2. - fClengthM1[ipla] - fClengthO1[ipla]/2.
            - GetChamberLengthO(ipla) / 2. + fgkCcthick;
    break;
  };

  // Position the cooling pipes
  for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {

    xpos   = 0.0;
    ypos   = (0.5 + iMCMrow) * xSize + x0 - 1.9;
    zpos   = fgkCheight - fgkSheight/2. - kcaframeOff/2. 
                        + ipla * (fgkCheight + fgkCspace);
    par[0] = GetChamberWidth(ipla) / 2. - fgkCcthick;
    par[1] = 0.3/2.;
    par[2] = 0.3/2.;
    gMC->Gsposp("UCOA",iCopy+iMCMrow,"TRD1",xpos,ypos,zpos,0,"ONLY",par,kNpar);

  }

}
