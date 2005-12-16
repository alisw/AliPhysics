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
//  Detailed TRD geometry for the spaceframe without holes                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TVirtualMC.h>

#include "AliTRDgeometryDetail.h"
#include "AliTRDCommonParam.h"

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
  // including the MCMs and the cooling pipes
  //
  //
  // Names of the TRD volumina (xx = detector number):
  //
  //      Lower part of the readout chambers (gas volume + radiator)
  //
  //        UAxx    Aluminum frames             (Al)
  //        UBxx    G10 frames                  (C)
  //        UCxx    Inner volumes               (Air)
  //
  //      Upper part of the readout chambers (readout plane + fee)
  //
  //        UDxx    G10 frames                  (C)
  //        UExx    Inner volumes of the G10    (Air)
  //        UFxx    Aluminum frames             (Al)
  //        UGxx    Inner volumes of the Al     (Air)
  //
  //      Inner material layers
  //
  //        UHxx    Radiator                    (Rohacell)
  //        UIxx    Entrance window             (Mylar)
  //        UJxx    Drift volume                (Xe/CO2)
  //        UKxx    Amplification volume        (Xe/CO2)
  //        ULxx    Pad plane                   (Cu)
  //        UMxx    Support structure           (Rohacell)
  //        UNxx    FEE + signal lines          (Cu)
  //

  const Int_t kNparTrd = 4;
  const Int_t kNparCha = 3;

  Float_t xpos, ypos, zpos;

  Float_t parTrd[kNparTrd];
  Float_t parCha[kNparCha];

  Char_t  cTagV[5];
  Char_t  cTagM[5];

  Int_t   idrotm;

  // Rotation matrix
  gMC->Matrix(idrotm,  0.0,  0.0, 90.0, 90.0, 90.0,  0.0);

  // The TRD mother volume for one sector (Air), full length in z-direction
  parTrd[0] = fgkSwidth1/2.;
  parTrd[1] = fgkSwidth2/2.;
  parTrd[2] = fgkSlenTR1/2.;
  parTrd[3] = fgkSheight/2.;
  gMC->Gsvolu("UTR1","TRD1",idtmed[1302-1],parTrd,kNparTrd);

  // Create the readout volumina
  CreateReadout(idtmed);

  // Create the volumina for the cooling
  CreateCooling(idtmed);

  for (Int_t icham = 0; icham < kNcham; icham++) {
    for (Int_t iplan = 0; iplan < kNplan; iplan++) {  

      Int_t iDet = GetDetectorSec(iplan,icham);

      // The lower part of the readout chambers (gas volume + radiator) 
      // The aluminum frames 
      sprintf(cTagV,"UA%02d",iDet);
      parCha[0] = fCwidth[iplan]/2.;
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.;
      parCha[2] = fgkCraH/2. + fgkCdrH/2.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The G10 frames 
      sprintf(cTagV,"UB%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. - fgkCalT; 
      parCha[1] = -1.;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The inner part (air)
      sprintf(cTagV,"UC%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. - fgkCalT - fgkCclsT; 
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.- fgkCclfT;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      // The upper part of the readout chambers (readout plane + fee)
      // The G10 frames
      sprintf(cTagV,"UD%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW;
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.;
      parCha[2] = fgkCamH/2.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1307-1],parCha,kNparCha);
      // The inner part of the G10 frame (air)
      sprintf(cTagV,"UE%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCcuT; 
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.- fgkCcuT;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);
      // The aluminum frames
      sprintf(cTagV,"UF%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW;
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.;
      parCha[2] = fgkCroH/2.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parCha,kNparCha);
      // The inner part of the aluminum frames
      sprintf(cTagV,"UG%02d",iDet);
      parCha[0] = fCwidth[iplan]/2. + fgkCroW - fgkCauT; 
      parCha[1] = fClength[iplan][icham]/2. - fgkHspace/2.- fgkCauT;
      parCha[2] = -1.;
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1302-1],parCha,kNparCha);

      // The material layers inside the chambers
      parCha[0] = -1.;
      parCha[1] = -1.;
      // Rohacell layer (radiator)
      parCha[2] = fgkRaThick/2;
      sprintf(cTagV,"UH%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1315-1],parCha,kNparCha);
      // Mylar layer (entrance window + HV cathode) 
      parCha[2] = fgkMyThick/2;
      sprintf(cTagV,"UI%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1308-1],parCha,kNparCha);
      // Xe/Isobutane layer (drift volume) 
      parCha[2] = fgkDrThick/2.;
      sprintf(cTagV,"UJ%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);
      // Xe/Isobutane layer (amplification volume)
      parCha[2] = fgkAmThick/2.;
      sprintf(cTagV,"UK%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1309-1],parCha,kNparCha);  
      // Cu layer (pad plane)
      parCha[2] = fgkCuThick/2;
      sprintf(cTagV,"UL%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1305-1],parCha,kNparCha);
      // G10 layer (support structure / honeycomb)
      parCha[2] = fgkSuThick/2;
      sprintf(cTagV,"UM%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1313-1],parCha,kNparCha);
      // Cu layer (FEE + signal lines)
      parCha[2] = fgkFeThick/2;
      sprintf(cTagV,"UN%02d",iDet);
      gMC->Gsvolu(cTagV,"BOX ",idtmed[1305-1],parCha,kNparCha);

      // Position the layers in the chambers
      xpos = 0;
      ypos = 0;
      // Lower part
      // Rohacell layer (radiator)
      zpos = fgkRaZpos;
      sprintf(cTagV,"UH%02d",iDet);
      sprintf(cTagM,"UC%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Mylar layer (entrance window + HV cathode)   
      zpos = fgkMyZpos;
      sprintf(cTagV,"UI%02d",iDet);
      sprintf(cTagM,"UC%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Xe/Isobutane layer (drift volume) 
      zpos = fgkDrZpos;
      sprintf(cTagV,"UJ%02d",iDet);
      sprintf(cTagM,"UC%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Upper part
      // Xe/Isobutane layer (amplification volume)
      zpos = fgkAmZpos;
      sprintf(cTagV,"UK%02d",iDet);
      sprintf(cTagM,"UE%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Readout part
      // Cu layer (pad plane)
      zpos = fgkCuZpos; 
      sprintf(cTagV,"UL%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // G10 layer (support structure)
      zpos = fgkSuZpos;
      sprintf(cTagV,"UM%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // Cu layer (FEE + signal lines)
      zpos = fgkFeZpos; 
      sprintf(cTagV,"UN%02d",iDet);
      sprintf(cTagM,"UG%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");

      // Position the inner volumes of the chambers in the frames
      xpos      = 0.0;
      ypos      = 0.0;
      zpos      = 0.0;
      // The inside of the lower G10 frame
      sprintf(cTagV,"UC%02d",iDet);
      sprintf(cTagM,"UB%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The lower G10 frame inside the aluminum frame
      sprintf(cTagV,"UB%02d",iDet);
      sprintf(cTagM,"UA%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The inside of the upper G10 frame
      sprintf(cTagV,"UE%02d",iDet);
      sprintf(cTagM,"UD%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");
      // The inside of the upper aluminum frame
      sprintf(cTagV,"UG%02d",iDet);
      sprintf(cTagM,"UF%02d",iDet);
      gMC->Gspos(cTagV,1,cTagM,xpos,ypos,zpos,0,"ONLY");      

      // Position the frames of the chambers in the TRD mother volume
      xpos  = 0.;
      ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.;
      for (Int_t ic = 0; ic < icham; ic++) {
        ypos += fClength[iplan][ic];        
      }
      ypos += fClength[iplan][icham]/2.;
      zpos  = fgkCraH/2. + fgkCdrH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace);
      // The lower aluminum frame, radiator + drift region
      sprintf(cTagV,"UA%02d",iDet);
      gMC->Gspos(cTagV,1,"UTR1",xpos,ypos,zpos,0,"ONLY");
      // The upper G10 frame, amplification region
      sprintf(cTagV,"UD%02d",iDet);
      zpos += fgkCamH/2. + fgkCraH/2. + fgkCdrH/2.;
      gMC->Gspos(cTagV,1,"UTR1",xpos,ypos,zpos,0,"ONLY");
      // The upper aluminum frame
      sprintf(cTagV,"UF%02d",iDet);
      zpos += fgkCroH/2. + fgkCamH/2.;
      gMC->Gspos(cTagV,1,"UTR1",xpos,ypos,zpos,0,"ONLY");
 
      // Position the MCM volumina
      PositionReadout(iplan,icham);

      // Position the volumina for the cooling
      PositionCooling(iplan,icham,idrotm);

    }
  }

  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTR1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UTR1",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UTR1",3,"BTR3",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::CreateReadout(Int_t *idtmed) const
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

  const Int_t   kNmcmChannel = 18;

  Int_t nMCMrow = AliTRDCommonParam::Instance()->GetRowMax(ipla,icha,0);
  Int_t nMCMcol = AliTRDCommonParam::Instance()->GetColMax(ipla) / kNmcmChannel;

  Float_t xSize = (GetChamberWidth(ipla)       - 2.*fgkCpadW) 
                / ((Float_t) nMCMcol);
  Float_t ySize = (GetChamberLength(ipla,icha) - 2.*fgkRpadW) 
                / ((Float_t) nMCMrow);
  Float_t x0    = AliTRDCommonParam::Instance()->GetCol0(ipla);
  Float_t y0    = AliTRDCommonParam::Instance()->GetRow0(ipla,icha,0);

  Int_t iCopy = GetDetector(ipla,icha,0) * 1000;
  for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {
    for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) {
      iCopy++;
      Float_t xpos = (0.5 + iMCMcol) * xSize + x0; 
      Float_t ypos = (0.5 + iMCMrow) * ySize + y0;
      Float_t zpos = fgkCH - fgkSheight/2. + 0.5/2.
                   + ipla * (fgkCH + fgkVspace);
      gMC->Gspos("UMCM",iCopy,"UTR1",xpos,ypos,zpos,0,"ONLY");    
    }
  }

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::CreateCooling(Int_t *idtmed) const
{
  //
  // Create the volumina of the cooling
  //

  const Int_t kNparTube = 3;

  Float_t parTube[kNparTube];
  Float_t xpos;
  Float_t ypos;
  Float_t zpos;

  // The aluminum pipe for the cooling
  parTube[0] = 0.0;
  parTube[1] = 0.0;
  parTube[2] = 0.0;
  gMC->Gsvolu("UCOA","TUBE",idtmed[1324-1],parTube,0);

  // The cooling water
  parTube[0] =  0.0;
  parTube[1] =  0.2/2.;
  parTube[2] = -1.;
  gMC->Gsvolu("UCOW","TUBE",idtmed[1314-1],parTube,kNparTube);

  // Water inside the cooling pipe
  xpos = 0.0;
  ypos = 0.0;
  zpos = 0.0;
  gMC->Gspos("UCOW",1,"UCOA",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDgeometryDetail::PositionCooling(Int_t ipla, Int_t icha, Int_t idrotm)
{
  //
  // Position the volumina of the cooling
  //

  const Int_t kNpar = 3;

  Float_t par[kNpar];
  Float_t xpos;
  Float_t ypos;
  Float_t zpos;

  Int_t   iCopy   = GetDetector(ipla,icha,0) * 100;
  Int_t   nMCMrow = AliTRDCommonParam::Instance()->GetRowMax(ipla,icha,0);

  Float_t ySize   = (GetChamberLength(ipla,icha) - 2.*fgkRpadW) 
                  / ((Float_t) nMCMrow);
  Float_t y0      = AliTRDCommonParam::Instance()->GetRow0(ipla,icha,0);

  // Position the cooling pipes
  for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) {

    xpos   = 0.0;
    ypos   = (0.5 + iMCMrow) * ySize + y0 - 1.9;
    zpos   = fgkCH - fgkSheight/2. + 0.5/2.
                   + ipla * (fgkCH + fgkVspace);
    par[0] = 0.0;
    par[1] = 0.3/2.;
    par[2] = GetChamberWidth(ipla)/2.+ fgkCroW;
    gMC->Gsposp("UCOA",iCopy+iMCMrow,"UTR1",xpos,ypos,zpos
                      ,idrotm,"ONLY",par,kNpar);

  }
}
