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
//  Transition Radiation Detector                                            //
//  This class contains the basic functions for the Transition Radiation     //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TMath.h>
#include <TNode.h>
#include <TPGON.h> 

#include "AliTRD.h"
#include "AliRun.h"

#include "AliConst.h"
 
ClassImp(AliTRD)
 
//_____________________________________________________________________________
AliTRD::AliTRD()
{
  //
  // Default constructor
  //

  fIshunt      = 0;
  fGasMix      = 0;
  fHits        = 0;
  fDigits      = 0;

  // The chamber dimensions
  for (Int_t iplan = 0; iplan < kNplan; iplan++) {
    fClengthI[iplan] = 0.;
    fClengthM[iplan] = 0.;
    fClengthO[iplan] = 0.;
  }

}
 
//_____________________________________________________________________________
AliTRD::AliTRD(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for the TRD
  //


  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) {
    Error("AliTRD","TRD needs FRAME to be present\n");
    exit(1);
  }

  // Allocate the hit array
  fHits   = new TClonesArray("AliTRDhit",  405);

  // Allocate the digits array
  fDigits = new TClonesArray("AliTRDdigit",10000);
   
  fIshunt = 0;
  fGasMix = 0;

  // The chamber dimensions
  for (Int_t iplan = 0; iplan < kNplan; iplan++) {
    fClengthI[iplan] = 0.;
    fClengthM[iplan] = 0.;
    fClengthO[iplan] = 0.;
    fCwidth[iplan]   = 0.;
  }
  
  SetMarkerColor(kWhite);   

}

//_____________________________________________________________________________
AliTRD::~AliTRD()
{
  //
  // TRD destructor
  //

  fIshunt = 0;

  delete fHits;
  delete fDigits;

}

//_____________________________________________________________________________
void AliTRD::AddDigit(Int_t *tracks, Int_t *digits)
{
  //
  // Add a digit for the TRD
  //

  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliTRDdigit(tracks,digits);

}

//_____________________________________________________________________________
void AliTRD::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a hit for the TRD
  //

  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTRDhit(fIshunt,track,vol,hits);

}

//_____________________________________________________________________________
void AliTRD::BuildGeometry()
{
  //
  // Create the ROOT TNode geometry for the TRD
  //

  TNode *Node, *Top;
  TPGON *pgon;
  const Int_t kColorTRD = 46;
  
  // Find the top node alice
  Top = gAlice->GetGeometry()->GetNode("alice");
  
  pgon = new TPGON("S_TRD","TRD","void",0,360,kNsect,4);
  Float_t ff    = TMath::Cos(kDegrad * 180 / kNsect);
  Float_t rrmin = kRmin / ff;
  Float_t rrmax = kRmax / ff;
  pgon->DefineSection(0,-kZmax1,rrmax,rrmax);
  pgon->DefineSection(1,-kZmax2,rrmin,rrmax);
  pgon->DefineSection(2, kZmax2,rrmin,rrmax);
  pgon->DefineSection(3, kZmax1,rrmax,rrmax);
  Top->cd();
  Node = new TNode("TRD","TRD","S_TRD",0,0,0,"");
  Node->SetLineColor(kColorTRD);
  fNodes->Add(Node);

}
 
//_____________________________________________________________________________
void AliTRD::CreateGeometry()
{
  //
  // Creates the volumes for the TRD chambers
  //
  // Author: Christoph Blume (C.Blume@gsi.de) 20/07/99
  //
  // The volumes:
  //    TRD        (Air)   --- The TRD mother volume for one sector. 
  //                           To be placed into the spaceframe.
  //
  //    UAFI(/M/O) (Al)    --- The aluminum frame of the inner(/middle/outer) chambers (readout)
  //    UCFI(/M/O) (C)     --- The carbon frame of the inner(/middle/outer) chambers 
  //                           (driftchamber + radiator)
  //    UAII(/M/O) (Air)   --- The inner part of the readout of the inner(/middle/outer) chambers
  //    UFII(/M/O) (Air)   --- The inner part of the chamner and radiator of the 
  //                           inner(/middle/outer) chambers
  //
  // The material layers in one chamber:
  //    UL01       (G10)   --- The gas seal of the radiator
  //    UL02       (CO2)   --- The gas in the radiator
  //    UL03       (PE)    --- The foil stack
  //    UL04       (Mylar) --- Entrance window to the driftvolume and HV-cathode
  //    UL05       (Xe)    --- The driftvolume
  //    UL06       (Xe)    --- The amplification region
  //    
  //    UL07       (Cu)    --- The pad plane
  //    UL08       (G10)   --- The Nomex honeycomb support structure
  //    UL09       (Cu)    --- FEE and signal lines
  //    UL10       (PE)    --- The cooling devices
  //    UL11       (Water) --- The cooling water

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) return;

  const Int_t npar_trd = 4;
  const Int_t npar_cha = 3;

  Float_t par_dum[3];
  Float_t par_trd[npar_trd];
  Float_t par_cha[npar_cha];
  Int_t iplan;

  Float_t xpos, ypos, zpos;

  Int_t *idtmed = fIdtmed->GetArray()-1299;

  // The length of the inner chambers
  for (iplan = 0; iplan < kNplan; iplan++) fClengthI[iplan] = 110.0;
  // The length of the middle chambers
  fClengthM[0] = 123.5;
  fClengthM[1] = 131.0;
  fClengthM[2] = 138.5;
  fClengthM[3] = 146.0;
  fClengthM[4] = 153.0;
  fClengthM[5] = 160.5;
  // The length of the outer chambers
  fClengthO[0] = 123.5;
  fClengthO[1] = 131.0;
  fClengthO[2] = 134.5;
  fClengthO[3] = 142.0;
  fClengthO[4] = 142.0;
  fClengthO[5] = 134.5;

  // The width of the chambers
  fCwidth[0] =  99.6;
  fCwidth[1] = 104.1;
  fCwidth[2] = 108.5;
  fCwidth[3] = 112.9;
  fCwidth[4] = 117.4;
  fCwidth[5] = 121.8;

  // The TRD mother volume for one sector (Air) (dimensions identical to BTR1-3)
  par_trd[0] = kSwidth1/2.;
  par_trd[1] = kSwidth2/2.;
  par_trd[2] = kSlength/2.;
  par_trd[3] = kSheight/2.;
  gMC->Gsvolu("TRD ","TRD1",idtmed[1302-1],par_trd,npar_trd);

  // The aluminum frames - readout + electronics (Al)
  // The inner chambers
  gMC->Gsvolu("UAFI","BOX ",idtmed[1301-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UAFM","BOX ",idtmed[1301-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UAFO","BOX ",idtmed[1301-1],par_dum,0);

  // The inner part of the aluminum frames (Air)
  // The inner chambers
  gMC->Gsvolu("UAII","BOX ",idtmed[1302-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UAIM","BOX ",idtmed[1302-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UAIO","BOX ",idtmed[1302-1],par_dum,0);

  // The carbon frames - radiator + driftchamber (C)
  // The inner chambers
  gMC->Gsvolu("UCFI","BOX ",idtmed[1307-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UCFM","BOX ",idtmed[1307-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UCFO","BOX ",idtmed[1307-1],par_dum,0);

  // The inner part of the carbon frames (Air)
  // The inner chambers
  gMC->Gsvolu("UCII","BOX ",idtmed[1302-1],par_dum,0);
  // The middle chambers
  gMC->Gsvolu("UCIM","BOX ",idtmed[1302-1],par_dum,0);
  // The outer chambers
  gMC->Gsvolu("UCIO","BOX ",idtmed[1302-1],par_dum,0);

  // The material layers inside the chambers
  par_cha[0] = -1.;
  par_cha[1] = -1.;
  // G10 layer (radiator seal)
  par_cha[2] = kSeThick/2;
  gMC->Gsvolu("UL01","BOX ",idtmed[1313-1],par_cha,npar_cha);
  // CO2 layer (radiator)
  par_cha[2] = kRaThick/2;
  gMC->Gsvolu("UL02","BOX ",idtmed[1312-1],par_cha,npar_cha);
  // PE layer (radiator)
  par_cha[2] = kPeThick/2;
  gMC->Gsvolu("UL03","BOX ",idtmed[1303-1],par_cha,npar_cha);
  // Mylar layer (entrance window + HV cathode) 
  par_cha[2] = kMyThick/2;
  gMC->Gsvolu("UL04","BOX ",idtmed[1308-1],par_cha,npar_cha);
  // Xe/Isobutane layer (drift volume, sensitive) 
  par_cha[2] = kDrThick/2.;
  gMC->Gsvolu("UL05","BOX ",idtmed[1309-1],par_cha,npar_cha);
  // Xe/Isobutane layer (amplification volume, not sensitive)
  par_cha[2] = kAmThick/2.;
  gMC->Gsvolu("UL06","BOX ",idtmed[1309-1],par_cha,npar_cha);
  
  // Cu layer (pad plane)
  par_cha[2] = kCuThick/2;
  gMC->Gsvolu("UL07","BOX ",idtmed[1305-1],par_cha,npar_cha);
  // G10 layer (support structure)
  par_cha[2] = kSuThick/2;
  gMC->Gsvolu("UL08","BOX ",idtmed[1313-1],par_cha,npar_cha);
  // Cu layer (FEE + signal lines)
  par_cha[2] = kFeThick/2;
  gMC->Gsvolu("UL09","BOX ",idtmed[1305-1],par_cha,npar_cha);
  // PE layer (cooling devices)
  par_cha[2] = kCoThick/2;
  gMC->Gsvolu("UL10","BOX ",idtmed[1303-1],par_cha,npar_cha);
  // Water layer (cooling)
  par_cha[2] = kWaThick/2;
  gMC->Gsvolu("UL11","BOX ",idtmed[1314-1],par_cha,npar_cha);

  // Position the layers in the chambers
  xpos = 0;
  ypos = 0;

  // G10 layer (radiator seal)
  zpos = kSeZpos;
  gMC->Gspos("UL01",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL01",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL01",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // CO2 layer (radiator)
  zpos = kRaZpos;
  gMC->Gspos("UL02",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL02",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL02",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // PE layer (radiator)
  zpos = 0;
  gMC->Gspos("UL03",1,"UL02",xpos,ypos,zpos,0,"ONLY");
  // Mylar layer (entrance window + HV cathode)   
  zpos = kMyZpos;
  gMC->Gspos("UL04",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL04",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (drift volume) 
  zpos = kDrZpos;
  gMC->Gspos("UL05",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL05",3,"UCIO",xpos,ypos,zpos,0,"ONLY");
  // Xe/Isobutane layer (amplification volume)
  zpos = kAmZpos;
  gMC->Gspos("UL06",1,"UCII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",2,"UCIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL06",3,"UCIO",xpos,ypos,zpos,0,"ONLY");

  // Cu layer (pad plane)
  zpos = kCuZpos;
  gMC->Gspos("UL07",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL07",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // G10 layer (support structure)
  zpos = kSuZpos;
  gMC->Gspos("UL08",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL08",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // Cu layer (FEE + signal lines)
  zpos = kFeZpos; 
  gMC->Gspos("UL09",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL09",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL09",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // PE layer (cooling devices)
  zpos = kCoZpos;
  gMC->Gspos("UL10",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL10",2,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL10",3,"UAIO",xpos,ypos,zpos,0,"ONLY");
  // Water layer (cooling)
  zpos = kWaZpos;
  gMC->Gspos("UL11",1,"UAII",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL11",1,"UAIM",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("UL11",1,"UAIO",xpos,ypos,zpos,0,"ONLY");

  // Position the chambers in the TRD mother volume
  for (iplan = 1; iplan <= kNplan; iplan++) {

    // The inner chambers ---------------------------------------------------------------

    // the aluminum frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2.;
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthI[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFI",iplan       ,"TRD ",xpos,ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2. - kCathick;
    par_cha[0] = fCwidth[iplan-1]/2.   - kCathick;
    par_cha[1] = fClengthI[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAII",iplan       ,"TRD ",xpos,ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2.;
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthI[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFI",iplan       ,"TRD ",xpos,ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2. - kCcthick;
    par_cha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    par_cha[1] = fClengthI[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCII",iplan       ,"TRD ",xpos,ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // The middle chambers --------------------------------------------------------------

    // the aluminum frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2.;
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthM[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFM",iplan       ,"TRD ",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UAFM",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2. - kCathick;
    par_cha[0] = fCwidth[iplan-1]/2.   - kCathick;
    par_cha[1] = fClengthM[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIM",iplan       ,"TRD ",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UAIM",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2.;
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthM[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]/2.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFM",iplan,       "TRD ",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UCFM",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2. - kCcthick;
    par_cha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    par_cha[1] = fClengthM[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]/2.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIM",iplan       ,"TRD ",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UCIM",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // The outer chambers ---------------------------------------------------------------

    // the aluminum frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2.;
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthO[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFO",iplan       ,"TRD ",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UAFO",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2. - kCathick;
    par_cha[0] = fCwidth[iplan-1]/2.   - kCathick;
    par_cha[1] = fClengthO[iplan-1]/2.  - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIO",iplan       ,"TRD ",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UAIO",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2.;
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthO[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFO",iplan,       "TRD ",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UCFO",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    //par_cha[0] = kSwidth1/2. + (iplan-1) * kCwidcha/2. - kCcthick;
    par_cha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    par_cha[1] = fClengthO[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM[iplan-1]    + fClengthO[iplan-1]/2.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIO",iplan       ,"TRD ",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UCIO",iplan+kNplan,"TRD ",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

  }

}
 
//_____________________________________________________________________________
void AliTRD::CreateMaterials()
{
  //
  // Create the materials for the TRD
  // Origin Y.Foka
  //

  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  // For polyethilene (CH2) 
  Float_t ape[2] = { 12., 1. };
  Float_t zpe[2] = {  6., 1. };
  Float_t wpe[2] = {  1., 2. };
  Float_t dpe    = 0.95;

  // For mylar (C5H4O2) 
  Float_t amy[3] = { 12., 1., 16. };
  Float_t zmy[3] = {  6., 1.,  8. };
  Float_t wmy[3] = {  5., 4.,  2. };
  Float_t dmy    = 1.39;

  // For CO2 
  Float_t aco[2] = { 12., 16. };
  Float_t zco[2] = {  6.,  8. };
  Float_t wco[2] = {  1.,  2. };
  Float_t dco    = 0.001977;

  // For water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;

  // For isobutane (C4H10)
  Float_t ais[2] = { 12.,  1. };
  Float_t zis[2] = {  6.,  1. };
  Float_t wis[2] = {  4., 10. };
  Float_t dis    = 0.00267;

  // For Xe/CO2-gas-mixture 
  // Xe-content of the Xe/CO2-mixture (90% / 10%) 
  Float_t fxc    = .90;
  // Xe-content of the Xe/Isobutane-mixture (97% / 3%) 
  Float_t fxi    = .97;
  Float_t dxe    = .005858;
  
  // General tracking parameter
  Float_t tmaxfd = -10.;
  Float_t stemax = -1e10;
  Float_t deemax = -0.1;
  Float_t epsil  =  1e-4;
  Float_t stmin  = -0.001;
  
  Float_t absl, radl, d, buf[1];
  Float_t agm[2], dgm, zgm[2], wgm[2];
  Int_t   nbuf;
  
  //////////////////////////////////////////////////////////////////////////
  //     Define Materials 
  //////////////////////////////////////////////////////////////////////////

  AliMaterial( 1, "Al $",  26.98, 13.0, 2.7     ,     8.9 ,    37.2);
  AliMaterial( 2, "Air$",  14.61,  7.3, 0.001205, 30420.0 , 67500.0);
  AliMaterial( 4, "Xe $", 131.29, 54.0, dxe     ,  1447.59,     0.0);
  AliMaterial( 5, "Cu $",  63.54, 29.0, 8.96    ,     1.43,    14.8);
  AliMaterial( 6, "C  $",  12.01,  6.0, 2.265   ,    18.8 ,    74.4);
  AliMaterial(12, "G10$",  20.00, 10.0, 1.7     ,    19.4 ,   999.0);

  // Mixtures 
  AliMixture(3, "Polyethilene$",   ape, zpe, dpe, -2, wpe);
  AliMixture(7, "Mylar$",          amy, zmy, dmy, -3, wmy);
  AliMixture(8, "CO2$",            aco, zco, dco, -2, wco);
  AliMixture(9, "Isobutane$",      ais, zis, dis, -2, wis);
  AliMixture(13,"Water$",          awa, zwa, dwa, -2, wwa);

  // Gas mixtures
  Char_t namate[21];
  // Xe/CO2-mixture
  // Get properties of Xe 
  gMC->Gfmate((*fIdmate)[4], namate, agm[0], zgm[0], d, radl, absl, buf, nbuf);
  // Get properties of CO2 
  gMC->Gfmate((*fIdmate)[8], namate, agm[1], zgm[1], d, radl, absl, buf, nbuf);
  // Create gas mixture 
  wgm[0] = fxc;
  wgm[1] = 1. - fxc;
  dgm    = wgm[0] * dxe + wgm[1] * dco;
  AliMixture(10, "Gas mixture 1$", agm, zgm, dgm,  2, wgm);
  // Xe/Isobutane-mixture
  // Get properties of Xe 
  gMC->Gfmate((*fIdmate)[4], namate, agm[0], zgm[0], d, radl, absl, buf, nbuf);
  // Get properties of Isobutane
  gMC->Gfmate((*fIdmate)[9], namate, agm[1], zgm[1], d, radl, absl, buf, nbuf);
  // Create gas mixture 
  wgm[0] = fxi;
  wgm[1] = 1. - fxi;
  dgm    = wgm[0] * dxe + wgm[1] * dis;
  AliMixture(11, "Gas mixture 2$", agm, zgm, dgm,  2, wgm);
 
  //////////////////////////////////////////////////////////////////////////
  //     Tracking Media Parameters 
  //////////////////////////////////////////////////////////////////////////

  // Al Frame 
  AliMedium(1, "Al Frame$",   1, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Air 
  AliMedium(2, "Air$",        2, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Polyethilene 
  AliMedium(3, "Radiator$",   3, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Xe 
  AliMedium(4, "Xe$",         4, 1, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Cu pads 
  AliMedium(5, "Padplane$",   5, 1, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Fee + cables 
  AliMedium(6, "Readout$",    1, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // C frame 
  AliMedium(7, "C Frame$",    6, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Mylar foils 
  AliMedium(8, "Mylar$",      7, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  if (fGasMix == 1) {
    // Gas-mixture (Xe/CO2) 
    AliMedium(9, "Gas-mix$",   10, 1, ISXFLD, SXMGMX
                  , tmaxfd, stemax, deemax, epsil, stmin);
  }
  else {
    // Gas-mixture (Xe/Isobutane) 
    AliMedium(9, "Gas-mix$",   11, 1, ISXFLD, SXMGMX
                  , tmaxfd, stemax, deemax, epsil, stmin);
  }
  // Nomex-honeycomb (use carbon for the time being) 
  AliMedium(10, "Nomex$",      6, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Kapton foils (use Mylar for the time being) 
  AliMedium(11, "Kapton$",     7, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Gas-filling of the radiator 
  AliMedium(12, "CO2$",        8, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // G10-plates
  AliMedium(13, "G10-plates$",12, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Cooling water
  AliMedium(14, "Water$",     13, 0, ISXFLD, SXMGMX
                , tmaxfd, stemax, deemax, epsil, stmin);

}

//_____________________________________________________________________________
void AliTRD::DrawModule()
{
  //
  // Draw a shaded view of the Transition Radiation Detector version 0
  //

  // Set everything unseen
  gMC->Gsatt("*"   ,"SEEN",-1);
  
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN", 0);
  
  // Set the volumes visible
  gMC->Gsatt("B032","SEEN", 0);
  gMC->Gsatt("B028","SEEN", 0);
  gMC->Gsatt("B029","SEEN", 0);
  gMC->Gsatt("B030","SEEN", 0);
  gMC->Gsatt("BTR1","SEEN", 0);
  gMC->Gsatt("BTR2","SEEN", 0);
  gMC->Gsatt("BTR3","SEEN", 0);
  gMC->Gsatt("TRD" ,"SEEN", 0);
  gMC->Gsatt("UCII","SEEN", 0);
  gMC->Gsatt("UCIM","SEEN", 0);
  gMC->Gsatt("UCIO","SEEN", 0);
  gMC->Gsatt("UL02","SEEN", 1);
  gMC->Gsatt("UL05","SEEN", 1);
  gMC->Gsatt("UL06","SEEN", 1);
  
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.4, .021, .021);
  gMC->Gdhead(1111, "Transition Radiation Detector");
  gMC->Gdman(18, 4, "MAN");

}

//_____________________________________________________________________________
Int_t AliTRD::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance between the mouse and the TRD detector on the screen
  // Dummy routine
  
  return 9999;

}
 
//_____________________________________________________________________________
void AliTRD::Init()
{
  //
  // Initialise the TRD detector after the geometry has been created
  //

  Int_t i;
  
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" TRD_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  
  // Here the TRD initialisation code (if any!)
  if (fGasMix == 1) 
    printf("          Gas Mixture: 90%% Xe + 10%% CO2\n");
  else
    printf("          Gas Mixture: 97%% Xe + 3%% Isobutane\n");

}

//_____________________________________________________________________________
void AliTRD::SetGasMix(Int_t imix)
{
  //
  // Defines the gas mixture (imix=0:  Xe/Isobutane imix=1: Xe/CO2)
  //
  
  if ((imix < 0) || (imix > 1)) {
    printf("Wrong input value: %d\n",imix);
    printf("Use standard setting\n");
    fGasMix = 0;
    return;
  }

  fGasMix = imix;

}

ClassImp(AliTRDhit)
 
//_____________________________________________________________________________
AliTRDhit::AliTRDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Create a TRD hit
  //

  // Store volume hierarchy
  fSector  = vol[0]; 
  fChamber = vol[1];
  fPlane   = vol[2];
  
  // Store position and charge
  fX       = hits[0];
  fY       = hits[1];
  fZ       = hits[2];
  fQ       = hits[3];

}

ClassImp(AliTRDdigit)

//_____________________________________________________________________________
AliTRDdigit::AliTRDdigit(Int_t *tracks, Int_t *digits)
            :AliDigit(tracks)
{
  //
  // Create a TRD digit
  //

  // Store the volume hierarchy
  fSector  = digits[0];
  fChamber = digits[1];
  fPlane   = digits[2];

  // Store the row, pad, and time bucket number
  fRow     = digits[3];
  fCol     = digits[4];
  fTime    = digits[5];

  // Store the signal amplitude
  fSignal  = digits[6];

}
