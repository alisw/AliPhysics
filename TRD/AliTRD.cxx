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
Revision 1.13  1999/11/02 16:35:56  fca
New version of TRD introduced

Revision 1.12  1999/11/01 20:41:51  fca
Added protections against using the wrong version of FRAME

Revision 1.11  1999/09/29 09:24:34  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector                                            //
//  This class contains the basic functions for the Transition Radiation     //
//  Detector, as well as the geometry.                                       //
//  Functions specific to one particular geometry are contained in the       // 
//  derived classes.                                                         //
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

  Int_t iplan;

  fIshunt      = 0;
  fGasMix      = 0;
  fHits        = 0;
  fDigits      = 0;
  fHole        = 0;

  fClusters    = 0;
  fNclusters   = 0;

  // The chamber dimensions
  for (iplan = 0; iplan < kNplan; iplan++) {
    fClengthI[iplan]  = 0.;
    fClengthM1[iplan] = 0.;
    fClengthM2[iplan] = 0.;
    fClengthO1[iplan] = 0.;
    fClengthO2[iplan] = 0.;
    fClengthO3[iplan] = 0.;
    fCwidth[iplan]    = 0.;
  }

  for (iplan = 0; iplan < kNplan; iplan++) {
    for (Int_t icham = 0; icham < kNcham; icham++) {
      for (Int_t isect = 0; isect < kNsect; isect++) {
        fRowMax[iplan][icham][isect] = 0;
      }
    }
    fColMax[iplan] = 0;
  }
  fTimeMax       = 0;

  fRowPadSize    = 0;
  fColPadSize    = 0;
  fTimeBinSize   = 0;

}
 
//_____________________________________________________________________________
AliTRD::AliTRD(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for the TRD
  //

  Int_t iplan;

  // Check that FRAME is there otherwise we have no place where to
  // put TRD
  AliModule* FRAME=gAlice->GetModule("FRAME");
  if (!FRAME) {
    Error("Ctor","TRD needs FRAME to be present\n");
    exit(1);
  } 

  // Define the TRD geometry according to the FRAME geometry
  if (FRAME->IsVersion() == 0) 
    // With hole
    fHole = 1;
  else 
    // Without hole
    fHole = 0; 

  // Allocate the hit array
  fHits      = new TClonesArray("AliTRDhit"    ,  405);

  // Allocate the digits array
  fDigits    = new TClonesArray("AliTRDdigit"  ,10000);

  // Allocate the cluster array
  fClusters  = new TClonesArray("AliTRDcluster",  400);
  fNclusters = 0;
   
  fIshunt = 0;
  fGasMix = 0;

  // The chamber dimensions
  for (iplan = 0; iplan < kNplan; iplan++) {
    fClengthI[iplan]  = 0.;
    fClengthM1[iplan] = 0.;
    fClengthM2[iplan] = 0.;
    fClengthO1[iplan] = 0.;
    fClengthO2[iplan] = 0.;
    fClengthO3[iplan] = 0.;
    fCwidth[iplan]    = 0.;
  }
  
  for (iplan = 0; iplan < kNplan; iplan++) {
    for (Int_t icham = 0; icham < kNcham; icham++) {
      for (Int_t isect = 0; isect < kNsect; isect++) {
        fRowMax[iplan][icham][isect] = 0;
      }
    }
    fColMax[iplan] = 0;
  }
  fTimeMax       = 0;

  fRowPadSize    = 0;
  fColPadSize    = 0;
  fTimeBinSize   = 0;

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
  delete fClusters;

}

//_____________________________________________________________________________
void AliTRD::AddCluster(Int_t *tracks, Int_t *clusters, Float_t *position)
{
  //
  // Add a cluster for the TRD
  // 

  TClonesArray &lclusters = *fClusters;
  new(lclusters[fNclusters++]) AliTRDcluster(tracks,clusters,position);

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
  //    TRD1-3     (Air)   --- The TRD mother volumes for one sector. 
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

  Int_t iplan;

  const Int_t npar_trd = 4;
  const Int_t npar_cha = 3;

  Float_t par_dum[3];
  Float_t par_trd[npar_trd];
  Float_t par_cha[npar_cha];

  Float_t xpos, ypos, zpos;

  Int_t *idtmed = fIdtmed->GetArray() - 1299;

  // The length of the inner chambers
  for (iplan = 0; iplan < kNplan; iplan++) 
    fClengthI[iplan] = 110.0;

  // The length of the middle chambers
  fClengthM1[0] = 123.5;
  fClengthM1[1] = 131.0;
  fClengthM1[2] = 138.5;
  fClengthM1[3] = 146.0;
  fClengthM1[4] = 153.0;
  fClengthM1[5] = 160.5;

  fClengthM2[0] = 123.5 - 7.0;
  fClengthM2[1] = 131.0 - 7.0;
  fClengthM2[2] = 138.5 - 7.0;
  fClengthM2[3] = 146.0 - 7.0;
  fClengthM2[4] = 153.0 - 7.0;
  fClengthM2[5] = 160.4 - 7.0;

  // The length of the outer chambers
  fClengthO1[0] = 123.5;
  fClengthO1[1] = 131.0;
  fClengthO1[2] = 134.5;
  fClengthO1[3] = 142.0;
  fClengthO1[4] = 142.0;
  fClengthO1[5] = 134.5;

  fClengthO2[0] = 123.5;
  fClengthO2[1] = 131.0;
  fClengthO2[2] = 134.5;
  fClengthO2[3] = 142.0;
  fClengthO2[4] = 142.0;
  fClengthO2[5] = 134.5;

  fClengthO3[0] =  86.5;
  fClengthO3[1] = 101.5;
  fClengthO3[2] = 112.5;
  fClengthO3[3] = 127.5;
  fClengthO3[4] = 134.5;
  fClengthO3[5] = 134.5;

  // The width of the chambers
  fCwidth[0]    =  99.6;
  fCwidth[1]    = 104.1;
  fCwidth[2]    = 108.5;
  fCwidth[3]    = 112.9;
  fCwidth[4]    = 117.4;
  fCwidth[5]    = 121.8;

  // The TRD mother volume for one sector (Air) (dimensions identical to BTR1)
  par_trd[0] = kSwidth1/2.;
  par_trd[1] = kSwidth2/2.;
  par_trd[2] = kSlenTR1/2.;
  par_trd[3] = kSheight/2.;
  gMC->Gsvolu("TRD1","TRD1",idtmed[1302-1],par_trd,npar_trd);
  
  // The TRD mother volume for one sector (Air) (dimensions identical to BTR2 + BTR3).
  // Only used for the geometry with holes.
  if (fHole) {

    par_trd[0] = kSwidth1/2.;
    par_trd[1] = kSwidth2/2.;
    par_trd[2] = kSlenTR2/2.;
    par_trd[3] = kSheight/2.;
    gMC->Gsvolu("TRD2","TRD1",idtmed[1302-1],par_trd,npar_trd);

    par_trd[0] = kSwidth1/2.;
    par_trd[1] = kSwidth2/2.;
    par_trd[2] = kSlenTR3/2.;
    par_trd[3] = kSheight/2.;
    gMC->Gsvolu("TRD3","TRD1",idtmed[1302-1],par_trd,npar_trd);

  }

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
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthI[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.   - kCathick;
    par_cha[1] = fClengthI[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthI[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFI",iplan       ,"TRD1",xpos,ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.   - kCcthick;
    par_cha[1] = fClengthI[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = 0.;
    zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCII",iplan       ,"TRD1",xpos,ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // The middle chambers --------------------------------------------------------------

    // the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthM1[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2.  + fClengthM1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UAFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCathick;
    par_cha[1] = fClengthM1[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2.  + fClengthM1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UAIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthM1[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFM",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UCFM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCcthick;
    par_cha[1] = fClengthM1[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIM",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UCIM",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // Only for the geometry with holes
    if (fHole) {

      // the aluminum frame
      par_cha[0] = fCwidth[iplan-1]/2.;
      par_cha[1] = fClengthM2[iplan-1]/2.;
      par_cha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAFM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);

      // the inner part of the aluminum frame
      par_cha[0] = fCwidth[iplan-1]/2.    - kCathick;
      par_cha[1] = fClengthM2[iplan-1]/2. - kCathick;
      par_cha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAIM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);

      // the carbon frame
      par_cha[0] = fCwidth[iplan-1]/2.;
      par_cha[1] = fClengthM2[iplan-1]/2.;
      par_cha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCFM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);

      // the inner part of the carbon frame
      par_cha[0] = fCwidth[iplan-1]/2.    - kCcthick;
      par_cha[1] = fClengthM2[iplan-1]/2. - kCcthick;
      par_cha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCIM",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
      
    }

    // The outer chambers ---------------------------------------------------------------

    // the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthO1[iplan-1]/2.;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAFO",iplan         ,"TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UAFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the aluminum frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCathick;
    par_cha[1] = fClengthO1[iplan-1]/2. - kCathick;
    par_cha[2] = kCaframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UAIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UAIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.;
    par_cha[1] = fClengthO1[iplan-1]/2.;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCFO",iplan,         "TRD1",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);
    gMC->Gsposp("UCFO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"MANY",par_cha,npar_cha);

    // the inner part of the carbon frame
    par_cha[0] = fCwidth[iplan-1]/2.    - kCcthick;
    par_cha[1] = fClengthO1[iplan-1]/2. - kCcthick;
    par_cha[2] = kCcframe/2.;
    xpos       = 0.;
    ypos       = fClengthI[iplan-1]/2. + fClengthM1[iplan-1]    + fClengthO1[iplan-1]/2.;
    zpos       = kCcframe/2.           - kSheight/2. + (iplan-1) * (kCheight + kCspace);
    gMC->Gsposp("UCIO",iplan         ,"TRD1",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);
    gMC->Gsposp("UCIO",iplan+  kNplan,"TRD1",xpos,-ypos,zpos,0,"ONLY",par_cha,npar_cha);

    // Only for the geometry with holes
    if (fHole) {

      // the aluminum frame
      par_cha[0] = fCwidth[iplan-1]/2.;
      par_cha[1] = fClengthO2[iplan-1]/2.;
      par_cha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]    + fClengthO2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAFO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);

      // the inner part of the aluminum frame
      par_cha[0] = fCwidth[iplan-1]/2.    - kCathick;
      par_cha[1] = fClengthO2[iplan-1]/2. - kCathick;
      par_cha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]    + fClengthO2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAIO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);

      // the carbon frame
      par_cha[0] = fCwidth[iplan-1]/2.;
      par_cha[1] = fClengthO2[iplan-1]/2.;
      par_cha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]    + fClengthO2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCFO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);

      // the inner part of the carbon frame
      par_cha[0] = fCwidth[iplan-1]/2.    - kCcthick;
      par_cha[1] = fClengthO2[iplan-1]/2. - kCcthick;
      par_cha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthM2[iplan-1]    + fClengthO2[iplan-1]/2. - kSlenTR2/2.;
      zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCIO",iplan+2*kNplan,"TRD2",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);

      // the aluminum frame
      par_cha[0] = fCwidth[iplan-1]/2.;
      par_cha[1] = fClengthO3[iplan-1]/2.;
      par_cha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthO3[iplan-1]/2. - kSlenTR3/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAFO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);

      // the inner part of the aluminum frame
      par_cha[0] = fCwidth[iplan-1]/2.    - kCathick;
      par_cha[1] = fClengthO3[iplan-1]/2. - kCathick;
      par_cha[2] = kCaframe/2.;
      xpos       = 0.;
      ypos       = fClengthO3[iplan-1]/2. - kSlenTR3/2.;
      zpos       = kCheight - kCaframe/2. - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UAIO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);

      // the carbon frame
      par_cha[0] = fCwidth[iplan-1]/2.;
      par_cha[1] = fClengthO3[iplan-1]/2.;
      par_cha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthO3[iplan-1]/2. - kSlenTR3/2.;
      zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCFO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"MANY",par_cha,npar_cha);

      // the inner part of the carbon frame
      par_cha[0] = fCwidth[iplan-1]/2.    - kCcthick;
      par_cha[1] = fClengthO3[iplan-1]/2. - kCcthick;
      par_cha[2] = kCcframe/2.;
      xpos       = 0.;
      ypos       = fClengthO3[iplan-1]/2. - kSlenTR3/2.;
      zpos       = kCcframe/2.            - kSheight/2. + (iplan-1) * (kCheight + kCspace);
      gMC->Gsposp("UCIO",iplan+4*kNplan,"TRD3",xpos, ypos,zpos,0,"ONLY",par_cha,npar_cha);

    }

  }

  if (fHole) {
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gspos("TRD1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
    gMC->Gspos("TRD2",1,"BTR2",xpos,ypos,zpos,0,"ONLY");
    gMC->Gspos("TRD3",1,"BTR3",xpos,ypos,zpos,0,"ONLY");
  }
  else {
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gspos("TRD1",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
    gMC->Gspos("TRD1",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
    gMC->Gspos("TRD1",3,"BTR3",xpos,ypos,zpos,0,"ONLY");
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
  if (fHole) {
    gMC->Gsatt("B071","SEEN", 0);
    gMC->Gsatt("B074","SEEN", 0);
    gMC->Gsatt("B075","SEEN", 0);
    gMC->Gsatt("B077","SEEN", 0);
    gMC->Gsatt("BTR1","SEEN", 0);
    gMC->Gsatt("BTR2","SEEN", 0);
    gMC->Gsatt("BTR3","SEEN", 0);
    gMC->Gsatt("TRD1","SEEN", 0);
    gMC->Gsatt("TRD2","SEEN", 0);
    gMC->Gsatt("TRD3","SEEN", 0);
  }
  else {
    gMC->Gsatt("B071","SEEN", 0);
    gMC->Gsatt("B074","SEEN", 0);
    gMC->Gsatt("B075","SEEN", 0);
    gMC->Gsatt("B077","SEEN", 0);
    gMC->Gsatt("BTR1","SEEN", 0);
    gMC->Gsatt("BTR2","SEEN", 0);
    gMC->Gsatt("BTR3","SEEN", 0);
    gMC->Gsatt("TRD1","SEEN", 0);
  }
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
  Int_t iplan;
  
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

  if (fHole)
    printf("          Geometry with holes\n");
  else
    printf("          Full geometry\n");

  // The default pad dimensions
  if (!(fRowPadSize))  fRowPadSize  = 4.5;
  if (!(fColPadSize))  fColPadSize  = 1.0;
  if (!(fTimeBinSize)) fTimeBinSize = 0.1;

  // The maximum number of pads
  // and the position of pad 0,0,0 
  // 
  // chambers seen from the top:
  //     +----------------------------+
  //     |                            |
  //     |                            |     ^
  //     |                            | rphi|
  //     |                            |     |
  //     |0                           |     | 
  //     +----------------------------+     +------>
  //                                             z 
  // chambers seen from the side:           ^
  //     +----------------------------+ time|
  //     |                            |     |
  //     |0                           |     |
  //     +----------------------------+     +------>
  //                                             z
  //                                             
  for (iplan = 0; iplan < kNplan; iplan++) {

    // The pad row (z-direction)
    for (Int_t isect = 0; isect < kNsect; isect++) {
      Float_t clengthI = fClengthI[iplan];
      Float_t clengthM = fClengthM1[iplan];
      Float_t clengthO = fClengthO1[iplan];
      if (fHole) {
        switch (isect) {
        case 12:
	case 13:
	case 14:
	case 15:
	case 16:
          clengthM = fClengthM2[iplan];
          clengthO = fClengthO2[iplan];
          break;
        case 4:
	case 5:
	case 6:
          clengthO = fClengthO3[iplan];
          break;
        };
      }
      fRowMax[iplan][0][isect] = 1 + TMath::Nint((clengthO - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][1][isect] = 1 + TMath::Nint((clengthM - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][2][isect] = 1 + TMath::Nint((clengthI - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][3][isect] = 1 + TMath::Nint((clengthM - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRowMax[iplan][4][isect] = 1 + TMath::Nint((clengthO - 2. * kCcthick) 
                                                           / fRowPadSize - 0.5);
      fRow0[iplan][0][isect]   = -clengthI/2. - clengthM - clengthO + kCcthick; 
      fRow0[iplan][1][isect]   = -clengthI/2. - clengthM            + kCcthick;
      fRow0[iplan][2][isect]   = -clengthI/2.                       + kCcthick;
      fRow0[iplan][3][isect]   =  clengthI/2.                       + kCcthick; 
      fRow0[iplan][4][isect]   =  clengthI/2. + clengthM            + kCcthick; 
    }

    // The pad column (rphi-direction)  
    fColMax[iplan]    = 1 + TMath::Nint((fCwidth[iplan] - 2. * kCcthick) 
                                                        / fColPadSize - 0.5);
    fCol0[iplan]      = -fCwidth[iplan]/2. + kCcthick;

  }

  // The time bucket
  fTimeMax = 1 + TMath::Nint(kDrThick / fTimeBinSize - 0.5);
  for (Int_t iplan = 0; iplan < kNplan; iplan++) {
    fTime0[iplan]   = kRmin + kCcframe/2. + kDrZpos - 0.5 * kDrThick
                            + iplan * (kCheight + kCspace);
  } 

}

//_____________________________________________________________________________
void AliTRD::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the TRD digits and cluster.
  //

  Int_t  buffersize = 4000;
  Char_t branchname[15];

  AliDetector::MakeBranch(option);

  Char_t *D = strstr(option,"D");
  sprintf(branchname,"%s",GetName());
  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits,  buffersize);
    printf("* AliTRD::MakeBranch * Making Branch %s for digits in TreeD\n",branchname);
  }

  sprintf(branchname,"%scluster",GetName());
  if (fClusters && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fClusters,buffersize);
    printf("* AliTRD::MakeBranch * Making Branch %s for cluster in TreeD\n",branchname);
  }

}

//_____________________________________________________________________________
void AliTRD::SetTreeAddress()
{
  //
  // Set the branch addresses for the trees.
  //

  Char_t branchname[15];

  AliDetector::SetTreeAddress();

  TBranch *branch;
  TTree   *treeD = gAlice->TreeD();

  if (treeD) {
    sprintf(branchname,"%scluster",GetName());    
    if (fClusters) {
      branch = treeD->GetBranch(branchname);
      if (branch) branch->SetAddress(&fClusters);
    }
  }

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

//______________________________________________________________________________
void AliTRD::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliTRD.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDetector::Streamer(R__b);
      R__b >> fGasMix;
      R__b.ReadStaticArray(fClengthI);
      R__b.ReadStaticArray(fClengthM1);
      R__b.ReadStaticArray(fClengthM2);
      R__b.ReadStaticArray(fClengthO1);
      R__b.ReadStaticArray(fClengthO2);
      R__b.ReadStaticArray(fClengthO3);
      R__b.ReadStaticArray(fCwidth);
      R__b.ReadStaticArray((int*)fRowMax);
      R__b.ReadStaticArray(fColMax);
      R__b >> fTimeMax;
      R__b.ReadStaticArray((float*)fRow0);
      R__b.ReadStaticArray(fCol0);
      R__b.ReadStaticArray(fTime0);
      R__b >> fRowPadSize;
      R__b >> fColPadSize;
      R__b >> fTimeBinSize;
      R__b >> fHole;
      // Stream the pointers but not the TClonesArray
      R__b >> fClusters;     // diff
      //R__b >> fNclusters;
   } else {
      R__b.WriteVersion(AliTRD::IsA());
      AliDetector::Streamer(R__b);
      R__b << fGasMix;
      R__b.WriteArray(fClengthI, 6);
      R__b.WriteArray(fClengthM1, 6);
      R__b.WriteArray(fClengthM2, 6);
      R__b.WriteArray(fClengthO1, 6);
      R__b.WriteArray(fClengthO2, 6);
      R__b.WriteArray(fClengthO3, 6);
      R__b.WriteArray(fCwidth, 6);
      R__b.WriteArray((int*)fRowMax, 540);
      R__b.WriteArray(fColMax, 6);
      R__b << fTimeMax;
      R__b.WriteArray((float*)fRow0, 540);
      R__b.WriteArray(fCol0, 6);
      R__b.WriteArray(fTime0, 6);
      R__b << fRowPadSize;
      R__b << fColPadSize;
      R__b << fTimeBinSize;
      R__b << fHole;
      // Stream the pointers but not the TClonesArrays
      R__b << fClusters;     // diff
      //R__b << fNclusters;
   }

}

ClassImp(AliTRDhit)
 
//_____________________________________________________________________________
AliTRDhit::AliTRDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
          :AliHit(shunt, track)
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

ClassImp(AliTRDcluster)

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(Int_t *tracks, Int_t *cluster, Float_t* position)
              :TObject()
{
  //
  // Create a TRD cluster
  //

  fSector    = cluster[0];
  fChamber   = cluster[1];
  fPlane     = cluster[2];

  fTimeSlice = cluster[3];
  fEnergy    = cluster[4];

  fX         = position[0];
  fY         = position[1];
  fZ         = position[2];

  fTracks[0] = tracks[0];
  fTracks[1] = tracks[1];
  fTracks[2] = tracks[2];

}

