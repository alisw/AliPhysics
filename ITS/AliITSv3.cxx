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
Revision 1.8  1999/10/05 08:05:09  fca
Minor corrections for uninitialised variables.

Revision 1.7  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System version 3                                           //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
// Authors: R. Barbera, A. Morsch.
// version 3.
// Created  1998.
//
//  NOTE: THIS IS THE OLD detailed TP-like geometry of the ITS. THIS WILL NOT 
// WORK with the geometry or module classes or any analysis classes. You are 
// strongly encouraged to uses AliITSv5.
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include "AliITShit.h"
#include "AliITSv3.h"
#include "AliRun.h"

#include "AliMC.h"
#include "AliConst.h"

ClassImp(AliITSv3)
 
//_____________________________________________________________________________
AliITSv3::AliITSv3() {
    //
    // Default constructor for ITS
    //
    fId3N = 6;
    fId3Name = new char*[fId3N];
    fId3Name[0] = "ITS1";
    fId3Name[1] = "ITS2";
    fId3Name[2] = "ITS3";
    fId3Name[3] = "ITS4";
    fId3Name[4] = "ITS5";
    fId3Name[5] = "ITS6";
    fMinorVersionV3=1;
}
 
//_____________________________________________________________________________
AliITSv3::~AliITSv3() {
    //
    // Standard destructor for the ITS
    //
  for (Int_t i=0;i<fId3N;++i) delete [] fId3Name[i];
  delete [] fId3Name;
  fId3Name = 0;
}

//_____________________________________________________________________________
AliITSv3::AliITSv3(const char *name, const char *title) : AliITS(name, title){
    //
    // Standard constructor for ITS
    // 
    fId3N = 6;
    fId3Name = new char*[fId3N];
    fId3Name[0] = "ITS1";
    fId3Name[1] = "ITS2";
    fId3Name[2] = "ITS3";
    fId3Name[3] = "ITS4";
    fId3Name[4] = "ITS5";
    fId3Name[5] = "ITS6";
    fMinorVersionV3=1;
}
 
//_____________________________________________________________________________
void AliITSv3::CreateGeometry(){
  //
  // Create ITS geometry for version 3
  //
  //

  
  const Float_t xx[14] = {  0.000,  0.000,-14.002, -6.288,-25.212,-16.292,
                          -35.713,-26.401,-45.340,-36.772,-18.740,-12.814,
                          -14.358,  0.000};
  const Float_t yy[14] = {  0.000, 27.056, 31.408, 25.019, 27.768, 22.664,
                           22.420, 18.727, 15.479, 13.680, -9.984, -6.175,
                           -3.775, 0.000 };
  const Float_t xbeg[13] = {  0.000, -0.352,-12.055, -8.755,-23.035,-19.085,
                            -33.362,-28.859,-42.774,-36.644,-18.352,-13.085,
                            -13.426 };
  const Float_t ybeg[13] = {  0.386, 27.165, 29.795, 25.377, 26.480, 22.632,
                             21.487, 18.305, 14.940, 13.509, -9.735, -5.755,
                             -3.53 };
  const Float_t xend[13] = {  0.000,-11.588, -8.208,-22.709,-18.738,-33.184,
                            -28.719,-42.756,-37.027,-19.002,-13.235,-13.837,
                              -.373 };
  const Float_t yend[13] = { 26.688, 30.658, 26.609, 27.405, 23.935, 22.452,
                             19.646, 15.922, 13.733, -9.639, -6.446, -4.585,
                              -.098 };
  const Float_t xarc[13] = { -0.500,-13.248,-13.505,-18.622,-37.171,-42.671,
                            -28.977,-33.178,-19.094,-22.781, -8.655,-11.736,
                             -0.500 };
  const Float_t yarc[13] = {  0.500, -4.093, -5.911, -9.200, 13.162, 15.543,
                             19.109, 22.066, 23.446, 27.024, 26.184, 30.294,
                             26.802 };
  const Float_t rarc[13] = { 0.5,0.7,0.5,0.5,0.7,0.5,0.7,
                             0.5,0.7,0.5,0.7,0.5,0.5 };
  const Float_t  rr     =   4.064516;
  const Float_t  tteta  =  63.00;
  const Float_t  pphi   = -35.00;
  const Float_t  gteta  =  87.78;
  const Double_t degrad = kPI/180.;
  const Double_t raddeg = 180./kPI;
  const Double_t twopi  = 2*kPI;
  
  Double_t biga, bigb;
  Float_t  dcei[3], dela[3], dchi[3], dpcb[3], darc[5], 
           dfra[10], dcer[3], dkap[3], dpla[3],
           xccc, yccc, aphi, dcop[3], dtra[3], dsil[3], 
           atheta1011, dbus[3], dtub[3], dwat[3],
           depx[3], dits[3], atheta1314, atheta1213, atheta1112, 
           dsup[3], xtra[8], ytra[8], ztra[8], dsrv[3];
  Double_t biga1, bigb1;
  Float_t  runo, xpos, ypos, zpos, rtwo, aphi1, aphi2, 
           dtra1[3], dtra2[3], dtra3[3],
           dtra4[3], dbox1[3], dbox2[3];
  Int_t    jbox1, jbox2;
  Float_t  xtra1[6], ytra1[6], ztra1[6];
  Int_t    i;
  Float_t  xpos1, ypos1;
  Int_t    j;
  Float_t  angle, dcone[5], dtube[3], dpgon[10];
  Float_t  rzero, xzero, yzero;
  Double_t coeffa, coeffb, coeffc;
  Int_t    idrotm[5250];
  Float_t  atheta, offset;
  Float_t  offset1, offset2, dgh[15];
  Float_t  xcc, ycc, sep, atheta12, atheta23, atheta34, atheta45, atheta56, 
           atheta67, atheta78, atheta89, xxm, dal1[3], dal2[3];
  //Float_t  yos;
  Float_t  r1, r2, r3;
  Double_t xcc1, ycc1, xcc2, ycc2;
  Float_t  atheta910;
  const char natra[][5] ={ "TR01","TR02","TR03","TR04",
                           "TR05","TR06","TR07","TR08"};
  const char natra1[][5] ={"TR11","TR12","TR13","TR14",
                           "TR15","TR16","TR17","TR18",
                           "TR19","TR20","TR21","TR22",
                           "TR23","TR24","TR25","TR26"};
  const char natra2[][5] ={"TR31","TR32","TR33","TR34","TR35","TR36"};
  const char natra3[][5] ={"TR41","TR42","TR43","TR44","TR45","TR46"};
  const char natra4[][5] ={"TR51","TR52","TR53","TR54","TR55","TR56",
                           "TR57","TR58","TR59","TR60","TR61","TR62",
                           "TR63","TR64","TR65","TR66"};
  
  Int_t *idtmed = fIdtmed->GetArray()-199;
  
  // --- Define a ghost volume containing the whole ITS and fill it with air
  //     or vacuum 
  
  dgh[0]  = 0.0;
  dgh[1]  = 360.0;
  dgh[2]  = 4.0;
  dgh[3]  = -70.0;
  dgh[4]  = 49.999;
  dgh[5]  = 49.999;
  dgh[6]  = -25.0;
  dgh[7]  = 3.0;
  dgh[8]  = 49.999;
  dgh[9]  = 25.0;
  dgh[10] = 3.0;
  dgh[11] = 49.999;
  dgh[12] = 70.0;
  dgh[13] = 49.999;
  dgh[14] = 49.999;
  gMC->Gsvolu("ITSV", "PCON", idtmed[275], dgh, 15);
  
  // --- Place the ghost volume in its mother volume (ALIC) and make it 
  //     invisible 
  
  gMC->Gspos("ITSV", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  gMC->Gsatt("ITSV", "SEEN", 0);
  
  //************************************************************************
  //*                                                                      *
  //*                               P I X E L S                            *
  //*                               ===========                            *
  //*                                                                      *
  //************************************************************************
  
  //        GOTO 2345             ! skip ITS layer no. 1 and 2 
  
  // --- Define a ghost volume containing the Silicon Pixel Detectors 
  //     (layer #1 and #2) and fill it with air or vacuum 
  
  xxm     = (49.999-3)/(70-25);
  dgh[0]  = 0.0;
  dgh[1]  = 360.0;
  dgh[2]  = 4.0;
  dgh[3]  = -25.-(9.-3.01)/xxm;
  dgh[4]  = 9.0;
  dgh[5]  = 9.0;
  dgh[6]  = -25.0;
  dgh[7]  = 3.01;
  dgh[8]  = 9.0;
  dgh[9]  = 25.0;
  dgh[10] = 3.01;
  dgh[11] = 9.0;
  dgh[12] = 25+(9-3.01)/xxm;
  dgh[13] = 9.0;
  dgh[14] = 9.0;
  gMC->Gsvolu("IT12", "PCON", idtmed[275], dgh, 15);
  
  // --- Place the ghost volume in its mother volume (ITSV) and make it 
  //     invisible 
  
  gMC->Gspos("IT12", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  gMC->Gsatt("IT12", "SEEN", 0);
  
  // --- Define a ghost volume containing a single element of layer #1 
  //     and fill it with air or vacuum 
  
  dbox1[0] = 0.005+0.01+0.0075;
  dbox1[1] = .79;
  dbox1[2] = 12.67;
  gMC->Gsvolu("IPV1", "BOX ", idtmed[203], dbox1, 3);
  
  //--Divide each element of layer #1 in three ladders along the beam direction

  gMC->Gsdvn("IPB1", "IPV1", 3, 3);
  
  // --- Make the ghost volumes invisible 
  
  gMC->Gsatt("IPV1", "SEEN", 0);
  gMC->Gsatt("IPB1", "SEEN", 0);
  
  // --- Define a volume containing the chip of pixels (silicon, layer #1)

  dchi[0] = 0.005;
  dchi[1] = 0.79;
  dchi[2] = dbox1[2] / 3.;
  gMC->Gsvolu("ICH1", "BOX ", idtmed[200], dchi, 3);
  
  // --- Define a volume containing the bus of pixels (silicon, layer #1) 
  
  dbus[0] = 0.01;
  dbus[1] = 0.64;
  dbus[2] = 4.19;
  gMC->Gsvolu("IBU1", "BOX ", idtmed[201], dbus, 3);
  
  // --- Define a volume containing the sensitive part of pixels 
  //     (silicon, layer #1) 
  
  dits[0] = 0.0075;
  dits[1] = 0.64;
  dits[2] = 4.19;
  gMC->Gsvolu("ITS1", "BOX ", idtmed[199], dits, 3);

  // --- Place the chip into its mother (IPB1) 
  
  xpos = dbox1[0] - dchi[0];
  ypos = 0.0;
  zpos = 0.0;
  gMC->Gspos("ICH1", 1, "IPB1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the sensitive volume into its mother (IPB1) 
  
  xpos = dbox1[0] - dchi[0] * 2. - dits[0];
  ypos = dchi[1] - dits[1];
  zpos = -(dchi[2] - dits[2]);
  gMC->Gspos("ITS1", 1, "IPB1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the bus into its mother (IPB1) 
  
  xpos = dbox1[0] - dchi[0] * 2. - dits[0] * 2. - dbus[0];
  ypos = dchi[1] - dbus[1];
  zpos = -(dchi[2] - dbus[2]);
  gMC->Gspos("IBU1", 1, "IPB1", xpos, ypos, zpos, 0, "ONLY");

  // --- Define a ghost volume containing a single element of layer #2 
  //     and fill it with air or vacuum 
  
  dbox2[0] = 0.005+0.01+0.0075;
  dbox2[1] = 0.79;
  dbox2[2] = 16.91;
  gMC->Gsvolu("IPV2", "BOX ", idtmed[203], dbox2, 3);

  //--Divide each element of layer #2 in four ladders along the beam direction
  
  gMC->Gsdvn("IPB2", "IPV2", 4, 3);
  
  // --- Make the ghost volumes invisible 
  
  gMC->Gsatt("IPV2", "SEEN", 0);
  gMC->Gsatt("IPB2", "SEEN", 0);
  
  // --- Define a volume containing the chip of pixels (silicon, layer #2) 
  
  dchi[0] = 0.005;
  dchi[1] = 0.79;
  dchi[2] = dbox2[2] / 4.;
  gMC->Gsvolu("ICH2", "BOX ", idtmed[200], dchi, 3);

  // --- Define a volume containing the bus of pixels (silicon, layer #2) 
  
  dbus[0] = 0.01;
  dbus[1] = 0.64;
  dbus[2] = 4.19;
  gMC->Gsvolu("IBU2", "BOX ", idtmed[201], dbus, 3);

  // --- Define a volume containing the sensitive part of pixels 
  //     (silicon, layer #2) 
  
  dits[0] = 0.0075;
  dits[1] = 0.64;
  dits[2] = 4.19;
  gMC->Gsvolu("ITS2", "BOX ", idtmed[199], dits, 3);

  // --- Place the chip into its mother (IPB2) 
  
  xpos = dbox1[0] - dbus[0] * 2. - dits[0] * 2. - dchi[0];
  ypos = 0.0;
  zpos = 0.0;
  gMC->Gspos("ICH2", 1, "IPB2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the sensitive volume into its mother (IPB2) 
  
  xpos = dbox1[0] - dbus[0] * 2. - dits[0];
  ypos = -(dchi[1] - dits[1]);
  zpos = -(dchi[2] - dits[2]);
  gMC->Gspos("ITS2", 1, "IPB2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the bus into its mother (IPB2) 
  
  xpos = dbox1[0] - dbus[0];
  ypos = -(dchi[1] - dbus[1]);
  zpos = -(dchi[2] - dbus[2]);
  gMC->Gspos("IBU2", 1, "IPB2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Define a generic segment of an element of the mechanical support 
  
  dsup[0] = 0.0;
  dsup[1] = 0.0;
  dsup[2] = 0.0;
  gMC->Gsvolu("SPIX", "BOX ", idtmed[202], dsup, 0);
  
  // --- Define a generic arc of an element of the mechanical support 
  
  darc[0] = 0.0;
  darc[1] = 0.0;
  darc[2] = 0.0;
  gMC->Gsvolu("SARC", "TUBS", idtmed[202], darc, 0);
  
  // --- Define the mechanical supports of layers #1 and #2 and place the 
  //     elements of the layers in it 
  
  jbox1 = 0;
  // counter over the number of elements of layer #1 ( 
  jbox2 = 0;
  
  // counter over the number of elements of layer #2 ( 
  for (i = 1; i <= 10; ++i) {
    
    // --- Place part # 1-2 (see sketch) 
    
    // number of carbon fiber supports (see sketch) 
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[0] - xbeg[0]) * (xend[0] - xbeg[0]) + 
                          (yend[0] - ybeg[0]) * (yend[0] - ybeg[0])   ) / 20.;
    dsup[2] = 25.0;
    xcc     = (  xx[0] +   xx[1]) / 20.;
    ycc     = (  yy[0] +   yy[1]) / 20.;
    xccc    = (xbeg[0] + xend[0]) / 20.;
    yccc    = (ybeg[0] + yend[0]) / 20.;
    if (xx[0] == xx[1]) {
      offset2 = 0.;
    } else {
      r1 = yy[1] - yy[0];
      r2 = xx[1] - xx[0];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    } // end if xx[0] == xx[1]
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) +
             ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.0;
    atheta12 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1100], 90., atheta12, 90.,
                                             atheta12 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 1, "IT12", xpos, ypos, zpos,
                         idrotm[(i-1) * 13 + 1100], "ONLY", dsup, 3);

    // --- Place part # 2-3 (see sketch) 
    
    offset1 = -35.0;
    dsup[0] = 0.01;
    dsup[1] = TMath::Sqrt((xend[1] - xbeg[1]) * (xend[1] - xbeg[1]) + 
                          (yend[1] - ybeg[1]) * (yend[1] - ybeg[1])) / 20.;
    dsup[2] = 25.0;
    xcc     = (  xx[1] +   xx[2]) / 20.;
    ycc     = (  yy[1] +   yy[2]) / 20.;
    xccc    = (xbeg[1] + xend[1]) / 20.;
    yccc    = (ybeg[1] + yend[1]) / 20.;
    if (xx[1] == xx[2]) {
      offset2 = 0.;
    } else {
      r1 = yy[2] - yy[1];
      r2 = xx[2] - xx[1];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    } // end if xx[1] == xx[2]
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) + 
             ypos1 * TMath::Sin(gteta * degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.0;
    atheta23 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1101], 90., atheta23, 90.,
                                             atheta23 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 2, "IT12", xpos, ypos, zpos,
		idrotm[(i-1) * 13 + 1101], "ONLY", dsup, 3);

    // --- Place an element of layer #2

    biga   = (yy[2] - yy[1]) / (xx[2] - xx[1]);
    bigb   = (xx[2] * yy[1] - xx[1] * yy[2]) / (xx[2] - xx[1]) / 10.;
    coeffa = biga * biga + 1.;
    coeffb = biga * bigb - biga * ycc - xcc;
    coeffc = xcc * xcc + ycc * ycc - ycc * 2. * bigb + 
            bigb * bigb - 0.08964*0.08964;
    xcc1   = (-coeffb + TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) /
                                                                   coeffa;
    ycc1   = biga * xcc1 + bigb;
    biga1  = -1. / biga;
    bigb1  = xcc1 / biga + ycc1;
    coeffa = biga1 * biga1 + 1.;
    coeffb = biga1 * bigb1 - biga1 * ycc1 - xcc1;
    coeffc = xcc1 * xcc1 + ycc1 * ycc1 - ycc1 * 2. * bigb1 + 
            bigb1 * bigb1 - (dsup[0] + dbox2[0]) * (dsup[0] + dbox2[0]);
    xcc2   = (-coeffb + TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / 
                                                                    coeffa;
    ycc2   =  biga1 * xcc2 + bigb1;
    xpos1  =   xcc2 * TMath::Cos(aphi) - ycc2 * TMath::Sin(aphi) + xzero;
    ypos1  =   xcc2 * TMath::Sin(aphi) + ycc2 * TMath::Cos(aphi) + yzero;
    xpos   =  xpos1 * TMath::Cos(gteta * degrad) + 
              ypos1 * TMath::Sin(gteta *degrad);
    ypos   = -xpos1 * TMath::Sin(gteta * degrad) + 
              ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.0;
    ++jbox2;
    gMC->Gspos("IPV2", jbox2, "IT12", xpos, ypos, zpos, 
                         idrotm[(i-1) * 13 + 1101], "ONLY");

    // --- Place part # 3-4 (see sketch)

    offset1 = -35.0;
    dsup[0] = 0.01;
    dsup[1] = TMath::Sqrt((xend[2] - xbeg[2]) * (xend[2] - xbeg[2]) + 
                          (yend[2] - ybeg[2]) * (yend[2] - ybeg[2])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[1] + xx[2]) / 20.;
    ycc     = (yy[1] + yy[2]) / 20.;
    xccc    = (xbeg[2] + xend[2]) / 20.;
    yccc    = (ybeg[2] + yend[2]) / 20.;
    if (xx[2] == xx[3]) {
      offset2 = 0.;
    } else {
      r1 = yy[3] - yy[2];
      r2 = xx[3] - xx[2];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    } // end if xx[2] == xx[3]
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) + 
             ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.0;
    atheta34 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1102], 90., atheta34, 90., 
                                             atheta34 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 3, "IT12", xpos, ypos, zpos, 
                         idrotm[(i-1) * 13 + 1102], "ONLY", dsup, 3);

    // --- Place part # 4-5 (see sketch)

    offset1 = -35.0;
    dsup[0] = 0.01;
    dsup[1] = TMath::Sqrt((xend[3] - xbeg[3]) * (xend[3] - xbeg[3]) + 
                          (yend[3] - ybeg[3]) * (yend[3] - ybeg[3])) / 20.;
    dsup[2] = 25.0;
    xcc     = (  xx[3] +   xx[4]) / 20.;
    ycc     = (  yy[3] +   yy[4]) / 20.;
    xccc    = (xbeg[3] + xend[3]) / 20.;
    yccc    = (ybeg[3] + yend[3]) / 20.;
    if (xx[3] == xx[4]) {
      offset2 = 0.;
    } else {
      r1 = yy[4] - yy[3];
      r2 = xx[4] - xx[3];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    } // end if xx[3] == xx[4]
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) + 
             ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.0;
    atheta45 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1103], 90., atheta45, 90., 
                                        atheta45 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 4, "IT12", xpos, ypos, zpos, 
                        idrotm[(i-1) * 13 + 1103], "ONLY", dsup, 3);

    // --- Place an element of layer #2

    biga   = (yy[4] - yy[3]) / (xx[4] - xx[3]);
    bigb   = (xx[4] * yy[3] - xx[3] * yy[4]) / (xx[4] - xx[3]) / 10.;
    coeffa = biga * biga + 1.;
    coeffb = biga * bigb - biga * ycc - xcc;
    coeffc = xcc * xcc + ycc * ycc - ycc * 2. * bigb + 
            bigb * bigb - .014285030400000001;
    xcc1   = (-coeffb - TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / 
                                                                      coeffa;
    ycc1   = biga * xcc1 + bigb;
    biga1  = -1. / biga;
    bigb1  = xcc1 / biga + ycc1;
    coeffa = biga1 * biga1 + 1.;
    coeffb = biga1 * bigb1 - biga1 * ycc1 - xcc1;
    coeffc = xcc1 * xcc1 + ycc1 * ycc1 - ycc1 * 2. * bigb1 + 
            bigb1 * bigb1 - (dsup[0] + dbox2[0]) * (dsup[0] + dbox2[0]);
    xcc2   = (-coeffb + TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / 
                                                                     coeffa;
    ycc2   =  biga1 * xcc2 + bigb1;
    xpos1  =   xcc2 * TMath::Cos(aphi) - ycc2 * TMath::Sin(aphi) + xzero;
    ypos1  =   xcc2 * TMath::Sin(aphi) + ycc2 * TMath::Cos(aphi) + yzero;
    xpos   =  xpos1 * TMath::Cos(gteta * degrad) + 
              ypos1 * TMath::Sin(gteta *degrad);
    ypos   = -xpos1 * TMath::Sin(gteta * degrad) + 
              ypos1 * TMath::Cos(gteta * degrad);
    zpos   = 0.0;
    ++jbox2;
    gMC->Gspos("IPV2", jbox2, "IT12", xpos, ypos, zpos, 
                           idrotm[(i-1) * 13 + 1103], "ONLY");

    // --- Place part # 5-6 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[4] - xbeg[4]) * (xend[4] - xbeg[4]) + 
                          (yend[4] - ybeg[4]) * (yend[4] - ybeg[4])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[4] + xx[5]) / 20.;
    ycc     = (yy[4] + yy[5]) / 20.;
    xccc    = (xbeg[4] + xend[4]) / 20.;
    yccc    = (ybeg[4] + yend[4]) / 20.;
    if (xx[4] == xx[5]) {
      offset2 = 0.;
    } else {
      r1 = yy[5] - yy[4];
      r2 = xx[5] - xx[4];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    }
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) + 
             ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta56 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1104], 90., atheta56, 90., 
                                              atheta56 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 5, "IT12", xpos, ypos, zpos,
                          idrotm[(i-1) * 13 + 1104], "ONLY", dsup, 3);

    // --- Place part # 6-7 (see sketch) 
    
    offset1 = -35.0;
    dsup[0] = 0.01;
    dsup[1] = TMath::Sqrt((xend[5] - xbeg[5]) * (xend[5] - xbeg[5]) + 
                          (yend[5] - ybeg[5]) * (yend[5] - ybeg[5])) / 20.;
    dsup[2] = 25.0;
    xcc     = (xx[5] + xx[6]) / 20.;
    ycc     = (yy[5] + yy[6]) / 20.;
    xccc    = (xbeg[5] + xend[5]) / 20.;
    yccc    = (ybeg[5] + yend[5]) / 20.;
    if (xx[5] == xx[6]) {
      offset2 = 0.;
    } else {
      r1 = yy[6] - yy[5];
      r2 = xx[6] - xx[5];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    } // end if xx[5] == xx[6]
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) + 
             ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta67 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1105], 90., atheta67, 90., 
                                             atheta67 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 6, "IT12", xpos, ypos, zpos, 
                         idrotm[(i-1) * 13 + 1105], "ONLY", dsup, 3);

    // --- Place an element of layer #2 
    
    biga   = (yy[6] - yy[5]) / (xx[6] - xx[5]);
    bigb   = (xx[6] * yy[5] - xx[5] * yy[6]) / (xx[6] - xx[5]) / 10.;
    coeffa = biga * biga + 1.;
    coeffb = biga * bigb - biga * ycc - xcc;
    coeffc = xcc * xcc + ycc * ycc - ycc * 2. * bigb + 
            bigb * bigb - .014285030400000001;
    xcc1   = (-coeffb - TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / 
                                                                      coeffa;
    ycc1   = biga * xcc1 + bigb;
    biga1  = -1. / biga;
    bigb1  = xcc1 / biga + ycc1;
    coeffa = biga1 * biga1 + 1.;
    coeffb = biga1 * bigb1 - biga1 * ycc1 - xcc1;
    coeffc = xcc1 * xcc1 + ycc1 * ycc1 - ycc1 * 2. * bigb1 +
            bigb1 * bigb1 - (dsup[0] + dbox2[0]) * (dsup[0] + dbox2[0]);
    xcc2   = (-coeffb - TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / 
                                                                   coeffa;
    ycc2   =  biga1 * xcc2 + bigb1;
    xpos1  =   xcc2 * TMath::Cos(aphi) - ycc2 * TMath::Sin(aphi) + xzero;
    ypos1  =   xcc2 * TMath::Sin(aphi) + ycc2 * TMath::Cos(aphi) + yzero;
    xpos   =  xpos1 * TMath::Cos(gteta * degrad) + 
              ypos1 * TMath::Sin(gteta *degrad);
    ypos   = -xpos1 * TMath::Sin(gteta * degrad) + 
              ypos1 * TMath::Cos(gteta * degrad);
    zpos   = 0.0;
    ++jbox2;
    gMC->Gspos("IPV2", jbox2, "IT12", xpos, ypos, zpos,
                       idrotm[(i-1) * 13 + 1105], "ONLY");
    
    // --- Place part # 7-8 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[6] - xbeg[6]) * (xend[6] - xbeg[6]) +
                          (yend[6] - ybeg[6]) * (yend[6] - ybeg[6])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[6] + xx[7]) / 20.;
    ycc     = (yy[6] + yy[7]) / 20.;
    xccc    = (xbeg[6] + xend[6]) / 20.;
    yccc    = (ybeg[6] + yend[6]) / 20.;
    if (xx[6] == xx[7]) {
      offset2 = 0.;
    } else {
      r1 = yy[7] - yy[6];
      r2 = xx[7] - xx[6];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    } // end if
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) + 
             ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta78 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1106], 90., atheta78, 90., 
                                             atheta78 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 7, "IT12", xpos, ypos, zpos, 
                           idrotm[(i-1) * 13 + 1106], "ONLY", dsup, 3);
    
    // --- Place part # 8-9 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[7] - xbeg[7]) * (xend[7] - xbeg[7]) + 
                          (yend[7] - ybeg[7]) * (yend[7] - ybeg[7])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[7] + xx[8]) / 20.;
    ycc     = (yy[7] + yy[8]) / 20.;
    xccc    = (xbeg[7] + xend[7]) / 20.;
    yccc    = (ybeg[7] + yend[7]) / 20.;
    if (xx[1] == xx[2]) {
      offset2 = 0.;
    } else {
      r1 = yy[8] - yy[7];
      r2 = xx[8] - xx[7];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    }
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero =     rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero =     rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 =   xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 =   xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  =  xpos1 * TMath::Cos(gteta * degrad) + 
             ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + 
             ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta89 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1107], 90., atheta89, 90., 
                                             atheta89 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 8, "IT12", xpos, ypos, zpos, 
                         idrotm[(i-1) * 13 + 1107], "ONLY", dsup, 3);
    
    // --- Place an element of layer #2 
    
    biga   = (yy[8] - yy[7]) / (xx[8] - xx[7]);
    bigb   = (xx[8] * yy[7] - xx[7] * yy[8]) / (xx[8] - xx[7]) / 10.;
    coeffa = biga * biga + 1.;
    coeffb = biga * bigb - biga * ycc - xcc;
    coeffc = xcc * xcc + ycc * ycc - ycc * 2. * bigb +
            bigb * bigb - .014285030400000001;
    xcc1   = (-coeffb - TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / 
                                                                   coeffa;
    ycc1   = biga * xcc1 + bigb;
    biga1  = -1. / biga;
    bigb1  = xcc1 / biga + ycc1;
    coeffa = biga1 * biga1 + 1.;
    coeffb = biga1 * bigb1 - biga1 * ycc1 - xcc1;
    coeffc = xcc1 * xcc1 + ycc1 * ycc1 - ycc1 * 2. * bigb1 + 
            bigb1 * bigb1 - (dsup[0] + dbox2[0]) * (dsup[0] + dbox2[0]);
    xcc2   = (-coeffb - TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / 
                                                                    coeffa;
    ycc2   =  biga1 * xcc2 + bigb1;
    xpos1  =   xcc2 * TMath::Cos(aphi) - ycc2 * TMath::Sin(aphi) + xzero;
    ypos1  =   xcc2 * TMath::Sin(aphi) + ycc2 * TMath::Cos(aphi) + yzero;
    xpos   =  xpos1 * TMath::Cos(gteta * degrad) + 
              ypos1 * TMath::Sin(gteta *degrad);
    ypos   = -xpos1 * TMath::Sin(gteta * degrad) + 
              ypos1 * TMath::Cos(gteta * degrad);
    zpos   = 0.0;
    ++jbox2;
    gMC->Gspos("IPV2", jbox2, "IT12", xpos, ypos, zpos, 
                      idrotm[(i-1) * 13 + 1107], "ONLY");

    // --- Place part # 9-10 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[8] - xbeg[8]) * (xend[8] - xbeg[8]) +
                          (yend[8] - ybeg[8]) * (yend[8] - ybeg[8])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[8] + xx[9]) / 20.;
    ycc     = (yy[8] + yy[9]) / 20.;
    xccc    = (xbeg[8] + xend[8]) / 20.;
    yccc    = (ybeg[8] + yend[8]) / 20.;
    if (xx[8] == xx[9]) {
      offset2 = 0.;
    } else {
      r1 = yy[9] - yy[8];
      r2 = xx[9] - xx[8];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    }
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 = xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 = xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta910 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1108], 90., atheta910, 90., atheta910 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 9, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1108], "ONLY", dsup, 3);
    
    // --- Place part # 10-11 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[9] - xbeg[9]) * (xend[9] - xbeg[9]) + (yend[9] - ybeg[9]) * (yend[9] - ybeg[9])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[9] + xx[10]) / 20.;
    ycc     = (yy[9] + yy[10]) / 20.;
    xccc    = (xbeg[9] + xend[9]) / 20.;
    yccc    = (ybeg[9] + yend[9]) / 20.;
    if (xx[9] == xx[10]) {
      offset2 = 0.;
    } else {
      r1 = yy[10] - yy[9];
      r2 = xx[10] - xx[9];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    }
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 = xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 = xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta1011 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1109], 90., atheta1011, 90.,atheta1011 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 10, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1109], "ONLY", dsup, 3);
    
    // --- Place part # 13-14 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[12] - xbeg[12]) * (xend[12] - xbeg[12]) + (yend[12] - ybeg[12]) * (yend[12] - ybeg[12])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[12] + xx[13]) / 20.;
    ycc     = (yy[12] + yy[13]) / 20.;
    xccc    = (xbeg[12] + xend[12]) / 20.;
    yccc    = (ybeg[12] + yend[12]) / 20.;
    if (xx[12] == xx[13]) {
      offset2 = 0.;
    } else {
      r1 = yy[12] - yy[13];
      r2 = xx[12] - xx[13];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    }
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 = xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 = xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta1314 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1112], 90., atheta1314, 90.,atheta1314 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 13, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1112], "ONLY", dsup, 3);
    
    // --- Place an element of layer #1 
    
    biga   = (yy[13] - yy[12]) / (xx[13] - xx[12]);
    bigb   = (xx[13] * yy[12] - xx[12] * yy[13]) / (xx[13] - xx[12]) / 10.;
    coeffa = biga * biga + 1.;
    coeffb = biga * bigb - biga * ycc - xcc;
    coeffc = xcc * xcc + ycc * ycc - ycc * 2. * bigb + bigb * bigb - .050216328100000006;
    xcc1   = (-coeffb + TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / coeffa;
    ycc1   = biga * xcc1 + bigb;
    biga1  = -1. / biga;
    bigb1  = xcc1 / biga + ycc1;
    coeffa = biga1 * biga1 + 1.;
    coeffb = biga1 * bigb1 - biga1 * ycc1 - xcc1;
    coeffc = xcc1 * xcc1 + ycc1 * ycc1 - ycc1 * 2. * bigb1 + bigb1 * bigb1 - (dsup[0] + dbox1[0]) * (dsup[0] + dbox1[0]);
    xcc2   = (-coeffb + TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / coeffa;
    ycc2   = biga1 * xcc2 + bigb1;
    xpos1  = xcc2 * TMath::Cos(aphi) - ycc2 * TMath::Sin(aphi) + xzero;
    ypos1  = xcc2 * TMath::Sin(aphi) + ycc2 * TMath::Cos(aphi) + yzero;
    xpos   = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos   = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos   = 0.;
    ++jbox1;
    gMC->Gspos("IPV1", jbox1, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1112], "ONLY");
    
    // --- Place part # 12-13 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[11] - xbeg[11]) * (xend[11] - xbeg[11]) + (yend[11] - ybeg[11]) * (yend[11] - ybeg[11])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[11] + xx[12]) / 20.;
    ycc     = (yy[11] + yy[12]) / 20.;
    xccc    = (xbeg[11] + xend[11]) / 20.;
    yccc    = (ybeg[11] + yend[11]) / 20.;
    if (xx[11] == xx[12]) {
      offset2 = 0.;
    } else {
      r1 = yy[12] - yy[11];
      r2 = xx[12] - xx[11];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    }
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 = xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 = xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta1213 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1111], 90., atheta1213, 90.,atheta1213 + 90., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 12, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1111], "ONLY", dsup, 3);
    
    // --- Place part # 11-12 (see sketch) 
    
    offset1 = -35.;
    dsup[0] = .01;
    dsup[1] = TMath::Sqrt((xend[10] - xbeg[10]) * (xend[10] - xbeg[10]) + (yend[10] - ybeg[10]) * (yend[10] - ybeg[10])) / 20.;
    dsup[2] = 25.;
    xcc     = (xx[10] + xx[11]) / 20.;
    ycc     = (yy[10] + yy[11]) / 20.;
    xccc    = (xbeg[10] + xend[10]) / 20.;
    yccc    = (ybeg[10] + yend[10]) / 20.;
    if (xx[10] == xx[11]) {
      offset2 = 0.;
    } else {
      r1 = yy[11] - yy[10];
      r2 = xx[11] - xx[10];
      offset2 = TMath::ATan2(r1, r2) * raddeg - 90.;
    }
    aphi  = (pphi + (i-1) * 36.) * degrad;
    xzero = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1 = xccc * TMath::Cos(aphi) - yccc * TMath::Sin(aphi) + xzero;
    ypos1 = xccc * TMath::Sin(aphi) + yccc * TMath::Cos(aphi) + yzero;
    xpos  = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos  = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos  = 0.;
    atheta1112 = (i-1) * 36. + offset1 + offset2 - gteta;
    AliMatrix(idrotm[(i-1) * 13 + 1110], 270., atheta1112, 90., atheta1112 + 270., 0., 0.);
    gMC->Gsposp("SPIX", (i-1) * 13 + 11, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1110], "ONLY", dsup, 3);

    // --- Place an element of layer #1 
    
    biga   = (yy[11] - yy[10]) / (xx[11] - xx[10]);
    bigb   = (xx[11] * yy[10] - xx[10] * yy[11]) / (xx[11] - xx[10]) / 10.;
    coeffa = biga * biga + 1.;
    coeffb = biga * bigb - biga * ycc - xcc;
    coeffc = xcc * xcc + ycc * ycc - ycc * 2. * bigb + bigb * bigb - .0035712576000000002;
    xcc1   = (-coeffb + TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / coeffa;
    ycc1   = biga * xcc1 + bigb;
    biga1  = -1. / biga;
    bigb1  = xcc1 / biga + ycc1;
    coeffa = biga1 * biga1 + 1.;
    coeffb = biga1 * bigb1 - biga1 * ycc1 - xcc1;
    coeffc = xcc1 * xcc1 + ycc1 * ycc1 - ycc1 * 2. * bigb1 + bigb1 * bigb1 - (dsup[0] + dbox1[0]) * (dsup[0] + dbox1[0]);
    xcc2   = (-coeffb + TMath::Sqrt(coeffb * coeffb - coeffa * coeffc)) / coeffa;
    ycc2   = biga1 * xcc2 + bigb1;
    xpos1  = xcc2 * TMath::Cos(aphi) - ycc2 * TMath::Sin(aphi) + xzero;
    ypos1  = xcc2 * TMath::Sin(aphi) + ycc2 * TMath::Cos(aphi) + yzero;
    xpos   = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos   = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos   = 0.;
    ++jbox1;
    gMC->Gspos("IPV1", jbox1, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1110], "ONLY");
    
    // --- Place arc # 13 (between part 1-2 and part 2-3) (see sketch) 
    
    darc[0] = rarc[12] / 10. - .02;
    darc[1] = rarc[12] / 10.;
    darc[2] = 25.;
    darc[3] = atheta12 - (i-1) * 36.;
    darc[4] = atheta23 - (i-1) * 36.;
    xcc     = xarc[12] / 10.;
    ycc     = yarc[12] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 13, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1112], "ONLY", darc, 5);
    
    // --- Place arc # 12 (between part 2-3 and part 3-4) (see sketch) 
    
    darc[0] = rarc[11] / 10. - .02;
    darc[1] = rarc[11] / 10.;
    darc[2] = 25.;
    darc[3] = atheta23 + 90. - (i-1) * 36.;
    darc[4] = atheta34 + 90. - (i-1) * 36.;
    xcc     = xarc[11] / 10.;
    ycc     = yarc[11] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 12, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1111], "ONLY", darc, 5);
    
    // --- Place arc # 11 (between part 3-4 and part 4-5) (see sketch) 
    
    darc[0] = rarc[10] / 10. - .02;
    darc[1] = rarc[10] / 10.;
    darc[2] = 25.;
    darc[3] = atheta45 + 180. - (i-1) * 36.;
    darc[4] = atheta34 + 180. - (i-1) * 36.;
    xcc     = xarc[10] / 10.;
    ycc     = yarc[10] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 11, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1110], "ONLY", darc, 5);
    
    // --- Place arc # 10 (between part 4-5 and part 5-6) (see sketch) 
    
    darc[0] = rarc[9] / 10. - .02;
    darc[1] = rarc[9] / 10.;
    darc[2] = 25.;
    darc[3] = atheta45 - 90. - (i-1) * 36.;
    darc[4] = atheta56 - 90. - (i-1) * 36.;
    xcc     = xarc[9] / 10.;
    ycc     = yarc[9] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 10, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1109], "ONLY", darc, 5);
    
    // --- Place arc # 9 (between part 5-6 and part) (see sketch) 
    
    darc[0] = rarc[8] / 10. - .02;
    darc[1] = rarc[8] / 10.;
    darc[2] = 25.;
    darc[3] = atheta67 + 45. - (i-1) * 36.;
    darc[4] = atheta56 + 45. - (i-1) * 36.;
    xcc     = xarc[8] / 10.;
    ycc     = yarc[8] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 9, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1108], "ONLY", darc, 5);
    
    // --- Place arc # 8 (between part 6-7 and part 7-8) (see sketch) 
    
    darc[0] = rarc[7] / 10. - .02;
    darc[1] = rarc[7] / 10.;
    darc[2] = 25.;
    darc[3] = atheta67 - (i-1) * 36.;
    darc[4] = atheta78 - (i-1) * 36.;
    xcc     = xarc[7] / 10.;
    ycc     = yarc[7] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 8, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1107], "ONLY", darc, 5);
    
    // --- Place arc # 7 (between part 7-8 and part 8-9) (see sketch) 
    
    darc[0] = rarc[6] / 10. - .02;
    darc[1] = rarc[6] / 10.;
    darc[2] = 25.;
    darc[3] = atheta89 + 45. - (i-1) * 36.;
    darc[4] = atheta78 + 45. - (i-1) * 36.;
    xcc     = xarc[6] / 10.;
    ycc     = yarc[6] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 7, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1106], "ONLY", darc, 5);
    
    // --- Place arc # 6 (between part 8-9 and part 9-10) (see sketch) 
    
    darc[0] = rarc[5] / 10. - .02;
    darc[1] = rarc[5] / 10.;
    darc[2] = 25.;
    darc[3] = atheta89 + 45. - (i-1) * 36.;
    darc[4] = atheta910 + 45. - (i-1) * 36.;
    xcc     = xarc[5] / 10.;
    ycc     = yarc[5] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 6, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1105], "ONLY", darc, 5);
    
    // --- Place arc # 5 (between part 9-10 and part 10-11) 
    //     (see sketch) 
    
    darc[0] = rarc[4] / 10. - .02;
    darc[1] = rarc[4] / 10.;
    darc[2] = 25.;
    darc[3] = atheta1011 + 45. - (i-1) * 36.;
    darc[4] = atheta910 + 45. - (i-1) * 36.;
    xcc     = xarc[4] / 10.;
    ycc     = yarc[4] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 5, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1104], "ONLY", darc, 5);
    
    // --- Place arc # 4 (between part 10-11 and part 11-12) 
    //     (see sketch) 
    
    darc[0] = rarc[3] / 10. - .02;
    darc[1] = rarc[3] / 10.;
    darc[2] = 25.;
    darc[3] = atheta1112 - 45. - (i-1) * 36.;
    darc[4] = atheta1011 - 225. - (i-1) * 36.;
    xcc     = xarc[3] / 10.;
    ycc     = yarc[3] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 4, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1103], "ONLY", darc, 5);
    
    // --- Place arc # 3 (between part 11-12 and part 12-13) 
    //     (see sketch) 
    
    darc[0] = rarc[2] / 10. - .02;
    darc[1] = rarc[2] / 10.;
    darc[2] = 25.;
    darc[3] = atheta1112 - 90. - (i-1) * 36.;
    darc[4] = atheta1213 - 90. - (i-1) * 36.;
    xcc     = xarc[2] / 10.;
    ycc     = yarc[2] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 3, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1102], "ONLY", darc, 5);
    
    // --- Place arc # 2 (between part 12-13 and part 13-14) 
    //     (see sketch) 
    
    darc[0] = rarc[1] / 10. - .02;
    darc[1] = rarc[1] / 10.;
    darc[2] = 25.;
    darc[3] = atheta1213 + 135. - (i-1) * 36.;
    darc[4] = atheta1314 + 165. - (i-1) * 36.;
    xcc     = xarc[1] / 10.;
    ycc     = yarc[1] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 2, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1101], "ONLY", darc, 5);
    
    // --- Place arc # 1 (between part 13-14 and part 1-2) 
    //     (see sketch) 
    
    darc[0] = rarc[0] / 10. - .02;
    darc[1] = rarc[0] / 10.;
    darc[2] = 25.;
    darc[3] = atheta12 + 45. - (i-1) * 36.;
    darc[4] = atheta1314 - (i-1) * 36.;
    xcc     = xarc[0] / 10.;
    ycc     = yarc[0] / 10.;
    aphi    = (pphi + (i-1) * 36.) * degrad;
    xzero   = rr * TMath::Cos((tteta + (i-1) * 36.) * degrad);
    yzero   = rr * TMath::Sin((tteta + (i-1) * 36.) * degrad);
    xpos1   = xcc * TMath::Cos(aphi) - ycc * TMath::Sin(aphi) + xzero;
    ypos1   = xcc * TMath::Sin(aphi) + ycc * TMath::Cos(aphi) + yzero;
    xpos    = xpos1 * TMath::Cos(gteta * degrad) + ypos1 * TMath::Sin(gteta *degrad);
    ypos    = -xpos1 * TMath::Sin(gteta * degrad) + ypos1 * TMath::Cos(gteta * degrad);
    zpos    = 0.;
    gMC->Gsposp("SARC", (i-1) * 13 + 1, "IT12", xpos, ypos, zpos, idrotm[(i-1) * 13 + 1100], "ONLY", darc, 5);
    
  }
  //************************************************************************
  //*                                                                      *
  //*                               D R I F T S                            *
  //*                               ===========                            *
  //*                                                                      *
  //************************************************************************
  
  // --- Define a ghost volume containing the Silicon Drift Detectors 
  //     (layer #3 and #4) and fill it with air or vacuum 
  
  xxm    = (49.999-3.)/(70.-25.);
  dgh[0] = 0;
  dgh[1] = 360;
  dgh[2] = 4;
  dgh[3] = -25.-(9.-3.01)/xxm-(9.01-9.)/xxm-(27.-9.01)/xxm;
  dgh[4] = 27.;
  dgh[5] = 27.;
  dgh[6] = -25.-(9.-3.01)/xxm-(9.01-9.)/xxm;
  dgh[7] = 9.01;
  dgh[8] = 27.;
  dgh[9] = 25.+(9.-3.01)/xxm+(9.01-9.)/xxm;
  dgh[10] = 9.01;
  dgh[11] = 27.;
  dgh[12] = 25.+(9.-3.01)/xxm+(9.01-9.)/xxm+(27.-9.01)/xxm;
  dgh[13] = 27.;
  dgh[14] = 27.;
  gMC->Gsvolu("IT34", "PCON", idtmed[275], dgh, 15);
  
  // --- Place the ghost volume in its mother volume (ITSV) and make it 
  //     invisible 
  
  gMC->Gspos("IT34", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  gMC->Gsatt("IT34", "SEEN", 0);
  
  // --- Layer #3 
  
  //        GOTO 3456           ! skip ITS layer no. 3 
  
  //--- Define a ghost volume containing a single ladder of layer #3 (with the
  //     smaller lenght of ribs) and fill it with air or vacuum 
  
  dbox1[0] = 0.5+(0.0172+0.03+0.0252+0.04+0.003);
  dbox1[1] = 3.85;
  // the widest element is the sensitive element 
  dbox1[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IDV1", "BOX ", idtmed[228], dbox1, 3);
  
  // --- Make the ghost volume invisible 
  
  gMC->Gsatt("IDV1", "SEEN", 0);
  
  // --- Define a volume containing the sensitive part of drifts 
  //     (silicon, layer #3) 
  
  dits[0] = .0172;
  // see material budget report by G. Feofilov 
  dits[1] = 3.85;
  dits[2] = 4.35;
  gMC->Gsvolu("ITS3", "BOX ", idtmed[224], dits, 3);
  
  //--- Define the part of the (smaller) rib between two sensitive parts made of
  //     carbon (layer #3) 
  
  dsup[0] = .5 - dits[0];
  dsup[1] = .01;
  dsup[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR11", "BOX ", idtmed[227], dsup, 3);
  
  //--- Define the first part of the (smaller) rib between two sensitive parts
  //     made of aluminum (layer #3) 
  
  dal1[0] = .5 - dits[0];
  dal1[1] = 0.00096/2.;
  dal1[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR12", "BOX ", idtmed[230], dal1, 3);
  
  //--- Define the part of the (smaller) rib between two sensitive parts made of
  //     kapton (layer #3) 
  
  dkap[0] = .5 - dits[0];
  dkap[1] = .01585;
  dkap[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR13", "BOX ", idtmed[236], dkap, 3);
  
  //--- Define the second part of the (smaller) rib between two sensitive parts
  //     made of aluminum (layer #3) 
  
  dal2[0] = .5 - dits[0];
  dal2[1] = 0.0027/2.;
  dal2[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR14", "BOX ", idtmed[230], dal2, 3);
  
  // --- Define the part of the (smaller) rib between two sensitive parts 
  //     made of silicon (the electronics) (layer #3) 
  
  dchi[0] = .5 - dits[0];
  dchi[1] = 0.0071/2.;
  dchi[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR15", "BOX ", idtmed[225], dal2, 3);
  
  // --- Define the part of the (smaller) rib between two sensitive parts 
  //     made of water (the cooler) (layer #3) 
  
  dwat[0] = .5 - dits[0];
  dwat[1] = 0.0093/2.;
  dwat[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR16", "BOX ", idtmed[231], dwat, 3);
  
  //--- Define the third part of the (smaller) rib between two sensitive parts
  //     made of aluminum (the cooling tubes) (layer #3) 
  
  dtub[0] = .5 - dits[0];
  dtub[1] = 0.00134/2.;
  dtub[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR17", "BOX ", idtmed[230], dtub, 3);
  
  // --- Define the part of the end-ladder stuff made of PCB (layer #3) 
  
  dpcb[0] = .03;
  // twice the foreseen thickness 
  dpcb[1] = 3.5;
  dpcb[2] = 7.5;
  gMC->Gsvolu("IEL1", "BOX ", idtmed[233], dpcb, 3);
  
  // --- Define the part of the end-ladder stuff made of copper (layer #3) 
  
  dcop[0] = .0252;
  // twice the foreseen thickness 
  dcop[1] = 3.5;
  dcop[2] = 7.5;
  gMC->Gsvolu("IEL2", "BOX ", idtmed[234], dcop, 3);
  
  // --- Define the part of the end-ladder stuff made of ceramics (layer #3)
  
  dcer[0] = .04;
  // twice the foreseen thickness 
  dcer[1] = 3.5;
  dcer[2] = 7.5;
  gMC->Gsvolu("IEL3", "BOX ", idtmed[235], dcer, 3);
  
  // --- Define the part of the end-ladder stuff made of silicon (layer #3) 
  
  dsil[0] = .003;
  // twice the foreseen thickness 
  dsil[1] = 3.5;
  dsil[2] = 7.5;
  gMC->Gsvolu("IEL4", "BOX ", idtmed[226], dsil, 3);
  
  //--- Place the sensitive part of the drifts (smaller ribs) into its mother
  //     (IDV1) 
  
  ypos = 0.;
  for (j = 1; j <= 5; ++j) {
    // odd elements are up and even elements are down 
    if (j == 1) {
      xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0. - dits[2] + 1. - dits[2] * 2. - .1 - dits[2];
    } else if (j == 2) {
      xpos = -dbox1[0] + dits[0];
      zpos = 0. - dits[2] + 1. - dits[2];
    } else if (j == 3) {
      xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0.;
    } else if (j == 4) {
      xpos = -dbox1[0] + dits[0];
      zpos = dits[2] + 0. - 1. + dits[2];
    } else if (j == 5) {
      xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = dits[2] + 0. - 1. + dits[2] * 2. + .1 + dits[2];
    }
    gMC->Gspos("ITS3", j, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  }
  
  // --- Place the smaller ribs into their mother (IDV1) 
  
  // --- Right ribs (just a matter of convention) 
  
  xpos = .5 - dbox1[0] + dits[0];
  zpos = 0.;
  
  // --- Carbon 
  
  ypos = 2.81;
  gMC->Gspos("IR11", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = dsup[1] + 2.81 + dal1[1];
  gMC->Gspos("IR12", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1];
  gMC->Gspos("IR13", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 

  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1];
  gMC->Gspos("IR14", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1];
  gMC->Gspos("IR15", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1];
  gMC->Gspos("IR16", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 2. 
    + dtub[1];
  gMC->Gspos("IR17", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Right ribs (just a matter of convention) 
  
  // --- Carbon 
  
  ypos = -2.81;
  gMC->Gspos("IR11", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = -(dsup[1] + 2.81 + dal1[1]);
  gMC->Gspos("IR12", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1]);
  gMC->Gspos("IR13", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1]);
  gMC->Gspos("IR14", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1]);
  gMC->Gspos("IR15", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1]);
  gMC->Gspos("IR16", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 2. + dtub[1]);
  gMC->Gspos("IR17", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the end-ladder stuff into its mother (IDV1) 
  
  
  // --- Negative-Z end-ladder 
  
  ypos = 0.;
  zpos = -(8.7*5.-2.*1.+2.*0.1)/2.-7.5;
  
  // --- PCB 
  
  xpos = dbox1[0] - dpcb[0];
  gMC->Gspos("IEL1", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL2", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL3", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL4", 1, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Positive-Z end-ladder 
  
  ypos = 0.;
  zpos = (8.7*5.-2.*1.+2.*0.1)/2.+7.5;
  
  // --- PCB 
  
  xpos = dbox1[0] - dpcb[0];
  gMC->Gspos("IEL1", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 

  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL2", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL3", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL4", 2, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  //--- Define a ghost volume containing a single ladder of layer #3 (with the
  //     larger lenght of ribs) and fill it with air or vacuum 
  
  dbox2[0] = 0.65+(0.0172+0.03+0.0252+0.04+0.003);
  dbox2[1] = 3.85;
  // the widest element is the sensitive element 
  dbox2[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IDV2", "BOX ", idtmed[228], dbox2, 3);
  
  // --- Make the ghost volume invisible 
  
  gMC->Gsatt("IDV2", "SEEN", 0);
  
  //--- Define the part of the (larger) rib between two sensitive parts madeof
  //     carbon (layer #3) 
  
  dsup[0] = .65 - dits[0];
  dsup[1] = .01;
  dsup[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR21", "BOX ", idtmed[227], dsup, 3);
  
  //--- Define the first part of the (larger) rib between two sensitive parts
  //     made of aluminum (layer #3) 
  
  dal1[0] = .65 - dits[0];
  dal1[1] = 0.00096/2.;
  dal1[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR22", "BOX ", idtmed[230], dal1, 3);
  
  //--- Define the part of the (larger) rib between two sensitive parts madeof
  //     kapton (layer #3) 
  
  dkap[0] = .65 - dits[0];
  dkap[1] = 0.0317/2.;
  dkap[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR23", "BOX ", idtmed[236], dkap, 3);
  
  //--- Define the second part of the (larger) rib between two sensitive parts
  //     made of aluminum (layer #3) 
  
  dal2[0] = .65 - dits[0];
  dal2[1] = 0.0027/2.;
  dal2[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR24", "BOX ", idtmed[230], dal2, 3);
  
  // --- Define the part of the (larger) rib between two sensitive parts 
  //     made of silicon (the electronics) (layer #3) 

  dchi[0] = .65 - dits[0];
  dchi[1] = 0.0071/2.;
  dchi[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR25", "BOX ", idtmed[225], dal2, 3);
  
  // --- Define the part of the (larger) rib between two sensitive parts 
  //     made of water (the cooler) (layer #3) 
  
  dwat[0] = .65 - dits[0];
  dwat[1] = 0.0093/2.;
  dwat[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR26", "BOX ", idtmed[231], dwat, 3);
  
  //--- Define the third part of the (larger) rib between two sensitive parts
  //     made of aluminum (the cooling tubes) (layer #3) 
  
  dtub[0] = .65 - dits[0];
  dtub[1] = 0.00134/2.;
  dtub[2] = (8.7*5.-2.*1.+2.*0.1)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR27", "BOX ", idtmed[230], dtub, 3);
  
  //--- Place the sensitive part of the drifts (smaller ribs) into its mother
  //     (IDV2) 
  
  ypos = 0.;
  for (j = 1; j <= 5; ++j) {
    // odd element are up and even elements are down 
    if (j == 1) {
      xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0. - dits[2] + 1. - dits[2] * 2. - .1 - dits[2];
    } else if (j == 2) {
      xpos = -dbox2[0] + dits[0];
      zpos = 0. - dits[2] + 1. - dits[2];
    } else if (j == 3) {
      xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0.;
    } else if (j == 4) {
      xpos = -dbox2[0] + dits[0];
      zpos = dits[2] + 0. - 1. + dits[2];
    } else if (j == 5) {
      xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = dits[2] + 0. - 1. + dits[2] * 2. + .1 + dits[2];
    }
    gMC->Gspos("ITS3", j, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  }
  
  // --- Place the larger ribs into their mother (IDV2) 
  
  
  // --- Right ribs (just a matter of convention) 
  
  xpos = .65 - dbox2[0] + dits[0];
  zpos = 0.;
  
  // --- Carbon 
  
  ypos = 2.81;
  gMC->Gspos("IR21", 1, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = dsup[1] + 2.81 + dal1[1];
  gMC->Gspos("IR22", 1, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1];
  gMC->Gspos("IR23", 1, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1];
  gMC->Gspos("IR24", 1, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1];
  gMC->Gspos("IR25", 1, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1];
  gMC->Gspos("IR26", 1, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 2. + dtub[1];
  gMC->Gspos("IR27", 1, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Right ribs (just a matter of convention) 
  
  // --- Carbon 
  
  ypos = -2.81;
  gMC->Gspos("IR21", 2, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = -(dsup[1] + 2.81 + dal1[1]);
  gMC->Gspos("IR22", 2, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1]);
  gMC->Gspos("IR23", 2, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1]);
  gMC->Gspos("IR24", 2, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1]);
  gMC->Gspos("IR25", 2, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1]);
  gMC->Gspos("IR26", 2, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 2. + dtub[1]);
  gMC->Gspos("IR27", 2, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the end-ladder stuff into its mother (IDV1) 
  
  
  // --- Negative-Z end-ladder 
  
  ypos = 0.;
  zpos = -(8.7*5.-2.*1.+2.*0.1)/2.-7.5;
  
  // --- PCB 
  
  xpos = dbox2[0] - dpcb[0];
  gMC->Gspos("IEL1", 3, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL2", 3, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL3", 3, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL4", 3, "IDV1", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Positive-Z end-ladder 
  
  //yos  = 0.;
  zpos = (8.7*5.-2.*1.+2.*0.1)/2.+7.5;
  
  // --- PCB 
  
  xpos = dbox2[0] - dpcb[0];
  gMC->Gspos("IEL1", 4, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL2", 4, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL3", 4, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL4", 4, "IDV2", xpos, ypos, zpos, 0, "ONLY");
  
  //--- Place the ghost volumes containing the drift ladders of layer #3 in their
  //     mother volume (IT34) 
  //     Odd elements have large ribs and even elements have small ribs 
  
  for (i = 1; i <= 12; ++i) {
    atheta = (i-1) * 30.;
    AliMatrix(idrotm[i+1299], 90., atheta, 90., atheta + 90., 0.,0.);
    if (i % 2 == 0) {
      rzero = 14.;
      xpos = rzero * TMath::Cos((i-1) * twopi / 12.);
      ypos = rzero * TMath::Sin((i-1) * twopi / 12.);
      zpos = 0.;
      gMC->Gspos("IDV1", i, "IT34", xpos, ypos, zpos, idrotm[i+1299], "ONLY");
    } else {
      rzero = 13.85;
      xpos = rzero * TMath::Cos((i-1) * twopi / 12.);
      ypos = rzero * TMath::Sin((i-1) * twopi / 12.);
      zpos = 0.;
      gMC->Gspos("IDV2", i, "IT34", xpos, ypos, zpos, idrotm[i+1299], "ONLY");
    }
  }
  
  
  // --- Layer #4 
  
  //        GOTO 4567           ! skip ITS layer no. 4 
  
  //--- Define a ghost volume containing a single ladder of layer #4 (with the
  //     smaller lenght of ribs) and fill it with air or vacuum 
  
  dbox1[0] = 0.5+(0.0172+0.03+0.0252+0.04+0.003);
  dbox1[1] = 3.5;
  // the widest element is the end-ladder stuff 
  dbox1[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IDV3", "BOX ", idtmed[228], dbox1, 3);
  
  // --- Make the ghost volume invisible 
  
  gMC->Gsatt("IDV3", "SEEN", 0);
  
  // --- Define a volume containing the sensitive part of drifts 
  //     (silicon, layer #4) 
  
  dits[0] = .0172;
  // see material budget report by G. Feofilov 
  dits[1] = 3.125;
  dits[2] = 4.35;
  gMC->Gsvolu("ITS4", "BOX ", idtmed[224], dits, 3);
  
  //--- Define the part of the (smaller) rib between two sensitive parts made of
  //     carbon (layer #4) 
  
  dsup[0] = .5 - dits[0];
  dsup[1] = .01;
  dsup[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR31", "BOX ", idtmed[227], dsup, 3);
  
  //--- Define the first part of the (smaller) rib between two sensitive parts
  //     made of aluminum (layer #4) 
  
  dal1[0] = .5 - dits[0];
  dal1[1] = 0.00096/2.;
  dal1[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR32", "BOX ", idtmed[230], dal1, 3);
  
  //--- Define the part of the (smaller) rib between two sensitive parts made of
  //     kapton (layer #4) 
  
  dkap[0] = .5 - dits[0];
  dkap[1] = 0.0317/2.;
  dkap[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR33", "BOX ", idtmed[236], dkap, 3);
  
  //--- Define the second part of the (smaller) rib between two sensitive parts
  //     made of aluminum (layer #4) 
  
  dal2[0] = .5 - dits[0];
  dal2[1] = 0.0027/2.;
  dal2[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR34", "BOX ", idtmed[230], dal2, 3);
  
  // --- Define the part of the (smaller) rib between two sensitive parts 
  //     made of silicon (the electronics) (layer #4) 
  
  dchi[0] = .5 - dits[0];
  dchi[1] = 0.0071/2.;
  dchi[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR35", "BOX ", idtmed[225], dal2, 3);
  
  // --- Define the part of the (smaller) rib between two sensitive parts 
  //     made of water (the cooler) (layer #4) 
  
  dwat[0] = .5 - dits[0];
  dwat[1] = 0.0093/2.;
  dwat[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IR36", "BOX ", idtmed[231], dwat, 3);
  
  //--- Define the third part of the (smaller) rib between two sensitive parts
  //     made of aluminum (the cooling tubes) (layer #4) 
  
  dtub[0] = .5 - dits[0];
  dtub[1] = 0.00134/2.;
  dtub[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR37", "BOX ", idtmed[230], dtub, 3);
  
  // --- Define the part of the end-ladder stuff made of PCB (layer #4) 
  
  dpcb[0] = .03;
  // twice the foreseen thickness 
  dpcb[1] = 3.5;
  dpcb[2] = 7.5;
  gMC->Gsvolu("IEL5", "BOX ", idtmed[233], dpcb, 3);
  
  // --- Define the part of the end-ladder stuff made of copper (layer #4) 
  
  dcop[0] = .0252;
  // twice the foreseen thickness 
  dcop[1] = 3.5;
  dcop[2] = 7.5;
  gMC->Gsvolu("IEL6", "BOX ", idtmed[234], dcop, 3);
  
  // --- Define the part of the end-ladder stuff made of ceramics (layer #4)
  
  dcer[0] = .04;
  // twice the foreseen thickness 
  dcer[1] = 3.5;
  dcer[2] = 7.5;
  gMC->Gsvolu("IEL7", "BOX ", idtmed[235], dcer, 3);
  
  // --- Define the part of the end-ladder stuff made of silicon (layer #4) 
  
  dsil[0] = .003;
  // twice the foreseen thickness 
  dsil[1] = 3.5;
  dsil[2] = 7.5;
  gMC->Gsvolu("IEL8", "BOX ", idtmed[226], dsil, 3);
  
  //--- Place the sensitive part of the drifts (smaller ribs) into its mother
  //     (IDV3) 
  
  ypos = 0.;
  for (j = 1; j <= 7; ++j) {
    // odd elements are down and even elements are up 
    if (j == 1) {
      xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0. - dits[2] + .7 - dits[2] * 2. + 0. - dits[2] * 2. + 1.3 - dits[2];
    } else if (j == 2) {
      xpos = -dbox1[0] + dits[0];
      zpos = 0. - dits[2] + .7 - dits[2] * 2. + 0. - dits[2];
    } else if (j == 3) {
      xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0. - dits[2] + .7 - dits[2];
    } else if (j == 4) {
      xpos = -dbox1[0] + dits[0];
      zpos = 0.;
    } else if (j == 5) {
      xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = dits[2] + 0. - .7 + dits[2];
    } else if (j == 6) {
      xpos = -dbox1[0] + dits[0];
      zpos = dits[2] + 0. - .7 + dits[2] * 2. + 0. + dits[2];
    } else if (j == 7) {
      xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = dits[2] + 0. - .7 + dits[2] * 2. + 0. + dits[2] * 2. - 1.3 + dits[2];
    }
    gMC->Gspos("ITS4", j, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  }
  
  // --- Place the smaller ribs into their mother (IDV3) 
  
  // --- Right ribs (just a matter of convention) 
  
  xpos = .5 - dbox1[0] + dits[0];
  zpos = 0.;
  
  // --- Carbon 
  
  ypos = 2.81;
  gMC->Gspos("IR31", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = dsup[1] + 2.81 + dal1[1];
  gMC->Gspos("IR32", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 

  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1];
  gMC->Gspos("IR33", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1];
  gMC->Gspos("IR34", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1];
  gMC->Gspos("IR35", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1];
  gMC->Gspos("IR36", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 2. 
    + dtub[1];
  gMC->Gspos("IR37", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Right ribs (just a matter of convention) 
  
  // --- Carbon 
  
  ypos = -2.81;
  gMC->Gspos("IR31", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = -(dsup[1] + 2.81 + dal1[1]);
  gMC->Gspos("IR32", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1]);
  gMC->Gspos("IR33", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1]);
  gMC->Gspos("IR34", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1] * 2. + dchi[1]);
  gMC->Gspos("IR35", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1]);
  gMC->Gspos("IR36", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 
	   2. + dtub[1]);
  gMC->Gspos("IR37", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the end-ladder stuff into its mother (IDV1) 
  
  
  // --- Negative-Z end-ladder 
  
  ypos = 0.;
  zpos = -(8.7*7.-2.*0.7-2.*1.3)/2.-7.5;
  
  // --- PCB 
  
  xpos = dbox1[0] - dpcb[0];
  gMC->Gspos("IEL5", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL6", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL7", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL8", 1, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Positive-Z end-ladder 
  
  ypos = 0.;
  zpos = (8.7*7.-2.*0.7-2.*1.3)/2.-7.5;
  
  // --- PCB 
  
  xpos = dbox1[0] - dpcb[0];
  gMC->Gspos("IEL5", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL6", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL7", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox1[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL8", 2, "IDV3", xpos, ypos, zpos, 0, "ONLY");
  
  //--- Define a ghost volume containing a single ladder of layer #4 (with the
  //     larger lenght of ribs) and fill it with air or vacuum 
  
  dbox2[0] = 0.65+(0.0172+0.03+0.0252+0.04+0.003);
  dbox2[1] = 3.5;
  // the widest element is the end-ladder stuff 
  dbox2[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lenght 
  gMC->Gsvolu("IDV4", "BOX ", idtmed[228], dbox2, 3);
  
  // --- Make the ghost volume invisible 
  
  gMC->Gsatt("IDV4", "SEEN", 0);
  
  //--- Define the part of the (larger) rib between two sensitive parts madeof
  //     carbon (layer #4) 
  
  dsup[0] = .65 - dits[0];
  dsup[1] = .01;
  dsup[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR41", "BOX ", idtmed[227], dsup, 3);
  
  //--- Define the first part of the (larger) rib between two sensitive parts
  //     made of aluminum (layer #4) 
  
  dal1[0] = .65 - dits[0];
  dal1[1] = 0.00096/2.;
  dal1[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR42", "BOX ", idtmed[230], dal1, 3);
  
  //--- Define the part of the (larger) rib between two sensitive parts madeof
  //     kapton (layer #4) 
  
  dkap[0] = .65 - dits[0];
  dkap[1] = 0.0317/2.;
  dkap[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR43", "BOX ", idtmed[236], dkap, 3);
  
  //--- Define the second part of the (larger) rib between two sensitive parts
  //     made of aluminum (layer #4) 
  
  dal2[0] = .65 - dits[0];
  dal2[1] = 0.0027/2.;
  dal2[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR44", "BOX ", idtmed[230], dal2, 3);
  
  // --- Define the part of the (larger) rib between two sensitive parts 
  //     made of silicon (the electronics) (layer #4) 
  
  dchi[0] = .65 - dits[0];
  dchi[1] = 0.0071/2.;
  dchi[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR45", "BOX ", idtmed[225], dal2, 3);
  
  // --- Define the part of the (larger) rib between two sensitive parts 
  //     made of water (the cooler) (layer #4) 
  
  dwat[0] = .65 - dits[0];
  dwat[1] = 0.0093/2.;
  dwat[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR46", "BOX ", idtmed[231], dwat, 3);
  
  //--- Define the third part of the (larger) rib between two sensitive parts
  //     made of aluminum (the cooling tubes) (layer #4) 
  
  dtub[0] = .65 - dits[0];
  dtub[1] = 0.00134/2.;
  dtub[2] = (8.7*7.-2.*0.7-2.*1.3)/2.+2.*7.5;
  // 7.5 cm is the lengh 
  gMC->Gsvolu("IR47", "BOX ", idtmed[230], dtub, 3);
  
  //--- Place the sensitive part of the drifts (smaller ribs) into its mother
  //     (IDV4) 
  
  ypos = 0.;
  for (j = 1; j <= 7; ++j) {
    // odd elements are down and even elements are up 
    if (j == 1) {
      xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0. - dits[2] + .7 - dits[2] * 2. + 0. - dits[2] * 2. + 1.3 - dits[2];
    } else if (j == 2) {
      xpos = -dbox2[0] + dits[0];
      zpos = 0. - dits[2] + .7 - dits[2] * 2. + 0. - dits[2];
    } else if (j == 3) {
      xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = 0. - dits[2] + .7 - dits[2];
    } else if (j == 4) {
      xpos = -dbox2[0] + dits[0];
      zpos = 0.;
    } else if (j == 5) {
      xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = dits[2] + 0. - .7 + dits[2];
    } else if (j == 6) {
      xpos = -dbox2[0] + dits[0];
      zpos = dits[2] + 0. - .7 + dits[2] * 2. + 0. + dits[2];
    } else if (j == 7) {
      xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0] * 2. - dits[0];
      zpos = dits[2] + 0. - .7 + dits[2] * 2. + 0. + dits[2] * 2. - 1.3 + dits[2];
    }
    gMC->Gspos("ITS4", j, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  }
  
  // --- Place the larger ribs into their mother (IDV4) 
  
  
  // --- Right ribs (just a matter of convention) 
  
  xpos = .65 - dbox2[0] + dits[0];
  zpos = 0.;
  
  // --- Carbon 
  
  ypos = 2.81;
  gMC->Gspos("IR41", 1, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = dsup[1] + 2.81 + dal1[1];
  gMC->Gspos("IR42", 1, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1];
  gMC->Gspos("IR43", 1, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1];
  gMC->Gspos("IR44", 1, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1];
  gMC->Gspos("IR45", 1, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1];
  gMC->Gspos("IR46", 1, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 2. 
    + dtub[1];
  gMC->Gspos("IR47", 1, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Right ribs (just a matter of convention) 
  
  // --- Carbon 
  
  ypos = -2.81;
  gMC->Gspos("IR41", 2, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #1 
  
  ypos = -(dsup[1] + 2.81 + dal1[1]);
  gMC->Gspos("IR42", 2, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Kapton 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1]);
  gMC->Gspos("IR43", 2, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #2 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1]);
  gMC->Gspos("IR44", 2, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (chip) 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1] * 2. + dchi[1]);
  gMC->Gspos("IR45", 2, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Water 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1]);
  gMC->Gspos("IR46", 2, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Aluminum #3 
  
  ypos = -(dsup[1] + 2.81 + dal1[1] * 2. + dkap[1] * 
	   2. + dal2[1] * 2. + dchi[1] * 2. + dwat[1] * 2. + dtub[1]);
  gMC->Gspos("IR47", 2, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Place the end-ladder stuff into its mother (IDV1) 
  
  
  // --- Negative-Z end-ladder 
  
  ypos = 0.;
  zpos = -(8.7*7.-2.*0.7-2.*1.3)/2.-7.5;
  
  // --- PCB 
  
  xpos = dbox2[0] - dpcb[0];
  gMC->Gspos("IEL5", 3, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL6", 3, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL7", 3, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL8", 3, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Positive-Z end-ladder 
  
  //yos  = 0.;
  zpos = (8.7*7.-2.*0.7-2.*1.3)/2.-7.5;
  
  // --- PCB 
  
  xpos = dbox2[0] - dpcb[0];
  gMC->Gspos("IEL5", 4, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Copper 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0];
  gMC->Gspos("IEL6", 4, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Ceramics 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0];
  gMC->Gspos("IEL7", 4, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  // --- Silicon (bus) 
  
  xpos = dbox2[0] - dpcb[0] * 2. - dcop[0] * 2. - dcer[0] * 2. - dsil[0];
  gMC->Gspos("IEL8", 4, "IDV4", xpos, ypos, zpos, 0, "ONLY");
  
  //--- Place the ghost volumes containing the drift ladders of layer #4 in their
  //     mother volume (IT34) 
  //     Odd elements have large ribs and even elements have small ribs 
  
  for (i = 1; i <= 24; ++i) {
    atheta = (i-1) * 15.;
    AliMatrix(idrotm[i+1399], 90., atheta, 90., atheta + 90., 0.,0.);
    if (i % 2 == 0) {
      rzero = 23.5;
      xpos = rzero * TMath::Cos((i-1) * twopi / 24.);
      ypos = rzero * TMath::Sin((i-1) * twopi / 24.);
      zpos = 0.;
      gMC->Gspos("IDV3", i, "IT34", xpos, ypos, zpos, idrotm[i+1399], "ONLY");
    } else {
      rzero = (24.0+22.8)/2.;
      xpos = rzero * TMath::Cos((i-1) * twopi / 24.);
      ypos = rzero * TMath::Sin((i-1) * twopi / 24.);
      zpos = 0.;
      gMC->Gspos("IDV4", i, "IT34", xpos, ypos, zpos, idrotm[i+1399], "ONLY");
    }
  }
  
  //************************************************************************
  //*                                                                      *
  //*                               S T R I P S                            *
  //*                               ===========                            *
  //*                                                                      *
  //************************************************************************
  
  // --- Define SSD with the 35+39 lay-out 
  
  if (fMinorVersionV3 < 3) {
    
    //--- Define ghost volume containing the Strip Detectors and fill it with air
    //     or vacuum 
    
    xxm    = (49.999-3.)/(70.-25.);
    dgh[0] = 0.;
    dgh[1] = 360.;
    dgh[2] = 4.;
    dgh[3] = -25.-(9.-3.01)/xxm-(9.01-9.)/xxm-(27.-9.01)/xxm-
      (37.-27)/xxm-(49.998-37.)/xxm;
    dgh[4] = 49.998;
    dgh[5] = 49.998;
    dgh[6] = -25.-(9.-3.01)/xxm-(9.01-9.)/xxm-(27.-9.01)/xxm-
      (37.-27)/xxm;
    dgh[7] = 37.;
    dgh[8] = 49.998;
    dgh[9] = 25.+(9.-3.01)/xxm+(9.01-9.)/xxm+(27.-9.01)/xxm+
      (37.-27)/xxm;
    dgh[10] = 37.;
    dgh[11] = 49.998;
    dgh[12] = 25.+(9.-3.01)/xxm+(9.01-9.)/xxm+(27.-9.01)/xxm+
      (37.-27)/xxm+(49.998-37.)/xxm;
    dgh[13] = 49.998;
    dgh[14] = 49.998;
    gMC->Gsvolu("IT56", "PCON", idtmed[275], dgh, 15);
    gMC->Gspos("IT56", 1, "ITSV", 0., 0., 0., 0, "ONLY");
    gMC->Gsatt("IT56", "SEEN", 0);
    
    // --- Layer #5 
    
    //        GOTO 5678           ! skip ITS layer no. 5 
    
    //--- Define a ghost volume containing a single ladder of layer #5 andfill
    //     it with air or vacuum 
    
    dbox1[0] = (0.0600+2.*0.0150)/2.;
    dbox1[1] = 3.75;
    dbox1[2] = 90.22/2.;
    gMC->Gsvolu("ISV1", "BOX ", idtmed[253], dbox1, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ISV1", "SEEN", 0);
    
    // --- Define a ghost volume containing the electronics and cooling of
    //     a single ladder of layer #5 and fill it with air or vacuum 
    
    dsrv[0] = (TMath::Sqrt(3.) / 2. * 4.2 + .47 + .05) / 2.;
    dsrv[1] = 3.75;
    dsrv[2] = 90.22/2.;
    gMC->Gsvolu("SSV1", "BOX ", idtmed[253], dsrv, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("SSV1", "SEEN", 0);
    
    // --- Define a ghost volume containing the end-ladder stuff of 
    //     a single ladder of layer #5 and fill it with air or vacuum 
    
    dela[0] = 2.;
    dela[1] = 3.5;
    dela[2] = 4.;
    gMC->Gsvolu("ELL5", "BOX ", idtmed[253], dela, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ELL5", "SEEN", 0);
    
    // --- Define a volume containing the sensitive part of the strips 
    //     (silicon, layer #5) 
    
    dits[0] = .015;
    dits[1] = 3.75;
    dits[2] = 2.1;
    gMC->Gsvolu("ITS5", "BOX ", idtmed[249], dits, 3);
    
    // --- Define a volume containing the electronics of the strips 
    //     (silicon, layer #5) 
    
    dchi[0] = .02;
    dchi[1] = 3.4;
    dchi[2] = .525;
    gMC->Gsvolu("SCH5", "BOX ", idtmed[250], dchi, 3);
    
    // --- Define the cooling tubes (aluminum, layer #5) 
    
    dtub[0] = .09;
    dtub[1] = dtub[0] + .01;
    dtub[2] = 90.22/2.;
    gMC->Gsvolu("STB5", "TUBE", idtmed[255], dtub, 3);
    
    // --- Define the cooling fluid (water or freon, layer #5) 
    
    dwat[0] = 0.;
    dwat[1] = .09;
    dwat[2] = 90.22/2.;
    gMC->Gsvolu("SWT5", "TUBE", idtmed[256], dwat, 3);
    //        CALL GSVOLU('SWT5','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
    
    //--- Define the (triangular) element of the heat bridge (carbon, layer #5)
    
    // water 
    dfra[0] = 120.;
    dfra[1] = 360.;
    dfra[2] = 3.;
    dfra[3] = 2.;
    dfra[4] = -.015;
    dfra[5] = TMath::Sqrt(3.) * 4.2 / 6.;
    dfra[6] = dfra[5] + .03;
    dfra[7] = .015;
    dfra[8] = dfra[5];
    dfra[9] = dfra[6];
    gMC->Gsvolu("SFR5", "PGON", idtmed[252], dfra, 10);
    
    // --- Define the element connecting the triangles of the heat bridge 
    //     (carbon, layer #5) 
    
    dcei[0] = 0.;
    dcei[1] = .03;
    dcei[2] = 90.22/2.;
    gMC->Gsvolu("SCE5", "TUBE", idtmed[252], dcei, 3);
    
    // --- Define the part of the end-ladder stuff made of plastic (G10FR4) 
    //     (layer #5) 
    
    dpla[0] = (10./(8.*7.))/2.;
    dpla[1] = 3.5;
    dpla[2] = 4.;
    gMC->Gsvolu("EPL5", "BOX ", idtmed[262], dpla, 3);
    
    // --- Define the part of the end-ladder stuff made of copper (layer #5) 
    
    dcop[0] = (2./(8.*7.))/2.;
    dcop[1] = 3.5;
    dcop[2] = 4.;
    gMC->Gsvolu("ECU5", "BOX ", idtmed[259], dcop, 3);
    
    // --- Define the part of the end-ladder stuff made of epoxy (layer #5) 
    
    depx[0] = (30./(8.*7.))/2.;
    depx[1] = 3.5;
    depx[2] = 4.;
    gMC->Gsvolu("EPX5", "BOX ", idtmed[262], depx, 3);
    
    // --- Define the part of the end-ladder stuff made of silicon (bus) 
    //     (layer #5) 
    
    dsil[0] = (20./(8.*7.))/2.;
    dsil[1] = 3.5;
    dsil[2] = 4.;
    gMC->Gsvolu("ESI5", "BOX ", idtmed[251], dsil, 3);
    
    // --- Place the end-ladder stuff into its mother (ELL5) 
    
    sep = (4. - (dpla[0] + dcop[0] + depx[0] + dsil[0]) * 2.) / 3.;
    ypos = 0.;
    zpos = 0.;
    
    // --- Plastic 
    
    xpos = -dela[0] + dpla[0];
    gMC->Gspos("EPL5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Copper 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0];
    gMC->Gspos("ECU5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Epoxy 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0];
    gMC->Gspos("EPX5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Silicon (bus) 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0] * 2. + sep + dsil[0];
    gMC->Gspos("ESI5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the sensitive part of the strips into its mother (ISV1) 
    
    ypos = 0.;
    for (j = 1; j <= 23; ++j) {
      if (j % 2 == 0) xpos = dbox1[0] - dits[0];
      else            xpos = -dbox1[0] + dits[0];
      zpos = ((j - 1) - 11.) * 3.91;
      gMC->Gspos("ITS5", j, "ISV1", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the electronics of the strips into its mother (SSV1) 
    
    ypos = 0.;
    for (j = 1; j <= 23; ++j) {
      if (j % 2 == 0) xpos = -dsrv[0] + .28;
      else            xpos = -dsrv[0] + .28 - dits[0] * 2. - .03;
      zpos = ((j - 1) - 11.) * 3.91 + .85;
      gMC->Gspos("SCH5", j, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    }
    
    //--- Place the cooling tubes and the cooling fluid into their mother (SSV1)
    
    xpos = -dsrv[0] + .41;
    zpos = 0.;
    
    // --- Left tube (just a matter of convention) 
    
    ypos = -2.25-0.1;
    gMC->Gspos("STB5", 1, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT5", 1, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right tube (just a matter of convention) 
    
    ypos = 2.25+0.1;
    gMC->Gspos("STB5", 2, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT5", 2, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the heat bridge elements into their mother (SSV1) 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 6. * 4.2;
    ypos = 0.;
    for (j = 1; j <= 24; ++j) {
      zpos = ((j - 1) - 11.) * 3.91 - -4.2/2.;
      gMC->Gspos("SFR5", j, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the elements connecting the triangles of the heat bridge 
    //     into their mother (SSV1) 
    
    zpos = 0.;
    
    // --- Left element (just a matter of convention) 
    
    xpos = -dsrv[0] + .47;
    ypos = -(2.1+0.015);
    gMC->Gspos("SCE5", 1, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right element 
    
    xpos = -dsrv[0] + .47;
    ypos = 2.1+0.015;
    gMC->Gspos("SCE5", 2, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Top element 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 2. * 4.2 + .015;
    ypos = 0.;
    gMC->Gspos("SCE5", 3, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the ghost volumes containing the strip ladders (ISV1), 
    //    electronics/cooling (SSV1) and end-ladder stuff (ELL5) of layer #5 in
    //     their mother volume (IT56) 
    
    offset1 = TMath::ATan2(.9, 40.);
    offset2 = 5.2;
    rzero   = dbox1[0] + 40.;
    runo    = dbox1[0] * 2. + 40. + dsrv[0];
    rtwo    = dbox1[0] * 2. + 40. + dela[0];
    for (i = 1; i <= 35; ++i) {
      atheta = (i-1) * twopi * raddeg / 35. + offset2;
      AliMatrix(idrotm[i+1499], 90., atheta, 90., atheta + 90., 0., 0.);
      
      // --- Strip ladders 
      
      xpos = rzero * TMath::Cos((i-1) * twopi / 35. + offset1);
      ypos = rzero * TMath::Sin((i-1) * twopi / 35. + offset1);
      zpos = 0.;
      gMC->Gspos("ISV1", i, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
      
      // --- Electronics/cooling 
      
      xpos = runo * TMath::Cos((i-1) * twopi / 35. + offset1);
      ypos = runo * TMath::Sin((i-1) * twopi / 35. + offset1);
      zpos = 0.;
      gMC->Gspos("SSV1", i, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
      
      // --- End-ladders (nagative-Z and positive-Z) 
      
      xpos = rtwo * TMath::Cos((i-1) * twopi / 35. + offset1);
      ypos = rtwo * TMath::Sin((i-1) * twopi / 35. + offset1);
      zpos = -(dbox1[2] + dela[2] + 6.);
      gMC->Gspos("ELL5", i, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
      zpos = dbox1[2] + dela[2] + 6.;
      gMC->Gspos("ELL5", i + 35, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
    }
    
    
    // --- Layer #6 
    
    //        GOTO 5778           ! skip ITS layer no. 6 
    
    //--- Define a ghost volume containing a single ladder of layer #6 andfill
    //     it with air or vacuum 
    
    dbox2[0] = (0.0600+2.*0.0150)/2.;
    dbox2[1] = 3.75;
    dbox2[2] = 101.95/2.;
    gMC->Gsvolu("ISV2", "BOX ", idtmed[253], dbox2, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ISV2", "SEEN", 0);
    
    // --- Define a ghost volume containing the electronics and cooling of
    //     a single ladder of layer #6 and fill it with air or vacuum 
    
    dsrv[0] = (TMath::Sqrt(3.) / 2. * 4.2 + .47 + .05) / 2.;
    dsrv[1] = 3.75;
    dsrv[2] = 101.95/2.;
    gMC->Gsvolu("SSV2", "BOX ", idtmed[253], dsrv, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("SSV2", "SEEN", 0);
    
    // --- Define a ghost volume containing the end-ladder stuff of 
    //     a single ladder of layer #6 and fill it with air or vacuum 
    
    dela[0] = 2.;
    dela[1] = 3.5;
    dela[2] = 4.;
    gMC->Gsvolu("ELL6", "BOX ", idtmed[253], dela, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ELL6", "SEEN", 0);
    
    // --- Define a volume containing the sensitive part of the strips 
    //     (silicon, layer #6) 
    
    dits[0] = .015;
    dits[1] = 3.75;
    dits[2] = 2.1;
    gMC->Gsvolu("ITS6", "BOX ", idtmed[249], dits, 3);
    
    // --- Define a volume containing the electronics of the strips 
    //     (silicon, layer #6) 
    
    dchi[0] = .02;
    dchi[1] = 3.4;
    dchi[2] = .525;
    gMC->Gsvolu("SCH6", "BOX ", idtmed[250], dchi, 3);
    
    // --- Define the cooling tubes (aluminum, layer #6) 
    
    dtub[0] = .09;
    dtub[1] = dtub[0] + .01;
    dtub[2] = 101.95/2.;
    gMC->Gsvolu("STB6", "TUBE", idtmed[255], dtub, 3);
    
    // --- Define the cooling fluid (water or freon, layer #6) 
    
    dwat[0] = 0.;
    dwat[1] = .09;
    dwat[2] = 101.95/2.;
    gMC->Gsvolu("SWT6", "TUBE", idtmed[256], dwat, 3);
    //        CALL GSVOLU('SWT6','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
    
    //--- Define the (triangular) element of the heat bridge (carbon, layer #6)
    
    // water 
    dfra[0] = 120.;
    dfra[1] = 360.;
    dfra[2] = 3.;
    dfra[3] = 2.;
    dfra[4] = -.015;
    dfra[5] = TMath::Sqrt(3.) * 4.2 / 6.;
    dfra[6] = dfra[5] + .03;
    dfra[7] = .015;
    dfra[8] = dfra[5];
    dfra[9] = dfra[6];
    gMC->Gsvolu("SFR6", "PGON", idtmed[252], dfra, 10);
    
    // --- Define the element connecting the triangles of the heat bridge 
    //     (carbon, layer #6) 
    
    dcei[0] = 0.;
    dcei[1] = .03;
    dcei[2] = 101.95/2.;
    gMC->Gsvolu("SCE6", "TUBE", idtmed[252], dcei, 3);
    
    // --- Define the part of the end-ladder stuff made of plastic (G10FR4) 
    //     (layer #6) 
    
    dpla[0] = (10./(8.*7.))/2.;
    dpla[1] = 3.5;
    dpla[2] = 4.;
    gMC->Gsvolu("EPL6", "BOX ", idtmed[262], dpla, 3);
    
    // --- Define the part of the end-ladder stuff made of copper (layer #6) 
    
    dcop[0] = (2./(8.*7.))/2.;
    dcop[1] = 3.5;
    dcop[2] = 4.;
    gMC->Gsvolu("ECU6", "BOX ", idtmed[259], dcop, 3);
    
    // --- Define the part of the end-ladder stuff made of epoxy (layer #6) 
    
    depx[0] = (30./(8.*7.))/2.;
    depx[1] = 3.5;
    depx[2] = 4.;
    gMC->Gsvolu("EPX6", "BOX ", idtmed[262], depx, 3);
    
    // --- Define the part of the end-ladder stuff made of silicon (bus) 
    //     (layer #6) 
    
    dsil[0] = (20./(8.*7.))/2.;
    dsil[1] = 3.5;
    dsil[2] = 4.;
    gMC->Gsvolu("ESI6", "BOX ", idtmed[251], dsil, 3);
    
    // --- Place the end-ladder stuff into its mother (ELL5) 
    
    sep  = (4. - (dpla[0] + dcop[0] + depx[0] + dsil[0]) * 2.) / 3.;
    ypos = 0.;
    zpos = 0.;
    
    // --- Plastic 
    
    xpos = -dela[0] + dpla[0];
    gMC->Gspos("EPL6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Copper 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0];
    gMC->Gspos("ECU6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Epoxy 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0];
    gMC->Gspos("EPX6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Silicon (bus) 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0] * 2. + sep + dsil[0];
    gMC->Gspos("ESI6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the sensitive part of the strips into its mother (ISV2) 
    
    ypos = 0.;
    for (j = 1; j <= 26; ++j) {
      if (j % 2 == 0) xpos = dbox2[0] - dits[0];
      else            xpos = -dbox2[0] + dits[0];
      zpos = ((j - 1) - 12.) * 3.91 - 1.96;
      gMC->Gspos("ITS6", j, "ISV2", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the electronics of the strips into its mother (SSV2) 
    
    ypos = 0.;
    for (j = 1; j <= 26; ++j) {
      if (j % 2 == 0) xpos = -dsrv[0] + .28;
      else            xpos = -dsrv[0] + .28 - dits[0] * 2. - .03;
      zpos = ((j - 1) - 12.) * 3.91 - 1.96 + .85;
      gMC->Gspos("SCH5", j, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    }

    //--- Place the cooling tubes and the cooling fluid into their mother (SSV1)
    
    xpos = -dsrv[0] + .41;
    zpos = 0.;
    
    // --- Left tube (just a matter of convention) 
    
    ypos = -2.25-0.1;
    gMC->Gspos("STB6", 1, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT6", 1, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right tube (just a matter of convention) 
    
    ypos = 2.25+0.;
    gMC->Gspos("STB6", 2, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT6", 2, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the heat bridge elements into their mother (SSV2) 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 6. * 4.2;
    ypos = 0.;
    for (j = 1; j <= 27; ++j) {
      zpos = ((j - 1) - 12.) * 3.91 - 1.96 - 4.2/2.;
      gMC->Gspos("SFR6", j, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the elements connecting the triangles of the heat bridge 
    //     into their mother (SSV2) 
    
    zpos = 0.;
    
    // --- Left element (just a matter of convention) 
    
    xpos = -dsrv[0] + .47;
    ypos = -(2.1+0.015);
    gMC->Gspos("SCE6", 1, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right element 
    
    xpos = -dsrv[0] + .47;
    ypos = 2.1+0.015;
    gMC->Gspos("SCE6", 2, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Top element 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 2. * 4.2 + .015;
    ypos = 0.;
    gMC->Gspos("SCE6", 3, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the ghost volumes containing the strip ladders (ISV2), 
    //    electronics/cooling (SSV2) and end-ladder stuff (ELL6) of layer #6 in
    //     their mother volume (IT56) 
    
    offset1 = TMath::ATan2(1., 45.);
    offset2 = 5.2;
    rzero   = dbox2[0] + 45.;
    runo    = dbox2[0] * 2. + 45. + dsrv[0];
    rtwo    = dbox2[0] * 2. + 45. + dela[0];
    for (i = 1; i <= 39; ++i) {
      atheta = (i-1) * twopi * raddeg / 39. + offset2;
      AliMatrix(idrotm[i+1599], 90., atheta, 90., atheta + 90., 0., 0.);
      
      // --- Strip ladders 
      
      xpos = rzero * TMath::Cos((i-1) * twopi / 39. + offset1);
      ypos = rzero * TMath::Sin((i-1) * twopi / 39. + offset1);
      zpos = 0.;
      gMC->Gspos("ISV2", i, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
      
      // --- Electronics/cooling 
      
      xpos = runo * TMath::Cos((i-1) * twopi / 39. + offset1);
      ypos = runo * TMath::Sin((i-1) * twopi / 39. + offset1);
      zpos = 0.;
      gMC->Gspos("SSV2", i, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
      
      // --- End-ladders (nagative-Z and positive-Z) 
      
      xpos = rtwo * TMath::Cos((i-1) * twopi / 39. + offset1);
      ypos = rtwo * TMath::Sin((i-1) * twopi / 39. + offset1);
      zpos = -(dbox2[2] + dela[2] + 6.);
      gMC->Gspos("ELL6", i, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
      zpos = dbox2[2] + dela[2] + 6.;
      gMC->Gspos("ELL6", i + 39, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
    }
    //#ifdef NEVER
  }
  
  // --- Define SSD with the 32+36 lay-out 
  
  if (fMinorVersionV3 >2 && fMinorVersionV3 < 6) {
    
    //--- Define ghost volume containing the Strip Detectors and fill it with air
    //     or vacuum 
    
    xxm    = (49.999-3.)/(70.-25.);
    dgh[0] = 0.;
    dgh[1] = 360.;
    dgh[2] = 4.;
    dgh[3] = -25. - (9.-3.01) / xxm - (9.01-9.) / xxm - 
      (27.-9.01) / xxm - (36.-27.) / xxm - (49.998-36.) / xxm;
    dgh[4] = 49.998;
    dgh[5] = 49.998;
    dgh[6] = -25. - (9.-3.01) / xxm - 
      (9.01-9.) / xxm - (27.-9.01) / xxm - (36.-27.) / xxm;
    dgh[7] = 36.;
    dgh[8] = 49.998;
    dgh[9] = (9.-3.01) / xxm + 25. + 
      (9.01-9.) / xxm + (27.-9.01) / xxm + (36.-27.) / xxm;
    dgh[10] = 36.;
    dgh[11] = 49.998;
    dgh[12] = (9.-3.01) / xxm + 25. + (9.01-9.) / xxm + 
      (27.-9.01) / xxm + (36.-27.) / xxm + (49.998-36.) / xxm;
    dgh[13] = 49.998;
    dgh[14] = 49.998;
    gMC->Gsvolu("IT56", "PCON", idtmed[275], dgh, 15);
    gMC->Gspos("IT56", 1, "ITSV", 0., 0., 0., 0, "ONLY");
    gMC->Gsatt("IT56", "SEEN", 0);
    
    // --- Layer #5 
    
    //        GOTO 6678           ! skip ITS layer no. 5 
    
    //--- Define a ghost volume containing a single ladder of layer #5 andfill
    //     it with air or vacuum 
    
    dbox1[0] = (0.0600+2.*0.0150)/2.;
    dbox1[1] = 3.75;
    dbox1[2] = 86.31/2.;
    gMC->Gsvolu("ISV1", "BOX ", idtmed[253], dbox1, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ISV1", "SEEN", 0);
    
    // --- Define a ghost volume containing the electronics and cooling of
    //     a single ladder of layer #5 and fill it with air or vacuum 
    
    dsrv[0] = (TMath::Sqrt(3.) / 2. * 4.2 + .47 + .05) / 2.;
    dsrv[1] = 3.75;
    dsrv[2] = 86.31/2.;
    gMC->Gsvolu("SSV1", "BOX ", idtmed[253], dsrv, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("SSV1", "SEEN", 0);
    
    // --- Define a ghost volume containing the end-ladder stuff of 
    //     a single ladder of layer #5 and fill it with air or vacuum 
    
    dela[0] = 2.;
    dela[1] = 3.5;
    dela[2] = 4.;
    gMC->Gsvolu("ELL5", "BOX ", idtmed[253], dela, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ELL5", "SEEN", 0);
    
    // --- Define a volume containing the sensitive part of the strips 
    //     (silicon, layer #5) 
    
    dits[0] = .015;
    dits[1] = 3.75;
    dits[2] = 2.1;
    gMC->Gsvolu("ITS5", "BOX ", idtmed[249], dits, 3);
    
    // --- Define a volume containing the electronics of the strips 
    //     (silicon, layer #5) 
    
    dchi[0] = .02;
    dchi[1] = 3.4;
    dchi[2] = .525;
    gMC->Gsvolu("SCH5", "BOX ", idtmed[250], dchi, 3);
    
    // --- Define the cooling tubes (aluminum, layer #5) 
    
    dtub[0] = .09;
    dtub[1] = dtub[0] + .01;
    dtub[2] = 86.31/2.;
    gMC->Gsvolu("STB5", "TUBE", idtmed[255], dtub, 3);
    
    // --- Define the cooling fluid (water or freon, layer #5) 
    
    dwat[0] = 0.;
    dwat[1] = .09;
    dwat[2] = 86.31/2.;
    gMC->Gsvolu("SWT5", "TUBE", idtmed[256], dwat, 3);
    //        CALL GSVOLU('SWT5','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
    
    //--- Define the (triangular) element of the heat bridge (carbon, layer #5)
    
    // water 
    dfra[0] = 120.;
    dfra[1] = 360.;
    dfra[2] = 3.;
    dfra[3] = 2.;
    dfra[4] = -.015;
    dfra[5] = TMath::Sqrt(3.) * 4.2 / 6.;
    dfra[6] = dfra[5] + .03;
    dfra[7] = .015;
    dfra[8] = dfra[5];
    dfra[9] = dfra[6];
    gMC->Gsvolu("SFR5", "PGON", idtmed[252], dfra, 10);
    
    // --- Define the element connecting the triangles of the heat bridge 
    //     (carbon, layer #5) 
    
    dcei[0] = 0.;
    dcei[1] = .03;
    dcei[2] = 86.31/2.;
    gMC->Gsvolu("SCE5", "TUBE", idtmed[252], dcei, 3);
    
    // --- Define the part of the end-ladder stuff made of plastic (G10FR4) 
    //     (layer #5) 
    
    dpla[0] = (10./(8.*7.))/2;
    dpla[1] = 3.5;
    dpla[2] = 4.;
    gMC->Gsvolu("EPL5", "BOX ", idtmed[262], dpla, 3);
    
    // --- Define the part of the end-ladder stuff made of copper (layer #5) 
    
    dcop[0] = (2./(8.*7.))/2;
    dcop[1] = 3.5;
    dcop[2] = 4.;
    gMC->Gsvolu("ECU5", "BOX ", idtmed[259], dcop, 3);
    
    // --- Define the part of the end-ladder stuff made of epoxy (layer #5) 
    
    depx[0] = (30./(8.*7.))/2.;
    depx[1] = 3.5;
    depx[2] = 4.;
    gMC->Gsvolu("EPX5", "BOX ", idtmed[262], depx, 3);
    
    // --- Define the part of the end-ladder stuff made of silicon (bus) 
    //     (layer #5) 
    
    dsil[0] = (20./(8.*7.))/2.;
    dsil[1] = 3.5;
    dsil[2] = 4.;
    gMC->Gsvolu("ESI5", "BOX ", idtmed[251], dsil, 3);
    
    // --- Place the end-ladder stuff into its mother (ELL5) 
    
    sep  = (4. - (dpla[0] + dcop[0] + depx[0] + dsil[0]) * 2.) / 3.;
    ypos = 0.;
    zpos = 0.;
    
    // --- Plastic 
    
    xpos = -dela[0] + dpla[0];
    gMC->Gspos("EPL5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Copper 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0];
    gMC->Gspos("ECU5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Epoxy 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0];
    gMC->Gspos("EPX5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Silicon (bus) 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0] * 2. + sep + dsil[0];
    gMC->Gspos("ESI5", 1, "ELL5", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the sensitive part of the strips into its mother (ISV1) 
    
    ypos = 0.;
    for (j = 1; j <= 22; ++j) {
      if (j % 2 == 0) xpos = dbox1[0] - dits[0];
      else            xpos = -dbox1[0] + dits[0];
      zpos = ((j - 1) - 10.) * 3.91 - 1.96;
      gMC->Gspos("ITS5", j, "ISV1", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the electronics of the strips into its mother (SSV1) 
    
    ypos = 0.;
    for (j = 1; j <= 22; ++j) {
      if (j % 2 == 0) xpos = -dsrv[0] + .28;
      else            xpos = -dsrv[0] + .28 - dits[0] * 2. - .03;
      zpos = ((j - 1) - 10.) * 3.91 - 1.96 + .85;
      gMC->Gspos("SCH5", j, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    }
    
    //--- Place the cooling tubes and the cooling fluid into their mother (SSV1)
    
    xpos = -dsrv[0] + .41;
    zpos = 0.;
    
    // --- Left tube (just a matter of convention) 
    
    ypos = -(2.25+0.1);
    gMC->Gspos("STB5", 1, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT5", 1, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right tube (just a matter of convention) 
    
    ypos = (2.25+0.1);
    gMC->Gspos("STB5", 2, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT5", 2, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the heat bridge elements into their mother (SSV1) 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 6. * 4.2;
    ypos = 0.;
    for (j = 1; j <= 23; ++j) {
      zpos = ((j - 1) - 10.) * 3.91 - 1.96 - 4.2/2.;
      gMC->Gspos("SFR5", j, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the elements connecting the triangles of the heat bridge 
    //     into their mother (SSV1) 
    
    zpos = 0.;
    
    // --- Left element (just a matter of convention) 
    
    xpos = -dsrv[0] + .47;
    ypos = -(2.1+0.015);
    gMC->Gspos("SCE5", 1, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right element 
    
    xpos = -dsrv[0] + .47;
    ypos = (2.1+0.015);
    gMC->Gspos("SCE5", 2, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Top element 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 2. * 4.2 + .015;
    ypos = 0.;
    gMC->Gspos("SCE5", 3, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the ghost volumes containing the strip ladders (ISV1), 
    //    electronics/cooling (SSV1) and end-ladder stuff (ELL5) of layer #5 in
    //     their mother volume (IT56) 
    
    offset1 = TMath::ATan2(.8, 36.6);
    offset2 = 5.2;
    rzero   = dbox1[0] + 36.6;
    runo    = dbox1[0] * 2. + 36.6 + dsrv[0];
    rtwo    = dbox1[0] * 2. + 36.6 + dela[0];
    for (i = 1; i <= 32; ++i) {
      atheta = (i-1) * twopi * raddeg / 32. + offset2;
      AliMatrix(idrotm[i+1499], 90., atheta, 90., atheta + 90., 0., 0.);
      
      // --- Strip ladders 
      
      xpos = rzero * TMath::Cos((i-1) * twopi / 32. + offset1);
      ypos = rzero * TMath::Sin((i-1) * twopi / 32. + offset1);
      zpos = 0.;
      gMC->Gspos("ISV1", i, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
      
      // --- Electronics/cooling 
      
      xpos = runo * TMath::Cos((i-1) * twopi / 32. + offset1);
      ypos = runo * TMath::Sin((i-1) * twopi / 32. + offset1);
      zpos = 0.;
      gMC->Gspos("SSV1", i, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
      
      // --- End-ladders (nagative-Z and positive-Z) 
      
      xpos = rtwo * TMath::Cos((i-1) * twopi / 32. + offset1);
      ypos = rtwo * TMath::Sin((i-1) * twopi / 32. + offset1);
      zpos = -(dbox1[2] + dela[2] + 6.);
      gMC->Gspos("ELL5", i, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
      zpos = dbox1[2] + dela[2] + 6.;
      gMC->Gspos("ELL5", i + 35, "IT56", xpos, ypos, zpos, idrotm[i+1499], "ONLY");
    }
    
    
    // --- Layer #6 
    
    //        GOTO 6778           ! skip ITS layer no. 6 
    
    //--- Define a ghost volume containing a single ladder of layer #6 andfill
    //     it with air or vacuum 
    
    dbox2[0] = (0.0600+2.*0.0150)/2.;
    dbox2[1] = 3.75;
    dbox2[2] = 94.13/2.;
    gMC->Gsvolu("ISV2", "BOX ", idtmed[253], dbox2, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ISV2", "SEEN", 0);
    
    // --- Define a ghost volume containing the electronics and cooling of
    //     a single ladder of layer #6 and fill it with air or vacuum 
    
    dsrv[0] = (TMath::Sqrt(3.) / 2. * 4.2 + .47 + .05) / 2.;
    dsrv[1] = 3.75;
    dsrv[2] = 94.13/2.;
    gMC->Gsvolu("SSV2", "BOX ", idtmed[253], dsrv, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("SSV2", "SEEN", 0);
    
    // --- Define a ghost volume containing the end-ladder stuff of 
    //     a single ladder of layer #6 and fill it with air or vacuum 
    
    dela[0] = 2.;
    dela[1] = 3.5;
    dela[2] = 4.;
    gMC->Gsvolu("ELL6", "BOX ", idtmed[253], dela, 3);
    
    // --- Make the ghost volume invisible 
    
    gMC->Gsatt("ELL6", "SEEN", 0);
    
    // --- Define a volume containing the sensitive part of the strips 
    //     (silicon, layer #6) 
    
    dits[0] = .015;
    dits[1] = 3.75;
    dits[2] = 2.1;
    gMC->Gsvolu("ITS6", "BOX ", idtmed[249], dits, 3);
    
    // --- Define a volume containing the electronics of the strips 
    //     (silicon, layer #6) 
    
    dchi[0] = .02;
    dchi[1] = 3.4;
    dchi[2] = .525;
    gMC->Gsvolu("SCH6", "BOX ", idtmed[250], dchi, 3);
    
    // --- Define the cooling tubes (aluminum, layer #6) 
    
    dtub[0] = .09;
    dtub[1] = dtub[0] + .01;
    dtub[2] = 94.13/2.;
    gMC->Gsvolu("STB6", "TUBE", idtmed[255], dtub, 3);
    
    // --- Define the cooling fluid (water or freon, layer #6) 
    
    dwat[0] = 0.;
    dwat[1] = .09;
    dwat[2] = 94.13/2.;
    gMC->Gsvolu("SWT6", "TUBE", idtmed[256], dwat, 3);
    //        CALL GSVOLU('SWT6','TUBE',IDTMED(258),DWAT,3,IOUT)   ! freon
    
    //--- Define the (triangular) element of the heat bridge (carbon, layer #6)
    
    // water 
    dfra[0] = 120.;
    dfra[1] = 360.;
    dfra[2] = 3.;
    dfra[3] = 2.;
    dfra[4] = -.015;
    dfra[5] = TMath::Sqrt(3.) * 4.2 / 6.;
    dfra[6] = dfra[5] + .03;
    dfra[7] = .015;
    dfra[8] = dfra[5];
    dfra[9] = dfra[6];
    gMC->Gsvolu("SFR6", "PGON", idtmed[252], dfra, 10);
    
    // --- Define the element connecting the triangles of the heat bridge 
    //     (carbon, layer #6) 
    
    dcei[0] = 0.;
    dcei[1] = .03;
    dcei[2] = 94.13/2.;
    gMC->Gsvolu("SCE6", "TUBE", idtmed[252], dcei, 3);
    
    // --- Define the part of the end-ladder stuff made of plastic (G10FR4) 
    //     (layer #6) 
    
    dpla[0] = (10./(8.*7.))/2;
    dpla[1] = 3.5;
    dpla[2] = 4.;
    gMC->Gsvolu("EPL6", "BOX ", idtmed[262], dpla, 3);
    
    // --- Define the part of the end-ladder stuff made of copper (layer #6) 
    
    dcop[0] = (2./(8.*7.))/2;
    dcop[1] = 3.5;
    dcop[2] = 4.;
    gMC->Gsvolu("ECU6", "BOX ", idtmed[259], dcop, 3);
    
    // --- Define the part of the end-ladder stuff made of epoxy (layer #6) 
    
    depx[0] = (30./(8.*7.))/2.;
    depx[1] = 3.5;
    depx[2] = 4.;
    gMC->Gsvolu("EPX6", "BOX ", idtmed[262], depx, 3);
    
    // --- Define the part of the end-ladder stuff made of silicon (bus) 
    //     (layer #6) 
    
    dsil[0] = (20./(8.*7.))/2.;
    dsil[1] = 3.5;
    dsil[2] = 4.;
    gMC->Gsvolu("ESI6", "BOX ", idtmed[251], dsil, 3);
    
    // --- Place the end-ladder stuff into its mother (ELL5) 
    
    sep  = (4. - (dpla[0] + dcop[0] + depx[0] + dsil[0]) * 2.) / 3.;
    ypos = 0.;
    zpos = 0.;
    
    // --- Plastic 
    
    xpos = -dela[0] + dpla[0];
    gMC->Gspos("EPL6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Copper 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0];
    gMC->Gspos("ECU6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Epoxy 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0];
    gMC->Gspos("EPX6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Silicon (bus) 
    
    xpos = -dela[0] + dpla[0] * 2. + sep + dcop[0] * 2. + sep + depx[0] * 2. + sep + dsil[0];
    gMC->Gspos("ESI6", 1, "ELL6", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the sensitive part of the strips into its mother (ISV2) 
    
    ypos = 0.;
    for (j = 1; j <= 24; ++j) {
      if (j % 2 == 0) xpos = -dbox2[0] + dits[0];
      else            xpos = dbox2[0] - dits[0];
      zpos = ((j - 1) - 11.) * 3.91 - 1.96;
      gMC->Gspos("ITS6", j, "ISV2", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the electronics of the strips into its mother (SSV2) 
    
    ypos = 0.;
    for (j = 1; j <= 24; ++j) {
      if (j % 2 == 0) xpos = -dsrv[0] + .28 - dits[0] * 2. - .03;
      else            xpos = -dsrv[0] + .28;
      zpos = ((j - 1) - 11.) * 3.91 - 1.96 + .85;
      gMC->Gspos("SCH5", j, "SSV1", xpos, ypos, zpos, 0, "ONLY");
    }
    
    //--- Place the cooling tubes and the cooling fluid into their mother (SSV2)
    
    xpos = -dsrv[0] + .41;
    zpos = 0.;
    
    // --- Left tube (just a matter of convention) 
    
    ypos = -(2.25+0.1);
    gMC->Gspos("STB6", 1, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT6", 1, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right tube (just a matter of convention) 
    
    ypos = (2.25+0.1);
    gMC->Gspos("STB6", 2, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    gMC->Gspos("SWT6", 2, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the heat bridge elements into their mother (SSV2) 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 6. * 4.2;
    ypos = 0.;
    for (j = 1; j <= 25; ++j) {
      zpos = ((j - 1) - 11.) * 3.91 - 1.96 - 4.2/2.;
      gMC->Gspos("SFR6", j, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    }
    
    // --- Place the elements connecting the triangles of the heat bridge 
    //     into their mother (SSV2) 
    
    zpos = 0.;
    
    // --- Left element (just a matter of convention) 
    
    xpos = -dsrv[0] + .47;
    ypos = -(2.1+0.015);
    gMC->Gspos("SCE6", 1, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Right element 
    
    xpos = -dsrv[0] + .47;
    ypos = (2.1+0.015);
    gMC->Gspos("SCE6", 2, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Top element 
    
    xpos = -dsrv[0] + .47 + TMath::Sqrt(3.) / 2. * 4.2 + .015;
    ypos = 0.;
    gMC->Gspos("SCE6", 3, "SSV2", xpos, ypos, zpos, 0, "ONLY");
    
    // --- Place the ghost volumes containing the strip ladders (ISV2), 
    //    electronics/cooling (SSV2) and end-ladder stuff (ELL6) of layer #6 in
    //     their mother volume (IT56) 
    
    offset1 = TMath::ATan2(.9, 41.2);
    offset2 = 5.2;
    rzero   = dbox2[0] + 41.2;
    runo    = dbox2[0] * 2. + 41.2 + dsrv[0];
    rtwo    = dbox2[0] * 2. + 41.2 + dela[0];
    for (i = 1; i <= 36; ++i) {
      atheta = (i-1) * twopi * raddeg / 36. + offset2;
      AliMatrix(idrotm[i+1599], 90., atheta, 90., atheta + 90., 0., 0.);
      
      // --- Strip ladders 
      
      xpos = rzero * TMath::Cos((i-1) * twopi / 36. + offset1);
      ypos = rzero * TMath::Sin((i-1) * twopi / 36. + offset1);
      zpos = 0.;
      gMC->Gspos("ISV2", i, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
      
      // --- Electronics/cooling 
      
      xpos = runo * TMath::Cos((i-1) * twopi / 36. + offset1);
      ypos = runo * TMath::Sin((i-1) * twopi / 36. + offset1);
      zpos = 0.;
      gMC->Gspos("SSV2", i, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
      
      // --- End-ladders (nagative-Z and positive-Z) 
      
      xpos = rtwo * TMath::Cos((i-1) * twopi / 36. + offset1);
      ypos = rtwo * TMath::Sin((i-1) * twopi / 36. + offset1);
      zpos = -(dbox2[2] + dela[2] + 6.);
      gMC->Gspos("ELL6", i, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
      zpos = dbox2[2] + dela[2] + 6.;
      gMC->Gspos("ELL6", i + 39, "IT56", xpos, ypos, zpos, idrotm[i+1599], "ONLY");
    }
    
    
  }
  
  //************************************************************************
  //*                                                                      *
  //*             E N D - C A P S     A N D     F R A M E S                *
  //*             =========================================                *
  //*                                                                      *
  //************************************************************************
  
  // --- Define a dummy cylinder for multiple scattering tests 
  
  //     GOTO 7890             ! skip dummy cylinder for multiple scatteringtests
  
  //      DITS(1)=49. 
  //      DITS(2)=DITS(1)+0.1 
  //      DITS(3)=60.3 
  //      CALL GSVOLU('ITST','TUBE',IDTMED(255),DITS,3,IOUT) 
  //      CALL GSPOS('ITST',1,'ITSV',0.,0.,0.,0,'ONLY') 
  // 7890  CONTINUE 
  
  // --- The 0.74% X0 outer wall (C) of the gas vessel at r=50cm --- 
  
  //      GOTO 8901                                    ! skip outer wall 
  
  if (fMinorVersionV3 == 0 || fMinorVersionV3 == 3) {
    
    dits[0] = 49.9;
    dits[1] = dits[0] + .06926;
    dits[2] = dpcb[2] * 2. + 62.7 - 10.5;
    // old value 60.3 
    gMC->Gsvolu("ITSG", "TUBE", idtmed[274], dits, 3);
    gMC->Gspos("ITSG", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  } else {
    goto L8901;
  }
 L8901:
  // --- The frame between the end-caps (octagonal lay-out) --- 
  
  //     GOTO 9012                                    ! skip octagonal frame
  
  if (fMinorVersionV3 == 1) {
    
    rzero    = 34.;
    dtra[0]  = .92;
    dtra[1]  = 1.;
    dtra[2]  = dpcb[2] * 2. + 50.5 - 10.5;
    dtra1[0] = .92;
    dtra1[1] = 1.;
    dtra1[2] = TMath::Sqrt(dtra[2] * dtra[2] + (55.4*55.4-50.5*50.5))/2.;
    angle    = 45.;
    offset   = angle / 2.;
    for (i = 0; i < 8; ++i) {
      xtra[i] = rzero * TMath::Cos(i * angle * degrad);
      ytra[i] = rzero * TMath::Sin(i * angle * degrad);
      ztra[i] = 0.;
      gMC->Gsvolu(natra[i], "TUBE", idtmed[274], dtra, 3);
      gMC->Gspos(natra[i], 1, "ITSV", xtra[i], ytra[i], ztra[i], 0, "ONLY");
    }
    
    atheta = 22.5;
    aphi1  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2  = 180. - aphi1;
    xpos   = (xtra[0] + xtra[1]) / 2.;
    ypos   = (ytra[0] + ytra[1]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[0], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5100], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[0], 1, "ITSV", xpos, ypos, zpos, idrotm[5100], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[1], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5101], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[1], 1, "ITSV", xpos, ypos, zpos, idrotm[5101], "ONLY");
    
    atheta = 67.5;
    aphi2  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1  = 180. - aphi2;
    xpos   = (xtra[1] + xtra[2]) / 2.;
    ypos   = (ytra[1] + ytra[2]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[2], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5102], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[2], 1, "ITSV", xpos, ypos, zpos, idrotm[5102], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[3], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5103], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[3], 1, "ITSV", xpos, ypos, zpos, idrotm[5103], "ONLY");
    
    atheta = 112.5;
    aphi1  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2  = 180. - aphi1;
    xpos   = (xtra[2] + xtra[3]) / 2.;
    ypos   = (ytra[2] + ytra[3]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[4], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5104], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[4], 1, "ITSV", xpos, ypos, zpos, idrotm[5104], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[5], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5105], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[5], 1, "ITSV", xpos, ypos, zpos, idrotm[5105], "ONLY");
    
    atheta = 157.5;
    aphi2  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1  = 180. - aphi2;
    xpos   = (xtra[3] + xtra[4]) / 2.;
    ypos   = (ytra[3] + ytra[4]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[6], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5106], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[6], 1, "ITSV", xpos, ypos, zpos, idrotm[5106], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[7], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5107], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[7], 1, "ITSV", xpos, ypos, zpos, idrotm[5107], "ONLY");
    
    atheta = 22.5;
    aphi2  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1  = 180. - aphi2;
    xpos   = (xtra[4] + xtra[5]) / 2.;
    ypos   = (ytra[4] + ytra[5]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[8], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5108], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[8], 1, "ITSV", xpos, ypos, zpos, idrotm[5108], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[9], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5109], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[9], 1, "ITSV", xpos, ypos, zpos, idrotm[5109], "ONLY");
    
    atheta = 67.5;
    aphi1  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2  = 180. - aphi1;
    xpos   = (xtra[5] + xtra[6]) / 2.;
    ypos   = (ytra[5] + ytra[6]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[10], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5110], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[10], 1, "ITSV", xpos, ypos, zpos, idrotm[5110], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[11], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5111], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[11], 1, "ITSV", xpos, ypos, zpos, idrotm[5111], "ONLY");
    
    atheta = 112.5;
    aphi2  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1  = 180. - aphi2;
    xpos   = (xtra[6] + xtra[7]) / 2.;
    ypos   = (ytra[6] + ytra[7]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[12], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5112], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[12], 1, "ITSV", xpos, ypos, zpos, idrotm[5112], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[13], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5113], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[13], 1, "ITSV", xpos, ypos, zpos, idrotm[5113], "ONLY");
    
    atheta = 157.5;
    aphi1  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2  = 180. - aphi1;
    xpos   = (xtra[7] + xtra[0]) / 2.;
    ypos   = (ytra[7] + ytra[0]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[14], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5114], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[14], 1, "ITSV", xpos, ypos, zpos, idrotm[5114], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[15], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5115], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[15], 1, "ITSV", xpos, ypos, zpos, idrotm[5115], "ONLY");
    
    
  } else if (fMinorVersionV3 == 4) {
    
    
    rzero    = 34.;
    dtra[0]  = .92;
    dtra[1]  = 1.;
    dtra[2]  = dpcb[2] * 2. + 50.5 - 10.5;
    dtra1[0] = .92;
    dtra1[1] = 1.;
    dtra1[2] = TMath::Sqrt(dtra[2] * dtra[2] + (55.4*55.4-50.5*50.5))/2.;
    angle    = 45.;
    offset   = angle / 2.;
    for (i = 0; i < 8; ++i) {
      xtra[i] = rzero * TMath::Cos(i * angle * degrad);
      ytra[i] = rzero * TMath::Sin(i * angle * degrad);
      ztra[i] = 0.;
      gMC->Gsvolu(natra[i], "TUBE", idtmed[274], dtra, 3);
      gMC->Gspos(natra[i], 1, "ITSV", xtra[i], ytra[i], ztra[i], 0, "ONLY");
    }
    
    atheta = 22.5;
    aphi1  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2  = 180. - aphi1;
    xpos   = (xtra[0] + xtra[1]) / 2.;
    ypos   = (ytra[0] + ytra[1]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[0], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5100], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[0], 1, "ITSV", xpos, ypos, zpos, idrotm[5100], "ONLY");
    zpos   = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[1], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5101], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[1], 1, "ITSV", xpos, ypos, zpos, idrotm[5101], "ONLY");
    
    atheta = 67.5;
    aphi2  = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1  = 180. - aphi2;
    xpos   = (xtra[1] + xtra[2]) / 2.;
    ypos   = (ytra[1] + ytra[2]) / 2.;
    zpos   = dtra[2] / 2.;
    gMC->Gsvolu(natra1[2], "TUBE", idtmed[274], dtra1, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5102], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[2], 1, "ITSV", xpos, ypos, zpos, idrotm[5102], "ONLY");
    zpos = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[3], "TUBE", idtmed[274], dtra1, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5103], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[3], 1, "ITSV", xpos, ypos, zpos, idrotm[5103], "ONLY");
    
    atheta = 112.5;
    aphi1 = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos  = (xtra[2] + xtra[3]) / 2.;
    ypos  = (ytra[2] + ytra[3]) / 2.;
    zpos  = dtra[2] / 2.;
    gMC->Gsvolu(natra1[4], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5104], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[4], 1, "ITSV", xpos, ypos, zpos, idrotm[5104], "ONLY");
    zpos  = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[5], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5105], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[5], 1, "ITSV", xpos, ypos, zpos, idrotm[5105], "ONLY");
    
    atheta = 157.5;
    aphi2 = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos  = (xtra[3] + xtra[4]) / 2.;
    ypos  = (ytra[3] + ytra[4]) / 2.;
    zpos  = dtra[2] / 2.;
    gMC->Gsvolu(natra1[6], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5106], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[6], 1, "ITSV", xpos, ypos, zpos, idrotm[5106], "ONLY");
    zpos  = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[7], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5107], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[7], 1, "ITSV", xpos, ypos, zpos, idrotm[5107], "ONLY");
    
    atheta = 22.5;
    aphi2 = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos  = (xtra[4] + xtra[5]) / 2.;
    ypos  = (ytra[4] + ytra[5]) / 2.;
    zpos  = dtra[2] / 2.;
    gMC->Gsvolu(natra1[8], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5108], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[8], 1, "ITSV", xpos, ypos, zpos, idrotm[5108], "ONLY");
    zpos  = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[9], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5109], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[9], 1, "ITSV", xpos, ypos, zpos, idrotm[5109], "ONLY");
    
    atheta = 67.5;
    aphi1 = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos = (xtra[5] + xtra[6]) / 2.;
    ypos = (ytra[5] + ytra[6]) / 2.;
    zpos = dtra[2] / 2.;
    gMC->Gsvolu(natra1[10], "TUBE", idtmed[274], dtra1, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5110], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[10], 1, "ITSV", xpos, ypos, zpos, idrotm[5110], "ONLY");
    zpos = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[11], "TUBE", idtmed[274], dtra1, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5111], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[11], 1, "ITSV", xpos, ypos, zpos, idrotm[5111], "ONLY");
    
    atheta = 112.5;
    aphi2 = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos  = (xtra[6] + xtra[7]) / 2.;
    ypos  = (ytra[6] + ytra[7]) / 2.;
    zpos  = dtra[2] / 2.;
    gMC->Gsvolu(natra1[12], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5112], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[12], 1, "ITSV", xpos, ypos, zpos, idrotm[5112], "ONLY");
    zpos  = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[13], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5113], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[13], 1, "ITSV", xpos, ypos, zpos, idrotm[5113], "ONLY");
    
    atheta = 157.5;
    aphi1 = TMath::ACos(dtra[2] / TMath::Sqrt(dtra[2] * dtra[2] + (50.5 / cos(28.*degrad) * (50.5 / cos(28.*degrad))- 50.5*50.5))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos  = (xtra[7] + xtra[0]) / 2.;
    ypos  = (ytra[7] + ytra[0]) / 2.;
    zpos  = dtra[2] / 2.;
    gMC->Gsvolu(natra1[14], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5114], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra1[14], 1, "ITSV", xpos, ypos, zpos, idrotm[5114], "ONLY");
    zpos  = -dtra[2] / 2.;
    gMC->Gsvolu(natra1[15], "TUBE", idtmed[274], dtra1, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5115], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra1[15], 1, "ITSV", xpos, ypos, zpos, idrotm[5115], "ONLY");
  } else {
    goto L9012;
  }
  
 L9012:
  // --- The frame between the end-caps (hexagonal lay-out) --- 
  
  //     GOTO 9123                                    ! skip hexagonal frame
  
  if (fMinorVersionV3 == 2) {
    
    rzero    = 33.5;
    dtra2[0] = .92;
    dtra2[1] = 1.;
    dtra2[2] = dpcb[2] * 2. + 50. - 10.5;
    dtra3[0] = .92;
    dtra3[1] = 1.;
    dtra3[2] = 16.75;
    dtra4[0] = .92;
    dtra4[1] = 1.;
    dtra4[2] = TMath::Sqrt(dtra2[2] * dtra2[2] + (59.9*59.9-50.*50.)) / 2.;
    angle = 60.;
    offset = angle / 2.;
    for (i = 0; i < 6; ++i) {
      xtra1[i] = rzero * TMath::Cos((i * angle + offset) *degrad);
      ytra1[i] = rzero * TMath::Sin((i * angle + offset) *degrad);
      ztra1[i] = 0.;
      gMC->Gsvolu(natra2[i], "TUBE", idtmed[274], dtra2, 3);
      gMC->Gspos(natra2[i], 1, "ITSV", xtra1[i], ytra1[i], ztra1[i], 0, "ONLY");
    }
    
    atheta = 60.;
    aphi = 90.;
    xpos = (xtra1[0] + xtra1[1]) / 2.;
    ypos = (ytra1[0] + ytra1[1]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[0], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5200], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[0], 1, "ITSV", xpos, ypos, zpos, idrotm[5200], "ONLY");
    
    atheta = 120.;
    aphi = 90.;
    xpos = (xtra1[1] + xtra1[2]) / 2.;
    ypos = (ytra1[1] + ytra1[2]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[1], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5201], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[1], 1, "ITSV", xpos, ypos, zpos, idrotm[5201], "ONLY");
    
    atheta = 180.;
    aphi = 90.;
    xpos = (xtra1[2] + xtra1[3]) / 2.;
    ypos = (ytra1[2] + ytra1[3]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[2], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5202], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[2], 1, "ITSV", xpos, ypos, zpos, idrotm[5202], "ONLY");
    
    atheta = 60.;
    aphi = 90.;
    xpos = (xtra1[3] + xtra1[4]) / 2.;
    ypos = (ytra1[3] + ytra1[4]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[3], "TUBE", idtmed[274], dtra3, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5203], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[3], 1, "ITSV", xpos, ypos, zpos, idrotm[5203], "ONLY");
    
    atheta = 120.;
    aphi = 90.;
    xpos = (xtra1[4] + xtra1[5]) / 2.;
    ypos = (ytra1[4] + ytra1[5]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[4], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5204], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[4], 1, "ITSV", xpos, ypos, zpos, idrotm[5204], "ONLY");
    
    atheta = 180.;
    aphi = 90.;
    xpos = (xtra1[5] + xtra1[0]) / 2.;
    ypos = (ytra1[5] + ytra1[0]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[5], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5205], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[5], 1, "ITSV", xpos, ypos, zpos, idrotm[5205], "ONLY");
    
    atheta = 60.;
    aphi2 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos  = (xtra1[0] + xtra1[1]) / 2.;
    ypos  = (ytra1[0] + ytra1[1]) / 2.;
    zpos  = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[0], "TUBE", idtmed[274], dtra4, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5210], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[0], 1, "ITSV", xpos, ypos, zpos, idrotm[5210], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[1], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5211], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[1], 1, "ITSV", xpos, ypos, zpos, idrotm[5211], "ONLY");
    
    atheta = 120.;
    aphi1 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos = (xtra1[1] + xtra1[2]) / 2.;
    ypos = (ytra1[1] + ytra1[2]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[2], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5212], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[2], 1, "ITSV", xpos, ypos, zpos, idrotm[5212], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[3], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5213], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[3], 1, "ITSV", xpos, ypos, zpos, idrotm[5213], "ONLY");
    
    atheta = 180.;
    aphi2 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos  = (xtra1[2] + xtra1[3]) / 2.;
    ypos  = (ytra1[2] + ytra1[3]) / 2.;
    zpos  = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[4], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5214], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[4], 1, "ITSV", xpos, ypos, zpos, idrotm[5214], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[5], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5215], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[5], 1, "ITSV", xpos, ypos, zpos, idrotm[5215], "ONLY");
    atheta = 180.;
    aphi1 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad)) 
								      - 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos = (xtra1[2] + xtra1[3]) / 2.;
    ypos = (ytra1[2] + ytra1[3]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[6], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5216], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[6], 1, "ITSV", xpos, ypos, zpos, idrotm[5216], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[7], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5217], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[7], 1, "ITSV", xpos, ypos, zpos, idrotm[5217], "ONLY");
    
    atheta = 60.;
    aphi2 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos = (xtra1[3] + xtra1[4]) / 2.;
    ypos = (ytra1[3] + ytra1[4]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[8], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5218], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[8], 1, "ITSV", xpos, ypos, zpos, idrotm[5218], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[9], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5219], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[9], 1, "ITSV", xpos, ypos, zpos, idrotm[5219], "ONLY");
    
    atheta = 120.;
    aphi1 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos  = (xtra1[4] + xtra1[5]) / 2.;
    ypos  = (ytra1[4] + ytra1[5]) / 2.;
    zpos  = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[10], "TUBE", idtmed[274], dtra4, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5220], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[10], 1, "ITSV", xpos, ypos, zpos, idrotm[5220], "ONLY");
    zpos  = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[11], "TUBE", idtmed[274], dtra4, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5221], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[11], 1, "ITSV", xpos, ypos, zpos, idrotm[5221], "ONLY");
    
    atheta = 180.;
    aphi2  = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1  = 180. - aphi2;
    xpos   = (xtra1[5] + xtra1[0]) / 2.;
    ypos   = (ytra1[5] + ytra1[0]) / 2.;
    zpos   = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[12], "TUBE", idtmed[274], dtra4, 3);
    r2     = atheta + 90.;
    r3     = atheta + 90.;
    AliMatrix(idrotm[5222], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[12], 1, "ITSV", xpos, ypos, zpos, idrotm[5222], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[13], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5223], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[13], 1, "ITSV", xpos, ypos, zpos, idrotm[5223], "ONLY");
    atheta = 180.;
    aphi1 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos  = (xtra1[5] + xtra1[0]) / 2.;
    ypos  = (ytra1[5] + ytra1[0]) / 2.;
    zpos  = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[14], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5224], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[14], 1, "ITSV", xpos, ypos, zpos, idrotm[5224], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[15], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5225], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[15], 1, "ITSV", xpos, ypos, zpos, idrotm[5225], "ONLY");
    
    
  } else if (fMinorVersionV3 == 5) {
    
    
    rzero    = 33.5;
    dtra2[0] = .92;
    dtra2[1] = 1.;
    dtra2[2] = dpcb[2] * 2. + 50. - 10.5;
    dtra3[0] = .92;
    dtra3[1] = 1.;
    dtra3[2] = 16.75;
    dtra4[0] = .92;
    dtra4[1] = 1.;
    dtra4[2] = TMath::Sqrt(dtra2[2] * dtra2[2] + (59.9*59.9-50.*50.)) / 2.;
    angle  = 60.;
    offset = angle / 2.;
    for (i = 0; i < 6; ++i) {
      xtra1[i] = rzero * TMath::Cos((i * angle + offset) *degrad);
      ytra1[i] = rzero * TMath::Sin((i * angle + offset) *degrad);
      ztra1[i] = 0.;
      gMC->Gsvolu(natra2[i], "TUBE", idtmed[274], dtra2, 3);
      gMC->Gspos(natra2[i], 1, "ITSV", xtra1[i], ytra1[i], ztra1[i], 0, "ONLY");
    }
    
    atheta = 60.;
    aphi = 90.;
    xpos = (xtra1[0] + xtra1[1]) / 2.;
    ypos = (ytra1[0] + ytra1[1]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[0], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5200], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[0], 1, "ITSV", xpos, ypos, zpos, idrotm[5200], "ONLY");
    
    atheta = 120.;
    aphi = 90.;
    xpos = (xtra1[1] + xtra1[2]) / 2.;
    ypos = (ytra1[1] + ytra1[2]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[1], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5201], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[1], 1, "ITSV", xpos, ypos, zpos, idrotm[5201], "ONLY");
    
    atheta = 180.;
    aphi = 90.;
    xpos = (xtra1[2] + xtra1[3]) / 2.;
    ypos = (ytra1[2] + ytra1[3]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[2], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5202], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[2], 1, "ITSV", xpos, ypos, zpos, idrotm[5202], "ONLY");
    
    atheta = 60.;
    aphi   = 90.;
    xpos   = (xtra1[3] + xtra1[4]) / 2.;
    ypos   = (ytra1[3] + ytra1[4]) / 2.;
    zpos   = 0.;
    gMC->Gsvolu(natra3[3], "TUBE", idtmed[274], dtra3, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5203], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[3], 1, "ITSV", xpos, ypos, zpos, idrotm[5203], "ONLY");
    
    atheta = 120.;
    aphi = 90.;
    xpos = (xtra1[4] + xtra1[5]) / 2.;
    ypos = (ytra1[4] + ytra1[5]) / 2.;
    zpos = 0.;
    gMC->Gsvolu(natra3[4], "TUBE", idtmed[274], dtra3, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5204], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[4], 1, "ITSV", xpos, ypos, zpos, idrotm[5204], "ONLY");
    
    atheta = 180.;
    aphi   = 90.;
    xpos   = (xtra1[5] + xtra1[0]) / 2.;
    ypos   = (ytra1[5] + ytra1[0]) / 2.;
    zpos   = 0.;
    gMC->Gsvolu(natra3[5], "TUBE", idtmed[274], dtra3, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5205], 90., atheta, aphi + 90., r2, aphi, r3);
    gMC->Gspos(natra3[5], 1, "ITSV", xpos, ypos, zpos, idrotm[5205], "ONLY");
    
    atheta = 60.;
    aphi2 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos = (xtra1[0] + xtra1[1]) / 2.;
    ypos = (ytra1[0] + ytra1[1]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[0], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5210], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[0], 1, "ITSV", xpos, ypos, zpos, idrotm[5210], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[1], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5211], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[1], 1, "ITSV", xpos, ypos, zpos, idrotm[5211], "ONLY");
    
    atheta = 120.;
    aphi1 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos = (xtra1[1] + xtra1[2]) / 2.;
    ypos = (ytra1[1] + ytra1[2]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[2], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5212], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[2], 1, "ITSV", xpos, ypos, zpos, idrotm[5212], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[3], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5213], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[3], 1, "ITSV", xpos, ypos, zpos, idrotm[5213], "ONLY");
    
    atheta = 180.;
    aphi2 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos = (xtra1[2] + xtra1[3]) / 2.;
    ypos = (ytra1[2] + ytra1[3]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[4], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5214], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[4], 1, "ITSV", xpos, ypos, zpos, idrotm[5214], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[5], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5215], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[5], 1, "ITSV", xpos, ypos, zpos, idrotm[5215], "ONLY");
    atheta = 180.;
    aphi1 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos = (xtra1[2] + xtra1[3]) / 2.;
    ypos = (ytra1[2] + ytra1[3]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[6], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5216], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[6], 1, "ITSV", xpos, ypos, zpos, idrotm[5216], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[7], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5217], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[7], 1, "ITSV", xpos, ypos, zpos, idrotm[5217], "ONLY");
    
    atheta = 60.;
    aphi2 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos = (xtra1[3] + xtra1[4]) / 2.;
    ypos = (ytra1[3] + ytra1[4]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[8], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5218], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[8], 1, "ITSV", xpos, ypos, zpos, idrotm[5218], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[9], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5219], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[9], 1, "ITSV", xpos, ypos, zpos, idrotm[5219], "ONLY");
    
    atheta = 120.;
    aphi1 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos  = (xtra1[4] + xtra1[5]) / 2.;
    ypos  = (ytra1[4] + ytra1[5]) / 2.;
    zpos  = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[10], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5220], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[10], 1, "ITSV", xpos, ypos, zpos, idrotm[5220], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[11], "TUBE", idtmed[274], dtra4, 3);
    r2 = atheta + 90.;
    r3 = atheta + 90.;
    AliMatrix(idrotm[5221], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[11], 1, "ITSV", xpos, ypos, zpos, idrotm[5221], "ONLY");
    
    atheta = 180.;
    aphi2 = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi1 = 180. - aphi2;
    xpos  = (xtra1[5] + xtra1[0]) / 2.;
    ypos  = (ytra1[5] + ytra1[0]) / 2.;
    zpos  = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[12], "TUBE", idtmed[274], dtra4, 3);
    r2    = atheta + 90.;
    r3    = atheta + 90.;
    AliMatrix(idrotm[5222], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[12], 1, "ITSV", xpos, ypos, zpos, idrotm[5222], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[13], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5223], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[13], 1, "ITSV", xpos, ypos, zpos, idrotm[5223], "ONLY");
    atheta = 180.;
    aphi1  = TMath::ACos(dtra2[2] / TMath::Sqrt(dtra2[2] * dtra2[2] + (50. / cos(34.*degrad) * (50. / cos(34.*degrad))- 50.*50.))) * raddeg;
    aphi2 = 180. - aphi1;
    xpos = (xtra1[5] + xtra1[0]) / 2.;
    ypos = (ytra1[5] + ytra1[0]) / 2.;
    zpos = dtra2[2] / 2.;
    gMC->Gsvolu(natra4[14], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5224], 90., atheta, aphi1 + 90., r2, aphi1, r3);
    gMC->Gspos(natra4[14], 1, "ITSV", xpos, ypos, zpos, idrotm[5224], "ONLY");
    zpos = -dtra2[2] / 2.;
    gMC->Gsvolu(natra4[15], "TUBE", idtmed[274], dtra4, 3);
    r2   = atheta + 90.;
    r3   = atheta + 90.;
    AliMatrix(idrotm[5225], 90., atheta, aphi2 + 90., r2, aphi2, r3);
    gMC->Gspos(natra4[15], 1, "ITSV", xpos, ypos, zpos, idrotm[5225], "ONLY");
  } else {
    goto L9123;
  }
  
 L9123:
  // --- Define the end-caps 
  
  //      GOTO 9234             ! skip both end-caps 
  
  // --- Define the Z>0 end-cap 
  
  //      GOTO 9345             ! skip the Z>0 end-cap 
  
  dcone[0] = 16.75;
  dcone[1] = 12.;
  dcone[2] = 12.02;
  dcone[3] = (338.-3.)*455./(338.-3.-10.)/10.;
  dcone[4] = .02 / TMath::Cos(45.*degrad) + (338.-3.)*455./(338.-3.-10.)/10.;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = dpcb[2] * 2. + (583.+(338.-3.))/2./10. - 10.5;
  // end-ladder electro 
  gMC->Gsvolu("RCON", "CONE", idtmed[274], dcone, 5);
  gMC->Gspos("RCON", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dtube[0] = .02 / TMath::Cos(45.*degrad) + (338.-3.)*455./(338.-3.-10.)/10.;
  dtube[1] = 49.9;
  // In the Simonetti's drawings 52. In the TP 50. 
  dtube[2] = .15;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = dpcb[2] * 2. + (583./2.+(338-1.5))/10. - 10.5;
  // end-ladder electro 
  gMC->Gsvolu("RTB1", "TUBE", idtmed[274], dtube, 3);
  gMC->Gspos("RTB1", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dtube[0] = 10.5;
  dtube[1] = 12.;
  dtube[2] = 26.8/2./10.;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = dpcb[2] * 2. + (583./2.-89.+26.8/2.)/10. - 10.5;
  // end-ladder elect 
  gMC->Gsvolu("RTB2", "TUBE", idtmed[274], dtube, 3);
  gMC->Gspos("RTB2", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dpgon[0] = 15.;
  dpgon[1] = 360.;
  dpgon[2] = 12.;
  dpgon[3] = 2.;
  dpgon[4] = dpcb[2] * 2. + (583./2.-62.2)/10. - 10.5;
  // end-ladder electron 
  dpgon[5] = 12.;
  dpgon[6] = 13.5;
  dpgon[7] = dpcb[2] * 2. + 583./2./10. - 10.5;
  // end-ladder electronics 
  dpgon[8] = 12.;
  dpgon[9] = 13.5;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gsvolu("RP03", "PGON", idtmed[274], dpgon, 10);
  gMC->Gspos("RP03", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dpgon[0] = 7.5;
  dpgon[1] = 360.;
  dpgon[2] = 24.;
  dpgon[3] = 2.;
  dpgon[4] = dpcb[2] * 2. + (583./2.+(338.-273.))/10. - 10.5;
  // end-ladder e 
  dpgon[5] = 21.;
  dpgon[6] = 23.;
  dpgon[7] = dpcb[2] * 2. + (583./2.+(338.-273.+15.))/10. - 10.5;
  // end-ladde 
  dpgon[8] = 21.;
  dpgon[9] = 23.;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gsvolu("RP04", "PGON", idtmed[274], dpgon, 10);
  gMC->Gspos("RP04", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  if (fMinorVersionV3 < 3 ) {
    offset2  = 5.2;
    dpgon[0] = offset2 + 360./(2.*35.);
    dpgon[1] = 360.;
    dpgon[2] = 35.;
    dpgon[3] = 2.;
    dpgon[4] = dpcb[2] * 2. + (583./2.+(338.-106.))/10. - 10.5;
    // end-ladde 
    dpgon[5] = 37.7;
    dpgon[6] = 40.;
    dpgon[7] = dpcb[2] * 2. + (583./2.+(338.-106.+15.))/10. - 10.5;
    // end-la 
    dpgon[8] = 37.7;
    dpgon[9] = 40.;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("RP05", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("RP05", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
    
    dpgon[0] = offset2 + 360./(2.*39.);
    dpgon[1] = 360.;
    dpgon[2] = 39.;
    dpgon[3] = 2.;
    dpgon[4] = dpcb[2] * 2. + (583./2.+(338.-56.))/10. - 10.5;
    // end-ladder 
    dpgon[5] = 42.7;
    dpgon[6] = 45.;
    dpgon[7] = dpcb[2] * 2. + (583./2.+(338.-56.+15.))/10. - 10.5;
    // end-la 
    dpgon[8] = 42.7;
    dpgon[9] = 45.;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("RP06", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("RP06", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  }
  if (fMinorVersionV3 > 2 && fMinorVersionV3 < 6) {
    offset2  = 5.2;
    dpgon[0] = offset2 + 5.625;
    dpgon[1] = 360.;
    dpgon[2] = 32.;
    dpgon[3] = 2.;
    dpgon[4] = (583./2.+(338.-106.))/10. - (40.-36.6) / TMath::Tan(45.*degrad) + dpcb[2] * 2. - 10.5;
    // end-ladder electronics 
    dpgon[5] = 34.3;
    dpgon[6] = 36.6;
    dpgon[7] = (583./2.+(338.-106.+15.))/10. - (40.-36.6) / TMath::Tan(45.*degrad) + dpcb[2] * 2. - 10.5;
    // end-ladder electr 
    dpgon[8] = 34.3;
    dpgon[9] = 36.6;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("RP05", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("RP05", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
    
    dpgon[0] = offset2 + 5.;
    dpgon[1] = 360.;
    dpgon[2] = 36.;
    dpgon[3] = 2.;
    dpgon[4] = (583./2.+(338.-56.))/10. - (45.-41.2) / TMath::Tan(45.*degrad) + dpcb[2] * 2. - 10.5;
    // end-ladder electronics 
    dpgon[5] = 38.9;
    dpgon[6] = 41.2;
    dpgon[7] = (583./2.+(338.-56.+15.))/10. - (45.-41.2) / TMath::Tan(45.*degrad) + dpcb[2] * 2. - 10.5;
    // end-ladder electr 
    dpgon[8] = 38.9;
    dpgon[9] = 41.2;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("RP06", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("RP06", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  }
  
  // 9345 CONTINUE 
  
  // --- Define the Z<0 end-cap 
  
  //      GOTO 9456             ! skip the Z<0 end-cap 
  
  dcone[0] = 16.75;
  dcone[1] = (338.-3.)*455./(338.-3.-10.)/10.;
  dcone[2] = .02 / TMath::Cos(45.*degrad) + (338.-3.)*455./(338.-3.-10.)/10.;
  dcone[3] = 12.;
  dcone[4] = 12.02;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = -(583.+(338.-3.))/2./10. - dpcb[2] * 2. + 10.5;
  // end-ladder electr 
  gMC->Gsvolu("LCON", "CONE", idtmed[274], dcone, 5);
  
  gMC->Gspos("LCON", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dtube[0] = .02 / TMath::Cos(45.*degrad) +  (338.-3.)*455./(338.-3.-10.)/10.;
  dtube[1] = 49.9;
  // In the Simonetti's drawings 52. In the TP 50. 
  dtube[2] = .15;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = -(583./2.+(338-1.5))/10. - dpcb[2] * 2. + 10.5;
  // end-ladder electr 
  gMC->Gsvolu("LTB1", "TUBE", idtmed[274], dtube, 3);
  
  gMC->Gspos("LTB1", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dtube[0] = 10.5;
  dtube[1] = 12.;
  dtube[2] = 26.8/2./10.;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = -(583./2.-89.+26.8/2.)/10. - dpcb[2] * 2. + 10.5;
  // end-ladder elec 
  gMC->Gsvolu("LTB2", "TUBE", idtmed[274], dtube, 3);
  ;
  gMC->Gspos("LTB2", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dpgon[0] = 15.;
  dpgon[1] = 360.;
  dpgon[2] = 12.;
  dpgon[3] = 2.;
  dpgon[4] = -583./2./10. - dpcb[2] * 2. + 10.5;
  // end-ladder electronics 
  dpgon[5] = 12.;
  dpgon[6] = 13.5;
  dpgon[7] = -(583./2.-62.2)/10. - dpcb[2] * 2. + 10.5;
  // end-ladder electro 
  dpgon[8] = 12.;
  dpgon[9] = 13.5;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gsvolu("LP03", "PGON", idtmed[274], dpgon, 10);
  gMC->Gspos("LP03", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  dpgon[0] = 7.5;
  dpgon[1] = 360.;
  dpgon[2] = 24.;
  dpgon[3] = 2.;
  dpgon[4] = -(583./2.+(338.-273.+15.))/10. - dpcb[2] * 2. + 10.5;
  // end-ladd 
  dpgon[5] = 21.;
  dpgon[6] = 23.;
  dpgon[7] = -(583./2.+(338.-273.))/10. - dpcb[2] * 2. + 10.5;
  // end-ladder 
  dpgon[8] = 21.;
  dpgon[9] = 23.;
  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gsvolu("LP04", "PGON", idtmed[274], dpgon, 10);
  gMC->Gspos("LP04", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  
  if (fMinorVersionV3 < 3) {
    offset2  = 5.2;
    dpgon[0] = offset2 + 360./(2.*35.);
    dpgon[1] = 360.;
    dpgon[2] = 35.;
    dpgon[3] = 2.;
    dpgon[4] = -(583./2.+(338.-106.))/10. - dpcb[2] * 2. + 10.5;
    // end-ladd 
    dpgon[5] = 37.7;
    dpgon[6] = 40.;
    dpgon[7] = -(583./2.+(338.-106.+15.))/10. - dpcb[2] * 2. + 10.5;
    // end-l 
    dpgon[8] = 37.7;
    dpgon[9] = 40.;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("LP05", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("LP05", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
    
    dpgon[0] = offset2 + 360./(2.*39.);
    dpgon[1] = 360.;
    dpgon[2] = 39.;
    dpgon[3] = 2.;
    dpgon[4] = -(583./2.+(338.-56.))/10. - dpcb[2] * 2. + 10.5;
    // end-ladde 
    dpgon[5] = 42.7;
    dpgon[6] = 45.;
    dpgon[7] = -(583./2.+(338.-56.+15.))/10. - dpcb[2] * 2. + 10.5;
    // end-l 
    dpgon[8] = 42.7;
    dpgon[9] = 45.;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("LP06", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("LP06", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  }
  if (fMinorVersionV3 > 2 && fMinorVersionV3 < 6) {
    offset2  = 5.2;
    dpgon[0] = offset2 + 5.625;
    dpgon[1] = 360.;
    dpgon[2] = 32.;
    dpgon[3] = 2.;
    dpgon[4] = (40.-36.6) / TMath::Tan(45.*degrad) - (583./2.+(338.-106.))/10. - dpcb[2] * 2. + 10.5;
    // end-ladder electronics 
    dpgon[5] = 34.3;
    dpgon[6] = 36.6;
    dpgon[7] = (40.-36.6) / TMath::Tan(45.*degrad) - (583./2.+(338.-106.+15.))/10. - dpcb[2] * 2. + 10.5;
    // end-ladder electr 
    dpgon[8] = 34.3;
    dpgon[9] = 36.6;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("LP05", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("LP05", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
    
    dpgon[0] = offset2 + 5.;
    dpgon[1] = 360.;
    dpgon[2] = 36.;
    dpgon[3] = 2.;
    dpgon[4] = (45.-41.2) / TMath::Tan(45.*degrad) - (583./2.+(338.-56.))/10. - dpcb[2] * 2. + 10.5;
    // end-ladder electronics 
    dpgon[5] = 38.9;
    dpgon[6] = 41.2;
    dpgon[7] = (45.-41.2) / TMath::Tan(45.*degrad) - (583./2.+(338.-56.+15.))/10. - dpcb[2] * 2. + 10.5;
    // end-ladder electr 
    dpgon[8] = 38.9;
    dpgon[9] = 41.2;
    xpos     = 0.;
    ypos     = 0.;
    zpos     = 0.;
    gMC->Gsvolu("LP06", "PGON", idtmed[274], dpgon, 10);
    gMC->Gspos("LP06", 1, "ITSV", xpos, ypos, zpos, 0, "ONLY");
  }
  
  // 9456 CONTINUE 
  
  
  // --- Outputs the geometry tree in the EUCLID/CAD format 
  
  if (fEuclidOut) {
    gMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
  }
  fMinorVersion = fMinorVersionV3;
}
 
//_____________________________________________________________________________
void AliITSv3::CreateMaterials()
{
  //
  // Create Materials for ITS
  //
   AliITS::CreateMaterials();
}

//_____________________________________________________________________________
void AliITSv3::Init(){
    //
    // Initialise its after it is built
    //
    Int_t i,j,l;

    fIdN       = fId3N;;
    fIdName    = new char*[fIdN];
    fIdSens    = new Int_t[fIdN];
    for(i=0;i<fId3N;i++) {
	l = strlen(fId3Name[i]);
	fIdName[i] = new char[l+1];
	for(j=0;j<l;j++) fIdName[i][j] = fId3Name[i][j];
	fIdName[i][l] = '\0'; // Null terminate this string.
    } // end for i
    //
    AliITS::Init();
    fMajorVersion = 3;
    fMinorVersion = fMinorVersionV3;
} 

//_____________________________________________________________________________
void AliITSv3::StepManager()
{
  //
  // Called at every step in ITS
  //
  Int_t         copy, id;
  Float_t       hits[8];
  Int_t         vol[4];
  TLorentzVector position, momentum;
  TClonesArray &lhits = *fHits;
  //
  // Track status
  vol[3] = 0;
  if(gMC->IsTrackInside())      vol[3] +=  1;
  if(gMC->IsTrackEntering())    vol[3] +=  2;
  if(gMC->IsTrackExiting())     vol[3] +=  4;
  if(gMC->IsTrackOut())         vol[3] +=  8;
  if(gMC->IsTrackDisappeared()) vol[3] += 16;
  if(gMC->IsTrackStop())        vol[3] += 32;
  if(gMC->IsTrackAlive())       vol[3] += 64;
  //
  // Fill hit structure.
  if(gMC->TrackCharge() && gMC->Edep()) {
    //
    // Only entering charged tracks
    if((id=gMC->CurrentVolID(copy))==fIdSens[0]) {  
      vol[0]=1;
      id=gMC->CurrentVolOffID(1,copy);      
      vol[1]=copy;
      id=gMC->CurrentVolOffID(2,copy);
      vol[2]=copy;                       
    } else if(id==fIdSens[1]) {
      vol[0]=2;
      id=gMC->CurrentVolOffID(1,copy);       
      vol[1]=copy;
      id=gMC->CurrentVolOffID(2,copy);
      vol[2]=copy;                    
    } else if(id==fIdSens[2]) {
      vol[0]=3;
      vol[1]=copy;
      id=gMC->CurrentVolOffID(1,copy);
      vol[2]=copy;             
    } else if(id==fIdSens[3]) {
      vol[0]=4;
      vol[1]=copy;
      id=gMC->CurrentVolOffID(1,copy);
      vol[2]=copy;                  
    } else if(id==fIdSens[4]) {
      vol[0]=5;
      vol[1]=copy;
      id=gMC->CurrentVolOffID(1,copy);
      vol[2]=copy;               
    } else if(id==fIdSens[5]) {
      vol[0]=6;
      vol[1]=copy;
      id=gMC->CurrentVolOffID(1,copy);
      vol[2]=copy;                      
    } else return;
    gMC->TrackPosition(position);
    gMC->TrackMomentum(momentum);
    hits[0]=position[0];
    hits[1]=position[1];
    hits[2]=position[2];          
    hits[3]=momentum[0];
    hits[4]=momentum[1];
    hits[5]=momentum[2];        
    hits[6]=gMC->Edep();
    hits[7]=gMC->TrackTime();
    new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }      
}

//____________________________________________________________________________
void AliITSv3::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliITSv3.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITS::Streamer(R__b);
      // This information does not need to be read. It is "hard wired"
      // into this class via its creators.
      //R__b >> fId3N;
      //R__b.ReadArray(fId3Name);
   } else {
      R__b.WriteVersion(AliITSv3::IsA());
      AliITS::Streamer(R__b);
      // This information does not need to be saved. It is "hard wired"
      // into this class via its creators.
      //R__b << fId3N;
      //R__b.WriteArray(fId3Name, __COUNTER__);
   }
}
