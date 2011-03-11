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

//-------------------------------------------------------------------------
//  Beam pipe class
//  Default version
//  Author: A.Morsch
//-------------------------------------------------------------------------

#include <Riostream.h>

#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TSystem.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliMagF.h"
#include "AliPIPEv0.h"
#include "AliRun.h"
#include "AliLog.h"
 
ClassImp(AliPIPEv0)

 
//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0():
    fPipeMaterial(kBe)   
{
// Constructor
}

//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0(const char *name, const char *title)
    : AliPIPE(name,title),
      fPipeMaterial(kBe)   
{
// Constructor
}

 
//___________________________________________
void AliPIPEv0::CreateGeometry()
{
//Begin_Html
/*
<img src="picts/pipe.gif">
*/
//End_Html


//Begin_Html
/*
<img src="picts/tree_pipe.gif">
*/
//End_Html

    AliDebugClass(1,"Create PIPEv0 geometry");
  

    Int_t *idtmed = fIdtmed->GetArray();
    Float_t ppcon[90], ptube[3], pbox[3];
    Int_t i=0;
    
    
    Int_t   idrotm[2099];  
    AliMatrix(idrotm[2001],90.,240.,  0.,  0., 90.,150.);
    AliMatrix(idrotm[2002],90.,  0.,  0.,  0., 90.,270.);
    AliMatrix(idrotm[2003],90.,120.,  0.,  0., 90., 30.);
    AliMatrix(idrotm[2004],90.,315., 90., 45.,  0.,  0.);
    AliMatrix(idrotm[2005],90.,270., 90.,  0.,  0.,  0.);
    AliMatrix(idrotm[2006],90.,225., 90.,315.,  0.,  0.);
    AliMatrix(idrotm[2007],90.,180., 90.,270.,  0.,  0.);
    AliMatrix(idrotm[2008],90.,135., 90.,225.,  0.,  0.);
    AliMatrix(idrotm[2009],90., 90., 90.,180.,  0.,  0.);
    AliMatrix(idrotm[2010],90., 45., 90.,135.,  0.,  0.);
    idrotm[2011] = 0;
    AliMatrix(idrotm[2012],90.,180., 90., 90.,180.,  0.);
    AliMatrix(idrotm[2013],90.,  0., 90., 90.,180.,  0.);
//
//  Bellow
//
//  distance between bellows
//  total size of bellow section
    const Float_t kdzb  = 14.6;
//  size of undulated region 
//  
//  Absorber side
//   
//  distance between bellows
    const Float_t kdzbbA =  4.6;
//  total size of bellow section
    const Float_t kdzbA  = 14.6;
//  size of undulated region 
    const Float_t kdzubA =  3.75;

// half-lengths of various beam pipe sections
// central Be-Pipe
    Float_t hlenQbbe1 = 40.;
    Float_t hlenQbbe2 = 36.5;
    Float_t hlenQbbe  = (hlenQbbe1+hlenQbbe2)/2.;
//
//
//    Float_t hlenQbt1 = 5.5/2.;
//
//  Pipe outside central region (non-absorber side)
    Float_t hlenQbab = 157.5 + 23./2.;
//
//  Flange non-absorber side
    Float_t hlenQb29 = 11.5/2.+1.75 + 5.0;
//
//  Bellow element 
    Float_t hlenQbe0 = kdzbA;
//
//  Inox pipe between Be and Bellow (absorber side)
    Float_t hlenQb24[3] = {11.3/2., 1.8, 3.3};
//
//
    Float_t hlenQb28 = (800.-hlenQbbe1-2.*hlenQbab-4.*hlenQb29-2.*hlenQbe0)/2.;
//
//  Position of the pump
    Float_t zPump = hlenQbbe1+2.*hlenQbab+2.*hlenQb29+kdzb;
//
//  Inner beam pipe radius
//  Be
    const Float_t kRinBe = 2.9;
//  Steel
    const Float_t kRinSt = 2.9;
//  Bellow
    const Float_t kRinSB = 2.92;
//
//  Outer beam pipe radius
//  Be
    const Float_t kRoutBe = 2.98;
//  Steel
    const Float_t kRoutSt = 2.98;
//  Bellow
    const Float_t kRoutSB = 3.00;

//
    Float_t dz;
    
//
// The peam pipe up to the Front Absorber
//
// Mother Volume QBPM
    ppcon[0]  =   0;
    ppcon[1]  = 360;
    ppcon[2]  =  20;
//  1 
    ppcon[3]  = -90.;
    ppcon[4]  =   0.;
    ppcon[5]  =   3.1;
//  2
    ppcon[6]  = -84.;
    ppcon[7]  =   0.;
    ppcon[8]  =   3.1;
//  3
    ppcon[9]  = -84.;
    ppcon[10] =   0.;
    ppcon[11] =   4.4;
//  4
    ppcon[12] = -90+2.*hlenQb24[2]+2.8+2.*hlenQb24[1];
    ppcon[13] =   0.;
    ppcon[14] =   4.4;
//  5
    ppcon[15]  = ppcon[12];
    ppcon[16] =   0.;
    ppcon[17] =   4.1;
//  6 
    ppcon[18] = ppcon[15] + 2.5 + 2.*kdzubA+0.2; 
    ppcon[19] =   0.;
    ppcon[20] =   4.1;
//  7 
    ppcon[21] = ppcon[18];
    ppcon[22] =   0.;
    ppcon[23] =   3.2;
//  8 
    ppcon[24] = ppcon[21] + 2.* kdzbbA-0.4; 
    ppcon[25] =   0.;
    ppcon[26] =   3.2;
//  9
    ppcon[27] = ppcon[24]; 
    ppcon[28] =   0.;
    ppcon[29] =   4.1;
//  10
    ppcon[30] = -44.;
    ppcon[31] =   0.;
    ppcon[32] =   4.1;
//  11
    ppcon[33] = -44.;
    ppcon[34] =    0;
    ppcon[35] =    3.06;
//  12
    ppcon[36] =  38.;
    ppcon[37] =    0;
    ppcon[38] =    3.06;
//  13
    ppcon[39] =  38.;
    ppcon[40] =    0;
    ppcon[41] =    4.1;
//  14
    ppcon[42] = hlenQbbe1+2.*hlenQbab-0.1;
    ppcon[43] =    0.;
    ppcon[44] =    4.1;
//  15
    ppcon[45] = ppcon[42];
    ppcon[46] =    0.;
    ppcon[47] =    4.1;
//  16
    ppcon[48] = ppcon[45]+2.*hlenQb29-5.;
    ppcon[49] =    0.;
    ppcon[50] =    4.1;
//  17
    ppcon[51] = ppcon[48];
    ppcon[52] =    0.;
    ppcon[53] =   56.;
//  18
    ppcon[54] = ppcon[51]+2.*kdzb+10.;
    ppcon[55] =    0.;
    ppcon[56] =   56.;
//  19
    ppcon[57] =   ppcon[54];
    ppcon[58] =    0.;
    ppcon[59] =    4.1;
//  20
    ppcon[60] =  800.;
    ppcon[61] =    0.;
    ppcon[62] =    4.1;
    
    gMC->Gsvolu("QBPM", "PCON", idtmed[kAir], ppcon,63);


//
// volume definitions of various sections
//

//
// The Vacuum 
    gMC->Gsvolu("QBVA","TUBE", idtmed[kVac], ptube, 0);
    ptube[0] =   0.0;
    ptube[1] =   kRinSt;
    ptube[2] =   (90.-hlenQbbe2)/2.;
    dz = -90. + ptube[2];
    gMC->Gsposp ("QBVA", 1, "QBPM", 0., 0., dz , 0, "ONLY", ptube, 3);
    dz = dz + ptube[2];

    ptube[1] =   kRinBe;
    ptube[2] =   hlenQbbe+hlenQbab;
    dz = dz + ptube[2];
    gMC->Gsposp ("QBVA", 2, "QBPM", 0., 0., dz , 0, "ONLY", ptube, 3);
    dz = dz + ptube[2];

    ptube[1] =   kRinSt;
    ptube[2] =   (800.-hlenQbbe1-2.*hlenQbab)/2.;
    dz = dz + ptube[2];
    gMC->Gsposp ("QBVA", 3, "QBPM", 0., 0., dz , 0, "ONLY", ptube, 3);

//
// Be Pipe in central Alice 
    ptube[0] = kRinBe;
    ptube[1] = kRoutBe;
    ptube[2] = hlenQbbe;
    
    gMC->Gsvolu("QBBE","TUBE", idtmed[kBe], ptube, 3);
    
//
//  Support Ring
//
    //  Mother
    ptube[0] = kRoutSB;
    ptube[1] = 4.0;
    ptube[2] = 0.6;
    gMC->Gsvolu("QBSR", "TUBE", idtmed[kAlu], ptube,3);
    // Inner support
    ptube[0] = kRoutSB;
    ptube[1] = 3.5;
    gMC->Gsvolu("QBSS", "TUBE", idtmed[kPA], ptube,3);
    gMC->Gspos("QBSS", 1, "QBSR", 0.0, 0.0,  0.0, 0, "ONLY");

    gMC->Gspos("QBSR", 1, "QBPM", 0.0, 0.0,  40., 0, "ONLY");
    gMC->Gspos("QBSR", 2, "QBPM", 0.0, 0.0, 150., 0, "ONLY");
    gMC->Gspos("QBSR", 3, "QBPM", 0.0, 0.0, 260., 0, "ONLY");
    gMC->Gspos("QBSR", 4, "QBPM", 0.0, 0.0,- 46., 0, "ONLY");
//
// Flange and Fixed Point: non absorber side
//
// ---------->
//
//  Mother
    ppcon[0]  =   0;
    ppcon[1]  = 360;
    ppcon[2]  =   4;
//  1: 
    ppcon[3]  = -hlenQb29;
    ppcon[4]  = kRinSt;
    ppcon[5]  = 5.8;
//  2
    ppcon[6]  = ppcon[3]+3.6;
    ppcon[7]  = kRinSt;
    ppcon[8]  = 5.8;
//  3
    ppcon[9]  = ppcon[6];
    ppcon[10] = kRinSt;
    ppcon[11] = 3.6;
//  4 
    ppcon[12] = hlenQb29;
    ppcon[13] = kRinSt;
    ppcon[14] = 3.6;
    
    gMC->Gsvolu("QB29", "PCON", idtmed[kAir], ppcon,15);
    

//    Flange
    ptube[0] = kRinSt;
    ptube[1] = 5.7;
    ptube[2] = 1.75;
    gMC->Gsvolu("QF29","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QF29", 1, "QB29", 0.0, 0.0, -hlenQb29+1.75, 0, "ONLY");
//    Pipe
    ptube[0] = kRinSt;
    ptube[1] = 3.06;
    ptube[2] = hlenQb29;
    gMC->Gsvolu("QS29","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QS29", 1, "QB29", 0.0, 0.0, 0., 0, "ONLY");
//    Fixed point
    ptube[0] = kRinSt;
    ptube[1] = 3.5;
    ptube[2] = 0.3;
    gMC->Gsvolu("QP29","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QP29", 1, "QB29", 0.0, 0.0, -hlenQb29+9.75+3., 0, "ONLY");
    
//
//
// Inox beam pipe: final section on non-absorber side

    ptube[0] =   kRinSt;
    ptube[1] =   kRoutSt;
    ptube[2] =   hlenQb28;    

    gMC->Gsvolu("QB28","TUBE", idtmed[kInox], ptube, 3);
//
//
// Air with high transport cuts outside QB28
    ptube[0] =   25.;
    ptube[1] =   100.;
    ptube[2] =   hlenQb28;    

    gMC->Gsvolu("QA28","TUBE", idtmed[kAirHigh], ptube, 3);
    gGeoManager->SetVolumeAttribute("QA28", "SEEN", 0);
//  Al-Be (40-60 wgt%, rho=2.7 g/cm**3) beam pipe
//
//  This section is under study (A.M. 1/2/2002)
//

    ptube[0] = kRinBe;
    if (fPipeMaterial == kAlu) {
	ptube[1] = 3.06;
    } else if (fPipeMaterial == kBe) {
	ptube[1] = kRoutBe;
    } else if (fPipeMaterial == kInox){
	ptube[1] = kRoutSt;
    }
    ptube[2] =   hlenQbab;    

    gMC->Gsvolu("QBAB","TUBE", idtmed[fPipeMaterial], ptube, 3);

// 2.5 mm thick SS tube for hanging pump
/*
    ptube[0] = Rin;
    ptube[1] = 3.15;
    ptube[2] = hlenQb26;
    
    gMC->Gsvolu("QB26","TUBE", idtmed[kInox], ptube, 3);
*/
//
// Bellows
//
//
// Mother Volume
    Float_t pconQBE0[33];
    pconQBE0[ 0]= 0;
    pconQBE0[ 1]= 360;
    pconQBE0[ 2]=  10;
//  1
    pconQBE0[ 3] = -kdzbA;
    pconQBE0[ 4] = kRinSB;
    pconQBE0[ 5] = kRoutSB;
//
    pconQBE0[ 6] = -kdzbA+2.5;
    pconQBE0[ 7] = kRinSB;
    pconQBE0[ 8] = kRoutSB;
//  2
    pconQBE0[ 9] = -kdzbA+2.5;
    pconQBE0[10] = kRinSB;
    pconQBE0[11] = 3.60;
//  3
    pconQBE0[12] = -kdzbA+2.5+2.*kdzubA;
    pconQBE0[13] = kRinSB;
    pconQBE0[14] = 3.60;
//  4
    pconQBE0[15] = -kdzbA+2.5+2.*kdzubA;
    pconQBE0[16] = kRinSB;
    pconQBE0[17] = kRoutSB;
//  5    
    pconQBE0[18] = -kdzbA+2.5+2.*kdzubA+2.*kdzbbA;
    pconQBE0[19] = kRinSB;
    pconQBE0[20] = kRoutSB;
//  6    
    pconQBE0[21] = -kdzbA+2.5+2.*kdzubA+2.*kdzbbA;
    pconQBE0[22] = kRinSB;
    pconQBE0[23] = 3.60;
//  7
    pconQBE0[24] = -kdzbA+2.5+4.*kdzubA+2.*kdzbbA;
    pconQBE0[25] = kRinSB;
    pconQBE0[26] = 3.60;
//  8
    pconQBE0[27] = -kdzbA+2.5+4.*kdzubA+2.*kdzbbA;
    pconQBE0[28] = kRinSB;
    pconQBE0[29] = kRoutSB;
// 9
    pconQBE0[30] = -kdzbA+5.0+4.*kdzubA+2.*kdzbbA;
    pconQBE0[31] = kRinSB;
    pconQBE0[32] = kRoutSB;

    gMC->Gsvolu("QBE0", "PCON", idtmed[kAir], pconQBE0, 33);
//
//  Undulated piece mother
    ptube[0] =  kRinSB;
    ptube[1] =  3.60;
    ptube[2] =  kdzubA;
    gMC->Gsvolu("QBEM","TUBE", idtmed[kAir], ptube, 3);
    dz = -kdzbA+kdzubA+2.5;
    gMC->Gspos("QBEM", 2 ,"QBE0", 0.0, 0.0,   dz, 0 , "ONLY");
    gMC->Gspos("QBEM", 1 ,"QBE0", 0.0, 0.0,  -dz, idrotm[2012], "ONLY");
//  
    Float_t pund[30];
    Float_t uw = 0.02;
    Float_t ur = 2.*kdzubA-36.*uw;
    Float_t uz = ur/37.;
    Float_t ut = uz+uw;
    
    pund[ 0] = 0;
    pund[ 1] = 360;
    pund[ 2] =  8;

    pund[ 3] = -ut;
    pund[ 4] = kRinSB;
    pund[ 5] = kRinSB+uw;
    
    pund[ 6] = -ut+uz;
    pund[ 7] = pund[4];
    pund[ 8] = pund[5];
    
    pund[ 9] = pund[6];
    pund[10] = pund[4];
    pund[11] = 3.6;
    
    pund[12] = pund[9]+uw;
    pund[13] = pund[10];
    pund[14] = pund[11];
    
    pund[15] = pund[12];
    pund[16] = 3.6-uw;
    pund[17] = pund[14];
    
    pund[18] = pund[12]+uz;
    pund[19] = pund[16];
    pund[20] = pund[17];
    
    pund[21] = pund[18];
    pund[22] = kRinSB;
    pund[23] = pund[20];
    
    pund[24] = pund[21]+uw;
    pund[25] = pund[22];
    pund[26] = pund[23];

    gMC->Gsvolu("QBEU", "PCON", idtmed[kInox], pund, 27);

    for (i = 0; i < 18; i++)
    {
	dz = -kdzubA+(1+2*i)*ut;
	gMC->Gspos("QBEU", i+1 ,"QBEM", 0.0, 0.0,   dz, 0 , "ONLY");
    }
    ptube[0] =  kRinSB;
    ptube[1] =  kRinSB+uw;
    ptube[2] =  uz;
    gMC->Gsvolu("QBEW","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QBEW", 1 ,"QBEM", 0.0, 0.0,   kdzubA-uz, 0 , "ONLY");

//
//  BeamPipe
    gMC->Gsvolu("QBEP","TUBE", idtmed[kInox], ptube, 0);
    ptube[0] =  kRinSB;
    ptube[1] =  kRoutSB;
    ptube[2] =  1.25;
    gMC->Gsposp("QBEP", 1 ,"QBE0", 0.0, 0.0, -kdzbA+1.25, 0 , "ONLY", ptube, 3);
    gMC->Gsposp("QBEP", 2 ,"QBE0", 0.0, 0.0,  kdzbA-1.25, 0 , "ONLY", ptube, 3);    
    ptube[2] = kdzbbA;
    gMC->Gsposp("QBEP", 3 ,"QBE0", 0.0, 0.0,  0., 0 , "ONLY", ptube, 3);    
//  
//
//  ----> End Bellow
//
// **** Placement of various sections on non-absorber side ****
//
    //
    // first the beryllium section
    Float_t zpos = -(hlenQbbe2-hlenQbbe1)/2;
    gMC->Gspos("QBBE", 1, "QBPM", 0., 0., zpos, 0, "ONLY");

    // next meta-metal transition QBT1 on on-absorber side
//    zpos = zpos + hlenQbbe + hlenQbt1;
//    gMC->Gspos("QBT1", 1, "QBPM", 0., 0.,  zpos, 0, "ONLY");
    
    // Aluminium OR Al-be alloy section
    zpos = hlenQbbe1+hlenQbab;
    gMC->Gspos("QBAB", 1, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");
    //
    // inox flange at the start of bellow
    zpos = zpos + hlenQbab + hlenQb29;
    gMC->Gspos("QB29", 1, "QBPM", 0.0, 0.0, zpos, idrotm[2012], "ONLY");
    //
    // bellow section
    zpos = zpos + hlenQb29 + hlenQbe0;
    gMC->Gspos("QBE0", 2 ,"QBPM", 0.0, 0.0, zpos, 0, "ONLY");
    // 
    // inox flange at the end of bellow and start of thick inox for pump
    zpos = zpos + hlenQbe0 + hlenQb29;
    gMC->Gspos("QB29", 2, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");
    //   
    //last inox section till 800 cm
    zpos = zpos + hlenQb29 + hlenQb28;
    gMC->Gspos("QB28", 1, "QBPM", 0.0, 0.0, zpos, 0, "ONLY"); 
    gMC->Gspos("QA28", 1, "ALIC", 0.0, 0.0, zpos, 0, "ONLY"); 
    
//******** end of placement on non-absorber side *********
    //
    // **** Absorber side *****   
    //
    //
//
// Beam pipes between elements
//    
    gMC->Gsvolu("QB24","TUBE", idtmed[kInox], ptube, 0);
    ptube[0] = kRinSt;
    ptube[1] = kRoutSt;
    ptube[2] = hlenQb24[0];
    Float_t bpbe[33];
    bpbe[ 0] =   0.;
    bpbe[ 1] = 360.;
    bpbe[ 2] =  10.;
    // 1
    bpbe[ 3] = -hlenQb24[0];
    bpbe[ 4] = kRinSB;
    bpbe[ 5] = kRoutSB;
    // 2
    bpbe[ 6] = hlenQb24[0] - 5.8;
    bpbe[ 7] = kRinSB;
    bpbe[ 8] = kRoutSB;
    // 3
    bpbe[ 9] = hlenQb24[0] - 5.8;
    bpbe[10] = kRinSt;
    bpbe[11] = 3.05;
    // 4
    bpbe[12] = hlenQb24[0] - 4.5;
    bpbe[13] = kRinSt;
    bpbe[14] = 3.05;
    // 5
    bpbe[15] = hlenQb24[0] - 4.5;
    bpbe[16] = kRinSt;
    bpbe[17] = kRoutSt ;
    // 6
    bpbe[18] = hlenQb24[0] - 3.5;
    bpbe[19] = kRinSt;
    bpbe[20] = kRoutSt ;
    // 7
    bpbe[21] = hlenQb24[0] - 3.5;
    bpbe[22] = kRinSt;
    bpbe[23] = 3.05;
    // 8
    bpbe[24] = hlenQb24[0] - 3.0;
    bpbe[25] = kRinSt;
    bpbe[26] = 3.05;
    // 9
    bpbe[27] = hlenQb24[0] - 3.0;
    bpbe[28] = kRinSt;
    bpbe[29] = kRoutSt;
    // 10
    bpbe[30] = hlenQb24[0];
    bpbe[31] = kRinSt;
    bpbe[32] = kRoutSt;

    gMC->Gsvolu("QA24","PCON", idtmed[kInox], bpbe, 33);

    dz = hlenQbbe2 + hlenQb24[0];
    
    gMC->Gspos("QA24", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY");
//
// Bellow on absorber side

    dz = dz+hlenQb24[0] + kdzbA;
    gMC->Gspos("QBE0", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY");
//
    ptube[2] = hlenQb24[1];
    dz = dz + kdzb + ptube[2];
    gMC->Gsposp("QB24", 2 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY", ptube, 3);
    dz = dz + ptube[2];
    
//
// Flange
// 
//  Mother Volume
    ptube[0] = kRinSB;
    ptube[1] = 4.300;
    ptube[2] = 1.4;
    
    gMC->Gsvolu("QFA0","TUBE", idtmed[kInox], ptube, 3);
    dz = dz + ptube[2];
    gMC->Gspos("QFA0", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY");
    dz = dz + ptube[2];
//
//
    ptube[0] = kRinSB;
    ptube[1] = kRoutSB;
    ptube[2] = hlenQb24[2];
    dz = dz + ptube[2];
    gMC->Gsposp("QB24", 3 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY", ptube, 3);

// --- Place the PIPE ghost volume (QBPM) in its mother volume (ALIC)
//    by rotating it to 180 deg. and make it invisible
// 
    gMC->Gspos("QBPM",1,"ALIC", 0, 0, 0, 0, "ONLY");
    gMC->Gsbool("QBPM", "L3DX");
    gMC->Gsbool("QBPM", "L3O3");
    gMC->Gsbool("QBPM", "L3O4");


//
// ******** Ion Pump volume description starts here ******
// 
    //
    // Getters ->
    pbox[0] =  6.50;
    pbox[1] =  6.75;
    pbox[2] = 15.60;
    gMC->Gsvolu("QI32","BOX", idtmed[kInox], pbox, 3);
    
    pbox[0] =  5.90;
    pbox[1] =  6.15;
    pbox[2] = 15.00;
    gMC->Gsvolu("QI42","BOX", idtmed[kGetter], pbox, 3);
    gMC->Gspos("QI42", 1, "QI32", 0.0, 0.0, 0.0, 0, "ONLY");
// <-

    ptube[0] =  0.0;
    ptube[1] = 19.0;
    ptube[2] =  2.5;
    gMC->Gsvolu("QI33","TUBE", idtmed[kInox], ptube, 3);


    ptube[0] =  0.0;
    ptube[1] = 15.0;
    ptube[2] =  2.5;
    gMC->Gsvolu("QI43","TUBE", idtmed[kAir], ptube, 3);
    gMC->Gspos("QI43", 1, "QI33", 0.0, 0.0, 0.0, 0, "ONLY");
// 
// Connecting tube ->
    ptube[0] =  0.0;
    ptube[1] =  4.5;
    ptube[2] = 14.6;
    gMC->Gsvolu("QI34","TUBE", idtmed[kInox], ptube, 3);
    
    ptube[0] =  0.0;
    ptube[1] =  3.9;
    ptube[2] = 14.6;
    gMC->Gsvolu("QI44","TUBE", idtmed[kAir], ptube, 3);
    gMC->Gspos("QI44", 1, "QI34", 0.0, 0.0, 0.0, 0, "ONLY");
// <-

  //
  // Flange ->
    ptube[0] =  4.6;
    ptube[1] =  7.30;
    ptube[2] =  2.15;
    gMC->Gsvolu("QI35","TUBE", idtmed[kInox], ptube, 3);
// <-
    gMC->Gspos("QI32", 1, "QBPM", 0.0, -44.25, zPump, 0, "ONLY");
    gMC->Gspos("QI33", 1, "QBPM", 0.0, -35.00, zPump,idrotm[2002], "ONLY");
    gMC->Gspos("QI34", 1, "QBPM", 0.0, -17.90, zPump,idrotm[2002], "ONLY");
    gMC->Gspos("QI35", 1, "QBPM", 0.0, -24.35, zPump,idrotm[2002], "ONLY");

    gGeoManager->SetVolumeAttribute("QBPM", "SEEN", 1);
    gGeoManager->SetVolumeAttribute("QBEM", "SEEN", 1);
}



//___________________________________________
void AliPIPEv0::CreateMaterials()
{
  //
  // Define materials for beam pipe
  //

  AliDebugClass(1,"Create PIPEv0 materials");
  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  // Steel (Inox)  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  // AlBe - alloy 
  Float_t aAlBe[2] = { 26.98, 9.01};
  Float_t zAlBe[2] = { 13.00, 4.00};
  Float_t wAlBe[2] = { 0.4, 0.6};
  //
  // Polyamid
  Float_t aPA[4] = {16., 14., 12.,  1.};
  Float_t zPA[4] = { 8.,  7.,  6.,  1.};
  Float_t wPA[4] = { 1.,  1.,  6., 11.};
  //
  // Air 
  //
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  Float_t dAir1 = 1.20479E-10;
  //
  // Kapton
  //
  Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
  Float_t zKapton[4]={1.,6.,7.,8.};
  Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
  Float_t dKapton = 1.42;
  //
  //     Berillium 
  AliMaterial(5, "BERILLIUM$", 9.01, 4., 1.848, 35.3, 36.7);
  //
  //     Carbon 
  AliMaterial(6,  "CARBON$   ", 12.01, 6., 2.265, 18.8, 49.9);
  //
  //     Aluminum 
  AliMaterial(9,  "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  //
  //     Air 
  AliMixture(15, "AIR$",      aAir, zAir, dAir, 4, wAir);
  AliMixture(35, "AIR_HIGH$", aAir, zAir, dAir, 4, wAir);
  //
  //     Vacuum 
  AliMixture(16, "VACUUM$ ", aAir, zAir, dAir1, 4, wAir);
  //
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  //
  //     reduced density steel to approximate pump getter material
  AliMixture(20, "GETTER$", asteel, zsteel, 1.00, 4, wsteel);
  //     Al-Be alloy
  //     
  AliMixture(21, "AlBe$", aAlBe, zAlBe, 2.07, 2, wAlBe);
  //     Polyamid
  //   
  AliMixture(22, "PA$", aPA, zPA, 1.14, -4, wPA);
  //
  //     Kapton
  AliMixture(23, "KAPTON", aKapton, zKapton, dKapton, 4, wKapton);
  //
  // **************** 
  //     Defines tracking media parameters. 
  //
  Float_t epsil  = .001;    // Tracking precision, 
  Float_t stemax = -0.01;   // Maximum displacement for multiple scat 
  Float_t tmaxfd = -20.;    // Maximum angle due to field deflection 
  Float_t deemax = -.3;     // Maximum fractional energy loss, DLS 
  Float_t stmin  = -.8;
  // *************** 
  //
  //    Beryllium 
  
  AliMedium(5, "BE",       5, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    Carbon 
  AliMedium(6, "C",        6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Aluminum 
  AliMedium(9, "ALU",      9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Air 
  AliMedium(15, "AIR",     15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_HIFG",35, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Vacuum 
  AliMedium(16, "VACUUM", 16, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Steel 
  AliMedium(19, "INOX",   19, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Getter 
  AliMedium(20, "GETTER", 20, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   AlBe - Aloy 
  AliMedium(21, "AlBe"  , 21, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   Polyamid
  AliMedium(22, "PA"  ,   22, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   KAPTON
  AliMedium(23, "KAPTON", 23, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

}










