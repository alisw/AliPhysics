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
Revision 1.19.4.1  2002/06/10 15:13:48  hristov
Merged with v3-08-02

Revision 1.21  2002/05/11 19:14:44  morsch
PCONE QBEU had last z defined twice.

Revision 1.20  2002/05/02 12:36:43  morsch
New version of beam-pipe geometry. Ref. AL502206PL
(used in production readiness report).

Revision 1.19  2002/03/25 12:34:40  morsch
Obsolete support rollers removed.

Revision 1.18  2002/02/01 18:02:41  morsch
Material of beam pipe between inner Be piece and forward detectors
can be set by SetPipeMaterial(mat), mat = kInox, kAlu, kBe ...

Revision 1.17  2001/09/24 13:11:50  morsch
Ion pump and bellows moved out by 15 cm to make space for forward
detectors.

Revision 1.16  2001/05/16 14:57:22  alibrary
New files for folders and Stack

Revision 1.15  2001/05/02 11:50:18  morsch
New layout of the non-absorber side provided by Y. Viyogi. Not the final design
but the prsent most realistic.

Revision 1.14  2001/01/20 16:56:33  morsch
Put air in connecting tubes and flanges of vacuum pump.

Revision 1.13  2001/01/20 16:35:27  morsch
Increase mother volume for bellows.

Revision 1.12  2000/12/21 16:41:06  morsch
Coding convention clean-up (RS3)

Revision 1.11  2000/11/28 16:06:57  morsch
Undulated beam-pipe replaced by Al-Be (40,60) pipe 1.5 mm thick.

Revision 1.10  2000/11/24 13:00:37  morsch
- Geometry and materials imported from euclid output
- include comments
- better struturing of volume tree
- improved version of flange close to front absorber
- more realistic pump materials
- undulated beam pipe imported from v3.

Revision 1.9  2000/10/02 21:28:15  fca
Removal of useless dependecies via forward declarations

Revision 1.8  2000/06/11 12:37:01  morsch
Coding rule violations corrected

Revision 1.7  2000/02/23 16:25:24  fca
AliVMC and AliGeant3 classes introduced
ReadEuclid moved from AliRun to AliModule

Revision 1.6  1999/09/29 09:24:30  fca
Introduction of the Copyright and cvs Log

*/

////////////////////////////////////////////////
//  Beam pipe class                            /
////////////////////////////////////////////////

#include "AliPIPEv0.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliMagF.h"
#include "AliMC.h"
#include "TSystem.h"

#include <iostream.h>
 
ClassImp(AliPIPEv0)
 
//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0()
{
// Constructor
    SetPipeMaterial();
}

//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0(const char *name, const char *title)
  : AliPIPE(name,title)
{
// Constructor
    SetPipeMaterial();
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

    if(fDebug) printf("%s: Create PIPEv0 geometry \n",ClassName());
  

    Int_t *idtmed = fIdtmed->GetArray();
    Float_t ppcon[84], ptube[3], pbox[3];
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
    const Float_t dzbb =  9.0;
//  total size of bellow section
    const Float_t dzb  = 15.0;
//  size of undulated region 
    const Float_t dzub =  2.0;

// half-lengths of various beam pipe sections
// central Be-Pipe
    Float_t hlenQbbe1 = 40.;
    Float_t hlenQbbe2 = 36.5;
    Float_t hlenQbbe  = (hlenQbbe1+hlenQbbe2)/2.;
//
//
//    Float_t hlenQbt1 = 5.5/2.;
//
//  Pipe outside central region (non-absober side)
    Float_t hlenQbab = 157.5;
//
//  Flange non-absorber side
    Float_t hlenQb29 = 11.5/2.+1.75 + 5.0;
//
//  Bellow element 
    Float_t hlenQbe0 = dzb;
//
//  Inox pipe between Be and Bellow (absorber side)
    Float_t hlenQb24[3] = {10.5/2., 1.8, 3.3};
//
//
    Float_t hlenQb28 = (800.-hlenQbbe1-2.*hlenQbab-4.*hlenQb29-2.*hlenQbe0)/2.;
//
//  Position of the pump
    Float_t zPump = hlenQbbe1+2.*hlenQbab+2.*hlenQb29+dzb;
//
//  Inner beam pipe radius
    Float_t RinBe = 2.9;
    Float_t RinSt = 2.92;
//
//
    Float_t RoutBe = 2.98;
    Float_t RoutSt = 3.00;


//
    Float_t dz;
    
//
// The peam pipe up to the Front Absorber
//
// Mother Volume QBPM
    ppcon[0]  =   0;
    ppcon[1]  = 360;
    ppcon[2]  =  18;
//  1: 
    ppcon[3]  = -90.;
    ppcon[4]  =   0.;
    ppcon[5]  =   4.4;
//  2
    ppcon[6]  = -90+2.*hlenQb24[2]+2.8+2.*hlenQb24[1];
    ppcon[7]  =   0.;
    ppcon[8]  =   4.4;
//  3
    ppcon[9]  = ppcon[6];
    ppcon[10] =   0.;
    ppcon[11] =   4.1;
//  4 
    ppcon[12] = ppcon[9] + 2. + 2.*dzub+0.2; 
    ppcon[13] =   0.;
    ppcon[14] =   4.1;
//  5 
    ppcon[15] = ppcon[12];
    ppcon[16] =   0.;
    ppcon[17] =   3.2;
//  6 
    ppcon[18] = ppcon[15] + 2.* dzbb-0.4; 
    ppcon[19] =   0.;
    ppcon[20] =   3.2;
//  7
    ppcon[21] = ppcon[18]; 
    ppcon[22] =   0.;
    ppcon[23] =   4.1;
//  8
    ppcon[24] = -44.;
    ppcon[25] =   0.;
    ppcon[26] =   4.1;
//  9
    ppcon[27] = -44.;
    ppcon[28] =    0;
    ppcon[29] =    3.0;
//  10
    ppcon[30] =  38.;
    ppcon[31] =    0;
    ppcon[32] =    3.0;
//  11
    ppcon[33] =  38.;
    ppcon[34] =    0;
    ppcon[35] =    3.6;
//  12
    ppcon[36] = hlenQbbe1+2.*hlenQbab-0.1;
    ppcon[37] =    0.;
    ppcon[38] =    3.6;
//  13
    ppcon[39] = ppcon[36];
    ppcon[40] =    0.;
    ppcon[41] =    3.6;
//  14
    ppcon[42] = ppcon[39]+2.*hlenQb29-5.;
    ppcon[43] =    0.;
    ppcon[44] =    3.6;
//  15
    ppcon[45] = ppcon[42];
    ppcon[46] =    0.;
    ppcon[47] =   56.;
//  16
    ppcon[48] = ppcon[45]+2.*dzb+10.;
    ppcon[49] =    0.;
    ppcon[50] =   56.;
//  17
    ppcon[51] =   ppcon[48];
    ppcon[52] =    0.;
    ppcon[53] =    3.6;
//  18
    ppcon[54] =  800.;
    ppcon[55] =    0.;
    ppcon[56] =    3.6;
    
    gMC->Gsvolu("QBPM", "PCON", idtmed[kAir], ppcon,57);


//
// volume definitions of various sections
//

//
// The Vacuum 
    gMC->Gsvolu("QBVA","TUBE", idtmed[kVac], ptube, 0);
    ptube[0] =   0.0;
    ptube[1] =   RinSt;
    ptube[2] =   (90.-hlenQbbe2)/2.;
    dz = -90. + ptube[2];
    gMC->Gsposp ("QBVA", 1, "QBPM", 0., 0., dz , 0, "ONLY", ptube, 3);
    dz = dz + ptube[2];

    ptube[1] =   RinBe;
    ptube[2] =   hlenQbbe+hlenQbab;
    dz = dz + ptube[2];
    gMC->Gsposp ("QBVA", 2, "QBPM", 0., 0., dz , 0, "ONLY", ptube, 3);
    dz = dz + ptube[2];

    ptube[1] =   RinSt;
    ptube[2] =   (800.-hlenQbbe1-2.*hlenQbab)/2.;
    dz = dz + ptube[2];
    gMC->Gsposp ("QBVA", 3, "QBPM", 0., 0., dz , 0, "ONLY", ptube, 3);

//
// Be Pipe in central Alice 
    ptube[0] = RinBe;
    ptube[1] = RoutBe;
    ptube[2] = hlenQbbe;
    
    gMC->Gsvolu("QBBE","TUBE", idtmed[kBe], ptube, 3);
    
//
//  Support Ring
//
    //  Mother
    ppcon[0]  =   0;
    ppcon[1]  = 360;
    ppcon[2]  =   6;
//  1: 
    ppcon[3]  = -1.;
    ppcon[4]  = 3.0;
    ppcon[5]  = 3.4;
//  2
    ppcon[6]  = -0.8;
    ppcon[7]  = 3.0;
    ppcon[8]  = 3.4;
//  3
    ppcon[9]  = -0.8;
    ppcon[10] = 3.0;
    ppcon[11] = 3.2;
//  4 
    ppcon[12] = 0.8;
    ppcon[13] = 3.0;
    ppcon[14] = 3.2;
//  5 
    ppcon[15] = 0.8;
    ppcon[16] = 3.0;
    ppcon[17] = 3.4;
//  6 
    ppcon[18] = 1.0;
    ppcon[19] = 3.0;
    ppcon[20] = 3.4;

    
    gMC->Gsvolu("QBSR", "PCON", idtmed[kC], ppcon,21);
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
    ppcon[4]  = RinSt;
    ppcon[5]  = 5.8;
//  2
    ppcon[6]  = ppcon[3]+3.6;
    ppcon[7]  = RinSt;
    ppcon[8]  = 5.8;
//  3
    ppcon[9]  = ppcon[6];
    ppcon[10] = RinSt;
    ppcon[11] = 3.6;
//  4 
    ppcon[12] = hlenQb29;
    ppcon[13] = RinSt;
    ppcon[14] = 3.6;
    
    gMC->Gsvolu("QB29", "PCON", idtmed[kAir], ppcon,15);
    

//    Flange
    ptube[0] = RinSt;
    ptube[1] = 5.7;
    ptube[2] = 1.75;
    gMC->Gsvolu("QF29","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QF29", 1, "QB29", 0.0, 0.0, -hlenQb29+1.75, 0, "ONLY");
//    Pipe
    ptube[0] = RinSt;
    ptube[1] = 3.0;
    ptube[2] = hlenQb29;
    gMC->Gsvolu("QS29","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QS29", 1, "QB29", 0.0, 0.0, 0., 0, "ONLY");
//    Fixed point
    ptube[0] = RinSt;
    ptube[1] = 3.5;
    ptube[2] = 0.3;
    gMC->Gsvolu("QP29","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QP29", 1, "QB29", 0.0, 0.0, -hlenQb29+9.75+3., 0, "ONLY");
    
//
//
// Inox beam pipe: final section on non-absorber side

    ptube[0] =   RinSt;
    ptube[1] =   RoutSt;
    ptube[2] =   hlenQb28;    

    gMC->Gsvolu("QB28","TUBE", idtmed[kInox], ptube, 3);


//  Al-Be (40-60 wgt%, rho=2.7 g/cm**3) beam pipe
//
//  This section is under study (A.M. 1/2/2002)
//

    ptube[0] = RinBe;
    if (fPipeMaterial == kAlu) {
	ptube[1] = 3.0;
    } else if (fPipeMaterial == kBe) {
	ptube[1] = RoutBe;
    } else if (fPipeMaterial == kInox){
	ptube[1] = RoutSt;
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
    pconQBE0[ 2]= 10;
//  1
    pconQBE0[ 3]= -dzb;
    pconQBE0[ 4]= RinSt;
    pconQBE0[ 5]= RoutSt;
//  2    
    pconQBE0[ 6]= -dzb+2.;
    pconQBE0[ 7]= RinSt;
    pconQBE0[ 8]= RoutSt;
//  3
    pconQBE0[ 9]= -dzb+2.;
    pconQBE0[10]= RinSt;
    pconQBE0[11]= 4.00;
//  4
    pconQBE0[12]= -dzb+2.+2.*dzub;
    pconQBE0[13]= RinSt;
    pconQBE0[14]= 4.00;
//  5
    pconQBE0[15]= -dzb+2.+2.*dzub;
    pconQBE0[16]= RinSt;
    pconQBE0[17]= RoutSt;
//  6    
    pconQBE0[18]= -dzb+2.+2.*dzub+2.*dzbb;
    pconQBE0[19]= RinSt;
    pconQBE0[20]= RoutSt;
//  7    
    pconQBE0[21]= -dzb+2.+2.*dzub+2.*dzbb;
    pconQBE0[22]= RinSt;
    pconQBE0[23]= 4.00;
//  8
    pconQBE0[24]= -dzb+2.+4.*dzub+2.*dzbb;
    pconQBE0[25]= RinSt;
    pconQBE0[26]= 4.00;
//  9
    pconQBE0[27]= -dzb+2.+4.*dzub+2.*dzbb;
    pconQBE0[28]= RinSt;
    pconQBE0[29]= RoutSt;
//  10 
    pconQBE0[30]= +dzb;
    pconQBE0[31]= RinSt;
    pconQBE0[32]= RoutSt;
    gMC->Gsvolu("QBE0", "PCON", idtmed[kAir], pconQBE0, 33);
//
//  Undulated piece mother
    ptube[0] =  RinSt;
    ptube[1] =  4.00;
    ptube[2] =  dzub;
    gMC->Gsvolu("QBEM","TUBE", idtmed[kAir], ptube, 3);
    dz = -dzb+2.+dzub;
    gMC->Gspos("QBEM", 2 ,"QBE0", 0.0, 0.0,   dz, 0 , "ONLY");
    gMC->Gspos("QBEM", 1 ,"QBE0", 0.0, 0.0,  -dz, idrotm[2012], "ONLY");
//  
    Float_t pund[30];
    Float_t uw = 0.02;
    Float_t ur = 2.*dzub-12.*uw;
    Float_t uz = ur/13.;
    Float_t ut = uz+uw;
    
    pund[ 0] = 0;
    pund[ 1] = 360;
    pund[ 2] =  8;

    pund[ 3] = -ut;
    pund[ 4] = RinSt;
    pund[ 5] = RinSt+uw;
    
    pund[ 6] = -ut+uz;
    pund[ 7] = pund[4];
    pund[ 8] = pund[5];
    
    pund[ 9] = pund[6];
    pund[10] = pund[4];
    pund[11] = 4.0;
    
    pund[12] = pund[9]+uw;
    pund[13] = pund[10];
    pund[14] = pund[11];
    
    pund[15] = pund[12];
    pund[16] = 4.0-uw;
    pund[17] = pund[14];
    
    pund[18] = pund[12]+uz;
    pund[19] = pund[16];
    pund[20] = pund[17];
    
    pund[21] = pund[18];
    pund[22] = RinSt;
    pund[23] = pund[20];
    
    pund[24] = pund[21]+uw;
    pund[25] = pund[22];
    pund[26] = pund[23];

    gMC->Gsvolu("QBEU", "PCON", idtmed[kInox], pund, 27);

    for (i = 0; i < 6; i++)
    {
	dz = -dzub+(1+2*i)*ut;
	gMC->Gspos("QBEU", i+1 ,"QBEM", 0.0, 0.0,   dz, 0 , "ONLY");
    }
    ptube[0] =  RinSt;
    ptube[1] =  RinSt+uw;
    ptube[2] =  uz;
    gMC->Gsvolu("QBEW","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QBEW", 1 ,"QBEM", 0.0, 0.0,   dzub-uz, 0 , "ONLY");
//
//  BeamPipe
    gMC->Gsvolu("QBEP","TUBE", idtmed[kInox], ptube, 0);
    ptube[0] =  RinSt;
    ptube[1] =  RoutSt;
    ptube[2] =  1.;
    gMC->Gsposp("QBEP", 1 ,"QBE0", 0.0, 0.0, -dzb+1., 0 , "ONLY", ptube, 3);
    gMC->Gsposp("QBEP", 2 ,"QBE0", 0.0, 0.0,  dzb-1., 0 , "ONLY", ptube, 3);    
    ptube[2] = dzbb;
    gMC->Gsposp("QBEP", 3 ,"QBE0", 0.0, 0.0,  0., 0 , "ONLY", ptube, 3);    
//  
//  End undulated part
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
    
//******** end of placement on no-absorber side *********

    //
    // **** Absorber side *****   
    //
    //
/*
    //  metal-metal transition :  Be-Alu on absorber side
    //  Mother Volume
    ptube[0] = 2.900;
    ptube[1] = 4.200;
    ptube[2] = 2.750;
    gMC->Gsvolu("QBT2","TUBE", idtmed[kAir], ptube, 3);
    // z  = 43.3 - 48.8
    gMC->Gspos("QBT2", 1, "QBPM", 0., 0., -hlenQbbe-ptube[2], idrotm[2012], "ONLY");

    ptube[0] = 2.900;
    ptube[1] = 3.150;
    ptube[2] = 0.375;
    //
    //  Be-part
    gMC->Gsvolu("QB02","TUBE", idtmed[kAlu], ptube, 3);

    ptube[1] = 3.000;
    gMC->Gsvolu("QBA2","TUBE", idtmed[kBe], ptube, 3);

    gMC->Gspos("QBA2", 1, "QB02", 0., 0., 0, 0, "ONLY");
//  z = -2.75 -> -2.00
    gMC->Gspos("QB02", 1, "QBT2", 0., 0.,-2.75+ptube[2], 0, "ONLY");

    // Alu part    
    ptube[0] = 2.900;
    ptube[1] = 3.150;
    ptube[2] = 2.375;
// z = -2.00 -> 2.75
    gMC->Gsvolu("QB04","TUBE", idtmed[kAlu], ptube, 3);
    gMC->Gspos("QB04", 1, "QBT2", 0., 0.,-2.+ptube[2], 0, "ONLY");
    
    
    ptube[0] = 3.15;
    ptube[1] = 3.50;
    ptube[2] = 0.10;
// z = 2.55 -> 2.75
    gMC->Gsvolu("QB06","TUBE", idtmed[kAlu], ptube, 3);
    gMC->Gspos("QB06", 1, "QBT2", 0., 0., 2.55+ptube[2], 0, "ONLY");
    
    
    // Fixation
    ptube[0] = 0.0;
    ptube[1] = 0.1;
    ptube[2] = 0.5;
    
    gMC->Gsvolu("QBA8","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QBA8", 1 ,"QBT2",  0.000,  3.650, -1.25, idrotm[2002], "ONLY");
    gMC->Gspos("QBA8", 2 ,"QBT2",  3.161, -1.825, -1.25, idrotm[2001], "ONLY");
    gMC->Gspos("QBA8", 3 ,"QBT2", -3.161, -1.825, -1.25, idrotm[2003], "ONLY");
    
    // Carbon ring
    ptube[0] = 3.15;
    ptube[1] = 4.10;
    ptube[2] = 0.55;
    
    gMC->Gsvolu("QB77","TUBE", idtmed[kC], ptube, 3);

    ptube[0] = 3.15;
    ptube[1] = 3.50;
    ptube[2] = 0.10;
    gMC->Gsvolu("QBB7","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QBB7", 1, "QB77", 0.0, 0.0, 0.55-0.2, 0, "ONLY");
    gMC->Gspos("QB77", 1, "QBT2", 0.0, 0.0, 2., 0, "ONLY");
 */
//
// Beam pipes between elements
//

    gMC->Gsvolu("QB24","TUBE", idtmed[kInox], ptube, 0);
    ptube[0] = RinSt;
    ptube[1] = RoutSt;
    ptube[2] = hlenQb24[0];
    dz = hlenQbbe2 + ptube[2];
    gMC->Gsposp("QB24", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY", ptube, 3);
//
// Bellow on absorber side
    dz = dz+hlenQb24[0] + dzb;
    gMC->Gspos("QBE0", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY");
//
    ptube[2] = hlenQb24[1];
    dz = dz + dzb + ptube[2];
    gMC->Gsposp("QB24", 2 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY", ptube, 3);
    dz = dz + ptube[2];
    
//
// Flange
// 
//  Mother Volume
    ptube[0] = RinSt;
    ptube[1] = 4.300;
    ptube[2] = 1.4;
    
    gMC->Gsvolu("QFA0","TUBE", idtmed[kInox], ptube, 3);
    dz = dz + ptube[2];
    gMC->Gspos("QFA0", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY");
    dz = dz + ptube[2];
//
//
    ptube[0] = RinSt;
    ptube[1] = RoutSt;
    ptube[2] = hlenQb24[2];
    dz = dz + ptube[2];
    gMC->Gsposp("QB24", 3 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY", ptube, 3);


// --- Place the PIPE ghost volume (QBPM) in its mother volume (ALIC)
//    by rotating it to 180 deg. and make it invisible
// 
    gMC->Gspos("QBPM",1,"ALIC",0,0,0,idrotm[2013], "ONLY");
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
    ptube[1] =  5.4;
    ptube[2] = 14.6;
    gMC->Gsvolu("QI34","TUBE", idtmed[kInox], ptube, 3);
    
    ptube[0] =  0.0;
    ptube[1] =  4.8;
    ptube[2] = 14.6;
    gMC->Gsvolu("QI44","TUBE", idtmed[kAir], ptube, 3);
    gMC->Gspos("QI44", 1, "QI34", 0.0, 0.0, 0.0, 0, "ONLY");
// <-

  //
  // Flange ->
    ptube[0] =  5.41;
    ptube[1] =  7.30;
    ptube[2] =  2.15;
    gMC->Gsvolu("QI35","TUBE", idtmed[kInox], ptube, 3);
// <-
    gMC->Gspos("QI32", 1, "QBPM", 0.0, -44.25, zPump, 0, "ONLY");
    gMC->Gspos("QI33", 1, "QBPM", 0.0, -35.00, zPump,idrotm[2002], "ONLY");
    gMC->Gspos("QI34", 1, "QBPM", 0.0, -17.90, zPump,idrotm[2002], "ONLY");
    gMC->Gspos("QI35", 1, "QBPM", 0.0, -24.35, zPump,idrotm[2002], "ONLY");

    gMC->Gsatt("QBPM", "SEEN", 1);
    gMC->Gsatt("QBEM", "SEEN", 1);
}



//___________________________________________
void AliPIPEv0::CreateMaterials()
{
  //
  // Define materials for beam pipe
  //

  if(fDebug) printf("%s: Create PIPEv0 materials \n",ClassName());
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  // Steel (Inox)  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  // AlBe - alloy 
  Float_t aAlBe[2] = { 26.98, 9.01};
  Float_t zAlBe[2] = { 13.00, 4.00};
  Float_t wAlBe[2] = { 0.4, 0.6};

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
  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  //
  //     Vacuum 
  AliMaterial(16, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  //
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  //
  //     reduced density steel to approximate pump getter material
  AliMixture(20, "GETTER$", asteel, zsteel, 1.00, 4, wsteel);
  //     Al-Be alloy
  //     
  AliMixture(21, "AlBe$", aAlBe, zAlBe, 2.07, 2, wAlBe);
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
  AliMedium(15, "AIR",    15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
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

}










