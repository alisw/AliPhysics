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
 
ClassImp(AliPIPEv0)
 
//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0()
{
// Constructor
}

//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0(const char *name, const char *title)
  : AliPIPE(name,title)
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

    printf("Create PIPEvTemp geometry \n");
  

    Int_t *idtmed = fIdtmed->GetArray();
    Float_t ppcon[48], ptube[3], pbox[3];
    Int_t i=0;
    
    enum {kC=6, kAlu=9, kInox=19, kGetter=20, kBe=5, kVac=16, kAir=15, kAlBe=21};
    
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
// The peam pipe up to the Front Absorber
//
// Mother Volume QBPM
    const Float_t dbe1 = 15.;
    const Float_t dbe2 = 15.;
    ppcon[0]  =    0;
    ppcon[1]  =  360;
    ppcon[2]  =   15;
//  1: 
    ppcon[3]  = - 90;
    ppcon[4]  =    0;
    ppcon[5]  =    5.8;
//  2
    ppcon[6]  = - 81.0;
    ppcon[7]  =    0.;
    ppcon[8]  =    5.8;
//  3
    ppcon[9]  = - 81.;
    ppcon[10]  =    0.;
    ppcon[11] =    4.22;
//  4
    ppcon[12] = - 28.00-dbe2;
    ppcon[13] =    0;
    ppcon[14] =    4.22;
//  5
    ppcon[15] = - 28.00-dbe2;
    ppcon[16] =    0;
    ppcon[17] =    3.2;
//  6
    ppcon[18] =    0;
    ppcon[19] =    0;
    ppcon[20] =    3.2;
//  7
    ppcon[21] =    28.+dbe1;
    ppcon[22] =    0;
    ppcon[23] =    3.2;
//  8
    ppcon[24] =   28.+dbe1;
    ppcon[25] =    0;
    ppcon[26] =    4.22;
//  9
    ppcon[27] =  250;
    ppcon[28] =    0;
    ppcon[29] =   4.22;
// 10
    ppcon[30] =  250;
    ppcon[31] =    0;
    ppcon[32] =    5;
// 11
    ppcon[33] =  365;
    ppcon[34] =    0;
    ppcon[35] =    5;
// 12
    ppcon[36] =  365;
    ppcon[37] =    0;
    ppcon[38] =   56;
// 13
    ppcon[39] =  405;
    ppcon[40] =    0;
    ppcon[41] =   56;
// 14
    ppcon[42] =  405;
    ppcon[43] =    0;
    ppcon[44] =   5;
// 15
    ppcon[45] =  800;
    ppcon[46] =    0;
    ppcon[47] =    5;
    
    gMC->Gsvolu("QBPM", "PCON", idtmed[kAir], ppcon, 48);

//
// The Vacuum 
    ptube[0] =   0.0;
    ptube[1] =   2.9;
    ptube[2] = 445.0;
    
    gMC->Gsvolu("QBVA","TUBE", idtmed[kVac], ptube, 3);
    gMC->Gspos ("QBVA", 1, "QBPM", 0., 0., 355., 0, "ONLY");
//
// volume definitions of various sections
//
//  distance between bellows
    const Float_t dzbb =  8.0;
//  size of bellow
    const Float_t dzb  = 11.4;

    //half-lengths of various beam pipe sections
    Float_t hlenQbbe=43.3;
    Float_t hlenQbt1=5.5/2.;
    Float_t hlenQbab=286.6/2.;
    Float_t hlenQb10=286.6/2.;
    Float_t hlenQb29=4.4/2.;
    Float_t hlenQbe0=(2.* dzb + dzbb)/2.;
    Float_t hlenQb26=20./2.;
    Float_t hlenQb28=400.6/2.;

// Be Pipe in central Alice : length 43.3 * 2
    ptube[0] =  2.90;
    ptube[1] =  3.00;
    ptube[2] = hlenQbbe;
    
    gMC->Gsvolu("QBBE","TUBE", idtmed[kBe], ptube, 3);
    
//
// Metal-Metal Transitions : non-absorber side
//
//  Be-Al (or Be to Be-Al alloy)
//  Mother Volume
    ptube[0] = 2.900;
    ptube[1] = 4.200;
    ptube[2] = hlenQbt1;
    gMC->Gsvolu("QBT1","TUBE", idtmed[kAir], ptube, 3);

    ptube[0] = 2.900;
    ptube[1] = 3.150;
    ptube[2] = 0.375;
    //
    //  Be-part
    gMC->Gsvolu("QB01","TUBE", idtmed[kAlu], ptube, 3);

    ptube[1] = 3.000;
    gMC->Gsvolu("QBA1","TUBE", idtmed[kBe], ptube, 3);

    gMC->Gspos("QBA1", 1, "QB01", 0., 0., 0, 0, "ONLY");
    gMC->Gspos("QB01", 1, "QBT1", 0., 0.,-2.75+ptube[2], 0, "ONLY");

    //  Al-part
    //
    ptube[0] = 2.900;
    ptube[1] = 3.150;
    ptube[2] = 2.375;

    gMC->Gsvolu("QB03","TUBE", idtmed[kAlu], ptube, 3);
    gMC->Gspos("QB03", 1, "QBT1", 0., 0.,-2.+ptube[2], 0, "ONLY");
    
    
    ptube[0] = 3.15;
    ptube[1] = 3.50;
    ptube[2] = 0.10;

    gMC->Gsvolu("QB05","TUBE", idtmed[kAlu], ptube, 3);
    gMC->Gspos("QB05", 1, "QBT1", 0., 0., 2.55+ptube[2], 0, "ONLY");
    
    
    // Fixations
    ptube[0] = 0.0;
    ptube[1] = 0.1;
    ptube[2] = 0.5;
    
    gMC->Gsvolu("QB08","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB08", 1 ,"QBT1",  0.000,  3.650, -1.25, idrotm[2002], "ONLY");
    gMC->Gspos("QB08", 2 ,"QBT1",  3.161, -1.825, -1.25, idrotm[2001], "ONLY");
    gMC->Gspos("QB08", 3 ,"QBT1", -3.161, -1.825, -1.25, idrotm[2003], "ONLY");
    
    // Carbon ring
    ptube[0] = 3.15;
    ptube[1] = 4.10;
    ptube[2] = 0.55;
    
    gMC->Gsvolu("QB07","TUBE", idtmed[kC], ptube, 3);

    ptube[0] = 3.15;
    ptube[1] = 3.50;
    ptube[2] = 0.10;
    gMC->Gsvolu("QBA7","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QBA7", 1, "QB07", 0.0, 0.0, 0.55-0.2, 0, "ONLY");
    gMC->Gspos("QB07", 1, "QBT1", 0.0, 0.0, 2., 0, "ONLY");


//
// 1st section Alu non-absorber side
    ptube[0] = 2.9;
    ptube[1] = 3.0;
    ptube[2] = hlenQb10;
    
    gMC->Gsvolu("QB10","TUBE", idtmed[kAlu], ptube, 3);
//
// Support rollers: non absorber side
//
//  Mother volume
    ptube[0] = 3.2;
    ptube[1] = 4.8;
    ptube[2] = 3.0;
    gMC->Gsvolu("QBRM","TUBE", idtmed[kAir], ptube, 3);

    ptube[0] = 0.0;
    ptube[1] = 0.7;
    ptube[2] = 3.0;
    
    gMC->Gsvolu("QB30","TUBE", idtmed[kInox], ptube, 3);
    
    for (i=0; i<8; i++) {
	Float_t phi = 45.+i*45.*kDegrad;
	Float_t xpos = 4.*TMath::Sin(phi);
	Float_t ypos = 4.*TMath::Cos(phi);
	gMC->Gspos("QB30", i+1, "QBRM", xpos, ypos, 0, idrotm[2004+i], "ONLY");
    }

//
// Flanges: non absorber side
    ptube[0] = 3.0;
    ptube[1] = 4.9;
    ptube[2] = hlenQb29;
    
    gMC->Gsvolu("QB29","TUBE", idtmed[kInox], ptube, 3);

//
// Inox beam pipe: final section on non-absorber side

    ptube[0] =   2.90;
    ptube[1] =   2.98;
    ptube[2] =   hlenQb28;    

    gMC->Gsvolu("QB28","TUBE", idtmed[kInox], ptube, 3);

//
//  Undulated beam pipe
// 
/*
    Float_t pitch=0.25;
    Float_t thick=0.015;
    Float_t zundul=171;
    Float_t rundul=3.0;
    char cn48[][5]={"QN21","QN22","QN23","QN24","QN25","QN26","QN27","QN28"};

    Undulation("QUND",pitch,thick,zundul,rundul,cn48);
    gMC->Gspos("QUND", 1, "QBPM", 0., 0., 335.+zundul, 0, "ONLY");
*/

//  Al-Be (40-60 wgt%, rho=2.7 g/cm**3) beam pipe
//
    ptube[0] =   2.90;
    ptube[1] =   3.05;
    ptube[2] =   hlenQbab;    

    gMC->Gsvolu("QBAB","TUBE", idtmed[kAlBe], ptube, 3);

    
//
//  last piece : inox pipe from pump till end 
//
    //ptube[0] =   2.90;
    //ptube[1] =   2.98;
    //ptube[2] =  61.55;    

    //gMC->Gsvolu("QB48","TUBE", idtmed[kInox], ptube, 3);

/*
    ptube[0] = 2.90;
    ptube[1] = 2.98;
    ptube[2] = 1.0;
    
    gMC->Gsvolu("QB27","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB27", 1, "QBPM", 0.0, 0.0, 208.1, 0, "ONLY");
*/
//
// Aluminium collar for support rollers
    
    gMC->Gsvolu("QB25","TUBE", idtmed[kAlu], ptube, 0);


// 2.5 mm thick SS tube for hanging pump
    ptube[0] = 2.90;
    ptube[1] = 3.15;
    ptube[2] = hlenQb26;
    
    gMC->Gsvolu("QB26","TUBE", idtmed[kInox], ptube, 3);

//
// Bellows
//
// Mother Volume
    ptube[0] =  2.90;
    ptube[1] =  3.75;
    ptube[2] = (2.*dzb+dzbb)/2.;
    gMC->Gsvolu("QBE0","TUBE", idtmed[kAir], ptube, 3);

    ptube[2] = dzb/2.;

    gMC->Gsvolu("QBEM","TUBE", idtmed[kAir], ptube, 3);
    Float_t dz = (dzb+dzbb)/2.;
    gMC->Gspos("QBEM", 2 ,"QBE0", 0.0, 0.0, -dz, 0 , "ONLY");
    gMC->Gspos("QBEM", 1 ,"QBE0", 0.0, 0.0,  dz, idrotm[2012], "ONLY");
    
    ptube[0] = 2.90;
    ptube[1] = 3.25;
    ptube[2] = 3.70;
    
    gMC->Gsvolu("QB19","TUBE", idtmed[kVac], ptube, 3);
    gMC->Gspos("QB19", 1 ,"QBEM", 0.0, 0.0, 0.5, 0 , "ONLY");
    
    ptube[0] = 3.25;
    ptube[1] = 3.74;
    ptube[2] = 0.095;
    
    gMC->Gsvolu("QB18","TUBE", idtmed[kVac], ptube, 3);
    for (i=0; i<15; i++) {
	gMC->Gspos("QB18", i+1, "QBEM", 0.0, 0.0, 3.3-i*0.4, 0, "ONLY");
    }
    
    ptube[0] = 2.90;
    ptube[1] = 3.00;
    ptube[2] = 1.20;
    
    gMC->Gsvolu("QB21","TUBE", idtmed[kVac], ptube, 3);
    gMC->Gspos("QB21", 1 ,"QBEM", 0.0, 0.0, -4.5, 0 , "ONLY");
    
    ptube[0] = 3.250;
    ptube[1] = 3.750;
    ptube[2] = 0.005;
    
    gMC->Gsvolu("QB15","TUBE", idtmed[kInox], ptube, 3);
    for (i=0; i<30; i++) {
	gMC->Gspos("QB15", i+1, "QBEM", 0.0, 0.0, 3.4-i*0.2, 0, "ONLY");
    }
    
    ptube[0] = 3.740;
    ptube[1] = 3.750;
    ptube[2] = 0.095;
    
    gMC->Gsvolu("QB16","TUBE", idtmed[kInox], ptube, 3);
    for (i=0; i<15; i++) {
	gMC->Gspos("QB16", i+1, "QBEM", 0.0, 0.0, 3.3-i*0.4, 0, "ONLY");
    }
    
    ptube[0] = 3.250;
    ptube[1] = 3.260;
    ptube[2] = 0.095;
    
    gMC->Gsvolu("QB17","TUBE", idtmed[kInox], ptube, 3);
    for (i=0; i<14; i++) {
	gMC->Gspos("QB17", i+1, "QBEM", 0.0, 0.0, 3.1-i*0.4, 0, "ONLY");
    }
    
    ptube[0] = 3.250;
    ptube[1] = 3.260;
    ptube[2] = 0.3975;

    gMC->Gsvolu("QB14","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB14", 2 ,"QBEM", 0.0, 0.0, -2.8025, 0 , "ONLY");
    gMC->Gspos("QB14", 1 ,"QBEM", 0.0, 0.0,  3.8025, 0 , "ONLY");
    
    ptube[0] = 2.900;
    ptube[1] = 3.260;
    ptube[2] = 0.050;
    
    gMC->Gsvolu("QB13","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB13", 2 ,"QBEM", 0.0, 0.0, -3.25, 0 , "ONLY");
    gMC->Gspos("QB13", 1 ,"QBEM", 0.0, 0.0,  4.25, 0 , "ONLY");
    
    ptube[0] = 2.900;
    ptube[1] = 3.000;
    ptube[2] = 0.700;
    
    gMC->Gsvolu("QB12","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB12", 1 ,"QBEM", 0.0, 0.0, 5.0, 0, "ONLY");


//
//  pipe between Bellows
    ptube[0] = 2.9;
    ptube[1] = 3.0;
    ptube[2] = dzbb/2.;
    gMC->Gsvolu("QB23","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB23", 1 ,"QBE0", 0.0, 0.0, 0.0, 0, "ONLY");
    
//
// End Bellow

// **** placement of various sections on non-absorber side ****
    //
    // first the beryllium section : 43.3 cm long on both sides of IP
    // zstart=0, zend=43.3, length=43.3 on one side
    Float_t zpos=0;
    gMC->Gspos("QBBE", 1, "QBPM", 0., 0., zpos, 0, "ONLY");

    // next meta-metal transition QBT1 on on-absorber side
    // zstart=43.3, zend=48.8, length=5.5
    zpos = zpos + hlenQbbe + hlenQbt1;
    gMC->Gspos("QBT1", 1, "QBPM", 0., 0.,  zpos, 0, "ONLY");
    
    // Aluminium OR Al-be alloy section
    //zpos=zpos + hlenQbt1 + hlenQb10;
    //gMC->Gspos("QB10", 1, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");
    // zstart= 48.8, zend=335.4, length=286.6
    zpos = zpos + hlenQbt1 + hlenQbab;
    gMC->Gspos("QBAB", 1, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");

    // inox flange at the start of bellow
    // zstart=335.4, zend=339.8, length=4.4
    zpos = zpos + hlenQbab + hlenQb29;
    gMC->Gspos("QB29", 1, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");

    // bellow section
    // zstart=339.8, zend=370.6, length=30.9
    zpos = zpos + hlenQb29 + hlenQbe0;
    gMC->Gspos("QBE0", 2 ,"QBPM", 0.0, 0.0, zpos, 0, "ONLY");

    // inox flange at the end of bellow and start of thick inox for pump
    // zstart=370.6, zend=375.0, length=4.4
    zpos = zpos + hlenQbe0 + hlenQb29;
    gMC->Gspos("QB29", 2, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");

    // 2.5mm thick inox section to hang the pump
    // zstart=375.0, zend=395.0, length=20.0
    zpos = zpos + hlenQb29 + hlenQb26;
    gMC->Gspos("QB26", 1, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");

    //inox flange at the end of pump section to join to next pipe
    // zstart=395.0, zend=399.4, length=4.4
    zpos = zpos + hlenQb26 + hlenQb29;
    gMC->Gspos("QB29", 3, "QBPM", 0.0, 0.0, zpos, 0, "ONLY");

    //last inox section till 800 cm
    // zstart=399.4, zend=800., length=400.6
    zpos = zpos + hlenQb29 + hlenQb28;
    gMC->Gspos("QB28", 1, "QBPM", 0.0, 0.0, zpos, 0, "ONLY"); 


    // aluminium collars and support rollers at two places

    ptube[0] = 3.0;
    ptube[1] = 3.15;
    ptube[2] = 2.75;

    gMC->Gsposp("QB25", 1, "QBPM", 0., 0., 654.8, 0, "ONLY", ptube, 3);
    gMC->Gspos("QBRM", 1, "QBPM", 0., 0., 654.8, 0, "ONLY");

    ptube[0] = 3.05;
    ptube[1] = 3.20;
    gMC->Gsposp("QB25", 2, "QBPM", 0., 0., 254.8, 0, "ONLY", ptube, 3);
    gMC->Gspos("QBRM", 2, "QBPM", 0., 0., 254.8, 0, "ONLY");



//******** end of placement on no-absorber side *********

    //
    // **** Absorber side *****   
    //
    //
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


// beam pipe between metal-metal transition and bellows
    ptube[0] = 2.9;
    ptube[1] = 3.0;
    ptube[2] = (81.7-(2.*dzb+dzbb)-(hlenQbbe+5.5))/2.;
    

    gMC->Gsvolu("QB24","TUBE", idtmed[kInox], ptube, 3);
//  z = 48.8 - 50.9    
    dz = (hlenQbbe+5.5)+ptube[2];
    gMC->Gspos("QB24", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY");
//
// beam pipe between flange and bellows    
    ptube[0] = 2.90;
    ptube[1] = 3.00;
    ptube[2] = 0.45;

    gMC->Gsvolu("QB22","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB22", 1 ,"QBPM", 0.0, 0.0, -82.15, 0, "ONLY");

// 
// Flange
// 
//  Mother Volume
    ptube[0] = 2.900;
    ptube[1] = 4.300;
    ptube[2] = 1.400;
    
    gMC->Gsvolu("QFA0","TUBE", idtmed[kAlu], ptube, 3);
    gMC->Gspos("QFA0", 1 ,"QBPM", 0.0, 0.0, -84.0, 0, "ONLY");
//
//  inner Inox piece
    ptube[0] = 2.900;
    ptube[1] = 3.500;
    ptube[2] = 0.450;
    gMC->Gsvolu("QFA1","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QFA1", 1 ,"QFA0", 0.0, 0.0, 0.225, 0, "ONLY");
//
//  8 x M5 Inox
    ptube[0] = 0.000;
    ptube[1] = 0.250;
    ptube[2] = 1.400;
    gMC->Gsvolu("QFA2","TUBE", idtmed[kInox], ptube, 3);
    for (i=0; i<8; i++) {
	Float_t phi = i*45.*kDegrad;
	Float_t xpos = 3.9*TMath::Sin(phi);
	Float_t ypos = 3.9*TMath::Cos(phi);
	gMC->Gspos("QFA2", i+1, "QFA0", xpos, ypos, 0., 0, "ONLY");
    }


    ptube[0] = 2.900;
    ptube[1] = 3.000;
    ptube[2] = 2.300;
    
    gMC->Gsvolu("QB32","TUBE", idtmed[kInox], ptube, 3);
    gMC->Gspos("QB32", 1 ,"QBPM", 0.0, 0.0, -90.+2.3, 0, "ONLY");

    //bellow on absorber side
    ptube[2] = (2.*dzb+dzbb)/2.;
    dz = (81.7-ptube[2]);
    gMC->Gspos("QBE0", 1 ,"QBPM", 0.0, 0.0, -dz, 0, "ONLY");


// --- Place the PIPE ghost volume (QBPM) in its mother volume (ALIC)
//    by rotating it to 180 deg. and make it invisible
// 

  
    gMC->Gspos("QBPM",1,"ALIC",0,0,0,idrotm[2013], "ONLY");


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
    //ptube[2] = 13.7;
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
    Float_t zPump = 385.;
    
    gMC->Gspos("QI32", 1, "QBPM", 0.0, -44.25, zPump, 0, "ONLY");
    gMC->Gspos("QI33", 1, "QBPM", 0.0, -35.00, zPump,idrotm[2002], "ONLY");
    gMC->Gspos("QI34", 1, "QBPM", 0.0, -17.90, zPump,idrotm[2002], "ONLY");
    gMC->Gspos("QI35", 1, "QBPM", 0.0, -24.35, zPump,idrotm[2002], "ONLY");

    gMC->Gsatt("QBPM", "SEEN", 0);
    gMC->Gsatt("QBEM", "SEEN", 0);
}



//___________________________________________
void AliPIPEv0::CreateMaterials()
{
  //
  // Define materials for beam pipe
  //

  printf("Create PIPEvTemp materials \n");
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










