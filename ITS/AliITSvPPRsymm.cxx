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
Revision 1.12  2000/12/10 16:00:44  barbera
Added last definition of special media like end-ladder boxes and cones

Revision 1.11  2000/10/30 08:02:25  barbera
PCON's changed into simpler CONS and TUBS. Services now allow for the rails to go through them.

Revision 1.3.2.7  2000/10/27 17:20:00  barbera
Position of rails w.r.t. the interaction point corrected.

Revision 1.9  2000/10/27 13:31:29  barbera
Rails between ITS and TPC added.

Revision 1.8  2000/10/27 13:03:08  barbera
Small changes in the SPD volumes and materials

Revision 1.6  2000/10/16 14:45:37  barbera
Mother volume ITSD modified to avoid some overlaps

Revision 1.5  2000/10/16 13:49:15  barbera
Services volumes slightly modified and material added following Pierluigi Barberis' information

Revision 1.4  2000/10/07 15:33:07  barbera
Small corrections to the ITSV mother volume

Revision 1.3  2000/10/07 13:06:50  barbera
Some new materials and media defined

Revision 1.2  2000/10/07 10:42:43  barbera
Mother volume ITSV corrected

Revision 1.1  2000/10/06 23:09:12  barbera
New  geometry (symmetric services

Revision 1.20  2000/10/02 21:28:08  fca
Removal of useless dependecies via forward declarations

Revision 1.19  2000/07/10 16:07:19  fca
Release version of ITS code

Revision 1.14.2.2  2000/05/19 10:09:21  nilsen
fix for bug with HP and Sun unix + fix for event display in ITS-working branch

Revision 1.14.2.1  2000/03/04 23:45:19  nilsen
Fixed up the comments/documentation.

Revision 1.14  1999/11/25 06:52:56  fca
Correct value of drca

Revision 1.13.2.1  1999/11/25 06:52:21  fca
Correct value of drca

Revision 1.13  1999/10/27 11:16:26  fca
Correction of problem in geometry

Revision 1.12  1999/10/22 08:25:25  fca
remove double definition of destructors

Revision 1.11  1999/10/22 08:16:49  fca
Correct destructors, thanks to I.Hrivnacova

Revision 1.10  1999/10/06 19:56:50  fca
Add destructor

Revision 1.9  1999/10/05 08:05:09  fca
Minor corrections for uninitialised variables.

Revision 1.8  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System version PPR  symmetric                                         //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
// Authors: R. Barbera
// version 6.
// Created  2000.
//
//  NOTE: THIS IS THE  SYMMETRIC PPR geometry of the ITS. 
// THIS WILL NOT WORK
// with the geometry or module classes or any analysis classes. You are 
// strongly encouraged to uses AliITSv5.
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TFile.h>    // only required for Tracking function?
#include <TCanvas.h>
#include <TObjArray.h>
#include <TClonesArray.h>


#include "AliMC.h"
#include "AliMagF.h"
#include "AliConst.h"

#include "AliITShit.h"
#include "AliITSvPPRsymm.h"
#include "AliRun.h"


ClassImp(AliITSvPPRsymm)
 
//_____________________________________________________________________________
AliITSvPPRsymm::AliITSvPPRsymm() {
////////////////////////////////////////////////////////////////////////
//    Standard default constructor for the ITS version 9.
////////////////////////////////////////////////////////////////////////

  fIdN = 6;
  fIdName = new TString[fIdN];
  fIdName[0] = "ITS1";
  fIdName[1] = "ITS2";
  fIdName[2] = "ITS3";
  fIdName[3] = "ITS4";
  fIdName[4] = "ITS5";
  fIdName[5] = "ITS6";
  fIdSens    = new Int_t[fIdN];
  for (Int_t i=0;i<fIdN;i++) fIdSens[i]=fIdName[i].Length();
}
//_____________________________________________________________________________
AliITSvPPRsymm::AliITSvPPRsymm(const char *name, const char *title) : AliITS(name, title){
////////////////////////////////////////////////////////////////////////
//    Standard constructor for the ITS version 9.
////////////////////////////////////////////////////////////////////////

  fIdN = 6;
  fIdName = new TString[fIdN];
  fIdName[0] = "ITS1";
  fIdName[1] = "ITS2";
  fIdName[2] = "ITS3";
  fIdName[3] = "ITS4";
  fIdName[4] = "ITS5";
  fIdName[5] = "ITS6";
  fIdSens    = new Int_t[fIdN];
  for (Int_t i=0;i<fIdN;i++) fIdSens[i]=fIdName[i].Length();

}
//____________________________________________________________________________
AliITSvPPRsymm::AliITSvPPRsymm(const AliITSvPPRsymm &source){
////////////////////////////////////////////////////////////////////////
//     Copy Constructor for ITS version 9.
////////////////////////////////////////////////////////////////////////
    if(&source == this) return;
    printf("Not allowed to copy AliITSvPPRsymm\n");
    return;
}
//_____________________________________________________________________________
AliITSvPPRsymm& AliITSvPPRsymm::operator=(const AliITSvPPRsymm &source){
////////////////////////////////////////////////////////////////////////
//    Assignment operator for the ITS version 7.
////////////////////////////////////////////////////////////////////////
  if(&source == this) return *this;
    printf("Not allowed to copy AliITSvPPRsymm\n");
  return *this;
}
//_____________________________________________________________________________
AliITSvPPRsymm::~AliITSvPPRsymm() {
////////////////////////////////////////////////////////////////////////
//    Standard destructor for the ITS version 7.
////////////////////////////////////////////////////////////////////////
}

//__________________________________________________________________________
void AliITSvPPRsymm::BuildGeometry(){
////////////////////////////////////////////////////////////////////////
//    Geometry builder for the ITS version 9.
////////////////////////////////////////////////////////////////////////
    TNode *node, *top;
    const int kColorITS=kYellow;
    //
    top = gAlice->GetGeometry()->GetNode("alice");

    new TTUBE("S_layer1","Layer1 of ITS","void",3.95,3.95+0.05475,12.25);
    top->cd();
    node = new TNode("Layer1","Layer1","S_layer1",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer2","Layer2 of ITS","void",7.,7.+0.05475,16.3);
    top->cd();
    node = new TNode("Layer2","Layer2","S_layer2",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer3","Layer3 of ITS","void",15.,15.+0.05288,21.1);
    top->cd();
    node = new TNode("Layer3","Layer3","S_layer3",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer4","Layer4 of ITS","void",24,24+0.05288,29.6);
    top->cd();
    node = new TNode("Layer4","Layer4","S_layer4",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer5","Layer5 of ITS","void",40,40+0.05382,45.1);
    top->cd();
    node = new TNode("Layer5","Layer5","S_layer5",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer6","Layer6 of ITS","void",45,45+0.05382,50.4);
    top->cd();
    node = new TNode("Layer6","Layer6","S_layer6",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);
}
//_____________________________________________________________________________
void AliITSvPPRsymm::CreateGeometry(){
////////////////////////////////////////////////////////////////////////
//    This routine defines and Creates the geometry for version 9 of the ITS.
////////////////////////////////////////////////////////////////////////
  
  //INNER RADII OF THE SILICON LAYERS 
  Float_t rl[6]    = { 3.8095,7.,15.,24.,38.1,43.5765 };   
  //THICKNESSES OF LAYERS (in % radiation length)
  Float_t drl[6]   = { 1.03,1.03,0.94,0.95,0.91,0.87 };   
  //HALF LENGTHS OF LAYERS  
  Float_t dzl[6]   = { 14.35,14.35,25.1,32.1,49.405,55.27 };
  //LENGTHS OF END-LADDER BOXES (ALL INCLUDED)
  Float_t dzb[6]   = { 12.4,12.4,13.5,15.,7.5,7.5 };    
  //THICKNESSES OF END-LADDER BOXES (ALL INCLUDED)
  Float_t drb[6]   = { rl[1]-rl[0],0.2,5.,5.,4.,4. };         

 
  Float_t dits[50], rlim, zmax;
  Float_t zpos;
  Float_t pcits[50], xltpc;
  Int_t idrotm[999], i;
  Float_t dgh[50];
  
  Int_t *idtmed = fIdtmed->GetArray()-199;
  
  // Rotation matrices
  
  // SPD - version 'a'
  
  AliMatrix(idrotm[201],90.0,90.0,90.0,180.0,0.0,0.0);
  AliMatrix(idrotm[202],90.0,90.0,90.0,0.0,0.0,0.0);
  AliMatrix(idrotm[203],90.0,350.0,90.0,260.0,0.0,0.0);
  AliMatrix(idrotm[204],90.0,170.0,90.0,80.0,0.0,0.0);
  AliMatrix(idrotm[205],90.0,10.0,90.0,100.0,0.0,0.0);
  AliMatrix(idrotm[206],90.0,190.0,90.0,280.0,0.0,0.0);
  AliMatrix(idrotm[207],90.0,342.0,90.0,72.0,0.0,0.0);
  AliMatrix(idrotm[208],90.0,156.999893,90.0,246.999893,0.0,0.0);
  AliMatrix(idrotm[209],90.0,147.999802,90.0,237.999893,0.0,0.0);
  AliMatrix(idrotm[210],90.0,138.999802,90.0,228.999802,0.0,0.0);
  AliMatrix(idrotm[211],90.0,129.999802,90.0,219.999802,0.0,0.0);
  AliMatrix(idrotm[212],90.0,36.7896,90.0,126.789597,0.0,0.0);
  AliMatrix(idrotm[213],90.0,343.579712,90.0,73.579697,0.0,0.0);
  AliMatrix(idrotm[214],90.0,95.413696,90.0,185.413696,0.0,0.0);
  AliMatrix(idrotm[215],90.0,5.4141,90.0,95.414101,0.0,0.0);
  AliMatrix(idrotm[216],90.0,318.296906,90.0,48.296902,0.0,0.0);
  AliMatrix(idrotm[217],90.0,67.000099,90.0,157.000107,0.0,0.0);
  AliMatrix(idrotm[218],90.0,337.003998,90.0,67.003998,0.0,0.0);
  AliMatrix(idrotm[219],90.0,247.000305,90.0,337.000305,0.0,0.0);
  AliMatrix(idrotm[220],90.0,305.633514,90.0,35.633499,0.0,0.0);
  AliMatrix(idrotm[221],90.0,58.000198,90.0,148.000198,0.0,0.0);
  AliMatrix(idrotm[222],90.0,327.997101,90.0,57.997101,0.0,0.0 );
  AliMatrix(idrotm[223],90.0,237.994202,90.0,327.994202,0.0,0.0);
  AliMatrix(idrotm[224],90.0,296.627502,90.0,26.627399,0.0,0.0);
  AliMatrix(idrotm[225],90.0,48.994099,90.0,138.994095,0.0,0.0);
  AliMatrix(idrotm[226],90.0,318.990997,90.0,48.991001,0.0,0.0);
  AliMatrix(idrotm[227],90.0,228.988205,90.0,318.98819,0.0,0.0);
  AliMatrix(idrotm[228],90.0,287.621399,90.0,17.621401,0.0,0.0);
  AliMatrix(idrotm[229],90.0,39.988098,90.0,129.988098,0.0,0.0);
  AliMatrix(idrotm[230],90.0,309.984985,90.0,39.985001,0.0,0.0);
  AliMatrix(idrotm[231],90.0,327.2612,90.0,57.2612,0.0,0.0);
  AliMatrix(idrotm[232],90.0,237.261398,90.0,327.261414,0.0,0.0);
  AliMatrix(idrotm[233],90.0,252.000504,90.0,342.000488,0.0,0.0 );
  AliMatrix(idrotm[234],90.0,71.9991,90.0,161.9991,0.0,0.0);
  AliMatrix(idrotm[235],90.0,270.0,90.0,0.0,0.0,0.0);
  AliMatrix(idrotm[236],90.0,180.013702,90.0,270.013702,0.0,0.0);
  AliMatrix(idrotm[237],90.0,0.0,90.0,90.0,180.0,0.0);
  AliMatrix(idrotm[238],90.0,144.0,90.0,234.0,0.0,0.0);
  AliMatrix(idrotm[239],90.0,216.0,90.0,306.0,0.0,0.0);
  AliMatrix(idrotm[240],90.0,288.0,90.0,18.0,0.0,0.0);
  AliMatrix(idrotm[241],90.0,324.0,90.0,54.0,0.0,0.0);
  AliMatrix(idrotm[242],90.0,36.0,90.0,126.0,0.0,0.0);
  AliMatrix(idrotm[243],90.0,108.0,90.0,198.0,0.0,0.0);
  AliMatrix(idrotm[244],90.0,0.0,90.0,270.0,180.0,0.0);
  AliMatrix(idrotm[245],90.0,342.0,90.0,252.0,180.0,0.0);
  AliMatrix(idrotm[246],90.0,130.0,90.0,40.0,180.0,0.0);
  AliMatrix(idrotm[247],90.0,139.0,90.0,49.0,180.0,0.0);
  AliMatrix(idrotm[248],90.0,148.0,90.0,58.0,180.0,0.0);
  AliMatrix(idrotm[249],90.0,157.0,90.0,67.0,180.0,0.0);
  
  // SDD
  
  AliMatrix(idrotm[301],0.0,0.0,90.0,90.0,90.0,180.0);  
  AliMatrix(idrotm[302],0.0,0.0,90.0,90.0,90.0,0.0);
  AliMatrix(idrotm[303],180.0,0.0,90.0,90.0,90.0,0.0); 
  AliMatrix(idrotm[304],180.0,0.0,90.0,90.0,90.0,180.0); 
  AliMatrix(idrotm[305],90.0,347.14,90.0,77.14,0.0,0.0); 
  AliMatrix(idrotm[306],90.0,321.43,90.0,51.43,0.0,0.0); 
  AliMatrix(idrotm[307],90.0,295.71,90.0,25.71,0.0,0.0);
  AliMatrix(idrotm[308],90.0,244.29,90.0,334.29,0.0,0.0);
  AliMatrix(idrotm[309],90.0,218.57,90.0,308.57,0.0,0.0);
  AliMatrix(idrotm[310],90.0,167.14,90.0,257.14,0.0,0.0);
  AliMatrix(idrotm[311],90.0,141.43,90.0,231.43,0.0,0.0);  
  AliMatrix(idrotm[312],90.0,0.0,0.0,0.0,90.0,270.0);
  AliMatrix(idrotm[313],90.0,115.71,90.0,205.71,0.0,0.0); 
  AliMatrix(idrotm[314],90.0,335.45,90.0,65.45,0.0,0.0); 
  AliMatrix(idrotm[315],90.0,319.09,90.0,49.09,0.0,0.0); 
  AliMatrix(idrotm[316],90.0,302.73,90.0,32.73,0.0,0.0); 
  AliMatrix(idrotm[317],90.0,286.36,90.0,16.36,0.0,0.0);
  AliMatrix(idrotm[318],90.0,270.0,90.0,360.0,0.0,0.0);
  AliMatrix(idrotm[319],90.0,253.64,90.0,343.64,0.0,0.0);
  AliMatrix(idrotm[320],90.0,237.27,90.0,327.27,0.0,0.0);
  AliMatrix(idrotm[321],90.0,12.86,90.0,102.86,0.0,0.0);  
  AliMatrix(idrotm[322],90.0,220.91,90.0,310.91,0.0,0.0);
  AliMatrix(idrotm[323],90.0,204.55,90.0,294.55,0.0,0.0); 
  AliMatrix(idrotm[324],90.0,188.18,90.0,278.18,0.0,0.0); 
  AliMatrix(idrotm[325],90.0,171.82,90.0,261.82,0.0,0.0); 
  AliMatrix(idrotm[326],90.0,155.45,90.0,245.45,0.0,0.0); 
  AliMatrix(idrotm[327],90.0,139.09,90.0,229.09,0.0,0.0);
  AliMatrix(idrotm[328],90.0,122.73,90.0,212.73,0.0,0.0);
  AliMatrix(idrotm[329],90.0,106.36,90.0,196.36,0.0,0.0);
  AliMatrix(idrotm[330],90.0,73.64,90.0,163.64,0.0,0.0);    
  AliMatrix(idrotm[331],90.0,40.91,90.0,130.91,0.0,0.0);  
  AliMatrix(idrotm[332],90.0,24.55,90.0,114.55,0.0,0.0);
  AliMatrix(idrotm[333],90.0,38.57,90.0,128.57,0.0,0.0); 
  AliMatrix(idrotm[334],90.0,351.82,90.0,81.82,0.0,0.0); 
  AliMatrix(idrotm[335],90.0,8.18,90.0,98.18,0.0,0.0); 
  AliMatrix(idrotm[336],90.0,64.29,90.0,154.29,0.0,0.0); 
  AliMatrix(idrotm[337],111.0,300.0,21.0,300.0,90.0,30.0);
  AliMatrix(idrotm[338],69.0,240.0,159.0,240.0,90.0,150.0);
  AliMatrix(idrotm[339],111.0,240.0,21.0,240.0,90.0,150.0);
  AliMatrix(idrotm[340],69.0,300.0,159.0,300.0,90.0,30.0);  
  AliMatrix(idrotm[341],128.0,0.0,38.0,0.0,90.0,270.0);  
  AliMatrix(idrotm[342],90.0,240.0,180.0,0.0,90.0,330.);
  AliMatrix(idrotm[343],90.0,120.0,180.0,0.0,90.0,210.0); 
  AliMatrix(idrotm[344],90.0,0.0,180.0,0.0,90.0,90.0); 
  AliMatrix(idrotm[345],90.0,180.0,90.0,90.0,0.0,0.0); 
  AliMatrix(idrotm[346],90.0,300.0,90.0,30.0,0.0,0.0); 
  AliMatrix(idrotm[347],90.0,240.0,90.0,150.0,0.0,0.0);
  AliMatrix(idrotm[348],90.0,180.0,0.0,0.0,90.0,270.0);
  AliMatrix(idrotm[349],90.0,235.0,90.0,145.0,0.0,0.0);
  AliMatrix(idrotm[350],90.0,90.0,90.0,180.0,0.0,0.0);  
  AliMatrix(idrotm[351],90.0,305.0,90.0,35.0,0.0,0.0);  
  AliMatrix(idrotm[352],0.0,0.0,90.0,0.0,90.0,90.0);
  AliMatrix(idrotm[353],90.0,60.0,90.0,150.0,0.0,0.0); 
  AliMatrix(idrotm[354],90.0,120.0,90.0,30.0,0.0,0.0); 
  AliMatrix(idrotm[355],90.0,180.0,90.0,90.0,180.0,0.0); 
  AliMatrix(idrotm[356],90.0,270.0,90.0,0.0,0.0,0.0); 
  AliMatrix(idrotm[366],90.0,57.27,90.0,147.27,0.0,0.0); 
  AliMatrix(idrotm[386],90.0,192.86,90.0,282.86,0.0,0.0);  
   
  // SSD
  
  AliMatrix(idrotm[501],90.0,148.24,90.0,238.24,0.0,0.0);
  AliMatrix(idrotm[503],90.0,137.65,90.0,227.65,0.0,0.0); 
  AliMatrix(idrotm[504],90.0,127.06,90.0,217.06,0.0,0.0);  
  AliMatrix(idrotm[505],90.0,116.47,90.0,206.47,0.0,0.0);  
  AliMatrix(idrotm[506],90.0,105.88,90.0,195.88,0.0,0.0);  
  AliMatrix(idrotm[507],90.0,95.29,90.0,185.29,0.0,0.0);  
  AliMatrix(idrotm[508],90.0,84.71,90.0,174.71,0.0,0.0);
  AliMatrix(idrotm[509],90.0,74.12,90.0,164.12,0.0,0.0);
  AliMatrix(idrotm[510],90.0,63.53,90.0,153.53,0.0,0.0);  
  AliMatrix(idrotm[511],90.0,52.94,90.0,142.94,0.0,0.0);
  AliMatrix(idrotm[512],90.0,42.35,90.0,132.35,0.0,0.0);
  AliMatrix(idrotm[513],90.0,31.76,90.0,121.76,0.0,0.0); 
  AliMatrix(idrotm[514],90.0,10.59,90.0,100.59,0.0,0.0);  
  AliMatrix(idrotm[515],90.0,349.41,90.0,79.41,0.0,0.0);  
  AliMatrix(idrotm[516],90.0,338.82,90.0,68.82,0.0,0.0);  
  AliMatrix(idrotm[517],90.0,328.24,90.0,58.24,0.0,0.0);  
  AliMatrix(idrotm[518],90.0,317.65,90.0,47.65,0.0,0.0);
  AliMatrix(idrotm[519],90.0,307.06,90.0,37.06,0.0,0.0);
  AliMatrix(idrotm[520],90.0,296.47,90.0,26.47,0.0,0.0);  
  AliMatrix(idrotm[521],90.0,285.88,90.0,15.88,0.0,0.0);
  AliMatrix(idrotm[522],90.0,275.29,90.0,5.29,0.0,0.0);
  AliMatrix(idrotm[523],90.0,264.71,90.0,354.71,0.0,0.0); 
  AliMatrix(idrotm[524],90.0,254.12,90.0,344.12,0.0,0.0);  
  AliMatrix(idrotm[525],90.0,243.53,90.0,333.53,0.0,0.0);  
  AliMatrix(idrotm[526],90.0,232.94,90.0,322.94,0.0,0.0);  
  AliMatrix(idrotm[527],90.0,222.35,90.0,312.35,0.0,0.0);  
  AliMatrix(idrotm[528],90.0,211.76,90.0,301.76,0.0,0.0);
  AliMatrix(idrotm[529],90.0,190.59,90.0,280.59,0.0,0.0);
  AliMatrix(idrotm[530],90.0,169.41,90.0,259.41,0.0,0.0);  
  AliMatrix(idrotm[531],90.0,158.82,90.0,248.82,0.0,0.0);
  AliMatrix(idrotm[532],90.0,360.0,90.0,90.0,0.0,0.0);
  AliMatrix(idrotm[533],90.0,180.0,90.0,270.0,0.0,0.0); 
  AliMatrix(idrotm[534],90.0,189.47,90.0,279.47,0.0,0.0);  
  AliMatrix(idrotm[535],90.0,198.95,90.0,288.95,0.0,0.0 );  
  AliMatrix(idrotm[537],90.0,217.89,90.0,307.89,0.0,0.0);  
  AliMatrix(idrotm[538],90.0,227.37,90.0,317.37,0.0,0.0);
  AliMatrix(idrotm[539],90.0,236.84,90.0,326.84,0.0,0.0);
  AliMatrix(idrotm[540],90.0,246.32,90.0,336.32,0.0,0.0);  
  AliMatrix(idrotm[541],90.0,255.79,90.0,345.79,0.0,0.0);
  AliMatrix(idrotm[542],90.0,265.26,90.0,355.26,0.0,0.0);
  AliMatrix(idrotm[543],90.0,274.74,90.0,4.74,0.0,0.0); 
  AliMatrix(idrotm[544],90.0,284.21,90.0,14.21,0.0,0.0);  
  AliMatrix(idrotm[545],90.0,293.68,90.0,23.68,0.0,0.0);  
  AliMatrix(idrotm[546],90.0,303.16,90.0,33.16,0.0,0.0);  
  AliMatrix(idrotm[547],90.0,312.63,90.0,42.63,0.0,0.0);  
  AliMatrix(idrotm[548],90.0,322.11,90.0,52.11,0.0,0.0);
  AliMatrix(idrotm[549],90.0,331.58,90.0,61.58,0.0,0.0);
  AliMatrix(idrotm[550],90.0,341.05,90.0,71.05,0.0,0.0);  
  AliMatrix(idrotm[551],90.0,350.53,90.0,80.53,0.0,0.0);
  AliMatrix(idrotm[552],90.0,9.47,90.0,99.47,0.0,0.0);
  AliMatrix(idrotm[553],90.0,18.95,90.0,108.95,0.0,0.0 ); 
  AliMatrix(idrotm[555],90.0,37.89,90.0,127.89,0.0,0.0);  
  AliMatrix(idrotm[556],90.0,47.37,90.0,137.37,0.0,0.0);  
  AliMatrix(idrotm[557],90.0,56.84,90.0,146.84,0.0,0.0);  
  AliMatrix(idrotm[558],90.0,66.32,90.0,156.32,0.0,0.0);
  AliMatrix(idrotm[559],90.0,75.79,90.0,165.79,0.0,0.0);
  AliMatrix(idrotm[560],90.0,85.26,90.0,175.26,0.0,0.0);  
  AliMatrix(idrotm[561],90.0,94.74,90.0,184.74,0.0,0.0);
  AliMatrix(idrotm[562],90.0,104.21,90.0,194.21,0.0,0.0);
  AliMatrix(idrotm[563],90.0,113.68,90.0,203.68,0.0,0.0); 
  AliMatrix(idrotm[564],90.0,123.16,90.0,213.16,0.0,0.0);  
  AliMatrix(idrotm[565],90.0,132.63,90.0,222.63,0.0,0.0);  
  AliMatrix(idrotm[566],90.0,142.11,90.0,232.11,0.0,0.0);  
  AliMatrix(idrotm[567],90.0,151.58,90.0,241.58,0.0,0.0);  
  AliMatrix(idrotm[568],90.0,161.05,90.0,251.05,0.0,0.0);
  AliMatrix(idrotm[569],90.0,170.53,90.0,260.53,0.0,0.0);
  AliMatrix(idrotm[570],90.0,180.0,90.0,90.0,180.0,0.0);  
  AliMatrix(idrotm[571],90.0,0.0,0.0,0.0,90.0,270.0);
  AliMatrix(idrotm[572],90.0,180.0,0.0,0.0,90.0,270.0);
  AliMatrix(idrotm[573],90.0,180.0,90.0,90.0,0.0,0.0); 
  AliMatrix(idrotm[575],90.0,120.0,180.0,0.0,90.0,210.0);  
  AliMatrix(idrotm[576],65.71,300.0,90.0,30.0,24.29,120.0);  
  AliMatrix(idrotm[577],114.29,300.0,90.0,30.0,155.71,120.0);  
  AliMatrix(idrotm[579],65.71,240.0,90.0,150.0,24.29,60.0);
  AliMatrix(idrotm[580],114.29,240.0,90.0,150.0,155.71,60.0);  
  AliMatrix(idrotm[581],90.0,240.0,180.0,0.0,90.0,330.0);
  AliMatrix(idrotm[583],90.0,0.0,180.0,0.0,90.0,90.0); 
  AliMatrix(idrotm[584],90.0,180.0,180.0,0.0,90.0,90.0);  
  AliMatrix(idrotm[586],180.0,0.0,90.0,90.0,90.0,0.0);  
  AliMatrix(idrotm[618],90.0,201.18,90.0,291.18,0.0,0.0);
  AliMatrix(idrotm[620],90.0,28.42,90.0,118.42,0.0,0.0);  
  AliMatrix(idrotm[623],90.0,208.42,90.0,298.42,0.0,0.0);
  AliMatrix(idrotm[633],132.46,0.0,90.0,90.0,42.46,360.0);
  AliMatrix(idrotm[653],90.0,21.18,90.0,111.18,0.0,0.0); 



  //     CONVERT INTO CM (RL(SI)=9.36 CM) 
  for (i = 0; i < 6; ++i) {
    drl[i] = drl[i] / 100. * 9.36;
  }
    
  //     FIELD CAGE HALF LENGTH 
  
  rlim  = 56.;
  zmax  = 76.708;
  xltpc = 284.;
    
  // --- Define ghost volume containing the whole ITS (including services) 
  //     and fill it with air 
  
  dgh[0] = 0.;
  dgh[1] = 360.;
  dgh[2] = 12.;
  dgh[3] = -xltpc-5.-0.1;
  dgh[4] = 62.4;
  dgh[5] = 85.;
  dgh[6] = -xltpc;
  dgh[7] = 61.5;
  dgh[8] = 85.;
  dgh[9] = -xltpc;
  dgh[10] = 61.5;
  dgh[11] = 61.5+4.;
  dgh[12] = -100.7;
  dgh[13] = 44.9;
  dgh[14] = 56.1;
  dgh[15] = -77.2;
  dgh[16] = 44.9;
  dgh[17] = 56.1;
  dgh[18] = -40.;
  dgh[19] = 3.295;
  dgh[20] = 56.1; 

/*
  dgh[21] = -35.;
  dgh[22] = 3.295;
  dgh[23] = 56.1;

  dgh[24] = -35.;
  dgh[25] = 5.;
  dgh[26] = 56.1;

  dgh[27] = -29.;
  dgh[28] = 5.;
  dgh[29] = 56.1;
  
  dgh[30] = -29.;
  dgh[31] = 3.295;
  dgh[32] = 56.1;

*/


  dgh[21] = 40.;
  dgh[22] = 3.295;
  dgh[23] = 56.1;
  dgh[24] = 77.2;
  dgh[25] = 44.9;
  dgh[26] = 56.1;
  dgh[27] = 100.7;
  dgh[28] = 44.9;
  dgh[29] = 56.1;
  dgh[30] = xltpc;
  dgh[31] = 61.5;
  dgh[32] = 61.5+4.;
  dgh[33] = xltpc;
  dgh[34] = 61.5;
  dgh[35] = 85.;
  dgh[36] = xltpc+4.+0.1;
  dgh[37] = 62.4;
  dgh[38] = 85.;

  gMC->Gsvolu("ITSV", "PCON", idtmed[275], dgh, 39);
  
  // --- Place the ghost volume in its mother volume (ALIC) and make it 
  //     invisible 
  
  gMC->Gspos("ITSV", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  //gMC->Gsatt("ITSV", "SEEN", 0); 
 
  
  // --- Define ghost volume containing the six layers and fill it with air 
  
  dgh[0] = 0.;
  dgh[1] = 360.;
  dgh[2] = 4.;
  dgh[3] = -77.2;
  dgh[4] = 45.;
  dgh[5] = 56.;
  dgh[6] = -40.;     
  dgh[7] = 3.3;
  dgh[8] = 56.;
  dgh[9] = 40.;
  dgh[10] = 3.3;
  dgh[11] = 56.;
  dgh[12] = 77.2;
  dgh[13] = 45.;
  dgh[14] = 56.;
  gMC->Gsvolu("ITSD", "PCON", idtmed[275], dgh, 15);
  
  // --- Place the ghost volume in its mother volume (ALIC) and make it 
  //     invisible 
  
  gMC->Gspos("ITSD", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  //gMC->Gsatt("ITSD", "SEEN", 0);
  
  // --- Define SPD (version 'a') volumes ----------------------------
  
  dits[0] = 3.7;
  dits[1] = 7.75;
  dits[2] = 24;
  gMC->Gsvolu("IT12", "TUBE", idtmed[200], dits, 3);   

  dits[0] = 3.7;
  dits[1] = 7.7;
  dits[2] = 24;
  dits[3] = 57;
  dits[4] = 100;
  gMC->Gsvolu("I12A", "TUBS", idtmed[200], dits, 5); 
  
  dits[0] = 0.843;
  dits[1] = 0.025;
  dits[2] = 19.344;
  gMC->Gsvolu("I10A", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.843;
  dits[1] = 0.025;
  dits[2] = 19.344;
  gMC->Gsvolu("I20A", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 1.3673;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I123", "BOX ", idtmed[200], dits, 3);
  
  dits[0] = 0.06;
  dits[1] = 0.08;
  dits[2] = 24;
  dits[3] = -36.79;
  dits[4] = 21.834;
  gMC->Gsvolu("I121", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0.1253;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I122", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.04;
  dits[1] = 0.06 ;
  dits[2] = 24;
  dits[3] = 126.79;
  dits[4] = 270;
  gMC->Gsvolu("I120", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0.1134;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I144", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.25;
  dits[1] = 0.06;
  dits[2] = 24;
  gMC->Gsvolu("I113", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.077;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I143", "BOX ", idtmed[200], dits, 3);   

  dits[0] = 0.04;
  dits[1] = 0.06;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 90;
  gMC->Gsvolu("I142", "TUBS", idtmed[200], dits, 5); 
  
  dits[0] = 0.0695;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I141", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.06;
  dits[1] = 0.08;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 108;
  gMC->Gsvolu("I140", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0.1835;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I139", "BOX ", idtmed[200], dits, 3);
  
  dits[0] = 0.1894 ;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I138", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.04;
  dits[1] = 0.06;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 75.261;
  gMC->Gsvolu("I137", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 1.3401;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I136", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.05;
  dits[1] = 0.07;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 72.739;
  gMC->Gsvolu("I135", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0.1193;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I134", "BOX ", idtmed[200], dits, 3);    
  
  dits[0] = 0.163;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I133", "BOX ", idtmed[200], dits, 3);   

  dits[0] = 0.04;
  dits[1] = 0.06;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 157.633;
  gMC->Gsvolu("I132", "TUBS", idtmed[200], dits, 5); 
  
  dits[0] = 0.2497;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I131", "BOX ", idtmed[200], dits, 3); 
  
  dits[0] = 0.06;
  dits[1] = 0.08;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 148.633;
  gMC->Gsvolu("I130", "TUBS", idtmed[200], dits, 5); 

  dits[0] = 0.292;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I129", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.163;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I128", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.04;
  dits[1] = 0.06;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 161.297;
  gMC->Gsvolu("I126", "TUBS", idtmed[200], dits, 5);
  
  dits[0] = 0.2433;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I125", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.06;
  dits[1] = 0.08;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 42.883;
  gMC->Gsvolu("I124", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0.793;
  dits[1] = 0.0125;
  dits[2] = 3.536;
  gMC->Gsvolu("I103", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.793;
  dits[1] = 0.015 ;
  dits[2] = 2.5;
  gMC->Gsvolu("I105", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.843;
  dits[1] = 0.01;
  dits[2] = 19.344;
  gMC->Gsvolu("I104", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.793;
  dits[1] = 0.0125;
  dits[2] = 3.536;
  gMC->Gsvolu("I1D3", "BOX ", idtmed[200], dits, 3);
  
  dits[0] = 0.06;
  dits[1] = 0.08;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 80;
  gMC->Gsvolu("I112", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0.04;
  dits[1] = 0.06;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 80;
  gMC->Gsvolu("I111", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0.15;
  dits[1] = 0.0146;
  dits[2] = 24;
  gMC->Gsvolu("I118", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.1315;
  dits[1] = 0.01;
  dits[2] = 24;
  gMC->Gsvolu("I110", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.025;
  dits[1] = 0.035;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 180;
  gMC->Gsvolu("I114", "TUBS", idtmed[200], dits, 5);  
  
  dits[0] = 0;
  dits[1] = 0.025;
  dits[2] = 24;
  dits[3] = 0;
  dits[4] = 180;
  gMC->Gsvolu("I115", "TUBS", idtmed[200], dits, 5);   

  dits[0] = 0.063;
  dits[1] = 0.035;
  dits[2] = 24;
  gMC->Gsvolu("I116", "BOX ", idtmed[200], dits, 3); 
  
  dits[0] = 0.705;
  dits[1] = 0.005;
  dits[2] = 3.536;
  gMC->Gsvolu("I101", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.793;
  dits[1] = 0.0075;
  dits[2] = 0.68;
  gMC->Gsvolu("I102", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.705;
  dits[1] = 0.005;
  dits[2] = 3.536;
  gMC->Gsvolu("I1D1", "BOX ", idtmed[200], dits, 3);
  
  dits[0] = 0.063;
  dits[1] = 0.025;
  dits[2] = 24;
  gMC->Gsvolu("I117", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.64;
  dits[1] = 0.005;
  dits[2] = 3.48;
  gMC->Gsvolu("ITS1", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.64;
  dits[1] = 0.005;
  dits[2] = 3.48;
  gMC->Gsvolu("ITS2", "BOX ", idtmed[200], dits, 3);  


  // --- Define SDD volumes ------------------------------------------

  dits[0] = 0;
  dits[1] = 360;
  dits[2] = 6;
  dits[3] = -34.6;
  dits[4] = 23.495;
  dits[5] = 28.5;
  dits[6] = -23.7;
  dits[7] = 23.495;
  dits[8] = 28.5;
  dits[9] = -23.7;
  dits[10] = 14.595; 
  dits[11] = 28.5;
  dits[12] = 23.7;
  dits[13] = 14.595;
  dits[14] = 28.5;
  dits[15] = 23.7;
  dits[16] = 23.495;
  dits[17] = 28.5;
  dits[18] = 34.65;
  dits[19] = 23.495;
  dits[20] = 28.5;
  gMC->Gsvolu("IT34", "PCON", idtmed[200], dits, 21);  

  dits[0] = 3.2;
  dits[1] = 2;
  dits[2] = 34.65;
  gMC->Gsvolu("I048", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.63;
  dits[1] = 0.135;
  dits[2] = 30.385;
  gMC->Gsvolu("I005", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.2;
  dits[1] = 2;
  dits[2] = 23.7;
  gMC->Gsvolu("I047", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.63;
  dits[1] = 0.135;
  dits[2] = 23.05;
  gMC->Gsvolu("I004", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.2;
  dits[1] = 2;
  dits[2] = 2.725;
  gMC->Gsvolu("I024", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.2;
  dits[1] = 2;
  dits[2] = 3.65;
  gMC->Gsvolu("I018", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.625;
  dits[1] = 0.015;
  dits[2] = 4.382;
  gMC->Gsvolu("I302", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.625;
  dits[1] = 0.015;
  dits[2] = 4.382;
  gMC->Gsvolu("I402", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0.2;
  dits[1] = 0.1815;
  dits[2] = 2.725;
  dits[3] = 0.015;
  gMC->Gsvolu("I025", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 0.183;
  dits[1] = 0.165;
  dits[2] = 2.725;
  dits[3] = 0.015;
  gMC->Gsvolu("I026", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 2.23;
  dits[1] = 2.1;
  dits[2] = 0.05;
  dits[3] = 0.03;
  gMC->Gsvolu("I021", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 2.615;
  dits[1] = 2.465;
  dits[2] = 0.06;
  dits[3] = 0.04;
  gMC->Gsvolu("I023", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 2.1;
  dits[1] = 2;
  dits[2] = 0.06;
  dits[3] = 0.04;
  gMC->Gsvolu("I022", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 2.15;
  dits[1] = 0.2;
  dits[2] = 0.85;
  gMC->Gsvolu("I028", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 1.25;
  dits[1] = 0.6;
  dits[2] = 0.075;
  gMC->Gsvolu("I029", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 1.25;
  dits[1] = 0.1;
  dits[2] = 0.075;
  dits[3] = 1;
  gMC->Gsvolu("I030", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 1.6;
  dits[1] = 7;
  dits[2] = 0;
  dits[3] = 0.075;
  dits[4] = 0.775;
  dits[5] = 0.775;
  dits[6] = 0;
  dits[7] = 0.075;
  dits[8] = 0.376;
  dits[9] = 0.376;
  dits[10] = 0;
  gMC->Gsvolu("I027", "TRAP", idtmed[200], dits, 11);  

  dits[0] = 0;
  dits[1] = 0.093;
  dits[2] = 2.725;
  gMC->Gsvolu("I032", "TUBE", idtmed[200], dits, 3);  

  dits[0] = 0.093;
  dits[1] = 0.1;
  dits[2] = 2.725;
  gMC->Gsvolu("I031", "TUBE", idtmed[200], dits, 3);  

  dits[0] = 0.7;
  dits[1] = 0.002;
  dits[2] = 2.725;
  gMC->Gsvolu("I046", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0.2;
  dits[1] = 0.182;
  dits[2] = 3.65;
  dits[3] = 0.015;
  gMC->Gsvolu("I019", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 0.183;
  dits[1] = 0.165;
  dits[2] = 3.65;
  dits[3] = 0.015;
  gMC->Gsvolu("I020", "TRD1", idtmed[200], dits, 4);  

  dits[0] = 0.3;
  dits[1] = 0.05;
  dits[2] = 0.15;
  gMC->Gsvolu("I033", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0.2;
  dits[1] = 0.01;
  dits[2] = 0.05;
  gMC->Gsvolu("I036", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0;
  dits[1] = 0.05;
  dits[2] = 0.225;
  gMC->Gsvolu("I034", "TUBE", idtmed[200], dits, 3);  

  dits[0] = 0.1;
  dits[1] = 0.15;
  dits[2] = 0.2;
  gMC->Gsvolu("I035", "TUBE", idtmed[200], dits, 3);

  dits[0] = 0.7;
  dits[1] = 0.002;
  dits[2] = 3.65;
  gMC->Gsvolu("I045", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0;
  dits[1] = 0.093;
  dits[2] = 3.65;
  gMC->Gsvolu("I038", "TUBE", idtmed[200], dits, 3);  

  dits[0] = 0.093;
  dits[1] = 0.1;
  dits[2] = 3.65;
  gMC->Gsvolu("I037", "TUBE", idtmed[200], dits, 3);
  
  dits[0] = 1;
  dits[1] = 0.01;
  dits[2] = 3.6;
  gMC->Gsvolu("I039", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0.25;
  dits[1] = 0.01;
  dits[2] = 3.4;
  gMC->Gsvolu("I040", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0.1;
  dits[1] = 0.12;
  dits[2] = 3.4;
  dits[3] = 90;
  dits[4] = 320;
  gMC->Gsvolu("I041", "TUBS", idtmed[200], dits, 5);  

  dits[0] = 0.4;
  dits[1] = 0.015;
  dits[2] = 0.4;
  gMC->Gsvolu("I042", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0.25;
  dits[1] = 0.015;
  dits[2] = 0.25;
  gMC->Gsvolu("I043", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 0.75;
  dits[1] = 0.002;
  dits[2] = 3.4;
  gMC->Gsvolu("I044", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 1.8125;
  dits[1] = 3.81;
  dits[2] = 0;
  dits[3] = 0.015;
  dits[4] = 0.242;
  dits[5] = 0.242;
  dits[6] = 0;
  dits[7] = 0.015;
  dits[8] = 1E-03;
  dits[9] = 1E-03;
  dits[10] = 0;
  gMC->Gsvolu("I303", "TRAP", idtmed[200], dits, 11);  

  dits[0] = 1.8125;
  dits[1] = 3.81;
  dits[2] = 0;
  dits[3] = 0.015;
  dits[4] = 0.242;
  dits[5] = 0.242;
  dits[6] = 0;
  dits[7] = 0.015;
  dits[8] = 1E-03;
  dits[9] = 1E-03;
  dits[10] = 0;
  gMC->Gsvolu("I403", "TRAP", idtmed[200], dits, 11);  

  dits[0] = 3.5;
  dits[1] = 0.014;
  dits[2] = 3.763;
  gMC->Gsvolu("ITS3", "BOX ", idtmed[200], dits, 3);  

  dits[0] = 3.5;
  dits[1] = 0.014;
  dits[2] = 3.763;
  gMC->Gsvolu("ITS4", "BOX ", idtmed[200], dits, 3);  


  // --- Define SSD volumes ------------------------------------------

    
  dits[0] = 0;
  dits[1] = 360;
  dits[2] = 6;
  dits[3] = -57.5;
  dits[4] = 43.5;
  dits[5] = 48;  
  dits[6] = -51.365;
  dits[7] = 43.5;
  dits[8] = 48;  
  dits[9] = -51.365;
  dits[10] = 36.7;
  dits[11] = 48;  
  dits[12] = 51.3651;
  dits[13] = 36.7;
  dits[14] = 48;  
  dits[15] = 51.3651;
  dits[16] = 43.5;
  dits[17] = 48;  
  dits[18] = 56.96;
  dits[19] = 43.5;
  dits[20] = 48;   
  gMC->Gsvolu("IT56", "PCON", idtmed[200], dits, 21);   
  
  dits[0] =  3.4;
  dits[1] = 1.955;
  dits[2] = 57.13;
  gMC->Gsvolu("I570", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.75;
  dits[1] = 0.045;
  dits[2] = 50.975;
  gMC->Gsvolu("I569", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.4;
  dits[1] = 1.955;
  dits[2] = 57.13;
  gMC->Gsvolu("I571", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.75;
  dits[1] = 0.045;
  dits[2] = 45.21;
  gMC->Gsvolu("I565", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.4;
  dits[1] = 1.955;
  dits[2] = 3.15;
  gMC->Gsvolu("I553", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.405;
  dits[1] = 1.955;
  dits[2] = 1.955;
  gMC->Gsvolu("I523", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.75;
  dits[1] = 0.015;
  dits[2] = 2.1;
  gMC->Gsvolu("I566", "BOX ", idtmed[200], dits, 3); 
  
  dits[0] = 3.4;
  dits[1] = 1.955;
  dits[2] = 3.15;
  gMC->Gsvolu("I544", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.41;
  dits[1] = 1.955;
  dits[2] = 1.955;
  gMC->Gsvolu("I516", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.75;
  dits[1] = 0.015;
  dits[2] = 2.1;
  gMC->Gsvolu("I562", "BOX ", idtmed[200], dits, 3);   
  
  dits[0] = 0;
  dits[1] = 0.07;
  dits[2] = 3.15;
  gMC->Gsvolu("I559", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 0.07;
  dits[1] = 0.1;
  dits[2] = 3.15;
  gMC->Gsvolu("I560", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 0.225;
  dits[1] = 0.195;
  dits[2] = 3.15;
  dits[3] = 0.025;
  gMC->Gsvolu("I558", "TRD1", idtmed[200], dits, 4);  
  
  dits[0] = 0.25;
  dits[1] = 0.22;
  dits[2] = 3.15;
  dits[3] = 0.025;
  gMC->Gsvolu("I557", "TRD1", idtmed[200], dits, 4);  
  
  dits[0] = 2.17;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I556", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 2 ;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I554", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 2.675;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I555", "BOX ", idtmed[200], dits, 3); 
  
  dits[0] = 0.3;
  dits[1] = 0.15;
  dits[2] = 0.15;
  gMC->Gsvolu("I561", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.025;
  dits[1] = 0.025;
  dits[2] = 0.05;
  gMC->Gsvolu("I519", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.304;
  dits[1] = 0.0275;
  dits[2] = 0.432;
  gMC->Gsvolu("I521", "BOX ", idtmed[200], dits, 3);   
  
  dits[0] = 0.16;
  dits[1] = 0.08;
  dits[2] = 0.08;
  gMC->Gsvolu("I520", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.4;
  dits[1] = 0.015;
  dits[2] = 0.525;
  gMC->Gsvolu("I518", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.15;
  dits[1] = 0.105;
  dits[2] = 0.29;
  dits[3] = 0.08;
  gMC->Gsvolu("I522", "TRD1", idtmed[200], dits, 4);  
  
  dits[0] = 0.07;
  dits[1] = 0.1;
  dits[2] = 1.955;
  gMC->Gsvolu("I542", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 0;
  dits[1] = 0.07;
  dits[2] = 1.955;
  gMC->Gsvolu("I541", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 0.3;
  dits[1] = 0.15;
  dits[2] = 0.15;
  gMC->Gsvolu("I543", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.25;
  dits[1] = 0.22;
  dits[2] = 1.955;
  dits[3] = 0.025;
  gMC->Gsvolu("I537", "TRD1", idtmed[200], dits, 4); 
  
  dits[0] = 0.225;
  dits[1] = 0.195;
  dits[2] = 1.955;
  dits[4] = 0.025;
  gMC->Gsvolu("I538", "TRD1", idtmed[200], dits, 4);  
  
  dits[0] = 2.17;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I536", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 2.675;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I535", "BOX ", idtmed[200], dits, 3);   
  
  dits[0] = 2;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I534", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0;
  dits[1] = 0.05;
  dits[2] = 0.17;
  gMC->Gsvolu("I540", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 0;
  dits[1] = 0.05;
  dits[2] = 0.205;
  gMC->Gsvolu("I539", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 3.65;
  dits[1] = 0.015;
  dits[2] = 2;
  gMC->Gsvolu("ITS6", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0;
  dits[1] = 0.07;
  dits[2] = 3.15;
  gMC->Gsvolu("I550", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 0.07;
  dits[1] = 0.1;
  dits[2] = 3.15;
  gMC->Gsvolu("I551", "TUBE", idtmed[200], dits, 3);  
  
  dits[0] = 0.225;
  dits[1] = 0.195;
  dits[2] = 3.15;
  dits[3] = 0.025;
  gMC->Gsvolu("I549", "TRD1", idtmed[200], dits, 4); 
  
  dits[0] = 0.25;
  dits[1] = 0.22;
  dits[2] = 3.15;
  dits[3] = 0.025;
  gMC->Gsvolu("I548", "TRD1", idtmed[200], dits, 4);  
  
  dits[0] = 2.17;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I547", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 2;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I545", "BOX ", idtmed[200], dits, 3);   
  
  dits[0] = 2.675;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I546", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.3;
  dits[1] = 0.15;
  dits[2] = 0.15;
  gMC->Gsvolu("I552", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.304;
  dits[1] = 0.0275;
  dits[2] = 0.4322;
  gMC->Gsvolu("I515", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.025;
  dits[1] = 0.025;
  dits[2] = 0.05;
  gMC->Gsvolu("I513", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.16;
  dits[1] = 0.08;
  dits[2] = 0.08;
  gMC->Gsvolu("I514", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.4;
  dits[1] = 0.015;
  dits[2] = 0.525;
  gMC->Gsvolu("I512", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 0.225;
  dits[1] = 0.195;
  dits[2] = 1.955;
  dits[3] = 0.025;
  gMC->Gsvolu("I528", "TRD1", idtmed[200], dits, 4); 
  
  dits[0] = 0.25;
  dits[1] = 0.22;
  dits[2] = 1.955;
  dits[3] = 0.025;
  gMC->Gsvolu("I527", "TRD1", idtmed[200], dits, 4);  
  
  dits[0] = 2.17;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I526", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 2.675;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I525", "BOX ", idtmed[200], dits, 3);  
   
  dits[0] = 2;
  dits[1] = 0.035;
  dits[2] = 0.05;
  gMC->Gsvolu("I524", "BOX ", idtmed[200], dits, 3);  
   
  dits[0] = 0;
  dits[1] = 0.05;
  dits[2] = 0.205;
  gMC->Gsvolu("I529", "TUBE", idtmed[200], dits, 3);  
   
  dits[0] = 0;
  dits[1] = 0.05;
  dits[2] = 0.17;
  gMC->Gsvolu("I530", "TUBE", idtmed[200], dits, 3);  
   
  dits[0] = 0.15;
  dits[1] = 0.105;
  dits[2] = 0.29;
  dits[3] = 0.08;
  gMC->Gsvolu("I517", "TRD1", idtmed[200], dits, 4);  
  
  dits[0] = 0;
  dits[1] = 0.07;
  dits[2] = 1.955;
  gMC->Gsvolu("I531", "TUBE", idtmed[200], dits, 3);  
     
  dits[0] = 0.07;
  dits[1] = 0.1;
  dits[2] = 1.955;
  gMC->Gsvolu("I532", "TUBE", idtmed[200], dits, 3);  
 
  dits[0] = 0.3;
  dits[1] = 0.15;
  dits[2] = 0.15;
  gMC->Gsvolu("I533", "BOX ", idtmed[200], dits, 3);  
  
  dits[0] = 3.65;
  dits[1] = 0.015;
  dits[2] = 2;
  gMC->Gsvolu("ITS5", "BOX ", idtmed[200], dits, 3);  


     
  // --- Place SPD (version 'a') volumes into their mother volume IT12

  gMC->Gspos("I12A",5,"IT12",0.0,0.0,0.0,idrotm[238],"MANY");
  gMC->Gspos("I12A",6,"IT12",0.0,0.0,0.0,idrotm[236],"MANY");
  gMC->Gspos("I12A",7,"IT12",0.0,0.0,0.0,idrotm[239],"MANY");
  gMC->Gspos("I12A",8,"IT12",0.0,0.0,0.0,idrotm[233],"MANY");
  gMC->Gspos("I12A",9,"IT12",0.0,0.0,0.0,idrotm[240],"MANY");
  gMC->Gspos("I12A",10,"IT12",0.0,0.0,0.0,idrotm[241],"MANY");
  gMC->Gspos("I12A",2,"IT12",0.0,0.0,0.0,idrotm[242],"MANY");
  gMC->Gspos("I12A",3,"IT12",0.0,0.0,0.0,idrotm[234],"MANY");
  gMC->Gspos("I12A",4,"IT12",0.0,0.0,0.0,idrotm[243],"MANY");
  gMC->Gspos("I12A",1,"IT12",0.0,0.0,0.0,0,"MANY");
  gMC->Gspos("I10A",2,"I12A",0.203,3.8206,0.0,idrotm[244],"ONLY");
  gMC->Gspos("I10A",1,"I12A",1.4531,3.8152,0.0,idrotm[245],"ONLY");
  gMC->Gspos("I20A",1,"I12A",3.0174,6.5143,0.0,idrotm[246],"ONLY");
  gMC->Gspos("I20A",2,"I12A",1.9612,6.9062,0.0,idrotm[247],"ONLY");
  gMC->Gspos("I20A",3,"I12A",0.8567,7.1279,0.0,idrotm[248],"ONLY");
  gMC->Gspos("I20A",4,"I12A",-0.2689,7.1742,0.0,idrotm[249],"ONLY");
  gMC->Gspos("I123",2,"I12A",-0.2978,5.5196,0.0,idrotm[214],"ONLY");
  gMC->Gspos("I121",2,"I12A",-0.2385,4.1518,0.0,idrotm[213],"ONLY");
  gMC->Gspos("I122",2,"I12A",-0.2968,4.0207,0.0,idrotm[212],"ONLY");
  gMC->Gspos("I120",2,"I12A",-0.3672,3.9056,0.0,0,"ONLY");
  gMC->Gspos("I144",1,"I12A",-0.2538,3.8556,0.0,0,"ONLY");
  gMC->Gspos("I113",3,"I12A",0.1095,3.9056,0.0,0,"ONLY");
  gMC->Gspos("I143",1,"I12A",0.4365,3.8556,0.0,idrotm[236],"ONLY");
  gMC->Gspos("I142",1,"I12A",0.5136,3.9056,0.0,idrotm[235],"ONLY");
  gMC->Gspos("I141",1,"I12A",0.5636,3.9752,0.0,idrotm[201],"ONLY");
  gMC->Gspos("I140",1,"I12A",0.6336,4.0447,0.0,idrotm[234],"ONLY");
  gMC->Gspos("I139",1,"I12A",0.8297,4.0545,0.0,idrotm[207],"ONLY");
  gMC->Gspos("I113",5,"I12A",1.2575,3.9681,0.0,idrotm[207],"ONLY");
  gMC->Gspos("I138",1,"I12A",1.66,3.7848,0.0,idrotm[207],"ONLY");
  gMC->Gspos("I137",1,"I12A",1.8556,3.7738,0.0,idrotm[233],"ONLY");
  gMC->Gspos("I136",1,"I12A",2.6224,4.874,0.0,idrotm[232],"ONLY");
  gMC->Gspos("I135",1,"I12A",3.2967,6.0337,0.0,idrotm[231],"ONLY");
  gMC->Gspos("I134",1,"I12A",3.266,6.1636,0.0,idrotm[230],"ONLY");
  gMC->Gspos("I113",1,"I12A",2.9903,6.4144,0.0,idrotm[211],"ONLY");
  gMC->Gspos("I133",3,"I12A",2.7631,6.7627,0.0,idrotm[230],"ONLY");
  gMC->Gspos("I132",3,"I12A",2.62,6.8555,0.0,idrotm[229],"ONLY");
  gMC->Gspos("I131",3,"I12A",2.648,6.6023,0.0,idrotm[228],"ONLY");
  gMC->Gspos("I130",3,"I12A",2.6569,6.3431,0.0,idrotm[227],"ONLY");
  gMC->Gspos("I129",3,"I12A",2.3906,6.4819,0.0,idrotm[226],"ONLY");
  gMC->Gspos("I113",2,"I12A",1.9488,6.7998,0.0,idrotm[210],"ONLY");
  gMC->Gspos("I133",2,"I12A",1.6699,7.1085,0.0,idrotm[226],"ONLY");
  gMC->Gspos("I132",2,"I12A",1.5142,7.1777,0.0,idrotm[225],"ONLY");
  gMC->Gspos("I131",2,"I12A",1.5814,6.932,0.0,idrotm[224],"ONLY");
  gMC->Gspos("I130",2,"I12A",1.6308,6.6774,0.0,idrotm[223],"ONLY");
  gMC->Gspos("I129",2,"I12A",1.346,6.7728,0.0,idrotm[222],"ONLY");
  gMC->Gspos("I113",6,"I12A",0.8599,7.0176,0.0,idrotm[209],"ONLY");
  gMC->Gspos("I133",1,"I12A",0.5362,7.2789,0.0,idrotm[222],"ONLY");
  gMC->Gspos("I132",1,"I12A",0.3715,7.3228,0.0,idrotm[221],"ONLY");
  gMC->Gspos("I131",1,"I12A",0.4763,7.0907,0.0,idrotm[220],"ONLY");
  gMC->Gspos("I130",1,"I12A",0.5649,6.8469,0.0,idrotm[219],"ONLY");
  gMC->Gspos("I129",1,"I12A",0.2688,6.8966,0.0,idrotm[218],"ONLY");
  gMC->Gspos("I113",4,"I12A",-0.2497,7.0624,0.0,idrotm[208],"ONLY");
  gMC->Gspos("I128",1,"I12A",-0.6103,7.2698,0.0,idrotm[218],"ONLY");
  gMC->Gspos("I126",2,"I12A",-0.7799,7.2874,0.0,idrotm[217],"ONLY");
  gMC->Gspos("I125",2,"I12A",-0.6315,7.0883,0.0,idrotm[216],"ONLY");
  gMC->Gspos("I124",2,"I12A",-0.4965,6.8742,0.0,idrotm[215],"ONLY");
  gMC->Gspos("I103",3,"I10A",-0.05,0.0075,-3.536,idrotm[237],"ONLY");
  gMC->Gspos("I103",4,"I10A",-0.05,0.0075,-10.708,idrotm[237],"ONLY");
  gMC->Gspos("I103",1,"I10A",-0.05,0.0075,10.708,0,"ONLY");
  gMC->Gspos("I103",2,"I10A",-0.05,0.0075,3.536,0,"ONLY");
  gMC->Gspos("I105",1,"I10A",-0.05,0.01,-16.844,idrotm[237],"ONLY");
  gMC->Gspos("I105",2,"I10A",-0.05,0.01,16.844,0,"ONLY");
  gMC->Gspos("I104",1,"I10A",0.0,-0.015,0.0,0,"ONLY");
  gMC->Gspos("I1D3",1,"I20A",-0.05,0.0075,-3.536,idrotm[237],"ONLY");
  gMC->Gspos("I1D3",2,"I20A",-0.05,0.0075,-10.708,idrotm[237],"ONLY");
  gMC->Gspos("I1D3",3,"I20A",-0.05,0.0075,10.708,0,"ONLY");
  gMC->Gspos("I1D3",4,"I20A",-0.05,0.0075,3.536,0,"ONLY");
  gMC->Gspos("I105",3,"I20A",-0.05,0.01,-16.844,idrotm[237],"ONLY");
  gMC->Gspos("I105",4,"I20A",-0.05,0.01,16.844,0,"ONLY");
  gMC->Gspos("I104",2,"I20A",0.0,-0.015,0.0,0,"ONLY");
  gMC->Gspos("I112",2,"I113",0.25,0.02,0.0,idrotm[206],"ONLY");
  gMC->Gspos("I111",2,"I113",0.1318,-0.0008,0.0,idrotm[205],"ONLY");
  gMC->Gspos("I118",1,"I113",0.0,-0.0454,0.0,0,"ONLY");
  gMC->Gspos("I110",1,"I113",0.0,0.0492,0.0,0,"ONLY");
  gMC->Gspos("I114",1,"I113",0.063,0.0042,0.0,idrotm[202],"ONLY");
  gMC->Gspos("I115",1,"I113",0.063,0.0042,0.0,idrotm[202],"ONLY");
  gMC->Gspos("I115",2,"I113",-0.063,0.0042,0.0,idrotm[201],"ONLY");
  gMC->Gspos("I114",2,"I113",-0.063,0.0042,0.0,idrotm[201],"ONLY");
  gMC->Gspos("I116",1,"I113",0.0,0.0042,0.0,0,"ONLY");
  gMC->Gspos("I111",1,"I113",-0.1318,-0.0008,0.0,idrotm[204],"ONLY");
  gMC->Gspos("I112",1,"I113",-0.25,0.02,0.0,idrotm[203],"ONLY");
  gMC->Gspos("I101",1,"I103",-0.088,0.0075,0.0,0,"ONLY");
  gMC->Gspos("I102",1,"I103",0.0,-0.005,-2.8,0,"ONLY");
  gMC->Gspos("I102",2,"I103",0.0,-0.005,-1.4,0,"ONLY");
  gMC->Gspos("I102",3,"I103",0.0,-0.005,0.0,0,"ONLY");
  gMC->Gspos("I102",4,"I103",0.0,-0.005,1.4,0,"ONLY");
  gMC->Gspos("I102",5,"I103",0.0,-0.005,2.8,0,"ONLY");
  gMC->Gspos("I1D1",1,"I1D3",-0.088,0.0075,0.0,0,"ONLY");
  gMC->Gspos("I102",6,"I1D3",0.0,-0.005,-2.8,0,"ONLY");
  gMC->Gspos("I102",7,"I1D3",0.0,-0.005,-1.4,0,"ONLY");
  gMC->Gspos("I102",8,"I1D3",0.0,-0.005,0.0,0,"ONLY");
  gMC->Gspos("I102",9,"I1D3",0.0,-0.005,1.4,0,"ONLY");
  gMC->Gspos("I102",10,"I1D3",0.0,-0.005,2.8,0,"ONLY");
  gMC->Gspos("I117",1,"I116",0.0,0.0,0.0,0,"ONLY");
  gMC->Gspos("ITS1",1,"I101",0.0,0.0,0.0,0,"ONLY");
  gMC->Gspos("ITS2",1,"I1D1",0.0,0.0,0.0,0,"ONLY");
    
  // --- Place SDD volumes into their mother volume IT34
  
  gMC->Gspos("I048",8,"IT34",-22.1376,-14.227,0.0,idrotm[328],"ONLY");
  gMC->Gspos("I048",7,"IT34",-24.7213,-7.2588,0.0,idrotm[329],"ONLY");
  gMC->Gspos("I048",6,"IT34",-26.315,0.0,0.0,idrotm[350],"ONLY");
  gMC->Gspos("I048",5,"IT34",-24.7213,7.2588,0.0,idrotm[330],"ONLY");
  gMC->Gspos("I048",4,"IT34",-22.1376,14.227,0.0,idrotm[366],"ONLY");
  gMC->Gspos("I048",3,"IT34",-16.8725,19.4719,0.0,idrotm[331],"ONLY");
  gMC->Gspos("I048",2,"IT34",-10.9317,23.937,0.0,idrotm[332],"ONLY");
  gMC->Gspos("I048",1,"IT34",-3.6667,25.5027,0.0,idrotm[335],"ONLY");
  gMC->Gspos("I048",22,"IT34",3.745,26.0472,0.0,idrotm[334],"ONLY");
  gMC->Gspos("I048",21,"IT34",10.7032,23.4367,0.0,idrotm[314],"ONLY");
  gMC->Gspos("I048",20,"IT34",17.2327,19.8876,0.0,idrotm[315],"ONLY");
  gMC->Gspos("I048",19,"IT34",21.6749,13.9296,0.0,idrotm[316],"ONLY");
  gMC->Gspos("I048",18,"IT34",25.2491,7.4138,0.0,idrotm[317],"ONLY");
  gMC->Gspos("I048",17,"IT34",25.765,0.0,0.0,idrotm[318],"ONLY");
  gMC->Gspos("I048",16,"IT34",25.2491,-7.4138,0.0,idrotm[319],"ONLY");
  gMC->Gspos("I048",15,"IT34",21.6749,-13.9296,0.0,idrotm[320],"ONLY");
  gMC->Gspos("I048",14,"IT34",17.2327,-19.8876,0.0,idrotm[322],"ONLY");
  gMC->Gspos("I048",13,"IT34",10.7032,-23.4367,0.0,idrotm[323],"ONLY");
  gMC->Gspos("I048",12,"IT34",3.745,-26.0472,0.0,idrotm[324],"ONLY");
  gMC->Gspos("I048",11,"IT34",-3.6667,-25.5027,0.0,idrotm[325],"ONLY");
  gMC->Gspos("I048",10,"IT34",-10.9316,-23.937,0.0,idrotm[326],"ONLY");
  gMC->Gspos("I048",9,"IT34",-16.8725,-19.4719,0.0,idrotm[327],"ONLY");
  gMC->Gspos("I005",9,"IT34",-15.4744,-17.8584,-0.15,idrotm[327],"ONLY");
  gMC->Gspos("I005",8,"IT34",-20.3415,-13.0727,-0.15,idrotm[328],"ONLY");
  gMC->Gspos("I005",7,"IT34",-22.6728,-6.6573,-0.15,idrotm[329],"ONLY");
  gMC->Gspos("I005",6,"IT34",-24.18,0.0,-0.15,idrotm[350],"ONLY");
  gMC->Gspos("I005",5,"IT34",-22.6728,6.6573,-0.15,idrotm[330],"ONLY");
  gMC->Gspos("I005",4,"IT34",-20.3415,13.0727,-0.15,idrotm[366],"ONLY");
  gMC->Gspos("I005",3,"IT34",-15.4744,17.8584,-0.15,idrotm[331],"ONLY");
  gMC->Gspos("I005",2,"IT34",-10.0447,21.9949,-0.15,idrotm[332],"ONLY");
  gMC->Gspos("I005",1,"IT34",-3.3629,23.3895,-0.15,idrotm[335],"ONLY");
  gMC->Gspos("I005",22,"IT34",3.4412,23.9339,-0.15,idrotm[334],"ONLY");
  gMC->Gspos("I005",21,"IT34",9.8163,21.4946,-0.15,idrotm[314],"ONLY");
  gMC->Gspos("I005",20,"IT34",15.8345,18.274,-0.15,idrotm[315],"ONLY");
  gMC->Gspos("I005",19,"IT34",19.8788,12.7753,-0.15,idrotm[316],"ONLY");
  gMC->Gspos("I005",18,"IT34",23.2005,6.8123,-0.15,idrotm[317],"ONLY");
  gMC->Gspos("I005",17,"IT34",23.63,0.0,-0.15,idrotm[318],"ONLY");
  gMC->Gspos("I005",16,"IT34",23.2005,-6.8123,-0.15,idrotm[319],"ONLY");
  gMC->Gspos("I005",15,"IT34",19.8788,-12.7753,-0.15,idrotm[320],"ONLY");
  gMC->Gspos("I005",14,"IT34",15.8345,-18.274,-0.15,idrotm[322],"ONLY");
  gMC->Gspos("I005",13,"IT34",9.8163,-21.4946,-0.15,idrotm[323],"ONLY");
  gMC->Gspos("I005",12,"IT34",3.4412,-23.9339,-0.15,idrotm[324],"ONLY");
  gMC->Gspos("I005",11,"IT34",-3.3629,-23.3895,-0.15,idrotm[325],"ONLY");
  gMC->Gspos("I005",10,"IT34",-10.0447,-21.9949,-0.15,idrotm[326],"ONLY");
  gMC->Gspos("I047",6,"IT34",-10.8893,-13.6547,0.0,idrotm[311],"ONLY");
  gMC->Gspos("I047",5,"IT34",-15.1948,-7.3174,0.0,idrotm[313],"ONLY");
  gMC->Gspos("I047",4,"IT34",-17.465,0.0,0.0,idrotm[350],"ONLY");
  gMC->Gspos("I047",3,"IT34",-15.1948,7.3175,0.0,idrotm[336],"ONLY");
  gMC->Gspos("I047",2,"IT34",-10.8892,13.6547,0.0,idrotm[333],"ONLY");
  gMC->Gspos("I047",1,"IT34",-3.7528,16.4422,0.0,idrotm[321],"ONLY");
  gMC->Gspos("I047",14,"IT34",3.8863,17.0271,0.0,idrotm[305],"ONLY");
  gMC->Gspos("I047",13,"IT34",10.5152,13.1856,0.0,idrotm[306],"ONLY");
  gMC->Gspos("I047",12,"IT34",15.7354,7.5778,0.0,idrotm[307],"ONLY");
  gMC->Gspos("I047",11,"IT34",16.865,0.0,0.0,idrotm[356],"ONLY");
  gMC->Gspos("I047",10,"IT34",15.7354,-7.5778,0.0,idrotm[308],"ONLY");
  gMC->Gspos("I047",9,"IT34",10.5152,-13.1856,0.0,idrotm[309],"ONLY");
  gMC->Gspos("I047",8,"IT34",3.8863,-17.0271,0.0,idrotm[386],"ONLY");
  gMC->Gspos("I047",7,"IT34",-3.7528,-16.4422,0.0,idrotm[310],"ONLY");
  gMC->Gspos("I004",6,"IT34",-9.5581,-11.9855,0.0,idrotm[311],"ONLY");
  gMC->Gspos("I004",5,"IT34",-13.2713,-6.3911,0.0,idrotm[313],"ONLY");
  gMC->Gspos("I004",4,"IT34",-15.33,0.0,0.0,idrotm[350],"ONLY");
  gMC->Gspos("I004",3,"IT34",-13.2713,6.3911,0.0,idrotm[336],"ONLY");
  gMC->Gspos("I004",2,"IT34",-9.5581,11.9855,0.0,idrotm[333],"ONLY");
  gMC->Gspos("I004",1,"IT34",-3.2777,14.3607,0.0,idrotm[321],"ONLY");
  gMC->Gspos("I004",14,"IT34",3.4113,14.9456,0.0,idrotm[305],"ONLY");
  gMC->Gspos("I004",13,"IT34",9.184,11.5164,0.0,idrotm[306],"ONLY");
  gMC->Gspos("I004",12,"IT34",13.8119,6.6514,0.0,idrotm[307],"ONLY");
  gMC->Gspos("I004",11,"IT34",14.73,0.0,0.0,idrotm[356],"ONLY");
  gMC->Gspos("I004",10,"IT34",13.8119,-6.6514,0.0,idrotm[308],"ONLY");
  gMC->Gspos("I004",9,"IT34",9.184,-11.5164,0.0,idrotm[309],"ONLY");
  gMC->Gspos("I004",8,"IT34",3.4112,-14.9456,0.0,idrotm[386],"ONLY");
  gMC->Gspos("I004",7,"IT34",-3.2777,-14.3607,0.0,idrotm[310],"ONLY");
  gMC->Gspos("I024",3,"I048",-0.0001,0.0,31.925,0,"ONLY");
  gMC->Gspos("I024",4,"I048",-0.0001,0.0,-31.925,idrotm[355],"ONLY");
  gMC->Gspos("I018",13,"I048",-0.0001,0.0,-25.55,0,"ONLY");
  gMC->Gspos("I018",12,"I048",-0.0001,0.0,-18.25,0,"ONLY");
  gMC->Gspos("I018",11,"I048",-0.0001,0.0,-10.95,0,"ONLY");
  gMC->Gspos("I018",10,"I048",-0.0001,0.0,25.55,0,"ONLY");
  gMC->Gspos("I018",9,"I048",-0.0001,0.0,18.25,0,"ONLY");
  gMC->Gspos("I018",8,"I048",-0.0001,0.0,10.95,0,"ONLY");
  gMC->Gspos("I018",7,"I048",-0.0001,0.0,3.65,0,"ONLY");
  gMC->Gspos("I018",6,"I048",-0.0001,0.0,3.65,0,"ONLY");
  gMC->Gspos("I402",5,"I005",0.0,-0.115,-3.55,0,"ONLY");
  gMC->Gspos("I402",4,"I005",0.0,0.115,3.85,0,"ONLY");
  gMC->Gspos("I402",2,"I005",0.0,0.115,18.75,0,"ONLY");
  gMC->Gspos("I402",3,"I005",0.0,-0.115,11.15,0,"ONLY");
  gMC->Gspos("I402",1,"I005",0.0,-0.115,25.9,0,"ONLY");
  gMC->Gspos("I402",6,"I005",0.0,0.115,-11.05,0,"ONLY");
  gMC->Gspos("I402",7,"I005",0.0,-0.115,-18.3,0,"ONLY");
  gMC->Gspos("I402",8,"I005",0.0,0.115,-25.9,0,"ONLY");
  gMC->Gspos("I024",1,"I047",0.0,0.0,20.975,0,"ONLY");
  gMC->Gspos("I018",4,"I047",0.0,0.0,7.3,0,"ONLY");
  gMC->Gspos("I018",5,"I047",0.0,0.0,14.6,0,"ONLY");
  gMC->Gspos("I018",1,"I047",0.0,0.0,0.0,0,"ONLY");
  gMC->Gspos("I018",3,"I047",0.0,0.0,-14.6,0,"ONLY");
  gMC->Gspos("I018",2,"I047",0.0,0.0,-7.3,0,"ONLY");
  gMC->Gspos("I024",2,"I047",0.0,0.0,-20.975,idrotm[355],"ONLY");
  gMC->Gspos("I302",4,"I004",0.0,-0.115,-3.7,0,"ONLY");
  gMC->Gspos("I302",3,"I004",0.0,0.115,3.7,0,"ONLY");
  gMC->Gspos("I302",6,"I004",0.0,-0.115,-18.35,0,"ONLY");
  gMC->Gspos("I302",5,"I004",0.0,0.115,-11.2,0,"ONLY");
  gMC->Gspos("I302",2,"I004",0.0,-0.115,10.95,0,"ONLY");
  gMC->Gspos("I302",1,"I004",0.0,0.115,18.55,0,"ONLY");
  gMC->Gspos("I025",2,"I024",1.987,-1.5842,0.0,idrotm[343],"ONLY");
  gMC->Gspos("I026",2,"I024",1.8824,-1.7349,0.0,idrotm[344],"ONLY");
  gMC->Gspos("I025",1,"I024",-1.9,-1.7349,0.0,idrotm[344],"ONLY");
  gMC->Gspos("I026",1,"I024",-1.9782,-1.5689,0.0,idrotm[342],"ONLY");
  gMC->Gspos("I026",3,"I024",0.0958,1.6914,0.0,idrotm[343],"ONLY");
  gMC->Gspos("I025",3,"I024",-0.087,1.7067,0.0,idrotm[342],"ONLY");
  gMC->Gspos("I021",10,"I024",1.0761,0.0836,1.7742,idrotm[337],"ONLY");
  gMC->Gspos("I021",9,"I024",-1.0761,0.0836,1.7742,idrotm[339],"ONLY");
  gMC->Gspos("I021",12,"I024",1.0761,0.0836,-0.1242,idrotm[340],"ONLY");
  gMC->Gspos("I021",11,"I024",-1.0761,0.0836,-0.1242,idrotm[338],"ONLY");
  gMC->Gspos("I021",13,"I024",-1.0761,0.0836,-1.8758,idrotm[339],"ONLY");
  gMC->Gspos("I021",14,"I024",1.0761,0.0836,-1.8758,idrotm[337],"ONLY");
  gMC->Gspos("I023",3,"I024",0.0,-1.7899,-1.0,idrotm[341],"ONLY");
  gMC->Gspos("I022",3,"I024",0.0,-1.7899,0.825,idrotm[312],"ONLY");
  gMC->Gspos("I028",1,"I024",0.0,-1.7999,1.875,0,"MANY");
  gMC->Gspos("I029",1,"I024",0.0,-0.9999,2.65,0,"ONLY");
  gMC->Gspos("I030",1,"I024",0.0,0.6001,2.65,idrotm[344],"ONLY");
  gMC->Gspos("I027",1,"I024",0.0,0.0001,1.9965,idrotm[352],"ONLY");
  gMC->Gspos("I032",1,"I024",1.7,-0.4999,0.0,0,"ONLY");
  gMC->Gspos("I031",1,"I024",1.7,-0.4999,0.0,0,"ONLY");
  gMC->Gspos("I031",2,"I024",-1.7,-0.4999,0.0,0,"ONLY");
  gMC->Gspos("I032",2,"I024",-1.7,-0.4999,0.0,0,"ONLY");
  gMC->Gspos("I046",6,"I024",-0.616,1.1702,0.0,idrotm[353],"ONLY");
  gMC->Gspos("I046",5,"I024",-0.566,1.1702,0.0,idrotm[353],"ONLY");
  gMC->Gspos("I046",4,"I024",0.616,1.1702,0.0,idrotm[354],"ONLY");
  gMC->Gspos("I046",3,"I024",0.566,1.1702,0.0,idrotm[354],"ONLY");
  gMC->Gspos("I046",2,"I024",0.516,1.1702,0.0,idrotm[354],"ONLY");
  gMC->Gspos("I046",1,"I024",-0.516,1.1702,0.0,idrotm[353],"ONLY");
  gMC->Gspos("I022",2,"I018",0.0,-1.79,-0.1,idrotm[312],"ONLY");
  gMC->Gspos("I021",8,"I018",1.0761,0.0835,0.8492,idrotm[337],"ONLY");
  gMC->Gspos("I021",7,"I018",-1.0761,0.0835,2.6008,idrotm[338],"ONLY");
  gMC->Gspos("I021",6,"I018",-1.0761,0.0835,0.8492,idrotm[339],"ONLY");
  gMC->Gspos("I021",5,"I018",1.0761,0.0835,-1.0492,idrotm[340],"ONLY");
  gMC->Gspos("I021",4,"I018",1.0761,0.0835,-2.8008,idrotm[337],"ONLY");
  gMC->Gspos("I021",3,"I018",-1.0761,0.0835,-1.0492,idrotm[338],"ONLY");
  gMC->Gspos("I021",2,"I018",-1.0761,0.0835,-2.8008,idrotm[339],"ONLY");
  gMC->Gspos("I023",2,"I018",0.0,-1.79,-1.925,idrotm[341],"ONLY");
  gMC->Gspos("I019",3,"I018",-0.087,1.7066,0.0,idrotm[342],"ONLY");
  gMC->Gspos("I020",3,"I018",0.0958,1.6913,0.0,idrotm[343],"ONLY");
  gMC->Gspos("I019",2,"I018",1.987,-1.5843,0.0,idrotm[343],"ONLY");
  gMC->Gspos("I020",2,"I018",1.8824,-1.735,0.0,idrotm[344],"ONLY");
  gMC->Gspos("I022",1,"I018",0.0,-1.79,3.55,idrotm[312],"ONLY");
  gMC->Gspos("I021",1,"I018",1.0761,0.0835,2.6008,idrotm[340],"ONLY");
  gMC->Gspos("I023",1,"I018",0.0,-1.79,1.725,idrotm[341],"ONLY");
  gMC->Gspos("I019",1,"I018",-1.9,-1.735,0.0,idrotm[344],"ONLY");
  gMC->Gspos("I020",1,"I018",-1.9782,-1.569,0.0,idrotm[342],"ONLY");
  gMC->Gspos("I033",1,"I018",1.8,-1.75,1.35,0,"MANY");
  gMC->Gspos("I033",4,"I018",1.8,-1.75,-2.65,0,"MANY");
  gMC->Gspos("I033",2,"I018",-1.8,-1.75,-2.65,idrotm[345],"MANY");
  gMC->Gspos("I033",3,"I018",-1.8,-1.75,1.35,idrotm[345],"MANY");
  gMC->Gspos("I036",1,"I018",0.3087,1.7191,3.56,idrotm[346],"ONLY");
  gMC->Gspos("I036",4,"I018",-0.3087,1.7191,3.56,idrotm[347],"ONLY");
  gMC->Gspos("I036",2,"I018",0.3087,1.7191,-0.11,idrotm[346],"ONLY");
  gMC->Gspos("I036",3,"I018",-0.3087,1.7191,-0.11,idrotm[347],"ONLY");
  gMC->Gspos("I034",1,"I018",1.6,-1.775,1.35,idrotm[312],"ONLY");
  gMC->Gspos("I034",4,"I018",1.6,-1.775,-2.65,idrotm[312],"ONLY");
  gMC->Gspos("I034",2,"I018",-1.6,-1.775,-2.65,idrotm[348],"ONLY");
  gMC->Gspos("I034",3,"I018",-1.6,-1.775,1.35,idrotm[348],"ONLY");
  gMC->Gspos("I035",2,"I018",-1.7,-0.55,2.8581,idrotm[345],"MANY");
  gMC->Gspos("I035",1,"I018",1.7,-0.55,2.8581,0,"MANY");
  gMC->Gspos("I045",1,"I018",0.7483,0.9337,0.0,idrotm[346],"ONLY");
  gMC->Gspos("I045",2,"I018",0.7065,0.9337,0.0,idrotm[346],"ONLY");
  gMC->Gspos("I045",3,"I018",-0.7483,0.9337,0.0,idrotm[347],"ONLY");
  gMC->Gspos("I045",4,"I018",-0.7065,0.9337,0.0,idrotm[347],"ONLY");
  gMC->Gspos("I038",1,"I018",1.7,-0.55,0.0,idrotm[346],"ONLY");
  gMC->Gspos("I037",1,"I018",1.7,-0.55,0.0,idrotm[346],"ONLY");
  gMC->Gspos("I037",2,"I018",-1.7,-0.55,0.0,idrotm[347],"ONLY");
  gMC->Gspos("I038",2,"I018",-1.7,-0.55,0.0,idrotm[347],"ONLY");
  gMC->Gspos("I039",1,"I018",1.8126,-0.485,0.0,idrotm[346],"ONLY");
  gMC->Gspos("I040",1,"I018",1.9204,-0.7118,0.0,idrotm[346],"ONLY");
  gMC->Gspos("I041",1,"I018",1.7,-0.55,0.0,idrotm[346],"ONLY");
  gMC->Gspos("I042",1,"I018",2.0342,-0.8189,3.12,idrotm[346],"ONLY");
  gMC->Gspos("I042",2,"I018",2.0342,-0.8189,2.28,idrotm[346],"ONLY");
  gMC->Gspos("I042",3,"I018",2.0342,-0.8189,1.38,idrotm[346],"ONLY");
  gMC->Gspos("I042",4,"I018",2.0342,-0.8189,0.48,idrotm[346],"ONLY");
  gMC->Gspos("I042",5,"I018",2.0342,-0.8189,-0.42,idrotm[346],"ONLY");
  gMC->Gspos("I042",6,"I018",2.0342,-0.8189,-1.32,idrotm[346],"ONLY");
  gMC->Gspos("I042",7,"I018",2.0342,-0.8189,-2.22,idrotm[346],"ONLY");
  gMC->Gspos("I042",8,"I018",2.0342,-0.8189,-3.12,idrotm[346],"ONLY");
  gMC->Gspos("I043",8,"I018",1.5592,0.0038,-3.15,idrotm[346],"ONLY");
  gMC->Gspos("I043",7,"I018",1.5592,0.0038,-2.25,idrotm[346],"ONLY");
  gMC->Gspos("I043",6,"I018",1.5592,0.0038,-1.35,idrotm[346],"ONLY");
  gMC->Gspos("I043",5,"I018",1.5592,0.0038,-0.45,idrotm[346],"ONLY");
  gMC->Gspos("I043",4,"I018",1.5592,0.0038,0.45,idrotm[346],"ONLY");
  gMC->Gspos("I043",3,"I018",1.5592,0.0038,1.35,idrotm[346],"ONLY");
  gMC->Gspos("I043",2,"I018",1.5592,0.0038,2.25,idrotm[346],"ONLY");
  gMC->Gspos("I043",1,"I018",1.5592,0.0038,3.15,idrotm[346],"ONLY");
  gMC->Gspos("I039",2,"I018",-1.8126,-0.485,0.0,idrotm[347],"ONLY");
  gMC->Gspos("I041",2,"I018",-1.7,-0.55,0.0,idrotm[347],"ONLY");
  gMC->Gspos("I040",2,"I018",-1.9204,-0.7118,0.0,idrotm[347],"ONLY");
  gMC->Gspos("I043",16,"I018",-1.5592,0.0038,-3.15,idrotm[347],"ONLY");
  gMC->Gspos("I042",9,"I018",-2.0342,-0.8189,-3.12,idrotm[347],"ONLY");
  gMC->Gspos("I043",15,"I018",-1.5592,0.0038,-2.25,idrotm[347],"ONLY");
  gMC->Gspos("I042",10,"I018",-2.0342,-0.8189,-2.22,idrotm[347],"ONLY");
  gMC->Gspos("I042",11,"I018",-2.0342,-0.8189,-1.32,idrotm[347],"ONLY");
  gMC->Gspos("I043",14,"I018",-1.5592,0.0038,-1.35,idrotm[347],"ONLY");
  gMC->Gspos("I042",12,"I018",-2.0342,-0.8189,-0.42,idrotm[347],"ONLY");
  gMC->Gspos("I043",13,"I018",-1.5592,0.0038,-0.45,idrotm[347],"ONLY");
  gMC->Gspos("I043",12,"I018",-1.5592,0.0038,0.45,idrotm[347],"ONLY");
  gMC->Gspos("I043",11,"I018",-1.5592,0.0038,1.35,idrotm[347],"ONLY");
  gMC->Gspos("I043",10,"I018",-1.5592,0.0038,2.25,idrotm[347],"ONLY");
  gMC->Gspos("I043",9,"I018",-1.5592,0.0038,3.15,idrotm[347],"ONLY");
  gMC->Gspos("I042",16,"I018",-2.0342,-0.8189,3.12,idrotm[347],"ONLY");
  gMC->Gspos("I042",15,"I018",-2.0342,-0.8189,2.28,idrotm[347],"ONLY");
  gMC->Gspos("I042",14,"I018",-2.0342,-0.8189,1.38,idrotm[347],"ONLY");
  gMC->Gspos("I042",13,"I018",-2.0342,-0.8189,0.48,idrotm[347],"ONLY");
  gMC->Gspos("I044",2,"I018",-2.7487,-1.3673,-0.2,idrotm[349],"ONLY");
  gMC->Gspos("I044",1,"I018",2.7487,-1.3673,-0.2,idrotm[351],"ONLY");
  gMC->Gspos("I303",1,"I302",1.8125,0.0,4.2605,idrotm[301],"ONLY");
  gMC->Gspos("I303",2,"I302",-1.8125,0.0,4.2605,idrotm[302],"ONLY");
  gMC->Gspos("I303",3,"I302",-1.8125,0.0,-4.2605,idrotm[303],"ONLY");
  gMC->Gspos("I303",4,"I302",1.8125,0.0,-4.2605,idrotm[304],"ONLY");
  gMC->Gspos("I403",1,"I402",1.8125,0.0,4.2605,idrotm[301],"ONLY");
  gMC->Gspos("I403",2,"I402",-1.8125,0.0,4.2605,idrotm[302],"ONLY");
  gMC->Gspos("I403",3,"I402",-1.8125,0.0,-4.2605,idrotm[303],"ONLY");
  gMC->Gspos("I403",4,"I402",1.8125,0.0,-4.2605,idrotm[304],"ONLY");
  gMC->Gspos("ITS3",1,"I302",0.0,0.0,0.0,0,"ONLY");
  gMC->Gspos("ITS4",1,"I402",0.0,0.0,0.0,0,"ONLY");
  
  // --- Place SSD volumes into their mother volume IT56  



  gMC->Gspos("I570",14,"IT56",-28.0681,-36.0619,-0.27,idrotm[566],"ONLY");
  gMC->Gspos("I570",15,"IT56",-21.677,-40.0556,-0.27,idrotm[567],"ONLY");
  gMC->Gspos("I570",16,"IT56",-14.838,-43.2217,-0.27,idrotm[568],"ONLY");
  gMC->Gspos("I570",17,"IT56",-7.4965,-44.9238,-0.27,idrotm[569],"ONLY");
  gMC->Gspos("I570",18,"IT56",0.0,-45.6977,-0.27,idrotm[533],"ONLY");
  gMC->Gspos("I570",19,"IT56",7.4965,-44.9238,-0.27,idrotm[534],"ONLY");
  gMC->Gspos("I570",20,"IT56",14.838,-43.2217,-0.27,idrotm[535],"ONLY");
  gMC->Gspos("I570",21,"IT56",21.677,-40.0556,-0.27,idrotm[623],"ONLY");
  gMC->Gspos("I570",22,"IT56",28.0681,-36.0619,-0.27,idrotm[537],"ONLY");
  gMC->Gspos("I570",23,"IT56",33.5085,-30.8468,-0.27,idrotm[538],"ONLY");
  gMC->Gspos("I570",24,"IT56",38.2566,-24.9943,-0.27,idrotm[539],"ONLY");
  gMC->Gspos("I570",25,"IT56",41.7089,-18.2952,-0.27,idrotm[540],"ONLY");
  gMC->Gspos("I570",26,"IT56",44.2994,-11.2181,-0.27,idrotm[541],"ONLY");
  gMC->Gspos("I570",27,"IT56",45.3894,-3.7611,-0.27,idrotm[542],"ONLY");
  gMC->Gspos("I570",28,"IT56",45.5416,3.7737,-0.27,idrotm[543],"ONLY");
  gMC->Gspos("I570",29,"IT56",44.1513,11.1806,-0.27,idrotm[544],"ONLY");
  gMC->Gspos("I570",30,"IT56",41.8487,18.3566,-0.27,idrotm[545],"ONLY");
  gMC->Gspos("I570",31,"IT56",38.1287,24.9107,-0.27,idrotm[546],"ONLY");
  gMC->Gspos("I570",32,"IT56",33.6209,30.9502,-0.27,idrotm[547],"ONLY");
  gMC->Gspos("I570",33,"IT56",27.9743,35.9414,-0.27,idrotm[548],"ONLY");
  gMC->Gspos("I570",34,"IT56",21.7497,40.1899,-0.27,idrotm[549],"ONLY");
  gMC->Gspos("I570",35,"IT56",14.7884,43.0772,-0.27,idrotm[550],"ONLY");
  gMC->Gspos("I570",36,"IT56",7.5216,45.0744,-0.27,idrotm[551],"ONLY");
  gMC->Gspos("I570",37,"IT56",0.0,45.545,-0.27,0,"ONLY");
  gMC->Gspos("I570",38,"IT56",-7.5216,45.0744,-0.27,idrotm[552],"ONLY");
  gMC->Gspos("I570",1,"IT56",-14.7884,43.0772,-0.27,idrotm[553],"ONLY");
  gMC->Gspos("I570",2,"IT56",-21.7497,40.1899,-0.27,idrotm[620],"ONLY");
  gMC->Gspos("I570",3,"IT56",-27.9743,35.9414,-0.27,idrotm[555],"ONLY");
  gMC->Gspos("I570",4,"IT56",-33.6209,30.9502,-0.27,idrotm[556],"ONLY");
  gMC->Gspos("I570",5,"IT56",-38.1287,24.9108,-0.27,idrotm[557],"ONLY");
  gMC->Gspos("I570",6,"IT56",-41.8487,18.3566,-0.27,idrotm[558],"ONLY");
  gMC->Gspos("I570",7,"IT56",-44.1513,11.1806,-0.27,idrotm[559],"ONLY");
  gMC->Gspos("I570",8,"IT56",-45.5416,3.7737,-0.27,idrotm[560],"ONLY");
  gMC->Gspos("I570",9,"IT56",-45.3894,-3.7611,-0.27,idrotm[561],"ONLY");
  gMC->Gspos("I570",10,"IT56",-44.2994,-11.2181,-0.27,idrotm[562],"ONLY");
  gMC->Gspos("I570",11,"IT56",-41.7089,-18.2952,-0.27,idrotm[563],"ONLY");
  gMC->Gspos("I570",12,"IT56",-38.2566,-24.9943,-0.27,idrotm[564],"ONLY");
  gMC->Gspos("I570",13,"IT56",-33.5086,-30.8468,-0.27,idrotm[565],"ONLY");
  gMC->Gspos("I569",8,"IT56",-43.5484,3.6085,0.0,idrotm[560],"ONLY");
  gMC->Gspos("I569",9,"IT56",-43.3963,-3.5959,0.0,idrotm[561],"ONLY");
  gMC->Gspos("I569",10,"IT56",-42.3606,-10.7271,0.0,idrotm[562],"ONLY");
  gMC->Gspos("I569",11,"IT56",-39.8773,-17.4918,0.0,idrotm[563],"ONLY");
  gMC->Gspos("I569",12,"IT56",-36.5823,-23.9004,0.0,idrotm[564],"ONLY");
  gMC->Gspos("I569",13,"IT56",-32.0371,-29.4922,0.0,idrotm[565],"ONLY");
  gMC->Gspos("I569",14,"IT56",-26.8397,-34.4836,0.0,idrotm[566],"ONLY");
  gMC->Gspos("I569",15,"IT56",-20.7251,-38.2967,0.0,idrotm[567],"ONLY");
  gMC->Gspos("I569",16,"IT56",-14.1886,-41.33,0.0,idrotm[568],"ONLY");
  gMC->Gspos("I569",17,"IT56",-7.1673,-42.9511,0.0,idrotm[569],"ONLY");
  gMC->Gspos("I569",18,"IT56",0.0,-43.6977,0.0,idrotm[533],"ONLY");
  gMC->Gspos("I569",19,"IT56",7.1673,-42.9511,0.0,idrotm[534],"ONLY");
  gMC->Gspos("I569",20,"IT56",14.1886,-41.33,0.0,idrotm[535],"ONLY");
  gMC->Gspos("I569",21,"IT56",20.7251,-38.2967,0.0,idrotm[623],"ONLY");
  gMC->Gspos("I569",22,"IT56",26.8397,-34.4836,0.0,idrotm[537],"ONLY");
  gMC->Gspos("I569",23,"IT56",32.0371,-29.4922,0.0,idrotm[538],"ONLY");
  gMC->Gspos("I569",24,"IT56",36.5822,-23.9004,0.0,idrotm[539],"ONLY");
  gMC->Gspos("I569",25,"IT56",39.8773,-17.4918,0.0,idrotm[540],"ONLY");
  gMC->Gspos("I569",26,"IT56",42.3606,-10.7272,0.0,idrotm[541],"ONLY");
  gMC->Gspos("I569",27,"IT56",43.3963,-3.5959,0.0,idrotm[542],"ONLY");
  gMC->Gspos("I569",28,"IT56",43.5484,3.6085,0.0,idrotm[543],"ONLY");
  gMC->Gspos("I569",29,"IT56",42.2125,10.6897,0.0,idrotm[544],"ONLY");
  gMC->Gspos("I569",30,"IT56",40.0172,17.5532,0.0,idrotm[545],"ONLY");
  gMC->Gspos("I569",31,"IT56",36.4544,23.8169,0.0,idrotm[546],"ONLY");
  gMC->Gspos("I569",32,"IT56",32.1494,29.5956,0.0,idrotm[547],"ONLY");
  gMC->Gspos("I569",33,"IT56",26.7459,34.3631,0.0,idrotm[548],"ONLY");
  gMC->Gspos("I569",34,"IT56",20.7978,38.431,0.0,idrotm[549],"ONLY");
  gMC->Gspos("I569",35,"IT56",14.139,41.1856,0.0,idrotm[550],"ONLY");
  gMC->Gspos("I569",36,"IT56",7.1924,43.1017,0.0,idrotm[551],"ONLY");
  gMC->Gspos("I569",37,"IT56",0.0,43.545,0.0,0,"ONLY");
  gMC->Gspos("I569",38,"IT56",-7.1924,43.1017,0.0,idrotm[552],"ONLY");
  gMC->Gspos("I569",1,"IT56",-14.139,41.1856,0.0,idrotm[553],"ONLY");
  gMC->Gspos("I569",2,"IT56",-20.7978,38.431,0.0,idrotm[620],"ONLY");
  gMC->Gspos("I569",3,"IT56",-26.7459,34.3631,0.0,idrotm[555],"ONLY");
  gMC->Gspos("I569",4,"IT56",-32.1494,29.5956,0.0,idrotm[556],"ONLY");
  gMC->Gspos("I569",5,"IT56",-36.4544,23.8169,0.0,idrotm[557],"ONLY");
  gMC->Gspos("I569",6,"IT56",-40.0172,17.5532,0.0,idrotm[558],"ONLY");
  gMC->Gspos("I569",7,"IT56",-42.2125,10.6897,0.0,idrotm[559],"ONLY");
  gMC->Gspos("I571",15,"IT56",-21.2916,-34.387,0.0,idrotm[501],"ONLY");
  gMC->Gspos("I571",14,"IT56",-27.351,-30.0026,0.0,idrotm[503],"ONLY");
  gMC->Gspos("I571",13,"IT56",-32.2758,-24.3735,0.0,idrotm[504],"ONLY");
  gMC->Gspos("I571",12,"IT56",-36.3422,-18.0963,0.0,idrotm[505],"ONLY");
  gMC->Gspos("I571",11,"IT56",-38.901,-11.0683,0.0,idrotm[506],"ONLY");
  gMC->Gspos("I571",10,"IT56",-40.4252,-3.7459,0.0,idrotm[507],"ONLY");
  gMC->Gspos("I571",9,"IT56",-40.2725,3.7318,0.0,idrotm[508],"ONLY");
  gMC->Gspos("I571",8,"IT56",-39.0486,11.1103,0.0,idrotm[509],"ONLY");
  gMC->Gspos("I571",7,"IT56",-36.2049,18.0279,0.0,idrotm[510],"ONLY");
  gMC->Gspos("I571",6,"IT56",-32.3982,24.466,0.0,idrotm[511],"ONLY");
  gMC->Gspos("I571",5,"IT56",-27.2476,29.8892,0.0,idrotm[512],"ONLY");
  gMC->Gspos("I571",4,"IT56",-21.3723,34.5175,0.0,idrotm[513],"ONLY");
  gMC->Gspos("I571",3,"IT56",-14.6104,37.7138,0.0,idrotm[653],"ONLY");
  gMC->Gspos("I571",2,"IT56",-7.4599,39.9072,0.0,idrotm[514],"ONLY");
  gMC->Gspos("I571",1,"IT56",0.0,40.445,0.0,0,"ONLY");
  gMC->Gspos("I571",34,"IT56",7.46,39.9071,0.0,idrotm[515],"ONLY");
  gMC->Gspos("I571",33,"IT56",14.6104,37.7138,0.0,idrotm[516],"ONLY");
  gMC->Gspos("I571",32,"IT56",21.3723,34.5175,0.0,idrotm[517],"ONLY");
  gMC->Gspos("I571",31,"IT56",27.2476,29.8892,0.0,idrotm[518],"ONLY");
  gMC->Gspos("I571",30,"IT56",32.3983,24.466,0.0,idrotm[519],"ONLY");
  gMC->Gspos("I571",29,"IT56",36.2049,18.0279,0.0,idrotm[520],"ONLY");
  gMC->Gspos("I571",28,"IT56",39.0486,11.1103,0.0,idrotm[521],"ONLY");
  gMC->Gspos("I571",27,"IT56",40.2725,3.7318,0.0,idrotm[522],"ONLY");
  gMC->Gspos("I571",26,"IT56",40.4252,-3.746,0.0,idrotm[523],"ONLY");
  gMC->Gspos("I571",25,"IT56",38.901,-11.0683,0.0,idrotm[524],"ONLY");
  gMC->Gspos("I571",24,"IT56",36.3422,-18.0963,0.0,idrotm[525],"ONLY");
  gMC->Gspos("I571",23,"IT56",32.2758,-24.3736,0.0,idrotm[526],"ONLY");
  gMC->Gspos("I571",22,"IT56",27.351,-30.0026,0.0,idrotm[527],"ONLY");
  gMC->Gspos("I571",21,"IT56",21.2915,-34.387,0.0,idrotm[528],"ONLY");
  gMC->Gspos("I571",20,"IT56",14.6658,-37.8569,0.0,idrotm[618],"ONLY");
  gMC->Gspos("I571",19,"IT56",7.4317,-39.7563,0.0,idrotm[529],"ONLY");
  gMC->Gspos("I571",18,"IT56",0.0,-40.5984,0.0,idrotm[533],"ONLY");
  gMC->Gspos("I571",17,"IT56",-7.4318,-39.7563,0.0,idrotm[530],"ONLY");
  gMC->Gspos("I571",16,"IT56",-14.6659,-37.8569,0.0,idrotm[531],"ONLY");
  gMC->Gspos("I565",13,"IT56",-30.6798,-23.1683,0.0,idrotm[504],"ONLY");
  gMC->Gspos("I565",12,"IT56",-34.5519,-17.2048,0.0,idrotm[505],"ONLY");
  gMC->Gspos("I565",11,"IT56",-36.9774,-10.521,0.0,idrotm[506],"ONLY");
  gMC->Gspos("I565",10,"IT56",-38.4338,-3.5614,0.0,idrotm[507],"ONLY");
  gMC->Gspos("I565",9,"IT56",-38.281,3.5473,0.0,idrotm[508],"ONLY");
  gMC->Gspos("I565",8,"IT56",-37.1249,10.563,0.0,idrotm[509],"ONLY");
  gMC->Gspos("I565",7,"IT56",-34.4146,17.1364,0.0,idrotm[510],"ONLY");
  gMC->Gspos("I565",6,"IT56",-30.8022,23.2608,0.0,idrotm[511],"ONLY");
  gMC->Gspos("I565",5,"IT56",-25.9002,28.4112,0.0,idrotm[512],"ONLY");
  gMC->Gspos("I565",4,"IT56",-20.3195,32.817,0.0,idrotm[513],"ONLY");
  gMC->Gspos("I565",3,"IT56",-13.8879,35.8489,0.0,idrotm[653],"ONLY");
  gMC->Gspos("I565",2,"IT56",-7.0924,37.9412,0.0,idrotm[514],"ONLY");
  gMC->Gspos("I565",1,"IT56",0.0,38.445,0.0,0,"ONLY");
  gMC->Gspos("I565",34,"IT56",7.0925,37.9412,0.0,idrotm[515],"ONLY");
  gMC->Gspos("I565",33,"IT56",13.888,35.8489,0.0,idrotm[516],"ONLY");
  gMC->Gspos("I565",32,"IT56",20.3195,32.817,0.0,idrotm[517],"ONLY");
  gMC->Gspos("I565",31,"IT56",25.9002,28.4112,0.0,idrotm[518],"ONLY");
  gMC->Gspos("I565",30,"IT56",30.8022,23.2607,0.0,idrotm[519],"ONLY");
  gMC->Gspos("I565",29,"IT56",34.4146,17.1364,0.0,idrotm[520],"ONLY");
  gMC->Gspos("I565",28,"IT56",37.125,10.5629,0.0,idrotm[521],"ONLY");
  gMC->Gspos("I565",27,"IT56",38.281,3.5472,0.0,idrotm[522],"ONLY");
  gMC->Gspos("I565",26,"IT56",38.4338,-3.5614,0.0,idrotm[523],"ONLY");
  gMC->Gspos("I565",25,"IT56",36.9774,-10.521,0.0,idrotm[524],"ONLY");
  gMC->Gspos("I565",24,"IT56",34.5519,-17.2048,0.0,idrotm[525],"ONLY");
  gMC->Gspos("I565",23,"IT56",30.6798,-23.1683,0.0,idrotm[526],"ONLY");
  gMC->Gspos("I565",22,"IT56",26.0036,-28.5246,0.0,idrotm[527],"ONLY");
  gMC->Gspos("I565",21,"IT56",20.2387,-32.6866,0.0,idrotm[528],"ONLY");
  gMC->Gspos("I565",20,"IT56",13.9433,-35.992,0.0,idrotm[618],"ONLY");
  gMC->Gspos("I565",19,"IT56",7.0642,-37.7904,0.0,idrotm[529],"ONLY");
  gMC->Gspos("I565",18,"IT56",0.0,-38.5984,0.0,idrotm[533],"ONLY");
  gMC->Gspos("I565",17,"IT56",-7.0643,-37.7904,0.0,idrotm[530],"ONLY");
  gMC->Gspos("I565",16,"IT56",-13.9434,-35.992,0.0,idrotm[531],"ONLY");
  gMC->Gspos("I565",15,"IT56",-20.2387,-32.6866,0.0,idrotm[501],"ONLY");
  gMC->Gspos("I565",14,"IT56",-26.0036,-28.5246,0.0,idrotm[503],"ONLY");
  gMC->Gspos("I553",1,"I570",0.005,0.0,53.98,0,"ONLY");
  gMC->Gspos("I523",2,"I570",0.0,0.0,48.875,0,"ONLY");
  gMC->Gspos("I523",3,"I570",0.0,0.0,44.965,0,"ONLY");
  gMC->Gspos("I523",4,"I570",0.0,0.0,41.055,0,"ONLY");
  gMC->Gspos("I523",5,"I570",0.0,0.0,33.235,0,"ONLY");
  gMC->Gspos("I523",6,"I570",0.0,0.0,37.145,0,"ONLY");
  gMC->Gspos("I523",7,"I570",0.0,0.0,29.325,0,"ONLY");
  gMC->Gspos("I523",8,"I570",0.0,0.0,25.415,0,"ONLY");
  gMC->Gspos("I523",9,"I570",0.0,0.0,21.505,0,"ONLY");
  gMC->Gspos("I523",10,"I570",0.0,0.0,13.685,0,"ONLY");
  gMC->Gspos("I523",11,"I570",0.0,0.0,17.595,0,"ONLY");
  gMC->Gspos("I523",12,"I570",0.0,0.0,9.775,0,"ONLY");
  gMC->Gspos("I523",13,"I570",0.0,0.0,5.865,0,"ONLY");
  gMC->Gspos("I523",14,"I570",0.0,0.0,1.955,0,"ONLY");
  gMC->Gspos("I523",15,"I570",0.0,0.0,-1.955,0,"ONLY");
  gMC->Gspos("I523",16,"I570",0.0,0.0,-9.775,0,"ONLY");
  gMC->Gspos("I523",17,"I570",0.0,0.0,-5.865,0,"ONLY");
  gMC->Gspos("I523",18,"I570",0.0,0.0,-13.685,0,"ONLY");
  gMC->Gspos("I523",19,"I570",0.0,0.0,-21.505,0,"ONLY");
  gMC->Gspos("I523",20,"I570",0.0,0.0,-17.595,0,"ONLY");
  gMC->Gspos("I523",21,"I570",0.0,0.0,-25.415,0,"ONLY");
  gMC->Gspos("I523",22,"I570",0.0,0.0,-29.325,0,"ONLY");
  gMC->Gspos("I523",23,"I570",0.0,0.0,-37.145,0,"ONLY");
  gMC->Gspos("I523",24,"I570",0.0,0.0,-33.235,0,"ONLY");
  gMC->Gspos("I523",25,"I570",0.0,0.0,-44.965,0,"ONLY");
  gMC->Gspos("I523",26,"I570",0.0,0.0,-41.055,0,"ONLY");
  gMC->Gspos("I553",2,"I570",-0.005,0.0,-53.98,idrotm[570],"ONLY");
  gMC->Gspos("I523",1,"I570",0.0,0.0,-48.875,0,"ONLY");
  gMC->Gspos("I566",1,"I569",0.0,-0.03,46.9203,idrotm[532],"ONLY");
  gMC->Gspos("I566",2,"I569",0.0,0.03,43.0103,0,"ONLY");
  gMC->Gspos("I566",3,"I569",0.0,-0.03,39.1003,idrotm[532],"ONLY");
  gMC->Gspos("I566",4,"I569",0.0,0.03,35.1903,0,"ONLY");
  gMC->Gspos("I566",5,"I569",0.0,-0.03,31.2803,idrotm[532],"ONLY");
  gMC->Gspos("I566",6,"I569",0.0,0.03,27.3703,0,"ONLY");
  gMC->Gspos("I566",7,"I569",0.0,-0.03,23.4603,idrotm[532],"ONLY");
  gMC->Gspos("I566",8,"I569",0.0,0.03,19.5503,0,"ONLY");
  gMC->Gspos("I566",9,"I569",0.0,-0.03,15.6403,idrotm[532],"ONLY");
  gMC->Gspos("I566",10,"I569",0.0,0.03,11.7303,0,"ONLY");
  gMC->Gspos("I566",11,"I569",0.0,-0.03,7.8203,idrotm[532],"ONLY");
  gMC->Gspos("I566",12,"I569",0.0,0.03,3.9103,0,"ONLY");
  gMC->Gspos("I566",13,"I569",0.0,-0.03,0.0003,0,"ONLY");
  gMC->Gspos("I566",14,"I569",0.0,0.03,-3.9097,0,"ONLY");
  gMC->Gspos("I566",15,"I569",0.0,-0.03,-7.8197,idrotm[532],"ONLY");
  gMC->Gspos("I566",16,"I569",0.0,0.03,-11.7297,0,"ONLY");
  gMC->Gspos("I566",17,"I569",0.0,-0.03,-15.6397,0,"ONLY");
  gMC->Gspos("I566",18,"I569",0.0,0.03,-19.5497,0,"ONLY");
  gMC->Gspos("I566",19,"I569",0.0,-0.03,-23.4597,idrotm[532],"ONLY");
  gMC->Gspos("I566",20,"I569",0.0,0.03,-27.3697,0,"ONLY");
  gMC->Gspos("I566",21,"I569",0.0,-0.03,-31.2797,idrotm[532],"ONLY");
  gMC->Gspos("I566",22,"I569",0.0,0.03,-35.1897,0,"ONLY");
  gMC->Gspos("I566",23,"I569",0.0,-0.03,-39.0997,0,"ONLY");
  gMC->Gspos("I566",24,"I569",0.0,0.03,-43.0097,0,"ONLY");
  gMC->Gspos("I566",25,"I569",0.0,-0.03,-46.9197,idrotm[532],"ONLY");
  gMC->Gspos("I544",1,"I571",0.0101,0.0,48.115,0,"ONLY");
  gMC->Gspos("I516",23,"I571",0.0001,0.0,43.01,0,"ONLY");
  gMC->Gspos("I516",22,"I571",0.0001,0.0,39.1,0,"ONLY");
  gMC->Gspos("I516",21,"I571",0.0001,0.0,35.19,0,"ONLY");
  gMC->Gspos("I516",20,"I571",0.0001,0.0,31.28,0,"ONLY");
  gMC->Gspos("I516",19,"I571",0.0001,0.0,27.37,0,"ONLY");
  gMC->Gspos("I516",18,"I571",0.0001,0.0,23.46,0,"ONLY");
  gMC->Gspos("I516",17,"I571",0.0001,0.0,19.55,0,"ONLY");
  gMC->Gspos("I516",16,"I571",0.0001,0.0,15.64,0,"ONLY");
  gMC->Gspos("I516",15,"I571",0.0001,0.0,11.73,0,"ONLY");
  gMC->Gspos("I516",14,"I571",0.0001,0.0,7.82,0,"ONLY");
  gMC->Gspos("I516",13,"I571",0.0001,0.0,3.91,0,"ONLY");
  gMC->Gspos("I516",12,"I571",0.0001,0.0,0.0,0,"ONLY");
  gMC->Gspos("I516",11,"I571",0.0001,0.0,-3.91,0,"ONLY");
  gMC->Gspos("I516",10,"I571",0.0001,0.0,-7.82,0,"ONLY");
  gMC->Gspos("I516",9,"I571",0.0001,0.0,-11.73,0,"ONLY");
  gMC->Gspos("I516",8,"I571",0.0001,0.0,-15.64,0,"ONLY");
  gMC->Gspos("I516",7,"I571",0.0001,0.0,-19.55,0,"ONLY");
  gMC->Gspos("I516",6,"I571",0.0001,0.0,-23.46,0,"ONLY");
  gMC->Gspos("I516",5,"I571",0.0001,0.0,-27.37,0,"ONLY");
  gMC->Gspos("I516",4,"I571",0.0001,0.0,-31.28,0,"ONLY");
  gMC->Gspos("I516",3,"I571",0.0001,0.0,-35.19,0,"ONLY");
  gMC->Gspos("I516",2,"I571",0.0001,0.0,-39.1,0,"ONLY");
  gMC->Gspos("I516",1,"I571",0.0001,0.0,-43.01,0,"ONLY");
  gMC->Gspos("I544",2,"I571",-0.0099,0.0,-48.115,idrotm[570],"ONLY");
  gMC->Gspos("I562",1,"I565",0.0,0.03,41.1546,0,"ONLY");
  gMC->Gspos("I562",2,"I565",0.0,-0.03,37.2246,0,"ONLY");
  gMC->Gspos("I562",3,"I565",0.0,0.03,33.3146,0,"ONLY");
  gMC->Gspos("I562",4,"I565",0.0,-0.03,29.3846,0,"ONLY");
  gMC->Gspos("I562",5,"I565",0.0,0.03,25.4746,0,"ONLY");
  gMC->Gspos("I562",6,"I565",0.0,-0.03,21.5446,0,"ONLY");
  gMC->Gspos("I562",7,"I565",0.0,0.03,17.6346,0,"ONLY");
  gMC->Gspos("I562",8,"I565",0.0,-0.03,13.7046,0,"ONLY");
  gMC->Gspos("I562",9,"I565",0.0,0.03,9.7946,0,"ONLY");
  gMC->Gspos("I562",10,"I565",0.0,-0.03,5.8645,0,"ONLY");
  gMC->Gspos("I562",11,"I565",0.0,0.03,1.9546,0,"ONLY");
  gMC->Gspos("I562",12,"I565",0.0,-0.03,-1.9754,0,"ONLY");
  gMC->Gspos("I562",13,"I565",0.0,0.03,-5.8855,0,"ONLY");
  gMC->Gspos("I562",14,"I565",0.0,-0.03,-9.8154,0,"ONLY");
  gMC->Gspos("I562",15,"I565",0.0,0.03,-13.7254,0,"ONLY");
  gMC->Gspos("I562",16,"I565",0.0,-0.03,-17.6555,0,"ONLY");
  gMC->Gspos("I562",17,"I565",0.0,0.03,-21.5655,0,"ONLY");
  gMC->Gspos("I562",18,"I565",0.0,-0.03,-25.4954,0,"ONLY");
  gMC->Gspos("I562",19,"I565",0.0,0.03,-29.4054,0,"ONLY");
  gMC->Gspos("I562",20,"I565",0.0,-0.03,-33.3354,0,"ONLY");
  gMC->Gspos("I562",21,"I565",0.0,0.03,-37.2454,0,"ONLY");
  gMC->Gspos("I562",22,"I565",0.0,-0.03,-41.1554,0,"ONLY");
  gMC->Gspos("I559",1,"I553",2.25,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I560",1,"I553",2.25,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I560",2,"I553",-2.25,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I558",1,"I553",-1.7167,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I557",1,"I553",-1.8533,-1.341,0.0,idrotm[581],"ONLY");
  gMC->Gspos("I558",2,"I553",1.8367,-1.3122,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I557",2,"I553",1.75,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I558",3,"I553",-0.12,1.6613,0.0,idrotm[581],"ONLY");
  gMC->Gspos("I557",3,"I553",0.1034,1.6901,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I556",3,"I553",-1.031,0.2033,-2.203,idrotm[580],"ONLY");
  gMC->Gspos("I556",1,"I553",1.0311,0.2033,-0.287,idrotm[576],"ONLY");
  gMC->Gspos("I554",1,"I553",0.0,-1.58,0.71,0,"ONLY");
  gMC->Gspos("I555",1,"I553",-0.0072,-1.58,-1.2311,idrotm[633],"ONLY");
  gMC->Gspos("I556",2,"I553",1.0311,0.2033,-2.203,idrotm[577],"ONLY");
  gMC->Gspos("I556",4,"I553",-1.031,0.2033,-0.287,idrotm[579],"ONLY");
  gMC->Gspos("I559",2,"I553",-2.25,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I561",1,"I553",2.1,-1.615,-0.24,0,"MANY");
  gMC->Gspos("I561",2,"I553",-2.1,-1.615,-0.24,idrotm[573],"MANY");
  gMC->Gspos("I519",37,"I523",0.0001,-1.79,-0.99,idrotm[586],"ONLY");
  gMC->Gspos("I519",36,"I523",-3.2986,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",35,"I523",-3.2986,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",34,"I523",-3.2286,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",33,"I523",-3.2286,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",32,"I523",-3.1586,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",31,"I523",-3.1586,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",30,"I523",-1.3436,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",29,"I523",-1.3436,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",28,"I523",-1.2736,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",27,"I523",-1.2736,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",26,"I523",-1.2036,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",25,"I523",-1.2036,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",24,"I523",-1.0458,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",23,"I523",-1.0458,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",22,"I523",-0.9758,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",21,"I523",-0.9758,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",20,"I523",-0.9058,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",19,"I523",-0.9058,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",18,"I523",0.9092,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",17,"I523",0.9092,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",16,"I523",0.9792,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",15,"I523",0.9792,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",14,"I523",1.0492,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",13,"I523",1.0492,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",12,"I523",1.207,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",11,"I523",1.207,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",10,"I523",1.277,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",9,"I523",1.277,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",8,"I523",1.347,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",7,"I523",1.347,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",6,"I523",3.162,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I519",5,"I523",3.162,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",4,"I523",3.232,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I519",3,"I523",3.232,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I521",12,"I523",-2.8209,-1.7925,-0.982,0,"ONLY");
  gMC->Gspos("I521",11,"I523",-1.6895,-1.7925,-0.982,0,"ONLY");
  gMC->Gspos("I521",10,"I523",-0.5631,-1.7925,-0.982,0,"ONLY");
  gMC->Gspos("I521",9,"I523",0.5633,-1.7925,-0.982,0,"ONLY");
  gMC->Gspos("I521",8,"I523",1.6861,-1.7925,-0.982,0,"ONLY");
  gMC->Gspos("I521",7,"I523",2.8161,-1.7925,-0.982,0,"ONLY");
  gMC->Gspos("I519",2,"I523",3.302,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I520",3,"I523",0.0001,-1.845,-1.19,0,"ONLY");
  gMC->Gspos("I520",2,"I523",-2.2499,-1.845,-1.19,0,"ONLY");
  gMC->Gspos("I521",6,"I523",-2.8209,-1.7075,-0.982,0,"ONLY");
  gMC->Gspos("I521",5,"I523",-1.6895,-1.7075,-0.982,0,"ONLY");
  gMC->Gspos("I521",4,"I523",-0.5631,-1.7075,-0.982,0,"ONLY");
  gMC->Gspos("I521",3,"I523",0.5633,-1.7075,-0.982,0,"ONLY");
  gMC->Gspos("I521",2,"I523",1.6861,-1.7075,-0.982,0,"ONLY");
  gMC->Gspos("I518",1,"I523",0.0001,-1.75,-1.065,0,"ONLY");
  gMC->Gspos("I519",1,"I523",3.302,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I520",1,"I523",2.2501,-1.845,-1.19,0,"ONLY");
  gMC->Gspos("I521",1,"I523",2.8161,-1.7075,-0.982,0,"ONLY");
  gMC->Gspos("I522",1,"I523",2.2501,-1.655,-1.3,idrotm[583],"MANY");
  gMC->Gspos("I522",2,"I523",-2.2499,-1.655,-1.3,idrotm[583],"MANY");
  gMC->Gspos("I542",2,"I523",-2.2499,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I541",2,"I523",-2.2499,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I541",1,"I523",2.2501,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I542",1,"I523",2.2501,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I543",1,"I523",2.1001,-1.615,0.955,0,"MANY");
  gMC->Gspos("I543",2,"I523",-2.0999,-1.615,0.955,idrotm[573],"MANY");
  gMC->Gspos("I537",2,"I523",1.7501,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I538",2,"I523",1.8368,-1.3122,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I537",3,"I523",0.1035,1.6901,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I538",3,"I523",-0.1199,1.6612,0.0,idrotm[581],"ONLY");
  gMC->Gspos("I538",1,"I523",-1.7166,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I537",1,"I523",-1.8532,-1.341,0.0,idrotm[581],"ONLY");
  gMC->Gspos("I536",3,"I523",-1.031,0.2033,-1.008,idrotm[580],"ONLY");
  gMC->Gspos("I536",4,"I523",-1.031,0.2033,0.908,idrotm[579],"ONLY");
  gMC->Gspos("I535",1,"I523",-0.0072,-1.58,-0.0361,idrotm[633],"ONLY");
  gMC->Gspos("I536",2,"I523",1.0312,0.2033,-1.008,idrotm[577],"ONLY");
  gMC->Gspos("I536",1,"I523",1.0312,0.2033,0.908,idrotm[576],"ONLY");
  gMC->Gspos("I534",1,"I523",0.0001,-1.58,1.905,0,"ONLY");
  gMC->Gspos("I540",1,"I523",0.0001,-1.785,1.905,idrotm[571],"ONLY");
  gMC->Gspos("I539",1,"I523",1.8001,-1.75,-0.195,idrotm[571],"ONLY");
  gMC->Gspos("I539",2,"I523",-1.7999,-1.75,-0.195,idrotm[572],"ONLY");
  gMC->Gspos("ITS6",1,"I566",0.0,0.0,0.0,0,"ONLY");
  gMC->Gspos("I550",1,"I544",2.25,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I551",1,"I544",2.25,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I551",2,"I544",-2.25,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I550",2,"I544",-2.25,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I549",1,"I544",1.7167,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I548",1,"I544",1.8533,-1.341,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I547",1,"I544",1.0311,0.2033,-0.287,idrotm[576],"ONLY");
  gMC->Gspos("I545",1,"I544",0.0,-1.58,0.71,0,"ONLY");
  gMC->Gspos("I547",2,"I544",1.0311,0.2033,-2.203,idrotm[577],"ONLY");
  gMC->Gspos("I546",1,"I544",-0.0073,-1.58,-1.2311,idrotm[633],"ONLY");
  gMC->Gspos("I547",4,"I544",-1.0311,0.2033,-0.287,idrotm[579],"ONLY");
  gMC->Gspos("I547",3,"I544",-1.0311,0.2033,-2.203,idrotm[580],"ONLY");
  gMC->Gspos("I548",2,"I544",-0.1033,1.6901,0.0,idrotm[581],"O]NLY");
  gMC->Gspos("I549",2,"I544",0.12,1.6613,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I549",3,"I544",-1.8367,-1.3122,0.0,idrotm[581],"ONLY");
  gMC->Gspos("I548",3,"I544",-1.75,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I552",1,"I544",2.1,-1.615,-0.24,0,"MANY");
  gMC->Gspos("I552",2,"I544",-2.1,-1.615,-0.24,idrotm[573],"MANY");
  gMC->Gspos("I515",12,"I516",-1.6896,-1.7075,-0.9822,0,"ONLY");
  gMC->Gspos("I515",11,"I516",-1.6896,-1.7925,-0.9822,0,"ONLY");
  gMC->Gspos("I513",37,"I516",0.0,-1.79,-1.035,idrotm[586],"ONLY");
  gMC->Gspos("I513",1,"I516",-3.2987,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I515",1,"I516",-2.816,-1.7075,-0.9822,0,"ONLY");
  gMC->Gspos("I514",1,"I516",-2.25,-1.845,-1.19,0,"ONLY");
  gMC->Gspos("I514",2,"I516",0.0,-1.845,-1.19,0,"ONLY");
  gMC->Gspos("I514",3,"I516",2.25,-1.845,-1.19,0,"ONLY");
  gMC->Gspos("I515",2,"I516",-2.816,-1.7925,-0.9822,0,"ONLY");
  gMC->Gspos("I513",2,"I516",-3.2987,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I515",3,"I516",-0.5632,-1.7075,-0.9822,0,"ONLY");
  gMC->Gspos("I515",4,"I516",-0.5632,-1.7925,-0.9822,0,"ONLY");
  gMC->Gspos("I515",5,"I516",0.5632,-1.7925,-0.9822,0,"ONLY");
  gMC->Gspos("I515",6,"I516",0.5632,-1.7075,-0.9822,0,"ONLY");
  gMC->Gspos("I515",7,"I516",1.6896,-1.7075,-0.9822,0,"ONLY");
  gMC->Gspos("I515",8,"I516",1.6896,-1.7925,-0.9822,0,"ONLY");
  gMC->Gspos("I515",9,"I516",2.816,-1.7925,-0.9822,0,"ONLY");
  gMC->Gspos("I515",10,"I516",2.816,-1.7075,-0.9822,0,"ONLY");
  gMC->Gspos("I513",3,"I516",-3.2287,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",4,"I516",-3.2287,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",5,"I516",-3.1587,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",6,"I516",-3.1587,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",7,"I516",-1.3437,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",8,"I516",-1.3437,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",9,"I516",-1.2737,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",10,"I516",-1.2737,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",11,"I516",-1.2037,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",12,"I516",-1.2037,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",13,"I516",-1.046,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",14,"I516",-1.046,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",15,"I516",-0.976,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",16,"I516",-0.976,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",17,"I516",-0.906,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",18,"I516",-0.906,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",19,"I516",0.9091,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",20,"I516",0.9091,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",21,"I516",0.9791,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",22,"I516",0.9791,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",23,"I516",1.0491,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",24,"I516",1.0491,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",25,"I516",1.2068,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",26,"I516",1.2068,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",27,"I516",1.2768,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",28,"I516",1.2768,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",29,"I516",1.3469,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",30,"I516",1.3469,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",31,"I516",3.1619,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",32,"I516",3.1619,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",33,"I516",3.2319,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I513",34,"I516",3.2319,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",35,"I516",3.3019,-1.79,-1.2943,0,"ONLY");
  gMC->Gspos("I513",36,"I516",3.3019,-1.71,-1.2943,0,"ONLY");
  gMC->Gspos("I512",1,"I516",0.0,-1.75,-1.065,0,"ONLY");
  gMC->Gspos("I528",1,"I516",1.7167,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I527",1,"I516",1.8534,-1.341,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I528",2,"I516",0.12,1.6613,0.0,idrotm[575],"ONLY");
  gMC->Gspos("I527",2,"I516",-0.1033,1.6901,0.0,idrotm[581],"ONLY");
  gMC->Gspos("I527",3,"I516",-1.75,-1.52,0.0,idrotm[583],"ONLY");
  gMC->Gspos("I528",3,"I516",-1.8367,-1.3122,0.0,idrotm[581],"ONLY");
  gMC->Gspos("I526",2,"I516",1.0311,0.2033,-1.008,idrotm[577],"ONLY");
  gMC->Gspos("I525",1,"I516",-0.0073,-1.58,-0.0361,idrotm[633],"ONLY");
  gMC->Gspos("I524",1,"I516",0.0,-1.58,1.905,0,"ONLY");
  gMC->Gspos("I526",1,"I516",1.0311,0.2033,0.908,idrotm[576],"ONLY");
  gMC->Gspos("I526",3,"I516",-1.0311,0.2033,0.908,idrotm[579],"ONLY");
  gMC->Gspos("I526",4,"I516",-1.0311,0.2033,-1.008,idrotm[580],"ONLY");
  gMC->Gspos("I529",1,"I516",1.8,-1.75,-0.195,idrotm[571],"ONLY");
  gMC->Gspos("I530",1,"I516",0.0,-1.785,1.905,idrotm[571],"ONLY");
  gMC->Gspos("I529",2,"I516",-1.8,-1.75,-0.195,idrotm[572],"ONLY");
  gMC->Gspos("I517",1,"I516",2.25,-1.655,-1.3,idrotm[583],"MANY");
  gMC->Gspos("I517",2,"I516",-2.25,-1.655,-1.3,idrotm[584],"MANY");
  gMC->Gspos("I531",2,"I516",-2.25,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I531",1,"I516",2.25,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I532",1,"I516",2.25,-1.615,0.0,0,"ONLY");
  gMC->Gspos("I532",2,"I516",-2.25,-1.615,0.0,idrotm[573],"ONLY");
  gMC->Gspos("I533",1,"I516",2.1,-1.615,0.955,0,"MANY");
  gMC->Gspos("I533",2,"I516",-2.1,-1.615,0.955,idrotm[573],"MANY");
  gMC->Gspos("ITS5",1,"I562",0.0,0.0,0.0,0,"ONLY");
  
  
  // --- Place subdetectors' mother volumes into ITS mother volume ITSD
    
  gMC->Gspos("IT12",1,"ITSD",0.0,0.0,0.0,0,"ONLY"); 
  gMC->Gspos("IT34",1,"ITSD",0.0,0.0,0.0,0,"ONLY"); 
  gMC->Gspos("IT56",1,"ITSD",0.0,0.0,0.0,0,"ONLY"); 
  //gMC->Gspos("IS01",1,"ITSD",0.0,0.0,0.0,0,"ONLY"); 
  //gMC->Gspos("IS02",1,"ITSD",0.0,0.0,0.0,0,"ONLY");         
  

  
  
  
  
  
  
  // ********************************************************************





  // SERVICES
  
  
  // --- DEFINE CABLES AT THE END OF THE ITS CONES - COPPER PART
  
  dgh[0] = 45.;
  dgh[1] = 45.+1.0;
  dgh[2] = 9.5;
  
  gMC->Gsvolu("ICCU", "TUBE", idtmed[279], dgh, 3);  
  gMC->Gspos("ICCU", 1, "ITSV", 0., 0., 86.7, 0, "ONLY");
  gMC->Gspos("ICCU", 2, "ITSV", 0., 0., -86.7, idrotm[200], "ONLY");
  
  // --- DEFINE CABLES AT THE END OF THE ITS CONES - CARBON PART
  
  dgh[0] = 45.+1.0;
  dgh[1] = 45.+1.0+1.5;
  dgh[2] = 9.5;
  
  gMC->Gsvolu("ICCC", "TUBE", idtmed[274], dgh, 3);  
  gMC->Gspos("ICCC", 1, "ITSV", 0., 0., 86.7, 0, "ONLY");
  gMC->Gspos("ICCC", 2, "ITSV", 0., 0., -86.7, idrotm[200], "ONLY");  
  
  // --- DEFINE PATCH PANELS AT THE END OF THE ITS CONES
  
  dgh[0] = 45.;
  dgh[1] = 56.;
  dgh[2] = 2.25;
  
  gMC->Gsvolu("IPAN", "TUBE", idtmed[285], dgh, 3);  
  gMC->Gspos("IPAN", 1, "ITSV", 0., 0., 98.45, 0, "ONLY");  
  gMC->Gspos("IPAN", 2, "ITSV", 0., 0., -98.45, idrotm[200], "ONLY"); 
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC - COPPER PART - UPPER PART
 
  dgh[0] = (xltpc-100.7)/2.;
  dgh[1] = 45.2;
  dgh[2] = 45.2+1.0;
  dgh[3] = 61.8;
  dgh[4] = 61.8+1.0;
  dgh[5] = 12.;    
  dgh[6] = 168.;
  gMC->Gsvolu("ICU1", "CONS", idtmed[279], dgh, 7);    
  gMC->Gspos("ICU1", 1, "ITSV", 0., 0., 100.7+dgh[0], 0, "ONLY");  
  gMC->Gspos("ICU1", 2, "ITSV", 0., 0., -(100.7+dgh[0]), idrotm[200], "ONLY");   
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC - COPPER PART - LOWER PART
  
  dgh[0] = (xltpc-100.7)/2.;
  dgh[1] = 45.2;
  dgh[2] = 45.2+1.0;
  dgh[3] = 61.8;
  dgh[4] = 61.8+1.0;
  dgh[5] = 192.;    
  dgh[6] = 348.;
  gMC->Gsvolu("ICU2", "CONS", idtmed[279], dgh, 7);    
  gMC->Gspos("ICU2", 1, "ITSV", 0., 0., 100.7+dgh[0], 0, "ONLY");  
  gMC->Gspos("ICU2", 2, "ITSV", 0., 0., -(100.7+dgh[0]), idrotm[200], "ONLY");     
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC - CARBON PART - UPPER PART
  
  dgh[0] = (xltpc-100.7)/2.;
  dgh[1] = 45.2+1.0;
  dgh[2] = 45.2+1.0+1.5;
  dgh[3] = 61.8+1.0;
  dgh[4] = 61.8+1.0+1.5;
  dgh[5] = 12.;    
  dgh[6] = 168.;  
  gMC->Gsvolu("ICC1", "CONS", idtmed[274], dgh, 7);    
  gMC->Gspos("ICC1", 1, "ITSV", 0., 0., 100.7+dgh[0], 0, "ONLY");  
  gMC->Gspos("ICC1", 2, "ITSV", 0., 0., -(100.7+dgh[0]), idrotm[200], "ONLY");   
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC - CARBON PART - LOWER PART
  
  dgh[0] = (xltpc-100.7)/2.;
  dgh[1] = 45.2+1.0;
  dgh[2] = 45.2+1.0+1.5;
  dgh[3] = 61.8+1.0;
  dgh[4] = 61.8+1.0+1.5;
  dgh[5] = 192.;    
  dgh[6] = 348.;  
  gMC->Gsvolu("ICC2", "CONS", idtmed[274], dgh, 7);    
  gMC->Gspos("ICC2", 1, "ITSV", 0., 0., 100.7+dgh[0], 0, "ONLY");  
  gMC->Gspos("ICC2", 2, "ITSV", 0., 0., -(100.7+dgh[0]), idrotm[200], "ONLY");     
    
  // --- DEFINE CABLES/COOLING BEHIND THE TPC - COPPER PART - UPPER PART
    
  dgh[0] = 62.5;
  dgh[1] = 74.5;
  dgh[2] = 0.5;
  dgh[3] = 12.;
  dgh[4] = 168.;
  gMC->Gsvolu("ICU3", "TUBS", idtmed[279], dgh, 5);    
  gMC->Gspos("ICU3", 1, "ITSV", 0., 0., xltpc+1.5+dgh[2], 0, "ONLY");  
  gMC->Gspos("ICU3", 2, "ITSV", 0., 0., -(xltpc+1.5+dgh[2]), idrotm[200], "ONLY");      
  
  // --- DEFINE CABLES/COOLING BEHIND THE TPC - COPPER PART - LOWER PART
  
  dgh[0] = 62.5;
  dgh[1] = 74.5;
  dgh[2] = 0.5;
  dgh[3] = 192.;
  dgh[4] = 348.;
  gMC->Gsvolu("ICU4", "TUBS", idtmed[279], dgh, 5);    
  gMC->Gspos("ICU4", 1, "ITSV", 0., 0., xltpc+1.5+dgh[2], 0, "ONLY");  
  gMC->Gspos("ICU4", 2, "ITSV", 0., 0., -(xltpc+1.5+dgh[2]), idrotm[200], "ONLY");      
     
  // --- DEFINE CABLES/COOLING BEHIND THE TPC - CARBON PART - UPPER PART

  dgh[0] = 62.5;
  dgh[1] = 74.5;
  dgh[2] = 0.75;
  dgh[3] = 12.;
  dgh[4] = 168.;
  gMC->Gsvolu("ICC3", "TUBS", idtmed[274], dgh, 5);    
  gMC->Gspos("ICC3", 1, "ITSV", 0., 0., xltpc+dgh[2], 0, "ONLY");  
  gMC->Gspos("ICC3", 2, "ITSV", 0., 0., -(xltpc+dgh[2]), idrotm[200], "ONLY"); 
    
  // --- DEFINE CABLES/COOLING BEHIND THE TPC - CARBON PART - LOWER PART

  dgh[0] = 62.5;
  dgh[1] = 74.5;
  dgh[2] = 0.75;
  dgh[3] = 192.;
  dgh[4] = 348.;
  gMC->Gsvolu("ICC4", "TUBS", idtmed[274], dgh, 5);    
  gMC->Gspos("ICC4", 1, "ITSV", 0., 0., xltpc+dgh[2], 0, "ONLY");  
  gMC->Gspos("ICC4", 2, "ITSV", 0., 0., -(xltpc+dgh[2]), idrotm[200], "ONLY"); 

  // --- DEFINE HOOK TO THE TPC ON OTHER SIDE W.R.T. THE ABSORBER - UPPER PART
  
  dgh[0] = 74.5;
  dgh[1] = 79.5;
  dgh[2] = 2.5;
  dgh[3] = 12.;
  dgh[4] = 168.;
  gMC->Gsvolu("IHK1", "TUBS", idtmed[284], dgh, 5);   
  gMC->Gspos("IHK1", 1, "ITSV", 0., 0., -xltpc-dgh[2], 0, "ONLY");      
  
  // --- DEFINE HOOK TO THE TPC ON OTHER SIDE W.R.T. THE ABSORBER - LOWER PART
  
  dgh[0] = 74.5;
  dgh[1] = 79.5;
  dgh[2] = 2.5;
  dgh[3] = 192.;
  dgh[4] = 348.;
  gMC->Gsvolu("IHK2", "TUBS", idtmed[284], dgh, 5);   
  gMC->Gspos("IHK2", 1, "ITSV", 0., 0., -xltpc-dgh[2], 0, "ONLY");        
  
  // --- DEFINE RAILS BETWEEN THE ITS AND THE TPC
  
  //dgh[0] = 0.85;
  //dgh[1] = 10.;
  //dgh[2] = 190.;  
  //gMC->Gsvolu("IRAI", "BOX ", idtmed[285], dgh, 3);   
  //gMC->Gspos("IRAI", 1, "ITSV", 53., 0., -69.5, 0, "ONLY");
  //gMC->Gspos("IRAI", 2, "ITSV", -53., 0., -69.5, 0, "ONLY");        
  
  // --- DEFINE CYLINDERS HOLDING RAILS BETWEEN THE ITS AND THE TPC
  

  dgh[0] = 58.;
  dgh[1] = 59.;
  dgh[2] = 0.6;   
  gMC->Gsvolu("ICYL", "TUBE", idtmed[285], dgh, 3);   
  gMC->Gspos("ICYL", 1, "ALIC", 0., 0., 74., 0, "ONLY");
  gMC->Gspos("ICYL", 2, "ALIC", 0., 0., -74., idrotm[200], "ONLY");     
  
  // --- Outputs the geometry tree in the EUCLID/CAD format 
  
  if (fEuclidOut) {
    gMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
  }
}
//_____________________________________________________________________________
void AliITSvPPRsymm::CreateMaterials(){
////////////////////////////////////////////////////////////////////////
  //
  // Create ITS materials
  //     This function defines the default materials used in the Geant
  // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
  // AliITSvPPRsymm, AliITSvPPRasymm.
  // In general it is automatically replaced by
  // the CreatMaterials routine defined in AliITSv?. Should the function
  // CreateMaterials not exist for the geometry version you are using this
  // one is used. See the definition found in AliITSv5 or the other routine
  // for a complete definition.
  //
  // Water H2O
  Float_t awat[2]  = { 1.00794,15.9994 };
  Float_t zwat[2]  = { 1.,8. };
  Float_t wwat[2]  = { 2.,1. };
  Float_t denswat  = 1.;
  // Freon
  Float_t afre[2]  = { 12.011,18.9984032 };
  Float_t zfre[2]  = { 6.,9. };
  Float_t wfre[2]  = { 5.,12. };
  Float_t densfre  = 1.5;
  // Ceramics
  //     94.4% Al2O3 , 2.8% SiO2 , 2.3% MnO , 0.5% Cr2O3 
  Float_t acer[5]  = { 26.981539,15.9994,28.0855,54.93805,51.9961 };
  Float_t zcer[5]  = { 13.,8.,14.,25.,	    24. };
  Float_t wcer[5]  = { .49976,1.01233,.01307,	    .01782,.00342 };
  Float_t denscer  = 3.6;
  //
  //     60% SiO2 , 40% G10FR4 
  // PC board
  Float_t apcb[3]  = { 28.0855,15.9994,17.749 };
  Float_t zpcb[3]  = { 14.,8.,8.875 };
  Float_t wpcb[3]  = { .28,.32,.4 };
  Float_t denspcb  = 1.8;
  // POLYETHYL
  Float_t apoly[2] = { 12.01,1. };
  Float_t zpoly[2] = { 6.,1. };
  Float_t wpoly[2] = { .33,.67 };
  // old SERVICES
  Float_t zserv[4] = { 1.,6.,26.,29. };
  Float_t aserv[4] = { 1.,12.,55.8,63.5 };
  Float_t wserv[4] = { .014,.086,.42,.48 };
  // Stainless steel
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  
  
  Int_t  isxfld  = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  
  // --- Define the various materials for GEANT --- 
  
  //  200-224 --> Silicon Pixel Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(0, "SPD Si$",      28.0855, 14., 2.33, 9.36, 999);
  AliMaterial(1, "SPD Si chip$", 28.0855, 14., 2.33, 9.36, 999);
  AliMaterial(2, "SPD Si bus$",  28.0855, 14., 2.33, 9.36, 999);
  AliMaterial(3, "SPD C$",       12.011,   6., 2.265,18.8, 999);
  // v. dens 
  AliMaterial(4, "SPD Air$",    14.61, 7.3, .001205, 30423., 999);
  AliMaterial(5, "SPD Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(6, "SPD Al$",     26.981539, 13., 2.6989, 8.9, 999);
  AliMixture( 7, "SPD Water $", awat, zwat, denswat, -2, wwat);
  AliMixture( 8, "SPD Freon$",  afre, zfre, densfre, -2, wfre);
  AliMaterial(9, "SPD End ladder$", 55.845, 26., 7.87/10., 1.76*10., 999); 
  //AliMaterial(9, "SPD End ladder$", 55.845, 26., -7.87/10., -1.76*10., 999); 
  AliMaterial(10, "SPD cone$",28.0855, 14., 2.33, 9.36, 999);       // check !!!!
  // ** 
  AliMedium(0, "SPD Si$",       0, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(1, "SPD Si chip$",  1, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(2, "SPD Si bus$",   2, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(3, "SPD C$",        3, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(4, "SPD Air$",      4, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(5, "SPD Vacuum$",   5, 0,isxfld,sxmgmx, 10.,1.00, .1, .100,10.00);
  AliMedium(6, "SPD Al$",       6, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(7, "SPD Water $",   7, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(8, "SPD Freon$",    8, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(9, "SPD End ladder",9, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(10, "SPD cone$",   10, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);   
  
  //  225-249 --> Silicon Drift Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(25, "SDD Si$",      28.0855, 14., 2.33,  9.36, 999);
  AliMaterial(26, "SDD Si chip$", 28.0855, 14., 2.33,  9.36, 999);
  AliMaterial(27, "SDD Si bus$",  28.0855, 14., 2.33,  9.36, 999);
  AliMaterial(28, "SDD C$",       12.011,   6., 2.265,18.8,  999);
  // v. dens 
  AliMaterial(29, "SDD Air$",     14.61, 7.3, .001205, 30423., 999);
  AliMaterial(30, "SDD Vacuum$",  1e-16, 1e-16, 1e-16, 1e16,  1e16);
  AliMaterial(31, "SDD Al$",      26.981539, 13., 2.6989, 8.9, 999);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(32, "SDD Water $", awat, zwat, denswat, 2, wwat);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture( 33, "SDD Freon$", afre, zfre, densfre, 2, wfre);
  AliMixture( 34, "SDD PCB$",   apcb, zpcb, denspcb, 3, wpcb);
  AliMaterial(35, "SDD Copper$", 63.546, 29., 8.96, 1.43, 999);
  AliMixture( 36, "SDD Ceramics$", acer, zcer, denscer, -5, wcer);
  AliMaterial(37, "SDD Kapton$", 12.011, 6., 1.3, 31.27, 999);
  AliMaterial(38, "SDD End ladder$", 69.9298, 29.8246, 0.3824, 36.5103, 999); 
  AliMaterial(39, "SDD cone$", 63.546, 29., 1.15, 1.265, 999);  
  //AliMaterial(38, "SDD End ladder$", 69.9298, 29.8246, -0.3824, -36.5103, 999); 
  //AliMaterial(39, "SDD cone$", 63.546, 29., -1.15, -1.265, 999);         
  // ** 
  // check A and Z 
  AliMedium(25, "SDD Si$",      25, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(26, "SDD Si chip$", 26, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(27, "SDD Si bus$",  27, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(28, "SDD C$",       28, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(29, "SDD Air$",     29, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(30, "SDD Vacuum$",  30, 0,isxfld,sxmgmx, 10.,1.00, .1, .100,10.00);
  AliMedium(31, "SDD Al$",      31, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(32, "SDD Water $",  32, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(33, "SDD Freon$",   33, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(34, "SDD PCB$",     34, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(35, "SDD Copper$",  35, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(36, "SDD Ceramics$",36, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(37, "SDD Kapton$",  37, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(38, "SDD End ladder",38, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(39, "SDD cone$",   39, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  
  //  250-274 --> Silicon Strip Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(50, "SSD Si$",      28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(51, "SSD Si chip$", 28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(52, "SSD Si bus$",  28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(53, "SSD C$",       12.011,   6., 2.265,18.8, 999.);
  // v. dens 
  AliMaterial(54, "SSD Air$",     14.61, 7.3, .001205, 30423., 999);
  AliMaterial(55, "SSD Vacuum$",  1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(56, "SSD Al$",      26.981539, 13., 2.6989, 8.9, 999);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(57, "SSD Water $", awat, zwat, denswat, 2, wwat);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(58, "SSD Freon$", afre, zfre, densfre, 2, wfre);
  AliMixture(59, "SSD PCB$",   apcb, zpcb, denspcb, 3, wpcb);
  AliMaterial(60, "SSD Copper$", 63.546, 29., 8.96, 1.43, 999.);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture( 61, "SSD Ceramics$", acer, zcer, denscer, 5, wcer);
  AliMaterial(62, "SSD Kapton$", 12.011, 6., 1.3, 31.27, 999.);
  // check A and Z 
  AliMaterial(63, "SSD G10FR4$", 17.749, 8.875, 1.8, 21.822, 999.);
  AliMaterial(64, "SSD End ladder$", 32.0988, 15.4021, 0.68, 35.3238, 999); 
  AliMaterial(65, "SSD cone$",63.546, 29., 1.15, 1.265, 999);   
  //AliMaterial(64, "SSD End ladder$", 32.0988, 15.4021, -0.68, -35.3238, 999); 
  //AliMaterial(65, "SSD cone$",63.546, 29., -1.15, -1.265, 999);     
  // ** 
  AliMedium(50, "SSD Si$",      50, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(51, "SSD Si chip$", 51, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(52, "SSD Si bus$",  52, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(53, "SSD C$",       53, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(54, "SSD Air$",     54, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(55, "SSD Vacuum$",  55, 0,isxfld,sxmgmx, 10.,1.00, .1, .100,10.00);
  AliMedium(56, "SSD Al$",      56, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(57, "SSD Water $",  57, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(58, "SSD Freon$",   58, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(59, "SSD PCB$",     59, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(60, "SSD Copper$",  60, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(61, "SSD Ceramics$",61, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(62, "SSD Kapton$",  62, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(63, "SSD G10FR4$",  63, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(64, "SPD End ladder",64, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(65, "SPD cone$",   65, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  
  //     275-299 --> General (end-caps, frames, cooling, cables, etc.) 
  
  AliMaterial(75, "GEN C$", 12.011, 6., 2.265, 18.8, 999.);
  // verify density 
  AliMaterial(76, "GEN Air$", 14.61, 7.3, .001205, 30423., 999);
  AliMaterial(77, "GEN Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMixture( 78, "GEN POLYETHYL$", apoly, zpoly, .95, -2, wpoly);
  AliMixture( 79, "GEN SERVICES$",  aserv, zserv, 4.68, 4, wserv);
  AliMaterial(80, "GEN Copper$", 63.546, 29., 8.96, 1.43, 999.);  
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(81, "GEN Water $", awat, zwat, denswat, 2, wwat);
  //AliMaterial(82, "GEN Cables$", 12.011, 6., 2.265, 18.8, 999.);  // check !!!
  //AliMaterial(83, "GEN patch pan$", 12.011, 6., 2.265, 18.8, 999.);  // check !!!  
  //AliMaterial(84, "GEN serv$", 12.011, 6., 2.265, 18.8, 999.);  // check !!!  
  AliMixture(85, "GEN Inox$", asteel, zsteel, 7.88, 4, wsteel);
  AliMaterial(86, "GEN Al$",      26.981539, 13., 2.6989, 8.9, 999);
  // ** 
  AliMedium(75,"GEN C$",        75, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(76,"GEN Air$",      76, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(77,"GEN Vacuum$",   77, 0,isxfld,sxmgmx, 10., .10, .1, .100,10.00);
  AliMedium(78,"GEN POLYETHYL$",78, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(79,"GEN SERVICES$", 79, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(80,"GEN Copper$",   80, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(81,"GEN Water $",   81, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  //AliMedium(82,"GEN Cables$",   82, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  //AliMedium(83,"GEN patch pan$",83, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);  
  //AliMedium(84,"GEN serv$",     84, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(85,"GEN Inox$",     85, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(86, "GEN Al$",      86, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);

}
//_____________________________________________________________________________
void AliITSvPPRsymm::Init(){
////////////////////////////////////////////////////////////////////////
//     Initialise the ITS after it has been created.
////////////////////////////////////////////////////////////////////////

    //
    AliITS::Init();
    fMajorVersion = 1;
    fMinorVersion = 0;
}  
 
//_____________________________________________________________________________
void AliITSvPPRsymm::DrawModule(){
////////////////////////////////////////////////////////////////////////
//     Draw a shaded view of the FMD version 9.
////////////////////////////////////////////////////////////////////////
  
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother visible
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("ITSD","SEEN",0);
  gMC->Gsatt("ITS1","SEEN",1);
  gMC->Gsatt("ITS2","SEEN",1);
  gMC->Gsatt("ITS3","SEEN",1);
  gMC->Gsatt("ITS4","SEEN",1);
  gMC->Gsatt("ITS5","SEEN",1);
  gMC->Gsatt("ITS6","SEEN",1);

  gMC->Gsatt("IPCB","SEEN",1);
  gMC->Gsatt("ICO2","SEEN",1);
  gMC->Gsatt("ICER","SEEN",0);
  gMC->Gsatt("ISI2","SEEN",0);
  gMC->Gsatt("IPLA","SEEN",0);
  gMC->Gsatt("ICO3","SEEN",0);
  gMC->Gsatt("IEPX","SEEN",0);
  gMC->Gsatt("ISI3","SEEN",1);
  gMC->Gsatt("ISUP","SEEN",0);
  gMC->Gsatt("ICHO","SEEN",0);
  gMC->Gsatt("ICMO","SEEN",0);
  gMC->Gsatt("ICMD","SEEN",0);
  gMC->Gsatt("ICCO","SEEN",1);
  gMC->Gsatt("ICCM","SEEN",0);
  gMC->Gsatt("ITMD","SEEN",0);
  gMC->Gsatt("ITTT","SEEN",1);

  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 300, -300, 300, -300, 300);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 11, 10, .07, .07);
  gMC->Gdhead(1111, "Inner Tracking System Version 1");
  gMC->Gdman(17, 6, "MAN");
}
//_____________________________________________________________________________
void AliITSvPPRsymm::StepManager(){
////////////////////////////////////////////////////////////////////////
//    Called for every step in the ITS, then calls the AliITShit class
// creator with the information to be recoreded about that hit.
////////////////////////////////////////////////////////////////////////

/*
  Int_t         copy, id;
  Float_t       hits[8];
  Int_t         vol[4];
  TLorentzVector position, momentum;
//  TClonesArray &lhits = *fHits;
//
// no hits for this  symmetric version.
//

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
  if( !(gMC->TrackCharge()) ) return;
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
//    new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->CurrentTrack(),vol,hits);
//
// no hits for this  symmetric version.
//
*/
}
/*
//____________________________________________________________________________
void AliITSvPPRsymm::Streamer(TBuffer &R__b){
////////////////////////////////////////////////////////////////////////
//    A dummy Streamer function for this class AliITSvPPRsymm. By default it
// only streams the AliITS class as it is required. Since this class
// dosen't contain any "real" data to be saved, it doesn't.
////////////////////////////////////////////////////////////////////////

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITS::Streamer(R__b);
   } else {
      R__b.WriteVersion(AliITSvPPRsymm::IsA());
      AliITS::Streamer(R__b);
   } // end if R__b.IsReading()
}
*/
