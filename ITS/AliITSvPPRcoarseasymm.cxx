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
//  Inner Traking System version PPR coarse asymmetric                                         //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
// Authors: R. Barbera
// version 6.
// Created  2000.
//
//  NOTE: THIS IS THE COARSE ASYMMETRIC PPR geometry of the ITS. 
// THIS WILL NOT WORK
// with the geometry or module classes or any analysis classes. You are 
// strongly encouraged to uses AliITSv5.
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <Riostream.h> 
#include <TMath.h>
#include <TNode.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TTUBE.h>
#include <TVector.h>
#include <TVirtualMC.h>
#include <TGeometry.h>


#include "AliMagF.h"
#include "AliConst.h"
#include "AliITShit.h"
#include "AliITSvPPRcoarseasymm.h"
#include "AliMagF.h"
#include "AliRun.h"

ClassImp(AliITSvPPRcoarseasymm)
 
//_____________________________________________________________________________
AliITSvPPRcoarseasymm::AliITSvPPRcoarseasymm():
fMajorVersion(9),
fMinorVersion(0),
fRails(0),
fSuppMat(0) {
////////////////////////////////////////////////////////////////////////
//    Standard default constructor for the ITS version 6.
////////////////////////////////////////////////////////////////////////

    fIdN = 0;
    fIdName = 0;
    fIdSens    = 0;
}
//_____________________________________________________________________________
AliITSvPPRcoarseasymm::AliITSvPPRcoarseasymm(const char *name, const char *title) : AliITS(name, title),
fMajorVersion(9),
fMinorVersion(0),
fRails(0),
fSuppMat(0) {
////////////////////////////////////////////////////////////////////////
//    Standard constructor for the ITS version 6.
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
    for (Int_t i=0;i<fIdN;i++) fIdSens[i]=0;
}
//____________________________________________________________________________
AliITSvPPRcoarseasymm::AliITSvPPRcoarseasymm(const AliITSvPPRcoarseasymm &s) :
 AliITS(s),
fMajorVersion(s.fMajorVersion),
fMinorVersion(s.fMinorVersion),
fRails(s.fRails),
fSuppMat(s.fSuppMat){
////////////////////////////////////////////////////////////////////////
//     Copy Constructor for ITS version 6.
////////////////////////////////////////////////////////////////////////
}
//_____________________________________________________________________________
AliITSvPPRcoarseasymm& AliITSvPPRcoarseasymm::operator=(const AliITSvPPRcoarseasymm &source){
////////////////////////////////////////////////////////////////////////
//    Assignment operator for the ITS version 6.
////////////////////////////////////////////////////////////////////////
  if(&source == this) return *this;
    Warning("= operator","Not allowed to copy AliITSvPPRcoarseasymm");
  return *this;
}
//_____________________________________________________________________________
AliITSvPPRcoarseasymm::~AliITSvPPRcoarseasymm() {
////////////////////////////////////////////////////////////////////////
//    Standard destructor for the ITS version 6.
////////////////////////////////////////////////////////////////////////
}

//__________________________________________________________________________
void AliITSvPPRcoarseasymm::BuildGeometry(){
////////////////////////////////////////////////////////////////////////
//    Geometry builder for the ITS version 6.
////////////////////////////////////////////////////////////////////////
    TNode *node, *top;
    const Int_t kColorITS=kYellow;
    //
    top = gAlice->GetGeometry()->GetNode("alice");

    new TTUBE("S_layer1","Layer1 of ITS","void",3.8095,3.8095+1.03*9.36/100.,14.35);
    top->cd();
    node = new TNode("Layer1","Layer1","S_layer1",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer2","Layer2 of ITS","void",7.,7.+1.03*9.36/100.,14.35);
    top->cd();
    node = new TNode("Layer2","Layer2","S_layer2",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer3","Layer3 of ITS","void",15.,15.+0.94*9.36/100.,25.1);
    top->cd();
    node = new TNode("Layer3","Layer3","S_layer3",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer4","Layer4 of ITS","void",24.1,24.1+0.95*9.36/100.,32.1);
    top->cd();
    node = new TNode("Layer4","Layer4","S_layer4",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer5","Layer5 of ITS","void",38.5,38.5+0.91*9.36/100.,49.405);
    top->cd();
    node = new TNode("Layer5","Layer5","S_layer5",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);

    new TTUBE("S_layer6","Layer6 of ITS","void",43.5765,43.5765+0.87*9.36/100.,55.27);
    top->cd();
    node = new TNode("Layer6","Layer6","S_layer6",0,0,0,"");
    node->SetLineColor(kColorITS);
    fNodes->Add(node);
}
//_____________________________________________________________________________
void AliITSvPPRcoarseasymm::CreateGeometry(){
////////////////////////////////////////////////////////////////////////
//    This routine defines and Creates the geometry for version 6 of the ITS.
////////////////////////////////////////////////////////////////////////
  
  //INNER RADII OF THE SILICON LAYERS 
  Float_t rl[6]    = { 3.8095,7.,15.,24.,38.5,43.5765 };   
  //THICKNESSES OF LAYERS (in % radiation length)
  Float_t drl[6]   = { 1.03,1.03,0.94,0.95,0.91,0.87 };   
  //HALF LENGTHS OF LAYERS  
  Float_t dzl[6]   = { 14.35,14.35,25.1,32.1,49.405,55.27 };
  //LENGTHS OF END-LADDER BOXES (ALL INCLUDED)
  Float_t dzb[6]   = { 12.4,12.4,13.5,15.,7.5,7.5 };   
  //THICKNESSES OF END-LADDER BOXES (ALL INCLUDED)
  Float_t drb[6]   = { rl[1]-rl[0],0.2,5.,5.,4.,4. };        

 
  Float_t dits[3], rlim, zmax;
  Float_t zpos;
  Float_t pcits[100], ztpc;
  Int_t idrotm[1999], i;
  Float_t dgh[100];
  
  Int_t rails = 1;       // flag for rails (1 --> rails in; 0 --> rails out)
  Int_t suppmat = 0;     // flag to change the material of the services
                         // supports (=0 copper, =1 aluminum, =2 carbon)  
  rails = GetRails();
  
  if(rails != 0 && rails != 1) {
     cout << "ITS - WARNING: the switch for rails is not set neither to 0 (rails out) nor to 1 (rails in)." 
     " The default value of 1 (rails in) will be used." << endl;
  }  
  
  if (rails == 0 ) {
     cout << "ITS: Rails are out." << endl; 
  } else {
     cout << "ITS: Rails are in." << endl;
  }      
  
  suppmat = GetSupportMaterial();   
  
  if (suppmat != 0 && suppmat != 1 && suppmat != 2) {
     cout << "ITS - WARNING: the flag for the material of services supports is not set neither to 0 (copper) nor to 1 (aluminum) nor to 2 (carbon)." 
     " The default value of 0 (copper) will be used." << endl;
  }  
  
  if (suppmat == 0) {
     cout << "ITS: The material of the services supports is copper." << endl; 
  } else if (suppmat == 1){
     cout << "ITS: The material of the services supports is aluminum." << endl;
  } else {
     cout << "ITS: The material of the services supports is carbon." << endl;
  }      
      

  Int_t *idtmed = fIdtmed->GetArray()-199;
  
  //     CONVERT INTO CM (RL(SI)=9.36 CM) 
  for (i = 0; i < 6; ++i) {
    drl[i] = drl[i] / 100. * 9.36;
  }
    
  //     FIELD CAGE HALF LENGTH 
  
  rlim  = 50.;
  zmax  = 74.;
  ztpc = 284.;
  
  // --- Define ghost volume containing the whole ITS (including services) 
  //     and fill it with air 
  
  dgh[0] = 0.;
  dgh[1] = 360.;
  dgh[2] = 16.;
  dgh[3] = -ztpc-5.-0.1;
  dgh[4] = 46;   
  dgh[5] = 85.;
  dgh[6] = -ztpc;
  dgh[7] = 46;   
  dgh[8] = 85.;
  dgh[9] = -ztpc;
  dgh[10] = 46;  
  dgh[11] = rlim+6;
  dgh[12] = -97.5;
  dgh[13] = 46;  
  dgh[14] = rlim+6;
  dgh[15] = -zmax;
  dgh[16] = 46;  
  dgh[17] = rlim+6;
  dgh[18] = -48;   
  dgh[19] = 6;
  dgh[20] = rlim+6;
  dgh[21] = -28.6;   
  dgh[22] = 6;
  dgh[23] = rlim+6;    
  dgh[24] = -27.6;  
  dgh[25] = 3.295;
  dgh[26] = rlim+6; 
  dgh[27] = 27.6;   
  dgh[28] = 3.295;
  dgh[29] = rlim+6;
  dgh[30] = 28.6;   
  dgh[31] = 6;
  dgh[32] = rlim+6;
  dgh[33] = 48;   
  dgh[34] = 6;
  dgh[35] = rlim+6;  
  dgh[36] = zmax;
  dgh[37] = 46;
  dgh[38] = rlim+6;
  dgh[39] = 97.5;
  dgh[40] = 46;  
  dgh[41] = rlim+6;
  dgh[42] = ztpc;
  dgh[43] = 62;     
  dgh[44] = 62+4.;  
  dgh[45] = ztpc;
  dgh[46] = 62;     
  dgh[47] = 85.;
  dgh[48] = ztpc+4.+0.1;
  dgh[49] = 62.4;
  dgh[50] = 85.;
  gMC->Gsvolu("ITSV", "PCON", idtmed[275], dgh, 51);
  
  // --- Place the ghost volume in its mother volume (ALIC) and make it 
  //     invisible 
  
  gMC->Gspos("ITSV", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  //gMC->Gsatt("ITSV", "SEEN", 0); 


  // --- Define ghost volume containing the six layers and fill it with air 
  
  dgh[0] = 0.;
  dgh[1] = 360.;
  dgh[2] = 8.;
  dgh[3] = -zmax;  
  dgh[4] = 46.;
  dgh[5] = rlim;
  dgh[6] = -47.5;    
  dgh[7] = 6.005;
  dgh[8] = rlim;
  dgh[9] = -28.5;    
  dgh[10] = 6.005;
  dgh[11] = rlim;  
  dgh[12] = -27.5;   
  dgh[13] = 3.3;
  dgh[14] = rlim;
  dgh[15] = 27.5;    
  dgh[16] = 3.3;
  dgh[17] = rlim;
  dgh[18] = 28.5;    
  dgh[19] = 6.005;
  dgh[20] = rlim;
  dgh[21] = 47.5;    
  dgh[22] = 6.005;
  dgh[23] = rlim;
  dgh[24] = zmax;    
  dgh[25] = 46.;
  dgh[26] = rlim;
  gMC->Gsvolu("ITSD", "PCON", idtmed[275], dgh, 27);
  
  // --- Place the ghost volume in its mother volume (ITSV) and make it 
  //     invisible 
  
  gMC->Gspos("ITSD", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  //gMC->Gsatt("ITSD", "SEEN", 0);
  

  //     ITS LAYERS (SILICON) 
  
  dits[0] = rl[0];
  dits[1] = rl[0] + drl[0];
  dits[2] = dzl[0];
  gMC->Gsvolu("ITS1", "TUBE", idtmed[199], dits, 3);
  gMC->Gspos("ITS1", 1, "ITSD", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[1];
  dits[1] = rl[1] + drl[1];
  dits[2] = dzl[1];
  gMC->Gsvolu("ITS2", "TUBE", idtmed[199], dits, 3);
  gMC->Gspos("ITS2", 1, "ITSD", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[2];
  dits[1] = rl[2] + drl[2];
  dits[2] = dzl[2];
  gMC->Gsvolu("ITS3", "TUBE", idtmed[224], dits, 3);
  gMC->Gspos("ITS3", 1, "ITSD", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[3];
  dits[1] = rl[3] + drl[3];
  dits[2] = dzl[3];
  gMC->Gsvolu("ITS4", "TUBE", idtmed[224], dits, 3);
  gMC->Gspos("ITS4", 1, "ITSD", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[4];
  dits[1] = rl[4] + drl[4];
  dits[2] = dzl[4];
  gMC->Gsvolu("ITS5", "TUBE", idtmed[249], dits, 3);
  gMC->Gspos("ITS5", 1, "ITSD", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[5];
  dits[1] = rl[5] + drl[5];
  dits[2] = dzl[5];
  gMC->Gsvolu("ITS6", "TUBE", idtmed[249], dits, 3);
  gMC->Gspos("ITS6", 1, "ITSD", 0., 0., 0., 0, "ONLY");
  
  // END-LADDER ELECTRONICS BOXES AND CABLES FOR SPD
  
  gMC->Gsvolu("IEL1", "TUBE", idtmed[208], dits, 0); 
  for (i = 0; i < 2; i++) {
    dits[0] = rl[i];
    dits[1] = dits[0] + drb[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    gMC->Gsposp("IEL1", i+1, "ITSD", 0., 0., zpos, 0, "ONLY", dits, 3);
    gMC->Gsposp("IEL1", i+1+6, "ITSD", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  // END-LADDER ELECTRONICS BOXES AND CABLES FOR SDD
  
  gMC->Gsvolu("IEL2", "TUBE", idtmed[237], dits, 0); 
  for (i = 2; i < 3; i++) {
    dits[0] = rl[i]-2.5;
    dits[1] = dits[0] + drb[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    gMC->Gsposp("IEL2", i+1, "ITSD", 0., 0., zpos, 0, "ONLY", dits, 3);
    gMC->Gsposp("IEL2", i+1+6, "ITSD", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  for (i = 3; i < 4; i++) {
    dits[0] = rl[i]-1.4;
    dits[1] = dits[0] + drb[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    gMC->Gsposp("IEL2", i+1, "ITSD", 0., 0., zpos, 0, "ONLY", dits, 3);
    gMC->Gsposp("IEL2", i+1+6, "ITSD", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  // END-LADDER ELECTRONICS BOXES AND CABLES FOR SSD
  
  gMC->Gsvolu("IEL3", "TUBE", idtmed[263], dits, 0); 
  for (i = 4; i < 5; i++) {
    dits[0] = rl[i]+1.4;
    dits[1] = dits[0] + drb[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    gMC->Gsposp("IEL3", i+1, "ITSD", 0., 0., zpos, 0, "ONLY", dits, 3);
    gMC->Gsposp("IEL3", i+1+6, "ITSD", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }    
  for (i = 5; i < 6; i++) {
    dits[0] = rl[i]+0.4235;
    dits[1] = dits[0] + drb[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    gMC->Gsposp("IEL3", i+1, "ITSD", 0., 0., zpos, 0, "ONLY", dits, 3);
    gMC->Gsposp("IEL3", i+1+6, "ITSD", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }        

  //    DEFINE THERMAL SCREEN FOR SPD
  
  pcits[0] = 8.3;
  pcits[1] = 8.5;
  pcits[2] = 42.5;
  gMC->Gsvolu("ICY1", "TUBE", idtmed[274], pcits, 3);   
  gMC->Gspos("ICY1", 1, "ITSD", 0., 0., 0., 0, "ONLY");

  //    DEFINE END CONES FOR SDD
  
  pcits[0] = (59.-42.5)/2.;
  pcits[1] = 8.5;
  pcits[2] = 8.5+0.1;
  pcits[3] = 28.;
  pcits[4] = 28.+0.1;  
  gMC->Gsvolu("ICO1", "CONE", idtmed[238], pcits, 5);    
  AliMatrix(idrotm[200], 90., 0., 90., 90., 180., 0.);
  gMC->Gspos("ICO1", 1, "ITSD", 0., 0., 42.5+pcits[0], 0, "ONLY");
  gMC->Gspos("ICO1", 2, "ITSD", 0., 0., -(42.5+pcits[0]), idrotm[200], "ONLY");

  //    DEFINE CYLINDER BETWEEN SDD AND SSD
  
  pcits[0] = (59.5-0.13/2.)/2.;
  pcits[1] = (59.5+0.13/2.)/2.;
  pcits[2] = 57.;
  gMC->Gsvolu("ICY2", "TUBE", idtmed[274], pcits, 3);   
  gMC->Gspos("ICY2", 1, "ITSD", 0., 0., 0., 0, "ONLY"); 

  //    DEFINE END CONES FOR SSD
  
  pcits[0] = (74.-59.)/2.;
  pcits[1] = 28.;
  pcits[2] = 28.+0.1;
  pcits[3] = 47.;
  pcits[4] = 47.+0.1;
  gMC->Gsvolu("ICO2", "CONE", idtmed[264], pcits, 5);    
  gMC->Gspos("ICO2", 1, "ITSD", 0., 0., 59.+pcits[0], 0, "ONLY");
  gMC->Gspos("ICO2", 2, "ITSD", 0., 0., -(59.+pcits[0]), idrotm[200], "ONLY");
  
    
  // ****************************  SERVICES  *********************************
  

  
  // --- DEFINE CABLES AT THE END OF THE ITS CONES - COPPER PART
  //     UPPER PART
  
  dgh[0] = 46.;    
  dgh[1] = 46.+1.0;  
  dgh[2] = 9.5;
  dgh[3] = 12.;
  dgh[4] = 168.;
  
  if (suppmat == 0) {
     gMC->Gsvolu("I1CU", "TUBS", idtmed[279], dgh, 5);    // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("I1CU", "TUBS", idtmed[285], dgh, 5);    // aluminum
  } else {
     gMC->Gsvolu("I1CU", "TUBS", idtmed[274], dgh, 5);    // carbon
  }     
  gMC->Gspos("I1CU", 1, "ITSV", 0., 0., 83.5, 0, "ONLY");
  gMC->Gspos("I1CU", 2, "ITSV", 0., 0., -83.5, idrotm[200], "ONLY");
  
  // --- DEFINE CABLES AT THE END OF THE ITS CONES - COPPER PART
  //     LOWER PART
  
  dgh[0] = 46.;    
  dgh[1] = 46.+1.0;  
  dgh[2] = 9.5;
  dgh[3] = 192.;
  dgh[4] = 348.;
  
  if (suppmat == 0) {
     gMC->Gsvolu("I2CU", "TUBS", idtmed[279], dgh, 5);    // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("I2CU", "TUBS", idtmed[285], dgh, 5);    // aluminum
  } else {
     gMC->Gsvolu("I2CU", "TUBS", idtmed[274], dgh, 5);    // carbon
  }     
  gMC->Gspos("I2CU", 1, "ITSV", 0., 0., 83.5, 0, "ONLY");
  gMC->Gspos("I2CU", 2, "ITSV", 0., 0., -83.5, idrotm[200], "ONLY");


  // --- DEFINE CABLES AT THE END OF THE ITS CONES - CARBON PART
  //     UPPER PART
  
  dgh[0] = 46.+1.0;  
  dgh[1] = 46.+1.0+1.5;   
  dgh[2] = 9.5;
  dgh[3] = 12.;
  dgh[4] = 168.;
  
  gMC->Gsvolu("I1CC", "TUBS", idtmed[274], dgh, 5);  
  gMC->Gspos("I1CC", 1, "ITSV", 0., 0., 83.5, 0, "ONLY");
  gMC->Gspos("I1CC", 2, "ITSV", 0., 0., -83.5, idrotm[200], "ONLY");  
  
  // --- DEFINE CABLES AT THE END OF THE ITS CONES - CARBON PART
  //     LOWER PART
  
  dgh[0] = 46.+1.0;  
  dgh[1] = 46.+1.0+1.5;   
  dgh[2] = 9.5;
  dgh[3] = 192.;
  dgh[4] = 348.;
  
  gMC->Gsvolu("I2CC", "TUBS", idtmed[274], dgh, 5);  
  gMC->Gspos("I2CC", 1, "ITSV", 0., 0., 83.5, 0, "ONLY");
  gMC->Gspos("I2CC", 2, "ITSV", 0., 0., -83.5, idrotm[200], "ONLY");  

  // --- DEFINE PATCH PANELS AT THE END OF THE ITS CONES
  //    UPPER PART
    
  dgh[0] = 46.;  
  dgh[1] = 56.;
  dgh[2] = 2.25;
  dgh[3] = 12.;
  dgh[4] = 168.;
  
  gMC->Gsvolu("IPA1", "TUBS", idtmed[285], dgh, 5);  
  gMC->Gspos("IPA1", 1, "ITSV", 0., 0., 95.25, 0, "ONLY");  
  gMC->Gspos("IPA1", 2, "ITSV", 0., 0., -95.25, idrotm[200], "ONLY"); 


  // --- DEFINE PATCH PANELS AT THE END OF THE ITS CONES
  //     LOWER PART
    
  dgh[0] = 46.;  
  dgh[1] = 56.;
  dgh[2] = 2.25;
  dgh[3] = 192.;
  dgh[4] = 348.;
  
  gMC->Gsvolu("IPA2", "TUBS", idtmed[285], dgh, 5);  
  gMC->Gspos("IPA2", 1, "ITSV", 0., 0., 95.25, 0, "ONLY");  
  gMC->Gspos("IPA2", 2, "ITSV", 0., 0., -95.25, idrotm[200], "ONLY"); 


  // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE ABSORBER SIDE - COPPER PART
  //     UPPER PART
  
  dgh[0] = (ztpc-97.5)/2.;
  dgh[1] = 46.2;     
  dgh[2] = 46.2+1.0;  
  dgh[3] = 62.3;     
  dgh[4] = 62.3+1.0;   
  dgh[5] = 12.;    
  dgh[6] = 168.;
  if (suppmat == 0) {
     gMC->Gsvolu("ICU1", "CONS", idtmed[279], dgh, 7);    // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU1", "CONS", idtmed[285], dgh, 7);    // aluminum  
  } else {
     gMC->Gsvolu("ICU1", "CONS", idtmed[274], dgh, 7);    // carbon
  }
  gMC->Gspos("ICU1", 1, "ITSV", 0., 0., 97.5+dgh[0], 0, "ONLY");   
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE ABSORBER SIDE - COPPER PART
  //     LOWER PART
  
  dgh[0] = (ztpc-97.5)/2.;
  dgh[1] = 46.2;      
  dgh[2] = 46.2+1.0;  
  dgh[3] = 62.3;      
  dgh[4] = 62.3+1.0;  
  dgh[5] = 192.;    
  dgh[6] = 348.;
  if (suppmat == 0) {
     gMC->Gsvolu("ICU2", "CONS", idtmed[279], dgh, 7);    // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU2", "CONS", idtmed[285], dgh, 7);    // aluminum  
  } else {
     gMC->Gsvolu("ICU2", "CONS", idtmed[274], dgh, 7);    // carbon
  }
  gMC->Gspos("ICU2", 1, "ITSV", 0., 0., 97.5+dgh[0], 0, "ONLY");  


   // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE ABSORBER SIDE - CARBON PART
   //     UPPER PART
  
  dgh[0] = (ztpc-97.5)/2.;
  dgh[1] = 46.2+1.0;      
  dgh[2] = 46.2+1.0+1.5;  
  dgh[3] = 62.3+1.0;      
  dgh[4] = 62.3+1.0+1.5;  
  dgh[5] = 12.;    
  dgh[6] = 168.;  
  gMC->Gsvolu("ICC1", "CONS", idtmed[274], dgh, 7);    
  gMC->Gspos("ICC1", 1, "ITSV", 0., 0., 97.5+dgh[0], 0, "ONLY");   
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE ABSORBER SIDE - CARBON PART
  //     LOWER PART
  
  dgh[0] = (ztpc-97.5)/2.;
  dgh[1] = 46.2+1.0;    
  dgh[2] = 46.2+1.0+1.5;
  dgh[3] = 62.3+1.0;    
  dgh[4] = 62.3+1.0+1.5;
  dgh[5] = 192.;    
  dgh[6] = 348.;  
  gMC->Gsvolu("ICC2", "CONS", idtmed[274], dgh, 7);    
  gMC->Gspos("ICC2", 1, "ITSV", 0., 0., 97.5+dgh[0], 0, "ONLY");  
   
  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON THE ABSORBER SIDE - COPPER PART
  //     UPPER PART
    
  dgh[0] = 62.1; 
  dgh[1] = 74.5;
  dgh[2] = 0.5;
  dgh[3] = 12.;
  dgh[4] = 168.;
  if (suppmat == 0) {
     gMC->Gsvolu("ICU3", "TUBS", idtmed[279], dgh, 5);    // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU3", "TUBS", idtmed[285], dgh, 5);    // aluminum
  } else {
     gMC->Gsvolu("ICU3", "TUBS", idtmed[274], dgh, 5);    // carbon
    }
  gMC->Gspos("ICU3", 1, "ITSV", 0., 0., ztpc+1.5+dgh[2], 0, "ONLY");  

  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON THE ABSORBER SIDE - COPPER PART
  //     LOWER PART
  
  dgh[0] = 62.1;  
  dgh[1] = 74.5;
  dgh[2] = 0.5;
  dgh[3] = 192.;
  dgh[4] = 348.;
  if (suppmat == 0) {
     gMC->Gsvolu("ICU4", "TUBS", idtmed[279], dgh, 5);    // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU4", "TUBS", idtmed[285], dgh, 5);    // aluminum
  } else {
     gMC->Gsvolu("ICU4", "TUBS", idtmed[274], dgh, 5);    // carbon
  }
  gMC->Gspos("ICU4", 1, "ITSV", 0., 0., ztpc+1.5+dgh[2], 0, "ONLY");     
     
  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON THE ABSORBER SIDE - CARBON PART
  //     UPPER PART

  dgh[0] = 62.1;  
  dgh[1] = 74.5;
  dgh[2] = 0.75;
  dgh[3] = 12.;
  dgh[4] = 168.;
  gMC->Gsvolu("ICC3", "TUBS", idtmed[274], dgh, 5);    
  gMC->Gspos("ICC3", 1, "ITSV", 0., 0., ztpc+dgh[2], 0, "ONLY");   
    
  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON THE ABSORBER SIDE - CARBON PART
  //     LOWER PART

  dgh[0] = 62.1;  
  dgh[1] = 74.5;
  dgh[2] = 0.75;
  dgh[3] = 192.;
  dgh[4] = 348.;
  gMC->Gsvolu("ICC4", "TUBS", idtmed[274], dgh, 5);    
  gMC->Gspos("ICC4", 1, "ITSV", 0., 0., ztpc+dgh[2], 0, "ONLY");  
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE OTHER SIDE W.R.T.
  //     THE ABSORBER - COPPER PART - UPPER PART
  
  dgh[0] = 46.;      
  dgh[1] = 46.+1.0;  
  dgh[2] = (ztpc-97.5+1.5)/2.;
  dgh[3] = 12.;
  dgh[4] = 168.;
  if (suppmat == 0) {
     gMC->Gsvolu("ICU5", "TUBS", idtmed[279], dgh, 5);   // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU5", "TUBS", idtmed[285], dgh, 5);   // aluminum
  } else {
     gMC->Gsvolu("ICU5", "TUBS", idtmed[274], dgh, 5);   // carbon
  }
  gMC->Gspos("ICU5", 1, "ITSV", 0., 0., -97.5-dgh[2], 0, "ONLY");  
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE OTHER SIDE W.R.T.
  //     THE ABSORBER - COPPER PART - LOWER PART
  
  dgh[0] = 46.;  
  dgh[1] = 46.+1.0;  
  dgh[2] = (ztpc-97.5+1.5)/2.;
  dgh[3] = 192.;
  dgh[4] = 348.;  
  if (suppmat == 0) {
     gMC->Gsvolu("ICU6", "TUBS", idtmed[279], dgh, 5);   // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU6", "TUBS", idtmed[285], dgh, 5);   // aluminum
  } else {
     gMC->Gsvolu("ICU6", "TUBS", idtmed[274], dgh, 5);   // carbon
  }
  gMC->Gspos("ICU6", 1, "ITSV", 0., 0., -97.5-dgh[2], 0, "ONLY");    
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE OTHER SIDE W.R.T.
  //     THE ABSORBER - CARBON PART - UPPER PART
  
  dgh[0] = 46.+1.0;  
  dgh[1] = 46.+1.0+1.5; 
  dgh[2] = (ztpc-97.5)/2.;
  dgh[3] = 12.;
  dgh[4] = 168.;  
  gMC->Gsvolu("ICC5", "TUBS", idtmed[274], dgh, 5);   
  gMC->Gspos("ICC5", 1, "ITSV", 0., 0., -97.5-dgh[2], 0, "ONLY");   
  
  // --- DEFINE CABLES/COOLING BELOW THE TPC ON THE OTHER SIDE W.R.T.
  //     THE ABSORBER - CARBON PART - LOWER PART
  
  dgh[0] = 46.+1.0;   
  dgh[1] = 46.+1.0+1.5;  
  dgh[2] = (ztpc-97.5)/2.;
  dgh[3] = 192.;
  dgh[4] = 348.;  
  gMC->Gsvolu("ICC6", "TUBS", idtmed[274], dgh, 5);   
  gMC->Gspos("ICC6", 1, "ITSV", 0., 0., -97.5-dgh[2], 0, "ONLY");      

  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON OTHER SIDE W.R.T. THE ABSORBER
  //     COPPER PART - UPPER PART
    
  dgh[0] = 46.;   
  dgh[1] = 74.5;
  dgh[2] = 0.5;
  dgh[3] = 12.;
  dgh[4] = 168.;  
  if (suppmat == 0) {
     gMC->Gsvolu("ICU7", "TUBS", idtmed[279], dgh, 5);   // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU7", "TUBS", idtmed[285], dgh, 5);   // aluminum
  } else {
     gMC->Gsvolu("ICU7", "TUBS", idtmed[274], dgh, 5);   // carbon
  }
  gMC->Gspos("ICU7", 1, "ITSV", 0., 0., -(ztpc+1.5+dgh[2]), 0, "ONLY");  
  
  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON OTHER SIDE W.R.T. THE ABSORBER
  //     COPPER PART - LOWER PART
    
  dgh[0] = 46.; 
  dgh[1] = 74.5;
  dgh[2] = 0.5;
  dgh[3] = 192.;
  dgh[4] = 348.;   
  if (suppmat == 0) {
     gMC->Gsvolu("ICU8", "TUBS", idtmed[279], dgh, 5);   // copper
  } else if (suppmat == 1) {
     gMC->Gsvolu("ICU8", "TUBS", idtmed[285], dgh, 5);   // aluminum
  } else {
     gMC->Gsvolu("ICU8", "TUBS", idtmed[274], dgh, 5);   // carbon
  }
  gMC->Gspos("ICU8", 1, "ITSV", 0., 0., -(ztpc+1.5+dgh[2]), 0, "ONLY");      
    
  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON OTHER SIDE W.R.T. THE ABSORBER
  //     CARBON PART - UPPER PART
  
  dgh[0] = 46.+1.0;  
  dgh[1] = 74.5;
  dgh[2] = 0.75;
  dgh[3] = 12.;
  dgh[4] = 168.;   
  gMC->Gsvolu("ICC7", "TUBS", idtmed[274], dgh, 5);   
  gMC->Gspos("ICC7", 1, "ITSV", 0., 0., -(ztpc+dgh[2]), 0, "ONLY"); 
  
  // --- DEFINE CABLES/COOLING BEHIND THE TPC ON OTHER SIDE W.R.T. THE ABSORBER
  //     CARBON PART - LOWER PART
  
  dgh[0] = 46.+1.0;  
  dgh[1] = 74.5;
  dgh[2] = 0.75;
  dgh[3] = 192.;
  dgh[4] = 348.;     
  gMC->Gsvolu("ICC8", "TUBS", idtmed[274], dgh, 5);   
  gMC->Gspos("ICC8", 1, "ITSV", 0., 0., -(ztpc+dgh[2]), 0, "ONLY");        
    
  // --- DEFINE HOOK TO THE TPC ON OTHER SIDE W.R.T. THE ABSORBER - UPPER PART
  
  dgh[0] = 74.5;
  dgh[1] = 79.5;
  dgh[2] = 2.5;
  dgh[3] = 12.;
  dgh[4] = 168.;    
  gMC->Gsvolu("IHK1", "TUBS", idtmed[284], dgh, 5);  
  gMC->Gspos("IHK1", 1, "ITSV", 0., 0., -ztpc-dgh[2], 0, "ONLY");   
  
  // --- DEFINE HOOK TO THE TPC ON OTHER SIDE W.R.T. THE ABSORBER - LOWER PART
  
  dgh[0] = 74.5;
  dgh[1] = 79.5;
  dgh[2] = 2.5;
  dgh[3] = 192.;
  dgh[4] = 348.;    
  gMC->Gsvolu("IHK2", "TUBS", idtmed[284], dgh, 5);  
  gMC->Gspos("IHK2", 1, "ITSV", 0., 0., -ztpc-dgh[2], 0, "ONLY");      
  
  // --- DEFINE RAILS BETWEEN THE ITS AND THE TPC
  
  if (rails == 1) {
  
     dgh[0] = 2.;          
     dgh[1] = 8.;           
     dgh[2] = 190.;         
     gMC->Gsvolu("IRA1", "BOX ", idtmed[239], dgh, 3);
     gMC->Gspos("IRA1", 1, "ITSV", 53.5, 0., -69.5, 0, "ONLY");   
     gMC->Gsvolu("IRA2", "BOX ", idtmed[239], dgh, 3);    
     gMC->Gspos("IRA2", 1, "ITSV", -53.5, 0., -69.5, 0, "ONLY");    

     dgh[0] = 2.-0.5;    // 0.5 was determined in such a way that the aluminum area is 20.9 cm^2      
     dgh[1] = 8.-0.5;    // 0.5 was determined in such a way that the aluminum area is 20.9 cm^2       
     dgh[2] = 190.;         
     gMC->Gsvolu("IRA3", "BOX ", idtmed[275], dgh, 3);   
     gMC->Gspos("IRA3", 1, "IRA1", 0., 0., 0., 0, "ONLY");   
     gMC->Gsvolu("IRA4", "BOX ", idtmed[275], dgh, 3);     
     gMC->Gspos("IRA4", 1, "IRA2", 0., 0., 0., 0, "ONLY");    

  }

  // --- DEFINE CYLINDERS HOLDING RAILS BETWEEN THE ITS AND THE TPC
  
  dgh[0] = 56.9;    
  dgh[1] = 59.;
  dgh[2] = 0.6;    
  gMC->Gsvolu("ICYL", "TUBE", idtmed[285], dgh, 3);   
  gMC->Gspos("ICYL", 1, "ALIC", 0., 0., 74.1, 0, "ONLY");       
  gMC->Gspos("ICYL", 2, "ALIC", 0., 0., -74.1, idrotm[200], "ONLY");  

  // --- DEFINE SUPPORTS FOR RAILS ATTACHED TO THE CYLINDERS

  dgh[0] = 0.;        
  dgh[1] = 3.;         
  dgh[2] = 5.;  // 5. comes from the fact that the volume has to be 567.6/2 cm^3       
  gMC->Gsvolu("ISR1", "TUBE", idtmed[286], dgh, 3);   
  gMC->Gspos("ISR1", 1, "ITSV", 53.4292, 10.7053, 79.75, 0, "ONLY");    
  gMC->Gspos("ISR1", 2, "ITSV", 53.4292, -10.7053, 79.75, 0, "ONLY");   
  gMC->Gspos("ISR1", 3, "ITSV", -53.4292, 10.7053, 79.75, 0, "ONLY"); 
  gMC->Gspos("ISR1", 4, "ITSV", -53.4292, -10.7053, 79.75, 0, "ONLY");  
  gMC->Gspos("ISR1", 5, "ITSV", 53.4292, 10.7053, -79.75, 0, "ONLY");   
  gMC->Gspos("ISR1", 6, "ITSV", 53.4292, -10.7053, -79.75, 0, "ONLY");   
  gMC->Gspos("ISR1", 7, "ITSV", -53.4292, 10.7053, -79.75, 0, "ONLY"); 
  gMC->Gspos("ISR1", 8, "ITSV", -53.4292, -10.7053, -79.75, 0, "ONLY");         
  
  // --- DEFINE SUPPORTS FOR RAILS ATTACHED TO THE ABSORBER

  dgh[0] = 5.;        
  dgh[1] = 12.;         
  dgh[2] = 5.;         
  gMC->Gsvolu("ISR2", "BOX ", idtmed[285], dgh, 3);   
  gMC->Gspos("ISR2", 1, "ALIC", 53.5, 0., 125.5, 0, "ONLY");
  gMC->Gsvolu("ISR3", "BOX ", idtmed[285], dgh, 3);   
  gMC->Gspos("ISR3", 1, "ALIC", -53.5, 0., 125.5, 0, "ONLY");  
  
  dgh[0] = 5.-2.;        
  dgh[1] = 12.-2.;         
  dgh[2] = 5.;         
  gMC->Gsvolu("ISR4", "BOX ", idtmed[275], dgh, 3);   
  gMC->Gspos("ISR4", 1, "ISR2", 0., 0., 0., 0, "ONLY");     
  gMC->Gsvolu("ISR5", "BOX ", idtmed[275], dgh, 3);   
  gMC->Gspos("ISR5", 1, "ISR3", 0., 0., 0., 0, "ONLY");
  
  // --- DEFINE SUPPORTS TO ATTACH THE ITS TO THE TPC
  
  dgh[0] = 0.;        
  dgh[1] = 5.;         
  dgh[2] = 2.;         
  gMC->Gsvolu("ISR6", "TUBE", idtmed[285], dgh, 3);   
  gMC->Gspos("ISR6", 1, "ALIC", 0., 54., 77., 0, "ONLY"); 
  gMC->Gspos("ISR6", 2, "ALIC", 0., 54., -77., 0, "ONLY"); 
  gMC->Gspos("ISR6", 3, "ALIC", 0., -54., -77., 0, "ONLY");                   
    
    
  // --- Outputs the geometry tree in the EUCLID/CAD format 
  
  if (fEuclidOut) {
    gMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
  }
}
//_____________________________________________________________________________
void AliITSvPPRcoarseasymm::CreateMaterials(){
////////////////////////////////////////////////////////////////////////
  //
  // Create ITS materials
  //     This function defines the default materials used in the Geant
  // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
  // AliITSvPPRcoarseasymm.
  // In general it is automatically replaced by
  // the CreateMaterials routine defined in AliITSv?. Should the function
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
  
  AliMaterial(0, "SPD Si$",      28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(1, "SPD Si chip$", 28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(2, "SPD Si bus$",  28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(3, "SPD C$",       12.011,   6., 2.265,18.8, 999.);
  // v. dens 
  AliMaterial(4, "SPD Air$",    14.61, 7.3, .001205, 30423., 999.);
  AliMaterial(5, "SPD Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(6, "SPD Al$",     26.981539, 13., 2.6989, 8.9, 999.);
  AliMixture( 7, "SPD Water $", awat, zwat, denswat, -2, wwat);
  AliMixture( 8, "SPD Freon$",  afre, zfre, densfre, -2, wfre);
  AliMaterial(9, "SPD End ladder$", 55.845, 26., 7.87/10., 1.76*10., 999.); 
  //AliMaterial(9, "SPD End ladder$", 55.845, 26., -7.87/10., -1.76*10., 999.);   
  AliMaterial(10, "SPD cone$",28.0855, 14., 2.33, 9.36, 999.);       // check !!!!
  // ** 
  AliMedium(0, "SPD Si$",        0, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(1, "SPD Si chip$",   1, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(2, "SPD Si bus$",    2, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(3, "SPD C$",         3, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(4, "SPD Air$",       4, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(5, "SPD Vacuum$",    5, 0,isxfld,sxmgmx, 10.,1.00, .1, .100,10.00);
  AliMedium(6, "SPD Al$",        6, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(7, "SPD Water $",    7, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(8, "SPD Freon$",     8, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(9, "SPD End ladder$",9, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(10, "SPD cone$",    10, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);   
  
  //  225-249 --> Silicon Drift Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(25, "SDD Si$",      28.0855, 14., 2.33,  9.36, 999.);
  AliMaterial(26, "SDD Si chip$", 28.0855, 14., 2.33,  9.36, 999.);
  AliMaterial(27, "SDD Si bus$",  28.0855, 14., 2.33,  9.36, 999.);
  AliMaterial(28, "SDD C$",       12.011,   6., 2.265,18.8,  999.);
  // v. dens 
  AliMaterial(29, "SDD Air$",     14.61, 7.3, .001205, 30423., 999.);
  AliMaterial(30, "SDD Vacuum$",  1e-16, 1e-16, 1e-16, 1e16,  1e16);
  AliMaterial(31, "SDD Al$",      26.981539, 13., 2.6989, 8.9, 999.);
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
  AliMaterial(35, "SDD Copper$", 63.546, 29., 8.96, 1.43, 999.);
  AliMixture( 36, "SDD Ceramics$", acer, zcer, denscer, -5, wcer);
  AliMaterial(37, "SDD Kapton$", 12.011, 6., 1.3, 31.27, 999.);
  AliMaterial(38, "SDD End ladder$", 69.9298, 29.8246, 0.3824, 36.5103, 999.); 
  AliMaterial(39, "SDD cone$",63.546, 29., 1.15, 1.265, 999.);   
  AliMaterial(40, "SDD M55J$",12.3565, 6.4561, 1.8097, 22.9570, 999.);         
  //AliMaterial(38, "SDD End ladder$", 69.9298, 29.8246, -0.3824, -36.5103, 999.); 
  //AliMaterial(39, "SDD cone$",63.546, 29., -1.15, -1.265, 999.);       

  // ** 
  // check A and Z 
  AliMedium(25, "SDD Si$",        25, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(26, "SDD Si chip$",   26, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(27, "SDD Si bus$",    27, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(28, "SDD C$",         28, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(29, "SDD Air$",       29, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(30, "SDD Vacuum$",    30, 0,isxfld,sxmgmx, 10.,1.00, .1, .100,10.00);
  AliMedium(31, "SDD Al$",        31, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(32, "SDD Water $",    32, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(33, "SDD Freon$",     33, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(34, "SDD PCB$",       34, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(35, "SDD Copper$",    35, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(36, "SDD Ceramics$",  36, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(37, "SDD Kapton$",    37, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(38, "SDD End ladder$",38, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(39, "SDD cone$",      39, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(40, "SDD M55J$",      40, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);  
  //  250-274 --> Silicon Strip Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(50, "SSD Si$",      28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(51, "SSD Si chip$", 28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(52, "SSD Si bus$",  28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(53, "SSD C$",       12.011,   6., 2.265,18.8, 999.);
  // v. dens 
  AliMaterial(54, "SSD Air$",     14.61, 7.3, .001205, 30423., 999.);
  AliMaterial(55, "SSD Vacuum$",  1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(56, "SSD Al$",      26.981539, 13., 2.6989, 8.9, 999.);
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
  AliMixture(61, "SSD Ceramics$", acer, zcer, denscer, 5, wcer);
  AliMaterial(62, "SSD Kapton$", 12.011, 6., 1.3, 31.27, 999.);
  // check A and Z 
  AliMaterial(63, "SSD G10FR4$", 17.749, 8.875, 1.8, 21.822, 999.);
  AliMaterial(64, "SSD End ladder$", 32.0988, 15.4021, 0.68, 35.3238, 999.); 
  AliMaterial(65, "SSD cone$",63.546, 29., 1.15, 1.265, 999.);  
  //AliMaterial(64, "SSD End ladder$", 32.0988, 15.4021, -0.68, -35.3238, 999.); 
  //AliMaterial(65, "SSD cone$",63.546, 29., -1.15, -1.265, 999.);    
  // ** 
  AliMedium(50, "SSD Si$",        50, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(51, "SSD Si chip$",   51, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(52, "SSD Si bus$",    52, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(53, "SSD C$",         53, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(54, "SSD Air$",       54, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(55, "SSD Vacuum$",    55, 0,isxfld,sxmgmx, 10.,1.00, .1, .100,10.00);
  AliMedium(56, "SSD Al$",        56, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(57, "SSD Water $",    57, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(58, "SSD Freon$",     58, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(59, "SSD PCB$",       59, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(60, "SSD Copper$",    60, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(61, "SSD Ceramics$",  61, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(62, "SSD Kapton$",    62, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(63, "SSD G10FR4$",    63, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(64, "SSD End ladder$",64, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(65, "SSD cone$",      65, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);

  //     275-299 --> General (end-caps, frames, cooling, cables, etc.) 
  
  AliMaterial(75, "GEN C$", 12.011, 6., 2.265, 18.8, 999.);
  // verify density 
  AliMaterial(76, "GEN Air$", 14.61, 7.3, .001205, 30423., 999.);
  AliMaterial(77, "GEN Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMixture( 78, "GEN POLYETHYL$", apoly, zpoly, .95, -2, wpoly);
  AliMixture( 79, "GEN SERVICES$",  aserv, zserv, 4.68, 4, wserv);
  AliMaterial(80, "GEN Copper$", 63.546, 29., 8.96, 1.43, 999.);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(81, "GEN Water $", awat, zwat, denswat, 2, wwat);
//  AliMaterial(82, "GEN Cables$", 12.011, 6., 2.265, 18.8, 999.);  // check !!!
//  AliMaterial(83, "GEN patch pan$", 12.011, 6., 2.265, 18.8, 999.);  // check !!!  
//  AliMaterial(84, "GEN serv$", 12.011, 6., 2.265, 18.8, 999.);  // check !!!  
  AliMixture(85, "GEN Inox$", asteel, zsteel, 7.88, 4, wsteel);
  AliMaterial(86, "GEN Al$",      26.981539, 13., 2.6989, 8.9, 999.);
  AliMaterial(87,"inox/alum$",    32.1502,15.3383,3.0705,6.9197,999.);
  // ** 
  AliMedium(75,"GEN C$",        75, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(76,"GEN Air$",      76, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(77,"GEN Vacuum$",   77, 0,isxfld,sxmgmx, 10., .10, .1, .100,10.00);
  AliMedium(78,"GEN POLYETHYL$",78, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(79,"GEN SERVICES$", 79, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(80,"GEN Copper$",   80, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(81,"GEN Water $",   81, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
//  AliMedium(82,"GEN Cables$",   82, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
//  AliMedium(83,"GEN patch pan$",83, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);  
//  AliMedium(84,"GEN serv$",     84, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(85,"GEN Inox$",     85, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(86, "GEN Al$",      86, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(87,"inox/alum$",    87, 0,isxfld,sxmgmx, 10., .01, .1, .003, .003);

}
//_____________________________________________________________________________
void AliITSvPPRcoarseasymm::Init(){
////////////////////////////////////////////////////////////////////////
//     Initialise the ITS after it has been created.
////////////////////////////////////////////////////////////////////////
    Int_t i;

    cout << endl;
    for(i=0;i<24;i++) cout << "*";cout << " ITSvPPRcoarseasymm_Init ";
    for(i=0;i<25;i++) cout << "*";cout << endl;
//
    AliITS::Init();
//
    for(i=0;i<72;i++) cout << "*";
    cout << endl;
}  
 
//_____________________________________________________________________________
void AliITSvPPRcoarseasymm::DrawModule() const {
////////////////////////////////////////////////////////////////////////
//     Draw a shaded view of the FMD version 6.
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
void AliITSvPPRcoarseasymm::StepManager(){
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
// no hits for this coarse asymmetric version.
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
//    new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
//
// no hits for this coarse asymmetric version.
//
*/
}

