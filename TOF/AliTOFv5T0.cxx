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
Revision 1.13  2006/10/12 16:35:43  arcelli
definition of the alignable volumes symbolic names added

Revision 1.12  2006/08/22 13:34:46  arcelli
removal of effective c++ warnings (C.Zampolli)

Revision 1.11  2006/07/12 16:03:44  arcelli
updates to match the new numbering of the TOF/TRD mother volumes in FRAME (ALICE convention)

Revision 1.10  2006/05/10 18:40:17  hristov
Larger strings for the names

Revision 1.9  2006/05/04 19:41:42  hristov
Possibility for partial TOF geometry (S.Arcelli)

Revision 1.8  2006/04/20 22:30:50  hristov
Coding conventions (Annalisa)

Revision 1.7  2006/04/16 22:29:05  hristov
Coding conventions (Annalisa)

Revision 1.6  2006/03/20 08:20:35  decaro
Al layer: positioning correction

Revision 1.5  2006/03/20 07:54:20  decaro
Correction of some layer thickness

Revision 1.4  2006/03/13 12:35:44  decaro
Suppression of fractional Z warning

Revision 1.3  2006/02/28 10:38:00  decaro
AliTOFGeometry::fAngles, AliTOFGeometry::fHeights,
AliTOFGeometry::fDistances arrays: dimension definition in the right
location

Revision 1.2  2006/02/27 18:12:14  decaro
Remove in StepManager the dependence of hit indexes from parametrized
TOF position

Revision 1.1  2005/12/15 08:55:33  decaro
New TOF geometry description (V5) -G. Cara Romeo and A. De Caro


Revision 0.1 2004 November G. Cara Romeo and A. De Caro
        Implement new TOF geometry version
	in order to
	   suppress few volume overlaps
	        (in the 4th TOF geometry version),
	   insert the realistic strip numbers and positions

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the functions for version 5 of the Time Of Flight    //
//  detector.                                                                //
//                                                                           //
//  VERSION WITH 5 MODULES AND TILTED STRIPS                                 //
//                                                                           //
//  FULL COVERAGE VERSION + OPTION for PHOS holes                            //
//                                                                           //
//                                                                           //
//Begin_Html                                                                 //
/*                                                                           //
<img src="picts/AliTOFv5T0Class.gif">                                        //
*/                                                                           //
//End_Html                                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TBRIK.h"
#include "TGeometry.h"
#include "TLorentzVector.h"
#include "TNode.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"

#include "AliConst.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMC.h"
#include "AliRun.h"

#include "AliTOFGeometry.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFv5T0.h"

extern TDirectory *gDirectory;
extern TVirtualMC *gMC;

extern AliRun *gAlice;

ClassImp(AliTOFv5T0)

//_____________________________________________________________________________
  AliTOFv5T0::AliTOFv5T0():
  fIdFTOA(-1),
  fIdFTOB(-1),
  fIdFTOC(-1),
  fIdFLTA(-1),
  fIdFLTB(-1),
  fIdFLTC(-1),
  fTOFHoles(kFALSE)
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv5T0::AliTOFv5T0(const char *name, const char *title):
  AliTOF(name,title,"tzero"),
  fIdFTOA(-1),
  fIdFTOB(-1),
  fIdFTOC(-1),
  fIdFLTA(-1),
  fIdFLTB(-1),
  fIdFLTC(-1),
  fTOFHoles(kFALSE)
{
  //
  // Standard constructor
  //
  //
  // Check that FRAME is there otherwise we have no place where to
  // put TOF


  AliModule* frame = (AliModule*)gAlice->GetModule("FRAME");
  if(!frame) {
    AliFatal("TOF needs FRAME to be present");
  } else{
    
    if (fTOFGeometry) delete fTOFGeometry;
    fTOFGeometry = new AliTOFGeometryV5();

    if(frame->IsVersion()==1) {
      AliDebug(1,Form("Frame version %d", frame->IsVersion())); 
      AliDebug(1,"Full Coverage for TOF");
      fTOFHoles=false;}    
    else {
      AliDebug(1,Form("Frame version %d", frame->IsVersion())); 
      AliDebug(1,"TOF with Holes for PHOS");
      fTOFHoles=true;}      
  }
  fTOFGeometry->SetHoles(fTOFHoles);

  //AliTOF::fTOFGeometry = fTOFGeometry;

  // Save the geometry
  TDirectory* saveDir = gDirectory;
  gAlice->GetRunLoader()->CdGAFile();
  fTOFGeometry->Write("TOFgeometry");
  saveDir->cd();

} 

//_____________________________________________________________________________
void AliTOFv5T0::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  //

  TString volPath;
  TString symName;

  TString vpL0  = "ALIC_1/B077_1/BSEGMO";
  TString vpL1 = "_1/BTOF";
  TString vpL2 = "_1";
  TString vpL3 = "/FTOA_0";
  TString vpL4 = "/FLTA_0/FSTR_";

  TString snSM  = "TOF/sm";
  TString snSTRIP = "/strip";

  Int_t nSectors=fTOFGeometry->NSectors();
  Int_t nStrips =fTOFGeometry->NStripA()+
                 2*fTOFGeometry->NStripB()+
                 2*fTOFGeometry->NStripC();

  //
  // The TOF MRPC Strips
  // The symbolic names are: TOF/sm00/strip01
  //                           ...
  //                         TOF/sm17/strip91
 
  Int_t imod=0;

  for (Int_t isect = 0; isect < nSectors; isect++) {
    for (Int_t istr = 1; istr <= nStrips; istr++) {
      
      volPath  = vpL0;
      volPath += isect;
      volPath += vpL1;
      volPath += isect;
      volPath += vpL2;
      volPath += vpL3;
      volPath += vpL4;
      volPath += istr;

      
      symName  = snSM;
      symName += Form("%02d",isect);
      symName += snSTRIP;
      symName += Form("%02d",istr);
            
      AliDebug(2,"--------------------------------------------"); 
      AliDebug(2,Form("Alignable object %d", imod)); 
      AliDebug(2,Form("volPath=%s\n",volPath.Data()));
      AliDebug(2,Form("symName=%s\n",symName.Data()));
      AliDebug(2,"--------------------------------------------"); 
	      
      gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data());
      imod++;
    }
  }


  //
  // The TOF supermodules
  // The symbolic names are: TOF/sm00
  //                           ...
  //                         TOF/sm17
  //
  for (Int_t isect = 0; isect < nSectors; isect++) {

    volPath  = vpL0;
    volPath += isect;
    volPath += vpL1;
    volPath += isect;
    volPath += vpL2;
    volPath += vpL3;

    symName  = snSM;
    symName += Form("%02d",isect);

      AliDebug(2,"--------------------------------------------"); 
      AliDebug(2,Form("Alignable object %d", isect+imod)); 
      AliDebug(2,Form("volPath=%s\n",volPath.Data()));
      AliDebug(2,Form("symName=%s\n",symName.Data()));
      AliDebug(2,"--------------------------------------------"); 
	      
    gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data());

  }
  
}
//____________________________________________________________________________
void AliTOFv5T0::BuildGeometry()
{
  //
  // Build TOF ROOT geometry for the ALICE event display
  //
  TNode *node, *top;
  const int kColorTOF  = 27;
  
  TGeometry *globalGeometry = (TGeometry*)gAlice->GetGeometry();

  // Find top TNODE
  top = globalGeometry->GetNode("alice");
  
  // Position the different copies
  const Float_t krTof  =(fTOFGeometry->Rmax()+fTOFGeometry->Rmin())/2.;
  const Float_t khTof  = fTOFGeometry->Rmax()-fTOFGeometry->Rmin();
  const Int_t   kNTof  = fTOFGeometry->NSectors();
  const Float_t kangle = k2PI/kNTof;

  const Float_t kInterCentrModBorder1 = 49.5;
  const Float_t kInterCentrModBorder2 = 57.5;

  Float_t ang;
  
  // define offset for nodes
  Float_t zOffsetB = (fTOFGeometry->ZlenA()*0.5 + (kInterCentrModBorder1+kInterCentrModBorder2)*0.5)*0.5;
  Float_t zOffsetA = 0.;
  // Define TOF basic volume
  
  char nodeName0[16], nodeName1[16], nodeName2[16];
  char nodeName3[16], nodeName4[16], rotMatNum[16];

  if (fTOFHoles) {
    new TBRIK("S_TOF_B","TOF box","void",
	      fTOFGeometry->StripLength()*0.5, khTof*0.5, fTOFGeometry->ZlenB()*0.5);
    new TBRIK("S_TOF_C","TOF box","void",
	      fTOFGeometry->StripLength()*0.5, khTof*0.5, fTOFGeometry->ZlenB()*0.5);
  }
  new TBRIK("S_TOF_A","TOF box","void",
            fTOFGeometry->StripLength()*0.5, khTof*0.5, fTOFGeometry->ZlenA()*0.5);
  
  for (Int_t nodeNum=1;nodeNum<19;nodeNum++){
    
    if (nodeNum<10) {
      sprintf(rotMatNum,"rot50%i",nodeNum);
      sprintf(nodeName0,"FTO00%i",nodeNum);
      sprintf(nodeName1,"FTO10%i",nodeNum);
      sprintf(nodeName2,"FTO20%i",nodeNum);
      sprintf(nodeName3,"FTO30%i",nodeNum);
      sprintf(nodeName4,"FTO40%i",nodeNum);
    }
    if (nodeNum>9) {
      sprintf(rotMatNum,"rot5%i",nodeNum);
      sprintf(nodeName0,"FTO0%i",nodeNum);
      sprintf(nodeName1,"FTO1%i",nodeNum);
      sprintf(nodeName2,"FTO2%i",nodeNum);
      sprintf(nodeName3,"FTO3%i",nodeNum);
      sprintf(nodeName4,"FTO4%i",nodeNum);
    }
    
    new TRotMatrix(rotMatNum,rotMatNum,90,-20*nodeNum,90,90-20*nodeNum,0,0);
    ang = (4.5-nodeNum) * kangle;

    if (fTOFHoles) {   
      top->cd();
      node = new TNode(nodeName2,nodeName2,"S_TOF_B", krTof*TMath::Cos(ang), krTof*TMath::Sin(ang), zOffsetB,rotMatNum);
      node->SetLineColor(kColorTOF);
      fNodes->Add(node);
      
      top->cd();
      node = new TNode(nodeName3,nodeName3,"S_TOF_C", krTof*TMath::Cos(ang), krTof*TMath::Sin(ang),-zOffsetB,rotMatNum);
      node->SetLineColor(kColorTOF);
      fNodes->Add(node);
    }

    top->cd();
    node = new TNode(nodeName4,nodeName4,"S_TOF_A", krTof*TMath::Cos(ang), krTof*TMath::Sin(ang), zOffsetA,rotMatNum);
    node->SetLineColor(kColorTOF);
    fNodes->Add(node);
  } // end loop on nodeNum

}

//_____________________________________________________________________________
void AliTOFv5T0::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv5T0.gif">
  */
  //End_Html
  //
  // Creates common geometry
  //
  AliTOF::CreateGeometry();
}
 

//_____________________________________________________________________________
void AliTOFv5T0::TOFpc(Float_t xtof,  Float_t ytof, Float_t zlenA,
		       Float_t zlenB)
{
  //
  // Definition of the Time Of Fligh Resistive Plate Chambers
  //

  const Float_t kPi = TMath::Pi();

  const Float_t kInterCentrModBorder1 = 49.5;
  const Float_t kInterCentrModBorder2 = 57.5;
  const Float_t kExterInterModBorder1 = 196.0;
  const Float_t kExterInterModBorder2 = 203.5;

  const Float_t kLengthExInModBorder  = 4.7;
  const Float_t kLengthInCeModBorder  = 7.0;

  // Al layers over all internal module walls (cm)
  const Float_t khAlWall = 0.03;

  // module wall thickness (cm)
  const Float_t kModuleWallThickness = 0.3;

  // Al honeycomb layer between strips and cards (cm)
  const Float_t kHoneycombLayerThickness = 1.5;

  AliDebug(2,Form("zlenA*0.5 = %d", zlenA*0.5));
  AliDebug(1, "************************* TOF geometry **************************");
  
  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xFLT, yFLT, zFLT - sizes of TOF modules (large)
  
  Float_t  xcoor, ycoor, zcoor;
  Float_t  par[3];
  Int_t    *idtmed = fIdtmed->GetArray()-499;
  Int_t    idrotm[100];

  par[0] =  xtof * 0.5;
  par[1] =  ytof * 0.5;
  par[2] = zlenA * 0.5;
  gMC->Gsvolu("FTOA", "BOX ", idtmed[503], par, 3);  // fibre glass
  
  if (fTOFHoles) {
    par[0] =  xtof * 0.5;
    par[1] =  ytof * 0.5;
    par[2] = (zlenA*0.5 - kInterCentrModBorder1)*0.5;
    gMC->Gsvolu("FTOB", "BOX ", idtmed[503], par, 3);  // fibre glass
    gMC->Gsvolu("FTOC", "BOX ", idtmed[503], par, 3);  // fibre glass
  }

  // Positioning of fibre glass modules (FTOA, FTOB and FTOC)
  
  //AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,-90.);
  AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,270.);

  xcoor = 0.;
  ycoor = 0.;
  zcoor = 0.;
  for(Int_t isec=0;isec<18;isec++){
    if(fTOFSectors[isec]==-1)continue;
    char name[16];
    sprintf(name, "BTOF%d",isec);
    if (fTOFHoles && (isec==11||isec==12)) {
    //    if (fTOFHoles && (isec==16||isec==17)) { \\Old 6h convention
      xcoor = 0.;
      ycoor = (zlenA*0.5 + kInterCentrModBorder1)*0.5;
      zcoor = 0.;
      gMC->Gspos("FTOB", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
      gMC->Gspos("FTOC", 0, name, xcoor,-ycoor, zcoor, idrotm[0], "ONLY");
    }
    else gMC->Gspos("FTOA", 0,name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  }
  // Large not sensitive volumes with Insensitive Freon (FLTA, FLTB and FLTC)
  
  Float_t xFLT, yFLT, zFLTA;
  
  xFLT  = xtof  - kModuleWallThickness*2.;
  yFLT  = ytof  - kModuleWallThickness*2.;
  zFLTA = zlenA - kModuleWallThickness*2.;
  
  par[0] = xFLT*0.5;
  par[1] = yFLT*0.5;
  par[2] = zFLTA*0.5;
  gMC->Gsvolu("FLTA", "BOX ", idtmed[507], par, 3); //  Freon mix

  xcoor = 0.;
  ycoor = 0.;
  zcoor = 0.;
  gMC->Gspos ("FLTA", 0, "FTOA", xcoor, ycoor, zcoor, 0, "ONLY");

  if (fTOFHoles) {
    par[0] = xFLT*0.5;
    par[1] = yFLT*0.5;
    par[2] = (zlenA*0.5 - kInterCentrModBorder1-kModuleWallThickness)*0.5;
    gMC->Gsvolu("FLTB", "BOX ", idtmed[507], par, 3); // Freon mix
    gMC->Gsvolu("FLTC", "BOX ", idtmed[507], par, 3); // Freon mix

    xcoor = 0.;
    ycoor = 0.;
    zcoor = kModuleWallThickness*0.5;
    gMC->Gspos ("FLTB", 0, "FTOB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos ("FLTC", 0, "FTOC", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  // Layer of Aluminum before detector (FALA, FALB and FALC)

  par[0] = xFLT*0.5;
  par[1] = khAlWall*0.5;
  par[2] = kInterCentrModBorder1 - (kModuleWallThickness + khAlWall);
  gMC->Gsvolu("FALA", "BOX ", idtmed[505], par, 3); // Alluminium

  xcoor = 0.;
  ycoor = (-yFLT + khAlWall)*0.5;
  zcoor = 0.;
  gMC->Gspos ("FALA", 0, "FLTA", xcoor, ycoor, zcoor, 0, "ONLY");

  par[0] = xFLT*0.5;
  par[1] = khAlWall*0.5;
  par[2] = (kExterInterModBorder2 - kInterCentrModBorder1 - 2.*(kModuleWallThickness + khAlWall))*0.5;
  gMC->Gsvolu("FALB", "BOX ", idtmed[505], par, 3); // Alluminium

  xcoor = 0.;
  ycoor = (-yFLT + khAlWall)*0.5;
  zcoor = (kExterInterModBorder2 + kInterCentrModBorder1)*0.5;
  gMC->Gspos ("FALB", 1, "FLTA", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos ("FALB", 2, "FLTA", xcoor, ycoor,-zcoor, 0, "ONLY");

  par[0] = xFLT*0.5;
  par[1] = khAlWall*0.5;
  par[2] = (zlenA*0.5 - kExterInterModBorder2 - 2.*(kModuleWallThickness + khAlWall))*0.5;
  gMC->Gsvolu("FALC", "BOX ", idtmed[505], par, 3); // Alluminium

  xcoor = 0.;
  ycoor = (-yFLT + khAlWall)*0.5;
  zcoor = (kExterInterModBorder2+zlenA*0.5)*0.5;
  gMC->Gspos ("FALC", 1, "FLTA", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos ("FALC", 2, "FLTA", xcoor, ycoor,-zcoor, 0, "ONLY");

  if (fTOFHoles) {
    xcoor = 0.;
    ycoor = (-yFLT + khAlWall)*0.5;
    zcoor = (zlenA*0.5 - kExterInterModBorder2)*0.5 - kModuleWallThickness*0.5;
    gMC->Gspos ("FALB", 1, "FLTB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos ("FALB", 2, "FLTC", xcoor, ycoor,-zcoor, 0, "ONLY");

    xcoor = 0.;
    ycoor = (-yFLT + khAlWall)*0.5;
    zcoor = (kExterInterModBorder2 - kInterCentrModBorder1)*0.5 + kModuleWallThickness*0.5;
    gMC->Gspos ("FALC", 1, "FLTB", xcoor, ycoor,-zcoor, 0, "ONLY");
    gMC->Gspos ("FALC", 2, "FLTC", xcoor, ycoor, zcoor, 0, "ONLY");
  }

  Float_t y0, alpha, tgal, beta, tgbe, trpa[11];

  // Fibre glass walls between central and intermediate modules (FWZ1 and FWZ2; holes -> FZ1B, FZ1C, FZ2B)

  tgal = (yFLT*0.5 - 2.*kLengthInCeModBorder)/(kInterCentrModBorder2 - kInterCentrModBorder1);
  alpha = TMath::ATan(tgal);
  beta = (kPi*0.5 - alpha)*0.5;
  tgbe = TMath::Tan(beta);
  trpa[0]  = xFLT*0.5;
  trpa[1]  = 0.;
  trpa[2]  = 0.;
  trpa[3]  = kModuleWallThickness;
  trpa[4]  = (kLengthInCeModBorder - kModuleWallThickness*tgbe)*0.5;
  trpa[5]  = (kLengthInCeModBorder + kModuleWallThickness*tgbe)*0.5;
  trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  trpa[7]  = kModuleWallThickness;
  trpa[8]  = (kLengthInCeModBorder - kModuleWallThickness*tgbe)*0.5;
  trpa[9]  = (kLengthInCeModBorder + kModuleWallThickness*tgbe)*0.5;
  trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  gMC->Gsvolu("FWZ1","TRAP", idtmed[503], trpa, 11);   // fibre glass

  AliMatrix (idrotm[1],90., 90.,180.,0.,90.,180.);
  AliMatrix (idrotm[4],90., 90.,  0.,0.,90.,  0.);

  xcoor = 0.;
  ycoor = -yFLT*0.5 + kLengthInCeModBorder*0.5;
  zcoor = kInterCentrModBorder1;
  gMC->Gspos("FWZ1", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[1],"ONLY");
  gMC->Gspos("FWZ1", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[4],"ONLY");

  if (fTOFHoles) {
    y0 = kLengthInCeModBorder - kModuleWallThickness*0.5*tgbe;
    trpa[0]  = xFLT*0.5;
    trpa[1]  = 0.;
    trpa[2]  = 0.;
    trpa[3]  = kModuleWallThickness*0.5;
    trpa[4]  = (y0 - kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[5]  = (y0 + kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    trpa[7]  = kModuleWallThickness*0.5;
    trpa[8]  = (y0 - kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[9]  = (y0 + kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    gMC->Gsvolu("FZ1B","TRAP", idtmed[503], trpa, 11);   // fibre glass

    xcoor = 0.;
    ycoor = -yFLT*0.5 + kLengthInCeModBorder*0.5 - kModuleWallThickness*0.25*tgbe;
    zcoor = -kInterCentrModBorder1 + (zlenA*0.5 + kInterCentrModBorder1)*0.5 - kModuleWallThickness;
    gMC->Gspos("FZ1B", 1,"FLTB", xcoor, ycoor, zcoor,idrotm[4],"ONLY");
    gMC->Gspos("FZ1B", 2,"FLTC", xcoor, ycoor,-zcoor,idrotm[1],"ONLY");
  }

  AliMatrix (idrotm[2],90.,270.,  0.,0.,90.,180.);
  AliMatrix (idrotm[5],90.,270.,180.,0.,90.,  0.);

  xcoor = 0.;
  ycoor = -kLengthInCeModBorder*0.5;
  zcoor = kInterCentrModBorder2;
  gMC->Gspos("FWZ1", 3,"FLTA", xcoor, ycoor, zcoor,idrotm[2],"ONLY");
  gMC->Gspos("FWZ1", 4,"FLTA", xcoor, ycoor,-zcoor,idrotm[5],"ONLY");

  if (fTOFHoles) {
    y0 = kLengthInCeModBorder + kModuleWallThickness*0.5*tgbe;
    trpa[0]  = xFLT*0.5;
    trpa[1]  = 0.;
    trpa[2]  = 0.;
    trpa[3]  = kModuleWallThickness*0.5;
    trpa[4]  = (y0 - kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[5]  = (y0 + kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    trpa[7]  = kModuleWallThickness*0.5;
    trpa[8]  = (y0 - kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[9]  = (y0 + kModuleWallThickness*0.5*tgbe)*0.5;
    trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    gMC->Gsvolu("FZ1C","TRAP", idtmed[503], trpa, 11);   // fibre glass

    xcoor = 0.;
    ycoor = -kLengthInCeModBorder*0.5 - kModuleWallThickness*0.25*tgbe;
    zcoor = -kInterCentrModBorder2 + (zlenA*0.5 + kInterCentrModBorder1)*0.5 - kModuleWallThickness;
    gMC->Gspos("FZ1C", 1,"FLTB", xcoor, ycoor, zcoor,idrotm[5],"ONLY");
    gMC->Gspos("FZ1C", 2,"FLTC", xcoor, ycoor,-zcoor,idrotm[2],"ONLY");
  }

  trpa[0] = 0.5*(kInterCentrModBorder2 - kInterCentrModBorder1)/TMath::Cos(alpha);
  trpa[1] = kModuleWallThickness;
  trpa[2] = xFLT*0.5;
  trpa[3] = -beta*kRaddeg;
  trpa[4] = 0.;
  trpa[5] = 0.;
  gMC->Gsvolu("FWZ2","PARA", idtmed[503], trpa, 6);    // fibre glass

  AliMatrix (idrotm[3],     alpha*kRaddeg,90.,90.+alpha*kRaddeg,90.,90.,180.);
  AliMatrix (idrotm[6],180.-alpha*kRaddeg,90.,90.-alpha*kRaddeg,90.,90.,  0.);

  xcoor = 0.;
  ycoor = -yFLT*0.25;
  zcoor = (kInterCentrModBorder2 + kInterCentrModBorder1)*0.5;
  gMC->Gspos("FWZ2", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[3],"ONLY");
  gMC->Gspos("FWZ2", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[6],"ONLY");

  if (fTOFHoles) {
    trpa[0] = 0.5*(kInterCentrModBorder2 - kInterCentrModBorder1)/TMath::Cos(alpha);
    trpa[1] = kModuleWallThickness*0.5;
    trpa[2] = xFLT*0.5;
    trpa[3] = -beta*kRaddeg;
    trpa[4] = 0.;
    trpa[5] = 0.;
    gMC->Gsvolu("FZ2B","PARA", idtmed[503], trpa, 6);    // fibre glass

    xcoor = 0.;
    ycoor = -yFLT*0.25 - kModuleWallThickness*0.5*tgbe;
    zcoor = -(kInterCentrModBorder2 + kInterCentrModBorder1)*0.5 + (zlenA*0.5 + kInterCentrModBorder1)*0.5 - kModuleWallThickness;
    gMC->Gspos("FZ2B", 1,"FLTB", xcoor, ycoor, zcoor,idrotm[6],"ONLY");
    gMC->Gspos("FZ2B", 2,"FLTC", xcoor, ycoor,-zcoor,idrotm[3],"ONLY");
  }

  // Fibre glass walls between intermediate and lateral modules (FWZ3 and FWZ4)

  tgal = (yFLT*0.5 - 2.*kLengthExInModBorder)/(kExterInterModBorder2 - kExterInterModBorder1);
  alpha = TMath::ATan(tgal);
  beta = (kPi*0.5 - alpha)*0.5;
  tgbe = TMath::Tan(beta);
  trpa[0]  = xFLT*0.5;
  trpa[1]  = 0.;
  trpa[2]  = 0.;
  trpa[3]  = kModuleWallThickness;
  trpa[4]  = (kLengthExInModBorder - kModuleWallThickness*tgbe)*0.5;
  trpa[5]  = (kLengthExInModBorder + kModuleWallThickness*tgbe)*0.5;
  trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  trpa[7]  = kModuleWallThickness;
  trpa[8]  = (kLengthExInModBorder - kModuleWallThickness*tgbe)*0.5;
  trpa[9]  = (kLengthExInModBorder + kModuleWallThickness*tgbe)*0.5;
  trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  gMC->Gsvolu("FWZ3","TRAP", idtmed[503], trpa, 11);    // fibre glass

  xcoor = 0.;
  ycoor = -kLengthExInModBorder*0.5;
  zcoor = kExterInterModBorder1;
  gMC->Gspos("FWZ3", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[5],"ONLY");
  gMC->Gspos("FWZ3", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[2],"ONLY");

  if (fTOFHoles) {
    xcoor = 0.;
    ycoor = -kLengthExInModBorder*0.5;
    zcoor = -kExterInterModBorder1 + (zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5;
    gMC->Gspos("FWZ3", 5,"FLTB", xcoor, ycoor, zcoor,idrotm[2],"ONLY");
    gMC->Gspos("FWZ3", 6,"FLTC", xcoor, ycoor,-zcoor,idrotm[5],"ONLY");
  }

  xcoor = 0.;
  ycoor = -yFLT*0.5 + kLengthExInModBorder*0.5;
  zcoor = kExterInterModBorder2;
  gMC->Gspos("FWZ3", 3,"FLTA", xcoor, ycoor, zcoor,idrotm[4],"ONLY");
  gMC->Gspos("FWZ3", 4,"FLTA", xcoor, ycoor,-zcoor,idrotm[1],"ONLY");

  if (fTOFHoles) {
    xcoor = 0.;
    ycoor = -yFLT*0.5 + kLengthExInModBorder*0.5;
    zcoor = -kExterInterModBorder2 + (zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5;
    gMC->Gspos("FWZ3", 7,"FLTB", xcoor, ycoor, zcoor,idrotm[1],"ONLY");
    gMC->Gspos("FWZ3", 8,"FLTC", xcoor, ycoor,-zcoor,idrotm[4],"ONLY");
  }

  trpa[0] = 0.5*(kExterInterModBorder2 - kExterInterModBorder1)/TMath::Cos(alpha);
  trpa[1] = kModuleWallThickness;
  trpa[2] = xFLT*0.5;
  trpa[3] = -beta*kRaddeg;
  trpa[4] = 0.;
  trpa[5] = 0.;
  gMC->Gsvolu("FWZ4","PARA", idtmed[503], trpa, 6);    // fibre glass

  AliMatrix (idrotm[13],alpha*kRaddeg,90.,90.+alpha*kRaddeg,90.,90.,180.);
  AliMatrix (idrotm[16],180.-alpha*kRaddeg,90.,90.-alpha*kRaddeg,90.,90.,0.);

  xcoor = 0.;
  ycoor = -yFLT*0.25;
  zcoor = (kExterInterModBorder2 + kExterInterModBorder1)*0.5;
  gMC->Gspos("FWZ4", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[16],"ONLY");
  gMC->Gspos("FWZ4", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[13],"ONLY");

  if (fTOFHoles) {
    xcoor = 0.;
    ycoor = -yFLT*0.25;
    zcoor = -(kExterInterModBorder2 + kExterInterModBorder1)*0.5 + (zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5;
    gMC->Gspos("FWZ4", 3,"FLTB", xcoor, ycoor, zcoor,idrotm[13],"ONLY");
    gMC->Gspos("FWZ4", 4,"FLTC", xcoor, ycoor,-zcoor,idrotm[16],"ONLY");
  }


  ///////////////// Detector itself //////////////////////

  const Int_t    knx   = fTOFGeometry->NpadX();  // number of pads along x
  const Int_t    knz   = fTOFGeometry->NpadZ();  // number of pads along z
  const Float_t  kPadX = fTOFGeometry->XPad();   // pad length along x
  const Float_t  kPadZ = fTOFGeometry->ZPad();   // pad length along z

  // new description for strip volume -double stack strip-
  // -- all constants are expressed in cm
  // heigth of different layers
  const Float_t khhony   = 1.0    ;   // heigth of HONY  Layer
  const Float_t khpcby   = 0.08   ;   // heigth of PCB   Layer
  const Float_t khrgly   = 0.055  ;   // heigth of RED GLASS  Layer

  const Float_t khfiliy  = 0.125  ;   // heigth of FISHLINE  Layer
  const Float_t khglassy = 0.160*0.5; // heigth of GLASS  Layer
  const Float_t khglfy   = khfiliy+2.*khglassy;// heigth of GLASS+FISHLINE  Layer

  const Float_t khcpcby  = 0.16   ;   // heigth of PCB  Central Layer
  const Float_t kwhonz   = 8.1    ;   // z dimension of HONEY  Layer
  const Float_t kwpcbz1  = 10.6   ;   // z dimension of PCB  Lower Layer
  const Float_t kwpcbz2  = 11.6   ;   // z dimension of PCB  Upper Layer
  const Float_t kwcpcbz  = 12.4   ;   // z dimension of PCB  Central Layer
  const Float_t kwrglz   = 8.     ;   // z dimension of RED GLASS  Layer
  const Float_t kwglfz   = 7.     ;   // z dimension of GLASS+FISHLN Layer
  const Float_t klsensmx = knx*kPadX; // length of Sensitive Layer
  const Float_t khsensmy = 0.05;//0.11;//0.16;// heigth of Sensitive Layer
  const Float_t kwsensmz = knz*kPadZ; // width of Sensitive Layer
  
  // heigth of the FSTR Volume (the strip volume)
  const Float_t khstripy = 2.*khhony+2.*khpcby+4.*khrgly+2.*khglfy+khcpcby;

  // width  of the FSTR Volume (the strip volume)
  const Float_t kwstripz = kwcpcbz;
  // length of the FSTR Volume (the strip volume)
  const Float_t klstripx = fTOFGeometry->StripLength();
  
  Float_t parfp[3]={klstripx*0.5,khstripy*0.5,kwstripz*0.5};
  // Coordinates of the strip center in the strip reference frame;
  // used for positioning internal strip volumes
  Float_t posfp[3]={0.,0.,0.};  

  // FSTR volume definition-filling this volume with non sensitive Gas Mixture
  gMC->Gsvolu("FSTR","BOX",idtmed[507],parfp,3);     // Freon mix

  //-- HONY Layer definition
  parfp[1] = khhony*0.5;
  parfp[2] = kwhonz*0.5;
  gMC->Gsvolu("FHON","BOX",idtmed[501],parfp,3);     // honeycomb (Nomex)
  // positioning 2 HONY Layers on FSTR volume
  posfp[1] =-khstripy*0.5+parfp[1];
  gMC->Gspos("FHON",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FHON",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  
  //-- PCB Layer definition
  parfp[1] = khpcby*0.5;
  parfp[2] = kwpcbz1*0.5;
  gMC->Gsvolu("FPC1","BOX",idtmed[502],parfp,3);     // G10
  parfp[2] = kwpcbz2*0.5;
  gMC->Gsvolu("FPC2","BOX",idtmed[502],parfp,3);     // G10
  // positioning 2 PCB Layers on FSTR volume
  posfp[1] =-khstripy*0.5+khhony+parfp[1];
  gMC->Gspos("FPC1",1,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FPC2",1,"FSTR",0., posfp[1],0.,0,"ONLY");

  //-- central PCB layer definition
  parfp[1] = khcpcby*0.5;
  parfp[2] = kwcpcbz*0.5;
  gMC->Gsvolu("FPCB","BOX",idtmed[502],parfp,3);     // G10
  // positioning the central PCB layer
  gMC->Gspos("FPCB",1,"FSTR",0.,0.,0.,0,"ONLY");
  
  //      Sensitive volume
  Float_t parfs[3] = {klsensmx*0.5, khsensmy*0.5, kwsensmz*0.5};
  gMC->Gsvolu("FSEN","BOX",idtmed[508],parfs,3);     // sensitive ...
  // dividing FSEN along z in knz=2 and along x in knx=48
  gMC->Gsdvn("FSEZ","FSEN",knz,3);
  gMC->Gsdvn("FPAD","FSEZ",knx,1);
  // positioning a Sensitive layer inside FPCB
  gMC->Gspos("FSEN",1,"FPCB",0.,0.,0.,0,"ONLY");

  //-- RED GLASS Layer definition
  parfp[1] = khrgly*0.5;
  parfp[2] = kwrglz*0.5;
  gMC->Gsvolu("FRGL","BOX",idtmed[509],parfp,3);     // glass
  // positioning 4 RED GLASS Layers on FSTR volume
  posfp[1] = -khstripy*0.5+khhony+khpcby+parfp[1];
  gMC->Gspos("FRGL",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FRGL",4,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  posfp[1] = (khcpcby+khrgly)*0.5;
  gMC->Gspos("FRGL",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FRGL",3,"FSTR",0., posfp[1],0.,0,"ONLY");

  //-- GLASS Layer definition
  parfp[1] = khglassy*0.5;
  parfp[2] = kwglfz*0.5;
  gMC->Gsvolu("FGLA","BOX",idtmed[509],parfp,3);     // glass

  // positioning 4 GLASS Layers on FSTR volume
  posfp[1] = -khstripy*0.5+khhony+khpcby+khrgly+parfp[1];
  gMC->Gspos("FGLA",1,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGLA",4,"FSTR",0., posfp[1],0.,0,"ONLY");
  posfp[1] = khcpcby*0.5+khrgly+khglassy*0.5;
  gMC->Gspos("FGLA",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGLA",3,"FSTR",0., posfp[1],0.,0,"ONLY");

  //-- FREON Layer definition
  parfp[1] = khfiliy*0.5;
  gMC->Gsvolu("FFIS","BOX",idtmed[507],parfp,3);     // freon

  // positioning 2 FREON Layers on FSTR volume
  posfp[1] = -khstripy*0.5+khhony+khpcby+khrgly+khglassy+parfp[1];
  gMC->Gspos("FFIS",1,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FFIS",2,"FSTR",0., posfp[1],0.,0,"ONLY");

  /*
  //-- GLASS+FISHLINE Layer definition
  parfp[1] = khglfy*0.5;
  parfp[2] = kwglfz*0.5;
  gMC->Gsvolu("FGLF","BOX",idtmed[504],parfp,3);

  // positioning 2 GLASS+FISHLINE Layers on FSTR volume
  posfp[1] = (khcpcby+khglfy)*0.5+khrgly;
  gMC->Gspos("FGLF",1,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGLF",2,"FSTR",0., posfp[1],0.,0,"ONLY");
  */

  //  Positioning the Strips  (FSTR) in the FLT volumes
  Int_t maxStripNumbers [5] ={fTOFGeometry->NStripC(),
			      fTOFGeometry->NStripB(),
			      fTOFGeometry->NStripA(),
			      fTOFGeometry->NStripB(),
			      fTOFGeometry->NStripC()};

  Int_t totalStrip = 0;
  Float_t xpos, zpos, ypos, ang;
  for(Int_t iplate =0; iplate < fTOFGeometry->NPlates(); iplate++){
    if (iplate>0) totalStrip += maxStripNumbers[iplate-1];
    for(Int_t istrip =0; istrip < maxStripNumbers[iplate]; istrip++){

      ang = fTOFGeometry->GetAngles(iplate,istrip);
      AliDebug(1, Form(" iplate = %1i, istrip = %2i ---> ang = %f", iplate, istrip, ang));
 
      if (ang>0.)       AliMatrix (idrotm[istrip+totalStrip+1],90.,0.,90.+ang,90., ang, 90.);
      else if (ang==0.) AliMatrix (idrotm[istrip+totalStrip+1],90.,0.,90.,90., 0., 0.);
      else if (ang<0.)  AliMatrix (idrotm[istrip+totalStrip+1],90.,0.,90.+ang,90.,-ang,270.);

      xpos = 0.;
      zpos = fTOFGeometry->GetDistances(iplate,istrip);
      ypos = fTOFGeometry->GetHeights(iplate,istrip);

      gMC->Gspos("FSTR",istrip+totalStrip+1,"FLTA", xpos, ypos,-zpos,idrotm[istrip+totalStrip+1],  "ONLY");

      if (fTOFHoles) {
	if (istrip+totalStrip+1>53) gMC->Gspos("FSTR",istrip+totalStrip+1,"FLTC", xpos, ypos,-zpos-(zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5,idrotm[istrip+totalStrip+1],"ONLY");
	if (istrip+totalStrip+1<39) gMC->Gspos("FSTR",istrip+totalStrip+1,"FLTB", xpos, ypos,-zpos+(zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5,idrotm[istrip+totalStrip+1],"ONLY");
      }
    }
  }

  //  1.5 cm Al honeycomb layer between strips and cards
  par[0] = xFLT*0.5;
  par[1] = kHoneycombLayerThickness*0.5;
  par[2] = zFLTA*0.5;
  gMC->Gsvolu("FPEA", "BOX ", idtmed[506], par, 3);   // Al honeycomb

  xcoor = 0.;
  ycoor = kHoneycombLayerThickness*0.5;
  zcoor = 0.;
  gMC->Gspos ("FPEA", 0, "FLTA", xcoor, ycoor, zcoor, 0, "ONLY");

  if (fTOFHoles) {
    par[0] = xFLT*0.5;
    par[1] = kHoneycombLayerThickness*0.5;
    par[2] = (zlenA*0.5 - kInterCentrModBorder2-kModuleWallThickness)*0.5;
    gMC->Gsvolu("FPEB", "BOX ", idtmed[506], par, 3);   // Al honeycomb

    xcoor = 0.;
    ycoor = kHoneycombLayerThickness*0.5;
    zcoor = (kInterCentrModBorder2-kInterCentrModBorder1)*0.5;
    gMC->Gspos ("FPEB", 1, "FLTB", xcoor, ycoor,-zcoor, 0, "ONLY");
    gMC->Gspos ("FPEB", 2, "FLTC", xcoor, ycoor, zcoor, 0, "ONLY");
  }

  // frame of Air
  par[0] = xFLT*0.5;
  par[1] = (yFLT*0.5 - kHoneycombLayerThickness)*0.5;
  par[2] = zFLTA *0.5;
  gMC->Gsvolu("FAIA", "BOX ", idtmed[500], par, 3); // Air

  xcoor = 0.;
  ycoor = kHoneycombLayerThickness + (yFLT*0.5 - kHoneycombLayerThickness)*0.5;
  zcoor = 0.;
  gMC->Gspos ("FAIA", 0, "FLTA", xcoor, ycoor, zcoor, 0, "ONLY");

  if (fTOFHoles) {
    par[0] = xFLT*0.5;
    par[1] = (yFLT*0.5 - kHoneycombLayerThickness)*0.5;
    par[2] = (zlenA*0.5 - kInterCentrModBorder2 - kModuleWallThickness)*0.5;
    gMC->Gsvolu("FAIB", "BOX ", idtmed[500], par, 3); // Air
    gMC->Gsvolu("FAIC", "BOX ", idtmed[500], par, 3); // Air

    xcoor = 0.;
    ycoor = kHoneycombLayerThickness + (yFLT*0.5 - kHoneycombLayerThickness)*0.5;
    zcoor = (kInterCentrModBorder2-kInterCentrModBorder1)*0.5;
    gMC->Gspos ("FAIB", 0, "FLTB", xcoor, ycoor,-zcoor, 0, "ONLY");
    gMC->Gspos ("FAIC", 0, "FLTC", xcoor, ycoor, zcoor, 0, "ONLY");
  }

  // start with cards and cooling tubes
  // finally, cards, cooling tubes and layer for thermal dispersion
  // 3 volumes
  
  // see GEOM200 in GEANT manual

  Float_t cardpar[3];

  // card volume definition
  cardpar[0]= xFLT*0.5;
  cardpar[1]= 5.;
  cardpar[2]= 0.1;
  gMC->Gsvolu("FCAR", "BOX ", idtmed[502], cardpar, 3); // PCB Card 

  //alu plate volume definition
  cardpar[1]= 3.5;
  cardpar[2]= 0.05;
  gMC->Gsvolu("FALP", "BOX ", idtmed[505], cardpar, 3); // Alu Plate

  // tube volume definition
  Float_t tubepar[3];
  tubepar[0]= 0.;
  tubepar[1]= 0.4;
  tubepar[2]= 61.;
  gMC->Gsvolu("FTUB", "TUBE", idtmed[511], tubepar, 3); // cooling tubes (steel)

  //tubepar[0]= 0.;
  tubepar[1]= 0.35;
  //tubepar[2]= 61.;
  gMC->Gsvolu("FITU", "TUBE", idtmed[510], tubepar, 3); // cooling water
  // positioning water tube into the steel one
  gMC->Gspos("FITU",1,"FTUB",0.,0.,0.,0,"ONLY");

  // rotation matrix
  AliMatrix(idrotm[99], 180., 90., 90., 90., 90., 0.);

  // central module positioning
  Float_t cardpos[3], aplpos2;
  Float_t stepforcardA = 6.625;
  Float_t tdis = 0.6;
  Float_t aplpos1 = -2.;

  cardpos[0]= 0.;
  cardpos[1]= -0.5;
  cardpos[2]= -53.;
  //  tubepos= -53.+tdis;
  Int_t icard;
  for (icard=39; icard<54; ++icard) {
    cardpos[2]= cardpos[2]+stepforcardA;
    aplpos2 = cardpos[2]+0.15;
    gMC->Gspos("FCAR",icard,"FAIA",cardpos[0],cardpos[1],     cardpos[2],         0,"ONLY"); 
    gMC->Gspos("FALP",icard,"FAIA",cardpos[0],   aplpos1,        aplpos2,         0,"ONLY");
    gMC->Gspos("FTUB",icard,"FAIA",        0.,cardpos[1],cardpos[2]+tdis,idrotm[99],"ONLY");
  }
  
  // intermediate module positioning
  Float_t stepforcardB= 7.05;
  Float_t offs = 53.;

  cardpos[2]= offs;
  for (icard=20; icard<39; ++icard) {
    cardpos[2]= cardpos[2]+stepforcardB;
    aplpos2 = cardpos[2]+0.15;

    gMC->Gspos("FCAR",icard+34,"FAIA",cardpos[0],cardpos[1],      cardpos[2],         0,"ONLY");
    gMC->Gspos("FALP",icard+34,"FAIA",cardpos[0],   aplpos1,         aplpos2,         0,"ONLY");
    gMC->Gspos("FTUB",icard+34,"FAIA",        0.,cardpos[1], cardpos[2]+tdis,idrotm[99],"ONLY");
    gMC->Gspos("FCAR",58-icard,"FAIA",cardpos[0],cardpos[1],     -cardpos[2],         0,"ONLY");
    gMC->Gspos("FALP",58-icard,"FAIA",cardpos[0],   aplpos1,        -aplpos2,         0,"ONLY");
    gMC->Gspos("FTUB",58-icard,"FAIA",        0.,cardpos[1],-cardpos[2]-tdis,idrotm[99],"ONLY");

    if (fTOFHoles) {
      gMC->Gspos("FCAR",icard+34+182,"FAIC",cardpos[0],cardpos[1],      cardpos[2]-(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FALP",icard+34+182,"FAIC",cardpos[0],   aplpos1,         aplpos2-(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FTUB",icard+34+182,"FAIC",        0.,cardpos[1], cardpos[2]+tdis-(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,idrotm[99],"ONLY");
      gMC->Gspos("FCAR",58-icard+ 91,"FAIB",cardpos[0],cardpos[1],     -cardpos[2]+(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FALP",58-icard+ 91,"FAIB",cardpos[0],   aplpos1,        -aplpos2+(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FTUB",58-icard+ 91,"FAIB",        0.,cardpos[1],-cardpos[2]-tdis+(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,idrotm[99],"ONLY");
    }
    
  }

  // outer module positioning
  Float_t stepforcardC= 8.45238;
  offs += zlenB;
  cardpos[2]= offs;
  for (icard=1; icard<20; ++icard) {
    cardpos[2]= cardpos[2]+stepforcardC;
    aplpos2 = cardpos[2]+0.15;

    gMC->Gspos("FCAR",icard+72,"FAIA",cardpos[0],cardpos[1],      cardpos[2],         0,"ONLY"); 
    gMC->Gspos("FALP",icard+72,"FAIA",cardpos[0],   aplpos1,         aplpos2,         0,"ONLY");
    gMC->Gspos("FTUB",icard+72,"FAIA",        0.,cardpos[1], cardpos[2]+tdis,idrotm[99],"ONLY");
    gMC->Gspos("FCAR",20-icard,"FAIA",cardpos[0],cardpos[1],     -cardpos[2],         0,"ONLY");
    gMC->Gspos("FALP",20-icard,"FAIA",cardpos[0],   aplpos1,        -aplpos2,         0,"ONLY");
    gMC->Gspos("FTUB",20-icard,"FAIA",        0.,cardpos[1],-cardpos[2]-tdis,idrotm[99],"ONLY");

    if (fTOFHoles) {
      gMC->Gspos("FCAR",icard+72+182,"FAIC",cardpos[0],cardpos[1],      cardpos[2]-(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FALP",icard+72+182,"FAIC",cardpos[0],   aplpos1,         aplpos2-(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FTUB",icard+72+182,"FAIC",        0.,cardpos[1], cardpos[2]+tdis-(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,idrotm[99],"ONLY");
      gMC->Gspos("FCAR",20-icard+ 91,"FAIB",cardpos[0],cardpos[1],     -cardpos[2]+(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FALP",20-icard+ 91,"FAIB",cardpos[0],   aplpos1,        -aplpos2+(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,         0,"ONLY");
      gMC->Gspos("FTUB",20-icard+ 91,"FAIB",        0.,cardpos[1],-cardpos[2]-tdis+(zlenA*0.5 + kInterCentrModBorder2 - kModuleWallThickness)*0.5,idrotm[99],"ONLY");
    }
  }

}
//_____________________________________________________________________________
void AliTOFv5T0::DrawModule() const
{
  //
  // Draw a shaded view of the Time Of Flight version 4
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);

  //
  //Set volumes visible
  // 

  //Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN", 0);

//=====> Level 1
  // Level 1 for TOF volumes
  gMC->Gsatt("B077","seen", 0);

//=====> Level 2
  // Level 2 for TOF volumes
  gMC->Gsatt("B071","seen", 0);
  gMC->Gsatt("B074","seen", 0);
  gMC->Gsatt("B075","seen", 0);
  gMC->Gsatt("B076","seen",-1); // all B076 sub-levels skipped -
  gMC->Gsatt("B080","seen", 0);  // B080 does not has sub-level                

  // Level 2 of B071
  gMC->Gsatt("B056","seen", 0);  // B056 does not has sub-levels  -
  gMC->Gsatt("B063","seen",-1); // all B063 sub-levels skipped   -
  gMC->Gsatt("B065","seen",-1); // all B065 sub-levels skipped   -
  gMC->Gsatt("B067","seen",-1); // all B067 sub-levels skipped   -
  gMC->Gsatt("B072","seen",-1); // all B072 sub-levels skipped   -

  gMC->Gsatt("BTR1","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO1","seen", 0);

  // Level 2 of B074
  gMC->Gsatt("BTR2","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO2","seen", 0);

  // Level 2 of B075
  gMC->Gsatt("BTR3","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO3","seen", 0);

  // Level 3 of B071, B074 and B075
  gMC->Gsatt("FTOA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FTOB","SEEN", 0);

  // Level 4 of B071, B074 and B075
  gMC->Gsatt("FLTA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FLTB","SEEN",0);
  if (fTOFHoles) gMC->Gsatt("FLTC","SEEN",0);

  // Level 5 of B071, B074 and B075
  gMC->Gsatt("FAIA","SEEN",-1);  // all FAIA sub-levels skipped   -
  if (fTOFHoles) gMC->Gsatt("FAIB","SEEN",-1);  // all FAIB sub-levels skipped   -
  if (fTOFHoles) gMC->Gsatt("FAIC","SEEN",-1);  // all FAIC sub-levels skipped   -

  gMC->Gsatt("FALA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FALB","SEEN", 0);

  gMC->Gsatt("FPEA","SEEN", 1);
  if (fTOFHoles) gMC->Gsatt("FPEB","SEEN", 1);

  gMC->Gsatt("FSTR","SEEN",-1);  // all FSTR sub-levels skipped   -

  gMC->Gsatt("FWZ1","SEEN", 0);
  gMC->Gsatt("FWZ2","SEEN", 0);
  gMC->Gsatt("FWZ3","SEEN", 0);
  gMC->Gsatt("FWZ4","SEEN", 0);
  if (fTOFHoles) {
    gMC->Gsatt("FZ1B","SEEN", 0);
    gMC->Gsatt("FZ1C","SEEN", 0);
    gMC->Gsatt("FZ2B","SEEN", 0);
  }

  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 100, 1000, 100, 1000, 100, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 10, 9.5, .018, .018);
  gMC->Gdhead(1111, "Time Of Flight");
  gMC->Gdman(18, 3, "MAN");
  gMC->Gdopt("hide","off");
}
//_____________________________________________________________________________
void AliTOFv5T0::DrawDetectorModules() const
{
  //
  // Draw a shaded view of the TOF detector version 4
  //
 
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);

  //
  //Set volumes visible
  // 

  //Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN", 0);

//=====> Level 1
  // Level 1 for TOF volumes
  gMC->Gsatt("B077","seen", 0);

//=====> Level 2
  // Level 2 for TOF volumes
  gMC->Gsatt("B071","seen", 0);
  gMC->Gsatt("B074","seen", 0);
  gMC->Gsatt("B075","seen", 0);
  gMC->Gsatt("B076","seen",-1); // all B076 sub-levels skipped -
  gMC->Gsatt("B080","seen", 0);  // B080 does not has sub-level                

  // Level 2 of B071
  gMC->Gsatt("B056","seen", 0);  // B056 does not has sub-levels  -
  gMC->Gsatt("B063","seen",-1); // all B063 sub-levels skipped   -
  gMC->Gsatt("B065","seen",-1); // all B065 sub-levels skipped   -
  gMC->Gsatt("B067","seen",-1); // all B067 sub-levels skipped   -
  gMC->Gsatt("B072","seen",-1); // all B072 sub-levels skipped   -

  gMC->Gsatt("BTR1","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO1","seen", 0);

  // Level 2 of B074
  gMC->Gsatt("BTR2","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO2","seen", 0);

  // Level 2 of B075
  gMC->Gsatt("BTR3","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO3","seen", 0);

  // Level 3 of B071, B075 and B074
  gMC->Gsatt("FTOA","seen",-2);  // all FTOA sub-levels skipped   -
  if (fTOFHoles) {
    gMC->Gsatt("FTOB","seen",-2);  // all FTOB sub-levels skipped   -
    gMC->Gsatt("FTOC","seen",-2);  // all FTOC sub-levels skipped   -
  }

  gMC->Gdopt("hide","on");
  gMC->Gdopt("shad","on");
  gMC->Gsatt("*", "fill", 5);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 100, 1000, 100, 1000, 0, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 10, 9.5, .018, .018);
  gMC->Gdhead(1111,"TOF detector");
  gMC->Gdman(18, 3, "MAN");
  gMC->Gdopt("hide","off");
}                                 

//_____________________________________________________________________________
void AliTOFv5T0::DrawDetectorStrips() const
{
  //
  // Draw a shaded view of the TOF strips for version 4
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);

  //
  //Set volumes visible
  // 
  
  //Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN", 0);
  
//=====> Level 1
  // Level 1 for TOF volumes
  gMC->Gsatt("B077","seen", 0);

//=====> Level 2
  // Level 2 for TOF volumes
  gMC->Gsatt("B071","seen", 0);
  gMC->Gsatt("B074","seen", 0);
  gMC->Gsatt("B075","seen", 0);
  gMC->Gsatt("B076","seen",-1); // all B076 sub-levels skipped -
  gMC->Gsatt("B080","seen", 0);  // B080 does not has sub-level                

  // Level 2 of B071
  gMC->Gsatt("B063","seen",-1); // all B063 sub-levels skipped   -
  gMC->Gsatt("B065","seen",-1); // all B065 sub-levels skipped   -
  gMC->Gsatt("B067","seen",-1); // all B067 sub-levels skipped   -
  gMC->Gsatt("B056","seen", 0);  // B056 does not has sub-levels  -
  gMC->Gsatt("B072","seen",-1); // all B072 sub-levels skipped   -

  gMC->Gsatt("BTR1","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO1","seen", 0);

  // Level 2 of B074
  gMC->Gsatt("BTR2","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO2","seen", 0);

  // Level 2 of B075
  gMC->Gsatt("BTR3","seen", 0);  // all BTR1 sub-levels skipped   -
  gMC->Gsatt("BTO3","seen", 0);

  // Level 3 of B071, B074 and B075
  gMC->Gsatt("FTOA","SEEN", 0);
  if (fTOFHoles) {
    gMC->Gsatt("FTOB","SEEN", 0);
    gMC->Gsatt("FTOC","SEEN", 0);
  }

  // Level 4 of B071, B074 and B075
  gMC->Gsatt("FLTA","SEEN", 0);
  if (fTOFHoles) {
    gMC->Gsatt("FLTB","SEEN", 0);
    gMC->Gsatt("FLTC","SEEN", 0);
  }

  // Level 5 of B071, B074 and B075
  gMC->Gsatt("FAIA","SEEN",-1);  // all FAIA sub-levels skipped   -
  if (fTOFHoles) {
    gMC->Gsatt("FAIB","SEEN",-1);  // all FAIB sub-levels skipped   -
    gMC->Gsatt("FAIC","SEEN",-1);  // all FAIC sub-levels skipped   -
  }

  gMC->Gsatt("FALA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FALB","SEEN", 0);

  gMC->Gsatt("FPEA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FPEB","SEEN", 0);

  gMC->Gsatt("FSTR","SEEN",-2);  // all FSTR sub-levels skipped   -

  gMC->Gsatt("FWZ1","SEEN", 0);
  gMC->Gsatt("FWZ2","SEEN", 0);
  gMC->Gsatt("FWZ3","SEEN", 0);
  gMC->Gsatt("FWZ4","SEEN", 0);
  if (fTOFHoles){
    gMC->Gsatt("FZ1B","SEEN", 0);
    gMC->Gsatt("FZ1C","SEEN", 0);
    gMC->Gsatt("FZ2B","SEEN", 0);
  }

  /*
  // Level 2 of FAIA
  // Level 2 of FAIB
  // Level 2 of FAIC
  gMC->Gsatt("FALP","SEEN",0);
  gMC->Gsatt("FCAR","SEEN",0);
  gMC->Gsatt("FTUB","SEEN",-1);  // all FTUB sub-levels skipped   -

  // Level 2 of FTUB
  gMC->Gsatt("FITU","SEEN",0);
  */

  /*
  // Level 2 of FSTR
  gMC->Gsatt("FGLF","SEEN",0);
  gMC->Gsatt("FHON","SEEN",0);
  gMC->Gsatt("FPC1","SEEN",0);
  gMC->Gsatt("FPC2","SEEN",0);
  gMC->Gsatt("FPCB","SEEN",0);
  gMC->Gsatt("FRGL","SEEN",0);

  // Level 2 of FPCB => Level 3 of FSTR
  gMC->Gsatt("FSEN","SEEN",0);
  gMC->Gsatt("FSEZ","SEEN",0);
  gMC->Gsatt("FPAD","SEEN",1);
  */

  gMC->Gdopt("hide","on");
  gMC->Gdopt("shad","on");
  gMC->Gsatt("*", "fill", 5);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, 0, 1000, 0, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 10, 9.5, .018, .018);
  gMC->Gdhead(1111,"TOF Strips");
  gMC->Gdman(18, 3, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTOFv5T0::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //

  //AliTOF::CreateMaterials();

  AliMagF *magneticField = (AliMagF*)gAlice->Field();

  Int_t   isxfld = magneticField->Integ();
  Float_t sxmgmx = magneticField->Max();

  Float_t we[7], ae[7], na[7], fr[7], vl[7];
  Int_t i;

  //--- Quartz (SiO2) to simulate float glass
  //    density tuned to have correct float glass 
  //    radiation length
  Float_t   aq[2] = { 28.0855,15.9994 };
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  Float_t   dq = 2.55; // std value: 2.2
  Int_t nq = -2;

  // --- Nomex
  Float_t anox[4] = {12.01,1.01,16.00,14.01};
  Float_t znox[4] = { 6.,  1.,  8.,  7.};
  Float_t wnox[4] = {14., 22., 2., 2.};
  Float_t dnox  = 0.048;
  Int_t nnox   = -4;

  //  { Si, C, H, O }
  Float_t ag10[4] = {28.09,12.01,1.01,16.00};
  Float_t zg10[4] = {14.,  6.,  1.,  8.};
  Float_t wmatg10[4];
  Int_t nlmatg10  = 4;
  for (i = 0; i < nlmatg10; ++i) {
    ae[i] = ag10[i];
    vl[i] = 1.;
  }
  ae[4] = 16.00;
  vl[4] = 1.;
  na[0] = 1.; 
  na[1] = 14.;
  na[2] = 20.;
  na[3] = 2.;
  na[4] = 3.;
  fr[0] = 0.6;
  fr[1] = 0.4;
  fr[2] = 0.4;
  fr[3] = 0.6;
  fr[4] = 0.4;
  MaterialMixer(we,ae,na,fr,vl,5);
  we[3] += we[4];
  wmatg10[0]= we[0];
  wmatg10[1]= we[1];
  wmatg10[2]= we[2];
  wmatg10[3]= we[3];
  Float_t densg10 = 1.7;

  // -- Water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;
  Int_t nwa = -2;

  // stainless steel
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };

  // AIR
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;

  // --- fibre glass
  Float_t afg[4] = {28.09,16.00,12.01,1.01};
  Float_t zfg[4] = {14., 8., 6., 1.};
  Float_t wfg[4] = {0.12906,0.29405,0.51502,0.06187};
  Float_t dfg  = 1.111;
  Int_t nfg   = 4;

  // --- Freon C2F4H2 + SF6
  Float_t afre[4]= {12.01,1.01,19.00,32.07};
  Float_t zfre[4]= { 6., 1., 9., 16.};
  Float_t wfre[4]= {0.21250,0.01787,0.74827,0.021355};
  Float_t densfre= 0.00375;
  Int_t nfre  = 4;

  //char namat[15] = "            ";
  //Float_t ama[2], zma[2], dma, radl, absl, buf[1];
  //Int_t nbuf;

  AliMixture ( 0, "Air$", aAir, zAir, dAir, 4, wAir);
  AliMixture ( 1, "Nomex$", anox, znox, dnox, nnox, wnox);
  AliMixture ( 2, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  AliMixture ( 3, "fibre glass$", afg, zfg, dfg, nfg, wfg);
  AliMaterial( 4, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial( 5, "Al honeycomb$", 26.98, 13., 0.0496, 483., 2483.);
  AliMixture ( 6, "Freon$",  afre, zfre, densfre, nfre, wfre);
  AliMixture ( 7, "Glass$", aq, zq, dq, nq, wq);
  /*
  // get freon and glass
  gMC->Gfmate((*fIdmate)[6],namat,ama[0],zma[0],dma,radl,absl,buf,nbuf);
  gMC->Gfmate((*fIdmate)[7],namat,ama[1],zma[1],dma,radl,absl,buf,nbuf);

  // --- glass-freon
  Float_t wgfr[2]= {0.0011,0.9989};
  Float_t dgfr = 1.434;
  Int_t ngfr  = 2;
  AliMixture ( 8, "glass-freon$", ama, zma, dgfr, ngfr, wgfr);
  */
  AliMixture ( 9, "Water$",  awa, zwa, dwa, nwa, wwa);
  AliMixture (10, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);

  Float_t epsil, stmin, deemax, stemax;
 
  //   Previous data
  //       EPSIL  = 0.1   ! Tracking precision,
  //       STEMAX = 0.1   ! Maximum displacement for multiple scattering
  //       DEEMAX = 0.1   ! Maximum fractional energy loss, DLS
  //       STMIN  = 0.1

  //   New data
  epsil  = .001;  // Tracking precision,
  stemax = -1.;   // Maximum displacement for multiple scattering
  deemax = -.3;   // Maximum fractional energy loss, DLS
  stmin  = -.8;

  AliMedium( 1, "Air$",         0, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 2,"Nomex$",        1, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 3,"G10$",          2, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 4,"fibre glass$",  3, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  //AliMedium( 5,"glass-freon$",  8, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 6,"Al Frame$",     4, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 7,"Al honeycomb$", 5, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 8,"Fre$",          6, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 9,"PCB-S$",        2, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(10,"Glass$",        7, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(11,"Water$",        9, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(12,"STEEL$",       10, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);

}
//_____________________________________________________________________________
void AliTOFv5T0::Init()
{
  //
  // Initialise the detector after the geometry has been defined
  //
  AliDebug(1, "**************************************"
           "  TOF  "
           "**************************************");
  AliDebug(1, "  Version 4 of TOF initialing, "
	   "symmetric TOF - Full Coverage version");
  
  AliTOF::Init();
  
  fIdFTOA = gMC->VolId("FTOA");
  if (fTOFHoles) {
    fIdFTOB = gMC->VolId("FTOB");
    fIdFTOC = gMC->VolId("FTOC");
  }
  fIdFLTA = gMC->VolId("FLTA");
  if (fTOFHoles) {
    fIdFLTB = gMC->VolId("FLTB");
    fIdFLTC = gMC->VolId("FLTC");
  }

  AliDebug(1, "**************************************"
           "  TOF  "
           "**************************************");
}
 
//_____________________________________________________________________________
void AliTOFv5T0::StepManager()
{

  //
  // Procedure called at each step in the Time Of Flight
  //

  TLorentzVector mom, pos;
  Float_t xm[3],pm[3],xpad[3],ppad[3];
  Float_t hits[14];
  Int_t   vol[5];
  Int_t   sector, plate, padx, padz, strip;
  Int_t   copy, padzid, padxid, stripid, i;
  Int_t   *idtmed = fIdtmed->GetArray()-499;
  Float_t incidenceAngle;

  const char* volpath;

  Int_t index = 0;

  if(
     gMC->IsTrackEntering()
     && gMC->TrackCharge()
     //&& gMC->GetMedium()==idtmed[508]
     && gMC->CurrentMedium()==idtmed[508]
     && gMC->CurrentVolID(copy)==fIdSens
     )
  {

    AliMC *mcApplication = (AliMC*)gAlice->GetMCApp();

    AddTrackReference(mcApplication->GetCurrentTrackNumber());
    //AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());

    // getting information about hit volumes
    
    padzid=gMC->CurrentVolOffID(1,copy);
    padz=copy;
    padz--;

    padxid=gMC->CurrentVolOffID(0,copy);
    padx=copy; 
    padx--;
    
    stripid=gMC->CurrentVolOffID(4,copy);
    strip=copy; 
    strip--;

    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);

    Double_t normMom=1./mom.Rho();

    //  getting the coordinates in pad ref system

    xm[0] = (Float_t)pos.X();
    xm[1] = (Float_t)pos.Y();
    xm[2] = (Float_t)pos.Z();

    pm[0] = (Float_t)mom.X()*normMom;
    pm[1] = (Float_t)mom.Y()*normMom;
    pm[2] = (Float_t)mom.Z()*normMom;
 
    gMC->Gmtod(xm,xpad,1); // from MRS to DRS: coordinates convertion
    gMC->Gmtod(pm,ppad,2); // from MRS to DRS: direction cosinus convertion


    if (TMath::Abs(ppad[1])>1) {
      AliWarning("Abs(ppad) > 1");
      ppad[1]=TMath::Sign((Float_t)1,ppad[1]);
    }
    incidenceAngle = TMath::ACos(ppad[1])*kRaddeg;

    plate = -1;
    if      (strip <  fTOFGeometry->NStripC()) {
      plate = 0;
      //strip = strip;
    }
    else if (strip >= fTOFGeometry->NStripC() && 
	     strip <  fTOFGeometry->NStripC() + fTOFGeometry->NStripB()) {
      plate = 1;
      strip = strip - fTOFGeometry->NStripC();
    }
    else if (strip >= fTOFGeometry->NStripC() + fTOFGeometry->NStripB() &&
	     strip <  fTOFGeometry->NStripC() + fTOFGeometry->NStripB() + fTOFGeometry->NStripA()) {
      plate = 2;
      strip = strip - fTOFGeometry->NStripC() - fTOFGeometry->NStripB();
    }
    else if (strip >= fTOFGeometry->NStripC() + fTOFGeometry->NStripB() + fTOFGeometry->NStripA() &&
	     strip <  fTOFGeometry->NStripC() + fTOFGeometry->NStripB() + fTOFGeometry->NStripA() + fTOFGeometry->NStripB()) {
      plate = 3;
      strip = strip - fTOFGeometry->NStripC() - fTOFGeometry->NStripB() - fTOFGeometry->NStripA();
    }
    else                                {
      plate = 4;
      strip = strip - fTOFGeometry->NStripC() - fTOFGeometry->NStripB() - fTOFGeometry->NStripA() - fTOFGeometry->NStripB();
    }

    volpath=gMC->CurrentVolOffName(7);
    index=atoi(&volpath[4]);
    sector=-1;
    sector=index;

    //Old 6h convention
    // if(index<5){
    //   sector=index+13;
    //	}
    // else{
    //   sector=index-5;
    // } 
 
    for(i=0;i<3;++i) {
      hits[i]   = pos[i];
      hits[i+3] = pm[i];
    }

    hits[6] = mom.Rho();
    hits[7] = pos[3];
    hits[8] = xpad[0];
    hits[9] = xpad[1];
    hits[10]= xpad[2];
    hits[11]= incidenceAngle;
    hits[12]= gMC->Edep();
    hits[13]= gMC->TrackLength();
    
    vol[0]= sector;
    vol[1]= plate;
    vol[2]= strip;
    vol[3]= padx;
    vol[4]= padz;    

    AddT0Hit(mcApplication->GetCurrentTrackNumber(),vol, hits);
    //AddT0Hit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol, hits);
  }
}
//-------------------------------------------------------------------
void AliTOFv5T0::MaterialMixer(Float_t* p,Float_t* a,Float_t* m,Float_t* d,Float_t* s,Int_t n) const
{
  // a[] atomic weights vector  (in)
  //     (atoms present in more compound appear separately)
  // m[] number of corresponding atoms in the mixture  (in)
  // d[] fraction of the compound relative to the corresponding atoms  (in)
  // s[] further possible weights     "        "        "        "     (in)
  Float_t t = 0.;
  for (Int_t i = 0; i < n; ++i) {
    p[i] = a[i]*m[i]*d[i]*s[i];
    t  += p[i];
  }
  for (Int_t i = 0; i < n; ++i) {
    p[i] = p[i]/t;
    //    AliInfo(Form((\n weight[%i] = %f (,i,p[i]));
  }
}
