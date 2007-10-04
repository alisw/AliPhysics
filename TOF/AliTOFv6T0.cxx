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
Revision 1.7  2007/10/03 18:07:26  arcelli
right handling of l2t matrices and alignable entries in case of TOF holes (Annalisa)

Revision 1.6  2007/10/03 10:41:16  arcelli
adding tracking-to-local matrices for new AliTOFcluster

Revision 1.5  2007/07/27 08:14:48  morsch
Write all track references into the same branch.

Revision 1.4  2007/05/29 16:51:05  decaro
Update of the front-end electronics and cooling system description

Revision 1.3.2  2007/05/29  decaro
FEA+cooling zone description: update
     FEA+cooling orientation (side A/ side C) -> correction
Revision 1.3.1  2007/05/24  decaro
Change the FEA+cooling zone description:
     - FCA1/FCA2, air boxes, contain:
                 FFEA volume, G10 box,
	         FAL1/FAL2/FAL3 volumes, aluminium boxes;
     - FRO1/FRO2/FRO3/FRO4/FBAR, aluminum boxes;
     - changed FTUB positions;

Revision 1.3  2007/05/04 14:05:42  decaro
Ineffective comment cleanup

Revision 1.2  2007/05/04 12:59:22  arcelli
Change the TOF SM paths for misalignment (one layer up)

Revision 1.1  2007/05/02 17:32:58  decaro
TOF geometry description as installed (G. Cara Romeo, A. De Caro)

Revision 0.1 2007 March G. Cara Romeo and A. De Caro
        Implemented a more realistic TOF geometry description,
	in terms of:
	   - material badget,
	   - services and front end electronics description,
	   - TOF crate readout modules
	     (added volume FTOS in ALIC_1/BBMO_1/BBCE_%i -for i=1,...,18-,
	      and in ALIC_1/BFMO_%i -for i=19,...,36- volumes)
	As the 5th version in terms of geometrical positioning of volumes.

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the functions for version 6 of the Time Of Flight    //
//  detector.                                                                //
//                                                                           //
//  VERSION WITH 6 MODULES AND TILTED STRIPS                                 //
//                                                                           //
//  FULL COVERAGE VERSION + OPTION for PHOS holes                            //
//                                                                           //
//                                                                           //
//Begin_Html                                                                 //
/*                                                                           //
<img src="picts/AliTOFv6T0Class.gif">                                        //
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
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TGeoVolume.h>

#include "AliConst.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliTrackReference.h"

#include "AliTOFGeometry.h"
#include "AliTOFv6T0.h"

extern TDirectory *gDirectory;
extern TVirtualMC *gMC;
extern TGeoManager *gGeoManager;

extern AliRun *gAlice;

ClassImp(AliTOFv6T0)

//_____________________________________________________________________________
  AliTOFv6T0::AliTOFv6T0():
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
AliTOFv6T0::AliTOFv6T0(const char *name, const char *title):
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
    fTOFGeometry = new AliTOFGeometry();

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
void AliTOFv6T0::AddAlignableVolumes() const
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

      if (fTOFHoles && (isect==11 || isect==12)) {
	if (istr<39) {
	  vpL3 = "/FTOB_0";
	  vpL4 = "/FLTB_0/FSTR_";
	}
	else if (istr>53) {
	  vpL3 = "/FTOC_0";
	  vpL4 = "/FLTC_0/FSTR_";
	}
	else continue;
      }
      else {
	vpL3 = "/FTOA_0";
	vpL4 = "/FLTA_0/FSTR_";
      }

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

      //T2L matrices for alignment
      TGeoPNEntry *e = gGeoManager->GetAlignableEntry(symName.Data());
      if (e) {
	const char *path = e->GetTitle();
	if (!gGeoManager->cd(path)) {
	  AliFatal(Form("Volume path %s not valid!",path));
	}
	TGeoHMatrix *globMatrix = gGeoManager->GetCurrentMatrix();
	Double_t phi = 20.0 * (isect % 18) + 10.0;
	TGeoHMatrix *t2l  = new TGeoHMatrix();
	t2l->RotateZ(phi);
	t2l->MultiplyLeft(&(globMatrix->Inverse()));
	e->SetMatrix(t2l);
      }
      else {
	AliError(Form("Alignable entry %s is not valid!",symName.Data()));
      }

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
void AliTOFv6T0::BuildGeometry()
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
  
  for (Int_t nodeNum=1;nodeNum<kNTof+1;nodeNum++){
    
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
void AliTOFv6T0::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv6T0.gif">
  */
  //End_Html
  //
  // Creates common geometry
  //
  AliTOF::CreateGeometry();
}
 

//_____________________________________________________________________________
void AliTOFv6T0::TOFpc(Float_t xtof, Float_t ytof, Float_t zlenA)
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

  // module wall thickness (cm)
  const Float_t kModuleWallThickness = 0.33;

  // honeycomb layer between strips and cards (cm)
  const Float_t kHoneycombLayerThickness = 2.;

  AliDebug(1, "************************* TOF geometry **************************");
  AliDebug(1,Form(" xtof   %d",  xtof));
  AliDebug(1,Form(" ytof   %d",  ytof));
  AliDebug(1,Form(" zlenA   %d", zlenA));
  AliDebug(2,Form(" zlenA*0.5 = %d", zlenA*0.5));
  
  // Definition of the of fibre glass modules (FTOA, FTOB and FTOC)
    
  Float_t  xcoor, ycoor, zcoor;
  Float_t  par[3];
  Int_t    *idtmed = fIdtmed->GetArray()-499;
  Int_t    idrotm[100];

  par[0] = xtof * 0.5;
  par[1] = ytof * 0.25;
  par[2] = zlenA * 0.5;
  gMC->Gsvolu("FTOA", "BOX ", idtmed[503], par, 3);  // fibre glass
   
  if (fTOFHoles) {
    par[0] =  xtof * 0.5;
    par[1] =  ytof * 0.25;
    par[2] = (zlenA*0.5 - kInterCentrModBorder1)*0.5;
    gMC->Gsvolu("FTOB", "BOX ", idtmed[503], par, 3);  // fibre glass
    gMC->Gsvolu("FTOC", "BOX ", idtmed[503], par, 3);  // fibre glass
  }

  // New supermodule card section description
  //  2 cm  honeycomb layer between strips and cards
  par[0] = xtof*0.5 + 2.;
  par[1] = kHoneycombLayerThickness*0.5;
  par[2] = zlenA*0.5 + 2.;
  gMC->Gsvolu("FPEA", "BOX ", idtmed[506], par, 3);    // Al + Cu honeycomb
  if (fTOFHoles) {
    //par[0] = xtof*0.5 + 2.;
    //par[1] = kHoneycombLayerThickness*0.5;
    par[2] = (zlenA*0.5 - kInterCentrModBorder1)*0.5 + 2.;
    gMC->Gsvolu("FPEB", "BOX ", idtmed[506], par, 3);  // Al + Cu honeycomb
  }

  // Definition of the air card containers (FAIA and FAIB)

  par[0] = xtof*0.5;
  par[1] = (ytof*0.5 - kHoneycombLayerThickness)*0.5;
  par[2] = zlenA*0.5;
  gMC->Gsvolu("FAIA", "BOX ", idtmed[500], par, 3);                // Air
  if (fTOFHoles) gMC->Gsvolu("FAIB", "BOX ", idtmed[500], par, 3); // Air

  // Positioning of fibre glass modules (FTOA, FTOB and FTOC) and
  // card containers (FPEA, FAIA and FAIB)

  //AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,-90.);
  AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,270.);
  xcoor = 0.;
  for(Int_t isec=0; isec<fTOFGeometry->NSectors(); isec++){
    if(fTOFSectors[isec]==-1)continue;
    char name[16];
    sprintf(name, "BTOF%d",isec);
    if (fTOFHoles && (isec==11||isec==12)) {
    //if (fTOFHoles && (isec==16||isec==17)) { \\Old 6h convention
      //xcoor = 0.;
      ycoor = (zlenA*0.5 + kInterCentrModBorder1)*0.5;
      zcoor = -ytof * 0.25;
      gMC->Gspos("FTOB", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
      gMC->Gspos("FTOC", 0, name, xcoor,-ycoor, zcoor, idrotm[0], "ONLY");
      //xcoor = 0.;
      //ycoor = (zlenA*0.5 + kInterCentrModBorder1)*0.5;
      zcoor = kHoneycombLayerThickness*0.5;
      gMC->Gspos("FPEB", 1, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
      gMC->Gspos("FPEB", 2, name, xcoor,-ycoor, zcoor, idrotm[0], "ONLY");
      //xcoor = 0.;
      ycoor = 0.;
      zcoor = kHoneycombLayerThickness + (ytof*0.5 - kHoneycombLayerThickness)*0.5;
      gMC->Gspos("FAIB", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
    }
    else {
      //xcoor = 0.;
      ycoor = 0.;
      zcoor = -ytof * 0.25;
      gMC->Gspos("FTOA", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
      //xcoor = 0.;
      //ycoor = 0.;
      zcoor = kHoneycombLayerThickness*0.5;
      gMC->Gspos("FPEA", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
      //xcoor = 0.;
      //ycoor = 0.;
      zcoor = kHoneycombLayerThickness + (ytof*0.5 - kHoneycombLayerThickness)*0.5;
      gMC->Gspos("FAIA", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
    }
  }

  // Definition and positioning
  // of the not sensitive volumes with Insensitive Freon (FLTA, FLTB and FLTC)

  Float_t xFLT, yFLT, zFLTA;
  
  xFLT  = xtof     - kModuleWallThickness*2.;
  yFLT  = ytof*0.5 - kModuleWallThickness;
  zFLTA = zlenA    - kModuleWallThickness*2.;
  
  par[0] = xFLT*0.5;
  par[1] = yFLT*0.5;
  par[2] = zFLTA*0.5;
  gMC->Gsvolu("FLTA", "BOX ", idtmed[507], par, 3); //  Freon mix

  xcoor = 0.;
  ycoor = kModuleWallThickness*0.5;
  zcoor = 0.;
  gMC->Gspos ("FLTA", 0, "FTOA", xcoor, ycoor, zcoor, 0, "ONLY");

  if (fTOFHoles) {
    par[2] = (zlenA*0.5 - kInterCentrModBorder1 - kModuleWallThickness)*0.5;
    gMC->Gsvolu("FLTB", "BOX ", idtmed[507], par, 3); // Freon mix
    gMC->Gsvolu("FLTC", "BOX ", idtmed[507], par, 3); // Freon mix

    //xcoor = 0.;
    //ycoor = kModuleWallThickness*0.5;
    //zcoor = 0.;
    gMC->Gspos ("FLTB", 0, "FTOB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos ("FLTC", 0, "FTOC", xcoor, ycoor, zcoor, 0, "ONLY");
  }

  Float_t alpha, tgal, beta, tgbe, trpa[11];

  // Definition and positioning
  // of the fibre glass walls between central and intermediate modules (FWZ1 and FWZ2)

  tgal = (yFLT - 2.*kLengthInCeModBorder)/(kInterCentrModBorder2 - kInterCentrModBorder1);
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
  ycoor = -(yFLT - kLengthInCeModBorder)*0.5;
  zcoor = kInterCentrModBorder1;
  gMC->Gspos("FWZ1", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[1],"ONLY");
  gMC->Gspos("FWZ1", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[4],"ONLY");

  AliMatrix (idrotm[2],90.,270.,  0.,0.,90.,180.);
  AliMatrix (idrotm[5],90.,270.,180.,0.,90.,  0.);

  xcoor = 0.;
  ycoor = (yFLT - kLengthInCeModBorder)*0.5;
  zcoor = kInterCentrModBorder2;
  gMC->Gspos("FWZ1", 3,"FLTA", xcoor, ycoor, zcoor,idrotm[2],"ONLY");
  gMC->Gspos("FWZ1", 4,"FLTA", xcoor, ycoor,-zcoor,idrotm[5],"ONLY");

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
  ycoor = 0.;
  zcoor = (kInterCentrModBorder2 + kInterCentrModBorder1)*0.5;
  gMC->Gspos("FWZ2", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[3],"ONLY");
  gMC->Gspos("FWZ2", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[6],"ONLY");

  // Definition and positioning
  // of the fibre glass walls between intermediate and lateral modules (FWZ3 and FWZ4)

  tgal = (yFLT - 2.*kLengthExInModBorder)/(kExterInterModBorder2 - kExterInterModBorder1);
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
  ycoor = (yFLT - kLengthExInModBorder)*0.5;
  zcoor = kExterInterModBorder1;
  gMC->Gspos("FWZ3", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[5],"ONLY");
  gMC->Gspos("FWZ3", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[2],"ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = (yFLT - kLengthExInModBorder)*0.5;
    zcoor = -kExterInterModBorder1 + (zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5;
    gMC->Gspos("FWZ3", 5,"FLTB", xcoor, ycoor, zcoor,idrotm[2],"ONLY");
    gMC->Gspos("FWZ3", 6,"FLTC", xcoor, ycoor,-zcoor,idrotm[5],"ONLY");
  }

  //xcoor = 0.;
  ycoor = -(yFLT - kLengthExInModBorder)*0.5;
  zcoor = kExterInterModBorder2;
  gMC->Gspos("FWZ3", 3,"FLTA", xcoor, ycoor, zcoor,idrotm[4],"ONLY");
  gMC->Gspos("FWZ3", 4,"FLTA", xcoor, ycoor,-zcoor,idrotm[1],"ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = -(yFLT - kLengthExInModBorder)*0.5;
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

  //xcoor = 0.;
  ycoor = 0.;
  zcoor = (kExterInterModBorder2 + kExterInterModBorder1)*0.5;
  gMC->Gspos("FWZ4", 1,"FLTA", xcoor, ycoor, zcoor,idrotm[16],"ONLY");
  gMC->Gspos("FWZ4", 2,"FLTA", xcoor, ycoor,-zcoor,idrotm[13],"ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = 0.;
    zcoor = -(kExterInterModBorder2 + kExterInterModBorder1)*0.5 +
      (zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5;
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
  const Float_t khhony   = 1.0;       // heigth of HONY  Layer
  const Float_t khpcby   = 0.08;      // heigth of PCB   Layer
  const Float_t khrgly   = 0.055;     // heigth of RED GLASS  Layer

  const Float_t khfiliy  = 0.125;     // heigth of FISHLINE  Layer
  const Float_t khglassy = 0.160*0.5; // heigth of GLASS  Layer
  const Float_t khglfy   = khfiliy+2.*khglassy; // heigth of GLASS+FISHLINE  Layer

  const Float_t khcpcby  = 0.16;      // heigth of PCB  Central Layer
  const Float_t kwhonz   = 8.1;       // z dimension of HONEY  Layer
  const Float_t kwpcbz1  = 10.6;      // z dimension of PCB  Lower Layer
  const Float_t kwpcbz2  = 11.6;      // z dimension of PCB  Upper Layer
  const Float_t kwcpcbz  = 13.;       // z dimension of PCB  Central Layer
  const Float_t kwrglz   = 8.;        // z dimension of RED GLASS  Layer
  const Float_t kwglfz   = 7.;        // z dimension of GLASS+FISHLN Layer
  const Float_t klsensmx = knx*kPadX; // length of Sensitive Layer
  const Float_t khsensmy = 0.05;      // heigth of Sensitive Layer
  const Float_t kwsensmz = knz*kPadZ; // width of Sensitive Layer

  // heigth of the FSTR Volume (the strip volume)
  const Float_t khstripy = 2.*khhony+2.*khpcby+4.*khrgly+2.*khglfy+khcpcby;

  // width  of the FSTR Volume (the strip volume)
  const Float_t kwstripz = kwcpcbz;
  // length of the FSTR Volume (the strip volume)
  const Float_t klstripx = fTOFGeometry->StripLength();
  
  Float_t parfp[3]={klstripx*0.5, khstripy*0.5, kwstripz*0.5};
  // Coordinates of the strip center in the strip reference frame;
  // used for positioning internal strip volumes
  Float_t posfp[3]={0.,0.,0.};

  // FSTR volume definition-filling this volume with non sensitive Gas Mixture
  gMC->Gsvolu("FSTR","BOX",idtmed[507],parfp,3); // Freon mix

  //-- HONY Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khhony*0.5;
  parfp[2] = kwhonz*0.5;
  gMC->Gsvolu("FHON","BOX",idtmed[501],parfp,3); // honeycomb (Nomex)
  // positioning 2 HONY Layers on FSTR volume
  //posfp[0] = 0.;
  posfp[1] =-khstripy*0.5+parfp[1];
  //posfp[2] = 0.;
  gMC->Gspos("FHON",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FHON",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  
  //-- PCB Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khpcby*0.5;
  parfp[2] = kwpcbz1*0.5;
  gMC->Gsvolu("FPC1","BOX",idtmed[502],parfp,3); // G10
  //parfp[0] = klstripx*0.5;
  //parfp[1] = khpcby*0.5;
  parfp[2] = kwpcbz2*0.5;
  gMC->Gsvolu("FPC2","BOX",idtmed[502],parfp,3); // G10
  // positioning 2 PCB Layers on FSTR volume
  //posfp[0] = 0.;
  posfp[1] =-khstripy*0.5+khhony+parfp[1];
  //posfp[2] = 0.;
  gMC->Gspos("FPC1",1,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FPC2",1,"FSTR",0., posfp[1],0.,0,"ONLY");

  //-- central PCB layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khcpcby*0.5;
  parfp[2] = kwcpcbz*0.5;
  gMC->Gsvolu("FPCB","BOX",idtmed[502],parfp,3); // G10
  // positioning the central PCB layer
  gMC->Gspos("FPCB",1,"FSTR",0.,0.,0.,0,"ONLY");

  //      Sensitive volume
  Float_t parfs[3] = {klsensmx*0.5, khsensmy*0.5, kwsensmz*0.5};
  gMC->Gsvolu("FSEN","BOX",idtmed[508],parfs,3); // sensitive
  // dividing FSEN along z in knz=2 and along x in knx=48
  gMC->Gsdvn("FSEZ","FSEN",knz,3);
  gMC->Gsdvn("FPAD","FSEZ",knx,1);
  // positioning a Sensitive layer inside FPCB
  gMC->Gspos("FSEN",1,"FPCB",0.,0.,0.,0,"ONLY");

  //-- RED GLASS Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khrgly*0.5;
  parfp[2] = kwrglz*0.5;
  gMC->Gsvolu("FRGL","BOX",idtmed[509],parfp,3); // glass
  // positioning 4 RED GLASS Layers on FSTR volume
  //posfp[0] = 0.;
  posfp[1] = -khstripy*0.5+khhony+khpcby+parfp[1];
  //posfp[2] = 0.;
  gMC->Gspos("FRGL",1,"FSTR",0., posfp[1],0.,0,"ONLY");
  gMC->Gspos("FRGL",4,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  //posfp[0] = 0.;
  posfp[1] = (khcpcby+khrgly)*0.5;
  //posfp[2] = 0.;
  gMC->Gspos("FRGL",2,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FRGL",3,"FSTR",0., posfp[1],0.,0,"ONLY");

  //-- GLASS+FISHLINE Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khglfy*0.5;
  parfp[2] = kwglfz*0.5;
  gMC->Gsvolu("FGLF","BOX",idtmed[504],parfp,3);

  // positioning 2 GLASS+FISHLINE Layers on FSTR volume
  //posfp[0] = 0.;
  posfp[1] = (khcpcby + khglfy)*0.5 + khrgly;
  //posfp[2] = 0.;
  gMC->Gspos("FGLF",1,"FSTR",0.,-posfp[1],0.,0,"ONLY");
  gMC->Gspos("FGLF",2,"FSTR",0., posfp[1],0.,0,"ONLY");

  //  Positioning the Strips (FSTR volumes) in the FLT volumes
  Int_t maxStripNumbers [5] ={fTOFGeometry->NStripC(),
			      fTOFGeometry->NStripB(),
			      fTOFGeometry->NStripA(),
			      fTOFGeometry->NStripB(),
			      fTOFGeometry->NStripC()};

  Int_t totalStrip = 0;
  Float_t xpos, zpos, ypos, ang;
  for(Int_t iplate = 0; iplate < fTOFGeometry->NPlates(); iplate++){
    if (iplate>0) totalStrip += maxStripNumbers[iplate-1];
    for(Int_t istrip = 0; istrip < maxStripNumbers[iplate]; istrip++){

      ang = fTOFGeometry->GetAngles(iplate,istrip);
      AliDebug(1, Form(" iplate = %1i, istrip = %2i ---> ang = %f", iplate, istrip, ang));
 
      if (ang>0.)       AliMatrix (idrotm[istrip+totalStrip+1],90.,0.,90.+ang,90., ang, 90.);
      else if (ang==0.) AliMatrix (idrotm[istrip+totalStrip+1],90.,0.,90.,90., 0., 0.);
      else if (ang<0.)  AliMatrix (idrotm[istrip+totalStrip+1],90.,0.,90.+ang,90.,-ang,270.);

      xpos = 0.;
      zpos = fTOFGeometry->GetDistances(iplate,istrip);
      ypos = fTOFGeometry->GetHeights(iplate,istrip) + yFLT*0.5;

      gMC->Gspos("FSTR",istrip+totalStrip+1,"FLTA", xpos, ypos,-zpos,idrotm[istrip+totalStrip+1],  "ONLY");

      if (fTOFHoles) {
	if (istrip+totalStrip+1>53)
	  gMC->Gspos("FSTR",istrip+totalStrip+1,"FLTC", xpos, ypos,-zpos-(zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5,idrotm[istrip+totalStrip+1],"ONLY");
	if (istrip+totalStrip+1<39)
	  gMC->Gspos("FSTR",istrip+totalStrip+1,"FLTB", xpos, ypos,-zpos+(zlenA*0.5 + kInterCentrModBorder1 - kModuleWallThickness)*0.5,idrotm[istrip+totalStrip+1],"ONLY");
      }
    }
  }

  // Definition of the cards, cooling tubes and layer for thermal dispersion
  // (3 volumes)

  // card volume definition
  //Float_t carpar[3] = {9.5, 5.75, 0.5};
  Float_t carpar[3] = {9.5, 5.6, 0.55};
  //gMC->Gsvolu("FCA1", "BOX ", idtmed[514], carpar, 3);   // PCB+Alu small Card 
  gMC->Gsvolu("FCA1", "BOX ", idtmed[500], carpar, 3);   // air
  carpar[0] = 19.25;
  //carpar[1] =  5.6;//5.75;
  //carpar[2] =  0.55;//0.5;
  //gMC->Gsvolu("FCA2", "BOX ", idtmed[514], carpar, 3);   // PCB+Alu long Card 
  gMC->Gsvolu("FCA2", "BOX ", idtmed[500], carpar, 3);   // air


  Float_t feaParam1[3] = {9.5, 5.6, 0.1};
  gMC->Gsvolu("FFEA", "BOX ", idtmed[502], feaParam1, 3);   // G10

  Float_t al1[3] = {9.5, 0.5, 0.25};
  gMC->Gsvolu("FAL1", "BOX ", idtmed[505], al1, 3);   // Aluminium
  Float_t al2[3] = {7.2, 0.8, 0.25};
  gMC->Gsvolu("FAL2", "BOX ", idtmed[505], al2, 3);   // Aluminium
  Float_t al3[3] = {3.35, 3.7, 0.1};
  gMC->Gsvolu("FAL3", "BOX ", idtmed[505], al3, 3);   // Aluminium

  gMC->Gspos("FFEA", 1, "FCA1", 0., 0., -carpar[2]+feaParam1[2], 0, "ONLY");
  gMC->Gspos("FAL1", 1, "FCA1", 0., carpar[1]-al1[1], -carpar[2]+2.*feaParam1[2]+al1[2], 0, "ONLY");
  gMC->Gspos("FAL3", 1, "FCA1", 0., carpar[1]-al3[1],  carpar[2]-al3[2], 0, "ONLY");
  gMC->Gspos("FAL2", 1, "FCA1", 0., carpar[1]-2.*al3[1],  carpar[2]-2.*al3[2]-al2[2], 0, "ONLY");


  gMC->Gspos("FFEA", 2, "FCA2", -(feaParam1[0]+0.25), 0., -carpar[2]+feaParam1[2], 0, "ONLY");
  gMC->Gspos("FAL1", 2, "FCA2", -(feaParam1[0]+0.25), carpar[1]-al1[1], -carpar[2]+2.*feaParam1[2]+al1[2], 0, "ONLY");
  gMC->Gspos("FAL3", 2, "FCA2", -(feaParam1[0]+0.25), carpar[1]-al3[1],  carpar[2]-al3[2], 0, "ONLY");
  gMC->Gspos("FAL2", 2, "FCA2", -(feaParam1[0]+0.25), carpar[1]-2.*al3[1],  carpar[2]-2.*al3[2]-al2[2], 0, "ONLY");

  gMC->Gspos("FFEA", 3, "FCA2",  (feaParam1[0]+0.25), 0., -carpar[2]+feaParam1[2], 0, "ONLY");
  gMC->Gspos("FAL1", 3, "FCA2",  (feaParam1[0]+0.25), carpar[1]-al1[1], -carpar[2]+2.*feaParam1[2]+al1[2], 0, "ONLY");
  gMC->Gspos("FAL3", 3, "FCA2",  (feaParam1[0]+0.25), carpar[1]-al3[1],  carpar[2]-al3[2], 0, "ONLY");
  gMC->Gspos("FAL2", 3, "FCA2",  (feaParam1[0]+0.25), carpar[1]-2.*al3[1],  carpar[2]-2.*al3[2]-al2[2], 0, "ONLY");

  Float_t feaRoof1[3] = {9.5, 0.25, 1.7};
  gMC->Gsvolu("FRO1", "BOX ", idtmed[505], feaRoof1, 3);   // Aluminium
  Float_t feaRoof2[3] = {3.35, 0.05, 1.5};
  gMC->Gsvolu("FRO2", "BOX ", idtmed[505], feaRoof2, 3);   // Aluminium
  Float_t feaRoof3[3] = {3.35, feaRoof1[1]+feaRoof2[1], 0.1};
  gMC->Gsvolu("FRO3", "BOX ", idtmed[505], feaRoof3, 3);   // Aluminium

  Float_t feaRoof4[3] = {3.35,
			 0.05,
			 carpar[2]-feaParam1[2]-al1[2]-al3[2]};
  gMC->Gsvolu("FRO4", "BOX ", idtmed[505], feaRoof4, 3);   // Aluminium

  Float_t bar[3] = {8.575, 0.6, 0.15};
  gMC->Gsvolu("FBAR", "BOX ", idtmed[505], bar, 3);   // Aluminium


  // tube volume definition
  Float_t tubepar[3] = {0., 0.4, xFLT*0.5-15.};
  gMC->Gsvolu("FTUB", "TUBE", idtmed[513], tubepar, 3);  // copper cooling tubes
  //tubepar[0]= 0.;
  tubepar[1]= 0.3;
  //tubepar[2]= xFLT*0.5 - 15.;
  gMC->Gsvolu("FITU", "TUBE", idtmed[510], tubepar, 3);  // cooling water
  // Positioning of the water tube into the steel one
  gMC->Gspos("FITU",1,"FTUB",0.,0.,0.,0,"ONLY");

  // cable
  Float_t cbpar[3] = {0., 0.5, tubepar[2]};
  gMC->Gsvolu("FCAB", "TUBE", idtmed[511], cbpar, 3);    // copper+alu

  // Alluminium components
  Float_t lonpar[3] = {tubepar[2], 6.15, 0.7};
  gMC->Gsvolu("FTLN", "BOX ", idtmed[505], lonpar, 3);   // alluminium
  lonpar[0] = 2.;
  lonpar[1] = 1.;
  lonpar[2] = zlenA*0.5;
  gMC->Gsvolu("FLON", "BOX ", idtmed[505], lonpar, 3);   // alluminium

  // rotation matrix
  AliMatrix(idrotm[99], 180., 90., 90., 90., 90., 0.);
  AliMatrix(idrotm[98],  90.,180., 90., 90.,180., 0.);

  // cards, tubes, cables  positioning
  Float_t carpos[3], rowstep = 6.66, ytub= 3.65, ycab= ytub-3.;
  Float_t rowgap[5] = {13.5, 22.9, 16.94, 23.8, 20.4};
  Int_t row, rowb[5] = {6, 7, 6, 19, 7}, nrow;
  carpos[0] = 25. - xtof*0.5;
  carpos[1] = (11.5 - (ytof*0.5 - kHoneycombLayerThickness))*0.5;
  row = 1;
  for (Int_t sg= -1; sg< 2; sg+= 2) {
    carpos[2] = sg*zlenA*0.5;
    for (Int_t nb=0; nb<5; ++nb) {
      carpos[2] = carpos[2] - sg*(rowgap[nb] - rowstep);
      nrow = row + rowb[nb];
      for ( ; row < nrow ; ++row) {

        carpos[2] -= sg*rowstep;

	if (nb==4) {
	  gMC->Gspos("FCA1",2*row,  "FAIA", carpos[0],carpos[1],carpos[2], 0,"ONLY");
	  gMC->Gspos("FCA1",2*row-1,"FAIA",-carpos[0],carpos[1],carpos[2], 0,"ONLY");
	  gMC->Gspos("FCA2", row,   "FAIA", 0., carpos[1], carpos[2], 0, "ONLY");

	  //gMC->Gspos("FTUB", row, "FAIA", 0., ytub, carpos[2]-sg, idrotm[99], "ONLY");
	  gMC->Gspos("FTUB", row, "FAIA", 0., carpos[1]+carpar[1]-bar[1], carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-bar[1]), idrotm[99], "ONLY");
	  gMC->Gspos("FCAB", row, "FAIA", 0., ycab, carpos[2]-1.1, idrotm[99], "ONLY");

	  gMC->Gspos("FRO1",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");

	  gMC->Gspos("FRO2",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");

	  gMC->Gspos("FRO3",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+(carpar[2]-feaRoof3[2]), 0,"ONLY");

	  gMC->Gspos("FRO4",4*row,  "FAIA", carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-3,"FAIA",-carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");

	  gMC->Gspos("FBAR",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]-bar[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]-bar[1],carpos[2]-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");

	}
	else {
	  switch (sg) {
	  case 1:
	    gMC->Gspos("FCA1",2*row,  "FAIA", carpos[0],carpos[1],carpos[2], 0,"ONLY");
	    gMC->Gspos("FCA1",2*row-1,"FAIA",-carpos[0],carpos[1],carpos[2], 0,"ONLY");
	    gMC->Gspos("FCA2", row,   "FAIA", 0., carpos[1], carpos[2], 0, "ONLY");
	    break;
	  case -1:
	    gMC->Gspos("FCA1",2*row,  "FAIA", carpos[0],carpos[1],carpos[2], idrotm[98],"ONLY");
	    gMC->Gspos("FCA1",2*row-1,"FAIA",-carpos[0],carpos[1],carpos[2], idrotm[98],"ONLY");
	    gMC->Gspos("FCA2", row,   "FAIA", 0., carpos[1], carpos[2], idrotm[98], "ONLY");
	    break;
	  }

	  //gMC->Gspos("FTUB", row, "FAIA", 0., ytub, carpos[2]-sg, idrotm[99], "ONLY");
	  gMC->Gspos("FTUB", row, "FAIA", 0., carpos[1]+carpar[1]-bar[1], carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-bar[1]), idrotm[99], "ONLY");
	  gMC->Gspos("FCAB", row, "FAIA", 0., ycab, carpos[2]-sg*1.1, idrotm[99], "ONLY");

	  gMC->Gspos("FRO1",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");

	  gMC->Gspos("FRO2",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");

	  gMC->Gspos("FRO3",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");

	  gMC->Gspos("FRO4",4*row,  "FAIA", carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-3,"FAIA",-carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");

	  gMC->Gspos("FBAR",4*row,  "FAIA", carpos[0],carpos[1]+carpar[1]-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-1,"FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-2,"FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-3,"FAIA",-carpos[0],carpos[1]+carpar[1]-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");

	}
      }
    }
    gMC->Gspos("FTLN", 5+4*sg, "FAIA", 0., -0.1, 369.9*sg, 0, "ONLY");
    gMC->Gspos("FTLN", 5+3*sg, "FAIA", 0., -0.1, 366.9*sg, 0, "ONLY");
    gMC->Gspos("FTLN", 5+2*sg, "FAIA", 0., -0.1, 198.8*sg, 0, "ONLY");
    gMC->Gspos("FTLN",   5+sg, "FAIA", 0., -0.1, 56.82*sg, 0, "ONLY");
  }
  gMC->Gspos("FCA1", 182, "FAIA", carpos[0],carpos[1],0., 0,"ONLY");
  gMC->Gspos("FCA1", 181, "FAIA",-carpos[0],carpos[1],0., 0,"ONLY");
  gMC->Gspos("FCA2",  91, "FAIA",  0., carpos[1], 0., 0, "ONLY");

  //gMC->Gspos("FTUB",  91, "FAIA",  0., ytub, -1., idrotm[99], "ONLY");
  gMC->Gspos("FTUB", 91, "FAIA", 0., carpos[1]+carpar[1]-bar[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-bar[1]), idrotm[99], "ONLY");
  gMC->Gspos("FCAB", 91, "FAIA",  0., ycab, -1.1, idrotm[99], "ONLY");

  gMC->Gspos("FRO1",364, "FAIA", carpos[0],carpos[1]+carpar[1]+feaRoof1[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
  gMC->Gspos("FRO1",363, "FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
  gMC->Gspos("FRO1",362, "FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
  gMC->Gspos("FRO1",361, "FAIA",-carpos[0],carpos[1]+carpar[1]+feaRoof1[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");

  gMC->Gspos("FRO2",364, "FAIA", carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
  gMC->Gspos("FRO2",363, "FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
  gMC->Gspos("FRO2",362, "FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
  gMC->Gspos("FRO2",361, "FAIA",-carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");

  gMC->Gspos("FRO3",364, "FAIA", carpos[0],carpos[1]+carpar[1]+feaRoof3[1],(carpar[2]-feaRoof3[2]), 0,"ONLY");
  gMC->Gspos("FRO3",363, "FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],(carpar[2]-feaRoof3[2]), 0,"ONLY");
  gMC->Gspos("FRO3",362, "FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],(carpar[2]-feaRoof3[2]), 0,"ONLY");
  gMC->Gspos("FRO3",361, "FAIA",-carpos[0],carpos[1]+carpar[1]+feaRoof3[1],(carpar[2]-feaRoof3[2]), 0,"ONLY");

  gMC->Gspos("FRO4",364, "FAIA", carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
  gMC->Gspos("FRO4",363, "FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
  gMC->Gspos("FRO4",362, "FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
  gMC->Gspos("FRO4",361, "FAIA",-carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");

  gMC->Gspos("FBAR",364, "FAIA", carpos[0],carpos[1]+carpar[1]-bar[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
  gMC->Gspos("FBAR",363, "FAIA", (feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
  gMC->Gspos("FBAR",362, "FAIA",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
  gMC->Gspos("FBAR",361, "FAIA",-carpos[0],carpos[1]+carpar[1]-bar[1],-(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");

  gMC->Gspos("FLON",  2, "FAIA",-24., ytub+1.4, 0., 0, "MANY");
  gMC->Gspos("FLON",  1, "FAIA", 24., ytub+1.4, 0., 0, "MANY");


  if (fTOFHoles) {
    row = 1;
    for (Int_t sg= -1; sg< 2; sg+= 2) {
      carpos[2] = sg*zlenA*0.5;
      for (Int_t nb=0; nb<4; ++nb) {
        carpos[2] = carpos[2] - sg*(rowgap[nb] - rowstep);
        nrow = row + rowb[nb];
        for ( ; row < nrow ; ++row) {
          carpos[2] -= sg*rowstep;

	  switch (sg) {
	  case 1:
	    gMC->Gspos("FCA1",2*row,  "FAIB", carpos[0],carpos[1],carpos[2], 0,"ONLY");
	    gMC->Gspos("FCA1",2*row-1,"FAIB",-carpos[0],carpos[1],carpos[2], 0,"ONLY");
	    gMC->Gspos("FCA2", row,   "FAIB", 0., carpos[1], carpos[2], 0, "ONLY");
	    break;
	  case -1:
	    gMC->Gspos("FCA1",2*row,  "FAIB", carpos[0],carpos[1],carpos[2], idrotm[98],"ONLY");
	    gMC->Gspos("FCA1",2*row-1,"FAIB",-carpos[0],carpos[1],carpos[2], idrotm[98],"ONLY");
	    gMC->Gspos("FCA2", row,   "FAIB", 0., carpos[1], carpos[2], idrotm[98], "ONLY");
	    break;
	  }

          //gMC->Gspos("FTUB", row, "FAIB", 0., ytub,carpos[2]-sg, idrotm[99], "ONLY");
	  gMC->Gspos("FTUB", row, "FAIB", 0., carpos[1]+carpar[1]-bar[1], carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-bar[1]), idrotm[99], "ONLY");
          gMC->Gspos("FCAB", row, "FAIB", 0., ycab,carpos[2]-sg*1.1, idrotm[99], "ONLY");

	  gMC->Gspos("FRO1",4*row,  "FAIB", carpos[0],carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-1,"FAIB", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-2,"FAIB",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");
	  gMC->Gspos("FRO1",4*row-3,"FAIB",-carpos[0],carpos[1]+carpar[1]+feaRoof1[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+feaRoof1[2]), 0,"ONLY");

	  gMC->Gspos("FRO2",4*row,  "FAIB", carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-1,"FAIB", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-2,"FAIB",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");
	  gMC->Gspos("FRO2",4*row-3,"FAIB",-carpos[0],carpos[1]+carpar[1]+2.*feaRoof1[1]+feaRoof2[1],carpos[2]+sg*(carpar[2]-2.*feaRoof3[2]-feaRoof2[2]), 0,"ONLY");

	  gMC->Gspos("FRO3",4*row,  "FAIB", carpos[0],carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-1,"FAIB", (feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-2,"FAIB",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");
	  gMC->Gspos("FRO3",4*row-3,"FAIB",-carpos[0],carpos[1]+carpar[1]+feaRoof3[1],carpos[2]+sg*(carpar[2]-feaRoof3[2]), 0,"ONLY");

	  gMC->Gspos("FRO4",4*row,  "FAIB", carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-1,"FAIB", (feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-2,"FAIB",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");
	  gMC->Gspos("FRO4",4*row-3,"FAIB",-carpos[0],          carpos[1]+carpar[1]+2.*feaRoof1[1]-feaRoof4[1],carpos[2]+sg*(carpar[2]-2.*al3[2]-feaRoof4[2]), 0,"ONLY");

	  gMC->Gspos("FBAR",4*row,  "FAIB", carpos[0],carpos[1]+carpar[1]-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-1,"FAIB", (feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-2,"FAIB",-(feaParam1[0]+0.25),carpos[1]+carpar[1]+-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");
	  gMC->Gspos("FBAR",4*row-3,"FAIB",-carpos[0],carpos[1]+carpar[1]-bar[1],carpos[2]-sg*(carpar[2]-2.*feaParam1[2]-2.*al1[2]+2.*feaRoof1[2]-2.*bar[1]), 0,"ONLY");

        }
      }
      gMC->Gspos("FTLN", 5+4*sg, "FAIB", 0., -0.1, 369.9*sg, 0, "ONLY");
      gMC->Gspos("FTLN", 5+3*sg, "FAIB", 0., -0.1, 366.9*sg, 0, "ONLY");
      gMC->Gspos("FTLN", 5+2*sg, "FAIB", 0., -0.1, 198.8*sg, 0, "ONLY");
      gMC->Gspos("FTLN",   5+sg, "FAIB", 0., -0.1, 56.82*sg, 0, "ONLY");
    }
  gMC->Gspos("FLON", 2, "FAIB",-24., ytub+1.4, 0., 0, "MANY");
  gMC->Gspos("FLON", 1, "FAIB", 24., ytub+1.4, 0., 0, "MANY");
  }

  // Cables and tubes on the side blocks
  const Float_t kcbll   = zlenA*0.5; // length of block
  const Float_t kcbllh  = zlenA*0.5 - kInterCentrModBorder2; // length  of block in case of hole
  const Float_t kcblw   = 13.5;      // width of block
  const Float_t kcblh1  = 2.;        // min. heigth of block
  const Float_t kcblh2  = 12.3;      // max. heigth of block
  // volume definition
  Float_t cblpar[11];
  tgal =  (kcblh2 - kcblh1)/(2.*kcbll);
  cblpar[0] = kcblw *0.5;
  cblpar[1] = 0.;
  cblpar[2] = 0.;
  cblpar[3] = kcbll *0.5;
  cblpar[4] = kcblh1 *0.5;
  cblpar[5] = kcblh2 *0.5;
  cblpar[6] = TMath::ATan(tgal)*kRaddeg;
  cblpar[7] = kcbll *0.5;
  cblpar[8] = kcblh1 *0.5;
  cblpar[9] = kcblh2 *0.5;
  cblpar[10]= cblpar[6];
  gMC->Gsvolu("FCBL", "TRAP", idtmed[512], cblpar, 11); // cables & tubes mix 
  Float_t sawpar[3] = {0.5, kcblh2*0.5, kcbll};
  gMC->Gsvolu("FSAW", "BOX ", idtmed[505], sawpar,  3); // Side Al walls
  // volume positioning
  AliMatrix(idrotm[7], 90., 90., 180., 0., 90., 180.);
  AliMatrix(idrotm[8], 90., 90., 0., 0., 90., 0.);
  xcoor = (xtof-kcblw)*0.5 - 2.*sawpar[0];
  ycoor = (kcblh1+kcblh2)*0.25 - (ytof*0.5 - kHoneycombLayerThickness)*0.5;
  zcoor = kcbll*0.5;
  gMC->Gspos("FCBL", 1, "FAIA", -xcoor, ycoor, -zcoor, idrotm[7], "ONLY");
  gMC->Gspos("FCBL", 2, "FAIA",  xcoor, ycoor, -zcoor, idrotm[7], "ONLY");
  gMC->Gspos("FCBL", 3, "FAIA", -xcoor, ycoor,  zcoor, idrotm[8], "ONLY");
  gMC->Gspos("FCBL", 4, "FAIA",  xcoor, ycoor,  zcoor, idrotm[8], "ONLY");
  xcoor = xtof*0.5-sawpar[0];
  ycoor = (kcblh2 - ytof*0.5 + kHoneycombLayerThickness)*0.5;
  gMC->Gspos("FSAW", 1, "FAIA", -xcoor, ycoor, 0., 0, "ONLY");
  gMC->Gspos("FSAW", 2, "FAIA",  xcoor, ycoor, 0., 0, "ONLY");
  if (fTOFHoles) {
    cblpar[3] = kcbllh *0.5;
    cblpar[5] = kcblh1*0.5 + kcbllh*tgal;
    cblpar[7] = kcbllh *0.5;
    cblpar[9] = cblpar[5];
    gMC->Gsvolu("FCBB", "TRAP", idtmed[512], cblpar, 11); // cables & tubes mix
    xcoor = (xtof - kcblw)*0.5 - 2.*sawpar[0];
    ycoor = (kcblh1 + 2.*cblpar[5])*0.25 - (ytof*0.5 - kHoneycombLayerThickness)*0.5;
    zcoor = kcbll-kcbllh*0.5;
    gMC->Gspos("FCBB", 1, "FAIB", -xcoor, ycoor, -zcoor, idrotm[7], "ONLY");
    gMC->Gspos("FCBB", 2, "FAIB",  xcoor, ycoor, -zcoor, idrotm[7], "ONLY");
    gMC->Gspos("FCBB", 3, "FAIB", -xcoor, ycoor,  zcoor, idrotm[8], "ONLY");
    gMC->Gspos("FCBB", 4, "FAIB",  xcoor, ycoor,  zcoor, idrotm[8], "ONLY");
    xcoor = xtof*0.5 - sawpar[0];
    ycoor = (kcblh2 - ytof*0.5 + kHoneycombLayerThickness)*0.5;
    gMC->Gspos("FSAW", 1, "FAIB", -xcoor, ycoor, 0., 0, "ONLY");
    gMC->Gspos("FSAW", 2, "FAIB",  xcoor, ycoor, 0., 0, "ONLY");
  }

  // TOF Supermodule cover definition and positioning
  Float_t covpar[3] = {xtof*0.5, 0.1, zlenA*0.5};
  gMC->Gsvolu("FCOV", "BOX ", idtmed[505], covpar, 3);    // Al cover
  xcoor = 0.;
  ycoor = 12.5*0.5 - 0.1;
  zcoor = 0.;
  gMC->Gspos("FCOV", 0, "FAIA", xcoor, ycoor, zcoor, 0, "ONLY");
  if (fTOFHoles) gMC->Gspos("FCOV", 0, "FAIB", xcoor, ycoor, zcoor, 0, "ONLY");

  // Services Volumes

  // Empty crate weight: 50 Kg, electronics cards + cables ~ 52 Kg.
  // Per each side (A and C) the total weight is: 2x102 ~ 204 Kg.
  // ... + weight of the connection pannel for the steel cooling system (Cr 18%, Ni 12%, Fe 70%)
  // + other remaining elements + various supports

  // Each FEA card weight + all supports
  // (including all bolts and not including the cable connectors)
  //  353.1 g.
  // Per each strip there are 4 FEA cards, then
  // the total weight of the front-end electonics section is: 353.1 g x 4 = 1412.4 g.

  Float_t serpar[3] = {29.*0.5, 121.*0.5, 90.*0.5};
  gMC->Gsvolu("FTOS", "BOX ", idtmed[515], serpar, 3); // Al + Cu + steel
  zcoor = (118.-90.)*0.5;
  Float_t phi = -10.,  ra = fTOFGeometry->Rmin() + ytof*0.5;
  for (Int_t i = 0; i < fTOFGeometry->NSectors(); i++) {
    phi += 20.;
    xcoor = ra * TMath::Cos(phi * kDegrad);
    ycoor = ra * TMath::Sin(phi * kDegrad);
    AliMatrix(idrotm[20+i], 90., phi, 90., phi + 270., 0., 0.);
    gMC->Gspos("FTOS", i, "BFMO", xcoor, ycoor, zcoor, idrotm[20+i], "ONLY");      
  }
  zcoor = (90. - 223.)*0.5;
  gMC->Gspos("FTOS", 1, "BBCE", ra, 0., zcoor, 0, "ONLY");

}
//_____________________________________________________________________________
void AliTOFv6T0::DrawModule() const
{
  //
  // Draw a shaded view of the Time Of Flight version 5
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

  char name[16];
  for (Int_t isec=0; isec<fTOFGeometry->NSectors(); isec++) {
    sprintf(name, "BREF%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BREF%d sub-levels skipped   -
    sprintf(name, "BTRD%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BTRD%d sub-levels skipped   -
    sprintf(name, "BTOF%d",isec);
    gMC->Gsatt(name,"seen",-2);  // all BTOF%d sub-levels skipped   -
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
void AliTOFv6T0::DrawDetectorModules() const
{
  //
  // Draw a shaded view of the TOF detector SuperModules version 5
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

  char name[16];
  for (Int_t isec=0; isec<fTOFGeometry->NSectors(); isec++) {
    sprintf(name, "BREF%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BREF%d sub-levels skipped   -
    sprintf(name, "BTRD%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BTRD%d sub-levels skipped   -
    sprintf(name, "BTOF%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BTOF%d sub-levels skipped   -
  }

  // Level 3 of B071, B075 and B074
  gMC->Gsatt("FTOA","seen",-2);  // all FTOA sub-levels skipped   -
  if (fTOFHoles) gMC->Gsatt("FTOB","seen",-2);  // all FTOB sub-levels skipped   -
  if (fTOFHoles) gMC->Gsatt("FTOC","seen",-2);  // all FTOC sub-levels skipped   -

  // Level 3 of B071, B075 and B074
  gMC->Gsatt("FAIA","seen",-1);  // all FAIA sub-levels skipped   -
  if (fTOFHoles) gMC->Gsatt("FAIB","seen",-1);  // all FAIB sub-levels skipped   -

  // Level 3 of B071, B075 and B074
  gMC->Gsatt("FPEA","seen",1);  // all FPEA sub-levels skipped   -
  if (fTOFHoles) gMC->Gsatt("FPEB","seen",1);  // all FPEB sub-levels skipped   -

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
void AliTOFv6T0::DrawDetectorStrips() const
{
  //
  // Draw a shaded view of the TOF strips for version 5
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

  char name[16];
  for (Int_t isec=0; isec<fTOFGeometry->NSectors(); isec++) {
    sprintf(name, "BREF%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BREF%d sub-levels skipped   -
    sprintf(name, "BTRD%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BTRD%d sub-levels skipped   -
    sprintf(name, "BTOF%d",isec);
    gMC->Gsatt(name,"seen", 0);  // all BTOF%d sub-levels skipped   -
  }

  // Level 3 of B071, B074 and B075
  gMC->Gsatt("FTOA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FTOB","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FTOC","SEEN", 0);

  // Level 4 of B071, B074 and B075
  gMC->Gsatt("FLTA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FLTB","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FLTC","SEEN", 0);

  // Level 5 of B071, B074 and B075
  gMC->Gsatt("FAIA","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FAIB","SEEN", 0);

  gMC->Gsatt("FPEA","SEEN", 1);
  if (fTOFHoles) gMC->Gsatt("FPEB","SEEN", 1);

  gMC->Gsatt("FSTR","SEEN",-2);  // all FSTR sub-levels skipped   -

  gMC->Gsatt("FWZ1","SEEN", 1);
  gMC->Gsatt("FWZ2","SEEN", 1);
  gMC->Gsatt("FWZ3","SEEN", 1);
  gMC->Gsatt("FWZ4","SEEN", 1);


  // Level 2 of FAIA
  // Level 2 of FAIB
  gMC->Gsatt("FCA1","SEEN", 0);
  gMC->Gsatt("FCA2","SEEN", 0);
  gMC->Gsatt("FCAB","SEEN", 0);
  gMC->Gsatt("FTUB","SEEN",-1);  // all FTUB sub-levels skipped   -
  gMC->Gsatt("FTLN","SEEN", 0);
  gMC->Gsatt("FLTN","SEEN", 0);
  gMC->Gsatt("FCBL","SEEN", 0);
  gMC->Gsatt("FSAW","SEEN", 0);
  gMC->Gsatt("FCOV","SEEN", 0);
  if (fTOFHoles) gMC->Gsatt("FCBB","SEEN", 0);

  // Level 2 of FTUB
  gMC->Gsatt("FITU","SEEN", 0);

  // Level 2 of FSTR
  gMC->Gsatt("FHON","SEEN", 1);
  gMC->Gsatt("FPC1","SEEN", 1);
  gMC->Gsatt("FPC2","SEEN", 1);
  gMC->Gsatt("FPCB","SEEN", 1);
  gMC->Gsatt("FRGL","SEEN", 1);
  gMC->Gsatt("FGLF","SEEN", 1);

  // Level 2 of FPCB => Level 3 of FSTR
  gMC->Gsatt("FSEN","SEEN", 0);
  gMC->Gsatt("FSEZ","SEEN", 0);
  gMC->Gsatt("FPAD","SEEN", 1);

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
void AliTOFv6T0::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //

  //AliTOF::CreateMaterials();

  AliMagF *magneticField = (AliMagF*)gAlice->Field();

  Int_t   isxfld = magneticField->Integ();
  Float_t sxmgmx = magneticField->Max();

  Float_t we[7], na[7];

  //--- Quartz (SiO2) to simulate float glass
  //    density tuned to have correct float glass 
  //    radiation length
  Float_t   aq[2] = { 28.09,16. };
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  //Float_t   dq = 2.55; // std value: 2.2
  Float_t   dq = 2.7;    // (+5.9%)
  Int_t nq = -2;

  // --- Nomex
  Float_t anox[4] = {12.01,1.01,16.00,14.01};
  Float_t znox[4] = { 6.,  1.,  8.,  7.};
  Float_t wnox[4] = {14., 22., 2., 2.};
  //Float_t dnox  = 0.048; //old value
  Float_t dnox  = 0.22;    // (x 4.6)
  Int_t nnox   = -4;

  // --- glass+freon { Si, O, C, F, H, S }
  Float_t agfr[6]= {28.09,16.00,12.01,19.00,1.01,32.065};
  Float_t zgfr[6]= {14.,  8.,  6.,  9.,  1.,  16.};
  Float_t wgfr[6]= {0.465, 0.530, 0.000484, 0.00383, 4.0e-05, 0.000646};
  Int_t ngfr  = 6;
  AliDebug(1,Form("wgfr: %d  %d  %d  %d  %d %d", wgfr[0], wgfr[1], wgfr[2], wgfr[3], wgfr[4], wgfr[5]));
  //Float_t dgfr = 1.35; // + FISHLINE (old value)
  Float_t dgfr = 1.6;    // + FISHLINE(+18.5 %)

  // --- G10  {Si, O, C, H, O}
  Float_t ag10[5] = {28.09,16.00,12.01,1.01,16.00};
  Float_t zg10[5] = {14., 8., 6., 1., 8.};
  Float_t wmatg10[5];
  Int_t nlmatg10 = 5;
  na[0]= 1. ,   na[1]= 2. ,   na[2]= 0. ,   na[3]= 0. ,   na[4]= 0.;
  MaterialMixer(we,ag10,na,5);
  wmatg10[0]= we[0]*0.6;
  wmatg10[1]= we[1]*0.6;
  na[0]= 0. ,   na[1]= 0. ,   na[2]= 14. ,   na[3]= 20. ,   na[4]= 3.;
  MaterialMixer(we,ag10,na,5);
  wmatg10[2]= we[2]*0.4;
  wmatg10[3]= we[3]*0.4;
  wmatg10[4]= we[4]*0.4;
  AliDebug(1,Form("wg10  %d  %d  %d  %d  %d", wmatg10[0], wmatg10[1], wmatg10[2], wmatg10[3], wmatg10[4]));
  //  Float_t densg10 = 1.7; //old value
  Float_t densg10 = 2.0; // (+17.8%)

  // -- Water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;
  Int_t nwa = -2;

  // AIR
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir   = 1.20479E-3;

  // --- fibre glass
  Float_t afg[4] = {28.09,16.00,12.01,1.01};
  Float_t zfg[4] = {14., 8., 6., 1.};
  Float_t wfg[4] = {0.12906,0.29405,0.51502,0.06187};
  //Float_t dfg    = 1.111;
  Float_t dfg    = 2.; // (+1.8%)
  Int_t nfg      = 4;

  // --- Freon C2F4H2 + SF6
  Float_t afre[4]= {12.01,1.01,19.00,32.07};
  Float_t zfre[4]= { 6., 1., 9., 16.};
  Float_t wfre[4]= {0.21250,0.01787,0.74827,0.021355};
  Float_t densfre= 0.00375;
  Int_t nfre     = 4;

  // --- Al + Cu + G10  {Al, Cu, Si, O, C, H, O}
  Float_t acar[10]= {26.98,
		     /*63.55,*/
		     ag10[0], ag10[1], ag10[2], ag10[3], ag10[4],
		     aAir[0], aAir[1], aAir[2], aAir[3]};
  Float_t zcar[10]= {13.,
		     /*29.,*/
		     zg10[0], zg10[1], zg10[2], zg10[3], zg10[4],
		     zAir[0], zAir[1], zAir[2], zAir[3]};
  Float_t wcar[10];
  wcar[0]= 0.4732;//0.7;
  //wcar[1]= 0.04;//0.05;
  wcar[1]= 0.2854*wmatg10[0];//0.25*wmatg10[0];
  wcar[2]= 0.2854*wmatg10[1];//0.25*wmatg10[1];
  wcar[3]= 0.2854*wmatg10[2];//0.25*wmatg10[2];
  wcar[4]= 0.2854*wmatg10[3];//0.25*wmatg10[3];
  wcar[5]= 0.2854*wmatg10[4];//0.25*wmatg10[4];
  wcar[6]= 0.2414*wAir[0];
  wcar[7]= 0.2414*wAir[1];
  wcar[8]= 0.2414*wAir[2];
  wcar[9]= 0.2414*wAir[3];

  AliDebug(1,Form("wcar  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f", wcar[0], wcar[1], wcar[2], wcar[3], wcar[4],
		  wcar[5], wcar[6], wcar[7], wcar[8], wcar[9]));
  Float_t dcar = 1.85;//1.9;

  // --- Cables, tubes {Al, Cu} ---
  Float_t acbt[2]= {26.98,63.55};
  Float_t zcbt[2]= {13., 29.};
  //Float_t wcbt[2]= {0.541,0.459};
  Float_t wcbt[2]= {0.407,0.593};
  //Float_t decbt  = 0.95;
  Float_t decbt  = 0.68;

  // --- Cable {Al, Cu}
  Float_t wcb[2] = {0.165,0.835};
  Float_t decb   = 0.962;

  // --- Honeycomb layer {Al, Cu}
  Float_t whon[2]= {0.9,0.1};
  //Float_t dhon   = 0.44;
  Float_t dhon   = 1.095; // (x 2.56)

  // --- Crates boxes {Al, Cu, Fe, Cr, Ni}
  Float_t acra[5]= {26.98,63.55,55.845,52.00,58.69};
  Float_t zcra[5]= {13., 29., 26., 24., 28.};
  Float_t wcra[5]= {0.7,0.2,0.07,0.018,0.012};
  Float_t dcra   = 0.77;

  AliMixture ( 0, "Air$", aAir, zAir, dAir, 4, wAir);
  AliMixture ( 1, "Nomex$", anox, znox, dnox, nnox, wnox);
  AliMixture ( 2, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  AliMixture ( 3, "fibre glass$", afg, zfg, dfg, nfg, wfg);
  AliMaterial( 4, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  AliMixture ( 5, "Al+Cu honeycomb$", acbt, zcbt, dhon, 2, whon);
  AliMixture ( 6, "Freon$", afre, zfre, densfre, nfre, wfre);
  AliMixture ( 7, "Glass$", aq, zq, dq, nq, wq);
  AliMixture ( 8, "glass-freon$", agfr, zgfr, dgfr, ngfr, wgfr);
  AliMixture ( 9, "Water$",  awa, zwa, dwa, nwa, wwa);
  AliMixture (10, "Al+Cu$", acbt, zcbt, decbt, 2, wcbt);
  AliMaterial(11, "Cu $", 63.54, 29., 8.96, 1.43, 10.);
  AliMixture (12, "Al+Cu (cable)$", acbt, zcbt, decb, 2, wcb);
  AliMixture (13, "Al+Cu+G10$", acar, zcar, dcar, 10/*7*/, wcar);
  AliMixture (14, "Al+Cu+steel$", acra, zcra, dcra, 5, wcra);
  AliMaterial(15, "Cu_sensitive$", 63.54, 29., 3.392, 1.43, 10.);

  Float_t epsil, stmin, deemax, stemax;

  //   STD data
  //  EPSIL  = 0.1   ! Tracking precision,
  //  STEMAX = 0.1   ! Maximum displacement for multiple scattering
  //  DEEMAX = 0.1   ! Maximum fractional energy loss, DLS
  //  STMIN  = 0.1

  // TOF data
  epsil  = .001;  // Tracking precision,
  stemax = -1.;   // Maximum displacement for multiple scattering
  deemax = -.3;   // Maximum fractional energy loss, DLS
  stmin  = -.8;

  AliMedium( 1, "Air$",         0, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 2,"Nomex$",        1, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 3,"G10$",          2, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 4,"fibre glass$",  3, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 5,"glass-freon$",  8, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 6,"Al Frame$",     4, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 7,"honeycomb$",    5, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 8,"Fre$",          6, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 9,"Cu-S$",        15, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(10,"Glass$",        7, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(11,"Water$",        9, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(12,"Cable$",       12, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(13,"Al+Cables$",   10, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(14,"Copper$",      11, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(15,"Cards$",       13, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(16,"Crates$",      14, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);

}
//_____________________________________________________________________________
void AliTOFv6T0::Init()
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
void AliTOFv6T0::StepManager()
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

    AddTrackReference(mcApplication->GetCurrentTrackNumber(), AliTrackReference::kTOF);
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
void AliTOFv6T0::MaterialMixer(Float_t* p,Float_t* a,Float_t* m,Int_t n) const
{
  // a[] atomic weights vector      (in)
  //     (atoms present in more compound appear separately)
  // m[] number of corresponding atoms in the compound  (in)
  Float_t t = 0.;
  for (Int_t i = 0; i < n; ++i) {
    p[i] = a[i]*m[i];
    t  += p[i];
  }
  for (Int_t i = 0; i < n; ++i) {
    p[i] = p[i]/t;
    //AliDebug(1,Form((\n weight[%i] = %f (,i,p[i]));
  }
}
