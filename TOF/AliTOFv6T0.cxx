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
Revision 1.11  2007/10/08 17:52:55  decaro
hole region in front of PHOS detector: update of sectors' numbers

Revision 1.10  2007/10/07 19:40:46  decaro
right handling of l2t matrices and alignable entries in case of TOF staging geometry

Revision 1.9  2007/10/07 19:36:29  decaro
TOF materials and volumes description: update

Revision 1.8  2007/10/04 13:15:37  arcelli
updates to comply with AliTOFGeometryV5 becoming AliTOFGeometry

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

#include <TDirectory.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TGeoVolume.h>
#include <TLorentzVector.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliTrackReference.h"

#include "AliTOFGeometry.h"
#include "AliTOFv6T0.h"

extern TVirtualMC *gMC;
extern TGeoManager *gGeoManager;

extern AliRun *gAlice;

ClassImp(AliTOFv6T0)

// TOF sectors with Nino masks: 0, 8, 9, 10, 16
const Bool_t AliTOFv6T0::fgkFEAwithMasks[18] = 
{kTRUE , kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
 kFALSE, kFALSE, kTRUE , kTRUE , kTRUE , kFALSE,
 kFALSE, kFALSE, kFALSE, kFALSE, kTRUE , kFALSE};
const Float_t AliTOFv6T0::fgkModuleWallThickness   =   0.33; // cm
const Float_t AliTOFv6T0::fgkInterCentrModBorder1  =  49.5 ; // cm
const Float_t AliTOFv6T0::fgkInterCentrModBorder2  =  57.5 ; // cm
const Float_t AliTOFv6T0::fgkExterInterModBorder1  = 196.0 ; // cm
const Float_t AliTOFv6T0::fgkExterInterModBorder2  = 203.5 ; // cm
//const Float_t AliTOFv6T0::fgkLengthInCeModBorder   =   7.2 ; // cm // it was 4.7 cm (AdC)
const Float_t AliTOFv6T0::fgkLengthInCeModBorderU  =   5.0 ; // cm
const Float_t AliTOFv6T0::fgkLengthInCeModBorderD  =   7.0 ; // cm
const Float_t AliTOFv6T0::fgkLengthExInModBorder   =   5.0 ; // cm // it was 7.0 cm (AdC)
const Float_t AliTOFv6T0::fgkModuleCoverThickness  =   2.0 ; // cm
const Float_t AliTOFv6T0::fgkFEAwidth1    = 19.0; // cm
const Float_t AliTOFv6T0::fgkFEAwidth2    = 39.5;//38.5; // cm
const Float_t AliTOFv6T0::fgkSawThickness =  1.0; // cm
const Float_t AliTOFv6T0::fgkCBLw  = 13.5; // cm
const Float_t AliTOFv6T0::fgkCBLh1 =  2.0; // cm
const Float_t AliTOFv6T0::fgkCBLh2 = 12.3; // cm
const Float_t AliTOFv6T0::fgkBetweenLandMask = 0.1; // cm
const Float_t AliTOFv6T0::fgkAl1parameters[3] = {fgkFEAwidth1*0.5, 0.4, 0.2}; // cm
const Float_t AliTOFv6T0::fgkAl2parameters[3] = {7.25, 0.75, 0.25}; // cm
const Float_t AliTOFv6T0::fgkAl3parameters[3] = {3., 4., 0.1}; // cm
const Float_t AliTOFv6T0::fgkRoof1parameters[3] = {fgkAl1parameters[0], fgkAl1parameters[2], 1.45}; // cm
const Float_t AliTOFv6T0::fgkRoof2parameters[3] = {fgkAl3parameters[0], 0.1, 1.15}; // cm
const Float_t AliTOFv6T0::fgkFEAparameters[3] = {fgkFEAwidth1*0.5, 5.6, 0.1}; // cm
const Float_t AliTOFv6T0::fgkBar[3] = {8.575, 0.6, 0.25}; // cm
const Float_t AliTOFv6T0::fgkBar1[3] = {fgkBar[0], fgkBar[1], 0.1}; // cm
const Float_t AliTOFv6T0::fgkBar2[3] = {fgkBar[0], 0.1, fgkBar[1] - 2.*fgkBar1[2]}; // cm
const Float_t AliTOFv6T0::fgkBarS[3] = {2., fgkBar[1], fgkBar[2]}; // cm
const Float_t AliTOFv6T0::fgkBarS1[3] = {fgkBarS[0], fgkBar1[1], fgkBar1[2]}; // cm
const Float_t AliTOFv6T0::fgkBarS2[3] = {fgkBarS[0], fgkBar2[1], fgkBar2[2]}; // cm

//_____________________________________________________________________________
  AliTOFv6T0::AliTOFv6T0():
  fIdFTOA(-1),
  fIdFTOB(-1),
  fIdFTOC(-1),
  fIdFLTA(-1),
  fIdFLTB(-1),
  fIdFLTC(-1)//,
//fTOFHoles(kFALSE)
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
  fIdFLTC(-1)//,
  //fTOFHoles(kFALSE)
{
  //
  // Standard constructor
  //

  //
  // Check that FRAME is there otherwise we have no place where to
  // put TOF

  /*
  AliModule* frame = (AliModule*)gAlice->GetModule("FRAME");

  if(!frame) {
    AliFatal("TOF needs FRAME to be present");
  } else {
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
  */

  if (fTOFGeometry) delete fTOFGeometry;
  fTOFGeometry = new AliTOFGeometry();
  fTOFGeometry->SetHoles(fTOFHoles);

  //AliTOF::fTOFGeometry = fTOFGeometry;

  // Save the geometry
  TDirectory* saveDir = gDirectory;
  AliRunLoader::Instance()->CdGAFile();
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

  AliGeomManager::ELayerID idTOF = AliGeomManager::kTOF;
  Int_t modUID, modnum=0;

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

      modUID = AliGeomManager::LayerToVolUID(idTOF, modnum++);
      if (fTOFSectors[isect]==-1) continue;

      if (fTOFHoles && (isect==13 || isect==14 || isect==15)) {
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
	      
      if(!gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data(),modUID))
	AliError(Form("Alignable entry %s not set",symName.Data()));

      //T2L matrices for alignment
      TGeoPNEntry *e = gGeoManager->GetAlignableEntryByUID(modUID);
      if (e) {
	TGeoHMatrix *globMatrix = e->GetGlobalOrig();
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

  AliDebug(1, "************************* TOF geometry **************************");
  AliDebug(1,Form(" xtof   %f",  xtof));
  AliDebug(1,Form(" ytof   %f",  ytof));
  AliDebug(1,Form(" zlenA   %f", zlenA));
  AliDebug(2,Form(" zlenA*0.5 = %f", zlenA*0.5));

  Float_t xFLT, yFLT, zFLTA;
  xFLT  = xtof     - 2.*fgkModuleWallThickness;
  yFLT  = ytof*0.5 -    fgkModuleWallThickness;
  zFLTA = zlenA    - 2.*fgkModuleWallThickness;

  CreateModules(xtof, ytof, zlenA, xFLT, yFLT, zFLTA);
  MakeStripsInModules(ytof, zlenA);

  CreateModuleCovers(xtof, zlenA);

  CreateBackZone(xtof, ytof, zlenA);
  MakeFrontEndElectronics(xtof);
  MakeFEACooling(xtof);
  MakeNinoMask(xtof);
  MakeSuperModuleCooling(xtof, ytof, zlenA);
  MakeSuperModuleServices(xtof, ytof, zlenA);

  MakeModulesInBTOFvolumes(ytof, zlenA);
  MakeCoversInBTOFvolumes();
  MakeBackInBTOFvolumes(ytof);

  MakeReadoutCrates(ytof);

}

//_____________________________________________________________________________
void AliTOFv6T0::CreateModules(Float_t xtof,  Float_t ytof, Float_t zlenA,
			       Float_t xFLT,  Float_t yFLT, Float_t zFLTA) const
{
  //
  // Create supermodule volume
  // and wall volumes to separate 5 modules
  //

  const Float_t kPi = TMath::Pi();

  Int_t *idtmed = fIdtmed->GetArray()-499;

  Int_t idrotm[8];

  // Definition of the of fibre glass modules (FTOA, FTOB and FTOC)
  Float_t  par[3];
  par[0] = xtof * 0.5;
  par[1] = ytof * 0.25;
  par[2] = zlenA * 0.5;
  gMC->Gsvolu("FTOA", "BOX ", idtmed[503], par, 3);  // Fibre glass

  if (fTOFHoles) {
    par[0] =  xtof * 0.5;
    par[1] =  ytof * 0.25;
    par[2] = (zlenA*0.5 - fgkInterCentrModBorder1)*0.5;
    gMC->Gsvolu("FTOB", "BOX ", idtmed[503], par, 3);  // Fibre glass
    gMC->Gsvolu("FTOC", "BOX ", idtmed[503], par, 3);  // Fibre glass
  }


  // Definition and positioning
  // of the not sensitive volumes with Insensitive Freon (FLTA, FLTB and FLTC)
  par[0] = xFLT*0.5;
  par[1] = yFLT*0.5;
  par[2] = zFLTA*0.5;
  gMC->Gsvolu("FLTA", "BOX ", idtmed[506], par, 3); // Freon mix

  Float_t xcoor, ycoor, zcoor;
  xcoor = 0.;
  ycoor = fgkModuleWallThickness*0.5;
  zcoor = 0.;
  gMC->Gspos ("FLTA", 0, "FTOA", xcoor, ycoor, zcoor, 0, "ONLY");

  if (fTOFHoles) {
    par[2] = (zlenA*0.5 - 2.*fgkModuleWallThickness - fgkInterCentrModBorder1)*0.5;
    gMC->Gsvolu("FLTB", "BOX ", idtmed[506], par, 3); // Freon mix
    gMC->Gsvolu("FLTC", "BOX ", idtmed[506], par, 3); // Freon mix

    //xcoor = 0.;
    //ycoor = fgkModuleWallThickness*0.5;
    zcoor = fgkModuleWallThickness;
    gMC->Gspos ("FLTB", 0, "FTOB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos ("FLTC", 0, "FTOC", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  // Definition and positioning
  // of the fibre glass walls between central and intermediate modules (FWZ1 and FWZ2)
  Float_t alpha, tgal, beta, tgbe, trpa[11];
  //tgal  = (yFLT - 2.*fgkLengthInCeModBorder)/(fgkInterCentrModBorder2 - fgkInterCentrModBorder1);
  tgal  = (yFLT - fgkLengthInCeModBorderU - fgkLengthInCeModBorderD)/(fgkInterCentrModBorder2 - fgkInterCentrModBorder1);
  alpha = TMath::ATan(tgal);
  beta  = (kPi*0.5 - alpha)*0.5;
  tgbe  = TMath::Tan(beta);
  trpa[0]  = xFLT*0.5;
  trpa[1]  = 0.;
  trpa[2]  = 0.;
  trpa[3]  = 2.*fgkModuleWallThickness;
  //trpa[4]  = (fgkLengthInCeModBorder - 2.*fgkModuleWallThickness*tgbe)*0.5;
  //trpa[5]  = (fgkLengthInCeModBorder + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[4]  = (fgkLengthInCeModBorderD - 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[5]  = (fgkLengthInCeModBorderD + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  trpa[7]  = 2.*fgkModuleWallThickness;
  trpa[8]  = (fgkLengthInCeModBorderD - 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[9]  = (fgkLengthInCeModBorderD + 2.*fgkModuleWallThickness*tgbe)*0.5;
  //trpa[8]  = (fgkLengthInCeModBorder - 2.*fgkModuleWallThickness*tgbe)*0.5;
  //trpa[9]  = (fgkLengthInCeModBorder + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  gMC->Gsvolu("FWZ1D", "TRAP", idtmed[503], trpa, 11); // Fibre glass

  AliMatrix (idrotm[0],90., 90.,180.,0.,90.,180.);
  AliMatrix (idrotm[1],90., 90.,  0.,0.,90.,  0.);

  //xcoor = 0.;
  //ycoor = -(yFLT - fgkLengthInCeModBorder)*0.5;
  ycoor = -(yFLT - fgkLengthInCeModBorderD)*0.5;
  zcoor = fgkInterCentrModBorder1;
  gMC->Gspos("FWZ1D", 1, "FLTA", xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FWZ1D", 2, "FLTA", xcoor, ycoor,-zcoor, idrotm[1], "ONLY");

  Float_t y0B, ycoorB, zcoorB;

  if (fTOFHoles) {
    //y0B = fgkLengthInCeModBorder - fgkModuleWallThickness*tgbe;
    y0B = fgkLengthInCeModBorderD - fgkModuleWallThickness*tgbe;
    trpa[0]  = xFLT*0.5;
    trpa[1]  = 0.;
    trpa[2]  = 0.;
    trpa[3]  = fgkModuleWallThickness;
    trpa[4]  = (y0B - fgkModuleWallThickness*tgbe)*0.5;
    trpa[5]  = (y0B + fgkModuleWallThickness*tgbe)*0.5;
    trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    trpa[7]  = fgkModuleWallThickness;
    trpa[8]  = (y0B - fgkModuleWallThickness*tgbe)*0.5;
    trpa[9]  = (y0B + fgkModuleWallThickness*tgbe)*0.5;
    trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    //xcoor = 0.;
    ycoorB = ycoor - fgkModuleWallThickness*0.5*tgbe;
    zcoorB = (zlenA*0.5 - 2.*fgkModuleWallThickness - fgkInterCentrModBorder1)*0.5 - 2.*fgkModuleWallThickness;
    gMC->Gsvolu("FWZAD", "TRAP", idtmed[503], trpa, 11); // Fibre glass
    gMC->Gspos("FWZAD", 1, "FLTB", xcoor, ycoorB, zcoorB, idrotm[1], "ONLY");
    gMC->Gspos("FWZAD", 2, "FLTC", xcoor, ycoorB,-zcoorB, idrotm[0], "ONLY");
  }



  tgal  = (yFLT - fgkLengthInCeModBorderU - fgkLengthInCeModBorderD)/(fgkInterCentrModBorder2 - fgkInterCentrModBorder1);
  alpha = TMath::ATan(tgal);
  beta  = (kPi*0.5 - alpha)*0.5;
  tgbe  = TMath::Tan(beta);
  trpa[0]  = xFLT*0.5;
  trpa[1]  = 0.;
  trpa[2]  = 0.;
  trpa[3]  = 2.*fgkModuleWallThickness;
  //trpa[4]  = (fgkLengthInCeModBorder - 2.*fgkModuleWallThickness*tgbe)*0.5;
  //trpa[5]  = (fgkLengthInCeModBorder + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[4]  = (fgkLengthInCeModBorderU - 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[5]  = (fgkLengthInCeModBorderU + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  trpa[7]  = 2.*fgkModuleWallThickness;
  trpa[8]  = (fgkLengthInCeModBorderU - 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[9]  = (fgkLengthInCeModBorderU + 2.*fgkModuleWallThickness*tgbe)*0.5;
  //trpa[8]  = (fgkLengthInCeModBorder - 2.*fgkModuleWallThickness*tgbe)*0.5;
  //trpa[9]  = (fgkLengthInCeModBorder + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  gMC->Gsvolu("FWZ1U", "TRAP", idtmed[503], trpa, 11); // Fibre glass


  AliMatrix (idrotm[2],90.,270.,  0.,0.,90.,180.);
  AliMatrix (idrotm[3],90.,270.,180.,0.,90.,  0.);

  //xcoor = 0.;
  //ycoor = (yFLT - fgkLengthInCeModBorder)*0.5;
  ycoor = (yFLT - fgkLengthInCeModBorderU)*0.5;
  zcoor = fgkInterCentrModBorder2;
  gMC->Gspos("FWZ1U", 1, "FLTA", xcoor, ycoor, zcoor,idrotm[2], "ONLY");
  gMC->Gspos("FWZ1U", 2, "FLTA", xcoor, ycoor,-zcoor,idrotm[3], "ONLY");

  if (fTOFHoles) {
    //y0B = fgkLengthInCeModBorder + fgkModuleWallThickness*tgbe;
    y0B = fgkLengthInCeModBorderU + fgkModuleWallThickness*tgbe;
    trpa[0]  = xFLT*0.5;
    trpa[1]  = 0.;
    trpa[2]  = 0.;
    trpa[3]  = fgkModuleWallThickness;
    trpa[4]  = (y0B - fgkModuleWallThickness*tgbe)*0.5;
    trpa[5]  = (y0B + fgkModuleWallThickness*tgbe)*0.5;
    trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    trpa[7]  = fgkModuleWallThickness;
    trpa[8]  = (y0B - fgkModuleWallThickness*tgbe)*0.5;
    trpa[9]  = (y0B + fgkModuleWallThickness*tgbe)*0.5;
    trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
    gMC->Gsvolu("FWZBU", "TRAP", idtmed[503], trpa, 11); // Fibre glass
    //xcoor = 0.;
    ycoorB = ycoor - fgkModuleWallThickness*0.5*tgbe;
    zcoorB = (zlenA*0.5 - 2.*fgkModuleWallThickness - fgkInterCentrModBorder1)*0.5 -
      (fgkInterCentrModBorder2 - fgkInterCentrModBorder1) - 2.*fgkModuleWallThickness;
    gMC->Gspos("FWZBU", 1, "FLTB", xcoor, ycoorB, zcoorB, idrotm[3], "ONLY");
    gMC->Gspos("FWZBU", 2, "FLTC", xcoor, ycoorB,-zcoorB, idrotm[2], "ONLY");
  }

  trpa[0] = 0.5*(fgkInterCentrModBorder2 - fgkInterCentrModBorder1)/TMath::Cos(alpha);
  trpa[1] = 2.*fgkModuleWallThickness;
  trpa[2] = xFLT*0.5;
  trpa[3] = -beta*kRaddeg;
  trpa[4] = 0.;
  trpa[5] = 0.;
  gMC->Gsvolu("FWZ2", "PARA", idtmed[503], trpa, 6); // Fibre glass

  AliMatrix (idrotm[4],     alpha*kRaddeg,90.,90.+alpha*kRaddeg,90.,90.,180.);
  AliMatrix (idrotm[5],180.-alpha*kRaddeg,90.,90.-alpha*kRaddeg,90.,90.,  0.);

  //xcoor = 0.;
  //ycoor = 0.;
  ycoor = (fgkLengthInCeModBorderD - fgkLengthInCeModBorderU)*0.5;
  zcoor = (fgkInterCentrModBorder2 + fgkInterCentrModBorder1)*0.5;
  gMC->Gspos("FWZ2", 1, "FLTA", xcoor, ycoor, zcoor, idrotm[4], "ONLY");
  gMC->Gspos("FWZ2", 2, "FLTA", xcoor, ycoor,-zcoor, idrotm[5], "ONLY");

  if (fTOFHoles) {
    trpa[0] = 0.5*(fgkInterCentrModBorder2 - fgkInterCentrModBorder1)/TMath::Cos(alpha);
    trpa[1] = fgkModuleWallThickness;
    trpa[2] = xFLT*0.5;
    trpa[3] = -beta*kRaddeg;
    trpa[4] = 0.;
    trpa[5] = 0.;
    gMC->Gsvolu("FWZC", "PARA", idtmed[503], trpa, 6); // Fibre glass
    //xcoor = 0.;
    ycoorB = ycoor - fgkModuleWallThickness*tgbe;
    zcoorB = (zlenA*0.5 - 2.*fgkModuleWallThickness - fgkInterCentrModBorder1)*0.5 -
      (fgkInterCentrModBorder2 - fgkInterCentrModBorder1)*0.5 - 2.*fgkModuleWallThickness;
    gMC->Gspos("FWZC", 1, "FLTB", xcoor, ycoorB, zcoorB, idrotm[5], "ONLY");
    gMC->Gspos("FWZC", 2, "FLTC", xcoor, ycoorB,-zcoorB, idrotm[4], "ONLY");
  }


  // Definition and positioning
  // of the fibre glass walls between intermediate and lateral modules (FWZ3 and FWZ4)
  tgal  = (yFLT - 2.*fgkLengthExInModBorder)/(fgkExterInterModBorder2 - fgkExterInterModBorder1);
  alpha = TMath::ATan(tgal);
  beta  = (kPi*0.5 - alpha)*0.5;
  tgbe  = TMath::Tan(beta);
  trpa[0]  = xFLT*0.5;
  trpa[1]  = 0.;
  trpa[2]  = 0.;
  trpa[3]  = 2.*fgkModuleWallThickness;
  trpa[4]  = (fgkLengthExInModBorder - 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[5]  = (fgkLengthExInModBorder + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[6]  = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  trpa[7]  = 2.*fgkModuleWallThickness;
  trpa[8]  = (fgkLengthExInModBorder - 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[9]  = (fgkLengthExInModBorder + 2.*fgkModuleWallThickness*tgbe)*0.5;
  trpa[10] = TMath::ATan(tgbe*0.5)*kRaddeg; //TMath::ATan((trpa[5] - trpa[4])/(2.*trpa[3]))*kRaddeg;
  gMC->Gsvolu("FWZ3", "TRAP", idtmed[503], trpa, 11); // Fibre glass

  //xcoor = 0.;
  ycoor = (yFLT - fgkLengthExInModBorder)*0.5;
  zcoor = fgkExterInterModBorder1;
  gMC->Gspos("FWZ3", 1, "FLTA", xcoor, ycoor, zcoor,idrotm[3], "ONLY");
  gMC->Gspos("FWZ3", 2, "FLTA", xcoor, ycoor,-zcoor,idrotm[2], "ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = (yFLT - fgkLengthExInModBorder)*0.5;
    zcoor = -fgkExterInterModBorder1 + (zlenA*0.5 + fgkInterCentrModBorder1 - 2.*fgkModuleWallThickness)*0.5;
    gMC->Gspos("FWZ3", 5, "FLTB", xcoor, ycoor, zcoor, idrotm[2], "ONLY");
    gMC->Gspos("FWZ3", 6, "FLTC", xcoor, ycoor,-zcoor, idrotm[3], "ONLY");
  }

  //xcoor = 0.;
  ycoor = -(yFLT - fgkLengthExInModBorder)*0.5;
  zcoor = fgkExterInterModBorder2;
  gMC->Gspos("FWZ3", 3, "FLTA", xcoor, ycoor, zcoor, idrotm[1], "ONLY");
  gMC->Gspos("FWZ3", 4, "FLTA", xcoor, ycoor,-zcoor, idrotm[0], "ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = -(yFLT - fgkLengthExInModBorder)*0.5;
    zcoor = -fgkExterInterModBorder2 + (zlenA*0.5 + fgkInterCentrModBorder1 - 2.*fgkModuleWallThickness)*0.5;
    gMC->Gspos("FWZ3", 7, "FLTB", xcoor, ycoor, zcoor, idrotm[0], "ONLY");
    gMC->Gspos("FWZ3", 8, "FLTC", xcoor, ycoor,-zcoor, idrotm[1], "ONLY");
  }

  trpa[0] = 0.5*(fgkExterInterModBorder2 - fgkExterInterModBorder1)/TMath::Cos(alpha);
  trpa[1] = 2.*fgkModuleWallThickness;
  trpa[2] = xFLT*0.5;
  trpa[3] = -beta*kRaddeg;
  trpa[4] = 0.;
  trpa[5] = 0.;
  gMC->Gsvolu("FWZ4", "PARA", idtmed[503], trpa, 6); // Fibre glass

  AliMatrix (idrotm[6],alpha*kRaddeg,90.,90.+alpha*kRaddeg,90.,90.,180.);
  AliMatrix (idrotm[7],180.-alpha*kRaddeg,90.,90.-alpha*kRaddeg,90.,90.,0.);

  //xcoor = 0.;
  ycoor = 0.;
  zcoor = (fgkExterInterModBorder2 + fgkExterInterModBorder1)*0.5;
  gMC->Gspos("FWZ4", 1, "FLTA", xcoor, ycoor, zcoor, idrotm[7], "ONLY");
  gMC->Gspos("FWZ4", 2, "FLTA", xcoor, ycoor,-zcoor, idrotm[6], "ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = 0.;
    zcoor = -(fgkExterInterModBorder2 + fgkExterInterModBorder1)*0.5 +
      (zlenA*0.5 + fgkInterCentrModBorder1 - 2.*fgkModuleWallThickness)*0.5;
    gMC->Gspos("FWZ4", 3, "FLTB", xcoor, ycoor, zcoor, idrotm[6], "ONLY");
    gMC->Gspos("FWZ4", 4, "FLTC", xcoor, ycoor,-zcoor, idrotm[7], "ONLY");
  }

}

//_____________________________________________________________________________
void AliTOFv6T0::CreateModuleCovers(Float_t xtof, Float_t zlenA) const
{
  //
  // Create covers for module:
  //   per each module zone, defined according to
  //   fgkInterCentrModBorder2, fgkExterInterModBorder1 and zlenA+2 values,
  //   there is a frame of thickness 2cm in Al
  //   and the contained zones in honeycomb of Al.
  //   There is also an interface layer (1.6mm thichness)
  //   and plastic and Cu corresponding to the flat cables.
  //

  Int_t  *idtmed = fIdtmed->GetArray()-499;

  Float_t par[3];
  par[0] = xtof*0.5 + 2.;
  par[1] = fgkModuleCoverThickness*0.5;
  par[2] = zlenA*0.5 + 2.;
  gMC->Gsvolu("FPEA", "BOX ", idtmed[500], par, 3); // Air
  if (fTOFHoles) gMC->Gsvolu("FPEB", "BOX ", idtmed[500], par, 3); // Air

  const Float_t kAlCoverThickness = 1.5;
  const Float_t kInterfaceCardThickness = 0.16;
  const Float_t kAlSkinThickness = 0.1;

  //par[0] = xtof*0.5 + 2.;
  par[1] = kAlCoverThickness*0.5;
  //par[2] = zlenA*0.5 + 2.;
  gMC->Gsvolu("FALT", "BOX ", idtmed[504], par, 3); // Al
  if (fTOFHoles) gMC->Gsvolu("FALB", "BOX ", idtmed[504], par, 3); // Al
  Float_t  xcoor, ycoor, zcoor;
  xcoor = 0.;
  ycoor = 0.;
  zcoor = 0.;
  gMC->Gspos("FALT", 0, "FPEA", xcoor, ycoor, zcoor, 0, "ONLY");
  if (fTOFHoles) gMC->Gspos("FALB", 0, "FPEB", xcoor, ycoor, zcoor, 0, "ONLY");

  par[0] = xtof*0.5;
  //par[1] = kAlCoverThickness*0.5;
  par[2] = fgkInterCentrModBorder2 - 2.;
  gMC->Gsvolu("FPE1", "BOX ", idtmed[505], par, 3); // Al honeycomb
  //xcoor = 0.;
  //ycoor = 0.;
  //zcoor = 0.;
  gMC->Gspos("FPE1", 0, "FALT", xcoor, ycoor, zcoor, 0, "ONLY");

  if (fTOFHoles) {
    //par[0] = xtof*0.5;
    par[1] = kAlCoverThickness*0.5 - kAlSkinThickness;
    //par[2] = fgkInterCentrModBorder2 - 2.;
    gMC->Gsvolu("FPE4", "BOX ", idtmed[515], par, 3); // Al honeycomb for holes
    //xcoor = 0.;
    //ycoor = 0.;
    //zcoor = 0.;
    gMC->Gspos("FPE4", 0, "FALB", xcoor, ycoor, zcoor, 0, "ONLY");
  }

  //par[0] = xtof*0.5;
  //par[1] = kAlCoverThickness*0.5;
  par[2] = (fgkExterInterModBorder1 - fgkInterCentrModBorder2)*0.5 - 2.;
  gMC->Gsvolu("FPE2", "BOX ", idtmed[505], par, 3); // Al honeycomb
  //xcoor = 0.;
  //ycoor = 0.;
  zcoor = (fgkExterInterModBorder1 + fgkInterCentrModBorder2)*0.5;
  gMC->Gspos("FPE2", 1, "FALT", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FPE2", 2, "FALT", xcoor, ycoor,-zcoor, 0, "ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = 0.;
    //zcoor = (fgkExterInterModBorder1 + fgkInterCentrModBorder2)*0.5;
    gMC->Gspos("FPE2", 1, "FALB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos("FPE2", 2, "FALB", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  //par[0] = xtof*0.5;
  //par[1] = kAlCoverThickness*0.5;
  par[2] = (zlenA*0.5 + 2. - fgkExterInterModBorder1)*0.5 - 2.;
  gMC->Gsvolu("FPE3", "BOX ", idtmed[505], par, 3); // Al honeycomb
  //xcoor = 0.;
  //ycoor = 0.;
  zcoor = (zlenA*0.5 + 2. + fgkExterInterModBorder1)*0.5;
  gMC->Gspos("FPE3", 1, "FALT", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FPE3", 2, "FALT", xcoor, ycoor,-zcoor, 0, "ONLY");

  if (fTOFHoles) {
    //xcoor = 0.;
    //ycoor = 0.;
    zcoor = (zlenA*0.5 + 2. + fgkExterInterModBorder1)*0.5;
    gMC->Gspos("FPE3", 1, "FALB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos("FPE3", 2, "FALB", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  // volumes for Interface cards
  par[0] = xtof*0.5;
  par[1] = kInterfaceCardThickness*0.5;
  par[2] = fgkInterCentrModBorder2 - 2.;
  gMC->Gsvolu("FIF1", "BOX ", idtmed[502], par, 3); // G10
  //xcoor = 0.;
  ycoor = kAlCoverThickness*0.5 + kInterfaceCardThickness*0.5;
  zcoor = 0.;
  gMC->Gspos("FIF1", 0, "FPEA", xcoor, ycoor, zcoor, 0, "ONLY");

  //par[0] = xtof*0.5;
  //par[1] = kInterfaceCardThickness*0.5;
  par[2] = (fgkExterInterModBorder1 - fgkInterCentrModBorder2)*0.5 - 2.;
  gMC->Gsvolu("FIF2", "BOX ", idtmed[502], par, 3); // G10
  //xcoor = 0.;
  //ycoor = kAlCoverThickness*0.5 + kInterfaceCardThickness*0.5;
  zcoor = (fgkExterInterModBorder1 + fgkInterCentrModBorder2)*0.5;
  gMC->Gspos("FIF2", 1, "FPEA", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FIF2", 2, "FPEA", xcoor, ycoor,-zcoor, 0, "ONLY");
  if (fTOFHoles) {
    gMC->Gspos("FIF2", 1, "FPEB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos("FIF2", 2, "FPEB", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  //par[0] = xtof*0.5;
  //par[1] = kInterfaceCardThickness*0.5;
  par[2] = (zlenA*0.5 + 2. - fgkExterInterModBorder1)*0.5 - 2.;
  gMC->Gsvolu("FIF3", "BOX ", idtmed[502], par, 3); // G10
  //xcoor = 0.;
  //ycoor = kAlCoverThickness*0.5 + kInterfaceCardThickness*0.5;
  zcoor = (zlenA*0.5 + 2. + fgkExterInterModBorder1)*0.5;
  gMC->Gspos("FIF3", 1, "FPEA", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FIF3", 2, "FPEA", xcoor, ycoor,-zcoor, 0, "ONLY");
  if (fTOFHoles) {
    gMC->Gspos("FIF3", 1, "FPEB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos("FIF3", 2, "FPEB", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  // volumes for flat cables
  // plastic
  const Float_t kPlasticFlatCableThickness = 0.25;
  par[0] = xtof*0.5;
  par[1] = kPlasticFlatCableThickness*0.5;
  par[2] = fgkInterCentrModBorder2 - 2.;
  gMC->Gsvolu("FFC1", "BOX ", idtmed[513], par, 3); // Plastic (CH2)
  //xcoor = 0.;
  ycoor = -kAlCoverThickness*0.5 - kPlasticFlatCableThickness*0.5;
  zcoor = 0.;
  gMC->Gspos("FFC1", 0, "FPEA", xcoor, ycoor, zcoor, 0, "ONLY");

  //par[0] = xtof*0.5;
  //par[1] = kPlasticFlatCableThickness*0.5;
  par[2] = (fgkExterInterModBorder1 - fgkInterCentrModBorder2)*0.5 - 2.;
  gMC->Gsvolu("FFC2", "BOX ", idtmed[513], par, 3); // Plastic (CH2)
  //xcoor = 0.;
  //ycoor = -kAlCoverThickness*0.5 - kPlasticFlatCableThickness*0.5;
  zcoor = (fgkExterInterModBorder1 + fgkInterCentrModBorder2)*0.5;
  gMC->Gspos("FFC2", 1, "FPEA", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FFC2", 2, "FPEA", xcoor, ycoor,-zcoor, 0, "ONLY");
  if (fTOFHoles) {
    gMC->Gspos("FFC2", 1, "FPEB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos("FFC2", 2, "FPEB", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  //par[0] = xtof*0.5;
  //par[1] = kPlasticFlatCableThickness*0.5;
  par[2] = (zlenA*0.5 + 2. - fgkExterInterModBorder1)*0.5 - 2.;
  gMC->Gsvolu("FFC3", "BOX ", idtmed[513], par, 3); // Plastic (CH2)
  //xcoor = 0.;
  //ycoor = -kAlCoverThickness*0.5 - kPlasticFlatCableThickness*0.5;
  zcoor = (zlenA*0.5 + 2. + fgkExterInterModBorder1)*0.5;
  gMC->Gspos("FFC3", 1, "FPEA", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FFC3", 2, "FPEA", xcoor, ycoor,-zcoor, 0, "ONLY");
  if (fTOFHoles) {
    gMC->Gspos("FFC3", 1, "FPEB", xcoor, ycoor, zcoor, 0, "ONLY");
    gMC->Gspos("FFC3", 2, "FPEB", xcoor, ycoor,-zcoor, 0, "ONLY");
  }

  // Cu
  const Float_t kCopperFlatCableThickness = 0.01;
  par[0] = xtof*0.5;
  par[1] = kCopperFlatCableThickness*0.5;
  par[2] = fgkInterCentrModBorder2 - 2.;
  gMC->Gsvolu("FCC1", "BOX ", idtmed[512], par, 3); // Cu
  gMC->Gspos("FCC1", 0, "FFC1", 0., 0., 0., 0, "ONLY");

  //par[0] = xtof*0.5;
  //par[1] = kCopperFlatCableThickness*0.5;
  par[2] = (fgkExterInterModBorder1 - fgkInterCentrModBorder2)*0.5 - 2.;
  gMC->Gsvolu("FCC2", "BOX ", idtmed[512], par, 3); // Cu
  gMC->Gspos("FCC2", 0, "FFC2", 0., 0., 0., 0, "ONLY");

  //par[0] = xtof*0.5;
  //par[1] = kCopperFlatCableThickness*0.5;
  par[2] = (zlenA*0.5 + 2. - fgkExterInterModBorder1)*0.5 - 2.;
  gMC->Gsvolu("FCC3", "BOX ", idtmed[512], par, 3); // Cu
  gMC->Gspos("FCC3", 0, "FFC3", 0., 0., 0., 0, "ONLY");

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeModulesInBTOFvolumes(Float_t ytof, Float_t zlenA) const
{
  //
  // Fill BTOF_%i (for i=0,...17) volumes
  // with volumes FTOA (MRPC strip container),
  // In case of TOF holes, three sectors (i.e. 13th, 14th and 15th)
  // are filled with volumes: FTOB and FTOC (MRPC containers),
  //

  const Int_t kSize=16;

  Int_t idrotm[1];

  //AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,-90.);
  AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,270.);

  Float_t xcoor, ycoor, zcoor;
  xcoor = 0.;

  // Positioning of fibre glass modules (FTOA, FTOB and FTOC)
  for(Int_t isec=0; isec<fTOFGeometry->NSectors(); isec++){
    if(fTOFSectors[isec]==-1)continue;

    char name[kSize];
    snprintf(name, kSize, "BTOF%d",isec);
    if (fTOFHoles && (isec==13 || isec==14 || isec==15)) {
      //xcoor = 0.;
      ycoor = (zlenA*0.5 + fgkInterCentrModBorder1)*0.5;
      zcoor = -ytof * 0.25;
      gMC->Gspos("FTOB", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
      gMC->Gspos("FTOC", 0, name, xcoor,-ycoor, zcoor, idrotm[0], "ONLY");
    }
    else {
      //xcoor = 0.;
      ycoor = 0.;
      zcoor = -ytof * 0.25;
      gMC->Gspos("FTOA", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
    }
  }

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeCoversInBTOFvolumes() const
{
  //
  // Fill BTOF_%i (for i=0,...17) volumes
  // with volumes FPEA (to separate strips from FEA cards)
  // In case of TOF holes, three sectors (i.e. 13th, 14th and 15th)
  // are filled with FPEB volumes
  // (to separate MRPC strips from FEA cards)
  //

  const Int_t kSize=16;

  Int_t idrotm[1];

  //AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,-90.);
  AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,270.);

  Float_t xcoor, ycoor, zcoor;
  xcoor = 0.;
  ycoor = 0.;
  zcoor = fgkModuleCoverThickness*0.5;

  char name[kSize];

  // Positioning of module covers (FPEA, FPEB)
  for(Int_t isec=0; isec<fTOFGeometry->NSectors(); isec++) {
    if(fTOFSectors[isec]==-1)continue;
    snprintf(name, kSize, "BTOF%d",isec);
    if (fTOFHoles && (isec==13 || isec==14 || isec==15))
      gMC->Gspos("FPEB", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
    else
      gMC->Gspos("FPEA", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  }

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeBackInBTOFvolumes(Float_t ytof) const
{
  //
  // Fill BTOF_%i (for i=0,...17) volumes with volumes called FAIA and
  // FAIC (FEA cards and services container).
  // In case of TOF holes, three sectors (i.e. 13th, 14th and 15th) are
  // filled with volumes FAIB (FEA cards and services container).
  //

  const Int_t kSize=16;

  Int_t idrotm[1];

  //AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,-90.);
  AliMatrix(idrotm[0], 90.,  0., 0., 0., 90.,270.);

  Float_t xcoor, ycoor, zcoor;
  xcoor = 0.;
  ycoor = 0.;
  zcoor = fgkModuleCoverThickness + (ytof*0.5 - fgkModuleCoverThickness)*0.5;

  char name[kSize];

  // Positioning of FEA cards and services containers (FAIA, FAIC and FAIB)
  for(Int_t isec=0; isec<fTOFGeometry->NSectors(); isec++) {
    if(fTOFSectors[isec]==-1)continue;
    snprintf(name, kSize, "BTOF%d",isec);
    if (fgkFEAwithMasks[isec])
      gMC->Gspos("FAIA", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
    else {
      if (fTOFHoles && (isec==13 || isec==14 || isec==15))
	gMC->Gspos("FAIB", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
      else
	gMC->Gspos("FAIC", 0, name, xcoor, ycoor, zcoor, idrotm[0], "ONLY");
    }
  }

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeStripsInModules(Float_t ytof, Float_t zlenA) const
{
  //
  // Define MRPC strip volume, called FSTR
  // Insert FSTR volume in FLTA/B/C volumes
  //

  Float_t yFLT  = ytof*0.5 - fgkModuleWallThickness;

  Int_t *idtmed = fIdtmed->GetArray()-499;

  ///////////////// Detector itself //////////////////////

  const Int_t    knx   = fTOFGeometry->NpadX();  // number of pads along x
  const Int_t    knz   = fTOFGeometry->NpadZ();  // number of pads along z
  const Float_t  kPadX = fTOFGeometry->XPad();   // pad length along x
  const Float_t  kPadZ = fTOFGeometry->ZPad();   // pad length along z

  // new description for strip volume -double stack strip-
  // -- all constants are expressed in cm
  // height of different layers
  const Float_t khhony   = 1.0;       // height of HONY Layer
  const Float_t khpcby   = 0.08;      // height of PCB Layer
  const Float_t khrgly   = 0.055;     // height of RED GLASS Layer

  const Float_t khfiliy  = 0.125;     // height of FISHLINE Layer
  const Float_t khglassy = 0.160*0.5; // semi-height of GLASS Layer
  const Float_t khglfy   = khfiliy+2.*khglassy; // height of GLASS Layer

  const Float_t khcpcby  = 0.16;      // height of PCB  Central Layer
  const Float_t kwhonz   = 8.1;       // z dimension of HONEY Layer
  const Float_t kwpcbz1  = 10.64;     // z dimension of PCB Lower Layer
  const Float_t kwpcbz2  = 11.6;      // z dimension of PCB Upper Layer
  const Float_t kwcpcbz  = 12.4;      // z dimension of PCB Central Layer

  const Float_t kwrglz   = 8.;        // z dimension of RED GLASS Layer
  const Float_t kwglfz   = 7.;        // z dimension of GLASS Layer
  const Float_t klsensmx = knx*kPadX; // length of Sensitive Layer
  const Float_t khsensmy = 0.0105;    // height of Sensitive Layer
  const Float_t kwsensmz = knz*kPadZ; // width of Sensitive Layer

  // height of the FSTR Volume (the strip volume)
  const Float_t khstripy = 2.*khhony+2.*khpcby+4.*khrgly+2.*khglfy+khcpcby;

  // width  of the FSTR Volume (the strip volume)
  const Float_t kwstripz = kwcpcbz;
  // length of the FSTR Volume (the strip volume)
  const Float_t klstripx = fTOFGeometry->StripLength();


  // FSTR volume definition-filling this volume with non sensitive Gas Mixture
  Float_t parfp[3]={klstripx*0.5, khstripy*0.5, kwstripz*0.5};
  gMC->Gsvolu("FSTR", "BOX", idtmed[506], parfp, 3); // Freon mix

  Float_t posfp[3]={0.,0.,0.};

  // NOMEX (HONEYCOMB) Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khhony*0.5;
  parfp[2] = kwhonz*0.5;
  gMC->Gsvolu("FHON", "BOX", idtmed[501], parfp, 3); // Nomex (Honeycomb)
  // positioning 2 NOMEX Layers on FSTR volume
  //posfp[0] = 0.;
  posfp[1] =-khstripy*0.5 + parfp[1];
  //posfp[2] = 0.;
  gMC->Gspos("FHON", 1, "FSTR", 0., posfp[1], 0., 0, "ONLY");
  gMC->Gspos("FHON", 2, "FSTR", 0.,-posfp[1], 0., 0, "ONLY");
  
  // Lower PCB Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khpcby*0.5;
  parfp[2] = kwpcbz1*0.5;
  gMC->Gsvolu("FPC1", "BOX", idtmed[502], parfp, 3); // G10

  // Upper PCB Layer definition
  //parfp[0] = klstripx*0.5;
  //parfp[1] = khpcby*0.5;
  parfp[2] = kwpcbz2*0.5;
  gMC->Gsvolu("FPC2", "BOX", idtmed[502], parfp, 3); // G10

  // positioning 2 external PCB Layers in FSTR volume
  //posfp[0] = 0.;
  posfp[1] =-khstripy*0.5+khhony+parfp[1];
  //posfp[2] = 0.;
  gMC->Gspos("FPC1", 1, "FSTR", 0.,-posfp[1], 0., 0, "ONLY");
  gMC->Gspos("FPC2", 1, "FSTR", 0., posfp[1], 0., 0, "ONLY");

  // Central PCB layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khcpcby*0.5;
  parfp[2] = kwcpcbz*0.5;
  gMC->Gsvolu("FPCB", "BOX", idtmed[502], parfp, 3); // G10
  gGeoManager->GetVolume("FPCB")->VisibleDaughters(kFALSE);
  // positioning the central PCB layer
  gMC->Gspos("FPCB", 1, "FSTR", 0., 0., 0., 0, "ONLY");

  // Sensitive volume definition
  Float_t parfs[3] = {klsensmx*0.5, khsensmy*0.5, kwsensmz*0.5};
  gMC->Gsvolu("FSEN", "BOX", idtmed[507], parfs, 3); // Cu sensitive
  // dividing FSEN along z in knz=2 and along x in knx=48
  gMC->Gsdvn("FSEZ", "FSEN", knz, 3);
  gMC->Gsdvn("FPAD", "FSEZ", knx, 1);
  // positioning sensitive layer inside FPCB
  gMC->Gspos("FSEN", 1, "FPCB", 0., 0., 0., 0, "ONLY");

  // RED GLASS Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khrgly*0.5;
  parfp[2] = kwrglz*0.5;
  gMC->Gsvolu("FRGL", "BOX", idtmed[508], parfp, 3); // red glass
  // positioning 4 RED GLASS Layers in FSTR volume
  //posfp[0] = 0.;
  posfp[1] = -khstripy*0.5+khhony+khpcby+parfp[1];
  //posfp[2] = 0.;
  gMC->Gspos("FRGL", 1, "FSTR", 0., posfp[1], 0., 0, "ONLY");
  gMC->Gspos("FRGL", 4, "FSTR", 0.,-posfp[1], 0., 0, "ONLY");
  //posfp[0] = 0.;
  posfp[1] = (khcpcby+khrgly)*0.5;
  //posfp[2] = 0.;
  gMC->Gspos("FRGL", 2, "FSTR", 0.,-posfp[1], 0., 0, "ONLY");
  gMC->Gspos("FRGL", 3, "FSTR", 0., posfp[1], 0., 0, "ONLY");

  // GLASS Layer definition
  //parfp[0] = klstripx*0.5;
  parfp[1] = khglassy;
  parfp[2] = kwglfz*0.5;
  gMC->Gsvolu("FGLF", "BOX", idtmed[508], parfp, 3); // glass
  // positioning 2 GLASS Layers in FSTR volume
  //posfp[0] = 0.;
  posfp[1] = (khcpcby + khglfy)*0.5 + khrgly;
  //posfp[2] = 0.;
  gMC->Gspos("FGLF", 1, "FSTR", 0.,-posfp[1], 0., 0, "ONLY");
  gMC->Gspos("FGLF", 2, "FSTR", 0., posfp[1], 0., 0, "ONLY");

  // Positioning the Strips (FSTR volumes) in the FLT volumes
  Int_t maxStripNumbers [5] ={fTOFGeometry->NStripC(),
			      fTOFGeometry->NStripB(),
			      fTOFGeometry->NStripA(),
			      fTOFGeometry->NStripB(),
			      fTOFGeometry->NStripC()};

  Int_t idrotm[91];

  Int_t totalStrip = 0;
  Float_t xpos, zpos, ypos, ang;
  for(Int_t iplate = 0; iplate < fTOFGeometry->NPlates(); iplate++){
    if (iplate>0) totalStrip += maxStripNumbers[iplate-1];
    for(Int_t istrip = 0; istrip < maxStripNumbers[iplate]; istrip++){

      ang = fTOFGeometry->GetAngles(iplate,istrip);
      AliDebug(1, Form(" iplate = %1i, istrip = %2i ---> ang = %f", iplate, istrip, ang));
 
      if (ang>0.)       AliMatrix (idrotm[istrip+totalStrip],90.,0.,90.+ang,90., ang, 90.);
      else if (ang==0.) AliMatrix (idrotm[istrip+totalStrip],90.,0.,90.,90., 0., 0.);
      else if (ang<0.)  AliMatrix (idrotm[istrip+totalStrip],90.,0.,90.+ang,90.,-ang,270.);

      xpos = 0.;
      ypos = fTOFGeometry->GetHeights(iplate,istrip) + yFLT*0.5;
      zpos = fTOFGeometry->GetDistances(iplate,istrip);
      gMC->Gspos("FSTR", istrip+totalStrip+1, "FLTA", xpos, ypos,-zpos, idrotm[istrip+totalStrip], "ONLY");

      if (fTOFHoles) {
	if (istrip+totalStrip+1>53)
	  gMC->Gspos("FSTR", istrip+totalStrip+1, "FLTC", xpos, ypos,-zpos-(zlenA*0.5 - 2.*fgkModuleWallThickness + fgkInterCentrModBorder1)*0.5, idrotm[istrip+totalStrip], "ONLY");
	if (istrip+totalStrip+1<39)
	  gMC->Gspos("FSTR", istrip+totalStrip+1, "FLTB", xpos, ypos,-zpos+(zlenA*0.5 - 2.*fgkModuleWallThickness + fgkInterCentrModBorder1)*0.5, idrotm[istrip+totalStrip], "ONLY");
      }
    }
  }

}

//_____________________________________________________________________________
void AliTOFv6T0::CreateBackZone(Float_t xtof, Float_t ytof, Float_t zlenA) const
{
  //
  // Define:
  //        - containers for FEA cards, cooling system
  //          signal cables and supermodule support structure
  //          (volumes called FAIA/B/C),
  //        - containers for FEA cards and some cooling
  //          elements for a FEA (volumes called FCA1/2).
  //

  Int_t *idtmed = fIdtmed->GetArray()-499;

  Int_t idrotm[1];

  // Definition of the air card containers (FAIA, FAIC and FAIB)

  Float_t  par[3];
  par[0] = xtof*0.5;
  par[1] = (ytof*0.5 - fgkModuleCoverThickness)*0.5;
  par[2] = zlenA*0.5;
  gMC->Gsvolu("FAIA", "BOX ", idtmed[500], par, 3); // Air
  if (fTOFHoles) gMC->Gsvolu("FAIB", "BOX ", idtmed[500], par, 3); // Air
  gMC->Gsvolu("FAIC", "BOX ", idtmed[500], par, 3); // Air

  Float_t feaParam[3] = {fgkFEAparameters[0], fgkFEAparameters[1], fgkFEAparameters[2]};
  Float_t feaRoof1[3] = {fgkRoof1parameters[0], fgkRoof1parameters[1], fgkRoof1parameters[2]};
  Float_t al3[3] = {fgkAl3parameters[0], fgkAl3parameters[1], fgkAl3parameters[2]};
  //Float_t feaRoof2[3] = {fgkRoof2parameters[0], fgkRoof2parameters[1], fgkRoof2parameters[2]};

  // FEA card mother-volume definition
  Float_t carpar[3] = {xtof*0.5 - fgkCBLw - fgkSawThickness,
		       feaParam[1] + feaRoof1[1] + fgkRoof2parameters[1]*0.5,
		       feaRoof1[2] + fgkBetweenLandMask*0.5 + al3[2]};
  gMC->Gsvolu("FCA1", "BOX ", idtmed[500], carpar, 3); // Air
  gMC->Gsvolu("FCA2", "BOX ", idtmed[500], carpar, 3); // Air

  // rotation matrix
  AliMatrix(idrotm[0],  90.,180., 90., 90.,180., 0.);

  // FEA card mother-volume positioning
  Float_t rowstep = 6.66;
  Float_t rowgap[5] = {13.5, 22.9, 16.94, 23.8, 20.4};
  Int_t rowb[5] = {6, 7, 6, 19, 7};
  Float_t carpos[3] = {0.,
		       -(ytof*0.5 - fgkModuleCoverThickness)*0.5 + carpar[1],
		       -0.8};
  gMC->Gspos("FCA1", 91, "FAIA", carpos[0], carpos[1], carpos[2], 0, "MANY");
  gMC->Gspos("FCA2", 91, "FAIC", carpos[0], carpos[1], carpos[2], 0, "MANY");

  Int_t row = 1;
  Int_t nrow = 0;
  for (Int_t sg= -1; sg< 2; sg+= 2) {
    carpos[2] = sg*zlenA*0.5 - 0.8;
    for (Int_t nb=0; nb<5; ++nb) {
      carpos[2] = carpos[2] - sg*(rowgap[nb] - rowstep);
      nrow = row + rowb[nb];
      for ( ; row < nrow ; ++row) {

        carpos[2] -= sg*rowstep;

	if (nb==4) {
	  gMC->Gspos("FCA1", row, "FAIA", carpos[0], carpos[1], carpos[2], 0, "ONLY");
	  gMC->Gspos("FCA2", row, "FAIC", carpos[0], carpos[1], carpos[2], 0, "ONLY");

	}
	else {
	  switch (sg) {
	  case 1:
	    gMC->Gspos("FCA1", row, "FAIA", carpos[0], carpos[1], carpos[2], 0, "ONLY");
	    gMC->Gspos("FCA2", row, "FAIC", carpos[0], carpos[1], carpos[2], 0, "ONLY");
	    break;
	  case -1:
	    gMC->Gspos("FCA1", row, "FAIA", carpos[0], carpos[1], carpos[2], idrotm[0], "ONLY");
	    gMC->Gspos("FCA2", row, "FAIC", carpos[0], carpos[1], carpos[2], idrotm[0], "ONLY");
	    break;
	  }

	}

      }
    }
  }

  if (fTOFHoles) {
    row = 1;
    for (Int_t sg= -1; sg< 2; sg+= 2) {
      carpos[2] = sg*zlenA*0.5 - 0.8;
      for (Int_t nb=0; nb<4; ++nb) {
        carpos[2] = carpos[2] - sg*(rowgap[nb] - rowstep);
        nrow = row + rowb[nb];
        for ( ; row < nrow ; ++row) {
          carpos[2] -= sg*rowstep;

	  switch (sg) {
	  case 1:
	    gMC->Gspos("FCA1", row, "FAIB", carpos[0], carpos[1], carpos[2], 0, "ONLY");
	    break;
	  case -1:
	    gMC->Gspos("FCA1", row, "FAIB", carpos[0], carpos[1], carpos[2], idrotm[0], "ONLY");
	    break;
	  }
	}
      }
    }
  }

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeFrontEndElectronics(Float_t xtof) const
{
  //
  // Fill FCA1/2 volumes with FEA cards (FFEA volumes).
  //

  Int_t *idtmed = fIdtmed->GetArray()-499;

  // FEA card volume definition
  Float_t feaParam[3] = {fgkFEAparameters[0], fgkFEAparameters[1], fgkFEAparameters[2]};
  gMC->Gsvolu("FFEA", "BOX ", idtmed[502], feaParam, 3); // G10

  Float_t al1[3] = {fgkAl1parameters[0], fgkAl1parameters[1], fgkAl1parameters[2]};
  Float_t al3[3] = {fgkAl3parameters[0], fgkAl3parameters[1], fgkAl3parameters[2]};
  Float_t feaRoof1[3] = {fgkRoof1parameters[0], fgkRoof1parameters[1], fgkRoof1parameters[2]};
  //Float_t feaRoof2[3] = {fgkRoof2parameters[0], fgkRoof2parameters[1], fgkRoof2parameters[2]};

  Float_t carpar[3] = {xtof*0.5 - fgkCBLw - fgkSawThickness,
		       feaParam[1] + feaRoof1[1] + fgkRoof2parameters[1]*0.5,
		       feaRoof1[2] + fgkBetweenLandMask*0.5 + al3[2]};

  // FEA card volume positioning
  Float_t xCoor = xtof*0.5 - 25.;
  Float_t yCoor =-carpar[1] + feaParam[1];
  Float_t zCoor =-carpar[2] + (2.*feaRoof1[2] - 2.*al1[2] - feaParam[2]);
  gMC->Gspos("FFEA", 1, "FCA1",-xCoor, yCoor, zCoor, 0, "ONLY");
  gMC->Gspos("FFEA", 4, "FCA1", xCoor, yCoor, zCoor, 0, "ONLY");
  gMC->Gspos("FFEA", 1, "FCA2",-xCoor, yCoor, zCoor, 0, "ONLY");
  gMC->Gspos("FFEA", 4, "FCA2", xCoor, yCoor, zCoor, 0, "ONLY");
  xCoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FFEA", 2, "FCA1",-xCoor, yCoor, zCoor, 0, "ONLY");
  gMC->Gspos("FFEA", 3, "FCA1", xCoor, yCoor, zCoor, 0, "ONLY");
  gMC->Gspos("FFEA", 2, "FCA2",-xCoor, yCoor, zCoor, 0, "ONLY");
  gMC->Gspos("FFEA", 3, "FCA2", xCoor, yCoor, zCoor, 0, "ONLY");

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeFEACooling(Float_t xtof) const
{
  //
  // Make cooling system attached to each FEA card
  // (FAL1, FRO1 and FBAR/1/2 volumes)
  // in FCA1/2 volume containers.
  //

  Int_t *idtmed = fIdtmed->GetArray()-499;

  // first FEA cooling element definition
  Float_t al1[3] = {fgkAl1parameters[0], fgkAl1parameters[1], fgkAl1parameters[2]};
  gMC->Gsvolu("FAL1", "BOX ", idtmed[504], al1, 3); // Al

  // second FEA cooling element definition
  Float_t feaRoof1[3] = {fgkRoof1parameters[0], fgkRoof1parameters[1], fgkRoof1parameters[2]};
  gMC->Gsvolu("FRO1", "BOX ", idtmed[504], feaRoof1, 3); // Al

  Float_t al3[3] = {fgkAl3parameters[0], fgkAl3parameters[1], fgkAl3parameters[2]};
  //Float_t feaRoof2[3] = {fgkRoof2parameters[0], fgkRoof2parameters[1], fgkRoof2parameters[2]};

  // definition and positioning of a small air groove in the FRO1 volume
  Float_t airHole[3] = {fgkRoof2parameters[0], fgkRoof2parameters[1]*0.5, feaRoof1[2]};
  gMC->Gsvolu("FREE", "BOX ", idtmed[500], airHole, 3); // Air
  gMC->Gspos("FREE", 1, "FRO1", 0., feaRoof1[1]-airHole[1], 0., 0, "ONLY");
  gGeoManager->GetVolume("FRO1")->VisibleDaughters(kFALSE);

  // third FEA cooling element definition
  Float_t bar[3] = {fgkBar[0], fgkBar[1], fgkBar[2]};
  gMC->Gsvolu("FBAR", "BOX ", idtmed[504], bar, 3); // Al

  Float_t feaParam[3] = {fgkFEAparameters[0], fgkFEAparameters[1], fgkFEAparameters[2]};

  Float_t carpar[3] = {xtof*0.5 - fgkCBLw - fgkSawThickness,
		       feaParam[1] + feaRoof1[1] + fgkRoof2parameters[1]*0.5,
		       feaRoof1[2] + fgkBetweenLandMask*0.5 + al3[2]};

  // fourth FEA cooling element definition
  Float_t bar1[3] = {fgkBar1[0], fgkBar1[1], fgkBar1[2]};
  gMC->Gsvolu("FBA1", "BOX ", idtmed[504], bar1, 3); // Al

  // fifth FEA cooling element definition
  Float_t bar2[3] = {fgkBar2[0], fgkBar2[1], fgkBar2[2]};
  gMC->Gsvolu("FBA2", "BOX ", idtmed[504], bar2, 3); // Al

  // first FEA cooling element positioning
  Float_t xcoor = xtof*0.5 - 25.;
  Float_t ycoor = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - al1[1];
  Float_t zcoor =-carpar[2] + 2.*feaRoof1[2] - al1[2];
  gMC->Gspos("FAL1", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL1", 4, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL1", 1, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL1", 4, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FAL1", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL1", 3, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL1", 2, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL1", 3, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");

  // second FEA cooling element positioning
  xcoor = xtof*0.5 - 25.;
  ycoor = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - feaRoof1[1];
  zcoor =-carpar[2] + feaRoof1[2];
  gMC->Gspos("FRO1", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "MANY"); // (AdC)
  gMC->Gspos("FRO1", 4, "FCA1", xcoor, ycoor, zcoor, 0, "MANY"); // (AdC)
  gMC->Gspos("FRO1", 1, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FRO1", 4, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FRO1", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "MANY"); // (AdC)
  gMC->Gspos("FRO1", 3, "FCA1", xcoor, ycoor, zcoor, 0, "MANY"); // (AdC)
  gMC->Gspos("FRO1", 2, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FRO1", 3, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");

  // third FEA cooling element positioning
  xcoor = xtof*0.5 - 25.;
  ycoor = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - bar[1];
  zcoor =-carpar[2] + bar[2];
  gMC->Gspos("FBAR", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAR", 4, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAR", 1, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAR", 4, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FBAR", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAR", 3, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAR", 2, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAR", 3, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");

  // fourth FEA cooling element positioning
  Float_t tubepar[3] = {0., 0.4, xtof*0.5 - fgkCBLw};
  xcoor = xtof*0.5 - 25.;
  ycoor = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - bar[1];
  zcoor =-carpar[2] + 2.*bar[2] + 2.*tubepar[1] + bar1[2];
  gMC->Gspos("FBA1", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA1", 4, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA1", 1, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA1", 4, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FBA1", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA1", 3, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA1", 2, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA1", 3, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");

  // fifth FEA cooling element positioning
  xcoor = xtof*0.5 - 25.;
  ycoor = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - bar2[1];
  zcoor =-carpar[2] + 2.*bar[2] + bar2[2];
  gMC->Gspos("FBA2", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 4, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 1, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 4, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FBA2", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 3, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 2, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 3, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");

  xcoor = xtof*0.5 - 25.;
  ycoor = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - 2.*bar2[1] - 2.*tubepar[1] - bar2[1];
  zcoor =-carpar[2] + 2.*bar[2] + bar2[2];
  gMC->Gspos("FBA2", 5, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 8, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 5, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 8, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FBA2", 6, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 7, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 6, "FCA2",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBA2", 7, "FCA2", xcoor, ycoor, zcoor, 0, "ONLY");

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeNinoMask(Float_t xtof) const
{
  //
  // Make cooling Nino mask
  // for each FEA card (FAL2/3 and FRO2 volumes)
  // in FCA1 volume container.
  //

  Int_t *idtmed = fIdtmed->GetArray()-499;

  // first Nino ASIC mask volume definition
  Float_t al2[3] = {fgkAl2parameters[0], fgkAl2parameters[1], fgkAl2parameters[2]};
  gMC->Gsvolu("FAL2", "BOX ", idtmed[504], al2, 3); // Al

  // second Nino ASIC mask volume definition
  Float_t al3[3] = {fgkAl3parameters[0], fgkAl3parameters[1], fgkAl3parameters[2]};
  gMC->Gsvolu("FAL3", "BOX ", idtmed[504], al3, 3); // Al

  // third Nino ASIC mask volume definition
  Float_t feaRoof2[3] = {fgkRoof2parameters[0], fgkRoof2parameters[1], fgkRoof2parameters[2]};
  gMC->Gsvolu("FRO2", "BOX ", idtmed[504], feaRoof2, 3); // Al

  Float_t feaRoof1[3] = {fgkRoof1parameters[0], fgkRoof1parameters[1], fgkRoof1parameters[2]};
  Float_t feaParam[3] = {fgkFEAparameters[0], fgkFEAparameters[1], fgkFEAparameters[2]};

  Float_t carpar[3] = {xtof*0.5 - fgkCBLw - fgkSawThickness,
		       feaParam[1] + feaRoof1[1] + fgkRoof2parameters[1]*0.5,
		       feaRoof1[2] + fgkBetweenLandMask*0.5 + al3[2]};

  // first Nino ASIC mask volume positioning
  Float_t xcoor = xtof*0.5 - 25.;
  Float_t ycoor = carpar[1] - 2.*al3[1];
  Float_t zcoor = carpar[2] - 2.*al3[2] - al2[2];
  gMC->Gspos("FAL2", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL2", 4, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FAL2", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL2", 3, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");

  // second Nino ASIC mask volume positioning
  xcoor = xtof*0.5 - 25.;
  ycoor = carpar[1] - al3[1];
  zcoor = carpar[2] - al3[2];
  gMC->Gspos("FAL3", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL3", 4, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FAL3", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FAL3", 3, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");

  // third Nino ASIC mask volume positioning
  xcoor = xtof*0.5 - 25.;
  ycoor = carpar[1] - fgkRoof2parameters[1];
  zcoor = carpar[2] - 2.*al3[2] - fgkRoof2parameters[2];
  gMC->Gspos("FRO2", 1, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FRO2", 4, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");
  xcoor = feaParam[0] + (fgkFEAwidth2*0.5 - fgkFEAwidth1);
  gMC->Gspos("FRO2", 2, "FCA1",-xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FRO2", 3, "FCA1", xcoor, ycoor, zcoor, 0, "ONLY");

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeSuperModuleCooling(Float_t xtof, Float_t ytof, Float_t zlenA) const
{
  //
  // Make cooling tubes (FTUB volume)
  // and cooling bars (FTLN and FLO1/2/3 volumes)
  // in FAIA/B/C volume containers.
  //

  Int_t *idtmed = fIdtmed->GetArray()-499;

  Int_t idrotm[1];

  // cooling tube volume definition
  Float_t tubepar[3] = {0., 0.4, xtof*0.5 - fgkCBLw - fgkSawThickness};
  gMC->Gsvolu("FTUB", "TUBE", idtmed[512], tubepar, 3); // Cu

  // water cooling tube volume definition
  Float_t tubeparW[3] = {0., 0.3, tubepar[2]};
  gMC->Gsvolu("FITU", "TUBE", idtmed[509], tubeparW, 3); // H2O

  // Positioning of the water tube into the steel one
  gMC->Gspos("FITU", 1, "FTUB", 0., 0., 0., 0, "ONLY");

  // definition of transverse components of SM cooling system
  Float_t trapar[3] = {tubepar[2], 6.175/*6.15*/, 0.7};
  gMC->Gsvolu("FTLN", "BOX ", idtmed[504], trapar, 3); // Al

  // rotation matrix
  AliMatrix(idrotm[0], 180., 90., 90., 90., 90., 0.);

  Float_t feaParam[3] = {fgkFEAparameters[0], fgkFEAparameters[1], fgkFEAparameters[2]};
  Float_t feaRoof1[3] = {fgkRoof1parameters[0], fgkRoof1parameters[1], fgkRoof1parameters[2]};
  Float_t bar[3] = {fgkBar[0], fgkBar[1], fgkBar[2]};
  Float_t bar2[3] = {fgkBar2[0], fgkBar2[1], fgkBar2[2]};
  Float_t al3[3] = {fgkAl3parameters[0], fgkAl3parameters[1], fgkAl3parameters[2]};
  //Float_t feaRoof2[3] = {fgkRoof2parameters[0], fgkRoof2parameters[1], fgkRoof2parameters[2]};

  Float_t carpar[3] = {xtof*0.5 - fgkCBLw - fgkSawThickness,
		       feaParam[1] + feaRoof1[1] + fgkRoof2parameters[1]*0.5,
		       feaRoof1[2] + fgkBetweenLandMask*0.5 + al3[2]};

  Float_t ytub =-(ytof*0.5 - fgkModuleCoverThickness)*0.5 + carpar[1] +
    carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - 2.*bar2[1] - tubepar[1];

  // Positioning of tubes for the SM cooling system
  Float_t ycoor = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - 2.*bar2[1] - tubepar[1];
  Float_t zcoor =-carpar[2] + 2.*bar[2] + tubepar[1];
  gMC->Gspos("FTUB", 1, "FCA1", 0., ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FTUB", 1, "FCA2", 0., ycoor, zcoor, idrotm[0], "ONLY");
  gGeoManager->GetVolume("FTUB")->VisibleDaughters(kFALSE);

  Float_t yFLTN = trapar[1] - (ytof*0.5 - fgkModuleCoverThickness)*0.5;
  for (Int_t sg= -1; sg< 2; sg+= 2) {
    // Positioning of transverse components for the SM cooling system
    gMC->Gspos("FTLN", 5+4*sg, "FAIA", 0., yFLTN, 369.9*sg, 0, "MANY");
    gMC->Gspos("FTLN", 5+3*sg, "FAIA", 0., yFLTN, 366.9*sg, 0, "MANY");
    gMC->Gspos("FTLN", 5+2*sg, "FAIA", 0., yFLTN, 198.8*sg, 0, "MANY");
    gMC->Gspos("FTLN",   5+sg, "FAIA", 0., yFLTN, 56.82*sg, 0, "MANY");
    gMC->Gspos("FTLN", 5+4*sg, "FAIC", 0., yFLTN, 369.9*sg, 0, "MANY");
    gMC->Gspos("FTLN", 5+3*sg, "FAIC", 0., yFLTN, 366.9*sg, 0, "MANY");
    gMC->Gspos("FTLN", 5+2*sg, "FAIC", 0., yFLTN, 198.8*sg, 0, "MANY");
    gMC->Gspos("FTLN",   5+sg, "FAIC", 0., yFLTN, 56.82*sg, 0, "MANY");
  }

  // definition of longitudinal components of SM cooling system
  Float_t lonpar1[3] = {2., 0.5, 56.82 - trapar[2]};
  Float_t lonpar2[3] = {lonpar1[0], lonpar1[1], (198.8 - 56.82)*0.5 - trapar[2]};
  Float_t lonpar3[3] = {lonpar1[0], lonpar1[1], (366.9 - 198.8)*0.5 - trapar[2]};
  gMC->Gsvolu("FLO1", "BOX ", idtmed[504], lonpar1, 3); // Al
  gMC->Gsvolu("FLO2", "BOX ", idtmed[504], lonpar2, 3); // Al
  gMC->Gsvolu("FLO3", "BOX ", idtmed[504], lonpar3, 3); // Al

  // Positioning of longitudinal components for the SM cooling system
  ycoor =  ytub + (tubepar[1] + 2.*bar2[1] + lonpar1[1]);
  gMC->Gspos("FLO1",  4, "FAIA",-24., ycoor, 0., 0, "MANY");
  gMC->Gspos("FLO1",  2, "FAIA", 24., ycoor, 0., 0, "MANY");
  gMC->Gspos("FLO1",  4, "FAIC",-24., ycoor, 0., 0, "MANY");
  gMC->Gspos("FLO1",  2, "FAIC", 24., ycoor, 0., 0, "MANY");

  zcoor = (198.8 + 56.82)*0.5;
  gMC->Gspos("FLO2",  4, "FAIA",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  2, "FAIA", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  4, "FAIC",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  2, "FAIC", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  8, "FAIA",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  6, "FAIA", 24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  8, "FAIC",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  6, "FAIC", 24., ycoor, zcoor, 0, "MANY");

  zcoor = (366.9 + 198.8)*0.5;
  gMC->Gspos("FLO3",  4, "FAIA",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  2, "FAIA", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  4, "FAIC",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  2, "FAIC", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  8, "FAIA",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  6, "FAIA", 24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  8, "FAIC",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  6, "FAIC", 24., ycoor, zcoor, 0, "MANY");

  ycoor =  ytub - (tubepar[1] + 2.*bar2[1] + lonpar1[1]);
  gMC->Gspos("FLO1",  3, "FAIA",-24., ycoor, 0., 0, "MANY");
  gMC->Gspos("FLO1",  1, "FAIA", 24., ycoor, 0., 0, "MANY");
  gMC->Gspos("FLO1",  3, "FAIC",-24., ycoor, 0., 0, "MANY");
  gMC->Gspos("FLO1",  1, "FAIC", 24., ycoor, 0., 0, "MANY");

  zcoor = (198.8 + 56.82)*0.5;
  gMC->Gspos("FLO2",  3, "FAIA",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  1, "FAIA", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  3, "FAIC",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  1, "FAIC", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  7, "FAIA",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  5, "FAIA", 24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  7, "FAIC",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO2",  5, "FAIC", 24., ycoor, zcoor, 0, "MANY");

  zcoor = (366.9 + 198.8)*0.5;
  gMC->Gspos("FLO3",  3, "FAIA",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  1, "FAIA", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  3, "FAIC",-24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  1, "FAIC", 24., ycoor,-zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  7, "FAIA",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  5, "FAIA", 24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  7, "FAIC",-24., ycoor, zcoor, 0, "MANY");
  gMC->Gspos("FLO3",  5, "FAIC", 24., ycoor, zcoor, 0, "MANY");


  Float_t carpos[3] = {25. - xtof*0.5,
		       (11.5 - (ytof*0.5 - fgkModuleCoverThickness))*0.5,
		       0.};
  if (fTOFHoles) {
    for (Int_t sg= -1; sg< 2; sg+= 2) {
      carpos[2] = sg*zlenA*0.5;
      gMC->Gspos("FTLN", 5+4*sg, "FAIB", 0., yFLTN, 369.9*sg, 0, "MANY");
      gMC->Gspos("FTLN", 5+3*sg, "FAIB", 0., yFLTN, 366.9*sg, 0, "MANY");
      gMC->Gspos("FTLN", 5+2*sg, "FAIB", 0., yFLTN, 198.8*sg, 0, "MANY");
      gMC->Gspos("FTLN",   5+sg, "FAIB", 0., yFLTN, 56.82*sg, 0, "MANY");
    }

    ycoor =  ytub + (tubepar[1] + 2.*bar2[1] + lonpar1[1]);
    zcoor = (198.8 + 56.82)*0.5;
    gMC->Gspos("FLO2", 2, "FAIB",-24., ycoor,-zcoor, 0, "MANY");
    gMC->Gspos("FLO2", 1, "FAIB",-24., ycoor, zcoor, 0, "MANY");
    zcoor = (366.9 + 198.8)*0.5;
    gMC->Gspos("FLO3", 2, "FAIB",-24., ycoor,-zcoor, 0, "MANY");
    gMC->Gspos("FLO3", 1, "FAIB",-24., ycoor, zcoor, 0, "MANY");
    ycoor =  ytub - (tubepar[1] + 2.*bar2[1] + lonpar1[1]);
    zcoor = (198.8 + 56.82)*0.5;
    gMC->Gspos("FLO2", 4, "FAIB", 24., ycoor,-zcoor, 0, "MANY");
    gMC->Gspos("FLO2", 3, "FAIB", 24., ycoor, zcoor, 0, "MANY");
    zcoor = (366.9 + 198.8)*0.5;
    gMC->Gspos("FLO3", 4, "FAIB", 24., ycoor,-zcoor, 0, "MANY");
    gMC->Gspos("FLO3", 3, "FAIB", 24., ycoor, zcoor, 0, "MANY");

  }

  Float_t barS[3] = {fgkBarS[0], fgkBarS[1], fgkBarS[2]};
  gMC->Gsvolu("FBAS", "BOX ", idtmed[504], barS, 3); // Al

  Float_t barS1[3] = {fgkBarS1[0], fgkBarS1[1], fgkBarS1[2]};
  gMC->Gsvolu("FBS1", "BOX ", idtmed[504], barS1, 3); // Al

  Float_t barS2[3] = {fgkBarS2[0], fgkBarS2[1], fgkBarS2[2]};
  gMC->Gsvolu("FBS2", "BOX ", idtmed[504], barS2, 3); // Al

  Float_t ytubBis = carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - 2.*barS2[1] - tubepar[1];
  ycoor = ytubBis;
  zcoor =-carpar[2] + barS[2];
  gMC->Gspos("FBAS", 1, "FCA1",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAS", 2, "FCA1", 24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAS", 1, "FCA2",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBAS", 2, "FCA2", 24., ycoor, zcoor, 0, "ONLY");

  zcoor =-carpar[2] + 2.*barS[2] + 2.*tubepar[1] + barS1[2];
  gMC->Gspos("FBS1", 1, "FCA1",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS1", 2, "FCA1", 24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS1", 1, "FCA2",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS1", 2, "FCA2", 24., ycoor, zcoor, 0, "ONLY");

  ycoor = ytubBis + (tubepar[1] + barS2[1]);
  zcoor =-carpar[2] + 2.*barS[2] + barS2[2];
  gMC->Gspos("FBS2", 1, "FCA1",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS2", 2, "FCA1", 24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS2", 1, "FCA2",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS2", 2, "FCA2", 24., ycoor, zcoor, 0, "ONLY");

  ycoor = ytubBis - (tubepar[1] + barS2[1]);
  //zcoor =-carpar[2] + 2.*barS[2] + barS2[2];
  gMC->Gspos("FBS2", 3, "FCA1",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS2", 4, "FCA1", 24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS2", 3, "FCA2",-24., ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FBS2", 4, "FCA2", 24., ycoor, zcoor, 0, "ONLY");

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeSuperModuleServices(Float_t xtof, Float_t ytof, Float_t zlenA) const
{
  //
  // Make signal cables (FCAB/L and FCBL/B volumes),
  // supemodule cover (FCOV volume) and wall (FSAW volume)
  // in FAIA/B/C volume containers.
  //

  Int_t *idtmed = fIdtmed->GetArray()-499;

  Int_t idrotm[3];

  Float_t tubepar[3] = {0., 0.4, xtof*0.5 - fgkCBLw - fgkSawThickness};
  Float_t al1[3] = {fgkAl1parameters[0], fgkAl1parameters[1], fgkAl1parameters[2]};
  Float_t al3[3] = {fgkAl3parameters[0], fgkAl3parameters[1], fgkAl3parameters[2]};
  Float_t feaRoof1[3] = {fgkRoof1parameters[0], fgkRoof1parameters[1], fgkRoof1parameters[2]};
  //Float_t feaRoof2[3] = {fgkRoof2parameters[0], fgkRoof2parameters[1], fgkRoof2parameters[2]};
  Float_t feaParam[3] = {fgkFEAparameters[0], fgkFEAparameters[1], fgkFEAparameters[2]};

  // FEA cables definition
  Float_t cbpar[3] = {0., 0.5, (tubepar[2] - (fgkFEAwidth2 - fgkFEAwidth1/6.)*0.5)*0.5};
  gMC->Gsvolu("FCAB", "TUBE", idtmed[510], cbpar, 3);    // copper+alu

  Float_t cbparS[3] = {cbpar[0], cbpar[1], (tubepar[2] - (xtof*0.5 - 25. + (fgkFEAwidth1 - fgkFEAwidth1/6.)*0.5))*0.5};
  gMC->Gsvolu("FCAL", "TUBE", idtmed[510], cbparS, 3);    // copper+alu

  // rotation matrix
  AliMatrix(idrotm[0], 180., 90., 90., 90., 90., 0.);

  Float_t carpar[3] = {xtof*0.5 - fgkCBLw - fgkSawThickness,
		       feaParam[1] + feaRoof1[1] + fgkRoof2parameters[1]*0.5,
		       feaRoof1[2] + fgkBetweenLandMask*0.5 + al3[2]};

  Float_t bar2[3] = {fgkBar2[0], fgkBar2[1], fgkBar2[2]};
  Float_t ytub =-(ytof*0.5 - fgkModuleCoverThickness)*0.5 + carpar[1] +
    carpar[1] - 2.*fgkRoof2parameters[1]*0.5 - 2.*feaRoof1[1] - 2.*bar2[1] - tubepar[1];

  // FEA cables positioning
  Float_t xcoor = (tubepar[2] + (fgkFEAwidth2 - fgkFEAwidth1/6.)*0.5)*0.5;
  Float_t ycoor = ytub - 3.;
  Float_t zcoor =-carpar[2] + (2.*feaRoof1[2] - 2.*al1[2] - 2.*feaParam[2] - cbpar[1]);
  gMC->Gspos("FCAB", 1, "FCA1",-xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FCAB", 2, "FCA1", xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FCAB", 1, "FCA2",-xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FCAB", 2, "FCA2", xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  xcoor = (tubepar[2] + (xtof*0.5 - 25. + (fgkFEAwidth1 - fgkFEAwidth1/6.)*0.5))*0.5;
  ycoor -= 2.*cbpar[1];
  gMC->Gspos("FCAL", 1, "FCA1",-xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FCAL", 2, "FCA1", xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FCAL", 1, "FCA2",-xcoor, ycoor, zcoor, idrotm[0], "ONLY");
  gMC->Gspos("FCAL", 2, "FCA2", xcoor, ycoor, zcoor, idrotm[0], "ONLY");


  // Cables and tubes on the side blocks
  // constants definition
  const Float_t kCBLl   = zlenA*0.5; // length of block
  const Float_t kCBLlh  = zlenA*0.5 - fgkInterCentrModBorder2; // length  of block in case of holes
  //const Float_t fgkCBLw   = 13.5;      // width of block
  //const Float_t fgkCBLh1  = 2.;        // min. height of block
  //const Float_t fgkCBLh2  = 12.3;      // max. height of block
  //const Float_t fgkSawThickness = 1.; // Al wall thickness

  // lateral cable and tube volume definition
  Float_t tgal =  (fgkCBLh2 - fgkCBLh1)/(2.*kCBLl);
  Float_t cblpar[11];
  cblpar[0] = fgkCBLw *0.5;
  cblpar[1] = 0.;
  cblpar[2] = 0.;
  cblpar[3] = kCBLl *0.5;
  cblpar[4] = fgkCBLh1 *0.5;
  cblpar[5] = fgkCBLh2 *0.5;
  cblpar[6] = TMath::ATan(tgal)*kRaddeg;
  cblpar[7] = kCBLl *0.5;
  cblpar[8] = fgkCBLh1 *0.5;
  cblpar[9] = fgkCBLh2 *0.5;
  cblpar[10]= cblpar[6];
  gMC->Gsvolu("FCBL", "TRAP", idtmed[511], cblpar, 11); // cables and tubes mix 

  // Side Al Walls definition
  Float_t sawpar[3] = {fgkSawThickness*0.5, fgkCBLh2*0.5, kCBLl};
  gMC->Gsvolu("FSAW", "BOX ", idtmed[504], sawpar,  3); // Al

  AliMatrix(idrotm[1], 90., 90., 180., 0., 90., 180.);
  AliMatrix(idrotm[2], 90., 90., 0., 0., 90., 0.);

  // lateral cable and tube volume positioning
  xcoor = (xtof - fgkCBLw)*0.5 - 2.*sawpar[0];
  ycoor = (fgkCBLh1 + fgkCBLh2)*0.25 - (ytof*0.5 - fgkModuleCoverThickness)*0.5;
  zcoor = kCBLl*0.5;
  gMC->Gspos("FCBL", 1, "FAIA", -xcoor, ycoor, -zcoor, idrotm[1], "ONLY");
  gMC->Gspos("FCBL", 2, "FAIA",  xcoor, ycoor, -zcoor, idrotm[1], "ONLY");
  gMC->Gspos("FCBL", 3, "FAIA", -xcoor, ycoor,  zcoor, idrotm[2], "ONLY");
  gMC->Gspos("FCBL", 4, "FAIA",  xcoor, ycoor,  zcoor, idrotm[2], "ONLY");
  gMC->Gspos("FCBL", 1, "FAIC", -xcoor, ycoor, -zcoor, idrotm[1], "ONLY");
  gMC->Gspos("FCBL", 2, "FAIC",  xcoor, ycoor, -zcoor, idrotm[1], "ONLY");
  gMC->Gspos("FCBL", 3, "FAIC", -xcoor, ycoor,  zcoor, idrotm[2], "ONLY");
  gMC->Gspos("FCBL", 4, "FAIC",  xcoor, ycoor,  zcoor, idrotm[2], "ONLY");

  if (fTOFHoles) {
    cblpar[3] = kCBLlh *0.5;
    cblpar[5] = fgkCBLh1*0.5 + kCBLlh*tgal;
    cblpar[7] = kCBLlh *0.5;
    cblpar[9] = cblpar[5];
    gMC->Gsvolu("FCBB", "TRAP", idtmed[511], cblpar, 11); // cables and tubes mix

    xcoor = (xtof - fgkCBLw)*0.5 - 2.*sawpar[0];
    ycoor = (fgkCBLh1 + 2.*cblpar[5])*0.25 - (ytof*0.5 - fgkModuleCoverThickness)*0.5;
    zcoor = kCBLl-kCBLlh*0.5;
    gMC->Gspos("FCBB", 1, "FAIB", -xcoor, ycoor, -zcoor, idrotm[1], "ONLY");
    gMC->Gspos("FCBB", 2, "FAIB",  xcoor, ycoor, -zcoor, idrotm[1], "ONLY");
    gMC->Gspos("FCBB", 3, "FAIB", -xcoor, ycoor,  zcoor, idrotm[2], "ONLY");
    gMC->Gspos("FCBB", 4, "FAIB",  xcoor, ycoor,  zcoor, idrotm[2], "ONLY");
  }

  // lateral cable and tube volume positioning
  xcoor = xtof*0.5 - sawpar[0];
  ycoor = (fgkCBLh2 - ytof*0.5 + fgkModuleCoverThickness)*0.5;
  zcoor = 0.;
  gMC->Gspos("FSAW", 1, "FAIA", -xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FSAW", 2, "FAIA",  xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FSAW", 1, "FAIC", -xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FSAW", 2, "FAIC",  xcoor, ycoor, zcoor, 0, "ONLY");

  if (fTOFHoles) {
    xcoor = xtof*0.5 - sawpar[0];
    ycoor = (fgkCBLh2 - ytof*0.5 + fgkModuleCoverThickness)*0.5;
    gMC->Gspos("FSAW", 1, "FAIB", -xcoor, ycoor, 0., 0, "ONLY");
    gMC->Gspos("FSAW", 2, "FAIB",  xcoor, ycoor, 0., 0, "ONLY");
  }

  // TOF Supermodule cover definition and positioning
  Float_t covpar[3] = {xtof*0.5, 0.075, zlenA*0.5};
  gMC->Gsvolu("FCOV", "BOX ", idtmed[504], covpar, 3); // Al
  if (fTOFHoles) {
    covpar[2] = (zlenA*0.5 - fgkInterCentrModBorder2)*0.5;
    gMC->Gsvolu("FCOB", "BOX ", idtmed[504], covpar, 3); // Al
    covpar[2] = fgkInterCentrModBorder2;
    gMC->Gsvolu("FCOP", "BOX ", idtmed[513], covpar, 3); // Plastic (CH2)
  }

  xcoor = 0.;
  ycoor = (ytof*0.5 - fgkModuleCoverThickness)*0.5 - covpar[1];
  zcoor = 0.;
  gMC->Gspos("FCOV", 0, "FAIA", xcoor, ycoor, zcoor, 0, "ONLY");
  gMC->Gspos("FCOV", 0, "FAIC", xcoor, ycoor, zcoor, 0, "ONLY");
  if (fTOFHoles) {
    zcoor = (zlenA*0.5 + fgkInterCentrModBorder2)*0.5;
    gMC->Gspos("FCOB", 1, "FAIB", xcoor, ycoor,  zcoor, 0, "ONLY");
    gMC->Gspos("FCOB", 2, "FAIB", xcoor, ycoor, -zcoor, 0, "ONLY");
    zcoor = 0.;
    gMC->Gspos("FCOP", 0, "FAIB", xcoor, ycoor,  zcoor, 0, "ONLY");
  }

}

//_____________________________________________________________________________
void AliTOFv6T0::MakeReadoutCrates(Float_t ytof) const
{
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
  //

  Int_t *idtmed = fIdtmed->GetArray()-499;

  Int_t idrotm[18];

  // volume definition
  Float_t serpar[3] = {29.*0.5, 121.*0.5, 90.*0.5};
  gMC->Gsvolu("FTOS", "BOX ", idtmed[514], serpar, 3); // Al + Cu + steel

  Float_t xcoor, ycoor, zcoor;
  zcoor = (118.-90.)*0.5;
  Float_t phi = -10.,  ra = fTOFGeometry->Rmin() + ytof*0.5;
  for (Int_t i = 0; i < fTOFGeometry->NSectors(); i++) {
    phi += 20.;
    xcoor = ra * TMath::Cos(phi * kDegrad);
    ycoor = ra * TMath::Sin(phi * kDegrad);
    AliMatrix(idrotm[i], 90., phi, 90., phi + 270., 0., 0.);
    gMC->Gspos("FTOS", i, "BFMO", xcoor, ycoor, zcoor, idrotm[i], "ONLY");
  }

  zcoor = (90. - 223.)*0.5;
  gMC->Gspos("FTOS", 1, "BBCE", ra, -3., zcoor, 0, "ONLY");

}

//_____________________________________________________________________________
void AliTOFv6T0::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //

  //AliTOF::CreateMaterials();

  AliMagF *magneticField = (AliMagF*)((AliMagF*)TGeoGlobalMagField::Instance()->GetField());

  Int_t   isxfld = magneticField->Integ();
  Float_t sxmgmx = magneticField->Max();

  //--- Quartz (SiO2) ---
  Float_t   aq[2] = { 28.0855,15.9994};
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  Float_t   dq = 2.7; // (+5.9%)
  Int_t nq = -2;

  // --- Nomex (C14H22O2N2) ---
  Float_t anox[4] = {12.011,1.00794,15.9994,14.00674};
  Float_t znox[4] = { 6.,  1.,  8.,  7.};
  Float_t wnox[4] = {14., 22., 2., 2.};
  //Float_t dnox  = 0.048; //old value
  Float_t dnox  = 0.22;    // (x 4.6)
  Int_t nnox   = -4;

  // --- G10  {Si, O, C, H, O} ---
  Float_t we[7], na[7];

  Float_t ag10[5] = {28.0855,15.9994,12.011,1.00794,15.9994};
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
  AliDebug(1,Form("wg10  %f  %f  %f  %f  %f", wmatg10[0], wmatg10[1], wmatg10[2], wmatg10[3], wmatg10[4]));
  //Float_t densg10 = 1.7; //old value
  Float_t densg10 = 2.0; // (+17.8%)

  // --- Water ---
  Float_t awa[2] = {  1.00794, 15.9994 };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;
  Int_t nwa = -2;

  // --- Air ---
  Float_t aAir[4]={12.011,14.00674,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir   = 1.20479E-3;

  // --- Fibre Glass ---
  Float_t afg[4] = {28.0855,15.9994,12.011,1.00794};
  Float_t zfg[4] = {14., 8., 6., 1.};
  Float_t wfg[4] = {0.12906,0.29405,0.51502,0.06187};
  //Float_t dfg    = 1.111;
  Float_t dfg    = 2.05; // (x1.845)
  Int_t nfg      = 4;

  // --- Freon C2F4H2 + SF6 ---
  Float_t afre[4] = {12.011,1.00794,18.9984032,32.0065};
  Float_t zfre[4] = { 6., 1., 9., 16.};
  Float_t wfre[4] = {0.21250,0.01787,0.74827,0.021355};
  Float_t densfre = 0.00375;
  Int_t nfre     = 4;

  // --- Cables and tubes {Al, Cu} ---
  Float_t acbt[2] = {26.981539,63.546};
  Float_t zcbt[2] = {13., 29.};
  Float_t wcbt[2] = {0.407,0.593};
  Float_t decbt   = 0.68;

  // --- Cable {CH2, Al, Cu} ---
  Float_t asc[4] = {12.011, 1.00794, 26.981539,63.546};
  Float_t zsc[4] = { 6., 1., 13., 29.};
  Float_t wsc[4];
  for (Int_t ii=0; ii<4; ii++) wsc[ii]=0.;

  Float_t wDummy[4], nDummy[4];
  for (Int_t ii=0; ii<4; ii++) wDummy[ii]=0.;
  for (Int_t ii=0; ii<4; ii++) nDummy[ii]=0.;
  nDummy[0] = 1.;
  nDummy[1] = 2.;
  MaterialMixer(wDummy,asc,nDummy,2);
  wsc[0] = 0.4375*wDummy[0];
  wsc[1] = 0.4375*wDummy[1];
  wsc[2] = 0.3244;
  wsc[3] = 0.2381;
  Float_t dsc = 1.223;

  // --- Crates boxes {Al, Cu, Fe, Cr, Ni} ---
  Float_t acra[5]= {26.981539,63.546,55.845,51.9961,58.6934};
  Float_t zcra[5]= {13., 29., 26., 24., 28.};
  Float_t wcra[5]= {0.7,0.2,0.07,0.018,0.012};
  Float_t dcra   = 0.77;

  // --- Polietilene CH2 ---
  Float_t aPlastic[2] = {12.011, 1.00794};
  Float_t zPlastic[2] = { 6., 1.};
  Float_t wPlastic[2] = { 1., 2.};
  //Float_t dPlastic = 0.92; // PDB value
  Float_t dPlastic = 0.93; // (~+1.1%)
  Int_t nwPlastic = -2;

  AliMixture ( 0, "Air$", aAir, zAir, dAir, 4, wAir);
  AliMixture ( 1, "Nomex$", anox, znox, dnox, nnox, wnox);
  AliMixture ( 2, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  AliMixture ( 3, "fibre glass$", afg, zfg, dfg, nfg, wfg);
  AliMaterial( 4, "Al $", 26.981539, 13., 2.7, -8.9, 999.);
  Float_t factor = 0.4/1.5*2./3.;
  AliMaterial( 5, "Al honeycomb$", 26.981539, 13., 2.7*factor, -8.9/factor, 999.);
  AliMixture ( 6, "Freon$", afre, zfre, densfre, nfre, wfre);
  AliMixture ( 7, "Glass$", aq, zq, dq, nq, wq);
  AliMixture ( 8, "Water$",  awa, zwa, dwa, nwa, wwa);
  AliMixture ( 9, "cables+tubes$", acbt, zcbt, decbt, 2, wcbt);
  AliMaterial(10, "Cu $", 63.546, 29., 8.96, -1.43, 999.);
  AliMixture (11, "cable$", asc, zsc, dsc, 4, wsc);
  AliMixture (12, "Al+Cu+steel$", acra, zcra, dcra, 5, wcra);
  AliMixture (13, "plastic$", aPlastic, zPlastic, dPlastic, nwPlastic, wPlastic);
  Float_t factorHoles = 1./36.5;
  AliMaterial(14, "Al honey for holes$", 26.981539, 13., 2.7*factorHoles, -8.9/factorHoles, 999.);

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

  AliMedium( 1,"Air$",          0, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 2,"Nomex$",        1, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 3,"G10$",          2, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 4,"fibre glass$",  3, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 5,"Al Frame$",     4, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 6,"honeycomb$",    5, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 7,"Fre$",          6, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 8,"Cu-S$",        10, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 9,"Glass$",        7, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(10,"Water$",        8, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(11,"Cable$",       11, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(12,"Cables+Tubes$", 9, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(13,"Copper$",      10, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(14,"Plastic$",     13, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(15,"Crates$",      12, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(16,"honey_holes$", 14, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);

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
     //&& gMC->GetMedium()==idtmed[507]
     && gMC->CurrentMedium()==idtmed[507]
     && gMC->CurrentVolID(copy)==fIdSens
     )
  {

    AliMC *mcApplication = (AliMC*)gAlice->GetMCApp();

    AddTrackReference(mcApplication->GetCurrentTrackNumber(), AliTrackReference::kTOF);
    //AddTrackReference(mcApplication->GetCurrentTrackNumber());

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
void AliTOFv6T0::MaterialMixer(Float_t * p, const Float_t * const a,
			       const Float_t * const m, Int_t n) const
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
