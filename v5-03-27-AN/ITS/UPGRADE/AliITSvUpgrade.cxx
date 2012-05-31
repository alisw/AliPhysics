/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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


/* $Id: AliITSvUpgrade.cxx */


//========================================================================
//
//        Geometry for the Upgrade of the Inner Tracking System
//
// Mario Sitta (sitta@to.infn.it)
//
//========================================================================



// $Log: AliITSvUpgrade.cxx,v $

#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TGeoVolume.h>
#include <TGeoXtru.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TVirtualMC.h>

#include "AliITSUpg.h"
#include "AliITSDetTypeSim.h"
#include "AliITShit.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliTrackReference.h"
#include "AliITSInitGeometryUpg.h"
#include "AliITSv11Geometry.h"
#include "AliITSv11GeometryUpgrade.h"
#include "AliITSv11GeomBeamPipe.h"
#include "AliITSvUpgrade.h"
#include "AliGeomManager.h"

const Double_t AliITSvUpgrade::fgkBeamPipeHalfZLen = 400;


ClassImp(AliITSvUpgrade)

//______________________________________________________________________
AliITSvUpgrade::AliITSvUpgrade():
  fMajorVersion(IsVersion()),
  fMinorVersion(-1),
  fNumberOfLayers(0),
  fLayTurbo(0),
  fLayRadii(0),
  fLayZLength(0),
  fLaddPerLay(0),
  fModPerLadd(0),
  fLadThick(0),
  fLadWidth(0),
  fLadTilt(0),
  fDetThick(0),
  fBeamPipe(0),
  fBeamPipeRmin(0),
  fBeamPipeRmax(0),
  fBeamPipeZlen(0),
  fInitGeom((AliITSVersion_t)fMajorVersion,fMinorVersion),
  fUpGeom(0),
  fBPGeom(0)
 {
    //    Standard default constructor
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
}

//______________________________________________________________________
AliITSvUpgrade::AliITSvUpgrade(const char *title):
  AliITSUpg("ITS", title),
  fMajorVersion(IsVersion()),
  fMinorVersion(-1),
  fNumberOfLayers(0),
  fLayTurbo(0),
  fLayRadii(0),
  fLayZLength(0),
  fLaddPerLay(0),
  fModPerLadd(0),
  fLadThick(0),
  fLadWidth(0),
  fLadTilt(0),
  fDetThick(0),
  fBeamPipe(0),
  fBeamPipeRmin(0),
  fBeamPipeRmax(0),
  fBeamPipeZlen(0),
  fInitGeom((AliITSVersion_t)fMajorVersion,fMinorVersion),
  fUpGeom(0),
  fBPGeom(0)
{
    //    Standard constructor for the Upgrade geometry.
    // Inputs:
    //   const char * title  Arbitrary title
    // Outputs:
    //   none.
    // Return:
    //   none.

  fIdN = 6;
  fIdName = new TString[fIdN];

  fIdName[0] = "ITS1";
  fIdName[1] = "ITS2";
  fIdName[2] = "ITS3";
  fIdName[3] = "ITS4";
  fIdName[4] = "ITS5";
  fIdName[5] = "ITS6";

  fIdSens    = new Int_t[fIdN];
  for(Int_t i=0; i<fIdN; i++) fIdSens[i] = 0;

  fLayTurbo   = new Bool_t[fIdN];
  fLayRadii   = new Double_t[fIdN];
  fLayZLength = new Double_t[fIdN];
  fLaddPerLay = new Int_t[fIdN];
  fModPerLadd = new Int_t[fIdN];
  fLadThick   = new Double_t[fIdN];
  fLadWidth   = new Double_t[fIdN];
  fLadTilt    = new Double_t[fIdN];
  fDetThick   = new Double_t[fIdN];

  fUpGeom = new AliITSv11GeometryUpgrade*[fIdN];

  if (fNumberOfLayers > 0)  // if not, we'll Fatal-ize in CreateGeometry
    for (Int_t j=0; j<fNumberOfLayers; j++) {
      fLayRadii[j] = 0.;
      fLayZLength[j] = 0.;
      fLaddPerLay[j] = 0;
      fModPerLadd[j] = 0;
      fLadWidth[j] = 0.;
      fDetThick[j] = 0.;
      fUpGeom[j] = 0;
    }
    
}

//______________________________________________________________________
AliITSvUpgrade::AliITSvUpgrade(const char *name, const char *title,
			       const Int_t nlay):
  AliITSUpg("ITS", title),
  fMajorVersion(IsVersion()),
  fMinorVersion(1),
  fNumberOfLayers(nlay),
  fLayTurbo(0),
  fLayRadii(0),
  fLayZLength(0),
  fLaddPerLay(0),
  fModPerLadd(0),
  fLadThick(0),
  fLadWidth(0),
  fLadTilt(0),
  fDetThick(0),
  fBeamPipe(0),
  fBeamPipeRmin(0),
  fBeamPipeRmax(0),
  fBeamPipeZlen(0),
  fInitGeom((AliITSVersion_t)fMajorVersion,fMinorVersion),
  fUpGeom(0),
  fBPGeom(0)
{
    //    Standard constructor for the Upgrade geometry.
    // Inputs:
    //   const char * name   Ignored, set to "ITS"
    //   const char * title  Arbitrary title
    //   const Int_t nlay    Number of layers
    // Outputs:
    //   none.
    // Return:
    //   none.
  
  fIdN = nlay;
  fIdName = new TString[fIdN];

  for (Int_t j=0; j<fIdN; j++)
    fIdName[j].Form("ITSupgSensor%d",j); // See AliITSv11GeometryUpgrade

  (void) name; // removes warning message

  fIdSens    = new Int_t[fIdN];
  for(Int_t i=0; i<fIdN; i++) fIdSens[i] = 0;

  fLayTurbo   = new Bool_t[fIdN];
  fLayRadii   = new Double_t[fIdN];
  fLayZLength = new Double_t[fIdN];
  fLaddPerLay = new Int_t[fIdN];
  fModPerLadd = new Int_t[fIdN];
  fLadThick   = new Double_t[fIdN];
  fLadWidth   = new Double_t[fIdN];
  fLadTilt    = new Double_t[fIdN];
  fDetThick   = new Double_t[fIdN];

  fUpGeom = new AliITSv11GeometryUpgrade*[fIdN];
  
  if (fNumberOfLayers > 0)  // if not, we'll Fatal-ize in CreateGeometry
    for (Int_t j=0; j<fNumberOfLayers; j++) {
      fLayRadii[j] = 0.;
      fLayZLength[j] = 0.;
      fLaddPerLay[j] = 0;
      fModPerLadd[j] = 0;
      fLadWidth[j] = 0.;
      fDetThick[j] = 0.;
      fUpGeom[j] = 0;
    }
    
}

//______________________________________________________________________
AliITSvUpgrade::~AliITSvUpgrade() {
    //    Standard destructor
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
  delete [] fLayTurbo;
  delete [] fLayRadii;
  delete [] fLayZLength;
  delete [] fLaddPerLay;
  delete [] fModPerLadd;
  delete [] fLadThick;
  delete [] fLadWidth;
  delete [] fLadTilt;
  delete [] fDetThick;
  delete [] fUpGeom;

}

//______________________________________________________________________
void AliITSvUpgrade::SetT2Lmatrix(Int_t uid, Double_t yShift, 
				  Bool_t yFlip, Bool_t yRot180) const {

  //
  // Creates the TGeo Local to Tracking transformation matrix
  // and sends it to the corresponding TGeoPNEntry 
  //
  // This function is used in AddAlignableVolumes()

  TGeoPNEntry *alignableEntry = gGeoManager->GetAlignableEntryByUID(uid);
  TGeoHMatrix* globMatrix = alignableEntry->GetGlobalOrig();

  Double_t *gtrans = globMatrix->GetTranslation(), rotMatrix[9];
  memcpy(&rotMatrix[0], globMatrix->GetRotationMatrix(), 9*sizeof(Double_t));
  Double_t al = TMath::ATan2(rotMatrix[1],rotMatrix[0]);
  if (yRot180) {
    al = TMath::ATan2(rotMatrix[1],-rotMatrix[0]);
  }
  Double_t xShift = gtrans[0]*TMath::Cos(al)+gtrans[1]*TMath::Sin(al);
  Double_t zShift = -gtrans[2];

  TGeoHMatrix *matLtoT = new TGeoHMatrix;
  matLtoT->SetDx( xShift ); // translation
  matLtoT->SetDy( yShift );
  matLtoT->SetDz( zShift );
  rotMatrix[0]= 0;  rotMatrix[1]= 1;  rotMatrix[2]= 0; // + rotation
  rotMatrix[3]= 1;  rotMatrix[4]= 0;  rotMatrix[5]= 0;
  rotMatrix[6]= 0;  rotMatrix[7]= 0;  rotMatrix[8]=-1;
  if (yFlip) rotMatrix[3] = -1;  // flipping in y  (for SPD1)
  if (yFlip) rotMatrix[1] = -1;  // flipping in y  (for SPD1)

  if (yRot180) { // rotation of pi around the axis perpendicular to the wafer
    if (yFlip) matLtoT->SetDx( -xShift ); // flipping in y  (for SPD1)
    matLtoT->SetDy( -yShift );
    matLtoT->SetDz( -zShift );
    rotMatrix[8]=1;
    rotMatrix[3] = -1;
    if (yFlip) rotMatrix[3] = 1;  // flipping in y  (for SPD1)
  }

  TGeoRotation rot;
  rot.SetMatrix(rotMatrix);
  matLtoT->MultiplyLeft(&rot);
  TGeoHMatrix *matTtoL = new TGeoHMatrix(matLtoT->Inverse());
  delete matLtoT;
  alignableEntry->SetMatrix(matTtoL);
}

//______________________________________________________________________
void AliITSvUpgrade::AddAlignableVolumes() const{
  // Creates entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path.
  // 
  // Records in the alignable entries the transformation matrices converting
  // TGeo local coordinates (in the RS of alignable volumes) to the tracking
  // system
  // For this, this function has to run before the misalignment because we
  // are using the ideal positions in the AliITSgeom object.
  // Inputs:
  //   none.
  // Outputs:
  //   none.
  // Return:
  //   none.

  AliInfo("Add ITS alignable volumes");

  if (!gGeoManager) {
    AliFatal("TGeoManager doesn't exist !");
    return;
  }
  // To be done - dummy for the time being
}

//______________________________________________________________________
void AliITSvUpgrade::AddBeamPipe(const Double_t rmin, const Double_t rmax,
				 const Double_t halfzlen) {

  // Define the parameters for the beam pipe
  if (fBeamPipe)
    AliWarning("Redefining beam pipe parameters");

  if (rmin <= 0) {
    AliError(Form("Beam pipe min radius (%f) wrong",rmin));
    return;
  } else
    fBeamPipeRmin = rmin;

  if (rmax <= 0) {
    AliError(Form("Beam pipe max radius (%f) wrong",rmax));
    return;
  } else
    fBeamPipeRmax = rmax;

  if (halfzlen < 0) {
    AliError(Form("Beam pipe half Zlength (%f) wrong",halfzlen));
    return;
  } else {
    if (halfzlen == 0) // Use default value
      fBeamPipeZlen = fgkBeamPipeHalfZLen;
    else
      fBeamPipeZlen = halfzlen;
  }

  fBeamPipe = kTRUE;
}

//______________________________________________________________________
void AliITSvUpgrade::CreateGeometry() {

  // Create the geometry and insert it in the mother volume ITSV


  TGeoManager *geoManager = gGeoManager;

  TGeoVolume *vALIC = geoManager->GetVolume("ALIC");

  new TGeoVolumeAssembly("ITSV");
  TGeoVolume *vITSV = geoManager->GetVolume("ITSV");
  vALIC->AddNode(vITSV, 2, 0);  // Copy number is 2 to cheat AliGeoManager::CheckSymNamesLUT

  //
  const Char_t *cvsDate="$Date: 2011-03-11 18:17:13 +0100 (Fri, 11 Mar 2011) $";
  const Char_t *cvsRevision="$Revision: 48336 $";
  const Int_t kLength=100;
  Char_t vstrng[kLength];
  if(fInitGeom.WriteVersionString(vstrng,kLength,(AliITSVersion_t)IsVersion(),
			     fMinorVersion,cvsDate,cvsRevision)) {
    vITSV->SetTitle(vstrng);
  }

  // Check that we have all needed parameters
  if (fNumberOfLayers <= 0)
    AliFatal(Form("Wrong number of layers (%d)",fNumberOfLayers));

  for (Int_t j=0; j<fNumberOfLayers; j++) {
    if (fLayRadii[j] <= 0)
      AliFatal(Form("Wrong layer radius for layer %d (%f)",j,fLayRadii[j]));
    if (fLayZLength[j] <= 0)
      AliFatal(Form("Wrong layer length for layer %d (%f)",j,fLayZLength[j]));
    if (fLaddPerLay[j] <= 0)
      AliFatal(Form("Wrong number of ladders for layer %d (%d)",j,
		    fLaddPerLay[j]));
    if (fModPerLadd[j] <= 0)
      AliFatal(Form("Wrong number of modules for layer %d (%d)",j,
		    fModPerLadd[j]));
    if (fLadThick[j] < 0)
      AliFatal(Form("Wrong ladder thickness for layer %d (%f)",j,
		    fLadThick[j]));
    if (fLayTurbo[j] && fLadWidth[j] <= 0)
      AliFatal(Form("Wrong ladder width for layer %d (%f)",j,
		    fLadWidth[j]));
    if (fDetThick[j] < 0)
      AliFatal(Form("Wrong module thickness for layer %d (%f)",j,
		    fDetThick[j]));

    if (j > 0) {
      if (fLayRadii[j] <= fLayRadii[j-1])
	AliError(Form("Layer %d radius (%f) is smaller than layer %d radius (%f)",
		      j,fLayRadii[j],j-1,fLayRadii[j-1]));
      if (fLayZLength[j] <= fLayZLength[j-1])
	AliWarning(Form("Layer %d length (%f) is smaller than layer %d length (%f)",
		      j,fLayZLength[j],j-1,fLayZLength[j-1]));
    } // if (j > 0)

    if (fLadThick[j] == 0)
      AliInfo(Form("Ladder thickness for layer %d not set, using default",j));

    if (fDetThick[j] == 0)
      AliInfo(Form("Module thickness for layer %d not set, using default",j));

  } // for (Int_t j=0; j<fNumberOfLayers; j++)

  // Now create the actual geometry
  for (Int_t j=0; j<fNumberOfLayers; j++) {
    if (fLayTurbo[j]) {
      fUpGeom[j] = new AliITSv11GeometryUpgrade(j,kTRUE,kFALSE);
      fUpGeom[j]->SetLadderWidth(fLadWidth[j]);
      fUpGeom[j]->SetLadderTilt(fLadTilt[j]);
    }
    else
      fUpGeom[j] = new AliITSv11GeometryUpgrade(j,kFALSE);
    fUpGeom[j]->SetRadius(fLayRadii[j]);
    fUpGeom[j]->SetZLength(fLayZLength[j]);
    fUpGeom[j]->SetNLadders(fLaddPerLay[j]);
    fUpGeom[j]->SetNModules(fModPerLadd[j]);
    if (fLadThick[j] != 0) fUpGeom[j]->SetLadderThick(fLadThick[j]);
    if (fDetThick[j] != 0) fUpGeom[j]->SetSensorThick(fDetThick[j]);
    fUpGeom[j]->CreateLayer(vITSV);
  }

  // Finally add the beam pipe
  if (fBeamPipe) {
    fBPGeom = new AliITSv11GeomBeamPipe(fBeamPipeRmin, fBeamPipeRmax,
					fBeamPipeZlen, kFALSE);
    fBPGeom->CreateBeamPipe(vALIC); // We put the BP in the ALIC volume
  }

}

//______________________________________________________________________
void AliITSvUpgrade::CreateMaterials(){
    // Create ITS materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSv11Hybrid.
    // In general it is automatically replaced by
    // the CreateMaterials routine defined in AliITSv?. Should the function
    // CreateMaterials not exist for the geometry version you are using this
    // one is used. See the definition found in AliITSv5 or the other routine
    // for a complete definition.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    Int_t   ifield = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
    Float_t fieldm = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();

    Float_t tmaxfd = 0.1; // 1.0; // Degree
    Float_t stemax = 1.0; // cm
    Float_t deemax = 0.1; // 30.0; // Fraction of particle's energy 0<deemax<=1
    Float_t epsil  = 1.0E-4; // 1.0; // cm
    Float_t stmin  = 0.0; // cm "Default value used"

    Float_t tmaxfdSi = 0.1; // .10000E+01; // Degree
    Float_t stemaxSi = 0.0075; //  .10000E+01; // cm
    Float_t deemaxSi = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilSi  = 1.0E-4;// .10000E+01;
    Float_t stminSi  = 0.0; // cm "Default value used"

    Float_t tmaxfdAir = 0.1; // .10000E+01; // Degree
    Float_t stemaxAir = .10000E+01; // cm
    Float_t deemaxAir = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilAir  = 1.0E-4;// .10000E+01;
    Float_t stminAir  = 0.0; // cm "Default value used"

    // AIR
    Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
    Float_t zAir[4]={6.,7.,8.,18.};
    Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
    Float_t dAir = 1.20479E-3;


    AliMaterial(1,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(1,"SI$",1,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(5,"AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(5,"AIR$",5,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMaterial(8,"BERILLIUM$",9.01, 4., 1.848, 35.3, 36.7);// From AliPIPEv3
    AliMedium(8,"BERILLIUM$",8,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

}

//______________________________________________________________________
void AliITSvUpgrade::DefineLayer(const Int_t nlay, const Double_t r,
				 const Double_t zlen, const Int_t nladd,
				 const Int_t nmod, const Double_t lthick,
				 const Double_t dthick){
    //     Sets the layer parameters
    // Inputs:
    //          nlay    layer number
    //          r       layer radius
    //          zlen    layer length
    //          nladd   number of ladders
    //          nmod    number of modules per ladder
    //          lthick  ladder thickness (if omitted, defaults to 0)
    //          dthick  detector thickness (if omitted, defaults to 0)
    // Outputs:
    //   none.
    // Return:
    //   none.

    if (nlay >= fNumberOfLayers || nlay < 0) {
      AliError(Form("Wrong layer number (%d)",nlay));
      return;
    }

    fLayTurbo[nlay] = kFALSE;
    fLayRadii[nlay] = r;
    fLayZLength[nlay] = zlen;
    fLaddPerLay[nlay] = nladd;
    fModPerLadd[nlay] = nmod;
    fLadThick[nlay] = lthick;
    fDetThick[nlay] = dthick;
}

//______________________________________________________________________
void AliITSvUpgrade::DefineLayerTurbo(const Int_t nlay, const Double_t r,
				      const Double_t zlen, const Int_t nladd,
				      const Int_t nmod, const Double_t width,
				      const Double_t tilt,
				      const Double_t lthick,
				      const Double_t dthick){
    //     Sets the layer parameters for a "turbo" layer
    //     (i.e. a layer whose ladders overlap in phi)
    // Inputs:
    //          nlay    layer number
    //          r       layer radius
    //          zlen    layer length
    //          nladd   number of ladders
    //          nmod    number of modules per ladder
    //          width   layer width
    //          tilt    layer tilt angle (degrees)
    //          lthick  ladder thickness (if omitted, defaults to 0)
    //          dthick  detector thickness (if omitted, defaults to 0)
    // Outputs:
    //   none.
    // Return:
    //   none.

    if (nlay >= fNumberOfLayers || nlay < 0) {
      AliError(Form("Wrong layer number (%d)",nlay));
      return;
    }

    fLayTurbo[nlay] = kTRUE;
    fLayRadii[nlay] = r;
    fLayZLength[nlay] = zlen;
    fLaddPerLay[nlay] = nladd;
    fModPerLadd[nlay] = nmod;
    fLadThick[nlay] = lthick;
    fLadWidth[nlay] = width;
    fLadTilt[nlay] = tilt;
    fDetThick[nlay] = dthick;
}

//______________________________________________________________________
void AliITSvUpgrade::GetBeamPipeParameters(Double_t &rmin, Double_t &rmax,
					   Double_t &hzlen){
    //     Gets the beam pipe parameters
    // Inputs:
    //   none.
    // Outputs:
    //          rmin    min radius
    //          rmax    max radius
    //          hzlen   half Z length
    // Return:
    //   none.

    rmin  = fBeamPipeRmin;
    rmax  = fBeamPipeRmax;
    hzlen = fBeamPipeZlen;

}

//______________________________________________________________________
void AliITSvUpgrade::GetLayerParameters(const Int_t nlay,
					Double_t &r, Double_t &zlen,
					Int_t &nladd, Int_t &nmod,
					Double_t &width, Double_t &tilt,
					Double_t &lthick, Double_t &dthick){
    //     Gets the layer parameters
    // Inputs:
    //          nlay    layer number
    // Outputs:
    //          r       layer radius
    //          zlen    layer length
    //          nladd   number of ladders
    //          nmod    number of modules per ladder
    //          width   ladder width
    //          tilt    ladder tilt angle
    //          lthick  ladder thickness
    //          dthick  detector thickness
    // Return:
    //   none.

    if (nlay >= fNumberOfLayers || nlay < 0) {
      AliError(Form("Wrong layer number (%d)",nlay));
      return;
    }

    r = fLayRadii[nlay];
    zlen = fLayZLength[nlay];
    nladd = fLaddPerLay[nlay];
    nmod = fModPerLadd[nlay];
    width = fLadWidth[nlay];
    tilt = fLadTilt[nlay];
    lthick = fLadThick[nlay];
    dthick = fDetThick[nlay];
}

//______________________________________________________________________
void AliITSvUpgrade::Init(){
    //     Initialise the ITS after it has been created.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    AliDebug(1,Form("Init: Major version %d Minor version %d",fMajorVersion,
		 fMinorVersion));
    UpdateInternalGeometry();
    AliITSUpg::Init();

}

//______________________________________________________________________
Bool_t AliITSvUpgrade::IsLayerTurbo(const Int_t nlay){
    //     Returns true if the layer is a "turbo" layer
    // Inputs:
    //          nlay    layer number
    // Outputs:
    //   none.
    // Return:
    //          kTRUE if the layer is a turbo layer

    if ( nlay < 0 || nlay > fNumberOfLayers ) {
      AliError(Form("Wrong layer number %d",nlay));
      return kFALSE;
    } else
      return fUpGeom[nlay]->IsTurbo();

}

//______________________________________________________________________
void AliITSvUpgrade::SetDefaults(){
    // sets the default segmentation, response, digit and raw cluster classes
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(!fDetTypeSim){
	Warning("SetDefaults","Error fDetTypeSim not defined");
	return;
    }

    fDetTypeSim->SetDefaults();

    return;
}

//______________________________________________________________________
void AliITSvUpgrade::StepManager(){
    //    Called for every step in the ITS, then calles the AliITShit class
    // creator with the information to be recoreded about that hit.
    //     The value of the macro ALIITSPRINTGEOM if set to 1 will allow the
    // printing of information to a file which can be used to create a .det
    // file read in by the routine CreateGeometry(). If set to 0 or any other
    // value except 1, the default behavior, then no such file is created nor
    // it the extra variables and the like used in the printing allocated.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(!(this->IsActive())) return;
    if(!(gMC->TrackCharge())) return;

    Int_t copy, lay = 0;
    Int_t id = gMC->CurrentVolID(copy);

    Bool_t notSens = kFALSE;
    while ((lay<fIdN)  && (notSens = id != fIdSens[lay])) ++lay;
    if (notSens) return;

    if(gMC->IsTrackExiting()) {
	AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kITS);
    } // if Outer ITS mother Volume

    static TLorentzVector position, momentum; // Saves on calls to construtors
    static AliITShit hit;// Saves on calls to constructors

    TClonesArray &lhits = *(Hits());
    Int_t   cpn0, cpn1, mod, status = 0;
    //
    // Track status
    if(gMC->IsTrackInside())      status +=  1;
    if(gMC->IsTrackEntering())    status +=  2;
    if(gMC->IsTrackExiting())     status +=  4;
    if(gMC->IsTrackOut())         status +=  8;
    if(gMC->IsTrackDisappeared()) status += 16;
    if(gMC->IsTrackStop())        status += 32;
    if(gMC->IsTrackAlive())       status += 64;

    //
    // retrieve the indices with the volume path
    //
    if (lay < 0 || lay > fIdN) {
      AliError(Form("Invalid value: lay=%d. Not an ITS sensitive volume",lay));
      return; // not an ITS sensitive volume.
    } else {
      copy = 1;
      gMC->CurrentVolOffID(1,cpn1);
      gMC->CurrentVolOffID(2,cpn0);
    } //

    fInitGeom.DecodeDetector(mod,lay,cpn0,cpn1,copy);
    // We should not need to pass by the switch !
    // This is time consuming...
    // therefore DecodeDetectorv11Upgrade(...) shouldn't be private !
    // and we should be able to use instead :
    //fInitGeom.DecodeDetectorv11Upgrade(mod,lay+1,cpn0,cpn1,copy);

    //
    // Fill hit structure.
    //
    hit.SetModule(mod);
    hit.SetTrack(gAlice->GetMCApp()->GetCurrentTrackNumber());
    gMC->TrackPosition(position);
    gMC->TrackMomentum(momentum);
    hit.SetPosition(position);
    hit.SetTime(gMC->TrackTime());
    hit.SetMomentum(momentum);
    hit.SetStatus(status);
    hit.SetEdep(gMC->Edep());
    hit.SetShunt(GetIshunt());
    if(gMC->IsTrackEntering()){
        hit.SetStartPosition(position);
        hit.SetStartTime(gMC->TrackTime());
        hit.SetStartStatus(status);
        return; // don't save entering hit.
    } // end if IsEntering
    // Fill hit structure with this new hit.
    //Info("StepManager","Calling Copy Constructor");
    new(lhits[fNhits++]) AliITShit(hit); // Use Copy Construtor.
    // Save old position... for next hit.
    hit.SetStartPosition(position);
    hit.SetStartTime(gMC->TrackTime());
    hit.SetStartStatus(status);

    return;
}

