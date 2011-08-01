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
#include "AliITSvUpgrade.h"
#include "AliGeomManager.h"

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
  fInitGeom((AliITSVersion_t)fMajorVersion,fMinorVersion),
  fUpGeom(0)
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
  fInitGeom((AliITSVersion_t)fMajorVersion,fMinorVersion),
  fUpGeom(0)
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
  fInitGeom((AliITSVersion_t)fMajorVersion,fMinorVersion),
  fUpGeom(0)
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

    Float_t tmaxfdServ = 1.0; // 10.0; // Degree
    Float_t stemaxServ = 1.0; // 0.01; // cm
    Float_t deemaxServ = 0.5; // 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilServ  = 1.0E-3; // 0.003; // cm
    Float_t stminServ  = 0.0; //0.003; // cm "Default value used"

    // Freon PerFluorobuthane C4F10 see 
    // http://st-support-cooling-electronics.web.cern.ch/
    //        st-support-cooling-electronics/default.htm
    Float_t afre[2]  = { 12.011,18.9984032 };
    Float_t zfre[2]  = { 6., 9. };
    Float_t wfre[2]  = { 4.,10. };
    Float_t densfre  = 1.52;


    //CM55J

    Float_t aCM55J[4]={12.0107,14.0067,15.9994,1.00794};
    Float_t zCM55J[4]={6.,7.,8.,1.};
    Float_t wCM55J[4]={0.908508078,0.010387573,0.055957585,0.025146765};
    Float_t dCM55J = 1.8;

    //ALCM55J

    Float_t aALCM55J[5]={12.0107,14.0067,15.9994,1.00794,26.981538};
    Float_t zALCM55J[5]={6.,7.,8.,1.,13.};
    Float_t wALCM55J[5]={0.817657902,0.0093488157,0.0503618265,0.0226320885,0.1};
    Float_t dALCM55J = 1.9866;

    //Si Chips

    Float_t aSICHIP[6]={12.0107,14.0067,15.9994,1.00794,28.0855,107.8682};
    Float_t zSICHIP[6]={6.,7.,8.,1.,14., 47.};
    Float_t wSICHIP[6]={0.039730642,0.001396798,0.01169634,0.004367771,0.844665,0.09814344903};
    Float_t dSICHIP = 2.36436;

    //Inox
    
    Float_t aINOX[9]={12.0107,54.9380, 28.0855,30.9738,32.066,58.6928,51.9961,95.94,55.845};
    Float_t zINOX[9]={6.,25.,14.,15.,16., 28.,24.,42.,26.};
    Float_t wINOX[9]={0.0003,0.02,0.01,0.00045,0.0003,0.12,0.17,0.025,0.654};
    Float_t dINOX = 8.03;

    //AISI 304 L (from F.Tosello's web page - M.S. 18 Oct 10)
    
    Float_t a304L[8]={12.0107,54.9380, 28.0855,30.9738,32.066,58.6928,51.9961,55.845};
    Float_t z304L[8]={6.,25.,14.,15.,16., 28.,24.,26.};
    Float_t w304L[8]={0.0003,0.02,0.01,0.00045,0.003,0.0925,0.19,0.6865};
    Float_t d304L = 8.03;

    //SDD HV microcable

    Float_t aHVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zHVm[5]={6.,1.,7.,8.,13.};
    Float_t wHVm[5]={0.520088819984,0.01983871336,0.0551367996,0.157399667056, 0.247536};
    Float_t dHVm = 1.6087;

    //SDD LV+signal cable

    Float_t aLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zLVm[5]={6.,1.,7.,8.,13.};
    Float_t wLVm[5]={0.21722436468,0.0082859922,0.023028867,0.06574077612, 0.68572};
    Float_t dLVm = 2.1035;

    //SDD hybrid microcab

    Float_t aHLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zHLVm[5]={6.,1.,7.,8.,13.};
    Float_t wHLVm[5]={0.24281879711,0.00926228815,0.02574224025,0.07348667449, 0.64869};
    Float_t dHLVm = 2.0502;

    //SDD anode microcab

    Float_t aALVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zALVm[5]={6.,1.,7.,8.,13.};
    Float_t wALVm[5]={0.392653705471,0.0128595919215,0.041626868025,0.118832707289, 0.431909};
    Float_t dALVm = 2.0502;

    //X7R capacitors - updated from F.Tosello's web page - M.S. 18 Oct 10

    Float_t aX7R[6]={137.327,47.867,15.9994,58.6928,63.5460,118.710};
    Float_t zX7R[6]={56.,22.,8.,28.,29.,50.};
    Float_t wX7R[6]={0.524732,0.176736,0.179282,0.079750,0.019750,0.019750};
    Float_t dX7R = 6.07914;

    //X7R weld, i.e. Sn 60% Pb 40% (from F.Tosello's web page - M.S. 15 Oct 10)

    Float_t aX7Rweld[2]={118.71 , 207.20};
    Float_t zX7Rweld[2]={ 50.   ,  82.  };
    Float_t wX7Rweld[2]={  0.60 ,   0.40};
    Float_t dX7Rweld   = 8.52358;

    // AIR

    Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
    Float_t zAir[4]={6.,7.,8.,18.};
    Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
    Float_t dAir = 1.20479E-3;

    // Water

    Float_t aWater[2]={1.00794,15.9994};
    Float_t zWater[2]={1.,8.};
    Float_t wWater[2]={0.111894,0.888106};
    Float_t dWater   = 1.0;

    // CERAMICS
  //     94.4% Al2O3 , 2.8% SiO2 , 2.3% MnO , 0.5% Cr2O3
    Float_t acer[5]  = { 26.981539,15.9994,28.0855,54.93805,51.9961 };
    Float_t zcer[5]  = {       13.,     8.,    14.,     25.,    24. };
    Float_t wcer[5]  = {.4443408,.5213375,.0130872,.0178135,.003421};
    Float_t denscer  = 3.6;

    //G10FR4

    Float_t zG10FR4[14] = {14.00,	20.00,	13.00,	12.00,	5.00,	22.00,	11.00,	19.00,	26.00,	9.00,	8.00,	6.00,	7.00,	1.00};
    Float_t aG10FR4[14] = {28.0855000,40.0780000,26.9815380,24.3050000,10.8110000,47.8670000,22.9897700,39.0983000,55.8450000,18.9984000,15.9994000,12.0107000,14.0067000,1.0079400};
    Float_t wG10FR4[14] = {0.15144894,0.08147477,0.04128158,0.00904554,0.01397570,0.00287685,0.00445114,0.00498089,0.00209828,0.00420000,0.36043788,0.27529426,0.01415852,0.03427566};
    Float_t densG10FR4= 1.8;
    
     //--- EPOXY  --- C18 H19 O3
      Float_t aEpoxy[3] = {15.9994, 1.00794, 12.0107} ; 
      Float_t zEpoxy[3] = {     8.,      1.,      6.} ; 
      Float_t wEpoxy[3] = {     3.,     19.,     18.} ; 
      Float_t dEpoxy = 1.8 ;

      // rohacell: C9 H13 N1 O2
    Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
    Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
    Float_t wrohac[4] = { 14.,   10.,    2.,     6.};
    Float_t drohac    = 0.052;

    // If he/she means stainless steel (inox) + Aluminium and Zeff=15.3383 then
//
// %Al=81.6164 %inox=100-%Al

    Float_t aInAl[5] = {27., 55.847,51.9961,58.6934,28.0855 };
    Float_t zInAl[5] = {13., 26.,24.,28.,14. };
    Float_t wInAl[5] = {.816164, .131443,.0330906,.0183836,.000919182};
    Float_t dInAl    = 3.075;

    // Aluminum alloy with 12% Copper - 21 Oct 10

    Float_t aAlCu12[2] = {26.9815, 63.546};
    Float_t zAlCu12[2] = {13.    , 29.   };
    Float_t wAlCu12[2] = { 0.88  ,  0.12 };
    Float_t dAlCu12    = 2.96;

    // Kapton

    Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
    Float_t zKapton[4]={1.,6.,7.,8.};
    Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
    Float_t dKapton   = 1.42;
    
    // Kapton + Cu (for Pixel Bus)

    Float_t aKaptonCu[5]={1.00794, 12.0107, 14.010, 15.9994, 63.5460};
    Float_t zKaptonCu[5]={1., 6., 7., 8., 29.};
    Float_t wKaptonCuBus[5];
    
    // Kapton + Cu (for Pixel MCM)

    Float_t wKaptonCuMCM[5];
    
    // Kapton + Cu (mix of two above)

    Float_t wKaptonCuMix[5];

    //SDD ruby sph.
    Float_t aAlOxide[2]  = { 26.981539,15.9994};
    Float_t zAlOxide[2]  = {       13.,     8.};
    Float_t wAlOxide[2]  = {0.4707, 0.5293};
    Float_t dAlOxide     = 3.97;

    // Silica for optical fibers: Si O2
    Float_t aoptfib[2] = { 28.0855, 15.9994};
    Float_t zoptfib[2] = { 14.,      8.    };
    Float_t woptfib[2] = {  1.,      2.    };
    Float_t doptfib    = 2.55;

    // Tetrafluorethylene-Perfluorpropylene (FEP) - 08 Mar 10
    Float_t aFEP[2] = { 12.0107, 18.9984};
    Float_t zFEP[2] = {  6.    ,  9.    };
    Float_t wFEP[2] = {  1.    ,  2.    };
    Float_t dFEP    = 2.15;

    // PVC (C2H3Cl)n - 08 Jul 10
    Float_t aPVC[3] = { 12.0107, 1.00794, 35.4527};
    Float_t zPVC[3] = {  6.    , 1.     , 35.   };
    Float_t wPVC[3] = {  2.    , 3.     ,  1.   };
    Float_t dPVC    = 1.3;

    // PBT (Polybutylene terephthalate = C12-H12-O4) - 01 Sep 10
    Float_t aPBT[3] = { 12.0107, 1.00794, 15.9994};
    Float_t zPBT[3] = {  6.    , 1.     ,  8.   };
    Float_t wPBT[3] = { 12.    ,12.     ,  4.   };
    Float_t dPBT    = 1.31;

    // POLYAX (POLYAX = C37-H24-O6-N2) - 03 Sep 10
    Float_t aPOLYAX[4] = { 12.0107, 1.00794, 15.9994, 14.00674};
    Float_t zPOLYAX[4] = {  6.    , 1.     ,  8.    ,  7.     };
    Float_t wPOLYAX[4] = { 37.    ,24.     ,  6.    ,  2.     };
    Float_t dPOLYAX    = 1.27;

    // PPS (PPS = C6-H4-S) - 05 Sep 10
    Float_t aPPS[3] = { 12.0107, 1.00794, 32.066};
    Float_t zPPS[3] = {  6.    , 1.     , 16.   };
    Float_t wPPS[3] = {  6.    , 4.     ,  1.   };
    Float_t dPPS    = 1.35;

    // Megolon (Polyolefin = (C-H2)n) - 20 Oct 10
    Float_t aMegolon[2] = { 12.0107, 1.00794};
    Float_t zMegolon[2] = {  6.    , 1.     };
    Float_t wMegolon[2] = {  1.    , 2.     };
    Float_t dMegolon    = 1.51; // Mean of various types

    // Standard glass (from glassproperties.com/glasses - M.S. 21 Oct 10)
    Float_t aStdGlass[7] = {15.9994  ,28.0855  ,22.98977 ,40.078   ,
			    24.305   ,26.981539,39.0983  };
    Float_t zStdGlass[7] = { 8.      ,14.      ,11.      ,20.      ,
			    12.      ,13.      ,19.      };
    Float_t wStdGlass[7] = { 0.468377, 0.348239, 0.096441, 0.071469,
			     0.006030, 0.005293, 0.004151};
    Float_t dStdGlass    = 2.53;

    // Glass Fiber (from F.Tosello's web page - M.S. 15 Oct 10)
    Float_t aGlass[11] = {15.9994  ,28.0855  ,40.078   ,26.981539,10.811   ,
		24.305   ,39.0983  ,22.98977 ,18.9984  ,47.867   ,55.845};
    Float_t zGlass[11] = { 8.      ,14.      ,20       ,13       , 5       ,
		12.      ,19       ,11       , 9       ,22       ,26    };
    Float_t wGlass[11] = { 0.473610, 0.252415, 0.135791, 0.068803, 0.023293,
		 0.015076, 0.008301, 0.007419, 0.007000, 0.004795, 0.003497};
    Float_t dGlass = 2.61;

    // Ryton R-4 04 (from F.Tosello's web page - M.S. 15 Oct 10)
    Float_t aRyton[14] = {15.9994  ,28.0855  ,40.078   ,26.981539,10.811   ,
			  24.305   ,39.0983  ,22.98977 ,18.9984  ,47.867   ,
			  55.845   ,12.0107  , 1.00794 ,32.066   };
    Float_t zRyton[14] = { 8.      ,14.      ,20.      ,13.      , 5.      ,
			  12.      ,19.      ,11.      , 9.      ,22.      ,
			  26.      , 6.      , 1.      ,16.      };
    Float_t wRyton[14] = { 0.189445, 0.100966, 0.054316, 0.027521, 0.009317,
			   0.006030, 0.003320, 0.002968, 0.002800, 0.001918,
			   0.001399, 0.399760, 0.022365, 0.177875};
    Float_t dRyton = 1.65;

    // Plexiglas (Poly(methyl methacrylate) (C5O2H8)n - M.S. 05 nov 10)
    Float_t aPlexy[3] = { 12.0107, 15.9994,  1.00794};
    Float_t zPlexy[3] = {  6.    , 8.     ,  1.   };
    Float_t wPlexy[3] = {  5.    , 2.     ,  8.   };
    Float_t dPlexy    = 1.18;

    //SSD NiSn capacitor ends
    Float_t aNiSn[2]  = { 56.6934,118.710};
    Float_t zNiSn[2]  = {     28.,     50.};
    Float_t wNiSn[2]  = {0.33, 0.67};
    Float_t dNiSn     = wNiSn[0]*8.908 + wNiSn[1]*7.310;

    AliMaterial(1,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(1,"SI$",1,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(2,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(2,"SPD SI CHIP$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(3,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(3,"SPD SI BUS$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(4,"C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(4,"C (M55J)$",4,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(5,"AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(5,"AIR$",5,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(6,"GEN AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(6,"GEN AIR$",6,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(7,"SDD SI CHIP$",aSICHIP,zSICHIP,dSICHIP,6,wSICHIP);
    AliMedium(7,"SDD SI CHIP$",7,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(9,"SDD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(9,"SDD C (M55J)$",9,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(10,"SDD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(10,"SDD AIR$",10,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMaterial(11,"AL$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
    AliMedium(11,"AL$",11,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(12, "Water$",aWater,zWater,dWater,2,wWater);
    AliMedium(12,"WATER$",12,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(13,"Freon$",afre,zfre,densfre,-2,wfre);
    AliMedium(13,"Freon$",13,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(14,"COPPER$",0.63546E+02,0.29000E+02,0.89600E+01,0.14300E+01,0.99900E+03);
    AliMedium(14,"COPPER$",14,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    AliMixture(15,"CERAMICS$",acer,zcer,denscer,5,wcer);
    AliMedium(15,"CERAMICS$",15,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(20,"SSD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(20,"SSD C (M55J)$",20,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(21,"SSD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(21,"SSD AIR$",21,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(25,"G10FR4$",aG10FR4,zG10FR4,densG10FR4,14,wG10FR4);
    AliMedium(25,"G10FR4$",25,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMixture(26,"GEN C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(26,"GEN C (M55J)$",26,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(27,"GEN Air$",aAir,zAir,dAir,4,wAir);
    AliMedium(27,"GEN Air$",27,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(35,"PLEXYGLAS$",aPlexy,zPlexy,dPlexy,-3,wPlexy);
    AliMedium(35,"PLEXYGLAS$",35,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(36,"STDGLASS$",aStdGlass,zStdGlass,dStdGlass,7,wStdGlass);
    AliMedium(36,"STDGLASS$",36,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(37,"ALCU12$",aAlCu12,zAlCu12,dAlCu12,2,wAlCu12);
    AliMedium(37,"ALCU12$",37,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(38,"MEGOLON$",aMegolon,zMegolon,dMegolon,-2,wMegolon);
    AliMedium(38,"MEGOLON$",38,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(39,"RYTON$",aRyton,zRyton,dRyton,14,wRyton);
    AliMedium(39,"RYTON$",39,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(40,"GLASS FIBER$",aGlass,zGlass,dGlass,11,wGlass);
    AliMedium(40,"GLASS FIBER$",40,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(41,"AISI304L$",a304L,z304L,d304L,8,w304L);
    AliMedium(41,"AISI304L$",41,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(42,"NICKEL$",0.58693E+02,0.28000E+02,0.89080E+01,0.14200E+01,0.99900E+03);
    AliMedium(42,"NICKEL$",42,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
               
    AliMixture(43,"SDD X7R weld$",aX7Rweld,zX7Rweld,dX7Rweld,2,wX7Rweld);
    AliMedium(43,"SDD X7R weld$",43,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(44,"PPS$",aPPS,zPPS,dPPS,-3,wPPS);
    AliMedium(44,"PPS$",44,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(45,"POLYAX$",aPOLYAX,zPOLYAX,dPOLYAX,-4,wPOLYAX);
    AliMedium(45,"POLYAX$",45,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(46,"PBT$",aPBT,zPBT,dPBT,-3,wPBT);
    AliMedium(46,"PBT$",46,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(47,"PVC$",aPVC,zPVC,dPVC,-3,wPVC);
    AliMedium(47,"PVC$",47,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    Double_t cuFrac = 0.56;
    Double_t kFrac  = 1.0 - cuFrac;
    Double_t cuDens = 8.96;
    Float_t dKaptonCuBus   = cuFrac * cuDens + kFrac * dKapton;
    for (Int_t j=0; j<4; j++)
      wKaptonCuBus[j] = wKapton[j]*kFrac;
    wKaptonCuBus[4] = cuFrac;
    AliMixture(48, "SPD-BUS CU KAPTON", aKaptonCu, zKaptonCu, dKaptonCuBus, 5, wKaptonCuBus);
    AliMedium(48,"SPD-BUS CU KAPTON$",48,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    
    cuFrac = 0.5;
    kFrac  = 1.0 - cuFrac;
    Float_t dKaptonCuMCM   = cuFrac * cuDens + kFrac * dKapton;
    for (Int_t j=0; j<4; j++)
      wKaptonCuMCM[j] = wKapton[j]*kFrac;
    wKaptonCuMCM[4] = cuFrac;
    AliMixture(49, "SPD-MCM CU KAPTON", aKaptonCu, zKaptonCu, dKaptonCuMCM, 5, wKaptonCuMCM);
    AliMedium(49,"SPD-MCM CU KAPTON$",49,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    
    cuFrac = (0.56 + 0.5) / 2.0;
    kFrac  = 1.0 - cuFrac;
    Float_t dKaptonCuMix   = cuFrac * cuDens + kFrac * dKapton;
    for (Int_t j=0; j<4; j++)
      wKaptonCuMix[j] = wKapton[j]*kFrac;
    wKaptonCuMix[4] = cuFrac;
    AliMixture(50, "SPD-MIX CU KAPTON", aKaptonCu, zKaptonCu, dKaptonCuMix, 5, wKaptonCuMix);
    AliMedium(50,"SPD-MIX CU KAPTON$",50,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(51,"SPD SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(51,"SPD SI$",51,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(52,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(52,"SPD SI CHIP$",52,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(53,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(53,"SPD SI BUS$",53,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(54,"SPD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(54,"SPD C (M55J)$",54,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(55,"SPD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(55,"SPD AIR$",55,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(56, "SPD KAPTON(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(56,"SPD KAPTON(POLYCH2)$",56,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // Gaseous Freon has same chemical composition but air density at 1.7 atm
    AliMixture(59,"GASEOUS FREON$",afre,zfre,1.7*dAir,-2,wfre);
    AliMedium(59,"GASEOUS FREON$",59,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(61,"EPOXY$",aEpoxy,zEpoxy,dEpoxy,-3,wEpoxy);
    AliMedium(61,"EPOXY$",61,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(62,"SILICON$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(62,"SILICON$",62,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(63, "KAPTONH(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(63,"KAPTONH(POLYCH2)$",63,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(64,"ALUMINUM$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
    AliMedium(64,"ALUMINUM$",64,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(65,"INOX$",aINOX,zINOX,dINOX,9,wINOX);
    AliMedium(65,"INOX$",65,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(66,"NiSn$",aNiSn,zNiSn,dNiSn,2,wNiSn);
    AliMedium(66,"NiSn$",66,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(67,"Sn$", 118.710, 50., 7.310, 1.206, 999.);
    AliMedium(67,"Sn$",67,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(68,"ROHACELL$",arohac,zrohac,drohac,-4,wrohac);
    AliMedium(68,"ROHACELL$",68,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMixture(69,"SDD C AL (M55J)$",aALCM55J,zALCM55J,dALCM55J,5,wALCM55J);
    AliMedium(69,"SDD C AL (M55J)$",69,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
  
    AliMixture(70, "SDDKAPTON (POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(70,"SDDKAPTON (POLYCH2)$",70,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMaterial(71,"ITS SANDW A$",0.12011E+02,0.60000E+01,0.2115E+00,0.17479E+03,0.99900E+03);
    AliMedium(71,"ITS SANDW A$",71,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(72,"ITS SANDW B$",0.12011E+02,0.60000E+01,0.27000E+00,0.18956E+03,0.99900E+03);
    AliMedium(72,"ITS SANDW B$",72,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(73,"ITS SANDW C$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    AliMedium(73,"ITS SANDW C$",73,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(74,"HEAT COND GLUE$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
    AliMedium(74,"HEAT COND GLUE$",74,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(75,"ELASTO SIL$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(75,"ELASTO SIL$",75,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // SPD bus (data from Petra Riedler)
    Float_t aSPDbus[5] = {1.00794,12.0107,14.01,15.9994,26.982 };
    Float_t zSPDbus[5] = {1.,6.,7.,8.,13.};
    Float_t wSPDbus[5] = {0.023523,0.318053,0.009776,0.078057,0.570591};
    Float_t dSPDbus    = 2.128505;

    //   AliMaterial(76,"SPDBUS(AL+KPT+EPOX)$",0.19509E+02,0.96502E+01,0.19060E+01,0.15413E+02,0.99900E+03);
    AliMixture(76,"SPDBUS(AL+KPT+EPOX)$",aSPDbus,zSPDbus,dSPDbus,5,wSPDbus);
    AliMedium(76,"SPDBUS(AL+KPT+EPOX)$",76,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
               
    AliMixture(77,"SDD X7R capacitors$",aX7R,zX7R,dX7R,6,wX7R);
    AliMedium(77,"SDD X7R capacitors$",77,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(78,"SDD ruby sph. Al2O3$",aAlOxide,zAlOxide,dAlOxide,2,wAlOxide);
    AliMedium(78,"SDD ruby sph. Al2O3$",78,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(79,"SDD SI insensitive$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(79,"SDD SI insensitive$",79,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(80,"SDD HV microcable$",aHVm,zHVm,dHVm,5,wHVm);
    AliMedium(80,"SDD HV microcable$",80,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(81,"SDD LV+signal cable$",aLVm,zLVm,dLVm,5,wLVm);
    AliMedium(81,"SDD LV+signal cable$",81,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(82,"SDD hybrid microcab$",aHLVm, zHLVm,dHLVm,5,wHLVm);
    AliMedium(82,"SDD hybrid microcab$",82,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(83,"SDD anode microcab$",aALVm,zALVm,dALVm,5,wALVm);
    AliMedium(83,"SDD anode microcab$",83,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    Float_t aDSring[4]={12.0107,      1.00794,     14.0067,      15.9994};
    Float_t zDSring[4]={ 6.,          1.,           7.,           8.};
    Float_t wDSring[4]={ 0.854323888, 0.026408778,  0.023050265,  0.096217069};
    Float_t dDSring = 0.2875;
    AliMixture(84,"SDD/SSD rings$",aDSring,zDSring,dDSring,4,wDSring);
    AliMedium(84,"SDD/SSD rings$",84,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(85,"inox/alum$",aInAl,zInAl,dInAl,5,wInAl);
    AliMedium(85,"inox/alum$",85,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // Special media

    AliMaterial(90,"SPD shield$", 12.011, 6., 1.93 , 22.36, 999);
    AliMedium(90,"SPD shield$",90,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    // SPD End Ladder (data from Petra Riedler)
    Float_t aSPDel[5] = {1.00794,12.0107,14.01,15.9994,63.54 };
    Float_t zSPDel[5] = {1.,6.,7.,8.,29.};
    Float_t wSPDel[5] = {0.004092,0.107274,0.011438,0.032476,0.844719};
    Float_t dSPDel    = 3.903403;

    //   AliMaterial(91, "SPD End ladder$", 47.0447, 21.7963, 3.6374, 4.4711, 999); 
    AliMixture(91,"SPD End ladder$",aSPDel,zSPDel,dSPDel,5,wSPDel);
    AliMedium(91,"SPD End ladder$",91,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    AliMaterial(92, "SPD cone$",28.0855, 14., 2.33, 9.36, 999);    
    AliMedium(92,"SPD cone$",92,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    /*  Material with fractional Z not actually used
    AliMaterial(93, "SDD End ladder$", 69.9298, 29.8246, 0.3824, 36.5103, 999);
    AliMedium(93,"SDD End ladder$",93,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    */
    AliMaterial(94, "SDD cone$",63.546, 29., 1.15, 1.265, 999);
    AliMedium(94,"SDD cone$",94,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    /* Material with fractional Z not actually used
    AliMaterial(95, "SSD End ladder$", 32.0988, 15.4021, 0.68, 35.3238, 999); 
    AliMedium(95,"SSD End ladder$",95,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    */
    AliMaterial(96, "SSD cone$",63.546, 29., 1.15, 1.265, 999);
    AliMedium(96,"SSD cone$",96,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    AliMixture(98,"SDD OPTICFIB$",aoptfib,zoptfib,doptfib,-2,woptfib);
    AliMedium(98,"SDD OPTICFIB$",98,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(95,"SSD FEP$",aFEP,zFEP,dFEP,-2,wFEP);
    AliMedium(95,"SSD FEP$",95,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    Float_t aA[13],zZ[13],wW[13],den;
    // Mean material for low-voltage cables on SPD trays Side A
    // (Copper + PolyEthylene (C2-H4)) (D.Elia for cable number and
    // cross-section area, M.Sitta for elemental computation) - 26 Feb 10
    wW[0] = 0.323024;//H
    wW[2] = 0.515464;//Cu
    wW[1] = 0.161512;//C
    wW[3] = 0.000000;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 5.078866;
    AliMixture(60,"SPD_LOWCABLES$",aA,zZ,den,+3,wW);
    AliMedium(60,"SPD_LOWCABLES$",60,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    // Mean material for high-voltage cables on SPD trays Side A & C
    // (Copper + HD PolyEthylene (C2-H2)) (D.Elia for cable number and
    // cross-section area, M.Sitta for elemental computation) - 10 Jun 10
    wW[0] = 0.083766;//H
    wW[2] = 0.417136;//Cu
    wW[1] = 0.499098;//C
    wW[3] = 0.000000;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 1.514930;
    AliMixture(58,"SPD_HICABLES$",aA,zZ,den,+3,wW);
    AliMedium(58,"SPD_HICABLES$",58,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    // PolyUrethane [C25-H42-N2-O6] - 07 Mar 10
    zZ[2] =  7.0; aA[2] =  14.0067; // Nitrogen - From Root TGeoElementTable

    wW[0] = 0.090724;//H
    wW[2] = 0.060035;//N
    wW[1] = 0.643513;//C
    wW[3] = 0.205728;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 1.158910;
    AliMixture(67,"POLYURETHANE$",aA,zZ,den,+4,wW);
    AliMedium(67,"POLYURETHANE$",67,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    //  POM (Polyoxymethylene = (CH2O)n ) - 02 May 10
    zZ[2] =  8.0; aA[2] =  15.9994; // Oxigen

    wW[0] = 0.067137;//H
    wW[1] = 0.400016;//C
    wW[2] = 0.532847;//O
    wW[3] = 0.000000;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 1.4200;
    AliMixture(57,"POLYOXYMETHYLENE$",aA,zZ,den,+3,wW);
    AliMedium(57,"POLYOXYMETHYLENE$",57,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);


    // Anticorodal (Aliminum alloy) - 08 nov 10
    // A,Z from Root TGeoElementTable, W from Web sites
    zZ[0] = 13.0; aA[0] =  26.9815; // Aluminium
    zZ[1] = 29.0; aA[1] =  63.546 ; // Copper
    zZ[2] = 26.0; aA[2] =  55.845 ; // Iron
    zZ[3] = 25.0; aA[3] =  54.938 ; // Manganese
    zZ[4] = 12.0; aA[4] =  24.305 ; // Magnesium
    zZ[5] = 14.0; aA[5] =  28.0855; // Silicon
    zZ[6] = 30.0; aA[6] =  65.39  ; // Zinc
    zZ[7] = 24.0; aA[7] =  51.9961; // Chromium
    zZ[8] = 22.0; aA[8] =  47.867 ; // Titanium

    wW[1] = 0.001000;//Cu
    wW[2] = 0.005000;//Fe
    wW[3] = 0.007000;//Mn - mean value
    wW[4] = 0.009000;//Mg - mean value
    wW[5] = 0.001000;//Si - mean value
    wW[6] = 0.002000;//Zn
    wW[7] = 0.002500;//Cr
    wW[8] = 0.001000;//Ti

    Double_t totFrac = 0;
    for (Int_t j=1; j<9; j++)
      totFrac += wW[j];
    wW[0] = 1. - totFrac;//Al - the remainder

    den = 2.69;
    AliMixture(93,"ANTICORODAL$",aA,zZ,den,+9,wW);
    AliMedium(93,"ANTICORODAL$",93,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // Hokotol (another Aluminium alloy) - 08 nov 10
    // A,Z from Root TGeoElementTable, W from Web sites
    zZ[0] = 13.0; aA[0] =  26.9815; // Aluminium
    zZ[1] = 29.0; aA[1] =  63.546 ; // Copper
    zZ[2] = 26.0; aA[2] =  55.845 ; // Iron
    zZ[3] = 25.0; aA[3] =  54.938 ; // Manganese
    zZ[4] = 12.0; aA[4] =  24.305 ; // Magnesium
    zZ[5] = 14.0; aA[5] =  28.0855; // Silicon
    zZ[6] = 30.0; aA[6] =  65.39  ; // Zinc
    zZ[7] = 24.0; aA[7] =  51.9961; // Chromium
    zZ[8] = 22.0; aA[8] =  47.867 ; // Titanium
    zZ[9] = 40.0; aA[9] =  91.224 ; // Zirconium

    wW[1] = 0.020500;//Cu - mean value
    wW[2] = 0.000300;//Fe
    wW[3] = 0.022000;//Mn - mean value
    wW[4] = 0.001000;//Mg - mean value
    wW[5] = 0.002000;//Si - mean value
    wW[6] = 0.066500;//Zn
    wW[7] = 0.005000;//Cr
    wW[8] = 0.000600;//Ti
    wW[9] = 0.001650;//Zr - mean value

    totFrac = 0;
    for (Int_t j=1; j<10; j++)
      totFrac += wW[j];
    wW[0] = 1. - totFrac;//Al - the remainder

    den = 2.69;
    AliMixture(34,"HOKOTOL$",aA,zZ,den,+10,wW);
    AliMedium(34,"HOKOTOL$",34,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // Ergal (7075) (yet another Aluminium alloy) - 09 nov 10
    // A,Z from Root TGeoElementTable, W from Web sites
    zZ[0] = 13.0; aA[0] =  26.9815; // Aluminium
    zZ[1] = 29.0; aA[1] =  63.546 ; // Copper
    zZ[2] = 26.0; aA[2] =  55.845 ; // Iron
    zZ[3] = 25.0; aA[3] =  54.938 ; // Manganese
    zZ[4] = 12.0; aA[4] =  24.305 ; // Magnesium
    zZ[5] = 14.0; aA[5] =  28.0855; // Silicon
    zZ[6] = 30.0; aA[6] =  65.39  ; // Zinc
    zZ[7] = 24.0; aA[7] =  51.9961; // Chromium
    zZ[8] = 22.0; aA[8] =  47.867 ; // Titanium

    wW[1] = 0.016000;//Cu - mean value
    wW[2] = 0.005000;//Fe
    wW[3] = 0.003000;//Mn
    wW[4] = 0.025000;//Mg - mean value
    wW[5] = 0.004000;//Si
    wW[6] = 0.056000;//Zn - mean value
    wW[7] = 0.002300;//Cr - mean value
    wW[8] = 0.002000;//Ti

    totFrac = 0;
    for (Int_t j=1; j<9; j++)
      totFrac += wW[j];
    wW[0] = 1. - totFrac;//Al - the remainder

    den = 2.69;
    AliMixture(33,"ERGAL$",aA,zZ,den,+9,wW);
    AliMedium(33,"ERGAL$",33,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    
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

