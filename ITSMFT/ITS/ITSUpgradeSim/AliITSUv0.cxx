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


/* $Id: AliITSUv0.cxx */


//========================================================================
//
//        Geometry for the Upgrade of the Inner Tracking System
//
// Mario Sitta (sitta@to.infn.it)
// Chinorat Kobdaj (kobdaj@g.sut.ac.th)
//
//========================================================================



// $Log: AliITSUv0.cxx,v $

#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoXtru.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TVirtualMC.h>

#include "AliITSU.h"
#include "AliITSUHit.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliTrackReference.h"
#include "AliITSUv0Layer.h"
#include "AliITSUv0.h"
#include "AliITSUGeomTGeo.h"
#include "AliGeomManager.h"
using namespace TMath;


ClassImp(AliITSUv0)

//______________________________________________________________________
AliITSUv0::AliITSUv0()
:  fNWrapVol(0)
  ,fWrapRMin(0)
  ,fWrapRMax(0)
  ,fWrapZSpan(0)
  ,fLay2WrapV(0)
  ,fLayTurbo(0)
  ,fLayPhi0(0)
  ,fLayRadii(0)
  ,fLayZLength(0)
  ,fStavPerLay(0)
  ,fModPerStav(0)
  ,fStaThick(0)
  ,fStaWidth(0)
  ,fStaTilt(0)
  ,fDetThick(0)
  ,fChipTypeID(0)
  ,fBuildLevel(0)
  ,fUpGeom(0)
  ,fStaveModel(kModel0)
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
AliITSUv0::AliITSUv0(const char *title, Int_t nlay)
  :AliITSU(title,nlay)
  ,fNWrapVol(0)
  ,fWrapRMin(0)
  ,fWrapRMax(0)
  ,fWrapZSpan(0)
  ,fLay2WrapV(0)
  ,fLayTurbo(0)
  ,fLayPhi0(0)
  ,fLayRadii(0)
  ,fLayZLength(0)
  ,fStavPerLay(0)
  ,fModPerStav(0)
  ,fStaThick(0)
  ,fStaWidth(0)
  ,fStaTilt(0)
  ,fDetThick(0)
  ,fChipTypeID(0)
  ,fBuildLevel(0)
  ,fUpGeom(0)
  ,fStaveModel(kModel0)
{
  //    Standard constructor for the Upgrade geometry.
  // Inputs:
  //   const char * name   Ignored, set to "ITS"
  //   const char * title  Arbitrary title
  //   const Int_t nlay    Number of layers
  //
  fLayerName = new TString[fNLayers];
  //
  for (Int_t j=0; j<fNLayers; j++)
    fLayerName[j].Form("%s%d",AliITSUGeomTGeo::GetITSSensorPattern(),j); // See AliITSUv0Layer
  //
  fLayTurbo   = new Bool_t[fNLayers];
  fLayPhi0    = new Double_t[fNLayers];
  fLayRadii   = new Double_t[fNLayers];
  fLayZLength = new Double_t[fNLayers];
  fStavPerLay = new Int_t[fNLayers];
  fModPerStav = new Int_t[fNLayers];
  fStaThick   = new Double_t[fNLayers];
  fStaWidth   = new Double_t[fNLayers];
  fStaTilt    = new Double_t[fNLayers];
  fDetThick   = new Double_t[fNLayers];
  fChipTypeID  = new UInt_t[fNLayers];
  fBuildLevel = new Int_t[fNLayers];


  fUpGeom = new AliITSUv0Layer*[fNLayers];
  
  if (fNLayers > 0) { // if not, we'll Fatal-ize in CreateGeometry
    for (Int_t j=0; j<fNLayers; j++) {
      fLayPhi0[j] = 0;
      fLayRadii[j] = 0.;
      fLayZLength[j] = 0.;
      fStavPerLay[j] = 0;
      fModPerStav[j] = 0;
      fStaWidth[j] = 0.;
      fDetThick[j] = 0.;
      fChipTypeID[j] = 0;
      fBuildLevel[j] = 0;
      fUpGeom[j] = 0;     
    }
  }
}

//______________________________________________________________________
AliITSUv0::~AliITSUv0() {
  //    Standard destructor
  // Inputs:
  //   none.
  // Outputs:
  //   none.
  // Return:
  //   none.
  delete [] fLayTurbo;
  delete [] fLayPhi0;
  delete [] fLayRadii;
  delete [] fLayZLength;
  delete [] fStavPerLay;
  delete [] fModPerStav;
  delete [] fStaThick;
  delete [] fStaWidth;
  delete [] fStaTilt;
  delete [] fDetThick;
  delete [] fChipTypeID;
  delete [] fBuildLevel;
  delete [] fUpGeom;
  delete [] fWrapRMin;
  delete [] fWrapRMax;
  delete [] fWrapZSpan;
  delete [] fLay2WrapV;
  //
}

//______________________________________________________________________
void AliITSUv0::AddAlignableVolumes() const{
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

  if (!gGeoManager) { AliFatal("TGeoManager doesn't exist !"); return;  }
  TString pth;
  //
  pth = Form("ALIC_1/%s_2",AliITSUGeomTGeo::GetITSVolPattern());
  // RS: to be checked with MS
  if( !gGeoManager->SetAlignableEntry(AliITSUGeomTGeo::ComposeSymNameITS(),pth.Data()) )
    AliFatal(Form("Unable to set alignable entry ! %s :: %s","ITS",pth.Data()));    
  //
  int chipNum = 0;
  //
  for (int lr=0; lr<fNLayers; lr++) {
    //
    TString wrpV = fLay2WrapV[lr]!=-1 ? Form("%s%d_1/",AliITSUGeomTGeo::GetITSWrapVolPattern(),fLay2WrapV[lr]) : "";
    pth = Form("ALIC_1/%s_2/%s%s%d_1",AliITSUGeomTGeo::GetITSVolPattern(),wrpV.Data(),AliITSUGeomTGeo::GetITSLayerPattern(),lr);
    //printf("SetAlignable: %s %s\n",snm.Data(),pth.Data());
    if( !gGeoManager->SetAlignableEntry(AliITSUGeomTGeo::ComposeSymNameLayer(lr),pth.Data()) ) {
      AliFatal(Form("Unable to set alignable entry ! %s :: %s",AliITSUGeomTGeo::ComposeSymNameLayer(lr),pth.Data()));
    }
    //
    for (int ld=0; ld<fStavPerLay[lr]; ld++) {
      //
      TString pthL = Form("%s/%s%d_%d",pth.Data(),AliITSUGeomTGeo::GetITSStavePattern(),lr,ld);
      //printf("SetAlignable: %s %s\n",snmL.Data(),pthL.Data());
      if ( !gGeoManager->SetAlignableEntry(AliITSUGeomTGeo::ComposeSymNameStave(lr,ld),pthL.Data()) ) {
	AliFatal(Form("Unable to set alignable entry ! %s :: %s",AliITSUGeomTGeo::ComposeSymNameStave(lr,ld),pthL.Data()));
      }
      //
      for (int md=0; md<fModPerStav[lr]; md++) {
	//
	TString pthM = Form("%s/%s%d_%d",pthL.Data(),AliITSUGeomTGeo::GetITSChipPattern(),lr,md);
	//
	// RS: Attention, this is a hack: AliGeomManager cannot accomodate all ITSU chips w/o
	// conflicts with TPC. For this reason we define the UID of the chip to be simply its ID
	//	int modUID = AliGeomManager::LayerToVolUID(lr+1,chipNum++); // here chipNum would be chip within the layer
	int modUID = AliITSUGeomTGeo::ChipVolUID( chipNum++ );
	// 
	if ( !gGeoManager->SetAlignableEntry(AliITSUGeomTGeo::ComposeSymNameChip(lr,ld,-1,-1,md),pthM.Data(),modUID) ) {
	  AliFatal(Form("Unable to set alignable entry ! %s :: %s",AliITSUGeomTGeo::ComposeSymNameChip(lr,ld,-1,-1,md),pthM.Data()));
	}
	//
      }
    }
  }
  //
}

//______________________________________________________________________
void AliITSUv0::SetNWrapVolumes(Int_t n)
{
  // book arrays for wrapper volumes
  if (fNWrapVol) AliFatal(Form("%d wrapper volumes already defined",fNWrapVol));
  if (n<1) return;
  fNWrapVol = n;
  fWrapRMin = new Double_t[fNWrapVol];
  fWrapRMax = new Double_t[fNWrapVol];
  fWrapZSpan= new Double_t[fNWrapVol];
  for (int i=fNWrapVol;i--;) fWrapRMin[i]=fWrapRMax[i]=fWrapZSpan[i]=-1;
  //
}

//______________________________________________________________________
void AliITSUv0::DefineWrapVolume(Int_t id, Double_t rmin,Double_t rmax, Double_t zspan)
{
  // set parameters of id-th wrapper volume
  if (id>=fNWrapVol||id<0) AliFatal(Form("id=%d of wrapper volume is not in 0-%d range",id,fNWrapVol-1));
  fWrapRMin[id] = rmin;
  fWrapRMax[id] = rmax;
  fWrapZSpan[id] = zspan;
}

//______________________________________________________________________
void AliITSUv0::CreateGeometry() {

  // Create the geometry and insert it in the mother volume ITSV
  TGeoManager *geoManager = gGeoManager;

  TGeoVolume *vALIC = geoManager->GetVolume("ALIC");

  new TGeoVolumeAssembly(AliITSUGeomTGeo::GetITSVolPattern());
  TGeoVolume *vITSV = geoManager->GetVolume(AliITSUGeomTGeo::GetITSVolPattern());
  vITSV->SetUniqueID(AliITSUGeomTGeo::GetUIDShift()); // store modID -> midUUID bitshift
  vALIC->AddNode(vITSV, 2, 0);  // Copy number is 2 to cheat AliGeoManager::CheckSymNamesLUT
  //
  const Int_t kLength=100;
  Char_t vstrng[kLength] = "xxxRS"; //?
  vITSV->SetTitle(vstrng);
  //
  // Check that we have all needed parameters
  if (fNLayers <= 0) AliFatal(Form("Wrong number of layers (%d)",fNLayers));
  //
  for (Int_t j=0; j<fNLayers; j++) {
    if (fLayRadii[j] <= 0)                 AliFatal(Form("Wrong layer radius for layer %d (%f)",j,fLayRadii[j]));
    if (fLayZLength[j] <= 0)               AliFatal(Form("Wrong layer length for layer %d (%f)",j,fLayZLength[j]));
    if (fStavPerLay[j] <= 0)               AliFatal(Form("Wrong number of staves for layer %d (%d)",j,fStavPerLay[j]));
    if (fModPerStav[j] <= 0)               AliFatal(Form("Wrong number of chips for layer %d (%d)",j,fModPerStav[j]));
    if (fStaThick[j] < 0)                  AliFatal(Form("Wrong stave thickness for layer %d (%f)",j,fStaThick[j]));
    if (fLayTurbo[j] && fStaWidth[j] <= 0) AliFatal(Form("Wrong stave width for layer %d (%f)",j,fStaWidth[j]));
    if (fDetThick[j] < 0)                  AliFatal(Form("Wrong chip thickness for layer %d (%f)",j,fDetThick[j]));
    //
    if (j > 0) {
      if (fLayRadii[j]<=fLayRadii[j-1])    AliFatal(Form("Layer %d radius (%f) is smaller than layer %d radius (%f)",
							 j,fLayRadii[j],j-1,fLayRadii[j-1]));
    } // if (j > 0)

    if (fStaThick[j] == 0) AliInfo(Form("Stave thickness for layer %d not set, using default",j));
    if (fDetThick[j] == 0) AliInfo(Form("Chip thickness for layer %d not set, using default",j));

  } // for (Int_t j=0; j<fNLayers; j++)

  // Create the wrapper volumes
  TGeoVolume **wrapVols = 0;
  if (fNWrapVol) {
    wrapVols = new TGeoVolume*[fNWrapVol];
    for (int id=0;id<fNWrapVol;id++) {
      wrapVols[id] = CreateWrapperVolume(id);
      vITSV->AddNode(wrapVols[id], 1, 0);
    }
  }
  //
  fLay2WrapV = new Int_t[fNLayers];

  // Now create the actual geometry
  for (Int_t j=0; j<fNLayers; j++) {
    TGeoVolume* dest = vITSV;
    //
    if (fLayTurbo[j]) {
      fUpGeom[j] = new AliITSUv0Layer(j,kTRUE,kFALSE);
      fUpGeom[j]->SetStaveWidth(fStaWidth[j]);
      fUpGeom[j]->SetStaveTilt(fStaTilt[j]);
    }
    else fUpGeom[j] = new AliITSUv0Layer(j,kFALSE);
    //
    fUpGeom[j]->SetPhi0(fLayPhi0[j]);
    fUpGeom[j]->SetRadius(fLayRadii[j]);
    fUpGeom[j]->SetZLength(fLayZLength[j]);
    fUpGeom[j]->SetNStaves(fStavPerLay[j]);
    fUpGeom[j]->SetNChips(fModPerStav[j]);
    fUpGeom[j]->SetChipType(fChipTypeID[j]);
    fUpGeom[j]->SetBuildLevel(fBuildLevel[j]);
    fUpGeom[j]->SetStaveModel(fStaveModel);
    AliDebug(1,Form("fBuildLevel: %d\n",fBuildLevel[j]));
    //
    if (fStaThick[j] != 0) fUpGeom[j]->SetStaveThick(fStaThick[j]);
    if (fDetThick[j] != 0) fUpGeom[j]->SetSensorThick(fDetThick[j]);
    //
    for (int iw=0;iw<fNWrapVol;iw++) {
      if (fLayRadii[j]>fWrapRMin[iw] && fLayRadii[j]<fWrapRMax[iw]) {
	AliInfo(Form("Will embed layer %d in wrapper volume %d",j,iw));
	if (fLayZLength[j]>=fWrapZSpan[iw]) AliFatal(Form("ZSpan %.3f of wrapper volume %d is less than ZSpan %.3f of layer %d",
							  fWrapZSpan[iw],iw,fLayZLength[j],j));
	dest = wrapVols[iw];
	fLay2WrapV[j] = iw;
	break;
      }
    }
    fUpGeom[j]->CreateLayer(dest);
  }
  delete[] wrapVols; // delete pointer only, not the volumes
  //
  ((TGeoVolumeAssembly*)vITSV)->GetShape()->ComputeBBox(); // RS enforce recomputing BBox
  //
}

//______________________________________________________________________
void AliITSUv0::CreateMaterials() {
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

  // Water
  Float_t aWater[2]={1.00794,15.9994};
  Float_t zWater[2]={1.,8.};
  Float_t wWater[2]={0.111894,0.888106};
  Float_t dWater   = 1.0;


  // Kapton
  Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
  Float_t zKapton[4]={1.,6.,7.,8.};
  Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
  Float_t dKapton   = 1.42;
 
  AliMixture(1,"AIR$",aAir,zAir,dAir,4,wAir);
  AliMedium(1, "AIR$",1,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

  AliMixture(2,"WATER$",aWater,zWater,dWater,2,wWater);
  AliMedium(2, "WATER$",2,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(3,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(3,  "SI$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(4,"BERILLIUM$",9.01, 4., 1.848, 35.3, 36.7);// From AliPIPEv3
  AliMedium(4,  "BERILLIUM$",4,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(5,"COPPER$",0.63546E+02,0.29000E+02,0.89600E+01,0.14300E+01,0.99900E+03);
  AliMedium(5,  "COPPER$",5,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    
  // needed for STAVE , Carbon, kapton, Epoxy, flexcable

  //AliMaterial(6,"CARBON$",12.0107,6,2.210,999,999);
  AliMaterial(6,"CARBON$",12.0107,6,2.210/1.3,999,999);
  AliMedium(6,  "CARBON$",6,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMixture(7,"KAPTON(POLYCH2)$", aKapton, zKapton, dKapton, 4, wKapton);
  AliMedium(7, "KAPTON(POLYCH2)$",7,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);


 
  // values below modified as compared to source AliITSv11 !

  //AliMaterial(7,"GLUE$",0.12011E+02,0.60000E+01,0.1930E+01/2.015,999,999); // original
  AliMaterial(7,"GLUE$",12.011,6,1.93/2.015,999,999);  // conform with ATLAS, Corrado, Stefan
  AliMedium(7,  "GLUE$",7,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  // All types of carbon
  // Unidirectional prepreg
 AliMaterial(8,"K13D2U2k$",12.0107,6,1.643,999,999);
 AliMedium(8,  "K13D2U2k$",8,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
 //Impregnated thread
 AliMaterial(9,"M60J3K$",12.0107,6,2.21,999,999);
 AliMedium(9,  "M60J3K$",9,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
 //Impregnated thread
 AliMaterial(10,"M55J6K$",12.0107,6,1.63,999,999);
 AliMedium(10,  "M55J6K$",10,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
 // Fabric(0/90)
 AliMaterial(11,"T300$",12.0107,6,1.725,999,999);
 AliMedium(11,  "T300$",11,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
 //AMEC Thermasol
 AliMaterial(12,"FGS003$",12.0107,6,1.6,999,999);
 AliMedium(12,  "FGS003$",12,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
 // Carbon fleece
 AliMaterial(13,"CarbonFleece$",12.0107,6,0.4,999,999);
 AliMedium(13,  "CarbonFleece$",13,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  // Flex cable
  Float_t aFCm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
  Float_t zFCm[5]={6.,1.,7.,8.,13.};
  Float_t wFCm[5]={0.520088819984,0.01983871336,0.0551367996,0.157399667056, 0.247536};
  //Float_t dFCm = 1.6087;  // original
  //Float_t dFCm = 2.55;   // conform with STAR
   Float_t dFCm = 2.595;   // conform with Corrado

  AliMixture(14,"FLEXCABLE$",aFCm,zFCm,dFCm,5,wFCm);
  AliMedium(14, "FLEXCABLE$",14,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

}

//______________________________________________________________________
void AliITSUv0::DefineLayer(Int_t nlay, double phi0, Double_t r,
			    Double_t zlen, Int_t nstav,
			    Int_t nmod, Double_t lthick,
			    Double_t dthick, UInt_t dettypeID)
{
  //     Sets the layer parameters
  // Inputs:
  //          nlay    layer number
  //          phi0    layer phi0
  //          r       layer radius
  //          zlen    layer length
  //          nstav   number of staves
  //          nmod    number of chips per stave
  //          lthick  stave thickness (if omitted, defaults to 0)
  //          dthick  detector thickness (if omitted, defaults to 0)
  // Outputs:
  //   none.
  // Return:
  //   none.
  
  if (nlay >= fNLayers || nlay < 0) {
    AliError(Form("Wrong layer number (%d)",nlay));
    return;
  }
  
  fLayTurbo[nlay] = kFALSE;
  fLayPhi0[nlay] = phi0;
  fLayRadii[nlay] = r;
  fLayZLength[nlay] = zlen;
  fStavPerLay[nlay] = nstav;
  fModPerStav[nlay] = nmod;
  fStaThick[nlay] = lthick;
  fDetThick[nlay] = dthick;
  fChipTypeID[nlay] = dettypeID;
    
}

//______________________________________________________________________
void AliITSUv0::DefineLayerTurbo(Int_t nlay, Double_t phi0, Double_t r, Double_t zlen, Int_t nstav,
				 Int_t nmod, Double_t width, Double_t tilt,
				 Double_t lthick,Double_t dthick,
				 UInt_t dettypeID, Int_t buildLevel)
{
  //     Sets the layer parameters for a "turbo" layer
  //     (i.e. a layer whose staves overlap in phi)
  // Inputs:
  //          nlay    layer number
  //          phi0    phi of 1st stave
  //          r       layer radius
  //          zlen    layer length
  //          nstav   number of staves
  //          nmod    number of chips per stave
  //          width   stave width
  //          tilt    layer tilt angle (degrees)
  //          lthick  stave thickness (if omitted, defaults to 0)
  //          dthick  detector thickness (if omitted, defaults to 0)
  //          dettypeID  ??
  //          buildLevel (if 0, all geometry is build, used for material budget studies)
  // Outputs:
  //   none.
  // Return:
  //   none.

  if (nlay >= fNLayers || nlay < 0) {
    AliError(Form("Wrong layer number (%d)",nlay));
    return;
  }

  fLayTurbo[nlay] = kTRUE;
  fLayPhi0[nlay] = phi0;
  fLayRadii[nlay] = r;
  fLayZLength[nlay] = zlen;
  fStavPerLay[nlay] = nstav;
  fModPerStav[nlay] = nmod;
  fStaThick[nlay] = lthick;
  fStaWidth[nlay] = width;
  fStaTilt[nlay] = tilt;
  fDetThick[nlay] = dthick;
  fChipTypeID[nlay] = dettypeID;
  fBuildLevel[nlay] = buildLevel;

}

//______________________________________________________________________
void AliITSUv0::GetLayerParameters(Int_t nlay, Double_t &phi0,
				   Double_t &r, Double_t &zlen,
				   Int_t &nstav, Int_t &nmod,
				   Double_t &width, Double_t &tilt,
				   Double_t &lthick, Double_t &dthick,
				   UInt_t &dettype) const
{
  //     Gets the layer parameters
  // Inputs:
  //          nlay    layer number
  // Outputs:
  //          phi0    phi of 1st stave
  //          r       layer radius
  //          zlen    layer length
  //          nstav   number of staves
  //          nmod    number of chips per stave
  //          width   stave width
  //          tilt    stave tilt angle
  //          lthick  stave thickness
  //          dthick  detector thickness
  //          dettype detector type
  // Return:
  //   none.

  if (nlay >= fNLayers || nlay < 0) {
    AliError(Form("Wrong layer number (%d)",nlay));
    return;
  }
  
  phi0   = fLayPhi0[nlay];
  r      = fLayRadii[nlay];
  zlen   = fLayZLength[nlay];
  nstav  = fStavPerLay[nlay];
  nmod   = fModPerStav[nlay];
  width  = fStaWidth[nlay];
  tilt   = fStaTilt[nlay];
  lthick = fStaThick[nlay];
  dthick = fDetThick[nlay];
  dettype= fChipTypeID[nlay];
}

//______________________________________________________________________
TGeoVolume* AliITSUv0::CreateWrapperVolume(Int_t id)
{
  //     Creates an air-filled wrapper cylindrical volume 
  // Inputs:
  //          volume id
  // Outputs:
  //          the wrapper volume

  if (fWrapRMin[id]<0 || fWrapRMax[id]<0 || fWrapZSpan[id]<0) AliFatal(Form("Wrapper volume %d was requested but not defined",id));
  // Now create the actual shape and volume
  //
  TGeoTube *tube = new TGeoTube(fWrapRMin[id], fWrapRMax[id], fWrapZSpan[id]/2.);

  TGeoMedium *medAir = gGeoManager->GetMedium("ITS_AIR$");

  char volnam[30];
  snprintf(volnam, 29, "%s%d", AliITSUGeomTGeo::GetITSWrapVolPattern(),id);

  TGeoVolume *wrapper = new TGeoVolume(volnam, tube, medAir);

  return wrapper;
}

//______________________________________________________________________
void AliITSUv0::Init()
{
  //     Initialise the ITS after it has been created.
  UpdateInternalGeometry();
  AliITSU::Init();
  //  
}

//______________________________________________________________________
Bool_t AliITSUv0::IsLayerTurbo(Int_t nlay)
{
  //     Returns true if the layer is a "turbo" layer
  if ( nlay < 0 || nlay > fNLayers ) {
    AliError(Form("Wrong layer number %d",nlay));
    return kFALSE;
  } 
  else return fUpGeom[nlay]->IsTurbo();
}

//______________________________________________________________________
void AliITSUv0::SetDefaults()
{
  // sets the default segmentation, response, digit and raw cluster classes
}

//______________________________________________________________________
void AliITSUv0::StepManager()
{
  //    Called for every step in the ITS, then calles the AliITSUHit class
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
  if(!(fMC->TrackCharge())) return;
  //
  Int_t copy, lay = 0;
  Int_t id = fMC->CurrentVolID(copy);

  Bool_t notSens = kFALSE;
  while ((lay<fNLayers)  && (notSens = (id!=fIdSens[lay]))) ++lay;
  //printf("R: %.1f | Lay: %d  NotSens: %d\n",positionRS.Pt(), lay, notSens);
	   
  if (notSens) return;

  if(fMC->IsTrackExiting()) {
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kITS);
  } // if Outer ITS mother Volume

  static TLorentzVector position, momentum; // Saves on calls to construtors
  static AliITSUHit hit;// Saves on calls to constructors

  TClonesArray &lhits = *(Hits());
  Int_t   cpn0, cpn1, mod, status = 0;
  //
  // Track status
  if(fMC->IsTrackInside())      status +=  1;
  if(fMC->IsTrackEntering())    status +=  2;
  if(fMC->IsTrackExiting())     status +=  4;
  if(fMC->IsTrackOut())         status +=  8;
  if(fMC->IsTrackDisappeared()) status += 16;
  if(fMC->IsTrackStop())        status += 32;
  if(fMC->IsTrackAlive())       status += 64;

  //
  // retrieve the indices with the volume path
  //
  if (lay < 0 || lay >= fNLayers) {
    AliError(Form("Invalid value: lay=%d. Not an ITS sensitive volume",lay));
    return; // not an ITS sensitive volume.
  } else {
    copy = 1;
    fMC->CurrentVolOffID(1,cpn1);
    fMC->CurrentVolOffID(2,cpn0);
  } //

  mod = fGeomTGeo->GetChipIndex(lay,cpn0,cpn1);
  //RS2DEL  fInitGeom.DecodeDetector(mod,lay+1,cpn0,cpn1,copy);
  //
  // Fill hit structure.
  //
  hit.SetChip(mod);
  hit.SetTrack(gAlice->GetMCApp()->GetCurrentTrackNumber());
  fMC->TrackPosition(position);
  fMC->TrackMomentum(momentum);
  hit.SetPosition(position);
  hit.SetTime(fMC->TrackTime());
  hit.SetMomentum(momentum);
  hit.SetStatus(status);
  hit.SetEdep(fMC->Edep());
  hit.SetShunt(GetIshunt());
  if(fMC->IsTrackEntering()){
    hit.SetStartPosition(position);
    hit.SetStartTime(fMC->TrackTime());
    hit.SetStartStatus(status);
    return; // don't save entering hit.
  } // end if IsEntering
    // Fill hit structure with this new hit.
    //Info("StepManager","Calling Copy Constructor");
  new(lhits[fNhits++]) AliITSUHit(hit); // Use Copy Construtor.
  // Save old position... for next hit.
  hit.SetStartPosition(position);
  hit.SetStartTime(fMC->TrackTime());
  hit.SetStartStatus(status);

  return;
}

//______________________________________________________________________
void AliITSUv0::SetLayerChipTypeID(Int_t lr, UInt_t id)
{
  // set det type
  if (!fChipTypeID || fNLayers<=lr) AliFatal(Form("Number of layers %d, %d is manipulated",fNLayers,lr));
  fChipTypeID[lr] = id;
}

//______________________________________________________________________
Int_t AliITSUv0::GetLayerChipTypeID(Int_t lr)
{
  // set det type
  if (!fChipTypeID || fNLayers<=lr) AliFatal(Form("Number of layers %d, %d is manipulated",fNLayers,lr));
  return fChipTypeID[lr];
}
