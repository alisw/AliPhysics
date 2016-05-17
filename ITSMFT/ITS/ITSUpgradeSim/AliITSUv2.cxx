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


/* $Id: AliITSUv2.cxx */


//========================================================================
//
//        Geometry for the Upgrade of the Inner Tracking System
//
// Mario Sitta (sitta@to.infn.it)
// Chinorat Kobdaj (kobdaj@g.sut.ac.th)
//
//========================================================================



// $Log: AliITSUv2.cxx,v $

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
#include "AliITSUv2Layer.h"
#include "AliITSUv2.h"
#include "AliITSUGeomTGeo.h"
#include "AliGeomManager.h"
using namespace TMath;


ClassImp(AliITSUv2)

//______________________________________________________________________
AliITSUv2::AliITSUv2()
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
  ,fUnitPerStave(0)
  ,fChipThick(0)
  ,fStaveWidth(0)
  ,fStaveTilt(0)
  ,fSensThick(0)
  ,fChipTypeID(0)
  ,fBuildLevel(0)
  ,fUpGeom(0)
  ,fStaveModelIB(kIBModel0)
  ,fStaveModelOB(kOBModel0)
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
AliITSUv2::AliITSUv2(const char *title, Int_t nlay)
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
  ,fUnitPerStave(0)
  ,fChipThick(0)
  ,fStaveWidth(0)
  ,fStaveTilt(0)
  ,fSensThick(0)
  ,fChipTypeID(0)
  ,fBuildLevel(0)
  ,fUpGeom(0)
  ,fStaveModelIB(kIBModel0)
  ,fStaveModelOB(kOBModel0)
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
    fLayerName[j].Form("%s%d",AliITSUGeomTGeo::GetITSSensorPattern(),j); // See AliITSUv2Layer
  //
  fLayTurbo     = new Bool_t[fNLayers];
  fLayPhi0      = new Double_t[fNLayers];
  fLayRadii     = new Double_t[fNLayers];
  fLayZLength   = new Double_t[fNLayers];
  fStavPerLay   = new Int_t[fNLayers];
  fUnitPerStave = new Int_t[fNLayers];
  fChipThick    = new Double_t[fNLayers];
  fStaveWidth   = new Double_t[fNLayers];
  fStaveTilt    = new Double_t[fNLayers];
  fSensThick    = new Double_t[fNLayers];
  fChipTypeID   = new UInt_t[fNLayers];
  fBuildLevel   = new Int_t[fNLayers];


  fUpGeom = new AliITSUv2Layer*[fNLayers];
  
  if (fNLayers > 0) { // if not, we'll Fatal-ize in CreateGeometry
    for (Int_t j=0; j<fNLayers; j++) {
      fLayPhi0[j]      = 0;
      fLayRadii[j]     = 0.;
      fLayZLength[j]   = 0.;
      fStavPerLay[j]   = 0;
      fUnitPerStave[j] = 0;
      fStaveWidth[j]   = 0.;
      fSensThick[j]    = 0.;
      fChipTypeID[j]   = 0;
      fBuildLevel[j]   = 0;
      fUpGeom[j]       = 0;
    }
  }
}

//______________________________________________________________________
AliITSUv2::~AliITSUv2() {
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
  delete [] fUnitPerStave;
  delete [] fChipThick;
  delete [] fStaveWidth;
  delete [] fStaveTilt;
  delete [] fSensThick;
  delete [] fChipTypeID;
  delete [] fBuildLevel;
  delete [] fUpGeom;
  delete [] fWrapRMin;
  delete [] fWrapRMax;
  delete [] fWrapZSpan;
  delete [] fLay2WrapV;
}

//______________________________________________________________________
void AliITSUv2::AddAlignableVolumes() const
{
  // Creates entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path.
  // 
  // Records in the alignable entries the transformation matrices converting
  // TGeo local coordinates (in the RS of alignable volumes) to the tracking
  // system
  // For this, this function has to run before the misalignment because we
  // are using the ideal positions in the AliITSgeom object.
  AliInfo("Add ITS alignable volumes");

  if (!gGeoManager) { AliFatal("TGeoManager doesn't exist !"); return;  }
  //
  TString path = Form("ALIC_1/%s_2",AliITSUGeomTGeo::GetITSVolPattern());
  TString sname = AliITSUGeomTGeo::ComposeSymNameITS();
  //
  AliDebug(1,Form("%s <-> %s",sname.Data(),path.Data()));
  if( !gGeoManager->SetAlignableEntry(sname.Data(),path.Data()) )
    AliFatal(Form("Unable to set alignable entry ! %s %s",sname.Data(),path.Data()));    
  //
  int lastUID = 0;
  for (int lr=0; lr<fNLayers; lr++) AddAlignableVolumesLayer(lr,path,lastUID);
  //
}

//______________________________________________________________________
void AliITSUv2::AddAlignableVolumesLayer(int lr, TString& parent,Int_t &lastUID) const
{
  // add alignable volumes for layer and its daughters
  //
  TString wrpV = fLay2WrapV[lr]!=-1 ? Form("%s%d_1/",AliITSUGeomTGeo::GetITSWrapVolPattern(),fLay2WrapV[lr]) : "";
  TString path = Form("%s/%s%s%d_1",parent.Data(),wrpV.Data(),AliITSUGeomTGeo::GetITSLayerPattern(),lr);
  TString sname = AliITSUGeomTGeo::ComposeSymNameLayer(lr);
  AliDebug(1,Form("Add %s <-> %s", sname.Data(),path.Data()));
  if ( !gGeoManager->SetAlignableEntry(sname.Data(),path.Data()) ) 
    AliFatal(Form("Unable to set alignable entry ! %s : %s",sname.Data(),path.Data()));
  //
  const AliITSUv2Layer* lrobj = fUpGeom[lr];
  int nstaves = lrobj->GetNStavesPerParent();
  for (int st=0; st<nstaves; st++) AddAlignableVolumesStave(lr,st,path,lastUID);
  //
}
    
//______________________________________________________________________
void AliITSUv2::AddAlignableVolumesStave(int lr, int st, TString& parent, Int_t &lastUID) const
{
  // add alignable volumes for stave and its daughters
  //
  TString path = Form("%s/%s%d_%d",parent.Data(),AliITSUGeomTGeo::GetITSStavePattern(),lr,st);
  TString sname = AliITSUGeomTGeo::ComposeSymNameStave(lr,st);
  AliDebug(1,Form("Add %s <-> %s", sname.Data(),path.Data()));
  if ( !gGeoManager->SetAlignableEntry(sname.Data(),path.Data()) ) 
    AliFatal(Form("Unable to set alignable entry ! %s : %s",sname.Data(),path.Data()));
  //
  const AliITSUv2Layer* lrobj = fUpGeom[lr];
  int nsstave = lrobj->GetNHalfStavesPerParent();
  int start = nsstave>0 ? 0:-1;
  //
  for (int sst=start; sst<nsstave; sst++) AddAlignableVolumesHalfStave(lr,st,sst,path,lastUID);
}

//______________________________________________________________________
void AliITSUv2::AddAlignableVolumesHalfStave(int lr, int st, int sst, TString& parent, Int_t &lastUID) const
{
  // add alignable volumes for halfstave (if any) and its daughters
  //
  TString path = parent;
  if (sst>=0) {
    path = Form("%s/%s%d_%d",parent.Data(),AliITSUGeomTGeo::GetITSHalfStavePattern(),lr,sst);
    TString sname = AliITSUGeomTGeo::ComposeSymNameHalfStave(lr,st,sst);
    AliDebug(1,Form("Add %s <-> %s", sname.Data(),path.Data()));
    if ( !gGeoManager->SetAlignableEntry(sname.Data(),path.Data()) ) 
      AliFatal(Form("Unable to set alignable entry ! %s : %s",sname.Data(),path.Data()));
    //
  }
  const AliITSUv2Layer* lrobj = fUpGeom[lr];
  int nmodules = lrobj->GetNModulesPerParent();
  int start = nmodules>0 ? 0:-1;
  //
  for (int md=start; md<nmodules; md++) AddAlignableVolumesModule(lr,st,sst,md,path,lastUID);
}

//______________________________________________________________________
void AliITSUv2::AddAlignableVolumesModule(int lr, int st, int sst, int md, TString& parent, Int_t &lastUID) const
{
  // add alignable volumes for module (if any) and its daughters
  //
  TString path = parent;
  if (md>=0) {
    path = Form("%s/%s%d_%d",parent.Data(),AliITSUGeomTGeo::GetITSModulePattern(),lr,md);
    TString sname = AliITSUGeomTGeo::ComposeSymNameModule(lr,st,sst,md);
    AliDebug(1,Form("Add %s <-> %s", sname.Data(),path.Data()));
    if ( !gGeoManager->SetAlignableEntry(sname.Data(),path.Data()) ) 
      AliFatal(Form("Unable to set alignable entry ! %s : %s",sname.Data(),path.Data()));
    //
  }
  //
  const AliITSUv2Layer* lrobj = fUpGeom[lr];
  int nchips = lrobj->GetNChipsPerParent();
  //
  for (int ic=0; ic<nchips; ic++) AddAlignableVolumesChip(lr,st,sst,md,ic,path,lastUID);
}  

//______________________________________________________________________
void AliITSUv2::AddAlignableVolumesChip(int lr, int st, int sst, int md, int ch, TString& parent, Int_t &lastUID) const
{
  // add alignable volumes for chip
  //
  TString path = Form("%s/%s%d_%d",parent.Data(),AliITSUGeomTGeo::GetITSChipPattern(),lr,ch);
  TString sname = AliITSUGeomTGeo::ComposeSymNameChip(lr,st,sst,md,ch);
  int modUID = AliITSUGeomTGeo::ChipVolUID( lastUID++ );
  //
  AliDebug(1,Form("Add %s <-> %s : ID=%d", sname.Data(),path.Data(),modUID));
  if ( !gGeoManager->SetAlignableEntry(sname,path.Data(),modUID) )
    AliFatal(Form("Unable to set alignable entry ! %s : %s | %d",sname.Data(),path.Data(),modUID));
  //
}

//______________________________________________________________________
void AliITSUv2::SetNWrapVolumes(Int_t n)
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
void AliITSUv2::DefineWrapVolume(Int_t id, Double_t rmin,Double_t rmax, Double_t zspan)
{
  // set parameters of id-th wrapper volume
  if (id>=fNWrapVol||id<0) AliFatal(Form("id=%d of wrapper volume is not in 0-%d range",id,fNWrapVol-1));
  fWrapRMin[id] = rmin;
  fWrapRMax[id] = rmax;
  fWrapZSpan[id] = zspan;
}

//______________________________________________________________________
void AliITSUv2::CreateGeometry() {

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
    if (fLayRadii[j] <= 0)                   AliFatal(Form("Wrong layer radius for layer %d (%f)",j,fLayRadii[j]));
    if (fLayZLength[j] <= 0)                 AliFatal(Form("Wrong layer length for layer %d (%f)",j,fLayZLength[j]));
    if (fStavPerLay[j] <= 0)                 AliFatal(Form("Wrong number of staves for layer %d (%d)",j,fStavPerLay[j]));
    if (fUnitPerStave[j] <= 0)               AliFatal(Form("Wrong number of chips for layer %d (%d)",j,fUnitPerStave[j]));
    if (fChipThick[j] < 0)                   AliFatal(Form("Wrong chip thickness for layer %d (%f)",j,fChipThick[j]));
    if (fLayTurbo[j] && fStaveWidth[j] <= 0) AliFatal(Form("Wrong stave width for layer %d (%f)",j,fStaveWidth[j]));
    if (fSensThick[j] < 0)                   AliFatal(Form("Wrong sensor thickness for layer %d (%f)",j,fSensThick[j]));
    //
    if (j > 0) {
      if (fLayRadii[j]<=fLayRadii[j-1])      AliFatal(Form("Layer %d radius (%f) is smaller than layer %d radius (%f)",
			   j,fLayRadii[j],j-1,fLayRadii[j-1]));
    } // if (j > 0)

    if (fChipThick[j] == 0) AliInfo(Form("Chip thickness for layer %d not set, using default",j));
    if (fSensThick[j] == 0) AliInfo(Form("Sensor thickness for layer %d not set, using default",j));

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
    fLay2WrapV[j] = -1;
    //
    if (fLayTurbo[j]) {
      fUpGeom[j] = new AliITSUv2Layer(j,kTRUE,kFALSE);
      fUpGeom[j]->SetStaveWidth(fStaveWidth[j]);
      fUpGeom[j]->SetStaveTilt(fStaveTilt[j]);
    }
    else fUpGeom[j] = new AliITSUv2Layer(j,kFALSE);
    //
    fUpGeom[j]->SetPhi0(fLayPhi0[j]);
    fUpGeom[j]->SetRadius(fLayRadii[j]);
    fUpGeom[j]->SetZLength(fLayZLength[j]);
    fUpGeom[j]->SetNStaves(fStavPerLay[j]);
    fUpGeom[j]->SetNUnits(fUnitPerStave[j]);
    fUpGeom[j]->SetChipType(fChipTypeID[j]);
    fUpGeom[j]->SetBuildLevel(fBuildLevel[j]);
    if (j < 3)
      fUpGeom[j]->SetStaveModel(fStaveModelIB);
    else
      fUpGeom[j]->SetStaveModel(fStaveModelOB);
    AliDebug(1,Form("fBuildLevel: %d\n",fBuildLevel[j]));
    //
    if (fChipThick[j] != 0) fUpGeom[j]->SetChipThick(fChipThick[j]);
    if (fSensThick[j] != 0) fUpGeom[j]->SetSensorThick(fSensThick[j]);
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
	fUpGeom[j]->CreateBarrelLayer(dest);	// #pnamwong
    
  }
  // CreateSuppCyl(kTRUE,wrapVols[0]);		// Disabled by pnamwong
  // CreateSuppCyl(kFALSE,wrapVols[2]);		// Disabled by pnamwong

  delete[] wrapVols; // delete pointer only, not the volumes
  //
  ((TGeoVolumeAssembly*)vITSV)->GetShape()->ComputeBBox(); // RS enforce recomputing BBox

  
  //-------------------CREATE Wrap Volume for Barrel--------------------
  // #pnamwong
	//TGeoVolume* barrVol = new TGeoVolumeAssembly(AliITSUGeomTGeo::GetITSWrapVolPattern() + TString("Barrel"));
	//AliITSUv2Layer::CreateBarrel(barrVol);
	//vALIC->AddNode(barrVol, 1, 0);

	//vALIC->AddNode(new TGeoVolumeAssembly("ITSVBarrelV3"), 1, 0);

  //--------------------------------------------------------------------

}

//____________________________________________________________
//Service Barrel
void AliITSUv2::CreateSuppCyl(const Bool_t innerBarrel,TGeoVolume *dest,const TGeoManager *mgr){
  // Creates the Service Barrel (as a simple cylinder) for IB and OB
  // Inputs:
  //         innerBarrel : if true, build IB service barrel, otherwise for OB
  //         dest        : the mother volume holding the service barrel
  //         mgr         : the gGeoManager pointer (used to get the material)
  //

  Double_t rminIB =  4.7;
  Double_t rminOB = 43.9;
  Double_t zLenOB ;
  Double_t cInt	= 0.22; //dimensioni cilindro di supporto interno
  Double_t cExt	= 1.00; //dimensioni cilindro di supporto esterno
//  Double_t phi1   =  180;
//  Double_t phi2   =  360;


  TGeoMedium *medCarbonFleece = mgr->GetMedium("ITS_CarbonFleece$");

  if (innerBarrel){
    zLenOB=((TGeoTube*)(dest->GetShape()))->GetDz();
//    TGeoTube*ibSuppSh = new TGeoTubeSeg(rminIB,rminIB+cInt,zLenOB,phi1,phi2);
    TGeoTube*ibSuppSh = new TGeoTube(rminIB,rminIB+cInt,zLenOB);
    TGeoVolume *ibSupp = new TGeoVolume("ibSuppCyl",ibSuppSh,medCarbonFleece);
    dest->AddNode(ibSupp,1);
  }
  else {
    zLenOB=((TGeoTube*)(dest->GetShape()))->GetDz();
    TGeoTube*obSuppSh=new TGeoTube(rminOB,rminOB+cExt,zLenOB);
    TGeoVolume *obSupp=new TGeoVolume("obSuppCyl",obSuppSh,medCarbonFleece);
    dest->AddNode(obSupp,1);
  }

  return;
}

//______________________________________________________________________
void AliITSUv2::CreateMaterials() {
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

  // PEEK CF30
  Float_t aPEEK[3]={12.0107,1.00794,15.9994};
  Float_t zPEEK[3]={6.,1.,8.};
  Float_t wPEEK[3]={19.,12.,3};
  Float_t dPEEK   = 1.32;

  // Kapton
  Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
  Float_t zKapton[4]={1.,6.,7.,8.};
  Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
  Float_t dKapton   = 1.42;

  // Tungsten Carbide
  Float_t aWC[2]={183.84, 12.0107};
  Float_t zWC[2]={74, 6};
  Float_t wWC[2]={0.5, 0.5};
  Float_t dWC   = 15.63;

  // Inox 304
  Float_t aInox304[4]={12.0107,51.9961,58.6928,55.845};
  Float_t zInox304[4]={6.,24.,28,26}; // C, Cr, Ni, Fe
  Float_t wInox304[4]={0.0003,0.18,0.10,0}; // [3] will be computed
  Float_t dInox304   = 7.85;

 
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
  AliMaterial(15,"GLUE$",12.011,6,1.93/2.015,999,999);  // conform with ATLAS, Corrado, Stefan
  AliMedium(15,  "GLUE$",15,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  // All types of carbon
  // Unidirectional prepreg
 AliMaterial(8,"K13D2U2k$",12.0107,6,1.643,999,999);
 AliMedium(8,  "K13D2U2k$",8,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
 AliMaterial(17,"K13D2U120$",12.0107,6,1.583,999,999);
 AliMedium(17,  "K13D2U120$",17,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
  // Carbon prepreg woven
 AliMaterial(18,"F6151B05M$",12.0107,6,2.133,999,999);
 AliMedium(18,  "F6151B05M$",18,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
 //Impregnated thread
 AliMaterial(9,"M60J3K$",12.0107,6,2.12,999,999);
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

 // PEEK CF30
 AliMixture(19,"PEEKCF30$",aPEEK,zPEEK,dPEEK,-3,wPEEK);
 AliMedium(19, "PEEKCF30$",19,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  // Flex cable
  Float_t aFCm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
  Float_t zFCm[5]={6.,1.,7.,8.,13.};
  Float_t wFCm[5]={0.520088819984,0.01983871336,0.0551367996,0.157399667056, 0.247536};
  //Float_t dFCm = 1.6087;  // original
  //Float_t dFCm = 2.55;   // conform with STAR
   Float_t dFCm = 2.595;   // conform with Corrado

  AliMixture(14,"FLEXCABLE$",aFCm,zFCm,dFCm,5,wFCm);
  AliMedium(14, "FLEXCABLE$",14,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(16,"ALUMINUM$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
  AliMedium(16,  "ALUMINUM$",16,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMixture(20,"TUNGCARB$",aWC,zWC,dWC,2,wWC);
  AliMedium(20, "TUNGCARB$",20,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  wInox304[3] = 1. - wInox304[0] - wInox304[1] - wInox304[2];
  AliMixture(21,"INOX304$",aInox304,zInox304,dInox304,4,wInox304);
  AliMedium(21, "INOX304$",21,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

	// --------------------------- 
	// needed for ITS IB Endwheels & DC-DC units.
	// #pnamwong
	
	// FR4 for PCB of DCDC, copied from AliITSv11
    Float_t zG10FR4[14] = {14.00,	20.00,	13.00,	12.00,	5.00,	22.00,	11.00,	19.00,	26.00,	9.00,	8.00,	6.00,	7.00,	1.00};
    Float_t aG10FR4[14] = {28.0855000,40.0780000,26.9815380,24.3050000,10.8110000,47.8670000,22.9897700,39.0983000,55.8450000,18.9984000,15.9994000,12.0107000,14.0067000,1.0079400};
    Float_t wG10FR4[14] = {0.15144894,0.08147477,0.04128158,0.00904554,0.01397570,0.00287685,0.00445114,0.00498089,0.00209828,0.00420000,0.36043788,0.27529426,0.01415852,0.03427566};
    Float_t densG10FR4= 1.8;
    AliMixture(22,"G10FR4$",aG10FR4,zG10FR4,densG10FR4,14,wG10FR4);
    AliMedium(22,"G10FR4$",22,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    
    // Polyeurethane (C3H8N20 for cooling pipe of DCDC, information from "www.chemnet.com"
    Float_t aPUR[4] = { 12.0107, 1.00794, 14.00674, 15.9994};
    Float_t zPUR[4] = { 6.     , 1.     ,  7.     ,  8.    };
    Float_t wPUR[4] = { 3.     , 8.     ,  2.     ,  1.    };
    Float_t dPUR    = 1.005;
    AliMixture(23,"PUR$",aPUR,zPUR,dPUR,-4,wPUR);
    AliMedium(23,"PUR$",23,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

	// FR4+CU for PCB of DCDC(DCDC-PCBCU), modified from G10FR4
    Float_t zG10FR4CU[15] = {14.00,	20.00,	13.00,	12.00,	5.00,	22.00,	11.00,	19.00,	26.00,	9.00,	8.00,	6.00,	7.00,	1.00, 29.00};
    Float_t aG10FR4CU[15] = {28.0855000,40.0780000,26.9815380,24.3050000,10.8110000,47.8670000,22.9897700,39.0983000,55.8450000,18.9984000,15.9994000,12.0107000,14.0067000,1.0079400,63.5460000};
    Float_t wG10FR4CU[15] = {0.137785237,0.074124127,0.037557161,0.008229453,0.012714814,0.002617301,0.004049559,0.004531515,0.001908973,0.003821077,0.327919224,0.250457249,0.01288114,0.031183315,0.090219864};
    Float_t densG10FR4CU= 2.527172093;
    AliMixture(24,"G10FR4CU$",aG10FR4CU,zG10FR4CU,densG10FR4CU,15,wG10FR4CU);
    AliMedium(24,"G10FR4CU$",24,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

	// Coil+Shield+Passive60%+Air for SHIELD of DCDC (DCDC-SHIELD).
	// Still not include Gap Pad (Bergquist GP3000S30).
    Float_t zDCDCSHIELD[20] = {29.00,		6.00,		 1.00,		  29.00,	   6.00,		1.00,		 29.00,		  13.00,	   8.00,		28.00,	  	 50.00,	  	  56.00,	   22.00,	    8.00,		 28.00,	  	  50.00,	   6.00,		7.00,		 8.00,		  18.00};
    Float_t aDCDCSHIELD[20] = {63.54600000, 12.01070000, 1.007940000, 63.54600000, 12.01070000, 1.007940000, 63.54600000, 26.98200000, 15.99940000, 58.69340000, 118.7100000, 137.3270000, 47.86700000, 15.99940000, 58.69340000, 118.7100000, 12.01070000, 14.00670000, 15.99940000, 39.94800000};
    Float_t wDCDCSHIELD[20] = {0.011184966, 0.160966032, 0.027016594, 0.013410143, 0.081330205, 0.013650490, 0.000891892, 0.000098852, 0.000087924, 0.000017375, 0.000013031, 0.009534679, 0.003323429, 0.003332538, 0.001506107, 0.001129580, 0.000083418, 0.508085744, 0.155924490, 0.008629022};
    Float_t densDCDCSHIELD= 0.643338921;
    AliMixture(25,"DCDCSHIELD$",aDCDCSHIELD,zDCDCSHIELD,densDCDCSHIELD,20,wDCDCSHIELD);
    AliMedium(25,"DCDCSHIELD$",25,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

	// Passive40%+Connector for SHIELD of DCDC (DCDC-PASSCNT).
    Float_t zDCDCPASSCNT[18] = {14.00,		29.00,		 13.00,		  8.00,	   	   28.00,		50.00,		 56.00,		  22.00,	   8.00,		28.00,	  	 50.00,	  	  6.00,	   	   1.00,	    29.00,		 6.00,		  7.00,		   8.00,		18.00};
    Float_t aDCDCPASSCNT[18] = {28.08600000, 63.54600000, 26.98200000, 15.99940000, 58.69340000, 118.7100000, 137.3270000, 47.86700000, 15.99940000, 58.69340000, 118.7100000, 12.01070000, 1.007940000, 63.54600000, 12.01070000, 14.00670000, 15.99940000, 39.94800000};
    Float_t wDCDCPASSCNT[18] = {0.005062500, 0.007129630, 0.000316083, 0.000281139, 0.000055556, 0.000041667, 0.030487431, 0.010626766, 0.010655893, 0.004815822, 0.003611867, 0.051287664, 0.008608169, 0.023958333, 0.000104540, 0.636736516, 0.195405633, 0.010813950};
    Float_t densDCDCPASSCNT= 0.759649185;
    AliMixture(26,"DCDCPASSCNT$",aDCDCPASSCNT,zDCDCPASSCNT,densDCDCPASSCNT,18,wDCDCPASSCNT);
    AliMedium(26,"DCDCPASSCNT$",26,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

	// Carbon Foam
	AliMaterial(27,"CarbonFoam$",12.0107,6,0.01,999,999);
	AliMedium(27,  "CarbonFoam$",27,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

	// --------------------------- 
}

//______________________________________________________________________
void AliITSUv2::DefineLayer(Int_t nlay, double phi0, Double_t r,
			    Double_t zlen, Int_t nstav,
			    Int_t nunit, Double_t lthick,
			    Double_t dthick, UInt_t dettypeID,
			    Int_t buildLevel)
{
  //     Sets the layer parameters
  // Inputs:
  //          nlay    layer number
  //          phi0    layer phi0
  //          r       layer radius
  //          zlen    layer length
  //          nstav   number of staves
  //          nunit   IB: number of chips per stave
  //                  OB: number of modules per half stave
  //          lthick  chip thickness (if omitted, defaults to 0)
  //          dthick  sensor thickness (if omitted, defaults to 0)
  //          dettypeID  ??
  //          buildLevel (if 0, all geometry is build, used for material budget studies)
  // Outputs:
  //   none.
  // Return:
  //   none.

  AliInfo(Form("L# %d Phi:%+.3f R:%+7.3f DZ:%7.2f Nst:%2d Nunit:%2d Lthick:%.4f Dthick:%.4f DetID:%d B:%d",
	       nlay,phi0,r,zlen,nstav,nunit,lthick,dthick,dettypeID,buildLevel));

  if (nlay >= fNLayers || nlay < 0) {
    AliError(Form("Wrong layer number (%d)",nlay));
    return;
  }
  
  fLayTurbo[nlay] = kFALSE;
  fLayPhi0[nlay] = phi0;
  fLayRadii[nlay] = r;
  fLayZLength[nlay] = zlen;
  fStavPerLay[nlay] = nstav;
  fUnitPerStave[nlay] = nunit;
  fChipThick[nlay] = lthick;
  fSensThick[nlay] = dthick;
  fChipTypeID[nlay] = dettypeID;
  fBuildLevel[nlay] = buildLevel;
    
}

//______________________________________________________________________
void AliITSUv2::DefineLayerTurbo(Int_t nlay, Double_t phi0, Double_t r, Double_t zlen, Int_t nstav,
				 Int_t nunit, Double_t width, Double_t tilt,
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
  //          nunit   IB: number of chips per stave
  //                  OB: number of modules per half stave
  //          width   stave width
  //          tilt    layer tilt angle (degrees)
  //          lthick  chip thickness (if omitted, defaults to 0)
  //          dthick  sensor thickness (if omitted, defaults to 0)
  //          dettypeID  ??
  //          buildLevel (if 0, all geometry is build, used for material budget studies)
  // Outputs:
  //   none.
  // Return:
  //   none.

  AliInfo(Form("L# %d Phi:%+.3f R:%+7.3f DZ:%7.2f Nst:%2d Nunit:%2d W:%7.4f Tilt:%+.3f Lthick:%.4f Dthick:%.4f DetID:%d B:%d",
	       nlay,phi0,r,zlen,nstav,nunit,width,tilt,lthick,dthick,dettypeID,buildLevel));

  if (nlay >= fNLayers || nlay < 0) {
    AliError(Form("Wrong layer number (%d)",nlay));
    return;
  }

  fLayTurbo[nlay] = kTRUE;
  fLayPhi0[nlay] = phi0;
  fLayRadii[nlay] = r;
  fLayZLength[nlay] = zlen;
  fStavPerLay[nlay] = nstav;
  fUnitPerStave[nlay] = nunit;
  fChipThick[nlay] = lthick;
  fStaveWidth[nlay] = width;
  fStaveTilt[nlay] = tilt;
  fSensThick[nlay] = dthick;
  fChipTypeID[nlay] = dettypeID;
  fBuildLevel[nlay] = buildLevel;

}

//______________________________________________________________________
void AliITSUv2::GetLayerParameters(Int_t nlay, Double_t &phi0,
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
  //          nmod    IB: number of chips per stave
  //                  OB: number of modules per half stave
  //          width   stave width
  //          tilt    stave tilt angle
  //          lthick  chip thickness
  //          dthick  sensor thickness
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
  nmod   = fUnitPerStave[nlay];
  width  = fStaveWidth[nlay];
  tilt   = fStaveTilt[nlay];
  lthick = fChipThick[nlay];
  dthick = fSensThick[nlay];
  dettype= fChipTypeID[nlay];
}

//______________________________________________________________________
TGeoVolume* AliITSUv2::CreateWrapperVolume(Int_t id)
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
void AliITSUv2::Init()
{
  //     Initialise the ITS after it has been created.
  UpdateInternalGeometry();
  AliITSU::Init();
  //  
}

//______________________________________________________________________
Bool_t AliITSUv2::IsLayerTurbo(Int_t nlay)
{
  //     Returns true if the layer is a "turbo" layer
  if ( nlay < 0 || nlay > fNLayers ) {
    AliError(Form("Wrong layer number %d",nlay));
    return kFALSE;
  } 
  else return fUpGeom[nlay]->IsTurbo();
}

//______________________________________________________________________
void AliITSUv2::SetDefaults()
{
  // sets the default segmentation, response, digit and raw cluster classes
}

//______________________________________________________________________
void AliITSUv2::StepManager()
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
  //
  if(!(this->IsActive())) return;
  if(!(fMC->TrackCharge())) return;
  //
  Int_t copy, lay = 0;
  Int_t id = fMC->CurrentVolID(copy);

  Bool_t notSens = kFALSE;
  while ((lay<fNLayers)  && (notSens = (id!=fIdSens[lay]))) ++lay;
  //printf("R: %.1f | Lay: %d  NotSens: %d\n",positionRS.Pt(), lay, notSens);
	   
  if (notSens) return;
  //
  if (lay < 0 || lay >= fNLayers) {
    AliError(Form("Invalid value: lay=%d. Not an ITS sensitive volume",lay));
    return; // not an ITS sensitive volume.
  } 
  //
  static TLorentzVector position, momentum; // Saves on calls to construtors
  static AliITSUHit hit;// Saves on calls to constructors
  //
  TClonesArray &lhits = *(Hits());
  Int_t chipID, status = 0;
  //
  // Track status
  if(fMC->IsTrackInside())      status +=  1;
  if(fMC->IsTrackEntering())    status +=  2;
  if(fMC->IsTrackExiting()) {
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kITS);
    status +=  4;
  } // if Outer ITS mother Volume
  if(fMC->IsTrackOut())         status +=  8;
  if(fMC->IsTrackDisappeared()) status += 16;
  if(fMC->IsTrackStop())        status += 32;
  if(fMC->IsTrackAlive())       status += 64;
  //
  // retrieve the indices with the volume path
  //
  fMC->TrackPosition(position);
  int chip=-1,module=-1,sstave=-1,stave=-1,level=0; // volume copies on different levels
  fMC->CurrentVolOffID(++level,chip);
  if (fGeomTGeo->GetNModules(lay)>0)    fMC->CurrentVolOffID(++level,module);
  if (fGeomTGeo->GetNHalfStaves(lay)>0) fMC->CurrentVolOffID(++level,sstave);
  fMC->CurrentVolOffID(++level,stave);
  //
  chipID = fGeomTGeo->GetChipIndex(lay,stave,sstave,module,chip);
  // Fill hit structure.
  //
  hit.SetChip(chipID);
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
void AliITSUv2::SetLayerChipTypeID(Int_t lr, UInt_t id)
{
  // set det type
  if (!fChipTypeID || fNLayers<=lr) AliFatal(Form("Number of layers %d, %d is manipulated",fNLayers,lr));
  fChipTypeID[lr] = id;
}

//______________________________________________________________________
Int_t AliITSUv2::GetLayerChipTypeID(Int_t lr)
{
  // set det type
  if (!fChipTypeID || fNLayers<=lr) AliFatal(Form("Number of layers %d, %d is manipulated",fNLayers,lr));
  return fChipTypeID[lr];
}
