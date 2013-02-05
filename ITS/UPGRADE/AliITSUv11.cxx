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


/* $Id: AliITSUv11.cxx */


//========================================================================
//
//        Geometry for the Upgrade of the Inner Tracking System
//
// Mario Sitta (sitta@to.infn.it)
//
//========================================================================



// $Log: AliITSUv11.cxx,v $

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

#include "AliITSU.h"
#include "AliITSUHit.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliTrackReference.h"
#include "AliITSv11Geometry.h"
#include "AliITSUv11Layer.h"
#include "AliITSUv11.h"
#include "AliITSUGeomTGeo.h"
#include "AliGeomManager.h"
using namespace TMath;


ClassImp(AliITSUv11)

//______________________________________________________________________
AliITSUv11::AliITSUv11()
:  fLayTurbo(0)
  ,fLayPhi0(0)
  ,fLayRadii(0)
  ,fLayZLength(0)
  ,fLaddPerLay(0)
  ,fModPerLadd(0)
  ,fLadThick(0)
  ,fLadWidth(0)
  ,fLadTilt(0)
  ,fDetThick(0)
  ,fDetTypeID(0)
  ,fUpGeom(0)
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
AliITSUv11::AliITSUv11(const char *title,const Int_t nlay)
  :AliITSU(title,nlay)
  ,fLayTurbo(0)
  ,fLayPhi0(0)
  ,fLayRadii(0)
  ,fLayZLength(0)
  ,fLaddPerLay(0)
  ,fModPerLadd(0)
  ,fLadThick(0)
  ,fLadWidth(0)
  ,fLadTilt(0)
  ,fDetThick(0)
  ,fDetTypeID(0)
  ,fUpGeom(0)
{
  //    Standard constructor for the Upgrade geometry.
  // Inputs:
  //   const char * name   Ignored, set to "ITS"
  //   const char * title  Arbitrary title
  //   const Int_t nlay    Number of layers
  //
  fLayerName = new TString[fNLayers];
  //
  for (Int_t j=0; j<fNLayers; j++) fLayerName[j].Form("%s%d",AliITSUGeomTGeo::GetITSSensorPattern(),j); // See AliITSUv11Layer
  //
  fLayTurbo   = new Bool_t[fNLayers];
  fLayPhi0    = new Double_t[fNLayers];
  fLayRadii   = new Double_t[fNLayers];
  fLayZLength = new Double_t[fNLayers];
  fLaddPerLay = new Int_t[fNLayers];
  fModPerLadd = new Int_t[fNLayers];
  fLadThick   = new Double_t[fNLayers];
  fLadWidth   = new Double_t[fNLayers];
  fLadTilt    = new Double_t[fNLayers];
  fDetThick   = new Double_t[fNLayers];
  fDetTypeID  = new UInt_t[fNLayers];

  fUpGeom = new AliITSUv11Layer*[fNLayers];
  
  if (fNLayers > 0) { // if not, we'll Fatal-ize in CreateGeometry
    for (Int_t j=0; j<fNLayers; j++) {
      fLayPhi0[j] = 0;
      fLayRadii[j] = 0.;
      fLayZLength[j] = 0.;
      fLaddPerLay[j] = 0;
      fModPerLadd[j] = 0;
      fLadWidth[j] = 0.;
      fDetThick[j] = 0.;
      fDetTypeID[j] = 0;
      fUpGeom[j] = 0;
    }
  }
}

//______________________________________________________________________
AliITSUv11::~AliITSUv11() {
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
  delete [] fLaddPerLay;
  delete [] fModPerLadd;
  delete [] fLadThick;
  delete [] fLadWidth;
  delete [] fLadTilt;
  delete [] fDetThick;
  delete [] fDetTypeID;
  delete [] fUpGeom;
}

//______________________________________________________________________
void AliITSUv11::AddAlignableVolumes() const{
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
  int modNum = 0;
  //
  for (int lr=0; lr<fNLayers; lr++) {
    //
    pth = Form("ALIC_1/%s_2/%s%d_1",AliITSUGeomTGeo::GetITSVolPattern(),AliITSUGeomTGeo::GetITSLayerPattern(),lr);
    //printf("SetAlignable: %s %s\n",snm.Data(),pth.Data());
    gGeoManager->SetAlignableEntry(AliITSUGeomTGeo::ComposeSymNameLayer(lr),pth.Data());
    //
    for (int ld=0; ld<fLaddPerLay[lr]; ld++) {
      //
      TString pthL = Form("%s/%s%d_%d",pth.Data(),AliITSUGeomTGeo::GetITSLadderPattern(),lr,ld);
      //printf("SetAlignable: %s %s\n",snmL.Data(),pthL.Data());
      gGeoManager->SetAlignableEntry(AliITSUGeomTGeo::ComposeSymNameLadder(lr,ld),pthL.Data());
      //
      for (int md=0; md<fModPerLadd[lr]; md++) {
	//
	TString pthM = Form("%s/%s%d_%d",pthL.Data(),AliITSUGeomTGeo::GetITSModulePattern(),lr,md);
	//
	// RS: Attention, this is a hack: AliGeomManager cannot accomodate all ITSU modules w/o
	// conflicts with TPC. For this reason we define the UID of the module to be simply its ID
	//	int modUID = AliGeomManager::LayerToVolUID(lr+1,modNum++); // here modNum would be module within the layer
	int modUID = AliITSUGeomTGeo::ModuleVolUID( modNum++ );
	// 
	gGeoManager->SetAlignableEntry(AliITSUGeomTGeo::ComposeSymNameModule(lr,ld,md),pthM.Data(),modUID);
	//
      }
    }
  }
  //
}

//______________________________________________________________________
void AliITSUv11::CreateGeometry() {

  // Create the geometry and insert it in the mother volume ITSV
  TGeoManager *geoManager = gGeoManager;

  TGeoVolume *vALIC = geoManager->GetVolume("ALIC");

  new TGeoVolumeAssembly(AliITSUGeomTGeo::GetITSVolPattern());
  TGeoVolume *vITSV = geoManager->GetVolume(AliITSUGeomTGeo::GetITSVolPattern());
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
    if (fLaddPerLay[j] <= 0)               AliFatal(Form("Wrong number of ladders for layer %d (%d)",j,fLaddPerLay[j]));
    if (fModPerLadd[j] <= 0)               AliFatal(Form("Wrong number of modules for layer %d (%d)",j,fModPerLadd[j]));
    if (fLadThick[j] < 0)                  AliFatal(Form("Wrong ladder thickness for layer %d (%f)",j,fLadThick[j]));
    if (fLayTurbo[j] && fLadWidth[j] <= 0) AliFatal(Form("Wrong ladder width for layer %d (%f)",j,fLadWidth[j]));
    if (fDetThick[j] < 0)                  AliFatal(Form("Wrong module thickness for layer %d (%f)",j,fDetThick[j]));
    //
    if (j > 0) {
      if (fLayRadii[j]<=fLayRadii[j-1])    AliFatal(Form("Layer %d radius (%f) is smaller than layer %d radius (%f)",
							 j,fLayRadii[j],j-1,fLayRadii[j-1]));
    } // if (j > 0)

    if (fLadThick[j] == 0) AliInfo(Form("Ladder thickness for layer %d not set, using default",j));
    if (fDetThick[j] == 0) AliInfo(Form("Module thickness for layer %d not set, using default",j));

  } // for (Int_t j=0; j<fNLayers; j++)

  // Now create the actual geometry
  for (Int_t j=0; j<fNLayers; j++) {
    if (fLayTurbo[j]) {
      fUpGeom[j] = new AliITSUv11Layer(j,kTRUE,kFALSE);
      fUpGeom[j]->SetLadderWidth(fLadWidth[j]);
      fUpGeom[j]->SetLadderTilt(fLadTilt[j]);
    }
    else fUpGeom[j] = new AliITSUv11Layer(j,kFALSE);
    //
    fUpGeom[j]->SetPhi0(fLayPhi0[j]);
    fUpGeom[j]->SetRadius(fLayRadii[j]);
    fUpGeom[j]->SetZLength(fLayZLength[j]);
    fUpGeom[j]->SetNLadders(fLaddPerLay[j]);
    fUpGeom[j]->SetNModules(fModPerLadd[j]);
    fUpGeom[j]->SetDetType(fDetTypeID[j]);
    //
    if (fLadThick[j] != 0) fUpGeom[j]->SetLadderThick(fLadThick[j]);
    if (fDetThick[j] != 0) fUpGeom[j]->SetSensorThick(fDetThick[j]);
    fUpGeom[j]->CreateLayer(vITSV);
  }
  //
}

//______________________________________________________________________
void AliITSUv11::CreateMaterials() {
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
void AliITSUv11::DefineLayer(const Int_t nlay, const double phi0, const Double_t r,
				 const Double_t zlen, const Int_t nladd,
				 const Int_t nmod, const Double_t lthick,
				 const Double_t dthick, const UInt_t dettypeID)
{
  //     Sets the layer parameters
  // Inputs:
  //          nlay    layer number
  //          phi0    layer phi0
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
  
  if (nlay >= fNLayers || nlay < 0) {
    AliError(Form("Wrong layer number (%d)",nlay));
    return;
  }
  
  fLayTurbo[nlay] = kFALSE;
  fLayPhi0[nlay] = phi0;
  fLayRadii[nlay] = r;
  fLayZLength[nlay] = zlen;
  fLaddPerLay[nlay] = nladd;
  fModPerLadd[nlay] = nmod;
  fLadThick[nlay] = lthick;
  fDetThick[nlay] = dthick;
  fDetTypeID[nlay] = dettypeID;
    
}

//______________________________________________________________________
void AliITSUv11::DefineLayerTurbo(Int_t nlay, Double_t phi0, Double_t r, Double_t zlen, Int_t nladd,
				  Int_t nmod, Double_t width, Double_t tilt,
				  Double_t lthick,Double_t dthick,
				  UInt_t dettypeID)
{
  //     Sets the layer parameters for a "turbo" layer
  //     (i.e. a layer whose ladders overlap in phi)
  // Inputs:
  //          nlay    layer number
  //          phi0    phi of 1st ladder
  //          r       layer radius
  //          zlen    layer length
  //          nladd   number of ladders
  //          nmod    number of modules per ladder
  //          width   ladder width
  //          tilt    layer tilt angle (degrees)
  //          lthick  ladder thickness (if omitted, defaults to 0)
  //          dthick  detector thickness (if omitted, defaults to 0)
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
  fLaddPerLay[nlay] = nladd;
  fModPerLadd[nlay] = nmod;
  fLadThick[nlay] = lthick;
  fLadWidth[nlay] = width;
  fLadTilt[nlay] = tilt;
  fDetThick[nlay] = dthick;
  fDetTypeID[nlay] = dettypeID;
  //
}

//______________________________________________________________________
void AliITSUv11::GetLayerParameters(Int_t nlay, Double_t &phi0,
				    Double_t &r, Double_t &zlen,
				    Int_t &nladd, Int_t &nmod,
				    Double_t &width, Double_t &tilt,
				    Double_t &lthick, Double_t &dthick) const
{
  //     Gets the layer parameters
  // Inputs:
  //          nlay    layer number
  // Outputs:
  //          phi0    phi of 1st ladder
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

  if (nlay >= fNLayers || nlay < 0) {
    AliError(Form("Wrong layer number (%d)",nlay));
    return;
  }
  
  phi0   = fLayPhi0[nlay];
  r      = fLayRadii[nlay];
  zlen   = fLayZLength[nlay];
  nladd  = fLaddPerLay[nlay];
  nmod   = fModPerLadd[nlay];
  width  = fLadWidth[nlay];
  tilt   = fLadTilt[nlay];
  lthick = fLadThick[nlay];
  dthick = fDetThick[nlay];
}

//______________________________________________________________________
void AliITSUv11::Init()
{
  //     Initialise the ITS after it has been created.
  UpdateInternalGeometry();
  AliITSU::Init();
  //  
}

//______________________________________________________________________
Bool_t AliITSUv11::IsLayerTurbo(Int_t nlay)
{
  //     Returns true if the layer is a "turbo" layer
  if ( nlay < 0 || nlay > fNLayers ) {
    AliError(Form("Wrong layer number %d",nlay));
    return kFALSE;
  } 
  else return fUpGeom[nlay]->IsTurbo();
}

//______________________________________________________________________
void AliITSUv11::SetDefaults()
{
  // sets the default segmentation, response, digit and raw cluster classes
}

//______________________________________________________________________
void AliITSUv11::StepManager()
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
  if(!(gMC->TrackCharge())) return;
  //
  Int_t copy, lay = 0;
  Int_t id = gMC->CurrentVolID(copy);

  Bool_t notSens = kFALSE;
  while ((lay<fNLayers)  && (notSens = (id!=fIdSens[lay]))) ++lay;
  //printf("R: %.1f | Lay: %d  NotSens: %d\n",positionRS.Pt(), lay, notSens);
	   
  if (notSens) return;

  if(gMC->IsTrackExiting()) {
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kITS);
  } // if Outer ITS mother Volume

  static TLorentzVector position, momentum; // Saves on calls to construtors
  static AliITSUHit hit;// Saves on calls to constructors

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
  if (lay < 0 || lay >= fNLayers) {
    AliError(Form("Invalid value: lay=%d. Not an ITS sensitive volume",lay));
    return; // not an ITS sensitive volume.
  } else {
    copy = 1;
    gMC->CurrentVolOffID(1,cpn1);
    gMC->CurrentVolOffID(2,cpn0);
  } //

  mod = fGeomTGeo->GetModuleIndex(lay,cpn0,cpn1);
  //RS2DEL  fInitGeom.DecodeDetector(mod,lay+1,cpn0,cpn1,copy);
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
  new(lhits[fNhits++]) AliITSUHit(hit); // Use Copy Construtor.
  // Save old position... for next hit.
  hit.SetStartPosition(position);
  hit.SetStartTime(gMC->TrackTime());
  hit.SetStartStatus(status);

  return;
}

//______________________________________________________________________
void AliITSUv11::SetLayerDetTypeID(Int_t lr, UInt_t id)
{
  // set det type
  if (!fDetTypeID || fNLayers<=lr) AliFatal(Form("Number of layers %d, %d is manipulated",fNLayers,lr));
  fDetTypeID[lr] = id;
}

//______________________________________________________________________
Int_t AliITSUv11::GetLayerDetTypeID(Int_t lr)
{
  // set det type
  if (!fDetTypeID || fNLayers<=lr) AliFatal(Form("Number of layers %d, %d is manipulated",fNLayers,lr));
  return fDetTypeID[lr];
}
