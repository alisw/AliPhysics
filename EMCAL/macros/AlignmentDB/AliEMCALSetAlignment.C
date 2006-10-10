/* $Id$*/

// Script to create alignment parameters and store them into CDB
// Three sets of alignment parameters can be created:
// 1) Ideal geometry
// 2) Geometry with disalignments and disorientations
// 3) Geometry small disalignments and disorientations

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"
#include "TRandom.h"
#include "TClonesArray.h"

#include "AliAlignObjAngles.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALSetAlignment()
{
  TControlBar *menu = new TControlBar("vertical","EMCAL alignment control");
  menu->AddButton("Help to run EMCAL alignment control","Help()",
		  "Explains how to use EMCAL alignment control menus");

  menu->AddButton("Ideal geometry","IdealAlignment()",
		  "Set ideal EMCAL geometry with zero displacement");
  menu->AddButton("Misaligned geometry","FullMisalignment()",
		  "Set EMCAL geometry with large displacement");
  menu->AddButton("Residual misaligned geometry","ResidualAlignment()",
		  "Set EMCAL geometry with small residual displacement");

  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\n\n\nSet EMCAL alignment parameters and write them into ALICE CDB.
Press button \"Ideal geometry\" to create EMCAL geometry with ideal geometry.
Press button \"Misaligned geometry\" to create EMCAL geometry with fully displaced and disorientated geometry.
Press button \"Residual misaligned geometry\" to create EMCAL geometry with infinitesimal displacement and disorientation\n\n\n";
  printf(string);
}

//------------------------------------------------------------------------
void IdealAlignment()
{
  // Create alignment objects for EMCAL with ideally aligned geometry,
  // i.e. with zero displacements and zero disorientations

  // *************************    1st step    ***************
  // Create TClonesArray of alignment objects for EMCAL

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",12);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  // null shifts and rotations

  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 

  TString basePath = "EMCAL/FullSupermodule"; 
  const Int_t nModules=10;

  for (Int_t iModule = 0; iModule < nModules; iModule++) {
    TString newPath = basePath;
    newPath += iModule+1;
    new(alobj[iModule]) AliAlignObjAngles(newPath.Data(),
					    dvoluid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  //half SMs
  new(alobj[10]) AliAlignObjAngles("EMCAL/HalfSupermodule1",
				   dvoluid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[11]) AliAlignObjAngles("EMCAL/HalfSupermodule2",
                                   dvoluid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);


  // *************************    2nd step    ***************
  // Make CDB storage and put TClonesArray in
  // 
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("EMCAL Expert");
  md->SetComment("Alignment objects for ideal geometry, i.e. applying them to TGeo has to leave geometry unchanged");
  AliCDBId id("EMCAL/Align/Data",0,999999);
  CDB->Put(array,id, md);
}

//------------------------------------------------------------------------
void ResidualAlignment()
{
  // Create alignment objects for EMCAL with residual alignment,
  // i.e. with infinitesimal displacement and disorientation

  // *************************    1st step    ***************
  // Create TClonesArray of alignment objects for EMCAL

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",12);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 

  Double_t sigmaTrans = 0.01; Double_t sigmaRot = 0.001;
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  TRandom *pRnd = new TRandom(4357);

  TString basePath = "EMCAL/FullSupermodule";
  for(Int_t iSM = 0; iSM < 10; iSM++) {
    dx = pRnd->Gaus(0,sigmaTrans); dy = pRnd->Gaus(0,sigmaTrans); dz = pRnd->Gaus(0,sigmaTrans);
    dpsi = pRnd->Gaus(0,sigmaRot); dtheta = pRnd->Gaus(0,sigmaRot); dphi = pRnd->Gaus(0,sigmaRot);
    TString newPath = basePath;
    newPath += iSM + 1;
    new(alobj[iSM]) AliAlignObjAngles(newPath.Data(),
				     dvoluid, dx,dy,dz,dpsi,dtheta,dphi, kTRUE);
  }
  
  dx = pRnd->Gaus(0,sigmaTrans); dy = pRnd->Gaus(0,sigmaTrans); dz = pRnd->Gaus(0,sigmaTrans);
  dpsi = pRnd->Gaus(0,sigmaRot); dtheta = pRnd->Gaus(0,sigmaRot); dphi = pRnd->Gaus(0,sigmaRot);
  new(alobj[10]) AliAlignObjAngles("EMCAL/HalfSupermodule1",
				   dvoluid, dx,dy,dz,dpsi,dtheta,dphi, kTRUE);
  
  dx = pRnd->Gaus(0,sigmaTrans); dy = pRnd->Gaus(0,sigmaTrans); dz = pRnd->Gaus(0,sigmaTrans);
  dpsi = pRnd->Gaus(0,sigmaRot); dtheta = pRnd->Gaus(0,sigmaRot); dphi = pRnd->Gaus(0,sigmaRot);
  new(alobj[11]) AliAlignObjAngles("EMCAL/HalfSupermodule2",
				   dvoluid, dx,dy,dz,dpsi,dtheta,dphi, kTRUE);

  // *************************    2nd step    ***************
  // Make CDB storage and put TClonesArray in
  // 
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("EMCAL Expert");
  md->SetComment("Alignment objects for slightly misaligned geometry, i.e. applying them to TGeo has to distirbes geometry very little (resisual misalignment");
  AliCDBId id("EMCAL/Align/Data",1000000,1999999);
  CDB->Put(array,id, md);
}

//------------------------------------------------------------------------
void FullMisalignment()
{
  // Create alignment objects for EMCAL with fully misaligned geometry

  // *************************    1st step    ***************
  // Create TClonesArray of alignment objects for EMCAL

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",12);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 


  Double_t sigmaTrans = 10.; Double_t sigmaRot = 0.1;
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  TRandom *pRnd = new TRandom(4357);

  TString basePath = "EMCAL/FullSupermodule";
  for(Int_t iSM = 0; iSM < 10; iSM++) {
    dx = pRnd->Gaus(0,sigmaTrans); dy = pRnd->Gaus(0,sigmaTrans); dz = pRnd->Gaus(0,sigmaTrans);
    dpsi = pRnd->Gaus(0,sigmaRot); dtheta = pRnd->Gaus(0,sigmaRot); dphi = pRnd->Gaus(0,sigmaRot);
    TString newPath = basePath;
    newPath += iSM + 1;
    new(alobj[iSM]) AliAlignObjAngles(newPath.Data(),
				     dvoluid, dx,dy,dz,dpsi,dtheta,dphi, kTRUE);
  }
  
  dx = pRnd->Gaus(0,sigmaTrans); dy = pRnd->Gaus(0,sigmaTrans); dz = pRnd->Gaus(0,sigmaTrans);
  dpsi = pRnd->Gaus(0,sigmaRot); dtheta = pRnd->Gaus(0,sigmaRot); dphi = pRnd->Gaus(0,sigmaRot);
  new(alobj[10]) AliAlignObjAngles("EMCAL/HalfSupermodule1",
				   dvoluid, dx,dy,dz,dpsi,dtheta,dphi, kTRUE);
  
  dx = pRnd->Gaus(0,sigmaTrans); dy = pRnd->Gaus(0,sigmaTrans); dz = pRnd->Gaus(0,sigmaTrans);
  dpsi = pRnd->Gaus(0,sigmaRot); dtheta = pRnd->Gaus(0,sigmaRot); dphi = pRnd->Gaus(0,sigmaRot);
  new(alobj[11]) AliAlignObjAngles("EMCAL/HalfSupermodule2",
				   dvoluid, dx,dy,dz,dpsi,dtheta,dphi, kTRUE);

  // *************************    2nd step    ***************
  // Make CDB storage and put TClonesArray in
  // 
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("EMCAL Expert");
  md->SetComment("Alignment objects for fully misaligned geometry, i.e. applying them to TGeo has to distirbes geometry very much");
  AliCDBId id("EMCAL/Align/Data",2000000,2999999);
  CDB->Put(array,id, md);
}
