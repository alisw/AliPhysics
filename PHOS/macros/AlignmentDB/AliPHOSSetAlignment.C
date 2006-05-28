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


void AliPHOSSetAlignment()
{
  TControlBar *menu = new TControlBar("vertical","PHOS alignment control");
  menu->AddButton("Help to run PHOS alignment control","Help()",
		  "Explains how to use PHOS alignment control menus");

  menu->AddButton("Ideal geometry","IdealAlignment()",
		  "Set ideal PHOS geometry with zero displacement");
  menu->AddButton("Misaligned geometry","FullMisalignment()",
		  "Set PHOS geometry with large displacement");
  menu->AddButton("Residual misaligned geometry","ResidualAlignment()",
		  "Set PHOS geometry with small residual displacement");

  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\n\n\nSet PHOS alignment parameters and write them into ALICE CDB.
Press button \"Ideal geometry\" to create PHOS geometry with ideal geometry.
Press button \"Misaligned geometry\" to create PHOS geometry with fully displaced and disorientated geometry.
Press button \"Residual misaligned geometry\" to create PHOS geometry with infinitesimal displacement and disorientation\n\n\n";
  printf(string);
}

//------------------------------------------------------------------------
void IdealAlignment()
{
  // Create alignment objects for PHOS with ideally aligned geometry,
  // i.e. with zero displacements and zero disorientations

  // *************************    1st step    ***************
  // Create TClonesArray of alignment objects for PHOS

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",11);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  // null shifts and rotations

  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 

  TString basePath = "/ALIC_1/PHOS_"; 
  const Int_t nModules=5;

  for (Int_t iModule = 1; iModule<=nModules; iModule++) {
    TString newPath = basePath;
    newPath += iModule;
    new(alobj[iModule-1]) AliAlignObjAngles(newPath.Data(),
					    dvoluid, dx, dy, dz, dpsi, dtheta, dphi);
  }

  // *************************    2nd step    ***************
  // Make CDB storage and put TClonesArray in
  // 
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("Alignment objects for ideal geometry, i.e. applying them to TGeo has to leave geometry unchanged");
  AliCDBId id("PHOS/Align/Data",0,0);
  CDB->Put(array,id, md);
}

//------------------------------------------------------------------------
void ResidualAlignment()
{
  // Create alignment objects for PHOS with residual alignment,
  // i.e. with infinitesimal displacement and disorientation

  // *************************    1st step    ***************
  // Create TClonesArray of alignment objects for PHOS

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",11);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  Double_t dpsi=0., dtheta=0., dphi=0.;
  Double_t displacement = 0.2;

  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 

  // Alignment for 5 PHOS modules
  new(alobj[0]) AliAlignObjAngles("/ALIC_1/PHOS_1",
				  dvoluid, -0.20, -0.1, +0.0, dpsi, dtheta, 0.2);
  new(alobj[1]) AliAlignObjAngles("/ALIC_1/PHOS_2",
				  dvoluid, -0.10, +0.0, -0.2, dpsi, dtheta, 0.2);
  new(alobj[2]) AliAlignObjAngles("/ALIC_1/PHOS_3",
				  dvoluid,  0.05, -0.1,  0.2, dpsi, dtheta, 0.0);
  new(alobj[3]) AliAlignObjAngles("/ALIC_1/PHOS_4",
				  dvoluid, +0.10, -0.0, -0.1, dpsi, dtheta, 0.1);
  new(alobj[4]) AliAlignObjAngles("/ALIC_1/PHOS_5",
				  dvoluid, +0.20, -0.1,  0.1, dpsi, dtheta, 0.2);

  // Alignment for PHOS cradle
  new(alobj[5]) AliAlignObjAngles("/ALIC_1/PCRA_0",
				  dvoluid, 0., 0., -displacement, dpsi, dtheta, dphi);
  new(alobj[6]) AliAlignObjAngles("/ALIC_1/PCRA_1",
				  dvoluid, 0., 0., +displacement, dpsi, dtheta, dphi);

  // Alignment for cradle wheels
  new(alobj[7])  AliAlignObjAngles("/ALIC_1/PWHE_0",
				   dvoluid, 0., 0., -displacement, dpsi, dtheta, dphi);
  new(alobj[8])  AliAlignObjAngles("/ALIC_1/PWHE_1",
				   dvoluid, 0., 0., -displacement, dpsi, dtheta, dphi);
  new(alobj[9])  AliAlignObjAngles("/ALIC_1/PWHE_2",
				   dvoluid, 0., 0., +displacement, dpsi, dtheta, dphi);
  new(alobj[10]) AliAlignObjAngles("/ALIC_1/PWHE_3",
				   dvoluid, 0., 0., +displacement, dpsi, dtheta, dphi);

  // *************************    2nd step    ***************
  // Make CDB storage and put TClonesArray in
  // 
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("Alignment objects for slightly misaligned geometry, i.e. applying them to TGeo has to distirbes geometry very little (resisual misalignment");
  AliCDBId id("PHOS/Align/Data",200,200);
  CDB->Put(array,id, md);
}

//------------------------------------------------------------------------
void FullMisalignment()
{
  // Create alignment objects for PHOS with fully misaligned geometry

  // *************************    1st step    ***************
  // Create TClonesArray of alignment objects for PHOS

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",11);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  Double_t dpsi=0., dtheta=0., dphi=0.;
  Double_t displacement = 10;

  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 

  // Alignment for 5 PHOS modules
  new(alobj[0]) AliAlignObjAngles("/ALIC_1/PHOS_1",
				  dvoluid, -20., -10.,   0., dpsi, dtheta, 5);
  new(alobj[1]) AliAlignObjAngles("/ALIC_1/PHOS_2",
				  dvoluid, -10.,   0., -10., dpsi, dtheta, 2);
  new(alobj[2]) AliAlignObjAngles("/ALIC_1/PHOS_3",
				  dvoluid,   5., -10.,  10., dpsi, dtheta, 0);
  new(alobj[3]) AliAlignObjAngles("/ALIC_1/PHOS_4",
				  dvoluid, +10.,  -0., -10., dpsi, dtheta, 2);
  new(alobj[4]) AliAlignObjAngles("/ALIC_1/PHOS_5",
				  dvoluid, +20., -10.,   0., dpsi, dtheta, 5);

  // Alignment for PHOS cradle
  new(alobj[5]) AliAlignObjAngles("/ALIC_1/PCRA_0",
				  dvoluid, 0., 0., -displacement, dpsi, dtheta, dphi);
  new(alobj[6]) AliAlignObjAngles("/ALIC_1/PCRA_1",
				  dvoluid, 0., 0., +displacement, dpsi, dtheta, dphi);

  // Alignment for cradle wheels
  new(alobj[7]) AliAlignObjAngles("/ALIC_1/PWHE_0",
				  dvoluid, 0., 0., -displacement, dpsi, dtheta, dphi);
  new(alobj[8]) AliAlignObjAngles("/ALIC_1/PWHE_1",
				  dvoluid, 0., 0., -displacement, dpsi, dtheta, dphi);
  new(alobj[9]) AliAlignObjAngles("/ALIC_1/PWHE_2",
				  dvoluid, 0., 0., +displacement, dpsi, dtheta, dphi);
  new(alobj[10]) AliAlignObjAngles("/ALIC_1/PWHE_3",
				   dvoluid, 0., 0., +displacement, dpsi, dtheta, dphi);

  // *************************    2nd step    ***************
  // Make CDB storage and put TClonesArray in
  // 
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("Alignment objects for fully misaligned geometry, i.e. applying them to TGeo has to distirbes geometry very much");
  AliCDBId id("PHOS/Align/Data",100,100);
  CDB->Put(array,id, md);
}
