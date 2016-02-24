// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveITSDigitsInfo.h"
#include <AliEveEventManager.h>

#include <AliITS.h>
#include <AliITSInitGeometry.h>
#include <AliITSgeomTGeo.h>
#include <AliITSsegmentationSPD.h>
#include <AliITSsegmentationSDD.h>
#include <AliITSsegmentationSSD.h>
#include <AliITSDDLModuleMapSDD.h>

#include <AliITSCalibrationSPD.h>
#include <AliITSCalibrationSDD.h>
#include <AliITSCalibrationSSD.h>
#include <AliITSdigit.h>
#include <AliITSdigitSPD.h>

#include <AliCDBEntry.h>
#include <AliCDBManager.h>

#include <AliRawReader.h>
#include <AliITSRawStreamSPD.h>
#include <AliITSRawStreamSDD.h>
#include <AliITSRawStreamSDDCompressed.h>
#include <AliITSRawStreamSSD.h>

#include <TGeoMatrix.h>
#include <TEveTrans.h>

#include <TMath.h>
#include <TVector3.h>
#include <TTree.h>

//==============================================================================
//==============================================================================
// AliEveITSModuleSelection
//==============================================================================

//______________________________________________________________________________
//
// Helper for selecting a range of ITS modules by type, layer, phi and
// theta. Taken as an argument to AliEveITSDigitsInfo::GetModuleIDs().

ClassImp(AliEveITSModuleSelection)

AliEveITSModuleSelection::AliEveITSModuleSelection():
  fType(-1),
  fLayer(-1),
  fMinPhi(-TMath::Pi()),
  fMaxPhi(TMath::Pi()),
  fMinTheta(-TMath::Pi()),
  fMaxTheta(TMath::Pi())
{
  // Constructor.
}


//==============================================================================
//==============================================================================
// AliEveITSDigitsInfo
//==============================================================================

//______________________________________________________________________________
//
// Stores ITS geometry information and event-data in format suitable
// for visualization.

ClassImp(AliEveITSDigitsInfo)

AliITSDDLModuleMapSDD *AliEveITSDigitsInfo::fgDDLMapSDD  = 0;
TObjArray             *AliEveITSDigitsInfo::fgDeadModSPD = 0;
TObjArray             *AliEveITSDigitsInfo::fgDeadModSDD = 0;
TObjArray             *AliEveITSDigitsInfo::fgDeadModSSD = 0;

/******************************************************************************/

AliEveITSDigitsInfo::AliEveITSDigitsInfo() :
  TObject(),
  TEveRefCnt(),
  fTree       (0),
  fSegSPD     (0), fSegSDD     (0), fSegSSD     (0),
  fSPDMinVal  (0), fSSDMinVal  (0), fSDDMinVal  (0),
  fSPDMaxVal  (0), fSSDMaxVal  (0), fSDDMaxVal  (0),
  fSPDHighLim (0), fSDDHighLim (0), fSSDHighLim (0),
  fSPDmap     (),  fSDDmap     (),  fSSDmap     ()
{
  // Default constructor.

  InitInternals();
}

void AliEveITSDigitsInfo::InitInternals()
{
  // Initialize internal geometry structures, in particular the
  // module-id to transformation-matrix mapping and segmentation
  // classes and data-structures.

  static const TEveException kEH("AliEveITSDigitsInfo::InitInternals ");

  AliEveEventManager::Instance()->AssertGeometry();

  SetITSSegmentation();

  // create tables for scaling
  fSPDMinVal = 0;
  fSDDMinVal = 5;
  fSSDMinVal = 2;

  fSPDMaxVal = 1;
  fSDDMaxVal = 80;
  fSSDMaxVal = 100;

  fSPDHighLim = 1;
  fSDDHighLim = 512;
  fSSDHighLim = 1024;

  // lowest scale factor refers to unscaled ITS module
  fSPDScaleX[0] = 1;
  fSPDScaleZ[0] = 1;
  fSDDScaleX[0] = 1;
  fSDDScaleZ[0] = 1;
  fSSDScale [0] = 1;

  // spd lowest resolution
  Int_t nx = 8; // fSegSPD->Npx()/8; // 32
  Int_t nz = 6; // fSegSPD->Npz()/2; // 128
  fSPDScaleX[1] = Int_t(nx);
  fSPDScaleZ[1] = Int_t(nz);
  fSPDScaleX[2] = Int_t(nx*2);
  fSPDScaleZ[2] = Int_t(nz*2);
  fSPDScaleX[3] = Int_t(nx*3);
  fSPDScaleZ[3] = Int_t(nz*3);
  fSPDScaleX[4] = Int_t(nx*4);
  fSPDScaleZ[4] = Int_t(nz*4);

  fSDDScaleX[1] = 2;
  fSDDScaleZ[1] = 2;
  fSDDScaleX[2] = 8;
  fSDDScaleZ[2] = 8;
  fSDDScaleX[3] = 16;
  fSDDScaleZ[3] = 16;
  fSDDScaleX[4] = 25;
  fSDDScaleZ[4] = 25;

  fSSDScale[1] = 3;
  fSSDScale[2] = 9;
  fSSDScale[3] = 20;
  fSSDScale[4] = 30;

  if (fgDDLMapSDD == 0)
  {
    fgDDLMapSDD = new AliITSDDLModuleMapSDD();
    AliCDBManager *man = AliCDBManager::Instance();
    Bool_t cacheStatus = man->GetCacheFlag();
    AliCDBEntry   *ddlMapSDD = man->Get("ITS/Calib/DDLMapSDD");

    if (!ddlMapSDD)
    {
      AliWarning("SDD DDL map file retrieval from OCDB failed! - Use default DDL map");
    }
    else
    {
      AliITSDDLModuleMapSDD *ddlsdd = (AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
      if (!ddlsdd)
      {
	AliWarning("SDD DDL map object not found in OCDB file! - Use default DDL map");
      }
      else
      {
	if (!cacheStatus)
	  ddlMapSDD->SetObject(0);
	ddlMapSDD->SetOwner(kTRUE);
	fgDDLMapSDD->SetDDLMap(ddlsdd);
      }
    }

    if (!cacheStatus)
      delete ddlMapSDD;
  }

  if (fgDeadModSPD == 0)
  {
    AliCDBManager *cdb = AliCDBManager::Instance();

    AliCDBEntry *deadSPD = cdb->Get("ITS/Calib/SPDDead");

    if (!deadSPD)
    {
      AliWarning("SPD Calibration object retrieval failed!");
    }
    else
    {
      fgDeadModSPD = (TObjArray*)deadSPD->GetObject();
      fgDeadModSPD->SetOwner(kTRUE);
    }
  }

  // if (fgDeadModSDD == 0)
  // {
  //   AliCDBManager *cdb = AliCDBManager::Instance();

  //   AliCDBEntry *deadSDD = cdb->Get("ITS/Calib/SDDDead");

  //   if (!deadSDD)
  //   {
  //     AliWarning("SDD Calibration object retrieval failed!");
  //   }
  //   else
  //   {
  //     fgDeadModSDD = (TObjArray*)deadSDD->GetObject();
  //     fgDeadModSDD->SetOwner(kTRUE);
  //   }
  // }

  // if (fgDeadModSSD == 0)
  // {
  //   AliCDBManager *cdb = AliCDBManager::Instance();

  //   AliCDBEntry *deadSSD = cdb->Get("ITS/Calib/SSDDead");

  //   if (!deadSSD)
  //   {
  //     AliWarning("SSD Calibration object retrieval failed!");
  //   }
  //   else
  //   {
  //     fgDeadModSSD = (TObjArray*)deadSSD->GetObject();
  //     fgDeadModSSD->SetOwner(kTRUE);
  //   }
  // }
}

/******************************************************************************/

AliEveITSDigitsInfo::~AliEveITSDigitsInfo()
{
  // Destructor.
  // Deletes the data-maps.
  // The tree, if it was passed into the object, is not touched.

  std::map<Int_t, TClonesArray*>::iterator j;
  for (j = fSPDmap.begin(); j != fSPDmap.end(); ++j)
    delete j->second;
  for (j = fSDDmap.begin(); j != fSDDmap.end(); ++j)
    delete j->second;
  for (j = fSSDmap.begin(); j != fSSDmap.end(); ++j)
    delete j->second;

  delete fSegSPD; delete fSegSDD; delete fSegSSD;
}

/******************************************************************************/

void AliEveITSDigitsInfo::SetTree(TTree* tree)
{
  // Set digit-tree to be used for digit retrieval. Data is loaded on
  // demand.

  fTree = tree;
}

void AliEveITSDigitsInfo::ReadRaw(AliRawReader* raw, Int_t mode)
{
  // Read raw-data into internal structures. AliITSdigit is used to
  // store raw-adata for all sub-detectors.

  if ((mode & 1) || (mode & 2))
  {
    AliITSRawStreamSPD inputSPD(raw);
    TClonesArray* digits = 0;
    while (inputSPD.Next())
    {
      Int_t module = inputSPD.GetModuleID();
      Int_t column = inputSPD.GetColumn();
      Int_t row    = inputSPD.GetRow();

      if (inputSPD.IsNewModule())
      {
	digits = fSPDmap[module];
	if (digits == 0)
	  fSPDmap[module] = digits = new TClonesArray("AliITSdigit", 16);
      } else if(!digits) {
		  AliFatal("No module found\n");
	  }

      AliITSdigit* d = new ((*digits)[digits->GetEntriesFast()]) AliITSdigit();
      d->SetCoord1(column);
      d->SetCoord2(row);
      d->SetSignal(1);

      // printf("SPD: %d %d %d\n",module,column,row);
    }
    raw->Reset();
  }

  if ((mode & 4) || (mode & 8))
  {
    AliITSRawStream* inputSDD = AliITSRawStreamSDD::CreateRawStreamSDD(raw);
    inputSDD->SetDDLModuleMap(fgDDLMapSDD);
    TClonesArray* digits = 0;
    while (inputSDD->Next())
    {
      Int_t module = inputSDD->GetModuleID();

      digits = fSDDmap[module];
      if (digits == 0)
	fSDDmap[module] = digits = new TClonesArray("AliITSdigit", 0);      

      if (inputSDD->IsCompletedModule()==kFALSE && inputSDD->IsCompletedDDL()==kFALSE)
      {
	Int_t anode  = inputSDD->GetCoord1()+inputSDD->GetChannel()*AliITSsegmentationSDD::GetNAnodesPerHybrid();
	Int_t time   = inputSDD->GetCoord2();
	Int_t signal = inputSDD->GetSignal();
	AliITSdigit* d = new ((*digits)[digits->GetEntriesFast()]) AliITSdigit();
	d->SetCoord1(anode);
	d->SetCoord2(time);
	d->SetSignal(signal);
      }
      // printf("SDD: %d %d %d %d\n",module,anode,time,signal);
    }
    delete inputSDD;
    raw->Reset();
  }

  if ((mode & 16) || (mode & 32))
  {
    AliITSRawStreamSSD input(raw);
    TClonesArray* digits = 0;
    while (input.Next())
    {
      Int_t module  = input.GetModuleID();
      Int_t side    = input.GetSideFlag();
      Int_t strip   = input.GetStrip();
      Int_t signal  = input.GetSignal();

      // !!!! IsNewModule() is false in the beginning of the stream, so we also
      // !!!! check for digits == 0. Should be fixed in SSD stream.
      if (input.IsNewModule() || digits == 0)
      {
	digits = fSSDmap[module];
	if (digits == 0)
	  fSSDmap[module] = digits = new TClonesArray("AliITSdigit", 0);
      }

      AliITSdigit* d = new ((*digits)[digits->GetEntriesFast()]) AliITSdigit();
      d->SetCoord1(side);
      d->SetCoord2(strip);
      d->SetSignal(signal);

      // printf("SSD: %d %d %d %d\n",module,side,strip,signal);
    }
    raw->Reset();
  }
}

/******************************************************************************/

void AliEveITSDigitsInfo::SetITSSegmentation()
{
  // Create the segmentation objects and fill internal
  // data-structures.

  // SPD
  fSegSPD = new AliITSsegmentationSPD("TGeo");

  Float_t fNzSPD = 160;
  Float_t fHlSPD = 3.48;
  Float_t fZ1pitchSPD = 0.0425, fZ2pitchSPD = 0.0625;

  fSPDZCoord[0] = fZ1pitchSPD - fHlSPD;
  for (Int_t m=1; m<fNzSPD; m++)
  {
    Float_t dz = fZ1pitchSPD;
    if (m==31 || m==32 || m==63 || m==64 || m==95 || m==96 || m==127 || m==128)
    {
      dz = fZ2pitchSPD;
    }
    fSPDZCoord[m] = fSPDZCoord[m-1] + dz;
  }

  for (Int_t m=0; m<fNzSPD; m++)
  {
    Float_t dz = fZ1pitchSPD;
    if (m==31 || m==32 || m==63 || m==64 || m==95 || m==96 || m==127 || m==128)
    {
      dz = fZ2pitchSPD;
    }
    fSPDZCoord[m] -= dz;
  }

  // SDD
  fSegSDD = new AliITSsegmentationSDD("TGeo");
  // !!!! Set default drift speed, eventually need to get it from CDB.
  fSegSDD->SetDriftSpeed(7.3);

  // SSD
  fSegSSD = new AliITSsegmentationSSD("TGeo");
}

/******************************************************************************/

TClonesArray* AliEveITSDigitsInfo::GetDigits(Int_t mod, Int_t subdet) const
{
  // Return TClonesArray of digits for specified module and sub-detector-id.

  switch(subdet)
  {
    case 0:
    {
      std::map<Int_t, TClonesArray*>::iterator i = fSPDmap.find(mod);
      if (i != fSPDmap.end())
      {
	return i->second;
      }
      else if (fTree)
      {
	TClonesArray *digitsSPD = 0;
	TBranch *br = fTree->GetBranch("ITSDigitsSPD");
	br->SetAddress(&digitsSPD);
	br->GetEntry(mod);
	fSPDmap[mod] = digitsSPD;
	return digitsSPD;
      }
      break;
    }

    case 1:
    {
      std::map<Int_t, TClonesArray*>::iterator i = fSDDmap.find(mod);
      if (i != fSDDmap.end())
      {
	return i->second;
      }
      else if (fTree)
      {
	TClonesArray *digitsSDD = 0;
	TBranch *br = fTree->GetBranch("ITSDigitsSDD");
	br->SetAddress(&digitsSDD);
	br->GetEntry(mod);
	fSDDmap[mod] = digitsSDD;
	return digitsSDD;
      }
      break;
    }

    case 2:
    {
      std::map<Int_t, TClonesArray*>::iterator i = fSSDmap.find(mod);
      if (i != fSSDmap.end())
      {
	return i->second;
      }
      else if (fTree)
      {
	TClonesArray *digitsSSD = 0;
	TBranch *br = fTree->GetBranch("ITSDigitsSSD");
	br->SetAddress(&digitsSSD);
	br->GetEntry(mod);
	fSSDmap[mod] = digitsSSD;
	return digitsSSD;
      }
      break;
    }
  }

  return 0;
}

/******************************************************************************/

void AliEveITSDigitsInfo::GetModuleIDs(AliEveITSModuleSelection* sel,
				       std::vector<UInt_t>& ids) const
{
  // Fill the id-vector with ids of modules that satisfy conditions
  // given by the AliEveITSModuleSelection object.

  Int_t idx0 = 0, idx1 = 0;
  switch(sel->GetType())
  {
    case 0:
      idx0 = 0;
      idx1 = AliITSgeomTGeo::GetModuleIndex(3, 1, 1) - 1;
      break;
    case 1:
      idx0 = AliITSgeomTGeo::GetModuleIndex(3, 1, 1);
      idx1 = AliITSgeomTGeo::GetModuleIndex(5, 1, 1) - 1;
      break;
    case 2:
      idx0 = AliITSgeomTGeo::GetModuleIndex(5, 1, 1);
      idx1 = AliITSgeomTGeo::GetNModules() - 1;
      break;
    default:
      idx1 = 0;
      idx1 = AliITSgeomTGeo::GetNModules() - 1;
      break;
  }

  TVector3  v;
  TEveTrans mx;
  Int_t     lay, lad, det;
  for (Int_t id = idx0; id < idx1; ++id)
  {
    AliITSgeomTGeo::GetModuleId(id, lay, lad, det);
    if (sel->GetLayer() == lay || sel->GetLayer() == -1)
    {
      // check data from matrix
      mx.SetFrom(*AliITSgeomTGeo::GetMatrix(id));
      mx.GetPos(v);
      if (v.Phi()   <= sel->GetMaxPhi()   && v.Phi()   >= sel->GetMinPhi()   &&
	  v.Theta() <= sel->GetMaxTheta() && v.Theta() >= sel->GetMinTheta())
      {
	ids.push_back(id);
      }
    }
  }
}

/******************************************************************************/

Bool_t AliEveITSDigitsInfo::HasData(Int_t module, Int_t det_id) const
{
  // Return true if given module has data.

  TClonesArray *digits = GetDigits(module, det_id);
  return digits && digits->GetEntriesFast() > 0;
}

Bool_t AliEveITSDigitsInfo::IsDead (Int_t module, Int_t det_id) const
{
  // Return true if given module is dead.
  // Only implemented for SPD.

  if (det_id == 0 && fgDeadModSPD)
    return ((AliITSCalibrationSPD*) fgDeadModSPD->At(module))->IsBad();
  if (det_id == 1 && fgDeadModSDD)
    return ((AliITSCalibrationSDD*) fgDeadModSDD->At(module))->IsBad();
  if (det_id == 2 && fgDeadModSSD)
    return ((AliITSCalibrationSSD*) fgDeadModSSD->At(module))->IsBad();
  return kFALSE;
}

/******************************************************************************/

void AliEveITSDigitsInfo::Print(Option_t* ) const
{
  // Print information about stored geometry and segmentation.

  printf("*********************************************************\n");
  printf("SPD module dimension (%f,%f)\n",           fSegSPD->Dx()*0.0001, fSegSPD->Dz()*0.0001);
  printf("SPD first,last module: %d, %d\n",
	 AliITSgeomTGeo::GetModuleIndex(1,1,1), 
	 AliITSgeomTGeo::GetModuleIndex(3,1,1) - 1);
  printf("SPD num cells per module (x::%d,z::%d)\n", fSegSPD->Npx(), fSegSPD->Npz());
  Int_t iz = 0, ix = 0;
  printf("SPD dimesion of (%d,%d) in pixel(%f,%f)\n",   ix, iz, fSegSPD->Dpx(ix), fSegSPD->Dpz(iz));
  iz = 32;
  printf("SPD dimesion of pixel (%d,%d) are (%f,%f)\n", ix, iz, fSegSPD->Dpx(ix)*0.001, fSegSPD->Dpz(iz)*0.001);

  printf("*********************************************************\n");
  printf("SDD module dimension (%f,%f)\n",           fSegSDD->Dx()*0.0001, fSegSDD->Dz()*0.0001);
  printf("SDD first,last module: %d, %d\n",
	 AliITSgeomTGeo::GetModuleIndex(3,1,1),
	 AliITSgeomTGeo::GetModuleIndex(5,1,1) - 1);
  printf("SDD num cells per module (x::%d,z::%d)\n", fSegSDD->Npx(), fSegSDD->Npz());
  printf("SDD dimesion of pixel are (%f,%f)\n",      fSegSDD->Dpx(1)*0.001, fSegSDD->Dpz(1)*0.001);

  Float_t ap, an;
  printf("*********************************************************\n");
  printf("SSD module dimension (%f,%f)\n",  fSegSSD->Dx()*0.0001, fSegSSD->Dz()*0.0001);
  printf("SSD first, last module: %d, %d\n",
	 AliITSgeomTGeo::GetModuleIndex(5,1,1),
	 AliITSgeomTGeo::GetNModules() - 1);
  printf("SSD strips in module %d\n",       fSegSSD->Npx());
  printf("SSD strip sizes are (%f,%f)\n",   fSegSSD->Dpx(1), fSegSSD->Dpz(1));
  fSegSSD->SetLayer(5);  fSegSSD->Angles(ap,an);
  printf("SSD layer 5 stereoP %f stereoN %f angle\n", ap, an);
  fSegSSD->SetLayer(6);  fSegSSD->Angles(ap,an);
  printf("SSD layer 6 stereoP %f stereoN %f angle\n", ap, an);
}


/*
  printf("num cells %d,%d scaled %d,%d \n",fSegSPD->Npz(),fSegSPD->Npx(),Nz,Nx);
  printf("%d digits in AliEveITSModule %d\n",ne, module);
  Float_t zn = i*(3.48*2)/Nz - 3.48 ;
  Float_t xo =  -fSegSPD->Dx()*0.00005 + fSegSPD->Dpx(0)*od->GetCoord2()*0.0001;
  Float_t xn =  -fSegSPD->Dx()*0.00005 + j*0.0001*fSegSPD->Dx()/Nx;
  Float_t dpx = 0.0001*fSegSPD->Dx()/Nx;
  Float_t dpz = 3.48*2/Nz;
  printf("Z::original (%3f) scaled (%3f, %3f) \n", zo, zn-dpz/2, zn+dpz/2);
  printf("X::original (%3f) scaled (%3f, %3f) \n", xo, xn-dpx/2, xn+dpx/2);
  printf("%d,%d maped to %d,%d \n", od->GetCoord1(), od->GetCoord2(), i,j );
*/
