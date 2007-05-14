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
//-------------------------------------------------------------------------
//   Implementation of AliGeomManager, the geometry manager class 
//   which interfaces to TGeo and the look-up table mapping unique
//   volume indices to symbolic volume names. For that it collects
//   several static methods.
//-------------------------------------------------------------------------

#include <TClass.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TObjString.h>
#include <TGeoPhysicalNode.h>
#include <TClonesArray.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>

#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliAlignObjAngles.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"

ClassImp(AliGeomManager)
  
Int_t AliGeomManager::fgLayerSize[kLastLayer - kFirstLayer] = {
  80, 160,  // ITS SPD first and second layer
  84, 176,  // ITS SDD first and second layer
  748, 950, // ITS SSD first and second layer
  36, 36,   // TPC inner and outer chambers
  90, 90, 90, 90, 90, 90,  // 6 TRD chambers' layers
  1638,     // TOF
  1, 1,     // PHOS ??
  7,        // HMPID ??
  1         // MUON ??
};

const char* AliGeomManager::fgLayerName[kLastLayer - kFirstLayer] = {
  "ITS inner pixels layer", "ITS outer pixels layer",
  "ITS inner drifts layer", "ITS outer drifts layer",
  "ITS inner strips layer", "ITS outer strips layer",
  "TPC inner chambers layer", "TPC outer chambers layer",
  "TRD chambers layer 1", "TRD chambers layer 2", "TRD chambers layer 3",
  "TRD chambers layer 4", "TRD chambers layer 5", "TRD chambers layer 6",
  "TOF layer",
  "?","?",
  "HMPID layer",
  "?"
};

TString* AliGeomManager::fgSymName[kLastLayer - kFirstLayer] = {
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,0x0,
  0x0,0x0,0x0,
  0x0,
  0x0,0x0,
  0x0,
  0x0
};

TGeoPNEntry** AliGeomManager::fgPNEntry[kLastLayer - kFirstLayer] = {
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,0x0,
  0x0,0x0,0x0,
  0x0,
  0x0,0x0,
  0x0,
  0x0
};

AliAlignObj** AliGeomManager::fgAlignObjs[kLastLayer - kFirstLayer] = {
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,0x0,
  0x0,0x0,0x0,
  0x0,
  0x0,0x0,
  0x0,
  0x0
};

AliGeomManager* AliGeomManager::fgInstance = 0x0;

TGeoManager* AliGeomManager::fgGeometry = 0x0;

//_____________________________________________________________________________
AliGeomManager* AliGeomManager::Instance()
{
// returns AliGeomManager instance (singleton)

	if (!fgInstance) {
	  fgInstance = new AliGeomManager();
	  fgInstance->Init();
	}
	return fgInstance;
}

//_____________________________________________________________________________
void AliGeomManager::Init()
{
// initialization
  if(!gGeoManager) AliFatal("Impossible to initialize AliGeomManager without an active geometry");
  fgGeometry = gGeoManager;
  InitSymNamesLUT();
  InitPNEntriesLUT();
}

//_____________________________________________________________________________
AliGeomManager::AliGeomManager():
  TObject(),

  //  fgGeometry(NULL),
  fAlignObjArray(NULL)
{
  // default constructor
}

//_____________________________________________________________________________
AliGeomManager::~AliGeomManager()
{
  // dummy destructor
  if(fAlignObjArray) fAlignObjArray->Delete();
  delete fAlignObjArray;
}

//_____________________________________________________________________________
Int_t AliGeomManager::LayerSize(Int_t layerId)
{
  // Get the layer size for layer corresponding to layerId.
  // Implemented only for ITS,TPC,TRD,TOF and HMPID
  //
  if (layerId < kFirstLayer || layerId >= kLastLayer) {
    AliErrorClass(Form("Invalid layer index %d ! Layer range is (%d -> %d) !",layerId,kFirstLayer,kLastLayer));
    return 0;
  }
  else {
    return fgLayerSize[layerId - kFirstLayer];
 }
}

//_____________________________________________________________________________
const char* AliGeomManager::LayerName(Int_t layerId)
{
  // Get the layer name corresponding to layerId.
  // Implemented only for ITS,TPC,TRD,TOF and HMPID
  //
  if (layerId < kFirstLayer || layerId >= kLastLayer) {
    AliErrorClass(Form("Invalid layer index %d ! Layer range is (%d -> %d) !",layerId,kFirstLayer,kLastLayer));
    return "Invalid Layer!";
  }
  else {
    return fgLayerName[layerId - kFirstLayer];
 }
}

//_____________________________________________________________________________
UShort_t AliGeomManager::LayerToVolUID(ELayerID layerId, Int_t modId)
{
  // From detector (layer) name and module number (according to detector
  // internal numbering) build the unique numerical identity of that volume
  // inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  // NO check for validity of given modId inside the layer for speed's sake.
  //
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
UShort_t AliGeomManager::LayerToVolUID(Int_t   layerId, Int_t modId)
{
  // From detector (layer) name and module number (according to detector
  // internal numbering) build the unique numerical identity of that volume
  // inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  // NO check for validity of given modId inside the layer for speed's sake.
  //
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
UShort_t AliGeomManager::LayerToVolUIDSafe(ELayerID layerId, Int_t modId)
{
  // From detector (layer) name and module number (according to detector
  // internal numbering) build the unique numerical identity of that volume
  // inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  // Check validity of given modId inside the layer.
  //
  if(modId < 0 || modId >= LayerSize(layerId)){
    AliErrorClass(Form("Invalid volume id %d ! Range of valid ids for layer \"%s\" is [0, %d] !",modId,LayerName(layerId),LayerSize(layerId)-1));
    return 0;
  }
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
UShort_t AliGeomManager::LayerToVolUIDSafe(Int_t layerId, Int_t modId)
{
  // From detector (layer) name and module number (according to detector
  // internal numbering) build the unique numerical identity of that volume
  // inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  // Check validity of given modId inside the layer.
  //
  if(modId < 0 || modId >= LayerSize(layerId)){
    AliErrorClass(Form("Invalid volume id %d ! Range of valid ids for layer \"%s\" is [0, %d] !",modId,LayerName(layerId),LayerSize(layerId)-1));
    return 0;
  }
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
AliGeomManager::ELayerID AliGeomManager::VolUIDToLayer(UShort_t voluid, Int_t &modId)
{
  // From voluid, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), return
  // the identity of the layer to which that volume belongs and sets the
  // argument modId to the identity of that volume internally to the layer.
  // NO check for validity of given voluid for speed's sake.
  //
  modId = voluid & 0x7ff;

  return VolUIDToLayer(voluid);
}

//_____________________________________________________________________________
AliGeomManager::ELayerID AliGeomManager::VolUIDToLayer(UShort_t voluid)
{
  // From voluid, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), return
  // the identity of the layer to which that volume belongs
  // NO check for validity of given voluid for speed's sake.
  //
  return ELayerID(voluid >> 11);
}

//_____________________________________________________________________________
AliGeomManager::ELayerID AliGeomManager::VolUIDToLayerSafe(UShort_t voluid, Int_t &modId)
{
  // From voluid, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), returns
  // the identity of the layer to which that volume belongs and sets the
  // argument modId to the identity of that volume internally to the layer.
  // Checks the validity of the given voluid
  //
  ELayerID layId = VolUIDToLayerSafe(voluid);
  if(layId){
    Int_t mId = Int_t(voluid & 0x7ff);
    if( mId>=0 && mId<LayerSize(layId)){
      modId = mId;
      return layId;
    }
  }

  AliErrorClass(Form("Invalid unique volume id: %d !",voluid));
  modId = -1;
  return kInvalidLayer;

}

//_____________________________________________________________________________
AliGeomManager::ELayerID AliGeomManager::VolUIDToLayerSafe(UShort_t voluid)
{
  // From voluid, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), returns
  // the identity of the layer to which that volume belongs
  // Checks the validity of the given voluid
  //
  if( (voluid >> 11) < kLastLayer)  return ELayerID(voluid >> 11);

  AliErrorClass(Form("Invalid layer id: %d !",(voluid >> 11)));
  return kInvalidLayer;

}

//_____________________________________________________________________________
Bool_t AliGeomManager::GetFromGeometry(const char *symname, AliAlignObj &alobj)
{
  // Get the alignment object which corresponds to the symbolic volume name
  // symname (in case equal to the TGeo volume path)
  // The method is extremely slow due to the searching by string,
  // therefore it should be used with great care!!
  // This method returns FALSE if the symname of the object was not
  // valid neither to get a TGeoPEntry nor as a volume path, or if the path
  // associated to the TGeoPNEntry was not valid.
  //

  // Reset the alignment object
  alobj.SetPars(0,0,0,0,0,0);
  alobj.SetSymName(symname);

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  if (!gGeoManager->GetListOfPhysicalNodes()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't contain any aligned nodes!");
    return kFALSE;
  }

  ReactIfChangedGeom();
  
  const char *path;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
  }else{
    AliWarningClass(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path = symname;
  }
  TObjArray* nodesArr = gGeoManager->GetListOfPhysicalNodes();
  TGeoPhysicalNode* node = NULL;
  for (Int_t iNode = 0; iNode < nodesArr->GetEntriesFast(); iNode++) {
    TGeoPhysicalNode* tempNode = (TGeoPhysicalNode*) nodesArr->UncheckedAt(iNode);
    const char *nodePath = tempNode->GetName();
    if (strcmp(path,nodePath) == 0) {
      node = tempNode;
      break;
    }
  }

  if (!node) {
    if (!gGeoManager->cd(path)) {
      AliErrorClass(Form("%s not valid neither as symbolic volume name nor as volume path!",path));
      return kFALSE;
    }
    else {
      AliWarningClass(Form("Volume (%s) has not been misaligned!",path));
      return kTRUE;
    }
  }

  TGeoHMatrix align,gprime,g,ginv,l;
  gprime = *node->GetMatrix();
  l = *node->GetOriginalMatrix();
  g = *node->GetMatrix(node->GetLevel()-1);
  g *= l;
  ginv = g.Inverse();
  align = gprime * ginv;

  return alobj.SetMatrix(align);
}


//_____________________________________________________________________________
void  AliGeomManager::InitAlignObjFromGeometry()
{
 // Loop over all alignable volumes and extract
 // the corresponding alignment objects from
 // the TGeo geometry

  if(fgAlignObjs[0]) return;
  
  ReactIfChangedGeom();
  
  for (Int_t iLayer = kFirstLayer; iLayer < AliGeomManager::kLastLayer; iLayer++) {
    fgAlignObjs[iLayer-kFirstLayer] = new AliAlignObj*[LayerSize(iLayer)];
    for (Int_t iModule = 0; iModule < LayerSize(iLayer); iModule++) {
      UShort_t volid = LayerToVolUID(iLayer,iModule);
      fgAlignObjs[iLayer-kFirstLayer][iModule] = new AliAlignObjAngles("",volid,0,0,0,0,0,0,kTRUE);
      const char *symname = SymName(volid);
      if (!GetFromGeometry(symname, *fgAlignObjs[iLayer-kFirstLayer][iModule]))
	AliErrorClass(Form("Failed to extract the alignment object for the volume (ID=%d and path=%s) !",volid,symname));
    }
  }
  
}

//_____________________________________________________________________________
AliAlignObj* AliGeomManager::GetAlignObj(UShort_t voluid) {
  // Returns the alignment object for given volume ID
  //
  Int_t modId;
  ELayerID layerId = VolUIDToLayer(voluid,modId);
  return GetAlignObj(layerId,modId);
}

//_____________________________________________________________________________
AliAlignObj* AliGeomManager::GetAlignObj(ELayerID layerId, Int_t modId)
{
  // Returns pointer to alignment object given its layer and module ID
  //
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }
  InitAlignObjFromGeometry();

  return fgAlignObjs[layerId-kFirstLayer][modId];
}

//_____________________________________________________________________________
const char* AliGeomManager::SymName(UShort_t voluid) {
  // Returns the symbolic volume name for given volume ID
  //
  Int_t modId;
  ELayerID layerId = VolUIDToLayer(voluid,modId);
  return SymName(layerId,modId);
}

//_____________________________________________________________________________
const char* AliGeomManager::SymName(ELayerID layerId, Int_t modId)
{
  // Returns the symbolic volume name given for a given layer
  // and module ID
  //
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }
  ReactIfChangedGeom();

  return fgSymName[layerId-kFirstLayer][modId].Data();
}

//_____________________________________________________________________________
void AliGeomManager::InitSymNamesLUT()
{
  // Initialize the look-up table which associates the unique
  // numerical identity of each alignable volume to the 
  // corresponding symbolic volume name
  // The LUTs are static; they are created at the creation of the
  // AliGeomManager instance and recreated if the geometry has changed
  //

  for (Int_t iLayer = 0; iLayer < (kLastLayer - kFirstLayer); iLayer++){
    if(!fgSymName[iLayer]) fgSymName[iLayer]=new TString[fgLayerSize[iLayer]];
  }

  TString symname;
  Int_t modnum; // in the following, set it to 0 at the start of each layer

  /*********************       ITS layers  ***********************/
  TString strSPD = "ITS/SPD";
  TString strSDD = "ITS/SDD";
  TString strSSD = "ITS/SSD";
  TString strStave = "/Stave";
  TString strLadder = "/Ladder";
  TString strSector = "/Sector";
  TString strSensor = "/Sensor";
  TString strEntryName1;
  TString strEntryName2;


  /*********************       SPD layer1  ***********************/
  {
    modnum = 0;

    for(Int_t c1 = 1; c1<=10; c1++){
      strEntryName1 = strSPD;
      strEntryName1 += 0;
      strEntryName1 += strSector;
      strEntryName1 += (c1-1);
      for(Int_t c2 =1; c2<=2; c2++){
	strEntryName2 = strEntryName1;
	strEntryName2 += strStave;
	strEntryName2 += (c2-1);
	for(Int_t c3 =1; c3<=4; c3++){
	  symname = strEntryName2;
	  symname += strLadder;
	  symname += (c3-1);
	  fgSymName[kSPD1-kFirstLayer][modnum] = symname.Data();
	  modnum++;
	}
      }
    }
  }
  
  /*********************       SPD layer2  ***********************/
  {
    modnum = 0;

    for(Int_t c1 = 1; c1<=10; c1++){
      strEntryName1 = strSPD;
      strEntryName1 += 1;
      strEntryName1 += strSector;
      strEntryName1 += (c1-1);
      for(Int_t c2 =1; c2<=4; c2++){
	strEntryName2 = strEntryName1;
	strEntryName2 += strStave;
	strEntryName2 += (c2-1);
	for(Int_t c3 =1; c3<=4; c3++){
	  symname = strEntryName2;
	  symname += strLadder;
	  symname += (c3-1);
	  fgSymName[kSPD2-kFirstLayer][modnum] = symname.Data();
	  modnum++;
	}
      }
    }
  }

  /*********************       SDD layer1  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=14; c1++){
      strEntryName1 = strSDD;
      strEntryName1 += 2;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 =1; c2<=6; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgSymName[kSDD1-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }

  /*********************       SDD layer2  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=22; c1++){
      strEntryName1 = strSDD;
      strEntryName1 += 3;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 = 1; c2<=8; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgSymName[kSDD2-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }

  /*********************       SSD layer1  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=34; c1++){
      strEntryName1 = strSSD;
      strEntryName1 += 4;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 = 1; c2<=22; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgSymName[kSSD1-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }

  /*********************       SSD layer2  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=38; c1++){
      strEntryName1 = strSSD;
      strEntryName1 += 5;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 = 1; c2<=25; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgSymName[kSSD2-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }


  /***************    TPC inner and outer layers    ****************/
  TString sAsector="TPC/EndcapA/Sector";
  TString sCsector="TPC/EndcapC/Sector";
  TString sInner="/InnerChamber";
  TString sOuter="/OuterChamber";
  
  /***************    TPC inner chambers' layer    ****************/
  {
    modnum = 0;
    
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sAsector;
      symname += cnt;
      symname += sInner;
      fgSymName[kTPC1-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sCsector;
      symname += cnt;
      symname += sInner;
      fgSymName[kTPC1-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
  }

  /***************    TPC outer chambers' layer    ****************/
  {
    modnum = 0;
    
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sAsector;
      symname += cnt;
      symname += sOuter;
      fgSymName[kTPC2-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sCsector;
      symname += cnt;
      symname += sOuter;
      fgSymName[kTPC2-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
  }    

  /*********************       TOF layer   ***********************/
  {
    modnum=0;
    
    Int_t nstrA=15;
    Int_t nstrB=19;
    Int_t nstrC=19;
    Int_t nSectors=18;
    Int_t nStrips=nstrA+2*nstrB+2*nstrC;
    
    TString snSM  = "TOF/sm";
    TString snSTRIP = "/strip";
    
    for (Int_t isect = 0; isect < nSectors; isect++) {
      for (Int_t istr = 1; istr <= nStrips; istr++) {	
	symname  = snSM;
	symname += Form("%02d",isect);
	symname += snSTRIP;
	symname += Form("%02d",istr);
	fgSymName[kTOF-kFirstLayer][modnum] = symname.Data();	
	modnum++;
      }
    }
  } 

  /*********************      HMPID layer   ***********************/
  {
    TString str = "/HMPID/Chamber";

    for (modnum=0; modnum < 7; modnum++) {
      symname = str;
      symname += modnum;
      fgSymName[kHMPID-kFirstLayer][modnum] = symname.Data();
    }
  }

  /*********************      TRD layers 1-6   *******************/
  //!! 6 layers with index increasing in outwards direction
  {
    Int_t arTRDlayId[6] = {kTRD1, kTRD2, kTRD3, kTRD4, kTRD5, kTRD6};

    TString snStr  = "TRD/sm";
    TString snApp1 = "/st";
    TString snApp2 = "/pl";
    
    for(Int_t layer=0; layer<6; layer++){
      modnum=0;
      for (Int_t isect = 0; isect < 18; isect++) {
	for (Int_t icham = 0; icham < 5; icham++) {
	  symname  = snStr;
	  symname += Form("%02d",isect);
	  symname += snApp1;
	  symname += icham;
	  symname += snApp2;
	  symname += layer;
	  fgSymName[arTRDlayId[layer]-kFirstLayer][modnum] = symname.Data();
	  modnum++;
	}
      }
    }
  }
}

//_____________________________________________________________________________
void AliGeomManager::InitPNEntriesLUT()
{
  // Initialize the look-up table which associates the unique
  // numerical identity of each alignable volume to the
  // corresponding TGeoPNEntry.
  // The LUTs are static; they are created at the creation of the
  // AliGeomManager instance and recreated if the geometry has changed
  //

  for (Int_t iLayer = 0; iLayer < (kLastLayer - kFirstLayer); iLayer++){
    if(!fgPNEntry[iLayer]) fgPNEntry[iLayer] = new TGeoPNEntry*[fgLayerSize[iLayer]];
  }

  for (Int_t iLayer = 0; iLayer < (kLastLayer-kFirstLayer); iLayer++){
    for(Int_t modnum=0; modnum<LayerSize(iLayer); modnum++){
      fgPNEntry[iLayer-kFirstLayer][modnum] = gGeoManager->GetAlignableEntry(fgSymName[iLayer-kFirstLayer][modnum].Data());
    }
  }
}

//______________________________________________________________________
TGeoHMatrix* AliGeomManager::GetMatrix(TGeoPNEntry* pne) 
{
  // Get the transformation matrix for a given PNEntry
  // by quering the TGeoManager

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if (pnode) return pnode->GetMatrix();

  const char* path = pne->GetTitle();
  if (!gGeoManager->cd(path)) {
    AliErrorClass(Form("Volume path %s not valid!",path));
    return NULL;
  }
  return gGeoManager->GetCurrentMatrix();
}

//______________________________________________________________________
TGeoHMatrix* AliGeomManager::GetMatrix(Int_t index) 
{
  // Get the global transformation matrix for a given alignable volume
  // identified by its unique ID 'index' by quering the TGeoManager

  ReactIfChangedGeom();

  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  return GetMatrix(pne);
}

//______________________________________________________________________
TGeoHMatrix* AliGeomManager::GetMatrix(const char* symname) 
{
  // Get the global transformation matrix for a given alignable volume
  //  identified by its symbolic name 'symname' by quering the TGeoManager

  ReactIfChangedGeom();
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if (!pne) return NULL;

  return GetMatrix(pne);
}

//______________________________________________________________________
Bool_t AliGeomManager::GetTranslation(Int_t index, Double_t t[3]) 
{
  // Get the translation vector for a given module 'index'
  // by quering the TGeoManager

  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *trans = m->GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliGeomManager::GetRotation(Int_t index, Double_t r[9]) 
{
  // Get the rotation matrix for a given module 'index'
  // by quering the TGeoManager

  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *rot = m->GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliGeomManager::GetOrigGlobalMatrix(const char *symname, TGeoHMatrix &m)
{
  // The method returns global matrix for the ideal detector geometry
  // Symname identifies either the corresponding TGeoPNEntry or directly
  // the volume path. The output global matrix is stored in 'm'.
  // Returns kFALSE in case TGeo has not been initialized or the symname
  // is invalid.
  //

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliErrorClass("Can't get the original global matrix! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }
  
  if (!gGeoManager->GetListOfPhysicalNodes()) {
    AliWarningClass("gGeoManager doesn't contain any aligned nodes!");
    if (!gGeoManager->cd(symname)) {
      AliErrorClass(Form("Volume path %s not valid!",symname));
      return kFALSE;
    }
    else {
      m = *gGeoManager->GetCurrentMatrix();
      return kTRUE;
    }
  }

  const char* path = NULL;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
  }else{
    AliWarningClass(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path=symname;
  }

  if (!gGeoManager->CheckPath(path)) {
    AliErrorClass(Form("Volume path %s not valid!",path));
    return kFALSE;
  }

  m.Clear();

  TIter next(gGeoManager->GetListOfPhysicalNodes());
  gGeoManager->cd(path);

  while(gGeoManager->GetLevel()){

    TGeoPhysicalNode *physNode = NULL;
    next.Reset();
    TGeoNode *node = gGeoManager->GetCurrentNode();
    while ((physNode=(TGeoPhysicalNode*)next())) 
      if (physNode->GetNode() == node) break;

    TGeoMatrix *lm = NULL;
    if (physNode) {
        lm = physNode->GetOriginalMatrix();
	if (!lm) lm = node->GetMatrix();
    } else
      lm = node->GetMatrix();

    m.MultiplyLeft(lm);

    gGeoManager->CdUp();
  }

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliGeomManager::GetOrigGlobalMatrix(Int_t index, TGeoHMatrix &m)
{
  // Get the original (ideal geometry) TGeo matrix for
  // a given module identified by 'index'.
  // The method is slow, so it should be used
  // with great care.

  m.Clear();

  ReactIfChangedGeom();
  const char *symname = SymName(index);
  if (!symname) return kFALSE;

  return GetOrigGlobalMatrix(symname,m);
}

//______________________________________________________________________
Bool_t AliGeomManager::GetOrigTranslation(Int_t index, Double_t t[3]) 
{
  // Get the original translation vector (ideal geometry)
  // for a given module 'index' by quering the TGeoManager

  TGeoHMatrix m;
  if (!GetOrigGlobalMatrix(index,m)) return kFALSE;

  Double_t *trans = m.GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliGeomManager::GetOrigRotation(Int_t index, Double_t r[9]) 
{
  // Get the original rotation matrix (ideal geometry)
  // for a given module 'index' by quering the TGeoManager

  TGeoHMatrix m;
  if (!GetOrigGlobalMatrix(index,m)) return kFALSE;

  Double_t *rot = m.GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
const TGeoHMatrix* AliGeomManager::GetTracking2LocalMatrix(Int_t index)
{
  // Get the matrix which transforms from the tracking to local r.s.
  // The method queries directly the TGeoPNEntry

  ReactIfChangedGeom();
  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  const TGeoHMatrix *m = pne->GetMatrix();
  if (!m)
    AliErrorClass(Form("TGeoPNEntry (%s) contains no matrix !",pne->GetName()));

  return m;
}

//______________________________________________________________________
Bool_t AliGeomManager::GetTrackingMatrix(Int_t index, TGeoHMatrix &m)
{
  // Get the matrix which transforms from the tracking r.s. to
  // the global one.
  // Returns kFALSE in case of error.

  m.Clear();

  TGeoHMatrix *m1 = GetMatrix(index);
  if (!m1) return kFALSE;

  const TGeoHMatrix *m2 = GetTracking2LocalMatrix(index);
  if (!m2) return kFALSE;

  m = *m1;
  m.Multiply(m2);

  return kTRUE;
}

//_____________________________________________________________________________
TGeoPNEntry* AliGeomManager::GetPNEntry(Int_t voluid) {
  // Returns the TGeoPNEntry for the given global volume ID "voluid"
  //
  Int_t modId;
  ELayerID layerId = VolUIDToLayer(voluid,modId);
  return GetPNEntry(layerId,modId);
}

//_____________________________________________________________________________
TGeoPNEntry* AliGeomManager::GetPNEntry(UShort_t voluid) {
  // Returns the TGeoPNEntry for the given global volume ID "voluid"
  //
  Int_t modId;
  ELayerID layerId = VolUIDToLayer(voluid,modId);
  return GetPNEntry(layerId,modId);
}

//_____________________________________________________________________________
TGeoPNEntry* AliGeomManager::GetPNEntry(ELayerID layerId, Int_t modId)
{
  // Returns the TGeoPNEntry for a given layer
  // and module ID
  //
  ReactIfChangedGeom();
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }

  return fgPNEntry[layerId-kFirstLayer][modId];
}

//_____________________________________________________________________________
Bool_t AliGeomManager::ApplyAlignObjsFromCDB(const char* AlignDetsList)
{
  // Calls AddAlignObjsFromCDBSingleDet for the detectors appearing in
  // the list passed as argument (called by AliSimulation and
  // AliReconstruction)
  // Read the alignment objects from CDB.
  // Each detector is supposed to have the
  // alignment objects in DET/Align/Data CDB path.
  // All the detector objects are then collected,
  // sorted by geometry level (starting from ALIC) and
  // then applied to the TGeo geometry.
  // Finally an overlaps check is performed.
  //
 
  if(!fAlignObjArray) fAlignObjArray = new TObjArray();
  fAlignObjArray->Clear();  	
  fAlignObjArray->SetOwner(0);

  TString alObjsNotLoaded="";
  TString alObjsLoaded="";

  TString AlignDetsString(AlignDetsList);
  TObjArray *detsarr = AlignDetsString.Tokenize(' ');
  TIter iter(detsarr);
  TObjString *str = 0;
  
  while((str = (TObjString*) iter.Next())){
    TString det(str->String());
    AliInfo(Form("Loading alignment objs for %s",det.Data()));
    if(!LoadAlignObjsFromCDBSingleDet(det.Data())){
      alObjsNotLoaded += det.Data();
      alObjsNotLoaded += " ";
    } else {
      alObjsLoaded += det.Data();
      alObjsLoaded += " ";
    }
  }

  if(!alObjsLoaded.IsNull()) AliInfo(Form("Alignment objects loaded for: %s",
					alObjsLoaded.Data()));
  if(!alObjsNotLoaded.IsNull()) AliInfo(Form("Didn't/couldn't load alignment objects for: %s",
					   alObjsNotLoaded.Data()));
 
  return(ApplyAlignObjsToGeom(fAlignObjArray));
}

//_____________________________________________________________________________
Bool_t AliGeomManager::LoadAlignObjsFromCDBSingleDet(const char* detName)
{
  // Adds the alignable objects found in the CDBEntry for the detector
  // passed as argument to the array of all alignment objects to be applyed
  // to geometry
  //
  // Fills array of single detector's alignable objects from CDB
  
  AliDebug(2, Form("Loading alignment objs for detector: %s",detName));
  
  AliCDBEntry *entry;
  	
  AliCDBPath path(detName,"Align","Data");
	
  entry=AliCDBManager::Instance()->Get(path.GetPath());
  if(!entry){ 
  	AliDebug(2,Form("Couldn't load alignment data for detector %s",detName));
	return kFALSE;
  }
  entry->SetOwner(1);
  TClonesArray *alignArray = (TClonesArray*) entry->GetObject();	
  alignArray->SetOwner(0);
  AliDebug(2,Form("Found %d alignment objects for %s",
			alignArray->GetEntries(),detName));

  AliAlignObj *alignObj=0;
  TIter iter(alignArray);
	
  // loop over align objects in detector
  while( ( alignObj=(AliAlignObj *) iter.Next() ) ){
  	fAlignObjArray->Add(alignObj);
  }
  // delete entry --- Don't delete, it is cached!
	
  AliDebug(2, Form("fAlignObjArray entries: %d",fAlignObjArray->GetEntries() ));
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliGeomManager::ApplyAlignObjsToGeom(TObjArray* alObjArray)
{
  // Read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray alObjArray and apply them to gGeoManager
  //
  ReactIfChangedGeom();

  alObjArray->Sort();
  Int_t nvols = alObjArray->GetEntriesFast();

  Bool_t flag = kTRUE;

  for(Int_t j=0; j<nvols; j++)
    {
      AliAlignObj* alobj = (AliAlignObj*) alObjArray->UncheckedAt(j);
      if (alobj->ApplyToGeometry() == kFALSE) flag = kFALSE;
    }

  if (AliDebugLevelClass() >= 1) {
    gGeoManager->GetTopNode()->CheckOverlaps(1);
    TObjArray* ovexlist = gGeoManager->GetListOfOverlaps();
    if(ovexlist->GetEntriesFast()){  
      AliError("The application of alignment objects to the geometry caused huge overlaps/extrusions!");
   }
  }

  return flag;

}

//_____________________________________________________________________________
Bool_t AliGeomManager::ApplyAlignObjsToGeom(const char* fileName, const char* clArrayName)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the file fileName and apply
  // them to the geometry
  //

  TFile* inFile = TFile::Open(fileName,"READ");
  if (!inFile || !inFile->IsOpen()) {
    AliErrorClass(Form("Could not open file %s !",fileName));
    return kFALSE;
  }

  TClonesArray* alObjArray = ((TClonesArray*) inFile->Get(clArrayName));
  inFile->Close();
  if (!alObjArray) {
    AliErrorClass(Form("Could not get array (%s) from file (%s) !",clArrayName,fileName));
    return kFALSE;
  }

  return ApplyAlignObjsToGeom(alObjArray);

}

//_____________________________________________________________________________
Bool_t AliGeomManager::ApplyAlignObjsToGeom(AliCDBParam* param, AliCDBId& Id)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the AliCDBEntry identified by
  // param (to get the AliCDBStorage) and Id; apply the alignment objects
  // to the geometry
  //

  AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage(param);
  AliCDBEntry* entry = storage->Get(Id);
  TClonesArray* AlObjArray = ((TClonesArray*) entry->GetObject());

  return ApplyAlignObjsToGeom(AlObjArray);

}

//_____________________________________________________________________________
Bool_t AliGeomManager::ApplyAlignObjsToGeom(const char* uri, const char* path, Int_t runnum, Int_t version, Int_t sversion)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the AliCDBEntry identified by
  // param (to get the AliCDBStorage) and Id; apply the alignment objects
  // to the geometry
  //

  AliCDBParam* param = AliCDBManager::Instance()->CreateParameter(uri);
  AliCDBId id(path, runnum, runnum, version, sversion);

  return ApplyAlignObjsToGeom(param, id);

}

//_____________________________________________________________________________
Bool_t AliGeomManager::ApplyAlignObjsToGeom(const char* detName, Int_t runnum, Int_t version, Int_t sversion)
{
  // read collection of alignment objects (AliAlignObj derived) saved
  // in the TClonesArray ClArrayName in the AliCDBEntry identified by
  // param (to get the AliCDBStorage) and Id; apply the alignment objects
  // to the geometry
  //

  AliCDBPath path(detName,"Align","Data");
  AliCDBEntry* entry = AliCDBManager::Instance()->Get(path.GetPath(),runnum,version,sversion);

  if(!entry) return kFALSE;
  TClonesArray* AlObjArray = ((TClonesArray*) entry->GetObject());

  return ApplyAlignObjsToGeom(AlObjArray);
}

//_____________________________________________________________________________
void AliGeomManager::ReactIfChangedGeom()
{
  // Check if the TGeo geometry has changed. In that case reinitialize the
  // look-up tables
  //
  if(HasGeomChanged())
    {
      fgGeometry = gGeoManager;
      InitSymNamesLUT();
      InitPNEntriesLUT();
    }
}
