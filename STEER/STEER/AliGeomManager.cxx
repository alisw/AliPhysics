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
#include <TSystem.h>
#include <TStopwatch.h>
#include <TGeoOverlap.h>
#include <TPluginManager.h>
#include <TROOT.h>

#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliAlignObjParams.h"
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
  5, 5,     // PHOS,CPV
  7,        // HMPID ??
  1,         // MUON ??
  12        // EMCAL
};

const char* AliGeomManager::fgLayerName[kLastLayer - kFirstLayer] = {
  "ITS inner pixels layer", "ITS outer pixels layer",
  "ITS inner drifts layer", "ITS outer drifts layer",
  "ITS inner strips layer", "ITS outer strips layer",
  "TPC inner chambers layer", "TPC outer chambers layer",
  "TRD chambers layer 1", "TRD chambers layer 2", "TRD chambers layer 3",
  "TRD chambers layer 4", "TRD chambers layer 5", "TRD chambers layer 6",
  "TOF layer",
  "PHOS EMC layer","PHOS CPV layer",
  "HMPID layer", 
  "MUON ?",
  "EMCAL layer"
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
  0x0,
  0x0
};

const char* AliGeomManager::fgkDetectorName[AliGeomManager::fgkNDetectors] = {"GRP","ITS","TPC","TRD","TOF","PHOS","HMPID","EMCAL","MUON","FMD","ZDC","PMD","T0","VZERO","ACORDE"
									      // #ifdef MFT_UPGRADE	
									      //  ,"MFT"
									      // #endif 
									      ,"MFT"   // AU
};
Int_t AliGeomManager::fgNalignable[fgkNDetectors] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
// #ifdef MFT_UPGRADE
// 						     , 0
// #endif 
						     , 0    // AU
};

TGeoManager* AliGeomManager::fgGeometry = 0x0;

//_____________________________________________________________________________
void AliGeomManager::LoadGeometry(const char *geomFileName)
{
  // initialization
  // Load geometry either from a file
  // or from the corresponding CDB entry

  if(fgGeometry->IsLocked()){
      AliErrorClass("Cannot load a new geometry, the current one being locked. Setting internal geometry to null!!");
      fgGeometry = NULL;
      return;
  }

  fgGeometry = NULL;
  if (geomFileName && (!gSystem->AccessPathName(geomFileName))) {
    fgGeometry = TGeoManager::Import(geomFileName);
    AliInfoClass(Form("From now on using geometry from custom geometry file \"%s\"",geomFileName));
  }

  if (!fgGeometry) {
    AliCDBPath path("GRP","Geometry","Data");
	
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(!entry) AliFatalClass("Couldn't load geometry data from CDB!");

    entry->SetOwner(0);
    fgGeometry = (TGeoManager*) entry->GetObject();
    if (!fgGeometry) AliFatalClass("Couldn't find TGeoManager in the specified CDB entry!");
    
    AliInfoClass(Form("From now on using geometry from CDB base folder \"%s\"",
		      AliCDBManager::Instance()->GetURI("GRP/Geometry/Data")));
  }
  ResetPNEntriesLUT();
  InitPNEntriesLUT();
  InitNalignable();
}

//_____________________________________________________________________________
void AliGeomManager::SetGeometry(TGeoManager * const geom)
{
  // Load already active geometry
  if (!geom) AliFatalClass("Pointer to the active geometry is 0x0!");
  ResetPNEntriesLUT();
  fgGeometry = geom;
  InitPNEntriesLUT();
  InitNalignable();
}

//_____________________________________________________________________________
AliGeomManager::AliGeomManager():
  TObject()
{
  // default constructor
}

//_____________________________________________________________________________
AliGeomManager::~AliGeomManager()
{
  // dummy destructor
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
  if(layId != AliGeomManager::kInvalidLayer){
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

  if (!fgGeometry || !fgGeometry->IsClosed()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  if (!fgGeometry->GetListOfPhysicalNodes()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't contain any aligned nodes!");
    return kFALSE;
  }

  const char *path;
  TGeoPNEntry* pne = fgGeometry->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
  }else{
    AliWarningClass(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path = symname;
  }
  TObjArray* nodesArr = fgGeometry->GetListOfPhysicalNodes();
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
    if (!fgGeometry->cd(path)) {
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
  //
  for (Int_t iLayer = kFirstLayer; iLayer < AliGeomManager::kLastLayer; iLayer++) {
    if (!fgAlignObjs[iLayer-kFirstLayer]) {
      fgAlignObjs[iLayer-kFirstLayer] = new AliAlignObj*[LayerSize(iLayer)];
    }
    for (Int_t iModule = 0; iModule < LayerSize(iLayer); iModule++) {
      UShort_t volid = LayerToVolUID(iLayer,iModule);
      fgAlignObjs[iLayer-kFirstLayer][iModule] = new AliAlignObjParams("",volid,0,0,0,0,0,0,kTRUE);
      const char *symname = SymName(volid);
      if (!GetFromGeometry(symname, *fgAlignObjs[iLayer-kFirstLayer][iModule]))
	AliErrorClass(Form("Failed to extract the alignment object for the volume (ID=%d and path=%s) !",volid,symname));
    }
  }

}

//_____________________________________________________________________________
AliAlignObj* AliGeomManager::GetAlignObj(UShort_t voluid)
{
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
const char* AliGeomManager::SymName(UShort_t voluid)
{
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
  if(!fgGeometry){
    AliErrorClass("No geometry instance loaded yet!");
    return NULL;
  }
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }

  TGeoPNEntry* pne = fgPNEntry[layerId-kFirstLayer][modId];
  if(!pne)
  {
    AliWarningClass(Form("Module %d of layer %s is not activated!",modId,LayerName(layerId)));
    return NULL;
  }
  return pne->GetName();

}

//_____________________________________________________________________________
Bool_t AliGeomManager::CheckSymNamesLUT(const char* /*detsToBeChecked*/)
{
  // Check the look-up table which associates the unique numerical identity of
  // each alignable volume to the corresponding symbolic volume name.
  // The LUT is now held inside the geometry and handled by TGeo.
  // The method is meant to be launched when loading a geometry to verify that
  // no changes in the symbolic names have been introduced, which would prevent
  // backward compatibility with alignment objects.
  // To accept both complete and partial geometry, this method skips the check
  // for TRD and TOF volumes which are missing in the partial geometry.
  // 

//  TString detsString(detsToBeChecked);
//  if(detsString.Contains("ALL")) detsString="ITS TPC TOF TRD HMPID PHOS EMCAL";

  // Temporary measure to face the case of reconstruction over detectors not present in the geometry
  TString detsString = "";
  if(fgGeometry->CheckPath("ALIC_1/ITSV_1")) detsString+="ITS ";
  if(fgGeometry->CheckPath("ALIC_1/TPC_M_1")) detsString+="TPC ";

  TString tofsm;
  TString baseTof("ALIC_1/B077_1/BSEGMO");
  TString middleTof("_1/BTOF");
  TString trailTof("_1/FTOA_0");
  Bool_t tofActive=kFALSE;
  Bool_t tofSMs[18];
  for(Int_t sm=0; sm<18; sm++)
  {
    tofSMs[sm]=kFALSE;
    tofsm=baseTof;
    tofsm += sm;
    tofsm += middleTof;
    tofsm += sm;
    tofsm += trailTof;
    if(fgGeometry->CheckPath(tofsm.Data()))
    {
      tofActive=kTRUE;
      tofSMs[sm]=kTRUE;
    }
  }
  if(tofActive) detsString+="TOF ";
  
  TString trdsm;
  TString baseTrd("ALIC_1/B077_1/BSEGMO");
  TString middleTrd("_1/BTRD");
  TString trailTrd("_1/UTR1_1");
  Bool_t trdActive=kFALSE;
  Bool_t trdSMs[18];
  for(Int_t sm=0; sm<18; sm++)
  {
    trdSMs[sm]=kFALSE;
    trdsm=baseTrd;
    trdsm += sm;
    trdsm += middleTrd;
    trdsm += sm;
    trdsm += trailTrd;
    if(fgGeometry->CheckPath(trdsm.Data()))
    {
      trdActive=kTRUE;
      trdSMs[sm]=kTRUE;
    }
  }
  if(trdActive) detsString+="TRD ";

  if(fgGeometry->CheckPath("ALIC_1/Hmp0_0")) detsString+="HMPID ";
  
  TString phosMod, cpvMod;
  TString basePhos("ALIC_1/PHOS_");
  Bool_t phosActive=kFALSE;
  Bool_t cpvActive=kFALSE;
  Bool_t phosMods[5];
  for(Int_t pmod=0; pmod<5; pmod++)
  {
    phosMods[pmod]=kFALSE;
    phosMod = basePhos;
    phosMod += (pmod+1);
    cpvMod = phosMod;
    cpvMod += "/PCPV_1";
    if(fgGeometry->CheckPath(phosMod.Data()))
    {
      phosActive=kTRUE;
      phosMods[pmod]=kTRUE;
      if(fgGeometry->CheckPath(cpvMod.Data())) cpvActive=kTRUE;
    }
  }
  if(phosActive) detsString+="PHOS ";

  // Check over the ten EMCAL full supermodules and the two EMCAL half supermodules
  TString emcalSM;
  TString baseEmcalSM("ALIC_1/XEN1_1/SM");
  Bool_t emcalActive=kFALSE;
  Bool_t emcalSMs[12] = {kFALSE};
  for(Int_t sm=0; sm<12; sm++)
  {
    emcalSM=baseEmcalSM;
    if(sm<10){
	emcalSM += "OD_";
	emcalSM += (sm+1);
    }else{
	emcalSM += "10_";
	emcalSM += (sm-9);
    }
    if(fgGeometry->CheckPath(emcalSM.Data()))
    {
      emcalActive=kTRUE;
      emcalSMs[sm]=kTRUE;
    }
  }
  if(emcalActive) detsString+="EMCAL ";
  

  TString symname;
  const char* sname;
  TGeoPNEntry* pne = 0x0;
  Int_t uid; // global unique identity
  Int_t modnum; // unique id inside layer; in the following, set it to 0 at the start of each layer

  if(detsString.Contains("ITS")){
  /*********************       ITS layers  ***********************/
    AliDebugClass(2,"Checking consistency of symbolic names for ITS layers");
    TString strSPD = "ITS/SPD";
    TString strSDD = "ITS/SDD";
    TString strSSD = "ITS/SSD";
    TString strStave = "/Stave";
    TString strHalfStave = "/HalfStave";
    TString strLadder = "/Ladder";
    TString strSector = "/Sector";
    TString strSensor = "/Sensor";
    TString strEntryName1;
    TString strEntryName2;
    TString strEntryName3;

    /*********************       SPD layer1  ***********************/
    {
      modnum = 0;

      for(Int_t cSect = 0; cSect<10; cSect++){
	strEntryName1 = strSPD;
	strEntryName1 += 0;
	strEntryName1 += strSector;
	strEntryName1 += cSect;

	for(Int_t cStave =0; cStave<2; cStave++){
	  strEntryName2 = strEntryName1;
	  strEntryName2 += strStave;
	  strEntryName2 += cStave;

	  for (Int_t cHS=0; cHS<2; cHS++) {
	    strEntryName3 = strEntryName2;
	    strEntryName3 += strHalfStave;
	    strEntryName3 += cHS;

	    for(Int_t cLad =0; cLad<2; cLad++){
	      symname = strEntryName3;
	      symname += strLadder;
	      symname += cLad+cHS*2;
	      uid = LayerToVolUID(kSPD1,modnum++);
	      pne = fgGeometry->GetAlignableEntryByUID(uid);
	      if(!pne)
	      {
		AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
		return kFALSE;
	      }
	      sname = pne->GetName();
	      if(symname.CompareTo(sname)) 
	      {
		AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d."
		      "Expected was %s, found was %s!", uid, symname.Data(), sname));
		return kFALSE;
	      }
	    }
	  }
	}
      }
    }

    /*********************       SPD layer2  ***********************/
    {
      modnum = 0;

      for(Int_t cSect = 0; cSect<10; cSect++){
	strEntryName1 = strSPD;
	strEntryName1 += 1;
	strEntryName1 += strSector;
	strEntryName1 += cSect;

	for(Int_t cStave =0; cStave<4; cStave++){
	  strEntryName2 = strEntryName1;
	  strEntryName2 += strStave;
	  strEntryName2 += cStave;

	  for (Int_t cHS=0; cHS<2; cHS++) {
	    strEntryName3 = strEntryName2;
	    strEntryName3 += strHalfStave;
	    strEntryName3 += cHS;

	    for(Int_t cLad =0; cLad<2; cLad++){
	      symname = strEntryName3;
	      symname += strLadder;
	      symname += cLad+cHS*2;
	      uid = LayerToVolUID(kSPD2,modnum++);
	      pne = fgGeometry->GetAlignableEntryByUID(uid);
	      if(!pne)
	      {
		AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
		return kFALSE;
	      }
	      sname = pne->GetName();
	      if(symname.CompareTo(sname)) 
	      {
		AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d."
		      "Expected was %s, found was %s!", uid, symname.Data(), sname));
		return kFALSE;
	      }
	    }
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
	  uid = LayerToVolUID(kSDD1,modnum++);
	  pne = fgGeometry->GetAlignableEntryByUID(uid);
	  if(!pne)
	  {
	    AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	    return kFALSE;
	  }
	  sname = pne->GetName();
	  if(symname.CompareTo(sname)) 
	  {
	    AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		  "Expected was %s, found was %s!", uid, symname.Data(), sname));
	    return kFALSE;
	  }
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
	  uid = LayerToVolUID(kSDD2,modnum++);
	  pne = fgGeometry->GetAlignableEntryByUID(uid);
	  if(!pne)
	  {
	    AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	    return kFALSE;
	  }
	  sname = pne->GetName();
	  if(symname.CompareTo(sname)) 
	  {
	    AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		  "Expected was %s, found was %s!", uid, symname.Data(), sname));
	    return kFALSE;
	  }
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
	  uid = LayerToVolUID(kSSD1,modnum++);
	  pne = fgGeometry->GetAlignableEntryByUID(uid);
	  if(!pne)
	  {
	    AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	    return kFALSE;
	  }
	  sname = pne->GetName();
	  if(symname.CompareTo(sname)) 
	  {
	    AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		  "Expected was %s, found was %s!", uid, symname.Data(), sname));
	    return kFALSE;
	  }
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
	  uid = LayerToVolUID(kSSD2,modnum++);
	  pne = fgGeometry->GetAlignableEntryByUID(uid);
	  if(!pne)
	  {
	    AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	    return kFALSE;
	  }
	  sname = pne->GetName();
	  if(symname.CompareTo(sname)) 
	  {
	    AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		  "Expected was %s, found was %s!", uid, symname.Data(), sname));
	    return kFALSE;
	  }
	}
      }
    }

    AliDebugClass(2,"Consistency check for ITS symbolic names finished successfully.");
  }

  if(detsString.Contains("TPC"))
  {
    /***************    TPC inner and outer layers    ****************/
    
    AliDebugClass(2,"Checking consistency of symbolic names for TPC layers");
    TString sAsector="TPC/EndcapA/Sector";
    TString sCsector="TPC/EndcapC/Sector";
    TString sInner="/InnerChamber";
    TString sOuter="/OuterChamber";

    /***************    TPC inner chambers' layer    ****************/
    {
      modnum = 0;

      for(Int_t cnt=1; cnt<=18; cnt++)
      {
	symname = sAsector;
	symname += cnt;
	symname += sInner;
	uid = LayerToVolUID(kTPC1,modnum++);
	pne = fgGeometry->GetAlignableEntryByUID(uid);
	if(!pne)
	{
	  AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	  return kFALSE;
	}
	sname = pne->GetName();
	if(symname.CompareTo(sname)) 
	{
	  AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		"Expected was %s, found was %s!", uid, symname.Data(), sname));
	  return kFALSE;
	}
      }

      for(Int_t cnt=1; cnt<=18; cnt++)
      {
	symname = sCsector;
	symname += cnt;
	symname += sInner;
	uid = LayerToVolUID(kTPC1,modnum++);
	pne = fgGeometry->GetAlignableEntryByUID(uid);
	if(!pne)
	{
	  AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	  return kFALSE;
	}
	sname = pne->GetName();
	if(symname.CompareTo(sname)) 
	{
	  AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		"Expected was %s, found was %s!", uid, symname.Data(), sname));
	  return kFALSE;
	}
      }
    }

    /***************    TPC outer chambers' layer    ****************/
    {
      modnum = 0;

      for(Int_t cnt=1; cnt<=18; cnt++)
      {
	symname = sAsector;
	symname += cnt;
	symname += sOuter;
	uid = LayerToVolUID(kTPC2,modnum++);
	pne = fgGeometry->GetAlignableEntryByUID(uid);
	if(!pne)
	{
	  AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	  return kFALSE;
	}
	sname = pne->GetName();
	if(symname.CompareTo(sname)) 
	{
	  AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		"Expected was %s, found was %s!", uid, symname.Data(), sname));
	  return kFALSE;
	}
      }

      for(Int_t cnt=1; cnt<=18; cnt++)
      {
	symname = sCsector;
	symname += cnt;
	symname += sOuter;
	uid = LayerToVolUID(kTPC2,modnum++);
	pne = fgGeometry->GetAlignableEntryByUID(uid);
	if(!pne)
	{
	  AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	  return kFALSE;
	}
	sname = pne->GetName();
	if(symname.CompareTo(sname)) 
	{
	  AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		"Expected was %s, found was %s!", uid, symname.Data(), sname));
	  return kFALSE;
	}
      }
    }

    AliDebugClass(2,"Consistency check for TPC symbolic names finished successfully.");
  }

  if(detsString.Contains("TOF"))
  {
    /*********************       TOF layer   ***********************/

    AliDebugClass(2,"Checking consistency of symbolic names for TOF layers");
    modnum=0;

    Int_t nstrA=15;
    Int_t nstrB=19;
    Int_t nstrC=19;
    Int_t nSectors=18;
    Int_t nStrips=nstrA+2*nstrB+2*nstrC;

    TString snSM  = "TOF/sm";
    TString snSTRIP = "/strip";
    
    for (Int_t isect = 0; isect < nSectors; isect++) {
	if(tofSMs[isect]) AliDebugClass(3,Form("Consistency check for symnames of TOF supermodule %d.",isect));
      for (Int_t istr = 1; istr <= nStrips; istr++) {	
	symname  = snSM;
	symname += Form("%02d",isect);
	symname += snSTRIP;
	symname += Form("%02d",istr);
	uid = LayerToVolUID(kTOF,modnum++);
	if(!tofSMs[isect]) continue; // taking possible missing TOF sectors (partial geometry) into account
	if ((isect==13 || isect==14 || isect==15) && (istr >= 39 && istr <= 53)) continue; //taking holes into account
	pne = fgGeometry->GetAlignableEntryByUID(uid);
	if(!pne)
	{
	  AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	  return kFALSE;
	}
	sname = pne->GetName();
	if(symname.CompareTo(sname)) 
	{
	  AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		"Expected was %s, found was %s!", uid, symname.Data(), sname));
	  return kFALSE;
	}
      }
    }
    
    AliDebugClass(2,"Consistency check for TOF symbolic names finished successfully.");
  } 

  if(detsString.Contains("HMPID"))
  {
    /*********************      HMPID layer   ***********************/

    AliDebugClass(2,"Checking consistency of symbolic names for HMPID layers");
    TString str = "/HMPID/Chamber";

    for (modnum=0; modnum < 7; modnum++) {
      symname = str;
      symname += modnum;
      uid = LayerToVolUID(kHMPID,modnum);
      pne = fgGeometry->GetAlignableEntryByUID(uid);
      if(!pne)
      {
	AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	return kFALSE;
      }
      sname = pne->GetName();
      if(symname.CompareTo(sname)) 
      {
	AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
	      "Expected was %s, found was %s!", uid, symname.Data(), sname));
	return kFALSE;
      }
    }

    AliDebugClass(2,"Consistency check for HMPID symbolic names finished successfully.");
  }

  if(detsString.Contains("TRD"))
  {
    /*********************      TRD layers 1-6   *******************/
    //!! 6 layers with index increasing in outwards direction
    
    AliDebugClass(2,"Checking consistency of symbolic names for TRD layers");
    Int_t arTRDlayId[6] = {kTRD1, kTRD2, kTRD3, kTRD4, kTRD5, kTRD6};

    TString snStr  = "TRD/sm";
    TString snApp1 = "/st";
    TString snApp2 = "/pl";
    
    for(Int_t layer=0; layer<6; layer++){
      modnum=0;
      AliDebugClass(3,Form("Consistency check for symnames of TRD layer %d.",layer));
      for (Int_t isect = 0; isect < 18; isect++) {
	for (Int_t icham = 0; icham < 5; icham++) {
	  symname  = snStr;
	  symname += Form("%02d",isect);
	  symname += snApp1;
	  symname += icham;
	  symname += snApp2;
	  symname += layer;
	  uid = LayerToVolUID(arTRDlayId[layer],modnum++);
	  if(!trdSMs[isect]) continue;
	  if ((isect==13 || isect==14 || isect==15) && icham==2) continue; //keeping holes into account
	  pne = fgGeometry->GetAlignableEntryByUID(uid);
	  if(!pne)
	  {
	    AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	    return kFALSE;
	  }
	  sname = pne->GetName();
	  if(symname.CompareTo(sname)) 
	  {
	    AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		  "Expected was %s, found was %s!", uid, symname.Data(), sname));
	    return kFALSE;
	  }
	}
      }
    }

    AliDebugClass(2,"Consistency check for TRD symbolic names finished successfully.");
  }

  if(detsString.Contains("PHOS"))
  {
    /*********************      PHOS EMC layer   ***********************/

    AliDebugClass(2,"Checking consistency of symbolic names for PHOS layers");
    
      TString str = "PHOS/Module";
      modnum=0;

      for (Int_t iModule=0; iModule < 5; iModule++) {
	if(!phosMods[iModule]) continue;
	symname = str;
	symname += (iModule+1);
	uid = LayerToVolUID(kPHOS1,iModule);
	pne = fgGeometry->GetAlignableEntryByUID(uid);
	if(!pne)
	{
	  AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	  return kFALSE;
	}
	sname = pne->GetName();
	if(symname.CompareTo(sname)) 
	{
	  AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		"Expected was %s, found was %s!", uid, symname.Data(), sname));
	  return kFALSE;
	}
	/*********************      PHOS CPV layer   ***********************/
	if(!cpvActive) continue;
	symname += "/CPV";
	uid = LayerToVolUID(kPHOS2,iModule);
	pne = fgGeometry->GetAlignableEntryByUID(uid);
	if(!pne)
	{
	  AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	  return kFALSE;
	}
	sname = pne->GetName();
	if(symname.CompareTo(sname)) 
	{
	  AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
		"Expected was %s, found was %s!", uid, symname.Data(), sname));
	  return kFALSE;
	}
      }
    AliDebugClass(2,"Consistency check for PHOS symbolic names finished successfully.");
  }

  if(detsString.Contains("EMCAL"))
  {
    /*********************      EMCAL layer   ***********************/

    AliDebugClass(2,"Checking consistency of symbolic names for EMCAL layers");
    TString str = "EMCAL/FullSupermodule";
    modnum=0;

    for (Int_t iModule=1; iModule <= 12; iModule++) {
      if(!emcalSMs[iModule-1]) continue;
      symname = str;
      symname += iModule;
      if(iModule >10) {
	symname = "EMCAL/HalfSupermodule";
	symname += iModule-10;
      }
      modnum = iModule-1;
      uid = LayerToVolUID(kEMCAL,modnum);
      pne = fgGeometry->GetAlignableEntryByUID(uid);
      if(!pne)
      {
	AliErrorClass(Form("In the currently loaded geometry there is no TGeoPNEntry with unique id %d",uid));
	return kFALSE;
      }
      sname = pne->GetName();
      if(symname.CompareTo(sname)) 
      {
	AliErrorClass(Form("Current loaded geometry differs in the definition of symbolic name for uid %d"
	      "Expected was %s, found was %s!", uid, symname.Data(), sname));
	return kFALSE;
      }
    }
  
    AliDebugClass(2,"Consistency check for EMCAL symbolic names finished successfully.");
  }
  
  return kTRUE;

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
  if(!fgGeometry) {
    AliErrorClass("Impossible to initialize PNEntries LUT without an active geometry");
    return;
  }

  for (Int_t iLayer = 0; iLayer < (kLastLayer - kFirstLayer); iLayer++){
    if (!fgPNEntry[iLayer]) fgPNEntry[iLayer] = new TGeoPNEntry*[fgLayerSize[iLayer]];
    for(Int_t modnum=0; modnum<fgLayerSize[iLayer]; modnum++){
      fgPNEntry[iLayer][modnum] = fgGeometry->GetAlignableEntryByUID(LayerToVolUID(iLayer+1,modnum));
    }
  }
}

//______________________________________________________________________
TGeoHMatrix* AliGeomManager::GetMatrix(TGeoPNEntry * const pne) 
{
  // Get the global transformation matrix for a given PNEntry
  // by quering the TGeoManager

  if (!fgGeometry || !fgGeometry->IsClosed()) {
    AliErrorClass("Can't get the global matrix! gGeoManager doesn't exist or it is still opened!");
    return NULL;
  }

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if (pnode) return pnode->GetMatrix();

  const char* path = pne->GetTitle();
  if (!fgGeometry->cd(path)) {
    AliErrorClass(Form("Volume path %s not valid!",path));
    return NULL;
  }
  return fgGeometry->GetCurrentMatrix();
}

//______________________________________________________________________
TGeoHMatrix* AliGeomManager::GetMatrix(Int_t index) 
{
  // Get the global transformation matrix for a given alignable volume
  // identified by its unique ID 'index' by quering the TGeoManager

  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  return GetMatrix(pne);
}

//______________________________________________________________________
TGeoHMatrix* AliGeomManager::GetMatrix(const char* symname) 
{
  // Get the global transformation matrix for a given alignable volume
  //  identified by its symbolic name 'symname' by quering the TGeoManager

  if (!fgGeometry || !fgGeometry->IsClosed()) {
    AliErrorClass("No active geometry or geometry not yet closed!");
    return NULL;
  }

  TGeoPNEntry* pne = fgGeometry->GetAlignableEntry(symname);
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
Bool_t AliGeomManager::GetDeltaForBranch(Int_t index, TGeoHMatrix &inclusiveD)
{
  // The method sets the matrix passed as argument as the global delta
  // (for the volume referred by the unique index) including the displacements
  // of all parent volumes in the branch.
  //

  TGeoHMatrix go,invgo;
  go = *GetOrigGlobalMatrix(index);
  invgo = go.Inverse();
  inclusiveD = *GetMatrix(index);
  inclusiveD.Multiply(&invgo);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliGeomManager::GetDeltaForBranch(AliAlignObj& aao, TGeoHMatrix &inclusiveD)
{
  // The method sets the matrix passed as argument as the global delta
  // (for the volume referred by the alignment object) including the displacements
  // of all parent volumes in the brach.
  //
  Int_t index = aao.GetVolUID();
  if(!index){
    AliErrorClass("Either the alignment object or its index are not valid");
    return kFALSE;
  }
  return GetDeltaForBranch(index, inclusiveD);
}

//______________________________________________________________________
Bool_t AliGeomManager::GetOrigGlobalMatrix(const char* symname, TGeoHMatrix &m) 
{
  // Get the global transformation matrix (ideal geometry) for a given alignable volume
  // The alignable volume is identified by 'symname' which has to be either a valid symbolic
  // name, the query being performed after alignment, or a valid volume path if the query is
  // performed before alignment.
  //
  m.Clear();

  if (!fgGeometry || !fgGeometry->IsClosed()) {
    AliErrorClass("No active geometry or geometry not yet closed!");
    return kFALSE;
  }
  if (!fgGeometry->GetListOfPhysicalNodes()) {
    AliWarningClass("gGeoManager doesn't contain any aligned nodes!");
    if (!fgGeometry->cd(symname)) {
      AliErrorClass(Form("Volume path %s not valid!",symname));
      return kFALSE;
    }
    else {
      m = *fgGeometry->GetCurrentMatrix();
      return kTRUE;
    }
  }

  TGeoPNEntry* pne = fgGeometry->GetAlignableEntry(symname);
  const char* path = NULL;
  if(pne){
    m = *pne->GetGlobalOrig();
    return kTRUE;
  }else{
    AliWarningClass(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path=symname;
  }

  return GetOrigGlobalMatrixFromPath(path,m);
}

//_____________________________________________________________________________
Bool_t AliGeomManager::GetOrigGlobalMatrixFromPath(const char *path, TGeoHMatrix &m)
{
  // The method returns the global matrix for the volume identified by 
  // 'path' in the ideal detector geometry.
  // The output global matrix is stored in 'm'.
  // Returns kFALSE in case TGeo has not been initialized or the volume
  // path is not valid.
  //
  m.Clear();

  if (!fgGeometry || !fgGeometry->IsClosed()) {
    AliErrorClass("Can't get the original global matrix! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  if (!fgGeometry->CheckPath(path)) {
    AliErrorClass(Form("Volume path %s not valid!",path));
    return kFALSE;
  }

  TIter next(fgGeometry->GetListOfPhysicalNodes());
  fgGeometry->cd(path);

  while(fgGeometry->GetLevel()){

    TGeoPhysicalNode *physNode = NULL;
    next.Reset();
    TGeoNode *node = fgGeometry->GetCurrentNode();
    while ((physNode=(TGeoPhysicalNode*)next())) 
      if (physNode->GetNode() == node) break;

    TGeoMatrix *lm = NULL;
    if (physNode) {
      lm = physNode->GetOriginalMatrix();
      if (!lm) lm = node->GetMatrix();
    } else
      lm = node->GetMatrix();

    m.MultiplyLeft(lm);

    fgGeometry->CdUp();
  }

  return kTRUE;
}

//_____________________________________________________________________________
TGeoHMatrix* AliGeomManager::GetOrigGlobalMatrix(TGeoPNEntry * const pne)
{
  // The method returns global matrix for the ideal detector geometry
  // using the corresponding TGeoPNEntry as an input.
  // The returned pointer should be copied by the user, since its content could
  // be overwritten by a following call to the method.
  // In case of missing TGeoManager the method returns NULL.
  //
  if (!fgGeometry || !fgGeometry->IsClosed()) {
    AliErrorClass("Can't get the global matrix! gGeoManager doesn't exist or it is still opened!");
    return NULL;
  }

  return pne->GetGlobalOrig();
}

//______________________________________________________________________
TGeoHMatrix* AliGeomManager::GetOrigGlobalMatrix(Int_t index)
{
  // The method returns global matrix from the ideal detector geometry
  // for the volume identified by its index.
  // The returned pointer should be copied by the user, since its content could
  // be overwritten by a following call to the method.
  // In case of missing TGeoManager the method returns NULL.
  // If possible, the method uses the LUT of original ideal matrices
  // for fast access. The LUT is reset in case a
  // new geometry is loaded.
  //
  TGeoPNEntry* pne = GetPNEntry(index);
  return pne->GetGlobalOrig();
}

//______________________________________________________________________
Bool_t AliGeomManager::GetOrigTranslation(Int_t index, Double_t t[3]) 
{
  // Get the original translation vector (ideal geometry)
  // for a given module 'index' by quering the TGeoManager

  TGeoHMatrix *m = GetOrigGlobalMatrix(index);
  if (!m) return kFALSE;

  Double_t *trans = m->GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliGeomManager::GetOrigRotation(Int_t index, Double_t r[9]) 
{
  // Get the original rotation matrix (ideal geometry)
  // for a given module 'index' by quering the TGeoManager

  TGeoHMatrix *m = GetOrigGlobalMatrix(index);
  if (!m) return kFALSE;

  Double_t *rot = m->GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
const TGeoHMatrix* AliGeomManager::GetTracking2LocalMatrix(Int_t index)
{
  // Get the matrix which transforms from the tracking to the local RS
  // The method queries directly the TGeoPNEntry

  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  const TGeoHMatrix *m = pne->GetMatrix();
  if (!m)
    AliErrorClass(Form("TGeoPNEntry (%s) contains no tracking-to-local matrix !",pne->GetName()));

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
TGeoPNEntry* AliGeomManager::GetPNEntry(ELayerID layerId, Int_t modId)
{
  // Returns the TGeoPNEntry for a given layer
  // and module ID
  //

  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }

  return fgPNEntry[layerId-kFirstLayer][modId];
}

//_____________________________________________________________________________
void AliGeomManager::CheckOverlapsOverPNs(Double_t threshold)
{
  // Check for overlaps/extrusions on physical nodes only;
  // this overlap-checker is meant to be used to check overlaps/extrusions
  // originated by the application of alignment objects.
  //

  TObjArray* ovexlist = 0x0;

  AliInfoClass("********* Checking overlaps/extrusions over physical nodes only *********");
  TObjArray* pnList = gGeoManager->GetListOfPhysicalNodes();
  TGeoVolume* mvol = 0;
  TGeoPhysicalNode* pn;
  TObjArray* overlaps = new TObjArray(64);
  overlaps->SetOwner();

  TStopwatch timer2;
  timer2.Start();
  for(Int_t pni=0; pni<pnList->GetEntriesFast(); pni++){
    pn = (TGeoPhysicalNode*) pnList->UncheckedAt(pni);
    // checking the volume of the mother (go upper in the tree in case it is an assembly)
    Int_t levup=1;
    while(((TGeoVolume*)pn->GetVolume(pn->GetLevel()-levup))->IsAssembly()) levup++;
      //Printf("Going to upper level");
    mvol = pn->GetVolume(pn->GetLevel()-levup);
    if(!mvol->IsSelected()){
      AliInfoClass(Form("Checking overlaps for volume %s",mvol->GetName()));
      mvol->CheckOverlaps(threshold);
      ovexlist = gGeoManager->GetListOfOverlaps();
      TIter next(ovexlist);
      TGeoOverlap *ov;
      while ((ov=(TGeoOverlap*)next())) overlaps->Add(ov->Clone());
      mvol->SelectVolume();
    }
  }
  mvol->SelectVolume(kTRUE); // clears the list of selected volumes

  AliInfoClass(Form("Number of overlapping/extruding PNs: %d",overlaps->GetEntriesFast()));
  timer2.Stop();
  timer2.Print();

  TIter nextN(overlaps);
  TGeoOverlap *ovlp;
  while ((ovlp=(TGeoOverlap*)nextN())) ovlp->PrintInfo();

  overlaps->Delete();
  delete overlaps;
}

//_____________________________________________________________________________
Int_t AliGeomManager::GetNalignable(const char* module)
{
  // Get number of declared alignable volumes in current geometry
  // for the given detector "module" passed as a vaild detector name
  // if the detector name is invalid return -1
  
  // return the detector index corresponding to detector
  Int_t index = -1 ;
  for (index = 0; index < fgkNDetectors ; index++) {
    if ( strcmp(module, fgkDetectorName[index]) == 0 )
      break ;
  }
  if(index==fgkNDetectors) return -1;
  return fgNalignable[index];
}
  
//_____________________________________________________________________________
void AliGeomManager::InitNalignable()
{
  // Set number of declared alignable volumes for given detector in current geometry
  // by looping on the list of PNEntries
  //
  
  Int_t nAlE = gGeoManager->GetNAlignable(); // total number of alignable entries
  TGeoPNEntry *pne = 0;
  const char* detName;
  
  for (Int_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
    detName = fgkDetectorName[iDet];
    Int_t nAlDet = 0;
    
    for(Int_t iE = 0; iE < nAlE; iE++)
    {
      pne = gGeoManager->GetAlignableEntry(iE);
      TString pneName = pne->GetName();
      if(pneName.Contains(detName)) nAlDet++;
      if(!strcmp(detName,"GRP")) if(pneName.Contains("ABSO")  || pneName.Contains("DIPO") || 
				    pneName.Contains("FRAME") || pneName.Contains("PIPE") || 
				    pneName.Contains("SHIL")) nAlDet++;
    }
    fgNalignable[iDet] = nAlDet;
  }

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
 
  TObjArray alignObjArray;
  alignObjArray.Clear();  	
  alignObjArray.SetOwner(0);

  TString alObjsNotLoaded="";
  TString alObjsLoaded="";

  TString AlignDetsString(AlignDetsList);
  TObjArray *detsarr = AlignDetsString.Tokenize(' ');
  TIter iter(detsarr);
  TObjString *str = 0;
  
  while((str = (TObjString*) iter.Next())){
    TString det(str->String());
    AliDebugClass(5,Form("Loading alignment objs for %s",det.Data()));
    if(!LoadAlignObjsFromCDBSingleDet(det.Data(),alignObjArray)){
      alObjsNotLoaded += det.Data();
      alObjsNotLoaded += " ";
    } else {
      alObjsLoaded += det.Data();
      alObjsLoaded += " ";
    }
  }
  detsarr->Delete();
  delete detsarr;

  if(!alObjsLoaded.IsNull()) AliInfoClass(Form("Alignment objects loaded for: %s",
					       alObjsLoaded.Data()));
  if(!alObjsNotLoaded.IsNull())
      AliFatalClass(Form("Could not load alignment objects from OCDB for: %s",
						  alObjsNotLoaded.Data()));
 
  return ApplyAlignObjsToGeom(alignObjArray);
}

//_____________________________________________________________________________
Bool_t AliGeomManager::LoadAlignObjsFromCDBSingleDet(const char* detName, TObjArray& alignObjArray)
{
  // Adds the alignable objects found in the CDBEntry for the detector
  // passed as argument to the array of all alignment objects to be applyed
  // to geometry
  //
  // Fills array of single detector's alignable objects from CDB
  
  AliDebugClass(2, Form("Loading alignment objs for detector: %s",detName));
  
  AliCDBEntry *entry;
  	
  AliCDBPath path(detName,"Align","Data");
	
  entry=AliCDBManager::Instance()->Get(path.GetPath());
  if(!entry){ 
  	AliDebugClass(2,Form("Couldn't load alignment data for detector %s",detName));
	return kFALSE;
  }
  entry->SetOwner(1);
  TClonesArray *alignArray = (TClonesArray*) entry->GetObject();	
  alignArray->SetOwner(0);
  Int_t nAlObjs = alignArray->GetEntries();
  AliDebugClass(2,Form("Found %d alignment objects for %s",nAlObjs,detName));
  Int_t nAlVols = GetNalignable(detName);
  if(nAlObjs!=nAlVols) AliWarningClass(Form("%d alignment objects loaded for %s, which has %d alignable volumes",nAlObjs,detName,GetNalignable(detName)));

  AliAlignObj *alignObj=0;
  TIter iter(alignArray);
	
  // loop over align objects in detector
  while( ( alignObj=(AliAlignObj *) iter.Next() ) ){
  	alignObjArray.Add(alignObj);
  }
  // delete entry --- Don't delete, it is cached!
	
  AliDebugClass(2, Form("fAlignObjArray entries: %d",alignObjArray.GetEntries() ));
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliGeomManager::ApplyAlignObjsToGeom(TObjArray& alignObjArray, Bool_t ovlpcheck)
{
    // Read collection of alignment objects (AliAlignObj derived) saved
    // in the TClonesArray alObjArray and apply them to gGeoManager
    //
    alignObjArray.Sort();
    Int_t nvols = alignObjArray.GetEntriesFast();

    Bool_t flag = kTRUE;

    for(Int_t j=0; j<nvols; j++)
    {
	AliAlignObj* alobj = (AliAlignObj*) alignObjArray.UncheckedAt(j);
	if(!alobj->ApplyToGeometry(ovlpcheck))
	{
	    flag = kFALSE;
	    AliDebugClass(5,Form("Error applying alignment object for volume %s !",alobj->GetSymName()));
	}else{
	    AliDebugClass(5,Form("Alignment object for volume %s applied successfully",alobj->GetSymName()));
	}

    }

    if (AliDebugLevelClass() > 5) {
	fgGeometry->CheckOverlaps(0.001);
	TObjArray* ovexlist = fgGeometry->GetListOfOverlaps();
	if(ovexlist->GetEntriesFast()){  
	    AliErrorClass("The application of alignment objects to the geometry caused huge overlaps/extrusions!");
	    fgGeometry->PrintOverlaps();
	}
    }

    // Update the TGeoPhysicalNodes
    fgGeometry->RefreshPhysicalNodes();

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

  TClonesArray* alignObjArray = ((TClonesArray*) inFile->Get(clArrayName));
  inFile->Close();
  if (!alignObjArray) {
    AliErrorClass(Form("Could not get array (%s) from file (%s) !",clArrayName,fileName));
    return kFALSE;
  }

  return ApplyAlignObjsToGeom(*alignObjArray);

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
  TClonesArray* alignObjArray = ((TClonesArray*) entry->GetObject());

  return ApplyAlignObjsToGeom(*alignObjArray);

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
  TClonesArray* alignObjArray = ((TClonesArray*) entry->GetObject());

  return ApplyAlignObjsToGeom(*alignObjArray);
}

//_____________________________________________________________________________
void AliGeomManager::ResetPNEntriesLUT()
{
  // cleans static arrays containing the information on currently loaded geometry
  //
  for (Int_t iLayer = 0; iLayer < (kLastLayer - kFirstLayer); iLayer++){
    if (!fgPNEntry[iLayer]) continue;
    for (Int_t modnum=0; modnum<fgLayerSize[iLayer]; modnum++) fgPNEntry[iLayer][modnum] = 0;
  }
  //
}
  
