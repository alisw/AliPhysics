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

///////////////////////////////////////////////////////////////////////////
//    AliITSMFTGeomTGeo is a simple interface class to TGeoManager         //
//    It is used in the simulation and reconstruction in order to        //
//    query the TGeo geometry                                            //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#include <TDatime.h>
#include <TMath.h>
#include <TSystem.h>
#include <TObjArray.h>

#include "AliITSMFTGeomTGeo.h"
#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliITSMFTSegmentationPix.h"
using namespace TMath;

ClassImp(AliITSMFTGeomTGeo)

UInt_t AliITSMFTGeomTGeo::fgUIDShift = 16;                // bit shift to go from mod.id to modUUID for TGeo
TString AliITSMFTGeomTGeo::fgITSVolName = "ITSV";
TString AliITSMFTGeomTGeo::fgITSLrName  = "ITSMFTLayer";
TString AliITSMFTGeomTGeo::fgITSStaveName = "ITSMFTStave";
TString AliITSMFTGeomTGeo::fgITSHalfStaveName = "ITSMFTHalfStave";
TString AliITSMFTGeomTGeo::fgITSModuleName = "ITSMFTModule";
TString AliITSMFTGeomTGeo::fgITSChipName = "ITSMFTChip";
TString AliITSMFTGeomTGeo::fgITSSensName = "ITSMFTSensor";
TString AliITSMFTGeomTGeo::fgITSWrapVolName = "ITSMFTWrapVol";
TString AliITSMFTGeomTGeo::fgITSChipTypeName[AliITSMFTAux::kNChipTypes] = {"Pix"};
//
TString AliITSMFTGeomTGeo::fgITSsegmFileName = "itsSegmentations.root";

//______________________________________________________________________
AliITSMFTGeomTGeo::AliITSMFTGeomTGeo(Bool_t /*build*/, Bool_t /*loadSegm*/)
  :fNLayers(0)
  ,fNChips(0)
  ,fNStaves(0)
  ,fNHalfStaves(0)
  ,fNModules(0)
  ,fNChipsPerModule(0)
  ,fNChipRowsPerModule(0)
  ,fNChipsPerHalfStave(0)
  ,fNChipsPerStave(0)
  ,fNChipsPerLayer(0)
  ,fLrChipType(0)
  ,fLastChipIndex(0)
  ,fMatSens(0)
  ,fMatT2L(0)
  ,fSegm(0)
{
  // default c-tor
  for (int i=AliITSMFTAux::kMaxLayers;i--;) fLr2Wrapper[i] = -1;
}

//______________________________________________________________________
AliITSMFTGeomTGeo::AliITSMFTGeomTGeo(const AliITSMFTGeomTGeo &src)
  :TObject(src)
  ,fNLayers(src.fNLayers)
  ,fNChips(src.fNChips)
  ,fNStaves(0)
  ,fNHalfStaves(0)
  ,fNModules(0)
  ,fNChipsPerModule(0)
  ,fNChipRowsPerModule(0)
  ,fNChipsPerHalfStave(0)
  ,fNChipsPerStave(0)
  ,fNChipsPerLayer(0)
  ,fLrChipType(0)
  ,fLastChipIndex(0)
  ,fMatSens(0)
  ,fMatT2L(0)
  ,fSegm(0)
{
  // copy c-tor
  if (fNLayers) {
    fNStaves   = new Int_t[fNLayers];
    fNChipsPerModule = new Int_t[fNLayers];
    fNChipRowsPerModule = new Int_t[fNLayers];
    fLrChipType  = new Int_t[fNLayers];
    fLastChipIndex   = new Int_t[fNLayers];
    fNChipsPerHalfStave = new Int_t[fNLayers];
    fNChipsPerStave = new Int_t[fNLayers];
    fNChipsPerLayer = new Int_t[fNLayers];
    //
    for (int i=fNLayers;i--;) {
      fNStaves[i] = src.fNStaves[i];
      fNHalfStaves[i] = src.fNHalfStaves[i];
      fNModules[i] = src.fNModules[i];
      fNChipsPerModule[i] = src.fNChipsPerModule[i];
      fNChipRowsPerModule[i] = src.fNChipRowsPerModule[i];
      fNChipsPerHalfStave[i] = src.fNChipsPerHalfStave[i];
      fNChipsPerStave[i] = src.fNChipsPerStave[i];
      fNChipsPerLayer[i] = src.fNChipsPerLayer[i];
      fLrChipType[i]  = src.fLrChipType[i];
      fLastChipIndex[i] = src.fLastChipIndex[i];
    }
    if (src.fMatSens) {
      fMatSens = new TObjArray(fNChips);
      fMatSens->SetOwner(kTRUE);
      for (int i=0;i<fNChips;i++) {
	const TGeoHMatrix* mat = (TGeoHMatrix*)src.fMatSens->At(i);
	fMatSens->AddAt(new TGeoHMatrix(*mat),i);
      }
    }
    if (src.fMatT2L) {
      fMatT2L = new TObjArray(fNChips);
      fMatT2L->SetOwner(kTRUE);
      for (int i=0;i<fNChips;i++) {
	const TGeoHMatrix* mat =(TGeoHMatrix*) src.fMatT2L->At(i);
	fMatT2L->AddAt(new TGeoHMatrix(*mat),i);
      }
    }
    if (src.fSegm) {
      int sz = src.fSegm->GetEntriesFast();
      fSegm = new TObjArray(sz);
      fSegm->SetOwner(kTRUE);
      for (int i=0;i<sz;i++) {
	AliITSMFTSegmentationPix* sg = (AliITSMFTSegmentationPix*)src.fSegm->UncheckedAt(i);
	if (!sg) continue;
	fSegm->AddAt(sg->Clone(),i);
      }
    }
  }
  for (int i=AliITSMFTAux::kMaxLayers;i--;) fLr2Wrapper[i] = src.fLr2Wrapper[i];
}

//______________________________________________________________________
AliITSMFTGeomTGeo::~AliITSMFTGeomTGeo()
{
  //d-tor
  delete[] fNStaves;
  delete[] fNHalfStaves;
  delete[] fNModules;
  delete[] fLrChipType;
  delete[] fNChipsPerModule;
  delete[] fNChipRowsPerModule;
  delete[] fNChipsPerHalfStave;
  delete[] fNChipsPerStave;
  delete[] fNChipsPerLayer;
  delete[] fLastChipIndex;
  delete fMatT2L;
  delete fMatSens;
  delete fSegm;
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetChipIndex(Int_t lay,Int_t sta,Int_t chipInStave) const
{
  // This routine computes the chip index number from the layer,
  // stave, and chip number in stave. 
  // Inputs:
  //    Int_t lay  The layer number. Starting from 0.
  //    Int_t sta  The stave number. Starting from 0
  //    Int_t chipInStave  The chip number in the stave. Starting from 0
  //
  return GetFirstChipIndex(lay) + fNChipsPerStave[lay]*sta + chipInStave;
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetChipIndex(Int_t lay,Int_t sta, Int_t substa, Int_t chipInSStave) const
{
  // This routine computes the chip index number from the layer,
  // stave, substave and chip number in substave. 
  // Inputs:
  //    Int_t lay  The layer number. Starting from 0.
  //    Int_t sta  The stave number. Starting from 0
  //    Int_t substa  The substave number. Starting from 0
  //    Int_t chipInSStave  The chip number in the sub stave. Starting from 0
  //
  int n = GetFirstChipIndex(lay) + fNChipsPerStave[lay]*sta + chipInSStave;
  if (fNHalfStaves[lay] && substa>0) n += fNChipsPerHalfStave[lay]*substa;
  return n;
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetChipIndex(Int_t lay,Int_t sta, Int_t substa, Int_t md, Int_t chipInMod) const
{
  // This routine computes the chip index number from the layer,
  // stave, substave module and chip number in module. 
  // Inputs:
  //    Int_t lay  The layer number. Starting from 0.
  //    Int_t sta  The stave number. Starting from 0
  //    Int_t substa  The substave number. Starting from 0
  //    Int_t module  The module number ...
  //    Int_t chipInSStave  The chip number in the module. Starting from 0
  //
  int n = GetFirstChipIndex(lay) + fNChipsPerStave[lay]*sta + chipInMod;
  if (fNHalfStaves[lay] && substa>0) n += fNChipsPerHalfStave[lay]*substa;
  if (fNModules[lay] && md>0)       n += fNChipsPerModule[lay]*md;
  return n;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetLayer(Int_t index,Int_t &lay,Int_t &indexInLr)  const
{
  // This routine computes the layer number a
  // given the chip index. The 
  // Inputs:
  //     Int_t index  The chip index number, starting from zero.
  // Outputs:
  //     Int_t indexInLr The chip index inside a layer, starting from zero.
  //     Int_t lay    The layer number. Starting from 0.
  //
  lay = GetLayer(index);
  indexInLr = index - GetFirstChipIndex(lay);
  return kTRUE;
  //
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetLayer(Int_t index) const
{
  // Get chip layer, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  return lay;
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetStave(Int_t index) const
{
  // Get chip stave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index/fNChipsPerStave[lay];
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetHalfStave(Int_t index) const
{
  // Get chip substave id in stave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  if (fNHalfStaves[lay]<0) return -1;
  index -= GetFirstChipIndex(lay);
  index %= fNChipsPerStave[lay];
  return index/fNChipsPerHalfStave[lay];
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetModule(Int_t index) const
{
  // Get chip module id in substave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  if (fNModules[lay]<0) return 0;
  index -= GetFirstChipIndex(lay);
  index %= fNChipsPerStave[lay];
  if (fNHalfStaves[lay]) index %= fNChipsPerHalfStave[lay];
  return index/fNChipsPerModule[lay];
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetChipIdInLayer(Int_t index) const
{
  // Get chip number within layer, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index;
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetChipIdInStave(Int_t index) const
{
  // Get chip number within stave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index%fNChipsPerStave[lay];
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetChipIdInHalfStave(Int_t index) const
{
  // Get chip number within stave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index%fNChipsPerHalfStave[lay];
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::GetChipIdInModule(Int_t index) const
{
  // Get chip number within module, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index%fNChipsPerModule[lay];
}

//______________________________________________________________________
const char* AliITSMFTGeomTGeo::GetSymName(Int_t index)  const
{
  // Get the TGeoPNEntry symbolic name
  // for a given chip identified by 'index'
  //
  Int_t lay, index2;
  if (!GetLayer(index,lay,index2)) return NULL;
  // return AliGeomManager::SymName((AliGeomManager::ELayerID)((lay-1)+AliGeomManager::kSPD1),index2);
  // RS: this is not optimal, but we cannod access directly AliGeomManager, since the latter has hardwired layers 
  //  TGeoPNEntry* pne = gGeoManager->GetAlignableEntryByUID( AliGeomManager::LayerToVolUID(lay+1,index2) );
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntryByUID( ChipVolUID(index) );
  if (!pne) {
    AliError(Form("Failed to find alignable entry with index %d: (Lr%d Chip:%d) !",index,lay,index2));
    return NULL;
  }
  return pne->GetName();
}

//______________________________________________________________________
const char* AliITSMFTGeomTGeo::ComposeSymNameITS()
{
  // sym name of the layer
  return "ITS";
}

//______________________________________________________________________
const char* AliITSMFTGeomTGeo::ComposeSymNameLayer(Int_t lr)
{
  // sym name of the layer
  return Form("%s/%s%d",ComposeSymNameITS(),GetITSLayerPattern(),lr);
}

//______________________________________________________________________
const char* AliITSMFTGeomTGeo::ComposeSymNameStave(Int_t lr, Int_t stave)
{
  // sym name of the stave at given layer
  return Form("%s/%s%d",ComposeSymNameLayer(lr),GetITSStavePattern(),stave);
}

//______________________________________________________________________
const char* AliITSMFTGeomTGeo::ComposeSymNameHalfStave(Int_t lr, Int_t stave, Int_t substave)
{
  // sym name of the stave at given layer
  return substave>=0 ? 
    Form("%s/%s%d",ComposeSymNameStave(lr,stave),GetITSHalfStavePattern(),substave) :
    ComposeSymNameStave(lr,stave);
}

//______________________________________________________________________
const char* AliITSMFTGeomTGeo::ComposeSymNameModule(Int_t lr, Int_t stave, Int_t substave, Int_t mod)
{
  // sym name of the substave at given layer/stave
  return mod>=0 ? 
    Form("%s/%s%d",ComposeSymNameHalfStave(lr,stave,substave),GetITSModulePattern(),mod) :
    ComposeSymNameHalfStave(lr,stave,substave);    
}

//______________________________________________________________________
const char* AliITSMFTGeomTGeo::ComposeSymNameChip(Int_t lr, Int_t sta, Int_t substave, Int_t mod, Int_t chip)
{
  // sym name of the chip in the given layer/stave/substave/module
  return Form("%s/%s%d",ComposeSymNameModule(lr,sta,substave,mod),GetITSChipPattern(),chip);
}

//______________________________________________________________________
TGeoHMatrix* AliITSMFTGeomTGeo::GetMatrix(Int_t index)  const
{
  // Get the transformation matrix for a given chip 'index'
  // by quering the TGeoManager
  static TGeoHMatrix matTmp;
  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if (pnode) return pnode->GetMatrix();

  const char* path = pne->GetTitle();
  gGeoManager->PushPath(); // Preserve the modeler state.
  if (!gGeoManager->cd(path)) {
    gGeoManager->PopPath();
    AliError(Form("Volume path %s not valid!",path));
    return NULL;
  }
  matTmp = *gGeoManager->GetCurrentMatrix();
  gGeoManager->PopPath();
  return &matTmp;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetTranslation(Int_t index, Double_t t[3])  const
{
  // Get the translation vector for a given chip 'index'
  // by quering the TGeoManager
  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *trans = m->GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetRotation(Int_t index, Double_t r[9])  const
{
  // Get the rotation matrix for a given chip 'index'
  // by quering the TGeoManager
  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *rot = m->GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetOrigMatrix(Int_t index, TGeoHMatrix &m) const
{
  // Get the original (ideal geometry) TGeo matrix for
  // a given chip identified by 'index'.
  // The method is slow, so it should be used
  // with great care.
  m.Clear();

  const char *symname = GetSymName(index);
  if (!symname) return kFALSE;

  return AliGeomManager::GetOrigGlobalMatrix(symname,m);
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetOrigTranslation(Int_t index, Double_t t[3])  const
{
  // Get the original translation vector (ideal geometry)
  // for a given chip 'index' by quering the TGeoManager
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *trans = m.GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetOrigRotation(Int_t index, Double_t r[9])  const
{
  // Get the original rotation matrix (ideal geometry)
  // for a given chip 'index' by quering the TGeoManager
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *rot = m.GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
TGeoHMatrix* AliITSMFTGeomTGeo::ExtractMatrixT2L(Int_t index) const
{
  // Get the matrix which transforms from the tracking to local r.s.
  // The method queries directly the TGeoPNEntry
  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  TGeoHMatrix *m = (TGeoHMatrix*) pne->GetMatrix();
  if (!m) AliError(Form("TGeoPNEntry (%s) contains no matrix !",pne->GetName()));

  return m;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetTrackingMatrix(Int_t index, TGeoHMatrix &m)
{
  // Get the matrix which transforms from the tracking r.s. to
  // the global one.
  // Returns kFALSE in case of error.
  m.Clear();

  TGeoHMatrix *m1 = GetMatrix(index);
  if (!m1) return kFALSE;

  const TGeoHMatrix *m2 = GetMatrixT2L(index);
  if (!m2) return kFALSE;

  m = *m1;
  m.Multiply(m2);

  return kTRUE;
}

//______________________________________________________________________
TGeoPNEntry* AliITSMFTGeomTGeo::GetPNEntry(Int_t index) const
{
  // Get a pointer to the TGeoPNEntry of a chip
  // identified by 'index'
  // Returns NULL in case of invalid index,
  // missing TGeoManager or invalid symbolic name
  //
  if (index >= fNChips) {
    AliError(Form("Invalid ITS chip index: %d (0 -> %d) !",index,fNChips));
    return NULL;
  }
  
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the matrix! gGeoManager doesn't exist or it is still opened!");
    return NULL;
  }
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntryByUID( ChipVolUID(index) );
  //  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(GetSymName(index));
  if (!pne) AliError(Form("The index %d does not correspond to a physical entry!",index));
  //
  return pne;
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::ExtractNumberOfLayers()
{
  // Determines the number of layers in the Upgrade Geometry
  //
  Int_t numberOfLayers = 0;
  //
  TGeoVolume *itsV = gGeoManager->GetVolume(GetITSVolPattern());
  if (!itsV) AliFatal(Form("ITS volume %s is not in the geometry",GetITSVolPattern()));
  SetUIDShift(itsV->GetUniqueID());
  //
  // Loop on all ITSV nodes, count Layer volumes by checking names
  // Build on the fly layer - wrapper correspondence
  TObjArray* nodes = itsV->GetNodes();
  Int_t nNodes = nodes->GetEntriesFast();
  //
  for (Int_t j=0; j<nNodes; j++) {
    int lrID = -1;
    TGeoNode* nd = (TGeoNode*)nodes->At(j);
    const char* name = nd->GetName();
    if (strstr(name,GetITSLayerPattern())) {
      numberOfLayers++;
      if ( (lrID=ExtractVolumeCopy(name,AliITSMFTGeomTGeo::GetITSLayerPattern()))<0 ) {
	AliFatal(Form("Failed to extract layer ID from the %s",name));
	exit(1);
      }
      //
      fLr2Wrapper[lrID] = -1; // not wrapped
    }
    else if (strstr(name,GetITSWrapVolPattern())) { // this is a wrapper volume, may cointain layers
      int wrID = -1;
      if ( (wrID=ExtractVolumeCopy(name,AliITSMFTGeomTGeo::GetITSWrapVolPattern()))<0 ) {
	AliFatal(Form("Failed to extract wrapper ID from the %s",name));
	exit(1);
      }
      //
      TObjArray* nodesW = nd->GetNodes();
      int nNodesW = nodesW->GetEntriesFast();
      for (Int_t jw=0; jw<nNodesW; jw++) {
	TGeoNode* ndW = (TGeoNode*)nodesW->At(jw);
	if (strstr(ndW->GetName(),GetITSLayerPattern())) {
	  if ( (lrID=ExtractVolumeCopy(ndW->GetName(),AliITSMFTGeomTGeo::GetITSLayerPattern()))<0 ) {
	    AliFatal(Form("Failed to extract layer ID from the %s",name));
	    exit(1);
	  }
	  numberOfLayers++;
	  fLr2Wrapper[lrID] = wrID;
	}
      }
    }
  }
  //  
  return numberOfLayers;
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::ExtractNumberOfStaves(Int_t lay) const
{
  // Determines the number of layers in the Upgrade Geometry
  //
  // Inputs:
  //   lay: layer number, starting from 0
  //
  // MS
  Int_t numberOfStaves = 0;
  char laynam[30];
  snprintf(laynam, 30, "%s%d",GetITSLayerPattern(),lay);
  TGeoVolume* volLr = gGeoManager->GetVolume(laynam);
  if (!volLr) { AliFatal(Form("can't find %s volume",laynam)); return -1; }
  //
  // Loop on all layer nodes, count Stave volumes by checking names
  Int_t nNodes = volLr->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) {
    //    AliInfo(Form("L%d %d of %d %s %s -> %d",lay,j,nNodes,volLr->GetNodes()->At(j)->GetName(),GetITSStavePattern(),numberOfStaves));
    if (strstr(volLr->GetNodes()->At(j)->GetName(),GetITSStavePattern())) numberOfStaves++;
  }
  //
  return numberOfStaves;
  //
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::ExtractNumberOfHalfStaves(Int_t lay) const
{
  // Determines the number of substaves in the stave of the layer
  //
  // Inputs:
  //   lay: layer number, starting from 0
  //
  // MS
  if (fgITSHalfStaveName.IsNull()) return 0; // for the setup w/o substave defined the stave and the substave is the same thing
  Int_t nSS = 0;
  char stavnam[30];
  snprintf(stavnam, 30, "%s%d", GetITSStavePattern(),lay);
  TGeoVolume* volLd = gGeoManager->GetVolume(stavnam);
  if (!volLd) AliFatal(Form("can't find %s volume",stavnam));
  //
  // Loop on all stave nodes, count Chip volumes by checking names
  Int_t nNodes = volLd->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(volLd->GetNodes()->At(j)->GetName(),GetITSHalfStavePattern())) nSS++;
  //
  return nSS;
  //
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::ExtractNumberOfModules(Int_t lay) const
{
  // Determines the number of modules in substave in the stave of the layer
  //
  // Inputs:
  //   lay: layer number, starting from 0
  //
  // for the setup w/o modules defined the module and the stave or the substave is the same thing
  if (fgITSModuleName.IsNull()) return 0;
  char stavnam[30];
  TGeoVolume* volLd = 0;
  if (!fgITSHalfStaveName.IsNull()) {
    snprintf(stavnam, 30, "%s%d", GetITSHalfStavePattern(),lay); 
    volLd = gGeoManager->GetVolume(stavnam);
  }
  if (!volLd) { // no substaves, check staves
    snprintf(stavnam, 30, "%s%d", GetITSStavePattern(),lay); 
    volLd = gGeoManager->GetVolume(stavnam);
  }
  if (!volLd) return 0;
  Int_t nMod = 0;
  //
  // Loop on all substave nodes, count module volumes by checking names
  Int_t nNodes = volLd->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(volLd->GetNodes()->At(j)->GetName(),GetITSModulePattern())) nMod++;
  //
  return nMod;
  //
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::ExtractNChipsPerModule(Int_t lay, int &nrow)  const
{
  // Determines the number of chips per module on the (sub)stave in the Upgrade Geometry
  // Also extract the layout: span of module centers in Z and X
  // Inputs:
  //   lay: layer number from 0
  // MS
  Int_t numberOfChips = 0;
  char stavnam[30];
  TGeoVolume* volLd = 0;
  if (!fgITSModuleName.IsNull()) {
    snprintf(stavnam, 30, "%s%d", GetITSModulePattern(),lay); 
    volLd = gGeoManager->GetVolume(stavnam);
  }
  if (!volLd) { // no modules on this layer, check substaves
    if (!fgITSHalfStaveName.IsNull()) {
      snprintf(stavnam, 30, "%s%d", GetITSHalfStavePattern(),lay); 
      volLd = gGeoManager->GetVolume(stavnam);
    }
  }
  if (!volLd) { // no substaves on this layer, check staves
    snprintf(stavnam, 30, "%s%d", GetITSStavePattern(),lay);
    volLd = gGeoManager->GetVolume(stavnam);
  }
  if (!volLd) AliFatal(Form("can't find volume containing chips on layer %d",lay));
  //
  // Loop on all stave nodes, count Chip volumes by checking names
  Int_t nNodes = volLd->GetNodes()->GetEntries();
  //
  double xmin=1e9,xmax=-1e9, zmin=1e9,zmax=-1e9;
  double lab[3],loc[3]={0,0,0};
  double dx=-1,dz=-1;
  for (Int_t j=0; j<nNodes; j++) {
    //    AliInfo(Form("L%d %d of %d %s %s -> %d",lay,j,nNodes,volLd->GetNodes()->At(j)->GetName(),GetITSChipPattern(),numberOfChips));
    TGeoNodeMatrix* node = (TGeoNodeMatrix*)volLd->GetNodes()->At(j);
    if (!strstr(node->GetName(),GetITSChipPattern())) continue;
    node->LocalToMaster(loc,lab);
    if (lab[0]>xmax) xmax=lab[0];
    if (lab[0]<xmin) xmin=lab[0];    
    if (lab[2]>zmax) zmax=lab[2];
    if (lab[2]<zmin) zmin=lab[2];    
    //
    numberOfChips++;
    //
    if (dx<0) {
      TGeoShape* chShape = node->GetVolume()->GetShape();
      TGeoBBox* bbox = dynamic_cast<TGeoBBox*>(chShape);
      if (!bbox) {
	AliFatal(Form("Chip %s volume is of unprocessed shape %s",node->GetName(),chShape->IsA()->GetName()));
      }
      else {
	dx = 2*bbox->GetDX();
	dz = 2*bbox->GetDZ();
      }
    }
  }
  //
  double spanX = xmax-xmin;
  double spanZ = zmax-zmin;  
  nrow = TMath::Nint(spanX/dx + 1);
  int ncol = TMath::Nint(spanZ/dz + 1);
  if (nrow*ncol != numberOfChips) 
    AliError(Form("Inconsistency between Nchips=%d and Nrow*Ncol=%d*%d->%d\n"
		  "Extracted chip dimensions (x,z): %.4f %.4f, Module Span: %.4f %.4f",
		  numberOfChips,nrow,ncol,nrow*ncol,
		  dx,dz,spanX,spanZ));
  return numberOfChips;
  //
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::ExtractLayerChipType(Int_t lay)  const
{
  // Determines the layer detector type the Upgrade Geometry
  //
  // Inputs:
  //   lay: layer number from 0
  // Outputs:
  //   none
  // Return:
  //   detector type id for the layer
  // MS
  char stavnam[30];
  snprintf(stavnam, 30, "%s%d", GetITSLayerPattern(),lay);
  TGeoVolume* volLd = gGeoManager->GetVolume(stavnam);
  if (!volLd) {AliFatal(Form("can't find %s volume",stavnam)); return -1;}
  //
  return volLd->GetUniqueID();
  //
}

//______________________________________________________________________
void AliITSMFTGeomTGeo::FetchMatrices()
{
  // store pointer on often used matrices for faster access
  if (!gGeoManager) AliFatal("Geometry is not loaded");
  fMatSens = new TObjArray(fNChips);
  fMatSens->SetOwner(kTRUE);
  for (int i=0;i<fNChips;i++) fMatSens->AddAt(new TGeoHMatrix(*ExtractMatrixSens(i)),i);
  CreateT2LMatrices();
}

//______________________________________________________________________
void AliITSMFTGeomTGeo::CreateT2LMatrices()
{
  // create tracking to local (Sensor!) matrices
  fMatT2L  = new TObjArray(fNChips);  
  fMatT2L->SetOwner(kTRUE);
  TGeoHMatrix matLtoT;
  double locA[3]={-100,0,0},locB[3]={100,0,0},gloA[3],gloB[3];
  for (int isn=0;isn<fNChips;isn++) {
    const TGeoHMatrix* matSens = GetMatrixSens(isn);
    if (!matSens) {AliFatal(Form("Failed to get matrix for sensor %d",isn)); return;}
    matSens->LocalToMaster(locA,gloA);
    matSens->LocalToMaster(locB,gloB);
    double dx = gloB[0]-gloA[0];
    double dy = gloB[1]-gloA[1];
    double t = (gloB[0]*dx+gloB[1]*dy)/(dx*dx+dy*dy),x=gloB[0]-dx*t,y=gloB[1]-dy*t;
    TGeoHMatrix* t2l = new TGeoHMatrix();
    t2l->RotateZ(ATan2(y,x)*RadToDeg()); // rotate in direction of normal to the sensor plane
    t2l->SetDx(x);
    t2l->SetDy(y);
    t2l->MultiplyLeft(&matSens->Inverse());
    fMatT2L->AddAt(t2l,isn);
    /*
    const double *gtrans = matSens->GetTranslation();
    memcpy(&rotMatrix[0], matSens->GetRotationMatrix(), 9*sizeof(Double_t));
    Double_t al = -ATan2(rotMatrix[1],rotMatrix[0]);
    Double_t rSens = Sqrt(gtrans[0]*gtrans[0] + gtrans[1]*gtrans[1]);
    Double_t tanAl = ATan2(gtrans[1],gtrans[0]) - Pi()/2; //angle of tangent
    Double_t alTr = tanAl - al;
    //
    // The X axis of tracking frame must always look outward
    loc[1] = rSens/2;
    matSens->LocalToMaster(loc,glo);
    double rPos = Sqrt(glo[0]*glo[0] + glo[1]*glo[1]);
    Bool_t rotOutward = rPos>rSens ? kFALSE : kTRUE;
    //
    // Transformation matrix
    matLtoT.Clear();
    matLtoT.SetDx(-rSens*Sin(alTr)); // translation
    matLtoT.SetDy(0.);
    matLtoT.SetDz(gtrans[2]);
    // Rotation matrix
    rotMatrix[0]= 0;  rotMatrix[1]= 1;  rotMatrix[2]= 0; // + rotation
    rotMatrix[3]=-1;  rotMatrix[4]= 0;  rotMatrix[5]= 0;
    rotMatrix[6]= 0;  rotMatrix[7]= 0;  rotMatrix[8]= 1;
    //
    TGeoRotation rot;
    rot.SetMatrix(rotMatrix);
    matLtoT.MultiplyLeft(&rot);
    if (rotOutward) matLtoT.RotateZ(180.);
    // Inverse transformation Matrix
    fMatT2L->AddAt(new TGeoHMatrix(matLtoT.Inverse()),isn);
    */
  }
  //
}

//______________________________________________________________________
Int_t AliITSMFTGeomTGeo::ExtractVolumeCopy(const char* name, const char* prefix) const
{
  // extract Number following the prefix in the name string
  TString nms = name;
  if (!nms.BeginsWith(prefix)) return -1;
  nms.Remove(0,strlen(prefix));
  if (!isdigit(nms.Data()[0])) return -1;
  return nms.Atoi();
  //
}
