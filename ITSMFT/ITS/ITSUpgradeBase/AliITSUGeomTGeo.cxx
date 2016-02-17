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
//    AliITSUGeomTGeo is a simple interface class to TGeoManager       //
//    It is used in the simulation and reconstruction in order to        //
//    query the TGeo ITS geometry                                        //
//                                                                       //
//    author - cvetan.cheshkov@cern.ch                                   //
//    15/02/2007                                                         //
//    adapted to ITSupg 18/07/2012 - ruben.shahoyan@cern.ch              //
//                                                                       //
//    ATTENTION: In opposite to ols AliITSgeomTGeo, all indices start    //
//    from 0, not from 1!!!                                              //
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

#include "AliITSUGeomTGeo.h"
#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliITSMFTSegmentationPix.h"
using namespace TMath;

ClassImp(AliITSUGeomTGeo)

// bit shift to go from mod.id to modUUID for TGeo
UInt_t AliITSUGeomTGeo::fgUIDShift = 16;
//

TString AliITSUGeomTGeo::fgITSVolName = "ITSV";
TString AliITSUGeomTGeo::fgITSLrName  = "ITSULayer";
TString AliITSUGeomTGeo::fgITSStaveName = "ITSUStave";
TString AliITSUGeomTGeo::fgITSHalfStaveName = "ITSUHalfStave";
TString AliITSUGeomTGeo::fgITSModuleName = "ITSUModule";
TString AliITSUGeomTGeo::fgITSChipName = "ITSUChip";
TString AliITSUGeomTGeo::fgITSSensName = "ITSUSensor";
TString AliITSUGeomTGeo::fgITSWrapVolName = "ITSUWrapVol";
TString AliITSUGeomTGeo::fgITSChipTypeName[AliITSMFTAux::kNChipTypes]={"Pix"};
//
TString AliITSUGeomTGeo::fgITSsegmFileName = "itsSegmentations.root";

//______________________________________________________________________
AliITSUGeomTGeo::AliITSUGeomTGeo(Bool_t build, Bool_t loadSegm)
  :AliITSMFTGeomTGeo()
  ,fVersion(kITSVNA)
  ,fNLayers(0)
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
{
  // default c-tor
  for (Int_t i=AliITSMFTAux::kMaxLayers;i--;) fLr2Wrapper[i] = -1;
  if (build) BuildITS(loadSegm);
}

//______________________________________________________________________
AliITSUGeomTGeo::AliITSUGeomTGeo(const AliITSUGeomTGeo &src)
  :AliITSMFTGeomTGeo(src)
  ,fVersion(src.fVersion)
  ,fNLayers(src.fNLayers)
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
  }
  for (int i=AliITSMFTAux::kMaxLayers;i--;) fLr2Wrapper[i] = src.fLr2Wrapper[i];
}

//______________________________________________________________________
AliITSUGeomTGeo::~AliITSUGeomTGeo()
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
}

//______________________________________________________________________
void AliITSUGeomTGeo::BuildITS(Bool_t loadSegm)
{
  // exract upg ITS parameters from TGeo
  if (fVersion!=kITSVNA) {AliWarning("Already built"); return;} // already initialized
  if (!gGeoManager) AliFatal("Geometry is not loaded");
  fNLayers    = ExtractNumberOfLayers();
  if (!fNLayers) return;
  //
  fNStaves         = new Int_t[fNLayers];
  fNHalfStaves     = new Int_t[fNLayers];
  fNModules        = new Int_t[fNLayers];
  fNChipsPerModule = new Int_t[fNLayers];
  fNChipRowsPerModule = new Int_t[fNLayers];
  fNChipsPerHalfStave = new Int_t[fNLayers];
  fNChipsPerStave  = new Int_t[fNLayers];
  fNChipsPerLayer  = new Int_t[fNLayers];
  fLrChipType      = new Int_t[fNLayers];
  fLastChipIndex   = new Int_t[fNLayers];
  fNChips = 0;
  
  for (int i=0;i<fNLayers;i++) {
    fLrChipType[i]      = ExtractLayerChipType(i);
    fNStaves[i]         = ExtractNumberOfStaves(i);
    fNHalfStaves[i]     = ExtractNumberOfHalfStaves(i);
    fNModules[i]        = ExtractNumberOfModules(i);
    fNChipsPerModule[i] = ExtractNChipsPerModule(i,fNChipRowsPerModule[i]);
    fNChipsPerHalfStave[i]= fNChipsPerModule[i]*Max(1,fNModules[i]);
    fNChipsPerStave[i]    = fNChipsPerHalfStave[i]*Max(1,fNHalfStaves[i]);
    fNChipsPerLayer[i]    = fNChipsPerStave[i]*fNStaves[i];
    fNChips               += fNChipsPerLayer[i];
    fLastChipIndex[i]     = fNChips-1;
  }
  //
  FetchMatrices();
  fVersion = kITSVUpg;
  //
  if (loadSegm) {  // fetch segmentations
    fSegm = new TObjArray();
    AliITSMFTSegmentationPix::LoadSegmentations(fSegm,GetITSsegmentationFileName());
  }
  //
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetChipIndex(Int_t lay,Int_t sta,Int_t chipInStave) const
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
Int_t AliITSUGeomTGeo::GetChipIndex(Int_t lay,Int_t sta, Int_t substa, Int_t chipInSStave) const
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
Int_t AliITSUGeomTGeo::GetChipIndex(Int_t lay,Int_t sta, Int_t substa, Int_t md, Int_t chipInMod) const
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
Bool_t AliITSUGeomTGeo::GetLayer(Int_t index,Int_t &lay,Int_t &indexInLr)  const
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
Int_t AliITSUGeomTGeo::GetLayer(Int_t index) const
{
  // Get chip layer, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  return lay;
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetStave(Int_t index) const
{
  // Get chip stave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index/fNChipsPerStave[lay];
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetHalfStave(Int_t index) const
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
Int_t AliITSUGeomTGeo::GetModule(Int_t index) const
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
Int_t AliITSUGeomTGeo::GetChipIdInLayer(Int_t index) const
{
  // Get chip number within layer, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index;
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetChipIdInStave(Int_t index) const
{
  // Get chip number within stave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index%fNChipsPerStave[lay];
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetChipIdInHalfStave(Int_t index) const
{
  // Get chip number within stave, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index%fNChipsPerHalfStave[lay];
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetChipIdInModule(Int_t index) const
{
  // Get chip number within module, from 0
  //
  int lay = 0;
  while(index>fLastChipIndex[lay]) lay++;
  index -= GetFirstChipIndex(lay);
  return index%fNChipsPerModule[lay];
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetChipId(Int_t index,Int_t &lay,Int_t &sta,Int_t &hsta, Int_t &mod, Int_t &chip)  const
{
  //
  // This routine computes the layer, stave, substave, module and chip number 
  // given the chip index number. 
  // Inputs:
  //     Int_t index  The chip index number, starting from zero.
  // Outputs:
  //     Int_t lay    The layer number. Starting from 0
  //     Int_t sta    The stave number. Starting from 0
  //     Int_t ssta   The halfstave number. Starting from 0
  //     Int_t mod    The module number. Starting from 0
  //     Int_t chip   The detector number. Starting from 0
  //
  lay  = GetLayer(index);
  index -= GetFirstChipIndex(lay);
  sta  = index/fNChipsPerStave[lay];
  index %= fNChipsPerStave[lay];
  hsta = fNHalfStaves[lay]>0 ? index/fNChipsPerHalfStave[lay] : -1;
  index %= fNChipsPerHalfStave[lay];
  mod  = fNModules[lay]>0 ? index/fNChipsPerModule[lay] : -1;
  chip = index%fNChipsPerModule[lay];
  //
  return kTRUE;
}

//______________________________________________________________________
TGeoHMatrix* AliITSUGeomTGeo::ExtractMatrixSens(Int_t index) const
{
  // Get the transformation matrix of the SENSOR (not necessary the same as the chip) 
  // for a given chip 'index' by quering the TGeoManager
  Int_t lay,stav,sstav,mod,chipInMod;
  GetChipId(index,lay,stav,sstav,mod,chipInMod);
  int wrID = fLr2Wrapper[lay];
  TString path = Form("/ALIC_1/%s_2/",GetITSVolPattern());
  if (wrID>=0) path += Form("%s%d_1/",GetITSWrapVolPattern(),wrID);
  path += Form("%s%d_1/%s%d_%d/",GetITSLayerPattern(),lay,GetITSStavePattern(),lay,stav);
  if (fNHalfStaves[lay]>0) path += Form("%s%d_%d/",GetITSHalfStavePattern(),lay,sstav);
  if (fNModules[lay]>0)   path += Form("%s%d_%d/",GetITSModulePattern(),lay,mod);
  path += Form("%s%d_%d/%s%d_1",GetITSChipPattern(),lay,chipInMod,GetITSSensorPattern(),lay);
  static TGeoHMatrix matTmp;
  gGeoManager->PushPath();
  if (!gGeoManager->cd(path.Data())) {
    gGeoManager->PopPath();
    AliError(Form("Error in cd-ing to %s",path.Data()));
    return 0;
  } // end if !gGeoManager
  matTmp = *gGeoManager->GetCurrentMatrix(); // matrix may change after cd
  //RSS
  //  printf("%d/%d/%d %s\n",lay,stav,detInSta,path.Data());
  //  mat->Print();
  // Retstore the modeler state.
  gGeoManager->PopPath();
  return &matTmp;
}

//______________________________________________________________________
void AliITSUGeomTGeo::Print(Option_t *) const
{
  // print
  printf("Geometry version %d, NLayers:%d NChips:%d\n",fVersion,fNLayers,fNChips);
  if (fVersion==kITSVNA) return;
  for (int i=0;i<fNLayers;i++) {
    printf("Lr%2d\tNStav:%2d\tNChips:%2d (%dx%-2d)\tNMod:%d\tNSubSt:%d\tNSt:%3d\tChipType:%3d\tChip#:%5d:%-5d\tWrapVol:%d\n",
	   i,fNStaves[i],fNChipsPerModule[i],fNChipRowsPerModule[i],
	   fNChipRowsPerModule[i] ? fNChipsPerModule[i]/fNChipRowsPerModule[i] : 0,
	   fNModules[i],fNHalfStaves[i],fNStaves[i],
	   fLrChipType[i],GetFirstChipIndex(i),GetLastChipIndex(i),fLr2Wrapper[i]);
  }
}

//______________________________________________________________________
const char* AliITSUGeomTGeo::GetSymName(Int_t index)  const
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
Bool_t AliITSUGeomTGeo::GetOrigMatrix(Int_t index, TGeoHMatrix &m) const
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
TGeoPNEntry* AliITSUGeomTGeo::GetPNEntry(Int_t index) const
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
const char* AliITSUGeomTGeo::ComposeSymNameITS()
{
  // sym name of the layer
  return "ITS";
}

//______________________________________________________________________
const char* AliITSUGeomTGeo::ComposeSymNameLayer(Int_t lr)
{
  // sym name of the layer
  return Form("%s/%s%d",ComposeSymNameITS(),GetITSLayerPattern(),lr);
}

//______________________________________________________________________
const char* AliITSUGeomTGeo::ComposeSymNameStave(Int_t lr, Int_t stave)
{
  // sym name of the stave at given layer
  return Form("%s/%s%d",ComposeSymNameLayer(lr),GetITSStavePattern(),stave);
}

//______________________________________________________________________
const char* AliITSUGeomTGeo::ComposeSymNameHalfStave(Int_t lr, Int_t stave, Int_t substave)
{
  // sym name of the stave at given layer
  return substave>=0 ? 
    Form("%s/%s%d",ComposeSymNameStave(lr,stave),GetITSHalfStavePattern(),substave) :
    ComposeSymNameStave(lr,stave);
}

//______________________________________________________________________
const char* AliITSUGeomTGeo::ComposeSymNameModule(Int_t lr, Int_t stave, Int_t substave, Int_t mod)
{
  // sym name of the substave at given layer/stave
  return mod>=0 ? 
    Form("%s/%s%d",ComposeSymNameHalfStave(lr,stave,substave),GetITSModulePattern(),mod) :
    ComposeSymNameHalfStave(lr,stave,substave);    
}

//______________________________________________________________________
const char* AliITSUGeomTGeo::ComposeSymNameChip(Int_t lr, Int_t sta, Int_t substave, Int_t mod, Int_t chip)
{
  // sym name of the chip in the given layer/stave/substave/module
  return Form("%s/%s%d",ComposeSymNameModule(lr,sta,substave,mod),GetITSChipPattern(),chip);
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::ExtractNumberOfLayers()
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
      if ( (lrID=ExtractVolumeCopy(name,GetITSLayerPattern()))<0 ) {
	AliFatal(Form("Failed to extract layer ID from the %s",name));
	exit(1);
      }
      //
      fLr2Wrapper[lrID] = -1; // not wrapped
    }
    else if (strstr(name,GetITSWrapVolPattern())) { // this is a wrapper volume, may cointain layers
      int wrID = -1;
      if ( (wrID=ExtractVolumeCopy(name,GetITSWrapVolPattern()))<0 ) {
	AliFatal(Form("Failed to extract wrapper ID from the %s",name));
	exit(1);
      }
      //
      TObjArray* nodesW = nd->GetNodes();
      int nNodesW = nodesW->GetEntriesFast();
      for (Int_t jw=0; jw<nNodesW; jw++) {
	TGeoNode* ndW = (TGeoNode*)nodesW->At(jw);
	if (strstr(ndW->GetName(),GetITSLayerPattern())) {
	  if ( (lrID=ExtractVolumeCopy(ndW->GetName(),GetITSLayerPattern()))<0 ) {
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
Int_t AliITSUGeomTGeo::ExtractNumberOfStaves(Int_t lay) const
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
Int_t AliITSUGeomTGeo::ExtractNumberOfHalfStaves(Int_t lay) const
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
Int_t AliITSUGeomTGeo::ExtractNumberOfModules(Int_t lay) const
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
Int_t AliITSUGeomTGeo::ExtractNChipsPerModule(Int_t lay, int &nrow)  const
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
Int_t AliITSUGeomTGeo::ExtractLayerChipType(Int_t lay)  const
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
Int_t AliITSUGeomTGeo::ExtractVolumeCopy(const char* name, const char* prefix) const
{
  // extract Number following the prefix in the name string
  TString nms = name;
  if (!nms.BeginsWith(prefix)) return -1;
  nms.Remove(0,strlen(prefix));
  if (!isdigit(nms.Data()[0])) return -1;
  return nms.Atoi();
  //
}

