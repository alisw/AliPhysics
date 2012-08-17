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
#include <TDatime.h>

#include "AliITSUGeomTGeo.h"
#include "AliLog.h"
#include "AliAlignObj.h"

ClassImp(AliITSUGeomTGeo)


const char* AliITSUGeomTGeo::fgkITSVolName = "ITSV";
const char* AliITSUGeomTGeo::fgkITSLrName  = "ITSULayer";
const char* AliITSUGeomTGeo::fgkITSLadName = "ITSULadder";
const char* AliITSUGeomTGeo::fgkITSModName = "ITSUModule";
const char* AliITSUGeomTGeo::fgkITSSensName ="ITSUSensor";
const char* AliITSUGeomTGeo::fgkITSDetTypeName[AliITSUGeomTGeo::kNDetTypes] = {"Pix"};
//
TString     AliITSUGeomTGeo::fgITSsegmFileName = "itsSegmentations.root";

//______________________________________________________________________
AliITSUGeomTGeo::AliITSUGeomTGeo(Bool_t build)
:  fVersion(kITSVNA)
  ,fNLayers(0)
  ,fNModules(0)
  ,fNLadders(0)
  ,fLrDetType(0)
  ,fNDetectors(0)
  ,fLastModIndex(0)
{
  // default c-tor
  if (build) BuildITS();
}

//______________________________________________________________________
AliITSUGeomTGeo::AliITSUGeomTGeo(const AliITSUGeomTGeo &src)
  :TObject(src)
  ,fVersion(src.fVersion)
  ,fNLayers(src.fNLayers)
  ,fNModules(src.fNModules)
  ,fNLadders(0)
  ,fLrDetType(0)
  ,fNDetectors(0)
  ,fLastModIndex(0)
{
  // copy c-tor
  if (fNLayers) {
    fNLadders   = new Int_t[fNLayers];
    fNDetectors = new Int_t[fNLayers];
    fLrDetType  = new Int_t[fNLayers];
    fLastModIndex   = new Int_t[fNLayers];
    for (int i=fNLayers;i--;) {
      fNLadders[i] = src.fNLadders[i];
      fNDetectors[i] = src.fNDetectors[i];
      fLrDetType[i]  = src.fLrDetType[i];
      fLastModIndex[i] = src.fLastModIndex[i];
    }
  }
}

//______________________________________________________________________
AliITSUGeomTGeo::~AliITSUGeomTGeo()
{
  //d-tor
  delete[] fNLadders;
  delete[] fLrDetType;
  delete[] fNDetectors;
  delete[] fLastModIndex;
}


//______________________________________________________________________
AliITSUGeomTGeo& AliITSUGeomTGeo::operator=(const AliITSUGeomTGeo &src)
{
  // cp op.
  if (this!=&src) {
    delete[] fNLadders;
    delete[] fLrDetType;
    delete[] fNDetectors;
    delete[] fLastModIndex;
    fNLadders = fLrDetType = fNDetectors = fLastModIndex = 0;
    fVersion = src.fVersion;
    fNLayers = src.fNLayers;
    fNModules = src.fNModules;
    if (fNLayers) {
      fNLadders   = new Int_t[fNLayers];
      fNDetectors = new Int_t[fNLayers];
      fLrDetType  = new Int_t[fNLayers];
      fLastModIndex   = new Int_t[fNLayers];
      for (int i=fNLayers;i--;) {
	fNLadders[i] = src.fNLadders[i];
	fNDetectors[i] = src.fNDetectors[i];
	fLrDetType[i]  = src.fLrDetType[i];
	fLastModIndex[i] = src.fLastModIndex[i];
      }
    }    
  }
  return *this;
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetModuleIndex(Int_t lay,Int_t lad,Int_t det) const
{
  // This routine computes the module index number from the layer,
  // ladder, and detector numbers. The number of ladders and detectors
  // per layer is set statically
  // see above for details.
  // Inputs:
  //    Int_t lay  The layer number. Starting from 0.
  //    Int_t lad  The ladder number. Starting from 0
  //    Int_t det  The detector number in the ladder. Starting from 0
  //
  return GetFirstModIndex(lay) + fNDetectors[lay]*lad + det;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetLayer(Int_t index,Int_t &lay,Int_t &index2)  const
{
  // This routine computes the layer number for a
  // given the module index. The 
  // Inputs:
  //     Int_t index  The module index number, starting from zero.
  // Outputs:
  //     Int_t index2 The module index inside a layer, starting from zero.
  //     Int_t lay    The layer number. Starting from 0.
  //
  lay = GetLayer(index);
  index2 = index - GetFirstModIndex(lay);
  return kTRUE;
  //
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetLayer(Int_t index) const
{
  // Get module layer, from 0
  //
  int lay = 0;
  while(index>fLastModIndex[lay]) lay++;
  return lay;
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetLadder(Int_t index) const
{
  // Get module ladder, from 0
  //
  int lay = 0;
  while(index>fLastModIndex[lay]) lay++;
  index -= GetFirstModIndex(lay);
  return index/fNDetectors[lay];
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetModIdInLayer(Int_t index) const
{
  // Get module number within layer, from 0
  //
  int lay = 0;
  while(index>fLastModIndex[lay]) lay++;
  index -= GetFirstModIndex(lay);
  return index;
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::GetModIdInLadder(Int_t index) const
{
  // Get module number within ladder, from 0
  //
  int lay = 0;
  while(index>fLastModIndex[lay]) lay++;
  index -= GetFirstModIndex(lay);
  return index%fNDetectors[lay];
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det)  const
{
  // The method is taken from the old AliITSgeom class by Bjorn Nilsen
  //
  // This routine computes the layer, ladder and detector number 
  // given the module index number. 
  // Inputs:
  //     Int_t index  The module index number, starting from zero.
  // Outputs:
  //     Int_t lay    The layer number. Starting from 0
  //     Int_t lad    The ladder number. Starting from 0
  //     Int_t det    The detector number. Starting from 0
  //
  lay  = GetLayer(index);
  index -= GetFirstModIndex(lay);
  lad  = index/fNDetectors[lay];
  det  = index%fNDetectors[lay];
  return kTRUE;
}

//______________________________________________________________________
const char* AliITSUGeomTGeo::GetSymName(Int_t index)  const
{
  // Get the TGeoPNEntry symbolic name
  // for a given module identified by 'index'
  //
  Int_t lay, index2;
  if (!GetLayer(index,lay,index2)) return NULL;
  // return AliGeomManager::SymName((AliGeomManager::ELayerID)((lay-1)+AliGeomManager::kSPD1),index2);
  // RS: this is not optimal, but we cannod access directly AliGeomManager, since the latter has hardwired layers 
  //  TGeoPNEntry* pne = gGeoManager->GetAlignableEntryByUID( AliGeomManager::LayerToVolUID(lay+1,index2) );
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntryByUID( ModuleVolUID(index) );
  if (!pne) {
    AliError(Form("Failed to find alignable entry with index %d: (Lr%d Mod:%d) !",index,lay,index2));
    return NULL;
  }
  return pne->GetName();
}

//______________________________________________________________________
TGeoHMatrix* AliITSUGeomTGeo::GetMatrix(Int_t index)  const
{
  // Get the transformation matrix for a given module 'index'
  // by quering the TGeoManager
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
  TGeoHMatrix *mat = gGeoManager->GetCurrentMatrix();
  gGeoManager->PopPath();
  return mat;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetTranslation(Int_t index, Double_t t[3])  const
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
Bool_t AliITSUGeomTGeo::GetRotation(Int_t index, Double_t r[9])  const
{
  // Get the rotation matrix for a given module 'index'
  // by quering the TGeoManager
  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *rot = m->GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetOrigMatrix(Int_t index, TGeoHMatrix &m) const
{
  // Get the original (ideal geometry) TGeo matrix for
  // a given module identified by 'index'.
  // The method is slow, so it should be used
  // with great care.
  m.Clear();

  const char *symname = GetSymName(index);
  if (!symname) return kFALSE;

  return AliGeomManager::GetOrigGlobalMatrix(symname,m);
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetOrigTranslation(Int_t index, Double_t t[3])  const
{
  // Get the original translation vector (ideal geometry)
  // for a given module 'index' by quering the TGeoManager
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *trans = m.GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetOrigRotation(Int_t index, Double_t r[9])  const
{
  // Get the original rotation matrix (ideal geometry)
  // for a given module 'index' by quering the TGeoManager
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *rot = m.GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
const TGeoHMatrix* AliITSUGeomTGeo::GetTracking2LocalMatrix(Int_t index) const
{
  // Get the matrix which transforms from the tracking to local r.s.
  // The method queries directly the TGeoPNEntry
  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  const TGeoHMatrix *m = pne->GetMatrix();
  if (!m)
    AliError(Form("TGeoPNEntry (%s) contains no matrix !",pne->GetName()));

  return m;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GetTrackingMatrix(Int_t index, TGeoHMatrix &m) const
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

//______________________________________________________________________
TGeoHMatrix* AliITSUGeomTGeo::GetMatrixSens(Int_t lay, Int_t ladd, Int_t detInLad)  const
{
  // Get the transformation matrix of the SENSOR (not ncessary the same as the module) for a given module 'index'
  // by quering the TGeoManager
  const TString kPathBase = Form("/ALIC_1/%s_2/",AliITSUGeomTGeo::GetITSVolPattern());
  const TString kNames = Form("%%s%s%%d_1/%s%%d_%%d/%s%%d_%%d/%s%%d_%%d"
			      ,AliITSUGeomTGeo::GetITSLayerPattern()
			      ,AliITSUGeomTGeo::GetITSLadderPattern()
			      ,AliITSUGeomTGeo::GetITSModulePattern()
			      ,AliITSUGeomTGeo::GetITSSensorPattern());
  TString path;
  //
  path.Form(kNames.Data(),kPathBase.Data(),lay,lay,ladd,lay,detInLad,lay,1);
  gGeoManager->PushPath();
  if (!gGeoManager->cd(path.Data())) {
    gGeoManager->PopPath();
    AliError(Form("Error in cd-ing to %s",path.Data()));
    return 0;
  } // end if !gGeoManager
  TGeoHMatrix* mat = gGeoManager->GetCurrentMatrix();
  // Retstore the modeler state.
  gGeoManager->PopPath();
  return mat;
}


//______________________________________________________________________
TGeoPNEntry* AliITSUGeomTGeo::GetPNEntry(Int_t index) const
{
  // Get a pointer to the TGeoPNEntry of a module
  // identified by 'index'
  // Returns NULL in case of invalid index,
  // missing TGeoManager or invalid symbolic name
  //
  if (index >= fNModules) {
    AliError(Form("Invalid ITS module index: %d (0 -> %d) !",index,fNModules));
    return NULL;
  }
  
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the matrix! gGeoManager doesn't exist or it is still opened!");
    return NULL;
  }

  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(GetSymName(index));
  if (!pne) AliError(Form("The symbolic volume name %s does not correspond to a physical entry!",GetSymName(index)));
  //
  return pne;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::LocalToGlobal(Int_t index,const Double_t *loc, Double_t *glob) const
{
  // Make the conversion from the local sensitive reference system to the global
  // reference system, for an arbitrary local position. The input is the pointer
  // to the array of local coordinates, the result is sent to the glob pointer.
  //
  // Please don't use this method to get the global coordinates of clusters, use
  // the direct method of AliCluster instead.
  //
  const TGeoHMatrix *m2 = GetTracking2LocalMatrix(index);
  if (!m2) return kFALSE;

  // The shift (in local y only) between alignable and sensitive volume
  // is extracted directly from the Tracking2Local matrix
  Double_t locSens[] = {loc[0], loc[1]+m2->GetTranslation()[1], loc[2]};

  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;
  ml->LocalToMaster(locSens,glob);
  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GlobalToLocal(Int_t index, const Double_t *glob, Double_t *loc) const
{
  // Make the conversion from the global reference system to the sensitive local
  // reference system, for an arbitrary global position. The input is the pointer
  // to the array of global coordinates, the result is sent to the loc pointer.
  //
  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;

  const TGeoHMatrix *m2 = GetTracking2LocalMatrix(index);
  if (!m2) return kFALSE;
  ml->MasterToLocal(glob,loc);
  // The shift (in local y only) between alignable and sensitive volume
  // is extracted directly from the Tracking2Local matrix
  loc[1] -= m2->GetTranslation()[1];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::LocalToGlobalVect(Int_t index, const Double_t *loc, Double_t *glob) const
{
  // Make the conversion from the local sensitive reference system to the global
  // reference system, for an arbitrary vector. The input is the pointer to the
  // array of local coordinates, the result is sent to the glob pointer.
  //
  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;
  ml->LocalToMasterVect(loc,glob);
  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSUGeomTGeo::GlobalToLocalVect(Int_t index, const Double_t *glob, Double_t *loc) const
{
  // Make the conversion from the global reference system to the sensitive local
  // reference system, for an arbitrary vector. The input is the pointer to the
  // array of global coordinates, the result is sent to the loc pointer.
  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;
  ml->MasterToLocalVect(glob,loc);
  return kTRUE;
}

//______________________________________________________________________
void AliITSUGeomTGeo::BuildITS()
{
  // exract upg ITS parameters from TGeo
  if (fVersion!=kITSVNA) {AliWarning("Already built"); return; // already initialized}
    if (!gGeoManager) AliFatal("Geometry is not loaded");
  }
  fNLayers    = ExtractNumberOfLayers();
  if (!fNLayers) return;
  //
  fNLadders   = new Int_t[fNLayers];
  fNDetectors = new Int_t[fNLayers];
  fLrDetType  = new Int_t[fNLayers];
  fLastModIndex   = new Int_t[fNLayers];
  fNModules = 0;
  for (int i=0;i<fNLayers;i++) {
    fNLadders[i]   = ExtractNumberOfLadders(i);
    fNDetectors[i] = ExtractNumberOfDetectors(i);
    fLrDetType[i]  = ExtractLayerDetType(i);
    fNModules     += fNLadders[i]*fNDetectors[i];
    fLastModIndex[i]   = fNModules-1;
  }
  //
  fVersion = kITSVUpg;
  //
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::ExtractNumberOfLayers() const
{
  // Determines the number of layers in the Upgrade Geometry
  //
  Int_t numberOfLayers = 0;
  //
  TGeoVolume *itsV = gGeoManager->GetVolume(fgkITSVolName);
  if (!itsV) AliFatal(Form("ITS volume %s is not in the geometry",fgkITSVolName));
  //
  // Loop on all ITSV nodes, count Layer volumes by checking names
  Int_t nNodes = itsV->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(itsV->GetNodes()->At(j)->GetName(),fgkITSLrName)) numberOfLayers++;
  //  
  return numberOfLayers;
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::ExtractNumberOfLadders(Int_t lay) const
{
  // Determines the number of layers in the Upgrade Geometry
  //
  // Inputs:
  //   lay: layer number, starting from 0
  //
  // MS
  Int_t numberOfLadders = 0;
  char laynam[30];
  snprintf(laynam, 30, "%s%d",fgkITSLrName,lay);
  TGeoVolume* volLr = gGeoManager->GetVolume(laynam);
  if (!volLr) AliFatal(Form("can't find %s volume",laynam));
  //
  // Loop on all layer nodes, count Ladder volumes by checking names
  Int_t nNodes = volLr->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(volLr->GetNodes()->At(j)->GetName(),fgkITSLadName)) numberOfLadders++;
  //
  return numberOfLadders;
  //
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::ExtractNumberOfDetectors(Int_t lay)  const
{
  // Determines the number of detectors per ladder in the Upgrade Geometry
  //
  // Inputs:
  //   lay: layer number from 0
  // MS
  Int_t numberOfModules = 0;
  char laddnam[30];
  snprintf(laddnam, 30, "%s%d", fgkITSLadName,lay);
  TGeoVolume* volLd = gGeoManager->GetVolume(laddnam);
  if (!volLd) AliFatal(Form("can't find %s volume",laddnam));
  //
  // Loop on all ladder nodes, count Module volumes by checking names
  Int_t nNodes = volLd->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(volLd->GetNodes()->At(j)->GetName(),fgkITSModName)) numberOfModules++;
  //
  return numberOfModules;
  //
}

//______________________________________________________________________
Int_t AliITSUGeomTGeo::ExtractLayerDetType(Int_t lay)  const
{
  // Determines the layer detector type the Upgrade Geometry
  //
  // Inputs:
  //   lay: layer number from 1
  // Outputs:
  //   none
  // Return:
  //   detector type id for the layer
  // MS
  char laddnam[30];
  snprintf(laddnam, 30, "%s%d", fgkITSLrName,lay);
  TGeoVolume* volLd = gGeoManager->GetVolume(laddnam);
  if (!volLd) {AliFatal(Form("can't find %s volume",laddnam)); return -1;}
  //
  return volLd->GetUniqueID();
  //
}

//______________________________________________________________________
UInt_t AliITSUGeomTGeo::ComposeDetTypeID(UInt_t segmId)
{
  if (segmId>=kMaxSegmPerDetType) AliFatalClass(Form("Id=%d is >= max.allowed %d",segmId,kMaxSegmPerDetType));
  return segmId + kDetTypePix*kMaxSegmPerDetType;
}

//______________________________________________________________________
void AliITSUGeomTGeo::Print(Option_t *) const
{
  // print
  printf("Geometry version %d, NLayers:%d NModules:%d\n",fVersion,fNLayers,fNModules);
  if (fVersion==kITSVNA) return;
  for (int i=0;i<fNLayers;i++) {
    printf("Lr%2d\tNLadd:%2d\tNDet:%2d\tDetType:%3d\tMod#:%4d:%4d\n",
	   i,fNLadders[i],fNDetectors[i],fLrDetType[i],GetFirstModIndex(i),GetLastModIndex(i));
  }
}
