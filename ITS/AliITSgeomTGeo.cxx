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
//    AliITSgeomTGeo is a simple interface class to TGeoManager          //
//    It is used in the simulation and reconstruction in order to        //
//    query the TGeo ITS geometry                                        //
//                                                                       //
//    author - cvetan.cheshkov@cern.ch                                   //
//    15/02/2007                                                         //
///////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>

#include "AliITSgeomTGeo.h"
#include "AliLog.h"
#include "AliAlignObj.h"

ClassImp(AliITSgeomTGeo)

const Int_t AliITSgeomTGeo::fgkNModules = 2198;
const Int_t AliITSgeomTGeo::fgkNLadders[kNLayers] = {20,40,14,22,34,38};
const Int_t AliITSgeomTGeo::fgkNDetectors[kNLayers] = {4,4,6,8,22,25};

//______________________________________________________________________
Int_t AliITSgeomTGeo::GetModuleIndex(Int_t lay,Int_t lad,Int_t det)
{
  // The method is taken from the old AliITSgeom class by Bjorn Nilsen
  //
  // This routine computes the module index number from the layer,
  // ladder, and detector numbers. The number of ladders and detectors
  // per layer is set statically
  // see above for details.
  // Inputs:
  //    Int_t lay  The layer number. Starting from 1.
  //    Int_t lad  The ladder number. Starting from 1.
  //    Int_t det  The detector number. Starting from 1.
  // Return:
  //    the module index number, starting from zero.
  //    -1 in case of error

  if (lay < 1 || lay > kNLayers) {
    AliErrorClass(Form("Invalid layer: %d (1 -> %d",lay,kNLayers));
    return -1;
  }

  if (lad < 1 || lad > fgkNLadders[lay-1] ||
      det < 1 || det > fgkNDetectors[lay-1]) {
    AliErrorClass(Form("Invalid layer,ladder,detector combination: %d, %d, %d",lay,lad,det));
    return -1;
  }

  Int_t index = fgkNDetectors[lay-1] * (lad-1) + (det-1);
  for(Int_t iLayer=0;iLayer < (lay-1); iLayer++)
    index += fgkNDetectors[iLayer]*fgkNLadders[iLayer];

  return index;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeo::GetLayer(Int_t index,Int_t &lay,Int_t &index2) 
{
  // The method is taken from the old AliITSgeom class by Bjorn Nilsen
  //
  // This routine computes the layer number for a
  // given the module index. The number of ladders and detectors
  // per layer is defined statically,
  // see above for details.
  // Inputs:
  //     Int_t index  The module index number, starting from zero.
  // Outputs:
  //     Int_t index2 The module index inside a layer, starting from zero.
  //     Int_t lay    The layer number. Starting from 1.
  // Return:
  //     kTRUE in case of valid index
  //     kFALSE in case of error

  if (index < 0 || index >= fgkNModules) {
    index2 = -1;
    AliErrorClass(Form("Invalid module index: %d (0 -> %d)",index,fgkNModules));
    return -1;
  }

  lay = 0;
  index2 = 0;
  do {
    index2 += fgkNLadders[lay]*fgkNDetectors[lay];
    lay++;
  } while(index2 <= index);
  index2 -= fgkNLadders[lay-1]*fgkNDetectors[lay-1];
  index2 = index - index2;

  return lay;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeo::GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det) 
{
  // The method is taken from the old AliITSgeom class by Bjorn Nilsen
  //
  // This routine computes the layer, ladder and detector number 
  // given the module index number. The number of ladders and detectors
  // per layer is defined statically,
  // see above for details.
  // Inputs:
  //     Int_t index  The module index number, starting from zero.
  // Outputs:
  //     Int_t lay    The layer number. Starting from 1.
  //     Int_t lad    The ladder number. Starting from 1.
  //     Int_t det    The detector number. Starting from 1.
  // Return:
  //     kTRUE in case of valid index
  //     kFALSE in case of error

  if (index < 0 || index >= fgkNModules) {
    lay = lad = det = -1;
    AliErrorClass(Form("Invalid module index: %d (0 -> %d)",index,fgkNModules));
    return kFALSE;
  }

  lay  = lad = det = 0;
  Int_t index2 = 0;
  do {
    index2 += fgkNLadders[lay]*fgkNDetectors[lay];
    lay++;
  } while(index2 <= index);
  index2 -= fgkNLadders[lay-1]*fgkNDetectors[lay-1];

  do {
    index2 += fgkNDetectors[lay-1];
    lad++;
  } while(index2 <= index);
  index2 -= fgkNDetectors[lay-1];

  det = index-index2+1;

  return kTRUE;
}

//______________________________________________________________________
const char* AliITSgeomTGeo::GetSymName(Int_t index) 
{
  // Get the TGeoPNEntry symbolic name
  // for a given module identified by 'index'

  if (index < 0 || index >= fgkNModules) {
    AliErrorClass(Form("Invalid ITS module index: %d (0 -> %d) !",index,fgkNModules));
    return NULL;
  }

  Int_t lay, index2;
  if (!GetLayer(index,lay,index2)) return NULL;

  return AliGeomManager::SymName((AliGeomManager::ELayerID)((lay-1)+AliGeomManager::kSPD1),index2);
}

//______________________________________________________________________
TGeoHMatrix* AliITSgeomTGeo::GetMatrix(Int_t index) 
{
  // Get the transformation matrix for a given module 'index'
  // by quering the TGeoManager

  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

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
Bool_t AliITSgeomTGeo::GetTranslation(Int_t index, Double_t t[3]) 
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
Bool_t AliITSgeomTGeo::GetRotation(Int_t index, Double_t r[9]) 
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
Bool_t AliITSgeomTGeo::GetOrigMatrix(Int_t index, TGeoHMatrix &m)
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
Bool_t AliITSgeomTGeo::GetOrigTranslation(Int_t index, Double_t t[3]) 
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
Bool_t AliITSgeomTGeo::GetOrigRotation(Int_t index, Double_t r[9]) 
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
const TGeoHMatrix* AliITSgeomTGeo::GetTracking2LocalMatrix(Int_t index)
{
  // Get the matrix which transforms from the tracking to local r.s.
  // The method queries directly the TGeoPNEntry

  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  const TGeoHMatrix *m = pne->GetMatrix();
  if (!m)
    AliErrorClass(Form("TGeoPNEntry (%s) contains no matrix !",pne->GetName()));

  return m;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeo::GetTrackingMatrix(Int_t index, TGeoHMatrix &m)
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
TGeoPNEntry* AliITSgeomTGeo::GetPNEntry(Int_t index)
{
  // Get a pointer to the TGeoPNEntry of a module
  // identified by 'index'
  // Returns NULL in case of invalid index,
  // missing TGeoManager or invalid symbolic name

  if (index < 0 || index >= fgkNModules) {
    AliErrorClass(Form("Invalid ITS module index: %d (0 -> %d) !",index,fgkNModules));
    return NULL;
  }
  
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliErrorClass("Can't get the matrix! gGeoManager doesn't exist or it is still opened!");
    return NULL;
  }

  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(GetSymName(index));
  if (!pne)
    AliErrorClass(Form("The symbolic volume name %s does not correspond to a physical entry!",
		       GetSymName(index)));

  return pne;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeo::LocalToGlobal(Int_t index,
				     const Double_t *loc, Double_t *glob)
{
  // Make the conversion from the local sensitive reference system to the global
  // reference system, for an arbitrary local position. The input is the pointer
  // to the array of local coordinates, the result is sent to the glob pointer.
  //
  // Please don't use this method to get the global coordinates of clusters, use
  // the direct method of AliCluster instead.

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
Bool_t AliITSgeomTGeo::GlobalToLocal(Int_t index,
				     const Double_t *glob, Double_t *loc)
{
  // Make the conversion from the global reference system to the sensitive local
  // reference system, for an arbitrary global position. The input is the pointer
  // to the array of global coordinates, the result is sent to the loc pointer.

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
Bool_t AliITSgeomTGeo::LocalToGlobalVect(Int_t index,
					 const Double_t *loc, Double_t *glob)
{
  // Make the conversion from the local sensitive reference system to the global
  // reference system, for an arbitrary vector. The input is the pointer to the
  // array of local coordinates, the result is sent to the glob pointer.

  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;
  ml->LocalToMasterVect(loc,glob);
  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeo::GlobalToLocalVect(Int_t index,
					 const Double_t *glob, Double_t *loc)
{
  // Make the conversion from the global reference system to the sensitive local
  // reference system, for an arbitrary vector. The input is the pointer to the
  // array of global coordinates, the result is sent to the loc pointer.

  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;
  ml->MasterToLocalVect(glob,loc);

  return kTRUE;
}
