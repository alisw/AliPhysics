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
//    AliITSgeomTGeoUpg is a simple interface class to TGeoManager       //
//    It is used in the simulation and reconstruction in order to        //
//    query the TGeo ITS geometry                                        //
//                                                                       //
//    author - cvetan.cheshkov@cern.ch                                   //
//    15/02/2007                                                         //
//    adapted to ITSupg 18/07/2012 - ruben.shahoyan@cern.ch              //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TDatime.h>

#include "AliITSgeomTGeoUpg.h"
#include "AliLog.h"
#include "AliAlignObj.h"

ClassImp(AliITSgeomTGeoUpg)


const char* AliITSgeomTGeoUpg::fgkITSVolName = "ITSV";
const char* AliITSgeomTGeoUpg::fgkITSLrName  = "ITSupgLayer";
const char* AliITSgeomTGeoUpg::fgkITSLadName = "ITSupgLadder";
const char* AliITSgeomTGeoUpg::fgkITSModName = "ITSupgModule";
const char* AliITSgeomTGeoUpg::fgkITSSensName ="ITSupgSensor";

const Int_t AliITSgeomTGeoUpg::fgkNModulesOld = 2198;
const Int_t AliITSgeomTGeoUpg::fgkNLaddersOld[AliITSgeomTGeoUpg::kNLayersOld] = {20,40,14,22,34,38};
const Int_t AliITSgeomTGeoUpg::fgkNDetectorsOld[AliITSgeomTGeoUpg::kNLayersOld] = {4,4,6,8,22,25};

Int_t  AliITSgeomTGeoUpg::fgVersion = 0;
Int_t  AliITSgeomTGeoUpg::fgNLayers = 0;
Int_t  AliITSgeomTGeoUpg::fgNModules = 0;
Int_t* AliITSgeomTGeoUpg::fgNLadders = 0;
Int_t* AliITSgeomTGeoUpg::fgNDetectors = 0;
Int_t* AliITSgeomTGeoUpg::fgLrDetType = 0;


//______________________________________________________________________
Int_t AliITSgeomTGeoUpg::GetModuleIndex(Int_t lay,Int_t lad,Int_t det)
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
  CheckInit();
  //
  if (lay < 1 || lay > fgNLayers) {
    AliErrorClass(Form("Invalid layer: %d (1 -> %d",lay,fgNLayers));
    return -1;
  }

  if (lad < 1 || lad > fgNLadders[lay-1] ||
      det < 1 || det > fgNDetectors[lay-1]) {
    AliErrorClass(Form("Invalid layer,ladder,detector combination: %d, %d, %d",lay,lad,det));
    return -1;
  }

  Int_t index = fgNDetectors[lay-1] * (lad-1) + (det-1);
  for(Int_t iLayer=0;iLayer < (lay-1); iLayer++)
    index += fgNDetectors[iLayer]*fgNLadders[iLayer];

  return index;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GetLayer(Int_t index,Int_t &lay,Int_t &index2) 
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
  CheckInit();
  //
  if (index < 0 || index >= fgNModules) {
    index2 = -1;
    AliErrorClass(Form("Invalid module index: %d (0 -> %d)",index,fgNModules));
    return -1;
  }

  lay = 0;
  index2 = 0;
  do {
    index2 += fgNLadders[lay]*fgNDetectors[lay];
    lay++;
  } while(index2 <= index);
  index2 -= fgNLadders[lay-1]*fgNDetectors[lay-1];
  index2 = index - index2;

  return lay;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det) 
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
  CheckInit();
  //
  if (index < 0 || index >= fgNModules) {
    lay = lad = det = -1;
    AliErrorClass(Form("Invalid module index: %d (0 -> %d)",index,fgNModules));
    return kFALSE;
  }

  lay  = lad = det = 0;
  Int_t index2 = 0;
  do {
    index2 += fgNLadders[lay]*fgNDetectors[lay];
    lay++;
  } while(index2 <= index);
  index2 -= fgNLadders[lay-1]*fgNDetectors[lay-1];

  do {
    index2 += fgNDetectors[lay-1];
    lad++;
  } while(index2 <= index);
  index2 -= fgNDetectors[lay-1];

  det = index-index2+1;

  return kTRUE;
}

//______________________________________________________________________
const char* AliITSgeomTGeoUpg::GetSymName(Int_t index) 
{
  // Get the TGeoPNEntry symbolic name
  // for a given module identified by 'index'
  CheckInit();
  //
  if (index < 0 || index >= fgNModules) {
    AliErrorClass(Form("Invalid ITS module index: %d (0 -> %d) !",index,fgNModules));
    return NULL;
  }

  Int_t lay, index2;
  if (!GetLayer(index,lay,index2)) return NULL;
  if (fgVersion == kITSVOld) return AliGeomManager::SymName((AliGeomManager::ELayerID)((lay-1)+AliGeomManager::kSPD1),index2);
  else {
    // RS: this is not optimal, but we cannod access directly AliGeomManager, since the latter has hardwired layers 
    TGeoPNEntry* pne = gGeoManager->GetAlignableEntryByUID( AliGeomManager::LayerToVolUID(lay,index2) );
    if (!pne) {
      AliErrorClass(Form("Failed to find alignable entry with index %d: (Lr%d Mod:%d) !",index,lay,index2));
      return NULL;
    }
    return pne->GetName();
  }
}

//______________________________________________________________________
TGeoHMatrix* AliITSgeomTGeoUpg::GetMatrix(Int_t index) 
{
  // Get the transformation matrix for a given module 'index'
  // by quering the TGeoManager
  CheckInit();
  //
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
Bool_t AliITSgeomTGeoUpg::GetTranslation(Int_t index, Double_t t[3]) 
{
  // Get the translation vector for a given module 'index'
  // by quering the TGeoManager
  CheckInit();
  //
  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *trans = m->GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GetRotation(Int_t index, Double_t r[9]) 
{
  // Get the rotation matrix for a given module 'index'
  // by quering the TGeoManager
  CheckInit();
  //
  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *rot = m->GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GetOrigMatrix(Int_t index, TGeoHMatrix &m)
{
  // Get the original (ideal geometry) TGeo matrix for
  // a given module identified by 'index'.
  // The method is slow, so it should be used
  // with great care.
  CheckInit();
  //
  m.Clear();

  const char *symname = GetSymName(index);
  if (!symname) return kFALSE;

  return AliGeomManager::GetOrigGlobalMatrix(symname,m);
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GetOrigTranslation(Int_t index, Double_t t[3]) 
{
  // Get the original translation vector (ideal geometry)
  // for a given module 'index' by quering the TGeoManager
  CheckInit();
  //
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *trans = m.GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GetOrigRotation(Int_t index, Double_t r[9]) 
{
  // Get the original rotation matrix (ideal geometry)
  // for a given module 'index' by quering the TGeoManager
  CheckInit();
  //
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *rot = m.GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
const TGeoHMatrix* AliITSgeomTGeoUpg::GetTracking2LocalMatrix(Int_t index)
{
  // Get the matrix which transforms from the tracking to local r.s.
  // The method queries directly the TGeoPNEntry
  CheckInit();
  //
  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  const TGeoHMatrix *m = pne->GetMatrix();
  if (!m)
    AliErrorClass(Form("TGeoPNEntry (%s) contains no matrix !",pne->GetName()));

  return m;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GetTrackingMatrix(Int_t index, TGeoHMatrix &m)
{
  // Get the matrix which transforms from the tracking r.s. to
  // the global one.
  // Returns kFALSE in case of error.
  CheckInit();
  //
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
TGeoPNEntry* AliITSgeomTGeoUpg::GetPNEntry(Int_t index)
{
  // Get a pointer to the TGeoPNEntry of a module
  // identified by 'index'
  // Returns NULL in case of invalid index,
  // missing TGeoManager or invalid symbolic name
  CheckInit();
  //
  if (index < 0 || index >= fgNModules) {
    AliErrorClass(Form("Invalid ITS module index: %d (0 -> %d) !",index,fgNModules));
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
Bool_t AliITSgeomTGeoUpg::LocalToGlobal(Int_t index,
				     const Double_t *loc, Double_t *glob)
{
  // Make the conversion from the local sensitive reference system to the global
  // reference system, for an arbitrary local position. The input is the pointer
  // to the array of local coordinates, the result is sent to the glob pointer.
  //
  // Please don't use this method to get the global coordinates of clusters, use
  // the direct method of AliCluster instead.
  CheckInit();
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
Bool_t AliITSgeomTGeoUpg::GlobalToLocal(Int_t index,
				     const Double_t *glob, Double_t *loc)
{
  // Make the conversion from the global reference system to the sensitive local
  // reference system, for an arbitrary global position. The input is the pointer
  // to the array of global coordinates, the result is sent to the loc pointer.
  CheckInit();
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
Bool_t AliITSgeomTGeoUpg::LocalToGlobalVect(Int_t index,
					 const Double_t *loc, Double_t *glob)
{
  // Make the conversion from the local sensitive reference system to the global
  // reference system, for an arbitrary vector. The input is the pointer to the
  // array of local coordinates, the result is sent to the glob pointer.
  CheckInit();
  //
  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;
  ml->LocalToMasterVect(loc,glob);
  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::GlobalToLocalVect(Int_t index,
					 const Double_t *glob, Double_t *loc)
{
  // Make the conversion from the global reference system to the sensitive local
  // reference system, for an arbitrary vector. The input is the pointer to the
  // array of global coordinates, the result is sent to the loc pointer.
  CheckInit();
  //
  TGeoHMatrix *ml = GetMatrix(index);
  if (!ml) return kFALSE;
  ml->MasterToLocalVect(glob,loc);

  return kTRUE;
}

//______________________________________________________________________
void AliITSgeomTGeoUpg::BuildITS()
{
  // build ITS info from TGeo 
  if (fgVersion!=kITSVNA) return; // already initialized
  // 
  if (!gGeoManager) AliFatalClass("Geometry is not loaded");
  //
  // get ITS version
  TGeoVolume *itsV = gGeoManager->GetVolume(fgkITSVolName);
  if (!gGeoManager) AliFatalClass(Form("ITS volume %s is not in the geometry",fgkITSVolName));
  //
  Int_t minor = 0;
  TDatime datetime;
  const Char_t *title = itsV->GetTitle();
  if(!ReadVersionString(title,fgVersion,minor,datetime))
    AliFatalClass(Form("Can't read title %s to extract version",title));
  //
  if (fgVersion==kITSVNA) AliFatalClass(Form("Failed to detrmine ITS version from title %s",title));
  //
  (fgVersion==kITSVOld) ? BuildITSOld() : BuildITSUpg();
  //
}

//______________________________________________________________________
void AliITSgeomTGeoUpg::BuildITSOld()
{
  // assign old ITS
  fgNLayers  = kNLayersOld;
  fgNModules = fgkNModulesOld;
  fgNLadders = (Int_t*)fgkNLaddersOld;
  fgNDetectors = (Int_t*)fgkNDetectorsOld;
  fgLrDetType = 0;
}

//______________________________________________________________________
void AliITSgeomTGeoUpg::BuildITSUpg()
{
  // exract upg ITS parameters from TGeo
  fgNLayers    = ExtractNumberOfLayers();
  if (!fgNLayers) return;
  //
  fgNLadders   = new Int_t[fgNLayers];
  fgNDetectors = new Int_t[fgNLayers];
  fgLrDetType  = new Int_t[fgNLayers];
  fgNModules = 0;
  for (int i=0;i<fgNLayers;i++) {
    fgNLadders[i]   = ExtractNumberOfLadders(i+1);
    fgNDetectors[i] = ExtractNumberOfDetectors(i+1);
    fgLrDetType[i]  = ExtractLayerDetType(i+1);
    fgNModules     += fgNLadders[i]*fgNDetectors[i];
  }
  //
}

//______________________________________________________________________
Bool_t AliITSgeomTGeoUpg::ReadVersionString(const Char_t *str,
					    Int_t &maj,Int_t &min,
					    TDatime &dt)
{
  // fills the string str with the major and minor version number
  // Inputs:
  //   Char_t *str   The character string to holding the major and minor
  //                 version numbers in
  // Outputs:
  //   Int_t           maj  The major number
  //   Int_t           min  The minor number
  //   TDatime         dt   The date and time of the cvs commit
  // Return:
  //   kTRUE if no errors
  enum {kv11=11,kv110=110,kvUpgrade=20}; // RS: to make consistent global numbering scheme
  
  Bool_t ok;
  Char_t cvsRevision[10],cvsDate[11],cvsTime[9];
  Int_t i,m,n=strlen(str),year,month,day,hours,minuts,seconds;
  memset(cvsRevision,0,10*sizeof(Char_t));
  memset(cvsDate,0,11*sizeof(Char_t));    
  memset(cvsTime,0,9*sizeof(Char_t));
  
  if(n<35) return kFALSE; // not enough space for numbers
  m = sscanf(str,"Major Version= %d  Minor Version= %d Revision: %9s "
	     "Date: %10s %8s",&i,&min,cvsRevision,cvsDate,cvsTime);
  ok = m==5;
  if(!ok) return !ok;
  m = sscanf(cvsDate,"%d/%d/%d",&year,&month,&day);
  ok = m==3;
  if(!ok) return !ok;
  m = sscanf(cvsTime,"%d:%d:%d",&hours,&minuts,&seconds);
  ok = m==3;
  if(!ok) return !ok;
  dt.Set(year,month,day,hours,minuts,seconds);
  //
  switch (i){
  case kv110:
  case kv11:{
    maj = kITSVOld;
  } break;
  case kvUpgrade:{
    maj = kITSVUpg;
  } break;
  default:{
    maj = kITSVNA;
  } break;
  } // end switch
  return ok;
}

//______________________________________________________________________
Int_t AliITSgeomTGeoUpg::ExtractNumberOfLayers()
{
  // Determines the number of layers in the Upgrade Geometry
  //
  // Inputs:
  //   none
  // Outputs:
  //   none
  // Return:
  //   the number of layers in the current geometry
  //   -1 if not Upgrade Geometry
  // MS
  //
  Int_t numberOfLayers = 0;
  //
  TGeoVolume *itsV = gGeoManager->GetVolume(fgkITSVolName);
  if (!itsV) AliFatalClass(Form("ITS volume %s is not in the geometry",fgkITSVolName));
  //
  // Loop on all ITSV nodes, count Layer volumes by checking names
  Int_t nNodes = itsV->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(itsV->GetNodes()->At(j)->GetName(),fgkITSLrName)) numberOfLayers++;
  //  
  return numberOfLayers;
}

//______________________________________________________________________
Int_t AliITSgeomTGeoUpg::ExtractNumberOfLadders(const Int_t lay)
{
  // Determines the number of layers in the Upgrade Geometry
  //
  // Inputs:
  //   lay: layer number, starting from 1
  // Outputs:
  //   none
  // Return:
  //   the number of ladders in layer lay
  //   -1 if not Upgrade Geometry
  // MS
  Int_t numberOfLadders = 0;
  char laynam[30];
  snprintf(laynam, 30, "%s%d",fgkITSLrName,lay);
  TGeoVolume* volLr = gGeoManager->GetVolume(laynam);
  if (!volLr) AliFatalClass(Form("can't find %s volume",laynam));
  //
  // Loop on all layer nodes, count Ladder volumes by checking names
  Int_t nNodes = volLr->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(volLr->GetNodes()->At(j)->GetName(),fgkITSLadName)) numberOfLadders++;
  //
  return numberOfLadders;
  //
}

//______________________________________________________________________
Int_t AliITSgeomTGeoUpg::ExtractNumberOfDetectors(const Int_t lay) 
{
  // Determines the number of detectors per ladder in the Upgrade Geometry
  //
  // Inputs:
  //   lay: layer number from 1
  // Outputs:
  //   none
  // Return:
  //   the number of modules per ladder in layer lay
  //   -1 if not Upgrade Geometry
  // MS
  Int_t numberOfModules = 0;
  char laddnam[30];
  snprintf(laddnam, 30, "%s%d", fgkITSLadName,lay);
  TGeoVolume* volLd = gGeoManager->GetVolume(laddnam);
  if (!volLd) AliFatalClass(Form("can't find %s volume",laddnam));
  //
  // Loop on all ladder nodes, count Module volumes by checking names
  Int_t nNodes = volLd->GetNodes()->GetEntries();
  for (Int_t j=0; j<nNodes; j++) if (strstr(volLd->GetNodes()->At(j)->GetName(),fgkITSModName)) numberOfModules++;
  //
  return numberOfModules;
  //
}

//______________________________________________________________________
Int_t AliITSgeomTGeoUpg::ExtractLayerDetType(const Int_t lay) 
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
  if (!volLd) {AliFatalClass(Form("can't find %s volume",laddnam)); return -1;}
  //
  return volLd->GetUniqueID();
  //
}
