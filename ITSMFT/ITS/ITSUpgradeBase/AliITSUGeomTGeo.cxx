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

//______________________________________________________________________
AliITSUGeomTGeo::AliITSUGeomTGeo(Bool_t build, Bool_t loadSegm)
  :AliITSMFTGeomTGeo(build,loadSegm)
  ,fVersion(kITSVNA)
{
  // default c-tor
  if (build) BuildITS(loadSegm);
}

//______________________________________________________________________
AliITSUGeomTGeo::AliITSUGeomTGeo(const AliITSUGeomTGeo &src)
  :AliITSMFTGeomTGeo(src)
  ,fVersion(src.fVersion)
{
  // copy c-tor
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
    fNChips              += fNChipsPerLayer[i];
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
  TString path = Form("/ALIC_1/%s_2/",AliITSMFTGeomTGeo::GetITSVolPattern());
  if (wrID>=0) path += Form("%s%d_1/",GetITSWrapVolPattern(),wrID);
  path += Form("%s%d_1/%s%d_%d/",AliITSMFTGeomTGeo::GetITSLayerPattern(),lay,AliITSMFTGeomTGeo::GetITSStavePattern(),lay,stav);
  if (fNHalfStaves[lay]>0) path += Form("%s%d_%d/",AliITSMFTGeomTGeo::GetITSHalfStavePattern(),lay,sstav);
  if (fNModules[lay]>0)   path += Form("%s%d_%d/",AliITSMFTGeomTGeo::GetITSModulePattern(),lay,mod);
  path += Form("%s%d_%d/%s%d_1",AliITSMFTGeomTGeo::GetITSChipPattern(),lay,chipInMod,AliITSMFTGeomTGeo::GetITSSensorPattern(),lay);
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
