// -*- C++ -*-

#include <TGeoManager.h>
#include <TGeoMatrix.h>

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliGeomManager.h"

#include "AliITSsegmentationSPD.h"
#include "AliITSAlignMille2Module.h"
#include "AliSPDUtils.h"

#include "AliLog.h"

#include "AliSPDGeometryUtils.h"

ClassImp(AliSPDGeometryUtils);

Bool_t AliSPDGeometryUtils::LoadGeom(Int_t rn) {
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man)
    return kFALSE;

  man->SetRun(rn);

  AliCDBEntry* entry = man->Get("GRP/Geometry/Data");
  if (!entry)
    return kFALSE;

  TGeoManager *geoManager = dynamic_cast<TGeoManager *>(entry->GetObject());
  if (!geoManager)
    return kFALSE;

  AliGeomManager::SetGeometry(geoManager);
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  return kTRUE;
}

Bool_t AliSPDGeometryUtils::GetChipXYZ(Int_t chipKey,     // 0-399 inner layer, 400-1199 outer layer
				       TVector3 pos[4]) { // position of the corners of the FO chip
  UInt_t module(999), offchip(999);
  if (!AliSPDUtils::GetOfflineFromOfflineChipKey(chipKey, module, offchip))
    return kFALSE;

  const UInt_t hs = AliSPDUtils::GetOnlineHSFromOffline(module);
  if (hs<2) offchip = 4 - offchip; // inversion  in the inner layer...

  const Int_t col[4] = {
    hs<2 ?  0 : 31,
    hs<2 ? 31 :  0,
    hs<2 ? 31 :  0,
    hs<2 ?  0 : 31
  };
  const Int_t aa[4] = { 0, 0, 255, 255 };

  const AliITSsegmentationSPD seg;
  for(Int_t ic=0; ic<4; ++ic) {
    Float_t localchip[3] = { 0., 0., 0. };
    seg.DetToLocal(aa[ic], col[ic] + 32*offchip, localchip[0], localchip[2]);

    const Double_t local[3] = { localchip[0], localchip[1], localchip[2] };
    Double_t glochip[3] = { 0., 0., 0. };
    if (!Loc2Glo(module, local, glochip)) {
      return kFALSE;
    }
    pos[ic].SetXYZ(glochip[0], glochip[1], glochip[2]);
  }
  return kTRUE;
}

Bool_t AliSPDGeometryUtils::Loc2Glo(Int_t id,
				    const Double_t *loc,
				    Double_t *glo) {
  static TGeoHMatrix mat;
  const Int_t vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(id);
  if (vid < 0) {
    AliErrorClassF("AliCDMesonUtils Did not find module with such ID %d\n", id);
    return kFALSE;
  }
  AliITSAlignMille2Module::SensVolMatrix(vid, &mat);
  mat.LocalToMaster(loc, glo);
  return kTRUE;
}
