// -*- C++ -*-

#include <TObject.h>

#include <TVector3.h>

// adapted from: AliPhysics/PWGUD/DIFFRACTIVE/xsAndTwoProng/AliCDMesonUtils.h
struct AliSPDGeometryUtils : public TObject {
  AliSPDGeometryUtils()
    : TObject() {}

  virtual ~AliSPDGeometryUtils() {}

  // to be called before GetChipXYZ
  static Bool_t LoadGeom(Int_t rn);

  static Bool_t GetChipXYZ(Int_t chipKey,    // 0-399 inner layer, 400-1199 outer layer
			   TVector3 pos[4]); // position of the corners of the FO chip

  static Bool_t Loc2Glo(Int_t id,
			const Double_t *loc,
			Double_t *glo);

  ClassDef(AliSPDGeometryUtils, 1);
};
