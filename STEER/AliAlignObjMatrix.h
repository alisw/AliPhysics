#ifndef ALIALIGNOBJMATRIX_H
#define ALIALIGNOBJMATRIX_H

/**************************************************************************
 * AliAlignObjMatrix: derived alignment class storing alignment           *
 *   information for a single volume in form of TGeoHMatrix, which        *
 *   includes the information for a translation, a rotation and a scale   *
 *************************************************************************/
#include "TObject.h"
#include "TString.h"
#include "TGeoMatrix.h"

#include "AliAlignObj.h"

class AliAlignObjMatrix : public AliAlignObj {
 public:
  AliAlignObjMatrix();
  AliAlignObjMatrix(const char* volpath, UShort_t volUId, Double_t x, Double_t y, Double_t z, Double_t psi, Double_t theta, Double_t phi, Bool_t global) throw (const Char_t *);
  AliAlignObjMatrix(const char* volpath, UShort_t volUId, TGeoMatrix& m, Bool_t global) throw (const Char_t *);
  AliAlignObjMatrix(const AliAlignObj& theAlignObj);
  AliAlignObjMatrix& operator= (const AliAlignObj& theAlignObj);
  virtual ~AliAlignObjMatrix();
  
  //Setters
  virtual void SetTranslation(Double_t x, Double_t y, Double_t z);
  virtual void SetTranslation(const TGeoMatrix& m);
  virtual void SetRotation(Double_t psi, Double_t theta, Double_t phi);
  virtual Bool_t SetRotation(const TGeoMatrix& m);

  //Getters
  virtual void GetTranslation(Double_t* tr)  const;
  virtual Bool_t GetAngles(Double_t* angles)  const;
  virtual void GetMatrix(TGeoHMatrix& m) const;

  virtual AliAlignObj& Inverse() const;
  
 protected:
  TGeoHMatrix fMatrix; // Transformation matrix
  
  ClassDef(AliAlignObjMatrix, 1)
};

#endif
