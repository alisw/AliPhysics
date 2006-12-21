#ifndef ALIALIGNOBJANGLES_H
#define ALIALIGNOBJANGLES_H

/*****************************************************************************
 * AliAlignObjAngles: derived alignment class storing alignment information  *
 *   for a single volume in form of three doubles for the translation        *
 *   and three doubles for the rotation expressed with the euler angles      *
 *   in the xyz-convention (http://mathworld.wolfram.com/EulerAngles.html),  *
 *   also known as roll, pitch, yaw. PLEASE NOTE THE ANGLES SIGNS ARE        *
 *   INVERSE WITH RESPECT TO THIS REFERENCE!!! In this way the representation*
 *   is fully consistent with the TGeo Rotation methods.                     *
 *****************************************************************************/
#include "TObject.h"
#include "TString.h"
#include "TGeoMatrix.h"

#include "AliAlignObj.h"

class AliAlignObjAngles : public AliAlignObj{
 public:
  AliAlignObjAngles();
  AliAlignObjAngles(const char* symname, UShort_t volUId, Double_t x, Double_t y, Double_t z, Double_t psi, Double_t theta, Double_t phi, Bool_t global) throw (const Char_t *);
  AliAlignObjAngles(const char* symname, UShort_t volUId, TGeoMatrix& m, Bool_t global) throw (const Char_t *);
  AliAlignObjAngles(const AliAlignObj& theAlignObj);
  AliAlignObjAngles& operator= (const AliAlignObj& theAlignObj);
  virtual ~AliAlignObjAngles();
  
  //Setters
  virtual void SetTranslation(Double_t x, Double_t y, Double_t z){
    fTranslation[0]=x; fTranslation[1]=y; fTranslation[2]=z;}
  virtual void SetTranslation(const TGeoMatrix& m);
  virtual void SetRotation(Double_t psi, Double_t theta, Double_t phi){
    fRotation[0]=psi; fRotation[1]=theta; fRotation[2]=phi;}
  virtual Bool_t SetRotation(const TGeoMatrix& m);
  
  //Getters
  virtual void GetTranslation(Double_t *tr)  const {
    tr[0] = fTranslation[0]; tr[1] = fTranslation[1]; tr[2] = fTranslation[2];}
  virtual Bool_t GetAngles(Double_t* angles)   const {
    angles[0] = fRotation[0]; angles[1] = fRotation[1];
    angles[2] = fRotation[2]; return kTRUE;}
  virtual void GetMatrix(TGeoHMatrix& m) const;

  virtual AliAlignObj& Inverse() const;
  
 protected:
  Double_t fTranslation[3]; // Translation vector
  Double_t fRotation[3]; // Roll-pitch-yaw angles
  
  ClassDef(AliAlignObjAngles, 1)
};

#endif
