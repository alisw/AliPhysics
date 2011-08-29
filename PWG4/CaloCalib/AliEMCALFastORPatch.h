#ifndef ALIEMCALFASTORPATCH_H
#define ALIEMCALFASTORPATCH_H

// $Id$

class TObjArray;
class AliESDCaloTrigger;

#include <TObject.h>

class AliEMCALFastORPatch : public TObject {
  
public:
  AliEMCALFastORPatch();
  AliEMCALFastORPatch(Int_t aid, Int_t size = 4);
  virtual ~AliEMCALFastORPatch();
  
  Float_t GetTotalAmplitude();
  Float_t GetFastORamplitude(Int_t i);
  Int_t GetNumberOfFastOR();
  Bool_t AddFastORat(AliESDCaloTrigger* f, Int_t i);
  Bool_t AddFastORat(Float_t amp, Int_t gCol, Int_t gRow, Int_t i);
  void RemoveFastORat(Int_t i);
  Int_t GetAbsId();
  void SetAbsId(Int_t aid);
  Int_t GetFastORrow(Int_t i);
  Int_t GetFastORcolumn(Int_t i);
  void Expand(Int_t size);
  Bool_t Contains(Int_t row, Int_t col);
  
private:
  Float_t       *fFastORamplitudes;
  Int_t         *fFastORcolumns;
  Int_t         *fFastORrows;
  Int_t          fSize;
  Int_t          fAbsId;
  
  ClassDef(AliEMCALFastORPatch, 1);
};
#endif //ALIEMCALFASTORPATCH_H