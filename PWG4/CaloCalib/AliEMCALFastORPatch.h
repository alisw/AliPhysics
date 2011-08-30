#ifndef ALIEMCALFASTORPATCH_H
#define ALIEMCALFASTORPATCH_H

// $Id$

class TArrayF;
class TArrayI;
class AliESDCaloTrigger;

#include <TObject.h>

class AliEMCALFastORPatch : public TObject {
  
public:
  AliEMCALFastORPatch();
  AliEMCALFastORPatch(Int_t aid, Int_t size = 4);
  virtual ~AliEMCALFastORPatch();
  
public:
  Int_t     GetAbsId()      const { return fAbsId;  }
  void      SetAbsId(Int_t aid)   { fAbsId = aid;   }
  
  Float_t   GetTotalAmplitude() const;
  Float_t   GetFastORamplitude(Int_t i) const;
  Int_t     GetNumberOfFastOR() const;
  Int_t     GetFastORrow(Int_t i) const;
  Int_t     GetFastORcolumn(Int_t i) const;
  Bool_t    Contains(Int_t row, Int_t col) const;
  Bool_t    AddFastORat(AliESDCaloTrigger* f, Int_t i);
  Bool_t    AddFastORat(Float_t amp, Int_t gCol, Int_t gRow, Int_t i);
  void      RemoveFastORat(Int_t i);
  void      Expand(Int_t size);
  
protected:
  Int_t          fAbsId;              // Unique ID of the FastOR patch (trigger)
  TArrayF       *fFastORamplitudes;   // Array containing amplitudes of the FastORs
  TArrayI       *fFastORcolumns;      // Array containing column position of the FastORs
  TArrayI       *fFastORrows;         // Array containing row position of the FastORs

private:
  AliEMCALFastORPatch(const AliEMCALFastORPatch&);            // not implemented
  AliEMCALFastORPatch &operator=(const AliEMCALFastORPatch&); // not implemented
  
  ClassDef(AliEMCALFastORPatch, 1);
};
#endif //ALIEMCALFASTORPATCH_H
