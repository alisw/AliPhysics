#ifndef _ALINANOAODSTORAGE_H_
#define _ALINANOAODSTORAGE_H_


//-------------------------------------------------------------------------
//  AliNanoAODStorage
//
//  Implements the storage for special AOD classes 
//
//-------------------------------------------------------------------------
#include <iostream>
#include "TObject.h"
#include "AliLog.h"


class AliNanoAODStorage  {

public:
  AliNanoAODStorage():fNVars(0), fVars(0) {;}
  virtual ~AliNanoAODStorage() {;};

  AliNanoAODStorage& operator=(const AliNanoAODStorage& sto);

  void AllocateInternalStorage(Int_t size);
  void SetVar(Int_t index, Double_t var) { 
    if(index>=0 && index < fNVars)  fVars[index] = var;  
    else AliFatal(Form("Variable %d not included in this special aod", index));
  }
  Double_t GetVar(Int_t index)  const {
    if(index>=0 && index < fNVars) return fVars[index]; 
    AliFatal(Form("Variable %d not included in this special aod", index));
    return 0;}
  virtual const char *ClassName() const {return "AliNanoAODStorage";} // Needed for alifatal. The class cannot inherit from TObject, otherwise we mix (and mess up) inheritance in the track, as it also inherits from AliVTrack which inherits from TObject.

protected:

  Int_t    fNVars;     // Number of kimematic variables, set by constructor
  std::vector<Double32_t> fVars; // Array of kinematic vars. Here we use an STL vector because it produces ~5% smaller files. It may be aslo splittable

  ClassDef(AliNanoAODStorage, 1)
};



#endif /* _ALINANOAODSTORAGE_H_ */
