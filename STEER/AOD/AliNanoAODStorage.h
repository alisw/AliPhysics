#ifndef _ALINANOAODSTORAGE_H_
#define _ALINANOAODSTORAGE_H_


//-------------------------------------------------------------------------
//  AliNanoAODStorage
//
//  Implements the storage for special AOD classes 
//
//-------------------------------------------------------------------------
#include "TObject.h"
#include "TString.h"


class AliNanoAODStorage  {

public:
  AliNanoAODStorage():fNVars(0), fNVarsString(0), fVars(0), fVarsString(0) {;}
  virtual ~AliNanoAODStorage() {fVars.clear();};

  AliNanoAODStorage& operator=(const AliNanoAODStorage& sto);

  void AllocateInternalStorage(Int_t size);
  void AllocateInternalStorage(Int_t size, Int_t sizeString);
  
  static Int_t GetStringParameters(const TString varListHeader);
  
  void SetVar(Int_t index, Double_t var) { 
    if(index>=0 && index < fNVars)  fVars[index] = var;
    else Complain(index);
  }
  void SetVarString(Int_t index, TString var) { 
    if(index>=0 && index < fNVarsString)  fVarsString[index] = var;
    else Complain(index);
  }
  Double_t GetVar(Int_t index)  const {
    if(index>=0 && index < fNVars) return fVars[index]; 
    Complain(index);
    return 0;}

  TString GetVarString(Int_t index)  const {
    if(index>=0 && index < fNVarsString) return fVarsString[index]; 
    Complain(index);
    return "";}
  
  virtual const char *ClassName() const {return "AliNanoAODStorage";} // Needed for alifatal. The class cannot inherit from TObject, otherwise we mix (and mess up) inheritance in the track, as it also inherits from AliVTrack which inherits from TObject.

protected:

  Int_t    fNVars;     // Number of kimematic variables, set by constructor
  Int_t    fNVarsString;     // Number of string variables, set by constructor
  std::vector<Double32_t> fVars; // Array of kinematic vars. Here we use an STL vector because it produces ~5% smaller files. It may be aslo splittable
  std::vector<TString> fVarsString; // Array of string vars. Use same structure as for fVars

  ClassDef(AliNanoAODStorage, 2)
private:
  void Complain(Int_t index) const;
};



#endif /* _ALINANOAODSTORAGE_H_ */
