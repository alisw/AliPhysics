#ifndef _ALINANOAODSTORAGE_H_
#define _ALINANOAODSTORAGE_H_


/// \class AliNanoAODStorage
/// \brief AliNanoAODStorage
///
/// Implements the storage for special AOD classes
///

#include "TObject.h"
#include "TString.h"


class AliNanoAODStorage  {

public:
  AliNanoAODStorage():fNVars(0), fNVarsInt(0), fVars(0), fVarsInt(0) {;}
  virtual ~AliNanoAODStorage() {fVars.clear();};

  AliNanoAODStorage& operator=(const AliNanoAODStorage& sto);

  void AllocateInternalStorage(Int_t size);
  void AllocateInternalStorage(Int_t size, Int_t sizeInt);
  
  static Int_t GetIntParameters(const TString varListHeader);
  
  void SetVar(Int_t index, Double_t var) { 
    if(index>=0 && index < fNVars)  fVars[index] = var;
    else Complain(index);
  }
  void SetVarInt(Int_t index, UInt_t var) { 
    if(index>=0 && index < fNVarsInt)  fVarsInt[index] = var;
    else Complain(index);
  }
  Double_t GetVar(Int_t index)  const {
    if(index>=0 && index < fNVars) return fVars[index]; 
    Complain(index);
    return 0;
  }

  UInt_t GetVarInt(Int_t index)  const {
    if(index>=0 && index < fNVarsInt) return fVarsInt[index]; 
    Complain(index);
    return 0;
  }
  
  virtual const char *ClassName() const {return "AliNanoAODStorage";} // Needed for alifatal. The class cannot inherit from TObject, otherwise we mix (and mess up) inheritance in the track, as it also inherits from AliVTrack which inherits from TObject.

protected:

  Int_t    fNVars;     ///< Number of kimematic variables, set by constructor
  Int_t    fNVarsInt;     ///< Number of int variables, set by constructor
  std::vector<Double32_t> fVars; // Array of kinematic vars. Here we use an STL vector because it produces ~5% smaller files. It may be aslo splittable
  std::vector<UInt_t> fVarsInt; // Array of int vars. Use same structure as for fVars, int has to be used for a bitwise variable describing the fired trigger particles, because this is bitwise coded

  ClassDef(AliNanoAODStorage, 4)
private:
  void Complain(Int_t index) const;
};



#endif /* _ALINANOAODSTORAGE_H_ */
