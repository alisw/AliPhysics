#ifndef ALI_MISALIGNER_H
#define ALI_MISALIGNER_H

#include "TObject.h"
#include "TString.h"
#include "AliCDBMetaData.h"

class TClonesArray;
class AliCDBManager;

// Base class for creating a TClonesArray of simulated misalignment objects
// for a given subdetector of type ideal,residual or full
//

class AliMisAligner : public TObject {

  public:
    AliMisAligner();
    virtual TClonesArray* MakeAlObjsArray() =0;
    virtual AliCDBMetaData* GetCDBMetaData() const =0;
    void SetMisalType(const char* misalType)
    {
      fMisalType=misalType;
    }
    const char* GetMisalType() const
    {
      return fMisalType.Data();
    }

  protected:
    TString fMisalType;

  private:
    ClassDef(AliMisAligner,0);
};

#endif

