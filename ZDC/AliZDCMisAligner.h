#ifndef ALI_ZDC_MISALIGNER_H
#define ALI_ZDC_MISALIGNER_H

#include "TNamed.h"
#include "TString.h"
#include "AliMisAligner.h"
#include "AliCDBMetaData.h"

class TClonesArray;
class AliCDBManager;

  // Create a TClonesArray of misalignment objects for ZDC
  // of ideal/residual/full type according to request by
  // the steering macro

class AliZDCMisAligner : public AliMisAligner {

  public:
    AliZDCMisAligner();
    TClonesArray* MakeAlObjsArray();
    AliCDBMetaData* GetCDBMetaData() const;

  private:
    ClassDef(AliZDCMisAligner,0);
};

#endif

