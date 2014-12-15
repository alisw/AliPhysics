#ifndef ALIPMDMISALIGNER_H
#define ALIPMDMISALIGNER_H

#include "AliMisAligner.h"

// Class building the alignment objects for PMD in the three
// canonical scenarios "ideal", "residual" and "full".
// It derives from AliMisAligner, thus providing the methods
// MakeAlObjsArray (builds and returns the array of alignment objects)
// and GetCDBMetaData (returns the metadata for the OCDB entry)
//

class TString;
class TNamed;
class TClonesArray;
class AliCDBManager;
class AliCDBMetaData;

class AliPMDMisAligner : public AliMisAligner {

    public:
	AliPMDMisAligner();
	TClonesArray* MakeAlObjsArray();
	AliCDBMetaData* GetCDBMetaData() const;

    private:
	ClassDef(AliPMDMisAligner,0);
};

#endif

