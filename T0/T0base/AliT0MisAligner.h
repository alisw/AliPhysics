#ifndef ALIT0MISALIGNER_H
#define ALIT0MISALIGNER_H

#include "AliMisAligner.h"

// Class building the alignment objects for T0 in the three
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

class AliT0MisAligner : public AliMisAligner {

    public:
	AliT0MisAligner();
	TClonesArray* MakeAlObjsArray();
	AliCDBMetaData* GetCDBMetaData() const;

    private:
	ClassDef(AliT0MisAligner,0);
};

#endif

