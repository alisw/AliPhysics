#ifndef ALITPCMISALIGNER_H
#define ALITPCMISALIGNER_H

#include "AliMisAligner.h"

// Class building the alignment objects for TPC in the three
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

class AliTPCMisAligner : public AliMisAligner {

    public:
	AliTPCMisAligner();
	TClonesArray* MakeAlObjsArray();
	AliCDBMetaData* GetCDBMetaData() const;

    private:
	ClassDef(AliTPCMisAligner,0);
};

#endif

