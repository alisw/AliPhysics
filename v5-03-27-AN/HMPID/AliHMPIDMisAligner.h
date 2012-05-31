#ifndef ALIHMPIDMISALIGNER_H
#define ALIHMPIDMISALIGNER_H

#include "AliMisAligner.h"

// Class building the alignment objects for HMPID in the three
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

class AliHMPIDMisAligner : public AliMisAligner {

    public:
	AliHMPIDMisAligner();
	TClonesArray* MakeAlObjsArray();
	AliCDBMetaData* GetCDBMetaData() const;

    private:
	ClassDef(AliHMPIDMisAligner,0);
};

#endif

