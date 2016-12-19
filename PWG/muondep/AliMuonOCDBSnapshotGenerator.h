#ifndef ALIMUONOCDBSNAPSHOTGENERATOR_H
#define ALIMUONOCDBSNAPSHOTGENERATOR_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMuonOCDBSnapshotGenerator : public TObject
{
    public:
        AliMuonOCDBSnapshotGenerator(Int_t runNumber, const char* localOCDBPath="local://./OCDB", const char* sourceOCDB="raw://");

        Bool_t CreateLocalOCDBWithDefaultObjects(Bool_t overwrite=kFALSE);

        Bool_t CreateSnapshot(Int_t mode, const char* snapshotName,
                const char* alignSpecificStorage);
        
        Bool_t PopulateLocalOCDBWithPatchedRecoParam();

        TString LocalOCDBPath() const { return fLocalOCDBPath; }
        TString SourceOCDBPath() const { return fSourceOCDBPath; }
        Int_t RunNumber() const { return fRunNumber; }

    private:

        void DefaultSpecificStorage(Int_t mode, const char* alignSpecificStorage);

        Bool_t ExcludeFromSnapshot(const char* path);

        Int_t fRunNumber;
        TString fLocalOCDBPath;
        TString fSourceOCDBPath;

        ClassDef(AliMuonOCDBSnapshotGenerator,0)
};

#endif

