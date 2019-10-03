#ifndef ALIMUONOCDBSNAPSHOTGENERATOR_H
#define ALIMUONOCDBSNAPSHOTGENERATOR_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

/// \ingroup pwg_muondep_submitter
/// \class AliMuonOCDBSnapshotGenerator
/// \brief Generate a custom OCDB snapshot for pseudo-ideal muon simulations
/// 
/// The snapshot generation is a two step process. First a local OCDB is generated
/// with default objects by CreateLocalOCDBWithDefaultObjects() and PopulateLocalOCDBWithPatchedRecoParam()
/// Then the snapshot itself is created using the objects from that local OCDB 
/// and from the source OCDB. The list of objects to be taken from the source OCDB
/// is given by what is available for the requested run in the requested source OCDB,
/// minus the objects specified by ExcludeFromSnapshot()
///
/// \author Laurent Aphecetche, Subatech

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

        Int_t fRunNumber; ///< reference run to be used for the objects taken from the sourceOCDB
        TString fLocalOCDBPath; ///< path of the temporary local OCDB that will be created as one source of the snapshot
        TString fSourceOCDBPath; ///< the reference OCDB for most of the objects of the snapshot

/// \cond CLASSIMP
        ClassDef(AliMuonOCDBSnapshotGenerator,0);
/// \endcond
};

#endif

