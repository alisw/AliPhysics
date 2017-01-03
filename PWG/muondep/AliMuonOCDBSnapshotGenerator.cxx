#include "AliMuonOCDBSnapshotGenerator.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliLog.h"
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMpCDB.h"
#include <cassert>
#include "Riostream.h"
#include "TObjArray.h"
#include "TSystem.h"


/// \cond CLASSIMP
ClassImp(AliMuonOCDBSnapshotGenerator)
/// \endcond

AliMuonOCDBSnapshotGenerator::AliMuonOCDBSnapshotGenerator(Int_t runNumber, const char* localOCDBPath, const char* sourceOCDBPath)
    : TObject(), fRunNumber(runNumber), fLocalOCDBPath(localOCDBPath), fSourceOCDBPath(sourceOCDBPath)
{

}

/// Generate a local OCDB containing some default MUON/*/* objects.
/// 
/// The actual creation of the default objects is delegated to the AliMUONCDB class.
/// 
/// \param overwrite whether or not existing local OCDB should be overwritten.
///
Bool_t AliMuonOCDBSnapshotGenerator::CreateLocalOCDBWithDefaultObjects(Bool_t overwrite)
{
    Bool_t shouldDoIt = kFALSE;

    if ( overwrite )
    {
        shouldDoIt = kTRUE;
    }
    else
    {
        if ( gSystem->AccessPathName("OCDB/MUON/Calib/HV/Run0_999999999_v0_s0.root") || 
                gSystem->AccessPathName("OCDB/MUON/Calib/Config/Run0_999999999_v0_s0.root") || 
                gSystem->AccessPathName("OCDB/MUON/Calib/LV/Run0_999999999_v0_s0.root") || 
                gSystem->AccessPathName("OCDB/MUON/Calib/OccupancyMap/Run0_999999999_v0_s0.root") || 
                gSystem->AccessPathName("OCDB/MUON/Calib/Pedestals/Run0_999999999_v0_s0.root") ||
                gSystem->AccessPathName("OCDB/MUON/Calib/RejectList/Run0_999999999_v0_s0.root") ||
                gSystem->AccessPathName("OCDB/MUON/Calib/RecoParam/Run0_999999999_v0_s0.root") )
        {
            // we don't have all our objects, should generate them
            // 
            shouldDoIt = kTRUE;
        }
    }

    if ( !shouldDoIt )
    {
        AliWarning(Form("Local default MUON Calib objects already there under directory %s/MUON/Calib : not redoing them. Erase them first if you'd prefer to re-create them.",gSystem->WorkingDirectory()));
        return kTRUE;
    }

    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage(SourceOCDBPath().Data());
    man->SetSpecificStorage("MUON/Calib/HV",LocalOCDBPath().Data());
    man->SetSpecificStorage("MUON/Calib/LV",LocalOCDBPath().Data());
    man->SetSpecificStorage("MUON/Calib/Config",LocalOCDBPath().Data());
    man->SetSpecificStorage("MUON/Calib/OccupancyMap",LocalOCDBPath().Data());
    man->SetSpecificStorage("MUON/Calib/Pedestals",LocalOCDBPath().Data());
    man->SetSpecificStorage("MUON/Calib/RejectList",LocalOCDBPath().Data());
    man->SetRun(RunNumber());

    AliMpCDB::LoadAll();

    AliMUONCDB::WriteHV(true,0);
    AliMUONCDB::WriteLV(true,0);
    AliMUONCDB::WriteOccupancyMap(true,0);
    AliMUONCDB::WriteConfig(0);
    AliMUONCDB::WritePedestals(true,0);
    AliMUONCDB::WriteRejectList(true,0);

    PopulateLocalOCDBWithPatchedRecoParam();

    return kTRUE; 
}

/// Create a snapshot 
///
/// \param mode @see DefaultSpecificStorage()
/// \param snapshotName name of the snapshot file (typicall OCDB_rec.root or OCDB_sim.root)
/// \param alignSpecificStorage the specific storage location to be used for the (mis)alignement file
///
Bool_t AliMuonOCDBSnapshotGenerator::CreateSnapshot(Int_t mode, const char* snapshotName,
        const char* alignSpecificStorage)
{
    if ( gSystem->AccessPathName(gSystem->DirName(snapshotName)) )
    {
       AliInfo(Form("Directory %s does not exist. Attempting to create it",
               gSystem->DirName(snapshotName)));
       if (gSystem->MakeDirectory(gSystem->DirName(snapshotName))) return kFALSE;
    }
  
  AliInfo(Form("Snapshot %s for mode %d with align=%s",snapshotName,mode,alignSpecificStorage));

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(SourceOCDBPath().Data());
  man->SetRun(RunNumber());

  DefaultSpecificStorage(mode,alignSpecificStorage);

  AliCDBStorage *defStorage = man->GetDefaultStorage();
  TObjArray* arrCDBID = defStorage->GetQueryCDBList();

  const TMap* stMap = man->GetStorageMap();
  man->SetCacheFlag(kTRUE);

  AliCDBId* cdbID = 0;
  TIter nxt(arrCDBID);
  while ((cdbID=static_cast<AliCDBId*>(nxt()))) { 
          // loop over default storage
    TString path = cdbID->GetPath();
    if (stMap->GetValue(path)) {
        // already defined in the specific storage
      continue; 
    }
    if (ExcludeFromSnapshot(path)) continue;
    man->Get(path.Data());
  }
  // 
  TIter nextSt(stMap);
  TObjString *str;
  while ((str=static_cast<TObjString*>(nextSt()))) { // exclusion is not applied to specific objects
    TString calType = str->GetString();
    if (calType=="default") continue;
    if (ExcludeFromSnapshot(calType)) continue;
    man->Get(calType.Data());
 
  }
  man->DumpToSnapshotFile(snapshotName,kFALSE);
  return kTRUE;
}

/// Set some specific storages, depending on the intended usage of the snapshot (simulation or reconstruction)
///
/// \param mode 0 for simulation and 1 for reconstruction
/// \param alignSpecificStorage the storage to be used for the (mis)alignement 
///
void AliMuonOCDBSnapshotGenerator::DefaultSpecificStorage(Int_t mode, const char* alignSpecificStorage)
{
   assert(mode==0 || mode==1);
   AliCDBManager* man = AliCDBManager::Instance();

   // copied from AliDPG/MC/OCDBConfig.C, but stripping out
   // non relevant parts (e.g. TPC)
   
  const Char_t *Ideal    = "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/";
  const Char_t *Residual = "alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/";
  const Char_t *Full     = "alien://Folder=/alice/simulation/2008/v4-15-Release/Full/";
  
  // DEFAULT SPECIFIC OBJECTS 
  const Char_t *SpecificStorageList[][3] = {
    // path                    sim       rec
    //
    // ITS
    "ITS/Align/Data",          Ideal,    Residual, // ok
    "ITS/Calib/SPDSparseDead", Full,     Residual, // ok ?
    // ZDC
    "ZDC/Align/Data",          Ideal,    Ideal     // ok
  };
  const Int_t nSpecificStorages = sizeof(SpecificStorageList) / (3 * sizeof(Char_t *));

  // set specific storages
  for (Int_t isto = 0; isto < nSpecificStorages; isto++) {
    if (SpecificStorageList[mode+1][isto]) {
      man->SetSpecificStorage(SpecificStorageList[isto][0], SpecificStorageList[isto][mode+1]);
    }
  }

  // special hacks to compute "ideal" simulations for Quick Acc x Eff

  TString sAlignSpecificStorage(alignSpecificStorage);

  sAlignSpecificStorage.ReplaceAll("\"","");

  AliInfo(Form("alignSpecificStorage=%s",sAlignSpecificStorage.Data()));

  man->SetSpecificStorage("MUON/Align/Data",sAlignSpecificStorage.Data());
  man->SetSpecificStorage("MUON/Calib/Pedestals",LocalOCDBPath().Data());
  man->SetSpecificStorage("MUON/Calib/OccupancyMap",LocalOCDBPath().Data());
  man->SetSpecificStorage("MUON/Calib/HV",LocalOCDBPath().Data());
  man->SetSpecificStorage("MUON/Calib/LV",LocalOCDBPath().Data());
  man->SetSpecificStorage("MUON/Calib/Config",LocalOCDBPath().Data());
}

/// Whether or not to exclude a given path from the snapshot
/// The idea here is to make the snapshot as slim as possible for muon simulations, so
/// we remove everything not strictly needed for the purpose, like TPC objects etc...
Bool_t AliMuonOCDBSnapshotGenerator::ExcludeFromSnapshot(const char* path)
{
    TString spath(path);

    // can not exclude align data
    if ( spath.Contains("/Align/Data") ) return kFALSE;
    if ( spath.Contains("MeanVertexTPC")) return kFALSE;
    if ( spath.Contains("TPC/Calib/RecoParam")) return kFALSE;
    if ( spath.Contains("HLT/Calib/esdLayout")) return kFALSE;
    if ( spath.Contains("GRP/Calib/MeanVertex")) return kFALSE;

    if ( spath.Contains("ACORDE") ) return kTRUE;
    if ( spath.Contains("EMCAL") ) return kTRUE;
    if ( spath.Contains("FMD") ) return kTRUE;
    if ( spath.Contains("HLT") ) return kTRUE;
    if ( spath.Contains("HMPID") ) return kTRUE;
    if ( spath.Contains("PHOS") ) return kTRUE;
    if ( spath.Contains("PMD") ) return kTRUE;
    if ( spath.Contains("TOF") ) return kTRUE;
    if ( spath.Contains("TPC") ) return kTRUE;
    if ( spath.Contains("TRD") ) return kTRUE;
    if ( spath.Contains("AD") ) return kTRUE;
    if ( spath.Contains("T0") ) return kTRUE;
    if ( spath.Contains("VZERO") ) return kTRUE;
    if ( spath.Contains("ZDC") ) return kTRUE;

    // MUON objects that are no longer in use
    if ( spath.Contains("MUON/Calib/Gains")) return kTRUE;
    if ( spath.Contains("MUON/Calib/Capacitances")) return kTRUE;
    if ( spath.Contains("MUON/Calib/Neighbours")) return kTRUE;
    if ( spath.Contains("MUON/Calib/GlobalTriggerBoardMasks")) return kTRUE;
    if ( spath.Contains("MUON/Calib/RegionalTriggerBoardMasks")) return kTRUE;

    return kFALSE;
}

/// Read the MUON/Calib/RecoParam object from the source OCDB, change the
/// PadGoodnessMask, and write it to the local OCDB.
Bool_t AliMuonOCDBSnapshotGenerator::PopulateLocalOCDBWithPatchedRecoParam()
{
    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage(SourceOCDBPath().Data());
    man->SetRun(RunNumber());
    // patch the recoparam to _not_ cut on any of the above

    AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam");

    TObject* o = entry->GetObject();

    if ( o->IsA() != TObjArray::Class() ) 
    {
        AliFatal("This code only works with TObjArray recoparams. Sorry");
        return kFALSE;
    }

    TObjArray* array = static_cast<TObjArray*>(o);
    for ( Int_t i = 0; i <= array->GetLast(); ++i ) 
    {
        AliDetectorRecoParam* p = static_cast<AliDetectorRecoParam*>(array->At(i));
        if (AliRecoParam::Convert(p->GetEventSpecie())==AliRecoParam::kLowMult ||
            AliRecoParam::Convert(p->GetEventSpecie())==AliRecoParam::kHighMult)
        {
            std::cout << Form("array[%d]=%s %s %s",i,
                    p ? p->ClassName() : "",
                    p ? AliRecoParam::GetEventSpecieName(AliRecoParam::Convert(p->GetEventSpecie())) :"",
                    p ? ( p->IsDefault() ? "default" : "") : "" ) << std::endl;
            p->Print("");
            AliMUONRecoParam* rp = dynamic_cast<AliMUONRecoParam*>(p);
            if (!rp) 
            {
                AliFatal("OUPS. OUPS");
                return kFALSE;
            }
            rp->SetPadGoodnessMask(0);
        }
    }

    man->SetSpecificStorage("MUON/Calib/RecoParam",LocalOCDBPath().Data());

    AliMUONCDB::WriteToCDB(array, "MUON/Calib/RecoParam", 0, AliCDBRunRange::Infinity(), "reconstruction parameters for MUON, patched for no cut", "AliMuonAccEffSubmitter");

    return kTRUE;
}


