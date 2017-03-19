#include "AliMuonCompactTreeMaker.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisMuonUtility.h"
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonCluster.h"
#include "AliESDMuonTrack.h"
#include "AliGeomManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMuonCompactEvent.h"
#include "AliMuonCompactMapping.h"
#include "AliVParticle.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include <cassert>

/// \ingroup compact
namespace {
    const Int_t kCANNOTSETUPOCDB = -1;
    const Int_t kCANNOTGETESDTREE = -2;
    const Int_t kCANNOTOPENESDFILE = -3; 
    UInt_t ENCODE(Int_t a16, Int_t b16)
    {
        return ( ( a16 & 0x0000FFFF ) << 16 ) | ( b16 & 0x0000FFFF );
    }

}

ClassImp(AliMuonCompactTreeMaker)

AliMuonCompactTreeMaker::AliMuonCompactTreeMaker(const char* ocdbPath): 
    AliAnalysisTaskSE("AliMuonCompactTreeMaker"), fOCDBPath(ocdbPath), fRunNumber(-1), fGeometryTransformer(nullptr), fOutputTree(nullptr), fCompactEvent()
{
    DefineOutput(1,TTree::Class());
}

AliMuonCompactTreeMaker::AliMuonCompactTreeMaker(): 
    AliAnalysisTaskSE("AliMuonCompactTreeMaker"), fOCDBPath(), fRunNumber(-1), fGeometryTransformer(nullptr), fOutputTree(nullptr), fCompactEvent()
{
}

void AliMuonCompactTreeMaker::CleanupOCDB()
{
    AliCDBManager::Instance()->ClearCache();
    AliGeomManager::Destroy();
    delete fGeometryTransformer;
    fGeometryTransformer = nullptr;
}

Bool_t AliMuonCompactTreeMaker::SetupOCDB(const char* alignSpecificStorage)
{
    AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.c_str());
    AliCDBManager::Instance()->SetRun(fRunNumber);
    AliMpCDB::LoadAll();
    
    AliGeomManager::LoadGeometry();

    AliCDBManager::Instance()->SetSpecificStorage("MUON/Align/Data",alignSpecificStorage);

    if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
        return kFALSE;
    }
    fGeometryTransformer = new AliMUONGeometryTransformer;
    fGeometryTransformer->LoadGeometryData();

    AliMuonCompactMapping::GetCompactMapping(fOCDBPath.c_str(),fRunNumber);

    return kTRUE;
}

void AliMuonCompactTreeMaker::GetClusterLocation(Int_t detElemId,
        Double_t xg, 
        Double_t yg, 
        Double_t zg,
        Int_t& bendingManuAbsIndex,
        Int_t& nonBendingManuAbsIndex)
{
    /// Get the pad under the center of the cluster
    AliMuonCompactMapping* cm = AliMuonCompactMapping::GetCompactMapping();

    Double_t x,y,z;

    fGeometryTransformer->Global2Local(detElemId,
            xg,yg,zg,x,y,z);

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    const AliMpVSegmentation* segB = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,de->GetCathodType(AliMp::kBendingPlane));
    const AliMpVSegmentation* segNB = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,de->GetCathodType(AliMp::kNonBendingPlane));

    AliMpPad padB = segB->PadByPosition(x,y);
    AliMpPad padNB = segNB->PadByPosition(x,y);

    Int_t manuBending = padB.GetManuId();
    Int_t manuNonBending = padNB.GetManuId();

    bendingManuAbsIndex=-1;
    nonBendingManuAbsIndex=-1;

    if ( padB.IsValid() )
    { 
        bendingManuAbsIndex = cm->FindManuAbsIndex(de->GetId(),manuBending);
    }

    if ( padNB.IsValid() )
    {
        nonBendingManuAbsIndex = cm->FindManuAbsIndex(de->GetId(),manuNonBending);
    }
}
void AliMuonCompactTreeMaker::ConvertEvent(AliESDEvent& esd, AliMuonCompactEvent& compactEvent)
{
    compactEvent.mTracks.clear();

    for ( Int_t i = 0; i < esd.GetNumberOfMuonTracks(); ++i )
    {
        AliESDMuonTrack* track = esd.GetMuonTrack(i);

        if  (!track->ContainTrackerData()) continue;

        if (track->GetMatchTrigger()<2) continue;

        if (track->Eta() < -4 || track->Eta() > -2.5 ) continue;

        if (track->GetRAtAbsorberEnd() < 17.5 || track->GetRAtAbsorberEnd() > 89.0 ) continue;

        AliMuonCompactTrack compactTrack(track->Px(),track->Py(),track->Pz());

        for ( Int_t j = 0; j < track->GetNClusters(); ++j )
        {
            UInt_t id = track->GetClusterId(j);
            AliESDMuonCluster* cluster = esd.FindMuonCluster(id);
            Int_t b,nb;
            GetClusterLocation(cluster->GetDetElemId(),
                    cluster->GetX(),
                    cluster->GetY(),
                    cluster->GetZ(),
                    b,
                    nb);
            if (b>=0 || nb>=0)
            {
                AliMuonCompactCluster cl(b,nb);
                compactTrack.mClusters.push_back(cl);
            }
            else
            {
                std::cout << "Got no manu for this cluster ?" << std::endl;
                cluster->Print();
            }
        }

        compactEvent.mTracks.push_back(compactTrack);
    }
}

void AliMuonCompactTreeMaker::UserExec(Option_t*)
{
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());

    assert(esd!=nullptr);

    fCompactEvent.Clear();

    // if ( esd.GetNumberOfMuonTracks() >= 2 )
    {
        ConvertEvent(*esd,fCompactEvent);
    }

    if (MCEvent())
    {
        Int_t nprimary(0);

        for ( Int_t i = 0; i < MCEvent()->GetNumberOfTracks(); ++i )
        {
            AliVParticle* part = MCEvent()->GetTrack(i);

            if ( AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) && part->GetMother() == -1 ) 
            {
                nprimary++;
                fCompactEvent.SetInput(part->Pt(),part->Y());
            }
        }

        if (nprimary!=1)
        {
            AliFatal(Form("Wrong nprimary=%d",nprimary));
        }
    }

    fOutputTree->Fill();
    PostData(1,fOutputTree);
}

void AliMuonCompactTreeMaker::UserCreateOutputObjects()
{
    OpenFile(1);

    fOutputTree = new TTree("compactevents","a tree with compacted tracks");
    fOutputTree->Branch("event",&fCompactEvent);

    PostData(1,fOutputTree);
}

void AliMuonCompactTreeMaker::NotifyRun()
{
    if ( fCurrentRunNumber != fRunNumber ) 
    {
        CleanupOCDB();
        fRunNumber = fCurrentRunNumber;

        // get the alignment (specific) storage
        TTree* tree = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetTree();
        TList* userInfo = tree->GetUserInfo();
        TMap* cdbMap = static_cast<TMap*>(userInfo->FindObject("cdbMap"));
        TObjString* alignSpecificStorage = static_cast<TObjString*>(cdbMap->GetValue("MUON/Align/Data"));
        if (!alignSpecificStorage)
        {
            // no specific storage, take it from the regular list
            alignSpecificStorage = static_cast<TObjString*>(cdbMap->FindObject("default"));
        }
        if ( !SetupOCDB(alignSpecificStorage->String().Data()) )
        {
            AliFatal("Cannot setup OCDB !");
        }
    }
}
