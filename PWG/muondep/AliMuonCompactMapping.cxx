#include "AliMuonCompactMapping.h"
#include "AliCDBManager.h"
#include "AliMpCDB.h"
#include "AliMpDEIterator.h"
#include <vector>
#include "AliMpDEManager.h"
#include "AliMpManuIterator.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include <cassert>
#include <algorithm>

/// \ingroup compact

namespace {
    Int_t DECODELOW(UInt_t e)
    {
        return e & 0x0000FFFF;
    }

    Int_t DECODEHIGH(UInt_t e)
    {
        return (e & 0xFFFF0000) >> 16;
    }
    UInt_t ENCODE(Int_t a16, Int_t b16)
    {
        return ( ( a16 & 0x0000FFFF ) << 16 ) | ( b16 & 0x0000FFFF );
    }
}

Int_t AliMuonCompactMapping::GetDetElemIdFromAbsManuIndex(Int_t index) const
{
    return GetDetElemIdFromAbsManuId(AbsManuId(index));
}

Int_t AliMuonCompactMapping::AbsManuId(UInt_t index) const 
{
    if ( index>mManuIds.size())
    {
        std::cout << index << " > " << mManuIds.size()
            << std::endl;
    }
    return mManuIds[index];
}

Int_t AliMuonCompactMapping::GetManuIdFromAbsManuId(UInt_t absManuId) const
{
    return DECODELOW(absManuId);
}

Int_t AliMuonCompactMapping::GetDetElemIdFromAbsManuId(UInt_t absManuId) const
{
    return DECODEHIGH(absManuId);
}

Int_t AliMuonCompactMapping::GetNofPadsFromAbsManuIndex(Int_t index) const
{
    return mNpads[index];
}

Int_t AliMuonCompactMapping::FindManuAbsIndex(Int_t detElemId, Int_t manuId) const
{
    UInt_t encoded = ENCODE(detElemId,manuId);
    auto p = mManuMap.find(encoded);
    return p->second;
}

AliMuonCompactMapping* AliMuonCompactMapping::GetCompactMapping(const char* ocdbPath, Int_t runNumber)
{
    static AliMuonCompactMapping* cm(0x0);

    if (!cm)
    {
        cm = new AliMuonCompactMapping;

        AliCDBManager* man = AliCDBManager::Instance();
        if (!man->IsDefaultStorageSet())
        {
            man->SetDefaultStorage(ocdbPath);
            man->SetRun(runNumber);
        }
        AliMpCDB::LoadAll();

        // first get a sorted list of the detElemId

        AliMpDEIterator deit;
        deit.First();
        std::vector<int> deids;
        while (!deit.IsDone())
        {
            Int_t detElemId = deit.CurrentDEId();

            if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger )
            {
                deids.push_back(detElemId);
            }
            deit.Next();
        }

        std::sort(deids.begin(),deids.end());

        // then for each DE get one sorted list of manuIds
        std::vector<int>::size_type totalNofManus(0);

        for ( std::vector<int>::size_type i = 0; i < deids.size(); ++i )
        {
            AliMpManuIterator it;
            Int_t detElemId, manuId;
            std::vector<int> bendingManuids;
            std::vector<int> nonBendingManuids;

            // get the list of manu ids on both planes
            // for this detection element
            while ( it.Next(detElemId,manuId) )
            {
                if ( detElemId == deids[i] )
                {
                    if ( manuId >= 1024 )
                    {
                        nonBendingManuids.push_back(manuId);
                    }
                    if ( manuId < 1024 )
                    {
                        bendingManuids.push_back(manuId);
                    }
                }
            }

            detElemId = deids[i];

            // insure manuids are sorted (should be the case
            // already, though)
            std::sort(bendingManuids.begin(),bendingManuids.end());
            std::sort(nonBendingManuids.begin(),nonBendingManuids.end());

            // add those manus to the global array of AbsManuIds
            // and update the relevant "links"
            

            std::vector<int>& allManuOfThisDE = bendingManuids;
            allManuOfThisDE.insert(allManuOfThisDE.end(),
                    nonBendingManuids.begin(),
                    nonBendingManuids.end());
          
            UInt_t ix = cm->mManuIds.size();

            for ( std::vector<int>::size_type i = 0; i < allManuOfThisDE.size(); ++i )
            {
                Int_t manuId = allManuOfThisDE[i];
                UInt_t encodedManu = ENCODE(detElemId,manuId);
                cm->mManuMap[encodedManu]=ix;
                cm->mManuIds.push_back(encodedManu);

                // get the number of pads per manu (quite usefull in 
                // some instances, e.g. to compute occupancies...)
                // 
                AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
                cm->mNpads.push_back(de->NofChannelsInManu(manuId));
                ++ix;
            }

            totalNofManus += allManuOfThisDE.size();

            // std::cout << Form("DE %04d (%3d) ",deids[i],
            //         allManuOfThisDE.size());
            //
            // for ( std::vector<int>::size_type j = 0; 
            //         j < allManuOfThisDE.size(); ++j )
            // {
            //     std::cout << Form("%04d ",allManuOfThisDE[j]);
            // }
            // std::cout << std::endl;
        }

        // std::cout << "Total number of manus : " << totalNofManus << std::endl;
        assert(totalNofManus==16828);
        assert(cm->mNpads.size()==totalNofManus);
    }
    return cm;
}

void AliMuonCompactMapping::GenerateStaticOffsets()
{
    // generate the detection element ids arrays
    // that are needed in some other part of the code
    
    Int_t currentDetElemId = -1;

    std::vector<int> de;
    std::vector<int> offset;

    for ( unsigned int i = 0; i < mManuIds.size(); ++i )
    {
        Int_t absManuId = AbsManuId(i);

        Int_t detElemId = GetDetElemIdFromAbsManuId(absManuId);

        if (detElemId != currentDetElemId)
        {
            currentDetElemId = detElemId;

            de.push_back(detElemId);
            offset.push_back(i);
        }
    }

    assert(de.size()==156);
    assert(offset.size()==156);

    std::cout << "std::vector<int> detectionElementIds = {";
    for ( std::vector<int>::size_type i = 0; i < de.size(); ++i )
    {
        std::cout << de[i];
        if ( i < de.size() - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << "};" << std::endl;

    std::cout << "std::vector<int> detectionElementIdOffsets = {";
    for ( std::vector<int>::size_type i = 0; i < offset.size(); ++i )
    {
        std::cout << offset[i];
        if ( i < offset.size() - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << "};" << std::endl;
}

std::ostream& operator<<(std::ostream& os,
        const AliMuonCompactMapping& cm)
{
    // consistency check
    for ( unsigned int i = 0; i < cm.mManuIds.size(); ++i )
    {
        Int_t absManuId = cm.AbsManuId(i);

        os << Form("ENCODEDMANUID[%6d]=%04d : DE %04d MANU %04d",
                i,absManuId,cm.GetDetElemIdFromAbsManuId(absManuId),
                cm.GetManuIdFromAbsManuId(absManuId))
            << std::endl;
    }

    return os;
}

