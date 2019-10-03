#ifndef ALIMUONCOMPACTMAPPING_H
#define ALIMUONCOMPACTMAPPING_H

#include <vector>
#include <map>
#include <iostream>

#include "Rtypes.h"

/**

@ingroup pwg_muondep_compact

@struct AliMuonCompactMapping

@brief A minimal mapping at Manu level

The base concept here is to get an absolute index 0..16827 for each manu.

*/

struct AliMuonCompactMapping
{
    AliMuonCompactMapping() : mManuIds(), mManuMap(), mNpads() {}

    // array containing the 32bits encoded
    // pair (detElemId,manuId) for each
    // absolute manu index
    // mManuIds[abs]=(de << 16) & local
    std::vector<UInt_t> mManuIds;

    // associate each encoded pair (detElemId,manuId)
    // to an index in mManuIds
    // (i.e. reverse structure of mManuIds
    std::map<int,int> mManuMap;

    // number of pads per manu
    std::vector<int> mNpads;
    
    friend std::ostream& operator<<(std::ostream& os, const AliMuonCompactMapping& cm);

    Int_t GetDetElemIdFromAbsManuIndex(Int_t index) const;

    Int_t AbsManuId(UInt_t index) const;

    Int_t GetManuIdFromAbsManuId(UInt_t absManuId) const;

    Int_t GetDetElemIdFromAbsManuId(UInt_t absManuId) const;

    Int_t GetNofPadsFromAbsManuIndex(Int_t index) const;

    Int_t FindManuAbsIndex(Int_t detElemId, Int_t manuId) const;

    static AliMuonCompactMapping* GetCompactMapping(const char* ocdbPath="raw://", Int_t runNumber=0);

    void GenerateStaticOffsets();
};

#endif
