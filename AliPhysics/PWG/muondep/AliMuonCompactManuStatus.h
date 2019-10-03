#ifndef ALIMUONCOMPACTMANUSTATUS_H
#define ALIMUONCOMPACTMANUSTATUS_H

#include <map>
#include <vector>
#include "Rtypes.h"

/**

@ingroup pwg_muondep_compact

@class AliMuonCompactManuStatus

@brief Utility class to compute status of MCH manus

*/

class AliMuonCompactManuStatus
{
public:

    static const UInt_t MANUBADPEDMASK;
    static const UInt_t MANUBADHVMASK;
    static const UInt_t MANUBADLVMASK;
    static const UInt_t MANUBADOCCMASK;
    static const UInt_t MANUOUTOFCONFIGMASK;
    static const UInt_t MANUREJECTMASK;

    AliMuonCompactManuStatus() {}

    static std::string CauseAsString(UInt_t cause);
    std::vector<UInt_t> BuildFromOCDB(Int_t runNumber, const char* ocdbPath="raw://");
    void Print(const std::vector<UInt_t>& manuStatus, bool all = false);
    void ReadManuStatus(const char* inputfile, std::map<int,std::vector<UInt_t> >& manuStatusForRuns);
    void WriteToBinaryFile(const char* runlist, const char* outputfile, const char* ocdbPath="raw://");
};

#endif

