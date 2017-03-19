#ifndef ALIMUONCOMPACTQUICKACCEFF_H
#define ALIMUONCOMPACTQUICKACCEFF_H

#include "Rtypes.h"
#include <vector>
#include <map>

class AliMuonCompactCluster;
class AliMuonCompactEvent;
class AliMuonCompactTrack;
class TH1;
class TTree;

/**

  @ingroup pwg_muondep_compact

  @class AliMuonCompactQuickAccEff

  @brief Quick Acc x Eff processor

  This class is meant to get a quick computation of
  the evolution of the Acc x Eff for some runs.

*/


class AliMuonCompactQuickAccEff
{
    public:

        AliMuonCompactQuickAccEff(int maxevents=0, bool rejectMonoCathodeClusters=false);

        void ComputeEvolution(const std::vector<AliMuonCompactEvent>& events, 
                std::vector<int>& vrunlist,
                const std::map<int,std::vector<UInt_t> >& manuStatusForRuns,
                const char* outputfile);

        Bool_t ValidateCluster(const AliMuonCompactCluster& cl,
                const std::vector<UInt_t>& manuStatus,
                UInt_t causeMask);


        Bool_t ValidateTrack(const AliMuonCompactTrack& track,
                const std::vector<UInt_t>& manuStatus,
                UInt_t causeMask);

        TH1* ComputeMinv(const std::vector<AliMuonCompactEvent>& events,
                const std::vector<UInt_t>& manustatus,
                UInt_t causeMask,
                Int_t& npairs);

        void ComputeEvolutionFromManuStatus(const char* treeFile,
                const char* runList,
                const char* outputfile,
                const char* manustatusfile,
                const char* ocdbPath="raw://",
                Int_t runNumber=0);
        
        UInt_t GetEvents(TTree* tree,std::vector<AliMuonCompactEvent>& events, Bool_t verbose=kFALSE);

        UInt_t GetEvents(const char* treeFile, std::vector<AliMuonCompactEvent>& events, Bool_t verbose=kFALSE);

    private:
        ULong64_t fMaxEvents;
        bool fRejectMonoCathodeClusters;
};

#endif

