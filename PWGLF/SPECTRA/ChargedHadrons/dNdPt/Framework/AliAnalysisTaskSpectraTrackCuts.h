/// \class AliAnalysisTaskTrackCuts
/// \brief AnalysisTask to study track cuts
///
/// Fills Histograms of AliESDTrackCuts class
///
/// \author Patrick Huhn <patrick.huhn@cern.ch>
/// \date Feb. 21, 2020

#ifndef AliAnalysisTaskTrackCuts_H
#define AliAnalysisTaskTrackCuts_H

#include "AliAnalysisTaskMKBase.h"

class AliESDtrackCuts;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliMCParticle;

class AliAnalysisTaskSpectraTrackCuts : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskSpectraTrackCuts();
                                AliAnalysisTaskSpectraTrackCuts(const char *name);
        virtual                 ~AliAnalysisTaskSpectraTrackCuts();

        virtual void           Terminate(Option_t* option);

        static AliAnalysisTaskSpectraTrackCuts* AddTaskTrackCuts(const char* name = "TrackCuts", const char* outfile = 0, int _cutMode=100);

        
    private:
        AliAnalysisTaskSpectraTrackCuts(const AliAnalysisTaskSpectraTrackCuts&); // not implemented
        AliAnalysisTaskSpectraTrackCuts& operator=(const AliAnalysisTaskSpectraTrackCuts&); // not implemented
        
    /// \cond CLASSIMP
        ClassDef(AliAnalysisTaskSpectraTrackCuts, 1);
    /// \endcond
};

#endif
