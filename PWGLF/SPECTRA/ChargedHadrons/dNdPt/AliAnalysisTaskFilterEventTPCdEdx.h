/// \class AliAnalysisTaskFilterEventTPCdEdx
/// \brief Task to filter out ESD events that have a very large dEdX track
///
/// TODO add proper description
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Apr 5, 2019

#ifndef AliAnalysisTaskFilterEventTPCdEdx_H
#define AliAnalysisTaskFilterEventTPCdEdx_H

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

class AliAnalysisTaskFilterEventTPCdEdx : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskFilterEventTPCdEdx();
                                AliAnalysisTaskFilterEventTPCdEdx(const char *name);
        virtual                 ~AliAnalysisTaskFilterEventTPCdEdx();

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event        
        virtual void            AnaTrack(Int_t flag = 0);        //called once for every track

        static AliAnalysisTaskFilterEventTPCdEdx* AddTaskFilterEventTPCdEdx(const char* name = "TaskFilterEventTPCdEdx", const char* treefile = 0, const char* outfile = 0);

    protected:    
        THnSparseD*             fHistdEdx;                  //-> dedx hist
        TTree*                  fesdTreeFiltered;           //! the tree with filtered esd events
        
    private:
        AliAnalysisTaskFilterEventTPCdEdx(const AliAnalysisTaskFilterEventTPCdEdx&); // not implemented
        AliAnalysisTaskFilterEventTPCdEdx& operator=(const AliAnalysisTaskFilterEventTPCdEdx&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskFilterEventTPCdEdx, 1);
    /// \endcond        
};

#endif