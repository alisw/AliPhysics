/// \class AliAnalysisTaskEffContStudy
/// \brief Task to study efficiency and contamination in MC
///
/// Fills THnSparse for corrections on tracking efficiency and 
/// contamination from secondary particles
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 8, 2019

#ifndef AliAnalysisTaskEffContStudy_H
#define AliAnalysisTaskEffContStudy_H

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

class AliAnalysisTaskEffContStudy : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskEffContStudy();
                                AliAnalysisTaskEffContStudy(const char *name);
        virtual                 ~AliAnalysisTaskEffContStudy();

        virtual void            AddOutput(); //called at the beginning
        virtual void            AnaTrack();  //called once for every track
        virtual void            AnaMCParticle();  //called once for every mc particle
        virtual void            AnaEvent();  //called once for every event        
        
        static AliAnalysisTaskEffContStudy* AddTaskEffContStudy(const char* name = "TaskEffContStudy", const char* outfile = 0);

    protected:    
        THnSparseD*             fHistEffCont;         //-> efficiency/contamination histogram                
        THnSparseD*             fHistEffContScaled;   //-> efficiency/contamination histogram including scaling, for testing
        
    private:
        AliAnalysisTaskEffContStudy(const AliAnalysisTaskEffContStudy&); // not implemented
        AliAnalysisTaskEffContStudy& operator=(const AliAnalysisTaskEffContStudy&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskEffContStudy, 1);
    /// \endcond        
};

#endif