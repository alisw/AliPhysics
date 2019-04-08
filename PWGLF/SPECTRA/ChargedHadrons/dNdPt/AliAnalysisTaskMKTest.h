/// \class AliAnalysisTaskMKTest
/// \brief Example of AnalysisTask derived from AliAnalysisTaskMKBase
///
/// Fills THnSparse with pt of all tracks that pass esd track cuts
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 8, 2019

#ifndef AliAnalysisTaskMKTest_H
#define AliAnalysisTaskMKTest_H

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

class AliAnalysisTaskMKTest : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskMKTest();
                                AliAnalysisTaskMKTest(const char *name);
        virtual                 ~AliAnalysisTaskMKTest();

        virtual void            AddOutput(); //called at the beginning
        virtual void            AnaTrack();  //called once for every track
        virtual void            AnaEvent();  //called once for every event        
        
        static AliAnalysisTaskMKTest* AddTaskMKTest(const char* name = "TaskMKTest", const char* outfile = 0);

    protected:    
        THnSparseD*             fHistPt;     //-> pt hist  
        
    private:
        AliAnalysisTaskMKTest(const AliAnalysisTaskMKTest&); // not implemented
        AliAnalysisTaskMKTest& operator=(const AliAnalysisTaskMKTest&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskMKTest, 1);
    /// \endcond        
};

#endif