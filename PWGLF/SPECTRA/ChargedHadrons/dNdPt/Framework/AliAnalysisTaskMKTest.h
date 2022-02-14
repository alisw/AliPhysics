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

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event        
        virtual void            AnaTrack(Int_t flag = 0);        //called once for every track

        static AliAnalysisTaskMKTest* AddTaskMKTest(const char* name = "TaskMKTest", const char* outfile = 0);

    protected:    
        THnSparseF*             fHistPt;     //-> pt hist  
        
    private:
        AliAnalysisTaskMKTest(const AliAnalysisTaskMKTest&); // not implemented
        AliAnalysisTaskMKTest& operator=(const AliAnalysisTaskMKTest&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskMKTest, 1);
    /// \endcond        
};

#endif
