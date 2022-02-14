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
        
        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event        
        virtual void            AnaTrackMC(Int_t flag = 0);    //called once for every track in DATA event
        virtual void            AnaParticleMC(Int_t flag = 0);      //called once for every track in MC event        
        
        static AliAnalysisTaskEffContStudy* AddTaskEffContStudy(const char* name = "TaskEffContStudy", const char* outfile = 0);

    protected:    
        THnSparseF*             fHistEffCont;         //-> efficiency/contamination histogram
        THnSparseF*             fHistEffContScaled;   //-> efficiency/contamination histogram including scaling, for testing
        
    private:
        AliAnalysisTaskEffContStudy(const AliAnalysisTaskEffContStudy&); // not implemented
        AliAnalysisTaskEffContStudy& operator=(const AliAnalysisTaskEffContStudy&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskEffContStudy, 1);
    /// \endcond        
};

#endif
