/// \class AliAnalysisTaskTPCMatchEff
/// \brief task to study the tpc-its matching efficiency with different cuts
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 8, 2019

#ifndef AliAnalysisTaskTPCMatchEff_H
#define AliAnalysisTaskTPCMatchEff_H

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

class AliAnalysisTaskTPCMatchEff : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskTPCMatchEff();
                                AliAnalysisTaskTPCMatchEff(const char *name);
        virtual                 ~AliAnalysisTaskTPCMatchEff();

        virtual void            AddOutput(); //called at the beginning
        virtual void            AnaTrack();  //called once for every track        
        virtual void            AnaEvent();  //called once for every event        
        
        static AliAnalysisTaskTPCMatchEff* AddTaskTPCMatchEff(const char* name = "TaskTPCMatchEff");

    protected:    
        THnSparseD*             fHistMCMatchEff;     //-> for mc matching efficiency
        THnSparseD*             fHistDATAMatchEff;   //-> for data matching efficiency
        
    private:
        AliAnalysisTaskTPCMatchEff(const AliAnalysisTaskTPCMatchEff&); // not implemented
        AliAnalysisTaskTPCMatchEff& operator=(const AliAnalysisTaskTPCMatchEff&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskTPCMatchEff, 1);
    /// \endcond        
};

#endif