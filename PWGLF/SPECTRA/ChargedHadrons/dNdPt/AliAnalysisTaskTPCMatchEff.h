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

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event        
        virtual void            AnaTrackDATA(Int_t flag = 0);    //called once for every track in data
        virtual void            AnaTrackMC(Int_t flag = 0);      //called once for every track in mc        
        
        static AliAnalysisTaskTPCMatchEff* AddTaskTPCMatchEff(const char* name = "TaskTPCMatchEff", const char* outfile = 0);

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