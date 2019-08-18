/// \class AliAnalysisTaskSpectraV0M
/// \brief Task for analysis of spectra vs. V0M in pp data
///
/// filling spectra histograms as well as v0 distributions
/// to be combined with a glauber to get peripheral raa
/// this task is intended for pp
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 25, 2019

#ifndef AliAnalysisTaskSpectraV0M_H
#define AliAnalysisTaskSpectraV0M_H

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

class AliAnalysisTaskSpectraV0M : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskSpectraV0M();
                                AliAnalysisTaskSpectraV0M(const char *name);
        virtual                 ~AliAnalysisTaskSpectraV0M();
        
        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event                
        virtual void            AnaTrack(Int_t flag = 0);        //called once for every track in DATA+MC event
        virtual void            AnaTrackMC(Int_t flag = 0);      //called once for every track in DATA event
        virtual void            AnaParticleMC(Int_t flag = 0);   //called once for every track in MC event        
        
        static AliAnalysisTaskSpectraV0M* AddTaskSpectraV0M(const char* name = "TaskSpectraV0M", const char* outfile = 0);

    protected:    
        THnSparseD*             fHistEffCont;         //-> efficiency/contamination histogram                
        THnSparseD*             fHistTrack;           //-> histogram of pt spectra vs. mult and cent
        THnSparseD*             fHistEvent;           //-> histogram of event numbers etc.
        
    private:
        AliAnalysisTaskSpectraV0M(const AliAnalysisTaskSpectraV0M&); // not implemented
        AliAnalysisTaskSpectraV0M& operator=(const AliAnalysisTaskSpectraV0M&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskSpectraV0M, 1);
    /// \endcond        
};

#endif