/// \class AliAnalysisTaskSpectra
/// \brief Simple task for quick analysis of pt spectra vs. mult and cent
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 25, 2019

#ifndef AliAnalysisTaskSpectra_H
#define AliAnalysisTaskSpectra_H

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

class AliAnalysisTaskSpectra : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskSpectra();
                                AliAnalysisTaskSpectra(const char *name);
        virtual                 ~AliAnalysisTaskSpectra();
        
        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event                
        virtual void            AnaTrack(Int_t flag = 0);        //called once for every track in DATA+MC event
        virtual void            AnaTrackMC(Int_t flag = 0);      //called once for every track in DATA event
        virtual void            AnaParticleMC(Int_t flag = 0);   //called once for every track in MC event        
        
        static AliAnalysisTaskSpectra* AddTaskSpectra(const char* name = "TaskSpectra", const char* outfile = 0);

    protected:    
        THnSparseF*             fHistEffCont;         //-> efficiency/contamination histogram
        THnSparseF*             fHistTrack;           //-> histogram of pt spectra vs. mult and cent
        THnSparseF*             fHistEvent;           //-> histogram of event numbers etc.
        
    private:
        AliAnalysisTaskSpectra(const AliAnalysisTaskSpectra&); // not implemented
        AliAnalysisTaskSpectra& operator=(const AliAnalysisTaskSpectra&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskSpectra, 1);
    /// \endcond        
};

#endif
