/// \class AliAnalysisTaskSpectraINEL0
/// \brief Simple task for study INEL0 normalisation
///
/// \author Patrick Huhn <patrick.huhn@cern.ch>, CERN
/// \date Mai 08, 2020

#ifndef AliAnalysisTaskSpectraINEL0_H
#define AliAnalysisTaskSpectraINEL0_H

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

class AliAnalysisTaskSpectraINEL0 : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskSpectraINEL0();
                                AliAnalysisTaskSpectraINEL0(const char *name);
        virtual                 ~AliAnalysisTaskSpectraINEL0();

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual Bool_t          oldEventCuts();               //it tests, whether old event cuts are fullfilled
        virtual void            FillDefaultHistograms(Int_t step); // called twice for each event
        virtual void            AnaTrack(Int_t flag = 0);        //called once for every track in DATA+MC event
        virtual void            AnaTrackMC(Int_t flag = 0);      //called once for every track in DATA event
        virtual void            AnaParticleMC(Int_t flag = 0);   //called once for every track in MC event

        static AliAnalysisTaskSpectraINEL0* AddTaskSpectraINEL0(const char* name = "TaskSpectra", const char* outfile = 0);

    protected:
        THnSparseF*             fHistEffCont;         //-> efficiency/contamination histogram
        THnSparseF*             fHistTrack;           //-> histogram of pt spectra vs. mult and cent
        THnSparseF*             fHistEvent;           //-> histogram of event numbers etc.
        THnSparseF*             fHistTrackINEL0;      //-> histogram for INEL0 study
        THnSparseF*             fHistVtxInfo;      //-> vertex reconstruction efficiency histogram


    private:
        AliAnalysisTaskSpectraINEL0(const AliAnalysisTaskSpectraINEL0&); // not implemented
        AliAnalysisTaskSpectraINEL0& operator=(const AliAnalysisTaskSpectraINEL0&); // not implemented

    /// \cond CLASSIMP
        ClassDef(AliAnalysisTaskSpectraINEL0, 1);
    /// \endcond
};

#endif /* AliAnalysisTaskSpectraINEL0_hpp */
