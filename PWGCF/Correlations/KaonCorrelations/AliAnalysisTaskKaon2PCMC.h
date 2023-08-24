/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskKaon2PCMC:
// Description: Analysis task to calculate Two-Particle Angular Correlation Functions of Neutral and Charged Kaons
// Author: Anjaly Sasikumar Menon
// (anjaly.sasikumar.menon@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskKaon2PCMC_H
#define AliAnalysisTaskKaon2PCMC_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "TMacro.h"

class AliAODEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliMCParticle;
class THnSparse;
class AliAODv0;
class AliAnalysisTaskKaon2PCMC : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskKaon2PCMC();
                                AliAnalysisTaskKaon2PCMC(const char *name);
        virtual                 ~AliAnalysisTaskKaon2PCMC();

        virtual void            UserCreateOutputObjects();
        Bool_t SelectK0TracksMC(AliMCParticle *mcTrack);
        Bool_t SelectKPosTracksMC(AliMCParticle *mcTrack);
        Bool_t SelectKNegTracksMC(AliMCParticle *mcTrack);
        Bool_t SelectKchTracksMC(AliMCParticle  *mcTrack);
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            Fill2DHistMCTruth(Double_t DPhi, Double_t DEta, TH2F* hist);
        virtual void            Fill2DHist(Double_t DPhi, Double_t DEta, TH3F* hist, Double_t fWeight);
        //virtual void ProcessMCParticles();

    private:
       AliMCEvent*             fmcEvent;       //! input event
       TList*                  fOutputList;    //! output list
       AliPIDResponse*         fPIDResponse;   //! pid response object’

       AliMCEvent*             fMCEvent;
       THnSparse*              fMCK0;
       THnSparse*              fMCKpos;
       THnSparse*              fMCKneg;
       TH2F*                   fHistKpKnMC;
       TH2F*                   fHistK0KchMC;
       TH1D*                   fHistGenMultiplicity;

       //Bool_t                        fAnalysisMC; // enable MC study
       //Bool_t                        fRejectEventPileUp; // enable to use Pile-up cuts
       
       Double_t        fPV[3];


       AliAnalysisTaskKaon2PCMC(const AliAnalysisTaskKaon2PCMC&); // not implemented
       AliAnalysisTaskKaon2PCMC& operator=(const AliAnalysisTaskKaon2PCMC&); // not implemented

       ClassDef(AliAnalysisTaskKaon2PCMC, 1);
};

#endif
