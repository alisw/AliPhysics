/// \class AliAnalysisTaskBaseWeights
/// \brief Task implementing MC weights for secondaries and contamination
///
/// This Tasks overrides the track and particle loops from AliAnalysisTaskMKBase
/// to call the Ana functions multiple times or less per mc particle
/// 
/// actual tasks doing analysis have to derive from this task in the same 
/// way as from AliAnalysisTaskMKBase
///
/// On data the task has no effect
/// On MC everything is transparent
/// 
/// by default, weigths are active
/// call SetUseMCWeights(kFALSE) to deactivate
/// there is randomness in the task, but by default it is fully deterministic
/// to switch on real randromness call SetUseRandomSeed(kTRUE)
/// otherwise seeds are based on the event number
/// 
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 24, 2019

#ifndef AliAnalysisTaskBaseWeights_H
#define AliAnalysisTaskBaseWeights_H

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
class TRandom;
class AliMCSpectraWeights;

class AliAnalysisTaskBaseWeights : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskBaseWeights();
                                AliAnalysisTaskBaseWeights(const char *name);
        virtual                 ~AliAnalysisTaskBaseWeights();

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event
        virtual void            LoopOverAllTracks(Int_t flag = 0);
        virtual void            LoopOverAllParticles(Int_t flag = 0);
        virtual void            AnaTrackMC(Int_t flag = 0);      //called once for every track in DATA event
        virtual void            AnaParticleMC(Int_t flag = 0);   //called once for every track in MC event
    
        // new functions added to steer the behavior of the MC re-weighting
        virtual void            SetUseMCWeights(Bool_t use=kTRUE) { fUseMCWeights = use;}   
        virtual Bool_t          GetUseMCWeights() { return fUseMCWeights;  }
    
        virtual void            SetUseRandomSeed(Bool_t use=kTRUE) { fUseRandomSeed = use; }
        virtual Bool_t          GetUseRandomSeed() { return fUseRandomSeed;  }        
    
        double                  GetRandomRoundDouble(double val);
        Double_t                MCScalingFactor();
        virtual void            FillDefaultHistograms(Int_t step=0);
    
        static AliAnalysisTaskBaseWeights* AddTaskBaseWeights(const char* name = "TaskBaseWeights", const char* outfile = 0, const char* collisionSystem = 0, Int_t sysFlag = 0, const char* prevTrainOutputPath = 0);
    protected:
        virtual UInt_t          GetSeed();
        
        Bool_t                  fUseMCWeights;      ///  use mc weights (default) or do nothing
        Bool_t                  fUseRandomSeed;     ///  use a random seed or a deterministic one (default)
        TRandom*                fRand;              //-> random generator to be used
        AliMCSpectraWeights*    fMCSpectraWeights;  //-> object to determine efficiency scaling
        Double_t                fMCweight;          //!<! MC weight of the current track/particle
        Double_t                fMCweightRandom;          //!<! MC weight of the current track/particle random rounded
        Double_t                fNch;
        Double_t                fNchWeighted;
        Double_t                fNchWeightedRandom;
        Double_t                fNacc;
        Double_t                fNaccWeighted;
        Double_t                fNaccWeightedRandom;
        THnSparseF*             fHistEffCont;         //-> efficiency/contamination histogram pure/weighted/randomWeight
        THnSparseF*             fHistMultCorrelation; //!<! N_acc vs N_ch pure/weighted/randomWeight
    private:
        AliAnalysisTaskBaseWeights(const AliAnalysisTaskBaseWeights&); // not implemented
        AliAnalysisTaskBaseWeights& operator=(const AliAnalysisTaskBaseWeights&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskBaseWeights, 4);
    /// \endcond        
};

#endif
