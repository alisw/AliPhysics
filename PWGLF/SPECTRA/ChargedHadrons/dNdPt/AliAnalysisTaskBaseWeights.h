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

        // new functions added to steer the behavior of the MC re-weighting
        virtual void            SetUseMCWeights(Bool_t use=kTRUE) { fUseMCWeights = use; }              
        virtual Bool_t          GetUseMCWeights() { return fUseMCWeights;  }
        virtual void            SetUseRandomSeed(Bool_t use=kTRUE) { fUseRandomSeed = use; }              
        virtual Bool_t          GetUseRandomSeed() { return fUseRandomSeed;  }        
             
        static AliAnalysisTaskBaseWeights* AddTaskBaseWeights(const char* name = "TaskBaseWeights", const char* outfile = 0, const char* collisionSystem = 0, Int_t sysFlag = 0, const char* prevTrainOutputPath = 0);

    protected:    
        // override the track and particle loops from the AliAnalysisTaskMKBase class
        virtual void            LoopOverAllTracks(Int_t flag = 0);    // loops over all tracks in the event, calls AnaTrack(), AnaTrackMC() and AnaTrackDATA() for each track
        virtual void            LoopOverAllParticles(Int_t flag = 0); // loops over all MC particles in the event, calls AnaParticleMC() for each particle
        virtual void            FillDefaultHistograms(Int_t step); // overwrite this function to include the filling of the MCSpectraWeights histograms
        virtual Double_t        MCScalingFactor();  //internal method to determine the proper scaling factor for primaries and secondaries
        virtual void            BaseAddOutput();   // add even more output
        
        virtual UInt_t          GetSeed();
        
        Bool_t               fUseMCWeights;      ///  use mc weights (default) or do nothing
        Bool_t               fUseRandomSeed;     ///  use a random seed or a deterministic one (default)        
        TRandom*             fRand;              //-> random generator to be used
        AliMCSpectraWeights* fMCSpectraWeights;  //-> object to determine efficiency scaling        
        Double_t             fMCweight;          //!<! MC weight of the current track/particle
        
    private:
        AliAnalysisTaskBaseWeights(const AliAnalysisTaskBaseWeights&); // not implemented
        AliAnalysisTaskBaseWeights& operator=(const AliAnalysisTaskBaseWeights&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskBaseWeights, 3);
    /// \endcond        
};

#endif