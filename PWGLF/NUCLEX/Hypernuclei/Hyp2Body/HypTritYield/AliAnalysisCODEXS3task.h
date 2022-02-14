// based on AliCODEX class
// contains additional trigger information (TRD and High Multiplicity Trigger)

#ifndef AliAnalysisCODEXS3task_H
#define AliAnalysisCODEXS3task_H

#include <Rtypes.h>
#include <TAxis.h>
#include <vector>

#include "AliAnalysisCODEX.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"

class AliPIDResponse;
class TList;
class TTree;
class TParticle;
class TH2I;

using AliAnalysisCODEX::Track;
using std::vector;

class AliAnalysisCODEXS3task : public AliAnalysisTaskSE {
  public:
    AliAnalysisCODEXS3task(const char *name = "AliCODEX");
    virtual ~AliAnalysisCODEXS3task();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(const Option_t*) {}
 AliInputEventHandler* handler;
    long GetParticleMask(TParticle* part);
    bool            mMCtrue;
    UInt_t          mCentralityMode;
    AliESDtrackCuts Cuts;
    AliEventCuts    mEventCuts;
    double          mPtCut;                 /// Minimum pt stored in the output trees
    double          mITSsaPtCut;            /// Maximum pt for the ITSsa track
    unsigned char   mPOI;                   /// Particles Of Interest (POI) to be stored in the output
    unsigned char   mEventPOI;              /// Events without POI with this PDG code can be rejected
    double          mNsigmaTPCselectionPOI; /// Maximum number of sigmas in the TPC from the expected signal of a POI
    double          mNsigmaTOFselectionPOI; /// Maximum number of sigmas in the TPC from the expected signal of a POI
    double          mStartingPtTOFselection;/// pt at which the TOF selection starts
    bool            mSkipEmptyEvents;       /// If true events without any tracks are not stored in the output tree
    bool            mITSstandalone;         /// If true ITS standalone PID is provided

    void Discard(const TString discard) { mToDiscard = discard;};
  	void SetTriggerMask(UInt_t triggerMask = 0) {mTriggerMask = triggerMask;};
  private:
    // Private methods
    AliAnalysisCODEXS3task(const AliAnalysisCODEXS3task&);            //! Not implemented
    AliAnalysisCODEXS3task& operator=(const AliAnalysisCODEXS3task&); //! Not implemented
    void Discard();                                               /// Method to discard a branch of the filtered tree

    TList*   mOutput; //!
    TTree*   mTree;   //!
    
TH1D * fHistEvents;
    
    AliPIDResponse* mPIDresponse; //!

    //
    AliAnalysisCODEX::Header        mHeader; /// Header
    vector<AliAnalysisCODEX::Track> mTracks; /// Tracks

    ///
    TH2I* mTimeChan;                                    //! 2D histogram with the Time/Channel correlation
    //
    TString mToDiscard; /// List of the branches to discard

		UInt_t mTriggerMask;
		
    ClassDef(AliAnalysisCODEXS3task,1)
};

#endif
