#ifndef AliAnalysisCODEXtask_H
#define AliAnalysisCODEXtask_H

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

class AliAnalysisCODEXtask : public AliAnalysisTaskSE {
  public:
    AliAnalysisCODEXtask(const char *name = "AliCODEX");
    virtual ~AliAnalysisCODEXtask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(const Option_t*) {}

    long GetParticleMask(TParticle* part);

    /// You know why I am putting these members as public?
    /// Because I am sick of writing down setter and getter
    /// for classes that I, and only I, use.
    /// If you don't like it, write 3x lines of code for each
    /// of these members...
    bool            mMCtrue;
    UInt_t          mCentralityMode;
    AliESDtrackCuts Cuts;
    AliEventCuts    mEventCuts;
    double          mPtCut;                 /// Minimum pt stored in the output trees
    unsigned char   mPOI;                   /// Particles Of Interest (POI) to be stored in the output
    double          mNsigmaTPCselectionPOI; /// Maximum number of sigmas in the TPC from the expected signal of a POI
  private:
    // Private methods
    AliAnalysisCODEXtask(const AliAnalysisCODEXtask&);            //! Not implemented
    AliAnalysisCODEXtask& operator=(const AliAnalysisCODEXtask&); //! Not implemented

    TList*   mOutput; //!
    TTree*   mTree;   //!

    AliPIDResponse* mPIDresponse; //!

    //
    AliAnalysisCODEX::Header        mHeader; /// Header
    vector<AliAnalysisCODEX::Track> mTracks; /// Tracks

    ///
    TH2I* mTimeChan;                                    /// 2D histogram with the Time/Channel correlation

    ClassDef(AliAnalysisCODEXtask,1)
};

#endif
