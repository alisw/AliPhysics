/// \class AliAnalysisTaskSpectraTrackCuts
/// \brief Simple task for quick look at track cuts
///
/// \author Patrick Huhn <patrick.huhn@cern.ch>, CERN
/// \date Mar 25, 2019

#ifndef AliAnalysisTaskSpectraTrackCuts_H
#define AliAnalysisTaskSpectraTrackCuts_H

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

class AliAnalysisTaskSpectraTrackCuts : public AliAnalysisTaskMKBase {
  public:
    AliAnalysisTaskSpectraTrackCuts();
    AliAnalysisTaskSpectraTrackCuts(const char* name);
    virtual ~AliAnalysisTaskSpectraTrackCuts();

    virtual void Terminate(Option_t *option="");

    virtual Bool_t IsEventSelected(); // called for each event
    virtual void AnaEvent();          // called once for every selected event
    virtual void
    AnaTrack(Int_t flag = 0); // called once for every track in DATA+MC event
    virtual void
    AnaTrackMC(Int_t flag = 0); // called once for every track in DATA event
    virtual void
    AnaParticleMC(Int_t flag = 0); // called once for every track in MC event

    static AliAnalysisTaskSpectraTrackCuts*
    AddTaskSpectraTrackCuts(const char* name = "TaskSpectraEtaPhi",
                            const char* outfile = 0, int _CutMode = 100);

  private:
    AliAnalysisTaskSpectraTrackCuts(
        const AliAnalysisTaskSpectraTrackCuts&); // not implemented
    AliAnalysisTaskSpectraTrackCuts&
    operator=(const AliAnalysisTaskSpectraTrackCuts&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSpectraTrackCuts, 1);
    /// \endcond
};

#endif
