/// \class AliAnalysisTaskSpectraEtaPhi
/// \brief Simple task for quick analysis of pt spectra vs. mult and cent

#ifndef AliAnalysisTaskSpectraEtaPhi_H
#define AliAnalysisTaskSpectraEtaPhi_H

#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisHelpersHist.h"
#include "THn.h"

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

class AliAnalysisTaskSpectraEtaPhi : public AliAnalysisTaskMKBase {
  public:
    AliAnalysisTaskSpectraEtaPhi();
    AliAnalysisTaskSpectraEtaPhi(const char* name);
    virtual ~AliAnalysisTaskSpectraEtaPhi();

    virtual void AddOutput();         // called at the beginning
    virtual Bool_t IsEventSelected(); // called for each event
    virtual void AnaEvent();          // called once for every selected event
    virtual void
    AnaTrack(Int_t flag = 0); // called once for every track in DATA+MC event
    virtual void
    AnaTrackMC(Int_t flag = 0); // called once for every track in DATA event
    virtual void
    AnaParticleMC(Int_t flag = 0); // called once for every track in MC event

    static AliAnalysisTaskSpectraEtaPhi*
    AddTaskSpectra(const char* name = "TaskSpectraEtaPhi",
                   const char* outfile = 0);

  protected:
    Hist::Hist<THnF> fHistEffContNCluster; //!<!   efficiency/contamination histogram
    Hist::Hist<THnF> fHistEffContZ;        //!<!   efficiency/contamination histogram
    Hist::Hist<THnF> fHistEffContEta;      //!<!   efficiency/contamination histogram
    Hist::Hist<THnF> fHistEffContPhi;      //!<!   efficiency/contamination histogram
    Hist::Hist<THnF> fHistTrackNCluster;   //!<!   histogram of pt spectra vs. mult and cent
    Hist::Hist<THnF> fHistTrackZ;          //!<!   histogram of pt spectra vs. mult and cent
    Hist::Hist<THnF> fHistTrackEta;        //!<!   histogram of pt spectra vs. mult and cent
    Hist::Hist<THnF> fHistTrackPhi;        //!<!   histogram of pt spectra vs. mult and cent
    Hist::Hist<TH3F> fHistEvent;           //!<!   histogram of event numbers etc.

  private:
    AliAnalysisTaskSpectraEtaPhi(
        const AliAnalysisTaskSpectraEtaPhi&); // not implemented
    AliAnalysisTaskSpectraEtaPhi&
    operator=(const AliAnalysisTaskSpectraEtaPhi&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSpectraEtaPhi, 1);
    /// \endcond
};

#endif
