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
class AliMCSpectraWeights;

class AliAnalysisTaskSpectraEtaPhi : public AliAnalysisTaskMKBase {
  public:
    AliAnalysisTaskSpectraEtaPhi();
    AliAnalysisTaskSpectraEtaPhi(const char* name);
    virtual ~AliAnalysisTaskSpectraEtaPhi();

    virtual void AddOutput();         // called at the beginning
    virtual Bool_t IsEventSelected(); // called for each event
    virtual void AnaEvent();          // called once for every selected event
    virtual void AnaEventMC();          // called once for every selected MC event
    virtual void AnaTrackMC(Int_t flag = 0);      //called once for every track in MC event
    virtual void
    AnaTrack(Int_t flag = 0); // called once for every track in DATA+MC event
    virtual void
    AnaParticleMC(Int_t flag = 0); // called once for every track in MC event

    static AliAnalysisTaskSpectraEtaPhi*
    AddTaskSpectra(const char* name = "TaskSpectraEtaPhi",
                   const char* outfile = 0);

  protected:
    Hist<THnF> fHistEffContNCluster; //!<!   efficiency/contamination histogram
    Hist<THnF> fHistEffContZ;        //!<!   efficiency/contamination histogram
    Hist<THnF> fHistEffContEta;      //!<!   efficiency/contamination histogram
    Hist<THnF> fHistEffContPhi;      //!<!   efficiency/contamination histogram
    Hist<THnF> fHistTrackNCluster;   //!<!   histogram of pt spectra vs. mult and cent
    Hist<THnF> fHistTrackZ;          //!<!   histogram of pt spectra vs. mult and cent
    Hist<THnF> fHistTrackEta;        //!<!   histogram of pt spectra vs. mult and cent
    Hist<THnF> fHistTrackPhi;        //!<!   histogram of pt spectra vs. mult and cent
    Hist<TH3F> fHistEvent;           //!<!   histogram of event numbers etc.
    AliMCSpectraWeights*    fMCSpectraWeights; //!<! object to determine efficiency scaling

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
