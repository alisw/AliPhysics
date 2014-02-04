#ifndef ALIANALYSISMUMUBASE_H
#define ALIANALYSISMUMUBASE_H

/**
 *
 * \class AliAnalysisMuMuBase
 *
 * \brief Base class of the sub-analysis for AliAnalysisTaskMuMu
 *
 * \author L. Aphecetche (Subatech)
 */

#include "TObject.h"
#include "TString.h"

class AliCounterCollection;
class AliAnalysisMuMuBinning;
class AliMergeableCollection;
class AliVParticle;
class AliVEvent;
class AliMCEvent;
class TH1;
class AliInputEventHandler;
class AliAnalysisMuMuCutRegistry;

class AliAnalysisMuMuBase : public TObject
{
public:

  AliAnalysisMuMuBase();
  virtual ~AliAnalysisMuMuBase() {}

  /** Define the histograms needed for the path starting at eventSelection/triggerClassName/centrality.
   * This method has to ensure the histogram creation is performed only once !
   */
  virtual void DefineHistogramCollection(const char* eventSelection,
                                         const char* triggerClassName,
                                         const char* centrality) = 0;

  /** Fill histograms for one event */
  virtual void FillHistosForEvent(const char* /*eventSelection*/,const char* /*triggerClassName*/,const char* /*centrality*/) {}

  /** Fill histograms for one MC event */
  virtual void FillHistosForMCEvent(const char* /*eventSelection*/,const char* /*triggerClassName*/,const char* /*centrality*/) {}
  
  /** Fill histograms for one track */
  virtual void FillHistosForTrack(const char* /*eventSelection*/,const char* /*triggerClassName*/,const char* /*centrality*/,
                                  const char* /*trackCutName*/,
                                  const AliVParticle& /*part*/) {}

  /** Fill histograms for one track pair */
  virtual void FillHistosForPair(const char* /*eventSelection*/,const char* /*triggerClassName*/,const char* /*centrality*/,
                                 const char* /*pairCutName*/,
                                 const AliVParticle& /*part1*/,
                                 const AliVParticle& /*part2*/) {}
  
  virtual void Init(AliCounterCollection& cc,
                    AliMergeableCollection& hc,
                    const AliAnalysisMuMuBinning& binning,
                    const AliAnalysisMuMuCutRegistry& cutRegister);

  virtual void SetEvent(AliVEvent* event, AliMCEvent* mcEvent=0x0);

  virtual Bool_t IsHistogramDisabled(const char* hname) const;
  
  virtual Bool_t IsHistogrammingDisabled() const;
  
  virtual void DisableHistograms(const char* pattern="*");
  
  AliVEvent* Event() const { return fEvent; }
  
  AliMCEvent* MCEvent() const { return fMCEvent; }
  
  /// Called at each new run
  void SetRun(const AliInputEventHandler* /*eventHandler*/) {}

  void Terminate(Option_t* /*opt*/="") {}
  
  enum EDataType
  {
    kHistoForMCInput = (1<<0),
    kHistoForData = (1<<1)
  };
  
  void SetMC() { fHasMC = kTRUE; }
  
  Bool_t HasMC() const { return fHasMC; }
  
  Bool_t AlwaysTrue(const AliVEvent& /*event*/) const { return kTRUE; }
  Bool_t AlwaysTrue(const AliVParticle& /*particle*/) const { return kTRUE; }
  Bool_t AlwaysTrue(const AliVParticle& /*particle*/, const AliVParticle& /*particle*/) const { return kTRUE; }
  void NameOfAlwaysTrue(TString& name) const { name = "ALL"; }

  Bool_t AlwaysFalse(const AliVEvent& /*event*/) const { return kFALSE; }
  Bool_t AlwaysFalse(const AliVParticle& /*particle*/) const { return kFALSE; }
  Bool_t AlwaysFalse(const AliVParticle& /*particle*/, const AliVParticle& /*particle*/) const { return kFALSE; }
  void NameOfAlwaysFalse(TString& name) const { name = "NONE"; }
  
protected:

  void CreateHistos(const TObjArray& paths,
                    const char* hname, const char* htitle,
                    Int_t nbinsx, Double_t xmin, Double_t xmax,
                    Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0) const;

  void CreateEventHistos(UInt_t dataType,
                         const char* what,
                         const char* hname, const char* htitle,
                         Int_t nbinsx, Double_t xmin, Double_t xmax,
                         Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0) const;

  void CreateEventHistos(UInt_t dataType,
                         const char* eventSelection,
                         const char* triggerClassName,
                         const char* centrality,
                         const char* hname, const char* htitle,
                         Int_t nbinsx, Double_t xmin, Double_t xmax,
                         Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0) const;
  
  
  void CreateTrackHistos(const char* eventSelection,
                         const char* triggerClassName,
                         const char* centrality,
                         const char* hname, const char* htitle,
                         Int_t nbinsx, Double_t xmin, Double_t xmax,
                         Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0) const;
  
  void CreatePairHistos(const char* eventSelection,
                        const char* triggerClassName,
                        const char* centrality,
                        const char* hname, const char* htitle,
                        Int_t nbinsx, Double_t xmin, Double_t xmax,
                        Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0) const;

  TH1* Histo(const char* eventSelection, const char* histoname);
  TH1* Histo(const char* eventSelection, const char* triggerClassName, const char* histoname);
  TH1* Histo(const char* eventSelection, const char* triggerClassName, const char* cent, const char* histoname);
  TH1* Histo(const char* eventSelection, const char* triggerClassName, const char* cent,
             const char* what, const char* histoname);

  TH1* MCHisto(const char* eventSelection, const char* histoname);
  TH1* MCHisto(const char* eventSelection, const char* triggerClassName, const char* histoname);
  TH1* MCHisto(const char* eventSelection, const char* triggerClassName, const char* cent, const char* histoname);
  TH1* MCHisto(const char* eventSelection, const char* triggerClassName, const char* cent,
             const char* what, const char* histoname);

  Int_t GetNbins(Double_t xmin, Double_t xmax, Double_t xstep);

  AliCounterCollection* CounterCollection() const { return fEventCounters; }
  AliMergeableCollection* HistogramCollection() const { return fHistogramCollection; }
  const AliAnalysisMuMuBinning* Binning() const { return fBinning; }
  const AliAnalysisMuMuCutRegistry* CutRegistry() const { return fCutRegistry; }
  
private:
  
  AliCounterCollection* fEventCounters; //! event counters
  AliMergeableCollection* fHistogramCollection; //! collection of histograms
  const AliAnalysisMuMuBinning* fBinning; //! binning for particles
  const AliAnalysisMuMuCutRegistry* fCutRegistry; //! registry of cut combinations
  AliVEvent* fEvent; //! current event
  AliMCEvent* fMCEvent; //! current MC event
  TList* fHistogramToDisable; // list of regexp of histo name to disable
  Bool_t fHasMC; // whether or not we're dealing with MC data
  
  ClassDef(AliAnalysisMuMuBase,1) // base class for a companion class to AliAnalysisMuMu
};

#endif

