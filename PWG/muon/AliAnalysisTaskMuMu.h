#ifndef ALIANALYSISTASKMUMU_H
#define ALIANALYSISTASKMUMU_H

/**
 * \defgroup pwg-muon-mumu pwg-muon-mumu
 *
 * \brief Small sub-framework to analyse muon pairs and more...
 *
 * Started as a simple invariant mass analysis and grew into a bit more general thing...
 *
 * Can now compute the charged particle multiplicy (from SPD tracklets only) in order
 * to be able to correlate it with e.g. J/psi or single mu.
 */

/**
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisTaskMuMu 
 *
 * \brief Steering class for mu pairs analysis (and more...)
 *
 * This class acts as a small sub-framework to steer various sub-analysis which 
 * share the same MergeableCollection and the same CounterCollection.
 *
 *  \author: L. Aphecetche (Subatech)
 */

#ifndef ALIANALYSISTASKSE_H
#  include "AliAnalysisTaskSE.h"
#endif

#ifndef ROOT_TMath
#  include "TMath.h"
#endif

class AliAnalysisMuMuBinning;
class AliCounterCollection;
class AliMergeableCollection;
class AliVParticle;
class TList;
class TObjArray;
class AliAnalysisMuMuBase;
class AliAnalysisMuMuCutRegistry;

class AliAnalysisTaskMuMu : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskMuMu();
  virtual ~AliAnalysisTaskMuMu();

  AliAnalysisMuMuCutRegistry* CutRegistry() const;
  
  AliAnalysisMuMuBinning* Binning() const;

  void AdoptSubAnalysis(AliAnalysisMuMuBase* analysis);
  
  virtual void DisableHistograms(const char* pattern="*");

  void SetBeamYear(const char* beamYear) { fBeamYear = beamYear; }
  
  virtual void FinishTaskOutput();
  
  virtual void NotifyRun();
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void Terminate(Option_t *);
  
  void UserCreateOutputObjects();

  virtual void UserExec(Option_t* opt);
  
private:
  
  void CreateTrackHisto(const char* eventSelection,
                        const char* triggerClassName,
                        const char* hname, const char* htitle,
                        Int_t nbinsx, Double_t xmin, Double_t xmax,
                        Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0,
                        Bool_t separatePlusAndMinus=kFALSE) const;
  
  void CreatePairHisto(const char* eventSelection,
                       const char* triggerClassName,
                       const char* hname, const char* htitle,
                       Int_t nbinsx, Double_t xmin, Double_t xmax,
                       Int_t nbinsy=-1, Double_t ymin=0.0, Double_t ymax=0.0) const;
  
  const char* DefaultCentralityName() const;

  AliVEvent* Event() const;
  
  void FillHistos(const char* eventSelection, const char* triggerClassName, const char* centrality);
  
  void Fill(const char* eventSelection, const char* triggerClassName);
  
  void FillMC();
  
  void GetSelectedTrigClassesInEvent(const AliVEvent* event, TObjArray& array);

  Bool_t IsHistogrammingDisabled() const;
  
  virtual Bool_t IsHistogramDisabled(const char* hname) const;
  
  Bool_t IsPP() const;
  
private:
  
  AliAnalysisTaskMuMu(const AliAnalysisTaskMuMu&); // not implemented (on purpose)
  AliAnalysisTaskMuMu& operator=(const AliAnalysisTaskMuMu&); // not implemented (on purpose)

private:
  
  AliMergeableCollection* fHistogramCollection; //! collection of histograms
  AliCounterCollection* fEventCounters; //! event counters
  mutable AliAnalysisMuMuBinning* fBinning; // binning for particles

  mutable AliAnalysisMuMuCutRegistry* fCutRegistry; // cuts (owner)
  
  TString fBeamYear; // beam and year
  
  TList* fHistogramToDisable; // list of regexp of histo name(s) to disable
  
  TObjArray* fSubAnalysisVector; // list of companion analysis
  
  ClassDef(AliAnalysisTaskMuMu,26) // a class to analyse muon pairs (and single also ;-) )
};

#endif

