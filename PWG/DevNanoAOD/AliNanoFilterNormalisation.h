#ifndef ALINANOFILTERNORMALISATION_H
#define ALINANOFILTERNORMALISATION_H

#include <TNamed.h>
#include <TString.h>

class TCollection;
class TH2D;

class AliNanoFilterNormalisation : public TNamed {
  public:
  enum NormBin {
    kAnyEvent = 0,
    kTriggeredEvent,
    kTriggeredEventWithQualityCuts,
    kTriggeredEventWithQualityCutsAndRecoVertex,
    kAnalysisEvent
  };

  AliNanoFilterNormalisation(TString name = "NanoFilterNormalisation", TString title = "NanoFilterNormalisation", int nMultBins = 101, float multBegin = -1, float multEnd = 100);
  AliNanoFilterNormalisation(TString name, TString title, int nMultBins, float* mBins);
  ~AliNanoFilterNormalisation();

  void FillCandidate(bool triggered, bool nonVertexRelatedSel, bool recoVertex, bool allCuts, float mult = -.5);
  void FillSelected(bool triggered, bool nonVertexRelatedSel, bool recoVertex, bool allCuts, float mult = -.5);

  Long64_t Merge(TCollection* col);

  double GetScalingFactor(NormBin bin, float mult = 0.5);
  double GetNcanditateEvents(NormBin bin, float mult = 0.5);
  double GetNselectedEvents(NormBin bin, float mult = 0.5);

  const TH2D*  GetCandidateEventsHistogram() const { return fCandidateEvents; }
  const TH2D*  GetSelectedEventsHistogram() const { return fSelectedEvents; }

  private:
  TH2D* fCandidateEvents; ///->
  TH2D* fSelectedEvents;  ///->

  AliNanoFilterNormalisation& operator=(const AliNanoFilterNormalisation&);
  AliNanoFilterNormalisation(const AliNanoFilterNormalisation&);

  ClassDef(AliNanoFilterNormalisation,1);

};

#endif