///
/// \file AliFemtoEventCutCentrality.cxx
///

#include "AliFemtoEventCutCentrality.h"
#include "TObjString.h"


AliFemtoEventCutCentrality::AliFemtoEventCutCentrality():
  AliFemtoEventCut()
  , fCentralityType(kV0)
  , fEventCentrality(0, 100)
  , fVertZPos(-100.0, 100.0)
  , fPsiEP(-1000.0, 1000.0)
  , fSelectTrigger(0)
  , fNEventsPassed(0)
  , fNEventsFailed(0)
{
}

bool AliFemtoEventCutCentrality::Pass(const AliFemtoEvent* event)
{
  /// Pass events if they fall within the multiplicity and z-vertex
  /// position range. Fail otherwise
  ///  int mult =  event->NumberOfTracks();
  const bool passes = PassCentrality(event)
                   && PassVertex(event)
                   && PassTrigger(event);

  passes ? fNEventsPassed++ : fNEventsFailed++;

  return passes;
}

TList* AliFemtoEventCutCentrality::AppendSettings(TList *settings,
                                                  const TString &prefix) const
{
  settings->AddVector(

    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.centrality.min=%f", fEventCentrality.first)),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.centrality.max=%f", fEventCentrality.second)),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.vertex.min=%f", fVertZPos.first)),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.vertex.max=%f", fVertZPos.second)),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.psiep.min=%f", fPsiEP.first)),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.psiep.max=%f", fPsiEP.second)),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.trigger=%d", fSelectTrigger)),

  NULL);

  return settings;
}

AliFemtoString AliFemtoEventCutCentrality::Report()
{
  /// Prepare report
  const TString
    report = TString::Format("Centraltiy:\t %f - %f\n", fEventCentrality.first, fEventCentrality.second)
           + TString::Format("Psi Ep:\t %E - %E\n", fPsiEP.first, fPsiEP.second)
           + TString::Format("Vertex Z-position:\t %E - %E\n", fVertZPos.first, fVertZPos.second)
           + TString::Format("Number of events which passed:\t%lu  Number which failed:\t%lu\n", fNEventsPassed, fNEventsFailed);

  return AliFemtoString(report);
}
