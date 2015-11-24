///
/// \file AliFemtoEventCutCentrality.cxx
///

#include "AliFemtoEventCutCentrality.h"
#include "TObjString.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoEventCutCentrality);
  /// \endcond
#endif

AliFemtoEventCutCentrality::AliFemtoEventCutCentrality():
  AliFemtoEventCut()
  , fCentralityType(kV0)
  , fSelectTrigger(0)
  , fNEventsPassed(0)
  , fNEventsFailed(0)
{
  /// Default constructor

  fEventCentrality[0] = 0;
  fEventCentrality[1] = 100;

  fVertZPos[0] = -100.0;
  fVertZPos[1] = 100.0;

  fPsiEP[0] = -1000.0;
  fPsiEP[1] = 1000.0;
}
//------------------------------
AliFemtoEventCutCentrality::~AliFemtoEventCutCentrality()
{ // Default destructor
}

//------------------------------
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

    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.centrality.min=%f", fEventCentrality[0])),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.centrality.max=%f", fEventCentrality[1])),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.vertex.min=%f", fVertZPos[0])),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.vertex.max=%f", fVertZPos[1])),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.psiep.min=%f", fPsiEP[0])),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.psiep.max=%f", fPsiEP[1])),
    new TObjString(prefix + TString::Format("AliFemtoEventCutCentrality.trigger=%d", fSelectTrigger)),

  NULL);

  return settings;
}

//------------------------------
AliFemtoString AliFemtoEventCutCentrality::Report()
{
  /// Prepare report
  TString report = TString::Format("Centraltiy:\t %f - %f\n", fEventCentrality[0], fEventCentrality[1]);

   report += TString::Format("Vertex Z-position:\t %E - %E\n", fVertZPos[0], fVertZPos[1])
           + TString::Format("Number of events which passed:\t%ld  Number which failed:\t%ld\n", fNEventsPassed, fNEventsFailed);

  return AliFemtoString(report);
}
