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

AliFemtoEventCutCentrality::AliFemtoEventCutCentrality(const Parameters &param):
  AliFemtoEventCut()
  , fCentralityType(param.centrality_type)
  , fEventCentrality(param.centrality_range)
  , fVertZPos(param.vertex_z_range)
  , fPsiEP(param.psi_ep_range)
  , fSelectTrigger(param.select_trigger)
  , fNEventsPassed(0)
  , fNEventsFailed(0)
{
}

AliFemtoConfigObject AliFemtoEventCutCentrality::GetConfigObject() const
{
  return AliFemtoConfigObject::BuildMap()
    ("_class", "AliFemtoEventCutCentrality")
    ("centrality_type", fCentralityType)
    ("centrality_range", fEventCentrality)
    ("vertex_z_range", fVertZPos)
    ("psi_ep_range", fPsiEP)
    ("select_trigger", fSelectTrigger);
}

AliFemtoConfigObject* AliFemtoEventCutCentrality::GetConfigObjectPtr() const
{
  auto cfg = new AliFemtoConfigObject(GetConfigObject());
  return cfg;
}


AliFemtoEventCutCentrality::Parameters::Parameters():
  centrality_type(kV0)
  , centrality_range(0, 100)
  , vertex_z_range(-100.0, 100.0)
  , psi_ep_range(-1000.0, 1000.0)
  , select_trigger(0)
{
}

AliFemtoEventCutCentrality::Parameters::Parameters(AliFemtoConfigObject &cfg)
  : Parameters()
{
  TString ctype;
  int ctype_int;
  if (cfg.pop_and_load("centrality_type", ctype)) {
    centrality_type = CentralityTypeFromName(ctype);
  }
  else if (cfg.pop_and_load("centrality_type", ctype_int)) {
    centrality_type = static_cast<CentralityType>(ctype_int);
  }
  cfg.pop_all()
    ("centrality_range", centrality_range)
    ("vertex_z_range", vertex_z_range)
    ("psi_ep_range", psi_ep_range)
    ("select_trigger", select_trigger);
}

// AliFemtoEventCutCentrality::Parameters::Parameters(const AliFemtoConfigObject& cfg)
//   : Parameters(AliFemtoConfigObject(cfg))
// {
// }

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

    new TObjString(prefix + Form("AliFemtoEventCutCentrality.centrality.min=%f", fEventCentrality.first)),
    new TObjString(prefix + Form("AliFemtoEventCutCentrality.centrality.max=%f", fEventCentrality.second)),
    new TObjString(prefix + Form("AliFemtoEventCutCentrality.vertex.min=%f", fVertZPos.first)),
    new TObjString(prefix + Form("AliFemtoEventCutCentrality.vertex.max=%f", fVertZPos.second)),
    new TObjString(prefix + Form("AliFemtoEventCutCentrality.psiep.min=%f", fPsiEP.first)),
    new TObjString(prefix + Form("AliFemtoEventCutCentrality.psiep.max=%f", fPsiEP.second)),
    new TObjString(prefix + Form("AliFemtoEventCutCentrality.trigger=%d", fSelectTrigger)),

  nullptr);

  return settings;
}

AliFemtoString AliFemtoEventCutCentrality::Report()
{
  /// Prepare report
  AliFemtoString
    report = std::string("AliFemtoEventCutCentrality Report\n")
           + Form("Centraltiy:\t %f - %f\n", fEventCentrality.first, fEventCentrality.second)
           + Form("Psi Ep:\t %E - %E\n", fPsiEP.first, fPsiEP.second)
           + Form("Vertex Z-position:\t %E - %E\n", fVertZPos.first, fVertZPos.second)
           + Form("Number of events which passed:\t%lu  Number which failed:\t%lu\n", fNEventsPassed, fNEventsFailed);

  return report;
}
