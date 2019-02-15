///
/// \file AliFemtoUser/AliFemtoAnalysisPionPionCuts.cxx
///

#include <TList.h>

#include "AliFemtoAnalysisPionPionCuts.h"

void
AliFemtoEventCutPionPionAK::AppendSettings(TCollection &settings) const
{
}

void
AliFemtoTrackCutPionPionAK::AppendSettings(TCollection &settings) const
{
}

void
AliFemtoPairCutPionPionAKDetaDphi::AppendSettings(TCollection &settings) const
{
  TString prefix = "AliFemtoPairCutPionPionAKDetaDphi.";
  settings.AddVector(
    new TObjString(prefix + Form("max_share_fraction=%g", max_share_fraction)),
    new TObjString(prefix + Form("max_share_quality=%g", max_share_quality)),
    new TObjString(prefix + Form("min_delta_eta=%g", min_delta_eta)),
    new TObjString(prefix + Form("min_delta_phi=%g", min_delta_phi)),
    new TObjString(prefix + Form("radius=%g", radius)),
    nullptr
  );
}

void
AliFemtoPairCutPionPionAKAvgSep::AppendSettings(TCollection &settings) const
{
  TString prefix = "AliFemtoPairCutPionPionAKAvgSep.";
  settings.AddVector(
    new TObjString(prefix + Form("max_share_fraction=%g", max_share_fraction)),
    new TObjString(prefix + Form("max_share_quality=%g", max_share_quality)),
    new TObjString(prefix + Form("min_avgsep=%g", min_avgsep)),
    nullptr
  );
}

/*
AliFemtoString
AliFemtoEventCutPionPionAK::Report()
{
  AliFemtoString report;
  return report;
}


AliFemtoString
AliFemtoTrackCutPionPionAK::Report()
{
  AliFemtoString report;
  return report;
}
*/


AliFemtoString
AliFemtoPairCutPionPionAKAvgSep::Report()
{
  AliFemtoString report;
  return report;
}


AliFemtoPairCutPionPionAKAvgSep::~AliFemtoPairCutPionPionAKAvgSep()
{}

AliFemtoString
AliFemtoPairCutPionPionAKDetaDphi::Report()
{
  AliFemtoString report;
  return report;
}

void
AliFemtoPairCutPionPionAKDetaDphi::EventBegin(const AliFemtoEvent *ev)
{
  fCurrentMagneticField = ev->MagneticField();

  // Correct AliFemto units error for magnetic field (back to Tesla)
  // TODO: Fix this bug in AliFemtoEventReaderAOD::CopyAODtoFemtoEvent
  if (fabs(fCurrentMagneticField) < 1e-10) {
    fCurrentMagneticField *= 1e13;
  }
  Super::EventBegin(ev);
}

/*

template <>
AliFemtoEventCutPionPionAK*
AliFemtoConfigObject::Into<AliFemtoEventCutPionPionAK>(bool)
{
  auto *cut = AliFemtoEventCutPionPionAK();
  return cut;
}

template <>
AliFemtoTrackCutPionPionAK*
AliFemtoConfigObject::Into<AliFemtoTrackCutPionPionAK>(bool)
{
  auto *cut = AliFemtoTrackCutPionPionAK();
  return cut;
}

template <>
AliFemtoConfigObject
AliFemtoConfigObject::From(const AliFemtoPairCutPionPionAK &cut)
{
  AliFemtoConfigObject cfg = AliFemtoConfigObject::BuildMap()
                              ("class", "AliFemtoPairCutPionPionAK");

  cut.StoreConfiguration(cfg);
  return cfg;
}
*/
