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
AliFemtoTrackCutPionPionAK::AppendSettings(TCollection &settings, TString prefix) const
{
  settings.AddVector(
    new TObjString(prefix + Form("eta_min=%f", eta_range.first)),
    new TObjString(prefix + Form("eta_max=%f", eta_range.second)),
    nullptr
  );
}

void
AliFemtoPairCutPionPionAKDetaDphi::AppendSettings(TCollection &settings) const
{
  TString prefix = "AliFemtoPairCutPionPionAKDetaDphi.";
  settings.AddVector(
    new TObjString(prefix + Form("max_share_fraction=%g", share_fraction_max)),
    new TObjString(prefix + Form("max_share_quality=%g", share_quality_max)),
    new TObjString(prefix + Form("min_delta_eta=%g", delta_eta_min)),
    new TObjString(prefix + Form("min_delta_phi=%g", delta_phistar_min)),
    new TObjString(prefix + Form("phistar_radius=%g", phistar_radius)),
    nullptr
  );
}

void
AliFemtoPairCutPionPionAKAvgSep::AppendSettings(TCollection &settings) const
{
  TString prefix = "AliFemtoPairCutPionPionAKAvgSep.";
  settings.AddVector(
    new TObjString(prefix + Form("max_share_fraction=%g", share_fraction_max)),
    new TObjString(prefix + Form("max_share_quality=%g", share_quality_max)),
    new TObjString(prefix + Form("min_avgsep=%g", avgsep_min)),
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

#define IMPL_INTO_CUT(__name) \
  template <> __name* AliFemtoConfigObject::Into<__name>(bool) { \
    auto *cut = new __name(*this); return cut; }

#define IMPL_FROM_CUT(__name) \
  template <> \
  AliFemtoConfigObject AliFemtoConfigObject::From(const __name &cut) { \
    return cut.GetConfiguration(); }

#define IMPL_TOFROM_CUT(T) \
  IMPL_INTO_CUT(T) \
  IMPL_FROM_CUT(T)


IMPL_TOFROM_CUT(AliFemtoEventCutPionPionAK)
IMPL_TOFROM_CUT(AliFemtoTrackCutPionPionAK)
IMPL_TOFROM_CUT(AliFemtoTrackCutPionPionIdealAK)
IMPL_TOFROM_CUT(AliFemtoTrackCutPionPionMisidentAK)
IMPL_TOFROM_CUT(AliFemtoPairCutPionPionAKAvgSep)
IMPL_TOFROM_CUT(AliFemtoPairCutPionPionAKDetaDphi)
