///
/// \class FemtoAnalysisPionPionObjectConstructor.cxx
///

#include "AliFemtoConfigObject.h"

#include "AliFemtoEventReader.h"
#include "AliFemtoEventReaderAODMultSelection.h"
// #include "AliFemtoEventReaderAlt.h"
#include "AliFemtoEventReaderESD.h"
#include "AliFemtoEventReaderESDChain.h"

#include "AliFemtoEventCut.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoEventCutCentrality.h"

#include "AliFemtoTrackCut.h"
#include "AliFemtoBasicTrackCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoAODTrackCut.h"
#include "AliFemtoV0TrackCut.h"

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutDetaDphi.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoDummyPairCut.h"
#include "AliFemtoPairCutPt.h"

#include "AliFemtoAnalysisPionPionCuts.h"

#include "AliFemtoCorrFctn.h"


#include "AliFemtoAnalysisPionPionObjectConstructor.h"


#if __cplusplus < 201103L

  #define RETURN_IF_CAST(__type)
  #define TRY_CONSTRUCTING_CLASS(__name)
  #define FORWARD_TO_BUILDER(__type, __name)

#else

  #define RETURN_IF_CAST(__type) \
    if (auto *ptr = dynamic_cast<const __type*>(&obj)) { \
      return Configuration<__type>::GetConfigurationOf(*ptr); }

  #define TRY_CONSTRUCTING_CLASS(__name) \
    if (classname == #__name) {          \
      return Configuration<__name>(*this); }

  #define FORWARD_TO_BUILDER(__type, __name) \
    if (classname == #__name) { return Into<__type>(); }

#endif

#define RETURN_IF_CAST_FROM(__type) \
  if (auto *ptr = dynamic_cast<const __type*>(&obj)) { \
    return AliFemtoConfigObject::From(*ptr); }

#define TRY_CONSTRUCTING_INTO(__type) \
    if (classname == #__type) {          \
      return Into<__type>(true); }


template<>
AliFemtoConfigObject
AliFemtoConfigObject::From<AliFemtoEventReader>(const AliFemtoEventReader &obj)
{
  RETURN_IF_CAST(AliFemtoEventReaderAODMultSelection);
  RETURN_IF_CAST(AliFemtoEventReaderAODChain);
  RETURN_IF_CAST(AliFemtoEventReaderAOD);
  // RETURN_IF_CAST(AliFemtoEventReaderAlt);
  // RETURN_IF_CAST(AliFemtoEventReaderESDChain);
  // RETURN_IF_CAST(AliFemtoEventReaderESD);

  return AliFemtoConfigObject::BuildMap()
            ("class", "AliFemtoEventCut");
}


template<>
AliFemtoEventReader*
AliFemtoConfigObject::Into<AliFemtoEventReader>(bool)
{
  std::string classname;
  if (!find_and_load("class", classname)) {
    TString msg = "Could not load string-property 'class' from object:\n" + Stringify(true);
    std::cerr << "[AliFemtoAnalysisPionPion::ConstructPairCut] " << msg;
    return nullptr;
  }

  TRY_CONSTRUCTING_CLASS(AliFemtoEventReaderAODMultSelection);
  TRY_CONSTRUCTING_CLASS(AliFemtoEventReaderAODChain);
  TRY_CONSTRUCTING_CLASS(AliFemtoEventReaderAOD);
  // TRY_CONSTRUCTING_CLASS(AliFemtoEventReaderESDChain)
  // TRY_CONSTRUCTING_CLASS(AliFemtoEventReaderESD)

  Warning("AliFemtoConfigObject::Construct<ConstructEventReader>",
          "Could not load class %s", classname.c_str());

  return nullptr;
}


//
//  EVENT CUTS
//

template<>
AliFemtoConfigObject
AliFemtoConfigObject::From<AliFemtoEventCut>(const AliFemtoEventCut &obj)
{
  RETURN_IF_CAST_FROM(AliFemtoEventCutPionPionAK)
  RETURN_IF_CAST(AliFemtoBasicEventCut);
  RETURN_IF_CAST(AliFemtoEventCutCentrality);

  return AliFemtoConfigObject::BuildMap()
            ("class", "AliFemtoEventCut");
}


template<>
AliFemtoEventCut*
AliFemtoConfigObject::Into<AliFemtoEventCut>(bool)
{
  std::string classname;
  if (!find_and_load("class", classname)) {
    std::cerr
      << "[AliFemtoAnalysisPionPion::ConstructPairCut] "
      << "Could not load string-property 'class' from object:\n"
      << Stringify(true);
    return nullptr;
  }

  TRY_CONSTRUCTING_INTO(AliFemtoEventCutPionPionAK)
  TRY_CONSTRUCTING_CLASS(AliFemtoBasicEventCut);
  TRY_CONSTRUCTING_CLASS(AliFemtoEventCutCentrality);

  Warning("AliFemtoConfigObject::Construct<AliFemtoEventCuut>",
          "Could not load class %s", classname.c_str());

  return nullptr;
}


//
//  TRACK CUTS
//

template<>
AliFemtoConfigObject
AliFemtoConfigObject::From<AliFemtoTrackCut>(const AliFemtoTrackCut &obj)
{
  RETURN_IF_CAST_FROM(AliFemtoTrackCutPionPionAK);
  RETURN_IF_CAST(AliFemtoESDTrackCut);
  RETURN_IF_CAST(AliFemtoAODTrackCut);
  // RETURN_IF_CAST(AliFemtoBasicTrackCut);
  // RETURN_IF_CAST(AliFemtoESDTrackCutNSigmaFilter);
      // AliFemtoKKTrackCut
      // AliFemtoKKTrackCutFull
      // AliFemtoKKTrackCutTest
      // AliFemtoKpm45TrackCut
      // AliFemtoKpmTrackCut
      // AliFemtoMCTrackCut
      // AliFemtoMJTrackCut
      // AliFemtoQATrackCut

  return AliFemtoConfigObject::BuildMap()
            ("class", "AliFemtoTrackCut");
}


template<>
AliFemtoTrackCut*
AliFemtoConfigObject::Into<AliFemtoTrackCut>(bool)
{
  std::string classname;
  if (!find_and_load("class", classname)) {
    std::cerr << "[AliFemtoAnalysisPionPion::Into] "
              << "Could not load string-property 'class' from object:\n"
              << Stringify(true)
              << "\n";
    return nullptr;
  }

  TRY_CONSTRUCTING_INTO(AliFemtoTrackCutPionPionAK);
  TRY_CONSTRUCTING_CLASS(AliFemtoAODTrackCut);
  TRY_CONSTRUCTING_CLASS(AliFemtoESDTrackCut);

  Warning("AliFemtoConfigObject::Into<ConstructAliFemtoTrackCut>",
          "Could not load class %s", classname.c_str());

  return nullptr;
}



//
//  V0 CUTS
//

// template<>
// AliFemtoConfigObject
// AliFemtoConfigObject::From<AliFemtoV0Cut>(const AliFemtoV0Cut &obj)
// {
// }


// template<>
// AliFemtoV0Cut*
// AliFemtoConfigObject::Into<AliFemtoV0Cut>(bool)
// {
//   return nullptr;
// }

// template<>
// AliFemtoV0TrackCut*
// AliFemtoConfigObject::Into<AliFemtoV0TrackCut>(bool)
// {
//   return nullptr;
// }


//
//  PARTICLE CUTS (forwards to Track/V0)
//


//
//  PARTICLE CUTS (forwards to Track/V0)
//

template<>
AliFemtoConfigObject
AliFemtoConfigObject::From<AliFemtoParticleCut>(const AliFemtoParticleCut &obj)
{
  if (auto *track = dynamic_cast<const AliFemtoTrackCut*>(&obj)) {
    return From(*track);
  }

  if (auto *v0 = dynamic_cast<const AliFemtoV0*>(&obj)) {
    return From(*v0);
  }

  return AliFemtoConfigObject::BuildMap()
            ("class", "AliFemtoParticleCut");
}

template<>
AliFemtoParticleCut*
AliFemtoConfigObject::Into<AliFemtoParticleCut>(bool)
{
  std::string classname;
  if (!find_and_load("class", classname)) {
    TString msg = "Could not load string-property 'class' from object:\n" + Stringify(true);
    std::cerr << "[AliFemtoAnalysisPionPion::AliFemtoParticleCut] " << msg;
    return nullptr;
  }

  FORWARD_TO_BUILDER(AliFemtoTrackCut, AliFemtoTrackCutPionPionAK);
  FORWARD_TO_BUILDER(AliFemtoTrackCut, AliFemtoESDTrackCut);
  FORWARD_TO_BUILDER(AliFemtoTrackCut, AliFemtoAODTrackCut);
  FORWARD_TO_BUILDER(AliFemtoV0TrackCut, AliFemtoV0TrackCut);
  // FORWARD_TO_BUILDER(AliFemtoV0TrackCut, AliFemtoXiTrackCut);

  Warning("AliFemtoConfigObject::Construct<AliFemtoParticleCut>",
          "Could not load class '%s'", classname.c_str());

  return nullptr;
}


//
// PAIR CUTS
//

template<>
AliFemtoConfigObject
AliFemtoConfigObject::From<AliFemtoPairCut>(const AliFemtoPairCut &obj)
{
  RETURN_IF_CAST_FROM(AliFemtoPairCutPionPionAKDetaDphi)
  RETURN_IF_CAST_FROM(AliFemtoPairCutPionPionAKAvgSep)
  RETURN_IF_CAST(AliFemtoPairCutAntiGamma)
  RETURN_IF_CAST(AliFemtoPairCutDetaDphi)
  RETURN_IF_CAST(AliFemtoShareQualityPairCut)
  RETURN_IF_CAST(AliFemtoDummyPairCut)

  return AliFemtoConfigObject::BuildMap()
            ("class", "AliFemtoPairCut");
}


template<>
AliFemtoPairCut*
AliFemtoConfigObject::Into<AliFemtoPairCut>(bool)
{
  std::string classname;
  if (!find_and_load("class", classname)) {
    TString msg = "Could not load string-property 'class' from object:\n" + Stringify(true);
    std::cerr << "[AliFemtoAnalysisPionPion::ConstructPairCut] " << msg;
    return nullptr;
  }

  TRY_CONSTRUCTING_INTO(AliFemtoPairCutPionPionAKDetaDphi)
  TRY_CONSTRUCTING_INTO(AliFemtoPairCutPionPionAKAvgSep)

  TRY_CONSTRUCTING_CLASS(AliFemtoPairCutAntiGamma)
  TRY_CONSTRUCTING_CLASS(AliFemtoPairCutDetaDphi)
  TRY_CONSTRUCTING_CLASS(AliFemtoShareQualityPairCut)
  TRY_CONSTRUCTING_CLASS(AliFemtoDummyPairCut)

  Warning("AliFemtoConfigObject::Construct<ConstructAliFemtoPairCut>",
          "Could not load class %s", classname.c_str());

  return nullptr;
}


// implement various standard AliFemtoConfigObj <-> Configuration<T> <-> Pointer
// functions
#define IMPL_CFG_INTO_OBJ(T) \
  template <>                \
  T* AliFemtoConfigObject::Into<T>(bool _debug) { \
    Configuration<T> cfg(*this);  \
    T *cut = new T(); cfg.Configure(*cut); return cut; }


IMPL_CFG_INTO_OBJ(AliFemtoEventCutCentrality)
IMPL_CFG_INTO_OBJ(AliFemtoBasicEventCut)
IMPL_CFG_INTO_OBJ(AliFemtoESDTrackCut)

#undef IMPL_CFG_INTO_OBJ


#if __cplusplus >= 201103L

void
AbstractConfiguration<AliFemtoCorrFctn>::Configure(AliFemtoCorrFctn &cf) const
{
  if (!pair_cut_cfg.is_empty()) {
    cf.SetPairSelectionCut(pair_cut_cfg.Clone().Into<AliFemtoPairCut>());
  }
}

#endif
