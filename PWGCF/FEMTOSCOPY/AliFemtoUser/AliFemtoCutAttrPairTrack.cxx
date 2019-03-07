///
/// \file AliFemtoUser/AliFemtoCutAttrPairTrack.h
///


#include "AliFemtoCutAttrPairTrack.h"



const std::pair<double, double>
  pwgfemto::PairCutTrackAttrPt::DEFAULT = {0, 10.0};

/*
template<>
pwgfemto::AliFemtoPairCutPionPionAK*
AliFemtoConfigObject::Into(bool)
{
  return nullptr;
}

template<>
AliFemtoConfigObject
AliFemtoConfigObject::From(const pwgfemto::AliFemtoPairCutPionPionAK &cut)
{
  AliFemtoConfigObject cfg = AliFemtoConfigObject::BuildMap()
                              ("_class", "AliFemtoPairCutPionPionAK");
  cut.FillConfiguration(cfg);
  return cfg;
}
*/
