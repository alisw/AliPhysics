///
/// \file AliFemtoUser/AliFemtoCutAttrTrack.cxx
///


#include "AliFemtoCutAttrTrack.h"


const std::pair<double, double>
  pwgfemto::TrackCutAttrEta::DEFAULT(-0.8, 0.8),
  pwgfemto::TrackCutAttrMomentum::DEFAULT(0, 100.0),
  pwgfemto::TrackCutAttrPt::DEFAULT(0, 100.0);

const int
  pwgfemto::TrackCutAttrCharge::DEFAULT = 1;
