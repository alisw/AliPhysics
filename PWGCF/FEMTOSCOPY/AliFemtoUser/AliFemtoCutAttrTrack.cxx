///
/// \file AliFemtoUser/AliFemtoCutAttrTrack.cxx
///


#include "AliFemtoCutAttrTrack.h"


const std::pair<double, double>
  pwgfemto::TrackCutAttrSigmaPion::DEFAULT(-1.0, 2.0),
  pwgfemto::TrackCutAttrSigmaProton::DEFAULT(-1.0, 2.0),
  pwgfemto::TrackCutAttrMostProbablePion::DEFAULT(-1.0, 2.0),
  pwgfemto::TrackCutAttrEta::DEFAULT(-0.8, 0.8),
  pwgfemto::TrackCutAttrRapidity::DEFAULT(-0.8, 0.8),
  pwgfemto::TrackCutAttrMomentum::DEFAULT(0, 100.0),
  pwgfemto::TrackCutAttrPt::DEFAULT(0, 100.0);

const int
  pwgfemto::TrackCutAttrCharge::DEFAULT = 1;
