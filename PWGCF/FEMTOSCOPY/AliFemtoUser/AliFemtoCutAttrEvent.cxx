///
/// \class AliFemtoUser/AliFemtoCutAttrEvent.cxx
///

#include "AliFemtoCutAttrEvent.h"

const std::pair<double, double>
  pwgfemto::EventCutAttrCentrality::DEFAULT = {0.0, 100.0},
  pwgfemto::EventCutAttrEpPsi::DEFAULT = {-1000.0, 1000.0};

const std::pair<int, int>
  pwgfemto::EventCutAttrMultiplicty::DEFAULT = {0, 100000};
